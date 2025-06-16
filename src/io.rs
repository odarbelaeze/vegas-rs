//! IO module for reading and writing data to and from files

use crate::error::IOResult as Result;
use crate::observables::Sensor;
use crate::state::{Spin, State};

use std::fs::File;
use std::iter::repeat;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{BooleanArray, Float64Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::{arrow::ArrowWriter, basic::Compression, file::properties::WriterProperties};

pub struct ObservableParquetIO {
    schema: Arc<Schema>,
    writer: ArrowWriter<File>,
    step: usize,
    stage: u64,
}

impl ObservableParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("step", DataType::UInt64, false),
            Field::new("stage", DataType::UInt64, false),
            Field::new("relax", DataType::Boolean, false),
            Field::new("beta", DataType::Float64, false),
            Field::new("energy", DataType::Float64, false),
            Field::new("magnetization", DataType::Float64, false),
        ]));
        let properties = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(properties))?;
        Ok(Self {
            schema: schema.clone(),
            writer,
            step: 0,
            stage: 0,
        })
    }

    pub fn write(&mut self, sensor: &Sensor, relax: bool) -> Result<()> {
        // Let's think that every step needs to be recorded here
        let step: UInt64Array = (self.step..self.step + sensor.len())
            .map(|i| i as u64)
            .collect();
        self.step += sensor.len();

        // Each time we call write we enter a different stage of the run
        let stage: UInt64Array = repeat(self.stage).take(sensor.len()).collect();
        self.stage += 1;

        // We should indicate if it's a relaxation step
        let relax: BooleanArray = repeat(Some(relax)).take(sensor.len()).collect();

        // Grab the beta value
        let beta: Float64Array = repeat(sensor.beta()).take(sensor.len()).collect();

        // Grab the measurement values
        let energy: Float64Array = Float64Array::from(sensor.energy().clone());
        let magnetization: Float64Array = Float64Array::from(sensor.magnetization().clone());

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(step),
                Arc::new(stage),
                Arc::new(relax),
                Arc::new(beta),
                Arc::new(energy),
                Arc::new(magnetization),
            ],
        )?;

        self.writer.write(&batch)?;

        Ok(())
    }

    pub fn close(self) -> Result<()> {
        self.writer.close()?;
        Ok(())
    }
}

pub struct StateParquetIO {
    schema: Arc<Schema>,
    writer: ArrowWriter<File>,
    step: usize,
    frequency: u64,
}

impl StateParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P, frequency: u64) -> Result<Self> {
        let file = File::create(path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("step", DataType::UInt64, false),
            Field::new("id", DataType::UInt64, false),
            Field::new("sx", DataType::Float64, false),
            Field::new("sy", DataType::Float64, false),
            Field::new("sz", DataType::Float64, false),
        ]));
        let properties = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(properties))?;
        Ok(Self {
            schema: schema.clone(),
            writer,
            step: 0,
            frequency,
        })
    }

    pub fn write<T: Spin>(&mut self, state: &State<T>) -> Result<()> {
        if (self.step as u64 % self.frequency) != 0 {
            self.step += 1;
            return Ok(());
        }
        let step = UInt64Array::from(
            repeat(self.step as u64)
                .take(state.len())
                .collect::<Vec<_>>(),
        );
        let id = UInt64Array::from((0..state.len()).map(|i| i as u64).collect::<Vec<_>>());
        let sx = Float64Array::from(state.spins().iter().map(|s| s.sx()).collect::<Vec<_>>());
        let sy = Float64Array::from(state.spins().iter().map(|s| s.sy()).collect::<Vec<_>>());
        let sz = Float64Array::from(state.spins().iter().map(|s| s.sz()).collect::<Vec<_>>());
        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(step),
                Arc::new(id),
                Arc::new(sx),
                Arc::new(sy),
                Arc::new(sz),
            ],
        )?;
        self.writer.write(&batch)?;
        Ok(())
    }
}
