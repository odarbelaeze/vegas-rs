//! IO module for reading and writing data to and from files

use crate::error::IOResult;
use crate::state::{Spin, State};
use crate::thermostat::Thermostat;

use std::fs::File;
use std::iter::repeat_n;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{BooleanArray, Float64Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::{arrow::ArrowWriter, basic::Compression, file::properties::WriterProperties};

pub struct ObservableParquetIO {
    schema: Arc<Schema>,
    writer: Option<ArrowWriter<File>>,
    step: usize,
    stage: u64,
}

impl ObservableParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P) -> IOResult<Self> {
        let file = File::create(path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("step", DataType::UInt64, false),
            Field::new("stage", DataType::UInt64, false),
            Field::new("relax", DataType::Boolean, false),
            Field::new("temperature", DataType::Float64, false),
            Field::new("field", DataType::Float64, false),
            Field::new("energy", DataType::Float64, false),
            Field::new("magnetization", DataType::Float64, false),
        ]));
        let properties = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(properties))?;
        Ok(Self {
            schema: schema.clone(),
            writer: Some(writer),
            step: 0,
            stage: 0,
        })
    }

    pub fn write(
        &mut self,
        thermostat: &Thermostat,
        energy: &[f64],
        magnetization: &[f64],
        relax: bool,
    ) -> IOResult<()> {
        debug_assert!(energy.len() == magnetization.len());
        // Let's think that every step needs to be recorded here
        let step: UInt64Array = (self.step..self.step + energy.len())
            .map(|i| i as u64)
            .collect();
        self.step += energy.len();

        // Each time we call write we enter a different stage of the run
        let stage: UInt64Array = repeat_n(self.stage, energy.len()).collect();
        self.stage += 1;

        let relax: BooleanArray = repeat_n(Some(relax), energy.len()).collect();
        let temp: Float64Array = repeat_n(thermostat.temperature(), energy.len()).collect();
        let field: Float64Array = repeat_n(thermostat.field(), energy.len()).collect();
        let energy: Float64Array = Float64Array::from(energy.to_owned());
        let magnetization: Float64Array = Float64Array::from(magnetization.to_owned());

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(step),
                Arc::new(stage),
                Arc::new(relax),
                Arc::new(temp),
                Arc::new(field),
                Arc::new(energy),
                Arc::new(magnetization),
            ],
        )?;

        if let Some(writer) = &mut self.writer {
            writer.write(&batch)?;
        } else {
            return Err(std::io::Error::other("Writer has been closed"))?;
        }
        Ok(())
    }
}

impl Drop for ObservableParquetIO {
    fn drop(&mut self) {
        if let Some(self_writer) = self.writer.take()
            && let Err(err) = self_writer.close()
        {
            eprintln!("error closing parquet writer: {}", err);
        }
    }
}

pub struct StateParquetIO {
    schema: Arc<Schema>,
    writer: Option<ArrowWriter<File>>,
}

impl StateParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P) -> IOResult<Self> {
        let file = File::create(path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("relax", DataType::Boolean, false),
            Field::new("stage", DataType::UInt64, false),
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
            writer: Some(writer),
        })
    }

    pub fn write<S: Spin>(
        &mut self,
        relax: bool,
        stage: usize,
        step: usize,
        state: &State<S>,
    ) -> IOResult<()> {
        let relax = BooleanArray::from(repeat_n(relax, state.len()).collect::<Vec<_>>());
        let stage = UInt64Array::from(repeat_n(stage as u64, state.len()).collect::<Vec<_>>());
        let step = UInt64Array::from(repeat_n(step as u64, state.len()).collect::<Vec<_>>());
        let id = UInt64Array::from((0..state.len()).map(|i| i as u64).collect::<Vec<_>>());
        let sx = Float64Array::from(state.spins().iter().map(|s| s.sx()).collect::<Vec<_>>());
        let sy = Float64Array::from(state.spins().iter().map(|s| s.sy()).collect::<Vec<_>>());
        let sz = Float64Array::from(state.spins().iter().map(|s| s.sz()).collect::<Vec<_>>());
        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(relax),
                Arc::new(stage),
                Arc::new(step),
                Arc::new(id),
                Arc::new(sx),
                Arc::new(sy),
                Arc::new(sz),
            ],
        )?;
        if let Some(writer) = &mut self.writer {
            writer.write(&batch)?;
        } else {
            return Err(std::io::Error::other("Writer has been closed"))?;
        }
        Ok(())
    }
}

impl Drop for StateParquetIO {
    fn drop(&mut self) {
        if let Some(writer) = self.writer.take()
            && let Err(err) = writer.close()
        {
            eprintln!("error closing parquet writer: {}", err);
        }
    }
}
