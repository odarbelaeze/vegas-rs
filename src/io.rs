//! IO module for reading and writing data to and from files

use crate::error::IOResult as Result;
use crate::observables::Reading;

use std::fs::File;
use std::iter::repeat;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{BooleanArray, Float64Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::{arrow::ArrowWriter, basic::Compression, file::properties::WriterProperties};

pub struct ParquetIO {
    schema: Arc<Schema>,
    writer: ArrowWriter<File>,
    step: usize,
    stage: u64,
}

impl ParquetIO {
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

    pub fn write(&mut self, readings: Vec<Reading>, relax: bool) -> Result<()> {
        // Let's think that every step needs to be recorded here
        let step: UInt64Array = (self.step..self.step + readings.len())
            .map(|i| i as u64)
            .collect();
        self.step += readings.len();

        // Each time we call write we enter a different stage of the run
        let stage: UInt64Array = repeat(self.stage).take(readings.len()).collect();
        self.stage += 1;

        // We should indicate if it's a relaxation step
        let relax: BooleanArray = repeat(Some(relax)).take(readings.len()).collect();

        // Let's build our data arrays
        let mut beta = Vec::<f64>::with_capacity(readings.len());
        let mut energy = Vec::<f64>::with_capacity(readings.len());
        let mut magnetization = Vec::<f64>::with_capacity(readings.len());
        readings.iter().for_each(|reading| {
            beta.push(reading.beta);
            energy.push(reading.energy);
            magnetization.push(reading.magnetization);
        });

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(step),
                Arc::new(stage),
                Arc::new(relax),
                Arc::new(Float64Array::from(beta)),
                Arc::new(Float64Array::from(energy)),
                Arc::new(Float64Array::from(magnetization)),
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
