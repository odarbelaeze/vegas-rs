//! Module for writing simulation data to Parquet files.
//!
//! This module provides functionality to write observable data and spin state data
//! to Parquet files using the Apache Arrow format.
//! It defines two main structs: `ObservableParquetIO` and `StateParquetIO`,
//! each responsible for writing different types of data.

use crate::{
    error::IoResult,
    state::{Spin, State},
    thermostat::Thermostat,
};
use arrow::{
    array::{BooleanArray, Float64Array, UInt64Array},
    datatypes::{DataType, Field, Schema},
    record_batch::RecordBatch,
};
use parquet::{arrow::ArrowWriter, basic::Compression, file::properties::WriterProperties};
use std::{
    fs::{File, rename},
    iter::repeat_n,
    path::{Path, PathBuf},
    sync::Arc,
};

pub struct ObservableParquetIO {
    path: PathBuf,
    temp_path: PathBuf,
    schema: Arc<Schema>,
    writer: Option<ArrowWriter<File>>,
}

impl ObservableParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P) -> IoResult<Self> {
        let temp_path = path.as_ref().with_extension("parquet.tmp");
        let file = File::create(&temp_path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("relax", DataType::Boolean, false),
            Field::new("stage", DataType::UInt64, false),
            Field::new("step", DataType::UInt64, false),
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
            path: path.as_ref().to_path_buf(),
            temp_path,
            schema: schema.clone(),
            writer: Some(writer),
        })
    }

    pub fn write<S: Spin>(
        &mut self,
        relax: bool,
        stage: usize,
        thermostat: &Thermostat<S>,
        energy: &[f64],
        magnetization: &[f64],
    ) -> IoResult<()> {
        debug_assert!(energy.len() == magnetization.len());
        let step: UInt64Array = (0..energy.len()).map(|i| i as u64).collect();
        let stage: UInt64Array = repeat_n(stage as u64, energy.len()).collect();
        let relax: BooleanArray = repeat_n(Some(relax), energy.len()).collect();
        let temperature: Float64Array = repeat_n(thermostat.temperature(), energy.len()).collect();
        let field: Float64Array = repeat_n(thermostat.field().magnitude(), energy.len()).collect();
        let energy: Float64Array = Float64Array::from(energy.to_owned());
        let magnetization: Float64Array = Float64Array::from(magnetization.to_owned());

        let batch = RecordBatch::try_new(
            self.schema.clone(),
            vec![
                Arc::new(relax),
                Arc::new(stage),
                Arc::new(step),
                Arc::new(temperature),
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
        if let Some(self_writer) = self.writer.take() {
            if let Err(err) = self_writer.close() {
                eprintln!("error closing parquet writer: {}", err);
                return;
            }
            if let Err(err) = rename(&self.temp_path, &self.path) {
                eprintln!("error renaming parquet file: {}", err);
            }
        }
    }
}

pub struct StateParquetIO {
    path: PathBuf,
    temp_path: PathBuf,
    schema: Arc<Schema>,
    writer: Option<ArrowWriter<File>>,
}

impl StateParquetIO {
    pub fn try_new<P: AsRef<Path>>(path: P) -> IoResult<Self> {
        let temp_path = path.as_ref().with_extension("parquet.tmp");
        let file = File::create(&temp_path)?;
        let schema = Arc::new(Schema::new(vec![
            Field::new("relax", DataType::Boolean, false),
            Field::new("stage", DataType::UInt64, false),
            Field::new("step", DataType::UInt64, false),
            Field::new("temperature", DataType::Float64, false),
            Field::new("field", DataType::Float64, false),
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
            path: path.as_ref().to_path_buf(),
            temp_path,
            schema: schema.clone(),
            writer: Some(writer),
        })
    }

    pub fn write<S: Spin>(
        &mut self,
        relax: bool,
        stage: usize,
        step: usize,
        thermostat: &Thermostat<S>,
        state: &State<S>,
    ) -> IoResult<()> {
        let relax = BooleanArray::from(repeat_n(relax, state.len()).collect::<Vec<_>>());
        let stage = UInt64Array::from(repeat_n(stage as u64, state.len()).collect::<Vec<_>>());
        let step = UInt64Array::from(repeat_n(step as u64, state.len()).collect::<Vec<_>>());
        let temperature: Float64Array = repeat_n(thermostat.temperature(), state.len()).collect();
        let field: Float64Array = repeat_n(thermostat.field().magnitude(), state.len()).collect();
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
                Arc::new(temperature),
                Arc::new(field),
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
        if let Some(writer) = self.writer.take() {
            if let Err(err) = writer.close() {
                eprintln!("error closing parquet writer: {}", err);
                return;
            }
            if let Err(err) = rename(&self.temp_path, &self.path) {
                eprintln!("error renaming parquet file: {}", err);
            }
        };
    }
}
