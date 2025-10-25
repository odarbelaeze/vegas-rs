//! Errors for the vegas package

use arrow::error::ArrowError;
use parquet::errors::ParquetError;
use std::io::Error as StdIOError;
use std::result::Result as StdResult;
use thiserror::Error;
use toml::{de::Error as TomlDeserializeError, ser::Error as TomlSerializeError};
use vegas_lattice::error::VegasLatticeError;

/// Error type for the vegas package
#[derive(Error, Debug)]
pub enum VegasError {
    #[error("program error: {0}")]
    ProgramError(#[from] ProgramError),
    #[error("io error: {0}")]
    IOError(#[from] IOError),
    #[error("lattice error: {0}")]
    LatticeError(#[from] VegasLatticeError),
    #[error("toml deserialization error: {0}")]
    TomlDeserializeError(#[from] TomlDeserializeError),
    #[error("toml serialization error: {0}")]
    TomlSerializeError(#[from] TomlSerializeError),
}

/// Error type for program misconfiguration
#[derive(Error, Debug)]
pub enum ProgramError {
    #[error("maximum temperature must be greater than minimum temperature")]
    TemperatureMaxLessThanMin,
    #[error("there should be at least one step")]
    NoSteps,
    #[error("temperature must be greater than zero")]
    ZeroTemperature,
    #[error("cooling rate must be greater than zero")]
    ZeroCoolRate,
    #[error("maximum field must be greater than zero")]
    ZeroField,
    #[error("field step must be greater than zero")]
    ZeroFieldStep,
    #[error("machine error: {0}")]
    MachineError(#[from] MachineError),
}

// Error type for IO operations
#[derive(Error, Debug)]
pub enum IOError {
    #[error("std io error: {0}")]
    StdIOError(#[from] StdIOError),
    #[error("parquet error: {0}")]
    ParquetError(#[from] ParquetError),
    #[error("arrow error: {0}")]
    ArrowError(#[from] ArrowError),
}

// Error type for machine operations
#[derive(Error, Debug)]
pub enum MachineError {
    #[error("instrument error: {0}")]
    InstrumentError(#[from] InstrumentError),
}

// Error type for instrument operations
#[derive(Error, Debug)]
pub enum InstrumentError {
    #[error("std io error: {0}")]
    StdIOError(#[from] StdIOError),
    #[error("parquet error: {0}")]
    ParquetError(#[from] ParquetError),
    #[error("arrow error: {0}")]
    ArrowError(#[from] ArrowError),
}

/// Result type for the vegas package
pub type Result<T> = StdResult<T, VegasError>;

/// Result type for program misconfiguration
pub type ProgramResult<T> = StdResult<T, ProgramError>;

/// Result type for IO operations
pub type IOResult<T> = StdResult<T, IOError>;

/// Result type for machine operations
pub type MachineResult<T> = StdResult<T, MachineError>;

/// Result type for instrument operations
pub type InstrumentResult<T> = StdResult<T, InstrumentError>;
