//! Errors for the vegas package

use std::io::Error as IoError;
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
    IoError(#[from] IoError),
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
}

/// Result type for the vegas package
pub type Result<T> = StdResult<T, VegasError>;
