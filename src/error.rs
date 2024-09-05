//! Errors for the vegas package

use std::io::Error as IoError;
use std::result::Result as StdResult;
use thiserror::Error;
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
}

/// Error type for program missconfiguration
#[derive(Error, Debug)]
pub enum ProgramError {
    #[error("maximum temperature must be greater than minimum temperature")]
    MaxTempLessThanMinTemp,
    #[error("there should be at least one step")]
    NoSteps,
    #[error("temperature must be greater than zero")]
    ZeroTemp,
    #[error("temperature delta must be greater than zero")]
    ZeroDelta,
}

/// Result type for the vegas package
pub type Result<T> = StdResult<T, VegasError>;
