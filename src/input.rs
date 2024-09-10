//! Input for a generic simulation.

use clap::ValueEnum;
use serde::{Deserialize, Serialize};

use crate::program::{CurieTemp, HysteresisLoop, Relax};

#[derive(Debug, Clone, ValueEnum, Serialize, Deserialize)]
pub enum Model {
    /// Ising model
    Ising,
    /// Heisenberg model
    Heisenberg,
}

impl Default for Model {
    fn default() -> Self {
        Model::Ising
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum UnitCellName {
    /// Simple cubic
    SC,
    /// Body-centered cubic
    BCC,
    /// Face-centered cubic
    FCC,
    /// Hexagonal close-packed
    HCP,
}

impl Default for UnitCellName {
    fn default() -> Self {
        UnitCellName::SC
    }
}

/// Unit cell to simulate.
#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub enum UnitCell {
    /// Unit cell by name
    Name(UnitCellName),
    /// Unit cell by path
    Path(String),
}

impl Default for UnitCell {
    fn default() -> Self {
        UnitCell::Name(UnitCellName::default())
    }
}

/// Size to expand the unit cell.
#[derive(Debug, Deserialize, Serialize)]
pub struct UnitCellSize {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

impl Default for UnitCellSize {
    fn default() -> Self {
        UnitCellSize { x: 1, y: 1, z: 1 }
    }
}

/// Periodic boundary conditions.
#[derive(Debug, Deserialize, Serialize)]
pub struct PeriodicBoundaryConditions {
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl Default for PeriodicBoundaryConditions {
    fn default() -> Self {
        PeriodicBoundaryConditions {
            x: true,
            y: true,
            z: true,
        }
    }
}

/// Sample to simulate.
#[derive(Debug, Default, Deserialize, Serialize)]
pub struct Sample {
    /// Unit cell to create
    pub unitcell: UnitCell,
    /// Size to expand the unit cell
    pub size: UnitCellSize,
    /// Periodic boundary conditions
    pub pbc: PeriodicBoundaryConditions,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(tag = "program")]
pub enum Program {
    /// Relaxation
    Relax(Relax),
    /// Curie temperature
    CurieTemp(CurieTemp),
    /// Hysteresis loop
    Hysteresis(HysteresisLoop),
}

impl Default for Program {
    fn default() -> Self {
        Program::Relax(Relax::default())
    }
}

/// Input for a generic simulation.
#[derive(Debug, Deserialize, Serialize)]
pub struct Input {
    /// Model to simulate
    pub sample: Sample,
    /// Steps to take
    pub steps: Vec<Program>,
}

impl Default for Input {
    fn default() -> Self {
        Input {
            sample: Sample::default(),
            steps: vec![
                Program::Relax(Relax::default()),
                Program::CurieTemp(CurieTemp::default()),
            ],
        }
    }
}
