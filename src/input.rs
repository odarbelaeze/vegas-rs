//! Input for a generic simulation.

use clap::ValueEnum;
use serde::{Deserialize, Serialize};

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

#[derive(Debug, Deserialize, Serialize)]
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

#[derive(Debug, Default, Deserialize, Serialize)]
pub struct Create {
    /// Unit cell to create
    pub unitcell: UnitCell,
    pub size: UnitCellSize,
    pub pbc: PeriodicBoundaryConditions,
}

#[derive(Debug, Default, Deserialize, Serialize)]
pub struct Input {
    /// Model to simulate
    pub create: Create,
}
