//! Input for a generic simulation.

use clap::ValueEnum;
use rand::Rng;
use serde::{Deserialize, Serialize};
use vegas_lattice::Lattice;

use crate::{
    error::Result,
    hamiltonian::Exchange,
    integrator::MetropolisIntegrator,
    machine::Machine,
    program::{CoolDown, HysteresisLoop, Program, Relax},
    state::{HeisenbergSpin, IsingSpin, Spin, State},
};

#[derive(Debug, Default, Clone, ValueEnum, Serialize, Deserialize)]
pub enum Model {
    /// Ising model
    #[default]
    Ising,
    /// Heisenberg model
    Heisenberg,
}

#[derive(Clone, Default, Debug, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum UnitCellName {
    /// Simple cubic
    #[default]
    SC,
    /// Body-centered cubic
    BCC,
    /// Face-centered cubic
    FCC,
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
pub enum Step {
    /// Relaxation
    Relax(Relax),
    /// Curie temperature
    CoolDown(CoolDown),
    /// Hysteresis loop
    Hysteresis(HysteresisLoop),
}

impl Default for Step {
    fn default() -> Self {
        Step::Relax(Relax::default())
    }
}

/// Input for a generic simulation.
#[derive(Debug, Deserialize, Serialize)]
pub struct Input {
    /// Model to simulate
    pub model: Model,
    /// Sample to simulate
    pub sample: Sample,
    /// Steps to take
    pub steps: Vec<Step>,
}

impl Default for Input {
    fn default() -> Self {
        Input {
            model: Model::default(),
            sample: Sample::default(),
            steps: vec![
                Step::Relax(Relax::default()),
                Step::CoolDown(CoolDown::default()),
            ],
        }
    }
}

impl Input {
    /// Create a new input.
    pub fn new(model: Model, sample: Sample, steps: Vec<Step>) -> Self {
        Input {
            model,
            sample,
            steps,
        }
    }

    fn run_with_spin<T: Spin, R: Rng>(&self, rng: &mut R) -> Result<()> {
        let lattice = self.lattice();
        let integrator = MetropolisIntegrator::new();
        let hamiltonian = Exchange::from_lattice(&lattice);
        let mut machine = Machine::new(
            2.8,
            0.0,
            hamiltonian,
            integrator,
            State::<T>::rand_with_size(rng, lattice.sites().len()),
        );
        for program in self.steps.iter() {
            match program {
                Step::Relax(relax) => {
                    relax.run(rng, &mut machine)?;
                }
                Step::CoolDown(curie) => {
                    curie.run(rng, &mut machine)?;
                }
                Step::Hysteresis(hysteresis) => {
                    hysteresis.run(rng, &mut machine)?;
                }
            }
        }
        Ok(())
    }

    /// Get the lattice.
    pub fn lattice(&self) -> Lattice {
        let unitcell = match &self.sample.unitcell {
            UnitCell::Name(name) => match name {
                UnitCellName::SC => Lattice::sc(1.0),
                UnitCellName::BCC => Lattice::bcc(1.0),
                UnitCellName::FCC => Lattice::fcc(1.0),
            },
            UnitCell::Path(_path) => todo!(),
        };
        let UnitCellSize { x, y, z } = self.sample.size;
        let PeriodicBoundaryConditions {
            x: pbc_x,
            y: pbc_y,
            z: pbc_z,
        } = self.sample.pbc;
        let mut lattice = unitcell.expand(x, y, z);
        if !pbc_x {
            lattice = lattice.drop_x();
        }
        if !pbc_y {
            lattice = lattice.drop_y();
        }
        if !pbc_z {
            lattice = lattice.drop_z();
        }
        lattice
    }

    pub fn run<R: Rng>(&self, rng: &mut R) -> Result<()> {
        match self.model {
            Model::Ising => self.run_with_spin::<IsingSpin, _>(rng),
            Model::Heisenberg => self.run_with_spin::<HeisenbergSpin, _>(rng),
        }
    }
}
