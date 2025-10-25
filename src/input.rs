//! Input for a generic simulation.

use std::{io::stdout, path::PathBuf};

use clap::ValueEnum;
use rand::Rng;
use serde::{Deserialize, Serialize};
use vegas_lattice::Lattice;

use crate::{
    error::Result,
    hamiltonian::{Exchange, Hamiltonian, ZeemanEnergy},
    instrument::{Instrument, RawStatSensor, StatSensor, StateSensor},
    integrator::MetropolisIntegrator,
    machine::Machine,
    program::{CoolDown, HysteresisLoop, Program, Relax},
    state::{HeisenbergSpin, IsingSpin, Spin, State},
    thermostat::Thermostat,
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

/// State output for a simulation.
#[derive(Debug, Deserialize, Serialize)]
pub struct StateOutput {
    /// Write the states to a parquet file
    pub path: PathBuf,
    /// Frequency of writing states
    pub frequency: usize,
}

/// Output for a generic simulation.
#[derive(Debug, Deserialize, Serialize)]
pub struct Output {
    /// Write the raw data into the given file
    pub raw: Option<PathBuf>,
    /// Write states to a parquet file
    pub state: Option<StateOutput>,
}

impl Default for Output {
    fn default() -> Self {
        Self {
            raw: Some("./output.parquet".into()),
            state: Some(StateOutput {
                path: "./state.parquet".into(),
                frequency: 1000,
            }),
        }
    }
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
    /// Output for the simulation
    pub output: Option<Output>,
}

impl Default for Input {
    fn default() -> Self {
        Input {
            model: Default::default(),
            sample: Default::default(),
            steps: vec![
                Step::Relax(Relax::default()),
                Step::CoolDown(CoolDown::default()),
            ],
            output: Some(Output::default()),
        }
    }
}

pub struct InputBuilder {
    model: Option<Model>,
    sample: Option<Sample>,
    steps: Option<Vec<Step>>,
    output: Option<Output>,
}

impl InputBuilder {
    pub fn new() -> Self {
        InputBuilder {
            model: None,
            sample: None,
            steps: None,
            output: None,
        }
    }

    pub fn model(mut self, model: Model) -> Self {
        self.model = Some(model);
        self
    }

    pub fn sample(mut self, sample: Sample) -> Self {
        self.sample = Some(sample);
        self
    }

    pub fn steps(mut self, steps: Vec<Step>) -> Self {
        self.steps = Some(steps);
        self
    }

    pub fn output(mut self, output: Output) -> Self {
        self.output = Some(output);
        self
    }

    pub fn build(self) -> Input {
        Input {
            model: self.model.unwrap_or_default(),
            sample: self.sample.unwrap_or_default(),
            steps: self.steps.unwrap_or_default(),
            output: self.output,
        }
    }
}

impl Default for InputBuilder {
    fn default() -> Self {
        InputBuilder::new()
    }
}

impl Input {
    fn run_with_spin<S: Spin + 'static, R: Rng>(&self, rng: &mut R) -> Result<()> {
        let lattice = self.lattice();
        let integrator = MetropolisIntegrator::new();
        let hamiltonian =
            hamiltonian!(Exchange::from_lattice(&lattice), ZeemanEnergy::new(S::up()));
        let instruments = self.instruments::<_, S>()?;
        let mut machine = Machine::new(
            Thermostat::new(2.8, 0.0),
            hamiltonian,
            integrator,
            instruments,
            State::<S>::rand_with_size(rng, lattice.sites().len()),
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

    fn lattice(&self) -> Lattice {
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

    fn instruments<H: Hamiltonian<S> + 'static, S: Spin + 'static>(
        &self,
    ) -> Result<Vec<Box<dyn Instrument<H, S>>>> {
        let mut instruments: Vec<Box<dyn Instrument<_, _>>> =
            vec![Box::new(StatSensor::<_, S>::new(Box::new(stdout())))];
        if let Some(output) = &self.output
            && let Some(raw_filename) = &output.raw
        {
            instruments.push(Box::new(RawStatSensor::<_, S>::try_new(raw_filename)?));
        }
        if let Some(output) = &self.output
            && let Some(state_output) = &output.state
        {
            instruments.push(Box::new(StateSensor::<_, S>::try_new(
                &state_output.path,
                state_output.frequency,
            )?));
        }
        Ok(instruments)
    }

    pub fn run<R: Rng>(&self, rng: &mut R) -> Result<()> {
        match self.model {
            Model::Ising => self.run_with_spin::<IsingSpin, _>(rng),
            Model::Heisenberg => self.run_with_spin::<HeisenbergSpin, _>(rng),
        }
    }
}
