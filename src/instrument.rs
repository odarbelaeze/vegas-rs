//! Instruments to monitor the simulation.
//!
//! This module defines instruments that can hook into the simulation process
//! to monitor and record various statistics and states during the simulation.
//! It includes instruments for recording statistical data and saving spin states
//! to Parquet files.

use crate::{
    accumulator::Accumulator,
    energy::Hamiltonian,
    error::{InstrumentResult, IoResult},
    io::{ObservableParquetIO, StateParquetIO},
    state::{Spin, State},
    thermostat::Thermostat,
};
use std::{io::Write, marker::PhantomData, path::Path};

/// An instrument allows to hook into the simulation at various points.
pub trait Instrument<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    /// Hook called when a relaxation starts.
    fn on_relax_start(
        &mut self,
        _thermostat: &Thermostat,
        _hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        Ok(())
    }

    /// Hook called when a relaxation ends.
    fn on_relax_end(&mut self) -> InstrumentResult<()> {
        Ok(())
    }

    /// Hook called when a measurement starts.
    fn on_measure_start(
        &mut self,
        _thermostat: &Thermostat,
        _hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        Ok(())
    }

    /// Hook called when a measurement ends.
    fn on_measure_end(&mut self) -> InstrumentResult<()> {
        Ok(())
    }

    /// Hook called after each integration step.
    fn after_step(&mut self, _state: &State<S>) -> InstrumentResult<()> {
        Ok(())
    }
}

/// An instrument that writes "statistics" to a given file.
pub struct StatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    output: Box<dyn Write>,
    energy_acc: Accumulator,
    magnetization_acc: Accumulator,
    thermostat: Option<Thermostat>,
    hamiltonian: Option<H>,
    phantom: PhantomData<S>,
}

impl<H, S> StatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    pub fn new(output: Box<dyn Write>) -> Self {
        Self {
            output,
            energy_acc: Accumulator::new(),
            magnetization_acc: Accumulator::new(),
            thermostat: None,
            hamiltonian: None,
            phantom: PhantomData,
        }
    }
}

impl<H, S> Instrument<H, S> for StatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    fn on_measure_start(
        &mut self,
        thermostat: &Thermostat,
        hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        self.thermostat = Some(thermostat.clone());
        self.hamiltonian = Some(hamiltonian.clone());
        Ok(())
    }

    fn on_measure_end(&mut self) -> InstrumentResult<()> {
        if let (Some(thermostat), Some(_)) = (&self.thermostat, &self.hamiltonian) {
            writeln!(
                self.output,
                "{:.16} {:.16} {:.16} {:.16} {:.16} {:.16} {:.16}",
                thermostat.temperature(),
                thermostat.field(),
                self.energy_acc.mean(),
                self.energy_acc.variance() / thermostat.temperature(),
                self.magnetization_acc.mean(),
                self.magnetization_acc.variance() / thermostat.temperature(),
                self.magnetization_acc.binder_cumulant(),
            )?;
        }
        self.thermostat = None;
        self.hamiltonian = None;
        self.energy_acc = Accumulator::new();
        self.magnetization_acc = Accumulator::new();
        Ok(())
    }

    fn after_step(&mut self, state: &State<S>) -> InstrumentResult<()> {
        if let (Some(thermostat), Some(hamiltonian)) = (&self.thermostat, &self.hamiltonian) {
            let energy = hamiltonian.total_energy(thermostat, state) / state.len() as f64;
            let magnetization = state.magnetization().magnitude() / state.len() as f64;
            self.energy_acc.collect(energy);
            self.magnetization_acc.collect(magnetization);
        }
        Ok(())
    }
}

/// An instrument that stores raw stats in a parquet file.
pub struct RawStatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    io: ObservableParquetIO,
    stage: usize,
    thermostat: Option<Thermostat>,
    hamiltonian: Option<H>,
    energy: Vec<f64>,
    magnetization: Vec<f64>,
    phantom: PhantomData<S>,
}

impl<H, S> RawStatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    pub fn try_new<P: AsRef<Path>>(path: P) -> IoResult<Self> {
        Ok(Self {
            io: ObservableParquetIO::try_new(path)?,
            stage: 0,
            thermostat: None,
            hamiltonian: None,
            energy: Vec::new(),
            magnetization: Vec::new(),
            phantom: PhantomData,
        })
    }
}

impl<H, S> Instrument<H, S> for RawStatSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    fn on_relax_start(
        &mut self,
        thermostat: &Thermostat,
        hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        self.thermostat = Some(thermostat.clone());
        self.hamiltonian = Some(hamiltonian.clone());
        self.energy.clear();
        self.magnetization.clear();
        Ok(())
    }

    fn on_relax_end(&mut self) -> InstrumentResult<()> {
        if let Some(thermostat) = &self.thermostat {
            self.io.write(
                true,
                self.stage,
                thermostat,
                &self.energy,
                &self.magnetization,
            )?;
        }
        self.stage += 1;
        self.hamiltonian = None;
        self.thermostat = None;
        self.energy.clear();
        self.magnetization.clear();
        Ok(())
    }

    fn on_measure_start(
        &mut self,
        thermostat: &Thermostat,
        hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        self.thermostat = Some(thermostat.clone());
        self.hamiltonian = Some(hamiltonian.clone());
        self.energy.clear();
        self.magnetization.clear();
        Ok(())
    }

    fn on_measure_end(&mut self) -> InstrumentResult<()> {
        if let Some(thermostat) = &self.thermostat {
            self.io.write(
                false,
                self.stage,
                thermostat,
                &self.energy,
                &self.magnetization,
            )?;
        }
        self.stage += 1;
        self.hamiltonian = None;
        self.thermostat = None;
        self.energy.clear();
        self.magnetization.clear();
        Ok(())
    }

    fn after_step(&mut self, state: &State<S>) -> InstrumentResult<()> {
        if let (Some(thermostat), Some(hamiltonian)) = (&self.thermostat, &self.hamiltonian) {
            let energy = hamiltonian.total_energy(thermostat, state) / state.len() as f64;
            let magnetization = state.magnetization().magnitude() / state.len() as f64;
            self.energy.push(energy);
            self.magnetization.push(magnetization);
        }
        Ok(())
    }
}

pub struct StateSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    io: StateParquetIO,
    frequency: usize,
    relax: Option<bool>,
    step: usize,
    stage: usize,
    phantom_h: PhantomData<H>,
    phantom_s: PhantomData<S>,
}

impl<H, S> StateSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    pub fn try_new<P: AsRef<Path>>(path: P, frequency: usize) -> IoResult<Self> {
        Ok(Self {
            io: StateParquetIO::try_new(path)?,
            frequency,
            relax: None,
            stage: 0,
            step: 0,
            phantom_h: PhantomData,
            phantom_s: PhantomData,
        })
    }
}

impl<H, S> Instrument<H, S> for StateSensor<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    fn on_relax_start(
        &mut self,
        _thermostat: &Thermostat,
        _hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        self.relax = Some(true);
        Ok(())
    }

    fn on_relax_end(&mut self) -> InstrumentResult<()> {
        self.relax = None;
        self.step = 0;
        self.stage += 1;
        Ok(())
    }

    fn on_measure_start(
        &mut self,
        _thermostat: &Thermostat,
        _hamiltonian: &H,
        _state: &State<S>,
    ) -> InstrumentResult<()> {
        self.relax = Some(false);
        Ok(())
    }

    fn on_measure_end(&mut self) -> InstrumentResult<()> {
        self.relax = None;
        self.step = 0;
        self.stage += 1;
        Ok(())
    }

    fn after_step(&mut self, state: &State<S>) -> InstrumentResult<()> {
        if self.step.is_multiple_of(self.frequency)
            && let Some(relax) = self.relax
        {
            self.io.write(relax, self.stage, self.step, state)?;
        }
        self.step += 1;
        Ok(())
    }
}
