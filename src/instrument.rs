//! Instrument module for hooking into the state of the simulation.

use std::{io::Write, marker::PhantomData};

use crate::{
    accumulator::Accumulator,
    error::InstrumentResult,
    hamiltonian::Hamiltonian,
    state::{Magnetization, Spin, State},
    thermostat::Thermostat,
};

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
    fn after_step(&mut self, _step: usize, _state: &State<S>) {}
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
                "{:.16} {:.16} {:.16} {:.16} {:.16} {:.16}",
                thermostat.temperature(),
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

    fn after_step(&mut self, _step: usize, state: &State<S>) {
        if let (Some(thermostat), Some(hamiltonian)) = (&self.thermostat, &self.hamiltonian) {
            let energy = hamiltonian.total_energy(thermostat, state) / state.len() as f64;
            let magnetization = state.magnetization().magnitude() / state.len() as f64;
            self.energy_acc.collect(energy);
            self.magnetization_acc.collect(magnetization);
        }
    }
}
