//! A machine to measure samples.

use rand::Rng;

use crate::{
    error::MachineResult,
    hamiltonian::Hamiltonian,
    instrument::Instrument,
    integrator::Integrator,
    state::{Spin, State},
    thermostat::Thermostat,
};

/// A box containing the sample with a given temperature and field.
pub struct Machine<H, I, S>
where
    H: Hamiltonian<S>,
    I: Integrator<S>,
    S: Spin,
{
    thermostat: Thermostat,
    hamiltonian: H,
    integrator: I,
    instruments: Vec<Box<dyn Instrument<H, S>>>,
    state: State<S>,
}

impl<H, I, S> Machine<H, I, S>
where
    H: Hamiltonian<S>,
    I: Integrator<S>,
    S: Spin,
{
    /// Create a new machine with a given temperature, field, hamiltonian,
    pub fn new(
        thermostat: Thermostat,
        hamiltonian: H,
        integrator: I,
        instruments: Vec<Box<dyn Instrument<H, S>>>,
        state: State<S>,
    ) -> Self {
        Machine {
            thermostat,
            hamiltonian,
            integrator,
            instruments,
            state,
        }
    }

    /// Get the current thermostat of the machine.
    pub fn thermostat(&self) -> &Thermostat {
        &self.thermostat
    }

    /// Set the thermostat of the machine.
    pub fn set_thermostat(&mut self, thermostat: Thermostat) {
        self.thermostat = thermostat;
    }

    /// Run and observe the machine for a given number of steps.
    fn run<R: Rng>(&mut self, rng: &mut R, steps: usize) -> MachineResult<()> {
        for _ in 0..steps {
            self.state =
                self.integrator
                    .step(rng, &self.thermostat, &self.hamiltonian, self.state.clone());
            for instrument in self.instruments.iter_mut() {
                instrument.after_step(&self.state)?;
            }
        }
        Ok(())
    }

    /// Relax the machine for a given number of steps.
    pub fn relax_for<R: Rng>(&mut self, rng: &mut R, relax_steps: usize) -> MachineResult<()> {
        for instrument in self.instruments.iter_mut() {
            instrument.on_relax_start(&self.thermostat, &self.hamiltonian, &self.state)?;
        }
        self.run(rng, relax_steps)?;
        for instrument in self.instruments.iter_mut() {
            instrument.on_relax_end()?;
        }
        Ok(())
    }

    /// Measure the machine for a given number of steps.
    pub fn measure_for<R: Rng>(&mut self, rng: &mut R, measure_steps: usize) -> MachineResult<()> {
        for instrument in self.instruments.iter_mut() {
            instrument.on_measure_start(&self.thermostat, &self.hamiltonian, &self.state)?;
        }
        self.run(rng, measure_steps)?;
        for instrument in self.instruments.iter_mut() {
            instrument.on_measure_end()?;
        }
        Ok(())
    }
}
