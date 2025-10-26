//! A machine to measure samples.
//!
//! A machine contains a sample with a given temperature and field, and allows
//! to relax or measure the sample using a given integrator and hamiltonian.
//! Instruments can be attached to the machine to measure various properties
//! of the sample during the simulation.
//!
//! # Example
//!
//! ```rust
//! use vegas::{
//!     energy::{Gauge, Hamiltonian},
//!     instrument::{Instrument, StatSensor},
//!     integrator::MetropolisIntegrator,
//!     machine::Machine,
//!     state::{IsingSpin, State},
//!     thermostat::Thermostat,
//! };
//! use rand::thread_rng;
//!
//! // Define a Hamiltonian (e.g., Gauge Hamiltonian).
//! let hamiltonian = Gauge::new(1.0);
//! let thermostat = Thermostat::new(2.5, 0.0);
//! let integrator = MetropolisIntegrator::new();
//! let mut rng = thread_rng();
//! let state: State<IsingSpin> = State::rand_with_size(&mut rng, 100);
//! let instruments = Vec::new();
//! let mut machine = Machine::new(thermostat, hamiltonian, integrator, instruments, state);
//! machine.relax_for(&mut rng, 10).unwrap();
//! machine.measure_for(&mut rng, 10).unwrap();
//! ```

use crate::{
    energy::Hamiltonian,
    error::MachineResult,
    instrument::Instrument,
    integrator::Integrator,
    state::{Spin, State},
    thermostat::Thermostat,
};
use rand::Rng;

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
    pub fn relax_for<R: Rng>(&mut self, rng: &mut R, steps: usize) -> MachineResult<()> {
        for instrument in self.instruments.iter_mut() {
            instrument.on_relax_start(&self.thermostat, &self.hamiltonian, &self.state)?;
        }
        self.run(rng, steps)?;
        for instrument in self.instruments.iter_mut() {
            instrument.on_relax_end()?;
        }
        Ok(())
    }

    /// Measure the machine for a given number of steps.
    pub fn measure_for<R: Rng>(&mut self, rng: &mut R, steps: usize) -> MachineResult<()> {
        for instrument in self.instruments.iter_mut() {
            instrument.on_measure_start(&self.thermostat, &self.hamiltonian, &self.state)?;
        }
        self.run(rng, steps)?;
        for instrument in self.instruments.iter_mut() {
            instrument.on_measure_end()?;
        }
        Ok(())
    }
}
