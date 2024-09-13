//! A machine to measure samples.

use rand::Rng;

use crate::{
    hamiltonian::{Hamiltonian, ZeemanEnergy},
    integrator::Integrator,
    observables::Sensor,
    state::{Spin, State},
};

/// A box containing the sample with a given temperature and field.
pub struct Machine<H, I, S>
where
    H: Hamiltonian<S>,
    I: Integrator<S>,
    S: Spin,
{
    temp: f64,
    field: f64,
    hamiltonian: H,
    integrator: I,
    state: State<S>,
}

impl<H, I, S> Machine<H, I, S>
where
    H: Hamiltonian<S>,
    I: Integrator<S>,
    S: Spin,
{
    /// Create a new machine with a given temperature, field, hamiltonian,
    pub fn new(temp: f64, field: f64, hamiltonian: H, integrator: I, state: State<S>) -> Self {
        Machine {
            temp,
            field,
            hamiltonian,
            integrator,
            state,
        }
    }

    /// Set the temperature of the machine.
    pub fn set_temperature(&mut self, temp: f64) {
        if temp < f64::EPSILON {
            return;
        }
        self.temp = temp;
    }

    /// Set the field of the machine.
    pub fn set_field(&mut self, field: f64) {
        self.field = field;
    }

    /// Run and observe the machine for a given number of steps.
    pub fn run<R: Rng>(&mut self, rng: &mut R, steps: usize) -> Sensor {
        let mut sensor = Sensor::new(self.temp);
        if self.field != 0.0 {
            let hamiltonian = hamiltonian!(
                self.hamiltonian.clone(),
                ZeemanEnergy::new(S::up(), self.field)
            );
            for _ in 0..steps {
                self.state = self
                    .integrator
                    .step(rng, self.temp, &hamiltonian, self.state.clone());
                sensor.observe(&hamiltonian, &self.state);
            }
            return sensor;
        }
        for _ in 0..steps {
            self.state =
                self.integrator
                    .step(rng, self.temp, &self.hamiltonian, self.state.clone());
            sensor.observe(&self.hamiltonian, &self.state);
        }
        sensor
    }

    pub fn state(&self) -> &State<S> {
        &self.state
    }
}
