//! Observables module.

use core::fmt;
use std::fmt::{Display, Formatter};

use crate::{
    hamiltonian::Hamiltonian,
    state::{Magnetization, Spin, State},
};

struct Accumulator {
    pub sum: f64,
    pub sum_sq: f64,
    pub sum_fourth: f64,
    pub count: usize,
}

impl Accumulator {
    pub fn new() -> Accumulator {
        Accumulator {
            sum: 0.0,
            sum_sq: 0.0,
            sum_fourth: 0.0,
            count: 0,
        }
    }

    pub fn add(&mut self, value: f64) {
        self.sum += value;
        self.sum_sq += value * value;
        self.sum_fourth += value * value * value * value;
        self.count += 1;
    }

    pub fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }

    pub fn variance(&self) -> f64 {
        self.sum_sq / self.count as f64 - self.mean().powi(2)
    }

    pub fn binder_cumulant(&self) -> f64 {
        1.0 - (self.sum_fourth / self.count as f64)
            / (3.0 * (self.sum_sq / self.count as f64).powi(2))
    }
}

impl Default for Accumulator {
    fn default() -> Self {
        Self::new()
    }
}

/// A sensor is a device that can measure the state of a system.
pub struct Sensor {
    energy: Accumulator,
    magnetization: Accumulator,
    beta: f64,
}

impl Sensor {
    /// Create a new sensor with a given beta.
    pub fn new(beta: f64) -> Sensor {
        Sensor {
            energy: Accumulator::new(),
            magnetization: Accumulator::new(),
            beta,
        }
    }

    /// Observe the state of the system.
    pub fn observe<H, S>(&mut self, hamiltonian: &H, state: &State<S>)
    where
        H: Hamiltonian<S>,
        S: Spin,
    {
        self.energy
            .add(hamiltonian.total_energy(state) / state.len() as f64);
        self.magnetization
            .add(state.magnetization().magnitude() / state.len() as f64);
    }

    /// Get the beta of the sensor.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Get the mean energy of the sensor.
    pub fn energy(&self) -> f64 {
        self.energy.mean()
    }

    /// Get the specific heat of the sensor.
    pub fn specific_heat(&self) -> f64 {
        self.energy.variance() / self.beta
    }

    /// Get the mean magnetization of the sensor.
    pub fn magnetization(&self) -> f64 {
        self.magnetization.mean()
    }

    /// Get the susceptibility of the sensor.
    pub fn susceptibility(&self) -> f64 {
        self.magnetization.variance() / self.beta
    }

    /// Get the Binder cumulant for the magnetization of the sensor.
    pub fn binder_cumulant(&self) -> f64 {
        self.magnetization.binder_cumulant()
    }
}

impl Display for Sensor {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.16} {:.16} {:.16} {:.16} {:.16} {:.16}",
            self.beta(),
            self.energy(),
            self.specific_heat(),
            self.magnetization(),
            self.susceptibility(),
            self.binder_cumulant()
        )
    }
}
