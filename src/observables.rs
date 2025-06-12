//! Observables module.

use core::fmt;
use std::fmt::{Display, Formatter};

use crate::{
    error::IOResult,
    hamiltonian::Hamiltonian,
    io::ParquetIO,
    state::{Magnetization, Spin, State},
};

/// An accumulator helps to compute statistical properties of a stream of measurements.
struct Accumulator {
    sum: f64,
    sum_sq: f64,
    sum_fourth: f64,
    count: usize,
}

impl Accumulator {
    /// Create an empty accumulator.
    pub fn new() -> Accumulator {
        Accumulator {
            sum: 0.0,
            sum_sq: 0.0,
            sum_fourth: 0.0,
            count: 0,
        }
    }

    /// Add a new measurement to the accumulator.
    pub fn collect(&mut self, value: f64) {
        self.sum += value;
        self.sum_sq += value * value;
        self.sum_fourth += value * value * value * value;
        self.count += 1;
    }

    /// Compute the mean of the measurements.
    pub fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }

    /// Compute the variance of the measurements.
    pub fn variance(&self) -> f64 {
        self.sum_sq / self.count as f64 - self.mean().powi(2)
    }

    /// Compute the Binder cumulant of the measurements.
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
    energy: Vec<f64>,
    magnetization: Vec<f64>,
    energy_acc: Accumulator,
    magnetization_acc: Accumulator,
    beta: f64,
}

impl Sensor {
    /// Create a new sensor with a given beta.
    pub fn new(beta: f64) -> Sensor {
        Sensor {
            energy: Vec::new(),
            magnetization: Vec::new(),
            energy_acc: Accumulator::new(),
            magnetization_acc: Accumulator::new(),
            beta,
        }
    }

    /// Observe the state of the system.
    pub fn observe<H, S>(&mut self, hamiltonian: &H, state: &State<S>)
    where
        H: Hamiltonian<S>,
        S: Spin,
    {
        let energy = hamiltonian.total_energy(state) / state.len() as f64;
        let magnetization = state.magnetization().magnitude() / state.len() as f64;
        self.energy.push(energy);
        self.magnetization.push(magnetization);
        self.energy_acc.collect(energy);
        self.magnetization_acc.collect(magnetization);
    }

    /// Get the beta of the sensor.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Get the mean energy of the sensor.
    pub fn mean_energy(&self) -> f64 {
        self.energy_acc.mean()
    }

    /// Get the specific heat of the sensor.
    pub fn specific_heat(&self) -> f64 {
        self.energy_acc.variance() / self.beta
    }

    /// Get the mean magnetization of the sensor.
    pub fn mean_magnetization(&self) -> f64 {
        self.magnetization_acc.mean()
    }

    /// Get the susceptibility of the sensor.
    pub fn susceptibility(&self) -> f64 {
        self.magnetization_acc.variance() / self.beta
    }

    /// Get the Binder cumulant for the magnetization of the sensor.
    pub fn binder_cumulant(&self) -> f64 {
        self.magnetization_acc.binder_cumulant()
    }

    /// Get the number of measurements taken by the sensor.
    pub fn len(&self) -> usize {
        self.energy.len()
    }

    /// Check if the sensor has no measurements.
    pub fn is_empty(&self) -> bool {
        self.energy.is_empty()
    }

    /// Get a reference to the energy measurements.
    pub fn energy(&self) -> &Vec<f64> {
        &self.energy
    }

    /// Get a reference to the magnetization measurements.
    pub fn magnetization(&self) -> &Vec<f64> {
        &self.magnetization
    }

    /// Write relax data to the output file.
    pub fn write_relax(&self, raw: &mut Option<ParquetIO>) -> IOResult<()> {
        if let Some(raw) = raw {
            raw.write(self, true)
        } else {
            Ok(())
        }
    }

    /// Write the sensor readings to an output file.
    pub fn write(&self, raw: &mut Option<ParquetIO>) -> IOResult<()> {
        println!("{}", self);
        if let Some(raw) = raw {
            raw.write(self, false)
        } else {
            Ok(())
        }
    }
}

impl Display for Sensor {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.16} {:.16} {:.16} {:.16} {:.16} {:.16}",
            self.beta(),
            self.mean_energy(),
            self.specific_heat(),
            self.mean_magnetization(),
            self.susceptibility(),
            self.binder_cumulant()
        )
    }
}
