//! Represents a thermal bath for spin systems.
//!
//! A thermostat is characterized by its temperature and external magnetic field.
//!
//! # Examples
//!
//! ```rust
//! use vegas::thermostat::Thermostat;
//!
//! // Create a thermostat with a reduced temperature of 3 and an external field of 0.5.
//! let thermostat = Thermostat::new(3.0, 0.5);
//! ```

/// A thermostat representing a thermal bath with a given temperature and field.
#[derive(Debug, Clone)]
pub struct Thermostat {
    temperature: f64,
    field: f64,
}

impl Thermostat {
    /// Create a new thermostat with a given temperature and field.
    pub fn new(temperature: f64, field: f64) -> Self {
        let _temp = if temperature < f64::EPSILON {
            f64::EPSILON
        } else {
            temperature
        };
        Thermostat {
            temperature: _temp,
            field,
        }
    }

    /// Create a new thermostat with temperature near zero and zero field.
    pub fn near_zero() -> Self {
        Thermostat {
            temperature: f64::EPSILON,
            field: 0.0,
        }
    }

    /// Returns a new thermostat with the given temperature.
    pub fn with_temperature(&self, temperature: f64) -> Self {
        let _temp = if temperature < f64::EPSILON {
            f64::EPSILON
        } else {
            temperature
        };
        Thermostat {
            temperature: _temp,
            field: self.field,
        }
    }

    /// Returns a new thermostat with the given field.
    pub fn with_field(&self, field: f64) -> Self {
        Thermostat {
            temperature: self.temperature,
            field,
        }
    }

    /// Get the temperature of the thermostat.
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Get the field of the thermostat.
    pub fn field(&self) -> f64 {
        self.field
    }
}
