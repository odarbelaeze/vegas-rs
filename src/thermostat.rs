//! Represents a thermal bath for spin systems.

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

    /// Set the temperature of the thermostat.
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

    /// Set the field of the thermostat.
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
