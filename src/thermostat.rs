//! A simple thermostat implementation.

/// A thermostat is a device that can heat or cool a system.
pub struct Thermostat {
    temp: f64,
}

impl Thermostat {
    /// Create a new thermostat with a given temperature.
    pub fn new(temp: f64) -> Self {
        Self { temp }
    }

    /// Get the temperature of the thermostat.
    pub fn temp(&self) -> f64 {
        self.temp
    }

    /// Heat the thermostat by a given amount.
    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    /// Cool the thermostat by a given amount.
    pub fn cool(&mut self, delta: f64) {
        self.heat(-delta);
    }
}
