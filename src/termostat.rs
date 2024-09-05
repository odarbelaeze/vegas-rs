//! A simple termostat implementation.

/// A termostat is a device that can heat or cool a system.
pub struct Termostat {
    temp: f64,
}

impl Termostat {
    /// Create a new termostat with a given temperature.
    pub fn new(temp: f64) -> Self {
        Self { temp }
    }

    /// Get the temperature of the termostat.
    pub fn temp(&self) -> f64 {
        self.temp
    }

    /// Heat the termostat by a given amount.
    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    /// Cool the termostat by a given amount.
    pub fn cool(&mut self, delta: f64) {
        self.heat(-delta);
    }
}
