/// An accumulator helps to compute statistical properties of a stream of measurements.
pub struct Accumulator {
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
