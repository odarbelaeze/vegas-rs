pub struct Termostat {
    temp: f64,
}

impl Termostat {
    pub fn new(temp: f64) -> Self {
        Self { temp }
    }

    pub fn temp(&self) -> f64 {
        self.temp
    }

    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    pub fn cool(&mut self, delta: f64) {
        self.heat(-delta);
    }
}
