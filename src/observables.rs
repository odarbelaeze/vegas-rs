use core::fmt;
use std::fmt::{Display, Formatter};

use crate::{
    energy::HamiltonianComponent,
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
        1.0 - self.sum_fourth / self.count as f64 / 3.0 * (self.sum_sq / self.count as f64).powi(2)
    }
}

impl Default for Accumulator {
    fn default() -> Self {
        Self::new()
    }
}

pub struct Sensor {
    energy: Accumulator,
    magnetization: Accumulator,
    beta: f64,
}

impl Sensor {
    pub fn new(beta: f64) -> Sensor {
        Sensor {
            energy: Accumulator::new(),
            magnetization: Accumulator::new(),
            beta,
        }
    }
    pub fn observe<H, S>(&mut self, hamiltonian: &H, state: &State<S>)
    where
        H: HamiltonianComponent<S>,
        S: Spin,
    {
        self.energy.add(hamiltonian.total_energy(state));
        self.magnetization
            .add(state.magnetization().magnitude() / state.len() as f64);
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn energy(&self) -> f64 {
        self.energy.mean()
    }

    pub fn specific_heat(&self) -> f64 {
        self.energy.variance() / self.beta
    }

    pub fn magnetization(&self) -> f64 {
        self.magnetization.mean()
    }

    pub fn susceptibility(&self) -> f64 {
        self.magnetization.variance() / self.beta
    }

    pub fn binder_cumulant(&self) -> f64 {
        self.magnetization.binder_cumulant()
    }
}

impl Display for Sensor {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.8} {:.8} {:.8} {:.8} {:.8} {:.8}",
            self.beta(),
            self.energy(),
            self.specific_heat(),
            self.magnetization(),
            self.susceptibility(),
            self.binder_cumulant()
        )
    }
}
