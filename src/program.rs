//! Programs to run on samples.

use crate::{
    energy::HamiltonianComponent,
    integrator::Integrator,
    observables::Sensor,
    state::{Spin, State},
    termostat::Termostat,
};

pub struct CurieTemp {
    max_temp: f64,
    min_temp: f64,
    cool_rate: f64,
    relax: usize,
    steps: usize,
}

impl CurieTemp {
    pub fn new(max_temp: f64, min_temp: f64, cool_rate: f64, relax: usize, steps: usize) -> Self {
        Self {
            max_temp,
            min_temp,
            cool_rate,
            relax,
            steps,
        }
    }

    pub fn with_max_temp(mut self, max_temp: f64) -> Self {
        self.max_temp = max_temp;
        self
    }

    pub fn with_min_temp(mut self, min_temp: f64) -> Self {
        self.min_temp = min_temp;
        self
    }

    pub fn with_cool_rate(mut self, cool_rate: f64) -> Self {
        self.cool_rate = cool_rate;
        self
    }

    pub fn with_relax(mut self, relax: usize) -> Self {
        self.relax = relax;
        self
    }

    pub fn with_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }

    pub fn run<I, H, S>(&self, integrator: &mut I, hamiltonian: &H, mut state: State<S>) -> State<S>
    where
        I: Integrator<S, H>,
        H: HamiltonianComponent<S>,
        S: Spin,
    {
        let mut termostat = Termostat::new(self.max_temp);
        loop {
            let mut sensor = Sensor::new(termostat.temp());
            for _ in 0..self.relax {
                state = integrator.step(hamiltonian, state, &termostat);
            }
            for _ in 0..self.steps {
                state = integrator.step(hamiltonian, state, &termostat);
                sensor.observe(hamiltonian, &state);
            }
            println!("{}", sensor);
            termostat.cool(self.cool_rate);
            if termostat.temp() < self.min_temp {
                break;
            }
        }
        state
    }
}

impl Default for CurieTemp {
    fn default() -> Self {
        Self::new(3.0, f64::EPSILON, 0.1, 1000, 20000)
    }
}
