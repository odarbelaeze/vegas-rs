extern crate rand;
use rand::distributions::{IndependentSample, Range};

use state::{Spin, State};
use energy::EnergyComponent;


pub trait Integrator<S: Spin, T: EnergyComponent<S>> {
    fn step(&self, energy: &T, state: &State<S>) -> State<S>;
}


pub struct MetropolisIntegrator {
    temp: f64,
}


impl MetropolisIntegrator {
    pub fn new(temp: f64) -> MetropolisIntegrator {
        MetropolisIntegrator { temp: temp }
    }

    pub fn temp(&self) -> f64 {
        self.temp
    }

    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    pub fn cool(&mut self, delta: f64) {
        self.heat( - delta);
    }
}


impl<S, T> Integrator<S, T> for MetropolisIntegrator where
    S: Spin + Clone,
    T: EnergyComponent<S>
{
    fn step(&self, energy: &T, state: &State<S>) -> State<S> {
        let mut new_state = (*state).clone();
        let sites = Range::new(0, new_state.len());
        let mut rng = rand::thread_rng();
        for _ in 0..new_state.len() {
            let site = sites.ind_sample(&mut rng);
            let old_energy = energy.energy(&new_state, site);
            new_state.set_at(site, Spin::rand());
            let new_energy = energy.energy(&new_state, site);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue
            }
            if rand::random::<f64>() < (- delta / self.temp).exp() {
                continue
            }
            new_state.set_at(site, state.at(site).clone());
        }
        new_state
    }
}