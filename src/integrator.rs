extern crate rand;

use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use energy::EnergyComponent;
use state::{Spin, State};

pub trait Integrator<S: Spin, T: EnergyComponent<S>> {
    fn step(&mut self, energy: &T, state: &State<S>) -> State<S>;
}

pub trait StateGenerator<S: Spin> {
    fn state(&mut self, nsites: usize) -> State<S>;
}

pub struct MetropolisIntegrator {
    rng: SmallRng,
    temp: f64,
}

impl MetropolisIntegrator {
    pub fn new(temp: f64) -> Self {
        Self {
            temp,
            rng: SmallRng::from_entropy(),
        }
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

impl<S, T> Integrator<S, T> for MetropolisIntegrator
where
    S: Spin + Clone,
    T: EnergyComponent<S>,
{
    fn step(&mut self, energy: &T, state: &State<S>) -> State<S> {
        let mut new_state = (*state).clone();
        let sites = Uniform::new(0, new_state.len());
        for _ in 0..new_state.len() {
            let site = sites.sample(&mut self.rng);
            let old_energy = energy.energy(&new_state, site);
            new_state.set_at(site, Spin::rand(&mut self.rng));
            let new_energy = energy.energy(&new_state, site);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue;
            }
            if self.rng.gen::<f64>() < (-delta / self.temp).exp() {
                continue;
            }
            new_state.set_at(site, state.at(site).clone());
        }
        new_state
    }
}

impl<S> StateGenerator<S> for MetropolisIntegrator
where
    S: Spin + Clone,
{
    fn state(&mut self, nsites: usize) -> State<S> {
        State::rand_with_size(nsites, &mut self.rng)
    }
}
