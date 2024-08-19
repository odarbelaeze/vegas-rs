//! Integrators for Monte Carlo simulations.

extern crate rand;

use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use crate::energy::HamiltonianComponent;
use crate::state::{Spin, State};

/// An integrator is a method that allows you to sample the phase space of a
/// system.
pub trait Integrator<S: Spin, T: HamiltonianComponent<S>> {
    /// Perform a single step of the integrator.
    fn step(&mut self, energy: &T, state: &State<S>) -> State<S>;
}

/// A state generator is a method that allows you to generate a random state.
pub trait StateGenerator<S: Spin> {
    /// Generate a random state.
    fn state(&mut self, nsites: usize) -> State<S>;
}

/// The most common integrator is the Metropolis integrator.
///
/// The Metropolis integrator is a Monte Carlo method that allows you to sample
/// the phase space of a system. It is based on the Metropolis algorithm, which
/// is a Markov chain Monte Carlo method.
pub struct MetropolisIntegrator {
    rng: SmallRng,
    temp: f64,
}

impl MetropolisIntegrator {
    /// Create a new Metropolis integrator with a given temperature.
    pub fn new(temp: f64) -> Self {
        Self {
            temp,
            rng: SmallRng::from_entropy(),
        }
    }

    /// Get the temperature of the integrator.
    pub fn temp(&self) -> f64 {
        self.temp
    }

    /// Increase the temperature of the integrator.
    pub fn heat(&mut self, delta: f64) {
        self.temp += delta;
    }

    /// Decrease the temperature of the integrator.
    pub fn cool(&mut self, delta: f64) {
        self.heat(-delta);
    }
}

impl<S, T> Integrator<S, T> for MetropolisIntegrator
where
    S: Spin + Clone,
    T: HamiltonianComponent<S>,
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
