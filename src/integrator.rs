//! Integrators for Monte Carlo simulations.

use rand::distributions::{Distribution, Uniform};
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
pub struct MetropolisIntegrator<R>
where
    R: Rng,
{
    rng: R,
    temp: f64,
}

impl<R> MetropolisIntegrator<R>
where
    R: Rng + SeedableRng,
{
    /// Create a new Metropolis integrator with a given temperature.
    pub fn new(temp: f64) -> Self {
        Self {
            temp,
            rng: R::from_entropy(),
        }
    }

    /// Create a new Metropolis integrator with a given random number generator
    pub fn new_with_rng(temp: f64, rng: R) -> Self {
        Self { temp, rng }
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

impl<S, T, R> Integrator<S, T> for MetropolisIntegrator<R>
where
    S: Spin + Clone,
    T: HamiltonianComponent<S>,
    R: Rng,
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

impl<S, R> StateGenerator<S> for MetropolisIntegrator<R>
where
    S: Spin + Clone,
    R: Rng,
{
    fn state(&mut self, nsites: usize) -> State<S> {
        State::rand_with_size(nsites, &mut self.rng)
    }
}
