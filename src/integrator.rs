//! Integrators for Monte Carlo simulations.

use rand::distributions::{Distribution, Uniform};
use rand::Rng;

use crate::energy::HamiltonianComponent;
use crate::state::{Spin, State};
use crate::termostat::Termostat;

/// An integrator is a method that allows you to sample the phase space of a
/// system.
pub trait Integrator<S: Spin, T: HamiltonianComponent<S>> {
    /// Perform a single step of the integrator.
    fn step<R: Rng>(
        &self,
        energy: &T,
        state: State<S>,
        termostat: &Termostat,
        rng: &mut R,
    ) -> State<S>;
}

/// The most common integrator is the Metropolis integrator.
///
/// The Metropolis integrator is a Monte Carlo method that allows you to sample
/// the phase space of a system. It is based on the Metropolis algorithm, which
/// is a Markov chain Monte Carlo method.
#[derive(Debug, Default)]
pub struct MetropolisIntegrator {}

impl MetropolisIntegrator {
    /// Create a new Metropolis integrator with a given temperature.
    pub fn new() -> Self {
        Self {}
    }
}

impl<S, H> Integrator<S, H> for MetropolisIntegrator
where
    S: Spin + Clone,
    H: HamiltonianComponent<S>,
{
    fn step<R: Rng>(
        &self,
        hamiltonian: &H,
        mut state: State<S>,
        termostat: &Termostat,
        rng: &mut R,
    ) -> State<S> {
        let sites = Uniform::new(0, state.len());
        for _ in 0..state.len() {
            let site = sites.sample(rng);
            let old_energy = hamiltonian.energy(&state, site);
            let old_spin = state.at(site).clone();
            state.set_at(site, Spin::rand(rng));
            let new_energy = hamiltonian.energy(&state, site);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue;
            }
            if rng.gen::<f64>() < (-delta / termostat.temp()).exp() {
                continue;
            }
            state.set_at(site, old_spin);
        }
        state
    }
}
