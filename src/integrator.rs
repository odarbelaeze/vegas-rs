//! Integrators for Monte Carlo simulations.

use rand::Rng;
use rand::distr::{Distribution, Uniform};

use crate::thermostat::Thermostat;
use crate::{
    hamiltonian::Hamiltonian,
    state::{Flip, Spin, State},
};

/// An integrator is a method that allows you to sample the phase space of a
/// system.
pub trait Integrator<S: Spin> {
    /// Perform a single step of the integrator.
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat,
        hamiltonian: &H,
        state: State<S>,
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

impl<S: Spin> Integrator<S> for MetropolisIntegrator {
    /// Perform a single step of the Metropolis integrator.
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat,
        hamiltonian: &H,
        mut state: State<S>,
    ) -> State<S> {
        let distribution = Uniform::new(0, state.len()).expect("should always be able to create");
        for _ in 0..state.len() {
            let site_index = distribution.sample(rng);
            let old_energy = hamiltonian.energy(thermostat, &state, site_index);
            let old_spin = state.at(site_index).clone();
            state.set_at(site_index, Spin::rand(rng));
            let new_energy = hamiltonian.energy(thermostat, &state, site_index);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue;
            }
            if rng.random::<f64>() < (-delta / thermostat.temperature()).exp() {
                continue;
            }
            state.set_at(site_index, old_spin);
        }
        state
    }
}

/// Metropolis integrator that flips spins instead of randomizing them.
///
/// The Metropolis integrator is a Monte Carlo method that allows you to sample
/// the phase space of a system. It is based on the Metropolis algorithm, which
/// is a Markov chain Monte Carlo method.
#[derive(Debug, Default)]
pub struct MetropolisFlipIntegrator {}

impl MetropolisFlipIntegrator {
    /// Create a new Metropolis integrator with a given temperature.
    pub fn new() -> Self {
        Self {}
    }
}

impl<S> Integrator<S> for MetropolisFlipIntegrator
where
    S: Spin + Flip,
{
    /// Perform a single step of the Metropolis integrator.
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat,
        hamiltonian: &H,
        mut state: State<S>,
    ) -> State<S> {
        let sites = Uniform::new(0, state.len()).expect("should always be able to create");
        for _ in 0..state.len() {
            let site = sites.sample(rng);
            let old_energy = hamiltonian.energy(thermostat, &state, site);
            let old_spin = state.at(site).clone();
            state.set_at(site, old_spin.flip());
            let new_energy = hamiltonian.energy(thermostat, &state, site);
            let delta = new_energy - old_energy;
            if delta < 0.0 {
                continue;
            }
            if rng.random::<f64>() < (-delta / thermostat.temperature()).exp() {
                continue;
            }
            state.set_at(site, old_spin);
        }
        state
    }
}
