//! Integrators for Monte Carlo simulations.
//!
//! This module contains various integrators that can be used to sample the
//! phase space of a system using Monte Carlo methods. It includes the
//! Metropolis integrator and a variant that flips spins instead of randomizing them.
//!
//! # Example
//!
//! ```rust
//! use vegas::{
//!     energy::{Gauge, Hamiltonian},
//!     integrator::{Integrator, MetropolisFlipIntegrator},
//!     state::{Field, Spin, IsingSpin, State},
//!     thermostat::Thermostat,
//! };
//! use rand::thread_rng;
//!
//! // Define a Hamiltonian (e.g., Gauge Hamiltonian).
//! let hamiltonian = Gauge::new(1.0);
//! let thermostat = Thermostat::new(2.5, Field::zero());
//! let integrator = MetropolisFlipIntegrator::new();
//! let mut rng = thread_rng();
//! let state: State<IsingSpin> = State::rand_with_size(&mut rng, 100);
//! let new_state = integrator.step(&mut rng, &thermostat, &hamiltonian, state);
//! ```

use std::collections::VecDeque;

use crate::{
    energy::Hamiltonian,
    state::{Flip, IsingSpin, Spin, State},
    thermostat::Thermostat,
};
use rand::Rng;
use rand::distr::{Distribution, Uniform};
use vegas_lattice::Lattice;

/// An integrator is a method that allows you to sample the phase space of a
/// system.
pub trait Integrator<S: Spin> {
    /// Perform a single step of the integrator.
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat<S>,
        hamiltonian: &H,
        state: State<S>,
    ) -> State<S>;
}

/// The most common integrator is the Metropolis integrator.
///
/// The Metropolis integrator is a Monte Carlo method that allows you to sample
/// the phase space of a system. It is based on the Metropolis algorithm, which
/// is a Markov Chain Monte Carlo method.
#[derive(Debug, Default)]
pub struct MetropolisIntegrator {}

impl MetropolisIntegrator {
    /// Create a new Metropolis integrator with a given temperature.
    pub fn new() -> Self {
        Self {}
    }
}

impl<S: Spin> Integrator<S> for MetropolisIntegrator {
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat<S>,
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
/// is a Markov Chain Monte Carlo method.
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
    fn step<R: Rng, H: Hamiltonian<S>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat<S>,
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

/// Wolff cluster integrator for Ising spins.
///
/// The Wolff integrator is a Monte Carlo method that allows you to sample
/// the phase space of a system using cluster updates. It is based on the
/// Wolff algorithm, which is a cluster Monte Carlo method.
#[derive(Debug)]
pub struct WolffIntegrator {
    neighbor_list: Vec<Vec<usize>>,
}

impl WolffIntegrator {
    /// Create a new Wolff integrator with a given neighbor list.
    pub fn new(neighbor_list: Vec<Vec<usize>>) -> Self {
        Self { neighbor_list }
    }

    /// Create a new Wolff integrator from a lattice.
    pub fn from_lattice(lattice: &Lattice) -> Self {
        let mut neighbor_list = vec![Vec::new(); lattice.sites().len()];
        for vertex in lattice.vertices() {
            neighbor_list[vertex.source()].push(vertex.target());
            neighbor_list[vertex.target()].push(vertex.source());
        }
        Self { neighbor_list }
    }
}

impl Integrator<IsingSpin> for WolffIntegrator {
    /// Perform a single step of the Wolff integrator.
    ///
    /// Even though the Hamiltonian is not used in this integrator, it is included
    /// in the function signature to comply with the `Integrator` trait. This method
    /// is only valid for Ising spins and the Exchange Hamiltonian.
    fn step<R: Rng, H: Hamiltonian<IsingSpin>>(
        &self,
        rng: &mut R,
        thermostat: &Thermostat<IsingSpin>,
        _hamiltonian: &H,
        state: State<IsingSpin>,
    ) -> State<IsingSpin> {
        // Make sure the neighbor list matches the state size
        debug_assert!(state.len() == self.neighbor_list.len());

        // Choose a random site to start the cluster
        let sites = Uniform::new(0, state.len()).expect("should always be able to create");
        let source = sites.sample(rng);

        // Build the cluster using a queue
        let mut queue = VecDeque::new();
        queue.push_back(source);
        let mut visited = vec![false; state.len()];
        let mut cluster = vec![];
        while let Some(site) = queue.pop_front() {
            if visited[site] {
                continue;
            }
            visited[site] = true;
            cluster.push(site);
            let spin = state.at(site);
            for &neighbor in &self.neighbor_list[site] {
                if !visited[neighbor] && state.at(neighbor) == spin {
                    let prob = 1.0 - (-2.0 / thermostat.temperature()).exp();
                    if rng.random::<f64>() < prob {
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        // Flip the spins in the cluster
        let mut state = state;
        for &site in &cluster {
            let old_spin = state.at(site).clone();
            state.set_at(site, old_spin.flip());
        }

        state
    }
}
