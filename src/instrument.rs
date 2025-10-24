//! Instrument module for hooking into the state of the simulation.

use crate::{
    hamiltonian::Hamiltonian,
    state::{Spin, State},
    thermostat::Thermostat,
};

/// An instrument allows to hook into the simulation at various points.
pub trait Instrument<H, S>
where
    H: Hamiltonian<S>,
    S: Spin,
{
    /// Hook called after each integration step.
    fn after_step(&mut self, _step: usize, _state: &State<S>) {}

    /// Hook called when a relaxation starts.
    fn before_relax(&mut self, _thermostat: &Thermostat, _hamiltonian: &H, _state: &State<S>) {}

    /// Hook called when a relaxation ends.
    fn on_relax_end(&mut self, _state: &State<S>) {}

    /// Hook called when a measurement starts.
    fn on_measure_start(&mut self, _thermostat: &Thermostat, _hamiltonian: &H, _state: &State<S>) {}

    /// Hook called when a measurement ends.
    fn on_measure_end(&mut self, _state: &State<S>) {}
}
