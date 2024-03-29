//! This module contains the energy components of the system.
//!
//! Energy components are parts of the Hamiltonian that can be aggregated.

use sprs::CsMat;
use state::{Spin, State};
use std::iter::Iterator;
use std::marker::PhantomData;

/// A trait that represents an energy component of the system.
///
/// An energy component is characterized by the fact that it can
/// compute the energy of a given site for a given state.
pub trait EnergyComponent<T: Spin> {
    /// Get the energy of a given site for a state.
    ///
    /// Panics:
    ///
    /// This function will panic of the index is greater than the size
    /// of the state vector (business as usual)
    fn energy(&self, state: &State<T>, index: usize) -> f64;

    /// Compute the total energy of a state.
    fn total_energy(&self, state: &State<T>) -> f64 {
        (0..state.len())
            .map(|i| self.energy(state, i))
            .fold(0f64, |s, i| s + i)
    }
}

/// Some constant energy that doesn't depend on the state.
pub struct Gauge {
    value: f64,
}

impl Gauge {
    /// Create a new gauge energy.
    pub fn new(val: f64) -> Self {
        Self { value: val }
    }
}

impl<T: Spin> EnergyComponent<T> for Gauge {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        self.value
    }
}

/// Strong preference for a given axis.
pub struct UniaxialAnisotropy<T: Spin> {
    reference: T,
    strength: f64,
}

impl<T: Spin> UniaxialAnisotropy<T> {
    pub fn new(s: T, k: f64) -> Self {
        Self {
            reference: s,
            strength: k,
        }
    }
}

impl<T: Spin> EnergyComponent<T> for UniaxialAnisotropy<T> {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        state.spins()[index].interact(&self.reference).powi(2) * self.strength
    }

    fn total_energy(&self, state: &State<T>) -> f64 {
        state
            .spins()
            .iter()
            .map(|s| s.interact(&self.reference).powi(2))
            .fold(0f64, |s, i| s + i)
    }
}

/// Energy resulting from a magnetic field.
pub struct ZeemanEnergy<T: Spin> {
    reference: T,
    strength: f64,
}

impl<T: Spin> ZeemanEnergy<T> {
    pub fn new(s: T, h: f64) -> Self {
        Self {
            reference: s,
            strength: h,
        }
    }
}

impl<T: Spin> EnergyComponent<T> for ZeemanEnergy<T> {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        -state.spins()[index].interact(&self.reference) * self.strength
    }

    fn total_energy(&self, state: &State<T>) -> f64 {
        -state
            .spins()
            .iter()
            .map(|s| s.interact(&self.reference))
            .fold(0f64, |s, i| s + i)
    }
}

/// Energy resulting from the exchange interaction.
pub struct ExchangeEnergy {
    exchange: CsMat<f64>,
}

impl ExchangeEnergy {
    pub fn new(exc: CsMat<f64>) -> Self {
        Self { exchange: exc }
    }
}

impl<T: Spin> EnergyComponent<T> for ExchangeEnergy {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        let site = state.at(index);
        if let Some(row) = self.exchange.outer_view(index) {
            row.iter()
                .map(|(nbi, exc)| (state.at(nbi), exc))
                .map(|(nb, exc)| -exc * site.interact(nb))
                .fold(0f64, |s, i| s + i)
        } else {
            // Just retun 0.0 for out of ranges.
            0.0
        }
    }

    fn total_energy(&self, state: &State<T>) -> f64 {
        (0..state.len())
            .map(|i| self.energy(state, i))
            .fold(0f64, |s, i| s + i)
            / 2.0
    }
}

/// A compound energy is the sum of two energy components.
///
/// The key point here is that you one of the energy components
/// can be a compound energy itself.
pub struct CompoundEnergy<T, U, V>
where
    T: Spin,
    U: EnergyComponent<T>,
    V: EnergyComponent<T>,
{
    a: U,
    b: V,
    phantom: PhantomData<T>,
}

impl<T, U, V> CompoundEnergy<T, U, V>
where
    T: Spin,
    U: EnergyComponent<T>,
    V: EnergyComponent<T>,
{
    pub fn new(a: U, b: V) -> Self {
        Self {
            a,
            b,
            phantom: PhantomData,
        }
    }
}

impl<T, U, V> EnergyComponent<T> for CompoundEnergy<T, U, V>
where
    T: Spin,
    U: EnergyComponent<T>,
    V: EnergyComponent<T>,
{
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        self.a.energy(state, index) + self.b.energy(state, index)
    }
}

/// A macro to easily build complex hamiltonians.
///
/// Examples:
///
/// ```
/// #[macro_use] extern crate vegas_rs;
/// use vegas_rs::state::{Spin, HeisenbergSpin};
/// use vegas_rs::energy::{Gauge, UniaxialAnisotropy};
/// fn main() {
///     let _hamiltonian =  hamiltonian!(
///         UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0),
///         Gauge::new(1.0)
///     );
/// }
/// ```
#[macro_export]
macro_rules! hamiltonian {
    (@flatten $I: expr,) => (
        $I
        );
    (@flatten $I: expr, $J: expr, $($K:expr,)*) => (
        hamiltonian!(@flatten hamiltonian!($I, $J), $($K,)*)
        );
    ($I: expr) => (
        $I
        );
    ($I: expr, $J: expr) => (
        $crate::energy::CompoundEnergy::new($I, $J)
        );
    ($I: expr, $J: expr, $($K: expr),+) => (
        hamiltonian!(@flatten hamiltonian!($I, $J), $($K,)+)
        );
}

#[cfg(test)]
mod tests {
    use super::{CompoundEnergy, EnergyComponent, Gauge, UniaxialAnisotropy, ZeemanEnergy};
    use state::{HeisenbergSpin, Spin, State};

    #[test]
    fn test_gauge_energy() {
        let ups = State::<HeisenbergSpin>::up_with_size(10);
        let gauge = Gauge::new(10.0);
        assert!(gauge.total_energy(&ups) - 100.0 < 1e-12)
    }

    #[test]
    fn test_anisotropy_energy() {
        let ups = State::<HeisenbergSpin>::up_with_size(10);
        let downs = State::<HeisenbergSpin>::down_with_size(10);
        let anisotropy = UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0);
        assert!(anisotropy.total_energy(&ups) - 100.0 < 1e-12);
        assert!(anisotropy.total_energy(&downs) - 100.0 < 1e-12)
    }

    #[test]
    fn test_zeeman_energy() {
        let ups = State::<HeisenbergSpin>::up_with_size(10);
        let downs = State::<HeisenbergSpin>::down_with_size(10);
        let anisotropy = ZeemanEnergy::new(HeisenbergSpin::up(), 1.0);
        assert!(anisotropy.total_energy(&ups) + 10.0 < 1e-12);
        assert!(anisotropy.total_energy(&downs) - 10.0 < 1e-12)
    }

    #[test]
    fn lets_try_a_simple_composition() {
        let ups = State::<HeisenbergSpin>::up_with_size(10);
        let gauge = Gauge::new(10.0);
        let anisotropy = UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0);
        let compound = CompoundEnergy::new(gauge, anisotropy);
        assert!(compound.total_energy(&ups) - 200.0 < 1e-12);
    }

    #[test]
    fn lets_try_with_a_macro() {
        let state = State::<HeisenbergSpin>::up_with_size(10);
        let hamiltonian = hamiltonian!(
            UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0),
            Gauge::new(1.0)
        );
        assert!(hamiltonian.total_energy(&state) - 200.0 < 1e-12);
    }
}
