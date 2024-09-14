//! This module contains the energy components of the system.
//!
//! Energy components are parts of the Hamiltonian that can be aggregated.

use crate::state::{Spin, State};
use sprs::{CsMat, TriMat};
use std::iter::Iterator;
use std::marker::PhantomData;
use vegas_lattice::Lattice;

/// A trait that represents an energy component of the system.
///
/// An energy component is characterized by the fact that it can
/// compute the energy of a given site for a given state.
pub trait Hamiltonian<S: Spin>: Clone {
    /// Get the energy of a given site for a state.
    ///
    /// Panics:
    ///
    /// This function will panic of the index is greater than the size
    /// of the state vector (business as usual)
    fn energy(&self, state: &State<S>, index: usize) -> f64;

    /// Compute the total energy of a state.
    fn total_energy(&self, state: &State<S>) -> f64 {
        (0..state.len())
            .map(|i| self.energy(state, i))
            .fold(0f64, |s, i| s + i)
    }
}

/// Some constant energy that doesn't depend on the state.
#[derive(Clone, Debug)]
pub struct Gauge {
    value: f64,
}

impl Gauge {
    /// Create a new gauge energy.
    pub fn new(val: f64) -> Self {
        Self { value: val }
    }
}

impl<T: Spin> Hamiltonian<T> for Gauge {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        self.value
    }
}

/// Strong preference for a given axis.
#[derive(Clone, Debug)]
pub struct UniaxialAnisotropy<S>
where
    S: Spin,
{
    reference: S,
    strength: f64,
}

impl<S> UniaxialAnisotropy<S>
where
    S: Spin,
{
    pub fn new(s: S, k: f64) -> Self {
        Self {
            reference: s,
            strength: k,
        }
    }
}

impl<S> Hamiltonian<S> for UniaxialAnisotropy<S>
where
    S: Spin,
{
    fn energy(&self, state: &State<S>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        let s = state.at(index);
        s.dot(&self.reference).powi(2) * self.strength
    }

    fn total_energy(&self, state: &State<S>) -> f64 {
        state
            .spins()
            .iter()
            .map(|s| (s.dot(&self.reference)).powi(2))
            .sum()
    }
}

/// Energy resulting from a magnetic field.
#[derive(Clone, Debug)]
pub struct ZeemanEnergy<S>
where
    S: Spin,
{
    reference: S,
    strength: f64,
}

impl<T> ZeemanEnergy<T>
where
    T: Spin,
{
    pub fn new(s: T, h: f64) -> Self {
        Self {
            reference: s,
            strength: h,
        }
    }
}

impl<S> Hamiltonian<S> for ZeemanEnergy<S>
where
    S: Spin,
{
    fn energy(&self, state: &State<S>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        let s = state.at(index);
        s.dot(&self.reference) * self.strength
    }

    fn total_energy(&self, state: &State<S>) -> f64 {
        -state
            .spins()
            .iter()
            .map(|s| s.dot(&self.reference))
            .sum::<f64>()
    }
}

/// Energy resulting from the exchange interaction.
#[derive(Clone, Debug)]
pub struct Exchange {
    exchange: CsMat<f64>,
}

impl Exchange {
    pub fn new(exc: CsMat<f64>) -> Self {
        Self { exchange: exc }
    }

    pub fn from_lattice(lattice: &Lattice) -> Self {
        let nsites = lattice.sites().len();
        let mut mat = TriMat::<f64>::new((nsites, nsites));
        for vertex in lattice.vertices() {
            if vertex.source() <= vertex.target() {
                mat.add_triplet(vertex.source(), vertex.target(), 1.0);
                mat.add_triplet(vertex.target(), vertex.source(), 1.0);
            }
        }
        let csr = mat.to_csr();
        Self::new(csr)
    }
}

impl<S> Hamiltonian<S> for Exchange
where
    S: Spin,
{
    fn energy(&self, state: &State<S>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        let site = state.at(index);
        if let Some(row) = self.exchange.outer_view(index) {
            row.iter()
                .map(|(nbi, exc)| (state.at(nbi), exc))
                .map(|(nb, exc)| -exc * (site.dot(nb)))
                .fold(0f64, |s, i| s + i)
        } else {
            // Just retun 0.0 for out of ranges.
            0.0
        }
    }

    fn total_energy(&self, state: &State<S>) -> f64 {
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
#[derive(Clone, Debug)]
pub struct Compound<S, U, V>
where
    S: Spin,
    U: Hamiltonian<S>,
    V: Hamiltonian<S>,
{
    a: U,
    b: V,
    phantom: PhantomData<S>,
}

impl<T, U, V> Compound<T, U, V>
where
    T: Spin,
    U: Hamiltonian<T>,
    V: Hamiltonian<T>,
{
    pub fn new(a: U, b: V) -> Self {
        Self {
            a,
            b,
            phantom: PhantomData,
        }
    }
}

impl<S, U, V> Hamiltonian<S> for Compound<S, U, V>
where
    S: Spin,
    U: Hamiltonian<S>,
    V: Hamiltonian<S>,
{
    fn energy(&self, state: &State<S>, index: usize) -> f64 {
        self.a.energy(state, index) + self.b.energy(state, index)
    }
}

/// A macro to easily build complex hamiltonians.
///
/// Examples:
///
/// ```
/// #[macro_use] extern crate vegas;
/// use vegas::state::{Spin, HeisenbergSpin};
/// use vegas::hamiltonian::{Gauge, UniaxialAnisotropy};
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
        $crate::hamiltonian::Compound::new($I, $J)
        );
    ($I: expr, $J: expr, $($K: expr),+) => (
        hamiltonian!(@flatten hamiltonian!($I, $J), $($K,)+)
        );
}

#[cfg(test)]
mod tests {
    use super::{Compound, Gauge, Hamiltonian, UniaxialAnisotropy, ZeemanEnergy};
    use crate::state::{HeisenbergSpin, Spin, State};

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
        let compound = Compound::new(gauge, anisotropy);
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
