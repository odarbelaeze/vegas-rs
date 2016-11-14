use std::iter::Iterator;
use std::marker::PhantomData;
use state::{Spin, State};


pub trait EnergyComponent<T: Spin> {
    /// Get the energy of a given site for a state.
    ///
    /// Panics:
    ///
    /// This function will panic of the index is greater than the size
    /// of the state vector (business as usual)
    fn energy(&self, state: &State<T>, index: usize) -> f64;

    fn total_energy(&self, state: &State<T>) -> f64 {
        (0..state.len())
            .map(|i| self.energy(state, i))
            .fold(0f64, |s, i| s + i)
    }
}


pub struct Gauge {
    value: f64,
}


impl Gauge {
    pub fn new(val: f64) -> Gauge {
        Gauge { value: val }
    }
}


impl<T: Spin> EnergyComponent<T> for Gauge {
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        debug_assert!(index < state.len());
        self.value
    }
}


pub struct UniaxialAnisotropy<T: Spin> {
    reference: T,
    strength: f64,
}

impl<T: Spin> UniaxialAnisotropy<T> {
    pub fn new(s: T, k: f64) -> UniaxialAnisotropy<T> {
        UniaxialAnisotropy {
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
        state.spins()
            .iter()
            .map(|s| s.interact(&self.reference).powi(2))
            .fold(0f64, |s, i| s + i)
    }
}


pub struct CompoundEnergy<T, U, V>
    where T: Spin,
          U: EnergyComponent<T>,
          V: EnergyComponent<T>
{
    a: U,
    b: V,
    phantom: PhantomData<T>,
}

impl<T, U, V> EnergyComponent<T> for CompoundEnergy<T, U, V>
    where T: Spin,
          U: EnergyComponent<T>,
          V: EnergyComponent<T>
{
    fn energy(&self, state: &State<T>, index: usize) -> f64 {
        self.a.energy(&state, index) + self.b.energy(&state, index)
    }
}

impl<T, U, V> CompoundEnergy<T, U, V>
    where T: Spin,
          U: EnergyComponent<T>,
          V: EnergyComponent<T>
{
    pub fn new(a: U, b: V) -> CompoundEnergy<T, U, V> {
        CompoundEnergy {
            a: a,
            b: b,
            phantom: PhantomData,
        }
    }
}



#[cfg(test)]
mod tests {
    use super::{EnergyComponent, Gauge, UniaxialAnisotropy, CompoundEnergy};
    use state::{Spin, State, HeisenbergSpin};

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
    fn lets_try_a_simple_composition() {
        let ups = State::<HeisenbergSpin>::up_with_size(10);
        let gauge = Gauge::new(10.0);
        let anisotropy = UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0);
        let compound = CompoundEnergy::new(gauge, anisotropy);
        assert!(compound.total_energy(&ups) - 200.0 < 1e-12);
    }
}