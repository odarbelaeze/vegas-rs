use std::ops::Mul;
use std::iter::Iterator;
use state::{Spin, State};


pub trait EnergyComponent<T: Spin>
    where for<'a, 'b> &'a T: Mul<&'b T, Output = f64>
{
    /// Get the energy of a given site for a state.
    ///
    /// Panics:
    ///
    /// This function will panic of the index is greater than the size
    /// of the state vector (business as usual)
    fn energy(&self, state: &State<T>, index: usize) -> f64;

    fn total_energy(&self, state: &State<T>) -> f64 {
        (0..state.spins().len())
            .map(|i| self.energy(state, i))
            .fold(0f64, |s, i| s + i)
    }
}
