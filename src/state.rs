//! This little library should have traits and routines to build a generic
//! atomistic Monte Carlo simulation program to deal with magnetic properties
//! of materials.

extern crate rand;

use rand::distributions::{IndependentSample, Range};


/// This trait specifies what a spin is for me.
pub trait Spin {
    /// New up an up Spin, this depends on what you're calling up.
    fn up() -> Self;

    /// New up a down Spin, this depends on what you're calling down, it should
    /// be anti parallel to `Spin::up()`.
    fn down() -> Self;

    /// New up a random spin.
    fn rand() -> Self;

    /// Interact with another spin
    fn interact(&self, other: &Self) -> f64;
}


/// This trait represents a spin which can be created as a perturbation of
/// another, useful for things like a Metropolis algorithm.
pub trait PerturbableSpin: Spin {
    /// New up a spin which is the perturbation of other.
    fn perturbation_of(other: &Self) -> Self;
}


pub enum IsingSpin {
    Up,
    Down,
}

impl Spin for IsingSpin {
    fn up() -> IsingSpin {
        IsingSpin::Up
    }

    fn down() -> IsingSpin {
        IsingSpin::Down
    }

    /// Randomly pick up or down for an Ising spin.
    fn rand() -> IsingSpin {
        let range = Range::new(0f64, 1f64);
        let mut rng = rand::thread_rng();
        let r = range.ind_sample(&mut rng);
        if r < 0.5f64 {
            IsingSpin::Up
        } else {
            IsingSpin::Down
        }
    }

    fn interact(&self, other: &IsingSpin) -> f64 {
        use self::IsingSpin::{Up, Down};
        match (self, other) {
            (&Up, &Up) | (&Down, &Down) => 1f64,
            _ => -1f64,
        }
    }
}

impl PerturbableSpin for IsingSpin {
    fn perturbation_of(other: &IsingSpin) -> IsingSpin {
        use self::IsingSpin::{Up, Down};
        match *other {
            Up => Down,
            Down => Up,
        }
    }
}


pub struct HeisenbergSpin([f64; 3]);

impl Spin for HeisenbergSpin {
    fn up() -> HeisenbergSpin {
        HeisenbergSpin([0f64, 0f64, 1f64])
    }

    fn down() -> HeisenbergSpin {
        HeisenbergSpin([0f64, 0f64, -1f64])
    }

    /// Gerate a random Heisenberg spin using the Marsaglia method for sphere
    /// point picking.
    fn rand() -> HeisenbergSpin {
        loop {
            let (a, b) = rand::random::<(f64, f64)>();
            let sum = a * a + b * b;
            if sum >= 1f64 {
                continue;
            }
            let dif = (1f64 - a * a - b * b).sqrt();
            return HeisenbergSpin([2f64 * a * dif, 2f64 * b * dif, 1f64 - 2f64 * sum]);
        }
    }

    fn interact(&self, other: &HeisenbergSpin) -> f64 {
        let &HeisenbergSpin(_self) = self;
        let &HeisenbergSpin(_other) = other;
        _self.iter()
            .zip(_other.iter())
            .map(|(a, b)| a * b)
            .fold(0f64, |sum, i| sum + i)
    }
}

impl PerturbableSpin for HeisenbergSpin {
    fn perturbation_of(_: &HeisenbergSpin) -> HeisenbergSpin {
        HeisenbergSpin::rand()
    }
}


pub struct State<T: Spin>(Vec<T>);

impl<T: Spin> State<T> {
    pub fn down_with_size(n: usize) -> State<T> {
        State::<T>((0..n).map(|_| T::down()).collect())
    }

    pub fn up_with_size(n: usize) -> State<T> {
        State::<T>((0..n).map(|_| T::up()).collect())
    }

    pub fn rand_with_size(n: usize) -> State<T> {
        State::<T>((0..n).map(|_| T::rand()).collect())
    }

    pub fn spins(&self) -> &Vec<T> {
        let &State::<T>(ref items) = self;
        items
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
}


#[cfg(test)]
mod tests {
    use super::{Spin, PerturbableSpin};
    use super::IsingSpin;
    use super::HeisenbergSpin;
    use super::State;

    fn real_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-15);
    }

    #[test]
    fn ising_spin_multiplies_correctly() {
        let up = IsingSpin::up();
        let down = IsingSpin::down();
        let rand = IsingSpin::rand();
        real_close(up.interact(&up), 1.0);
        real_close(up.interact(&down), -1.0);
        real_close(down.interact(&up), -1.0);
        real_close(down.interact(&down), 1.0);
        real_close(rand.interact(&rand), 1.0);
        real_close(rand.interact(&up) + rand.interact(&down), 0.0);
    }

    #[test]
    fn heisemberg_spin_multiplies_correctly() {
        let up = HeisenbergSpin::up();
        let down = HeisenbergSpin::down();
        let rand = HeisenbergSpin::rand();
        real_close(up.interact(&up), 1.0);
        real_close(up.interact(&down), -1.0);
        real_close(rand.interact(&rand), 1.0);
    }

    #[test]
    fn heisenberg_spins_are_random() {
        let HeisenbergSpin(a) = HeisenbergSpin::rand();
        let HeisenbergSpin(b) = HeisenbergSpin::rand();
        assert!(a != b);
    }

    #[test]
    fn random_heisenberg_spins_are_unit() {
        for _ in 0..100 {
            let HeisenbergSpin(a) = HeisenbergSpin::rand();
            let norm = a.iter().map(|i| i * i).fold(0f64, |s, i| s + i);
            real_close(norm, 1.0);
        }
    }

    #[test]
    fn perturbation_of_heisenberg_spins() {
        let a = HeisenbergSpin::rand();
        let b = HeisenbergSpin::perturbation_of(&a);
        let HeisenbergSpin(aitems) = a;
        let HeisenbergSpin(bitems) = b;
        assert!(aitems != bitems);
    }

    #[test]
    fn lengths_of_states() {
        let State(items) = State::<HeisenbergSpin>::up_with_size(10);
        assert_eq!(items.len(), 10);
    }
}