//! This little library should have traits and routines to build a generic
//! atomistic Monte Carlo simulation program to deal with magnetic properties
//! of materials.

extern crate rand;

use std::ops::Mul;
use rand::distributions::{IndependentSample, Range};


/// This trait specifies what a spin is for me.
pub trait Spin
    where for<'a, 'b> &'a Self: Mul<&'b Self, Output = f64>
{
    /// New up an up Spin, this depends on what you're calling up.
    fn up() -> Self;

    /// New up a down Spin, this depends on what you're calling down, it should
    /// be anti parallel to `Spin::up()`.
    fn down() -> Self;

    /// New up a random spin.
    fn rand() -> Self;
}


/// This trait represents a spin which can be created as a perturbation of
/// another, useful for things like a Metropolis algorithm.
pub trait PerturbableSpin: Spin
    where for<'a, 'b> &'a Self: Mul<&'b Self, Output = f64>
{
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
}

impl PerturbableSpin for IsingSpin {
    fn perturbation_of(other: &IsingSpin) -> IsingSpin {
        use IsingSpin::{Up, Down};
        match *other {
            Up => Down,
            Down => Up,
        }
    }
}

impl<'a, 'b> Mul<&'a IsingSpin> for &'b IsingSpin {
    type Output = f64;

    fn mul(self, other: &'a IsingSpin) -> f64 {
        use IsingSpin::{Up, Down};
        match (self, other) {
            (&Up, &Up) | (&Down, &Down) => 1f64,
            _ => -1f64,
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
}

impl PerturbableSpin for HeisenbergSpin {
    fn perturbation_of(_: &HeisenbergSpin) -> HeisenbergSpin {
        HeisenbergSpin::rand()
    }
}

impl<'a, 'b> Mul<&'a HeisenbergSpin> for &'b HeisenbergSpin {
    type Output = f64;

    fn mul(self, other: &'a HeisenbergSpin) -> f64 {
        let &HeisenbergSpin(_self) = self;
        let &HeisenbergSpin(_other) = other;
        _self.iter()
            .zip(_other.iter())
            .map(|(a, b)| a * b)
            .fold(0f64, |sum, i| sum + i)
    }
}


pub struct State<T: Spin>(Vec<T>) where for<'b, 'a> &'a T: std::ops::Mul<&'b T, Output = f64>;


#[cfg(test)]
mod tests {
    use super::{Spin, PerturbableSpin};
    use super::IsingSpin;
    use super::HeisenbergSpin;

    #[test]
    fn ising_spin_multiplies_correctly() {
        assert_eq!(&IsingSpin::Up * &IsingSpin::Up, 1.0);
        assert_eq!(&IsingSpin::Down * &IsingSpin::Up, -1.0);
        assert_eq!(&IsingSpin::Up * &IsingSpin::Down, -1.0);
        assert_eq!(&IsingSpin::Down * &IsingSpin::Down, 1.0);
    }

    #[test]
    fn heisemberg_spin_multiplies_correctly() {
        assert_eq!(&HeisenbergSpin([1.0, 1.0, 1.0]) * &HeisenbergSpin([1.0, 1.0, 1.0]),
                   3.0);
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
            assert!((norm - 1f64).abs() < 1e-15);
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
}