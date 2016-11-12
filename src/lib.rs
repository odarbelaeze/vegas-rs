//! This little library should have traits and routines to build a generic
//! atomistic Monte Carlo simulation program to deal with magnetic properties
//! of materials.

use std::ops::{Mul};


/// This trait specifies what a spin is for me.
pub trait Spin where
    for<'a, 'b> &'a Self: Mul<&'b Self, Output = f64>
{
    fn up() -> Self;
    fn down() -> Self;

    /// This method should return a spin based on an old spin reference
    ///
    /// Examples:
    ///
    /// ```rust,ignore
    /// let old = Spin::up();
    /// let new = Spin::perturbation_of(old);
    /// ```
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

    fn perturbation_of(other: &HeisenbergSpin) -> HeisenbergSpin {
        let &HeisenbergSpin(elems) = other;
        HeisenbergSpin([ - elems[0], - elems[1], - elems[2], ])
    }
}


impl<'a, 'b> Mul<&'a HeisenbergSpin> for &'b HeisenbergSpin {
    type Output = f64;

    fn mul(self, other: &'a HeisenbergSpin) -> f64 {
        let &HeisenbergSpin(_self) = self;
        let &HeisenbergSpin(_other) = other;
        _self
            .iter()
            .zip(_other.iter())
            .map(|(a, b)| a * b)
            .fold(0f64, |sum, i| sum + i)
    }
}


// struct State<T: Spin> where  for<'a, 'b> &'a T: Mul<&'b T, Output = f64> (Vec<T>);


#[cfg(test)]
mod tests {
    use super::IsingSpin;
    use super::HeisenbergSpin;

    #[test]
    fn ising_spin_multiplies_correctly() {
        assert_eq!(&IsingSpin::Up  * &IsingSpin::Up, 1.0);
        assert_eq!(&IsingSpin::Down  * &IsingSpin::Up, -1.0);
        assert_eq!(&IsingSpin::Up * &IsingSpin::Down, -1.0);
        assert_eq!(&IsingSpin::Down * &IsingSpin::Down, 1.0);
    }

    #[test]
    fn heisemberg_spin_multiplies_correctly() {
        assert_eq!(&HeisenbergSpin([1.0, 1.0, 1.0]) * &HeisenbergSpin([1.0, 1.0, 1.0]), 3.0);
    }
}
