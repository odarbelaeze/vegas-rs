//! Defines spins, magnetizations, and states.
//!
//! Spins are the basic building blocks of the library. They represent the spin
//! of an atom in a magnetic material. The library provides a `Spin` trait that
//! you can implement for your own spin types.
//!
//! # Examples
//!
//! ```rust
//! use rand::SeedableRng;
//! use rand_pcg::Pcg64;
//! use vegas::state::{HeisenbergSpin, IsingSpin, State, Spin};
//! let mut rng = Pcg64::from_rng(&mut rand::rng());
//!
//! let ising_state = State::<IsingSpin>::rand_with_size(&mut rng, 100);
//! let heisenberg_state = State::<HeisenbergSpin>::rand_with_size(&mut rng, 100);
//! let ising_spin = IsingSpin::rand(&mut rng);
//! let heisenberg_spin = HeisenbergSpin::up();
//! ```

use crate::util::marsaglia;
use rand::{
    Rng,
    distr::{Distribution, Uniform},
};
use std::iter::Sum;

/// This trait specifies what a spin is.
pub trait Spin: Clone {
    /// New up an up Spin, this depends on what you're calling up.
    fn up() -> Self;

    /// New up a down Spin, this depends on what you're calling down, it should
    /// be anti parallel to `Spin::up()`.
    fn down() -> Self;

    /// New up a random spin.
    fn rand<R: Rng>(rng: &mut R) -> Self;

    /// Create a spin from its projections along the x, y, and z axes.
    fn from_projections(sx: f64, sy: f64, sz: f64) -> Field<Self>;

    /// Dot product of two spins.
    fn dot(&self, other: &Self) -> f64;

    /// Projection of the spin along the x-axis.
    fn sx(&self) -> f64;

    /// Projection of the spin along the y-axis.
    fn sy(&self) -> f64;

    /// Projection of the spin along the z-axis.
    fn sz(&self) -> f64;
}

/// This trait represents a spin which can be flipped.
pub trait Flip {
    /// Flip the spin.
    fn flip(&self) -> Self;
}

/// This enum represents an Ising spin.
#[derive(Debug, Clone, PartialEq)]
pub enum IsingSpin {
    Up,
    Down,
}

impl Spin for IsingSpin {
    fn up() -> Self {
        IsingSpin::Up
    }

    fn down() -> Self {
        IsingSpin::Down
    }

    /// Randomly pick up or down for an Ising spin.
    fn rand<R: Rng>(rng: &mut R) -> Self {
        let range = Uniform::new(0f64, 1f64).expect("should always be able to create");
        let r = range.sample(rng);
        if r < 0.5f64 {
            IsingSpin::Up
        } else {
            IsingSpin::Down
        }
    }

    fn from_projections(_sx: f64, _sy: f64, sz: f64) -> Field<Self> {
        if sz >= 0f64 {
            Field::new(IsingSpin::Up, sz.abs())
        } else {
            Field::new(IsingSpin::Down, sz.abs())
        }
    }

    #[inline]
    fn dot(&self, other: &Self) -> f64 {
        use self::IsingSpin::{Down, Up};
        match (self, other) {
            (&Up, &Up) | (&Down, &Down) => 1f64,
            _ => -1f64,
        }
    }

    #[inline]
    fn sx(&self) -> f64 {
        0.0
    }

    #[inline]
    fn sy(&self) -> f64 {
        0.0
    }

    #[inline]
    fn sz(&self) -> f64 {
        use self::IsingSpin::{Down, Up};
        match self {
            Up => 1f64,
            Down => -1f64,
        }
    }
}

impl Flip for IsingSpin {
    fn flip(&self) -> Self {
        use self::IsingSpin::{Down, Up};
        match self {
            Up => Down,
            Down => Up,
        }
    }
}

/// Heisenberg spin.
#[derive(Debug, Clone, PartialEq)]
pub struct HeisenbergSpin([f64; 3]);

impl Spin for HeisenbergSpin {
    fn up() -> Self {
        HeisenbergSpin([0f64, 0f64, 1f64])
    }

    fn down() -> Self {
        HeisenbergSpin([0f64, 0f64, -1f64])
    }

    fn rand<R: Rng>(rng: &mut R) -> Self {
        let (x, y, z) = marsaglia(rng);
        HeisenbergSpin([x, y, z])
    }

    fn from_projections(sx: f64, sy: f64, sz: f64) -> Field<Self> {
        let magnitude = (sx * sx + sy * sy + sz * sz).sqrt();
        if magnitude.abs() < f64::EPSILON {
            Field::zero()
        } else {
            Field::new(
                HeisenbergSpin([sx / magnitude, sy / magnitude, sz / magnitude]),
                magnitude,
            )
        }
    }

    #[inline]
    fn dot(&self, other: &Self) -> f64 {
        let &HeisenbergSpin(_self) = self;
        let &HeisenbergSpin(_other) = other;
        _self
            .iter()
            .zip(_other.iter())
            .map(|(a, b)| a * b)
            .fold(0f64, |sum, i| sum + i)
    }

    #[inline]
    fn sx(&self) -> f64 {
        self.0[0]
    }

    #[inline]
    fn sy(&self) -> f64 {
        self.0[1]
    }

    #[inline]
    fn sz(&self) -> f64 {
        self.0[2]
    }
}

/// Field represents a magnetic field for the given spin type.
#[derive(Debug, Clone)]
pub struct Field<S: Spin> {
    orientation: S,
    magnitude: f64,
}

impl<S: Spin> Field<S> {
    /// Create a new field with given direction and magnitude.
    pub fn new(orientation: S, magnitude: f64) -> Self {
        Field {
            orientation,
            magnitude,
        }
    }

    /// Create a zero field.
    pub fn zero() -> Self {
        Field {
            orientation: S::up(),
            magnitude: 0.0,
        }
    }

    /// Get the magnitude of the field.
    pub fn magnitude(&self) -> f64 {
        self.magnitude.abs()
    }

    /// Get the orientation of the field.
    pub fn orientation(&self) -> &S {
        &self.orientation
    }
}

impl<S: Spin> Default for Field<S> {
    fn default() -> Self {
        Field::zero()
    }
}

impl<S: Spin> Sum<S> for Field<S> {
    fn sum<I: Iterator<Item = S>>(iter: I) -> Self {
        let (px, py, pz) = iter.fold((0f64, 0f64, 0f64), |(accx, accy, accz), i| {
            (accx + i.sx(), accy + i.sy(), accz + i.sz())
        });
        S::from_projections(px, py, pz)
    }
}

/// A state of spins.
#[derive(Clone)]
pub struct State<S: Spin>(Vec<S>);

impl<S: Spin> State<S> {
    /// Create a new state with a given number of spins pointing down.
    pub fn down_with_size(n: usize) -> Self {
        State::<S>((0..n).map(|_| S::down()).collect())
    }

    /// Create a new state with a given number of spins pointing up.
    pub fn up_with_size(n: usize) -> Self {
        State::<S>((0..n).map(|_| S::up()).collect())
    }

    /// Create a new state with a given number of random spins.
    pub fn rand_with_size<R: Rng>(rng: &mut R, n: usize) -> Self {
        State::<S>((0..n).map(|_| S::rand(rng)).collect())
    }

    /// View the spins in a state.
    pub fn spins(&self) -> &Vec<S> {
        let State::<S>(items) = self;
        items
    }

    /// View a particular spin in a state.
    pub fn at(&self, index: usize) -> &S {
        &self.spins()[index]
    }

    /// Set a particular spin in a state.
    pub fn set_at(&mut self, index: usize, spin: S) {
        self.0[index] = spin;
    }

    /// Get the length of a state.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if a state is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Get the magnetization of a state.
    pub fn magnetization(&self) -> Field<S>
    where
        S: Spin,
    {
        self.spins().iter().cloned().sum()
    }
}

#[cfg(test)]
mod tests {
    use crate::state::{HeisenbergSpin, IsingSpin, Spin, State};
    use rand::SeedableRng;
    use rand_pcg::Pcg64;

    fn assert_real_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-15);
    }

    #[test]
    fn ising_spin_multiplies_correctly() {
        let up = IsingSpin::up();
        let down = IsingSpin::down();
        let mut rng = Pcg64::from_rng(&mut rand::rng());
        let rand = IsingSpin::rand(&mut rng);
        assert_real_close(up.dot(&up), 1.0);
        assert_real_close(up.dot(&down), -1.0);
        assert_real_close(down.dot(&up), -1.0);
        assert_real_close(down.dot(&down), 1.0);
        assert_real_close(rand.dot(&rand), 1.0);
        assert_real_close(rand.dot(&up) + rand.dot(&down), 0.0);
    }

    #[test]
    fn ising_magnetization_can_be_added_with_spin() {
        let state = State::<IsingSpin>::up_with_size(10);
        let mag = state.magnetization();
        assert_eq!(mag.magnitude(), 10.0);
        assert_eq!(mag.orientation(), &IsingSpin::up());
    }

    #[test]
    fn heisemberg_spin_multiplies_correctly() {
        let up = HeisenbergSpin::up();
        let down = HeisenbergSpin::down();
        let rand = HeisenbergSpin::rand(&mut rand::rng());
        assert_real_close(up.dot(&up), 1.0);
        assert_real_close(up.dot(&down), -1.0);
        assert_real_close(rand.dot(&rand), 1.0);
    }

    #[test]
    fn heisenberg_spins_are_random() {
        let HeisenbergSpin(a) = HeisenbergSpin::rand(&mut rand::rng());
        let HeisenbergSpin(b) = HeisenbergSpin::rand(&mut rand::rng());
        assert!(a != b);
    }

    #[test]
    fn random_heisenberg_spins_are_unit() {
        for _ in 0..100 {
            let HeisenbergSpin(a) = HeisenbergSpin::rand(&mut rand::rng());
            let norm = a.iter().map(|i| i * i).fold(0f64, |s, i| s + i);
            assert_real_close(norm, 1.0);
        }
    }

    #[test]
    fn lengths_of_states() {
        let State(items) = State::<HeisenbergSpin>::up_with_size(10);
        assert_eq!(items.len(), 10);
    }
}
