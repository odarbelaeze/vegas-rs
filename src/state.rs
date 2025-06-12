//! This little library should have traits and routines to build a generic
//! atomistic Monte Carlo simulation program to deal with magnetic properties
//! of materials.

use std::iter::Sum;
use std::ops::Add;

use rand::Rng;
use rand::distr::{Distribution, Uniform};

use super::util::marsaglia;

/// This trait specifies what a spin is.
pub trait Spin: Clone + Add<Self, Output = Self::MagnetizationType> {
    type MagnetizationType: Magnetization<SpinType = Self>;

    /// New up an up Spin, this depends on what you're calling up.
    fn up() -> Self;

    /// New up a down Spin, this depends on what you're calling down, it should
    /// be anti parallel to `Spin::up()`.
    fn down() -> Self;

    /// New up a random spin.
    fn rand<T: Rng>(rng: &mut T) -> Self;

    /// Dot product of two spins.
    fn dot(&self, other: &Self) -> f64;
}

/// This trait represents a spin which can be flipped.
pub trait Flip {
    /// Flip the spin.
    fn flip(&self) -> Self;
}

/// This trait represents a magnetization.
pub trait Magnetization: Default + Clone + Sum<Self::SpinType> {
    type SpinType: Spin;

    /// Create a new magnetization.
    fn new() -> Self
    where
        Self: Sized,
    {
        Default::default()
    }
    /// Get the magnitude of the magnetization.
    fn magnitude(&self) -> f64;
    /// Get the orientation of the magnetization.
    fn orientation(&self) -> Self::SpinType;
}

/// This enum represents an Ising spin.
#[derive(Debug, Clone, PartialEq)]
pub enum IsingSpin {
    Up,
    Down,
}

impl Spin for IsingSpin {
    type MagnetizationType = IsingMagnetization;

    fn up() -> Self {
        IsingSpin::Up
    }

    fn down() -> Self {
        IsingSpin::Down
    }

    /// Randomly pick up or down for an Ising spin.
    fn rand<T: Rng>(rng: &mut T) -> Self {
        let range = Uniform::new(0f64, 1f64).expect("should always be able to create");
        let r = range.sample(rng);
        if r < 0.5f64 {
            IsingSpin::Up
        } else {
            IsingSpin::Down
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

/// Ising magnetization.
#[derive(Debug, Clone)]
pub struct IsingMagnetization {
    magnitude: usize,
    reference: IsingSpin,
}

impl IsingMagnetization {
    /// Create a new Ising magnetization.
    pub fn new() -> Self {
        IsingMagnetization {
            magnitude: 0,
            reference: IsingSpin::up(),
        }
    }
}

impl Magnetization for IsingMagnetization {
    type SpinType = IsingSpin;

    fn magnitude(&self) -> f64 {
        self.magnitude as f64
    }

    fn orientation(&self) -> Self::SpinType {
        self.reference.clone()
    }
}

impl Default for IsingMagnetization {
    fn default() -> Self {
        IsingMagnetization::new()
    }
}

impl Sum<IsingSpin> for IsingMagnetization {
    fn sum<I: Iterator<Item = IsingSpin>>(iter: I) -> Self {
        iter.fold(IsingMagnetization::new(), |acc, i| acc + i)
    }
}

impl Add for IsingSpin {
    type Output = IsingMagnetization;

    fn add(self, other: IsingSpin) -> IsingMagnetization {
        use IsingSpin::{Down, Up};
        match (self, other) {
            (Up, Up) => IsingMagnetization {
                magnitude: 2,
                reference: IsingSpin::up(),
            },
            (Down, Down) => IsingMagnetization {
                magnitude: 2,
                reference: IsingSpin::down(),
            },
            _ => IsingMagnetization {
                magnitude: 0,
                reference: IsingSpin::up(),
            },
        }
    }
}

impl Add<IsingSpin> for IsingMagnetization {
    type Output = IsingMagnetization;

    fn add(self, rhs: IsingSpin) -> Self::Output {
        if self.magnitude == 0 {
            return IsingMagnetization {
                magnitude: 1,
                reference: rhs,
            };
        }
        if self.reference == rhs {
            IsingMagnetization {
                magnitude: self.magnitude + 1,
                reference: self.reference,
            }
        } else {
            IsingMagnetization {
                magnitude: self.magnitude - 1,
                reference: self.reference,
            }
        }
    }
}

/// Heisenberg spin.
#[derive(Debug, Clone, PartialEq)]
pub struct HeisenbergSpin([f64; 3]);

impl Spin for HeisenbergSpin {
    type MagnetizationType = HeisenbergMagnetization;

    fn up() -> Self {
        HeisenbergSpin([0f64, 0f64, 1f64])
    }

    fn down() -> Self {
        HeisenbergSpin([0f64, 0f64, -1f64])
    }

    fn rand<T: Rng>(rng: &mut T) -> Self {
        let (x, y, z) = marsaglia(rng);
        HeisenbergSpin([x, y, z])
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
}

impl Flip for HeisenbergSpin {
    fn flip(&self) -> Self {
        let HeisenbergSpin(a) = self;
        HeisenbergSpin([-a[0], -a[1], -a[2]])
    }
}

/// Heisenberg magnetization.
#[derive(Debug, Clone)]
pub struct HeisenbergMagnetization([f64; 3]);

impl HeisenbergMagnetization {
    pub fn new() -> Self {
        HeisenbergMagnetization([0f64, 0f64, 0f64])
    }
}

impl Magnetization for HeisenbergMagnetization {
    type SpinType = HeisenbergSpin;

    fn magnitude(&self) -> f64 {
        let HeisenbergMagnetization(a) = self;
        a.iter().map(|i| i * i).sum::<f64>().sqrt()
    }

    fn orientation(&self) -> Self::SpinType {
        let HeisenbergMagnetization(a) = self;
        let magnitude = self.magnitude();
        HeisenbergSpin([a[0] / magnitude, a[1] / magnitude, a[2] / magnitude])
    }
}

impl Default for HeisenbergMagnetization {
    fn default() -> Self {
        HeisenbergMagnetization::new()
    }
}

impl Sum<HeisenbergSpin> for HeisenbergMagnetization {
    fn sum<I: Iterator<Item = HeisenbergSpin>>(iter: I) -> Self {
        iter.fold(HeisenbergMagnetization::new(), |acc, i| acc + i)
    }
}

impl Add for HeisenbergSpin {
    type Output = HeisenbergMagnetization;

    fn add(self, other: HeisenbergSpin) -> HeisenbergMagnetization {
        let HeisenbergSpin(a) = self;
        let HeisenbergSpin(b) = other;
        HeisenbergMagnetization([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }
}

impl Add<HeisenbergSpin> for HeisenbergMagnetization {
    type Output = HeisenbergMagnetization;

    fn add(self, rhs: HeisenbergSpin) -> Self::Output {
        let HeisenbergMagnetization(a) = self;
        let HeisenbergSpin(b) = rhs;
        HeisenbergMagnetization([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }
}

/// A state of spins.
#[derive(Clone)]
pub struct State<T: Spin>(Vec<T>);

impl<T: Spin> State<T> {
    /// Create a new state with a given number of spins pointing down.
    pub fn down_with_size(n: usize) -> Self {
        State::<T>((0..n).map(|_| T::down()).collect())
    }

    /// Create a new state with a given number of spins pointing up.
    pub fn up_with_size(n: usize) -> Self {
        State::<T>((0..n).map(|_| T::up()).collect())
    }

    /// Create a new state with a given number of random spins.
    pub fn rand_with_size<R: Rng>(rng: &mut R, n: usize) -> Self {
        State::<T>((0..n).map(|_| T::rand(rng)).collect())
    }

    /// View the spins in a state.
    pub fn spins(&self) -> &Vec<T> {
        let State::<T>(items) = self;
        items
    }

    /// View a particular spin in a state.
    pub fn at(&self, index: usize) -> &T {
        &self.spins()[index]
    }

    /// Set a particular spin in a state.
    pub fn set_at(&mut self, index: usize, spin: T) {
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
    pub fn magnetization(&self) -> T::MagnetizationType
    where
        T: Spin + Add<T, Output = T::MagnetizationType> + Clone,
        T::MagnetizationType: Default + Sum<T>,
    {
        self.spins().iter().cloned().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::HeisenbergSpin;
    use super::IsingMagnetization;
    use super::IsingSpin;
    use super::Spin;
    use super::State;
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
        let mag = IsingMagnetization::new();
        let up = IsingSpin::up();
        let result = mag + up;
        assert_eq!(result.magnitude, 1);
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
