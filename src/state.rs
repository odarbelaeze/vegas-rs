//! This little library should have traits and routines to build a generic
//! atomistic Monte Carlo simulation program to deal with magnetic properties
//! of materials.

use std::iter::Sum;
use std::ops::Add;

use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use rand_distr::Normal;

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

    #[deprecated(since = "0.0.4", note = "Use `&s1.dot(&s2)` instead")]
    /// Interact with another spin
    fn interact(&self, other: &Self) -> f64;

    fn dot(&self, other: &Self) -> f64;
}

pub trait Magnetization: Default + Clone + Add<Self, Output = Self> + Sum<Self::SpinType> {
    type SpinType: Spin;

    fn new() -> Self
    where
        Self: Sized,
    {
        Default::default()
    }
    fn magnitude(&self) -> f64;
    fn orientation(&self) -> Self::SpinType;
}

/// This trait represents a spin which can be created as a perturbation of
/// another, useful for things like a Metropolis algorithm.
pub trait PerturbableSpin: Spin {
    /// New up a spin which is the perturbation of other.
    fn perturbation_of<R: Rng>(other: &Self, rng: &mut R) -> Self;
}

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
        let range = Uniform::new(0f64, 1f64);
        let r = range.sample(rng);
        if r < 0.5f64 {
            IsingSpin::Up
        } else {
            IsingSpin::Down
        }
    }

    fn interact(&self, other: &Self) -> f64 {
        self.dot(other)
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

impl PerturbableSpin for IsingSpin {
    fn perturbation_of<T>(other: &Self, _: &mut T) -> Self {
        use self::IsingSpin::{Down, Up};
        match *other {
            Up => Down,
            Down => Up,
        }
    }
}

#[derive(Debug, Clone)]
pub struct IsingMagnetization {
    magnitude: usize,
    reference: IsingSpin,
}

impl IsingMagnetization {
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
        match (self, other) {
            (IsingSpin::Up, IsingSpin::Up) => IsingMagnetization {
                magnitude: 2,
                reference: IsingSpin::up(),
            },
            (IsingSpin::Down, IsingSpin::Down) => IsingMagnetization {
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

impl Add for IsingMagnetization {
    type Output = IsingMagnetization;

    fn add(self, other: IsingMagnetization) -> IsingMagnetization {
        match (&self.reference, &other.reference) {
            (IsingSpin::Up, IsingSpin::Up) | (IsingSpin::Down, IsingSpin::Down) => {
                IsingMagnetization {
                    magnitude: self.magnitude + other.magnitude,
                    reference: self.reference,
                }
            }
            _ => {
                if self.magnitude > other.magnitude {
                    IsingMagnetization {
                        magnitude: self.magnitude - other.magnitude,
                        reference: self.reference,
                    }
                } else {
                    IsingMagnetization {
                        magnitude: other.magnitude - self.magnitude,
                        reference: other.reference,
                    }
                }
            }
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

    /// Gerate a random Heisenberg spin using the Marsaglia method for sphere
    /// point picking.
    fn rand<T: Rng>(rng: &mut T) -> Self {
        let distribution = Normal::new(0f64, 1f64).unwrap();
        let x = distribution.sample(rng);
        let y = distribution.sample(rng);
        let z = distribution.sample(rng);
        let sum = x * x + y * y + z * z;
        if sum == 0f64 {
            return HeisenbergSpin::up();
        }
        let norm = 1f64 / sum.sqrt();
        HeisenbergSpin([x * norm, y * norm, z * norm])
    }

    fn interact(&self, other: &Self) -> f64 {
        self.dot(other)
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

impl PerturbableSpin for HeisenbergSpin {
    fn perturbation_of<R: Rng>(_: &Self, rng: &mut R) -> Self {
        Self::rand(rng)
    }
}

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

impl Add for HeisenbergMagnetization {
    type Output = HeisenbergMagnetization;

    fn add(self, other: HeisenbergMagnetization) -> HeisenbergMagnetization {
        let HeisenbergMagnetization(a) = self;
        let HeisenbergMagnetization(b) = other;
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

#[derive(Clone)]
pub struct State<T: Spin>(Vec<T>);

impl<T: Spin> State<T> {
    pub fn down_with_size(n: usize) -> Self {
        State::<T>((0..n).map(|_| T::down()).collect())
    }

    pub fn up_with_size(n: usize) -> Self {
        State::<T>((0..n).map(|_| T::up()).collect())
    }

    pub fn rand_with_size<R: Rng>(n: usize, rng: &mut R) -> Self {
        State::<T>((0..n).map(|_| T::rand(rng)).collect())
    }

    pub fn spins(&self) -> &Vec<T> {
        let State::<T>(items) = self;
        items
    }

    pub fn at(&self, index: usize) -> &T {
        &self.spins()[index]
    }

    pub fn set_at(&mut self, index: usize, spin: T) {
        self.0[index] = spin;
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn magnetization(&self) -> T::MagnetizationType
    where
        T: Spin + Add<T, Output = T::MagnetizationType> + Clone,
        T::MagnetizationType:
            Default + Add<T::MagnetizationType, Output = T::MagnetizationType> + Sum<T>,
    {
        self.spins().iter().cloned().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::HeisenbergSpin;
    use super::IsingMagnetization;
    use super::IsingSpin;
    use super::State;
    use super::{PerturbableSpin, Spin};
    use rand::thread_rng;

    fn real_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-15);
    }

    #[test]
    fn ising_spin_multiplies_correctly() {
        let up = IsingSpin::up();
        let down = IsingSpin::down();
        let rand = IsingSpin::rand(&mut thread_rng());
        real_close(up.dot(&up), 1.0);
        real_close(up.dot(&down), -1.0);
        real_close(down.dot(&up), -1.0);
        real_close(down.dot(&down), 1.0);
        real_close(rand.dot(&rand), 1.0);
        real_close(rand.dot(&up) + rand.dot(&down), 0.0);
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
        let rand = HeisenbergSpin::rand(&mut thread_rng());
        real_close(up.dot(&up), 1.0);
        real_close(up.dot(&down), -1.0);
        real_close(rand.dot(&rand), 1.0);
    }

    #[test]
    fn heisenberg_spins_are_random() {
        let HeisenbergSpin(a) = HeisenbergSpin::rand(&mut thread_rng());
        let HeisenbergSpin(b) = HeisenbergSpin::rand(&mut thread_rng());
        assert!(a != b);
    }

    #[test]
    fn random_heisenberg_spins_are_unit() {
        for _ in 0..100 {
            let HeisenbergSpin(a) = HeisenbergSpin::rand(&mut thread_rng());
            let norm = a.iter().map(|i| i * i).fold(0f64, |s, i| s + i);
            real_close(norm, 1.0);
        }
    }

    #[test]
    fn perturbation_of_heisenberg_spins() {
        let a = HeisenbergSpin::rand(&mut thread_rng());
        let b = HeisenbergSpin::perturbation_of(&a, &mut thread_rng());
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
