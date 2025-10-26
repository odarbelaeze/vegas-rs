//! Utility functions for magnetic field simulations.

use rand::{
    Rng,
    distr::{Distribution, Uniform},
};

/// Marsaglia's method for generating random points on a unit sphere.
///
/// Returns a tuple (x, y, z) representing a point in the surface of the unit sphere.
///
/// # Examples
///
/// ```rust
/// use rand::thread_rng;
/// use vegas::util::marsaglia;
/// let mut rng = thread_rng();
/// let (x, y, z) = marsaglia(&mut rng);
/// println!("Random point on unit sphere: ({}, {}, {})", x, y, z);
/// ```
pub fn marsaglia<R: Rng>(rng: &mut R) -> (f64, f64, f64) {
    loop {
        let distribution = Uniform::new(-1.0, 1.0).expect("should always be able to create");
        let x1 = distribution.sample(rng);
        let x2 = distribution.sample(rng);
        if x1 * x1 + x2 * x2 >= 1f64 {
            continue;
        }
        let x = 2f64 * x1 * (1f64 - x1 * x1 - x2 * x2).sqrt();
        let y = 2f64 * x2 * (1f64 - x1 * x1 - x2 * x2).sqrt();
        let z = 1f64 - 2f64 * (x1 * x1 + x2 * x2);
        return (x, y, z);
    }
}
