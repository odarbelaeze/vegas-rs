//! Describes states for spin systems, for now only the Heisenberg-like
//! state is implemented.

use std::ops::Mul;

extern crate rand;

use rand::distributions::normal::StandardNormal;


pub trait StateConstructors {
    fn up(size: usize) -> Self;
    fn rand(size: usize) -> Self;
}


pub trait SpinConstructors {
    fn up() -> Self;
    fn rand() -> Self;
}


#[derive(Copy, Clone)]
pub struct Spin {
    x: f64,
    y: f64,
    z: f64,
}


impl SpinConstructors for Spin {
    fn up() -> Spin {
        Spin { x: 0.0f64, y: 0.0f64, z: 1.0f64,  }
    }

    fn rand() -> Spin {
        let StandardNormal(x) = rand::random();
        let StandardNormal(y) = rand::random();
        let StandardNormal(z) = rand::random();
        let norm = (x * x + y * y + z * z).sqrt();
        Spin { x: x / norm, y: y / norm, z: z / norm, }
    }
}


impl Mul for Spin {

    type Output = f64;

    fn mul(self, other: Spin) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

}


pub type State = Vec<Spin>;


impl StateConstructors for State {

    fn up(size: usize) -> State {
        vec![Spin::up(); size]
    }

    fn rand(size: usize) -> State {
        (0..size).map(|_| { Spin::rand() }).collect()
    }

}
