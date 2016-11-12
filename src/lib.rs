use std::ops::{Mul};


trait Spin:
// I want to write this condition here so that implementing mul for references
// is defined and outputs f64
    // Mul<&'a Self, Output = f64> + Sized
{
}


enum IsingSpin {
    Up,
    Down,
}


impl Spin for IsingSpin {
}


impl<'a, 'b> Mul<&'a IsingSpin> for &'b IsingSpin {
    type Output = f64;

    fn mul(self, other: &'a IsingSpin) -> f64 {
        use IsingSpin::{Up, Down};
        match (self, other) {
            (&Up, &Up) => 1f64,
            (&Down, &Down) => 1f64,
            _ => -1f64,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::IsingSpin;

    #[test]
    fn ising_spin_multiplies_correctly() {
        assert_eq!(&IsingSpin::Up  * &IsingSpin::Up, 1.0);
        assert_eq!(&IsingSpin::Down  * &IsingSpin::Up, -1.0);
        assert_eq!(&IsingSpin::Up * &IsingSpin::Down, -1.0);
        assert_eq!(&IsingSpin::Down * &IsingSpin::Down, 1.0);
    }
}
