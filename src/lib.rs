use std::ops::{Mul};


trait Spin where
    for<'a, 'b> &'a Self: Mul<&'b Self, Output = f64>
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


struct HeisenbergSpin([f64; 3]);


impl Spin for HeisenbergSpin {}


impl<'a, 'b> Mul<&'a HeisenbergSpin> for &'b HeisenbergSpin {
    type Output = f64;

    fn mul(self, other: &'a HeisenbergSpin) -> f64 {
        let &HeisenbergSpin(_self) = self;
        let &HeisenbergSpin(_other) = other;
        _self.iter().zip(_other.iter()).map(|(a, b)| a * b).fold(0f64, |sum, i| sum + i)
    }
}


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
