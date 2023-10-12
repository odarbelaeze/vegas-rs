#![feature(test)]
extern crate rand;
extern crate test;
extern crate vegas_rs;

use rand::XorShiftRng;
use vegas_rs::state::{HeisenbergSpin, IsingSpin, Spin, State};

#[bench]
fn create_1k_ising_spins(b: &mut test::Bencher) {
    let mut rng = XorShiftRng::new_unseeded();
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = IsingSpin::rand(&mut rng);
        }
    });
}

#[bench]
fn create_1k_heisenberg_spins(b: &mut test::Bencher) {
    let mut rng = XorShiftRng::new_unseeded();
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = HeisenbergSpin::rand(&mut rng);
        }
    });
}

#[bench]
fn create_a_1k_spin_state(b: &mut test::Bencher) {
    let mut rng = XorShiftRng::new_unseeded();
    b.iter(|| {
        let _state = State::<HeisenbergSpin>::rand_with_size(1_000, &mut rng);
    })
}
