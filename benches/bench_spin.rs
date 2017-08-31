#![feature(test)]
extern crate vegas_rs;
extern crate test;
extern crate rand;

use vegas_rs::state::{Spin, IsingSpin, HeisenbergSpin, State};
use rand::XorShiftRng;


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
