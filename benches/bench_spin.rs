#![feature(test)]
extern crate vegas_rs;
extern crate test;

use vegas_rs::state::{Spin, IsingSpin, HeisenbergSpin, State};


#[bench]
fn create_1k_ising_spins(b: &mut test::Bencher) {
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = IsingSpin::rand();
        }
    });
}

#[bench]
fn create_1k_heisenberg_spins(b: &mut test::Bencher) {
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = HeisenbergSpin::rand();
        }
    });
}

#[bench]
fn create_a_1k_spin_state(b: &mut test::Bencher) {
    b.iter(|| {
        let _state = State::<HeisenbergSpin>::rand_with_size(1_000);
    })
}
