#![feature(test)]
extern crate vegas_rs;
extern crate test;

use vegas_rs::state::{Spin, IsingSpin, HeisenbergSpin, State};


#[bench]
fn create_10k_ising_spins(b: &mut test::Bencher) {
    b.iter(|| {
        for _ in 0..10_000 {
            let _spin = IsingSpin::rand();
        }
    });
}

#[bench]
fn create_10k_heisenberg_spins(b: &mut test::Bencher) {
    b.iter(|| {
        for _ in 0..10_000 {
            let _spin = HeisenbergSpin::rand();
        }
    });
}

#[bench]
fn create_a_10k_spin_state(b: &mut test::Bencher) {
    b.iter(|| {
        let _state = State::<HeisenbergSpin>::rand_with_size(10_000);
    })
}
