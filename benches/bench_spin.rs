#![feature(test)]
extern crate vegas_rs;
extern crate test;

use vegas_rs::state::{Spin, IsingSpin, HeisenbergSpin};


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
