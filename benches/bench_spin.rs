#![feature(test)]
extern crate rand;
extern crate test;
extern crate vegas;

use rand::{rngs::SmallRng, SeedableRng};
use rand_pcg::Pcg64;
use vegas::state::{HeisenbergSpin, IsingSpin, Spin, State};

#[bench]
fn create_1k_ising_spins(b: &mut test::Bencher) {
    let mut rng = Pcg64::from_entropy();
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = IsingSpin::rand(&mut rng);
        }
    });
}

#[bench]
fn create_1k_heisenberg_spins(b: &mut test::Bencher) {
    let mut rng = Pcg64::from_entropy();
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = HeisenbergSpin::rand(&mut rng);
        }
    });
}

#[bench]
fn create_1k_heisenberg_spins_with_small_rng(b: &mut test::Bencher) {
    let mut rng = SmallRng::from_entropy();
    b.iter(|| {
        for _ in 0..1_000 {
            let _spin = HeisenbergSpin::rand(&mut rng);
        }
    });
}

#[bench]
fn create_a_1k_spin_state(b: &mut test::Bencher) {
    let mut rng = Pcg64::from_entropy();
    b.iter(|| {
        let _state = State::<HeisenbergSpin>::rand_with_size(&mut rng, 1_000);
    })
}
