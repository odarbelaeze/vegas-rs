#![feature(test)]
extern crate rand;
extern crate test;
extern crate vegas;

use vegas::hamiltonian::{Compound, Gauge, Hamiltonian, UniaxialAnisotropy};
use vegas::state::{HeisenbergSpin, Spin, State};
use vegas::thermostat::Thermostat;

#[bench]
fn gauge_energy_of_a_1_k_heisenberg_state(b: &mut test::Bencher) {
    let state = State::<HeisenbergSpin>::rand_with_size(&mut rand::rng(), 1_000);
    let gauge = Gauge::new(1.0);
    let thermostat = Thermostat::near_zero();
    b.iter(|| {
        let energy = gauge.total_energy(&thermostat, &state);
        assert!(energy == 1_000.0);
    })
}

#[bench]
fn anisotropy_energy_of_a_1k_heisenberg_state(b: &mut test::Bencher) {
    let state = State::<HeisenbergSpin>::up_with_size(1_000);
    let anisotropy = UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0);
    let thermostat = Thermostat::near_zero();
    b.iter(|| {
        let energy = anisotropy.total_energy(&thermostat, &state);
        assert!(energy == 1_000.0);
    })
}

#[bench]
fn compound_energy_hit_1k_heisenberg(b: &mut test::Bencher) {
    let state = State::<HeisenbergSpin>::up_with_size(1_000);
    let gauge = Gauge::new(1.0);
    let anisotropy = UniaxialAnisotropy::new(HeisenbergSpin::up(), 1.0);
    let compound = Compound::new(gauge, anisotropy);
    let thermostat = Thermostat::near_zero();
    b.iter(|| {
        let energy = compound.total_energy(&thermostat, &state);
        assert!(energy == 2_000.0);
    })
}
