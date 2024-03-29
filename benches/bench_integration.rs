#![feature(test)]

#[macro_use]
extern crate vegas_rs;
extern crate test;

use vegas_rs::energy::{Gauge, UniaxialAnisotropy};
use vegas_rs::integrator::{Integrator, MetropolisIntegrator, StateGenerator};
use vegas_rs::state::{HeisenbergSpin, IsingSpin, Spin, State};

#[bench]
fn integration_of_1k_heisenberg_spin_with_gauge(b: &mut test::Bencher) {
    let gauge = Gauge::new(1.0);
    let mut integrator = MetropolisIntegrator::new(3.0);
    let mut state: State<HeisenbergSpin> = integrator.state(1_000);
    b.iter(|| state = integrator.step(&gauge, &state))
}

#[bench]
fn integration_of_1k_ising_spin_with_gauge(b: &mut test::Bencher) {
    let gauge = Gauge::new(1.0);
    let mut integrator = MetropolisIntegrator::new(3.0);
    let mut state: State<IsingSpin> = integrator.state(1_000);
    b.iter(|| state = integrator.step(&gauge, &state))
}

#[bench]
fn integration_of_1k_heisenberg_spin_with_compound_energy(b: &mut test::Bencher) {
    let hamiltonian = hamiltonian!(
        Gauge::new(1.0),
        UniaxialAnisotropy::new(HeisenbergSpin::up(), 10.0)
    );
    let mut integrator = MetropolisIntegrator::new(3.0);
    let mut state: State<HeisenbergSpin> = integrator.state(1_000);
    b.iter(|| state = integrator.step(&hamiltonian, &state))
}

