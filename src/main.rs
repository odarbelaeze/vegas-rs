#[macro_use] extern crate vegas_rs;

use vegas_rs::state::{State, HeisenbergSpin};
use vegas_rs::energy::{EnergyComponent, Gauge};
use vegas_rs::integrator::{Integrator, StateGenerator, MetropolisIntegrator};


pub fn main() {


    let hamiltonian = hamiltonian!(
        Gauge::new(10.0)
    );

    let mut integrator = MetropolisIntegrator::new(3.0);
    let mut state: State<HeisenbergSpin> = integrator.state(100);

    loop {
        let steps = 1000;
        let mut energy_sum = 0.0;
        for _ in 0..steps {
            state = integrator.step(&hamiltonian, &state);
            energy_sum += hamiltonian.total_energy(&state)
        }
        println!("{} {}", integrator.temp(), energy_sum / steps as f64);
        if integrator.temp() < 0.1 { break }
        integrator.cool(0.1);
    }

}
