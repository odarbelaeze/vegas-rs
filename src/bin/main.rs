#[macro_use] extern crate vegas;

use vegas::lattice::Adjacency;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;
use vegas::state::State;
use vegas::state::StateConstructors;
use vegas::state::CommonObservables;
use vegas::energy::EnergyComponent;
use vegas::energy::ExchangeComponent;
use vegas::integrator::Integrator;
use vegas::integrator::MetropolisIntegrator;


pub fn main() {

    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((10, 10, 10))
        .vertices(Vertex::list_for_cubic())
        .natoms(1)
        .finalize();

    let mut state = State::rand(latt.nsites());

    let hamiltonian = hamiltonian!(
        ExchangeComponent::new(Adjacency::new(&latt))
    );

    let mut integrator = MetropolisIntegrator::new(3.0);

    loop {
        for _ in 0..1000 {
            state = integrator.step(&hamiltonian, &state);
            println!("{} {}", hamiltonian.total_energy(&state), state.mag_len());
        }
        if integrator.temp() < 0.1 {
            break
        }
        integrator.cool(0.1);
        println!("");
    }

}
