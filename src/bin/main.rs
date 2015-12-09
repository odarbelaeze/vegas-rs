extern crate rand;
use rand::distributions::{IndependentSample, Range};


#[macro_use] extern crate vegas;
use vegas::lattice::Adjacency;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;
use vegas::state::Spin;
use vegas::state::State;
use vegas::state::SpinConstructors;
use vegas::state::StateConstructors;
use vegas::energy::EnergyComponent;
use vegas::energy::ExchangeComponent;


fn step<T: EnergyComponent>(energy: &T, state: &State, temp: f64) -> State {
    let mut new_state = (*state).clone();
    let sites = Range::new(0, new_state.len());
    let mut rng = rand::thread_rng();
    for _ in 0..new_state.len() {
        let site = sites.ind_sample(&mut rng);
        let old_energy = energy.energy(&new_state, site);
        new_state[site] = Spin::rand();
        let new_energy = energy.energy(&new_state, site);
        let delta = new_energy - old_energy;
        if delta < 0.0 {
            continue
        }
        if rand::random::<f64>() < (- delta / temp).exp() {
            continue
        }
        new_state[site] = state[site];
    }
    new_state
}


fn main() {

    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((10, 10, 10))
        .vertices(Vertex::list_for_cubic())
        .natoms(1)
        .finalize();

    let mut state = State::rand(latt.nsites());

    let exchange = hamiltonian!(
        ExchangeComponent::new(Adjacency::new(&latt))
    );

    let mut temp = 3.0;

    loop {
        for _ in 0..1000 {
            state = step(&exchange, &state, temp);
            println!("{}", exchange.total_energy(&state));
        }
        if temp < 0.1 {
            break
        }
        temp -= 0.1;
        println!("");
    }

}
