
extern crate rand;

use rand::distributions::{IndependentSample, Range};

extern crate vegas;

use vegas::lattice::Adjacency;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;
use vegas::state::Spin;
use vegas::state::State;
use vegas::state::SpinConstructors;
use vegas::state::StateConstructors;



trait EnergyComponent {
    fn energy(&self, state: &State, index: usize) -> f64;
    fn total_energy(&self, state: &State) -> f64 {
        let mut total = 0f64;
        for i in 0..state.len() {
            total += self.energy(state, i);
        }
        total
    }
}


struct ExchangeComponent {
    adjacency: Adjacency,
}


impl EnergyComponent for ExchangeComponent {
    fn energy(&self, state: &State, index: usize) -> f64 {
        let mut ene = 0f64;
        let s = state[index];
        for nb in self.adjacency.nbhs_of(index).unwrap() {
            ene -= s * state[*nb];
        }
        ene
    }
}


fn step(energy: &ExchangeComponent, state: &State, temp: f64) -> State {
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
        .pbc((true, true, false))
        .shape((30, 30, 1))
        .vertices(Vertex::list_for_cubic())
        .natoms(1)
        .finalize();

    let mut state = State::rand(latt.nsites());
    let exchange = ExchangeComponent { adjacency: Adjacency::new(latt) };
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
