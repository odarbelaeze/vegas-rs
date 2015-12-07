use std::ops::Mul;

extern crate rand;
use rand::distributions::normal::StandardNormal;
use rand::distributions::{IndependentSample, Range};

extern crate vegas;

use vegas::lattice::Lattice;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;


struct Adjacency {
    lims: Vec<usize>,
    nbhs: Vec<usize>,
}


impl Adjacency {
    pub fn new(lattice: Lattice) -> Adjacency
    {
        let mut lims = vec![0];
        let mut nbhs = vec![];
        for site in lattice.sites() {
            let pnbhs: Vec<usize>  = lattice.targets(&site)
                .expect("No sites ma frien").iter().map(|site| {
                    lattice.index(&site).expect("Site outside lattice")
                }).collect();
            let last = lims.last().unwrap().clone();
            lims.push(last + pnbhs.len());
            nbhs.extend(pnbhs.iter());
        }
        Adjacency { lims: lims, nbhs: nbhs, }
    }

    fn nbhs_of<'a>(&'a self, item: usize) -> Option<&'a[usize]> {
        if item >= self.lims.len() - 1 {
            return None
        }
        let low = self.lims[item] as usize;
        let hi = self.lims[item + 1] as usize;
        Some(&self.nbhs[low..hi])
    }
}


#[derive(Copy, Clone)]
struct Spin {
    x: f64,
    y: f64,
    z: f64,
}


impl Spin {
    pub fn up() -> Spin {
        Spin { x: 0.0f64, y: 0.0f64, z: 1.0f64,  }
    }

    pub fn rand() -> Spin {
        let StandardNormal(x) = rand::random();
        let StandardNormal(y) = rand::random();
        let StandardNormal(z) = rand::random();
        let norm = (x * x + y * y + z * z).sqrt();
        Spin { x: x / norm, y: y / norm, z: z / norm, }
    }
}


impl Mul for Spin {

    type Output = f64;

    fn mul(self, other: Spin) ->  f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}


trait EnergyComponent {
    fn energy(&self, state: &Vec<Spin>, index: usize) ->  f64;
    fn total_energy(&self, state: &Vec<Spin>) -> f64 {
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
    fn energy(&self, state: &Vec<Spin>, index: usize) -> f64 {
        let mut ene = 0f64;
        let s = state[index];
        for nb in self.adjacency.nbhs_of(index).unwrap() {
            ene -= s * state[*nb];
        }
        ene
    }
}

fn step(energy: &ExchangeComponent, state: &Vec<Spin>, temp: f64) -> Vec<Spin> {
    let mut new_state = state.clone();
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

    let exchange = ExchangeComponent { adjacency: Adjacency::new(latt) };
    let mut state = (0..1000).map(|_| {
        Spin::rand()
    }).collect();

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
