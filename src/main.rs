#[macro_use] extern crate vegas_rs;
extern crate docopt;
extern crate vegas_lattice;
extern crate sprs;


use std::error::Error;
use std::fs::File;
use std::io::Read;

use docopt::Docopt;
use sprs::TriMat;
use vegas_lattice::Lattice;

use vegas_rs::state::{State, HeisenbergSpin};
use vegas_rs::energy::{EnergyComponent, Gauge, ExchangeEnergy};
use vegas_rs::integrator::{Integrator, StateGenerator, MetropolisIntegrator};


const USAGE: &'static str = "
Vegas rust.

Usage:
  vegas bench
  vegas lattice <lattice>
  vegas (-h | --help)
  vegas --version

Options:
  -h --help     Show this screen.
  --version     Show version.
";

const VERSION: &'static str = "
Vegas rust, version: 0.1.0
";


fn cool_down<T: EnergyComponent<HeisenbergSpin>>(hamiltonian: T, len: usize) {
    let mut integrator = MetropolisIntegrator::new(3.0);
    let mut state: State<HeisenbergSpin> = integrator.state(len);
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



fn bench() {
    let hamiltonian = hamiltonian!(
        Gauge::new(10.0)
    );
    cool_down(hamiltonian, 100);
}


fn bench_lattice(input: &str) -> Result<(), Box<Error>> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let lattice: Lattice = data.parse()?;
    println!("# Successfuly read the lattice!");

    let nsites = lattice.sites().len();

    let mut mat = TriMat::new((nsites, nsites));
    for vertex in lattice.vertices() {
        if vertex.source() <= vertex.target() {
            mat.add_triplet(vertex.source(), vertex.target(), 1.0);
            mat.add_triplet(vertex.target(), vertex.source(), 1.0);
        }
    }

    let csr = mat.to_csr();

    println!("# Simulating with {} sites", nsites);
    println!("# Simulating with {} exchanges", lattice.vertices().len());

    let hamiltonian = hamiltonian!(
        ExchangeEnergy::new(csr)
    );

    cool_down(hamiltonian, nsites);
    Ok(())
}


fn check_error(res: Result<(), Box<Error>>) {
    match res {
        Err(e) => {
            eprintln!("Error: {}", e.description());
            match e.cause() {
                Some(cause) => eprintln!("Cause: {}", cause),
                _ => ()
            };
            std::process::exit(1);
        },
        _ => {},
    }
}


pub fn main() {
    let version = VERSION.trim().to_string();
    let args = Docopt::new(USAGE)
        .and_then(|doc| doc.version(Some(version)).parse())
        .unwrap_or_else(|e| e.exit());
    if args.get_bool("bench") {
        bench()
    } else if args.get_bool("lattice") {
        check_error(bench_lattice(args.get_str("<lattice>")))
    }
}
