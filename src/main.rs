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
use vegas_rs::energy::{EnergyComponent, Gauge};
use vegas_rs::integrator::{Integrator, StateGenerator, MetropolisIntegrator};


const USAGE: &'static str = "
Vegas rust.

Usage:
  vegas bench
  vegas read <lattice>
  vegas (-h | --help)
  vegas --version

Options:
  -h --help     Show this screen.
  --version     Show version.
";

const VERSION: &'static str = "
Vegas rust, version: 0.1.0
";


fn bench() {
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


fn read(input: &str) -> Result<(), Box<Error>> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let lattice: Lattice = data.parse()?;
    println!("Successfuly read the lattice!");
    let mut mat = TriMat::new((lattice.sites().len(), lattice.sites().len()));
    for vertex in lattice.vertices() {
        println!("{:?}", vertex);
        mat.add_triplet(vertex.source(), vertex.target(), 1.0);
    }
    
    println!("{:?}", mat);
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
    } else if args.get_bool("read") {
        check_error(read(args.get_str("<lattice>")))
    }
}
