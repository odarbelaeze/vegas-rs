#[macro_use]
extern crate vegas;
extern crate clap;
extern crate sprs;
extern crate vegas_lattice;

use std::error::Error;
use std::fs::File;
use std::io::Read;

use clap::{Parser, Subcommand};
use sprs::TriMat;
use vegas_lattice::Lattice;

use vegas::energy::{Exchage, Gauge, HamiltonianComponent};
use vegas::integrator::{Integrator, MetropolisIntegrator, StateGenerator};
use vegas::state::{HeisenbergSpin, State};

fn cool_down<T: HamiltonianComponent<HeisenbergSpin>>(hamiltonian: T, len: usize) {
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
        if integrator.temp() < 0.1 {
            break;
        }
        integrator.cool(0.1);
    }
}

fn bench() {
    let hamiltonian = hamiltonian!(Gauge::new(10.0));
    cool_down(hamiltonian, 100);
}

fn bench_lattice(input: &str) -> Result<(), Box<dyn Error>> {
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

    let hamiltonian = hamiltonian!(Exchage::new(csr));

    cool_down(hamiltonian, nsites);
    Ok(())
}

fn check_error(res: Result<(), Box<dyn Error>>) {
    if let Err(e) = res {
        eprintln!("Error: {}", e);
        if let Some(source) = e.source() {
            eprintln!("source: {}", source);
        }
        std::process::exit(1);
    }
}

#[derive(Debug, Subcommand)]
enum SubCommand {
    #[command(about = "Run benchmark")]
    Bench {},
    #[command(about = "Simulate the given lattice")]
    Lattice { lattice: String },
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about=None)]
struct Cli {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

fn main() {
    let cli = Cli::parse();
    match cli.subcmd {
        SubCommand::Bench {} => bench(),
        SubCommand::Lattice { lattice } => check_error(bench_lattice(&lattice)),
    }
}
