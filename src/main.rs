#[macro_use]
extern crate vegas;
extern crate clap;
extern crate rand_pcg;
extern crate sprs;
extern crate vegas_lattice;

use std::error::Error;
use std::fs::File;
use std::io::Read;

use clap::{Parser, Subcommand};
use rand::SeedableRng;
use rand_pcg::Pcg64;
use vegas::program::CurieTemp;
use vegas_lattice::{Axis, Lattice};

use vegas::energy::Exchage;
use vegas::integrator::MetropolisIntegrator;
use vegas::state::{HeisenbergSpin, IsingSpin, State};

fn bench(lattice: Lattice, model: &str) {
    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));
    match model {
        "ising" => {
            let program = CurieTemp::default().with_max_temp(5.0);
            let mut rng = Pcg64::from_entropy();
            let state = State::<IsingSpin>::rand_with_size(lattice.sites().len(), &mut rng);
            let mut integrator = MetropolisIntegrator::new(rng);
            program.run(&mut integrator, &hamiltonian, state);
        }
        "heisenberg" => {
            let program = CurieTemp::default().with_max_temp(2.5);
            let mut rng = Pcg64::from_entropy();
            let state = State::<HeisenbergSpin>::rand_with_size(lattice.sites().len(), &mut rng);
            let mut integrator = MetropolisIntegrator::new(rng);
            program.run(&mut integrator, &hamiltonian, state);
        }
        _ => eprintln!("Unknown model: {}", model),
    }
}

fn bench_model(length: usize, model: &str) {
    let lattice = Lattice::sc(1.0)
        .expand_along(Axis::X, length)
        .expand_along(Axis::Y, length)
        .expand_along(Axis::Z, length);
    bench(lattice, model);
}

fn bench_lattice(model: &str, input: &str) -> Result<(), Box<dyn Error>> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let lattice: Lattice = data.parse()?;

    println!("# Successfuly read the lattice!");
    println!("# Simulating with {} sites", lattice.sites().len());
    println!("# Simulating with {} exchanges", lattice.vertices().len());

    bench(lattice, model);

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
    Bench { model: String, length: usize },
    #[command(about = "Simulate the given lattice")]
    Lattice { model: String, lattice: String },
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
        SubCommand::Bench { length, model } => bench_model(length, &model),
        SubCommand::Lattice { lattice, model } => check_error(bench_lattice(&model, &lattice)),
    }
}
