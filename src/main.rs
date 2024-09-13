#[macro_use]
extern crate vegas;
extern crate clap;
extern crate rand_pcg;
extern crate sprs;
extern crate vegas_lattice;

use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use clap::{Parser, Subcommand};
use rand::SeedableRng;
use rand_pcg::Pcg64;
use vegas_lattice::Lattice;

use vegas::{
    error::Result,
    hamiltonian::Exchage,
    input::{Input, Model},
    integrator::{MetropolisFlipIntegrator, MetropolisIntegrator},
    machine::Machine,
    program::{CurieTemp, Program, Relax},
    state::{HeisenbergSpin, IsingSpin, State},
};

fn bench(lattice: Lattice, model: Model) -> Result<()> {
    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));
    match model {
        Model::Ising => {
            let program = CurieTemp::default().set_max_temp(5.0);
            let mut rng = Pcg64::from_entropy();
            let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
            program.run(&mut rng, &mut machine)
        }
        Model::Heisenberg => {
            let program = CurieTemp::default().set_max_temp(2.5);
            let mut rng = Pcg64::from_entropy();
            let state = State::<HeisenbergSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
            program.run(&mut rng, &mut machine)
        }
    }
}

fn bench_ising(length: usize) -> Result<()> {
    let lattice = Lattice::sc(1.0).expand_x(length).expand_y(length).drop_z();
    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));
    let cool_rate = 0.05;
    let relax = Relax::default().set_steps(500000).set_temp(2.8);
    let curie = CurieTemp::default()
        .set_steps(500000)
        .set_max_temp(2.8)
        .set_min_temp(1.8)
        .set_cool_rate(cool_rate);
    let mut rng = Pcg64::from_entropy();
    let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
    let integrator = MetropolisFlipIntegrator::new();
    let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
    relax.run(&mut rng, &mut machine)?;
    curie.run(&mut rng, &mut machine)
}

fn bench_model(model: Model, length: usize) -> Result<()> {
    let lattice = Lattice::sc(1.0).expand_all(length);
    bench(lattice, model)
}

fn bench_lattice(model: Model, input: &str) -> Result<()> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let lattice: Lattice = data.parse()?;

    println!("# Successfuly read the lattice!");
    println!("# Simulating with {} sites", lattice.sites().len());
    println!("# Simulating with {} exchanges", lattice.vertices().len());

    bench(lattice, model)
}

fn run_input(input: PathBuf) -> Result<()> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let input: Input = toml::from_str(&data)?;
    let mut rng = Pcg64::from_entropy();
    input.run(&mut rng)
}

fn print_default_input() -> Result<()> {
    let input = Input::default();
    let input = toml::to_string_pretty(&input)?;
    println!("{}", input);
    Ok(())
}

fn check_error(res: Result<()>) {
    if let Err(e) = res {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

#[derive(Debug, Subcommand)]
enum SubCommand {
    #[command(about = "Run 2D Ising model")]
    Ising {
        /// Length of the side lattice
        length: usize,
    },
    #[command(about = "Run benchmark")]
    Bench {
        /// Model to run
        model: Model,
        /// Length of the side lattice
        length: usize,
    },
    #[command(about = "Simulate the given lattice")]
    Lattice {
        /// Model to run
        model: Model,
        /// Path to the lattice file
        lattice: String,
    },
    #[command(about = "Generate a sample input file")]
    Input,
    #[command(about = "Run the program")]
    Run { input: PathBuf },
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about)]
struct Cli {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

fn main() {
    let cli = Cli::parse();
    match cli.subcmd {
        SubCommand::Ising { length } => check_error(bench_ising(length)),
        SubCommand::Bench { length, model } => check_error(bench_model(model, length)),
        SubCommand::Lattice { lattice, model } => check_error(bench_lattice(model, &lattice)),
        SubCommand::Input => check_error(print_default_input()),
        SubCommand::Run { input } => check_error(run_input(input)),
    }
}
