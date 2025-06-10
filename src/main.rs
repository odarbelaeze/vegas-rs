use std::{
    fs::File,
    io::{stdin, Read},
    path::{Path, PathBuf},
};

use clap::{Parser, Subcommand};
use rand::SeedableRng;
use rand_pcg::Pcg64;
use vegas_lattice::Lattice;

use vegas::{
    error::{IOError, Result},
    hamiltonian::Exchange,
    input::{Input, Model},
    integrator::{MetropolisFlipIntegrator, MetropolisIntegrator},
    machine::Machine,
    program::{CoolDown, Program, Relax},
    state::{HeisenbergSpin, IsingSpin, State},
};

fn bench(lattice: Lattice, model: Model) -> Result<()> {
    let hamiltonian = Exchange::from_lattice(&lattice);
    match model {
        Model::Ising => {
            let program = CoolDown::default().set_max_temperature(5.0);
            let mut rng = Pcg64::from_rng(&mut rand::rng());
            let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
            program.run(&mut rng, &mut machine, &mut None)
        }
        Model::Heisenberg => {
            let program = CoolDown::default().set_max_temperature(2.5);
            let mut rng = Pcg64::from_rng(&mut rand::rng());
            let state = State::<HeisenbergSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
            program.run(&mut rng, &mut machine, &mut None)
        }
    }
}

fn bench_ising(length: usize) -> Result<()> {
    let lattice = Lattice::sc(1.0).expand_x(length).expand_y(length).drop_z();
    let hamiltonian = Exchange::from_lattice(&lattice);
    let cool_rate = 0.05;
    let relax = Relax::default().set_steps(500000).set_temperature(2.8);
    let curie = CoolDown::default()
        .set_steps(500000)
        .set_max_temperature(2.8)
        .set_min_temperature(1.8)
        .set_cool_rate(cool_rate);
    let mut rng = Pcg64::from_rng(&mut rand::rng());
    let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
    let integrator = MetropolisFlipIntegrator::new();
    let mut machine = Machine::new(2.8, 0.0, hamiltonian, integrator, state);
    relax.run(&mut rng, &mut machine, &mut None)?;
    curie.run(&mut rng, &mut machine, &mut None)
}

fn bench_model(model: Model, length: usize) -> Result<()> {
    let lattice = Lattice::sc(1.0).expand_all(length);
    bench(lattice, model)
}

fn bench_lattice(model: Model, input: &Path) -> Result<()> {
    let mut data = String::new();
    let mut file = File::open(input).map_err(IOError::from)?;
    file.read_to_string(&mut data).map_err(IOError::from)?;
    let lattice: Lattice = data.parse()?;

    println!("# Successfuly read the lattice!");
    println!("# Simulating with {} sites", lattice.sites().len());
    println!("# Simulating with {} exchanges", lattice.vertices().len());

    bench(lattice, model)
}

fn run_input(input: PathBuf) -> Result<()> {
    let mut data = String::new();
    if input == PathBuf::from("-") {
        stdin().read_to_string(&mut data).map_err(IOError::from)?;
    } else {
        let mut file = File::open(input).map_err(IOError::from)?;
        file.read_to_string(&mut data).map_err(IOError::from)?;
    };
    let input: Input = toml::from_str(&data)?;
    let mut rng = Pcg64::from_rng(&mut rand::rng());
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
    /// Run a 2D Ising model
    Ising {
        /// Length of the side lattice
        length: usize,
    },
    /// Run benchmarks
    Bench {
        /// Model to run
        model: Model,
        /// Length of the side lattice
        length: usize,
    },
    /// Run a lattice
    Lattice {
        /// Model to run
        model: Model,
        /// Path to the lattice file
        lattice: PathBuf,
    },
    /// Print the default input
    Input,
    /// Run an input file
    Run { input: PathBuf },
}

#[derive(Debug, Parser)]
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
