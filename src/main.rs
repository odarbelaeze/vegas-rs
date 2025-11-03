//! A command line interface for running Vegas simulations and benchmarks.

use clap::{Parser, Subcommand};
use rand::SeedableRng;
use rand_pcg::Pcg64;
use std::{
    fs::File,
    io::{Read, stdin, stdout},
    path::PathBuf,
};
use vegas::{
    energy::Exchange,
    error::{IoError, VegasResult},
    input::{MetropolisInput, Model, WolffInput},
    instrument::{Instrument, StatSensor},
    integrator::MetropolisIntegrator,
    machine::Machine,
    program::{CoolDown, Program},
    state::{Field, HeisenbergSpin, IsingSpin, State},
    thermostat::Thermostat,
};
use vegas_lattice::Lattice;

fn bench(lattice: Lattice, model: Model) -> VegasResult<()> {
    let hamiltonian = Exchange::from_lattice(&lattice);
    match model {
        Model::Ising => {
            let program = CoolDown::default()
                .set_max_temperature(5.0)
                .set_cool_rate(0.05);
            let mut rng = Pcg64::from_rng(&mut rand::rng());
            let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let thermostat = Thermostat::new(2.8, Field::zero());
            let instruments: Vec<Box<dyn Instrument<_, _>>> =
                vec![Box::new(StatSensor::<_, _>::new(Box::new(stdout())))];
            let mut machine = Machine::new(thermostat, hamiltonian, integrator, instruments, state);
            program.run(&mut rng, &mut machine)?;
            Ok(())
        }
        Model::Heisenberg => {
            let program = CoolDown::default()
                .set_max_temperature(2.5)
                .set_cool_rate(0.05);
            let mut rng = Pcg64::from_rng(&mut rand::rng());
            let state = State::<HeisenbergSpin>::rand_with_size(&mut rng, lattice.sites().len());
            let integrator = MetropolisIntegrator::new();
            let thermostat = Thermostat::new(2.8, Field::zero());
            let instruments: Vec<Box<dyn Instrument<_, _>>> =
                vec![Box::new(StatSensor::<_, _>::new(Box::new(stdout())))];
            let mut machine = Machine::new(thermostat, hamiltonian, integrator, instruments, state);
            program.run(&mut rng, &mut machine)?;
            Ok(())
        }
    }
}

fn bench_model(model: Model, length: usize) -> VegasResult<()> {
    let lattice = Lattice::sc(1.0).expand_all(length);
    bench(lattice, model)
}

fn run_metropolis(input: PathBuf) -> VegasResult<()> {
    let mut data = String::new();
    if input == PathBuf::from("-") {
        stdin().read_to_string(&mut data).map_err(IoError::from)?;
    } else {
        let mut file = File::open(input).map_err(IoError::from)?;
        file.read_to_string(&mut data).map_err(IoError::from)?;
    };
    let input: MetropolisInput = toml::from_str(&data)?;
    let mut rng = Pcg64::from_rng(&mut rand::rng());
    input.run(&mut rng)
}

fn run_wolff(input: PathBuf) -> VegasResult<()> {
    let mut data = String::new();
    if input == PathBuf::from("-") {
        stdin().read_to_string(&mut data).map_err(IoError::from)?;
    } else {
        let mut file = File::open(input).map_err(IoError::from)?;
        file.read_to_string(&mut data).map_err(IoError::from)?;
    };
    let input: WolffInput = toml::from_str(&data)?;
    let mut rng = Pcg64::from_rng(&mut rand::rng());
    input.run(&mut rng)
}

fn print_default_input() -> VegasResult<()> {
    let input = MetropolisInput::default();
    let input = toml::to_string_pretty(&input)?;
    println!("{}", input);
    Ok(())
}

fn check_error(res: VegasResult<()>) {
    if let Err(e) = res {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

#[derive(Debug, Subcommand)]
enum SubCommand {
    /// Run benchmarks
    Bench {
        /// Model to run
        model: Model,
        /// Length of the side lattice
        length: usize,
    },
    /// Print the default input
    Input,
    /// Run a Metropolis simulation from an input file
    Metropolis { input: PathBuf },
    /// Run a Wolf simulation from an input file
    Wolff { input: PathBuf },
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
        SubCommand::Bench { length, model } => check_error(bench_model(model, length)),
        SubCommand::Input => check_error(print_default_input()),
        SubCommand::Metropolis { input } => check_error(run_metropolis(input)),
        SubCommand::Wolff { input } => check_error(run_wolff(input)),
    }
}
