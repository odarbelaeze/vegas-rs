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
    input::{Input, Model},
    instrument::{Instrument, StatSensor},
    integrator::MetropolisIntegrator,
    machine::Machine,
    program::{CoolDown, Program},
    state::{Field, HeisenbergSpin, IsingSpin, State},
    thermostat::Thermostat,
};
use vegas_lattice::Lattice;

fn bench_model(model: Model, length: usize, seed: Option<u64>) -> VegasResult<()> {
    let lattice = Lattice::sc(1.0).expand_all(length);
    let hamiltonian = Exchange::from_lattice(1.0, &lattice);
    let mut rng = match seed {
        Some(s) => Pcg64::seed_from_u64(s),
        None => Pcg64::from_rng(&mut rand::rng()),
    };
    match model {
        Model::Ising => {
            let program = CoolDown::default()
                .set_max_temperature(5.0)
                .set_cool_rate(0.05);
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

fn run_input(input: PathBuf, seed: Option<u64>) -> VegasResult<()> {
    let mut data = String::new();
    if &input == "-" {
        stdin().read_to_string(&mut data).map_err(IoError::from)?;
    } else {
        let mut file = File::open(input).map_err(IoError::from)?;
        file.read_to_string(&mut data).map_err(IoError::from)?;
    };
    let input: Input = toml::from_str(&data)?;
    let mut rng = match seed {
        Some(s) => Pcg64::seed_from_u64(s),
        None => Pcg64::from_rng(&mut rand::rng()),
    };
    input.run(&mut rng)
}

fn print_default_input() -> VegasResult<()> {
    let input = Input::default();
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
    /// Print the default input
    Input,
    /// Run benchmarks
    Bench {
        /// Model to run
        model: Model,
        /// Length of the side lattice
        length: usize,
        /// Seed for RNG, random if omitted
        #[arg(short, long)]
        seed: Option<u64>,
    },
    /// Run a simulation from an input file
    Run {
        /// Input file, or "-" for stdin
        input: PathBuf,
        /// Seed for RNG, random if omitted
        #[arg(short, long)]
        seed: Option<u64>,
    },
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
        SubCommand::Input => check_error(print_default_input()),
        SubCommand::Bench {
            length,
            model,
            seed,
        } => check_error(bench_model(model, length, seed)),
        SubCommand::Run { input, seed } => check_error(run_input(input, seed)),
    }
}
