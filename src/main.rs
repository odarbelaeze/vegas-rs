#[macro_use]
extern crate vegas;
extern crate clap;
extern crate sprs;
extern crate vegas_lattice;

use std::error::Error;
use std::fs::File;
use std::io::Read;

use clap::{Parser, Subcommand};
use vegas_lattice::{Axis, Lattice};

use vegas::energy::{Exchage, HamiltonianComponent};
use vegas::integrator::{Integrator, MetropolisIntegrator, StateGenerator};
use vegas::state::{IsingMagnetization, IsingSpin, Magnetization};

fn cool_down<T>(hamiltonian: T, len: usize)
where
    T: HamiltonianComponent<IsingSpin>,
{
    let mut integrator = MetropolisIntegrator::new(5.0);
    let mut state = integrator.state(len);
    loop {
        let steps = 5000;
        let mut energy_sum = 0.0;
        let mut magnetization_sum = IsingMagnetization::new();
        for _ in 0..steps {
            state = integrator.step(&hamiltonian, &state);
            energy_sum += hamiltonian.total_energy(&state);
            magnetization_sum = magnetization_sum
                + state
                    .spins()
                    .clone()
                    .into_iter()
                    .sum::<IsingMagnetization>();
        }
        println!(
            "{} {} {}",
            integrator.temp(),
            energy_sum / steps as f64,
            magnetization_sum.magnitude() / steps as f64
        );
        if integrator.temp() < 1.0 {
            break;
        }
        integrator.cool(0.1);
    }
}

fn bench() {
    let lattice = Lattice::sc(1.0)
        .expand_along(Axis::X, 10)
        .expand_along(Axis::Y, 10)
        .expand_along(Axis::Z, 10);
    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));
    cool_down(hamiltonian, lattice.sites().len());
}

fn bench_lattice(input: &str) -> Result<(), Box<dyn Error>> {
    let mut data = String::new();
    let mut file = File::open(input)?;
    file.read_to_string(&mut data)?;
    let lattice: Lattice = data.parse()?;
    println!("# Successfuly read the lattice!");

    println!("# Simulating with {} sites", lattice.sites().len());
    println!("# Simulating with {} exchanges", lattice.vertices().len());

    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));

    cool_down(hamiltonian, lattice.sites().len());
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
