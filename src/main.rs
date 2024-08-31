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
use rand_pcg::Pcg64;
use vegas_lattice::{Axis, Lattice};

use vegas::energy::{Exchage, HamiltonianComponent};
use vegas::integrator::{Integrator, MetropolisIntegrator, StateGenerator};
use vegas::state::{HeisenbergSpin, IsingSpin, Magnetization, Spin};

fn cool_down<T, S>(hamiltonian: T, len: usize)
where
    S: Spin,
    T: HamiltonianComponent<S>,
{
    let mut integrator = MetropolisIntegrator::<Pcg64>::new(5.0);
    let mut state = integrator.state(len);
    loop {
        let relax = 1000;
        let steps = 100000;
        let mut energy_sum = 0.0;
        let mut magnetization_sum = 0.0;
        let mut mag_square_sum = 0.0;
        let mut mag_to_the_fourth_sum = 0.0;
        for _ in 0..relax {
            state = integrator.step(&hamiltonian, &state);
        }
        for _ in 0..steps {
            state = integrator.step(&hamiltonian, &state);
            energy_sum += hamiltonian.total_energy(&state);
            let magnetization = state.magnetization().magnitude() / len as f64;
            magnetization_sum += magnetization;
            mag_square_sum += magnetization.powi(2);
            mag_to_the_fourth_sum += magnetization.powi(4);
        }
        println!(
            "{} {} {} {} {}",
            integrator.temp(),
            energy_sum / steps as f64,
            magnetization_sum / steps as f64,
            (mag_square_sum / steps as f64 - (magnetization_sum / steps as f64).powi(2))
                / integrator.temp(),
            1.0 - (mag_to_the_fourth_sum / steps as f64) / 3.0
                * (mag_square_sum / steps as f64).powi(2)
        );
        if integrator.temp() < 0.1 {
            break;
        }
        integrator.cool(0.05);
    }
}

fn bench(length: usize, model: &str) {
    let lattice = Lattice::sc(1.0)
        .expand_along(Axis::X, length)
        .expand_along(Axis::Y, length)
        .expand_along(Axis::Z, length);
    let len = lattice.sites().len();
    let hamiltonian = hamiltonian!(Exchage::from_lattice(&lattice));
    match model {
        "ising" => cool_down::<_, IsingSpin>(hamiltonian, len),
        "heisenberg" => cool_down::<_, HeisenbergSpin>(hamiltonian, len),
        _ => eprintln!("Unknown model: {}", model),
    }
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

    cool_down::<_, HeisenbergSpin>(hamiltonian, lattice.sites().len());
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
        SubCommand::Bench { length, model } => bench(length, &model),
        SubCommand::Lattice { lattice } => check_error(bench_lattice(&lattice)),
    }
}
