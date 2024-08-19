//! `vegas` is a library that allows you to create atomistic simulations
//! of magnetic materials.
//!
//! # Table of Contentsble
//!
//! * [Installation](#installation)
//! * [Spins](#spins)
//! * [Hamiltonians](#hamiltonians)
//! * [Example](#example)
//!
//! # Installation
//!
//! Install the `vegas` library by adding the following to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! vegas = "*"
//! ```
//!
//! # Spins
//!
//! Spins are the basic building blocks of the library. They represent the spin
//! of an atom in a magnetic material. The library provides a `Spin` trait that
//! you can implement for your own spin types.
//!
//! The library provides a `HeisenbergSpin` type that represents a spin in a
//! Heisenberg model. The `HeisenbergSpin` type is a three-dimensional vector
//! that represents the spin of an atom.
//!
//! Furthermore, the library provides an `IsingSpin` type that represents a spin
//! in an Ising model. The `IsingSpin` type implemented as an enum that can take
//! the up or down variants.
//!
//! # Hamiltonians
//!
//! A hamiltonian is a function that calculates the energy of a spin system.
//!
//! This library provides a `EnergyComponent` trait that you can implement for
//! your own hamiltonians.
//!
//! Among others this library provides the following hamiltonians:
//!
//! * `Exchange` - A hamiltonian that calculates the exchange energy of a spin system.
//! * `Gauge` - A hamiltonian that calculates the gauge energy of a spin system.
//! * `UniaxialAnisotropy` - A hamiltonian that calculates the uniaxial anisotropy energy of a spin system.
//! * `ZeemanEnergy` - A hamiltonian that calculates the Zeeman energy of a spin system.
//! * `Compound` - A hamiltonian that combines multiple hamiltonians.
//!
//! # Example
//!
//! ```rust
//! use vegas::hamiltonian;
//! use vegas::integrator::{Integrator, MetropolisIntegrator, StateGenerator};
//! use vegas::energy::{HamiltonianComponent, Gauge};
//! use vegas::state::{HeisenbergSpin, State};
//!
//! let hamiltonian = hamiltonian!(Gauge::new(10.0));
//! let mut integrator = MetropolisIntegrator::new(3.0);
//! let mut state: State<HeisenbergSpin> = integrator.state(1000);
//! loop {
//!     let steps = 1000;
//!     let mut energy_sum = 0.0;
//!     for _ in 0..steps {
//!         state = integrator.step(&hamiltonian, &state);
//!         energy_sum += hamiltonian.total_energy(&state)
//!     }
//!     println!("{} {}", integrator.temp(), energy_sum / steps as f64);
//!     if integrator.temp() < 0.1 {
//!         break;
//!     }
//!     integrator.cool(0.1);
//! }
//! ```

extern crate rand;
extern crate sprs;

pub mod energy;
pub mod integrator;
pub mod state;
