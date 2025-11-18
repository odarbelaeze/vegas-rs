//! `vegas` is a library that allows you to run atomistic simulations of
//! magnetic materials.
//!
//! # Table of Contents
//!
//! * [Installation](#installation)
//! * [Spins](#spins)
//! * [Hamiltonians](#hamiltonians)
//! * [Instruments](#instruments)
//! * [Machine](#machine)
//! * [Programs](#programs)
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
//! # Conceptual Overview
//!
//! ## Spins
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
//! ## Hamiltonians
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
//! * `Zeeman` - A hamiltonian that calculates the Zeeman energy of a spin system.
//! * `Compound` - A hamiltonian that combines multiple hamiltonians.
//!
//! ## Instruments
//!
//! Instruments allow you to measure properties of the spin system during the simulation.
//! The library provides an `Instrument` trait that you can implement for your own instruments.
//!
//! Among others this library provides the following instruments:
//!
//! * `StatSensor` - An instrument that measures the statistical properties of the spin system.
//! * `ObservableSensor` - An instrument that measures the observables of the spin system.
//! * `StateSensor` - An instrument that writes the state of the spin system.
//!
//! ## Machine
//!
//! Think of the `Machine` struct as a box that contains the sample you want to
//! measure. It holds the current state of the system, the hamiltonian that
//! describes the interactions, the integrator that updates the state, the
//! thermostat that controls the temperature, and the instruments that measure
//! properties of the system.
//!
//! A machine can relax or measure the system for a given number of steps during which
//! the thermal conditions don't change.
//!
//! The machine also provides hooks for changing instruments to measure different
//! properties during the simulation.
//!
//! ## Programs
//!
//! Programs are high-level routines that control the simulation process. They define sequences of
//! operations to be performed on the machine, such as cooling down the system or relaxing it at a
//! certain temperature.
//!
//! Programs implement the `Program` trait, which requires a `run` method. This method takes a
//! random number generator and a mutable reference to a `Machine`, allowing the program to
//! manipulate the machine's state and behavior.
//!
//! Some of the provided programs include:
//!
//! * `Relax` - A program that relaxes the system at a specified temperature for a given number of steps.
//! * `CoolDown` - A program that gradually cools down the system from a maximum temperature to a minimum temperature over a series of steps.
//! * `HysteresisLoop` - A program that simulates a hysteresis loop by varying the external magnetic field and measuring the system's response.
//!
//! # Example
//!
//! ```rust
//! use rand::SeedableRng;
//! use rand_pcg::Pcg64;
//! use vegas_lattice::Lattice;
//!
//! use vegas::{
//!     energy::Exchange,
//!     instrument::{Instrument, StatSensor},
//!     integrator::MetropolisIntegrator,
//!     machine::Machine,
//!     program::{CoolDown, Program, Relax},
//!     state::{Field, IsingSpin, State},
//!     thermostat::Thermostat,
//! };
//!
//! let lattice = Lattice::sc(1.0).expand_x(10).expand_y(10).drop_z();
//! let hamiltonian = Exchange::from_lattice(1.0, &lattice);
//! let program = CoolDown::default()
//!     .set_max_temperature(5.0)
//!     .set_steps(10)
//!     .set_relax(10);
//! let mut rng = Pcg64::from_rng(&mut rand::rng());
//! let state = State::<IsingSpin>::rand_with_size(&mut rng, lattice.sites().len());
//! let integrator = MetropolisIntegrator::new();
//! let thermostat = Thermostat::new(2.8, Field::zero());
//! let instruments: Vec<Box<dyn Instrument<_, _>>> =
//!     vec![Box::new(StatSensor::<_, _>::new(Box::new(std::io::stdout())))];
//! let mut machine = Machine::new(thermostat, hamiltonian, integrator, instruments, state);
//! let _ = program.run(&mut rng, &mut machine);
//! ```

#[macro_use]
pub mod energy;

pub mod accumulator;
pub mod error;
pub mod input;
pub mod instrument;
pub mod integrator;
pub mod machine;
pub mod output;
pub mod program;
pub mod state;
pub mod thermostat;
pub mod util;
