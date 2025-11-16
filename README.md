# vegas

[![Crates.io](https://img.shields.io/crates/v/vegas.svg)](https://crates.io/crates/vegas)
[![Build Status](https://github.com/odarbelaeze/vegas-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/odarbelaeze/vegas-rs/actions/workflows/rust.yml)
[![Documentation](https://docs.rs/vegas/badge.svg)](https://docs.rs/vegas)
[![DOI](https://zenodo.org/badge/48195121.svg)](https://zenodo.org/badge/latestdoi/48195121)

## Introduction

**vegas** is a feature rich atomistic magnetic material simulation platform
written in [rust](https://rust-lang.org/). It supports Ising and Heisenberg
spins, as well as a couple of Monte Carlo algorithms, namely Metropolis and
Wolff.

Vegas is meant to be used as a library to build your custom magnetic material
simulation programs. That said, there's an included program that can handle
generic input and can be used as a reference implementation for your own
programs.

## Installation

To install **vegas** you need to have the rust installed. Then, you can install
**vegas** by running the following command:

```bash
cargo install vegas
```

If you want to install **vegas** as a library, you can add it via cargo:

```bash
cargo add vegas
```

## Features

### As a library

- Statistical metrics accumulators.
- Static dispatching of compund Hamiltonians.
- Pre-defined energy components: Gauge, Exchange, Anisotropy, Zeeman.
- Powerful error handling via the `thiserror` crate.
- Flexible instrumentation system, using dynamic dispatching.
- Support for different integration algorithms such as Metropolis.
- Parquet input output support via the `parquet` crate.
- Pre-defined programs: Relax, CoolDown, HysteresisLoop.

### As a command line tool

You can use the toml input file format to run simulations. An example of input
file is given below:

```toml
# Model definition can be Ising or Heisenberg
model = "Ising"
algorithm = "Metropolis"

# You can create unit cells of different lattice types.
[sample.unitcell]
name = "sc"

# You can expand your unit cell to create larger samples.
[sample.size]
x = 10
y = 10
z = 1

# You can set periodic boundary conditions in each direction.
[sample.pbc]
x = true
y = true
z = false


# You can control the steps of the simulation.
[[steps]]
program = "Relax"
steps = 1000
temperature = 4.0

[[steps]]
program = "CoolDown"
max_temperature = 4.0
min_temperature = 0.1
cool_rate = 0.1
relax = 1000
steps = 20000

# You can define outputs to be written during the simulation.
[output]
raw = "./output.parquet"

[output.state]
path = "./state.parquet"
frequency = 1000
```

You can run the simulation by executing the following command:

```bash
vegas run input.toml
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on
GitHub. There are currently some missing features that would benefit the
package, such as:

- Custom exchange interaction values, we currently support only one value.
- More Hamiltonian terms (Dzyaloshinskii-Moriya, Dipolar, etc).
- More integration algorithms (Wolff, Swendsen-Wang, etc).
