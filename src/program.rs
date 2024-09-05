//! Programs to run on samples.

use rand::Rng;

use crate::{
    energy::HamiltonianComponent,
    error::{ProgramError, Result},
    integrator::Integrator,
    observables::Sensor,
    state::{Spin, State},
    thermostat::Thermostat,
};

/// A program is a sequence of steps that can be run on a system.
pub trait Program {
    /// Run the program on a system returning the last state.
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        integrator: &I,
        hamiltonian: &H,
        state: State<S>,
    ) -> Result<State<S>>
    where
        S: Spin,
        H: HamiltonianComponent<S>,
        I: Integrator<S, H>,
        R: Rng;
}

/// A program that cools the system to find the Curie temperature.
pub struct CurieTemp {
    max_temp: f64,
    min_temp: f64,
    cool_rate: f64,
    relax: usize,
    steps: usize,
}

impl CurieTemp {
    /// Create a new Curie temperature program.
    pub fn new(max_temp: f64, min_temp: f64, cool_rate: f64, relax: usize, steps: usize) -> Self {
        Self {
            max_temp,
            min_temp,
            cool_rate,
            relax,
            steps,
        }
    }

    /// Set the maximum temperature.
    pub fn set_max_temp(mut self, max_temp: f64) -> Self {
        self.max_temp = max_temp;
        self
    }

    /// Set the minimum temperature.
    pub fn set_min_temp(mut self, min_temp: f64) -> Self {
        self.min_temp = min_temp;
        self
    }

    /// Set the cooling rate.
    pub fn set_cool_rate(mut self, cool_rate: f64) -> Self {
        self.cool_rate = cool_rate;
        self
    }

    /// Set the number of relaxation steps.
    pub fn set_relax(mut self, relax: usize) -> Self {
        self.relax = relax;
        self
    }

    /// Set the number of steps.
    pub fn set_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }
}

impl Default for CurieTemp {
    fn default() -> Self {
        Self::new(3.0, f64::EPSILON, 0.1, 1000, 20000)
    }
}

impl Program for CurieTemp {
    /// Run the program.
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        integrator: &I,
        hamiltonian: &H,
        mut state: State<S>,
    ) -> Result<State<S>>
    where
        I: Integrator<S, H>,
        H: HamiltonianComponent<S>,
        S: Spin,
        R: Rng,
    {
        if self.max_temp < self.min_temp {
            return Err(ProgramError::MaxTempLessThanMinTemp.into());
        }
        if self.steps == 0 {
            return Err(ProgramError::NoSteps.into());
        }
        if self.min_temp < f64::EPSILON {
            return Err(ProgramError::ZeroTemp.into());
        }
        if self.cool_rate < f64::EPSILON {
            return Err(ProgramError::ZeroDelta.into());
        }
        let mut thermostat = Thermostat::new(self.max_temp);
        loop {
            let mut sensor = Sensor::new(thermostat.temp());
            for _ in 0..self.relax {
                state = integrator.step(rng, &thermostat, hamiltonian, state);
            }
            for _ in 0..self.steps {
                state = integrator.step(rng, &thermostat, hamiltonian, state);
                sensor.observe(hamiltonian, &state);
            }
            println!("{}", sensor);
            thermostat.cool(self.cool_rate);
            if thermostat.temp() < self.min_temp {
                break;
            }
        }
        Ok(state)
    }
}

/// A program that relaxes the system.
pub struct Relax {
    steps: usize,
    temp: f64,
}

impl Relax {
    /// Create a new relaxation program.
    pub fn new(steps: usize, temp: f64) -> Self {
        Self { steps, temp }
    }

    /// Set the number of steps.
    pub fn set_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }

    /// Set the temperature.
    pub fn set_temp(mut self, temp: f64) -> Self {
        self.temp = temp;
        self
    }
}

impl Default for Relax {
    fn default() -> Self {
        Self::new(1000, 3.0)
    }
}

impl Program for Relax {
    /// Run the program.
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        integrator: &I,
        hamiltonian: &H,
        mut state: State<S>,
    ) -> Result<State<S>>
    where
        I: Integrator<S, H>,
        H: HamiltonianComponent<S>,
        S: Spin,
        R: Rng,
    {
        if self.steps == 0 {
            return Err(ProgramError::NoSteps.into());
        }
        if self.temp < f64::EPSILON {
            return Err(ProgramError::ZeroTemp.into());
        }
        let thermostat = Thermostat::new(self.temp);
        let mut sensor = Sensor::new(thermostat.temp());
        for _ in 0..self.steps {
            state = integrator.step(rng, &thermostat, hamiltonian, state);
            sensor.observe(hamiltonian, &state);
        }
        Ok(state)
    }
}
