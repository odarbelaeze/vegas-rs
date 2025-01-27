//! Programs to run on samples.

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{
    error::{ProgramError, Result},
    hamiltonian::Hamiltonian,
    integrator::Integrator,
    io::RawIO,
    machine::Machine,
    state::Spin,
};

/// A program is a sequence of steps that can be run on a system.
pub trait Program {
    /// Run the program on a system returning the last state.
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        machine: &mut Machine<H, I, S>,
        output: &mut Option<RawIO>,
    ) -> Result<()>
    where
        S: Spin,
        H: Hamiltonian<S>,
        I: Integrator<S>,
        R: Rng;
}

/// A program that relaxes the system.
#[derive(Debug, Deserialize, Serialize)]
pub struct Relax {
    steps: usize,
    temperature: f64,
}

impl Relax {
    /// Create a new relaxation program.
    pub fn new(steps: usize, temperature: f64) -> Self {
        Self { steps, temperature }
    }

    /// Set the number of steps.
    pub fn set_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }

    /// Set the temperature.
    pub fn set_temperature(mut self, temperature: f64) -> Self {
        self.temperature = temperature;
        self
    }
}

impl Default for Relax {
    fn default() -> Self {
        Self::new(1000, 3.0)
    }
}

impl Program for Relax {
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        machine: &mut Machine<H, I, S>,
        output: &mut Option<RawIO>,
    ) -> Result<()>
    where
        I: Integrator<S>,
        H: Hamiltonian<S>,
        S: Spin,
        R: Rng,
    {
        if self.steps == 0 {
            return Err(ProgramError::NoSteps.into());
        }
        if self.temperature < f64::EPSILON {
            return Err(ProgramError::ZeroTemperature.into());
        }
        machine.set_temperature(self.temperature);
        let sensor = machine.run(rng, self.steps);
        sensor.write(output)?;
        Ok(())
    }
}

/// A program that cools the system to find the Curie temperature.
#[derive(Debug, Deserialize, Serialize)]
pub struct CoolDown {
    max_temperature: f64,
    min_temperature: f64,
    cool_rate: f64,
    relax: usize,
    steps: usize,
}

impl CoolDown {
    /// Create a new Curie temperature program.
    pub fn new(
        max_temperature: f64,
        min_temperature: f64,
        cool_rate: f64,
        relax: usize,
        steps: usize,
    ) -> Self {
        Self {
            max_temperature,
            min_temperature,
            cool_rate,
            relax,
            steps,
        }
    }

    /// Set the maximum temperature.
    pub fn set_max_temperature(mut self, max_temp: f64) -> Self {
        self.max_temperature = max_temp;
        self
    }

    /// Set the minimum temperature.
    pub fn set_min_temperature(mut self, min_temp: f64) -> Self {
        self.min_temperature = min_temp;
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

impl Default for CoolDown {
    fn default() -> Self {
        Self::new(3.0, 0.1, 0.1, 1000, 20000)
    }
}

impl Program for CoolDown {
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        machine: &mut Machine<H, I, S>,
        output: &mut Option<RawIO>,
    ) -> Result<()>
    where
        I: Integrator<S>,
        H: Hamiltonian<S>,
        S: Spin,
        R: Rng,
    {
        if self.max_temperature < self.min_temperature {
            return Err(ProgramError::TemperatureMaxLessThanMin.into());
        }
        if self.steps == 0 {
            return Err(ProgramError::NoSteps.into());
        }
        if self.min_temperature < f64::EPSILON {
            return Err(ProgramError::ZeroTemperature.into());
        }
        if self.cool_rate < f64::EPSILON {
            return Err(ProgramError::ZeroCoolRate.into());
        }
        let mut temperature = self.max_temperature;
        loop {
            machine.set_temperature(temperature);
            let sensor = machine.run(rng, self.relax);
            sensor.write(output)?;
            let sensor = machine.run(rng, self.steps);
            sensor.print();
            sensor.write(output)?;
            temperature -= self.cool_rate;
            if temperature < self.min_temperature {
                break;
            }
        }
        Ok(())
    }
}

/// A program that runs a hysteresis loop.
#[derive(Debug, Deserialize, Serialize)]
pub struct HysteresisLoop {
    steps: usize,
    relax: usize,
    temperature: f64,
    max_field: f64,
    field_step: f64,
}

impl HysteresisLoop {
    /// Create a new hysteresis loop program.
    pub fn new(
        steps: usize,
        relax: usize,
        temperature: f64,
        max_field: f64,
        field_step: f64,
    ) -> Self {
        Self {
            steps,
            relax,
            temperature,
            max_field,
            field_step,
        }
    }

    /// Set the number of steps.
    pub fn set_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }

    /// Set the number of relaxation steps.
    pub fn set_relax(mut self, relax: usize) -> Self {
        self.relax = relax;
        self
    }

    /// Set the temperature.
    pub fn set_temperature(mut self, temperature: f64) -> Self {
        self.temperature = temperature;
        self
    }

    /// Set the maximum field.
    pub fn set_max_field(mut self, max_field: f64) -> Self {
        self.max_field = max_field;
        self
    }

    /// Set the field step.
    pub fn set_field_step(mut self, field_step: f64) -> Self {
        self.field_step = field_step;
        self
    }
}

impl Default for HysteresisLoop {
    fn default() -> Self {
        Self::new(1000, 1000, 3.0, 1.0, 0.1)
    }
}

impl Program for HysteresisLoop {
    fn run<R, I, H, S>(
        &self,
        rng: &mut R,
        machine: &mut Machine<H, I, S>,
        output: &mut Option<RawIO>,
    ) -> Result<()>
    where
        R: Rng,
        I: Integrator<S>,
        H: Hamiltonian<S>,
        S: Spin,
    {
        if self.steps == 0 {
            return Err(ProgramError::NoSteps.into());
        }
        if self.temperature < f64::EPSILON {
            return Err(ProgramError::ZeroTemperature.into());
        }
        if self.max_field < f64::EPSILON {
            return Err(ProgramError::ZeroField.into());
        }
        if self.field_step < f64::EPSILON {
            return Err(ProgramError::ZeroFieldStep.into());
        }
        machine.set_temperature(self.temperature);
        let mut field = 0.0;
        loop {
            machine.set_field(field);
            let sensor = machine.run(rng, self.relax);
            sensor.write(output)?;
            let sensor = machine.run(rng, self.steps);
            sensor.print();
            sensor.write(output)?;
            field += self.field_step;
            if field > self.max_field {
                break;
            }
        }
        loop {
            machine.set_field(field);
            let sensor = machine.run(rng, self.relax);
            sensor.write(output)?;
            let sensor = machine.run(rng, self.steps);
            sensor.print();
            sensor.write(output)?;
            field -= self.field_step;
            if field < -self.max_field {
                break;
            }
        }
        loop {
            machine.set_field(field);
            let sensor = machine.run(rng, self.relax);
            sensor.write(output)?;
            let sensor = machine.run(rng, self.steps);
            sensor.print();
            sensor.write(output)?;
            field += self.field_step;
            if field < self.max_field {
                break;
            }
        }

        Ok(())
    }
}
