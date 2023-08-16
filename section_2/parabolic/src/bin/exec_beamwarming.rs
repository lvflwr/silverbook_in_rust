//! Solve the diffusion equation by the [parabolic::solver::beamwarming_solver].
//!
//! # Formulation
//! The diffusion equation is given by
//! ```math
//! \frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2} (x \in [-1, 1]),
//! ```
//! where `u` is the diffusion quantity and `\alpha` is the diffusion coefficient.
//!
//! The initial condition is given by
//! ```math
//! u(x, 0) = -x + 1 (x \ge 0), u(x, 0) = x + 1 (x < 0).
//! ```
//!
//! For the boundary condition, see [parabolic::solver::beamwarming_solver].
//!
//! # Scheme
//! See [parabolic::solver::beamwarming_solver].
//!
//! # Input Format
//! Input should be a YAML file in the following format:
//! ```yaml
//! n_x: 100
//! step_max: 10000
//! mu: 0.5
//! lambda: 0.5
//! ncycle_out: 1000
//! ```
//!
//! For the meaning of each parameter, see [ExecBeamwarmingInputParams].
//!
//! # Output Format
//! See [parabolic::output::output].

use ndarray::prelude::*;
use parabolic::input;
use parabolic::input::InputParams;
use parabolic::solver::beamwarming_solver::{BeamwarmingSolver, BeamwarmingSolverNewParams};
use serde_derive::{Deserialize, Serialize};
use std::fs::{self, File};
use std::process;

/// Solve the diffusion equation with the given input parameters and output the results to a file.
fn main() {
    // read input parameters
    let mut inputfile = File::open("inputs/section_2/parabolic/exec_beamwarming.yml")
        .unwrap_or_else(|err| {
            eprintln!("Problem opening input file: {}", err);
            process::exit(1);
        });
    let input_params: ExecBeamwarmingInputParams = input::read_input_params(&mut inputfile)
        .unwrap_or_else(|err| {
            eprintln!("Problem reading input parameters: {}", err);
            process::exit(1);
        });

    // setup output files
    let dir_str = "outputs/section_2/parabolic";
    fs::create_dir_all(dir_str).unwrap_or_else(|err| {
        eprintln!("Problem creating output directory: {}", err);
        process::exit(1);
    });
    let mut outputfile =
        File::create(format!("{}/exec_beamwarming.dat", dir_str)).unwrap_or_else(|err| {
            eprintln!("Problem creating output files: {}", err);
            process::exit(1);
        });

    // setup coordinates
    let x: Array1<f64> = Array1::linspace(-1.0, 1.0, input_params.n_x + 1);

    // initialize the solver
    let new_params = BeamwarmingSolverNewParams {
        u: x.map(|x| if *x < 0.0 { *x + 1.0 } else { -(*x) + 1.0 }),
        step_max: input_params.step_max,
        mu: input_params.mu,
        lambda: input_params.lambda,
    };
    let mut solver = BeamwarmingSolver::new(new_params).unwrap_or_else(|err| {
        eprintln!("Problem creating solver: {}", err);
        process::exit(1);
    });

    // run
    parabolic::run(&x, &mut solver, &mut outputfile, input_params.ncycle_out).unwrap_or_else(
        |err| {
            eprintln!("Application error: {}", err);
            process::exit(1);
        },
    );
}

/// Input parameters.
#[derive(Debug, Serialize, Deserialize)]
pub struct ExecBeamwarmingInputParams {
    /// Number of cells.
    pub n_x: usize,
    /// Maximum number of time steps.
    pub step_max: usize,
    /// diffusion coefficient * dt / dx^2.
    pub mu: f64,
    /// Weighting factor in differencing scheme.
    pub lambda: f64,
    /// Number of cycles between outputs.
    pub ncycle_out: usize,
}

impl InputParams for ExecBeamwarmingInputParams {
    fn validate_params(&self) -> Result<(), &'static str> {
        if self.n_x == 0 {
            return Err("n_x must be positive");
        }
        if self.step_max == 0 {
            return Err("step_max must be positive");
        }
        if self.mu <= 0.0 {
            return Err("mu must be positive");
        }
        if self.lambda < 0.0 || self.lambda > 1.0 {
            return Err("lambda must be between 0 and 1");
        }
        if self.ncycle_out == 0 {
            return Err("ncycle_out must be positive");
        }

        Ok(())
    }
}
