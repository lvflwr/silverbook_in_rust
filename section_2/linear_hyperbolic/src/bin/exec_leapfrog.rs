//! Solve the transport equation by the [LeapfrogSolver].
//!
//! # Formulation
//! The transport equation is given by
//! ```math
//! \frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0 (x \in [-1, 1]),
//! ```
//! where `u` is the transported quantity and `c` (`> 0`) is the advection velocity.
//!
//! The initial condition is given by
//! ```math
//! u(x, 0) = 0 (x \ge 0), u(x, 0) = 1 (x < 0).
//! ```
//!
//! For the boundary condition, see [LeapfrogSolver].
//!
//! # Scheme
//! See [LeapfrogSolver].
//!
//! # Input Format
//! Input should be a YAML file in the following format:
//! ```yaml
//! n_x: 20
//! step_max: 6
//! n_cfl: 1.0
//! ncycle_out: 2
//! ```
//!
//! For the meaning of each parameter, see [ExecLeapfrogInputParams].
//!
//! # Output Format
//! See [linear_hyperbolic::output::output].

use linear_hyperbolic::input;
use linear_hyperbolic::input::InputParams;
use linear_hyperbolic::solver::leapfrog_solver::{LeapfrogSolver, LeapfrogSolverNewParams};
use ndarray::prelude::*;
use serde_derive::{Deserialize, Serialize};
use std::fs::{self, File};
use std::process;

/// Solve the transport equation with the given input parameters and output the results to a file.
fn main() {
    // read input parameters
    let mut inputfile = File::open("inputs/section_2/linear_hyperbolic/exec_leapfrog.yml")
        .unwrap_or_else(|err| {
            eprintln!("Problem opening input file: {}", err);
            process::exit(1);
        });
    let input_params: ExecLeapfrogInputParams = input::read_input_params(&mut inputfile)
        .unwrap_or_else(|err| {
            eprintln!("Problem reading input parameters: {}", err);
            process::exit(1);
        });

    // setup output files
    let dir_str = "outputs/section_2/linear_hyperbolic";
    fs::create_dir_all(dir_str).unwrap_or_else(|err| {
        eprintln!("Problem creating output directory: {}", err);
        process::exit(1);
    });
    let mut outputfile =
        File::create(format!("{}/exec_leapfrog.dat", dir_str)).unwrap_or_else(|err| {
            eprintln!("Problem creating output files: {}", err);
            process::exit(1);
        });

    // setup coordinates
    let x: Array1<f64> = Array1::linspace(-1.0, 1.0, input_params.n_x + 1);

    // initialize the solver
    let new_params = LeapfrogSolverNewParams {
        u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
        step_max: input_params.step_max,
        n_cfl: input_params.n_cfl,
    };
    let mut solver = LeapfrogSolver::new(new_params).unwrap_or_else(|err| {
        eprintln!("Problem creating solver: {}", err);
        process::exit(1);
    });

    // run
    linear_hyperbolic::run(&x, &mut solver, &mut outputfile, input_params.ncycle_out)
        .unwrap_or_else(|err| {
            eprintln!("Application error: {}", err);
            process::exit(1);
        });
}

/// Input parameters.
#[derive(Debug, Serialize, Deserialize)]
pub struct ExecLeapfrogInputParams {
    /// Number of cells.
    pub n_x: usize,
    /// Maximum number of time steps.
    pub step_max: usize,
    /// CFL number.
    pub n_cfl: f64,
    /// Number of cycles between outputs.
    pub ncycle_out: usize,
}

impl InputParams for ExecLeapfrogInputParams {
    fn validate_params(&self) -> Result<(), &'static str> {
        if self.n_x == 0 {
            return Err("n_x must be positive");
        }
        if self.step_max == 0 {
            return Err("step_max must be positive");
        }
        if self.n_cfl <= 0.0 {
            return Err("n_cfl must be positive");
        }
        if self.ncycle_out == 0 {
            return Err("ncycle_out must be positive");
        }

        Ok(())
    }
}
