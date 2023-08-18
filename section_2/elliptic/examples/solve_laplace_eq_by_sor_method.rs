//! Solve the diffusion equation by the [elliptic::solver::sor_solver].
//!
//! # Formulation
//! The diffusion equation is given by
//! ```math
//! \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0,
//! ```
//! where `u` is the diffusion quantity.
//!
//! The boundary condition is given by
//! ```math
//! u(x, y) = 1 (y = y_{+}), u(x, y) = 0 (x = x_{\pm} or y = y_{-}).
//! ```
//! See also [elliptic::solver::sor_solver] for the boundary condition.
//!
//! # Scheme
//! See [elliptic::solver::sor_solver].
//!
//! # Input Format
//! Input should be a YAML file in the following format:
//! ```yaml
//! n_x: 20
//! n_y: 20
//! n_iter_max: 10000
//! omega: 1.5
//! ```
//!
//! For the meaning of each parameter, see [ExecSorInputParams].
//!
//! # Output Format
//! See [elliptic::output::output].

use elliptic::input;
use elliptic::input::InputParams;
use elliptic::solver::sor_solver::{SorSolver, SorSolverNewParams};
use ndarray::prelude::*;
use serde_derive::{Deserialize, Serialize};
use std::fs::{self, File};
use std::process;

/// Solve the diffusion equation with the given input parameters and output the results to a file.
fn main() {
    // read input parameters
    let mut inputfile =
        File::open("inputs/section_2/elliptic/solve_laplace_eq_by_sor_method/input.yml")
            .unwrap_or_else(|err| {
                eprintln!("Problem opening input file: {}", err);
                process::exit(1);
            });
    let input_params: ExecSorInputParams =
        input::read_input_params(&mut inputfile).unwrap_or_else(|err| {
            eprintln!("Problem reading input parameters: {}", err);
            process::exit(1);
        });

    // setup output files
    let dir_str = "outputs/section_2/elliptic/solve_laplace_eq_by_sor_method";
    fs::create_dir_all(dir_str).unwrap_or_else(|err| {
        eprintln!("Problem creating output directory: {}", err);
        process::exit(1);
    });
    let mut outputfile = File::create(format!("{}/solution.dat", dir_str)).unwrap_or_else(|err| {
        eprintln!("Problem creating output files: {}", err);
        process::exit(1);
    });

    // setup initial and boundary conditions
    let mut u_init: Array2<f64> = Array::zeros((input_params.n_x + 1, input_params.n_y + 1));
    u_init
        .slice_mut(s![.., input_params.n_y])
        .assign(&Array::ones(input_params.n_x + 1));

    // initialize the solver
    let new_params = SorSolverNewParams {
        u_init,
        n_iter_max: input_params.n_iter_max,
        omega: input_params.omega,
    };
    let mut solver = SorSolver::new(new_params).unwrap_or_else(|err| {
        eprintln!("Problem creating solver: {}", err);
        process::exit(1);
    });

    // run
    elliptic::run(&mut solver, &mut outputfile).unwrap_or_else(|err| {
        eprintln!("Application error: {}", err);
        process::exit(1);
    });
}

/// Input parameters.
#[derive(Debug, Serialize, Deserialize)]
pub struct ExecSorInputParams {
    /// Number of grids in x direction.
    pub n_x: usize,
    /// Number of grids in y direction.
    pub n_y: usize,
    /// Maximum number of iterations.
    pub n_iter_max: usize,
    /// Relaxation parameter.
    pub omega: f64,
}

impl InputParams for ExecSorInputParams {
    fn validate_params(&self) -> Result<(), &'static str> {
        if self.n_x == 0 {
            return Err("n_x must be positive");
        }
        if self.n_y == 0 {
            return Err("n_y must be positive");
        }
        if self.n_iter_max == 0 {
            return Err("n_iter_max must be positive");
        }
        if self.omega < 1.0 || self.omega > 2.0 {
            return Err("omega must be between 1 and 2");
        }

        Ok(())
    }
}
