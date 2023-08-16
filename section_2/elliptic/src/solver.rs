//! Solvers for the diffusion equation.

pub mod point_jacobi_solver;
pub mod sor_solver;

use ndarray::prelude::*;
use std::error::Error;

/// Solver for the diffusion equation.
pub trait Solver {
    /// Execute solving the diffusion equation.
    fn exec(&mut self) -> Result<(), Box<dyn Error>>;
    /// Return a reference to `u`.
    fn borrow_u(&self) -> &Array2<f64>;
    /// Return the number of iterations.
    fn get_n_iter(&self) -> usize;
}

/// Parameters for creating a new solver.
pub trait NewParams {
    /// Validate the parameters for creating a new solver.
    fn validate_new_params(&self) -> Result<(), &'static str>;
}
