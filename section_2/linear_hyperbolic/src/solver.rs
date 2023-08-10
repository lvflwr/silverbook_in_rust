//! Solvers for the transport equation.

pub mod beamwarming_solver;
pub mod ftcs_solver;
pub mod lax_solver;
pub mod laxwendroff_solver;
pub mod leapfrog_solver;
pub mod maccormack_solver;
pub mod upwind_solver;

use ndarray::prelude::*;
use std::error::Error;

/// Solver for the transport equation.
pub trait Solver {
    /// Return a reference to the current `u`.
    fn borrow_u(&self) -> &Array1<f64>;
    /// Return the current `step`.
    fn get_step(&self) -> usize;
    /// Return `true` if the calculation has been completed.
    fn is_completed(&self) -> bool;
    /// Integrate the transport equation by one step.
    fn integrate(&mut self) -> Result<(), Box<dyn Error>>;
}

/// Parameters for creating a new solver.
pub trait NewParams {
    /// Validate the parameters for creating a new solver.
    fn validate_new_params(&self) -> Result<(), &'static str>;
}
