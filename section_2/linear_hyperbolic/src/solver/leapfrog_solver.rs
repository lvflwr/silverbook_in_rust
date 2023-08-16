//! Solver for the transport equation using the Leap-Frog method.
//!
//! # Scheme
//! The Leap-Frog method is given by
//! ```math
//! u_j^{n+1} = u_j^{n-1} - \nu (u_{j+1}^n - u_{j-1}^n),
//! ```
//! where `\nu = c \frac{\Delta t}{\Delta x}`.
//!
//! # Boundary Condition
//! The boundary condition is fixed as
//! ```math
//! u(x_{\pm}, t) = u(x_{\pm}, 0).
//! ```

use super::{NewParams, Solver};
use ndarray::prelude::*;
use std::error::Error;

/// Solver for the transport equation using the Leap-Frog method.
#[derive(Debug)]
pub struct LeapfrogSolver {
    u: Array1<f64>,
    step_max: usize,
    n_cfl: f64,
    u_prev: Array1<f64>,
    step: usize,
    completed: bool,
}

impl LeapfrogSolver {
    /// Create a new `LeapfrogSolver` instance.
    pub fn new(new_params: LeapfrogSolverNewParams) -> Result<Self, &'static str> {
        new_params.validate_new_params()?;

        Ok(Self {
            u: new_params.u.clone(),
            step_max: new_params.step_max,
            n_cfl: new_params.n_cfl,
            u_prev: new_params.u,
            step: 0,
            completed: false,
        })
    }

    fn calculate_u_next(&self) -> Array1<f64> {
        self.u
            .indexed_iter()
            .map(|(i, _)| {
                if i == 0 || i == self.u.len() - 1 {
                    return self.u[i];
                }

                self.u_prev[i] - 0.5 * self.n_cfl * (self.u[i + 1] - self.u[i - 1])
            })
            .collect()
    }
}

impl Solver for LeapfrogSolver {
    fn borrow_u(&self) -> &Array1<f64> {
        &self.u
    }

    fn get_step(&self) -> usize {
        self.step
    }

    fn is_completed(&self) -> bool {
        self.completed
    }

    fn integrate(&mut self) -> Result<(), Box<dyn Error>> {
        if self.completed {
            return Err(Box::<dyn Error>::from(
                "calculation has already been completed",
            ));
        }

        let next_u = self.calculate_u_next();
        self.u_prev = self.u.clone();
        self.u = next_u;
        self.step += 1;

        if self.step >= self.step_max {
            self.completed = true;
        }

        Ok(())
    }
}

/// Parameters for creating a new `LeapfrogSolver` instance.
pub struct LeapfrogSolverNewParams {
    /// Initial value of `u`.
    pub u: Array1<f64>,
    /// Maximum number of time steps.
    pub step_max: usize,
    /// CFL number.
    pub n_cfl: f64,
}

impl NewParams for LeapfrogSolverNewParams {
    fn validate_new_params(&self) -> Result<(), &'static str> {
        if self.u.is_empty() {
            return Err("u must not be empty");
        }
        if self.step_max == 0 {
            return Err("step_max must be positive");
        }
        if self.n_cfl <= 0.0 {
            return Err("n_cfl must be positive");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_leapfrog_integrate_works() {
        // setup leapfrog solver and run integrate()
        let u_init = array![1.0, 1.0, 0.0, 0.0, 0.0];
        let new_params = LeapfrogSolverNewParams {
            u: u_init,
            step_max: 6,
            n_cfl: 1.0,
        };
        let mut leapfrog_solver = LeapfrogSolver::new(new_params).unwrap();
        leapfrog_solver.integrate().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![1.0, 1.5, 0.5, 0.0, 0.0];
        let is_u_correctly_updated = (leapfrog_solver.u - u_exact)
            .iter()
            .all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
        assert_eq!(leapfrog_solver.step, 1);
    }
}
