//! Solver for the transport equation using the upwind method.
//!
//! # Scheme
//! The upwind method is given by
//! ```math
//! u_j^{n+1} = u_j^n - \nu (u_j^n - u_{j-1}^n),
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

/// Solver for the transport equation using the upwind method.
#[derive(Debug)]
pub struct UpwindSolver {
    u: Array1<f64>,
    step_max: usize,
    n_cfl: f64,
    step: usize,
    completed: bool,
}

impl UpwindSolver {
    /// Create a new `UpwindSolver` instance.
    pub fn new(new_params: UpwindSolverNewParams) -> Result<Self, &'static str> {
        new_params.validate_new_params()?;

        Ok(Self {
            u: new_params.u,
            step_max: new_params.step_max,
            n_cfl: new_params.n_cfl,
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

                self.u[i] - self.n_cfl * (self.u[i] - self.u[i - 1])
            })
            .collect()
    }
}

impl Solver for UpwindSolver {
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

        self.u = self.calculate_u_next();
        self.step += 1;

        if self.step >= self.step_max {
            self.completed = true;
        }

        Ok(())
    }
}

/// Parameters for creating a new `UpwindSolver` instance.
pub struct UpwindSolverNewParams {
    pub u: Array1<f64>,
    pub step_max: usize,
    pub n_cfl: f64,
}

impl NewParams for UpwindSolverNewParams {
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
    fn fn_upwind_integrate_works() {
        // setup upwind solver and run integrate()
        let u_init = array![1.0, 1.0, 0.0, 0.0, 0.0];
        let new_params = UpwindSolverNewParams {
            u: u_init,
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut upwind_solver = UpwindSolver::new(new_params).unwrap();
        upwind_solver.integrate().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![1.0, 1.0, 0.5, 0.0, 0.0];
        let is_u_correctly_updated = (upwind_solver.u - u_exact).iter().all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
        assert_eq!(upwind_solver.step, 1);
    }
}
