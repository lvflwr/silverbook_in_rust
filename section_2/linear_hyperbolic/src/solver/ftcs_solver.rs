//! Solver for the transport equation using the FTCS method.
//!
//! # Scheme
//! The FTCS (Forward in Time and Central Difference in Space) method is given by
//! ```math
//! u_j^{n+1} = u_j^n - \frac{1}{2} \nu (u_{j+1}^n - u_{j-1}^n),
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

/// Solver for the transport equation using the FTCS method.
#[derive(Debug)]
pub struct FtcsSolver {
    u: Array1<f64>,
    step_max: usize,
    n_cfl: f64,
    step: usize,
    completed: bool,
}

impl FtcsSolver {
    /// Create a new `FtcsSolver` instance.
    pub fn new(new_params: FtcsSolverNewParams) -> Result<Self, &'static str> {
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

                self.u[i] - 0.5 * self.n_cfl * (self.u[i + 1] - self.u[i - 1])
            })
            .collect()
    }
}

impl Solver for FtcsSolver {
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

/// Parameters for creating a new `FtcsSolver` instance.
pub struct FtcsSolverNewParams {
    /// Initial value of `u`.
    pub u: Array1<f64>,
    /// Maximum number of time steps.
    pub step_max: usize,
    /// CFL number.
    pub n_cfl: f64,
}

impl NewParams for FtcsSolverNewParams {
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
    fn fn_ftcs_integrate_works() {
        // setup ftcs solver and run integrate()
        let u_init = array![1.0, 1.0, 0.0, 0.0, 0.0];
        let new_params = FtcsSolverNewParams {
            u: u_init,
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut ftcs_solver = FtcsSolver::new(new_params).unwrap();
        ftcs_solver.integrate().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![1.0, 1.25, 0.25, 0.0, 0.0];
        let is_u_correctly_updated = (ftcs_solver.u - u_exact).iter().all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
        assert_eq!(ftcs_solver.step, 1);
    }
}
