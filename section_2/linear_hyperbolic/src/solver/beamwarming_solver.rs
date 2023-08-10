//! Solver for the transport equation using the Beam-Warming method.
//!
//! # Scheme
//! The Beam-Warming method is given by
//! ```math
//! -\frac{\nu}{2} \lambda u_{j-1}^{n+1} + u_j^{n+1} + \frac{\nu}{2} \lambda u_{j+1}^{n+1} =
//! \frac{\nu}{2} (1 - \lambda) u_{j-1}^n + u_j^n - \frac{\nu}{2} (1 - \lambda) u_{j+1}^n,
//! ```
//! where `\nu = c \frac{\Delta t}{\Delta x}` and `\lambda \in (0, 1)` is the weighting factor.
//!
//! The Beam-Warming method is equivalent to the Crank-Nicolson method when `\lambda = 0.5`,
//! explicit euler method when `\lambda = 0` and implicit euler method when `\lambda = 1`.
//!
//! # Boundary Condition
//! The boundary condition is fixed as
//! ```math
//! u(x_{\pm}, t) = u(x_{\pm}, 0).
//! ```

use super::{NewParams, Solver};
use crate::math::trinomial_eq::TrinomialEq;
use ndarray::prelude::*;
use std::error::Error;

/// Solver for the transport equation using the Beam-Warming method.
#[derive(Debug)]
pub struct BeamwarmingSolver {
    u: Array1<f64>,
    step_max: usize,
    n_cfl: f64,
    lambda: f64,
    trinomial_eq: TrinomialEq,
    step: usize,
    completed: bool,
}

impl BeamwarmingSolver {
    /// Create a new `BeamwarmingSolver` instance.
    pub fn new(new_params: BeamwarmingSolverNewParams) -> Result<Self, &'static str> {
        new_params.validate_new_params()?;

        let u_len = new_params.u.len();

        Ok(Self {
            u: new_params.u,
            step_max: new_params.step_max,
            n_cfl: new_params.n_cfl,
            lambda: new_params.lambda,
            trinomial_eq: TrinomialEq::new(Self::create_mat_coef(
                u_len,
                new_params.n_cfl,
                new_params.lambda,
            )),
            step: 0,
            completed: false,
        })
    }

    fn calculate_u_next(&self) -> Result<Array1<f64>, Box<dyn Error>> {
        let coef_lower_rhs = 0.5 * self.n_cfl * (1.0 - self.lambda);
        let coef_diag_rhs = 1.0;
        let coef_upper_rhs = -coef_lower_rhs;

        let mut u_next: Array1<f64> = (0..self.u.len())
            .map(|i| {
                if i == 0 {
                    return coef_diag_rhs * self.u[i] + coef_upper_rhs * self.u[i + 1];
                }
                if i == self.u.len() - 1 {
                    return coef_lower_rhs * self.u[i - 1] + coef_diag_rhs * self.u[i];
                }

                coef_lower_rhs * self.u[i - 1]
                    + coef_diag_rhs * self.u[i]
                    + coef_upper_rhs * self.u[i + 1]
            })
            .collect();

        self.trinomial_eq.solve(&mut u_next)?;

        Ok(u_next
            .indexed_iter()
            .map(|(i, v)| {
                if i == 0 || i == u_next.len() - 1 {
                    return self.u[i];
                }

                *v
            })
            .collect())
    }

    fn create_mat_coef(n_dim: usize, n_cfl: f64, lambda: f64) -> Array1<(f64, f64, f64)> {
        let coef_lower = -0.5 * n_cfl * lambda;
        let coef_diag = 1.0;
        let coef_upper = -coef_lower;

        Array::from_elem(n_dim, (coef_lower, coef_diag, coef_upper))
    }
}

impl Solver for BeamwarmingSolver {
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

        self.u = self.calculate_u_next()?;
        self.step += 1;

        if self.step >= self.step_max {
            self.completed = true;
        }

        Ok(())
    }
}

/// Parameters for creating a new `BeamwarmingSolver` instance.
pub struct BeamwarmingSolverNewParams {
    pub u: Array1<f64>,
    pub step_max: usize,
    pub n_cfl: f64,
    pub lambda: f64,
}

impl NewParams for BeamwarmingSolverNewParams {
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
        if self.lambda < 0.0 || self.lambda > 1.0 {
            return Err("lambda must be between 0 and 1");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_beamwarming_integrate_works() {
        // setup beamwarming solver and run integrate()
        let u_init = array![1.0, 1.0, 0.0, 0.0, 0.0];
        let new_params = BeamwarmingSolverNewParams {
            u: u_init,
            step_max: 3,
            n_cfl: 1.0,
            lambda: 0.5,
        };
        let mut beamwarming_solver = BeamwarmingSolver::new(new_params).unwrap();
        beamwarming_solver.integrate().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![1.0, 1.22910216718, 0.52631578947, 0.12383900929, 0.0];
        let is_u_correctly_updated = (beamwarming_solver.u - u_exact)
            .iter()
            .all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
        assert_eq!(beamwarming_solver.step, 1);
    }
}
