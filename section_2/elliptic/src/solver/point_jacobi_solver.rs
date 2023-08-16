//! Solver for the diffusion equation using the Point Jacobi method.
//!
//! # Scheme
//! The Point Jacobi method is given by
//! ```math
//! u_{j,k}^{n+1} = \frac{1}{4} (u_{j-1,k}^n + u_{j+1,k}^n + u_{j,k-1}^n + u_{j,k+1}^n).
//! ```
//!
//! # Boundary Condition
//! The boundary condition is fixed as
//! ```math
//! u(x_{\pm}, y_{\pm}) = u_init(x_{\pm}, y_{\pm}).
//! ```

use super::{NewParams, Solver};
use ndarray::prelude::*;
use std::error::Error;

/// Solver for the diffusion equation using the Point Jacobi method.
#[derive(Debug)]
pub struct PointJacobiSolver {
    u: Array2<f64>,
    n_iter_max: usize,
    epsilon: f64,
    n_iter: usize,
    executed: bool,
    converged: bool,
}

impl PointJacobiSolver {
    /// Create a new `PointJacobiSolver` instance.
    pub fn new(new_params: PointJacobiSolverNewParams) -> Result<Self, &'static str> {
        new_params.validate_new_params()?;

        Ok(Self {
            u: new_params.u_init,
            n_iter_max: new_params.n_iter_max,
            epsilon: 1.0e-10,
            n_iter: 0,
            executed: false,
            converged: false,
        })
    }

    fn iterate(&mut self) {
        let u_next = self.calculate_u_next();

        self.converged = (&u_next - &self.u).iter().all(|u| u.abs() <= self.epsilon);
        self.u = u_next;
        self.n_iter += 1;
    }

    fn calculate_u_next(&self) -> Array2<f64> {
        let mut u_next = self.u.clone();
        for i_x in 1..self.u.shape()[0] - 1 {
            for i_y in 1..self.u.shape()[1] - 1 {
                if i_x == 0
                    || i_x == self.u.shape()[0] - 1
                    || i_y == 0
                    || i_y == self.u.shape()[1] - 1
                {
                    continue;
                }

                u_next[[i_x, i_y]] = 0.25
                    * (self.u[[i_x - 1, i_y]]
                        + self.u[[i_x + 1, i_y]]
                        + self.u[[i_x, i_y - 1]]
                        + self.u[[i_x, i_y + 1]]);
            }
        }

        u_next
    }
}

impl Solver for PointJacobiSolver {
    fn exec(&mut self) -> Result<(), Box<dyn Error>> {
        if self.executed {
            return Err(Box::<dyn Error>::from("solver has already been executed"));
        }
        self.executed = true;

        while !self.converged {
            if self.n_iter >= self.n_iter_max {
                return Err(Box::<dyn Error>::from(
                    "maximum number of iterations reached",
                ));
            }

            self.iterate();
        }

        Ok(())
    }

    fn borrow_u(&self) -> &Array2<f64> {
        &self.u
    }

    fn get_n_iter(&self) -> usize {
        self.n_iter
    }
}

/// Parameters for creating a new `PointJacobiSolver` instance.
pub struct PointJacobiSolverNewParams {
    /// Initial values of `u`.
    pub u_init: Array2<f64>,
    /// Maximum number of iterations.
    pub n_iter_max: usize,
}

impl NewParams for PointJacobiSolverNewParams {
    fn validate_new_params(&self) -> Result<(), &'static str> {
        if self.u_init.is_empty() {
            return Err("u must not be empty");
        }
        if self.n_iter_max == 0 {
            return Err("n_iter_max must be positive");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_point_jacobi_exec_works() {
        // setup ftcs solver and run integrate()
        let u_init = array![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let new_params = PointJacobiSolverNewParams {
            u_init,
            n_iter_max: 100,
        };
        let mut solver = PointJacobiSolver::new(new_params).unwrap();
        solver.exec().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.12499999994, 0.37499999994, 1.0],
            [0.0, 0.12499999994, 0.37499999994, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let is_u_correctly_updated = (solver.u - u_exact).iter().all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
    }
}
