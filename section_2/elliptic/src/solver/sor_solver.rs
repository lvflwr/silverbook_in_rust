//! Solver for the diffusion equation using the SOR method.
//!
//! # Scheme
//! The SOR method is given by
//! ```math
//! u_{j,k}^{n+1} = u_{j,k}^n + \frac{1}{4} \omega (u_{j-1,k}^{n+1} + u_{j+1,k}^n - u_{j,k}^n + u_{j,k-1}^{n+1} + u_{j,k+1}^n),
//! ```
//! where `\omega \in [1, 2]` is the relaxation parameter.
//!
//! # Boundary Condition
//! The boundary condition is fixed as
//! ```math
//! u(x_{\pm}, y_{\pm}) = u_init(x_{\pm}, y_{\pm}).
//! ```

use super::{NewParams, Solver};
use ndarray::prelude::*;
use std::error::Error;

/// Solver for the diffusion equation using the SOR method.
#[derive(Debug)]
pub struct SorSolver {
    u: Array2<f64>,
    n_iter_max: usize,
    omega: f64,
    epsilon: f64,
    n_iter: usize,
    executed: bool,
    converged: bool,
}

impl SorSolver {
    /// Create a new `SorSolver` instance.
    pub fn new(new_params: SorSolverNewParams) -> Result<Self, &'static str> {
        new_params.validate_new_params()?;

        Ok(Self {
            u: new_params.u_init,
            n_iter_max: new_params.n_iter_max,
            omega: new_params.omega,
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

                u_next[[i_x, i_y]] = (1.0 - self.omega) * u_next[[i_x, i_y]]
                    + 0.25
                        * self.omega
                        * (u_next[[i_x - 1, i_y]]
                            + u_next[[i_x + 1, i_y]]
                            + u_next[[i_x, i_y - 1]]
                            + u_next[[i_x, i_y + 1]]);
            }
        }

        u_next
    }
}

impl Solver for SorSolver {
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

/// Parameters for creating a new `SorSolver` instance.
pub struct SorSolverNewParams {
    /// Initial values of `u`.
    pub u_init: Array2<f64>,
    /// Maximum number of iterations.
    pub n_iter_max: usize,
    /// Relaxation parameter.
    pub omega: f64,
}

impl NewParams for SorSolverNewParams {
    fn validate_new_params(&self) -> Result<(), &'static str> {
        if self.u_init.is_empty() {
            return Err("u must not be empty");
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_sor_exec_works() {
        // setup ftcs solver and run integrate()
        let u_init = array![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let new_params = SorSolverNewParams {
            u_init,
            n_iter_max: 100,
            omega: 1.5,
        };
        let mut solver = SorSolver::new(new_params).unwrap();
        solver.exec().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.12500000000, 0.37499999998, 1.0],
            [0.0, 0.12500000000, 0.37499999998, 1.0],
            [0.0, 0.0, 0.0, 1.0]
        ];
        let is_u_correctly_updated = (solver.u - u_exact).iter().all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
    }
}
