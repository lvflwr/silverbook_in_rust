//! Solver for the transport equation using upwind method.

use ndarray::prelude::*;

/// Solver for the transport equation using upwind method.
#[derive(Debug)]
pub struct UpwindSolver {
    u: Array1<f64>,
    v_adv: f64,
    dx: f64,
    dt: f64,
    t_max: f64,
    t: f64,
    step: usize,
    diff_method: DiffMethod,
    completed: bool,
}

impl UpwindSolver {
    /// Create a new `UpwindSolver` instance.
    pub fn new(
        u: Array1<f64>,
        v_adv: f64,
        dx: f64,
        dt: f64,
        t_max: f64,
        diff_method: DiffMethod,
    ) -> Self {
        Self {
            u,
            v_adv,
            dx,
            dt,
            t_max,
            t: 0.0,
            step: 0,
            diff_method,
            completed: false,
        }
    }

    /// Return a reference to the current `u`.
    pub fn borrow_u(&self) -> &Array1<f64> {
        &self.u
    }

    /// Return the current `t`.
    pub fn get_t(&self) -> f64 {
        self.t
    }

    /// Return the current `step`.
    pub fn get_step(&self) -> usize {
        self.step
    }

    /// Return `true` if the calculation has been completed.
    pub fn is_completed(&self) -> bool {
        self.completed
    }

    /// Integrate the transport equation by one time step.
    ///
    /// # Errors
    /// Returns an error if the calculation has already been completed.
    pub fn integrate(&mut self) -> Result<(), &'static str> {
        if self.completed {
            return Err("calculation has already been completed");
        }

        self.u = self
            .diff_method
            .calculate_next_u(&self.u, self.v_adv, self.dx, self.dt);
        self.t += self.dt;
        self.step += 1;

        if self.t >= self.t_max {
            self.completed = true;
        }

        Ok(())
    }
}

/// Difference methods.
#[derive(Debug)]
pub enum DiffMethod {
    /// Forward difference method.
    Forward,
    /// Backward difference method.
    Backward,
}

impl DiffMethod {
    fn calculate_next_u(&self, u: &Array1<f64>, v_adv: f64, dx: f64, dt: f64) -> Array1<f64> {
        match self {
            DiffMethod::Forward => self.calculate_next_u_by_forward(u, v_adv, dx, dt),
            DiffMethod::Backward => self.calculate_next_u_by_backward(u, v_adv, dx, dt),
        }
    }

    fn calculate_next_u_by_forward(
        &self,
        u: &Array1<f64>,
        v_adv: f64,
        dx: f64,
        dt: f64,
    ) -> Array1<f64> {
        u.indexed_iter()
            .map(|(i, _)| {
                if i == 0 || i == u.len() - 1 {
                    u[i]
                } else {
                    u[i] - v_adv * dt / dx * (u[i + 1] - u[i])
                }
            })
            .collect()
    }

    fn calculate_next_u_by_backward(
        &self,
        u: &Array1<f64>,
        v_adv: f64,
        dx: f64,
        dt: f64,
    ) -> Array1<f64> {
        u.indexed_iter()
            .map(|(i, _)| {
                if i == 0 || i == u.len() - 1 {
                    u[i]
                } else {
                    u[i] - v_adv * dt / dx * (u[i] - u[i - 1])
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_upwind_integrate_works() {
        // setup upwind solver and run integrate()
        let u_init = array![1.0, 1.0, 0.0, 0.0, 0.0];
        let mut upwind_solver = UpwindSolver::new(u_init, 1.0, 0.1, 0.1, 0.5, DiffMethod::Backward);
        upwind_solver.integrate().unwrap();

        // check if u, t and step are correctly updated
        let u_exact = array![1.0, 1.0, 1.0, 0.0, 0.0];
        let is_u_correctly_updated = (&upwind_solver.u - u_exact).iter().all(|u| u.abs() < 1e-10);
        assert!(is_u_correctly_updated);
        assert!((upwind_solver.t - 0.1).abs() < 1e-10);
        assert_eq!(upwind_solver.step, 1);
    }
}
