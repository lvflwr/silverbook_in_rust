//! Module for solving the trinomial equations.

use ndarray::prelude::*;

/// Solver for the trinomial equations.
#[derive(Debug)]
pub struct TrinomialEq {
    mat_coef: Array1<(f64, f64, f64)>,
}

impl TrinomialEq {
    /// Create a new `TrinomialEq` instance.
    ///
    /// # Arguments
    /// * `mat_coef` - coefficient matrix of the trinomial equation.
    /// The 1st component of each element is the diagonal component of the coefficient matrix
    /// and the 0th and 2nd components are the lower and upper components, respectively.
    pub fn new(mut mat_coef: Array1<(f64, f64, f64)>) -> Self {
        Self::decompose_mat_coef(&mut mat_coef);

        Self { mat_coef }
    }

    /// Solve the trinomial equation.
    ///
    /// # Arguments
    /// * `vec_rhs` - right-hand side vector of the trinomial equation.
    ///
    /// # Examples
    /// ```
    /// use ndarray::prelude::*;
    /// use linear_hyperbolic::math::trinomial_eq::TrinomialEq;
    ///
    /// let mat_coef = array![
    ///   (0.0, 1.0, 2.0),
    ///   (3.0, 4.0, 5.0),
    ///   (6.0, 7.0, 0.0),
    /// ];
    /// let trinomial_eq = TrinomialEq::new(mat_coef);
    /// let mut vec_rhs = array![8.0, 9.0, 10.0];
    /// trinomial_eq.solve(&mut vec_rhs).unwrap();
    ///
    /// let exact_solution = array![21.0 / 22.0, 155.0 / 44.0, -35.0 / 22.0];
    /// let is_correctly_solved = (&vec_rhs - exact_solution).iter().all(|x| x.abs() < 1e-10);
    /// assert!(is_correctly_solved);
    /// ```
    ///
    /// # Errors
    /// Returns an error if the length of `vec_rhs` is not equal to the length of `mat_coef`.
    pub fn solve(&self, vec_rhs: &mut Array1<f64>) -> Result<(), &'static str> {
        if vec_rhs.len() != self.mat_coef.len() {
            return Err("The length of vec_rhs must be equal to the length of mat_coef");
        }

        // Forward elimination
        for i in 1..vec_rhs.len() {
            vec_rhs[i] -= self.mat_coef[i].0 * vec_rhs[i - 1];
        }

        // Back substitution
        for i in (0..vec_rhs.len()).rev() {
            if i == vec_rhs.len() - 1 {
                vec_rhs[i] /= self.mat_coef[i].1;
                continue;
            }

            vec_rhs[i] = (vec_rhs[i] - self.mat_coef[i].2 * vec_rhs[i + 1]) / self.mat_coef[i].1;
        }

        Ok(())
    }

    fn decompose_mat_coef(mat_coef: &mut Array1<(f64, f64, f64)>) {
        // Forward elimination
        for i in 1..mat_coef.len() {
            mat_coef[i].0 /= mat_coef[i - 1].1;
            mat_coef[i].1 -= mat_coef[i].0 * mat_coef[i - 1].2;
        }
    }
}
