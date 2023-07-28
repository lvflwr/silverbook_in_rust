//! The `bad_upwind` crate provides a means to find out what results can occur when the upwind method is improperly used.
//!
//! Section 1.3 of the book states that when solving the transport equation by the finite difference method, the upwind side
//! of the advection must be used to calculate the value of the next time step in order to get a good answer.
//!
//! We define the upwind method that uses information on the upwind side of the advection as a good upwind method and the one
//! that uses information on the downwind side as a bad upwind method.
//!
//! And we calculate the solution of the transport equation by the good and bad upwind methods, respectively, to see the difference
//! in the results.
//!
//! # Formulation and Scheme
//! See [run].
//!
//! # Input Format
//! See [input::read_input_params].
//!
//! # Output Format
//! See [output::output].

pub mod input;
pub mod output;
pub mod upwind_solver;

use input::InputParams;
use ndarray::prelude::*;
use std::error::Error;
use std::io::Write;
use upwind_solver::DiffMethod;
use upwind_solver::UpwindSolver;

/// Solve the transport equation by the upwind method.
///
/// # Formulation
/// The transport equation is given by
/// ```math
/// \frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0 (x \in [-1, 1])),
/// ```
/// where `u` is the transported quantity and `c` (`> 0`) is the advection velocity.
///
/// The initial condition is given by
/// ```math
/// u(x, 0) = 0 (x \ge 0), u(x, 0) = 1 (x < 0).
/// ```
///
/// The boundary condition is given by
/// ```math
/// u(-1, t) = 1, u(1, t) = 0.
/// ```
///
/// # Scheme
/// The upwind method is used to solve the transport equation.
///
/// The good upwind method, in this case, [DiffMethod::Forward], is given by
/// ```math
/// u_j^{n+1} = u_j^n -  c \frac{\Delta t}{\Delta x} (u_j^n - u_{j-1}^n).
/// ```
///
/// The bad upwind method, in this case, [DiffMethod::Forward], is given by
/// ```math
/// u_j^{n+1} = u_j^n -  c \frac{\Delta t}{\Delta x} (u_{j+1}^n - u_{j}^n).
/// ```
pub fn run(
    input_params: &InputParams,
    diff_method: DiffMethod,
    outputstream: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    // setup coordinates
    let x: Array1<f64> = Array1::linspace(-1.0, 1.0, input_params.n_x + 1);
    let dx = x[1] - x[0];

    // initialize the upwind solver
    let mut upwind_solver = UpwindSolver::new(
        x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
        input_params.v_adv,
        dx,
        input_params.dt,
        input_params.t_max,
        diff_method,
    );

    // calculate and output
    output::output(outputstream, 0.0, &x, upwind_solver.borrow_u())?;
    while !upwind_solver.is_completed() {
        upwind_solver.integrate()?;

        if upwind_solver.get_step() % input_params.ncycle_out == 0 {
            output::output(
                outputstream,
                upwind_solver.get_t(),
                &x,
                upwind_solver.borrow_u(),
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use input::InputParams;

    #[test]
    fn fn_run_works_with_good_upwind_method() {
        // setup input parameters and output stream, and execute run()
        let input_params = InputParams {
            v_adv: 1.0,
            n_x: 20,
            t_max: 0.5,
            dt: 0.1,
            ncycle_out: 5,
        };
        let mut outputstream: Vec<u8> = Vec::new();
        run(&input_params, DiffMethod::Backward, &mut outputstream).unwrap();

        // check if the output is correct
        let output_expected = "\
0.00 -1.0000000000 1.0000000000
0.00 -0.9000000000 1.0000000000
0.00 -0.8000000000 1.0000000000
0.00 -0.7000000000 1.0000000000
0.00 -0.6000000000 1.0000000000
0.00 -0.5000000000 1.0000000000
0.00 -0.4000000000 1.0000000000
0.00 -0.3000000000 1.0000000000
0.00 -0.2000000000 1.0000000000
0.00 -0.1000000000 1.0000000000
0.00 0.0000000000 0.0000000000
0.00 0.1000000000 0.0000000000
0.00 0.2000000000 0.0000000000
0.00 0.3000000000 0.0000000000
0.00 0.4000000000 0.0000000000
0.00 0.5000000000 0.0000000000
0.00 0.6000000000 0.0000000000
0.00 0.7000000000 0.0000000000
0.00 0.8000000000 0.0000000000
0.00 0.9000000000 0.0000000000
0.00 1.0000000000 0.0000000000


0.50 -1.0000000000 1.0000000000
0.50 -0.9000000000 1.0000000000
0.50 -0.8000000000 1.0000000000
0.50 -0.7000000000 1.0000000000
0.50 -0.6000000000 1.0000000000
0.50 -0.5000000000 1.0000000000
0.50 -0.4000000000 1.0000000000
0.50 -0.3000000000 1.0000000000
0.50 -0.2000000000 1.0000000000
0.50 -0.1000000000 1.0000000000
0.50 0.0000000000 1.0000000000
0.50 0.1000000000 1.0000000000
0.50 0.2000000000 1.0000000000
0.50 0.3000000000 1.0000000000
0.50 0.4000000000 1.0000000000
0.50 0.5000000000 0.0000000000
0.50 0.6000000000 0.0000000000
0.50 0.7000000000 0.0000000000
0.50 0.8000000000 0.0000000000
0.50 0.9000000000 0.0000000000
0.50 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_bad_upwind_method() {
        // setup input parameters and output stream, and execute run()
        let input_params = InputParams {
            v_adv: 1.0,
            n_x: 20,
            t_max: 0.5,
            dt: 0.1,
            ncycle_out: 5,
        };
        let mut outputstream: Vec<u8> = Vec::new();
        run(&input_params, DiffMethod::Forward, &mut outputstream).unwrap();

        // check if the output is correct
        let output_expected = "\
0.00 -1.0000000000 1.0000000000
0.00 -0.9000000000 1.0000000000
0.00 -0.8000000000 1.0000000000
0.00 -0.7000000000 1.0000000000
0.00 -0.6000000000 1.0000000000
0.00 -0.5000000000 1.0000000000
0.00 -0.4000000000 1.0000000000
0.00 -0.3000000000 1.0000000000
0.00 -0.2000000000 1.0000000000
0.00 -0.1000000000 1.0000000000
0.00 0.0000000000 0.0000000000
0.00 0.1000000000 0.0000000000
0.00 0.2000000000 0.0000000000
0.00 0.3000000000 0.0000000000
0.00 0.4000000000 0.0000000000
0.00 0.5000000000 0.0000000000
0.00 0.6000000000 0.0000000000
0.00 0.7000000000 0.0000000000
0.00 0.8000000000 0.0000000000
0.00 0.9000000000 0.0000000000
0.00 1.0000000000 0.0000000000


0.50 -1.0000000000 1.0000000000
0.50 -0.9000000000 1.0000000000
0.50 -0.8000000000 1.0000000000
0.50 -0.7000000000 1.0000000000
0.50 -0.6000000000 1.0000000000
0.50 -0.5000000000 2.0000000000
0.50 -0.4000000000 -8.0000000000
0.50 -0.3000000000 32.0000000000
0.50 -0.2000000000 -48.0000000000
0.50 -0.1000000000 32.0000000000
0.50 0.0000000000 0.0000000000
0.50 0.1000000000 0.0000000000
0.50 0.2000000000 0.0000000000
0.50 0.3000000000 0.0000000000
0.50 0.4000000000 0.0000000000
0.50 0.5000000000 0.0000000000
0.50 0.6000000000 0.0000000000
0.50 0.7000000000 0.0000000000
0.50 0.8000000000 0.0000000000
0.50 0.9000000000 0.0000000000
0.50 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }
}
