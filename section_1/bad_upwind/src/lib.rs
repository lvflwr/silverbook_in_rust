//! The `bad_upwind` crate provides a means to find out what results can occur when the upwind method is improperly used.
//!
//! Section 1.3 of the book states that when solving the transport equation by the finite difference method, the upwind side
//! of the advection must be used to calculate the value of the next time step in order to get a good answer.
//!
//! We define the upwind method that uses information on the upwind side of the advection as a good upwind method and the one
//! that uses information on the downwind side as a bad upwind method.
//!
//! Both the good and bad upwind methods are implemented in this crate.
//!
//! Using this crate, you can actually compute and see the difference between the good and bad upwind methods.

pub mod input;
pub mod output;
pub mod upwind_solver;

use ndarray::prelude::*;
use std::error::Error;
use std::io::Write;
use upwind_solver::UpwindSolver;

/// Run the solver and output the results.
pub fn run(
    x: &Array1<f64>,
    upwind_solver: &mut UpwindSolver,
    outputstream: &mut impl Write,
    ncycle_out: usize,
) -> Result<(), Box<dyn Error>> {
    // calculate and output
    output::output(outputstream, 0.0, x, upwind_solver.borrow_u())?;
    while !upwind_solver.is_completed() {
        upwind_solver.integrate()?;

        if upwind_solver.get_step() % ncycle_out == 0 {
            output::output(
                outputstream,
                upwind_solver.get_t(),
                x,
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
    use upwind_solver::DiffMethod;

    #[test]
    fn fn_run_works_with_good_upwind_method() {
        // setup input parameters
        let input_params = InputParams {
            v_adv: 1.0,
            n_x: 20,
            t_max: 0.5,
            dt: 0.1,
            ncycle_out: 5,
        };

        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, input_params.n_x + 1);

        // initialize the upwind solver
        let mut upwind_solver = UpwindSolver::new(
            x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            input_params.v_adv,
            x[1] - x[0],
            input_params.dt,
            input_params.t_max,
            DiffMethod::Backward,
        );

        // execute run()
        run(
            &x,
            &mut upwind_solver,
            &mut outputstream,
            input_params.ncycle_out,
        )
        .unwrap();

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
        // setup input parameters
        let input_params = InputParams {
            v_adv: 1.0,
            n_x: 20,
            t_max: 0.5,
            dt: 0.1,
            ncycle_out: 5,
        };

        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, input_params.n_x + 1);

        // initialize the upwind solver
        let mut upwind_solver = UpwindSolver::new(
            x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            input_params.v_adv,
            x[1] - x[0],
            input_params.dt,
            input_params.t_max,
            DiffMethod::Forward,
        );

        // execute run()
        run(
            &x,
            &mut upwind_solver,
            &mut outputstream,
            input_params.ncycle_out,
        )
        .unwrap();

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
