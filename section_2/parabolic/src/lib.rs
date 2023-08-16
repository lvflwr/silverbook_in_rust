//! This crate provides a means of comparing the stability of several different schemes for the diffusion equations.
//!
//! Section 2.3 of the book introduces two schemes for solving the diffusion equation: the FTCS method and the Beam-Warming method,
//! and discusses the stability of these schemes.
//!
//! All of the schemes mentioned in the book are implemented in this crate.
//!
//! Using this crate, you can actually compute and check the stability of each scheme.

pub mod input;
pub mod math;
pub mod output;
pub mod solver;

use ndarray::prelude::*;
use solver::Solver;
use std::error::Error;
use std::io::Write;

/// Run the solver and output the results.
pub fn run(
    x: &Array1<f64>,
    solver: &mut impl Solver,
    outputstream: &mut impl Write,
    ncycle_out: usize,
) -> Result<(), Box<dyn Error>> {
    // calculate and output
    output::output(outputstream, 0, x, solver.borrow_u())?;
    while !solver.is_completed() {
        solver.integrate()?;

        if solver.get_step() % ncycle_out == 0 {
            output::output(outputstream, solver.get_step(), x, solver.borrow_u())?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use solver::beamwarming_solver::{BeamwarmingSolver, BeamwarmingSolverNewParams};
    use solver::ftcs_solver::{FtcsSolver, FtcsSolverNewParams};

    #[test]
    fn fn_run_works_with_ftcs_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = FtcsSolverNewParams {
            u: x.map(|x| if *x < 0.0 { *x + 1.0 } else { -(*x) + 1.0 }),
            step_max: 500,
            mu: 0.5,
        };
        let mut solver = FtcsSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 500).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 0.0000000000
0 -0.9000000000 0.1000000000
0 -0.8000000000 0.2000000000
0 -0.7000000000 0.3000000000
0 -0.6000000000 0.4000000000
0 -0.5000000000 0.5000000000
0 -0.4000000000 0.6000000000
0 -0.3000000000 0.7000000000
0 -0.2000000000 0.8000000000
0 -0.1000000000 0.9000000000
0 0.0000000000 1.0000000000
0 0.1000000000 0.9000000000
0 0.2000000000 0.8000000000
0 0.3000000000 0.7000000000
0 0.4000000000 0.6000000000
0 0.5000000000 0.5000000000
0 0.6000000000 0.4000000000
0 0.7000000000 0.3000000000
0 0.8000000000 0.2000000000
0 0.9000000000 0.1000000000
0 1.0000000000 0.0000000000


500 -1.0000000000 0.0000000000
500 -0.9000000000 0.0002577989
500 -0.8000000000 0.0005155977
500 -0.7000000000 0.0007481615
500 -0.6000000000 0.0009807252
500 -0.5000000000 0.0011652888
500 -0.4000000000 0.0013498524
500 -0.3000000000 0.0014683496
500 -0.2000000000 0.0015868467
500 -0.1000000000 0.0016276780
500 0.0000000000 0.0016685094
500 0.1000000000 0.0016276780
500 0.2000000000 0.0015868467
500 0.3000000000 0.0014683496
500 0.4000000000 0.0013498524
500 0.5000000000 0.0011652888
500 0.6000000000 0.0009807252
500 0.7000000000 0.0007481615
500 0.8000000000 0.0005155977
500 0.9000000000 0.0002577989
500 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_beamwarming_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = BeamwarmingSolverNewParams {
            u: x.map(|x| if *x < 0.0 { *x + 1.0 } else { -(*x) + 1.0 }),
            step_max: 500,
            mu: 0.5,
            lambda: 0.5,
        };
        let mut solver = BeamwarmingSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 500).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 0.0000000000
0 -0.9000000000 0.1000000000
0 -0.8000000000 0.2000000000
0 -0.7000000000 0.3000000000
0 -0.6000000000 0.4000000000
0 -0.5000000000 0.5000000000
0 -0.4000000000 0.6000000000
0 -0.3000000000 0.7000000000
0 -0.2000000000 0.8000000000
0 -0.1000000000 0.9000000000
0 0.0000000000 1.0000000000
0 0.1000000000 0.9000000000
0 0.2000000000 0.8000000000
0 0.3000000000 0.7000000000
0 0.4000000000 0.6000000000
0 0.5000000000 0.5000000000
0 0.6000000000 0.4000000000
0 0.7000000000 0.3000000000
0 0.8000000000 0.2000000000
0 0.9000000000 0.1000000000
0 1.0000000000 0.0000000000


500 -1.0000000000 0.0000000000
500 -0.9000000000 0.0003963585
500 -0.8000000000 0.0007172735
500 -0.7000000000 0.0010212068
500 -0.6000000000 0.0013009629
500 -0.5000000000 0.0015499185
500 -0.4000000000 0.0017621794
500 -0.3000000000 0.0019327205
500 -0.2000000000 0.0020575040
500 -0.1000000000 0.0021335757
500 0.0000000000 0.0021591347
500 0.1000000000 0.0021335757
500 0.2000000000 0.0020575040
500 0.3000000000 0.0019327205
500 0.4000000000 0.0017621794
500 0.5000000000 0.0015499185
500 0.6000000000 0.0013009629
500 0.7000000000 0.0010212068
500 0.8000000000 0.0007172735
500 0.9000000000 0.0003963585
500 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }
}
