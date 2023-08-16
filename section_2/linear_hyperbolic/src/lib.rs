//! This crate provides a means of comparing the stability of several different schemes for the linear wave equations.
//!
//! Section 2.2 of the book introduces several schemes for the linear wave equations and discusses the stability of each.
//!
//! All of the schemes mentioned in the book are implemented in this crate.
//!
//! Using this crate, you can actually compute and see how the dissipative and dispersive errors arise for each scheme.

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
    use solver::lax_solver::{LaxSolver, LaxSolverNewParams};
    use solver::laxwendroff_solver::{LaxwendroffSolver, LaxwendroffSolverNewParams};
    use solver::leapfrog_solver::{LeapfrogSolver, LeapfrogSolverNewParams};
    use solver::maccormack_solver::{MaccormackSolver, MaccormackSolverNewParams};
    use solver::upwind_solver::{UpwindSolver, UpwindSolverNewParams};

    #[test]
    fn fn_run_works_with_ftcs_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = FtcsSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut solver = FtcsSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 0.9997558594
6 -0.5000000000 1.0056152344
6 -0.4000000000 0.9484863281
6 -0.3000000000 1.2316894531
6 -0.2000000000 0.5249023438
6 -0.1000000000 1.1459960938
6 0.0000000000 1.6743164062
6 0.1000000000 1.0532226562
6 0.2000000000 0.3464355469
6 0.3000000000 0.0632324219
6 0.4000000000 0.0061035156
6 0.5000000000 0.0002441406
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_lax_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = LaxSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut solver = LaxSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 0.9997558594
6 -0.5000000000 0.9997558594
6 -0.4000000000 0.9953613281
6 -0.3000000000 0.9953613281
6 -0.2000000000 0.9624023438
6 -0.1000000000 0.9624023438
6 0.0000000000 0.8305664062
6 0.1000000000 0.8305664062
6 0.2000000000 0.5339355469
6 0.3000000000 0.5339355469
6 0.4000000000 0.1779785156
6 0.5000000000 0.1779785156
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_leapfrog_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = LeapfrogSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 1.0,
        };
        let mut solver = LeapfrogSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 0.9843750000
6 -0.5000000000 1.0156250000
6 -0.4000000000 0.7968750000
6 -0.3000000000 1.1406250000
6 -0.2000000000 0.6562500000
6 -0.1000000000 0.9687500000
6 0.0000000000 1.4062500000
6 0.1000000000 1.0937500000
6 0.2000000000 0.6093750000
6 0.3000000000 0.2656250000
6 0.4000000000 0.0468750000
6 0.5000000000 0.0156250000
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_laxwendroff_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = LaxwendroffSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut solver = LaxwendroffSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 0.9999961853
6 -0.5000000000 1.0001335144
6 -0.4000000000 0.9981422424
6 -0.3000000000 1.0125617981
6 -0.2000000000 0.9626083374
6 -0.1000000000 1.0046310425
6 0.0000000000 1.1624221802
6 0.1000000000 1.0363540649
6 0.2000000000 0.5867729187
6 0.3000000000 0.1974449158
6 0.4000000000 0.0361518860
6 0.5000000000 0.0027809143
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_maccormack_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = MaccormackSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut solver = MaccormackSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 0.9999961853
6 -0.5000000000 1.0001335144
6 -0.4000000000 0.9981422424
6 -0.3000000000 1.0125617981
6 -0.2000000000 0.9626083374
6 -0.1000000000 1.0046310425
6 0.0000000000 1.1624221802
6 0.1000000000 1.0363540649
6 0.2000000000 0.5867729187
6 0.3000000000 0.1974449158
6 0.4000000000 0.0361518860
6 0.5000000000 0.0027809143
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_upwind_solver() {
        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup coordinates
        let x: Array1<f64> = Array1::linspace(-1.0, 1.0, 20 + 1);

        // initialize the solver
        let new_params = UpwindSolverNewParams {
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 6,
            n_cfl: 0.5,
        };
        let mut solver = UpwindSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 6).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


6 -1.0000000000 1.0000000000
6 -0.9000000000 1.0000000000
6 -0.8000000000 1.0000000000
6 -0.7000000000 1.0000000000
6 -0.6000000000 1.0000000000
6 -0.5000000000 1.0000000000
6 -0.4000000000 1.0000000000
6 -0.3000000000 1.0000000000
6 -0.2000000000 1.0000000000
6 -0.1000000000 1.0000000000
6 0.0000000000 0.9843750000
6 0.1000000000 0.8906250000
6 0.2000000000 0.6562500000
6 0.3000000000 0.3437500000
6 0.4000000000 0.1093750000
6 0.5000000000 0.0156250000
6 0.6000000000 0.0000000000
6 0.7000000000 0.0000000000
6 0.8000000000 0.0000000000
6 0.9000000000 0.0000000000
6 1.0000000000 0.0000000000


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
            u: x.map(|x| if *x < 0.0 { 1.0 } else { 0.0 }),
            step_max: 3,
            n_cfl: 1.0,
            lambda: 0.5,
        };
        let mut solver = BeamwarmingSolver::new(new_params).unwrap();

        // execute run()
        run(&x, &mut solver, &mut outputstream, 3).unwrap();

        // check if the output is correct
        let output_expected = "\
0 -1.0000000000 1.0000000000
0 -0.9000000000 1.0000000000
0 -0.8000000000 1.0000000000
0 -0.7000000000 1.0000000000
0 -0.6000000000 1.0000000000
0 -0.5000000000 1.0000000000
0 -0.4000000000 1.0000000000
0 -0.3000000000 1.0000000000
0 -0.2000000000 1.0000000000
0 -0.1000000000 1.0000000000
0 0.0000000000 0.0000000000
0 0.1000000000 0.0000000000
0 0.2000000000 0.0000000000
0 0.3000000000 0.0000000000
0 0.4000000000 0.0000000000
0 0.5000000000 0.0000000000
0 0.6000000000 0.0000000000
0 0.7000000000 0.0000000000
0 0.8000000000 0.0000000000
0 0.9000000000 0.0000000000
0 1.0000000000 0.0000000000


3 -1.0000000000 1.0000000000
3 -0.9000000000 0.7769564522
3 -0.8000000000 0.8338211969
3 -0.7000000000 0.9186703690
3 -0.6000000000 0.9522185033
3 -0.5000000000 1.0207447161
3 -0.4000000000 0.9086533220
3 -0.3000000000 1.1851764733
3 -0.2000000000 0.7183036043
3 -0.1000000000 1.1268431186
3 0.0000000000 1.3398428513
3 0.1000000000 0.9316474270
3 0.2000000000 0.4637756569
3 0.3000000000 0.1903174095
3 0.4000000000 0.0695039387
3 0.5000000000 0.0235061476
3 0.6000000000 0.0075308953
3 0.7000000000 0.0023180371
3 0.8000000000 0.0006913478
3 0.9000000000 0.0002031287
3 1.0000000000 0.0000000000


";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }
}
