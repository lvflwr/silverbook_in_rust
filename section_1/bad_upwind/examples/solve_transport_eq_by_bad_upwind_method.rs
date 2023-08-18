//! Solve the transport equation by the bad upwind method, in this case, [DiffMethod::Forward].
//!
//! # Formulation
//! The transport equation is given by
//! ```math
//! \frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0 (x \in [-1, 1])),
//! ```
//! where `u` is the transported quantity and `c` (`> 0`) is the advection velocity.
//!
//! The initial condition is given by
//! ```math
//! u(x, 0) = 0 (x \ge 0), u(x, 0) = 1 (x < 0).
//! ```
//!
//! For the boundary condition, see [bad_upwind::upwind_solver].
//!
//! # Scheme
//! See [DiffMethod::Forward].
//!
//! # Input Format
//! See [input::read_input_params].
//!
//! # Output Format
//! See [bad_upwind::output::output].

use bad_upwind::input;
use bad_upwind::upwind_solver::{DiffMethod, UpwindSolver};
use ndarray::prelude::*;
use std::fs::{self, File};
use std::process;

/// Solve the equation with the given input parameters and output the result to a file.
fn main() {
    // read input parameters
    let mut inputfile =
        File::open("inputs/section_1/bad_upwind/solve_transport_eq_by_bad_upwind_method/input.yml")
            .unwrap_or_else(|err| {
                eprintln!("Problem opening input file: {}", err);
                process::exit(1);
            });
    let input_params = input::read_input_params(&mut inputfile).unwrap_or_else(|err| {
        eprintln!("Problem reading input parameters: {}", err);
        process::exit(1);
    });

    // setup output files
    let dir_str = "outputs/section_1/bad_upwind/solve_transport_eq_by_bad_upwind_method";
    fs::create_dir_all(dir_str).unwrap_or_else(|err| {
        eprintln!("Problem creating output directory: {}", err);
        process::exit(1);
    });
    let mut outputfile = File::create(format!("{}/solution.dat", dir_str)).unwrap_or_else(|err| {
        eprintln!("Problem creating output files: {}", err);
        process::exit(1);
    });

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

    // run
    bad_upwind::run(
        &x,
        &mut upwind_solver,
        &mut outputfile,
        input_params.ncycle_out,
    )
    .unwrap_or_else(|err| {
        eprintln!("Application error: {}", err);
        process::exit(1);
    });
}
