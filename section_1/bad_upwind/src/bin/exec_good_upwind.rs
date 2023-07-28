//! Solve the transport equation by the good upwind method, in this case, [DiffMethod::Backward].
//!
//! See [bad_upwind] for details.

use bad_upwind::input;
use bad_upwind::upwind_solver::DiffMethod;
use std::fs::{self, File};
use std::process;

/// Solve the equation with the given input parameters and output the result to a file.
fn main() {
    // read input parameters
    let mut inputfile = File::open("inputs/section_1/bad_upwind/exec_good_upwind.yml")
        .unwrap_or_else(|err| {
            eprintln!("Problem opening input file: {}", err);
            process::exit(1);
        });
    let input_params = input::read_input_params(&mut inputfile).unwrap_or_else(|err| {
        eprintln!("Problem reading input parameters: {}", err);
        process::exit(1);
    });

    // setup output files
    let dir_str = "outputs/section_1/bad_upwind";
    fs::create_dir_all(dir_str).unwrap_or_else(|err| {
        eprintln!("Problem creating output directory: {}", err);
        process::exit(1);
    });
    let mut outputfile =
        File::create(format!("{}/exec_good_upwind.dat", dir_str)).unwrap_or_else(|err| {
            eprintln!("Problem creating output files: {}", err);
            process::exit(1);
        });

    // run
    bad_upwind::run(&input_params, DiffMethod::Backward, &mut outputfile).unwrap_or_else(|err| {
        eprintln!("Application error: {}", err);
        process::exit(1);
    });
}
