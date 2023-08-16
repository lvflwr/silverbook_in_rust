//! Module to output the results.

use ndarray::prelude::*;
use std::io::{Error, Write};

/// Output the results.
///
/// # Output Format
/// The output is formatted as follows:
/// ```text
/// step_0 x_0 u_0
/// step_0 x_1 u_1
/// step_0 x_2 u_2
/// ...
/// step_0 x_n u_n
///
///
/// step_1 x_0 u_0
/// step_1 x_1 u_1
/// step_1 x_2 u_2
/// ...
/// step_1 x_n u_n
///
///
/// ...
/// step_m x_0 u_0
/// step_m x_1 u_1
/// step_m x_2 u_2
/// ...
/// step_m x_n u_n
/// ```
///
/// # Examples
/// ```
/// use ndarray::prelude::*;
/// use linear_hyperbolic::output;
///
/// let mut outputstream: Vec<u8> = Vec::new();
/// let step = 3;
/// let x = array![-1.0, 0.0, 1.0];
/// let u = array![0.0, 1.0, 2.0];
/// output::output(&mut outputstream, step, &x, &u).unwrap();
///
/// let output_expected = "\
/// 3 -1.0000000000 0.0000000000
/// 3 0.0000000000 1.0000000000
/// 3 1.0000000000 2.0000000000
///
///
/// ";
/// assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
/// ```
///
/// # Errors
/// Returns an error if the output fails.
pub fn output(
    outputstream: &mut impl Write,
    step: usize,
    x: &Array1<f64>,
    u: &Array1<f64>,
) -> Result<(), Error> {
    for (x, u) in x.iter().zip(u.iter()) {
        writeln!(outputstream, "{} {:.10} {:.10}", step, x, u)?;
    }
    writeln!(outputstream)?;
    writeln!(outputstream)?;

    Ok(())
}
