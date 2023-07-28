//! Module to output the results.

use ndarray::prelude::*;
use std::io::{Error, Write};

/// Output the results.
///
/// # Output Format
/// The output is formatted as follows:
/// ```text
/// t_0 x_0 u_0
/// t_0 x_1 u_1
/// t_0 x_2 u_2
/// ...
/// t_0 x_n u_n
///
///
/// t_1 x_0 u_0
/// t_1 x_1 u_1
/// t_1 x_2 u_2
/// ...
/// t_1 x_n u_n
///
///
/// ...
/// t_m x_0 u_0
/// t_m x_1 u_1
/// t_m x_2 u_2
/// ...
/// t_m x_n u_n
/// ```
///
/// # Examples
/// ```
/// let t = 3.0;
/// let x = ndarray::array![-1.0, 0.0, 1.0];
/// let u = ndarray::array![0.0, 1.0, 2.0];
/// let mut outputstream: Vec<u8> = Vec::new();
/// bad_upwind::output::output(&mut outputstream, t, &x, &u).unwrap();
///
/// let output_expected = "\
/// 3.00 -1.0000000000 0.0000000000
/// 3.00 0.0000000000 1.0000000000
/// 3.00 1.0000000000 2.0000000000
///
///
/// ";
/// assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
/// ```
///
/// # Errors
/// Returns an error if output fails.
pub fn output(
    outputstream: &mut impl Write,
    t: f64,
    x: &Array1<f64>,
    u: &Array1<f64>,
) -> Result<(), Error> {
    for (x, u) in x.iter().zip(u.iter()) {
        writeln!(outputstream, "{:.2} {:.10} {:.10}", t, x, u)?;
    }
    writeln!(outputstream)?;
    writeln!(outputstream)?;

    Ok(())
}
