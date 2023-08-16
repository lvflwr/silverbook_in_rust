//! Module to output the results.

use ndarray::prelude::*;
use std::io::{Error, Write};

/// Output the results.
///
/// # Output Format
/// The output is formatted as follows:
/// ```text
/// x0 y0 u_x0_y0
/// x0 y1 u_x0_y1
/// x0 y2 u_x0_y2
/// ...
/// x0 ym u_x0_ym
///
/// x1 y0 u_x1_y0
/// x1 y1 u_x1_y1
/// x1 y2 u_x1_y2
/// ...
/// x1 ym u_x1_ym
///
/// ...
/// xn y0 u_xn_y0
/// xn y1 u_xn_y1
/// xn y2 u_xn_y2
/// ...
/// xn ym u_xn_ym
/// ```
///
/// # Examples
/// ```
/// use ndarray::prelude::*;
/// use elliptic::output;
///
/// let mut outputstream: Vec<u8> = Vec::new();
/// let u = array![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]];
/// output::output(&mut outputstream, &u).unwrap();
///
/// let output_expected = "\
/// 0 0 0.0000000000
/// 0 1 1.0000000000
/// 0 2 2.0000000000
///
/// 1 0 3.0000000000
/// 1 1 4.0000000000
/// 1 2 5.0000000000
///
/// 2 0 6.0000000000
/// 2 1 7.0000000000
/// 2 2 8.0000000000
///
/// ";
/// assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
/// ```
///
/// # Errors
/// Returns an error if the output fails.
pub fn output(outputstream: &mut impl Write, u: &Array2<f64>) -> Result<(), Error> {
    for (i_x, u_at_x) in u.outer_iter().enumerate() {
        for (i_y, u_val) in u_at_x.iter().enumerate() {
            writeln!(outputstream, "{} {} {:.10}", i_x, i_y, u_val)?;
        }
        writeln!(outputstream)?;
    }

    Ok(())
}
