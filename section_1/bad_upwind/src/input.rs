//! Module to read the input parameters.

use serde_derive::{Deserialize, Serialize};
use std::error::Error;
use std::io::prelude::*;

/// Input parameters.
#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct InputParams {
    /// Advection velocity.
    pub v_adv: f64,
    /// Number of cells.
    pub n_x: usize,
    /// Maximum time.
    pub t_max: f64,
    /// Time step.
    pub dt: f64,
    /// Number of cycles between outputs.
    pub ncycle_out: usize,
}

impl InputParams {
    fn validate_params(&self) -> Result<(), &'static str> {
        if self.v_adv <= 0.0 {
            return Err("v_adv must be positive");
        }
        if self.n_x == 0 {
            return Err("n_x must be positive");
        }
        if self.t_max < self.dt {
            return Err("t_max must be greater than or equal to dt");
        }
        if self.dt <= 0.0 {
            return Err("dt must be positive");
        }
        if self.ncycle_out == 0 {
            return Err("ncycle_out must be positive");
        }

        Ok(())
    }
}

/// Read the input parameters from the input in YAML format.
///
/// # Input Format
/// The input must be formatted as follows:
/// ```yaml
/// v_adv: 1.0
/// n_x: 100
/// t_max: 1.0
/// dt: 0.01
/// ncycle_out: 1
/// ```
///
/// For the meaning of each parameter, see [InputParams].
///
/// # Examples
/// ```
/// use bad_upwind::input::{self, InputParams};
///
/// let input_params = InputParams {
///   v_adv: 1.0,
///   n_x: 100,
///   t_max: 1.0,
///   dt: 0.01,
///   ncycle_out: 1,
/// };
/// let input_str = serde_yaml::to_string(&input_params).unwrap();
/// let input_params_read = input::read_input_params(&mut input_str.as_bytes()).unwrap();
///
/// assert_eq!(input_params_read, input_params);
/// ```
///
/// # Errors
/// Returns an error if the input is invalid.
pub fn read_input_params(inputstream: &mut impl Read) -> Result<InputParams, Box<dyn Error>> {
    let mut contents = String::new();
    inputstream.read_to_string(&mut contents)?;
    let input_params: InputParams = serde_yaml::from_str(&contents)?;
    input_params.validate_params()?;

    Ok(input_params)
}
