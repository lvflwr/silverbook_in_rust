//! Module to read the input parameters.

use serde::de::DeserializeOwned;
use serde::Serialize;
use std::error::Error;
use std::io::prelude::*;

/// Read the input parameters from the input.
///
/// The format of the input should be defined by a struct that implements [InputParams], [Serialize] and [DeserializeOwned].
///
/// # Examples
/// ```
/// use serde_derive::{Deserialize, Serialize};
/// use linear_hyperbolic::input::{self, InputParams};
///
/// #[derive(Debug, Serialize, Deserialize, PartialEq)]
/// pub struct SpecificInputParams {
///    pub a: usize,
///    pub b: f64,
///    pub c: f64,
/// }
///
/// impl InputParams for SpecificInputParams {
///     fn validate_params(&self) -> Result<(), &'static str> {
///         if self.b <= 0.0 {
///             return Err("b must be positive");
///         }
///
///         Ok(())
///     }
/// }
///
/// let input_params = SpecificInputParams {
///   a: 3,
///   b: 100.0,
///   c: 1.0,
/// };
/// let input_str = serde_yaml::to_string(&input_params).unwrap();
/// let input_params_read: SpecificInputParams = input::read_input_params(&mut input_str.as_bytes()).unwrap();
///
/// assert_eq!(input_params_read, input_params);
/// ```
///
/// # Errors
/// Returns an error if the input is invalid.
pub fn read_input_params<T: InputParams + Serialize + DeserializeOwned>(
    inputstream: &mut impl Read,
) -> Result<T, Box<dyn Error>> {
    let mut contents = String::new();
    inputstream.read_to_string(&mut contents)?;
    let input_params: T = serde_yaml::from_str(&contents)?;
    input_params.validate_params()?;

    Ok(input_params)
}

/// Input parameters.
pub trait InputParams {
    /// Validate the input parameters.
    fn validate_params(&self) -> Result<(), &'static str>;
}
