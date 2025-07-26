//! Error types for the qbessel-rs library

use core::fmt;

/// Error type for Bessel function computations
#[derive(Debug, Clone)]
pub enum BesselError {
    /// Invalid argument provided to a function
    InvalidArgument(String),
    /// Convergence error in iterative algorithms
    ConvergenceError(String),
    /// Numerical error due to precision limits
    NumericalError(String),
}

impl fmt::Display for BesselError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BesselError::InvalidArgument(msg) => write!(f, "Invalid argument: {}", msg),
            BesselError::ConvergenceError(msg) => write!(f, "Convergence error: {}", msg),
            BesselError::NumericalError(msg) => write!(f, "Numerical error: {}", msg),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for BesselError {}

/// Result type for Bessel function computations
pub type BesselResult<T> = Result<T, BesselError>;
