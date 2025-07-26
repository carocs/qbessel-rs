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

/// Error type for quantum Bessel function applications
#[derive(Debug, Clone)]
pub enum QuantumBesselError {
    /// Invalid input parameter
    InvalidInput(String),
    /// Computation error
    ComputationError(String),
    /// Bessel function error
    BesselError(BesselError),
}

impl fmt::Display for QuantumBesselError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QuantumBesselError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
            QuantumBesselError::ComputationError(msg) => write!(f, "Computation error: {}", msg),
            QuantumBesselError::BesselError(err) => write!(f, "Bessel error: {}", err),
        }
    }
}

impl From<BesselError> for QuantumBesselError {
    fn from(err: BesselError) -> Self {
        QuantumBesselError::BesselError(err)
    }
}

impl From<crate::core::func::BesselError> for QuantumBesselError {
    fn from(err: crate::core::func::BesselError) -> Self {
        // Convert from func::BesselError to error::BesselError
        let converted_err = match err {
            crate::core::func::BesselError::InvalidArgument(msg) => BesselError::InvalidArgument(msg),
            crate::core::func::BesselError::ConvergenceError(msg) => BesselError::ConvergenceError(msg),
            crate::core::func::BesselError::NumericalError(msg) => BesselError::NumericalError(msg),
        };
        QuantumBesselError::BesselError(converted_err)
    }
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

#[cfg(feature = "std")]
impl std::error::Error for QuantumBesselError {}

/// Result type for Bessel function computations
pub type BesselResult<T> = Result<T, BesselError>;
