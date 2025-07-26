//! # Quantum-Bessel-RS (`qbessel-rs`)
//! 
//! A high-performance Rust library for Bessel/Neumann functions with quantum cryptography applications.
//! 
//! ## Core Functions
//! 
//! - `Jν(x)` - Bessel functions of the first kind
//! - `Yν(x)` - Neumann functions (Bessel functions of the second kind)  
//! - `Hν(x)` - Hankel functions
//! - `Iν(x), Kν(x)` - Modified Bessel functions
//! 
//! ## Example Usage
//! 
//! ```rust
//! use qbessel_rs::core::{bessel_j, neumann_y, hankel_h1, modified_bessel_i};
//! 
//! // Compute J₀(2.0)
//! let j0_result = bessel_j(0.0, 2.0).unwrap();
//! 
//! // Compute Y₁(1.5)
//! let y1_result = neumann_y(1.0, 1.5).unwrap();
//! 
//! // Compute H₁₀(1.0)
//! let h1_result = hankel_h1(0.0, 1.0).unwrap();
//! 
//! // Compute I₀(0.5)
//! let i0_result = modified_bessel_i(0.0, 0.5).unwrap();
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod core;
pub mod error;

// Re-export core functions for convenience
pub use core::*;
pub use error::{BesselError, BesselResult};
