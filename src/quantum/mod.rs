//! Quantum cryptography modules using Bessel functions
//! 
//! This module provides specialized implementations for quantum cryptographic
//! applications, including Quantum Key Distribution (QKD) and lattice-based
//! post-quantum cryptography.

pub mod qkd;
pub mod lattice;

pub use qkd::{BesselBeam, BeamProfile};
pub use lattice::{BesselNoiseSampler, BesselFunctionType, NeumannNoiseSampler, LatticeProblemGenerator};
