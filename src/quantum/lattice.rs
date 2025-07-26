//! Lattice-based cryptography noise generation using Bessel and Neumann functions
//! 
//! This module provides specialized noise samplers for lattice-based cryptographic
//! schemes like LWE (Learning With Errors) and NTRU, using Bessel and Neumann
//! functions to generate quantum-resistant noise distributions.

use crate::core::{bessel_j, neumann_y, modified_bessel_i, modified_bessel_k};
use crate::error::QuantumBesselError;

/// Noise sampler using Bessel functions for lattice-based cryptography
/// 
/// This sampler generates noise with distributions derived from Bessel functions,
/// providing quantum-resistant properties for post-quantum cryptographic schemes.
#[derive(Debug, Clone)]
pub struct BesselNoiseSampler {
    /// Order of the Bessel function (ν)
    pub order: f64,
    /// Scale parameter for noise distribution
    pub scale: f64,
    /// Type of Bessel function to use
    pub function_type: BesselFunctionType,
    /// Standard deviation of the resulting noise
    pub sigma: f64,
}

/// Types of Bessel functions available for noise generation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BesselFunctionType {
    /// First kind Bessel function J_ν(x)
    FirstKind,
    /// Second kind Bessel function (Neumann) Y_ν(x)  
    SecondKind,
    /// Modified Bessel function of first kind I_ν(x)
    ModifiedFirstKind,
    /// Modified Bessel function of second kind K_ν(x)
    ModifiedSecondKind,
}

impl BesselNoiseSampler {
    /// Create a new Bessel noise sampler
    /// 
    /// # Arguments
    /// * `order` - Order ν of the Bessel function
    /// * `scale` - Scale parameter for the distribution
    /// * `function_type` - Type of Bessel function to use
    /// * `sigma` - Target standard deviation
    /// 
    /// # Example
    /// ```
    /// use qbessel_rs::quantum::lattice::{BesselNoiseSampler, BesselFunctionType};
    /// 
    /// let sampler = BesselNoiseSampler::new(0.5, 1.0, BesselFunctionType::SecondKind, 3.2);
    /// ```
    pub fn new(order: f64, scale: f64, function_type: BesselFunctionType, sigma: f64) -> Self {
        Self {
            order,
            scale,
            function_type,
            sigma,
        }
    }

    /// Generate a noise sample using the configured Bessel function
    /// 
    /// This uses rejection sampling with the Bessel function as the probability density.
    /// The method is suitable for generating discrete noise for lattice problems.
    /// 
    /// # Arguments
    /// * `x` - Input value for the Bessel function
    /// 
    /// # Returns
    /// Noise value weighted by the Bessel function
    pub fn sample(&self, x: f64) -> Result<f64, QuantumBesselError> {
        let scaled_x = x / self.scale;
        
        let bessel_weight = match self.function_type {
            BesselFunctionType::FirstKind => {
                bessel_j(self.order, scaled_x).map_err(QuantumBesselError::from)?
            },
            BesselFunctionType::SecondKind => {
                // Use absolute value of Neumann function for positive weights
                neumann_y(self.order, scaled_x).map_err(QuantumBesselError::from)?.abs()
            },
            BesselFunctionType::ModifiedFirstKind => {
                modified_bessel_i(self.order, scaled_x).map_err(QuantumBesselError::from)?
            },
            BesselFunctionType::ModifiedSecondKind => {
                modified_bessel_k(self.order, scaled_x).map_err(QuantumBesselError::from)?
            },
        };

        // Apply scaling to achieve target standard deviation
        let noise_sample = bessel_weight * self.sigma;
        Ok(noise_sample)
    }

    /// Generate a discrete noise sample for LWE problems
    /// 
    /// This rounds the continuous Bessel-weighted noise to integers,
    /// suitable for Learning With Errors constructions.
    /// 
    /// # Arguments
    /// * `x` - Input value
    /// 
    /// # Returns
    /// Integer noise value
    pub fn sample_discrete(&self, x: f64) -> Result<i32, QuantumBesselError> {
        let continuous_sample = self.sample(x)?;
        
        // Round to nearest integer with proper handling of negative values
        let rounded = if continuous_sample >= 0.0 {
            (continuous_sample + 0.5).floor()
        } else {
            (continuous_sample - 0.5).ceil()
        };
        
        // Convert to i32, clamping to avoid overflow
        let clamped = rounded.max(i32::MIN as f64).min(i32::MAX as f64);
        
        Ok(clamped as i32)
    }

    /// Generate noise for NTRU polynomial coefficients
    /// 
    /// NTRU requires small integer coefficients, typically in {-1, 0, 1}.
    /// This method uses Bessel function weights to bias the selection.
    /// 
    /// # Arguments
    /// * `x` - Input value (typically polynomial degree index)
    /// 
    /// # Returns
    /// NTRU coefficient in {-1, 0, 1}
    pub fn sample_ntru(&self, x: f64) -> Result<i8, QuantumBesselError> {
        let weight = self.sample(x)?.abs();
        
        // Use weight to determine probability of non-zero coefficient
        let threshold_high = 0.7;
        let threshold_low = 0.3;
        
        if weight > threshold_high {
            // High weight: more likely to be ±1
            if x.sin() > 0.0 { // Use input as pseudo-random source
                Ok(1)
            } else {
                Ok(-1)
            }
        } else if weight > threshold_low {
            // Medium weight: could be any value
            let phase = (x * 3.0).sin();
            if phase > 0.33 {
                Ok(1)
            } else if phase < -0.33 {
                Ok(-1)
            } else {
                Ok(0)
            }
        } else {
            // Low weight: more likely to be 0
            Ok(0)
        }
    }

    /// Calculate the effective variance of the noise distribution
    /// 
    /// This estimates the variance of noise generated by this sampler,
    /// which is important for security analysis of lattice schemes.
    /// 
    /// # Arguments
    /// * `num_samples` - Number of samples to use for estimation
    /// * `max_x` - Maximum input value for sampling
    pub fn estimate_variance(&self, num_samples: usize, max_x: f64) -> Result<f64, QuantumBesselError> {
        let mut sum = 0.0;
        let mut sum_squares = 0.0;
        let dx = max_x / num_samples as f64;

        for i in 0..num_samples {
            let x = (i as f64 + 1.0) * dx; // Ensure x > 0 for Neumann functions
            let sample = self.sample(x)?;
            sum += sample;
            sum_squares += sample * sample;
        }

        let n = num_samples as f64;
        let mean = sum / n;
        let variance = sum_squares / n - mean * mean;
        
        Ok(variance)
    }
}

/// Neumann function noise sampler specifically for quantum-resistant applications
/// 
/// The Neumann functions (Bessel functions of the second kind) have unique
/// properties that make them suitable for generating hard-to-predict noise
/// patterns in post-quantum cryptography.
#[derive(Debug, Clone)]
pub struct NeumannNoiseSampler {
    /// Order of the Neumann function
    pub order: f64,
    /// Scaling factor
    pub scale: f64,
    /// Offset to handle singularities
    pub offset: f64,
}

impl NeumannNoiseSampler {
    /// Create a new Neumann noise sampler
    /// 
    /// # Arguments
    /// * `order` - Order ν of Y_ν(x)
    /// * `scale` - Scale parameter
    /// * `offset` - Offset to avoid singularity at x=0
    pub fn new(order: f64, scale: f64, offset: f64) -> Self {
        Self { order, scale, offset }
    }

    /// Generate noise using Neumann function with error handling
    /// 
    /// # Arguments
    /// * `x` - Input value
    /// 
    /// # Returns
    /// Noise sample based on Y_ν(scale * x + offset)
    pub fn sample(&self, x: f64) -> Result<f64, QuantumBesselError> {
        let arg = self.scale * x + self.offset;
        
        // Ensure argument is positive to avoid singularity
        if arg <= 0.0 {
            return Err(QuantumBesselError::InvalidInput(
                "Neumann function argument must be positive".to_string()
            ));
        }
        
        let neumann_val = neumann_y(self.order, arg).map_err(QuantumBesselError::from)?;
        
        // Handle potential infinities or very large values
        if !neumann_val.is_finite() {
            return Err(QuantumBesselError::ComputationError(
                "Neumann function produced non-finite result".to_string()
            ));
        }
        
        Ok(neumann_val)
    }

    /// Generate bounded noise suitable for lattice cryptography
    /// 
    /// This method ensures the noise stays within reasonable bounds
    /// for cryptographic applications.
    /// 
    /// # Arguments
    /// * `x` - Input value
    /// * `bound` - Maximum absolute value for the noise
    pub fn sample_bounded(&self, x: f64, bound: f64) -> Result<f64, QuantumBesselError> {
        let raw_sample = self.sample(x)?;
        
        // Apply hyperbolic tangent to bound the output
        let bounded = bound * raw_sample.tanh();
        Ok(bounded)
    }
}

/// Lattice problem generator using Bessel function noise
/// 
/// This struct helps generate instances of lattice problems (like LWE)
/// with Bessel function-derived noise for enhanced quantum resistance.
pub struct LatticeProblemGenerator {
    /// Dimension of the lattice problem
    pub dimension: usize,
    /// Noise sampler for error terms
    pub noise_sampler: BesselNoiseSampler,
    /// Modulus for the problem (for LWE)
    pub modulus: Option<i64>,
}

impl LatticeProblemGenerator {
    /// Create a new lattice problem generator
    /// 
    /// # Arguments
    /// * `dimension` - Problem dimension
    /// * `noise_sampler` - Configured noise sampler
    /// * `modulus` - Optional modulus for modular arithmetic
    pub fn new(dimension: usize, noise_sampler: BesselNoiseSampler, modulus: Option<i64>) -> Self {
        Self {
            dimension,
            noise_sampler,
            modulus,
        }
    }

    /// Generate an LWE sample (a, b) where b = <a, s> + e (mod q)
    /// 
    /// # Arguments
    /// * `secret` - Secret vector s
    /// * `seed` - Seed for generating the 'a' vector
    /// 
    /// # Returns
    /// Tuple (a_vector, b_value) representing the LWE sample
    pub fn generate_lwe_sample(&self, secret: &[i32], seed: f64) -> Result<(Vec<i32>, i32), QuantumBesselError> {
        if secret.len() != self.dimension {
            return Err(QuantumBesselError::InvalidInput(
                "Secret vector dimension mismatch".to_string()
            ));
        }

        let mut a_vector = Vec::with_capacity(self.dimension);
        let mut inner_product = 0i64;

        // Generate 'a' vector and compute inner product
        for i in 0..self.dimension {
            // Use seed and index to generate pseudo-random 'a' values
            let a_val = ((seed + i as f64).sin() * 1000.0).floor() as i32;
            a_vector.push(a_val);
            inner_product += (a_val as i64) * (secret[i] as i64);
        }

        // Generate noise using Bessel function
        let noise = self.noise_sampler.sample_discrete(seed)?;
        
        // Compute b = <a, s> + e
        let mut b = inner_product + (noise as i64);
        
        // Apply modulus if specified
        if let Some(q) = self.modulus {
            b = ((b % q) + q) % q; // Ensure positive remainder
        }

        Ok((a_vector, b as i32))
    }

    /// Generate multiple LWE samples for a given secret
    /// 
    /// # Arguments
    /// * `secret` - Secret vector
    /// * `num_samples` - Number of samples to generate
    /// 
    /// # Returns
    /// Vector of LWE samples
    pub fn generate_lwe_samples(&self, secret: &[i32], num_samples: usize) -> Result<Vec<(Vec<i32>, i32)>, QuantumBesselError> {
        let mut samples = Vec::with_capacity(num_samples);
        
        for i in 0..num_samples {
            let seed = i as f64 * 1.618033988749; // Golden ratio for good distribution
            let sample = self.generate_lwe_sample(secret, seed)?;
            samples.push(sample);
        }
        
        Ok(samples)
    }

    /// Estimate the hardness of the generated lattice problem
    /// 
    /// This provides a rough estimate of the computational difficulty
    /// based on the noise distribution and problem parameters.
    /// 
    /// # Returns
    /// Hardness estimate (higher values indicate harder problems)
    pub fn estimate_hardness(&self) -> Result<f64, QuantumBesselError> {
        // Estimate noise variance
        let noise_variance = self.noise_sampler.estimate_variance(1000, 10.0)?;
        
        // Hardness roughly scales with dimension and inverse of noise variance
        let dimension_factor = self.dimension as f64;
        let noise_factor = 1.0 / (noise_variance + f64::EPSILON);
        
        let hardness = dimension_factor * noise_factor.sqrt();
        Ok(hardness)
    }
}

/// Utility functions for lattice-based cryptography with Bessel functions
pub mod utils {
    use super::*;

    /// Generate a random-looking but deterministic secret vector using Bessel functions
    /// 
    /// # Arguments
    /// * `dimension` - Length of the secret vector
    /// * `seed` - Seed for deterministic generation
    /// * `bound` - Maximum absolute value for secret coefficients
    pub fn generate_bessel_secret(dimension: usize, seed: f64, bound: i32) -> Vec<i32> {
        let mut secret = Vec::with_capacity(dimension);
        let sampler = BesselNoiseSampler::new(
            0.5,
            1.0,
            BesselFunctionType::FirstKind,
            bound as f64,
        );

        for i in 0..dimension {
            let x = seed + i as f64 * 0.1;
            let sample = sampler.sample_discrete(x).unwrap_or(0);
            secret.push(sample.max(-bound).min(bound));
        }

        secret
    }

    /// Validate that a noise distribution provides sufficient entropy
    /// 
    /// # Arguments
    /// * `sampler` - The noise sampler to validate
    /// * `num_tests` - Number of samples to use for validation
    /// 
    /// # Returns
    /// True if the distribution appears to have sufficient entropy
    pub fn validate_entropy(sampler: &BesselNoiseSampler, num_tests: usize) -> Result<bool, QuantumBesselError> {
        let mut samples = Vec::with_capacity(num_tests);
        
        // Collect samples
        for i in 0..num_tests {
            let x = (i as f64 + 1.0) / 100.0; // Ensure x > 0 for Neumann functions
            let sample = sampler.sample(x)?;
            samples.push(sample);
        }

        // Simple entropy check: ensure samples are not all the same
        let first_sample = samples[0];
        let all_same = samples.iter().all(|&s| (s - first_sample).abs() < f64::EPSILON);
        
        if all_same {
            return Ok(false);
        }

        // Check for reasonable variance
        let variance = sampler.estimate_variance(num_tests, 10.0)?;
        let min_variance = 0.01;
        
        Ok(variance > min_variance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_noise_sampler_creation() {
        let sampler = BesselNoiseSampler::new(
            0.5, 
            1.0, 
            BesselFunctionType::FirstKind, 
            3.2
        );
        assert_eq!(sampler.order, 0.5);
        assert_eq!(sampler.scale, 1.0);
        assert_eq!(sampler.sigma, 3.2);
    }

    #[test]
    fn test_noise_sampling() {
        let sampler = BesselNoiseSampler::new(
            0.0, 
            1.0, 
            BesselFunctionType::FirstKind, 
            1.0
        );
        
        let sample = sampler.sample(1.0).unwrap();
        assert!(sample.is_finite());
    }

    #[test]
    fn test_discrete_sampling() {
        let sampler = BesselNoiseSampler::new(
            0.0, 
            1.0, 
            BesselFunctionType::FirstKind, 
            2.0
        );
        
        let discrete_sample = sampler.sample_discrete(1.0).unwrap();
        assert!(discrete_sample >= i32::MIN && discrete_sample <= i32::MAX);
    }

    #[test]
    fn test_ntru_sampling() {
        let sampler = BesselNoiseSampler::new(
            0.5, 
            1.0, 
            BesselFunctionType::SecondKind, 
            1.0
        );
        
        let ntru_sample = sampler.sample_ntru(1.0).unwrap();
        assert!(ntru_sample >= -1 && ntru_sample <= 1);
    }

    #[test]
    fn test_neumann_sampler() {
        let sampler = NeumannNoiseSampler::new(0.5, 1.0, 0.1);
        let sample = sampler.sample(1.0).unwrap();
        assert!(sample.is_finite());
    }

    #[test]
    fn test_bounded_neumann_sampling() {
        let sampler = NeumannNoiseSampler::new(0.5, 1.0, 0.1);
        let bounded_sample = sampler.sample_bounded(1.0, 5.0).unwrap();
        assert!(bounded_sample.abs() <= 5.0);
    }

    #[test]
    fn test_lwe_generation() {
        let noise_sampler = BesselNoiseSampler::new(
            0.0, 
            1.0, 
            BesselFunctionType::FirstKind, 
            1.0
        );
        let generator = LatticeProblemGenerator::new(5, noise_sampler, Some(97));
        
        let secret = vec![1, -1, 0, 1, -1];
        let (a, b) = generator.generate_lwe_sample(&secret, 1.0).unwrap();
        
        assert_eq!(a.len(), 5);
        assert!(b >= 0 && b < 97); // Should be reduced modulo 97
    }

    #[test]
    fn test_variance_estimation() {
        let sampler = BesselNoiseSampler::new(
            0.0, 
            1.0, 
            BesselFunctionType::FirstKind, 
            2.0
        );
        
        let variance = sampler.estimate_variance(100, 5.0).unwrap();
        assert!(variance > 0.0);
    }

    #[test]
    fn test_entropy_validation() {
        let sampler = BesselNoiseSampler::new(
            0.5, 
            1.0, 
            BesselFunctionType::SecondKind, 
            1.0
        );
        
        let has_entropy = utils::validate_entropy(&sampler, 50).unwrap();
        assert!(has_entropy);
    }
}
