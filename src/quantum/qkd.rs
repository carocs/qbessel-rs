//! Quantum Key Distribution (QKD) utilities using Bessel functions
//! 
//! This module provides tools for modeling Bessel beams in quantum cryptography
//! applications, particularly for QKD systems where beam propagation and 
//! intensity profiles are critical for security analysis.

use crate::core::{bessel_j, BesselResult};
use crate::error::QuantumBesselError;

#[cfg(feature = "std")]
use std::f64::consts::PI;
#[cfg(not(feature = "std"))]
use libm::PI;

/// Represents a Bessel beam for QKD applications
/// 
/// Bessel beams are non-diffracting solutions to the wave equation,
/// making them ideal for long-distance quantum communication.
#[derive(Debug, Clone)]
pub struct BesselBeam {
    /// Transverse wave vector component (k_ρ)
    pub k_rho: f64,
    /// Beam order (typically 0 for fundamental mode)
    pub order: f64,
    /// Wavelength of the beam
    pub wavelength: f64,
    /// Maximum propagation distance before significant diffraction
    pub max_distance: Option<f64>,
}

impl BesselBeam {
    /// Create a new Bessel beam with specified transverse wave vector
    /// 
    /// # Arguments
    /// * `k_rho` - Transverse wave vector component
    /// * `order` - Beam order (0 for fundamental Bessel beam)
    /// * `wavelength` - Wavelength of the beam
    /// 
    /// # Example
    /// ```
    /// use qbessel_rs::quantum::qkd::BesselBeam;
    /// 
    /// let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
    /// ```
    pub fn new(k_rho: f64, order: f64, wavelength: f64) -> Self {
        let max_distance = if k_rho > 0.0 {
            // Approximate non-diffracting distance: z_max ≈ R/tan(θ) where θ = k_rho/k
            let k_total = 2.0 * PI / wavelength;
            Some(1e-3 / (k_rho / k_total)) // Assuming 1mm aperture
        } else {
            None
        };

        Self {
            k_rho,
            order,
            wavelength,
            max_distance,
        }
    }

    /// Calculate the intensity profile at radius r
    /// 
    /// For a Bessel beam, the intensity is proportional to |J_n(k_ρ * r)|²
    /// where n is the beam order.
    /// 
    /// # Arguments
    /// * `radius` - Radial distance from beam axis
    /// 
    /// # Returns
    /// Normalized intensity at the given radius
    pub fn intensity_profile(&self, radius: f64) -> BesselResult<f64> {
        let arg = self.k_rho * radius;
        let bessel_val = bessel_j(self.order, arg)?;
        Ok(bessel_val * bessel_val) // |J_n(k_ρ * r)|²
    }

    /// Calculate the beam's central intensity (at r = 0)
    /// 
    /// For J_0 beams, this is always 1.0
    /// For higher order beams, this is 0.0
    pub fn central_intensity(&self) -> f64 {
        if self.order == 0.0 {
            1.0
        } else {
            0.0
        }
    }

    /// Check for significant diffraction effects
    /// 
    /// Returns true if the propagation distance exceeds the non-diffracting range,
    /// which could compromise QKD security by introducing beam spreading.
    /// 
    /// # Arguments
    /// * `distance` - Propagation distance to check
    pub fn diffraction_check(&self, distance: f64) -> bool {
        match self.max_distance {
            Some(max_dist) => distance > max_dist,
            None => true, // Conservative: assume diffraction if max_distance unknown
        }
    }

    /// Calculate the beam's effective radius containing a given fraction of power
    /// 
    /// This is useful for determining the required detector size in QKD systems.
    /// 
    /// # Arguments
    /// * `power_fraction` - Fraction of total power to contain (0.0 to 1.0)
    /// * `max_radius` - Maximum radius to search within
    /// * `steps` - Number of integration steps for numerical calculation
    pub fn effective_radius(&self, power_fraction: f64, max_radius: f64, steps: usize) -> Result<f64, QuantumBesselError> {
        if power_fraction <= 0.0 || power_fraction > 1.0 {
            return Err(QuantumBesselError::InvalidInput("Power fraction must be between 0 and 1".to_string()));
        }

        let dr = max_radius / steps as f64;
        let mut total_power = 0.0;
        let mut enclosed_power = 0.0;
        let target_power = power_fraction;

        // Numerical integration using trapezoidal rule
        // P(r) = 2π ∫₀ʳ I(r') r' dr' for cylindrical symmetry
        for i in 0..=steps {
            let r = i as f64 * dr;
            let intensity = self.intensity_profile(r).map_err(QuantumBesselError::from)?;
            let power_element = if i == 0 || i == steps {
                intensity * r * dr / 2.0 // Trapezoidal rule endpoints
            } else {
                intensity * r * dr
            };
            
            total_power += power_element;
            
            if i < steps {
                enclosed_power += power_element;
                let current_fraction = enclosed_power / (total_power + f64::EPSILON);
                
                if current_fraction >= target_power {
                    return Ok(r);
                }
            }
        }

        // If we reach here, the power fraction requires radius > max_radius
        Ok(max_radius)
    }

    /// Calculate the quantum bit error rate (QBER) contribution from beam imperfections
    /// 
    /// This estimates the error rate introduced by non-ideal Bessel beam profiles
    /// in QKD systems, considering detector positioning and beam alignment.
    /// 
    /// # Arguments
    /// * `detector_radius` - Radius of the detector aperture
    /// * `misalignment` - Lateral misalignment of detector from beam center
    pub fn qber_contribution(&self, detector_radius: f64, misalignment: f64) -> Result<f64, QuantumBesselError> {
        // Calculate power collected by detector considering misalignment
        let steps = 100;
        let dr = detector_radius / steps as f64;
        let mut collected_power = 0.0;
        let mut total_power = 0.0;

        // Integrate over detector area, accounting for misalignment
        for i in 0..=steps {
            let r = i as f64 * dr;
            // Consider detector positioned at distance 'misalignment' from beam center
            let effective_r = (r * r + misalignment * misalignment).sqrt();
            let intensity = self.intensity_profile(effective_r).map_err(QuantumBesselError::from)?;
            
            let power_element = if i == 0 || i == steps {
                intensity * r * dr / 2.0
            } else {
                intensity * r * dr
            };
            
            collected_power += power_element;
        }

        // Calculate total beam power (reference)
        let reference_radius = detector_radius * 10.0; // Large reference
        let dr_ref = reference_radius / steps as f64;
        
        for i in 0..=steps {
            let r = i as f64 * dr_ref;
            let intensity = self.intensity_profile(r).map_err(QuantumBesselError::from)?;
            let power_element = if i == 0 || i == steps {
                intensity * r * dr_ref / 2.0
            } else {
                intensity * r * dr_ref
            };
            total_power += power_element;
        }

        // QBER contribution is roughly proportional to power loss
        let collection_efficiency = collected_power / (total_power + f64::EPSILON);
        let qber_contrib = (1.0 - collection_efficiency) / 4.0; // Simplified model
        
        Ok(qber_contrib.min(0.5)) // Cap at 50% error rate
    }
}

/// Beam profile analysis for QKD security assessment
pub struct BeamProfile {
    /// Radial positions for profile sampling
    pub radii: Vec<f64>,
    /// Intensity values at corresponding radii
    pub intensities: Vec<f64>,
    /// Peak intensity value
    pub peak_intensity: f64,
    /// Full width at half maximum (FWHM)
    pub fwhm: Option<f64>,
}

impl BeamProfile {
    /// Generate a beam profile for analysis
    /// 
    /// # Arguments
    /// * `beam` - The Bessel beam to analyze
    /// * `max_radius` - Maximum radius for profile generation
    /// * `num_points` - Number of sampling points
    pub fn generate(beam: &BesselBeam, max_radius: f64, num_points: usize) -> Result<Self, QuantumBesselError> {
        let mut radii = Vec::with_capacity(num_points);
        let mut intensities = Vec::with_capacity(num_points);
        let mut peak_intensity = 0.0;

        let dr = max_radius / (num_points - 1) as f64;

        for i in 0..num_points {
            let r = i as f64 * dr;
            let intensity = beam.intensity_profile(r).map_err(QuantumBesselError::from)?;
            
            radii.push(r);
            intensities.push(intensity);
            
            if intensity > peak_intensity {
                peak_intensity = intensity;
            }
        }

        // Calculate FWHM (Full Width at Half Maximum)
        let half_max = peak_intensity / 2.0;
        let mut fwhm = None;
        
        // Find first crossing of half-maximum
        for i in 1..intensities.len() {
            if intensities[i-1] >= half_max && intensities[i] < half_max {
                // Linear interpolation for more accurate FWHM
                let r1 = radii[i-1];
                let r2 = radii[i];
                let i1 = intensities[i-1];
                let i2 = intensities[i];
                
                let r_half = r1 + (r2 - r1) * (half_max - i1) / (i2 - i1);
                fwhm = Some(r_half * 2.0); // Full width
                break;
            }
        }

        Ok(BeamProfile {
            radii,
            intensities,
            peak_intensity,
            fwhm,
        })
    }

    /// Assess beam quality for QKD applications
    /// 
    /// Returns a quality score between 0.0 and 1.0, where 1.0 indicates
    /// an ideal Bessel beam suitable for QKD.
    pub fn quality_assessment(&self) -> f64 {
        let mut score = 1.0;

        // Penalize if FWHM is too large (beam too spread out)
        if let Some(fwhm) = self.fwhm {
            let ideal_fwhm = 1.0; // Normalized ideal FWHM
            let fwhm_penalty = (fwhm / ideal_fwhm - 1.0).abs() / 10.0;
            score -= fwhm_penalty.min(0.3);
        }

        // Penalize low peak intensity
        let intensity_penalty = (1.0 - self.peak_intensity) / 5.0;
        score -= intensity_penalty.min(0.2);

        score.max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bessel_beam_creation() {
        let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
        assert_eq!(beam.k_rho, 5.0);
        assert_eq!(beam.order, 0.0);
        assert!(beam.max_distance.is_some());
    }

    #[test]
    fn test_central_intensity() {
        let beam_j0 = BesselBeam::new(5.0, 0.0, 632.8e-9);
        let beam_j1 = BesselBeam::new(5.0, 1.0, 632.8e-9);
        
        assert_eq!(beam_j0.central_intensity(), 1.0);
        assert_eq!(beam_j1.central_intensity(), 0.0);
    }

    #[test]
    fn test_intensity_profile() {
        let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
        
        // At r=0, J_0(0) = 1, so intensity should be 1
        let intensity_center = beam.intensity_profile(0.0).unwrap();
        assert_relative_eq!(intensity_center, 1.0, epsilon = 1e-10);
        
        // At some radius, intensity should be positive but less than 1
        let intensity_r = beam.intensity_profile(0.1).unwrap();
        assert!(intensity_r > 0.0);
        assert!(intensity_r < 1.0);
    }

    #[test]
    fn test_diffraction_check() {
        let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
        
        // Short distance should not show diffraction
        assert!(!beam.diffraction_check(1e-6));
        
        // Very long distance should show diffraction
        assert!(beam.diffraction_check(1.0));
    }

    #[test]
    fn test_beam_profile_generation() {
        let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
        let profile = BeamProfile::generate(&beam, 2.0, 100).unwrap();
        
        assert_eq!(profile.radii.len(), 100);
        assert_eq!(profile.intensities.len(), 100);
        assert!(profile.peak_intensity > 0.0);
        assert!(profile.fwhm.is_some());
    }

    #[test]
    fn test_quality_assessment() {
        let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
        let profile = BeamProfile::generate(&beam, 2.0, 100).unwrap();
        let quality = profile.quality_assessment();
        
        assert!(quality >= 0.0);
        assert!(quality <= 1.0);
    }
}
