//! Example demonstrating quantum cryptography applications using Bessel functions
//! 
//! This example shows how to use the qbessel-rs library for:
//! 1. QKD (Quantum Key Distribution) beam modeling
//! 2. Lattice-based cryptography noise generation

use qbessel_rs::quantum::qkd::{BesselBeam, BeamProfile};
use qbessel_rs::quantum::lattice::{BesselNoiseSampler, BesselFunctionType, LatticeProblemGenerator};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Quantum Cryptography with Bessel Functions ===\n");

    // === QKD Bessel Beam Example ===
    println!("1. QKD Bessel Beam Analysis");
    println!("---------------------------");

    // Create a Bessel beam for QKD (wavelength = 632.8 nm, HeNe laser)
    let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);
    println!("Created Bessel beam: k_ρ = {}, order = {}, λ = {} nm", 
             beam.k_rho, beam.order, beam.wavelength * 1e9);

    // Calculate intensity at different radii
    println!("\nIntensity profile:");
    for r in [0.0, 0.1, 0.2, 0.3, 0.5] {
        let intensity = beam.intensity_profile(r)?;
        println!("  r = {:.1} → I = {:.4}", r, intensity);
    }

    // Check diffraction effects
    let distances = [1e-6, 1e-3, 1e-2, 1e-1];
    println!("\nDiffraction analysis:");
    for &dist in &distances {
        let has_diffraction = beam.diffraction_check(dist);
        println!("  Distance {:.0e} m → Diffraction: {}", dist, has_diffraction);
    }

    // Calculate effective radius for 90% power containment
    let eff_radius = beam.effective_radius(0.9, 2.0, 1000)?;
    println!("Effective radius (90% power): {:.3}", eff_radius);

    // Estimate QBER contribution
    let qber = beam.qber_contribution(0.5, 0.1)?;
    println!("QBER contribution (detector r=0.5, misalignment=0.1): {:.4}", qber);

    // Generate beam profile for analysis
    let profile = BeamProfile::generate(&beam, 2.0, 200)?;
    let quality = profile.quality_assessment();
    println!("Beam quality score: {:.3}", quality);
    if let Some(fwhm) = profile.fwhm {
        println!("Full Width at Half Maximum: {:.3}", fwhm);
    }

    println!("\n");

    // === Lattice-Based Cryptography Example ===
    println!("2. Lattice-Based Cryptography Noise Generation");
    println!("----------------------------------------------");

    // Create Bessel noise sampler for LWE
    let lwe_sampler = BesselNoiseSampler::new(
        0.0,                              // J_0 Bessel function
        1.0,                              // Scale parameter
        BesselFunctionType::FirstKind,    // Use J_ν(x)
        3.2                               // Target standard deviation
    );

    println!("Created LWE noise sampler: order = {}, σ = {}", 
             lwe_sampler.order, lwe_sampler.sigma);

    // Generate some noise samples
    println!("\nLWE noise samples:");
    for i in 0..5 {
        let x = i as f64 * 0.5;
        let continuous = lwe_sampler.sample(x)?;
        let discrete = lwe_sampler.sample_discrete(x)?;
        println!("  x = {:.1} → continuous: {:.4}, discrete: {}", 
                 x, continuous, discrete);
    }

    // Estimate noise variance
    let variance = lwe_sampler.estimate_variance(1000, 10.0)?;
    println!("Estimated noise variance: {:.4}", variance);

    // Create NTRU noise sampler
    let ntru_sampler = BesselNoiseSampler::new(
        0.5,                              // J_{0.5} Bessel function
        1.0,                              // Scale parameter
        BesselFunctionType::SecondKind,   // Use Y_ν(x) (Neumann)
        1.0                               // Standard deviation
    );

    println!("\nNTRU coefficient generation:");
    for i in 0..10 {
        let x = (i as f64 + 1.0) * 0.3; // Ensure x > 0 for Neumann functions
        let coeff = ntru_sampler.sample_ntru(x)?;
        print!("{:2} ", coeff);
    }
    println!();

    // === LWE Problem Generation ===
    println!("\n3. LWE Problem Instance Generation");
    println!("----------------------------------");

    let dimension = 5;
    let modulus = 97;
    let generator = LatticeProblemGenerator::new(
        dimension, 
        lwe_sampler, 
        Some(modulus)
    );

    // Generate a secret vector
    let secret = vec![1, -1, 0, 1, -1];
    println!("Secret vector: {:?}", secret);

    // Generate LWE samples
    let samples = generator.generate_lwe_samples(&secret, 3)?;
    println!("Generated {} LWE samples:", samples.len());
    for (i, (a, b)) in samples.iter().enumerate() {
        println!("  Sample {}: a = {:?}, b = {}", i + 1, a, b);
    }

    // Estimate problem hardness
    let hardness = generator.estimate_hardness()?;
    println!("Estimated problem hardness: {:.2}", hardness);

    println!("\n");

    // === Security Analysis ===
    println!("4. Security Analysis");
    println!("-------------------");

    // Validate entropy of noise distribution
    let has_entropy = qbessel_rs::quantum::lattice::utils::validate_entropy(&ntru_sampler, 100)?;
    println!("NTRU sampler has sufficient entropy: {}", has_entropy);

    // Generate Bessel-based secret
    let bessel_secret = qbessel_rs::quantum::lattice::utils::generate_bessel_secret(
        8, 1.234, 2
    );
    println!("Bessel-generated secret: {:?}", bessel_secret);

    println!("\n=== Analysis Complete ===");
    println!("The quantum cryptography modules demonstrate:");
    println!("• QKD beam modeling with intensity profiles and diffraction analysis");
    println!("• Lattice-based noise generation for post-quantum cryptography");
    println!("• LWE problem instance generation with Bessel function noise");
    println!("• Security analysis tools for quantum-resistant systems");

    Ok(())
}
