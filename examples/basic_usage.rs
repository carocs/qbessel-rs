//! Basic usage example for qbessel-rs
//! 
//! This example demonstrates the core Bessel function implementations.

use qbessel_rs::core::{bessel_j, neumann_y, hankel_h1, modified_bessel_i, modified_bessel_k};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Quantum-Bessel-RS Basic Usage Example ===\n");

    // Bessel functions of the first kind J_ν(x)
    println!("Bessel Functions J_ν(x):");
    println!("J₀(0.0) = {:.10}", bessel_j(0.0, 0.0)?);
    println!("J₀(1.0) = {:.10}", bessel_j(0.0, 1.0)?);
    println!("J₁(1.0) = {:.10}", bessel_j(1.0, 1.0)?);
    println!("J₂(2.0) = {:.10}", bessel_j(2.0, 2.0)?);
    println!();

    // Neumann functions (Bessel functions of the second kind) Y_ν(x)
    println!("Neumann Functions Y_ν(x):");
    println!("Y₀(1.0) = {:.10}", neumann_y(0.0, 1.0)?);
    println!("Y₁(2.0) = {:.10}", neumann_y(1.0, 2.0)?);
    println!("Y₂(3.0) = {:.10}", neumann_y(2.0, 3.0)?);
    println!();

    // Hankel functions H_ν(x) = J_ν(x) ± iY_ν(x)
    println!("Hankel Functions H_ν(x):");
    let h1 = hankel_h1(0.0, 1.0)?;
    println!("H₁₀(1.0) = {:.6} + {:.6}i", h1.re, h1.im);
    
    let h1_1 = hankel_h1(1.0, 2.0)?;
    println!("H₁₁(2.0) = {:.6} + {:.6}i", h1_1.re, h1_1.im);
    println!();

    // Modified Bessel functions I_ν(x) and K_ν(x)
    println!("Modified Bessel Functions:");
    println!("I₀(0.0) = {:.10}", modified_bessel_i(0.0, 0.0)?);
    println!("I₀(1.0) = {:.10}", modified_bessel_i(0.0, 1.0)?);
    println!("I₁(1.0) = {:.10}", modified_bessel_i(1.0, 1.0)?);
    println!();
    
    println!("K₀(1.0) = {:.10}", modified_bessel_k(0.0, 1.0)?);
    println!("K₁(1.0) = {:.10}", modified_bessel_k(1.0, 1.0)?);
    println!();

    // Demonstrate quantum cryptography applications
    println!("=== Quantum Cryptography Applications ===");
    
    // QKD beam propagation example
    println!("\nQKD Beam Propagation:");
    println!("For a Bessel beam with k_ρ = 5.0:");
    let k_rho = 5.0;
    let radius = 2.0;
    let intensity = bessel_j(0.0, k_rho * radius)?.powi(2);
    println!("Intensity at r={}: I(r) = |J₀({:.1} × {:.1})|² = {:.6}", 
             radius, k_rho, radius, intensity);
    
    // Fiber eavesdrop detection example
    println!("\nFiber Eavesdrop Detection:");
    println!("Neumann function for fiber mode analysis:");
    let fiber_param = 2.4;
    let y_val = neumann_y(1.0, fiber_param)?;
    println!("Y₁({:.1}) = {:.6} (used in fiber mode calculations)", fiber_param, y_val);
    
    // Thermal noise in QRNG
    println!("\nQRNG Thermal Noise Analysis:");
    println!("Modified Bessel functions for thermal noise modeling:");
    let thermal_param = 1.5;
    let i_val = modified_bessel_i(0.0, thermal_param)?;
    let k_val = modified_bessel_k(0.0, thermal_param)?;
    println!("I₀({:.1}) = {:.6}", thermal_param, i_val);
    println!("K₀({:.1}) = {:.6}", thermal_param, k_val);
    println!("Ratio I₀/K₀ = {:.6} (thermal noise characteristic)", i_val / k_val);

    println!("\n=== Implementation Features ===");
    println!("✓ Miller's algorithm for small arguments");
    println!("✓ Asymptotic expansions for large arguments");
    println!("✓ Wronskian recurrence for Neumann functions");
    println!("✓ Temme's method for modified Bessel functions");
    println!("✓ Complex Hankel function support");
    println!("✓ Error handling for invalid arguments");
    println!("✓ no_std compatibility (with libm feature)");

    Ok(())
}
