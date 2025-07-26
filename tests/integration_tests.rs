//! Integration tests for qbessel-rs
//! 
//! These tests verify the core Bessel function implementations against known values.

use qbessel_rs::core::{bessel_j, neumann_y, hankel_h1, hankel_h2, modified_bessel_i, modified_bessel_k};

#[test]
fn test_bessel_j_known_values() {
    // Test J₀(0) = 1
    let result = bessel_j(0.0, 0.0).unwrap();
    assert!((result - 1.0).abs() < 1e-10, "J₀(0) should be 1, got {}", result);
    
    // Test J₁(0) = 0
    let result = bessel_j(1.0, 0.0).unwrap();
    assert!(result.abs() < 1e-10, "J₁(0) should be 0, got {}", result);
    
    // Test J₀(1) ≈ 0.7651976865579666
    let result = bessel_j(0.0, 1.0).unwrap();
    let expected = 0.7651976865579666;
    assert!((result - expected).abs() < 1e-10, "J₀(1) should be {}, got {}", expected, result);
    
    // Test J₁(1) ≈ 0.44005058574493355
    let result = bessel_j(1.0, 1.0).unwrap();
    let expected = 0.44005058574493355;
    assert!((result - expected).abs() < 1e-10, "J₁(1) should be {}, got {}", expected, result);
}

#[test]
fn test_neumann_y_known_values() {
    // Test Y₀(1) ≈ 0.08825696421567696
    let result = neumann_y(0.0, 1.0).unwrap();
    let expected = 0.08825696421567696;
    assert!((result - expected).abs() < 1e-6, "Y₀(1) should be approximately {}, got {}", expected, result);
    
    // Test Y₁(1) ≈ -0.7812128213002887
    let result = neumann_y(1.0, 1.0).unwrap();
    let expected = -0.7812128213002887;
    assert!((result - expected).abs() < 1e-6, "Y₁(1) should be approximately {}, got {}", expected, result);
}

#[test]
fn test_hankel_functions() {
    let x = 1.0;
    let nu = 0.0;
    
    let h1 = hankel_h1(nu, x).unwrap();
    let h2 = hankel_h2(nu, x).unwrap();
    
    let j_val = bessel_j(nu, x).unwrap();
    let y_val = neumann_y(nu, x).unwrap();
    
    // H₁ = J + iY
    assert!((h1.re - j_val).abs() < 1e-10, "H₁ real part should equal J");
    assert!((h1.im - y_val).abs() < 1e-10, "H₁ imaginary part should equal Y");
    
    // H₂ = J - iY
    assert!((h2.re - j_val).abs() < 1e-10, "H₂ real part should equal J");
    assert!((h2.im + y_val).abs() < 1e-10, "H₂ imaginary part should equal -Y");
}

#[test]
fn test_modified_bessel_i_known_values() {
    // Test I₀(0) = 1
    let result = modified_bessel_i(0.0, 0.0).unwrap();
    assert!((result - 1.0).abs() < 1e-10, "I₀(0) should be 1, got {}", result);
    
    // Test I₁(0) = 0
    let result = modified_bessel_i(1.0, 0.0).unwrap();
    assert!(result.abs() < 1e-10, "I₁(0) should be 0, got {}", result);
    
    // Test I₀(1) ≈ 1.2660658777520084
    let result = modified_bessel_i(0.0, 1.0).unwrap();
    let expected = 1.2660658777520084;
    assert!((result - expected).abs() < 1e-10, "I₀(1) should be {}, got {}", expected, result);
}

#[test]
fn test_modified_bessel_k_known_values() {
    // Test K₀(1) ≈ 0.42102443824070834
    let result = modified_bessel_k(0.0, 1.0).unwrap();
    let expected = 0.42102443824070834;
    assert!((result - expected).abs() < 1e-6, "K₀(1) should be approximately {}, got {}", expected, result);
    
    // Test K₁(1) ≈ 0.6019072301972346
    let result = modified_bessel_k(1.0, 1.0).unwrap();
    let expected = 0.6019072301972346;
    assert!((result - expected).abs() < 1e-6, "K₁(1) should be approximately {}, got {}", expected, result);
}

#[test]
fn test_error_conditions() {
    // Test negative x for regular Bessel functions
    assert!(bessel_j(0.0, -1.0).is_err(), "J₀(-1) should return an error");
    
    // Test non-positive x for Neumann functions
    assert!(neumann_y(0.0, 0.0).is_err(), "Y₀(0) should return an error");
    assert!(neumann_y(0.0, -1.0).is_err(), "Y₀(-1) should return an error");
    
    // Test negative x for modified Bessel I functions
    assert!(modified_bessel_i(0.0, -1.0).is_err(), "I₀(-1) should return an error");
    
    // Test non-positive x for modified Bessel K functions
    assert!(modified_bessel_k(0.0, 0.0).is_err(), "K₀(0) should return an error");
    assert!(modified_bessel_k(0.0, -1.0).is_err(), "K₀(-1) should return an error");
}

#[test]
fn test_asymptotic_behavior() {
    // Test large x behavior - should use asymptotic approximations
    let large_x = 50.0;
    
    // These should not panic and should return reasonable values
    let j0_large = bessel_j(0.0, large_x).unwrap();
    let y0_large = neumann_y(0.0, large_x).unwrap();
    let i0_large = modified_bessel_i(0.0, large_x).unwrap();
    let k0_large = modified_bessel_k(0.0, large_x).unwrap();
    
    // Basic sanity checks
    assert!(j0_large.is_finite(), "J₀({}) should be finite", large_x);
    assert!(y0_large.is_finite(), "Y₀({}) should be finite", large_x);
    assert!(i0_large.is_finite(), "I₀({}) should be finite", large_x);
    assert!(k0_large.is_finite(), "K₀({}) should be finite", large_x);
    
    // I₀ should be large and positive for large x
    assert!(i0_large > 1e10, "I₀({}) should be very large", large_x);
    
    // K₀ should be small and positive for large x
    assert!(k0_large > 0.0 && k0_large < 1e-10, "K₀({}) should be very small", large_x);
}
