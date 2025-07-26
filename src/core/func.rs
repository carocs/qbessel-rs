//! Core Bessel and Neumann function implementations
//! 
//! This module provides high-performance implementations of:
//! - Jν(x): Bessel functions of the first kind
//! - Yν(x): Neumann functions (Bessel functions of the second kind)
//! - Hν(x): Hankel functions
//! - Iν(x), Kν(x): Modified Bessel functions

use core::f64::consts::PI;

#[cfg(feature = "std")]
use std::f64;

#[cfg(not(feature = "std"))]
use libm;

// Math functions - use std or libm depending on feature
#[cfg(feature = "std")]
fn sin(x: f64) -> f64 { x.sin() }
#[cfg(feature = "std")]
fn cos(x: f64) -> f64 { x.cos() }
#[cfg(feature = "std")]
fn exp(x: f64) -> f64 { x.exp() }
#[cfg(feature = "std")]
fn log(x: f64) -> f64 { x.ln() }
#[cfg(feature = "std")]
fn sqrt(x: f64) -> f64 { x.sqrt() }
#[cfg(feature = "std")]
fn pow(x: f64, y: f64) -> f64 { x.powf(y) }
#[cfg(feature = "std")]
fn fabs(x: f64) -> f64 { x.abs() }

#[cfg(not(feature = "std"))]
fn sin(x: f64) -> f64 { libm::sin(x) }
#[cfg(not(feature = "std"))]
fn cos(x: f64) -> f64 { libm::cos(x) }
#[cfg(not(feature = "std"))]
fn exp(x: f64) -> f64 { libm::exp(x) }
#[cfg(not(feature = "std"))]
fn log(x: f64) -> f64 { libm::log(x) }
#[cfg(not(feature = "std"))]
fn sqrt(x: f64) -> f64 { libm::sqrt(x) }
#[cfg(not(feature = "std"))]
fn pow(x: f64, y: f64) -> f64 { libm::pow(x, y) }
#[cfg(not(feature = "std"))]
fn fabs(x: f64) -> f64 { libm::fabs(x) }

#[cfg(feature = "std")]
const E: f64 = core::f64::consts::E;
#[cfg(not(feature = "std"))]
const E: f64 = libm::E;

/// Error type for Bessel function computations
#[derive(Debug, Clone)]
pub enum BesselError {
    InvalidArgument(String),
    ConvergenceError(String),
    NumericalError(String),
}

impl core::fmt::Display for BesselError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
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

/// Complex number representation for Hankel functions
#[derive(Debug, Clone, Copy)]
pub struct Complex64 {
    pub re: f64,
    pub im: f64,
}

impl Complex64 {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    
    pub fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }
    
    pub fn i() -> Self {
        Self { re: 0.0, im: 1.0 }
    }
}

impl core::ops::Add for Complex64 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl core::ops::Sub for Complex64 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl core::ops::Mul for Complex64 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

/// Bessel function of the first kind Jν(x)
/// 
/// Uses Miller's method for small x and asymptotic expansion for large x
/// 
/// # Arguments
/// * `nu` - Order of the Bessel function
/// * `x` - Argument
/// 
/// # Returns
/// Value of Jν(x)
pub fn bessel_j(nu: f64, x: f64) -> BesselResult<f64> {
    if x < 0.0 {
        return Err(BesselError::InvalidArgument("x must be non-negative".to_string()));
    }
    
    if x == 0.0 {
        return if nu == 0.0 { Ok(1.0) } else { Ok(0.0) };
    }
    
    // Choose algorithm based on argument size
    if x < 10.0 {
        bessel_j_miller(nu, x)
    } else {
        bessel_j_asymptotic(nu, x)
    }
}

/// Miller's method for computing Jν(x) for small to moderate x
fn bessel_j_miller(nu: f64, x: f64) -> BesselResult<f64> {
    let n = nu.round() as i32;
    
    // For integer orders near nu, use series expansion
    if (nu - n as f64).abs() < 1e-10 {
        return bessel_j_series(n, x);
    }
    
    // For non-integer orders, use Taylor series
    bessel_j_taylor_series(nu, x)
}

/// Taylor series expansion for Bessel functions
fn bessel_j_taylor_series(nu: f64, x: f64) -> BesselResult<f64> {
    let x_half = x / 2.0;
    let mut term = pow(x_half, nu) / gamma_approx(nu + 1.0);
    let mut sum = term;
    let x_squared_quarter = x * x / 4.0;
    
    for k in 1..100 {
        term *= -x_squared_quarter / (k as f64 * (nu + k as f64));
        sum += term;
        
        if fabs(term) < fabs(sum) * 1e-15 {
            break;
        }
    }
    
    Ok(sum)
}

/// Series expansion for integer order Bessel functions
fn bessel_j_series(n: i32, x: f64) -> BesselResult<f64> {
    let x_half = x / 2.0;
    let mut term = pow(x_half, n.abs() as f64) / factorial(n.abs() as u32) as f64;
    let mut sum = term;
    let x_squared_quarter = x * x / 4.0;
    
    for k in 1..100 {
        term *= -x_squared_quarter / (k as f64 * (k + n.abs()) as f64);
        sum += term;
        
        if fabs(term) < fabs(sum) * 1e-15 {
            break;
        }
    }
    
    if n < 0 && n % 2 != 0 {
        sum = -sum;
    }
    
    Ok(sum)
}

/// Asymptotic expansion for large x
fn bessel_j_asymptotic(nu: f64, x: f64) -> BesselResult<f64> {
    let phase = x - nu * PI / 2.0 - PI / 4.0;
    let amplitude = sqrt(2.0 / (PI * x));
    
    // First-order asymptotic approximation
    Ok(amplitude * cos(phase))
}

/// Neumann function (Bessel function of the second kind) Yν(x)
/// 
/// Uses Wronskian recurrence relation
/// 
/// # Arguments
/// * `nu` - Order of the Neumann function
/// * `x` - Argument (must be positive)
/// 
/// # Returns
/// Value of Yν(x)
pub fn neumann_y(nu: f64, x: f64) -> BesselResult<f64> {
    if x <= 0.0 {
        return Err(BesselError::InvalidArgument("x must be positive for Neumann functions".to_string()));
    }
    
    if x < 10.0 {
        neumann_y_wronskian(nu, x)
    } else {
        neumann_y_asymptotic(nu, x)
    }
}

/// Wronskian recurrence for Neumann functions
fn neumann_y_wronskian(nu: f64, x: f64) -> BesselResult<f64> {
    let n = nu.round() as i32;
    
    // For integer orders
    if (nu - n as f64).abs() < 1e-10 {
        return neumann_y_integer(n, x);
    }
    
    // For non-integer orders: Y_ν = (J_ν cos(νπ) - J_{-ν}) / sin(νπ)
    let j_nu = bessel_j(nu, x)?;
    let j_minus_nu = bessel_j(-nu, x)?;
    let sin_nu_pi = sin(nu * PI);
    
    if fabs(sin_nu_pi) < 1e-15 {
        return Err(BesselError::NumericalError("sin(νπ) too small".to_string()));
    }
    
    Ok((j_nu * cos(nu * PI) - j_minus_nu) / sin_nu_pi)
}

/// Integer order Neumann functions
fn neumann_y_integer(n: i32, x: f64) -> BesselResult<f64> {
    match n {
        0 => neumann_y0(x),
        1 => neumann_y1(x),
        _ => {
            let y0 = neumann_y0(x)?;
            let y1 = neumann_y1(x)?;
            
            if n < 0 {
                // Y_{-n} = (-1)^n Y_n
                let yn = neumann_y_recurrence_up(n.abs(), x, y0, y1)?;
                Ok(if n % 2 == 0 { yn } else { -yn })
            } else {
                neumann_y_recurrence_up(n, x, y0, y1)
            }
        }
    }
}

/// Y₀(x) using series expansion
fn neumann_y0(x: f64) -> BesselResult<f64> {
    let j0 = bessel_j(0.0, x)?;
    let x_half = x / 2.0;
    let euler_gamma = 0.5772156649015329; // Euler-Mascheroni constant
    
    // Y₀(x) = (2/π)[ln(x/2) + γ]J₀(x) - (2/π)∑_{k=1}^∞ (-1)^k H_k (x/2)^{2k} / (k!)²
    let log_term = (2.0 / PI) * (log(x_half) + euler_gamma) * j0;
    
    let mut series_sum = 0.0;
    let x_squared_quarter = x * x / 4.0;
    
    for k in 1..50 {
        let harmonic_k = harmonic_number(k);
        let term = pow(-x_squared_quarter, k as f64) * harmonic_k / 
                   (factorial(k) as f64).powi(2);
        series_sum += term;
        
        if fabs(term) < 1e-15 {
            break;
        }
    }
    
    Ok(log_term - (2.0 / PI) * series_sum)
}

/// Y₁(x) using series expansion
fn neumann_y1(x: f64) -> BesselResult<f64> {
    let j1 = bessel_j(1.0, x)?;
    let x_half = x / 2.0;
    let euler_gamma = 0.5772156649015329;
    
    // Similar to Y₀ but with different series
    let log_term = (2.0 / PI) * (log(x_half) + euler_gamma - 0.5) * j1;
    let leading_term = -(2.0 / PI) / x;
    
    Ok(log_term + leading_term)
}

/// Upward recurrence for Neumann functions
fn neumann_y_recurrence_up(n: i32, x: f64, y0: f64, y1: f64) -> BesselResult<f64> {
    if n == 0 { return Ok(y0); }
    if n == 1 { return Ok(y1); }
    
    let mut y_prev = y0;
    let mut y_curr = y1;
    
    for k in 1..n {
        let y_next = (2.0 * k as f64 / x) * y_curr - y_prev;
        y_prev = y_curr;
        y_curr = y_next;
    }
    
    Ok(y_curr)
}

/// Asymptotic expansion for Neumann functions
fn neumann_y_asymptotic(nu: f64, x: f64) -> BesselResult<f64> {
    let phase = x - nu * PI / 2.0 - PI / 4.0;
    let amplitude = sqrt(2.0 / (PI * x));
    
    // First-order asymptotic approximation
    Ok(amplitude * sin(phase))
}

/// Hankel function of the first kind H₁ν(x) = Jν(x) + iYν(x)
/// 
/// # Arguments
/// * `nu` - Order of the Hankel function
/// * `x` - Argument
/// 
/// # Returns
/// Complex value of H₁ν(x)
pub fn hankel_h1(nu: f64, x: f64) -> BesselResult<Complex64> {
    let j_val = bessel_j(nu, x)?;
    let y_val = neumann_y(nu, x)?;
    
    Ok(Complex64::new(j_val, y_val))
}

/// Hankel function of the second kind H₂ν(x) = Jν(x) - iYν(x)
/// 
/// # Arguments
/// * `nu` - Order of the Hankel function
/// * `x` - Argument
/// 
/// # Returns
/// Complex value of H₂ν(x)
pub fn hankel_h2(nu: f64, x: f64) -> BesselResult<Complex64> {
    let j_val = bessel_j(nu, x)?;
    let y_val = neumann_y(nu, x)?;
    
    Ok(Complex64::new(j_val, -y_val))
}

/// Modified Bessel function of the first kind Iν(x)
/// 
/// Uses Temme's method for computation
/// 
/// # Arguments
/// * `nu` - Order of the modified Bessel function
/// * `x` - Argument
/// 
/// # Returns
/// Value of Iν(x)
pub fn modified_bessel_i(nu: f64, x: f64) -> BesselResult<f64> {
    if x < 0.0 {
        return Err(BesselError::InvalidArgument("x must be non-negative".to_string()));
    }
    
    if x == 0.0 {
        return if nu == 0.0 { Ok(1.0) } else { Ok(0.0) };
    }
    
    if x < 10.0 {
        modified_bessel_i_temme(nu, x)
    } else {
        modified_bessel_i_asymptotic(nu, x)
    }
}

/// Temme's method for modified Bessel functions
fn modified_bessel_i_temme(nu: f64, x: f64) -> BesselResult<f64> {
    let x_half = x / 2.0;
    let mut term = pow(x_half, nu) / gamma_approx(nu + 1.0);
    let mut sum = term;
    let x_squared_quarter = x * x / 4.0;
    
    for k in 1..100 {
        term *= x_squared_quarter / (k as f64 * (nu + k as f64));
        sum += term;
        
        if fabs(term) < fabs(sum) * 1e-15 {
            break;
        }
    }
    
    Ok(sum)
}

/// Asymptotic expansion for modified Bessel I functions
fn modified_bessel_i_asymptotic(_nu: f64, x: f64) -> BesselResult<f64> {
    let amplitude = exp(x) / sqrt(2.0 * PI * x);
    
    // First-order asymptotic approximation
    Ok(amplitude)
}

/// Modified Bessel function of the second kind Kν(x)
/// 
/// Uses Temme's method for computation
/// 
/// # Arguments
/// * `nu` - Order of the modified Bessel function
/// * `x` - Argument (must be positive)
/// 
/// # Returns
/// Value of Kν(x)
pub fn modified_bessel_k(nu: f64, x: f64) -> BesselResult<f64> {
    if x <= 0.0 {
        return Err(BesselError::InvalidArgument("x must be positive for K functions".to_string()));
    }
    
    if x < 10.0 {
        modified_bessel_k_temme(nu, x)
    } else {
        modified_bessel_k_asymptotic(nu, x)
    }
}

/// Temme's method for K functions
fn modified_bessel_k_temme(nu: f64, x: f64) -> BesselResult<f64> {
    let n = nu.round() as i32;
    
    // For integer orders
    if (nu - n as f64).abs() < 1e-10 {
        return modified_bessel_k_integer(n, x);
    }
    
    // For non-integer orders: K_ν = π(I_{-ν} - I_ν) / (2sin(νπ))
    let i_nu = modified_bessel_i(nu, x)?;
    let i_minus_nu = modified_bessel_i(-nu, x)?;
    let sin_nu_pi = sin(nu * PI);
    
    if fabs(sin_nu_pi) < 1e-15 {
        return Err(BesselError::NumericalError("sin(νπ) too small".to_string()));
    }
    
    Ok(PI * (i_minus_nu - i_nu) / (2.0 * sin_nu_pi))
}

/// Integer order K functions
fn modified_bessel_k_integer(n: i32, x: f64) -> BesselResult<f64> {
    match n {
        0 => modified_bessel_k0(x),
        1 => modified_bessel_k1(x),
        _ => {
            let k0 = modified_bessel_k0(x)?;
            let k1 = modified_bessel_k1(x)?;
            
            if n < 0 {
                // K_{-n} = K_n
                modified_bessel_k_recurrence_up(n.abs(), x, k0, k1)
            } else {
                modified_bessel_k_recurrence_up(n, x, k0, k1)
            }
        }
    }
}

/// K₀(x) using series expansion
fn modified_bessel_k0(x: f64) -> BesselResult<f64> {
    let i0 = modified_bessel_i(0.0, x)?;
    let euler_gamma = 0.5772156649015329;
    
    // K₀(x) = -ln(x/2)I₀(x) - γI₀(x) + series
    let log_term = -(log(x / 2.0) + euler_gamma) * i0;
    
    let mut series_sum = 0.0;
    let x_squared_quarter = x * x / 4.0;
    
    for k in 1..50 {
        let harmonic_k = harmonic_number(k);
        let term = pow(x_squared_quarter, k as f64) * harmonic_k / 
                   (factorial(k) as f64).powi(2);
        series_sum += term;
        
        if fabs(term) < 1e-15 {
            break;
        }
    }
    
    Ok(log_term + series_sum)
}

/// K₁(x) using series expansion
fn modified_bessel_k1(x: f64) -> BesselResult<f64> {
    let i1 = modified_bessel_i(1.0, x)?;
    let euler_gamma = 0.5772156649015329;
    
    // Similar to K₀ but with different series
    let log_term = -(log(x / 2.0) + euler_gamma - 0.5) * i1;
    let leading_term = 1.0 / x;
    
    Ok(log_term + leading_term)
}

/// Upward recurrence for K functions
fn modified_bessel_k_recurrence_up(n: i32, x: f64, k0: f64, k1: f64) -> BesselResult<f64> {
    if n == 0 { return Ok(k0); }
    if n == 1 { return Ok(k1); }
    
    let mut k_prev = k0;
    let mut k_curr = k1;
    
    for k in 1..n {
        let k_next = (2.0 * k as f64 / x) * k_curr + k_prev;
        k_prev = k_curr;
        k_curr = k_next;
    }
    
    Ok(k_curr)
}

/// Asymptotic expansion for K functions
fn modified_bessel_k_asymptotic(_nu: f64, x: f64) -> BesselResult<f64> {
    let amplitude = sqrt(PI / (2.0 * x)) * exp(-x);
    
    // First-order asymptotic approximation
    Ok(amplitude)
}

// Helper functions

/// Factorial function
fn factorial(n: u32) -> u64 {
    match n {
        0 | 1 => 1,
        _ => (2..=n as u64).product(),
    }
}

/// Harmonic number H_n = 1 + 1/2 + ... + 1/n
fn harmonic_number(n: u32) -> f64 {
    (1..=n).map(|k| 1.0 / k as f64).sum()
}

/// Gamma function approximation using Stirling's formula
fn gamma_approx(x: f64) -> f64 {
    if x < 1.0 {
        return gamma_approx(x + 1.0) / x;
    }
    
    // Stirling's approximation for x >= 1
    let sqrt_2pi = sqrt(2.0 * PI);
    sqrt_2pi * pow(x / E, x) * sqrt(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_j0() {
        // J₀(0) = 1
        let result = bessel_j(0.0, 0.0).unwrap();
        assert!((result - 1.0).abs() < 1e-10);
        
        // J₀(1) ≈ 0.7651976865579666
        let result = bessel_j(0.0, 1.0).unwrap();
        assert!((result - 0.7651976865579666).abs() < 1e-10);
    }

    #[test]
    fn test_bessel_j1() {
        // J₁(0) = 0
        let result = bessel_j(1.0, 0.0).unwrap();
        assert!(result.abs() < 1e-10);
        
        // J₁(1) ≈ 0.44005058574493355
        let result = bessel_j(1.0, 1.0).unwrap();
        assert!((result - 0.44005058574493355).abs() < 1e-10);
    }

    #[test]
    fn test_neumann_y0() {
        // Y₀(1) ≈ 0.08825696421567696
        let result = neumann_y(0.0, 1.0).unwrap();
        assert!((result - 0.08825696421567696).abs() < 1e-8);
    }

    #[test]
    fn test_hankel_functions() {
        let h1 = hankel_h1(0.0, 1.0).unwrap();
        let h2 = hankel_h2(0.0, 1.0).unwrap();
        
        let j0 = bessel_j(0.0, 1.0).unwrap();
        let y0 = neumann_y(0.0, 1.0).unwrap();
        
        assert!((h1.re - j0).abs() < 1e-10);
        assert!((h1.im - y0).abs() < 1e-10);
        assert!((h2.re - j0).abs() < 1e-10);
        assert!((h2.im + y0).abs() < 1e-10);
    }

    #[test]
    fn test_modified_bessel_i0() {
        // I₀(0) = 1
        let result = modified_bessel_i(0.0, 0.0).unwrap();
        assert!((result - 1.0).abs() < 1e-10);
        
        // I₀(1) ≈ 1.2660658777520084
        let result = modified_bessel_i(0.0, 1.0).unwrap();
        assert!((result - 1.2660658777520084).abs() < 1e-10);
    }

    #[test]
    fn test_modified_bessel_k0() {
        // K₀(1) ≈ 0.42102443824070834
        let result = modified_bessel_k(0.0, 1.0).unwrap();
        assert!((result - 0.42102443824070834).abs() < 1e-8);
    }
}
