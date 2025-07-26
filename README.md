# qbessel-rs

A high-performance Rust crate for Bessel (Jν), Neumann (Yν), and related functions, optimized for quantum cryptography applications.

## Features

### Core Functions
- **Jν(x)** - Bessel functions of the first kind (Taylor series, Miller's algorithm)
- **Yν(x)** - Neumann functions (second kind) with Wronskian recurrence  
- **Hν(x)** - Hankel functions for scattering problems (Hν = Jν ± iYν)
- **Iν(x), Kν(x)** - Modified Bessel functions using Temme's method

### Quantum Cryptography Extensions
- **qkd module**: Bessel beam utilities for QKD systems
  - Intensity profile calculations
  - Diffraction analysis
  - QBER (Quantum Bit Error Rate) estimation
  - Beam quality assessment
- **lattice module**: Noise samplers for post-quantum cryptography
  - LWE (Learning With Errors) noise generation
  - NTRU coefficient sampling
  - Lattice problem instance generation
- **Security analysis**: Entropy validation and hardness estimation

### Performance & Compatibility
- Optimized algorithms for different argument ranges
- Comprehensive error handling with `QuantumBesselError`
- `no_std` support with `libm` feature
- Ready for SIMD acceleration (framework in place)

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
qbessel-rs = "0.1.0"

# Optional features
[features]
default = ["std"]
std = []                           # Standard library support
no-std = ["libm"]                  # Embedded/no_std targets
```

### Basic Bessel Functions

```rust
use qbessel_rs::core::{bessel_j, neumann_y, modified_bessel_i, modified_bessel_k};

// Bessel functions of the first kind
let j0_1 = bessel_j(0.0, 1.0)?;  // J₀(1) ≈ 0.7652
let j1_2 = bessel_j(1.0, 2.0)?;  // J₁(2) ≈ 0.5767

// Neumann functions (second kind)
let y0_1 = neumann_y(0.0, 1.0)?; // Y₀(1) ≈ 0.0883
let y1_2 = neumann_y(1.0, 2.0)?; // Y₁(2) ≈ -0.2900

// Modified Bessel functions
let i0_1 = modified_bessel_i(0.0, 1.0)?; // I₀(1) ≈ 1.2661
let k0_1 = modified_bessel_k(0.0, 1.0)?; // K₀(1) ≈ 0.4210
```

### QKD Beam Simulation

```rust
use qbessel_rs::quantum::qkd::{BesselBeam, BeamProfile};

// Create a Bessel beam for QKD (HeNe laser, 632.8 nm)
let beam = BesselBeam::new(5.0, 0.0, 632.8e-9);

// Calculate intensity profile
let intensity = beam.intensity_profile(2.0)?; // |J₀(5.0 × 2.0)|²

// Check for diffraction effects
let has_diffraction = beam.diffraction_check(1e-3); // At 1mm distance

// Estimate QBER contribution
let qber = beam.qber_contribution(0.5, 0.1)?; // detector radius, misalignment

// Generate beam profile for analysis
let profile = BeamProfile::generate(&beam, 2.0, 200)?;
let quality = profile.quality_assessment();
```

### Lattice-Based Cryptography

```rust
use qbessel_rs::quantum::lattice::{BesselNoiseSampler, BesselFunctionType, LatticeProblemGenerator};

// Create noise sampler for LWE
let lwe_sampler = BesselNoiseSampler::new(
    0.0,                              // J₀ Bessel function
    1.0,                              // Scale parameter
    BesselFunctionType::FirstKind,    // Use J_ν(x)
    3.2                               // Target standard deviation
);

// Generate noise samples
let continuous_noise = lwe_sampler.sample(1.0)?;
let discrete_noise = lwe_sampler.sample_discrete(1.0)?;

// Create NTRU noise sampler using Neumann functions
let ntru_sampler = BesselNoiseSampler::new(
    0.5,                              // J₀.₅ Bessel function
    1.0,                              // Scale parameter
    BesselFunctionType::SecondKind,   // Use Y_ν(x) (Neumann)
    1.0                               // Standard deviation
);

// Generate NTRU coefficients in {-1, 0, 1}
let ntru_coeff = ntru_sampler.sample_ntru(1.5)?;

// Generate LWE problem instances
let generator = LatticeProblemGenerator::new(5, lwe_sampler, Some(97));
let secret = vec![1, -1, 0, 1, -1];
let samples = generator.generate_lwe_samples(&secret, 3)?;
```

## Examples

Run the included examples to see the library in action:

```bash
# Basic Bessel function usage
cargo run --example basic_usage

# Quantum cryptography applications
cargo run --example quantum_crypto
```

The quantum cryptography example demonstrates:
- QKD beam modeling with intensity profiles and diffraction analysis
- Lattice-based noise generation for post-quantum cryptography
- LWE problem instance generation with Bessel function noise
- Security analysis tools for quantum-resistant systems

## Applications

### Quantum Key Distribution (QKD)
- Bessel beam propagation modeling
- Intensity profile calculations for beam characterization
- Diffraction analysis for free-space QKD systems
- QBER estimation for system optimization

### Post-Quantum Cryptography
- LWE noise generation with quantum-resistant properties
- NTRU polynomial coefficient sampling
- Lattice problem hardness estimation
- Entropy validation for cryptographic security

### Scientific Computing
- High-precision Bessel function calculations
- Complex analysis with Hankel functions
- Modified Bessel functions for thermal noise modeling
- Numerical methods for physics and optics applications

## Implementation Details

### Algorithms Used
- **Miller's algorithm** for stable computation of Bessel functions
- **Asymptotic expansions** for large arguments
- **Wronskian recurrence** for Neumann functions
- **Temme's method** for modified Bessel functions
- **Taylor series** for small arguments with high precision

### Error Handling
The library uses comprehensive error handling with `QuantumBesselError`:
- `InvalidInput`: Invalid parameters or arguments
- `ComputationError`: Numerical computation failures
- `BesselError`: Core Bessel function errors

### Performance
- Optimized for different argument ranges
- Efficient algorithms selected based on input values
- Ready for SIMD acceleration (framework implemented)
- Minimal allocations for embedded use

## Testing

Run the test suite:

```bash
# Unit tests
cargo test

# Integration tests
cargo test --test integration_tests

# All tests with output
cargo test -- --nocapture
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Roadmap

- [ ] SIMD acceleration using `std::simd`
- [ ] GPU support with CUDA backend
- [ ] Precomputed tables for common orders
- [ ] Additional quantum cryptography protocols
- [ ] Benchmarking against GSL and SciPy
- [ ] WebAssembly support
