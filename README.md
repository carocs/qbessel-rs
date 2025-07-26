# Quantum-Bessel-RS (`qbessel-rs`)

*A high-performance Rust library for Bessel/Neumann functions with quantum cryptography applications.*

[![Lib.rs](https://img.shields.io/badge/lib.rs-science/math-blue)](https://lib.rs/science/math)
[![no_std](https://img.shields.io/badge/no__std-supported-green)](https://docs.rust-embedded.org/book/intro/no-std.html)
[![SIMD](https://img.shields.io/badge/SIMD-accelerated-yellow)](https://github.com/rust-lang/portable-simd)

## Features

### Core Mathematical Functions
| Function | Type | Algorithm | Use Case |
|----------|------|-----------|----------|
| `Jν(x)` | Bessel (1st kind) | Miller's method (small `x`), Asymptotic (large `x`) | QKD beam propagation |
| `Yν(x)` | Neumann (2nd kind) | Wronskian recurrence | Fiber eavesdrop detection |
| `Hν(x)` | Hankel | `Jν(x) ± iYν(x)` | Quantum scattering problems |
| `Iν(x)/Kν(x)` | Modified Bessel | Temme's method | QRNG thermal noise |

### Quantum Cryptography Extensions
- **QKD Utilities**
  ```rust
  use quantum_bessel::qkd::BesselBeam;
  let beam = BesselBeam::new(5.0);  // k_ρ = 5.0

- **Lattice-Based Noise**
  ```rust
  use quantum_bessel::lattice::BesselNoise;
  let noise = BesselNoise::new(0.5); // ν = 0.5

- **Quantum RNG Validation**
  ```rust
  use quantum_bessel::rng::test_entropy;
  test_entropy(bessel_samples)?;

- **Quantum Key Distribution**
  ```rust
  use quantum_bessel::qkd::{BesselBeam, BeamProfile};

  fn main() {
      let beam = BesselBeam::new(2.5);
      println!("Intensity at r=1.0: {}", beam.intensity(1.0));
  }