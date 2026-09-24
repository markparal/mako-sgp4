<p align="center">
  <img src="https://raw.githubusercontent.com/markparal/mako-sgp4/HEAD/assets/logo_dark.png" alt="mako-sgp4 logo" width="240">
</p>

# mako-sgp4

[![Build status](https://github.com/markparal/mako-sgp4/actions/workflows/ci.yml/badge.svg)](https://github.com/markparal/mako-sgp4/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/mako-sgp4.svg)](https://crates.io/crates/mako-sgp4)
[![Documentation](https://docs.rs/mako-sgp4/badge.svg)](https://docs.rs/mako-sgp4)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust crate to parse and propagate General Perturbation Element Sets (GPs) using Simplified Perturbations Models (SGP4 / SDP4). Both Two-Line Elements (TLEs) and Orbit Mean-Elements Messages (OMMs) are supported. This code is implemented using the theory and equations found in *History of Analytical Orbit Modeling in the U.S. Space Surveillance System* by Hoots et al. Practical implementation adjustments were made based on *Revisiting Spacetrack Report #3: Rev 3* by Vallado et al.

The name mako-sgp4 pays tribute to the Shortfin Mako Shark, the fastest shark species. The speed and efficiency of the SGP4 propagator make this an apt name.

## Accuracy

The SGP4 propagator is a mean orbital elements propagator. At any point in time, its accuracy to the true position of an orbiting body is typically on the order of hundreds of meters to kilometers. Thus, it is not intended for high-precision operations.

mako-sgp4 is verified against the standard Vallado test cases as well as custom test cases (found in `test/`). In all propagation cases, mako-sgp4 agrees with Vallado's SGP4 propagator to within 1 meter.

## Documentation

The full API is on [docs.rs](https://docs.rs/mako-sgp4). The crate requires Rust 1.85 or newer (`edition = "2024"`). The style guide and math spec live in the `docs/` directory.

```bash
# Unit tests, doctests, and default features (XML / JSON / CSV)
cargo test

# TLE and OMM KVN only
cargo test --no-default-features

# Rustdoc
cargo doc --no-deps --open
```

## Usage

Add mako-sgp4 to your projects:

```bash
cargo add mako-sgp4
```

Propagate a TLE to its epoch with this example:

```rust
use mako_sgp4::{from_tle_string, sgp4_prop_delta};

fn main() {
    // Define the TLE string
    let tle_string = "\
    ISS (ZARYA)
    1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
    ";

    // Parse the TLE string into an SGP4 propagator. Returns a vector of SGP4 propagators
    let tle_sgp4s = from_tle_string(tle_string).unwrap();

    // Propagate the TLE propagator to the epoch time. Returns a StateVector struct.
    let state = sgp4_prop_delta(&tle_sgp4s[0], 0.0).unwrap();

    // Print the result in TEME coordinates.
    println!(
        "{}\nr_TEME = [{:.3}, {:.3}, {:.3}] km\nv_TEME = [{:.3}, {:.3}, {:.3}] km/s",
        tle_sgp4s[0].gp.common_name,
        state.r_x,
        state.r_y,
        state.r_z,
        state.v_x,
        state.v_y,
        state.v_z
    );
}
```

SGP4 propagation is accomplished with one of the following functions:
- `sgp4_prop_delta` - Propagates delta_t minutes from epoch
- `sgp4_prop_datetime` - Propagates to a specified UTC datetime

Example code lives in [`examples/`](examples/) and can be run with:

```bash
cargo run --example propagate_tle
```

### Features

TLE and OMM KVN parsing are always available. XML, JSON, and CSV OMM support are optional Cargo features but **on by default**.

| Feature | Formats |
| --- | --- |
| *(none / core)* | TLE, OMM KVN |
| `xml` | OMM XML |
| `json` | OMM JSON |
| `csv` | OMM CSV |

```bash
cargo add mako-sgp4                                           # TLE, KVN, XML, JSON, CSV
cargo add mako-sgp4 --no-default-features                     # TLE and KVN only, no 3rd party dependencies
cargo add mako-sgp4 --no-default-features --features json,csv # Exclude XML
```

## Future Work
- Add additional test cases
- Finish math spec
- Fit state data to GP
- Python wrapper

## References
- [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
- [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
- [Fundamentals of Astrodynamics Github Repository by Vallado et al](https://github.com/CelesTrak/fundamentals-of-astrodynamics?tab=readme-ov-file)
- [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
- [Space-Track](https://www.space-track.org/auth/login)
- [Celestrak](https://celestrak.org/)
