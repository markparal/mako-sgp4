//! Example demonstrating how to export a TLE to a file, import it, and propagate it using the mako-sgp4 crate.

// ------------------
// External Libraries
// ------------------
use std::env;

// ------------------
// Internal Libraries
// ------------------
use mako_sgp4::{from_tle_file, from_tle_string, sgp4_prop_delta, to_tle_file};

// -------
// Structs
// -------

// -----
// Enums
// -----

// ------
// Traits
// ------

// ---------
// Constants
// ---------

// ---------
// Functions
// ---------
fn main() {
    // Define the TLE string
    let tle_string = "\
    ISS (ZARYA)
    1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
    ";

    // Parse the TLE string into an SGP4 propagator. Returns a vector of SGP4 propagators
    let tle_sgp4s_original = from_tle_string(tle_string).unwrap();

    // Export the GP data from the SGP4 propagator to a TLE file
    let path = env::temp_dir().join("tle_example.tle");
    let path_str = path.to_str().unwrap();
    to_tle_file(&tle_sgp4s_original, path_str).unwrap();

    // Import the TLE file into an SGP4 propagator. Returns a vector of SGP4 propagators
    let tle_sgp4s_from_file = from_tle_file(path_str).unwrap();

    // Propagate the TLE propagator from the file to the epoch time. Returns a StateVector struct
    let state = sgp4_prop_delta(&tle_sgp4s_from_file[0], 0.0).unwrap();

    // Print the result in TEME coordinates
    println!(
        "{}\nr_TEME = [{:.3}, {:.3}, {:.3}] km\nv_TEME = [{:.3}, {:.3}, {:.3}] km/s",
        tle_sgp4s_from_file[0].gp.common_name,
        state.r_x,
        state.r_y,
        state.r_z,
        state.v_x,
        state.v_y,
        state.v_z
    );
}

// ----------
// Unit Tests
// ----------
