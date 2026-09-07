//! Example demonstrating how to export a TLE to OMM-CSV, import it, and propagate to epoch using the mako-sgp4 crate.

// ------------------
// External Libraries
// ------------------
use std::env;

// ------------------
// Internal Libraries
// ------------------
use mako_sgp4::{from_omm_csv_file, from_tle_string, sgp4_prop_delta, to_omm_csv_file};

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
    let tle_sgp4s = from_tle_string(tle_string).unwrap();

    // Export the GP data from the SGP4 propagator to an OMM-CSV file
    let path = env::temp_dir().join("omm_example.csv");
    let path_str = path.to_str().unwrap();
    to_omm_csv_file(&tle_sgp4s, path_str).unwrap();

    // Import the OMM-CSV file into an SGP4 propagator. Returns a vector of SGP4 propagators
    let omm_sgp4s = from_omm_csv_file(path_str).unwrap();

    // Propagate the OMM-CSV propagator to the epoch time. Returns a StateVector struct
    let state = sgp4_prop_delta(&omm_sgp4s[0], 0.0).unwrap();

    // Print the result in TEME coordinates
    println!(
        "{}\nr_TEME = [{:.3}, {:.3}, {:.3}] km\nv_TEME = [{:.3}, {:.3}, {:.3}] km/s",
        omm_sgp4s[0].gp.common_name,
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
