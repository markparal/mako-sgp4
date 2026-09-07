//! Example demonstrating how to export a TLE to an OMM-KVN string using the mako-sgp4 crate.

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------
use mako_sgp4::{from_tle_string, to_omm_kvn_string};

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

    // Export the GP data from the SGP4 propagator to an OMM-KVN string
    let omm_kvn = to_omm_kvn_string(&tle_sgp4s);

    // Print the result
    println!("{}", omm_kvn);
}

// ----------
// Unit Tests
// ----------
