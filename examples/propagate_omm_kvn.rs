//! Example demonstrating how to propagate an OMM-KVN to a desired datetime using the mako-sgp4 crate.

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------
use mako_sgp4::{DateTime, Timezone, from_omm_kvn_string, sgp4_prop_datetime};

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
    // Define the OMM-KVN string
    let omm_kvn = "\
    CCSDS_OMM_VERS = 2.0
    CREATION_DATE  = 
    ORIGINATOR     = 

    OBJECT_NAME    = 2026-106A
    OBJECT_ID      = 2026-106A
    CENTER_NAME    = EARTH
    REF_FRAME      = TEME
    TIME_SYSTEM    = UTC
    MEAN_ELEMENT_THEORY = SGP/SGP4

    EPOCH          = 2026-06-14T15:07:48.259488
    MEAN_MOTION    = 15.11169557
    ECCENTRICITY   = .00147468
    INCLINATION    = 97.5103
    RA_OF_ASC_NODE = 247.7605
    ARG_OF_PERICENTER = 169.6213
    MEAN_ANOMALY   = 190.5325

    EPHEMERIS_TYPE = 0
    CLASSIFICATION_TYPE = U
    NORAD_CAT_ID   = 69097
    ELEMENT_SET_NO = 999
    REV_AT_EPOCH   = 459
    BSTAR          = .39221734E-3
    MEAN_MOTION_DOT = .6535E-4
    MEAN_MOTION_DDOT = 0
    ";

    // Parse the OMM-KVN string into an SGP4 propagator. Returns a vector of SGP4 propagators
    let omm_kvn_sgp4s = from_omm_kvn_string(omm_kvn).unwrap();

    // Define a desired datetime to propagate to
    let desired_datetime = DateTime {
        year: 2026,
        month: 9,
        day: 5,
        hour: 15,
        minute: 7,
        second: 48.259488,
        timezone: Timezone::UTC,
    };

    // Propagate the OMM-KVN propagator to the epoch time. Returns a StateVector struct
    let state = sgp4_prop_datetime(&omm_kvn_sgp4s[0], &desired_datetime).unwrap();

    // Print the result in TEME coordinates
    println!(
        "{}\nr_TEME = [{:.3}, {:.3}, {:.3}] km\nv_TEME = [{:.3}, {:.3}, {:.3}] km/s",
        omm_kvn_sgp4s[0].gp.common_name,
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
