//! mako-sgp4: parse and propagate general perturbation element sets with SGP4/SDP4.
//!
//! The crate is `#![no_std]` and needs only `alloc`. The `std` feature (on by
//! default) adds the `*_file` read and write functions, and the `csv` feature
//! turns it on.
//!
//! # Examples
//! ```rust
//! use mako_sgp4::{from_tle_lines, sgp4_prop_delta};
//!
//! // Parse a TLE and propagate to epoch
//! let line0 = "ISS (ZARYA)";
//! let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
//! let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
//! let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
//! let state = sgp4_prop_delta(&sgp4, 0.0).unwrap();
//!
//! assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
//! assert!(state.position[0].is_finite());
//! ```

#![no_std]

extern crate alloc;

#[cfg(any(feature = "std", test))]
extern crate std;

pub mod common;
pub mod gp;
// Not yet called outside its own tests; remove this allow once the GP fit
// uses all of it (implementation step 6)
#[cfg_attr(not(test), allow(dead_code))]
mod linalg;
mod math;
pub mod sgp4;
pub mod time;

#[doc(inline)]
pub use common::{CoordinateFrame, StateVector, WGS72, WGS84, Wgs, calc_period, deg2rad};
#[doc(inline)]
pub use gp::{
    GenPerturbElementSet, GpError, from_omm_kvn_lines, from_omm_kvn_string, from_tle_lines,
    from_tle_string, to_omm_kvn_string, to_tle_string,
};
#[cfg(feature = "csv")]
#[doc(inline)]
pub use gp::{from_omm_csv_file, from_omm_csv_string, to_omm_csv_file, to_omm_csv_string};
#[cfg(all(feature = "json", feature = "std"))]
#[doc(inline)]
pub use gp::{from_omm_json_file, to_omm_json_file};
#[cfg(feature = "json")]
#[doc(inline)]
pub use gp::{from_omm_json_string, to_omm_json_string};
#[cfg(feature = "std")]
#[doc(inline)]
pub use gp::{from_omm_kvn_file, from_tle_file, to_omm_kvn_file, to_tle_file};
#[cfg(all(feature = "xml", feature = "std"))]
#[doc(inline)]
pub use gp::{from_omm_xml_file, to_omm_xml_file};
#[cfg(feature = "xml")]
#[doc(inline)]
pub use gp::{from_omm_xml_string, to_omm_xml_string};
#[doc(inline)]
pub use sgp4::{Sgp4, Sgp4Error, calc_theta_g, init_sgp4, sgp4_prop_datetime, sgp4_prop_delta};
#[doc(inline)]
pub use time::{
    DateError, DateTime, Timezone, dayofyr2utc, utc2jday, utc2mjday, validate_datetime,
};
