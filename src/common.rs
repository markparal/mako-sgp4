//! Module to handle common constants, coordinate types, and utility functions shared across
//! GP parsing and SGP4 propagation.

// ------------------
// External Libraries
// ------------------
use std::f64::consts::PI;

// ------------------
// Internal Libraries
// ------------------

// -------
// Structs
// -------

/// World Geodetic System (WGS) parameters
///
/// This struct contains the important Earth parameters defined by different WGS
/// standards (e.g. WGS-72, WGS-84).
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::WGS72;
///
/// // WGS-72 is the default Earth model for TLE / SGP4
/// let wgs = WGS72;
/// assert!((wgs.mu - 398600.8).abs() < 1e-9);
/// assert!((wgs.r_earth_eq - 6378.135).abs() < 1e-9);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
#[derive(Default, Clone, Copy)]
pub struct Wgs {
    /// Earth's standard gravitational parameter \[km^3/s^2\]
    pub mu: f64,

    /// Earth's equatorial radius \[km\]
    pub r_earth_eq: f64,

    /// Earth's J2 harmonic \[\]
    pub j2: f64,

    /// k2 constant \[Earth Radii^2\]
    pub k2: f64,

    /// Earth's J3 harmonic \[\]
    pub j3: f64,

    /// Earth's J4 harmonic \[\]
    pub j4: f64,

    /// k4 constant \[Earth Radii^4\]
    pub k4: f64,

    /// The square root of the standard gravitational parameter \[Earth radii^1.5 / min\]
    pub ke: f64,
}

/// Satellite state vector
///
/// Position \[km\] and velocity \[km/s\] of a satellite in a specified coordinate frame.
/// SGP4 propagation returns this type in [`CoordinateFrame::TEME`].
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::{CoordinateFrame, StateVector};
///
/// // Define a TEME state (position in km, velocity in km/s)
/// let state = StateVector {
///     r_x: 1.0,
///     r_y: 0.0,
///     r_z: 0.0,
///     v_x: 0.0,
///     v_y: 7.5,
///     v_z: 0.0,
///     coordinate_frame: CoordinateFrame::TEME,
/// };
///
/// // Assert the frame used by SGP4
/// assert_eq!(state.coordinate_frame, CoordinateFrame::TEME);
/// ```
///
/// # References
#[derive(Default, Clone, Copy)]
pub struct StateVector {
    /// Position X component \[km\]
    pub r_x: f64,

    /// Position Y component \[km\]
    pub r_y: f64,

    /// Position Z component \[km\]
    pub r_z: f64,

    /// Velocity X component \[km/s\]
    pub v_x: f64,

    /// Velocity Y component \[km/s\]
    pub v_y: f64,

    /// Velocity Z component \[km/s\]
    pub v_z: f64,

    /// Coordinate frame of the state vector
    pub coordinate_frame: CoordinateFrame,
}

// -----
// Enums
// -----

/// Coordinate frames
///
/// Represents the coordinate frame used for a [`StateVector`].
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::CoordinateFrame;
///
/// // SGP4 state vectors are TEME; J2000 is the enum default
/// let frame_teme = CoordinateFrame::TEME;
/// let frame_j2000 = CoordinateFrame::J2000;
///
/// assert_eq!(frame_j2000, CoordinateFrame::default());
/// assert_ne!(frame_teme, frame_j2000);
/// ```
///
/// # References
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateFrame {
    /// J2000, an Earth-centered inertial (ECI) coordinate frame
    #[default]
    J2000,

    /// True Equator Mean Equinox (TEME), an Earth-centered inertial (ECI) coordinate frame
    TEME,
}

// ------
// Traits
// ------

// ---------
// Constants
// ---------

/// Fundamental and derived constants for WGS-72
///
/// Earth model parameters used as the default for TLE/GP processing with SGP4.
///
/// - `mu`: 398600.8 - standard gravitational parameter \[km^3 / s^2\]
/// - `r_earth_eq`: 6378.135 - Earth's equatorial radius \[km\]
/// - `j2`: 0.001082616 - second zonal harmonic (Earth's oblateness)
/// - `k2`: 0.000541308 - `0.5 * j2` \[Earth radii^2\]
/// - `j3`: -0.00000253881 - third zonal harmonic (pear-shaped component)
/// - `j4`: -0.00000165597 - fourth zonal harmonic
/// - `k4`: 0.00000062098875 - `-3/8 * j4` \[Earth radii^4\]
/// - `ke`: 0.07436691613317 - square root of `mu` \[Earth radii^1.5 / min\]
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::WGS72;
///
/// // TLE / SGP4 default Earth model
/// assert!((WGS72.mu - 398600.8).abs() < 1e-9);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub const WGS72: Wgs = Wgs {
    mu: 398600.8,
    r_earth_eq: 6378.135,
    j2: 0.001082616,
    k2: 0.000541308,
    j3: -0.00000253881,
    j4: -0.00000165597,
    k4: 0.00000062098875,
    ke: 0.07436691613317,
};

/// Fundamental and derived constants for WGS-84
///
/// Earth model parameters for the WGS-84 geodetic system.
///
/// - `mu`: 398600.5 - standard gravitational parameter \[km^3 / s^2\]
/// - `r_earth_eq`: 6378.137 - Earth's equatorial radius \[km\]
/// - `j2`: 0.00108262998905 - second zonal harmonic (Earth's oblateness)
/// - `k2`: 0.000541314994525 - `0.5 * j2` \[Earth radii^2\]
/// - `j3`: -0.00000253215306 - third zonal harmonic (pear-shaped component)
/// - `j4`: -0.00000161098761 - fourth zonal harmonic
/// - `k4`: 0.0000006041203538 - `-3/8 * j4` \[Earth radii^4\]
/// - `ke`: 0.07436685316871 - square root of `mu` \[Earth radii^1.5 / min\]
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::{WGS72, WGS84};
///
/// // WGS-84 differs from the TLE default (WGS-72)
/// assert!((WGS84.mu - 398600.5).abs() < 1e-9);
/// assert_ne!(WGS84.mu, WGS72.mu);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub const WGS84: Wgs = Wgs {
    mu: 398600.5,
    r_earth_eq: 6378.137,
    j2: 0.00108262998905,
    k2: 0.000541314994525,
    j3: -0.00000253215306,
    j4: -0.00000161098761,
    k4: 0.0000006041203538,
    ke: 0.07436685316871,
};

// ---------
// Functions
// ---------

/// Convert an angle from degrees to radians
///
/// Multiplies the input angle by `pi / 180`.
///
/// # Arguments
/// * `theta` - The angle in degrees
///
/// # Returns
/// * `theta_rad` - The angle in radians
///
/// # Examples
/// ```rust
/// use std::f64::consts::PI;
/// use mako_sgp4::common::deg2rad;
///
/// // Define a right angle in degrees
/// let theta = 90.0;
///
/// // Convert to radians
/// let theta_rad = deg2rad(theta);
///
/// // Assert 90 deg is pi/2
/// assert!((theta_rad - PI / 2.0).abs() < 1e-12);
/// ```
pub fn deg2rad(theta: f64) -> f64 {
    // Convert to radians
    PI / 180. * theta
}

/// Calculate the orbital period from semi-major axis and gravitational parameter
///
/// Uses Kepler's third law: `T = 2 pi sqrt(a^3 / mu)`, returned in minutes.
///
/// # Arguments
/// * `a` - The semi-major axis \[km\]
/// * `mu` - The standard gravitational parameter \[km^3 / s^2\]
///
/// # Returns
/// * `period` - The period in minutes \[min\]
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::{WGS72, calc_period};
///
/// // Circular orbit at Earth's equatorial radius (WGS-72)
/// let period = calc_period(WGS72.r_earth_eq, WGS72.mu);
///
/// // Period is about 84.5 minutes
/// assert!((period - 84.489).abs() < 1e-3);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
pub fn calc_period(a: f64, mu: f64) -> f64 {
    // Calculate time period in minutes
    2. * PI * (a.powi(3) / mu).sqrt() / 60.
}

// ----------
// Unit Tests
// ----------
