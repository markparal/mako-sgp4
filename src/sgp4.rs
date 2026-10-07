//! Module for propagating GP element sets with SGP4

// ------------------
// External Libraries
// ------------------
use std::f64::consts::PI;

// ------------------
// Internal Libraries
// ------------------
use crate::common::{CoordinateFrame, StateVector, WGS72, Wgs, calc_period, deg2rad};
use crate::gp::GenPerturbElementSet;
use crate::time::{DateError, DateTime, utc2jday};

// -------
// Structs
// -------

/// Simplified General Perturbations 4 (SGP4) parameters
///
/// This struct contains the epoch parameters which are necessary to propagate the state vectors of a satellite with SGP4.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a TLE; from_tle_lines initializes SGP4 parameters
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // ISS is a near-Earth satellite
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// assert!(!sgp4.deep_space);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone)]
pub struct Sgp4 {
    /// WGS model
    pub wgs: Wgs,

    /// General Perturbation Element Set
    pub gp: GenPerturbElementSet,

    /// Julian date at epoch \[days\]
    pub jd0: f64,

    /// Fractional Julian date at epoch \[days\]
    pub jdfrac0: f64,

    /// Deep space satellite
    pub deep_space: bool,

    /// Brouwer mean elements at epoch
    pub brouwer0: BrouwerMeanElements,

    /// Atmospheric Drag Parameters
    pub atm_params: AtmDragParams,

    /// Earth Zonal Harmonics Parameters
    pub zonal_params: EarthZonalParams,

    /// Solar 3rd Body Parameters
    pub solar_params: ThirdBodyParams,

    /// Lunar 3rd Body Parameters
    pub lunar_params: ThirdBodyParams,

    /// Account for whole day resonance effects of Earth's gravity
    pub whole_day_resonance: bool,

    /// Whole day resonance parameters of Earth's gravity
    pub whole_day_resonance_params: WholeDayResonanceParams,

    /// Account for half day resonance effects of Earth's gravity
    pub half_day_resonance: bool,

    /// Half day resonance parameters of Earth's gravity
    pub half_day_resonance_params: HalfDayResonanceParams,
}

/// Brouwer Mean Orbital Elements
///
/// This struct contains the mean orbital elements of a TLE converted to Brouwer convention. TLEs report mean orbital elements
/// in Kozai convention.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a TLE to recover Brouwer mean elements
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // Mean motion and semi-major axis are physical
/// assert!(sgp4.brouwer0.n > 0.0);
/// assert!(sgp4.brouwer0.a > 1.0);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct BrouwerMeanElements {
    /// Orbital inclination \[rad\]
    pub i: f64,

    /// The cosine of the orbital inclination
    pub theta: f64,

    /// Right ascension of the ascending node (RAAN) \[rad\]
    pub raan: f64,

    /// Orbital eccentricity \[\]
    pub e: f64,

    /// The square root of 1 minus the orbital eccentricity squared \[\]
    pub beta: f64,

    /// Argument of perigee \[rad\]
    pub omega: f64,

    /// Mean anomaly \[rad\]
    pub m: f64,

    /// Mean motion \[rad/min\]
    pub n: f64,

    /// Semi-major axis \[Earth Radii\]
    pub a: f64,

    /// The orbital period \[mins\]
    pub period: f64,
}

/// Atmospheric Drag Effects
///
/// This struct contains the parameters necessary to account for the impacts of atmospheric drag on an orbit.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a TLE to recover atmospheric drag parameters
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // Perigee height is above the Earth
/// assert!(sgp4.atm_params.hp > 0.0);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct AtmDragParams {
    /// Perigee height \[km\]
    pub hp: f64,

    /// q0 parameter of power-law density function \[Earth Radii\]
    pub q0: f64,

    /// s parameter of power-law density function \[Earth Radii\]
    pub s: f64,

    /// Zeta constant \[1 / Earth Radii\]
    pub zeta: f64,

    /// Eta constant \[\]
    pub eta: f64,

    /// C1 constant \[\]
    pub c1: f64,

    /// C3 constant \[\]
    pub c3: f64,

    /// C4 constant \[\]
    pub c4: f64,

    /// C5 constant \[\]
    pub c5: f64,

    /// D2 constant \[\]
    pub d2: f64,

    /// D3 constant \[\]
    pub d3: f64,

    /// D4 constant \[\]
    pub d4: f64,
}

/// Earth Zonal Harmonics
///
/// This struct contains the parameters necessary to account for the impacts of Earth's zonal harmonics on an orbit.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a TLE to recover Earth zonal harmonic rates
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // RAAN precession is non-zero for an inclined LEO
/// assert!(sgp4.zonal_params.raan_dot.is_finite());
/// assert!(sgp4.zonal_params.raan_dot != 0.0);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct EarthZonalParams {
    /// Rate of change of mean anomaly \[rad / min\]
    pub m_dot: f64,

    /// Rate of change of the argument of perigee \[rad / min\]
    pub omega_dot: f64,

    /// Rate of change of the right ascension of the ascending node \[rad / min\]
    pub raan_dot: f64,
}

/// Solar and Lunar 3rd Body Effects
///
/// This struct contains the parameters necessary to account for the impacts of the Sun and Moon on an orbit.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a near-Earth TLE; third-body rates stay at the default zeros
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // Near-Earth satellites do not apply lunar/solar secular rates
/// assert_eq!(sgp4.solar_params.n, 0.0);
/// assert_eq!(sgp4.lunar_params.n, 0.0);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct ThirdBodyParams {
    /// Third body orbital inclination cosine \[\]
    pub cos_i: f64,

    /// Third body orbital inclination sine \[\]
    pub sin_i: f64,

    /// Third body eccentricity \[\]
    pub e: f64,

    /// Third body mean motion \[rad/min\]
    pub n: f64,

    /// Third body argument of perigee cosine \[\]
    pub cos_omega: f64,

    /// Third body argument of perigee sine \[\]
    pub sin_omega: f64,

    /// Third body right ascension of the ascending node (RAAN) \[rad\]
    pub raan: f64,

    /// Third body mean anomaly \[rad\]
    pub m: f64,

    /// The square root of 1 minus the orbital eccentricity squared \[\]
    pub beta: f64,

    /// Third body perturbation coefficient \[rad/min\]
    pub c: f64,

    /// x1 constant
    pub x1: f64,

    /// x2 constant
    pub x2: f64,

    /// x3 constant
    pub x3: f64,

    /// x4 constant
    pub x4: f64,

    /// x5 constant
    pub x5: f64,

    /// x6 constant
    pub x6: f64,

    /// x7 constant
    pub x7: f64,

    /// x8 constant
    pub x8: f64,

    /// z1 constant
    pub z1: f64,

    /// z2 constant
    pub z2: f64,

    /// z3 constant
    pub z3: f64,

    /// z11 constant
    pub z11: f64,

    /// z13 constant
    pub z13: f64,

    /// z21 constant
    pub z21: f64,

    /// z23 constant
    pub z23: f64,

    /// z22 constant
    pub z22: f64,

    /// z12 constant
    pub z12: f64,

    /// z31 constant
    pub z31: f64,

    /// z32 constant
    pub z32: f64,

    /// z33 constant
    pub z33: f64,

    /// Rate of change of the orbital eccentricity \[1 / min\]
    pub e_dot: f64,

    /// Rate of change of the orbital inclination \[rad / min\]
    pub i_dot: f64,

    /// Rate of change of the mean anomaly \[rad / min\]
    pub m_dot: f64,

    /// Rate of change of the argument of perigee \[rad / min\]
    pub omega_dot: f64,

    /// Rate of change of the right ascension of the ascending node \[rad / min\]
    pub raan_dot: f64,
}

/// Half day resonance effects of Earth's gravity
///
/// This struct contains the parameters necessary to account for the impacts of half day resonance effects on an orbit.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a near-Earth TLE; 12-hour resonance is not applied
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// assert!(!sgp4.half_day_resonance);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct HalfDayResonanceParams {
    /// Greenwich sidereal time at epoch \[rad\]
    pub theta_g: f64,

    /// lam0 constant
    pub lam0: f64,

    /// lam0 rate of change
    pub lam0_dot: f64,

    /// d2201 constant
    pub d2201: f64,

    /// d2211 constant
    pub d2211: f64,

    /// d3210 constant
    pub d3210: f64,

    /// d3222 constant
    pub d3222: f64,

    /// d5220 constant
    pub d5220: f64,

    /// d5232 constant
    pub d5232: f64,

    /// d4422 constant
    pub d4422: f64,

    /// d5421 constant
    pub d5421: f64,

    /// d5433 constant
    pub d5433: f64,

    /// d4410 constant
    pub d4410: f64,
}

/// Whole day resonance effects of Earth's gravity
///
/// This struct contains the parameters necessary to account for the impacts of whole day resonance effects on an orbit.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a near-Earth TLE; 24-hour resonance is not applied
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// assert!(!sgp4.whole_day_resonance);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
#[derive(Default, Clone, Copy)]
pub struct WholeDayResonanceParams {
    /// Greenwich sidereal time at epoch \[rad\]
    pub theta_g: f64,

    /// lam0 constant
    pub lam0: f64,

    /// lam0 rate of change
    pub lam0_dot: f64,

    /// lam31 constant
    pub lam31: f64,

    /// lam22 constant
    pub lam22: f64,

    /// lam33 constant
    pub lam33: f64,

    /// delta1 constant
    pub delta1: f64,

    /// delta2 constant
    pub delta2: f64,

    /// delta3 constant
    pub delta3: f64,
}

// -----
// Enums
// -----

/// SGP4 errors
///
/// Failures that can occur while initializing or propagating an SGP4 model.
/// The element-set variants match Vallado's non-physical orbit checks. The
/// date variant wraps [`DateError`] when Julian-day conversion fails.
///
/// # Examples
/// ```rust
/// use mako_sgp4::sgp4::Sgp4Error;
/// use mako_sgp4::time::DateError;
///
/// // A non-UTC datetime is a recoverable SGP4 error
/// let err = Sgp4Error::InvalidDateTime(DateError::DateNotUTC);
/// assert_eq!(err, Sgp4Error::InvalidDateTime(DateError::DateNotUTC));
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
#[derive(Debug, Clone, PartialEq)]
pub enum Sgp4Error {
    /// Mean motion is less than or equal to zero
    InvalidMeanMotion,

    /// Mean eccentricity is outside the range 0.0 to 1.0
    InvalidMeanEccentricity,

    /// Perturbed eccentricity is outside the range 0.0 to 1.0
    InvalidPerturbedEccentricity,

    /// Semilatus rectum is less than zero
    InvalidSemilatusRectum,

    /// Satellite has decayed (position radius is less than 1 Earth radius)
    SatelliteDecayed,

    /// Epoch or propagation datetime could not be converted to Julian date
    InvalidDateTime(DateError),

    /// Propagated position or velocity is not finite (NaN or infinite)
    NonFiniteState,
}

// ------
// Traits
// ------

// ---------
// Constants
// ---------

/// Conversion factor from revolutions per day to radians per minute
///
/// Used to convert GP mean motion (`revs/day`) into the SGP4 internal
/// mean-motion unit of `rad/min`.
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
const XPDOTP: f64 = 1440.0 / (2.0 * PI);

/// Earth's rotational rate \[rad/min\]
///
/// Sidereal Earth rotation rate used when evaluating Greenwich sidereal time
/// and half-/whole-day resonance terms during deep-space propagation.
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
const RPTIM: f64 = 4.375_269_088_011_3e-3;

// ---------
// Functions
// ---------

/// Build an [`Sgp4`] struct for state propagation from a [`GenPerturbElementSet`] struct
///
/// Given a [`GenPerturbElementSet`] struct, calculate the time-independent parameters necessary
/// to propagate a satellite's states in time. These parameters include
/// - Brouwer mean orbital elements
/// - Atmospheric drag parameters
/// - Earth zonal harmonics parameters
/// - Solar and Lunar 3rd body effects
/// - Resonance effects of Earth's gravity
///
/// # Arguments
/// * `gp` - The General Perturbation Element Set parameters
/// * `wgs` - Optional, specify World Geodetic System (WGS) parameters (defaults to WGS-72, the standard for TLEs)
///
/// # Returns
/// * `Ok(Sgp4)` - The time-independent parameters for the SGP4 propagator
/// * `Err(Sgp4Error)` - If Julian-day conversion or epoch propagation fails
///
/// # Errors
/// * [`Sgp4Error::InvalidDateTime`] - If the GP epoch datetime is not UTC or Julian-day conversion fails
/// * [`Sgp4Error::InvalidMeanMotion`] - If mean motion is less than or equal to zero
/// * [`Sgp4Error::InvalidMeanEccentricity`] - If mean eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidPerturbedEccentricity`] - If perturbed eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidSemilatusRectum`] - If the semilatus rectum is less than zero
/// * [`Sgp4Error::SatelliteDecayed`] - If the satellite has decayed
/// * [`Sgp4Error::NonFiniteState`] - If the propagated state is NaN or infinite
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::WGS72;
/// use mako_sgp4::gp::from_tle_lines;
/// use mako_sgp4::sgp4::init_sgp4;
///
/// // Parse a TLE to obtain a General Perturbation (GP) Element Set
/// let tle_line0 = "ISS (ZARYA)";
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let parsed = from_tle_lines(tle_line1, tle_line2, Some(tle_line0)).unwrap();
///
/// // Initialize the SGP4 propagator (None uses WGS-72)
/// let sgp4 = init_sgp4(&parsed.gp, Some(&WGS72)).unwrap();
/// assert!(!sgp4.deep_space);
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn init_sgp4(gp: &GenPerturbElementSet, wgs: Option<&Wgs>) -> Result<Sgp4, Sgp4Error> {
    // Use WGS72 or custom WGS models if provided
    let wgs_sgp4 = if let Some(wgs_passed) = wgs {
        *wgs_passed
    } else {
        WGS72
    };

    // Extract General Perturbation (GP) Element Set contents in proper units
    let i0 = deg2rad(gp.inclination); // [rad]
    let n0_kozai = gp.mean_motion / XPDOTP; // [rad/min]
    let e0 = gp.eccentricity; // []
    let omega0 = deg2rad(gp.argument_of_perigee); // [rad]
    let raan0 = deg2rad(gp.right_ascension_of_ascending_node); // [rad]
    let m0 = deg2rad(gp.mean_anomaly); // [rad]

    // Extract GP epoch in Julian day format
    let (jd0, jdfrac0) = utc2jday(&gp.epoch_datetime).map_err(Sgp4Error::InvalidDateTime)?;

    // Recover Brouwer mean motion from Kozai mean motion (mean motion in GP)
    let theta0 = i0.cos();
    let beta0 = (1. - e0.powi(2)).sqrt();
    let a1 = (wgs_sgp4.ke / n0_kozai).powf(2. / 3.);
    let delta1 = (3. / 2.) * (wgs_sgp4.k2 / a1.powf(2.)) * (3. * i0.cos().powf(2.) - 1.)
        / (1. - e0.powf(2.)).powf(3. / 2.);
    let a2 = a1 * (1. - (1. / 3.) * delta1 - delta1.powf(2.) - (134. / 81.) * delta1.powf(3.));
    let delta0 = (3. / 2.) * (wgs_sgp4.k2 / a2.powf(2.)) * (3. * i0.cos().powf(2.) - 1.)
        / (1. - e0.powf(2.)).powf(3. / 2.);
    let n0 = n0_kozai / (1. + delta0); // [rad/min]
    let a0 = (wgs_sgp4.ke / n0).powf(2. / 3.); // [Earth radii]
    let a0_km = a0 * wgs_sgp4.r_earth_eq; // [km]
    let period0 = calc_period(a0_km, wgs_sgp4.mu); // [min]

    // Store Brouwer mean elements
    let brouwer0 = BrouwerMeanElements {
        i: i0,
        theta: theta0,
        raan: raan0,
        e: e0,
        beta: beta0,
        omega: omega0,
        m: m0,
        n: n0,
        a: a0,
        period: period0,
    };

    // Initialize atmospheric drag parameters
    let atm_params = init_atm_effects(&wgs_sgp4, gp, &brouwer0);

    // Initialize Earth zonal harmonics parameters
    let zonal_params = init_zonal_effects(&wgs_sgp4, &brouwer0);

    // Check for deep space satellite
    let mut deep_space = false;
    if period0 >= 225. {
        deep_space = true;
    }

    // Lunar and solar gravity effects
    let (lunar_params, solar_params) =
        init_lunar_solar_effects(deep_space, jd0, jdfrac0, &brouwer0);

    // Earth gravity resonance effects (use Vallado criteria instead of Hoots)
    let mut whole_day_resonance = false;
    let mut half_day_resonance = false;
    let mut whole_day_resonance_params = WholeDayResonanceParams::default();
    let mut half_day_resonance_params = HalfDayResonanceParams::default();
    if (n0 > 0.0034906585) && (n0 < 0.0052359877) {
        whole_day_resonance = true;
        whole_day_resonance_params = init_earth_gravity_resonance_wholeday(
            jd0,
            jdfrac0,
            &brouwer0,
            &zonal_params,
            &lunar_params,
            &solar_params,
        );
    }
    if (8.26e-3..=9.24e-3).contains(&n0) && (e0 >= 0.5) {
        half_day_resonance = true;
        half_day_resonance_params = init_earth_gravity_resonance_halfday(
            jd0,
            jdfrac0,
            &brouwer0,
            &zonal_params,
            &lunar_params,
            &solar_params,
        );
    }

    // Construct SGP4 propagator
    let sgp4 = Sgp4 {
        wgs: wgs_sgp4,
        gp: gp.clone(),
        jd0,
        jdfrac0,
        deep_space,
        brouwer0,
        atm_params,
        zonal_params,
        lunar_params,
        solar_params,
        whole_day_resonance,
        whole_day_resonance_params,
        half_day_resonance,
        half_day_resonance_params,
    };

    // Propagate to epoch so initialization failures surface through the same checks as propagation
    sgp4_prop_delta(&sgp4, 0.0)?;

    Ok(sgp4)
}

/// Initialize the atmospheric drag effects
///
/// Computes the power-law density constants (`q0`, `s`, `zeta`, `eta`) and the
/// drag coefficients (`C1`-`C5`, `D2`-`D4`) used during SGP4 secular updates.
///
/// # Arguments
/// * `wgs` - The WGS model
/// * `gp` - The General Perturbation (GP) Element Set
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `AtmDragParams` - The atmospheric drag parameters
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn init_atm_effects(
    wgs: &Wgs,
    gp: &GenPerturbElementSet,
    brouwer0: &BrouwerMeanElements,
) -> AtmDragParams {
    // Define initial constants
    let a30 = -wgs.j3; // [Earth Radii^3]
    let q0 = (120. + wgs.r_earth_eq) / wgs.r_earth_eq; // [Earth radii]

    // Determine parameter s based on perigee height
    let rp = brouwer0.a * (1. - brouwer0.e); // Radius of perigee [Earth Radii]
    let hp = (rp - 1.) * wgs.r_earth_eq; // Perigee height [km]

    let s: f64 = if hp >= 156. {
        (78. + wgs.r_earth_eq) / wgs.r_earth_eq
    } else if hp >= 98. {
        (hp - 78. + wgs.r_earth_eq) / wgs.r_earth_eq
    } else {
        (20. + wgs.r_earth_eq) / wgs.r_earth_eq
    }; // [Earth radii]

    // Calculate atmospheric drag parameters
    let zeta = 1. / (brouwer0.a - s);
    let eta = brouwer0.a * brouwer0.e * zeta;
    let psisq = (1. - eta.powi(2)).abs(); // Vallado uses abs to handle the case when eta > 1 (sub-orbital / decayed orbits)

    let c2_1 = (q0 - s).powi(4) * zeta.powi(4) * brouwer0.n * psisq.powf(-3.5);
    let c2_2 = brouwer0.a
        * (1. + (3. / 2.) * eta.powi(2) + 4. * brouwer0.e * eta + brouwer0.e * eta.powi(3));
    let c2_3 = (3. / 2.)
        * (wgs.k2 * zeta / psisq)
        * (-(1. / 2.) + (3. / 2.) * brouwer0.theta.powi(2))
        * (8. + 24. * eta.powi(2) + 3. * eta.powi(4));
    let c2 = c2_1 * (c2_2 + c2_3);

    let c1 = gp.bstar * c2;
    // Vallado drop C3 when eccentricity is too small (avoids /e blow-up)
    let c3 = if brouwer0.e > 1.0e-4 {
        ((q0 - s).powf(4.) * zeta.powf(5.) * a30 * brouwer0.n * brouwer0.i.sin())
            / (wgs.k2 * brouwer0.e)
    } else {
        0.0
    };

    let c4_1 = 2.
        * brouwer0.n
        * (q0 - s).powi(4)
        * zeta.powi(4)
        * brouwer0.a
        * brouwer0.beta.powi(2)
        * psisq.powf(-3.5);
    let c4_2 = 2. * eta * (1. + brouwer0.e * eta) + 0.5 * brouwer0.e + 0.5 * eta.powi(3);
    let c4_3 = 2. * wgs.k2 * zeta / (brouwer0.a * psisq);
    let c4_4 = 3.
        * (1. - 3. * brouwer0.theta.powi(2))
        * (1. + 3. / 2. * eta.powi(2) - 2. * brouwer0.e * eta - 0.5 * brouwer0.e * eta.powi(3));
    let c4_5 = 3. / 4.
        * (1. - brouwer0.theta.powi(2))
        * (2. * eta.powi(2) - brouwer0.e * eta - brouwer0.e * eta.powi(3))
        * (2. * brouwer0.omega).cos();
    let c4 = c4_1 * (c4_2 - c4_3 * (c4_4 + c4_5));

    let c5_1 = 2.
        * (q0 - s).powi(4)
        * zeta.powi(4)
        * brouwer0.a
        * brouwer0.beta.powi(2)
        * psisq.powf(-3.5);
    let c5_2 = 1. + 11. / 4. * eta * (eta + brouwer0.e) + brouwer0.e * eta.powi(3);
    let c5 = c5_1 * c5_2;

    let d2 = 4. * brouwer0.a * zeta * c1.powi(2);
    let d3 = 4. / 3. * brouwer0.a * zeta.powi(2) * (17. * brouwer0.a + s) * c1.powi(3);
    let d4 =
        2. / 3. * brouwer0.a.powi(2) * zeta.powi(3) * (221. * brouwer0.a + 31. * s) * c1.powi(4);

    // Store atmospheric drag parameters
    AtmDragParams {
        hp,
        q0,
        s,
        zeta,
        eta,
        c1,
        c3,
        c4,
        c5,
        d2,
        d3,
        d4,
    }
}

/// Initialize the Earth zonal harmonics effects
///
/// Computes the secular rates of mean anomaly, argument of perigee, and RAAN
/// due to Earth's J2 / J4 zonal harmonics at the GP epoch.
///
/// # Arguments
/// * `wgs` - The WGS model
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `EarthZonalParams` - The Earth zonal harmonics parameters
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn init_zonal_effects(wgs: &Wgs, brouwer0: &BrouwerMeanElements) -> EarthZonalParams {
    // Calculate orbital element rates of change due to zonal harmonics
    let m_dot_1 = 3. * wgs.k2 * (-1. + 3. * brouwer0.theta.powi(2))
        / (2. * brouwer0.a.powi(2) * brouwer0.beta.powi(3));
    let m_dot_2 =
        3. * wgs.k2.powi(2) * (13. - 78. * brouwer0.theta.powi(2) + 137. * brouwer0.theta.powi(4))
            / (16. * brouwer0.a.powi(4) * brouwer0.beta.powi(7));
    let m_dot = (m_dot_1 + m_dot_2) * brouwer0.n;

    let omega_dot_1 = -3. * wgs.k2 * (1. - 5. * brouwer0.theta.powi(2))
        / (2. * brouwer0.a.powi(2) * brouwer0.beta.powi(4));
    let omega_dot_2 =
        3. * wgs.k2.powi(2) * (7. - 114. * brouwer0.theta.powi(2) + 395. * brouwer0.theta.powi(4))
            / (16. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let omega_dot_3 =
        5. * wgs.k4 * (3. - 36. * brouwer0.theta.powi(2) + 49. * brouwer0.theta.powi(4))
            / (4. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let omega_dot = (omega_dot_1 + omega_dot_2 + omega_dot_3) * brouwer0.n;

    let raan_dot_1 = -3. * wgs.k2 * brouwer0.theta / (brouwer0.a.powi(2) * brouwer0.beta.powi(4));
    let raan_dot_2 = 3. * wgs.k2.powi(2) * (4. * brouwer0.theta - 19. * brouwer0.theta.powi(3))
        / (2. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let raan_dot_3 = 5. * wgs.k4 * brouwer0.theta * (3. - 7. * brouwer0.theta.powi(2))
        / (2. * brouwer0.a.powi(4) * brouwer0.beta.powi(8));
    let raan_dot = (raan_dot_1 + raan_dot_2 + raan_dot_3) * brouwer0.n;

    // Store Earth zonal parameters
    EarthZonalParams {
        m_dot,
        omega_dot,
        raan_dot,
    }
}

/// Initialize the Lunar and Solar third body effects
///
/// For deep-space satellites (period >= 225 min), evaluates Sun/Moon geometry at
/// the TLE epoch and returns secular third-body rates. Near-Earth satellites
/// receive default (zero) parameters.
///
/// # Arguments
/// * `deep_space` - Is this a deep space satellite
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
///
/// # Returns
/// * `(ThirdBodyParams, ThirdBodyParams)` - The Lunar and Solar third body parameters
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn init_lunar_solar_effects(
    deep_space: bool,
    jd0: f64,
    jdfrac0: f64,
    brouwer0: &BrouwerMeanElements,
) -> (ThirdBodyParams, ThirdBodyParams) {
    // Check if the satellite is not in deep space
    if !deep_space {
        return (ThirdBodyParams::default(), ThirdBodyParams::default());
    }

    // Lunar/Solar element epochs (12/31/1899 12:00:00 UTC) \[Julian date\]
    let epoch_sm = 2415020.0;

    // Lunar constants
    // Lunar eccentricity
    let e_m = 0.05490;

    // Lunar mean motion \[rad/min\] (Spacetrack/Vallado digits; extra precision drifts periodics)
    let n_m = 1.5835218e-4;

    // Lunar perturbation coefficient \[rad/min\]
    let c_m = 4.7968065e-7;

    // Lunar right ascension of the ascending node (RAAN) with respect to the ecliptic plane at epoch \[rad\]
    let raan_me0 = 4.5236020;

    // Lunar right ascension of the ascending node (RAAN) with respect to the ecliptic plane time rate of change at epoch \[rad/day\]
    let raan_me0_dot = -9.2422029e-4;

    // Lunar longitude of perigee with respect to the ecliptic plane at epoch \[rad\]
    let u_me0 = 5.8351514;

    // Lunar longitude of perigee with respect to the ecliptic plane time rate of change at epoch \[rad/day\]
    let u_me0_dot = 0.0019443680;

    // Lunar mean anomaly at epoch \[rad\]
    let m_m0 = 4.7199672;

    // Lunar mean anomaly time rate of change at epoch \[rad/day\]
    let m_m0_dot = 0.22997150;

    // Solar constants
    // Solar inclination sin and cos
    let sin_i_s = 0.39785416;
    let cos_i_s = 0.91744867;

    // Solar eccentricity
    let e_s = 0.01675;

    // Solar mean motion \[rad/min\]
    let n_s = 1.19459e-5;

    // Solar right ascension of the ascending node (RAAN) \[rad\]
    let raan_s = 0.0;

    // Solar argument of periapsis cos and sin
    let sin_omega_s = -0.98088458;
    let cos_omega_s = 0.1945905;

    // Solar perturbation coefficient \[rad/min\]
    let c_s = 2.9864797e-6;

    // Solar mean anomaly at epoch \[rad\]
    let m_s0 = 6.2565837;

    // Solar mean anomaly time rate of change at epoch \[rad/day\]
    let m_s0_dot = 0.017201977;

    // Find the difference in time between the Solar / Lunar epoch and the TLE epoch
    let delta_t = jd0 + jdfrac0 - epoch_sm;

    // Calculate the Lunar RAAN wrt to the ecliptic plane at TLE epoch
    let raan_me = (raan_me0 + raan_me0_dot * delta_t).rem_euclid(2.0 * PI);
    let sin_raan_me = raan_me.sin();
    let cos_raan_me = raan_me.cos();

    // Lunar equatorial inclination (Spacetrack/Hoots legacy linearization of the
    // ecliptic-to-equatorial transform - not the exact acos form)
    let cos_i_m = 0.91375164 - 0.03568096 * cos_raan_me;
    let sin_i_m = (1.0 - cos_i_m * cos_i_m).sqrt();

    // Calculate the Lunar longitude of perigee referred to the ecliptic
    let gamma_m = u_me0 + u_me0_dot * delta_t;

    // Lunar ascending-node trig in the equatorial frame (legacy SDP4)
    let sin_raan_m = 0.089683511 * sin_raan_me / sin_i_m;
    let cos_raan_m = (1.0 - sin_raan_m * sin_raan_m).sqrt();
    let raan_m = sin_raan_m.atan2(cos_raan_m);

    // Lunar argument of periapsis (legacy SDP4 zx formulation)
    let mut zx = 0.39785416 * sin_raan_me / sin_i_m;
    let zy = cos_raan_m * cos_raan_me + 0.91744867 * sin_raan_m * sin_raan_me;
    zx = zx.atan2(zy);
    let omega_m = gamma_m + zx - raan_me;
    let sin_omega_m = omega_m.sin();
    let cos_omega_m = omega_m.cos();

    // Calculate the Lunar mean anomaly \[rad\]
    let m_m = (m_m0 + m_m0_dot * delta_t - gamma_m).rem_euclid(2.0 * PI);

    // Calculate the Solar mean anomaly \[rad\]
    let m_s = (m_s0 + m_s0_dot * delta_t).rem_euclid(2.0 * PI);

    // Calculate the Lunar secular rates
    let lunar_params = calc_lunar_solar_secular_rates(
        &ThirdBodyParams {
            cos_i: cos_i_m,
            sin_i: sin_i_m,
            e: e_m,
            n: n_m,
            cos_omega: cos_omega_m,
            sin_omega: sin_omega_m,
            raan: raan_m,
            m: m_m,
            c: c_m,
            ..Default::default()
        },
        brouwer0,
    );

    // Calculate the Solar secular rates
    let solar_params = calc_lunar_solar_secular_rates(
        &ThirdBodyParams {
            cos_i: cos_i_s,
            sin_i: sin_i_s,
            e: e_s,
            n: n_s,
            cos_omega: cos_omega_s,
            sin_omega: sin_omega_s,
            raan: raan_s,
            m: m_s,
            c: c_s,
            ..Default::default()
        },
        brouwer0,
    );

    (lunar_params, solar_params)
}

/// Calculate the secular rates of a third body's orbital elements
///
/// Builds the frozen geometric coefficients (`x*`, `z*`) and secular element
/// rates for one third body (Sun or Moon) relative to the satellite orbit.
/// The input [`ThirdBodyParams`] should contain the third-body geometry
/// (`cos_i`, `sin_i`, `e`, `n`, `cos_omega`, `sin_omega`, `raan`, `m`, `c`);
/// remaining fields are filled in by this function.
///
/// # Arguments
/// * `body` - Third-body geometry at the satellite epoch
/// * `brouwer0` - Satellite Brouwer mean elements at epoch
///
/// # Returns
/// * `ThirdBodyParams` - Secular rates and frozen geometric coefficients (`x*`, `z*`)
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn calc_lunar_solar_secular_rates(
    body: &ThirdBodyParams,
    brouwer0: &BrouwerMeanElements,
) -> ThirdBodyParams {
    // Precompute common quantities
    let cos_raan_diff = (brouwer0.raan - body.raan).cos();
    let sin_raan_diff = (brouwer0.raan - body.raan).sin();
    let cos_omega0 = brouwer0.omega.cos();
    let sin_omega0 = brouwer0.omega.sin();
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();
    let beta_x = (1. - body.e.powi(2)).sqrt();

    // Calculate 3rd body constants
    let a1 = body.cos_omega * cos_raan_diff + body.sin_omega * body.cos_i * sin_raan_diff;
    let a3 = -body.sin_omega * cos_raan_diff + body.cos_omega * body.cos_i * sin_raan_diff;
    let a7 = -body.cos_omega * sin_raan_diff + body.sin_omega * body.cos_i * cos_raan_diff;
    let a8 = body.sin_omega * body.sin_i;
    let a9 = body.sin_omega * sin_raan_diff + body.cos_omega * body.cos_i * cos_raan_diff;
    let a10 = body.cos_omega * body.sin_i;
    let a2 = a7 * cos_i0 + a8 * sin_i0;
    let a4 = a9 * cos_i0 + a10 * sin_i0;
    let a5 = -a7 * sin_i0 + a8 * cos_i0;
    let a6 = -a9 * sin_i0 + a10 * cos_i0;

    let x1 = a1 * cos_omega0 + a2 * sin_omega0;
    let x2 = a3 * cos_omega0 + a4 * sin_omega0;
    let x3 = -a1 * sin_omega0 + a2 * cos_omega0;
    let x4 = -a3 * sin_omega0 + a4 * cos_omega0;
    let x5 = a5 * sin_omega0;
    let x6 = a6 * sin_omega0;
    let x7 = a5 * cos_omega0;
    let x8 = a6 * cos_omega0;

    let z31 = 12. * x1.powi(2) - 3. * x3.powi(2);
    let z32 = 24. * x1 * x2 - 6. * x3 * x4;
    let z33 = 12. * x2.powi(2) - 3. * x4.powi(2);
    let z1 = 6. * (a1.powi(2) + a2.powi(2)) + (1. + brouwer0.e.powi(2)) * z31;
    let z2 = 12. * (a1 * a3 + a2 * a4) + (1. + brouwer0.e.powi(2)) * z32;
    let z3 = 6. * (a3.powi(2) + a4.powi(2)) + (1. + brouwer0.e.powi(2)) * z33;
    let z11 = -6. * a1 * a5 + brouwer0.e.powi(2) * (-24. * x1 * x7 - 6. * x3 * x5);
    let z13 = -6. * a3 * a6 + brouwer0.e.powi(2) * (-24. * x2 * x8 - 6. * x4 * x6);
    let z21 = 6. * a2 * a5 + brouwer0.e.powi(2) * (24. * x1 * x5 - 6. * x3 * x7);
    let z23 = 6. * a4 * a6 + brouwer0.e.powi(2) * (24. * x2 * x6 - 6. * x4 * x8);
    let z22 = 6. * a4 * a5
        + 6. * a2 * a6
        + brouwer0.e.powi(2) * (24. * x2 * x5 + 24. * x1 * x6 - 6. * x4 * x7 - 6. * x3 * x8);
    let z12 = -6. * a1 * a6
        - 6. * a3 * a5
        - brouwer0.e.powi(2) * (24. * x2 * x7 + 24. * x1 * x8 + 6. * x3 * x6 + 6. * x4 * x5);

    // Calculate secular rates
    let e_x_dot =
        -15. * body.c * body.n * (brouwer0.e * brouwer0.beta / brouwer0.n) * (x1 * x3 + x2 * x4);

    let i_x_dot = (-body.c * body.n / (2. * brouwer0.n * brouwer0.beta)) * (z11 + z13);

    let m_x_dot = (-body.c * body.n / brouwer0.n) * (z1 + z3 - 14. - 6. * brouwer0.e.powi(2));

    // Neglect the RAAN rate within 3 deg of 0 or 180 deg inclination to avoid dividing by sin(i) (Vallado sgp4fix for 180 deg incl)
    let raan_rate_valid = brouwer0.i >= deg2rad(3.) && brouwer0.i <= PI - deg2rad(3.);

    let mut raan_x_dot = 0.;
    if raan_rate_valid {
        raan_x_dot = body.c * body.n / (2. * brouwer0.n * brouwer0.beta * sin_i0) * (z21 + z23);
    }

    let mut omega_x_dot = body.c * body.n * brouwer0.beta / brouwer0.n * (z31 + z33 - 6.);
    if raan_rate_valid {
        omega_x_dot -= raan_x_dot * cos_i0;
    }

    // Store the 3rd body parameters
    ThirdBodyParams {
        beta: beta_x,
        x1,
        x2,
        x3,
        x4,
        x5,
        x6,
        x7,
        x8,
        z1,
        z2,
        z3,
        z11,
        z13,
        z21,
        z23,
        z22,
        z12,
        z31,
        z32,
        z33,
        e_dot: e_x_dot,
        i_dot: i_x_dot,
        m_dot: m_x_dot,
        raan_dot: raan_x_dot,
        omega_dot: omega_x_dot,
        ..*body
    }
}

/// Initialize the half day resonance effects of Earth's gravity
///
/// Initializes 12-hour resonance coefficients (`d2201` ... `d5433`) and the
/// auxiliary longitude `lambda0` used by the Euler-Maclaurin deep-space integrator.
///
/// # Arguments
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
/// * `zonal_params` - The Earth zonal harmonics parameters
/// * `lunar_params` - The Lunar third body parameters
/// * `solar_params` - The Solar third body parameters
///
/// # Returns
/// * `HalfDayResonanceParams` - The half day resonance parameters
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn init_earth_gravity_resonance_halfday(
    jd0: f64,
    jdfrac0: f64,
    brouwer0: &BrouwerMeanElements,
    zonal_params: &EarthZonalParams,
    lunar_params: &ThirdBodyParams,
    solar_params: &ThirdBodyParams,
) -> HalfDayResonanceParams {
    // Precompute common quantities
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();

    // Define constants
    let c22s22 = 1.7891679e-6;
    let c32s32 = 3.7393792e-7;
    let c44s44 = 7.3636953e-9;
    let c52s52 = 1.1428639e-7;
    let c54s54 = 2.1765803e-9;

    // Calculate functions of inclination
    let f220 = (3. / 4.) * (1. + cos_i0).powi(2);
    let f221 = (3. / 2.) * sin_i0.powi(2);
    let f321 = (15. / 8.) * sin_i0 * (1. - 2. * cos_i0 - 3. * cos_i0.powi(2));
    let f322 = (-15. / 8.) * sin_i0 * (1. + 2. * cos_i0 - 3. * cos_i0.powi(2));
    let f441 = (105. / 4.) * sin_i0.powi(2) * (1. + cos_i0).powi(2);
    let f442 = (315. / 8.) * sin_i0.powi(4);
    let f522 = (315. / 32.)
        * (sin_i0.powi(3) - 2. * sin_i0.powi(3) * cos_i0 - 5. * sin_i0.powi(3) * cos_i0.powi(2)
            + sin_i0 * ((-2. / 3.) + (4. / 3.) * cos_i0 + 2. * cos_i0.powi(2)));
    let f523 = (105. / 16.)
        * sin_i0
        * (1. + 2. * cos_i0
            - 3. * cos_i0.powi(2)
            - (3. / 2.) * sin_i0.powi(2) * (1. + 2. * cos_i0 - 5. * cos_i0.powi(2)));
    let f542 = (945. / 32.)
        * sin_i0
        * (2. - 8. * cos_i0 + cos_i0.powi(2) * (-12. + 8. * cos_i0 + 10. * cos_i0.powi(2)));
    let f543 = (945. / 32.)
        * sin_i0
        * (cos_i0.powi(2) * (12. + 8. * cos_i0 - 10. * cos_i0.powi(2)) - 2. - 8. * cos_i0);

    // Calculate functions of eccentricity
    let g211: f64;
    let g201 = -0.306 - 0.44 * (brouwer0.e - 0.64);
    let g310: f64;
    let g322: f64;
    let g410: f64;
    let g422: f64;
    let g520: f64;
    let g521: f64;
    let g532: f64;
    let g533: f64;
    if brouwer0.e <= 0.65 {
        g211 = 3.616 - 13.247 * brouwer0.e + 16.29 * brouwer0.e.powi(2);
        g310 = -19.302 + 117.39 * brouwer0.e - 228.419 * brouwer0.e.powi(2)
            + 156.591 * brouwer0.e.powi(3);
        g322 = -18.9068 + 109.7927 * brouwer0.e - 214.6334 * brouwer0.e.powi(2)
            + 146.5816 * brouwer0.e.powi(3);
        g410 = -41.122 + 242.694 * brouwer0.e - 471.094 * brouwer0.e.powi(2)
            + 313.953 * brouwer0.e.powi(3);
        g422 = -146.407 + 841.88 * brouwer0.e - 1629.014 * brouwer0.e.powi(2)
            + 1083.435 * brouwer0.e.powi(3);
        g520 = -532.114 + 3017.977 * brouwer0.e - 5740.032 * brouwer0.e.powi(2)
            + 3708.276 * brouwer0.e.powi(3);
    } else {
        g211 = -72.099 + 331.819 * brouwer0.e - 508.738 * brouwer0.e.powi(2)
            + 266.724 * brouwer0.e.powi(3);
        g310 = -346.844 + 1582.851 * brouwer0.e - 2415.925 * brouwer0.e.powi(2)
            + 1246.113 * brouwer0.e.powi(3);
        g322 = -342.585 + 1554.908 * brouwer0.e - 2366.899 * brouwer0.e.powi(2)
            + 1215.972 * brouwer0.e.powi(3);
        g410 = -1052.797 + 4758.686 * brouwer0.e - 7193.992 * brouwer0.e.powi(2)
            + 3651.957 * brouwer0.e.powi(3);
        g422 = -3581.69 + 16178.11 * brouwer0.e - 24462.77 * brouwer0.e.powi(2)
            + 12422.52 * brouwer0.e.powi(3);
        if brouwer0.e < 0.715 {
            g520 = 1464.74 - 4664.75 * brouwer0.e + 3763.64 * brouwer0.e.powi(2);
        } else {
            g520 = -5149.66 + 29936.92 * brouwer0.e - 54087.36 * brouwer0.e.powi(2)
                + 31324.56 * brouwer0.e.powi(3);
        }
    }
    if brouwer0.e < 0.7 {
        g521 = -822.71072 + 4568.6173 * brouwer0.e - 8491.4146 * brouwer0.e.powi(2)
            + 5337.524 * brouwer0.e.powi(3);
        g532 = -853.666 + 4690.25 * brouwer0.e - 8624.77 * brouwer0.e.powi(2)
            + 5341.4 * brouwer0.e.powi(3);
        g533 = -919.2277 + 4988.61 * brouwer0.e - 9064.77 * brouwer0.e.powi(2)
            + 5542.21 * brouwer0.e.powi(3);
    } else {
        g521 = -51752.104 + 218913.95 * brouwer0.e - 309468.16 * brouwer0.e.powi(2)
            + 146349.42 * brouwer0.e.powi(3);
        g532 = -40023.88 + 170470.89 * brouwer0.e - 242699.48 * brouwer0.e.powi(2)
            + 115605.82 * brouwer0.e.powi(3);
        g533 = -37995.78 + 161616.52 * brouwer0.e - 229838.2 * brouwer0.e.powi(2)
            + 109377.94 * brouwer0.e.powi(3);
    }

    // Calculate the quadruples
    let d2201 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(2) * (c22s22 * f220 * g201);
    let d2211 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(2) * (c22s22 * f221 * g211);
    let d3210 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(3) * (c32s32 * f321 * g310);
    let d3222 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(3) * (c32s32 * f322 * g322);
    let d5220 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c52s52 * f522 * g520);
    let d5232 = 3. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c52s52 * f523 * g532);
    let d4422 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(4) * (c44s44 * f442 * g422); // 2x typo in Hoots et al 2004
    let d5421 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c54s54 * f542 * g521); // 2x typo in Hoots et al 2004
    let d5433 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(5) * (c54s54 * f543 * g533); // Typo in Hoots et al 2004
    let d4410 = 6. * brouwer0.n.powi(2) / brouwer0.a.powi(4) * (c44s44 * f441 * g410); // 2x typo in Hoots et al 2004

    // Calculate the initial value for the auxilary variable lam0
    let theta_g = calc_theta_g(jd0, jdfrac0);
    let lam0 = (brouwer0.m + 2. * brouwer0.raan - 2. * theta_g).rem_euclid(2.0 * PI);
    let lam0_dot = zonal_params.m_dot
        + (lunar_params.m_dot + solar_params.m_dot)
        + 2. * zonal_params.raan_dot
        + 2. * (lunar_params.raan_dot + solar_params.raan_dot)
        - 2. * RPTIM;

    // Store resonance parameters
    HalfDayResonanceParams {
        theta_g,
        lam0,
        lam0_dot,
        d2201,
        d2211,
        d3210,
        d3222,
        d5220,
        d5232,
        d4422,
        d5421,
        d5433,
        d4410,
    }
}

/// Initialize the whole day resonance effects of Earth's gravity
///
/// Initializes 24-hour resonance coefficients (`delta1`-`delta3`, `lambda31`/`lambda22`/`lambda33`)
/// and the auxiliary longitude `lambda0` used by the Euler-Maclaurin deep-space
/// integrator.
///
/// # Arguments
/// * `jd0` - The Julian date at epoch \[days\]
/// * `jdfrac0` - The fractional Julian date at epoch \[days\]
/// * `brouwer0` - The Brouwer mean elements at epoch
/// * `zonal_params` - The Earth zonal harmonics parameters
/// * `lunar_params` - The Lunar third body parameters
/// * `solar_params` - The Solar third body parameters
///
/// # Returns
/// * `WholeDayResonanceParams` - The whole day resonance parameters
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn init_earth_gravity_resonance_wholeday(
    jd0: f64,
    jdfrac0: f64,
    brouwer0: &BrouwerMeanElements,
    zonal_params: &EarthZonalParams,
    lunar_params: &ThirdBodyParams,
    solar_params: &ThirdBodyParams,
) -> WholeDayResonanceParams {
    // Precompute common quantities
    let cos_i0 = brouwer0.i.cos();
    let sin_i0 = brouwer0.i.sin();

    // Define constants
    let q31 = 2.1460748e-6;
    let q22 = 1.7891679e-6;
    let q33 = 2.2123015e-7;
    let lam31 = 0.13130908;
    let lam22 = 2.88431980;
    let lam33 = 0.37448087;

    // Calculate functions of inclination
    let f220 = (3. / 4.) * (1. + cos_i0).powi(2);
    let f311 = (15. / 16.) * sin_i0.powi(2) * (1. + 3. * cos_i0) - (3. / 4.) * (1. + cos_i0);
    let f330 = (15. / 8.) * (1. + cos_i0).powi(3);

    // Calculate functions of eccentricity
    let g200 = 1. - (5. / 2.) * brouwer0.e.powi(2) + (13. / 16.) * brouwer0.e.powi(4);
    let g310 = 1. + 2. * brouwer0.e.powi(2);
    let g300 = 1. - 6. * brouwer0.e.powi(2) + (423. / 64.) * brouwer0.e.powi(4);

    // Calculate coefficients of the resonance terms
    let delta1 = (3. * brouwer0.n.powi(2) / brouwer0.a.powi(3)) * f311 * g310 * q31;
    let delta2 = (6. * brouwer0.n.powi(2) / brouwer0.a.powi(2)) * f220 * g200 * q22;
    let delta3 = (9. * brouwer0.n.powi(2) / brouwer0.a.powi(3)) * f330 * g300 * q33;

    // Calculate the initial value for the auxilary variable lam0
    let theta_g = calc_theta_g(jd0, jdfrac0);
    let lam0 = brouwer0.m + brouwer0.raan + brouwer0.omega - theta_g;
    let lam0_dot_1 = zonal_params.m_dot
        + (lunar_params.m_dot + solar_params.m_dot)
        + zonal_params.raan_dot
        + (lunar_params.raan_dot + solar_params.raan_dot);
    let lam0_dot_2 =
        zonal_params.omega_dot + (lunar_params.omega_dot + solar_params.omega_dot) - RPTIM;
    let lam0_dot = lam0_dot_1 + lam0_dot_2;

    // Store resonance parameters
    WholeDayResonanceParams {
        theta_g,
        lam0,
        lam0_dot,
        lam31,
        lam22,
        lam33,
        delta1,
        delta2,
        delta3,
    }
}

/// Calculate Greenwich mean sidereal time (GMST) / longitude of Greenwich at a Julian date.
///
/// Used as `theta_g` when initializing 12 h / 24 h resonance terms. Also gives the rotation
/// angle between TEME and an Earth-fixed frame (neglecting polar motion).
///
/// # Arguments
/// * `jd0` - Julian day (integer part) \[days\]
/// * `jdfrac0` - Julian day fraction \[days\]
///
/// # Returns
/// * `theta_g` - GMST \[rad\], wrapped to \[0, 2 * pi)
///
/// # Examples
/// ```rust
/// use mako_sgp4::sgp4::calc_theta_g;
///
/// // GMST at the J2000.0 epoch (2000-01-01 12:00 UT1)
/// let theta_g = calc_theta_g(2451545.0, 0.0);
///
/// // GMST is about 280.46 degrees
/// assert!((theta_g.to_degrees() - 280.46061837).abs() < 1e-8);
/// ```
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
pub fn calc_theta_g(jd0: f64, jdfrac0: f64) -> f64 {
    // Calculate the Julian centuries since J2000.0
    let tut1 = (jd0 + jdfrac0 - 2451545.0) / 36525.0; // [centuries]

    // Calculate the Greenwich sidereal time in seconds
    let temp = -6.2e-6 * tut1.powi(3)
        + 0.093104 * tut1.powi(2)
        + (876600.0 * 3600.0 + 8640184.812866) * tut1
        + 67310.54841; // [seconds]

    // Calculate the Greenwich sidereal time in radians
    (deg2rad(temp / 240.0)).rem_euclid(2.0 * PI) // [radians] 360/86400 = 1/240 degrees per second
}

/// Propagate an initialized [`Sgp4`] model to a UTC [`DateTime`].
///
/// Converts `datetime` to Julian date, forms minutes since the TLE epoch, then
/// calls [`sgp4_prop_delta`].
///
/// # Arguments
/// * `sgp4` - Initialized propagator
/// * `datetime` - Propagation epoch in UTC
///
/// # Returns
/// * `Ok(StateVector)` - TEME position \[km\] and velocity \[km/s\]
/// * `Err(Sgp4Error)` - If Julian-day conversion fails or intermediate elements are non-physical
///
/// # Errors
/// * [`Sgp4Error::InvalidDateTime`] - If `datetime` is not UTC or Julian-day conversion fails
/// * [`Sgp4Error::InvalidMeanMotion`] - If mean motion is less than or equal to zero
/// * [`Sgp4Error::InvalidMeanEccentricity`] - If mean eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidPerturbedEccentricity`] - If perturbed eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidSemilatusRectum`] - If the semilatus rectum is less than zero
/// * [`Sgp4Error::SatelliteDecayed`] - If the satellite has decayed
/// * [`Sgp4Error::NonFiniteState`] - If the propagated state is NaN or infinite
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::CoordinateFrame;
/// use mako_sgp4::gp::from_tle_lines;
/// use mako_sgp4::sgp4::sgp4_prop_datetime;
///
/// // Parse a TLE and initialize the SGP4 propagator
/// let tle_line0 = "ISS (ZARYA)";
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(tle_line1, tle_line2, Some(tle_line0)).unwrap();
///
/// // Propagate to the TLE epoch
/// let state_vector = sgp4_prop_datetime(&sgp4, &sgp4.gp.epoch_datetime).unwrap();
/// assert_eq!(state_vector.coordinate_frame, CoordinateFrame::TEME);
/// assert!(state_vector.r_x.is_finite());
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn sgp4_prop_datetime(sgp4: &Sgp4, datetime: &DateTime) -> Result<StateVector, Sgp4Error> {
    // Convert datetime to Julian day format
    let (jd_prop, jdfrac_prop) = utc2jday(datetime).map_err(Sgp4Error::InvalidDateTime)?;

    // Get minutes since epoch. Subtract whole days and day-fractions separately
    // to avoid floating-point rounding errors
    let delta_t = ((jd_prop - sgp4.jd0) + (jdfrac_prop - sgp4.jdfrac0)) * 1440.;

    // Propagate the state vector
    sgp4_prop_delta(sgp4, delta_t)
}

/// Propagate an initialized [`Sgp4`] model by `delta_t` minutes from the TLE epoch.
///
/// Applies secular drag/zonal updates (and deep-space lunar/solar + resonance when
/// applicable), long-period periodics, then the short-period Kepler / J2 solution.
/// Returns [`Sgp4Error`] on non-physical intermediate elements (e.g. mean motion <= 0,
/// eccentricity out of range, decayed radius), matching Vallado-style checks.
///
/// # Arguments
/// * `sgp4` - Initialized propagator
/// * `delta_t` - Minutes since epoch (may be negative)
///
/// # Returns
/// * `Ok(StateVector)` - TEME position \[km\] and velocity \[km/s\]
/// * `Err(Sgp4Error)` - If intermediate orbital elements become non-physical
///
/// # Errors
/// * [`Sgp4Error::InvalidMeanMotion`] - If mean motion is less than or equal to zero
/// * [`Sgp4Error::InvalidMeanEccentricity`] - If mean eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidPerturbedEccentricity`] - If perturbed eccentricity is outside the range 0.0 to 1.0
/// * [`Sgp4Error::InvalidSemilatusRectum`] - If the semilatus rectum is less than zero
/// * [`Sgp4Error::SatelliteDecayed`] - If the satellite has decayed
/// * [`Sgp4Error::NonFiniteState`] - If the propagated state is NaN or infinite
///
/// # Examples
/// ```rust
/// use mako_sgp4::common::CoordinateFrame;
/// use mako_sgp4::gp::from_tle_lines;
/// use mako_sgp4::sgp4::sgp4_prop_delta;
///
/// // Parse a TLE and initialize the SGP4 propagator
/// let tle_line0 = "ISS (ZARYA)";
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(tle_line1, tle_line2, Some(tle_line0)).unwrap();
///
/// // Propagate 6 hours past epoch
/// let state_vector = sgp4_prop_delta(&sgp4, 360.0).unwrap();
/// assert_eq!(state_vector.coordinate_frame, CoordinateFrame::TEME);
/// assert!(state_vector.r_x.is_finite());
/// ```
///
/// # References
/// - [Revisiting Spacetrack Report #3: Rev 3 by Vallado et al](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev3.pdf)
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
pub fn sgp4_prop_delta(sgp4: &Sgp4, delta_t: f64) -> Result<StateVector, Sgp4Error> {
    // Create mutable variables for the orbital elements
    let mut m: f64;
    let mut omega: f64;
    let mut raan: f64;
    let mut e = sgp4.brouwer0.e;
    let mut i = sgp4.brouwer0.i;
    let mut n = sgp4.brouwer0.n;

    // Account for Earth zonal gravity and partial atmospheric drag effects
    let m_df = sgp4.brouwer0.m + sgp4.brouwer0.n * delta_t + sgp4.zonal_params.m_dot * delta_t;
    let omega_df = sgp4.brouwer0.omega + sgp4.zonal_params.omega_dot * delta_t;
    let raan_df = sgp4.brouwer0.raan + sgp4.zonal_params.raan_dot * delta_t;

    // Neglect delta_omega and delta_m if deep space or perigee height is less than 220 km, or e <= 1e-4
    let delta_omega: f64;
    let delta_m: f64;
    if sgp4.deep_space || sgp4.atm_params.hp < 220. || sgp4.brouwer0.e <= 1.0e-4 {
        delta_omega = 0.;
        delta_m = 0.;
    } else {
        delta_omega = sgp4.gp.bstar * sgp4.atm_params.c3 * sgp4.brouwer0.omega.cos() * delta_t;
        delta_m = (-2. / 3.)
            * (sgp4.atm_params.q0 - sgp4.atm_params.s).powi(4)
            * sgp4.gp.bstar
            * sgp4.atm_params.zeta.powi(4)
            * (1. / (sgp4.brouwer0.e * sgp4.atm_params.eta))
            * ((1. + sgp4.atm_params.eta * m_df.cos()).powi(3)
                - (1. + sgp4.atm_params.eta * sgp4.brouwer0.m.cos()).powi(3));
    }

    m = m_df + delta_omega + delta_m;
    omega = omega_df - delta_omega - delta_m;
    raan = raan_df
        - (21. / 2.)
            * (sgp4.brouwer0.n * sgp4.wgs.k2 * sgp4.brouwer0.theta
                / (sgp4.brouwer0.a.powi(2) * sgp4.brouwer0.beta.powi(2)))
            * sgp4.atm_params.c1
            * delta_t.powi(2);

    // Account for Lunar and Solar third body effects
    if sgp4.deep_space {
        m += (sgp4.lunar_params.m_dot + sgp4.solar_params.m_dot) * delta_t;
        omega += (sgp4.lunar_params.omega_dot + sgp4.solar_params.omega_dot) * delta_t;
        raan += (sgp4.lunar_params.raan_dot + sgp4.solar_params.raan_dot) * delta_t;
        e += (sgp4.lunar_params.e_dot + sgp4.solar_params.e_dot) * delta_t;
        i += (sgp4.lunar_params.i_dot + sgp4.solar_params.i_dot) * delta_t;
    }

    // Account for the whole and half day resonance effects of Earth's gravity
    // Vallado dspace: +/-720 min Euler-Maclaurin steps (works for negative tsince)
    if sgp4.half_day_resonance {
        let mut lami = sgp4.half_day_resonance_params.lam0;
        let mut ni = sgp4.brouwer0.n;
        let mut lami_dot: f64;
        let mut ni_dot: f64;
        let mut lami_ddot: f64;
        let mut ni_ddot: f64;
        let step = if delta_t >= 0.0 { 720.0 } else { -720.0 };
        let em_steps = (delta_t / step).floor() as i32;
        let t_em = delta_t - em_steps as f64 * step;

        (lami_dot, ni_dot, lami_ddot, ni_ddot) = half_day_euler_maclaurin_step(
            lami,
            ni,
            sgp4.brouwer0.omega,
            &sgp4.half_day_resonance_params,
        );

        for em_step in 0..em_steps {
            // h^2/2 term uses |720|^2; linear term uses signed step (Vallado step2 = 259200)
            lami += lami_dot * step + 0.5 * lami_ddot * 518400.;
            ni += ni_dot * step + 0.5 * ni_ddot * 518400.;

            let omegai =
                sgp4.brouwer0.omega + sgp4.zonal_params.omega_dot * (em_step + 1) as f64 * step;
            (lami_dot, ni_dot, lami_ddot, ni_ddot) =
                half_day_euler_maclaurin_step(lami, ni, omegai, &sgp4.half_day_resonance_params);
        }

        lami = lami + (lami_dot * t_em) + (0.5 * lami_ddot * t_em.powi(2));
        ni = ni + (ni_dot * t_em) + (0.5 * ni_ddot * t_em.powi(2));

        let theta_t =
            (sgp4.half_day_resonance_params.theta_g + RPTIM * delta_t).rem_euclid(2.0 * PI);
        n = ni;
        m = lami - 2. * raan + 2. * theta_t;
    } else if sgp4.whole_day_resonance {
        let mut lami = sgp4.whole_day_resonance_params.lam0;
        let mut ni = sgp4.brouwer0.n;
        let mut lami_dot: f64;
        let mut ni_dot: f64;
        let mut lami_ddot: f64;
        let mut ni_ddot: f64;
        let step = if delta_t >= 0.0 { 720.0 } else { -720.0 };
        let em_steps = (delta_t / step).floor() as i32;
        let t_em = delta_t - em_steps as f64 * step;

        (lami_dot, ni_dot, lami_ddot, ni_ddot) =
            whole_day_euler_maclaurin_step(lami, ni, &sgp4.whole_day_resonance_params);

        for _ in 0..em_steps {
            lami += lami_dot * step + 0.5 * lami_ddot * 518400.;
            ni += ni_dot * step + 0.5 * ni_ddot * 518400.;

            (lami_dot, ni_dot, lami_ddot, ni_ddot) =
                whole_day_euler_maclaurin_step(lami, ni, &sgp4.whole_day_resonance_params);
        }

        lami = lami + (lami_dot * t_em) + (0.5 * lami_ddot * t_em.powi(2));
        ni = ni + (ni_dot * t_em) + (0.5 * ni_ddot * t_em.powi(2));

        let theta_t =
            (sgp4.whole_day_resonance_params.theta_g + RPTIM * delta_t).rem_euclid(2.0 * PI);
        n = ni;
        m = lami - raan - omega + theta_t;
    }

    // Mean motion must remain positive before the drag semi-major-axis update
    if n <= 0.0 {
        return Err(Sgp4Error::InvalidMeanMotion);
    }

    // Account for remaining atmospheric drag effects
    let a: f64;
    let il_atm: f64 = if sgp4.deep_space || sgp4.atm_params.hp < 220. {
        e += -sgp4.gp.bstar * (sgp4.atm_params.c4 * delta_t);
        let a_1 = 1. - sgp4.atm_params.c1 * delta_t; // Drop quadratic term, different from Hoots et al 2004
        a = (sgp4.wgs.ke / n).powf(2. / 3.) * a_1.powi(2);
        n = sgp4.wgs.ke / a.powf(1.5);
        let il_1 = 3. / 2. * sgp4.atm_params.c1 * delta_t.powi(2);
        sgp4.brouwer0.n * il_1
    } else {
        e += -sgp4.gp.bstar
            * (sgp4.atm_params.c4 * delta_t
                + sgp4.atm_params.c5 * (m.sin() - sgp4.brouwer0.m.sin()));
        let a_1 = 1. - sgp4.atm_params.c1 * delta_t - sgp4.atm_params.d2 * delta_t.powi(2);
        let a_2 = -sgp4.atm_params.d3 * delta_t.powi(3) - sgp4.atm_params.d4 * delta_t.powi(4);
        a = (sgp4.wgs.ke / n).powf(2. / 3.) * (a_1 + a_2).powi(2);
        n = sgp4.wgs.ke / a.powf(1.5);
        let il_1 = 3. / 2. * sgp4.atm_params.c1 * delta_t.powi(2);
        let il_2 = (sgp4.atm_params.d2 + 2. * sgp4.atm_params.c1.powi(2)) * delta_t.powi(3);
        let il_3 = 1. / 4.
            * (3. * sgp4.atm_params.d3
                + 12. * sgp4.atm_params.c1 * sgp4.atm_params.d2
                + 10. * sgp4.atm_params.c1.powi(3))
            * delta_t.powi(4);
        let il_4 = 1. / 5.
            * (3. * sgp4.atm_params.d4
                + 12. * sgp4.atm_params.c1 * sgp4.atm_params.d3
                + 6. * sgp4.atm_params.d2.powi(2)
                + 30. * sgp4.atm_params.c1.powi(2) * sgp4.atm_params.d2
                + 15. * sgp4.atm_params.c1.powi(4))
            * delta_t.powi(5);
        sgp4.brouwer0.n * (il_1 + il_2 + il_3 + il_4)
    };

    // Mean eccentricity check before the near-zero floor (Vallado)
    if !(-0.001..1.0).contains(&e) {
        return Err(Sgp4Error::InvalidMeanEccentricity);
    }

    // Vallado eccentricity guard after atmospheric drag
    if e < 1.0e-6 {
        e = 1.0e-6;
    }

    // Vallado mean-element recovery before lunar-solar periodics. Use `%` (C fmod, keeps sign)
    // so the Lyddane node continuity check below sees the same wrapped angles as Vallado
    m += il_atm;
    let lm = (m + omega + raan) % (2.0 * PI);
    raan %= 2.0 * PI;
    omega %= 2.0 * PI;
    m = (lm - omega - raan) % (2.0 * PI);

    // Account for long-period periodic effects of lunar and solar gravity
    if sgp4.deep_space {
        let m_m = sgp4.lunar_params.m + sgp4.lunar_params.n * delta_t;
        let m_s = sgp4.solar_params.m + sgp4.solar_params.n * delta_t;
        let f_m = m_m + 2. * sgp4.lunar_params.e * m_m.sin();
        let f_s = m_s + 2. * sgp4.solar_params.e * m_s.sin();
        let f2_m = 0.5 * f_m.sin().powi(2) - 0.25;
        let f2_s = 0.5 * f_s.sin().powi(2) - 0.25;
        let f3_m = -0.5 * f_m.sin() * f_m.cos();
        let f3_s = -0.5 * f_s.sin() * f_s.cos();
        let delta_e_m = -(30. * sgp4.brouwer0.beta * sgp4.lunar_params.c * sgp4.brouwer0.e
            / sgp4.brouwer0.n)
            * (f2_m
                * (sgp4.lunar_params.x2 * sgp4.lunar_params.x3
                    + sgp4.lunar_params.x1 * sgp4.lunar_params.x4)
                + f3_m
                    * (sgp4.lunar_params.x2 * sgp4.lunar_params.x4
                        - sgp4.lunar_params.x1 * sgp4.lunar_params.x3));
        let delta_e_s = -(30. * sgp4.brouwer0.beta * sgp4.solar_params.c * sgp4.brouwer0.e
            / sgp4.brouwer0.n)
            * (f2_s
                * (sgp4.solar_params.x2 * sgp4.solar_params.x3
                    + sgp4.solar_params.x1 * sgp4.solar_params.x4)
                + f3_s
                    * (sgp4.solar_params.x2 * sgp4.solar_params.x4
                        - sgp4.solar_params.x1 * sgp4.solar_params.x3));
        let delta_i_m = -(sgp4.lunar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta)
            * (f2_m * sgp4.lunar_params.z12
                + f3_m * (sgp4.lunar_params.z13 - sgp4.lunar_params.z11));
        let delta_i_s = -(sgp4.solar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta)
            * (f2_s * sgp4.solar_params.z12
                + f3_s * (sgp4.solar_params.z13 - sgp4.solar_params.z11));
        let delta_m_m = -(2. * sgp4.lunar_params.c / sgp4.brouwer0.n)
            * (f2_m * sgp4.lunar_params.z2 + f3_m * (sgp4.lunar_params.z3 - sgp4.lunar_params.z1)
                - 3. * sgp4.lunar_params.e * f_m.sin() * (7. + 3. * sgp4.brouwer0.e.powi(2)));
        let delta_m_s = -(2. * sgp4.solar_params.c / sgp4.brouwer0.n)
            * (f2_s * sgp4.solar_params.z2 + f3_s * (sgp4.solar_params.z3 - sgp4.solar_params.z1)
                - 3. * sgp4.solar_params.e * f_s.sin() * (7. + 3. * sgp4.brouwer0.e.powi(2)));
        let delta_raan_m = (sgp4.lunar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta)
            * (f2_m * sgp4.lunar_params.z22
                + f3_m * (sgp4.lunar_params.z23 - sgp4.lunar_params.z21)); // / sgp4.lunar_params.i.sin();
        let delta_raan_s = (sgp4.solar_params.c / sgp4.brouwer0.n / sgp4.brouwer0.beta)
            * (f2_s * sgp4.solar_params.z22
                + f3_s * (sgp4.solar_params.z23 - sgp4.solar_params.z21)); // / sgp4.solar_params.i.sin();
        let delta_omega_m = (2. * sgp4.brouwer0.beta * sgp4.lunar_params.c / sgp4.brouwer0.n)
            * (f2_m * sgp4.lunar_params.z32
                + f3_m * (sgp4.lunar_params.z33 - sgp4.lunar_params.z31)
                - 9. * sgp4.lunar_params.e * f_m.sin()); // - delta_raan_m * sgp4.lunar_params.i.cos();
        let delta_omega_s = (2. * sgp4.brouwer0.beta * sgp4.solar_params.c / sgp4.brouwer0.n)
            * (f2_s * sgp4.solar_params.z32
                + f3_s * (sgp4.solar_params.z33 - sgp4.solar_params.z31)
                - 9. * sgp4.solar_params.e * f_s.sin()); // - delta_raan_s * sgp4.solar_params.i.cos();
        let delta_e_ls = delta_e_m + delta_e_s;
        let delta_i_ls = delta_i_m + delta_i_s;
        let delta_m_ls = delta_m_m + delta_m_s;
        let delta_raan_ls = delta_raan_m + delta_raan_s;
        let delta_omega_ls = delta_omega_m + delta_omega_s;

        e += delta_e_ls;
        if !(0.0..1.0).contains(&e) {
            return Err(Sgp4Error::InvalidPerturbedEccentricity);
        }
        i += delta_i_ls;

        // Vallado threshold (i >= 0.2 on the perturbed inclination), differs from Hoots et al 2004
        if i >= 0.2 {
            // Notation is confusing in paper, delta_omega_ls stores more than just delta_omega
            // Same for delta_raan_ls, stores more than just delta_raan
            raan += delta_raan_ls / i.sin();
            omega += delta_omega_ls - delta_raan_ls * i.cos() / i.sin();
            m += delta_m_ls;
        } else {
            // Lyddane modification for inclinations below 0.2 rad (legacy)
            let alpha = i.sin() * raan.sin()
                + raan.cos() * delta_raan_ls
                + i.cos() * raan.sin() * delta_i_ls;
            let beta = i.sin() * raan.cos() - raan.sin() * delta_raan_ls
                + i.cos() * raan.cos() * delta_i_ls;
            let raan_old = raan;

            // Use Vallado's modification which is more numerically stable when i ~ 0
            let m_omega_raan =
                m + omega + delta_m_ls + delta_omega_ls + (i.cos() - delta_i_ls * i.sin()) * raan;

            // Calculate RAAN
            raan = alpha.atan2(beta);
            // Maintain RAAN continuity across atan2 branch cuts
            if (raan_old - raan).abs() > PI {
                if raan < raan_old {
                    raan += 2.0 * PI;
                } else {
                    raan -= 2.0 * PI;
                }
            }

            // Calculate omega
            m += delta_m_ls;
            omega = m_omega_raan - m - i.cos() * raan;
        }
    }

    // Vallado keep inclination in [0, pi] after lunisolar periodics
    if sgp4.deep_space && i < 0.0 {
        i = -i;
        raan += PI;
        omega -= PI;
    }

    let il = m + omega + raan;

    // Account for long-period periodic effects of Earth's gravity
    let beta_update = (1. - e.powi(2)).sqrt();
    let a30 = -sgp4.wgs.j3; // [Earth Radii^3]
    let axn = e * omega.cos();
    // Vallado floor on 1 + cos(i) so retrograde equatorial orbits (i = 180 deg) do not divide by zero
    let mut one_plus_cos_i = 1. + i.cos();
    if one_plus_cos_i.abs() < 1.5e-12 {
        one_plus_cos_i = 1.5e-12;
    }
    let ill = a30 * i.sin() / (8. * sgp4.wgs.k2 * a * beta_update.powi(2))
        * e
        * omega.cos()
        * (3. + 5. * i.cos())
        / one_plus_cos_i;
    let aynl = a30 * i.sin() / (4. * sgp4.wgs.k2 * a * beta_update.powi(2));
    let ilt = il + ill;
    let ayn = e * omega.sin() + aynl;

    // Account for short-period periodic effects of Earth's gravity (solve Kepler's equation)
    let u = (ilt - raan).rem_euclid(2.0 * PI);
    let mut e_omega = u;
    let mut delta_e_omega: f64;

    // Newton-Raphson iteration to solve Kepler's equation (10 iterations max per Vallado)
    for _ in 0..10 {
        delta_e_omega = (u - ayn * e_omega.cos() + axn * e_omega.sin() - e_omega)
            / (1. - ayn * e_omega.sin() - axn * e_omega.cos());

        // Protect against oversized steps
        if delta_e_omega.abs() >= 0.95 {
            if delta_e_omega > 0.0 {
                delta_e_omega = 0.95;
            } else {
                delta_e_omega = -0.95;
            }
        }

        e_omega += delta_e_omega;

        // Verify convergence
        if delta_e_omega.abs() < 1e-12 {
            break;
        }
    }

    // Return position and velocity vectors in the TEME frame
    e = (axn.powi(2) + ayn.powi(2)).sqrt();
    let pl = a * (1. - e.powi(2));
    if pl < 0.0 {
        return Err(Sgp4Error::InvalidSemilatusRectum);
    }
    let cos_ecc_anomaly = (axn * e_omega.cos() + ayn * e_omega.sin()) / e;
    let sin_ecc_anomaly = (axn * e_omega.sin() - ayn * e_omega.cos()) / e;
    let r = a * (1. - e * cos_ecc_anomaly);
    let r_dot = sgp4.wgs.ke * a.sqrt() * e * sin_ecc_anomaly / r;
    let r_f_dot = sgp4.wgs.ke * pl.sqrt() / r;
    let cos_u = a / r
        * (e_omega.cos() - axn + ayn * (e * sin_ecc_anomaly) / (1. + (1. - e.powi(2)).sqrt()));
    let sin_u = a / r
        * (e_omega.sin() - ayn - axn * (e * sin_ecc_anomaly) / (1. + (1. - e.powi(2)).sqrt()));
    let u = sin_u.atan2(cos_u);
    let delta_r = sgp4.wgs.k2 / (2. * pl) * (1. - i.cos().powi(2)) * (2. * u).cos();
    let delta_u = -sgp4.wgs.k2 / (4. * pl.powi(2)) * (7. * i.cos().powi(2) - 1.) * (2. * u).sin();
    let delta_raan = 3. * sgp4.wgs.k2 * i.cos() / (2. * pl.powi(2)) * (2. * u).sin();
    let delta_i = 3. * sgp4.wgs.k2 * i.cos() / (2. * pl.powi(2)) * i.sin() * (2. * u).cos();
    let delta_r_dot = -sgp4.wgs.k2 * n / pl * (1. - i.cos().powi(2)) * (2. * u).sin();
    let delta_r_f_dot = sgp4.wgs.k2 * n / pl
        * ((1. - i.cos().powi(2)) * (2. * u).cos() - 3. / 2. * (1. - 3. * i.cos().powi(2)));
    let rk = r
        * (1.
            - 3. / 2. * sgp4.wgs.k2 * (1. - e.powi(2)).sqrt() / pl.powi(2)
                * (3. * i.cos().powi(2) - 1.))
        + delta_r;
    if rk < 1.0 {
        return Err(Sgp4Error::SatelliteDecayed);
    }
    let uk = u + delta_u;
    let raan_k = raan + delta_raan;
    let i_k = i + delta_i;
    let r_dot_k = r_dot + delta_r_dot;
    let r_f_dot_k = r_f_dot + delta_r_f_dot;

    let mx = -raan_k.sin() * i_k.cos();
    let my = raan_k.cos() * i_k.cos();
    let mz = i_k.sin();

    let nx = raan_k.cos();
    let ny = raan_k.sin();
    let nz = 0.;

    let ux = mx * uk.sin() + nx * uk.cos();
    let uy = my * uk.sin() + ny * uk.cos();
    let uz = mz * uk.sin() + nz * uk.cos();

    let vx = mx * uk.cos() - nx * uk.sin();
    let vy = my * uk.cos() - ny * uk.sin();
    let vz = mz * uk.cos() - nz * uk.sin();

    let rx = rk * ux * sgp4.wgs.r_earth_eq;
    let ry = rk * uy * sgp4.wgs.r_earth_eq;
    let rz = rk * uz * sgp4.wgs.r_earth_eq;

    let r_dot_x = (r_dot_k * ux + r_f_dot_k * vx) * sgp4.wgs.r_earth_eq / 60.;
    let r_dot_y = (r_dot_k * uy + r_f_dot_k * vy) * sgp4.wgs.r_earth_eq / 60.;
    let r_dot_z = (r_dot_k * uz + r_f_dot_k * vz) * sgp4.wgs.r_earth_eq / 60.;

    // Never return a non-finite state as a success
    if ![rx, ry, rz, r_dot_x, r_dot_y, r_dot_z]
        .iter()
        .all(|x| x.is_finite())
    {
        return Err(Sgp4Error::NonFiniteState);
    }

    Ok(StateVector {
        r_x: rx,
        r_y: ry,
        r_z: rz,
        v_x: r_dot_x,
        v_y: r_dot_y,
        v_z: r_dot_z,
        coordinate_frame: CoordinateFrame::TEME,
    })
}

/// Half day Euler-Maclaurin integration step
///
/// Evaluates `lambda_dot`, `n_dot`, `lambda_ddot`, and `n_ddot` for one 12-hour resonance integrator step
/// at the current auxiliary longitude, mean motion, and argument of perigee.
///
/// # Arguments
/// * `lami` - The auxilary variable at time i
/// * `ni` - The mean motion at time i
/// * `omegai` - The argument of perigee at time i
/// * `half_day_resonance_params` - The half day resonance parameters
///
/// # Returns
/// * `lami_dot` - The rate of change of the auxilary variable at time i+1
/// * `ni_dot` - The rate of change of the mean motion at time i+1
/// * `lami_ddot` - The 2nd derivative of the auxilary variable at time i+1
/// * `ni_ddot` - The 2nd derivative of the mean motion at time i+1
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn half_day_euler_maclaurin_step(
    lami: f64,
    ni: f64,
    omegai: f64,
    half_day_resonance_params: &HalfDayResonanceParams,
) -> (f64, f64, f64, f64) {
    // Define constants
    let g22 = 5.7686396;
    let g32 = 0.95240898;
    let g44 = 1.8014998;
    let g52 = 1.0508330;
    let g54 = 4.4108898;

    // Calculate the rate of change of the auxilary variable
    let lami_dot = ni + half_day_resonance_params.lam0_dot;

    // Calculate the rate of change of the mean motion
    let ni_dot_2201 = half_day_resonance_params.d2201 * (2. * omegai + lami - g22).sin();
    let ni_dot_2211 = half_day_resonance_params.d2211 * (lami - g22).sin();
    let ni_dot_3210 = half_day_resonance_params.d3210 * (omegai + lami - g32).sin();
    let ni_dot_3222 = half_day_resonance_params.d3222 * (-omegai + lami - g32).sin();
    let ni_dot_5220 = half_day_resonance_params.d5220 * (omegai + lami - g52).sin();
    let ni_dot_5232 = half_day_resonance_params.d5232 * (-omegai + lami - g52).sin();
    let ni_dot_4422 = half_day_resonance_params.d4422 * (2. * lami - g44).sin();
    let ni_dot_5421 = half_day_resonance_params.d5421 * (omegai + 2. * lami - g54).sin();
    let ni_dot_5433 = half_day_resonance_params.d5433 * (-omegai + 2. * lami - g54).sin();
    let ni_dot_4410 = half_day_resonance_params.d4410 * (2. * omegai + 2. * lami - g44).sin();
    let ni_dot = ni_dot_2201
        + ni_dot_2211
        + ni_dot_3210
        + ni_dot_3222
        + ni_dot_5220
        + ni_dot_5232
        + ni_dot_4422
        + ni_dot_5421
        + ni_dot_5433
        + ni_dot_4410;

    // Calculate the 2nd derivative of the auxilary variable
    let lami_ddot = ni_dot;

    // Calculate the 2nd derivative of the mean motion
    let ni_ddot_2201 = 1. * half_day_resonance_params.d2201 * (2. * omegai + lami - g22).cos();
    let ni_ddot_2211 = 1. * half_day_resonance_params.d2211 * (lami - g22).cos();
    let ni_ddot_3210 = 1. * half_day_resonance_params.d3210 * (omegai + lami - g32).cos();
    let ni_ddot_3222 = 1. * half_day_resonance_params.d3222 * (-omegai + lami - g32).cos();
    let ni_ddot_5220 = 1. * half_day_resonance_params.d5220 * (omegai + lami - g52).cos();
    let ni_ddot_5232 = 1. * half_day_resonance_params.d5232 * (-omegai + lami - g52).cos();
    let ni_ddot_4422 = 2. * half_day_resonance_params.d4422 * (2. * lami - g44).cos();
    let ni_ddot_5421 = 2. * half_day_resonance_params.d5421 * (omegai + 2. * lami - g54).cos();
    let ni_ddot_5433 = 2. * half_day_resonance_params.d5433 * (-omegai + 2. * lami - g54).cos();
    let ni_ddot_4410 = 2. * half_day_resonance_params.d4410 * (2. * omegai + 2. * lami - g44).cos();
    let ni_ddot = lami_dot
        * (ni_ddot_2201
            + ni_ddot_2211
            + ni_ddot_3210
            + ni_ddot_3222
            + ni_ddot_5220
            + ni_ddot_5232
            + ni_ddot_4422
            + ni_ddot_5421
            + ni_ddot_5433
            + ni_ddot_4410);

    (lami_dot, ni_dot, lami_ddot, ni_ddot)
}

/// Whole day Euler-Maclaurin integration step
///
/// Evaluates `lambda_dot`, `n_dot`, `lambda_ddot`, and `n_ddot` for one 24-hour resonance integrator step
/// at the current auxiliary longitude and mean motion.
///
/// # Arguments
/// * `lami` - The auxilary variable at time i
/// * `ni` - The mean motion at time i
/// * `whole_day_resonance_params` - The whole day resonance parameters
///
/// # Returns
/// * `lami_dot` - The rate of change of the auxilary variable at time i+1
/// * `ni_dot` - The rate of change of the mean motion at time i+1
/// * `lami_ddot` - The 2nd derivative of the auxilary variable at time i+1
/// * `ni_ddot` - The 2nd derivative of the mean motion at time i+1
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [History of Analytical Orbit Modeling in the U.S. Space Surveillance System by Hoots et al](https://arc.aiaa.org/doi/abs/10.2514/1.9161?journalCode=jgcd)
fn whole_day_euler_maclaurin_step(
    lami: f64,
    ni: f64,
    whole_day_resonance_params: &WholeDayResonanceParams,
) -> (f64, f64, f64, f64) {
    // Calculate the rate of change of the auxilary variable
    let lami_dot = ni + whole_day_resonance_params.lam0_dot;

    // Calculate the rate of change of the mean motion
    let ni_dot_1 =
        whole_day_resonance_params.delta1 * (lami - whole_day_resonance_params.lam31).sin();
    let ni_dot_2 =
        whole_day_resonance_params.delta2 * (2. * (lami - whole_day_resonance_params.lam22)).sin();
    let ni_dot_3 =
        whole_day_resonance_params.delta3 * (3. * (lami - whole_day_resonance_params.lam33)).sin();
    let ni_dot = ni_dot_1 + ni_dot_2 + ni_dot_3;

    // Calculate the 2nd derivative of the auxilary variable
    let lami_ddot = ni_dot;

    // Calculate the 2nd derivative of the mean motion
    let ni_ddot_1 =
        whole_day_resonance_params.delta1 * (lami - whole_day_resonance_params.lam31).cos();
    let ni_ddot_2 = 2.
        * whole_day_resonance_params.delta2
        * (2. * (lami - whole_day_resonance_params.lam22)).cos();
    let ni_ddot_3 = 3.
        * whole_day_resonance_params.delta3
        * (3. * (lami - whole_day_resonance_params.lam33)).cos();
    let ni_ddot = lami_dot * (ni_ddot_1 + ni_ddot_2 + ni_ddot_3);

    (lami_dot, ni_dot, lami_ddot, ni_ddot)
}

// ----------
// Unit Tests
// ----------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{from_omm_kvn_string, from_tle_string};
    use crate::time::Timezone;
    use serde::Deserialize;
    use std::collections::HashMap;
    use toml::from_str;

    // -------------------------------------------------------
    // Structs for deserializing the Vallado test case TOML
    // -------------------------------------------------------

    #[derive(Deserialize)]
    struct ValladoCases {
        test: HashMap<String, ValladoCase>,
    }

    #[derive(Deserialize)]
    struct ValladoCase {
        name: String,
        tle: String,
        #[serde(default)]
        #[allow(dead_code)]
        start_mins_from_epoch: f64,
        #[serde(default)]
        #[allow(dead_code)]
        end_mins_from_epoch: f64,
        #[serde(default)]
        #[allow(dead_code)]
        delta_time_mins: f64,
        ephem: String,
        #[serde(default)]
        exception: bool,
    }

    // One ephemeris row: minutes from epoch, TEME position and velocity, and UTC time.
    // A row without a calendar stamp uses the TLE epoch.
    struct EphemRow {
        t_mins: f64,
        rx: f64,
        ry: f64,
        rz: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        datetime: DateTime,
    }

    /// Parse the optional Vallado UTC stamp after the state columns
    ///
    /// The stamp follows the seven ephemeris numbers. Fortran spacing can split
    /// the clock time across tokens.
    ///
    /// # Arguments
    /// * `cols` - Whitespace-separated fields from one ephemeris line
    ///
    /// # Returns
    /// * `Some(DateTime)` - Parsed UTC stamp
    /// * `None` - If the line has no stamp or the stamp cannot be parsed
    fn parse_vallado_datetime(cols: &[&str]) -> Option<DateTime> {
        if cols.len() < 11 {
            return None;
        }

        let year = cols[7].parse().ok()?;
        let month = cols[8].parse().ok()?;
        let day = cols[9].parse().ok()?;
        let time_str = cols[10..].join(" ");
        let parts: Vec<&str> = time_str
            .split(':')
            .map(str::trim)
            .filter(|part| !part.is_empty())
            .collect();
        if parts.len() != 3 {
            return None;
        }

        Some(DateTime {
            year,
            month,
            day,
            hour: parts[0].parse().ok()?,
            minute: parts[1].parse().ok()?,
            second: parts[2].parse().ok()?,
            timezone: Timezone::UTC,
        })
    }

    /// Parse a Vallado ephemeris block into rows
    ///
    /// Each line holds minutes from epoch, TEME position in km, and TEME velocity
    /// in km/s. A line without a calendar stamp uses the TLE epoch.
    ///
    /// # Arguments
    /// * `ephem` - Ephemeris text from the test case
    /// * `epoch` - TLE epoch, used when a line has no calendar stamp
    ///
    /// # Returns
    /// * `Vec<EphemRow>` - Parsed rows
    ///
    /// # Panics
    /// * If a line has a calendar stamp that cannot be parsed
    fn parse_vallado_ephem(ephem: &str, epoch: DateTime) -> Vec<EphemRow> {
        ephem
            .lines()
            .filter_map(|line| {
                let cols: Vec<&str> = line.split_whitespace().collect();
                if cols.len() < 7 {
                    return None;
                }
                let datetime = if cols.len() >= 11 {
                    parse_vallado_datetime(&cols).unwrap_or_else(|| {
                        panic!("failed to parse Vallado calendar stamp: {line}");
                    })
                } else {
                    epoch
                };
                Some(EphemRow {
                    t_mins: cols[0].parse().ok()?,
                    rx: cols[1].parse().ok()?,
                    ry: cols[2].parse().ok()?,
                    rz: cols[3].parse().ok()?,
                    vx: cols[4].parse().ok()?,
                    vy: cols[5].parse().ok()?,
                    vz: cols[6].parse().ok()?,
                    datetime,
                })
            })
            .collect()
    }

    /// Shift a UTC stamp so it matches an exact minute offset
    ///
    /// Vallado prints the calendar time with limited precision. The reference
    /// state is at exact `t_mins`, and a fraction of a millisecond is enough to
    /// miss 1 m on a near-decay orbit. The printed stamp must still be within
    /// 1 ms of `t_mins`.
    ///
    /// # Arguments
    /// * `datetime` - Parsed Vallado calendar stamp
    /// * `jd0` - Integer Julian day of the TLE epoch
    /// * `jdfrac0` - Fractional Julian day of the TLE epoch
    /// * `t_mins` - Minutes from epoch for this row
    ///
    /// # Returns
    /// * `DateTime` - Stamp adjusted onto `t_mins`
    ///
    /// # Panics
    /// * If the stamp is more than 1 ms from `t_mins`
    fn align_datetime_to_t_mins(
        datetime: DateTime,
        jd0: f64,
        jdfrac0: f64,
        t_mins: f64,
    ) -> DateTime {
        let (jd, jdfrac) = utc2jday(&datetime).expect("propagation datetime must be UTC");
        let recovered_mins = ((jd - jd0) + (jdfrac - jdfrac0)) * 1440.0;
        let delta_mins = t_mins - recovered_mins;
        assert!(
            delta_mins.abs() * 60.0 < 1.0e-3,
            "Vallado UTC stamp is more than 1 ms from t_mins={t_mins} (recovered {recovered_mins})"
        );
        DateTime {
            second: datetime.second + delta_mins * 60.0,
            ..datetime
        }
    }

    // Agreement with Vallado reference ephemerides to 1e-6 km (1 mm) and 1e-6 km/s (1 mm/s).
    const VALLADO_STATE_TOL_KM: f64 = 1e-6;

    /// Compare a propagated state with one Vallado ephemeris row
    ///
    /// Position and velocity must agree to within `VALLADO_STATE_TOL_KM`.
    ///
    /// # Arguments
    /// * `key` - TOML case key, included in the failure message
    /// * `name` - Case description, included in the failure message
    /// * `t_mins` - Minutes from epoch for this row
    /// * `got` - State returned by mako-sgp4
    /// * `row` - Reference ephemeris row
    ///
    /// # Panics
    /// * If any component differs by at least the Vallado tolerance
    fn assert_state_near(key: &str, name: &str, t_mins: f64, got: &StateVector, row: &EphemRow) {
        let checks = [
            ("Position x", got.r_x, row.rx),
            ("Position y", got.r_y, row.ry),
            ("Position z", got.r_z, row.rz),
            ("Velocity x", got.v_x, row.vx),
            ("Velocity y", got.v_y, row.vy),
            ("Velocity z", got.v_z, row.vz),
        ];
        for (label, value, expected) in checks {
            assert!(
                (value - expected).abs() < VALLADO_STATE_TOL_KM,
                "{key}: {name}\n- {label} mismatch {value} vs {expected} (t_mins = {t_mins})"
            );
        }
    }

    /// Build an initialized SGP4 model from the ISS test TLE
    ///
    /// # Returns
    /// * `Sgp4` - Propagator for the ISS element set
    ///
    /// # Panics
    /// * If the TLE does not parse
    fn iss_sgp4() -> Sgp4 {
        from_tle_string(
            "\
1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
",
        )
        .expect("ISS TLE should parse")
        .into_iter()
        .next()
        .expect("ISS TLE should yield one propagator")
    }

    /// Build the J2000 epoch as a UTC datetime
    ///
    /// # Returns
    /// * `DateTime` - 2000-01-01 12:00:00 UTC
    fn utc_j2000() -> DateTime {
        DateTime {
            year: 2000,
            month: 1,
            day: 1,
            hour: 12,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        }
    }

    /// Reject a non-UTC or pre-Gregorian propagation time
    ///
    /// Checks both initialization and `sgp4_prop_datetime`.
    ///
    /// # Panics
    /// * If either path accepts an invalid datetime
    #[test]
    fn test_invalid_date_time() {
        let mut gp = GenPerturbElementSet {
            epoch_datetime: DateTime {
                timezone: Timezone::UT1,
                ..utc_j2000()
            },
            mean_motion: 15.0,
            ..GenPerturbElementSet::default()
        };
        assert!(matches!(
            init_sgp4(&gp, None),
            Err(Sgp4Error::InvalidDateTime(DateError::DateNotUTC))
        ));

        gp.epoch_datetime = DateTime {
            year: 1500,
            ..utc_j2000()
        };
        assert!(matches!(
            init_sgp4(&gp, None),
            Err(Sgp4Error::InvalidDateTime(DateError::DateTooEarly))
        ));

        let sgp4 = iss_sgp4();
        let mut datetime = sgp4.gp.epoch_datetime;
        datetime.timezone = Timezone::UT1;
        assert!(matches!(
            sgp4_prop_datetime(&sgp4, &datetime),
            Err(Sgp4Error::InvalidDateTime(DateError::DateNotUTC))
        ));

        datetime = DateTime {
            year: 1500,
            ..utc_j2000()
        };
        assert!(matches!(
            sgp4_prop_datetime(&sgp4, &datetime),
            Err(Sgp4Error::InvalidDateTime(DateError::DateTooEarly))
        ));
    }

    /// Reject a mean motion that is not positive
    ///
    /// Checks both initialization and `sgp4_prop_delta`.
    ///
    /// # Panics
    /// * If either path accepts a non-positive mean motion
    #[test]
    fn test_invalid_mean_motion() {
        let gp = GenPerturbElementSet {
            epoch_datetime: utc_j2000(),
            mean_motion: 0.0,
            ..GenPerturbElementSet::default()
        };
        assert!(matches!(
            init_sgp4(&gp, None),
            Err(Sgp4Error::InvalidMeanMotion)
        ));

        let mut sgp4 = iss_sgp4();
        sgp4.brouwer0.n = 0.0;
        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::InvalidMeanMotion)
        ));
    }

    /// Reject a mean eccentricity outside 0 to 1
    ///
    /// # Panics
    /// * If propagation accepts an eccentricity of 1.5 or -0.5
    #[test]
    fn test_invalid_mean_eccentricity() {
        let mut sgp4 = iss_sgp4();
        sgp4.brouwer0.e = 1.5;
        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::InvalidMeanEccentricity)
        ));

        sgp4.brouwer0.e = -0.5;
        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::InvalidMeanEccentricity)
        ));
    }

    /// Reject a deep-space perturbed eccentricity outside 0 to 1
    ///
    /// # Panics
    /// * If propagation accepts the constructed lunar terms
    #[test]
    fn test_invalid_perturbed_eccentricity() {
        let mut sgp4 = iss_sgp4();
        sgp4.deep_space = true;
        sgp4.lunar_params.c = 1.0e6;
        sgp4.lunar_params.x1 = 1.0;
        sgp4.lunar_params.x2 = 1.0;
        sgp4.lunar_params.x3 = 1.0;
        sgp4.lunar_params.x4 = 1.0;

        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::InvalidPerturbedEccentricity)
        ));
    }

    /// Reject a semilatus rectum that is not positive
    ///
    /// # Panics
    /// * If propagation accepts an eccentricity of 0.9999 on the ISS model
    #[test]
    fn test_invalid_semilatus_rectum() {
        let mut sgp4 = iss_sgp4();
        sgp4.brouwer0.e = 0.9999;

        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::InvalidSemilatusRectum)
        ));
    }

    /// Reject a decayed satellite
    ///
    /// # Panics
    /// * If propagation accepts a mean motion of 1 rad/min on the ISS model
    #[test]
    fn test_satellite_decayed() {
        let mut sgp4 = iss_sgp4();
        sgp4.brouwer0.n = 1.0;

        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::SatelliteDecayed)
        ));
    }

    /// Non-finite states are returned as errors
    ///
    /// # Panics
    /// * If a NaN state is returned as Ok
    #[test]
    fn test_non_finite_state() {
        let mut sgp4 = iss_sgp4();
        sgp4.brouwer0.i = f64::NAN;

        assert!(matches!(
            sgp4_prop_delta(&sgp4, 0.0),
            Err(Sgp4Error::NonFiniteState)
        ));
    }

    /// Compare propagation with the Vallado reference ephemerides
    ///
    /// Each non-exception TLE is propagated with both minutes from epoch and
    /// the printed UTC stamp. Exception cases must fail to parse.
    ///
    /// # Panics
    /// * If a case fails to parse, fails to propagate, or misses the tolerance
    #[test]
    fn test_sgp4_vallado_cases() {
        let content = std::fs::read_to_string("test/vallado_cases.toml")
            .expect("could not read test/vallado_cases.toml");
        let cases: ValladoCases =
            from_str(&content).expect("could not parse test/vallado_cases.toml");

        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];

            if case.exception {
                assert!(
                    from_tle_string(&case.tle).is_err(),
                    "case {key}: expected initialization to return Err"
                );
                continue;
            }

            let sgp4s = from_tle_string(&case.tle).expect("Vallado TLE should parse");
            assert!(!sgp4s.is_empty(), "case {key}: TLE failed to parse");
            let sgp4 = &sgp4s[0];

            for row in parse_vallado_ephem(&case.ephem, sgp4.gp.epoch_datetime) {
                // Propagate using the delta time
                let state_delta = sgp4_prop_delta(sgp4, row.t_mins).unwrap_or_else(|err| {
                    panic!("{key}: {} at t_mins={}: {err:?}", case.name, row.t_mins)
                });
                assert_state_near(key, &case.name, row.t_mins, &state_delta, &row);

                // Propagate using the datetime. Align the printed stamp to t_mins
                // so the datetime path is checked at the same instant as the
                // reference ephemeris
                let datetime =
                    align_datetime_to_t_mins(row.datetime, sgp4.jd0, sgp4.jdfrac0, row.t_mins);
                let state_datetime = sgp4_prop_datetime(sgp4, &datetime).unwrap_or_else(|err| {
                    panic!(
                        "{key}: {} datetime at t_mins={}: {err:?}",
                        case.name, row.t_mins
                    )
                });
                assert_state_near(key, &case.name, row.t_mins, &state_datetime, &row);
            }
        }
    }

    /// Parse compact, split, and missing Vallado calendar stamps
    ///
    /// # Panics
    /// * If a stamp is parsed as the wrong calendar time
    #[test]
    fn test_parse_vallado_datetime() {
        let epoch = DateTime {
            year: 2000,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };
        let compact = "     360.00000000   -7154.03120202   -3783.17682504   -3536.19412294  4.741887409 -4.151817765 -2.093935425    2000  6 28  0:50:19.733571";
        let spaced_seconds = "0.0 8827.0 -41223.0 3.63 3.00 0.64 0.00 2004 2 9 10:20: 1.494254";
        let spaced_hms =
            "1844040.0 -31652.0 -66335.0 12774.0 1.71 1.91 -0.60 2009 7 2 9: 0: 0.000282";
        let no_stamp = "       0.00000000    7022.46529266   -1400.08296755       0.03995155  1.893841015  6.405893759  4.534807250";

        let compact_row = &parse_vallado_ephem(compact, epoch)[0];
        let compact_dt = compact_row.datetime;
        assert_eq!(compact_dt.year, 2000);
        assert_eq!(compact_dt.month, 6);
        assert_eq!(compact_dt.day, 28);
        assert_eq!(compact_dt.hour, 0);
        assert_eq!(compact_dt.minute, 50);
        assert!((compact_dt.second - 19.733571).abs() < 1e-9);
        assert_eq!(compact_dt.timezone, Timezone::UTC);

        let spaced_row = &parse_vallado_ephem(spaced_seconds, epoch)[0];
        let spaced_dt = spaced_row.datetime;
        assert_eq!(spaced_dt.hour, 10);
        assert_eq!(spaced_dt.minute, 20);
        assert!((spaced_dt.second - 1.494254).abs() < 1e-9);

        let hms_row = &parse_vallado_ephem(spaced_hms, epoch)[0];
        let hms_dt = hms_row.datetime;
        assert_eq!(hms_dt.year, 2009);
        assert_eq!(hms_dt.hour, 9);
        assert_eq!(hms_dt.minute, 0);
        assert!((hms_dt.second - 0.000282).abs() < 1e-12);

        let epoch_row = &parse_vallado_ephem(no_stamp, epoch)[0];
        assert_eq!(epoch_row.datetime, epoch);
    }

    // -------------------------------------------------------
    // Structs for deserializing python-sgp4 test cases
    // -------------------------------------------------------

    // Root of test/python-sgp4_cases.toml. Each entry is one OMM and its states.
    #[derive(Deserialize)]
    struct PythonSgp4Cases {
        test: HashMap<String, PythonSgp4Case>,
    }

    // One KVN OMM and parallel sample arrays.
    // time_delta is minutes from epoch. time_utc is the matching UTC timestamp.
    // Position is TEME meters. Velocity is TEME meters per second.
    #[derive(Deserialize)]
    struct PythonSgp4Case {
        name: String,
        omm: String,
        time_delta: Vec<f64>,
        time_utc: Vec<String>,
        x_teme_m: Vec<f64>,
        y_teme_m: Vec<f64>,
        z_teme_m: Vec<f64>,
        vx_teme_m_per_s: Vec<f64>,
        vy_teme_m_per_s: Vec<f64>,
        vz_teme_m_per_s: Vec<f64>,
    }

    // A position component passes when it differs by less than 1e-6 km (1 millimeter).
    const PYTHON_SGP4_POS_TOL_M: f64 = 1e-3;

    // A velocity component passes when it differs by less than 1e-6 km/s (1 millimeter per second).
    const PYTHON_SGP4_VEL_TOL_M_S: f64 = 1e-3;

    /// Parse an ISO-8601 UTC timestamp with a fractional second
    ///
    /// The calendar date and the clock time are separated by `T`. Fractional
    /// seconds stay in the `second` field.
    ///
    /// # Arguments
    /// * `text` - Timestamp text from a python-sgp4 case
    ///
    /// # Returns
    /// * `DateTime` - Parsed UTC timestamp
    ///
    /// # Panics
    /// * If the timestamp is missing a date, a time, or a numeric field
    fn parse_iso_utc(text: &str) -> DateTime {
        // Split the calendar date from the clock time
        let (date, time) = text
            .split_once('T')
            .unwrap_or_else(|| panic!("UTC timestamp missing T: {text}"));
        let date_parts: Vec<&str> = date.split('-').collect();
        let time_parts: Vec<&str> = time.split(':').collect();
        assert!(
            date_parts.len() == 3 && time_parts.len() == 3,
            "UTC timestamp has the wrong shape: {text}"
        );

        // Build the UTC datetime. Fractional seconds stay in `second`.
        DateTime {
            year: date_parts[0].parse().expect("year"),
            month: date_parts[1].parse().expect("month"),
            day: date_parts[2].parse().expect("day"),
            hour: time_parts[0].parse().expect("hour"),
            minute: time_parts[1].parse().expect("minute"),
            second: time_parts[2].parse().expect("second"),
            timezone: Timezone::UTC,
        }
    }

    /// Compare one propagated state with a python-sgp4 reference sample
    ///
    /// mako-sgp4 reports kilometers, so each component is scaled to meters
    /// before the tolerance check.
    ///
    /// # Arguments
    /// * `key` - TOML case key, included in the failure message
    /// * `name` - Case description, included in the failure message
    /// * `path` - Which propagator entry point was used, `delta` or `utc`
    /// * `t_mins` - Minutes from epoch for this sample
    /// * `got` - State returned by mako-sgp4
    /// * `expected` - Reference position in meters and velocity in meters per second
    ///
    /// # Panics
    /// * If any component differs by at least the python-sgp4 tolerance
    fn assert_python_sgp4_state(
        key: &str,
        name: &str,
        path: &str,
        t_mins: f64,
        got: &StateVector,
        expected: &[f64; 6],
    ) {
        // Scale mako kilometers into the reference meter units
        let checks = [
            (
                "Position x",
                got.r_x * 1000.0,
                expected[0],
                PYTHON_SGP4_POS_TOL_M,
            ),
            (
                "Position y",
                got.r_y * 1000.0,
                expected[1],
                PYTHON_SGP4_POS_TOL_M,
            ),
            (
                "Position z",
                got.r_z * 1000.0,
                expected[2],
                PYTHON_SGP4_POS_TOL_M,
            ),
            (
                "Velocity x",
                got.v_x * 1000.0,
                expected[3],
                PYTHON_SGP4_VEL_TOL_M_S,
            ),
            (
                "Velocity y",
                got.v_y * 1000.0,
                expected[4],
                PYTHON_SGP4_VEL_TOL_M_S,
            ),
            (
                "Velocity z",
                got.v_z * 1000.0,
                expected[5],
                PYTHON_SGP4_VEL_TOL_M_S,
            ),
        ];
        for (label, value, reference, tol) in checks {
            assert!(
                (value - reference).abs() < tol,
                "{key}: {name} {path}\n- {label} mismatch {value} vs {reference} (t_mins = {t_mins})"
            );
        }
    }

    /// Compare OMM propagation with the python-sgp4 reference states
    ///
    /// Each case is propagated with both minutes from epoch and the stored
    /// UTC timestamp.
    ///
    /// # Panics
    /// * If an OMM fails to parse, a sample fails to propagate, or a state misses the tolerance
    #[test]
    fn test_sgp4_python_sgp4_cases() {
        // Load and parse the reference file
        let content = std::fs::read_to_string("test/python-sgp4_cases.toml")
            .expect("could not read test/python-sgp4_cases.toml");
        let cases: PythonSgp4Cases =
            from_str(&content).expect("could not parse test/python-sgp4_cases.toml");

        // Stable order so a failure names the same case on every run
        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];
            let n = case.time_delta.len();

            // Every sample array must describe the same set of times
            assert!(n > 0, "{key}: no propagation samples");
            assert_eq!(case.time_utc.len(), n, "{key}: time_utc length");
            assert_eq!(case.x_teme_m.len(), n, "{key}: x_teme_m length");
            assert_eq!(case.y_teme_m.len(), n, "{key}: y_teme_m length");
            assert_eq!(case.z_teme_m.len(), n, "{key}: z_teme_m length");
            assert_eq!(case.vx_teme_m_per_s.len(), n, "{key}: vx length");
            assert_eq!(case.vy_teme_m_per_s.len(), n, "{key}: vy length");
            assert_eq!(case.vz_teme_m_per_s.len(), n, "{key}: vz length");

            // One OMM record initializes one propagator
            let sgp4s = from_omm_kvn_string(&case.omm)
                .unwrap_or_else(|err| panic!("{key}: {} failed to parse: {err:?}", case.name));
            assert_eq!(sgp4s.len(), 1, "{key}: expected one OMM");
            let sgp4 = &sgp4s[0];

            for i in 0..n {
                let expected = [
                    case.x_teme_m[i],
                    case.y_teme_m[i],
                    case.z_teme_m[i],
                    case.vx_teme_m_per_s[i],
                    case.vy_teme_m_per_s[i],
                    case.vz_teme_m_per_s[i],
                ];
                let t_mins = case.time_delta[i];

                // Minutes from epoch
                let state_delta = sgp4_prop_delta(sgp4, t_mins).unwrap_or_else(|err| {
                    panic!("{key}: {} delta at t_mins={t_mins}: {err:?}", case.name)
                });
                assert_python_sgp4_state(key, &case.name, "delta", t_mins, &state_delta, &expected);

                // Same instant as a UTC timestamp
                let datetime = parse_iso_utc(&case.time_utc[i]);
                let state_utc = sgp4_prop_datetime(sgp4, &datetime).unwrap_or_else(|err| {
                    panic!("{key}: {} utc at t_mins={t_mins}: {err:?}", case.name)
                });
                assert_python_sgp4_state(key, &case.name, "utc", t_mins, &state_utc, &expected);
            }
        }
    }
}
