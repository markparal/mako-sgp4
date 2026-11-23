//! Module for fitting a general perturbations element set to timestamped
//! states and positions.
//!
//! The public types describe the fit's input and output: observations
//! ([`FitObservation`]) with optional measurement noise, options
//! ([`FitOptions`]), and the result ([`FitResult`]) or error ([`FitError`]).
//! Each observation pairs a time with a [`FitMeasurement`], so one fit can
//! mix kinds of measurement. Observations are built with
//! [`new_fit_observation`] and options from [`FitOptions::default`] with
//! fields assigned, because the fit structs are `#[non_exhaustive]`.
//!
//! The fit works in equinoctial elements, which stay well defined for
//! circular and equatorial orbits where the classical argument of perigee and
//! node are not. This module also holds the conversions between equinoctial
//! elements, the classical elements a TLE or OMM stores, and TEME states, the
//! time helpers that place observations relative to the fit epoch, and the
//! validation of the input.

// ------------------
// External Libraries
// ------------------
use alloc::boxed::Box;
use alloc::string::String;
use alloc::vec::Vec;
use core::f64::consts::PI;
use core::fmt;

// ------------------
// Internal Libraries
// ------------------
use crate::common::{CoordinateFrame, Wgs};
use crate::gp::GenPerturbElementSet;
use crate::linalg::{householder_qr, jacobi_svd, solve_lower};
use crate::math::{abs, atan2, cos, rem_euclid, sin, sqrt};
use crate::sgp4::{Sgp4, Sgp4Error};
use crate::time::{DateTime, Timezone, round_to_tle_epoch, utc2jday};

// -------
// Structs
// -------

/// One timestamped observation for the GP fit
///
/// A time and a [`FitMeasurement`]. Every kind of measurement goes through
/// this one type, so kinds can be mixed freely in one list, in any order,
/// with repeated times. Built with [`new_fit_observation`].
///
/// # Examples
/// ```rust
/// use mako_sgp4::{CoordinateFrame, DateTime, FitMeasurement, Timezone, new_fit_observation};
///
/// // A TEME state at 2024-03-01 12:00:00 UTC, with no measurement noise
/// let datetime = DateTime {
///     year: 2024,
///     month: 3,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UTC,
/// };
/// let measurement = FitMeasurement::State {
///     frame: CoordinateFrame::TEME,
///     position: [6778.0, 0.0, 0.0],
///     velocity: [0.0, 4.8, 5.9],
///     noise: None,
/// };
/// let obs = new_fit_observation(datetime, measurement);
///
/// assert_eq!(obs.datetime, datetime);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct FitObservation {
    /// Observation time; UTC is the only time scale today
    pub datetime: DateTime,

    /// What was measured at that time
    pub measurement: FitMeasurement,
}

/// Identifying fields for the fitted element set
///
/// Each field that is set overrides the same field of the initial GP set, if
/// one is given; each field that is `None` comes from the initial GP set, or
/// else from the default: empty name and designator, catalog number 0,
/// classification `U`, and element set and rev numbers 0. The fit doesn't
/// check these against any output format; the TLE and OMM writers do.
///
/// # Examples
/// ```rust
/// use mako_sgp4::GpMetadata;
///
/// // Name the result and give it a catalog number; leave the rest unset
/// let mut metadata = GpMetadata::default();
/// metadata.name = Some("ISS (ZARYA)".into());
/// metadata.catalog_number = Some(25544);
///
/// assert_eq!(metadata.classification, None);
/// ```
#[derive(Debug, Default, Clone, PartialEq)]
#[non_exhaustive]
pub struct GpMetadata {
    /// Object name
    pub name: Option<String>,

    /// Satellite catalog number
    pub catalog_number: Option<i32>,

    /// International designator
    pub intl_designator: Option<String>,

    /// Classification, such as `U` for unclassified
    pub classification: Option<char>,

    /// Element set number; never incremented by the fit
    pub element_set_number: Option<i32>,

    /// Revolution number at the fit epoch, used as given
    pub rev_number: Option<i64>,
}

/// Options for the GP fit
///
/// Built from [`FitOptions::default`] with fields assigned. Every value is
/// checked before the fit starts, and the first one out of range is named in
/// [`FitError::InvalidOptions`]; the metadata is not checked.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitOptions;
///
/// // Center the B* prior on a ballistic-coefficient estimate and name the result
/// let mut opts = FitOptions::default();
/// opts.bstar_prior_center = Some(2.5e-4);
/// opts.metadata.catalog_number = Some(25544);
///
/// assert_eq!(opts.outlier_cutoff, Some(10.0));
/// ```
#[derive(Clone)]
#[non_exhaustive]
pub struct FitOptions {
    /// Fit epoch; `None` uses the latest observation time. Either way it is
    /// rounded to the nearest 1e-8 day, the TLE epoch resolution
    pub epoch: Option<DateTime>,

    /// Element set that gives the starting orbit, ndot and nddot, the prior
    /// center unless `bstar_prior_center` is set, and any metadata not set in
    /// `metadata`
    pub initial_guess: Option<GenPerturbElementSet>,

    /// Identifying fields for the result
    pub metadata: GpMetadata,

    /// Center of the B* prior \[1/ER\]; `None` uses the initial GP set's B*,
    /// or else 0. Must be finite if set
    pub bstar_prior_center: Option<f64>,

    /// One-sigma width of the B* prior \[1/ER\], finite and positive; default 1e-2
    pub bstar_prior_width: f64,

    /// Outlier cutoff as a multiple of the median normalized residual,
    /// finite and at least 3; default `Some(10.0)`, and `None` turns
    /// rejection off
    pub outlier_cutoff: Option<f64>,

    /// Largest fraction of the input observations that may be rejected,
    /// above 0 and below 0.5; default 0.25
    pub max_rejected_fraction: f64,

    /// Jacobian builds per Levenberg-Marquardt run, at least 1
    pub max_iterations: usize,

    /// Convergence threshold on the Gauss-Newton step, in units of its own
    /// one-sigma size, finite and positive
    pub tolerance: f64,

    /// Starting Levenberg-Marquardt damping, finite, positive and below
    /// `max_damping`; default 1e-3
    pub initial_damping: f64,

    /// Damping above which the fit is stuck, finite and positive, and must
    /// exceed `initial_damping`; default 1e10
    pub max_damping: f64,

    /// Earth model; `None` uses WGS-72
    pub wgs: Option<Wgs>,
}

/// Result of a GP fit
///
/// The fitted element set ready to propagate, its residuals against every
/// input observation, and its covariance. Residuals and RMS values come from
/// the returned `sgp4`, so they describe exactly what is propagated.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitResult;
///
/// // Summarize a fit: position RMS in meters and the number kept
/// fn summarize(result: &FitResult) -> (f64, usize) {
///     let kept = result.residuals.len() - result.rejected.len();
///     (result.rms_position * 1000.0, kept)
/// }
/// ```
#[derive(Clone)]
#[non_exhaustive]
pub struct FitResult {
    /// The fitted element set, initialized and ready to propagate or write out
    pub sgp4: Sgp4,

    /// Unweighted RMS of the 3D position residuals over kept observations \[km\]
    pub rms_position: f64,

    /// Unweighted RMS of the 3D velocity residuals over kept observations
    /// \[km/s\]; `None` if no kept observation has a velocity
    pub rms_velocity: Option<f64>,

    /// One residual per input observation, in input order, including
    /// rejected ones
    pub residuals: Vec<FitResidual>,

    /// Input indices of the rejected observations, ascending
    pub rejected: Vec<usize>,

    /// Covariance of the fitted unknowns, with the data's information scaled
    /// by `variance_factor` and the B* prior's left as given (unscaled for a
    /// single state). Order and units:
    /// n (Kozai) \[rad/min\], a_f, a_g, chi, psi \[\], lambda \[rad\], B* \[1/ER\]
    pub covariance: [[f64; 7]; 7],

    /// Retrograde factor, +1 or -1, which defines a_f, a_g, chi and psi in
    /// the covariance
    pub retrograde_factor: f64,

    /// Variance factor of the fit; `None` with no degrees of freedom (a
    /// single state)
    pub variance_factor: Option<f64>,

    /// Jacobian builds, summed over all outlier passes
    pub iterations: usize,
}

/// Residual of one observation against the fitted element set
///
/// Observed minus computed, from the returned [`Sgp4`]. The RIC axes come
/// from the fitted state at the observation time: radial along r, cross-track
/// along r x v, and in-track completing the frame.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitResidual;
///
/// // Radial position error of one observation, in meters
/// fn radial_meters(residual: &FitResidual) -> f64 {
///     residual.position_ric[0] * 1000.0
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct FitResidual {
    /// Position residual in TEME \[km\]
    pub position: [f64; 3],

    /// Velocity residual in TEME \[km/s\]; `None` for a position-only observation
    pub velocity: Option<[f64; 3]>,

    /// Position residual in radial, in-track, cross-track axes \[km\]
    pub position_ric: [f64; 3],

    /// TEME velocity residual rotated into the RIC axes \[km/s\]; not the rate
    /// of change of the RIC position error. `None` for a position-only
    /// observation
    pub velocity_ric: Option<[f64; 3]>,
}

/// A validated observation, ready for the noise model
///
/// Built by [`validate_observation`]: the time as a split Julian date, the
/// measurement in TEME, and a factor of the measurement noise (all zero
/// without noise).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CheckedObservation {
    /// Time as (Julian day, day fraction) from [`crate::utc2jday`]
    pub time: (f64, f64),

    /// Measurement and the factor of its measurement noise
    pub measurement: CheckedMeasurement,
}

/// Equinoctial orbital elements
///
/// The six orbital unknowns of the fit, in the form of equation 1 of the fit
/// plan. With retrograde factor f_r (+1 for inclination up to 90 degrees, -1
/// above), they relate to the classical elements by
///
/// - a_f = e*cos(w + f_r*O), a_g = e*sin(w + f_r*O)
/// - chi = tan(i/2)^f_r * sin(O), psi = tan(i/2)^f_r * cos(O)
/// - lambda = M + w + f_r*O
///
/// where w is the argument of perigee and O the right ascension of the
/// ascending node. For mean elements `n` is the Kozai mean motion; for an
/// osculating state it is the two-body mean motion.
///
/// # References
/// - D. A. Danielson et al., *Semianalytic Satellite Theory*, NPS-MA-95-002, 1995, Section 2.1
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EquinoctialElements {
    /// Mean motion \[rad/min\]
    pub n: f64,

    /// Eccentricity vector component along the equinoctial f axis \[\]
    pub a_f: f64,

    /// Eccentricity vector component along the equinoctial g axis \[\]
    pub a_g: f64,

    /// Node vector component, tan(i/2)^f_r * sin(O) \[\]
    pub chi: f64,

    /// Node vector component, tan(i/2)^f_r * cos(O) \[\]
    pub psi: f64,

    /// Mean longitude, M + w + f_r*O \[rad\]
    pub lambda: f64,
}

/// Classical orbital elements in the units of a GP element set
///
/// Angles are in degrees in [0, 360) and the mean motion in revolutions per
/// day, as [`crate::GenPerturbElementSet`] stores them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ClassicalElements {
    /// Mean motion \[revs/day\]
    pub mean_motion: f64,

    /// Eccentricity \[\]
    pub eccentricity: f64,

    /// Inclination \[degrees\]
    pub inclination: f64,

    /// Right ascension of the ascending node \[degrees\]
    pub right_ascension_of_ascending_node: f64,

    /// Argument of perigee \[degrees\]
    pub argument_of_perigee: f64,

    /// Mean anomaly \[degrees\]
    pub mean_anomaly: f64,
}

// -----
// Enums
// -----

/// What an observation measured
///
/// A full state or a position only, each in a stated frame and with optional
/// measurement noise. The noise is added to the fit's internal model-error
/// term, so `None` (or a sigma of 0) still gives the observation a
/// well-defined weight. More kinds, such as angles or range, can be added as
/// new variants. The variants are built directly, as struct literals.
///
/// # Examples
/// ```rust
/// use mako_sgp4::{CoordinateFrame, FitMeasurement, FitPositionNoise, FitStateNoise};
///
/// // A GPS state with 10 m and 1 cm/s measurement noise
/// let state = FitMeasurement::State {
///     frame: CoordinateFrame::TEME,
///     position: [6778.0, 0.0, 0.0],
///     velocity: [0.0, 4.8, 5.9],
///     noise: Some(FitStateNoise::Sigmas {
///         position: 0.010,
///         velocity: 1e-5,
///     }),
/// };
///
/// // A radar position with 100 m isotropic measurement noise
/// let position = FitMeasurement::Position {
///     frame: CoordinateFrame::TEME,
///     position: [6778.0, 0.0, 0.0],
///     noise: Some(FitPositionNoise::Sigma(0.1)),
/// };
///
/// assert_ne!(state, position);
/// ```
// A covariance variant is about 0.3 kB against 16 bytes for the sigma forms;
// boxing it would cost an allocation per observation for no real saving
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FitMeasurement {
    /// Position and velocity
    State {
        /// Frame of the position and velocity; TEME is the only frame today
        frame: CoordinateFrame,

        /// Position \[km\]
        position: [f64; 3],

        /// Velocity \[km/s\]
        velocity: [f64; 3],

        /// Measurement noise; `None` uses the model-error term only
        noise: Option<FitStateNoise>,
    },

    /// Position only, as from a sensor such as radar
    Position {
        /// Frame of the position; TEME is the only frame today
        frame: CoordinateFrame,

        /// Position \[km\]
        position: [f64; 3],

        /// Measurement noise; `None` uses the model-error term only
        noise: Option<FitPositionNoise>,
    },
}

/// Measurement noise of a full state
///
/// Added to the fit's model-error term, never replacing it, so a sigma of 0
/// or a singular covariance is accepted. A covariance need only be symmetric
/// to roundoff, so one rotated from RIC sigmas into TEME passes.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitStateNoise;
///
/// // 10 m position noise, 1 cm/s velocity noise, no correlation
/// let noise = FitStateNoise::Sigmas {
///     position: 0.010,
///     velocity: 1e-5,
/// };
///
/// assert!(matches!(noise, FitStateNoise::Sigmas { .. }));
/// ```
// A covariance variant is about 0.3 kB against 16 bytes for the sigma forms;
// boxing it would cost an allocation per observation for no real saving
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FitStateNoise {
    /// Isotropic one-sigma noise, each finite and at least 0
    Sigmas {
        /// Position sigma \[km\]
        position: f64,

        /// Velocity sigma \[km/s\]
        velocity: f64,
    },

    /// Full covariance of position then velocity, including position-velocity
    /// correlation \[km^2, km^2/s, km^2/s^2\]; symmetric and positive semidefinite
    Covariance([[f64; 6]; 6]),
}

/// Measurement noise of a position
///
/// Added to the fit's model-error term, never replacing it, so a sigma of 0
/// or a singular covariance is accepted.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitPositionNoise;
///
/// // 10 m range noise along x, 5 km cross-range noise along y and z
/// let noise = FitPositionNoise::Covariance([
///     [1e-4, 0.0, 0.0],
///     [0.0, 25.0, 0.0],
///     [0.0, 0.0, 25.0],
/// ]);
///
/// assert!(matches!(noise, FitPositionNoise::Covariance(_)));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FitPositionNoise {
    /// Isotropic one-sigma noise \[km\], finite and at least 0
    Sigma(f64),

    /// Full covariance \[km^2\]; symmetric and positive semidefinite
    Covariance([[f64; 3]; 3]),
}

/// GP fit errors
///
/// Input errors are found before any propagation, in the order options,
/// observations, then the set as a whole.
///
/// # Examples
/// ```rust
/// use mako_sgp4::FitError;
///
/// // Name the failed option in a message
/// let err = FitError::InvalidOptions {
///     field: "bstar_prior_width",
/// };
///
/// assert!(err.to_string().contains("bstar_prior_width"));
/// ```
#[derive(Clone)]
#[non_exhaustive]
pub enum FitError {
    /// The observation at this index has a non-finite value, a time that is
    /// not a valid date, a negative or non-finite noise sigma, or a noise
    /// covariance that is not finite, not symmetric, or not positive
    /// semidefinite
    InvalidObservation {
        /// Index of the first invalid observation
        index: usize,
    },

    /// An option is out of range, or the epoch is not a valid date
    InvalidOptions {
        /// Name of the first invalid option, such as "epoch"
        field: &'static str,
    },

    /// No full state and fewer than 3 distinct position times. Outlier
    /// rejection never causes it
    TooFewObservations,

    /// No starting guess could be built: without an initial GP set, no
    /// candidate could be built from the data or SGP4 failed on every one;
    /// with one, the fixed-point loop at the epoch couldn't evaluate any
    /// iterate
    StartingGuessFailed,

    /// SGP4 rejected the initial GP set or the epoch
    Sgp4(Sgp4Error),

    /// The observations don't determine every unknown, so no covariance and
    /// no solution can be given
    Underdetermined,

    /// The fit hit the iteration cap or the damping ceiling, or SGP4 couldn't
    /// evaluate a Jacobian column on either side. Holds the best solution,
    /// with its covariance
    NotConverged(Box<FitResult>),
}

/// A validated measurement and the factor of its measurement noise
///
/// In TEME. The noise is kept as the factor F from [`checked_noise_factor`],
/// with F*F^T the measurement-noise covariance up to its roundoff; F is all
/// zero without measurement noise.
// A state is about 0.3 kB against 0.1 kB for a position; boxing it would
// cost an allocation per observation for no real saving
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum CheckedMeasurement {
    /// Full state
    State {
        /// Position \[km\]
        position: [f64; 3],

        /// Velocity \[km/s\]
        velocity: [f64; 3],

        /// Factor of the measurement noise of position then velocity
        /// \[km, km/s\]
        sensor_factor: [[f64; 6]; 6],
    },

    /// Position only
    Position {
        /// Position \[km\]
        position: [f64; 3],

        /// Factor of the measurement noise \[km\]
        sensor_factor: [[f64; 3]; 3],
    },
}

/// Noise factor of one observation
///
/// Equation 3 of the fit plan: the lower-triangular L_j, with a positive
/// diagonal, of the total noise covariance R_j = L_j*L_j^T. The variant
/// matches the observation's measurement: position then velocity for a full
/// state, position alone for a position.
// A state factor is 288 bytes against 72 for a position; boxing it would
// cost an allocation per observation for no real saving
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum NoiseFactor {
    /// Full state, 6 x 6 \[km, km/s\]
    State([[f64; 6]; 6]),

    /// Position, 3 x 3 \[km\]
    Position([[f64; 3]; 3]),
}

// ------
// Traits
// ------

/// Default fit options
///
/// No epoch override, initial GP set or metadata; a B* prior of width 1e-2
/// centered on the initial GP set's B* or 0; outlier rejection at 10 times the
/// median normalized residual, rejecting at most a quarter of the
/// observations; Levenberg-Marquardt damping from 1e-3 up to 1e10; WGS-72.
impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            epoch: None,
            initial_guess: None,
            metadata: GpMetadata::default(),
            bstar_prior_center: None,
            bstar_prior_width: 1e-2,
            outlier_cutoff: Some(10.0),
            max_rejected_fraction: 0.25,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
            initial_damping: 1e-3,
            max_damping: 1e10,
            wgs: None,
        }
    }
}

/// Debug output for [`FitError`]
///
/// Written by hand because [`FitResult`] holds an [`Sgp4`], which has no
/// `Debug`; `NotConverged` prints without its result.
impl fmt::Debug for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitError::InvalidObservation { index } => f
                .debug_struct("InvalidObservation")
                .field("index", index)
                .finish(),
            FitError::InvalidOptions { field } => f
                .debug_struct("InvalidOptions")
                .field("field", field)
                .finish(),
            FitError::TooFewObservations => f.write_str("TooFewObservations"),
            FitError::StartingGuessFailed => f.write_str("StartingGuessFailed"),
            FitError::Sgp4(err) => f.debug_tuple("Sgp4").field(err).finish(),
            FitError::Underdetermined => f.write_str("Underdetermined"),
            FitError::NotConverged(_) => f.write_str("NotConverged(..)"),
        }
    }
}

/// Human-readable [`FitError`] messages
///
/// Messages start lowercase and have no trailing period, so they compose
/// into longer messages.
impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitError::InvalidObservation { index } => write!(
                f,
                "observation {index} has a non-finite value, an invalid time, or invalid measurement noise"
            ),
            FitError::InvalidOptions { field } => {
                write!(
                    f,
                    "fit option `{field}` is out of range or not a valid date"
                )
            }
            FitError::TooFewObservations => {
                f.write_str("too few observations: need a full state or 3 distinct position times")
            }
            FitError::StartingGuessFailed => f.write_str("no starting guess could be built"),
            FitError::Sgp4(_) => f.write_str("SGP4 rejected the initial GP set or the epoch"),
            FitError::Underdetermined => {
                f.write_str("the observations don't determine every element")
            }
            FitError::NotConverged(_) => f.write_str("the fit did not converge"),
        }
    }
}

/// [`FitError`] is a standard error, usable with `?` into `Box<dyn Error>`
impl core::error::Error for FitError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            FitError::Sgp4(err) => Some(err),
            _ => None,
        }
    }
}

// ---------
// Constants
// ---------

/// Minutes in a day
///
/// Converts day counts to the minutes SGP4 propagates in.
const MINUTES_PER_DAY: f64 = 1440.0;

/// Default Jacobian builds per Levenberg-Marquardt run
///
/// Provisional; set by the LM convergence tests (implementation step 6).
const DEFAULT_MAX_ITERATIONS: usize = 25;

/// Default convergence threshold, in units of the step's own one sigma
///
/// Provisional; set by the LM convergence tests (implementation step 6).
const DEFAULT_TOLERANCE: f64 = 1e-3;

/// Smallest allowed outlier cutoff, as a multiple of the median normalized residual
///
/// Below about 3 the test starts rejecting good observations.
const MIN_OUTLIER_CUTOFF: f64 = 3.0;

/// Relative tolerance on the asymmetry of a noise covariance
///
/// Entry (i, j) may differ from (j, i) by this times d_i*d_j, where d_i is
/// the square root of the mean variance of the 3-vector holding component i.
/// Measuring against each 3-vector's own scale keeps the test fair across
/// the position and velocity blocks of a state covariance, and passes a
/// covariance rotated from RIC sigmas into TEME, whose roundoff is about
/// 1e-16 of the largest variance in its block.
const SYMMETRY_TOLERANCE: f64 = 1e-12;

/// Tolerance of the positive-semidefinite test of a noise covariance
///
/// The covariance, scaled as for [`SYMMETRY_TOLERANCE`], must be reproduced
/// to this in every entry by the factor built from its non-negative
/// eigenvalues, so a singular covariance and one with a roundoff-sized
/// negative eigenvalue pass. Roundoff in a rotated covariance is about 1e-16
/// of its block's variance.
const PSD_TOLERANCE: f64 = 1e-12;

/// SGP4 model error \[km\]
///
/// The one-sigma position error every observation carries from SGP4 itself,
/// added to its measurement noise. An internal constant with no user option;
/// the velocity term is n0 times this.
const MODEL_ERROR_SIGMA: f64 = 1.0;

/// Seconds in a minute
///
/// Converts a mean motion in rad/min to the rad/s that scales a velocity
/// sigma in km/s.
const SECONDS_PER_MINUTE: f64 = 60.0;

// ---------
// Functions
// ---------

/// Build an observation from a time and a measurement
///
/// # Arguments
/// * `datetime` - The observation time
/// * `measurement` - What was measured at that time
///
/// # Returns
/// * `obs` - The observation
///
/// # Examples
/// ```rust
/// use mako_sgp4::{CoordinateFrame, DateTime, FitMeasurement, FitPositionNoise, Timezone, new_fit_observation};
///
/// // A radar position with 100 m noise at 2024-03-01 12:00:00 UTC
/// let datetime = DateTime {
///     year: 2024,
///     month: 3,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UTC,
/// };
/// let position = FitMeasurement::Position {
///     frame: CoordinateFrame::TEME,
///     position: [6778.0, 0.0, 0.0],
///     noise: Some(FitPositionNoise::Sigma(0.1)),
/// };
/// let obs = new_fit_observation(datetime, position);
///
/// assert_eq!(obs.measurement, position);
/// ```
pub fn new_fit_observation(datetime: DateTime, measurement: FitMeasurement) -> FitObservation {
    FitObservation {
        datetime,
        measurement,
    }
}

/// Validate the fit input
///
/// Algorithm step 1 of the fit plan: the options, then each observation in
/// order, then the set as a whole, stopping at the first failure.
///
/// # Arguments
/// * `observations` - The input observations
/// * `opts` - The fit options
///
/// # Returns
/// * `checked` - One checked observation per input observation, in input order
///
/// # Errors
/// * [`FitError::InvalidOptions`] - If an option is out of range
/// * [`FitError::InvalidObservation`] - If an observation is invalid
/// * [`FitError::TooFewObservations`] - If there is no full state and fewer
///   than 3 distinct position times
pub(crate) fn validate_input(
    observations: &[FitObservation],
    opts: &FitOptions,
) -> Result<Vec<CheckedObservation>, FitError> {
    validate_options(opts)?;
    let checked = observations
        .iter()
        .enumerate()
        .map(|(index, obs)| validate_observation(index, obs))
        .collect::<Result<Vec<_>, _>>()?;
    check_observation_count(&checked)?;
    Ok(checked)
}

/// Validate the fit options
///
/// Checks run in the order the fields are listed in [`FitOptions`], so with
/// two invalid fields the first is named. The metadata and the Earth model
/// are not checked.
///
/// # Arguments
/// * `opts` - The fit options
///
/// # Returns
/// * `()` - If every option is in range
///
/// # Errors
/// * [`FitError::InvalidOptions`] - Naming the first invalid field.
///   `initial_damping` at or above `max_damping` names `max_damping`
pub(crate) fn validate_options(opts: &FitOptions) -> Result<(), FitError> {
    let fraction = opts.max_rejected_fraction;
    let checks: [(bool, &'static str); 10] = [
        (!is_positive(opts.bstar_prior_width), "bstar_prior_width"),
        (
            opts.bstar_prior_center.is_some_and(|c| !c.is_finite()),
            "bstar_prior_center",
        ),
        (
            opts.outlier_cutoff
                .is_some_and(|c| !c.is_finite() || c < MIN_OUTLIER_CUTOFF),
            "outlier_cutoff",
        ),
        (
            !fraction.is_finite() || fraction <= 0.0 || fraction >= 0.5,
            "max_rejected_fraction",
        ),
        (opts.max_iterations == 0, "max_iterations"),
        (!is_positive(opts.tolerance), "tolerance"),
        (!is_positive(opts.initial_damping), "initial_damping"),
        (!is_positive(opts.max_damping), "max_damping"),
        (opts.initial_damping >= opts.max_damping, "max_damping"),
        (
            opts.epoch
                .as_ref()
                .is_some_and(|e| julian_date(e).is_none()),
            "epoch",
        ),
    ];
    match checks.iter().find(|(failed, _)| *failed) {
        Some((_, field)) => Err(FitError::InvalidOptions { field }),
        None => Ok(()),
    }
}

/// Validate one observation
///
/// The time must convert to a Julian date, the position and any velocity
/// must be finite, and any measurement noise must pass
/// [`checked_variance`] and [`checked_noise_factor`]; a sigma form is checked
/// and factored as the equivalent diagonal covariance. There is no range check on
/// magnitudes. The measurement is matched on its frame and the time on its
/// time scale exhaustively, so a new frame or scale forces a conversion
/// decision here.
///
/// # Arguments
/// * `index` - The observation's input index, reported on failure
/// * `obs` - The observation
///
/// # Returns
/// * `checked` - The time as a split Julian date and the TEME measurement,
///   with the factor of its measurement noise (all zero without noise)
///
/// # Errors
/// * [`FitError::InvalidObservation`] - With `index`, if any check fails
pub(crate) fn validate_observation(
    index: usize,
    obs: &FitObservation,
) -> Result<CheckedObservation, FitError> {
    let invalid = || FitError::InvalidObservation { index };
    let time = julian_date(&obs.datetime).ok_or_else(invalid)?;
    let measurement = match obs.measurement {
        FitMeasurement::State {
            frame,
            position,
            velocity,
            noise,
        } => {
            let (position, velocity) = match frame {
                CoordinateFrame::TEME => (position, velocity),
            };
            if !all_finite(&position) || !all_finite(&velocity) {
                return Err(invalid());
            }
            let cov = match noise {
                None => [[0.0; 6]; 6],
                Some(FitStateNoise::Sigmas {
                    position: sigma_r,
                    velocity: sigma_v,
                }) => {
                    let var_r = checked_variance(sigma_r).ok_or_else(invalid)?;
                    let var_v = checked_variance(sigma_v).ok_or_else(invalid)?;
                    diagonal(core::array::from_fn(|i| if i < 3 { var_r } else { var_v }))
                }
                Some(FitStateNoise::Covariance(cov)) => cov,
            };
            CheckedMeasurement::State {
                position,
                velocity,
                sensor_factor: checked_noise_factor(&cov).ok_or_else(invalid)?,
            }
        }
        FitMeasurement::Position {
            frame,
            position,
            noise,
        } => {
            let position = match frame {
                CoordinateFrame::TEME => position,
            };
            if !all_finite(&position) {
                return Err(invalid());
            }
            let cov = match noise {
                None => [[0.0; 3]; 3],
                Some(FitPositionNoise::Sigma(sigma)) => {
                    diagonal([checked_variance(sigma).ok_or_else(invalid)?; 3])
                }
                Some(FitPositionNoise::Covariance(cov)) => cov,
            };
            CheckedMeasurement::Position {
                position,
                sensor_factor: checked_noise_factor(&cov).ok_or_else(invalid)?,
            }
        }
    };
    Ok(CheckedObservation { time, measurement })
}

/// Check that the observations can support a fit
///
/// A fit needs at least one full state or 3 distinct position times; with
/// the B* prior row that always gives at least 7 rows. Positions at the same
/// time count once.
///
/// # Arguments
/// * `observations` - The checked observations
///
/// # Returns
/// * `()` - If there are enough observations
///
/// # Errors
/// * [`FitError::TooFewObservations`] - If there is no full state and fewer
///   than 3 distinct position times
pub(crate) fn check_observation_count(observations: &[CheckedObservation]) -> Result<(), FitError> {
    let mut distinct = [(0.0, 0.0); 3];
    let mut count = 0;
    for obs in observations {
        match obs.measurement {
            CheckedMeasurement::State { .. } => return Ok(()),
            CheckedMeasurement::Position { .. } => {
                let seen = distinct[..count]
                    .iter()
                    .any(|t| minutes_since_epoch(*t, obs.time) == 0.0);
                if !seen && count < 3 {
                    distinct[count] = obs.time;
                    count += 1;
                }
            }
        }
    }
    if count == 3 {
        Ok(())
    } else {
        Err(FitError::TooFewObservations)
    }
}

/// Variance of a measurement-noise sigma
///
/// # Arguments
/// * `sigma` - The one-sigma noise
///
/// # Returns
/// * `Some(variance)` - sigma^2, if sigma is finite, at least 0, and its
///   square is finite
/// * `None` - Otherwise
pub(crate) fn checked_variance(sigma: f64) -> Option<f64> {
    let variance = sigma * sigma;
    (sigma.is_finite() && sigma >= 0.0 && variance.is_finite()).then_some(variance)
}

/// Check and factor a measurement-noise covariance
///
/// The matrix holds one or two 3-vectors (a position, or a position then a
/// velocity). Each component i gets the scale d_i, the square root of the
/// mean variance of its 3-vector, so the position and velocity blocks of a
/// state covariance are each judged in their own units. The matrix must be
/// finite, have no negative variance, and be symmetric to
/// [`SYMMETRY_TOLERANCE`] times d_i*d_j. A block with zero scale must be
/// exactly zero, including its correlations, and then counts as no noise; an
/// all-zero matrix is accepted.
///
/// The symmetrized matrix, with each entry scaled by 1/(d_i*d_j), is
/// decomposed by [`jacobi_svd`]; for a symmetric matrix the singular values
/// are the absolute eigenvalues, and an eigenvalue is negative where its left
/// and right singular vectors point opposite ways. The factor keeps the
/// non-negative eigenvalues only, so F*F^T is the covariance projected onto
/// the positive semidefinite matrices. The matrix is positive semidefinite if
/// that projection reproduces it to [`PSD_TOLERANCE`] in every scaled entry.
/// Unlike a diagonal shift, the projection changes no direction by more than
/// the covariance's own roundoff, so an accurate axis beside a huge one keeps
/// its weight.
///
/// # Arguments
/// * `cov` - The covariance \[N x N\], with N = 3 or 6
///
/// # Returns
/// * `Some(f)` - A factor F, with F*F^T the symmetrized covariance up to its
///   roundoff; its rows are zero for a zero-scale block
/// * `None` - If any check fails
///
/// # Panics
/// * If N is not a multiple of 3
pub(crate) fn checked_noise_factor<const N: usize>(cov: &[[f64; N]; N]) -> Option<[[f64; N]; N]> {
    assert!(N % 3 == 0, "covariance size must be a multiple of 3");
    if cov.iter().flatten().any(|x| !x.is_finite()) {
        return None;
    }
    // A variance is never negative, even with roundoff in a rotated matrix
    if (0..N).any(|i| cov[i][i] < 0.0) {
        return None;
    }

    // Scale of each component: the root mean variance of its 3-vector
    let scale: [f64; N] = core::array::from_fn(|i| {
        let b = 3 * (i / 3);
        sqrt(cov[b][b] / 3.0 + cov[b + 1][b + 1] / 3.0 + cov[b + 2][b + 2] / 3.0)
    });

    // Symmetric to roundoff, then symmetrized and scaled
    let mut scaled = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..N {
            let limit = SYMMETRY_TOLERANCE * scale[i] * scale[j];
            if abs(cov[i][j] - cov[j][i]) > limit {
                return None;
            }
            let sym = 0.5 * cov[i][j] + 0.5 * cov[j][i];
            if scale[i] > 0.0 && scale[j] > 0.0 {
                scaled[i][j] = sym / (scale[i] * scale[j]);
            } else if sym != 0.0 {
                // A zero-scale block must have no variance and no correlation
                return None;
            }
        }
    }

    // Keep the non-negative eigenvalues: column k of the scaled factor is
    // sqrt(sigma_k)*v_k where u_k and v_k agree
    let svd = jacobi_svd(&scaled);
    let mut factor = [[0.0; N]; N];
    for k in 0..N {
        let agree: f64 = (0..N).map(|i| svd.u[i][k] * svd.v[i][k]).sum();
        if agree > 0.0 {
            let root = sqrt(svd.sigma[k]);
            for (row, v_row) in factor.iter_mut().zip(&svd.v) {
                row[k] = root * v_row[k];
            }
        }
    }

    // Positive semidefinite if the projection reproduces the matrix
    for i in 0..N {
        for j in 0..N {
            let fft: f64 = (0..N).map(|k| factor[i][k] * factor[j][k]).sum();
            if abs(fft - scaled[i][j]) > PSD_TOLERANCE {
                return None;
            }
        }
    }
    Some(core::array::from_fn(|i| {
        core::array::from_fn(|j| scale[i] * factor[i][j])
    }))
}

/// Noise factor of one observation
///
/// Equation 3 of the fit plan: L_j with L_j*L_j^T = R_j = R_sensor,j +
/// sigma_model^2 * M_j, where M_j = diag(I3, n0^2 * I3) for a state and I3
/// for a position, and R_sensor,j = F*F^T from [`checked_noise_factor`]. The
/// velocity term n0 * sigma_model keeps position and velocity weights
/// consistent in every orbit regime, because a position error sigma on an
/// orbit goes with about n * sigma of velocity error.
///
/// L_j comes from [`stacked_factor`], not from a Cholesky of the sum, which
/// can fail in floating point when the measurement noise is many orders of
/// magnitude above the model error. The model-error term has full rank, so
/// L_j exists for every observation that passed validation, whatever n0 is.
///
/// # Arguments
/// * `measurement` - The checked measurement, with its measurement noise
/// * `n0` - The reference mean motion fixed for the fit \[rad/min\]
///
/// # Returns
/// * `l` - The noise factor L_j
///
/// # Panics
/// * In debug builds, if `n0` is not finite and positive
pub(crate) fn noise_factor(measurement: &CheckedMeasurement, n0: f64) -> NoiseFactor {
    debug_assert!(n0.is_finite() && n0 > 0.0, "n0 must be finite and positive");
    let sigma_v = n0 / SECONDS_PER_MINUTE * MODEL_ERROR_SIGMA;
    match measurement {
        CheckedMeasurement::State { sensor_factor, .. } => {
            let model = core::array::from_fn(|i| if i < 3 { MODEL_ERROR_SIGMA } else { sigma_v });
            NoiseFactor::State(stacked_factor(sensor_factor, &model))
        }
        CheckedMeasurement::Position { sensor_factor, .. } => {
            NoiseFactor::Position(stacked_factor(sensor_factor, &[MODEL_ERROR_SIGMA; 3]))
        }
    }
}

/// Noise factors of every observation
///
/// Builds L_j for each observation once, for the n0 fixed at the start of the
/// fit, so the weights stay the same through every iteration and outlier
/// pass.
///
/// # Arguments
/// * `observations` - The checked observations
/// * `n0` - The reference mean motion fixed for the fit \[rad/min\]
///
/// # Returns
/// * `factors` - L_j for each observation, in input order
///
/// # Panics
/// * In debug builds, if `n0` is not finite and positive
pub(crate) fn noise_factors(observations: &[CheckedObservation], n0: f64) -> Vec<NoiseFactor> {
    observations
        .iter()
        .map(|obs| noise_factor(&obs.measurement, n0))
        .collect()
}

/// Normalized residual of one observation
///
/// Equation 3 of the fit plan: rho_j = |L_j^-1 * y_j| / sqrt(m_j), the
/// whitened residual's norm per measured component, so a full state and a
/// position with the same error in every component score the same.
///
/// # Arguments
/// * `l` - The noise factor L_j \[N x N\]
/// * `y` - The residual y_j, observed minus computed \[N\]
///
/// # Returns
/// * `rho` - The normalized residual
pub(crate) fn normalized_residual<const N: usize>(l: &[[f64; N]; N], y: &[f64; N]) -> f64 {
    let w = solve_lower(l, y);
    sqrt(w.iter().map(|x| x * x).sum::<f64>() / N as f64)
}

/// Factor of a factored covariance plus a positive diagonal
///
/// Finds the lower-triangular L, with a positive diagonal, such that
/// L*L^T = F*F^T + diag(m)^2. A Householder QR of the 2N x N matrix
/// [F^T; diag(m)] gives an upper-triangular U with U^T*U = F*F^T + diag(m)^2,
/// and L is U^T with each column's sign set to make the diagonal positive.
/// With every m_i positive the stacked matrix has full column rank, so the
/// factor always exists. The matrix is divided by its largest entry first, so
/// no square in the QR can overflow.
///
/// # Arguments
/// * `f` - Any factor F \[N x N\]
/// * `m` - The square roots of the added diagonal, each positive \[N\]
///
/// # Returns
/// * `l` - The lower-triangular factor L
///
/// # Panics
/// * If N is above 6
fn stacked_factor<const N: usize>(f: &[[f64; N]; N], m: &[f64; N]) -> [[f64; N]; N] {
    assert!(N <= 6, "noise factors are at most 6 x 6");
    let scale = f
        .iter()
        .flatten()
        .chain(m)
        .fold(0.0, |acc: f64, x| acc.max(abs(*x)));

    // Rows 0..N hold F^T and rows N..2N hold diag(m), both divided by scale
    let mut a = [[0.0; N]; 12];
    for k in 0..N {
        for j in 0..N {
            a[k][j] = f[j][k] / scale;
        }
        a[N + k][k] = m[k] / scale;
    }
    let mut b = [0.0; 12];
    let u = householder_qr(&mut a[..2 * N], &mut b[..2 * N]).r;

    // L = U^T, scaled back, with row j of U negated where its diagonal is
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if j > i {
                0.0
            } else if u[j][j] < 0.0 {
                -u[j][i] * scale
            } else {
                u[j][i] * scale
            }
        })
    })
}

/// Convert a time to a split Julian date
///
/// Matches on the time scale exhaustively, so adding a scale forces a
/// conversion decision here.
///
/// # Arguments
/// * `datetime` - The time
///
/// # Returns
/// * `Some((jd, frac))` - The Julian day and day fraction from [`crate::utc2jday`]
/// * `None` - If the time is not a valid date
pub(crate) fn julian_date(datetime: &DateTime) -> Option<(f64, f64)> {
    match datetime.timezone {
        Timezone::UTC => utc2jday(datetime).ok(),
    }
}

/// Retrograde factor for an inclination
///
/// The equinoctial node vector is singular at one inclination: i = 180 degrees
/// for f_r = +1 and i = 0 for f_r = -1. Choosing f_r = +1 up to 90 degrees and
/// -1 above keeps every orbit at least 90 degrees from its singularity.
///
/// # Arguments
/// * `inclination` - The inclination \[rad\]
///
/// # Returns
/// * `f_r` - +1.0 for an inclination up to 90 degrees, -1.0 above
pub(crate) fn retrograde_factor(inclination: f64) -> f64 {
    if inclination <= PI / 2.0 { 1.0 } else { -1.0 }
}

/// Convert classical elements to equinoctial elements
///
/// Chooses the retrograde factor from the inclination with
/// [`retrograde_factor`].
///
/// # Arguments
/// * `classical` - The classical elements, in GP units
///
/// # Returns
/// * `(elements, f_r)` - The equinoctial elements and the retrograde factor
///   they use
pub(crate) fn classical_to_equinoctial(
    classical: &ClassicalElements,
) -> (EquinoctialElements, f64) {
    let deg = PI / 180.0;
    let i = classical.inclination * deg;
    let raan = classical.right_ascension_of_ascending_node * deg;
    let argp = classical.argument_of_perigee * deg;
    let ma = classical.mean_anomaly * deg;
    let f_r = retrograde_factor(i);

    // tan(i/2) = sin(i)/(1 + cos(i)), and its reciprocal sin(i)/(1 - cos(i)),
    // stay finite at both ends of the range
    let tan_half = if f_r > 0.0 {
        sin(i) / (1.0 + cos(i))
    } else {
        sin(i) / (1.0 - cos(i))
    };
    let lon_peri = argp + f_r * raan;
    let e = classical.eccentricity;

    let elements = EquinoctialElements {
        n: classical.mean_motion * 2.0 * PI / MINUTES_PER_DAY,
        a_f: e * cos(lon_peri),
        a_g: e * sin(lon_peri),
        chi: tan_half * sin(raan),
        psi: tan_half * cos(raan),
        lambda: rem_euclid(ma + lon_peri, 2.0 * PI),
    };
    (elements, f_r)
}

/// Convert equinoctial elements to classical elements
///
/// Implements equation 9 of the fit plan. The special cases are exact zeros
/// only: e = 0 sets the argument of perigee to 0, and chi = psi = 0 sets the
/// node to 0, which covers both degenerate inclinations (i = 0 for f_r = +1,
/// i = 180 degrees for f_r = -1). A small nonzero value needs no threshold,
/// because any argument of perigee gives the same orbit when the mean anomaly
/// is consistent with it.
///
/// # Arguments
/// * `elements` - The equinoctial elements
/// * `f_r` - The retrograde factor the elements use, +1.0 or -1.0
///
/// # Returns
/// * `classical` - The classical elements, angles in degrees in [0, 360) and
///   mean motion in revolutions per day
pub(crate) fn equinoctial_to_classical(
    elements: &EquinoctialElements,
    f_r: f64,
) -> ClassicalElements {
    let e = sqrt(elements.a_f * elements.a_f + elements.a_g * elements.a_g);

    // Node and inclination from the node vector; atan2 of two zeros is not
    // reliably 0 for negative zeros, so the exact-zero case is set directly
    let node_norm = sqrt(elements.chi * elements.chi + elements.psi * elements.psi);
    let raan = if node_norm == 0.0 {
        0.0
    } else {
        atan2(elements.chi, elements.psi)
    };
    let half_angle = 2.0 * atan2(node_norm, 1.0);
    let i = if f_r > 0.0 {
        half_angle
    } else {
        PI - half_angle
    };

    // Argument of perigee from the longitude of perigee w + f_r*O
    let argp = if e == 0.0 {
        0.0
    } else {
        atan2(elements.a_g, elements.a_f) - f_r * raan
    };
    let ma = elements.lambda - argp - f_r * raan;

    let rad = 180.0 / PI;
    ClassicalElements {
        mean_motion: elements.n * MINUTES_PER_DAY / (2.0 * PI),
        eccentricity: e,
        inclination: i * rad,
        right_ascension_of_ascending_node: wrap_degrees(raan * rad),
        argument_of_perigee: wrap_degrees(argp * rad),
        mean_anomaly: wrap_degrees(ma * rad),
    }
}

/// Retrograde factor for a state vector
///
/// Uses the sign of the angular momentum's z component, so an orbit counts as
/// prograde up to and including 90 degrees, as in [`retrograde_factor`].
///
/// # Arguments
/// * `r` - The position \[km\]
/// * `v` - The velocity \[km/s\]
///
/// # Returns
/// * `f_r` - +1.0 for a prograde or polar orbit, -1.0 for a retrograde one
pub(crate) fn state_retrograde_factor(r: &[f64; 3], v: &[f64; 3]) -> f64 {
    if cross(r, v)[2] >= 0.0 { 1.0 } else { -1.0 }
}

/// Convert a state vector to osculating equinoctial elements
///
/// Computes the two-body elements of the state directly in equinoctial form,
/// without passing through the classical angles, so circular and equatorial
/// states convert without loss. The mean longitude comes from the eccentric
/// longitude F by Kepler's equation in equinoctial form,
/// lambda = F + a_g*cos(F) - a_f*sin(F).
///
/// # Arguments
/// * `r` - The position \[km\]
/// * `v` - The velocity \[km/s\]
/// * `mu` - The gravitational parameter \[km^3/s^2\]
/// * `f_r` - The retrograde factor to use, +1.0 or -1.0
///
/// # Returns
/// * `Some(elements)` - The osculating elements, with `n` the two-body mean
///   motion and `lambda` in [0, 2*pi)
/// * `None` - If the state isn't an ellipse (zero angular momentum, or a
///   parabolic or hyperbolic orbit), or lies exactly on the singular
///   inclination of `f_r`
///
/// # References
/// - R. A. Broucke and P. J. Cefola, "On the Equinoctial Orbit Elements,"
///   *Celestial Mechanics* 5, 1972
/// - D. A. Danielson et al., *Semianalytic Satellite Theory*, NPS-MA-95-002, 1995, Section 2.1
pub(crate) fn state_to_equinoctial(
    r: &[f64; 3],
    v: &[f64; 3],
    mu: f64,
    f_r: f64,
) -> Option<EquinoctialElements> {
    let r_norm = norm(r);
    let h = cross(r, v);
    let h_norm = norm(&h);
    if r_norm == 0.0 || h_norm == 0.0 {
        return None;
    }

    // Semi-major axis from the vis-viva equation; only ellipses have 1/a > 0
    let inv_a = 2.0 / r_norm - dot(v, v) / mu;
    if inv_a.is_nan() || inv_a <= 0.0 {
        return None;
    }
    let a = 1.0 / inv_a;

    // Node vector from the orbit normal w: chi = w_x/(1 + f_r*w_z),
    // psi = -w_y/(1 + f_r*w_z)
    let w = [h[0] / h_norm, h[1] / h_norm, h[2] / h_norm];
    let denom = 1.0 + f_r * w[2];
    if denom == 0.0 {
        return None;
    }
    let chi = w[0] / denom;
    let psi = -w[1] / denom;

    // Equinoctial frame axes f and g in the orbit plane
    let k = 1.0 / (1.0 + chi * chi + psi * psi);
    let f_hat = [
        k * (1.0 - chi * chi + psi * psi),
        k * 2.0 * chi * psi,
        k * -2.0 * f_r * chi,
    ];
    let g_hat = [
        k * 2.0 * f_r * chi * psi,
        k * f_r * (1.0 + chi * chi - psi * psi),
        k * 2.0 * psi,
    ];

    // Eccentricity vector (v x h)/mu - r/|r|, resolved on f and g
    let vxh = cross(v, &h);
    let e_vec = [
        vxh[0] / mu - r[0] / r_norm,
        vxh[1] / mu - r[1] / r_norm,
        vxh[2] / mu - r[2] / r_norm,
    ];
    let a_f = dot(&e_vec, &f_hat);
    let a_g = dot(&e_vec, &g_hat);

    // Eccentric longitude F from the position on f and g
    let x1 = dot(r, &f_hat);
    let y1 = dot(r, &g_hat);
    // Roundoff can push a near-parabolic state's eccentricity to 1 or above
    let one_minus_e_sq = 1.0 - a_f * a_f - a_g * a_g;
    if one_minus_e_sq.is_nan() || one_minus_e_sq <= 0.0 {
        return None;
    }
    let b = sqrt(one_minus_e_sq);
    let beta = 1.0 / (1.0 + b);
    let cos_f = a_f + ((1.0 - a_f * a_f * beta) * x1 - a_f * a_g * beta * y1) / (a * b);
    let sin_f = a_g + ((1.0 - a_g * a_g * beta) * y1 - a_f * a_g * beta * x1) / (a * b);
    let ecc_lon = atan2(sin_f, cos_f);
    let lambda = ecc_lon + a_g * cos(ecc_lon) - a_f * sin(ecc_lon);

    Some(EquinoctialElements {
        n: sqrt(mu * inv_a * inv_a * inv_a) * 60.0,
        a_f,
        a_g,
        chi,
        psi,
        lambda: rem_euclid(lambda, 2.0 * PI),
    })
}

/// Minutes from the fit epoch to a time, from split Julian dates
///
/// Implements tau = ((jd - jd0) + (frac - frac0)) * 1440. Subtracting the
/// whole days and the fractions separately keeps microsecond precision; a
/// single summed Julian date resolves only about 40 microseconds.
///
/// # Arguments
/// * `epoch` - The fit epoch as (Julian day, day fraction) from [`crate::utc2jday`]
/// * `time` - The time as (Julian day, day fraction) from [`crate::utc2jday`]
///
/// # Returns
/// * `tau` - Minutes from the epoch to the time, negative before the epoch \[min\]
pub(crate) fn minutes_since_epoch(epoch: (f64, f64), time: (f64, f64)) -> f64 {
    ((time.0 - epoch.0) + (time.1 - epoch.1)) * MINUTES_PER_DAY
}

/// Index of the latest of a set of times
///
/// The default fit epoch is the latest input observation time. Ties go to the
/// first of the equal times, so the choice is reproducible.
///
/// # Arguments
/// * `times` - The times as (Julian day, day fraction) from [`crate::utc2jday`]
///
/// # Returns
/// * `Some(index)` - The index of the latest time
/// * `None` - If `times` is empty
pub(crate) fn latest_time_index(times: &[(f64, f64)]) -> Option<usize> {
    let mut latest: Option<usize> = None;
    for (index, time) in times.iter().enumerate() {
        match latest {
            Some(best) if minutes_since_epoch(times[best], *time) <= 0.0 => {}
            _ => latest = Some(index),
        }
    }
    latest
}

/// Fit epoch from the options and the observation times
///
/// The epoch is the user's epoch if given, otherwise the latest observation
/// time, even if that observation is later rejected, so the epoch can't move
/// during a fit. Either way it is then rounded to the nearest TLE epoch tick
/// (1e-8 day, at most 0.43 ms away), so writing the result as a TLE doesn't
/// shift the orbit.
///
/// # Arguments
/// * `observation_times` - The observation times, each passing
///   [`crate::validate_datetime`]
/// * `epoch` - The user's epoch, if set
///
/// # Returns
/// * `Some(epoch)` - The fit epoch, on the TLE grid
/// * `None` - If there are no observations and no user epoch, a time can't be
///   converted to a Julian date, or rounding overflows the year
pub(crate) fn fit_epoch(
    observation_times: &[DateTime],
    epoch: Option<&DateTime>,
) -> Option<DateTime> {
    let chosen = match epoch {
        Some(user) => *user,
        None => {
            let jds = observation_times
                .iter()
                .map(utc2jday)
                .collect::<Result<Vec<_>, _>>()
                .ok()?;
            observation_times[latest_time_index(&jds)?]
        }
    };
    round_to_tle_epoch(&chosen)
}

/// Wrap an angle into [0, 360)
///
/// A tiny negative angle wraps to 360 in floating point, so that case is
/// folded back to 0.
///
/// # Arguments
/// * `angle` - The angle \[degrees\]
///
/// # Returns
/// * `wrapped` - The angle in [0, 360) \[degrees\]
fn wrap_degrees(angle: f64) -> f64 {
    let wrapped = rem_euclid(angle, 360.0);
    if wrapped >= 360.0 { 0.0 } else { wrapped }
}

/// Whether a value is finite and positive
///
/// # Arguments
/// * `x` - The value
///
/// # Returns
/// * `positive` - True if x is finite and above 0
fn is_positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// Whether every component of a 3-vector is finite
///
/// # Arguments
/// * `a` - The vector
///
/// # Returns
/// * `finite` - True if no component is NaN or infinite
fn all_finite(a: &[f64; 3]) -> bool {
    a.iter().all(|x| x.is_finite())
}

/// Diagonal matrix from its diagonal
///
/// # Arguments
/// * `d` - The diagonal entries
///
/// # Returns
/// * `m` - The N x N matrix with `d` on the diagonal and zeros elsewhere
fn diagonal<const N: usize>(d: [f64; N]) -> [[f64; N]; N] {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}

/// Dot product of two 3-vectors
///
/// # Arguments
/// * `a` - The first vector
/// * `b` - The second vector
///
/// # Returns
/// * `d` - The dot product
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors
///
/// # Arguments
/// * `a` - The first vector
/// * `b` - The second vector
///
/// # Returns
/// * `c` - The cross product a x b
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector
///
/// # Arguments
/// * `a` - The vector
///
/// # Returns
/// * `n` - The norm
fn norm(a: &[f64; 3]) -> f64 {
    sqrt(dot(a, a))
}

// ----------
// Unit Tests
// ----------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::WGS72;
    use crate::math::{abs, powf};

    /// Classical elements in GP units
    ///
    /// # Arguments
    /// * `n` - Mean motion \[revs/day\]
    /// * `e` - Eccentricity
    /// * `i` - Inclination \[degrees\]
    /// * `raan` - Right ascension of the ascending node \[degrees\]
    /// * `argp` - Argument of perigee \[degrees\]
    /// * `ma` - Mean anomaly \[degrees\]
    ///
    /// # Returns
    /// * `classical` - The elements
    fn classical(n: f64, e: f64, i: f64, raan: f64, argp: f64, ma: f64) -> ClassicalElements {
        ClassicalElements {
            mean_motion: n,
            eccentricity: e,
            inclination: i,
            right_ascension_of_ascending_node: raan,
            argument_of_perigee: argp,
            mean_anomaly: ma,
        }
    }

    /// Smallest difference between two angles, allowing for the 360-degree wrap
    ///
    /// # Arguments
    /// * `a` - The first angle \[degrees\]
    /// * `b` - The second angle \[degrees\]
    ///
    /// # Returns
    /// * `d` - The absolute difference, at most 180 \[degrees\]
    fn angle_diff(a: f64, b: f64) -> f64 {
        let d = rem_euclid(a - b, 360.0);
        if d > 180.0 { 360.0 - d } else { d }
    }

    /// Two-body TEME state from classical elements, by the perifocal frame
    ///
    /// # Arguments
    /// * `ce` - The elements; the mean motion sets the semi-major axis
    /// * `mu` - The gravitational parameter \[km^3/s^2\]
    ///
    /// # Returns
    /// * `(r, v)` - The position \[km\] and velocity \[km/s\]
    fn two_body_state(ce: &ClassicalElements, mu: f64) -> ([f64; 3], [f64; 3]) {
        let deg = PI / 180.0;
        let n = ce.mean_motion * 2.0 * PI / 86400.0;
        let a = powf(mu / (n * n), 1.0 / 3.0);
        let e = ce.eccentricity;
        let (i, raan, argp) = (
            ce.inclination * deg,
            ce.right_ascension_of_ascending_node * deg,
            ce.argument_of_perigee * deg,
        );

        // Solve Kepler's equation for the eccentric anomaly
        let ma = ce.mean_anomaly * deg;
        let mut ea = ma;
        for _ in 0..50 {
            ea -= (ea - e * sin(ea) - ma) / (1.0 - e * cos(ea));
        }
        let nu = 2.0 * atan2(sqrt(1.0 + e) * sin(ea / 2.0), sqrt(1.0 - e) * cos(ea / 2.0));

        // Perifocal position and velocity, then rotate by O, i and w
        let p = a * (1.0 - e * e);
        let r_pf = [
            p * cos(nu) / (1.0 + e * cos(nu)),
            p * sin(nu) / (1.0 + e * cos(nu)),
            0.0,
        ];
        let v_pf = [-sqrt(mu / p) * sin(nu), sqrt(mu / p) * (e + cos(nu)), 0.0];
        let (co, so, ci, si, cw, sw) = (cos(raan), sin(raan), cos(i), sin(i), cos(argp), sin(argp));
        let rot = [
            [co * cw - so * sw * ci, -co * sw - so * cw * ci, so * si],
            [so * cw + co * sw * ci, -so * sw + co * cw * ci, -co * si],
            [sw * si, cw * si, ci],
        ];
        let apply = |x: &[f64; 3]| -> [f64; 3] {
            core::array::from_fn(|row| (0..3).map(|k| rot[row][k] * x[k]).sum())
        };
        (apply(&r_pf), apply(&v_pf))
    }

    /// Prograde and retrograde orbits covering near-circular, near-equatorial and polar cases
    const CASES: [(f64, f64, f64, f64, f64, f64); 8] = [
        (15.5, 0.0005, 51.6, 247.4, 130.5, 325.0),
        (15.1, 0.0021, 2.47, 183.9, 236.4, 60.4),
        (1.0027, 0.0005, 0.057, 356.8, 344.5, 218.8),
        (2.006, 0.74, 63.4, 10.0, 270.0, 5.0),
        (14.2, 0.01, 90.0, 120.0, 45.0, 300.0),
        (14.6, 0.002, 98.7, 30.0, 80.0, 200.0),
        (13.0, 0.05, 140.0, 300.0, 10.0, 100.0),
        (15.0, 0.001, 179.9, 75.0, 200.0, 20.0),
    ];

    /// The retrograde factor is +1 up to 90 degrees and -1 above, for elements and states
    ///
    /// # Panics
    /// * If a factor has the wrong sign
    #[test]
    fn test_retrograde_factor() {
        let deg = PI / 180.0;
        assert_eq!(retrograde_factor(0.0), 1.0);
        assert_eq!(retrograde_factor(90.0 * deg), 1.0);
        assert_eq!(retrograde_factor(90.0001 * deg), -1.0);
        assert_eq!(retrograde_factor(PI), -1.0);
        for case in CASES {
            let ce = classical(case.0, case.1, case.2, case.3, case.4, case.5);
            let (r, v) = two_body_state(&ce, WGS72.mu);
            let expected = if case.2 <= 90.0 { 1.0 } else { -1.0 };
            assert_eq!(state_retrograde_factor(&r, &v), expected, "i = {}", case.2);
        }
    }

    /// Classical -> equinoctial -> classical and equinoctial -> classical ->
    /// equinoctial both round-trip to roundoff, prograde and retrograde
    ///
    /// # Panics
    /// * If an element or the retrograde factor isn't recovered
    #[test]
    fn test_classical_equinoctial_round_trip() {
        for case in CASES {
            let ce = classical(case.0, case.1, case.2, case.3, case.4, case.5);
            let (eq, f_r) = classical_to_equinoctial(&ce);
            assert_eq!(f_r, if case.2 <= 90.0 { 1.0 } else { -1.0 });

            let back = equinoctial_to_classical(&eq, f_r);
            assert!(abs(back.mean_motion - ce.mean_motion) < 1e-13 * ce.mean_motion);
            assert!(abs(back.eccentricity - ce.eccentricity) < 1e-15);
            assert!(
                abs(back.inclination - ce.inclination) < 1e-11,
                "i = {}",
                case.2
            );
            assert!(
                angle_diff(
                    back.right_ascension_of_ascending_node,
                    ce.right_ascension_of_ascending_node
                ) < 1e-10
            );
            assert!(angle_diff(back.argument_of_perigee, ce.argument_of_perigee) < 1e-8);
            assert!(angle_diff(back.mean_anomaly, ce.mean_anomaly) < 1e-8);

            let again = classical_to_equinoctial(&back).0;
            for (x, y) in [
                (again.n, eq.n),
                (again.a_f, eq.a_f),
                (again.a_g, eq.a_g),
                (again.chi, eq.chi),
                (again.psi, eq.psi),
            ] {
                assert!(
                    abs(x - y) < 1e-14 * (1.0 + abs(y)),
                    "i = {}: {x} vs {y}",
                    case.2
                );
            }
            assert!(abs(rem_euclid(again.lambda - eq.lambda + PI, 2.0 * PI) - PI) < 1e-13);
        }
    }

    /// Exact zeros give w = 0 and O = 0, for i = 0 with f_r = +1 and i = 180 with f_r = -1
    ///
    /// # Panics
    /// * If an angle isn't exactly zero or the mean anomaly isn't the mean longitude
    #[test]
    fn test_exact_zero_special_cases() {
        for (f_r, inclination) in [(1.0, 0.0), (-1.0, 180.0)] {
            for (a_f, a_g, chi, psi) in [(0.0, 0.0, 0.0, 0.0), (-0.0, -0.0, -0.0, -0.0)] {
                let eq = EquinoctialElements {
                    n: 0.06,
                    a_f,
                    a_g,
                    chi,
                    psi,
                    lambda: 1.0,
                };
                let ce = equinoctial_to_classical(&eq, f_r);
                assert_eq!(ce.eccentricity, 0.0);
                assert_eq!(ce.inclination, inclination);
                assert_eq!(ce.argument_of_perigee, 0.0);
                assert_eq!(ce.right_ascension_of_ascending_node, 0.0);
                assert!(abs(ce.mean_anomaly - 180.0 / PI) < 1e-12);
            }
        }

        // e = 0 with a node: w = 0 and M takes up the rest of the longitude
        let eq = EquinoctialElements {
            n: 0.06,
            a_f: 0.0,
            a_g: 0.0,
            chi: 0.3,
            psi: 0.4,
            lambda: 2.0,
        };
        let ce = equinoctial_to_classical(&eq, 1.0);
        assert_eq!(ce.argument_of_perigee, 0.0);
        let raan = atan2(0.3, 0.4);
        assert!(abs(ce.mean_anomaly - wrap_degrees((2.0 - raan) * 180.0 / PI)) < 1e-12);
    }

    /// Exact-zero elements propagate as the conversion intends
    ///
    /// With chi = psi = 0 the node is set to 0 and the longitude of perigee
    /// w + f_r*O is kept, so any other node with the same longitude of perigee
    /// gives the same orbit. With e = 0, w = 0 is a convention: SGP4 floors the
    /// eccentricity at 1e-6 along w, so a different w moves the orbit by up to
    /// about 2*a*1e-6.
    ///
    /// # Panics
    /// * If a node split moves the orbit by 1 mm or more, or a perigee split at
    ///   e = 0 moves it by more than the eccentricity floor allows
    #[test]
    fn test_exact_zero_propagation() {
        use crate::gp::GenPerturbElementSet;
        use crate::sgp4::{init_sgp4, sgp4_prop_delta};
        use crate::time::{DateTime, Timezone};

        let epoch = DateTime {
            year: 2024,
            month: 3,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };
        let build = |ce: &ClassicalElements, raan: f64, argp: f64, ma: f64| {
            let gp = GenPerturbElementSet {
                epoch_datetime: epoch,
                inclination: ce.inclination,
                right_ascension_of_ascending_node: wrap_degrees(raan),
                eccentricity: ce.eccentricity,
                argument_of_perigee: wrap_degrees(argp),
                mean_anomaly: wrap_degrees(ma),
                mean_motion: ce.mean_motion,
                bstar: 1e-4,
                classification: 'U',
                ..Default::default()
            };
            init_sgp4(&gp, None).expect("valid elements")
        };
        let distance = |a: &crate::sgp4::Sgp4, b: &crate::sgp4::Sgp4, tsince: f64| {
            let s1 = sgp4_prop_delta(a, tsince).expect("propagates");
            let s2 = sgp4_prop_delta(b, tsince).expect("propagates");
            norm(&core::array::from_fn(|k| s1.position[k] - s2.position[k]))
        };
        let n = 15.2 * 2.0 * PI / MINUTES_PER_DAY;

        for (f_r, inclination) in [(1.0, 0.0), (-1.0, 180.0)] {
            // chi = psi = 0 with e = 0.001: O = 0, and moving O by 50 degrees
            // with w moved to keep w + f_r*O gives the same orbit
            let eq = EquinoctialElements {
                n,
                a_f: 0.0008,
                a_g: 0.0006,
                chi: 0.0,
                psi: 0.0,
                lambda: 1.0,
            };
            let ce = equinoctial_to_classical(&eq, f_r);
            assert_eq!(ce.inclination, inclination);
            assert_eq!(ce.right_ascension_of_ascending_node, 0.0);
            let zeroed = build(&ce, 0.0, ce.argument_of_perigee, ce.mean_anomaly);
            let node_split = build(
                &ce,
                50.0,
                ce.argument_of_perigee - f_r * 50.0,
                ce.mean_anomaly,
            );
            for tsince in [0.0, 45.0, 720.0, 1440.0] {
                let d = distance(&zeroed, &node_split, tsince);
                assert!(
                    d < 1e-6,
                    "i = {inclination}, t = {tsince}: node split {d} km"
                );
            }

            // e = 0: w = 0, and any other w differs only through the 1e-6 floor
            let eq = EquinoctialElements {
                a_f: 0.0,
                a_g: 0.0,
                ..eq
            };
            let ce = equinoctial_to_classical(&eq, f_r);
            assert_eq!(ce.argument_of_perigee, 0.0);
            let zeroed = build(&ce, 0.0, 0.0, ce.mean_anomaly);
            let perigee_split = build(&ce, 0.0, 37.0, ce.mean_anomaly - 37.0);
            let a_km = powf(WGS72.mu / (n / 60.0 * n / 60.0), 1.0 / 3.0);
            let d = distance(&zeroed, &perigee_split, 0.0);
            assert!(
                d < 2.2 * a_km * 1e-6,
                "i = {inclination}: perigee split {d} km"
            );
        }
    }

    /// The state conversion matches the classical conversion on two-body states
    ///
    /// # Panics
    /// * If an equinoctial element from the state differs from the one from the
    ///   classical elements
    #[test]
    fn test_state_to_equinoctial() {
        for case in CASES {
            let ce = classical(case.0, case.1, case.2, case.3, case.4, case.5);
            let (expected, f_r) = classical_to_equinoctial(&ce);
            let (r, v) = two_body_state(&ce, WGS72.mu);
            let got = state_to_equinoctial(&r, &v, WGS72.mu, f_r).expect("elliptical state");

            assert!(
                abs(got.n - expected.n) < 1e-12 * expected.n,
                "i = {}: n",
                case.2
            );
            for (name, x, y) in [
                ("a_f", got.a_f, expected.a_f),
                ("a_g", got.a_g, expected.a_g),
                ("chi", got.chi, expected.chi),
                ("psi", got.psi, expected.psi),
            ] {
                assert!(
                    abs(x - y) < 1e-11 * (1.0 + abs(y)),
                    "i = {}: {name} {x} vs {y}",
                    case.2
                );
            }
            let dl = abs(rem_euclid(got.lambda - expected.lambda + PI, 2.0 * PI) - PI);
            assert!(dl < 1e-10, "i = {}: lambda off by {dl}", case.2);
        }
    }

    /// Non-elliptical and degenerate states are rejected
    ///
    /// # Panics
    /// * If a hyperbolic, parabolic, radial, zero or singular state converts
    #[test]
    fn test_state_to_equinoctial_rejects() {
        let mu = WGS72.mu;
        let r = [7000.0, 0.0, 0.0];
        let v_esc = sqrt(2.0 * mu / 7000.0);
        assert!(state_to_equinoctial(&r, &[0.0, v_esc * 1.1, 0.0], mu, 1.0).is_none());
        assert!(state_to_equinoctial(&r, &[0.0, v_esc, 0.0], mu, 1.0).is_none());
        assert!(state_to_equinoctial(&r, &[3.0, 0.0, 0.0], mu, 1.0).is_none());
        // Nearly radial and bound: e^2 = 1 - 1e-20 rounds to 1
        assert!(state_to_equinoctial(&r, &[5.0, 1e-9, 0.0], mu, 1.0).is_none());
        assert!(state_to_equinoctial(&[0.0; 3], &[0.0, 7.5, 0.0], mu, 1.0).is_none());
        // Exactly retrograde equatorial with f_r = +1 is the singular inclination
        assert!(state_to_equinoctial(&r, &[0.0, -7.5, 0.0], mu, 1.0).is_none());
        assert!(state_to_equinoctial(&r, &[0.0, -7.5, 0.0], mu, -1.0).is_some());
    }

    /// Minutes from the epoch keep microsecond precision, and the latest time is found
    ///
    /// # Panics
    /// * If a 1-microsecond offset isn't resolved, or the wrong index is chosen
    #[test]
    fn test_minutes_since_epoch_and_latest() {
        use crate::time::{DateTime, Timezone, utc2jday};
        let t = |second: f64| {
            utc2jday(&DateTime {
                year: 2025,
                month: 6,
                day: 30,
                hour: 23,
                minute: 59,
                second,
                timezone: Timezone::UTC,
            })
            .expect("valid date")
        };
        let epoch = t(10.0);
        let later = t(10.000001);
        let tau = minutes_since_epoch(epoch, later);
        assert!(
            abs(tau * 60.0 - 1e-6) < 1e-10,
            "1 us came out as {} s",
            tau * 60.0
        );
        assert_eq!(minutes_since_epoch(later, epoch), -tau);

        assert_eq!(latest_time_index(&[]), None);
        assert_eq!(latest_time_index(&[t(5.0), later, t(1.0), epoch]), Some(1));
        // Equal times: the first one wins
        assert_eq!(latest_time_index(&[t(1.0), later, later]), Some(1));
    }

    /// The fit epoch is the user's epoch or the latest observation, rounded to the TLE grid
    ///
    /// # Panics
    /// * If the wrong time is chosen, the result isn't on the grid, or a
    ///   missing epoch isn't reported
    #[test]
    fn test_fit_epoch() {
        use crate::time::{DateTime, Timezone};
        let dt = |hour, second| DateTime {
            year: 2025,
            month: 6,
            day: 30,
            hour,
            minute: 0,
            second,
            timezone: Timezone::UTC,
        };
        let times = [dt(3, 1.0), dt(9, 0.0001), dt(6, 0.0)];

        // Default: the latest observation, 09:00:00.0001 rounds to 09:00:00
        let epoch = fit_epoch(&times, None).expect("an epoch");
        assert_eq!((epoch.hour, epoch.minute, epoch.second), (9, 0, 0.0));

        // A user epoch wins, even outside the observation span, and is rounded
        // to the nearest tick: 23:00:59.9997 is 95,902,777.43 ticks, so tick
        // 95,902,777 at 23:00:59.999328
        let user = dt(23, 59.9997);
        let epoch = fit_epoch(&times, Some(&user)).expect("an epoch");
        assert_eq!((epoch.day, epoch.hour, epoch.minute), (30, 23, 0));
        assert!(
            abs(epoch.second - 59.999328) < 1e-9,
            "second {}",
            epoch.second
        );

        // No observations and no user epoch
        assert!(fit_epoch(&[], None).is_none());
        assert!(fit_epoch(&[], Some(&user)).is_some());
    }

    /// A UTC time on 2024-03-01 during the 12:00 hour
    ///
    /// # Arguments
    /// * `minute` - The minute
    /// * `second` - The second
    ///
    /// # Returns
    /// * `datetime` - The time
    fn time_at(minute: i32, second: f64) -> DateTime {
        DateTime {
            year: 2024,
            month: 3,
            day: 1,
            hour: 12,
            minute,
            second,
            timezone: Timezone::UTC,
        }
    }

    /// A LEO-like TEME state observation with the given noise
    ///
    /// # Arguments
    /// * `minute` - The minute of the observation time
    /// * `noise` - The measurement noise
    ///
    /// # Returns
    /// * `obs` - The observation
    fn state_at(minute: i32, noise: Option<FitStateNoise>) -> FitObservation {
        new_fit_observation(
            time_at(minute, 0.0),
            FitMeasurement::State {
                frame: CoordinateFrame::TEME,
                position: [6778.0, 0.0, 0.0],
                velocity: [0.0, 4.8, 5.9],
                noise,
            },
        )
    }

    /// A LEO-like TEME position observation with the given noise
    ///
    /// # Arguments
    /// * `minute` - The minute of the observation time
    /// * `noise` - The measurement noise
    ///
    /// # Returns
    /// * `obs` - The observation
    fn position_at(minute: i32, noise: Option<FitPositionNoise>) -> FitObservation {
        new_fit_observation(
            time_at(minute, 0.0),
            FitMeasurement::Position {
                frame: CoordinateFrame::TEME,
                position: [6778.0, 0.0, 0.0],
                noise,
            },
        )
    }

    /// A rotation matrix from three angles
    ///
    /// # Arguments
    /// * `a` - Rotation about z \[rad\]
    /// * `b` - Rotation about x \[rad\]
    /// * `c` - Rotation about z \[rad\]
    ///
    /// # Returns
    /// * `rot` - Rz(a)*Rx(b)*Rz(c)
    fn rotation(a: f64, b: f64, c: f64) -> [[f64; 3]; 3] {
        let (ca, sa, cb, sb, cc, sc) = (cos(a), sin(a), cos(b), sin(b), cos(c), sin(c));
        [
            [ca * cc - sa * cb * sc, -ca * sc - sa * cb * cc, sa * sb],
            [sa * cc + ca * cb * sc, -sa * sc + ca * cb * cc, -ca * sb],
            [sb * sc, sb * cc, cb],
        ]
    }

    /// Rotate a diagonal covariance: rot * diag(variances) * rot^T
    ///
    /// # Arguments
    /// * `rot` - The rotation matrix
    /// * `variances` - The variances along the rotated axes
    ///
    /// # Returns
    /// * `cov` - The rotated covariance, symmetric only to roundoff
    fn rotate_covariance(rot: &[[f64; 3]; 3], variances: [f64; 3]) -> [[f64; 3]; 3] {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| (0..3).map(|k| rot[i][k] * variances[k] * rot[j][k]).sum())
        })
    }

    /// Check that a factor reproduces a covariance
    ///
    /// Entry (i, j) of F*F^T must match the covariance to 1e-13 times d_i*d_j,
    /// with d_i the root mean variance of the 3-vector holding component i,
    /// which covers roundoff in every block.
    ///
    /// # Arguments
    /// * `f` - The factor \[N x N\]
    /// * `cov` - The covariance it should reproduce \[N x N\]
    ///
    /// # Panics
    /// * If F*F^T differs from `cov`
    fn assert_factor<const N: usize>(f: &[[f64; N]; N], cov: &[[f64; N]; N]) {
        let scale: [f64; N] = core::array::from_fn(|i| {
            let b = 3 * (i / 3);
            sqrt((cov[b][b] + cov[b + 1][b + 1] + cov[b + 2][b + 2]) / 3.0)
        });
        for i in 0..N {
            for j in 0..N {
                let fft: f64 = (0..N).map(|k| f[i][k] * f[j][k]).sum();
                let tol = 1e-13 * scale[i] * scale[j];
                assert!(
                    abs(fft - cov[i][j]) <= tol,
                    "entry ({i}, {j}): {fft} against {}",
                    cov[i][j]
                );
            }
        }
    }

    /// The default options have the documented values and validate
    ///
    /// # Panics
    /// * If a default differs from the plan or fails validation
    #[test]
    fn test_default_options() {
        let opts = FitOptions::default();
        assert!(opts.epoch.is_none());
        assert!(opts.initial_guess.is_none());
        assert_eq!(opts.metadata, GpMetadata::default());
        assert_eq!(opts.metadata.name, None);
        assert_eq!(opts.bstar_prior_center, None);
        assert_eq!(opts.bstar_prior_width, 1e-2);
        assert_eq!(opts.outlier_cutoff, Some(10.0));
        assert_eq!(opts.max_rejected_fraction, 0.25);
        assert!(opts.max_iterations >= 1);
        assert_eq!(opts.tolerance, 1e-3);
        assert_eq!(opts.initial_damping, 1e-3);
        assert_eq!(opts.max_damping, 1e10);
        assert!(opts.wgs.is_none());
        assert!(validate_options(&opts).is_ok());
    }

    /// Each out-of-range option is named, and the first in order wins
    ///
    /// # Panics
    /// * If an invalid option passes or the wrong field is named
    #[test]
    fn test_invalid_options() {
        type Mutate = fn(&mut FitOptions);
        let cases: [(Mutate, &str); 25] = [
            (|o| o.bstar_prior_width = 0.0, "bstar_prior_width"),
            (|o| o.bstar_prior_width = -1e-3, "bstar_prior_width"),
            (|o| o.bstar_prior_width = f64::NAN, "bstar_prior_width"),
            (|o| o.bstar_prior_width = f64::INFINITY, "bstar_prior_width"),
            (
                |o| o.bstar_prior_center = Some(f64::NAN),
                "bstar_prior_center",
            ),
            (
                |o| o.bstar_prior_center = Some(f64::NEG_INFINITY),
                "bstar_prior_center",
            ),
            (|o| o.outlier_cutoff = Some(2.99), "outlier_cutoff"),
            (|o| o.outlier_cutoff = Some(f64::NAN), "outlier_cutoff"),
            (|o| o.outlier_cutoff = Some(f64::INFINITY), "outlier_cutoff"),
            (|o| o.max_rejected_fraction = 0.0, "max_rejected_fraction"),
            (|o| o.max_rejected_fraction = 0.5, "max_rejected_fraction"),
            (
                |o| o.max_rejected_fraction = f64::NAN,
                "max_rejected_fraction",
            ),
            (|o| o.max_iterations = 0, "max_iterations"),
            (|o| o.tolerance = 0.0, "tolerance"),
            (|o| o.tolerance = f64::INFINITY, "tolerance"),
            (|o| o.initial_damping = 0.0, "initial_damping"),
            (|o| o.initial_damping = f64::NAN, "initial_damping"),
            (|o| o.max_damping = -1.0, "max_damping"),
            (|o| o.max_damping = f64::INFINITY, "max_damping"),
            (|o| o.initial_damping = 1e10, "max_damping"),
            (|o| o.initial_damping = 1e11, "max_damping"),
            (|o| o.epoch = Some(time_at(60, 0.0)), "epoch"),
            (|o| o.epoch = Some(time_at(0, 61.0)), "epoch"),
            (
                |o| {
                    let mut early = time_at(0, 0.0);
                    early.year = 1500;
                    o.epoch = Some(early)
                },
                "epoch",
            ),
            (
                |o| {
                    o.bstar_prior_width = 0.0;
                    o.epoch = Some(time_at(60, 0.0))
                },
                "bstar_prior_width",
            ),
        ];
        for (mutate, expected) in cases {
            let mut opts = FitOptions::default();
            mutate(&mut opts);
            match validate_options(&opts) {
                Err(FitError::InvalidOptions { field }) => assert_eq!(field, expected),
                other => panic!("{expected}: got {other:?}"),
            }
        }

        // Legal edge values pass: negative B*, the minimum cutoff, rejection off
        let mut opts = FitOptions {
            bstar_prior_center: Some(-1e-3),
            outlier_cutoff: Some(3.0),
            max_rejected_fraction: 0.499,
            max_iterations: 1,
            epoch: Some(time_at(0, 60.5)),
            ..FitOptions::default()
        };
        assert!(validate_options(&opts).is_ok());
        opts.outlier_cutoff = None;
        assert!(validate_options(&opts).is_ok());
    }

    /// Invalid values, times and noise give `InvalidObservation` with the right index
    ///
    /// # Panics
    /// * If an invalid observation passes or the wrong index is reported
    #[test]
    fn test_invalid_observation() {
        let nan = f64::NAN;
        let mut bad: Vec<FitObservation> = Vec::new();

        // Non-finite position or velocity
        let mut obs = state_at(0, None);
        if let FitMeasurement::State { position, .. } = &mut obs.measurement {
            position[1] = nan;
        }
        bad.push(obs);
        let mut obs = state_at(0, None);
        if let FitMeasurement::State { velocity, .. } = &mut obs.measurement {
            velocity[2] = f64::INFINITY;
        }
        bad.push(obs);
        let mut obs = position_at(0, None);
        if let FitMeasurement::Position { position, .. } = &mut obs.measurement {
            position[0] = f64::NEG_INFINITY;
        }
        bad.push(obs);

        // Invalid times: month 13, second 61, NaN second, before the Gregorian calendar
        for change in [
            |t: &mut DateTime| t.month = 13,
            |t: &mut DateTime| t.second = 61.0,
            |t: &mut DateTime| t.second = f64::NAN,
            |t: &mut DateTime| t.year = 1500,
        ] {
            let mut obs = position_at(0, None);
            change(&mut obs.datetime);
            bad.push(obs);
        }

        // Invalid sigmas: negative, NaN, infinite, and one whose square overflows
        for sigma in [-1e-3, nan, f64::INFINITY, 1e200] {
            bad.push(position_at(0, Some(FitPositionNoise::Sigma(sigma))));
            bad.push(state_at(
                0,
                Some(FitStateNoise::Sigmas {
                    position: 0.01,
                    velocity: sigma,
                }),
            ));
            bad.push(state_at(
                0,
                Some(FitStateNoise::Sigmas {
                    position: sigma,
                    velocity: 1e-5,
                }),
            ));
        }

        // Invalid position covariances: non-finite, asymmetric beyond the
        // tolerance, a negative eigenvalue, and a zero block with a nonzero
        // entry
        let id = diagonal([1.0; 3]);
        let mut non_finite = id;
        non_finite[0][2] = nan;
        let mut asymmetric = id;
        asymmetric[0][1] = 1e-9;
        let indefinite = [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let cancelling = diagonal([1.0, -1.0, 0.0]);
        for cov in [
            non_finite,
            asymmetric,
            indefinite,
            cancelling,
            diagonal([-1.0; 3]),
        ] {
            bad.push(position_at(0, Some(FitPositionNoise::Covariance(cov))));
        }

        // A negative variance inside the tolerance of a huge block is still
        // rejected
        bad.push(position_at(
            0,
            Some(FitPositionNoise::Covariance(diagonal([
                -2.5, 4.5e12, 4.5e12,
            ]))),
        ));

        // Invalid state covariances. A negative eigenvalue confined to the
        // velocity block is caught although it is tiny beside the position
        // variances; a test against the whole trace would hide it
        let mut velocity_indefinite = diagonal([1.0, 1.0, 1.0, 1e-6, 1e-6, 1e-6]);
        velocity_indefinite[3][4] = 2e-6;
        velocity_indefinite[4][3] = 2e-6;
        let mut zero_velocity_correlated = diagonal([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);
        zero_velocity_correlated[0][3] = 1e-3;
        zero_velocity_correlated[3][0] = 1e-3;
        let mut rv_asymmetric = diagonal([1.0, 1.0, 1.0, 1e-6, 1e-6, 1e-6]);
        rv_asymmetric[1][4] = 1e-4;
        rv_asymmetric[4][1] = 1e-4 + 1e-14;
        for cov in [velocity_indefinite, zero_velocity_correlated, rv_asymmetric] {
            bad.push(state_at(0, Some(FitStateNoise::Covariance(cov))));
        }

        // Each bad observation, placed second, is reported at index 1
        for (case, obs) in bad.iter().enumerate() {
            let observations = [state_at(1, None), *obs, state_at(2, None)];
            match validate_input(&observations, &FitOptions::default()) {
                Err(FitError::InvalidObservation { index }) => assert_eq!(index, 1, "case {case}"),
                other => panic!("case {case}: got {other:?}"),
            }
        }
    }

    /// Valid noise forms are accepted, factored and equal where equivalent
    ///
    /// # Panics
    /// * If valid noise is rejected or equivalent forms differ
    #[test]
    fn test_noise_accepted() {
        let checked_noise = |obs: &FitObservation| -> CheckedMeasurement {
            validate_observation(0, obs).expect("valid").measurement
        };

        // A sigma of 0 and an all-zero covariance equal no noise
        let none = checked_noise(&position_at(0, None));
        assert_eq!(
            none,
            checked_noise(&position_at(0, Some(FitPositionNoise::Sigma(0.0))))
        );
        assert_eq!(
            none,
            checked_noise(&position_at(
                0,
                Some(FitPositionNoise::Covariance([[0.0; 3]; 3]))
            ))
        );
        let none = checked_noise(&state_at(0, None));
        let zero_sigmas = FitStateNoise::Sigmas {
            position: 0.0,
            velocity: -0.0,
        };
        assert_eq!(none, checked_noise(&state_at(0, Some(zero_sigmas))));

        // A diagonal covariance equals the sigma form
        assert_eq!(
            checked_noise(&position_at(0, Some(FitPositionNoise::Sigma(0.1)))),
            checked_noise(&position_at(
                0,
                Some(FitPositionNoise::Covariance(diagonal([0.1 * 0.1; 3])))
            ))
        );
        let sigmas = FitStateNoise::Sigmas {
            position: 0.01,
            velocity: 1e-5,
        };
        let diag = FitStateNoise::Covariance(diagonal(core::array::from_fn(|i| {
            if i < 3 { 0.01 * 0.01 } else { 1e-5 * 1e-5 }
        })));
        assert_eq!(
            checked_noise(&state_at(0, Some(sigmas))),
            checked_noise(&state_at(0, Some(diag)))
        );

        // Rotated from RIC sigmas into TEME: anisotropic (10 m against 10 km)
        // and singular (one zero sigma), asymmetric at roundoff. All pass, and
        // the factor reproduces the covariance up to the shift
        let rot = rotation(0.7, 1.1, -2.3);
        for variances in [[1e-4, 100.0, 100.0], [0.0, 1.0, 1.0], [0.0, 0.0, 4.0]] {
            let cov = rotate_covariance(&rot, variances);
            let CheckedMeasurement::Position { sensor_factor, .. } =
                checked_noise(&position_at(0, Some(FitPositionNoise::Covariance(cov))))
            else {
                panic!("expected a position")
            };
            assert_factor(&sensor_factor, &cov);
        }

        // A state covariance with position-velocity correlation, the velocity
        // block in its own units, and a zero velocity block
        let mut correlated = diagonal([1.0, 1.0, 1.0, 1e-6, 1e-6, 1e-6]);
        correlated[0][3] = 0.9e-3;
        correlated[3][0] = 0.9e-3;
        let zero_velocity = diagonal([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);
        for cov in [correlated, zero_velocity] {
            let CheckedMeasurement::State { sensor_factor, .. } =
                checked_noise(&state_at(0, Some(FitStateNoise::Covariance(cov))))
            else {
                panic!("expected a state")
            };
            assert_factor(&sensor_factor, &cov);
        }
    }

    /// The set needs a full state or 3 distinct position times
    ///
    /// # Panics
    /// * If a set with too few observations passes or a sufficient one fails
    #[test]
    fn test_too_few_observations() {
        let too_few: [&[FitObservation]; 4] = [
            &[],
            &[position_at(0, None), position_at(1, None)],
            &[
                position_at(0, None),
                position_at(0, None),
                position_at(0, None),
            ],
            &[
                position_at(0, None),
                position_at(1, None),
                position_at(0, None),
            ],
        ];
        for (case, observations) in too_few.iter().enumerate() {
            assert!(
                matches!(
                    validate_input(observations, &FitOptions::default()),
                    Err(FitError::TooFewObservations)
                ),
                "case {case}"
            );
        }

        let enough: [&[FitObservation]; 5] = [
            &[state_at(0, None)],
            &[
                position_at(0, None),
                position_at(1, None),
                position_at(2, None),
            ],
            &[
                position_at(0, None),
                position_at(0, None),
                position_at(1, None),
                position_at(2, None),
            ],
            &[state_at(0, None), position_at(1, None)],
            &[state_at(0, None), state_at(0, None)],
        ];
        for (case, observations) in enough.iter().enumerate() {
            let checked = validate_input(observations, &FitOptions::default()).expect("enough");
            assert_eq!(checked.len(), observations.len(), "case {case}");
        }

        // Positions 1 microsecond apart are distinct times
        let mut close = [
            position_at(0, None),
            position_at(0, None),
            position_at(0, None),
        ];
        close[1].datetime.second = 1e-6;
        close[2].datetime.second = 2e-6;
        assert!(validate_input(&close, &FitOptions::default()).is_ok());

        // Options are checked first, then observations, then the count
        let opts = FitOptions {
            max_iterations: 0,
            ..FitOptions::default()
        };
        assert!(matches!(
            validate_input(&[], &opts),
            Err(FitError::InvalidOptions {
                field: "max_iterations"
            })
        ));
        let mut invalid = position_at(0, None);
        invalid.datetime.month = 0;
        assert!(matches!(
            validate_input(&[invalid], &FitOptions::default()),
            Err(FitError::InvalidObservation { index: 0 })
        ));
    }

    /// A mixed list of measurement kinds validates in input order, values kept
    ///
    /// # Panics
    /// * If a kind, value, time or the order is not kept
    #[test]
    fn test_mixed_measurements() {
        let sigma = FitPositionNoise::Sigma(0.1);
        let observations = [
            position_at(3, Some(sigma)),
            state_at(1, None),
            position_at(2, None),
        ];
        let checked = validate_input(&observations, &FitOptions::default()).expect("valid");

        assert_eq!(checked.len(), 3);
        for (obs, check) in observations.iter().zip(&checked) {
            assert_eq!(check.time, utc2jday(&obs.datetime).unwrap());
        }
        assert_eq!(
            checked[0].measurement,
            CheckedMeasurement::Position {
                position: [6778.0, 0.0, 0.0],
                sensor_factor: checked_noise_factor(&diagonal([0.1 * 0.1; 3])).unwrap(),
            }
        );
        assert_eq!(
            checked[1].measurement,
            CheckedMeasurement::State {
                position: [6778.0, 0.0, 0.0],
                velocity: [0.0, 4.8, 5.9],
                sensor_factor: [[0.0; 6]; 6],
            }
        );
        assert!(matches!(
            checked[2].measurement,
            CheckedMeasurement::Position { .. }
        ));
    }

    /// Error messages, the hand-written Debug and the error source
    ///
    /// # Panics
    /// * If a message, Debug output or source is wrong
    #[test]
    fn test_fit_error_traits() {
        use alloc::format;
        use alloc::string::ToString;
        use core::error::Error;

        let err = FitError::InvalidObservation { index: 4 };
        assert!(err.to_string().starts_with("observation 4 "));
        assert_eq!(format!("{err:?}"), "InvalidObservation { index: 4 }");
        assert!(err.source().is_none());

        let err = FitError::InvalidOptions { field: "epoch" };
        assert!(err.to_string().contains("`epoch`"));

        let err = FitError::Sgp4(Sgp4Error::SatelliteDecayed);
        assert_eq!(format!("{err:?}"), "Sgp4(SatelliteDecayed)");
        assert!(err.source().is_some());

        for err in [
            FitError::TooFewObservations,
            FitError::StartingGuessFailed,
            FitError::Underdetermined,
        ] {
            let message = err.to_string();
            assert!(
                message.starts_with(|c: char| c.is_ascii_lowercase()),
                "{message}"
            );
            assert!(!message.ends_with('.'), "{message}");
        }
    }

    /// Mean motion in rad/min from rev/day
    ///
    /// # Arguments
    /// * `revs_per_day` - The mean motion \[revs/day\]
    ///
    /// # Returns
    /// * `n` - The mean motion \[rad/min\]
    fn rad_per_min(revs_per_day: f64) -> f64 {
        revs_per_day * 2.0 * PI / MINUTES_PER_DAY
    }

    /// The checked form of an observation that must validate
    ///
    /// # Arguments
    /// * `obs` - The observation
    ///
    /// # Returns
    /// * `checked` - The checked observation
    ///
    /// # Panics
    /// * If the observation fails validation
    fn checked(obs: &FitObservation) -> CheckedObservation {
        validate_observation(0, obs).expect("valid observation")
    }

    /// L_j reproduces the measurement noise plus the model-error term
    ///
    /// # Panics
    /// * If L_j*L_j^T differs from equation 3 or a velocity sigma is off
    #[test]
    fn test_noise_factor_model() {
        // No measurement noise: velocity sigma n0 * sigma_model, about 1.1 m/s
        // in LEO and 0.07 m/s at GEO
        for (revs_per_day, sigma_v) in [(15.5, 1.127e-3), (1.0027, 7.29e-5)] {
            let n0 = rad_per_min(revs_per_day);
            let NoiseFactor::State(l) = noise_factor(&checked(&state_at(0, None)).measurement, n0)
            else {
                panic!("expected a state")
            };
            let var_v = (n0 / 60.0) * (n0 / 60.0);
            assert_lower_positive(&l);
            assert_factor(&l, &diagonal([1.0, 1.0, 1.0, var_v, var_v, var_v]));
            assert!(abs(l[3][3] / sigma_v - 1.0) < 1e-3);
        }

        // A position's factor doesn't depend on n0
        let position = checked(&position_at(0, None)).measurement;
        let NoiseFactor::Position(l) = noise_factor(&position, rad_per_min(15.5)) else {
            panic!("expected a position")
        };
        assert_factor(&l, &diagonal([1.0; 3]));
        assert_eq!(
            noise_factor(&position, rad_per_min(1.0)),
            NoiseFactor::Position(l)
        );

        // Measurement noise is added to the model term, correlations included
        let n0 = rad_per_min(15.5);
        let var_v = (n0 / 60.0) * (n0 / 60.0);
        let mut cov = diagonal([0.04, 0.04, 0.04, 1e-8, 1e-8, 1e-8]);
        cov[1][4] = 1e-5;
        cov[4][1] = 1e-5;
        let obs = checked(&state_at(0, Some(FitStateNoise::Covariance(cov))));
        let NoiseFactor::State(l) = noise_factor(&obs.measurement, n0) else {
            panic!("expected a state")
        };
        let mut total = cov;
        for (i, row) in total.iter_mut().enumerate() {
            row[i] += if i < 3 { 1.0 } else { var_v };
        }
        assert_lower_positive(&l);
        assert_factor(&l, &total);
    }

    /// Equivalent noise forms give the same L_j, and L_j whitens
    ///
    /// # Panics
    /// * If equivalent forms differ or a whitened norm differs from the
    ///   rotated weights
    #[test]
    fn test_noise_forms_agree() {
        let n0 = rad_per_min(15.5);
        let model = |obs: &FitObservation| noise_factor(&checked(obs).measurement, n0);

        // Sigma, diagonal covariance, and zero noise against none
        let sigma = model(&position_at(0, Some(FitPositionNoise::Sigma(0.3))));
        let cov = FitPositionNoise::Covariance(diagonal([0.3 * 0.3; 3]));
        assert_eq!(sigma, model(&position_at(0, Some(cov))));
        assert_eq!(
            model(&position_at(0, None)),
            model(&position_at(0, Some(FitPositionNoise::Sigma(0.0))))
        );
        let sigmas = FitStateNoise::Sigmas {
            position: 0.3,
            velocity: 2e-4,
        };
        let cov = FitStateNoise::Covariance(diagonal(core::array::from_fn(|i| {
            if i < 3 { 0.3 * 0.3 } else { 2e-4 * 2e-4 }
        })));
        assert_eq!(
            model(&state_at(0, Some(sigmas))),
            model(&state_at(0, Some(cov)))
        );
        let zero = FitStateNoise::Sigmas {
            position: 0.0,
            velocity: 0.0,
        };
        assert_eq!(model(&state_at(0, None)), model(&state_at(0, Some(zero))));

        // A rotated, anisotropic covariance (10 m against 10 km) and a
        // singular one: |L^-1*y|^2 = sum_k (rot^T*y)_k^2 / (variance_k + 1),
        // up to the shift of about 1e-12 of the mean variance
        let rot = rotation(0.7, 1.1, -2.3);
        let y = [1.3, -0.4, 2.2];
        for variances in [[1e-4, 100.0, 100.0], [0.0, 0.0, 4.0]] {
            let noise = FitPositionNoise::Covariance(rotate_covariance(&rot, variances));
            let NoiseFactor::Position(l) = model(&position_at(0, Some(noise))) else {
                panic!("expected a position")
            };
            let expected: f64 = (0..3)
                .map(|k| {
                    let along: f64 = (0..3).map(|i| rot[i][k] * y[i]).sum();
                    along * along / (variances[k] + 1.0)
                })
                .sum();
            let rho = normalized_residual(&l, &y);
            assert!(abs(3.0 * rho * rho / expected - 1.0) < 1e-10);
        }
    }

    /// The same error in every component scores the same rho for a position
    /// and a state
    ///
    /// # Panics
    /// * If the two normalized residuals differ
    #[test]
    fn test_outlier_normalization() {
        let n0 = rad_per_min(15.5);
        let sigma_v = n0 / 60.0;
        let factors = noise_factors(
            &[checked(&position_at(0, None)), checked(&state_at(0, None))],
            n0,
        );
        let (NoiseFactor::Position(l_pos), NoiseFactor::State(l_state)) = (factors[0], factors[1])
        else {
            panic!("expected a position then a state")
        };

        // A 40 km glitch in every component, velocity in its own sigma units
        let glitch = 40.0;
        let y_pos = [glitch; 3];
        let y_state: [f64; 6] =
            core::array::from_fn(|i| if i < 3 { glitch } else { glitch * sigma_v });
        let rho_pos = normalized_residual(&l_pos, &y_pos);
        let rho_state = normalized_residual(&l_state, &y_state);
        assert!(abs(rho_pos / glitch - 1.0) < 1e-12);
        assert!(abs(rho_state / rho_pos - 1.0) < 1e-12);

        // Without the sqrt(m_j) normalization the state would score sqrt(2) higher
        let raw = |rho: f64, m: f64| rho * sqrt(m);
        assert!(abs(raw(rho_state, 6.0) / raw(rho_pos, 3.0) - sqrt(2.0)) < 1e-12);
    }

    /// Every observation that passes validation has a usable noise factor
    ///
    /// # Panics
    /// * If a factor is not finite or not lower triangular with a positive
    ///   diagonal, or an accurate axis loses its weight
    #[test]
    fn test_noise_factor_extremes() {
        let n0 = rad_per_min(1.0027);
        let position_factor = |cov: [[f64; 3]; 3]| {
            let obs = checked(&position_at(0, Some(FitPositionNoise::Covariance(cov))));
            let NoiseFactor::Position(l) = noise_factor(&obs.measurement, n0) else {
                panic!("expected a position")
            };
            assert!(l.iter().flatten().all(|x| x.is_finite()));
            assert_lower_positive(&l);
            l
        };

        // A 10 m range axis beside 1e6 km cross-range axes, rotated into
        // TEME: the range axis keeps its weight, 1/(1e-4 + 1) per km^2. A
        // diagonal shift of 1e-12 of the mean variance would have added
        // about 0.7 km^2 to it
        let rot = rotation(0.7, 1.1, -2.3);
        let l = position_factor(rotate_covariance(&rot, [1e-4, 1e12, 1e12]));
        let range_axis = [rot[0][0], rot[1][0], rot[2][0]];
        let rho = normalized_residual(&l, &range_axis);
        assert!(abs(3.0 * rho * rho * (1e-4 + 1.0) - 1.0) < 1e-3);

        // Variances near the top of the f64 range don't overflow, including
        // a singular, rotated one
        let l = position_factor(diagonal([1e300; 3]));
        assert!(abs(l[0][0] / 1e150 - 1.0) < 1e-12);
        position_factor(rotate_covariance(&rot, [0.0, 1e300, 1e300]));
    }

    /// Check that a factor is lower triangular with a positive diagonal
    ///
    /// # Arguments
    /// * `l` - The factor \[N x N\]
    ///
    /// # Panics
    /// * If an entry above the diagonal is nonzero or a diagonal entry is not
    ///   positive
    fn assert_lower_positive<const N: usize>(l: &[[f64; N]; N]) {
        for (i, row) in l.iter().enumerate() {
            assert!(row[i] > 0.0, "diagonal entry {i} is {}", row[i]);
            for (j, x) in row.iter().enumerate().skip(i + 1) {
                assert_eq!(*x, 0.0, "entry ({i}, {j}) above the diagonal");
            }
        }
    }
}
