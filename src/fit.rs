//! Module for fitting a general perturbations element set to timestamped
//! states and positions.
//!
//! The fit works in equinoctial elements, which stay well defined for
//! circular and equatorial orbits where the classical argument of perigee and
//! node are not. This module holds the conversions between equinoctial
//! elements, the classical elements a TLE or OMM stores, and TEME states, and
//! the time helpers that place observations relative to the fit epoch.

// ------------------
// External Libraries
// ------------------
use alloc::vec::Vec;
use core::f64::consts::PI;

// ------------------
// Internal Libraries
// ------------------
use crate::math::{atan2, cos, rem_euclid, sin, sqrt};
use crate::time::{DateTime, round_to_tle_epoch, utc2jday};

// -------
// Structs
// -------

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

// ---------
// Constants
// ---------

/// Minutes in a day
///
/// Converts day counts to the minutes SGP4 propagates in.
const MINUTES_PER_DAY: f64 = 1440.0;

// ---------
// Functions
// ---------

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
    let b = sqrt(1.0 - a_f * a_f - a_g * a_g);
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
}
