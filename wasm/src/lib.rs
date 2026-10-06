//! WebAssembly bindings for mako-sgp4
//!
//! Exposes a [`Satellite`] class to JavaScript that parses a GP element set and propagates it
//! with SGP4/SDP4.

// ------------------
// External Libraries
// ------------------
use wasm_bindgen::prelude::*;

// ------------------
// Internal Libraries
// ------------------
use mako_sgp4::{
    DateTime, Sgp4, Timezone, WGS84, calc_theta_g, from_omm_kvn_string, from_tle_string,
    sgp4_prop_delta, utc2jday,
};

// -------
// Structs
// -------

/// SGP4 propagator for a single satellite
///
/// Wraps an initialized [`Sgp4`] model for use from JavaScript. All times are minutes since
/// the element set epoch, and all positions are in kilometers.
///
/// # Examples
/// ```rust
/// use mako_sgp4_wasm::Satellite;
///
/// // Parse a TLE
/// let tle = "\
/// ISS (ZARYA)
/// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
/// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sat = Satellite::new(tle).ok().unwrap();
///
/// // Propagate one orbit at 1 minute steps
/// let track = sat.track_teme(0.0, 92.0, 1.0);
/// assert_eq!(track.len(), 93 * 3);
/// ```
#[wasm_bindgen]
pub struct Satellite {
    /// Initialized SGP4 model
    sgp4: Sgp4,

    /// Element set epoch as Unix time \[ms\]
    epoch_unix_ms: f64,
}

// ---------
// Constants
// ---------

/// Maximum number of points in a track
///
/// Caps the memory a single track call can allocate, so a tiny step or a huge time span
/// cannot exhaust the browser's memory. One million points is about 24 MB.
const MAX_TRACK_POINTS: usize = 1_000_000;

// ---------
// Functions
// ---------

#[wasm_bindgen]
impl Satellite {
    /// Parse the first element set in a string
    ///
    /// Accepts a TLE (2 or 3 lines) or an OMM in KVN format. Text starting with
    /// `CCSDS_OMM_VERS` is parsed as OMM KVN, otherwise as TLE.
    ///
    /// # Arguments
    /// * `text` - Element set text
    ///
    /// # Returns
    /// * `Ok(Satellite)` - Initialized propagator for the first element set
    /// * `Err(JsError)` - If no element set could be parsed
    ///
    /// # Errors
    /// * If `text` is empty, cannot be parsed, or contains no element sets
    ///
    /// # Examples
    /// ```rust
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE without a name line
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    /// assert_eq!(sat.norad_id(), 25544);
    /// ```
    #[wasm_bindgen(constructor)]
    pub fn new(text: &str) -> Result<Satellite, JsError> {
        // Reject empty input
        let trimmed = text.trim();
        if trimmed.is_empty() {
            return Err(JsError::new("No element set provided"));
        }

        // Pick a parser from the header
        let parsed = if trimmed.starts_with("CCSDS_OMM_VERS") {
            from_omm_kvn_string(trimmed)
        } else {
            from_tle_string(trimmed)
        };

        // Keep the first element set
        let sgp4 = parsed
            .map_err(|e| JsError::new(&format!("Could not parse element set: {e:?}")))?
            .into_iter()
            .next()
            .ok_or_else(|| JsError::new("No element set found"))?;

        // Julian date of the Unix epoch
        let unix_epoch = DateTime {
            year: 1970,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.,
            timezone: Timezone::UTC,
        };
        let (jd_unix, jdfrac_unix) = utc2jday(&unix_epoch)
            .map_err(|e| JsError::new(&format!("Could not convert Unix epoch: {e:?}")))?;

        // Subtract whole days and day-fractions separately to limit rounding error
        let epoch_unix_ms = ((sgp4.jd0 - jd_unix) + (sgp4.jdfrac0 - jdfrac_unix)) * 86_400_000.;

        Ok(Satellite {
            sgp4,
            epoch_unix_ms,
        })
    }

    /// Object name
    ///
    /// # Returns
    /// * `name` - Object name, or an empty string if none was given
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.sgp4.gp.common_name.clone()
    }

    /// NORAD satellite catalog number
    ///
    /// # Returns
    /// * `norad_id` - Catalog number
    #[wasm_bindgen(getter, js_name = noradId)]
    pub fn norad_id(&self) -> i32 {
        self.sgp4.gp.satellite_catalog_number
    }

    /// Element set epoch as Unix time
    ///
    /// Pass to `new Date()` in JavaScript.
    ///
    /// # Returns
    /// * `epoch_unix_ms` - Milliseconds since 1970-01-01 00:00:00 UTC \[ms\]
    #[wasm_bindgen(getter, js_name = epochUnixMs)]
    pub fn epoch_unix_ms(&self) -> f64 {
        self.epoch_unix_ms
    }

    /// Orbital inclination
    ///
    /// # Returns
    /// * `inclination` - Inclination \[deg\]
    #[wasm_bindgen(getter)]
    pub fn inclination(&self) -> f64 {
        self.sgp4.gp.inclination
    }

    /// Orbital eccentricity
    ///
    /// # Returns
    /// * `eccentricity` - Eccentricity \[\]
    #[wasm_bindgen(getter)]
    pub fn eccentricity(&self) -> f64 {
        self.sgp4.gp.eccentricity
    }

    /// Orbital period from the element set mean motion
    ///
    /// # Returns
    /// * `period` - Period \[min\]
    #[wasm_bindgen(getter, js_name = periodMinutes)]
    pub fn period_minutes(&self) -> f64 {
        1440. / self.sgp4.gp.mean_motion
    }

    /// Deep space flag
    ///
    /// # Returns
    /// * `deep_space` - True if propagated with SDP4 (period of 225 minutes or more)
    #[wasm_bindgen(getter, js_name = deepSpace)]
    pub fn deep_space(&self) -> bool {
        self.sgp4.deep_space
    }

    /// Convert Unix time to minutes since epoch
    ///
    /// # Arguments
    /// * `unix_ms` - Milliseconds since 1970-01-01 00:00:00 UTC, e.g. `Date.now()` \[ms\]
    ///
    /// # Returns
    /// * `t` - Minutes since epoch \[min\]
    ///
    /// # Examples
    /// ```rust
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    ///
    /// // One hour after epoch is 60 minutes
    /// let t = sat.minutes_since_epoch(sat.epoch_unix_ms() + 3_600_000.0);
    /// assert!((t - 60.0).abs() < 1e-6);
    /// ```
    #[wasm_bindgen(js_name = minutesSinceEpoch)]
    pub fn minutes_since_epoch(&self, unix_ms: f64) -> f64 {
        (unix_ms - self.epoch_unix_ms()) / 60_000.
    }

    /// Greenwich mean sidereal time (GMST)
    ///
    /// Rotate an Earth model by this angle about its polar axis to align it with the TEME frame.
    ///
    /// # Arguments
    /// * `t` - Minutes since epoch \[min\]
    ///
    /// # Returns
    /// * `theta_g` - GMST \[rad\], wrapped to \[0, 2 * pi)
    ///
    /// # Examples
    /// ```rust
    /// use std::f64::consts::PI;
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    ///
    /// // GMST is always within one revolution
    /// let theta_g = sat.gmst(0.0);
    /// assert!((0.0..2.0 * PI).contains(&theta_g));
    /// ```
    pub fn gmst(&self, t: f64) -> f64 {
        calc_theta_g(self.sgp4.jd0, self.sgp4.jdfrac0 + t / 1440.)
    }

    /// Propagate to a single time
    ///
    /// # Arguments
    /// * `t` - Minutes since epoch \[min\]
    ///
    /// # Returns
    /// * `Ok(state)` - `[x, y, z, vx, vy, vz]` in TEME \[km, km/s\]
    /// * `Err(JsError)` - If propagation fails
    ///
    /// # Errors
    /// * If intermediate orbital elements become non-physical (e.g. the satellite decays)
    ///
    /// # Examples
    /// ```rust
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    ///
    /// // Propagate to epoch
    /// let state = sat.propagate(0.0).ok().unwrap();
    /// assert_eq!(state.len(), 6);
    /// ```
    pub fn propagate(&self, t: f64) -> Result<Vec<f64>, JsError> {
        let sv = sgp4_prop_delta(&self.sgp4, t)
            .map_err(|e| JsError::new(&format!("Propagation failed: {e:?}")))?;
        Ok(vec![sv.r_x, sv.r_y, sv.r_z, sv.v_x, sv.v_y, sv.v_z])
    }

    /// Propagate a TEME position track
    ///
    /// Stops early if propagation fails (e.g. the satellite decays) or after 1,000,000 points.
    ///
    /// # Arguments
    /// * `start` - First time, minutes since epoch \[min\]
    /// * `end` - Last time, minutes since epoch \[min\]
    /// * `step` - Time step \[min\]
    ///
    /// # Returns
    /// * `track` - `[x0, y0, z0, x1, y1, z1, ...]` in TEME \[km\], empty if `step` is not positive,
    ///   `end` is before `start`, or any input is not finite
    ///
    /// # Examples
    /// ```rust
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    ///
    /// // Propagate 90 minutes at 1 minute steps (91 points)
    /// let track = sat.track_teme(0.0, 90.0, 1.0);
    /// assert_eq!(track.len(), 91 * 3);
    /// ```
    #[wasm_bindgen(js_name = trackTeme)]
    pub fn track_teme(&self, start: f64, end: f64, step: f64) -> Vec<f64> {
        let mut track = Vec::new();
        for t in sample_times(start, end, step) {
            match sgp4_prop_delta(&self.sgp4, t) {
                Ok(sv) => track.extend_from_slice(&[sv.r_x, sv.r_y, sv.r_z]),
                Err(_) => break,
            }
        }
        track
    }

    /// Propagate a geodetic ground track
    ///
    /// Positions are rotated into an Earth-fixed frame by GMST (neglecting polar motion) and
    /// converted to geodetic coordinates on the WGS-84 ellipsoid. Stops early if propagation
    /// fails (e.g. the satellite decays) or after 1,000,000 points.
    ///
    /// # Arguments
    /// * `start` - First time, minutes since epoch \[min\]
    /// * `end` - Last time, minutes since epoch \[min\]
    /// * `step` - Time step \[min\]
    ///
    /// # Returns
    /// * `track` - `[lat0, lon0, alt0, ...]` \[deg, deg, km\], empty if `step` is not positive,
    ///   `end` is before `start`, or any input is not finite
    ///
    /// # Examples
    /// ```rust
    /// use mako_sgp4_wasm::Satellite;
    ///
    /// // Parse a TLE
    /// let tle = "\
    /// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
    /// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    /// let sat = Satellite::new(tle).ok().unwrap();
    ///
    /// // Altitude at epoch is a few hundred kilometers
    /// let track = sat.track_geodetic(0.0, 0.0, 1.0);
    /// assert!((300.0..450.0).contains(&track[2]));
    /// ```
    #[wasm_bindgen(js_name = trackGeodetic)]
    pub fn track_geodetic(&self, start: f64, end: f64, step: f64) -> Vec<f64> {
        let mut track = Vec::new();
        for t in sample_times(start, end, step) {
            match sgp4_prop_delta(&self.sgp4, t) {
                Ok(sv) => {
                    let (lat, lon, alt) = teme2geodetic(sv.r_x, sv.r_y, sv.r_z, self.gmst(t));
                    track.extend_from_slice(&[lat, lon, alt]);
                }
                Err(_) => break,
            }
        }
        track
    }
}

/// Generate evenly spaced sample times
///
/// The number of samples is capped at [`MAX_TRACK_POINTS`].
///
/// # Arguments
/// * `start` - First time \[min\]
/// * `end` - Last time, included if it falls on a step \[min\]
/// * `step` - Time step \[min\]
///
/// # Returns
/// * `times` - Iterator of times, empty if `step` is not positive, `end` is before `start`,
///   or any input is not finite
fn sample_times(start: f64, end: f64, step: f64) -> impl Iterator<Item = f64> {
    // Reject non-finite inputs, non-positive steps, and reversed spans (NaN fails every check)
    let span = end - start;
    let valid = step.is_finite() && step > 0. && span.is_finite() && span >= 0.;

    // Count samples from the span to avoid accumulating floating-point error
    let n = if valid {
        ((span / step).floor() as usize).min(MAX_TRACK_POINTS - 1) + 1
    } else {
        0
    };
    (0..n).map(move |k| start + k as f64 * step)
}

/// Convert a TEME position to geodetic coordinates
///
/// Rotates about the polar axis by GMST into an Earth-fixed frame (neglecting polar motion),
/// then solves for geodetic latitude on the WGS-84 ellipsoid by fixed-point iteration.
///
/// # Arguments
/// * `x` - TEME position X component \[km\]
/// * `y` - TEME position Y component \[km\]
/// * `z` - TEME position Z component \[km\]
/// * `theta_g` - Greenwich mean sidereal time \[rad\]
///
/// # Returns
/// * `(lat, lon, alt)` - Geodetic latitude \[deg\], longitude in \[-180, 180\] \[deg\], and altitude \[km\]
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
fn teme2geodetic(x: f64, y: f64, z: f64, theta_g: f64) -> (f64, f64, f64) {
    // Rotate about Z by -GMST into the Earth-fixed frame
    let (s, c) = theta_g.sin_cos();
    let x_ef = c * x + s * y;
    let y_ef = -s * x + c * y;

    // Longitude and distance from the polar axis
    let lon = y_ef.atan2(x_ef);
    let p = x_ef.hypot(y_ef);

    // Iterate on geodetic latitude, starting from the spherical estimate
    let e2 = WGS84.flattening * (2. - WGS84.flattening);
    let mut lat = z.atan2(p * (1. - e2));
    let mut n = WGS84.r_earth_eq;
    for _ in 0..6 {
        n = WGS84.r_earth_eq / (1. - e2 * lat.sin().powi(2)).sqrt();
        lat = (z + n * e2 * lat.sin()).atan2(p);
    }

    // Altitude above the ellipsoid, using the polar form near the poles
    let alt = if lat.cos().abs() > 1e-10 {
        p / lat.cos() - n
    } else {
        z.abs() - n * (1. - e2)
    };

    (lat.to_degrees(), lon.to_degrees(), alt)
}

// ----------
// Unit Tests
// ----------

#[cfg(test)]
mod tests {
    use super::*;

    const ISS: &str = "\
ISS (ZARYA)
1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
";

    #[test]
    fn parses_tle() {
        let sat = Satellite::new(ISS).ok().unwrap();
        assert_eq!(sat.norad_id(), 25544);
        assert_eq!(sat.name(), "ISS (ZARYA)");
        assert!(!sat.deep_space());
        assert!((sat.period_minutes() - 91.6).abs() < 0.1);
    }

    #[test]
    fn parses_omm_kvn() {
        let omm = "\
CCSDS_OMM_VERS = 2.0
OBJECT_NAME    = 2026-106A
OBJECT_ID      = 2026-106A
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
        let sat = Satellite::new(omm).ok().unwrap();
        assert_eq!(sat.norad_id(), 69097);
        assert_eq!(sat.name(), "2026-106A");
    }

    #[test]
    fn epoch_matches_unix_time() {
        // 2008 day 264.51782528 = 2008-09-20 12:25:40.104 UTC
        let sat = Satellite::new(ISS).ok().unwrap();
        assert!((sat.epoch_unix_ms() - 1_221_913_540_104.0).abs() < 1.0);
    }

    #[test]
    fn empty_track_for_bad_step() {
        let sat = Satellite::new(ISS).ok().unwrap();
        assert!(sat.track_teme(0.0, 90.0, 0.0).is_empty());
        assert!(sat.track_teme(90.0, 0.0, 1.0).is_empty());
    }

    #[test]
    fn track_is_capped() {
        let sat = Satellite::new(ISS).ok().unwrap();
        assert_eq!(
            sat.track_teme(0.0, 1440.0, 1e-9).len(),
            MAX_TRACK_POINTS * 3
        );
    }

    #[test]
    fn empty_track_for_non_finite_inputs() {
        assert_eq!(sample_times(0.0, f64::INFINITY, 1.0).count(), 0);
        assert_eq!(sample_times(f64::NAN, 90.0, 1.0).count(), 0);
        assert_eq!(sample_times(0.0, 90.0, f64::NAN).count(), 0);
        assert_eq!(sample_times(0.0, 90.0, f64::INFINITY).count(), 0);
    }

    #[test]
    fn geodetic_track_is_plausible() {
        let sat = Satellite::new(ISS).ok().unwrap();
        let track = sat.track_geodetic(0.0, 1440.0, 5.0);
        assert_eq!(track.len(), 289 * 3);
        for p in track.chunks(3) {
            // Geodetic latitude can exceed the inclination by up to about 0.2 deg
            assert!(p[0].abs() <= 51.9, "latitude {}", p[0]);
            assert!((-180.0..=180.0).contains(&p[1]), "longitude {}", p[1]);
            assert!((300.0..450.0).contains(&p[2]), "altitude {}", p[2]);
        }
    }

    #[test]
    fn geodetic_on_equator() {
        let (lat, lon, alt) = teme2geodetic(WGS84.r_earth_eq + 400.0, 0.0, 0.0, 0.0);
        assert!(lat.abs() < 1e-9 && lon.abs() < 1e-9);
        assert!((alt - 400.0).abs() < 1e-9);
    }

    #[test]
    fn geodetic_on_pole() {
        let b = WGS84.r_earth_eq * (1.0 - WGS84.flattening);
        let (lat, _, alt) = teme2geodetic(0.0, 0.0, b + 400.0, 0.0);
        assert!((lat - 90.0).abs() < 1e-9);
        assert!((alt - 400.0).abs() < 1e-6);
    }
}
