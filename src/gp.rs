//! Module to handle the input and processing of GP (General Perturbation)
//! elements. Element sets can be parsed from TLE (including Alpha-5 catalog
//! numbers) and OMM formats, and written back to TLE and OMM KVN, XML, JSON,
//! and CSV.

// ------------------
// External Libraries
// ------------------
use std::fs;

#[cfg(any(feature = "xml", feature = "json", feature = "csv"))]
use std::collections::HashMap;

#[cfg(feature = "xml")]
use roxmltree::{Document, Node};

#[cfg(feature = "json")]
use serde_json::Value;

#[cfg(feature = "csv")]
use csv::{ReaderBuilder, WriterBuilder};

// ------------------
// Internal Libraries
// ------------------
use crate::sgp4::{Sgp4, Sgp4Error, init_sgp4};
use crate::time::{DateTime, Timezone, dayofyr2utc};

// -------
// Structs
// -------

/// General Perturbation Element Set for an Earth-orbiting satellite.
///
/// This struct represents the parsed contents of a standard General Perturbation Element Set.
/// The General Perturbation Element Set is the standard set of orbital elements used with the
/// SGP4 propagator. Elements are commonly distributed as Two-Line Elements (TLEs) or Orbit
/// Mean-Elements Messages (OMMs).
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Parse a TLE; the GP element set is stored on the SGP4 struct
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// assert_eq!(sgp4.gp.common_name, "ISS (ZARYA)");
/// assert!((sgp4.gp.inclination - 51.6416).abs() < 1e-4);
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
#[derive(Default, Clone)]
pub struct GenPerturbElementSet {
    /// Common name of the satellite (e.g., "ISS (ZARYA)")
    pub common_name: String,

    /// NORAD satellite catalog number
    pub satellite_catalog_number: i32,

    /// Classification (`U` = Unclassified, `C` = Classified, `S` = Secret)
    pub classification: char,

    /// International designator in `Y-NP` form (launch year (4+ digits), launch number (3+ digits), piece (1+ characters))
    pub international_designator: String,

    /// Epoch UTC datetime
    pub epoch_datetime: DateTime,

    /// First time derivative of mean motion \[revs/day^2\]
    pub first_derivative_of_mean_motion: f64,

    /// Second time derivative of mean motion \[revs/day^3\]
    pub second_derivative_of_mean_motion: f64,

    /// B* drag term \[1/Earth radii\]
    pub bstar: f64,

    /// Ephemeris type (always zero)
    pub ephemeris_type: i32,

    /// Element set number
    pub element_set_number: i32,

    /// Orbital inclination \[degrees\]
    pub inclination: f64,

    /// Right ascension of the ascending node (RAAN) \[degrees\]
    pub right_ascension_of_ascending_node: f64,

    /// Orbital eccentricity \[\]
    pub eccentricity: f64,

    /// Argument of perigee \[degrees\]
    pub argument_of_perigee: f64,

    /// Mean anomaly \[degrees\]
    pub mean_anomaly: f64,

    /// Mean motion \[revs/day\]
    pub mean_motion: f64,

    /// Revolution number at epoch \[revs\]
    pub revolution_number_at_epoch: i64,
}

// -----
// Enums
// -----

/// GP errors
#[derive(Debug, Clone, PartialEq)]
pub enum GpError {
    /// TLE line 0 is invalid
    InvalidTleLine0,

    /// TLE line 1 is invalid
    InvalidTleLine1,

    /// TLE line 2 is invalid
    InvalidTleLine2,

    /// TLE epoch is invalid
    InvalidTleEpoch,

    /// TLE line 1 and line 2 have different NORAD catalog numbers
    MismatchedTleCatalog,

    /// A TLE file could not be read
    Io(String),

    /// TLE catalog number is negative or outside the Alpha-5 range
    InvalidTLECatalogNumber,

    /// TLE epoch datetime is before 1957 or after 2056
    InvalidTLEDateTime,

    /// A TLE data line is not 68 characters before the checksum
    InvalidTLELine,

    /// TLE international designator launch number or piece exceeds 3 characters
    InvalidTLEInternationalDesignator,

    /// An OMM field cannot be parsed as the requested type
    InvalidOmmField,

    /// OMM EPOCH cannot be parsed as a UTC datetime
    InvalidOmmEpoch,

    /// OMM XML document cannot be parsed
    InvalidOmmXml,

    /// OMM JSON cannot be parsed or is not an object or array of objects
    InvalidOmmJson,

    /// OMM CSV cannot be parsed
    InvalidOmmCsv,

    /// SGP4 initialization failed
    Sgp4(Sgp4Error),
}

// ------
// Traits
// ------

/// Convert an OMM field value into a typed GP field
///
/// Used by KVN and XML OMM parsing to turn a present field string into String,
/// char, integer, floating-point, or DateTime values, and to supply defaults
/// when a field is missing or empty.
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_string;
///
/// // OMM numbers may use Fortran D exponents; missing classification is U
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// INCLINATION    = 97.5103
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734D-3
/// ";
/// let sgp4 = &from_omm_kvn_string(omm).unwrap()[0];
///
/// // Assert typed field conversion and the unclassified default
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4.gp.classification, 'U');
/// assert!((sgp4.gp.bstar - 0.00039221734).abs() < 1e-12);
/// ```
trait FromOmm: Sized {
    /// Default used when the OMM field is missing or empty
    fn omm_default() -> Self;

    /// Parse a present, non-empty OMM field value
    fn from_omm(field: &str, value: &str) -> Result<Self, GpError>;
}

impl FromOmm for String {
    fn omm_default() -> Self {
        String::new()
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        Ok(value.to_string())
    }
}

impl FromOmm for char {
    fn omm_default() -> Self {
        'U'
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        Ok(value.chars().next().unwrap_or('U'))
    }
}

impl FromOmm for f64 {
    fn omm_default() -> Self {
        0.0
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        let normalized = value.replace(['D', 'd'], "E");
        normalized
            .parse::<f64>()
            .map_err(|_| GpError::InvalidOmmField)
    }
}

impl FromOmm for i32 {
    fn omm_default() -> Self {
        0
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        value.parse::<i32>().map_err(|_| GpError::InvalidOmmField)
    }
}

impl FromOmm for i64 {
    fn omm_default() -> Self {
        0
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        value.parse::<i64>().map_err(|_| GpError::InvalidOmmField)
    }
}

impl FromOmm for DateTime {
    fn omm_default() -> Self {
        DateTime::default()
    }

    fn from_omm(_field: &str, value: &str) -> Result<Self, GpError> {
        parse_omm_epoch(value)
    }
}

// ---------
// Constants
// ---------

/// Celestrak GP CSV column order
///
/// Header names written by [`to_omm_csv_string`] and [`to_omm_csv_file`],
/// matching Celestrak GP CSV.
///
/// # References
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "csv")]
const OMM_CSV_HEADERS: [&str; 17] = [
    "OBJECT_NAME",
    "OBJECT_ID",
    "EPOCH",
    "MEAN_MOTION",
    "ECCENTRICITY",
    "INCLINATION",
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY",
    "EPHEMERIS_TYPE",
    "CLASSIFICATION_TYPE",
    "NORAD_CAT_ID",
    "ELEMENT_SET_NO",
    "REV_AT_EPOCH",
    "BSTAR",
    "MEAN_MOTION_DOT",
    "MEAN_MOTION_DDOT",
];

// ---------
// Functions
// ---------

/// Builds a [`Sgp4`] struct from the lines of a Two-Line Element (TLE) set.
///
/// Given the two required TLE lines (line 1 and line 2), and an optional
/// name line (line 0), this function parses the input into a [`Sgp4`] struct.
///
/// # Arguments
/// * `line1` - The first TLE data line (TLE line 1)
/// * `line2` - The second TLE data line (TLE line 2)
/// * `line0` - Optional name line (TLE line 0)
///
/// # Returns
/// * `Ok(Sgp4)` - Struct containing the parsed SGP4 parameters
/// * `Err(GpError)` - If a TLE line or field cannot be parsed, or SGP4 initialization fails
///
/// # Errors
/// * [`GpError::InvalidTleLine0`] - If the optional name line is empty or longer than 24 characters
/// * [`GpError::InvalidTleLine1`] - If line 1 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleLine2`] - If line 2 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleEpoch`] - If the epoch cannot be converted to a UTC datetime
/// * [`GpError::MismatchedTleCatalog`] - If line 1 and line 2 have different NORAD catalog numbers
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Define a 3-line TLE (optional name, then line 1 and line 2)
/// let line0 = "ISS (ZARYA)";
/// let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
///
/// // Parse the TLE lines into an SGP4 propagator
/// let sgp4 = from_tle_lines(line1, line2, Some(line0)).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4.gp.satellite_catalog_number, 25544);
/// assert_eq!(sgp4.gp.common_name, "ISS (ZARYA)");
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_lines(line1: &str, line2: &str, line0: Option<&str>) -> Result<Sgp4, GpError> {
    // Create mutable General Perturbation Element Set struct
    let mut gp = GenPerturbElementSet::default();

    // Extract the common name of the satellite from line 0
    if let Some(name_line) = line0 {
        if name_line.is_empty() || name_line.len() > 24 {
            return Err(GpError::InvalidTleLine0);
        }
        gp.common_name = name_line.to_string();
    }

    // Line 1 must be 69 characters before field slices or checksum
    if line1.len() != 69 {
        return Err(GpError::InvalidTleLine1);
    }

    // Line 2 must be 69 characters before field slices or checksum
    if line2.len() != 69 {
        return Err(GpError::InvalidTleLine2);
    }

    // Validate the TLE checksum. A mismatch is a warning, not an error.
    if !tle_checksum(line1) || !tle_checksum(line2) {
        eprintln!("warning: TLE checksum failed; continuing with parse");
    }

    // Satellite catalog number from line 1
    gp.satellite_catalog_number =
        parse_tle_catalog(&line1[2..7]).ok_or(GpError::InvalidTleLine1)?;

    // Line 2 must carry the same NORAD catalog number
    let catalog2 = parse_tle_catalog(&line2[2..7]).ok_or(GpError::InvalidTleLine2)?;
    if catalog2 != gp.satellite_catalog_number {
        return Err(GpError::MismatchedTleCatalog);
    }

    // Classification
    gp.classification = line1[7..8]
        .trim()
        .parse::<char>()
        .map_err(|_| GpError::InvalidTleLine1)?;

    // International designator (COSPAR YYYY-NNNP)
    let intl_des = line1[9..17].trim();
    if intl_des.is_empty() {
        // Handle case where international designator is not present
        gp.international_designator = String::new();
    } else if intl_des.len() < 2 {
        return Err(GpError::InvalidTleLine1);
    } else {
        // Expand two-digit launch year and insert the COSPAR dash
        let launch_year = tle_full_year(
            intl_des[..2]
                .parse::<i32>()
                .map_err(|_| GpError::InvalidTleLine1)?,
        );
        gp.international_designator = format!("{}-{}", launch_year, &intl_des[2..]);
    }

    // Epoch year (last two numbers)
    let yr_two_digit: i32 = line1[18..20]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine1)?;
    let epoch_year = tle_full_year(yr_two_digit);

    // Epoch day of year
    let epoch_day: f64 = line1[20..32]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine1)?;

    // Epoch UTC datetime
    gp.epoch_datetime = dayofyr2utc(epoch_year, epoch_day).map_err(|_| GpError::InvalidTleEpoch)?;

    // 1st derivative of mean motion [revs/day^2]
    gp.first_derivative_of_mean_motion = line1[33..43]
        .trim()
        .parse::<f64>()
        .map_err(|_| GpError::InvalidTleLine1)?
        * 2.0;

    // 2nd derivative of mean motion [revs/days^3]
    // Sign is a single column (space, +, or -) and must not be trimmed
    let nddot_mantissa: f64 = if line1.as_bytes()[44] == b'-' {
        format!("-0.{}", line1[45..50].trim())
            .parse()
            .map_err(|_| GpError::InvalidTleLine1)?
    } else {
        format!("0.{}", line1[45..50].trim())
            .parse()
            .map_err(|_| GpError::InvalidTleLine1)?
    };
    let nddot_exp: i32 = line1[50..52]
        .parse()
        .map_err(|_| GpError::InvalidTleLine1)?;
    gp.second_derivative_of_mean_motion = nddot_mantissa * 10.0_f64.powi(nddot_exp) * 6.0;

    // B* [1/Earth Radii]
    // Sign is a single column (space, +, or -) and must not be trimmed
    let bstar_mantissa: f64 = if line1.as_bytes()[53] == b'-' {
        format!("-0.{}", line1[54..59].trim())
            .parse()
            .map_err(|_| GpError::InvalidTleLine1)?
    } else {
        format!("0.{}", line1[54..59].trim())
            .parse()
            .map_err(|_| GpError::InvalidTleLine1)?
    };
    let bstar_exp: i32 = line1[59..61]
        .parse()
        .map_err(|_| GpError::InvalidTleLine1)?;
    gp.bstar = bstar_mantissa * 10.0_f64.powi(bstar_exp);

    // Ephemeris type
    if line1[62..63].trim().is_empty() {
        // Handle case where ephemeris type is not present
        gp.ephemeris_type = 0;
    } else {
        gp.ephemeris_type = line1[62..63]
            .parse()
            .map_err(|_| GpError::InvalidTleLine1)?;
    }

    // Element set number
    gp.element_set_number = line1[64..68]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine1)?;

    // Inclination [degs]
    gp.inclination = line2[8..16]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Right ascension of ascending node [degs]
    gp.right_ascension_of_ascending_node = line2[17..25]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Eccentricity
    gp.eccentricity = format!("0.{}", line2[26..33].trim())
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Argument of perigee [degs]
    gp.argument_of_perigee = line2[34..42]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Mean anomaly [degs]
    gp.mean_anomaly = line2[43..51]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Mean motion [revs/day]
    gp.mean_motion = line2[52..63]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Revolution number at epoch
    gp.revolution_number_at_epoch = line2[63..68]
        .trim()
        .parse()
        .map_err(|_| GpError::InvalidTleLine2)?;

    // Initialize the SGP4 parameters
    init_sgp4(&gp, None).map_err(GpError::Sgp4)
}

/// Builds a vector of [`Sgp4`] structs from a string containing Two-Line Element (TLE) sets.
///
/// This function parses a string containing one or more TLEs in either
/// 2-line or 3-line (name + 2 lines) format and returns all successfully
/// parsed entries.
///
/// # Arguments
/// * `tle_string` - A string containing one or more Two-Line Element (TLE) sets
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 parameters
/// * `Err(GpError)` - If a TLE in the string cannot be parsed or SGP4 initialization fails
///
/// # Errors
/// * [`GpError::InvalidTleLine0`] - If an optional name line is empty or longer than 24 characters
/// * [`GpError::InvalidTleLine1`] - If a line 1 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleLine2`] - If a line 2 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleEpoch`] - If an epoch cannot be converted to a UTC datetime
/// * [`GpError::MismatchedTleCatalog`] - If a TLE has different NORAD catalog numbers on line 1 and line 2
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_string;
///
/// // Define one or more 2-line or 3-line TLEs, separated by newlines
/// let tle = "\
/// ISS (ZARYA)
/// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
/// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
/// ";
///
/// // Parse the TLE string into SGP4 propagators
/// let sgp4s = from_tle_string(tle).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 25544);
/// assert_eq!(sgp4s[0].gp.common_name, "ISS (ZARYA)");
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_string(tle_string: &str) -> Result<Vec<Sgp4>, GpError> {
    // Parse the string into lines, removing spaces
    let lines: Vec<&str> = tle_string
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect();

    // Create the tles vector
    let mut sgp4s = Vec::new();
    let mut i = 0;

    // Iterate through the lines building TLE structs when possible
    while i < lines.len() {
        // Find TLEs within the string, either 2 or 3 line entries
        if lines[i].starts_with('1') {
            // This is likely a 2 line entry, check length
            if i + 1 >= lines.len() {
                break;
            }
            // Check that next line starts with '2'
            if lines[i + 1].starts_with('2') {
                let sgp4 = from_tle_lines(lines[i], lines[i + 1], None)?;
                sgp4s.push(sgp4);
                i += 2;
            } else {
                i += 1;
            }
        } else {
            // This is likely a 3 line entry, check length
            if i + 2 >= lines.len() {
                break;
            }
            // Check that next line 2 lines starts with '1' and '2'
            if lines[i + 1].starts_with('1') && lines[i + 2].starts_with('2') {
                let sgp4 = from_tle_lines(lines[i + 1], lines[i + 2], Some(lines[i]))?;
                sgp4s.push(sgp4);
                i += 3;
            } else {
                i += 1;
            }
        }
    }
    // Return vector of SGP4 structs
    Ok(sgp4s)
}

/// Builds a vector of [`Sgp4`] structs from a file containing Two-Line Element (TLE) sets.
///
/// This function parses a file containing one or more TLEs in either
/// 2-line or 3-line (name + 2 lines) format and returns all successfully
/// parsed entries.
///
/// # Arguments
/// * `file_path` - A path to a file containing one or more Two-Line Element (TLE) sets
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 structs
/// * `Err(GpError)` - If the file cannot be read, a TLE cannot be parsed, or SGP4 initialization fails
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be read
/// * [`GpError::InvalidTleLine0`] - If an optional name line is empty or longer than 24 characters
/// * [`GpError::InvalidTleLine1`] - If a line 1 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleLine2`] - If a line 2 is not 69 characters or a field cannot be parsed
/// * [`GpError::InvalidTleEpoch`] - If an epoch cannot be converted to a UTC datetime
/// * [`GpError::MismatchedTleCatalog`] - If a TLE has different NORAD catalog numbers on line 1 and line 2
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_file;
///
/// // Define a 3-line TLE
/// let tle = "\
/// ISS (ZARYA)
/// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
/// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
/// ";
///
/// // Write the TLE to a temporary file
/// let path = std::env::temp_dir().join("mako_sgp4_from_tle_file.tle");
/// std::fs::write(&path, tle).unwrap();
///
/// // Parse the TLE file into SGP4 propagators
/// let sgp4s = from_tle_file(path.to_str().unwrap()).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 25544);
/// assert_eq!(sgp4s[0].gp.common_name, "ISS (ZARYA)");
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
pub fn from_tle_file(file_path: &str) -> Result<Vec<Sgp4>, GpError> {
    // Open the TLE file
    let tle_string = fs::read_to_string(file_path).map_err(|err| GpError::Io(err.to_string()))?;

    // Parse tle string into a vector of SGP4 structs
    from_tle_string(&tle_string)
}

/// Calculate the checksum of the Two-Line Element (TLE) line.
///
/// Given a TLE line, calculate the checksum of that line. Follow the following rules:
/// - Ignore alpha characters
/// - Sum digits 0-9 as integer values
/// - '-' is treated as 1
/// - Return checksum % 10
///
/// # Arguments
/// * `line` - The TLE line to calculate the checksum of
///
/// # Returns
/// * `checksum` - The checksum of the TLE line (integer 0-9)
fn calc_checksum(line: &str) -> i32 {
    // Initialize checksum to 0
    let mut checksum = 0;

    // Loop through the line and calculate the checksum
    for c in line.chars().take(68) {
        match c {
            '0'..='9' => checksum += (c as u8 - b'0') as i32,
            '-' => checksum += 1,
            _ => {}
        }
    }

    // Calculate the checksum
    checksum %= 10;

    checksum
}

/// Check if the TLE line has been corrupted by running a checksum test.
///
/// Given a TLE line, check if the checksum of that line is valid.
///
/// # Arguments
/// * `line` - The TLE line to check the checksum of
///
/// # Returns
/// * `bool` - True if the checksum of the line is valid, false if otherwise
fn tle_checksum(line: &str) -> bool {
    // The checksum digit is the 69th character
    if line.len() < 69 {
        return false;
    }

    // Compare the checksum to the last character of the line
    let Ok(digit) = line[68..69].parse::<i32>() else {
        return false;
    };
    calc_checksum(line) == digit
}

/// Convert a two-digit TLE year to a four-digit Gregorian year.
///
/// TLE years follow the convention that 57-99 map to 1957-1999 and 00-56 map
/// to 2000-2056, matching the original SGP4 epoch encoding.
///
/// # Arguments
/// * `two_digit_year` - Year modulo 100 as stored in a TLE
///
/// # Returns
/// * Four-digit year in the range 1957-2056
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_tle_lines;
///
/// // Two-digit TLE year 08 maps to 2008 (00-56 -> 2000-2056)
/// let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
/// let tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let sgp4 = from_tle_lines(tle_line1, tle_line2, None).unwrap();
///
/// // Assert the epoch year is the four-digit Gregorian year
/// assert_eq!(sgp4.gp.epoch_datetime.year, 2008);
/// ```
fn tle_full_year(two_digit_year: i32) -> i32 {
    if two_digit_year < 57 {
        2000 + two_digit_year
    } else {
        1900 + two_digit_year
    }
}

/// Parse a TLE NORAD catalog field (numeric or Alpha-5).
///
/// # Arguments
/// * `field` - Five-character catalog substring from TLE line 1 or line 2
///
/// # Returns
/// * `Some(n)` - The catalog number
/// * `None` - If the field is empty or not a valid catalog number
fn parse_tle_catalog(field: &str) -> Option<i32> {
    let catalog = field.trim();
    let first_char = catalog.chars().next()?;
    if first_char.is_ascii_digit() {
        // This is a classic numeric catalog number
        catalog.parse().ok()
    } else if first_char.is_ascii_alphabetic() {
        // This is an Alpha-5 catalog number
        let alpha = alpha5_digit(first_char)?;
        let rest: i32 = catalog[1..].parse().ok()?;
        Some(alpha * 10000 + rest)
    } else {
        None
    }
}

/// Converts an Alpha-5 digit to an integer
///
/// Given a character, map it to the corresponding integer value as defined in the Alpha-5 standard.
///
/// # Arguments
/// * `c` - The character to convert
///
/// # Returns
/// * `Some(i)` - The integer value of the character
/// * `None` - If the character is not a valid Alpha-5 digit
///
/// # References
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
fn alpha5_digit(c: char) -> Option<i32> {
    // Convert the character to uppercase to match the standard
    let c_upper = c.to_ascii_uppercase();

    // Match the character to the corresponding integer value
    match c_upper {
        'A'..='H' => Some((c_upper as i32) - ('A' as i32) + 10),
        'J'..='N' => Some((c_upper as i32) - ('J' as i32) + 18),
        'P'..='Z' => Some((c_upper as i32) - ('P' as i32) + 23),
        _ => None, // Standard skips 'I' and 'O'
    }
}

/// Converts an integer to an Alpha-5 digit
///
/// Inverse of [`alpha5_digit`]. Integers 10-17 map to A-H, 18-22 to J-N,
/// and 23-33 to P-Z. The letters I and O are skipped.
///
/// # Arguments
/// * `digit` - The integer value to convert
///
/// # Returns
/// * `Some(c)` - The Alpha-5 character
/// * `None` - If the integer is outside 10-33 or would map to I or O
///
/// # References
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
fn alpha5_letter(digit: i32) -> Option<char> {
    // Match the integer to the corresponding Alpha-5 character
    match digit {
        10..=17 => Some((b'A' + (digit as u8 - 10)) as char),
        18..=22 => Some((b'J' + (digit as u8 - 18)) as char),
        23..=33 => Some((b'P' + (digit as u8 - 23)) as char),
        _ => None, // Standard skips 'I' and 'O'
    }
}

/// Format a NORAD catalog number for a TLE
///
/// Numbers below 100000 are written as a 5-character right-justified
/// field. Numbers 100000-339999 use Alpha-5 encoding.
///
/// # Arguments
/// * `catalog` - The NORAD catalog number
///
/// # Returns
/// * `Ok(String)` - A 5-character catalog field
/// * `Err(GpError)` - If the catalog number cannot be written as a TLE field
///
/// # Errors
/// * [`GpError::InvalidTLECatalogNumber`] - If the catalog number is negative or greater than 339999
///
/// # References
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
fn format_tle_catalog_number(catalog: i32) -> Result<String, GpError> {
    // Reject negative catalog numbers
    if catalog < 0 {
        return Err(GpError::InvalidTLECatalogNumber);
    }

    if catalog < 100000 {
        // This is a classic numeric catalog number
        return Ok(format!("{:>5}", catalog));
    }

    // This is an Alpha-5 catalog number
    let first = catalog / 10000;
    let rest = catalog % 10000;
    let Some(letter) = alpha5_letter(first) else {
        return Err(GpError::InvalidTLECatalogNumber);
    };
    Ok(format!("{}{:04}", letter, rest))
}

/// Format an international designator for a TLE
///
/// Converts YYYY-NNNP form back to the 8-character TLE field YYNNNPPP.
/// Missing designators are written as eight spaces. The TLE field allows
/// a 3-digit launch number and a 3-character piece; longer values stored
/// on the GP cannot be written as a TLE.
///
/// # Arguments
/// * `intl` - The international designator stored on the GP
///
/// # Returns
/// * `Ok(String)` - An 8-character international designator field
/// * `Err(GpError)` - If the designator cannot be written as a TLE field
///
/// # Errors
/// * [`GpError::InvalidTLEInternationalDesignator`] - If the launch number or piece exceeds 3 characters
fn format_tle_intl_des(intl: &str) -> Result<String, GpError> {
    let intl = intl.trim();
    if intl.is_empty() {
        // Handle case where international designator is not present
        return Ok("        ".to_string());
    }

    // Collapse YYYY-NNNP back to the two-digit TLE year plus the launch piece
    let field = if let Some((year_str, rest)) = intl.split_once('-') {
        // TLE allows 3 launch-number digits and 3 piece characters
        let launch_len = rest.bytes().take_while(u8::is_ascii_digit).count();
        let part_code = &rest[launch_len..];
        if launch_len > 3 || part_code.len() > 3 {
            return Err(GpError::InvalidTLEInternationalDesignator);
        }

        let year = year_str.parse::<i32>().unwrap_or(0);
        let yy = ((year % 100) + 100) % 100;
        format!("{:02}{}", yy, rest)
    } else if intl.len() > 8 {
        return Err(GpError::InvalidTLEInternationalDesignator);
    } else {
        intl.to_string()
    };

    // TLE international designator field is 8 characters, left-justified
    Ok(format!("{:<8}", field))
}

/// Format a TLE epoch from a UTC datetime
///
/// Writes YYDDD.DDDDDDDD: two-digit year and day of year with eight
/// fractional digits.
///
/// # Arguments
/// * `epoch` - The epoch datetime
///
/// # Returns
/// * `Ok(String)` - A 14-character TLE epoch field
/// * `Err(GpError)` - If the datetime cannot be written as a TLE epoch
///
/// # Errors
/// * [`GpError::InvalidTLEDateTime`] - If the year is before 1957 or after 2056
fn format_tle_epoch(epoch: &DateTime) -> Result<String, GpError> {
    // TLE two-digit years cover 1957-2056
    if !(1957..=2056).contains(&epoch.year) {
        return Err(GpError::InvalidTLEDateTime);
    }

    // TLE years are stored modulo 100
    let two_digit_year = ((epoch.year % 100) + 100) % 100;

    // Check for leap year
    let is_leap = (epoch.year % 4 == 0 && epoch.year % 100 != 0) || (epoch.year % 400 == 0);

    // Days per month (non-leap year)
    let days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    // Sum complete months to get the integer day of year
    let mut day_int = epoch.day;
    for m in 1..epoch.month {
        let days_this_month = if m == 2 && is_leap {
            29 // February in leap year
        } else {
            days_per_month[(m - 1) as usize]
        };
        day_int += days_this_month;
    }

    // Convert clock time to a fractional day
    let mut second = epoch.second;
    if second < 0.0 {
        second = 0.0;
    }
    let day_frac = (epoch.hour as f64 * 3600.0 + epoch.minute as f64 * 60.0 + second) / 86400.0;
    let dayofyr = day_int as f64 + day_frac;

    Ok(format!("{:02}{:012.8}", two_digit_year, dayofyr))
}

/// Format the TLE-printed first derivative of mean motion
///
/// Writes the 10-character field as a signed decimal without a leading
/// zero, for example `-.00002182` or ` .00000234`.
///
/// # Arguments
/// * `value` - The TLE-printed first derivative (stored value divided by 2) \[revs/day^2\]
///
/// # Returns
/// * A 10-character mean-motion-dot field
fn format_tle_ndot(value: f64) -> String {
    // Eight digits after the implied 0. prefix
    let mut digits = (value.abs() * 1e8).round() as i64;
    if digits > 99999999 {
        digits = 99999999;
    }

    // Negative values use '-', positive values use a leading space
    if value < 0.0 && digits != 0 {
        return format!("-.{:08}", digits);
    }
    format!(" .{:08}", digits)
}

/// Format a TLE exponential field (second derivative of mean motion or BSTAR)
///
/// Writes the 8-character assumed-decimal form, for example `-11606-4`
/// meaning -0.11606 times 10 to the -4. Zero is written as ` 00000+0`.
///
/// # Arguments
/// * `value` - The TLE-printed value
///
/// # Returns
/// * An 8-character TLE exponential field
fn format_tle_exp(value: f64) -> String {
    if !value.is_finite() || value == 0.0 {
        return " 00000+0".to_string();
    }

    let sign = if value < 0.0 { '-' } else { ' ' };
    let abs = value.abs();

    // Mantissa in [0.1, 1) so the five digits sit after an implied decimal
    let mut exp = abs.log10().floor() as i32 + 1;
    let mantissa = abs / 10f64.powi(exp);
    let mut digits = (mantissa * 100000.0).round() as i32;

    if digits >= 100000 {
        digits = 10000;
        exp += 1;
    }
    if digits == 0 {
        return " 00000+0".to_string();
    }

    // TLE exponent is a single digit
    if exp > 9 {
        exp = 9;
        digits = 99999;
    } else if exp < -9 {
        return " 00000+0".to_string();
    }

    let exp_sign = if exp >= 0 { '+' } else { '-' };
    format!("{}{:05}{}{}", sign, digits, exp_sign, exp.abs())
}

/// Format eccentricity for a TLE
///
/// Writes seven digits with the decimal point assumed. For example,
/// 0.0006703 becomes `0006703`.
///
/// # Arguments
/// * `eccentricity` - Orbital eccentricity
///
/// # Returns
/// * A 7-character eccentricity field
fn format_tle_eccentricity(eccentricity: f64) -> String {
    // Seven digits after the implied 0. prefix
    let mut digits = (eccentricity.abs() * 1e7).round() as i64;
    digits = digits.clamp(0, 9999999);
    format!("{:07}", digits)
}

/// Append a TLE checksum digit to a 68-character line
///
/// Computes the checksum of the first 68 characters with [`calc_checksum`]
/// and appends it so the returned line is 69 characters.
///
/// # Arguments
/// * `line68` - The first 68 characters of a TLE data line
///
/// # Returns
/// * `Ok(String)` - A 69-character TLE line including checksum
/// * `Err(GpError)` - If `line68` is not exactly 68 characters
///
/// # Errors
/// * [`GpError::InvalidTLELine`] - If `line68` is not exactly 68 characters
fn tle_line_with_checksum(line68: &str) -> Result<String, GpError> {
    if line68.len() != 68 {
        return Err(GpError::InvalidTLELine);
    }

    // Append checksum as the 69th character
    Ok(format!("{}{}", line68, calc_checksum(line68)))
}

/// Build TLE line 1 from a general perturbation element set
///
/// MEAN_MOTION_DOT and MEAN_MOTION_DDOT are written as the TLE-printed
/// values (stored derivatives divided by 2 and 6).
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * `Ok(String)` - TLE line 1, including checksum
/// * `Err(GpError)` - If a TLE field cannot be written
fn format_tle_line1(gp: &GenPerturbElementSet) -> Result<String, GpError> {
    // Line 1
    // Satellite catalog number
    let catalog = format_tle_catalog_number(gp.satellite_catalog_number)?;

    // Classification
    let classification = if gp.classification.is_ascii_graphic() {
        gp.classification
    } else {
        'U'
    };

    // International designator
    let intl_des = format_tle_intl_des(&gp.international_designator)?;

    // Epoch year and day of year
    let epoch = format_tle_epoch(&gp.epoch_datetime)?;

    // 1st derivative of mean motion (TLE prints the stored value divided by 2)
    let ndot = format_tle_ndot(gp.first_derivative_of_mean_motion / 2.0);

    // 2nd derivative of mean motion (TLE prints the stored value divided by 6)
    let nddot = format_tle_exp(gp.second_derivative_of_mean_motion / 6.0);

    // BSTAR drag term
    let bstar = format_tle_exp(gp.bstar);

    // Ephemeris type
    let ephem_type = if (0..=9).contains(&gp.ephemeris_type) {
        gp.ephemeris_type
    } else {
        0
    };

    // Element set number (4 digits)
    let elset = gp.element_set_number.rem_euclid(10000);

    let line68 = format!(
        "1 {}{} {} {} {} {} {} {} {:>4}",
        catalog, classification, intl_des, epoch, ndot, nddot, bstar, ephem_type, elset
    );
    tle_line_with_checksum(&line68)
}

/// Build TLE line 2 from a general perturbation element set
///
/// Writes inclination, RAAN, eccentricity, argument of perigee, mean
/// anomaly, mean motion, and revolution number in TLE column layout.
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * `Ok(String)` - TLE line 2, including checksum
/// * `Err(GpError)` - If a TLE field cannot be written
fn format_tle_line2(gp: &GenPerturbElementSet) -> Result<String, GpError> {
    // Line 2
    // Satellite catalog number
    let catalog = format_tle_catalog_number(gp.satellite_catalog_number)?;

    // Inclination [degs]
    let inclination = format!("{:8.4}", gp.inclination);

    // Right ascension of ascending node [degs]
    let right_ascension_of_ascending_node = format!("{:8.4}", gp.right_ascension_of_ascending_node);

    // Eccentricity (decimal point assumed)
    let ecc = format_tle_eccentricity(gp.eccentricity);

    // Argument of perigee [degs]
    let argument_of_perigee = format!("{:8.4}", gp.argument_of_perigee);

    // Mean anomaly [degs]
    let mean_anomaly = format!("{:8.4}", gp.mean_anomaly);

    // Mean motion [revs/day]
    let mean_motion = format!("{:11.8}", gp.mean_motion);

    // Revolution number at epoch (5 digits)
    let rev = format!("{:>5}", gp.revolution_number_at_epoch.rem_euclid(100000));

    let line68 = format!(
        "2 {} {} {} {} {} {} {}{}",
        catalog,
        inclination,
        right_ascension_of_ascending_node,
        ecc,
        argument_of_perigee,
        mean_anomaly,
        mean_motion,
        rev
    );
    tle_line_with_checksum(&line68)
}

/// Build one TLE from a general perturbation element set
///
/// A name line (line 0) is written when [`GenPerturbElementSet::common_name`]
/// is non-empty. Names longer than 24 characters are truncated. MEAN_MOTION
/// derivatives are written as the TLE-printed values (stored derivatives
/// divided by 2 and 6).
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * `Ok(String)` - One TLE in 2-line or 3-line format, including a trailing newline
/// * `Err(GpError)` - If a TLE field cannot be written
fn gp_to_tle(gp: &GenPerturbElementSet) -> Result<String, GpError> {
    let mut tle = String::new();

    // Optional name line (TLE line 0)
    let name = gp.common_name.trim();
    if !name.is_empty() {
        let mut name = name.to_string();
        if name.len() > 24 {
            name.truncate(24);
        }
        tle.push_str(&name);
        tle.push('\n');
    }

    // TLE line 1 and line 2
    tle.push_str(&format_tle_line1(gp)?);
    tle.push('\n');
    tle.push_str(&format_tle_line2(gp)?);
    tle.push('\n');
    Ok(tle)
}

/// Builds a Two-Line Element (TLE) string from a slice of [`Sgp4`] structs.
///
/// Each element set is written as a 2-line TLE, or a 3-line TLE when a
/// common name is present. Catalog numbers of 100000 and above use
/// Alpha-5 encoding. MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the
/// TLE-printed values (the stored derivatives divided by 2 and 6).
///
/// An empty slice produces an empty string.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
///
/// # Returns
/// * `Ok(String)` - One or more TLEs
/// * `Err(GpError)` - If a TLE field cannot be written
///
/// # Errors
/// * [`GpError::InvalidTLECatalogNumber`] - If a catalog number is negative or greater than 339999
/// * [`GpError::InvalidTLEDateTime`] - If an epoch year is before 1957 or after 2056
/// * [`GpError::InvalidTLEInternationalDesignator`] - If a launch number or piece exceeds 3 characters
/// * [`GpError::InvalidTLELine`] - If a formatted TLE data line is not 68 characters before the checksum
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_tle_string, to_tle_string};
///
/// // Define a 3-line TLE
/// let tle = "\
/// ISS (ZARYA)
/// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
/// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
/// ";
///
/// // Parse, export, and parse again
/// let sgp4s = from_tle_string(tle).unwrap();
/// let exported = to_tle_string(&sgp4s).unwrap();
/// let reparsed = from_tle_string(&exported).unwrap();
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 25544);
/// assert_eq!(reparsed[0].gp.common_name, "ISS (ZARYA)");
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
pub fn to_tle_string(sgp4s: &[Sgp4]) -> Result<String, GpError> {
    // Serialize each element set as one TLE
    let mut records = String::new();
    for sgp4 in sgp4s {
        records.push_str(&gp_to_tle(&sgp4.gp)?);
    }

    Ok(records)
}

/// Writes a Two-Line Element (TLE) file from a slice of [`Sgp4`] structs.
///
/// This function serializes the provided element sets with
/// [`to_tle_string`] and writes the result to `tle_file_path`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
/// * `tle_file_path` - Destination path for the TLE file
///
/// # Returns
/// * `Ok(())` - If the file was written
/// * `Err(GpError)` - If a TLE field cannot be written or the file cannot be written
///
/// # Errors
/// * [`GpError::InvalidTLECatalogNumber`] - If a catalog number is negative or greater than 339999
/// * [`GpError::InvalidTLEDateTime`] - If an epoch year is before 1957 or after 2056
/// * [`GpError::InvalidTLEInternationalDesignator`] - If a launch number or piece exceeds 3 characters
/// * [`GpError::InvalidTLELine`] - If a formatted TLE data line is not 68 characters before the checksum
/// * [`GpError::Io`] - If the file cannot be written
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_tle_file, from_tle_string, to_tle_file};
///
/// // Define a 3-line TLE and parse it
/// let tle = "\
/// ISS (ZARYA)
/// 1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921
/// 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
/// ";
/// let sgp4s = from_tle_string(tle).unwrap();
///
/// // Write the element sets to a temporary TLE file
/// let path = std::env::temp_dir().join("mako_sgp4_to_tle_file.tle");
/// let path_str = path.to_str().unwrap();
/// to_tle_file(&sgp4s, path_str).unwrap();
///
/// // Parse the written file
/// let reparsed = from_tle_file(path_str).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 25544);
/// assert_eq!(reparsed[0].gp.common_name, "ISS (ZARYA)");
/// ```
///
/// # References
/// - [Celestrak TLE Format](https://celestrak.org/columns/v04n03/#FAQ01)
/// - [Alpha-5 Standard](https://www.space-track.org/documentation#tle-alpha5)
pub fn to_tle_file(sgp4s: &[Sgp4], tle_file_path: &str) -> Result<(), GpError> {
    // Serialize the element sets and write the TLE file
    let tle_string = to_tle_string(sgp4s)?;
    fs::write(tle_file_path, tle_string).map_err(|err| GpError::Io(err.to_string()))?;
    Ok(())
}

/// Builds a [`Sgp4`] struct from the lines of a single OMM KVN record.
///
/// Given the key-value lines of one Orbit Mean-Elements Message, this function
/// parses the input into a [`GenPerturbElementSet`] and initializes SGP4,
/// matching [`from_tle_lines`]. Missing numeric fields default to 0 and missing
/// string fields default to an empty string.
///
/// # Arguments
/// * `lines` - The KVN key-value lines of one OMM record
///
/// # Returns
/// * `Ok(Sgp4)` - Struct containing the parsed SGP4 parameters
/// * `Err(GpError)` - If a field cannot be parsed or SGP4 initialization fails
///
/// # Errors
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_lines;
///
/// // Define one OMM record as key-value lines (missing fields use defaults)
/// let omm = "\
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
///
/// // Split the record into trimmed, non-empty lines
/// let lines: Vec<&str> = omm.lines().map(str::trim).filter(|l| !l.is_empty()).collect();
///
/// // Parse the OMM lines into an SGP4 propagator
/// let sgp4 = from_omm_kvn_lines(&lines).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4.gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4.gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_lines(lines: &[&str]) -> Result<Sgp4, GpError> {
    sgp4_from_omm_lookup(|field| kvn_lookup(lines, field))
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS KVN (key-value notation) format.
///
/// This function parses a string containing one or more OMM records. A new record
/// starts at each CCSDS_OMM_VERS line. Each record is parsed by
/// [`from_omm_kvn_lines`].
///
/// # Arguments
/// * `omm_kvn_string` - A string containing one or more OMM records in KVN format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 parameters
/// * `Err(GpError)` - If a field cannot be parsed or SGP4 initialization fails
///
/// # Errors
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_string;
///
/// // Define an OMM KVN record. Concatenated records start at each CCSDS_OMM_VERS line
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// EPHEMERIS_TYPE = 0
/// CLASSIFICATION_TYPE = U
/// NORAD_CAT_ID   = 69097
/// ELEMENT_SET_NO = 999
/// REV_AT_EPOCH   = 459
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
///
/// // Parse the OMM string into SGP4 propagators
/// let sgp4s = from_omm_kvn_string(omm).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_string(omm_kvn_string: &str) -> Result<Vec<Sgp4>, GpError> {
    // Split the string into individual OMM records
    let mut sgp4s = Vec::new();
    let mut current_lines: Vec<&str> = Vec::new();

    // Iterate through the lines, starting a new record at each CCSDS_OMM_VERS
    for line in omm_kvn_string.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // A new CCSDS_OMM_VERS line starts the next concatenated record
        if trimmed.to_ascii_uppercase().starts_with("CCSDS_OMM_VERS") && !current_lines.is_empty() {
            let sgp4 = from_omm_kvn_lines(&current_lines)?;
            sgp4s.push(sgp4);
            current_lines.clear();
        }

        current_lines.push(trimmed);
    }

    // Parse the final record
    if !current_lines.is_empty() {
        let sgp4 = from_omm_kvn_lines(&current_lines)?;
        sgp4s.push(sgp4);
    }

    // Return vector of SGP4 structs
    Ok(sgp4s)
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS KVN format.
///
/// This function parses a file containing one or more OMM records in KVN
/// format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_kvn_file_path` - A path to a file containing one or more OMM records in KVN format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 structs
/// * `Err(GpError)` - If the file cannot be read or an OMM cannot be parsed
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be read
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_kvn_file;
///
/// // Define an OMM KVN record
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
///
/// // Write the OMM to a temporary file
/// let path = std::env::temp_dir().join("mako_sgp4_from_omm_kvn_file.txt");
/// std::fs::write(&path, omm).unwrap();
///
/// // Parse the OMM file into SGP4 propagators
/// let sgp4s = from_omm_kvn_file(path.to_str().unwrap()).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn from_omm_kvn_file(omm_kvn_file_path: &str) -> Result<Vec<Sgp4>, GpError> {
    // Open the OMM KVN file
    let omm_kvn_string =
        fs::read_to_string(omm_kvn_file_path).map_err(|err| GpError::Io(err.to_string()))?;

    // Parse OMM string into a vector of SGP4 structs
    from_omm_kvn_string(&omm_kvn_string)
}

/// Parse a KVN field from OMM record lines
///
/// Search the provided lines for a key-value entry matching the given field
/// name. Missing or empty fields return a type-specific default: 0 for numbers
/// and an empty string for strings. Classification defaults to U.
///
/// # Arguments
/// * `lines` - The KVN lines of one OMM record
/// * `field` - The KVN keyword to look up
///
/// # Returns
/// * `Ok(T)` - The parsed value as the requested type
/// * `Err(GpError)` - If a present field cannot be parsed
#[cfg(test)]
fn kvn_parse<T: FromOmm>(lines: &[&str], field: &str) -> Result<T, GpError> {
    omm_typed_value(kvn_lookup(lines, field), field)
}

/// Look up a KVN keyword and return the cleaned value, if present
///
/// # Arguments
/// * `lines` - The KVN lines of one OMM record
/// * `field` - The KVN keyword to look up
///
/// # Returns
/// * `Some(value)` - The cleaned field value
/// * `None` - If the field is not present
fn kvn_lookup(lines: &[&str], field: &str) -> Option<String> {
    // Search each line for a matching KVN keyword
    for line in lines {
        let line = line.trim();

        // Skip empty lines and COMMENT lines
        if line.is_empty() || line.to_ascii_uppercase().starts_with("COMMENT") {
            continue;
        }

        // Split KEY = VALUE
        let Some((keyword, value)) = line.split_once('=') else {
            continue;
        };
        if !keyword.trim().eq_ignore_ascii_case(field) {
            continue;
        }

        return Some(clean_kvn_value(value));
    }

    None
}

/// Convert an optional OMM field string into type T
///
/// Missing or empty values use the type default.
///
/// # Arguments
/// * `value` - The raw field text, if the keyword was present
/// * `field` - The OMM keyword name, used in parse error messages
///
/// # Returns
/// * `Ok(T)` - The parsed value as the requested type
/// * `Err(GpError)` - If a present field cannot be parsed
fn omm_typed_value<T: FromOmm>(value: Option<String>, field: &str) -> Result<T, GpError> {
    match value {
        Some(v) if !v.is_empty() => T::from_omm(field, &v),
        _ => Ok(T::omm_default()),
    }
}

/// Build an [`Sgp4`] struct from an OMM field lookup
///
/// Shared by KVN and XML parsers. Missing numeric fields default to 0 and
/// missing string fields default to an empty string. OMM MEAN_MOTION_DOT and
/// MEAN_MOTION_DDOT are the TLE-printed values (already divided by 2 and 6)
/// and are scaled so [`GenPerturbElementSet`] stores the true derivatives.
///
/// # Arguments
/// * `lookup` - Returns the text of one OMM keyword, if present
///
/// # Returns
/// * `Ok(Sgp4)` - Struct containing the parsed SGP4 parameters
/// * `Err(GpError)` - If a field cannot be parsed or SGP4 initialization fails
fn sgp4_from_omm_lookup<F>(lookup: F) -> Result<Sgp4, GpError>
where
    F: Fn(&str) -> Option<String>,
{
    // Create a General Perturbation Element Set struct
    let gp = GenPerturbElementSet {
        common_name: omm_typed_value(lookup("OBJECT_NAME"), "OBJECT_NAME")?,
        satellite_catalog_number: omm_typed_value(lookup("NORAD_CAT_ID"), "NORAD_CAT_ID")?,
        classification: omm_typed_value(lookup("CLASSIFICATION_TYPE"), "CLASSIFICATION_TYPE")?,
        international_designator: omm_typed_value(lookup("OBJECT_ID"), "OBJECT_ID")?,
        epoch_datetime: omm_typed_value(lookup("EPOCH"), "EPOCH")?,
        first_derivative_of_mean_motion: omm_typed_value::<f64>(
            lookup("MEAN_MOTION_DOT"),
            "MEAN_MOTION_DOT",
        )? * 2.0,
        second_derivative_of_mean_motion: omm_typed_value::<f64>(
            lookup("MEAN_MOTION_DDOT"),
            "MEAN_MOTION_DDOT",
        )? * 6.0,
        bstar: omm_typed_value(lookup("BSTAR"), "BSTAR")?,
        ephemeris_type: omm_typed_value(lookup("EPHEMERIS_TYPE"), "EPHEMERIS_TYPE")?,
        element_set_number: omm_typed_value(lookup("ELEMENT_SET_NO"), "ELEMENT_SET_NO")?,
        inclination: omm_typed_value(lookup("INCLINATION"), "INCLINATION")?,
        right_ascension_of_ascending_node: omm_typed_value(
            lookup("RA_OF_ASC_NODE"),
            "RA_OF_ASC_NODE",
        )?,
        eccentricity: omm_typed_value(lookup("ECCENTRICITY"), "ECCENTRICITY")?,
        argument_of_perigee: omm_typed_value(lookup("ARG_OF_PERICENTER"), "ARG_OF_PERICENTER")?,
        mean_anomaly: omm_typed_value(lookup("MEAN_ANOMALY"), "MEAN_ANOMALY")?,
        mean_motion: omm_typed_value(lookup("MEAN_MOTION"), "MEAN_MOTION")?,
        revolution_number_at_epoch: omm_typed_value(lookup("REV_AT_EPOCH"), "REV_AT_EPOCH")?,
    };

    // Initialize the SGP4 parameters
    init_sgp4(&gp, None).map_err(GpError::Sgp4)
}

/// Clean a KVN value string
///
/// Trim whitespace, strip matching double quotes, and remove optional trailing
/// units written in square brackets.
///
/// # Arguments
/// * `value` - The raw value text to the right of the equals sign
///
/// # Returns
/// * The cleaned value string
fn clean_kvn_value(value: &str) -> String {
    let mut v = value.trim();

    // Strip matching double quotes
    if v.len() >= 2 && v.starts_with('"') && v.ends_with('"') {
        v = &v[1..v.len() - 1];
        v = v.trim();
    }

    // Strip optional units in square brackets
    if let Some(idx) = v.find('[') {
        v = v[..idx].trim();
    }

    v.to_string()
}

/// Format one CCSDS KVN keyword and value line
///
/// Keywords of 14 characters or fewer are padded so the equals sign aligns
/// with Celestrak OMM KVN files.
///
/// # Arguments
/// * `keyword` - The KVN keyword
/// * `value` - The KVN value
///
/// # Returns
/// * One `KEYWORD = VALUE` line
fn format_kvn_line(keyword: &str, value: &str) -> String {
    if keyword.len() >= 14 {
        format!("{} = {}", keyword, value)
    } else {
        format!("{:<14} = {}", keyword, value)
    }
}

/// Trim trailing zeros from a decimal string
///
/// # Arguments
/// * `value` - A decimal number already formatted as text
///
/// # Returns
/// * The number text without trailing fractional zeros
fn trim_decimal_zeros(value: &str) -> String {
    if !value.contains('.') {
        return value.to_string();
    }

    let trimmed = value.trim_end_matches('0');
    if trimmed.ends_with('.') {
        return trimmed.trim_end_matches('.').to_string();
    }
    trimmed.to_string()
}

/// Strip a leading zero before a decimal point
///
/// Celestrak OMM text uses `.00147468` rather than `0.00147468`. Whole
/// numbers such as `0` are left unchanged.
///
/// # Arguments
/// * `value` - A decimal number already formatted as text
///
/// # Returns
/// * The number text without a leading zero before the decimal point
fn strip_leading_decimal_zero(value: &str) -> String {
    if let Some(rest) = value.strip_prefix("-0.") {
        return format!("-.{}", rest);
    }
    if let Some(rest) = value.strip_prefix("0.") {
        return format!(".{}", rest);
    }
    value.to_string()
}

/// Format an OMM decimal number
///
/// Rounds to 12 decimal places to suppress binary float noise, then drops
/// trailing zeros. Values between -1 and 1 are written without a leading
/// zero (`.###`). Zero is written as `0`.
///
/// # Arguments
/// * `value` - The number to format
///
/// # Returns
/// * Decimal text suitable for a KVN, XML, or CSV value
fn format_omm_decimal(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }

    let rounded = (value * 1e12).round() / 1e12;
    strip_leading_decimal_zero(&trim_decimal_zeros(&format!("{:.12}", rounded)))
}

/// Format an OMM scientific-notation number
///
/// Matches Celestrak style: a leading decimal mantissa and an explicit
/// exponent sign, for example `.39221734E-3`. Zero is written as `0`.
///
/// # Arguments
/// * `value` - The number to format
///
/// # Returns
/// * Scientific-notation text suitable for BSTAR, MEAN_MOTION_DOT, and
///   MEAN_MOTION_DDOT
fn format_omm_sci(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }

    let sign = if value < 0.0 { "-" } else { "" };
    let abs = value.abs();

    // Mantissa in [0.1, 1) so the leading digit sits after the decimal point
    let mut exp = abs.log10().floor() as i32 + 1;
    let mut mantissa = abs / 10f64.powi(exp);

    if mantissa >= 1.0 {
        mantissa /= 10.0;
        exp += 1;
    } else if mantissa < 0.1 {
        mantissa *= 10.0;
        exp -= 1;
    }

    // 10 digits after the decimal, then drop trailing zeros
    const SCALE: f64 = 10_000_000_000.0;
    let mut digits = (mantissa * SCALE).round() as i64;
    if digits >= 10_000_000_000 {
        exp += 1;
        digits = 1_000_000_000;
    }

    let mut frac = format!("{:010}", digits);
    frac = frac.trim_end_matches('0').to_string();
    if frac.is_empty() {
        frac = "0".to_string();
    }

    let exp_sign = if exp >= 0 { "+" } else { "-" };
    format!("{}.{}E{}{}", sign, frac, exp_sign, exp.abs())
}

/// Format an OMM EPOCH string from a UTC datetime
///
/// Writes `YYYY-MM-DDTHH:MM:SS.ssssss` with six fractional-second digits.
///
/// # Arguments
/// * `epoch` - The epoch datetime
///
/// # Returns
/// * The EPOCH value in CCSDS calendar form
fn format_omm_epoch(epoch: &DateTime) -> String {
    let mut second = epoch.second;
    if second < 0.0 {
        second = 0.0;
    }
    if second >= 60.0 {
        second = 59.999999;
    }

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:09.6}",
        epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, second
    )
}

/// Format a classification character for OMM export
///
/// Non-printable values default to unclassified (`U`).
///
/// # Arguments
/// * `classification` - The GP classification character
///
/// # Returns
/// * A single-character CLASSIFICATION_TYPE value
fn format_omm_classification(classification: char) -> String {
    if classification.is_ascii_graphic() {
        return classification.to_string();
    }
    "U".to_string()
}

/// Build one OMM KVN record from a general perturbation element set
///
/// Metadata that is not stored on [`GenPerturbElementSet`] is filled with
/// SGP4 defaults. MEAN_MOTION_DOT and MEAN_MOTION_DDOT are written as the
/// TLE-printed values (stored derivatives divided by 2 and 6).
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * One OMM record in KVN format, including a trailing newline
fn gp_to_omm_kvn(gp: &GenPerturbElementSet) -> String {
    // OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
    let mean_motion_dot = gp.first_derivative_of_mean_motion / 2.0;
    let mean_motion_ddot = gp.second_derivative_of_mean_motion / 6.0;

    vec![
        // Header
        format_kvn_line("CCSDS_OMM_VERS", "2.0"),
        format_kvn_line("CREATION_DATE", ""),
        format_kvn_line("ORIGINATOR", ""),
        String::new(),
        // Metadata
        format_kvn_line("OBJECT_NAME", &gp.common_name),
        format_kvn_line("OBJECT_ID", &gp.international_designator),
        format_kvn_line("CENTER_NAME", "EARTH"),
        format_kvn_line("REF_FRAME", "TEME"),
        format_kvn_line("TIME_SYSTEM", "UTC"),
        format_kvn_line("MEAN_ELEMENT_THEORY", "SGP/SGP4"),
        String::new(),
        // Mean elements
        format_kvn_line("EPOCH", &format_omm_epoch(&gp.epoch_datetime)),
        format_kvn_line("MEAN_MOTION", &format_omm_decimal(gp.mean_motion)),
        format_kvn_line("ECCENTRICITY", &format_omm_decimal(gp.eccentricity)),
        format_kvn_line("INCLINATION", &format_omm_decimal(gp.inclination)),
        format_kvn_line(
            "RA_OF_ASC_NODE",
            &format_omm_decimal(gp.right_ascension_of_ascending_node),
        ),
        format_kvn_line(
            "ARG_OF_PERICENTER",
            &format_omm_decimal(gp.argument_of_perigee),
        ),
        format_kvn_line("MEAN_ANOMALY", &format_omm_decimal(gp.mean_anomaly)),
        String::new(),
        // TLE / SGP4 parameters
        format_kvn_line("EPHEMERIS_TYPE", &gp.ephemeris_type.to_string()),
        format_kvn_line(
            "CLASSIFICATION_TYPE",
            &format_omm_classification(gp.classification),
        ),
        format_kvn_line("NORAD_CAT_ID", &gp.satellite_catalog_number.to_string()),
        format_kvn_line("ELEMENT_SET_NO", &gp.element_set_number.to_string()),
        format_kvn_line("REV_AT_EPOCH", &gp.revolution_number_at_epoch.to_string()),
        format_kvn_line("BSTAR", &format_omm_sci(gp.bstar)),
        format_kvn_line("MEAN_MOTION_DOT", &format_omm_sci(mean_motion_dot)),
        format_kvn_line("MEAN_MOTION_DDOT", &format_omm_sci(mean_motion_ddot)),
        // Trailing newline so concatenated records split on CCSDS_OMM_VERS
        String::new(),
    ]
    .join("\n")
}

/// Builds a CCSDS KVN OMM string from a slice of [`Sgp4`] structs.
///
/// Each element set is written as one OMM record starting with
/// `CCSDS_OMM_VERS`. Metadata fields that are not stored on
/// [`GenPerturbElementSet`] are filled with SGP4 defaults: Earth-centered,
/// TEME, UTC, and SGP/SGP4. CREATION_DATE and ORIGINATOR are left empty.
/// OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
/// (the stored derivatives divided by 2 and 6).
///
/// An empty slice produces an empty string.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
///
/// # Returns
/// * `String` - One or more OMM records in KVN format
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_kvn_string, to_omm_kvn_string};
///
/// // Define an OMM KVN record
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
///
/// // Parse, export, and parse again
/// let sgp4s = from_omm_kvn_string(omm).unwrap();
/// let exported = to_omm_kvn_string(&sgp4s);
/// let reparsed = from_omm_kvn_string(&exported).unwrap();
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn to_omm_kvn_string(sgp4s: &[Sgp4]) -> String {
    // Serialize each element set as one KVN record
    let mut records = String::new();
    for sgp4 in sgp4s {
        records.push_str(&gp_to_omm_kvn(&sgp4.gp));
    }

    records
}

/// Writes a CCSDS KVN OMM file from a slice of [`Sgp4`] structs.
///
/// This function serializes the provided element sets with
/// [`to_omm_kvn_string`] and writes the result to `omm_kvn_file_path`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
/// * `omm_kvn_file_path` - Destination path for the KVN file
///
/// # Returns
/// * `Ok(())` - If the file is written
/// * `Err(GpError)` - If the file cannot be written
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be written
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_kvn_file, from_omm_kvn_string, to_omm_kvn_file};
///
/// // Define an OMM KVN record and parse it
/// let omm = "\
/// CCSDS_OMM_VERS = 2.0
/// OBJECT_NAME    = 2026-106A
/// OBJECT_ID      = 2026-106A
/// EPOCH          = 2026-06-14T15:07:48.259488
/// MEAN_MOTION    = 15.11169557
/// ECCENTRICITY   = .00147468
/// INCLINATION    = 97.5103
/// RA_OF_ASC_NODE = 247.7605
/// ARG_OF_PERICENTER = 169.6213
/// MEAN_ANOMALY   = 190.5325
/// NORAD_CAT_ID   = 69097
/// BSTAR          = .39221734E-3
/// MEAN_MOTION_DOT = .6535E-4
/// MEAN_MOTION_DDOT = 0
/// ";
/// let sgp4s = from_omm_kvn_string(omm).unwrap();
///
/// // Write the element sets to a temporary KVN file
/// let path = std::env::temp_dir().join("mako_sgp4_to_omm_kvn_file.txt");
/// let path_str = path.to_str().unwrap();
/// to_omm_kvn_file(&sgp4s, path_str).unwrap();
///
/// // Parse the written file
/// let reparsed = from_omm_kvn_file(path_str).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
pub fn to_omm_kvn_file(sgp4s: &[Sgp4], omm_kvn_file_path: &str) -> Result<(), GpError> {
    // Serialize the element sets and write the KVN file
    let omm_kvn_string = to_omm_kvn_string(sgp4s);
    fs::write(omm_kvn_file_path, omm_kvn_string).map_err(|err| GpError::Io(err.to_string()))?;
    Ok(())
}

/// Collect leaf OMM keyword values from one XML omm element
///
/// Wrapper tags (metadata, meanElements, tleParameters) and COMMENT elements
/// are skipped. Keys are stored in uppercase so lookup is case-insensitive.
///
/// # Arguments
/// * `omm` - The omm element to flatten
///
/// # Returns
/// * A map of OMM keyword to text value
#[cfg(feature = "xml")]
fn xml_omm_fields(omm: Node) -> HashMap<String, String> {
    let mut fields = HashMap::new();

    for node in omm.descendants() {
        if !node.is_element() {
            continue;
        }

        let name = node.tag_name().name();

        // Skip the record root, comments, and wrapper elements with children
        if name.eq_ignore_ascii_case("omm") || name.eq_ignore_ascii_case("COMMENT") {
            continue;
        }
        if node.children().any(|child| child.is_element()) {
            continue;
        }

        let text = node.text().unwrap_or("").trim();
        fields.insert(name.to_ascii_uppercase(), text.to_string());
    }

    fields
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS XML format.
///
/// This function parses a string containing one or more OMM records. Each
/// omm element is flattened to keyword and value pairs by local tag name, then
/// initialized with the same GP mapping as [`from_omm_kvn_lines`]. COMMENT
/// elements and wrapper tags are ignored. Empty tags use the same defaults
/// as missing KVN fields.
///
/// # Arguments
/// * `omm_xml_string` - A string containing one or more OMM records in XML format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 parameters
/// * `Err(GpError)` - If the XML document or a field cannot be parsed
///
/// # Errors
/// * [`GpError::InvalidOmmXml`] - If the XML document cannot be parsed
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_xml_string;
///
/// // Define an OMM XML document with one omm record
/// let omm = r#"<ndm>
/// <omm id="CCSDS_OMM_VERS" version="2.0">
/// <body><segment>
/// <metadata>
/// <OBJECT_NAME>2026-106A</OBJECT_NAME>
/// <OBJECT_ID>2026-106A</OBJECT_ID>
/// </metadata>
/// <data><meanElements>
/// <EPOCH>2026-06-14T15:07:48.259488</EPOCH>
/// <MEAN_MOTION>15.11169557</MEAN_MOTION>
/// <ECCENTRICITY>.00147468</ECCENTRICITY>
/// <INCLINATION>97.5103</INCLINATION>
/// <RA_OF_ASC_NODE>247.7605</RA_OF_ASC_NODE>
/// <ARG_OF_PERICENTER>169.6213</ARG_OF_PERICENTER>
/// <MEAN_ANOMALY>190.5325</MEAN_ANOMALY>
/// </meanElements><tleParameters>
/// <EPHEMERIS_TYPE>0</EPHEMERIS_TYPE>
/// <CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE>
/// <NORAD_CAT_ID>69097</NORAD_CAT_ID>
/// <ELEMENT_SET_NO>999</ELEMENT_SET_NO>
/// <REV_AT_EPOCH>459</REV_AT_EPOCH>
/// <BSTAR>.39221734E-3</BSTAR>
/// <MEAN_MOTION_DOT>.6535E-4</MEAN_MOTION_DOT>
/// <MEAN_MOTION_DDOT>0</MEAN_MOTION_DDOT>
/// </tleParameters></data>
/// </segment></body></omm></ndm>"#;
///
/// // Parse the OMM XML into SGP4 propagators
/// let sgp4s = from_omm_xml_string(omm).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "xml")]
pub fn from_omm_xml_string(omm_xml_string: &str) -> Result<Vec<Sgp4>, GpError> {
    // Parse the XML document
    let doc = Document::parse(omm_xml_string).map_err(|_| GpError::InvalidOmmXml)?;

    // Each omm element is one GP record
    let mut sgp4s = Vec::new();
    for node in doc.descendants() {
        if !node.is_element() {
            continue;
        }
        if !node.tag_name().name().eq_ignore_ascii_case("omm") {
            continue;
        }

        let fields = xml_omm_fields(node);
        let sgp4 = sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned())?;
        sgp4s.push(sgp4);
    }

    // Return vector of SGP4 structs
    Ok(sgp4s)
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CCSDS XML format.
///
/// This function parses a file containing one or more OMM records in XML
/// format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_xml_file_path` - A path to a file containing one or more OMM records in XML format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 structs
/// * `Err(GpError)` - If the file cannot be read or the XML cannot be parsed
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be read
/// * [`GpError::InvalidOmmXml`] - If the XML document cannot be parsed
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_xml_file;
///
/// // Define an OMM XML document with one omm record
/// let omm = r#"<ndm>
/// <omm id="CCSDS_OMM_VERS" version="2.0">
/// <body><segment>
/// <metadata>
/// <OBJECT_NAME>2026-106A</OBJECT_NAME>
/// <OBJECT_ID>2026-106A</OBJECT_ID>
/// </metadata>
/// <data><meanElements>
/// <EPOCH>2026-06-14T15:07:48.259488</EPOCH>
/// <MEAN_MOTION>15.11169557</MEAN_MOTION>
/// <ECCENTRICITY>.00147468</ECCENTRICITY>
/// <INCLINATION>97.5103</INCLINATION>
/// <RA_OF_ASC_NODE>247.7605</RA_OF_ASC_NODE>
/// <ARG_OF_PERICENTER>169.6213</ARG_OF_PERICENTER>
/// <MEAN_ANOMALY>190.5325</MEAN_ANOMALY>
/// </meanElements><tleParameters>
/// <NORAD_CAT_ID>69097</NORAD_CAT_ID>
/// <BSTAR>.39221734E-3</BSTAR>
/// <MEAN_MOTION_DOT>.6535E-4</MEAN_MOTION_DOT>
/// <MEAN_MOTION_DDOT>0</MEAN_MOTION_DDOT>
/// </tleParameters></data>
/// </segment></body></omm></ndm>"#;
///
/// // Write the OMM to a temporary file
/// let path = std::env::temp_dir().join("mako_sgp4_from_omm_xml_file.xml");
/// std::fs::write(&path, omm).unwrap();
///
/// // Parse the OMM file into SGP4 propagators
/// let sgp4s = from_omm_xml_file(path.to_str().unwrap()).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "xml")]
pub fn from_omm_xml_file(omm_xml_file_path: &str) -> Result<Vec<Sgp4>, GpError> {
    // Open the OMM XML file
    let omm_xml_string =
        fs::read_to_string(omm_xml_file_path).map_err(|err| GpError::Io(err.to_string()))?;

    // Parse OMM string into a vector of SGP4 structs
    from_omm_xml_string(&omm_xml_string)
}

/// Escape text for inclusion in an XML element
///
/// Replaces the five XML markup characters so satellite names and other
/// free-text fields remain well-formed.
///
/// # Arguments
/// * `value` - The raw text to escape
///
/// # Returns
/// * The escaped text
#[cfg(feature = "xml")]
fn escape_xml_text(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&apos;"),
            _ => escaped.push(c),
        }
    }
    escaped
}

/// Format one XML element from a keyword and value
///
/// Empty values are written as a self-closing tag, matching Celestrak OMM XML.
///
/// # Arguments
/// * `name` - The XML element name
/// * `value` - The element text
///
/// # Returns
/// * One XML element
#[cfg(feature = "xml")]
fn format_xml_tag(name: &str, value: &str) -> String {
    if value.is_empty() {
        format!("<{} />", name)
    } else {
        format!("<{}>{}</{}>", name, escape_xml_text(value), name)
    }
}

/// Build one OMM XML record from a general perturbation element set
///
/// Metadata that is not stored on [`GenPerturbElementSet`] is filled with
/// SGP4 defaults. MEAN_MOTION_DOT and MEAN_MOTION_DDOT are written as the
/// TLE-printed values (stored derivatives divided by 2 and 6). Layout matches
/// Celestrak: the `omm` opening tag on its own line, then the record body
/// packed onto the next line.
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * One `omm` element, including a trailing newline
#[cfg(feature = "xml")]
fn gp_to_omm_xml(gp: &GenPerturbElementSet) -> String {
    // OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
    let mean_motion_dot = gp.first_derivative_of_mean_motion / 2.0;
    let mean_motion_ddot = gp.second_derivative_of_mean_motion / 6.0;

    let mut xml = String::from("<omm id=\"CCSDS_OMM_VERS\" version=\"2.0\">\n");
    xml.push_str("<header><CREATION_DATE /><ORIGINATOR /></header>");
    xml.push_str("<body><segment><metadata>");
    xml.push_str(&format_xml_tag("OBJECT_NAME", &gp.common_name));
    xml.push_str(&format_xml_tag("OBJECT_ID", &gp.international_designator));
    xml.push_str(&format_xml_tag("CENTER_NAME", "EARTH"));
    xml.push_str(&format_xml_tag("REF_FRAME", "TEME"));
    xml.push_str(&format_xml_tag("TIME_SYSTEM", "UTC"));
    xml.push_str(&format_xml_tag("MEAN_ELEMENT_THEORY", "SGP4"));
    xml.push_str("</metadata><data><meanElements>");
    xml.push_str(&format_xml_tag(
        "EPOCH",
        &format_omm_epoch(&gp.epoch_datetime),
    ));
    xml.push_str(&format_xml_tag(
        "MEAN_MOTION",
        &format_omm_decimal(gp.mean_motion),
    ));
    xml.push_str(&format_xml_tag(
        "ECCENTRICITY",
        &format_omm_decimal(gp.eccentricity),
    ));
    xml.push_str(&format_xml_tag(
        "INCLINATION",
        &format_omm_decimal(gp.inclination),
    ));
    xml.push_str(&format_xml_tag(
        "RA_OF_ASC_NODE",
        &format_omm_decimal(gp.right_ascension_of_ascending_node),
    ));
    xml.push_str(&format_xml_tag(
        "ARG_OF_PERICENTER",
        &format_omm_decimal(gp.argument_of_perigee),
    ));
    xml.push_str(&format_xml_tag(
        "MEAN_ANOMALY",
        &format_omm_decimal(gp.mean_anomaly),
    ));
    xml.push_str("</meanElements><tleParameters>");
    xml.push_str(&format_xml_tag(
        "EPHEMERIS_TYPE",
        &gp.ephemeris_type.to_string(),
    ));
    xml.push_str(&format_xml_tag(
        "CLASSIFICATION_TYPE",
        &format_omm_classification(gp.classification),
    ));
    xml.push_str(&format_xml_tag(
        "NORAD_CAT_ID",
        &gp.satellite_catalog_number.to_string(),
    ));
    xml.push_str(&format_xml_tag(
        "ELEMENT_SET_NO",
        &gp.element_set_number.to_string(),
    ));
    xml.push_str(&format_xml_tag(
        "REV_AT_EPOCH",
        &gp.revolution_number_at_epoch.to_string(),
    ));
    xml.push_str(&format_xml_tag("BSTAR", &format_omm_sci(gp.bstar)));
    xml.push_str(&format_xml_tag(
        "MEAN_MOTION_DOT",
        &format_omm_sci(mean_motion_dot),
    ));
    xml.push_str(&format_xml_tag(
        "MEAN_MOTION_DDOT",
        &format_omm_sci(mean_motion_ddot),
    ));
    xml.push_str("</tleParameters></data></segment></body></omm>\n");

    xml
}

/// Builds a CCSDS XML OMM string from a slice of [`Sgp4`] structs.
///
/// Each element set is written as one `omm` record inside an `ndm` document.
/// The `ndm` root includes the CCSDS NDM XML schema location used by
/// Celestrak. Metadata fields that are not stored on
/// [`GenPerturbElementSet`] are filled with SGP4 defaults: Earth-centered,
/// TEME, UTC, and SGP4. CREATION_DATE and ORIGINATOR are left empty.
/// OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
/// (the stored derivatives divided by 2 and 6).
///
/// An empty slice produces an `ndm` document with no `omm` records.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
///
/// # Returns
/// * `String` - An NDM XML document containing one or more OMM records
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_xml_string, to_omm_xml_string};
///
/// // Define an OMM XML document with one omm record
/// let omm = r#"<ndm>
/// <omm id="CCSDS_OMM_VERS" version="2.0">
/// <body><segment>
/// <metadata>
/// <OBJECT_NAME>2026-106A</OBJECT_NAME>
/// <OBJECT_ID>2026-106A</OBJECT_ID>
/// </metadata>
/// <data><meanElements>
/// <EPOCH>2026-06-14T15:07:48.259488</EPOCH>
/// <MEAN_MOTION>15.11169557</MEAN_MOTION>
/// <ECCENTRICITY>.00147468</ECCENTRICITY>
/// <INCLINATION>97.5103</INCLINATION>
/// <RA_OF_ASC_NODE>247.7605</RA_OF_ASC_NODE>
/// <ARG_OF_PERICENTER>169.6213</ARG_OF_PERICENTER>
/// <MEAN_ANOMALY>190.5325</MEAN_ANOMALY>
/// </meanElements><tleParameters>
/// <NORAD_CAT_ID>69097</NORAD_CAT_ID>
/// <BSTAR>.39221734E-3</BSTAR>
/// <MEAN_MOTION_DOT>.6535E-4</MEAN_MOTION_DOT>
/// <MEAN_MOTION_DDOT>0</MEAN_MOTION_DDOT>
/// </tleParameters></data>
/// </segment></body></omm></ndm>"#;
///
/// // Parse, export, and parse again
/// let sgp4s = from_omm_xml_string(omm).unwrap();
/// let exported = to_omm_xml_string(&sgp4s);
/// let reparsed = from_omm_xml_string(&exported).unwrap();
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "xml")]
pub fn to_omm_xml_string(sgp4s: &[Sgp4]) -> String {
    // Wrap all records in an NDM document with the CCSDS schema location
    let mut xml = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <ndm xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
         xsi:noNamespaceSchemaLocation=\"https://sanaregistry.org/r/ndmxml_unqualified/ndmxml-2.0.0-master-2.0.xsd\">\n",
    );
    for sgp4 in sgp4s {
        xml.push_str(&gp_to_omm_xml(&sgp4.gp));
    }
    xml.push_str("</ndm>\n");

    xml
}

/// Writes a CCSDS XML OMM file from a slice of [`Sgp4`] structs.
///
/// This function serializes the provided element sets with
/// [`to_omm_xml_string`] and writes the result to `omm_xml_file_path`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
/// * `omm_xml_file_path` - Destination path for the XML file
///
/// # Returns
/// * `Ok(())` - If the file is written
/// * `Err(GpError)` - If the file cannot be written
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be written
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_xml_file, from_omm_xml_string, to_omm_xml_file};
///
/// // Define an OMM XML document and parse it
/// let omm = r#"<ndm>
/// <omm id="CCSDS_OMM_VERS" version="2.0">
/// <body><segment>
/// <metadata>
/// <OBJECT_NAME>2026-106A</OBJECT_NAME>
/// <OBJECT_ID>2026-106A</OBJECT_ID>
/// </metadata>
/// <data><meanElements>
/// <EPOCH>2026-06-14T15:07:48.259488</EPOCH>
/// <MEAN_MOTION>15.11169557</MEAN_MOTION>
/// <ECCENTRICITY>.00147468</ECCENTRICITY>
/// <INCLINATION>97.5103</INCLINATION>
/// <RA_OF_ASC_NODE>247.7605</RA_OF_ASC_NODE>
/// <ARG_OF_PERICENTER>169.6213</ARG_OF_PERICENTER>
/// <MEAN_ANOMALY>190.5325</MEAN_ANOMALY>
/// </meanElements><tleParameters>
/// <NORAD_CAT_ID>69097</NORAD_CAT_ID>
/// <BSTAR>.39221734E-3</BSTAR>
/// <MEAN_MOTION_DOT>.6535E-4</MEAN_MOTION_DOT>
/// <MEAN_MOTION_DDOT>0</MEAN_MOTION_DDOT>
/// </tleParameters></data>
/// </segment></body></omm></ndm>"#;
/// let sgp4s = from_omm_xml_string(omm).unwrap();
///
/// // Write the element sets to a temporary XML file
/// let path = std::env::temp_dir().join("mako_sgp4_to_omm_xml_file.xml");
/// let path_str = path.to_str().unwrap();
/// to_omm_xml_file(&sgp4s, path_str).unwrap();
///
/// // Parse the written file
/// let reparsed = from_omm_xml_file(path_str).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "xml")]
pub fn to_omm_xml_file(sgp4s: &[Sgp4], omm_xml_file_path: &str) -> Result<(), GpError> {
    // Serialize the element sets and write the XML file
    let omm_xml_string = to_omm_xml_string(sgp4s);
    fs::write(omm_xml_file_path, omm_xml_string).map_err(|err| GpError::Io(err.to_string()))?;
    Ok(())
}

/// Build an [`Sgp4`] struct from one JSON OMM object
///
/// # Arguments
/// * `record` - One OMM record as a JSON object
///
/// # Returns
/// * `Ok(Sgp4)` - Struct containing the parsed SGP4 parameters
/// * `Err(GpError)` - If the JSON value is not an object or a field cannot be parsed
#[cfg(feature = "json")]
fn sgp4_from_json_record(record: &Value) -> Result<Sgp4, GpError> {
    let fields = json_omm_fields(record)?;
    sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned())
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in JSON format.
///
/// This function parses a string containing one OMM object or an array of OMM
/// objects. Each object is flattened to keyword and value pairs, then
/// initialized with the same GP mapping as [`from_omm_kvn_lines`]. JSON numbers
/// are converted to text. Null, missing, and empty string fields use the same
/// defaults as missing KVN fields.
///
/// # Arguments
/// * `omm_json_string` - A string containing one or more OMM records in JSON format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 parameters
/// * `Err(GpError)` - If the JSON document or a field cannot be parsed
///
/// # Errors
/// * [`GpError::InvalidOmmJson`] - If the JSON cannot be parsed or is not an object or array of objects
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_json_string;
///
/// // Define OMM JSON as an array of objects (a single object is also accepted)
/// let omm = r#"[
/// {
/// "OBJECT_NAME": "2026-106A",
/// "OBJECT_ID": "2026-106A",
/// "EPOCH": "2026-06-14T15:07:48.259488",
/// "MEAN_MOTION": 15.11169557,
/// "ECCENTRICITY": 0.00147468,
/// "INCLINATION": 97.5103,
/// "RA_OF_ASC_NODE": 247.7605,
/// "ARG_OF_PERICENTER": 169.6213,
/// "MEAN_ANOMALY": 190.5325,
/// "EPHEMERIS_TYPE": 0,
/// "CLASSIFICATION_TYPE": "U",
/// "NORAD_CAT_ID": 69097,
/// "ELEMENT_SET_NO": 999,
/// "REV_AT_EPOCH": 459,
/// "BSTAR": 0.00039221734,
/// "MEAN_MOTION_DOT": 6.535e-05,
/// "MEAN_MOTION_DDOT": 0
/// }
/// ]"#;
///
/// // Parse the OMM JSON into SGP4 propagators
/// let sgp4s = from_omm_json_string(omm).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "json")]
pub fn from_omm_json_string(omm_json_string: &str) -> Result<Vec<Sgp4>, GpError> {
    // Parse the JSON document
    let value: Value =
        serde_json::from_str(omm_json_string).map_err(|_| GpError::InvalidOmmJson)?;

    // A file may be one object or an array of objects
    let mut sgp4s = Vec::new();
    match &value {
        Value::Array(records) => {
            for record in records {
                sgp4s.push(sgp4_from_json_record(record)?);
            }
        }
        Value::Object(_) => {
            sgp4s.push(sgp4_from_json_record(&value)?);
        }
        _ => return Err(GpError::InvalidOmmJson),
    }

    // Return vector of SGP4 structs
    Ok(sgp4s)
}

/// Collect OMM keyword values from one JSON object
///
/// Keys are stored in uppercase so lookup is case-insensitive. Strings are
/// used as-is, numbers are converted to text, and null becomes an empty
/// string. Nested arrays and objects are ignored.
///
/// # Arguments
/// * `record` - The JSON object to flatten
///
/// # Returns
/// * `Ok(HashMap)` - A map of OMM keyword to text value
/// * `Err(GpError)` - If the JSON value is not an object
#[cfg(feature = "json")]
fn json_omm_fields(record: &Value) -> Result<HashMap<String, String>, GpError> {
    let Some(object) = record.as_object() else {
        return Err(GpError::InvalidOmmJson);
    };

    let mut fields = HashMap::new();
    for (key, value) in object {
        if key.eq_ignore_ascii_case("COMMENT") {
            continue;
        }

        // Convert JSON scalars to the same text form used by KVN and XML
        let text = match value {
            Value::Null => String::new(),
            Value::String(s) => s.trim().to_string(),
            Value::Number(n) => n.to_string(),
            Value::Bool(b) => b.to_string(),
            Value::Array(_) | Value::Object(_) => continue,
        };
        fields.insert(key.to_ascii_uppercase(), text);
    }

    Ok(fields)
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in JSON format.
///
/// This function parses a file containing one OMM object or an array of OMM
/// objects in JSON format and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_json_file_path` - A path to a file containing one or more OMM records in JSON format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 structs
/// * `Err(GpError)` - If the file cannot be read or the JSON cannot be parsed
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be read
/// * [`GpError::InvalidOmmJson`] - If the JSON cannot be parsed or is not an object or array of objects
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_json_file;
///
/// // Define OMM JSON as an array of objects
/// let omm = r#"[
/// {
/// "OBJECT_NAME": "2026-106A",
/// "OBJECT_ID": "2026-106A",
/// "EPOCH": "2026-06-14T15:07:48.259488",
/// "MEAN_MOTION": 15.11169557,
/// "ECCENTRICITY": 0.00147468,
/// "INCLINATION": 97.5103,
/// "RA_OF_ASC_NODE": 247.7605,
/// "ARG_OF_PERICENTER": 169.6213,
/// "MEAN_ANOMALY": 190.5325,
/// "NORAD_CAT_ID": 69097,
/// "BSTAR": 0.00039221734,
/// "MEAN_MOTION_DOT": 6.535e-05,
/// "MEAN_MOTION_DDOT": 0
/// }
/// ]"#;
///
/// // Write the OMM to a temporary file
/// let path = std::env::temp_dir().join("mako_sgp4_from_omm_json_file.json");
/// std::fs::write(&path, omm).unwrap();
///
/// // Parse the OMM file into SGP4 propagators
/// let sgp4s = from_omm_json_file(path.to_str().unwrap()).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "json")]
pub fn from_omm_json_file(omm_json_file_path: &str) -> Result<Vec<Sgp4>, GpError> {
    // Open the OMM JSON file
    let omm_json_string =
        fs::read_to_string(omm_json_file_path).map_err(|err| GpError::Io(err.to_string()))?;

    // Parse OMM string into a vector of SGP4 structs
    from_omm_json_string(&omm_json_string)
}

/// Format an OMM number as a JSON numeric literal
///
/// JSON requires a digit before the decimal point, so this does not use
/// Celestrak KVN scientific notation or a leading-dot fraction. Zero is
/// written as `0`.
///
/// # Arguments
/// * `value` - The number to format
///
/// # Returns
/// * A JSON number literal
#[cfg(feature = "json")]
fn format_omm_json_number(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }

    let rounded = (value * 1e12).round() / 1e12;
    match serde_json::Number::from_f64(rounded) {
        Some(n) => n.to_string(),
        None => "0".to_string(),
    }
}

/// Build one OMM JSON object from a general perturbation element set
///
/// Field names and types match Celestrak GP JSON: strings for names, epoch,
/// and classification; numbers for the remaining GP fields. MEAN_MOTION_DOT
/// and MEAN_MOTION_DDOT are the TLE-printed values (stored derivatives
/// divided by 2 and 6).
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * One pretty-printed JSON object
#[cfg(feature = "json")]
fn gp_to_omm_json(gp: &GenPerturbElementSet) -> String {
    // OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
    let mean_motion_dot = gp.first_derivative_of_mean_motion / 2.0;
    let mean_motion_ddot = gp.second_derivative_of_mean_motion / 6.0;

    let fields = [
        format!(
            "        \"OBJECT_NAME\": {}",
            serde_json::to_string(&gp.common_name).unwrap()
        ),
        format!(
            "        \"OBJECT_ID\": {}",
            serde_json::to_string(&gp.international_designator).unwrap()
        ),
        format!(
            "        \"EPOCH\": {}",
            serde_json::to_string(&format_omm_epoch(&gp.epoch_datetime)).unwrap()
        ),
        format!(
            "        \"MEAN_MOTION\": {}",
            format_omm_json_number(gp.mean_motion)
        ),
        format!(
            "        \"ECCENTRICITY\": {}",
            format_omm_json_number(gp.eccentricity)
        ),
        format!(
            "        \"INCLINATION\": {}",
            format_omm_json_number(gp.inclination)
        ),
        format!(
            "        \"RA_OF_ASC_NODE\": {}",
            format_omm_json_number(gp.right_ascension_of_ascending_node)
        ),
        format!(
            "        \"ARG_OF_PERICENTER\": {}",
            format_omm_json_number(gp.argument_of_perigee)
        ),
        format!(
            "        \"MEAN_ANOMALY\": {}",
            format_omm_json_number(gp.mean_anomaly)
        ),
        format!("        \"EPHEMERIS_TYPE\": {}", gp.ephemeris_type),
        format!(
            "        \"CLASSIFICATION_TYPE\": {}",
            serde_json::to_string(&format_omm_classification(gp.classification)).unwrap()
        ),
        format!("        \"NORAD_CAT_ID\": {}", gp.satellite_catalog_number),
        format!("        \"ELEMENT_SET_NO\": {}", gp.element_set_number),
        format!(
            "        \"REV_AT_EPOCH\": {}",
            gp.revolution_number_at_epoch
        ),
        format!("        \"BSTAR\": {}", format_omm_json_number(gp.bstar)),
        format!(
            "        \"MEAN_MOTION_DOT\": {}",
            format_omm_json_number(mean_motion_dot)
        ),
        format!(
            "        \"MEAN_MOTION_DDOT\": {}",
            format_omm_json_number(mean_motion_ddot)
        ),
    ];

    format!("    {{\n{}\n    }}", fields.join(",\n"))
}

/// Builds a Celestrak-style JSON OMM string from a slice of [`Sgp4`] structs.
///
/// Each element set is written as one JSON object in a top-level array.
/// Strings are used for OBJECT_NAME, OBJECT_ID, EPOCH, and
/// CLASSIFICATION_TYPE. Numeric GP fields are JSON numbers.
/// OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
/// (the stored derivatives divided by 2 and 6).
///
/// An empty slice produces `[]`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
///
/// # Returns
/// * `String` - A JSON array of OMM objects
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_json_string, to_omm_json_string};
///
/// // Define OMM JSON as an array of objects
/// let omm = r#"[
/// {
/// "OBJECT_NAME": "2026-106A",
/// "OBJECT_ID": "2026-106A",
/// "EPOCH": "2026-06-14T15:07:48.259488",
/// "MEAN_MOTION": 15.11169557,
/// "ECCENTRICITY": 0.00147468,
/// "INCLINATION": 97.5103,
/// "RA_OF_ASC_NODE": 247.7605,
/// "ARG_OF_PERICENTER": 169.6213,
/// "MEAN_ANOMALY": 190.5325,
/// "NORAD_CAT_ID": 69097,
/// "BSTAR": 0.00039221734,
/// "MEAN_MOTION_DOT": 6.535e-05,
/// "MEAN_MOTION_DDOT": 0
/// }
/// ]"#;
///
/// // Parse, export, and parse again
/// let sgp4s = from_omm_json_string(omm).unwrap();
/// let exported = to_omm_json_string(&sgp4s);
/// let reparsed = from_omm_json_string(&exported).unwrap();
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "json")]
pub fn to_omm_json_string(sgp4s: &[Sgp4]) -> String {
    // An empty slice is an empty JSON array
    if sgp4s.is_empty() {
        return "[]\n".to_string();
    }

    // Serialize each element set as one JSON object
    let mut objects = Vec::new();
    for sgp4 in sgp4s {
        objects.push(gp_to_omm_json(&sgp4.gp));
    }

    format!("[\n{}\n]\n", objects.join(",\n"))
}

/// Writes a Celestrak-style JSON OMM file from a slice of [`Sgp4`] structs.
///
/// This function serializes the provided element sets with
/// [`to_omm_json_string`] and writes the result to `omm_json_file_path`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
/// * `omm_json_file_path` - Destination path for the JSON file
///
/// # Returns
/// * `Ok(())` - If the file is written
/// * `Err(GpError)` - If the file cannot be written
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be written
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_json_file, from_omm_json_string, to_omm_json_file};
///
/// // Define OMM JSON and parse it
/// let omm = r#"[
/// {
/// "OBJECT_NAME": "2026-106A",
/// "OBJECT_ID": "2026-106A",
/// "EPOCH": "2026-06-14T15:07:48.259488",
/// "MEAN_MOTION": 15.11169557,
/// "ECCENTRICITY": 0.00147468,
/// "INCLINATION": 97.5103,
/// "RA_OF_ASC_NODE": 247.7605,
/// "ARG_OF_PERICENTER": 169.6213,
/// "MEAN_ANOMALY": 190.5325,
/// "NORAD_CAT_ID": 69097,
/// "BSTAR": 0.00039221734,
/// "MEAN_MOTION_DOT": 6.535e-05,
/// "MEAN_MOTION_DDOT": 0
/// }
/// ]"#;
/// let sgp4s = from_omm_json_string(omm).unwrap();
///
/// // Write the element sets to a temporary JSON file
/// let path = std::env::temp_dir().join("mako_sgp4_to_omm_json_file.json");
/// let path_str = path.to_str().unwrap();
/// to_omm_json_file(&sgp4s, path_str).unwrap();
///
/// // Parse the written file
/// let reparsed = from_omm_json_file(path_str).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "json")]
pub fn to_omm_json_file(sgp4s: &[Sgp4], omm_json_file_path: &str) -> Result<(), GpError> {
    // Serialize the element sets and write the JSON file
    let omm_json_string = to_omm_json_string(sgp4s);
    fs::write(omm_json_file_path, omm_json_string).map_err(|err| GpError::Io(err.to_string()))?;
    Ok(())
}

/// Build an [`Sgp4`] struct from one CSV OMM row
///
/// # Arguments
/// * `headers` - Uppercase OMM keywords from the CSV header row
/// * `record` - One OMM data row
///
/// # Returns
/// * `Ok(Sgp4)` - Struct containing the parsed SGP4 parameters
/// * `Err(GpError)` - If a field cannot be parsed
#[cfg(feature = "csv")]
fn sgp4_from_csv_record(headers: &[String], record: &csv::StringRecord) -> Result<Sgp4, GpError> {
    let fields = csv_omm_fields(headers, record);
    sgp4_from_omm_lookup(|field| fields.get(&field.to_ascii_uppercase()).cloned())
}

/// Builds a vector of [`Sgp4`] structs from a string containing Orbit Mean-Elements Message
/// (OMM) sets in CSV format.
///
/// This function parses a string containing a header row of OMM keywords and
/// one data row per record. Each row is flattened to keyword and value pairs,
/// then initialized with the same GP mapping as [`from_omm_kvn_lines`]. Empty
/// cells use the same defaults as missing KVN fields.
///
/// # Arguments
/// * `omm_csv_string` - A string containing one or more OMM records in CSV format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 parameters
/// * `Err(GpError)` - If the CSV document or a field cannot be parsed
///
/// # Errors
/// * [`GpError::InvalidOmmCsv`] - If the CSV document cannot be parsed
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_csv_string;
///
/// // Define an OMM CSV with a header row and one data row
/// let omm = "\
/// OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
/// 2026-106A,2026-106A,2026-06-14T15:07:48.259488,15.11169557,.00147468,97.5103,247.7605,169.6213,190.5325,0,U,69097,999,459,.39221734E-3,.6535E-4,0
/// ";
///
/// // Parse the OMM CSV into SGP4 propagators
/// let sgp4s = from_omm_csv_string(omm).unwrap();
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "csv")]
pub fn from_omm_csv_string(omm_csv_string: &str) -> Result<Vec<Sgp4>, GpError> {
    // Parse the CSV document, using the first row as OMM keywords
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(omm_csv_string.as_bytes());

    let headers: Vec<String> = reader
        .headers()
        .map_err(|_| GpError::InvalidOmmCsv)?
        .iter()
        .map(|header| header.trim().to_ascii_uppercase())
        .collect();

    // Each data row is one GP record
    let mut sgp4s = Vec::new();
    for result in reader.records() {
        let record = result.map_err(|_| GpError::InvalidOmmCsv)?;
        sgp4s.push(sgp4_from_csv_record(&headers, &record)?);
    }

    // Return vector of SGP4 structs
    Ok(sgp4s)
}

/// Collect OMM keyword values from one CSV row
///
/// Keys are stored in uppercase so lookup is case-insensitive. COMMENT columns
/// are skipped. Empty cells become empty strings so missing-field defaults apply.
///
/// # Arguments
/// * `headers` - Uppercase OMM keywords from the CSV header row
/// * `record` - One OMM data row
///
/// # Returns
/// * A map of OMM keyword to text value
#[cfg(feature = "csv")]
fn csv_omm_fields(headers: &[String], record: &csv::StringRecord) -> HashMap<String, String> {
    let mut fields = HashMap::new();

    for (index, header) in headers.iter().enumerate() {
        if header.is_empty() || header.eq_ignore_ascii_case("COMMENT") {
            continue;
        }

        let text = record.get(index).unwrap_or("").trim().to_string();
        fields.insert(header.clone(), text);
    }

    fields
}

/// Builds a vector of [`Sgp4`] structs from a file containing Orbit Mean-Elements Message
/// (OMM) sets in CSV format.
///
/// This function parses a file containing a header row of OMM keywords and
/// one data row per record and returns all successfully parsed entries.
///
/// # Arguments
/// * `omm_csv_file_path` - A path to a file containing one or more OMM records in CSV format
///
/// # Returns
/// * `Ok(Vec<Sgp4>)` - All successfully parsed SGP4 structs
/// * `Err(GpError)` - If the file cannot be read or the CSV cannot be parsed
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be read
/// * [`GpError::InvalidOmmCsv`] - If the CSV document cannot be parsed
/// * [`GpError::InvalidOmmField`] - If a present field cannot be parsed as the requested type
/// * [`GpError::InvalidOmmEpoch`] - If EPOCH cannot be parsed as a UTC datetime
/// * [`GpError::Sgp4`] - If SGP4 initialization fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::from_omm_csv_file;
///
/// // Define an OMM CSV with a header row and one data row
/// let omm = "\
/// OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,NORAD_CAT_ID,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
/// 2026-106A,2026-106A,2026-06-14T15:07:48.259488,15.11169557,.00147468,97.5103,247.7605,169.6213,190.5325,69097,.39221734E-3,.6535E-4,0
/// ";
///
/// // Write the OMM to a temporary file
/// let path = std::env::temp_dir().join("mako_sgp4_from_omm_csv_file.csv");
/// std::fs::write(&path, omm).unwrap();
///
/// // Parse the OMM file into SGP4 propagators
/// let sgp4s = from_omm_csv_file(path.to_str().unwrap()).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the catalog number and name
/// assert_eq!(sgp4s.len(), 1);
/// assert_eq!(sgp4s[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(sgp4s[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "csv")]
pub fn from_omm_csv_file(omm_csv_file_path: &str) -> Result<Vec<Sgp4>, GpError> {
    // Open the OMM CSV file
    let omm_csv_string =
        fs::read_to_string(omm_csv_file_path).map_err(|err| GpError::Io(err.to_string()))?;

    // Parse OMM string into a vector of SGP4 structs
    from_omm_csv_string(&omm_csv_string)
}

/// Build one OMM CSV data row from a general perturbation element set
///
/// MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values (stored
/// derivatives divided by 2 and 6). BSTAR and the mean-motion derivatives
/// use Celestrak scientific notation.
///
/// # Arguments
/// * `gp` - The element set to serialize
///
/// # Returns
/// * One CSV data row, in [`OMM_CSV_HEADERS`] order
#[cfg(feature = "csv")]
fn gp_to_omm_csv_row(gp: &GenPerturbElementSet) -> [String; 17] {
    // OMM MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed values
    let mean_motion_dot = gp.first_derivative_of_mean_motion / 2.0;
    let mean_motion_ddot = gp.second_derivative_of_mean_motion / 6.0;

    [
        gp.common_name.clone(),
        gp.international_designator.clone(),
        format_omm_epoch(&gp.epoch_datetime),
        format_omm_decimal(gp.mean_motion),
        format_omm_decimal(gp.eccentricity),
        format_omm_decimal(gp.inclination),
        format_omm_decimal(gp.right_ascension_of_ascending_node),
        format_omm_decimal(gp.argument_of_perigee),
        format_omm_decimal(gp.mean_anomaly),
        gp.ephemeris_type.to_string(),
        format_omm_classification(gp.classification),
        gp.satellite_catalog_number.to_string(),
        gp.element_set_number.to_string(),
        gp.revolution_number_at_epoch.to_string(),
        format_omm_sci(gp.bstar),
        format_omm_sci(mean_motion_dot),
        format_omm_sci(mean_motion_ddot),
    ]
}

/// Builds a Celestrak-style CSV OMM string from a slice of [`Sgp4`] structs.
///
/// The first row is the OMM keyword header. Each following row is one
/// element set. MEAN_MOTION_DOT and MEAN_MOTION_DDOT are the TLE-printed
/// values (the stored derivatives divided by 2 and 6).
///
/// An empty slice produces a header-only CSV.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
///
/// # Returns
/// * `String` - A CSV document with a header row and one row per record
///
/// # Panics
/// * If the in-memory CSV writer fails
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_csv_string, to_omm_csv_string};
///
/// // Define an OMM CSV with a header row and one data row
/// let omm = "\
/// OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
/// 2026-106A,2026-106A,2026-06-14T15:07:48.259488,15.11169557,.00147468,97.5103,247.7605,169.6213,190.5325,0,U,69097,999,459,.39221734E-3,.6535E-4,0
/// ";
///
/// // Parse, export, and parse again
/// let sgp4s = from_omm_csv_string(omm).unwrap();
/// let exported = to_omm_csv_string(&sgp4s);
/// let reparsed = from_omm_csv_string(&exported).unwrap();
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "csv")]
pub fn to_omm_csv_string(sgp4s: &[Sgp4]) -> String {
    // Write the header row, then one data row per element set
    let mut writer = WriterBuilder::new().from_writer(Vec::new());
    writer
        .write_record(OMM_CSV_HEADERS)
        .expect("Cannot write OMM CSV");
    for sgp4 in sgp4s {
        writer
            .write_record(gp_to_omm_csv_row(&sgp4.gp))
            .expect("Cannot write OMM CSV");
    }

    let bytes = writer.into_inner().expect("Cannot write OMM CSV");
    String::from_utf8(bytes).expect("Cannot write OMM CSV")
}

/// Writes a Celestrak-style CSV OMM file from a slice of [`Sgp4`] structs.
///
/// This function serializes the provided element sets with
/// [`to_omm_csv_string`] and writes the result to `omm_csv_file_path`.
///
/// # Arguments
/// * `sgp4s` - The SGP4 structs to export
/// * `omm_csv_file_path` - Destination path for the CSV file
///
/// # Returns
/// * `Ok(())` - If the file is written
/// * `Err(GpError)` - If the file cannot be written
///
/// # Errors
/// * [`GpError::Io`] - If the file cannot be written
///
/// # Examples
/// ```rust
/// use mako_sgp4::gp::{from_omm_csv_file, from_omm_csv_string, to_omm_csv_file};
///
/// // Define an OMM CSV and parse it
/// let omm = "\
/// OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,NORAD_CAT_ID,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
/// 2026-106A,2026-106A,2026-06-14T15:07:48.259488,15.11169557,.00147468,97.5103,247.7605,169.6213,190.5325,69097,.39221734E-3,.6535E-4,0
/// ";
/// let sgp4s = from_omm_csv_string(omm).unwrap();
///
/// // Write the element sets to a temporary CSV file
/// let path = std::env::temp_dir().join("mako_sgp4_to_omm_csv_file.csv");
/// let path_str = path.to_str().unwrap();
/// to_omm_csv_file(&sgp4s, path_str).unwrap();
///
/// // Parse the written file
/// let reparsed = from_omm_csv_file(path_str).unwrap();
///
/// // Remove the temporary file
/// let _ = std::fs::remove_file(&path);
///
/// // Assert the round-trip catalog number and name
/// assert_eq!(reparsed.len(), 1);
/// assert_eq!(reparsed[0].gp.satellite_catalog_number, 69097);
/// assert_eq!(reparsed[0].gp.common_name, "2026-106A");
/// ```
///
/// # References
/// - [CCSDS Orbit Data Messages Specification](https://ccsds.org/Pubs/502x0b3e1.pdf)
/// - [CCSDS XML Specification for Navigation Data Messages](https://ccsds.org/Pubs/505x0b3e2.pdf)
/// - [Celestrak GP Data Formats](https://celestrak.org/NORAD/documentation/gp-data-formats.php)
#[cfg(feature = "csv")]
pub fn to_omm_csv_file(sgp4s: &[Sgp4], omm_csv_file_path: &str) -> Result<(), GpError> {
    // Serialize the element sets and write the CSV file
    let omm_csv_string = to_omm_csv_string(sgp4s);
    fs::write(omm_csv_file_path, omm_csv_string).map_err(|err| GpError::Io(err.to_string()))?;
    Ok(())
}

/// Parse an OMM EPOCH string into a UTC DateTime
///
/// Accepts an ISO-8601 date and time with optional fractional seconds and an
/// optional trailing Z. The date and time may be separated by T or a space.
///
/// # Arguments
/// * `epoch` - The EPOCH string from the OMM record
///
/// # Returns
/// * `Ok(DateTime)` - The epoch as a UTC datetime
/// * `Err(GpError)` - If the EPOCH string cannot be parsed
fn parse_omm_epoch(epoch: &str) -> Result<DateTime, GpError> {
    // Strip a trailing Z timezone marker
    let mut s = epoch.trim();
    if s.ends_with('Z') || s.ends_with('z') {
        s = &s[..s.len() - 1];
    }

    // Split date and time
    let (date, time) = if let Some((date, time)) = s.split_once('T') {
        (date, time)
    } else if let Some((date, time)) = s.split_once(' ') {
        (date, time)
    } else {
        return Err(GpError::InvalidOmmEpoch);
    };

    let date_parts: Vec<&str> = date.split('-').collect();
    let time_parts: Vec<&str> = time.split(':').collect();
    if date_parts.len() != 3 || time_parts.len() != 3 {
        return Err(GpError::InvalidOmmEpoch);
    }

    // Parse calendar date
    let year = date_parts[0]
        .parse::<i32>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;
    let month = date_parts[1]
        .parse::<i32>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;
    let day = date_parts[2]
        .parse::<i32>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;

    // Parse clock time
    let hour = time_parts[0]
        .parse::<i32>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;
    let minute = time_parts[1]
        .parse::<i32>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;
    let second = time_parts[2]
        .parse::<f64>()
        .map_err(|_| GpError::InvalidOmmEpoch)?;

    Ok(DateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
        timezone: Timezone::UTC,
    })
}

// ----------
// Unit Tests
// ----------
#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use std::collections::HashMap;
    use toml::from_str;

    /// Reject a TLE whose two data lines name different catalog numbers
    ///
    /// # Panics
    /// * If parsing returns anything other than a catalog mismatch
    #[test]
    fn test_tle_catalog_mismatch() {
        let line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
        let line2 = "2 25545  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
        let err = match from_tle_lines(line1, line2, None) {
            Err(err) => err,
            Ok(_) => panic!("expected MismatchedTleCatalog"),
        };
        assert_eq!(err, GpError::MismatchedTleCatalog);
    }

    /// Reject a TLE path that cannot be read
    ///
    /// # Panics
    /// * If the missing file returns anything other than an I/O error
    #[test]
    fn test_tle_file_io_error() {
        let err = match from_tle_file("test/this_tle_file_does_not_exist.txt") {
            Err(err) => err,
            Ok(_) => panic!("expected GpError::Io"),
        };
        assert!(matches!(err, GpError::Io(_)));
    }

    /// Reject catalog numbers outside the TLE range
    ///
    /// # Panics
    /// * If a negative number or 340000 is accepted
    #[test]
    fn test_invalid_tle_catalog_number() {
        let err = match format_tle_catalog_number(-1) {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLECatalogNumber"),
        };
        assert_eq!(err, GpError::InvalidTLECatalogNumber);

        let err = match format_tle_catalog_number(340000) {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLECatalogNumber"),
        };
        assert_eq!(err, GpError::InvalidTLECatalogNumber);
    }

    /// Reject TLE epochs outside 1957 through 2056
    ///
    /// # Panics
    /// * If 1956 or 2057 is accepted
    #[test]
    fn test_invalid_tle_datetime() {
        let mut epoch = DateTime {
            year: 1956,
            month: 12,
            day: 31,
            hour: 0,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };
        let err = match format_tle_epoch(&epoch) {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLEDateTime"),
        };
        assert_eq!(err, GpError::InvalidTLEDateTime);

        epoch.year = 2057;
        let err = match format_tle_epoch(&epoch) {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLEDateTime"),
        };
        assert_eq!(err, GpError::InvalidTLEDateTime);
    }

    /// Reject a TLE data line that is too short for a checksum
    ///
    /// # Panics
    /// * If the short line is accepted
    #[test]
    fn test_invalid_tle_line() {
        let err = match tle_line_with_checksum("1 25544U") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLELine"),
        };
        assert_eq!(err, GpError::InvalidTLELine);
    }

    /// Check TLE international designator length limits
    ///
    /// Launch numbers and piece designators that do not fit the TLE columns are
    /// rejected. Designators that fit are formatted into the eight-character field.
    ///
    /// # Panics
    /// * If an overlong designator is accepted or a fitting one is formatted wrong
    #[test]
    fn test_invalid_tle_international_designator() {
        let err = match format_tle_intl_des("2026-1234A") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLEInternationalDesignator"),
        };
        assert_eq!(err, GpError::InvalidTLEInternationalDesignator);

        let err = match format_tle_intl_des("2026-106ABCD") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidTLEInternationalDesignator"),
        };
        assert_eq!(err, GpError::InvalidTLEInternationalDesignator);

        assert_eq!(format_tle_intl_des("1998-067A").unwrap(), "98067A  ");
        assert_eq!(format_tle_intl_des("2025-286AD").unwrap(), "25286AD ");
        assert_eq!(format_tle_intl_des("2025-001ABC").unwrap(), "25001ABC");
    }

    /// Reject OMM fields that are not the expected type
    ///
    /// # Panics
    /// * If a non-numeric mean motion or catalog id is accepted
    #[test]
    fn test_invalid_omm_field() {
        let err = match from_omm_kvn_string("MEAN_MOTION = not_a_number") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmField"),
        };
        assert_eq!(err, GpError::InvalidOmmField);

        let err = match from_omm_kvn_string("NORAD_CAT_ID = abc") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmField"),
        };
        assert_eq!(err, GpError::InvalidOmmField);
    }

    /// Reject OMM epochs that are not a full UTC timestamp
    ///
    /// # Panics
    /// * If a non-date, a date without a time, or a bad clock time is accepted
    #[test]
    fn test_invalid_omm_epoch() {
        let err = match from_omm_kvn_string("EPOCH = not-an-epoch") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmEpoch"),
        };
        assert_eq!(err, GpError::InvalidOmmEpoch);

        let err = match from_omm_kvn_string("EPOCH = 2026-06-14") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmEpoch"),
        };
        assert_eq!(err, GpError::InvalidOmmEpoch);

        let err = match from_omm_kvn_string("EPOCH = 2026-06-14Txx:00:00") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmEpoch"),
        };
        assert_eq!(err, GpError::InvalidOmmEpoch);
    }

    /// Reject OMM XML that is not a document
    ///
    /// # Panics
    /// * If the bad text is accepted
    #[cfg(feature = "xml")]
    #[test]
    fn test_invalid_omm_xml() {
        let err = match from_omm_xml_string("<<<not xml") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmXml"),
        };
        assert_eq!(err, GpError::InvalidOmmXml);
    }

    /// Reject OMM JSON that is not an object or an array of objects
    ///
    /// # Panics
    /// * If invalid JSON, a number, or an array of numbers is accepted
    #[cfg(feature = "json")]
    #[test]
    fn test_invalid_omm_json() {
        let err = match from_omm_json_string("not json") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmJson"),
        };
        assert_eq!(err, GpError::InvalidOmmJson);

        let err = match from_omm_json_string("1") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmJson"),
        };
        assert_eq!(err, GpError::InvalidOmmJson);

        let err = match from_omm_json_string("[1]") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmJson"),
        };
        assert_eq!(err, GpError::InvalidOmmJson);
    }

    /// Reject OMM CSV whose rows do not match the header
    ///
    /// # Panics
    /// * If the short row is accepted
    #[cfg(feature = "csv")]
    #[test]
    fn test_invalid_omm_csv() {
        let err = match from_omm_csv_string("OBJECT_NAME,MEAN_MOTION\nonly_one_field") {
            Err(err) => err,
            Ok(_) => panic!("expected InvalidOmmCsv"),
        };
        assert_eq!(err, GpError::InvalidOmmCsv);
    }

    /// Reject an OMM path that cannot be read
    ///
    /// # Panics
    /// * If the missing file returns anything other than an I/O error
    #[test]
    fn test_omm_file_io_error() {
        let err = match from_omm_kvn_file("test/this_omm_file_does_not_exist.txt") {
            Err(err) => err,
            Ok(_) => panic!("expected GpError::Io"),
        };
        assert!(matches!(err, GpError::Io(_)));
    }

    /// Expand a two-digit TLE year into 1957 through 2056
    ///
    /// # Panics
    /// * If a year maps to the wrong century
    #[test]
    fn test_tle_full_year() {
        assert_eq!(tle_full_year(0), 2000);
        assert_eq!(tle_full_year(25), 2025);
        assert_eq!(tle_full_year(56), 2056);
        assert_eq!(tle_full_year(57), 1957);
        assert_eq!(tle_full_year(98), 1998);
        assert_eq!(tle_full_year(99), 1999);
    }

    /// Convert between Alpha-5 digits and their numeric values
    ///
    /// Letters I and O are skipped. Values outside 10 through 33 have no letter.
    ///
    /// # Panics
    /// * If a letter or value maps to the wrong digit
    #[test]
    fn test_alpha5_digit() {
        assert_eq!(alpha5_digit('A'), Some(10));
        assert_eq!(alpha5_digit('H'), Some(17));
        assert_eq!(alpha5_digit('J'), Some(18));
        assert_eq!(alpha5_digit('N'), Some(22));
        assert_eq!(alpha5_digit('P'), Some(23));
        assert_eq!(alpha5_digit('Z'), Some(33));
        assert_eq!(alpha5_digit('a'), Some(10));
        assert_eq!(alpha5_digit('I'), None);
        assert_eq!(alpha5_digit('O'), None);
        assert_eq!(alpha5_digit('0'), None);

        assert_eq!(alpha5_letter(10), Some('A'));
        assert_eq!(alpha5_letter(17), Some('H'));
        assert_eq!(alpha5_letter(18), Some('J'));
        assert_eq!(alpha5_letter(22), Some('N'));
        assert_eq!(alpha5_letter(23), Some('P'));
        assert_eq!(alpha5_letter(33), Some('Z'));
        assert_eq!(alpha5_letter(9), None);
        assert_eq!(alpha5_letter(34), None);
    }

    /// Compute the checksum of a known TLE line
    ///
    /// # Panics
    /// * If the checksum is not 1
    #[test]
    fn test_checksum_calculation() {
        // Define the TLE line
        let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";

        // Calculate the checksum of the TLE line
        let checksum = calc_checksum(tle_line1);

        // Assert the checksum is correct
        assert_eq!(checksum, 1);
    }

    /// Accept a TLE line with the right checksum and reject a wrong one
    ///
    /// # Panics
    /// * If the valid line fails or the altered line passes
    #[test]
    fn test_checksum_comparison() {
        // Define the TLE line
        let tle_line1 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2921";
        let tle_line2 = "1 25544U 98067A   08264.51782528 -.00002182 -00100-2 -11606-4 0  2922";

        // Calculate the checksum of the TLE line
        let checksum = tle_checksum(tle_line1);
        let checksum2 = tle_checksum(tle_line2);

        // Assert the checksum is correct
        assert!(checksum);
        assert!(!checksum2);
    }

    // -------------------------------------------------------
    // Structs for deserializing the TLE / OMM parsing test TOML
    // -------------------------------------------------------

    #[derive(Deserialize)]
    struct ParsingCases {
        test: HashMap<String, ParsingCase>,
    }

    #[derive(Deserialize)]
    struct ParsingCase {
        name: String,
        #[serde(default)]
        tle: String,
        #[serde(default)]
        omm_kvn: String,
        #[serde(default)]
        exception: bool,
        #[serde(default)]
        common_name: String,
        #[serde(default)]
        satellite_catalog_number: i32,
        #[serde(default)]
        classification: String,
        #[serde(default)]
        international_designator: String,
        #[serde(default)]
        epoch_datetime_year: i32,
        #[serde(default)]
        epoch_datetime_month: i32,
        #[serde(default)]
        epoch_datetime_day: i32,
        #[serde(default)]
        epoch_datetime_hour: i32,
        #[serde(default)]
        epoch_datetime_minute: i32,
        #[serde(default)]
        epoch_datetime_second: f64,
        #[serde(default)]
        first_derivative_of_mean_motion: f64,
        #[serde(default)]
        second_derivative_of_mean_motion: f64,
        #[serde(default)]
        bstar: f64,
        #[serde(default)]
        ephemeris_type: i32,
        #[serde(default)]
        element_set_number: i32,
        #[serde(default)]
        inclination: f64,
        #[serde(default)]
        right_ascension_of_ascending_node: f64,
        #[serde(default)]
        eccentricity: f64,
        #[serde(default)]
        argument_of_perigee: f64,
        #[serde(default)]
        mean_anomaly: f64,
        #[serde(default)]
        mean_motion: f64,
        #[serde(default)]
        revolution_number_at_epoch: i64,
    }

    // Absolute tolerance for parsed GP numeric fields.
    const TLE_PARSE_TOL: f64 = 1e-9;

    /// Compare one parsed numeric field with its expected value
    ///
    /// # Arguments
    /// * `key` - TOML case key, included in the failure message
    /// * `name` - Case description, included in the failure message
    /// * `source` - Parser that produced the value, included in the failure message
    /// * `label` - Field name, included in the failure message
    /// * `value` - Parsed number
    /// * `expected` - Number stored in the test case
    ///
    /// # Panics
    /// * If the absolute difference is at least the parse tolerance
    fn assert_near(key: &str, name: &str, source: &str, label: &str, value: f64, expected: f64) {
        assert!(
            (value - expected).abs() < TLE_PARSE_TOL,
            "{key} ({name}) via {source}: {label} mismatch {value} vs {expected}"
        );
    }

    /// Compare a parsed element set with one TOML case
    ///
    /// # Arguments
    /// * `key` - TOML case key, included in the failure message
    /// * `name` - Case description, included in the failure message
    /// * `source` - Parser that produced the element set
    /// * `gp` - Parsed element set
    /// * `case` - Expected fields
    ///
    /// # Panics
    /// * If any field differs from the case
    fn assert_gp_matches(
        key: &str,
        name: &str,
        source: &str,
        gp: &GenPerturbElementSet,
        case: &ParsingCase,
    ) {
        assert_eq!(
            gp.common_name, case.common_name,
            "{key} ({name}) via {source}: common_name"
        );
        assert_eq!(
            gp.satellite_catalog_number, case.satellite_catalog_number,
            "{key} ({name}) via {source}: satellite_catalog_number"
        );
        assert_eq!(
            gp.classification.to_string(),
            case.classification,
            "{key} ({name}) via {source}: classification"
        );
        assert_eq!(
            gp.international_designator, case.international_designator,
            "{key} ({name}) via {source}: international_designator"
        );
        assert_eq!(
            gp.epoch_datetime.year, case.epoch_datetime_year,
            "{key} ({name}) via {source}: epoch year"
        );
        assert_eq!(
            gp.epoch_datetime.month, case.epoch_datetime_month,
            "{key} ({name}) via {source}: epoch month"
        );
        assert_eq!(
            gp.epoch_datetime.day, case.epoch_datetime_day,
            "{key} ({name}) via {source}: epoch day"
        );
        assert_eq!(
            gp.epoch_datetime.hour, case.epoch_datetime_hour,
            "{key} ({name}) via {source}: epoch hour"
        );
        assert_eq!(
            gp.epoch_datetime.minute, case.epoch_datetime_minute,
            "{key} ({name}) via {source}: epoch minute"
        );
        assert_near(
            key,
            name,
            source,
            "epoch second",
            gp.epoch_datetime.second,
            case.epoch_datetime_second,
        );
        assert_near(
            key,
            name,
            source,
            "n-dot",
            gp.first_derivative_of_mean_motion,
            case.first_derivative_of_mean_motion,
        );
        assert_near(
            key,
            name,
            source,
            "n-ddot",
            gp.second_derivative_of_mean_motion,
            case.second_derivative_of_mean_motion,
        );
        assert_near(key, name, source, "bstar", gp.bstar, case.bstar);
        assert_eq!(
            gp.ephemeris_type, case.ephemeris_type,
            "{key} ({name}) via {source}: ephemeris_type"
        );
        assert_eq!(
            gp.element_set_number, case.element_set_number,
            "{key} ({name}) via {source}: element_set_number"
        );
        assert_near(
            key,
            name,
            source,
            "inclination",
            gp.inclination,
            case.inclination,
        );
        assert_near(
            key,
            name,
            source,
            "raan",
            gp.right_ascension_of_ascending_node,
            case.right_ascension_of_ascending_node,
        );
        assert_near(
            key,
            name,
            source,
            "eccentricity",
            gp.eccentricity,
            case.eccentricity,
        );
        assert_near(
            key,
            name,
            source,
            "argp",
            gp.argument_of_perigee,
            case.argument_of_perigee,
        );
        assert_near(
            key,
            name,
            source,
            "mean anomaly",
            gp.mean_anomaly,
            case.mean_anomaly,
        );
        assert_near(
            key,
            name,
            source,
            "mean motion",
            gp.mean_motion,
            case.mean_motion,
        );
        assert_eq!(
            gp.revolution_number_at_epoch, case.revolution_number_at_epoch,
            "{key} ({name}) via {source}: revolution_number_at_epoch"
        );
    }

    /// Parse one TLE case, with or without a name line
    ///
    /// # Arguments
    /// * `tle` - TLE text from the test case
    ///
    /// # Returns
    /// * `Ok(Sgp4)` - Parsed propagator
    /// * `Err(GpError)` - If the TLE is invalid
    fn sgp4_from_case_tle(tle: &str) -> Result<Sgp4, GpError> {
        let lines: Vec<&str> = tle
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .collect();
        if lines[0].starts_with('1') {
            from_tle_lines(lines[0], lines[1], None)
        } else {
            from_tle_lines(lines[1], lines[2], Some(lines[0]))
        }
    }

    /// Compare TLE parsing with the TOML cases and the fixture file
    ///
    /// Exception cases must fail. The other cases must match from lines, from a
    /// string, and from the fixture file.
    ///
    /// # Panics
    /// * If a case fails to parse or a field misses the expected value
    #[test]
    fn test_tle_parsing_cases() {
        let content = std::fs::read_to_string("test/tle_parsing_cases.toml")
            .expect("could not read test/tle_parsing_cases.toml");
        let cases: ParsingCases =
            from_str(&content).expect("could not parse test/tle_parsing_cases.toml");

        let from_file = from_tle_file("test/tle_parsing_cases.txt")
            .expect("could not read test/tle_parsing_cases.txt");
        let expected_file_count = cases.test.values().filter(|c| !c.exception).count();
        assert_eq!(
            from_file.len(),
            expected_file_count,
            "tle_parsing_cases.txt should contain one entry per non-error TOML case"
        );

        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();

        for key in keys {
            let case = &cases.test[key];

            if case.exception {
                assert!(
                    sgp4_from_case_tle(&case.tle).is_err(),
                    "case {key} ({}): expected from_tle_lines to return Err",
                    case.name
                );
                assert!(
                    from_tle_string(&case.tle).is_err(),
                    "case {key} ({}): expected from_tle_string to return Err",
                    case.name
                );
                continue;
            }

            let from_lines = sgp4_from_case_tle(&case.tle).unwrap_or_else(|err| {
                panic!("case {key} ({}): from_tle_lines failed: {err:?}", case.name)
            });
            assert_gp_matches(key, &case.name, "from_tle_lines", &from_lines.gp, case);

            let from_string = from_tle_string(&case.tle).unwrap_or_else(|err| {
                panic!(
                    "case {key} ({}): from_tle_string failed: {err:?}",
                    case.name
                )
            });
            assert_eq!(
                from_string.len(),
                1,
                "case {key} ({}): expected one TLE from string",
                case.name
            );
            assert_gp_matches(key, &case.name, "from_tle_string", &from_string[0].gp, case);

            let file_match = from_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from tle_parsing_cases.txt parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_tle_file", &file_match.gp, case);
        }
    }

    /// Load the OMM parsing cases from TOML
    ///
    /// # Returns
    /// * `ParsingCases` - Parsed case file
    ///
    /// # Panics
    /// * If the file cannot be read or parsed
    fn load_omm_parsing_cases() -> ParsingCases {
        let content = std::fs::read_to_string("test/omm_parsing_cases.toml")
            .expect("could not read test/omm_parsing_cases.toml");
        from_str(&content).expect("could not parse test/omm_parsing_cases.toml")
    }

    /// Return OMM case keys in a stable order
    ///
    /// # Arguments
    /// * `cases` - Parsed OMM case file
    ///
    /// # Returns
    /// * `Vec<&String>` - Sorted case keys
    fn sorted_omm_keys(cases: &ParsingCases) -> Vec<&String> {
        let mut keys: Vec<&String> = cases.test.keys().collect();
        keys.sort();
        keys
    }

    /// Compare an OMM fixture parse with the non-error TOML cases
    ///
    /// # Arguments
    /// * `cases` - Parsed OMM case file
    /// * `from_file` - Propagators parsed from the fixture
    /// * `source` - Parser name, included in failure messages
    /// * `fixture` - Fixture file name, included in failure messages
    ///
    /// # Panics
    /// * If the record count differs or a catalog number is missing or mismatched
    #[cfg(any(feature = "xml", feature = "json", feature = "csv"))]
    fn assert_omm_file_matches_cases(
        cases: &ParsingCases,
        from_file: &[Sgp4],
        source: &str,
        fixture: &str,
    ) {
        let expected_file_count = cases.test.values().filter(|c| !c.exception).count();
        assert_eq!(
            from_file.len(),
            expected_file_count,
            "{fixture} should contain one entry per non-error TOML case"
        );

        for key in sorted_omm_keys(cases) {
            let case = &cases.test[key];
            if case.exception {
                continue;
            }

            let file_match = from_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from {fixture} parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, source, &file_match.gp, case);
        }
    }

    /// Parse one KVN OMM case into a propagator
    ///
    /// # Arguments
    /// * `omm_kvn` - KVN text from the test case
    ///
    /// # Returns
    /// * `Ok(Sgp4)` - Parsed propagator
    /// * `Err(GpError)` - If the OMM is invalid
    fn sgp4_from_case_omm(omm_kvn: &str) -> Result<Sgp4, GpError> {
        let lines: Vec<&str> = omm_kvn
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .collect();
        from_omm_kvn_lines(&lines)
    }

    /// Compare KVN OMM parsing with the TOML cases and the fixture file
    ///
    /// Exception cases must fail. The other cases must match from lines, from a
    /// string, and from the fixture file.
    ///
    /// # Panics
    /// * If a case fails to parse or a field misses the expected value
    #[test]
    fn test_omm_kvn_parsing_cases() {
        let cases = load_omm_parsing_cases();
        let from_file = from_omm_kvn_file("test/omm_parsing_cases.txt")
            .expect("could not parse test/omm_parsing_cases.txt");
        let expected_file_count = cases.test.values().filter(|c| !c.exception).count();
        assert_eq!(
            from_file.len(),
            expected_file_count,
            "omm_parsing_cases.txt should contain one entry per non-error TOML case"
        );

        for key in sorted_omm_keys(&cases) {
            let case = &cases.test[key];

            if case.exception {
                assert!(
                    sgp4_from_case_omm(&case.omm_kvn).is_err(),
                    "case {key} ({}): expected from_omm_kvn_lines to return Err",
                    case.name
                );
                assert!(
                    from_omm_kvn_string(&case.omm_kvn).is_err(),
                    "case {key} ({}): expected from_omm_kvn_string to return Err",
                    case.name
                );
                continue;
            }

            let from_lines = sgp4_from_case_omm(&case.omm_kvn).unwrap_or_else(|err| {
                panic!(
                    "case {key} ({}): from_omm_kvn_lines failed: {err:?}",
                    case.name
                )
            });
            assert_gp_matches(key, &case.name, "from_omm_kvn_lines", &from_lines.gp, case);

            let from_string = from_omm_kvn_string(&case.omm_kvn).unwrap_or_else(|err| {
                panic!(
                    "case {key} ({}): from_omm_kvn_string failed: {err:?}",
                    case.name
                )
            });
            assert_eq!(
                from_string.len(),
                1,
                "case {key} ({}): expected one OMM from string",
                case.name
            );
            assert_gp_matches(
                key,
                &case.name,
                "from_omm_kvn_string",
                &from_string[0].gp,
                case,
            );

            let file_match = from_file
                .iter()
                .find(|s| s.gp.satellite_catalog_number == case.satellite_catalog_number)
                .unwrap_or_else(|| {
                    panic!(
                        "case {key} ({}): catalog {} missing from omm_parsing_cases.txt parse",
                        case.name, case.satellite_catalog_number
                    )
                });
            assert_gp_matches(key, &case.name, "from_omm_kvn_file", &file_match.gp, case);
        }
    }

    /// Compare XML OMM parsing with the TOML cases and the fixture file
    ///
    /// # Panics
    /// * If the fixture fails to parse or a field misses the expected value
    #[cfg(feature = "xml")]
    #[test]
    fn test_omm_xml_parsing_cases() {
        let cases = load_omm_parsing_cases();
        let from_xml_file = from_omm_xml_file("test/omm_parsing_cases.xml")
            .expect("could not parse test/omm_parsing_cases.xml");
        assert_omm_file_matches_cases(
            &cases,
            &from_xml_file,
            "from_omm_xml_file",
            "omm_parsing_cases.xml",
        );
    }

    /// Compare JSON OMM parsing with the TOML cases and the fixture file
    ///
    /// # Panics
    /// * If the fixture fails to parse or a field misses the expected value
    #[cfg(feature = "json")]
    #[test]
    fn test_omm_json_parsing_cases() {
        let cases = load_omm_parsing_cases();
        let from_json_file = from_omm_json_file("test/omm_parsing_cases.json")
            .expect("could not parse test/omm_parsing_cases.json");
        assert_omm_file_matches_cases(
            &cases,
            &from_json_file,
            "from_omm_json_file",
            "omm_parsing_cases.json",
        );
    }

    /// Compare CSV OMM parsing with the TOML cases and the fixture file
    ///
    /// # Panics
    /// * If the fixture fails to parse or a field misses the expected value
    #[cfg(feature = "csv")]
    #[test]
    fn test_omm_csv_parsing_cases() {
        let cases = load_omm_parsing_cases();
        let from_csv_file = from_omm_csv_file("test/omm_parsing_cases.csv")
            .expect("could not parse test/omm_parsing_cases.csv");
        assert_omm_file_matches_cases(
            &cases,
            &from_csv_file,
            "from_omm_csv_file",
            "omm_parsing_cases.csv",
        );
    }

    /// Compare two element sets field by field
    ///
    /// # Arguments
    /// * `label` - Record label, included in the failure message
    /// * `original` - Element set before export
    /// * `exported` - Element set parsed from the export
    ///
    /// # Panics
    /// * If any field differs
    fn assert_gp_eq(label: &str, original: &GenPerturbElementSet, exported: &GenPerturbElementSet) {
        assert_eq!(
            original.common_name, exported.common_name,
            "{label}: common_name"
        );
        assert_eq!(
            original.satellite_catalog_number, exported.satellite_catalog_number,
            "{label}: satellite_catalog_number"
        );
        assert_eq!(
            original.classification, exported.classification,
            "{label}: classification"
        );
        assert_eq!(
            original.international_designator, exported.international_designator,
            "{label}: international_designator"
        );
        assert_eq!(
            original.epoch_datetime.year, exported.epoch_datetime.year,
            "{label}: epoch year"
        );
        assert_eq!(
            original.epoch_datetime.month, exported.epoch_datetime.month,
            "{label}: epoch month"
        );
        assert_eq!(
            original.epoch_datetime.day, exported.epoch_datetime.day,
            "{label}: epoch day"
        );
        assert_eq!(
            original.epoch_datetime.hour, exported.epoch_datetime.hour,
            "{label}: epoch hour"
        );
        assert_eq!(
            original.epoch_datetime.minute, exported.epoch_datetime.minute,
            "{label}: epoch minute"
        );
        assert!(
            (original.epoch_datetime.second - exported.epoch_datetime.second).abs() < TLE_PARSE_TOL,
            "{label}: epoch second {} vs {}",
            original.epoch_datetime.second,
            exported.epoch_datetime.second
        );
        assert!(
            (original.first_derivative_of_mean_motion - exported.first_derivative_of_mean_motion)
                .abs()
                < TLE_PARSE_TOL,
            "{label}: n-dot"
        );
        assert!(
            (original.second_derivative_of_mean_motion - exported.second_derivative_of_mean_motion)
                .abs()
                < TLE_PARSE_TOL,
            "{label}: n-ddot"
        );
        assert!(
            (original.bstar - exported.bstar).abs() < TLE_PARSE_TOL,
            "{label}: bstar"
        );
        assert_eq!(
            original.ephemeris_type, exported.ephemeris_type,
            "{label}: ephemeris_type"
        );
        assert_eq!(
            original.element_set_number, exported.element_set_number,
            "{label}: element_set_number"
        );
        assert!(
            (original.inclination - exported.inclination).abs() < TLE_PARSE_TOL,
            "{label}: inclination"
        );
        assert!(
            (original.right_ascension_of_ascending_node
                - exported.right_ascension_of_ascending_node)
                .abs()
                < TLE_PARSE_TOL,
            "{label}: raan"
        );
        assert!(
            (original.eccentricity - exported.eccentricity).abs() < TLE_PARSE_TOL,
            "{label}: eccentricity"
        );
        assert!(
            (original.argument_of_perigee - exported.argument_of_perigee).abs() < TLE_PARSE_TOL,
            "{label}: argp"
        );
        assert!(
            (original.mean_anomaly - exported.mean_anomaly).abs() < TLE_PARSE_TOL,
            "{label}: mean anomaly"
        );
        assert!(
            (original.mean_motion - exported.mean_motion).abs() < TLE_PARSE_TOL,
            "{label}: mean motion"
        );
        assert_eq!(
            original.revolution_number_at_epoch, exported.revolution_number_at_epoch,
            "{label}: revolution_number_at_epoch"
        );
    }

    /// Create the export directory and return a path inside it
    ///
    /// # Arguments
    /// * `file_name` - File name under test/export
    ///
    /// # Returns
    /// * `String` - Path of the export file
    ///
    /// # Panics
    /// * If the directory cannot be created or the path is not valid UTF-8
    fn export_test_path(file_name: &str) -> String {
        let dir = std::path::Path::new("test/export");
        std::fs::create_dir_all(dir).expect("could not create test/export");
        dir.join(file_name)
            .to_str()
            .expect("export path is valid UTF-8")
            .to_string()
    }

    /// Export no propagators as an empty KVN string
    ///
    /// # Panics
    /// * If the export is not empty
    #[test]
    fn test_omm_kvn_export_empty() {
        // An empty slice should produce an empty KVN string
        let exported = to_omm_kvn_string(&[]);
        assert_eq!(exported, "");
    }

    /// Export KVN text and parse it back to the same element sets
    ///
    /// # Panics
    /// * If a record count or field differs
    #[test]
    fn test_omm_kvn_export_string_roundtrip() {
        // Parse the KVN test file, export, and parse the export
        let original = from_omm_kvn_file("test/omm_parsing_cases.txt")
            .expect("could not parse test/omm_parsing_cases.txt");
        let exported = to_omm_kvn_string(&original);
        let reparsed = from_omm_kvn_string(&exported).unwrap();

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Write KVN to a file and parse it back to the same element sets
    ///
    /// # Panics
    /// * If the file cannot be written or a field differs
    #[test]
    fn test_omm_kvn_export_file_roundtrip() {
        // Parse the KVN test file and write it back out
        let original = from_omm_kvn_file("test/omm_parsing_cases.txt")
            .expect("could not parse test/omm_parsing_cases.txt");
        let out_path = export_test_path("omm_kvn.txt");
        to_omm_kvn_file(&original, &out_path).expect("could not write exported KVN file");

        // Parse the written file and compare GP fields
        let reparsed = from_omm_kvn_file(&out_path).expect("could not parse exported KVN file");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Export no propagators as an XML document with no OMM records
    ///
    /// # Panics
    /// * If the exported document parses as a non-empty list
    #[cfg(feature = "xml")]
    #[test]
    fn test_omm_xml_export_empty() {
        // An empty slice should produce an NDM document with no OMM records
        let exported = to_omm_xml_string(&[]);
        let reparsed = from_omm_xml_string(&exported).unwrap();
        assert!(reparsed.is_empty());
    }

    /// Export XML text and parse it back to the same element sets
    ///
    /// # Panics
    /// * If a record count or field differs
    #[cfg(feature = "xml")]
    #[test]
    fn test_omm_xml_export_string_roundtrip() {
        // Parse the XML test file, export, and parse the export
        let original = from_omm_xml_file("test/omm_parsing_cases.xml")
            .expect("could not parse test/omm_parsing_cases.xml");
        let exported = to_omm_xml_string(&original);
        let reparsed = from_omm_xml_string(&exported).unwrap();

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Write XML to a file and parse it back to the same element sets
    ///
    /// # Panics
    /// * If the file cannot be written or a field differs
    #[cfg(feature = "xml")]
    #[test]
    fn test_omm_xml_export_file_roundtrip() {
        // Parse the XML test file and write it back out
        let original = from_omm_xml_file("test/omm_parsing_cases.xml")
            .expect("could not parse test/omm_parsing_cases.xml");
        let out_path = export_test_path("omm_xml.xml");
        to_omm_xml_file(&original, &out_path).expect("could not write exported XML file");

        // Parse the written file and compare GP fields
        let reparsed = from_omm_xml_file(&out_path).expect("could not parse exported XML file");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Export no propagators as an empty JSON array
    ///
    /// # Panics
    /// * If the exported array parses as a non-empty list
    #[cfg(feature = "json")]
    #[test]
    fn test_omm_json_export_empty() {
        // An empty slice should produce an empty JSON array
        let exported = to_omm_json_string(&[]);
        let reparsed = from_omm_json_string(&exported).unwrap();
        assert!(reparsed.is_empty());
    }

    /// Export JSON text and parse it back to the same element sets
    ///
    /// # Panics
    /// * If a record count or field differs
    #[cfg(feature = "json")]
    #[test]
    fn test_omm_json_export_string_roundtrip() {
        // Parse the JSON test file, export, and parse the export
        let original = from_omm_json_file("test/omm_parsing_cases.json")
            .expect("could not parse test/omm_parsing_cases.json");
        let exported = to_omm_json_string(&original);
        let reparsed = from_omm_json_string(&exported).unwrap();

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Write JSON to a file and parse it back to the same element sets
    ///
    /// # Panics
    /// * If the file cannot be written or a field differs
    #[cfg(feature = "json")]
    #[test]
    fn test_omm_json_export_file_roundtrip() {
        // Parse the JSON test file and write it back out
        let original = from_omm_json_file("test/omm_parsing_cases.json")
            .expect("could not parse test/omm_parsing_cases.json");
        let out_path = export_test_path("omm_json.json");
        to_omm_json_file(&original, &out_path).expect("could not write exported JSON file");

        // Parse the written file and compare GP fields
        let reparsed = from_omm_json_file(&out_path).expect("could not parse exported JSON file");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Export no propagators as a header-only CSV
    ///
    /// # Panics
    /// * If the exported CSV parses as a non-empty list
    #[cfg(feature = "csv")]
    #[test]
    fn test_omm_csv_export_empty() {
        // An empty slice should produce a header-only CSV
        let exported = to_omm_csv_string(&[]);
        let reparsed = from_omm_csv_string(&exported).unwrap();
        assert!(reparsed.is_empty());
    }

    /// Export CSV text and parse it back to the same element sets
    ///
    /// # Panics
    /// * If a record count or field differs
    #[cfg(feature = "csv")]
    #[test]
    fn test_omm_csv_export_string_roundtrip() {
        // Parse the CSV test file, export, and parse the export
        let original = from_omm_csv_file("test/omm_parsing_cases.csv")
            .expect("could not parse test/omm_parsing_cases.csv");
        let exported = to_omm_csv_string(&original);
        let reparsed = from_omm_csv_string(&exported).unwrap();

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Write CSV to a file and parse it back to the same element sets
    ///
    /// # Panics
    /// * If the file cannot be written or a field differs
    #[cfg(feature = "csv")]
    #[test]
    fn test_omm_csv_export_file_roundtrip() {
        // Parse the CSV test file and write it back out
        let original = from_omm_csv_file("test/omm_parsing_cases.csv")
            .expect("could not parse test/omm_parsing_cases.csv");
        let out_path = export_test_path("omm_csv.csv");
        to_omm_csv_file(&original, &out_path).expect("could not write exported CSV file");

        // Parse the written file and compare GP fields
        let reparsed = from_omm_csv_file(&out_path).expect("could not parse exported CSV file");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Export no propagators as an empty TLE string
    ///
    /// # Panics
    /// * If the export fails or is not empty
    #[test]
    fn test_tle_export_empty() {
        // An empty slice should produce an empty TLE string
        let exported = to_tle_string(&[]).expect("empty TLE export should succeed");
        assert_eq!(exported, "");
    }

    /// Export TLE text and parse it back to the same element sets
    ///
    /// # Panics
    /// * If a record count or field differs
    #[test]
    fn test_tle_export_string_roundtrip() {
        // Parse the TLE test file, export, and parse the export
        let original = from_tle_file("test/tle_parsing_cases.txt")
            .expect("could not read test/tle_parsing_cases.txt");
        let exported = to_tle_string(&original).expect("TLE export should succeed");
        let reparsed = from_tle_string(&exported).expect("exported TLE string should parse");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Write a TLE file and parse it back to the same element sets
    ///
    /// # Panics
    /// * If the file cannot be written or a field differs
    #[test]
    fn test_tle_export_file_roundtrip() {
        // Parse the TLE test file and write it back out
        let original = from_tle_file("test/tle_parsing_cases.txt")
            .expect("could not read test/tle_parsing_cases.txt");
        let out_path = export_test_path("tle.txt");
        to_tle_file(&original, &out_path).expect("TLE file export should succeed");

        // Parse the written file and compare GP fields
        let reparsed = from_tle_file(&out_path).expect("could not parse exported TLE file");

        assert_eq!(original.len(), reparsed.len());
        for (i, (a, b)) in original.iter().zip(reparsed.iter()).enumerate() {
            assert_gp_eq(&format!("record {i}"), &a.gp, &b.gp);
        }
    }

    /// Check exported TLE line length, checksum, and name length
    ///
    /// Data lines must be 69 characters and pass the checksum. Name lines must be
    /// non-empty and at most 24 characters.
    ///
    /// # Panics
    /// * If a line has the wrong length or a bad checksum
    #[test]
    fn test_tle_export_line_format() {
        // Exported data lines must be 69 characters and pass the checksum
        let original = from_tle_file("test/tle_parsing_cases.txt")
            .expect("could not read test/tle_parsing_cases.txt");
        let exported = to_tle_string(&original).expect("TLE export should succeed");

        for line in exported.lines() {
            if line.starts_with('1') || line.starts_with('2') {
                assert_eq!(line.len(), 69, "TLE data line length: {line}");
                assert!(tle_checksum(line), "TLE checksum: {line}");
            } else {
                assert!(!line.is_empty() && line.len() <= 24, "TLE name: {line}");
            }
        }
    }

    /// Use type defaults when a KVN field is missing
    ///
    /// Numbers default to 0, strings to empty, and classification to U.
    ///
    /// # Panics
    /// * If a missing field returns another default
    #[test]
    fn test_kvn_parse_defaults() {
        // Define an empty OMM record
        let lines: [&str; 0] = [];

        // Parse missing fields
        let name: String = kvn_parse(&lines, "OBJECT_NAME").unwrap();
        let catalog: i32 = kvn_parse(&lines, "NORAD_CAT_ID").unwrap();
        let rev: i64 = kvn_parse(&lines, "REV_AT_EPOCH").unwrap();
        let motion: f64 = kvn_parse(&lines, "MEAN_MOTION").unwrap();
        let class: char = kvn_parse(&lines, "CLASSIFICATION_TYPE").unwrap();

        // Assert numeric fields default to 0, strings to empty, classification to U
        assert_eq!(name, "");
        assert_eq!(catalog, 0);
        assert_eq!(rev, 0);
        assert_eq!(motion, 0.0);
        assert_eq!(class, 'U');
    }
}
