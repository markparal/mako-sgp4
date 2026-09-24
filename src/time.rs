//! Module for datetime representation and calendar conversions used by GP
//! parsing and SGP4 propagation (UTC Julian / Modified Julian dates).

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------

// -------
// Structs
// -------

/// A datetime structure
///
/// Represents a point in time with year, month, day, hour, minute, second
/// components, and an associated timezone.
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::{DateTime, Timezone};
///
/// // Define a UTC datetime (required by Julian-date conversions)
/// let datetime = DateTime {
///     year: 2000,
///     month: 1,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UTC,
/// };
///
/// // Assert the calendar components
/// assert_eq!(datetime.year, 2000);
/// assert_eq!(datetime.timezone, Timezone::UTC);
/// ```
///
/// # References
#[derive(Default, Debug, Clone, PartialEq, Copy)]
pub struct DateTime {
    /// The year
    pub year: i32,

    /// The month (1-12)
    pub month: i32,

    /// The day of month (1-31)
    pub day: i32,

    /// The hour (0-23)
    pub hour: i32,

    /// The minute (0-59)
    pub minute: i32,

    /// The second with fractional component (0.0-59.999...)
    pub second: f64,

    /// The timezone associated with this datetime
    pub timezone: Timezone,
}

// -----
// Enums
// -----

/// Date conversion errors
///
/// Errors that can occur during Julian-date and day-of-year conversions.
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::{DateError, DateTime, Timezone, utc2jday};
///
/// // Julian conversion requires UTC
/// let datetime = DateTime {
///     year: 2000,
///     month: 1,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UT1,
/// };
///
/// // Assert the non-UTC error
/// let err = utc2jday(&datetime).unwrap_err();
/// assert_eq!(err, DateError::DateNotUTC);
/// ```
///
/// # References
#[derive(Debug, Clone, PartialEq)]
pub enum DateError {
    /// The provided date is before October 10th, 1582 (Gregorian calendar adoption)
    DateTooEarly,

    /// The day of year is invalid (less than 1 or greater than 365/366)
    InvalidDayOfYear,

    /// The datetime is not in UTC
    DateNotUTC,
}

/// Timezone options for datetime representation
///
/// Represents the time scale used for the datetime.
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::Timezone;
///
/// // UTC is the default and the time scale used by GP / SGP4
/// let tz_utc = Timezone::UTC;
/// let tz_ut1 = Timezone::UT1;
///
/// assert_eq!(tz_utc, Timezone::default());
/// assert_ne!(tz_utc, tz_ut1);
/// ```
///
/// # References
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Timezone {
    /// Coordinated Universal Time (UTC)
    ///
    /// UTC is the primary time standard by which the world regulates clocks and time.
    /// It is within about 1 second of mean solar time at 0 deg longitude.
    #[default]
    UTC,

    /// Universal Time 1 (UT1)
    ///
    /// UT1 is a form of Universal Time that is directly related to the rotation of the Earth.
    /// It is based on the Earth's rotation and is used in astronomical calculations.
    /// UT1 differs from UTC by up to 0.9 seconds due to variations in Earth's rotation.
    UT1,
}

// ------
// Traits
// ------

// ---------
// Constants
// ---------

// ---------
// Functions
// ---------

/// Convert a datetime in UTC to Julian date (JD) format.
///
/// The Julian date is a continuous count of days since 4713-01-01 12:00:00 BCE (Julian calendar).
/// This function converts a UTC datetime to Julian date format, valid for any date after
/// October 10th, 1582 (Gregorian calendar adoption).
///
/// # Arguments
/// * `utc_datetime` - The datetime as a [`DateTime`] structure (in UTC)
///
/// # Returns
/// * `Ok((jd, jdfrac))` - Julian day (integer part) and day fraction
/// * `Err(DateError)` - If the datetime is invalid for conversion
///
/// # Errors
/// * `DateError::DateTooEarly` - If the date is before October 10th, 1582
/// * `DateError::DateNotUTC` - If the datetime is not in UTC
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::utc2jday;
/// use mako_sgp4::time::{DateTime, Timezone};
///
/// // Define J2000.0: 2000-01-01 12:00:00 UTC
/// let datetime = DateTime {
///     year: 2000,
///     month: 1,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UTC,
/// };
///
/// // Convert to Julian date (integer day and day fraction)
/// let (jd, jdfrac) = utc2jday(&datetime).unwrap();
///
/// // Assert JD 2451545.0
/// assert!(((jd + jdfrac) - 2_451_545.0).abs() < 1e-8);
/// assert!((0.0..1.0).contains(&jdfrac));
/// ```
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [Satellite Orbits by Montenbruck et al](https://link.springer.com/book/10.1007/978-3-642-58351-3)
pub fn utc2jday(utc_datetime: &DateTime) -> Result<(f64, f64), DateError> {
    // Calculate the MJD
    let (mjd, mjdfrac) = utc2mjday(utc_datetime)?;

    // Modify MJD to be JD
    let mut jd: f64 = mjd + 2400000.5;
    let mut jdfrac: f64 = mjdfrac;

    // Make JD whole
    if !(0.0..1.0).contains(&jdfrac) {
        jd += jdfrac.floor();
        jdfrac = jdfrac - jdfrac.floor();
    }

    Ok((jd, jdfrac))
}

/// Convert a datetime in UTC to Modified Julian date (MJD) format.
///
/// The Modified Julian date is a continuous count of days since 1858-11-17 00:00:00 CE.
/// MJD is related to Julian Date (JD) by: MJD = JD - 2400000.5
/// This function is valid for any date after October 10th, 1582 (Gregorian calendar adoption).
///
/// # Arguments
/// * `utc_datetime` - The datetime as a [`DateTime`] structure (in UTC)
///
/// # Returns
/// * `Ok((mjd, mjdfrac))` - Modified Julian day (integer part) and day fraction
/// * `Err(DateError)` - If the datetime is invalid for conversion
///
/// # Errors
/// * `DateError::DateTooEarly` - If the date is before October 10th, 1582
/// * `DateError::DateNotUTC` - If the datetime is not in UTC
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::utc2mjday;
/// use mako_sgp4::time::{DateTime, Timezone};
///
/// // Define J2000.0: 2000-01-01 12:00:00 UTC
/// let datetime = DateTime {
///     year: 2000,
///     month: 1,
///     day: 1,
///     hour: 12,
///     minute: 0,
///     second: 0.0,
///     timezone: Timezone::UTC,
/// };
///
/// // Convert to Modified Julian date (MJD = JD - 2400000.5)
/// let (mjd, mjdfrac) = utc2mjday(&datetime).unwrap();
///
/// // Assert MJD 51544.5
/// assert!(((mjd + mjdfrac) - 51_544.5).abs() < 1e-8);
/// assert!((0.0..1.0).contains(&mjdfrac));
/// ```
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [Satellite Orbits by Montenbruck et al](https://link.springer.com/book/10.1007/978-3-642-58351-3)
pub fn utc2mjday(utc_datetime: &DateTime) -> Result<(f64, f64), DateError> {
    // Verify that datetime is UTC
    if utc_datetime.timezone != Timezone::UTC {
        return Err(DateError::DateNotUTC);
    }

    // Verify date is after Oct 10th, 1582
    if utc_datetime.year < 1582
        || (utc_datetime.year == 1582 && utc_datetime.month < 10)
        || (utc_datetime.year == 1582 && utc_datetime.month == 10 && utc_datetime.day < 10)
    {
        return Err(DateError::DateTooEarly);
    }

    // Cast inputs as f64
    let year = utc_datetime.year as f64;
    let month = utc_datetime.month as f64;
    let day = utc_datetime.day as f64;
    let hour = utc_datetime.hour as f64;
    let minute = utc_datetime.minute as f64;
    let second = utc_datetime.second;

    // Modify month and year to account for leap years, start year in March instead of January
    let year_leap: f64;
    let month_leap: f64;
    if month <= 2. {
        year_leap = year - 1.;
        month_leap = month + 12.;
    } else {
        year_leap = year;
        month_leap = month;
    }

    // Account for leap days with B auxilary quantity
    let b_leap: f64 =
        (year_leap / 400.).floor() - (year_leap / 100.).floor() + (year_leap / 4.).floor();

    // Calculate the modified Julian date
    let mut mjd = 365. * year_leap - 679004. + b_leap + (30.6001 * (month_leap + 1.)).floor() + day;
    let mut mjdfrac = (second + minute * 60. + hour * 3600.) / 86400.;

    // Validate mjdfrac
    if !(0.0..1.0).contains(&mjdfrac) {
        mjd += mjdfrac.floor();
        mjdfrac = mjdfrac - mjdfrac.floor();
    }

    Ok((mjd, mjdfrac))
}

/// Convert a year and day of year to a UTC datetime
///
/// Converts a year and day of year (with fractional day) into a full UTC datetime.
/// The day of year is 1-based (1 = January 1st, 365/366 = December 31st).
///
/// # Arguments
/// * `year` - The year
/// * `dayofyr` - The day of year with fractional component (e.g., 123.5 = day 123 at 12:00:00 UTC)
///
/// # Returns
/// * `Ok(DateTime)` - The datetime as a [`DateTime`] structure (in UTC)
/// * `Err(DateError)` - If the day of year is out of range
///
/// # Errors
/// * `DateError::InvalidDayOfYear` - If the day of year is less than 1 or exceeds the number of days in the year
///
/// # Examples
/// ```rust
/// use mako_sgp4::time::{Timezone, dayofyr2utc};
///
/// // Define a TLE-style epoch: day 123.5 of 2024 is May 2 at 12:00:00 UTC
/// let datetime = dayofyr2utc(2024, 123.5).unwrap();
///
/// // Assert the calendar date, noon, and UTC
/// assert_eq!(datetime.year, 2024);
/// assert_eq!(datetime.month, 5);
/// assert_eq!(datetime.day, 2);
/// assert_eq!(datetime.hour, 12);
/// assert_eq!(datetime.timezone, Timezone::UTC);
/// ```
///
/// # References
/// - [Fundamentals of Astrodynamics and Applications by Vallado et al](https://celestrak.org/software/vallado-sw.php)
/// - [Satellite Orbits by Montenbruck et al](https://link.springer.com/book/10.1007/978-3-642-58351-3)
pub fn dayofyr2utc(year: i32, dayofyr: f64) -> Result<DateTime, DateError> {
    // Validate day of year is positive
    if dayofyr < 1.0 {
        return Err(DateError::InvalidDayOfYear);
    }

    // Check for leap year
    let is_leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);

    // Days per month (non-leap year)
    let days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    // Adjust February for leap year
    let max_days = if is_leap { 366 } else { 365 };

    // Extract integer and fractional parts
    let day_int = dayofyr.floor() as i32;

    // Validate day of year doesn't exceed days in year (check integer part)
    if day_int > max_days {
        return Err(DateError::InvalidDayOfYear);
    }
    let day_frac = dayofyr - day_int as f64;

    // Find which month the day falls in and calculate day of month
    let mut day_count = 0;
    let mut month = 1;
    let mut day = 1;

    for (idx, &days_in_month) in days_per_month.iter().enumerate() {
        let days_this_month = if idx == 1 && is_leap {
            29 // February in leap year
        } else {
            days_in_month
        };

        if day_int <= day_count + days_this_month {
            month = (idx + 1) as i32;
            day = day_int - day_count;
            break;
        }
        day_count += days_this_month;
    }

    // Convert fractional day to hours, minutes, seconds
    let total_seconds = day_frac * 86400.0;
    let hour = (total_seconds / 3600.0).floor() as i32;
    let remaining_seconds = total_seconds - (hour as f64 * 3600.0);
    let minute = (remaining_seconds / 60.0).floor() as i32;
    let second = remaining_seconds - (minute as f64 * 60.0);

    // Protect against rounding errors and handle overflow
    let mut final_second = second;
    let mut final_minute = minute;
    let mut final_hour = hour;
    let mut final_day = day;
    let mut final_month = month;
    let mut final_year = year;

    // Handle second overflow
    if final_second >= 60.0 {
        final_second -= 60.0;
        final_minute += 1;
    }

    // Handle minute overflow
    if final_minute >= 60 {
        final_minute -= 60;
        final_hour += 1;
    }

    // Handle hour overflow
    if final_hour >= 24 {
        final_hour -= 24;
        final_day += 1;
    }

    // Handle day overflow (check if day exceeds days in current month)
    let days_in_current_month = if final_month == 2 && is_leap {
        29 // February in leap year
    } else {
        days_per_month[(final_month - 1) as usize]
    };

    if final_day > days_in_current_month {
        final_day -= days_in_current_month;
        final_month += 1;
    }

    // Handle month overflow
    if final_month > 12 {
        final_month -= 12;
        final_year += 1;
    }

    // Store as DateTime
    let datetime = DateTime {
        year: final_year,
        month: final_month,
        day: final_day,
        hour: final_hour,
        minute: final_minute,
        second: final_second,
        timezone: Timezone::UTC,
    };

    Ok(datetime)
}

// ----------
// Unit Tests
// ----------

#[cfg(test)]
mod tests {
    use super::*;

    /// Convert known UTC datetimes to Julian dates
    ///
    /// Covers a 20th century date, a 21st century date, the year 2000, a date
    /// before Gregorian adoption, and a non-UTC timezone.
    ///
    /// # Panics
    /// * If a conversion misses the expected Julian date or the expected error
    #[test]
    fn test_utc2jday() {
        // Make a test date in 20th century
        let datetime1 = DateTime {
            year: 1959,
            month: 3,
            day: 25,
            hour: 12,
            minute: 34,
            second: 49.123,
            timezone: Timezone::UTC,
        };

        let (jd1, jdfrac1) = utc2jday(&datetime1).unwrap();
        let jd1_total = jd1 + jdfrac1;
        let jd1_expect = 2_436_653.024_179_664_4;
        assert!(
            (jd1_total - jd1_expect).abs() < 1e-8,
            "Julian Date Test failed: expected {}, got {}",
            jd1_expect,
            jd1_total
        );

        // Make a test date in 21st century
        let datetime2 = DateTime {
            year: 2026,
            month: 1,
            day: 23,
            hour: 0,
            minute: 10,
            second: 32.999,
            timezone: Timezone::UTC,
        };

        let (jd2, jdfrac2) = utc2jday(&datetime2).unwrap();
        let jd2_total = jd2 + jdfrac2;
        let jd2_expect = 2_461_063.507_326_377;
        assert!(
            (jd2_total - jd2_expect).abs() < 1e-8,
            "Julian Date Test failed: expected {}, got {}",
            jd2_expect,
            jd2_total
        );

        // Make a test date around year 2000
        let datetime3 = DateTime {
            year: 2000,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };

        let (jd3, jdfrac3) = utc2jday(&datetime3).unwrap();
        let jd3_total = jd3 + jdfrac3;
        let jd3_expect = 2451544.5;
        assert!(
            (jd3_total - jd3_expect).abs() < 1e-8,
            "Julian Date Test failed: expected {}, got {}",
            jd3_expect,
            jd3_total
        );

        // Test date too early (before Oct 10, 1582)
        let datetime4 = DateTime {
            year: 1582,
            month: 10,
            day: 9,
            hour: 12,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };

        let result = utc2jday(&datetime4);
        assert!(
            result.is_err(),
            "Should return error for date before Oct 10, 1582"
        );
        assert_eq!(result.unwrap_err(), DateError::DateTooEarly);

        // Test non-UTC datetime
        let datetime5 = DateTime {
            year: 1990,
            month: 10,
            day: 9,
            hour: 12,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UT1,
        };

        let result = utc2jday(&datetime5);
        assert!(result.is_err(), "Should return error for non-UTC date");
        assert_eq!(result.unwrap_err(), DateError::DateNotUTC);
    }

    /// Convert known UTC datetimes to Modified Julian dates
    ///
    /// Covers a 20th century date, a 21st century date, the year 2000, a date
    /// before Gregorian adoption, and a non-UTC timezone.
    ///
    /// # Panics
    /// * If a conversion misses the expected Modified Julian date or the expected error
    #[test]
    fn test_utc2mjday() {
        // Make a test date in 20th century
        let datetime1 = DateTime {
            year: 1997,
            month: 4,
            day: 2,
            hour: 16,
            minute: 12,
            second: 35.505,
            timezone: Timezone::UTC,
        };

        let (mjd1, mjdfrac1) = utc2mjday(&datetime1).unwrap();
        let mjd1_total = mjd1 + mjdfrac1;
        let mjd1_expect = 50_540.675_410_937_5;
        assert!(
            (mjd1_total - mjd1_expect).abs() < 1e-8,
            "Modified Julian Date Test failed: expected {}, got {}",
            mjd1_expect,
            mjd1_total
        );

        // Make a test date in 21st century
        let datetime2 = DateTime {
            year: 2013,
            month: 8,
            day: 12,
            hour: 2,
            minute: 49,
            second: 57.623,
            timezone: Timezone::UTC,
        };

        let (mjd2, mjdfrac2) = utc2mjday(&datetime2).unwrap();
        let mjd2_total = mjd2 + mjdfrac2;
        let mjd2_expect = 56_516.118_028_043_97;
        assert!(
            (mjd2_total - mjd2_expect).abs() < 1e-8,
            "Modified Julian Date Test failed: expected {}, got {}",
            mjd2_expect,
            mjd2_total
        );

        // Make a test date around year 2000
        let datetime3 = DateTime {
            year: 2000,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };

        let (mjd3, mjdfrac3) = utc2mjday(&datetime3).unwrap();
        let mjd3_total = mjd3 + mjdfrac3;
        let mjd3_expect = 51544.0;
        assert!(
            (mjd3_total - mjd3_expect).abs() < 1e-8,
            "Modified Julian Date Test failed: expected {}, got {}",
            mjd3_expect,
            mjd3_total
        );

        // Test date too early (before Oct 10, 1582)
        let datetime4 = DateTime {
            year: 1582,
            month: 10,
            day: 9,
            hour: 12,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UTC,
        };
        let result = utc2mjday(&datetime4);
        assert!(
            result.is_err(),
            "Should return error for date before Oct 10, 1582"
        );
        assert_eq!(result.unwrap_err(), DateError::DateTooEarly);

        // Test non-UTC datetime
        let datetime5 = DateTime {
            year: 1990,
            month: 10,
            day: 9,
            hour: 12,
            minute: 0,
            second: 0.0,
            timezone: Timezone::UT1,
        };

        let result = utc2mjday(&datetime5);
        assert!(result.is_err(), "Should return error for non-UTC date");
        assert_eq!(result.unwrap_err(), DateError::DateNotUTC);
    }

    /// Convert a year and day-of-year into a UTC datetime
    ///
    /// Checks a fractional day, a time near the next year, and day 366 of a leap year.
    ///
    /// # Panics
    /// * If a converted calendar field is wrong
    #[test]
    fn test_dayofyr_rounding() {
        // Test date in 20th century - Day 100.5 of 1959 (April 10, 1959 at 12:00:00)
        let datetime1 = dayofyr2utc(1959, 100.5).unwrap();
        assert_eq!(
            datetime1.year, 1959,
            "Day 100.5 of 1959: year should be 1959, got {}",
            datetime1.year
        );
        assert_eq!(
            datetime1.month, 4,
            "Day 100.5 of 1959: month should be 4 (April), got {}",
            datetime1.month
        );
        assert_eq!(
            datetime1.day, 10,
            "Day 100.5 of 1959: day should be 10, got {}",
            datetime1.day
        );
        assert_eq!(
            datetime1.hour, 12,
            "Day 100.5 of 1959: hour should be 12, got {}",
            datetime1.hour
        );
        assert_eq!(
            datetime1.minute, 0,
            "Day 100.5 of 1959: minute should be 0, got {}",
            datetime1.minute
        );
        assert!(
            (datetime1.second - 0.0).abs() < 1e-6,
            "Day 100.5 of 1959: second should be 0.0, got {}",
            datetime1.second
        );

        // Test date in 21st century - Day 200.75 of 2024 (July 18, 2024 at 18:00:00)
        let datetime2 = dayofyr2utc(2024, 200.75).unwrap();
        assert_eq!(
            datetime2.year, 2024,
            "Day 200.75 of 2024: year should be 2024, got {}",
            datetime2.year
        );
        assert_eq!(
            datetime2.month, 7,
            "Day 200.75 of 2024: month should be 7 (July), got {}",
            datetime2.month
        );
        assert_eq!(
            datetime2.day, 18,
            "Day 200.75 of 2024: day should be 18, got {}",
            datetime2.day
        );
        assert_eq!(
            datetime2.hour, 18,
            "Day 200.75 of 2024: hour should be 18, got {}",
            datetime2.hour
        );
        assert_eq!(
            datetime2.minute, 0,
            "Day 200.75 of 2024: minute should be 0, got {}",
            datetime2.minute
        );
        assert!(
            (datetime2.second - 0.0).abs() < 1e-6,
            "Day 200.75 of 2024: second should be 0.0, got {}",
            datetime2.second
        );

        // Test year rollover - Day 365.9999999999 of 2023 (very close to midnight of 2024)
        let datetime3 = dayofyr2utc(2023, 365.9999999999).unwrap();
        assert_eq!(
            datetime3.year, 2023,
            "Day 365.9999999999 of 2023: year should be 2023, got {}",
            datetime3.year
        );
        assert_eq!(
            datetime3.month, 12,
            "Day 365.9999999999 of 2023: month should be 12 (December), got {}",
            datetime3.month
        );
        assert_eq!(
            datetime3.day, 31,
            "Day 365.9999999999 of 2023: day should be 31, got {}",
            datetime3.day
        );
        assert_eq!(
            datetime3.hour, 23,
            "Day 365.9999999999 of 2023: hour should be 23, got {}",
            datetime3.hour
        );
        assert_eq!(
            datetime3.minute, 59,
            "Day 365.9999999999 of 2023: minute should be 59, got {}",
            datetime3.minute
        );
        assert!(
            datetime3.second >= 59.0 && datetime3.second < 60.0,
            "Day 365.9999999999 of 2023: second should be between 59.0 and 60.0, got {}",
            datetime3.second
        );

        // Test leap year day of year - Day 366 of 2024 (December 31, 2024)
        let datetime4 = dayofyr2utc(2024, 366.0).unwrap();
        assert_eq!(
            datetime4.year, 2024,
            "Day 366 of 2024: year should be 2024, got {}",
            datetime4.year
        );
        assert_eq!(
            datetime4.month, 12,
            "Day 366 of 2024: month should be 12 (December), got {}",
            datetime4.month
        );
        assert_eq!(
            datetime4.day, 31,
            "Day 366 of 2024: day should be 31, got {}",
            datetime4.day
        );
        assert_eq!(
            datetime4.hour, 0,
            "Day 366 of 2024: hour should be 0, got {}",
            datetime4.hour
        );
        assert_eq!(
            datetime4.minute, 0,
            "Day 366 of 2024: minute should be 0, got {}",
            datetime4.minute
        );
        assert!(
            (datetime4.second - 0.0).abs() < 1e-6,
            "Day 366 of 2024: second should be 0.0, got {}",
            datetime4.second
        );

        // Test leap year with fractional day - Day 60.5 of 2024 (February 29, 2024 at 12:00:00)
        let datetime5 = dayofyr2utc(2024, 60.5).unwrap();
        assert_eq!(
            datetime5.year, 2024,
            "Day 60.5 of 2024: year should be 2024, got {}",
            datetime5.year
        );
        assert_eq!(
            datetime5.month, 2,
            "Day 60.5 of 2024: month should be 2 (February), got {}",
            datetime5.month
        );
        assert_eq!(
            datetime5.day, 29,
            "Day 60.5 of 2024: day should be 29 (leap day), got {}",
            datetime5.day
        );
        assert_eq!(
            datetime5.hour, 12,
            "Day 60.5 of 2024: hour should be 12, got {}",
            datetime5.hour
        );
        assert_eq!(
            datetime5.minute, 0,
            "Day 60.5 of 2024: minute should be 0, got {}",
            datetime5.minute
        );
        assert!(
            (datetime5.second - 0.0).abs() < 1e-6,
            "Day 60.5 of 2024: second should be 0.0, got {}",
            datetime5.second
        );

        // Test day of year too high - Day 366.1 of non-leap year 2023
        let result = dayofyr2utc(2023, 366.1);
        assert!(
            result.is_err(),
            "Day 366.1 of 2023 (non-leap year): should return error for day of year exceeding 365, got Ok({:?})",
            result
        );
        assert_eq!(
            result.unwrap_err(),
            DateError::InvalidDayOfYear,
            "Day 366.1 of 2023: error should be InvalidDayOfYear"
        );

        // Test day of year too high - Day 367 of leap year 2024
        let result = dayofyr2utc(2024, 367.0);
        assert!(
            result.is_err(),
            "Day 367 of 2024 (leap year): should return error for day of year exceeding 366, got Ok({:?})",
            result
        );
        assert_eq!(
            result.unwrap_err(),
            DateError::InvalidDayOfYear,
            "Day 367 of 2024: error should be InvalidDayOfYear"
        );

        // Test day of year too low - Day 0.5
        let result = dayofyr2utc(2024, 0.5);
        assert!(
            result.is_err(),
            "Day 0.5 of 2024: should return error for day of year less than 1, got Ok({:?})",
            result
        );
        assert_eq!(
            result.unwrap_err(),
            DateError::InvalidDayOfYear,
            "Day 0.5 of 2024: error should be InvalidDayOfYear"
        );

        // Test day of year too low - Day 0.0
        let result = dayofyr2utc(2024, 0.0);
        assert!(
            result.is_err(),
            "Day 0.0 of 2024: should return error for day of year equal to 0, got Ok({:?})",
            result
        );
        assert_eq!(
            result.unwrap_err(),
            DateError::InvalidDayOfYear,
            "Day 0.0 of 2024: error should be InvalidDayOfYear"
        );
    }
}
