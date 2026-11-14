//! Module of floating-point math functions backed by `libm`.
//!
//! `f64` methods such as `sin` and `sqrt` live in `std`, not `core`. These free
//! functions call `libm` instead, so the crate builds without `std` and gives
//! bit-identical results on every platform.

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------

// ---------
// Functions
// ---------

/// Sine of `x` \[rad\]
#[inline]
pub(crate) fn sin(x: f64) -> f64 {
    libm::sin(x)
}

/// Cosine of `x` \[rad\]
#[inline]
pub(crate) fn cos(x: f64) -> f64 {
    libm::cos(x)
}

/// Four-quadrant arctangent of `y / x` \[rad\]
#[inline]
pub(crate) fn atan2(y: f64, x: f64) -> f64 {
    libm::atan2(y, x)
}

/// Square root of `x`
#[inline]
pub(crate) fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}

/// `x` raised to the floating-point power `y`
#[inline]
pub(crate) fn powf(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

/// `x` raised to the integer power `n`
///
/// Uses square-and-multiply, the same algorithm as the compiler's `powi` builtin.
#[inline]
pub(crate) fn powi(x: f64, n: i32) -> f64 {
    let reciprocal = n < 0;
    let mut exponent = n.unsigned_abs();
    let mut base = x;
    let mut result = 1.0;
    loop {
        if exponent & 1 != 0 {
            result *= base;
        }
        exponent >>= 1;
        if exponent == 0 {
            break;
        }
        base *= base;
    }
    if reciprocal { 1.0 / result } else { result }
}

/// Base-10 logarithm of `x`
#[inline]
pub(crate) fn log10(x: f64) -> f64 {
    libm::log10(x)
}

/// Absolute value of `x`
#[inline]
pub(crate) fn abs(x: f64) -> f64 {
    libm::fabs(x)
}

/// Largest integer less than or equal to `x`
#[inline]
pub(crate) fn floor(x: f64) -> f64 {
    libm::floor(x)
}

/// Nearest integer to `x`, rounding half-way cases away from zero
#[inline]
pub(crate) fn round(x: f64) -> f64 {
    libm::round(x)
}

/// Least non-negative remainder of `x / y`
///
/// Matches `f64::rem_euclid` in `std`.
#[inline]
pub(crate) fn rem_euclid(x: f64, y: f64) -> f64 {
    let r = x % y;
    if r < 0.0 { r + abs(y) } else { r }
}

// ----------
// Unit Tests
// ----------

#[cfg(test)]
mod tests {
    use super::*;

    /// `powi` matches repeated multiplication and handles negative and zero exponents
    ///
    /// # Panics
    /// * If a result differs from the expected value
    #[test]
    fn test_powi() {
        let x = 1.234_567_f64;
        assert_eq!(powi(x, 0), 1.0);
        assert_eq!(powi(x, 1), x);
        assert_eq!(powi(x, 2), x * x);
        assert_eq!(powi(x, 3), x * (x * x));
        assert_eq!(powi(x, -2), 1.0 / (x * x));
        assert_eq!(powi(2.0, 10), 1024.0);
    }

    /// `rem_euclid` returns a non-negative remainder like `std`
    ///
    /// # Panics
    /// * If a remainder is negative or wrong
    #[test]
    fn test_rem_euclid() {
        assert_eq!(rem_euclid(7.0, 4.0), 3.0);
        assert_eq!(rem_euclid(-1.0, 4.0), 3.0);
        assert_eq!(rem_euclid(-1.0, -4.0), 3.0);
        assert_eq!(rem_euclid(0.0, 4.0), 0.0);
    }
}
