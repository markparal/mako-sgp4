//! Module of small dense linear algebra for least-squares fitting.
//!
//! Householder QR reduces a tall weighted Jacobian to a square triangular
//! factor without forming its normal matrix, a one-sided Jacobi SVD factors
//! that triangle, and Cholesky factors noise covariances for whitening. All
//! sizes are compile-time constants, so every factor lives on the stack.

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------
use crate::math::{abs, sqrt};

// -------
// Structs
// -------

/// Householder QR of a least-squares problem
///
/// For the problem min ||b - A*x|| with A of size m x N, the factorization
/// A = Q*R gives the same solution as R*x = (Q^T*b)[..N]. Q itself is not kept.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Qr<const N: usize> {
    /// Upper-triangular factor R \[N x N\]
    pub r: [[f64; N]; N],

    /// First N entries of Q^T*b
    pub qtb: [f64; N],

    /// Sum of squares of the remaining entries of Q^T*b, which is the residual
    /// sum of squares of the least-squares solution
    pub residual_sq: f64,
}

/// Singular value decomposition of a square matrix
///
/// Holds the factors of A = U*S*V^T, with S the diagonal matrix of singular
/// values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Svd<const N: usize> {
    /// Left singular vectors, one per column; a column is zero where its
    /// singular value is zero
    pub u: [[f64; N]; N],

    /// Singular values, sorted from largest to smallest
    pub sigma: [f64; N],

    /// Right singular vectors, one per column
    pub v: [[f64; N]; N],
}

// ---------
// Constants
// ---------

/// Sweep cap for the Jacobi SVD
///
/// A 7 x 7 matrix converges in well under 10 sweeps; the cap only guards
/// against a matrix that never settles.
const MAX_JACOBI_SWEEPS: usize = 60;

// ---------
// Functions
// ---------

/// Reduce a least-squares problem with Householder reflections
///
/// Factors the m x N matrix `a` as Q*R and applies Q^T to `b` at the same
/// time, so Q is never stored. Both inputs are overwritten. A column that is
/// already zero on and below the diagonal is left as it is, giving a zero on
/// the diagonal of R.
///
/// # Arguments
/// * `a` - The matrix, one array per row \[m x N\]; overwritten
/// * `b` - The right-hand side \[m\]; overwritten with Q^T*b
///
/// # Returns
/// * `qr` - The triangular factor R, the first N entries of Q^T*b, and the
///   residual sum of squares
///
/// # Panics
/// * If `a` and `b` have different numbers of rows
///
/// # References
/// - G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Section 5.2
pub(crate) fn householder_qr<const N: usize>(a: &mut [[f64; N]], b: &mut [f64]) -> Qr<N> {
    assert_eq!(
        a.len(),
        b.len(),
        "matrix and right-hand side row counts differ"
    );
    let m = a.len();

    for k in 0..N.min(m) {
        // Norm of the column on and below the diagonal
        let norm = sqrt((k..m).map(|i| a[i][k] * a[i][k]).sum());
        if norm == 0.0 {
            continue;
        }

        // Reflect x onto alpha*e0 with alpha = -sign(x0)*||x||, so
        // v0 = x0 - alpha never cancels. Then v = (v0, x1, x2, ...) and
        // ||v||^2 = 2*||x||*(||x|| + |x0|).
        let x0 = a[k][k];
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        let v0 = x0 - alpha;
        let vtv = 2.0 * norm * (norm + abs(x0));

        // Apply I - 2*v*v^T/||v||^2 to the remaining columns and to b
        for j in (k + 1)..N {
            let s = v0 * a[k][j] + ((k + 1)..m).map(|i| a[i][k] * a[i][j]).sum::<f64>();
            let f = 2.0 * s / vtv;
            a[k][j] -= f * v0;
            for row in &mut a[(k + 1)..] {
                row[j] -= f * row[k];
            }
        }
        let s = v0 * b[k] + ((k + 1)..m).map(|i| a[i][k] * b[i]).sum::<f64>();
        let f = 2.0 * s / vtv;
        b[k] -= f * v0;
        for (bi, row) in b[(k + 1)..].iter_mut().zip(&a[(k + 1)..]) {
            *bi -= f * row[k];
        }

        a[k][k] = alpha;
    }

    // Copy out the upper triangle and split Q^T*b into its kept part and residual
    let n_kept = N.min(m);
    let mut r = [[0.0; N]; N];
    for (i, row) in r.iter_mut().enumerate().take(n_kept) {
        row[i..].copy_from_slice(&a[i][i..]);
    }
    let mut qtb = [0.0; N];
    qtb[..n_kept].copy_from_slice(&b[..n_kept]);
    let residual_sq = b[n_kept..].iter().map(|x| x * x).sum();

    Qr {
        r,
        qtb,
        residual_sq,
    }
}

/// Singular value decomposition of a square matrix by one-sided Jacobi rotations
///
/// Rotates pairs of columns of `a` until every pair is orthogonal to working
/// precision; the column norms are then the singular values. One-sided Jacobi
/// gives the small singular values to high relative accuracy, which the rank
/// test and the covariance rely on.
///
/// # Arguments
/// * `a` - The matrix \[N x N\]
///
/// # Returns
/// * `svd` - U, S and V with A = U*S*V^T, singular values in descending order
///
/// # References
/// - J. Demmel and K. Veselic, "Jacobi's Method is More Accurate than QR,"
///   *SIAM J. Matrix Anal. Appl.* 13(4), 1992
/// - G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Section 8.6
pub(crate) fn jacobi_svd<const N: usize>(a: &[[f64; N]; N]) -> Svd<N> {
    let mut w = *a;
    let mut v = identity::<N>();

    for _ in 0..MAX_JACOBI_SWEEPS {
        let mut rotated = false;
        for p in 0..N {
            for q in (p + 1)..N {
                let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
                for row in &w {
                    alpha += row[p] * row[p];
                    beta += row[q] * row[q];
                    gamma += row[p] * row[q];
                }
                // Columns already orthogonal to working precision, or one is zero
                if abs(gamma) <= f64::EPSILON * sqrt(alpha * beta) {
                    continue;
                }
                rotated = true;

                // Rotation angle that makes columns p and q orthogonal
                let zeta = (beta - alpha) / (2.0 * gamma);
                let sign = if zeta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (abs(zeta) + sqrt(1.0 + zeta * zeta));
                let c = 1.0 / sqrt(1.0 + t * t);
                let s = c * t;
                rotate_columns(&mut w, p, q, c, s);
                rotate_columns(&mut v, p, q, c, s);
            }
        }
        if !rotated {
            break;
        }
    }

    // Column norms are the singular values; normalized columns are U
    let mut sigma = [0.0; N];
    for (j, s) in sigma.iter_mut().enumerate() {
        *s = sqrt(w.iter().map(|row| row[j] * row[j]).sum());
    }
    let mut u = [[0.0; N]; N];
    for (urow, wrow) in u.iter_mut().zip(&w) {
        for ((uij, wij), &sj) in urow.iter_mut().zip(wrow).zip(&sigma) {
            if sj > 0.0 {
                *uij = wij / sj;
            }
        }
    }

    // Sort descending, carrying the columns of U and V
    let mut order = [0usize; N];
    for (i, o) in order.iter_mut().enumerate() {
        *o = i;
    }
    order.sort_by(|&i, &j| sigma[j].total_cmp(&sigma[i]));
    let mut sorted = Svd {
        u: [[0.0; N]; N],
        sigma: [0.0; N],
        v: [[0.0; N]; N],
    };
    for (dst, &src) in order.iter().enumerate() {
        sorted.sigma[dst] = sigma[src];
        for i in 0..N {
            sorted.u[i][dst] = u[i][src];
            sorted.v[i][dst] = v[i][src];
        }
    }
    sorted
}

/// Project a right-hand side onto the left singular vectors
///
/// Computes c = U^T*y, the coordinates used by [`damped_step`].
///
/// # Arguments
/// * `svd` - The decomposition
/// * `y` - The vector, usually the `qtb` of a [`Qr`]
///
/// # Returns
/// * `c` - The coordinates of `y` along each left singular vector
pub(crate) fn svd_project<const N: usize>(svd: &Svd<N>, y: &[f64; N]) -> [f64; N] {
    let mut c = [0.0; N];
    for (j, cj) in c.iter_mut().enumerate() {
        *cj = svd.u.iter().zip(y).map(|(urow, yi)| urow[j] * yi).sum();
    }
    c
}

/// Damped least-squares step from an SVD
///
/// Computes delta(mu) = sum over i of sigma_i/(sigma_i^2 + mu)*c_i*v_i, which
/// solves (A^T*A + mu*I)*delta = A^T*b from the factors alone, so a new mu
/// after a rejected step costs no new factorization. With mu = 0 this is the
/// Gauss-Newton (pseudo-inverse) step, and a zero singular value contributes
/// nothing rather than an infinite term.
///
/// # Arguments
/// * `svd` - The decomposition
/// * `c` - The projected right-hand side from [`svd_project`]
/// * `mu` - The damping parameter, 0 or greater
///
/// # Returns
/// * `delta` - The step
pub(crate) fn damped_step<const N: usize>(svd: &Svd<N>, c: &[f64; N], mu: f64) -> [f64; N] {
    let mut delta = [0.0; N];
    for (j, (&s, &cj)) in svd.sigma.iter().zip(c).enumerate() {
        let denom = s * s + mu;
        if denom == 0.0 {
            continue;
        }
        let weight = s / denom * cj;
        for (d, vrow) in delta.iter_mut().zip(&svd.v) {
            *d += weight * vrow[j];
        }
    }
    delta
}

/// Cholesky factorization of a symmetric positive definite matrix
///
/// Computes the lower-triangular L with A = L*L^T, reading only the lower
/// triangle of `a`.
///
/// # Arguments
/// * `a` - The matrix \[N x N\]
///
/// # Returns
/// * `Some(l)` - The lower-triangular factor
/// * `None` - If `a` isn't positive definite or holds a non-finite value
pub(crate) fn cholesky<const N: usize>(a: &[[f64; N]; N]) -> Option<[[f64; N]; N]> {
    let mut l = [[0.0; N]; N];
    for j in 0..N {
        let pivot = a[j][j] - (0..j).map(|k| l[j][k] * l[j][k]).sum::<f64>();
        // A NaN pivot isn't finite, so this rejects it too
        if !pivot.is_finite() || pivot <= 0.0 {
            return None;
        }
        let d = sqrt(pivot);
        l[j][j] = d;
        for i in (j + 1)..N {
            let s = a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if !s.is_finite() {
                return None;
            }
            l[i][j] = s / d;
        }
    }
    Some(l)
}

/// Solve a lower-triangular system by forward substitution
///
/// Solves L*x = y. Whitening a residual y with noise covariance R = L*L^T is
/// x = L^-1*y.
///
/// # Arguments
/// * `l` - The lower-triangular factor from [`cholesky`]
/// * `y` - The right-hand side
///
/// # Returns
/// * `x` - The solution
pub(crate) fn solve_lower<const N: usize>(l: &[[f64; N]; N], y: &[f64; N]) -> [f64; N] {
    let mut x = [0.0; N];
    for i in 0..N {
        let s = y[i] - (0..i).map(|k| l[i][k] * x[k]).sum::<f64>();
        x[i] = s / l[i][i];
    }
    x
}

/// Identity matrix
///
/// # Returns
/// * `m` - The N x N identity matrix
fn identity<const N: usize>() -> [[f64; N]; N] {
    let mut m = [[0.0; N]; N];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

/// Rotate two columns of a matrix in place
///
/// Applies p = c*p - s*q and q = s*p + c*q to columns p and q, using the old
/// values on the right of both.
///
/// # Arguments
/// * `m` - The matrix \[N x N\]; overwritten
/// * `p` - Index of the first column
/// * `q` - Index of the second column
/// * `c` - Cosine of the rotation angle
/// * `s` - Sine of the rotation angle
fn rotate_columns<const N: usize>(m: &mut [[f64; N]; N], p: usize, q: usize, c: f64, s: f64) {
    for row in m.iter_mut() {
        let (mp, mq) = (row[p], row[q]);
        row[p] = c * mp - s * mq;
        row[q] = s * mp + c * mq;
    }
}

// ----------
// Unit Tests
// ----------

// Matrix checks read most clearly with explicit indices
#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;
    use alloc::vec::Vec;

    /// Next value of a SplitMix64 generator, scaled to [-1, 1)
    ///
    /// A fixed generator in the test module keeps fixtures from shifting with a
    /// crate version.
    ///
    /// # Arguments
    /// * `state` - The generator state; advanced
    ///
    /// # Returns
    /// * `x` - A uniform value in [-1, 1)
    fn uniform(state: &mut u64) -> f64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }

    /// Random matrix with entries in [-1, 1)
    ///
    /// # Arguments
    /// * `rng` - The generator state
    /// * `m` - Number of rows
    ///
    /// # Returns
    /// * `a` - The m x N matrix
    fn random_matrix<const N: usize>(rng: &mut u64, m: usize) -> Vec<[f64; N]> {
        (0..m)
            .map(|_| core::array::from_fn(|_| uniform(rng)))
            .collect()
    }

    /// Random orthogonal matrix
    ///
    /// Builds Q = H0*H1*... from the Householder reflections that reduce a
    /// random square matrix.
    ///
    /// # Arguments
    /// * `rng` - The generator state
    ///
    /// # Returns
    /// * `q` - The N x N orthogonal matrix
    fn random_orthogonal<const N: usize>(rng: &mut u64) -> [[f64; N]; N] {
        let mut a = random_matrix::<N>(rng, N);
        let mut q = identity::<N>();
        for k in 0..N {
            let norm = sqrt((k..N).map(|i| a[i][k] * a[i][k]).sum());
            let x0 = a[k][k];
            let alpha = if x0 >= 0.0 { -norm } else { norm };
            let mut v = [0.0; N];
            v[k] = x0 - alpha;
            for i in (k + 1)..N {
                v[i] = a[i][k];
            }
            let vtv: f64 = v.iter().map(|x| x * x).sum();
            for j in 0..N {
                let s: f64 = (0..N).map(|i| v[i] * a[i][j]).sum();
                for i in 0..N {
                    a[i][j] -= 2.0 * s / vtv * v[i];
                }
            }
            for row in q.iter_mut() {
                let s: f64 = (0..N).map(|i| row[i] * v[i]).sum();
                for i in 0..N {
                    row[i] -= 2.0 * s / vtv * v[i];
                }
            }
        }
        q
    }

    /// Random tall matrix with chosen singular values
    ///
    /// Builds U*diag(sigma)*V^T, with U the orthonormal columns Q = A*R^-1 of a
    /// random tall matrix and V a random orthogonal matrix.
    ///
    /// # Arguments
    /// * `rng` - The generator state
    /// * `m` - Number of rows, at least N
    /// * `sigma` - The singular values
    ///
    /// # Returns
    /// * `a` - The m x N matrix
    fn matrix_with_singular_values<const N: usize>(
        rng: &mut u64,
        m: usize,
        sigma: &[f64; N],
    ) -> Vec<[f64; N]> {
        let tall = random_matrix::<N>(rng, m);
        let qr = householder_qr(&mut tall.clone(), &mut alloc::vec![0.0; m]);
        let u: Vec<[f64; N]> = tall
            .iter()
            .map(|row| {
                // Solve x*R = row for x, the row times R^-1
                let mut x = [0.0; N];
                for j in 0..N {
                    let s = row[j] - (0..j).map(|k| x[k] * qr.r[k][j]).sum::<f64>();
                    x[j] = s / qr.r[j][j];
                }
                x
            })
            .collect();
        let v = random_orthogonal::<N>(rng);
        u.iter()
            .map(|urow| {
                core::array::from_fn(|j| (0..N).map(|k| urow[k] * sigma[k] * v[j][k]).sum())
            })
            .collect()
    }

    /// Solve an upper-triangular system by back substitution
    ///
    /// # Arguments
    /// * `r` - The upper-triangular matrix
    /// * `y` - The right-hand side
    ///
    /// # Returns
    /// * `x` - The solution of R*x = y
    fn solve_upper<const N: usize>(r: &[[f64; N]; N], y: &[f64; N]) -> [f64; N] {
        let mut x = [0.0; N];
        for i in (0..N).rev() {
            let s = y[i] - ((i + 1)..N).map(|k| r[i][k] * x[k]).sum::<f64>();
            x[i] = s / r[i][i];
        }
        x
    }

    /// Largest absolute difference between two vectors
    ///
    /// # Arguments
    /// * `x` - The first vector
    /// * `y` - The second vector
    ///
    /// # Returns
    /// * `d` - The largest elementwise absolute difference
    fn max_diff<const N: usize>(x: &[f64; N], y: &[f64; N]) -> f64 {
        x.iter().zip(y).map(|(a, b)| abs(a - b)).fold(0.0, f64::max)
    }

    /// Euclidean norm of a vector
    ///
    /// # Arguments
    /// * `x` - The vector
    ///
    /// # Returns
    /// * `n` - The norm
    fn norm<const N: usize>(x: &[f64; N]) -> f64 {
        sqrt(x.iter().map(|v| v * v).sum())
    }

    /// QR reproduces A^T*A and the least-squares solution and residual of a tall system
    ///
    /// # Panics
    /// * If R^T*R differs from A^T*A, R isn't upper triangular, or the solution
    ///   or residual differs from the direct computation
    #[test]
    fn test_qr_least_squares() {
        let mut rng = 1;
        let a = random_matrix::<7>(&mut rng, 40);
        let b: Vec<f64> = (0..40).map(|_| uniform(&mut rng)).collect();
        let qr = householder_qr(&mut a.clone(), &mut b.clone());

        // R^T*R = A^T*A, and R is upper triangular
        for i in 0..7 {
            for j in 0..7 {
                let ata: f64 = a.iter().map(|row| row[i] * row[j]).sum();
                let rtr: f64 = (0..7).map(|k| qr.r[k][i] * qr.r[k][j]).sum();
                assert!(abs(ata - rtr) < 1e-12 * 40.0, "A^T*A[{i}][{j}]");
            }
            for j in 0..i {
                assert_eq!(qr.r[i][j], 0.0, "R below the diagonal");
            }
        }

        // The solution satisfies the normal equations, and the residual matches
        let x = solve_upper(&qr.r, &qr.qtb);
        let resid: Vec<f64> = a
            .iter()
            .zip(&b)
            .map(|(row, bi)| bi - (0..7).map(|k| row[k] * x[k]).sum::<f64>())
            .collect();
        for j in 0..7 {
            let g: f64 = a.iter().zip(&resid).map(|(row, r)| row[j] * r).sum();
            assert!(abs(g) < 1e-12, "A^T*(b - A*x)[{j}] = {g}");
        }
        let rss: f64 = resid.iter().map(|r| r * r).sum();
        assert!(abs(rss - qr.residual_sq) < 1e-12 * rss.max(1.0));
    }

    /// A square system, as from a single state with the prior row, has zero residual
    ///
    /// # Panics
    /// * If the square solve doesn't reproduce the exact solution
    #[test]
    fn test_qr_square() {
        let mut rng = 2;
        let a = random_matrix::<7>(&mut rng, 7);
        let x_true: [f64; 7] = core::array::from_fn(|_| uniform(&mut rng));
        let b: Vec<f64> = a
            .iter()
            .map(|row| (0..7).map(|k| row[k] * x_true[k]).sum())
            .collect();
        let qr = householder_qr(&mut a.clone(), &mut b.clone());
        assert_eq!(qr.residual_sq, 0.0);
        assert!(max_diff(&solve_upper(&qr.r, &qr.qtb), &x_true) < 1e-12);
    }

    /// QR with SVD solves an ill-conditioned system to about kappa*eps, while the
    /// normal equations lose about twice the digits
    ///
    /// # Panics
    /// * If the QR/SVD solution misses its accuracy bound, or the normal
    ///   equations come within a factor of 1000 of it
    #[test]
    fn test_qr_svd_ill_conditioned() {
        let mut rng = 3;
        // Condition number 1e7, roughly a long GEO arc after column scaling
        let sigma = [1.0, 0.3, 0.1, 1e-2, 1e-4, 1e-6, 1e-7];
        let a = matrix_with_singular_values(&mut rng, 60, &sigma);
        let x_true: [f64; 7] = core::array::from_fn(|_| uniform(&mut rng));
        let b: Vec<f64> = a
            .iter()
            .map(|row| (0..7).map(|k| row[k] * x_true[k]).sum())
            .collect();

        let qr = householder_qr(&mut a.clone(), &mut b.clone());
        let svd = jacobi_svd(&qr.r);
        let x = damped_step(&svd, &svd_project(&svd, &qr.qtb), 0.0);
        let qr_err = max_diff(&x, &x_true) / norm(&x_true);
        assert!(qr_err < 1e-7, "QR/SVD relative error {qr_err:e}");

        // Normal equations A^T*A*x = A^T*b, solved by Cholesky
        let mut ata = [[0.0; 7]; 7];
        let mut atb = [0.0; 7];
        for (row, bi) in a.iter().zip(&b) {
            for i in 0..7 {
                atb[i] += row[i] * bi;
                for j in 0..7 {
                    ata[i][j] += row[i] * row[j];
                }
            }
        }
        let l = cholesky(&ata).expect("A^T*A is positive definite at condition number 1e7");
        let y = solve_lower(&l, &atb);
        let mut lt = [[0.0; 7]; 7];
        for i in 0..7 {
            for j in 0..7 {
                lt[i][j] = l[j][i];
            }
        }
        let x_ne = solve_upper(&lt, &y);
        let ne_err = max_diff(&x_ne, &x_true) / norm(&x_true);
        // Measured: QR/SVD about 1e-10, normal equations about 8e-4
        assert!(
            ne_err > 1e3 * qr_err,
            "normal equations {ne_err:e} vs QR/SVD {qr_err:e}"
        );
    }

    /// The SVD reconstructs its input, has orthonormal factors, and recovers
    /// known singular values in descending order
    ///
    /// # Panics
    /// * If U*S*V^T differs from R, U or V isn't orthonormal, or a singular
    ///   value is wrong or out of order
    #[test]
    fn test_svd_factors() {
        let mut rng = 4;
        let sigma_true = [5.0, 2.0, 1.0, 0.5, 1e-2, 1e-4, 1e-6];
        let a = matrix_with_singular_values(&mut rng, 30, &sigma_true);
        let qr = householder_qr(&mut a.clone(), &mut alloc::vec![0.0; 30]);
        let svd = jacobi_svd(&qr.r);

        for i in 0..7 {
            for j in 0..7 {
                let usv: f64 = (0..7)
                    .map(|k| svd.u[i][k] * svd.sigma[k] * svd.v[j][k])
                    .sum();
                assert!(abs(usv - qr.r[i][j]) < 1e-14 * 5.0, "U*S*V^T[{i}][{j}]");
                let utu: f64 = (0..7).map(|k| svd.u[k][i] * svd.u[k][j]).sum();
                let vtv: f64 = (0..7).map(|k| svd.v[k][i] * svd.v[k][j]).sum();
                let id = if i == j { 1.0 } else { 0.0 };
                assert!(abs(utu - id) < 1e-12, "U^T*U[{i}][{j}]");
                assert!(abs(vtv - id) < 1e-12, "V^T*V[{i}][{j}]");
            }
        }
        for k in 0..7 {
            let rel = abs(svd.sigma[k] - sigma_true[k]) / sigma_true[k];
            assert!(
                rel < 1e-8,
                "sigma {k} = {} vs {}",
                svd.sigma[k],
                sigma_true[k]
            );
            if k > 0 {
                assert!(svd.sigma[k] <= svd.sigma[k - 1]);
            }
        }
    }

    /// A zero column gives a zero singular value whose right singular vector
    /// names that column, and the pseudo-inverse step stays finite
    ///
    /// # Panics
    /// * If the smallest singular value isn't zero, its vector doesn't point
    ///   at the zero column, or the step isn't finite
    #[test]
    fn test_svd_rank_deficient() {
        let mut rng = 5;
        let mut a = random_matrix::<7>(&mut rng, 20);
        for row in a.iter_mut() {
            row[3] = 0.0;
        }
        let mut b: Vec<f64> = (0..20).map(|_| uniform(&mut rng)).collect();
        let qr = householder_qr(&mut a, &mut b);
        let svd = jacobi_svd(&qr.r);

        assert_eq!(svd.sigma[6], 0.0);
        assert!(
            abs(abs(svd.v[3][6]) - 1.0) < 1e-12,
            "last right singular vector doesn't name column 3"
        );
        let step = damped_step(&svd, &svd_project(&svd, &qr.qtb), 0.0);
        assert!(step.iter().all(|x| x.is_finite()));
        assert_eq!(step[3], 0.0);
    }

    /// The damped step from one SVD matches a fresh QR of [A; sqrt(mu)*I] for every mu
    ///
    /// # Panics
    /// * If any re-solved step differs from the fresh factorization
    #[test]
    fn test_damped_step_resolve() {
        let mut rng = 6;
        let sigma = [3.0, 1.0, 0.4, 0.1, 1e-2, 1e-3, 1e-4];
        let a = matrix_with_singular_values(&mut rng, 25, &sigma);
        let b: Vec<f64> = (0..25).map(|_| uniform(&mut rng)).collect();
        let qr = householder_qr(&mut a.clone(), &mut b.clone());
        let svd = jacobi_svd(&qr.r);
        let c = svd_project(&svd, &qr.qtb);

        for mu in [0.0, 1e-8, 1e-4, 1e-2, 1.0, 1e3] {
            let step = damped_step(&svd, &c, mu);

            // Fresh factorization of the augmented system [A; sqrt(mu)*I]*delta = [b; 0]
            let mut aug = a.clone();
            let mut baug = b.clone();
            for k in 0..7 {
                let mut row = [0.0; 7];
                row[k] = sqrt(mu);
                aug.push(row);
                baug.push(0.0);
            }
            let fresh_qr = householder_qr(&mut aug, &mut baug);
            let fresh = solve_upper(&fresh_qr.r, &fresh_qr.qtb);
            let rel = max_diff(&step, &fresh) / norm(&fresh);
            assert!(rel < 1e-9, "mu = {mu}: relative difference {rel:e}");
        }
    }

    /// Cholesky reproduces a positive definite matrix and forward substitution inverts it
    ///
    /// # Panics
    /// * If L*L^T differs from A, L isn't lower triangular, or L*x differs from y
    #[test]
    fn test_cholesky_factor_and_solve() {
        let mut rng = 7;
        let g = random_matrix::<6>(&mut rng, 6);
        // A = G^T*G + I is symmetric positive definite
        let mut a = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] = (0..6).map(|k| g[k][i] * g[k][j]).sum::<f64>();
            }
            a[i][i] += 1.0;
        }
        let l = cholesky(&a).expect("positive definite");
        for i in 0..6 {
            for j in 0..6 {
                let llt: f64 = (0..6).map(|k| l[i][k] * l[j][k]).sum();
                assert!(abs(llt - a[i][j]) < 1e-13 * 10.0, "L*L^T[{i}][{j}]");
            }
            for j in (i + 1)..6 {
                assert_eq!(l[i][j], 0.0);
            }
        }
        let y: [f64; 6] = core::array::from_fn(|_| uniform(&mut rng));
        let x = solve_lower(&l, &y);
        let lx: [f64; 6] = core::array::from_fn(|i| (0..6).map(|k| l[i][k] * x[k]).sum());
        assert!(max_diff(&lx, &y) < 1e-14);
    }

    /// Cholesky rejects matrices that aren't positive definite
    ///
    /// # Panics
    /// * If an indefinite, singular, zero, negative or non-finite matrix is accepted
    #[test]
    fn test_cholesky_rejects() {
        let indefinite = [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let singular = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let zero = [[0.0; 3]; 3];
        let negative = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let nan = [[1.0, 0.0, 0.0], [f64::NAN, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inf = [[f64::INFINITY, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for (name, m) in [
            ("indefinite", indefinite),
            ("singular", singular),
            ("zero", zero),
            ("negative", negative),
            ("nan", nan),
            ("inf", inf),
        ] {
            assert!(cholesky(&m).is_none(), "{name} accepted");
        }
    }
}
