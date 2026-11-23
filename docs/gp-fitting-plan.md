# GP Fitting Plan (v0.4.0)

Working design notes for fitting a General Perturbations (GP) element set to timestamped state and position data. This is a planning document, not part of the math spec. The final 0.4.0 design moves into `docs/mathspec.md` as Section 9 after the fit passes validation and before 0.4.0 is published.

## Goal
Given timestamped TEME states or positions from GPS, a high-fidelity propagator, another ephemeris, or a sensor such as radar, find the GP element set (TLE/OMM) whose SGP4 propagation best matches them in a weighted least-squares sense. This is the reverse of propagation.

The method builds on Vallado & Crawford, *SGP4 Orbit Determination*, AIAA 2008-6770. The B* prior, outlier cutoff, mode-boundary handling and noise handling were designed in this plan and have not been compared with that paper.

## Noise
**Measurement noise and model error are separate.** Each residual has two independent error sources: the sensor's measurement noise and SGP4's model error. Model error is a property of SGP4, not of the data, so it is internal to the fit and not a user option. Measurement noise is the only noise a user supplies, optionally and per observation, and it is added to the model-error term rather than replacing it. Sensor sigmas used alone as weights would mis-state relative accuracy: a 10 m GPS state would outweigh a 100 m radar position 100 to 1, though both residuals are about 1 km of model error, and that source would then dominate the fit or, as a minority, be rejected wholesale by the outlier test.

**What this means in practice.** SGP4's model error, about 1 km in LEO, dominates every residual from a good source (GPS, a high-fidelity propagator), so those observations are weighted about equally whether or not their measurement noise is given, which is close to optimal for them. Measurement noise matters when it is comparable to or larger than the model error (radar, optical angles, poor receivers), and then it down-weights those observations in the directions where they are poor.

## Implementation steps
Each step is one reviewable change with its tests, built on the steps before it. The fit types are public from step 3; `fit_gp` is public from step 7, where it first fits end to end. Every step's "done when" also includes `cargo build --no-default-features` succeeding (see the no_std decision).

| Step | Builds | Done when |
|---|---|---|
| 1. Linear algebra | Householder QR, 7×7 Jacobi SVD and Cholesky, in the crate with no new dependency (Solve accuracy). | The QR/SVD accuracy and re-solve tests pass; Cholesky rejects matrices that aren't positive definite. |
| 2. Elements and time | Equinoctial ↔ classical elements with f_r (equations 1 and 9), state → osculating equinoctial, τ from the split Julian date, epoch default and rounding. | The element round-trip and exact-zero rows of the Output tests pass. |
| 3. API and validation | Every public type, with the observation constructor and `Default`, and algorithm step 1. | The Options tests, the `InvalidObservation` checks and the `TooFewObservations` edge cases pass. |
| 4. Noise model | Rⱼ, Lⱼ, whitening and ρⱼ for states and positions, for a given n₀ (equation 3). | Equivalent noise forms give the same Rⱼ and Lⱼ, and the outlier-normalization test passes on ρⱼ alone. |
| 5. Starting guess | Fixed-point loop (single and two-loop), nearest-state candidates, Gibbs/Herrick–Gibbs triples, n₀, scoring and the mode record (algorithm steps 2–3). | The candidate-selection parts of the starting-guess and starting-triple tests and the n₀ part of the velocity-term check pass; the Orekit cases' guesses meet the Orekit fixed-point tolerance. |
| 6. LM solve | Residuals, the finite-difference Jacobian with its one-sided rules, the B* prior row, column scaling, the LM step, convergence and held flags (algorithm steps 4–5, equations 4–7). | The finite-difference, LM-convergence and held-flag tests pass, and noise-free recovery holds on the internal elements. |
| 7. Finish | Covariance and rank test, output elements, metadata and rev number, residuals and RMS from the returned `Sgp4`, `NotConverged` (algorithm step 7, equation 8). `fit_gp` now fits end to end, without outlier rejection. | Every Orekit criterion and the Output, edge-case, deterministic, Monte Carlo, positions-and-noise, starting-guess and API-stability tests pass, except the parts that need outlier rejection. |
| 8. Outlier rejection | Algorithm step 6. | The Outliers tests and the rejection parts of earlier tests pass. |
| 9. Real data | The Swarm, GPS and QZSS fixtures; final values for everything set by testing. | The real-data criteria pass and the Values set by testing table is filled in. |
| 10. Docs and math spec | Every "the docs must say" item, the fit-window guidance (open item 13), the design moved into `docs/mathspec.md` Section 9, and the Nielsen citation checked. | The docs and math spec are reviewed. |

---

# 0.4.0 design

## Decisions
| Topic | Decision |
|---|---|
| Observations | Each `FitObservation` is a time plus a `FitMeasurement`: a full state (position and velocity) or a position only, mixed freely in one list. A set needs at least one full state or 3 distinct position times; with the prior row that always gives at least p = 7 rows. `FitMeasurement` is `#[non_exhaustive]`, so other kinds (angles, range) can be added later without a breaking change. |
| Frame | Each measurement states its frame as a `CoordinateFrame`, so other frames can be added later without changing the API. `CoordinateFrame` currently has only TEME, so there is no frame check and no frame error. The fit matches on the frame exhaustively inside the crate, so adding a variant forces a decision about its conversion to TEME at compile time. ECEF input (GPS receivers) is deferred. |
| Time scale | Each time carries its scale in `DateTime.timezone`. `Timezone` currently has only UTC, so there is no time-scale check and no time-scale error; matching exhaustively forces a conversion decision when a scale is added. Every time must pass `validate_datetime`, or the fit returns `InvalidObservation` or `InvalidOptions { field }`. The output epoch is always UTC, because TLE and OMM epochs are. τⱼ is computed from the split Julian date that `utc2jday` returns, so time differences keep microsecond precision. Leap seconds are not handled: an arc that spans one is 1 s off after it, as SGP4 itself is. The docs must say so. |
| Input order | Unsorted input and duplicate timestamps are accepted. `residuals[j]` always matches input j. |
| Noise model | Rⱼ = R_sensor,j + σ²_model·Mⱼ, with Mⱼ = diag(I₃, n₀²·I₃) for a state and I₃ for a position. R_sensor,j is the observation's optional measurement noise, given as isotropic sigmas or a full covariance (3×3 for a position, 6×6 for a state, with r–v correlation), and 0 without it. Adding the two, rather than taking a threshold or the larger, keeps weights continuous and works per direction: radar range accurate to metres with kilometres of cross-range error exceeds the model error in some directions and not others. A sigma of 0 or a singular covariance is accepted, because σ²_model·Mⱼ keeps every Rⱼ positive definite, and a covariance need only be symmetric to roundoff, so one rotated from RIC sigmas into TEME passes. Each Rⱼ is built once, when n₀ is set, and stays fixed through every outlier pass. Each residual is whitened by Lⱼ⁻¹, where Rⱼ = Lⱼ·Lⱼᵀ. |
| Model error | σ_model = 1 km, an internal constant with no user option. The velocity term is n₀·σ_model, with n₀ the median mean motion of the starting candidates, fixed for the fit: a position error σ on an orbit goes with about n·σ of velocity error, so the two stay consistent in every regime (about 1.1 m/s in LEO, 0.07 m/s at GEO, where a fixed 1 m/s would give velocities almost no weight). Where measurement noise is small beside it, σ_model only sets the absolute level: the balance against the B* prior, the covariance of a single state and the units of the convergence floor. Other covariances are rescaled by σ̂₀² and outlier rejection is scale-free, so fit quality shows in σ̂₀² and the RMS values. Where measurement noise is comparable, σ_model also sets how those observations balance against model-dominated ones. The term is isotropic, though SGP4's real error is mostly along-track (Known limitations). σ_model is fixed rather than estimated: on nearly noise-free data an estimate goes to 0, the data weight grows without bound and B* escapes its prior to absorb roundoff and unmodeled forces. An estimate would help only when measurement-dominated and model-dominated observations are mixed and the true model error is far from 1 km; it can be added later, with a floor. |
| Epoch | Defaults to the latest input observation time, even if that observation is later rejected, so the epoch can't move during a fit. A user override may fall anywhere, including outside the observation span, with no warning. Either way, the epoch is rounded to the nearest 1e-8 day, the TLE epoch resolution (about 0.86 ms), before anything else runs; without this, writing the TLE would shift the orbit by up to about 3 m along-track in LEO. The rounding is always applied, with no option: it isn't a format check (the writers still do those), it costs OMM users at most 0.43 ms of epoch placement and nothing in accuracy, and it makes every result TLE-safe by default. |
| Unknowns | p = 7: the six orbital elements in equinoctial form plus B*. `n` is Kozai, the value `init_sgp4` takes and the TLE/OMM publishes. B* is always fitted; there is no option to hold it fixed (a tight prior does that). |
| B* prior | A weak prior row, (B* − B*₀)/σ_B, is always in the fit. B*₀ is the starting B* and stays fixed for the whole fit. It comes from `bstar_prior_center` if set, otherwise from `initial_guess`, otherwise 0, so a user can supply a B* without also forcing the starting orbit. σ_B is `bstar_prior_width`, default 1e-2. It is a soft penalty, not a bound: "B* ≈ center ± width", where a tight width effectively holds B* fixed and a loose one lets the data decide where it can. It matters most for a single state at the epoch and for short arcs, where the data says little or nothing about B*: without a center, a LEO state → TLE comes out with B* = 0 and no drag. The docs must give the ballistic-coefficient estimate B* = ρ₀·B/2, with B = C_d·A/m [m²/kg] and ρ₀ = 2.461e-5 kg/m²/ER (Vallado). |
| Initial guess and metadata | Two separate options. `initial_guess` is an optional element set that gives the starting orbit, ndot/nddot, and B*₀ unless `bstar_prior_center` is set. `metadata` sets any of name, catalog number, designator, classification, element set number and rev number. Each metadata field comes from `metadata` if set, otherwise from `initial_guess`, otherwise the default: empty name and designator, catalog number 0, classification `U`, element set and rev numbers 0. The element set number is never incremented. `fit_gp` doesn't check metadata against any output format; the TLE/OMM writers do, as for any other element set. `ephemeris_type` is always 0. |
| Rev number | When taken from `initial_guess`, advanced by the whole revolutions between its epoch and the fit epoch, as Orekit does. When set in `metadata`, used as given, as the value at the fit epoch. |
| ndot/2 and nddot/6 | Copied from `initial_guess` if given, otherwise 0. SGP4 doesn't use them, so they never affect the fit. |
| Starting orbit | If `initial_guess` is given, it always supplies the starting orbit, even when the data could, because an option the user sets is honored. A stale guess (e.g. from before a maneuver) isn't overruled; it shows up as a large `rms_position` and σ̂₀², or as `NotConverged`. Without `initial_guess`, the full states nearest the epoch give up to 3 candidates, topped up to 3 with independent Gibbs/Herrick–Gibbs triples of positions, and they are scored by median ρ over all observations. No observation is used by more than one candidate, so one glitched observation can spoil at most one candidate. |
| Single state | There is one entry point, `fit_gp`. A single full state goes through it like any other set: 6 state rows plus the prior row make a square 7×7 system, and LM becomes Newton's method on SGP4(x, τ) = z. When the state is at the epoch and the fixed-point guess has converged, LM returns that guess after one Jacobian build, which the convergence test needs. At τ = 0, B* has no effect on the state, so B* stays at B*₀; set `bstar_prior_center` to give the result a realistic B*. |
| Mode boundaries | Only the three flags stored in the GP set matter: deep space (period 225 min), whole-day resonance and half-day resonance (n window and e = 0.5). Each fit is held to the flags of its starting guess, so its elements always behave the same in ordinary SGP4: a trial step that would change a flag is rejected, and finite differences beside a gate are one-sided. A fit whose answer lies across an edge stalls against it and returns `NotConverged`; the docs must say so. `init_sgp4` is unchanged. There is no second fit on the other side of an edge (Known limitations). |
| Cost function | Weighted least squares with the B* prior row, solved by Levenberg–Marquardt. |
| Jacobian | Central differences, recomputed each iteration with steps h_k = c·s_k from the current x_k (equation 5). c ≈ 6e-6 for n, a_f, a_g, χ, ψ, λ; c_B ≈ 1e-3 for B*, because position depends on B* through a low-order polynomial, so a large step costs almost no accuracy and keeps the column clear of roundoff. Both are tuned by a step-halving test. Next to a gate, or where one side makes SGP4 fail, the 3-point one-sided formula keeps O(h²) accuracy. |
| Solve accuracy | Each Jacobian column is scaled to unit norm. The scaled Ã is factored by Householder QR and the 7×7 R by SVD; ÃᵀÃ is never formed, because that squares the condition number and loses about twice the digits on long arcs, GEO and near-circular orbits. Every LM step, the Gauss–Newton test step and the covariance come from the singular values, so re-solving after a rejected step costs almost nothing. QR and a 7×7 Jacobi SVD are written in the crate, about 150 lines, with no new dependency. |
| LM loop | μ starts at `initial_damping` (default 1e-3, nearly Gauss–Newton) at the start of every LM run, including each pass of the outlier loop. It moves by Nielsen's gain-ratio rule. Converged when the undamped Gauss–Newton step is smaller than `tolerance` × its own 1σ, measured with the scaled covariance and floored for noise-free fits. μ above `max_damping` (default 1e10) means the fit is stuck: `NotConverged`. `max_iterations` counts Jacobian builds per LM run; the ceiling bounds the rejected steps in between. |
| Underdetermined | Tested on the singular values of the scaled Ã: σ_min/σ_max < κ_tol returns `FitError::Underdetermined` with no solution. κ_tol is about 100× the Jacobian error the step-halving test measures (about 1e-4), so a direction that differs from the others only by finite-difference noise is never reported with a finite variance. The test doesn't depend on column order, and the smallest singular vector names the unobservable combination. |
| Outliers | After a fit converges, reject observation j if ρⱼ > c·median(ρ), the median over all input observations, and its residual is above numerical noise; then refit. c is `outlier_cutoff`, default 10 (about 9σ), minimum 3; `None` turns rejection off. The scale comes from the fit's residuals, not the noise model, so the test is scale-free; because measurement noise is added to the model-error term and never replaces it, no source can be swamped by another's stated noise. At most a fraction `max_rejected_fraction` of input observations is rejected, default 0.25, below the median's 50% breakdown point. There is no separate minimum count: the user supplies the data, `None` turns rejection off, and the cap limits small sets. A pass that would exceed the cap, or leave no full state and fewer than 3 distinct position times, keeps the previous rejected set, and a refit that fails to converge is discarded in favour of the previous pass's fit, with no error. No Huber weighting. |
| Covariance | Exactly one covariance, always valid: σ̂₀²·P, the standard least-squares covariance scaled by how well the fit actually matched, as Orekit returns; for a single state, where σ̂₀² is `None`, P from the noise model. There is no optional or NaN covariance. It is a diagnostic: a TLE carries no covariance, and no covariance from the fit describes prediction beyond the arc. The docs must say it is the formal covariance, which assumes independent errors; that it is optimistic for dense data, because SGP4's model error is correlated in time; and that the honest measures of a TLE's accuracy are its residual RMS and its prediction error against held-out data. There is no robust (autocorrelation-consistent) estimate. |
| Output | Residuals and RMS values come from the returned `Sgp4`, re-initialized from the final GP set, so they describe exactly what the user propagates. Residuals are observed − computed, with RIC axes from the fitted state. RMS values are unweighted 3D norms over kept observations; `rms_velocity` and the velocity residuals are `None` where there is no velocity. `FitResult` includes `retrograde_factor`, without which the covariance can't be interpreted. `rejected` is ascending. The docs must say that `velocity_ric` is a rotated TEME residual, not the rate of change of the RIC position error. |
| no_std | The crate is `#![no_std]` with `alloc`, and the fit stays within that: `Vec`, `String` and `Box` come from `alloc`, and all math (QR, SVD, Cholesky, Gibbs/Herrick–Gibbs, element conversions) uses `libm` (`sqrt`, `atan2`, …), never the `f64` methods that need `std`. No new dependency. Every fit type and function builds with `--no-default-features`. |
| API stability | Every public struct and enum in the fit API is `#[non_exhaustive]`, so later releases add fields, variants and types without a breaking change. The `FitMeasurement` variants carry their fields inline and are built as struct literals, so a new field on an existing variant would be breaking; new noise forms, frames and measurement kinds go into their own non-exhaustive enums instead. Fields that position-only observations leave empty (the velocity residual, `rms_velocity`) are `Option`. Callers build options with `FitOptions::default()` and assign fields, and build observations with `new_fit_observation`, because struct literals don't work on non-exhaustive structs outside the crate. 0.4.0 is the cheap moment: adding the attribute later would itself be breaking. The `Fit` prefix is kept because everything is re-exported at the crate root. |

## API
Field names and types are provisional.

```rust
#[non_exhaustive]
pub struct FitObservation {
    pub datetime: DateTime,               // UTC is the only time scale today
    pub measurement: FitMeasurement,
}

#[non_exhaustive]                         // angles, range, ... later as new variants
pub enum FitMeasurement {                 // variants are built as struct literals
    State {
        frame: CoordinateFrame,           // TEME is the only frame today
        position: [f64; 3],               // [km]
        velocity: [f64; 3],               // [km/s]
        noise: Option<FitStateNoise>,     // measurement noise; None: model term only
    },
    Position {
        frame: CoordinateFrame,
        position: [f64; 3],               // [km]
        noise: Option<FitPositionNoise>,  // measurement noise; None: model term only
    },
}

#[non_exhaustive]                         // e.g. RIC sigmas can be added later
pub enum FitStateNoise {
    Sigmas { position: f64, velocity: f64 },  // [km], [km/s]; each ≥ 0, 0 meaning no noise
    Covariance([[f64; 6]; 6]),            // [km², km²/s, km²/s²], includes r–v correlation
}

#[non_exhaustive]
pub enum FitPositionNoise {
    Sigma(f64),                           // isotropic [km], ≥ 0
    Covariance([[f64; 3]; 3]),            // [km²]
}

// The one constructor, for every kind of measurement (a free function, as
// everywhere in the crate):
//
//     let obs = new_fit_observation(t, FitMeasurement::Position {
//         frame: CoordinateFrame::TEME, position: r, noise: None,
//     });
pub fn new_fit_observation(datetime: DateTime, measurement: FitMeasurement) -> FitObservation;

// GpMetadata and FitOptions implement Default. Being #[non_exhaustive], they are
// built from the default with fields assigned:
//
//     let mut opts = FitOptions::default();
//     opts.bstar_prior_center = Some(2.5e-4);
//     opts.metadata.catalog_number = Some(25544);

#[non_exhaustive]
pub struct GpMetadata {               // default: all None
    pub name: Option<String>,
    pub catalog_number: Option<i32>,
    pub intl_designator: Option<String>,
    pub classification: Option<char>,
    pub element_set_number: Option<i32>,
    pub rev_number: Option<i64>,      // at the fit epoch
}

#[non_exhaustive]
pub struct FitOptions {
    pub epoch: Option<DateTime>,             // default: latest observation
    pub initial_guess: Option<GenPerturbElementSet>, // default None
    pub metadata: GpMetadata,                // default: all None
    pub bstar_prior_center: Option<f64>,     // B*₀ [1/ER]; None: initial_guess's B*, else 0
    pub bstar_prior_width: f64,              // 1σ of the B* prior [1/ER], default 1e-2
    pub outlier_cutoff: Option<f64>,         // multiple of median ρ, default Some(10.0); None: off
    pub max_rejected_fraction: f64,          // of input observations, default 0.25, in (0, 0.5)
    pub max_iterations: usize,               // Jacobian builds per LM run, default set by testing
    pub tolerance: f64,                      // step size in units of 1σ, default about 1e-3
    pub initial_damping: f64,                // starting μ (unit-norm columns), default 1e-3
    pub max_damping: f64,                    // μ ceiling, default 1e10
    pub wgs: Option<Wgs>,                    // default WGS72
}

#[non_exhaustive]
pub struct FitResult {
    pub sgp4: Sgp4,                    // ready to propagate or write out
    pub rms_position: f64,             // √(Σ‖Δr‖²/N_kept), unweighted [km]
    pub rms_velocity: Option<f64>,     // same for velocity [km/s]; None if no kept observation has velocity
    pub residuals: Vec<FitResidual>,   // one per input observation, including rejected ones
    pub rejected: Vec<usize>,          // input indices of rejected observations, ascending
    pub covariance: [[f64; 7]; 7],     // σ̂₀²·P (P for a single state); always valid; order, units:
                                       // n (Kozai, rad/min), a_f, a_g, χ, ψ, λ (rad), B* (1/ER)
    pub retrograde_factor: f64,        // f_r, +1 or −1; defines a_f, a_g, χ, ψ in the covariance
    pub variance_factor: Option<f64>,  // σ̂₀², None with no degrees of freedom (a single state)
    pub iterations: usize,             // Jacobian builds, summed over all outlier passes
}

#[non_exhaustive]
pub struct FitResidual {                // observed − computed, from the returned Sgp4
    pub position: [f64; 3],             // TEME [km]
    pub velocity: Option<[f64; 3]>,     // TEME [km/s]; None for a position-only observation
    pub position_ric: [f64; 3],         // radial, in-track, cross-track [km], axes from the fitted state
    pub velocity_ric: Option<[f64; 3]>, // TEME velocity residual rotated into RIC axes [km/s]; None likewise
}

#[non_exhaustive]
pub enum FitError {
    /// Observation at this index has a non-finite value, a time that is not a valid date,
    /// a negative noise sigma, or a noise covariance that is not symmetric or not positive semidefinite
    InvalidObservation { index: usize },
    /// A FitOptions value is out of range, or the epoch is not a valid date.
    /// `field` names the option, e.g. "epoch" or "bstar_prior_width".
    InvalidOptions { field: &'static str },
    /// No full state and fewer than 3 distinct position times. Outlier rejection never causes it.
    TooFewObservations,
    /// No starting guess could be built: without an initial GP set, no candidate could be
    /// built from the data (e.g. positions with no pass of 3 points inside the span cap)
    /// or SGP4 failed on every one; with one, the fixed-point loop at the
    /// epoch couldn't evaluate any iterate
    StartingGuessFailed,
    /// SGP4 rejected the initial GP set or the epoch
    Sgp4(Sgp4Error),
    /// The covariance can't be computed, so no solution is returned
    Underdetermined,
    /// Iteration cap, μ ceiling, or a Jacobian column SGP4 couldn't evaluate on either side.
    /// Holds the best solution, with its covariance; never used as a fit result.
    NotConverged(Box<FitResult>),
}

pub fn fit_gp(obs: &[FitObservation], opts: &FitOptions) -> Result<FitResult, FitError>;
```

## Algorithm
1. **Validate input.** Checks run in this order and stop at the first failure.
   1. Options: `bstar_prior_width` finite and positive; `bstar_prior_center`, if set, finite (no range check: negative B* is legal, and format limits belong to the writers); `outlier_cutoff`, if set, finite and ≥ 3; `max_rejected_fraction` finite, > 0 and < 0.5; `max_iterations` ≥ 1; `tolerance` finite and positive; `initial_damping` and `max_damping` finite and positive, with `initial_damping` < `max_damping` (a violation names `max_damping`); `epoch`, if set, a valid date. Otherwise `InvalidOptions { field }`, naming the first field that failed. Metadata is not checked.
   2. Each observation: time a valid date; position and any velocity finite; any measurement-noise sigma finite and ≥ 0; any measurement-noise covariance finite, symmetric to a relative tolerance of about 1e-12 (then symmetrized), and positive semidefinite. Each 3-vector block (position; velocity) is judged in its own units: component i gets the scale dᵢ = √(mean variance of its block), entry (i, j) may be asymmetric by at most 1e-12·dᵢ·dⱼ, and the matrix scaled by 1/(dᵢ·dⱼ) plus ε·I, ε ≈ 1e-12, must have a Cholesky factor. One shift by ε·tr(R_sensor) would hide a negative eigenvalue in a velocity block (km²/s²) beside km² position variances. A block with zero scale must be exactly zero, correlations included; an all-zero matrix counts as no noise. Otherwise `InvalidObservation { index }`. There is no range check on magnitudes; getting units right is the caller's job.
   3. The set: at least one full state or 3 distinct position times, or `TooFewObservations`. With the prior row this always gives at least p = 7 rows.
2. **Starting guess.**
   1. **B*₀** is `bstar_prior_center` if set, otherwise the initial GP set's B*, otherwise 0, for every candidate. It is also the prior center.
   2. **Sources.** With an initial GP set: propagate it to the fit epoch with SGP4 and use that state as the only candidate; if SGP4 can't, return `FitError::Sgp4(err)`. Without one: the full states nearest the epoch, at distinct times, up to 3, then triples (step 2.3) until there are 3 candidates. A triple gives a full state at its middle point's time. If no candidate can be built, return `StartingGuessFailed`.
   3. **Triples from three positions.**
      1. The middle point r₂ needs positions on both sides, since both methods return the velocity at r₂. Add one position on each side, spread as widely as the span cap allows, because wider spread reduces the effect of noise. Triples are tried in order of the middle point's distance from the epoch.
      2. Separations are θ₁₂ = atan2(‖r₁ × r₂‖, r₁ · r₂) and likewise θ₂₃; `atan2` keeps precision at small angles, where `acos` does not.
      3. Four checks: coplanar (r₁ within about 1° of the plane of r₂ and r₃); consistent direction (r₁ × r₂ and r₂ × r₃ point the same way, so a wrap past 180° is caught); a span cap (θ₁₂ + θ₂₃ below about 60–90°); and timing.
      4. **Timing.** The separations are measured from geometry alone, so two positions one orbit plus 10° apart look exactly like positions 10° apart, and Gibbs or Herrick–Gibbs would return a velocity wrong by km/s. Each leg's Δt must be within a factor of about 2 of θ/√(μ/r̄³), where r̄ is the mean radius of its two points. The factor of 2 covers the spread in angular rate up to e ≈ 0.75, and the check rules out any triple whose points come from different passes.
      5. Use Gibbs if both θ₁₂ and θ₂₃ are at least θ_switch (about 3°), otherwise Herrick–Gibbs.
      6. A triple that fails a check is skipped, and the next middle point is tried. Observations already used by an accepted candidate, including the full states, are not available to later triples. If the data doesn't allow 3 independent candidates, fewer are built.
   4. **The fixed-point loop** finds mean elements whose SGP4 state at τ = 0 reproduces a target state, in the same equinoctial elements as LM:
      1. Convert the target state to osculating equinoctial elements e_target = (n, a_f, a_g, χ, ψ, λ), with f_r chosen from its inclination and held fixed.
      2. Start from x₀ = e_target, with B* at B*₀.
      3. Each iteration: initialize SGP4 from x_k, convert its state at τ = 0 to osculating equinoctial e_k, form Δ = e_target − e_k (Δλ wrapped into (−π, π]), and set x_{k+1} = x_k + s·Δ, with s = 1.
      4. **Mismatch.** q = max(‖Δr‖/tol_r, ‖Δv‖/tol_v) over the full state. The loop has converged when q < 1. It is checked before the first step, so a target that x₀ already reproduces exits at once.
      5. **Damping.** s is halved when a step doesn't reduce q (as in Orekit's O5 case) or when SGP4 can't initialize or evaluate the iterate (e.g. e ≥ 1). After a step that reduces q, s ← min(1, 2s), so brief damping doesn't leave the loop crawling: at a fixed s = 0.25, going from 1 km to 1 mm takes about 48 iterations instead of about 3.
      6. **Iteration cap.** Set by testing. Every trial step counts, including rejected ones, so damping can't loop forever. If the cap is reached, the iterate with the smallest q is kept.
      7. **Failure.** If no iterate can be evaluated at all, the loop fails. A failed data candidate is dropped; if none is left, or the loop fails on the initial GP set, return `StartingGuessFailed`.
      - Classical elements are not used, because osculating ω, M and Ω are badly defined near e = 0 and i = 0, which makes the update oscillate (O1, O2, O5).
      - Each iterate uses whatever SGP4 branch its own elements select; nothing is frozen inside the loop.
      - **A source state not at the epoch** (the normal case, because the epoch is rounded to the TLE grid): run the loop at the state's time, propagate the result to the epoch, and run the loop again there. The second loop is the same code path as the initial GP set.
   5. **Noise model.** Set n₀ to the median mean motion (in rad/s) of the candidates at the epoch, the only candidate's with an initial GP set, and build every Rⱼ from it (equation 3). Both stay fixed for the whole fit, so every candidate is scored with the same weights, one glitched candidate among 3 can't move them, and J stays comparable between iterations.
   6. **Scoring.** Propagate each candidate to every observation time and compute ρⱼ (equation 3). If SGP4 fails at an observation time, ρⱼ = ∞ there. The candidate with the lowest median(ρ) starts LM; on an exact tie, the one whose source (a state's time, or a triple's middle point) is nearer the epoch wins, so results are reproducible. The median keeps a few glitched or failed observations from deciding the score. All candidates use the same B*₀. There is no early exit when a candidate scores well: scoring costs about N SGP4 calls per candidate, about a fifth of one LM iteration for all 3. With an initial GP set there is one candidate and no scoring.
3. **Record the mode.** Store the `deep_space`, `whole_day_resonance` and `half_day_resonance` flags and f_r from the winning candidate's final iterate at the fit epoch, not from `initial_guess` at its own epoch, because the elements at the fit epoch can fall on the other side of a boundary.
4. **Iterate** (equations 5–7):
   1. Evaluate the residuals and cost.
   2. Build the Jacobian by central differences, with steps h_k = c·s_k.
   3. Apply the weights and scale the columns.
   4. Factor Ã by QR, then R by SVD, and solve for the damped step.
   5. Accept or reject the step by the gain ratio r, and update μ by Nielsen's rule.
   - Rejecting a step reuses the same Jacobian.
   - A trial step that makes SGP4 fail (e ≥ 1, n ≤ 0, decay) counts as a rejected step, not an abort.
   - A trial step that changes any recorded mode flag counts as a rejected step.
   - If a central difference would cross a gate edge, use the 3-point one-sided difference on the current side.
   - If SGP4 fails at one perturbed point, use the 3-point one-sided difference on the other side. If both sides fail, stop iterating as for the iteration cap: keep the best solution so far, which leads to `NotConverged`. The current point always evaluates, because it was an accepted step.
5. **Converge.** Before each step, compute the undamped Gauss–Newton step from the same SVD and test it (equation 7). If it passes, the run has converged. If σ_min/σ_max < κ_tol at this iterate, skip the test this iteration. The run also stops on `max_iterations` Jacobian builds or μ > `max_damping`; both mark it as not converged and keep the best solution so far. The roundoff case, where steps fail only because J is at noise level, passes the step test before μ climbs, so reaching the ceiling means the fit is really stuck. A fit stalled against a gate edge ends here and returns `NotConverged`.
6. **Reject outliers.** Skip if `outlier_cutoff` is `None` or the first LM run didn't converge: residuals from an unconverged fit can't be trusted to identify outliers. Each pass works from the current converged fit:
   1. Compute ρⱼ for every input observation, including ones rejected earlier, so every pass uses the same population.
   2. Flag j if ρⱼ > c·median(ρ), the median over all input observations, and its residual is above numerical noise: ‖Δr‖ > 1 mm or, for an observation with velocity, ‖Δv‖ > 1 µm/s. The floors stop roundoff being flagged in noise-free fits; the velocity floor lets a state with a good position and a bad velocity be rejected. Both are fixed constants.
   3. If the flagged set would hold more than ⌊`max_rejected_fraction`·N_in⌋ observations, where N_in is the number of input observations, or would leave no full state and fewer than 3 distinct position times, keep the previous rejected set. The fit returns normally; `rejected` and `rms_position` show what happened.
   4. Stop when the rejected set is unchanged.
   5. Otherwise refit (steps 4–5) without the flagged observations, starting from the current solution. Every point is rechecked each pass, so a wrongly flagged point can return. If the refit doesn't converge, discard it, keep the previous pass's fit and rejected set, and stop.
   6. Stop after 10 passes (set by testing). If the limit is hit, the set is oscillating: keep the pass with the lowest median(ρ) over all input observations.
7. **Finish.**
   1. Compute a fresh, undamped Jacobian at the solution.
   2. Compute P from its SVD (equation 8). If σ_min/σ_max < κ_tol, return `FitError::Underdetermined`. κ_tol is set by testing: degenerate cases must fail and the shortest legitimate arcs must pass. Legitimate correlations are far from the threshold: an n–λ correlation of 0.99999 gives σ_min/σ_max of about 3e-3.
   3. The returned covariance is σ̂₀²·P, or P for a single state, where σ̂₀² is `None`.
   4. Convert the equinoctial elements to classical elements (equation 9), in degrees wrapped into [0, 360) and n in rev/day. The special cases are exact zeros only: e = 0 sets ω = 0, and χ = ψ = 0 sets Ω = 0, which covers both degenerate inclinations (i = 0 when f_r = +1, i = 180° when f_r = −1). Small nonzero values need no threshold. For the node this is exact: any Ω gives the same orbit as long as ω + f_r·Ω is kept. For the eccentricity it is not: SGP4 floors e at 1e-6 along ω, so below that the orbit depends on ω, and ω = 0 at e = 0 is a convention that moves the orbit by up to about 2·a·1e-6 (14 m in LEO) against other splits (open item 19).
   5. Copy ndot/2 and nddot/6 from `initial_guess`, or set them to 0. Set `tle_checksum_valid` to `None`.
   6. Re-initialize SGP4 from the final GP set and compute every residual and RMS from it, not from the internal x, because the round trip through classical elements loses a little precision. Residuals are observed − computed; velocity residuals exist only for observations with velocity. RIC axes come from the fitted state at each time: R = r̂, C = (r × v)/‖r × v‖, I = C × R.
   7. Return the result, or `FitError::NotConverged` if the fit didn't converge. `NotConverged` holds the last accepted iterate, which has the lowest J because accepted iterates always lower it, with all fields filled as for a normal result; if its covariance fails, return `Underdetermined` instead.

## Equations

```
 GP FIT: WEIGHTED NONLINEAR LEAST SQUARES
 ═══════════════════════════════════════════════════════════════════════

 1. UNKNOWNS                                        p = 7

        x = [ n   a_f   a_g   χ   ψ   λ   B* ]ᵀ  ∈ ℝᵖ

        a_f = e·cos(ω + f_r·Ω)         χ = tanᶠʳ(i/2)·sin Ω
        a_g = e·sin(ω + f_r·Ω)         ψ = tanᶠʳ(i/2)·cos Ω
        λ   = M + ω + f_r·Ω            f_r = +1 (i ≤ 90°),  −1 (i > 90°)


 2. OBSERVATIONS                      j = 1 … N,    τⱼ = tⱼ − t₀  (min)

        τⱼ = ( (jdⱼ − jd₀) + (fracⱼ − frac₀) )·1440      split JD from utc2jday

        zⱼ = [ rⱼ ; vⱼ ] ∈ ℝ⁶,   mⱼ = 6,   Sⱼ = I₆               full state
        zⱼ = rⱼ ∈ ℝ³,          mⱼ = 3,   Sⱼ = [ I₃  0₃ ]        position only

            ⎡ z₁ ⎤             ⎡ S₁·SGP4(x, τ₁)  ⎤
        z = ⎢ ⋮  ⎥ ,    h(x) = ⎢        ⋮        ⎥   ∈ ℝᵐ,   m = Σⱼ mⱼ
            ⎣ z_N⎦             ⎣ S_N·SGP4(x, τ_N)⎦

        y(x) = z − h(x)                                    residuals


 3. WEIGHTS AND OUTLIERS

        Rⱼ = R_sensor,j + σ²_model · Mⱼ     R_sensor,j: the observation's measurement
                                             noise (0 if none)

        Mⱼ = diag( I₃ , n₀²·I₃ ),  or I₃ for a position

        n₀ = median mean motion of the starting candidates [rad/s],
             fixed for the fit (not the unknown n)

        σ_model = 1 km                       internal constant, SGP4 model error

        Rⱼ = Lⱼ·Lⱼᵀ   (Cholesky)        W^½ = blockdiag( L₁⁻¹ … L_N⁻¹ )

        Rejected observations are removed from z, h and R, so N and m
        count only kept observations.

        ρⱼ = ‖ Lⱼ⁻¹·yⱼ ‖ / √mⱼ          normalized residual

        reject j  if  ρⱼ > c · median(ρ)    median over all input observations
                 and  ( ‖Δrⱼ‖ > 1 mm  or  ‖Δvⱼ‖ > 1 µm/s  (with velocity) )

        at most ⌊ max_rejected_fraction · N_in ⌋ rejected,  N_in = input observations,
        and the kept set keeps a full state or 3 distinct position times


 4. NONLINEAR PROBLEM

        min  J(x) = y(x)ᵀ · W · y(x)  +  ( (B* − B*₀) / σ_B )²
         x


 5. LINEARIZED SUBPROBLEM  (each iteration, about x_k)

        h(x_k + δ) ≈ h(x_k) + H·δ          H = ∂h/∂x |ₓₖ   ∈ ℝ^(m×p)

                   h(x_k + h_k·e_k) − h(x_k − h_k·e_k)
        H₍:,k₎  =  ─────────────────────────────────────    central diff
                                  2·h_k

                   −3·h(x_k) + 4·h(x_k ± h_k·e_k) − h(x_k ± 2h_k·e_k)
        H₍:,k₎  =  ± ───────────────────────────────────────────────    one-sided
                                         2·h_k

        h_k = c·s_k       recomputed from the current x_k each iteration

            element            s_k                   c
            n                  |n|                   ≈ 6e-6  (ε^⅓, tuned)
            a_f a_g χ ψ λ      1                     ≈ 6e-6
            B*                 max( |B*|, σ_B )      c_B ≈ 1e-3  (tuned)

            ⎡ W^½ · H    ⎤          ⎡ W^½ · y(x_k)       ⎤
        A = ⎢            ⎥      b = ⎢                    ⎥
            ⎣ e_Bᵀ / σ_B ⎦          ⎣ (B*₀ − B*_k) / σ_B ⎦

        e_B = unit vector selecting B*

        min  ‖ b − A·δ ‖²
         δ


 6. COLUMN SCALING

        D = diag( ‖A₍:,1₎‖ … ‖A₍:,p₎‖ )

        Ã = A·D⁻¹          δ̃ = D·δ          ⇒   diag(ÃᵀÃ) = I


 7. LEVENBERG–MARQUARDT STEP

        Ã = Q·R            (Householder QR, m × 7 → 7 × 7)
        R = U·Σ·Vᵀ         (Jacobi SVD, σ₁ ≥ … ≥ σ₇)        ÃᵀÃ is never formed

        c = Uᵀ·Qᵀ·b        g̃ = Ãᵀb = V·Σ·c

        δ̃(μ) = Σᵢ  σᵢ / (σᵢ² + μ) · cᵢ · vᵢ      solves ( ÃᵀÃ + μ·I )·δ̃ = g̃

        δ = D⁻¹·δ̃          x_(k+1) = x_k + δ

        gain ratio      r = ( J − J' ) / δ̃ᵀ·( μ·δ̃ + g̃ )

        r > 0  →  accept,  μ ← μ·max( ⅓, 1 − (2r − 1)³ ),  ν ← 2
        r ≤ 0  →  reject,  μ ← μ·ν,  ν ← 2ν,  re-solve from the same SVD

        μ₀ = initial_damping,  ν₀ = 2,  reset at the start of every LM run
        μ > max_damping  →  stop, not converged

        convergence  (undamped step δ̃_GN = δ̃(0) = Σᵢ cᵢ/σᵢ · vᵢ)

            δ̃_GNᵀ · g̃  ≤  tolerance² · s²

            s² = max( J / max(m + 1 − p, 1),  s²_floor )

        δ̃_GNᵀ·g̃ = δ_GNᵀ·P⁻¹·δ_GN, the step's squared Mahalanobis length.
        s² ≈ σ̂₀², so the test uses the scaled covariance. s²_floor is about
        the σ̂₀² that 1 µm RMS residuals give, set by testing, because
        noise-free fits drive J to roundoff.


 8. COVARIANCE & FIT QUALITY       (at x̂,  μ = 0,  fresh Jacobian)

        P = ( HᵀWH + e_B·e_Bᵀ/σ_B² )⁻¹  =  D⁻¹ · V · Σ⁻² · Vᵀ · D⁻¹

        Underdetermined  if  σ₇ / σ₁ < κ_tol        κ_tol ≈ 100 × Jacobian error ≈ 1e-4

                    J(x̂)
        σ̂₀²  =  ───────────                  ≈ 1 if the noise model is realistic
                 m + 1 − p

        σ̂₀² is None when m + 1 − p = 0, which only a single full state gives

        returned covariance:  σ̂₀²·P,  or  P  when σ̂₀² is None

        σ_k  = √P_kk          corr_ij = P_ij / √( P_ii · P_jj )


 9. OUTPUT ELEMENTS                (at x̂, f_r fixed for the fit)

        e = √( a_f² + a_g² )                ω + f_r·Ω = atan2( a_g, a_f )

        Ω = atan2( χ, ψ )                   i = 2·atan( √(χ² + ψ²) ),  or π − that if f_r = −1

        M = λ − ω − f_r·Ω                   n [rev/day] = n [rad/min] · 1440 / 2π

        e = 0  →  ω = 0          χ = ψ = 0  →  Ω = 0  (i = 0, or i = 180° if f_r = −1)

        angles in degrees, wrapped into [0, 360)
```

Because the columns of Ã have unit norm, Marquardt's diag(ÃᵀÃ) damping reduces to μ·I: column scaling and Marquardt scaling are the same operation.

The prior row keeps the B* column norm at or above 1/σ_B, so column scaling can't divide by zero on that column, and a B* column that is exactly zero (GEO over a few hours) is harmless: P then reports σ_B for B*. The prior counts in J and the degrees of freedom, and it is excluded from `rms_position` and `rms_velocity`.

The median is used because outliers inflate an RMS scale: one point can never exceed √N times the RMS, so a 10× RMS cutoff can't fire below N = 100. For clean Gaussian residuals, median(ρ) is about 0.94 for full states and 0.89 for positions, so c = 10 means about 9σ. With residuals at about σ_model = 1 km, as in LEO, that is a position glitch above about 10·0.94·√6 ≈ 23 km in a full state, or 10·0.89·√3 ≈ 15 km in a position; the threshold scales with the actual residuals, not with σ_model; smaller glitches stay in at full weight. Rejected observations stay in `residuals` and are excluded from J, the RMS values and the degrees of freedom.

## Validation

### Orekit cross-check cases
Orekit's TLE generation tests (`FixedPointTleGenerationAlgorithmTest` and `LeastSquaresTleGenerationAlgorithmTest`, Apache 2.0) are single-state round trips: parse a TLE, propagate it to its epoch, convert that state back to a TLE, and compare. Only the TLEs are reused; no Orekit code is copied. Each case runs through `fit_gp` on the TEME state at the TLE epoch, a square 7×7 system with the prior row, without `initial_guess`, because passing the TLE would start LM at the answer and test nothing. The metadata comes from `metadata`.
- **Elements:** the orbital elements are recovered from the state alone. The fit also runs once without `bstar_prior_center`, where B* comes out as 0, because at τ = 0 the state doesn't depend on it.
- **Text round trip:** the test passes the TLE's B* as `bstar_prior_center`, so the result's B* equals it exactly, and copies ndot and nddot onto the result before writing (legitimate, because SGP4 doesn't use them).

| Case | Object | Why it is useful |
|---|---|---|
| O1 | 31135 (LEO, i = 2.47°, e = 0.0021) | Near-equatorial, low drag |
| O2 | 27508 (GEO, i = 0.057°, e = 0.0005) | Deep space, near-equatorial, near-circular |
| O3 | 25544 (ISS) | Standard LEO; Orekit reproduces both TLE lines exactly |
| O4 | 05709 (n ≈ 0.996 rev/day, i = 10.8°) | Deep space, 24-hour resonance; Orekit reproduces line 2 exactly |
| O5 | 33153 (INTELSAT 25, i = 0.0042°, e = 0.0003) | Hard convergence case: Orekit's default fixed-point settings fail, and it needs damping (scale 0.5) and more iterations |

```
O1  1 31135U 07013A   11003.00000000  .00000816  00000+0  47577-4 0    11
    2 31135   2.4656 183.9084 0021119 236.4164  60.4567 15.10546832    15

O2  1 27508U 02040A   12021.25695307 -.00000113  00000-0  10000-3 0  7326
    2 27508   0.0571 356.7800 0005033 344.4621 218.7816  1.00271798 34501

O3  1 25544U 98067A   21035.14486477  .00001026  00000-0  26816-4 0  9998
    2 25544  51.6455 280.7636 0002243 335.6496 186.1723 15.48938788267977

O4  1 05709U 71116A   21105.62692147  .00000088  00000-0  00000-0 0  9999
    2 05709  10.8207 310.3659 0014139  71.9531 277.0561  0.99618926100056

O5  1 33153U 08034A   21327.46310733 -.00000207  00000+0  00000+0 0  9990
    2 33153   0.0042  20.7353 0003042 213.9370 323.2156  1.00270917 48929
```

Orekit's tolerances, as targets:

| Case | Orekit fixed point (relative, per element) | Orekit least squares (relative, per element) | Orekit least-squares RMS |
|---|---|---|---|
| O1 | 5.2e-9 | 1.0e-12 | 3.755e-9 |
| O2 | 9.2e-8 | 1.6e-12 | 8.8e-11 |
| O3, O4 | TLE lines reproduced exactly | n/a | n/a |
| O5 | Epoch and metadata reproduced; converges only with damping | n/a | n/a |

**Passing criteria:**
- **Elements:** n, e, i, ω, Ω, M recovered within the Orekit fixed-point tolerance for O1 and O2; the least-squares columns are the stretch target.
- **Text round trip:** O3 and O4 reproduce both TLE lines character for character after `to_tle_*`.
- **Hard case:** O5 converges with the default `FitOptions`. Needing user tuning, as Orekit does, counts as a failure.
- **Metadata:** catalog number, classification, international designator and element set number carry through `initial_guess` exactly, and a field set in `metadata` overrides it.
- **Rev number:** moving the state forward or back by one orbital period before conversion moves the rev number by one, as in Orekit's test.
- **Stale guess:** an `initial_guess` from before a simulated maneuver is still used as the start, and the bad fit shows in `rms_position` and σ̂₀², or as `NotConverged`.

### Synthetic tests
Truth comes from a known TLE, so SGP4 model error is zero and exact recovery is possible. The layers run from cheapest and most deterministic to statistical.

**Noise-free recovery.** Full states from a truth TLE recover the six orbital elements to near machine precision (RMS < 1 mm), and `residuals[j]` lines up with observation j. With B*₀ equal to the truth, B* is recovered too; otherwise the prior biases it by about (σ_marginal/σ_B)²·(B*₀ − B*_true), which sets its tolerance.

**Deterministic checks** (no Monte Carlo):
- **Velocity term:** the internal velocity sigma equals n₀·σ_model, with n₀ the median candidate mean motion, for a LEO and a GEO fixture, and stays fixed through the fit and every outlier pass. A large velocity glitch on one of the 3 nearest states doesn't change n₀ beyond the spread of the other two.
- **B* prior center:** with `bstar_prior_center` set and an `initial_guess` with a different B*, B*₀ is the option's value. A single LEO state at the epoch returns B* equal to the center exactly; a tight width on a 1-day ISS arc holds B* at the center within the bias formula.
- **B* prior width:** an ISS fit over 1 day gives the same B* with `bstar_prior_width` at 1e-2 and at 1e6, within the bias formula. A GEO fit stays near B*₀ and reports σ_B* ≈ σ_B. A 6-hour GEO fit completes without a divide by zero.
- **Two sources:** two full-state sources, both with model-error-sized residuals, fit together with equal weights; the solution matches the correctly weighted combined fit, neither is rejected, and the result is the same whichever is the majority.

**Positions and measurement noise.**
- **Bookkeeping:** noise-free fits with all states, all positions, and mixed sets (alternating, random subsets, a velocity block at the start or end, a single state among many positions) recover the orbit to RMS < 1 mm and agree with each other; position-only observations report no velocity residual.
- **Information:** P_pos-only − P_mixed is positive semidefinite; trace(P) falls monotonically as the fraction of full states rises; with the states' measurement σ_v very large the mixed fit converges to the position-only result.
- **Noise forms agree:** a diagonal `Covariance` equal to the `Sigma` or `Sigmas` form gives the same solution and P; a rotated, anisotropic covariance changes P as the rotated weights predict; a TEME covariance rotated from RIC sigmas (asymmetric at roundoff) and a singular covariance are accepted; a sigma of 0 gives the same fit as `None`; a covariance with a negative eigenvalue, or asymmetric beyond the tolerance, returns `InvalidObservation { index }`.
- **Two sensors:** GPS-like full states every 10 min (σ_r = 10 m, σ_v = 1 cm/s) and radar-like positions (σ_r = 100 m). With simulated model error at σ_model added to both, combined beats either alone (smaller trace(P), NEES ≈ p) and σ̂₀² ≈ 1.
- **Measurement noise below model error:** states with 10 m measurement noise give the same solution and P as with no measurement noise, within about 0.01%; with measurement noise 10× σ_model, those observations are down-weighted about 100×.
- **Outlier normalization:** the same large outlier in all components, once in a position and once in a state, reaches about the same ρⱼ and is rejected; if the state scores about √2 higher, the √mⱼ normalization is missing.

**Monte Carlo** (fixed seeds, K fits). Add Gaussian noise with exactly σ_model in position and n₀·σ_model in velocity.

| Check | Statistic | Pass if |
|---|---|---|
| Unbiased | Mean of x̂ − x_true | Within ~3·std/√K of 0 |
| Covariance correct | NEES: (x̂ − x)ᵀ·C⁻¹·(x̂ − x) with the returned C, averaged | ≈ p, within χ² bounds for K·p degrees of freedom |
| Noise level honored | σ̂₀² | Mean ≈ 1; (m + 1 − p)·σ̂₀² follows χ²(m + 1 − p) |

Run these on an orbit and arc where the data determines B*, such as LEO over several days. Where the prior dominates, as at GEO, the expected NEES is about p − 1. With few rows σ̂₀² runs slightly low, because the prior row adds almost nothing to J but counts as a degree of freedom, so the χ² check needs m ≫ p.

**Outliers.**

| Case | Expected |
|---|---|
| Clean Gaussian noise, K Monte Carlo fits | `rejected` is empty, and σ̂₀² and NEES are unchanged |
| Noise-free data | `rejected` is empty |
| 1–2 glitches among about 100 observations | Exactly the glitches are rejected, and the result matches a fit with them removed |
| Noise 100× smaller or larger than σ_model, same glitches | Same rejected set as with noise at σ_model |
| A full state with a good position and a glitched velocity | It is rejected (velocity floor) |
| A wrongly flagged point in an early pass | It returns once the fit improves |
| 3 observations, one glitched | ⌊0.25·3⌋ = 0, so nothing is rejected |
| More than `max_rejected_fraction` of observations would be flagged | The previous pass's set is kept, and the fit returns normally |
| A pass's refit doesn't converge | The previous pass's fit and rejected set are returned |
| A rejected set that oscillates between passes | Stops at the pass limit with the lowest-median(ρ) pass |
| The first fit doesn't converge | Rejection is skipped; `NotConverged` |
| Rejection that would leave no full state and fewer than 3 position times | The previous pass's set is kept |

**Edge cases.**

| Observations | Rows | Expected |
|---|---|---|
| None | 1 | `TooFewObservations` |
| 2 positions, or 3 at the same time | 7, 10 | `TooFewObservations` |
| 1 full state | 7 | Fits; `variance_factor` is `None`; covariance is P |
| 3 positions | 10 | Fits |
| 1 full state + 1 position | 10 | Fits |
| 2 full states, same time | 13 | Fits, with a covariance |
| 3 positions microseconds apart, with an initial GP set | 10 | `Underdetermined` |
| A synthetic rank-deficient Jacobian (unit test of the rank check) | – | `Underdetermined`, with no solution |

**Starting guess.**
- A glitch on the state nearest the epoch: another candidate wins the scoring, and the fit matches one with the glitch removed.
- Glitches on the 2 nearest states: the third candidate wins.
- Noise-free data where two candidates score exactly the same: the one nearer the epoch is chosen.
- The two-loop path (loop at the source state's time, propagate, loop again at the epoch), which is the normal path in real use; the Orekit cases only exercise the single-loop path. Each case starts from a noise-free truth TLE and checks that the starting guess reproduces SGP4's propagated source state at the epoch and that LM converges to the truth:
  - **0.43 ms offset:** half a TLE-grid step, as with the default epoch; the guess matches the single-loop result to well below 1 mm.
  - **Several hours:** the B*₀ and re-initialization errors are visible in the guess but removed by LM.
  - **Epoch outside the observation span:** a user epoch days before or after every observation.
  - **Deep space:** a GEO or 12-hour orbit, where re-initializing at the epoch restarts the resonance integration.

**Starting triples.** Velocities only far from the epoch; no velocities at all; the nearest middle point failing coplanarity (the next triple is used); positions from several passes (the timing check rejects straddling triples); an HEO (e ≈ 0.74) pass near perigee and near apogee (the timing check accepts both); a glitch on one position (it appears in at most one triple, and another wins); a glitched single state among many positions (a triple wins); too few positions for 3 independent triples (fewer candidates, the fit still runs); positions where no triple passes the checks and there is no state (`StartingGuessFailed`).

**Held flags.** O1–O5 don't sit near an edge, so these cases are separate, with the truth just inside each gate (225 min, the whole-day gates, the half-day gates, and e = 0.5).
- Truth far from every edge: the flags are the starting flags, and the fit converges.
- Truth just across an edge from the start: the fit stays on its side and returns `NotConverged`.
- A finite-difference step beside a gate stays on the current side, so the Jacobian has no km-over-step spike, and the 3-point one-sided column matches the central column away from the gate.

**Finite-difference steps.** At LEO, near-circular LEO, GEO, HEO and 12-hour fixtures, compare the Jacobian built with steps h, h/2 and 2h; each column should agree to about 1e-6 relative. Where the B* column is below roundoff (GEO over short arcs), check instead that the B* entries of P agree. A perturbed point where SGP4 fails falls back to the one-sided formula on the other side. This test sets c and c_B.

**LM convergence.** Noise-free fits (the Orekit cases and noise-free recovery) must stop by the step test, not by the μ ceiling or the iteration cap; this sets s²_floor. Noisy fits must stop by the step test with `iterations` well below `max_iterations`; this sets `tolerance` and `max_iterations`.

**Solver and rank test.**
- **QR/SVD accuracy:** on ill-conditioned fits (a 10-day GEO arc, a near-circular LEO arc), the step and P agree with a reference built from a Jacobian checked by Richardson extrapolation, to the Jacobian's own accuracy. The same data solved by normal equations is kept as a comparison, to show the digits saved.
- **Rank test:** degenerate cases return `Underdetermined`, and the shortest legitimate arcs pass, among them a single 10-minute LEO pass of positions only, with and without 100 m measurement noise; this sets κ_tol. The smallest singular vector of a degenerate case points at the expected element combination.
- **Re-solve:** after a rejected step, the new μ reuses the SVD, and the step matches a fresh factorization of [Ã; √μ·I] to roundoff.

**Output.**

| Case | Expected |
|---|---|
| Equinoctial → classical → equinoctial | Recovered to roundoff, prograde and retrograde (f_r = ±1) |
| e = 0 exactly; i = 0 (f_r = +1) and i = 180° (f_r = −1) exactly | ω = 0 and Ω = 0, also for negative zeros; at i = 0 and 180°, moving Ω while keeping ω + f_r·Ω propagates the same to 1 mm; at e = 0, another ω differs only within SGP4's 1e-6 eccentricity floor (about 2·a·1e-6) |
| Residual sign | A known +1 km radial offset in an observation gives `position_ric[0]` ≈ +1 |
| Residuals versus the returned `Sgp4` | Propagating `sgp4` and subtracting from each observation reproduces `residuals` exactly |
| `NotConverged` | All fields filled; `residuals` and `rejected` match the last accepted iterate |

**Options.** Each invalid value returns `InvalidOptions { field }` naming that field: `bstar_prior_width` ≤ 0 or not finite; `bstar_prior_center` set but not finite; `outlier_cutoff` below 3 or not finite; `max_rejected_fraction` ≤ 0, ≥ 0.5 or not finite; `max_iterations` = 0; `tolerance` ≤ 0 or not finite; `initial_damping` or `max_damping` ≤ 0 or not finite, or `initial_damping` ≥ `max_damping`, which names `max_damping`; `epoch` not a valid date. With two invalid fields, the first in validation order is named.

**API stability.** Doctests, which compile as an outside crate, build `FitOptions` from `default()` with assigned fields, build mixed observations through `new_fit_observation` and `FitMeasurement` literals, and read every `FitResult` field.

**Practical notes.**
- **Random numbers:** fixed seeds throughout, with a small SplitMix64 + Box–Muller generator in the test module, so fixtures don't shift with `rand` versions.
- **Runtime:** everything except the Monte Carlo runs in the default `cargo test`. Run the Monte Carlo at K ≈ 50 by default, with a K ≈ 500 version behind `#[ignore]`.

### Real-data tests
Every test above uses SGP4-generated truth, so model error never appears. These fixtures test behavior under real model error, and they are required before 0.4.0 is published.

| Fixture | Source | Orbit | What it adds |
|---|---|---|---|
| Swarm | ESA Swarm reduced-dynamic orbit (SP3) | LEO | Observable B* and real drag mismatch |
| GPS | IGS orbits (SP3) | MEO, about 2 rev/day | Deep-space path, with near-resonance error that SGP4 doesn't model at low e |
| QZSS | IGS MGEX orbits (SP3) | GEO altitude, about 1 rev/day | Whole-day resonance path |

**Passing criteria:**
- The outlier cutoff rejects no good data in a maneuver-free window.
- σ̂₀² is stable across neighboring arcs.
- **Covariance realism:** compare prediction error against held-out data within the fit arc with σ̂₀²·P, and record how optimistic it is. The result sets the numbers in the docs' covariance warning.
- **Density:** thin the same data from 1 Hz to one state per 20 min and record how much the formal covariance shrinks.
- **Fit-window guidance:** fit arcs of increasing length and record where the held-out prediction error stops improving. This sets the guidance in open item 13.
- **Position only:** refit every fixture window with each velocity dropped, so each observation becomes a position at the same time. With no full state, the start must come from Gibbs/Herrick–Gibbs triples, so this also exercises the triple checks on real data, with the timing check passing for every accepted triple. The fit converges with default options. Along the arc, its propagated positions agree with the full-state fit's to within that fit's `rms_position`. Record the position-only fit's held-out prediction error and σ̂₀² next to the full-state values.
- **Measurement noise on real data:** repeat the position-only fit with Gaussian noise of about 100 m added to each position and stated as its measurement noise. The elements agree with the clean position-only fit within its covariance, and σ̂₀² stays close to the clean value, because the stated noise is well below the model error.

**Building the fixtures:**
- Confirm each source's current access and license terms before committing data. Commit only a small converted extract.
- Pull a real QZSS file first, and confirm it holds the satellite needed; current availability is unverified.
- Convert from Earth-fixed to TEME offline with Earth orientation data, using a separate tool such as Orekit, and take velocities from that tool if the SP3 files hold only positions. The conversion error must stay well below SGP4's error.
- Choose windows without maneuvers.

### Values set by testing
| Value | Where | Set by |
|---|---|---|
| c, c_B | Finite-difference steps | Step-halving test |
| `tolerance`, `max_iterations`, s²_floor | LM convergence | LM convergence tests |
| κ_tol | Rank test | Rank test |
| tol_r, tol_v, iteration cap | Fixed-point loop | Orekit cases and two-loop tests |
| Pass limit (10) | Outlier loop | Outlier tests |
| θ_switch, span cap, coplanarity tolerance, timing factor | Starting triples | Starting-triple tests |

## Known limitations
These are documented for users rather than handled in code.

- **Fixed model error.** σ_model = 1 km is right for LEO over days but not for every orbit class and arc length. It affects the B* prior balance, the single-state covariance, and how observations with large measurement noise balance against model-dominated ones; σ̂₀² rescales every other covariance.
- **Unstated measurement noise.** An observation without measurement noise is weighted as if its error were SGP4's model error alone, so a source much noisier than that (a poor receiver) gets more weight than it should unless its noise is given.
- **Small glitches.** Glitches beyond the outlier cutoff (in LEO with the default cutoff and 1 km residuals, about 23 km in a full state and 15 km in a position) are rejected; smaller ones stay in at full weight.
- **Fits across a gate edge.** A fit whose answer lies across a deep-space or resonance edge returns `NotConverged`, after running to the iteration cap or the μ ceiling. There is no reliable workaround: the flags come from the fixed-point result at the epoch, not directly from `initial_guess`. Crossing an edge moves the propagated orbit discontinuously (measured between orbits about 1e-15 rev/day apart: 2–6 km at the 225-min boundary, including at epoch; zero at epoch at the resonance gates, growing to 4–152 km by 14 days), so the fit can't simply ignore the flags. Such orbits are rare: GEO objects sit well inside the whole-day window, GPS isn't half-day resonant (e ≪ 0.5), and few objects have a period near 225 min.
- **The covariance is optimistic on real data.** σ̂₀²·P assumes independent residuals; real SGP4 model error is correlated in time, so it is too small, more so the denser the sampling (about 35× for a day of 1 Hz GPS). The docs must say so and recommend thinning dense data to about one state per minute. No covariance describes prediction beyond the fit arc.
- **Maneuvers inside the arc.** A single GP set can't describe the orbit before and after a burn, so the fit returns a compromise with a large `rms_position` and no error. Outlier rejection doesn't catch it, because median(ρ) is itself large and the 25% cap stops rejection long before half the data could go. An automatic detector isn't provided: a runs test on residual signs would false-alarm on real data, because SGP4 model error is itself smooth and systematic. This matches practice elsewhere: CelesTrak fits operator ephemerides and issues a separate `[PM]` set after a known maneuver; the Naval Space Command AUTODC process checks RMS over the fit span and its subsets and refers failures to analysts; Orekit models maneuvers as estimated parameters, which SGP4 can't. The docs must say: with known maneuver times, split the observations at each one and fit each segment separately; without them, an RMS well above what is typical for the orbit class is the warning sign, and a plot of along-track residuals against time shows a kink at the maneuver. Possible later additions, outside `fit_gp`: an RMS over subsets of the arc from the existing residuals, or a split-arc helper that fits each half separately and compares.
- **Sparse position-only data.** Without a full state, the start needs a triple: three positions within the span cap (about 60–90° of arc) that pass the timing check. Sparse tracking, such as one or two optical points per night on a GEO object, never gives one, so without `initial_guess` the fit returns `StartingGuessFailed`. That is the intended response; the docs must say that such data needs an `initial_guess`.
- **Isotropic model error.** Mⱼ adds σ_model equally in every direction, but SGP4's real error is mostly along-track and much smaller radially. A sensor that is accurate in one direction, such as radar range, is still weighted as if its error there were σ_model, so the fit uses less of that sensor's information than it could. Sources far better than the model error (the goal's GPS and ephemeris sources) are unaffected. A model-error term shaped in the candidate orbit's RIC axes could be added later without an API change.
- **No progressive arc extension.** A start more than half a revolution-per-arc valley out in n (about 30 m/s in velocity for a 3-day LEO arc) can converge to the wrong revolution count. Starts built from the data are far better than that; a stale `initial_guess` is the user's choice.
- **No real-data test of the half-day resonance path.** The fixtures cover LEO, MEO and the whole-day resonance, but no Molniya-type orbit. Precise public HEO ephemerides are hard to find; look for one when building the fixtures, and otherwise rely on the SGP4-generated held-flag tests.

---

## Open items
- **13. Fit-window guidance (docs).** Recommended arc lengths per orbit class (LEO, MEO, 12-hour, HEO, GEO), set from the real-data fixtures: held-out prediction error shows where a longer arc stops helping. LEO 1–3 days and GEO 5–10 days are placeholders.
- **14. Deferred: covariance extras.** Conversion to Kozai elements (P_K = G·P·Gᵀ) and a 6×6 position and velocity covariance at time t (C = B·P·Bᵀ). The 6×6 helper leaves out most SGP4 model error and is likely overconfident past the fit span, so it waits until a real-data fixture shows how far to trust it. Adding it later breaks nothing.
- **18. Measure the remaining branch jumps.** The perigee-220 km (simplified drag), e = 1e-4 (C3 and δω/δM dropped) and Lyddane (i < 0.2 rad in deep space) branches aren't stored flags, so the fit doesn't hold them. Measure the jump across each once; if it is large, extend the one-sided finite-difference rule to that edge. The Lyddane switch is tested at each propagation time, not at the epoch: in an arc near i = 0.2 rad (O4 sits at 10.8°, 0.66° below it, and old GEO satellites drift through it), some observations can use one branch and some the other, and a finite-difference step can flip a single observation's branch. Measure that jump per observation, and if it is large, apply the one-sided rule to the affected rows, not to whole columns.
- **19. Near-zero eccentricity.** SGP4 replaces e < 1e-6 by 1e-6 along ω ([sgp4.rs](../src/sgp4.rs), the eccentricity floor in the secular update). Below that, the propagated orbit depends on the direction of (a_f, a_g) but not on its length, so the Jacobian columns for a_f and a_g lose their length direction and become large in the angular one as e → 0. A truth or iterate with e < 1e-6 is rare in real data (GEO and GPS sit near 1e-4), but a synthetic e = 0 truth hits it. Measure how LM behaves on an e = 0 truth when building step 6; if it stalls or reports `Underdetermined`, decide between documenting it and holding e at the floor.

## Decision log
Closed items, in the order they were settled. The current design is in the sections above; these record why.

- **1. Initial guess and metadata.** Replaced a `FitStart` enum with separate `initial_guess` and `metadata` options, so a user can supply either, both, or neither, at the cost of a per-field precedence rule.
- **2. Sigma weighting.** Began as global user sigmas (1 km, σ_v derived) and became an internal model-error term plus optional per-observation measurement noise (items 22, 24, 26).
- **3. Mode boundaries.** Each fit is held to one flag set. A second complete fit per nearby edge was designed and later dropped (item 27).
- **4. RIC residuals.** `FitResidual` reports both TEME and RIC.
- **5. Covariance output.** One covariance, σ̂₀²·P (P for a single state), plus σ̂₀². A robust covariance was designed and later dropped (item 27).
- **6. B\* always fitted.** The covariance is always the full 7×7 with no zeroed row or column.
- **7. Outlier handling.** Median-scaled cutoff, no Huber. Errors between about 3σ and 10σ stay in at full weight; Huber can be added later if real data shows a need.
- **8. Finite-difference steps.** Per-element scales with a larger B* coefficient. A 1e-7 B* step gave an exactly zero column at GEO over 6 hours; the prior row makes that harmless.
- **9. Convergence criteria.** Mahalanobis step test with a noise-free floor; `initial_damping` and `max_damping` are options.
- **10. Input validation.** Unsorted input and duplicate timestamps are accepted; zero degrees of freedom occurs only for a single state, where σ̂₀² is `None`.
- **11. Single state.** `gp_from_state` was dropped; a single state goes through `fit_gp` as a square system.
- **12. Observation type.** A time plus a `#[non_exhaustive]` `FitMeasurement` enum, so new kinds don't break the API. Its variants carry their fields inline (decision 31).
- **15. Rev number at a moved epoch.** Advanced by whole revolutions, as Orekit does.
- **16. Starting guess.** Up to 3 independent candidates scored by median ρ: the nearest states, topped up with Gibbs/Herrick–Gibbs triples.
- **17. Weakly observable B\*.** Solved by the B* prior. B*'s effect on position falls steadily with altitude and grows about as t². At GEO, B* = 1 shifts position by about 3 mm at 6 hours and 6 m at 14 days, and the 1σ B* over 14 days is about 4 (position-only, σ_r = 100 m). The default prior also limits how much B* can absorb unmodeled forces at high altitude; users who want that can loosen `bstar_prior_width`.
- **19. Progressive arc extension: not adopted.** It would guard against a start more than half a revolution-per-arc valley out in n. Data-built starts are far better than that, and a stale `initial_guess` is the user's choice. Recorded as a known limitation.
- **21. 0.4.0 scope.** The plan had grown into a general multi-sensor orbit-determination design and was split into a core 0.4.0 and roadmap steps; the roadmap was later cut (item 27) and its remaining step merged back (item 28).
- **22. Derived velocity sigma.** σ_v = n·σ_r (now n₀·σ_model, item 26), because a position error σ_r on an orbit goes with about n·σ_r of velocity error; a fixed 1 m/s gave GEO velocities almost no weight.
- **23. Outlier rejection in 0.4.0.** In its simple form: global median, 10× cutoff, 1 mm and 1 µm/s floors, 25% cap. It is scale-free and small, real GPS receiver data does glitch, and because measurement noise is added to the model-error term (item 26), it can't swamp a source.
- **24. No per-observation noise in 0.4.0** (reversed by item 28). With model error dominating, per-observation sigmas used as weights could only mis-state relative accuracy; item 26 resolved this by adding measurement noise to the model-error term instead.
- **25. Plan cleanup.** Reorganized into a design section and per-step roadmap sections, with no decision changed.
- **26. Measurement noise and model error separated.** `sigma_position` and `sigma_velocity` were documented as SGP4 model error but set by the user, with a warning not to set them to sensor noise. Both were removed: σ_model = 1 km is internal, with velocity n₀·σ_model, and measurement noise is added to it per observation. σ_model is fixed because an estimate collapses on noise-free data and frees B* from its prior. Removing the options before release is free; adding one back later isn't breaking.
- **27. Roadmap cut to step 1.** Removed: alternative flag sets with the stall stop (complex, for a rare case, with a proximity window that relied on the optimistic covariance); the Newey–West robust covariance (already conditional on real data, and the covariance is only a diagnostic); and noise groups (superseded by item 26, taking item 20 with them). The gate-edge limitation and its measured jumps were kept, and the real-data realism test now sets the docs' covariance warning. Any of them can return without a breaking change.
- **28. Roadmap 1 merged into 0.4.0.** Positions, Gibbs/Herrick–Gibbs triples and measurement noise were the only roadmap step left, and the 0.4.0 API was already shaped around them. They moved in unchanged except for the `noise` argument on `FitStateObs::new`, `StartingGuessFailed` covering data that yields no candidate, and the position outlier cutoff. The Scope table became implementation steps.
- **29. Red team after the merge.** Noise validation relaxed: sigmas ≥ 0, covariances symmetric to roundoff and positive semidefinite. Sparse position-only data returning `StartingGuessFailed` is intended and documented. The isotropic model-error term is a recorded limitation. Added a real-data measurement-noise test and a short position-only pass to the rank-test fixtures.
- **30. no_std.** The fit follows the crate's `#![no_std]` + `alloc` setup: `libm` for all math, no new dependency, and a `--no-default-features` build in every implementation step.
- **31. Step 3 (API and validation).** One constructor, the free function `new_fit_observation(datetime, measurement)`, for every kind of measurement. The `FitStateObs` and `FitPositionObs` payload structs were dropped: `FitMeasurement::State` and `::Position` carry their fields inline and are built as struct literals, which a `#[non_exhaustive]` enum allows outside the crate. The cost is that a new field on an existing variant would be breaking, which is unlikely because noise forms, frames and new kinds each have their own non-exhaustive enum. The per-kind shortcuts (`from_state`, `from_position`) were dropped too. The covariance check scales each position or velocity block by its own variance, because one trace-based shift hides velocity-block errors. `FitError` has a hand-written `Debug`, because `Sgp4` has none. The noise enums allow `clippy::large_enum_variant`: boxing a 6×6 covariance would cost an allocation per observation. The fit types are re-exported at the crate root now; `fit_gp` follows in step 7.

## References
- D. A. Vallado and P. Crawford, *SGP4 Orbit Determination*, AIAA 2008-6770.
- K. Levenberg, "A Method for the Solution of Certain Non-Linear Problems in Least Squares," 1944.
- D. W. Marquardt, "An Algorithm for Least-Squares Estimation of Nonlinear Parameters," 1963.
- H. B. Nielsen, "Damping Parameter in Marquardt's Method," IMM-REP-1999-05, Technical University of Denmark, 1999.
- Orekit test suite: `src/test/java/org/orekit/propagation/analytical/tle/generation/` ([GitHub mirror](https://github.com/CS-SI/Orekit)).
- Citation details for Nielsen were written from memory and must be checked before the design moves into the math spec.
