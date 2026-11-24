//! A bound-constrained minimiser by sequential quadratic programming (SLSQP), the engine of
//! the sketch solver ([`super::sketch::solve`]).
//!
//! The method is Kraft's SLSQP as published (D. Kraft, *A software package for sequential
//! quadratic programming*, DFVLR-FB 88-28, 1988; ACM TOMS 20, 1994), restricted to simple
//! bounds `lower ≤ x ≤ upper`:
//!
//! - the Hessian approximation `B` is kept as `L D Lᵀ` (unit lower `L`, positive diagonal
//!   `D`), started at the identity and updated by BFGS with Powell's damping (θ from the
//!   0.2 / 0.8 rule), applied to the factors as two rank-one modifications (Fletcher & Powell;
//!   Gill, Golub, Murray & Saunders), the second, negative one with the backward recurrence
//!   that keeps the factors positive definite;
//! - the quadratic subproblem `min ½ dᵀBd + gᵀd, lower − x ≤ d ≤ upper − x` is solved as a
//!   least-squares problem with inequality constraints in those factors,
//!   `min ‖D^½ Lᵀ d + D^-½ L⁻¹ g‖`, reduced to a least-distance problem and solved by
//!   Lawson & Hanson's NNLS (*Solving Least Squares Problems*, ch. 23); only finite bounds
//!   become constraints;
//! - the step length starts at 1 and is accepted by Armijo's test with η = 0.1 (the merit
//!   function is the objective: there are no general constraints); otherwise it is shrunk by
//!   safeguarded quadratic interpolation, `α ← max(0.1 α, α_q)`, or halved after a
//!   non-finite value; after 10 finite trials the last one is accepted;
//! - when the subproblem gives no descent direction the Hessian approximation is restarted
//!   from the identity (at most [`MAX_RESTARTS`] times in a run); if there is still none
//!   the run ends as converged (`FtolReached` / `XtolReached` by tolerances relaxed tenfold,
//!   else `RoundoffLimited`).
//!
//! Stopping tests (NLopt's conventions, which CadQuery's settings are written for) are made
//! between successive accepted iterates, from the second one on. Every point evaluated lies
//! within the bounds. The result is the best point evaluated.
//!
//! Safety caps, which report `Failure` and are far beyond anything a sound run reaches:
//! [`MAX_ITERATIONS`] iterations, [`MAX_TRIALS_TOTAL`] trials in one line search (only
//! non-finite values can make a search that long), and NNLS's own iteration limit
//! (Lawson & Hanson's 3 × the number of constraints, here at least 100).

/// When to stop (all tests between successive accepted iterates; 0 / empty disables).
#[derive(Debug, Clone)]
pub(crate) struct Stopping {
    pub ftol_rel: f64,
    pub ftol_abs: f64,
    pub xtol_rel: f64,
    pub xtol_abs: Vec<f64>,
    pub max_evaluations: usize,
}

/// How a run ended (NLopt's public result codes, which CadQuery reports).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Result {
    Success = 1,
    FtolReached = 3,
    XtolReached = 4,
    MaxevalReached = 5,
    Failure = -1,
    InvalidArgs = -2,
    RoundoffLimited = -4,
}

/// What [`minimize`] found: the best point evaluated and its value.
#[derive(Debug, Clone)]
pub(crate) struct Outcome {
    pub x: Vec<f64>,
    pub minimum: f64,
    pub evaluations: usize,
    pub result: Result,
}

/// Armijo's sufficient-decrease constant.
const ARMIJO: f64 = 0.1;
/// The least factor a step is shrunk by when interpolating.
const SHRINK_LIMIT: f64 = 0.1;
/// Finite trials per line search; the last is accepted.
const LINE_SEARCH_TRIALS: usize = 10;
/// Restarts of the Hessian approximation from the identity in one run.
const MAX_RESTARTS: usize = 5;
/// The factor on the tolerances when no descent direction is left.
const RELAXED: f64 = 10.0;
/// Safety cap on iterations (reports `Failure`).
const MAX_ITERATIONS: usize = 1_000_000;
/// Safety cap on trials in one line search, non-finite ones included (reports `Failure`).
const MAX_TRIALS_TOTAL: usize = 2_000;

/// `objective` made least within `lower ≤ x ≤ upper` from `x0`, by SLSQP (module doc).
/// `objective(x, Some(g))` returns f(x) and fills g with ∇f(x); with `None` only f(x).
pub(crate) fn minimize(
    x0: &[f64],
    lower: &[f64],
    upper: &[f64],
    mut objective: impl FnMut(&[f64], Option<&mut [f64]>) -> f64,
    stop: &Stopping,
) -> Outcome {
    let n = x0.len();
    assert!(lower.len() == n && upper.len() == n, "slsqp: bounds of the wrong length");
    assert!(stop.xtol_abs.is_empty() || stop.xtol_abs.len() == n, "slsqp: xtol_abs of the wrong length");
    if n == 0 {
        let value = objective(&[], None);
        return Outcome { x: Vec::new(), minimum: value, evaluations: 1, result: Result::Success };
    }
    for i in 0..n {
        if lower[i] > upper[i] || x0[i] < lower[i] || x0[i] > upper[i] {
            return Outcome { x: x0.to_vec(), minimum: f64::INFINITY, evaluations: 0, result: Result::InvalidArgs };
        }
    }
    Solver::new(x0, lower, upper, stop).run(&mut objective)
}

/// The state of a run.
struct Solver<'a> {
    n: usize,
    lower: &'a [f64],
    upper: &'a [f64],
    stop: &'a Stopping,
    /// The accepted iterate, its value and gradient.
    x: Vec<f64>,
    f: f64,
    g: Vec<f64>,
    /// `B = L D Lᵀ`, `L` unit lower triangular stored by rows (n × n).
    l: Vec<f64>,
    d: Vec<f64>,
    identity: bool,
    restarts: usize,
    evaluations: usize,
    best_x: Vec<f64>,
    best_f: f64,
    last_x: Vec<f64>,
    last_f: f64,
}

/// How the evaluation bookkeeping says to go on.
enum Next {
    Go,
    End(Result),
}

impl<'a> Solver<'a> {
    fn new(x0: &[f64], lower: &'a [f64], upper: &'a [f64], stop: &'a Stopping) -> Self {
        let n = x0.len();
        let mut s = Solver {
            n,
            lower,
            upper,
            stop,
            x: x0.to_vec(),
            f: f64::INFINITY,
            g: vec![0.0; n],
            l: vec![0.0; n * n],
            d: vec![1.0; n],
            identity: true,
            restarts: 0,
            evaluations: 0,
            best_x: x0.to_vec(),
            best_f: f64::INFINITY,
            last_x: x0.to_vec(),
            last_f: f64::INFINITY,
        };
        s.reset_hessian();
        s
    }

    fn reset_hessian(&mut self) {
        let n = self.n;
        self.l.iter_mut().for_each(|v| *v = 0.0);
        for i in 0..n {
            self.l[i * n + i] = 1.0;
        }
        self.d.iter_mut().for_each(|v| *v = 1.0);
        self.identity = true;
    }

    /// One call of the objective at `x` with the gradient into `g`, recorded.
    fn evaluate(&mut self, objective: &mut impl FnMut(&[f64], Option<&mut [f64]>) -> f64, x: &[f64], g: &mut [f64]) -> f64 {
        let value = objective(x, Some(g));
        self.evaluations += 1;
        if value < self.best_f {
            self.best_f = value;
            self.best_x.copy_from_slice(x);
        }
        self.last_x.copy_from_slice(x);
        self.last_f = value;
        value
    }

    fn limit_reached(&self) -> bool {
        self.stop.max_evaluations > 0 && self.evaluations >= self.stop.max_evaluations
    }

    fn outcome(&self, result: Result) -> Outcome {
        if self.best_f.is_finite() {
            Outcome { x: self.best_x.clone(), minimum: self.best_f, evaluations: self.evaluations, result }
        } else {
            Outcome { x: self.last_x.clone(), minimum: self.last_f, evaluations: self.evaluations, result }
        }
    }

    fn run(mut self, objective: &mut impl FnMut(&[f64], Option<&mut [f64]>) -> f64) -> Outcome {
        let n = self.n;
        let x0 = self.x.clone();
        let mut g = vec![0.0; n];
        self.f = self.evaluate(objective, &x0, &mut g);
        self.g.copy_from_slice(&g);
        if self.limit_reached() {
            return self.outcome(Result::MaxevalReached);
        }
        if !self.f.is_finite() || self.g.iter().any(|v| !v.is_finite()) {
            return self.outcome(Result::RoundoffLimited);
        }
        let mut subproblem = Subproblem::new(n);
        let mut direction = vec![0.0; n];
        let mut trial = vec![0.0; n];
        let mut trial_g = vec![0.0; n];
        let mut first = true;
        for _ in 0..MAX_ITERATIONS {
            // A descent direction from the quadratic model, restarting B if there is none.
            let slope = loop {
                let solved = subproblem.solve(&self, &mut direction);
                let slope = match solved {
                    Ok(()) => dot(&self.g, &direction),
                    Err(e) => {
                        if !self.identity && self.restarts < MAX_RESTARTS {
                            self.restarts += 1;
                            self.reset_hessian();
                            continue;
                        }
                        return self.outcome(e);
                    }
                };
                if slope.is_nan() {
                    return self.outcome(Result::RoundoffLimited);
                }
                if slope < 0.0 {
                    break slope;
                }
                if !self.identity && self.restarts < MAX_RESTARTS {
                    self.restarts += 1;
                    self.reset_hessian();
                    continue;
                }
                // No progress possible: converged, by tolerances relaxed, between the
                // accepted iterate and itself.
                let x = self.x.clone();
                return self.outcome(self.tolerances(&x, self.f, &x, self.f, RELAXED).unwrap_or(Result::RoundoffLimited));
            };

            // Line search on f along the direction.
            let mut alpha = 1.0;
            let mut finite_trials = 0;
            let mut trials = 0;
            let f_trial = loop {
                trials += 1;
                if trials > MAX_TRIALS_TOTAL {
                    return self.outcome(Result::Failure);
                }
                for i in 0..n {
                    trial[i] = (self.x[i] + alpha * direction[i]).max(self.lower[i]).min(self.upper[i]);
                }
                let value = self.evaluate(objective, &trial, &mut trial_g);
                if !value.is_finite() {
                    if self.limit_reached() {
                        return self.outcome(Result::MaxevalReached);
                    }
                    alpha *= 0.5;
                    continue;
                }
                finite_trials += 1;
                let decrease = value - self.f;
                if decrease <= ARMIJO * alpha * slope || finite_trials >= LINE_SEARCH_TRIALS {
                    break value;
                }
                if self.limit_reached() {
                    return self.outcome(Result::MaxevalReached);
                }
                // The minimum of the quadratic through f(0), f'(0) and f(α).
                let q = 0.5 * slope * alpha * alpha / (slope * alpha - decrease);
                alpha = if q.is_finite() { q.max(SHRINK_LIMIT * alpha) } else { SHRINK_LIMIT * alpha };
            };

            // Accepted: the stopping tests, then the evaluation limit.
            let stopped = if first { None } else { self.tolerances(&trial, f_trial, &self.x, self.f, 1.0) };
            first = false;
            if self.limit_reached() {
                return self.outcome(Result::MaxevalReached);
            }
            if let Some(result) = stopped {
                return self.outcome(result);
            }
            if trial_g.iter().any(|v| !v.is_finite()) {
                return self.outcome(Result::RoundoffLimited);
            }
            let s: Vec<f64> = (0..n).map(|i| trial[i] - self.x[i]).collect();
            let y: Vec<f64> = (0..n).map(|i| trial_g[i] - self.g[i]).collect();
            self.x.copy_from_slice(&trial);
            self.f = f_trial;
            self.g.copy_from_slice(&trial_g);
            self.damped_bfgs(&s, &y);
        }
        self.outcome(Result::Failure)
    }

    /// The tolerance tests between a new iterate and the one before, the tolerances scaled
    /// by `factor` (`xtol_abs` excepted).
    fn tolerances(&self, x_new: &[f64], f_new: f64, x_old: &[f64], f_old: f64, factor: f64) -> Option<Result> {
        let st = self.stop;
        if f_old.is_finite() && f_new.is_finite() {
            let change = (f_new - f_old).abs();
            let abs = st.ftol_abs > 0.0 && change < factor * st.ftol_abs;
            let rel = st.ftol_rel > 0.0 && (f_new == f_old || change < factor * st.ftol_rel * 0.5 * (f_new.abs() + f_old.abs()));
            if abs || rel {
                return Some(Result::FtolReached);
            }
        }
        if st.xtol_rel > 0.0 {
            let step: f64 = x_new.iter().zip(x_old).map(|(a, b)| (a - b).abs()).sum();
            let size: f64 = x_new.iter().map(|a| a.abs()).sum();
            if step < factor * st.xtol_rel * size {
                return Some(Result::XtolReached);
            }
        }
        if !st.xtol_abs.is_empty() && x_new.iter().zip(x_old).zip(&st.xtol_abs).all(|((a, b), t)| (a - b).abs() < *t) {
            return Some(Result::XtolReached);
        }
        None
    }

    /// `B s` from the factors.
    fn times_b(&self, s: &[f64]) -> Vec<f64> {
        let n = self.n;
        // u = D Lᵀ s
        let mut u = vec![0.0; n];
        for j in 0..n {
            let mut sum = s[j];
            for i in j + 1..n {
                sum += self.l[i * n + j] * s[i];
            }
            u[j] = self.d[j] * sum;
        }
        // L u
        let mut out = vec![0.0; n];
        for i in 0..n {
            let mut sum = u[i];
            for j in 0..i {
                sum += self.l[i * n + j] * u[j];
            }
            out[i] = sum;
        }
        out
    }

    /// BFGS with Powell's damping on the factors: `B + ηηᵀ/sᵀη − BssᵀB/sᵀBs`.
    fn damped_bfgs(&mut self, s: &[f64], y: &[f64]) {
        let bs = self.times_b(s);
        let sbs = dot(s, &bs);
        if !(sbs > 0.0) || !sbs.is_finite() {
            return;
        }
        let sy = dot(s, y);
        let eta: Vec<f64> = if sy < 0.2 * sbs {
            let theta = 0.8 * sbs / (sbs - sy);
            (0..self.n).map(|i| theta * y[i] + (1.0 - theta) * bs[i]).collect()
        } else {
            y.to_vec()
        };
        let seta = dot(s, &eta);
        if !(seta > 0.0) || !seta.is_finite() {
            return;
        }
        self.rank_one(&eta, 1.0 / seta);
        self.rank_one(&bs, -1.0 / sbs);
        self.identity = false;
    }

    /// `L D Lᵀ + σ z zᵀ` refactored in place (Gill, Golub, Murray & Saunders; for σ < 0
    /// Fletcher & Powell's backward recurrence, which keeps D positive).
    fn rank_one(&mut self, z: &[f64], sigma: f64) {
        let n = self.n;
        // v = L⁻¹ z
        let mut v = z.to_vec();
        for i in 0..n {
            let mut sum = v[i];
            for j in 0..i {
                sum -= self.l[i * n + j] * v[j];
            }
            v[i] = sum;
        }
        let mut t = vec![0.0; n + 1];
        if sigma > 0.0 {
            t[0] = 1.0 / sigma;
            for j in 0..n {
                t[j + 1] = t[j] + v[j] * v[j] / self.d[j];
            }
        } else {
            let mut last = 1.0 / sigma;
            for j in 0..n {
                last += v[j] * v[j] / self.d[j];
            }
            let floor = f64::EPSILON / sigma;
            if last > floor {
                last = floor;
            }
            t[n] = last;
            for j in (0..n).rev() {
                t[j] = t[j + 1] - v[j] * v[j] / self.d[j];
            }
        }
        let mut w = z.to_vec();
        for j in 0..n {
            let dj = self.d[j];
            let beta = v[j] / (dj * t[j + 1]);
            self.d[j] = dj * t[j + 1] / t[j];
            for r in j + 1..n {
                w[r] -= v[j] * self.l[r * n + j];
                self.l[r * n + j] += beta * w[r];
            }
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// The quadratic subproblem as a least-distance problem (work space kept between solves).
struct Subproblem {
    n: usize,
    /// Rows of E⁻¹ = L⁻ᵀ D^-½ for constrained variables, and the free step −B⁻¹g.
    row: Vec<f64>,
    free: Vec<f64>,
}

impl Subproblem {
    fn new(n: usize) -> Self {
        Subproblem { n, row: vec![0.0; n], free: vec![0.0; n] }
    }

    /// The step `d` minimising `½ dᵀBd + gᵀd` within `lower − x ≤ d ≤ upper − x`.
    fn solve(&mut self, s: &Solver, d: &mut [f64]) -> std::result::Result<(), Result> {
        let n = self.n;
        if s.d.iter().any(|&v| !(v > 0.0) || !v.is_finite()) {
            return Err(Result::RoundoffLimited);
        }
        // free = −L⁻ᵀ D⁻¹ L⁻¹ g
        let mut u = s.g.clone();
        forward(&s.l, n, &mut u);
        for j in 0..n {
            u[j] /= s.d[j];
        }
        backward(&s.l, n, &mut u);
        for j in 0..n {
            self.free[j] = -u[j];
        }
        // The finite bounds as rows of G̃ z ≥ h̃, with z = E d − E·free (E = D^½ Lᵀ).
        let mut columns: Vec<Vec<f64>> = Vec::new();
        for i in 0..n {
            let lo = s.lower[i] - s.x[i];
            let hi = s.upper[i] - s.x[i];
            let has_lo = s.lower[i].is_finite();
            let has_hi = s.upper[i].is_finite();
            if !has_lo && !has_hi {
                continue;
            }
            // Row i of E⁻¹: (L⁻¹ e_i)_k / √d_k.
            self.row.iter_mut().for_each(|v| *v = 0.0);
            self.row[i] = 1.0;
            for r in i + 1..n {
                let mut sum = 0.0;
                for k in i..r {
                    sum -= s.l[r * n + k] * self.row[k];
                }
                self.row[r] = sum;
            }
            for k in i..n {
                self.row[k] /= s.d[k].sqrt();
            }
            if has_lo {
                let mut c: Vec<f64> = self.row.clone();
                c.push(lo - self.free[i]);
                columns.push(c);
            }
            if has_hi {
                let mut c: Vec<f64> = self.row.iter().map(|v| -v).collect();
                c.push(-(hi - self.free[i]));
                columns.push(c);
            }
        }
        let z = if columns.iter().all(|c| c[n] <= 0.0) { vec![0.0; n] } else { least_distance(&columns, n)? };
        // d = E⁻¹ z + free = L⁻ᵀ D^-½ z + free
        let mut e = z;
        for j in 0..n {
            e[j] /= s.d[j].sqrt();
        }
        backward(&s.l, n, &mut e);
        for j in 0..n {
            d[j] = e[j] + self.free[j];
        }
        if d.iter().any(|v| !v.is_finite()) {
            return Err(Result::RoundoffLimited);
        }
        Ok(())
    }
}

/// u ← L⁻¹ u (unit lower, by rows).
fn forward(l: &[f64], n: usize, u: &mut [f64]) {
    for i in 0..n {
        let mut sum = u[i];
        for j in 0..i {
            sum -= l[i * n + j] * u[j];
        }
        u[i] = sum;
    }
}

/// u ← L⁻ᵀ u.
fn backward(l: &[f64], n: usize, u: &mut [f64]) {
    for i in (0..n).rev() {
        let mut sum = u[i];
        for j in i + 1..n {
            sum -= l[j * n + i] * u[j];
        }
        u[i] = sum;
    }
}

/// Least distance programming (Lawson & Hanson ch. 23): the least `‖z‖` with
/// `Σ_k c_k[..n]·z ≥ c_k[n]` for each column `c_k`, by NNLS on `[G̃ᵀ; h̃ᵀ] u ≈ e_{n+1}`.
fn least_distance(columns: &[Vec<f64>], n: usize) -> std::result::Result<Vec<f64>, Result> {
    let mut target = vec![0.0; n + 1];
    target[n] = 1.0;
    let u = nnls(columns, &target)?;
    let mut r = vec![0.0; n + 1];
    for (c, &uk) in columns.iter().zip(&u) {
        if uk != 0.0 {
            for i in 0..=n {
                r[i] += c[i] * uk;
            }
        }
    }
    r[n] -= 1.0;
    let norm: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
    if !(norm > 1e3 * f64::EPSILON) || r[n].abs() <= f64::EPSILON * norm {
        // The constraints are inconsistent (cannot happen with a start within the bounds).
        return Err(Result::Failure);
    }
    Ok((0..n).map(|i| -r[i] / r[n]).collect())
}

/// Lawson & Hanson's NNLS: `u ≥ 0` least `‖A u − b‖`, `A` given by columns.
fn nnls(columns: &[Vec<f64>], b: &[f64]) -> std::result::Result<Vec<f64>, Result> {
    let m = columns.len();
    let rows = b.len();
    let mut u = vec![0.0; m];
    let mut passive: Vec<usize> = Vec::new();
    let mut in_passive = vec![false; m];
    let mut blocked = vec![false; m];
    let max_iterations = (3 * m).max(100);
    let mut iterations = 0;
    let residual = |u: &[f64]| -> Vec<f64> {
        let mut r = b.to_vec();
        for (c, &uk) in columns.iter().zip(u) {
            if uk != 0.0 {
                for i in 0..rows {
                    r[i] -= c[i] * uk;
                }
            }
        }
        r
    };
    loop {
        let r = residual(&u);
        // The dual vector w = Aᵀ(b − A u); the largest positive one among the zero set.
        let mut best: Option<(usize, f64)> = None;
        for k in 0..m {
            if in_passive[k] || blocked[k] {
                continue;
            }
            let w = dot(&columns[k], &r);
            if w > 0.0 && best.map_or(true, |(_, bw)| w > bw) {
                best = Some((k, w));
            }
        }
        let Some((t, _)) = best else { return Ok(u) };
        // A column dependent on the passive ones cannot enter.
        let mut trial_set = passive.clone();
        trial_set.push(t);
        let Some(z) = least_squares_on(columns, &trial_set, b) else {
            blocked[t] = true;
            continue;
        };
        if !(z[z.len() - 1] > 0.0) {
            blocked[t] = true;
            continue;
        }
        passive.push(t);
        in_passive[t] = true;
        blocked.iter_mut().for_each(|v| *v = false);
        let mut z = z;
        loop {
            iterations += 1;
            if iterations > max_iterations {
                return Err(Result::Failure);
            }
            if z.iter().all(|&v| v > 0.0) {
                for (p, &k) in passive.iter().enumerate() {
                    u[k] = z[p];
                }
                break;
            }
            let mut alpha = f64::INFINITY;
            for (p, &k) in passive.iter().enumerate() {
                if z[p] <= 0.0 {
                    let a = u[k] / (u[k] - z[p]);
                    if a < alpha {
                        alpha = a;
                    }
                }
            }
            for (p, &k) in passive.iter().enumerate() {
                u[k] += alpha * (z[p] - u[k]);
            }
            let mut kept = Vec::with_capacity(passive.len());
            for (p, &k) in passive.iter().enumerate() {
                if z[p] <= 0.0 && u[k] <= 1e-15 * (1.0 + u[k].abs()) || u[k] <= 0.0 {
                    u[k] = 0.0;
                    in_passive[k] = false;
                } else {
                    kept.push(k);
                }
            }
            passive = kept;
            if passive.is_empty() {
                break;
            }
            match least_squares_on(columns, &passive, b) {
                Some(nz) => z = nz,
                None => return Err(Result::RoundoffLimited),
            }
        }
    }
}

/// The least-squares coefficients of `b` on the columns `set` (Householder QR), or `None`
/// when the columns are (numerically) dependent.
fn least_squares_on(columns: &[Vec<f64>], set: &[usize], b: &[f64]) -> Option<Vec<f64>> {
    let rows = b.len();
    let p = set.len();
    if p > rows {
        return None;
    }
    let mut a: Vec<Vec<f64>> = set.iter().map(|&k| columns[k].clone()).collect();
    let norms: Vec<f64> = a.iter().map(|c| c.iter().map(|v| v * v).sum::<f64>().sqrt()).collect();
    let mut rhs = b.to_vec();
    for j in 0..p {
        let norm: f64 = a[j][j..].iter().map(|v| v * v).sum::<f64>().sqrt();
        if !(norm > 1e-12 * norms[j]) {
            return None;
        }
        let alpha = if a[j][j] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = a[j][j..].to_vec();
        v[0] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv > 0.0 {
            for col in a.iter_mut().skip(j) {
                let s: f64 = v.iter().zip(&col[j..]).map(|(x, y)| x * y).sum::<f64>() * 2.0 / vv;
                for (i, vi) in v.iter().enumerate() {
                    col[j + i] -= s * vi;
                }
            }
            let s: f64 = v.iter().zip(&rhs[j..]).map(|(x, y)| x * y).sum::<f64>() * 2.0 / vv;
            for (i, vi) in v.iter().enumerate() {
                rhs[j + i] -= s * vi;
            }
        }
    }
    let mut z = vec![0.0; p];
    for j in (0..p).rev() {
        let mut sum = rhs[j];
        for k in j + 1..p {
            sum -= a[k][j] * z[k];
        }
        z[j] = sum / a[j][j];
    }
    if z.iter().all(|v| v.is_finite()) {
        Some(z)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INF: f64 = f64::INFINITY;

    fn cadquery(n: usize) -> Stopping {
        Stopping { ftol_rel: 0.0, ftol_abs: 0.0, xtol_rel: 1e-9, xtol_abs: vec![1e-12; n], max_evaluations: 0 }
    }

    /// Runs with every evaluated point checked against the bounds.
    fn run(x0: &[f64], lower: &[f64], upper: &[f64], f: impl Fn(&[f64], &mut [f64]) -> f64, stop: &Stopping) -> (Outcome, Vec<(Vec<f64>, f64)>) {
        let mut trace = Vec::new();
        let out = minimize(
            x0,
            lower,
            upper,
            |x, g| {
                for i in 0..x.len() {
                    assert!(x[i] >= lower[i] && x[i] <= upper[i], "evaluated outside the bounds: {x:?}");
                }
                let mut scratch = vec![0.0; x.len()];
                let v = match g {
                    Some(g) => f(x, g),
                    None => f(x, &mut scratch),
                };
                trace.push((x.to_vec(), v));
                v
            },
            stop,
        );
        (out, trace)
    }

    fn rosenbrock(x: &[f64], g: &mut [f64]) -> f64 {
        let (a, b) = (x[1] - x[0] * x[0], 1.0 - x[0]);
        g[0] = -400.0 * x[0] * a - 2.0 * b;
        g[1] = 200.0 * a;
        100.0 * a * a + b * b
    }

    fn close(a: &[f64], b: &[f64], tol: f64) {
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() <= tol, "{a:?} vs {b:?}");
        }
    }

    fn best_of(trace: &[(Vec<f64>, f64)]) -> (Vec<f64>, f64) {
        let mut best = (trace[0].0.clone(), trace[0].1);
        for (x, v) in trace {
            if *v < best.1 {
                best = (x.clone(), *v);
            }
        }
        best
    }

    #[test]
    fn rosenbrock_unbounded() {
        let (out, trace) = run(&[-1.2, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[1.0, 1.0], 1e-8);
        assert!(out.minimum < 1e-12);
        assert_eq!(out.evaluations, trace.len());
    }

    #[test]
    fn rosenbrock_bounded() {
        let (out, _) = run(&[-1.2, 1.0], &[-INF; 2], &[0.5, INF], rosenbrock, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[0.5, 0.25], 1e-8);
        assert!(out.x[0] <= 0.5);
        assert!((out.minimum - 0.25).abs() < 1e-12);
    }

    #[test]
    fn rosenbrock_evaluation_limit() {
        let mut stop = cadquery(2);
        stop.max_evaluations = 10;
        let (out, trace) = run(&[-1.2, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &stop);
        assert_eq!(out.result, Result::MaxevalReached);
        assert_eq!(out.evaluations, 10);
        assert_eq!(trace.len(), 10);
        let (bx, bv) = best_of(&trace);
        assert!(bv <= trace[0].1);
        assert_eq!(out.x, bx);
        assert_eq!(out.minimum, bv);
    }

    #[test]
    fn one_evaluation_allowed() {
        let mut stop = cadquery(2);
        stop.max_evaluations = 1;
        let (out, _) = run(&[-1.2, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &stop);
        assert_eq!((out.result, out.evaluations), (Result::MaxevalReached, 1));
        assert_eq!(out.x, vec![-1.2, 1.0]);
    }

    /// A chain of squares, (x₀ + 3)² + Σ (x_{i+1} − x_i − 1)².
    fn chain(x: &[f64], g: &mut [f64]) -> f64 {
        let n = x.len();
        g.iter_mut().for_each(|v| *v = 0.0);
        let mut f = 0.0;
        let r0 = x[0] + 3.0;
        f += r0 * r0;
        g[0] += 2.0 * r0;
        for i in 0..n - 1 {
            let r = x[i + 1] - x[i] - 1.0;
            f += r * r;
            g[i + 1] += 2.0 * r;
            g[i] -= 2.0 * r;
        }
        f
    }

    #[test]
    fn chain_with_an_active_lower_bound() {
        // Unconstrained optimum x_i = i − 3 (x₂ = −1); x₂ ≥ 0 makes it a bound optimum, the
        // variables after x₂ at x_i = i − 2.
        let n = 7;
        let mut lower = vec![-INF; n];
        lower[2] = 0.0;
        let x0 = vec![0.5; n];
        let (out, _) = run(&x0, &lower, &[INF; 7], chain, &cadquery(n));
        assert_eq!(out.result, Result::XtolReached);
        assert_eq!(out.x[2], 0.0);
        // The first two: least squares of x₀ = −3, x₁ − x₀ = 1, −x₁ = 1 → x₀ = −8/3, x₁ = −4/3,
        // each residual 1/3.
        let expected = [-8.0 / 3.0, -4.0 / 3.0, 0.0, 1.0, 2.0, 3.0, 4.0];
        close(&out.x, &expected, 1e-8);
        assert!((out.minimum - 3.0 * (1.0f64 / 3.0).powi(2)).abs() < 1e-12);
    }

    #[test]
    fn started_on_a_bound() {
        // f = (x₀ + 1)² + (x₁ − 1.5)² + x₀x₁ with x₀ ≥ 0, start (0, 3): x₀ stays at 0.
        let f = |x: &[f64], g: &mut [f64]| {
            g[0] = 2.0 * (x[0] + 1.0) + x[1];
            g[1] = 2.0 * (x[1] - 1.5) + x[0];
            (x[0] + 1.0).powi(2) + (x[1] - 1.5).powi(2) + x[0] * x[1]
        };
        let (out, _) = run(&[0.0, 3.0], &[0.0, -INF], &[INF; 2], f, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[0.0, 1.5], 1e-8);
    }

    #[test]
    fn against_a_bound() {
        // f = (x₀ + 1)² + (x₁ − 2)², x₀ ≥ 0: minimum (0, 2), value 1; ends with no descent left.
        let f = |x: &[f64], g: &mut [f64]| {
            g[0] = 2.0 * (x[0] + 1.0);
            g[1] = 2.0 * (x[1] - 2.0);
            (x[0] + 1.0).powi(2) + (x[1] - 2.0).powi(2)
        };
        let (out, _) = run(&[1.0, 0.0], &[0.0, -INF], &[INF; 2], f, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        assert_eq!(out.x[0], 0.0);
        close(&out.x, &[0.0, 2.0], 1e-12);
        assert!((out.minimum - 1.0).abs() < 1e-15);
        assert!(out.evaluations <= 6, "{}", out.evaluations);
    }

    #[test]
    fn started_at_the_minimum() {
        let (out, _) = run(&[1.0, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &cadquery(2));
        assert_eq!((out.result, out.evaluations), (Result::XtolReached, 1));
        assert_eq!(out.x, vec![1.0, 1.0]);
    }

    #[test]
    fn underdetermined() {
        // (x₀ + x₁ − 1)² + (x₂x₃ − 2)²: a family of minima.
        let f = |x: &[f64], g: &mut [f64]| {
            let a = x[0] + x[1] - 1.0;
            let b = x[2] * x[3] - 2.0;
            g[0] = 2.0 * a;
            g[1] = 2.0 * a;
            g[2] = 2.0 * b * x[3];
            g[3] = 2.0 * b * x[2];
            a * a + b * b
        };
        let (out, _) = run(&[3.0, -1.0, 0.5, 1.0], &[-INF; 4], &[INF; 4], f, &cadquery(4));
        assert_eq!(out.result, Result::XtolReached);
        assert!(out.minimum < 1e-20, "{}", out.minimum);
        assert!((out.x[0] + out.x[1] - 1.0).abs() < 1e-10);
        assert!((out.x[2] * out.x[3] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn infinite_values() {
        // f = 5x + 1/x, +∞ for x ≤ 0: minimum at 1/√5; the first full step from 2 lands at
        // −2.75.
        let f = |x: &[f64], g: &mut [f64]| {
            if x[0] <= 0.0 {
                return INF;
            }
            g[0] = 5.0 - 1.0 / (x[0] * x[0]);
            5.0 * x[0] + 1.0 / x[0]
        };
        let (out, trace) = run(&[2.0], &[-INF], &[INF], f, &cadquery(1));
        assert!(trace.iter().any(|(_, v)| v.is_infinite()));
        assert_eq!(out.result, Result::XtolReached);
        assert!(out.minimum.is_finite());
        close(&out.x, &[1.0 / 5f64.sqrt()], 1e-8);
    }

    #[test]
    fn mostly_on_bounds() {
        // Σ cᵢ xᵢ + (x₂ − 0.1)² over [−0.6, 0.6]⁸: all but x₂ at a bound.
        let c = [1.0, -2.0, 0.0, 0.5, -0.5, 3.0, -1.0, 0.25];
        let f = |x: &[f64], g: &mut [f64]| {
            let mut v = (x[2] - 0.1).powi(2);
            for i in 0..8 {
                v += c[i] * x[i];
                g[i] = c[i];
            }
            g[2] += 2.0 * (x[2] - 0.1);
            v
        };
        let (out, _) = run(&[0.0; 8], &[-0.6; 8], &[0.6; 8], f, &cadquery(8));
        assert_eq!(out.result, Result::XtolReached);
        for i in 0..8 {
            if i != 2 {
                let bound = if c[i] > 0.0 { -0.6 } else { 0.6 };
                assert!((out.x[i] - bound).abs() < 4.0 * f64::EPSILON, "{:?}", out.x);
            }
        }
        assert!((out.x[2] - 0.1).abs() < 1e-8);
        let expected: f64 = c.iter().enumerate().filter(|(i, _)| *i != 2).map(|(_, ci)| -0.6 * ci.abs()).sum();
        assert!((out.minimum - expected).abs() < 1e-12);
    }

    #[test]
    fn start_outside_the_bounds() {
        let mut called = false;
        let out = minimize(&[2.0], &[0.0], &[1.0], |_, _| { called = true; 0.0 }, &cadquery(1));
        assert_eq!((out.result, out.evaluations, out.minimum), (Result::InvalidArgs, 0, INF));
        assert_eq!(out.x, vec![2.0]);
        assert!(!called);
        let out = minimize(&[0.5], &[1.0], &[0.0], |_, _| 0.0, &cadquery(1));
        assert_eq!((out.result, out.evaluations), (Result::InvalidArgs, 0));
    }

    #[test]
    fn no_variables() {
        let mut asked_gradient = false;
        let out = minimize(&[], &[], &[], |_, g| { asked_gradient = g.is_some(); 7.0 }, &cadquery(0));
        assert_eq!((out.result, out.evaluations, out.minimum), (Result::Success, 1, 7.0));
        assert!(out.x.is_empty() && !asked_gradient);
    }

    /// Hock–Schittkowski problems with bounds only (starts moved into the bounds where the
    /// collection's lie outside).
    #[test]
    fn hock_schittkowski_1() {
        let (out, _) = run(&[-2.0, 1.0], &[-INF, -1.5], &[INF; 2], rosenbrock, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[1.0, 1.0], 1e-6);
        assert!(out.minimum < 1e-12);
    }

    #[test]
    fn hock_schittkowski_2() {
        // From the collection's (−2, 1) moved onto the bound, the path ends at the other local
        // minimum on x₁ = 1.5, x₀ ≈ −1.2210; from (2, 2) at the collection's optimum.
        let (out, _) = run(&[-2.0, 1.5], &[-INF, 1.5], &[INF; 2], rosenbrock, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        assert!(out.x[1] == 1.5 && (out.x[0] + 1.2210262421).abs() < 1e-8, "{:?}", out.x);
        let (out, _) = run(&[2.0, 2.0], &[-INF, 1.5], &[INF; 2], rosenbrock, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        // On x₁ = 1.5: d/dx₀ = −400x₀(1.5 − x₀²) − 2(1 − x₀) = 0, by Newton from 1.2.
        let mut t: f64 = 1.2;
        for _ in 0..50 {
            let h = -400.0 * t * (1.5 - t * t) - 2.0 * (1.0 - t);
            let dh = -600.0 + 1200.0 * t * t + 2.0;
            t -= h / dh;
        }
        close(&out.x, &[t, 1.5], 1e-6);
        assert!((out.minimum - 0.0504261879).abs() < 1e-6);
    }

    #[test]
    fn hock_schittkowski_3() {
        let f = |x: &[f64], g: &mut [f64]| {
            let a = x[1] - x[0];
            g[0] = -2e-5 * a;
            g[1] = 1.0 + 2e-5 * a;
            x[1] + 1e-5 * a * a
        };
        let (out, _) = run(&[10.0, 1.0], &[-INF, 0.0], &[INF; 2], f, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[0.0, 0.0], 1e-6);
        assert!(out.minimum.abs() < 1e-6);
    }

    #[test]
    fn hock_schittkowski_4() {
        let f = |x: &[f64], g: &mut [f64]| {
            g[0] = (x[0] + 1.0).powi(2);
            g[1] = 1.0;
            (x[0] + 1.0).powi(3) / 3.0 + x[1]
        };
        let (out, _) = run(&[1.125, 0.125], &[1.0, 0.0], &[INF; 2], f, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[1.0, 0.0], 1e-12);
        assert!((out.minimum - 8.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn hock_schittkowski_5() {
        let f = |x: &[f64], g: &mut [f64]| {
            let c = (x[0] + x[1]).cos();
            g[0] = c + 2.0 * (x[0] - x[1]) - 1.5;
            g[1] = c - 2.0 * (x[0] - x[1]) + 2.5;
            (x[0] + x[1]).sin() + (x[0] - x[1]).powi(2) - 1.5 * x[0] + 2.5 * x[1] + 1.0
        };
        let (out, _) = run(&[0.0, 0.0], &[-1.5, -3.0], &[4.0, 3.0], f, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        let pi3 = std::f64::consts::PI / 3.0;
        close(&out.x, &[0.5 - pi3, -0.5 - pi3], 1e-6);
        assert!((out.minimum - (-(3.0f64).sqrt() / 2.0 - pi3)).abs() < 1e-10);
    }

    #[test]
    fn hock_schittkowski_38() {
        let f = |x: &[f64], g: &mut [f64]| {
            let a = x[1] - x[0] * x[0];
            let b = x[3] - x[2] * x[2];
            g[0] = -400.0 * x[0] * a - 2.0 * (1.0 - x[0]);
            g[1] = 200.0 * a + 20.2 * (x[1] - 1.0) + 19.8 * (x[3] - 1.0);
            g[2] = -360.0 * x[2] * b - 2.0 * (1.0 - x[2]);
            g[3] = 180.0 * b + 20.2 * (x[3] - 1.0) + 19.8 * (x[1] - 1.0);
            100.0 * a * a + (1.0 - x[0]).powi(2) + 90.0 * b * b + (1.0 - x[2]).powi(2)
                + 10.1 * ((x[1] - 1.0).powi(2) + (x[3] - 1.0).powi(2))
                + 19.8 * (x[1] - 1.0) * (x[3] - 1.0)
        };
        let (out, _) = run(&[-3.0, -1.0, -3.0, -1.0], &[-10.0; 4], &[10.0; 4], f, &cadquery(4));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[1.0; 4], 1e-6);
        assert!(out.minimum < 1e-12);
    }

    #[test]
    fn hock_schittkowski_45() {
        let f = |x: &[f64], g: &mut [f64]| {
            let p: f64 = x.iter().product();
            for i in 0..5 {
                g[i] = -x.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, v)| v).product::<f64>() / 120.0;
            }
            2.0 - p / 120.0
        };
        let upper = [1.0, 2.0, 3.0, 4.0, 5.0];
        let (out, _) = run(&[0.5, 1.0, 1.5, 2.0, 2.5], &[0.0; 5], &upper, f, &cadquery(5));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &upper, 1e-12);
        assert!((out.minimum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hock_schittkowski_110() {
        let f = |x: &[f64], g: &mut [f64]| {
            let p: f64 = x.iter().product();
            let q = p.powf(0.2);
            let mut v = -q;
            for i in 0..10 {
                let (a, b) = ((x[i] - 2.0).ln(), (10.0 - x[i]).ln());
                v += a * a + b * b;
                g[i] = 2.0 * a / (x[i] - 2.0) - 2.0 * b / (10.0 - x[i]) - 0.2 * q / x[i];
            }
            v
        };
        let (out, _) = run(&[9.0; 10], &[2.001; 10], &[9.999; 10], f, &cadquery(10));
        assert_eq!(out.result, Result::XtolReached);
        // The symmetric stationary point by bisection (the collection prints 9.35025655,
        // 1e-5 off; its value agrees).
        let h = |t: f64| 2.0 * (t - 2.0).ln() / (t - 2.0) - 2.0 * (10.0 - t).ln() / (10.0 - t) - 0.2 * t;
        let (mut a, mut b) = (9.0, 9.9);
        for _ in 0..100 {
            let m = 0.5 * (a + b);
            if h(a) * h(m) <= 0.0 { b = m } else { a = m }
        }
        close(&out.x, &[0.5 * (a + b); 10], 1e-6);
        assert!((out.minimum + 45.77846971).abs() < 1e-6);
    }

    #[test]
    fn fixed_variables() {
        // x₁ fixed at 0.3; the others of a coupled quadratic optimised.
        let f = |x: &[f64], g: &mut [f64]| {
            g[0] = 2.0 * (x[0] - 1.0) + x[1];
            g[1] = x[0] + 2.0 * (x[1] - 2.0) + x[2];
            g[2] = x[1] + 2.0 * (x[2] + 1.0);
            (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2) + (x[2] + 1.0).powi(2) + x[0] * x[1] + x[1] * x[2]
        };
        let (out, trace) = run(&[0.0, 0.3, 0.0], &[-INF, 0.3, -INF], &[INF, 0.3, INF], f, &cadquery(3));
        assert_eq!(out.result, Result::XtolReached);
        assert!(trace.iter().all(|(x, _)| x[1] == 0.3));
        close(&out.x, &[1.0 - 0.15, 0.3, -1.0 - 0.15], 1e-8);
    }

    #[test]
    fn mixed_bounds_and_a_start_pushed_out() {
        // (x₀ − 5)² + (x₁ + 5)² + (x₂ − 0.5)² in [−1, 2] × [−2, ∞) × (−∞, 1], started at the
        // upper bound of x₀ with the gradient pushing out of it.
        let f = |x: &[f64], g: &mut [f64]| {
            g[0] = 2.0 * (x[0] - 5.0);
            g[1] = 2.0 * (x[1] + 5.0);
            g[2] = 2.0 * (x[2] - 0.5);
            (x[0] - 5.0).powi(2) + (x[1] + 5.0).powi(2) + (x[2] - 0.5).powi(2)
        };
        let (out, _) = run(&[2.0, 0.0, -3.0], &[-1.0, -2.0, -INF], &[2.0, INF, 1.0], f, &cadquery(3));
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[2.0, -2.0, 0.5], 1e-10);
        assert!(out.x[0] <= 2.0 && out.x[1] >= -2.0);
    }

    #[test]
    fn a_quadratic_in_few_iterations() {
        // ½ xᵀAx − bᵀx with A positive definite (n = 5): quasi-Newton steps end it in a few
        // more evaluations than n (the line search is not exact).
        let a = [
            [4.0, 1.0, 0.0, 0.5, 0.0],
            [1.0, 3.0, 0.2, 0.0, 0.0],
            [0.0, 0.2, 2.0, 0.3, 0.1],
            [0.5, 0.0, 0.3, 5.0, 1.0],
            [0.0, 0.0, 0.1, 1.0, 1.5],
        ];
        let b = [1.0, -2.0, 0.5, 3.0, -1.0];
        let f = |x: &[f64], g: &mut [f64]| {
            let mut v = 0.0;
            for i in 0..5 {
                let ax: f64 = (0..5).map(|j| a[i][j] * x[j]).sum();
                g[i] = ax - b[i];
                v += 0.5 * x[i] * ax - b[i] * x[i];
            }
            v
        };
        let (out, _) = run(&[0.0; 5], &[-INF; 5], &[INF; 5], f, &cadquery(5));
        assert_eq!(out.result, Result::XtolReached);
        assert!(out.evaluations <= 3 * 5, "{}", out.evaluations);
        let mut g = [0.0; 5];
        f(&out.x, &mut g);
        assert!(g.iter().all(|v| v.abs() < 1e-10), "{g:?}");
    }

    #[test]
    fn a_noisy_gradient() {
        // A sum of squares with forward differences of step 1e-10, as the sketch solver's.
        let cost = |x: &[f64]| (x[0] + x[1] - 3.0).powi(2) + (x[0] * x[1] - 2.0).powi(2) + (x[2] - x[0]).powi(2);
        let f = |x: &[f64], g: &mut [f64]| {
            let v = cost(x);
            let mut y = x.to_vec();
            for i in 0..3 {
                y[i] += 1e-10;
                g[i] = (cost(&y) - v) / 1e-10;
                y[i] = x[i];
            }
            v
        };
        let (out, _) = run(&[0.2, 3.0, 0.0], &[-INF; 3], &[INF; 3], f, &cadquery(3));
        assert_eq!(out.result, Result::XtolReached);
        assert!(out.minimum < 1e-16, "{}", out.minimum);
    }

    #[test]
    fn deterministic() {
        let (a, ta) = run(&[-1.2, 1.0], &[-INF, -INF], &[0.8, INF], rosenbrock, &cadquery(2));
        let (b, tb) = run(&[-1.2, 1.0], &[-INF, -INF], &[0.8, INF], rosenbrock, &cadquery(2));
        assert_eq!(ta.len(), tb.len());
        for ((xa, va), (xb, vb)) in ta.iter().zip(&tb) {
            assert_eq!(xa.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), xb.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            assert_eq!(va.to_bits(), vb.to_bits());
        }
        assert_eq!((a.result, a.evaluations, a.minimum.to_bits()), (b.result, b.evaluations, b.minimum.to_bits()));
    }

    #[test]
    fn relative_tolerance_only() {
        let stop = Stopping { ftol_rel: 0.0, ftol_abs: 0.0, xtol_rel: 1e-10, xtol_abs: Vec::new(), max_evaluations: 0 };
        let (out, _) = run(&[-1.2, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &stop);
        assert_eq!(out.result, Result::XtolReached);
        close(&out.x, &[1.0, 1.0], 1e-8);
    }

    #[test]
    fn absolute_f_tolerance() {
        let stop = Stopping { ftol_rel: 0.0, ftol_abs: 1e-6, xtol_rel: 0.0, xtol_abs: Vec::new(), max_evaluations: 0 };
        let (out, _) = run(&[-1.2, 1.0], &[-INF; 2], &[INF; 2], rosenbrock, &stop);
        assert_eq!(out.result, Result::FtolReached);
        assert!(out.minimum < 1e-4);
    }

    #[test]
    fn flat_and_singular_problems_end() {
        // Constant: no descent at the start.
        let (out, _) = run(&[1.0, 2.0], &[-INF; 2], &[INF; 2], |_, g| { g.iter_mut().for_each(|v| *v = 0.0); 3.0 }, &cadquery(2));
        assert_eq!((out.result, out.evaluations), (Result::XtolReached, 1));
        // Linear, falling towards a bound on one variable and unbounded in no other.
        let (out, _) = run(&[1.0, 2.0], &[0.0, -INF], &[INF; 2], |x, g| { g[0] = 1.0; g[1] = 0.0; x[0] }, &cadquery(2));
        assert_eq!(out.result, Result::XtolReached);
        assert_eq!(out.x, vec![0.0, 2.0]);
        // A gradient that is NaN: RoundoffLimited, not a hang or panic.
        let (out, _) = run(&[1.0], &[-INF], &[INF], |x, g| { g[0] = f64::NAN; x[0] * x[0] }, &cadquery(1));
        assert_eq!(out.result, Result::RoundoffLimited);
        // All tolerances off and no progress possible: RoundoffLimited.
        let stop = Stopping { ftol_rel: 0.0, ftol_abs: 0.0, xtol_rel: 0.0, xtol_abs: Vec::new(), max_evaluations: 0 };
        let (out, _) = run(&[1.0], &[-INF], &[INF], |_, g| { g[0] = 0.0; 1.0 }, &stop);
        assert_eq!(out.result, Result::RoundoffLimited);
    }
}
