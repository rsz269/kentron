//! Levenberg–Marquardt least squares: the variables that make the sum of squares of some
//! residuals least, from a start, within bounds on each variable.
//!
//! Each step solves the normal equations of the residuals' linearization, damped by μ:
//! (JᵀJ + μI) h = −Jᵀr. μ shrinks after a step that gains about what the linearization
//! foretold and grows after one that does not (Nielsen's rule), so the method runs as
//! Gauss–Newton near a solution (quadratically, where the residuals vanish there) and as
//! steepest descent far from one. Steps leaving the bounds are cut back onto them.

/// When to stop.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Settings {
    /// Most steps tried (taken or not).
    pub max_steps: usize,
    /// Stop once a step moves the variables less than `step_relative` times their size plus
    /// `step_absolute`.
    pub step_relative: f64,
    pub step_absolute: f64,
    /// Stop once the sum of squares is at most this.
    pub cost: f64,
    /// Stop once no component of the gradient (of half the sum of squares) exceeds this.
    pub gradient: f64,
}

/// Why the minimization stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stop {
    /// The sum of squares fell to `Settings::cost`.
    Cost,
    /// The gradient vanished (within `Settings::gradient`).
    Gradient,
    /// Steps became too small to move the variables (within `Settings::step_*`).
    Step,
    /// `Settings::max_steps` steps were tried.
    Steps,
}

#[derive(Debug, Clone)]
pub(crate) struct Outcome {
    pub variables: Vec<f64>,
    /// The sum of squares of the residuals at `variables`.
    pub cost: f64,
    /// The steps tried and the residual evaluations made.
    pub steps: usize,
    pub evaluations: usize,
    pub stop: Stop,
}

/// The residuals at some variables, and the Jacobian's rows: each row the nonzero partial
/// derivatives of one residual, by variable index.
pub(crate) struct Evaluation {
    pub residuals: Vec<f64>,
    pub rows: Vec<Vec<(usize, f64)>>,
}

/// The variables making the sum of squares of `evaluate`'s residuals least, from `start`,
/// each kept within `lower[i]..=upper[i]`.
pub(crate) fn minimize(start: &[f64], lower: &[f64], upper: &[f64], mut evaluate: impl FnMut(&[f64]) -> Evaluation, settings: Settings) -> Outcome {
    let n = start.len();
    let clamp = |x: &mut [f64]| {
        for (i, v) in x.iter_mut().enumerate() {
            *v = v.clamp(lower[i], upper[i]);
        }
    };
    let mut x = start.to_vec();
    clamp(&mut x);
    let mut at = evaluate(&x);
    let mut evaluations = 1;
    let half_cost = |e: &Evaluation| 0.5 * e.residuals.iter().map(|r| r * r).sum::<f64>();
    let mut cost = half_cost(&at);
    let (mut normal, mut gradient) = normal_equations(&at, n);
    let mut damping = 1e-3 * (0..n).map(|i| normal[i * n + i]).fold(0.0, f64::max).max(1e-12);
    let mut growth = 2.0;
    let mut steps = 0;
    let stop = loop {
        if 2.0 * cost <= settings.cost {
            break Stop::Cost;
        }
        // (The gradient within the bounds: what pushes a variable on a bound past it counts not.)
        let held = |i: usize| (x[i] <= lower[i] && gradient[i] > 0.0) || (x[i] >= upper[i] && gradient[i] < 0.0);
        if (0..n).all(|i| held(i) || gradient[i].abs() <= settings.gradient) {
            break Stop::Gradient;
        }
        if steps >= settings.max_steps {
            break Stop::Steps;
        }
        steps += 1;
        // Variables on a bound that descent would push past it are held there (an active
        // set): the step is solved for the others.
        let free: Vec<usize> = (0..n).filter(|&i| !held(i)).collect();
        let f = free.len();
        let mut damped = vec![0.0; f * f];
        for (a, &i) in free.iter().enumerate() {
            for (b, &j) in free.iter().enumerate() {
                damped[a * f + b] = normal[i * n + j];
            }
            damped[a * f + a] += damping;
        }
        let negative: Vec<f64> = free.iter().map(|&i| -gradient[i]).collect();
        let Some(reduced) = cholesky_solve(&mut damped, &negative, f) else {
            damping *= growth;
            growth *= 2.0;
            continue;
        };
        let mut h = vec![0.0; n];
        for (a, &i) in free.iter().enumerate() {
            h[i] = reduced[a];
        }
        // A step crossing a bound is shortened to stop short of it (a large step from far
        // away would otherwise put variables on bounds where they need not end), unless that
        // leaves little of it: then the variables crossing stop on their bounds.
        let mut fraction: f64 = 1.0;
        for i in 0..n {
            let target = x[i] + h[i];
            if target < lower[i] && x[i] > lower[i] {
                fraction = fraction.min(0.995 * (lower[i] - x[i]) / h[i]);
            } else if target > upper[i] && x[i] < upper[i] {
                fraction = fraction.min(0.995 * (upper[i] - x[i]) / h[i]);
            }
        }
        if fraction < 0.1 {
            fraction = 1.0;
        }
        let mut trial: Vec<f64> = x.iter().zip(&h).map(|(a, b)| a + b * fraction).collect();
        clamp(&mut trial);
        let taken: Vec<f64> = trial.iter().zip(&x).map(|(a, b)| a - b).collect();
        let size = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        let length = taken.iter().map(|v| v * v).sum::<f64>().sqrt();
        if length <= settings.step_relative * size + settings.step_absolute {
            break Stop::Step;
        }
        // What the linearization foretells the step gains: −gᵀh − ½ hᵀ(JᵀJ)h.
        let mut quadratic = 0.0;
        for i in 0..n {
            let row: f64 = (0..n).map(|j| normal[i * n + j] * taken[j]).sum();
            quadratic += taken[i] * row;
        }
        let foretold = -gradient.iter().zip(&taken).map(|(g, h)| g * h).sum::<f64>() - 0.5 * quadratic;
        let next = evaluate(&trial);
        evaluations += 1;
        let next_cost = half_cost(&next);
        let gain = cost - next_cost;
        if foretold > 0.0 && gain > 0.0 {
            let ratio = gain / foretold;
            damping *= (1.0 / 3.0_f64).max(1.0 - (2.0 * ratio - 1.0).powi(3));
            growth = 2.0;
            x = trial;
            at = next;
            cost = next_cost;
            (normal, gradient) = normal_equations(&at, n);
        } else {
            damping *= growth;
            growth *= 2.0;
            if !damping.is_finite() || damping > 1e300 {
                break Stop::Step;
            }
        }
    };
    Outcome { variables: x, cost: 2.0 * cost, steps, evaluations, stop }
}

/// JᵀJ (dense, row-major) and Jᵀr.
fn normal_equations(e: &Evaluation, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut normal = vec![0.0; n * n];
    let mut gradient = vec![0.0; n];
    for (row, &r) in e.rows.iter().zip(&e.residuals) {
        for &(i, a) in row {
            gradient[i] += a * r;
            for &(j, b) in row {
                normal[i * n + j] += a * b;
            }
        }
    }
    (normal, gradient)
}

/// Solves `a x = b` for a symmetric positive definite `a` (row-major, overwritten by its
/// Cholesky factor); `None` where `a` is not positive definite.
fn cholesky_solve(a: &mut [f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j {
            d -= a[j * n + k] * a[j * n + k];
        }
        if !(d > 0.0) {
            return None;
        }
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = s / d;
        }
    }
    let mut y = b.to_vec();
    for i in 0..n {
        for k in 0..i {
            y[i] -= a[i * n + k] * y[k];
        }
        y[i] /= a[i * n + i];
    }
    for i in (0..n).rev() {
        for k in i + 1..n {
            y[i] -= a[k * n + i] * y[k];
        }
        y[i] /= a[i * n + i];
    }
    Some(y)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings() -> Settings {
        Settings { max_steps: 1000, step_relative: 1e-15, step_absolute: 1e-15, cost: 1e-30, gradient: 0.0 }
    }

    /// Rosenbrock's valley as residuals (10(y − x²), 1 − x), from (−1.2, 1): its least at (1, 1).
    #[test]
    fn rosenbrock_valley() {
        let out = minimize(&[-1.2, 1.0], &[f64::NEG_INFINITY; 2], &[f64::INFINITY; 2], |x| Evaluation {
            residuals: vec![10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]],
            rows: vec![vec![(0, -20.0 * x[0]), (1, 10.0)], vec![(0, -1.0)]],
        }, settings());
        assert!((out.variables[0] - 1.0).abs() < 1e-10 && (out.variables[1] - 1.0).abs() < 1e-10, "{out:?}");
    }

    /// A residual x + 1 with x kept at or above 0: the least within the bound is at 0.
    #[test]
    fn a_bound_holds_the_least_on_it() {
        let out = minimize(&[2.0], &[0.0], &[f64::INFINITY], |x| Evaluation { residuals: vec![x[0] + 1.0], rows: vec![vec![(0, 1.0)]] }, settings());
        assert_eq!(out.variables[0], 0.0, "{out:?}");
        assert!((out.cost - 1.0).abs() < 1e-12);
    }
}
