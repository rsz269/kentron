//! Newton's method for an unconstrained minimum, step by step as IPOPT takes it on a problem
//! with neither constraints nor bounds (as CadQuery's assembly solver poses its problem), so
//! that where the minimum is not unique (a part free to turn about an axis) it ends where
//! IPOPT does.
//!
//! Each iteration solves (H + δI) d = −g. δ is 0 if H is positive definite;
//! else IPOPT's inertia correction: the first time 1e-4, multiplied by 100 until H + δI is
//! positive definite; afterwards a third of the last δ, multiplied by 8 until it is. The step
//! is d, halved until the objective falls by at least 1e-8 of what the gradient foretells
//! (Armijo). It stops once no gradient component exceeds `tolerance`, or, IPOPT's "acceptable"
//! level, once none exceeds `acceptable` and the objective changed by at most
//! `acceptable_change` (relative to it, or to 1 if less) over the last iteration.

/// When to stop (CadQuery's IPOPT options: tol 1e-14, acceptable_tol 1e-6 by default,
/// acceptable_obj_change_tol 1e-12, acceptable_iter 1).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Settings {
    pub tolerance: f64,
    pub acceptable: f64,
    pub acceptable_change: f64,
    pub max_iterations: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stop {
    /// The gradient vanished within `tolerance` (IPOPT's "Solve_Succeeded").
    Optimal,
    /// Within `acceptable` and no longer changing ("Solved_To_Acceptable_Level").
    Acceptable,
    /// `max_iterations` taken.
    Iterations,
    /// No step along the direction lowered the objective.
    LineSearch,
}

#[derive(Debug, Clone)]
pub(crate) struct Outcome {
    pub variables: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub stop: Stop,
    /// The regularization δ of each iteration (for comparing with IPOPT's log).
    pub regularizations: Vec<f64>,
}

/// The objective, its gradient and its Hessian (dense, row-major) at some variables.
pub(crate) struct Evaluation {
    pub objective: f64,
    pub gradient: Vec<f64>,
    pub hessian: Vec<f64>,
}

pub(crate) fn minimize(start: &[f64], mut evaluate: impl FnMut(&[f64]) -> Evaluation, settings: Settings) -> Outcome {
    let n = start.len();
    let mut x = start.to_vec();
    let mut at = evaluate(&x);
    let mut evaluations = 1;
    let mut last_delta = 0.0_f64;
    let mut previous: Option<f64> = None;
    let mut regularizations = Vec::new();
    let mut iterations = 0;
    let stop = loop {
        let error = at.gradient.iter().fold(0.0_f64, |m, g| m.max(g.abs()));
        if error <= settings.tolerance {
            break Stop::Optimal;
        }
        if let Some(p) = previous {
            if error <= settings.acceptable && (at.objective - p).abs() / at.objective.abs().max(1.0) <= settings.acceptable_change {
                break Stop::Acceptable;
            }
        }
        if iterations >= settings.max_iterations {
            break Stop::Iterations;
        }
        // The direction, with IPOPT's inertia correction.
        let negative: Vec<f64> = at.gradient.iter().map(|g| -g).collect();
        let solve = |delta: f64| -> Option<Vec<f64>> {
            let mut a = at.hessian.clone();
            for i in 0..n {
                a[i * n + i] += delta;
            }
            cholesky_solve(&mut a, &negative, n)
        };
        let mut delta = 0.0;
        let mut direction = solve(0.0);
        if direction.is_none() {
            let first = last_delta == 0.0;
            delta = if first { 1e-4 } else { (last_delta / 3.0).max(1e-20) };
            loop {
                direction = solve(delta);
                if direction.is_some() || delta > 1e40 {
                    break;
                }
                delta *= if first { 100.0 } else { 8.0 };
            }
            last_delta = delta;
        }
        regularizations.push(delta);
        let Some(d) = direction else { break Stop::LineSearch };
        // Armijo backtracking.
        let slope: f64 = at.gradient.iter().zip(&d).map(|(g, d)| g * d).sum();
        let mut alpha = 1.0;
        let mut accepted = None;
        for _ in 0..60 {
            let trial: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
            let next = evaluate(&trial);
            evaluations += 1;
            if next.objective <= at.objective + 1e-8 * alpha * slope {
                accepted = Some((trial, next));
                break;
            }
            alpha *= 0.5;
        }
        let Some((trial, next)) = accepted else { break Stop::LineSearch };
        previous = Some(at.objective);
        x = trial;
        at = next;
        iterations += 1;
    };
    Outcome { variables: x, objective: at.objective, iterations, evaluations, stop, regularizations }
}

/// Solves `a x = b` by Cholesky; `None` where `a` is not positive definite (a pivot not above
/// 0: IPOPT's inertia test; one barely above, as where a part is free to turn and only
/// CadQuery's penalty holds it, passes, as it does for IPOPT).
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

    /// Rosenbrock's function (1 − x)² + 100(y − x²)² from (−1.2, 1): its minimum (1, 1).
    #[test]
    fn rosenbrock() {
        let out = minimize(&[-1.2, 1.0], |v| {
            let (x, y) = (v[0], v[1]);
            Evaluation {
                objective: (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2),
                gradient: vec![-2.0 * (1.0 - x) - 400.0 * x * (y - x * x), 200.0 * (y - x * x)],
                hessian: vec![2.0 - 400.0 * y + 1200.0 * x * x, -400.0 * x, -400.0 * x, 200.0],
            }
        }, Settings { tolerance: 1e-12, acceptable: 0.0, acceptable_change: 0.0, max_iterations: 200 });
        assert!(out.stop == Stop::Optimal && (out.variables[0] - 1.0).abs() < 1e-9 && (out.variables[1] - 1.0).abs() < 1e-9, "{out:?}");
    }
}
