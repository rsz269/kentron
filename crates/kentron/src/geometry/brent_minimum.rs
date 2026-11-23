//! A local minimum of a function of one variable in a bracket, without derivatives: Brent's
//! method (parabolic interpolation through the three best points, safeguarded by
//! golden-section steps), R. P. Brent, *Algorithms for Minimization without Derivatives*
//! (1973), ch. 5.

/// The golden-section fraction `(3 − √5)/2`.
const GOLDEN: f64 = 0.381_966_011_250_105_1;

/// Brent's minimiser with a location tolerance `τ(x) = relative·|x| + absolute`.
pub(crate) struct BrentMinimum {
    relative: f64,
    absolute: f64,
    max_iterations: usize,
    done: bool,
    /// The best point and value of the last search, if it evaluated anything.
    best: Option<(f64, f64)>,
    iterations: usize,
}

impl BrentMinimum {
    pub(crate) fn new(relative_tolerance: f64, max_iterations: usize, absolute_tolerance: f64) -> Self {
        BrentMinimum { relative: relative_tolerance, absolute: absolute_tolerance, max_iterations, done: false, best: None, iterations: 0 }
    }

    /// Searches `[min(ax, cx), max(ax, cx)]` from `bx` (inside it) for a local minimum of
    /// `function` (`None` where it cannot be evaluated: the search then stops, not done).
    /// Done when the bracket around the best point `x` is within `2·τ(x)` of it on both sides.
    /// Evaluates only inside the bracket, at least `τ(x)` from the best point, at most
    /// `max_iterations + 1` times.
    pub(crate) fn perform(&mut self, function: &mut dyn FnMut(f64) -> Option<f64>, ax: f64, bx: f64, cx: f64) {
        self.done = false;
        self.best = None;
        self.iterations = 0;
        let (mut a, mut b) = (ax.min(cx), ax.max(cx));
        let start = bx.clamp(a, b);
        let Some(f_start) = function(start) else { return };
        // x the best point, w the second best, v the previous w.
        let (mut x, mut w, mut v) = (start, start, start);
        let (mut fx, mut fw, mut fv) = (f_start, f_start, f_start);
        self.best = Some((x, fx));
        // d the last step, e the one before.
        let (mut d, mut e): (f64, f64) = (0.0, 0.0);
        for iteration in 0..=self.max_iterations {
            let middle = 0.5 * (a + b);
            let tol1 = self.relative * x.abs() + self.absolute;
            let tol2 = 2.0 * tol1;
            if (x - middle).abs() <= tol2 - 0.5 * (b - a) {
                self.done = true;
                return;
            }
            if iteration == self.max_iterations {
                return;
            }
            let mut golden = true;
            if e.abs() > tol1 {
                // The parabola through (v, w, x).
                let r = (x - w) * (fx - fv);
                let mut q = (x - v) * (fx - fw);
                let mut p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if q > 0.0 {
                    p = -p;
                }
                q = q.abs();
                let previous = e;
                if !(p.abs() >= (0.5 * q * previous).abs() || p <= q * (a - x) || p >= q * (b - x)) {
                    e = d;
                    d = p / q;
                    let u = x + d;
                    if u - a < tol2 || b - u < tol2 {
                        d = tol1.copysign(middle - x);
                    }
                    golden = false;
                }
            }
            if golden {
                e = if x >= middle { a - x } else { b - x };
                d = GOLDEN * e;
            }
            let u = (if d.abs() >= tol1 { x + d } else { x + tol1.copysign(d) }).clamp(a, b);
            self.iterations += 1;
            let Some(fu) = function(u) else { return };
            if fu <= fx {
                if u >= x {
                    a = x;
                } else {
                    b = x;
                }
                (v, w, x) = (w, x, u);
                (fv, fw, fx) = (fw, fx, fu);
                self.best = Some((x, fx));
            } else {
                if u < x {
                    a = u;
                } else {
                    b = u;
                }
                if fu <= fw || w == x {
                    (v, w) = (w, u);
                    (fv, fw) = (fw, fu);
                } else if fu <= fv || v == x || v == w {
                    v = u;
                    fv = fu;
                }
            }
        }
    }

    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// The minimum's abscissa.
    pub(crate) fn location(&self) -> f64 {
        assert!(self.done, "no minimum located");
        self.best.expect("a minimum").0
    }

    /// The value at [`location`](Self::location).
    pub(crate) fn minimum(&self) -> f64 {
        assert!(self.done, "no minimum located");
        self.best.expect("a minimum").1
    }

    /// The best point and value of the last search, done or not.
    pub(crate) fn best(&self) -> Option<(f64, f64)> {
        self.best
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parabola_minimum() {
        let mut m = BrentMinimum::new(1e-10, 100, 1e-12);
        m.perform(&mut |x| Some((x - 0.3) * (x - 0.3) + 2.0), -1.0, 0.0, 1.0);
        assert!(m.is_done());
        assert!((m.location() - 0.3).abs() < 1e-7);
        assert_eq!(m.minimum(), (m.location() - 0.3).powi(2) + 2.0);
        assert!(m.iterations < 50, "{}", m.iterations);
    }

    #[test]
    fn local_minima_within_the_tolerance() {
        // cos 3x + 0.2x: local minima near 1.0250 and −1.0694.
        let f = |x: f64| (3.0 * x).cos() + 0.2 * x;
        for (a, b, c, rho, alpha, expected) in [(0.0f64, 0.5f64, 2.0f64, 1e-8f64, 1e-10f64, 1.0250f64), (0.5, 1.0, 1.5, 1e-3, 1e-6, 1.0250), (-2.0, -1.0, 0.0, 1e-6, 1e-9, -1.0694), (0.0, 2.0, 2.0, 1e-10, 1e-12, 1.0250)] {
            let mut calls = Vec::new();
            let mut m = BrentMinimum::new(rho, 100, alpha);
            m.perform(&mut |x| {
                calls.push(x);
                (x >= a.min(c) && x <= a.max(c)).then(|| f(x))
            }, a, b, c);
            assert!(m.is_done());
            let x = m.location();
            // The true minimum, from the derivative −3 sin 3x + 0.2 = 0.
            let mut t: f64 = expected;
            for _ in 0..50 {
                t -= (-3.0 * (3.0 * t).sin() + 0.2) / (-9.0 * (3.0 * t).cos());
            }
            // (Values locate a minimum to about √ε at best: f − f* ~ f″δ²/2 is lost in rounding below.)
            assert!((x - t).abs() <= 2.0 * (rho * x.abs() + alpha) + 1e-8, "{x} {t}");
            assert_eq!(m.minimum(), f(x));
            assert!(calls.iter().all(|&u| u >= a.min(c) && u <= a.max(c)));
        }
    }

    #[test]
    fn kinks_and_monotone_functions() {
        let mut m = BrentMinimum::new(1e-9, 200, 1e-12);
        m.perform(&mut |x| Some((x - 0.123).abs()), -1.0, 0.7, 1.0);
        assert!(m.is_done() && (m.location() - 0.123).abs() < 1e-8);
        // Increasing: next to the lower end.
        m.perform(&mut |x| Some(x), 2.0, 2.5, 3.0);
        assert!(m.is_done() && m.location() - 2.0 <= 2.0 * (1e-9 * 2.0 + 1e-12) + 1e-15);
    }

    #[test]
    fn failed_evaluation_stops_undone() {
        let mut m = BrentMinimum::new(1e-9, 100, 1e-12);
        m.perform(&mut |x| Some(x * x), -1.0, 0.5, 1.0);
        assert!(m.is_done());
        m.perform(&mut |_| None, -1.0, 0.5, 1.0);
        assert!(!m.is_done());
        let mut n = BrentMinimum::new(1e-9, 100, 1e-12);
        n.perform(&mut |x| (x > 0.2).then_some(x * x), -1.0, 0.5, 1.0);
        assert!(!n.is_done());
    }

    #[test]
    fn iteration_cap() {
        let mut m = BrentMinimum::new(1e-15, 20, 1e-300);
        let mut count = 0;
        m.perform(&mut |x| {
            count += 1;
            Some((1.0 / x).sin())
        }, 1e-4, 0.5e-3, 1e-3);
        assert!(count <= 21);
    }
}
