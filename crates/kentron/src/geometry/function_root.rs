//! One root of a function of one variable near a guess, within bounds: a bounded local
//! minimisation of `½F²`. Newton steps, cut to the bounds and to half their width; a step
//! across a sign change is finished by safeguarded Newton in that bracket; a step that gains
//! little (as towards a double root) is extended while `|F|` keeps dropping, and one that
//! overshoots goes back to the minimum of `|F|` along it, found as a root of `F′` (signed by
//! `F`) by Brent's method, so that a tangency is located as precisely as a simple root.

use super::function_roots::FunctionWithDerivative;

/// The result of [`function_root`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct FunctionRoot {
    pub done: bool,
    /// The point reached (the root when done).
    pub root: f64,
    /// `F′` at `root`.
    pub derivative: f64,
    /// Iterations performed (0 when the guess already satisfies the stopping test).
    pub iterations: usize,
}

/// Values at least this large stand for infinity.
const HUGE: f64 = 0.5 * crate::tolerance::INFINITE;

/// A value or gradient this small is zero.
fn negligible(v: f64) -> bool {
    v.abs() <= 1e-200
}

/// The function with its evaluation bookkeeping.
struct Solver<'a> {
    function: &'a mut dyn FunctionWithDerivative,
    /// The abscissa of the last evaluation.
    last: f64,
    /// Every evaluation `(x, F)` since the last line search began.
    seen: Vec<(f64, f64)>,
}

impl Solver<'_> {
    /// `(F, F′)` at `x`, `None` where it fails or is infinite.
    fn eval(&mut self, x: f64) -> Option<(f64, f64)> {
        let (mut f, mut d) = (0.0, 0.0);
        self.last = x;
        if !self.function.values(x, &mut f, &mut d) || !f.is_finite() || f.abs() >= HUGE || !(f * d).is_finite() {
            return None;
        }
        self.seen.push((x, f));
        Some((f, d))
    }

    /// The root in the bracket `[x0, x1]` (values of opposite signs): safeguarded Newton, then
    /// Newton steps while `|F|` drops. Returns the root, F and F′ there.
    fn bracketed(&mut self, (x0, f0): (f64, f64), (x1, f1): (f64, f64), tolerance: f64) -> Option<(f64, f64, f64)> {
        let (lo, hi) = (x0.min(x1), x0.max(x1));
        let (mut xl, mut xh) = if f0 < 0.0 { (x0, x1) } else { (x1, x0) };
        let mut x = x0 - f0 * (x1 - x0) / (f1 - f0);
        if !(x > lo && x < hi) {
            x = 0.5 * (x0 + x1);
        }
        let mut dx_old = hi - lo;
        let mut dx = dx_old;
        let (mut f, mut d) = self.eval(x)?;
        let tol = tolerance.min(1e-3 * (hi - lo)).max(0.0);
        for _ in 0..200 {
            if f == 0.0 {
                return Some((x, f, d));
            }
            if f < 0.0 {
                xl = x;
            } else {
                xh = x;
            }
            let leaves = ((x - xh) * d - f) * ((x - xl) * d - f) > 0.0;
            let slow = (2.0 * f).abs() > (dx_old * d).abs();
            dx_old = dx;
            if leaves || slow || d == 0.0 {
                dx = 0.5 * (xh - xl);
                x = xl + dx;
            } else {
                dx = f / d;
                x -= dx;
            }
            (f, d) = self.eval(x)?;
            if dx.abs() <= tol || xh == xl {
                break;
            }
        }
        for _ in 0..8 {
            if f == 0.0 || d == 0.0 {
                break;
            }
            let next = x - f / d;
            if !(next >= lo && next <= hi) || next == x {
                break;
            }
            let (fn_, dn) = self.eval(next)?;
            if fn_.abs() >= f.abs() {
                break;
            }
            (x, f, d) = (next, fn_, dn);
        }
        Some((x, f, d))
    }

    /// The minimum of `|F|` between `x0` and `x1` where `F` has the sign `sign` at `x0`: a
    /// root of `sign·F′` (towards `x1`) where it changes sign, otherwise Brent's minimiser on
    /// `|F|`. Where the search meets a value of the other sign, the root between instead (the
    /// flag set). Returns the point, F and F′ there.
    fn line_minimum(&mut self, x0: f64, x1: f64, sign: f64, tolerance: f64) -> Option<(f64, f64, f64, bool)> {
        self.seen.clear();
        let towards = if x1 >= x0 { 1.0 } else { -1.0 };
        let (f0, d0) = self.eval(x0)?;
        let (f1, d1) = self.eval(x1)?;
        let (s0, s1) = (sign * d0 * towards, sign * d1 * towards);
        let mut found = None;
        if s0 < 0.0 && s1 > 0.0 {
            let mut failed = false;
            let mut slope = |x: f64| -> Option<f64> {
                match self.eval(x) {
                    Some((_, d)) => Some(sign * d * towards),
                    None => {
                        failed = true;
                        None
                    }
                }
            };
            let tol = 1e-15 * (x0.abs() + x1.abs()) + 1e-300;
            let r = brent_zero(&mut slope, x0, s0, x1, s1, tol, 400);
            if !failed {
                found = Some(r);
            }
        }
        if found.is_none() {
            let mut minimum = super::brent_minimum::BrentMinimum::new(1e-12, 200, tolerance.max(1e-300));
            let mid = 0.5 * (x0 + x1);
            minimum.perform(&mut |x| self.eval(x).map(|(f, _)| f.abs()), x0, mid, x1);
            let (x, _) = minimum.best()?;
            found = Some(x);
        }
        // A value of the other sign met: the root between it and the nearest one of this sign.
        let crossing = self.seen.iter().copied().filter(|&(_, f)| f * sign < 0.0).min_by(|a, b| (a.0 - x0).abs().total_cmp(&(b.0 - x0).abs()));
        if let Some((z, fz)) = crossing {
            let (w, fw) = self.seen.iter().copied().filter(|&(x, f)| f * sign > 0.0 && (x - x0) * towards <= (z - x0) * towards).min_by(|a, b| (a.0 - z).abs().total_cmp(&(b.0 - z).abs())).unwrap_or((x0, f0));
            let (r, fr, dr) = self.bracketed((w, fw), (z, fz), tolerance)?;
            return Some((r, fr, dr, true));
        }
        let x = found.expect("a point");
        let (best_x, best_f) = [(x0, f0.abs()), (x1, f1.abs())].into_iter().fold((x, f64::INFINITY), |b, c| if c.1 < b.1 { c } else { b });
        let x = if self.eval(x)?.0.abs() <= best_f { x } else { best_x };
        let (f, d) = self.eval(x)?;
        Some((x, f, d, false))
    }
}

/// Brent's zero finder on `[a, b]` (values of opposite signs given), to `tol`.
fn brent_zero(g: &mut dyn FnMut(f64) -> Option<f64>, a: f64, fa: f64, b: f64, fb: f64, tol: f64, limit: usize) -> f64 {
    let (mut a, mut fa, mut b, mut fb) = (a, fa, b, fb);
    let (mut c, mut fc) = (a, fa);
    let (mut d, mut e) = (b - a, b - a);
    for _ in 0..limit {
        if (fb > 0.0) == (fc > 0.0) && fb != 0.0 && fc != 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            (a, b, c) = (b, c, b);
            (fa, fb, fc) = (fb, fc, fb);
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return b;
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let (qq, r) = (fa / fc, fb / fc);
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol1 { d } else { tol1.copysign(xm) };
        match g(b) {
            Some(v) => fb = v,
            None => return a,
        }
    }
    b
}

/// A root of `function` in `[a, b]` near `guess` (clamped into the bounds), or where `½F²` is
/// locally least in them: a tangency, or a bound towards which `|F|` drops. Done when `F` or
/// `F·F′` vanishes, or an iteration moves by at most `tolerance` (a root is then refined to
/// full precision); not done when an evaluation fails or after `iteration_limit` iterations.
/// On success the function was last evaluated at the root, and then told so
/// (`get_state_number`, once).
pub(crate) fn function_root(function: &mut dyn FunctionWithDerivative, guess: f64, tolerance: f64, a: f64, b: f64, iteration_limit: usize) -> FunctionRoot {
    let mut s = Solver { function, last: f64::NAN, seen: Vec::new() };
    let failed = |x: f64, iterations: usize| FunctionRoot { done: false, root: x, derivative: 0.0, iterations };
    let mut x = guess.clamp(a.min(b), a.max(b));
    let (lo, hi) = (a.min(b), a.max(b));
    let Some((mut f, mut d)) = s.eval(x) else { return failed(x, 0) };
    let finish = |s: &mut Solver, x: f64, d: f64, iterations: usize| {
        let mut d = d;
        if s.last != x {
            match s.eval(x) {
                Some((_, dx)) => d = dx,
                None => return FunctionRoot { done: false, root: x, derivative: d, iterations },
            }
        }
        s.function.get_state_number();
        FunctionRoot { done: true, root: x, derivative: d, iterations }
    };
    if negligible(f) || negligible(f * d) {
        return finish(&mut s, x, d, 0);
    }
    let extent = (hi - lo).min(1e100);
    let max_step = 0.5 * extent;
    for iteration in 1..=iteration_limit {
        let sign = if f > 0.0 { 1.0 } else { -1.0 };
        let newton = -f / d;
        let mut step = newton;
        if !step.is_finite() || step.abs() > max_step {
            step = max_step.copysign(if step.is_nan() { -f * d } else { step });
        }
        let y = (x + step).clamp(lo, hi);
        if y == x {
            // At a bound with the descent pointing out, or a step below rounding.
            return finish(&mut s, x, d, iteration);
        }
        let Some((fy, dy)) = s.eval(y) else { return failed(x, iteration) };
        // (A step cut short, at a bound or to half the range, is searched along: |F| may
        // have dropped further inside, or changed sign there.)
        let cut = !(newton.abs() <= (y - x).abs());
        let (x_new, f_new, d_new);
        if fy == 0.0 {
            return finish(&mut s, y, dy, iteration);
        } else if (fy > 0.0) != (f > 0.0) {
            // A root between: finished in the bracket.
            let Some((r, fr, dr)) = s.bracketed((x, f), (y, fy), tolerance) else { return failed(x, iteration) };
            let _ = fr;
            return finish(&mut s, r, dr, iteration);
        } else if !cut && fy.abs() < f.abs() && {
            // (The linear model's predicted drop of ½F² against the actual one: a step that
            // gains much less than predicted is searched along instead.)
            let predicted = f * f - (f + d * (y - x)).powi(2);
            f * f - fy * fy >= 0.25 * predicted
        } {
            let (mut prev, mut best, mut fbest, mut dbest) = (x, y, fy, dy);
            // Slow, as towards a multiple root: |F| dropping by a constant factor, the next
            // Newton step on the same way and not much shorter.
            let next = -fy / dy;
            let linear = fy.abs() > 0.1 * f.abs() && next * (y - x) > 0.0 && next.abs() >= 0.3 * (y - x).abs();
            if linear && y > lo && y < hi {
                // On along the step while |F| drops.
                let mut length = y - x;
                for _ in 0..30 {
                    length *= 2.0;
                    let z = (x + length).clamp(lo, hi);
                    if z == best {
                        break;
                    }
                    let Some((fz, dz)) = s.eval(z) else { return failed(x, iteration) };
                    if fz == 0.0 {
                        return finish(&mut s, z, dz, iteration);
                    }
                    if (fz > 0.0) != (f > 0.0) {
                        let Some((r, _, dr)) = s.bracketed((best, fbest), (z, fz), tolerance) else { return failed(x, iteration) };
                        return finish(&mut s, r, dr, iteration);
                    }
                    if fz.abs() >= fbest.abs() {
                        // The least |F| lies about `best`: on the side where F′ says.
                        let towards = if z > x { 1.0 } else { -1.0 };
                        let (from, to) = if sign * dbest * towards > 0.0 { (prev, best) } else { (best, z) };
                        if let Some((r, fr, dr, root)) = s.line_minimum(from, to, sign, tolerance) {
                            if root {
                                return finish(&mut s, r, dr, iteration);
                            }
                            if fr.abs() <= fbest.abs() {
                                (best, fbest, dbest) = (r, fr, dr);
                            }
                        }
                        break;
                    }
                    (prev, best, fbest, dbest) = (best, z, fz, dz);
                    if z == lo || z == hi {
                        break;
                    }
                }
            }
            (x_new, f_new, d_new) = (best, fbest, dbest);
        } else {
            // Overshot without crossing, or gained too little: the least |F| along the step.
            match s.line_minimum(x, y, sign, tolerance) {
                Some((r, _, dr, true)) => return finish(&mut s, r, dr, iteration),
                Some((r, fr, dr, _)) if fr.abs() < f.abs() => (x_new, f_new, d_new) = (r, fr, dr),
                Some(_) => return finish(&mut s, x, d, iteration),
                None => return failed(x, iteration),
            }
        }
        let moved = (x_new - x).abs();
        (x, f, d) = (x_new, f_new, d_new);
        if negligible(f) || negligible(f * d) {
            return finish(&mut s, x, d, iteration);
        }
        if moved <= tolerance {
            return finish(&mut s, x, d, iteration);
        }
    }
    failed(x, iteration_limit)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Test<F: FnMut(f64) -> Option<(f64, f64)>> {
        f: F,
        last: f64,
        evaluations: Vec<f64>,
        states: Vec<f64>,
    }

    impl<F: FnMut(f64) -> Option<(f64, f64)>> FunctionWithDerivative for Test<F> {
        fn value(&mut self, x: f64, value: &mut f64) -> bool {
            self.last = x;
            self.evaluations.push(x);
            (self.f)(x).map(|(v, _)| *value = v).is_some()
        }
        fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool {
            self.last = x;
            self.evaluations.push(x);
            (self.f)(x).map(|(_, d)| *derivative = d).is_some()
        }
        fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool {
            self.last = x;
            self.evaluations.push(x);
            (self.f)(x).map(|(v, d)| (*value, *derivative) = (v, d)).is_some()
        }
        fn get_state_number(&mut self) -> i32 {
            self.states.push(self.last);
            0
        }
    }

    fn test<F: FnMut(f64) -> Option<(f64, f64)>>(f: F) -> Test<F> {
        Test { f, last: f64::NAN, evaluations: Vec::new(), states: Vec::new() }
    }

    #[test]
    fn cube_root() {
        let mut f = test(|x| Some((x * x * x - 2.0, 3.0 * x * x)));
        let r = function_root(&mut f, 0.5, 1e-9, 0.0, 3.0, 100);
        assert!(r.done);
        assert!((r.root - 2f64.cbrt()).abs() < 1e-14, "{}", r.root);
        assert_eq!(f.states, [r.root]);
        assert_eq!(*f.evaluations.last().unwrap(), r.root);
        assert!(r.iterations <= 100);
    }

    #[test]
    fn no_root_stationary() {
        let mut f = test(|x| Some((x * x + 1.0, 2.0 * x)));
        let r = function_root(&mut f, 0.7, 1e-9, -2.0, 2.0, 100);
        assert!(r.done);
        assert!(r.root.abs() < 1e-8, "{}", r.root);
    }

    #[test]
    fn root_beyond_a_bound() {
        let mut f = test(|x| Some((x - 5.0, 1.0)));
        let r = function_root(&mut f, 0.5, 1e-9, 0.0, 1.0, 100);
        assert!(r.done && r.root == 1.0);
        // From a guess outside, clamped.
        let mut g = test(|x| Some((x - 0.25, 1.0)));
        let r = function_root(&mut g, -3.0, 1e-9, 0.0, 1.0, 100);
        assert!(r.done && (r.root - 0.25).abs() < 1e-15);
        assert_eq!(g.evaluations[0], 0.0);
    }

    #[test]
    fn double_root() {
        let mut f = test(|t| Some(((2.0 * t - 1.0) * (2.0 * t - 1.0), 4.0 * (2.0 * t - 1.0))));
        let r = function_root(&mut f, 0.2, 1e-9, -0.5, 1.5, 50);
        assert!(r.done);
        assert!((r.root - 0.5).abs() < 1e-12, "{}", r.root);
        // Near miss: the least |F|.
        let mut g = test(|t| Some(((2.0 * t - 1.0) * (2.0 * t - 1.0) + 1e-3, 4.0 * (2.0 * t - 1.0))));
        let r = function_root(&mut g, 0.2, 1e-9, -0.5, 1.5, 50);
        assert!(r.done && (r.root - 0.5).abs() < 1e-12, "{}", r.root);
    }

    #[test]
    fn failures_and_degenerate_bounds() {
        let mut f = test(|_| None);
        let r = function_root(&mut f, 0.2, 1e-9, 0.0, 1.0, 50);
        assert!(!r.done && f.states.is_empty());
        let mut g = test(|x| Some((x - 0.3, 1.0)));
        let r = function_root(&mut g, 0.2, 1e-9, 0.6, 0.6, 50);
        assert!(r.done && r.root == 0.6);
        let mut h = test(|x| Some((1e101 * x + 1.0, 1e101)));
        assert!(!function_root(&mut h, 0.2, 1e-9, 0.0, 1.0, 50).done);
    }

    #[test]
    fn transcendental() {
        let mut f = test(|x| Some((x.cos() - x, -x.sin() - 1.0)));
        let r = function_root(&mut f, 0.0, 1e-6, -1e100, 2e100, 50);
        assert!(r.done);
        assert!((r.root.cos() - r.root).abs() < 1e-15);
    }
}
