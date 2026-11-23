//! Roots of a function of one variable on an interval: every root of `F − k` from a uniform
//! sampling (sign changes, exact zeros, ends within a value tolerance, and tangential roots
//! where `|F − k|` dips below the tolerance without changing sign), stretches where `F` is
//! null, and Brent's bracketed root finder.
//!
//! Also the callback traits the scalar solvers of the crate take.

/// A function of one variable.
pub(crate) trait Function {
    /// `F(x)` into `value`; `false` where it cannot be evaluated (`value` is then not to be
    /// trusted).
    fn value(&mut self, x: f64, value: &mut f64) -> bool;
}

/// A function of one variable with its derivative.
///
/// The solvers call [`get_state_number`](Self::get_state_number) each time they record a
/// root, right after evaluating the function at exactly that root (and nowhere else in
/// between), so that an implementor can keep whatever it computed there.
pub(crate) trait FunctionWithDerivative {
    fn value(&mut self, x: f64, value: &mut f64) -> bool;
    fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool;
    fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool;
    /// Notification that the last evaluation was at a root; the number returned is kept with
    /// the root.
    fn get_state_number(&mut self) -> i32;
}

/// The result of [`bracketed_root`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct BracketedRoot {
    pub done: bool,
    pub root: f64,
    pub value: f64,
    pub iterations: i32,
}

/// Brent's root finder (inverse quadratic and secant steps safeguarded by bisection) on a
/// bracket `[bound1, bound2]` (either order) across which `F` changes sign. Converged when the
/// bracket around the best point is within `2·zeps·|x| + tolerance/2` of it, or `F` there is
/// exactly zero. At most `iteration_limit` evaluations after the two at the ends; no
/// evaluation outside the bracket. Without a sign change: not done, at `bound2`.
pub(crate) fn bracketed_root(function: &mut dyn Function, bound1: f64, bound2: f64, tolerance: f64, iteration_limit: i32, zeps: f64) -> BracketedRoot {
    let mut eval = |x: f64| -> Option<f64> {
        let mut v = 0.0;
        (function.value(x, &mut v) && !v.is_nan()).then_some(v)
    };
    let (Some(fa), Some(fb)) = (eval(bound1), eval(bound2)) else {
        return BracketedRoot { done: false, root: bound2, value: f64::NAN, iterations: 0 };
    };
    brent_zero(&mut eval, bound1, fa, bound2, fb, tolerance, zeps, iteration_limit)
}

/// Brent's zero finder on `[a, b]` with the end values given.
#[allow(clippy::too_many_arguments)]
fn brent_zero(eval: &mut dyn FnMut(f64) -> Option<f64>, a: f64, fa: f64, b: f64, fb: f64, tolerance: f64, zeps: f64, iteration_limit: i32) -> BracketedRoot {
    let (mut a, mut fa, mut b, mut fb) = (a, fa, b, fb);
    if fa * fb > 0.0 {
        return BracketedRoot { done: false, root: b, value: fb, iterations: 0 };
    }
    // `b` the best point, `a` the previous one, `c` the other end of the bracket.
    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    let mut iterations = 0;
    loop {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * zeps * b.abs() + 0.5 * tolerance;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return BracketedRoot { done: true, root: b, value: fb, iterations };
        }
        if iterations >= iteration_limit {
            return BracketedRoot { done: false, root: b, value: fb, iterations };
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            // Inverse quadratic interpolation (secant with two points).
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < min1.min(min2) {
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
        iterations += 1;
        match eval(b) {
            Some(v) => fb = v,
            None => return BracketedRoot { done: false, root: b, value: f64::NAN, iterations },
        }
    }
}

/// The result of [`function_roots`].
#[derive(Clone, Debug, Default)]
pub(crate) struct FunctionRoots {
    pub done: bool,
    /// `|F − k|` within the null tolerance at every sample.
    pub all_null: bool,
    /// The roots, ascending.
    pub solutions: Vec<f64>,
    /// The state number of each root.
    pub state_numbers: Vec<i32>,
}

/// The roots being collected, with the protocol of recording one.
struct Recorder<'a> {
    function: &'a mut dyn FunctionWithDerivative,
    k: f64,
    out: FunctionRoots,
    /// Candidates closer than this to a recorded root are the same root.
    duplicate: f64,
}

impl Recorder<'_> {
    /// `F − k` and its derivative at `x`.
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let (mut v, mut d) = (0.0, 0.0);
        (self.function.values(x, &mut v, &mut d) && !v.is_nan()).then_some((v - self.k, d))
    }

    fn value(&mut self, x: f64) -> Option<f64> {
        let mut v = 0.0;
        (self.function.value(x, &mut v) && !v.is_nan()).then_some(v - self.k)
    }

    fn derivative(&mut self, x: f64) -> Option<f64> {
        let mut d = 0.0;
        (self.function.derivative(x, &mut d) && !d.is_nan()).then_some(d)
    }

    /// Records the root `x` unless one is recorded within the duplicate distance: evaluates the
    /// function there and keeps its state number.
    fn record(&mut self, x: f64) {
        if self.out.solutions.iter().any(|&s| (s - x).abs() <= self.duplicate) {
            return;
        }
        let mut v = 0.0;
        let _ = self.function.value(x, &mut v);
        let state = self.function.get_state_number();
        let at = self.out.solutions.partition_point(|&s| s < x);
        self.out.solutions.insert(at, x);
        self.out.state_numbers.insert(at, state);
    }

    /// The root of `F − k` in `[lo, hi]`, where its values `glo`, `ghi` have strictly opposite
    /// signs: Newton steps kept in a shrinking bracket (bisection where a step leaves it or
    /// gains too little), then Newton steps kept while they lower `|F − k|`.
    fn root_in(&mut self, lo: f64, glo: f64, hi: f64, ghi: f64, eps_x: f64) -> Option<f64> {
        // Oriented so that g(xl) < 0 < g(xh).
        let (mut xl, mut xh) = if glo < 0.0 { (lo, hi) } else { (hi, lo) };
        let (outer_lo, outer_hi) = (lo.min(hi), lo.max(hi));
        // From the secant point.
        let mut x = lo - glo * (hi - lo) / (ghi - glo);
        if !(x > outer_lo && x < outer_hi) {
            x = 0.5 * (lo + hi);
        }
        let mut dx_old = (hi - lo).abs();
        let mut dx = dx_old;
        let (mut g, mut dg) = self.values(x)?;
        let mut iterations = 0;
        while g != 0.0 {
            if g < 0.0 {
                xl = x;
            } else {
                xh = x;
            }
            iterations += 1;
            if iterations > 200 {
                return None;
            }
            let leaves = ((x - xh) * dg - g) * ((x - xl) * dg - g) > 0.0;
            let slow = (2.0 * g).abs() > (dx_old * dg).abs();
            dx_old = dx;
            if leaves || slow || !dg.is_finite() || dg == 0.0 {
                dx = 0.5 * (xh - xl);
                x = xl + dx;
            } else {
                dx = g / dg;
                x -= dx;
            }
            (g, dg) = self.values(x)?;
            if dx.abs() <= eps_x {
                break;
            }
        }
        // Polish: Newton steps inside the bracket while |g| drops.
        for _ in 0..8 {
            if g == 0.0 || dg == 0.0 || !dg.is_finite() {
                break;
            }
            let next = x - g / dg;
            if !(next >= outer_lo && next <= outer_hi) || next == x {
                break;
            }
            match self.values(next) {
                Some((gn, dn)) if gn.abs() < g.abs() => {
                    x = next;
                    g = gn;
                    dg = dn;
                }
                _ => break,
            }
        }
        Some(x)
    }

    /// The extremum of `F − k` in `[lo, hi]` around a sampled one at `mid` (`sign` +1 for a
    /// positive minimum, −1 for a negative maximum), recording the roots of any sign change met
    /// on the way and then the extremum itself when `|F − k|` there is below `eps_f`.
    fn tangential(&mut self, lo: f64, mid: f64, hi: f64, sign: f64, values: [f64; 3], eps_x: f64, eps_f: f64) {
        // h = sign·(F − k) is positive at the three samples; its minimum is sought.
        let mut estimates: Vec<(f64, f64, f64)> = Vec::new(); // (x, h, |h'|)
        // From the derivative changing sign across the span.
        let d_lo = self.derivative(lo);
        let d_hi = self.derivative(hi);
        if let (Some(d_lo), Some(d_hi)) = (d_lo, d_hi) {
            let (s_lo, s_hi) = (sign * d_lo, sign * d_hi);
            if s_lo < 0.0 && s_hi > 0.0 {
                let mut slope = |x: f64| self.derivative(x).map(|d| sign * d);
                let r = brent_zero(&mut slope, lo, s_lo, hi, s_hi, eps_x.min(1e-15 * (lo.abs() + hi.abs())) + 1e-300, 1e-16, 200);
                if r.root >= lo && r.root <= hi {
                    if let Some((g, d)) = self.values(r.root) {
                        estimates.push((r.root, sign * g, d.abs()));
                    }
                }
            }
        }
        // Without derivative, or where that lands on something no lower than the sample.
        let best_sample = values[1];
        if estimates.first().is_none_or(|e| e.1 > best_sample) {
            let rho = 1e-10;
            let mut minimum = super::brent_minimum::BrentMinimum::new(rho, 100, eps_x.max(1e-300));
            {
                let k = self.k;
                let f = &mut *self.function;
                let mut h = |x: f64| {
                    let mut v = 0.0;
                    (f.value(x, &mut v) && !v.is_nan()).then_some(sign * (v - k))
                };
                minimum.perform(&mut h, lo, mid, hi);
            }
            if let Some((x, h)) = minimum.best() {
                if let Some(d) = self.derivative(x) {
                    estimates.push((x, h, d.abs()));
                }
            }
        }
        let Some(&best) = estimates.iter().min_by(|a, b| {
            if (a.1 - b.1).abs() <= eps_f { a.2.total_cmp(&b.2) } else { a.1.total_cmp(&b.1) }
        }) else {
            return;
        };
        let (x, h, _) = best;
        if h < 0.0 {
            // Two close roots the sampling missed, either side of the extremum.
            let g_lo = sign * values[0];
            let g_hi = sign * values[2];
            let g_x = sign * h;
            if let Some(r) = self.root_in(lo, g_lo, x, g_x, eps_x) {
                self.record(r);
            }
            if let Some(r) = self.root_in(x, g_x, hi, g_hi, eps_x) {
                self.record(r);
            }
            return;
        }
        if h < eps_f {
            self.record(x);
        }
    }
}

/// Every root of `F − k` on `[a, b]` (either order): `F − k` sampled at
/// `max(2·sample_count, 20) + 1` evenly spaced points; then, in this order, the roots of the
/// sign changes between samples, samples where `F − k` is exactly zero, the ends where
/// `|F − k| ≤ eps_f`, and tangential roots (an extremum of `F − k` with `|F − k| < eps_f`
/// where the samples show a positive minimum or a negative maximum). Roots are located to
/// `max(eps_x_given, 1e-10·(|a| + |b|))` and refined to full precision; a root within
/// `1e-10·(|a| + |b|)` of one already recorded is the same. Not done when a sample cannot be
/// evaluated; all null (no roots) when `|F − k| ≤ eps_null` at every sample.
#[allow(clippy::too_many_arguments)]
pub(crate) fn function_roots(function: &mut dyn FunctionWithDerivative, a: f64, b: f64, sample_count: i32, eps_x_given: f64, eps_f: f64, eps_null: f64, k: f64) -> FunctionRoots {
    let (x0, xn) = (a.min(b), a.max(b));
    let n = (2 * sample_count.max(0) as usize).max(20);
    let scale = x0.abs() + xn.abs();
    let duplicate = 1e-10 * scale;
    let eps_x = eps_x_given.max(duplicate);
    let mut r = Recorder { function, k, out: FunctionRoots::default(), duplicate };
    let step = (xn - x0) / n as f64;
    // (Each sample a step on from the one before: the last one may miss `xn` by rounding,
    // which decides on which side of a jump of F at `xn` it falls, and with it whether a
    // sign change is seen there; the reference results depend on it.)
    let mut xs: Vec<f64> = Vec::with_capacity(n + 1);
    let mut x = x0;
    for _ in 0..=n {
        xs.push(x);
        x += step;
    }
    let mut gs = Vec::with_capacity(n + 1);
    for &x in &xs {
        match r.value(x) {
            Some(g) => gs.push(g),
            None => return FunctionRoots::default(),
        }
    }
    if gs.iter().all(|g| g.abs() <= eps_null) {
        return FunctionRoots { done: true, all_null: true, ..FunctionRoots::default() };
    }
    r.out.done = true;
    // Sign changes.
    for i in 0..n {
        let (g0, g1) = (gs[i], gs[i + 1]);
        if (g0 < 0.0 && g1 > 0.0) || (g0 > 0.0 && g1 < 0.0) {
            if let Some(x) = r.root_in(xs[i], g0, xs[i + 1], g1, eps_x) {
                r.record(x);
            }
        }
    }
    // Exact zeros at samples (not where F − k stays zero half a step on).
    for i in 0..=n {
        if gs[i] == 0.0 {
            let ahead = (xs[i] + 0.5 * step).min(xn);
            if r.value(ahead) != Some(0.0) {
                r.record(xs[i]);
            }
        }
    }
    // The ends.
    if gs[0].abs() <= eps_f {
        r.record(x0);
    }
    if gs[n].abs() <= eps_f {
        r.record(xn);
    }
    // Tangential roots.
    for i in 1..n {
        let (gp, g, gn) = (gs[i - 1], gs[i], gs[i + 1]);
        let sign = if g > 0.0 && gp > g && gn > g {
            1.0
        } else if g < 0.0 && gp < g && gn < g {
            -1.0
        } else {
            continue;
        };
        r.tangential(xs[i - 1], xs[i], xs[i + 1], sign, [sign * gp, sign * g, sign * gn], eps_x, eps_f);
    }
    r.out
}

/// Equally spaced samples `a`…`b` (`n` of them, 1-based).
#[derive(Clone, Copy, Debug)]
pub(crate) struct FunctionSample {
    pub a: f64,
    pub b: f64,
    pub n: i32,
}

impl FunctionSample {
    /// Sample `index`, `1 ≤ index ≤ n`: `a` first, `b` last.
    pub(crate) fn parameter(&self, index: i32) -> f64 {
        assert!(index >= 1 && index <= self.n, "sample {index} of {}", self.n);
        let (n, i) = (self.n as f64, index as f64);
        ((n - i) * self.a + (i - 1.0) * self.b) / (n - 1.0)
    }
}

/// A search for roots or null stretches failed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct NumericError;

/// The result of [`function_all_roots`].
#[derive(Clone, Debug, Default)]
pub(crate) struct FunctionAllRoots {
    /// Stretches where `F` is null.
    pub intervals: Vec<(f64, f64)>,
    pub interval_states: Vec<(i32, i32)>,
    /// Isolated roots outside them, ascending.
    pub points: Vec<f64>,
    pub point_states: Vec<i32>,
}

/// The state number at `x` (after evaluating there).
fn state_at(function: &mut dyn FunctionWithDerivative, x: f64) -> i32 {
    let mut v = 0.0;
    let _ = function.value(x, &mut v);
    function.get_state_number()
}

/// The largest (`last`) or smallest root of `F = ±c` on `[lo, hi]` with its state: that of
/// `F = c`, which must exist, or that of `F = −c` beyond it.
fn null_end(function: &mut dyn FunctionWithDerivative, lo: f64, hi: f64, c: f64, eps_x: f64, eps_f: f64, last: bool) -> Result<(f64, i32), NumericError> {
    let pick = |r: &FunctionRoots| -> Option<(f64, i32)> {
        let i = if last { r.solutions.len().checked_sub(1)? } else { (!r.solutions.is_empty()).then_some(0)? };
        Some((r.solutions[i], r.state_numbers[i]))
    };
    let plus = function_roots(function, lo, hi, 10, eps_x, eps_f, 0.0, c);
    if !plus.done || plus.all_null {
        return Err(NumericError);
    }
    let mut best = pick(&plus).ok_or(NumericError)?;
    let minus = function_roots(function, lo, hi, 10, eps_x, eps_f, 0.0, -c);
    if !minus.done || minus.all_null {
        return Err(NumericError);
    }
    if let Some(other) = pick(&minus) {
        if (last && other.0 > best.0) || (!last && other.0 < best.0) {
            best = other;
        }
    }
    Ok(best)
}

/// The stretches of the samples `s` where `F` is null (`|F| ≤ eps_nul` at two or more
/// consecutive samples; their ends where `|F|` reaches `eps_nul`), and the isolated roots of
/// `F` outside them (see [`function_roots`], with `eps_x`, `eps_f`).
pub(crate) fn function_all_roots(function: &mut dyn FunctionWithDerivative, s: &FunctionSample, eps_x: f64, eps_f: f64, eps_nul: f64) -> Result<FunctionAllRoots, NumericError> {
    let n = s.n;
    let mut out = FunctionAllRoots::default();
    let mut values = Vec::with_capacity(n.max(0) as usize);
    for i in 1..=n {
        let mut v = 0.0;
        if !function.value(s.parameter(i), &mut v) {
            return Err(NumericError);
        }
        values.push(v);
    }
    let null = |v: f64| v.abs() <= eps_nul;
    let mut v_last = 0.0;
    let mut open: Option<(f64, i32)> = None;
    for i in 1..=n {
        let v = values[(i - 1) as usize];
        if null(v) {
            if open.is_none() && i >= 2 && null(values[(i - 2) as usize]) {
                // Opens between samples i − 1 and i.
                let start = if i == 2 {
                    let x = s.parameter(1);
                    (x, state_at(function, x))
                } else {
                    let c = eps_nul * if v_last < 0.0 { -1.0 } else { 1.0 };
                    null_end(function, s.parameter(i - 2), s.parameter(i - 1), c, eps_x, eps_f, true)?
                };
                open = Some(start);
            }
        } else {
            if let Some((start, start_state)) = open.take() {
                let c = eps_nul * if v < 0.0 { -1.0 } else { 1.0 };
                let (end, end_state) = null_end(function, s.parameter(i - 1), s.parameter(i), c, eps_x, eps_f, false)?;
                out.intervals.push((start, end));
                out.interval_states.push((start_state, end_state));
            }
            v_last = v;
        }
    }
    if let Some((start, start_state)) = open {
        let x = s.parameter(n);
        out.intervals.push((start, x));
        out.interval_states.push((start_state, state_at(function, x)));
    }
    // Isolated roots between the stretches.
    let (s1, sn) = (s.parameter(1), s.parameter(n));
    let mut ranges: Vec<(f64, f64, i32)> = Vec::new();
    if out.intervals.is_empty() {
        ranges.push((s1, sn, n));
    } else {
        let length = (sn - s1).abs();
        let count = |a: f64, b: f64| -> i32 { (((n as f64) * (b - a).abs() / length) as i32).max(3) };
        let mut from = s1;
        let mut skip_first = out.intervals[0].0 == s1;
        for &(start, end) in &out.intervals {
            if !skip_first {
                ranges.push((from, start, count(from, start)));
            }
            skip_first = false;
            from = end;
        }
        if out.intervals[out.intervals.len() - 1].1 != sn {
            ranges.push((from, sn, count(from, sn)));
        }
    }
    for (a, b, count) in ranges {
        let r = function_roots(function, a, b, count, eps_x, eps_f, 0.0, 0.0);
        if !r.done || r.all_null {
            return Err(NumericError);
        }
        out.points.extend(r.solutions);
        out.point_states.extend(r.state_numbers);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A function from closures, counting state numbers.
    struct Test<F: FnMut(f64) -> Option<(f64, f64)>> {
        f: F,
        last: f64,
        states: Vec<f64>,
    }

    impl<F: FnMut(f64) -> Option<(f64, f64)>> FunctionWithDerivative for Test<F> {
        fn value(&mut self, x: f64, value: &mut f64) -> bool {
            self.last = x;
            (self.f)(x).map(|(v, _)| *value = v).is_some()
        }
        fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool {
            self.last = x;
            (self.f)(x).map(|(_, d)| *derivative = d).is_some()
        }
        fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool {
            self.last = x;
            (self.f)(x).map(|(v, d)| (*value, *derivative) = (v, d)).is_some()
        }
        fn get_state_number(&mut self) -> i32 {
            self.states.push(self.last);
            self.states.len() as i32
        }
    }

    fn test<F: FnMut(f64) -> Option<(f64, f64)>>(f: F) -> Test<F> {
        Test { f, last: f64::NAN, states: Vec::new() }
    }

    #[test]
    fn roots_of_a_cubic_in_order() {
        let mut f = test(|x| Some(((x - 0.13) * (x - 0.47) * (x - 0.83), (x - 0.47) * (x - 0.83) + (x - 0.13) * (x - 0.83) + (x - 0.13) * (x - 0.47))));
        let r = function_roots(&mut f, 0.0, 1.0, 10, 1e-9, 1e-12, 1e-12, 0.0);
        assert!(r.done && !r.all_null);
        assert_eq!(r.solutions.len(), 3);
        for (x, e) in r.solutions.iter().zip([0.13, 0.47, 0.83]) {
            assert!((x - e).abs() < 1e-14, "{x} {e}");
        }
        assert_eq!(r.state_numbers, [1, 2, 3]);
        // Each state recorded right after an evaluation at the root.
        assert_eq!(f.states, r.solutions);
    }

    #[test]
    fn tangential_root_found_by_rediscretization() {
        let mut f = test(|x| Some(((x - 0.31) * (x - 0.31), 2.0 * (x - 0.31))));
        let r = function_roots(&mut f, 0.0, 1.0, 10, 1e-9, 1e-8, 1e-12, 0.0);
        assert_eq!(r.solutions.len(), 1);
        assert!((r.solutions[0] - 0.31).abs() < 1e-12, "{:?}", r.solutions);
    }

    #[test]
    fn exact_double_root() {
        let mut f = test(|t| Some(((2.0 * t - 1.0).powi(2), 4.0 * (2.0 * t - 1.0))));
        let r = function_roots(&mut f, 0.03, 1.07, 7, 1e-14, 1e-8, 1e-8, 0.0);
        assert_eq!(r.solutions.len(), 1, "{:?}", r.solutions);
        assert!((r.solutions[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn sine_roots_ascending() {
        let mut f = test(|x| Some((x.sin(), x.cos())));
        let r = function_roots(&mut f, 10.0, 0.0, 10, 1e-12, 1e-12, 1e-12, 0.0);
        let pi = std::f64::consts::PI;
        assert_eq!(r.solutions.len(), 4, "{:?}", r.solutions);
        for (i, x) in r.solutions.iter().enumerate() {
            assert!((x - pi * i as f64).abs() < 1e-14, "{x}");
        }
        // Recorded in the order: sign changes (π, 2π, 3π), then the end 0.
        assert_eq!(f.states, [pi, 2.0 * pi, 3.0 * pi, 0.0].map(|x| r.solutions.iter().copied().find(|s| (s - x).abs() < 1e-12).unwrap()));
    }

    #[test]
    fn zero_at_a_sample_and_at_the_ends() {
        // Zero exactly at the eleventh of the 21 samples on [0, 1] (each a step on from the
        // last), and within eps_f at both ends.
        let s = (0..10).fold(0.0, |x, _| x + 0.05);
        let mut f = test(move |x| Some(((x - s) * (x - s) * (x - s), 3.0 * (x - s) * (x - s))));
        let r = function_roots(&mut f, 0.0, 1.0, 10, 1e-12, 1e-9, 1e-12, 0.0);
        assert!(r.solutions.contains(&s), "{:?}", r.solutions);
        let mut g = test(|x| Some((x * (1.0 - x) * 1e-10, (1.0 - 2.0 * x) * 1e-10)));
        let r = function_roots(&mut g, 0.0, 1.0, 10, 1e-12, 1e-9, 1e-12, 0.0);
        assert_eq!(r.solutions, [0.0, 1.0]);
    }

    #[test]
    fn close_pair_and_near_miss() {
        let mut f = test(|x| Some((x * x - 1e-12, 2.0 * x)));
        let r = function_roots(&mut f, -0.97, 1.0, 10, 1e-14, 1e-8, 1e-14, 0.0);
        assert_eq!(r.solutions.len(), 2, "{:?}", r.solutions);
        assert!((r.solutions[0] + 1e-6).abs() < 1e-15 && (r.solutions[1] - 1e-6).abs() < 1e-15);
        let mut g = test(|x| Some((x * x + 1e-9, 2.0 * x)));
        let r = function_roots(&mut g, -0.97, 1.0, 10, 1e-14, 1e-8, 1e-14, 0.0);
        assert_eq!(r.solutions.len(), 1);
        assert!(r.solutions[0].abs() < 1e-12);
        // Above the tolerance: nothing.
        let mut h = test(|x| Some((x * x + 1e-7, 2.0 * x)));
        assert!(function_roots(&mut h, -0.97, 1.0, 10, 1e-14, 1e-8, 1e-14, 0.0).solutions.is_empty());
    }

    #[test]
    fn constant_and_failing_functions() {
        let mut f = test(|_| Some((1e-13, 0.0)));
        let r = function_roots(&mut f, 0.0, 1.0, 10, 1e-9, 1e-12, 1e-12, 0.0);
        assert!(r.done && r.all_null && r.solutions.is_empty());
        let mut g = test(|x| ((x - 0.5).abs() > 0.01).then_some((x - 0.3, 1.0)));
        let r = function_roots(&mut g, 0.0, 1.0, 10, 1e-9, 1e-12, 1e-12, 0.0);
        assert!(!r.done && r.solutions.is_empty() && g.states.is_empty());
    }

    #[test]
    fn brent_on_a_bracket() {
        struct Square;
        impl Function for Square {
            fn value(&mut self, x: f64, value: &mut f64) -> bool {
                *value = x * x - 2.0;
                true
            }
        }
        let r = bracketed_root(&mut Square, 0.0, 2.0, 1e-12, 100, 1e-16);
        assert!(r.done);
        assert!((r.root - 2f64.sqrt()).abs() < 1e-11);
        assert!(r.iterations < 15, "{}", r.iterations);
        let r = bracketed_root(&mut Square, 2.0, 3.0, 1e-12, 100, 1e-16);
        assert!(!r.done && r.iterations == 0 && r.root == 3.0 && r.value == 7.0);
    }

    #[test]
    fn samples() {
        let s = FunctionSample { a: 0.1, b: 0.7, n: 7 };
        assert!((s.parameter(1) - 0.1).abs() < 1e-16);
        assert!((s.parameter(7) - 0.7).abs() < 1e-15);
        assert!((s.parameter(4) - 0.4).abs() < 1e-15);
    }

    #[test]
    fn null_stretch_and_isolated_roots() {
        // Zero on [−1, 1], (x ∓ 1)(2.5 − x) beyond: roots at 2.5 too.
        let f = |x: f64| -> (f64, f64) {
            if x.abs() <= 1.0 {
                (0.0, 0.0)
            } else if x > 0.0 {
                ((x - 1.0) * (2.5 - x), 3.5 - 2.0 * x)
            } else {
                ((x + 1.0) * (2.5 - x), 1.5 - 2.0 * x)
            }
        };
        let mut t = test(|x| Some(f(x)));
        let s = FunctionSample { a: -3.0, b: 3.0, n: 25 };
        let r = function_all_roots(&mut t, &s, 1e-12, 1e-8, 1e-8).unwrap();
        assert_eq!(r.intervals.len(), 1, "{:?}", r.intervals);
        let (a, b) = r.intervals[0];
        assert!((a + 1.0).abs() < 1e-6 && (b - 1.0).abs() < 1e-6, "{a} {b}");
        assert!(r.points.iter().any(|p| (p - 2.5).abs() < 1e-12), "{:?}", r.points);
        assert_eq!(r.points.len(), r.point_states.len());
    }
}
