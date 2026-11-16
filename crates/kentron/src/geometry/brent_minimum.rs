//! OCCT's one-variable minimizer `math_BrentMinimum` (Brent's method: parabolic
//! interpolation with golden-section steps), ported.

/// `CGOLD`: 0.5 * (3 - sqrt(5)), as OCCT writes it.
const GOLDEN_SECTION_RATIO: f64 = 0.3819660;

/// `math_BrentMinimum`: the minimum of a function of one variable inside a bracketing
/// triplet. The state persists between [`perform`](Self::perform) calls as in OCCT (in
/// particular `IsDone()` is not reset when a call stops on a failed evaluation).
#[doc(alias = "math_BrentMinimum")]
pub(crate) struct BrentMinimum {
    /// `a`, `b`: the current bracket.
    lower: f64,
    upper: f64,
    /// `x`, `fx`: the best point so far and its value.
    location: f64,
    value_at_location: f64,
    /// `fv`, `fw`: the values at the previous two best points.
    value_at_v: f64,
    value_at_w: f64,
    /// `XTol`: the relative tolerance on the location.
    relative_tolerance: f64,
    /// `EPSZ`: the absolute tolerance added to it.
    absolute_tolerance: f64,
    done: bool,
    iterations: usize,
    max_iterations: usize,
}

impl BrentMinimum {
    /// `math_BrentMinimum(TolX, NbIterations, ZEPS)`.
    pub(crate) fn new(relative_tolerance: f64, max_iterations: usize, absolute_tolerance: f64) -> Self {
        BrentMinimum {
            lower: 0.0,
            upper: 0.0,
            location: 0.0,
            value_at_location: 0.0,
            value_at_v: 0.0,
            value_at_w: 0.0,
            relative_tolerance,
            absolute_tolerance,
            done: false,
            iterations: 0,
            max_iterations,
        }
    }

    /// `math_BrentMinimum::IsSolutionReached`.
    fn is_solution_reached(&self) -> bool {
        let two_tolerance = 2.0 * (self.relative_tolerance * self.location.abs() + self.absolute_tolerance);
        self.location <= two_tolerance + self.lower && self.location >= self.upper - two_tolerance
    }

    /// `math_BrentMinimum::Perform(F, Ax, Bx, Cx)`: `function` returns `None` where OCCT's
    /// `math_Function::Value` returns false (the search then stops where it is).
    #[doc(alias = "math_BrentMinimum::Perform")]
    pub(crate) fn perform(&mut self, function: &mut dyn FnMut(f64) -> Option<f64>, ax: f64, bx: f64, cx: f64) {
        let mut e: f64 = 0.0;
        let mut d: f64 = f64::MAX;
        self.lower = if ax < cx { ax } else { cx };
        self.upper = if ax > cx { ax } else { cx };
        let mut x = bx;
        let mut w = bx;
        let mut v = bx;
        self.location = x;
        match function(x) {
            Some(value) => self.value_at_location = value,
            None => return,
        }
        self.value_at_w = self.value_at_location;
        self.value_at_v = self.value_at_location;
        self.iterations = 1;
        while self.iterations <= self.max_iterations {
            let xm = 0.5 * (self.lower + self.upper);
            let tol1 = self.relative_tolerance * x.abs() + self.absolute_tolerance;
            let tol2 = 2.0 * tol1;
            if self.is_solution_reached() {
                self.done = true;
                return;
            }
            let fx = self.value_at_location;
            if e.abs() > tol1 {
                let r = (x - w) * (fx - self.value_at_v);
                let mut q = (x - v) * (fx - self.value_at_w);
                let mut p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if q > 0.0 {
                    p = -p;
                }
                q = q.abs();
                let etemp = e;
                e = d;
                if p.abs() >= (0.5 * q * etemp).abs() || p <= q * (self.lower - x) || p >= q * (self.upper - x) {
                    e = if x >= xm { self.lower - x } else { self.upper - x };
                    d = GOLDEN_SECTION_RATIO * e;
                } else {
                    d = p / q;
                    let u = x + d;
                    if u - self.lower < tol2 || self.upper - u < tol2 {
                        d = tol1.copysign(xm - x);
                    }
                }
            } else {
                e = if x >= xm { self.lower - x } else { self.upper - x };
                d = GOLDEN_SECTION_RATIO * e;
            }
            let u = if d.abs() >= tol1 { x + d } else { x + tol1.copysign(d) };
            let fu = match function(u) {
                Some(value) => value,
                None => return,
            };
            if fu <= fx {
                if u >= x {
                    self.lower = x;
                } else {
                    self.upper = x;
                }
                // SHFT(v, w, x, u); SHFT(fv, fw, fx, fu)
                v = w;
                w = x;
                x = u;
                self.location = x;
                self.value_at_v = self.value_at_w;
                self.value_at_w = fx;
                self.value_at_location = fu;
            } else {
                if u < x {
                    self.lower = u;
                } else {
                    self.upper = u;
                }
                if fu <= self.value_at_w || w == x {
                    v = w;
                    w = u;
                    self.value_at_v = self.value_at_w;
                    self.value_at_w = fu;
                } else if fu <= self.value_at_v || v == x || v == w {
                    v = u;
                    self.value_at_v = fu;
                }
            }
            self.iterations += 1;
        }
        self.done = false;
    }

    #[doc(alias = "math_BrentMinimum::IsDone")]
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// `Location()`; panics where OCCT raises `StdFail_NotDone`.
    #[doc(alias = "math_BrentMinimum::Location")]
    pub(crate) fn location(&self) -> f64 {
        assert!(self.done, "math_BrentMinimum not done");
        self.location
    }

    /// `Minimum()`; panics where OCCT raises `StdFail_NotDone`.
    #[doc(alias = "math_BrentMinimum::Minimum")]
    pub(crate) fn minimum(&self) -> f64 {
        assert!(self.done, "math_BrentMinimum not done");
        self.value_at_location
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parabola_minimum() {
        let mut brent = BrentMinimum::new(1.0e-10, 100, 1.0e-12);
        brent.perform(&mut |x| Some((x - 0.3) * (x - 0.3) + 2.0), -1.0, 0.0, 1.0);
        assert!(brent.is_done());
        // The value is flat to double precision within ~1e-8 of the minimum.
        assert!((brent.location() - 0.3).abs() < 1.0e-7);
    }

    /// Against OCCT 8.1's `math_BrentMinimum` on `cos(3x) + 0.2x`, undefined outside the
    /// bracket: `(ax, bx, cx, TolX, ZEPS)` to `(location, minimum, iterations)`.
    #[test]
    fn matches_occt() {
        let cases = [
            ([-1.0, 0.1, 2.0, 1e-3, 1e-9], (1.0256949747715178, -0.7927811028765742, 7)),
            ([0.0, 1.0, 3.0, 1e-6, 1e-12], (1.0249593071700436, -0.7927835361267458, 9)),
            ([-2.0, -1.9, 4.0, 1e-10, 1e-12], (-1.0694362655111846, -1.2116625566063857, 18)),
            ([0.5, 1.5, 1.6, 1e-3, 1e-5], (1.0252880508800348, -0.792783049478366, 7)),
        ];
        for ([ax, bx, cx, tol_x, zeps], (location, minimum, iterations)) in cases {
            let mut brent = BrentMinimum::new(tol_x, 100, zeps);
            brent.perform(&mut |x| if x < ax || x > cx { None } else { Some((3.0 * x).cos() + 0.2 * x) }, ax, bx, cx);
            assert!(brent.is_done());
            assert_eq!((brent.location(), brent.minimum(), brent.iterations), (location, minimum, iterations));
        }
    }

    #[test]
    fn failed_evaluation_stops_undone() {
        let mut brent = BrentMinimum::new(1.0e-10, 100, 1.0e-12);
        brent.perform(&mut |x| if x > 0.5 { None } else { Some(-x) }, 0.0, 0.25, 1.0);
        assert!(!brent.is_done());
    }
}
