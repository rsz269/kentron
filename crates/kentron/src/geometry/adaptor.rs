//! A read-only view of a curve over a parameter range, for the algorithms that take any
//! curve (point–curve extrema, deflection sampling, arc length, curve–plane intersection,
//! sweep trihedra, composite curves), and small helpers on parameters.
//!
//! The view evaluates a B-spline one-sidedly at the ends of its range (the piece inside the
//! range, not the one beyond a knot at the end) and on the piece just ahead of a knot that a
//! parameter misses by a rounding step, and lists the pieces of a requested continuity.

use std::sync::Arc;

use crate::geometry::curve::Curve;
use crate::geometry::{Continuity, ParametricCurve};
use crate::math::{Point, Vector};

/// A curve viewed over a parameter range `[first_parameter, last_parameter]`.
///
/// `curve_type` and `basis_curve` tell an algorithm which elementary curve it has, so that it
/// can use the curve's own data (a line's direction, a circle's radius, a B-spline's knots);
/// views that are not one of those keep the defaults (`None`).
pub(crate) trait CurveAdaptor {
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    /// Increasing breakpoints from the first to the last parameter (at least those two)
    /// between which the curve is at least `continuity`-continuous.
    fn intervals(&self, continuity: Continuity) -> Vec<f64>;
    fn is_periodic(&self) -> bool;
    /// The period (asked of periodic curves only).
    fn period(&self) -> f64;
    fn value(&self, u: f64) -> Point;
    fn d1(&self, u: f64) -> (Point, Vector);
    fn d2(&self, u: f64) -> (Point, Vector, Vector);
    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector);
    /// The `n`-th derivative, `n ≥ 1`.
    fn dn(&self, u: f64, n: u32) -> Vector;
    /// A parameter step that moves the point along the curve by at most `r3d`.
    fn resolution(&self, r3d: f64) -> f64;
    /// Which elementary curve this is; `None` for any other.
    fn curve_type(&self) -> Option<CurveKind> {
        None
    }
    /// The curve `curve_type` names.
    fn basis_curve(&self) -> Option<&Curve> {
        None
    }
    /// The same curve viewed over `[first, last]`.
    fn trimmed(&self, first: f64, last: f64, tolerance: f64) -> Box<dyn CurveAdaptor> {
        let _ = (first, last, tolerance);
        panic!("this curve view cannot be trimmed")
    }
}

/// The kinds of curve a [`GeomAdaptor`] can view.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    Bezier,
    BSpline,
    Offset,
}

/// One of kentron's curves (never a trimmed one: its basis is taken) over a range, which may
/// reach beyond the curve's own (evaluation then extends it as the curve itself does).
/// Cloning shares the curve.
#[derive(Clone, Debug)]
pub(crate) struct GeomAdaptor {
    curve: Arc<Curve>,
    first: f64,
    last: f64,
}

impl GeomAdaptor {
    pub(crate) fn new(curve: &Curve, first: f64, last: f64) -> GeomAdaptor {
        let mut basis = curve;
        while let Curve::Trimmed(t) = basis {
            basis = t.basis();
        }
        GeomAdaptor { curve: Arc::new(basis.clone()), first, last }
    }

    /// The curve viewed (a trimmed curve's basis).
    pub(crate) fn curve(&self) -> &Curve {
        &self.curve
    }

    /// The same curve over `[first, last]`.
    pub(crate) fn trim(&self, first: f64, last: f64) -> GeomAdaptor {
        GeomAdaptor { curve: Arc::clone(&self.curve), first, last }
    }

    pub(crate) fn kind(&self) -> CurveKind {
        match &*self.curve {
            Curve::Line(_) => CurveKind::Line,
            Curve::Circle(_) => CurveKind::Circle,
            Curve::Ellipse(_) => CurveKind::Ellipse,
            Curve::Hyperbola(_) => CurveKind::Hyperbola,
            Curve::Parabola(_) => CurveKind::Parabola,
            Curve::Bezier(_) => CurveKind::Bezier,
            Curve::BSpline(_) => CurveKind::BSpline,
            Curve::Offset(_) => CurveKind::Offset,
            Curve::Trimmed(_) => unreachable!("a view keeps a trimmed curve's basis"),
        }
    }

    /// The continuity of the curve inside the open range.
    pub(crate) fn continuity(&self) -> Continuity {
        match &*self.curve {
            Curve::BSpline(b) => {
                let p = b.degree();
                let knots = b.knots();
                let mults = b.multiplicities();
                let nk = knots.len();
                let tol = crate::tolerance::PARAMETRIC;
                let inside = |k: f64, lo: f64, hi: f64| k > lo + tol && k < hi - tol;
                let mut highest: Option<usize> = None;
                if b.is_periodic() {
                    let (k0, kn) = (knots[0], knots[nk - 1]);
                    let period = kn - k0;
                    let lo = in_period(self.first, k0, kn);
                    let hi = lo + (self.last - self.first);
                    let mut first_after: Option<(f64, usize)> = None;
                    // Knots of every period the range crosses (one period's knots: all but the
                    // last, which is the next period's first).
                    let periods = if period > 0.0 { ((hi - k0) / period).floor().max(0.0) as i64 + 1 } else { 1 };
                    for j in 0..periods.min(1_000_000) {
                        for i in 0..nk - 1 {
                            let k = knots[i] + j as f64 * period;
                            if inside(k, lo, hi) {
                                highest = Some(highest.map_or(mults[i], |m| m.max(mults[i])));
                            }
                            if k > lo + tol && first_after.is_none_or(|(f, _)| k < f) {
                                first_after = Some((k, mults[i]));
                            }
                        }
                    }
                    if highest.is_none() {
                        highest = first_after.map(|(_, m)| m).or(Some(mults[0]));
                    }
                } else {
                    for i in 1..nk.saturating_sub(1) {
                        if inside(knots[i], self.first, self.last) {
                            highest = Some(highest.map_or(mults[i], |m| m.max(mults[i])));
                        }
                    }
                }
                match highest {
                    None => Continuity::CN,
                    Some(m) => match p as i64 - m as i64 {
                        r if r <= 0 => Continuity::C0,
                        1 => Continuity::C1,
                        2 => Continuity::C2,
                        3 => Continuity::C3,
                        _ => Continuity::CN,
                    },
                }
            }
            Curve::Offset(o) => match o.basis().continuity() {
                Continuity::CN => Continuity::CN,
                Continuity::C3 => Continuity::C2,
                Continuity::C2 => Continuity::C1,
                Continuity::C1 => Continuity::C0,
                Continuity::G1 => Continuity::G1,
                Continuity::G2 => Continuity::G2,
                Continuity::C0 => Continuity::C0,
            },
            _ => Continuity::CN,
        }
    }

    /// The basis curve's own resolution.
    pub(crate) fn resolution(&self, r3d: f64) -> f64 {
        self.curve.resolution(r3d)
    }

    pub(crate) fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// The derivatives of orders `0..=order` at `u` as coordinates, on the polynomial piece the
    /// rules for B-splines choose (see `CurveAdaptor for GeomAdaptor`); `None` where the curve's
    /// own evaluation applies.
    fn piece_derivatives(&self, u: f64, order: usize) -> Option<Vec<[f64; 3]>> {
        let Curve::BSpline(b) = &*self.curve else { return None };
        let knots = b.knots();
        let nk = knots.len();
        if nk < 2 {
            return None;
        }
        let spans = nk - 1;
        let (k0, kn) = (knots[0], knots[nk - 1]);
        // (A bound counts as a knot it misses by rounding only.)
        let near = |a: f64, k: f64| (a - k).abs() <= 1e-13 * k.abs().max(1.0);
        let at_first = u == self.first;
        if at_first || u == self.last {
            if b.is_periodic() && (u < k0 || u > kn) {
                return None;
            }
            let span = if at_first {
                // The span starting at or after the bound.
                let i = knots.partition_point(|&k| k <= u);
                let start = if i < nk && near(u, knots[i]) { i } else { i.saturating_sub(1) };
                start.min(spans - 1)
            } else {
                // The span ending at or before the bound.
                let i = knots.partition_point(|&k| k < u);
                let end = if i < nk && near(u, knots[i]) {
                    i
                } else if i > 0 && near(u, knots[i - 1]) {
                    i - 1
                } else {
                    i
                };
                end.clamp(1, spans) - 1
            };
            return Some(b.span_derivatives(u, span, order));
        }
        // At a knot, or just before it by a rounding step: the piece after it.
        let w = if b.is_periodic() { in_period(u, k0, kn) } else { u };
        let i = knots.partition_point(|&k| k < w);
        if i >= 1 && i + 1 < nk {
            let k = knots[i];
            if k - w <= epsilon(kn.abs().min(w.abs())) {
                return Some(b.span_derivatives(w, i, order));
            }
        }
        None
    }
}

fn point_of(c: [f64; 3]) -> Point {
    Point::new(c[0], c[1], c[2])
}

fn vector_of(c: [f64; 3]) -> Vector {
    Vector::new(c[0], c[1], c[2])
}

impl CurveAdaptor for GeomAdaptor {
    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let (first, last) = (self.first, self.last);
        match &*self.curve {
            Curve::BSpline(b) => {
                if continuity == Continuity::C0
                    || (!b.is_periodic() && continuity <= self.continuity())
                    || matches!(continuity, Continuity::G1 | Continuity::G2)
                {
                    return vec![first, last];
                }
                let p = b.degree();
                let order = match continuity {
                    Continuity::C1 => 1,
                    Continuity::C2 => 2,
                    Continuity::C3 => 3,
                    _ => p,
                };
                let knots = b.knots();
                let mults = b.multiplicities();
                let nk = knots.len();
                let tol = self.resolution(1e-7).min(1e-9);
                let breaks_at = |m: usize| m + order > p;
                let mut out = vec![first];
                let mut push = |k: f64| {
                    if k > first + tol && k < last - tol {
                        out.push(k);
                    }
                };
                if b.is_periodic() {
                    let (k0, kn) = (knots[0], knots[nk - 1]);
                    let period = kn - k0;
                    if period > 0.0 {
                        let j0 = ((first - k0) / period).floor() as i64 - 1;
                        let j1 = ((last - k0) / period).floor() as i64 + 1;
                        let mut all = Vec::new();
                        for j in j0..=j1.min(j0 + 1_000_000) {
                            for i in 0..nk - 1 {
                                if breaks_at(mults[i]) {
                                    all.push(knots[i] + j as f64 * period);
                                }
                            }
                        }
                        all.sort_by(f64::total_cmp);
                        all.into_iter().for_each(&mut push);
                    }
                } else {
                    for i in 1..nk.saturating_sub(1) {
                        if breaks_at(mults[i]) {
                            push(knots[i]);
                        }
                    }
                }
                out.push(last);
                out
            }
            Curve::Offset(o) => {
                let higher = match continuity {
                    Continuity::C0 => Continuity::C1,
                    Continuity::C1 => Continuity::C2,
                    Continuity::C2 => Continuity::C3,
                    _ => Continuity::CN,
                };
                let basis = GeomAdaptor::new(o.basis(), first, last);
                let mut out = vec![first];
                out.extend(basis.intervals(higher).into_iter().filter(|&t| t > first && t < last));
                out.push(last);
                out
            }
            _ => vec![first, last],
        }
    }

    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.curve.period().unwrap_or(0.0)
    }

    fn value(&self, u: f64) -> Point {
        match self.piece_derivatives(u, 0) {
            Some(d) => point_of(d[0]),
            None => self.curve.point(u),
        }
    }

    fn d1(&self, u: f64) -> (Point, Vector) {
        match self.piece_derivatives(u, 1) {
            Some(d) => (point_of(d[0]), vector_of(d[1])),
            None => self.curve.d1(u),
        }
    }

    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        match self.piece_derivatives(u, 2) {
            Some(d) => (point_of(d[0]), vector_of(d[1]), vector_of(d[2])),
            None => self.curve.d2(u),
        }
    }

    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
        match self.piece_derivatives(u, 3) {
            Some(d) => (point_of(d[0]), vector_of(d[1]), vector_of(d[2]), vector_of(d[3])),
            None => self.curve.d3(u),
        }
    }

    fn dn(&self, u: f64, n: u32) -> Vector {
        if let Some(d) = self.piece_derivatives(u, n as usize) {
            return vector_of(d[n as usize]);
        }
        if let (Curve::Offset(o), 3) = (&*self.curve, n) {
            if let Some(v) = crate::geometry::eval::offset_adaptor_third_derivative(o, u) {
                return v;
            }
        }
        self.curve.dn(u, n)
    }

    fn resolution(&self, r3d: f64) -> f64 {
        GeomAdaptor::resolution(self, r3d)
    }

    fn curve_type(&self) -> Option<CurveKind> {
        Some(self.kind())
    }

    fn basis_curve(&self) -> Option<&Curve> {
        Some(&self.curve)
    }

    fn trimmed(&self, first: f64, last: f64, _tolerance: f64) -> Box<dyn CurveAdaptor> {
        Box::new(self.trim(first, last))
    }
}

/// The gap between `x` and the next double away from zero (the smallest subnormal at zero).
pub(crate) fn epsilon(x: f64) -> f64 {
    if x < 0.0 { x - x.next_down() } else { x.next_up() - x }
}

/// Whether `v` stands for an infinite value (at least half of `tolerance::INFINITE`).
fn infinite(v: f64) -> bool {
    v.abs() >= 0.5 * crate::tolerance::INFINITE
}

/// `u` shifted by whole periods `last − first` into `[first, last)`; unchanged for infinite
/// values or a degenerate period.
pub(crate) fn in_period(u: f64, first: f64, last: f64) -> f64 {
    let period = last - first;
    if infinite(u) || infinite(first) || infinite(last) || period < epsilon(last) {
        return u;
    }
    let k = ((first - u) / period).ceil();
    first.max(u + period * k)
}

/// Puts `u1` into the period starting at `first` (just below its end counting as before the
/// start) and `u2` into the period following `u1` (one that would end within `preci` of `u1`
/// going on a whole period). Returns the bounds themselves for an infinite or degenerate
/// period.
pub(crate) fn adjust_periodic(first: f64, last: f64, preci: f64, u1: f64, u2: f64) -> (f64, f64) {
    let period = last - first;
    if infinite(first) || infinite(last) || period < epsilon(last) {
        return (first, last);
    }
    let mut a = u1 - ((u1 - first) / period).floor() * period;
    if last - a < preci {
        a -= period;
    }
    let mut b = u2 - ((u2 - a) / period).floor() * period;
    if b - a < preci {
        b += period;
    }
    (a, b)
}

/// The 1-based interval of the non-decreasing table `values` containing `x`: 0 below the
/// table, `len + 1` above it, otherwise the largest `i ≤ len − 1` with `values[i − 1] < x`
/// (1 if there is none).
pub(crate) fn hunt(values: &[f64], x: f64) -> usize {
    let n = values.len();
    assert!(n > 0, "hunt in an empty table");
    if x < values[0] {
        return 0;
    }
    if x > values[n - 1] {
        return n + 1;
    }
    let below = values.partition_point(|&v| v < x);
    below.min(n - 1).max(1)
}

/// The 1-based index of the first knot of a B-spline's useful range: the smallest `i` with
/// `mults[0] + … + mults[i − 1] > degree`.
pub(super) fn first_u_knot_index(degree: usize, mults: &[usize]) -> usize {
    let mut sum = 0;
    for (i, &m) in mults.iter().enumerate() {
        sum += m;
        if sum > degree {
            return i + 1;
        }
    }
    panic!("multiplicities never exceed the degree")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::BSplineCurve;
    use std::f64::consts::{PI, TAU};

    #[test]
    fn spacing_of_doubles() {
        assert_eq!(epsilon(1.0), 2f64.powi(-52));
        assert_eq!(epsilon(-1.0), 2f64.powi(-52));
        assert_eq!(epsilon(2.0), 2f64.powi(-51));
        assert_eq!(epsilon(0.0), 5e-324);
        assert_eq!(epsilon(-0.0), 5e-324);
    }

    #[test]
    fn parameters_into_a_period() {
        assert_eq!(in_period(7.0, 0.0, TAU), 7.0 - TAU);
        assert_eq!(in_period(-0.5, 0.0, TAU), TAU - 0.5);
        for u in [0.0, 1.0, 3.0, 6.0] {
            assert_eq!(in_period(u, 0.0, TAU), u);
        }
        assert_eq!(in_period(5.0, -1e100, 1e100), 5.0);
        assert_eq!(in_period(5.0, 1.0, 1.0), 5.0);
        let (a, b) = adjust_periodic(0.0, TAU, 1e-9, TAU - 1e-12, 0.5);
        assert!((a + 1e-12).abs() < 1e-15, "{a}");
        assert_eq!(b, 0.5);
        assert_eq!(adjust_periodic(0.0, 1e100, 1e-9, 3.0, 4.0), (0.0, 1e100));
        // A range ending on its start goes on a whole period.
        let (a, b) = adjust_periodic(0.0, TAU, 1e-9, 1.0, 1.0);
        assert_eq!((a, b), (1.0, 1.0 + TAU));
    }

    #[test]
    fn hunting_in_a_table() {
        let t = [0.0, 1.0, 2.0, 3.0];
        let got: Vec<usize> = [-1.0, 0.0, 0.5, 1.0, 2.5, 3.0, 4.0].iter().map(|&x| hunt(&t, x)).collect();
        assert_eq!(got, [0, 1, 1, 1, 3, 3, 5]);
        assert_eq!(hunt(&[2.0], 2.0), 1);
    }

    #[test]
    fn first_knot_of_the_range() {
        assert_eq!(first_u_knot_index(3, &[4, 1, 1, 4]), 1);
        assert_eq!(first_u_knot_index(2, &[1, 1, 1, 1, 1, 1]), 3);
    }

    /// A degree-2 B-spline with a C0 knot at 0.5 (a corner).
    fn corner() -> Curve {
        let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 2.0, 0.0), Point::new(4.0, 0.0, 0.0)];
        Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], 2))
    }

    #[test]
    fn pieces_at_a_corner() {
        let c = corner();
        let Curve::BSpline(b) = &c else { unreachable!() };
        let left = b.span_derivatives(0.5, 0, 1)[1];
        let right = b.span_derivatives(0.5, 1, 1)[1];
        assert!((left[1] - right[1]).abs() > 1.0, "a corner: {left:?} {right:?}");
        let v = |a: Vector| [a.x, a.y, a.z];
        assert_eq!(v(GeomAdaptor::new(&c, 0.0, 0.5).d1(0.5).1), left);
        assert_eq!(v(GeomAdaptor::new(&c, 0.5, 1.0).d1(0.5).1), right);
        let just_below = 0.5f64.next_down();
        let d = GeomAdaptor::new(&c, 0.0, 1.0).d1(just_below).1;
        let expected = b.span_derivatives(just_below, 1, 1)[1];
        assert_eq!(v(d), expected);
        // Further below, the left piece.
        let d = GeomAdaptor::new(&c, 0.0, 1.0).d1(0.49).1;
        assert!((d.y - b.span_derivatives(0.49, 0, 1)[1][1]).abs() < 1e-12);
    }

    #[test]
    fn breakpoints_and_continuity() {
        let c = corner();
        let whole = GeomAdaptor::new(&c, 0.0, 1.0);
        assert_eq!(whole.intervals(Continuity::C1), [0.0, 0.5, 1.0]);
        assert_eq!(whole.intervals(Continuity::C0), [0.0, 1.0]);
        assert_eq!(whole.continuity(), Continuity::C0);
        let half = GeomAdaptor::new(&c, 0.5, 1.0);
        assert_eq!(half.intervals(Continuity::C1), [0.5, 1.0]);
        assert_eq!(GeomAdaptor::new(&c, 0.0, 0.5).continuity(), Continuity::CN);
        assert_eq!(whole.kind(), CurveKind::BSpline);
        assert_eq!(whole.curve_type(), Some(CurveKind::BSpline));
    }

    #[test]
    fn periodic_breakpoints_over_two_periods() {
        let poles: Vec<Point> = (0..6).map(|i| {
            let a = TAU * i as f64 / 6.0;
            Point::new(a.cos(), a.sin(), 0.0)
        }).collect();
        let knots: Vec<f64> = (0..=6).map(|i| i as f64).collect();
        let b = BSplineCurve::new_periodic(poles, knots, vec![1; 7], 3);
        let c = Curve::BSpline(b);
        let view = GeomAdaptor::new(&c, 0.0, 12.0);
        let got = view.intervals(Continuity::CN);
        let expected: Vec<f64> = (0..=12).map(|i| i as f64).collect();
        assert_eq!(got, expected);
        assert!(view.is_periodic());
        assert_eq!(view.period(), 6.0);
        // Within the knot period, as the curve itself.
        let p = view.value(1.3);
        let q = c.point(1.3);
        assert!((p - q).length() < 1e-14);
        let _ = PI;
    }

    #[test]
    fn trimmed_curves_unwrap() {
        let line = Curve::Line(crate::geometry::curve::Line::new(Point::ORIGIN, crate::math::Direction::new(1.0, 0.0, 0.0)));
        let a = GeomAdaptor::new(&line, -1.0, 2.0);
        assert_eq!(a.intervals(Continuity::C2), [-1.0, 2.0]);
        assert_eq!(a.kind(), CurveKind::Line);
        let t = a.trimmed(0.0, 1.0, 1e-7);
        assert_eq!((t.first_parameter(), t.last_parameter()), (0.0, 1.0));
        assert_eq!(t.period(), 0.0);
    }
}
