//! The part of OCCT's `Adaptor3d_Curve` interface that the algorithms ported for any curve
//! (OCCT's `GeomAbs_OtherCurve` paths: `Extrema_GGExtPC`, `GCPnts_TangentialDeflection`)
//! call, for curves that are not one of kentron's own [`Curve`](super::curve::Curve)s (such
//! as `GeomFill_SnglrFunc`, the cross product of a curve's first two derivatives).

use crate::math::{Point, Vector};

use super::Continuity;

/// A curve as OCCT's adaptors present one to the algorithms that take any curve.
#[doc(alias = "Adaptor3d_Curve")]
pub(crate) trait CurveAdaptor {
    fn first_parameter(&self) -> f64;

    fn last_parameter(&self) -> f64;

    /// `NbIntervals(S)` and `Intervals(T, S)`: the bounds of the intervals over which the
    /// curve is of continuity `continuity`, from the first parameter to the last (two values
    /// at least).
    fn intervals(&self, continuity: Continuity) -> Vec<f64>;

    fn is_periodic(&self) -> bool;

    /// `Period()` (only called for periodic curves).
    fn period(&self) -> f64;

    /// `Value`, `D0`.
    fn value(&self, u: f64) -> Point;

    fn d1(&self, u: f64) -> (Point, Vector);

    fn d2(&self, u: f64) -> (Point, Vector, Vector);

    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector);

    /// `DN(U, N)`, `n` at least 1.
    fn dn(&self, u: f64, n: u32) -> Vector;

    /// `Resolution(R3d)`: the parameter step that moves a point along the curve by at most
    /// `r3d`.
    fn resolution(&self, r3d: f64) -> f64;

    /// `GetType()`: the kind of curve, `None` for `GeomAbs_OtherCurve`.
    fn curve_type(&self) -> Option<CurveKind> {
        None
    }

    /// The curve itself, for the kinds `curve_type` names (`Line()`, `Circle()`, `BSpline()`,
    /// …).
    fn basis_curve(&self) -> Option<&Curve> {
        None
    }

    /// `Trim(First, Last, Tol)`: the same curve over `[first, last]`. Panics by default, as
    /// the base `Adaptor3d_Curve::Trim` raises `Standard_NotImplemented`.
    fn trimmed(&self, _first: f64, _last: f64, _tolerance: f64) -> Box<dyn CurveAdaptor> {
        panic!("Adaptor3d_Curve::Trim: Standard_NotImplemented")
    }
}

use std::sync::Arc;

use super::curve::Curve;
use super::ParametricCurve;

/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1e-9;
/// `Precision::Confusion()`.
const CONFUSION: f64 = 1e-7;

/// OCCT's `Epsilon(x)`: the gap from `x` to the next double away from zero.
pub(crate) fn epsilon(x: f64) -> f64 {
    if x >= 0.0 { x.next_up() - x } else { x - x.next_down() }
}

/// `ElCLib::InPeriod(u, first, last)`: `u` brought into `[first, last]` by whole periods.
pub(crate) fn in_period(u: f64, first: f64, last: f64) -> f64 {
    // (Precision::IsInfinite.)
    let infinite = |x: f64| x.abs() >= 0.5 * 2e100;
    if infinite(u) || infinite(first) || infinite(last) {
        return u;
    }
    let period = last - first;
    if period < epsilon(last) {
        return u;
    }
    first.max(u + period * ((first - u) / period).ceil())
}

/// `ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)`: `u1` into `[first, last)`, `u2`
/// into `(u1, u1 + period]`.
pub(crate) fn adjust_periodic(first: f64, last: f64, preci: f64, u1: f64, u2: f64) -> (f64, f64) {
    // (Precision::IsInfinite.)
    let infinite = |x: f64| x.abs() >= 0.5 * 2e100;
    if infinite(first) || infinite(last) {
        return (first, last);
    }
    let period = last - first;
    if period < epsilon(last) {
        return (first, last);
    }
    let mut u1 = u1 - ((u1 - first) / period).floor() * period;
    if last - u1 < preci {
        u1 -= period;
    }
    let mut u2 = u2 - ((u2 - u1) / period).floor() * period;
    if u2 - u1 < preci {
        u2 += period;
    }
    (u1, u2)
}

/// `BSplCLib::Hunt`: 1-based, `0` below the first value and `len + 1` above the last.
pub(crate) fn hunt(values: &[f64], x: f64) -> usize {
    if values[0] > x {
        return 0;
    }
    if values[values.len() - 1] < x {
        return values.len() + 1;
    }
    let mut position = 1;
    if values.len() <= 1 {
        return position;
    }
    let mut high = values.len();
    while high - position != 1 {
        let middle = (high + position) / 2;
        if values[middle - 1] < x {
            position = middle;
        } else {
            high = middle;
        }
    }
    position
}

/// `BSplCLib::LocateParameter(Knots, U, IsPeriodic, FromK1, ToK2, KnotIndex, NewU, UFirst,
/// ULast)`, 1-based.
fn locate_parameter(knots: &[f64], u: f64, periodic: bool, from_k1: usize, to_k2: usize, u_first: f64, u_last: f64) -> (usize, f64) {
    let (first, last) = if from_k1 < to_k2 { (from_k1, to_k2) } else { (to_k2, from_k1) };
    let last1 = last - 1;
    let mut new_u = u;
    if periodic && (new_u < u_first || new_u > u_last) {
        new_u = in_period(new_u, u_first, u_last);
    }
    let mut index = hunt(knots, new_u);
    let at = |i: usize| knots[i - 1];
    let eps = epsilon(at(knots.len()).abs().min(u.abs()));
    if index < knots.len() && (new_u - at(index + 1)).abs() <= eps {
        index += 1;
    }
    index = index.max(first);
    if index > last1 {
        index = last1;
    }
    if index != last1 {
        let (mut k1, mut k2) = (at(index), at(index + 1));
        while (k2 - k1).abs() <= eps {
            index += 1;
            if index >= knots.len() {
                break;
            }
            k1 = k2;
            k2 = at(index + 1);
        }
    }
    (index, new_u)
}

/// `BSplCLib::FirstUKnotIndex`, 1-based.
fn first_u_knot_index(degree: usize, mults: &[usize]) -> usize {
    let mut index = 1;
    let mut sigma = mults[0];
    while sigma <= degree {
        index += 1;
        sigma += mults[index - 1];
    }
    index
}

/// `BSplCLib::LastUKnotIndex`, 1-based.
fn last_u_knot_index(degree: usize, mults: &[usize]) -> usize {
    let mut index = mults.len();
    let mut sigma = mults[index - 1];
    while sigma <= degree {
        index -= 1;
        sigma += mults[index - 1];
    }
    index
}

/// `BSplCLib::Intervals`: the bounds of the intervals of `[first, last]` where the B-spline
/// is `continuity` times differentiable.
#[allow(clippy::too_many_arguments)]
fn bspline_intervals(knots: &[f64], mults: &[usize], degree: usize, periodic: bool, continuity: usize, first: f64, last: f64, tolerance: f64) -> Vec<f64> {
    let first_index = if periodic { 1 } else { first_u_knot_index(degree, mults) };
    let last_index = if periodic { knots.len() } else { last_u_knot_index(degree, mults) };
    let mut new_knots = Vec::with_capacity(last_index - first_index + 1);
    for index in first_index..=last_index {
        if mults[index - 1] as isize > degree as isize - continuity as isize || index == first_index || index == last_index {
            new_knots.push(knots[index - 1]);
        }
    }
    let n = new_knots.len();
    let (mut cur_first, mut cur_last) = (first, last);
    let mut period = 0.0;
    let (mut first_period, mut last_period) = (0i64, 0i64);
    if periodic {
        let (lower, upper) = (knots[0], knots[knots.len() - 1]);
        period = upper - lower;
        while cur_first < lower {
            cur_first += period;
            first_period -= 1;
        }
        while cur_last < lower {
            cur_last += period;
            last_period -= 1;
        }
        while cur_first >= upper {
            cur_first -= period;
            first_period += 1;
        }
        while cur_last >= upper {
            cur_last -= period;
            last_period += 1;
        }
    }
    // (The overload without multiplicities: not periodic, from 1 to the count.)
    let (mut index1, _) = locate_parameter(&new_knots, cur_first, false, 1, n, 0.0, 1.0);
    let (mut index2, _) = locate_parameter(&new_knots, cur_last, false, 1, n, 0.0, 1.0);
    let at = |i: usize| new_knots[i - 1];
    if index1 < n && (at(index1 + 1) - cur_first).abs() < tolerance {
        index1 += 1;
    }
    if n > 0 && (at(index2) - cur_last).abs() < tolerance {
        index2 -= 1;
    }
    let count = (index2 as i64 - index1 as i64 + 1 + (last_period - first_period) * (n as i64 - 1)) as usize;
    let mut out = vec![0.0; count + 1];
    if periodic && last_period != first_period {
        let mut k = 0;
        for i in index1..n {
            out[k] = at(i) + first_period as f64 * period;
            k += 1;
        }
        for period_number in first_period + 1..last_period {
            for i in 1..n {
                out[k] = at(i) + period_number as f64 * period;
                k += 1;
            }
        }
        for i in 1..=index2 {
            out[k] = at(i) + last_period as f64 * period;
            k += 1;
        }
    } else {
        for (k, i) in (index1..=index2).enumerate() {
            out[k] = at(i) + first_period as f64 * period;
        }
    }
    out[0] = first;
    out[count] = last;
    out
}

/// The kinds of curve `GeomAdaptor_Curve::GetType` tells apart, for kentron's curves.
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

/// OCCT's `GeomAdaptor_Curve`: a curve over `[first, last]` (a trimmed curve's basis, as
/// `Load` takes it), evaluating a B-spline exactly at either bound on the span inside.
#[doc(alias = "GeomAdaptor_Curve")]
#[derive(Clone, Debug)]
pub(crate) struct GeomAdaptor {
    curve: Arc<Curve>,
    first: f64,
    last: f64,
}

impl GeomAdaptor {
    /// `GeomAdaptor_Curve(C, First, Last)` (a trimmed curve loaded as its basis).
    pub(crate) fn new(curve: &Curve, first: f64, last: f64) -> GeomAdaptor {
        let mut curve = curve;
        while let Curve::Trimmed(t) = curve {
            curve = t.basis();
        }
        GeomAdaptor { curve: Arc::new(curve.clone()), first, last }
    }

    pub(crate) fn curve(&self) -> &Curve {
        &self.curve
    }

    /// `Trim(First, Last)`: the same curve over `[first, last]`.
    pub(crate) fn trim(&self, first: f64, last: f64) -> GeomAdaptor {
        GeomAdaptor { curve: self.curve.clone(), first, last }
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
            Curve::Trimmed(_) => unreachable!("trimmed curves are loaded as their basis"),
        }
    }

    /// `LocalContinuity(U1, U2)` of a B-spline.
    fn local_continuity(&self, u1: f64, u2: f64) -> Continuity {
        let Curve::BSpline(b) = &*self.curve else { unreachable!() };
        let (knots, mults) = (b.knots(), b.multiplicities());
        let nb = knots.len();
        let periodic = b.is_periodic();
        let (uf, ul) = if periodic { (knots[0], knots[nb - 1]) } else { (0.0, 1.0) };
        let (mut index1, new_first) = locate_parameter(knots, u1, periodic, 1, nb, uf, ul);
        let (mut index2, new_last) = locate_parameter(knots, u2, periodic, 1, nb, uf, ul);
        if (new_first - knots[index1]).abs() < P_CONFUSION && index1 < nb {
            index1 += 1;
        }
        if (new_last - knots[index2 - 1]).abs() < P_CONFUSION {
            index2 -= 1;
        }
        if periodic && index1 == nb {
            index1 = 1;
        }
        let mult_max = if index2 as isize - index1 as isize <= 0 && !periodic {
            100
        } else {
            let mut m = mults[index1];
            for i in index1 + 1..=index2 {
                m = m.max(mults[i - 1]);
            }
            b.degree() as isize - m as isize
        };
        match mult_max {
            m if m <= 0 => Continuity::C0,
            1 => Continuity::C1,
            2 => Continuity::C2,
            3 => Continuity::C3,
            _ => Continuity::CN,
        }
    }

    /// `Continuity()`.
    pub(crate) fn continuity(&self) -> Continuity {
        match &*self.curve {
            Curve::BSpline(_) => self.local_continuity(self.first, self.last),
            Curve::Offset(o) => match o.basis().continuity() {
                Continuity::CN => Continuity::CN,
                Continuity::C3 => Continuity::C2,
                Continuity::C2 => Continuity::C1,
                Continuity::C1 => Continuity::C0,
                Continuity::G1 => Continuity::G1,
                Continuity::G2 => Continuity::G2,
                // (OCCT throws Standard_NoSuchObject.)
                Continuity::C0 => Continuity::C0,
            },
            _ => Continuity::CN,
        }
    }

    /// `Resolution(R3d)` (kentron's B-spline and other curves' own).
    pub(crate) fn resolution(&self, r3d: f64) -> f64 {
        self.curve.resolution(r3d)
    }

    /// For a B-spline evaluated inside, the span OCCT evaluates on where kentron's own
    /// evaluation would take another: a parameter short of an interior knot by no more than
    /// `Epsilon` is taken to the knot (`BSplCLib::LocateParameter`, which places both the
    /// adaptor's span cache, `BSplCLib_CacheParams`, and `Geom_BSplineCurve::DN`), and
    /// evaluated on the span after it. The knots (0-based) of that span.
    fn snapped_span(&self, u: f64) -> Option<(usize, usize)> {
        let Curve::BSpline(b) = &*self.curve else { return None };
        let knots = b.knots();
        let n = knots.len();
        let u = if b.is_periodic() { cache_normalization(u, knots[0], knots[n - 1]) } else { u };
        let eps = epsilon(knots[n - 1].abs().min(u.abs()));
        let i = knots.partition_point(|&k| k <= u);
        // (The next knot just past `u`, an interior one.)
        (i > 0 && i + 1 < n && knots[i] - u <= eps).then_some((i, i + 1))
    }

    /// `IsBoundary`: for a B-spline evaluated exactly at a bound, the knots (0-based) of
    /// the span inside.
    fn boundary_span(&self, u: f64) -> Option<(usize, usize)> {
        let Curve::BSpline(b) = &*self.curve else { return None };
        if u != self.first && u != self.last {
            return None;
        }
        let nb = b.knots().len();
        let (mut start, mut finish) = locate_u(b.knots(), b.is_periodic(), u, P_CONFUSION / 2.0);
        if u == self.first {
            start = start.max(1);
            if start >= finish {
                finish = start + 1;
            }
        } else {
            finish = finish.min(nb);
            if start >= finish {
                start = finish - 1;
            }
        }
        if start == 0 || finish > nb {
            // (A periodic curve's bound just past an end knot, brought round to the other
            // end: OCCT's `LocalD*` from knot 0 or to knot NbKnots + 1 then evaluate the end
            // span at the parameter brought into the period, as the curve's own evaluation
            // does.)
            return None;
        }
        Some((start - 1, finish - 1))
    }

    pub(crate) fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }
}

fn point(c: [f64; 3]) -> Point {
    Point::new(c[0], c[1], c[2])
}

fn vector(c: [f64; 3]) -> Vector {
    Vector::new(c[0], c[1], c[2])
}

/// `Geom_BSplineCurve::PeriodicNormalization` of a periodic B-spline over `[first, last]`:
/// `u` shifted by whole periods into it (`last` itself kept).
fn periodic_normalization(u: f64, first: f64, last: f64) -> f64 {
    if u >= first && u <= last {
        return u;
    }
    let period = last - first;
    let periods = ((u - first) / period).floor();
    let mut u = u;
    if periods != 0.0 {
        u -= periods * period;
    }
    if u > last {
        let shifted = u - period;
        u = if shifted == u { shifted.next_down() } else { shifted };
    } else if u < first {
        let shifted = u + period;
        u = if shifted == u { shifted.next_up() } else { shifted };
    }
    u
}

/// `BSplCLib_CacheParams::PeriodicNormalization`.
fn cache_normalization(u: f64, first: f64, last: f64) -> f64 {
    let period = last - first;
    if u < first {
        u + period * (((first - u) / period).trunc() + 1.0)
    } else if u > last {
        u - period * (((u - last) / period).trunc() + 1.0)
    } else {
        u
    }
}

/// `Geom_BSplineCurve::LocateU` over the knots (without repetition), 1-based.
fn locate_u(knots: &[f64], periodic: bool, u: f64, tolerance: f64) -> (usize, usize) {
    let n = knots.len();
    let (u_first, u_last) = (knots[0], knots[n - 1]);
    let new_u = if periodic { periodic_normalization(u, u_first, u_last) } else { u };
    let tolerance = tolerance.abs();
    if (new_u - u_first).abs() <= tolerance {
        (1, 1)
    } else if (new_u - u_last).abs() <= tolerance {
        (n, n)
    } else if new_u < u_first {
        (0, 1)
    } else if new_u > u_last {
        (n, n + 1)
    } else {
        let mut i1 = hunt(knots, new_u).clamp(1, n);
        while i1 < n && (knots[i1] - new_u).abs() <= tolerance {
            i1 += 1;
        }
        if (knots[i1 - 1] - new_u).abs() <= tolerance { (i1, i1) } else { (i1, i1 + 1) }
    }
}

impl CurveAdaptor for GeomAdaptor {
    fn curve_type(&self) -> Option<CurveKind> {
        Some(self.kind())
    }

    fn basis_curve(&self) -> Option<&Curve> {
        Some(&self.curve)
    }

    /// `GeomAdaptor_Curve::Trim` (the tolerance unused, as in OCCT).
    fn trimmed(&self, first: f64, last: f64, _tolerance: f64) -> Box<dyn CurveAdaptor> {
        Box::new(self.trim(first, last))
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `NbIntervals(S)` and `Intervals(T, S)`.
    fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        match &*self.curve {
            Curve::BSpline(b) => {
                if (!b.is_periodic() && continuity <= self.continuity()) || continuity == Continuity::C0 {
                    return vec![self.first, self.last];
                }
                let degree = b.degree();
                let order = match continuity {
                    Continuity::C1 => 1,
                    Continuity::C2 => 2,
                    Continuity::C3 => 3,
                    Continuity::CN => degree,
                    // (OCCT throws Standard_DomainError.)
                    _ => return vec![self.first, self.last],
                };
                let eps = self.resolution(CONFUSION).min(P_CONFUSION);
                bspline_intervals(b.knots(), b.multiplicities(), degree, b.is_periodic(), order, self.first, self.last, eps)
            }
            Curve::Offset(o) => {
                let base = match continuity {
                    Continuity::C0 => Continuity::C1,
                    Continuity::C1 => Continuity::C2,
                    Continuity::C2 => Continuity::C3,
                    _ => Continuity::CN,
                };
                let basis = GeomAdaptor::new(o.basis(), self.first, self.last);
                let inner = basis.intervals(base);
                let mut out = vec![self.first];
                if inner.len() > 2 {
                    out.extend(inner.iter().copied().filter(|&t| t > self.first && t < self.last));
                }
                out.push(self.last);
                out
            }
            _ => vec![self.first, self.last],
        }
    }

    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.curve.period().unwrap_or(0.0)
    }

    fn value(&self, u: f64) -> Point {
        match (self.boundary_span(u).or_else(|| self.snapped_span(u)), &*self.curve) {
            (Some((a, _)), Curve::BSpline(c)) => point(c.span_derivatives(u, a, 0)[0]),
            _ => self.curve.point(u),
        }
    }

    fn d1(&self, u: f64) -> (Point, Vector) {
        match (self.boundary_span(u).or_else(|| self.snapped_span(u)), &*self.curve) {
            (Some((a, _)), Curve::BSpline(c)) => {
                let d = c.span_derivatives(u, a, 1);
                (point(d[0]), vector(d[1]))
            }
            _ => self.curve.d1(u),
        }
    }

    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        match (self.boundary_span(u).or_else(|| self.snapped_span(u)), &*self.curve) {
            (Some((a, _)), Curve::BSpline(c)) => {
                let d = c.span_derivatives(u, a, 2);
                (point(d[0]), vector(d[1]), vector(d[2]))
            }
            // (Not the second derivative `d3` returns with the third: an offset curve's
            // differs, as OCCT's `EvalD3` does from its `EvalD2`.)
            _ => self.curve.d2(u),
        }
    }

    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
        match (self.boundary_span(u).or_else(|| self.snapped_span(u)), &*self.curve) {
            (Some((a, _)), Curve::BSpline(c)) => {
                let d = c.span_derivatives(u, a, 3);
                (point(d[0]), vector(d[1]), vector(d[2]), vector(d[3]))
            }
            _ => self.curve.d3(u),
        }
    }

    fn dn(&self, u: f64, n: u32) -> Vector {
        match (self.boundary_span(u).or_else(|| self.snapped_span(u)), &*self.curve) {
            (Some((a, _)), Curve::BSpline(c)) => vector(c.span_derivatives(u, a, n as usize)[n as usize]),
            // (OCCT throws where the third derivative is undefined; kentron's own then.)
            (_, Curve::Offset(c)) if n == 3 => super::eval::offset_adaptor_third_derivative(c, u).unwrap_or_else(|| self.curve.dn(u, n)),
            _ => self.curve.dn(u, n),
        }
    }

    fn resolution(&self, r3d: f64) -> f64 {
        GeomAdaptor::resolution(self, r3d)
    }
}
