//! The Frenet frame of a spine piece: tangent `T`, principal normal `N` and binormal `B`,
//! the basis of every pipe frame law.
//!
//! A piece is a curve travelled from `first` to `last` (backwards when `last < first`); the
//! frame belongs to the curve as travelled, so tangents point the way of travel. Where the
//! curvature vanishes the textbook formula `B = C′ × C″ / |C′ × C″|` has no value:
//!
//! - along lines and exactly straight stretches the binormal and normal are the default ones
//!   of the tangent ([`Frame::from_z`]: `N` its X axis, `B` its Y axis);
//! - at isolated zeros of `C′ × C″` (inflections of planar curves, and the rarer points where
//!   a space curve's curvature touches zero) the frame is the one-sided limit from the
//!   evaluation's side, computed from the first higher derivative that leaves the tangent;
//!   such points are found once at construction and split the curve's intervals, so that a
//!   frame asked for *at* one is the limit from the interval it is asked in;
//! - where `C′` itself vanishes the tangent comes from the first non-vanishing derivative.

use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::{Continuity, ParametricCurve};
use crate::math::{Direction, Frame, Point, Vector};

/// Parameters closer than this are one point (singular points, and their distance to
/// interval bounds).
const SAME_PARAMETER: f64 = 1e-9;

/// `C′ × C″` counts as zero below this multiple of `ε·|C′|·|C″|` (rounding level).
const ROUNDING: f64 = 64.0 * f64::EPSILON;

/// The highest derivative tried for the limit at a singular point.
const HIGHEST_DERIVATIVE: u32 = 8;

/// A curve view of the piece as travelled, with the range of its parameter: the curve
/// itself from `first` to `last`, or, backwards, a reversed copy over the reversed range.
pub(super) fn travelled(curve: &Curve, first: f64, last: f64) -> (GeomAdaptor, f64, f64) {
    if last >= first {
        return (GeomAdaptor::new(curve, first, last), first, last);
    }
    let mut reversed = curve.clone();
    reversed.reverse();
    let (a, b) = (curve.reversed_parameter(first), curve.reversed_parameter(last));
    (GeomAdaptor::new(&reversed, a, b), a, b)
}

/// The unit vector along `v`.
pub(super) fn unit(v: Vector) -> Vector {
    v * (1.0 / v.length())
}

/// `v` turned about the unit `axis` by `angle` (Rodrigues).
pub(super) fn rotate_about(v: Vector, axis: Vector, angle: f64) -> Vector {
    let (s, c) = angle.sin_cos();
    v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c))
}

/// `v` turned by the least rotation taking the unit vector `from` onto the unit vector
/// `to` (about `from × to`); unchanged where the two are parallel to the last bit.
pub(super) fn least_rotation(v: Vector, from: Vector, to: Vector) -> Vector {
    let axis = from.cross(to);
    let s = axis.length();
    if s == 0.0 {
        return v;
    }
    let c = from.dot(to);
    rotate_about(v, axis * (1.0 / s), s.atan2(c))
}

/// The angle turning `a` onto `b` about the unit `axis` (both perpendicular to it), in
/// `(−π, π]`.
pub(super) fn signed_angle(a: Vector, b: Vector, axis: Vector) -> f64 {
    axis.dot(a.cross(b)).atan2(a.dot(b))
}

/// The default normal and binormal of the unit tangent `t`: the X and Y axes of
/// [`Frame::from_z`].
pub(super) fn default_normal(t: Vector) -> (Vector, Vector) {
    let frame = Frame::from_z(Point::ORIGIN, Direction::from_vector(t));
    (Vector::from(frame.x()), Vector::from(frame.y()))
}

/// The unit tangent of `view` at `u`, along the first non-vanishing derivative where `C′`
/// vanishes (see [`null_derivative_tangent`]).
pub(super) fn unit_tangent(view: &GeomAdaptor, u: f64) -> Vector {
    let (_, d1) = view.d1(u);
    if d1.length() > 0.0 { unit(d1) } else { null_derivative_tangent(view, u) }
}

/// The tangent of `view` at `u` where `C′` vanishes: along the first non-vanishing higher
/// derivative (second or third), pointing forward (an acute angle with the chord to a point
/// a thousandth of the range ahead, or behind at the end); else along a one-sided
/// three-point difference of points.
pub(super) fn null_derivative_tangent(view: &GeomAdaptor, u: f64) -> Vector {
    let (a, b) = (view.first_parameter(), view.last_parameter());
    let h = 1e-3 * (b - a);
    let ahead = u + h <= b;
    let p = |x: f64| view.value(x);
    let chord = if ahead { p(u + h) - p(u) } else { p(u) - p(u - h) };
    for j in 2..=3 {
        let d = view.dn(u, j);
        if d.length() > 0.0 {
            let d = unit(d);
            return if d.dot(chord) < 0.0 { -d } else { d };
        }
    }
    let d = if ahead { (p(u + h) - p(u)) * 4.0 - (p(u + 2.0 * h) - p(u)) } else { (p(u) - p(u - h)) * 4.0 - (p(u) - p(u - 2.0 * h)) };
    unit(d)
}

/// A Frenet frame: unit tangent, normal and binormal.
#[derive(Clone, Copy, Debug)]
pub(super) struct Trihedron {
    pub t: Vector,
    pub n: Vector,
    pub b: Vector,
}

/// The Frenet frame along a travelled curve.
#[derive(Clone, Debug)]
pub(super) struct Frenet {
    curve: GeomAdaptor,
    /// The C2 breaks merged with the singular points, from the first parameter to the last.
    intervals: Vec<f64>,
    /// The interior singular points (zeros of `C′ × C″`).
    singular: Vec<f64>,
    /// A line: the default frame everywhere.
    line: bool,
    /// The size of `C″` on each interval (against which rounding is judged).
    second: Vec<f64>,
}

impl Frenet {
    /// The frame along `curve` (a travelled view, see [`travelled`]), its singular points
    /// searched for unless it is a line or a conic.
    pub(super) fn new(curve: GeomAdaptor) -> Frenet {
        Self::build(curve, true)
    }

    /// The frame without the search for singular points: for laws that ask only for
    /// tangents, or for the frame at the curve's ends.
    pub(super) fn without_singular_points(curve: GeomAdaptor) -> Frenet {
        Self::build(curve, false)
    }

    fn build(curve: GeomAdaptor, search: bool) -> Frenet {
        let kind = curve.kind();
        let line = kind == CurveKind::Line;
        let breaks = curve.intervals(Continuity::C2);
        let analytic = matches!(kind, CurveKind::Line | CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola);
        let mut singular = Vec::new();
        if search && !analytic {
            for w in breaks.windows(2) {
                singular.extend(singular_points(&curve, w[0], w[1]));
            }
        }
        let mut intervals = breaks.clone();
        intervals.extend(singular.iter().copied());
        intervals.sort_by(f64::total_cmp);
        let second = intervals.windows(2).map(|w| derivative_scale(&curve.trim(w[0], w[1]), 2)).collect();
        Frenet { curve, intervals, singular, line, second }
    }

    pub(super) fn curve(&self) -> &GeomAdaptor {
        &self.curve
    }

    pub(super) fn first(&self) -> f64 {
        self.curve.first_parameter()
    }

    pub(super) fn last(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// The bounds of the intervals the frame is smooth on.
    pub(super) fn intervals(&self) -> &[f64] {
        &self.intervals
    }

    /// The interval holding `u`: the first of two sharing it as a bound.
    pub(super) fn interval_from_left(&self, u: f64) -> usize {
        let n = self.intervals.len() - 1;
        self.intervals[1..].partition_point(|&b| b < u).min(n - 1)
    }

    /// The interval holding `u`: the one starting at it where it is a bound (the last at the
    /// end).
    pub(super) fn interval_from_right(&self, u: f64) -> usize {
        let n = self.intervals.len() - 1;
        self.intervals.partition_point(|&a| a <= u).clamp(1, n) - 1
    }

    /// The parameter at fraction `s` of the way.
    pub(super) fn at_fraction(&self, s: f64) -> f64 {
        let (a, b) = (self.first(), self.last());
        if s == 1.0 { b } else { a + (b - a) * s }
    }

    /// The frame at `u` on the whole curve: the curve's own side at interior knots (the span
    /// starting there), inwards at the ends.
    pub(super) fn frame(&self, u: f64) -> Trihedron {
        self.frame_in(u, self.interval_from_right(u))
    }

    /// The unit tangent at `u` on the whole curve (as [`Frenet::frame`]).
    pub(super) fn tangent(&self, u: f64) -> Vector {
        let i = self.interval_from_right(u);
        self.tangent_in(u, i)
    }

    /// The unit tangent at `u` seen from interval `i`.
    pub(super) fn tangent_in(&self, u: f64, i: usize) -> Vector {
        let view = self.view(i);
        let (_, d1) = view.d1(u);
        if d1.length() > 0.0 { unit(d1) } else { null_derivative_tangent(&view, u) }
    }

    /// The curve restricted to interval `i` (one-sided evaluation at its bounds).
    fn view(&self, i: usize) -> GeomAdaptor {
        if self.intervals.len() == 2 {
            return self.curve.clone();
        }
        self.curve.trim(self.intervals[i], self.intervals[i + 1])
    }

    /// The frame at `u` seen from interval `i` (one-sided at its bounds).
    pub(super) fn frame_in(&self, u: f64, i: usize) -> Trihedron {
        let view = self.view(i);
        let (_, d1, d2) = view.d2(u);
        if !(d1.length() > 0.0) {
            return self.null_derivative_frame(&view, u, i);
        }
        let t = unit(d1);
        if self.line {
            let (n, b) = default_normal(t);
            return Trihedron { t, n, b };
        }
        let x = d1.cross(d2);
        let (a, b) = (self.intervals[i], self.intervals[i + 1]);
        let near_singular = |bound: f64| (u - bound).abs() <= SAME_PARAMETER && self.singular.iter().any(|&s| s == bound);
        if x.length() > ROUNDING * d1.length() * d2.length().max(self.second[i]) && !near_singular(a) && !near_singular(b) {
            let bn = unit(x);
            return Trihedron { t, n: bn.cross(t), b: bn };
        }
        self.limit_frame(&view, u, t, i)
    }

    /// The frame at `u` where `C′ × C″` vanishes: the binormal along `T × C⁽ʲ⁾` for the
    /// first derivative leaving the tangent, signed as the regular frame a little way into
    /// the interval; the default frame where none does (a straight stretch).
    fn limit_frame(&self, view: &GeomAdaptor, u: f64, t: Vector, i: usize) -> Trihedron {
        let mut along = None;
        for j in 2..=HIGHEST_DERIVATIVE {
            let dj = view.dn(u, j);
            let v = t.cross(dj);
            if v.length() > ROUNDING * dj.length().max(derivative_scale(view, j)) && v.length() > 0.0 {
                along = Some(unit(v));
                break;
            }
        }
        let Some(mut bn) = along else {
            let (n, b) = default_normal(t);
            return Trihedron { t, n, b };
        };
        if let Some(inner) = self.inner_binormal(view, u, i) {
            if bn.dot(inner) < 0.0 {
                bn = -bn;
            }
        }
        Trihedron { t, n: bn.cross(t), b: bn }
    }

    /// The binormal of the nearest clearly regular point from `u` into interval `i`.
    fn inner_binormal(&self, view: &GeomAdaptor, u: f64, i: usize) -> Option<Vector> {
        let (a, b) = (self.intervals[i], self.intervals[i + 1]);
        let width = b - a;
        let direction = if u - a <= b - u { 1.0 } else { -1.0 };
        let mut step = 1e-7 * width;
        while step <= 0.5 * width {
            let (_, d1, d2) = view.d2(u + direction * step);
            let x = d1.cross(d2);
            if x.length() > 1e-8 * d1.length() * d2.length() {
                return Some(unit(x));
            }
            step *= 10.0;
        }
        None
    }

    /// The frame where `C′` vanishes: the tangent of [`null_derivative_tangent`], the
    /// normal and binormal of a nearby point inwards turned by the least rotation between
    /// the tangents.
    fn null_derivative_frame(&self, view: &GeomAdaptor, u: f64, i: usize) -> Trihedron {
        let t = null_derivative_tangent(view, u);
        let (a, b) = (self.intervals[i], self.intervals[i + 1]);
        let h = 1e-2 * (b - a);
        let near = if u - a <= b - u { u + h } else { u - h };
        let (_, d1) = view.d1(near);
        if !(d1.length() > 0.0) || h == 0.0 {
            let (n, b) = default_normal(t);
            return Trihedron { t, n, b };
        }
        let there = self.frame_in(near, i);
        let n = least_rotation(there.n, there.t, t);
        let n = unit(n - t * t.dot(n));
        Trihedron { t, n, b: t.cross(n) }
    }
}

/// The size of the `j`-th derivative over a view: the largest at its ends and middle.
fn derivative_scale(view: &GeomAdaptor, j: u32) -> f64 {
    let (a, b) = (view.first_parameter(), view.last_parameter());
    [a, 0.5 * (a + b), b].iter().map(|&u| view.dn(u, j).length()).fold(0.0, f64::max)
}

/// The interior zeros of `C′ × C″` on the C2 stretch `[a, b]` of `curve`, increasing.
///
/// `|C′ × C″|` is sampled densely (per knot span); at each local minimum of the samples the
/// cross product is projected on its change across the minimum, and a sign change of that
/// projection is solved for by bisection. A root where `|C′ × C″|` is not below `1e-9` of
/// its largest sampled value is not a zero (a space curve passing near one).
fn singular_points(curve: &GeomAdaptor, a: f64, b: f64) -> Vec<f64> {
    let view = curve.trim(a, b);
    let spans = view.intervals(Continuity::CN);
    let per_span = match curve.curve() {
        Curve::BSpline(s) => (4 * s.degree()).max(12),
        Curve::Bezier(s) => (4 * s.degree()).max(16),
        _ => 24,
    };
    let mut us = Vec::with_capacity(per_span * (spans.len() - 1) + 1);
    for w in spans.windows(2) {
        for k in 0..per_span {
            us.push(w[0] + (w[1] - w[0]) * k as f64 / per_span as f64);
        }
    }
    us.push(b);
    let cross = |u: f64| {
        let (_, d1, d2) = view.d2(u);
        d1.cross(d2)
    };
    let xs: Vec<Vector> = us.iter().map(|&u| cross(u)).collect();
    let rs: Vec<f64> = xs.iter().map(|x| x.length()).collect();
    let largest = rs.iter().copied().fold(0.0, f64::max);
    if !(largest > 0.0) {
        return Vec::new();
    }
    let mut found: Vec<f64> = Vec::new();
    for k in 1..us.len() - 1 {
        if !(rs[k] <= rs[k - 1] && rs[k] <= rs[k + 1]) {
            continue;
        }
        let root = if rs[k] == 0.0 {
            Some(us[k])
        } else {
            let along = xs[k + 1] - xs[k - 1];
            let f = |u: f64| cross(u).dot(along);
            bisect(&f, us[k - 1], us[k + 1])
        };
        let Some(s) = root else { continue };
        if cross(s).length() > 1e-9 * largest {
            continue;
        }
        if s - a <= SAME_PARAMETER || b - s <= SAME_PARAMETER || found.last().is_some_and(|&p| s - p <= SAME_PARAMETER) {
            continue;
        }
        found.push(s);
    }
    found
}

/// A sign change of `f` in `[lo, hi]`, by bisection to the last bit; `None` without one.
fn bisect(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> Option<f64> {
    let (mut lo, mut hi) = (lo, hi);
    let (mut flo, fhi) = (f(lo), f(hi));
    if flo == 0.0 {
        return Some(lo);
    }
    if fhi == 0.0 {
        return Some(hi);
    }
    if flo.signum() == fhi.signum() {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        let fm = f(mid);
        if fm == 0.0 {
            return Some(mid);
        }
        if fm.signum() == flo.signum() {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    Some(if flo.abs() <= f(hi).abs() { lo } else { hi })
}

/// The Frenet frame law of one spine piece (a pipe's Frenet mode).
#[derive(Clone, Debug)]
pub(super) struct FrenetLaw {
    frenet: Frenet,
}

impl FrenetLaw {
    /// The law along `curve` travelled from `first` to `last`.
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> FrenetLaw {
        let (view, _, _) = travelled(curve, first, last);
        FrenetLaw { frenet: Frenet::new(view) }
    }

    /// The tangent and normal at fraction `fraction` of the way, seen from the interval
    /// holding it (the one before where it is a bound).
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let u = self.frenet.at_fraction(fraction);
        let f = self.frenet.frame_in(u, self.frenet.interval_from_left(u));
        (f.t, f.n)
    }
}

#[cfg(test)]
#[path = "frenet_tests.rs"]
mod tests;
