//! Distance extrema: the points where the distance between two objects is locally
//! minimal or maximal (where the connecting segment is perpendicular to both).

use std::fmt;
use std::ops::Index;

use crate::math::{Point, Point2};

use super::curve::Curve;
use super::control::ControlPoint;
use super::curve2::Curve2;
use super::parametric::ParametricCurve;
use super::surface::Surface;

/// One extremum between two curves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveExtremum<P = Point> {
    pub distance: f64,
    /// Parameter and point on the first curve.
    pub t1: f64,
    pub p1: P,
    /// Parameter and point on the second curve.
    pub t2: f64,
    pub p2: P,
}

/// The extrema between two curves.
#[doc(alias = "GeomAPI_ExtremaCurveCurve")]
#[doc(alias = "Extrema_ExtCC")]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveCurveExtrema<P = Point> {
    /// Isolated extrema (where the connecting segment is perpendicular to both curves, and
    /// for parallel pieces that only touch, the touching pair), sorted by increasing
    /// distance.
    Isolated(Vec<CurveExtremum<P>>),
    /// The curves are parallel (or concentric) over an interval: infinitely many extrema,
    /// all at `distance`.
    Parallel { distance: f64 },
}

impl<P: Copy> CurveCurveExtrema<P> {
    /// The smallest extremal distance, if any.
    #[doc(alias = "LowerDistance")]
    pub fn min_distance(&self) -> Option<f64> {
        match self {
            CurveCurveExtrema::Isolated(e) => e.first().map(|e| e.distance),
            CurveCurveExtrema::Parallel { distance } => Some(*distance),
        }
    }

    /// The isolated extrema (empty for parallel curves).
    pub fn isolated(&self) -> &[CurveExtremum<P>] {
        match self {
            CurveCurveExtrema::Isolated(e) => e,
            CurveCurveExtrema::Parallel { .. } => &[],
        }
    }

    pub fn is_parallel(&self) -> bool {
        matches!(self, CurveCurveExtrema::Parallel { .. })
    }
}

/// Extrema between two 3D curves over their whole domains (use trimmed curves to restrict
/// them).
pub fn curve_curve(a: &Curve, b: &Curve) -> CurveCurveExtrema {
    curves_extrema(a, b, &kind_of(a), &kind_of(b))
}

/// Extrema between two 2D curves.
#[doc(alias = "Geom2dAPI_ExtremaCurveCurve")]
pub fn curve_curve2(a: &Curve2, b: &Curve2) -> CurveCurveExtrema<Point2> {
    curves_extrema(a, b, &kind_of2(a), &kind_of2(b))
}

/// The closest pair of points of two curves, including their end points (unlike the
/// extrema, which only include perpendicular pairs).
#[doc(alias = "GeomAPI_ExtremaCurveCurve::NearestPoints")]
pub fn closest_points(a: &Curve, b: &Curve) -> Option<CurveExtremum> {
    let interior = curve_curve(a, b).isolated().first().copied();
    // Or a pair with one end of one curve, nearest on the other.
    match (interior, closest_end_pair(a, b)) {
        (Some(i), Some(e)) => Some(if e.distance < i.distance { e } else { i }),
        (i, e) => i.or(e),
    }
}

/// What [`curves_extrema`] needs to know of a curve's basis: a line, a circle (with its
/// frame), or anything else.
enum Kind {
    Line { origin: [f64; 3], direction: [f64; 3] },
    Circle { center: [f64; 3], x: [f64; 3], y: [f64; 3], normal: [f64; 3], radius: f64 },
    Other,
}

fn kind_of(c: &Curve) -> Kind {
    let v = |d: crate::math::Direction| [d.x(), d.y(), d.z()];
    match c {
        Curve::Trimmed(t) => kind_of(t.basis()),
        Curve::Line(l) => Kind::Line { origin: l.origin.coords(), direction: v(l.direction) },
        Curve::Circle(c) => {
            let f = c.frame();
            Kind::Circle { center: f.origin().coords(), x: v(f.x()), y: v(f.y()), normal: v(f.z()), radius: c.radius() }
        }
        _ => Kind::Other,
    }
}

fn kind_of2(c: &Curve2) -> Kind {
    let v = |d: crate::math::Direction2| [d.x(), d.y(), 0.0];
    match c {
        Curve2::Trimmed(t) => kind_of2(t.basis()),
        Curve2::Line(l) => Kind::Line { origin: l.origin.coords(), direction: v(l.direction) },
        Curve2::Circle(c) => {
            let f = c.frame();
            Kind::Circle { center: f.origin().coords(), x: v(f.x()), y: v(f.y()), normal: [0.0, 0.0, 1.0], radius: c.radius() }
        }
        _ => Kind::Other,
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// An interval `[lo, hi]` (empty when `lo > hi`), as OCCT's `Bnd_Range`.
#[derive(Clone, Copy)]
struct Range {
    lo: f64,
    hi: f64,
}

impl Range {
    const EMPTY: Range = Range { lo: f64::INFINITY, hi: f64::NEG_INFINITY };
    fn of(a: f64, b: f64) -> Range {
        Range { lo: a.min(b), hi: a.max(b) }
    }
    fn add(&mut self, x: f64) {
        self.lo = self.lo.min(x);
        self.hi = self.hi.max(x);
    }
    fn common(self, other: Range) -> Range {
        Range { lo: self.lo.max(other.lo), hi: self.hi.min(other.hi) }
    }
    fn is_void(self) -> bool {
        self.lo > self.hi
    }
    fn delta(self) -> f64 {
        if self.is_void() { 0.0 } else { self.hi - self.lo }
    }
    fn middle(self) -> f64 {
        0.5 * (self.lo + self.hi)
    }
}

/// The extrema of two curves (2D or 3D): parallel and concentric lines and circles by
/// OCCT's rules (`Extrema_ExtCC`), other pairs by Newton steps from a grid of parameters.
fn curves_extrema<C>(a: &C, b: &C, ka: &Kind, kb: &Kind) -> CurveCurveExtrema<C::Point>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let infinite = |t: f64| t.abs() >= 0.5 * crate::tolerance::INFINITE;
    let (u11, u12, u21, u22) = (a.first_parameter(), a.last_parameter(), b.first_parameter(), b.last_parameter());
    let make = |t1: f64, t2: f64| {
        let (p1, p2) = (a.point(t1), b.point(t2));
        CurveExtremum { distance: norm(minus(p1.coords(), p2.coords())), t1, p1, t2, p2 }
    };
    // The closest pair of ends.
    let ends = || -> Option<CurveExtremum<C::Point>> {
        let mut best: Option<CurveExtremum<C::Point>> = None;
        for t1 in [u11, u12] {
            for t2 in [u21, u22] {
                if infinite(t1) || infinite(t2) {
                    continue;
                }
                let e = make(t1, t2);
                if best.is_none_or(|b| e.distance < b.distance) {
                    best = Some(e);
                }
            }
        }
        best
    };
    let confusion = crate::tolerance::LINEAR;
    let isolated = |mut v: Vec<CurveExtremum<C::Point>>| {
        v.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        CurveCurveExtrema::Isolated(v)
    };
    match (ka, kb) {
        (Kind::Line { origin: oa, direction: da }, Kind::Line { origin: ob, direction: db }) if norm(cross(*da, *db)) <= 1e-12 => {
            // Parallel lines: as many solutions where the ranges face each other.
            let distance = norm(cross(minus(*oa, *ob), *db));
            if (infinite(u11) && infinite(u12)) || (infinite(u21) && infinite(u22)) {
                return CurveCurveExtrema::Parallel { distance };
            }
            let opposite = dot(*da, *db) < 0.0;
            let on_b = |t: f64| dot(minus(a.point(t).coords(), *ob), *db);
            let mut projected = Range::EMPTY;
            let big = crate::tolerance::INFINITE;
            projected.add(if infinite(u11) { if opposite { big } else { -big } } else { on_b(u11) });
            projected.add(if infinite(u12) { if opposite { -big } else { big } } else { on_b(u12) });
            let common = Range::of(u21, u22).common(projected);
            if common.delta() > confusion {
                CurveCurveExtrema::Parallel { distance }
            } else if !common.is_void() {
                let t2 = common.middle();
                let t1 = dot(minus(b.point(t2).coords(), *oa), *da);
                isolated(vec![make(t1, t2)])
            } else {
                isolated(ends().into_iter().collect())
            }
        }
        (Kind::Line { origin, direction }, Kind::Circle { center, normal, radius, .. })
        | (Kind::Circle { center, normal, radius, .. }, Kind::Line { origin, direction })
            if norm(cross(*direction, *normal)) <= 1e-12 && norm(cross(minus(*center, *origin), *direction)) <= confusion =>
        {
            // A circle round the line: parallel if the circle's plane crosses the line's
            // range.
            let line_first = matches!(ka, Kind::Line { .. });
            let (lo, hi) = if line_first { (u11, u12) } else { (u21, u22) };
            let on_line = dot(minus(*center, *origin), *direction);
            match on_line >= lo - confusion && on_line <= hi + confusion {
                true => CurveCurveExtrema::Parallel { distance: *radius },
                false => isolated(Vec::new()),
            }
        }
        (
            Kind::Circle { center: ca, normal: na, radius: ra, .. },
            Kind::Circle { center: cb, x: xb, y: yb, normal: nb, radius: rb },
        ) if norm(cross(*na, *nb)) <= 1e-12 && norm(cross(minus(*cb, *ca), *na)) <= confusion => {
            let height = dot(minus(*cb, *ca), *na);
            let square_distance = (ra - rb).powi(2) + height * height;
            concentric(a, b, (u11, u12, u21, u22), (*cb, *xb, *yb), square_distance, *ra, *rb, &make, ends)
        }
        _ => {
            let mut found: Vec<CurveExtremum<C::Point>> = stationary_pairs(a, b, ka, kb).into_iter().map(|(t1, t2)| make(t1, t2)).collect();
            // Beyond lines and conics, the global minimum too where it is at an end (as
            // OCCT's global search finds it).
            if matches!(ka, Kind::Other) || matches!(kb, Kind::Other) {
                if let Some(e) = closest_end_pair(a, b) {
                    if found.iter().all(|f| f.distance > e.distance + confusion) {
                        found.push(e);
                    }
                }
            }
            isolated(found)
        }
    }
}

/// The closest pair of points of two curves with at least one at an end of its curve.
fn closest_end_pair<C>(a: &C, b: &C) -> Option<CurveExtremum<C::Point>>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let mut best: Option<CurveExtremum<C::Point>> = None;
    for (x, y, swapped) in [(a, b, false), (b, a, true)] {
        for t in [x.first_parameter(), x.last_parameter()] {
            if t.abs() >= 0.5 * crate::tolerance::INFINITE {
                continue;
            }
            let q = x.point(t);
            let Ok(found) = PointCurveSearch::new(y).find_with(q, crate::tolerance::LINEAR, SearchMode::Min, true) else { continue };
            if let Some(m) = found.min() {
                let (t1, p1, t2, p2) = if swapped { (m.parameter, m.point, t, q) } else { (t, q, m.parameter, m.point) };
                let e = CurveExtremum { distance: m.distance(), t1, p1, t2, p2 };
                if best.is_none_or(|b| e.distance < b.distance) {
                    best = Some(e);
                }
            }
        }
    }
    best
}

/// Coaxial circles (OCCT's circle-circle `PrepareParallelResult`): parallel where the
/// arcs face each other, else the touching or closest pairs.
#[allow(clippy::too_many_arguments)]
fn concentric<C>(
    a: &C,
    b: &C,
    (u11, u12, u21, u22): (f64, f64, f64, f64),
    (cb, xb, yb): ([f64; 3], [f64; 3], [f64; 3]),
    square_distance: f64,
    ra: f64,
    rb: f64,
    make: &dyn Fn(f64, f64) -> CurveExtremum<C::Point>,
    ends: impl Fn() -> Option<CurveExtremum<C::Point>>,
) -> CurveCurveExtrema<C::Point>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    use std::f64::consts::PI;
    let period = 2.0 * PI;
    let in_period = |x: f64, lo: f64| lo + (x - lo).rem_euclid(period);
    // The angle of a point about circle b.
    let angle_on_b = |p: [f64; 3]| {
        let r = minus(p, cb);
        dot(r, yb).atan2(dot(r, xb)).rem_euclid(period)
    };
    let precision = (ra.abs() * f64::EPSILON).max(rb.abs() * f64::EPSILON);
    let angular = 1e-12;
    let confusion = crate::tolerance::LINEAR;
    let (p11, tangent1) = a.d1(u11);
    let p12 = a.point(u12);
    let range = Range::of(u21, u22);
    let par1 = in_period(angle_on_b(p11.coords()), u21);
    let tangent2 = b.d1(par1).1;
    let opposite = dot(C::Point::vector_coords(tangent1), C::Point::vector_coords(tangent2)) < 0.0;
    let mut par2 = in_period(angle_on_b(p12.coords()), u21);
    if opposite {
        if range.delta() > angular && par1 - par2 < angular {
            par2 -= period;
        }
    } else if range.delta() > angular && par2 - par1 < angular {
        par2 += period;
    }
    let mut projected = Range::EMPTY;
    projected.add(par1 - period);
    projected.add(par2 - period);
    // The extrema from a point of circle b to arc a.
    let from_b = |t2: f64| -> Vec<PointExtremum<C::Point>> {
        PointCurveSearch::with_domain(a, u11, u12).find(b.point(t2), confusion).map(|r| r.extrema).unwrap_or_default()
    };
    let mut parallel: Option<f64> = None;
    let mut found: Vec<CurveExtremum<C::Point>> = Vec::new();
    let mut min_square = f64::INFINITY;
    for _ in 0..3 {
        let common = projected.common(range);
        if common.delta() > confusion {
            let extrema = from_b(common.middle());
            if let Some(min) = extrema.iter().map(|e| e.square_distance).reduce(f64::min) {
                if min <= min_square + (1.0 + min) * precision {
                    found.clear();
                    parallel = Some(min);
                    let delta = min - square_distance;
                    if delta * delta < 4.0 * min.max(square_distance) * confusion * confusion {
                        break;
                    }
                }
            }
        } else if !common.is_void() {
            let t2 = common.middle();
            let extrema = from_b(t2);
            let nearer = parallel.is_none() || extrema.iter().any(|e| e.square_distance < min_square);
            if nearer {
                found.clear();
                parallel = None;
                for e in &extrema {
                    found.push(make(e.parameter, t2));
                    min_square = min_square.min(e.square_distance);
                }
            }
        } else if let Some(e) = ends() {
            let d = e.distance * e.distance;
            if d <= min_square + (1.0 + d) * precision {
                parallel = None;
                found = vec![e];
                min_square = min_square.min(d);
            }
        }
        projected.lo += PI;
        projected.hi += PI;
    }
    match parallel {
        Some(d) => CurveCurveExtrema::Parallel { distance: d.sqrt() },
        None => {
            found.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            CurveCurveExtrema::Isolated(found)
        }
    }
}

/// The parameters of the stationary pairs of the distance between two curves (minima,
/// maxima and saddles), within their ranges: Newton steps on `(A - B)·A' = 0` and
/// `(A - B)·B' = 0` from the most promising nodes of a grid.
fn stationary_pairs<C>(a: &C, b: &C, ka: &Kind, kb: &Kind) -> Vec<(f64, f64)>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let infinite = |t: f64| t.abs() >= 0.5 * crate::tolerance::INFINITE;
    // Lines crossing at a single closest pair, in closed form.
    if let (Kind::Line { origin: oa, direction: da }, Kind::Line { origin: ob, direction: db }) = (ka, kb) {
        let w = minus(*oa, *ob);
        let (dd, e, f) = (dot(*da, *db), dot(*da, w), dot(*db, w));
        let det = 1.0 - dd * dd;
        let (s, t) = ((dd * f - e) / det, (f - dd * e) / det);
        let inside = |x: f64, lo: f64, hi: f64| x >= lo - 1e-12 * (1.0 + x.abs()) && x <= hi + 1e-12 * (1.0 + x.abs());
        return match inside(s, a.first_parameter(), a.last_parameter()) && inside(t, b.first_parameter(), b.last_parameter()) {
            true => vec![(s, t)],
            false => Vec::new(),
        };
    }
    // Finite ranges: an unbounded line is cut to where the other curve projects onto it.
    let range = |c: &C, k: &Kind, other: &C| -> (f64, f64) {
        let (lo, hi) = (c.first_parameter(), c.last_parameter());
        if !infinite(lo) && !infinite(hi) {
            return (lo, hi);
        }
        let (olo, ohi) = (other.first_parameter().max(-1e6), other.last_parameter().min(1e6));
        let mut seen = Range::EMPTY;
        if let Kind::Line { origin, direction } = k {
            for i in 0..=64 {
                let q = other.point(olo + (ohi - olo) * i as f64 / 64.0).coords();
                seen.add(dot(minus(q, *origin), *direction));
            }
            let margin = 1.0 + seen.delta();
            (if infinite(lo) { seen.lo - margin } else { lo }, if infinite(hi) { seen.hi + margin } else { hi })
        } else {
            (lo.max(-1e3), hi.min(1e3))
        }
    };
    let (ra, rb) = (range(a, ka, b), range(b, kb, a));
    let samples = |c: &C, (lo, hi): (f64, f64)| -> Vec<f64> {
        let n = match c.period() {
            Some(p) if c.is_periodic() => ((hi - lo) / p * 36.0).ceil().max(8.0) as usize,
            _ => 32,
        };
        let mut v: Vec<f64> = (0..=n).map(|i| if i == n { hi } else { lo + (hi - lo) * i as f64 / n as f64 }).collect();
        v.extend(c.breaks().into_iter().filter(|k| *k > lo && *k < hi));
        v.sort_by(f64::total_cmp);
        v.dedup();
        v
    };
    let (sa, sb) = (samples(a, ra), samples(b, rb));
    let gradient = |s: f64, t: f64| -> ([f64; 2], [[f64; 2]; 2], f64) {
        let (pa, a1, a2) = a.d2(s);
        let (pb, b1, b2) = b.d2(t);
        let (a1, a2, b1, b2) = (C::Point::vector_coords(a1), C::Point::vector_coords(a2), C::Point::vector_coords(b1), C::Point::vector_coords(b2));
        let r = minus(pa.coords(), pb.coords());
        let g = [dot(r, a1), -dot(r, b1)];
        let h = [[dot(a1, a1) + dot(r, a2), -dot(a1, b1)], [-dot(a1, b1), dot(b1, b1) - dot(r, b2)]];
        (g, h, (dot(a1, a1) * dot(b1, b1)).sqrt())
    };
    // Grid nodes where the gradient is smallest among their neighbours.
    let size: Vec<Vec<f64>> = sa
        .iter()
        .map(|&s| {
            sb.iter()
                .map(|&t| {
                    let (g, _, scale) = gradient(s, t);
                    (g[0] * g[0] + g[1] * g[1]).sqrt() / scale.max(1e-300)
                })
                .collect()
        })
        .collect();
    let mut seeds = Vec::new();
    for i in 0..sa.len() {
        for j in 0..sb.len() {
            let here = size[i][j];
            let lowest = (i.saturating_sub(1)..=(i + 1).min(sa.len() - 1))
                .all(|k| (j.saturating_sub(1)..=(j + 1).min(sb.len() - 1)).all(|l| size[k][l] >= here));
            if lowest {
                seeds.push((sa[i], sb[j]));
            }
        }
    }
    let wrap = |c: &C, x: f64, (lo, hi): (f64, f64)| match c.period().filter(|_| c.is_periodic()) {
        Some(p) => lo + (x - lo).rem_euclid(p),
        None => x,
    }
    .min(if c.is_periodic() { f64::INFINITY } else { hi })
    .max(if c.is_periodic() { f64::NEG_INFINITY } else { lo });
    let mut found: Vec<(f64, f64)> = Vec::new();
    for (mut s, mut t) in seeds {
        let mut converged = false;
        for _ in 0..60 {
            let (g, h, scale) = gradient(s, t);
            let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
            if det.abs() <= 1e-300 {
                break;
            }
            let ds = (-g[0] * h[1][1] + g[1] * h[0][1]) / det;
            let dt = (-g[1] * h[0][0] + g[0] * h[1][0]) / det;
            let (ns, nt) = (wrap(a, s + ds, ra), wrap(b, t + dt, rb));
            let small = (ns - s).abs() <= 1e-14 * (1.0 + s.abs()) && (nt - t).abs() <= 1e-14 * (1.0 + t.abs());
            (s, t) = (ns, nt);
            let (g, _, _) = gradient(s, t);
            if small || (g[0].abs() + g[1].abs()) <= 1e-13 * scale.max(1.0) * (1.0 + norm(minus(a.point(s).coords(), b.point(t).coords()))) {
                converged = true;
                break;
            }
        }
        if !converged {
            let (g, _, scale) = gradient(s, t);
            let size = (g[0].abs() + g[1].abs()) / scale.max(1e-300);
            if size > 1e-7 {
                continue;
            }
        }
        let inside = |c: &C, x: f64| c.is_periodic() || (x >= c.first_parameter() - 1e-12 && x <= c.last_parameter() + 1e-12);
        if !inside(a, s) || !inside(b, t) {
            continue;
        }
        let same = |x: &(f64, f64)| {
            let near = |p: f64, q: f64, c: &C| match c.period().filter(|_| c.is_periodic()) {
                Some(per) => {
                    let d = (p - q).rem_euclid(per);
                    d.min(per - d) <= 1e-7
                }
                None => (p - q).abs() <= 1e-7 * (1.0 + p.abs()),
            };
            near(x.0, s, a) && near(x.1, t, b)
        };
        if !found.iter().any(same) {
            found.push((s, t));
        }
    }
    found
}

/// The foot of a perpendicular from a point to a curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Projection<P = Point> {
    pub parameter: f64,
    pub point: P,
    pub distance: f64,
}

/// The orthogonal projections of `p` onto `curve`, sorted by increasing distance.
#[doc(alias = "GeomAPI_ProjectPointOnCurve")]
#[doc(alias = "Extrema_ExtPC")]
pub fn point_curve(p: Point, curve: &Curve) -> Vec<Projection> {
    projections(&PointCurveSearch::new(curve), p)
}

/// Every extremum the search finds (minima and maxima), as projections by increasing
/// distance.
fn projections<C>(search: &PointCurveSearch<C>, p: C::Point) -> Vec<Projection<C::Point>>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let Ok(found) = search.find(p, crate::tolerance::LINEAR) else { return Vec::new() };
    let mut out: Vec<Projection<C::Point>> =
        found.iter().map(|e| Projection { parameter: e.parameter, point: e.point, distance: e.distance() }).collect();
    out.sort_by(|a, b| a.distance.total_cmp(&b.distance));
    out
}

/// The orthogonal projections of `p` onto a 2D curve, sorted by increasing distance.
#[doc(alias = "Geom2dAPI_ProjectPointOnCurve")]
pub fn point_curve2(p: Point2, curve: &Curve2) -> Vec<Projection<Point2>> {
    projections(&PointCurveSearch::new(curve), p)
}

/// Which local extrema of the distance a [`PointCurveSearch`] looks for.
#[doc(alias = "ExtremaPC::SearchMode")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SearchMode {
    /// Minima and maxima.
    #[default]
    MinMax,
    Min,
    Max,
}

/// A local extremum of the distance from a point to a curve.
#[doc(alias = "ExtremaPC::ExtremumResult")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointExtremum<P = Point> {
    pub parameter: f64,
    pub point: P,
    pub square_distance: f64,
    pub is_minimum: bool,
    /// Both flags are set where the domain is a single parameter.
    pub is_maximum: bool,
}

impl<P> PointExtremum<P> {
    pub fn distance(&self) -> f64 {
        self.square_distance.sqrt()
    }
}

/// The result of a [`PointCurveSearch`].
#[doc(alias = "ExtremaPC::Result")]
#[derive(Debug, Clone, PartialEq)]
pub struct PointCurveExtrema<P = Point> {
    /// The isolated extrema, sorted by increasing parameter (empty if none were found or
    /// the solutions are not isolated).
    pub extrema: Vec<PointExtremum<P>>,
    /// `Some(d²)` if every point of the curve in the domain is at the same distance `d`
    /// (a point on the axis of a circle).
    pub infinite_square_distance: Option<f64>,
}

impl<P> PointCurveExtrema<P> {
    /// True if the solutions are not isolated (see `infinite_square_distance`).
    #[doc(alias = "IsInfinite")]
    pub fn is_infinite(&self) -> bool {
        self.infinite_square_distance.is_some()
    }

    #[doc(alias = "NbExt")]
    pub fn len(&self) -> usize {
        self.extrema.len()
    }

    pub fn is_empty(&self) -> bool {
        self.extrema.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, PointExtremum<P>> {
        self.extrema.iter()
    }

    /// The index of the closest extremum (the first one on ties).
    #[doc(alias = "MinIndex")]
    pub fn min_index(&self) -> Option<usize> {
        (0..self.extrema.len()).min_by(|&a, &b| {
            self.extrema[a].square_distance.total_cmp(&self.extrema[b].square_distance)
        })
    }

    /// The index of the farthest extremum (the first one on ties).
    #[doc(alias = "MaxIndex")]
    pub fn max_index(&self) -> Option<usize> {
        // `max_by` keeps the last of equal elements; compare reversed with `min_by`.
        (0..self.extrema.len()).min_by(|&a, &b| {
            self.extrema[b].square_distance.total_cmp(&self.extrema[a].square_distance)
        })
    }

    /// The closest extremum.
    pub fn min(&self) -> Option<&PointExtremum<P>> {
        self.min_index().map(|i| &self.extrema[i])
    }

    /// The farthest extremum.
    pub fn max(&self) -> Option<&PointExtremum<P>> {
        self.max_index().map(|i| &self.extrema[i])
    }

    /// The smallest squared distance of an extremum (infinity if there is none).
    #[doc(alias = "MinSquareDistance")]
    pub fn min_square_distance(&self) -> f64 {
        self.min().map_or(f64::INFINITY, |e| e.square_distance)
    }

    /// The largest squared distance of an extremum (0 if there is none).
    #[doc(alias = "MaxSquareDistance")]
    pub fn max_square_distance(&self) -> f64 {
        self.max().map_or(0.0, |e| e.square_distance)
    }
}

impl<P> Index<usize> for PointCurveExtrema<P> {
    type Output = PointExtremum<P>;
    fn index(&self, i: usize) -> &PointExtremum<P> {
        &self.extrema[i]
    }
}

impl<'a, P> IntoIterator for &'a PointCurveExtrema<P> {
    type Item = &'a PointExtremum<P>;
    type IntoIter = std::slice::Iter<'a, PointExtremum<P>>;
    fn into_iter(self) -> Self::IntoIter {
        self.extrema.iter()
    }
}

/// Why a [`PointCurveSearch`] failed.
#[doc(alias = "ExtremaPC::Status")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExtremaError {
    /// A non-positive or non-finite tolerance, a reversed or NaN domain, or a point with
    /// a coordinate that is NaN or at or beyond ±[`INFINITE`](crate::tolerance::INFINITE).
    InvalidInput,
    /// The computation did not converge.
    NumericalError,
}

impl fmt::Display for ExtremaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "invalid extrema input",
            Self::NumericalError => "extrema computation failed to converge",
        })
    }
}

impl std::error::Error for ExtremaError {}

/// Local extrema of the distance from points to a curve over a parameter domain.
///
/// Works on [`Curve`] and [`Curve2`] (any [`ParametricCurve`]); lines and conics are
/// solved analytically, other curves numerically. [`find`](Self::find) returns the
/// interior critical points; [`find_with_endpoints`](Self::find_with_endpoints) adds the
/// ends of a bounded domain where they are local extrema (an end of a closed curve's full
/// period is not an end). Ends at or beyond ±[`INFINITE`](crate::tolerance::INFINITE) are
/// unbounded. A root within [`PARAMETRIC`](crate::tolerance::PARAMETRIC) of a domain end
/// is reported exactly at that end.
#[doc(alias = "ExtremaPC_Curve")]
#[doc(alias = "ExtremaPC2d_Curve")]
#[doc(alias = "ExtremaPC_Line")]
#[doc(alias = "ExtremaPC_Circle")]
#[doc(alias = "ExtremaPC_Ellipse")]
#[doc(alias = "ExtremaPC_Hyperbola")]
#[doc(alias = "ExtremaPC_Parabola")]
#[doc(alias = "ExtremaPC_BezierCurve")]
#[doc(alias = "ExtremaPC_BSplineCurve")]
#[doc(alias = "ExtremaPC_OffsetCurve")]
#[derive(Debug, Clone, Copy)]
pub struct PointCurveSearch<'a, C> {
    curve: &'a C,
    first: f64,
    last: f64,
}

impl<'a, C: ParametricCurve> PointCurveSearch<'a, C>
where
    C::Point: ControlPoint<Vector = C::Vector>,
{
    /// Searches over the curve's whole domain.
    pub fn new(curve: &'a C) -> Self {
        Self { curve, first: curve.first_parameter(), last: curve.last_parameter() }
    }

    /// Searches over `[first, last]`. For a periodic curve the range may span several
    /// periods; for a bounded non-periodic one it is intersected with the curve's domain.
    #[doc(alias = "Domain1D")]
    pub fn with_domain(curve: &'a C, first: f64, last: f64) -> Self {
        Self { curve, first, last }
    }

    pub fn curve(&self) -> &C {
        self.curve
    }

    pub fn domain(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// The interior extrema (minima and maxima) of the distance from `p`; extrema closer
    /// than `tolerance` (in distance) are merged.
    #[doc(alias = "Perform")]
    pub fn find(&self, p: C::Point, tolerance: f64) -> Result<PointCurveExtrema<C::Point>, ExtremaError> {
        self.find_with(p, tolerance, SearchMode::MinMax, false)
    }

    /// [`find`](Self::find), plus the domain ends that are local extrema.
    #[doc(alias = "PerformWithEndpoints")]
    pub fn find_with_endpoints(
        &self,
        p: C::Point,
        tolerance: f64,
    ) -> Result<PointCurveExtrema<C::Point>, ExtremaError> {
        self.find_with(p, tolerance, SearchMode::MinMax, true)
    }

    /// The general form: only the extrema `mode` asks for, with or without the domain
    /// ends.
    pub fn find_with(
        &self,
        p: C::Point,
        tolerance: f64,
        mode: SearchMode,
        endpoints: bool,
    ) -> Result<PointCurveExtrema<C::Point>, ExtremaError> {
        use crate::tolerance::{INFINITE, PARAMETRIC};
        let q = p.coords();
        let valid_point = q.iter().all(|x| x.is_finite() && x.abs() < INFINITE);
        if !(tolerance > 0.0 && tolerance.is_finite()) || !valid_point || !(self.first <= self.last) {
            return Err(ExtremaError::InvalidInput);
        }
        let curve = self.curve;
        let period = curve.period().filter(|_| curve.is_periodic());
        // The domain, within a bounded non-periodic curve's own.
        let (mut first, mut last) = (self.first, self.last);
        if period.is_none() {
            first = first.max(curve.first_parameter());
            last = last.min(curve.last_parameter());
            if first > last {
                return Ok(PointCurveExtrema { extrema: Vec::new(), infinite_square_distance: None });
            }
        }
        if let Some(period) = period {
            // Steps of a tiny fraction of the period must be representable.
            let ulp = |x: f64| f64::EPSILON * x.abs().max(1.0);
            if ulp(first.abs().max(last.abs())) * 1e6 > period {
                return Err(ExtremaError::NumericalError);
            }
        }
        let geometry = Geometry { curve, q };
        let mut found: Vec<PointExtremum<C::Point>> = Vec::new();
        let flags = |is_minimum: bool, is_maximum: bool| match mode {
            SearchMode::MinMax => (is_minimum, is_maximum),
            SearchMode::Min => (is_minimum, false),
            SearchMode::Max => (false, is_maximum),
        };
        if first == last {
            if endpoints {
                let (is_minimum, is_maximum) = flags(true, true);
                found.push(geometry.extremum(first, is_minimum, is_maximum));
            }
            return Ok(PointCurveExtrema { extrema: found, infinite_square_distance: None });
        }
        // Ends at or beyond half of INFINITE are infinite (as OCCT's Precision::IsInfinite).
        let (lo_unbounded, hi_unbounded) = (first <= -0.5 * INFINITE, last >= 0.5 * INFINITE);
        // A domain of whole periods has no ends (its last point is its first).
        let whole_periods = period.is_some_and(|p| {
            let n = ((last - first) / p).round();
            n >= 1.0 && ((last - first) - n * p).abs() <= PARAMETRIC * n
        });
        let (a, b) = geometry.search_range(first, last, lo_unbounded, hi_unbounded);
        let samples = geometry.samples(a, b, period);
        // The distance is the same everywhere: no isolated extrema.
        let distances: Vec<f64> = samples.iter().map(|&t| geometry.square_distance(t).sqrt()).collect();
        let (dmin, dmax) = distances.iter().fold((f64::INFINITY, 0.0f64), |(lo, hi), &d| (lo.min(d), hi.max(d)));
        if dmax - dmin <= tolerance.min(1e-7 * dmax.max(1.0)) {
            let mean = distances.iter().map(|d| d * d).sum::<f64>() / distances.len() as f64;
            return Ok(PointCurveExtrema { extrema: Vec::new(), infinite_square_distance: Some(mean) });
        }
        // Roots of the slope: where it changes sign (from - to + at a minimum, also across
        // a corner), is zero, or touches zero.
        let mut roots: Vec<(f64, Option<bool>)> = Vec::new();
        let values: Vec<f64> = samples.iter().map(|&t| geometry.slope(t)).collect();
        for i in 0..samples.len() - 1 {
            let (t0, t1, f0, f1) = (samples[i], samples[i + 1], values[i], values[i + 1]);
            if f0 == 0.0 {
                roots.push((t0, None));
            } else if f0 * f1 < 0.0 {
                roots.push((geometry.root(t0, t1, f0, f1), Some(f0 < 0.0)));
            } else if let Some(t) = geometry.touching_root(t0, t1, f0, f1, tolerance) {
                roots.push((t, None));
            }
        }
        if values[values.len() - 1] == 0.0 {
            roots.push((samples[samples.len() - 1], None));
        }
        // A root just outside a bounded end (within the parametric tolerance) is at the end.
        for (end, outward, bounded) in [(first, -1.0, !lo_unbounded), (last, 1.0, !hi_unbounded)] {
            if !bounded || whole_periods {
                continue;
            }
            let (inside, outside) = (geometry.slope(end), geometry.slope(end + outward * PARAMETRIC));
            if inside != 0.0 && inside * outside < 0.0 {
                // Going up the parameter: from - to + is a minimum.
                let rising = if outward > 0.0 { inside < 0.0 } else { outside < 0.0 };
                roots.push((end, Some(rising)));
            }
        }
        roots.sort_by(|a, b| a.0.total_cmp(&b.0));
        // Roots at the ends of the domain, exactly there; outside it, dropped; the last end
        // of whole periods is the first.
        let mut kept: Vec<(f64, Option<bool>)> = Vec::new();
        for (t, kind) in roots {
            let t = if (t - first).abs() <= PARAMETRIC { first } else if (t - last).abs() <= PARAMETRIC { last } else { t };
            if t < first || t > last || (whole_periods && t == last) {
                continue;
            }
            // Merged with the previous one if at the same point (not a period away).
            if let Some(&(previous, _)) = kept.last() {
                let close = geometry.point_distance(previous, t) <= tolerance && period.is_none_or(|p| (t - previous).abs() < 0.5 * p);
                if close || t - previous <= PARAMETRIC {
                    continue;
                }
            }
            kept.push((t, kind));
        }
        // A closed curve's two ends are one point: one representative (the first).
        let seam = !whole_periods
            && !lo_unbounded
            && !hi_unbounded
            && first == curve.first_parameter()
            && last == curve.last_parameter()
            && geometry.point_distance(first, last) <= tolerance;
        if seam && kept.first().is_some_and(|k| k.0 == first) {
            kept.retain(|k| k.0 != last);
        }
        for (t, kind) in kept {
            let curvature = geometry.slope_derivative(t);
            let (is_minimum, is_maximum) = if let Some(minimum) = kind {
                (minimum, !minimum)
            } else if curvature != 0.0 {
                (curvature > 0.0, curvature < 0.0)
            } else {
                // Degenerate: by the distances on both sides.
                let h = 1e-4 * (b - a).min(1.0);
                let d = geometry.square_distance(t);
                let (l, r) = (geometry.square_distance(t - h), geometry.square_distance(t + h));
                (l >= d && r >= d, l <= d && r <= d)
            };
            let (is_minimum, is_maximum) = flags(is_minimum, is_maximum);
            if is_minimum || is_maximum {
                found.push(geometry.extremum(t, is_minimum, is_maximum));
            }
        }
        // The ends, where the distance is extremal there and they are not roots already.
        if endpoints && !whole_periods {
            let ends = [(first, !lo_unbounded, 1.0), (last, !hi_unbounded, -1.0)];
            for (t, bounded, inward) in ends {
                if !bounded || found.iter().any(|e| e.parameter == t) {
                    continue;
                }
                // The seam of a closed curve once only.
                if seam && t == last && found.iter().any(|e| e.parameter == first) {
                    continue;
                }
                // Moving inward, the distance grows (a minimum) or shrinks (a maximum).
                let f = inward * geometry.slope(t);
                // A zero slope makes the end a root: by the curvature of the distance.
                let f = if f == 0.0 { geometry.slope_derivative(t) } else { f };
                let (is_minimum, is_maximum) = flags(f > 0.0, f < 0.0);
                if is_minimum || is_maximum {
                    found.push(geometry.extremum(t, is_minimum, is_maximum));
                }
            }
        }
        found.sort_by(|x, y| x.parameter.total_cmp(&y.parameter));
        Ok(PointCurveExtrema { extrema: found, infinite_square_distance: None })
    }
}

/// A curve and a point, for [`PointCurveSearch`]: the squared distance `d²(t)/2` has slope
/// `f(t) = (C(t) - q)·C'(t)`, zero at the extrema, and `f'(t) = C'·C' + (C - q)·C''`.
struct Geometry<'a, C> {
    curve: &'a C,
    q: [f64; 3],
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn minus(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

impl<C> Geometry<'_, C>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    fn square_distance(&self, t: f64) -> f64 {
        let r = minus(self.curve.point(t).coords(), self.q);
        dot(r, r)
    }

    fn point_distance(&self, s: f64, t: f64) -> f64 {
        let r = minus(self.curve.point(s).coords(), self.curve.point(t).coords());
        dot(r, r).sqrt()
    }

    fn slope(&self, t: f64) -> f64 {
        let (c, d) = self.curve.d1(t);
        dot(minus(c.coords(), self.q), C::Point::vector_coords(d))
    }

    fn slope_derivative(&self, t: f64) -> f64 {
        let (c, d1, d2) = self.curve.d2(t);
        let d1 = C::Point::vector_coords(d1);
        dot(d1, d1) + dot(minus(c.coords(), self.q), C::Point::vector_coords(d2))
    }

    fn extremum(&self, t: f64, is_minimum: bool, is_maximum: bool) -> PointExtremum<C::Point> {
        let point = self.curve.point(t);
        let r = minus(point.coords(), self.q);
        PointExtremum { parameter: t, point, square_distance: dot(r, r), is_minimum, is_maximum }
    }

    /// The finite range searched: the domain, an unbounded side cut where the curve is far
    /// beyond the point (no extrema lie past it on a line, parabola or hyperbola).
    fn search_range(&self, first: f64, last: f64, lo_unbounded: bool, hi_unbounded: bool) -> (f64, f64) {
        if !lo_unbounded && !hi_unbounded {
            return (first, last);
        }
        let center = if lo_unbounded && hi_unbounded { 0.0 } else if lo_unbounded { last } else { first };
        let base = self.square_distance(center).sqrt();
        let origin = self.curve.point(center).coords();
        let far = |t: f64| {
            let r = minus(self.curve.point(t).coords(), origin);
            dot(r, r).sqrt() > 4.0 * (base + 1.0)
        };
        let (mut a, mut b) = (if lo_unbounded { center - 1.0 } else { first }, if hi_unbounded { center + 1.0 } else { last });
        for _ in 0..200 {
            let (a_far, b_far) = (!lo_unbounded || (far(a) && self.slope(a) < 0.0), !hi_unbounded || (far(b) && self.slope(b) > 0.0));
            if a_far && b_far {
                break;
            }
            // Not past where the curve can be evaluated (cosh overflows).
            let finite = |t: f64| self.curve.point(t).coords().iter().all(|x| x.is_finite());
            let (next_a, next_b) = (center - 2.0 * (center - a), center + 2.0 * (b - center));
            let (grow_a, grow_b) = (!a_far && finite(next_a), !b_far && finite(next_b));
            if !grow_a && !grow_b {
                break;
            }
            if grow_a {
                a = next_a;
            }
            if grow_b {
                b = next_b;
            }
        }
        (a, b)
    }

    /// Parameters splitting `[a, b]` finely enough that each piece holds at most one
    /// extremum: evenly, then halved where the tangent turns much.
    fn samples(&self, a: f64, b: f64, period: Option<f64>) -> Vec<f64> {
        let base = match period {
            Some(p) => ((b - a) / p * 48.0).ceil().max(8.0) as usize,
            None => 64,
        };
        let tangent = |t: f64| {
            let d = C::Point::vector_coords(self.curve.d1(t).1);
            let n = dot(d, d).sqrt();
            if n > 0.0 { [d[0] / n, d[1] / n, d[2] / n] } else { [0.0; 3] }
        };
        // Even steps, and the curve's breaks (each period's) where it may have corners.
        let mut grid: Vec<f64> = (0..=base).map(|i| if i == base { b } else { a + (b - a) * i as f64 / base as f64 }).collect();
        for k in self.curve.breaks() {
            match period {
                Some(p) => {
                    let mut t = k + ((a - k) / p).ceil() * p;
                    while t < b {
                        grid.push(t);
                        t += p;
                    }
                }
                None if a < k && k < b => grid.push(k),
                None => {}
            }
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup();
        let mut out = vec![a];
        let mut stack: Vec<(f64, f64, usize)> = grid.windows(2).rev().filter(|w| w[0] < w[1]).map(|w| (w[0], w[1], 0)).collect();
        while let Some((s, t, depth)) = stack.pop() {
            // Turning by more than about 25° (or a vanishing tangent): halved.
            if depth < 8 && dot(tangent(s), tangent(t)) < 0.9 {
                let m = 0.5 * (s + t);
                stack.push((m, t, depth + 1));
                stack.push((s, m, depth + 1));
                continue;
            }
            out.push(t);
        }
        out
    }

    /// The root of the slope in `[s, t]`, where it changes sign: Newton steps kept inside
    /// the bracket, halving it where they leave.
    fn root(&self, mut s: f64, mut t: f64, mut fs: f64, _ft: f64) -> f64 {
        let mut x = 0.5 * (s + t);
        for _ in 0..100 {
            let f = self.slope(x);
            if f == 0.0 {
                return x;
            }
            if (f < 0.0) == (fs < 0.0) {
                s = x;
                fs = f;
            } else {
                t = x;
            }
            let df = self.slope_derivative(x);
            let newton = x - f / df;
            let next = if df != 0.0 && newton > s && newton < t { newton } else { 0.5 * (s + t) };
            if (next - x).abs() <= 1e-16 * (1.0 + x.abs()) || t - s <= 1e-16 * (1.0 + x.abs()) {
                return next;
            }
            x = next;
        }
        x
    }

    /// A root where the slope touches zero without changing sign (a degenerate extremum):
    /// where its size has a minimum near zero inside `[s, t]`.
    fn touching_root(&self, s: f64, t: f64, fs: f64, ft: f64, tolerance: f64) -> Option<f64> {
        // Only worth a look where the slope's size dips in the middle.
        let m = 0.5 * (s + t);
        let fm = self.slope(m);
        if fm.abs() >= fs.abs().min(ft.abs()) || (fm < 0.0) != (fs < 0.0) {
            return None;
        }
        // Golden-section search for the smallest |f|.
        let (mut a, mut b) = (s, t);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..80 {
            let (x1, x2) = (b - g * (b - a), a + g * (b - a));
            if self.slope(x1).abs() < self.slope(x2).abs() {
                b = x2;
            } else {
                a = x1;
            }
        }
        let x = 0.5 * (a + b);
        let speed = dot(C::Point::vector_coords(self.curve.d1(x).1), C::Point::vector_coords(self.curve.d1(x).1)).sqrt();
        (self.slope(x).abs() <= tolerance * speed.max(1e-300)).then_some(x)
    }
}

/// One extremum between a curve and a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceExtremum {
    pub distance: f64,
    pub t: f64,
    pub curve_point: Point,
    pub u: f64,
    pub v: f64,
    pub surface_point: Point,
}

/// The extrema between a curve and a surface.
#[doc(alias = "GeomAPI_ExtremaCurveSurface")]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveSurfaceExtrema {
    /// Isolated extrema, sorted by increasing distance.
    Isolated(Vec<CurveSurfaceExtremum>),
    /// The curve runs at constant distance from the surface (e.g. a line parallel to a
    /// plane): infinitely many extrema at `distance`.
    Parallel { distance: f64 },
}

/// The local minima of the distance from a curve to a surface (each curve point's
/// distance to its nearest surface point), or [`Parallel`](CurveSurfaceExtrema::Parallel)
/// where that distance is the same all along. An unbounded curve is searched over
/// `[-1000, 1000]` of its parameter.
pub fn curve_surface(curve: &Curve, surface: &Surface) -> CurveSurfaceExtrema {
    let (lo, hi) = (curve.first_parameter().max(-1e3), curve.last_parameter().min(1e3));
    let n = 96;
    let mut samples: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(n + 1);
    let mut start = None;
    for i in 0..=n {
        let t = if i == n { hi } else { lo + (hi - lo) * i as f64 / n as f64 };
        let (u, v, d) = nearest_on_surface(curve.point(t), surface, start);
        start = Some((u, v));
        samples.push((t, u, v, d));
    }
    let (dmin, dmax) = samples.iter().fold((f64::INFINITY, 0.0f64), |(a, b), s| (a.min(s.3), b.max(s.3)));
    if dmax - dmin <= crate::tolerance::LINEAR {
        return CurveSurfaceExtrema::Parallel { distance: dmin };
    }
    let at = |t: f64, start: (f64, f64)| nearest_on_surface(curve.point(t), surface, Some(start));
    let mut found: Vec<CurveSurfaceExtremum> = Vec::new();
    for i in 0..samples.len() {
        let d = samples[i].3;
        let left = if i > 0 { samples[i - 1].3 } else { f64::INFINITY };
        let right = if i + 1 < samples.len() { samples[i + 1].3 } else { f64::INFINITY };
        if !(d <= left && d <= right) {
            continue;
        }
        // Golden-section search for the smallest distance around the sample.
        let (mut a, mut b) = (samples[i.saturating_sub(1)].0, samples[(i + 1).min(samples.len() - 1)].0);
        let mut guess = (samples[i].1, samples[i].2);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..90 {
            let (x1, x2) = (b - g * (b - a), a + g * (b - a));
            let (f1, f2) = (at(x1, guess), at(x2, guess));
            if f1.2 < f2.2 {
                b = x2;
                guess = (f1.0, f1.1);
            } else {
                a = x1;
                guess = (f2.0, f2.1);
            }
        }
        let t = 0.5 * (a + b);
        let (u, v, distance) = at(t, guess);
        if found.iter().any(|e| (e.t - t).abs() <= 1e-9 * (1.0 + t.abs())) {
            continue;
        }
        use super::parametric::ParametricSurface;
        found.push(CurveSurfaceExtremum { distance, t, curve_point: curve.point(t), u, v, surface_point: surface.point(u, v) });
    }
    found.sort_by(|x, y| x.distance.total_cmp(&y.distance));
    CurveSurfaceExtrema::Isolated(found)
}

/// Newton steps from `start` towards the point of `surface` nearest to `p` (a local
/// minimum of the distance, within the surface's bounds): its parameters and the distance.
pub(crate) fn newton_on_surface(p: Point, surface: &Surface, start: (f64, f64)) -> (f64, f64, f64) {
    let foot = newton_foot(p, surface, start);
    (foot.u, foot.v, foot.point.distance(p))
}

/// Where Newton steps from `start` end on the surface: the parameters, the point and its
/// first derivatives there.
pub(crate) struct SurfaceFoot {
    pub u: f64,
    pub v: f64,
    pub point: Point,
    pub du: crate::math::Vector,
    pub dv: crate::math::Vector,
}

/// A Newton step `e` at most a millionth of the range (`small` a trillionth) and a hundredth
/// of the step `last` before it.
fn close_to(e: f64, small: f64, last: f64) -> bool {
    e <= 1e6 * small && e <= 1e-2 * last
}

/// [`newton_on_surface`] with the point's derivatives (to spare evaluating them again).
pub(crate) fn newton_foot(p: Point, surface: &Surface, start: (f64, f64)) -> SurfaceFoot {
    newton_foot_to(p, surface, start, true, false)
}

/// [`newton_foot`] with every step taken to the rounding noise (none spared as converged):
/// for distances whose signs decide where surfaces cross, as booleans' do.
pub(crate) fn newton_foot_to_the_last(p: Point, surface: &Surface, start: (f64, f64)) -> SurfaceFoot {
    newton_foot_to(p, surface, start, false, false)
}

/// [`newton_foot_to_the_last`] from a sample of a grid over the surface (a fresh search, no
/// earlier foot to stay near): where the steps circle the foot without closing in, they stop
/// at the nearest point met.
pub(crate) fn newton_foot_from_sample(p: Point, surface: &Surface, start: (f64, f64)) -> SurfaceFoot {
    newton_foot_to(p, surface, start, false, true)
}

fn newton_foot_to(p: Point, surface: &Surface, start: (f64, f64), spare_last: bool, stall: bool) -> SurfaceFoot {
    use super::parametric::ParametricSurface;
    // An offset surface's normal lines are its basis's: a point's foot on it is where its
    // foot on the basis is, found with the basis's second derivatives (the offset's need the
    // basis's third).
    if let Surface::Offset(o) = surface {
        let foot = newton_foot_to(p, o.basis(), start, spare_last, stall);
        let (point, du, dv) = surface.d1(foot.u, foot.v);
        return SurfaceFoot { u: foot.u, v: foot.v, point, du, dv };
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let (u_periodic, v_periodic) = (surface.is_u_periodic(), surface.is_v_periodic());
    let clamp = |u: f64, v: f64| (if u_periodic { u } else { u.clamp(u0, u1) }, if v_periodic { v } else { v.clamp(v0, v1) });
    let (mut u, mut v) = clamp(start.0, start.1);
    // Converged once a step moves by a trillionth of the range (the next would be at the
    // rounding noise, where steps no longer shrink).
    let span = |a: f64, b: f64, t: f64| if (b - a).is_finite() { b - a } else { 1.0 + t.abs() };
    let (u_small, v_small) = (1e-12 * span(u0, u1, u), 1e-12 * span(v0, v1, v));
    let mut before = (f64::INFINITY, f64::INFINITY);
    // (With `stall`: where the distance is not convex (a point far beyond the surface's
    // centres of curvature), steepest-descent steps can circle the foot without closing in,
    // to the step limit; the nearest point met is kept, and once STALLED such steps in a row
    // bring the point no nearer, the steps stop there.)
    const STALLED: usize = 4;
    let (mut nearest, mut stalled) = ((f64::INFINITY, u, v), 0);
    for step in 0..50 {
        let d = surface.d2(u, v);
        let r = d.point - p;
        let f = r.dot(r);
        if f < nearest.0 {
            nearest = (f, u, v);
            stalled = 0;
        }
        let (gu, gv) = (r.dot(d.du), r.dot(d.dv));
        let (huu, huv, hvv) = (d.du.dot(d.du) + r.dot(d.duu), d.du.dot(d.dv) + r.dot(d.duv), d.dv.dot(d.dv) + r.dot(d.dvv));
        let det = huu * hvv - huv * huv;
        let convex = det.abs() > 1e-300 && huu > 0.0 && det > 0.0;
        if stall && !convex && f > nearest.0 {
            stalled += 1;
            if stalled >= STALLED {
                (u, v) = (nearest.1, nearest.2);
                break;
            }
        }
        let (du, dv) = if convex {
            ((-gu * hvv + gv * huv) / det, (-gv * huu + gu * huv) / det)
        } else {
            // Not convex here: a gradient step instead.
            let len2 = d.du.dot(d.du).max(d.dv.dot(d.dv)).max(1e-300);
            (-gu / len2, -gv / len2)
        };
        // Steps not settling (as on a bound, the step in both cut short there): a parameter on
        // its bound that the step would take past it stays there, and the other alone moves,
        // to the nearest point along that edge.
        const SETTLING: usize = 12;
        let held = |x: f64, dx: f64, lo: f64, hi: f64, periodic: bool| step >= SETTLING && !periodic && ((x <= lo && dx < 0.0) || (x >= hi && dx > 0.0));
        let (du, dv) = match (held(u, du, u0, u1, u_periodic), held(v, dv, v0, v1, v_periodic)) {
            (true, true) => (0.0, 0.0),
            (true, false) => (0.0, if hvv > 0.0 { -gv / hvv } else { -gv / d.dv.dot(d.dv).max(1e-300) }),
            (false, true) => (if huu > 0.0 { -gu / huu } else { -gu / d.du.dot(d.du).max(1e-300) }, 0.0),
            (false, false) => (du, dv),
        };
        let (nu, nv) = clamp(u + du, v + dv);
        let done = (nu - u).abs() <= u_small && (nv - v).abs() <= v_small;
        (u, v) = (nu, nv);
        // Close (where sparing it), the steps shrinking as Newton's do near a foot (each within about the square
        // of the one before: this one a millionth of the range, a hundredth of the last at
        // most): the next would be at the rounding noise, and is not taken. (The point and
        // its derivatives are evaluated where it ends, not carried on from those before: a
        // point far off the surface magnifies their errors in its distance along the normal.)
        let (eu, ev) = (du.abs(), dv.abs());
        if done || spare_last && convex && close_to(eu, u_small, before.0) && close_to(ev, v_small, before.1) {
            break;
        }
        before = (eu.max(u_small), ev.max(v_small));
    }
    let (point, du, dv) = surface.d1(u, v);
    SurfaceFoot { u, v, point, du, dv }
}

/// What stands for a surface in caches of what is found on it: a B-spline or Bézier
/// surface's identity, else three of its points.
#[derive(PartialEq)]
pub(crate) enum SurfaceKey {
    Shared(crate::geometry::SurfaceId),
    Points([u64; 9]),
}

impl SurfaceKey {
    pub(crate) fn of(surface: &Surface) -> SurfaceKey {
        match surface.shared_id() {
            Some(id) => SurfaceKey::Shared(id),
            None => SurfaceKey::Points(fingerprint(surface)),
        }
    }
}

/// Three points of a surface, standing for it in the cache.
fn fingerprint(surface: &Surface) -> [u64; 9] {
    use super::parametric::ParametricSurface;
    let (u0, u1, v0, v1) = surface.bounds();
    let mut out = [0; 9];
    for (k, (u, v)) in [(0.31, 0.27), (0.5, 0.5), (0.77, 0.69)].into_iter().enumerate() {
        let q = surface.point(u0 + (u1 - u0) * u, v0 + (v1 - v0) * v);
        out[3 * k..3 * k + 3].copy_from_slice(&[q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]);
    }
    out
}

/// For a B-spline surface winding round six times or more along u or v (a helical sweep, its
/// turns close together, more than its usual samples tell apart), sixteen samples a turn
/// that way, so that a point's nearest sample is on the turn it is on; else 0. The turning is measured along the middle iso-lines.
/// Kept for the surface.
pub(crate) fn windings(surface: &Surface) -> (usize, usize) {
    windings_from(surface, 6.0)
}

/// [`windings`] for a surface winding round `turns` times or more; under six, only where its
/// turns come close together (not one only bending to and fro, as a fillet along a wavy
/// edge).
pub(crate) fn windings_from(surface: &Surface, turns: f64) -> (usize, usize) {
    let ((u_turn, u_close), (v_turn, v_close)) = turning(surface);
    let winding = |turn: f64, close: bool| if turn >= turns * std::f64::consts::TAU && (close || turn >= 6.0 * std::f64::consts::TAU) { ((16.0 * turn / std::f64::consts::TAU).ceil() as usize).min(1024) } else { 0 };
    (winding(u_turn, u_close), winding(v_turn, v_close))
}

/// For a B-spline surface coming back near itself along its middle iso-line in u or in v
/// (three quarters of it near another part, as a thread's flank near its next turn: a sixth
/// of the way along or more, within a twelfth of its length), samples enough that way for two in each gap between its turns (at most 512);
/// else 0. Kept for the surface.
pub(crate) fn close_samples(surface: &Surface) -> (usize, usize) {
    static FOUND: std::sync::Mutex<Vec<(SurfaceKey, (usize, usize))>> = std::sync::Mutex::new(Vec::new());
    use super::parametric::ParametricSurface;
    let Surface::BSpline(_) = surface else { return (0, 0) };
    let key = SurfaceKey::of(surface);
    let lock = || FOUND.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Some(&(_, found)) = lock().iter().find(|(k, _)| *k == key) {
        return found;
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let samples = |points: Vec<Point>| -> usize {
        let n = points.len() - 1;
        let length: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
        // Each point's nearest of those apart from it along the line by a sixth or more,
        // either way round (a closed line's ends are one point); the largest of those of three
        // quarters of the points: most of the line near another part of it, as a thread's
        // turns are, not only some of it (a line pinched to a pole, or passing near itself
        // once, as a fillet along curves crossing nearly).
        let mut nearest: Vec<f64> = (0..=n)
            .map(|i| {
                (0..=n)
                    .filter(|&j| {
                        let apart = i.abs_diff(j);
                        apart >= n / 6 && n - apart >= n / 6
                    })
                    .map(|j| points[i].distance(points[j]))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        nearest.sort_by(f64::total_cmp);
        let gap = nearest[3 * n / 4];
        if !(length > 0.0) || !(gap < length / 12.0) {
            return 0;
        }
        ((2.0 * length / gap.max(1e-300)).ceil() as usize).min(512)
    };
    let n = 128;
    let found = (
        samples((0..=n).map(|i| surface.point(u0 + (u1 - u0) * i as f64 / n as f64, 0.5 * (v0 + v1))).collect()),
        samples((0..=n).map(|j| surface.point(0.5 * (u0 + u1), v0 + (v1 - v0) * j as f64 / n as f64)).collect()),
    );
    let mut cache = lock();
    if cache.len() >= 256 {
        cache.remove(0);
    }
    cache.push((key, found));
    found
}

/// How far a B-spline surface turns (in radians) along its middle iso-lines in u and in v,
/// and whether each comes back near itself all along (most of its points within a few of
/// its chords of a point a sixth of it along or more away: the next turn, not only its ends
/// overlapping, as a fillet's carried past the ends of a closed edge); none for others.
/// Kept for the surface.
fn turning(surface: &Surface) -> ((f64, bool), (f64, bool)) {
    static FOUND: std::sync::Mutex<Vec<(SurfaceKey, ((f64, bool), (f64, bool)))>> = std::sync::Mutex::new(Vec::new());
    use super::parametric::ParametricSurface;
    use crate::math::Vector;
    let Surface::BSpline(_) = surface else { return ((0.0, false), (0.0, false)) };
    let key = SurfaceKey::of(surface);
    let lock = || FOUND.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Some(&(_, found)) = lock().iter().find(|(k, _)| *k == key) {
        return found;
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let angle = |x: Vector, y: Vector| {
        let (lx, ly) = (x.length(), y.length());
        if lx > 0.0 && ly > 0.0 { (x.dot(y) / (lx * ly)).clamp(-1.0, 1.0).acos() } else { 0.0 }
    };
    let turned = |points: Vec<Point>| {
        let turn = points.windows(3).map(|w| angle(w[1] - w[0], w[2] - w[1])).sum::<f64>();
        let n = points.len();
        let chord = points.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>() / (n - 1) as f64;
        // (Counted only until the half is reached, or can no longer be.)
        let mut returns = 0;
        for i in 0..n {
            if 2 * returns >= n || 2 * (returns + n - i) < n {
                break;
            }
            if (0..n).any(|j| i.abs_diff(j) >= n / 6 && points[i].distance(points[j]) <= 4.0 * chord) {
                returns += 1;
            }
        }
        let close = 2 * returns >= n;
        (turn, close)
    };
    // (Over chords long enough that a fitted surface's ripples between its knots do not
    // add up: up to some forty turns told.)
    let (nu, nv) = (128, 128);
    let u_turn = turned((0..=nu).map(|i| surface.point(u0 + (u1 - u0) * i as f64 / nu as f64, 0.5 * (v0 + v1))).collect());
    let v_turn = turned((0..=nv).map(|j| surface.point(0.5 * (u0 + u1), v0 + (v1 - v0) * j as f64 / nv as f64)).collect());
    let found = (u_turn, v_turn);
    let mut cache = lock();
    if cache.len() >= 256 {
        cache.remove(0);
    }
    cache.push((key, found));
    found
}

/// The parameters `(u, v)` of the point of `surface` closest to `p` within its bounds, and the
/// distance: in closed form on elementary surfaces, otherwise by Newton steps from `start`
/// (a nearby solution, as when following a curve) or from the best of a grid of samples.
#[doc(alias = "GeomAPI_ProjectPointOnSurf")]
pub fn nearest_on_surface(p: Point, surface: &Surface, start: Option<(f64, f64)>) -> (f64, f64, f64) {
    use super::parametric::ParametricSurface;
    if let Some((u, v)) = surface.parameters(p) {
        return (u, v, surface.point(u, v).distance(p));
    }
    let (u0, u1, v0, v1) = surface.bounds();
    // Infinite directions (an extrusion's) are searched near the origin.
    let finite = |a: f64, b: f64| if a.abs() < 1e9 && b.abs() < 1e9 { (a, b) } else { (-1e3, 1e3) };
    let ((u0, u1), (v0, v1)) = (finite(u0, u1), finite(v0, v1));
    let newton = |start: (f64, f64)| newton_on_surface(p, surface, start);
    if let Some(start) = start {
        let found = newton(start);
        // Accepted where the point lies on the surface (following a curve on it).
        if found.2 <= 1e-6 * (1.0 + (p - Point::ORIGIN).length()) {
            return found;
        }
    }
    // A grid (finer on B-splines of many spans), then Newton steps from its few nearest
    // points (the nearest alone may be a pole, where the steps stall, or a local minimum).
    // (And on one winding round many times, enough to tell its turns apart.)
    let (nu, nv) = match surface {
        Surface::BSpline(b) => {
            let (wu, wv) = windings(surface);
            ((2 * b.u_knots().len()).clamp(16, 64).max(wu), (2 * b.v_knots().len()).clamp(16, 64).max(wv))
        }
        _ => (16, 16),
    };
    let points = grid_points(surface, (u0, u1, v0, v1), nu, nv);
    let mut grid: Vec<(f64, f64, f64)> = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            let (u, v) = (u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
            grid.push((u, v, points[i * (nv + 1) + j].distance(p)));
        }
    }
    // By distance, ties in grid order (as a stable sort): only the nearest few are sorted,
    // the rest only if they run out before the starts are found.
    let order = |a: &(f64, f64, f64), b: &(f64, f64, f64)| a.2.total_cmp(&b.2).then(a.0.total_cmp(&b.0)).then(a.1.total_cmp(&b.1));
    let nearest = 32.min(grid.len());
    if nearest < grid.len() {
        grid.select_nth_unstable_by(nearest, order);
    }
    grid[..nearest].sort_unstable_by(order);
    let mut best = grid[0];
    // Starts at distinct points (a collapsed edge, as a pole, is one point for many
    // parameters), each moved a little into the surface where it collapses (the steps stall
    // on it).
    let size = surface.point(u0, v0).distance(surface.point(u1, v1)) + surface.point(u0, v1).distance(surface.point(u1, v0));
    let same = 1e-9 * (1.0 + size);
    let mut starts: Vec<Point> = Vec::new();
    for k in 0..grid.len() {
        if starts.len() >= 6 {
            break;
        }
        if k == nearest {
            grid[nearest..].sort_unstable_by(order);
        }
        let (u, v, _) = grid[k];
        let q = surface.point(u, v);
        if starts.iter().any(|s| s.distance(q) <= same) {
            continue;
        }
        starts.push(q);
        let (du, dv) = (1e-3 * (u1 - u0), 1e-3 * (v1 - v0));
        let inward = |x: f64, lo: f64, hi: f64, d: f64| if x - lo < hi - x { x + d } else { x - d };
        let collapsed_u = q.distance(surface.point(inward(u, u0, u1, du), v)) <= same;
        let collapsed_v = q.distance(surface.point(u, inward(v, v0, v1, dv))) <= same;
        let start = (if collapsed_v { inward(u, u0, u1, du) } else { u }, if collapsed_u { inward(v, v0, v1, dv) } else { v });
        let found = newton(start);
        if found.2 < best.2 {
            best = found;
        }
    }
    best
}

/// What stands for a surface in the cache of grids: a B-spline or Bézier surface's
/// identity, else three of its points (by their bits); with the grid's steps.
#[derive(PartialEq)]
enum GridKey {
    Shared(super::SurfaceId, usize, usize),
    Points([u64; 9], usize, usize),
}

/// Grids of points of the surfaces searched lately (shared by the threads of an
/// operation): a search from scratch evaluates many.
static GRIDS: std::sync::Mutex<Vec<(GridKey, std::sync::Arc<Vec<Point>>)>> = std::sync::Mutex::new(Vec::new());

/// The points of `surface` at the `(nu + 1) × (nv + 1)` grid over the bounds given, row by
/// row along u (cached).
fn grid_points(surface: &Surface, (u0, u1, v0, v1): (f64, f64, f64, f64), nu: usize, nv: usize) -> std::sync::Arc<Vec<Point>> {
    use super::parametric::ParametricSurface;
    let key = match surface.shared_id() {
        Some(id) => GridKey::Shared(id, nu, nv),
        None => {
            let mut bits = [0; 9];
            for (k, (a, b)) in [(0.31, 0.27), (0.5, 0.5), (0.77, 0.69)].into_iter().enumerate() {
                let q = surface.point(u0 + (u1 - u0) * a, v0 + (v1 - v0) * b);
                bits[3 * k..3 * k + 3].copy_from_slice(&[q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]);
            }
            GridKey::Points(bits, nu, nv)
        }
    };
    let lock = || GRIDS.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Some((_, grid)) = lock().iter().find(|(k, _)| *k == key) {
        return grid.clone();
    }
    // (Made outside the lock; another thread may make it too.)
    let mut points = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            points.push(surface.point(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64));
        }
    }
    let grid = std::sync::Arc::new(points);
    let mut cache = lock();
    if !cache.iter().any(|(k, _)| *k == key) {
        if cache.len() >= 64 {
            cache.remove(0);
        }
        cache.push((key, grid.clone()));
    }
    grid
}

/// The point of `curve` over `[first, last]` closest to `p` (its ends included).
#[doc(alias = "GeomAPI_ProjectPointOnCurve::LowerDistance")]
pub fn nearest_on_curve(p: Point, curve: &Curve, first: f64, last: f64) -> Projection {
    NearestOnCurve::new(curve, first, last).nearest(p)
}

/// [`nearest_on_curve`] for many points onto one curve over one range: its coarse samples
/// taken once.
pub struct NearestOnCurve<'a> {
    curve: &'a Curve,
    lo: f64,
    hi: f64,
    period: Option<f64>,
    samples: Vec<(f64, Point)>,
}

impl<'a> NearestOnCurve<'a> {
    pub fn new(curve: &'a Curve, first: f64, last: f64) -> NearestOnCurve<'a> {
        use crate::geometry::ParametricCurve;
        let (lo, hi) = (first.min(last), first.max(last));
        // Over a whole period of a periodic curve the parameter wraps round (so that a point
        // just before the curve's origin is not stuck at the range's start).
        let period = curve.period().filter(|p| hi - lo >= p - 1e-12);
        // Coarse samples, then (for each point) Newton steps on (C(t) - p)·C'(t) = 0 from the
        // best.
        let n = if matches!(curve, Curve::Line(_)) { 1 } else { 64 };
        let samples = (0..=n)
            .map(|i| {
                let t = lo + (hi - lo) * i as f64 / n as f64;
                (t, curve.point(t))
            })
            .collect();
        NearestOnCurve { curve, lo, hi, period, samples }
    }

    /// The point of the curve closest to `p`.
    pub fn nearest(&self, p: Point) -> Projection {
        use crate::geometry::ParametricCurve;
        let (curve, lo, hi, period) = (self.curve, self.lo, self.hi, self.period);
        let clamp = |t: f64| match period {
            Some(p) => lo + (t - lo).rem_euclid(p),
            None => t.clamp(lo, hi),
        };
        let mut best = (lo, f64::INFINITY);
        for &(t, q) in &self.samples {
            let d = q.distance(p);
            if d < best.1 {
                best = (t, d);
            }
        }
        let mut t = best.0;
        for _ in 0..30 {
            let (c, d1, d2) = curve.d2(t);
            let f = (c - p).dot(d1);
            let df = d1.dot(d1) + (c - p).dot(d2);
            if df.abs() <= 1e-300 {
                break;
            }
            let step = t - f / df;
            let next = if period.is_some() { step } else { clamp(step) };
            let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
            t = next;
            if done {
                break;
            }
        }
        let t = clamp(t);
        let d = curve.point(t).distance(p);
        let t = if d < best.1 { t } else { best.0 };
        let point = curve.point(t);
        Projection { parameter: t, point, distance: point.distance(p) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::{BSplineSurface, KnotVector};
    use crate::geometry::parametric::ParametricSurface;

    /// A point near a bicubic saddle is found where it is nearest, from a grid or from a
    /// close start.
    #[test]
    fn nearest_point_on_a_bspline_surface() {
        let poles: Vec<Vec<Point>> = (0..4)
            .map(|i| (0..4).map(|j| Point::new(i as f64, j as f64, ((i as f64 - 1.5) * (j as f64 - 1.5)) * 0.3)).collect())
            .collect();
        let knots = || KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
        let surface = Surface::BSpline(BSplineSurface::new(poles, knots(), knots()));
        let (u, v) = (0.3, 0.8);
        let d = surface.d1(u, v);
        let normal = d.1.cross(d.2).normalized();
        let p = d.0 + normal * 0.05;
        for start in [None, Some((0.35, 0.75))] {
            let (fu, fv, dist) = nearest_on_surface(p, &surface, start);
            assert!((fu - u).abs() < 1e-9 && (fv - v).abs() < 1e-9, "{fu} {fv}");
            assert!((dist - 0.05).abs() < 1e-12);
        }
    }
}
