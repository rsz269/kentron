//! Intersections of curves and surfaces.

use std::f64::consts::TAU;

use crate::math::{Axis2, Direction, Frame, Point, Point2, Vector, Vector2};
use crate::tolerance::ANGULAR;

use super::curve::{Circle, Ellipse, Hyperbola, Line, Parabola};
use super::curve2::{Circle2, Curve2};
use super::surface::{Cylinder, Plane, Surface};
use super::{extrema, ParametricCurve, ParametricSurface};

/// A point where two 2D curves meet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveIntersection2 {
    pub point: Point2,
    /// Parameter on the first curve.
    pub t1: f64,
    /// Parameter on the second curve.
    pub t2: f64,
}

/// The intersection points of two 2D curves (within `tolerance`), sorted by parameter on
/// the first curve. Overlapping stretches are reported by their end points.
#[doc(alias = "Geom2dAPI_InterCurveCurve")]
pub fn curves2(a: &Curve2, b: &Curve2, tolerance: f64) -> Vec<CurveIntersection2> {
    let tolerance = tolerance.max(1e-12);
    let Some(((a0, a1), (b0, b1))) = clipped_domains(a, b, tolerance) else { return Vec::new() };
    let pa = Polyline2::new(a, a0, a1);
    let pb = Polyline2::new(b, b0, b1);
    let mut found: Vec<CurveIntersection2> = Vec::new();
    let mut accept = |s: f64, t: f64| {
        let (p, q) = (a.point(s), b.point(t));
        if p.distance(q) <= tolerance {
            found.push(CurveIntersection2 { point: p, t1: s, t2: t });
        }
    };
    // Close pairs of polyline segments, refined on the curves.
    for i in 0..pa.points.len() - 1 {
        let (p0, p1) = (pa.points[i], pa.points[i + 1]);
        let margin_a = pa.sag[i];
        for j in 0..pb.points.len() - 1 {
            let (q0, q1) = (pb.points[j], pb.points[j + 1]);
            let margin = margin_a + pb.sag[j] + tolerance;
            if p0.x.min(p1.x) - margin > q0.x.max(q1.x)
                || q0.x.min(q1.x) - margin > p0.x.max(p1.x)
                || p0.y.min(p1.y) - margin > q0.y.max(q1.y)
                || q0.y.min(q1.y) - margin > p0.y.max(p1.y)
            {
                continue;
            }
            let (u, v, d) = segment_closest(p0, p1, q0, q1);
            if d > margin {
                continue;
            }
            let s = pa.params[i] + u * (pa.params[i + 1] - pa.params[i]);
            let t = pb.params[j] + v * (pb.params[j + 1] - pb.params[j]);
            // The refined pair, unless refining went astray (at a cusp, say).
            let (rs, rt) = refine_crossing(a, b, s, t, (a0, a1), (b0, b1));
            if a.point(rs).distance(b.point(rt)) <= a.point(s).distance(b.point(t)) {
                accept(rs, rt);
            } else {
                accept(s, t);
            }
        }
    }
    // Ends of either curve lying on the other (touching ends, ends of overlaps).
    let ends = |c: &Curve2, lo: f64, hi: f64| {
        let mut e = Vec::new();
        if !c.is_closed() || !c.is_periodic() {
            for t in [lo, hi] {
                if t.is_finite() && t.abs() < 1e50 {
                    e.push(t);
                }
            }
        }
        e
    };
    for s in ends(a, a.first_parameter(), a.last_parameter()) {
        if let Some(x) = extrema::point_curve2(a.point(s), b).first() {
            accept(s, x.parameter);
        }
    }
    for t in ends(b, b.first_parameter(), b.last_parameter()) {
        if let Some(x) = extrema::point_curve2(b.point(t), a).first() {
            accept(x.parameter, t);
        }
    }
    simplify(a, b, found, tolerance)
}

/// The parameter windows of two curves that can meet: unbounded ends are cut where the
/// curve leaves the other curve's bounding box. `None` if the boxes are apart.
fn clipped_domains(a: &Curve2, b: &Curve2, tolerance: f64) -> Option<((f64, f64), (f64, f64))> {
    let bounded = |c: &Curve2| {
        let (lo, hi) = (c.first_parameter(), c.last_parameter());
        lo.is_finite() && hi.is_finite() && lo.abs() < 1e50 && hi.abs() < 1e50
    };
    let window = |c: &Curve2| {
        let (mut lo, mut hi) = (c.first_parameter(), c.last_parameter());
        if !(lo.is_finite() && lo.abs() < 1e50) {
            lo = -100.0;
        }
        if !(hi.is_finite() && hi.abs() < 1e50) {
            hi = 100.0;
        }
        (lo, hi)
    };
    let sampled_box = |c: &Curve2, (lo, hi): (f64, f64)| {
        let mut bb = crate::math::BoundingBox2::EMPTY;
        for i in 0..=256 {
            bb.add_point(c.point(lo + (hi - lo) * i as f64 / 256.0));
        }
        bb
    };
    let box_of = |c: &Curve2| {
        if bounded(c) {
            super::bounds::curve2(c, 0.0)
        } else {
            sampled_box(c, window(c))
        }
    };
    // Both unbounded: the box round both windows.
    let (box_a, box_b) = match (bounded(a), bounded(b)) {
        (false, false) => {
            let mut bb = box_of(a);
            bb.add_box(&box_of(b));
            (bb, bb)
        }
        _ => (box_of(a), box_of(b)),
    };
    let clip = |c: &Curve2, other: crate::math::BoundingBox2| -> Option<(f64, f64)> {
        let other = other.enlarged(tolerance + 1e-9 * (other.max - other.min).length());
        let (lo, hi) = (c.first_parameter(), c.last_parameter());
        let finite = |t: f64| t.is_finite() && t.abs() < 1e50;
        if finite(lo) && finite(hi) {
            return Some((lo, hi));
        }
        if let Curve2::Line(line) = c {
            // The stretch of the line inside the box.
            let Axis2 { origin: o, direction: d } = line.axis();
            let (mut t0, mut t1) = (f64::NEG_INFINITY, f64::INFINITY);
            for (oc, dc, mn, mx) in [(o.x, d.x(), other.min.x, other.max.x), (o.y, d.y(), other.min.y, other.max.y)] {
                if dc.abs() < 1e-300 {
                    if oc < mn || oc > mx {
                        return None;
                    }
                } else {
                    let (u, v) = ((mn - oc) / dc, (mx - oc) / dc);
                    t0 = t0.max(u.min(v));
                    t1 = t1.min(u.max(v));
                }
            }
            let (t0, t1) = (if finite(lo) { t0.max(lo) } else { t0 }, if finite(hi) { t1.min(hi) } else { t1 });
            return (t0 <= t1).then_some((t0, t1));
        }
        // Step outwards until the curve is well past the box.
        let center = other.min + 0.5 * (other.max - other.min);
        let radius = 0.5 * (other.max - other.min).length() + tolerance;
        let reach = |from: f64, sign: f64| {
            let mut step = 1.0;
            for _ in 0..80 {
                let t = from + sign * step;
                if c.point(t).distance(center) > 4.0 * radius + 1.0 {
                    return t;
                }
                step *= 2.0;
            }
            from + sign * step
        };
        let t0 = if finite(lo) { lo } else { reach(if finite(hi) { hi } else { 0.0 }, -1.0) };
        let t1 = if finite(hi) { hi } else { reach(if finite(lo) { lo } else { 0.0 }, 1.0) };
        Some((t0, t1))
    };
    Some((clip(a, box_b)?, clip(b, box_a)?))
}

/// A curve sampled densely enough that each chord stays within `sag` of its arc.
struct Polyline2 {
    params: Vec<f64>,
    points: Vec<Point2>,
    /// Per chord: how far the arc may stray from it.
    sag: Vec<f64>,
}

impl Polyline2 {
    fn new(c: &Curve2, lo: f64, hi: f64) -> Self {
        let mut cuts: Vec<f64> = vec![lo, hi];
        cuts.extend(c.breaks().into_iter().filter(|&k| lo < k && k < hi));
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let per_piece = (512 / (cuts.len() - 1)).clamp(16, 64);
        let mut params = vec![lo];
        for w in cuts.windows(2) {
            for i in 1..=per_piece {
                params.push(if i == per_piece { w[1] } else { w[0] + (w[1] - w[0]) * i as f64 / per_piece as f64 });
            }
        }
        let points: Vec<Point2> = params.iter().map(|&t| c.point(t)).collect();
        // The arc's distance from its chord, estimated from the midpoint (doubled, as a
        // margin for the parts between samples).
        let sag = (0..points.len() - 1)
            .map(|i| {
                let m = c.point(0.5 * (params[i] + params[i + 1]));
                let (_, _, d) = segment_closest(m, m, points[i], points[i + 1]);
                2.0 * d
            })
            .collect();
        Self { params, points, sag }
    }
}

/// The closest points of segments `p0p1` and `q0q1`, as fractions along each, and their
/// distance.
fn segment_closest(p0: Point2, p1: Point2, q0: Point2, q1: Point2) -> (f64, f64, f64) {
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (mut s, mut t);
    if a <= 1e-300 && e <= 1e-300 {
        return (0.0, 0.0, p0.distance(q0));
    }
    if a <= 1e-300 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = d1.dot(r);
        if e <= 1e-300 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            s = if denom > 1e-14 * a * e { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
            t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
        }
    }
    let d = (p0 + s * d1).distance(q0 + t * d2);
    (s, t, d)
}

/// Newton steps towards `a(s) = b(t)`, falling back to minimising the distance where the
/// curves are tangent; parameters stay inside the windows.
fn refine_crossing(a: &Curve2, b: &Curve2, mut s: f64, mut t: f64, wa: (f64, f64), wb: (f64, f64)) -> (f64, f64) {
    let clamp = |x: f64, (lo, hi): (f64, f64)| x.clamp(lo, hi);
    let gap = |s: f64, t: f64| a.point(s).distance(b.point(t));
    let mut best = (s, t, gap(s, t));
    for _ in 0..60 {
        let (p, da, dda) = a.d2(s);
        let (q, db, ddb) = b.d2(t);
        let r = p - q;
        let det = da.x * (-db.y) - da.y * (-db.x);
        let scale = da.length() * db.length();
        let (ns, nt) = if det.abs() > 1e-10 * scale {
            // Solve [da, -db] (ds, dt) = -r.
            let ds = (-r.x * (-db.y) + r.y * (-db.x)) / det;
            let dt = (da.x * (-r.y) - da.y * (-r.x)) / det;
            (s + ds, t + dt)
        } else {
            // Tangent: Newton on the gradient of |r|^2 / 2.
            let g = (r.dot(da), -r.dot(db));
            let h11 = da.dot(da) + r.dot(dda);
            let h22 = db.dot(db) - r.dot(ddb);
            let h12 = -da.dot(db);
            let hdet = h11 * h22 - h12 * h12;
            if hdet.abs() <= 1e-300 {
                break;
            }
            (s - (h22 * g.0 - h12 * g.1) / hdet, t - (h11 * g.1 - h12 * g.0) / hdet)
        };
        let (ns, nt) = (clamp(ns, wa), clamp(nt, wb));
        let d = gap(ns, nt);
        if d < best.2 {
            best = (ns, nt, d);
        }
        let moved = (ns - s).abs() + (nt - t).abs();
        (s, t) = (ns, nt);
        if moved <= 1e-15 * (1.0 + s.abs() + t.abs()) || d == 0.0 {
            break;
        }
    }
    (best.0, best.1)
}

/// Merges duplicates and reduces overlapping stretches to their ends.
fn simplify(a: &Curve2, b: &Curve2, mut found: Vec<CurveIntersection2>, tolerance: f64) -> Vec<CurveIntersection2> {
    let gap = |x: &CurveIntersection2| a.point(x.t1).distance(b.point(x.t2));
    // Best first, so that a duplicate keeps the most accurate point.
    found.sort_by(|x, y| gap(x).total_cmp(&gap(y)));
    let mut kept: Vec<CurveIntersection2> = Vec::new();
    for x in found {
        if !kept.iter().any(|k| k.point.distance(x.point) <= tolerance) {
            kept.push(x);
        }
    }
    kept.sort_by(|x, y| x.t1.total_cmp(&y.t1));
    // Consecutive points joined by a stretch of `a` lying on `b` bound an overlap: only
    // the stretch's ends stay.
    let on_b = |s: f64, hint: f64| {
        let p = a.point(s);
        let (lo, hi) = (b.first_parameter(), b.last_parameter());
        let mut t = hint;
        for _ in 0..30 {
            let (q, d1, d2) = b.d2(t);
            let r = q - p;
            let f = r.dot(d1);
            let df = d1.dot(d1) + r.dot(d2);
            if df.abs() <= 1e-300 {
                break;
            }
            let next = (t - f / df).clamp(lo, hi);
            let done = (next - t).abs() <= 1e-14 * (1.0 + t.abs());
            t = next;
            if done {
                break;
            }
        }
        b.point(t).distance(p) <= tolerance
    };
    let joined = |x: &CurveIntersection2, y: &CurveIntersection2| {
        (1..8).all(|k| {
            let f = k as f64 / 8.0;
            on_b(x.t1 + f * (y.t1 - x.t1), x.t2 + f * (y.t2 - x.t2))
        })
    };
    let mut out: Vec<CurveIntersection2> = Vec::new();
    let mut i = 0;
    while i < kept.len() {
        let mut j = i;
        while j + 1 < kept.len() && joined(&kept[j], &kept[j + 1]) {
            j += 1;
        }
        out.push(kept[i]);
        if j > i {
            out.push(kept[j]);
        }
        i = j + 1;
    }
    out
}

/// The intersection points of two circles, computed exactly (no tolerance), sorted by
/// parameter on the first circle. Tangent circles give one point; coincident circles
/// none.
#[doc(alias = "IntAna2d_AnaIntersection")]
pub fn circles2(a: &Circle2, b: &Circle2) -> Vec<CurveIntersection2> {
    let (c1, c2) = (a.center(), b.center());
    let (r1, r2) = (a.radius(), b.radius());
    let between = c2 - c1;
    let d = between.length();
    if d == 0.0 {
        return Vec::new();
    }
    let u = between * (1.0 / d);
    let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let h2 = r1 * r1 - along * along;
    let scale = r1.max(r2).max(d);
    let eps = (1e-15 * scale) * (1e-15 * scale) * 64.0 + 4.0 * f64::EPSILON * r1 * r1;
    if h2 < -eps {
        return Vec::new();
    }
    let mut points = Vec::new();
    let base = c1 + along * u;
    if h2 <= eps {
        points.push(base);
    } else {
        let h = h2.sqrt();
        let n = Vector2::new(-u.y, u.x);
        points.push(base + h * n);
        points.push(base + (-h) * n);
    }
    let angle = |c: &Circle2, p: Point2| {
        let f = c.frame();
        let v = p - f.origin();
        let t = v.dot(Vector2::from(f.y())).atan2(v.dot(Vector2::from(f.x())));
        let t = t.rem_euclid(TAU);
        if t >= TAU { 0.0 } else { t }
    };
    let mut out: Vec<CurveIntersection2> =
        points.into_iter().map(|p| CurveIntersection2 { point: p, t1: angle(a, p), t2: angle(b, p) }).collect();
    out.sort_by(|x, y| x.t1.total_cmp(&y.t1));
    out
}

/// The intersection of two elementary surfaces (planes, cylinders, cones, spheres, tori),
/// when it has a closed form.
#[doc(alias = "IntAna_QuadQuadGeo")]
#[derive(Debug, Clone, PartialEq)]
pub enum QuadricIntersection {
    Empty,
    /// The surfaces coincide.
    Same,
    /// Isolated (tangency) points.
    Points(Vec<Point>),
    Lines(Vec<Line>),
    Circles(Vec<Circle>),
    Ellipses(Vec<Ellipse>),
    Parabolas(Vec<Parabola>),
    Hyperbolas(Vec<Hyperbola>),
    /// The intersection has no closed form; see [`cylinder_quadric`].
    NotAnalytic,
}

/// The intersection of two elementary surfaces; distances below `tolerance` count as
/// zero (so nearly tangent cylinders meet along a single line).
#[doc(alias = "IntAna_QuadQuadGeo")]
pub fn quadrics(a: &Surface, b: &Surface, tolerance: f64) -> QuadricIntersection {
    use QuadricIntersection as Q;
    let tol = tolerance;
    let parallel = |u: Vector, v: Vector| u.cross(v).length() <= ANGULAR.sqrt();
    match (a, b) {
        (Surface::Plane(p1), Surface::Plane(p2)) => {
            let (n1, n2) = (Vector::from(p1.normal()), Vector::from(p2.normal()));
            if parallel(n1, n2) {
                return if p1.distance(p2.origin()) <= tol { Q::Same } else { Q::Empty };
            }
            let n = n1.cross(n2);
            let (d1, d2) = (n1.dot(p1.origin() - Point::ORIGIN), n2.dot(p2.origin() - Point::ORIGIN));
            let origin = Point::ORIGIN + (n2.cross(n) * d1 + n.cross(n1) * d2) * (1.0 / n.length_squared());
            Q::Lines(vec![Line::new(origin, Direction::from_vector(n))])
        }
        (Surface::Plane(plane), Surface::Cylinder(cylinder)) | (Surface::Cylinder(cylinder), Surface::Plane(plane)) => {
            plane_cylinder(plane, cylinder, tol)
        }
        (Surface::Plane(plane), Surface::Sphere(sphere)) | (Surface::Sphere(sphere), Surface::Plane(plane)) => {
            let n = Vector::from(plane.normal());
            let c = sphere.frame().origin();
            let d = n.dot(c - plane.origin());
            let r = sphere.radius();
            if d.abs() > r + tol {
                Q::Empty
            } else if d.abs() >= r - tol {
                Q::Points(vec![c - n * d])
            } else {
                let center = c - n * d;
                Q::Circles(vec![Circle::new(Frame::from_z(center, plane.normal()), (r * r - d * d).sqrt())])
            }
        }
        (Surface::Sphere(s1), Surface::Sphere(s2)) => {
            let (c1, c2) = (s1.frame().origin(), s2.frame().origin());
            let (r1, r2) = (s1.radius(), s2.radius());
            let d = c1.distance(c2);
            if d <= tol {
                return if (r1 - r2).abs() <= tol { Q::Same } else { Q::Empty };
            }
            if d > r1 + r2 + tol || d < (r1 - r2).abs() - tol {
                return Q::Empty;
            }
            let axis = (c2 - c1) * (1.0 / d);
            let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
            let center = c1 + axis * along;
            let radius2 = r1 * r1 - along * along;
            if radius2 <= tol * tol {
                return Q::Points(vec![center]);
            }
            Q::Circles(vec![Circle::new(Frame::from_z(center, Direction::from_vector(axis)), radius2.sqrt())])
        }
        (Surface::Cylinder(c1), Surface::Cylinder(c2)) => {
            let (z1, z2) = (Vector::from(c1.frame().z()), Vector::from(c2.frame().z()));
            if !parallel(z1, z2) {
                // Skew axes: apart, or touching at one point on their common perpendicular.
                let (o1, o2) = (c1.frame().origin(), c2.frame().origin());
                let normal = z1.cross(z2);
                let w = o2 - o1;
                let gap = w.dot(normal) / normal.length();
                let (r1, r2) = (c1.radius(), c2.radius());
                if gap.abs() > r1 + r2 + tol {
                    return Q::Empty;
                }
                if (gap.abs() - (r1 + r2)).abs() <= tol {
                    // The feet of the common perpendicular on each axis.
                    let (a, b, c) = (z1.dot(z1), z1.dot(z2), z2.dot(z2));
                    let (d, e) = (z1.dot(w), z2.dot(w));
                    let denom = a * c - b * b;
                    let s = (d * c - b * e) / denom;
                    let t = (b * d - a * e) / denom;
                    let (p1, p2) = (o1 + z1 * s, o2 + z2 * t);
                    let across = p2 - p1;
                    return Q::Points(vec![p1 + across * (r1 / (r1 + r2))]);
                }
                return equal_cylinders(c1, c2, tol).unwrap_or(Q::NotAnalytic);
            }
            // Parallel axes: intersect the two circles in a plane across the axes.
            let o1 = c1.frame().origin();
            let mut o2 = c2.frame().origin() - o1;
            o2 = o2 - z1 * o2.dot(z1);
            let d = o2.length();
            let (r1, r2) = (c1.radius(), c2.radius());
            if d <= tol {
                return if (r1 - r2).abs() <= tol { Q::Same } else { Q::Empty };
            }
            if d > r1 + r2 + tol || d < (r1 - r2).abs() - tol {
                return Q::Empty;
            }
            let x = o2 * (1.0 / d);
            let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
            let across2 = r1 * r1 - along * along;
            let y = z1.cross(x);
            let direction = c1.frame().z();
            let foot = o1 + x * along;
            if across2 <= tol * tol {
                return Q::Lines(vec![Line::new(foot, direction)]);
            }
            let across = across2.sqrt();
            Q::Lines(vec![Line::new(foot + y * across, direction), Line::new(foot - y * across, direction)])
        }
        (Surface::Cylinder(cylinder), Surface::Sphere(sphere)) | (Surface::Sphere(sphere), Surface::Cylinder(cylinder)) => {
            // Only a sphere centered on the axis meets the cylinder in circles.
            let f = cylinder.frame();
            let z = Vector::from(f.z());
            let d = sphere.frame().origin() - f.origin();
            let height = d.dot(z);
            if (d - z * height).length() > tol {
                return Q::NotAnalytic;
            }
            let (r, big) = (cylinder.radius(), sphere.radius());
            if r > big + tol {
                return Q::Empty;
            }
            let center = f.origin() + z * height;
            if r >= big - tol {
                return Q::Circles(vec![Circle::new(Frame::new(center, f.z(), f.x()), r)]);
            }
            let h = (big * big - r * r).sqrt();
            Q::Circles(vec![
                Circle::new(Frame::new(center + z * h, f.z(), f.x()), r),
                Circle::new(Frame::new(center - z * h, f.z(), f.x()), r),
            ])
        }
        (Surface::Plane(plane), Surface::Cone(cone)) | (Surface::Cone(cone), Surface::Plane(plane)) => {
            plane_cone(plane, cone, tol)
        }
        (Surface::Cone(c1), Surface::Cone(c2)) if coaxial_cones(c1, c2, tol).is_some() => coaxial_cones(c1, c2, tol).expect("checked"),
        (Surface::Cone(c1), Surface::Cone(c2)) if c1.apex().distance(c2.apex()) <= tol => common_apex_cones(c1, c2, tol),
        (Surface::Cone(_), _) | (_, Surface::Cone(_)) => {
            // Coincident surfaces have proportional equations; other meetings are curves
            // of higher degree.
            match (quadric_of(a), quadric_of(b)) {
                (Some(qa), Some(qb)) if proportional(&qa, &qb, tol) => Q::Same,
                _ => Q::NotAnalytic,
            }
        }
        (a, b) => todo!("intersecting a {} and a {}", surface_kind(a), surface_kind(b)),
    }
}

/// Two cones on one axis meet in circles: where their meridians (lines in the plane of
/// height along the axis and signed radius) give the same distance from the axis. `None`
/// if the axes differ.
fn coaxial_cones(c1: &crate::geometry::surface::Cone, c2: &crate::geometry::surface::Cone, tol: f64) -> Option<QuadricIntersection> {
    let (f1, f2) = (c1.frame(), c2.frame());
    let z = Vector::from(f1.z());
    let z2 = Vector::from(f2.z());
    let offset = f2.origin() - f1.origin();
    if z.cross(z2).length() > ANGULAR || (offset - z * offset.dot(z)).length() > tol {
        return None;
    }
    // Signed radius as a function of height h along z from f1's origin: r = a + b h.
    let line = |c: &crate::geometry::surface::Cone, h0: f64, along: f64| {
        let (sin, cos) = c.semi_angle_sin_cos();
        let b = along * sin / cos;
        (c.radius() - b * h0, b)
    };
    let (a1, b1) = line(c1, 0.0, 1.0);
    let (a2, b2) = line(c2, offset.dot(z), z2.dot(z).signum());
    let mut heights: Vec<f64> = Vec::new();
    // The same radius, or opposite ones (the other nappe).
    for (db, da) in [(b1 - b2, a2 - a1), (b1 + b2, -a2 - a1)] {
        if db.abs() <= 1e-12 * (b1.abs() + b2.abs() + 1.0) {
            // Parallel meridians: the same cone, or none of it.
            if da.abs() <= tol {
                return Some(QuadricIntersection::Same);
            }
            continue;
        }
        let h = da / db;
        if !heights.iter().any(|&x| (x - h).abs() <= tol) {
            heights.push(h);
        }
    }
    let mut circles = Vec::new();
    let mut points = Vec::new();
    for h in heights {
        let center = f1.origin() + z * h;
        let r = (a1 + b1 * h).abs();
        if r <= tol {
            points.push(center);
        } else {
            circles.push(Circle::new(Frame::new(center, f1.z(), f1.x()), r));
        }
    }
    Some(match (circles.is_empty(), points.is_empty()) {
        (true, true) => QuadricIntersection::Empty,
        (false, true) => QuadricIntersection::Circles(circles),
        (true, false) => QuadricIntersection::Points(points),
        (false, false) => QuadricIntersection::Circles(circles),
    })
}

/// Two cones (on different axes) with one apex meet in lines through it: their common
/// generators, at the angle of each cone from its axis (on either nappe), at most four.
fn common_apex_cones(c1: &crate::geometry::surface::Cone, c2: &crate::geometry::surface::Cone, tol: f64) -> QuadricIntersection {
    let apex = c1.apex();
    let (z1, z2) = (Vector::from(c1.frame().z()), Vector::from(c2.frame().z()));
    let k = z1.dot(z2);
    let across = z1.cross(z2).normalized();
    let p = c1.semi_angle_sin_cos().1.abs();
    let mut lines: Vec<Line> = Vec::new();
    // A line's direction d (one of ±d, so d·z1 = cos α1): d·z2 = ±cos α2, its part in the
    // axes' plane a z1 + b z2 from those, the rest across.
    for q in [c2.semi_angle_sin_cos().1.abs(), -c2.semi_angle_sin_cos().1.abs()] {
        let det = 1.0 - k * k;
        let (a, b) = ((p - q * k) / det, (q - p * k) / det);
        let in_plane = z1 * a + z2 * b;
        let rest = 1.0 - in_plane.length_squared();
        // (Nearly none across: cones touching along one generator.)
        let spread = if rest.abs() <= tol * tol { vec![0.0] } else if rest > 0.0 { vec![rest.sqrt(), -rest.sqrt()] } else { Vec::new() };
        for g in spread {
            let d = in_plane + across * g;
            let d = Direction::from_vector(d);
            if !lines.iter().any(|l| Vector::from(l.direction).cross(Vector::from(d)).length() <= ANGULAR) {
                lines.push(Line::new(apex, d));
            }
        }
    }
    if lines.is_empty() { QuadricIntersection::Points(vec![apex]) } else { QuadricIntersection::Lines(lines) }
}

fn surface_kind(s: &Surface) -> String {
    format!("{s:?}").split('(').next().unwrap_or_default().to_lowercase()
}

/// Cylinders of equal radius whose axes cross meet in two ellipses, in the planes through
/// the crossing point bisecting the axes; `None` for other crossing cylinders.
fn equal_cylinders(c1: &Cylinder, c2: &Cylinder, tol: f64) -> Option<QuadricIntersection> {
    let r = c1.radius();
    if (r - c2.radius()).abs() > tol {
        return None;
    }
    let (a, b) = (Vector::from(c1.frame().z()), Vector::from(c2.frame().z()));
    let (o1, o2) = (c1.frame().origin(), c2.frame().origin());
    // The axes' closest points.
    let n = a.cross(b);
    let w = o1 - o2;
    if n.dot(w).abs() / n.length() > tol {
        return None;
    }
    let (aa, ab, bb, aw, bw) = (a.dot(a), a.dot(b), b.dot(b), a.dot(w), b.dot(w));
    let den = aa * bb - ab * ab;
    let s = (ab * bw - bb * aw) / den;
    let center = o1 + a * s;
    let across = Direction::from_vector(n);
    let cos = a.dot(b);
    let mut ellipses = Vec::new();
    for (normal, cut) in [(a - b, (0.5 * (1.0 - cos)).sqrt()), (a + b, (0.5 * (1.0 + cos)).sqrt())] {
        let normal = Direction::from_vector(normal);
        // The major axis, in the plane across the minor one (along a × b).
        let major = Vector::from(normal).cross(Vector::from(across));
        let frame = Frame::new(center, normal, Direction::from_vector(major));
        ellipses.push(Ellipse::new(frame, r / cut, r));
    }
    Some(QuadricIntersection::Ellipses(ellipses))
}

/// True if two quadric equations are multiples of each other.
fn proportional(a: &[f64; 10], b: &[f64; 10], tol: f64) -> bool {
    let (ia, ib) = (
        (0..10).max_by(|&i, &j| a[i].abs().total_cmp(&a[j].abs())).expect("coefficients"),
        (0..10).max_by(|&i, &j| b[i].abs().total_cmp(&b[j].abs())).expect("coefficients"),
    );
    if ia != ib || a[ia] == 0.0 {
        return false;
    }
    let k = b[ib] / a[ia];
    (0..10).all(|i| (a[i] * k - b[i]).abs() <= tol * b[ib].abs().max(1.0))
}

/// A plane meets a cone in a circle (across the axis), in generators (through the
/// apex), or in an ellipse, parabola or hyperbola (left to [`ruled_quadric`]).
fn plane_cone(plane: &Plane, cone: &crate::geometry::surface::Cone, tol: f64) -> QuadricIntersection {
    use QuadricIntersection as Q;
    let n = Vector::from(plane.normal());
    let f = cone.frame();
    let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
    let (sin, cos) = cone.semi_angle_sin_cos();
    let apex = cone.apex();
    if n.dot(apex - plane.origin()).abs() <= tol {
        // Generators lying in the plane: n·D(θ) = 0 with D = sin α (cos θ x + sin θ y) + cos α z.
        let (a, b, c) = (n.dot(x) * sin, n.dot(y) * sin, -n.dot(z) * cos);
        let r = a.hypot(b);
        if r <= ANGULAR || c.abs() > r + ANGULAR {
            return Q::Points(vec![apex]);
        }
        let base = b.atan2(a);
        let spread = (c / r).clamp(-1.0, 1.0).acos();
        let mut lines = Vec::new();
        for theta in [base + spread, base - spread] {
            let d = (x * theta.cos() + y * theta.sin()) * sin + z * cos;
            if lines.iter().all(|l: &Line| Vector::from(l.direction).cross(d).length() > ANGULAR.sqrt()) {
                lines.push(Line::new(apex, Direction::from_vector(d)));
            }
        }
        return Q::Lines(lines);
    }
    if n.cross(z).length() <= ANGULAR.sqrt() {
        // Across the axis: a circle at the height where the plane meets the axis.
        let h = n.dot(plane.origin() - f.origin()) / n.dot(z);
        let radius = (cone.radius() + h * sin / cos).abs();
        let center = f.origin() + z * h;
        return Q::Circles(vec![Circle::new(Frame::new(center, f.z(), f.x()), radius)]);
    }
    Q::NotAnalytic
}

/// A plane meets a cylinder in a circle (axis across the plane), lines (axis parallel to
/// it) or an ellipse.
fn plane_cylinder(plane: &Plane, cylinder: &Cylinder, tol: f64) -> QuadricIntersection {
    use QuadricIntersection as Q;
    let n = Vector::from(plane.normal());
    let f = cylinder.frame();
    let z = Vector::from(f.z());
    let r = cylinder.radius();
    let cos = n.dot(z);
    if cos.abs() <= ANGULAR.sqrt() {
        // The axis is parallel to the plane: generators at the plane's distance.
        let d = n.dot(f.origin() - plane.origin());
        if d.abs() > r + tol {
            return Q::Empty;
        }
        let foot = f.origin() - n * d;
        if d.abs() >= r - tol {
            return Q::Lines(vec![Line::new(foot, f.z())]);
        }
        let across = z.cross(n) * (r * r - d * d).sqrt();
        return Q::Lines(vec![Line::new(foot + across, f.z()), Line::new(foot - across, f.z())]);
    }
    // Where the axis crosses the plane.
    let center = f.origin() + z * (n.dot(plane.origin() - f.origin()) / cos);
    if (cos.abs() - 1.0).abs() <= ANGULAR.sqrt() {
        return Q::Circles(vec![Circle::new(Frame::new(center, f.z(), f.x()), r)]);
    }
    // The minor axis is across the cylinder's axis, the major one along the plane's
    // steepest slope with respect to it.
    let minor = z.cross(n).normalized();
    let major = n.cross(minor);
    let frame = Frame::new(center, plane.normal(), Direction::from_vector(major));
    Q::Ellipses(vec![Ellipse::new(frame, r / cos.abs(), r)])
}

/// A branch of the intersection of a ruled elementary surface (cylinder or cone: points
/// `A(u) + v·D(u)`) with a quadric, parametrized by the angle `u` around its axis.
#[doc(alias = "IntAna_Curve")]
#[derive(Debug, Clone, PartialEq)]
pub struct QuadricCurve {
    base: Surface,
    /// The other quadric in world coordinates: `A x² + B y² + C z² + 2D xy + 2E xz +
    /// 2F yz + 2G x + 2H y + 2I z + J = 0`.
    quadric: [f64; 10],
    first: f64,
    last: f64,
    upper: bool,
}

/// The value of a quadric (as in [`QuadricCurve`]) at `p`.
fn quadric_value(q: &[f64; 10], p: Point) -> f64 {
    let (x, y, z) = (p.x, p.y, p.z);
    q[0] * x * x + q[1] * y * y + q[2] * z * z + 2.0 * (q[3] * x * y + q[4] * x * z + q[5] * y * z)
        + 2.0 * (q[6] * x + q[7] * y + q[8] * z)
        + q[9]
}

/// The implicit equation of an elementary quadric surface in world coordinates.
fn quadric_of(surface: &Surface) -> Option<[f64; 10]> {
    // From a symmetric matrix M, a vector g and a constant: pᵀMp + 2g·p + c.
    let from = |m: [[f64; 3]; 3], g: Vector, c: f64| {
        [m[0][0], m[1][1], m[2][2], m[0][1], m[0][2], m[1][2], g.x, g.y, g.z, c]
    };
    // |p - o|² - k (z·(p - o))², with the identity minus k z zᵀ as M.
    let axial = |o: Point, z: Vector, k: f64, r2: f64| {
        let m = [
            [1.0 - k * z.x * z.x, -k * z.x * z.y, -k * z.x * z.z],
            [-k * z.y * z.x, 1.0 - k * z.y * z.y, -k * z.y * z.z],
            [-k * z.z * z.x, -k * z.z * z.y, 1.0 - k * z.z * z.z],
        ];
        let o = o - Point::ORIGIN;
        let mo = Vector::new(
            m[0][0] * o.x + m[0][1] * o.y + m[0][2] * o.z,
            m[1][0] * o.x + m[1][1] * o.y + m[1][2] * o.z,
            m[2][0] * o.x + m[2][1] * o.y + m[2][2] * o.z,
        );
        from(m, -mo, o.dot(mo) - r2)
    };
    Some(match surface {
        Surface::Plane(plane) => {
            let n = Vector::from(plane.normal());
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, n.x / 2.0, n.y / 2.0, n.z / 2.0, -n.dot(plane.origin() - Point::ORIGIN)]
        }
        Surface::Sphere(s) => axial(s.frame().origin(), Vector::from(s.frame().z()), 0.0, s.radius() * s.radius()),
        Surface::Cylinder(c) => axial(c.frame().origin(), Vector::from(c.frame().z()), 1.0, c.radius() * c.radius()),
        Surface::Cone(c) => {
            // |d|² = (d·z)² / cos² α for d from the apex.
            let cos = c.semi_angle_sin_cos().1;
            axial(c.apex(), Vector::from(c.frame().z()), 1.0 / (cos * cos), 0.0)
        }
        _ => return None,
    })
}

impl QuadricCurve {
    /// The angular range `(first, last)` of the branch.
    pub fn domain(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// Whether this is the branch with the larger generator parameter `v`.
    pub fn is_upper(&self) -> bool {
        self.upper
    }

    /// The coefficients `(a, b, c)` of `a v² + b v + c = 0` for the generator parameter at
    /// `angle` (the quadric along the generator is quadratic in `v`).
    fn height_equation(&self, angle: f64) -> (f64, f64, f64) {
        let q = |v: f64| quadric_value(&self.quadric, self.base.point(angle, v));
        let (minus, zero, plus) = (q(-1.0), q(0.0), q(1.0));
        // A plane is met along a generator linearly: the rounding left in `a` is noise (where
        // the generator runs nearly along the plane it can outweigh `b`).
        let planar = self.quadric[..6].iter().all(|&c| c == 0.0);
        (if planar { 0.0 } else { 0.5 * (plus + minus) - zero }, 0.5 * (plus - minus), zero)
    }

    /// The branch's generator parameter at `angle`; `None` where it does not exist.
    fn height(&self, angle: f64) -> Option<f64> {
        let (a, b, c) = self.height_equation(angle);
        let scale = a.abs().max(b.abs()).max(c.abs()).max(1e-300);
        if a.abs() <= 1e-12 * scale {
            return (b.abs() > 1e-12 * scale).then(|| -c / b);
        }
        // (Rounding in the discriminant goes with the size of its terms.)
        let disc = b * b - 4.0 * a * c;
        if disc < -1e-12 * (b * b).max(4.0 * (a * c).abs()) {
            return None;
        }
        let root = disc.max(0.0).sqrt();
        let sign = if self.upper == (a > 0.0) { 1.0 } else { -1.0 };
        Some((-b + sign * root) / (2.0 * a))
    }

    /// The point of the branch at `angle`, `None` outside its domain.
    pub fn point(&self, angle: f64) -> Option<Point> {
        let v = self.height(angle)?;
        Some(self.base.point(angle, v))
    }

    /// The point at `angle` within the domain, robust at its ends (where the two
    /// branches meet and rounding may make the discriminant slightly negative).
    pub fn point_clamped(&self, angle: f64) -> Point {
        let angle = angle.clamp(self.first, self.last);
        let (a, b, c) = self.height_equation(angle);
        let scale = a.abs().max(b.abs()).max(c.abs()).max(1e-300);
        let v = if a.abs() <= 1e-12 * scale {
            -c / b
        } else {
            let root = (b * b - 4.0 * a * c).max(0.0).sqrt();
            let sign = if self.upper == (a > 0.0) { 1.0 } else { -1.0 };
            (-b + sign * root) / (2.0 * a)
        };
        self.base.point(angle, v)
    }
}

/// The intersection of a cylinder or cone with a quadric as parametrized branches.
#[doc(alias = "IntAna_IntQuadQuad")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QuadricCurves {
    pub curves: Vec<QuadricCurve>,
    /// For each curve, the curve continuing it at its end, if any, and whether that one
    /// runs in the opposite direction.
    pub next: Vec<Option<(usize, bool)>>,
    /// For each curve, the curves continuing it at its start and at its end (running the
    /// opposite way, from the same angle); `None` at an end where the branch stops (beyond
    /// the allowed range of `v`) or for a branch all the way round.
    pub joins: Vec<[Option<usize>; 2]>,
    /// True if the surfaces coincide (then `curves` is empty).
    pub same: bool,
}

/// The intersection of `cylinder` with a plane, cylinder, cone or sphere as parametrized
/// branches, linked end to end. `None` for other surfaces.
#[doc(alias = "IntAna_IntQuadQuad")]
pub fn cylinder_quadric(cylinder: &Cylinder, other: &Surface, tolerance: f64) -> Option<QuadricCurves> {
    ruled_quadric(&Surface::Cylinder(*cylinder), other, tolerance, None)
}

/// The intersection of a cylinder or cone (`base`) with a plane, cylinder, cone or sphere
/// as branches parametrized by the angle around the base's axis, keeping only points with
/// the generator parameter `v` in `v_range` (if given: a cone meets a plane parallel to one
/// of its generators in a curve running off to infinity). `None` for other surfaces.
#[doc(alias = "IntAna_IntQuadQuad")]
pub fn ruled_quadric(base: &Surface, other: &Surface, tolerance: f64, v_range: Option<(f64, f64)>) -> Option<QuadricCurves> {
    if !matches!(base, Surface::Cylinder(_) | Surface::Cone(_)) {
        return None;
    }
    let quadric = quadric_of(other)?;
    if let QuadricIntersection::Same = quadrics(base, other, tolerance) {
        return Some(QuadricCurves { curves: Vec::new(), next: Vec::new(), joins: Vec::new(), same: true });
    }
    let branch = |first: f64, last: f64, upper: bool| QuadricCurve { base: base.clone(), quadric, first, last, upper };
    let probe = branch(0.0, TAU, true);
    let linear = (0..8).all(|i| {
        let (a, b, c) = probe.height_equation(i as f64);
        a.abs() <= 1e-12 * a.abs().max(b.abs()).max(c.abs()).max(1e-300)
    });
    let (v_lo, v_hi) = v_range.unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
    // Where a branch exists, and why it ends: the discriminant vanishing (the branches
    // join) or `v` leaving its range.
    let exists = |upper: bool, angle: f64| branch(0.0, TAU, upper).height(angle).filter(|v| *v >= v_lo && *v <= v_hi).is_some();
    let joins_here = |angle: f64| {
        let (a, b, c) = probe.height_equation(angle);
        !linear && b * b - 4.0 * a * c <= 1e-6 * (b * b).max(4.0 * (a * c).abs())
    };
    const SAMPLES: usize = 1440;
    let angles: Vec<f64> = (0..SAMPLES).map(|i| TAU * i as f64 / SAMPLES as f64).collect();
    let refine = |upper: bool, inside: f64, outside: f64| {
        let (mut a, mut b) = (inside, outside);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            if exists(upper, mid) {
                a = mid;
            } else {
                b = mid;
            }
        }
        a
    };
    // Between the samples at `angle` and `angle + step`, both out of range, where `v` runs
    // across the whole range: the angles where it enters and leaves it.
    let crossing = |upper: bool, angle: f64, step: f64| -> Option<(f64, f64)> {
        let height = |t: f64| branch(0.0, TAU, upper).height(t);
        let (h0, h1) = (height(angle)?, height(angle + step)?);
        if !((h0 < v_lo && h1 > v_hi) || (h0 > v_hi && h1 < v_lo)) {
            return None;
        }
        let root = |level: f64| -> Option<f64> {
            let (mut a, mut b) = (angle, angle + step);
            let below = |t: f64| height(t).map(|h| (h < level) == (h0 < level));
            for _ in 0..100 {
                let mid = 0.5 * (a + b);
                if below(mid)? {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            Some(0.5 * (a + b))
        };
        let (enter, leave) = if h0 < v_lo { (root(v_lo)?, root(v_hi)?) } else { (root(v_hi)?, root(v_lo)?) };
        // (Not across an asymptote, where `v` runs off to infinity and back from the other.)
        let lands = |t: f64, level: f64| height(t).is_some_and(|h| (h - level).abs() <= 1e-6 * (1.0 + level.abs() + (v_hi - v_lo)));
        let (enter_level, leave_level) = if h0 < v_lo { (v_lo, v_hi) } else { (v_hi, v_lo) };
        (leave > enter && lands(enter, enter_level) && lands(leave, leave_level)).then_some((enter, leave))
    };
    // Each branch's domains: (start, end, start joins, end joins).
    let mut domains: Vec<(bool, f64, f64, bool, bool)> = Vec::new();
    for upper in if linear { vec![true] } else { vec![true, false] } {
        let valid: Vec<bool> = angles.iter().map(|&t| exists(upper, t)).collect();
        if valid.iter().all(|&v| v) {
            domains.push((upper, 0.0, TAU, false, false));
            continue;
        }
        // Walk round from an invalid sample (where none is valid, for domains between them).
        let start = valid.iter().position(|&v| !v).expect("an invalid sample");
        let step = TAU / SAMPLES as f64;
        let mut i = 0;
        while i < SAMPLES {
            let k = (start + i) % SAMPLES;
            let next = (k + 1) % SAMPLES;
            if !valid[k] && valid[next] {
                // A domain starts between k and next; find where it ends.
                let base_angle = angles[start] + step * i as f64;
                let a = refine(upper, base_angle + step, base_angle);
                let mut j = i + 1;
                while valid[(start + j) % SAMPLES] {
                    j += 1;
                }
                let end_angle = angles[start] + step * j as f64;
                let b = refine(upper, end_angle - step, end_angle);
                domains.push((upper, a, b, joins_here(a), joins_here(b)));
                i = j;
            } else if let Some(domain) = crossing(upper, angles[start] + step * i as f64, step) {
                // A domain narrower than the samples: `v` running across its whole range
                // between two of them (a plane crossing a narrow band nearly along a
                // generator).
                domains.push((upper, domain.0, domain.1, false, false));
                i += 1;
            } else {
                i += 1;
            }
        }
    }
    let mut curves = Vec::new();
    let mut joins: Vec<[Option<usize>; 2]> = Vec::new();
    for &(upper, a, b, _, _) in &domains {
        curves.push(branch(a, b, upper));
        joins.push([None, None]);
    }
    // Pair joining ends: the other branch ending at the same angle.
    let same_angle = |x: f64, y: f64| ((x - y).rem_euclid(TAU)).min((y - x).rem_euclid(TAU)) <= 1e-6;
    for i in 0..domains.len() {
        for (end, joined) in [(0usize, domains[i].3), (1usize, domains[i].4)] {
            if !joined || joins[i][end].is_some() {
                continue;
            }
            let angle = if end == 0 { domains[i].1 } else { domains[i].2 };
            let partner = (0..domains.len()).find_map(|j| {
                if j == i || domains[j].0 == domains[i].0 {
                    return None;
                }
                [(0usize, domains[j].1, domains[j].3), (1usize, domains[j].2, domains[j].4)]
                    .into_iter()
                    .find(|&(_, t, joined)| joined && same_angle(t, angle))
                    .map(|(e, _, _)| (j, e))
            });
            if let Some((j, e)) = partner {
                joins[i][end] = Some(j);
                joins[j][e] = Some(i);
            }
        }
    }
    let next = joins.iter().map(|j| j[1].map(|k| (k, true))).collect();
    Some(QuadricCurves { curves, next, joins, same: false })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::Sphere;
    use crate::geometry::{ParametricCurve, ParametricSurface};

    /// Unsigned distance from `p` to an elementary surface.
    fn distance(s: &Surface, p: Point) -> f64 {
        let (u, v) = s.parameters(p).unwrap();
        s.point(u, v).distance(p)
    }

    fn cylinder(origin: Point, axis: Direction, radius: f64) -> Surface {
        Surface::Cylinder(Cylinder::new(Frame::from_z(origin, axis), radius))
    }

    #[test]
    fn a_plane_across_a_narrow_band_of_a_cone() {
        // A plane along the axis of a steep cone (radius 150, semi-angle -0.17), 23 from the
        // axis, crosses the band 0 <= v <= 3.46 within angles narrower than the samples, and
        // nearly along the generators (where rounding in the linear equation for `v` once
        // made the quadratic formula's far root of it).
        use crate::geometry::surface::Cone;
        let frame = Frame::new(Point::new(0.0, 250.0, 14.4), Direction::new(0.0, 0.0, -1.0), Direction::new(0.73, -0.68, 0.0));
        let cone = Surface::Cone(Cone::new(frame, -0.17, 150.6));
        let plane = Surface::Plane(Plane::new(Point::new(-108.6, 108.6, 153.6), Direction::new(-1.0, 1.0, 0.0)));
        let curves = ruled_quadric(&cone, &plane, 1e-7, Some((0.0, 3.4641))).unwrap();
        assert_eq!(curves.curves.len(), 2);
        for c in &curves.curves {
            let (a, b) = c.domain();
            assert!(b > a && b - a < 0.01, "{a} {b}");
            for i in 0..=4 {
                let p = c.point_clamped(a + (b - a) * i as f64 / 4.0);
                assert!(distance(&cone, p) < 1e-9 && distance(&plane, p) < 1e-9, "{p:?}");
                let (_, v) = cone.parameters(p).unwrap();
                assert!((-1e-6..=3.4641 + 1e-6).contains(&v), "{v}");
            }
        }
    }

    #[test]
    fn plane_meets_cylinder() {
        let cyl = cylinder(Point::new(1.0, 2.0, 0.0), Direction::Z, 2.0);
        let oblique = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 3.0), Direction::new(0.3, -0.2, 1.0)));
        let QuadricIntersection::Ellipses(e) = quadrics(&oblique, &cyl, 1e-7) else { panic!() };
        for i in 0..16 {
            let p = e[0].point(i as f64 * 0.4);
            assert!(distance(&cyl, p) < 1e-12 && distance(&oblique, p) < 1e-12);
        }
        let across = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 3.0), Direction::Z));
        assert!(matches!(quadrics(&across, &cyl, 1e-7), QuadricIntersection::Circles(_)));
        let along = Surface::Plane(Plane::new(Point::new(2.0, 0.0, 0.0), Direction::X));
        let QuadricIntersection::Lines(lines) = quadrics(&along, &cyl, 1e-7) else { panic!() };
        assert_eq!(lines.len(), 2);
        for l in lines {
            assert!(distance(&cyl, l.origin) < 1e-12 && distance(&along, l.origin) < 1e-12);
        }
    }

    #[test]
    fn plane_meets_sphere() {
        let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0));
        let plane = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 1.0), Direction::new(0.0, 1.0, 1.0)));
        let QuadricIntersection::Circles(c) = quadrics(&plane, &sphere, 1e-7) else { panic!() };
        for i in 0..8 {
            let p = c[0].point(i as f64);
            assert!(distance(&sphere, p) < 1e-12 && distance(&plane, p) < 1e-12);
        }
    }

    #[test]
    fn equal_crossing_cylinders_meet_in_ellipses() {
        let a = cylinder(Point::new(0.0, 0.0, -3.0), Direction::Z, 1.0);
        let b = cylinder(Point::new(-3.0, 0.0, 0.0), Direction::new(1.0, 0.0, 0.4), 1.0);
        let QuadricIntersection::Ellipses(es) = quadrics(&a, &b, 1e-7) else { panic!() };
        assert_eq!(es.len(), 2);
        for e in es {
            for i in 0..12 {
                let p = e.point(i as f64 * 0.5);
                assert!(distance(&a, p) < 1e-12 && distance(&b, p) < 1e-12, "{p:?}");
            }
        }
    }

    #[test]
    fn coaxial_cones_meet_in_a_circle_or_not_at_all() {
        use crate::geometry::surface::Cone;
        let cone = |z: f64, angle: f64, radius: f64| Surface::Cone(Cone::new(Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X), angle, radius));
        // The same angle, apexes apart: parallel on each nappe; across nappes they meet where
        // the radii are opposite (h = -radius over tan, here far below).
        match quadrics(&cone(0.0, 0.3, 2.0), &cone(1.0, 0.3, 2.0), 1e-7) {
            QuadricIntersection::Circles(c) => assert_eq!(c.len(), 1),
            other => panic!("{other:?}"),
        }
        // Different angles: one circle on the same nappe (and one across), where
        // 2 + h tan 0.3 = 3 + h tan 0.2.
        let QuadricIntersection::Circles(circles) = quadrics(&cone(0.0, 0.3, 2.0), &cone(0.0, 0.2, 3.0), 1e-7) else { panic!() };
        let h = 1.0 / (0.3f64.tan() - 0.2f64.tan());
        assert!(circles.iter().any(|c| (c.center().z - h).abs() < 1e-9 && (c.radius() - (2.0 + h * 0.3f64.tan())).abs() < 1e-9));
    }

    /// Cones on crossing axes with one apex (OCCT bug10604) meet in generators through it,
    /// each on both cones; cones apart by more than their angles only at the apex.
    #[test]
    fn cones_with_one_apex_meet_in_lines() {
        use crate::geometry::surface::Cone;
        let cone = |axis: Vector, angle: f64| {
            let z = Direction::from_vector(axis);
            let apex_to_origin = 10.0;
            Surface::Cone(Cone::new(Frame::from_z(Point::ORIGIN + Vector::from(z) * apex_to_origin, z), angle, apex_to_origin * angle.tan()))
        };
        let (a, b) = (cone(Vector::new(1.0, 0.0, 0.0), 0.785), cone(Vector::new(0.995, 0.0995, 0.0), 0.733));
        let QuadricIntersection::Lines(lines) = quadrics(&a, &b, 1e-7) else { panic!("{:?}", quadrics(&a, &b, 1e-7)) };
        eprintln!("{lines:?}");
        assert_eq!(lines.len(), 2);
        for line in &lines {
            for t in [-3.0, 1.0, 7.0] {
                let p = line.origin + Vector::from(line.direction) * t;
                assert!(distance(&a, p) < 1e-9 && distance(&b, p) < 1e-9, "{p:?}");
            }
        }
        let narrow = cone(Vector::new(0.0, 1.0, 0.0), 0.1);
        assert!(matches!(quadrics(&cone(Vector::new(1.0, 0.0, 0.0), 0.1), &narrow, 1e-7), QuadricIntersection::Points(_)));
    }

    #[test]
    fn cones_meet_quadrics() {
        use crate::geometry::surface::Cone;
        let cone = Surface::Cone(Cone::new(Frame::WORLD, 0.4, 2.0));
        let across = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 1.5), Direction::Z));
        let QuadricIntersection::Circles(c) = quadrics(&across, &cone, 1e-7) else { panic!() };
        assert!(distance(&cone, c[0].point(1.0)) < 1e-12);
        let Surface::Cone(k) = &cone else { unreachable!() };
        let through_apex = Surface::Plane(Plane::new(k.apex(), Direction::new(1.0, 0.0, 0.2)));
        let QuadricIntersection::Lines(lines) = quadrics(&through_apex, &cone, 1e-7) else { panic!() };
        eprintln!("{lines:?}");
        assert_eq!(lines.len(), 2);
        for l in &lines {
            let p = l.origin + Vector::from(l.direction) * 3.0;
            assert!(distance(&cone, p) < 1e-12 && distance(&through_apex, p) < 1e-12);
        }
        // A plane parallel to a generator: a parabola, kept within the range of v.
        let (sin, cos) = 0.4f64.sin_cos();
        let parallel = Surface::Plane(Plane::new(Point::new(1.0, 0.0, 0.0), Direction::new(cos, 0.0, -sin)));
        assert!(matches!(quadrics(&parallel, &cone, 1e-7), QuadricIntersection::NotAnalytic));
        let curves = ruled_quadric(&cone, &parallel, 1e-7, Some((-5.0, 5.0))).unwrap();
        assert!(!curves.curves.is_empty());
        for c in &curves.curves {
            let (a, b) = c.domain();
            for i in 1..10 {
                let p = c.point(a + (b - a) * i as f64 / 10.0).unwrap();
                assert!(distance(&cone, p) < 1e-9 && distance(&parallel, p) < 1e-9, "{p:?}");
            }
        }
        // A cylinder through the cone.
        let pipe = cylinder(Point::new(-5.0, 0.0, 1.0), Direction::X, 0.5);
        let curves = ruled_quadric(&cone, &pipe, 1e-7, Some((-5.0, 5.0))).unwrap();
        for c in &curves.curves {
            let (a, b) = c.domain();
            for i in 1..10 {
                let p = c.point(a + (b - a) * i as f64 / 10.0).unwrap();
                assert!(distance(&cone, p) < 1e-9 && distance(&pipe, p) < 1e-9, "{p:?}");
            }
        }
    }

    #[test]
    fn crossing_cylinders() {
        // A thin cylinder through a thick one: two closed branches all the way round.
        let thick = Cylinder::new(Frame::WORLD, 2.0);
        let thin = cylinder(Point::new(-5.0, 0.3, 0.5), Direction::X, 1.0);
        let curves = cylinder_quadric(&Cylinder::new(Frame::from_z(Point::new(-5.0, 0.3, 0.5), Direction::X), 1.0), &Surface::Cylinder(thick), 1e-7).unwrap();
        assert_eq!(curves.curves.len(), 2);
        assert!(curves.next.iter().all(|n| n.is_none()));
        for c in &curves.curves {
            for i in 0..36 {
                let p = c.point(i as f64 * 0.17).unwrap();
                assert!(distance(&Surface::Cylinder(thick), p) < 1e-9 && distance(&thin, p) < 1e-9);
            }
        }
        // Offset so that the thin cylinder only bites: one loop of two linked branches.
        let biting = Cylinder::new(Frame::from_z(Point::new(-5.0, 1.8, 0.0), Direction::X), 1.0);
        let curves = cylinder_quadric(&biting, &Surface::Cylinder(thick), 1e-7).unwrap();
        assert_eq!(curves.curves.len(), 2);
        assert_eq!(curves.next, vec![Some((1, true)), Some((0, true))]);
        let (a, b) = curves.curves[0].domain();
        assert!(a < b && b - a < TAU);
        let (lo, hi) = (curves.curves[0].point(a).unwrap(), curves.curves[1].point(a).unwrap());
        assert!(lo.distance(hi) < 1e-4, "branches meet at the domain's ends");
    }
}
