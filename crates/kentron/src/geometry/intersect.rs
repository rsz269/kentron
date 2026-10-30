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
        (Surface::Torus(torus), other) | (other, Surface::Torus(torus)) => torus_elementary(torus, other, tol),
        _ => Q::NotAnalytic,
    }
}

/// A torus meets a plane across its axis, or a cylinder, cone, sphere or torus about the
/// same axis, in circles (where its meridian circle meets the other's meridian, in the
/// plane of distance from the axis and height along it), and a plane through its axis in
/// two circles. Otherwise not in closed form.
fn torus_elementary(torus: &crate::geometry::surface::Torus, other: &Surface, tol: f64) -> QuadricIntersection {
    use QuadricIntersection as Q;
    let f = torus.frame();
    let (c, z) = (f.origin(), Vector::from(f.z()));
    let (big, small) = (torus.major_radius(), torus.minor_radius());
    let parallel = |u: Vector| u.cross(z).length() <= ANGULAR.sqrt();
    let height = |p: Point| z.dot(p - c);
    let on_axis = |p: Point| {
        let d = p - c;
        (d - z * d.dot(z)).length() <= tol
    };
    // Meridian lines `a ρ + b h = k` and circles (centre, radius) in the (ρ, h) plane.
    enum Meridian {
        Line(f64, f64, f64),
        Circle((f64, f64), f64),
    }
    let meridian = match other {
        Surface::Plane(plane) => {
            let n = Vector::from(plane.normal());
            if parallel(n) {
                Meridian::Line(0.0, 1.0, height(plane.origin()))
            } else if n.dot(z).abs() <= ANGULAR.sqrt() && n.dot(c - plane.origin()).abs() <= tol {
                // Through the axis: the meridian circles on either side.
                let across = Vector::from(plane.normal()).cross(z).normalized();
                return Q::Circles(
                    [1.0, -1.0]
                        .iter()
                        .map(|&side| Circle::new(Frame::new(c + across * (side * big), plane.normal(), Direction::from_vector(across * side)), small))
                        .collect(),
                );
            } else {
                return Q::NotAnalytic;
            }
        }
        Surface::Cylinder(cyl) if parallel(Vector::from(cyl.frame().z())) && on_axis(cyl.frame().origin()) => {
            Meridian::Line(1.0, 0.0, cyl.radius())
        }
        Surface::Cone(cone) if parallel(Vector::from(cone.frame().z())) && on_axis(cone.frame().origin()) => {
            // ρ = R + s (h - h0) tan α.
            let (sin, cos) = cone.semi_angle_sin_cos();
            let s = Vector::from(cone.frame().z()).dot(z).signum();
            let slope = s * sin / cos;
            Meridian::Line(1.0, -slope, cone.radius() - slope * height(cone.frame().origin()))
        }
        Surface::Sphere(sphere) if on_axis(sphere.frame().origin()) => Meridian::Circle((0.0, height(sphere.frame().origin())), sphere.radius()),
        Surface::Torus(t) if parallel(Vector::from(t.frame().z())) && on_axis(t.frame().origin()) => {
            Meridian::Circle((t.major_radius(), height(t.frame().origin())), t.minor_radius())
        }
        _ => return Q::NotAnalytic,
    };
    // Where the torus's meridian circle (centre (R, 0), radius r) meets it.
    let mut found: Vec<(f64, f64)> = Vec::new();
    let line = |a: f64, b: f64, k: f64, found: &mut Vec<(f64, f64)>| {
        let len = (a * a + b * b).sqrt();
        let (a, b, k) = (a / len, b / len, k / len);
        let d = a * big - k;
        let foot = (big - a * d, -b * d);
        if (d.abs() - small).abs() <= tol {
            found.push(foot);
        } else if d.abs() < small {
            let half = (small * small - d * d).sqrt();
            found.push((foot.0 - b * half, foot.1 + a * half));
            found.push((foot.0 + b * half, foot.1 - a * half));
        }
    };
    match meridian {
        Meridian::Line(a, b, k) => line(a, b, k, &mut found),
        Meridian::Circle((p, q), radius) => {
            let gap = ((p - big).powi(2) + q * q).sqrt();
            if gap <= tol {
                return if (radius - small).abs() <= tol { Q::Same } else { Q::Empty };
            }
            // The radical line: 2(p - R) ρ + 2q h = p² + q² - radius² - R² + r².
            line(2.0 * (p - big), 2.0 * q, p * p + q * q - radius * radius - big * big + small * small, &mut found);
        }
    }
    let circles: Vec<Circle> = found
        .into_iter()
        .filter(|&(rho, _)| rho > tol)
        .map(|(rho, h)| Circle::new(Frame::new(c + z * h, f.z(), f.x()), rho))
        .collect();
    if circles.is_empty() { Q::Empty } else { Q::Circles(circles) }
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
    // Otherwise a conic: with the plane's points `q + s e1 + t e2` (q the apex's foot, e1 the
    // axis's direction in the plane, e2 across it), the double cone `(w·z)² = cos²α |w|²`
    // (w from the apex) reads `A s² + B s + C - cos²α t² = 0`: an ellipse (A < 0), a parabola
    // (A = 0) or a hyperbola of two branches (A > 0).
    let d = n.dot(apex - plane.origin());
    let foot = apex - n * d;
    let e1 = (z - n * n.dot(z)).normalized();
    let e2 = n.cross(e1);
    let (nz, ez) = (n.dot(z), e1.dot(z));
    let c2 = cos * cos;
    let (a2, b, c0) = (ez * ez - c2, -2.0 * d * nz * ez, d * d * (nz * nz - c2));
    let plane_z = Direction::from_vector(n);
    let at = |s: f64| foot + e1 * s;
    if a2.abs() <= 1e-12 {
        // A parabola, `s - s_v = (cos²α / B) t²`, opening the way of B along e1.
        if b.abs() <= 1e-300 {
            return Q::NotAnalytic;
        }
        let vertex = at(-c0 / b);
        let opening = if b > 0.0 { e1 } else { -e1 };
        return Q::Parabolas(vec![Parabola::new(Frame::new(vertex, plane_z, Direction::from_vector(opening)), b.abs() / (4.0 * c2))]);
    }
    let s0 = -b / (2.0 * a2);
    let k = a2 * s0 * s0 - c0;
    let centre = at(s0);
    if k.abs() <= tol * tol {
        return Q::NotAnalytic;
    }
    if a2 < 0.0 {
        // `(s - s0)² / (k / A) + t² / (-k / cos²α) = 1`.
        if k >= 0.0 {
            return Q::Empty;
        }
        let (rs, rt) = ((k / a2).sqrt(), (-k / c2).sqrt());
        let (x, major, minor) = if rs >= rt { (e1, rs, rt) } else { (e2, rt, rs) };
        return Q::Ellipses(vec![Ellipse::new(Frame::new(centre, plane_z, Direction::from_vector(x)), major, minor)]);
    }
    // A hyperbola across s (k > 0) or across t (k < 0); a branch either way along it.
    let (x, major, minor) = if k > 0.0 { (e1, (k / a2).sqrt(), (k / c2).sqrt()) } else { (e2, (-k / c2).sqrt(), (-k / a2).sqrt()) };
    Q::Hyperbolas(vec![
        Hyperbola::new(Frame::new(centre, plane_z, Direction::from_vector(x)), major, minor),
        Hyperbola::new(Frame::new(centre, plane_z, Direction::from_vector(-x)), major, minor),
    ])
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
    /// The other surface itself (its quadric's value found from its own frame where the
    /// world coefficients' rounding matters).
    other: Surface,
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

/// The value of the quadric of `surface` (as [`quadric_of`] gives it) at `p`, from the
/// surface's own frame: without the cancelling of large terms that the world coefficients
/// have far from the origin. `None` for other surfaces.
fn local_value(surface: &Surface, p: Point) -> Option<f64> {
    // (A difference of squares as a product, both factors exact enough.)
    let axial = |origin: Point, z: Vector, rho_on: &dyn Fn(f64) -> f64| {
        let d = p - origin;
        let h = d.dot(z);
        let rho = (d - z * h).length();
        let on = rho_on(h);
        (rho - on) * (rho + on)
    };
    Some(match surface {
        Surface::Plane(plane) => Vector::from(plane.normal()).dot(p - plane.origin()),
        Surface::Sphere(s) => {
            let d = (p - s.frame().origin()).length();
            (d - s.radius()) * (d + s.radius())
        }
        Surface::Cylinder(c) => axial(c.frame().origin(), Vector::from(c.frame().z()), &|_| c.radius()),
        Surface::Cone(c) => {
            let (sin, cos) = c.semi_angle_sin_cos();
            axial(c.apex(), Vector::from(c.frame().z()), &|h| h.abs() * sin / cos)
        }
        _ => return None,
    })
}

/// The size of the terms [`quadric_value`] sums at `p`, which its rounding goes with.
fn quadric_size(q: &[f64; 10], p: Point) -> f64 {
    let (x, y, z) = (p.x, p.y, p.z);
    (q[0] * x * x).abs() + (q[1] * y * y).abs() + (q[2] * z * z).abs() + 2.0 * ((q[3] * x * y).abs() + (q[4] * x * z).abs() + (q[5] * y * z).abs())
        + 2.0 * ((q[6] * x).abs() + (q[7] * y).abs() + (q[8] * z).abs())
        + q[9].abs()
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
        let points = self.base.points_along_v(angle, [-1.0, 0.0, 1.0]);
        let [minus, zero, plus] = points.map(|p| quadric_value(&self.quadric, p));
        // A plane is met along a generator linearly: the rounding left in `a` is noise (where
        // the generator runs nearly along the plane it can outweigh `b`).
        let planar = self.quadric[..6].iter().all(|&c| c == 0.0);
        // (So is an `a` within the rounding of the quadric's terms: a generator parallel to
        // one of the other quadric's, as of two cones with opposite half-angles, whose terms
        // far from the origin cancel to small values.)
        let a = 0.5 * (plus + minus) - zero;
        let rounding = 8.0 * f64::EPSILON * points.map(|p| quadric_size(&self.quadric, p)).into_iter().fold(0.0, f64::max);
        (if planar || a.abs() <= rounding { 0.0 } else { a }, 0.5 * (plus - minus), zero)
    }

    /// The branch's generator parameter at `angle`; `None` where it does not exist.
    fn height(&self, angle: f64) -> Option<f64> {
        let (a, b, c) = self.height_equation(angle);
        let scale = a.abs().max(b.abs()).max(c.abs()).max(1e-300);
        if a.abs() <= 1e-12 * scale {
            return (b.abs() > 1e-12 * scale).then(|| -c / b);
        }
        let disc = b * b - 4.0 * a * c;
        // (Rounding in the discriminant goes with the size of its terms.)
        let within = disc >= -1e-12 * (b * b).max(4.0 * (a * c).abs());
        if let Some(square) = self.near_double_root(angle, a, b, disc) {
            // (Off the quadric by more than the tolerance's square, and by the discriminant
            // too: no root.)
            return (within || square >= -1e-14).then(|| -b / (2.0 * a) + if self.upper { 1.0 } else { -1.0 } * square.max(0.0).sqrt());
        }
        within.then(|| self.root(a, b, disc))
    }

    /// Near a double root of `a v² + b v + c` (where the surfaces touch, or the branches
    /// meet) the discriminant `disc` is mostly rounding: the square of the roots' distance
    /// from the parabola's vertex, from the generator's depth into the quadric there, found
    /// from the other surface's own frame. `None` elsewhere.
    fn near_double_root(&self, angle: f64, a: f64, b: f64, disc: f64) -> Option<f64> {
        // (Small against `b`, or within the rounding of `c`'s terms where `b` is small too.)
        let vertex = self.base.point(angle, -b / (2.0 * a));
        let noise = 4.0 * a.abs() * 1e3 * f64::EPSILON * quadric_size(&self.quadric, vertex);
        if disc.abs() > (1e-4 * b * b).max(noise) {
            return None;
        }
        let depth = local_value(&self.other, vertex)?;
        Some(-depth / a)
    }

    /// The branch's root of `a v² + b v + c`, of discriminant `disc`.
    fn root(&self, a: f64, b: f64, disc: f64) -> f64 {
        let sign = if self.upper == (a > 0.0) { 1.0 } else { -1.0 };
        (-b + sign * disc.max(0.0).sqrt()) / (2.0 * a)
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
            let disc = b * b - 4.0 * a * c;
            match self.near_double_root(angle, a, b, disc) {
                Some(square) => -b / (2.0 * a) + if self.upper { 1.0 } else { -1.0 } * square.max(0.0).sqrt(),
                None => self.root(a, b, disc),
            }
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
    let branch = |first: f64, last: f64, upper: bool| QuadricCurve { base: base.clone(), quadric, other: other.clone(), first, last, upper };
    let probe = branch(0.0, TAU, true);
    let linear = (0..8).all(|i| {
        let (a, b, c) = probe.height_equation(i as f64);
        a.abs() <= 1e-12 * a.abs().max(b.abs()).max(c.abs()).max(1e-300)
    });
    let (v_lo, v_hi) = v_range.unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
    // Where a branch exists, and why it ends: the discriminant vanishing (the branches
    // join) or `v` leaving its range.
    let whole = [branch(0.0, TAU, false), branch(0.0, TAU, true)];
    let exists = |upper: bool, angle: f64| whole[upper as usize].height(angle).filter(|v| *v >= v_lo && *v <= v_hi).is_some();
    // (Or, where `v` is still well within its range there, the branch ended for its roots
    // meeting: the discriminant only nears zero as the square root of the angle's distance
    // from the join, too slowly to vanish against the terms of `c` alone where `b` is zero.)
    let joins_here = |upper: bool, angle: f64| {
        let (a, b, c) = probe.height_equation(angle);
        let within = |v: f64| {
            let margin = 1e-6 * (1.0 + v.abs()) + if v_hi > v_lo && v_hi.is_finite() && v_lo.is_finite() { 1e-6 * (v_hi - v_lo) } else { 0.0 };
            v > v_lo + margin && v < v_hi - margin
        };
        !linear && (b * b - 4.0 * a * c <= 1e-6 * (b * b).max(4.0 * (a * c).abs()) || whole[upper as usize].height(angle).is_some_and(within))
    };
    const SAMPLES: usize = 1440;
    let angles: Vec<f64> = (0..SAMPLES).map(|i| TAU * i as f64 / SAMPLES as f64).collect();
    let refine = |upper: bool, inside: f64, outside: f64| {
        let (mut a, mut b) = (inside, outside);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            // (Done once the halves are neighbouring numbers.)
            if mid == a || mid == b {
                break;
            }
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
        let height = |t: f64| whole[upper as usize].height(t);
        let (h0, h1) = (height(angle)?, height(angle + step)?);
        if !((h0 < v_lo && h1 > v_hi) || (h0 > v_hi && h1 < v_lo)) {
            return None;
        }
        let root = |level: f64| -> Option<f64> {
            let (mut a, mut b) = (angle, angle + step);
            let below = |t: f64| height(t).map(|h| (h < level) == (h0 < level));
            for _ in 0..100 {
                let mid = 0.5 * (a + b);
                if mid == a || mid == b {
                    break;
                }
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
    let step = TAU / SAMPLES as f64;
    // Stretches of samples where no branch can start or end, nor run across the range of
    // `v`: there the first sample stands for the rest.
    let (same_as, quiet_after) = quiet_stretches(&probe, linear, (v_lo, v_hi), SAMPLES);
    // Each branch's domains: (start, end, start joins, end joins).
    let mut domains: Vec<(bool, f64, f64, bool, bool)> = Vec::new();
    for upper in if linear { vec![true] } else { vec![true, false] } {
        let mut valid: Vec<bool> = Vec::with_capacity(SAMPLES);
        for (k, &t) in angles.iter().enumerate() {
            let v = if same_as[k] == k { exists(upper, t) } else { valid[same_as[k]] };
            valid.push(v);
        }
        if valid.iter().all(|&v| v) {
            domains.push((upper, 0.0, TAU, false, false));
            continue;
        }
        // Walk round from an invalid sample (where none is valid, for domains between them).
        let start = valid.iter().position(|&v| !v).expect("an invalid sample");
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
                domains.push((upper, a, b, joins_here(upper, a), joins_here(upper, b)));
                i = j;
            } else if let Some(domain) = (!quiet_after[k]).then(|| crossing(upper, angles[start] + step * i as f64, step)).flatten() {
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

/// A trigonometric polynomial `c₀ + Σ cₖ cos kθ + sₖ sin kθ` (k up to 7), fitted to values
/// at 16 equally spaced angles.
struct Trigonometric {
    cos: [f64; 8],
    sin: [f64; 8],
}

impl Trigonometric {
    const ANGLES: usize = 16;

    fn fit(values: &[f64]) -> Self {
        let n = Self::ANGLES as f64;
        let (mut cos, mut sin) = ([0.0; 8], [0.0; 8]);
        for (j, f) in values.iter().enumerate() {
            let t = TAU * j as f64 / n;
            for k in 0..8 {
                let (s, c) = (k as f64 * t).sin_cos();
                cos[k] += f * c * if k == 0 { 1.0 / n } else { 2.0 / n };
                sin[k] += f * s * 2.0 / n;
            }
        }
        Self { cos, sin }
    }

    /// The value and the derivative at `t`.
    fn at(&self, t: f64) -> (f64, f64) {
        let (s1, c1) = t.sin_cos();
        let (mut s, mut c) = (0.0, 1.0);
        let (mut f, mut d) = (self.cos[0], 0.0);
        for k in 1..8 {
            (s, c) = (s * c1 + c * s1, c * c1 - s * s1);
            f += self.cos[k] * c + self.sin[k] * s;
            d += k as f64 * (self.sin[k] * c - self.cos[k] * s);
        }
        (f, d)
    }

    /// The sum of the coefficients' sizes, bounding the value.
    fn size(&self) -> f64 {
        self.cos.iter().chain(&self.sin).map(|c| c.abs()).sum()
    }

    /// True if the polynomial stays further than `slack` from zero over `[t, t + span]`
    /// (by its second derivative's bound).
    fn clear(&self, t: f64, span: f64, slack: f64) -> bool {
        let (f, d) = self.at(t);
        let curvature: f64 = (1..8).map(|k| (k * k) as f64 * (self.cos[k].abs() + self.sin[k].abs())).sum();
        f.abs() > d.abs() * span + 0.5 * curvature * span * span + slack
    }
}

/// For `samples` angles round the base of a ruled quadric's section (`probe`, either
/// branch), the sample each one surely agrees with on whether a branch exists there, and
/// whether a branch surely neither starts, ends nor runs across `range` (of the generator
/// parameter `v`) from each sample to the next. A branch starts or ends only where the
/// discriminant of `a v² + b v + c`, `a` (`b`, for a linear equation) or the quadric's value at
/// either end of the range changes sign, all trigonometric polynomials of the angle of degree
/// at most 4 (2 for a cylinder): fitted, they show stretches clear of zero. (Where the fit
/// fails, from a surface that is not one, every sample stands alone.)
fn quiet_stretches(probe: &QuadricCurve, linear: bool, (v_lo, v_hi): (f64, f64), samples: usize) -> (Vec<usize>, Vec<bool>) {
    let mut same_as: Vec<usize> = (0..samples).collect();
    let mut quiet_after = vec![false; samples];
    let n = Trigonometric::ANGLES;
    // The functions, and the slack each is clear of zero by: generous against the rounding
    // the samples allow for (and the rounding of a near double root).
    let functions = |t: f64| {
        let (a, b, c) = probe.height_equation(t);
        let level = |v: f64| if v.is_finite() { a * v * v + b * v + c } else { 0.0 };
        [a, b, c, b * b - 4.0 * a * c, level(v_lo), level(v_hi)]
    };
    let values: Vec<[f64; 6]> = (0..n).map(|j| functions(TAU * j as f64 / n as f64)).collect();
    let fits: Vec<Trigonometric> = (0..6).map(|i| Trigonometric::fit(&values.iter().map(|v| v[i]).collect::<Vec<_>>())).collect();
    let [a, b, c, disc, lo, hi] = [0, 1, 2, 3, 4, 5].map(|i| &fits[i]);
    let (a_size, b_size, c_size) = (a.size(), b.size(), c.size());
    let coefficients = a_size + b_size + c_size;
    // (Checked halfway between the fitted angles.)
    let fitted = (0..n).all(|j| {
        let t = TAU * (j as f64 + 0.5) / n as f64;
        functions(t).iter().zip(&fits).all(|(f, fit)| (f - fit.at(t).0).abs() <= 1e-9 * (fit.size() + coefficients))
    });
    if !fitted || !(coefficients > 0.0 && coefficients.is_finite()) {
        return (same_as, quiet_after);
    }
    let noise = values
        .iter()
        .enumerate()
        .filter(|(_, v)| v[0] != 0.0)
        .map(|(j, v)| {
            let vertex = probe.base.point(TAU * j as f64 / n as f64, -v[1] / (2.0 * v[0]));
            4.0 * v[0].abs() * 1e3 * f64::EPSILON * quadric_size(&probe.quadric, vertex)
        })
        .fold(0.0, f64::max);
    let mut guards: Vec<(&Trigonometric, f64)> = Vec::new();
    if linear {
        guards.push((b, 1e-9 * coefficients));
    } else {
        guards.push((a, 1e-9 * coefficients));
        guards.push((disc, 1e-3 * b_size * b_size + 100.0 * noise + 1e-9 * (disc.size() + 4.0 * a_size * c_size)));
    }
    for (fit, v) in [(lo, v_lo), (hi, v_hi)] {
        if v.is_finite() {
            guards.push((fit, 1e-9 * (fit.size() + a_size * v * v + b_size * v.abs() + c_size)));
        }
    }
    let step = TAU / samples as f64;
    let quiet = |k: usize, len: usize| guards.iter().all(|(g, slack)| g.clear(k as f64 * step, len as f64 * step, *slack));
    let mut stack: Vec<(usize, usize)> = (0..samples).step_by(32).map(|k| (k, 32.min(samples - k))).collect();
    while let Some((k, len)) = stack.pop() {
        if quiet(k, len) {
            for i in k..k + len {
                same_as[i] = k;
                quiet_after[i] = true;
            }
        } else if len > 1 {
            stack.push((k, len / 2));
            stack.push((k + len / 2, len - len / 2));
        }
    }
    (same_as, quiet_after)
}

/// Where an extrusion meets a surface unchanged along the extrusion's direction (a plane
/// or cylinder along it, or another extrusion the same way): lines along it, through the
/// points where the extrusion's basis curve meets the surface. `None` for another surface
/// (or where the basis curve runs along the surface, the surfaces overlapping).
pub fn extrusion_lines(extrusion: &super::surface::ExtrusionSurface, other: &Surface, tolerance: f64) -> Option<Vec<Line>> {
    let d = extrusion.direction();
    let along = |u: Vector| u.cross(Vector::from(d)).length() <= ANGULAR.sqrt();
    let invariant = match other {
        Surface::Plane(p) => Vector::from(p.normal()).dot(Vector::from(d)).abs() <= ANGULAR.sqrt(),
        Surface::Cylinder(c) => along(Vector::from(c.frame().z())),
        Surface::Extrusion(e) => along(Vector::from(e.direction())),
        _ => false,
    };
    if !invariant {
        return None;
    }
    let found = curve_surface(extrusion.basis(), other, tolerance);
    if !found.segments.is_empty() {
        return None;
    }
    Some(found.points.iter().map(|p| Line::new(p.point, d)).collect())
}

/// A point where a curve meets a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfacePoint {
    pub point: Point,
    /// Parameter on the curve.
    pub t: f64,
    /// Parameters on the surface.
    pub u: f64,
    pub v: f64,
}

/// Where a curve meets a surface: isolated points, and the stretches of the curve (as
/// parameter ranges) lying on the surface.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CurveSurfaceIntersection {
    pub points: Vec<CurveSurfacePoint>,
    pub segments: Vec<(f64, f64)>,
}

/// The points where a curve meets a surface (within its bounds), within `tolerance`, sorted
/// by parameter on the curve; and the stretches of the curve on the surface. A line meets a
/// plane, cylinder or sphere in closed form; otherwise the curve is followed in fine steps
/// (an unbounded one over the extent of a bounded surface, else near the surface's
/// parametric centre).
#[doc(alias = "GeomAPI_IntCS")]
pub fn curve_surface(curve: &super::curve::Curve, surface: &Surface, tolerance: f64) -> CurveSurfaceIntersection {
    use super::curve::Curve;
    let tol = tolerance.max(1e-9);
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    let on_surface = |t: f64, start: Option<(f64, f64)>| {
        let p = curve.point(t);
        let (u, v, d) = extrema::nearest_on_surface(p, surface, start);
        (CurveSurfacePoint { point: p, t, u, v }, d)
    };
    let mut out = CurveSurfaceIntersection::default();
    if let Curve::Line(line) = curve {
        match line_elementary(line, surface, tol) {
            Some(None) => {
                out.segments.push((first, last));
                return out;
            }
            Some(Some(ts)) => {
                let slack = |t: f64| 1e-12 * (1.0 + t.abs());
                for t in ts.into_iter().filter(|&t| t >= first - slack(t) && t <= last + slack(t)) {
                    out.points.push(on_surface(t.clamp(first, last), None).0);
                }
                return out;
            }
            None => {}
        }
    }
    let (lo, hi) = search_range(curve, surface);
    if !(hi > lo) {
        return out;
    }
    let n = match curve {
        Curve::BSpline(b) => (32 * b.knots().len()).clamp(512, 8192),
        _ => 512,
    };
    let mut samples: Vec<(CurveSurfacePoint, f64)> = Vec::with_capacity(n + 1);
    let mut start = None;
    for i in 0..=n {
        let t = if i == n { hi } else { lo + (hi - lo) * i as f64 / n as f64 };
        let s = on_surface(t, start);
        start = Some((s.0.u, s.0.v));
        samples.push(s);
    }
    // Runs of samples on the surface are stretches of the curve on it.
    let mut i = 0;
    let mut in_segment = vec![false; samples.len()];
    while i < samples.len() {
        let mut j = i;
        while j < samples.len() && samples[j].1 <= tol {
            j += 1;
        }
        if j >= i + 3 {
            out.segments.push((samples[i].0.t, samples[j - 1].0.t));
            in_segment[i..j].iter_mut().for_each(|x| *x = true);
        }
        i = j.max(i + 1);
    }
    // Each local minimum of the distance: found by golden sections, then made exact by
    // Newton steps on curve(t) = surface(u, v).
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for i in 0..samples.len() {
        let d = samples[i].1;
        let left = if i > 0 { samples[i - 1].1 } else { f64::INFINITY };
        let right = if i + 1 < samples.len() { samples[i + 1].1 } else { f64::INFINITY };
        // (Within a step of a crossing the distance is at most the step's length.)
        let chord = |j: usize, k: usize| samples[j].0.point.distance(samples[k].0.point);
        let reach = chord(i.saturating_sub(1), i).max(chord(i, (i + 1).min(samples.len() - 1)));
        if in_segment[i] || !(d <= left && d <= right) || d > reach + 1e3 * tol {
            continue;
        }
        let (mut a, mut b) = (samples[i.saturating_sub(1)].0.t, samples[(i + 1).min(samples.len() - 1)].0.t);
        let mut best = samples[i];
        for _ in 0..80 {
            let (x1, x2) = (b - g * (b - a), a + g * (b - a));
            let guess = Some((best.0.u, best.0.v));
            let (f1, f2) = (on_surface(x1, guess), on_surface(x2, guess));
            if f1.1 < f2.1 {
                b = x2;
                if f1.1 < best.1 {
                    best = f1;
                }
            } else {
                a = x1;
                if f2.1 < best.1 {
                    best = f2;
                }
            }
        }
        let found = polish_curve_surface(curve, surface, best.0, lo, hi);
        if found.point.distance(surface.point(found.u, found.v)) > tol {
            continue;
        }
        if out.points.iter().any(|q| q.point.distance(found.point) <= tol || (q.t - found.t).abs() <= 1e-9 * (1.0 + found.t.abs())) {
            continue;
        }
        out.points.push(found);
    }
    out.points.sort_by(|a, b| a.t.total_cmp(&b.t));
    out
}

/// Where a line meets a plane, cylinder or sphere: `Some(None)` if it lies on it,
/// `Some(Some(parameters))` else; `None` for another surface.
fn line_elementary(line: &Line, surface: &Surface, tol: f64) -> Option<Option<Vec<f64>>> {
    let (o, d) = (line.origin, Vector::from(line.direction));
    // Where a line at `h` from a centre, nearest it at parameter `t0`, meets a circle of
    // radius `r` about that centre (`speed` the line's parametric speed across it).
    let circle = |t0: f64, h: f64, r: f64, speed: f64| {
        if (h - r).abs() <= tol {
            vec![t0]
        } else if h < r {
            let half = (r * r - h * h).sqrt() / speed;
            vec![t0 - half, t0 + half]
        } else {
            Vec::new()
        }
    };
    Some(match surface {
        Surface::Plane(plane) => {
            let n = Vector::from(plane.normal());
            let (along, off) = (n.dot(d), n.dot(o - plane.origin()));
            if along.abs() <= 1e-12 {
                return Some(if off.abs() <= tol { None } else { Some(Vec::new()) });
            }
            Some(vec![-off / along])
        }
        Surface::Sphere(sphere) => {
            let w = o - sphere.frame().origin();
            let t0 = -w.dot(d);
            Some(circle(t0, (w + d * t0).length(), sphere.radius(), 1.0))
        }
        Surface::Cylinder(cylinder) => {
            let f = cylinder.frame();
            let z = Vector::from(f.z());
            let across = |v: Vector| v - z * v.dot(z);
            let (w, dd) = (across(o - f.origin()), across(d));
            let a = dd.length_squared();
            if a <= 1e-24 {
                return Some(if (w.length() - cylinder.radius()).abs() <= tol { None } else { Some(Vec::new()) });
            }
            let t0 = -w.dot(dd) / a;
            Some(circle(t0, (w + dd * t0).length(), cylinder.radius(), a.sqrt()))
        }
        _ => return None,
    })
}

/// The parameter range of `curve` to search for `surface`: its own where finite; along an
/// unbounded line, across the extent of a bounded surface (or around its parametric
/// centre).
fn search_range(curve: &super::curve::Curve, surface: &Surface) -> (f64, f64) {
    use super::curve::Curve;
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    let finite = |x: f64| x.abs() < 1e100;
    if finite(first) && finite(last) {
        return (first, last);
    }
    let Curve::Line(line) = curve else { return (first.max(-1e3), last.min(1e3)) };
    let (u0, u1, v0, v1) = surface.bounds();
    let along = |p: Point| (p - line.origin).dot(Vector::from(line.direction));
    let (lo, hi) = if [u0, u1, v0, v1].into_iter().all(finite) {
        let n = 16;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in 0..=n {
            for j in 0..=n {
                let t = along(surface.point(u0 + (u1 - u0) * i as f64 / n as f64, v0 + (v1 - v0) * j as f64 / n as f64));
                (lo, hi) = (lo.min(t), hi.max(t));
            }
        }
        let margin = 0.05 * (hi - lo) + 1.0;
        (lo - margin, hi + margin)
    } else {
        let mid = |a: f64, b: f64| if finite(a) && finite(b) { 0.5 * (a + b) } else if finite(a) { a } else if finite(b) { b } else { 0.0 };
        let t = along(surface.point(mid(u0, u1), mid(v0, v1)));
        (t - 1e4, t + 1e4)
    };
    (lo.max(first), hi.min(last))
}

/// Newton steps on curve(t) = surface(u, v) from `start` (keeping within the ranges).
fn polish_curve_surface(curve: &super::curve::Curve, surface: &Surface, start: CurveSurfacePoint, lo: f64, hi: f64) -> CurveSurfacePoint {
    let (u0, u1, v0, v1) = surface.bounds();
    let clamp_u = |u: f64| if surface.is_u_periodic() { u } else { u.clamp(u0, u1) };
    let clamp_v = |v: f64| if surface.is_v_periodic() { v } else { v.clamp(v0, v1) };
    let (mut t, mut u, mut v) = (start.t, start.u, start.v);
    let gap = |t: f64, u: f64, v: f64| curve.point(t).distance(surface.point(u, v));
    let mut best = (gap(t, u, v), t, u, v);
    for _ in 0..30 {
        let (c, dc) = curve.d1(t);
        let (p, du, dv) = surface.d1(u, v);
        let r = p - c;
        // dc·dt − du·du_ − dv·dv_ = r
        let (a, b, e) = (dc, -du, -dv);
        let det = a.dot(b.cross(e));
        if det.abs() <= 1e-300 {
            break;
        }
        let (x, y, z) = (r.dot(b.cross(e)) / det, a.dot(r.cross(e)) / det, a.dot(b.cross(r)) / det);
        (t, u, v) = ((t + x).clamp(lo, hi), clamp_u(u + y), clamp_v(v + z));
        let g = gap(t, u, v);
        if g < best.0 {
            best = (g, t, u, v);
        } else if g >= best.0 {
            break;
        }
    }
    let (_, t, u, v) = best;
    if best.0 < gap(start.t, start.u, start.v) {
        CurveSurfacePoint { point: curve.point(t), t, u, v }
    } else {
        start
    }
}

/// The point where two surfaces meet nearest `near`, in the plane through `near` across
/// their curve of intersection there: Newton steps on `s1(u1, v1) = s2(u2, v2)` from the
/// points nearest `near` (or the parameters `start`), with the parameters `[u1, v1, u2,
/// v2]`. `None` where the steps do not converge (the surfaces do not meet near there, or
/// meet tangentially).
pub fn surfaces_point(s1: &Surface, s2: &Surface, near: Point, start: Option<[f64; 4]>) -> Option<(Point, [f64; 4])> {
    let mut x = match start {
        Some(x) => x,
        None => {
            let (u1, v1, _) = extrema::nearest_on_surface(near, s1, None);
            let (u2, v2, _) = extrema::nearest_on_surface(near, s2, None);
            [u1, v1, u2, v2]
        }
    };
    let clamp = |s: &Surface, u: f64, v: f64| {
        let (u0, u1, v0, v1) = s.bounds();
        (if s.is_u_periodic() { u } else { u.clamp(u0, u1) }, if s.is_v_periodic() { v } else { v.clamp(v0, v1) })
    };
    // The direction of the curve, fixed at the start.
    let (_, a1, b1) = s1.d1(x[0], x[1]);
    let (_, a2, b2) = s2.d1(x[2], x[3]);
    let along = a1.cross(b1).cross(a2.cross(b2));
    if along.length() <= 1e-300 {
        return None;
    }
    let along = along.normalized();
    let scale = 1.0 + (near - Point::ORIGIN).length();
    for _ in 0..40 {
        let (p1, du1, dv1) = s1.d1(x[0], x[1]);
        let (p2, du2, dv2) = s2.d1(x[2], x[3]);
        let gap = p1 - p2;
        let f = [gap.x, gap.y, gap.z, (p1 - near).dot(along)];
        if f[..3].iter().all(|c| c.abs() <= 1e-14 * scale) && f[3].abs() <= 1e-12 * scale {
            break;
        }
        let m = [
            [du1.x, dv1.x, -du2.x, -dv2.x],
            [du1.y, dv1.y, -du2.y, -dv2.y],
            [du1.z, dv1.z, -du2.z, -dv2.z],
            [du1.dot(along), dv1.dot(along), 0.0, 0.0],
        ];
        let step = solve4(m, [-f[0], -f[1], -f[2], -f[3]])?;
        let (u1, v1) = clamp(s1, x[0] + step[0], x[1] + step[1]);
        let (u2, v2) = clamp(s2, x[2] + step[2], x[3] + step[3]);
        x = [u1, v1, u2, v2];
    }
    let (p1, p2) = (s1.point(x[0], x[1]), s2.point(x[2], x[3]));
    (p1.distance(p2) <= 1e-9 * scale).then(|| (p1 + (p2 - p1) * 0.5, x))
}

/// The solution of a 4 × 4 linear system (Gaussian elimination with partial pivoting);
/// `None` if singular.
fn solve4(mut m: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    let size = m.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs()));
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))?;
        if m[pivot][col].abs() <= 1e-14 * size {
            return None;
        }
        m.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..4 {
            let k = m[row][col] / m[col][col];
            for c in col..4 {
                m[row][c] -= k * m[col][c];
            }
            b[row] -= k * b[col];
        }
    }
    let mut x = [0.0; 4];
    for row in (0..4).rev() {
        let rest: f64 = (row + 1..4).map(|c| m[row][c] * x[c]).sum();
        x[row] = (b[row] - rest) / m[row][row];
    }
    Some(x)
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
    fn lines_and_curves_meeting_surfaces() {
        use crate::geometry::curve::{Curve, TrimmedCurve};
        use crate::geometry::surface::Torus;
        let x_axis = Curve::Line(Line::new(Point::new(-10.0, 0.0, 0.0), Direction::X));
        // A sphere and a cylinder in closed form.
        let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0));
        let found = curve_surface(&x_axis, &sphere, 1e-7);
        let xs: Vec<f64> = found.points.iter().map(|p| p.point.x).collect();
        assert_eq!(xs.len(), 2);
        assert!((xs[0] + 2.0).abs() < 1e-12 && (xs[1] - 2.0).abs() < 1e-12, "{xs:?}");
        let along = cylinder(Point::ORIGIN, Direction::X, 1.0);
        let on = Curve::Line(Line::new(Point::new(0.0, 1.0, 0.0), Direction::X));
        assert_eq!(curve_surface(&on, &along, 1e-7).segments.len(), 1);
        // A torus (R 3, r 1) four times along a diameter, followed in steps; a bounded
        // piece of the line only twice.
        let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0));
        let found = curve_surface(&x_axis, &torus, 1e-7);
        let xs: Vec<f64> = found.points.iter().map(|p| p.point.x).collect();
        assert_eq!(xs.len(), 4, "{xs:?}");
        for (x, want) in xs.iter().zip([-4.0, -2.0, 2.0, 4.0]) {
            assert!((x - want).abs() < 1e-9, "{xs:?}");
        }
        for p in &found.points {
            assert!(torus.point(p.u, p.v).distance(p.point) < 1e-9);
        }
        let piece = Curve::Trimmed(TrimmedCurve::new(x_axis.clone(), 5.0, 9.0));
        assert_eq!(curve_surface(&piece, &torus, 1e-7).points.len(), 2);
    }

    #[test]
    fn a_torus_meets_planes_and_coaxial_surfaces_in_circles() {
        use crate::geometry::surface::Torus;
        let torus = Surface::Torus(Torus::new(Frame::WORLD, 10.0, 2.0));
        let radii = |q: QuadricIntersection| match q {
            QuadricIntersection::Circles(c) => {
                let mut r: Vec<f64> = c.iter().map(|c| c.radius()).collect();
                r.sort_by(f64::total_cmp);
                r
            }
            other => panic!("{other:?}"),
        };
        // Across the axis at height 1: radii 10 ± √3.
        let across = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 1.0), Direction::Z));
        let r = radii(quadrics(&torus, &across, 1e-7));
        assert!((r[0] - (10.0 - 3f64.sqrt())).abs() < 1e-12 && (r[1] - (10.0 + 3f64.sqrt())).abs() < 1e-12, "{r:?}");
        // Through the axis: the two meridian circles.
        let through = Surface::Plane(Plane::new(Point::ORIGIN, Direction::Y));
        let QuadricIntersection::Circles(c) = quadrics(&torus, &through, 1e-7) else { panic!() };
        assert_eq!(c.len(), 2);
        for c in &c {
            assert!((c.radius() - 2.0).abs() < 1e-12 && (c.center().distance(Point::ORIGIN) - 10.0).abs() < 1e-12);
        }
        // A coaxial cylinder of radius 11 at heights ±√3, a sphere about the centre.
        let r = radii(quadrics(&cylinder(Point::ORIGIN, Direction::Z, 11.0), &torus, 1e-7));
        assert!(r.len() == 2 && r.iter().all(|r| (r - 11.0).abs() < 1e-12));
        let ball = Surface::Sphere(Sphere::new(Frame::WORLD, 10.0));
        let QuadricIntersection::Circles(c) = quadrics(&torus, &ball, 1e-7) else { panic!() };
        for c in &c {
            for i in 0..8 {
                let p = c.point(i as f64);
                assert!(distance(&torus, p) < 1e-9 && distance(&ball, p) < 1e-9);
            }
        }
        // Off the axis: not in closed form.
        let off = cylinder(Point::new(1.0, 0.0, 0.0), Direction::Z, 11.0);
        assert_eq!(quadrics(&torus, &off, 1e-7), QuadricIntersection::NotAnalytic);
    }

    #[test]
    fn an_extrusion_meets_a_plane_along_lines_and_points_lie_on_both_surfaces() {
        use crate::geometry::curve::Curve;
        use crate::geometry::surface::ExtrusionSurface;
        let arc = Curve::Circle(Circle::new(Frame::WORLD, 5.0));
        let tube = ExtrusionSurface::new(arc, Direction::Z);
        let plane = Surface::Plane(Plane::new(Point::new(3.0, 0.0, 0.0), Direction::X));
        let lines = extrusion_lines(&tube, &plane, 1e-7).unwrap();
        assert_eq!(lines.len(), 2);
        for l in &lines {
            assert!((l.origin.x - 3.0).abs() < 1e-9 && (l.origin.y.abs() - 4.0).abs() < 1e-9);
        }
        let across = Surface::Plane(Plane::new(Point::ORIGIN, Direction::Z));
        assert!(extrusion_lines(&tube, &across, 1e-7).is_none());
        // Two crossing cylinders: a point near their curve moved onto it.
        let (a, b) = (cylinder(Point::ORIGIN, Direction::Z, 2.0), cylinder(Point::ORIGIN, Direction::X, 2.0));
        let (p, _) = surfaces_point(&a, &b, Point::new(1.0, 1.5, 1.4), None).unwrap();
        assert!(distance(&a, p) < 1e-12 && distance(&b, p) < 1e-12, "{p:?}");
    }

    #[test]
    fn planes_meet_a_cone_in_conics() {
        // Planes at every slant to a cone's axis (half-angle 30°): each conic found lies on
        // both, and is of the kind the slant makes.
        use crate::geometry::surface::Cone;
        let cone = Surface::Cone(Cone::new(Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.2, 0.1, 1.0), Direction::X), std::f64::consts::FRAC_PI_6, 4.0));
        let on_cone = |p: Point| {
            let Surface::Cone(c) = &cone else { unreachable!() };
            let z = Vector::from(c.frame().z());
            let w = p - c.apex();
            let (sin, cos) = c.semi_angle_sin_cos();
            (w.dot(z).abs() * sin.abs() - (w - z * w.dot(z)).length() * cos).abs()
        };
        for (tilt, kind) in [(1.2, "ellipse"), (std::f64::consts::FRAC_PI_6, "parabola"), (0.4, "hyperbola"), (0.0, "hyperbola")] {
            // The plane's normal at `tilt` from the axis, through a point off the apex.
            let Surface::Cone(c) = &cone else { unreachable!() };
            let f = c.frame();
            let n = Vector::from(f.z()) * tilt.sin() + Vector::from(f.x()) * tilt.cos();
            let plane = Surface::Plane(Plane::new(c.apex() + Vector::from(f.z()) * -5.0 + Vector::from(f.x()) * 0.7, Direction::from_vector(n)));
            let curves: Vec<crate::geometry::curve::Curve> = match quadrics(&plane, &cone, 1e-7) {
                QuadricIntersection::Ellipses(v) if kind == "ellipse" => v.into_iter().map(crate::geometry::curve::Curve::Ellipse).collect(),
                QuadricIntersection::Parabolas(v) if kind == "parabola" => v.into_iter().map(crate::geometry::curve::Curve::Parabola).collect(),
                QuadricIntersection::Hyperbolas(v) if kind == "hyperbola" => v.into_iter().map(crate::geometry::curve::Curve::Hyperbola).collect(),
                other => panic!("tilt {tilt}: {other:?}"),
            };
            for curve in &curves {
                for i in -10..=10 {
                    let p = curve.point(0.3 * i as f64);
                    assert!(on_cone(p) < 1e-9 * (1.0 + (p - Point::ORIGIN).length()), "tilt {tilt}: {p:?} off the cone by {}", on_cone(p));
                    assert!(distance(&plane, p) < 1e-9 * (1.0 + (p - Point::ORIGIN).length()), "tilt {tilt}: {p:?} off the plane");
                }
            }
        }
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
        let QuadricIntersection::Parabolas(parabolas) = quadrics(&parallel, &cone, 1e-7) else { panic!() };
        for p in parabolas.iter().flat_map(|c| (-5..=5).map(move |i| crate::geometry::curve::Curve::Parabola(c.clone()).point(0.5 * i as f64))) {
            assert!(distance(&cone, p) < 1e-9 && distance(&parallel, p) < 1e-9, "{p:?}");
        }
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
    fn cylinders_touching_far_from_the_axis_origin() {
        // Square cylinders tangent at one point (two fillets' bottoms on one plane), the
        // base one's axis starting far off: both branches pass through that point exactly
        // (the double root found from the other's own frame, not from the discriminant).
        let base = cylinder(Point::new(-249.2, -2.0, 7.0), Direction::X, 8.0);
        let other = cylinder(Point::new(74.9, 70.0, 9.0), -Direction::Y, 10.0);
        let touch = Point::new(74.9, -2.0, -1.0);
        let curves = ruled_quadric(&base, &other, 1e-7, None).unwrap();
        assert_eq!(curves.curves.len(), 2);
        for c in &curves.curves {
            let (a, b) = c.domain();
            let gap = |t: f64| c.point(t).map_or(f64::INFINITY, |p| p.distance(touch));
            let n = 4096;
            let step = (b - a) / n as f64;
            let i = (0..=n).min_by(|&i, &j| gap(a + step * i as f64).total_cmp(&gap(a + step * j as f64))).unwrap();
            // (Golden sections about the nearest sample.)
            let (mut lo, mut hi) = (a + step * (i as f64 - 1.0), a + step * (i as f64 + 1.0));
            for _ in 0..80 {
                let (x1, x2) = (hi - 0.618 * (hi - lo), lo + 0.618 * (hi - lo));
                if gap(x1) < gap(x2) {
                    hi = x2;
                } else {
                    lo = x1;
                }
            }
            let nearest = gap(0.5 * (lo + hi));
            assert!(nearest < 1e-7, "{nearest}");
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
