//! The continuity of two faces across an edge and an edge's valid range, as OCCT's BRepLib
//! finds them (DRAW's `getedgeregularity` and `validrange`): ports of
//! `BRepLib::ContinuityOfFaces` (with GeomLProp's principal curvatures) and
//! `BRepLib::FindValidRange`.

use kentron::geometry::curve2::Curve2;
use kentron::geometry::surface::Surface;
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Point, Transformable, Vector};
use kentron::topology::{CoEdgeId, Model, NodeId, Orientation, Shape};

/// OCCT's Precision::Confusion.
const CONFUSION: f64 = 1e-7;

/// The coedges of `edge` in `face` (two where it is the face's seam).
fn coedges(model: &Model, face: Shape, edge: Shape) -> Vec<CoEdgeId> {
    let forward = Shape { orientation: Orientation::Forward, ..face };
    let mut out = Vec::new();
    for wire in model.children(forward) {
        for coedge in model.children(wire) {
            if let Ok(id) = CoEdgeId::try_from(coedge.node) {
                if NodeId::from(model[id].edge) == edge.node {
                    out.push(id);
                }
            }
        }
    }
    out
}

fn elementary(s: &Surface) -> bool {
    matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_))
}

/// A face's surface at a point of an edge's pcurve on it (OCCT's SurfaceProperties).
struct Side<'a> {
    surface: &'a Surface,
    pcurve: Curve2,
    reversed: bool,
}

impl Side<'_> {
    /// The point, the derivative across the curve (along the normal of its tangent in the
    /// parameters, turned for a reversed face) and the unit normal (turned too).
    fn at(&self, t: f64) -> (Point, Vector, Vector, kentron::geometry::SurfaceD2) {
        let (uv, tangent) = self.pcurve.d1(t);
        let d = self.surface.d2(uv.x, uv.y);
        let length = (tangent.x * tangent.x + tangent.y * tangent.y).sqrt();
        let across = if length < CONFUSION {
            Vector::ZERO
        } else {
            let (ox, oy) = (-tangent.y / length, tangent.x / length);
            let (ox, oy) = if self.reversed { (-ox, -oy) } else { (ox, oy) };
            d.du * ox + d.dv * oy
        };
        let n = d.du.cross(d.dv);
        let n = n * (1.0 / n.length().max(1e-300));
        (d.point, across, if self.reversed { -n } else { n }, d)
    }
}

/// The principal curvatures (largest first) and their directions at a point (GeomLProp's
/// ComputeSurfCurvatures), signed for a reversed face; `None` where they are undefined.
fn curvatures(d: &kentron::geometry::SurfaceD2, reversed: bool) -> Option<[(Vector, f64); 2]> {
    let n = d.du.cross(d.dv);
    if n.length() < 1e-300 {
        return None;
    }
    let n = n * (1.0 / n.length());
    let (e, f, g) = (d.du.dot(d.du), d.du.dot(d.dv), d.dv.dot(d.dv));
    let (l, m, nn) = (n.dot(d.duu), n.dot(d.duv), n.dot(d.dvv));
    let (a0, b0, c0) = (e * m - f * l, e * nn - g * l, f * nn - g * m);
    let largest = a0.abs().max(b0.abs()).max(c0.abs());
    let sign = if reversed { -1.0 } else { 1.0 };
    let unit = |v: Vector| v * (1.0 / v.length().max(1e-300));
    if largest < f64::EPSILON {
        // (An umbilic point: one curvature, the directions along u and across it.)
        if g < f64::EPSILON {
            return None;
        }
        let k = nn / g * sign;
        return Some([(unit(d.du.cross(n)), k), (unit(d.du), k)]);
    }
    let (a, b, c) = (a0 / largest, b0 / largest, c0 / largest);
    let roots = |a: f64, b: f64, c: f64| -> Option<(f64, f64)> {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let q = -0.5 * (b + b.signum() * disc.sqrt());
        Some((q / a, if q != 0.0 { c / q } else { q / a }))
    };
    let (k1, v1, k2, v2) = if a.abs() > f64::EPSILON {
        let (r1, r2) = roots(a, b, c)?;
        let k = |r: f64| ((l * r + 2.0 * m) * r + nn) / ((e * r + 2.0 * f) * r + g);
        (k(r1), d.du * r1 + d.dv, k(r2), d.du * r2 + d.dv)
    } else if c.abs() > f64::EPSILON {
        let (r1, r2) = roots(c, b, a)?;
        let k = |r: f64| ((nn * r + 2.0 * m) * r + l) / ((g * r + 2.0 * f) * r + e);
        (k(r1), d.du + d.dv * r1, k(r2), d.du + d.dv * r2)
    } else {
        (l / e, d.du, nn / g, d.dv)
    };
    let (max, min) = if k1 < k2 { ((unit(v2), k2), (unit(v1), k1)) } else { ((unit(v1), k1), (unit(v2), k2)) };
    Some([(max.0, max.1 * sign), (min.0, min.1 * sign)])
}

fn parallel(a: Vector, b: Vector) -> bool {
    a.cross(b).length_squared() <= CONFUSION * CONFUSION
}

/// The continuity of `f1` and `f2` across `edge` ("C0", "G1", "C1", "G2", "C2" or "CN"), as
/// OCCT's BRepLib::ContinuityOfFaces with angular tolerance `angle`: sampled at 21 points
/// along the edge (its ends kept clear), the lowest found.
pub fn continuity_of_faces(model: &Model, edge: Shape, f1: Shape, f2: Shape, angle: f64) -> &'static str {
    let seam = f1.node == f2.node && f1.orientation == f2.orientation;
    let (c1, c2) = (coedges(model, f1, edge), coedges(model, f2, edge));
    let (first, second) = if f1.node == f2.node {
        match c1[..] {
            [a, b] => (a, b),
            _ => return "C0",
        }
    } else {
        match (c1.first(), c2.first()) {
            (Some(&a), Some(&b)) => (a, b),
            _ => return "C0",
        }
    };
    let (Some((p1, t0, t1)), Some((p2, _, _))) = (model.effective_pcurve(first), model.effective_pcurve(second)) else { return "C0" };
    let surface = |f: Shape| model.surface(f).map(|s| s.transformed(&f.placement));
    let (Some(s1), Some(s2)) = (surface(f1), surface(f2)) else { return "C0" };
    let both_elementary = elementary(&s1) && elementary(&s2);
    if seam && both_elementary {
        return "CN";
    }
    let a = Side { surface: &s1, pcurve: p1, reversed: f1.orientation == Orientation::Reversed };
    let b = Side { surface: &s2, pcurve: p2, reversed: f2.orientation == Orientation::Reversed };
    // C0 < G1 < C1 < G2 < C2 < CN.
    let mut cont = if both_elementary { 5 } else { 4 };
    let eps = (t1 - t0) / 100.0;
    let (f, l) = (t0 + eps, t1 - eps);
    for i in 0..=20 {
        if cont == 0 {
            break;
        }
        let t = f + (l - f) * i as f64 / 20.0;
        let (_, d1, n1, dd1) = a.at(t);
        let (_, d2, n2, dd2) = b.at(t);
        let (q1, q2) = (d1.length_squared(), d2.length_squared());
        if !(d1.cross(d2).length_squared() <= angle * angle * q1 * q2) {
            return "C0";
        }
        if n1.dot(n2) < 0.0 {
            return "C0";
        }
        let mut here = 1;
        if (q1.sqrt() - q2.sqrt()).abs() < CONFUSION && d1.dot(d2) > CONFUSION * CONFUSION {
            here = 2;
        }
        if cont >= 3 {
            let (Some(k1), Some(k2)) = (curvatures(&dd1, a.reversed), curvatures(&dd2, b.reversed)) else { return "C0" };
            for step in 0..2 {
                let (x0, y0, x1, y1) = (k1[0], k2[step], k1[1], k2[1 - step]);
                if parallel(x0.0, y0.0) && (x0.1 - y0.1).abs() < CONFUSION && parallel(x1.0, y1.0) && (x1.1 - y1.1).abs() < CONFUSION {
                    here = if here == 2 && x0.0.dot(y0.0) > CONFUSION && x1.0.dot(y1.0) > CONFUSION { 4 } else { 3 };
                    break;
                }
            }
        }
        cont = cont.min(here);
    }
    if both_elementary && cont == 4 {
        cont = 5;
    }
    ["C0", "G1", "C1", "G2", "C2", "CN"][cont]
}

/// The parameter change a 3D distance `tol` can need at most on `curve` (OCCT's
/// Adaptor3d_Curve::Resolution): for a B-spline, `tol` over the bound its poles set on its
/// derivative (BSplCLib::Resolution); for a circle, over its radius; else over its derivative's
/// largest length found.
fn resolution(curve: &kentron::geometry::curve::Curve, first: f64, last: f64, tol: f64) -> f64 {
    use kentron::geometry::curve::Curve;
    match curve {
        Curve::BSpline(b) => {
            let (poles, flat, p) = (b.poles(), b.flat_knots(), b.degree() as f64);
            let mut bound: f64 = 0.0;
            for i in 0..poles.len().saturating_sub(1) {
                let span = (flat[i + b.degree() + 1] - flat[i + 1]).abs();
                if span > 1e-300 {
                    bound = bound.max(poles[i + 1].distance(poles[i]) * p / span);
                }
            }
            if bound > 0.0 { tol / bound } else { last - first }
        }
        Curve::Circle(c) => tol / c.radius(),
        Curve::Line(_) => tol,
        _ => {
            let largest = (0..=64).map(|i| curve.d1(first + (last - first) * i as f64 / 64.0).1.length()).fold(0.0, f64::max);
            if largest > 0.0 { tol / largest } else { last - first }
        }
    }
}

/// The range of `edge`'s curve outside its vertices' tolerances (OCCT's
/// BRepLib::FindValidRange, stepping and bisecting as it does, its result its): `None` where
/// there is none (or an end of the curve lies outside its vertex's tolerance).
pub fn valid_range(model: &Model, edge: Shape) -> Option<(f64, f64)> {
    let (curve, first, last) = model.edge_curve(edge)?;
    let (v1, v2) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
    let (p1, p2) = (model.point(v1), model.point(v2));
    // (Widened by the confusion tolerance, as OCCT's are, "to have correspondence with
    // intersection precision".)
    let (tol1, tol2) = (model.tolerance(v1) + CONFUSION, model.tolerance(v2) + CONFUSION);
    if last - first < 1e-9 {
        return None;
    }
    let largest = first.abs().max(last.abs());
    let eps = (resolution(&curve, first, last, model.tolerance(edge)) * 0.1).max(f64::EPSILON * largest).max(1e-9);
    let spline = matches!(curve, kentron::geometry::curve::Curve::BSpline(_) | kentron::geometry::curve::Curve::Bezier(_));
    // (Where the curve's derivative is small, steps double, as OCCT's do on B-splines.)
    let small = if spline { (0.01 / resolution(&curve, first, last, 1.0)).powi(2) } else { 0.0 };
    // From one end, the parameter where the curve leaves the vertex's sphere (OCCT's
    // findNearestValidPoint).
    let leave = |from: f64, to: f64, p: Point, tol: f64| -> Option<f64> {
        let out = |t: f64| curve.point(t).distance_squared(p) > tol * tol;
        if out(from) {
            return None;
        }
        let forward = to > from;
        let mut step = (resolution(&curve, first, last, tol) * 1.01).max(eps);
        if !forward {
            step = -step;
        }
        let past = |t: f64| if forward { t > to } else { t < to };
        let (mut inside, mut outside) = (from, from);
        loop {
            inside = outside;
            outside += step;
            if past(outside) {
                if !out(to) {
                    return None;
                }
                outside = to;
                break;
            }
            if small > 0.0 {
                let mut local = step;
                loop {
                    let (_, d1) = curve.d1(outside);
                    if !out(outside) && d1.length_squared() < small {
                        local *= 2.0;
                        outside += local;
                        if !past(outside) {
                            continue;
                        }
                        outside = to;
                        if !out(outside) {
                            return None;
                        }
                    }
                    break;
                }
            }
            if out(outside) {
                break;
            }
        }
        while (outside - inside).abs() > eps {
            let mid = 0.5 * (inside + outside);
            if out(mid) {
                outside = mid;
            } else {
                inside = mid;
            }
        }
        Some(0.5 * (inside + outside))
    };
    let u1 = leave(first, last, p1, tol1)?;
    if last - u1 < eps {
        return None;
    }
    let u2 = leave(last, first, p2, tol2)?;
    if u2 - first < eps || u1 > u2 {
        return None;
    }
    Some((u1, u2))
}
