//! Projecting curves onto planes and surfaces.

use crate::math::{AffineTransform, AffineTransformable, Direction, Direction2, Frame, Frame2, Matrix3, Point, Point2, Vector, Vector2};
use crate::tolerance::ANGULAR;

use super::curve::{Circle, Curve, Ellipse, Hyperbola, Line, Parabola, TrimmedCurve};
use super::{ParametricCurve, ParametricSurface};
use super::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Line2};
use super::surface::{Cone, Cylinder, Plane, Sphere, Surface};

/// Projects `curve` onto `plane` along `direction`. Conics stay conics of the same kind
/// where possible (a projected trimmed curve stays trimmed); other curves become
/// B-splines. With `keep_parametrization`, the result has the same parametrization as
/// `curve`. `None` if `direction` is parallel to the plane.
#[doc(alias = "GeomProjLib::ProjectOnPlane")]
pub fn curve_on_plane(
    curve: &Curve,
    plane: &Plane,
    direction: Direction,
    keep_parametrization: bool,
) -> Option<Curve> {
    let n = Vector::from(plane.normal());
    let d = Vector::from(direction);
    let dn = d.dot(n);
    if dn.abs() <= ANGULAR {
        return None;
    }
    // P -> P - ((P - O)·n / (d·n)) d.
    let column = |e: Vector| e - d * (e.dot(n) / dn);
    let (cx, cy, cz) = (column(Vector::new(1.0, 0.0, 0.0)), column(Vector::new(0.0, 1.0, 0.0)), column(Vector::new(0.0, 0.0, 1.0)));
    let matrix = Matrix3::from_rows([[cx.x, cy.x, cz.x], [cx.y, cy.y, cz.y], [cx.z, cy.z, cz.z]]);
    let translation = d * ((plane.origin() - Point::ORIGIN).dot(n) / dn);
    let map = AffineTransform::new(matrix, translation);
    let (projected, reparametrization) = project_affine(curve, &map)?;
    if keep_parametrization && reparametrization.0 == 1.0 && reparametrization.1 != 0.0 {
        // A circle takes the shift of the parameter by turning its frame.
        if let Some(turned) = turned_circle(&projected, reparametrization.1) {
            return Some(turned);
        }
    }
    if keep_parametrization && reparametrization != (1.0, 0.0) {
        // The same parameters: a B-spline of the curve, projected pole by pole.
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        if !(a.is_finite() && b.is_finite() && a.abs() < 1e50 && b.abs() < 1e50) {
            return None;
        }
        let options = super::convert::CurveApproximationOptions { tolerance: crate::tolerance::LINEAR, ..Default::default() };
        let bspline = super::convert::approximate_curve(curve, &options)?.curve;
        return Some(Curve::BSpline(bspline.affine_transformed(&map)));
    }
    Some(projected)
}

/// A circle (or a trimmed one) `c` as the circle whose point at `t` is `c`'s at `t + shift`
/// (its frame turned by `shift`); `None` for another curve.
fn turned_circle(c: &Curve, shift: f64) -> Option<Curve> {
    match c {
        Curve::Circle(circle) => {
            let f = circle.frame();
            let (sin, cos) = shift.sin_cos();
            let x = Vector::from(f.x()) * cos + Vector::from(f.y()) * sin;
            Some(Curve::Circle(Circle::new(Frame::new(f.origin(), f.z(), Direction::try_from_vector(x)?), circle.radius())))
        }
        Curve::Trimmed(t) => {
            let (a, b) = t.range();
            Some(Curve::Trimmed(TrimmedCurve::new(turned_circle(t.basis(), shift)?, a - shift, b - shift)))
        }
        _ => None,
    }
}

/// The image of `curve` under the affine map `map`, as a curve of the same kind where
/// possible, with the parameter change `(scale, shift)`: the image of the point at `t` is
/// at `scale·t + shift`.
fn project_affine(curve: &Curve, map: &AffineTransform) -> Option<(Curve, (f64, f64))> {
    let point = |p: Point| map.transform_point(p);
    let vector = |v: Vector| map.transform_vector(v);
    let frame_of = |o: Point, x: Vector, y: Vector| {
        let z = Direction::try_from_vector(x.cross(y))?;
        Some(Frame::new(o, z, Direction::try_from_vector(x)?))
    };
    match curve {
        Curve::Line(l) => {
            let v = vector(Vector::from(l.direction));
            let scale = v.length();
            let dir = Direction::try_from_vector(v)?;
            Some((Curve::Line(Line::new(point(l.origin), dir)), (scale, 0.0)))
        }
        Curve::Circle(c) => {
            let f = c.frame();
            ellipse_from(point(f.origin()), vector(Vector::from(f.x()) * c.radius()), vector(Vector::from(f.y()) * c.radius()))
        }
        Curve::Ellipse(c) => {
            let f = c.frame();
            ellipse_from(point(f.origin()), vector(Vector::from(f.x()) * c.major_radius()), vector(Vector::from(f.y()) * c.minor_radius()))
        }
        Curve::Hyperbola(h) => {
            // O + a cosh t + b sinh t, with conjugate semi-diameters a and b made
            // perpendicular by a shift of the parameter.
            let f = h.frame();
            let a = vector(Vector::from(f.x()) * h.major_radius());
            let b = vector(Vector::from(f.y()) * h.minor_radius());
            let delta = 0.5 * (-2.0 * a.dot(b) / (a.dot(a) + b.dot(b))).atanh();
            let (sh, ch) = (delta.sinh(), delta.cosh());
            let (a2, b2) = (a * ch + b * sh, a * sh + b * ch);
            let frame = frame_of(point(f.origin()), a2, b2)?;
            Some((Curve::Hyperbola(Hyperbola::new(frame, a2.length(), b2.length())), (1.0, -delta)))
        }
        Curve::Parabola(pa) => {
            // O + t²/(4F) X + t Y, mapped: O' + t² x / (4F) + t y. With y = alpha x^ + beta m
            // (m perpendicular to x), s = |beta| t, completing the square.
            let f = pa.frame();
            let x = vector(Vector::from(f.x()));
            let y = vector(Vector::from(f.y()));
            let xl = x.length();
            let axis = Direction::try_from_vector(x)?;
            let alpha = y.dot(Vector::from(axis));
            let side = y - Vector::from(axis) * alpha;
            let beta = side.length();
            if beta <= 1e-15 * y.length().max(1.0) {
                return None;
            }
            let m = side * (1.0 / beta);
            // Along the axis: c s² + (alpha / beta) s, with c = |x| / (4 F beta²).
            let c = xl / (4.0 * pa.focal() * beta * beta);
            let k = alpha / beta;
            let shift = k / (2.0 * c);
            let vertex = point(f.origin()) + Vector::from(axis) * (-k * k / (4.0 * c)) + m * (-shift);
            let frame = frame_of(vertex, Vector::from(axis), m)?;
            Some((Curve::Parabola(Parabola::new(frame, 1.0 / (4.0 * c))), (beta, shift)))
        }
        Curve::Trimmed(t) => {
            let (projected, (scale, shift)) = project_affine(t.basis(), map)?;
            let (a, b) = t.range();
            let (a, b) = (scale * a + shift, scale * b + shift);
            Some((Curve::Trimmed(TrimmedCurve::new(projected, a.min(b), a.max(b))), (scale, shift)))
        }
        Curve::BSpline(b) => Some((Curve::BSpline(b.affine_transformed(map)), (1.0, 0.0))),
        Curve::Bezier(b) => Some((Curve::Bezier(b.affine_transformed(map)), (1.0, 0.0))),
        other => {
            let (a, b) = (other.first_parameter(), other.last_parameter());
            if !(a.is_finite() && b.is_finite() && a.abs() < 1e50 && b.abs() < 1e50) {
                return None;
            }
            let bspline = super::convert::curve_to_bspline(other).or_else(|| {
                let options = super::convert::CurveApproximationOptions { tolerance: crate::tolerance::LINEAR, ..Default::default() };
                super::convert::approximate_curve(other, &options).map(|x| x.curve)
            })?;
            Some((Curve::BSpline(bspline.affine_transformed(map)), (1.0, 0.0)))
        }
    }
}

/// The ellipse (or circle) `o + a cos t + b sin t`, and the parameter change to its
/// principal frame.
fn ellipse_from(o: Point, a: Vector, b: Vector) -> Option<(Curve, (f64, f64))> {
    // Principal axes at t = delta: tan 2·delta = 2 a·b / (|a|² - |b|²).
    let mut delta = 0.5 * (2.0 * a.dot(b)).atan2(a.dot(a) - b.dot(b));
    let (sin, cos) = delta.sin_cos();
    let (mut a2, mut b2) = (a * cos + b * sin, b * cos - a * sin);
    if a2.length() < b2.length() {
        // The major axis a quarter turn on.
        (a2, b2) = (b2, -a2);
        delta += std::f64::consts::FRAC_PI_2;
    }
    let z = Direction::try_from_vector(a.cross(b))?;
    let frame = Frame::new(o, z, Direction::try_from_vector(a2)?);
    // The frame's y must be along b2 (a2 × b2 has the orientation of a × b).
    let (major, minor) = (a2.length(), b2.length());
    let curve = if (major - minor).abs() <= 1e-12 * major {
        Curve::Circle(Circle::new(frame, major))
    } else {
        Curve::Ellipse(Ellipse::new(frame, major, minor))
    };
    Some((curve, (1.0, -delta)))
}

/// The exact curve in `surface`'s parameter space of a curve lying on an elementary
/// surface, when it is analytic: e.g. a circle around a cone's or cylinder's axis maps to
/// a line `v = const` (with direction `-u` if the circle runs the other way round). `None`
/// if there is no such representation.
#[doc(alias = "ProjLib_Cone")]
#[doc(alias = "ProjLib_Cylinder")]
#[doc(alias = "ProjLib_Plane")]
#[doc(alias = "ProjLib_Sphere")]
#[doc(alias = "ProjLib_Torus")]
pub fn analytic_pcurve(curve: &Curve, surface: &Surface) -> Option<Curve2> {
    match surface {
        Surface::Plane(plane) => curve_on_plane_parameters(curve, plane),
        Surface::Cylinder(cylinder) => curve_on_cylinder_parameters(curve, cylinder),
        Surface::Sphere(sphere) => curve_on_sphere_parameters(curve, sphere),
        Surface::Cone(cone) => curve_on_cone_parameters(curve, cone),
        // Others (tori, surfaces of revolution, ...) have no closed form here yet.
        _ => None,
    }
}

/// The direction about which a circle turns as its parameter grows: its frame's Z, or the
/// opposite for a left-handed (reflected) frame.
fn circle_axis(c: &crate::geometry::curve::Circle) -> Vector {
    Vector::from(c.frame().x()).cross(Vector::from(c.frame().y()))
}

/// `v = const` lines of circles about the axis, `u = const` lines of generators.
fn curve_on_cylinder_parameters(curve: &Curve, cylinder: &Cylinder) -> Option<Curve2> {
    let f = cylinder.frame();
    let z = Vector::from(f.z());
    let tol = crate::tolerance::LINEAR;
    let line = |o: (f64, f64), d: (f64, f64)| Some(Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1))));
    match curve {
        Curve::Line(l) => {
            let along = Vector::from(l.direction).dot(z);
            if (along.abs() - 1.0).abs() > ANGULAR.sqrt() || (cylinder_distance(cylinder, l.origin)).abs() > tol {
                return None;
            }
            let (u, v) = cylinder.parameters(l.origin);
            line((u, v), (0.0, along.signum()))
        }
        Curve::Circle(c) => {
            let normal = circle_axis(c).dot(z);
            let v = (c.center() - f.origin()).dot(z);
            let on_axis = (c.center() - f.origin()) - z * v;
            if (normal.abs() - 1.0).abs() > ANGULAR.sqrt() || on_axis.length() > tol || (c.radius() - cylinder.radius()).abs() > tol {
                return None;
            }
            let (u0, _) = cylinder.parameters(c.center() + Vector::from(c.frame().x()) * c.radius());
            line((u0, v), (normal.signum(), 0.0))
        }
        _ => None,
    }
}

/// `v = const` lines of parallels, `u = const` lines of meridians whose frame starts on
/// the equator and turns towards a pole (valid over latitudes in `[-π/2, π/2]`).
fn curve_on_sphere_parameters(curve: &Curve, sphere: &Sphere) -> Option<Curve2> {
    let f = sphere.frame();
    let z = Vector::from(f.z());
    let tol = crate::tolerance::LINEAR;
    let Curve::Circle(c) = curve else { return None };
    let line = |o: (f64, f64), d: (f64, f64)| Some(Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1))));
    let (cx, cy) = (Vector::from(c.frame().x()), Vector::from(c.frame().y()));
    let cz = circle_axis(c);
    let offset = c.center() - f.origin();
    let normal = cz.dot(z);
    if (normal.abs() - 1.0).abs() <= ANGULAR.sqrt() {
        // A parallel: centered on the axis, at the latitude of its points.
        let height = offset.dot(z);
        if (offset - z * height).length() > tol || (height.hypot(c.radius()) - sphere.radius()).abs() > tol {
            return None;
        }
        let (u0, v) = sphere.parameters(c.center() + cx * c.radius());
        return line((u0, v), (normal.signum(), 0.0));
    }
    // A meridian: a great circle through the poles, starting on the equator.
    let starts_on_equator = cx.dot(z).abs() <= ANGULAR.sqrt() && (cy.dot(z).abs() - 1.0).abs() <= ANGULAR.sqrt();
    if offset.length() > tol || (c.radius() - sphere.radius()).abs() > tol || !starts_on_equator {
        return None;
    }
    let (u0, _) = sphere.parameters(c.center() + cx * c.radius());
    line((u0, 0.0), (0.0, cy.dot(z).signum()))
}

/// `v = const` lines of circles about the axis, `u = const` lines of generators.
fn curve_on_cone_parameters(curve: &Curve, cone: &Cone) -> Option<Curve2> {
    let f = cone.frame();
    let z = Vector::from(f.z());
    let tol = crate::tolerance::LINEAR;
    let (sin, cos) = cone.semi_angle_sin_cos();
    let line = |o: (f64, f64), d: (f64, f64)| Some(Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1))));
    let on_cone = |p: Point| {
        let (u, v) = cone.parameters(p);
        let _ = u;
        let d = p - f.origin();
        let h = d.dot(z);
        let radial = (d - z * h).length();
        ((radial - cone.radius()) * cos - h * sin).abs() <= tol && v.is_finite()
    };
    match curve {
        Curve::Line(l) => {
            // A generator: away from the apex, its direction is ±D(u).
            let (u, v) = cone.parameters(l.origin);
            let d = Vector::from(l.direction);
            let generator = (Vector::from(f.x()) * u.cos() + Vector::from(f.y()) * u.sin()) * sin + z * cos;
            let along = d.dot(generator);
            if (along.abs() - 1.0).abs() > ANGULAR.sqrt() || !on_cone(l.origin) || cone.apex().distance(l.origin) <= tol {
                return None;
            }
            line((u, v), (0.0, along.signum()))
        }
        Curve::Circle(c) => {
            let normal = circle_axis(c).dot(z);
            let offset = c.center() - f.origin();
            let h = offset.dot(z);
            let start = c.center() + Vector::from(c.frame().x()) * c.radius();
            if (normal.abs() - 1.0).abs() > ANGULAR.sqrt() || (offset - z * h).length() > tol || !on_cone(start) {
                return None;
            }
            let (u0, v) = cone.parameters(start);
            line((u0, v), (normal.signum(), 0.0))
        }
        _ => None,
    }
}

/// Signed distance from `p` to the cylinder (positive outside).
fn cylinder_distance(cylinder: &Cylinder, p: Point) -> f64 {
    let f = cylinder.frame();
    let d = p - f.origin();
    let z = Vector::from(f.z());
    (d - z * d.dot(z)).length() - cylinder.radius()
}

/// The parameters of `p` (or of the nearest point) on `surface`: in closed form on elementary
/// surfaces, otherwise by inversion, from `near` (the previous point's) if given.
fn parameters_near(surface: &Surface, p: Point, near: Option<(f64, f64)>) -> Option<(f64, f64)> {
    if let Some(uv) = surface.parameters(p) {
        return Some(uv);
    }
    // Along a curve the next point's foot is near the last one's: Newton steps from there
    // (the point may lie off the surface, as in files accurate to 1e-3). From a foot on the
    // edge of a closed surface, also from the same place on its other edge (a curve
    // leaving the seam runs on from either side).
    let Some(start) = near else {
        let (u, v, _) = super::extrema::nearest_on_surface(p, surface, None);
        return Some((u, v));
    };
    let mut starts = vec![start];
    let (u0, u1, v0, v1) = surface.bounds();
    let other = |x: f64, lo: f64, hi: f64| {
        let at = 1e-9 * (hi - lo).abs().max(1.0);
        if (x - lo).abs() <= at { Some(hi) } else if (x - hi).abs() <= at { Some(lo) } else { None }
    };
    if surface.u_period().is_none() && surface.is_u_closed() && u0.is_finite() && u1.is_finite() {
        if let Some(u) = other(start.0, u0, u1) {
            starts.push((u, start.1));
        }
    }
    if surface.v_period().is_none() && surface.is_v_closed() && v0.is_finite() && v1.is_finite() {
        if let Some(v) = other(start.1, v0, v1) {
            starts.push((start.0, v));
        }
    }
    let mut found: Vec<(f64, f64, f64)> = starts.into_iter().map(|start| super::extrema::newton_on_surface(p, surface, start)).collect();
    // (From a collapsed edge, where one parameter is arbitrary, or where the steps led off
    // the surface, as from too far away, the whole surface searched.)
    let (u_free, v_free) = collapsed_at(surface, start.0, start.1);
    let off = found.iter().map(|f| f.2).fold(f64::INFINITY, f64::min) > 1e-6 * (1.0 + (p - Point::ORIGIN).length());
    if u_free || v_free || off {
        found.push(super::extrema::nearest_on_surface(p, surface, None));
    }
    let best = found.into_iter().min_by(|a, b| a.2.total_cmp(&b.2)).expect("a start");
    Some((best.0, best.1))
}

/// Samples of a curve in a surface's parameters where they lie on the edge of a closed (not
/// periodic) surface: put on the side of it their neighbours are on.
fn sides_of_seams(surface: &Surface, points: &mut [Point2]) {
    use super::ParametricSurface;
    let (u0, u1, v0, v1) = surface.bounds();
    let fix = |points: &mut [Point2], get: fn(&Point2) -> f64, set: fn(&mut Point2, f64), lo: f64, hi: f64| {
        let at = 1e-9 * (hi - lo).abs().max(1.0);
        let on_edge = |x: f64| (x - lo).abs() <= at || (x - hi).abs() <= at;
        for i in 0..points.len() {
            if !on_edge(get(&points[i])) {
                continue;
            }
            // The nearest sample off the edge, before or after.
            let neighbour = (1..points.len()).find_map(|k| {
                [i.checked_sub(k), Some(i + k)].into_iter().flatten().filter(|&j| j < points.len()).map(|j| get(&points[j])).find(|&x| !on_edge(x))
            });
            // (All on the edge, as a curve along the seam: on the first one's side.)
            let x = neighbour.unwrap_or(if i > 0 { get(&points[i - 1]) } else { get(&points[0]) });
            set(&mut points[i], if (x - lo).abs() < (x - hi).abs() { lo } else { hi });
        }
    };
    if surface.u_period().is_none() && surface.is_u_closed() && u0.is_finite() && u1.is_finite() {
        fix(points, |p| p.x, |p, x| p.x = x, u0, u1);
    }
    if surface.v_period().is_none() && surface.is_v_closed() && v0.is_finite() && v1.is_finite() {
        fix(points, |p| p.y, |p, y| p.y = y, v0, v1);
    }
}

/// Where `(u, v)` lies on a collapsed edge of a free-form surface (all of an iso-line one
/// point, as a sphere's pole made a B-spline): whether `u` is arbitrary there (the line of
/// constant `v` collapses), and whether `v` is.
fn collapsed_at(surface: &Surface, u: f64, v: f64) -> (bool, bool) {
    use super::ParametricSurface;
    if !matches!(surface, Surface::BSpline(_) | Surface::Bezier(_)) {
        return (false, false);
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let size = surface.point(u0, v0).distance(surface.point(u1, v1)) + surface.point(u0, v1).distance(surface.point(u1, v0));
    let tiny = 1e-9 * (1.0 + size);
    let same = |a: Point, b: Point, c: Point| a.distance(b) <= tiny && a.distance(c) <= tiny;
    let u_free = same(surface.point(u0, v), surface.point(u1, v), surface.point(0.5 * (u0 + u1), v));
    let v_free = same(surface.point(u, v0), surface.point(u, v1), surface.point(u, 0.5 * (v0 + v1)));
    (u_free, v_free)
}

/// Samples on a collapsed edge of a free-form surface ([`collapsed_at`]) given the arbitrary
/// parameter of their neighbour along the curve.
fn off_collapsed_edges(surface: &Surface, points: &mut [Point2]) {
    let n = points.len();
    let free: Vec<(bool, bool)> = points.iter().map(|p| collapsed_at(surface, p.x, p.y)).collect();
    for i in 0..n {
        let (u_free, v_free) = free[i];
        if !(u_free || v_free) {
            continue;
        }
        let neighbour = (1..n).find_map(|k| [i.checked_sub(k), Some(i + k)].into_iter().flatten().find(|&j| j < n && free[j] == (false, false)));
        if let Some(j) = neighbour {
            if u_free {
                points[i].x = points[j].x;
            }
            if v_free {
                points[i].y = points[j].y;
            }
        }
    }
}

/// The curve in `surface`'s parameter space of `curve` over `[first, last]` (lying on the
/// surface), as a B-spline with the same parametrization, accurate to `tolerance` in 3D.
/// The periodic `u` parameter is unwrapped so that the result is continuous, starting in
/// `[0, 2π)`. `None` if the surface has no closed-form parameters or the accuracy cannot
/// be reached.
#[doc(alias = "ProjLib_ComputeApprox")]
pub fn pcurve_on_range(curve: &Curve, surface: &Surface, first: f64, last: f64, tolerance: f64) -> Option<BSplineCurve2> {
    pcurve_on_range_within(curve, surface, first, last, &[tolerance])
}

/// Like [`pcurve_on_range`] to the first of `tolerances` (tightest first) that can be
/// reached: what trying them in turn gives, with the curve sampled and its points inverted
/// once for all of them.
pub fn pcurve_on_range_within(curve: &Curve, surface: &Surface, first: f64, last: f64, tolerances: &[f64]) -> Option<BSplineCurve2> {
    use super::{ParametricCurve, ParametricSurface};
    let u_period = surface.u_period();
    let v_period = surface.v_period();
    // Where u is undefined (at a sphere's poles) it is taken from the neighbouring sample.
    let singular = |v: f64| match surface {
        Surface::Sphere(_) => (v.abs() - std::f64::consts::FRAC_PI_2).abs() <= 1e-9,
        Surface::Cone(c) => (v + c.radius() / c.semi_angle_sin_cos().0).abs() <= 1e-9,
        _ => false,
    };
    let loosest = tolerances.iter().copied().fold(0.0, f64::max);
    // The first fit reaching each tolerance.
    let mut reached: Vec<Option<BSplineCurve2>> = vec![None; tolerances.len()];
    let first_reached = |reached: Vec<Option<BSplineCurve2>>| reached.into_iter().flatten().next();
    // Each doubling keeps the samples so far (their parameters come out the same) and adds
    // the middles, where the last fit was checked.
    let mut raw: Vec<Point2> = Vec::new();
    let mut middles: Vec<Point> = Vec::new();
    // How far the curve is off the surface at the samples.
    let mut offset: f64 = 0.0;
    let mut intervals = 8;
    while intervals <= 1024 {
        let params: Vec<f64> = (0..=intervals).map(|i| first + (last - first) * i as f64 / intervals as f64).collect();
        let mut next: Vec<Point2> = Vec::with_capacity(params.len());
        for (i, &t) in params.iter().enumerate() {
            let uv = match raw.is_empty() {
                false if i % 2 == 0 => raw[i / 2],
                _ => {
                    let p = if raw.is_empty() { curve.point(t) } else { middles[i / 2] };
                    let Some((u, v)) = parameters_near(surface, p, next.last().map(|q| (q.x, q.y))) else {
                        return first_reached(reached);
                    };
                    let off = surface.point(u, v).distance(p);
                    offset = if off.is_nan() { f64::INFINITY } else { offset.max(off) };
                    Point2::new(u, v)
                }
            };
            next.push(uv);
        }
        raw = next;
        // A fit is off by as much as the curve is off the surface at its samples, and finer
        // fits check near them: tolerances well below that are out of reach.
        let in_reach = |tolerance: f64| 2.0 * tolerance >= offset;
        let Some(open) = tolerances.iter().zip(&reached).position(|(&tolerance, r)| r.is_none() && in_reach(tolerance)) else {
            return first_reached(reached);
        };
        if reached[..open].iter().any(Option::is_some) {
            return first_reached(reached);
        }
        let mut points = raw.clone();
        let n = points.len();
        for i in 0..n {
            if singular(points[i].y) {
                let neighbour = if i + 1 < n && !singular(points[i + 1].y) { i + 1 } else { i.saturating_sub(1) };
                points[i].x = points[neighbour].x;
            }
        }
        off_collapsed_edges(surface, &mut points);
        sides_of_seams(surface, &mut points);
        if let Some(period) = u_period {
            for i in 1..n {
                points[i].x += ((points[i - 1].x - points[i].x) / period).round() * period;
            }
        }
        if let Some(period) = v_period {
            for i in 1..n {
                points[i].y += ((points[i - 1].y - points[i].y) / period).round() * period;
            }
        }
        let Ok(pcurve) = super::fit::interpolate_with_parameters(&points, &params) else {
            return first_reached(reached);
        };
        middles.clear();
        let mut error: f64 = 0.0;
        for w in params.windows(2) {
            let t = 0.5 * (w[0] + w[1]);
            let p = pcurve.point(t);
            let q = curve.point(t);
            middles.push(q);
            let d = surface.point(p.x, p.y).distance(q);
            error = if d.is_nan() { f64::INFINITY } else { error.max(d) };
        }
        if error <= loosest {
            for (k, &tolerance) in tolerances.iter().enumerate() {
                if reached[k].is_none() && error <= tolerance {
                    reached[k] = Some(pcurve.clone());
                }
            }
        }
        // Done once the tightest tolerance still in reach is met.
        if tolerances.iter().zip(&reached).find(|(&tolerance, _)| in_reach(tolerance)).is_some_and(|(_, r)| r.is_some()) {
            break;
        }
        intervals *= 2;
    }
    first_reached(reached)
}

/// Like [`pcurve_on_range`], but a degree-1 B-spline through samples placed adaptively
/// (where the curve's parameters change fast, e.g. passing close to a pole), each chord
/// within `tolerance` in 3D: a last resort where smooth fits fail.
pub fn polyline_pcurve_on_range(curve: &Curve, surface: &Surface, first: f64, last: f64, tolerance: f64) -> Option<BSplineCurve2> {
    use super::{ParametricCurve, ParametricSurface};
    let singular = |v: f64| match surface {
        Surface::Sphere(_) => (v.abs() - std::f64::consts::FRAC_PI_2).abs() <= 1e-9,
        Surface::Cone(c) => (v + c.radius() / c.semi_angle_sin_cos().0).abs() <= 1e-9,
        _ => false,
    };
    let u_period = surface.u_period();
    let previous = std::cell::Cell::new(None);
    // (None where the curve itself lies off the surface by more than the tolerance: no
    // chords through its feet come within it.)
    let uv = |t: f64| -> Option<Point2> {
        let p = curve.point(t);
        let (u, v) = parameters_near(surface, p, previous.get())?;
        if !singular(v) && !(surface.point(u, v).distance(p) <= tolerance) {
            return None;
        }
        previous.set(Some((u, v)));
        Some(Point2::new(u, v))
    };
    let v_period = surface.v_period();
    let unwrap = |p: Point2, near: Point2| {
        let x = u_period.map_or(p.x, |period| p.x + ((near.x - p.x) / period).round() * period);
        let y = v_period.map_or(p.y, |period| p.y + ((near.y - p.y) / period).round() * period);
        Point2::new(x, y)
    };
    // Bisect each stretch until the chord's middle matches the curve.
    let mut params = vec![first];
    let mut points = vec![uv(first)?];
    let mut stack = vec![(last, 0usize)];
    while let Some((t1, depth)) = stack.pop() {
        let (t0, p0) = (*params.last().expect("start"), *points.last().expect("start"));
        let p1 = unwrap(uv(t1)?, p0);
        let tm = 0.5 * (t0 + t1);
        let mid = Point2::new(0.5 * (p0.x + p1.x), 0.5 * (p0.y + p1.y));
        let good = singular(p0.y) || singular(p1.y) || surface.point(mid.x, mid.y).distance(curve.point(tm)) <= tolerance;
        if good || depth >= 18 {
            params.push(t1);
            points.push(p1);
            if points.len() > 20_000 {
                return None;
            }
        } else {
            stack.push((t1, depth + 1));
            stack.push((tm, depth + 1));
        }
    }
    // At a pole u is that of the neighbouring sample.
    let n = points.len();
    for i in 0..n {
        if singular(points[i].y) {
            let neighbour = if i + 1 < n && !singular(points[i + 1].y) { i + 1 } else { i.saturating_sub(1) };
            points[i].x = points[neighbour].x;
        }
    }
    let knots: Vec<f64> = params.clone();
    let mut multiplicities = vec![1; n];
    multiplicities[0] = 2;
    multiplicities[n - 1] = 2;
    BSplineCurve2::try_new(points, knots, multiplicities, 1).ok()
}

/// The curve in the plane's (u, v) space, with the same parametrization, of a line or
/// circle lying in (a plane parallel to) `plane`; `None` if it does not.
fn curve_on_plane_parameters(curve: &Curve, plane: &Plane) -> Option<Curve2> {
    let (x, y) = (Vector::from(plane.frame().x()), Vector::from(plane.frame().y()));
    // A unit vector parallel to the plane maps to a unit 2D direction.
    let in_plane = |d: Direction| {
        let d = Vector::from(d);
        let projected = Vector2::new(d.dot(x), d.dot(y));
        ((projected.length() - 1.0).abs() <= ANGULAR.sqrt()).then(|| Direction2::from_vector(projected))
    };
    let uv = |p: Point| {
        let (u, v) = plane.parameters(p);
        Point2::new(u, v)
    };
    match curve {
        Curve::Line(line) => Some(Curve2::Line(Line2::new(uv(line.origin), in_plane(line.direction)?))),
        Curve::Circle(circle) => {
            let f = circle.frame();
            let frame = Frame2::from_axes(uv(f.origin()), in_plane(f.x())?, in_plane(f.y())?)?;
            Some(Curve2::Circle(Circle2::new(frame, circle.radius())))
        }
        Curve::Ellipse(ellipse) => {
            let f = ellipse.frame();
            let frame = Frame2::from_axes(uv(f.origin()), in_plane(f.x())?, in_plane(f.y())?)?;
            Some(Curve2::Ellipse(Ellipse2::new(frame, ellipse.major_radius(), ellipse.minor_radius())))
        }
        Curve::BSpline(b) => {
            // An affine image of a B-spline is the B-spline of the images of its poles.
            if b.poles().iter().any(|&p| plane.distance(p) > crate::tolerance::LINEAR) {
                return None;
            }
            let poles: Vec<Point2> = b.poles().iter().map(|&p| uv(p)).collect();
            let (knots, multiplicities) = (b.knots().to_vec(), b.multiplicities().to_vec());
            let curve = if crate::geometry::ParametricCurve::is_periodic(b) {
                BSplineCurve2::try_new_periodic(poles, knots, multiplicities, b.degree()).ok()?
            } else {
                BSplineCurve2::try_new(poles, knots, multiplicities, b.degree()).ok()?
            };
            Some(Curve2::BSpline(match b.weights() {
                Some(w) => curve.with_weights(w.to_vec()),
                None => curve,
            }))
        }
        Curve::Bezier(b) => {
            if b.poles().iter().any(|&p| plane.distance(p) > crate::tolerance::LINEAR) {
                return None;
            }
            let poles: Vec<Point2> = b.poles().iter().map(|&p| uv(p)).collect();
            Some(Curve2::Bezier(match b.weights() {
                Some(w) => crate::geometry::curve2::BezierCurve2::rational(poles, w.to_vec()),
                None => crate::geometry::curve2::BezierCurve2::new(poles),
            }))
        }
        Curve::Parabola(c) => {
            let f = c.frame();
            let frame = Frame2::from_axes(uv(f.origin()), in_plane(f.x())?, in_plane(f.y())?)?;
            Some(Curve2::Parabola(crate::geometry::curve2::Parabola2::new(frame, c.focal())))
        }
        Curve::Hyperbola(c) => {
            let f = c.frame();
            let frame = Frame2::from_axes(uv(f.origin()), in_plane(f.x())?, in_plane(f.y())?)?;
            Some(Curve2::Hyperbola(crate::geometry::curve2::Hyperbola2::new(frame, c.major_radius(), c.minor_radius())))
        }
        _ => None,
    }
}

/// Approximates the curve in `surface`'s parameter space of a curve lying on it (or close
/// to it), within `tolerance`, unwrapping periodic parameters so that the result is
/// continuous. The result has the same parametrization as `curve`.
#[doc(alias = "ProjLib_ComputeApproxOnPolarSurface")]
#[doc(alias = "ProjLib_ComputeApprox")]
pub fn pcurve(curve: &Curve, surface: &Surface, tolerance: f64) -> Option<BSplineCurve2> {
    fitted_pcurve(curve, surface, tolerance, false)
}

/// The orthogonal projection of `curve` onto `surface` (each point to its nearest on the
/// surface) in the surface's parameter space, within `tolerance` of the projected points:
/// as [`pcurve`], but for a curve off the surface too. The result has the same
/// parametrization as `curve`.
#[doc(alias = "GeomProjLib::Curve2d")]
#[doc(alias = "ProjLib_ProjectedCurve")]
pub fn projected_pcurve(curve: &Curve, surface: &Surface, tolerance: f64) -> Option<BSplineCurve2> {
    fitted_pcurve(curve, surface, tolerance, true)
}

/// [`pcurve`], or with `projecting`, [`projected_pcurve`].
fn fitted_pcurve(curve: &Curve, surface: &Surface, tolerance: f64, projecting: bool) -> Option<BSplineCurve2> {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if !(a.is_finite() && b.is_finite() && a.abs() < 1e50 && b.abs() < 1e50 && a < b) {
        return None;
    }
    // An exact pcurve where there is one (of a curve on the surface).
    if let Some(exact) = analytic_pcurve(curve, surface).filter(|_| !projecting) {
        let trimmed = Curve2::Trimmed(super::curve2::TrimmedCurve2::new(exact, a, b));
        if let Some(bspline) = super::convert::curve2_to_bspline(&trimmed) {
            return Some(bspline);
        }
    }
    // Otherwise samples projected onto the surface (each from the last), unwrapped across
    // periodic seams, and fitted on more and more spans.
    let mut cuts: Vec<f64> = vec![a, b];
    cuts.extend(curve.breaks().into_iter().filter(|&k| a < k && k < b));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let unwrap = |x: f64, previous: f64, period: Option<f64>| match period {
        Some(p) => x + ((previous - x) / p).round() * p,
        None => x,
    };
    let (u_period, v_period) = (surface.u_period().filter(|_| surface.is_u_periodic()), surface.v_period().filter(|_| surface.is_v_periodic()));
    // The parameters of the curve's points at `ts` (increasing): each projected from the
    // last, unwrapped across periodic seams. On a closed surface that is not periodic, the
    // first point may be on the seam, on either side: the second is projected afresh, and an
    // end across the seam from its neighbour moved to the other side.
    let closed_open = [surface.is_u_closed() && u_period.is_none(), surface.is_v_closed() && v_period.is_none()];
    let project_along = |ts: &[f64]| -> Vec<[f64; 3]> {
        let mut uv: Vec<[f64; 3]> = Vec::with_capacity(ts.len());
        let mut hint: Option<(f64, f64)> = None;
        let mut last: Option<(Point, f64)> = None;
        for (k, &t) in ts.iter().enumerate() {
            let seed = if k == 1 && closed_open.iter().any(|&c| c) { None } else { hint };
            let p = curve.point(t);
            let (mut u, mut v, mut d) = super::extrema::nearest_on_surface(p, surface, seed);
            // (The distance changes by no more than the point moves: further, the search from
            // the last foot was caught (at a pole, say) and the whole surface is searched.)
            if let Some((q, before)) = last {
                if seed.is_some() && d > before + q.distance(p) + 1e-9 * (1.0 + before) {
                    let (gu, gv, gd) = super::extrema::nearest_on_surface(p, surface, None);
                    if gd < d {
                        (u, v, d) = (gu, gv, gd);
                    }
                }
            }
            last = Some((p, d));
            if let Some((pu, pv)) = hint {
                u = unwrap(u, pu, u_period);
                v = unwrap(v, pv, v_period);
            }
            hint = Some((u, v));
            uv.push([u, v, 0.0]);
        }
        // Where the surface collapses (a surface of revolution's axis, a pole), u (or v) is
        // arbitrary there: the neighbouring sample's.
        let collapsed: Vec<[bool; 2]> = uv
            .iter()
            .map(|q| {
                // (Next to nothing beside the other derivative.)
                let (_, du, dv) = surface.d1(q[0], q[1]);
                let (lu, lv) = (du.length(), dv.length());
                [lu <= 1e-6 * lv, lv <= 1e-6 * lu]
            })
            .collect();
        for k in 0..2 {
            for i in 0..uv.len() {
                if !collapsed[i][k] {
                    continue;
                }
                let neighbour = (1..uv.len()).find_map(|d| [i.checked_sub(d), Some(i + d)].into_iter().flatten().find(|&j| j < uv.len() && !collapsed[j][k]));
                if let Some(j) = neighbour {
                    uv[i][k] = uv[j][k];
                }
            }
        }
        let (su0, su1, sv0, sv1) = surface.bounds();
        for (is_closed, lo, hi, k) in [(closed_open[0], su0, su1, 0usize), (closed_open[1], sv0, sv1, 1)] {
            if !is_closed || uv.len() < 3 {
                continue;
            }
            let last = uv.len() - 1;
            for (end, next) in [(0, 1), (last, last - 1)] {
                if (uv[end][k] - uv[next][k]).abs() > 0.5 * (hi - lo) {
                    uv[end][k] = lo + hi - uv[end][k];
                }
            }
        }
        uv
    };
    // Where the curve crosses a knot line of a B-spline surface its parameters there need not
    // be smooth (a rational circle's parametrization is not): the fit is cut there too.
    let knot_lines = match surface {
        Surface::BSpline(s) => Some((s.u_knots().to_vec(), s.v_knots().to_vec())),
        Surface::Trimmed(t) => match t.basis() {
            Surface::BSpline(s) => Some((s.u_knots().to_vec(), s.v_knots().to_vec())),
            _ => None,
        },
        _ => None,
    };
    if let Some((uk, vk)) = knot_lines {
        let span = |knots: &[f64], x: f64| knots.partition_point(|&k| k <= x + 1e-12 * (1.0 + x.abs()));
        const N: usize = 256;
        let ts: Vec<f64> = (0..=N).map(|i| a + (b - a) * i as f64 / N as f64).collect();
        let uvs = project_along(&ts);
        let marks: Vec<(usize, usize)> = uvs.iter().map(|p| (span(&uk, p[0]), span(&vk, p[1]))).collect();
        // (The end samples, on a closed surface's seam, are left out.)
        for i in 1..N - 1 {
            if marks[i] == marks[i + 1] {
                continue;
            }
            let (mut lo, mut hi) = (ts[i], ts[i + 1]);
            let mut seed = (uvs[i][0], uvs[i][1]);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                let (u, v, _) = super::extrema::nearest_on_surface(curve.point(mid), surface, Some(seed));
                if (span(&uk, u), span(&vk, v)) == marks[i] {
                    lo = mid;
                    seed = (u, v);
                } else {
                    hi = mid;
                }
            }
            cuts.push(0.5 * (lo + hi));
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|x, y| (*x - *y).abs() <= 1e-9 * (b - a));
    }
    let degree = 3;
    let mut spans_per_piece = 2;
    while spans_per_piece * (cuts.len() - 1) <= 512 {
        let mut flat = vec![a; degree + 1];
        for (w, piece) in cuts.windows(2).zip(0..) {
            if piece > 0 {
                flat.extend(std::iter::repeat(w[0]).take(degree));
            }
            for i in 1..spans_per_piece {
                flat.push(w[0] + (w[1] - w[0]) * i as f64 / spans_per_piece as f64);
            }
        }
        flat.extend(std::iter::repeat(b).take(degree + 1));
        let mut spans = flat.clone();
        spans.dedup();
        let mut params = Vec::new();
        for w in spans.windows(2) {
            for i in 0..12 {
                params.push(w[0] + (w[1] - w[0]) * i as f64 / 12.0);
            }
        }
        params.push(b);
        let uv = project_along(&params);
        let fitted = crate::geometry::fit::least_squares_with_knots::<Point2>(&uv, &params, degree, flat)?;
        let mut worst = 0.0f64;
        for w in spans.windows(2) {
            for i in 0..=24 {
                let t = w[0] + (w[1] - w[0]) * i as f64 / 24.0;
                let q = fitted.point(t);
                // (Against the curve's point, or projecting, its foot on the surface.)
                let target = if projecting {
                    let (u, v, _) = super::extrema::nearest_on_surface(curve.point(t), surface, Some((q.x, q.y)));
                    surface.point(u, v)
                } else {
                    curve.point(t)
                };
                worst = worst.max(surface.point(q.x, q.y).distance(target));
            }
        }
        if worst <= tolerance.max(1e-9) {
            return Some(fitted);
        }
        spans_per_piece *= 2;
    }
    None
}
