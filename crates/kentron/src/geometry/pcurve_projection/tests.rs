//! Pcurves against closed forms: exact lines where the feet move affinely, accuracy within
//! ten times the tolerance elsewhere.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, TAU};

use super::*;
use crate::geometry::bspline::{BSplineSurface, KnotVector};
use crate::geometry::curve::{Circle, Ellipse, Line};
use crate::geometry::surface::{Cone, Cylinder, Sphere, Torus};
use crate::math::{Direction, Frame};

/// The nine check parameters of `[first, last]`.
fn checks(first: f64, last: f64) -> impl Iterator<Item = f64> {
    (0..=8).map(move |i| first + (last - first) * i as f64 / 8.0)
}

fn project(surface: &Surface, curve: &Curve, first: f64, last: f64, preci: f64) -> Curve2 {
    PCurveProjector::new(surface, preci).perform(curve, first, last, preci, preci).expect("a pcurve")
}

/// The pcurve is a line through `expected(t)` at the check points (to 1e-9).
fn assert_line(pc: &Curve2, first: f64, last: f64, expected: impl Fn(f64) -> (f64, f64)) {
    assert!(matches!(pc, Curve2::Line(_)), "not a line: {pc:?}");
    for t in checks(first, last) {
        let (q, (u, v)) = (pc.point(t), expected(t));
        assert!((q.x - u).abs() <= 1e-9 && (q.y - v).abs() <= 1e-9, "t = {t}: {q:?} vs ({u}, {v})");
    }
}

/// The pcurve is a B-spline within `tol` (in each parameter) of `expected(t)` at the check
/// points and densely between.
fn assert_close(pc: &Curve2, first: f64, last: f64, tol: f64, expected: impl Fn(f64) -> (f64, f64)) {
    assert!(matches!(pc, Curve2::BSpline(_)), "not a B-spline: {pc:?}");
    for i in 0..=400 {
        let t = first + (last - first) * i as f64 / 400.0;
        let (q, (u, v)) = (pc.point(t), expected(t));
        assert!((q.x - u).abs() <= tol && (q.y - v).abs() <= tol, "t = {t}: {q:?} vs ({u}, {v})");
    }
}

fn cylinder() -> Surface {
    Surface::Cylinder(Cylinder::new(Frame::WORLD, 2.0))
}

fn frame(origin: (f64, f64, f64), z: (f64, f64, f64), x: (f64, f64, f64)) -> Frame {
    Frame::new(Point::new(origin.0, origin.1, origin.2), Direction::new(z.0, z.1, z.2), Direction::new(x.0, x.1, x.2))
}

#[test]
fn circle_on_cylinder() {
    let circle = Curve::Circle(Circle::new(frame((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 2.0));
    let pc = project(&cylinder(), &circle, 0.0, TAU, 1e-7);
    assert_line(&pc, 0.0, TAU, |t| (t, 1.0));
}

/// The oblique section of the cylinder by a plane tilted 0.4 about the x axis: at `t` the
/// point of angle `π/2 + t`, height `2 tan(0.4) cos t`.
fn oblique_ellipse() -> Curve {
    let (s, c) = 0.4f64.sin_cos();
    Curve::Ellipse(Ellipse::new(frame((0.0, 0.0, 0.0), (0.0, -s, c), (0.0, c, s)), 2.0 / c, 2.0))
}

#[test]
fn ellipse_on_cylinder() {
    let k = 2.0 * 0.4f64.tan();
    for (first, last) in [(0.0, TAU), (1.0, 4.0)] {
        for preci in [1e-7, 1e-5] {
            let pc = project(&cylinder(), &oblique_ellipse(), first, last, preci);
            assert_close(&pc, first, last, 10.0 * preci, |t| (FRAC_PI_2 + t, k * t.cos()));
        }
    }
}

#[test]
fn arc_across_seam() {
    let arc = Curve::Circle(Circle::new(frame((0.0, 0.0, 0.5), (0.0, 0.0, 1.0), (5f64.cos(), 5f64.sin(), 0.0)), 2.0));
    let pc = project(&cylinder(), &arc, 0.0, 3.0, 1e-7);
    assert_line(&pc, 0.0, 3.0, |t| (5.0 + t, 0.5));
}

#[test]
fn circle_backwards_from_the_seam() {
    let circle = Curve::Circle(Circle::new(frame((0.0, 0.0, -0.5), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0)), 2.0));
    let pc = project(&cylinder(), &circle, 0.0, 4.0, 1e-7);
    assert_line(&pc, 0.0, 4.0, |t| (TAU - t, -0.5));
}

fn sphere() -> Surface {
    Surface::Sphere(Sphere::new(Frame::WORLD, 3.0))
}

#[test]
fn latitude_on_sphere() {
    let circle = Curve::Circle(Circle::new(frame((0.0, 0.0, 1.5), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 3.0 * FRAC_PI_6.cos()));
    let pc = project(&sphere(), &circle, 0.0, TAU, 1e-7);
    assert_line(&pc, 0.0, TAU, |t| (t, FRAC_PI_6));
}

/// The great circle through the poles at longitude `u`, starting on the equator.
fn meridian(u: f64) -> Curve {
    Curve::Circle(Circle::new(frame((0.0, 0.0, 0.0), (u.sin(), -u.cos(), 0.0), (u.cos(), u.sin(), 0.0)), 3.0))
}

#[test]
fn meridian_over_pole() {
    let s = sphere();
    let curve = meridian(0.0);
    let pc = project(&s, &curve, 0.5, 2.5, 1e-7);
    for t in checks(0.5, 2.5) {
        let q = pc.point(t);
        assert!(q.y.abs() <= FRAC_PI_2, "t = {t}: {q:?}");
        if (t - FRAC_PI_2).abs() >= 0.3 {
            assert!(s.point(q.x, q.y).distance(curve.point(t)) <= 3e-3, "t = {t}: {q:?}");
        }
    }
}

#[test]
fn meridian_to_pole() {
    let pc = project(&sphere(), &meridian(0.0), 0.0, FRAC_PI_2, 1e-7);
    assert_line(&pc, 0.0, FRAC_PI_2, |t| (0.0, t));
    // At another longitude the pole's longitude follows the meridian's.
    let pc = project(&sphere(), &meridian(1.0), 0.0, FRAC_PI_2, 1e-7);
    assert_line(&pc, 0.0, FRAC_PI_2, |t| (1.0, t));
}

#[test]
fn meridian_on_torus() {
    let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0));
    let tube = Curve::Circle(Circle::new(frame((3.0, 0.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0)), 1.0));
    let pc = project(&torus, &tube, 0.0, TAU, 1e-7);
    assert_line(&pc, 0.0, TAU, |t| (0.0, t));
}

#[test]
fn circle_on_cone() {
    let (alpha, r) = (0.3, 1.5);
    let cone = Surface::Cone(Cone::new(Frame::WORLD, alpha, r));
    let circle = Curve::Circle(Circle::new(frame((0.0, 0.0, 2.0 * alpha.cos()), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), r + 2.0 * alpha.sin()));
    let pc = project(&cone, &circle, 0.0, TAU, 1e-7);
    assert_line(&pc, 0.0, TAU, |t| (t, 2.0));
}

#[test]
fn generator_into_the_apex() {
    let cone_s = Cone::new(Frame::WORLD, 0.3, 1.5);
    let cone = Surface::Cone(cone_s);
    let v_apex = -1.5 / 0.3f64.sin();
    let start = cone.point(0.7, 1.0);
    let length = 1.0 - v_apex;
    let line = Curve::Line(Line::new(start, Direction::from_vector(cone_s.apex() - start)));
    let pc = project(&cone, &line, 0.0, length, 1e-7);
    assert_line(&pc, 0.0, length, |t| (0.7, 1.0 - t));
}

/// A bicubic patch bent both ways.
fn patch() -> Surface {
    let poles: Vec<Vec<Point>> = (0..5)
        .map(|i| (0..4).map(|j| Point::new(i as f64, j as f64 * 1.2, 0.3 * ((i as f64 - 2.0).powi(2) - (j as f64 - 1.5).powi(2)) + 0.1 * (i * j) as f64)).collect())
        .collect();
    let u = KnotVector::new(vec![0.0, 0.4, 1.0], vec![4, 1, 4], 3);
    let v = KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
    Surface::BSpline(BSplineSurface::new(poles, u, v))
}

/// A cubic curve whose control points lie on the patch (so the curve lies slightly off it).
fn curve_near(surface: &Surface) -> Curve {
    let uv = [(0.05, 0.1), (0.3, 0.5), (0.5, 0.2), (0.7, 0.8), (0.95, 0.6)];
    let poles: Vec<Point> = uv.iter().map(|&(u, v)| surface.point(u, v)).collect();
    Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 0.5, 1.0], vec![4, 1, 4], 3))
}

#[test]
fn curve_near_bspline() {
    let s = patch();
    let curve = curve_near(&s);
    for (first, last) in [(0.0, 1.0), (0.2, 0.9)] {
        let pc = project(&s, &curve, first, last, 1e-4);
        assert!(matches!(pc, Curve2::BSpline(_)));
        let mut start = None;
        for t in checks(first, last) {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(curve.point(t), &s, start);
            start = Some((u, v));
            let q = pc.point(t);
            assert!(s.point(q.x, q.y).distance(s.point(u, v)) <= 1e-4, "t = {t}");
        }
    }
}

#[test]
fn isos_of_a_bspline() {
    let s = patch();
    let iso = s.u_iso(0.0);
    let pc = project(&s, &iso, 0.0, 1.0, 1e-7);
    assert_line(&pc, 0.0, 1.0, |t| (0.0, t));
    let Curve2::Line(line) = &pc else { unreachable!() };
    assert_eq!(line.origin.x, 0.0);
    let iso = s.v_iso(0.37);
    let pc = project(&s, &iso, 0.1, 0.8, 1e-7);
    assert_line(&pc, 0.1, 0.8, |t| (t, 0.37));
}

/// A rational quadratic circle (its speed jumps at its double knots) on a cylinder: accurate
/// across the knots.
#[test]
fn rational_circle_on_cylinder() {
    let circle = Circle::new(frame((0.0, 0.0, 0.25), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 2.0);
    let curve = Curve::BSpline(crate::geometry::convert::curve_to_bspline(&Curve::Circle(circle)).expect("a B-spline circle"));
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    for preci in [1e-7, 1e-4] {
        let pc = project(&cylinder(), &curve, first, last, preci);
        for i in 0..=997 {
            let t = first + (last - first) * i as f64 / 997.0;
            let q = pc.point(t);
            let p = curve.point(t);
            let exact = p.y.atan2(p.x).rem_euclid(TAU);
            let exact = if i == 997 { exact + TAU } else { exact };
            let du = (q.x - exact).abs().min((q.x - exact - TAU).abs());
            assert!(du * 2.0 <= 10.0 * preci && (q.y - 0.25).abs() <= 10.0 * preci, "t = {t}: {q:?} vs {exact}");
        }
    }
}

/// A curve going round a torus, crossing both seams: continuous, and within the tolerance of
/// the feet.
#[test]
fn torus_crossing_both_seams() {
    let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0));
    let n = 60;
    let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let points: Vec<Point> = params.iter().map(|&t| torus.point(5.0 + 2.5 * t, 4.0 + 4.0 * t)).collect();
    let curve = Curve::BSpline(crate::geometry::fit::interpolate_with_parameters(&points, &params).expect("a curve"));
    let preci = 1e-5;
    let pc = project(&torus, &curve, 0.0, 1.0, preci);
    let mut previous: Option<Point2> = None;
    for i in 0..=500 {
        let t = i as f64 / 500.0;
        let q = pc.point(t);
        if let Some(p) = previous {
            assert!((q.x - p.x).abs() < 0.1 && (q.y - p.y).abs() < 0.1, "a jump at {t}");
        }
        previous = Some(q);
        let (u, v) = torus.parameters(curve.point(t)).unwrap();
        assert!(torus.point(q.x, q.y).distance(torus.point(u, v)) <= 10.0 * preci, "t = {t}");
    }
    // (Past both seams: u from 5 beyond 2π, v from 4 beyond 2π.)
    let end = pc.point(1.0);
    assert!((end.x - 7.5).abs() < 1e-3 && (end.y - 8.0).abs() < 1e-3, "{end:?}");
}

/// A line held a millimetre above a flat patch parametrised unevenly: the feet are its
/// shadow, at the gap.
#[test]
fn line_off_a_patch() {
    let xs = [0.0, 0.2, 1.1, 3.0];
    let poles: Vec<Vec<Point>> = (0..4).map(|i| (0..4).map(|j| Point::new(xs[i], xs[j] * 0.5 + 0.1 * xs[i], 0.0)).collect()).collect();
    let kv = || KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
    let s = Surface::BSpline(BSplineSurface::new(poles, kv(), kv()));
    let h = 1e-3;
    let line = Curve::Line(Line::new(Point::new(0.3, 0.2, h), Direction::new(1.0, 0.4, 0.0)));
    let preci = 1e-5;
    let pc = project(&s, &line, 0.0, 2.0, preci);
    for i in 0..=200 {
        let t = 2.0 * i as f64 / 200.0;
        let q = pc.point(t);
        let p = line.point(t);
        let foot = s.point(q.x, q.y);
        assert!((foot.distance(p) - h).abs() <= 10.0 * preci, "t = {t}");
        assert!(foot.distance(Point::new(p.x, p.y, 0.0)) <= 10.0 * preci, "t = {t}");
    }
}

/// A patch whose side u = 0 is one point, and its iso v = 0.6 from that point: the start's
/// v continues the curve's.
#[test]
fn iso_from_a_collapsed_side() {
    let poles: Vec<Vec<Point>> = (0..4)
        .map(|i| (0..4).map(|j| if i == 0 { Point::new(0.0, 0.0, 1.0) } else { Point::new(i as f64, j as f64 - 1.5, 0.2 * j as f64 + 1.0 - 0.3 * i as f64) }).collect())
        .collect();
    let kv = || KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
    let s = Surface::BSpline(BSplineSurface::new(poles, kv(), kv()));
    let iso = s.v_iso(0.6);
    let pc = project(&s, &iso, 0.0, 1.0, 1e-7);
    let start = pc.point(0.0);
    assert!(start.x.abs() <= 1e-9 && (start.y - 0.6).abs() <= 1e-6, "{start:?}");
    for t in checks(0.0, 1.0) {
        let q = pc.point(t);
        assert!(s.point(q.x, q.y).distance(iso.point(t)) <= 1e-6, "t = {t}: {q:?}");
    }
}

#[test]
fn tiny_edge() {
    let s = patch();
    let p = s.point(0.3, 0.4);
    let line = Curve::Line(Line::new(p, Direction::new(1.0, 0.0, 0.0)));
    let pc = project(&s, &line, 0.0, 1e-9, 1e-7);
    let (a, b) = (pc.point(0.0), pc.point(1e-9));
    assert!(a.distance(b) > 0.0);
    assert!(s.point(a.x, a.y).distance(p) <= 1e-7);
    // (Degenerate ranges are refused.)
    assert!(PCurveProjector::new(&s, 1e-7).perform(&line, 1.0, 1.0, 1e-7, 1e-7).is_none());
}

/// A circle once round a cylinder written as a B-spline (closed, not periodic), from its
/// seam: within the bounds, accurate.
#[test]
fn circle_on_a_bspline_cylinder() {
    let full = Surface::Trimmed(crate::geometry::surface::TrimmedSurface::new(cylinder(), 0.0, TAU, -1.0, 1.0));
    let s = Surface::BSpline(crate::geometry::convert::surface_to_bspline(&full).expect("a tube"));
    let circle = Curve::Circle(Circle::new(frame((0.0, 0.0, 0.3), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 2.0));
    let preci = 1e-6;
    let pc = project(&s, &circle, 0.0, TAU, preci);
    let (u0, u1, _, _) = s.bounds();
    for i in 0..=300 {
        let t = TAU * i as f64 / 300.0;
        let q = pc.point(t);
        assert!(q.x >= u0 - 1e-9 && q.x <= u1 + 1e-9, "t = {t}: {q:?}");
        assert!(s.point(q.x, q.y).distance(circle.point(t)) <= 10.0 * preci, "t = {t}");
    }
}
