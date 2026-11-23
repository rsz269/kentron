//! Curve–plane meetings against closed forms and against a dense independent sampling of
//! the distance to the plane.

use std::f64::consts::PI;

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Ellipse, Hyperbola, Line, Parabola};
use crate::geometry::ParametricCurve;

/// The plane through `origin` with normal `normal` (not necessarily unit).
pub(crate) fn plane(origin: Point, normal: Vector) -> Plane {
    plane_through(origin, Direction::from_vector(normal))
}

/// A cubic B-spline over `[0, 6]` waving up and down in z along x (single knots: C2).
pub(crate) fn wavy_bspline() -> Curve {
    let poles: Vec<Point> = (0..9).map(|i| Point::new(i as f64 * 0.75, 0.1 * i as f64, if i % 2 == 0 { -1.0 } else { 1.0 })).collect();
    Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], vec![4, 1, 1, 1, 1, 1, 4], 3))
}

/// A quadratic B-spline winding about the z axis while rising (a rough helix, two turns).
pub(crate) fn helix_bspline() -> Curve {
    let n = 17;
    let poles: Vec<Point> = (0..n)
        .map(|i| {
            let a = i as f64 * PI / 4.0;
            Point::new(a.cos(), a.sin(), 0.15 * i as f64)
        })
        .collect();
    let knots: Vec<f64> = (0..=(n - 2)).map(|i| i as f64).collect();
    let mut mults = vec![1; knots.len()];
    mults[0] = 3;
    *mults.last_mut().unwrap() = 3;
    Curve::BSpline(BSplineCurve::new(poles, knots, mults, 2))
}

/// A quadratic arc `x = w`, `y = 0`, `z = (w − ½)²` over `[0, 1]`, touching z = 0 at ½.
pub(crate) fn tangent_bspline() -> Curve {
    let poles = vec![Point::new(0.0, 0.0, 0.25), Point::new(0.5, 0.0, -0.25), Point::new(1.0, 0.0, 0.25)];
    Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 1.0], vec![3, 3], 2))
}

/// The view of `curve` over its whole range.
pub(crate) fn whole(curve: &Curve) -> GeomAdaptor {
    GeomAdaptor::new(curve, curve.first_parameter(), curve.last_parameter())
}

/// The crossings of the signed distance sampled at `n` steps (sign changes only).
fn sampled_crossings(a: &dyn CurveAdaptor, pl: &Plane, n: usize) -> usize {
    let pd = PlaneData::new(pl);
    let (w0, w1) = (a.first_parameter(), a.last_parameter());
    let mut count = 0;
    let mut prev = pd.distance(a.value(w0));
    for i in 1..=n {
        let f = pd.distance(a.value(w0 + (w1 - w0) * i as f64 / n as f64));
        if (prev < 0.0) != (f < 0.0) {
            count += 1;
        }
        prev = f;
    }
    count
}

/// Every reported point lies on the plane and on the curve at its parameter; points ascend;
/// `(u, v)` are the point's plane coordinates; transitions follow the sign of `N·C′`.
fn check_points(a: &dyn CurveAdaptor, pl: &Plane, found: &CurvePlaneIntersection) {
    let pd = PlaneData::new(pl);
    for w in found.points.windows(2) {
        assert!(w[0].w <= w[1].w, "not ascending: {} {}", w[0].w, w[1].w);
    }
    for p in &found.points {
        assert!(pd.distance(p.point).abs() < 1e-11, "distance {} at {}", pd.distance(p.point), p.w);
        assert!(p.point.distance(a.value(p.w)) < 1e-15 + 1e-15 * p.point.distance(Point::ORIGIN));
        let rel = p.point - pd.origin;
        assert!((rel.dot(pd.x) - p.u).abs() < 1e-11 && (rel.dot(pd.y) - p.v).abs() < 1e-11);
        let d = a.d1(p.w).1;
        let c = pd.n.dot(d) / d.length();
        let expected = if c < -1e-12 { Transition::In } else if c > 1e-12 { Transition::Out } else { Transition::Tangent };
        assert_eq!(p.transition, expected);
    }
}

#[test]
fn plane_axes_follow_the_default_frame() {
    for n in [Vector::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 0.0), Vector::new(0.3, -2.0, 0.7), Vector::new(-1.0, 1.0, 1.0)] {
        let d = Direction::from_vector(n);
        let p = plane_through(Point::new(1.0, 2.0, 3.0), d);
        let f = Frame::from_z(Point::new(1.0, 2.0, 3.0), d);
        assert!(p.frame().x().approx_eq(f.x(), 1e-15) && p.frame().y().approx_eq(f.y(), 1e-15));
        assert_eq!(p.frame().z(), d);
    }
    let p = plane_through(Point::ORIGIN, Direction::Z);
    assert_eq!(p.frame().x(), Direction::X);
    assert_eq!(p.frame().y(), Direction::Y);
}

#[test]
fn lines() {
    let line = Curve::Line(Line::new(Point::new(0.0, 0.0, -2.0), Direction::new(0.0, 0.6, 0.8)));
    let pl = plane(Point::new(1.0, 1.0, 1.0), Vector::new(0.0, 0.0, 1.0));
    let found = intersect_curve_plane(&GeomAdaptor::new(&line, -10.0, 10.0), &pl).unwrap();
    assert!(!found.parallel);
    assert_eq!(found.points.len(), 1);
    let p = found.points[0];
    assert!((p.w - 3.0 / 0.8).abs() < 1e-14);
    assert!((p.point.z - 1.0).abs() < 1e-14 && (p.point.y - 0.6 * 3.75).abs() < 1e-14);
    assert!((p.u - -1.0).abs() < 1e-14 && (p.v - (0.6 * 3.75 - 1.0)).abs() < 1e-14);
    assert_eq!(p.transition, Transition::Out);
    // Beyond a short range: nothing.
    let found = intersect_curve_plane(&GeomAdaptor::new(&line, 0.0, 3.0), &pl).unwrap();
    assert!(found.points.is_empty() && !found.parallel);
    // Against the normal: In.
    let found = intersect_curve_plane(&GeomAdaptor::new(&line, -10.0, 10.0), &plane(Point::new(0.0, 0.0, 1.0), Vector::new(0.0, 0.0, -1.0))).unwrap();
    assert_eq!(found.points[0].transition, Transition::In);
    // Parallel, off the plane or in it.
    let flat = Curve::Line(Line::new(Point::new(0.0, 0.0, 1.0), Direction::new(1.0, 1.0, 0.0)));
    for z in [1.0, 2.0] {
        let found = intersect_curve_plane(&GeomAdaptor::new(&flat, -10.0, 10.0), &plane(Point::new(0.0, 0.0, z), Vector::new(0.0, 0.0, 1.0))).unwrap();
        assert!(found.parallel && found.points.is_empty());
    }
}

#[test]
fn circles() {
    let r = 2.0;
    let circle = Curve::Circle(Circle::new(Frame::WORLD, r));
    let full = GeomAdaptor::new(&circle, 0.0, 2.0 * PI);
    // The plane x = d (normal in the circle's plane): w = ±acos(d/r).
    let d = 0.7;
    let pl = plane(Point::new(d, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    let found = intersect_curve_plane(&full, &pl).unwrap();
    check_points(&full, &pl, &found);
    let alpha = (d / r).acos();
    assert_eq!(found.points.len(), 2);
    assert!((found.points[0].w - alpha).abs() < 1e-14);
    assert!((found.points[1].w - (2.0 * PI - alpha)).abs() < 1e-14);
    assert_eq!(found.points[0].transition, Transition::In);
    assert_eq!(found.points[1].transition, Transition::Out);
    // Trimmed to one of them, and far away.
    let found = intersect_curve_plane(&GeomAdaptor::new(&circle, 0.0, PI), &pl).unwrap();
    assert_eq!(found.points.len(), 1);
    let found = intersect_curve_plane(&full, &plane(Point::new(3.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))).unwrap();
    assert!(found.points.is_empty() && !found.parallel);
    // In the plane, or parallel to it.
    for z in [0.0, 0.5] {
        let found = intersect_curve_plane(&full, &plane(Point::new(0.0, 0.0, z), Vector::new(0.0, 0.0, 1.0))).unwrap();
        assert!(found.parallel && found.points.is_empty());
    }
    // Tangent: one point.
    let found = intersect_curve_plane(&full, &plane(Point::new(0.0, r, 0.0), Vector::new(0.0, 1.0, 0.0))).unwrap();
    assert_eq!(found.points.len(), 1);
    assert!((found.points[0].w - PI / 2.0).abs() < 1e-14);
    assert_eq!(found.points[0].transition, Transition::Tangent);
    // A range starting mid-period: parameters in [w0, w0 + 2π).
    let shifted = GeomAdaptor::new(&circle, 4.0, 4.0 + 2.0 * PI);
    let found = intersect_curve_plane(&shifted, &pl).unwrap();
    assert_eq!(found.points.len(), 2);
    assert!((found.points[0].w - (2.0 * PI - alpha)).abs() < 1e-14);
    assert!((found.points[1].w - (2.0 * PI + alpha)).abs() < 1e-13);
}

#[test]
fn ellipses_parabolas_hyperbolas() {
    let frame = Frame::new(Point::new(0.5, -0.25, 0.1), Direction::new(0.1, 0.2, 1.0), Direction::new(1.0, 0.3, 0.0));
    let cases: Vec<(GeomAdaptor, usize)> = vec![
        (GeomAdaptor::new(&Curve::Ellipse(Ellipse::new(frame, 3.0, 1.0)), 0.0, 2.0 * PI), 2),
        (GeomAdaptor::new(&Curve::Parabola(Parabola::new(frame, 0.5)), -5.0, 5.0), 2),
        (GeomAdaptor::new(&Curve::Hyperbola(Hyperbola::new(frame, 1.0, 0.5)), -3.0, 3.0), 2),
    ];
    let planes = [plane(Point::new(1.5, 0.0, 0.0), Vector::new(1.0, 0.1, 0.0)), plane(Point::new(1.5, 0.3, 0.2), Vector::new(1.0, -0.4, 2.0))];
    for (a, expected) in &cases {
        for pl in &planes {
            let found = intersect_curve_plane(a, pl).unwrap();
            check_points(a, pl, &found);
            assert!(!found.parallel);
            assert_eq!(found.points.len(), sampled_crossings(a, pl, 20000), "{:?}", a.kind());
            let _ = expected;
        }
    }
    // A parabola and an ellipse in the cutting plane: parallel; an ellipse off it: nothing.
    let flat = Frame::WORLD;
    let par = GeomAdaptor::new(&Curve::Parabola(Parabola::new(flat, 1.0)), -2.0, 2.0);
    let found = intersect_curve_plane(&par, &plane(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0))).unwrap();
    assert!(found.parallel && found.points.is_empty());
    let ell = GeomAdaptor::new(&Curve::Ellipse(Ellipse::new(flat, 2.0, 1.0)), 0.0, 2.0 * PI);
    let found = intersect_curve_plane(&ell, &plane(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0))).unwrap();
    assert!(found.parallel && found.points.is_empty());
    let found = intersect_curve_plane(&ell, &plane(Point::new(0.0, 0.0, 1.0), Vector::new(0.0, 0.0, 1.0))).unwrap();
    assert!(!found.parallel && found.points.is_empty());
    // Closed forms: the parabola x = t²/4, y = t against x = 1 crosses at t = ±2, ascending;
    // the hyperbola (cosh t, sinh t) against x = cosh 1 at t = ±1.
    let found = intersect_curve_plane(&par, &plane(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))).unwrap();
    assert_eq!(found.points.len(), 2);
    assert!((found.points[0].w + 2.0).abs() < 1e-14 && (found.points[1].w - 2.0).abs() < 1e-14);
    let hyp = GeomAdaptor::new(&Curve::Hyperbola(Hyperbola::new(flat, 1.0, 1.0)), -3.0, 3.0);
    let found = intersect_curve_plane(&hyp, &plane(Point::new(1f64.cosh(), 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))).unwrap();
    assert_eq!(found.points.len(), 2);
    assert!((found.points[0].w + 1.0).abs() < 1e-14 && (found.points[1].w - 1.0).abs() < 1e-14);
    assert_eq!(found.points[0].transition, Transition::In);
    assert_eq!(found.points[1].transition, Transition::Out);
}

#[test]
fn bsplines() {
    let wavy = wavy_bspline();
    let a = whole(&wavy);
    for pl in [
        plane(Point::new(0.0, 0.0, 0.1), Vector::new(0.0, 0.0, 1.0)),
        plane(Point::new(2.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.2)),
        plane(Point::new(0.0, 0.0, 0.0), Vector::new(0.2, 0.0, 1.0)),
        plane(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0)),
    ] {
        let found = intersect_curve_plane(&a, &pl).unwrap();
        check_points(&a, &pl, &found);
        assert_eq!(found.points.len(), sampled_crossings(&a, &pl, 20000));
        assert!(!found.parallel);
        // Transitions alternate.
        for w in found.points.windows(2) {
            assert_ne!(w[0].transition, w[1].transition);
        }
    }
    // A z = 0.1 plane crosses the wave five or more times.
    let found = intersect_curve_plane(&a, &plane(Point::new(0.0, 0.0, 0.1), Vector::new(0.0, 0.0, 1.0))).unwrap();
    assert!(found.points.len() >= 5);
    // A range beyond the curve's own (the B-spline extended): crossings out there are found.
    let ext = GeomAdaptor::new(&wavy, -0.5, 6.5);
    let pl = plane(Point::new(-0.2, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    let found = intersect_curve_plane(&ext, &pl).unwrap();
    check_points(&ext, &pl, &found);
    assert_eq!(found.points.len(), 1);
    assert!(found.points[0].w < 0.0);
}

#[test]
fn helix_and_periodic_bsplines() {
    let helix = helix_bspline();
    let a = whole(&helix);
    for pl in [
        plane(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        plane(Point::new(0.0, 0.0, 1.2), Vector::new(0.0, 0.0, 1.0)),
        plane(Point::new(0.2, 0.1, 0.0), Vector::new(1.0, 1.0, 0.3)),
    ] {
        let found = intersect_curve_plane(&a, &pl).unwrap();
        check_points(&a, &pl, &found);
        assert_eq!(found.points.len(), sampled_crossings(&a, &pl, 20000));
    }
    // A closed periodic B-spline over a range starting mid-period.
    let poles: Vec<Point> = (0..6).map(|i| {
        let t = i as f64 * PI / 3.0;
        Point::new(t.cos(), t.sin(), 0.3 * (2.0 * t).sin())
    }).collect();
    let periodic = Curve::BSpline(BSplineCurve::new_periodic(poles, (0..=6).map(|i| i as f64).collect(), vec![1; 7], 3));
    let shifted = GeomAdaptor::new(&periodic, 2.5, 8.5);
    let pl = plane(Point::new(0.1, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    let found = intersect_curve_plane(&shifted, &pl).unwrap();
    check_points(&shifted, &pl, &found);
    assert_eq!(found.points.len(), 2);
    assert!(found.points.iter().all(|p| p.w >= 2.5 && p.w < 8.5));
}

#[test]
fn tangent_and_lying_in_the_plane() {
    let t = tangent_bspline();
    let a = whole(&t);
    let pl = plane(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0));
    let found = intersect_curve_plane(&a, &pl).unwrap();
    assert_eq!(found.points.len(), 1);
    assert!((found.points[0].w - 0.5).abs() < 1e-12);
    assert_eq!(found.points[0].transition, Transition::Tangent);
    // Touching at an interior sample point and at a parameter between samples.
    for shift in [0.0f64, 0.0123] {
        let poles = vec![Point::new(0.0, 0.0, (0.5 + shift).powi(2)), Point::new(0.5, 0.0, (0.5 + shift).powi(2) - (0.5 + shift)), Point::new(1.0, 0.0, (0.5 - shift).powi(2))];
        let c = Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 1.0], vec![3, 3], 2));
        let found = intersect_curve_plane(&whole(&c), &pl).unwrap();
        assert_eq!(found.points.len(), 1, "{shift}");
        assert!((found.points[0].w - (0.5 + shift)).abs() < 1e-12);
    }
    // A flat curve in the plane: nothing.
    let flat = Curve::BSpline(BSplineCurve::new(vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 1.0, 0.0)], vec![0.0, 1.0], vec![4, 4], 3));
    let found = intersect_curve_plane(&whole(&flat), &pl).unwrap();
    assert!(found.points.is_empty() && !found.parallel);
    // In the plane on [1, 2] only (degree 1, C0 knots): points only there, tangent.
    let partly = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, 1.0), Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 0.0, 1.0)],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![2, 1, 1, 2],
        1,
    ));
    let a = whole(&partly);
    let found = intersect_curve_plane(&a, &pl).unwrap();
    assert!(!found.points.is_empty());
    for p in &found.points {
        assert!(p.w >= 1.0 && p.w <= 2.0 && p.point.z.abs() <= 1e-8, "{p:?}");
    }
    // A Bezier crossing four times.
    // z = 100 (t − 0.1)(t − 0.4)(t − 0.6)(t − 0.9) in Bernstein form.
    let bez = Curve::Bezier(BezierCurve::new(vec![
        Point::new(0.0, 0.0, 2.16),
        Point::new(1.0, 0.0, -6.09),
        Point::new(2.0, 0.0, 7.826666666666666),
        Point::new(3.0, 0.0, -6.09),
        Point::new(4.0, 0.0, 2.16),
    ]));
    let a = whole(&bez);
    let found = intersect_curve_plane(&a, &pl).unwrap();
    check_points(&a, &pl, &found);
    assert_eq!(found.points.len(), sampled_crossings(&a, &pl, 20000));
    assert_eq!(found.points.len(), 4);
    for (p, w) in found.points.iter().zip([0.1, 0.4, 0.6, 0.9]) {
        assert!((p.w - w).abs() < 1e-12, "{} vs {w}", p.w);
    }
}

/// A crossing exactly at an interior knot (C1 there: degree 2, double knot) is reported once.
#[test]
fn a_root_at_a_knot_is_reported_once() {
    let c = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, -1.0), Point::new(1.0, 0.0, -0.5), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.5), Point::new(4.0, 0.0, 2.0)],
        vec![0.0, 1.0, 2.0],
        vec![3, 2, 3],
        2,
    ));
    let a = whole(&c);
    let found = intersect_curve_plane(&a, &plane(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0))).unwrap();
    assert_eq!(found.points.len(), 1);
    assert!((found.points[0].w - 1.0).abs() < 1e-12);
    assert_eq!(found.points[0].transition, Transition::Out);
}

#[test]
fn records_for_the_line_surface_intersection() {
    let frame = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 1.0, 0.0));
    let ax = Ax3::from_frame(&frame);
    let (m, t) = ax.transformation();
    let p = [2.0, 3.0, 5.0];
    let local: Vec<f64> = (0..3).map(|i| dot3(m[i], p) + t[i]).collect();
    assert!((local[0] - 2f64.sqrt()).abs() < 1e-14 && local[1].abs() < 1e-14 && (local[2] - 2.0).abs() < 1e-14);
    let r = line_plane([0.0, 0.0, 0.0], [0.0, 0.0, 2.0], &ax, TOLERANCE_ANGULAIRE, 0.0, 0.0);
    assert!(r.done && !r.parallel && !r.inquadric);
    assert_eq!(r.points.len(), 1);
    assert!((r.points[0].1 - 1.5).abs() < 1e-15 && (r.points[0].0[2] - 3.0).abs() < 1e-15);
    let r = line_plane([0.0, 0.0, 3.0], [1.0, 0.0, 0.0], &ax, TOLERANCE_ANGULAIRE, 0.0, 0.0);
    assert!(r.done && r.parallel && r.inquadric && r.points.is_empty());
    let r = line_plane([0.0, 0.0, 4.0], [1.0, 0.0, 0.0], &ax, TOLERANCE_ANGULAIRE, 0.0, 0.0);
    assert!(r.done && r.parallel && !r.inquadric);
    let _ = ConicQuad::default();
}
