//! Geom2d_BSplineCurve_Test.cxx.
//!
//! Not ported: copy tests and `WeightsArray_*` (see geom2d_bezier_curve).

use kentron::geometry::curve2::BSplineCurve2;
use kentron::geometry::ParametricCurve;
use kentron::math::{Point2, Transform2, Transformable2, Vector2};

use crate::{assert_point2_near, p2};

/// The cubic Bezier-like B-spline of the OCCT fixture.
fn cubic() -> BSplineCurve2 {
    BSplineCurve2::new(
        vec![p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 1.0), p2(3.0, 0.0)],
        vec![0.0, 1.0],
        vec![4, 4],
        3,
    )
}

fn rational_quadratic(poles: [Point2; 3], weights: [f64; 3]) -> BSplineCurve2 {
    BSplineCurve2::new(poles.to_vec(), vec![0.0, 1.0], vec![3, 3], 2).with_weights(weights.to_vec())
}

fn pentagon() -> BSplineCurve2 {
    let poles = vec![p2(1.0, 0.0), p2(0.309, 0.951), p2(-0.809, 0.588), p2(-0.809, -0.588), p2(0.309, -0.951)];
    BSplineCurve2::new_periodic(poles, (0..6).map(|i| i as f64 * 0.2).collect(), vec![1; 6], 3)
}

#[test]
fn evaluation() {
    let c = cubic();
    assert_point2_near(c.point(0.0), p2(0.0, 0.0), 1e-10);
    assert_point2_near(c.point(1.0), p2(3.0, 0.0), 1e-10);
    assert!(c.d1(0.5).1.length() > 0.0);
    assert!(c.d2(0.5).1.length() > 0.0);
    assert!(c.dn(0.5, 1).length() > 0.0);
    assert!(c.dn(0.5, 4).length() <= 1e-10);
    assert_point2_near(c.start_point(), p2(0.0, 0.0), 1e-10);
    assert_point2_near(c.end_point(), p2(3.0, 0.0), 1e-10);
}

#[test]
fn properties() {
    let c = cubic();
    assert_eq!(c.degree(), 3);
    assert!(!c.is_periodic() && !c.is_rational() && !c.is_closed());
    assert!(c.is_cn(3));
    assert_eq!(c.weights(), None);
    assert_eq!((c.knot(0), c.knot(1)), (0.0, 1.0));
    assert_eq!(c.knots().len(), 2);
    assert_eq!(c.multiplicities(), &[4, 4]);
}

#[test]
fn set_pole_and_weight() {
    let mut c = cubic();
    c.set_pole(1, p2(1.0, 5.0));
    assert_eq!(c.pole(1), p2(1.0, 5.0));
    assert_point2_near(c.start_point(), p2(0.0, 0.0), 1e-10);
    let mut r = rational_quadratic([p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 0.0)], [1.0, 1.0, 1.0]);
    let before = r.point(0.5);
    r.set_weight(1, 5.0);
    assert_eq!(r.weight(1), 5.0);
    assert!(r.is_rational());
    assert!(r.point(0.5).distance(before) > 1e-6);
}

#[test]
fn knot_operations_keep_shape() {
    let mut c = cubic();
    let before = c.point(0.5);
    c.increase_degree(5);
    assert_eq!(c.degree(), 5);
    assert_point2_near(c.point(0.5), before, 1e-10);

    let mut c = cubic();
    c.insert_knot(0.5, 1);
    assert_eq!(c.knots().len(), 3);
    assert_point2_near(c.point(0.5), before, 1e-10);
    let quarter = c.point(0.25);
    assert!(c.remove_knot(1, 0, 1e-6));
    assert_eq!(c.knots().len(), 2);
    assert_point2_near(c.point(0.25), quarter, 1e-6);

    let mut c = cubic();
    c.insert_knots(&[(0.25, 1), (0.75, 1)]);
    assert_eq!(c.knots().len(), 4);
    assert_point2_near(c.point(0.5), before, 1e-10);
}

#[test]
fn segment_and_reverse() {
    let mut c = cubic();
    let (a, b) = (c.point(0.25), c.point(0.75));
    c.segment(0.25, 0.75);
    assert_point2_near(c.start_point(), a, 1e-6);
    assert_point2_near(c.end_point(), b, 1e-6);
    let mut c = cubic();
    let (s, e) = (c.start_point(), c.end_point());
    c.reverse();
    assert_point2_near(c.start_point(), e, 1e-10);
    assert_point2_near(c.end_point(), s, 1e-10);
    let mut r = rational_quadratic([p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)], [1.0, 0.5_f64.sqrt(), 1.0]);
    let mid = r.point(0.5);
    r.segment(0.25, 0.75);
    assert_point2_near(r.point(0.5), mid, 1e-6);
    assert!(r.is_rational());
}

#[test]
fn resolution_and_transform() {
    assert!(cubic().resolution(1.0) > 0.0);
    let c = cubic();
    let moved = c.transformed(&Transform2::translation(Vector2::new(10.0, 20.0)));
    let (a, b) = (c.point(0.5), moved.point(0.5));
    assert!((b.x - a.x - 10.0).abs() <= 1e-10 && (b.y - a.y - 20.0).abs() <= 1e-10);
}

#[test]
fn periodic_curve() {
    let mut c = pentagon();
    assert!(c.is_periodic());
    let before = c.point(0.5);
    c.set_not_periodic();
    assert!(!c.is_periodic());
    assert_point2_near(c.point(0.5), before, 1e-10);
}

#[test]
fn move_point() {
    let mut c = cubic();
    let target = p2(1.5, 0.8);
    let last = c.poles().len() - 1;
    let modified = c.move_point(0.5, target, 0, last).expect("move failed");
    assert!(modified.0 <= modified.1);
    assert_point2_near(c.point(0.5), target, 1e-6);
}

#[test]
fn locate() {
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    let (i1, i2) = c.locate(0.25, 1e-10);
    assert!(i2 < c.knots().len());
    assert!(c.knot(i1) <= 0.25 + 1e-10 && c.knot(i2) >= 0.25 - 1e-10);
}

#[test]
fn local_d1_matches_global() {
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    let (p, d) = c.d1(0.25);
    let (pl, dl) = c.local_d1(0.25, 0, 1);
    assert_point2_near(p, pl, 1e-10);
    assert!((d.x - dl.x).abs() <= 1e-10 && (d.y - dl.y).abs() <= 1e-10);
}

#[test]
fn closed_curve() {
    let c = BSplineCurve2::new(vec![p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 0.0), p2(0.0, 0.0)], vec![0.0, 1.0], vec![4, 4], 3);
    assert!(c.is_closed());
}

#[test]
fn invalid_data_is_rejected() {
    // 4 poles need multiplicities summing to 4 + degree + 1.
    assert!(BSplineCurve2::try_new(vec![p2(0.0, 0.0); 4], vec![0.0, 1.0], vec![3, 3], 3).is_err());
    assert!(BSplineCurve2::try_new(vec![p2(0.0, 0.0); 4], vec![1.0, 0.0], vec![4, 4], 3).is_err());
}
