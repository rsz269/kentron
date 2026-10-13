//! Geom2d_BezierCurve_Test.cxx.
//!
//! Not ported: copy tests, and `WeightsArray_*` (OCCT array ownership; `weights()` returns
//! `None` for polynomial curves and `weight(i)` is 1).

use kentron::geometry::curve2::BezierCurve2;
use kentron::geometry::ParametricCurve;
use kentron::math::{Transform2, Transformable2, Vector2};

use crate::{assert_point2_near, p2};

fn cubic() -> BezierCurve2 {
    BezierCurve2::new(vec![p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 1.0), p2(3.0, 0.0)])
}

fn quadratic(weights: [f64; 3]) -> BezierCurve2 {
    BezierCurve2::rational(vec![p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 0.0)], weights.to_vec())
}

/// A rational quadratic quarter of the unit circle.
fn quarter_circle() -> BezierCurve2 {
    BezierCurve2::rational(vec![p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)], vec![1.0, 0.5_f64.sqrt(), 1.0])
}

#[test]
fn evaluation() {
    let c = cubic();
    assert_point2_near(c.point(0.0), p2(0.0, 0.0), 1e-10);
    assert_point2_near(c.point(1.0), p2(3.0, 0.0), 1e-10);
    assert!(c.d1(0.5).1.length() > 0.0);
    assert!(c.d2(0.5).1.length() > 0.0);
    // The third derivative of a cubic is constant.
    let (a, b) = (c.d3(0.5).3, c.d3(0.25).3);
    assert!((a.x - b.x).abs() <= 1e-8 && (a.y - b.y).abs() <= 1e-8);
    assert!(c.dn(0.5, 1).length() > 0.0);
    assert!(c.dn(0.5, 4).length() <= 1e-10);
}

#[test]
fn end_points_and_properties() {
    let c = cubic();
    assert_point2_near(c.start_point(), p2(0.0, 0.0), 1e-10);
    assert_point2_near(c.end_point(), p2(3.0, 0.0), 1e-10);
    assert_eq!(c.degree(), 3);
    assert!(!c.is_periodic() && !c.is_rational() && !c.is_closed());
    assert!(c.is_cn(100));
    assert_eq!((c.first_parameter(), c.last_parameter()), (0.0, 1.0));
    assert_eq!(c.poles().len(), 4);
    assert_eq!(c.weights(), None);
    assert!(BezierCurve2::MAX_DEGREE >= 25);
}

#[test]
fn set_pole() {
    let mut c = cubic();
    c.set_pole(1, p2(1.0, 5.0));
    assert_eq!(c.pole(1), p2(1.0, 5.0));
    assert_point2_near(c.start_point(), p2(0.0, 0.0), 1e-10);
}

#[test]
fn set_weight_pulls_the_curve() {
    let mut c = quadratic([1.0, 1.0, 1.0]);
    let before = c.point(0.5);
    c.set_weight(1, 10.0);
    assert_eq!(c.weight(1), 10.0);
    assert!(c.is_rational());
    assert!(c.point(0.5).y > before.y);
}

#[test]
fn insert_and_remove_poles() {
    let mut c = cubic();
    c.insert_pole_after(1, p2(1.5, 2.0));
    assert_eq!(c.poles().len(), 5);
    assert_eq!(c.degree(), 4);
    c.remove_pole(2);
    assert_eq!(c.poles().len(), 4);
    let mut r = quadratic([1.0, 2.0, 1.0]);
    r.insert_pole_before_with_weight(1, p2(0.5, 0.5), 3.0);
    assert_eq!(r.poles().len(), 4);
    assert_eq!(r.weight(1), 3.0);
    let mut r = quadratic([1.0, 1.0, 1.0]);
    r.set_pole_with_weight(1, p2(1.0, 2.0), 5.0);
    assert_eq!(r.pole(1), p2(1.0, 2.0));
    assert_eq!(r.weight(1), 5.0);
    assert!(r.is_rational());
}

#[test]
fn increase_degree_keeps_shape() {
    let mut c = cubic();
    let before = c.point(0.5);
    c.increase_degree(5);
    assert_eq!(c.degree(), 5);
    assert_point2_near(c.point(0.5), before, 1e-10);
    let mut r = quadratic([1.0, 2.0, 1.0]);
    let before = r.point(0.5);
    r.increase_degree(5);
    assert_eq!(r.degree(), 5);
    assert!(r.is_rational());
    assert_point2_near(r.point(0.5), before, 1e-10);
}

#[test]
fn segment() {
    let mut c = cubic();
    let (a, b) = (c.point(0.25), c.point(0.75));
    c.segment(0.25, 0.75);
    assert_point2_near(c.start_point(), a, 1e-6);
    assert_point2_near(c.end_point(), b, 1e-6);
    let mut r = quarter_circle();
    let mid = r.point(0.5);
    r.segment(0.25, 0.75);
    assert_point2_near(r.point(0.5), mid, 1e-6);
    assert!(r.is_rational());
}

#[test]
fn reverse() {
    let mut c = cubic();
    let (s, e) = (c.start_point(), c.end_point());
    c.reverse();
    assert_point2_near(c.start_point(), e, 1e-10);
    assert_point2_near(c.end_point(), s, 1e-10);
    let mut r = quadratic([1.0, 3.0, 2.0]);
    r.reverse();
    assert_eq!((r.weight(0), r.weight(2)), (2.0, 1.0));
    assert_eq!(cubic().reversed_parameter(0.3), 0.7);
    assert_eq!(cubic().reversed_parameter(0.0), 1.0);
}

#[test]
fn resolution_is_positive() {
    assert!(cubic().resolution(1.0) > 0.0);
}

#[test]
fn transform() {
    let c = cubic();
    let moved = c.transformed(&Transform2::translation(Vector2::new(10.0, 20.0)));
    let (a, b) = (c.point(0.5), moved.point(0.5));
    assert!((b.x - a.x - 10.0).abs() <= 1e-10 && (b.y - a.y - 20.0).abs() <= 1e-10);
}

#[test]
fn rational_quarter_circle() {
    let c = quarter_circle();
    assert!(c.is_rational());
    assert!((c.point(0.5).distance(p2(0.0, 0.0)) - 1.0).abs() <= 1e-6);
}

#[test]
fn closed_curve() {
    let c = BezierCurve2::new(vec![p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 0.0), p2(0.0, 0.0)]);
    assert!(c.is_closed());
}

#[test]
fn linear_curve() {
    let c = BezierCurve2::new(vec![p2(0.0, 0.0), p2(3.0, 4.0)]);
    assert_eq!(c.degree(), 1);
    assert_point2_near(c.point(0.5), p2(1.5, 2.0), 1e-10);
}
