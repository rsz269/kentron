//! Geom_BezierCurve_Test.cxx.
//!
//! Not ported: copy tests and `WeightsArray_*` (see GTEST_PORTING.md).

use kentron::geometry::curve::BezierCurve;
use kentron::geometry::{Continuity, ParametricCurve};
use kentron::math::{Point, Transform, Transformable, Vector};

use crate::{assert_point_near, p};

fn cubic() -> BezierCurve {
    BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0), p(3.0, 0.0, 0.0)])
}

fn quadratic(weights: [f64; 3]) -> BezierCurve {
    BezierCurve::rational(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)], weights.to_vec())
}

fn quarter_circle() -> BezierCurve {
    BezierCurve::rational(vec![p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)], vec![1.0, 0.5_f64.sqrt(), 1.0])
}

#[test]
fn evaluation() {
    let c = cubic();
    assert_point_near(c.point(0.0), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.point(1.0), p(3.0, 0.0, 0.0), 1e-10);
    assert!(c.d1(0.5).1.length() > 0.0);
    assert!(c.d2(0.5).1.length() > 0.0);
    let (a, b) = (c.d3(0.5).3, c.d3(0.25).3);
    assert!((a - b).length() <= 1e-8);
    assert!(c.dn(0.5, 1).length() > 0.0);
    assert!(c.dn(0.5, 4).length() <= 1e-10);
}

#[test]
fn properties() {
    let c = cubic();
    assert_point_near(c.start_point(), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.end_point(), p(3.0, 0.0, 0.0), 1e-10);
    assert_eq!(c.degree(), 3);
    assert!(!c.is_periodic() && !c.is_rational() && !c.is_closed());
    assert!(c.is_cn(100));
    assert_eq!(c.continuity(), Continuity::CN);
    assert_eq!((c.first_parameter(), c.last_parameter()), (0.0, 1.0));
    assert_eq!(c.poles().len(), 4);
    assert_eq!(c.pole(0), p(0.0, 0.0, 0.0));
    assert_eq!(c.weights(), None);
    assert!(BezierCurve::<Point>::MAX_DEGREE >= 25);
}

#[test]
fn set_pole_and_weight() {
    let mut c = cubic();
    c.set_pole(1, p(1.0, 5.0, 0.0));
    assert_eq!(c.pole(1), p(1.0, 5.0, 0.0));
    assert_point_near(c.start_point(), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.end_point(), p(3.0, 0.0, 0.0), 1e-10);
    let mut r = quadratic([1.0, 1.0, 1.0]);
    let before = r.point(0.5);
    r.set_weight(1, 10.0);
    assert_eq!(r.weight(1), 10.0);
    assert!(r.is_rational());
    assert!(r.point(0.5).y > before.y);
    let mut r = quadratic([1.0, 1.0, 1.0]);
    r.set_pole_with_weight(1, p(1.0, 2.0, 0.0), 5.0);
    assert_eq!((r.pole(1), r.weight(1)), (p(1.0, 2.0, 0.0), 5.0));
}

#[test]
fn insert_and_remove_poles() {
    let mut c = cubic();
    c.insert_pole_after(1, p(1.5, 2.0, 0.0));
    assert_eq!(c.poles().len(), 5);
    assert_eq!(c.degree(), 4);
    assert_eq!(c.pole(2), p(1.5, 2.0, 0.0));
    c.remove_pole(2);
    assert_eq!(c.poles().len(), 4);
    assert_eq!(c.degree(), 3);
    let mut c = cubic();
    c.insert_pole_before(1, p(0.5, 2.0, 0.0));
    assert_eq!(c.poles().len(), 5);
    assert_eq!(c.pole(1), p(0.5, 2.0, 0.0));
    let mut r = quadratic([1.0, 2.0, 1.0]);
    r.insert_pole_after_with_weight(1, p(1.5, 0.5, 0.0), 3.0);
    assert_eq!(r.poles().len(), 4);
    assert_eq!(r.weight(2), 3.0);
    assert!(r.is_rational());
}

#[test]
fn increase_degree_keeps_shape() {
    let mut c = cubic();
    let before = c.point(0.5);
    c.increase_degree(5);
    assert_eq!((c.degree(), c.poles().len()), (5, 6));
    assert_point_near(c.point(0.5), before, 1e-10);
    let mut r = quadratic([1.0, 2.0, 1.0]);
    let before = r.point(0.5);
    r.increase_degree(5);
    assert!(r.is_rational());
    assert_point_near(r.point(0.5), before, 1e-10);
}

#[test]
fn segment() {
    let mut c = cubic();
    let (a, b) = (c.point(0.25), c.point(0.75));
    c.segment(0.25, 0.75);
    assert_point_near(c.start_point(), a, 1e-6);
    assert_point_near(c.end_point(), b, 1e-6);
    let mut r = quarter_circle();
    let mid = r.point(0.5);
    r.segment(0.25, 0.75);
    assert_point_near(r.point(0.5), mid, 1e-6);
    assert!(r.is_rational());
}

#[test]
fn reverse() {
    let mut c = cubic();
    let (s, e) = (c.start_point(), c.end_point());
    c.reverse();
    assert_point_near(c.start_point(), e, 1e-10);
    assert_point_near(c.end_point(), s, 1e-10);
    let mut r = BezierCurve::rational(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)], vec![1.0, 3.0, 2.0]);
    r.reverse();
    assert_eq!((r.weight(0), r.weight(2)), (2.0, 1.0));
    assert_eq!(cubic().reversed_parameter(0.3), 0.7);
    assert_eq!(cubic().reversed_parameter(0.0), 1.0);
}

#[test]
fn resolution_and_transform() {
    assert!(cubic().resolution(1.0) > 0.0);
    let c = cubic();
    let moved = c.transformed(&Transform::translation(Vector::new(10.0, 20.0, 30.0)));
    assert!(((moved.point(0.5) - c.point(0.5)) - Vector::new(10.0, 20.0, 30.0)).length() <= 1e-10);
}

#[test]
fn rational_quarter_circle() {
    let c = quarter_circle();
    assert!(c.is_rational());
    assert!((c.point(0.5).distance(p(0.0, 0.0, 0.0)) - 1.0).abs() <= 1e-6);
}

#[test]
fn closed_and_linear_curves() {
    assert!(BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0), p(0.0, 0.0, 0.0)]).is_closed());
    let line = BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(3.0, 4.0, 0.0)]);
    assert_eq!(line.degree(), 1);
    assert_point_near(line.point(0.5), p(1.5, 2.0, 0.0), 1e-10);
    assert!(line.d2(0.5).2.length() <= 1e-10);
}

#[test]
fn degree_limit() {
    let poles = |n: usize| (1..=n).map(|i| p(i as f64 + 10.0, i as f64 * 2.0 + 20.0, i as f64 * 3.0 + 45.0)).collect::<Vec<_>>();
    assert_eq!(BezierCurve::new(poles(26)).degree(), 25);
    assert!(BezierCurve::try_new(poles(29)).is_none());
    occt_tests::assert_panics!(BezierCurve::new(poles(29)));
}
