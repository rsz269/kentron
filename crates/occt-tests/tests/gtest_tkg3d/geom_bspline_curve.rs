//! Geom_BSplineCurve_Test.cxx.
//!
//! Not ported: copy tests and `WeightsArray_*`. `FirstUKnotIndex`/`LastUKnotIndex` are the
//! first and last indices of `knots()`.

use kentron::geometry::curve::BSplineCurve;
use kentron::geometry::{Continuity, KnotDistribution, ParametricCurve};
use kentron::math::{Point, Transform, Transformable, Vector};

use crate::{assert_point_near, p};

fn cubic() -> BSplineCurve {
    BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0), p(3.0, 0.0, 0.0)], vec![0.0, 1.0], vec![4, 4], 3)
}

fn rational_quadratic(poles: [Point; 3], weights: [f64; 3]) -> BSplineCurve {
    BSplineCurve::new(poles.to_vec(), vec![0.0, 1.0], vec![3, 3], 2).with_weights(weights.to_vec())
}

fn pentagon() -> BSplineCurve {
    let poles = vec![p(1.0, 0.0, 0.0), p(0.309, 0.951, 0.0), p(-0.809, 0.588, 0.0), p(-0.809, -0.588, 0.0), p(0.309, -0.951, 0.0)];
    BSplineCurve::new_periodic(poles, (0..6).map(|i| i as f64 * 0.2).collect(), vec![1; 6], 3)
}

#[test]
fn evaluation() {
    let c = cubic();
    assert_point_near(c.point(0.0), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.point(1.0), p(3.0, 0.0, 0.0), 1e-10);
    assert!(c.d1(0.5).1.length() > 0.0);
    assert!(c.d2(0.5).1.length() > 0.0);
    assert!((c.d3(0.5).3 - c.d3(0.25).3).length() <= 1e-8);
    assert!(c.dn(0.5, 1).length() > 0.0);
    assert!(c.dn(0.5, 4).length() <= 1e-10);
    assert_point_near(c.start_point(), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.end_point(), p(3.0, 0.0, 0.0), 1e-10);
}

#[test]
fn properties_and_access() {
    let c = cubic();
    assert_eq!(c.degree(), 3);
    assert!(!c.is_periodic() && !c.is_rational() && !c.is_closed());
    assert!(c.is_cn(3));
    assert_eq!(c.continuity(), Continuity::CN);
    assert_eq!(c.knots(), &[0.0, 1.0]);
    assert_eq!(c.flat_knots().len(), 8);
    assert_eq!(c.multiplicities(), &[4, 4]);
    assert_eq!(c.poles().len(), 4);
    assert_eq!(c.pole(0), p(0.0, 0.0, 0.0));
    assert_eq!(c.weights(), None);
    assert_eq!(c.knot_distribution(), KnotDistribution::PiecewiseBezier);
}

#[test]
fn rational_weights() {
    let r = rational_quadratic([p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)], [1.0, 2.0, 1.0]);
    assert_eq!(r.weights(), Some(&[1.0, 2.0, 1.0][..]));
    let mut r = rational_quadratic([p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)], [1.0, 1.0, 1.0]);
    let before = r.point(0.5);
    r.set_weight(1, 5.0);
    assert_eq!(r.weight(1), 5.0);
    assert!(r.is_rational());
    assert!(r.point(0.5).distance(before) > 1e-6);
}

#[test]
fn set_pole() {
    let mut c = cubic();
    c.set_pole(1, p(1.0, 2.0, 3.0));
    assert_eq!(c.pole(1), p(1.0, 2.0, 3.0));
    assert_point_near(c.start_point(), p(0.0, 0.0, 0.0), 1e-10);
}

#[test]
fn degree_and_knot_insertion_keep_shape() {
    let before = cubic().point(0.5);
    let mut c = cubic();
    c.increase_degree(5);
    assert_eq!(c.degree(), 5);
    assert_point_near(c.point(0.5), before, 1e-10);
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    assert_eq!(c.knots().len(), 3);
    assert_point_near(c.point(0.5), before, 1e-10);
    let mut c = cubic();
    c.insert_knots(&[(0.25, 1), (0.75, 1)]);
    assert_eq!(c.knots().len(), 4);
    assert_point_near(c.point(0.5), before, 1e-10);
    let mut r = rational_quadratic([p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)], [1.0, 2.0, 1.0]);
    let before = r.point(0.5);
    r.increase_degree(4);
    assert_eq!(r.degree(), 4);
    assert!(r.is_rational());
    assert_point_near(r.point(0.5), before, 1e-10);
}

#[test]
fn remove_knot() {
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    let before = c.point(0.25);
    assert!(c.remove_knot(1, 0, 1e-6));
    assert_eq!(c.knots().len(), 2);
    assert_point_near(c.point(0.25), before, 1e-6);
}

#[test]
fn increase_multiplicity() {
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    assert_eq!(c.multiplicity(1), 1);
    let before = c.point(0.5);
    c.increase_multiplicity(1, 2);
    assert_eq!(c.multiplicity(1), 2);
    assert_point_near(c.point(0.5), before, 1e-10);
}

#[test]
fn segment_keeps_parametrization() {
    let mut c = cubic();
    let (a, b) = (c.point(0.25), c.point(0.75));
    c.segment(0.25, 0.75);
    assert!((c.first_parameter() - 0.25).abs() <= 1e-10);
    assert!((c.last_parameter() - 0.75).abs() <= 1e-10);
    assert_point_near(c.start_point(), a, 1e-6);
    assert_point_near(c.end_point(), b, 1e-6);
    let mut r = rational_quadratic([p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)], [1.0, 0.5_f64.sqrt(), 1.0]);
    let mid = r.point(0.5);
    r.segment(0.25, 0.75);
    assert_point_near(r.point(0.5), mid, 1e-6);
    assert!(r.is_rational());
}

#[test]
fn reverse_resolution_transform() {
    let mut c = cubic();
    let (s, e) = (c.start_point(), c.end_point());
    c.reverse();
    assert_point_near(c.start_point(), e, 1e-10);
    assert_point_near(c.end_point(), s, 1e-10);
    assert!(cubic().resolution(1.0) > 0.0);
    let moved = cubic().transformed(&Transform::translation(Vector::new(10.0, 20.0, 30.0)));
    assert!(((moved.point(0.5) - cubic().point(0.5)) - Vector::new(10.0, 20.0, 30.0)).length() <= 1e-10);
}

#[test]
fn periodic_curve() {
    let mut c = pentagon();
    assert!(c.is_periodic());
    let before = c.point(0.5);
    c.set_not_periodic();
    assert!(!c.is_periodic());
    assert_point_near(c.point(0.5), before, 1e-10);
    let mut c = pentagon();
    c.set_origin(2);
    assert!(c.is_periodic());
    let q = c.point(c.first_parameter());
    assert!(q.x * q.x + q.y * q.y > 0.0);
}

#[test]
fn multi_knot_span() {
    let poles = vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 0.0, 0.0), p(3.0, 2.0, 0.0), p(4.0, 0.0, 0.0)];
    let c = BSplineCurve::new(poles, vec![0.0, 1.0, 2.0, 3.0], vec![3, 1, 1, 3], 2);
    assert_eq!((c.knots().len(), c.poles().len()), (4, 5));
    assert_point_near(c.point(0.0), p(0.0, 0.0, 0.0), 1e-10);
    assert_point_near(c.point(3.0), p(4.0, 0.0, 0.0), 1e-10);
}

#[test]
fn local_evaluation_matches_global() {
    let mut c = cubic();
    c.insert_knot(0.5, 1);
    assert_point_near(c.local_point(0.25, 0, 1), c.point(0.25), 1e-10);
    let ((q, d), (ql, dl)) = (c.d1(0.25), c.local_d1(0.25, 0, 1));
    assert_point_near(q, ql, 1e-10);
    assert!((d - dl).length() <= 1e-10);
}

#[test]
fn set_knot() {
    let poles = vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0), p(3.0, 1.0, 0.0)];
    let mut c = BSplineCurve::new(poles, vec![0.0, 0.5, 1.0], vec![3, 1, 3], 2);
    c.set_knot(1, 0.6);
    assert_eq!(c.knot(1), 0.6);
}

#[test]
fn move_point() {
    let mut c = cubic();
    let target = p(1.5, 0.8, 0.0);
    let last = c.poles().len() - 1;
    assert!(c.move_point(0.5, target, 0, last).is_some());
    assert_point_near(c.point(0.5), target, 1e-6);
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
fn is_equal() {
    let mut other = cubic();
    assert!(cubic().is_equal(&other, 1e-10));
    other.set_pole(1, p(10.0, 10.0, 10.0));
    assert!(!cubic().is_equal(&other, 1e-10));
}

#[test]
fn closed_curve() {
    let c = BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0), p(0.0, 0.0, 0.0)], vec![0.0, 1.0], vec![4, 4], 3);
    assert!(c.is_closed());
}
