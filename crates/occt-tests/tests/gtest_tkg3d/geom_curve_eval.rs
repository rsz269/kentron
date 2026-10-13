//! Geom_CurveEval_Test.cxx: curve evaluation, including degenerate offset curves.

use std::f64::consts::{FRAC_PI_3, FRAC_PI_4, FRAC_PI_6, PI};

use kentron::geometry::curve::{Circle, Curve, Line, OffsetCurve, Parabola, TrimmedCurve};
use kentron::geometry::{fit, ParametricCurve};
use kentron::math::{Direction, Frame};
use kentron::tolerance;

use crate::{assert_point_near, p};

const TOL: f64 = tolerance::LINEAR;

fn circle(radius: f64) -> Circle {
    Circle::new(Frame::WORLD, radius)
}

#[test]
fn circle_point_and_derivatives() {
    let c = circle(5.0);
    assert!((c.point(FRAC_PI_4).distance(p(0.0, 0.0, 0.0)) - 5.0).abs() <= TOL);
    let (q, d1) = c.d1(FRAC_PI_3);
    assert!((q.distance(p(0.0, 0.0, 0.0)) - 5.0).abs() <= TOL);
    assert!((q - p(0.0, 0.0, 0.0)).dot(d1).abs() <= 1.0e-10);
    let (q2, _, d2) = c.d2(FRAC_PI_3);
    assert_point_near(q2, q, TOL);
    assert!((d2.length() - 5.0).abs() <= 1.0e-10);
    let (q3, _, _, d3) = c.d3(FRAC_PI_3);
    assert_point_near(q3, q, TOL);
    assert!(d3.length() > 0.0);
}

#[test]
fn line_second_derivative_is_zero() {
    let (_, d1, d2) = Line::new(p(1.0, 2.0, 3.0), Direction::X).d2(5.0);
    assert!(d2.length() <= TOL);
    assert!((d1.x - 1.0).abs() <= TOL && d1.y.abs() <= TOL && d1.z.abs() <= TOL);
}

#[test]
fn bspline_d1_consistent_with_dn() {
    let points = [p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(3.0, 1.0, 0.0), p(5.0, 3.0, 0.0), p(7.0, 0.0, 0.0)];
    let c = fit::approximate(&points).unwrap();
    let mid = 0.5 * (c.first_parameter() + c.last_parameter());
    let (q, d1) = c.d1(mid);
    assert_point_near(q, c.point(mid), TOL);
    assert!((d1 - c.dn(mid, 1)).length() <= TOL);
}

#[test]
fn dn_of_order_one_matches_d1() {
    let c = circle(3.0);
    assert!((c.d1(FRAC_PI_6).1 - c.dn(FRAC_PI_6, 1)).length() <= TOL);
}

#[test]
fn offset_collapsed_to_the_center() {
    // Offsetting a radius-5 circle about +Z by -5 lands on its center.
    let collapsed = OffsetCurve::new(Curve::Circle(circle(5.0)), -5.0, Direction::Z);
    assert_point_near(collapsed.point(0.0), p(0.0, 0.0, 0.0), TOL);
    assert_point_near(collapsed.d1(0.0).0, p(0.0, 0.0, 0.0), TOL);
    let regular = OffsetCurve::new(Curve::Circle(circle(5.0)), 2.0, Direction::Z);
    let _ = regular.point(0.0);
}

#[test]
fn parabola_third_derivative_is_zero() {
    let parabola = Parabola::new(Frame::WORLD, 2.0);
    assert!(parabola.d3(1.0).3.length() <= TOL);
}

#[test]
fn derivative_of_order_zero_panics() {
    let c = fit::approximate(&[p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 0.0)]).unwrap();
    let mid = 0.5 * (c.first_parameter() + c.last_parameter());
    occt_tests::assert_panics!(c.dn(mid, 0));
}

#[test]
fn trimmed_curve_delegates_to_basis() {
    let basis = Curve::Circle(circle(5.0));
    let trimmed = TrimmedCurve::new(basis.clone(), 0.0, PI);
    let ((bp, bd), (tp, td)) = (basis.d1(FRAC_PI_4), trimmed.d1(FRAC_PI_4));
    assert_point_near(bp, tp, TOL);
    assert!((bd - td).length() <= TOL);
}
