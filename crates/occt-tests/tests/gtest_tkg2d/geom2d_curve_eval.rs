//! Geom2d_CurveEval_Test.cxx: evaluation through the curve trait, including degenerate
//! offset points.

use std::f64::consts::FRAC_PI_4;

use kentron::geometry::curve2::{Circle2, Curve2, Ellipse2, Line2, OffsetCurve2};
use kentron::geometry::{fit, ParametricCurve};
use kentron::math::{Direction2, Frame2};
use kentron::tolerance;

use crate::{assert_point2_near, p2};

const TOL: f64 = tolerance::LINEAR;

#[test]
fn circle_d0_d1() {
    let c = Circle2::new(Frame2::WORLD, 4.0);
    let p = c.point(FRAC_PI_4);
    assert!((p.distance(p2(0.0, 0.0)) - 4.0).abs() <= TOL);
    let (q, d1) = c.d1(FRAC_PI_4);
    assert_point2_near(q, p, TOL);
    assert!((q - p2(0.0, 0.0)).dot(d1).abs() <= 1e-10);
}

#[test]
fn offset_to_the_center_evaluates_to_the_center() {
    // Offsetting a radius-3 circle by 3 inwards (to the left of a counter-clockwise
    // circle, so negative, as in OCCT's test) collapses it to its center.
    let collapsed = OffsetCurve2::new(Curve2::Circle(Circle2::new(Frame2::WORLD, 3.0)), -3.0);
    assert_point2_near(collapsed.point(0.0), p2(0.0, 0.0), TOL);
    assert_point2_near(collapsed.d1(0.0).0, p2(0.0, 0.0), TOL);
}

#[test]
fn bspline_d1_consistent_with_derivative_order_one() {
    let points = [p2(0.0, 0.0), p2(1.0, 2.0), p2(3.0, 1.0), p2(5.0, 3.0), p2(7.0, 0.0)];
    let c = fit::approximate(&points).unwrap();
    let mid = 0.5 * (c.first_parameter() + c.last_parameter());
    let (p, d1) = c.d1(mid);
    assert_point2_near(p, c.point(mid), TOL);
    assert!((d1 - c.dn(mid, 1)).length() <= TOL);
}

#[test]
fn ellipse_higher_derivatives() {
    let e = Ellipse2::new(Frame2::WORLD, 5.0, 3.0);
    assert!(e.d2(FRAC_PI_4).2.length() > 0.0);
    assert!(e.d3(FRAC_PI_4).3.length() > 0.0);
}

#[test]
fn line_dn() {
    let l = Line2::new(p2(1.0, 2.0), Direction2::X);
    let d1 = l.dn(5.0, 1);
    assert!((d1.x - 1.0).abs() <= TOL && d1.y.abs() <= TOL);
    assert!(l.dn(5.0, 2).length() <= TOL);
}

#[test]
fn derivative_of_order_zero_panics() {
    let c = fit::approximate(&[p2(0.0, 0.0), p2(1.0, 1.0), p2(2.0, 0.0)]).unwrap();
    let mid = 0.5 * (c.first_parameter() + c.last_parameter());
    occt_tests::assert_panics!(c.dn(mid, 0));
}
