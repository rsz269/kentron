//! Geom2d_Line_Test.cxx, plus the line cases of Adaptor2d_Line2d_Test.cxx.

use kentron::geometry::curve2::Line2;
use kentron::geometry::{Continuity, ParametricCurve};
use kentron::math::{Axis2, Direction2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

use crate::{assert_point2_near, assert_vector2_near, p2};

const TOL: f64 = tolerance::LINEAR;

fn x_line() -> Line2 {
    Line2::new(p2(0.0, 0.0), Direction2::X)
}

#[test]
fn construct_from_point_and_direction() {
    let line = x_line();
    assert_eq!(line.direction, Direction2::X);
    assert_eq!(line.origin, p2(0.0, 0.0));
}

#[test]
fn construct_from_axis() {
    let line = Line2::from(Axis2::new(p2(3.0, 4.0), Direction2::X));
    assert_eq!(line.origin, p2(3.0, 4.0));
}

#[test]
fn parameter_bounds_are_infinite() {
    assert!(x_line().first_parameter() < -1.0e90);
    assert!(x_line().last_parameter() > 1.0e90);
}

#[test]
fn not_closed_not_periodic() {
    assert!(!x_line().is_closed());
    assert!(!x_line().is_periodic());
    assert_eq!(x_line().period(), None);
}

#[test]
fn continuity() {
    assert_eq!(x_line().continuity(), Continuity::CN);
    assert!(x_line().is_cn(100));
}

#[test]
fn evaluation() {
    let line = x_line();
    assert_point2_near(line.point(5.0), p2(5.0, 0.0), TOL);
    let (_, d1) = line.d1(5.0);
    assert_vector2_near(d1, Vector2::new(1.0, 0.0), TOL);
    let (_, _, d2) = line.d2(5.0);
    assert!(d2.length() <= TOL);
    let (_, _, _, d3) = line.d3(5.0);
    assert!(d3.length() <= TOL);
    assert!(line.dn(0.0, 4).length() <= TOL);
    assert_vector2_near(line.dn(0.0, 1), Vector2::new(1.0, 0.0), TOL);
}

#[test]
fn distance_to_point() {
    assert!((x_line().distance_to_point(p2(5.0, 3.0)) - 3.0).abs() <= TOL);
    assert!(x_line().distance_to_point(p2(0.0, 0.0)) <= TOL);
}

#[test]
fn transform_translation_and_rotation() {
    let moved = x_line().transformed(&Transform2::translation(Vector2::new(5.0, 10.0)));
    assert_point2_near(moved.origin, p2(5.0, 10.0), TOL);
    let turned = x_line().transformed(&Transform2::rotation(p2(0.0, 0.0), std::f64::consts::FRAC_PI_2));
    assert!(turned.direction.x().abs() <= TOL);
    assert!((turned.direction.y().abs() - 1.0).abs() <= TOL);
}

#[test]
fn reversed_parameter() {
    assert_eq!(x_line().reversed_parameter(5.0), -5.0);
}

#[test]
fn diagonal_line_value() {
    let diagonal = Line2::new(p2(0.0, 0.0), Direction2::new(1.0, 1.0));
    assert_point2_near(diagonal.point(2.0_f64.sqrt()), p2(1.0, 1.0), TOL);
}

#[test]
fn resolution_of_a_unit_speed_line_is_the_tolerance() {
    assert!((x_line().resolution(0.01) - 0.01).abs() <= 1.0e-12);
}
