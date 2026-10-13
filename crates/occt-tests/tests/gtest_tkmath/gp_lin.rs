//! gp_Lin_Test.cxx (gp_Lin and gp_Lin2d)

use std::f64::consts::FRAC_PI_2;

use occt_tests::assert_near;
use kentron::geometry::curve::Line;
use kentron::geometry::curve2::Line2;
use kentron::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use kentron::tolerance;

#[test]
fn constructor_from_ax1() {
    let line = Line::from(Axis::new(Point::new(1.0, 2.0, 3.0), Direction::X));
    assert!(line.origin.approx_eq(Point::new(1.0, 2.0, 3.0), tolerance::LINEAR));
    assert!(line.direction.approx_eq(Direction::X, tolerance::ANGULAR));
}

#[test]
fn constructor_from_point_and_dir() {
    let line = Line::new(Point::ORIGIN, Direction::Z);
    assert!(line.origin.approx_eq(Point::ORIGIN, tolerance::LINEAR));
    assert!(line.direction.approx_eq(Direction::Z, tolerance::ANGULAR));
}

#[test]
fn distance_to_point() {
    let line = Line::new(Point::ORIGIN, Direction::X);
    // distance from the point to the X axis: sqrt(3^2 + 4^2) = 5
    assert_near!(line.distance_to_point(Point::new(5.0, 3.0, 4.0)), 5.0, tolerance::LINEAR);
}

#[test]
fn distance_to_line_parallel() {
    let l1 = Line::new(Point::ORIGIN, Direction::X);
    let l2 = Line::new(Point::new(0.0, 3.0, 4.0), Direction::X);
    assert_near!(l1.distance_to_line(&l2), 5.0, tolerance::LINEAR);
}

#[test]
fn distance_to_line_intersecting() {
    let l1 = Line::new(Point::ORIGIN, Direction::X);
    let l2 = Line::new(Point::ORIGIN, Direction::Y);
    assert_near!(l1.distance_to_line(&l2), 0.0, tolerance::LINEAR);
}

#[test]
fn contains() {
    let line = Line::new(Point::ORIGIN, Direction::X);
    assert!(line.contains(Point::new(5.0, 0.0, 0.0), tolerance::LINEAR));
    assert!(!line.contains(Point::new(5.0, 1.0, 0.0), tolerance::LINEAR));
}

#[test]
fn angle() {
    let l1 = Line::new(Point::ORIGIN, Direction::X);
    let l2 = Line::new(Point::ORIGIN, Direction::Y);
    assert_near!(l1.angle(&l2), FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn translate() {
    let translated = Line::new(Point::ORIGIN, Direction::X).translated(Vector::new(0.0, 5.0, 0.0));
    assert_near!(translated.origin.y, 5.0, tolerance::LINEAR);
}

#[test]
fn rotate() {
    let rotated = Line::new(Point::new(1.0, 0.0, 0.0), Direction::X).rotated(Axis::Z, FRAC_PI_2);
    assert_near!(rotated.origin.x, 0.0, tolerance::LINEAR);
    assert_near!(rotated.origin.y, 1.0, tolerance::LINEAR);
    assert!(rotated.direction.approx_eq(Direction::Y, tolerance::ANGULAR));
}

#[test]
fn mirror() {
    let mirrored = Line::new(Point::new(1.0, 1.0, 0.0), Direction::X).mirrored_axis(Axis::X);
    assert_near!(mirrored.origin.y, -1.0, tolerance::LINEAR);
}

#[test]
fn transform() {
    let t = Transform::translation(Vector::new(0.0, 0.0, 10.0));
    let transformed = Line::new(Point::ORIGIN, Direction::X).transformed(&t);
    assert_near!(transformed.origin.z, 10.0, tolerance::LINEAR);
}

/// OCC15489: a 2D line from the equation a*x + b*y + c = 0 needs a non-zero (a, b).
#[test]
fn lin2d_construct_from_equation_zero_direction_throws_exception() {
    assert!(Line2::from_equation(0.0, 0.0, 1.0).is_none());
}

#[test]
fn lin2d_construct_from_equation_valid_coefficients_correct_origin() {
    // 1e-20 * x - y + 2 = 0: the origin (closest point to (0, 0)) must be
    // x ~= -1.9999999999999999e-20, y ~= 2
    let line = Line2::from_equation(1e-20, -1.0, 2.0).unwrap();
    assert_near!(line.origin.x, -1.9999999999999999e-20, 1e-25);
    assert_near!(line.origin.y, 2.0, 0.001);
}
