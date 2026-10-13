//! gp_Circ_Test.cxx

use std::f64::consts::{FRAC_PI_2, PI};

use occt_tests::assert_near;
use kentron::geometry::curve::Circle;
use kentron::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};
use kentron::tolerance;

fn circle_at(center: Point, radius: f64) -> Circle {
    Circle::new(Frame::from_z(center, Direction::Z), radius)
}

#[test]
fn constructor() {
    assert_near!(circle_at(Point::ORIGIN, 5.0).radius(), 5.0, tolerance::LINEAR);
}

#[test]
fn area() {
    assert_near!(circle_at(Point::ORIGIN, 1.0).area(), PI, tolerance::LINEAR);
}

#[test]
fn length() {
    assert_near!(circle_at(Point::ORIGIN, 1.0).circumference(), 2.0 * PI, tolerance::LINEAR);
}

#[test]
fn center() {
    let center = Point::new(1.0, 2.0, 3.0);
    assert!(circle_at(center, 5.0).center().approx_eq(center, tolerance::LINEAR));
}

#[test]
fn translate() {
    let translated = circle_at(Point::ORIGIN, 5.0).translated(Vector::new(10.0, 20.0, 30.0));
    assert_near!(translated.center().x, 10.0, tolerance::LINEAR);
    assert_near!(translated.center().y, 20.0, tolerance::LINEAR);
    assert_near!(translated.center().z, 30.0, tolerance::LINEAR);
    assert_near!(translated.radius(), 5.0, tolerance::LINEAR);
}

#[test]
fn rotate() {
    let rotated = circle_at(Point::new(1.0, 0.0, 0.0), 2.0).rotated(Axis::Z, FRAC_PI_2);
    assert_near!(rotated.center().x, 0.0, tolerance::LINEAR);
    assert_near!(rotated.center().y, 1.0, tolerance::LINEAR);
    assert_near!(rotated.radius(), 2.0, tolerance::LINEAR);
}

#[test]
fn mirror() {
    let mirrored = circle_at(Point::new(1.0, 2.0, 3.0), 5.0).mirrored_point(Point::ORIGIN);
    assert_near!(mirrored.center().x, -1.0, tolerance::LINEAR);
    assert_near!(mirrored.center().y, -2.0, tolerance::LINEAR);
    assert_near!(mirrored.center().z, -3.0, tolerance::LINEAR);
}

#[test]
fn transform() {
    let t = Transform::scaling(Point::ORIGIN, 2.0);
    let transformed = circle_at(Point::ORIGIN, 3.0).transformed(&t);
    assert_near!(transformed.radius(), 6.0, tolerance::LINEAR);
}
