//! gp_Pln_Test.cxx

use std::f64::consts::FRAC_PI_2;

use occt_tests::assert_near;
use kentron::geometry::surface::Plane;
use kentron::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};
use kentron::tolerance;

#[test]
fn constructor_from_ax3() {
    let plane = Plane::from_frame(Frame::from_z(Point::ORIGIN, Direction::Z));
    assert!(plane.origin().approx_eq(Point::ORIGIN, tolerance::LINEAR));
    assert!(plane.normal().approx_eq(Direction::Z, tolerance::ANGULAR));
}

#[test]
fn constructor_from_point_and_dir() {
    let origin = Point::new(1.0, 2.0, 3.0);
    let plane = Plane::new(origin, Direction::Z);
    assert!(plane.origin().approx_eq(origin, tolerance::LINEAR));
    assert!(plane.normal().approx_eq(Direction::Z, tolerance::ANGULAR));
}

#[test]
fn constructor_from_coefficients() {
    // plane z = 5: 0*x + 0*y + 1*z - 5 = 0
    let (a, b, c, d) = Plane::from_equation(0.0, 0.0, 1.0, -5.0).unwrap().equation();
    assert_near!(a, 0.0, tolerance::LINEAR);
    assert_near!(b, 0.0, tolerance::LINEAR);
    assert_near!(c, 1.0, tolerance::LINEAR);
    assert_near!(d, -5.0, tolerance::LINEAR);
}

#[test]
fn distance_to_point() {
    let plane = Plane::new(Point::ORIGIN, Direction::Z);
    assert_near!(plane.distance(Point::new(3.0, 4.0, 5.0)), 5.0, tolerance::LINEAR);
}

#[test]
fn contains() {
    let plane = Plane::new(Point::ORIGIN, Direction::Z);
    assert!(plane.contains(Point::new(3.0, 4.0, 0.0), tolerance::LINEAR));
    assert!(!plane.contains(Point::new(3.0, 4.0, 1.0), tolerance::LINEAR));
}

#[test]
fn coefficients() {
    let (a, b, c, d) = Plane::new(Point::ORIGIN, Direction::Z).equation();
    assert_near!(a, 0.0, tolerance::LINEAR);
    assert_near!(b, 0.0, tolerance::LINEAR);
    assert_near!(c, 1.0, tolerance::LINEAR);
    assert_near!(d, 0.0, tolerance::LINEAR);
}

#[test]
fn mirror() {
    let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Direction::Z);
    let mirrored = plane.mirrored_plane(Point::ORIGIN, Direction::Z);
    assert_near!(mirrored.origin().z, -5.0, tolerance::LINEAR);
}

#[test]
fn rotate() {
    let rotated = Plane::new(Point::new(1.0, 0.0, 0.0), Direction::Z).rotated(Axis::Z, FRAC_PI_2);
    assert_near!(rotated.origin().x, 0.0, tolerance::LINEAR);
    assert_near!(rotated.origin().y, 1.0, tolerance::LINEAR);
}

#[test]
fn scale() {
    let scaled = Plane::new(Point::new(1.0, 0.0, 0.0), Direction::Z).scaled(Point::ORIGIN, 3.0);
    assert_near!(scaled.origin().x, 3.0, tolerance::LINEAR);
}

#[test]
fn translate() {
    let translated = Plane::new(Point::ORIGIN, Direction::Z).translated(Vector::new(1.0, 2.0, 3.0));
    assert_near!(translated.origin().x, 1.0, tolerance::LINEAR);
    assert_near!(translated.origin().y, 2.0, tolerance::LINEAR);
    assert_near!(translated.origin().z, 3.0, tolerance::LINEAR);
}

#[test]
fn transform() {
    let t = Transform::translation(Vector::new(0.0, 0.0, 10.0));
    let transformed = Plane::new(Point::ORIGIN, Direction::Z).transformed(&t);
    assert_near!(transformed.origin().z, 10.0, tolerance::LINEAR);
}
