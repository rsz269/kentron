//! gp_Dir2d_Test.cxx

use std::f64::consts::{FRAC_PI_2, PI};

use occt_tests::assert_near;
use kentron::math::{Axis2, Direction2, Point2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

#[test]
fn coordinate_constructor() {
    let d = Direction2::new(1.0, 0.0);
    assert_near!(d.x(), 1.0, tolerance::LINEAR);
    assert_near!(d.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn xy_constructor() {
    let d = Direction2::new(0.0, 3.0);
    assert_near!(d.x(), 0.0, tolerance::LINEAR);
    assert_near!(d.y(), 1.0, tolerance::LINEAR);
}

#[test]
fn vec2d_constructor() {
    let d = Direction2::from_vector(Vector2::new(3.0, 0.0));
    assert_near!(d.x(), 1.0, tolerance::LINEAR);
    assert_near!(d.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn angle() {
    assert_near!(Direction2::X.angle(Direction2::Y), FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn is_equal() {
    assert!(Direction2::new(1.0, 0.0).approx_eq(Direction2::X, tolerance::ANGULAR));
}

#[test]
fn is_normal() {
    assert!(Direction2::X.is_perpendicular(Direction2::Y, tolerance::ANGULAR));
}

#[test]
fn is_opposite() {
    assert!(Direction2::X.is_opposite(Direction2::NEG_X, tolerance::ANGULAR));
}

#[test]
fn is_parallel() {
    assert!(Direction2::X.is_parallel(Direction2::NEG_X, tolerance::ANGULAR));
}

#[test]
fn dot() {
    assert_near!(Direction2::X.dot(Direction2::Y), 0.0, tolerance::LINEAR);
    assert_near!(Direction2::X.dot(Direction2::new(1.0, 0.0)), 1.0, tolerance::LINEAR);
}

#[test]
fn crossed() {
    // 1*1 - 0*0 = 1
    assert_near!(Direction2::X.cross(Direction2::Y), 1.0, tolerance::LINEAR);
}

#[test]
fn rotate() {
    let rotated = Direction2::X.rotated(Point2::ORIGIN, FRAC_PI_2);
    assert_near!(rotated.x(), 0.0, tolerance::LINEAR);
    assert_near!(rotated.y(), 1.0, tolerance::LINEAR);
}

#[test]
fn mirror() {
    // mirror about the Y axis: 2*(d.a)*a - d = 2*0*(0,1) - (1,0) = (-1,0)
    let mirrored = Direction2::X.mirrored_axis(Axis2::Y);
    assert_near!(mirrored.x(), -1.0, tolerance::LINEAR);
    assert_near!(mirrored.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn transform() {
    let t = Transform2::rotation(Point2::ORIGIN, PI);
    let transformed = Direction2::X.transformed(&t);
    assert_near!(transformed.x(), -1.0, tolerance::LINEAR);
    assert_near!(transformed.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn zero_magnitude_throws() {
    assert!(Direction2::try_new(0.0, 0.0).is_none());
}
