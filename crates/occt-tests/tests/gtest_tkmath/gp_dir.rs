//! gp_Dir_Test.cxx

use std::f64::consts::{FRAC_PI_2, PI};

use occt_tests::assert_near;
use kentron::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use kentron::tolerance;

#[test]
fn coordinate_constructor() {
    let d = Direction::new(0.0, 0.0, 1.0);
    assert_near!(d.x(), 0.0, tolerance::LINEAR);
    assert_near!(d.y(), 0.0, tolerance::LINEAR);
    assert_near!(d.z(), 1.0, tolerance::LINEAR);
}

#[test]
fn xyz_constructor() {
    let d = Direction::new(0.0, 3.0, 4.0);
    assert_near!(d.x(), 0.0, tolerance::LINEAR);
    assert_near!(d.y(), 0.6, tolerance::LINEAR);
    assert_near!(d.z(), 0.8, tolerance::LINEAR);
}

#[test]
fn vec_constructor() {
    let d = Direction::from_vector(Vector::new(3.0, 0.0, 0.0));
    assert_near!(d.x(), 1.0, tolerance::LINEAR);
    assert_near!(d.y(), 0.0, tolerance::LINEAR);
    assert_near!(d.z(), 0.0, tolerance::LINEAR);
}

#[test]
fn predefined_directions() {
    assert!(Direction::new(1.0, 0.0, 0.0).approx_eq(Direction::X, tolerance::ANGULAR));
    assert!(Direction::new(0.0, 1.0, 0.0).approx_eq(Direction::Y, tolerance::ANGULAR));
    assert!(Direction::new(0.0, 0.0, 1.0).approx_eq(Direction::Z, tolerance::ANGULAR));
}

#[test]
fn angle() {
    assert_near!(Direction::X.angle(Direction::Y), FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn angle_with_ref() {
    let angle = Direction::X.signed_angle(Direction::Y, Direction::Z);
    assert_near!(angle, FRAC_PI_2, tolerance::ANGULAR);

    // reversing the reference direction negates the angle
    let reversed = Direction::X.signed_angle(Direction::Y, Direction::NEG_Z);
    assert_near!(reversed, -FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn is_equal() {
    let d1 = Direction::new(1.0, 0.0, 0.0);
    let d2 = Direction::new(1.0, 0.0, 0.0);
    assert!(d1.approx_eq(d2, tolerance::ANGULAR));
}

#[test]
fn is_normal() {
    assert!(Direction::X.is_perpendicular(Direction::Y, tolerance::ANGULAR));
}

#[test]
fn is_opposite() {
    assert!(Direction::X.is_opposite(Direction::NEG_X, tolerance::ANGULAR));
}

#[test]
fn is_parallel() {
    assert!(Direction::X.is_parallel(Direction::NEG_X, tolerance::ANGULAR));
}

#[test]
fn cross() {
    let cross = Direction::from_vector(Direction::X.cross(Direction::Y));
    assert_near!(cross.x(), 0.0, tolerance::LINEAR);
    assert_near!(cross.y(), 0.0, tolerance::LINEAR);
    assert_near!(cross.z(), 1.0, tolerance::LINEAR);
}

#[test]
fn cross_cross() {
    let a = Vector::from(Direction::new(1.0, 0.0, 0.0));
    let b = Vector::from(Direction::new(1.0, 1.0, 0.0));
    let c = Vector::from(Direction::new(0.0, 0.0, 1.0));
    let cc = Direction::from_vector(a.cross(b.cross(c)));
    // the result must be a valid normalized direction
    let mag = (cc.x() * cc.x() + cc.y() * cc.y() + cc.z() * cc.z()).sqrt();
    assert_near!(mag, 1.0, tolerance::LINEAR);
}

#[test]
fn dot() {
    assert_near!(Direction::X.dot(Direction::Y), 0.0, tolerance::LINEAR);
    assert_near!(Direction::X.dot(Direction::new(1.0, 0.0, 0.0)), 1.0, tolerance::LINEAR);
}

#[test]
fn mirror_point() {
    // mirror about the Y axis: 2*(d.a)*a - d = 2*0*(0,1,0) - (1,0,0) = (-1,0,0)
    let mirrored = Direction::X.mirrored_axis(Axis::Y);
    assert_near!(mirrored.x(), -1.0, tolerance::LINEAR);
    assert_near!(mirrored.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn rotate() {
    let rotated = Direction::X.rotated(Axis::Z, FRAC_PI_2);
    assert_near!(rotated.x(), 0.0, tolerance::LINEAR);
    assert_near!(rotated.y(), 1.0, tolerance::LINEAR);
    assert_near!(rotated.z(), 0.0, tolerance::LINEAR);
}

#[test]
fn transform() {
    let t = Transform::rotation(Axis::Z, PI);
    let transformed = Direction::X.transformed(&t);
    assert_near!(transformed.x(), -1.0, tolerance::LINEAR);
    assert_near!(transformed.y(), 0.0, tolerance::LINEAR);
}

#[test]
fn zero_magnitude_input_throws_exception() {
    assert!(Direction::try_new(0.0, 0.0, 0.0).is_none());
}

/// OCC22558: mirroring (1,1,1) in the plane through the origin with normal (1,0,0) gives
/// (-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)).
#[test]
fn occ22558_mirror_through_ax2_correct_result() {
    let mirrored = Direction::new(1.0, 1.0, 1.0).mirrored_plane(Point::ORIGIN, Direction::X);

    let expected = 1.0 / 3.0_f64.sqrt();
    assert_near!(mirrored.x(), -expected, tolerance::LINEAR);
    assert_near!(mirrored.y(), expected, tolerance::LINEAR);
    assert_near!(mirrored.z(), expected, tolerance::LINEAR);
}
