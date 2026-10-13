//! gp_Vec2d_Test.cxx

use std::f64::consts::FRAC_PI_2;

use occt_tests::{assert_double_eq, assert_near};
use kentron::math::{Direction2, Point2, Vector2};
use kentron::tolerance;

#[test]
fn default_constructor() {
    let v = Vector2::default();
    assert_double_eq!(v.x, 0.0);
    assert_double_eq!(v.y, 0.0);
}

#[test]
fn coordinate_constructor() {
    let v = Vector2::new(1.0, 2.0);
    assert_double_eq!(v.x, 1.0);
    assert_double_eq!(v.y, 2.0);
}

#[test]
fn two_point_constructor() {
    let v = Point2::new(4.0, 6.0) - Point2::new(1.0, 2.0);
    assert_double_eq!(v.x, 3.0);
    assert_double_eq!(v.y, 4.0);
}

#[test]
fn dir2d_constructor() {
    let v = Vector2::from(Direction2::Y);
    assert_near!(v.length(), 1.0, tolerance::LINEAR);
    assert_near!(v.y, 1.0, tolerance::LINEAR);
}

#[test]
fn add_subtract() {
    let v1 = Vector2::new(1.0, 2.0);
    let v2 = Vector2::new(4.0, 5.0);

    let sum = v1 + v2;
    assert_double_eq!(sum.x, 5.0);
    assert_double_eq!(sum.y, 7.0);

    let diff = v2 - v1;
    assert_double_eq!(diff.x, 3.0);
    assert_double_eq!(diff.y, 3.0);
}

#[test]
fn multiply_divide() {
    let v = Vector2::new(2.0, 4.0);

    let scaled = v * 3.0;
    assert_double_eq!(scaled.x, 6.0);
    assert_double_eq!(scaled.y, 12.0);

    let divided = v / 2.0;
    assert_double_eq!(divided.x, 1.0);
    assert_double_eq!(divided.y, 2.0);
}

#[test]
fn dot() {
    assert_double_eq!(Vector2::new(1.0, 0.0).dot(Vector2::new(0.0, 1.0)), 0.0);
    // 1*4 + 2*5 = 14
    assert_double_eq!(Vector2::new(1.0, 2.0).dot(Vector2::new(4.0, 5.0)), 14.0);
}

#[test]
fn crossed() {
    // 1*1 - 0*0 = 1
    assert_near!(Vector2::new(1.0, 0.0).cross(Vector2::new(0.0, 1.0)), 1.0, tolerance::LINEAR);
    // 3*5 - 4*2 = 7
    assert_near!(Vector2::new(3.0, 4.0).cross(Vector2::new(2.0, 5.0)), 7.0, tolerance::LINEAR);
}

#[test]
fn magnitude() {
    let v = Vector2::new(3.0, 4.0);
    assert_near!(v.length(), 5.0, tolerance::LINEAR);
    assert_near!(v.length_squared(), 25.0, tolerance::LINEAR);
}

#[test]
fn normalize() {
    let n = Vector2::new(3.0, 4.0).normalized();
    assert_near!(n.length(), 1.0, tolerance::LINEAR);
    assert_near!(n.x, 0.6, tolerance::LINEAR);
    assert_near!(n.y, 0.8, tolerance::LINEAR);
}

#[test]
fn angle() {
    let a = Vector2::new(1.0, 0.0).angle(Vector2::new(0.0, 1.0));
    assert_near!(a, FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn is_parallel() {
    assert!(Vector2::new(1.0, 0.0).is_parallel(Vector2::new(5.0, 0.0), tolerance::ANGULAR));
}

#[test]
fn is_normal() {
    assert!(Vector2::new(1.0, 0.0).is_perpendicular(Vector2::new(0.0, 1.0), tolerance::ANGULAR));
}

#[test]
fn is_opposite() {
    assert!(Vector2::new(1.0, 0.0).is_opposite(Vector2::new(-3.0, 0.0), tolerance::ANGULAR));
}

#[test]
fn reverse() {
    let r = -Vector2::new(1.0, 2.0);
    assert_double_eq!(r.x, -1.0);
    assert_double_eq!(r.y, -2.0);
}

#[test]
fn set_linear_form() {
    let result = 2.0 * Vector2::new(1.0, 0.0) + 3.0 * Vector2::new(0.0, 1.0);
    assert_near!(result.x, 2.0, tolerance::LINEAR);
    assert_near!(result.y, 3.0, tolerance::LINEAR);
}

/// OCC26750: perpendicularity must also hold for an angle of -PI/2.
#[test]
fn occ26750_is_normal_negative_angle() {
    let v1 = Vector2::new(1.0, 0.0);
    let v2 = Vector2::new(0.0, -1.0);
    assert!(v1.is_perpendicular(v2, tolerance::ANGULAR), "vectors at -90 degrees are perpendicular");
    assert!(
        Direction2::X.is_perpendicular(Direction2::NEG_Y, tolerance::ANGULAR),
        "X and -Y directions are perpendicular"
    );
}
