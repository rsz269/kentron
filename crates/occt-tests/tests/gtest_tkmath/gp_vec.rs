//! gp_Vec_Test.cxx

use std::f64::consts::FRAC_PI_2;

use occt_tests::{assert_double_eq, assert_near};
use kentron::math::{Direction, Point, Vector};
use kentron::tolerance;

#[test]
fn default_constructor() {
    let v = Vector::default();
    assert_double_eq!(v.x, 0.0);
    assert_double_eq!(v.y, 0.0);
    assert_double_eq!(v.z, 0.0);
}

#[test]
fn coordinate_constructor() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_double_eq!(v.x, 1.0);
    assert_double_eq!(v.y, 2.0);
    assert_double_eq!(v.z, 3.0);
}

#[test]
fn two_point_constructor() {
    let v = Point::new(4.0, 6.0, 8.0) - Point::new(1.0, 2.0, 3.0);
    assert_double_eq!(v.x, 3.0);
    assert_double_eq!(v.y, 4.0);
    assert_double_eq!(v.z, 5.0);
}

#[test]
fn dir_constructor() {
    let v = Vector::from(Direction::Z);
    assert_near!(v.length(), 1.0, tolerance::LINEAR);
    assert_near!(v.z, 1.0, tolerance::LINEAR);
}

#[test]
fn add_subtract() {
    let v1 = Vector::new(1.0, 2.0, 3.0);
    let v2 = Vector::new(4.0, 5.0, 6.0);

    let sum = v1 + v2;
    assert_double_eq!(sum.x, 5.0);
    assert_double_eq!(sum.y, 7.0);
    assert_double_eq!(sum.z, 9.0);

    let diff = v2 - v1;
    assert_double_eq!(diff.x, 3.0);
    assert_double_eq!(diff.y, 3.0);
    assert_double_eq!(diff.z, 3.0);
}

#[test]
fn multiply_divide() {
    let v = Vector::new(2.0, 4.0, 6.0);

    let scaled = v * 3.0;
    assert_double_eq!(scaled.x, 6.0);
    assert_double_eq!(scaled.y, 12.0);
    assert_double_eq!(scaled.z, 18.0);

    let divided = v / 2.0;
    assert_double_eq!(divided.x, 1.0);
    assert_double_eq!(divided.y, 2.0);
    assert_double_eq!(divided.z, 3.0);
}

#[test]
fn dot() {
    assert_double_eq!(Vector::new(1.0, 0.0, 0.0).dot(Vector::new(0.0, 1.0, 0.0)), 0.0);
    // 1*4 + 2*5 + 3*6 = 32
    assert_double_eq!(Vector::new(1.0, 2.0, 3.0).dot(Vector::new(4.0, 5.0, 6.0)), 32.0);
}

#[test]
fn cross() {
    let cross = Vector::new(1.0, 0.0, 0.0).cross(Vector::new(0.0, 1.0, 0.0));
    assert_near!(cross.x, 0.0, tolerance::LINEAR);
    assert_near!(cross.y, 0.0, tolerance::LINEAR);
    assert_near!(cross.z, 1.0, tolerance::LINEAR);
}

#[test]
fn cross_magnitude() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(0.0, 1.0, 0.0);
    assert_near!(v1.cross(v2).length(), 1.0, tolerance::LINEAR);
}

#[test]
fn dot_cross_triple_scalar_product() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(0.0, 1.0, 0.0);
    let v3 = Vector::new(0.0, 0.0, 1.0);
    // triple scalar product of orthogonal unit vectors = 1
    assert_near!(v1.dot(v2.cross(v3)), 1.0, tolerance::LINEAR);
}

#[test]
fn magnitude() {
    let v = Vector::new(3.0, 4.0, 0.0);
    assert_near!(v.length(), 5.0, tolerance::LINEAR);
    assert_near!(v.length_squared(), 25.0, tolerance::LINEAR);
}

#[test]
fn normalize() {
    let n = Vector::new(3.0, 4.0, 0.0).normalized();
    assert_near!(n.length(), 1.0, tolerance::LINEAR);
    assert_near!(n.x, 0.6, tolerance::LINEAR);
    assert_near!(n.y, 0.8, tolerance::LINEAR);
}

#[test]
fn angle() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(0.0, 1.0, 0.0);
    assert_near!(v1.angle(v2), FRAC_PI_2, tolerance::ANGULAR);
}

#[test]
fn is_parallel() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(5.0, 0.0, 0.0);
    assert!(v1.is_parallel(v2, tolerance::ANGULAR));
}

#[test]
fn is_normal() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(0.0, 1.0, 0.0);
    assert!(v1.is_perpendicular(v2, tolerance::ANGULAR));
}

#[test]
fn is_opposite() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(-3.0, 0.0, 0.0);
    assert!(v1.is_opposite(v2, tolerance::ANGULAR));
}

#[test]
fn reverse() {
    let r = -Vector::new(1.0, 2.0, 3.0);
    assert_double_eq!(r.x, -1.0);
    assert_double_eq!(r.y, -2.0);
    assert_double_eq!(r.z, -3.0);
}

#[test]
fn set_linear_form() {
    let v1 = Vector::new(1.0, 0.0, 0.0);
    let v2 = Vector::new(0.0, 1.0, 0.0);
    let result = 2.0 * v1 + 3.0 * v2;
    assert_near!(result.x, 2.0, tolerance::LINEAR);
    assert_near!(result.y, 3.0, tolerance::LINEAR);
    assert_near!(result.z, 0.0, tolerance::LINEAR);
}
