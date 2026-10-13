//! gp_XYZ_Test.cxx (gp_XYZ is folded into `Vector`)

use occt_tests::{assert_double_eq, assert_near};
use kentron::math::Vector;
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
fn set_coord() {
    let mut v = Vector::default();
    (v.x, v.y, v.z) = (7.0, 8.0, 9.0);
    assert_double_eq!(v.x, 7.0);
    assert_double_eq!(v.y, 8.0);
    assert_double_eq!(v.z, 9.0);
}

#[test]
fn accessors() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_double_eq!(v[0], 1.0);
    assert_double_eq!(v[1], 2.0);
    assert_double_eq!(v[2], 3.0);
}

#[test]
fn add() {
    let sum = Vector::new(1.0, 2.0, 3.0) + Vector::new(4.0, 5.0, 6.0);
    assert_double_eq!(sum.x, 5.0);
    assert_double_eq!(sum.y, 7.0);
    assert_double_eq!(sum.z, 9.0);
}

#[test]
fn subtract() {
    let diff = Vector::new(4.0, 5.0, 6.0) - Vector::new(1.0, 2.0, 3.0);
    assert_double_eq!(diff.x, 3.0);
    assert_double_eq!(diff.y, 3.0);
    assert_double_eq!(diff.z, 3.0);
}

#[test]
fn multiply() {
    let scaled = Vector::new(1.0, 2.0, 3.0) * 2.0;
    assert_double_eq!(scaled.x, 2.0);
    assert_double_eq!(scaled.y, 4.0);
    assert_double_eq!(scaled.z, 6.0);
}

#[test]
fn divide() {
    let divided = Vector::new(4.0, 6.0, 8.0) / 2.0;
    assert_double_eq!(divided.x, 2.0);
    assert_double_eq!(divided.y, 3.0);
    assert_double_eq!(divided.z, 4.0);
}

#[test]
fn cross() {
    let cross = Vector::new(1.0, 0.0, 0.0).cross(Vector::new(0.0, 1.0, 0.0));
    assert_near!(cross.x, 0.0, tolerance::LINEAR);
    assert_near!(cross.y, 0.0, tolerance::LINEAR);
    assert_near!(cross.z, 1.0, tolerance::LINEAR);
}

#[test]
fn dot() {
    // 1*4 + 2*5 + 3*6 = 32
    assert_double_eq!(Vector::new(1.0, 2.0, 3.0).dot(Vector::new(4.0, 5.0, 6.0)), 32.0);
}

#[test]
fn cross_magnitude() {
    let a = Vector::new(1.0, 0.0, 0.0);
    let b = Vector::new(0.0, 1.0, 0.0);
    assert_near!(a.cross(b).length(), 1.0, tolerance::LINEAR);
}

#[test]
fn cross_square_magnitude() {
    let a = Vector::new(1.0, 0.0, 0.0);
    let b = Vector::new(0.0, 1.0, 0.0);
    assert_near!(a.cross(b).length_squared(), 1.0, tolerance::LINEAR);
}

#[test]
fn dot_cross() {
    let a = Vector::new(1.0, 0.0, 0.0);
    let b = Vector::new(0.0, 1.0, 0.0);
    let c = Vector::new(0.0, 0.0, 1.0);
    assert_near!(a.dot(b.cross(c)), 1.0, tolerance::LINEAR);
}

#[test]
fn modulus() {
    assert_near!(Vector::new(3.0, 4.0, 0.0).length(), 5.0, tolerance::LINEAR);
}

#[test]
fn square_modulus() {
    assert_near!(Vector::new(3.0, 4.0, 0.0).length_squared(), 25.0, tolerance::LINEAR);
}

#[test]
fn set_linear_form() {
    let a = Vector::new(1.0, 0.0, 0.0);
    let b = Vector::new(0.0, 1.0, 0.0);
    let result = 2.0 * a + 3.0 * b;
    assert_near!(result.x, 2.0, tolerance::LINEAR);
    assert_near!(result.y, 3.0, tolerance::LINEAR);
    assert_near!(result.z, 0.0, tolerance::LINEAR);
}

#[test]
fn is_equal() {
    let a = Vector::new(1.0, 2.0, 3.0);
    let b = Vector::new(1.0, 2.0, 3.0);
    let c = Vector::new(1.0, 2.0, 4.0);
    assert!(a.approx_eq(b, tolerance::LINEAR));
    assert!(!a.approx_eq(c, tolerance::LINEAR));
}
