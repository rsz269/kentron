//! gp_Pnt_Test.cxx

use std::f64::consts::FRAC_PI_2;

use occt_tests::{assert_double_eq, assert_near};
use kentron::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use kentron::tolerance;

fn assert_point(p: Point, x: f64, y: f64, z: f64) {
    assert_near!(p.x, x, tolerance::LINEAR);
    assert_near!(p.y, y, tolerance::LINEAR);
    assert_near!(p.z, z, tolerance::LINEAR);
}

#[test]
fn default_constructor() {
    let p = Point::default();
    assert_double_eq!(p.x, 0.0);
    assert_double_eq!(p.y, 0.0);
    assert_double_eq!(p.z, 0.0);
}

#[test]
fn coordinate_constructor() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_double_eq!(p.x, 1.0);
    assert_double_eq!(p.y, 2.0);
    assert_double_eq!(p.z, 3.0);
}

#[test]
fn xyz_constructor() {
    let p = Point::ORIGIN + Vector::new(4.0, 5.0, 6.0);
    assert_double_eq!(p.x, 4.0);
    assert_double_eq!(p.y, 5.0);
    assert_double_eq!(p.z, 6.0);
}

#[test]
fn set_coord() {
    let mut p = Point::default();
    p.x = 10.0;
    p.y = 20.0;
    p.z = 30.0;
    assert_double_eq!(p.x, 10.0);
    assert_double_eq!(p.y, 20.0);
    assert_double_eq!(p.z, 30.0);

    p = Point::new(7.0, 8.0, 9.0);
    assert_double_eq!(p.x, 7.0);
    assert_double_eq!(p.y, 8.0);
    assert_double_eq!(p.z, 9.0);
}

#[test]
fn distance() {
    let p1 = Point::new(0.0, 0.0, 0.0);
    let p2 = Point::new(3.0, 4.0, 0.0);
    assert_near!(p1.distance(p2), 5.0, tolerance::LINEAR);
}

#[test]
fn square_distance() {
    let p1 = Point::new(1.0, 2.0, 3.0);
    let p2 = Point::new(4.0, 6.0, 3.0);
    // (4-1)^2 + (6-2)^2 + (3-3)^2 = 9 + 16 + 0 = 25
    assert_near!(p1.distance_squared(p2), 25.0, tolerance::LINEAR);
}

#[test]
fn is_equal() {
    let p1 = Point::new(1.0, 2.0, 3.0);
    let p2 = Point::new(1.0, 2.0, 3.0);
    let p3 = Point::new(1.0, 2.0, 4.0);
    assert!(p1.approx_eq(p2, tolerance::LINEAR));
    assert!(!p1.approx_eq(p3, tolerance::LINEAR));
}

#[test]
fn bary_center() {
    // equal weights: midpoint
    let p = Point::barycenter(&[(Point::new(0.0, 0.0, 0.0), 1.0), (Point::new(10.0, 0.0, 0.0), 1.0)]);
    assert_point(p, 5.0, 0.0, 0.0);
}

#[test]
fn translate_by_vec() {
    let result = Point::new(1.0, 2.0, 3.0).translated(Vector::new(10.0, 20.0, 30.0));
    assert_point(result, 11.0, 22.0, 33.0);
}

#[test]
fn translate_by_two_points() {
    let from = Point::new(0.0, 0.0, 0.0);
    let to = Point::new(5.0, 5.0, 5.0);
    let result = Point::new(1.0, 2.0, 3.0).translated(to - from);
    assert_point(result, 6.0, 7.0, 8.0);
}

#[test]
fn scale() {
    let result = Point::new(2.0, 4.0, 6.0).scaled(Point::ORIGIN, 2.0);
    assert_point(result, 4.0, 8.0, 12.0);
}

#[test]
fn rotate() {
    let result = Point::new(1.0, 0.0, 0.0).rotated(Axis::Z, FRAC_PI_2);
    assert_point(result, 0.0, 1.0, 0.0);
}

#[test]
fn mirror_point() {
    let result = Point::new(1.0, 0.0, 0.0).mirrored_point(Point::ORIGIN);
    assert_point(result, -1.0, 0.0, 0.0);
}

#[test]
fn mirror_axis() {
    let result = Point::new(1.0, 1.0, 0.0).mirrored_axis(Axis::X);
    assert_point(result, 1.0, -1.0, 0.0);
}

#[test]
fn mirror_plane() {
    // mirror through the XY plane (Z=0)
    let result = Point::new(1.0, 2.0, 3.0).mirrored_plane(Point::ORIGIN, Direction::Z);
    assert_point(result, 1.0, 2.0, -3.0);
}

#[test]
fn transformed() {
    let t = Transform::translation(Vector::new(0.0, 0.0, 5.0));
    let result = Point::new(1.0, 0.0, 0.0).transformed(&t);
    assert_point(result, 1.0, 0.0, 5.0);
}
