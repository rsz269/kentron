//! gp_Pnt2d_Test.cxx

use std::f64::consts::FRAC_PI_2;

use occt_tests::{assert_double_eq, assert_near};
use kentron::math::{Axis2, Direction2, Point2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

fn assert_point2(p: Point2, x: f64, y: f64) {
    assert_near!(p.x, x, tolerance::LINEAR);
    assert_near!(p.y, y, tolerance::LINEAR);
}

#[test]
fn default_constructor() {
    let p = Point2::default();
    assert_double_eq!(p.x, 0.0);
    assert_double_eq!(p.y, 0.0);
}

#[test]
fn coordinate_constructor() {
    let p = Point2::new(1.0, 2.0);
    assert_double_eq!(p.x, 1.0);
    assert_double_eq!(p.y, 2.0);
}

#[test]
fn xy_constructor() {
    let p = Point2::ORIGIN + Vector2::new(3.0, 4.0);
    assert_double_eq!(p.x, 3.0);
    assert_double_eq!(p.y, 4.0);
}

#[test]
fn set_coord() {
    let mut p = Point2::default();
    p.x = 10.0;
    p.y = 20.0;
    assert_double_eq!(p.x, 10.0);
    assert_double_eq!(p.y, 20.0);

    p = Point2::new(7.0, 8.0);
    assert_double_eq!(p.x, 7.0);
    assert_double_eq!(p.y, 8.0);
}

#[test]
fn distance() {
    assert_near!(Point2::new(0.0, 0.0).distance(Point2::new(3.0, 4.0)), 5.0, tolerance::LINEAR);
}

#[test]
fn square_distance() {
    // (4-1)^2 + (6-2)^2 = 9 + 16 = 25
    let d = Point2::new(1.0, 2.0).distance_squared(Point2::new(4.0, 6.0));
    assert_near!(d, 25.0, tolerance::LINEAR);
}

#[test]
fn is_equal() {
    let p1 = Point2::new(1.0, 2.0);
    assert!(p1.approx_eq(Point2::new(1.0, 2.0), tolerance::LINEAR));
    assert!(!p1.approx_eq(Point2::new(1.0, 3.0), tolerance::LINEAR));
}

#[test]
fn translate_by_vec2d() {
    assert_point2(Point2::new(1.0, 2.0).translated(Vector2::new(10.0, 20.0)), 11.0, 22.0);
}

#[test]
fn scale() {
    assert_point2(Point2::new(2.0, 4.0).scaled(Point2::ORIGIN, 2.0), 4.0, 8.0);
}

#[test]
fn rotate() {
    assert_point2(Point2::new(1.0, 0.0).rotated(Point2::ORIGIN, FRAC_PI_2), 0.0, 1.0);
}

#[test]
fn mirror_point() {
    assert_point2(Point2::new(1.0, 0.0).mirrored_point(Point2::ORIGIN), -1.0, 0.0);
}

#[test]
fn mirror_axis() {
    // mirror about the X axis
    assert_point2(Point2::new(1.0, 1.0).mirrored_axis(Axis2::X), 1.0, -1.0);
}

#[test]
fn transformed() {
    let t = Transform2::translation(Vector2::new(0.0, 5.0));
    assert_point2(Point2::new(1.0, 0.0).transformed(&t), 1.0, 5.0);
}

/// OCC22736: mirroring twice about the same axis is the identity, whether applied in turn
/// or composed.
#[test]
fn occ22736_two_identical_mirror_composition_is_identity() {
    let first = Point2::new(2.0, 1.0);
    let second = Point2::new(3.0, 1.0);
    let axis = Axis2::new(first, Direction2::from_vector(second - first));

    let p1 = Point2::new(1.0, 0.0);
    let p2 = Point2::new(1.0, 2.0);

    let m1 = Transform2::mirror_axis(axis);
    let m2 = Transform2::mirror_axis(axis);

    // M1 maps p1 to p2
    let p1_m1 = p1.transformed(&m1);
    assert_near!(p1_m1.x, p2.x, tolerance::LINEAR);
    assert_near!(p1_m1.y, p2.y, tolerance::LINEAR);

    // M2 maps it back to p1
    let p1_m1_m2 = p1_m1.transformed(&m2);
    assert_near!(p1_m1_m2.x, p1.x, tolerance::LINEAR);
    assert_near!(p1_m1_m2.y, p1.y, tolerance::LINEAR);

    // M2*M1 is the identity on p1
    let p1_comp = p1.transformed(&(m2 * m1));
    assert_near!(p1_comp.x, p1.x, tolerance::LINEAR);
    assert_near!(p1_comp.y, p1.y, tolerance::LINEAR);
}
