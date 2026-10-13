//! Geom2d_Transformation_Test.cxx: 2D transformations (`math::Transform2`).
//!
//! Not ported: `Form()` classification checks (kentron has no transformation "form" enum;
//! `is_identity` covers the identity case), `Copy`.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use kentron::math::{Axis2, Direction2, Point2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

use crate::{assert_point2_near, p2};

const TOL: f64 = tolerance::LINEAR;

fn apply(t: &Transform2, x: f64, y: f64) -> Point2 {
    p2(x, y).transformed(t)
}

#[test]
fn identity() {
    let t = Transform2::IDENTITY;
    assert!(t.is_identity());
    assert!((t.scale_factor() - 1.0).abs() <= TOL);
    assert!(!t.is_negative());
    assert_eq!(t.matrix(), [[1.0, 0.0], [0.0, 1.0]]);
}

#[test]
fn translation() {
    let t = Transform2::translation(Vector2::new(5.0, 10.0));
    assert_point2_near(apply(&t, 1.0, 2.0), p2(6.0, 12.0), TOL);
    let between = Transform2::translation(p2(4.0, 5.0) - p2(1.0, 1.0));
    assert_point2_near(apply(&between, 0.0, 0.0), p2(3.0, 4.0), TOL);
    assert!(!t.is_identity());
}

#[test]
fn rotation() {
    let quarter = Transform2::rotation(p2(0.0, 0.0), FRAC_PI_2);
    assert_point2_near(apply(&quarter, 1.0, 0.0), p2(0.0, 1.0), TOL);
    let half = Transform2::rotation(p2(0.0, 0.0), PI);
    assert_point2_near(apply(&half, 3.0, 4.0), p2(-3.0, -4.0), TOL);
    let about = Transform2::rotation(p2(1.0, 1.0), FRAC_PI_2);
    assert_point2_near(apply(&about, 1.0, 1.0), p2(1.0, 1.0), TOL);
    assert_point2_near(apply(&about, 2.0, 1.0), p2(1.0, 2.0), TOL);
}

#[test]
fn scaling() {
    let t = Transform2::scaling(p2(0.0, 0.0), 2.0);
    assert!((t.scale_factor() - 2.0).abs() <= TOL);
    assert_point2_near(apply(&t, 3.0, 4.0), p2(6.0, 8.0), TOL);
    let negative = Transform2::scaling(p2(0.0, 0.0), -1.0);
    assert!(!negative.is_negative());
    assert!((negative.scale_factor() + 1.0).abs() <= TOL);
    assert_point2_near(apply(&negative, 3.0, 4.0), p2(-3.0, -4.0), TOL);
    let about = Transform2::scaling(p2(2.0, 3.0), 2.0);
    assert_point2_near(apply(&about, 2.0, 3.0), p2(2.0, 3.0), TOL);
    assert_point2_near(apply(&about, 4.0, 3.0), p2(6.0, 3.0), TOL);
}

#[test]
fn mirrors() {
    let point = Transform2::mirror_point(p2(0.0, 0.0));
    assert!(!point.is_negative());
    assert_point2_near(apply(&point, 3.0, 4.0), p2(-3.0, -4.0), TOL);
    let axis = Transform2::mirror_axis(Axis2::new(p2(0.0, 0.0), Direction2::X));
    assert!(axis.is_negative());
    assert_point2_near(apply(&axis, 3.0, 4.0), p2(3.0, -4.0), TOL);
}

#[test]
fn inverse() {
    let t = Transform2::translation(Vector2::new(5.0, 10.0));
    let back = apply(&t, 7.0, 3.0).transformed(&t.inverse());
    assert_point2_near(back, p2(7.0, 3.0), TOL);
}

#[test]
fn composition() {
    let translate = Transform2::translation(Vector2::new(1.0, 0.0));
    let rotate = Transform2::rotation(p2(0.0, 0.0), FRAC_PI_2);
    assert_point2_near(apply(&(rotate * translate), 0.0, 0.0), p2(0.0, 1.0), TOL);
}

#[test]
fn powers() {
    assert!(Transform2::rotation(p2(0.0, 0.0), FRAC_PI_4).powi(0).is_identity());
    let twice = Transform2::rotation(p2(0.0, 0.0), FRAC_PI_2).powi(2);
    assert_point2_near(apply(&twice, 1.0, 0.0), p2(-1.0, 0.0), TOL);
}
