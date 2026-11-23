use std::f64::consts::PI;

use super::*;
use crate::math::{Direction, Frame};
use crate::primitives::{Cuboid, Cylinder, Sphere, Torus};
use crate::{properties, validity};

fn measure(model: &Model, shape: Shape) -> (f64, f64) {
    (properties::volume(model, shape).mass, properties::area(model, shape).mass)
}

fn assert_close(what: &str, value: f64, expected: f64, relative: f64) {
    assert!((value - expected).abs() <= relative * expected.abs(), "{what} {value} against {expected} ({:e} off)", (value - expected) / expected);
}

fn check(model: &Model, shape: Shape) {
    let report = validity::check(model, shape);
    assert!(report.is_valid(), "{report:?}");
}

/// A ball slid is a capsule: a cylinder between two half balls.
#[test]
fn a_ball_slid_is_a_capsule() {
    let mut model = Model::new();
    let ball = Sphere::new(1.5).with_frame(Frame::from_z(Point::new(1.0, 2.0, 3.0), Direction::Z)).build(&mut model).unwrap();
    let by = Vector::new(3.0, -2.0, 4.0);
    let swept = solid(&mut model, ball, &RigidMotion::Slide { by }).unwrap();
    let (volume, area) = measure(&model, swept);
    let (r, l) = (1.5f64, by.length());
    assert_close("volume", volume, 4.0 / 3.0 * PI * r.powi(3) + PI * r * r * l, 1e-8);
    assert_close("area", area, 4.0 * PI * r * r + 2.0 * PI * r * l, 1e-8);
    check(&model, swept);
}

/// A box slid aslant fills the box and, for each axis, its face square to it times the
/// slide's reach along it.
#[test]
fn a_box_slid_aslant() {
    let mut model = Model::new();
    let (a, b, c) = (2.0, 3.0, 1.5);
    let block = Cuboid::new(a, b, c).build(&mut model).unwrap();
    let by = Vector::new(1.2, -0.7, 2.5);
    let swept = solid(&mut model, block, &RigidMotion::Slide { by }).unwrap();
    let (volume, _) = measure(&model, swept);
    assert_close("volume", volume, a * b * c + by.x.abs() * b * c + by.y.abs() * a * c + by.z.abs() * a * b, 1e-10);
    check(&model, swept);
}

/// A torus slid along its axis: itself and, for the reach, the ring between its outer and
/// inner equators.
#[test]
fn a_torus_slid_along_its_axis() {
    let mut model = Model::new();
    let (big, small) = (3.0, 1.0);
    let axis = Direction::new(0.3, -0.2, 0.9);
    let ring = Torus::new(big, small).with_frame(Frame::from_z(Point::new(0.5, 0.0, -1.0), axis)).build(&mut model).unwrap();
    let length = 5.0;
    let swept = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::from(axis) * length }).unwrap();
    let (volume, area) = measure(&model, swept);
    assert_close("volume", volume, 2.0 * PI * PI * big * small * small + length * 4.0 * PI * big * small, 1e-8);
    assert_close("area", area, 4.0 * PI * PI * big * small + 4.0 * PI * big * length, 1e-8);
    check(&model, swept);
}

/// A drum turned about a line beside it and parallel to its axis sweeps the annular sector
/// between its nearest and farthest circles, and itself at each end.
#[test]
fn a_drum_turned_about_a_parallel_line() {
    let mut model = Model::new();
    let (r, h, big) = (1.0, 2.0, 4.0);
    let drum = Cylinder::new(r, h).with_frame(Frame::from_z(Point::new(big, 0.0, 0.0), Direction::Z)).build(&mut model).unwrap();
    let angle = 1.3;
    let swept = solid(&mut model, drum, &RigidMotion::Turn { axis: Axis::new(Point::ORIGIN, Direction::Z), from: -0.4, to: angle - 0.4 }).unwrap();
    let (volume, area) = measure(&model, swept);
    let section = PI * r * r + 2.0 * big * r * angle;
    assert_close("volume", volume, section * h, 1e-8);
    assert_close("area", area, 2.0 * section + h * (2.0 * PI * r + 2.0 * big * angle), 1e-8);
    check(&model, swept);
}

/// A ball screwed sweeps a tube round the helix its centre runs, closed by half balls:
/// `πr²·L + 4/3·πr³` (its turns apart, its curvature gentler than the ball).
#[test]
fn a_ball_screwed_is_a_helical_tube() {
    let mut model = Model::new();
    let (r, big, lead) = (0.8, 4.0, 3.0);
    let ball = Sphere::new(r).with_frame(Frame::from_z(Point::new(big, 0.0, 0.0), Direction::new(0.2, 0.3, 1.0))).build(&mut model).unwrap();
    let (from, to) = (-0.5, 4.0);
    let swept = solid(&mut model, ball, &RigidMotion::Screw { axis: Axis::new(Point::ORIGIN, Direction::Z), lead, from, to }).unwrap();
    let (volume, area) = measure(&model, swept);
    let length = (to - from) * (big * big + (lead / TAU).powi(2)).sqrt();
    assert_close("volume", volume, PI * r * r * length + 4.0 / 3.0 * PI * r.powi(3), 1e-7);
    assert_close("area", area, 2.0 * PI * r * length + 4.0 * PI * r * r, 1e-7);
    check(&model, swept);
}
