use std::f64::consts::PI;

use super::*;
use crate::math::{Direction, Frame, Transform};
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

/// A ball carried straight, ever faster: still a capsule (the motion given by its places).
#[test]
fn a_ball_carried_straight_unevenly_is_a_capsule() {
    let mut model = Model::new();
    let ball = Sphere::new(1.0).build(&mut model).unwrap();
    let by = Vector::new(2.0, 1.0, -2.0);
    let carried = RigidMotion::Carried(Carriage(std::sync::Arc::new(move |t: f64| Transform::translation(by * (t * t)))));
    let swept = solid(&mut model, ball, &carried).unwrap();
    let (volume, area) = measure(&model, swept);
    assert_close("volume", volume, 4.0 / 3.0 * PI + PI * 3.0, 1e-8);
    assert_close("area", area, 4.0 * PI + 2.0 * PI * 3.0, 1e-8);
    check(&model, swept);
}

/// A bar carried by a screw given by its places sweeps what the screw sweeps it through.
#[test]
fn a_screw_given_by_its_places_sweeps_as_the_screw() {
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let (lead, from, to) = (5.0, -0.3, 2.9);
    let volumes: Vec<f64> = (0..2)
        .map(|k| {
            let mut model = Model::new();
            let bar = Cuboid::new(3.0, 1.0, 1.5).with_frame(Frame::from_z(Point::new(4.0, -0.5, 0.0), Direction::Z)).build(&mut model).unwrap();
            let screw = RigidMotion::Screw { axis, lead, from, to };
            let motion = if k == 0 { screw } else { RigidMotion::Carried(Carriage(std::sync::Arc::new(move |t: f64| screw.at(t)))) };
            let swept = solid(&mut model, bar, &motion).unwrap_or_else(|e| panic!("motion {k}: {e:?}"));
            check(&model, swept);
            measure(&model, swept).0
        })
        .collect();
    assert_close("volume", volumes[1], volumes[0], 1e-8);
}

/// A box carried round a circle, turning with it, runs along its top and bottom beside faces
/// it carries out: refused; as is a ring whose contact moves over it (tumbling about one line
/// while carried along another).
#[test]
fn motions_beyond_the_product_are_refused() {
    let mut model = Model::new();
    let key = Cuboid::new(3.0, 1.6, 4.0).build(&mut model).unwrap();
    let round = RigidMotion::Carried(Carriage(std::sync::Arc::new(|t: f64| Transform::rotation(Axis::new(Point::new(0.0, 10.0, 0.0), Direction::Z), t))));
    assert!(solid(&mut model, key, &round).is_err());
    let torus = Torus::new(3.0, 1.0).build(&mut model).unwrap();
    let tumble = RigidMotion::Carried(Carriage(std::sync::Arc::new(|t: f64| {
        Transform::translation(Vector::new(0.0, 10.0 * t, 0.0)) * Transform::rotation(Axis::new(Point::ORIGIN, Direction::X), 1.2 * t)
    })));
    assert!(solid(&mut model, torus, &tumble).is_err());
}

/// A torus leaning 35° pushed straight down 8, against the volume its slide covers summed
/// column by column (each vertical line meets the torus in intervals, each reaching 8 lower),
/// on a grid fine enough for 1e-3.
#[test]
fn a_leaning_torus_pushed_down_against_its_columns() {
    let mut model = Model::new();
    let lean = 35f64.to_radians();
    let a = Vector::new(lean.sin(), 0.0, lean.cos());
    let ring = Torus::new(3.0, 1.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, 4.0), Direction::from_vector(a))).build(&mut model).unwrap();
    let swept = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::new(0.0, 0.0, -8.0) }).unwrap();
    check(&model, swept);
    let (volume, _) = measure(&model, swept);
    // (The torus as (|q⊥| − 3)² + (q·a)² ≤ 1, q from its centre.)
    let inside = |p: Point| {
        let q = p - Point::new(0.0, 0.0, 4.0);
        let along = q.dot(a);
        let across = (q - a * along).length();
        (across - 3.0).powi(2) + along * along - 1.0
    };
    let n = 300;
    let h = 9.0 / n as f64;
    let mut total = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (-4.5 + h * (i as f64 + 0.5), -4.5 + h * (j as f64 + 0.5));
            let f = |z: f64| inside(Point::new(x, y, z));
            // The column's crossings, by sign changes on fine steps, each narrowed by bisection.
            let steps = 300;
            let z = |k: usize| -0.5 + 9.0 * k as f64 / steps as f64;
            let mut crossings = Vec::new();
            for k in 0..steps {
                let f0 = f(z(k));
                if (f0 < 0.0) != (f(z(k + 1)) < 0.0) {
                    let (mut l, mut r) = (z(k), z(k + 1));
                    for _ in 0..50 {
                        let m = 0.5 * (l + r);
                        if (f(m) < 0.0) == (f0 < 0.0) { l = m } else { r = m }
                    }
                    crossings.push(0.5 * (l + r));
                }
            }
            let mut covered: Vec<(f64, f64)> = crossings.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0] - 8.0, c[1])).collect();
            covered.sort_by(|p, q| p.0.total_cmp(&q.0));
            let mut end = f64::NEG_INFINITY;
            for (lo, hi) in covered {
                let start = lo.max(end);
                if hi > start {
                    total += hi - start;
                }
                end = end.max(hi);
            }
        }
    }
    assert_close("volume", volume, total * h * h, 1e-3);
}
