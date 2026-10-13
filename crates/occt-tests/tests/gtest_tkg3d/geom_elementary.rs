//! Geom_Line_Test.cxx, Geom_Circle_Test.cxx, Geom_Plane_Test.cxx, and the trimming checks
//! of GeomAdaptor_Curve_Test.cxx.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use kentron::geometry::curve::{Circle, Curve, Line, TrimmedCurve};
use kentron::geometry::surface::Plane;
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Axis, Direction, Frame, Transform, Transformable, Vector};
use kentron::tolerance;

use crate::{assert_point_near, assert_vector_near, p};

const TOL: f64 = tolerance::LINEAR;

// ---- line ----

#[test]
fn line_construction() {
    let line = Line::from(Axis::new(p(0.0, 0.0, 0.0), Direction::X));
    assert_eq!(line.origin, p(0.0, 0.0, 0.0));
    let line = Line::new(p(1.0, 2.0, 3.0), Direction::Z);
    assert_eq!(line.origin, p(1.0, 2.0, 3.0));
    assert!(line.direction.approx_eq(Direction::Z, tolerance::ANGULAR));
}

#[test]
fn line_evaluation() {
    assert_point_near(Line::new(p(0.0, 0.0, 0.0), Direction::X).point(5.0), p(5.0, 0.0, 0.0), TOL);
    let (q, d1) = Line::new(p(0.0, 0.0, 0.0), Direction::Y).d1(3.0);
    assert!((q.y - 3.0).abs() <= TOL);
    assert_vector_near(d1, Vector::new(0.0, 1.0, 0.0), TOL);
    assert!(Line::new(p(0.0, 0.0, 0.0), Direction::X).d2(1.0).2.length() <= TOL);
}

#[test]
fn line_infinite_parameters() {
    let line = Line::new(p(0.0, 0.0, 0.0), Direction::X);
    assert_eq!(line.first_parameter(), -tolerance::INFINITE);
    assert_eq!(line.last_parameter(), tolerance::INFINITE);
}

#[test]
fn line_reverse_and_transform() {
    let mut line = Line::new(p(0.0, 0.0, 0.0), Direction::X);
    line.reverse();
    assert!(line.direction.approx_eq(Direction::NEG_X, tolerance::ANGULAR));
    let moved = Line::new(p(0.0, 0.0, 0.0), Direction::X).transformed(&Transform::translation(Vector::new(0.0, 0.0, 5.0)));
    assert!((moved.origin.z - 5.0).abs() <= TOL);
}

// ---- circle ----

fn unit_circle() -> Circle {
    Circle::new(Frame::WORLD, 1.0)
}

#[test]
fn circle_construction() {
    assert!((Circle::new(Frame::WORLD, 5.0).radius() - 5.0).abs() <= TOL);
    let c = Circle::new(Frame::from_z(p(1.0, 2.0, 3.0), Direction::Z), 5.0);
    assert_point_near(c.center(), p(1.0, 2.0, 3.0), TOL);
}

#[test]
fn circle_evaluation() {
    let c = unit_circle();
    assert_point_near(c.point(0.0), p(1.0, 0.0, 0.0), TOL);
    assert_point_near(c.point(FRAC_PI_2), p(0.0, 1.0, 0.0), TOL);
    assert_point_near(c.point(PI), p(-1.0, 0.0, 0.0), TOL);
    assert_point_near(c.point(3.0 * FRAC_PI_2), p(0.0, -1.0, 0.0), TOL);
    assert_vector_near(c.d1(0.0).1, Vector::new(0.0, 1.0, 0.0), TOL);
    assert_vector_near(c.d2(0.0).2, Vector::new(-1.0, 0.0, 0.0), TOL);
}

#[test]
fn circle_reversed_parameter_and_scaling() {
    assert!((unit_circle().reversed_parameter(FRAC_PI_4) - (2.0 * PI - FRAC_PI_4)).abs() <= TOL);
    let scaled = Circle::new(Frame::WORLD, 5.0).transformed(&Transform::scaling(p(0.0, 0.0, 0.0), 2.0));
    assert!((scaled.radius() - 10.0).abs() <= TOL);
}

// ---- plane ----

#[test]
fn plane_construction_and_evaluation() {
    let plane = Plane::new(p(1.0, 2.0, 3.0), Direction::Z);
    assert_point_near(plane.origin(), p(1.0, 2.0, 3.0), TOL);
    let xy = Plane::new(p(0.0, 0.0, 0.0), Direction::Z);
    assert_point_near(xy.point(3.0, 4.0), p(3.0, 4.0, 0.0), TOL);
}

#[test]
fn plane_equation() {
    let (a, b, c, d) = Plane::new(p(0.0, 0.0, 5.0), Direction::Z).equation();
    assert!(a.abs() <= TOL && b.abs() <= TOL && (c - 1.0).abs() <= TOL && (d + 5.0).abs() <= TOL);
}

#[test]
fn plane_isolines_are_lines() {
    let plane = Plane::new(p(0.0, 0.0, 0.0), Direction::Z);
    let Curve::Line(u_iso) = plane.u_iso(2.0) else { panic!("u-iso of a plane is not a line") };
    assert!((u_iso.point(0.0).x - 2.0).abs() <= TOL);
    let Curve::Line(v_iso) = plane.v_iso(3.0) else { panic!("v-iso of a plane is not a line") };
    assert!((v_iso.point(0.0).y - 3.0).abs() <= TOL);
}

#[test]
fn plane_transform() {
    let moved = Plane::new(p(0.0, 0.0, 0.0), Direction::Z).transformed(&Transform::translation(Vector::new(0.0, 0.0, 10.0)));
    assert!((moved.origin().z - 10.0).abs() <= TOL);
}

// ---- trimming (GeomAdaptor_Curve_Test) ----

fn x_line() -> Curve {
    Curve::Line(Line::new(p(0.0, 0.0, 0.0), Direction::X))
}

#[test]
fn trimming_validation() {
    TrimmedCurve::new(x_line(), 0.0, 1.0);
    TrimmedCurve::new(x_line(), 1.0, 1.0);
    TrimmedCurve::new(x_line(), 1.0 + 0.5 * TOL, 1.0);
    TrimmedCurve::new(x_line(), 0.0, 1.0e-9);
    occt_tests::assert_panics!(TrimmedCurve::new(x_line(), 2.0, 1.0));
    occt_tests::assert_panics!(TrimmedCurve::new(x_line(), 1.0 + 10.0 * TOL, 1.0));
}

#[test]
fn degenerate_curves_can_be_trimmed() {
    for center in [p(0.0, 0.0, 0.0), p(5.0, -3.0, 2.0)] {
        let point_circle = Curve::Circle(Circle::new(Frame::from_z(center, Direction::Z), 0.0));
        let t = TrimmedCurve::new(point_circle, 0.0, PI);
        assert_point_near(t.point(1.0), center, TOL);
    }
}
