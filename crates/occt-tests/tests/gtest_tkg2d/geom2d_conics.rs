//! Geom2d_Circle_Test.cxx, Geom2d_Ellipse_Test.cxx, Geom2d_Hyperbola_Test.cxx,
//! Geom2d_Parabola_Test.cxx.

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve2::{Circle2, Ellipse2, Hyperbola2, Line2, Parabola2};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction2, Frame2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

use crate::{assert_point2_near, assert_vector2_near, p2};

const TOL: f64 = tolerance::LINEAR;

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= TOL
}

fn world() -> Frame2 {
    Frame2::new(p2(0.0, 0.0), Direction2::X, true)
}

fn translation(x: f64, y: f64) -> Transform2 {
    Transform2::translation(Vector2::new(x, y))
}

// ---- circle ----

fn circle() -> Circle2 {
    Circle2::new(world(), 5.0)
}

#[test]
fn circle_construction() {
    assert!(near(circle().radius(), 5.0));
    let c = Circle2::new(Frame2::new(p2(1.0, 2.0), Direction2::X, true), 3.0);
    assert!(near(c.radius(), 3.0));
    assert_eq!(c.center(), p2(1.0, 2.0));
    let c = Circle2::new(Frame2::from_axes(p2(3.0, 4.0), Direction2::X, Direction2::Y).unwrap(), 7.0);
    assert!(near(c.radius(), 7.0));
    assert_eq!(c.center(), p2(3.0, 4.0));
}

#[test]
fn circle_domain() {
    assert_eq!(circle().first_parameter(), 0.0);
    assert_eq!(circle().last_parameter(), 2.0 * PI);
    assert!(circle().is_closed());
    assert!(circle().is_periodic());
    assert_eq!(circle().period(), Some(2.0 * PI));
}

#[test]
fn circle_set_radius() {
    let mut c = circle();
    c.set_radius(10.0);
    assert!(near(c.radius(), 10.0));
}

#[test]
fn circle_evaluation() {
    let c = circle();
    assert_point2_near(c.point(0.0), p2(5.0, 0.0), TOL);
    assert_point2_near(c.point(FRAC_PI_2), p2(0.0, 5.0), TOL);
    assert_point2_near(c.point(PI), p2(-5.0, 0.0), TOL);
    let (_, d1) = c.d1(0.0);
    assert_vector2_near(d1, Vector2::new(0.0, 5.0), TOL);
    let (_, _, d2) = c.d2(0.0);
    assert_vector2_near(d2, Vector2::new(-5.0, 0.0), TOL);
    assert_vector2_near(c.dn(0.0, 1), Vector2::new(0.0, 5.0), TOL);
    let mut u = 0.0;
    while u < 2.0 * PI {
        assert!(near(c.point(u).distance(c.center()), 5.0));
        u += PI / 6.0;
    }
}

#[test]
fn circle_reversed_parameter() {
    assert!(near(circle().reversed_parameter(PI / 3.0), 2.0 * PI - PI / 3.0));
}

#[test]
fn circle_transform() {
    let c = circle().transformed(&translation(10.0, 20.0));
    assert_point2_near(c.center(), p2(10.0, 20.0), TOL);
    assert!(near(c.radius(), 5.0));
}

// ---- ellipse ----

fn ellipse() -> Ellipse2 {
    Ellipse2::new(world(), 10.0, 5.0)
}

#[test]
fn ellipse_radii_and_domain() {
    let e = ellipse();
    assert!(near(e.major_radius(), 10.0) && near(e.minor_radius(), 5.0));
    assert_eq!(e.first_parameter(), 0.0);
    assert_eq!(e.last_parameter(), 2.0 * PI);
    assert!(e.is_closed() && e.is_periodic());
}

#[test]
fn ellipse_focal_properties() {
    let e = ellipse();
    assert!(near(e.eccentricity(), (1.0_f64 - 25.0 / 100.0).sqrt()));
    assert!(near(e.focal_distance(), 2.0 * 75.0_f64.sqrt()));
    let (f1, f2) = e.foci();
    let c = 75.0_f64.sqrt();
    assert_point2_near(f1, p2(c, 0.0), TOL);
    assert_point2_near(f2, p2(-c, 0.0), TOL);
    assert!(near(e.semi_latus_rectum(), 2.5));
    let mut u = 0.0;
    while u < 2.0 * PI {
        let p = e.point(u);
        assert!(near(p.distance(f1) + p.distance(f2), 20.0));
        u += PI / 6.0;
    }
}

#[test]
fn ellipse_evaluation() {
    let e = ellipse();
    assert_point2_near(e.point(0.0), p2(10.0, 0.0), TOL);
    assert_point2_near(e.point(FRAC_PI_2), p2(0.0, 5.0), TOL);
    assert_vector2_near(e.d1(0.0).1, Vector2::new(0.0, 5.0), TOL);
    assert_vector2_near(e.d2(0.0).2, Vector2::new(-10.0, 0.0), TOL);
}

#[test]
fn ellipse_setters_transform_and_reverse() {
    let mut e = ellipse();
    e.set_major_radius(20.0);
    assert!(near(e.major_radius(), 20.0));
    e.set_minor_radius(3.0);
    assert!(near(e.minor_radius(), 3.0));
    let moved = ellipse().transformed(&translation(5.0, 10.0));
    assert_point2_near(moved.center(), p2(5.0, 10.0), TOL);
    assert!(near(moved.major_radius(), 10.0));
    assert!(near(ellipse().reversed_parameter(PI / 4.0), 2.0 * PI - PI / 4.0));
}

// ---- hyperbola ----

fn hyperbola() -> Hyperbola2 {
    Hyperbola2::new(world(), 5.0, 3.0)
}

#[test]
fn hyperbola_properties() {
    let h = hyperbola();
    assert!(near(h.major_radius(), 5.0) && near(h.minor_radius(), 3.0));
    assert!(!h.is_closed() && !h.is_periodic());
    assert!(near(h.eccentricity(), 34.0_f64.sqrt() / 5.0));
    assert!(near(h.focal_distance(), 2.0 * 34.0_f64.sqrt()));
    let (f1, f2) = h.foci();
    assert!(near(f1.x, 34.0_f64.sqrt()) && near(f2.x, -(34.0_f64.sqrt())));
    assert!(near(h.semi_latus_rectum(), 1.8));
    let (a1, a2) = h.asymptotes();
    assert!(near(a1.origin.x, 0.0) && near(a2.origin.x, 0.0));
}

#[test]
fn hyperbola_evaluation() {
    let h = hyperbola();
    assert_point2_near(h.point(0.0), p2(5.0, 0.0), TOL);
    assert_vector2_near(h.d1(0.0).1, Vector2::new(0.0, 3.0), TOL);
    let (f1, f2) = h.foci();
    let mut u = -2.0;
    while u <= 2.0 {
        let p = h.point(u);
        assert!(near((p.distance(f1) - p.distance(f2)).abs(), 10.0));
        u += 0.5;
    }
}

#[test]
fn hyperbola_setters_and_transform() {
    let mut h = hyperbola();
    h.set_major_radius(8.0);
    h.set_minor_radius(4.0);
    assert!(near(h.major_radius(), 8.0) && near(h.minor_radius(), 4.0));
    let moved = hyperbola().transformed(&translation(10.0, 20.0));
    assert_point2_near(moved.center(), p2(10.0, 20.0), TOL);
}

// ---- parabola ----

fn parabola() -> Parabola2 {
    Parabola2::new(world(), 4.0)
}

#[test]
fn parabola_properties() {
    let p = parabola();
    assert!(near(p.focal(), 4.0));
    assert!(!p.is_closed() && !p.is_periodic());
    assert_point2_near(p.focus(), p2(4.0, 0.0), TOL);
    assert!(near(p.semi_latus_rectum(), 8.0));
}

#[test]
fn parabola_evaluation() {
    let p = parabola();
    assert_point2_near(p.point(0.0), p2(0.0, 0.0), TOL);
    let (a, b) = (p.point(3.0), p.point(-3.0));
    assert!(near(a.x, b.x) && near(a.y, -b.y));
    let (_, d1) = p.d1(0.0);
    assert!(near(d1.x, 0.0) && d1.length() > TOL);
    let (_, _, d2) = p.d2(0.0);
    assert!(d2.x > 0.0);
}

#[test]
fn parabola_focus_directrix_property() {
    let p = parabola();
    let directrix = Line2::from(p.directrix());
    let mut u = -5.0;
    while u <= 5.0 {
        let q = p.point(u);
        assert!((q.distance(p.focus()) - directrix.distance_to_point(q)).abs() <= 1.0e-6);
        u += 1.0;
    }
}

#[test]
fn parabola_set_focal_and_transform() {
    let mut p = parabola();
    p.set_focal(10.0);
    assert!(near(p.focal(), 10.0));
    let moved = parabola().transformed(&translation(5.0, 10.0));
    assert_point2_near(moved.apex(), p2(5.0, 10.0), TOL);
}
