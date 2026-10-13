//! Geom2d_TrimmedCurve_Test.cxx, plus the trimming-range checks of
//! Geom2dAdaptor_Curve_Test.cxx (kentron trims through `TrimmedCurve2`, not adaptors).

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve2::{Circle2, Curve2, Line2, TrimmedCurve2};
use kentron::geometry::{Continuity, ParametricCurve};
use kentron::math::{Direction2, Frame2, Transform2, Transformable2, Vector2};
use kentron::tolerance;

use crate::{assert_point2_near, p2};

const TOL: f64 = tolerance::LINEAR;

fn line() -> Curve2 {
    Curve2::Line(Line2::new(p2(0.0, 0.0), Direction2::X))
}

fn circle() -> Curve2 {
    Curve2::Circle(Circle2::new(Frame2::WORLD, 5.0))
}

#[test]
fn construct_from_line_and_circle() {
    let t = TrimmedCurve2::new(line(), 0.0, 10.0);
    assert_eq!((t.first_parameter(), t.last_parameter()), (0.0, 10.0));
    let t = TrimmedCurve2::new(circle(), 0.0, PI);
    assert_eq!((t.first_parameter(), t.last_parameter()), (0.0, PI));
    assert_eq!(t.basis(), &circle());
}

#[test]
fn start_and_end_points() {
    let t = TrimmedCurve2::new(line(), 2.0, 8.0);
    assert_point2_near(t.start_point(), p2(2.0, 0.0), TOL);
    assert_point2_near(t.end_point(), p2(8.0, 0.0), TOL);
    let arc = TrimmedCurve2::new(circle(), 0.0, FRAC_PI_2);
    assert_point2_near(arc.start_point(), p2(5.0, 0.0), TOL);
    assert_point2_near(arc.end_point(), p2(0.0, 5.0), TOL);
}

#[test]
fn evaluation_matches_basis() {
    let arc = TrimmedCurve2::new(circle(), 0.0, PI);
    let mut u = 0.0;
    while u <= PI {
        assert_point2_near(arc.point(u), circle().point(u), TOL);
        u += PI / 4.0;
    }
    let segment = TrimmedCurve2::new(line(), 0.0, 10.0);
    let ((p, d), (bp, bd)) = (segment.d1(5.0), line().d1(5.0));
    assert_point2_near(p, bp, TOL);
    assert!((d - bd).length() <= TOL);
}

#[test]
fn continuity() {
    let t = TrimmedCurve2::new(line(), 0.0, 10.0);
    assert_eq!(t.continuity(), Continuity::CN);
    assert!(t.is_cn(10));
}

#[test]
fn closed_and_periodic() {
    let half = TrimmedCurve2::new(circle(), 0.0, PI);
    assert!(!half.is_closed());
    assert!(!half.is_periodic());
    let full = TrimmedCurve2::new(circle(), 0.0, 2.0 * PI);
    assert!(full.is_closed());
    assert!(full.is_periodic());
    assert!(!TrimmedCurve2::new(line(), 0.0, 10.0).is_periodic());
}

#[test]
fn set_trim() {
    let mut t = TrimmedCurve2::new(line(), 0.0, 10.0);
    t.set_trim(3.0, 7.0);
    assert_eq!((t.first_parameter(), t.last_parameter()), (3.0, 7.0));
}

#[test]
fn transform() {
    let t = TrimmedCurve2::new(line(), 0.0, 10.0).transformed(&Transform2::translation(Vector2::new(5.0, 5.0)));
    assert_point2_near(t.start_point(), p2(5.0, 5.0), TOL);
}

#[test]
fn trimming_validation() {
    // Valid, equal, and slightly inverted (within LINEAR) bounds are accepted.
    TrimmedCurve2::new(line(), 0.0, 1.0);
    TrimmedCurve2::new(line(), 1.0, 1.0);
    TrimmedCurve2::new(line(), 1.0 + 0.5 * TOL, 1.0);
    // Clearly inverted bounds are rejected.
    occt_tests::assert_panics!(TrimmedCurve2::new(line(), 2.0, 1.0));
    occt_tests::assert_panics!(TrimmedCurve2::new(line(), 1.0 + 10.0 * TOL, 1.0));
}

#[test]
fn degenerate_curves_can_be_trimmed() {
    let point_circle = Curve2::Circle(Circle2::new(Frame2::WORLD, 0.0));
    let t = TrimmedCurve2::new(point_circle, 0.0, PI);
    assert_point2_near(t.point(1.0), p2(0.0, 0.0), TOL);
}
