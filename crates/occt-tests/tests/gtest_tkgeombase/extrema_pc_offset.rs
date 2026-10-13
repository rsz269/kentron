//! ExtremaPC_OffsetCurve_Test.cxx (closely related OCCT cases are grouped into one test).
//!
//! Not ported: `ScaledTransformedOffsetUsesParameterSafeFallback` (OCCT's transformed-curve
//! adaptor, which keeps the untransformed parametrization).

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve::{Circle, Curve, Ellipse, Line, OffsetCurve};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Point};

use crate::{assert_near, p};

const TOL: f64 = 1.0e-6;

fn offset(basis: Curve, distance: f64) -> Curve {
    Curve::Offset(OffsetCurve::new(basis, distance, Direction::Z))
}

fn circle(radius: f64) -> Curve {
    Curve::Circle(Circle::new(Frame::WORLD, radius))
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

#[track_caller]
fn assert_min_distance(c: &Curve, q: Point, d: f64) {
    let r = with_endpoints(c, q);
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), d, TOL);
}

#[test]
fn offset_circle() {
    let c = offset(circle(10.), 5.);
    assert_min_distance(&c, p(30., 0., 0.), 15.);
    assert_min_distance(&c, p(12., 0., 0.), 3.);
    assert_min_distance(&c, p(15., 0., 0.), 0.);
    assert_min_distance(&c, p(15., 0., 10.), 10.);
    assert_min_distance(&offset(circle(10.), -5.), p(20., 0., 0.), 15.);
}

#[test]
fn offset_ellipse() {
    let c = offset(Curve::Ellipse(Ellipse::new(Frame::WORLD, 20., 10.)), 5.);
    assert_min_distance(&c, p(40., 0., 0.), 15.);
    assert_min_distance(&c, p(0., 25., 0.), 10.);
}

#[test]
fn offset_line_segment() {
    let c = offset(Curve::Line(Line::new(p(0., 0., 0.), Direction::X)), 5.);
    let on = c.point(50.);
    let r = PointCurveSearch::with_domain(&c, 0., 100.).find_with_endpoints(p(50., on.y + 5., 0.), TOL).unwrap();
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), 5., TOL);
}

#[test]
fn partial_range() {
    let c = offset(circle(10.), 5.);
    let r = PointCurveSearch::with_domain(&c, 0., FRAC_PI_2).find_with_endpoints(p(20., 20., 0.), TOL).unwrap();
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (-TOL..=FRAC_PI_2 + TOL).contains(&e.parameter)));
}

#[test]
fn stored_distances_match_points() {
    let c = offset(circle(10.), 5.);
    for q in [p(25., 0., 0.), p(18., 12., 3.)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        for e in &r {
            assert_near(e.square_distance, q.distance(e.point).powi(2), TOL);
        }
    }
}

#[test]
fn offset_direction_sign_is_kept() {
    let basis = Curve::Circle(Circle::new(Frame::new(p(0., 0., 2.), Direction::Z, Direction::X), 3.));
    let c = Curve::Offset(OffsetCurve::new(basis, 1., -Direction::Z));
    let q = p(8., 0., 7.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find(q, TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert_near(r.min().unwrap().parameter, 0., TOL);
    assert_near(r.min_square_distance(), q.distance(c.point(0.)).powi(2), TOL);
}
