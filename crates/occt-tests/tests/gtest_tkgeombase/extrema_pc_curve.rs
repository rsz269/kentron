//! ExtremaPC_Curve_Test.cxx (closely related OCCT cases are grouped into one test).
//!
//! OCCT's transformed-curve adaptor is replaced by transforming the curve itself; the
//! `TransformedCurve_CompoundScale_*PreservesParameter` cases are not ported because a
//! scaled line or parabola is reparametrized (OCCT's adaptor keeps the old parameter).
//! `IsInitialized` checks have no counterpart (a search always has a curve).

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4, PI};

use kentron::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola,
    TrimmedCurve,
};
use kentron::geometry::extrema::{
    PointCurveExtrema, PointCurveSearch, PointExtremum, SearchMode,
};
use kentron::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};

use crate::{assert_near, assert_point_near, p};

const TOL: f64 = 1.0e-6;

fn x_line() -> Curve {
    Curve::Line(Line::new(p(0., 0., 0.), Direction::X))
}

fn trimmed(c: Curve, t1: f64, t2: f64) -> Curve {
    Curve::Trimmed(TrimmedCurve::new(c, t1, t2))
}

fn circle10() -> Curve {
    Curve::Circle(Circle::new(Frame::WORLD, 10.))
}

fn test_bspline() -> Curve {
    Curve::BSpline(BSplineCurve::new(
        vec![p(0., 0., 0.), p(1., 2., 0.), p(2., 2., 0.), p(3., 0., 0.)],
        vec![0., 1.],
        vec![4, 4],
        3,
    ))
}

fn test_bezier() -> Curve {
    Curve::Bezier(BezierCurve::new(vec![p(0., 0., 0.), p(1., 3., 0.), p(3., 3., 0.), p(4., 0., 0.)]))
}

fn find(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find(q, TOL).unwrap()
}

fn find_mode(c: &Curve, q: Point, mode: SearchMode) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with(q, TOL, mode, false).unwrap()
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

/// Every extremum's stored squared distance matches its point.
#[track_caller]
fn assert_distances_consistent(q: Point, r: &PointCurveExtrema) {
    for e in r {
        assert_near(e.square_distance, q.distance(e.point).powi(2), TOL);
    }
}

#[test]
fn line_segment_extrema() {
    let c = x_line();
    let s = PointCurveSearch::with_domain(&c, 0., 100.);
    for (q, d2) in [(p(5., 0., 0.), 0.), (p(5., 3., 4.), 25.)] {
        let r = s.find_with_endpoints(q, TOL).unwrap();
        assert_eq!(r.len(), 3, "1 minimum + 2 maxima at the ends");
        let min = r.min().unwrap();
        assert_near(min.parameter, 5., TOL);
        assert_near(min.square_distance, d2, TOL);
        assert!(min.is_minimum);
    }
    let r = s.find_with_endpoints(p(-50., 3., 4.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert_near(r.min().unwrap().parameter, 0., TOL);
}

#[test]
fn circle_extrema() {
    let c = circle10();
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find_with_endpoints(p(20., 0., 0.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    let (min, max) = (r.min().unwrap(), r.max().unwrap());
    assert_near(min.distance(), 10., TOL);
    assert!(min.is_minimum);
    assert_near(max.distance(), 30., TOL);
    assert!(!max.is_minimum);

    let r = with_endpoints(&c, p(3., 0., 0.));
    assert_eq!(r.len(), 2);
    assert_near(r.min_square_distance().sqrt(), 7., TOL);

    let r = with_endpoints(&c, p(0., 0., 0.));
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap().sqrt(), 10., TOL);
}

#[test]
fn every_curve_kind() {
    let c = trimmed(x_line(), -100., 100.);
    let r = with_endpoints(&c, p(5., 3., 4.));
    assert_eq!(r.len(), 3);
    assert_near(r.min().unwrap().parameter, 5., TOL);
    assert_near(r.min_square_distance(), 25., TOL);

    let r = find(&circle10(), p(20., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), 10., TOL);

    find(&test_bspline(), p(1.5, 3., 0.));

    let r = find(&Curve::Ellipse(Ellipse::new(Frame::WORLD, 20., 10.)), p(30., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), 10., TOL);

    let r = find(&trimmed(Curve::Parabola(Parabola::new(Frame::WORLD, 2.)), -10., 10.), p(2., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), 2., TOL);

    let q = p(10., 5., 0.);
    let r = find(&trimmed(Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 10., 5.)), -2., 2.), q);
    assert!(!r.is_empty());
    assert_distances_consistent(q, &r);

    let q = p(2., 5., 0.);
    let r = find(&test_bezier(), q);
    assert!(!r.is_empty());
    assert_distances_consistent(q, &r);

    let r = find(&Curve::Offset(OffsetCurve::new(circle10(), 5., Direction::Z)), p(25., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance().sqrt(), 10., TOL);

    let r = find(&trimmed(circle10(), 0., FRAC_PI_2), p(15., 15., 0.));
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (-TOL..=FRAC_PI_2 + TOL).contains(&e.parameter)));

    let r = find(&circle10(), p(10., 0., 0.));
    assert_near(r.min_square_distance().sqrt(), 0., TOL);

    let r = find(&trimmed(x_line(), 0., 100.), p(50., 5., 0.));
    assert!(!r.is_empty());
    assert_near(r[0].parameter, 50., TOL);
}

fn extremum(square_distance: f64, is_minimum: bool) -> PointExtremum {
    PointExtremum { parameter: 0., point: Point::ORIGIN, square_distance, is_minimum, is_maximum: !is_minimum }
}

#[test]
fn result_accessors() {
    let r = PointCurveExtrema {
        extrema: vec![extremum(100., false), extremum(50., true), extremum(75., false)],
        infinite_square_distance: None,
    };
    assert_near(r.min_square_distance(), 50., TOL);
    assert_eq!(r.min_index(), Some(1));
    assert_near(r.max_square_distance(), 100., TOL);
    assert_eq!(r.max_index(), Some(0));

    let empty = PointCurveExtrema::<Point> { extrema: vec![], infinite_square_distance: None };
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.min_index(), None);
    assert!(empty.min_square_distance().is_infinite());
}

#[test]
fn search_modes() {
    let c = circle10();
    let q = p(20., 0., 0.);
    let s = PointCurveSearch::new(&c);
    assert_eq!(s.find_with(q, TOL, SearchMode::MinMax, true).unwrap().len(), 2);
    for endpoints in [true, false] {
        let r = s.find_with(q, TOL, SearchMode::Min, endpoints).unwrap();
        assert_eq!(r.len(), 1);
        assert!(r[0].is_minimum);
        assert_near(r[0].distance(), 10., TOL);
        let r = s.find_with(q, TOL, SearchMode::Max, endpoints).unwrap();
        assert_eq!(r.len(), 1);
        assert!(!r[0].is_minimum);
        assert_near(r[0].distance(), 30., TOL);
    }
    let r = find_mode(&test_bspline(), p(1.5, 3., 0.), SearchMode::Min);
    assert!(r.iter().all(|e| e.is_minimum));
}

fn translation(x: f64, y: f64, z: f64) -> Transform {
    Transform::translation(Vector::new(x, y, z))
}

fn z_rotation(angle: f64) -> Transform {
    Transform::rotation(Axis::Z, angle)
}

/// The search on a moved curve, for the moved point, matches the original search.
#[track_caller]
fn assert_moves_with(c: &Curve, q: Point, t: &Transform, with_parameter: bool) {
    let r = find(c, q);
    let rt = find(&c.transformed(t), q.transformed(t));
    assert_eq!(r.len(), rt.len());
    assert_near(r.min_square_distance(), rt.min_square_distance(), TOL);
    let (min, min_t) = (r.min().unwrap(), rt.min().unwrap());
    if with_parameter {
        assert_near(min.parameter, min_t.parameter, TOL);
    }
    assert_point_near(min_t.point, min.point.transformed(t), TOL);
}

#[test]
fn identity_transform() {
    assert_moves_with(&trimmed(x_line(), -100., 100.), p(5., 3., 4.), &Transform::IDENTITY, true);
    assert_moves_with(&test_bspline(), p(1.5, 3., 0.), &Transform::IDENTITY, true);
}

#[test]
fn translated_curves() {
    let c = trimmed(x_line(), -100., 100.).transformed(&translation(0., 10., 0.));
    let r = find(&c, p(5., 13., 4.));
    let min = r.min().expect("no extrema");
    assert_near(min.parameter, 5., TOL);
    assert_near(min.square_distance, 25., TOL);
    assert_near(min.point.y, 10., TOL);

    let c = circle10().transformed(&translation(100., 0., 0.));
    assert_near(find(&c, p(120., 0., 0.)).min_square_distance().sqrt(), 10., TOL);

    assert_moves_with(&test_bspline(), p(1.5, 3., 0.), &translation(10., 20., 30.), true);
    assert_moves_with(&test_bezier(), p(2., 5., 0.), &translation(10., 20., 30.), true);

    let c = Curve::Offset(OffsetCurve::new(circle10(), 5., Direction::Z)).transformed(&translation(50., 0., 0.));
    assert_near(find(&c, p(75., 0., 0.)).min_square_distance().sqrt(), 10., TOL);
}

#[test]
fn rotated_curves() {
    let c = trimmed(x_line(), 0., 100.).transformed(&z_rotation(FRAC_PI_2));
    let r = find(&c, p(3., 5., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.distance(), 3., TOL);
    assert_point_near(min.point, p(0., 5., 0.), TOL);

    assert_moves_with(&test_bspline(), p(1.5, 3., 0.), &z_rotation(FRAC_PI_2), true);
    assert_moves_with(&test_bspline(), p(1.5, 3., 0.), &(z_rotation(FRAC_PI_4) * translation(5., 5., 0.)), true);
    assert_moves_with(&test_bezier(), p(2., 5., 0.), &(z_rotation(FRAC_PI_3) * translation(10., 0., 0.)), false);
}

#[test]
fn scaled_free_form_curves() {
    let scale = |f| Transform::scaling(Point::ORIGIN, f);
    let c = test_bspline();
    let r = find(&c, p(1.5, 3., 0.));
    let rt = find(&c.transformed(&scale(2.)), p(3., 6., 0.));
    assert_near(rt.min_square_distance(), 4. * r.min_square_distance(), TOL);
    assert_near(r.min().unwrap().parameter, rt.min().unwrap().parameter, TOL);

    let r = with_endpoints(&c, p(1.5, 3., 0.));
    let rt = with_endpoints(&c.transformed(&scale(2.)), p(3., 6., 0.));
    assert_near(rt.min_square_distance(), 4. * r.min_square_distance(), TOL);

    let c = test_bezier();
    let r = find(&c, p(2., 5., 0.));
    let rt = find(&c.transformed(&scale(3.)), p(6., 15., 0.));
    assert_near(rt.min_square_distance(), 9. * r.min_square_distance(), 1.0e-4);
}

#[test]
fn moved_and_restricted_domains() {
    let c = trimmed(circle10(), 0., FRAC_PI_2).transformed(&translation(0., 0., 50.));
    for e in &find(&c, p(15., 15., 50.)) {
        assert!((-TOL..=FRAC_PI_2 + TOL).contains(&e.parameter));
        assert_near(e.point.z, 50., TOL);
    }

    let c = x_line().transformed(&translation(0., 5., 0.));
    let r = PointCurveSearch::with_domain(&c, 0., 10.).find_with_endpoints(p(5., 8., 0.), TOL).unwrap();
    let min = r.min().expect("no extrema");
    assert_near(min.parameter, 5., TOL);
    assert_near(min.square_distance, 9., TOL);
    assert_near(min.point.y, 5., TOL);

    let c = test_bspline();
    let r = with_endpoints(&c, p(1.5, 3., 0.));
    let rt = with_endpoints(&c.transformed(&translation(0., 0., 10.)), p(1.5, 3., 10.));
    assert_eq!(r.len(), rt.len());
    assert_near(r.min_square_distance(), rt.min_square_distance(), TOL);
}

#[test]
fn moved_curve_search_modes() {
    let c = circle10().transformed(&translation(0., 0., 5.));
    let q = p(20., 0., 5.);
    let r = find_mode(&c, q, SearchMode::Min);
    assert_eq!(r.len(), 1);
    assert!(r[0].is_minimum);
    assert_near(r[0].distance(), 10., TOL);
    let r = find_mode(&c, q, SearchMode::Max);
    assert_eq!(r.len(), 1);
    assert!(!r[0].is_minimum);
    assert_near(r[0].distance(), 30., TOL);
}
