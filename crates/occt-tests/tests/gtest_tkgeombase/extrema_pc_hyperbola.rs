//! ExtremaPC_Hyperbola_Test.cxx (closely related OCCT cases are grouped into one test).

use kentron::geometry::curve::{Curve, Hyperbola};
use kentron::geometry::curve2::{Curve2, Hyperbola2};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Frame2, Point};

use crate::{assert_near, assert_planar_matches_2d, p};

const TOL: f64 = 1.0e-8;

fn hyperbola(frame: Frame, major: f64, minor: f64) -> Curve {
    Curve::Hyperbola(Hyperbola::new(frame, major, minor))
}

/// The 20 x 10 hyperbola in the XY plane.
fn h20() -> Curve {
    hyperbola(Frame::WORLD, 20., 10.)
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

/// The closest extremum lies on the curve at its parameter, at its reported distance.
#[track_caller]
fn assert_min_consistent(c: &Curve, q: Point, r: &PointCurveExtrema) {
    let min = r.min().expect("no extrema");
    let on = c.point(min.parameter);
    assert_near(min.point.distance(on), 0., TOL);
    assert_near(min.square_distance, q.distance(on).powi(2), TOL);
}

#[test]
fn point_on_major_axis() {
    let c = h20();
    // Beyond the vertex the closest points are at cosh(u) = 2, not at the vertex.
    let r = with_endpoints(&c, p(50., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance(), 400., 1.);
    assert_near(with_endpoints(&c, p(20., 0., 0.)).min_square_distance(), 0., TOL);
    assert_near(with_endpoints(&c, p(0., 0., 0.)).min_square_distance().sqrt(), 20., TOL);
}

#[test]
fn point_off_axis() {
    let c = h20();
    for q in [p(30., 15., 0.), p(30., -15., 0.), p(25., 5., 0.)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_min_consistent(&c, q, &r);
    }
    assert!(with_endpoints(&c, p(30., -15., 0.)).min().unwrap().parameter < 0.);
    assert!(with_endpoints(&c, p(25., 5., 0.)).min_square_distance().sqrt() < 26.);
}

#[test]
fn point_out_of_plane() {
    let c = h20();
    let r = with_endpoints(&c, p(30., 0., 5.));
    assert!(!r.is_empty());
    assert!(r.min_square_distance() >= 25.);
    assert_near(with_endpoints(&c, p(20., 0., -8.)).min_square_distance(), 64., TOL);
}

#[test]
fn bounded_domains() {
    let c = h20();
    let search = |t1, t2, q| PointCurveSearch::with_domain(&c, t1, t2).find_with_endpoints(q, TOL).unwrap();
    let r = search(0., 10., p(30., -20., 0.));
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| e.parameter >= -TOL));
    let r = search(-10., 0., p(30., 20., 0.));
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| e.parameter <= TOL));
    let r = search(0.5, 1.5, p(30., 10., 0.));
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (0.5 - TOL..=1.5 + TOL).contains(&e.parameter)));
}

#[test]
fn other_semi_axes() {
    for (c, q) in [
        (hyperbola(Frame::WORLD, 2., 1.), p(5., 2., 0.)),
        (hyperbola(Frame::WORLD, 10., 10.), p(20., 15., 0.)),
    ] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_min_consistent(&c, q, &r);
    }
    let c = hyperbola(Frame::WORLD, 200., 100.);
    let r = with_endpoints(&c, p(300., 100., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
}

#[test]
fn point_on_hyperbola() {
    let c = h20();
    for t in [0., 1.5, -2.] {
        assert_near(with_endpoints(&c, c.point(t)).min_square_distance(), 0., TOL);
    }
}

#[test]
fn offset_placement() {
    let c = hyperbola(Frame::from_z(p(100., 200., 50.), Direction::Z), 20., 10.);
    assert_near(with_endpoints(&c, p(120., 200., 50.)).min_square_distance(), 0., TOL);
    let q = p(150., 210., 50.);
    let r = with_endpoints(&c, q);
    assert!(!r.is_empty());
    assert_min_consistent(&c, q, &r);
}

#[test]
fn other_planes() {
    for (normal, q) in [(Direction::X, p(5., 30., 10.)), (Direction::Y, p(30., 5., 10.))] {
        let c = hyperbola(Frame::from_z(p(0., 0., 0.), normal), 20., 10.);
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        let min = r.min().unwrap();
        assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
        assert!(r.min_square_distance() >= 25., "the out-of-plane offset is 5");
    }
}

#[test]
fn projected_point_and_distance_are_consistent() {
    let c = h20();
    let r = with_endpoints(&c, p(35., 20., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);

    let q = p(40., 25., 3.);
    let r = with_endpoints(&c, q);
    let min = r.min().expect("no extrema");
    assert_near(min.square_distance, q.distance(min.point).powi(2), TOL);
}

#[test]
fn planar_case_matches_2d() {
    let c2 = Curve2::Hyperbola(Hyperbola2::new(Frame2::WORLD, 3., 1.5));
    assert_planar_matches_2d(&hyperbola(Frame::WORLD, 3., 1.5), &c2, TOL);
}
