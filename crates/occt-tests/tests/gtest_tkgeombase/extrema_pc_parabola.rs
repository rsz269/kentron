//! ExtremaPC_Parabola_Test.cxx (closely related OCCT cases are grouped into one test).

use kentron::geometry::curve::{Curve, Parabola};
use kentron::geometry::curve2::{Curve2, Parabola2};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Frame2, Point};

use crate::{assert_near, assert_planar_matches_2d, p};

const TOL: f64 = 1.0e-8;

fn parabola(frame: Frame, focal: f64) -> Curve {
    Curve::Parabola(Parabola::new(frame, focal))
}

/// Focal length 5 in the XY plane, opening along X.
fn p5() -> Curve {
    parabola(Frame::WORLD, 5.)
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

#[track_caller]
fn assert_min_on_curve(c: &Curve, r: &PointCurveExtrema) {
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
}

#[test]
fn point_on_axis() {
    let c = p5();
    for (q, d2) in [(p(0., 0., 0.), 0.), (p(10., 0., 0.), 100.), (p(-10., 0., 0.), 100.)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_near(r.min_square_distance(), d2, TOL);
    }
}

#[test]
fn zero_focal_degenerates_to_line() {
    let c = parabola(Frame::new(p(1., 2., 3.), Direction::Z, Direction::X), 0.);
    let r = PointCurveSearch::with_domain(&c, -10., 10.).find(p(5., 7., 3.), TOL).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, 4., TOL);
    assert_near(r[0].point.distance(p(5., 2., 3.)), 0., TOL);
}

#[test]
fn point_off_axis() {
    let c = p5();
    let q = p(5., 10., 0.);
    let r = with_endpoints(&c, q);
    let min = r.min().expect("no extrema");
    assert_near(min.square_distance, q.distance(c.point(min.parameter)).powi(2), TOL);

    let q = p(5., -10., 0.);
    let r = with_endpoints(&c, q);
    assert_min_on_curve(&c, &r);
    let min = r.min().unwrap();
    assert_near(min.square_distance, q.distance(c.point(min.parameter)).powi(2), TOL);
    assert!(min.parameter < 0.);

    let r = with_endpoints(&c, p(0., 50., 0.));
    assert_min_on_curve(&c, &r);
    assert!(r.min_square_distance().sqrt() < 50.);
    assert!(r.min().unwrap().parameter > 0.);
}

#[test]
fn point_out_of_plane() {
    let c = p5();
    let r = with_endpoints(&c, p(5., 10., 7.));
    assert!(!r.is_empty());
    assert!(r.min_square_distance() >= 49.);
    assert_near(with_endpoints(&c, p(0., 0., -5.)).min_square_distance(), 25., TOL);
}

#[test]
fn bounded_domains() {
    let c = p5();
    let search = |t1, t2, q| PointCurveSearch::with_domain(&c, t1, t2).find_with_endpoints(q, TOL).unwrap();
    assert!(search(0., 10., p(0., -10., 0.)).min().unwrap().parameter >= 0.);
    assert!(search(-10., 0., p(0., 10., 0.)).min().unwrap().parameter <= 0.);
    let r = search(4., 6., p(5., 5., 0.));
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (4. - TOL..=6. + TOL).contains(&e.parameter)));
}

#[test]
fn other_focal_lengths() {
    let c = parabola(Frame::WORLD, 0.5);
    let q = p(1., 2., 0.);
    let r = with_endpoints(&c, q);
    assert_min_on_curve(&c, &r);
    let min = r.min().unwrap();
    assert_near(min.square_distance, q.distance(c.point(min.parameter)).powi(2), TOL);

    let c = parabola(Frame::WORLD, 100.);
    assert_min_on_curve(&c, &with_endpoints(&c, p(50., 100., 0.)));
}

#[test]
fn point_on_parabola() {
    let c = p5();
    for t in [0., 10., -15.] {
        assert_near(with_endpoints(&c, c.point(t)).min_square_distance(), 0., TOL);
    }
}

#[test]
fn offset_placement() {
    let c = parabola(Frame::new(p(100., 200., 50.), Direction::Z, Direction::X), 5.);
    assert_near(with_endpoints(&c, p(100., 200., 50.)).min_square_distance(), 0., TOL);
    let q = p(110., 210., 50.);
    let r = with_endpoints(&c, q);
    assert_min_on_curve(&c, &r);
    let min = r.min().unwrap();
    assert_near(min.square_distance, q.distance(c.point(min.parameter)).powi(2), TOL);
}

#[test]
fn other_planes() {
    for (z, x, q) in [(Direction::X, Direction::Y, p(3., 10., 5.)), (Direction::Y, Direction::X, p(10., 3., 5.))] {
        let c = parabola(Frame::new(p(0., 0., 0.), z, x), 5.);
        let r = with_endpoints(&c, q);
        assert_min_on_curve(&c, &r);
        assert!(r.min_square_distance() >= 9., "the out-of-plane offset is 3");
    }
}

#[test]
fn projected_point_and_distance_are_consistent() {
    let c = p5();
    assert_min_on_curve(&c, &with_endpoints(&c, p(10., 15., 0.)));
    let q = p(8., 12., 3.);
    let r = with_endpoints(&c, q);
    let min = r.min().expect("no extrema");
    assert_near(min.square_distance, q.distance(min.point).powi(2), TOL);
}

#[test]
fn planar_case_matches_2d() {
    assert_planar_matches_2d(
        &parabola(Frame::WORLD, 2.),
        &Curve2::Parabola(Parabola2::new(Frame2::WORLD, 2.)),
        TOL,
    );
}
