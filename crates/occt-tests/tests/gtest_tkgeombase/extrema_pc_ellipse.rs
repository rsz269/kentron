//! ExtremaPC_Ellipse_Test.cxx (closely related OCCT cases are grouped into one test).

use std::f64::consts::PI;

use kentron::geometry::curve::{Curve, Ellipse};
use kentron::geometry::curve2::{Curve2, Ellipse2};
use kentron::geometry::extrema::{ExtremaError, PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Frame2, Point};

use crate::{assert_near, assert_planar_matches_2d, p};

const TOL: f64 = 1.0e-8;

fn ellipse(frame: Frame, major: f64, minor: f64) -> Curve {
    Curve::Ellipse(Ellipse::new(frame, major, minor))
}

/// The 20 x 10 ellipse in the XY plane.
fn e20() -> Curve {
    ellipse(Frame::WORLD, 20., 10.)
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

/// Every extremum's point is the curve's point at its parameter.
#[track_caller]
fn assert_points_on_curve(c: &Curve, r: &PointCurveExtrema, tolerance: f64) {
    for e in r {
        assert!(e.point.distance(c.point(e.parameter)) <= tolerance, "{e:?}");
    }
}

#[test]
fn shifted_full_period_gives_unique_representatives() {
    let c = e20();
    let u_min = -10.;
    let r = PointCurveSearch::with_domain(&c, u_min, u_min + 2. * PI).find_with_endpoints(p(30., 0., 0.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert!(r.iter().all(|e| (u_min..=u_min + 2. * PI).contains(&e.parameter)));
}

#[test]
fn zero_tolerance_is_invalid() {
    let c = e20();
    let s = PointCurveSearch::with_domain(&c, 0.5, 0.5);
    assert_eq!(s.find_with_endpoints(p(30., 0., 0.), 0.), Err(ExtremaError::InvalidInput));
}

#[test]
fn point_on_major_axis() {
    let c = e20();
    for q in [p(30., 0., 0.), p(-30., 0., 0.)] {
        let r = with_endpoints(&c, q);
        assert!(r.len() >= 2);
        assert_near(r.min_square_distance().sqrt(), 10., TOL);
    }
    // Inside: the closest points are at cos(u) = 2/3, not at the vertex.
    let r = with_endpoints(&c, p(10., 0., 0.));
    assert!(r.len() >= 2);
    assert_near(r.min_square_distance(), 600. / 9., 0.01);
}

#[test]
fn point_on_minor_axis() {
    let c = e20();
    for (q, d) in [(p(0., 20., 0.), 10.), (p(0., 5., 0.), 5.)] {
        let r = with_endpoints(&c, q);
        assert!(r.len() >= 2);
        assert_near(r.min_square_distance().sqrt(), d, TOL);
    }
}

#[test]
fn point_in_quadrants() {
    let c = e20();
    for q in [p(15., 15., 0.), p(-15., 8., 0.), p(-12., -7., 0.), p(18., -5., 0.)] {
        let r = with_endpoints(&c, q);
        assert!(r.len() >= 2);
        let min = r.min().unwrap();
        assert!(min.square_distance > 0.);
        assert!(min.point.distance(c.point(min.parameter)) <= TOL);
        assert!(min.distance() < q.distance(p(0., 0., 0.)));
    }
}

#[test]
fn point_at_center() {
    let c = e20();
    let r = with_endpoints(&c, p(0., 0., 0.));
    assert!(r.len() >= 4, "minima at the minor vertices, maxima at the major ones");
    assert_near(r.min_square_distance().sqrt(), 10., TOL);
    assert_near(r.max_square_distance().sqrt(), 20., TOL);
    assert_points_on_curve(&c, &r, TOL);
}

#[test]
fn point_out_of_plane() {
    let c = e20();
    assert_near(with_endpoints(&c, p(30., 0., 5.)).min_square_distance(), 125., TOL);
    assert_near(with_endpoints(&c, p(0., 20., -8.)).min_square_distance(), 164., TOL);
}

/// OCCT's PartialArc_* cases (which use the full ellipse).
#[test]
fn closest_vertex() {
    let c = e20();
    for q in [p(30., 0., 0.), p(-30., 0., 0.)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_near(r.min_square_distance().sqrt(), 10., TOL);
    }
}

#[test]
fn other_placements() {
    let c = ellipse(Frame::from_z(p(100., 200., 50.), Direction::Z), 20., 10.);
    assert_near(with_endpoints(&c, p(130., 200., 50.)).min_square_distance().sqrt(), 10., TOL);

    let c = ellipse(Frame::from_z(p(0., 0., 0.), Direction::X), 20., 10.);
    let r = with_endpoints(&c, p(5., 30., 0.));
    assert!(!r.is_empty());
    assert!(r.min_square_distance() > 0.);
    assert_points_on_curve(&c, &r, 1e-7);
}

#[test]
fn nearly_circular_and_very_flat() {
    let r = with_endpoints(&ellipse(Frame::WORLD, 10., 9.9), p(20., 0., 0.));
    assert!(r.len() >= 2);
    assert_near(r.min_square_distance().sqrt(), 10., TOL);
    let r = with_endpoints(&ellipse(Frame::WORLD, 20., 2.), p(0., 10., 0.));
    assert!(r.len() >= 2);
    assert_near(r.min_square_distance().sqrt(), 8., TOL);
}

#[test]
fn point_on_ellipse() {
    let c = e20();
    for q in [p(20., 0., 0.), p(0., 10., 0.), c.point(PI / 4.)] {
        assert_near(with_endpoints(&c, q).min_square_distance(), 0., TOL);
    }
}

#[test]
fn extremum_counts() {
    let c = e20();
    let r = with_endpoints(&c, p(5., 3., 0.));
    assert!(r.len() >= 4, "a point inside has 2 minima and 2 maxima");
    assert!(r.iter().all(|e| e.square_distance > 0.));
    assert_eq!(with_endpoints(&c, p(30., 0., 0.)).len(), 2);
}

#[test]
fn narrow_ranges() {
    let c = e20();
    let r = PointCurveSearch::with_domain(&c, 0., 0.05).find_with_endpoints(p(22., 2., 0.), TOL).unwrap();
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (-TOL..=0.05 + TOL).contains(&e.parameter)));

    let q = p(25., 5., 0.);
    let r = PointCurveSearch::with_domain(&c, PI / 6. - 0.1, PI / 6. + 0.1).find_with_endpoints(q, TOL).unwrap();
    for e in &r {
        let on = c.point(e.parameter);
        assert_near(on.distance(e.point), 0., TOL);
        assert_near(e.square_distance, q.distance(on).powi(2), TOL);
    }

    let r = PointCurveSearch::with_domain(&c, 0., 0.001).find_with_endpoints(p(21., 0.5, 0.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert_eq!(r[0].parameter, 0.);
    assert!(r[0].is_maximum);
    assert_eq!(r[1].parameter, 0.001);
    assert!(r[1].is_minimum);
}

#[test]
fn planar_case_matches_2d() {
    assert_planar_matches_2d(&ellipse(Frame::WORLD, 5., 2.), &Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 5., 2.)), TOL);
}
