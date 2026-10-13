//! ExtremaPC_Circle_Test.cxx (closely related OCCT cases are grouped into one test).

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve::{Circle, Curve};
use kentron::geometry::curve2::{Circle2, Curve2};
use kentron::geometry::extrema::{ExtremaError, PointCurveExtrema, PointCurveSearch, PointExtremum};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Frame2, Point};
use kentron::tolerance::INFINITE;

use crate::{assert_near, assert_planar_matches_2d, p};

const TOL: f64 = 1.0e-10;

fn circle(frame: Frame, radius: f64) -> Curve {
    Curve::Circle(Circle::new(frame, radius))
}

fn xy_circle(radius: f64) -> Curve {
    circle(Frame::WORLD, radius)
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

/// Expects exactly a minimum and a maximum; returns them.
#[track_caller]
fn min_max(r: &PointCurveExtrema) -> (PointExtremum, PointExtremum) {
    assert_eq!(r.len(), 2, "{r:?}");
    (*r.min().unwrap(), *r.max().unwrap())
}

#[test]
fn point_outside_on_x_axis() {
    let c = xy_circle(10.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find_with_endpoints(p(20., 0., 0.), TOL).unwrap();
    let (min, max) = min_max(&r);
    assert_near(min.distance(), 10., TOL);
    assert!(min.is_minimum);
    assert_near(min.parameter, 0., TOL);
    assert_near(max.distance(), 30., TOL);
    assert!(!max.is_minimum);
    assert_near(max.parameter, PI, TOL);
}

#[test]
fn shifted_full_period_gives_unique_representatives() {
    let c = xy_circle(10.);
    let u_min = 7.;
    let r = PointCurveSearch::with_domain(&c, u_min, u_min + 2. * PI).find_with_endpoints(p(20., 0., 0.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert!(r.iter().all(|e| (u_min..=u_min + 2. * PI).contains(&e.parameter)));
}

#[test]
fn large_tolerance_preserves_opposite_extrema() {
    let c = xy_circle(1.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find(p(2., 0., 0.), 4.).unwrap();
    assert_eq!(r.len(), 2);
    assert_near(r.min_square_distance(), 1., TOL);
    assert_near(r.max_square_distance(), 9., TOL);
}

#[test]
fn more_than_one_period_keeps_constrained_endpoint() {
    let c = xy_circle(10.);
    let end = 2. * PI + 0.1;
    let r = PointCurveSearch::with_domain(&c, 0., end).find_with_endpoints(c.point(end), TOL).unwrap();
    assert!(r.iter().any(|e| (e.parameter - end).abs() <= TOL));
}

#[test]
fn point_outside() {
    let c = xy_circle(10.);
    let (min, max) = min_max(&with_endpoints(&c, p(0., 25., 0.)));
    assert_near(min.distance(), 15., TOL);
    assert_near(min.parameter, FRAC_PI_2, TOL);
    assert_near(max.distance(), 35., TOL);

    let (min, _) = min_max(&with_endpoints(&c, p(20., 20., 0.)));
    assert_near(min.distance(), 800f64.sqrt() - 10., TOL);
}

#[test]
fn point_inside() {
    let c = xy_circle(10.);
    for (q, near, far) in [(p(3., 0., 0.), 7., 13.), (p(0., 5., 0.), 5., 15.)] {
        let (min, max) = min_max(&with_endpoints(&c, q));
        assert_near(min.distance(), near, TOL);
        assert_near(max.distance(), far, TOL);
    }
    let (min, _) = min_max(&with_endpoints(&c, p(0.1, 0.1, 0.)));
    assert_near(min.distance(), 10. - 0.02f64.sqrt(), TOL);
}

#[test]
fn point_at_center_is_degenerate() {
    let r = with_endpoints(&xy_circle(10.), p(0., 0., 0.));
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap().sqrt(), 10., TOL);
}

#[test]
fn point_near_center_is_not_degenerate() {
    let c = xy_circle(10.);
    let r = PointCurveSearch::new(&c).find_with_endpoints(p(1e-8, 0., 0.), 1e-12).unwrap();
    let (min, max) = min_max(&r);
    assert_near(min.distance(), 10. - 1e-8, 1e-6);
    assert_near(max.distance(), 10. + 1e-8, 1e-6);
    assert_near(min.parameter, 0., TOL);
    assert_near(max.parameter, PI, TOL);
}

#[test]
fn point_out_of_plane() {
    let c = xy_circle(10.);
    assert_near(min_max(&with_endpoints(&c, p(20., 0., 5.))).0.square_distance, 125., TOL);
    assert_near(min_max(&with_endpoints(&c, p(0., 15., -7.))).0.square_distance, 74., TOL);
}

#[test]
fn point_on_axis_is_degenerate() {
    let r = with_endpoints(&xy_circle(10.), p(0., 0., 15.));
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap(), 325., TOL);
}

/// OCCT's PartialArc_* cases (which use the full circle).
#[test]
fn closest_and_farthest_points() {
    let c = xy_circle(10.);
    for q in [p(20., 0., 0.), p(-20., 0., 0.)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_near(r.min_square_distance().sqrt(), 10., TOL);
    }
    let r = with_endpoints(&c, p(20., 0., 0.));
    let at_pi = r.iter().find(|e| (e.parameter - PI).abs() < 0.1).expect("extremum at PI");
    assert_near(at_pi.distance(), 30., TOL);
}

#[test]
fn other_placements() {
    let c = circle(Frame::from_z(p(100., 200., 50.), Direction::Z), 15.);
    assert_near(min_max(&with_endpoints(&c, p(130., 200., 50.))).0.distance(), 15., TOL);
    let c = circle(Frame::from_z(p(0., 0., 0.), Direction::X), 10.);
    assert_near(min_max(&with_endpoints(&c, p(5., 20., 0.))).0.square_distance, 125., TOL);
    let c = circle(Frame::from_z(p(0., 0., 0.), Direction::Y), 10.);
    assert_near(min_max(&with_endpoints(&c, p(20., 3., 0.))).0.square_distance, 109., TOL);
}

#[test]
fn small_and_large_radius() {
    assert_near(min_max(&with_endpoints(&xy_circle(0.001), p(0.01, 0., 0.))).0.distance(), 0.009, TOL);
    assert_near(min_max(&with_endpoints(&xy_circle(1000.), p(1500., 0., 0.))).0.distance(), 500., TOL);
}

#[test]
fn point_on_circle() {
    let c = xy_circle(10.);
    for (q, t) in [(p(10., 0., 0.), 0.), (p(-10., 0., 0.), PI)] {
        let (min, _) = min_max(&with_endpoints(&c, q));
        assert_near(min.square_distance, 0., TOL);
        assert_near(min.parameter, t, TOL);
    }
}

#[test]
fn non_finite_query_is_rejected() {
    let c = xy_circle(10.);
    let s = PointCurveSearch::new(&c);
    assert_eq!(s.find(p(INFINITE, 0., 0.), TOL), Err(ExtremaError::InvalidInput));
    assert_eq!(s.find_with_endpoints(p(INFINITE, 0., 0.), TOL), Err(ExtremaError::InvalidInput));
}

#[test]
fn planar_case_matches_2d() {
    assert_planar_matches_2d(&xy_circle(4.), &Curve2::Circle(Circle2::new(Frame2::WORLD, 4.)), TOL);
}
