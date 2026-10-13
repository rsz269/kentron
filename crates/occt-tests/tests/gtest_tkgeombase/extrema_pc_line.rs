//! ExtremaPC_Line_Test.cxx (closely related OCCT cases are grouped into one test).

use kentron::geometry::curve::{Curve, Line};
use kentron::geometry::curve2::{Curve2, Line2};
use kentron::geometry::extrema::{ExtremaError, PointCurveSearch};
use kentron::math::{Direction, Direction2, Point};

use crate::{assert_min, assert_near, assert_point_near, p, p2};

const TOL: f64 = 1.0e-10;

fn line(origin: Point, direction: Direction) -> Curve {
    Curve::Line(Line::new(origin, direction))
}

fn x_line() -> Curve {
    line(p(0., 0., 0.), Direction::X)
}

/// `find_with_endpoints` on `curve` over `[t1, t2]`: expects `count` extrema and the given
/// closest parameter and squared distance.
#[track_caller]
fn check(curve: &Curve, (t1, t2): (f64, f64), q: Point, count: usize, t: Option<f64>, d2: Option<f64>) {
    let r = PointCurveSearch::with_domain(curve, t1, t2).find_with_endpoints(q, TOL).unwrap();
    assert_eq!(r.len(), count, "{r:?}");
    assert_min(&r, t, d2, TOL);
}

const WIDE: (f64, f64) = (-100., 100.);

#[test]
fn singleton_domain_is_minimum_and_maximum() {
    let c = x_line();
    for q in [p(5., 3., 0.), p(2., 3., 0.)] {
        let r = PointCurveSearch::with_domain(&c, 2., 2.).find_with_endpoints(q, TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert!(r[0].is_minimum && r[0].is_maximum);
        assert_near(r[0].parameter, 2., TOL);
    }
}

#[test]
fn effectively_unbounded_domain_does_not_add_endpoints() {
    check(&x_line(), (-1.0e100, 1.0e100), p(2., 3., 0.), 1, Some(2.), None);
}

#[test]
fn half_bounded_domains_keep_the_finite_end() {
    let c = x_line();
    for (domain, q) in [((0., 1.0e100), p(-2., 3., 0.)), ((-1.0e100, 0.), p(2., 3., 0.))] {
        let r = PointCurveSearch::with_domain(&c, domain.0, domain.1).find_with_endpoints(q, TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
        assert!(r[0].is_minimum);
    }
}

#[test]
fn reversed_domain_is_invalid_input() {
    let c = x_line();
    assert_eq!(PointCurveSearch::with_domain(&c, 1., -1.).find(p(0., 0., 0.), TOL), Err(ExtremaError::InvalidInput));
}

#[test]
fn point_on_line() {
    let r = PointCurveSearch::with_domain(&x_line(), -100., 100.).find_with_endpoints(p(0., 0., 0.), TOL).unwrap();
    assert_eq!(r.len(), 3, "1 minimum + 2 maxima at the ends");
    assert_min(&r, Some(0.), Some(0.), TOL);
    assert!(r.min().unwrap().is_minimum);
    check(&x_line(), WIDE, p(5., 0., 0.), 3, Some(5.), Some(0.));
    check(&x_line(), WIDE, p(-7.5, 0., 0.), 3, Some(-7.5), Some(0.));
}

#[test]
fn point_off_line() {
    check(&x_line(), WIDE, p(5., 3., 0.), 3, Some(5.), Some(9.));
    check(&x_line(), WIDE, p(5., 0., 4.), 3, Some(5.), Some(16.));
    check(&x_line(), WIDE, p(5., 3., 4.), 3, Some(5.), Some(25.));
}

#[test]
fn projection_outside_bounds_clamps_to_an_end() {
    check(&x_line(), (0., 100.), p(-50., 3., 4.), 2, Some(0.), None);
    check(&x_line(), (0., 100.), p(150., 3., 4.), 2, Some(100.), None);
}

#[test]
fn projection_at_a_bound() {
    check(&x_line(), (10., 100.), p(10., 5., 0.), 2, Some(10.), None);
    check(&x_line(), WIDE, p(100., 5., 0.), 2, Some(100.), None);
}

#[test]
fn narrow_bounds() {
    check(&x_line(), WIDE, p(50., 1., 0.), 3, Some(50.), None);
}

#[test]
fn other_axes() {
    check(&line(p(0., 0., 0.), Direction::Y), WIDE, p(3., 5., 4.), 3, Some(5.), Some(25.));
    check(&line(p(0., 0., 0.), Direction::Z), WIDE, p(3., 4., 10.), 3, Some(10.), Some(25.));
}

#[test]
fn diagonal_lines() {
    let c = line(p(0., 0., 0.), Direction::new(1., 1., 0.));
    check(&c, WIDE, p(5., 5., 3.), 3, Some(5. * 2f64.sqrt()), Some(9.));
    let c = line(p(0., 0., 0.), Direction::new(1., 1., 1.));
    check(&c, WIDE, p(3., 3., 3.), 3, None, Some(0.));
}

#[test]
fn offset_origin() {
    let c = line(p(10., 20., 30.), Direction::X);
    check(&c, WIDE, p(15., 20., 30.), 3, Some(5.), Some(0.));
    check(&c, WIDE, p(15., 23., 34.), 3, Some(5.), Some(25.));
}

#[test]
fn negative_parameters() {
    check(&x_line(), WIDE, p(-25., 0., 0.), 3, Some(-25.), None);
    check(&x_line(), (-50., -10.), p(-5., 1., 0.), 2, Some(-10.), None);
}

#[test]
fn large_coordinates() {
    let c = line(p(1e6, 1e6, 1e6), Direction::X);
    let r = PointCurveSearch::with_domain(&c, -100., 200.)
        .find_with_endpoints(p(1e6 + 100., 1e6 + 3., 1e6 + 4.), TOL)
        .unwrap();
    assert_eq!(r.len(), 3);
    assert_min(&r, Some(100.), Some(25.), 1e-6);
}

#[test]
fn small_coordinates() {
    let c = line(p(1e-6, 1e-6, 1e-6), Direction::X);
    check(&c, WIDE, p(1e-6 + 1e-3, 1e-6 + 3e-4, 1e-6 + 4e-4), 3, Some(1e-3), Some(25e-8));
}

#[test]
fn projected_point_and_distance_are_consistent() {
    let c = x_line();
    let r = PointCurveSearch::with_domain(&c, -100., 100.).find_with_endpoints(p(7., 5., 0.), TOL).unwrap();
    assert_eq!(r.len(), 3);
    assert_point_near(r.min().unwrap().point, p(7., 0., 0.), TOL);

    let q = p(10., 6., 8.);
    let r = PointCurveSearch::with_domain(&c, -100., 100.).find_with_endpoints(q, TOL).unwrap();
    assert_eq!(r.len(), 3);
    let min = r.min().unwrap();
    assert_near(min.square_distance, q.distance(min.point).powi(2), TOL);
    assert_near(r.min_square_distance(), 100., TOL);
}

#[test]
fn ends_are_maxima() {
    let c = x_line();
    let r = PointCurveSearch::with_domain(&c, 0., 10.).find_with_endpoints(p(5., 1., 0.), TOL).unwrap();
    assert_eq!(r.len(), 3);
    assert_min(&r, Some(5.), Some(1.), TOL);
    let t = r.max().unwrap().parameter;
    assert!(t.abs() < TOL || (t - 10.).abs() < TOL);
}

#[test]
fn planar_case_matches_2d() {
    let c3 = x_line();
    let c2 = Curve2::Line(Line2::new(p2(0., 0.), Direction2::X));
    let r3 = PointCurveSearch::new(&c3).find(p(8., 3., 6.), TOL).unwrap();
    let r2 = PointCurveSearch::new(&c2).find(p2(8., 3.), TOL).unwrap();
    assert_eq!(r3.len(), 1);
    assert_eq!(r2.len(), 1);
    assert_near(r3[0].parameter, r2[0].parameter, TOL);
    assert_near(r3[0].square_distance, r2[0].square_distance + 36., TOL);
}
