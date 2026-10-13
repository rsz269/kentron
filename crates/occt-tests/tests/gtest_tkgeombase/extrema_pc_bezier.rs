//! ExtremaPC_BezierCurve_Test.cxx (closely related OCCT cases are grouped into one test).
//!
//! `SmallNonConstantCurveIsNotPromotedToInfinite` drives OCCT's internal grid evaluator;
//! it is ported through the public search.

use kentron::geometry::curve::{BezierCurve, Curve};
use kentron::geometry::curve2::{BezierCurve2, Curve2};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::Point;

use crate::{assert_near, p, p2};

const TOL: f64 = 1.0e-6;

fn bezier(poles: &[(f64, f64, f64)]) -> Curve {
    Curve::Bezier(BezierCurve::new(poles.iter().map(|&(x, y, z)| p(x, y, z)).collect()))
}

fn rational_quadratic(p0: Point, p1: Point, p2_: Point, w1: f64) -> Curve {
    Curve::Bezier(BezierCurve::rational(vec![p0, p1, p2_], vec![1., w1, 1.]))
}

/// The arch (0,0) (1,2) (3,2) (4,0).
fn arch() -> Curve {
    bezier(&[(0., 0., 0.), (1., 2., 0.), (3., 2., 0.), (4., 0., 0.)])
}

fn with_endpoints(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with_endpoints(q, TOL).unwrap()
}

/// Each extremum's point and squared distance match the curve at its parameter.
#[track_caller]
fn assert_consistent(c: &Curve, q: Point, r: &PointCurveExtrema) {
    assert!(!r.is_empty());
    for e in r {
        let on = c.point(e.parameter);
        assert_near(e.point.distance(on), 0., TOL);
        assert_near(e.square_distance, q.distance(on).powi(2), TOL);
    }
}

#[test]
fn point_on_curve() {
    let c = arch();
    for t in [0., 1., 0.5] {
        let r = with_endpoints(&c, c.point(t));
        assert!(!r.is_empty());
        assert_near(r.min_square_distance(), 0., TOL);
    }
}

#[test]
fn linear_bezier() {
    let c = bezier(&[(0., 0., 0.), (10., 0., 0.)]);
    assert_near(with_endpoints(&c, p(5., 3., 0.)).min_square_distance(), 9., TOL);
    assert_near(with_endpoints(&c, p(0., 0., 0.)).min_square_distance(), 0., TOL);
    assert_near(with_endpoints(&c, p(-5., 0., 0.)).min_square_distance(), 25., TOL);
}

#[test]
fn quadratic_bezier() {
    let c = bezier(&[(-2., 0., 0.), (0., 2., 0.), (2., 0., 0.)]);
    let r = with_endpoints(&c, p(0., 3., 0.));
    assert!(r.iter().any(|e| (e.parameter - 0.5).abs() < 0.01), "extremum at the apex");

    let r = with_endpoints(&c, p(0., 0., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    assert!(r.min_square_distance().sqrt() < 2.);

    let c = bezier(&[(0., 0., 0.), (2., 4., 0.), (6., 0., 0.)]);
    let q = p(3., 3., 0.);
    assert_consistent(&c, q, &with_endpoints(&c, q));
}

#[test]
fn cubic_and_quartic_shapes() {
    let cases = [
        (bezier(&[(0., 0., 0.), (0., 4., 0.), (4., -4., 0.), (4., 0., 0.)]), p(2., 0., 0.)),
        (bezier(&[(0., 0., 0.), (4., 3., 0.), (-2., 3., 0.), (2., 0., 0.)]), p(1., 1.5, 0.)),
        (bezier(&[(0., 0., 0.), (1., 3., 0.), (3., 3., 0.), (4., 0., 0.)]), p(2., -1., 0.)),
        (bezier(&[(0., 0., 0.), (1., 2., 0.), (2., -1., 0.), (3., 2., 0.), (4., 0., 0.)]), p(2., 1., 0.)),
        (bezier(&[(0., 0., 0.), (1., 1., 1.), (2., 1., 2.), (3., 0., 3.)]), p(1.5, 0.5, 1.5)),
    ];
    for (c, q) in &cases {
        assert_consistent(c, *q, &with_endpoints(c, *q));
    }
    let quartic = &cases[3].0;
    assert_near(with_endpoints(quartic, quartic.point(0.25)).min_square_distance(), 0., TOL);
    let c = bezier(&[(0., 0., 0.), (1., 2., 1.), (2., 2., 2.), (3., 0., 3.)]);
    assert_near(with_endpoints(&c, c.point(0.3)).min_square_distance(), 0., TOL);
}

#[test]
fn rational_quarter_circle_from_center_is_infinite() {
    let c = rational_quadratic(p(1., 0., 0.), p(1., 1., 0.), p(0., 1., 0.), 1. / 2f64.sqrt());
    let r = with_endpoints(&c, p(0., 0., 0.));
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap(), 1., 0.01);
}

#[test]
fn rational_high_weight() {
    let c = rational_quadratic(p(0., 0., 0.), p(1., 2., 0.), p(2., 0., 0.), 3.);
    let q = p(1., 1.5, 0.);
    let r = with_endpoints(&c, q);
    let min = r.min().expect("no extrema");
    let on = c.point(min.parameter);
    assert_near(min.point.distance(on), 0., TOL);
    assert_near(min.square_distance, q.distance(on).powi(2), TOL);
}

#[test]
fn bounded_ranges() {
    let c = arch();
    let r = with_endpoints(&c, c.point(0.5));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance(), 0., TOL);

    let r = PointCurveSearch::with_domain(&c, 0., 0.3).find_with_endpoints(c.point(0.5), TOL).unwrap();
    assert!(!r.is_empty());
    assert!(r.iter().all(|e| (-0.001..=0.301).contains(&e.parameter)));
}

#[test]
fn constant_curve_is_infinite() {
    let c = bezier(&[(1., 1., 1.); 4]);
    let r = with_endpoints(&c, p(2., 1., 1.));
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap(), 1., TOL);
}

#[test]
fn singleton_domain_returns_the_sole_point() {
    let c = bezier(&[(0., 0., 0.), (1., 2., 0.), (3., 0., 0.)]);
    let r = PointCurveSearch::with_domain(&c, 0.4, 0.4).find_with_endpoints(p(5., 5., 0.), TOL).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, 0.4, TOL);
    assert!(r[0].is_minimum && r[0].is_maximum);
}

#[test]
fn repeated_end_pole_does_not_fail() {
    let c = bezier(&[(0., 0., 0.), (0., 0., 0.), (2., 1., 0.), (3., 0., 0.)]);
    let r = with_endpoints(&c, p(0., 0., 0.));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance(), 0., TOL);
}

#[test]
fn closed_curve_seam_has_a_single_representative() {
    let c = bezier(&[(1., 0., 0.), (1., 1., 0.), (1., -1., 0.), (1., 0., 0.)]);
    let r = with_endpoints(&c, p(2., 0., 0.));
    let at_seam = r
        .iter()
        .filter(|e| (e.parameter - c.first_parameter()).abs() <= TOL || (e.parameter - c.last_parameter()).abs() <= TOL)
        .count();
    assert_eq!(at_seam, 1);
}

#[test]
fn far_and_near_points() {
    let c = bezier(&[(0., 0., 0.), (1., 1., 0.), (2., 1., 0.), (3., 0., 0.)]);
    let r = with_endpoints(&c, p(100., 100., 100.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    assert!(r.min_square_distance().sqrt() > 170.);

    let on = c.point(0.5);
    let r = with_endpoints(&c, p(on.x, on.y + 1e-8, on.z));
    assert!(!r.is_empty());
    assert_near(r.min_square_distance(), 0., 1e-12);

    let q = p(1.5, 2., 0.);
    assert_consistent(&c, q, &with_endpoints(&c, q));
}

#[test]
fn interior_extrema_are_perpendicular() {
    let c = arch();
    let q = p(2., 3., 0.);
    let r = with_endpoints(&c, q);
    let interior: Vec<_> = r.iter().filter(|e| (0.01..=0.99).contains(&e.parameter)).collect();
    let perpendicular = interior
        .iter()
        .filter(|e| {
            let (on, tangent) = c.d1(e.parameter);
            (q - on).dot(tangent).abs() < 0.1
        })
        .count();
    if !interior.is_empty() {
        assert!(perpendicular >= 1);
    }
}

#[test]
fn planar_case_matches_2d() {
    let ys = [0., 3., -1., 0.];
    let c3 = Curve::Bezier(BezierCurve::new((0..4).map(|i| p(i as f64, ys[i], 2.)).collect()));
    let c2 = Curve2::Bezier(BezierCurve2::new((0..4).map(|i| p2(i as f64, ys[i])).collect()));
    let r3 = PointCurveSearch::new(&c3).find(p(1.3, 1., 7.), TOL).unwrap();
    let r2 = PointCurveSearch::new(&c2).find(p2(1.3, 1.), TOL).unwrap();
    assert_eq!(r3.len(), r2.len());
    for (e3, e2) in r3.iter().zip(&r2) {
        assert_near(e3.parameter, e2.parameter, TOL);
        assert_near(e3.square_distance, e2.square_distance + 25., TOL);
    }
}

#[test]
fn result_follows_an_edited_pole() {
    let mut bezier = BezierCurve::new(vec![p(0., 0., 1.), p(2., 3., 1.), p(4., 0., 1.)]);
    let q = p(2., 2., 4.);
    let before = PointCurveSearch::new(&Curve::Bezier(bezier.clone())).find(q, TOL).unwrap();
    let t_before = before.min().unwrap().parameter;
    bezier.set_pole(1, p(3.5, -4., 1.));
    let c = Curve::Bezier(bezier);
    let after = PointCurveSearch::new(&c).find(q, TOL).unwrap();
    let min = after.min().expect("no extrema");
    assert!((min.parameter - t_before).abs() > 1e-3);
    assert_near(min.point.distance(c.point(min.parameter)), 0., 1e-12);
}

#[test]
fn small_non_constant_curve_is_not_infinite() {
    let length = 1.0e-3;
    let c = bezier(&(0..12).map(|i| (i as f64 / 11. * length, 0., 0.)).collect::<Vec<_>>());
    let r = PointCurveSearch::with_domain(&c, 0., 1.).find(p(0.5 * length, 1., 0.), 1.0e-10).unwrap();
    assert!(!r.is_infinite());
    assert_near(r.min().expect("no extrema").parameter, 0.5, 1.0e-7);
}
