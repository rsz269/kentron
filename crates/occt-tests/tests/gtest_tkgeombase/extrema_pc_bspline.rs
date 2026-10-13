//! ExtremaPC_BSplineCurve_Test.cxx (closely related OCCT cases are grouped into one test).
//!
//! Not ported: `GridUtilitiesRejectInvalidSamplesAndShortOpenClosure` (OCCT's internal grid
//! evaluator).

use kentron::geometry::curve::{BSplineCurve, Curve};
use kentron::geometry::curve2::{BSplineCurve2, Curve2};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch, SearchMode};
use kentron::geometry::ParametricCurve;
use kentron::math::Point;

use crate::{assert_near, p, p2};

const TOL: f64 = 1.0e-6;

fn bspline(poles: &[(f64, f64, f64)], knots: &[f64], mults: &[usize], degree: usize) -> Curve {
    Curve::BSpline(BSplineCurve::new(
        poles.iter().map(|&(x, y, z)| p(x, y, z)).collect(),
        knots.to_vec(),
        mults.to_vec(),
        degree,
    ))
}

/// A single-span cubic arch from (0,0) to (3,0).
fn cubic() -> Curve {
    bspline(&[(0., 0., 0.), (1., 2., 0.), (2., 2., 0.), (3., 0., 0.)], &[0., 1.], &[4, 4], 3)
}

fn linear() -> Curve {
    bspline(&[(0., 0., 0.), (10., 0., 0.)], &[0., 1.], &[2, 2], 1)
}

fn quadratic() -> Curve {
    bspline(&[(0., 0., 0.), (5., 10., 0.), (10., 0., 0.)], &[0., 1.], &[3, 3], 2)
}

fn multi_span() -> Curve {
    bspline(
        &[(0., 0., 0.), (1., 2., 0.), (2., 0., 0.), (3., 2., 0.), (4., 0., 0.)],
        &[0., 0.33, 0.67, 1.],
        &[3, 1, 1, 3],
        2,
    )
}

fn spatial() -> Curve {
    bspline(&[(0., 0., 0.), (1., 1., 1.), (2., 1., 2.), (3., 0., 3.)], &[0., 1.], &[4, 4], 3)
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
fn point_on_curve() {
    let c = cubic();
    for q in [p(0., 0., 0.), p(3., 0., 0.), c.point(0.5)] {
        let r = with_endpoints(&c, q);
        assert!(!r.is_empty());
        assert_near(r.min_square_distance(), 0., TOL);
    }
    let c = quadratic();
    assert_near(with_endpoints(&c, c.point(0.5)).min_square_distance(), 0., TOL);
    let c = spatial();
    assert_near(with_endpoints(&c, c.point(0.5)).min_square_distance(), 0., TOL);
}

#[test]
fn singleton_domain_returns_the_sole_point() {
    let c = cubic();
    let r = PointCurveSearch::with_domain(&c, 0.35, 0.35).find_with_endpoints(p(5., 5., 0.), TOL).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, 0.35, TOL);
    assert!(r[0].is_minimum && r[0].is_maximum);
}

#[test]
fn point_near_curve() {
    let c = cubic();
    let d2 = with_endpoints(&c, p(1.5, 3., 0.)).min_square_distance();
    assert!(d2 > 0. && d2 < 4.);

    let r = with_endpoints(&c, p(1.5, -1., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    assert!(r.min_square_distance() >= 1.);

    assert!(with_endpoints(&c, p(1.5, 2., 5.)).min_square_distance() >= 25.);
}

#[test]
fn linear_bspline() {
    let c = linear();
    assert_near(with_endpoints(&c, p(5., 0., 0.)).min_square_distance(), 0., TOL);
    assert_near(with_endpoints(&c, p(5., 3., 4.)).min_square_distance(), 25., TOL);
}

#[test]
fn quadratic_and_spatial() {
    let c = quadratic();
    let q = p(5., 8., 0.);
    assert_min_consistent(&c, q, &with_endpoints(&c, q));
    let c = spatial();
    let q = p(1.5, 2., 1.5);
    assert_min_consistent(&c, q, &with_endpoints(&c, q));
}

#[test]
fn multi_span_curve() {
    let c = multi_span();
    let r = with_endpoints(&c, p(2.5, 3., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    assert!(r.min_square_distance().sqrt() < 2.);

    let r = with_endpoints(&c, p(2.5, 1., 0.));
    assert!(!r.is_empty());
    for e in &r {
        assert_near(e.point.distance(c.point(e.parameter)), 0., TOL);
    }
}

#[test]
fn partial_ranges() {
    let c = cubic();
    for ((t1, t2), q) in [((0., 0.5), p(2.5, 1., 0.)), ((0.5, 1.), p(0.5, 1., 0.))] {
        let r = PointCurveSearch::with_domain(&c, t1, t2).find_with_endpoints(q, TOL).unwrap();
        assert!(!r.is_empty());
        assert!(r.iter().all(|e| (t1 - TOL..=t2 + TOL).contains(&e.parameter)));
    }
}

#[test]
fn projected_point_and_distance_are_consistent() {
    let c = cubic();
    let r = with_endpoints(&c, p(1.5, 3., 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);

    let q = p(2., 2.5, 1.);
    let r = with_endpoints(&c, q);
    let min = r.min().expect("no extrema");
    assert_near(min.square_distance, q.distance(min.point).powi(2), TOL);

    assert!(!with_endpoints(&c, p(1.5, 2., 0.)).is_empty());
}

#[test]
fn high_degree_spline() {
    let c = bspline(
        &[
            (0., 0., 0.),
            (0.5, 1., 0.),
            (1., -0.5, 0.),
            (1.5, 1.5, 0.),
            (2., 0., 0.),
            (2.5, 1., 0.),
            (3., -0.5, 0.),
            (3.5, 0.5, 0.),
        ],
        &[0., 1.],
        &[8, 8],
        7,
    );
    let r = with_endpoints(&c, p(1.75, 0.8, 0.));
    let min = r.min().expect("no extrema");
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    assert!(r.min_square_distance().sqrt() < 2.);
}

#[test]
fn c1_junction_point_on_curve() {
    let c = bspline(&[(0., 0., 0.), (1., 1., 0.), (2., 1., 0.), (3., 0., 0.)], &[0., 0.5, 1.], &[3, 1, 3], 2);
    let r = PointCurveSearch::new(&c).find(c.point(0.5), TOL).unwrap();
    assert_near(r.min().expect("no extrema").parameter, 0.5, TOL);
    assert_near(r.min_square_distance(), 0., TOL);
}

#[test]
fn clustered_c0_junctions() {
    let c = bspline(&[(0., 0., 0.), (1., 1., 0.), (2., 0., 0.), (3., 1., 0.)], &[0., 0.499, 0.501, 1.], &[2, 1, 1, 2], 1);
    let s = PointCurveSearch::new(&c);
    let q = p(1.5, 1.5, 0.);
    let at = |r: &PointCurveExtrema, t: f64| r.iter().find(|e| (e.parameter - t).abs() <= TOL).copied();

    let r = s.find(q, TOL).unwrap();
    assert!(at(&r, 0.499).expect("corner at 0.499").is_minimum);
    assert!(at(&r, 0.501).expect("corner at 0.501").is_maximum);

    let r = s.find_with(q, TOL, SearchMode::Min, false).unwrap();
    assert!(at(&r, 0.499).is_some() && at(&r, 0.501).is_none());

    let r = s.find_with(q, TOL, SearchMode::Max, false).unwrap();
    assert!(at(&r, 0.499).is_none() && at(&r, 0.501).is_some());
}

#[test]
fn planar_rational_case_matches_2d() {
    let xs: Vec<f64> = (1..=4).map(|i| 10. + i as f64).collect();
    let ys: Vec<f64> = (1..=4).map(|i| if i % 2 == 0 { 2. } else { -1. }).collect();
    let weights: Vec<f64> = (1..=4).map(|i| 0.5 + 0.25 * i as f64).collect();
    let c3 = Curve::BSpline(
        BSplineCurve::new((0..4).map(|i| p(xs[i], -4., ys[i])).collect(), vec![2., 6.], vec![4, 4], 3)
            .with_weights(weights.clone()),
    );
    let c2 = Curve2::BSpline(
        BSplineCurve2::new((0..4).map(|i| p2(xs[i], ys[i])).collect(), vec![2., 6.], vec![4, 4], 3)
            .with_weights(weights),
    );
    let r3 = PointCurveSearch::with_domain(&c3, 2.5, 5.5).find_with_endpoints(p(12., 3., 0.5), TOL).unwrap();
    let r2 = PointCurveSearch::with_domain(&c2, 2.5, 5.5).find_with_endpoints(p2(12., 0.5), TOL).unwrap();
    assert_eq!(r3.len(), r2.len());
    for (e3, e2) in r3.iter().zip(&r2) {
        assert_near(e3.parameter, e2.parameter, TOL);
        assert_near(e3.square_distance, e2.square_distance + 49., TOL);
    }
}

#[test]
fn non_planar_curve() {
    let c = bspline(&[(0., 0., 0.), (1., 2., 1.), (2., -1., 3.), (4., 0., -2.)], &[0., 1.], &[4, 4], 3);
    assert!(!PointCurveSearch::new(&c).find(p(1., 0.5, 2.), TOL).unwrap().is_empty());
}

#[test]
fn shifted_multi_period_domain_finds_every_minimum() {
    let poles = vec![p(1., 0., 0.), p(0.3, 1., 0.4), p(-0.8, 0.6, -0.2), p(-0.7, -0.7, 0.8), p(0.4, -0.9, -0.5)];
    let knots = (0..6).map(|i| 0.2 * i as f64).collect();
    let c = Curve::BSpline(BSplineCurve::new_periodic(poles, knots, vec![1; 6], 3));
    let period = c.period().expect("periodic");
    let target = c.first_parameter() + 0.17 * period;
    let s = PointCurveSearch::with_domain(&c, target - 0.1 * period, target + 3.1 * period);
    let r = s.find_with(c.point(target), 1.0e-9, SearchMode::Min, false).unwrap();
    let zero: Vec<_> = r.iter().filter(|e| e.square_distance <= 1.0e-16).collect();
    assert_eq!(zero.len(), 4);
    for (k, e) in zero.iter().enumerate() {
        assert_near(e.parameter, target + k as f64 * period, 1.0e-7);
    }
}
