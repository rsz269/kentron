//! ExtremaPC_Comparison_Test.cxx: `PointCurveSearch` against the projection query
//! `extrema::point_curve` (OCCT compares ExtremaPC with its older Extrema_ExtPC).

use kentron::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola,
    TrimmedCurve,
};
use kentron::geometry::extrema::PointCurveSearch;
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Point};

use crate::{assert_near, p, projection_distances, Rng};

const TOL: f64 = 1.0e-6;
const TOL_REL: f64 = 1.0e-4;

fn trimmed(c: Curve, t1: f64, t2: f64) -> Curve {
    Curve::Trimmed(TrimmedCurve::new(c, t1, t2))
}

fn circle10() -> Curve {
    Curve::Circle(Circle::new(Frame::WORLD, 10.))
}

fn bspline(poles: &[(f64, f64, f64)], knots: &[f64], mults: &[usize], degree: usize) -> Curve {
    Curve::BSpline(BSplineCurve::new(
        poles.iter().map(|&(x, y, z)| p(x, y, z)).collect(),
        knots.to_vec(),
        mults.to_vec(),
        degree,
    ))
}

fn cubic_arch() -> Curve {
    bspline(&[(0., 0., 0.), (1., 2., 0.), (2., 2., 0.), (3., 0., 0.)], &[0., 1.], &[4, 4], 3)
}

/// The search's closest distance matches the projections' (relatively, when not tiny).
#[track_caller]
fn compare_min_distance(c: &Curve, q: Point) -> f64 {
    let r = PointCurveSearch::new(c).find(q, TOL).unwrap();
    let reference = projection_distances(q, c).map(|(near, _)| near);
    let found = r.min().map(|e| e.distance());
    match (found, reference) {
        (None, None) => f64::INFINITY,
        (Some(new), Some(old)) => {
            if old > TOL {
                assert!((new - old).abs() / old < TOL_REL, "{new} vs {old}");
            } else {
                assert_near(new, old, TOL);
            }
            new
        }
        (found, reference) => panic!("search: {found:?}, projections: {reference:?}"),
    }
}

#[test]
fn line() {
    let c = trimmed(Curve::Line(Line::new(p(0., 0., 0.), Direction::X)), -100., 100.);
    compare_min_distance(&c, p(5., 0., 0.));
    compare_min_distance(&c, p(5., 3., 4.));
}

#[test]
fn circle() {
    let c = circle10();
    for q in [p(20., 0., 0.), p(3., 0., 0.), p(15., 0., 5.)] {
        compare_min_distance(&c, q);
    }
    assert!(compare_min_distance(&c, p(10., 0., 0.)) < TOL);
}

#[test]
fn ellipse() {
    let c = Curve::Ellipse(Ellipse::new(Frame::WORLD, 20., 10.));
    compare_min_distance(&c, p(30., 0., 0.));
    compare_min_distance(&c, p(15., 12., 0.));
}

#[test]
fn parabola_and_hyperbola() {
    compare_min_distance(&trimmed(Curve::Parabola(Parabola::new(Frame::WORLD, 2.)), -10., 10.), p(1., 2., 0.));
    compare_min_distance(&trimmed(Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 10., 5.)), -2., 2.), p(15., 3., 0.));
}

#[test]
fn bspline_and_bezier() {
    let c = cubic_arch();
    compare_min_distance(&c, p(1.5, 3., 0.));
    assert!(compare_min_distance(&c, c.point(0.5)) < TOL);
    let bezier = Curve::Bezier(BezierCurve::new(vec![p(0., 0., 0.), p(1., 3., 0.), p(3., 3., 0.), p(4., 0., 0.)]));
    compare_min_distance(&bezier, p(2., 4., 0.));
}

#[test]
fn multi_span_bspline_finds_a_point_on_the_curve() {
    let c = bspline(
        &[(0., 0., 0.), (1., 3., 0.), (2., 0., 0.), (3., 3., 0.), (4., 0., 0.), (5., 3., 0.)],
        &[0., 0.33, 0.67, 1.],
        &[3, 1, 2, 3],
        2,
    );
    let r = PointCurveSearch::new(&c).find(p(3., 4., 0.), TOL).unwrap();
    let min = r.min().expect("no extrema");
    assert_near(c.point(min.parameter).distance(min.point), 0., 1.0e-6);
}

#[test]
fn offset_circle() {
    compare_min_distance(&Curve::Offset(OffsetCurve::new(circle10(), 5., Direction::Z)), p(25., 0., 0.));
}

/// The search is never worse than the projections.
#[track_caller]
fn assert_not_worse(c: &Curve, rng: &mut Rng, (lo, hi): (f64, f64), count: usize) {
    let mut failures = 0;
    for _ in 0..count {
        let q = p(rng.uniform(lo, hi), rng.uniform(lo, hi), rng.uniform(lo, hi));
        let r = PointCurveSearch::new(c).find(q, TOL).unwrap();
        if let (Some(new), Some((old, _))) = (r.min(), projection_distances(q, c)) {
            if new.distance() > old * (1. + TOL_REL) + TOL {
                failures += 1;
            }
        }
    }
    assert_eq!(failures, 0, "{failures} of {count} points gave a worse minimum");
}

#[test]
fn random_points() {
    assert_not_worse(&circle10(), &mut Rng::new(12345), (-30., 30.), 100);
    let c = bspline(
        &[(0., 0., 0.), (1., 2., 1.), (2., 1., 2.), (3., 3., 1.), (4., 0., 0.)],
        &[0., 0.33, 0.67, 1.],
        &[3, 1, 1, 3],
        2,
    );
    assert_not_worse(&c, &mut Rng::new(54321), (-2., 7.), 100);
}
