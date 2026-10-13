//! ExtremaPC_ExtendedGeometry_Test.cxx: `PointCurveSearch` on varied placements and shapes,
//! checked against the projection query `extrema::point_curve` (OCCT compares with its
//! older Extrema_ExtPC).

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola,
    TrimmedCurve,
};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Point};

use crate::{assert_near, p, projection_distances, Rng};

const TOL: f64 = 1.0e-6;

fn trimmed(c: Curve, t1: f64, t2: f64) -> Curve {
    Curve::Trimmed(TrimmedCurve::new(c, t1, t2))
}

fn frame(origin: Point, normal: Direction) -> Frame {
    Frame::from_z(origin, normal)
}

fn circle(frame: Frame, radius: f64) -> Curve {
    Curve::Circle(Circle::new(frame, radius))
}

fn ellipse(frame: Frame, major: f64, minor: f64) -> Curve {
    Curve::Ellipse(Ellipse::new(frame, major, minor))
}

fn parabola(frame: Frame, focal: f64, half_range: f64) -> Curve {
    trimmed(Curve::Parabola(Parabola::new(frame, focal)), -half_range, half_range)
}

fn hyperbola(frame: Frame, major: f64, minor: f64, half_range: f64) -> Curve {
    trimmed(Curve::Hyperbola(Hyperbola::new(frame, major, minor)), -half_range, half_range)
}

fn bspline(poles: &[(f64, f64, f64)], knots: &[f64], mults: &[usize], degree: usize) -> Curve {
    Curve::BSpline(BSplineCurve::new(
        poles.iter().map(|&(x, y, z)| p(x, y, z)).collect(),
        knots.to_vec(),
        mults.to_vec(),
        degree,
    ))
}

fn bezier(poles: &[(f64, f64, f64)]) -> Curve {
    Curve::Bezier(BezierCurve::new(poles.iter().map(|&(x, y, z)| p(x, y, z)).collect()))
}

fn find(c: &Curve, q: Point) -> PointCurveExtrema {
    PointCurveSearch::new(c).find(q, TOL).unwrap()
}

/// The closest extremum matches the projections (when there are any).
#[track_caller]
fn compare_min(c: &Curve, q: Point, tolerance: f64) -> PointCurveExtrema {
    let r = find(c, q);
    if let Some((near, _)) = projection_distances(q, c) {
        assert!(!r.is_empty());
        assert_near(r.min_square_distance().sqrt(), near, tolerance);
    }
    r
}

/// Closest and farthest extrema match the projections.
#[track_caller]
fn compare_min_max(c: &Curve, q: Point, tolerance: f64) {
    let r = compare_min(c, q, tolerance);
    if let Some((_, far)) = projection_distances(q, c) {
        assert_near(r.max_square_distance().sqrt(), far, tolerance);
    }
}

#[test]
fn circles() {
    compare_min_max(&circle(frame(p(100., 200., 50.), Direction::Z), 25.), p(150., 200., 50.), TOL);
    let diagonal = Direction::new(1., 1., 0.);
    compare_min_max(&circle(frame(p(0., 0., 0.), diagonal), 15.), p(20., 20., 5.), TOL);
    compare_min_max(&circle(frame(p(0., 0., 0.), Direction::Y), 10.), p(15., 5., 0.), TOL);
    compare_min_max(&circle(frame(p(0., 0., 0.), Direction::X), 10.), p(5., 15., 0.), TOL);
    let skew = Direction::new(1., 1., 1.);
    compare_min_max(&circle(frame(p(10., 20., 30.), skew), 8.), p(15., 25., 35.), TOL);
    compare_min_max(&circle(Frame::WORLD, 0.001), p(0.002, 0., 0.), 1.0e-9);
    compare_min_max(&circle(Frame::WORLD, 1000.), p(1500., 0., 0.), TOL);
}

#[test]
fn ellipses() {
    compare_min_max(&ellipse(frame(p(50., 100., 25.), Direction::Z), 30., 15.), p(90., 100., 25.), TOL);
    compare_min_max(&ellipse(frame(p(0., 0., 0.), Direction::Y), 20., 8.), p(25., 10., 5.), TOL);
    compare_min_max(&ellipse(Frame::WORLD, 50., 2.), p(40., 5., 0.), TOL);
    compare_min_max(&ellipse(Frame::WORLD, 10.1, 10.), p(15., 8., 0.), TOL);
}

#[test]
fn lines() {
    let diagonal = trimmed(Curve::Line(Line::new(p(0., 0., 0.), Direction::new(1., 1., 1.))), -100., 100.);
    compare_min(&diagonal, p(10., 0., 0.), TOL);
    let skew = trimmed(Curve::Line(Line::new(p(10., 20., 30.), Direction::new(2., 3., 1.))), -50., 50.);
    compare_min(&skew, p(25., 35., 40.), TOL);
    let x = trimmed(Curve::Line(Line::new(p(0., 0., 0.), Direction::X)), -100., 100.);
    assert_near(find(&x, p(50., 0., 0.)).min_square_distance().sqrt(), 0., TOL);
}

#[test]
fn parabolas() {
    compare_min(&parabola(frame(p(10., 20., 0.), Direction::Z), 3., 15.), p(20., 30., 0.), TOL);
    compare_min(&parabola(frame(p(0., 0., 0.), Direction::Y), 5., 10.), p(5., 5., 10.), TOL);
    compare_min(&parabola(Frame::WORLD, 0.5, 5.), p(2., 3., 0.), TOL);
    compare_min(&parabola(Frame::WORLD, 50., 100.), p(50., 25., 0.), TOL);
}

#[test]
fn hyperbolas() {
    compare_min(&hyperbola(frame(p(50., 50., 0.), Direction::Z), 15., 8., 1.5), p(80., 50., 0.), TOL);
    compare_min(&hyperbola(frame(p(0., 0., 0.), Direction::Y), 10., 5., 2.), p(20., 5., 0.), TOL);
    let q = p(50., 10., 0.);
    let r = find(&hyperbola(Frame::WORLD, 30., 3., 1.), q);
    let min = r.min().expect("no extrema");
    assert_near(q.distance(min.point).powi(2), min.square_distance, TOL);
}

#[test]
fn bsplines() {
    let c = bspline(&[(0., 0., 0.), (1., 1., 0.), (2., 0., 0.), (3., 1., 0.)], &[0., 0.33, 0.67, 1.], &[2, 1, 1, 2], 1);
    compare_min(&c, p(1.5, 1.5, 0.), 0.01);
    let c = bspline(
        &[(0., 0., 0.), (1., 2., 0.), (2., 0., 0.), (3., 2., 0.), (4., 0., 0.)],
        &[0., 0.5, 1.],
        &[3, 2, 3],
        2,
    );
    compare_min(&c, p(2., 2.5, 0.), 0.01);
    let helix = bspline(
        &[(10., 0., 0.), (7., 7., 5.), (0., 10., 10.), (-7., 7., 15.), (-10., 0., 20.), (-7., -7., 25.), (0., -10., 30.)],
        &[0., 1.],
        &[7, 7],
        6,
    );
    compare_min(&helix, p(0., 0., 15.), 0.05);
}

#[test]
fn beziers() {
    compare_min(&bezier(&[(0., 0., 0.), (5., 10., 0.), (10., 0., 0.)]), p(5., 8., 0.), 0.01);
    let wavy = bezier(&[(0., 0., 0.), (1., 3., 0.), (2., -1., 0.), (3., 4., 0.), (4., 0., 0.), (5., 2., 0.), (6., 1., 0.)]);
    compare_min(&wavy, p(3., 3., 0.), 0.01);
    compare_min(&bezier(&[(0., 0., 0.), (3., 5., 2.), (6., 0., 4.), (9., 5., 6.)]), p(4.5, 3., 3.), 0.01);
}

#[test]
fn offset_curves() {
    let offset = |basis, d| Curve::Offset(OffsetCurve::new(basis, d, Direction::Z));
    compare_min(&offset(ellipse(Frame::WORLD, 20., 10.), 3.), p(30., 5., 0.), 0.05);
    compare_min(&offset(bezier(&[(0., 0., 0.), (2., 4., 0.), (4., 4., 0.), (6., 0., 0.)]), 1.), p(3., 5., 0.), 0.05);
    let r = compare_min(&offset(circle(Frame::WORLD, 15.), -5.), p(20., 0., 0.), 0.01);
    assert_near(r.min_square_distance().sqrt(), 10., 0.01);
}

#[test]
fn arcs() {
    compare_min(&trimmed(circle(Frame::WORLD, 10.), 0., FRAC_PI_2), p(15., 15., 0.), TOL);
    compare_min(&trimmed(circle(Frame::WORLD, 10.), PI, 1.5 * PI), p(-15., -5., 0.), TOL);
    compare_min(&trimmed(ellipse(Frame::WORLD, 20., 10.), 0., PI), p(0., 20., 0.), TOL);
}

#[test]
fn special_points() {
    let c = circle(Frame::WORLD, 10.);
    assert_near(find(&c, p(10., 0., 0.)).min_square_distance().sqrt(), 0., TOL);
    let b = bspline(&[(0., 0., 0.), (1., 2., 0.), (3., 2., 0.), (4., 0., 0.)], &[0., 1.], &[4, 4], 3);
    assert_near(find(&b, b.point(0.5)).min_square_distance().sqrt(), 0., 0.01);
    compare_min_max(&c, p(10000., 10000., 0.), TOL);
    assert_near(find(&c, p(10.0001, 0., 0.)).min_square_distance().sqrt(), 0.0001, 1.0e-5);
}

/// Counts the random points (in `[lo, hi)³`) whose closest extremum matches the
/// projections within `tolerance`; returns (matching, compared).
fn random_agreement(c: &Curve, seed: u64, (lo, hi): (f64, f64), tolerance: f64) -> (usize, usize) {
    let mut rng = Rng::new(seed);
    let (mut passed, mut compared) = (0, 0);
    for _ in 0..100 {
        let q = p(rng.uniform(lo, hi), rng.uniform(lo, hi), rng.uniform(lo, hi));
        let r = find(c, q);
        let Some((near, _)) = projection_distances(q, c) else { continue };
        if r.is_empty() {
            continue;
        }
        compared += 1;
        if (r.min_square_distance().sqrt() - near).abs() < tolerance {
            passed += 1;
        }
    }
    (passed, compared)
}

#[test]
fn random_points() {
    let (passed, _) = random_agreement(&circle(Frame::WORLD, 10.), 12345, (-100., 100.), TOL);
    assert!(passed >= 95, "{passed}/100");
    let (passed, _) = random_agreement(&ellipse(Frame::WORLD, 20., 10.), 54321, (-100., 100.), TOL);
    assert!(passed >= 95, "{passed}/100");
    let (passed, compared) = random_agreement(&parabola(Frame::WORLD, 5., 20.), 11111, (-50., 50.), 0.01);
    if compared > 0 {
        assert!(passed as f64 / compared as f64 >= 0.8, "{passed}/{compared}");
    }
}
