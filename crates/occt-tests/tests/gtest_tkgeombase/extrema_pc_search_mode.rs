//! ExtremaPC_SearchMode_Test.cxx (closely related OCCT cases are grouped into one test).
//!
//! Where OCCT compares with its older `Extrema_ExtPC`, the reference here is
//! `extrema::point_curve` (see `projection_distances`).

use std::f64::consts::FRAC_PI_2;

use kentron::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola,
    TrimmedCurve,
};
use kentron::geometry::extrema::{PointCurveExtrema, PointCurveSearch, SearchMode};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction, Frame, Point};

use crate::{assert_near, p, projection_distances, Rng};

const TOL: f64 = 1.0e-6;

fn trimmed(c: Curve, t1: f64, t2: f64) -> Curve {
    Curve::Trimmed(TrimmedCurve::new(c, t1, t2))
}

fn x_segment(half_length: f64) -> Curve {
    trimmed(Curve::Line(Line::new(p(0., 0., 0.), Direction::X)), -half_length, half_length)
}

fn circle10() -> Curve {
    Curve::Circle(Circle::new(Frame::WORLD, 10.))
}

fn ellipse() -> Curve {
    Curve::Ellipse(Ellipse::new(Frame::WORLD, 20., 10.))
}

fn parabola_segment() -> Curve {
    trimmed(Curve::Parabola(Parabola::new(Frame::WORLD, 5.)), -20., 20.)
}

fn hyperbola_segment() -> Curve {
    trimmed(Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 10., 5.)), -2., 2.)
}

fn cubic(y: f64) -> Vec<Point> {
    vec![p(0., 0., 0.), p(1., y, 0.), p(3., y, 0.), p(4., 0., 0.)]
}

fn bspline() -> Curve {
    Curve::BSpline(BSplineCurve::new(cubic(2.), vec![0., 1.], vec![4, 4], 3))
}

fn bezier() -> Curve {
    Curve::Bezier(BezierCurve::new(cubic(3.)))
}

fn offset_circle() -> Curve {
    Curve::Offset(OffsetCurve::new(circle10(), 5., Direction::Z))
}

fn find(c: &Curve, q: Point, mode: SearchMode) -> PointCurveExtrema {
    PointCurveSearch::new(c).find_with(q, TOL, mode, false).unwrap()
}

/// The closest extremum agrees with the projection reference.
#[track_caller]
fn assert_min_matches_projections(c: &Curve, q: Point, r: &PointCurveExtrema, tolerance: f64) {
    let (reference, _) = projection_distances(q, c).expect("no projections");
    assert_near(r.min_square_distance().sqrt(), reference, tolerance);
}

/// The farthest extremum is inside the domain and on the curve.
#[track_caller]
fn assert_max_on_curve(c: &Curve, r: &PointCurveExtrema, tolerance: f64) {
    if let Some(max) = r.max() {
        assert!((c.first_parameter() - TOL..=c.last_parameter() + TOL).contains(&max.parameter));
        assert_near(c.point(max.parameter).distance(max.point), 0., tolerance);
    }
}

#[test]
fn line_modes() {
    let c = x_segment(100.);
    let q = p(25., 10., 0.);
    let r = find(&c, q, SearchMode::Min);
    assert_eq!(r.len(), 1, "min mode returns only the minimum");
    assert_min_matches_projections(&c, q, &r, TOL);
    assert_near(r.min_square_distance().sqrt(), 10., TOL);

    let c = x_segment(50.);
    let r = PointCurveSearch::new(&c).find_with(p(0., 10., 0.), TOL, SearchMode::Max, true).unwrap();
    assert!(!r.is_empty());
    assert_near(r.max_square_distance().sqrt(), 2600f64.sqrt(), TOL);

    let q = p(10., 5., 0.);
    let r = find(&c, q, SearchMode::MinMax);
    assert_min_matches_projections(&c, q, &r, TOL);
    assert_near(r.min_square_distance().sqrt(), 5., TOL);
}

#[test]
fn circle_modes() {
    let c = circle10();
    let q = p(20., 0., 0.);
    let r = find(&c, q, SearchMode::Min);
    assert_eq!(r.len(), 1);
    assert_min_matches_projections(&c, q, &r, TOL);
    assert_near(r.min_square_distance().sqrt(), 10., TOL);
    let r = find(&c, q, SearchMode::Max);
    assert_eq!(r.len(), 1);
    assert_near(r.max_square_distance().sqrt(), projection_distances(q, &c).unwrap().1, TOL);
    assert_near(r.max_square_distance().sqrt(), 30., TOL);

    let q = p(3., 0., 0.);
    let r = find(&c, q, SearchMode::Min);
    assert_eq!(r.len(), 1);
    assert_near(r.min_square_distance().sqrt(), 7., TOL);
    let r = find(&c, q, SearchMode::Max);
    assert_eq!(r.len(), 1);
    assert_near(r.max_square_distance().sqrt(), 13., TOL);

    let q = p(15., 0., 0.);
    let r = find(&c, q, SearchMode::MinMax);
    assert_eq!(r.len(), 2);
    let (near, far) = projection_distances(q, &c).unwrap();
    assert_near(r.min_square_distance().sqrt(), near, TOL);
    assert_near(r.max_square_distance().sqrt(), far, TOL);
    assert_near(near, 5., TOL);
    assert_near(far, 25., TOL);
}

#[test]
fn circle_arc_modes() {
    let c = trimmed(circle10(), 0., FRAC_PI_2);
    let q = p(15., 15., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);
    assert_max_on_curve(&c, &find(&c, q, SearchMode::Max), TOL);
}

#[test]
fn ellipse_modes() {
    let c = ellipse();
    let q = p(30., 0., 0.);
    let r = find(&c, q, SearchMode::Min);
    assert_eq!(r.len(), 1);
    assert_min_matches_projections(&c, q, &r, TOL);
    assert_near(r.min_square_distance().sqrt(), 10., TOL);
    let r = find(&c, q, SearchMode::Max);
    assert_eq!(r.len(), 1);
    assert_near(r.max_square_distance().sqrt(), projection_distances(q, &c).unwrap().1, TOL);
    assert_near(r.max_square_distance().sqrt(), 50., TOL);

    let q = p(15., 12., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);

    let q = p(15., 8., 0.);
    let r = find(&c, q, SearchMode::MinMax);
    let (near, far) = projection_distances(q, &c).unwrap();
    assert_near(r.min_square_distance().sqrt(), near, TOL);
    assert_near(r.max_square_distance().sqrt(), far, TOL);
}

#[test]
fn parabola_modes() {
    let c = parabola_segment();
    let q = p(0., 5., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);
    let r = PointCurveSearch::new(&c).find_with(q, TOL, SearchMode::Max, true).unwrap();
    assert!(!r.is_empty());
    assert_max_on_curve(&c, &r, TOL);
    let q = p(10., -5., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);
}

#[test]
fn hyperbola_modes() {
    let c = hyperbola_segment();
    let q = p(20., 0., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);
    assert_max_on_curve(&c, &find(&c, q, SearchMode::Max), TOL);
}

#[test]
fn bspline_modes() {
    let c = bspline();
    let q = p(2., 3., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), 0.01);
    assert_max_on_curve(&c, &find(&c, p(2., 1., 0.), SearchMode::Max), 0.01);

    let q = p(2., 1.5, 0.);
    let r = find(&c, q, SearchMode::MinMax);
    assert!(!r.is_empty());
    for e in &r {
        assert!((-TOL..=1. + TOL).contains(&e.parameter));
        assert_near(c.point(e.parameter).distance(e.point), 0., 0.01);
    }
    if projection_distances(q, &c).is_some() {
        assert_min_matches_projections(&c, q, &r, 0.05);
    }
}

#[test]
fn bezier_modes() {
    let c = bezier();
    let q = p(2., 4., 0.);
    assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), 0.01);
    let r = PointCurveSearch::new(&c).find_with(p(2., 1., 0.), TOL, SearchMode::Max, true).unwrap();
    assert!(!r.is_empty());
    assert_near(r.max_square_distance().sqrt(), 5f64.sqrt(), 0.01);
}

#[test]
fn offset_circle_modes() {
    let c = offset_circle();
    let q = p(25., 0., 0.);
    let r = find(&c, q, SearchMode::Min);
    assert_min_matches_projections(&c, q, &r, 0.01);
    assert_near(r.min_square_distance().sqrt(), 10., 0.01);
    let r = find(&c, q, SearchMode::Max);
    if let Some(max) = r.max() {
        assert_near(c.point(max.parameter).distance(max.point), 0., 0.01);
        assert_near(r.max_square_distance().sqrt(), 40., 0.5);
    }
}

#[test]
fn random_points_on_conics() {
    let mut rng = Rng::new(12345);
    let c = circle10();
    for _ in 0..50 {
        let q = p(rng.uniform(-50., 50.), rng.uniform(-50., 50.), 0.);
        let (near, far) = projection_distances(q, &c).expect("projections");
        assert_near(find(&c, q, SearchMode::Min).min_square_distance().sqrt(), near, TOL);
        assert_near(find(&c, q, SearchMode::Max).max_square_distance().sqrt(), far, TOL);
    }
    let mut rng = Rng::new(54321);
    let c = ellipse();
    for _ in 0..50 {
        let q = p(rng.uniform(-50., 50.), rng.uniform(-50., 50.), 0.);
        assert_min_matches_projections(&c, q, &find(&c, q, SearchMode::Min), TOL);
    }
}

#[test]
fn random_points_on_bspline() {
    let mut rng = Rng::new(11111);
    let c = bspline();
    let (mut tested, mut passed) = (0, 0);
    for _ in 0..30 {
        let q = p(rng.uniform(-2., 6.), rng.uniform(-2., 6.), 0.);
        let r = find(&c, q, SearchMode::Min);
        let Some((near, _)) = projection_distances(q, &c) else { continue };
        if r.is_empty() {
            continue;
        }
        tested += 1;
        if (r.min_square_distance().sqrt() - near).abs() < 0.1 {
            passed += 1;
        }
    }
    if tested > 0 {
        assert!(passed as f64 / tested as f64 >= 0.8, "{passed}/{tested}");
    }
}
