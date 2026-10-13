//! GeomConvert_Test.cxx, GeomConvert_CompCurveToBSplineCurve_Test.cxx,
//! Geom2dConvert_CompCurveToBSplineCurve_Test.cxx.

use std::f64::consts::PI;

use kentron::geometry::convert::{self, Concatenation};
use kentron::geometry::curve::{BSplineCurve, Circle, Curve, Line, TrimmedCurve};
use kentron::geometry::curve2::{Circle2, Curve2, TrimmedCurve2};
use kentron::geometry::surface::{Plane, Surface, TrimmedSurface};
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Direction, Frame, Frame2, Point};
use kentron::tolerance::{LINEAR, PARAMETRIC};

use crate::{assert_point_near, p};

// ---- GeomConvert ----

#[test]
fn circle_to_bspline() {
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 5.0));
    let bspline = convert::curve_to_bspline(&circle).unwrap();
    assert_point_near(bspline.point(0.0), circle.point(0.0), LINEAR);
}

#[test]
fn line_to_bspline() {
    let line = Curve::Line(Line::new(p(0.0, 0.0, 0.0), Direction::X));
    let segment = Curve::Trimmed(TrimmedCurve::new(line, 0.0, 10.0));
    let bspline = convert::curve_to_bspline(&segment).unwrap();
    assert_point_near(bspline.point(5.0), p(5.0, 0.0, 0.0), LINEAR);
}

#[test]
fn plane_to_bspline_surface() {
    let patch = Surface::Trimmed(TrimmedSurface::new(Surface::Plane(Plane::XY), 0.0, 10.0, 0.0, 10.0));
    let bspline = convert::surface_to_bspline(&patch).unwrap();
    assert_point_near(bspline.point(5.0, 5.0), p(5.0, 5.0, 0.0), LINEAR);
}

// ---- GeomConvert_CompCurveToBSplineCurve ----

fn cubic(poles: [(f64, f64); 4]) -> BSplineCurve {
    BSplineCurve::new(
        poles.iter().map(|&(x, y)| p(x, y, 0.0)).collect(),
        vec![0.0, 1.0],
        vec![4, 4],
        3,
    )
}

fn concat(first: BSplineCurve) -> Concatenation {
    Concatenation::from_curve(&Curve::BSpline(first))
}

#[track_caller]
fn assert_ends(chain: &Concatenation, start: Point, end: Point) {
    let result = chain.result().expect("result curve");
    assert_point_near(result.start_point(), start, LINEAR);
    assert_point_near(result.end_point(), end, LINEAR);
}

#[test]
fn concatenate_clamped_bsplines() {
    let mut chain = concat(cubic([(0., 0.), (1., 1.), (2., 1.), (3., 0.)]));
    let second = cubic([(3., 0.), (4., -1.), (5., -1.), (6., 0.)]);
    assert!(chain.add(&Curve::BSpline(second), LINEAR));
    assert_ends(&chain, p(0., 0., 0.), p(6., 0., 0.));
}

#[test]
fn concatenate_trimmed_circle_arcs() {
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 5.0));
    let arc1 = TrimmedCurve::new(circle.clone(), 0.0, PI / 2.0);
    let arc2 = TrimmedCurve::new(circle, PI / 2.0, PI);
    assert_point_near(arc1.end_point(), arc2.start_point(), LINEAR);

    let mut chain = Concatenation::from_curve(&Curve::Trimmed(arc1));
    assert!(chain.add(&Curve::Trimmed(arc2), LINEAR));
    assert_ends(&chain, p(5., 0., 0.), p(-5., 0., 0.));
}

#[test]
fn concatenate_with_reversal() {
    let mut chain = concat(cubic([(0., 0.), (1., 1.), (2., 1.), (3., 0.)]));
    let second = cubic([(6., 0.), (5., -1.), (4., -1.), (3., 0.)]);
    assert!(chain.add(&Curve::BSpline(second), LINEAR));
    assert_ends(&chain, p(0., 0., 0.), p(6., 0., 0.));
}

#[test]
fn fails_for_disjoint_curves() {
    let mut chain = concat(cubic([(0., 0.), (1., 1.), (2., 1.), (3., 0.)]));
    let second = cubic([(10., 0.), (11., 1.), (12., 1.), (13., 0.)]);
    assert!(!chain.add(&Curve::BSpline(second), LINEAR));
}

/// Bug30007: the curves must be matched by their actual end points, not their poles.
#[test]
fn concatenate_non_clamped_bspline() {
    let first = BSplineCurve::new(
        [(0., 0.), (1., 2.), (2., 2.), (3., 2.), (4., 2.), (5., 0.)]
            .iter()
            .map(|&(x, y)| p(x, y, 0.0))
            .collect(),
        vec![0.0, 0.33, 0.67, 1.0],
        vec![4, 1, 1, 4],
        3,
    );
    let first_end = first.end_point();
    let first_start = first.start_point();
    let second = BSplineCurve::new(
        vec![first_end, p(6., -1., 0.), p(7., -1., 0.), p(8., 0., 0.)],
        vec![0.0, 1.0],
        vec![4, 4],
        3,
    );
    let mut chain = concat(first);
    assert!(chain.add(&Curve::BSpline(second), LINEAR));
    assert_ends(&chain, first_start, p(8., 0., 0.));
}

#[test]
fn prepend_curve() {
    let mut chain = concat(cubic([(3., 0.), (4., 1.), (5., 1.), (6., 0.)]));
    let second = cubic([(0., 0.), (1., -1.), (2., -1.), (3., 0.)]);
    assert!(chain.add(&Curve::BSpline(second), LINEAR));
    assert_ends(&chain, p(0., 0., 0.), p(6., 0., 0.));
}

#[test]
fn empty_initial_curve() {
    let mut chain = Concatenation::new();
    assert!(chain.add(&Curve::BSpline(cubic([(0., 0.), (1., 1.), (2., 1.), (3., 0.)])), LINEAR));
    assert!(chain.result().is_some());
}

#[test]
fn clear_and_reuse() {
    let curve = Curve::BSpline(cubic([(0., 0.), (1., 1.), (2., 1.), (3., 0.)]));
    let mut chain = Concatenation::from_curve(&curve);
    assert!(chain.result().is_some());
    chain.clear();
    assert!(chain.result().is_none());
    assert!(chain.add(&curve, LINEAR));
    assert!(chain.result().is_some());
}

// ---- Geom2dConvert_CompCurveToBSplineCurve ----

/// OCC30747: a closed contour assembled from circle arcs gives a closed curve.
#[test]
fn closed_contour_from_circle_arcs_2d() {
    let circle = Curve2::Circle(Circle2::new(Frame2::WORLD, 50.0));
    let (first, last) = (circle.first_parameter(), circle.last_parameter());
    let count = 10;
    let delta = (first + last) / count as f64;
    let mut chain =
        Concatenation::from_curve(&Curve2::Trimmed(TrimmedCurve2::new(circle.clone(), first, delta)));
    for i in 1..count {
        let end = if i == count - 1 { first } else { (i + 1) as f64 * delta };
        let arc = TrimmedCurve2::new(circle.clone(), i as f64 * delta, end);
        chain.add(&Curve2::Trimmed(arc), PARAMETRIC);
    }
    assert!(chain.result().expect("result curve").is_closed());
}
