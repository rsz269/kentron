//! ExtremaPC2d_OffsetCurve_Test.cxx, ExtremaPC2d_Curve_Test.cxx, ExtremaPC2d_Test.cxx.
//!
//! Not ported: comparisons of OCCT's dispatcher with its per-type evaluators (kentron has
//! one search), null handles, result-storage identity, the evaluator caches
//! (`EndpointPathRefreshesPartitionAfterMutation`, `AnalyticalDispatchCopiesPrimitiveGeometry`),
//! and the internal helpers of ExtremaPC2d.hxx (`IsValidTolerance`, `IsFiniteParameter`,
//! `IsClosedPeriodicDomain`, `AddEndpointExtrema`).

use std::f64::consts::{FRAC_PI_2, PI};

use kentron::geometry::curve2::{
    BezierCurve2, Circle2, Curve2, Line2, OffsetCurve2, TrimmedCurve2,
};
use kentron::geometry::extrema::{
    ExtremaError, PointCurveExtrema, PointCurveSearch, PointExtremum, SearchMode,
};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction2, Frame2, Point2, Vector2};

use crate::{assert_near, p2};

const TOL: f64 = 1.0e-8;

fn circle(radius: f64) -> Curve2 {
    Curve2::Circle(Circle2::new(Frame2::WORLD, radius))
}

fn offset_circle(offset: f64) -> Curve2 {
    Curve2::Offset(OffsetCurve2::new(circle(3.), offset))
}

// ---- offset curves ----

#[test]
fn offset_circle_min_and_max() {
    let c = offset_circle(1.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find(p2(8., 0.), TOL).unwrap();
    assert_eq!(r.len(), 2);
    assert_near(r.min_square_distance(), 16., TOL);
    assert_near(r.max_square_distance(), 144., TOL);
}

#[test]
fn offset_circle_center_is_infinite() {
    let c = offset_circle(1.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find(p2(0., 0.), TOL).unwrap();
    assert!(r.is_infinite());
    assert_near(r.infinite_square_distance.unwrap(), 16., TOL);
}

#[test]
fn negative_offset_uses_the_other_side() {
    let c = offset_circle(-1.);
    let r = PointCurveSearch::with_domain(&c, 0., 2. * PI).find(p2(8., 0.), TOL).unwrap();
    assert_near(r.min_square_distance(), 36., TOL);
}

#[test]
fn offset_line_projection_is_exact() {
    let c = Curve2::Offset(OffsetCurve2::new(Curve2::Line(Line2::new(p2(0., 0.), Direction2::X)), 2.));
    let q = c.point(2.) + Vector2::new(0., 3.);
    let r = PointCurveSearch::with_domain(&c, 0., 4.).find(q, TOL).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, 2., TOL);
    assert_near(r[0].square_distance, 9., TOL);
}

#[test]
fn offset_arc_search_modes() {
    let c = offset_circle(1.);
    let s = PointCurveSearch::with_domain(&c, 0., FRAC_PI_2);
    let r = s.find_with(p2(8., -1.), TOL, SearchMode::Min, true).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, 0., TOL);
    let r = s.find_with(p2(8., -1.), TOL, SearchMode::Max, true).unwrap();
    assert_eq!(r.len(), 1);
    assert_near(r[0].parameter, FRAC_PI_2, TOL);
}

#[test]
fn offset_invalid_domain_and_tolerance() {
    let c = offset_circle(1.);
    assert_eq!(PointCurveSearch::with_domain(&c, 1., 0.).find(p2(0., 0.), TOL), Err(ExtremaError::InvalidInput));
    assert_eq!(PointCurveSearch::with_domain(&c, 0., 2. * PI).find(p2(0., 0.), 0.), Err(ExtremaError::InvalidInput));
}

// ---- curve kinds ----

#[test]
fn offset_search_modes() {
    let c = offset_circle(1.);
    let s = PointCurveSearch::with_domain(&c, 0., 2. * PI);
    let r = s.find_with(p2(8., 0.), TOL, SearchMode::Min, false).unwrap();
    assert_eq!(r.len(), 1);
    assert!(r[0].is_minimum);
    let r = s.find_with(p2(8., 0.), TOL, SearchMode::Max, false).unwrap();
    assert_eq!(r.len(), 1);
    assert!(r[0].is_maximum);
}

#[test]
fn trimmed_curve_uses_its_trim() {
    let c = Curve2::Trimmed(TrimmedCurve2::new(circle(3.), 0.2, 1.));
    let r = PointCurveSearch::new(&c).find_with_endpoints(p2(6., 0.), TOL).unwrap();
    assert!(r.iter().all(|e| (0.2 - TOL..=1. + TOL).contains(&e.parameter)));
    let r = PointCurveSearch::with_domain(&c, 0.5, 2.).find_with_endpoints(p2(6., 0.), TOL).unwrap();
    assert!(r.iter().all(|e| (0.5 - TOL..=1. + TOL).contains(&e.parameter)));
}

#[test]
fn edited_bezier_changes_the_result() {
    let mut bezier = BezierCurve2::new(vec![p2(0., 0.), p2(2., 3.), p2(4., 0.)]);
    let q = p2(2., 2.);
    let before = PointCurveSearch::new(&Curve2::Bezier(bezier.clone())).find(q, TOL).unwrap().min_square_distance();
    bezier.set_pole(1, p2(3.5, -4.));
    let c = Curve2::Bezier(bezier);
    let r = PointCurveSearch::new(&c).find(q, TOL).unwrap();
    assert!((r.min_square_distance() - before).abs() > 1.0e-4);
    let min = r.min().unwrap();
    assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
}

// ---- result ----

fn extremum(parameter: f64, square_distance: f64, is_minimum: bool) -> PointExtremum<Point2> {
    PointExtremum { parameter, point: p2(0., 0.), square_distance, is_minimum, is_maximum: !is_minimum }
}

#[test]
fn result_accessors() {
    let mut r = PointCurveExtrema::<Point2> { extrema: vec![], infinite_square_distance: None };
    assert_eq!(r.len(), 0);
    assert_eq!(r.min_index(), None);
    assert_eq!(r.max_index(), None);
    assert!(r.min_square_distance().is_infinite());
    assert_eq!(r.max_square_distance(), 0.);
    r.extrema = vec![extremum(0., 9., true), extremum(1., 4., true), extremum(2., 9., false)];
    assert_eq!(r.min_index(), Some(1));
    assert_eq!(r.max_index(), Some(0), "the first of equally far extrema");
}
