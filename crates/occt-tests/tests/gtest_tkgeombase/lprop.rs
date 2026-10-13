//! GeomLProp_CLProps2d_Test.cxx, GeomLProp_CurAndInf2d_Test.cxx.
//!
//! The stateful `SetCurve`/`SetParameter` API becomes a plain function call per point.

use std::f64::consts::PI;

use kentron::geometry::curve2::{Circle2, Curve2, Ellipse2, Hyperbola2, Line2, Parabola2};
use kentron::geometry::props::{
    curvature_extrema2, curvature_features2, curve_point2, inflections2, CurvatureFeatureKind,
};
use kentron::math::{Direction2, Frame2};
use kentron::tolerance::{LINEAR, PARAMETRIC};

use crate::{assert_near, p2};

fn circle() -> Curve2 {
    Curve2::Circle(Circle2::new(Frame2::WORLD, 5.0))
}

fn line() -> Curve2 {
    Curve2::Line(Line2::new(p2(0.0, 0.0), Direction2::X))
}

// ---- GeomLProp_CLProps2d ----

#[test]
fn circle_value() {
    let q = curve_point2(&circle(), 0.0, LINEAR).point;
    assert_near(q.x, 5.0, LINEAR);
    assert_near(q.y, 0.0, LINEAR);
    let q = curve_point2(&circle(), PI / 2.0, LINEAR).point;
    assert_near(q.x, 0.0, LINEAR);
    assert_near(q.y, 5.0, LINEAR);
}

#[test]
fn circle_tangent() {
    let tangent = curve_point2(&circle(), 0.0, LINEAR).tangent.expect("tangent");
    assert_near(tangent.x(), 0.0, LINEAR);
    assert_near(tangent.y(), 1.0, LINEAR);
}

#[test]
fn circle_curvature_constant() {
    for i in 0..8 {
        let t = i as f64 * PI / 4.0;
        assert_near(curve_point2(&circle(), t, LINEAR).curvature, 0.2, LINEAR);
    }
}

#[test]
fn circle_normal_perpendicular_to_tangent() {
    let props = curve_point2(&circle(), 0.0, LINEAR);
    let (n, t) = (props.normal.expect("normal"), props.tangent.expect("tangent"));
    assert_near(n.dot(t).abs(), 0.0, LINEAR);
}

#[test]
fn circle_center_of_curvature() {
    for t in [0.0, PI / 2.0] {
        let center = curve_point2(&circle(), t, LINEAR).center_of_curvature.expect("center");
        assert_near(center.x, 0.0, LINEAR);
        assert_near(center.y, 0.0, LINEAR);
    }
}

#[test]
fn line_tangent() {
    assert!(curve_point2(&line(), 0.0, LINEAR).tangent.is_some());
    let tangent = curve_point2(&line(), 5.0, LINEAR).tangent.unwrap();
    assert_near(tangent.x(), 1.0, LINEAR);
    assert_near(tangent.y(), 0.0, LINEAR);
}

#[test]
fn line_curvature_is_zero_and_center_undefined() {
    let props = curve_point2(&line(), 5.0, LINEAR);
    assert_near(props.curvature, 0.0, LINEAR);
    assert!(props.center_of_curvature.is_none());
}

#[test]
fn ellipse_curvature_at_vertices() {
    let ellipse = Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 10.0, 5.0));
    assert_near(curve_point2(&ellipse, 0.0, LINEAR).curvature, 10.0 / 25.0, 1e-6);
    assert_near(curve_point2(&ellipse, PI / 2.0, LINEAR).curvature, 5.0 / 100.0, 1e-6);
}

#[test]
fn circle_derivatives() {
    let props = curve_point2(&circle(), 0.0, LINEAR);
    assert_near(props.d1.x, 0.0, LINEAR);
    assert_near(props.d1.y, 5.0, LINEAR);
    assert_near(props.d2.x, -5.0, LINEAR);
    assert_near(props.d2.y, 0.0, LINEAR);
}

// ---- GeomLProp_CurAndInf2d ----

fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 10.0, 3.0))
}

fn hyperbola() -> Curve2 {
    Curve2::Hyperbola(Hyperbola2::new(Frame2::WORLD, 6.0, 3.0))
}

#[test]
fn circle_has_no_features() {
    assert!(curvature_features2(&circle()).is_empty());
    assert!(inflections2(&circle()).is_empty());
    assert!(curvature_extrema2(&circle()).is_empty());
}

#[test]
fn ellipse_curvature_extrema() {
    let extrema = curvature_extrema2(&ellipse());
    assert_eq!(extrema.len(), 4);
    let count = |kind| extrema.iter().filter(|f| f.kind == kind).count();
    assert_eq!(count(CurvatureFeatureKind::MinRadius), 2);
    assert_eq!(count(CurvatureFeatureKind::MaxRadius), 2);
    for (feature, expected) in extrema.iter().zip([0.0, PI / 2.0, PI, 3.0 * PI / 2.0]) {
        assert_near(feature.parameter, expected, 1e-6);
    }
    assert!(extrema.windows(2).all(|w| w[0].parameter < w[1].parameter));
}

#[test]
fn ellipse_has_no_inflections() {
    assert!(inflections2(&ellipse()).is_empty());
    assert_eq!(curvature_features2(&ellipse()).len(), 4);
}

#[test]
fn hyperbola_vertex_only() {
    let extrema = curvature_extrema2(&hyperbola());
    assert_eq!(extrema.len(), 1);
    assert_near(extrema[0].parameter, 0.0, PARAMETRIC);
    assert_eq!(extrema[0].kind, CurvatureFeatureKind::MinRadius);
    assert!(inflections2(&hyperbola()).is_empty());
}

#[test]
fn parabola_vertex_only() {
    let parabola = Curve2::Parabola(Parabola2::new(Frame2::WORLD, 2.0));
    let extrema = curvature_extrema2(&parabola);
    assert_eq!(extrema.len(), 1);
    assert_near(extrema[0].parameter, 0.0, PARAMETRIC);
    assert_eq!(extrema[0].kind, CurvatureFeatureKind::MinRadius);
}
