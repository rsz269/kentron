//! Geom2d_OffsetCurve_Test.cxx and Adaptor2d_OffsetCurve_Test.cxx.

use kentron::geometry::convert::{self, CurveApproximationOptions};
use kentron::geometry::curve2::{Circle2, Curve2, Line2, OffsetCurve2};
use kentron::geometry::{measure, Continuity, ParametricCurve};
use kentron::math::{Direction2, Frame2};
use kentron::tolerance;

use crate::{assert_point2_near, p2};

const TOL: f64 = tolerance::LINEAR;

fn line() -> Curve2 {
    Curve2::Line(Line2::new(p2(0.0, 0.0), Direction2::X))
}

#[test]
fn offset_value_and_setter() {
    let mut c = OffsetCurve2::new(Curve2::Circle(Circle2::new(Frame2::WORLD, 5.0)), 2.0);
    assert_eq!(c.offset(), 2.0);
    c.set_offset(10.0);
    assert_eq!(c.offset(), 10.0);
}

#[test]
fn offset_circle_properties() {
    let c = OffsetCurve2::new(Curve2::Circle(Circle2::new(Frame2::WORLD, 5.0)), 2.0);
    assert!(c.is_periodic() && c.is_closed());
    // Offsetting to the right of a counter-clockwise circle grows it (OCCT's
    // Geom2d_OffsetCurve: P + offset·(T.y, -T.x)/|T|).
    assert!((c.point(0.0).distance(p2(0.0, 0.0)) - 7.0).abs() <= TOL);
}

#[test]
fn offset_line_points() {
    // Right of a line along +X is -Y.
    assert_point2_near(OffsetCurve2::new(line(), 2.0).point(3.0), p2(3.0, -2.0), TOL);
    assert_point2_near(OffsetCurve2::new(line(), -2.0).point(3.0), p2(3.0, 2.0), TOL);
    assert_point2_near(OffsetCurve2::new(line(), 0.0).point(3.0), p2(3.0, 0.0), TOL);
    let (_, d1) = OffsetCurve2::new(line(), 2.0).d1(1.0);
    assert!((d1.x - 1.0).abs() <= TOL && d1.y.abs() <= TOL);
}

#[test]
fn offset_line_domain() {
    let c = OffsetCurve2::new(line(), 1.0);
    assert!(!c.is_closed() && !c.is_periodic());
    assert_eq!(c.basis(), &line());
    assert!((c.resolution(0.01) - 0.01).abs() <= 1.0e-12);
}

/// fclasses/bug31381: approximating an offset of a rational B-spline circle keeps its
/// length.
#[test]
fn rational_bspline_offset_approximation_keeps_length() {
    let circle = Curve2::Circle(Circle2::new(Frame2::WORLD, 1.0));
    let bspline = convert::curve2_to_bspline(&circle).unwrap();
    let offset = Curve2::Offset(OffsetCurve2::new(Curve2::BSpline(bspline), -0.5));
    let options = CurveApproximationOptions {
        tolerance: 1.0e-4,
        continuity: Continuity::C2,
        max_segments: 25,
        max_degree: 5,
    };
    let approximation = convert::approximate_curve2(&offset, &options).unwrap();
    let length = measure::length2(&Curve2::BSpline(approximation.curve));
    // The offset is a circle of radius 0.5. OCCT expects its own approximation's length,
    // 3.1415850499856011, within 1e-6; see GTEST_PORTING.md.
    assert!((length - std::f64::consts::PI).abs() <= 1.0e-5, "{length}");
}
