//! GCPnts_AbscissaPoint_Test.cxx, GCPnts_QuasiUniformDeflection_Test.cxx.

use std::f64::consts::PI;

use kentron::geometry::curve::{BSplineCurve, Circle, Curve, Line, TrimmedCurve};
use kentron::geometry::measure;
use kentron::geometry::sample::quasi_uniform_deflection;
use kentron::math::{Direction, Frame};
use kentron::tolerance::LINEAR;

use crate::{assert_near, p};

fn x_line() -> Curve {
    Curve::Line(Line::new(p(0.0, 0.0, 0.0), Direction::X))
}

fn circle(radius: f64) -> Curve {
    Curve::Circle(Circle::new(Frame::WORLD, radius))
}

#[test]
fn line_length() {
    assert_near(measure::length_between(&x_line(), 0.0, 10.0), 10.0, LINEAR);
}

#[test]
fn circle_arc_length() {
    assert_near(measure::length_between(&circle(1.0), 0.0, PI), PI, LINEAR);
}

#[test]
fn full_circle_length() {
    assert_near(measure::length_between(&circle(5.0), 0.0, 2.0 * PI), 10.0 * PI, LINEAR);
}

#[test]
fn parameter_at_length_line() {
    let segment = Curve::Trimmed(TrimmedCurve::new(x_line(), 0.0, 10.0));
    assert_near(measure::parameter_at_length(&segment, 0.0, 5.0), 5.0, LINEAR);
}

#[test]
fn parameter_at_length_circle() {
    assert_near(measure::parameter_at_length(&circle(1.0), 0.0, PI / 2.0), PI / 2.0, LINEAR);
}

/// ModDataBug_33828: stable point counts for decreasing deflections.
#[test]
fn quasi_uniform_deflection_stable_point_counts() {
    let poles = [
        (0.163, 0.233),
        (0.158, 0.204),
        (0.139, 0.180),
        (0.086, 0.159),
        (0.055, 0.163),
        (0.009, 0.196),
        (-0.004, 0.225),
        (0.002, 0.281),
        (0.019, 0.307),
        (0.070, 0.332),
        (0.101, 0.331),
        (0.149, 0.301),
        (0.164, 0.274),
        (0.163, 0.246),
    ]
    .map(|(x, y)| p(x, y, 0.0));
    let curve = Curve::BSpline(BSplineCurve::new(
        poles.to_vec(),
        vec![0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0],
        vec![4, 2, 2, 2, 2, 2, 4],
        3,
    ));
    for (deflection, expected) in [
        (0.5, 2),
        (0.1, 3),
        (0.05, 5),
        (0.025, 5),
        (0.007, 9),
        (0.005, 17),
        (0.0005, 33),
        (0.0003, 65),
        (0.0002, 65),
        (0.0001, 73),
    ] {
        let sampling = quasi_uniform_deflection(&curve, deflection).unwrap();
        assert_eq!(sampling.points.len(), expected, "deflection {deflection}");
        assert!(sampling.deflection <= deflection, "deflection {deflection}");
    }
}
