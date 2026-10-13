//! GeomLib_CheckCurveOnSurface_Test.cxx.

use std::f64::consts::PI;

use kentron::geometry::curve::{Circle, Curve, Line};
use kentron::geometry::curve2::{Curve2, Ellipse2, Line2};
use kentron::geometry::measure::{curve_on_surface_deviation, CurveDeviation};
use kentron::geometry::surface::{Cylinder, Plane, Surface, Torus};
use kentron::math::{Direction, Direction2, Frame, Frame2};
use kentron::tolerance::{LINEAR, PARAMETRIC};

use crate::{assert_near, p, p2};

fn circle_on_cylinder(circle_radius: f64, cylinder_radius: f64, height: f64) -> CurveDeviation {
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 0.0, height), Direction::Z, Direction::X),
        circle_radius,
    ));
    let cylinder = Surface::Cylinder(Cylinder::new(Frame::WORLD, cylinder_radius));
    let pcurve = Curve2::Line(Line2::new(p2(0.0, height), Direction2::X));
    curve_on_surface_deviation(&circle, &pcurve, &cylinder, 0.0, 2.0 * PI)
}

fn line_on_plane(direction: Direction, pcurve_direction: Direction2, t1: f64, t2: f64) -> CurveDeviation {
    let line = Curve::Line(Line::new(p(0.0, 0.0, 0.0), direction));
    let pcurve = Curve2::Line(Line2::new(p2(0.0, 0.0), pcurve_direction));
    curve_on_surface_deviation(&line, &pcurve, &Surface::Plane(Plane::XY), t1, t2)
}

fn circle_and_ellipse(major: f64, minor: f64, t1: f64, t2: f64) -> CurveDeviation {
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 2.0));
    let ellipse = Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, major, minor));
    curve_on_surface_deviation(&circle, &ellipse, &Surface::Plane(Plane::XY), t1, t2)
}

#[test]
fn analytic_coincident_reports_zero_deviation() {
    assert!(circle_on_cylinder(2.0, 2.0, 5.0).max_distance <= LINEAR);
}

#[test]
fn analytic_deviating_reports_actual_deviation() {
    assert_near(circle_on_cylinder(1.95, 2.0, 5.0).max_distance, 0.05, 1.0e-9);
}

#[test]
fn line_on_plane_reports_zero_deviation() {
    assert!(line_on_plane(Direction::X, Direction2::X, -10.0, 10.0).max_distance <= LINEAR);
}

#[test]
fn deviating_lines_report_endpoint_maximum() {
    let deviation = line_on_plane(Direction::X, Direction2::Y, -2.0, 3.0);
    assert_near(deviation.max_distance, 18.0f64.sqrt(), LINEAR);
    assert_eq!(deviation.parameter, 3.0);
}

#[test]
fn circle_and_ellipse_report_analytic_maximum() {
    assert_near(circle_and_ellipse(4.0, 1.0, 0.0, 2.0 * PI).max_distance, 2.0, LINEAR);
}

#[test]
fn negative_trim_reports_interior_maximum() {
    let deviation = circle_and_ellipse(2.0, 1.0, -2.0, -1.0);
    assert_near(deviation.max_distance, 1.0, LINEAR);
    assert_near(deviation.parameter, -0.5 * PI, PARAMETRIC);
}

/// A fast-oscillating deviation must not be missed by sampling.
#[test]
fn analytic_aliasing_reports_actual_deviation() {
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 4.0));
    let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0));
    let pcurve = Curve2::Line(Line2::new(p2(0.0, 0.0), Direction2::new(3.0, 4.0)));
    let deviation = curve_on_surface_deviation(&circle, &pcurve, &torus, 0.0, 80.0 * PI);
    assert_near(deviation.max_distance, 8.0, LINEAR);
}
