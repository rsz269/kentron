//! GeomProjLib_Test.cxx, ProjLib_Cone_Test.cxx, ProjLib_ComputeApproxOnPolarSurface_Test.cxx.

use std::f64::consts::PI;

use kentron::geometry::curve::{Circle, Curve, Hyperbola, Parabola, TrimmedCurve};
use kentron::geometry::curve2::Curve2;
use kentron::geometry::surface::{Cone, Cylinder, Plane, Surface};
use kentron::geometry::{measure, project, ParametricCurve, ParametricSurface};
use kentron::math::{Direction, Frame};
use kentron::tolerance::ANGULAR;

use crate::{assert_near, p};

// ---- GeomProjLib ----

fn occ31661_frame() -> Frame {
    Frame::new(p(0.0, 0.0, 0.0), Direction::new(1.0, 1.0, 1.0), Direction::new(2.0, 0.0, -2.0))
}

/// OCC31661_1.
#[test]
fn project_parabola_on_plane() {
    let parabola = Curve::Parabola(Parabola::new(occ31661_frame(), 10.0));
    let projected = project::curve_on_plane(&parabola, &Plane::XY, Direction::Z, false).unwrap();
    assert!(matches!(projected, Curve::Parabola(_)), "{projected:?}");

    let trimmed = Curve::Trimmed(TrimmedCurve::new(parabola, -100.0, 100.0));
    let projected = project::curve_on_plane(&trimmed, &Plane::XY, Direction::Z, false).unwrap();
    assert!(matches!(projected, Curve::Trimmed(_)), "{projected:?}");
    assert_near(measure::length(&projected), 408.40363195229503, 2.0e-5);
    assert_near(projected.first_parameter(), -91.077748943768597, 2.0e-6);
    assert_near(projected.last_parameter(), 72.221567418462357, 2.0e-6);
}

/// OCC31661_2.
#[test]
fn project_hyperbola_on_plane() {
    let hyperbola = Curve::Hyperbola(Hyperbola::new(occ31661_frame(), 10.0, 10.0));
    let projected = project::curve_on_plane(&hyperbola, &Plane::XY, Direction::Z, false).unwrap();
    assert!(matches!(projected, Curve::Hyperbola(_)), "{projected:?}");

    let trimmed = Curve::Trimmed(TrimmedCurve::new(hyperbola, -5.0, 5.0));
    let projected = project::curve_on_plane(&trimmed, &Plane::XY, Direction::Z, false).unwrap();
    assert!(matches!(projected, Curve::Trimmed(_)), "{projected:?}");
    assert_near(measure::length(&projected), 1664.3732976598988, 2.0e-5);
    assert_near(projected.first_parameter(), -5.23179933356147, 1.0e-7);
    assert_near(projected.last_parameter(), 4.76820064934972, 1.0e-7);
}

// ---- ProjLib_Cone ----

fn cone() -> Surface {
    Surface::Cone(Cone::new(Frame::WORLD, PI / 6.0, 5.0))
}

/// The cone's radius at height `z`: OCCT's tests use radius 2, off the cone, which
/// `analytic_pcurve` (for curves on the surface) rejects; see GTEST_PORTING.md.
fn cone_radius_at(z: f64) -> f64 {
    5.0 + z * (PI / 6.0).tan()
}

#[test]
fn cone_circle_parallel_axes_gives_line() {
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 0.0, 10.0), Direction::Z, Direction::X),
        cone_radius_at(10.0),
    ));
    let Some(Curve2::Line(line)) = project::analytic_pcurve(&circle, &cone()) else {
        panic!("expected a line");
    };
    assert_near(line.direction.y(), 0.0, ANGULAR);
}

#[test]
fn cone_circle_non_parallel_axes_not_analytic() {
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 0.0, 10.0), Direction::X, Direction::Y),
        2.0,
    ));
    assert!(project::analytic_pcurve(&circle, &cone()).is_none());
}

#[test]
fn cone_circle_opposite_normal_gives_reversed_line() {
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 0.0, 6.0), -Direction::Z, Direction::X),
        cone_radius_at(6.0),
    ));
    let Some(Curve2::Line(line)) = project::analytic_pcurve(&circle, &cone()) else {
        panic!("expected a line");
    };
    assert!(line.direction.x() < 0.0);
}

// ---- ProjLib_ComputeApproxOnPolarSurface ----

#[test]
fn project_circle_on_cylinder() {
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 0.0, 2.0), Direction::Z, Direction::X),
        5.0,
    ));
    let surface = Surface::Cylinder(Cylinder::new(Frame::WORLD, 5.0));
    let pcurve = project::pcurve(&circle, &surface, 1.0e-6).unwrap();

    let first = circle.first_parameter().max(pcurve.first_parameter());
    let last = circle.last_parameter().min(pcurve.last_parameter());
    assert!(first < last);
    for t in [first, 0.5 * (first + last), last] {
        let uv = pcurve.point(t);
        assert!(surface.point(uv.x, uv.y).distance(circle.point(t)) <= 1.0e-3, "t = {t}");
    }
}

