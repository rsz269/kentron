//! IntAna_IntQuadQuad_Test.cxx (including the IntAna_QuadQuadGeo and
//! IntAna2d_AnaIntersection cases).
//!
//! The out-of-range index checks are not ported: the curve links are an
//! `Option<(usize, bool)>` per curve.

use std::f64::consts::PI;

use kentron::geometry::curve::Line;
use kentron::geometry::curve2::Circle2;
use kentron::geometry::intersect::{self, QuadricCurves, QuadricIntersection};
use kentron::geometry::surface::{Cylinder, Sphere, Surface};
use kentron::math::{Axis2, Direction, Frame, Frame2};

use crate::{assert_near, p, p2};

fn sphere5() -> Surface {
    Surface::Sphere(Sphere::new(Frame::WORLD, 5.0))
}

fn cylinder3() -> Cylinder {
    Cylinder::new(Frame::WORLD, 3.0)
}

#[track_caller]
fn assert_links_valid(result: &QuadricCurves) {
    assert_eq!(result.next.len(), result.curves.len());
    for (i, next) in result.next.iter().enumerate() {
        if let Some((j, _)) = *next {
            assert!(j < result.curves.len());
            assert_ne!(j, i, "curve {i} continues into itself");
        }
    }
}

#[test]
fn cylinder_vs_sphere() {
    let result = intersect::cylinder_quadric(&cylinder3(), &sphere5(), 1e-7).unwrap();
    assert!(!result.same);
    assert!(!result.curves.is_empty());
    for curve in &result.curves {
        let (first, _) = curve.domain();
        assert!(curve.point(first).is_some());
    }
    assert_links_valid(&result);
}

#[test]
fn connected_curves() {
    let sphere = Surface::Sphere(Sphere::new(
        Frame::new(p(3.0, 0.0, 0.0), Direction::Z, Direction::X),
        2.0,
    ));
    let result =
        intersect::cylinder_quadric(&Cylinder::new(Frame::WORLD, 2.0), &sphere, 1e-7).unwrap();
    assert!(!result.same);
    assert_links_valid(&result);
}

fn z_cylinder(x: f64, radius: f64) -> Surface {
    Surface::Cylinder(Cylinder::new(Frame::new(p(x, 0.0, 0.0), Direction::Z, Direction::X), radius))
}

#[test]
fn cylinder_cylinder_near_tangent_collapses_to_single_line() {
    let result = intersect::quadrics(&z_cylinder(0.0, 1.0), &z_cylinder(2.0 - 1.0e-9, 1.0), 1.0e-4);
    let QuadricIntersection::Lines(lines) = result else { panic!("expected lines: {result:?}") };
    assert_eq!(lines.len(), 1);
}

#[test]
fn cylinder_cylinder_non_tangent_has_two_lines() {
    let result = intersect::quadrics(&z_cylinder(0.0, 1.0), &z_cylinder(1.8, 1.0), 1.0e-7);
    let QuadricIntersection::Lines(lines) = result else { panic!("expected lines: {result:?}") };
    assert_eq!(lines.len(), 2);
}

#[test]
fn cylinder_cylinder_skew_externally_tangent_has_point() {
    let (r1, r2) = (3.0, 2.0);
    let c1 = Cylinder::new(Frame::WORLD, r1);
    let c2 = Cylinder::new(Frame::new(p(0.0, r1 + r2, 0.0), Direction::X, Direction::Y), r2);
    let result = intersect::quadrics(&Surface::Cylinder(c1), &Surface::Cylinder(c2), 1.0e-7);
    let QuadricIntersection::Points(points) = result else { panic!("expected a point: {result:?}") };
    assert_eq!(points.len(), 1);
    assert_near(Line::from(c1.frame().z_axis()).distance_to_point(points[0]), r1, 1.0e-7);
    assert_near(Line::from(c2.frame().z_axis()).distance_to_point(points[0]), r2, 1.0e-7);
}

/// ModDataBug_23939: internally tangent circles meet at parameter 0 on both.
#[test]
fn tangent_circle_parameters() {
    let c1 = Circle2::new(Frame2::WORLD, 10.0);
    let c2 = Circle2::new(Frame2::from_axis(Axis2::new(p2(7.0, 0.0), Axis2::X.direction), true), 3.0);
    let points = intersect::circles2(&c1, &c2);
    assert_eq!(points.len(), 1);
    assert_near(points[0].t1, 0.0, 1.0e-15);
    assert_near(points[0].t2, 0.0, 1.0e-15);
    assert_near(points[0].point.x, 10.0, 1.0e-12);
    assert_near(points[0].point.y, 0.0, 1.0e-12);

    let reversed = intersect::circles2(&c2, &c1);
    assert_eq!(reversed.len(), 1);
    assert_near(reversed[0].t1, 0.0, 1.0e-15);
    assert_near(reversed[0].t2, 0.0, 1.0e-15);
}

/// ModDataBug_24375: nearly coincident circles still give their two crossings.
#[test]
fn nearly_coincident_circle_parameters() {
    let c1 = Circle2::new(Frame2::from_axis(Axis2::new(p2(10.0, 0.0), Axis2::X.direction), true), 100.0);
    let c2 = Circle2::new(
        Frame2::from_axis(Axis2::new(p2(10.00000000000001, 0.0), Axis2::X.direction), true),
        100.0,
    );
    let points = intersect::circles2(&c1, &c2);
    assert_eq!(points.len(), 2);
    for (point, expected) in points.iter().zip([PI / 2.0, 3.0 * PI / 2.0]) {
        assert_near(point.t1, expected, 1.0e-7);
        assert_near(point.t2, expected, 1.0e-7);
    }
}
