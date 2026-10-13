//! Geom_SurfaceEval_Test.cxx: surface evaluation.

use std::f64::consts::{FRAC_PI_3, FRAC_PI_4, FRAC_PI_6};

use kentron::geometry::curve::{Circle, Curve, Line};
use kentron::geometry::surface::{
    BSplineSurface, BezierSurface, Cylinder, ExtrusionSurface, KnotVector, OffsetSurface, Plane,
    RevolutionSurface, Sphere, Surface,
};
use kentron::geometry::ParametricSurface;
use kentron::math::{Axis, Direction, Frame};
use kentron::tolerance;

use crate::{assert_point_near, p};

const TOL: f64 = tolerance::LINEAR;

fn xy_plane() -> Plane {
    Plane::new(p(0.0, 0.0, 0.0), Direction::Z)
}

#[test]
fn plane_point_and_partials() {
    assert_point_near(xy_plane().point(1.0, 2.0), p(1.0, 2.0, 0.0), TOL);
    let (_, du, dv) = xy_plane().d1(1.0, 2.0);
    assert!((du.length() - 1.0).abs() <= TOL && (dv.length() - 1.0).abs() <= TOL);
    assert!(du.dot(dv).abs() <= TOL);
}

#[test]
fn sphere_partials() {
    let (q, du, dv) = Sphere::new(Frame::WORLD, 5.0).d1(FRAC_PI_4, FRAC_PI_6);
    assert!((q.distance(p(0.0, 0.0, 0.0)) - 5.0).abs() <= TOL);
    assert!(du.length() > TOL && dv.length() > TOL);
}

fn assert_d2_consistent(s: &impl ParametricSurface, u: f64, v: f64) {
    let d2 = s.d2(u, v);
    let (q, du, dv) = s.d1(u, v);
    assert_point_near(d2.point, q, TOL);
    assert!((d2.du - du).length() <= TOL && (d2.dv - dv).length() <= TOL);
    assert!((d2.duu - s.dn(u, v, 2, 0)).length() <= TOL);
    assert!((d2.dvv - s.dn(u, v, 0, 2)).length() <= TOL);
    assert!((d2.duv - s.dn(u, v, 1, 1)).length() <= TOL);
}

#[test]
fn bspline_d2_consistent() {
    let poles: Vec<Vec<_>> = (0..3)
        .map(|i| (0..3).map(|j| p(i as f64, j as f64, 0.1 * ((i as f64 - 1.0) * (j as f64 - 1.0)))).collect())
        .collect();
    let knots = || KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2);
    assert_d2_consistent(&BSplineSurface::new(poles, knots(), knots()), 0.5, 0.5);
}

#[test]
fn bezier_d2_consistent() {
    let poles = vec![
        vec![p(0.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 2.0, 0.0)],
        vec![p(1.0, 0.0, 0.5), p(1.0, 1.0, 1.0), p(1.0, 2.0, 0.5)],
        vec![p(2.0, 0.0, 0.0), p(2.0, 1.0, 0.0), p(2.0, 2.0, 0.0)],
    ];
    assert_d2_consistent(&BezierSurface::new(poles), 0.3, 0.7);
}

#[test]
fn offset_plane_point() {
    let offset = OffsetSurface::new(Surface::Plane(xy_plane()), 3.0);
    assert_point_near(offset.point(1.0, 2.0), p(1.0, 2.0, 3.0), TOL);
}

#[test]
fn cylinder_partials() {
    let (q, du, dv) = Cylinder::new(Frame::WORLD, 4.0).d1(FRAC_PI_4, 3.0);
    assert!(((q.x * q.x + q.y * q.y).sqrt() - 4.0).abs() <= TOL);
    assert!((q.z - 3.0).abs() <= TOL);
    assert!(du.length() > TOL);
    assert!(dv.x.abs() <= TOL && dv.y.abs() <= TOL && (dv.z - 1.0).abs() <= TOL);
}

#[test]
fn extrusion_point() {
    let s = ExtrusionSurface::new(Curve::Circle(Circle::new(Frame::WORLD, 3.0)), Direction::Z);
    let q = s.point(FRAC_PI_3, 5.0);
    assert!(((q.x * q.x + q.y * q.y).sqrt() - 3.0).abs() <= TOL);
    assert!((q.z - 5.0).abs() <= TOL);
}

#[test]
fn revolution_partials() {
    let s = RevolutionSurface::new(Curve::Line(Line::new(p(1.0, 0.0, 0.0), Direction::Z)), Axis::Z);
    let (q, du, dv) = s.d1(FRAC_PI_4, 2.0);
    assert!(((q.x * q.x + q.y * q.y).sqrt() - 1.0).abs() <= TOL);
    assert!((q.z - 2.0).abs() <= TOL);
    assert!(du.length() > TOL);
    assert!(dv.x.abs() <= TOL && dv.y.abs() <= TOL && (dv.z - 1.0).abs() <= TOL);
}
