//! Lines against surfaces, closed forms and sampled ones, against the exact crossings.

use std::f64::consts::{PI, TAU};

use super::*;
use crate::geometry::bspline::{BSplineSurface, KnotVector};
use crate::geometry::surface::{Cone, Cylinder, Sphere, Torus, TrimmedSurface};
use crate::math::{Direction, Frame};

fn ws(face: &LineFace, origin: Point, d: Vector) -> Vec<f64> {
    let (points, lo, hi) = face.perform(origin, d, -INFINITE, INFINITE).expect("ready");
    assert!(lo <= -1e100 && hi >= 1e100);
    for p in &points {
        assert!(p.point.distance(origin + d * p.w) <= 1e-12 * (1.0 + p.w.abs()));
    }
    let mut w: Vec<f64> = points.iter().map(|p| p.w).collect();
    w.sort_by(f64::total_cmp);
    w
}

fn assert_ws(got: &[f64], want: &[f64], tol: f64) {
    assert_eq!(got.len(), want.len(), "{got:?} vs {want:?}");
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() <= tol, "{got:?} vs {want:?}");
    }
}

#[test]
fn line_meets_a_tilted_cylinder() {
    let frame = Frame::new(Point::new(1.0, 0.0, 0.0), Direction::new(0.0, 1.0, 1.0), Direction::new(1.0, 0.0, 0.0));
    let face = LineFace::new(Surface::Cylinder(Cylinder::new(frame, 2.0)), (0.0, TAU, -10.0, 10.0));
    assert_ws(&ws(&face, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0)), &[-1.0, 3.0], 1e-12);
    // A line along the axis: none.
    assert_ws(&ws(&face, Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)), &[], 0.0);
    // A generator: none either.
    assert_ws(&ws(&face, Point::new(3.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)), &[], 0.0);
}

#[test]
fn line_meets_a_torus_four_times() {
    let face = LineFace::new(Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0)), (0.0, TAU, 0.0, TAU));
    assert_ws(&ws(&face, Point::new(-10.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)), &[6.0, 8.0, 12.0, 14.0], 1e-9);
    // Along the axis: none.
    assert_ws(&ws(&face, Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0)), &[], 0.0);
}

#[test]
fn sphere_through_its_centre() {
    let centre = Point::new(3.0, -2.0, 7.0);
    let face = LineFace::new(Surface::Sphere(Sphere::new(Frame::new(centre, Direction::Z, Direction::X), 1.5)), (0.0, TAU, -PI / 2.0, PI / 2.0));
    for k in 0..20 {
        let (a, b) = (0.7 * k as f64, 0.37 * k as f64 - 1.0);
        let d = Vector::new(a.cos() * b.cos(), a.sin() * b.cos(), b.sin());
        assert_ws(&ws(&face, centre, d), &[-1.5, 1.5], 1e-12);
    }
}

#[test]
fn tangent_and_periodic_cases() {
    // A line touching a cylinder: one point (once), tangent.
    let face = LineFace::new(Surface::Cylinder(Cylinder::new(Frame::WORLD, 2.0)), (0.0, TAU, -1.0, 1.0));
    let (points, _, _) = face.perform(Point::new(-5.0, 2.0, 0.0), Vector::new(1.0, 0.0, 0.0), 0.0, INFINITE).unwrap();
    assert_eq!(points.len(), 1, "{points:?}");
    assert!((points[0].w - 5.0).abs() <= 1e-6);
    assert_eq!(points[0].transition, Transition::Tangent);
    // A cylinder face over [π, 3π]: u comes out in it.
    let face = LineFace::new(Surface::Cylinder(Cylinder::new(Frame::WORLD, 2.0)), (PI, 3.0 * PI, -1.0, 1.0));
    let (points, _, _) = face.perform(Point::new(-5.0, 0.5, 0.0), Vector::new(1.0, 0.0, 0.0), 0.0, INFINITE).unwrap();
    assert_eq!(points.len(), 2);
    assert!(points.iter().all(|p| p.u >= PI && p.u < 3.0 * PI));
    // In and out of the cylinder: against the normal (outward) first.
    let mut sorted = points.clone();
    sorted.sort_by(|a, b| a.w.total_cmp(&b.w));
    assert_eq!((sorted[0].transition, sorted[1].transition), (Transition::In, Transition::Out));
    // A cone face on one nappe: one point where the line crosses both.
    let cone = Surface::Cone(Cone::new(Frame::WORLD, 0.5, 0.0));
    let face = LineFace::new(cone, (0.0, TAU, 0.0, 10.0));
    let (points, _, _) = face.perform(Point::new(0.5, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0), -INFINITE, INFINITE).unwrap();
    assert_eq!(points.len(), 1, "{points:?}");
    assert!((points[0].point.z - 0.5 / 0.5f64.tan()).abs() <= 1e-12);
}

#[test]
fn torus_touched_on_its_inner_equator() {
    let face = LineFace::new(Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0)), (0.0, TAU, 0.0, TAU));
    let (points, _, _) = face.perform(Point::new(-10.0, 2.0, 0.0), Vector::new(1.0, 0.0, 0.0), -INFINITE, INFINITE).unwrap();
    assert!(!points.is_empty());
    for p in &points {
        // (Each found point on the torus.)
        let (u, v) = (p.u, p.v);
        assert!(Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0)).point(u, v).distance(p.point) <= 1e-6);
    }
    // (The outer crossings, at x = ±√(16 − 4) = ±√12.)
    assert!(points.iter().any(|p| (p.w - (10.0 - 12f64.sqrt())).abs() <= 1e-9));
    assert!(points.iter().any(|p| (p.w - (10.0 + 12f64.sqrt())).abs() <= 1e-9));
}

/// A wave `z = 0.01 sin(20 x)` as a fine B-spline over the unit square.
fn wave() -> Surface {
    let n = 160;
    let xs: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let points: Vec<Point> = xs.iter().map(|&x| Point::new(x, 0.0, 0.01 * (20.0 * x).sin())).collect();
    let row = crate::geometry::fit::interpolate_with_parameters(&points, &xs).expect("a row");
    let poles: Vec<Vec<Point>> = row.poles().iter().map(|p| vec![*p, Point::new(p.x, 1.0, p.z)]).collect();
    let u = KnotVector::new(row.knots().to_vec(), row.multiplicities().to_vec(), row.degree());
    let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
    Surface::BSpline(BSplineSurface::new(poles, u, v))
}

#[test]
fn line_along_a_wave() {
    let s = wave();
    let face = LineFace::new(s.clone(), (0.0, 1.0, 0.0, 1.0));
    let (points, _, _) = face.perform(Point::new(-1.0, 0.5, 0.005), Vector::new(1.0, 0.0, 0.0), -INFINITE, INFINITE).unwrap();
    let mut w: Vec<f64> = points.iter().map(|p| p.w - 1.0).collect();
    w.sort_by(f64::total_cmp);
    // sin(20 x) = 0.5 at x = (π/6 + 2kπ)/20 and (5π/6 + 2kπ)/20 within [0, 1].
    let mut want = Vec::new();
    for k in 0..4 {
        for a in [PI / 6.0, 5.0 * PI / 6.0] {
            let x = (a + TAU * k as f64) / 20.0;
            if x <= 1.0 {
                want.push(x);
            }
        }
    }
    want.sort_by(f64::total_cmp);
    // (The interpolated wave differs from the sine by little; its own crossings are checked
    // on it.)
    assert_eq!(w.len(), want.len(), "{w:?} vs {want:?}");
    for (g, x) in w.iter().zip(&want) {
        assert!((g - x).abs() <= 1e-4, "{w:?} vs {want:?}");
    }
    for p in &points {
        assert!(s.point(p.u, p.v).distance(p.point) <= 1e-12);
    }
}

/// A bilinear square patch.
fn square() -> Surface {
    let poles = vec![vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)], vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)]];
    let kv = || KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
    Surface::BSpline(BSplineSurface::new(poles, kv(), kv()))
}

#[test]
fn crossing_at_a_patch_edge() {
    let face = LineFace::new(square(), (0.0, 1.0, 0.0, 1.0));
    for x in [1.0 + 5e-10, 1.0, 1.0 - 1e-12, 0.0, -5e-10] {
        let (points, _, _) = face.perform(Point::new(x, 0.3, 1.0), Vector::new(0.0, 0.0, -1.0), 0.0, INFINITE).unwrap();
        assert_eq!(points.len(), 1, "{x}: {points:?}");
        assert!((points[0].w - 1.0).abs() <= 1e-12);
    }
    let (points, _, _) = face.perform(Point::new(1.0 + 1e-6, 0.3, 1.0), Vector::new(0.0, 0.0, -1.0), 0.0, INFINITE).unwrap();
    assert!(points.is_empty());
}

#[test]
fn not_ready_without_extent() {
    let face = LineFace::new(square(), (0.0, 1.0, 0.5, 0.5));
    assert!(face.perform(Point::new(0.5, 0.5, 1.0), Vector::new(0.0, 0.0, -1.0), 0.0, INFINITE).is_none());
}

#[test]
fn sampled_sphere_and_trim() {
    // A sphere as a B-spline: a line through its middle, off the axis.
    let sphere = Surface::Trimmed(TrimmedSurface::new(Surface::Sphere(Sphere::new(Frame::WORLD, 2.0)), 0.0, TAU, -PI / 2.0, PI / 2.0));
    let s = Surface::BSpline(crate::geometry::convert::surface_to_bspline(&sphere).expect("a B-spline sphere"));
    let (u0, u1, v0, v1) = s.bounds();
    let face = LineFace::new(s, (u0, u1, v0, v1));
    let (y, z) = (0.3, 0.4);
    let h = (4.0 - y * y - z * z) as f64;
    assert_ws(&ws(&face, Point::new(-5.0, y, z), Vector::new(1.0, 0.0, 0.0)), &[5.0 - h.sqrt(), 5.0 + h.sqrt()], 1e-9);
    // Through the poles, and through the seam: each crossing once.
    assert_ws(&ws(&face, Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0)), &[3.0, 7.0], 1e-9);
    assert_ws(&ws(&face, Point::new(-5.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)), &[3.0, 7.0], 1e-9);
}
