//! Point inversion, singular places, closure and curve feet against closed forms.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use super::*;
use crate::geometry::bspline::{BSplineCurve, BSplineSurface, KnotVector};
use crate::geometry::curve::Circle;
use crate::geometry::surface::{Cone, Cylinder, Sphere, Torus, TrimmedSurface};
use crate::math::{Direction, Frame};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn period_shift_values() {
    assert!(close(period_shift(7.0, 0.5, TAU), -TAU, 1e-15));
    assert!(close(period_shift(-3.5, 3.0, TAU), TAU, 1e-15));
    assert_eq!(period_shift(1.0, 1.0 + PI, TAU), 0.0);
    assert_eq!(period_shift(1.0, 1.2, TAU), 0.0);
    assert!(close(period_shift(1.0, 3.5, 0.0), 2.5, 1e-15));
    assert!(close(period_shift(20.0, 0.0, -TAU), -3.0 * TAU, 1e-12));
}

fn cylinder(r: f64) -> Surface {
    Surface::Cylinder(Cylinder::new(Frame::WORLD, r))
}

#[test]
fn nearest_on_elementary_surfaces() {
    // A cylinder face over [π, 3π]: an angle of 0.5 comes out a turn on.
    let s = Surface::Trimmed(TrimmedSurface::new(cylinder(2.0), PI, 3.0 * PI, -5.0, 5.0));
    let mut inv = SurfaceInverse::new(&s);
    let p = Point::new(3.0 * 0.5f64.cos(), 3.0 * 0.5f64.sin(), 1.25);
    let (uv, gap) = inv.nearest(p, 1e-7);
    assert!(close(uv.x, 0.5 + TAU, 1e-12) && close(uv.y, 1.25, 1e-12), "{uv:?}");
    assert!(close(gap, 1.0, 1e-12));
    // A sphere: the latitude asin(z/|p|).
    let s = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0));
    let mut inv = SurfaceInverse::new(&s);
    let p = Point::new(1.0, -2.0, 1.5);
    let (uv, gap) = inv.nearest(p, 1e-7);
    let len = (p - Point::ORIGIN).length();
    assert!(close(uv.y, (1.5 / len).asin(), 1e-12));
    assert!(close(uv.x, (-2.0f64).atan2(1.0) + TAU, 1e-12), "{uv:?}");
    assert!(close(gap, len - 2.0, 1e-12));
    // A torus: both angles.
    let s = Surface::Torus(Torus::new(Frame::WORLD, 5.0, 1.0));
    let mut inv = SurfaceInverse::new(&s);
    let (u, v) = (2.0, 4.0);
    let q = s.point(u, v);
    let (_, du, dv) = s.d1(u, v);
    let p = q + du.cross(dv).normalized() * 0.25;
    let (uv, gap) = inv.nearest(p, 1e-7);
    assert!(close(uv.x, u, 1e-12) && close(uv.y, v, 1e-12), "{uv:?}");
    assert!(close(gap, 0.25, 1e-12));
}

/// A bicubic patch with a fold: a wave along u.
fn wavy() -> Surface {
    let poles: Vec<Vec<Point>> = (0..6)
        .map(|i| (0..5).map(|j| Point::new(i as f64, j as f64, if i % 2 == 0 { 0.0 } else { 1.5 } + 0.1 * (j as f64 - 2.0).powi(2))).collect())
        .collect();
    let u = KnotVector::new(vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], 3);
    let v = KnotVector::new(vec![0.0, 1.0, 2.0], vec![4, 1, 4], 3);
    Surface::BSpline(BSplineSurface::new(poles, u, v))
}

#[test]
fn nearest_on_a_bspline_patch() {
    let s = wavy();
    let mut inv = SurfaceInverse::new(&s);
    for &(u, v, h) in &[(0.3, 0.4, 0.01), (1.7, 1.2, -0.02), (2.5, 0.9, 0.005), (1.1, 1.9, 1e-4)] {
        let (q, du, dv) = s.d1(u, v);
        let p = q + du.cross(dv).normalized() * h;
        let (uv, gap) = inv.nearest(p, 1e-7);
        assert!(close(uv.x, u, 1e-9) && close(uv.y, v, 1e-9), "{u} {v}: {uv:?}");
        assert!(close(gap, h.abs(), 1e-12 * (1.0 + h.abs()) + 1e-15), "{gap} {h}");
        // The local search from a nearby hint finds the same foot.
        let (uv2, gap2) = inv.nearest_from(Point2::new(u + 0.05, v - 0.05), p, 1e-7, 0.0);
        assert!(close(uv2.x, u, 1e-9) && close(uv2.y, v, 1e-9), "{uv2:?}");
        assert!(close(gap2, gap, 1e-12));
    }
    // Beyond the bound u = 0: the foot on the boundary.
    let (q, du, _) = s.d1(0.0, 1.0);
    let p = q - du.normalized() * 0.3;
    let (uv, _) = inv.nearest(p, 1e-7);
    assert!(close(uv.x, 0.0, 1e-12), "{uv:?}");
    let brute = brute_nearest(&s, p);
    assert!(close(s.point(uv.x, uv.y).distance(p), brute, 1e-7));
    // Above a crest of the wave, from a hint in the trough: the global search still finds the
    // nearest point.
    for p in [Point::new(1.0, 2.0, 2.2), Point::new(3.0, 0.5, 2.0), Point::new(2.2, 3.5, 0.9)] {
        let (uv, gap) = inv.nearest(p, 1e-7);
        let brute = brute_nearest(&s, p);
        assert!(gap <= brute + 1e-7, "{p:?}: {gap} vs {brute}");
        assert!(close(s.point(uv.x, uv.y).distance(p), gap, 1e-12));
    }
}

/// The distance from `p` to a dense grid of the surface, refined round its best point.
fn brute_nearest(s: &Surface, p: Point) -> f64 {
    let (u0, u1, v0, v1) = s.bounds();
    let n = 400;
    let mut best = (f64::INFINITY, 0.0, 0.0);
    for i in 0..=n {
        for j in 0..=n {
            let (u, v) = (u0 + (u1 - u0) * i as f64 / n as f64, v0 + (v1 - v0) * j as f64 / n as f64);
            let d = s.point(u, v).distance(p);
            if d < best.0 {
                best = (d, u, v);
            }
        }
    }
    // (A finer grid round it.)
    let (hu, hv) = ((u1 - u0) / n as f64, (v1 - v0) / n as f64);
    for i in -50..=50 {
        for j in -50..=50 {
            let u = (best.1 + hu * i as f64 / 50.0).clamp(u0, u1);
            let v = (best.2 + hv * j as f64 / 50.0).clamp(v0, v1);
            best.0 = best.0.min(s.point(u, v).distance(p));
        }
    }
    best.0
}

#[test]
fn nearest_from_unwraps_on_a_cylinder() {
    let s = cylinder(1.5);
    let mut inv = SurfaceInverse::new(&s);
    let mut hint = Point2::new(6.0, 0.0);
    let mut last = 6.0;
    for k in 1..=10 {
        let angle = 6.0 + 0.2 * k as f64;
        let p = Point::new(1.5 * angle.cos(), 1.5 * angle.sin(), 0.5);
        let (uv, _) = inv.nearest_from(hint, p, 1e-7, 0.0);
        assert!(close(uv.x, angle, 1e-12), "{uv:?} {angle}");
        assert!(uv.x > last);
        last = uv.x;
        hint = uv;
    }
    assert!(last > TAU);
}

#[test]
fn singular_places_of_surfaces() {
    let sphere_s = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0));
    let mut sphere = SurfaceInverse::new(&sphere_s);
    let places = sphere.singular_places().to_vec();
    assert_eq!(places.len(), 2);
    assert!(places.iter().all(|p| p.spread == 0.0 && !p.fixed_is_u));
    assert!(close(places[0].fixed, FRAC_PI_2, 1e-15) && close(places[1].fixed, -FRAC_PI_2, 1e-15));
    assert!(places[0].point.distance(Point::new(0.0, 0.0, 2.0)) < 1e-12);

    let cone = Cone::new(Frame::WORLD, 0.3, 2.0);
    let s = Surface::Cone(cone);
    let mut inv = SurfaceInverse::new(&s);
    let places = inv.singular_places().to_vec();
    assert_eq!(places.len(), 1);
    assert!(close(places[0].fixed, -2.0 / 0.3f64.sin(), 1e-12));
    assert!(s.point(1.0, places[0].fixed).distance(cone.apex()) < 1e-12);

    let spindle = Surface::Torus(Torus::new(Frame::WORLD, 1.0, 2.0));
    let mut inv = SurfaceInverse::new(&spindle);
    let places = inv.singular_places().to_vec();
    assert_eq!(places.len(), 2);
    assert!(places.iter().all(|p| p.spread == 0.0));
    for p in &places {
        // (On the axis.)
        let q = spindle.point(0.7, p.fixed);
        assert!(q.x.hypot(q.y) < 1e-12, "{q:?}");
    }
    let ring = Surface::Torus(Torus::new(Frame::WORLD, 5.0, 1.0));
    let places = SurfaceInverse::new(&ring).singular_places().to_vec();
    assert_eq!(places.len(), 1);
    assert!(close(places[0].spread, 4.0, 1e-12) && close(places[0].fixed, PI, 1e-12));

    // A patch whose row u = 0 is one point: that place first, spread 0.
    let poles: Vec<Vec<Point>> = (0..4)
        .map(|i| (0..4).map(|j| if i == 0 { Point::new(0.0, 0.0, 1.0) } else { Point::new(i as f64, j as f64 - 1.5, 0.2 * j as f64) }).collect())
        .collect();
    let kv = || KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
    let s = Surface::BSpline(BSplineSurface::new(poles, kv(), kv()));
    let mut inv = SurfaceInverse::new(&s);
    let places = inv.singular_places().to_vec();
    assert_eq!(places.len(), 4);
    assert!(places[0].fixed_is_u && places[0].fixed == 0.0 && places[0].spread == 0.0);
    assert!(places[1].spread > 1.0);
    assert!(places[0].point.distance(Point::new(0.0, 0.0, 1.0)) < 1e-15);
}

#[test]
fn closure() {
    let full = Surface::Trimmed(TrimmedSurface::new(cylinder(1.0), 0.0, TAU, 0.0, 2.0));
    let mut inv = SurfaceInverse::new(&full);
    assert!(inv.closed_in_u(1e-7));
    assert!(!inv.closed_in_v(1e-7));
    let half = Surface::Trimmed(TrimmedSurface::new(cylinder(1.0), 0.0, PI, 0.0, 2.0));
    assert!(!SurfaceInverse::new(&half).closed_in_u(1e-7));
    // A B-spline tube (rational circles): its first and last rows of poles coincide.
    let tube = Surface::BSpline(crate::geometry::convert::surface_to_bspline(&full).expect("tube"));
    let mut inv = SurfaceInverse::new(&tube);
    assert!(inv.closed_in_u(1e-7) || inv.closed_in_v(1e-7));
    // A flat patch.
    let poles: Vec<Vec<Point>> = (0..4).map(|i| (0..4).map(|j| Point::new(i as f64, j as f64, 0.0)).collect()).collect();
    let kv = || KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
    let flat = Surface::BSpline(BSplineSurface::new(poles, kv(), kv()));
    let mut inv = SurfaceInverse::new(&flat);
    assert!(!inv.closed_in_u(1e-3) && !inv.closed_in_v(1e-3));
}

#[test]
fn feet_on_curves() {
    // A circle and a point off its plane: the foot at the point's angle.
    let frame = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 0.0, 0.0));
    let circle = Curve::Circle(Circle::new(frame, 2.0));
    let p = Point::new(1.0 + 3.0 * 2.5f64.cos(), 2.0 + 3.0 * 2.5f64.sin(), 4.0);
    let foot = foot_on_curve(&circle, p, 1e-7, 0.0, TAU, false);
    assert!(close(foot.parameter, 2.5, 1e-12));
    assert!(close(foot.distance, 2.0f64.sqrt(), 1e-12));
    // Over a range past the turn: the angle brought into it.
    let foot = foot_on_curve(&circle, p, 1e-7, 4.0, 10.0, false);
    assert!(close(foot.parameter, 2.5 + TAU, 1e-12), "{}", foot.parameter);

    // A cubic arc and a point near its end: snapped to the end within preci only.
    let arc = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(2.0, 1.0, 0.5), Point::new(3.0, 0.0, 0.0)],
        vec![0.0, 1.0],
        vec![4, 4],
        3,
    ));
    let p = Point::new(3.0, 0.0, 0.0) + Vector::new(0.0, 1e-5, 0.0);
    let snapped = foot_on_curve(&arc, p, 1e-4, 0.0, 1.0, true);
    assert_eq!(snapped.parameter, 1.0);
    let free = foot_on_curve(&arc, p, 1e-4, 0.0, 1.0, false);
    let (q, d) = arc.d1(free.parameter);
    assert!((p - q).dot(d).abs() <= 1e-12 * d.length(), "{}", free.parameter);
    assert!(free.distance < 1e-5);
    // Following from a previous foot.
    let p = arc.point(0.4) + Vector::new(0.0, 0.0, 0.01);
    let next = foot_on_curve_from(0.35, &arc, p, 1e-7, 0.0, 1.0);
    let brute = (0..=100_000).map(|i| arc.point(i as f64 / 100_000.0).distance(p)).fold(f64::INFINITY, f64::min);
    assert!(next.distance <= brute + 1e-9);

    // An open curve whose foot lies slightly beyond its end: the range widened.
    let line = Curve::Line(crate::geometry::curve::Line::new(Point::ORIGIN, Direction::new(1.0, 0.0, 0.0)));
    let p = Point::new(1.0 + 5e-8, 1.0, 0.0);
    let foot = foot_on_curve(&line, p, 1e-6, 0.0, 1.0, false);
    assert!(close(foot.parameter, 1.0 + 5e-8, 1e-15), "{}", foot.parameter);
    let p = Point::new(2.0, 1.0, 0.0);
    let foot = foot_on_curve(&line, p, 1e-6, 0.0, 1.0, false);
    assert!(foot.parameter <= 1.0 + 1e-6 + 1e-15 && foot.parameter > 1.0);
}
