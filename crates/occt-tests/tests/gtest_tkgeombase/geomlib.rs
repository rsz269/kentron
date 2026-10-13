//! GeomLib_Test.cxx.
//!
//! Not ported: CancelDenominatorDerivative (an approximation-internal reparametrization),
//! the array helpers (RemovePointsFromArray, DensifyArray1OfReal, FuseIntervals), the
//! internal isIsoLine, and GTransform (kentron has no 2D affine transformation).

use kentron::geometry::curve::{BezierCurve, Curve};
use kentron::geometry::curve2::{BezierCurve2, Curve2};
use kentron::geometry::extend::{self, SurfaceSide};
use kentron::geometry::props::{self, SurfaceNormal};
use kentron::geometry::surface::{BSplineSurface, BezierSurface, KnotVector, Plane, Sphere, Surface};
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Direction, Frame, Point, Vector};

use crate::{p, p2};

/// A biquadratic patch over [0, 1]², rational with one weight 1.75 if `rational`.
fn make_surface(rational: bool) -> BSplineSurface {
    let poles: Vec<Vec<Point>> = (0..3)
        .map(|i| {
            let x = 0.5 * i as f64;
            (0..3).map(|j| p(x, 0.5 * j as f64, x * x + 0.2 * x * (0.5 * j as f64))).collect()
        })
        .collect();
    let knots = || KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2);
    let surface = BSplineSurface::new(poles, knots(), knots());
    if rational {
        let weights = (0..3)
            .map(|i| (0..3).map(|j| if i == 1 && j == 1 { 1.75 } else { 1.0 }).collect())
            .collect();
        surface.with_weights(weights)
    } else {
        surface
    }
}

/// Degree (1, 2) surface whose `u = 0` and `u = 1` boundaries have the same points, or
/// projectively equal homogeneous poles if `projectively_equal`.
fn make_rational_boundary_surface(projectively_equal: bool) -> BSplineSurface {
    let first_weights = [1.0, 2.0, 1.0];
    let last_weights = [2.0, if projectively_equal { 4.0 } else { 1.0 }, 2.0];
    let first: Vec<Point> = (0..3).map(|j| p(1.0, 1.0 + j as f64, 0.2 * j as f64)).collect();
    let last: Vec<Point> = first
        .iter()
        .enumerate()
        .map(|(j, &q)| {
            if projectively_equal {
                q
            } else {
                Point::ORIGIN + (q - Point::ORIGIN) * (first_weights[j] / last_weights[j])
            }
        })
        .collect();
    BSplineSurface::new(
        vec![first, last],
        KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1),
        KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2),
    )
    .with_weights(vec![first_weights.to_vec(), last_weights.to_vec()])
}

fn check_extended_surface(rational: bool, in_u: bool, after: bool) {
    let context = format!("rational={rational}, in_u={in_u}, after={after}");
    let source = make_surface(rational);
    let boundary = if after { 1.0 } else { 0.0 };
    let length = 0.4;
    let continuity_tolerance = 1.0e-10;
    let length_tolerance = length * 5.0e-3;
    let uv = |t: f64, s: f64| if in_u { (s, t) } else { (t, s) };

    let params = [0.25, 0.5, 0.75];
    let reference: Vec<(Point, Vector)> = params
        .iter()
        .map(|&t| {
            let (u, v) = uv(t, boundary);
            let (q, du, dv) = source.d1(u, v);
            (q, if in_u { du } else { dv })
        })
        .collect();

    let side = match (in_u, after) {
        (true, true) => SurfaceSide::UMax,
        (true, false) => SurfaceSide::UMin,
        (false, true) => SurfaceSide::VMax,
        (false, false) => SurfaceSide::VMin,
    };
    let result = extend::surface_by_length(&Surface::BSpline(source), length, 1, side).unwrap();
    let (u_min, u_max, v_min, v_max) = result.bounds();
    let end = match side {
        SurfaceSide::UMin => u_min,
        SurfaceSide::UMax => u_max,
        SurfaceSide::VMin => v_min,
        SurfaceSide::VMax => v_max,
    };
    assert!(if after { end > boundary } else { end < boundary }, "{context}");

    for (&t, &(q0, d0)) in params.iter().zip(&reference) {
        let (u, v) = uv(t, boundary);
        let (q, du, dv) = result.d1(u, v);
        assert!(q.distance(q0) <= continuity_tolerance, "{context}");
        let d = if in_u { du } else { dv };
        assert!((d - d0).length() <= continuity_tolerance, "{context}");
        let (u, v) = uv(t, end);
        let distance = result.point(u, v).distance(q0);
        assert!((distance - length).abs() <= length_tolerance, "{context}: {distance}");
    }
    if rational {
        assert!(result.is_u_rational() || result.is_v_rational(), "{context}");
        assert!(result.weights().unwrap().iter().flatten().all(|&w| w > 0.0), "{context}");
    }
}

#[test]
fn extend_surface_by_length_all_directions_and_rationality() {
    for rational in [false, true] {
        for in_u in [false, true] {
            for after in [false, true] {
                check_extended_surface(rational, in_u, after);
            }
        }
    }
}

fn unit_segment() -> Curve {
    Curve::Bezier(BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)]))
}

#[test]
fn extend_curve_to_point_reaches_target_and_preserves_start() {
    let target = p(2.0, 0.0, 0.0);
    let curve = extend::curve_to_point(&unit_segment(), target, 1, true).unwrap();
    assert!(curve.point(curve.first_parameter()).distance(p(0.0, 0.0, 0.0)) <= 1.0e-12);
    assert!(curve.point(curve.last_parameter()).distance(target) <= 1.0e-12);
}

#[test]
fn adjust_ends_reaches_end_points_and_tangents() {
    let (first, last) = (p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0));
    let tangent = Vector::new(1.0, 0.0, 0.0);
    let curve = extend::adjust_ends(&unit_segment(), first, last, tangent, tangent).unwrap();
    for (t, expected) in [(curve.first_parameter(), first), (curve.last_parameter(), last)] {
        let (q, d) = curve.d1(t);
        assert!(q.distance(expected) <= 1.0e-12);
        assert!(d.cross(tangent).length() <= 1.0e-12 * d.length());
    }
}

#[test]
fn to_3d_preserves_rational_bezier_geometry() {
    let curve2 = BezierCurve2::rational(vec![p2(0.0, 0.0), p2(0.5, 1.0), p2(1.0, 0.0)], vec![1.0, 2.0, 1.0]);
    let Curve::Bezier(curve3) = Curve2::Bezier(curve2.clone()).to_3d(&Frame::WORLD) else {
        panic!("expected a Bezier curve");
    };
    assert!(curve3.is_rational());
    for t in [0.0, 0.3, 1.0] {
        let q = curve2.point(t);
        assert!(curve3.point(t).distance(p(q.x, q.y, 0.0)) <= 1.0e-12);
    }
}

#[test]
fn rational_bspline_closedness_uses_homogeneous_poles() {
    let different = make_rational_boundary_surface(false);
    assert!(different.point(0.0, 0.25).distance(different.point(1.0, 0.25)) > 0.1);
    assert!(!different.is_u_closed_between(0.0, 1.0, 1.0e-10));

    let equivalent = make_rational_boundary_surface(true);
    assert!(equivalent.point(0.0, 0.25).distance(equivalent.point(1.0, 0.25)) <= 1.0e-12);
    assert!(equivalent.is_u_closed_between(0.0, 1.0, 1.0e-10));
}

#[test]
fn rational_bezier_closedness_accounts_for_weights() {
    let poles: Vec<Vec<Point>> =
        (0..2).map(|_| (0..3).map(|j| p(1.0, 1.0 + j as f64, 0.0)).collect()).collect();
    let weights = (0..2)
        .map(|i| (0..3).map(|j| if i == 0 && j == 1 { 2.0 } else { 1.0 }).collect())
        .collect();
    let surface = BezierSurface::rational(poles, weights);
    assert!(surface.point(0.0, 0.25).distance(surface.point(1.0, 0.25)) > 0.1);
    assert!(!surface.is_u_closed_between(0.0, 1.0, 1.0e-10));
}

#[test]
fn surface_normal_regular_and_at_pole() {
    let plane = Surface::Plane(Plane::XY);
    let Some(SurfaceNormal::Regular(n)) = props::surface_normal(&plane, 2.0, -3.0, 1.0e-9) else {
        panic!("expected a regular normal");
    };
    assert!(n.is_parallel(Direction::Z, 1.0e-12));

    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0));
    let half_pi = 0.5 * (-1.0f64).acos();
    let Some(SurfaceNormal::Singular(n)) = props::surface_normal(&sphere, 0.4, half_pi, 1.0e-9) else {
        panic!("expected an estimated normal");
    };
    assert!(n.is_parallel(Direction::Z, 1.0e-7));
}
