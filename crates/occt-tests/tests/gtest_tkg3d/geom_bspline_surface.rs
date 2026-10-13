//! Geom_BSplineSurface_Test.cxx.
//!
//! Not ported: copy tests and `WeightsArray_*`.

use std::f64::consts::TAU;

use kentron::geometry::surface::{BSplineSurface, KnotVector};
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Point, Transform, Transformable, Vector};

use crate::{assert_point_near, p};

fn quadratic_knots() -> KnotVector {
    KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2)
}

fn grid(f: impl Fn(usize, usize) -> Point) -> Vec<Vec<Point>> {
    (1..=3).map(|i| (1..=3).map(|j| f(i, j)).collect()).collect()
}

/// Poles `(i, j, (i + j) / 10)`, single quadratic patch.
fn quadratic() -> BSplineSurface {
    BSplineSurface::new(grid(|i, j| p(i as f64, j as f64, (i + j) as f64 * 0.1)), quadratic_knots(), quadratic_knots())
}

/// A bump in the middle with weights all 1 except `center` at (1, 1).
fn bump(center: f64) -> BSplineSurface {
    let mut weights = vec![vec![1.0; 3]; 3];
    weights[1][1] = center;
    BSplineSurface::new(
        grid(|i, j| p(i as f64, j as f64, if i == 2 && j == 2 { 1.0 } else { 0.0 })),
        quadratic_knots(),
        quadratic_knots(),
    )
    .with_weights(weights)
}

fn periodic_knots() -> KnotVector {
    KnotVector::periodic(vec![0.0, 1.0, 2.0, 3.0], vec![1; 4], 2)
}

#[test]
fn evaluation() {
    let s = quadratic();
    assert_point_near(s.point(0.0, 0.0), p(1.0, 1.0, 0.2), 1e-10);
    assert_point_near(s.point(1.0, 1.0), p(3.0, 3.0, 0.6), 1e-10);
    let (_, du, dv) = s.d1(0.5, 0.5);
    assert!(du.length() > 0.0 && dv.length() > 0.0);
    assert!(s.d2(0.5, 0.5).du.length() > 0.0);
    assert!(s.dn(0.5, 0.5, 1, 0).length() > 0.0 && s.dn(0.5, 0.5, 0, 1).length() > 0.0);
    let d3 = s.d3(0.5, 0.5);
    assert!(d3.duuu.length() <= 1e-10 && d3.dvvv.length() <= 1e-10);
}

#[test]
fn properties_and_access() {
    let s = quadratic();
    assert_eq!((s.u_degree(), s.v_degree()), (2, 2));
    assert!(!s.is_u_periodic() && !s.is_v_periodic());
    assert!(!s.is_u_rational() && !s.is_v_rational());
    assert_eq!((s.u_pole_count(), s.v_pole_count()), (3, 3));
    assert_eq!((s.u_knots().len(), s.v_knots().len()), (2, 2));
    assert_eq!(s.bounds(), (0.0, 1.0, 0.0, 1.0));
    assert_eq!(s.u_knots(), &[0.0, 1.0]);
    assert_eq!(s.u_multiplicities(), &[3, 3]);
    assert_eq!(s.v_multiplicities().len(), 2);
    assert_eq!(s.weights(), None);
}

#[test]
fn set_pole_and_weights() {
    let mut s = quadratic();
    s.set_pole(1, 1, p(5.0, 5.0, 5.0));
    assert_eq!(s.pole(1, 1), p(5.0, 5.0, 5.0));
    let mut r = bump(2.0);
    assert!(r.is_u_rational() || r.is_v_rational());
    let before = r.point(0.5, 0.5);
    r.set_weight(1, 1, 10.0);
    assert_eq!(r.weight(1, 1), 10.0);
    assert!(r.point(0.5, 0.5).z > before.z);
    let mut r = bump(2.0);
    r.set_weight_row(1, &[3.0, 4.0, 5.0]);
    assert_eq!((r.weight(1, 0), r.weight(1, 1), r.weight(1, 2)), (3.0, 4.0, 5.0));
    let mut r = bump(2.0);
    r.set_weight_col(1, &[6.0, 7.0, 8.0]);
    assert_eq!((r.weight(0, 1), r.weight(1, 1), r.weight(2, 1)), (6.0, 7.0, 8.0));
}

#[test]
fn pole_rows_and_columns() {
    let mut s = quadratic();
    s.set_pole_row(1, &[p(10.0, 1.0, 0.0), p(10.0, 2.0, 0.0), p(10.0, 3.0, 0.0)]);
    assert!((0..3).all(|j| (s.pole(1, j).x - 10.0).abs() <= 1e-10));
    let mut s = quadratic();
    s.set_pole_col(1, &[p(1.0, 10.0, 0.0), p(2.0, 10.0, 0.0), p(3.0, 10.0, 0.0)]);
    assert!((0..3).all(|i| (s.pole(i, 1).y - 10.0).abs() <= 1e-10));
}

#[test]
fn rational_surface() {
    let mut weights = vec![vec![1.0; 3]; 3];
    weights[1][1] = 3.0;
    let s = BSplineSurface::new(grid(|i, j| p(i as f64, j as f64, 0.0)), quadratic_knots(), quadratic_knots())
        .with_weights(weights);
    assert!(s.is_u_rational() && s.is_v_rational());
    assert_eq!(s.weights().unwrap()[1][1], 3.0);
    let mid = s.point(0.5, 0.5);
    assert!((mid.x - 2.0).abs() <= 0.1 && (mid.y - 2.0).abs() <= 0.1);
}

#[test]
fn knot_operations_keep_shape() {
    let before = quadratic().point(0.5, 0.5);
    let mut s = quadratic();
    s.increase_degree(4, 3);
    assert_eq!((s.u_degree(), s.v_degree()), (4, 3));
    assert_point_near(s.point(0.5, 0.5), before, 1e-10);
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    assert_eq!(s.u_knots().len(), 3);
    assert_point_near(s.point(0.5, 0.5), before, 1e-10);
    let mut s = quadratic();
    s.insert_v_knot(0.5, 1);
    assert_eq!(s.v_knots().len(), 3);
    assert_point_near(s.point(0.5, 0.5), before, 1e-10);
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    assert_eq!(s.u_multiplicity(1), 1);
    s.increase_u_multiplicity(1, 2);
    assert_eq!(s.u_multiplicity(1), 2);
    assert_point_near(s.point(0.5, 0.5), before, 1e-10);
}

#[test]
fn remove_knots() {
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    let before = s.point(0.25, 0.5);
    assert!(s.remove_u_knot(1, 0, 1e-6));
    assert_eq!(s.u_knots().len(), 2);
    assert_point_near(s.point(0.25, 0.5), before, 1e-6);
    let mut s = quadratic();
    s.insert_v_knot(0.5, 1);
    let before = s.point(0.5, 0.25);
    assert!(s.remove_v_knot(1, 0, 1e-6));
    assert_eq!(s.v_knots().len(), 2);
    assert_point_near(s.point(0.5, 0.25), before, 1e-6);
}

#[test]
fn set_knots() {
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    s.set_u_knot(1, 0.6);
    assert_eq!(s.u_knot(1), 0.6);
    let mut s = quadratic();
    s.insert_v_knot(0.5, 1);
    s.set_v_knot(1, 0.6);
    assert_eq!(s.v_knot(1), 0.6);
}

#[test]
fn segment_keeps_parametrization() {
    let mut s = quadratic();
    let before = s.point(0.5, 0.5);
    s.segment(0.25, 0.75, 0.25, 0.75);
    let (u1, u2, v1, v2) = s.bounds();
    assert!((u1 - 0.25).abs() <= 1e-10 && (u2 - 0.75).abs() <= 1e-10);
    assert!((v1 - 0.25).abs() <= 1e-10 && (v2 - 0.75).abs() <= 1e-10);
    assert_point_near(s.point(0.5, 0.5), before, 1e-6);
}

#[test]
fn exchange_and_reverse() {
    let mut s = quadratic();
    let before = s.point(0.3, 0.7);
    s.exchange_uv();
    assert_point_near(s.point(0.7, 0.3), before, 1e-10);
    let mut s = quadratic();
    let before = s.point(0.3, 0.5);
    s.u_reverse();
    assert_point_near(s.point(0.7, 0.5), before, 1e-10);
    let mut s = quadratic();
    let before = s.point(0.5, 0.3);
    s.v_reverse();
    assert_point_near(s.point(0.5, 0.7), before, 1e-10);
}

#[test]
fn isolines() {
    let s = quadratic();
    let (u_iso, v_iso) = (s.u_iso(0.5), s.v_iso(0.5));
    for i in 0..=4 {
        let t = i as f64 * 0.25;
        assert_point_near(u_iso.point(t), s.point(0.5, t), 1e-10);
        assert_point_near(v_iso.point(t), s.point(t, 0.5), 1e-10);
    }
    let weights = (0..3).map(|i| (0..3).map(|j| 1.0 + 0.5 * (i + j) as f64).collect()).collect();
    let r = quadratic().with_weights(weights);
    let u_iso = r.u_iso(0.5);
    for i in 0..=4 {
        let t = i as f64 * 0.25;
        assert_point_near(u_iso.point(t), r.point(0.5, t), 1e-10);
    }
}

#[test]
fn resolution_and_transform() {
    let (ur, vr) = quadratic().resolution(1.0);
    assert!(ur > 0.0 && vr > 0.0);
    let moved = quadratic().transformed(&Transform::translation(Vector::new(10.0, 20.0, 30.0)));
    assert!(((moved.point(0.5, 0.5) - quadratic().point(0.5, 0.5)) - Vector::new(10.0, 20.0, 30.0)).length() <= 1e-10);
}

fn ring_in_u(rational: bool) -> BSplineSurface {
    let s = BSplineSurface::new(
        grid(|i, j| {
            let a = TAU * (i - 1) as f64 / 3.0;
            p(a.cos(), a.sin(), j as f64 * 0.5)
        }),
        periodic_knots(),
        quadratic_knots(),
    );
    if rational {
        s.with_weights((0..3).map(|i| vec![1.0 + 0.5 * i as f64; 3]).collect())
    } else {
        s
    }
}

fn ring_in_v(rational: bool) -> BSplineSurface {
    let s = BSplineSurface::new(
        grid(|i, j| {
            let a = TAU * (j - 1) as f64 / 3.0;
            p(i as f64 * 0.5, a.cos(), a.sin())
        }),
        quadratic_knots(),
        periodic_knots(),
    );
    if rational {
        s.with_weights(vec![(0..3).map(|j| 1.0 + 0.5 * j as f64).collect(); 3])
    } else {
        s
    }
}

#[test]
fn unwrapping_periodic_surfaces_keeps_shape() {
    for rational in [false, true] {
        let mut s = ring_in_u(rational);
        assert!(s.is_u_periodic());
        let before = s.point(0.5, 0.5);
        s.set_u_not_periodic();
        assert!(!s.is_u_periodic());
        assert_point_near(s.point(0.5, 0.5), before, 1e-10);
        let mut s = ring_in_v(rational);
        assert!(s.is_v_periodic());
        let before = s.point(0.5, 0.5);
        s.set_v_not_periodic();
        assert!(!s.is_v_periodic());
        assert_point_near(s.point(0.5, 0.5), before, 1e-10);
    }
}

#[test]
fn local_evaluation_matches_global() {
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    assert_point_near(s.local_point(0.25, 0.5, (0, 1), (0, 1)), s.point(0.25, 0.5), 1e-10);
}

#[test]
fn move_point() {
    let mut s = quadratic();
    let target = p(2.0, 2.0, 1.0);
    assert!(s.move_point(0.5, 0.5, target, (0, 2), (0, 2)).is_some());
    assert_point_near(s.point(0.5, 0.5), target, 1e-6);
}

#[test]
fn locate() {
    let mut s = quadratic();
    s.insert_u_knot(0.5, 1);
    let (_, i2) = s.locate_u(0.25, 1e-10);
    assert!(i2 < s.u_knots().len());
    let mut s = quadratic();
    s.insert_v_knot(0.5, 1);
    let (_, i2) = s.locate_v(0.25, 1e-10);
    assert!(i2 < s.v_knots().len());
}
