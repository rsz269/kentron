//! Geom_BezierSurface_Test.cxx.
//!
//! Not ported: copy tests and `WeightsArray_*`.

use kentron::geometry::surface::BezierSurface;
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Point, Transform, Transformable, Vector};

use crate::{assert_point_near, p};

/// Poles `(i, j, (i + j) / 10)` for i, j in 1..=3.
fn quadratic() -> BezierSurface {
    BezierSurface::new(
        (1..=3).map(|i| (1..=3).map(|j| p(i as f64, j as f64, (i + j) as f64 * 0.1)).collect()).collect(),
    )
}

fn unit_square() -> Vec<Vec<Point>> {
    vec![vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)], vec![p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0)]]
}

fn rational_square(weights: [[f64; 2]; 2]) -> BezierSurface {
    BezierSurface::rational(unit_square(), weights.iter().map(|r| r.to_vec()).collect())
}

#[test]
fn evaluation() {
    let s = quadratic();
    assert_point_near(s.point(0.0, 0.0), p(1.0, 1.0, 0.2), 1e-10);
    assert_point_near(s.point(1.0, 1.0), p(3.0, 3.0, 0.6), 1e-10);
    let (_, du, dv) = s.d1(0.5, 0.5);
    assert!(du.length() > 0.0 && dv.length() > 0.0);
    assert!(s.d2(0.5, 0.5).du.length() > 0.0);
    assert!(s.dn(0.5, 0.5, 1, 0).length() > 0.0);
    assert!(s.dn(0.5, 0.5, 0, 1).length() > 0.0);
    let d3 = s.d3(0.5, 0.5);
    assert!(d3.duuu.length() <= 1e-10 && d3.dvvv.length() <= 1e-10);
}

#[test]
fn properties() {
    let s = quadratic();
    assert_eq!((s.u_degree(), s.v_degree()), (2, 2));
    assert!(!s.is_u_periodic() && !s.is_v_periodic());
    assert!(!s.is_u_rational() && !s.is_v_rational());
    assert_eq!((s.u_pole_count(), s.v_pole_count()), (3, 3));
    assert!(s.is_cn_u(100) && s.is_cn_v(100));
    assert_eq!(s.bounds(), (0.0, 1.0, 0.0, 1.0));
    assert_eq!(s.weights(), None);
    assert_eq!(s.poles().len(), 3);
    assert_eq!(s.pole(0, 0), p(1.0, 1.0, 0.2));
    assert!(BezierSurface::MAX_DEGREE >= 25);
}

#[test]
fn set_pole_and_weights() {
    let mut s = quadratic();
    s.set_pole(1, 1, p(5.0, 5.0, 5.0));
    assert_eq!(s.pole(1, 1), p(5.0, 5.0, 5.0));
    let mut r = rational_square([[1.0, 1.0], [1.0, 1.0]]);
    let before = r.point(0.5, 0.5);
    r.set_weight(1, 1, 10.0);
    assert_eq!(r.weight(1, 1), 10.0);
    assert!(r.is_u_rational() || r.is_v_rational());
    assert!(r.point(0.5, 0.5).distance(before) > 1e-6);
    let mut r = rational_square([[1.0, 1.0], [1.0, 2.0]]);
    r.set_pole_with_weight(0, 0, p(0.5, 0.5, 1.0), 5.0);
    assert_eq!((r.pole(0, 0), r.weight(0, 0)), (p(0.5, 0.5, 1.0), 5.0));
}

#[test]
fn rows_and_columns() {
    let mut s = quadratic();
    s.set_pole_row(1, &[p(10.0, 1.0, 0.0), p(10.0, 2.0, 0.0), p(10.0, 3.0, 0.0)]);
    assert!((0..3).all(|j| (s.pole(1, j).x - 10.0).abs() <= 1e-10));
    let mut s = quadratic();
    s.set_pole_col(1, &[p(1.0, 10.0, 0.0), p(2.0, 10.0, 0.0), p(3.0, 10.0, 0.0)]);
    assert!((0..3).all(|i| (s.pole(i, 1).y - 10.0).abs() <= 1e-10));
    let mut r = rational_square([[1.0, 2.0], [3.0, 4.0]]);
    r.set_weight_row(0, &[5.0, 6.0]);
    assert_eq!((r.weight(0, 0), r.weight(0, 1)), (5.0, 6.0));
    let mut r = rational_square([[1.0, 2.0], [3.0, 4.0]]);
    r.set_weight_col(0, &[7.0, 8.0]);
    assert_eq!((r.weight(0, 0), r.weight(1, 0)), (7.0, 8.0));
    let mut r = rational_square([[1.0, 2.0], [2.0, 1.0]]);
    r.set_pole_col_with_weights(1, &[p(5.0, 0.0, 0.0), p(5.0, 1.0, 0.0)], &[3.0, 4.0]);
    assert_eq!(r.pole(0, 1), p(5.0, 0.0, 0.0));
    assert_eq!((r.weight(0, 1), r.weight(1, 1)), (3.0, 4.0));
}

#[test]
fn insert_and_remove_rows_and_columns() {
    let row = [p(1.5, 1.0, 0.0), p(1.5, 2.0, 0.0), p(1.5, 3.0, 0.0)];
    let mut s = quadratic();
    s.insert_pole_row_after(0, &row);
    assert_eq!((s.u_pole_count(), s.u_degree()), (4, 3));
    s.remove_pole_row(1);
    assert_eq!(s.u_pole_count(), 3);
    let col = [p(1.0, 1.5, 0.0), p(2.0, 1.5, 0.0), p(3.0, 1.5, 0.0)];
    let mut s = quadratic();
    s.insert_pole_col_after(0, &col);
    assert_eq!((s.v_pole_count(), s.v_degree()), (4, 3));
    s.remove_pole_col(1);
    assert_eq!(s.v_pole_count(), 3);
    let mut s = quadratic();
    s.insert_pole_row_before(1, &[p(0.5, 1.0, 0.0), p(0.5, 2.0, 0.0), p(0.5, 3.0, 0.0)]);
    assert_eq!(s.u_pole_count(), 4);
    assert_eq!(s.pole(1, 0), p(0.5, 1.0, 0.0));
    let mut s = quadratic();
    s.insert_pole_col_before(1, &[p(1.0, 0.5, 0.0), p(2.0, 0.5, 0.0), p(3.0, 0.5, 0.0)]);
    assert_eq!(s.v_pole_count(), 4);
    assert_eq!(s.pole(0, 1), p(1.0, 0.5, 0.0));
    let mut r = rational_square([[1.0, 2.0], [2.0, 1.0]]);
    r.insert_pole_row_after_with_weights(0, &[p(0.0, 0.5, 0.0), p(1.0, 0.5, 0.0)], &[3.0, 4.0]);
    assert_eq!(r.u_pole_count(), 3);
    assert_eq!((r.weight(1, 0), r.weight(1, 1)), (3.0, 4.0));
}

#[test]
fn increase_degree_keeps_shape() {
    let mut s = quadratic();
    let before = s.point(0.5, 0.5);
    s.increase_degree(4, 3);
    assert_eq!((s.u_degree(), s.v_degree(), s.u_pole_count(), s.v_pole_count()), (4, 3, 5, 4));
    assert_point_near(s.point(0.5, 0.5), before, 1e-10);
    let mut poles = unit_square();
    poles[1][1] = p(1.0, 1.0, 1.0);
    let mut r = BezierSurface::rational(poles, vec![vec![1.0, 2.0], vec![2.0, 1.0]]);
    let before = r.point(0.5, 0.5);
    r.increase_degree(3, 3);
    assert_eq!((r.u_degree(), r.v_degree()), (3, 3));
    assert_point_near(r.point(0.5, 0.5), before, 1e-10);
}

#[test]
fn segment() {
    let mut s = quadratic();
    let before = s.point(0.5, 0.5);
    s.segment(0.25, 0.75, 0.25, 0.75);
    assert_point_near(s.point(0.5, 0.5), before, 1e-6);
    let poles = (1..=3).map(|i| (1..=3).map(|j| p(i as f64, j as f64, (i + j) as f64 * 0.1)).collect()).collect();
    let weights = (0..3).map(|i| (0..3).map(|j| 1.0 + 0.3 * (i + j) as f64).collect()).collect();
    let mut r = BezierSurface::rational(poles, weights);
    let mid = r.point(0.5, 0.5);
    r.segment(0.25, 0.75, 0.25, 0.75);
    assert_point_near(r.point(0.5, 0.5), mid, 1e-5);
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
    let r = rational_square([[1.0, 2.0], [2.0, 1.0]]);
    let (u_iso, v_iso) = (r.u_iso(0.5), r.v_iso(0.5));
    for i in 0..=4 {
        let t = i as f64 * 0.25;
        assert_point_near(u_iso.point(t), r.point(0.5, t), 1e-10);
        assert_point_near(v_iso.point(t), r.point(t, 0.5), 1e-10);
    }
}

#[test]
fn resolution_and_transform() {
    let (ur, vr) = quadratic().resolution(1.0);
    assert!(ur > 0.0 && vr > 0.0);
    let moved = quadratic().transformed(&Transform::translation(Vector::new(10.0, 20.0, 30.0)));
    assert!(((moved.point(0.5, 0.5) - quadratic().point(0.5, 0.5)) - Vector::new(10.0, 20.0, 30.0)).length() <= 1e-10);
}
