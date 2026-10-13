//! GProp_PEquation_Test.cxx, GProp_PGProps_Test.cxx.

use kentron::properties::{self, PointSetShape};
use kentron::tolerance::LINEAR;

use crate::{assert_near, assert_point_near, p};

// ---- GProp_PEquation ----

#[test]
fn coincident_points() {
    let q = p(1.0, 2.0, 3.0);
    let PointSetShape::Point(c) = properties::point_set_shape(&[q, q, q], 1e-6) else {
        panic!("expected a point");
    };
    assert_point_near(c, q, LINEAR);
}

#[test]
fn collinear_points() {
    let points = [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0)];
    let PointSetShape::Line(line) = properties::point_set_shape(&points, 1e-6) else {
        panic!("expected a line");
    };
    assert_near(line.direction.x().abs(), 1.0, LINEAR);
}

#[test]
fn coplanar_points() {
    let points = [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0)];
    let PointSetShape::Plane(plane) = properties::point_set_shape(&points, 1e-6) else {
        panic!("expected a plane");
    };
    assert_near(plane.normal().z().abs(), 1.0, LINEAR);
}

#[test]
fn space_points() {
    let points = [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 0.0, 1.0)];
    let PointSetShape::Space { edges, .. } = properties::point_set_shape(&points, 1e-6) else {
        panic!("expected a box");
    };
    assert!(edges.iter().all(|e| e.length() > LINEAR));
}

// ---- GProp_PGProps ----

#[test]
fn empty_set() {
    assert_near(properties::points(&[]).mass, 0.0, LINEAR);
}

#[test]
fn single_point() {
    let props = properties::points(&[p(1.0, 2.0, 3.0)]);
    assert_near(props.mass, 1.0, LINEAR);
    assert_point_near(props.center_of_mass, p(1.0, 2.0, 3.0), LINEAR);
}

#[test]
fn two_points() {
    let props = properties::points(&[p(0.0, 0.0, 0.0), p(2.0, 4.0, 6.0)]);
    assert_near(props.mass, 2.0, LINEAR);
    assert_point_near(props.center_of_mass, p(1.0, 2.0, 3.0), LINEAR);
}

#[test]
fn weighted_points() {
    let props = properties::point_masses(&[(p(0.0, 0.0, 0.0), 1.0), (p(4.0, 0.0, 0.0), 3.0)]);
    assert_near(props.mass, 4.0, LINEAR);
    assert_point_near(props.center_of_mass, p(3.0, 0.0, 0.0), LINEAR);
}

fn cross_points() -> [kentron::math::Point; 4] {
    [p(1.0, 0.0, 0.0), p(-1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, -1.0, 0.0)]
}

#[test]
fn symmetric_points() {
    let props = properties::points(&cross_points());
    assert_near(props.mass, 4.0, LINEAR);
    assert_point_near(props.center_of_mass, p(0.0, 0.0, 0.0), LINEAR);
}

#[test]
fn barycenter() {
    let g = properties::barycenter(&[p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(0.0, 6.0, 0.0)]);
    assert_point_near(g, p(1.0, 2.0, 0.0), LINEAR);
}

#[test]
fn matrix_of_inertia_symmetric_points() {
    let m = properties::points(&cross_points()).inertia;
    assert_near(m[(0, 0)], 2.0, LINEAR);
    assert_near(m[(1, 1)], 2.0, LINEAR);
    assert_near(m[(2, 2)], 4.0, LINEAR);
    assert_near(m[(0, 1)], 0.0, LINEAR);
    assert_near(m[(0, 2)], 0.0, LINEAR);
    assert_near(m[(1, 2)], 0.0, LINEAR);
}
