//! BRep_Tool_Test.cxx

use kentron::geometry::ParametricCurve;
use occt_tests::{assert_double_eq, assert_near};
use kentron::geometry::curve::{BezierCurve, Circle, Curve};
use kentron::geometry::surface::Plane;
use kentron::make;
use kentron::math::{Direction, Frame, Point};
use kentron::tolerance;
use kentron::topology::{FaceId, Model, NodeKind};

use crate::{first, make_box, make_edge};

fn unit_circle() -> Circle {
    Circle::new(Frame::from_z(Point::ORIGIN, Direction::Z), 1.0)
}

#[test]
fn pnt_from_vertex() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let p = model.point(first(&model, solid, NodeKind::Vertex));

    // within the box [0,10] x [0,20] x [0,30]
    let tol = tolerance::LINEAR;
    assert!(-tol <= p.x && p.x <= 10.0 + tol);
    assert!(-tol <= p.y && p.y <= 20.0 + tol);
    assert!(-tol <= p.z && p.z <= 30.0 + tol);
}

#[test]
fn tolerance_vertex() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert!(model.tolerance(first(&model, solid, NodeKind::Vertex)) >= 0.0);
}

#[test]
fn tolerance_edge() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert!(model.tolerance(first(&model, solid, NodeKind::Edge)) >= 0.0);
}

#[test]
fn tolerance_face() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert!(model.tolerance(first(&model, solid, NodeKind::Face)) >= 0.0);
}

#[test]
fn curve_from_edge() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let (_, first_param, last_param) =
        model.edge_curve(first(&model, solid, NodeKind::Edge)).expect("a box edge has a curve");
    assert!(first_param < last_param);
}

#[test]
fn surface_from_face() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let face = first(&model, solid, NodeKind::Face);
    assert!(model.face_surface(FaceId::try_from(face.node).unwrap()).is_some());
}

#[test]
fn is_closed_circle_edge() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::from_z(Point::ORIGIN, Direction::Z), 5.0);
    let edge = make::edge_from_curve(&mut model, Curve::Circle(circle)).unwrap();
    assert!(model.is_closed_edge(edge), "a full circle edge is closed");
}

#[test]
fn is_closed_line_edge() {
    let mut model = Model::new();
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
    assert!(!model.is_closed_edge(edge), "a segment is not closed");
}

#[test]
fn degenerated() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert!(!model.is_degenerate(first(&model, solid, NodeKind::Edge)));
}

#[test]
fn curve_on_plane_rejects_edge_range_outside_bounded_curve() {
    let mut model = Model::new();
    let bezier = BezierCurve::new(vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.5, 0.25, 0.0),
        Point::new(1.0, 0.0, 0.0),
    ]);
    let edge = make::edge_from_curve(&mut model, Curve::Bezier(bezier)).unwrap();

    // slightly outside the curve's domain [0, 1]
    const FIRST: f64 = -0.000002673149646;
    const LAST: f64 = 0.999997326850354;
    model.set_edge_range(edge, FIRST, LAST).unwrap();

    assert!(model.pcurve_on_plane(edge, &Plane::XY).is_none());
    let (first_param, last_param) = model.edge_range(edge);
    assert_double_eq!(first_param, FIRST);
    assert_double_eq!(last_param, LAST);
}

#[test]
fn curve_on_plane_projects_shifted_periodic_range() {
    let mut model = Model::new();
    let circle = unit_circle();
    let edge = make::edge_from_curve(&mut model, Curve::Circle(circle)).unwrap();

    let (range_first, range_last) = (circle.period().unwrap(), 2.0 * circle.period().unwrap());
    model.set_edge_range(edge, range_first, range_last).unwrap();

    let (_, first_param, last_param) =
        model.pcurve_on_plane(edge, &Plane::XY).expect("projection failed");
    assert_double_eq!(first_param, range_first);
    assert_double_eq!(last_param, range_last);
}

#[test]
fn curve_on_plane_projects_periodic_edge_built_from_equal_parameters() {
    let mut model = Model::new();
    let circle = unit_circle();
    let edge = make::edge_from_curve_range(&mut model, Curve::Circle(circle), 0.0, 0.0).unwrap();

    let (edge_first, edge_last) = model.edge_range(edge);
    assert_near!(edge_last - edge_first, circle.period().unwrap(), tolerance::PARAMETRIC);

    let (_, first_param, last_param) =
        model.pcurve_on_plane(edge, &Plane::XY).expect("projection failed");
    assert_double_eq!(first_param, edge_first);
    assert_double_eq!(last_param, edge_last);
}

#[test]
fn curve_on_plane_rejects_equal_periodic_range() {
    let mut model = Model::new();
    let edge = make::edge_from_curve(&mut model, Curve::Circle(unit_circle())).unwrap();
    model.set_edge_range(edge, 0.0, 0.0).unwrap();

    assert!(model.pcurve_on_plane(edge, &Plane::XY).is_none());
    let (first_param, last_param) = model.edge_range(edge);
    assert_double_eq!(first_param, 0.0);
    assert_double_eq!(last_param, 0.0);
}
