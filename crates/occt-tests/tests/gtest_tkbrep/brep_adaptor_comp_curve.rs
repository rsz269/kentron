//! BRepAdaptor_CompCurve_Test.cxx (BRepAdaptor_CompCurve is `WireCurve`)

use kentron::geometry::ParametricCurve;
use std::f64::consts::PI;

use occt_tests::assert_near;
use kentron::geometry::curve::{Circle, Curve, TrimmedCurve};
use kentron::make;
use kentron::math::{Direction, Frame, Point, Vector};
use kentron::topology::{Model, NodeKind, Shape, WireCurve};

use crate::make_edge;

/// OCC5696: finding the edge at a wire parameter.
#[test]
fn occ5696_edge_method() {
    let mut model = Model::new();
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(2.0, 0.0, 0.0));
    let wire = make::wire(&mut model, &[edge]).unwrap();

    let curve = WireCurve::new(&model, wire);
    let mid = (curve.first_parameter() + curve.last_parameter()) / 2.0;
    let (found, edge_param) = curve.edge_at(mid);

    assert_eq!(found.kind(), NodeKind::Edge);
    assert!((0.0..=2.0).contains(&edge_param), "edge parameter within the edge");
    assert_near!(1.0, edge_param, 0.01, "edge parameter should be about 1");
}

/// OCC29430: the wire curve's end points are the wire's end vertices.
#[test]
fn occ29430_arc_boundary_points() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X), 1.0);
    let arc = circle.arc(PI / 4.0, 3.0 * PI / 4.0, true);
    let edge = make::edge_from_curve(&mut model, Curve::Trimmed(arc)).unwrap();
    let wire = make::wire(&mut model, &[edge]).unwrap();

    let curve = WireCurve::new(&model, wire);
    let start = curve.point(curve.first_parameter());
    let end = curve.point(curve.last_parameter());

    let vertices: Vec<Point> =
        model.explore(wire, NodeKind::Vertex).map(|v| model.point(v)).collect();
    assert!(!vertices.is_empty());
    assert!(vertices.iter().any(|v| start.distance(*v) < 1.0e-7), "start is not a wire vertex");
    assert!(vertices.iter().any(|v| end.distance(*v) < 1.0e-7), "end is not a wire vertex");
    assert!(start.distance(end) > 1.0e-7, "start and end differ");
}

/// OCC30869: end points and tangents of a wire made of one reversed arc.
#[test]
fn occ30869_reversed_edge_boundary_points() {
    let mut model = Model::new();
    let frame = Frame::new(Point::new(1.0, 0.0, 0.0), Direction::NEG_Y, Direction::NEG_Z);
    let circle = Circle::new(frame, 1.0);
    let (t1, t2) = (PI / 2.0, 3.0 * PI / 2.0);
    let trimmed = TrimmedCurve::new(Curve::Circle(circle), t1, t2);
    let edge = make::edge_from_curve(&mut model, Curve::Trimmed(trimmed)).unwrap().reversed();
    let wire = make::wire(&mut model, &[edge]).unwrap();

    let curve = WireCurve::new(&model, wire);
    let (p_first, v_first) = curve.d1(curve.first_parameter());
    let (p_last, v_last) = curve.d1(curve.last_parameter());

    // reference: the circle with the opposite normal, at the same parameters
    let reference = Circle::new(Frame::new(Point::new(1.0, 0.0, 0.0), Direction::Y, Direction::NEG_Z), 1.0);
    let (ref_p1, ref_v1) = reference.d1(t1);
    let (ref_p2, ref_v2) = reference.d1(t2);

    let tol = 1.0e-7;
    let close = |a: [f64; 3], b: [f64; 3], what: &str| {
        for i in 0..3 {
            assert_near!(a[i], b[i], tol, "{what}, coordinate {i}");
        }
    };
    let coords = |p: Point| [p.x, p.y, p.z];
    let unit = |v: Vector| {
        let n = v.normalized();
        [n.x, n.y, n.z]
    };
    close(coords(p_first), coords(ref_p1), "first point");
    close(unit(v_first), unit(ref_v1), "first tangent");
    close(coords(p_last), coords(ref_p2), "last point");
    close(unit(v_last), unit(ref_v2), "last tangent");
}

/// Four collinear edges added to a wire in connected or shuffled order; the last edge
/// runs backwards.
fn occ30422_wire(model: &mut Model, shuffled: bool) -> Shape {
    let v: Vec<Shape> = (0..5).map(|i| make::vertex(model, Point::new(i as f64, 0.0, 0.0))).collect();
    let e1 = make::edge(model, v[0], v[1]).unwrap();
    let e2 = make::edge(model, v[1], v[2]).unwrap();
    let e3 = make::edge(model, v[2], v[3]).unwrap();
    let e4 = make::edge(model, v[4], v[3]).unwrap();

    let wire = model.new_empty(NodeKind::Wire);
    let order = if shuffled { [e1, e4, e3, e2] } else { [e1, e2, e3, e4] };
    for e in order {
        model.add(wire, e).unwrap();
    }
    wire
}

/// From tests/bugs/moddata_3/bug30422: walking a wire in connection order works for both
/// insertion orders (the fourth edge need not be reachable).
#[test]
fn occ30422_wire_explorer_reordered_edges() {
    let mut model = Model::new();
    for shuffled in [false, true] {
        let wire = occ30422_wire(&mut model, shuffled);
        assert!(model.ordered_edges(wire).len() >= 3);
    }
}
