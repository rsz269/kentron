//! BRepFilletAPI_MakeChamfer_Test.cxx
//!
//! OCCT picks "the first edge" and "its first face" by exploration order; kentron's order
//! differs, so the edges are chosen by position instead (any box edge is alike). Where OCCT
//! only checks validity, the volumes the chamfers leave are checked too.

use kentron::blend::{chamfer, chamfer_edges, ChamferSize, ChamferedEdge};
use kentron::math::Point;
use kentron::primitives::{Cuboid, Cylinder};
use kentron::topology::{Model, NodeKind, Shape};
use kentron::{boolean, properties, validity};

fn assert_valid(model: &Model, shape: Shape) {
    let report = validity::check(model, shape);
    assert!(report.problems.is_empty(), "{report:?}");
}

fn volume(model: &Model, shape: Shape) -> f64 {
    properties::volume(model, shape).mass
}

/// The box edge along Z at x = y = 0.
fn vertical_edge_at_origin(model: &Model, shape: Shape) -> Shape {
    model
        .unique_subshapes(shape, NodeKind::Edge)
        .into_iter()
        .find(|&e| {
            let (a, b) = model.edge_vertices(e);
            let (a, b) = (model.point(a), model.point(b));
            a.x == 0.0 && a.y == 0.0 && b.x == 0.0 && b.y == 0.0
        })
        .expect("a vertical edge at the origin")
}

#[test]
fn symmetric_chamfer() {
    let mut model = Model::new();
    let block = Cuboid::new(20.0, 20.0, 20.0).build(&mut model).unwrap();
    let edge = vertical_edge_at_origin(&model, block);
    let result = chamfer(&mut model, block, &[edge], 2.0).unwrap();
    assert_valid(&model, result);
    assert!((volume(&model, result) - (8000.0 - 2.0 * 2.0 / 2.0 * 20.0)).abs() < 1e-9);
}

#[test]
fn asymmetric_chamfer() {
    let mut model = Model::new();
    let block = Cuboid::new(20.0, 20.0, 20.0).build(&mut model).unwrap();
    let edge = vertical_edge_at_origin(&model, block);
    // A face at the edge: the one at x = 0, 1 along it, 3 along the other.
    let face = model
        .unique_subshapes(block, NodeKind::Face)
        .into_iter()
        .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| model.point(v).x == 0.0))
        .unwrap();
    let edges = [ChamferedEdge { edge, size: ChamferSize::Distances(1.0, 3.0), face: Some(face) }];
    let result = chamfer_edges(&mut model, block, &edges).unwrap();
    assert_valid(&model, result);
    assert!((volume(&model, result) - (8000.0 - 1.0 * 3.0 / 2.0 * 20.0)).abs() < 1e-9);
    // Cut back 1 along the face at x = 0 (its edge at the origin moved to y = 1), 3 along
    // the face at y = 0 (to x = 3).
    let points: Vec<Point> = model.unique_subshapes(result, NodeKind::Vertex).into_iter().map(|v| model.point(v)).collect();
    assert!(points.iter().any(|p| p.x == 0.0 && (p.y - 1.0).abs() < 1e-12));
    assert!(points.iter().any(|p| p.y == 0.0 && (p.x - 3.0).abs() < 1e-12));
}

/// OCCT issue 1177: every edge of a 50 x 50 x 10 slab chamfered by 5 (the chamfers of the
/// top and bottom edges meet halfway up the sides, taking the side faces away).
#[test]
fn issue1177_chamfer_all_edges_flat_box_succeeds_without_crash() {
    let mut model = Model::new();
    let slab = Cuboid::new(50.0, 50.0, 10.0).build(&mut model).unwrap();
    let edges = model.unique_subshapes(slab, NodeKind::Edge);
    let result = chamfer(&mut model, slab, &edges, 5.0).unwrap();
    assert_valid(&model, result);
    let tolerance = model
        .unique_subshapes(result, NodeKind::Vertex)
        .into_iter()
        .chain(model.unique_subshapes(result, NodeKind::Edge))
        .map(|s| model.tolerance(s))
        .fold(0.0f64, f64::max);
    assert!(tolerance <= 1.01e-4, "{tolerance}");
}

#[test]
fn chamfer_more_faces() {
    let mut model = Model::new();
    let block = Cuboid::new(20.0, 20.0, 20.0).build(&mut model).unwrap();
    let edge = vertical_edge_at_origin(&model, block);
    let result = chamfer(&mut model, block, &[edge], 2.0).unwrap();
    assert!(model.unique_subshapes(result, NodeKind::Face).len() > 6);
}

/// Edges at vertices with three or more faces of a box fused with a cylinder: the chamfer
/// may succeed (valid) or fail, but must not crash.
#[test]
fn chamfer_after_boolean_fusion() {
    let mut model = Model::new();
    let block = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let frame = kentron::math::Frame::from_z(Point::new(5.0, 0.0, 5.0), kentron::math::Direction::Y);
    let cylinder = Cylinder::new(3.0, 10.0).with_frame(frame).build(&mut model).unwrap();
    let fused = boolean::union(&mut model, block, cylinder).unwrap();
    let faces = model.unique_subshapes(fused, NodeKind::Face);
    let faces_at = |model: &Model, v: Shape| faces.iter().filter(|&&f| model.unique_subshapes(f, NodeKind::Vertex).iter().any(|w| w.node == v.node)).count();
    let edges: Vec<ChamferedEdge> = model
        .unique_subshapes(fused, NodeKind::Edge)
        .into_iter()
        .filter(|&e| model.unique_subshapes(e, NodeKind::Vertex).into_iter().any(|v| faces_at(&model, v) >= 3))
        .take(4)
        .map(|edge| ChamferedEdge { edge, size: ChamferSize::Distances(0.5, 0.5), face: None })
        .collect();
    assert!(!edges.is_empty());
    if let Ok(result) = chamfer_edges(&mut model, fused, &edges) {
        assert_valid(&model, result);
    }
}

#[test]
fn sequential_chamfer_no_crash() {
    let mut model = Model::new();
    let mut shape = Cuboid::new(20.0, 20.0, 20.0).build(&mut model).unwrap();
    let (mut previous, mut successes) = (0, 0);
    for i in 0..3 {
        let edges = model.unique_subshapes(shape, NodeKind::Edge);
        if i > 0 && previous > 0 {
            assert_ne!(edges.len(), previous);
        }
        previous = edges.len();
        let edge = edges[(i * 3 + 1) % edges.len()];
        match chamfer(&mut model, shape, &[edge], 1.0) {
            Ok(result) => {
                shape = result;
                successes += 1;
            }
            Err(_) => break,
        }
    }
    assert!(successes >= 1);
}
