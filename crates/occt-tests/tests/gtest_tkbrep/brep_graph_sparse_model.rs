//! BRepGraph_SparseModel_Test.cxx: compaction, deduplication and removal on partial
//! topology (free vertices, edges without curves, empty containers).
//!
//! OCCT's `CleanupRemovedReferences` is part of `Model::remove_node` in kentron.

use kentron::geometry::curve::{Curve, Line};
use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Point};
use kentron::topology::{
    CoEdge, DeduplicateOptions, Edge, EdgeId, FaceId, Model, NodeKind, Orientation, Product,
    Shape, VertexId, WireId,
};
use kentron::validity;

use crate::graph_support::id_of;

fn vertex(model: &mut Model, x: f64, y: f64) -> VertexId {
    model.add_vertex(Point::new(x, y, 0.0), 1.0e-7)
}

/// A curveless edge from (0, y) to (10, y).
fn bare_edge(model: &mut Model, y: f64) -> EdgeId {
    let a = vertex(model, 0.0, y);
    let b = vertex(model, 10.0, y);
    model.add_edge(a, b, 1.0e-7).unwrap()
}

/// An edge on a line from (0, y) to (10, y).
fn line_edge(model: &mut Model, y: f64) -> EdgeId {
    let edge = bare_edge(model, y);
    let line = Curve::Line(Line::new(Point::new(0.0, y, 0.0), Direction::X));
    model.set_edge_curve(edge.into(), line, 1.0e-7).unwrap();
    model.set_edge_range(edge.into(), 0.0, 10.0).unwrap();
    edge
}

fn single_coedge_wire(model: &mut Model, edge: EdgeId) -> WireId {
    let coedge = model.add_coedge(edge, Orientation::Forward).unwrap();
    model.add_wire(&[coedge]).unwrap()
}

/// A face on the XY plane bounded by `wire`.
fn face_on_xy(model: &mut Model, wire: WireId) -> FaceId {
    let face = model.new_empty(NodeKind::Face);
    model
        .set_face_surface(face, Surface::Plane(Plane::new(Point::ORIGIN, Direction::Z)), 1.0e-7)
        .unwrap();
    model.add(face, wire.into()).unwrap();
    id_of(face)
}

/// vertex - edge - coedge - wire - face - shell - solid, plus a part holding the solid.
fn canonical_chain(model: &mut Model, edge: EdgeId) -> (Shape, kentron::topology::ProductId) {
    let wire = single_coedge_wire(model, edge);
    let face = face_on_xy(model, wire);
    let shell = model.new_empty(NodeKind::Shell);
    model.add(shell, face.into()).unwrap();
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell).unwrap();
    let part = model.new_part(solid).unwrap();
    (solid, part)
}

fn assert_compact_is_noop(model: &mut Model) {
    let counts: Vec<usize> = NodeKind::ALL.iter().map(|&k| model.count(k)).collect();
    let report = model.compact();
    assert_eq!(report.nodes_before, report.nodes_after);
    let after: Vec<usize> = NodeKind::ALL.iter().map(|&k| model.count(k)).collect();
    assert_eq!(after, counts);
    assert!(validity::check_model_structure(model).is_valid());
}

// ---- compaction of partial topology ----

#[test]
fn compact_bare_vertices_noop() {
    let mut model = Model::new();
    for x in [0.0, 1.0, 2.0] {
        vertex(&mut model, x, 0.0);
    }
    assert_eq!(model.count(NodeKind::Vertex), 3);
    assert_compact_is_noop(&mut model);
}

#[test]
fn compact_free_edge_without_curve_noop() {
    let mut model = Model::new();
    bare_edge(&mut model, 0.0);
    assert_eq!(model.count(NodeKind::Vertex), 2);
    assert_eq!(model.count(NodeKind::Edge), 1);
    assert_compact_is_noop(&mut model);
}

#[test]
fn compact_free_wire_single_coedge_noop() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    single_coedge_wire(&mut model, edge);
    assert_eq!(model.count(NodeKind::Wire), 1);
    assert!(model.count(NodeKind::CoEdge) >= 1);
    assert_compact_is_noop(&mut model);
}

#[test]
fn compact_empty_containers_noop() {
    for kind in [NodeKind::Compound, NodeKind::Shell, NodeKind::Solid, NodeKind::Face] {
        let mut model = Model::new();
        model.new_empty(kind);
        assert_eq!(model.count(kind), 1);
        assert_compact_is_noop(&mut model);
    }
}

#[test]
fn compact_remove_vertex_from_free_edge() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let (start, _) = model.edge_end_vertices(edge);
    model.remove_node(start);
    assert_eq!(model.count(NodeKind::Vertex), 1);
    assert_eq!(model.count(NodeKind::Edge), 1);
    model.compact();
    assert_eq!(model.count(NodeKind::Vertex), 1);
    assert_eq!(model.count(NodeKind::Edge), 1);
}

#[test]
fn compact_remove_edge_from_free_wire_retires_coedge() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let wire = single_coedge_wire(&mut model, edge);
    let coedges = model.count(NodeKind::CoEdge);
    model.remove_node(edge);
    assert_eq!(model.count(NodeKind::Edge), 0);
    assert_eq!(model.count(NodeKind::CoEdge), coedges - 1);
    assert!(model[wire].coedges.is_empty());
    model.compact();
    assert_eq!(model.count(NodeKind::CoEdge), 0);
    assert_eq!(model.count(NodeKind::Edge), 0);
    assert_eq!(model.count(NodeKind::Wire), 1);
    assert_eq!(model.count(NodeKind::Vertex), 2);
}

#[test]
fn compact_remove_face_free_wire_coedges_survive() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let wire = single_coedge_wire(&mut model, edge);
    let face = face_on_xy(&mut model, wire);
    let coedges = model.count(NodeKind::CoEdge);
    model.remove_node(face);
    assert_eq!(model.count(NodeKind::Face), 0);
    assert_eq!(model.count(NodeKind::CoEdge), coedges);
    assert_eq!(model.count(NodeKind::Wire), 1);
    model.compact();
    assert_eq!(model.count(NodeKind::Face), 0);
    assert_eq!(model.count(NodeKind::CoEdge), coedges);
    assert_eq!(model.count(NodeKind::Wire), 1);
    assert_eq!(model.count(NodeKind::Edge), 1);
    assert_eq!(model.count(NodeKind::Vertex), 2);
}

#[test]
fn compact_remove_part_shape_empties_part() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let (solid, part) = canonical_chain(&mut model, edge);
    assert_eq!(model.count(NodeKind::Product), 1);
    assert_eq!(model.count(NodeKind::Solid), 1);
    model.remove_node(solid.node);
    assert_eq!(model.count(NodeKind::Solid), 0);
    assert_eq!(model[part].shape, None);
    assert_eq!(model.count(NodeKind::Product), 1);
    model.compact();
    assert_eq!(model.count(NodeKind::Solid), 0);
    assert_eq!(model.count(NodeKind::Product), 1);
    let (_, product) = model.table::<Product>().iter().next().unwrap();
    assert_eq!(product.shape, None);
}

#[test]
fn compact_multiple_removals_validates() {
    let mut model = Model::new();
    let extra = vertex(&mut model, 100.0, 0.0);
    let free = bare_edge(&mut model, 0.0);
    let edge = bare_edge(&mut model, 1.0);
    let wire = single_coedge_wire(&mut model, edge);
    let face = face_on_xy(&mut model, wire);
    let vertices = model.count(NodeKind::Vertex);
    model.remove_node(extra);
    model.remove_node(free);
    model.remove_node(face);
    assert_eq!(model.count(NodeKind::Vertex), vertices - 1);
    assert_eq!(model.count(NodeKind::Edge), 1);
    assert_eq!(model.count(NodeKind::Face), 0);
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

// ---- deduplication of partial topology ----

#[test]
fn dedup_bare_vertices_no_rewrite() {
    let mut model = Model::new();
    for x in [0.0, 1.0, 2.0] {
        vertex(&mut model, x, 0.0);
    }
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 0);
    assert_eq!(report.distinct_curves, 0);
    assert_eq!(report.shared_surfaces, 0);
    assert_eq!(report.shared_curves, 0);
    assert!(!report.merge_applied);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn dedup_two_identical_vertices_merge_when_asked() {
    let mut model = Model::new();
    vertex(&mut model, 0.0, 0.0);
    vertex(&mut model, 0.0, 0.0);
    assert_eq!(model.count(NodeKind::Vertex), 2);
    let report =
        model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    assert!(report.merged(NodeKind::Vertex) >= 1);
    assert_eq!(model.count(NodeKind::Vertex), 1);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn dedup_two_identical_vertices_no_merge_by_default() {
    let mut model = Model::new();
    vertex(&mut model, 0.0, 0.0);
    vertex(&mut model, 0.0, 0.0);
    let report = model.deduplicate();
    assert_eq!(report.merged(NodeKind::Vertex), 0);
    assert!(!report.merge_applied);
    assert_eq!(model.count(NodeKind::Vertex), 2);
}

#[test]
fn dedup_free_edges_with_curve_share_curves() {
    let mut model = Model::new();
    let a = line_edge(&mut model, 0.0);
    let b = bare_edge(&mut model, 1.0);
    // Same (infinite) line as `a`, used over another range: equal geometry.
    let line = Curve::Line(Line::new(Point::ORIGIN, Direction::X));
    model.set_edge_curve(b.into(), line, 1.0e-7).unwrap();
    model.set_edge_range(b.into(), 0.0, 10.0).unwrap();
    assert_ne!(model[a].curve, model[b].curve);
    let report = model.deduplicate();
    assert_eq!(report.distinct_curves, 1);
    assert_eq!(report.shared_curves, 1);
    assert_eq!(model[a].curve, model[b].curve);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn dedup_free_wire_no_geometry_rewrite() {
    let mut model = Model::new();
    let edge = line_edge(&mut model, 0.0);
    single_coedge_wire(&mut model, edge);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 0);
    assert_eq!(report.distinct_curves, 1);
    assert_eq!(report.shared_curves, 0);
}

#[test]
fn dedup_empty_compound_noop() {
    let mut model = Model::new();
    model.new_empty(NodeKind::Compound);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 0);
    assert_eq!(report.distinct_curves, 0);
    assert!(validity::check_model_structure(&model).is_valid());
}

// ---- removal ----

#[test]
fn remove_face_clears_coedge_face_and_keeps_coedge() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let wire = single_coedge_wire(&mut model, edge);
    let face = face_on_xy(&mut model, wire);
    model.remove_node(face);
    for (_, coedge) in model.table::<CoEdge>().iter() {
        assert!(!coedge.face.is_valid() || model.get(coedge.face).is_none());
        assert_eq!(coedge.edge, edge);
    }
    assert_eq!(model.count(NodeKind::CoEdge), 1);
}

#[test]
fn remove_edge_retires_free_wire_coedge() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    single_coedge_wire(&mut model, edge);
    model.remove_node(edge);
    assert_eq!(model.count(NodeKind::CoEdge), 0);
    assert!(model.slot_count(NodeKind::CoEdge) >= 1);
}

#[test]
fn remove_vertex_clears_edge_vertex_use() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let (start, _) = model.edge_end_vertices(edge);
    model.remove_node(start);
    assert_eq!(model.count(NodeKind::Vertex), 1);
    assert_eq!(model.count(NodeKind::Edge), 1);
    for (_, e) in model.table::<Edge>().iter() {
        assert!(model.get(e.start).is_none(), "start vertex use should be gone");
        assert!(model.get(e.end).is_some(), "end vertex use must persist");
    }
}

#[test]
fn remove_wire_prunes_face_wire_use() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let wire = single_coedge_wire(&mut model, edge);
    let face = face_on_xy(&mut model, wire);
    model.remove_node(wire);
    assert_eq!(model.count(NodeKind::Wire), 0);
    assert_eq!(model.count(NodeKind::Face), 1);
    assert!(model[face].wires.is_empty());
}

// ---- products and pipelines ----

#[test]
fn empty_product_compact_and_audit() {
    let mut model = Model::new();
    model.new_assembly();
    assert_eq!(model.count(NodeKind::Product), 1);
    assert_eq!(model.count(NodeKind::Instance), 0);
    assert_compact_is_noop(&mut model);
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn full_chain_compact_and_dedup_no_regression() {
    let mut model = Model::new();
    let edge = line_edge(&mut model, 0.0);
    canonical_chain(&mut model, edge);
    model.compact();
    assert!(validity::check_model(&model).is_valid());
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.distinct_curves, 1);
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn pipeline_remove_compact_valid() {
    let mut model = Model::new();
    let edge = bare_edge(&mut model, 0.0);
    let extra = vertex(&mut model, 100.0, 0.0);
    let wire = single_coedge_wire(&mut model, edge);
    let face = face_on_xy(&mut model, wire);
    model.remove_node(face);
    model.remove_node(extra);
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn pipeline_remove_dedup_compact_valid() {
    let mut model = Model::new();
    let a = line_edge(&mut model, 0.0);
    line_edge(&mut model, 1.0);
    let (start, _) = model.edge_end_vertices(a);
    model.remove_node(start);
    model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    model.compact();
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn compound_with_empty_shells_compact_and_audit() {
    let mut model = Model::new();
    let a = model.new_empty(NodeKind::Shell);
    let b = model.new_empty(NodeKind::Shell);
    let compound = model.new_empty(NodeKind::Compound);
    model.add_all(compound, &[a, b]).unwrap();
    model.remove_node(a.node);
    assert_eq!(model.count(NodeKind::Shell), 1);
    model.compact();
    assert_eq!(model.count(NodeKind::Shell), 1);
    assert_eq!(model.count(NodeKind::Compound), 1);
    assert!(validity::check_model(&model).is_valid());
}

// ---- BRepGraph_EdgeCases_Test.cxx ----
//
// Not ported: TopoDS import of null shapes, rebuild-in-place tests (a rebuilt model is a
// new `Model`), parallel import.

#[test]
fn new_model_is_empty_with_zero_counts() {
    let model = Model::new();
    assert!(model.is_empty());
    for kind in NodeKind::ALL {
        assert_eq!(model.count(kind), 0);
        assert_eq!(model.slot_count(kind), 0);
    }
}

#[test]
fn empty_compound_counts() {
    let mut model = Model::new();
    model.new_empty(NodeKind::Compound);
    assert!(!model.is_empty());
    for kind in NodeKind::ALL {
        assert_eq!(model.count(kind), usize::from(kind == NodeKind::Compound));
    }
}

#[test]
fn invalid_ids_resolve_to_nothing() {
    let model = Model::new();
    assert!(model.get(kentron::topology::FaceId::from_index(7)).is_none());
    assert!(model.get(kentron::topology::FaceId::INVALID).is_none());
    assert_eq!(model.stamp(kentron::topology::FaceId::from_index(7)), None);
}
