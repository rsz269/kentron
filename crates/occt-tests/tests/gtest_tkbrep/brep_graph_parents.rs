//! BRepGraph_ReverseIterator_Test.cxx (`Model::parents`, `Model::ancestors_of`)
//!
//! Not ported: `OccurrencesOfChild_ProductAndTopologyChildren` (OCCT wraps a part's
//! topology in an occurrence; a kentron part refers to its topology root directly).

use kentron::make;
use kentron::math::Point;
use kentron::topology::{
    CoEdge, Edge, Face, Model, NodeId, NodeKind, Shape, Shell, Vertex, Wire,
};

use crate::graph_support::{box_model, compound, first, id_of};

fn parents_of_kind(model: &Model, node: impl Into<NodeId>, kind: NodeKind) -> Vec<NodeId> {
    model.parents(node.into()).into_iter().filter(|p| p.kind == kind).collect()
}

/// Every box edge is shared by 2 faces.
#[test]
fn faces_of_edge_box_edge_shared_by_two_faces() {
    let (model, _) = box_model();
    assert_eq!(model.ancestors_of::<Face>(first::<Edge>(&model)).len(), 2);
}

/// Every box vertex is shared by 3 edges.
#[test]
fn edges_of_vertex_box_vertex_shared_by_three_edges() {
    let (model, _) = box_model();
    assert_eq!(model.ancestors_of::<Edge>(first::<Vertex>(&model)).len(), 3);
}

#[test]
fn solids_of_shell_box_shell_has_one_solid() {
    let (model, _) = box_model();
    assert_eq!(parents_of_kind(&model, first::<Shell>(&model), NodeKind::Solid).len(), 1);
}

#[test]
fn shells_of_face_box_face_has_one_shell() {
    let (model, _) = box_model();
    assert_eq!(parents_of_kind(&model, first::<Face>(&model), NodeKind::Shell).len(), 1);
}

#[test]
fn faces_of_wire_box_wire_has_one_face() {
    let (model, _) = box_model();
    assert_eq!(parents_of_kind(&model, first::<Wire>(&model), NodeKind::Face).len(), 1);
}

#[test]
fn wires_of_edge_box_edge_belongs_to_two_wires() {
    let (model, _) = box_model();
    assert_eq!(model.ancestors_of::<Wire>(first::<Edge>(&model)).len(), 2);
}

#[test]
fn sequential_iteration_skips_removed_parent() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let faces = model.ancestors_of::<Face>(edge);
    assert_eq!(faces.len(), 2);
    model.remove_node(faces[0]);
    assert_eq!(model.ancestors_of::<Face>(edge).len(), 1);
}

#[test]
fn indexed_access_tracks_active_reverse_bucket() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let faces = model.ancestors_of::<Face>(edge);
    model.remove_node(faces[0]);
    let remaining = model.ancestors_of::<Face>(edge);
    assert_eq!(remaining, [faces[1]]);
    assert!(model.get(remaining[0]).is_some());
}

#[test]
fn range_for_works_correctly() {
    let (model, _) = box_model();
    let faces = model.ancestors_of::<Face>(first::<Edge>(&model));
    assert!(faces.iter().all(|&f| model.get(f).is_some()));
    assert_eq!(faces.len(), 2);
}

#[test]
fn refs_shells_of_face_returns_valid_ref_id() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let shell = first::<Shell>(&model);
    assert_eq!(parents_of_kind(&model, face, NodeKind::Shell), [NodeId::from(shell)]);
    assert!(model[shell].faces.iter().any(|&u| model[u].face == face));
}

#[test]
fn refs_edges_of_vertex_returns_valid_ref_id() {
    let (model, _) = box_model();
    let vertex = first::<Vertex>(&model);
    let edges = model.ancestors_of::<Edge>(vertex);
    assert_eq!(edges.len(), 3);
    for edge in edges {
        let record = &model[edge];
        assert!(model[record.start].vertex == vertex || model[record.end].vertex == vertex);
    }
}

/// Each box edge has a coedge on each of its 2 faces.
#[test]
fn co_edges_of_edge_box_edge_has_co_edges() {
    let (model, _) = box_model();
    assert_eq!(parents_of_kind(&model, first::<Edge>(&model), NodeKind::CoEdge).len(), 2);
}

#[test]
fn wires_of_co_edge_box_co_edge_belongs_to_one_wire() {
    let (model, _) = box_model();
    assert!(model.get(model[first::<CoEdge>(&model)].wire).is_some());
}

#[test]
fn definition_returns_face_definition() {
    let (model, _) = box_model();
    let face = model.ancestors_of::<Face>(first::<Edge>(&model))[0];
    assert!(!model[face].wires.is_empty());
}

#[test]
fn definition_returns_edge_definition() {
    let (model, _) = box_model();
    let edge = model.ancestors_of::<Edge>(first::<Vertex>(&model))[0];
    assert!(model.get(model[edge].start).is_some());
    assert!(model.get(model[edge].end).is_some());
}

#[test]
fn skips_removed_parent_edges_of_vertex() {
    let (mut model, _) = box_model();
    let vertex = first::<Vertex>(&model);
    let edges = model.ancestors_of::<Edge>(vertex);
    assert_eq!(edges.len(), 3);
    model.remove_node(edges[0]);
    assert_eq!(model.ancestors_of::<Edge>(vertex).len(), 2);
}

#[test]
fn skips_removed_parent_all_removed() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    for face in model.ancestors_of::<Face>(edge) {
        model.remove_node(face);
    }
    assert!(model.ancestors_of::<Face>(edge).is_empty());
}

#[test]
fn starting_index_skips_to_position() {
    let (model, _) = box_model();
    let edges = model.ancestors_of::<Edge>(first::<Vertex>(&model));
    assert_eq!(edges.iter().skip(1).count(), 2);
}

/// A compound holding one node of every kind is a parent of each of them.
#[test]
fn compounds_of_child_all_topology_kinds() {
    let mut model = Model::new();
    let v0 = make::vertex(&mut model, Point::ORIGIN);
    let v1 = make::vertex(&mut model, Point::new(1.0, 0.0, 0.0));
    let edge = make::edge(&mut model, v0, v1).unwrap();
    let wire = model.new_empty(NodeKind::Wire);
    model.add(wire, edge).unwrap();
    let coedge = Shape::from(model[id_of::<Wire>(wire)].coedges[0]);
    let face = model.new_empty(NodeKind::Face);
    model.add(face, wire).unwrap();
    let shell = model.new_empty(NodeKind::Shell);
    model.add(shell, face).unwrap();
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell).unwrap();
    let comp_solid = model.new_empty(NodeKind::CompSolid);
    model.add(comp_solid, solid).unwrap();
    let nested = compound(&mut model, &[solid]);
    let c = compound(&mut model, &[v0, coedge, edge, wire, face, shell, solid, comp_solid, nested]);

    let compounds_of = |node: Shape| parents_of_kind(&model, node.node, NodeKind::Compound);
    for child in [v0, coedge, edge, wire, face, shell, comp_solid, nested] {
        assert_eq!(compounds_of(child), [c.node], "{:?}", child.kind());
    }
    assert_eq!(compounds_of(solid), [nested.node, c.node], "in creation order");
    assert_eq!(parents_of_kind(&model, wire.node, NodeKind::Face), [face.node]);
}

