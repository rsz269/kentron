//! BRepGraph_RelatedIterator_Test.cxx (`Model::related`)

use kentron::math::{Transform, Vector};
use kentron::topology::{CoEdge, Edge, Face, Model, NodeId, Relation, Shell, Solid, Vertex, Wire};

use crate::graph_support::{box_model, compound, first, id_of};
use crate::make_box;

fn count(model: &Model, node: impl Into<NodeId>, relation: Relation) -> usize {
    model.related(node.into()).iter().filter(|(r, _)| *r == relation).count()
}

fn related_of(model: &Model, node: impl Into<NodeId>, relation: Relation) -> Vec<NodeId> {
    model.related(node.into()).into_iter().filter(|(r, _)| *r == relation).map(|(_, n)| n).collect()
}

#[test]
fn face_of_box_returns_boundary_edges_and_outer_wire() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert_eq!(count(&model, face, Relation::BoundaryEdge), 4);
    assert_eq!(count(&model, face, Relation::AdjacentFace), 4);
    assert_eq!(count(&model, face, Relation::OuterWire), 1);
    let outer = model.outer_wire(face).expect("outer wire");
    assert_eq!(related_of(&model, face, Relation::OuterWire), [NodeId::from(outer)]);
}

#[test]
fn edge_of_box_returns_incident_vertices_and_faces() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert_eq!(count(&model, edge, Relation::IncidentVertex), 2);
    assert_eq!(count(&model, edge, Relation::ReferencedByFace), 2);
}

/// Products and instances have no topological relations.
#[test]
fn assembly_nodes_yield_no_relations() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid.node).unwrap();
    let assembly = model.new_assembly();
    let instance = model
        .add_instance(assembly, part, Transform::translation(Vector::new(1.0, 2.0, 3.0)))
        .unwrap();
    assert!(model.related(assembly.into()).is_empty());
    assert!(model.related(instance.into()).is_empty());
}

#[test]
fn edge_referenced_by_face_removed_face_is_skipped() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let faces = model.ancestors_of::<Face>(edge);
    assert_eq!(faces.len(), 2);
    model.remove_node(faces[0]);
    assert_eq!(related_of(&model, edge, Relation::ReferencedByFace), [NodeId::from(faces[1])]);
}

#[test]
fn container_nodes_yield_no_relations() {
    let (model, _) = box_model();
    assert!(model.related(first::<Solid>(&model).into()).is_empty());
    assert!(model.related(first::<Shell>(&model).into()).is_empty());
}

#[test]
fn wire_of_box_returns_co_edges() {
    let (model, _) = box_model();
    assert_eq!(count(&model, first::<Wire>(&model), Relation::WireCoEdge), 4);
}

#[test]
fn wire_of_box_returns_owning_face() {
    let (model, _) = box_model();
    assert_eq!(count(&model, first::<Wire>(&model), Relation::OwningFace), 1);
}

/// Each box vertex touches 3 edges.
#[test]
fn vertex_of_box_returns_incident_edges() {
    let (model, _) = box_model();
    assert_eq!(count(&model, first::<Vertex>(&model), Relation::IncidentEdge), 3);
}

#[test]
fn co_edge_of_box_returns_parent_edge_and_owning_face() {
    let (model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    assert_eq!(count(&model, coedge, Relation::ParentEdge), 1);
    assert_eq!(count(&model, coedge, Relation::OwningFace), 1);
}

#[test]
fn compound_yields_no_relations() {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 10.0, 10.0);
    let box2 = make_box(&mut model, 20.0, 20.0, 20.0);
    let c = compound(&mut model, &[box1, box2]);
    assert!(model.related(id_of::<kentron::topology::Compound>(c).into()).is_empty());
}

/// 4 boundary edges + 4 adjacent faces + 1 outer wire at least.
#[test]
fn range_for_works_correctly() {
    let (model, _) = box_model();
    assert!(model.related(first::<Face>(&model).into()).len() >= 9);
}

#[test]
fn edge_of_box_all_relations_sequential() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    let faces = related_of(&model, edge, Relation::ReferencedByFace);
    let vertices = related_of(&model, edge, Relation::IncidentVertex);
    assert_eq!(faces.len(), 2);
    assert_eq!(vertices.len(), 2);
    assert_ne!(faces[0], faces[1]);
    assert_ne!(vertices[0], vertices[1]);
}

#[test]
fn edge_of_box_removed_face_correct_transition() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let faces = model.ancestors_of::<Face>(edge);
    assert_eq!(faces.len(), 2);
    model.remove_node(faces[0]);
    assert_eq!(count(&model, edge, Relation::ReferencedByFace), 1);
    assert_eq!(count(&model, edge, Relation::IncidentVertex), 2);
}

#[test]
fn vertex_of_box_all_parent_edges_yielded() {
    let (model, _) = box_model();
    let edges = related_of(&model, first::<Vertex>(&model), Relation::IncidentEdge);
    assert_eq!(edges.len(), 3);
    for i in 0..edges.len() {
        for j in i + 1..edges.len() {
            assert_ne!(edges[i], edges[j]);
        }
    }
}
