//! BRepGraph_Views_Test.cxx: the graph's query "views".
//!
//! Most of that suite checks that OCCT's grouped accessors agree with each other (defs vs
//! topo vs refs views, UID tables, cache registry, mutation guards, TopoDS shape views);
//! kentron has one query API. Ported: the topological queries and how they treat removed
//! records.

use std::collections::HashSet;

use kentron::primitives::Sphere;
use kentron::topology::{
    CoEdge, CoEdgeId, DeduplicateOptions, Edge, EdgeId, Face, FaceId, Model, Relation,
};
use kentron::validity;

use crate::graph_support::{box_model, compound, first};

fn faces_of(model: &Model, edge: EdgeId) -> HashSet<FaceId> {
    model.ancestors_of::<Face>(edge).into_iter().collect()
}

#[test]
fn box_face_has_four_adjacent_faces() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let adjacent = model
        .related(face.into())
        .into_iter()
        .filter(|(r, _)| *r == Relation::AdjacentFace)
        .count();
    assert_eq!(adjacent, 4);
}

#[test]
fn box_edge_has_two_faces() {
    let (model, _) = box_model();
    assert_eq!(faces_of(&model, first::<Edge>(&model)).len(), 2);
}

#[test]
fn find_coedge_on_non_adjacent_face_is_none() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    let adjacent = faces_of(&model, edge);
    let (other, _) = model.table::<Face>().iter().find(|(f, _)| !adjacent.contains(f)).unwrap();
    assert_eq!(model.find_coedge(edge, other), None);
}

#[test]
fn find_coedge_skips_removed_coedges() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let wire = model.outer_wire(face).unwrap();
    let coedge = model[wire].coedges[0];
    let edge = model[coedge].edge;
    assert_eq!(model.find_coedge(edge, face), Some(coedge));
    model.remove_node(coedge);
    assert_ne!(model.find_coedge(edge, face), Some(coedge));
}

#[test]
fn outer_wire_skips_removed_wire() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let wire = model.outer_wire(face).unwrap();
    model.remove_node(wire);
    assert_ne!(model.outer_wire(face), Some(wire));
}

#[test]
fn boundary_and_manifold_follow_face_removal() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert!(!model.is_boundary_edge(edge));
    assert!(model.is_manifold_edge(edge));
    let face = *faces_of(&model, edge).iter().next().unwrap();
    model.remove_node(face);
    assert!(model.is_boundary_edge(edge));
    assert!(!model.is_manifold_edge(edge));
}

#[test]
fn removed_edge_is_neither_boundary_nor_manifold() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    model.remove_node(edge);
    assert!(!model.is_boundary_edge(edge));
    assert!(!model.is_manifold_edge(edge));
}

#[test]
fn sphere_seam_pair_is_symmetric() {
    let mut model = Model::new();
    Sphere::new(10.0).build(&mut model).unwrap();
    let (coedge, pair) = model
        .table::<CoEdge>()
        .iter()
        .find_map(|(c, _)| model.seam_pair(c).map(|p| (c, p)))
        .expect("no seam pair on a sphere");
    assert_eq!(model.seam_pair(pair), Some(coedge));
}

#[test]
fn seam_pair_of_out_of_range_coedge_is_none() {
    let (model, _) = box_model();
    assert_eq!(model.seam_pair(CoEdgeId::from_index(10_000)), None);
}

#[test]
fn merged_boxes_stay_valid_after_compact() {
    let mut model = Model::new();
    let a = crate::make_box(&mut model, 10.0, 20.0, 30.0);
    let b = crate::make_box(&mut model, 10.0, 20.0, 30.0);
    compound(&mut model, &[a, b]);
    model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    model.compact();
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn removing_a_solid_from_a_compound_then_compacting_stays_valid() {
    let (mut model, solid) = box_model();
    let c = compound(&mut model, &[solid]);
    let before = model.revision(c.node).unwrap().subtree;
    model.remove_node(solid.node);
    assert!(model.revision(c.node).unwrap().subtree > before);
    model.compact();
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn split_with_removed_vertex_fails() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, _) = model.edge_end_vertices(edge);
    model.remove_node(start);
    assert!(model.split_edge(edge, start, 0.5).is_err());
}

#[test]
fn replace_wire_edge_in_removed_wire_fails() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let wire = model.outer_wire(face).unwrap();
    let old = model[model[wire].coedges[0]].edge;
    let (a, b) = model.edge_end_vertices(old);
    let new = model.add_edge(a, b, 1.0e-7).unwrap();
    model.remove_node(wire);
    assert!(model.replace_wire_edge(wire, old, new, false).is_err());
}

#[test]
fn edge_has_two_live_vertex_uses() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    for vertex_use in [model[edge].start, model[edge].end] {
        let u = model.get(vertex_use).expect("removed vertex use");
        assert!(model.get(u.vertex).is_some());
    }
}

#[test]
fn removed_node_has_no_stable_lookup() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let sid = model.stable_id(face).unwrap();
    assert_eq!(model.find_stable::<Face>(sid), Some(face));
    model.remove_node(face);
    assert_eq!(model.find_stable::<Face>(sid), None);
    assert_eq!(model.find_stable_node(sid), None);
    // The slot keeps its stable id (for history), but the face is gone.
    assert_eq!(model.stable_id(face), Some(sid));
    assert_eq!(model.stamp(face), None);
}
