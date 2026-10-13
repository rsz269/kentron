//! BRepGraph_Iterator_Test.cxx (per-kind iteration is `model.table::<T>().iter()`)
//!
//! Not ported: `RootProductIterator_MatchesStoredRoots` (OCCT auto-creates a root product
//! per imported shape; kentron only has the products it is asked to build).

use kentron::topology::{CoEdge, Edge, Face, Model, Shell, Solid, Vertex, Wire};

use crate::graph_support::{box_model, first};

#[test]
fn face_iterator_count_matches_topology() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Face>().iter().count(), 6);
}

#[test]
fn edge_iterator_count_matches_topology() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Edge>().iter().count(), 12);
}

#[test]
fn vertex_iterator_count_matches_topology() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Vertex>().iter().count(), 8);
}

#[test]
fn solid_iterator_box_has_one_solid() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Solid>().iter().count(), 1);
}

#[test]
fn shell_iterator_box_has_one_shell() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Shell>().iter().count(), 1);
}

#[test]
fn wire_iterator_box_has_six_wires() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Wire>().iter().count(), 6);
}

#[test]
fn current_id_returns_valid_typed_ids() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert!(model.get(face).is_some());
}

#[test]
fn current_returns_definition() {
    let (model, _) = box_model();
    let (_, vertex) = model.table::<Vertex>().iter().next().unwrap();
    assert!(vertex.point.x.is_finite() && vertex.point.y.is_finite() && vertex.point.z.is_finite());
}

#[test]
fn removed_face_skipped_by_default_iterator() {
    let (mut model, _) = box_model();
    let before = model.table::<Face>().iter().count();
    model.remove_node(first::<Face>(&model));
    assert_eq!(model.table::<Face>().iter().count(), before - 1);
}

#[test]
fn full_traverse_includes_removed_face() {
    let (mut model, _) = box_model();
    let before = model.table::<Face>().iter().count();
    model.remove_node(first::<Face>(&model));
    assert_eq!(model.table::<Face>().iter_all().count(), before);
}

#[test]
fn range_for_works_correctly() {
    let (model, _) = box_model();
    let mut count = 0;
    for (_, _face) in model.table::<Face>().iter() {
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
fn empty_graph_iterator_is_empty() {
    let model = Model::new();
    assert!(model.table::<Face>().iter().next().is_none());
    assert!(model.table::<Edge>().iter().next().is_none());
    assert!(model.table::<Vertex>().iter().next().is_none());
    assert!(model.table::<Solid>().iter().next().is_none());
}

/// 6 faces x 4 edges = 24 coedges.
#[test]
fn co_edge_iterator_count_matches_topology() {
    let (model, _) = box_model();
    assert_eq!(model.table::<CoEdge>().iter().count(), 24);
}
