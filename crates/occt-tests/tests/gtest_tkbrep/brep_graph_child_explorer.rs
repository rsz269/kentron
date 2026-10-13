//! BRepGraph_ChildExplorer_Test.cxx (`Model::explore`, `Model::explore_all`,
//! `Model::children`)
//!
//! Not ported: `InvalidRoot_Empty` (no invalid shapes), `RootEqualsTarget_LinkKindNone`
//! and the parent/reference link kinds (explorer internals; direct children's use records
//! are read from the parent), `DirectChildren_ProductShapeRoot_ViaOccurrenceRef` (OCCT wraps
//! a part's root in an occurrence; a kentron part refers to it directly).

use std::collections::{HashMap, HashSet};

use occt_tests::assert_near;
use kentron::math::{Transform, Vector};
use kentron::tolerance;
use kentron::topology::{
    Face, Model, NodeKind, Orientation, Shape, Shell, Wire,
};

use crate::graph_support::{box_model, compound, first, id_of};
use crate::make_box;

#[test]
fn box_edge_occurrences_count24() {
    let (model, solid) = box_model();
    assert_eq!(model.explore(solid, NodeKind::Edge).count(), 24);
}

#[test]
fn box_face_occurrences_count6() {
    let (model, solid) = box_model();
    assert_eq!(model.explore(solid, NodeKind::Face).count(), 6);
}

#[test]
fn box_vertex_occurrences() {
    let (model, solid) = box_model();
    let vertices: Vec<Shape> = model.explore(solid, NodeKind::Vertex).collect();
    assert!(!vertices.is_empty());
    let distinct: HashSet<_> = vertices.iter().map(|v| v.node).collect();
    assert_eq!(distinct.len(), 8);
}

#[test]
fn face_edge_occurrences_4() {
    let (model, _) = box_model();
    assert_eq!(model.explore(first::<Face>(&model), NodeKind::Edge).count(), 4);
}

#[test]
fn root_equals_target_returns_self() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let found: Vec<Shape> = model.explore(face, NodeKind::Face).collect();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].node, face);
}

#[test]
fn avoid_kind_shell_skips_contained_faces() {
    let (model, solid) = box_model();
    assert_eq!(model.explore(solid, NodeKind::Face).skip_below(NodeKind::Shell).count(), 0);
}

#[test]
fn avoid_kind_emit_boundary_returns_faces_instead_of_edges() {
    let (model, solid) = box_model();
    let found: Vec<Shape> = model.explore(solid, NodeKind::Edge).stop_at(NodeKind::Face).collect();
    assert!(found.iter().all(|s| s.kind() == NodeKind::Face));
    assert_eq!(found.len(), 6);
}

#[test]
fn avoid_kind_same_as_target_is_ignored() {
    let (model, solid) = box_model();
    let found: Vec<Shape> = model.explore(solid, NodeKind::Face).skip_below(NodeKind::Face).collect();
    assert!(found.iter().all(|s| s.kind() == NodeKind::Face));
    assert_eq!(found.len(), 6);
}

#[test]
fn all_descendants_recursive_yields_all_kinds() {
    let (model, solid) = box_model();
    let mut counts: HashMap<NodeKind, usize> = HashMap::new();
    for s in model.explore_all(solid) {
        *counts.entry(s.kind()).or_default() += 1;
    }
    assert_eq!(counts.len(), 6, "unexpected kinds: {counts:?}");
    assert_eq!(counts[&NodeKind::Shell], 1);
    assert_eq!(counts[&NodeKind::Face], 6);
    assert_eq!(counts[&NodeKind::Wire], 6);
    assert_eq!(counts[&NodeKind::CoEdge], 24, "4 coedges per face");
    assert_eq!(counts[&NodeKind::Edge], 24, "one edge occurrence per coedge");
    assert_eq!(counts[&NodeKind::Vertex], 48, "2 vertices per edge occurrence");
}

#[test]
fn all_descendants_avoid_face_boundary_stops_below_faces() {
    let (model, solid) = box_model();
    let found: Vec<Shape> = model.explore_all(solid).stop_at(NodeKind::Face).collect();
    assert_eq!(found.iter().filter(|s| s.kind() == NodeKind::Shell).count(), 1);
    assert_eq!(found.iter().filter(|s| s.kind() == NodeKind::Face).count(), 6);
    assert_eq!(found.len(), 7);
}

#[test]
fn no_cum_loc_identity_location() {
    let (model, solid) = box_model();
    for edge in model.explore(solid, NodeKind::Edge).compose_placement(false) {
        assert!(edge.placement.is_identity());
    }
}

#[test]
fn no_cum_ori_forward_orientation() {
    let (model, solid) = box_model();
    for edge in model.explore(solid, NodeKind::Edge).compose_orientation(false) {
        assert_eq!(edge.orientation, Orientation::Forward);
    }
}

#[test]
fn global_location_box_identity() {
    let (model, solid) = box_model();
    for edge in model.explore(solid, NodeKind::Edge) {
        assert!(edge.placement.is_identity());
    }
}

#[test]
fn global_orientation_box_edges_forward_or_reversed() {
    let (model, solid) = box_model();
    for edge in model.explore(solid, NodeKind::Edge) {
        assert!(matches!(edge.orientation, Orientation::Forward | Orientation::Reversed));
    }
}

#[test]
fn compound_face_count() {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 10.0, 10.0);
    let box2 = make_box(&mut model, 20.0, 20.0, 20.0);
    let c = compound(&mut model, &[box1, box2]);
    assert_eq!(model.explore(c, NodeKind::Face).count(), 12);
}

#[test]
fn node_of_kind_face() {
    let (model, solid) = box_model();
    let mut explorer = model.explore(solid, NodeKind::Edge);
    explorer.next().expect("an edge");
    let face = explorer.ancestor(NodeKind::Face).expect("the edge's face");
    assert_eq!(face.kind(), NodeKind::Face);
}

#[test]
fn deep_compound_no_stack_overflow() {
    let mut model = Model::new();
    let mut inner = make_box(&mut model, 1.0, 1.0, 1.0);
    for _ in 0..100 {
        inner = compound(&mut model, &[inner]);
    }
    assert_eq!(model.explore(inner, NodeKind::Face).count(), 6);
}

#[test]
fn recreate_reset_and_reexplore() {
    let (model, solid) = box_model();
    assert_eq!(model.explore(solid, NodeKind::Face).count(), 6);
    assert_eq!(model.explore(solid, NodeKind::Edge).count(), 24);
}

#[test]
fn co_edge_target_reachable() {
    let (model, solid) = box_model();
    let found: Vec<Shape> = model.explore(solid, NodeKind::CoEdge).collect();
    assert!(found.iter().all(|s| s.kind() == NodeKind::CoEdge));
    assert_eq!(found.len(), 24);
}

#[test]
fn co_edge_target_from_face_count4() {
    let (model, _) = box_model();
    assert_eq!(model.explore(first::<Face>(&model), NodeKind::CoEdge).count(), 4);
}

/// Direct children come in the order of the parent's use records.
#[test]
fn direct_children_shell_faces_count_and_order() {
    let (model, _) = box_model();
    let shell = first::<Shell>(&model);
    let expected: Vec<_> = model[shell].faces.iter().map(|&u| model[u].face).collect();
    let actual: Vec<_> = model.children(shell).map(|s| id_of::<Face>(s)).collect();
    assert_eq!(actual, expected);
}

#[test]
fn direct_children_shell_faces_expose_parent_and_ref() {
    let (model, _) = box_model();
    let shell = first::<Shell>(&model);
    for (child, &face_use) in model.children(shell).zip(&model[shell].faces) {
        assert_eq!(model[face_use].shell, shell);
        assert_eq!(model[face_use].face, id_of::<Face>(child));
    }
}

#[test]
fn direct_children_product_occurrences_expose_occurrence_refs() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let i0 = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    let i1 = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();

    let uses = &model[assembly].children;
    assert_eq!(uses.len(), 2);
    assert_eq!([model[uses[0]].instance, model[uses[1]].instance], [i0, i1]);
    assert!(uses.iter().all(|&u| model[u].product == assembly));
    assert_eq!(model.children(assembly).count(), 2);
}

#[test]
fn direct_children_removed_face_ref_is_skipped() {
    let (mut model, _) = box_model();
    let shell = first::<Shell>(&model);
    let before = model[shell].faces.len();
    let removed_use = model[shell].faces[0];
    let removed_face = model[removed_use].face;
    model.remove_use(removed_use);

    let children: Vec<_> = model.children(shell).collect();
    assert!(children.iter().all(|c| id_of::<Face>(*c) != removed_face));
    assert_eq!(children.len(), before - 1);
}

#[test]
fn direct_children_wire_children_are_co_edges() {
    let (model, _) = box_model();
    let wire = first::<Wire>(&model);
    let children: Vec<_> = model.children(wire).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::CoEdge));
    assert_eq!(children.len(), model[wire].coedges.len());
}

#[test]
fn direct_children_compound_children_basic() {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 10.0, 10.0);
    let box2 = make_box(&mut model, 20.0, 20.0, 20.0);
    let c = compound(&mut model, &[box1, box2]);
    assert_eq!(model.children(c).filter(|s| s.kind() == NodeKind::Solid).count(), 2);
}

/// Walking down level by level, composing as we go, matches the recursive walk.
#[test]
fn direct_children_chained_traversal_parity_with_recursive_traversal() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0)
        .moved(&Transform::translation(Vector::new(100.0, 0.0, 0.0)));
    let c = compound(&mut model, &[solid]);

    let mut expected: HashMap<_, Shape> = HashMap::new();
    for face in model.explore(c, NodeKind::Face) {
        expected.entry(face.node).or_insert(face);
    }

    let mut visited = 0;
    for solid in model.children(c).filter(|s| s.kind() == NodeKind::Solid) {
        for shell in model.children(solid) {
            for face in model.children(shell) {
                let e = expected[&face.node];
                assert_eq!(face.placement, e.placement);
                assert_eq!(face.orientation, e.orientation);
                visited += 1;
            }
        }
    }
    assert_eq!(visited, expected.len());
}

/// One part instanced twice: its solid is reached twice, at different placements.
#[test]
fn recursive_shared_product_children_have_distinct_contexts() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::translation(Vector::new(10.0, 0.0, 0.0))).unwrap();
    model.add_instance(assembly, part, Transform::translation(Vector::new(20.0, 0.0, 0.0))).unwrap();

    let solids: Vec<Shape> = model.explore(assembly, NodeKind::Solid).collect();
    assert_eq!(solids.len(), 2);
    assert!(solids.iter().all(|s| s.node == solid.node));
    assert_ne!(solids[0].placement, solids[1].placement);
}

#[test]
fn recursive_product_occurrence_chain_reaches_nested_occurrences() {
    let mut model = Model::new();
    let part = model.new_assembly();
    let sub = model.new_assembly();
    let top = model.new_assembly();
    let sub_instance = model.add_instance(top, sub, Transform::IDENTITY).unwrap();
    let part_instance = model.add_instance(sub, part, Transform::IDENTITY).unwrap();

    let instances: Vec<_> = model.explore(top, NodeKind::Instance).map(|s| s.node).collect();
    assert_eq!(instances.len(), 2);
    assert!(instances.contains(&sub_instance.into()));
    assert!(instances.contains(&part_instance.into()));
}

/// The placement of a part's root composes with the instance placement.
#[test]
fn recursive_product_part_root_context_composed_with_occurrence() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid.with_placement(Transform::translation(Vector::new(0.0, 20.0, 0.0)))).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::translation(Vector::new(10.0, 0.0, 0.0))).unwrap();

    let found: Vec<Shape> = model.explore(assembly, NodeKind::Solid).collect();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].node, solid.node);
    let t = found[0].placement.translation_part();
    assert_near!(t.x, 10.0, tolerance::LINEAR);
    assert_near!(t.y, 20.0, tolerance::LINEAR);
    assert_near!(t.z, 0.0, tolerance::LINEAR);
}

#[test]
fn direct_children_high_fanout_direct_children_complete() {
    let mut model = Model::new();
    let boxes: Vec<Shape> = (0..80).map(|i| make_box(&mut model, 1.0 + i as f64, 2.0, 3.0)).collect();
    let c = compound(&mut model, &boxes);
    assert_eq!(model.children(c).count(), 80);
}

#[test]
fn high_fanout_completes_all_children() {
    let mut model = Model::new();
    let boxes: Vec<Shape> = (0..500).map(|_| make_box(&mut model, 1.0, 1.0, 1.0)).collect();
    let c = compound(&mut model, &boxes);
    assert_eq!(model.explore(c, NodeKind::Face).count(), 500 * 6);
}

