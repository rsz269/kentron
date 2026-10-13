//! BRepGraph_Sharing_Test.cxx

use kentron::math::{Transform, Vector};
use kentron::tolerance;
use kentron::topology::{
    CoEdge, Edge, Face, Model, NodeKind, Shell, Solid, Vertex, Wire,
};
use kentron::validity;

use crate::graph_support::{box_model, compound, id_of};
use crate::make_box;

/// In a box every edge is shared by exactly two faces.
#[test]
fn edge_def_each_shared_by_two_faces() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Edge), 12);
    for (edge, _) in model.table::<Edge>().iter() {
        assert_eq!(model.ancestors_of::<Face>(edge).len(), 2, "{edge:?}");
    }
}

#[test]
fn face_def_each_has_valid_surface() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Face), 6);
    for (face, record) in model.table::<Face>().iter() {
        assert!(model.get(record.surface).is_some(), "{face:?} has no surface");
    }
}

#[test]
fn solid_def_has_one_shell_ref() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Solid), 1);
    let (_, solid) = model.table::<Solid>().iter().next().unwrap();
    assert_eq!(solid.shells.len(), 1);
}

#[test]
fn shell_def_has_six_face_refs() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Shell), 1);
    let (_, shell) = model.table::<Shell>().iter().next().unwrap();
    assert_eq!(shell.faces.len(), 6);
}

#[test]
fn face_def_outer_wire_idx_valid() {
    let (model, _) = box_model();
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.outer_wire(face).is_some(), "{face:?} has no outer wire");
    }
}

#[test]
fn wire_def_co_edges_count_four_per_box_face() {
    let (model, _) = box_model();
    for (wire, record) in model.table::<Wire>().iter() {
        assert_eq!(record.coedges.len(), 4, "{wire:?}");
    }
}

#[test]
fn edge_def_vertex_defs_both_valid() {
    let (model, _) = box_model();
    for (edge, record) in model.table::<Edge>().iter() {
        assert!(model.get(record.start).is_some(), "{edge:?} start");
        assert!(model.get(record.end).is_some(), "{edge:?} end");
    }
}

/// Shared edges have coedges on different faces.
#[test]
fn shared_edge_incidence_refs_different_orientation() {
    let (model, _) = box_model();
    let shared = model
        .table::<Edge>()
        .iter()
        .filter(|(edge, _)| {
            let faces: Vec<_> = model
                .table::<CoEdge>()
                .iter()
                .filter(|(_, c)| c.edge == *edge)
                .map(|(_, c)| c.face)
                .collect();
            faces.windows(2).any(|w| w[0] != w[1])
        })
        .count();
    assert!(shared > 0, "expected edges with coedges on different faces");
}

/// Box edges are not closed: their end vertices differ.
#[test]
fn non_closed_edge_start_end_different() {
    let (model, _) = box_model();
    for (edge, _) in model.table::<Edge>().iter() {
        if model.is_degenerate(edge) {
            continue;
        }
        let (start, end) = model.edge_end_vertices(edge);
        assert_ne!(start, end, "{edge:?}");
    }
}

#[test]
fn vertex_def_points_match_expected_box_corners() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Vertex), 8);
    let tol = tolerance::LINEAR;
    for (_, vertex) in model.table::<Vertex>().iter() {
        let p = vertex.point;
        assert!(-tol <= p.x && p.x <= 10.0 + tol);
        assert!(-tol <= p.y && p.y <= 20.0 + tol);
        assert!(-tol <= p.z && p.z <= 30.0 + tol);
    }
}

/// The same solid used twice by a compound is stored once, with two child entries.
#[test]
fn compound_two_identical_boxes() {
    let (mut model, solid) = box_model();
    let c = compound(&mut model, &[solid, solid]);

    assert_eq!(model.count(NodeKind::Solid), 1);
    assert_eq!(model.count(NodeKind::Face), 6);
    assert_eq!(model.count(NodeKind::Edge), 12);
    assert_eq!(model.count(NodeKind::Vertex), 8);

    let children = &model[id_of::<kentron::topology::Compound>(c)].children;
    assert_eq!(children.len(), 2);
    assert_ne!(children[0], children[1], "each use has its own child entry");
    assert_eq!(model[children[0]].child, model[children[1]].child);
}

/// One solid used by two compounds: one child entry per parent.
#[test]
fn same_reference_representation_across_different_parents_uses_distinct_child_refs() {
    let (mut model, solid) = box_model();
    let parent_a = compound(&mut model, &[solid]);
    let parent_b = compound(&mut model, &[solid]);
    compound(&mut model, &[parent_a, parent_b]);

    let a_children = &model[id_of::<kentron::topology::Compound>(parent_a)].children;
    let b_children = &model[id_of::<kentron::topology::Compound>(parent_b)].children;
    assert_eq!(a_children.len(), 1);
    assert_eq!(b_children.len(), 1);
    assert_ne!(a_children[0], b_children[0]);
    assert_eq!(model[a_children[0]].child, model[b_children[0]].child);
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn compound_two_distinct_boxes() {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 20.0, 30.0);
    let box2 = make_box(&mut model, 5.0, 15.0, 25.0);
    compound(&mut model, &[box1, box2]);

    assert_eq!(model.count(NodeKind::Solid), 2);
    assert_eq!(model.count(NodeKind::Face), 12);
    assert_eq!(model.count(NodeKind::Edge), 24);
    assert_eq!(model.count(NodeKind::Vertex), 16);
}

/// A moved use of a solid shares its definition; the placement lives on the child entry.
fn assert_moved_use_shares_definition(offset: Vector) {
    let (mut model, solid) = box_model();
    let moved = solid.moved(&Transform::translation(offset));
    let c = compound(&mut model, &[solid, moved]);

    assert_eq!(model.count(NodeKind::Solid), 1);
    assert_eq!(model.count(NodeKind::Face), 6);
    assert_eq!(model.count(NodeKind::Edge), 12);
    assert_eq!(model.count(NodeKind::Vertex), 8);

    let children = &model[id_of::<kentron::topology::Compound>(c)].children;
    assert_eq!(children.len(), 2);
    assert!(!model[children[1]].placement.is_identity());
}

#[test]
fn compound_with_location_more_usages_than_defs() {
    assert_moved_use_shares_definition(Vector::new(100.0, 0.0, 0.0));
}

#[test]
fn translated_copy_same_t_shape_shared_defs() {
    assert_moved_use_shares_definition(Vector::new(50.0, 50.0, 50.0));
}

