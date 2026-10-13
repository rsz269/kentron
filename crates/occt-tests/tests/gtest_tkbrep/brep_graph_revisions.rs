//! BRepGraph_VersionStamp_Test.cxx and BRepGraph_MutationGen_Test.cxx: revision counters
//! and version stamps.

use kentron::math::Transform;
use kentron::topology::{
    CoEdge, Edge, Face, Model, NodeId, NodeKind, Revision, Shell, Solid, Wire,
};

use crate::graph_support::{box_model, first};

// ---- VersionStamp ----

#[test]
fn stamp_of_valid_node_returns_valid_stamp() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let stamp = model.stamp(face).expect("no stamp for a live face");
    assert_eq!(Some(stamp.stable_id), model.stable_id(face));
    assert_eq!(stamp.revision, 0);
    assert_eq!(stamp.model, model.uid());
}

#[test]
fn stamp_of_invalid_node_returns_none() {
    let (model, _) = box_model();
    assert_eq!(model.stamp(NodeId { kind: NodeKind::Face, index: u32::MAX }), None);
}

#[test]
fn is_stale_unmutated_node_returns_false() {
    let (model, _) = box_model();
    let stamp = model.stamp(first::<Face>(&model)).unwrap();
    assert!(!model.is_stale(&stamp));
}

#[test]
fn is_stale_mutated_node_returns_true() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let stamp = model.stamp(face).unwrap();
    let tolerance = model.tolerance(face);
    model.set_face_tolerance(face, tolerance + 0.01).unwrap();
    assert!(model.is_stale(&stamp));
}

#[test]
fn is_stale_removed_node_returns_true() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let stamp = model.stamp(face).unwrap();
    model.remove_node(face);
    assert!(model.is_stale(&stamp));
}

#[test]
fn is_stale_other_model_returns_true() {
    let (model, _) = box_model();
    let stamp = model.stamp(first::<Face>(&model)).unwrap();
    // A rebuilt model is another model, even with identical content.
    let (rebuilt, _) = box_model();
    assert_ne!(rebuilt.uid(), model.uid());
    assert!(rebuilt.is_stale(&stamp));
}

#[test]
fn stamp_identity_same_node_equal() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert_eq!(model.stamp(face), model.stamp(face));
}

#[test]
fn stamp_identity_different_nodes_not_equal() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert_ne!(model.stamp(face), model.stamp(kentron::topology::Id::<Face>::from_index(1)));
}

#[test]
fn is_same_item_same_version_returns_true() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert!(model.stamp(face).unwrap().is_same_item(&model.stamp(face).unwrap()));
}

#[test]
fn is_same_item_different_version_still_same_item() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let before = model.stamp(face).unwrap();
    let tolerance = model.tolerance(face);
    model.set_face_tolerance(face, tolerance + 0.01).unwrap();
    let after = model.stamp(face).unwrap();
    assert_ne!(before, after);
    assert!(before.is_same_item(&after));
}

#[test]
fn stamp_of_assembly_nodes_works_for_products_and_instances() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    for node in [NodeId::from(part), assembly.into(), instance.into()] {
        let stamp = model.stamp(node).expect("no stamp");
        assert!(!model.is_stale(&stamp));
    }
}

// ---- MutationGen ----

#[test]
fn own_revision_incremented_on_mutation() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert_eq!(model.revision(edge), Some(Revision { own: 0, subtree: 0 }));
    model.set_edge_tolerance(edge, 0.5).unwrap();
    assert_eq!(model.revision(edge), Some(Revision { own: 1, subtree: 1 }));
}

#[test]
fn adding_an_instance_bumps_the_assembly() {
    let mut model = Model::new();
    let parent = model.new_assembly();
    let child = model.new_assembly();
    let parent_stamp = model.stamp(parent).unwrap();
    let parent_revision = model.revision(parent).unwrap().own;
    let instance = model.add_instance(parent, child, Transform::IDENTITY).unwrap();
    assert!(model.is_stale(&parent_stamp));
    assert_eq!(model.revision(parent).unwrap().own, parent_revision + 1);
    assert_eq!(model.revision(instance).unwrap().own, 0);
}

#[test]
fn own_revision_multiple_increments() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    model.set_edge_tolerance(edge, 0.1).unwrap();
    model.set_edge_tolerance(edge, 0.2).unwrap();
    assert_eq!(model.revision(edge).unwrap().own, 2);
}

fn revisions<T: kentron::topology::Node>(model: &Model) -> Vec<Revision> {
    model.table::<T>().iter().map(|(id, _)| model.revision(id).unwrap()).collect()
}

#[test]
fn subtree_revision_propagated_to_ancestors() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    model.set_edge_tolerance(edge, 0.5).unwrap();
    assert_eq!(model.revision(edge), Some(Revision { own: 1, subtree: 1 }));
    for (kind, revs) in [
        ("wire", revisions::<Wire>(&model)),
        ("face", revisions::<Face>(&model)),
        ("shell", revisions::<Shell>(&model)),
        ("solid", revisions::<Solid>(&model)),
    ] {
        assert!(revs.iter().any(|r| r.subtree > 0), "no {kind} saw the edit");
        assert!(revs.iter().all(|r| r.own == 0), "a {kind} counted the edit as its own");
    }
}

#[test]
fn subtree_revision_solid_edit_propagates_through_instances() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    let instance_before = model.revision(instance).unwrap().subtree;
    let assembly_before = model.revision(assembly).unwrap().subtree;
    let shell = model.new_empty(NodeKind::Shell);
    model.add(solid, shell).unwrap();
    assert!(model.revision(instance).unwrap().subtree > instance_before);
    assert!(model.revision(assembly).unwrap().subtree > assembly_before);
}

#[test]
fn clearing_a_surface_bumps_the_face() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    assert_eq!(model.revision(face), Some(Revision::default()));
    model.clear_face_surface(face.into()).unwrap();
    let r = model.revision(face).unwrap();
    assert!(r.own > 0 && r.subtree > 0);
}

#[test]
fn clearing_a_curve_bumps_the_edge() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert_eq!(model.revision(edge), Some(Revision::default()));
    model.clear_edge_curve(edge.into()).unwrap();
    let r = model.revision(edge).unwrap();
    assert!(r.own > 0 && r.subtree > 0);
}

#[test]
fn clearing_a_pcurve_bumps_the_coedge() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    assert_eq!(model.revision(coedge), Some(Revision::default()));
    model.clear_pcurve(coedge).unwrap();
    let r = model.revision(coedge).unwrap();
    assert!(r.own > 0 && r.subtree > 0);
}

#[test]
fn clearing_a_mesh_bumps_the_face() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    assert_eq!(model.revision(face), Some(Revision::default()));
    model.set_face_mesh(face.into(), None).unwrap();
    let r = model.revision(face).unwrap();
    assert!(r.own > 0 && r.subtree > 0);
}

#[test]
fn clearing_a_polyline_bumps_the_edge() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert_eq!(model.revision(edge), Some(Revision::default()));
    model.set_edge_polyline(edge.into(), None).unwrap();
    let r = model.revision(edge).unwrap();
    assert!(r.own > 0 && r.subtree > 0);
}
