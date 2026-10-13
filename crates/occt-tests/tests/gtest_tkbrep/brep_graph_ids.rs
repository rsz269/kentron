//! BRepGraph_NodeId_Test.cxx, BRepGraph_RefId_Test.cxx, BRepGraph_RepId_Test.cxx
//!
//! Not ported (C++ id mechanics with no Rust counterpart): `++`/`+`/`-` arithmetic,
//! hash equality across id types, out-of-range enum kinds, BRepGraph_ItemId and
//! BRepGraph_TypedIdDispatch; UID / version-stamp / GUID checks of the RefId test (ported
//! with BRepGraph_VersionStamp_Test).

use kentron::topology::{
    EdgeCurve, EdgeId, FaceId, FaceMesh, FaceSurface, FaceUseId, Id, NodeId, NodeKind, PCurve,
    Relation, Vertex, VertexId,
};

use crate::graph_support::box_model;

#[test]
fn construction_default_invalid() {
    let face = FaceId::default();
    assert!(!face.is_valid());
    assert_eq!(face, FaceId::INVALID);
}

#[test]
fn construction_from_index() {
    let face = FaceId::from_index(5);
    assert!(face.is_valid());
    assert_eq!(face.index(), 5);
}

#[test]
fn implicit_conversion_to_node_id() {
    let node = NodeId::from(FaceId::from_index(5));
    assert_eq!(node.kind, NodeKind::Face);
    assert_eq!(node.index, 5);
}

/// Typed ids work where node ids are expected: the faces adjacent to a box face.
#[test]
fn implicit_conversion_pass_to_function() {
    let (model, _) = box_model();
    let face = FaceId::from_index(0);
    let adjacent = model.related(face.into()).iter().filter(|(r, _)| *r == Relation::AdjacentFace).count();
    assert!(adjacent > 0);
}

#[test]
fn from_node_id_correct_kind() {
    let edge = EdgeId::try_from(NodeId { kind: NodeKind::Edge, index: 3 }).unwrap();
    assert_eq!(edge.index(), 3);
}

#[test]
fn comparison_typed_vs_typed() {
    assert_eq!(FaceId::from_index(3), FaceId::from_index(3));
    assert_ne!(FaceId::from_index(3), FaceId::from_index(5));
}

#[test]
fn comparison_typed_vs_node_id() {
    let face = FaceId::from_index(3);
    let matching = NodeId { kind: NodeKind::Face, index: 3 };
    assert_eq!(face, matching);
    assert_eq!(matching, face);
    assert_ne!(face, NodeId { kind: NodeKind::Edge, index: 3 });
    assert_ne!(face, NodeId { kind: NodeKind::Face, index: 5 });
}

#[test]
fn from_node_id_wrong_kind_returns_invalid() {
    assert!(EdgeId::try_from(NodeId { kind: NodeKind::Face, index: 3 }).is_err());
}

#[test]
fn out_of_range_metadata_queries_return_false() {
    let (model, _) = box_model();
    let out_of_range = VertexId::from_index(model.table::<Vertex>().capacity_len() as u32);
    assert!(model.get(out_of_range).is_none());
    assert!(model.stable_id(out_of_range).is_none());
}

// RefId_Test (use ids)

#[test]
fn default_ref_id_is_invalid() {
    assert!(!FaceUseId::default().is_valid());
}

// RepId_Test (geometry slot ids)

#[test]
fn typed_default_rep_id_is_invalid() {
    assert!(!Id::<EdgeCurve>::default().is_valid());
    assert!(!Id::<FaceSurface>::default().is_valid());
    assert!(!Id::<PCurve>::default().is_valid());
    assert!(!Id::<FaceMesh>::default().is_valid());
}

#[test]
fn start_id_is_zero() {
    let start = Id::<EdgeCurve>::from_index(0);
    assert!(start.is_valid());
    assert_eq!(start.index(), 0);
}

#[test]
fn equality_and_comparison() {
    let a1 = Id::<EdgeCurve>::from_index(1);
    let a2 = Id::<EdgeCurve>::from_index(2);
    let a1_copy = Id::<EdgeCurve>::from_index(1);
    assert_eq!(a1, a1_copy);
    assert_ne!(a1, a2);
    assert!(a1 < a2);
    assert!(a1 <= a1_copy);
    assert!(a2 > a1);
    assert!(a1_copy >= a1);
}

// ---- stable ids of use records (BRepGraph_RefId_Test.cxx, UID part) ----
//
// Not ported: `RefUID_Default_IsInvalid`, `RefUID_Equality_IgnoresGeneration` (no default
// or generation-tagged stable ids), stamps of use records (`MutFaceRef_*`,
// `RefDomain_Stamp*`: kentron stamps nodes only), and rebuild staleness (a rebuilt model is
// a new `Model`).

fn assert_use_stable_ids_round_trip<T: kentron::topology::Record>(model: &kentron::topology::Model) {
    let table = model.table::<T>();
    assert!(!table.is_empty());
    for (id, _) in table.iter() {
        let sid = model.stable_id(id).unwrap();
        assert_eq!(model.find_stable::<T>(sid), Some(id));
    }
}

#[test]
fn use_stable_ids_round_trip() {
    use kentron::topology::{FaceUse, ShellUse, VertexUse, WireUse};
    let (model, _) = crate::graph_support::box_model();
    assert_use_stable_ids_round_trip::<VertexUse>(&model);
    assert_use_stable_ids_round_trip::<WireUse>(&model);
    assert_use_stable_ids_round_trip::<FaceUse>(&model);
    assert_use_stable_ids_round_trip::<ShellUse>(&model);
}

#[test]
fn removed_use_has_no_stable_lookup() {
    use kentron::topology::{FaceUse, Shell};
    let (mut model, _) = crate::graph_support::box_model();
    let shell = crate::graph_support::first::<Shell>(&model);
    let face_use = model[shell].faces[0];
    let sid = model.stable_id(face_use).unwrap();
    model.remove_use(face_use);
    assert_eq!(model.find_stable::<FaceUse>(sid), None);
}

#[test]
fn stable_id_of_a_new_use_resolves() {
    use kentron::topology::{CompoundChild, NodeKind};
    let (mut model, solid) = crate::graph_support::box_model();
    let compound = model.new_empty(NodeKind::Compound);
    model.add(compound, solid).unwrap();
    let child = model[crate::graph_support::id_of::<kentron::topology::Compound>(compound)].children[0];
    let sid = model.stable_id(child).unwrap();
    assert_eq!(model.find_stable::<CompoundChild>(sid), Some(child));
}
