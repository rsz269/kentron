//! BRepGraph_LayerHistory_Test.cxx: the modification history.
//!
//! OCCT creates the history layer on first use (`LayerRegistry().Ensure<...>()`); in kentron
//! it is opted into with `Model::enable_history`, which the fixture here does.

use kentron::math::Point;
use kentron::topology::{ChangeKind, Edge, EdgeId, History, Id, NodeId, NodeKind, StableId, Wire};

use crate::graph_support;

/// The box fixture with history recording turned on.
fn box_model() -> (kentron::topology::Model, kentron::topology::Shape) {
    let (mut model, solid) = graph_support::box_model();
    model.enable_history();
    (model, solid)
}

fn node(kind: NodeKind, index: u32) -> NodeId {
    NodeId { kind, index }
}

fn edge(index: u32) -> NodeId {
    node(NodeKind::Edge, index)
}

/// The id the next edge would get.
fn next_edge(model: &kentron::topology::Model) -> NodeId {
    edge(model.slot_count(NodeKind::Edge) as u32)
}

// ---- derivation chains ----

/// edge0 -> edge1 -> edge2 through two modifications that each name a fresh edge id.
fn chain(model: &mut kentron::topology::Model, kind: NodeKind) -> (NodeId, NodeId, NodeId) {
    let a = node(kind, 0);
    let mut b = a;
    model.apply_modification(a, "Step1", |m, _| {
        b = node(kind, m.slot_count(kind) as u32);
        vec![b]
    });
    let mut c = b;
    model.apply_modification(b, "Step2", |m, _| {
        c = node(kind, m.slot_count(kind) as u32 + 1);
        vec![c]
    });
    (a, b, c)
}

#[test]
fn find_original_chain_abc_returns_a() {
    let (mut model, _) = box_model();
    let (a, _, c) = chain(&mut model, NodeKind::Edge);
    assert_eq!(model.history().unwrap().find_original(c), a);
}

#[test]
fn find_derived_chain_abc_contains_b_and_c() {
    let (mut model, _) = box_model();
    let (a, b, c) = chain(&mut model, NodeKind::Edge);
    let derived = model.history().unwrap().find_derived(a);
    assert!(derived.contains(&b));
    assert!(derived.contains(&c));
}

#[test]
fn find_original_unmodified_node_returns_self() {
    let (model, _) = box_model();
    let face = node(NodeKind::Face, 0);
    assert_eq!(model.history().unwrap().find_original(face), face);
}

#[test]
fn find_derived_unmodified_node_returns_empty() {
    let (model, _) = box_model();
    assert!(model.history().unwrap().find_derived(node(NodeKind::Face, 0)).is_empty());
}

#[test]
fn find_original_two_apply_transitive_trace() {
    let (mut model, _) = box_model();
    let (a, b, c) = chain(&mut model, NodeKind::Vertex);
    assert_eq!(model.history().unwrap().find_original(c), a);
    assert_eq!(model.history().unwrap().find_original(b), a);
}

// ---- enabling ----

#[test]
fn disabled_record_history_no_record_stored() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().set_enabled(false);
    model.history_mut().unwrap().record("Disabled", edge(0), &[edge(1)]);
    assert_eq!(model.history().unwrap().len(), before);
}

#[test]
fn disabled_apply_modification_modifier_still_runs() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().set_enabled(false);
    let mut ran = false;
    model.apply_modification(edge(0), "DisabledOp", |_, _| {
        ran = true;
        Vec::new()
    });
    assert!(ran);
}

#[test]
fn re_enabled_records_after_re_enable() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().set_enabled(false);
    assert!(!model.history().unwrap().is_enabled());
    model.history_mut().unwrap().set_enabled(true);
    assert!(model.history().unwrap().is_enabled());
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().record("ReEnabled", edge(0), &[edge(1)]);
    assert_eq!(model.history().unwrap().len(), before + 1);
}

// ---- recording ----

#[test]
fn apply_modification_empty_replacements() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.apply_modification(edge(0), "Delete", |_, _| Vec::new());
    assert_eq!(model.history().unwrap().len(), before + 1);
    assert_eq!(model.history().unwrap().events().last().unwrap().operation, "Delete");
}

#[test]
fn apply_modification_multiple_replacements() {
    let (mut model, _) = box_model();
    let (mut n1, mut n2) = (edge(0), edge(0));
    model.apply_modification(edge(0), "Split", |m, _| {
        n1 = next_edge(m);
        n2 = edge(n1.index + 1);
        vec![n1, n2]
    });
    let derived = model.history().unwrap().find_derived(edge(0));
    assert!(derived.contains(&n1));
    assert!(derived.contains(&n2));
}

#[test]
fn record_history_empty_replacements_stored() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().record("Erase", edge(0), &[]);
    assert_eq!(model.history().unwrap().len(), before + 1);
}

#[test]
fn history_record_sequence_number_monotonic() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().record("Op1", edge(0), &[edge(1)]);
    model.history_mut().unwrap().record("Op2", edge(1), &[edge(2)]);
    let n = model.history().unwrap().len();
    assert!(n >= 2);
    assert!(model.history().unwrap().event(n - 2).sequence < model.history().unwrap().event(n - 1).sequence);
}

#[test]
fn history_record_operation_name_stored() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().record("MyCustomOp", edge(0), &[edge(1)]);
    assert_eq!(model.history().unwrap().events().last().unwrap().operation, "MyCustomOp");
}

#[test]
fn nb_history_records_after_multiple_ops_correct() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    for op in ["A", "B", "C"] {
        model.history_mut().unwrap().record(op, edge(0), &[edge(1)]);
    }
    assert_eq!(model.history().unwrap().len(), before + 3);
}

#[test]
fn apply_modification_when_modifier_panics_does_not_record_history() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    occt_tests::assert_panics!(model.apply_modification(edge(0), "ThrowingModification", |_, _| {
        panic!("Synthetic failure")
    }));
    assert_eq!(model.history().unwrap().len(), before);
}

// ---- edge split ----

fn wires_of_edge(model: &kentron::topology::Model, edge: EdgeId) -> Vec<Id<Wire>> {
    model.ancestors_of::<Wire>(edge)
}

#[test]
fn split_edge_rewrites_all_containing_wires() {
    let (mut model, _) = box_model();
    let old = EdgeId::from_index(0);
    let (first, last) = model.edge_range(old);
    let vertices_before = model.count(NodeKind::Vertex);
    let vertex = model.add_vertex(Point::new(1.0, 2.0, 3.0), 1.0e-7);
    assert_eq!(model.count(NodeKind::Vertex), vertices_before + 1);
    let wires = wires_of_edge(&model, old);
    assert!(!wires.is_empty());
    let slots_before = model.slot_count(NodeKind::Edge);
    let active_before = model.count(NodeKind::Edge);
    let (a, b) = model.split_edge(old, vertex, 0.5 * (first + last)).unwrap();
    assert_eq!(model.slot_count(NodeKind::Edge), slots_before + 2);
    assert_eq!(model.count(NodeKind::Edge), active_before + 1);
    for wire in wires {
        let coedges = &model[wire].coedges;
        let position = |e: EdgeId| coedges.iter().position(|&c| model[c].edge == e);
        assert_eq!(position(old), None);
        let pa = position(a).expect("first half missing");
        let pb = position(b).expect("second half missing");
        let distance = pa.abs_diff(pb);
        assert!(distance == 1 || distance + 1 == coedges.len());
    }
}

#[test]
fn split_edge_ignores_removed_coedge_entries() {
    let (mut model, _) = box_model();
    let old = EdgeId::from_index(0);
    let (first, last) = model.edge_range(old);
    let wires = wires_of_edge(&model, old);
    assert!(wires.len() > 1);
    let wire = wires[0];
    let removed = *model[wire].coedges.iter().find(|&&c| model[c].edge == old).unwrap();
    model.remove_node(removed);
    assert!(model.table::<kentron::topology::CoEdge>().is_removed(removed));
    let wire_len = model[wire].coedges.len();
    let vertex = model.add_vertex(Point::new(4.0, 5.0, 6.0), 1.0e-7);
    let (a, b) = model.split_edge(old, vertex, 0.5 * (first + last)).unwrap();
    assert!(model.table::<kentron::topology::CoEdge>().is_removed(removed));
    assert_eq!(model[wire].coedges.len(), wire_len);
    // The removed coedge's record is left as it was.
    let (record, _) = model
        .table::<kentron::topology::CoEdge>()
        .iter_all()
        .find(|(id, _, _)| *id == removed)
        .map(|(_, c, r)| (c.edge, r))
        .unwrap();
    assert_eq!(record, old);
    let coedges = model.table::<kentron::topology::CoEdge>();
    assert!(coedges.iter().any(|(_, c)| c.edge == a));
    assert!(coedges.iter().any(|(_, c)| c.edge == b));
}

#[test]
fn apply_modification_split_edge_records_both_derived_nodes() {
    let (mut model, _) = box_model();
    let old = EdgeId::from_index(0);
    let (first, last) = model.edge_range(old);
    let vertex = model.add_vertex(Point::new(4.0, 5.0, 6.0), 1.0e-7);
    let before = model.history().unwrap().len();
    model.apply_modification(old, "Split", |m, target| {
        let target = Id::<Edge>::try_from(target).unwrap();
        let (a, b) = m.split_edge(target, vertex, 0.5 * (first + last)).unwrap();
        vec![a.into(), b.into()]
    });
    assert_eq!(model.history().unwrap().len(), before + 1);
    assert!(model.history().unwrap().find_derived(old.into()).len() >= 2);
}

// ---- per-kind lookup ----

#[test]
fn find_modified_per_kind_lookup_returns_only_modified() {
    let (mut model, _) = box_model();
    let mut image = edge(0);
    model.apply_modification(edge(0), "Modify", |m, _| {
        image = next_edge(m);
        vec![image]
    });
    assert_eq!(model.history().unwrap().modified(edge(0)), Some(&[image][..]));
    assert_eq!(model.history().unwrap().generated(edge(0)), None);
}

#[test]
fn find_generated_per_kind_lookup_returns_only_generated() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().record_as("Generate", edge(0), &[edge(1)], ChangeKind::Generated);
    assert_eq!(model.history().unwrap().generated(edge(0)), Some(&[edge(1)][..]));
    assert_eq!(model.history().unwrap().modified(edge(0)), None);
}

#[test]
fn record_deleted_marks_nodes_as_deleted() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().record_deleted("DeleteEdges", &[edge(0), edge(1)]);
    assert_eq!(model.history().unwrap().len(), before + 1);
    assert!(model.history().unwrap().is_deleted(edge(0)));
    assert!(model.history().unwrap().is_deleted(edge(1)));
    assert!(model.history().unwrap().deleted().contains(&edge(0)));
    assert!(model.history().unwrap().deleted().contains(&edge(1)));
}

#[test]
fn record_deleted_empty_list_is_noop() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().record_deleted("Noop", &[]);
    assert_eq!(model.history().unwrap().len(), before);
}

#[test]
fn record_empty_replacements_auto_deleted() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().record_as("Consumed", edge(0), &[], ChangeKind::Generated);
    assert!(model.history().unwrap().is_deleted(edge(0)));
    assert_eq!(model.history().unwrap().len(), 1);
    assert_eq!(model.history().unwrap().event(0).kind, ChangeKind::Deleted);
}

#[test]
fn find_derived_modified_and_generated_collects_both() {
    let (mut model, _) = box_model();
    let mut modified = edge(0);
    model.apply_modification(edge(0), "FirstModified", |m, _| {
        modified = next_edge(m);
        vec![modified]
    });
    let generated = edge(2);
    model.history_mut().unwrap().record_as("ThenGenerated", modified, &[generated], ChangeKind::Generated);
    let all = model.history().unwrap().find_derived(edge(0));
    assert!(!all.is_empty());
    assert!(all.contains(&generated));
    assert!(model.history().unwrap().generated(modified).is_some());
    assert!(model.history().unwrap().modified(edge(0)).is_some());
}

#[test]
fn record_batch_explicit_generated_kind_stored_correctly() {
    let (mut model, _) = box_model();
    model.history_mut().unwrap().record_pairs("BatchGen", &[edge(0)], &[edge(1)], ChangeKind::Generated);
    assert_eq!(model.history().unwrap().len(), 1);
    assert_eq!(model.history().unwrap().event(0).kind, ChangeKind::Generated);
    assert!(model.history().unwrap().generated(edge(0)).is_some());
    assert_eq!(model.history().unwrap().modified(edge(0)), None);
}

#[test]
fn record_explicit_modified_kind_default_matches() {
    let (mut model, _) = box_model();
    let before = model.history().unwrap().len();
    model.history_mut().unwrap().record_as("ExplicitMod", edge(0), &[edge(1)], ChangeKind::Modified);
    assert_eq!(model.history().unwrap().len(), before + 1);
    assert_eq!(model.history().unwrap().event(before).kind, ChangeKind::Modified);
}

#[test]
fn is_deleted_unrelated_node_returns_false() {
    let (model, _) = box_model();
    assert!(!model.history().unwrap().is_deleted(edge(0)));
}

#[test]
fn record_replaced_maps_image_and_marks_deleted() {
    let (mut model, _) = box_model();
    let original = node(NodeKind::Face, 0);
    let replacement = node(NodeKind::Face, 1);
    model.history_mut().unwrap().record_replaced("ReplaceFace", original, replacement);
    assert_eq!(model.history().unwrap().len(), 1);
    assert_eq!(model.history().unwrap().event(0).kind, ChangeKind::Replaced);
    assert!(model.history().unwrap().is_deleted(original));
    assert_eq!(model.history().unwrap().modified(original), Some(&[replacement][..]));
    assert_eq!(model.history().unwrap().originals(replacement), Some(&[original][..]));
}

#[test]
fn structural_replacement_does_not_emit_semantic_history() {
    let (mut model, _) = box_model();
    assert!(model.count(NodeKind::Face) > 1);
    model.history_mut().unwrap().clear();
    let removed = node(NodeKind::Face, 1);
    model.replace_node(removed, node(NodeKind::Face, 0)).unwrap();
    assert_eq!(model.history().unwrap().len(), 0);
    assert!(!model.history().unwrap().is_deleted(removed));
    assert_eq!(model.history().unwrap().modified(removed), None);
}

#[test]
fn find_originals_derived_with_two_parents_returns_both() {
    let (mut model, _) = box_model();
    let (a, b) = (node(NodeKind::Face, 0), node(NodeKind::Face, 1));
    let derived = edge(0);
    model.history_mut().unwrap().record("FirstParent", a, &[derived]);
    model.history_mut().unwrap().record("SecondParent", b, &[derived]);
    let originals = model.history().unwrap().originals(derived).unwrap();
    assert_eq!(originals.len(), 2);
    assert!(originals.contains(&a));
    assert!(originals.contains(&b));
}

// ---- histories keyed by stable id ----

#[test]
fn record_stable_ids_modified_appends_audit_record() {
    let (original, image) = (StableId(10), StableId(11));
    let mut history = History::<StableId>::new();
    history.record_as("UidMod", original, &[image], ChangeKind::Modified);
    assert_eq!(history.len(), 1);
    let event = history.event(0);
    assert_eq!(event.operation, "UidMod");
    assert_eq!(event.kind, ChangeKind::Modified);
    assert_eq!(event.mapping, vec![(original, vec![image])]);
    assert!(history.is_known_input(original));
    assert_eq!(history.modified(original), Some(&[image][..]));
}

#[test]
fn record_deleted_stable_ids_appends_audit_record() {
    let (a, b) = (StableId(20), StableId(21));
    let mut history = History::<StableId>::new();
    history.record_deleted("UidDelete", &[a, b]);
    assert_eq!(history.len(), 1);
    let event = history.event(0);
    assert_eq!(event.operation, "UidDelete");
    assert_eq!(event.kind, ChangeKind::Deleted);
    assert_eq!(event.mapping, vec![(a, vec![]), (b, vec![])]);
    assert!(history.is_deleted(a) && history.is_deleted(b));
    assert!(history.is_known_input(a) && history.is_known_input(b));
}

#[test]
fn clear_releases_history_and_allows_fresh_recording() {
    let mut history = History::<StableId>::new();
    for i in 0..1000 {
        history.record("Stress", StableId(i + 1), &[StableId(i + 1001)]);
    }
    assert_eq!(history.len(), 1000);
    history.clear();
    assert_eq!(history.len(), 0);
    let (original, image) = (StableId(5020), StableId(5021));
    history.record("AfterClear", original, &[image]);
    assert_eq!(history.len(), 1);
    assert_eq!(history.modified(original), Some(&[image][..]));
    assert!(!history.is_deleted(original));
}

// ---- opting in ----

#[test]
fn new_model_records_no_history() {
    let (mut model, _) = graph_support::box_model();
    assert!(model.history().is_none());
    let mut ran = false;
    model.apply_modification(edge(0), "NotRecorded", |_, _| {
        ran = true;
        vec![edge(1)]
    });
    assert!(ran);
    assert!(model.history().is_none());
}

#[test]
fn enable_history_keeps_existing_events() {
    let mut model = kentron::topology::Model::new();
    model.enable_history().record("A", edge(0), &[edge(1)]);
    assert_eq!(model.enable_history().len(), 1);
}

#[test]
fn disable_history_returns_recorded_events() {
    let (mut model, _) = box_model();
    model.apply_modification(edge(0), "Recorded", |_, _| vec![edge(1)]);
    let history = model.disable_history().expect("history was enabled");
    assert_eq!(history.len(), 1);
    assert!(model.history().is_none());
    model.apply_modification(edge(1), "NotRecorded", |_, _| vec![edge(2)]);
    assert!(model.history().is_none());
}
