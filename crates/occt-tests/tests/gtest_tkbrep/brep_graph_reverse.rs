//! BRepGraph_Reverse_Test.cxx: reversing and reordering wires and edges.
//!
//! OCCT's dry-run `Check*` queries are ported by trying the edit on a copy of the model.
//! Not ported: the `LayerTopoSupplement_*` and `Relations_SupplementLayerClearedOnGraphClear`
//! tests (supplement layer) and `ReversedShell_Reconstruct*` (TopoDS round trips).

use kentron::math::Point;
use kentron::topology::{
    CoEdgeId, Edge, EdgeId, Model, NodeKind, Orientation, Shape, VertexId, Wire, WireId,
};
use kentron::validity;

use crate::graph_support::{box_model, first, id_of};

fn segment(model: &mut Model, start: VertexId, end: VertexId) -> EdgeId {
    model.add_edge(start, end, 1.0e-7).unwrap()
}

/// Vertices at x = 0, 1, 2, ... on the X axis.
fn vertices(model: &mut Model, n: usize) -> Vec<VertexId> {
    (0..n).map(|i| model.add_vertex(Point::new(i as f64, 0.0, 0.0), 1.0e-7)).collect()
}

/// A wire of forward coedges along `edges`.
fn wire_of(model: &mut Model, edges: &[EdgeId]) -> (WireId, Vec<CoEdgeId>) {
    let coedges: Vec<CoEdgeId> = edges.iter().map(|&e| model.add_coedge(e, Orientation::Forward).unwrap()).collect();
    (model.add_wire(&coedges).unwrap(), coedges)
}

#[test]
fn wire_order_and_orientation_flipped() {
    let (mut model, _) = box_model();
    let wire = first::<Wire>(&model);
    let before = model[wire].coedges.clone();
    let orientations: Vec<Orientation> = before.iter().map(|&c| model[c].orientation).collect();
    assert!(before.len() > 1);

    model.reverse_wire(wire).unwrap();
    let after = &model[wire].coedges;
    let n = before.len();
    for (i, &coedge) in after.iter().enumerate() {
        assert_eq!(coedge, before[n - 1 - i], "order reversed at {i}");
        assert_eq!(model[coedge].orientation, orientations[n - 1 - i].reversed(), "flipped at {i}");
    }
}

#[test]
fn wire_order_set_co_edge_order_canonicalizes_unordered_permutation() {
    let (mut model, _) = box_model();
    let wire = first::<Wire>(&model);
    let c = model[wire].coedges.clone();
    assert!(c.len() >= 4);
    model.set_coedge_order(wire, &[c[0], c[2], c[1], c[3]]).unwrap();
    assert_eq!(model[wire].coedges.len(), c.len());
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn wire_order_check_co_edge_order_reports_current_and_reordered() {
    let (model, _) = box_model();
    let wire = first::<Wire>(&model);
    let c = model[wire].coedges.clone();

    let mut copy = model.clone();
    assert!(!copy.set_coedge_order(wire, &c).unwrap(), "the current order changes nothing");
    let mut copy = model.clone();
    assert!(copy.set_coedge_order(wire, &[c[0], c[2], c[1], c[3]]).unwrap(), "reordered");
    let mut copy = model.clone();
    assert!(copy.set_coedge_order(WireId::INVALID, &c).is_err(), "invalid wire");
}

#[test]
fn wire_order_set_co_edge_order_canonicalizes_open_wire_permutation() {
    let mut model = Model::new();
    let v = vertices(&mut model, 4);
    let edges = [segment(&mut model, v[0], v[1]), segment(&mut model, v[1], v[2]), segment(&mut model, v[2], v[3])];
    let (wire, c) = wire_of(&mut model, &edges);

    model.set_coedge_order(wire, &[c[2], c[0], c[1]]).unwrap();
    assert_eq!(model[wire].coedges, c, "an open chain starts at its free end");
}

#[test]
fn wire_order_check_append_co_edge_reports_ready_and_already_contained() {
    let mut model = Model::new();
    let v = vertices(&mut model, 3);
    let first_edge = segment(&mut model, v[0], v[1]);
    let (wire, initial) = wire_of(&mut model, &[first_edge]);
    let second_edge = segment(&mut model, v[1], v[2]);
    let second = model.add_coedge(second_edge, Orientation::Forward).unwrap();

    assert!(model.clone().append_coedge(wire, second).is_ok(), "ready");
    assert!(model.clone().append_coedge(wire, initial[0]).is_err(), "already contained");
}

#[test]
fn wire_order_check_replace_edge_reports_ready_and_disconnected() {
    let mut model = Model::new();
    let v = vertices(&mut model, 5);
    let old = segment(&mut model, v[0], v[1]);
    let next = segment(&mut model, v[1], v[2]);
    let (wire, _) = wire_of(&mut model, &[old, next]);
    let compatible = segment(&mut model, v[0], v[1]);
    let disconnected = segment(&mut model, v[3], v[4]);

    assert_eq!(model.clone().replace_wire_edge(wire, old, compatible, false).ok(), Some(true));
    assert_eq!(model.clone().replace_wire_edge(wire, old, old, false).ok(), Some(false));
    assert!(model.clone().replace_wire_edge(wire, old, disconnected, false).is_err());
}

#[test]
fn wire_order_remove_internal_co_edge_from_closed_wire_keeps_connected_order() {
    let (mut model, _) = box_model();
    let wire = first::<Wire>(&model);
    let before = model[wire].coedges.clone();
    model.remove_coedge(wire, before[1]).unwrap();
    assert_eq!(model[wire].coedges.len(), before.len() - 1);
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn edge_start_end_vertex_refs_swapped() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, end) = (model[edge].start, model[edge].end);
    assert_ne!(start, end);
    model.reverse_edge(edge).unwrap();
    assert_eq!((model[edge].start, model[edge].end), (end, start));
}

#[test]
fn edge_round_trip_restores_original() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, end) = (model[edge].start, model[edge].end);
    model.reverse_edge(edge).unwrap();
    model.reverse_edge(edge).unwrap();
    assert_eq!((model[edge].start, model[edge].end), (start, end));
}

/// A new shell made of a box's faces, each used with orientation `Reversed`.
fn reversed_face_shell() -> (Model, kentron::topology::ShellId) {
    let (mut model, solid) = box_model();
    let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
    let shell = model.new_empty(NodeKind::Shell);
    for face in faces {
        model.add(shell, Shape { orientation: Orientation::Reversed, ..face }).unwrap();
    }
    (model, id_of(shell))
}

#[test]
fn reversed_shell_graph_stores_face_ref_orientation() {
    let (model, shell) = reversed_face_shell();
    let uses = &model[shell].faces;
    assert_eq!(uses.len(), 6);
    for &u in uses {
        assert_eq!(model[u].orientation, Orientation::Reversed, "{:?}", model[u].face);
    }
}

#[test]
fn reversed_shell_parent_explorer_reports_shell() {
    let (model, shell) = reversed_face_shell();
    for &face_use in &model[shell].faces {
        let face = model[face_use].face;
        let mut explorer = model.explore_up(face).direct_only();
        let mut found = false;
        while let Some(parent) = explorer.next() {
            if parent.node == shell {
                found = true;
                assert_eq!(explorer.leaf_orientation(), Orientation::Reversed, "{face:?}");
            }
        }
        assert!(found, "{face:?} should have the new shell as a parent");
    }
}
