//! BRepGraph_Test.cxx: the general BRepGraph suite. Most of it repeats cases ported in the
//! more specific modules (build counts, history, revisions, stable ids, removal); this
//! module holds the rest: topology analyses written on top of the query API, wire edge
//! replacement and use removal.
//!
//! Not ported: reconstruction (`Reconstruct*`, `Shape_*`, `HasOriginalShape_*`),
//! allocator and parallel-build tests, UID generation across rebuilds, and the
//! "prunes orphan" removal tests (`Remove{Face,Shell,Wire,CoEdge}_PrunesOrphan*`,
//! `RemoveShapeRoot_*`): in kentron detaching a child never removes it implicitly;
//! `AddPCurve_*` (a coedge on a face outside any wire) and
//! `DetectToleranceConflicts_ManualConflict_Detected` (vacuous on a box: it shares no
//! curve).

use std::collections::{HashMap, HashSet};

use kentron::topology::{
    Edge, EdgeId, Face, FaceId, Model, NodeKind, Shape, Shell, Wire, WireId,
};

use crate::graph_support::{box_model, compound, first, id_of};

fn faces_of_edge(model: &Model, edge: EdgeId) -> usize {
    model.ancestors_of::<Face>(edge).into_iter().collect::<HashSet<_>>().len()
}

fn free_edges(model: &Model) -> Vec<EdgeId> {
    model
        .table::<Edge>()
        .iter()
        .filter(|&(e, _)| !model.is_degenerate(e) && faces_of_edge(model, e) == 1)
        .map(|(e, _)| e)
        .collect()
}

/// Groups of faces connected through shared edges.
fn face_components(model: &Model) -> usize {
    let faces: Vec<FaceId> = model.table::<Face>().iter().map(|(f, _)| f).collect();
    let mut component: HashMap<FaceId, usize> = HashMap::new();
    let mut count = 0;
    for &start in &faces {
        if component.contains_key(&start) {
            continue;
        }
        let mut stack = vec![start];
        component.insert(start, count);
        while let Some(face) = stack.pop() {
            for edge in model.descendants_of::<Edge>(face) {
                for next in model.ancestors_of::<Face>(edge) {
                    if component.insert(next, count).is_none() {
                        stack.push(next);
                    }
                }
            }
        }
        count += 1;
    }
    count
}

/// A model with the first face of each of the given boxes, in one compound.
fn first_faces_of_boxes(sizes: &[f64]) -> Model {
    let mut source = Model::new();
    let faces: Vec<Shape> = sizes
        .iter()
        .map(|&s| {
            let solid = crate::make_box(&mut source, s, s, s);
            crate::first(&source, solid, NodeKind::Face)
        })
        .collect();
    let mut model = Model::new();
    let copies: Vec<Shape> = faces.iter().map(|&f| model.import(&source, f).unwrap()).collect();
    compound(&mut model, &copies);
    model
}

/// A copy of `edge` between the same vertices (reversed if `reversed`).
fn compatible_replacement(model: &mut Model, edge: EdgeId, reversed: bool) -> EdgeId {
    let (start, end) = model.edge_end_vertices(edge);
    let (a, b) = if reversed { (end, start) } else { (start, end) };
    let tolerance = model.tolerance(edge);
    let copy = model.add_edge(a, b, tolerance).unwrap();
    if let Some((curve, first, last)) = model.edge_curve(edge) {
        model.set_edge_curve(copy.into(), curve, tolerance).unwrap();
        model.set_edge_range(copy.into(), first, last).unwrap();
    }
    copy
}

fn first_wire(model: &Model) -> WireId {
    first::<Wire>(model)
}

#[test]
fn box_counts() {
    let (model, _) = box_model();
    let expected = [
        (NodeKind::Solid, 1),
        (NodeKind::Shell, 1),
        (NodeKind::Face, 6),
        (NodeKind::Wire, 6),
        (NodeKind::Edge, 12),
        (NodeKind::Vertex, 8),
    ];
    for (kind, n) in expected {
        assert_eq!(model.count(kind), n, "{kind:?}");
    }
    let total: usize = NodeKind::ALL.iter().map(|&k| model.count(k)).sum();
    assert_eq!(total, 1 + 1 + 6 + 6 + 24 + 12 + 8);
}

#[test]
fn box_edges_have_two_faces() {
    let (model, _) = box_model();
    for (edge, _) in model.table::<Edge>().iter() {
        assert_eq!(faces_of_edge(&model, edge), 2, "{edge:?}");
    }
    assert!(free_edges(&model).is_empty());
}

#[test]
fn single_face_edges_are_free() {
    let model = first_faces_of_boxes(&[10.0]);
    assert_eq!(free_edges(&model).len(), 4);
}

#[test]
fn separate_faces_are_separate_components() {
    let two = first_faces_of_boxes(&[10.0, 20.0]);
    assert_eq!(two.count(NodeKind::Face), 2);
    assert_eq!(face_components(&two), 2);
    let three = first_faces_of_boxes(&[10.0, 20.0, 30.0]);
    assert_eq!(face_components(&three), 3);
    let (solid, _) = box_model();
    assert_eq!(face_components(&solid), 1);
}

#[test]
fn box_has_no_same_domain_faces() {
    let (model, _) = box_model();
    let faces: Vec<FaceId> = model.table::<Face>().iter().map(|(f, _)| f).collect();
    for (i, &a) in faces.iter().enumerate() {
        for &b in &faces[i + 1..] {
            assert_ne!(model.surface(a), model.surface(b), "{a:?} and {b:?} share a surface");
        }
    }
}

#[test]
fn box_has_no_missing_pcurves() {
    let (model, _) = box_model();
    for (face, _) in model.table::<Face>().iter() {
        for edge in model.descendants_of::<Edge>(face) {
            if model.is_degenerate(edge) {
                continue;
            }
            let coedge = model.find_coedge(edge, face).expect("edge without coedge on its face");
            assert!(model.pcurve(coedge).is_some());
        }
    }
}

#[test]
fn box_has_no_degenerate_wires() {
    let (model, _) = box_model();
    for (wire, w) in model.table::<Wire>().iter() {
        assert!(w.coedges.len() >= 2);
        if let Some(face) = model.face_of_wire(wire) {
            if model.outer_wire(face) == Some(wire) {
                assert!(model.is_closed_wire(wire));
            }
        }
        assert_eq!(w.coedges.len(), 4);
    }
}

#[test]
fn replace_wire_edge_substitutes() {
    let (mut model, _) = box_model();
    let wire = first_wire(&model);
    let coedge = model[wire].coedges[0];
    let old = model[coedge].edge;
    let new = compatible_replacement(&mut model, old, false);
    assert!(model.replace_wire_edge(wire, old, new, false).unwrap());
    let coedge = model[wire].coedges[0];
    assert_eq!(model[coedge].edge, new);
    // The old edge no longer has a coedge in this wire; the new one does.
    assert!(!model.ancestors_of::<Wire>(old).contains(&wire));
    assert!(model.ancestors_of::<Wire>(new).contains(&wire));
}

#[test]
fn replace_wire_edge_reversed_flips_orientation() {
    let (mut model, _) = box_model();
    let wire = first_wire(&model);
    let coedge = model[wire].coedges[0];
    let old = model[coedge].edge;
    let orientation = model[coedge].orientation;
    let new = compatible_replacement(&mut model, old, true);
    assert!(model.replace_wire_edge(wire, old, new, true).unwrap());
    let coedge = model[wire].coedges[0];
    assert_eq!(model[coedge].edge, new);
    assert_eq!(model[coedge].orientation, orientation.reversed());
}

#[test]
fn remove_coedge_detaches_it_from_wire_and_edge() {
    let (mut model, _) = box_model();
    let wire = first_wire(&model);
    let before = model[wire].coedges.len();
    let coedge = model[wire].coedges[0];
    let edge = model[coedge].edge;
    model.remove_coedge(wire, coedge).unwrap();
    assert_eq!(model[wire].coedges.len(), before - 1);
    assert!(model.get(coedge).is_none());
    assert!(!model.ancestors_of::<Wire>(edge).contains(&wire));
    assert!(validity_ok(&model));
}

fn validity_ok(model: &Model) -> bool {
    kentron::validity::check_model_structure(model).is_valid()
}

#[test]
fn removing_a_compound_child_keeps_the_child() {
    let (mut model, solid) = box_model();
    let c = compound(&mut model, &[solid]);
    model.remove(c, solid).unwrap();
    assert_eq!(model.children(c).count(), 0);
    assert!(model.get(id_of::<kentron::topology::Solid>(solid)).is_some());
    assert!(validity_ok(&model));
}

#[test]
fn shared_child_has_one_use_per_compound() {
    let (mut model, solid) = box_model();
    let a = compound(&mut model, &[solid]);
    let b = compound(&mut model, &[]);
    model.add(b, solid).unwrap();
    let use_a = model[id_of::<kentron::topology::Compound>(a)].children[0];
    let use_b = model[id_of::<kentron::topology::Compound>(b)].children[0];
    assert_ne!(use_a, use_b);
    let parents = |m: &Model| m.parents(solid.node);
    assert!(parents(&model).contains(&a.node) && parents(&model).contains(&b.node));
    model.remove_use(use_a);
    assert!(model.get(use_a).is_none() && model.get(use_b).is_some());
    assert!(!parents(&model).contains(&a.node) && parents(&model).contains(&b.node));
    assert_eq!(model.children(a).count(), 0);
    assert_eq!(model.children(b).count(), 1);
    model.remove_use(use_b);
    assert!(!parents(&model).contains(&b.node));
    assert!(validity_ok(&model));
}

#[test]
fn replace_node_reparents_all_coedges() {
    let (mut model, _) = box_model();
    let wire = first_wire(&model);
    let old = model[model[wire].coedges[0]].edge;
    let new = EdgeId::from_index((old.index() + 1) % model.count(NodeKind::Edge) as u32);
    model.replace_node(old, new).unwrap();
    let coedges = &model.table::<kentron::topology::CoEdge>();
    assert!(coedges.iter().all(|(_, c)| c.edge != old));
    assert!(coedges.iter().filter(|(_, c)| c.edge == new).count() > 0);
}

#[test]
fn box_shell_and_solid_structure() {
    let (model, _) = box_model();
    let shell = first::<Shell>(&model);
    assert_eq!(model[shell].faces.len(), 6);
    let solid = first::<kentron::topology::Solid>(&model);
    assert_eq!(model[solid].shells.len(), 1);
    for (face, f) in model.table::<Face>().iter() {
        assert_eq!(f.wires.len(), 1, "box faces have no inner wires");
        assert!(f.tolerance >= 0.0);
        assert_eq!(model.outer_wire(face), Some(model[f.wires[0]].wire));
    }
    for (_, v) in model.table::<kentron::topology::Vertex>().iter() {
        assert!(v.tolerance > 0.0);
    }
    for (edge, e) in model.table::<Edge>().iter() {
        assert!(e.tolerance > 0.0);
        let (first, last) = model.edge_range(edge);
        assert!(first < last);
    }
}
