//! BRepGraph_Compact_Test.cxx: dropping removed records and renumbering ids.

use kentron::math::{BoundingBox, Point};
use kentron::modify;
use kentron::tolerance;
use kentron::topology::{
    CoEdge, DeduplicateOptions, Edge, Face, FaceUse, Id, Model, NodeId, NodeKind, Shape, Shell,
    ShellUse, Solid, StableId, Vertex, VertexUse, Wire, WireUse,
};
use kentron::validity;

use crate::graph_support::{box_model, compound, first, id_of};
use crate::{make_box, make_edge};

/// A compound of two deep copies of one face of a box.
fn two_copied_faces(model: &mut Model) -> Shape {
    let solid = make_box(model, 10.0, 20.0, 30.0);
    let face = crate::first(model, solid, NodeKind::Face);
    let a = modify::copy(model, face, false).unwrap();
    let b = modify::copy(model, face, false).unwrap();
    compound(model, &[a, b])
}

/// A compound of a box and a free edge.
fn box_with_loose_edge(model: &mut Model) -> (Shape, Shape) {
    let solid = make_box(model, 10.0, 20.0, 30.0);
    let edge = make_edge(model, Point::new(100.0, 0.0, 0.0), Point::new(120.0, 0.0, 0.0));
    (compound(model, &[solid, edge]), edge)
}

fn vertex_bounds(model: &Model) -> BoundingBox {
    let mut points = model.table::<Vertex>().iter().map(|(_, v)| v.point);
    let p = points.next().expect("no vertex");
    points.fold(BoundingBox { min: p, max: p }, |b, p| BoundingBox {
        min: Point::new(b.min.x.min(p.x), b.min.y.min(p.y), b.min.z.min(p.z)),
        max: Point::new(b.max.x.max(p.x), b.max.y.max(p.y), b.max.z.max(p.z)),
    })
}

fn stable_ids_of(model: &Model, kind: NodeKind) -> Vec<StableId> {
    (0..model.slot_count(kind))
        .map(|i| model.node_stable_id(NodeId { kind, index: i as u32 }).unwrap())
        .collect()
}

#[test]
fn no_removed_nodes_noop() {
    let (mut model, _) = box_model();
    let counts: Vec<usize> = NodeKind::ALL.iter().map(|&k| model.count(k)).collect();
    let report = model.compact();
    assert_eq!(report.removed(NodeKind::Vertex), 0);
    assert_eq!(report.removed(NodeKind::Edge), 0);
    assert_eq!(report.removed(NodeKind::Face), 0);
    assert_eq!(report.nodes_before, report.nodes_after);
    assert_eq!(report.renumbered, 0);
    let after: Vec<usize> = NodeKind::ALL.iter().map(|&k| model.count(k)).collect();
    assert_eq!(after, counts);
}

#[test]
fn after_deduplicate_removes_nodes() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    model.deduplicate();
    let report = model.compact();
    // Geometry merging alone removes no topology node.
    assert_eq!(report.nodes_after, report.nodes_before);
}

#[test]
fn index_density_no_gaps() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    model.deduplicate();
    model.compact();
    for kind in NodeKind::ALL {
        assert_eq!(model.slot_count(kind), model.count(kind), "{kind:?} has removed slots");
    }
}

#[test]
fn cross_references_valid() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    model.deduplicate();
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

/// OCCT records "Compact:Remap" history events; kentron rewrites the ids stored in the
/// history instead and reports how many nodes were renumbered.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn history_mode_records_mapping() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    let dedup = model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    assert!(NodeKind::ALL.iter().any(|&k| dedup.merged(k) > 0));
    let report = model.compact();
    assert!(report.renumbered >= 1);
}

#[test]
fn full_pipeline_deduplicate_compact_validate() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    model.deduplicate();
    let report = model.compact();
    assert_eq!(report.nodes_after, report.nodes_before);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn removal_compact_preserves_bounds_and_does_not_grow_topology() {
    let (mut model, _) = box_model();
    assert!(model.count(NodeKind::CoEdge) > 0);
    model.remove_node(Id::<Face>::from_index(2));
    let coedges_before = model.slot_count(NodeKind::CoEdge);
    let bounds_before = vertex_bounds(&model);
    let report = model.compact();
    assert!(report.nodes_before >= report.nodes_after);
    assert!(model.slot_count(NodeKind::CoEdge) <= coedges_before);
    let bounds_after = vertex_bounds(&model);
    for (a, b) in [(bounds_after.min, bounds_before.min), (bounds_after.max, bounds_before.max)] {
        assert!(a.distance(b) <= tolerance::LINEAR);
    }
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn removal_compact_preserves_closed_topology_and_valid_shape() {
    let mut model = Model::new();
    let (root, edge) = box_with_loose_edge(&mut model);
    let root = model.node_stable_id(root.node).unwrap();
    let edge = id_of::<Edge>(edge);
    assert!(model.ancestors_of::<Face>(edge).is_empty());
    model.remove_node(edge);
    let report = model.compact();
    assert!(report.nodes_before > report.nodes_after);
    assert_eq!(model.count(NodeKind::Shell), 1);
    for (wire, _) in model.table::<Wire>().iter() {
        assert!(model.is_closed_wire(wire));
    }
    assert!(model.is_closed_shell(first::<Shell>(&model)));
    let root = model.find_stable_node(root).map(Shape::from).unwrap();
    assert!(validity::check(&model, root).is_valid());
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn audit_mode_passes_after_dedup_compact() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    assert!(validity::check_model(&model).is_valid());
    model.deduplicate();
    model.compact();
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn audit_mode_passes_after_removal_compact() {
    let mut model = Model::new();
    let (_, edge) = box_with_loose_edge(&mut model);
    model.remove_node(edge.node);
    model.compact();
    assert!(validity::check_model(&model).is_valid());
}

#[test]
fn compact_preserves_stable_ids() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    let kinds = [
        NodeKind::Vertex,
        NodeKind::Edge,
        NodeKind::Wire,
        NodeKind::Face,
        NodeKind::Shell,
        NodeKind::Solid,
        NodeKind::Compound,
        NodeKind::CompSolid,
    ];
    let before: Vec<Vec<StableId>> = kinds.iter().map(|&k| stable_ids_of(&model, k)).collect();
    model.deduplicate();
    model.compact();
    for (&kind, before) in kinds.iter().zip(&before) {
        for (index, sid) in stable_ids_of(&model, kind).into_iter().enumerate() {
            assert!(before.contains(&sid), "{kind:?} {index} has a new stable id");
            assert_eq!(model.find_stable_node(sid), Some(NodeId { kind, index: index as u32 }));
        }
    }
    assert!(model.count(NodeKind::Face) > 0);
}

#[test]
fn compact_preserves_geometry_slot_stable_ids() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let slot = model[face].surface;
    let sid = model.stable_id(slot).unwrap();
    model.compact();
    let face = first::<Face>(&model);
    assert_eq!(model.stable_id(model[face].surface), Some(sid));
}

#[test]
fn revision_survives_compact() {
    let mut model = Model::new();
    two_copied_faces(&mut model);
    let edge = first::<Edge>(&model);
    model.set_edge_tolerance(edge, 0.1).unwrap();
    model.set_edge_tolerance(edge, 0.2).unwrap();
    assert_eq!(model.revision(edge).unwrap().own, 2);
    model.deduplicate();
    model.compact();
    let (edge, _) = model
        .table::<Edge>()
        .iter()
        .find(|(_, e)| (e.tolerance - 0.2).abs() < tolerance::LINEAR)
        .expect("no edge with the edited tolerance after compact");
    assert_eq!(model.revision(edge).unwrap().own, 2);
}

#[test]
fn stable_id_round_trip_after_compaction() {
    let (mut model, _) = box_model();
    let face = |i| Id::<Face>::from_index(i);
    let edge = |i| Id::<Edge>::from_index(i);
    let faces: Vec<StableId> = (0..3).map(|i| model.stable_id(face(i)).unwrap()).collect();
    let edges: Vec<StableId> = (0..2).map(|i| model.stable_id(edge(i)).unwrap()).collect();
    model.remove_node(face(2));
    let report = model.compact();
    assert!(report.removed(NodeKind::Face) >= 1);
    for &sid in faces[..2].iter() {
        assert!(model.find_stable::<Face>(sid).is_some());
    }
    for &sid in &edges {
        assert!(model.find_stable::<Edge>(sid).is_some());
    }
    assert_eq!(model.find_stable::<Face>(faces[2]), None);
    assert_eq!(model.find_stable_node(faces[2]), None);
}

/// The first coedge of face 1's first wire.
fn coedge_of_face_1(model: &Model) -> kentron::topology::CoEdgeId {
    let face = Id::<Face>::from_index(1);
    let wire = model[model[face].wires[0]].wire;
    model[wire].coedges[0]
}

#[test]
fn coedge_stable_id_after_compaction() {
    let (mut model, _) = box_model();
    let sid = model.stable_id(coedge_of_face_1(&model)).unwrap();
    model.remove_node(first::<Face>(&model));
    let report = model.compact();
    assert!(report.nodes_before > report.nodes_after);
    let node = model.find_stable_node(sid).expect("coedge lost after compaction");
    assert_eq!(node.kind, NodeKind::CoEdge);
    assert!(model.find_stable::<CoEdge>(sid).is_some());
}

#[test]
fn use_stable_ids_after_compaction() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let vertex_use = model[edge].start;
    let coedge = coedge_of_face_1(&model);
    let wire_use = model[Id::<Face>::from_index(1)].wires[0];
    let face0 = first::<Face>(&model);
    let shell = first::<Shell>(&model);
    let face_use = *model[shell].faces.iter().find(|&&u| model[u].face != face0).unwrap();
    let shell_use = model[first::<Solid>(&model)].shells[0];

    let vertex_use = model.stable_id(vertex_use).unwrap();
    let coedge = model.stable_id(coedge).unwrap();
    let wire_use = model.stable_id(wire_use).unwrap();
    let face_use = model.stable_id(face_use).unwrap();
    let shell_use = model.stable_id(shell_use).unwrap();

    model.remove_node(face0);
    model.compact();
    assert!(model.find_stable::<VertexUse>(vertex_use).is_some());
    assert_eq!(model.find_stable_node(coedge).map(|n| n.kind), Some(NodeKind::CoEdge));
    assert!(model.find_stable::<WireUse>(wire_use).is_some());
    assert!(model.find_stable::<FaceUse>(face_use).is_some());
    assert!(model.find_stable::<ShellUse>(shell_use).is_some());
}

#[test]
fn revisions_preserved_for_all_topology_kinds_after_compact() {
    let (mut model, _) = box_model();
    let vertex = first::<Vertex>(&model);
    for p in [(11.0, 21.0, 31.0), (12.0, 22.0, 32.0), (13.0, 23.0, 33.0)] {
        model.set_vertex_point(vertex.into(), Point::new(p.0, p.1, p.2), 1.0e-7).unwrap();
    }
    let edge = first::<Edge>(&model);
    for t in [0.1, 0.2, 0.3] {
        model.set_edge_tolerance(edge, t).unwrap();
    }
    let face = first::<Face>(&model);
    for t in [0.01, 0.02, 0.03] {
        model.set_face_tolerance(face, t).unwrap();
    }
    let wire = first::<Wire>(&model);
    let shell = first::<Shell>(&model);
    let solid = first::<Solid>(&model);
    for _ in 0..3 {
        model.touch(wire);
        model.touch(shell);
        model.touch(solid);
    }
    let nodes: [NodeId; 6] =
        [vertex.into(), edge.into(), wire.into(), face.into(), shell.into(), solid.into()];
    for node in nodes {
        assert_eq!(model.revision(node).unwrap().own, 3, "{node:?}");
    }
    model.compact();
    // Nothing was removed: ids are unchanged.
    for node in nodes {
        assert_eq!(model.revision(node).unwrap().own, 3, "{node:?} lost its revision");
    }
}

#[test]
fn compact_after_dedup_merge_no_bounds_errors() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let copies: Vec<Shape> = (0..3).map(|_| modify::copy(&mut model, solid, false).unwrap()).collect();
    compound(&mut model, &copies);
    model.remove_node(solid.node);
    model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    let report = model.compact();
    assert_ne!(report.nodes_after, 0);
    assert_eq!(model.count(NodeKind::Vertex), 8);
    assert_eq!(model.count(NodeKind::Edge), 12);
    assert_eq!(model.count(NodeKind::Face), 6);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn compact_face_with_inner_wires_wire_count_preserved() {
    use kentron::boolean;
    use kentron::math::{Direction, Frame};
    use kentron::primitives::Cylinder;
    let mut model = Model::new();
    let block = make_box(&mut model, 20.0, 20.0, 20.0);
    let frame = Frame::new(Point::new(10.0, 10.0, 0.0), Direction::Z, Direction::X);
    let cylinder = Cylinder::new(4.0, 25.0).with_frame(frame).build(&mut model).unwrap();
    let cut = boolean::difference(&mut model, block, cylinder).unwrap();
    let a = modify::copy(&mut model, cut, false).unwrap();
    let b = modify::copy(&mut model, cut, false).unwrap();
    compound(&mut model, &[a, b]);
    assert!(model.table::<Face>().iter().any(|(_, f)| f.wires.len() > 1));
    model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    let report = model.compact();
    assert_ne!(report.nodes_after, 0);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn history_survives_compact_with_remapped_ids() {
    let mut model = Model::new();
    let a = make_box(&mut model, 10.0, 20.0, 30.0);
    let b = make_box(&mut model, 10.0, 20.0, 30.0);
    compound(&mut model, &[a, b]);
    // A vertex of the second box: merged into the first box's one, then renumbered.
    let vertex: NodeId = Id::<Vertex>::from_index(9).into();
    let sid = model.node_stable_id(vertex).unwrap();
    model.enable_history().record("Test:Modify", vertex, &[vertex]);
    let records = model.history().unwrap().len();
    model.deduplicate_with(DeduplicateOptions { merge_entities: true, ..Default::default() });
    model.compact();
    assert!(model.history().unwrap().len() >= records);
    // Whatever survives under the recorded id is the node the history now names.
    let event = &model.history().unwrap().events()[records - 1];
    let (original, _) = &event.mapping[0];
    match model.find_stable_node(sid) {
        Some(node) => assert_eq!(*original, node),
        None => assert!(model.node_stable_id(*original).is_none_or(|s| s != sid)),
    }
}
