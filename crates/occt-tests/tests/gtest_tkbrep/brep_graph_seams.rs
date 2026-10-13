//! BRepGraph_SeamRedesign_Test.cxx and BRepGraph_ReplaceVertex_Test.cxx
//!
//! Not ported: `CylinderLateralWire_NbCoEdges_MatchesTopoDSIterator` and
//! `Reconstruct_CylinderWire_TopoDSIteratorOrder` (TopoDS round trips); the validator's
//! issue wording in `Validate_DetectsAsymmetricSeamPair`.

use kentron::make;
use kentron::math::Point;
use kentron::primitives;
use kentron::topology::{CoEdge, CoEdgeId, Edge, EdgeId, Model, Wire};
use kentron::validity;

use crate::graph_support::{box_model, first};

fn build(f: impl FnOnce(&mut Model)) -> Model {
    let mut model = Model::new();
    f(&mut model);
    model
}

fn cylinder() -> Model {
    build(|m| {
        primitives::Cylinder::new(5.0, 10.0).build(m).unwrap();
    })
}

fn find_seam(model: &Model) -> CoEdgeId {
    model.table::<CoEdge>().iter().map(|(id, _)| id).find(|&c| model.is_seam(c)).expect("a seam coedge")
}

/// Both halves of the cylinder's seam are in its lateral wire.
#[test]
fn cylinder_lateral_wire_both_seam_halves_in_co_edge_ids() {
    let model = cylinder();
    let seam = find_seam(&model);
    let pair = model.seam_pair(seam).expect("seam pair");
    let wire = &model[model[seam].wire];
    assert!(wire.coedges.contains(&seam));
    assert!(wire.coedges.contains(&pair));
}

#[test]
fn seam_pair_derived_query_is_symmetric() {
    let model = build(|m| {
        primitives::Sphere::new(10.0).build(m).unwrap();
    });
    let mut seams = 0;
    for (coedge, _) in model.table::<CoEdge>().iter() {
        if let Some(pair) = model.seam_pair(coedge) {
            seams += 1;
            assert_eq!(model.seam_pair(pair), Some(coedge));
        }
    }
    assert!(seams > 0, "a sphere has a seam");
}

#[test]
fn seam_pair_box_has_no_seams() {
    let (model, _) = box_model();
    for (coedge, _) in model.table::<CoEdge>().iter() {
        assert!(!model.is_seam(coedge), "{coedge:?}");
    }
}

#[test]
fn seam_pair_free_wire_has_no_seams() {
    let model = build(|m| {
        let edge = make::edge_between(m, Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).unwrap();
        make::wire(m, &[edge]).unwrap();
    });
    for (coedge, _) in model.table::<CoEdge>().iter() {
        assert!(!model.is_seam(coedge));
    }
}

#[test]
fn edge_ops_is_seam_on_face_derived_from_connectivity() {
    let model = cylinder();
    let seam = &model[find_seam(&model)];
    assert!(model.is_seam_on_face(seam.edge, seam.face));

    let (box_model, _) = box_model();
    for (_, coedge) in box_model.table::<CoEdge>().iter() {
        assert!(!box_model.is_seam_on_face(coedge.edge, coedge.face));
    }
}

/// A seam edge appears twice in its wire.
#[test]
fn nb_distinct_edges_accounts_for_seam_halves() {
    let model = cylinder();
    let mut wires_with_seams = 0;
    for (_, wire) in model.table::<Wire>().iter() {
        let mut edges: Vec<EdgeId> = wire.coedges.iter().map(|&c| model[c].edge).collect();
        let raw = edges.len();
        edges.sort();
        edges.dedup();
        assert!(edges.len() <= raw);
        if raw > edges.len() {
            wires_with_seams += 1;
        }
    }
    assert!(wires_with_seams >= 1);
}

/// Giving both seam halves the same orientation breaks the seam and the model check.
#[test]
fn validate_detects_asymmetric_seam_pair() {
    let mut model = cylinder();
    let seam = find_seam(&model);
    let mate = model.seam_pair(seam).expect("pair");
    let orientation = model[seam].orientation;
    model.set_coedge_orientation(mate, orientation).unwrap();

    let (edge, face) = (model[seam].edge, model[seam].face);
    assert!(!model.is_seam_on_face(edge, face));
    assert!(!validity::check_model(&model).is_valid());
}

// ReplaceVertex

#[test]
fn start_vertex_swapped_to_fresh_vertex_audit_clean() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let old = model[edge].start;
    let vertex = model.add_vertex(Point::new(99.0, 99.0, 99.0), 1.0e-7);

    let new = model.replace_edge_vertex(edge, old, vertex).unwrap();
    assert_ne!(new, old);
    assert_eq!(model[edge].start, new);
    assert!(model.get(old).is_none(), "the old use is retired");
    assert_eq!(model[new].vertex, vertex);
    // only structural checks: replacing one edge's end breaks wire connectivity
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn end_vertex_swapped_to_fresh_vertex_preserves_orientation() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let old = model[edge].end;
    let orientation = model[old].orientation;
    let vertex = model.add_vertex(Point::new(1.0, 2.0, 3.0), 1.0e-7);

    let new = model.replace_edge_vertex(edge, old, vertex).unwrap();
    assert_eq!(model[new].orientation, orientation);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn same_vertex_idempotent() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let old = model[edge].start;
    let same = model[old].vertex;
    assert_eq!(model.replace_edge_vertex(edge, old, same).unwrap(), old);
    assert!(model.get(old).is_some());
}

#[test]
fn inactive_edge_rejected() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let old = model[edge].start;
    let vertex = model.add_vertex(Point::new(5.0, 5.0, 5.0), 1.0e-7);
    model.remove_node(edge);
    assert!(model.replace_edge_vertex(edge, old, vertex).is_err());
}

#[test]
fn wrong_parent_rejected() {
    let (mut model, _) = box_model();
    let edges: Vec<EdgeId> = model.table::<Edge>().iter().map(|(id, _)| id).take(2).collect();
    let other_start = model[edges[1]].start;
    let vertex = model.add_vertex(Point::new(0.0, 0.0, 100.0), 1.0e-7);
    assert!(model.replace_edge_vertex(edges[0], other_start, vertex).is_err());
}

#[test]
fn relations_rebuilt_vertex_to_edge() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let old = model[edge].start;
    let old_vertex = model[old].vertex;
    let vertex = model.add_vertex(Point::new(1000.0, 0.0, 0.0), 1.0e-7);
    model.replace_edge_vertex(edge, old, vertex).unwrap();

    assert_eq!(model.ancestors_of::<Edge>(vertex), [edge]);
    assert!(!model.ancestors_of::<Edge>(old_vertex).contains(&edge));
    assert!(validity::check_model_structure(&model).is_valid());
}
