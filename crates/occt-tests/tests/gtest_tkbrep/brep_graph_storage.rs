//! BRepGraphInc_Test.cxx: the incidence storage underneath BRepGraph.
//!
//! Most of that suite tests OCCT's import (`BRepGraphInc_Populate`), reconstruction to
//! TopoDS, its side layer for internal/external sub-shapes, and the incremental reverse
//! index ("relations") kept in sync by every setter. kentron has no TopoDS, and its reverse
//! queries (`ancestors_of`, `parents`) are part of the model's contract, checked by the
//! other modules. Ported here: entity counts and seam/degenerate structure of primitives,
//! volume/area preservation through `Model::import`, and wire coedge-order normalization.
//! Not ported: `Populate_AppliesNormalChildLocationsToDefinitions` and the stacked variant
//! (a placed edge in a wire; kentron wires hold no placements), and
//! `Relations_RemoveWireFromSharedWire*` (a wire shared by two faces, which kentron's
//! coedges, bound to one face, cannot express).

use std::collections::HashSet;

use kentron::math::{Point, Transform, Vector};
use kentron::primitives::{Cuboid, Cylinder, Sphere};
use kentron::properties;
use kentron::topology::{CoEdge, Edge, Face, Model, NodeKind, Orientation, Shape, VertexId};

use crate::graph_support::{box_model, compound};

#[test]
fn box_coedge_count() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::CoEdge), 24);
}

#[test]
fn cylinder_has_a_seam() {
    let mut model = Model::new();
    Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
    let seams = model.table::<CoEdge>().iter().filter(|&(c, _)| model.is_seam(c)).count();
    assert!(seams >= 2, "a seam is a pair of coedges");
}

#[test]
fn cylinder_seam_edge_lists_its_face_once() {
    let mut model = Model::new();
    Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
    for (edge, _) in model.table::<Edge>().iter() {
        let faces = model.ancestors_of::<Face>(edge);
        let unique: HashSet<_> = faces.iter().collect();
        assert_eq!(unique.len(), faces.len(), "{edge:?} lists a face twice");
    }
}

#[test]
fn sphere_degenerate_edges_have_no_curve() {
    let mut model = Model::new();
    Sphere::new(8.0).build(&mut model).unwrap();
    let degenerate: Vec<_> =
        model.table::<Edge>().iter().filter(|&(e, _)| model.is_degenerate(e)).collect();
    assert!(degenerate.len() >= 2, "a sphere has degenerate edges at its poles");
    for (edge, _) in degenerate {
        assert!(model.edge_curve(edge).is_none(), "{edge:?} is degenerate but has a curve");
    }
}

fn assert_import_preserves(build: impl Fn(&mut Model) -> Shape, volume: bool) {
    let mut source = Model::new();
    let shape = build(&mut source);
    let mut model = Model::new();
    let copy = model.import(&source, shape).unwrap();
    let a = properties::area(&source, shape).mass;
    assert!((properties::area(&model, copy).mass - a).abs() <= 1.0e-7 * a.max(1.0));
    if volume {
        let v = properties::volume(&source, shape).mass;
        assert!((properties::volume(&model, copy).mass - v).abs() <= 1.0e-7 * v.max(1.0));
    }
    for kind in [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
        assert_eq!(model.explore(copy, kind).count(), source.explore(shape, kind).count());
    }
}

#[test]
fn import_preserves_box() {
    assert_import_preserves(|m| Cuboid::new(10.0, 20.0, 30.0).build(m).unwrap(), true);
}

#[test]
fn import_preserves_cylinder() {
    assert_import_preserves(|m| Cylinder::new(5.0, 15.0).build(m).unwrap(), true);
}

#[test]
fn import_preserves_sphere() {
    assert_import_preserves(|m| Sphere::new(8.0).build(m).unwrap(), false);
}

#[test]
fn compound_with_translated_child_keeps_volume() {
    let mut model = Model::new();
    let a = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let moved = kentron::modify::transform(
        &mut model,
        b,
        &Transform::translation(Vector::new(100.0, 0.0, 0.0)),
        Default::default(),
    )
    .unwrap();
    let c = compound(&mut model, &[a, moved]);
    assert!((properties::volume(&model, c).mass - 2000.0).abs() <= 1.0e-7);
    let mut copy = Model::new();
    let c2 = copy.import(&model, c).unwrap();
    assert!((properties::volume(&copy, c2).mass - 2000.0).abs() <= 1.0e-7);
    assert_eq!(copy.explore(c2, NodeKind::Solid).count(), 2);
    assert_eq!(copy.explore(c2, NodeKind::Face).count(), 12);
}

// ---- wire coedge order ----

fn vertex(model: &mut Model, x: f64, tolerance: f64) -> VertexId {
    model.add_vertex(Point::new(x, 0.0, 0.0), tolerance)
}

fn coedge(model: &mut Model, a: VertexId, b: VertexId) -> kentron::topology::CoEdgeId {
    let edge = model.add_edge(a, b, 1.0e-7).unwrap();
    model.add_coedge(edge, Orientation::Forward).unwrap()
}

#[test]
fn coedge_order_normalizes_connected_chain() {
    let mut m = Model::new();
    let (a, b, c) = (vertex(&mut m, 0.0, 1.0e-7), vertex(&mut m, 1.0, 1.0e-7), vertex(&mut m, 2.0, 1.0e-7));
    let bc = coedge(&mut m, b, c);
    let ab = coedge(&mut m, a, b);
    let wire = m.add_wire(&[bc, ab]).unwrap();
    assert!(m.set_coedge_order(wire, &[bc, ab]).unwrap());
    assert_eq!(m[wire].coedges, vec![ab, bc]);
    assert!(!m.set_coedge_order(wire, &[ab, bc]).unwrap());
}

#[test]
fn coedge_order_connects_within_vertex_tolerance() {
    let mut m = Model::new();
    let a = vertex(&mut m, 0.0, 1.0e-7);
    let b = vertex(&mut m, 1.0, 1.0e-2);
    let c = vertex(&mut m, 1.005, 1.0e-2);
    let d = vertex(&mut m, 2.0, 1.0e-7);
    let ab = coedge(&mut m, a, b);
    let cd = coedge(&mut m, c, d);
    let wire = m.add_wire(&[cd, ab]).unwrap();
    m.set_coedge_order(wire, &[cd, ab]).unwrap();
    assert_eq!(m[wire].coedges, vec![ab, cd]);
}

#[test]
fn coedge_order_keeps_disconnected_runs() {
    let mut m = Model::new();
    let v: Vec<VertexId> = [0.0, 1.0, 2.0, 10.0, 11.0].iter().map(|&x| vertex(&mut m, x, 1.0e-7)).collect();
    let bc = coedge(&mut m, v[1], v[2]);
    let xy = coedge(&mut m, v[3], v[4]);
    let ab = coedge(&mut m, v[0], v[1]);
    let wire = m.add_wire(&[bc, xy, ab]).unwrap();
    m.set_coedge_order(wire, &[bc, xy, ab]).unwrap();
    assert_eq!(m[wire].coedges, vec![ab, bc, xy]);
}
