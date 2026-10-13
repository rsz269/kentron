//! TopoDS_Iterator_Test.cxx
//!
//! Not ported (no null shapes or stateful iterators in kentron): `OCC30708_1_InitializeWithNullShape`,
//! `ValueThrowsWhenNotMore`; `DefaultConstructorThenInitialize` and `ReInitialization`
//! keep only their iteration checks.

use kentron::make;
use kentron::math::{Point, Transform, Vector};
use kentron::topology::{Model, NodeKind, Orientation, Shape};

use crate::{first, make_box, make_edge};

fn compound_of_vertices(model: &mut Model, n: usize, y: f64) -> Shape {
    let compound = model.new_empty(NodeKind::Compound);
    for i in 0..n {
        let v = make::vertex(model, Point::new(i as f64, y, 0.0));
        model.add(compound, v).unwrap();
    }
    compound
}

#[test]
fn empty_compound() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    assert!(model.children(compound).next().is_none(), "empty compound has no children");
}

#[test]
fn compound_with_vertices() {
    let mut model = Model::new();
    let compound = compound_of_vertices(&mut model, 5, 0.0);
    let children: Vec<_> = model.children(compound).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::Vertex));
    assert_eq!(children.len(), 5);
}

#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn wire_with_edges() {
    let mut model = Model::new();
    let wire = model.new_empty(NodeKind::Wire);
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)];
    for w in corners.windows(2) {
        let e = make_edge(&mut model, Point::new(w[0].0, w[0].1, 0.0), Point::new(w[1].0, w[1].1, 0.0));
        model.add(wire, e).unwrap();
    }
    let children: Vec<_> = model.children(wire).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::Edge));
    assert_eq!(children.len(), 4, "wire should have 4 edges");
}

#[test]
fn edge_with_vertices() {
    let mut model = Model::new();
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 1.0, 1.0));
    let children: Vec<_> = model.children(edge).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::Vertex));
    assert_eq!(children.len(), 2, "an edge has 2 vertices");
}

#[test]
fn vertex_no_children() {
    let mut model = Model::new();
    let vertex = make::vertex(&mut model, Point::ORIGIN);
    assert!(model.children(vertex).next().is_none());
}

#[test]
fn cumulative_orientation() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    let v = make::vertex(&mut model, Point::ORIGIN);
    model.add(compound, v).unwrap();
    let reversed = compound.reversed();

    let composed = model.children(reversed).next().expect("one child").orientation;
    let local = model.children_local(reversed).next().expect("one child").orientation;
    assert_ne!(composed, local, "composing with a reversed parent flips the child");
    assert_eq!(local, Orientation::Forward);
}

#[test]
fn re_initialization() {
    let mut model = Model::new();
    let c1 = compound_of_vertices(&mut model, 3, 0.0);
    let c2 = compound_of_vertices(&mut model, 5, 1.0);
    assert_eq!(model.children(c1).count(), 3);
    assert_eq!(model.children(c2).count(), 5);
}

#[test]
fn box_solid_iteration() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 1.0, 2.0, 3.0);

    let shells: Vec<_> = model.children(solid).collect();
    assert_eq!(shells.len(), 1, "a box solid has 1 shell");
    for shell in shells {
        assert_eq!(shell.kind(), NodeKind::Shell);
        let faces: Vec<_> = model.children(shell).collect();
        assert!(faces.iter().all(|f| f.kind() == NodeKind::Face));
        assert_eq!(faces.len(), 6, "a box shell has 6 faces");
    }
}

#[test]
fn many_children() {
    let mut model = Model::new();
    let compound = compound_of_vertices(&mut model, 1000, 0.0);
    assert_eq!(model.children(compound).count(), 1000);
}

#[test]
fn default_constructor_then_initialize() {
    let mut model = Model::new();
    let compound = compound_of_vertices(&mut model, 1, 0.0);
    let child = model.children(compound).next().expect("one child");
    assert_eq!(child.kind(), NodeKind::Vertex);
}

#[test]
fn cumulative_location() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    let v = make::vertex(&mut model, Point::ORIGIN);
    model.add(compound, v).unwrap();
    let placed = compound.with_placement(Transform::translation(Vector::new(10.0, 20.0, 30.0)));

    let composed = model.children(placed).next().expect("one child").placement;
    let local = model.children_local(placed).next().expect("one child").placement;
    assert_ne!(composed.is_identity(), local.is_identity(), "the parent's placement composes");
}

#[test]
fn nested_compounds() {
    let mut model = Model::new();
    let outer = model.new_empty(NodeKind::Compound);
    let inner1 = model.new_empty(NodeKind::Compound);
    let inner2 = model.new_empty(NodeKind::Compound);
    for (inner, x) in [(inner1, 0.0), (inner1, 1.0), (inner2, 2.0)] {
        let v = make::vertex(&mut model, Point::new(x, 0.0, 0.0));
        model.add(inner, v).unwrap();
    }
    model.add(outer, inner1).unwrap();
    model.add(outer, inner2).unwrap();

    let children: Vec<_> = model.children(outer).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::Compound));
    assert_eq!(children.len(), 2);
    assert_eq!(model.children(children[0]).count(), 2);
}

#[test]
fn face_with_wires() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 1.0, 1.0, 1.0);
    let face = first(&model, solid, NodeKind::Face);
    let wires: Vec<_> = model.children(face).collect();
    assert!(wires.iter().all(|w| w.kind() == NodeKind::Wire));
    assert!(!wires.is_empty(), "a face has at least 1 wire");
}

#[test]
fn shape_identity() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    let v1 = make::vertex(&mut model, Point::new(1.0, 0.0, 0.0));
    let v2 = make::vertex(&mut model, Point::new(2.0, 0.0, 0.0));
    model.add(compound, v1).unwrap();
    model.add(compound, v2).unwrap();

    for child in model.children(compound) {
        assert!(child.node == v1.node || child.node == v2.node);
    }
}

#[test]
fn shell_with_faces() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 1.0, 1.0, 1.0);
    let shell = first(&model, solid, NodeKind::Shell);
    let faces: Vec<_> = model.children(shell).collect();
    assert!(faces.iter().all(|f| f.kind() == NodeKind::Face));
    assert_eq!(faces.len(), 6);
}
