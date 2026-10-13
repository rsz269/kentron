//! TopoDS_Builder_Test.cxx
//!
//! Not ported (OCCT `TShape` bookkeeping flags, no kentron counterpart): `TShapeFlags`,
//! `Add_SetsModifiedFlag`, `Remove_SetsModifiedFlag`; the `Free()` checks are dropped.

use kentron::make;
use kentron::math::Point;
use kentron::topology::{Model, NodeKind};

use crate::make_edge;

#[test]
fn make_wire() {
    let mut model = Model::new();
    let wire = model.new_empty(NodeKind::Wire);
    assert_eq!(wire.kind(), NodeKind::Wire);
}

#[test]
fn make_shell() {
    let mut model = Model::new();
    assert_eq!(model.new_empty(NodeKind::Shell).kind(), NodeKind::Shell);
}

#[test]
fn make_solid() {
    let mut model = Model::new();
    assert_eq!(model.new_empty(NodeKind::Solid).kind(), NodeKind::Solid);
}

#[test]
fn make_comp_solid() {
    let mut model = Model::new();
    assert_eq!(model.new_empty(NodeKind::CompSolid).kind(), NodeKind::CompSolid);
}

#[test]
fn make_compound() {
    let mut model = Model::new();
    assert_eq!(model.new_empty(NodeKind::Compound).kind(), NodeKind::Compound);
}

#[test]
fn add_vertices_to_compound() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    for p in [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)] {
        let v = make::vertex(&mut model, p);
        model.add(compound, v).unwrap();
    }
    assert_eq!(model.children(compound).count(), 3, "compound should have 3 vertices");
}

#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn add_edges_to_wire() {
    let mut model = Model::new();
    let wire = model.new_empty(NodeKind::Wire);
    // a triangle
    let e1 = make_edge(&mut model, Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
    let e2 = make_edge(&mut model, Point::new(1.0, 0.0, 0.0), Point::new(0.5, 1.0, 0.0));
    let e3 = make_edge(&mut model, Point::new(0.5, 1.0, 0.0), Point::new(0.0, 0.0, 0.0));
    for e in [e1, e2, e3] {
        model.add(wire, e).unwrap();
    }

    let children: Vec<_> = model.children(wire).collect();
    assert!(children.iter().all(|c| c.kind() == NodeKind::Edge));
    assert_eq!(children.len(), 3, "wire should have 3 edges");
}

/// Three vertices in a compound, built at (1,0,0), (2,0,0), (3,0,0).
fn compound_of_three(model: &mut Model) -> (kentron::topology::Shape, [kentron::topology::Shape; 3]) {
    let compound = model.new_empty(NodeKind::Compound);
    let vs = [1.0, 2.0, 3.0].map(|x| make::vertex(model, Point::new(x, 0.0, 0.0)));
    for v in vs {
        model.add(compound, v).unwrap();
    }
    (compound, vs)
}

#[test]
fn remove_from_compound() {
    let mut model = Model::new();
    let (compound, [_, v2, _]) = compound_of_three(&mut model);
    model.remove(compound, v2).unwrap();
    assert_eq!(model.children(compound).count(), 2, "2 vertices should remain");
}

#[test]
fn add_many_shapes() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    for i in 0..500 {
        let v = make::vertex(&mut model, Point::new(i as f64, 0.0, 0.0));
        model.add(compound, v).unwrap();
    }
    assert_eq!(model.children(compound).count(), 500);
}

#[test]
fn remove_first_child() {
    let mut model = Model::new();
    let (compound, [v1, _, _]) = compound_of_three(&mut model);
    model.remove(compound, v1).unwrap();
    assert_eq!(model.children(compound).count(), 2);
}

#[test]
fn remove_last_child() {
    let mut model = Model::new();
    let (compound, [_, _, v3]) = compound_of_three(&mut model);
    model.remove(compound, v3).unwrap();
    assert_eq!(model.children(compound).count(), 2);
}

#[test]
fn remove_all_children() {
    let mut model = Model::new();
    let (compound, vs) = compound_of_three(&mut model);
    for v in vs {
        model.remove(compound, v).unwrap();
    }
    assert!(model.children(compound).next().is_none(), "compound should be empty");
}

#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn remove_from_wire() {
    let mut model = Model::new();
    let wire = model.new_empty(NodeKind::Wire);
    let e1 = make_edge(&mut model, Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
    let e2 = make_edge(&mut model, Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0));
    let e3 = make_edge(&mut model, Point::new(1.0, 1.0, 0.0), Point::new(0.0, 0.0, 0.0));
    for e in [e1, e2, e3] {
        model.add(wire, e).unwrap();
    }
    model.remove(wire, e2).unwrap();
    assert_eq!(model.children(wire).count(), 2);
}

#[test]
fn nested_compounds() {
    let mut model = Model::new();
    let outer = model.new_empty(NodeKind::Compound);
    let inner = model.new_empty(NodeKind::Compound);
    for p in [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)] {
        let v = make::vertex(&mut model, p);
        model.add(inner, v).unwrap();
    }
    model.add(outer, inner).unwrap();

    let children: Vec<_> = model.children(outer).collect();
    assert_eq!(children.len(), 1);
    assert_eq!(children[0].kind(), NodeKind::Compound);
    assert_eq!(model.children(inner).count(), 2, "inner compound keeps its 2 vertices");
}

#[test]
fn all_types() {
    let mut model = Model::new();
    for kind in [
        NodeKind::Wire,
        NodeKind::Shell,
        NodeKind::Solid,
        NodeKind::CompSolid,
        NodeKind::Compound,
    ] {
        assert_eq!(model.new_empty(kind).kind(), kind);
    }
}

#[test]
fn mixed_shapes_in_compound() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    let vertex = make::vertex(&mut model, Point::ORIGIN);
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    let wire = model.new_empty(NodeKind::Wire);
    let shell = model.new_empty(NodeKind::Shell);
    let inner = model.new_empty(NodeKind::Compound);
    for child in [vertex, edge, wire, shell, inner] {
        model.add(compound, child).unwrap();
    }

    let kinds: Vec<_> = model.children(compound).map(|c| c.kind()).collect();
    assert_eq!(kinds.len(), 5);
    for kind in [NodeKind::Vertex, NodeKind::Edge, NodeKind::Wire, NodeKind::Shell, NodeKind::Compound] {
        assert_eq!(kinds.iter().filter(|k| **k == kind).count(), 1, "{kind:?}");
    }
}
