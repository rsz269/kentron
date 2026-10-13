//! TopoDS_TShape_Test.cxx and TopoDS_Edge_Test.cxx
//!
//! Not ported (OCCT `TShape` flags, no kentron counterpart): `FlagSettersGetters`,
//! `Orientable_DifferentTypes`, `FlagsIndependent`, and TopoDS_Edge_Test's
//! `BUC60828_InfiniteFlag`. `NbChildren` is the number of [`Model::children`].
//!
//! [`Model::children`]: kentron::topology::Model::children

use kentron::make;
use kentron::math::Point;
use kentron::topology::{Model, NodeKind};

use crate::make_edge;

#[test]
fn shape_type_all_types() {
    let mut model = Model::new();
    assert_eq!(make::vertex(&mut model, Point::ORIGIN).kind(), NodeKind::Vertex);
    assert_eq!(
        make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).kind(),
        NodeKind::Edge
    );
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
fn nb_children() {
    let mut model = Model::new();

    let vertex = make::vertex(&mut model, Point::ORIGIN);
    assert_eq!(model.children(vertex).count(), 0, "a vertex has no children");

    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    assert_eq!(model.children(edge).count(), 2, "an edge has its 2 vertices");

    let wire = model.new_empty(NodeKind::Wire);
    assert_eq!(model.children(wire).count(), 0);
    model.add(wire, edge).unwrap();
    assert_eq!(model.children(wire).count(), 1);

    let compound = model.new_empty(NodeKind::Compound);
    for i in 0..10 {
        let v = make::vertex(&mut model, Point::new(i as f64, 0.0, 0.0));
        model.add(compound, v).unwrap();
    }
    assert_eq!(model.children(compound).count(), 10);
}

#[test]
fn empty_copy() {
    let mut model = Model::new();
    let compound = model.new_empty(NodeKind::Compound);
    for p in [Point::ORIGIN, Point::new(1.0, 0.0, 0.0)] {
        let v = make::vertex(&mut model, p);
        model.add(compound, v).unwrap();
    }
    assert_eq!(model.children(compound).count(), 2);

    let copy = model.empty_copy(compound);
    assert_eq!(copy.kind(), NodeKind::Compound);
    assert_eq!(model.children(copy).count(), 0, "an empty copy has no children");
}

#[test]
fn empty_copy_all_types() {
    let mut model = Model::new();

    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    let copy = model.empty_copy(edge);
    assert_eq!(copy.kind(), NodeKind::Edge);
    assert_eq!(model.children(copy).count(), 0);

    for kind in [NodeKind::Wire, NodeKind::Shell, NodeKind::Solid, NodeKind::CompSolid] {
        let shape = model.new_empty(kind);
        assert_eq!(model.empty_copy(shape).kind(), kind);
    }
}

#[test]
fn nb_children_consistency_with_iterator() {
    let mut model = Model::new();
    let wire = model.new_empty(NodeKind::Wire);
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    model.add(wire, edge).unwrap();
    assert_eq!(model.children(wire).count(), 1);
}
