//! Helpers shared by the BRepGraph_* ports.

use kentron::topology::{Id, Model, Node, NodeKind, Record, Shape};

use crate::make_box;

/// A 10 x 20 x 30 box (the fixture of most BRepGraph tests).
pub fn box_model() -> (Model, Shape) {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    (model, solid)
}

/// A compound of `children`.
pub fn compound(model: &mut Model, children: &[Shape]) -> Shape {
    let c = model.new_empty(NodeKind::Compound);
    for &child in children {
        model.add(c, child).unwrap();
    }
    c
}

/// The first (non-removed) `T` record.
pub fn first<T: Record>(model: &Model) -> Id<T> {
    model.table::<T>().iter().next().expect("no record of that kind").0
}

/// The typed id of a shape's node; panics on a kind mismatch.
pub fn id_of<T: Node>(shape: Shape) -> Id<T> {
    Id::try_from(shape.node).expect("unexpected node kind")
}
