//! BVH_ShapeSelect_Test.cxx (selecting sub-shapes by box with a tree over their boxes)
//!
//! The boxes of a box's vertices, edges and faces (`properties::bounding_box`, OCCT's
//! `BRepBndLib::Add`) in a `math::Bvh`; OCCT's typed and untyped selectors are one
//! `Selector` here, checked against `Bvh::overlapping`.

use kentron::math::{BoundingBox, Bvh, PairSelector, Point, Selector};
use kentron::primitives::Cuboid;
use kentron::properties;
use kentron::topology::{Model, NodeKind, Shape};

/// A 10 × 10 × 10 box, its sub-shapes (vertices, edges, faces) and their boxes.
fn box_set() -> (Model, Vec<(NodeKind, Shape)>, Vec<BoundingBox>) {
    let mut model = Model::new();
    let b = Cuboid::from_corners(Point::ORIGIN, Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
    let mut shapes = Vec::new();
    for kind in [NodeKind::Vertex, NodeKind::Edge, NodeKind::Face] {
        shapes.extend(model.unique_subshapes(b, kind).into_iter().map(|s| (kind, s)));
    }
    let boxes = shapes.iter().map(|&(_, s)| properties::bounding_box(&model, s)).collect();
    (model, shapes, boxes)
}

/// Counts the boxes overlapping a query box.
struct ShapeSelector {
    query: BoundingBox,
    accepted: usize,
}

impl Selector for ShapeSelector {
    fn enter(&mut self, b: &BoundingBox) -> Option<f64> {
        self.query.overlaps(b).then_some(0.0)
    }
    fn accept(&mut self, _: usize, _: f64) -> bool {
        self.accepted += 1;
        true
    }
}

/// Counts the pairs of overlapping boxes.
#[derive(Default)]
struct PairCounter {
    accepted: usize,
}

impl PairSelector for PairCounter {
    fn enter(&mut self, a: &BoundingBox, b: &BoundingBox) -> bool {
        a.overlaps(b)
    }
    fn accept(&mut self, _: usize, _: usize) -> bool {
        self.accepted += 1;
        true
    }
}

/// A box has 8 vertices, 12 edges and 6 faces: selecting by a vertex's, an edge's or a face's
/// box finds 7, 11 or 17 sub-shapes' boxes (itself and its neighbours, touching included).
#[test]
fn bug30655_shape_selection_typed_and_void() {
    let (_, shapes, boxes) = box_set();
    let tree = Bvh::new(&boxes);
    for (kind, expected, count) in [(NodeKind::Vertex, 7, 8), (NodeKind::Edge, 11, 12), (NodeKind::Face, 17, 6)] {
        let of_kind: Vec<usize> = (0..shapes.len()).filter(|&i| shapes[i].0 == kind).collect();
        assert_eq!(of_kind.len(), count);
        for i in of_kind {
            let mut selector = ShapeSelector { query: boxes[i], accepted: 0 };
            assert_eq!(tree.select(&mut selector), expected);
            assert_eq!(selector.accepted, expected);
            assert_eq!(tree.overlapping(&boxes[i]).len(), expected);
        }
    }
}

/// All pairs of overlapping sub-shape boxes: 8·7 + 12·11 + 6·17 = 290.
#[test]
fn bug30655_pair_selection_typed_and_void() {
    let (_, _, boxes) = box_set();
    let tree = Bvh::new(&boxes);
    let mut selector = PairCounter::default();
    assert_eq!(tree.select_pairs(&tree, &mut selector), 290);
    assert_eq!(selector.accepted, 290);
    assert_eq!(tree.overlapping_pairs(&tree).len(), 290);
}
