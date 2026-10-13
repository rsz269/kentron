//! BVH_SpatialMedianBuilder_Test.cxx (BVH_SpatialMedianBuilder over BVH_Triangulation:
//! `math::Bvh` over the triangles' boxes)
//!
//! See bvh_binned_builder.rs: the tests run against kentron's one builder (whose default leaf
//! size is 4, OCCT's 5).
//!
//! Not ported: `UseMainAxisFlag` (an option of OCCT's builder), `CompareWithBinnedBuilder`
//! (one builder), `FloatPrecision` (no single-precision boxes).

use crate::bvh_support::{leaves_at_most, tree, triangle, triangles};
use kentron::math::{BoundingBox, Bvh, BvhOptions};

/// Unit triangles one unit apart.
fn strip(n: usize) -> Vec<BoundingBox> {
    triangles(n, 1.0, 1.0)
}

#[test]
fn default_constructor() {
    assert_eq!(BvhOptions::default(), BvhOptions { leaf_size: 4, max_depth: 32 });
}

#[test]
fn custom_parameters() {
    assert_eq!(tree(&[], 10, 20).options(), BvhOptions { leaf_size: 10, max_depth: 20 });
}

#[test]
fn build_empty_set() {
    assert_eq!(Bvh::new(&[]).node_count(), 0);
}

#[test]
fn build_single_triangle() {
    let t = Bvh::new(&[triangle((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))]);
    assert_eq!(t.node_count(), 1);
    let root = t.root().expect("a root");
    assert!(root.is_leaf());
    assert_eq!(root.len(), 1);
}

#[test]
fn build_multiple_triangles() {
    let t = Bvh::new(&strip(10));
    assert!(t.node_count() > 0);
    assert!(t.depth() >= 1);
}

#[test]
fn leaf_node_size_respected() {
    assert!(leaves_at_most(&tree(&strip(20), 3, 32), 3));
}

#[test]
fn max_depth_respected() {
    assert!(tree(&strip(100), 1, 5).depth() <= 5);
}

#[test]
fn splits_at_spatial_median() {
    let t = tree(&strip(8), 2, 32);
    let (left, right) = t.root().expect("a root").children().expect("an inner root");
    for side in [left, right] {
        assert!((2..=6).contains(&side.len()));
    }
}

#[test]
fn build2_d() {
    // Triangles in a plane.
    let t = Bvh::new(&strip(10));
    assert!(t.node_count() > 0);
    assert!(t.depth() >= 1);
}

#[test]
fn clustered_data() {
    let boxes: Vec<BoundingBox> = (0..30)
        .map(|i| {
            let (x, y) = (if i < 15 { 0.0 } else { 100.0 }, (i % 15) as f64);
            triangle((x, y, 0.0), (x + 0.5, y + 1.0, 0.0), (x + 1.0, y, 0.0))
        })
        .collect();
    let t = Bvh::new(&boxes);
    assert!(t.node_count() > 0);
    assert!(t.depth() >= 1);
}

#[test]
fn identical_boxes() {
    let boxes = vec![triangle((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)); 10];
    assert!(Bvh::new(&boxes).node_count() > 0);
}

#[test]
fn large_data_set() {
    let boxes: Vec<BoundingBox> = (0..500)
        .map(|i| {
            let (x, y, z) = (((i * 7) % 100) as f64, ((i * 11) % 100) as f64, ((i * 13) % 100) as f64);
            triangle((x, y, z), (x + 0.5, y + 1.0, z), (x + 1.0, y, z + 1.0))
        })
        .collect();
    let t = Bvh::new(&boxes);
    assert!(t.node_count() > 0);
    assert!(t.depth() > 5);
    assert!(t.sah_cost() > 0.0);
}
