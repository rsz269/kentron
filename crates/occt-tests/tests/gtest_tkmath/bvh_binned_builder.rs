//! BVH_BinnedBuilder_Test.cxx (BVH_BinnedBuilder, with BVH_BoxSet: `math::Bvh`)
//!
//! kentron has one builder (a median split along the longest side of the boxes' centres);
//! the builder tests run against it with OCCT's leaf sizes and depth limits (`BvhOptions`).
//! `Length` is `node_count`, leaves' `NbPrimitives` are `leaf_sizes`, `EstimateSAH` is
//! `sah_cost`.
//!
//! Not ported: `BuildWithDifferentBinCounts` (binning is OCCT's algorithm), `Build2D` (boxes
//! are 3D; a flat set is tested by the spatial median builder's `build2_d`).

use crate::bvh_support::{cube, leaves_at_most, row, total_in_leaves, tree};
use kentron::math::{BoundingBox, Bvh, BvhOptions};

#[test]
fn default_constructor() {
    let options = BvhOptions::default();
    assert!(options.leaf_size > 0);
    assert!(options.max_depth > 0);
}

#[test]
fn custom_parameters() {
    let t = tree(&[], 5, 20);
    assert_eq!(t.options(), BvhOptions { leaf_size: 5, max_depth: 20 });
}

#[test]
fn build_empty_set() {
    assert_eq!(Bvh::new(&[]).node_count(), 0);
}

#[test]
fn build_single_element() {
    let t = Bvh::new(&[cube((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))]);
    assert_eq!(t.node_count(), 1);
    assert!(t.root().expect("a root").is_leaf());
}

#[test]
fn build_multiple_elements() {
    let t = tree(&row(10), 1, 32);
    assert!(t.node_count() > 1);
    assert!(t.depth() > 0);
}

#[test]
fn sah_optimization() {
    // Two clusters far apart.
    let mut boxes: Vec<BoundingBox> = (0..5).map(|i| cube((0.1 * i as f64, 0.0, 0.0), (0.1 * i as f64 + 0.1, 0.1, 0.1))).collect();
    boxes.extend((0..5).map(|i| cube((100.0 + 0.1 * i as f64, 0.0, 0.0), (100.1 + 0.1 * i as f64, 0.1, 0.1))));
    assert!(tree(&boxes, 1, 32).sah_cost() > 0.0);
}

#[test]
fn leaf_node_size_respected() {
    assert!(leaves_at_most(&tree(&row(10), 3, 32), 3));
}

#[test]
fn random_distribution() {
    let mut boxes = Vec::new();
    for x in 0..5 {
        for y in 0..5 {
            for z in 0..5 {
                let (x, y, z) = (2.0 * x as f64, 2.0 * y as f64, 2.0 * z as f64);
                boxes.push(cube((x, y, z), (x + 1.0, y + 1.0, z + 1.0)));
            }
        }
    }
    let t = tree(&boxes, 1, 32);
    assert!(t.node_count() > 1);
    assert_eq!(total_in_leaves(&t), 125);
}

#[test]
fn compare_tree_quality() {
    // Smaller leaves: a deeper tree.
    assert!(tree(&row(50), 1, 32).depth() >= tree(&row(50), 4, 32).depth());
}

#[test]
fn max_depth_respected() {
    assert!(tree(&row(100), 1, 5).depth() <= 5);
}

#[test]
fn overlapping_boxes() {
    let boxes: Vec<BoundingBox> = (0..10).map(|i| cube((0.5 * i as f64, 0.0, 0.0), (0.5 * i as f64 + 2.0, 1.0, 1.0))).collect();
    assert!(tree(&boxes, 1, 32).node_count() > 1);
}

#[test]
fn identical_boxes() {
    let boxes = vec![cube((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)); 10];
    assert!(tree(&boxes, 1, 32).node_count() >= 1);
}
