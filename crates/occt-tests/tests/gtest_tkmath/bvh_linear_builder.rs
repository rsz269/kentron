//! BVH_LinearBuilder_Test.cxx (BVH_LinearBuilder: `math::Bvh`)
//!
//! See bvh_binned_builder.rs: the tests run against kentron's one builder.
//!
//! Not ported: `MortonCodeSorting` (Morton codes are the linear builder's algorithm),
//! `Build2D` (boxes are 3D), `CompareWithBinnedBuilder` (one builder).

use crate::bvh_support::{cube, row, total_in_leaves, tree};
use kentron::math::{BoundingBox, Bvh, BvhOptions};

#[test]
fn default_constructor() {
    let options = BvhOptions::default();
    assert!(options.leaf_size > 0);
    assert!(options.max_depth > 0);
}

#[test]
fn custom_parameters() {
    assert_eq!(tree(&[], 5, 20).options(), BvhOptions { leaf_size: 5, max_depth: 20 });
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
fn leaf_node_size_respected() {
    assert!(crate::bvh_support::leaves_at_most(&tree(&row(10), 3, 32), 3));
}

#[test]
fn large_data_set() {
    let mut boxes = Vec::new();
    for x in 0..10 {
        for y in 0..10 {
            for z in 0..10 {
                let (x, y, z) = (2.0 * x as f64, 2.0 * y as f64, 2.0 * z as f64);
                boxes.push(cube((x, y, z), (x + 1.0, y + 1.0, z + 1.0)));
            }
        }
    }
    let t = tree(&boxes, 4, 32);
    assert!(t.node_count() > 1);
    assert_eq!(total_in_leaves(&t), 1000);
}

#[test]
fn max_depth_parameter() {
    let depth = tree(&row(100), 1, 10).depth();
    assert!(depth > 0);
    assert!(depth < 20);
}

#[test]
fn clustered_data() {
    let mut boxes: Vec<BoundingBox> = Vec::new();
    for i in 0..10 {
        let d = 0.1 * i as f64;
        boxes.push(cube((d, d, d), (d + 0.1, d + 0.1, d + 0.1)));
        boxes.push(cube((100.0 + d, 100.0 + d, 100.0 + d), (100.1 + d, 100.1 + d, 100.1 + d)));
    }
    let t = tree(&boxes, 1, 32);
    assert!(t.node_count() > 1);
    assert!(t.sah_cost() > 0.0);
}

#[test]
fn overlapping_boxes() {
    let boxes: Vec<BoundingBox> = (0..10).map(|i| cube((0.5 * i as f64, 0.0, 0.0), (0.5 * i as f64 + 2.0, 1.0, 1.0))).collect();
    assert!(tree(&boxes, 1, 32).node_count() > 1);
}

#[test]
fn identical_boxes() {
    assert!(tree(&vec![cube((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)); 10], 1, 32).node_count() >= 1);
}

#[test]
fn negative_coordinates() {
    let boxes: Vec<BoundingBox> =
        (0..10).map(|i| cube((-10.0 + 2.0 * i as f64, -5.0, -5.0), (-9.0 + 2.0 * i as f64, -4.0, -4.0))).collect();
    assert!(tree(&boxes, 1, 32).node_count() > 1);
}
