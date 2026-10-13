//! BVH_SweepPlaneBuilder_Test.cxx (BVH_SweepPlaneBuilder over BVH_Triangulation: `math::Bvh`
//! over the triangles' boxes)
//!
//! See bvh_binned_builder.rs: the tests run against kentron's one builder.

use crate::bvh_support::{leaves_at_most, tree, triangle};
use kentron::math::{BoundingBox, Bvh, BvhOptions};

#[test]
fn default_constructor() {
    let options = BvhOptions::default();
    assert!(options.leaf_size > 0);
    assert_eq!(options.max_depth, 32);
}

#[test]
fn custom_parameters() {
    assert_eq!(tree(&[], 5, 20).options(), BvhOptions { leaf_size: 5, max_depth: 20 });
}

#[test]
fn build_empty_triangulation() {
    assert_eq!(Bvh::new(&[]).node_count(), 0);
}

#[test]
fn build_single_triangle() {
    let t = tree(&[triangle((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))], 1, 32);
    assert_eq!(t.node_count(), 1);
    let root = t.root().expect("a root");
    assert!(root.is_leaf());
    assert_eq!(root.boxes(), &[0]);
}

#[test]
fn build_two_triangles() {
    // The two halves of a unit square: the same box.
    let boxes = [triangle((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)), triangle((0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0))];
    let t = tree(&boxes, 1, 32);
    assert!(t.node_count() > 2);
    assert!(t.depth() >= 1);
}

#[test]
fn build_multiple_triangles_grid() {
    let p = |i: usize, j: usize| (i as f64, j as f64, 0.0);
    let mut boxes = Vec::new();
    for i in 0..2 {
        for j in 0..2 {
            boxes.push(triangle(p(i, j), p(i, j + 1), p(i + 1, j + 1)));
            boxes.push(triangle(p(i, j), p(i + 1, j + 1), p(i + 1, j)));
        }
    }
    let t = tree(&boxes, 2, 32);
    assert!(t.node_count() > 0);
    assert!(t.depth() >= 1);
    let root = t.bbox();
    assert!(boxes.iter().all(|b| root.overlaps(b)));
}

#[test]
fn split_along_different_axes() {
    let boxes: Vec<BoundingBox> = (0..10)
        .map(|i| {
            let x = 2.0 * i as f64;
            triangle((x, 0.0, 0.0), (x + 1.0, 0.0, 0.0), (x + 0.5, 1.0, 0.0))
        })
        .collect();
    let t = tree(&boxes, 1, 32);
    assert!(t.node_count() > 1);
    assert!(t.depth() >= 1);
}

#[test]
fn degenerate_case_all_same_position() {
    let boxes = vec![triangle((0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0)); 5];
    assert!(tree(&boxes, 1, 32).node_count() > 0);
}

#[test]
fn sah_quality_verify_better() {
    // Two clusters: the root's children should barely overlap.
    let cluster = |x0: f64| (0..5).map(move |i| {
        let x = x0 + 0.5 * i as f64;
        triangle((x, 0.0, 0.0), (x + 0.4, 0.0, 0.0), (x + 0.2, 0.4, 0.0))
    });
    let boxes: Vec<BoundingBox> = cluster(0.0).chain(cluster(10.0)).collect();
    let t = tree(&boxes, 2, 32);
    assert!(t.depth() >= 1);
    let (left, right) = t.root().expect("a root").children().expect("an inner root");
    let (a, b) = (left.bbox(), right.bbox());
    let overlap: f64 = [(a.min.x, a.max.x, b.min.x, b.max.x), (a.min.y, a.max.y, b.min.y, b.max.y), (a.min.z, a.max.z, b.min.z, b.max.z)]
        .iter()
        .map(|&(l1, h1, l2, h2)| (h1.min(h2) - l1.max(l2)).max(0.0))
        .sum();
    assert!(overlap < 0.1 * t.bbox().size().x);
}

#[test]
fn leaf_size_respect_max_size() {
    let boxes: Vec<BoundingBox> = (0..20)
        .map(|i| {
            let x = i as f64;
            triangle((x, 0.0, 0.0), (x + 0.9, 0.0, 0.0), (x + 0.5, 0.9, 0.0))
        })
        .collect();
    assert!(leaves_at_most(&tree(&boxes, 5, 32), 5));
}
