//! Shared helpers of the BVH test modules.

use kentron::math::{BoundingBox, Bvh, BvhOptions, Point};

/// The box of a triangle.
pub fn triangle(a: (f64, f64, f64), b: (f64, f64, f64), c: (f64, f64, f64)) -> BoundingBox {
    BoundingBox::from_points([Point::new(a.0, a.1, a.2), Point::new(b.0, b.1, b.2), Point::new(c.0, c.1, c.2)])
}

/// An axis-aligned box from its corners.
pub fn cube(min: (f64, f64, f64), max: (f64, f64, f64)) -> BoundingBox {
    BoundingBox::new(Point::new(min.0, min.1, min.2), Point::new(max.0, max.1, max.2))
}

/// A tree with the given leaf size and depth limit (OCCT's `Builder(leafSize, maxDepth)`).
pub fn tree(boxes: &[BoundingBox], leaf_size: usize, max_depth: usize) -> Bvh {
    Bvh::with_options(boxes, BvhOptions { leaf_size, max_depth })
}

/// `n` unit boxes along X, one every 2 units (the builders' tests' common set).
pub fn row(n: usize) -> Vec<BoundingBox> {
    (0..n).map(|i| cube((2.0 * i as f64, 0.0, 0.0), (2.0 * i as f64 + 1.0, 1.0, 1.0))).collect()
}

/// `n` triangles along X from (x, 0, 0), (x + w/2, 1, 0), (x + w, 0, 0), x stepping by `step`.
pub fn triangles(n: usize, step: f64, width: f64) -> Vec<BoundingBox> {
    (0..n)
        .map(|i| {
            let x = step * i as f64;
            triangle((x, 0.0, 0.0), (x + 0.5 * width, 1.0, 0.0), (x + width, 0.0, 0.0))
        })
        .collect()
}

/// Every leaf respects the leaf size.
pub fn leaves_at_most(tree: &Bvh, size: usize) -> bool {
    tree.leaf_sizes().into_iter().all(|n| n <= size)
}

/// The leaves hold every box once.
pub fn total_in_leaves(tree: &Bvh) -> usize {
    tree.leaf_sizes().into_iter().sum()
}
