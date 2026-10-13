//! BVH_Traverse_Test.cxx (BVH_Traverse, BVH_PairTraverse: `math::Selector` and `PairSelector`
//! run by `Bvh::select` and `Bvh::select_pairs`)
//!
//! A selector's `RejectNode` is `enter` (returning the metric, or `None` to reject),
//! `IsMetricBetter`, `RejectMetric` and `Stop` keep their meaning. Unlike OCCT's, kentron's
//! traversal also offers each box of a leaf to `enter` before `accept`.
//!
//! Not ported: `NullTree`, `NullTrees` (there are no null trees; empty ones are tested).

use kentron::math::{BoundingBox, Bvh, PairSelector, Point, Selector};

/// The boxes of `n` triangles along X, as OCCT's `CreateSimpleTriangulationBVH`: triangle `i`
/// has corners (2i, 0, 0), (2i + 1, 1, 0), (2i + 2, 0, 0).
fn triangles(n: usize, x0: f64, step: f64, width: f64) -> Vec<BoundingBox> {
    (0..n)
        .map(|i| {
            let x = x0 + step * i as f64;
            BoundingBox::from_points([Point::new(x, 0.0, 0.0), Point::new(x + 0.5 * width, 1.0, 0.0), Point::new(x + width, 0.0, 0.0)])
        })
        .collect()
}

fn simple_tree(n: usize) -> Bvh {
    Bvh::new(&triangles(n, 0.0, 2.0, 2.0))
}

/// Three triangles from `x0`, one unit apart.
fn three_triangles(x0: f64) -> Vec<BoundingBox> {
    triangles(3, x0, 1.0, 1.0)
}

#[derive(Default)]
struct CountAll {
    accepted: usize,
}

impl Selector for CountAll {
    fn enter(&mut self, _: &BoundingBox) -> Option<f64> {
        Some(0.0)
    }
    fn accept(&mut self, _: usize, _: f64) -> bool {
        self.accepted += 1;
        true
    }
}

struct BoxSelector {
    query: BoundingBox,
    accepted: usize,
}

impl Selector for BoxSelector {
    fn enter(&mut self, b: &BoundingBox) -> Option<f64> {
        self.query.overlaps(b).then_some(0.0)
    }
    fn accept(&mut self, _: usize, _: f64) -> bool {
        self.accepted += 1;
        true
    }
}

/// Boxes within `max` of a point, the closest remembered.
struct DistanceSelector {
    point: Point,
    max: f64,
    min: f64,
    closest: Option<usize>,
    accepted: usize,
}

impl DistanceSelector {
    fn new(point: Point, max: f64) -> Self {
        DistanceSelector { point, max, min: f64::INFINITY, closest: None, accepted: 0 }
    }
}

impl Selector for DistanceSelector {
    fn enter(&mut self, b: &BoundingBox) -> Option<f64> {
        let d = b.distance_to_point(self.point);
        (d <= self.max).then_some(d)
    }
    fn accept(&mut self, index: usize, metric: f64) -> bool {
        self.accepted += 1;
        if metric < self.min {
            self.min = metric;
            self.closest = Some(index);
        }
        true
    }
    fn reject_metric(&self, metric: f64) -> bool {
        metric > self.max
    }
}

struct Limited {
    max: usize,
    accepted: usize,
}

impl Selector for Limited {
    fn enter(&mut self, _: &BoundingBox) -> Option<f64> {
        Some(0.0)
    }
    fn accept(&mut self, _: usize, _: f64) -> bool {
        self.accepted += 1;
        true
    }
    fn stop(&self) -> bool {
        self.accepted >= self.max
    }
}

#[derive(Default)]
struct CountAllPairs {
    accepted: usize,
}

impl PairSelector for CountAllPairs {
    fn enter(&mut self, _: &BoundingBox, _: &BoundingBox) -> bool {
        true
    }
    fn accept(&mut self, _: usize, _: usize) -> bool {
        self.accepted += 1;
        true
    }
}

#[derive(Default)]
struct OverlapDetector {
    overlaps: usize,
    rejects: usize,
}

impl PairSelector for OverlapDetector {
    fn enter(&mut self, a: &BoundingBox, b: &BoundingBox) -> bool {
        self.rejects += 1;
        a.overlaps(b)
    }
    fn accept(&mut self, _: usize, _: usize) -> bool {
        self.overlaps += 1;
        true
    }
}

#[test]
fn count_all_elements() {
    let mut selector = CountAll::default();
    assert_eq!(simple_tree(10).select(&mut selector), 10);
    assert_eq!(selector.accepted, 10);
}

#[test]
fn empty_tree() {
    assert_eq!(Bvh::new(&[]).select(&mut CountAll::default()), 0);
}

#[test]
fn box_selection() {
    let query = BoundingBox::new(Point::new(0.0, 0.0, -1.0), Point::new(10.0, 2.0, 1.0));
    let count = simple_tree(10).select(&mut BoxSelector { query, accepted: 0 });
    assert!(count > 0);
    assert!(count < 10);
}

#[test]
fn empty_box_selection() {
    let query = BoundingBox::new(Point::new(100.0, 100.0, 100.0), Point::new(200.0, 200.0, 200.0));
    assert_eq!(simple_tree(10).select(&mut BoxSelector { query, accepted: 0 }), 0);
}

#[test]
fn full_box_selection() {
    let query = BoundingBox::new(Point::new(-100.0, -100.0, -100.0), Point::new(100.0, 100.0, 100.0));
    assert_eq!(simple_tree(10).select(&mut BoxSelector { query, accepted: 0 }), 10);
}

#[test]
fn distance_based_selection() {
    let mut selector = DistanceSelector::new(Point::new(5.0, 0.5, 0.0), 5.0);
    let count = simple_tree(10).select(&mut selector);
    assert!(count > 0);
    assert!(count <= 10);
    assert!(selector.closest.is_some_and(|i| i < 10));
    assert_eq!(selector.accepted, count);
}

#[test]
fn early_termination() {
    let mut selector = Limited { max: 5, accepted: 0 };
    assert_eq!(simple_tree(100).select(&mut selector), 5);
    assert_eq!(selector.accepted, 5);
}

#[test]
fn large_data_set() {
    assert_eq!(simple_tree(1000).select(&mut CountAll::default()), 1000);
}

#[test]
fn metric_based_pruning() {
    let mut selector = DistanceSelector::new(Point::new(1000.0, 1000.0, 1000.0), 1.0);
    assert_eq!(simple_tree(50).select(&mut selector), 0);
}

#[test]
fn is_out_verification() {
    let b1 = BoundingBox::new(Point::new(0.0, 0.0, 0.0), Point::new(3.0, 1.0, 0.0));
    let b2 = BoundingBox::new(Point::new(100.0, 0.0, 0.0), Point::new(103.0, 1.0, 0.0));
    assert!(!b1.overlaps(&b2));
    assert!(!b2.overlaps(&b1));
}

#[test]
fn triangulation_box_verification() {
    let union = |boxes: &[BoundingBox]| boxes.iter().fold(BoundingBox::EMPTY, |mut u, b| {
        u.add_box(b);
        u
    });
    let (b1, b2) = (union(&three_triangles(0.0)), union(&three_triangles(100.0)));
    assert!(!b1.is_empty());
    assert!(!b2.is_empty());
    assert!((b1.min.x - 0.0).abs() <= 0.01 && (b1.max.x - 3.0).abs() <= 0.01);
    assert!((b2.min.x - 100.0).abs() <= 0.01 && (b2.max.x - 103.0).abs() <= 0.01);
    assert!(!b1.overlaps(&b2));
    assert!(!b2.overlaps(&b1));
}

#[test]
fn bvh_root_box_verification() {
    let (t1, t2) = (Bvh::new(&three_triangles(0.0)), Bvh::new(&three_triangles(100.0)));
    let (b1, b2) = (t1.bbox(), t2.bbox());
    assert!((b1.min.x - 0.0).abs() <= 0.01 && (b1.max.x - 3.0).abs() <= 0.01);
    assert!((b2.min.x - 100.0).abs() <= 0.01 && (b2.max.x - 103.0).abs() <= 0.01);
    assert!(!b1.overlaps(&b2));
    // Three boxes fit in one leaf: the root is a leaf.
    assert_eq!(t1.node_count(), 1);
    assert_eq!(t2.node_count(), 1);
}

#[test]
fn count_all_pairs() {
    assert_eq!(simple_tree(5).select_pairs(&simple_tree(5), &mut CountAllPairs::default()), 25);
}

#[test]
fn empty_first_tree() {
    assert_eq!(Bvh::new(&[]).select_pairs(&simple_tree(5), &mut CountAllPairs::default()), 0);
}

#[test]
fn empty_second_tree() {
    assert_eq!(simple_tree(5).select_pairs(&Bvh::new(&[]), &mut CountAllPairs::default()), 0);
}

#[test]
fn overlap_detection_same_trees() {
    let tree = simple_tree(10);
    assert!(tree.select_pairs(&tree, &mut OverlapDetector::default()) >= 10);
}

#[test]
fn overlap_detection_non_overlapping() {
    let (t1, t2) = (Bvh::new(&three_triangles(0.0)), Bvh::new(&three_triangles(100.0)));
    let mut selector = OverlapDetector::default();
    let count = t1.select_pairs(&t2, &mut selector);
    assert!(selector.rejects > 0);
    assert_eq!(count, 0);
}

#[test]
fn asymmetric_pairs() {
    assert_eq!(simple_tree(3).select_pairs(&simple_tree(7), &mut CountAllPairs::default()), 21);
}

#[test]
fn large_data_sets() {
    assert_eq!(simple_tree(50).select_pairs(&simple_tree(50), &mut CountAllPairs::default()), 2500);
}
