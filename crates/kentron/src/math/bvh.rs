//! A bounding volume hierarchy: a binary tree of axis-aligned boxes over a set of boxes, for
//! finding those overlapping a box, containing a point or crossed by a ray, and the
//! overlapping pairs of two sets, without testing every one.

use super::{BoundingBox, Point, Vector};

/// How a [`Bvh`] is built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BvhOptions {
    /// The largest number of boxes in a leaf (unless the depth limit is reached).
    #[doc(alias = "LeafNodeSize")]
    pub leaf_size: usize,
    /// The largest depth of a node (the root's is 0).
    #[doc(alias = "MaxTreeDepth")]
    pub max_depth: usize,
}

impl Default for BvhOptions {
    fn default() -> Self {
        BvhOptions { leaf_size: 4, max_depth: 32 }
    }
}

/// A query run over a [`Bvh`] by [`Bvh::select`]: it decides which nodes to enter (with a
/// metric ordering them, such as a distance) and which boxes to accept.
#[doc(alias = "BVH_Traverse")]
pub trait Selector {
    /// The metric of a node, or of a box before it is offered to `accept`; `None` skips it.
    #[doc(alias = "RejectNode")]
    fn enter(&mut self, bbox: &BoundingBox) -> Option<f64>;

    /// Offers box `index` (with its metric); returns whether it is accepted (counted).
    fn accept(&mut self, index: usize, metric: f64) -> bool;

    /// Whether metric `a` is better than `b`: of two children, the better is visited first.
    #[doc(alias = "IsMetricBetter")]
    fn better(&self, a: f64, b: f64) -> bool {
        a < b
    }

    /// Whether a node waiting to be visited is no longer worth it (its metric is worse than
    /// what was found since).
    #[doc(alias = "RejectMetric")]
    fn reject_metric(&self, _metric: f64) -> bool {
        false
    }

    /// Whether to stop the traversal.
    fn stop(&self) -> bool {
        false
    }
}

/// A query run over two trees by [`Bvh::select_pairs`].
#[doc(alias = "BVH_PairTraverse")]
pub trait PairSelector {
    /// Whether to look into a pair of nodes, or of boxes before offering them to `accept`.
    #[doc(alias = "RejectNode")]
    fn enter(&mut self, a: &BoundingBox, b: &BoundingBox) -> bool;

    /// Offers the pair of boxes `(i, j)`; returns whether it is accepted (counted).
    fn accept(&mut self, i: usize, j: usize) -> bool;

    fn stop(&self) -> bool {
        false
    }
}

#[derive(Debug, Clone)]
struct Node {
    bbox: BoundingBox,
    /// A leaf holds `order[first..first + count]`; an inner node (`count == 0`) has its
    /// children at `first` and `first + 1`.
    first: usize,
    count: usize,
}

/// A tree over boxes given by index; queries return those indices.
#[doc(alias = "BVH_Tree")]
#[doc(alias = "BVH_BinnedBuilder")]
#[derive(Debug, Clone, Default)]
pub struct Bvh {
    nodes: Vec<Node>,
    order: Vec<usize>,
    boxes: Vec<BoundingBox>,
    options: BvhOptions,
}

impl Bvh {
    /// The tree over `boxes`, split at the median of the boxes' centres along the longest
    /// side. Empty boxes are left out of every query.
    pub fn new(boxes: &[BoundingBox]) -> Self {
        Self::with_options(boxes, BvhOptions::default())
    }

    /// [`Bvh::new`] with other leaf size and depth limits.
    pub fn with_options(boxes: &[BoundingBox], options: BvhOptions) -> Self {
        let options = BvhOptions { leaf_size: options.leaf_size.max(1), ..options };
        let mut tree = Bvh { nodes: Vec::new(), order: Vec::new(), boxes: boxes.to_vec(), options };
        let mut items: Vec<usize> = (0..boxes.len()).filter(|&i| !boxes[i].is_empty()).collect();
        if !items.is_empty() {
            tree.nodes.push(Node { bbox: BoundingBox::EMPTY, first: 0, count: 0 });
            tree.build(0, &mut items, 0, 0);
            tree.order = items;
        }
        tree
    }

    /// The centre of a box, kept finite for boxes with infinite sides.
    fn centre(b: &BoundingBox) -> Point {
        let mid = |lo: f64, hi: f64| {
            let m = 0.5 * (lo + hi);
            if m.is_nan() { 0.0 } else { m.clamp(-1e300, 1e300) }
        };
        Point::new(mid(b.min.x, b.max.x), mid(b.min.y, b.max.y), mid(b.min.z, b.max.z))
    }

    /// Fills node `at` with `items`, the part of the final order starting at `offset`
    /// (sorted in place, in sub-slices for the children).
    fn build(&mut self, at: usize, items: &mut [usize], offset: usize, depth: usize) {
        let mut bbox = BoundingBox::EMPTY;
        let mut centres = BoundingBox::EMPTY;
        for &i in items.iter() {
            bbox.add_box(&self.boxes[i]);
            centres.add_point(Self::centre(&self.boxes[i]));
        }
        self.nodes[at].bbox = bbox;
        let extent = centres.max - centres.min;
        if items.len() <= self.options.leaf_size || depth >= self.options.max_depth {
            self.nodes[at].first = offset;
            self.nodes[at].count = items.len();
            return;
        }
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let half = items.len() / 2;
        // Boxes all centred alike are split in halves as they come.
        if extent.length() > 0.0 {
            // (The keys worked out once, moved with their items.)
            let mut keyed: Vec<(f64, usize)> = items
                .iter()
                .map(|&i| {
                    let c = Self::centre(&self.boxes[i]);
                    ([c.x, c.y, c.z][axis], i)
                })
                .collect();
            keyed.select_nth_unstable_by(half, |a, b| a.0.total_cmp(&b.0));
            for (item, (_, i)) in items.iter_mut().zip(keyed) {
                *item = i;
            }
        }
        let children = self.nodes.len();
        self.nodes.push(Node { bbox: BoundingBox::EMPTY, first: 0, count: 0 });
        self.nodes.push(Node { bbox: BoundingBox::EMPTY, first: 0, count: 0 });
        self.nodes[at].first = children;
        self.nodes[at].count = 0;
        let (left, right) = items.split_at_mut(half);
        self.build(children, left, offset, depth + 1);
        self.build(children + 1, right, offset + half, depth + 1);
    }

    /// The number of boxes in the tree (empty boxes are left out).
    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// The number of nodes.
    #[doc(alias = "Length")]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// The depth of the deepest leaf (0 for a root alone or an empty tree).
    pub fn depth(&self) -> usize {
        let mut deepest = 0;
        let mut stack = if self.nodes.is_empty() { Vec::new() } else { vec![(0usize, 0usize)] };
        while let Some((n, d)) = stack.pop() {
            deepest = deepest.max(d);
            if self.nodes[n].count == 0 {
                stack.extend([(self.nodes[n].first, d + 1), (self.nodes[n].first + 1, d + 1)]);
            }
        }
        deepest
    }

    /// The number of boxes in each leaf.
    #[doc(alias = "NbPrimitives")]
    pub fn leaf_sizes(&self) -> Vec<usize> {
        self.nodes.iter().filter(|n| n.count > 0).map(|n| n.count).collect()
    }

    /// The box around every box of the tree.
    pub fn bbox(&self) -> BoundingBox {
        self.nodes.first().map_or(BoundingBox::EMPTY, |n| n.bbox)
    }

    pub fn options(&self) -> BvhOptions {
        self.options
    }

    /// The root node, if the tree is not empty.
    pub fn root(&self) -> Option<BvhNode<'_>> {
        (!self.nodes.is_empty()).then_some(BvhNode { tree: self, index: 0 })
    }

    /// The tree's cost by the surface area heuristic: the expected number of node visits and
    /// box tests for a random ray through the root's box (each node weighted by its area
    /// relative to the root's: 1 for visiting it, and 1 per box of a leaf).
    #[doc(alias = "EstimateSAH")]
    pub fn sah_cost(&self) -> f64 {
        let Some(root) = self.nodes.first() else { return 0.0 };
        let total = root.bbox.area();
        if total <= 0.0 {
            return self.nodes.len() as f64 + self.len() as f64;
        }
        self.nodes.iter().map(|n| n.bbox.area() / total * (1.0 + n.count as f64)).sum()
    }

    /// Runs `selector` over the tree, visiting better nodes first; returns the number of
    /// boxes it accepted.
    #[doc(alias = "Select")]
    pub fn select(&self, selector: &mut impl Selector) -> usize {
        let mut accepted = 0;
        if self.nodes.is_empty() {
            return 0;
        }
        let Some(root) = selector.enter(&self.nodes[0].bbox) else { return 0 };
        let mut stack = vec![(0usize, root)];
        while let Some((n, metric)) = stack.pop() {
            if selector.reject_metric(metric) {
                continue;
            }
            let node = &self.nodes[n];
            if node.count > 0 {
                for &i in &self.order[node.first..node.first + node.count] {
                    if let Some(m) = selector.enter(&self.boxes[i]) {
                        if selector.accept(i, m) {
                            accepted += 1;
                        }
                    }
                    if selector.stop() {
                        return accepted;
                    }
                }
                continue;
            }
            let (l, r) = (node.first, node.first + 1);
            let (ml, mr) = (selector.enter(&self.nodes[l].bbox), selector.enter(&self.nodes[r].bbox));
            if selector.stop() {
                return accepted;
            }
            // The better child on top of the stack, visited first.
            match (ml, mr) {
                (Some(a), Some(b)) if selector.better(a, b) => stack.extend([(r, b), (l, a)]),
                (Some(a), Some(b)) => stack.extend([(l, a), (r, b)]),
                (Some(a), None) => stack.push((l, a)),
                (None, Some(b)) => stack.push((r, b)),
                (None, None) => {}
            }
        }
        accepted
    }

    /// Runs `selector` over the pairs of nodes of `self` and `other`; returns the number of
    /// pairs of boxes it accepted.
    #[doc(alias = "BVH_PairTraverse::Select")]
    pub fn select_pairs(&self, other: &Bvh, selector: &mut impl PairSelector) -> usize {
        let mut accepted = 0;
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return 0;
        }
        let mut stack = vec![(0usize, 0usize)];
        while let Some((a, b)) = stack.pop() {
            let (na, nb) = (&self.nodes[a], &other.nodes[b]);
            if !selector.enter(&na.bbox, &nb.bbox) {
                continue;
            }
            match (na.count > 0, nb.count > 0) {
                (true, true) => {
                    for &i in &self.order[na.first..na.first + na.count] {
                        for &j in &other.order[nb.first..nb.first + nb.count] {
                            if selector.enter(&self.boxes[i], &other.boxes[j]) && selector.accept(i, j) {
                                accepted += 1;
                            }
                            if selector.stop() {
                                return accepted;
                            }
                        }
                    }
                }
                (false, true) => stack.extend([(na.first, b), (na.first + 1, b)]),
                (true, false) => stack.extend([(a, nb.first), (a, nb.first + 1)]),
                (false, false) => {
                    let size = |n: &Node| (n.bbox.max - n.bbox.min).length();
                    if size(na) >= size(nb) {
                        stack.extend([(na.first, b), (na.first + 1, b)]);
                    } else {
                        stack.extend([(a, nb.first), (a, nb.first + 1)]);
                    }
                }
            }
            if selector.stop() {
                return accepted;
            }
        }
        accepted
    }

    /// The box closest to `p` by `distance` (of a box index to `p`, never less than the
    /// distance to its box), and that distance; `None` for an empty tree.
    #[doc(alias = "BVH_Distance")]
    pub fn nearest(&self, p: Point, mut distance: impl FnMut(usize) -> f64) -> Option<(usize, f64)> {
        struct Nearest<'a, F> {
            p: Point,
            distance: &'a mut F,
            best: Option<(usize, f64)>,
        }
        impl<F: FnMut(usize) -> f64> Selector for Nearest<'_, F> {
            fn enter(&mut self, bbox: &BoundingBox) -> Option<f64> {
                let d = bbox.distance_to_point(self.p);
                self.best.is_none_or(|(_, best)| d < best).then_some(d)
            }
            fn accept(&mut self, index: usize, _: f64) -> bool {
                let d = (self.distance)(index);
                let better = self.best.is_none_or(|(_, best)| d < best);
                if better {
                    self.best = Some((index, d));
                }
                better
            }
            fn reject_metric(&self, metric: f64) -> bool {
                self.best.is_some_and(|(_, best)| metric >= best)
            }
        }
        let mut search = Nearest { p, distance: &mut distance, best: None };
        self.select(&mut search);
        search.best
    }

    /// The closest pair of a box of `self` and a box of `other` by `distance` (between the
    /// two indices, never less than the distance between their boxes), and that distance;
    /// `None` if a tree is empty. Pairs of nodes are visited nearest first and skipped once
    /// no closer than the best pair found.
    #[doc(alias = "BVH_PairDistance")]
    pub fn closest_pair(&self, other: &Bvh, mut distance: impl FnMut(usize, usize) -> f64) -> Option<(usize, usize, f64)> {
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return None;
        }
        let mut best: Option<(usize, usize, f64)> = None;
        // Pairs of nodes by increasing box distance.
        let mut queue = std::collections::BinaryHeap::new();
        let gap = |a: usize, b: usize| self.nodes[a].bbox.distance_to_box(&other.nodes[b].bbox);
        queue.push(Pending(gap(0, 0), 0, 0));
        while let Some(Pending(d, a, b)) = queue.pop() {
            if best.is_some_and(|(_, _, found)| d >= found) {
                break;
            }
            let (na, nb) = (&self.nodes[a], &other.nodes[b]);
            match (na.count > 0, nb.count > 0) {
                (true, true) => {
                    for &i in &self.order[na.first..na.first + na.count] {
                        for &j in &other.order[nb.first..nb.first + nb.count] {
                            if best.is_some_and(|(_, _, found)| self.boxes[i].distance_to_box(&other.boxes[j]) >= found) {
                                continue;
                            }
                            let d = distance(i, j);
                            if best.is_none_or(|(_, _, found)| d < found) {
                                best = Some((i, j, d));
                            }
                        }
                    }
                }
                _ => {
                    let split_a = nb.count > 0 || (na.count == 0 && (na.bbox.max - na.bbox.min).length() >= (nb.bbox.max - nb.bbox.min).length());
                    let pairs = if split_a { [(na.first, b), (na.first + 1, b)] } else { [(a, nb.first), (a, nb.first + 1)] };
                    for (x, y) in pairs {
                        queue.push(Pending(gap(x, y), x, y));
                    }
                }
            }
        }
        best
    }

    /// Calls `visit` with the index of every box `accept`ed, descending into the nodes whose
    /// box `enter` accepts.
    fn walk(&self, enter: impl Fn(&BoundingBox) -> bool, mut visit: impl FnMut(usize)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n];
            if !enter(&node.bbox) {
                continue;
            }
            if node.count > 0 {
                for &i in &self.order[node.first..node.first + node.count] {
                    if enter(&self.boxes[i]) {
                        visit(i);
                    }
                }
            } else {
                stack.push(node.first);
                stack.push(node.first + 1);
            }
        }
    }

    /// The boxes overlapping `query` (touching counts).
    pub fn overlapping(&self, query: &BoundingBox) -> Vec<usize> {
        let mut out = Vec::new();
        self.walk(|b| overlap(b, query), |i| out.push(i));
        out
    }

    /// The boxes containing `p`.
    pub fn containing(&self, p: Point) -> Vec<usize> {
        let mut out = Vec::new();
        self.walk(|b| b.contains(p), |i| out.push(i));
        out
    }

    /// The boxes crossed by the ray from `origin` along `direction` (forward only).
    pub fn along_ray(&self, origin: Point, direction: Vector) -> Vec<usize> {
        let mut out = Vec::new();
        self.walk(|b| ray_meets(b, origin, direction), |i| out.push(i));
        out
    }

    /// The pairs `(i, j)` of a box of `self` and a box of `other` that overlap.
    pub fn overlapping_pairs(&self, other: &Bvh) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return out;
        }
        let mut stack = vec![(0usize, 0usize)];
        while let Some((a, b)) = stack.pop() {
            let (na, nb) = (&self.nodes[a], &other.nodes[b]);
            if !overlap(&na.bbox, &nb.bbox) {
                continue;
            }
            match (na.count > 0, nb.count > 0) {
                (true, true) => {
                    for &i in &self.order[na.first..na.first + na.count] {
                        for &j in &other.order[nb.first..nb.first + nb.count] {
                            if overlap(&self.boxes[i], &other.boxes[j]) {
                                out.push((i, j));
                            }
                        }
                    }
                }
                // Descend into the inner node (the larger one when both are).
                (false, true) => stack.extend([(na.first, b), (na.first + 1, b)]),
                (true, false) => stack.extend([(a, nb.first), (a, nb.first + 1)]),
                (false, false) => {
                    let size = |n: &Node| (n.bbox.max - n.bbox.min).length();
                    if size(na) >= size(nb) {
                        stack.extend([(na.first, b), (na.first + 1, b)]);
                    } else {
                        stack.extend([(a, nb.first), (a, nb.first + 1)]);
                    }
                }
            }
        }
        out
    }
}

impl Bvh {
    /// Pairs `(i, j)` of a box of this tree and one of `other` at most `reach` apart
    /// (overlapping ones included): both trees walked together, pairs of nodes farther apart
    /// passed over.
    pub fn pairs_within(&self, other: &Bvh, reach: f64) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return out;
        }
        let near = |a: &BoundingBox, b: &BoundingBox| a.distance_to_box(b) <= reach;
        let mut stack = vec![(0usize, 0usize)];
        while let Some((a, b)) = stack.pop() {
            let (na, nb) = (&self.nodes[a], &other.nodes[b]);
            if !near(&na.bbox, &nb.bbox) {
                continue;
            }
            match (na.count > 0, nb.count > 0) {
                (true, true) => {
                    for &i in &self.order[na.first..na.first + na.count] {
                        for &j in &other.order[nb.first..nb.first + nb.count] {
                            if near(&self.boxes[i], &other.boxes[j]) {
                                out.push((i, j));
                            }
                        }
                    }
                }
                (false, true) => stack.extend([(na.first, b), (na.first + 1, b)]),
                (true, false) => stack.extend([(a, nb.first), (a, nb.first + 1)]),
                (false, false) => {
                    let size = |n: &Node| (n.bbox.max - n.bbox.min).length();
                    if size(na) >= size(nb) {
                        stack.extend([(na.first, b), (na.first + 1, b)]);
                    } else {
                        stack.extend([(a, nb.first), (a, nb.first + 1)]);
                    }
                }
            }
        }
        out
    }
}

/// A pair of nodes waiting in [`Bvh::closest_pair`], the nearest first.
struct Pending(f64, usize, usize);

impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for Pending {}

impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Pending {
    // Reversed: the heap pops the smallest distance.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.0.total_cmp(&self.0)
    }
}

/// A node of a [`Bvh`], to inspect the tree.
#[derive(Debug, Clone, Copy)]
pub struct BvhNode<'a> {
    tree: &'a Bvh,
    index: usize,
}

impl<'a> BvhNode<'a> {
    /// The box around the node's boxes.
    pub fn bbox(&self) -> BoundingBox {
        self.tree.nodes[self.index].bbox
    }

    pub fn is_leaf(&self) -> bool {
        self.tree.nodes[self.index].count > 0
    }

    /// The node's two children, unless it is a leaf.
    pub fn children(&self) -> Option<(BvhNode<'a>, BvhNode<'a>)> {
        let n = &self.tree.nodes[self.index];
        (n.count == 0).then(|| (BvhNode { tree: self.tree, index: n.first }, BvhNode { tree: self.tree, index: n.first + 1 }))
    }

    /// The number of boxes below the node.
    pub fn len(&self) -> usize {
        match self.children() {
            Some((a, b)) => a.len() + b.len(),
            None => self.tree.nodes[self.index].count,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The indices of the boxes of a leaf (none for an inner node).
    pub fn boxes(&self) -> &'a [usize] {
        let n = &self.tree.nodes[self.index];
        &self.tree.order[n.first..n.first + n.count]
    }
}

fn overlap(a: &BoundingBox, b: &BoundingBox) -> bool {
    a.min.x <= b.max.x && b.min.x <= a.max.x && a.min.y <= b.max.y && b.min.y <= a.max.y && a.min.z <= b.max.z && b.min.z <= a.max.z
}

/// True if the ray meets the box (slab test).
fn ray_meets(b: &BoundingBox, origin: Point, direction: Vector) -> bool {
    let (mut t0, mut t1) = (0.0f64, f64::INFINITY);
    for (o, d, lo, hi) in [
        (origin.x, direction.x, b.min.x, b.max.x),
        (origin.y, direction.y, b.min.y, b.max.y),
        (origin.z, direction.z, b.min.z, b.max.z),
    ] {
        if d == 0.0 {
            if o < lo || o > hi {
                return false;
            }
            continue;
        }
        let (a, c) = ((lo - o) / d, (hi - o) / d);
        let (near, far) = if a <= c { (a, c) } else { (c, a) };
        t0 = t0.max(near);
        t1 = t1.min(far);
        if t0 > t1 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(x: f64, y: f64, z: f64, size: f64) -> BoundingBox {
        BoundingBox::new(Point::new(x, y, z), Point::new(x + size, y + size, z + size))
    }

    /// A grid of unit cubes, one every 2 units.
    fn grid(n: usize) -> Vec<BoundingBox> {
        let mut boxes = Vec::new();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    boxes.push(cube(2.0 * i as f64, 2.0 * j as f64, 2.0 * k as f64, 1.0));
                }
            }
        }
        boxes
    }

    fn brute(boxes: &[BoundingBox], test: impl Fn(&BoundingBox) -> bool) -> Vec<usize> {
        (0..boxes.len()).filter(|&i| test(&boxes[i])).collect()
    }

    fn sorted(mut v: Vec<usize>) -> Vec<usize> {
        v.sort_unstable();
        v
    }

    #[test]
    fn queries_match_testing_every_box() {
        let boxes = grid(6);
        let tree = Bvh::new(&boxes);
        assert_eq!(tree.len(), boxes.len());
        let query = cube(2.5, 0.5, 3.5, 3.0);
        assert_eq!(sorted(tree.overlapping(&query)), brute(&boxes, |b| overlap(b, &query)));
        let p = Point::new(4.5, 6.5, 8.0);
        assert_eq!(sorted(tree.containing(p)), brute(&boxes, |b| b.contains(p)));
        let (o, d) = (Point::new(-1.0, 0.5, 0.5), Vector::new(1.0, 0.3, 0.2));
        assert_eq!(sorted(tree.along_ray(o, d)), brute(&boxes, |b| ray_meets(b, o, d)));
        // Along an axis, through one row of cubes.
        assert_eq!(tree.along_ray(Point::new(-1.0, 4.5, 6.5), Vector::new(1.0, 0.0, 0.0)).len(), 6);
    }

    #[test]
    fn pairs_match_testing_every_pair() {
        let a = grid(5);
        let b: Vec<BoundingBox> = grid(4).iter().map(|x| BoundingBox::new(x.min + Vector::new(0.7, 0.7, 0.7), x.max + Vector::new(0.7, 0.7, 0.7))).collect();
        let (ta, tb) = (Bvh::new(&a), Bvh::new(&b));
        let mut pairs = ta.overlapping_pairs(&tb);
        pairs.sort_unstable();
        let mut expected = Vec::new();
        for i in 0..a.len() {
            for j in 0..b.len() {
                if overlap(&a[i], &b[j]) {
                    expected.push((i, j));
                }
            }
        }
        assert_eq!(pairs, expected);
    }

    #[test]
    fn empty_and_unbounded_boxes() {
        let boxes = vec![BoundingBox::EMPTY, BoundingBox::WHOLE, cube(0.0, 0.0, 0.0, 1.0)];
        let tree = Bvh::new(&boxes);
        assert_eq!(tree.len(), 2);
        assert_eq!(sorted(tree.containing(Point::new(0.5, 0.5, 0.5))), vec![1, 2]);
        assert_eq!(tree.containing(Point::new(5.0, 5.0, 5.0)), vec![1]);
        assert!(Bvh::new(&[]).overlapping(&BoundingBox::WHOLE).is_empty());
    }
}
