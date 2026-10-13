//! Walking down the topology graph.

use super::id::{NodeId, NodeKind};
use crate::math::Transform;
use super::shape::Shape;
use super::Model;

/// Iterator over the nodes below a root shape, visited depth-first (children in order). A
/// child used by several parents is visited once per use, with the orientation and
/// placement accumulated along that path. The walk goes through every link: solids to
/// shells, faces to wires, wires to coedges, coedges to edges, edges to vertices,
/// compounds to children, assemblies to instances to products to part roots.
///
/// Created by [`Model::explore`] (nodes of one kind) or [`Model::explore_all`]; refine
/// it with the builder methods before iterating.
#[doc(alias = "TopExp_Explorer")]
#[doc(alias = "BRepGraph_ChildExplorer")]
#[derive(Debug, Clone)]
pub struct Explorer<'a> {
    model: &'a Model,
    root: Shape,
    kind: Option<NodeKind>,
    skip_below: Option<NodeKind>,
    stop_at: Option<NodeKind>,
    compose_orientation: bool,
    compose_placement: bool,
    /// Shapes still to visit.
    stack: Vec<Step>,
    /// The shapes from the root to the last visited one.
    path: Vec<Step>,
    /// The placements met (relative to the root's frame): the root's first. Steps name
    /// them by index, so that the stack and path hold small entries and only children
    /// placed in their parents (of compounds and assemblies) need a product.
    placements: Vec<Transform>,
}

/// A shape on the way: its node, composed orientation, placement (an index into
/// [`Explorer::placements`]) and depth below the root.
#[derive(Debug, Clone, Copy)]
struct Step {
    node: NodeId,
    orientation: super::Orientation,
    placement: u32,
    depth: u32,
}

impl Explorer<'_> {
    /// Does not enter nodes of `kind` (nor yield them). Ignored if `kind` is the target.
    pub fn skip_below(self, kind: NodeKind) -> Self {
        Self { skip_below: Some(kind), ..self }
    }

    /// Yields nodes of `kind` instead of descending into them. Ignored if `kind` is the
    /// target.
    pub fn stop_at(self, kind: NodeKind) -> Self {
        Self { stop_at: Some(kind), ..self }
    }

    /// When false, every yielded shape is forward.
    pub fn compose_orientation(self, compose: bool) -> Self {
        Self { compose_orientation: compose, ..self }
    }

    /// When false, every yielded shape has an identity placement.
    pub fn compose_placement(self, compose: bool) -> Self {
        Self { compose_placement: compose, ..self }
    }

    /// The node of `kind` on the path from the root to the last yielded shape (the root
    /// counts), if any.
    pub fn ancestor(&self, kind: NodeKind) -> Option<Shape> {
        let ancestors = &self.path[..self.path.len().saturating_sub(1)];
        ancestors.iter().rev().find(|s| s.node.kind == kind).map(|s| self.shape(s))
    }

    fn shape(&self, step: &Step) -> Shape {
        Shape { node: step.node, orientation: step.orientation, placement: self.placements[step.placement as usize] }
    }
}

impl Iterator for Explorer<'_> {
    type Item = Shape;

    fn next(&mut self) -> Option<Shape> {
        while let Some(step) = self.stack.pop() {
            let depth = step.depth as usize;
            self.path.truncate(depth);
            self.path.push(step);
            let kind = step.node.kind;
            let is_target = self.kind.is_none_or(|k| k == kind);
            let excluded = |option: Option<NodeKind>| option == Some(kind) && self.kind != Some(kind);
            if excluded(self.skip_below) && depth > 0 {
                continue;
            }
            let stop = excluded(self.stop_at);
            // A target of a given kind is not entered: nothing of the same kind lies below,
            // except for instances, which assemblies nest.
            let nests = matches!(kind, NodeKind::Instance | NodeKind::Product);
            let descend = !stop && (self.kind.is_none() || !is_target || nests);
            if descend {
                // The children onto the stack, then turned round in place (the first on top).
                let base = self.stack.len();
                // Only compounds and assemblies place their children.
                let placing = matches!(kind, NodeKind::Compound | NodeKind::Product);
                let mut i = 0;
                while let Some(child) = self.model.child_at(step.node, i) {
                    i += 1;
                    let placement = if placing && child.placement != Transform::IDENTITY {
                        let parent = self.placements[step.placement as usize];
                        self.placements.push(if parent == Transform::IDENTITY { child.placement } else { parent * child.placement });
                        (self.placements.len() - 1) as u32
                    } else {
                        step.placement
                    };
                    self.stack.push(Step { node: child.node, orientation: child.orientation.compose(step.orientation), placement, depth: step.depth + 1 });
                }
                self.stack[base..].reverse();
            }
            // `explore_all` lists what is below the root, not the root itself.
            let is_root_of_all = self.kind.is_none() && depth == 0;
            if (is_target || stop) && !is_root_of_all {
                let mut out = self.shape(&step);
                if !self.compose_orientation {
                    out.orientation = super::Orientation::Forward;
                }
                if !self.compose_placement {
                    out.placement = Transform::IDENTITY;
                }
                return Some(out);
            }
        }
        None
    }
}

impl Model {
    /// Every occurrence of a `kind` node below `shape` (see [`Explorer`]); `shape` itself if
    /// it is of that kind. A box has 6 faces but 24 edge occurrences: each of its 12 edges
    /// is reached through 2 faces.
    pub fn explore(&self, shape: impl Into<Shape>, kind: NodeKind) -> Explorer<'_> {
        self.explorer(shape.into(), Some(kind))
    }

    /// Every node occurrence below `shape` (not `shape` itself), of any kind.
    pub fn explore_all(&self, shape: impl Into<Shape>) -> Explorer<'_> {
        self.explorer(shape.into(), None)
    }

    fn explorer(&self, root: Shape, kind: Option<NodeKind>) -> Explorer<'_> {
        Explorer {
            model: self,
            root,
            kind,
            skip_below: None,
            stop_at: None,
            compose_orientation: true,
            compose_placement: true,
            stack: vec![Step { node: root.node, orientation: root.orientation, placement: 0, depth: 0 }],
            path: Vec::new(),
            placements: vec![root.placement],
        }
    }

    /// The distinct `kind` sub-shapes below `shape`, ignoring orientation, in order of
    /// first occurrence.
    #[doc(alias = "TopExp::MapShapes")]
    pub fn unique_subshapes(&self, shape: impl Into<Shape>, kind: NodeKind) -> Vec<Shape> {
        let mut unique: Vec<Shape> = Vec::new();
        let mut seen = crate::hashing::FastSet::default();
        for found in self.explore(shape, kind) {
            // Same node at the same placement, whatever the orientation.
            let placed_elsewhere = found.placement != Transform::IDENTITY;
            if placed_elsewhere {
                if !unique.iter().any(|u| u.node == found.node && u.placement == found.placement) {
                    unique.push(found);
                }
            } else if seen.insert(found.node) {
                unique.push(found);
            }
        }
        unique
    }
}

/// Iterator over the ancestors of a node, walking up every path to the top separately: a
/// node reachable along several paths (a solid used twice in an assembly) yields its
/// ancestors once per path. Each ancestor is yielded nearest first, as a [`Shape`] whose
/// placement is the ancestor's placement relative to the top of the current path. The
/// walk goes up through every link, including part roots, instances and assemblies.
///
/// Created by [`Model::explore_up`].
#[doc(alias = "BRepGraph_ParentExplorer")]
#[derive(Debug, Clone)]
pub struct AncestorExplorer<'a> {
    model: &'a Model,
    leaf: NodeId,
    kind: Option<NodeKind>,
    direct_only: bool,
    skip_above: Option<NodeKind>,
    stop_at: Option<NodeKind>,
    /// The paths up, found on the first step: each from the leaf's parent to the top, with
    /// each ancestor's placement and orientation relative to that top, and the leaf's.
    paths: Option<Vec<AncestorPath>>,
    /// The path and position in it of the next ancestor.
    at: (usize, usize),
    /// The path of the last yielded ancestor.
    current: Option<usize>,
}

#[derive(Debug, Clone)]
struct AncestorPath {
    /// Ancestors nearest first, each with whether it is yielded.
    nodes: Vec<(Shape, bool)>,
    leaf: Shape,
}

impl AncestorExplorer<'_> {
    /// Only yields ancestors of `kind`.
    pub fn of_kind(self, kind: NodeKind) -> Self {
        Self { kind: Some(kind), ..self }
    }

    /// Only yields the direct parents.
    pub fn direct_only(self) -> Self {
        Self { direct_only: true, ..self }
    }

    /// Does not go above nodes of `kind` (nor yield them). Ignored if `kind` is the
    /// target kind.
    pub fn skip_above(self, kind: NodeKind) -> Self {
        Self { skip_above: Some(kind), ..self }
    }

    /// Yields nodes of `kind` instead of going above them. Ignored if `kind` is the target
    /// kind.
    pub fn stop_at(self, kind: NodeKind) -> Self {
        Self { stop_at: Some(kind), ..self }
    }

    /// The placement of the starting node relative to the top of the current path.
    pub fn leaf_placement(&self) -> Transform {
        self.current_leaf().map_or(Transform::IDENTITY, |l| l.placement)
    }

    /// The orientation of the starting node relative to the top of the current path.
    pub fn leaf_orientation(&self) -> super::Orientation {
        self.current_leaf().map_or(super::Orientation::Forward, |l| l.orientation)
    }

    fn current_leaf(&self) -> Option<Shape> {
        Some(self.paths.as_ref()?.get(self.current?)?.leaf)
    }

    /// Every path up from the leaf, cut where `skip_above` and `stop_at` say.
    fn find_paths(&self) -> Vec<AncestorPath> {
        use super::Orientation;
        let target = self.kind;
        let excluded = |option: Option<NodeKind>, kind: NodeKind| option == Some(kind) && target != Some(kind);
        // Each partial path: the links from the leaf up (parent, orientation, placement).
        let mut complete: Vec<Vec<(NodeId, Orientation, Transform, bool)>> = Vec::new();
        let mut stack: Vec<Vec<(NodeId, Orientation, Transform, bool)>> = vec![Vec::new()];
        while let Some(path) = stack.pop() {
            let last = path.last().map_or(self.leaf, |l| l.0);
            let stopped = path.last().is_some_and(|l| l.3);
            let links = if stopped || (self.direct_only && !path.is_empty()) { Vec::new() } else { self.model.parent_links(last) };
            let mut grown = false;
            // In reverse, so that the first parent's paths come first.
            for link in links.into_iter().rev() {
                let kind = link.parent.kind;
                if excluded(self.skip_above, kind) || path.iter().any(|l| l.0 == link.parent) || link.parent == self.leaf {
                    continue;
                }
                let stop = excluded(self.stop_at, kind);
                let mut next = path.clone();
                next.push((link.parent, link.orientation, link.placement, stop));
                stack.push(next);
                grown = true;
            }
            if !grown && !path.is_empty() {
                complete.push(path);
            }
        }
        complete
            .into_iter()
            .map(|links| {
                // Placements relative to the top: composed from the top down.
                let n = links.len();
                let mut placement = Transform::IDENTITY;
                let mut orientation = Orientation::Forward;
                let mut nodes = vec![(Shape::new(links[0].0), false); n];
                for k in (0..n).rev() {
                    let (node, _, _, stop) = links[k];
                    let yielded = target.is_none_or(|t| t == node.kind) || stop;
                    nodes[k] = (Shape { node, orientation, placement }, yielded);
                    // The link of the node below within this one.
                    let (_, o, p, _) = links[k];
                    placement = placement * p;
                    orientation = o.compose(orientation);
                }
                AncestorPath { nodes, leaf: Shape { node: self.leaf, orientation, placement } }
            })
            .collect()
    }
}

impl Iterator for AncestorExplorer<'_> {
    type Item = Shape;

    fn next(&mut self) -> Option<Shape> {
        if self.paths.is_none() {
            self.paths = Some(self.find_paths());
        }
        let paths = self.paths.as_ref().expect("the paths");
        loop {
            let (p, i) = self.at;
            let path = paths.get(p)?;
            if i >= path.nodes.len() {
                self.at = (p + 1, 0);
                continue;
            }
            self.at = (p, i + 1);
            let (shape, yielded) = path.nodes[i];
            if yielded {
                self.current = Some(p);
                return Some(shape);
            }
        }
    }
}

impl Model {
    /// The ancestors of `node`, path by path (see [`AncestorExplorer`]).
    pub fn explore_up(&self, node: impl Into<NodeId>) -> AncestorExplorer<'_> {
        AncestorExplorer {
            model: self,
            leaf: node.into(),
            kind: None,
            direct_only: false,
            skip_above: None,
            stop_at: None,
            paths: None,
            at: (0, 0),
            current: None,
        }
    }
}
