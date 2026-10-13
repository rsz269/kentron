//! Recording replacements of sub-shapes and rebuilding shapes with them applied.

use super::id::NodeKind;
use super::shape::Shape;
use super::Model;

/// A set of replacements and removals of shapes, applied to a shape by rebuilding every
/// ancestor of a replaced shape. Shapes are matched by node (and, optionally, placement),
/// regardless of orientation.
///
/// Replacement chains (`a -> b`, `b -> c`) are followed by [`Substitution::resolve`];
/// cycles (`a -> b -> a`, or a replacement that contains the replaced shape) are allowed,
/// and [`Substitution::apply`] stops at them instead of recursing forever.
#[doc(alias = "BRepTools_ReShape")]
#[derive(Debug, Clone, Default)]
pub struct Substitution {
    replacements: Vec<(Shape, Option<Shape>)>,
    match_placement: bool,
}

impl Substitution {
    pub fn new() -> Self {
        Self::default()
    }

    /// Also distinguish shapes by placement when matching (by default only the node
    /// counts).
    pub fn matching_placement(self) -> Self {
        Self { match_placement: true, ..self }
    }

    /// The index of the record for `shape`, if any.
    fn find(&self, shape: Shape) -> Option<usize> {
        self.replacements
            .iter()
            .position(|(old, _)| old.node == shape.node && (!self.match_placement || old.placement == shape.placement))
    }

    /// Records that `old` is to be replaced by `new`, overwriting any earlier record for
    /// `old`. Replacing a shape by itself records nothing.
    pub fn replace(&mut self, old: Shape, new: Shape) {
        if old.node == new.node && (!self.match_placement || old.placement == new.placement) {
            return;
        }
        match self.find(old) {
            Some(i) => self.replacements[i] = (old, Some(new)),
            None => self.replacements.push((old, Some(new))),
        }
    }

    /// Records that `old` is to be removed.
    pub fn remove(&mut self, old: Shape) {
        match self.find(old) {
            Some(i) => self.replacements[i] = (old, None),
            None => self.replacements.push((old, None)),
        }
    }

    pub fn is_recorded(&self, shape: Shape) -> bool {
        self.find(shape).is_some()
    }

    /// The direct replacement of `shape`: itself if unrecorded, `None` if removed.
    pub fn get(&self, shape: Shape) -> Option<Shape> {
        let Some(i) = self.find(shape) else { return Some(shape) };
        let (old, new) = self.replacements[i];
        // Oriented and placed as `shape` is relative to the recorded shape.
        new.map(|n| {
            let orientation = if shape.orientation == old.orientation { n.orientation } else { n.orientation.reversed() };
            let placement = if self.match_placement { n.placement } else { shape.placement * old.placement.inverse() * n.placement };
            Shape { node: n.node, orientation, placement }
        })
    }

    /// The end of `shape`'s replacement chain: itself if unrecorded, `None` if the chain
    /// ends in a removal.
    pub fn resolve(&self, shape: Shape) -> Option<Shape> {
        let mut current = shape;
        let mut seen = vec![current.node];
        while self.is_recorded(current) {
            let next = self.get(current)?;
            // A cycle ends where it would come round.
            if seen.contains(&next.node) {
                return Some(current);
            }
            seen.push(next.node);
            current = next;
        }
        Some(current)
    }

    pub fn clear(&mut self) {
        self.replacements.clear();
    }

    /// `shape` with all recorded replacements applied below it, rebuilding changed
    /// ancestors as new nodes in `model`; `shape` itself if nothing below it changed,
    /// `None` if it was removed.
    pub fn apply(&self, model: &mut Model, shape: Shape) -> Option<Shape> {
        self.apply_until(model, shape, NodeKind::Vertex)
    }

    /// Like [`Substitution::apply`], but only looks for recorded shapes down to the
    /// `until` level (e.g. `NodeKind::Compound` to replace sub-compounds only).
    pub fn apply_until(&self, model: &mut Model, shape: Shape, until: NodeKind) -> Option<Shape> {
        let mut done = std::collections::HashMap::new();
        self.rebuild(model, shape, until, &mut done)
    }

    /// `shape` rebuilt; `done` maps the nodes already rebuilt (as used forward and
    /// unplaced) to their results.
    fn rebuild(
        &self,
        model: &mut Model,
        shape: Shape,
        until: NodeKind,
        done: &mut std::collections::HashMap<super::NodeId, Option<Shape>>,
    ) -> Option<Shape> {
        if shape.kind() > until {
            return Some(shape);
        }
        if self.is_recorded(shape) {
            // Replaced as recorded (not rebuilt again: a replacement may contain the shape).
            return self.resolve(shape);
        }
        // The node rebuilt once, then used as `shape` uses it.
        let local = |s: Shape| Shape { orientation: super::Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..s };
        let rebuilt = match done.get(&shape.node) {
            Some(r) => *r,
            None => {
                let r = self.rebuild_node(model, local(shape), until, done);
                done.insert(shape.node, r);
                r
            }
        };
        rebuilt.map(|r| Shape {
            node: r.node,
            orientation: r.orientation.compose(shape.orientation),
            placement: shape.placement * r.placement,
        })
    }

    fn rebuild_node(
        &self,
        model: &mut Model,
        shape: Shape,
        until: NodeKind,
        done: &mut std::collections::HashMap<super::NodeId, Option<Shape>>,
    ) -> Option<Shape> {
        let children: Vec<Shape> = model.children_local(shape).collect();
        let mut rebuilt = Vec::with_capacity(children.len());
        let mut changed = false;
        for child in &children {
            let r = self.rebuild(model, *child, until, done);
            changed |= r != Some(*child);
            rebuilt.push(r);
        }
        if !changed {
            return Some(shape);
        }
        let copy = model.empty_copy(shape);
        if shape.kind() == NodeKind::CoEdge {
            // A coedge holds its edge: a new coedge on the replacement edge.
            let edge = rebuilt.into_iter().flatten().next()?;
            let id = super::EdgeId::try_from(edge.node).ok()?;
            let coedge = model.add_coedge(id, edge.orientation).ok()?;
            let old = super::CoEdgeId::try_from(shape.node).ok()?;
            if let Some((pcurve, first, last)) = model.pcurve(old).cloned().map(|p| (p, model.coedge_range(old).0, model.coedge_range(old).1)) {
                let _ = model.set_pcurve_with_range(coedge, pcurve, first, last);
            }
            model.remove_node(copy.node);
            return Some(Shape::from(coedge));
        }
        for child in rebuilt.into_iter().flatten() {
            if model.add(copy, child).is_err() {
                return None;
            }
        }
        Some(copy)
    }
}
