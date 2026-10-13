//! A placed, oriented reference to a topology node.

use crate::math::Transform;

use super::entities::Orientation;
use super::id::{Id, Node, NodeId, NodeKind};

/// A use of a topology node: which node, how oriented, where placed. Cheap to copy;
/// every query about it takes the [`Model`](super::Model) that owns the node.
#[doc(alias = "TopoDS_Shape")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shape {
    pub node: NodeId,
    pub orientation: Orientation,
    pub placement: Transform,
}

impl Shape {
    /// `node`, forward, unplaced.
    pub fn new(node: NodeId) -> Self {
        Self { node, orientation: Orientation::Forward, placement: Transform::IDENTITY }
    }

    pub fn kind(&self) -> NodeKind {
        self.node.kind
    }
}

impl Shape {
    /// The same use with the opposite orientation.
    #[doc(alias = "TopoDS_Shape::Reversed")]
    pub fn reversed(&self) -> Shape {
        Shape { orientation: self.orientation.reversed(), ..*self }
    }

    /// The same node and orientation, with `placement` replacing the current one.
    #[doc(alias = "TopoDS_Shape::Located")]
    pub fn with_placement(&self, placement: Transform) -> Shape {
        Shape { placement, ..*self }
    }

    /// The same use moved by `t` (applied after the current placement).
    #[doc(alias = "TopoDS_Shape::Moved")]
    pub fn moved(&self, t: &Transform) -> Shape {
        Shape { placement: *t * self.placement, ..*self }
    }

    /// True if both refer to the same node at the same placement, whatever their
    /// orientations.
    #[doc(alias = "TopoDS_Shape::IsSame")]
    pub fn is_same(&self, other: &Shape) -> bool {
        self.node == other.node && self.placement == other.placement
    }
}

/// An unplaced, forward use of the node.
impl From<NodeId> for Shape {
    fn from(node: NodeId) -> Self {
        Shape::new(node)
    }
}

/// An unplaced, forward use of the node.
impl<T: Node> From<Id<T>> for Shape {
    fn from(id: Id<T>) -> Self {
        Shape::new(id.into())
    }
}

/// Values kept per occurrence of a node (its node and placement, whatever its orientation),
/// for algorithms meeting shared sub-shapes several times.
#[derive(Debug, Clone)]
pub(crate) struct ByOccurrence<V>(std::collections::HashMap<NodeId, Vec<(Transform, V)>>);

impl<V> Default for ByOccurrence<V> {
    fn default() -> Self {
        ByOccurrence(std::collections::HashMap::new())
    }
}

impl<V> ByOccurrence<V> {
    pub fn get(&self, s: &Shape) -> Option<&V> {
        self.0.get(&s.node)?.iter().find(|(t, _)| *t == s.placement).map(|(_, v)| v)
    }

    /// Sets the value of `s`'s occurrence, replacing any.
    pub fn insert(&mut self, s: &Shape, v: V) {
        let entry = self.0.entry(s.node).or_default();
        match entry.iter_mut().find(|(t, _)| *t == s.placement) {
            Some(slot) => slot.1 = v,
            None => entry.push((s.placement, v)),
        }
    }
}
