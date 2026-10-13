//! Typed ids into [`Model`](super::Model) tables.

use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

/// Index of a `T` in its model table. Cheap (`u32`), `Copy`, and only meaningful for the
/// model that issued it. Not stable across [`Model::compact`](super::Model::compact); use
/// [`StableId`] for identity that must survive compaction.
pub struct Id<T> {
    index: u32,
    _kind: PhantomData<fn() -> T>,
}

impl<T> Id<T> {
    pub const INVALID: Self = Self::from_index(u32::MAX);

    pub const fn from_index(index: u32) -> Self {
        Self { index, _kind: PhantomData }
    }

    pub const fn index(self) -> u32 {
        self.index
    }

    pub const fn is_valid(self) -> bool {
        self.index != u32::MAX
    }
}

// Manual impls: derives would needlessly require `T: Clone` etc.
impl<T> Clone for Id<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Id<T> {}
impl<T> PartialEq for Id<T> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}
impl<T> Eq for Id<T> {}
impl<T> PartialOrd for Id<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<T> Ord for Id<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index)
    }
}
impl<T> Hash for Id<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state)
    }
}
impl<T> Default for Id<T> {
    fn default() -> Self {
        Self::INVALID
    }
}
impl<T> fmt::Debug for Id<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = std::any::type_name::<T>().rsplit("::").next().unwrap_or("?");
        match self.is_valid() {
            true => write!(f, "{name}#{}", self.index),
            false => write!(f, "{name}#invalid"),
        }
    }
}

/// The kinds of topology node, from the most complex to the simplest.
#[doc(alias = "TopAbs_ShapeEnum")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NodeKind {
    /// A placed instance of a product.
    Instance,
    /// A reusable part or assembly definition.
    Product,
    Compound,
    CompSolid,
    Solid,
    Shell,
    Face,
    Wire,
    CoEdge,
    Edge,
    Vertex,
}

impl NodeKind {
    /// Every kind, in declaration order.
    pub const ALL: [NodeKind; 11] = [
        NodeKind::Instance,
        NodeKind::Product,
        NodeKind::Compound,
        NodeKind::CompSolid,
        NodeKind::Solid,
        NodeKind::Shell,
        NodeKind::Face,
        NodeKind::Wire,
        NodeKind::CoEdge,
        NodeKind::Edge,
        NodeKind::Vertex,
    ];
}

/// An id of a topology node of any kind.
#[doc(alias = "BRepGraph_NodeId")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId {
    pub kind: NodeKind,
    pub index: u32,
}

/// Identity of a model entry that survives [`Model::compact`](super::Model::compact):
/// a counter issued once per model and never reused. Use it for history and for
/// references kept outside the model.
#[doc(alias = "BRepGraph_UID")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct StableId(pub u32);

/// A topology entity type: one that a [`NodeId`] can refer to.
pub trait Node: super::Record {
    const KIND: NodeKind;
}

macro_rules! nodes {
    ($($ty:ident => $kind:ident),* $(,)?) => {
        $(
            impl Node for super::$ty {
                const KIND: NodeKind = NodeKind::$kind;
            }
        )*
    };
}

nodes! {
    Instance => Instance, Product => Product, Compound => Compound, CompSolid => CompSolid,
    Solid => Solid, Shell => Shell, Face => Face, Wire => Wire, CoEdge => CoEdge,
    Edge => Edge, Vertex => Vertex,
}

impl<T: Node> From<Id<T>> for NodeId {
    fn from(id: Id<T>) -> Self {
        NodeId { kind: T::KIND, index: id.index() }
    }
}

/// Fails (returning the node back) if the node is of another kind.
impl<T: Node> TryFrom<NodeId> for Id<T> {
    type Error = NodeId;

    fn try_from(node: NodeId) -> Result<Self, NodeId> {
        match node.kind == T::KIND {
            true => Ok(Id::from_index(node.index)),
            false => Err(node),
        }
    }
}

impl<T: Node> PartialEq<NodeId> for Id<T> {
    fn eq(&self, node: &NodeId) -> bool {
        node.kind == T::KIND && node.index == self.index()
    }
}

impl<T: Node> PartialEq<Id<T>> for NodeId {
    fn eq(&self, id: &Id<T>) -> bool {
        id == self
    }
}
