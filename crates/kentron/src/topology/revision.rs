//! Change tracking: per-entry revision counters and version stamps.

use std::sync::atomic::{AtomicU64, Ordering};

use super::id::{NodeId, NodeKind, StableId};
use super::Model;

/// Change counters of one model entry.
///
/// `own` counts edits of the entry itself (its geometry, tolerance, children list, ...);
/// `subtree` also counts edits anywhere below it, so a cache computed from a node's whole
/// subtree is fresh while `subtree` is unchanged.
#[doc(alias = "OwnGen")]
#[doc(alias = "SubtreeGen")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Revision {
    pub own: u32,
    pub subtree: u32,
}

/// Identity of one [`Model`] instance, unique within the process. A clone of a model keeps
/// the identity (its ids mean the same things).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ModelUid(u64);

impl Default for ModelUid {
    fn default() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        ModelUid(NEXT.fetch_add(1, Ordering::Relaxed))
    }
}

/// A node's identity and version at one point in time: compare against the model later to
/// tell whether the node changed ([`Model::is_stale`]).
#[doc(alias = "BRepGraph_VersionStamp")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Stamp {
    pub model: ModelUid,
    pub kind: NodeKind,
    pub stable_id: StableId,
    /// The node's own revision when stamped.
    pub revision: u32,
}

impl Stamp {
    /// True if both stamps name the same node, whatever its version.
    pub fn is_same_item(&self, other: &Stamp) -> bool {
        self.model == other.model && self.stable_id == other.stable_id
    }
}

impl Model {
    /// This model's identity.
    pub fn uid(&self) -> ModelUid {
        self.uid
    }

    /// The change counters of `node` (also for removed nodes); `None` if it is out of range.
    /// The subtree counter adds the edits of everything below the node to its own.
    pub fn revision(&self, node: impl Into<NodeId>) -> Option<Revision> {
        let node = node.into();
        let (_, own, _) = self.node_slot(node)?;
        let mut subtree = own.own.wrapping_add(own.subtree);
        let mut seen = std::collections::HashSet::new();
        let mut stack = self.child_nodes(node);
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            if let Some((_, r, false)) = self.node_slot(n) {
                subtree = subtree.wrapping_add(r.own);
                stack.extend(self.child_nodes(n));
            }
        }
        Some(Revision { own: own.own, subtree })
    }

    fn reset_revision(&mut self, node: NodeId) {
        self.bump(node, false);
    }

    /// Resets the change counters of `shape` and everything below it: what a builder made
    /// starts unedited. (What values derived from the nodes go by is unchanged.)
    pub(crate) fn settle(&mut self, shape: super::Shape) {
        let mut seen = std::collections::HashSet::new();
        let mut stack = vec![shape.node];
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            self.reset_revision(n);
            stack.extend(self.child_nodes(n));
        }
    }

    /// A stamp of `node`'s current version; `None` if it is invalid or removed.
    pub fn stamp(&self, node: impl Into<NodeId>) -> Option<Stamp> {
        let node = node.into();
        let (stable_id, revision, removed) = self.node_slot(node)?;
        (!removed).then_some(Stamp {
            model: self.uid,
            kind: node.kind,
            stable_id,
            revision: revision.own,
        })
    }

    /// True if the stamped node was edited or removed since, or `stamp` comes from another
    /// model.
    pub fn is_stale(&self, stamp: &Stamp) -> bool {
        if stamp.model != self.uid {
            return true;
        }
        let Some(node) = self.find_stable_node(stamp.stable_id) else { return true };
        node.kind != stamp.kind || self.node_slot(node).is_none_or(|(_, r, _)| r.own != stamp.revision)
    }

    /// Marks `node` as edited: bumps its own revision and the subtree revision of it and
    /// all its ancestors.
    #[doc(alias = "MarkDirty")]
    pub fn touch(&mut self, node: impl Into<NodeId>) {
        let node = node.into();
        if self.is_live(node) {
            self.bump(node, true);
        }
    }

    /// Counts an edit of `node` (if `own`) or, without, sets its counters back to 0.
    fn bump(&mut self, node: NodeId, own: bool) {
        macro_rules! bump {
            ($table:ident) => {{
                let id = super::Id::from_index(node.index);
                if own {
                    let _ = self.$table.get_mut(id);
                } else if let Some(r) = self.$table.revision_mut(id) {
                    *r = Revision::default();
                }
            }};
        }
        match node.kind {
            NodeKind::Vertex => bump!(vertices),
            NodeKind::Edge => bump!(edges),
            NodeKind::CoEdge => bump!(coedges),
            NodeKind::Wire => bump!(wires),
            NodeKind::Face => bump!(faces),
            NodeKind::Shell => bump!(shells),
            NodeKind::Solid => bump!(solids),
            NodeKind::Compound => bump!(compounds),
            NodeKind::CompSolid => bump!(comp_solids),
            NodeKind::Product => bump!(products),
            NodeKind::Instance => bump!(instances),
        }
    }
}
