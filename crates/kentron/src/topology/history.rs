//! The record of which entries an operation derived from which.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;

use super::id::NodeId;

/// How the outputs of a [`HistoryEvent`] relate to its inputs.
#[doc(alias = "BRepGraph_LayerHistory::Kind")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChangeKind {
    /// The input lives on as the outputs (same entity, changed).
    Modified,
    /// The outputs are new entities produced from the input.
    Generated,
    /// The input was consumed and has no image.
    Deleted,
    /// The input was detached and continued by the output.
    Replaced,
}

/// One recorded operation.
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryEvent<K = NodeId> {
    /// Increases with every recorded event, also across [`History::clear`].
    pub sequence: u64,
    pub operation: String,
    pub kind: ChangeKind,
    /// Each input with its images (empty for deleted inputs).
    pub mapping: Vec<(K, Vec<K>)>,
}

/// Modification history: an append-only list of [`HistoryEvent`]s with lookup maps in both
/// directions. A model keeps one keyed by [`NodeId`] once [`Model::enable_history`] is
/// called ([`Model::history`]); keyed by
/// [`StableId`](super::StableId) it can describe changes across models or compactions.
///
/// [`Model::history`]: super::Model::history
/// [`Model::enable_history`]: super::Model::enable_history
#[doc(alias = "BRepGraph_LayerHistory")]
#[doc(alias = "BRepTools_History")]
#[derive(Debug, Clone)]
pub struct History<K = NodeId> {
    enabled: bool,
    next_sequence: u64,
    events: Vec<HistoryEvent<K>>,
    modified: HashMap<K, Vec<K>>,
    generated: HashMap<K, Vec<K>>,
    originals: HashMap<K, Vec<K>>,
    deleted: HashSet<K>,
}

impl<K> Default for History<K> {
    fn default() -> Self {
        Self {
            enabled: true,
            next_sequence: 0,
            events: Vec::new(),
            modified: HashMap::new(),
            generated: HashMap::new(),
            originals: HashMap::new(),
            deleted: HashSet::new(),
        }
    }
}

impl<K: Copy + Eq + Hash> History<K> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_enabled(&self) -> bool {
        self.enabled
    }

    /// While disabled, recording calls do nothing (the recorded events are kept). A new
    /// history is enabled.
    pub fn set_enabled(&mut self, enabled: bool) {
        self.enabled = enabled;
    }

    /// The number of recorded events.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// The `i`-th event, oldest first.
    pub fn event(&self, i: usize) -> &HistoryEvent<K> {
        &self.events[i]
    }

    pub fn events(&self) -> &[HistoryEvent<K>] {
        &self.events
    }

    /// Records that `original` became `images`: a [`ChangeKind::Modified`] event, or a
    /// [`ChangeKind::Deleted`] one if `images` is empty.
    pub fn record(&mut self, operation: &str, original: K, images: &[K]) {
        self.record_as(operation, original, images, ChangeKind::Modified);
    }

    /// Like [`record`](Self::record) with an explicit kind (`Modified` or `Generated`);
    /// empty `images` still record a deletion.
    pub fn record_as(&mut self, operation: &str, original: K, images: &[K], kind: ChangeKind) {
        self.push(operation, kind, vec![(original, images.to_vec())]);
    }

    /// One event mapping each of `originals` to the image at the same position.
    pub fn record_pairs(&mut self, operation: &str, originals: &[K], images: &[K], kind: ChangeKind) {
        assert_eq!(originals.len(), images.len(), "record_pairs: length mismatch");
        let mapping = originals.iter().zip(images).map(|(&o, &i)| (o, vec![i])).collect();
        self.push(operation, kind, mapping);
    }

    /// One event marking `nodes` as deleted; nothing is recorded if `nodes` is empty.
    pub fn record_deleted(&mut self, operation: &str, nodes: &[K]) {
        if nodes.is_empty() {
            return;
        }
        let mapping = nodes.iter().map(|&n| (n, Vec::new())).collect();
        self.push(operation, ChangeKind::Deleted, mapping);
    }

    /// Records that `original` was detached and continued by `replacement`: `original`
    /// counts as deleted and `replacement` as its modified image.
    pub fn record_replaced(&mut self, operation: &str, original: K, replacement: K) {
        self.push(operation, ChangeKind::Replaced, vec![(original, vec![replacement])]);
    }

    fn push(&mut self, operation: &str, kind: ChangeKind, mapping: Vec<(K, Vec<K>)>) {
        if !self.enabled || mapping.is_empty() {
            return;
        }
        let kind = match kind {
            ChangeKind::Modified | ChangeKind::Generated
                if mapping.iter().all(|(_, images)| images.is_empty()) =>
            {
                ChangeKind::Deleted
            }
            kind => kind,
        };
        for (original, images) in &mapping {
            if images.is_empty() || matches!(kind, ChangeKind::Deleted | ChangeKind::Replaced) {
                self.deleted.insert(*original);
            }
            if images.is_empty() {
                continue;
            }
            let forward = match kind {
                ChangeKind::Generated => &mut self.generated,
                _ => &mut self.modified,
            };
            push_unique(forward.entry(*original).or_default(), images);
            for image in images {
                push_unique(self.originals.entry(*image).or_default(), &[*original]);
            }
        }
        self.events.push(HistoryEvent {
            sequence: self.next_sequence,
            operation: operation.to_owned(),
            kind,
            mapping,
        });
        self.next_sequence += 1;
    }

    /// The modified images of `node`, if it has any.
    pub fn modified(&self, node: K) -> Option<&[K]> {
        self.modified.get(&node).map(Vec::as_slice)
    }

    /// The entries generated from `node`, if any.
    pub fn generated(&self, node: K) -> Option<&[K]> {
        self.generated.get(&node).map(Vec::as_slice)
    }

    /// The direct inputs `node` was derived from, if any.
    pub fn originals(&self, node: K) -> Option<&[K]> {
        self.originals.get(&node).map(Vec::as_slice)
    }

    /// The root of `node`'s derivation chain (following the first input at each step), or
    /// `node` itself if it was not derived from anything.
    pub fn find_original(&self, node: K) -> K {
        let mut current = node;
        let mut seen = HashSet::from([node]);
        while let Some(&previous) = self.originals.get(&current).and_then(|o| o.first()) {
            if !seen.insert(previous) {
                break;
            }
            current = previous;
        }
        current
    }

    /// Every entry transitively derived from `node` (modified or generated), excluding
    /// `node` itself.
    pub fn find_derived(&self, node: K) -> Vec<K> {
        let mut result = Vec::new();
        let mut seen = HashSet::from([node]);
        let mut stack = vec![node];
        while let Some(current) = stack.pop() {
            let images = self.modified.get(&current).into_iter().chain(self.generated.get(&current));
            for &image in images.flatten() {
                if seen.insert(image) {
                    result.push(image);
                    stack.push(image);
                }
            }
        }
        result
    }

    pub fn is_deleted(&self, node: K) -> bool {
        self.deleted.contains(&node)
    }

    pub fn deleted(&self) -> &HashSet<K> {
        &self.deleted
    }

    /// True if `node` is an input of some event.
    pub fn is_known_input(&self, node: K) -> bool {
        self.events.iter().any(|e| e.mapping.iter().any(|(o, _)| *o == node))
    }

    /// Forgets all events (the enabled flag is kept).
    pub fn clear(&mut self) {
        self.events.clear();
        self.modified.clear();
        self.generated.clear();
        self.originals.clear();
        self.deleted.clear();
    }
}

impl History<NodeId> {
    /// Rewrites every stored node id with `map` (after the model renumbered its nodes).
    pub(crate) fn remap(&mut self, map: &dyn Fn(NodeId) -> NodeId) {
        for event in &mut self.events {
            for (original, images) in &mut event.mapping {
                *original = map(*original);
                for image in images.iter_mut() {
                    *image = map(*image);
                }
            }
        }
        let remap_lists = |m: &mut HashMap<NodeId, Vec<NodeId>>| {
            *m = std::mem::take(m).into_iter().map(|(k, v)| (map(k), v.into_iter().map(map).collect())).collect();
        };
        remap_lists(&mut self.modified);
        remap_lists(&mut self.generated);
        remap_lists(&mut self.originals);
        self.deleted = std::mem::take(&mut self.deleted).into_iter().map(map).collect();
    }
}

fn push_unique<K: Copy + PartialEq>(list: &mut Vec<K>, items: &[K]) {
    for item in items {
        if !list.contains(item) {
            list.push(*item);
        }
    }
}
