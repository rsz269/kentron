//! Values computed from parts of a model (sampled boundaries, boxes, ...), kept while the
//! parts they were computed from are unchanged.
//!
//! Each value is stored with a fingerprint of the revisions of the records it depends on
//! ([`Model::face_fingerprint`], [`Model::edge_fingerprint`]): every edit of a record moves
//! its revision on, so a changed part gives another fingerprint and the value is computed
//! again.

use std::any::{Any, TypeId};
use std::sync::{Arc, Mutex, MutexGuard};

use super::entities::{EdgeId, FaceId};
use super::id::{Id, NodeId};
use super::Model;

type Entries = crate::hashing::FastMap<(NodeId, TypeId), (u64, Arc<dyn Any + Send + Sync>)>;

/// The cache of derived values of a [`Model`], shared by the threads reading the model
/// (locked only to look values up and to keep them, not while computing them).
#[derive(Default)]
pub struct Derived {
    entries: Mutex<Entries>,
}

impl Derived {
    fn lock(&self) -> MutexGuard<'_, Entries> {
        // (A panic elsewhere while it was held leaves the map whole.)
        self.entries.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
    }
}

impl Derived {
    /// Forgets every kept value.
    pub(crate) fn clear(&self) {
        self.lock().clear();
    }

    /// Forgets the values derived from nodes for which `keep` is false.
    pub(crate) fn retain_nodes(&self, keep: impl Fn(NodeId) -> bool) {
        let dropped: Vec<_> = {
            let mut entries = self.lock();
            let dropped = entries.extract_if(|key, _| !keep(key.0)).collect();
            entries.shrink_to_fit();
            dropped
        };
        // (Freed with the lock released.)
        drop(dropped);
    }

    pub(crate) fn len(&self) -> usize {
        self.lock().len()
    }
}

impl Clone for Derived {
    fn clone(&self) -> Self {
        Derived { entries: Mutex::new(self.lock().clone()) }
    }
}

impl std::fmt::Debug for Derived {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Derived({} values)", self.lock().len())
    }
}

/// Values kept at most; beyond, the cache starts afresh.
const CAPACITY: usize = 1 << 18;

/// A step of the fingerprint hash.
fn mix(hash: u64, value: u64) -> u64 {
    (hash.rotate_left(5) ^ value).wrapping_mul(0x517c_c1b7_2722_0a95)
}

impl Model {
    /// How many derived values are kept.
    #[cfg(test)]
    pub(crate) fn derived_count(&self) -> usize {
        self.derived.len()
    }

    /// The value of type `T` derived from `node`: computed by `compute` unless a value
    /// computed under the same `fingerprint` is kept.
    pub(crate) fn derived<T: Any + Send + Sync>(&self, node: NodeId, fingerprint: u64, compute: impl FnOnce() -> T) -> Arc<T> {
        let key = (node, TypeId::of::<T>());
        if let Some((kept, value)) = self.derived.lock().get(&key) {
            if *kept == fingerprint {
                return value.clone().downcast::<T>().expect("a derived value of its own type");
            }
        }
        let value = Arc::new(compute());
        let old = {
            let mut entries = self.derived.lock();
            let old = if entries.len() >= CAPACITY { std::mem::take(&mut *entries) } else { Default::default() };
            entries.insert(key, (fingerprint, value.clone()));
            old
        };
        // (Dropped with the lock released, and a full cache's many values on a thread of
        // their own: freeing them all takes long, and nothing waits for it.)
        if !old.is_empty() {
            std::thread::spawn(move || drop(old));
        }
        value
    }

    /// A slot's index and how often what it is was edited (not its meshes), mixed into
    /// `hash`.
    fn mix_slot<T: super::Record>(&self, hash: u64, id: Id<T>) -> u64 {
        let content = T::table(self).content(id).map_or(u64::MAX, u64::from);
        mix(mix(hash, id.index() as u64), content)
    }

    /// A slot's index and own revision (every edit), mixed into `hash`.
    fn mix_slot_revision<T: super::Record>(&self, hash: u64, id: Id<T>) -> u64 {
        let revision = T::table(self).revision(id).map_or(u64::MAX, |r| r.own as u64);
        mix(mix(hash, id.index() as u64), revision)
    }

    /// A fingerprint of the meshes and polylines kept with a shape's faces, coedges and
    /// edges: another one after any is set, changed or cleared (which
    /// [`Self::shape_fingerprint`] does not see).
    pub(crate) fn mesh_fingerprint(&self, shape: super::Shape) -> u64 {
        use super::NodeKind;
        let mut hash = 0;
        for face in self.explore(shape, NodeKind::Face) {
            let Ok(id) = FaceId::try_from(face.node) else { continue };
            let Some(record) = self.faces.get(id) else { continue };
            hash = self.mix_slot_revision(hash, record.mesh);
            for &wire_use in &record.wires {
                let Some(wire) = self.wire_uses.get(wire_use).and_then(|u| self.wires.get(u.wire)) else { continue };
                for &coedge in &wire.coedges {
                    if let Some(c) = self.coedges.get(coedge) {
                        hash = self.mix_slot_revision(hash, c.mesh_polyline);
                    }
                }
            }
        }
        for edge in self.explore(shape, NodeKind::Edge) {
            let Ok(id) = EdgeId::try_from(edge.node) else { continue };
            if let Some(record) = self.edges.get(id) {
                hash = self.mix_slot_revision(hash, record.polyline);
            }
        }
        hash
    }

    /// A fingerprint of an edge and what it is made of (its curve, its vertices): another
    /// one after any of them is edited.
    pub(crate) fn edge_fingerprint(&self, edge: EdgeId) -> u64 {
        let mut hash = self.mix_slot(0, edge);
        let Some(record) = self.edges.get(edge) else { return hash };
        hash = self.mix_slot(hash, record.curve);
        for vertex_use in [record.start, record.end] {
            hash = self.mix_slot(hash, vertex_use);
            if let Some(u) = self.vertex_uses.get(vertex_use) {
                hash = self.mix_slot(hash, u.vertex);
            }
        }
        hash
    }

    /// A fingerprint of a shape and everything in it (its faces, edges and vertices, as
    /// placed and oriented): another one after any of them is edited, or once the shape is
    /// made of others.
    pub(crate) fn shape_fingerprint(&self, shape: super::Shape) -> u64 {
        use super::{NodeKind, Shape};
        let placed = |hash: u64, s: &Shape| {
            let rows = s.placement.rows();
            let hash = rows.iter().flatten().fold(hash, |h, v| mix(h, v.to_bits()));
            mix(mix(hash, s.orientation as u64), s.node.index as u64)
        };
        let mut hash = placed(0, &shape);
        for face in self.explore(shape, NodeKind::Face) {
            hash = placed(hash, &face);
            if let Ok(id) = FaceId::try_from(face.node) {
                hash = mix(hash, self.face_fingerprint(id));
            }
        }
        for edge in self.explore(shape, NodeKind::Edge) {
            hash = placed(hash, &edge);
            if let Ok(id) = EdgeId::try_from(edge.node) {
                hash = mix(hash, self.edge_fingerprint(id));
            }
        }
        for vertex in self.explore(shape, NodeKind::Vertex) {
            hash = placed(hash, &vertex);
            if let Ok(id) = super::VertexId::try_from(vertex.node) {
                hash = self.mix_slot(hash, id);
            }
        }
        hash
    }

    /// A fingerprint of a face and what it is made of (its surface, wires, coedges with
    /// their pcurves, and edges): another one after any of them is edited.
    pub(crate) fn face_fingerprint(&self, face: FaceId) -> u64 {
        let mut hash = self.mix_slot(0, face);
        let Some(record) = self.faces.get(face) else { return hash };
        hash = self.mix_slot(hash, record.surface);
        for &wire_use in &record.wires {
            hash = self.mix_slot(hash, wire_use);
            let Some(u) = self.wire_uses.get(wire_use) else { continue };
            hash = self.mix_slot(hash, u.wire);
            let Some(wire) = self.wires.get(u.wire) else { continue };
            for &coedge in &wire.coedges {
                hash = self.mix_slot(hash, coedge);
                if let Some(c) = self.coedges.get(coedge) {
                    hash = self.mix_slot(hash, c.pcurve);
                    hash = mix(hash, self.edge_fingerprint(c.edge));
                }
            }
        }
        hash
    }
}

#[cfg(test)]
mod tests {
    use crate::math::Point;
    use crate::primitives::Cuboid;
    use crate::topology::{EdgeId, FaceId, Model, NodeKind, Shape};

    #[test]
    fn edits_below_a_face_change_its_fingerprint_and_its_derived_values() {
        let mut model = Model::new();
        let cuboid = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let face = model.explore(cuboid, NodeKind::Face).next().unwrap();
        let id = FaceId::try_from(face.node).unwrap();
        let edge = model.explore(face, NodeKind::Edge).next().unwrap();
        let edge_id = EdgeId::try_from(edge.node).unwrap();
        let vertex = model.explore(edge, NodeKind::Vertex).next().unwrap();
        let computed = std::cell::Cell::new(0);
        let value = |model: &Model| {
            *model.derived(face.node, model.face_fingerprint(id), || {
                computed.set(computed.get() + 1);
                computed.get()
            })
        };
        let unchanged = model.face_fingerprint(id);
        assert_eq!((value(&model), value(&model)), (1, 1), "kept while unchanged");
        // Reading changes nothing; editing a vertex, an edge's curve or a pcurve does.
        let _ = model.surface(face);
        assert_eq!(model.face_fingerprint(id), unchanged);
        let p = model.point(vertex);
        model.set_vertex_point(Shape::from(vertex.node), Point::new(p.x, p.y, p.z), 1e-6).unwrap();
        assert_ne!(model.face_fingerprint(id), unchanged);
        assert_eq!(value(&model), 2, "computed again after an edit below the face");
        let edge_print = model.edge_fingerprint(edge_id);
        let (first, last) = model.edge_range(edge);
        model.set_edge_range(edge, first, last).unwrap();
        assert_ne!(model.edge_fingerprint(edge_id), edge_print);
        assert_eq!(value(&model), 3);
        let wire = model.children(face).next().unwrap();
        let coedge = crate::topology::CoEdgeId::try_from(model.children(wire).next().unwrap().node).unwrap();
        let before = model.face_fingerprint(id);
        let (pcurve, ..) = model.effective_pcurve(coedge).unwrap();
        model.set_pcurve(coedge, pcurve).unwrap();
        assert_ne!(model.face_fingerprint(id), before);
    }
}
