//! Navigating between topology entities by id: upward (parents, ancestors), sideways
//! (related entities) and per-entity derived facts.
//!
//! Downward navigation needs no API: entity records hold their children's use ids
//! (`model[shell].faces`, then `model[face_use].face`). These queries work on definitions
//! and ignore placements; use [`Shape`](super::Shape) queries for placed geometry.

use super::entities::{CoEdgeId, EdgeId, FaceId, Orientation, ShellId, VertexId, WireId};
use super::id::{Id, Node, NodeId};
use super::Model;

/// How a node returned by [`Model::related`] relates to the queried node.
#[doc(alias = "BRepGraph_RelatedIterator::RelationKind")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Relation {
    /// Face -> an edge bounding it.
    BoundaryEdge,
    /// Face -> a face sharing an edge with it.
    AdjacentFace,
    /// Face -> its outer boundary wire.
    OuterWire,
    /// Edge -> a face using it.
    ReferencedByFace,
    /// Edge -> its start or end vertex.
    IncidentVertex,
    /// Wire -> one of its coedges.
    WireCoEdge,
    /// Wire or coedge -> the face it bounds.
    OwningFace,
    /// Vertex -> an edge ending at it.
    IncidentEdge,
    /// Coedge -> its edge.
    ParentEdge,
    /// Coedge -> the other coedge of a seam edge on the same face.
    SeamPair,
}

impl Model {
    /// The nodes directly containing `node` (a face's shells, an edge's coedges, a
    /// coedge's wire, ...), skipping removed ones, in creation order.
    #[doc(alias = "BRepGraph_ReverseIterator")]
    pub fn parents(&self, node: NodeId) -> Vec<NodeId> {
        let mut out: Vec<NodeId> = Vec::new();
        for link in self.parent_links(node) {
            if !out.contains(&link.parent) {
                out.push(link.parent);
            }
        }
        out
    }

    /// The distinct `T` nodes containing `node` at any depth, nearest first (e.g. the
    /// faces using an edge, the solids containing a face).
    #[doc(alias = "BRepGraph_ParentExplorer")]
    pub fn ancestors_of<T: Node>(&self, node: impl Into<NodeId>) -> Vec<Id<T>> {
        let mut found: Vec<Id<T>> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut queue: std::collections::VecDeque<NodeId> = self.parents(node.into()).into();
        while let Some(n) = queue.pop_front() {
            if !seen.insert(n) {
                continue;
            }
            if let Ok(id) = Id::<T>::try_from(n) {
                found.push(id);
            }
            queue.extend(self.parents(n));
        }
        found
    }

    /// For each of `edges`, the distinct faces using it (as [`Self::ancestors_of`] finds
    /// them, in no particular order), found in one pass over the coedges and one over the
    /// wire uses instead of a walk up from every edge.
    pub fn faces_around(&self, edges: &[EdgeId]) -> Vec<Vec<FaceId>> {
        let mut wires_of: crate::hashing::FastMap<EdgeId, Vec<WireId>> = Default::default();
        for (_, c) in self.coedges.iter() {
            if edges.contains(&c.edge) && self.is_live(c.wire.into()) {
                wires_of.entry(c.edge).or_default().push(c.wire);
            }
        }
        let mut faces_of: crate::hashing::FastMap<WireId, Vec<FaceId>> = Default::default();
        for (_, u) in self.wire_uses.iter() {
            if wires_of.values().any(|w| w.contains(&u.wire)) && self.is_live(u.face.into()) {
                faces_of.entry(u.wire).or_default().push(u.face);
            }
        }
        edges
            .iter()
            .map(|e| {
                let mut faces: Vec<FaceId> = Vec::new();
                for w in wires_of.get(e).into_iter().flatten() {
                    for &f in faces_of.get(w).into_iter().flatten() {
                        if !faces.contains(&f) {
                            faces.push(f);
                        }
                    }
                }
                faces
            })
            .collect()
    }

    /// The distinct `T` nodes below `node`, in depth-first order.
    pub fn descendants_of<T: Node>(&self, node: impl Into<NodeId>) -> Vec<Id<T>> {
        let node = node.into();
        let mut found: Vec<Id<T>> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        if !self.is_live(node) {
            return found;
        }
        let mut stack: Vec<NodeId> = self.child_nodes(node).into_iter().rev().collect();
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            if let Ok(id) = Id::<T>::try_from(n) {
                found.push(id);
            }
            stack.extend(self.child_nodes(n).into_iter().rev());
        }
        found
    }

    /// The nodes topologically related to a face, wire, coedge, edge or vertex, with the
    /// kind of relation (not deduplicated). Containers (shells, solids, compounds,
    /// products, instances) have no relations.
    #[doc(alias = "BRepGraph_RelatedIterator")]
    pub fn related(&self, node: NodeId) -> Vec<(Relation, NodeId)> {
        use super::id::NodeKind;
        let mut out: Vec<(Relation, NodeId)> = Vec::new();
        if !self.is_live(node) {
            return out;
        }
        let edge_faces = |e: EdgeId| -> Vec<FaceId> {
            self.coedges.iter().filter(|(_, c)| c.edge == e && self.get(c.face).is_some()).map(|(_, c)| c.face).collect()
        };
        match node.kind {
            NodeKind::Face => {
                let face = FaceId::from_index(node.index);
                if let Some(w) = self.outer_wire(face) {
                    out.push((Relation::OuterWire, w.into()));
                }
                let mut adjacent: Vec<FaceId> = Vec::new();
                for c in self.face_coedges(face) {
                    let edge = self.coedges[c].edge;
                    out.push((Relation::BoundaryEdge, edge.into()));
                    for f in edge_faces(edge) {
                        if f != face && !adjacent.contains(&f) {
                            adjacent.push(f);
                        }
                    }
                }
                out.extend(adjacent.into_iter().map(|f| (Relation::AdjacentFace, f.into())));
            }
            NodeKind::Edge => {
                let edge = EdgeId::from_index(node.index);
                for f in edge_faces(edge) {
                    out.push((Relation::ReferencedByFace, f.into()));
                }
                let record = &self.edges[edge];
                for u in [record.start, record.end] {
                    if let Some(u) = self.get(u) {
                        out.push((Relation::IncidentVertex, u.vertex.into()));
                    }
                }
            }
            NodeKind::Wire => {
                let wire = WireId::from_index(node.index);
                for &c in &self.wires[wire].coedges {
                    out.push((Relation::WireCoEdge, c.into()));
                }
                if let Some(f) = self.face_of_wire(wire) {
                    out.push((Relation::OwningFace, f.into()));
                }
            }
            NodeKind::CoEdge => {
                let coedge = CoEdgeId::from_index(node.index);
                let record = &self.coedges[coedge];
                out.push((Relation::ParentEdge, record.edge.into()));
                if self.get(record.face).is_some() {
                    out.push((Relation::OwningFace, record.face.into()));
                }
                if let Some(pair) = self.seam_pair(coedge) {
                    out.push((Relation::SeamPair, pair.into()));
                }
            }
            NodeKind::Vertex => {
                let vertex = VertexId::from_index(node.index);
                for (_, u) in self.vertex_uses.iter() {
                    if u.vertex == vertex && self.get(u.edge).is_some() {
                        out.push((Relation::IncidentEdge, u.edge.into()));
                    }
                }
            }
            _ => {}
        }
        out
    }

    /// The outer boundary wire of a face (its first wire), if it has wires.
    #[doc(alias = "BRepGraph_Tool::Face::OuterWire")]
    pub fn outer_wire(&self, face: FaceId) -> Option<WireId> {
        let face = self.get(face)?;
        face.wires.iter().filter_map(|&u| self.get(u)).map(|u| u.wire).find(|&w| self.get(w).is_some())
    }

    /// The (u, v) bounds of a face, from its boundary or, without one, its surface's domain.
    #[doc(alias = "BRepGraph_Tool::Face::Bounds")]
    pub fn face_bounds(&self, face: FaceId) -> (f64, f64, f64, f64) {
        use crate::geometry::ParametricCurve;
        let mut bounds = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for coedge in self.face_coedges(face) {
            let Some((pcurve, first, last)) = self.effective_pcurve(coedge) else { continue };
            // The ends, and enough samples to catch a curved pcurve's extent.
            for i in 0..=32 {
                let p = pcurve.point(first + (last - first) * i as f64 / 32.0);
                bounds = (bounds.0.min(p.x), bounds.1.max(p.x), bounds.2.min(p.y), bounds.3.max(p.y));
            }
        }
        if bounds.0 > bounds.1 {
            if let Some(surface) = self.face_surface(face) {
                return crate::geometry::ParametricSurface::bounds(surface);
            }
        }
        bounds
    }

    /// The face a wire bounds, if any.
    #[doc(alias = "BRepGraph_Tool::Wire::FaceOf")]
    pub fn face_of_wire(&self, wire: WireId) -> Option<FaceId> {
        self.wire_uses.iter().find(|(_, u)| u.wire == wire && self.get(u.face).is_some()).map(|(_, u)| u.face)
    }

    /// The start and end vertices of an edge (in the edge's own direction).
    #[doc(alias = "BRepGraph_Tool::Edge::StartVertexId")]
    #[doc(alias = "BRepGraph_Tool::Edge::EndVertexId")]
    pub fn edge_end_vertices(&self, edge: EdgeId) -> (VertexId, VertexId) {
        let record = &self.edges[edge];
        let vertex = |u| self.get(u).map_or(Id::INVALID, |u: &super::entities::VertexUse| u.vertex);
        (vertex(record.start), vertex(record.end))
    }

    /// An edge from `start` to `end` (or, with `either_direction`, from `end` to `start`).
    #[doc(alias = "BRepGraph_Tool::Edge::FindByVertices")]
    pub fn find_edge(&self, start: VertexId, end: VertexId, either_direction: bool) -> Option<EdgeId> {
        self.edges.iter().map(|(id, _)| id).find(|&e| {
            let (a, b) = self.edge_end_vertices(e);
            (a == start && b == end) || (either_direction && a == end && b == start)
        })
    }

    /// The coedge of `edge` on `face` (the first of the two for a seam edge), if any.
    #[doc(alias = "BRepGraph_Tool::Edge::FindPCurveCoEdgeId")]
    pub fn find_coedge(&self, edge: EdgeId, face: FaceId) -> Option<CoEdgeId> {
        self.face_coedges(face).into_iter().find(|&c| self.get(c).is_some_and(|c| c.edge == edge))
    }

    /// The coedge of `edge` on `face` with the given orientation: on a seam edge, the two
    /// coedges differ by orientation.
    pub fn find_coedge_oriented(
        &self,
        edge: EdgeId,
        face: FaceId,
        orientation: Orientation,
    ) -> Option<CoEdgeId> {
        self.face_coedges(face)
            .into_iter()
            .find(|&c| self.get(c).is_some_and(|c| c.edge == edge && c.orientation == orientation))
    }

    /// True if exactly one live face uses `edge` (a free edge of an open shell). False for
    /// a removed edge.
    #[doc(alias = "BRepGraph_Tool::Edge::IsBoundary")]
    pub fn is_boundary_edge(&self, edge: EdgeId) -> bool {
        self.get(edge).is_some() && self.ancestors_of::<super::entities::Face>(edge).len() == 1
    }

    /// True if exactly two live faces use `edge`. False for a removed edge.
    #[doc(alias = "BRepGraph_Tool::Edge::IsManifold")]
    pub fn is_manifold_edge(&self, edge: EdgeId) -> bool {
        self.get(edge).is_some() && self.ancestors_of::<super::entities::Face>(edge).len() == 2
    }

    /// True if the wire's coedges form a closed loop: each coedge ends where the next one
    /// starts, and the last where the first one starts.
    #[doc(alias = "BRepGraph_Tool::Wire::IsClosed")]
    pub fn is_closed_wire(&self, wire: WireId) -> bool {
        // The vertices each coedge runs from and to, following its orientation.
        let ends = |c: CoEdgeId| -> Option<(VertexId, VertexId)> {
            let coedge = self.get(c)?;
            let edge = self.get(coedge.edge)?;
            let (a, b) = (self.get(edge.start)?.vertex, self.get(edge.end)?.vertex);
            Some(if coedge.orientation == Orientation::Reversed { (b, a) } else { (a, b) })
        };
        let Some(wire) = self.get(wire) else { return false };
        let Some(ends) = wire.coedges.iter().map(|&c| ends(c)).collect::<Option<Vec<_>>>() else { return false };
        !ends.is_empty() && (0..ends.len()).all(|i| ends[i].1 == ends[(i + 1) % ends.len()].0)
    }

    /// True if each face of the shell meets other faces along all of its edges, enclosing
    /// a volume.
    #[doc(alias = "BRepGraph_Tool::Shell::IsClosed")]
    pub fn is_closed_shell(&self, shell: ShellId) -> bool {
        // Every (non-degenerate) edge is used as often in one direction as in the other.
        if self.get(shell).is_none() {
            return false;
        }
        let mut uses: std::collections::HashMap<NodeId, i64> = std::collections::HashMap::new();
        for edge in self.explore(crate::topology::Shape::from(shell), crate::topology::NodeKind::Edge) {
            if self.is_degenerate(edge) {
                continue;
            }
            let count = uses.entry(edge.node).or_default();
            match edge.orientation {
                Orientation::Forward => *count += 1,
                Orientation::Reversed => *count -= 1,
                _ => {}
            }
        }
        !uses.is_empty() && uses.values().all(|&c| c == 0)
    }

    /// True if `face` uses `edge` twice, with opposite orientations (a seam).
    #[doc(alias = "BRepGraph_Tool::Edge::IsSeamOnFace")]
    pub fn is_seam_on_face(&self, edge: EdgeId, face: FaceId) -> bool {
        let uses: Vec<_> = self
            .face_coedges(face)
            .into_iter()
            .filter_map(|c| self.get(c))
            .filter(|c| c.edge == edge)
            .map(|c| c.orientation)
            .collect();
        uses.len() == 2 && uses[0] == uses[1].reversed() && uses[0] != uses[1]
    }

    /// True if the coedge's edge is a seam of its face (used twice by it).
    #[doc(alias = "BRepGraph_Tool::CoEdge::IsSeam")]
    pub fn is_seam(&self, coedge: CoEdgeId) -> bool {
        self.get(coedge).is_some_and(|c| self.is_seam_on_face(c.edge, c.face))
    }
}
