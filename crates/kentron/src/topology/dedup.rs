//! Sharing equal geometry between entities and merging identical topology
//! ([`Model::deduplicate_with`]).

use std::collections::HashMap;

use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point2, Transform2, Transformable2, Vector2};

use super::entities::*;
use super::id::{Id, NodeId, NodeKind};
use super::slots::{EdgeCurve, FaceSurface};
use super::{DeduplicateOptions, DeduplicateReport, Model};

/// Parameters sampled to compare geometry: the range, or `[-10, 10]` where it is infinite.
fn samples(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let (lo, hi) = if lo.is_finite() && hi.is_finite() && lo.abs() < 1e50 && hi.abs() < 1e50 { (lo, hi) } else { (-10.0, 10.0) };
    (0..n).map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64).collect()
}

/// True if two edge curves are the same curve over the same range, with the same
/// parametrization.
fn same_curve(a: &EdgeCurve, b: &EdgeCurve, linear: f64) -> bool {
    let close = |x: f64, y: f64| (x - y).abs() <= crate::tolerance::PARAMETRIC * (1.0 + x.abs());
    if !close(a.first, b.first) || !close(a.last, b.last) {
        return false;
    }
    samples(a.first, a.last, 9).into_iter().all(|t| a.curve.point(t).distance(b.curve.point(t)) <= linear)
}

/// How a face's parameters map onto those of an equal surface: `Some(identity)` for the
/// same parametrization, the in-plane rigid motion between two frames of one plane, or
/// `None` if the surfaces differ.
fn same_surface(a: &Surface, b: &Surface, linear: f64, angular: f64) -> Option<Transform2> {
    if let (Surface::Plane(p), Surface::Plane(q)) = (a, b) {
        let (fp, fq) = (p.frame(), q.frame());
        let normal = crate::math::Vector::from(fp.z()).dot(crate::math::Vector::from(fq.z()));
        if normal < 1.0 - angular.max(1e-12) || q.distance(fp.origin()) > linear {
            return None;
        }
        // (u, v) on p to (u', v') on q: a rotation by the angle between the X axes and a
        // shift by p's origin seen from q.
        let (xp, xq, yq) = (crate::math::Vector::from(fp.x()), crate::math::Vector::from(fq.x()), crate::math::Vector::from(fq.y()));
        let angle = xp.dot(yq).atan2(xp.dot(xq));
        let shift = fp.origin() - fq.origin();
        let t = Transform2::translation(Vector2::new(shift.dot(xq), shift.dot(yq))) * Transform2::rotation(Point2::ORIGIN, angle);
        return Some(t);
    }
    let (au, av) = (a.is_u_periodic(), a.is_v_periodic());
    if au != b.is_u_periodic() || av != b.is_v_periodic() {
        return None;
    }
    let ((u0, u1, v0, v1), (bu0, bu1, bv0, bv1)) = (a.bounds(), b.bounds());
    let close = |x: f64, y: f64| x == y || (x - y).abs() <= crate::tolerance::PARAMETRIC * (1.0 + x.abs());
    if !(close(u0, bu0) && close(u1, bu1) && close(v0, bv0) && close(v1, bv1)) {
        return None;
    }
    for u in samples(u0, u1, 5) {
        for v in samples(v0, v1, 5) {
            if a.point(u, v).distance(b.point(u, v)) > linear {
                return None;
            }
        }
    }
    Some(Transform2::IDENTITY)
}

impl Model {
    pub(super) fn deduplicate_impl(&mut self, options: DeduplicateOptions) -> DeduplicateReport {
        let mut report = DeduplicateReport::default();
        let (linear, angular) = (options.linear_tolerance, options.angular_tolerance);
        let apply = !options.analyze_only;

        // Surfaces: each face's slot and its class's (the first face's) with the map of
        // parameters onto it.
        let mut classes: Vec<Id<FaceSurface>> = Vec::new();
        let mut surface_of: Vec<(FaceId, Id<FaceSurface>, Id<FaceSurface>, Transform2)> = Vec::new();
        let faces: Vec<(FaceId, Id<FaceSurface>)> = self.faces.iter().filter(|(_, f)| self.get(f.surface).is_some()).map(|(id, f)| (id, f.surface)).collect();
        for &(face, slot) in &faces {
            let surface = &self.face_surfaces[slot].surface;
            let found = classes.iter().find_map(|&c| {
                if c == slot {
                    return Some((c, Transform2::IDENTITY));
                }
                same_surface(surface, &self.face_surfaces[c].surface, linear, angular).map(|t| (c, t))
            });
            let (class, map) = found.unwrap_or_else(|| {
                classes.push(slot);
                (slot, Transform2::IDENTITY)
            });
            surface_of.push((face, slot, class, map));
        }
        report.distinct_surfaces = classes.len();
        // (Rewrites made: none when only analyzing.)
        report.shared_surfaces = if apply { surface_of.iter().filter(|(_, slot, class, _)| slot != class).count() } else { 0 };

        // Curves.
        let mut curve_classes: Vec<Id<EdgeCurve>> = Vec::new();
        let mut curve_of: Vec<(EdgeId, Id<EdgeCurve>, Id<EdgeCurve>)> = Vec::new();
        let edges: Vec<(EdgeId, Id<EdgeCurve>)> = self.edges.iter().filter(|(_, e)| self.get(e.curve).is_some()).map(|(id, e)| (id, e.curve)).collect();
        for &(edge, slot) in &edges {
            let curve = &self.edge_curves[slot];
            let class = curve_classes.iter().copied().find(|&c| c == slot || same_curve(curve, &self.edge_curves[c], linear));
            let class = class.unwrap_or_else(|| {
                curve_classes.push(slot);
                slot
            });
            curve_of.push((edge, slot, class));
        }
        report.distinct_curves = curve_classes.len();
        report.shared_curves = if apply { curve_of.iter().filter(|(_, slot, class)| slot != class).count() } else { 0 };

        if apply {
            for &(face, slot, class, map) in &surface_of {
                if slot == class {
                    continue;
                }
                if !map.is_identity() {
                    self.map_face_parameters(face, &map);
                }
                self.faces[face].surface = class;
            }
            for &(edge, slot, class) in &curve_of {
                if slot != class {
                    self.edges[edge].curve = class;
                }
            }
            // Slots no live record uses any more.
            let used: std::collections::HashSet<_> = self.faces.iter().map(|(_, f)| f.surface).collect();
            for &(_, slot, _, _) in &surface_of {
                if !used.contains(&slot) {
                    self.face_surfaces.remove(slot);
                }
            }
            let used: std::collections::HashSet<_> = self.edges.iter().map(|(_, e)| e.curve).collect();
            for &(_, slot, _) in &curve_of {
                if !used.contains(&slot) {
                    self.edge_curves.remove(slot);
                }
            }
        }

        if options.merge_entities {
            let surface_class: HashMap<FaceId, Id<FaceSurface>> = surface_of.iter().map(|&(f, _, c, _)| (f, c)).collect();
            let curve_class: HashMap<EdgeId, Id<EdgeCurve>> = curve_of.iter().map(|&(e, _, c)| (e, c)).collect();
            let merges = self.entity_merges(linear, &surface_class, &curve_class);
            for (kind, pairs) in &merges {
                report.merged[*kind as usize] += pairs.len();
            }
            if apply && merges.iter().any(|(_, p)| !p.is_empty()) {
                report.merge_applied = true;
                let record = options.record_history && self.history.as_ref().is_some_and(|h| h.is_enabled());
                for (_, pairs) in merges {
                    for (old, new) in pairs {
                        self.merge_node(old, new);
                        if record {
                            if let Some(history) = &mut self.history {
                                history.record_replaced("deduplicate", old, new);
                                report.history_events += 1;
                            }
                        }
                    }
                }
            }
        }
        report
    }

    /// Re-expresses the stored pcurves of a face's coedges through `map` (onto an equal
    /// plane with another frame).
    fn map_face_parameters(&mut self, face: FaceId, map: &Transform2) {
        for coedge in self.face_coedges(face) {
            let slot = self.coedges[coedge].pcurve;
            if let Some(p) = self.pcurves.get_mut(slot) {
                p.curve = p.curve.transformed(map);
            }
            let slot = self.coedges[coedge].pcurve_polyline;
            if let Some(p) = self.pcurve_polylines.get_mut(slot) {
                let nodes: Vec<Point2> = p.polyline.nodes().iter().map(|q| q.transformed(map)).collect();
                p.polyline = crate::mesh::Polyline2::new(nodes);
            }
        }
    }

    /// The (merged, kept) pairs of identical vertices, edges, wires and faces, each kind
    /// judged with the earlier kinds merged.
    fn entity_merges(
        &self,
        linear: f64,
        surface_class: &HashMap<FaceId, Id<FaceSurface>>,
        curve_class: &HashMap<EdgeId, Id<EdgeCurve>>,
    ) -> Vec<(NodeKind, Vec<(NodeId, NodeId)>)> {
        // Vertices at the same point.
        let mut vertex: HashMap<VertexId, VertexId> = HashMap::new();
        let mut kept: Vec<VertexId> = Vec::new();
        let mut vertex_pairs = Vec::new();
        for (id, v) in self.vertices.iter() {
            match kept.iter().copied().find(|&k| self.vertices[k].point.distance(v.point) <= linear) {
                Some(k) => {
                    vertex.insert(id, k);
                    vertex_pairs.push((id.into(), k.into()));
                }
                None => {
                    vertex.insert(id, id);
                    kept.push(id);
                }
            }
        }
        let canon_v = |v: VertexId| vertex.get(&v).copied().unwrap_or(v);
        // Edges on the same curve between the same vertices.
        let mut edge: HashMap<EdgeId, EdgeId> = HashMap::new();
        let mut edge_keys: HashMap<(Id<EdgeCurve>, VertexId, VertexId), EdgeId> = HashMap::new();
        let mut edge_pairs = Vec::new();
        for (id, _) in self.edges.iter() {
            let Some(&curve) = curve_class.get(&id) else {
                edge.insert(id, id);
                continue;
            };
            let (a, b) = self.edge_end_vertices(id);
            let key = (curve, canon_v(a), canon_v(b));
            match edge_keys.get(&key) {
                Some(&k) => {
                    edge.insert(id, k);
                    edge_pairs.push((id.into(), k.into()));
                }
                None => {
                    edge_keys.insert(key, id);
                    edge.insert(id, id);
                }
            }
        }
        let canon_e = |e: EdgeId| edge.get(&e).copied().unwrap_or(e);
        // Wires of the same coedges.
        let mut wire: HashMap<WireId, WireId> = HashMap::new();
        let mut wire_keys: HashMap<Vec<(EdgeId, Orientation)>, WireId> = HashMap::new();
        let mut wire_pairs = Vec::new();
        for (id, w) in self.wires.iter() {
            let key: Vec<(EdgeId, Orientation)> = w.coedges.iter().filter_map(|&c| self.get(c)).map(|c| (canon_e(c.edge), c.orientation)).collect();
            if key.is_empty() {
                wire.insert(id, id);
                continue;
            }
            match wire_keys.get(&key) {
                Some(&k) => {
                    wire.insert(id, k);
                    wire_pairs.push((id.into(), k.into()));
                }
                None => {
                    wire_keys.insert(key, id);
                    wire.insert(id, id);
                }
            }
        }
        let canon_w = |w: WireId| wire.get(&w).copied().unwrap_or(w);
        // Faces on the same surface with the same wires.
        let mut face_keys: HashMap<(Id<FaceSurface>, Vec<(WireId, Orientation)>), FaceId> = HashMap::new();
        let mut face_pairs = Vec::new();
        for (id, f) in self.faces.iter() {
            let Some(&surface) = surface_class.get(&id) else { continue };
            let key = (surface, f.wires.iter().filter_map(|&u| self.get(u)).map(|u| (canon_w(u.wire), u.orientation)).collect::<Vec<_>>());
            match face_keys.get(&key) {
                Some(&k) => face_pairs.push((id.into(), k.into())),
                None => {
                    face_keys.insert(key, id);
                }
            }
        }
        vec![(NodeKind::Vertex, vertex_pairs), (NodeKind::Edge, edge_pairs), (NodeKind::Wire, wire_pairs), (NodeKind::Face, face_pairs)]
    }

    /// Merges `old` into the identical `new`: its uses go to `new` and it is removed (with
    /// the coedges of a merged wire), leaving any geometry it shares with `new` alone.
    fn merge_node(&mut self, old: NodeId, new: NodeId) {
        match old.kind {
            NodeKind::Edge => {
                let id = EdgeId::from_index(old.index);
                self.edges[id].curve = Id::INVALID;
            }
            NodeKind::Face => {
                let id = FaceId::from_index(old.index);
                self.faces[id].surface = Id::INVALID;
            }
            _ => {}
        }
        let coedges = match old.kind {
            NodeKind::Wire => self.wires[WireId::from_index(old.index)].coedges.clone(),
            _ => Vec::new(),
        };
        if self.replace_node_impl(old, new).is_ok() {
            for c in coedges {
                self.remove_node_impl(c.into());
            }
        }
    }
}
