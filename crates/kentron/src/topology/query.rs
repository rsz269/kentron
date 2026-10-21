//! Navigating topology and reading the geometry attached to it.
//!
//! Queries returning owned geometry ([`Model::point`], [`Model::edge_curve`]) apply the
//! shape's placement: the point of a placed vertex is where the vertex is in the assembly.
//! Queries borrowing stored geometry ([`Model::surface`], [`Model::face_mesh`],
//! [`Model::edge_polyline`]) return it in the node's own frame; apply `shape.placement` to
//! get assembly coordinates.

use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Plane;
use crate::math::{Point, Transform, Transformable};
use crate::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};

use super::entities::{
    CoEdgeId, CompSolidId, CompoundId, EdgeId, FaceId, InstanceId, Orientation, ProductId, ShellId,
    SolidId, VertexId, WireId,
};

use super::id::{NodeId, NodeKind};
use super::shape::Shape;
use super::Model;

/// Iterator over the direct children of a shape, in order. Created by
/// [`Model::children`] (orientation and placement composed with the parent's) or
/// [`Model::children_local`] (as stored, relative to the parent).
#[doc(alias = "TopoDS_Iterator")]
#[derive(Debug, Clone)]
pub struct Children<'a> {
    model: &'a Model,
    parent: Shape,
    compose: bool,
    /// The parent is placed at the identity (as most are): its children's placements are
    /// their own, without a product each.
    unplaced: bool,
    next: usize,
}

impl Iterator for Children<'_> {
    type Item = Shape;

    fn next(&mut self) -> Option<Shape> {
        let child = self.model.child_at(self.parent.node, self.next)?;
        self.next += 1;
        Some(if self.compose {
            Shape {
                node: child.node,
                orientation: child.orientation.compose(self.parent.orientation),
                placement: if self.unplaced { child.placement } else { self.parent.placement * child.placement },
            }
        } else {
            child
        })
    }
}

impl Model {
    /// The coedges of all wires of `face`, in order.
    pub(super) fn face_coedges(&self, face: FaceId) -> Vec<CoEdgeId> {
        let Some(face) = self.get(face) else { return Vec::new() };
        face.wires
            .iter()
            .filter_map(|&w| self.get(w))
            .filter_map(|w| self.get(w.wire))
            .flat_map(|w| w.coedges.iter().copied())
            .collect()
    }

    /// The `i`-th child of `node` as stored (relative to it); `None` past the last one or
    /// for a missing node.
    pub(super) fn child_at(&self, node: NodeId, i: usize) -> Option<Shape> {
        let local = |node: NodeId, orientation: Orientation| Shape {
            node,
            orientation,
            placement: Transform::IDENTITY,
        };
        match node.kind {
            NodeKind::Vertex => None,
            NodeKind::Edge => {
                let edge = self.get(EdgeId::from_index(node.index))?;
                let use_id = *[edge.start, edge.end].get(i)?;
                let vertex_use = self.get(use_id)?;
                Some(local(vertex_use.vertex.into(), vertex_use.orientation))
            }
            NodeKind::CoEdge => {
                let coedge = self.get(CoEdgeId::from_index(node.index))?;
                (i == 0).then(|| local(coedge.edge.into(), coedge.orientation))
            }
            NodeKind::Wire => {
                let wire = self.get(WireId::from_index(node.index))?;
                wire.coedges.get(i).map(|&c| local(c.into(), Orientation::Forward))
            }
            NodeKind::Face => {
                let face = self.get(FaceId::from_index(node.index))?;
                let used = self.get(*face.wires.get(i)?)?;
                Some(local(used.wire.into(), used.orientation))
            }
            NodeKind::Shell => {
                let shell = self.get(ShellId::from_index(node.index))?;
                let used = self.get(*shell.faces.get(i)?)?;
                Some(local(used.face.into(), used.orientation))
            }
            NodeKind::Solid => {
                let solid = self.get(SolidId::from_index(node.index))?;
                let used = self.get(*solid.shells.get(i)?)?;
                Some(local(used.shell.into(), used.orientation))
            }
            NodeKind::CompSolid => {
                let comp_solid = self.get(CompSolidId::from_index(node.index))?;
                let used = self.get(*comp_solid.solids.get(i)?)?;
                Some(local(used.solid.into(), used.orientation))
            }
            NodeKind::Compound => {
                let compound = self.get(CompoundId::from_index(node.index))?;
                let child = self.get(*compound.children.get(i)?)?;
                Some(Shape { node: child.child, orientation: child.orientation, placement: child.placement })
            }
            NodeKind::Product => {
                // A part's root shape, then its instances.
                let product = self.get(ProductId::from_index(node.index))?;
                let i = match product.shape {
                    Some(shape) if i == 0 => return Some(shape),
                    Some(_) => i - 1,
                    None => i,
                };
                let child = self.get(*product.children.get(i)?)?;
                Some(Shape {
                    node: child.instance.into(),
                    orientation: Orientation::Forward,
                    placement: child.placement,
                })
            }
            NodeKind::Instance => {
                let instance = self.get(InstanceId::from_index(node.index))?;
                (i == 0).then(|| local(instance.target, Orientation::Forward))
            }
        }
    }
}

impl Model {
    /// The children of `node` it does not place (those of every kind of node but compounds
    /// and products), each with its orientation relative to `node`, passed to `f` in order,
    /// as [`Self::child_at`] has them; false, and none passed, for a compound or a product.
    #[inline]
    pub(super) fn local_children(&self, node: NodeId, mut f: impl FnMut(NodeId, Orientation)) -> bool {
        match node.kind {
            NodeKind::Vertex => {}
            NodeKind::Edge => {
                if let Some(edge) = self.get(EdgeId::from_index(node.index)) {
                    for use_id in [edge.start, edge.end] {
                        let Some(vertex_use) = self.get(use_id) else { break };
                        f(vertex_use.vertex.into(), vertex_use.orientation);
                    }
                }
            }
            NodeKind::CoEdge => {
                if let Some(coedge) = self.get(CoEdgeId::from_index(node.index)) {
                    f(coedge.edge.into(), coedge.orientation);
                }
            }
            NodeKind::Wire => {
                if let Some(wire) = self.get(WireId::from_index(node.index)) {
                    for &c in &wire.coedges {
                        f(c.into(), Orientation::Forward);
                    }
                }
            }
            NodeKind::Face => {
                if let Some(face) = self.get(FaceId::from_index(node.index)) {
                    for &w in &face.wires {
                        let Some(used) = self.get(w) else { break };
                        f(used.wire.into(), used.orientation);
                    }
                }
            }
            NodeKind::Shell => {
                if let Some(shell) = self.get(ShellId::from_index(node.index)) {
                    for &u in &shell.faces {
                        let Some(used) = self.get(u) else { break };
                        f(used.face.into(), used.orientation);
                    }
                }
            }
            NodeKind::Solid => {
                if let Some(solid) = self.get(SolidId::from_index(node.index)) {
                    for &u in &solid.shells {
                        let Some(used) = self.get(u) else { break };
                        f(used.shell.into(), used.orientation);
                    }
                }
            }
            NodeKind::CompSolid => {
                if let Some(comp_solid) = self.get(CompSolidId::from_index(node.index)) {
                    for &u in &comp_solid.solids {
                        let Some(used) = self.get(u) else { break };
                        f(used.solid.into(), used.orientation);
                    }
                }
            }
            NodeKind::Instance => {
                if let Some(instance) = self.get(InstanceId::from_index(node.index)) {
                    f(instance.target, Orientation::Forward);
                }
            }
            NodeKind::Compound | NodeKind::Product => return false,
        }
        true
    }

    /// The direct children of `shape`, with orientation and placement composed with
    /// `shape`'s: a solid's shells, a face's wires, an edge's two vertices, ...
    pub fn children(&self, shape: impl Into<Shape>) -> Children<'_> {
        let parent = shape.into();
        Children { model: self, parent, compose: true, unplaced: parent.placement == Transform::IDENTITY, next: 0 }
    }

    /// Like [`Model::children`], but with each child's orientation and placement relative
    /// to `shape`.
    pub fn children_local(&self, shape: impl Into<Shape>) -> Children<'_> {
        Children { model: self, parent: shape.into(), compose: false, unplaced: true, next: 0 }
    }

    /// The distinct `kind` sub-shapes below `root`, each with the distinct
    /// `ancestor_kind` shapes below `root` that contain it, in order of first occurrence.
    #[doc(alias = "TopExp::MapShapesAndAncestors")]
    pub fn subshape_ancestors(
        &self,
        root: impl Into<Shape>,
        kind: NodeKind,
        ancestor_kind: NodeKind,
    ) -> Vec<(Shape, Vec<Shape>)> {
        let root = root.into();
        let mut out: Vec<(Shape, Vec<Shape>)> = self.unique_subshapes(root, kind).into_iter().map(|s| (s, Vec::new())).collect();
        let same = |a: &Shape, b: &Shape| a.node == b.node && a.placement == b.placement;
        for ancestor in self.unique_subshapes(root, ancestor_kind) {
            for found in self.unique_subshapes(ancestor, kind) {
                if let Some(entry) = out.iter_mut().find(|(s, _)| same(s, &found)) {
                    if !entry.1.iter().any(|a| same(a, &ancestor)) {
                        entry.1.push(ancestor);
                    }
                }
            }
        }
        out
    }

    /// The point of a vertex. Panics if `vertex` is not a vertex.
    #[doc(alias = "BRep_Tool::Pnt")]
    pub fn point(&self, vertex: impl Into<Shape>) -> Point {
        let vertex = vertex.into();
        let id = VertexId::try_from(vertex.node).unwrap_or_else(|n| panic!("{n:?} is not a vertex"));
        let point = self.get(id).unwrap_or_else(|| panic!("no {id:?}")).point;
        point.transformed(&vertex.placement)
    }

    /// The tolerance of a vertex, edge or face. Panics for other kinds.
    #[doc(alias = "BRep_Tool::Tolerance")]
    pub fn tolerance(&self, shape: impl Into<Shape>) -> f64 {
        let node = shape.into().node;
        let tolerance = match node.kind {
            NodeKind::Vertex => self.get(VertexId::from_index(node.index)).map(|v| v.tolerance),
            NodeKind::Edge => self.get(EdgeId::from_index(node.index)).map(|e| e.tolerance),
            NodeKind::Face => self.get(FaceId::from_index(node.index)).map(|f| f.tolerance),
            kind => panic!("a {kind:?} has no tolerance"),
        };
        tolerance.unwrap_or_else(|| panic!("no {node:?}"))
    }

    /// The 3D curve of an edge with its parameter range; `None` for an edge without one
    /// (degenerate edges).
    #[doc(alias = "BRep_Tool::Curve")]
    pub fn edge_curve(&self, edge: impl Into<Shape>) -> Option<(Curve, f64, f64)> {
        let edge = edge.into();
        let id = EdgeId::try_from(edge.node).unwrap_or_else(|n| panic!("{n:?} is not an edge"));
        let slot = self.get(self.get(id)?.curve)?;
        let curve = if edge.placement == Transform::IDENTITY {
            slot.curve.clone()
        } else {
            slot.curve.transformed(&edge.placement)
        };
        Some((curve, slot.first, slot.last))
    }

    /// Whether the edge has a curve (as [`Self::edge_curve`] would give, without it).
    pub(crate) fn has_edge_curve(&self, edge: impl Into<Shape>) -> bool {
        let edge = edge.into();
        let id = EdgeId::try_from(edge.node).unwrap_or_else(|n| panic!("{n:?} is not an edge"));
        self.get(id).and_then(|record| self.get(record.curve)).is_some()
    }

    /// The parameter range of an edge.
    #[doc(alias = "BRep_Tool::Range")]
    pub fn edge_range(&self, edge: impl Into<Shape>) -> (f64, f64) {
        let edge = edge.into();
        let id = EdgeId::try_from(edge.node).unwrap_or_else(|n| panic!("{n:?} is not an edge"));
        let record = self.get(id).unwrap_or_else(|| panic!("no {id:?}"));
        let slot = self.get(record.curve).unwrap_or_else(|| panic!("{id:?} has no 3D curve"));
        (slot.first, slot.last)
    }

    /// The start and end vertices of an edge, in the direction given by the edge's
    /// orientation (swapped for a reversed edge).
    #[doc(alias = "TopExp::FirstVertex")]
    #[doc(alias = "TopExp::LastVertex")]
    #[doc(alias = "TopExp::Vertices")]
    pub fn edge_vertices(&self, edge: impl Into<Shape>) -> (Shape, Shape) {
        let edge = edge.into();
        let id = EdgeId::try_from(edge.node).unwrap_or_else(|n| panic!("{n:?} is not an edge"));
        let record = self.get(id).unwrap_or_else(|| panic!("no {id:?}"));
        let vertex = |u| {
            let u: &super::VertexUse = self.get(u).unwrap_or_else(|| panic!("{id:?} lacks a vertex"));
            Shape::new(u.vertex.into()).with_placement(edge.placement)
        };
        let (start, end) = match edge.orientation {
            Orientation::Reversed => (vertex(record.end), vertex(record.start)),
            _ => (vertex(record.start), vertex(record.end)),
        };
        // As TopExp::Vertices: the first vertex forward, the last reversed.
        (start, end.reversed())
    }

    /// A vertex shared by two edges, if any.
    #[doc(alias = "TopExp::CommonVertex")]
    pub fn common_vertex(&self, edge1: impl Into<Shape>, edge2: impl Into<Shape>) -> Option<Shape> {
        let (a, b) = (edge1.into(), edge2.into());
        let vertices = |e: Shape| self.unique_subshapes(e, NodeKind::Vertex);
        let others = vertices(b);
        vertices(a).into_iter().find(|v| others.iter().any(|w| w.node == v.node && w.placement == v.placement))
    }

    /// True if the edge starts and ends at the same vertex.
    #[doc(alias = "BRep_Tool::IsClosed")]
    pub fn is_closed_edge(&self, edge: impl Into<Shape>) -> bool {
        let (start, end) = self.edge_vertices(edge);
        start.node == end.node
    }

    /// True if the edge is collapsed to a point (e.g. at a sphere's pole).
    #[doc(alias = "BRep_Tool::Degenerated")]
    pub fn is_degenerate(&self, edge: impl Into<Shape>) -> bool {
        let edge = edge.into();
        let id = EdgeId::try_from(edge.node).unwrap_or_else(|n| panic!("{n:?} is not an edge"));
        let record = self.get(id).unwrap_or_else(|| panic!("no {id:?}"));
        let vertex = |u| self.get(u).map(|u: &super::VertexUse| u.vertex);
        self.get(record.curve).is_none() && vertex(record.start).is_some() && vertex(record.start) == vertex(record.end)
    }

    /// The (u, v) bounds of a face's boundary in its surface's parameter space:
    /// `(u_min, u_max, v_min, v_max)`.
    #[doc(alias = "BRepTools::UVBounds")]
    pub fn face_uv_bounds(&self, face: impl Into<Shape>) -> (f64, f64, f64, f64) {
        let face = face.into();
        let id = FaceId::try_from(face.node).unwrap_or_else(|n| panic!("{n:?} is not a face"));
        self.face_bounds(id)
    }

    /// The surface of a face (in the face's own coordinates: apply `face.placement` to
    /// get assembly coordinates), if it has one.
    pub fn surface(&self, face: impl Into<Shape>) -> Option<&crate::geometry::surface::Surface> {
        let face = face.into();
        let id = FaceId::try_from(face.node).unwrap_or_else(|n| panic!("{n:?} is not a face"));
        self.face_surface(id)
    }

    /// The mesh stored on a face, if any (in the face's own frame).
    #[doc(alias = "BRep_Tool::Triangulation")]
    pub fn face_mesh(&self, face: impl Into<Shape>) -> Option<&TriangleMesh> {
        let id = FaceId::try_from(face.into().node).ok()?;
        self.get(self.get(id)?.mesh).map(|slot| &slot.mesh)
    }

    /// The polyline stored on an edge, if any (in the edge's own frame).
    #[doc(alias = "BRep_Tool::Polygon3D")]
    pub fn edge_polyline(&self, edge: impl Into<Shape>) -> Option<&Polyline> {
        let id = EdgeId::try_from(edge.into().node).ok()?;
        self.get(self.get(id)?.polyline).map(|slot| &slot.polyline)
    }

    /// The stored pcurve of a coedge, if any.
    pub fn pcurve(&self, coedge: CoEdgeId) -> Option<&Curve2> {
        self.get(self.get(coedge)?.pcurve).map(|slot| &slot.curve)
    }

    /// The parameter range of a coedge on its pcurve.
    pub fn coedge_range(&self, coedge: CoEdgeId) -> (f64, f64) {
        let record = self.get(coedge).unwrap_or_else(|| panic!("no {coedge:?}"));
        match self.get(record.pcurve) {
            Some(slot) => (slot.first, slot.last),
            None => self.edge_range(record.edge),
        }
    }

    /// The pcurve of a coedge with its range: the stored one or, on a planar face without
    /// one, the projection of the edge's 3D curve.
    pub fn effective_pcurve(&self, coedge: CoEdgeId) -> Option<(Curve2, f64, f64)> {
        let record = self.get(coedge)?;
        if let Some(slot) = self.get(record.pcurve) {
            return Some((slot.curve.clone(), slot.first, slot.last));
        }
        let surface = self.get(self.get(record.face)?.surface)?;
        if !matches!(surface.surface, crate::geometry::surface::Surface::Plane(_)) {
            return None;
        }
        let edge = self.get(self.get(record.edge)?.curve)?;
        let curve = crate::geometry::project::analytic_pcurve(&edge.curve, &surface.surface)?;
        Some((curve, edge.first, edge.last))
    }

    /// A coedge's 2D polyline, if any.
    pub fn pcurve_polyline(&self, coedge: CoEdgeId) -> Option<&Polyline2> {
        self.get(self.get(coedge)?.pcurve_polyline).map(|slot| &slot.polyline)
    }

    /// True if the coedge's range equals its edge's range.
    #[doc(alias = "BRep_Tool::SameRange")]
    pub fn is_same_range(&self, coedge: CoEdgeId) -> bool {
        let Some(record) = self.get(coedge) else { return false };
        if self.get(record.edge).and_then(|e| self.get(e.curve)).is_none() {
            return true;
        }
        let (a, b) = self.coedge_range(coedge);
        let (c, d) = self.edge_range(record.edge);
        (a - c).abs() <= 1e-9 * (1.0 + c.abs()) && (b - d).abs() <= 1e-9 * (1.0 + d.abs())
    }

    /// True if the coedge's pcurve, mapped onto its face's surface, follows the edge's 3D
    /// curve at equal parameters (within the edge tolerance).
    #[doc(alias = "BRep_Tool::SameParameter")]
    pub fn is_same_parameter(&self, coedge: CoEdgeId) -> bool {
        use crate::geometry::{ParametricCurve, ParametricSurface};
        let Some(record) = self.get(coedge) else { return false };
        let Some(edge) = self.get(record.edge) else { return false };
        let Some(curve) = self.get(edge.curve) else { return true };
        // Same parameter at the same points needs the same range (as OCCT's flag does).
        if !self.is_same_range(coedge) {
            return false;
        }
        let Some((pcurve, _, _)) = self.effective_pcurve(coedge) else { return true };
        let Some(surface) = self.face_surface(record.face) else { return true };
        const SAMPLES: usize = 23;
        (0..=SAMPLES).all(|i| {
            let t = curve.first + (curve.last - curve.first) * i as f64 / SAMPLES as f64;
            let uv = pcurve.point(t);
            surface.point(uv.x, uv.y).distance(curve.curve.point(t)) <= edge.tolerance
        })
    }

    /// The other coedge of a seam edge on the same face.
    #[doc(alias = "BRepGraph_Tool::CoEdge::SeamPair")]
    pub fn seam_pair(&self, coedge: CoEdgeId) -> Option<CoEdgeId> {
        let record = self.get(coedge)?;
        self.face_coedges(record.face).into_iter().find(|&other| {
            other != coedge
                && self.get(other).is_some_and(|o| o.edge == record.edge && o.orientation != record.orientation)
        })
    }

    /// A coedge's discretization as nodes of its face's mesh, if any.
    #[doc(alias = "BRep_Tool::PolygonOnTriangulation")]
    pub fn coedge_mesh_polyline(&self, coedge: CoEdgeId) -> Option<&MeshPolyline> {
        self.get(self.get(coedge)?.mesh_polyline).map(|slot| &slot.polyline)
    }

    /// The edges of a wire in the order they connect, each oriented along the wire. Edges
    /// that do not connect to the chain are skipped.
    #[doc(alias = "BRepTools_WireExplorer")]
    pub fn ordered_edges(&self, wire: impl Into<Shape>) -> Vec<Shape> {
        let wire = wire.into();
        let Ok(id) = WireId::try_from(wire.node) else { return Vec::new() };
        let Some(record) = self.get(id) else { return Vec::new() };
        // Each use: the edge oriented along the wire, and its ends that way.
        let uses: Vec<(Shape, NodeId, NodeId)> = record
            .coedges
            .iter()
            .filter_map(|&c| {
                let coedge = self.get(c)?;
                let orientation = coedge.orientation.compose(wire.orientation);
                let edge = Shape { node: coedge.edge.into(), orientation, placement: wire.placement };
                let (a, b) = self.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
                let (start, end) = if orientation == Orientation::Reversed { (b.node, a.node) } else { (a.node, b.node) };
                Some((edge, start, end))
            })
            .collect();
        let Some(first) = uses.first() else { return Vec::new() };
        let mut used = vec![false; uses.len()];
        used[0] = true;
        let mut chain = std::collections::VecDeque::from([0usize]);
        let (head, mut tail) = (first.1, first.2);
        // Onwards while a use goes on from the end (a loop may pass its start again, as a
        // seam loop does); closed edges first where several go on from one vertex.
        let next_from = |at: NodeId, used: &[bool]| {
            (0..uses.len()).filter(|&j| !used[j] && uses[j].1 == at).min_by_key(|&j| uses[j].1 != uses[j].2)
        };
        while let Some(next) = next_from(tail, &used) {
            used[next] = true;
            chain.push_back(next);
            tail = uses[next].2;
        }
        // Back from the start, for a chain that does not close.
        if tail != head {
            let mut front = head;
            while let Some(prev) = (0..uses.len()).find(|&j| !used[j] && uses[j].2 == front) {
                used[prev] = true;
                chain.push_front(prev);
                front = uses[prev].1;
            }
        }
        chain.into_iter().map(|i| uses[i].0).collect()
    }

    /// The curve of `edge` in the parameter space of `plane`, obtained by projecting its 3D
    /// curve, with the edge's parameter range. `None` if the projection is not possible
    /// (e.g. the range lies outside a bounded curve's domain).
    #[doc(alias = "BRep_Tool::CurveOnPlane")]
    pub fn pcurve_on_plane(&self, edge: Shape, plane: &Plane) -> Option<(Curve2, f64, f64)> {
        use crate::geometry::ParametricCurve;
        let (curve, first, last) = self.edge_curve(edge)?;
        if !(first < last) {
            return None;
        }
        // Within a bounded, non-periodic curve's domain.
        let basis = match &curve {
            Curve::Trimmed(t) => t.basis().clone(),
            c => c.clone(),
        };
        if !basis.is_periodic() && (first < basis.first_parameter() || last > basis.last_parameter()) {
            return None;
        }
        let surface = crate::geometry::surface::Surface::Plane(*plane);
        let pcurve = crate::geometry::project::analytic_pcurve(&basis, &surface)?;
        Some((pcurve, first, last))
    }
}
