//! Low-level construction and editing of topology.
//!
//! These operations keep the model's records consistent (parents, children, uses) but do
//! not check geometric validity; modeling operations ([`crate::make`],
//! [`crate::primitives`], ...) build valid shapes on top of them.

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::math::Point;
use crate::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};

use super::entities::{
    CoEdge, CoEdgeId, CompSolid, CompSolidId, Compound, CompoundChildId, CompoundId, Edge, EdgeId, Face, FaceId, FaceUse, InstanceId,
    Orientation, Product, ProductId, Shell, ShellId, ShellUse, Solid, SolidId, SolidUse, Vertex, VertexId, VertexUse,
    VertexUseId, Wire, WireId, WireUse,
};
use super::slots::{EdgeCurve, EdgeMeshPolyline, EdgePolyline, FaceMesh, FaceSurface, PCurve, PCurvePolyline};
use crate::geometry::ParametricCurve;
use super::id::{Id, NodeId, NodeKind};
use super::shape::Shape;
use super::{Model, Record};
use crate::math::Transform;

impl Model {
    /// The live edge named by `shape`.
    fn edge_id(&self, shape: Shape) -> Result<EdgeId> {
        let id = EdgeId::try_from(shape.node)
            .map_err(|n| Error::InvalidInput(format!("{n:?} is not an edge")))?;
        self.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?;
        Ok(id)
    }

    /// A new node of `kind` with no children and no geometry: an empty wire, shell,
    /// compound, ..., or a vertex, edge or face whose geometry is yet to be set.
    #[doc(alias = "TopoDS_Builder::MakeCompound")]
    #[doc(alias = "BRep_Builder::MakeVertex")]
    pub fn new_empty(&mut self, kind: NodeKind) -> Shape {
        let node: NodeId = match kind {
            NodeKind::Vertex => self.insert(Vertex { point: Point::ORIGIN, tolerance: 0.0 }).into(),
            NodeKind::Edge => self
                .insert(Edge {
                    curve: Id::INVALID,
                    polyline: Id::INVALID,
                    tolerance: 0.0,
                    start: Id::INVALID,
                    end: Id::INVALID,
                })
                .into(),
            NodeKind::Wire => self.insert(Wire::default()).into(),
            NodeKind::Face => self
                .insert(Face { surface: Id::INVALID, mesh: Id::INVALID, tolerance: 0.0, wires: Vec::new() })
                .into(),
            NodeKind::Shell => self.insert(Shell::default()).into(),
            NodeKind::Solid => self.insert(Solid::default()).into(),
            NodeKind::Compound => self.insert(Compound::default()).into(),
            NodeKind::CompSolid => self.insert(CompSolid::default()).into(),
            NodeKind::Product => self.insert(Product::default()).into(),
            NodeKind::CoEdge | NodeKind::Instance => {
                todo!("new_empty({kind:?}): needs the node it refers to")
            }
        };
        Shape::new(node)
    }

    /// Appends `child` (with its orientation and placement) to `parent`. Fails if `parent`
    /// cannot contain that kind of node (a wire holds edges, a shell faces, a compound
    /// anything, ...).
    #[doc(alias = "TopoDS_Builder::Add")]
    pub fn add(&mut self, parent: Shape, child: Shape) -> Result<()> {
        fn invalid(msg: String) -> Error {
            Error::InvalidInput(msg)
        }
        if child.placement != Transform::IDENTITY && parent.kind() != NodeKind::Compound {
            todo!("add: placed children of a {:?}", parent.kind());
        }
        match (parent.kind(), child.kind()) {
            (NodeKind::Face, NodeKind::Wire) => {
                let face = FaceId::from_index(parent.node.index);
                let wire = WireId::from_index(child.node.index);
                let coedges = self.get(wire).ok_or_else(|| invalid(format!("no {wire:?}")))?.coedges.clone();
                self.get(face).ok_or_else(|| invalid(format!("no {face:?}")))?;
                let used = self.insert(WireUse { face, wire, orientation: child.orientation });
                self.faces[face].wires.push(used);
                for coedge in coedges {
                    self.coedges[coedge].face = face;
                }
            }
            (NodeKind::Shell, NodeKind::Face) => {
                let shell = ShellId::from_index(parent.node.index);
                let face = FaceId::from_index(child.node.index);
                self.get(shell).ok_or_else(|| invalid(format!("no {shell:?}")))?;
                self.get(face).ok_or_else(|| invalid(format!("no {face:?}")))?;
                let used = self.insert(FaceUse { shell, face, orientation: child.orientation });
                self.shells[shell].faces.push(used);
            }
            (NodeKind::Solid, NodeKind::Shell) => {
                let solid = SolidId::from_index(parent.node.index);
                let shell = ShellId::from_index(child.node.index);
                self.get(solid).ok_or_else(|| invalid(format!("no {solid:?}")))?;
                self.get(shell).ok_or_else(|| invalid(format!("no {shell:?}")))?;
                let used = self.insert(ShellUse { solid, shell, orientation: child.orientation });
                self.solids[solid].shells.push(used);
            }
            (NodeKind::CompSolid, NodeKind::Solid) => {
                let comp_solid = CompSolidId::from_index(parent.node.index);
                let solid = SolidId::from_index(child.node.index);
                self.get(comp_solid).ok_or_else(|| invalid(format!("no {comp_solid:?}")))?;
                self.get(solid).ok_or_else(|| invalid(format!("no {solid:?}")))?;
                let used = self.insert(SolidUse { comp_solid, solid, orientation: child.orientation });
                self.comp_solids[comp_solid].solids.push(used);
            }
            (NodeKind::Compound, _) => {
                let compound = CompoundId::from_index(parent.node.index);
                self.get(compound).ok_or_else(|| invalid(format!("no {compound:?}")))?;
                let used = self.insert(super::CompoundChild {
                    compound,
                    child: child.node,
                    orientation: child.orientation,
                    placement: child.placement,
                });
                self.compounds[compound].children.push(used);
            }
            (NodeKind::Wire, NodeKind::Edge) => {
                let wire = WireId::from_index(parent.node.index);
                let edge = EdgeId::from_index(child.node.index);
                self.get(wire).ok_or_else(|| invalid(format!("no {wire:?}")))?;
                let coedge = self.add_coedge(edge, child.orientation)?;
                let face = self.wire_face(wire);
                self.wires[wire].coedges.push(coedge);
                let record = &mut self.coedges[coedge];
                (record.wire, record.face) = (wire, face);
            }
            (NodeKind::Wire, NodeKind::CoEdge) => {
                let wire = WireId::from_index(parent.node.index);
                let coedge = CoEdgeId::from_index(child.node.index);
                self.get(wire).ok_or_else(|| invalid(format!("no {wire:?}")))?;
                let record = self.get(coedge).ok_or_else(|| invalid(format!("no {coedge:?}")))?;
                if self.get(record.wire).is_some() {
                    return Err(invalid(format!("{coedge:?} is already in a wire")));
                }
                let face = self.wire_face(wire);
                self.wires[wire].coedges.push(coedge);
                let record = &mut self.coedges[coedge];
                (record.wire, record.face) = (wire, face);
            }
            (NodeKind::Edge, NodeKind::Vertex) => {
                // A forward vertex starts the edge, a reversed one ends it.
                let edge = EdgeId::from_index(parent.node.index);
                let vertex = VertexId::from_index(child.node.index);
                self.get(edge).ok_or_else(|| invalid(format!("no {edge:?}")))?;
                self.get(vertex).ok_or_else(|| invalid(format!("no {vertex:?}")))?;
                let at_start = match child.orientation {
                    Orientation::Forward => true,
                    Orientation::Reversed => false,
                    o => todo!("add: a {o:?} vertex of an edge"),
                };
                let slot = if at_start { self.edges[edge].start } else { self.edges[edge].end };
                if self.get(slot).is_some() {
                    return Err(invalid(format!("{edge:?} already has that end")));
                }
                let used = self.insert(VertexUse { edge, vertex, orientation: child.orientation });
                if at_start {
                    self.edges[edge].start = used;
                } else {
                    self.edges[edge].end = used;
                }
            }
            (p, c) => todo!("add a {c:?} to a {p:?}"),
        }
        Ok(())
    }

    /// Fails unless `parent` and `child` are live and `child` can be added to `parent`.
    fn check_addable(&self, parent: Shape, child: Shape) -> Result<()> {
        let invalid = |msg: String| Err(Error::InvalidInput(msg));
        if !self.is_live(parent.node) {
            return invalid(format!("no {:?}", parent.node));
        }
        if !self.is_live(child.node) {
            return invalid(format!("no {:?}", child.node));
        }
        let fits = matches!(
            (parent.kind(), child.kind()),
            (NodeKind::Face, NodeKind::Wire)
                | (NodeKind::Shell, NodeKind::Face)
                | (NodeKind::Solid, NodeKind::Shell)
                | (NodeKind::CompSolid, NodeKind::Solid)
                | (NodeKind::Compound, _)
                | (NodeKind::Wire, NodeKind::Edge | NodeKind::CoEdge)
                | (NodeKind::Edge, NodeKind::Vertex)
        );
        match fits {
            true => Ok(()),
            false => invalid(format!("a {:?} cannot hold a {:?}", parent.kind(), child.kind())),
        }
    }

    /// The vertices a coedge runs from and to, following its orientation.
    fn coedge_vertices(&self, coedge: CoEdgeId) -> (VertexId, VertexId) {
        let record = &self.coedges[coedge];
        let (a, b) = self.edge_end_vertices(record.edge);
        if record.orientation == Orientation::Reversed { (b, a) } else { (a, b) }
    }

    /// The face a wire bounds (invalid if none).
    fn wire_face(&self, wire: WireId) -> FaceId {
        self.face_of_wire(wire).unwrap_or(Id::INVALID)
    }

    /// Fails unless `product` can be placed in `assembly` (both live, `assembly` not a part,
    /// no cycle, also counting the `pending` placements about to be made).
    fn check_instance(&self, assembly: ProductId, product: ProductId, pending: &[ProductId]) -> Result<()> {
        let invalid = |msg: String| Err(Error::InvalidInput(msg));
        let (Some(a), Some(_)) = (self.get(assembly), self.get(product)) else {
            return invalid(format!("no {assembly:?} or {product:?}"));
        };
        if a.is_part() {
            return invalid(format!("{assembly:?} is a part"));
        }
        // A cycle: the assembly is the product or lies below it.
        let mut stack = vec![product];
        let mut seen = std::collections::HashSet::new();
        while let Some(p) = stack.pop() {
            if p == assembly {
                return invalid(format!("placing {product:?} in {assembly:?} makes a cycle"));
            }
            if !seen.insert(p) {
                continue;
            }
            let Some(record) = self.get(p) else { continue };
            for u in record.children.iter().filter_map(|&u| self.get(u)) {
                if let Some(i) = self.get(u.instance) {
                    if let Ok(t) = ProductId::try_from(i.target) {
                        stack.push(t);
                    }
                }
            }
            if p == assembly {
                stack.extend(pending);
            }
        }
        Ok(())
    }

    fn place_instance(&mut self, assembly: ProductId, product: ProductId, placement: Transform) -> InstanceId {
        let instance = self.insert(super::Instance { target: product.into() });
        let used = self.insert(super::ProductChild { product: assembly, instance, placement });
        self.products[assembly].children.push(used);
        instance
    }

    /// Appends several children at once: either all are added or, if any cannot be (wrong
    /// kind, missing node), none is.
    pub fn add_all(&mut self, parent: Shape, children: &[Shape]) -> Result<()> {
        for &child in children {
            self.check_addable(parent, child)?;
        }
        for &child in children {
            self.add(parent, child)?;
        }
        Ok(())
    }

    /// Replaces the first use of `old` in `parent` by `new`, keeping its position. Fails if
    /// `parent` has no such child or cannot contain `new`.
    pub fn replace_child(&mut self, parent: Shape, old: Shape, new: Shape) -> Result<()> {
        self.check_addable(parent, new)?;
        let position = self
            .children_local(parent)
            .position(|c| c == Shape { placement: c.placement, ..old } && (parent.kind() != NodeKind::Compound || c.placement == old.placement))
            .ok_or_else(|| Error::InvalidInput(format!("{:?} is not a child of {:?}", old.node, parent.node)))?;
        let index = parent.node.index;
        let n = new.node.index;
        match parent.kind() {
            NodeKind::Face => {
                let u = self.faces[FaceId::from_index(index)].wires[position];
                let record = &mut self.wire_uses[u];
                (record.wire, record.orientation) = (WireId::from_index(n), new.orientation);
            }
            NodeKind::Shell => {
                let u = self.shells[ShellId::from_index(index)].faces[position];
                let record = &mut self.face_uses[u];
                (record.face, record.orientation) = (FaceId::from_index(n), new.orientation);
            }
            NodeKind::Solid => {
                let u = self.solids[SolidId::from_index(index)].shells[position];
                let record = &mut self.shell_uses[u];
                (record.shell, record.orientation) = (ShellId::from_index(n), new.orientation);
            }
            NodeKind::CompSolid => {
                let u = self.comp_solids[CompSolidId::from_index(index)].solids[position];
                let record = &mut self.solid_uses[u];
                (record.solid, record.orientation) = (SolidId::from_index(n), new.orientation);
            }
            NodeKind::Compound => {
                let u = self.compounds[CompoundId::from_index(index)].children[position];
                let record = &mut self.compound_children[u];
                (record.child, record.orientation, record.placement) = (new.node, new.orientation, new.placement);
            }
            NodeKind::Wire => {
                let wire = WireId::from_index(index);
                let (old_c, new_c) = (CoEdgeId::from_index(old.node.index), CoEdgeId::from_index(n));
                self.wires[wire].coedges[position] = new_c;
                self.coedges[old_c].wire = Id::INVALID;
                let face = self.coedges[old_c].face;
                let record = &mut self.coedges[new_c];
                (record.wire, record.face) = (wire, face);
            }
            NodeKind::Edge => {
                let edge = EdgeId::from_index(index);
                let u = if position == 0 { self.edges[edge].start } else { self.edges[edge].end };
                let record = &mut self.vertex_uses[u];
                (record.vertex, record.orientation) = (VertexId::from_index(n), new.orientation);
            }
            kind => return Err(Error::InvalidInput(format!("cannot replace the children of a {kind:?}"))),
        }
        Ok(())
    }

    /// Removes the use of `child` (same node, orientation and placement) from `parent`
    /// (its first such use, if there are several). Fails if `parent` has no such child.
    #[doc(alias = "TopoDS_Builder::Remove")]
    pub fn remove(&mut self, parent: Shape, child: Shape) -> Result<()> {
        let position = self
            .children_local(parent)
            .position(|c| c == child)
            .ok_or_else(|| Error::InvalidInput(format!("{:?} is not a child of {:?}", child.node, parent.node)))?;
        let index = parent.node.index;
        match parent.kind() {
            NodeKind::Face => self.detach_use(self.faces[FaceId::from_index(index)].wires[position]),
            NodeKind::Shell => self.detach_use(self.shells[ShellId::from_index(index)].faces[position]),
            NodeKind::Solid => self.detach_use(self.solids[SolidId::from_index(index)].shells[position]),
            NodeKind::CompSolid => self.detach_use(self.comp_solids[CompSolidId::from_index(index)].solids[position]),
            NodeKind::Compound => self.detach_use(self.compounds[CompoundId::from_index(index)].children[position]),
            NodeKind::Edge => {
                let edge = EdgeId::from_index(index);
                let u = if position == 0 { self.edges[edge].start } else { self.edges[edge].end };
                self.detach_use(u);
            }
            NodeKind::Wire => {
                let wire = WireId::from_index(index);
                let coedge = self.wires[wire].coedges.remove(position);
                let record = &mut self.coedges[coedge];
                (record.wire, record.face) = (Id::INVALID, Id::INVALID);
            }
            NodeKind::Product => {
                let product = ProductId::from_index(index);
                match (self.products[product].shape, position) {
                    (Some(_), 0) => self.products[product].shape = None,
                    (shape, i) => {
                        let u = self.products[product].children[i - usize::from(shape.is_some())];
                        self.detach_use(u);
                    }
                }
            }
            kind => return Err(Error::InvalidInput(format!("cannot remove the children of a {kind:?}"))),
        }
        Ok(())
    }

    /// A new vertex at `point`, not used by any edge.
    #[doc(alias = "BRepGraph::EditorView::Vertices::Add")]
    pub fn add_vertex(&mut self, point: Point, tolerance: f64) -> VertexId {
        self.insert(Vertex { point, tolerance })
    }

    /// A new edge from `start` to `end` without geometry (see [`Model::set_edge_curve`]).
    /// Fails, adding nothing, if a vertex is missing or removed.
    #[doc(alias = "BRepGraph::EditorView::Edges::Add")]
    pub fn add_edge(&mut self, start: VertexId, end: VertexId, tolerance: f64) -> Result<EdgeId> {
        for v in [start, end] {
            if self.get(v).is_none() {
                return Err(Error::InvalidInput(format!("no {v:?}")));
            }
        }
        let edge = self.insert(Edge {
            curve: Id::INVALID,
            polyline: Id::INVALID,
            tolerance,
            start: Id::INVALID,
            end: Id::INVALID,
        });
        // OCCT's convention: the first vertex is used forward, the last reversed.
        let start = self.insert(VertexUse { edge, vertex: start, orientation: Orientation::Forward });
        let end = self.insert(VertexUse { edge, vertex: end, orientation: Orientation::Reversed });
        let record = &mut self.edges[edge];
        record.start = start;
        record.end = end;
        Ok(edge)
    }

    /// A new coedge using `edge` with `orientation`, not yet in a wire. Fails, adding
    /// nothing, if the edge is missing or removed.
    #[doc(alias = "BRepGraph::EditorView::CoEdges::Add")]
    pub fn add_coedge(&mut self, edge: EdgeId, orientation: Orientation) -> Result<CoEdgeId> {
        if self.get(edge).is_none() {
            return Err(Error::InvalidInput(format!("no {edge:?}")));
        }
        Ok(self.insert(CoEdge {
            wire: Id::INVALID,
            edge,
            face: Id::INVALID,
            orientation,
            pcurve: Id::INVALID,
            pcurve_polyline: Id::INVALID,
            mesh_polyline: Id::INVALID,
        }))
    }

    /// A new wire made of `coedges` (which must not belong to another wire), in the given
    /// order, connected or not (see [`Model::set_coedge_order`]). Fails, adding nothing, if
    /// a coedge is missing, removed or already in a wire.
    #[doc(alias = "BRepGraph::EditorView::Wires::Add")]
    pub fn add_wire(&mut self, coedges: &[CoEdgeId]) -> Result<WireId> {
        for (i, &c) in coedges.iter().enumerate() {
            let record = self.get(c).ok_or_else(|| Error::InvalidInput(format!("no {c:?}")))?;
            if record.wire.is_valid() || coedges[..i].contains(&c) {
                return Err(Error::InvalidInput(format!("{c:?} is already in a wire")));
            }
        }
        let wire = self.insert(Wire { coedges: coedges.to_vec() });
        for &c in coedges {
            self.coedges[c].wire = wire;
        }
        Ok(wire)
    }

    /// Appends a free coedge to a wire. Fails if the wire already contains it or if it does
    /// not connect to the wire's end.
    pub fn append_coedge(&mut self, wire: WireId, coedge: CoEdgeId) -> Result<()> {
        let record = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?;
        let current = self.get(wire).ok_or_else(|| Error::InvalidInput(format!("no {wire:?}")))?;
        if record.wire.is_valid() && self.get(record.wire).is_some() {
            return Err(Error::InvalidInput(format!("{coedge:?} is already in a wire")));
        }
        if let Some(&last) = current.coedges.last() {
            if self.coedge_vertices(last).1 != self.coedge_vertices(coedge).0 {
                return Err(Error::InvalidInput(format!("{coedge:?} does not start where {wire:?} ends")));
            }
        }
        let face = self.wire_face(wire);
        self.wires[wire].coedges.push(coedge);
        let record = &mut self.coedges[coedge];
        (record.wire, record.face) = (wire, face);
        Ok(())
    }

    /// Removes a coedge from a wire (and the model), keeping the remaining coedges in
    /// connected order: removing one from a closed wire leaves an open chain starting after
    /// it.
    pub fn remove_coedge(&mut self, wire: WireId, coedge: CoEdgeId) -> Result<()> {
        let current = self.get(wire).ok_or_else(|| Error::InvalidInput(format!("no {wire:?}")))?;
        let position = current
            .coedges
            .iter()
            .position(|&c| c == coedge)
            .ok_or_else(|| Error::InvalidInput(format!("{coedge:?} is not in {wire:?}")))?;
        let closed = self.is_closed_wire(wire);
        let mut coedges = self.wires[wire].coedges.clone();
        coedges.remove(position);
        if closed {
            // The open chain starts after the removed coedge.
            let count = coedges.len().max(1);
            coedges.rotate_left(position % count);
        }
        self.wires[wire].coedges = coedges;
        self.coedges[coedge].wire = Id::INVALID;
        self.remove_node(coedge);
        Ok(())
    }

    /// Reorders a wire's coedges: `order` must be a permutation of them, and is normalized
    /// into connected order (an open chain starts at its free end). Coedges connect where
    /// their vertices are the same or within the vertices' tolerances. Runs that do not
    /// connect to each other are kept, each in connected order, in the order their
    /// coedges first appear in `order`. Returns whether `order` needed reordering or the
    /// wire's order changed.
    pub fn set_coedge_order(&mut self, wire: WireId, order: &[CoEdgeId]) -> Result<bool> {
        let current = self.get(wire).ok_or_else(|| Error::InvalidInput(format!("no {wire:?}")))?.coedges.clone();
        let (mut a, mut b) = (current.clone(), order.to_vec());
        a.sort();
        b.sort();
        if a != b {
            return Err(Error::InvalidInput(format!("not a permutation of the coedges of {wire:?}")));
        }
        // The vertices a coedge runs from and to.
        let ends = |c: CoEdgeId| -> (VertexId, VertexId) {
            let coedge = &self.coedges[c];
            let edge = &self.edges[coedge.edge];
            let (start, end) = (self.vertex_uses[edge.start].vertex, self.vertex_uses[edge.end].vertex);
            if coedge.orientation == Orientation::Reversed { (end, start) } else { (start, end) }
        };
        let meets = |u: VertexId, v: VertexId| {
            let (p, q) = (&self.vertices[u], &self.vertices[v]);
            u == v || p.point.distance(q.point) <= p.tolerance + q.tolerance
        };
        let mut left = order.to_vec();
        let mut sorted = Vec::with_capacity(order.len());
        while !left.is_empty() {
            // The run of the first coedge left, from its free end if it is an open chain.
            let mut first = left[0];
            let mut seen = vec![first];
            while let Some(&p) = left.iter().find(|&&p| p != first && meets(ends(p).1, ends(first).0)) {
                if seen.contains(&p) {
                    first = left[0];
                    break;
                }
                seen.push(p);
                first = p;
            }
            let mut c = first;
            left.retain(|&x| x != c);
            sorted.push(c);
            // Where several coedges continue the chain (around a seam), the one whose pcurve
            // starts where the current one's ends; first, one from the very vertex (where
            // vertices lie within one another's tolerances, as the ends of edges shorter than
            // those).
            let uv_gap = |a: CoEdgeId, b: CoEdgeId| -> f64 {
                let at = |c: CoEdgeId, start: bool| {
                    let record = &self.coedges[c];
                    let p = self.get(record.pcurve)?;
                    let forward = record.orientation != Orientation::Reversed;
                    let t = if start == forward { p.first } else { p.last };
                    Some(p.curve.point(t))
                };
                match (at(a, false), at(b, true)) {
                    (Some(p), Some(q)) => (p - q).length(),
                    _ => 0.0,
                }
            };
            loop {
                let next = left
                    .iter()
                    .enumerate()
                    .filter(|&(_, &n)| meets(ends(c).1, ends(n).0))
                    .min_by(|(_, &x), (_, &y)| (ends(c).1 != ends(x).0).cmp(&(ends(c).1 != ends(y).0)).then(uv_gap(c, x).total_cmp(&uv_gap(c, y))))
                    .map(|(i, _)| i);
                let Some(i) = next else { break };
                c = left.remove(i);
                sorted.push(c);
            }
        }
        let changed = sorted != current || sorted != order;
        self.wires[wire].coedges = sorted;
        Ok(changed)
    }

    /// Makes the wire's coedges that use `old` use `new` instead. Fails if `new` does not
    /// connect to the neighbouring coedges like `old` did (`either_direction` allows it to
    /// run backwards). Returns whether anything changed.
    pub fn replace_wire_edge(
        &mut self,
        wire: WireId,
        old: EdgeId,
        new: EdgeId,
        either_direction: bool,
    ) -> Result<bool> {
        let current = self.get(wire).ok_or_else(|| Error::InvalidInput(format!("no {wire:?}")))?.coedges.clone();
        let ends = |m: &Model, e: EdgeId| -> Option<(VertexId, VertexId)> {
            let r = m.get(e)?;
            Some((m.get(r.start)?.vertex, m.get(r.end)?.vertex))
        };
        if old == new {
            return Ok(false);
        }
        let (old_ends, new_ends) = (ends(self, old), ends(self, new));
        let new_ends = new_ends.ok_or_else(|| Error::InvalidInput(format!("no {new:?}")))?;
        let mut changed = false;
        for c in current {
            if self.coedges[c].edge != old {
                continue;
            }
            let flip = match old_ends {
                Some(o) if o == new_ends => false,
                Some((a, b)) if either_direction && (b, a) == new_ends => true,
                _ => return Err(Error::InvalidInput(format!("{new:?} does not connect like {old:?}"))),
            };
            let record = &mut self.coedges[c];
            record.edge = new;
            if flip {
                record.orientation = record.orientation.reversed();
            }
            changed = true;
        }
        Ok(changed)
    }

    /// Reverses a wire: its coedges in the opposite order, each with flipped orientation.
    pub fn reverse_wire(&mut self, wire: WireId) -> Result<()> {
        self.get(wire).ok_or_else(|| Error::InvalidInput(format!("no {wire:?}")))?;
        let mut coedges = self.wires[wire].coedges.clone();
        coedges.reverse();
        for &c in &coedges {
            let record = &mut self.coedges[c];
            record.orientation = record.orientation.reversed();
        }
        self.wires[wire].coedges = coedges;
        Ok(())
    }

    /// Reverses an edge's direction: swaps its start and end uses.
    pub fn reverse_edge(&mut self, edge: EdgeId) -> Result<()> {
        let record = self.get(edge).ok_or_else(|| Error::InvalidInput(format!("no {edge:?}")))?;
        let (start, end) = (record.start, record.end);
        // The vertices change places; each end keeps its orientation.
        // The end uses change places, each turned round (a start is used forward).
        for u in [start, end] {
            if let Some(used) = self.vertex_uses.get_mut(u) {
                used.orientation = used.orientation.reversed();
            }
        }
        let record = &mut self.edges[edge];
        (record.start, record.end) = (end, start);
        Ok(())
    }

    /// Soft-removes a node: its record and geometry slots are flagged as removed and every
    /// use of it is detached from its parents (its own children are left alone). Records
    /// that cannot exist without it go too: the coedges of a removed edge (and their places
    /// in wires); a removed face leaves its wires free; a removed vertex leaves its edges
    /// open at that end; a part whose shape is removed becomes empty; a removed product
    /// takes its own instances and every instance of it along. Removing an already removed
    /// node does nothing.
    #[doc(alias = "BRepGraph::EditorView::RemoveNode")]
    pub fn remove_node(&mut self, node: impl Into<NodeId>) {
        self.remove_node_impl(node.into());
    }

    /// Removes `node` and every node below it that nothing outside the removed set still
    /// uses (as [`Model::remove_node`] does for each). Shared sub-shapes, such as the edges
    /// a removed face shares with its neighbours, are kept.
    #[doc(alias = "BRepGraph::EditorView::GenOps::RemoveSubgraph")]
    pub fn remove_subtree(&mut self, node: impl Into<NodeId>) {
        self.remove_subtree_impl(node.into());
    }

    /// [`Model::replace_node`] for many nodes at once (each key of `map` by its value), in
    /// one pass over the uses: the replaced nodes are left unused rather than removed.
    pub fn replace_nodes(&mut self, map: &std::collections::HashMap<NodeId, NodeId>) -> Result<()> {
        self.replace_nodes_impl(map)
    }

    /// Redirects every use of `old` to `new` (same kind) and removes `old`. A structural
    /// edit: nothing is recorded in the history.
    #[doc(alias = "BRepGraph::EditorView::GenOps::ReplaceNode")]
    pub fn replace_node(&mut self, old: impl Into<NodeId>, new: impl Into<NodeId>) -> Result<()> {
        self.replace_node_impl(old.into(), new.into())
    }

    /// Splits `edge` at curve parameter `parameter` into two new edges meeting at `vertex`
    /// (start..vertex, vertex..end) and makes every wire that used `edge` use both halves,
    /// in its place; `edge` is removed. Returns the two halves.
    #[doc(alias = "BRepGraph::EditorView::Edges::Split")]
    pub fn split_edge(
        &mut self,
        edge: EdgeId,
        vertex: VertexId,
        parameter: f64,
    ) -> Result<(EdgeId, EdgeId)> {
        let record = self.get(edge).ok_or_else(|| Error::InvalidInput(format!("no {edge:?}")))?.clone();
        self.get(vertex).ok_or_else(|| Error::InvalidInput(format!("no {vertex:?}")))?;
        let curve = self.get(record.curve).cloned().ok_or_else(|| Error::InvalidInput(format!("{edge:?} has no 3D curve")))?;
        if !(curve.first < parameter && parameter < curve.last) {
            return Err(Error::InvalidInput(format!("{parameter} is not inside the range of {edge:?}")));
        }
        let (start, end) = self.edge_end_vertices(edge);
        let halves = [(start, vertex, curve.first, parameter), (vertex, end, parameter, curve.last)];
        let mut ids = Vec::with_capacity(2);
        for (a, b, first, last) in halves {
            let half = self.add_edge(a, b, record.tolerance)?;
            self.set_edge_curve(half.into(), curve.curve.clone(), record.tolerance)?;
            self.set_edge_range(half.into(), first, last)?;
            ids.push(half);
        }
        // Every wire using the edge uses the two halves in its place, in its direction.
        let coedges: Vec<CoEdgeId> = self.coedges.iter().filter(|(_, c)| c.edge == edge).map(|(c, _)| c).collect();
        for c in coedges {
            let old = self.coedges[c].clone();
            let pcurve = self.get(old.pcurve).cloned();
            let mut new = Vec::with_capacity(2);
            for (k, &half) in ids.iter().enumerate() {
                let coedge = self.add_coedge(half, old.orientation)?;
                if let Some(p) = &pcurve {
                    let (first, last) = if k == 0 { (p.first, parameter) } else { (parameter, p.last) };
                    self.set_pcurve_with_range(coedge, p.curve.clone(), first, last)?;
                }
                let r = &mut self.coedges[coedge];
                (r.wire, r.face) = (old.wire, old.face);
                new.push(coedge);
            }
            if old.orientation == Orientation::Reversed {
                new.reverse();
            }
            if let Some(wire) = self.wires.get_mut(old.wire) {
                if let Some(i) = wire.coedges.iter().position(|&x| x == c) {
                    wire.coedges.splice(i..=i, new);
                }
            }
        }
        self.remove_node(edge);
        Ok((ids[0], ids[1]))
    }

    /// Soft-removes a use entry and detaches it from its parent's child list.
    #[doc(alias = "BRepGraph::EditorView::RemoveRef")]
    pub fn remove_use<U: UseRecord>(&mut self, id: Id<U>) {
        if U::table(self).get(id).is_some() {
            self.detach_use(id);
        }
    }

    /// Removes several use entries: all of them or, if any is missing, none.
    pub fn remove_uses<U: UseRecord>(&mut self, ids: &[Id<U>]) -> Result<()> {
        if let Some(missing) = ids.iter().find(|&&id| U::table(self).get(id).is_none()) {
            return Err(Error::InvalidInput(format!("no {missing:?}")));
        }
        for &id in ids {
            self.detach_use(id);
        }
        Ok(())
    }

    /// A new part: a product whose content is the topology rooted at `shape`. Fails if
    /// `shape` is not a live topology node (a product or instance cannot be a part's shape).
    pub fn new_part(&mut self, shape: impl Into<Shape>) -> Result<ProductId> {
        let shape = shape.into();
        if matches!(shape.kind(), NodeKind::Product | NodeKind::Instance) || !self.is_live(shape.node) {
            return Err(Error::InvalidInput(format!("{:?} cannot be the shape of a part", shape.node)));
        }
        Ok(self.insert(Product { shape: Some(shape), children: Vec::new() }))
    }

    /// A deep copy of `shape` from another model (or another part of this one, through a
    /// clone): its nodes, uses, geometry and meshes are duplicated into this model. Shared
    /// subshapes stay shared within the copy.
    pub fn import(&mut self, source: &Model, shape: Shape) -> Result<Shape> {
        let mut copied = std::collections::HashMap::new();
        let node = self.import_node(source, shape.node, &mut copied)?;
        Ok(Shape { node, ..shape })
    }

    /// The assembly an instance is placed in, and its placement there.
    pub fn instance_placement(&self, instance: InstanceId) -> Option<(ProductId, Transform)> {
        self.get(instance)?;
        self.product_children.iter().find(|(_, u)| u.instance == instance).map(|(_, u)| (u.product, u.placement))
    }

    /// The instances of `product`, in any assembly.
    pub fn instances_of(&self, product: ProductId) -> Vec<InstanceId> {
        self.instances.iter().filter(|(_, i)| i.target == product).map(|(id, _)| id).collect()
    }

    /// Moves a compound's child: its placement becomes `t` applied after the current one.
    #[doc(alias = "BRepGraph_Transform::MoveRef")]
    pub fn move_compound_child(&mut self, child: CompoundChildId, t: &Transform) -> Result<()> {
        let record = self.compound_children.get_mut(child).ok_or_else(|| Error::InvalidInput(format!("no {child:?}")))?;
        record.placement = *t * record.placement;
        Ok(())
    }

    /// Moves an instance within its assembly: its placement becomes `t` applied after the
    /// current one. Fails (changing nothing) if `t` scales: instances are placed rigidly.
    pub fn move_instance(&mut self, instance: InstanceId, t: &Transform) -> Result<()> {
        if (t.scale_factor() - 1.0).abs() > crate::tolerance::PARAMETRIC || t.is_negative() {
            return Err(Error::InvalidInput("instances are placed rigidly: no scaling".into()));
        }
        let (u, _) = self
            .product_children
            .iter()
            .find(|(_, u)| u.instance == instance)
            .ok_or_else(|| Error::InvalidInput(format!("{instance:?} is not placed in an assembly")))?;
        let record = &mut self.product_children[u];
        record.placement = *t * record.placement;
        Ok(())
    }

    /// A new, empty assembly product.
    pub fn new_assembly(&mut self) -> ProductId {
        self.insert(Product::default())
    }

    /// Places an instance of `product` in `assembly`. Fails if either product is missing
    /// or removed, if `assembly` is a part, or if the instance would make the product graph
    /// cyclic (including `product == assembly`).
    pub fn add_instance(
        &mut self,
        assembly: ProductId,
        product: ProductId,
        placement: Transform,
    ) -> Result<InstanceId> {
        self.check_instance(assembly, product, &[])?;
        Ok(self.place_instance(assembly, product, placement))
    }

    /// Places several instances at once: all of them or none.
    pub fn add_instances(
        &mut self,
        assembly: ProductId,
        children: &[(ProductId, Transform)],
    ) -> Result<Vec<InstanceId>> {
        let mut pending = Vec::new();
        for &(product, _) in children {
            self.check_instance(assembly, product, &pending)?;
            pending.push(product);
        }
        Ok(children.iter().map(|&(product, placement)| self.place_instance(assembly, product, placement)).collect())
    }

    /// Removes an instance (and its use in its assembly).
    pub fn remove_instance(&mut self, instance: InstanceId) -> Result<()> {
        self.get(instance).ok_or_else(|| Error::InvalidInput(format!("no {instance:?}")))?;
        self.remove_node(instance);
        Ok(())
    }

    /// Products not instanced by any other product, in creation order.
    #[doc(alias = "BRepGraph::RootProductIds")]
    pub fn root_products(&self) -> Vec<ProductId> {
        let placed: std::collections::HashSet<u32> = self
            .product_children
            .iter()
            .filter_map(|(_, u)| self.get(u.instance))
            .filter(|i| i.target.kind == NodeKind::Product)
            .map(|i| i.target.index)
            .collect();
        self.products.iter().map(|(id, _)| id).filter(|id| !placed.contains(&id.index())).collect()
    }

    /// A new node of the same kind as `shape`, with the same own geometry (point, curve,
    /// surface, tolerance) but no children.
    #[doc(alias = "TopoDS_TShape::EmptyCopy")]
    pub fn empty_copy(&mut self, shape: Shape) -> Shape {
        let index = shape.node.index;
        let node: NodeId = match shape.kind() {
            NodeKind::Vertex => {
                let v = self.vertices[VertexId::from_index(index)].clone();
                self.insert(v).into()
            }
            NodeKind::Edge => {
                let e = self.edges[EdgeId::from_index(index)].clone();
                let edge = self.insert(Edge { curve: Id::INVALID, polyline: Id::INVALID, tolerance: e.tolerance, start: Id::INVALID, end: Id::INVALID });
                if let Some(c) = self.get(e.curve).cloned() {
                    let slot = self.insert(EdgeCurve { edge, ..c });
                    self.edges[edge].curve = slot;
                }
                edge.into()
            }
            NodeKind::Face => {
                let f = self.faces[FaceId::from_index(index)].clone();
                let face = self.insert(Face { surface: Id::INVALID, mesh: Id::INVALID, tolerance: f.tolerance, wires: Vec::new() });
                if let Some(s) = self.get(f.surface).cloned() {
                    let slot = self.insert(FaceSurface { face, ..s });
                    self.faces[face].surface = slot;
                }
                face.into()
            }
            NodeKind::CoEdge => {
                let c = self.coedges[CoEdgeId::from_index(index)].clone();
                self.add_coedge(c.edge, c.orientation).expect("the coedge's edge").into()
            }
            NodeKind::Instance => {
                let target = self.instances[InstanceId::from_index(index)].target;
                self.insert(super::Instance { target }).into()
            }
            kind => self.new_empty(kind).node,
        };
        Shape { node, ..shape }
    }

    /// Sets the point and tolerance of a vertex.
    #[doc(alias = "BRep_Builder::UpdateVertex")]
    pub fn set_vertex_point(&mut self, vertex: Shape, point: Point, tolerance: f64) -> Result<()> {
        let id = VertexId::try_from(vertex.node)
            .map_err(|n| Error::InvalidInput(format!("{n:?} is not a vertex")))?;
        self.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?;
        if vertex.placement != crate::math::Transform::IDENTITY {
            todo!("setting the point of a placed vertex");
        }
        let record = &mut self.vertices[id];
        (record.point, record.tolerance) = (point, tolerance);
        Ok(())
    }

    /// Sets the tolerance of a vertex.
    pub fn set_vertex_tolerance(&mut self, vertex: VertexId, tolerance: f64) -> Result<()> {
        self.get(vertex).ok_or_else(|| Error::InvalidInput(format!("no {vertex:?}")))?;
        self.vertices[vertex].tolerance = tolerance;
        Ok(())
    }

    /// Sets the tolerance of an edge.
    #[doc(alias = "BRepGraph::EditorView::Edges::SetTolerance")]
    pub fn set_edge_tolerance(&mut self, edge: EdgeId, tolerance: f64) -> Result<()> {
        self.get(edge).ok_or_else(|| Error::InvalidInput(format!("no {edge:?}")))?;
        self.edges[edge].tolerance = tolerance;
        Ok(())
    }

    /// Sets the tolerance of a face.
    #[doc(alias = "BRepGraph::EditorView::Faces::SetTolerance")]
    pub fn set_face_tolerance(&mut self, face: FaceId, tolerance: f64) -> Result<()> {
        self.get(face).ok_or_else(|| Error::InvalidInput(format!("no {face:?}")))?;
        self.faces[face].tolerance = tolerance;
        Ok(())
    }

    /// Sets the 3D curve and tolerance of an edge (its range becomes the curve's domain).
    #[doc(alias = "BRep_Builder::UpdateEdge")]
    pub fn set_edge_curve(&mut self, edge: Shape, curve: Curve, tolerance: f64) -> Result<()> {
        let id = self.edge_id(edge)?;
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let slot = self.edges[id].curve;
        let record = EdgeCurve { edge: id, curve, first, last };
        if slot.is_valid() && (slot.index() as usize) < self.edge_curves.capacity_len() {
            // Reuse the edge's slot, reviving it if it was cleared.
            self.edge_curves.replace(slot, record);
        } else {
            let slot = self.insert(record);
            self.edges[id].curve = slot;
        }
        self.edges[id].tolerance = tolerance;
        Ok(())
    }

    /// Sets or clears the polyline of an edge.
    pub fn set_edge_polyline(&mut self, edge: Shape, polyline: Option<Polyline>) -> Result<()> {
        let id = self.edge_id(edge)?;
        let slot = self.edges[id].polyline;
        let live = slot.is_valid() && (slot.index() as usize) < self.edge_polylines.capacity_len();
        match polyline {
            Some(polyline) if live => self.edge_polylines.replace(slot, EdgePolyline { edge: id, polyline }),
            Some(polyline) => {
                let slot = self.insert(EdgePolyline { edge: id, polyline });
                self.edges.get_mut_keeping_content(id).expect("the edge").polyline = slot;
            }
            None => {
                if live {
                    self.edge_polylines.remove(slot);
                }
                self.edges.get_mut_keeping_content(id);
            }
        }
        Ok(())
    }

    /// Sets the surface and tolerance of a face.
    #[doc(alias = "BRep_Builder::UpdateFace")]
    pub fn set_face_surface(&mut self, face: Shape, surface: Surface, tolerance: f64) -> Result<()> {
        let id = FaceId::try_from(face.node)
            .map_err(|n| Error::InvalidInput(format!("{n:?} is not a face")))?;
        let slot = self.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?.surface;
        if self.get(slot).is_some() {
            self.face_surfaces[slot].surface = surface;
        } else if slot.is_valid() && (slot.index() as usize) < self.face_surfaces.capacity_len() {
            // Reuse the face's slot, reviving it if it was cleared.
            self.face_surfaces.replace(slot, FaceSurface { face: id, surface });
        } else {
            let slot = self.insert(FaceSurface { face: id, surface });
            self.faces[id].surface = slot;
        }
        self.faces[id].tolerance = tolerance;
        Ok(())
    }

    /// Sets or clears the mesh of a face. A face may have a mesh and no surface.
    pub fn set_face_mesh(&mut self, face: Shape, mesh: Option<TriangleMesh>) -> Result<()> {
        let id = FaceId::try_from(face.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a face")))?;
        let slot = self.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?.mesh;
        let live = slot.is_valid() && (slot.index() as usize) < self.face_meshes.capacity_len();
        match mesh {
            Some(mesh) if live => self.face_meshes.replace(slot, FaceMesh { face: id, mesh }),
            Some(mesh) => {
                let slot = self.insert(FaceMesh { face: id, mesh });
                self.faces.get_mut_keeping_content(id).expect("the face").mesh = slot;
            }
            None => {
                if live {
                    self.face_meshes.remove(slot);
                }
                self.faces.get_mut_keeping_content(id);
            }
        }
        Ok(())
    }

    /// Makes the end of `edge` given by `old` (its start or end use) refer to `vertex`,
    /// keeping its orientation. Returns the new use (`old` itself if `vertex` is already
    /// there). Fails if the edge is removed or `old` is not one of its ends. Only this edge
    /// changes: edges sharing the old vertex keep it.
    #[doc(alias = "BRepGraph::EditorView::EdgeOps::ReplaceVertex")]
    pub fn replace_edge_vertex(
        &mut self,
        edge: EdgeId,
        old: VertexUseId,
        vertex: VertexId,
    ) -> Result<VertexUseId> {
        let record = self.get(edge).ok_or_else(|| Error::InvalidInput(format!("no {edge:?}")))?;
        let at_start = record.start == old;
        if !at_start && record.end != old {
            return Err(Error::InvalidInput(format!("{old:?} is not an end of {edge:?}")));
        }
        if self.get(vertex).is_none() {
            return Err(Error::InvalidInput(format!("no {vertex:?}")));
        }
        let previous = self.get(old).ok_or_else(|| Error::InvalidInput(format!("no {old:?}")))?.clone();
        if previous.vertex == vertex {
            return Ok(old);
        }
        let new = self.insert(VertexUse { edge, vertex, orientation: previous.orientation });
        self.vertex_uses.remove(old);
        let record = &mut self.edges[edge];
        if at_start {
            record.start = new;
        } else {
            record.end = new;
        }
        Ok(new)
    }

    /// Sets a coedge's orientation relative to its edge.
    pub fn set_coedge_orientation(&mut self, coedge: CoEdgeId, orientation: Orientation) -> Result<()> {
        let record = self.coedges.get_mut(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?;
        record.orientation = orientation;
        Ok(())
    }

    /// Removes the 3D curve of an edge. The edge keeps its slot, which a later
    /// [`Model::set_edge_curve`] reuses.
    pub fn clear_edge_curve(&mut self, edge: Shape) -> Result<()> {
        let id = self.edge_id(edge)?;
        let slot = self.edges[id].curve;
        self.edge_curves.remove(slot);
        self.edges.get_mut(id);
        Ok(())
    }

    /// Removes the surface of a face, keeping its slot (see [`Model::clear_edge_curve`]).
    pub fn clear_face_surface(&mut self, face: Shape) -> Result<()> {
        let id = FaceId::try_from(face.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a face")))?;
        let slot = self.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?.surface;
        if slot.is_valid() {
            self.face_surfaces.remove(slot);
        }
        self.faces.get_mut(id);
        Ok(())
    }

    /// Sets the pcurve of a coedge, keeping its current parameter range (the curve's domain
    /// if it had none).
    pub fn set_pcurve(&mut self, coedge: CoEdgeId, curve: Curve2) -> Result<()> {
        let record = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?;
        let (first, last) = match self.get(record.pcurve) {
            Some(old) => (old.first, old.last),
            None => (curve.first_parameter(), curve.last_parameter()),
        };
        self.set_pcurve_with_range(coedge, curve, first, last)
    }

    /// Sets the pcurve of a coedge restricted to `[first, last]`.
    pub fn set_pcurve_with_range(
        &mut self,
        coedge: CoEdgeId,
        curve: Curve2,
        first: f64,
        last: f64,
    ) -> Result<()> {
        let slot = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?.pcurve;
        let record = PCurve { coedge, curve, first, last };
        if slot.is_valid() && (slot.index() as usize) < self.pcurves.capacity_len() {
            self.pcurves.replace(slot, record);
        } else {
            let slot = self.insert(record);
            self.coedges[coedge].pcurve = slot;
        }
        Ok(())
    }

    /// Removes the pcurve of a coedge, keeping its slot (see [`Model::clear_edge_curve`]).
    pub fn clear_pcurve(&mut self, coedge: CoEdgeId) -> Result<()> {
        let slot = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?.pcurve;
        if slot.is_valid() {
            self.pcurves.remove(slot);
        }
        self.coedges.get_mut(coedge);
        Ok(())
    }

    /// Sets the parameter range of a coedge on its pcurve.
    pub fn set_coedge_range(&mut self, coedge: CoEdgeId, first: f64, last: f64) -> Result<()> {
        let slot = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?.pcurve;
        let pcurve = self.pcurves.get_mut(slot).ok_or_else(|| Error::InvalidInput(format!("{coedge:?} has no pcurve")))?;
        (pcurve.first, pcurve.last) = (first, last);
        Ok(())
    }

    /// Sets or clears a coedge's 2D polyline in its face's parameter space.
    pub fn set_pcurve_polyline(&mut self, coedge: CoEdgeId, polyline: Option<Polyline2>) -> Result<()> {
        let slot = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?.pcurve_polyline;
        let live = slot.is_valid() && (slot.index() as usize) < self.pcurve_polylines.capacity_len();
        match polyline {
            Some(polyline) if live => self.pcurve_polylines.replace(slot, PCurvePolyline { coedge, polyline }),
            Some(polyline) => {
                let slot = self.insert(PCurvePolyline { coedge, polyline });
                self.coedges[coedge].pcurve_polyline = slot;
            }
            None if live => {
                self.pcurve_polylines.remove(slot);
            }
            None => {}
        }
        Ok(())
    }

    /// Sets or clears a coedge's discretization as nodes of its face's mesh.
    pub fn set_coedge_mesh_polyline(
        &mut self,
        coedge: CoEdgeId,
        polyline: Option<MeshPolyline>,
    ) -> Result<()> {
        let slot = self.get(coedge).ok_or_else(|| Error::InvalidInput(format!("no {coedge:?}")))?.mesh_polyline;
        let live = slot.is_valid() && (slot.index() as usize) < self.edge_mesh_polylines.capacity_len();
        match polyline {
            Some(polyline) if live => self.edge_mesh_polylines.replace(slot, EdgeMeshPolyline { coedge, polyline }),
            Some(polyline) => {
                let slot = self.insert(EdgeMeshPolyline { coedge, polyline });
                self.coedges.get_mut_keeping_content(coedge).expect("the coedge").mesh_polyline = slot;
            }
            None if live => {
                self.edge_mesh_polylines.remove(slot);
            }
            None => {}
        }
        Ok(())
    }

    /// Removes the meshes and polylines stored on `shape` and everything below it.
    #[doc(alias = "BRepTools::Clean")]
    pub fn clear_meshes(&mut self, shape: Shape) {
        let faces: Vec<Shape> = self.unique_subshapes(shape, NodeKind::Face);
        let edges: Vec<Shape> = self.unique_subshapes(shape, NodeKind::Edge);
        for face in faces {
            let _ = self.set_face_mesh(face, None);
            let Ok(id) = FaceId::try_from(face.node) else { continue };
            for coedge in self.face_coedges(id) {
                let _ = self.set_coedge_mesh_polyline(coedge, None);
                let _ = self.set_pcurve_polyline(coedge, None);
            }
        }
        for edge in edges {
            let _ = self.set_edge_polyline(edge, None);
        }
    }

    /// Sets the parameter range of an edge on its 3D curve.
    #[doc(alias = "BRep_Builder::Range")]
    pub fn set_edge_range(&mut self, edge: Shape, first: f64, last: f64) -> Result<()> {
        let id = self.edge_id(edge)?;
        let slot = self.edges[id].curve;
        let curve = self
            .edge_curves
            .get_mut(slot)
            .ok_or_else(|| Error::InvalidInput(format!("{id:?} has no 3D curve")))?;
        curve.first = first;
        curve.last = last;
        Ok(())
    }
}

/// A record placing a child in a parent (see [`super::entities`]).
pub trait UseRecord: Record {
    /// Removes the use and takes it out of its parent's list.
    #[doc(hidden)]
    fn detach(model: &mut Model, id: Id<Self>);
}

impl UseRecord for super::VertexUse {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(edge) = model.get(id).map(|u| u.edge) else { return };
        if let Some(record) = model.edges.get_mut(edge) {
            if record.start == id {
                record.start = Id::INVALID;
            }
            if record.end == id {
                record.end = Id::INVALID;
            }
        }
        model.vertex_uses.remove(id);
    }
}

impl UseRecord for super::WireUse {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(u) = model.get(id).cloned() else { return };
        if let Some(face) = model.faces.get_mut(u.face) {
            face.wires.retain(|&w| w != id);
        }
        // The wire's coedges no longer bound the face.
        if let Some(wire) = model.get(u.wire).cloned() {
            for c in wire.coedges {
                if let Some(coedge) = model.coedges.get_mut(c) {
                    if coedge.face == u.face {
                        coedge.face = Id::INVALID;
                    }
                }
            }
        }
        model.wire_uses.remove(id);
    }
}

impl UseRecord for super::FaceUse {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(shell) = model.get(id).map(|u| u.shell) else { return };
        if let Some(record) = model.shells.get_mut(shell) {
            record.faces.retain(|&f| f != id);
        }
        model.face_uses.remove(id);
    }
}

impl UseRecord for super::ShellUse {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(solid) = model.get(id).map(|u| u.solid) else { return };
        if let Some(record) = model.solids.get_mut(solid) {
            record.shells.retain(|&s| s != id);
        }
        model.shell_uses.remove(id);
    }
}

impl UseRecord for super::SolidUse {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(comp_solid) = model.get(id).map(|u| u.comp_solid) else { return };
        if let Some(record) = model.comp_solids.get_mut(comp_solid) {
            record.solids.retain(|&s| s != id);
        }
        model.solid_uses.remove(id);
    }
}

impl UseRecord for super::CompoundChild {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(compound) = model.get(id).map(|u| u.compound) else { return };
        if let Some(record) = model.compounds.get_mut(compound) {
            record.children.retain(|&c| c != id);
        }
        model.compound_children.remove(id);
    }
}

impl UseRecord for super::ProductChild {
    fn detach(model: &mut Model, id: Id<Self>) {
        let Some(product) = model.get(id).map(|u| u.product) else { return };
        if let Some(record) = model.products.get_mut(product) {
            record.children.retain(|&c| c != id);
        }
        model.product_children.remove(id);
    }
}
