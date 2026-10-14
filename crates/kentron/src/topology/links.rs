//! The upward links of the model (which records use a node), and the edits that follow
//! them: removing nodes, redirecting uses, copying between models.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::math::Transform;

use super::entities::*;
use super::id::{Id, NodeId, NodeKind};
use super::shape::Shape;
use super::slots::{EdgeCurve, EdgeMeshPolyline, EdgePolyline, FaceMesh, FaceSurface, PCurve, PCurvePolyline};
use super::Model;

/// One use of a node by a parent: the parent, and the node's orientation and placement in
/// it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct ParentLink {
    pub parent: NodeId,
    pub orientation: Orientation,
    pub placement: Transform,
}

impl Model {
    /// True if `node` names a live record.
    pub(crate) fn is_live(&self, node: NodeId) -> bool {
        self.node_slot(node).is_some_and(|(_, _, removed)| !removed)
    }

    /// The live uses of `node` by live parents, in the parents' creation order (by stable
    /// id), the uses of one parent in their order.
    pub(super) fn parent_links(&self, node: NodeId) -> Vec<ParentLink> {
        let mut links: Vec<ParentLink> = Vec::new();
        let link = |parent: NodeId, orientation: Orientation| ParentLink { parent, orientation, placement: Transform::IDENTITY };
        let index = node.index;
        match node.kind {
            NodeKind::Vertex => {
                for (_, u) in self.vertex_uses.iter() {
                    if u.vertex.index() == index {
                        links.push(link(u.edge.into(), u.orientation));
                    }
                }
            }
            NodeKind::Edge => {
                for (id, c) in self.coedges.iter() {
                    if c.edge.index() == index {
                        links.push(link(id.into(), c.orientation));
                    }
                }
            }
            NodeKind::CoEdge => {
                if let Some(c) = self.get(CoEdgeId::from_index(index)) {
                    if self.get(c.wire).is_some() {
                        links.push(link(c.wire.into(), Orientation::Forward));
                    }
                }
            }
            NodeKind::Wire => {
                for (_, u) in self.wire_uses.iter() {
                    if u.wire.index() == index {
                        links.push(link(u.face.into(), u.orientation));
                    }
                }
            }
            NodeKind::Face => {
                for (_, u) in self.face_uses.iter() {
                    if u.face.index() == index {
                        links.push(link(u.shell.into(), u.orientation));
                    }
                }
            }
            NodeKind::Shell => {
                for (_, u) in self.shell_uses.iter() {
                    if u.shell.index() == index {
                        links.push(link(u.solid.into(), u.orientation));
                    }
                }
            }
            NodeKind::Solid => {
                for (_, u) in self.solid_uses.iter() {
                    if u.solid.index() == index {
                        links.push(link(u.comp_solid.into(), u.orientation));
                    }
                }
            }
            NodeKind::Instance => {
                for (_, u) in self.product_children.iter() {
                    if u.instance.index() == index {
                        links.push(ParentLink { parent: u.product.into(), orientation: Orientation::Forward, placement: u.placement });
                    }
                }
            }
            NodeKind::Compound | NodeKind::CompSolid | NodeKind::Product => {}
        }
        // Any node: compounds holding it, parts rooted at it, instances of it.
        for (_, c) in self.compound_children.iter() {
            if c.child == node {
                links.push(ParentLink { parent: c.compound.into(), orientation: c.orientation, placement: c.placement });
            }
        }
        for (id, p) in self.products.iter() {
            if let Some(shape) = p.shape.filter(|s| s.node == node) {
                links.push(ParentLink { parent: id.into(), orientation: shape.orientation, placement: shape.placement });
            }
        }
        for (id, i) in self.instances.iter() {
            if i.target == node {
                links.push(link(id.into(), Orientation::Forward));
            }
        }
        links.retain(|l| self.is_live(l.parent));
        // Creation order of the parents (a stable sort keeps each parent's uses in order).
        links.sort_by_key(|l| self.node_stable_id(l.parent));
        links
    }

    /// Soft-removes a record of a geometry or mesh slot, if it is live.
    fn drop_slot<T: super::Record>(&mut self, id: Id<T>) {
        if id.is_valid() {
            T::table_mut(self).remove(id);
        }
    }

    /// Takes a coedge out of its wire's list (the coedge itself stays).
    fn unlink_coedge(&mut self, coedge: CoEdgeId) {
        let Some(wire) = self.get(coedge).map(|c| c.wire) else { return };
        if self.get(wire).is_some() {
            self.wires[wire].coedges.retain(|&c| c != coedge);
        }
    }

    /// Removes a coedge record with its slots, out of its wire.
    fn drop_coedge(&mut self, coedge: CoEdgeId) {
        let Some(record) = self.get(coedge).cloned() else { return };
        self.unlink_coedge(coedge);
        self.drop_slot(record.pcurve);
        self.drop_slot(record.pcurve_polyline);
        self.drop_slot(record.mesh_polyline);
        self.coedges.remove(coedge);
    }

    /// Removes a use record and takes it out of its parent's list.
    pub(super) fn detach_use<U: super::UseRecord>(&mut self, id: Id<U>) {
        U::detach(self, id);
    }

    /// Soft-removes the record of `node` itself.
    fn remove_record(&mut self, node: NodeId) {
        macro_rules! remove {
            ($table:ident) => {{
                self.$table.remove(Id::from_index(node.index));
            }};
        }
        match node.kind {
            NodeKind::Vertex => remove!(vertices),
            NodeKind::Edge => remove!(edges),
            NodeKind::CoEdge => remove!(coedges),
            NodeKind::Wire => remove!(wires),
            NodeKind::Face => remove!(faces),
            NodeKind::Shell => remove!(shells),
            NodeKind::Solid => remove!(solids),
            NodeKind::Compound => remove!(compounds),
            NodeKind::CompSolid => remove!(comp_solids),
            NodeKind::Product => remove!(products),
            NodeKind::Instance => remove!(instances),
        }
    }

    pub(super) fn remove_node_impl(&mut self, node: NodeId) {
        if !self.is_live(node) {
            return;
        }
        let index = node.index;
        match node.kind {
            NodeKind::Vertex => {
                let uses: Vec<VertexUseId> = self.vertex_uses.iter().filter(|(_, u)| u.vertex.index() == index).map(|(id, _)| id).collect();
                for u in uses {
                    // The edge stays, open at that end.
                    self.detach_use(u);
                }
            }
            NodeKind::Edge => {
                let id = EdgeId::from_index(index);
                let record = self.edges[id].clone();
                // A curve shared with other edges stays theirs.
                let sharer = self.edges.iter().find(|&(e, r)| e != id && r.curve == record.curve && record.curve.is_valid()).map(|(e, _)| e);
                match sharer {
                    Some(other) => self.edge_curves.get_mut(record.curve).expect("the curve").edge = other,
                    None => self.drop_slot(record.curve),
                }
                self.drop_slot(record.polyline);
                for u in [record.start, record.end] {
                    if self.get(u).is_some() {
                        self.vertex_uses.remove(u);
                    }
                }
                let coedges: Vec<CoEdgeId> = self.coedges.iter().filter(|(_, c)| c.edge == id).map(|(c, _)| c).collect();
                for c in coedges {
                    self.drop_coedge(c);
                }
            }
            NodeKind::CoEdge => {
                let id = CoEdgeId::from_index(index);
                let record = self.coedges[id].clone();
                self.unlink_coedge(id);
                self.drop_slot(record.pcurve);
                self.drop_slot(record.pcurve_polyline);
                self.drop_slot(record.mesh_polyline);
            }
            NodeKind::Wire => {
                let id = WireId::from_index(index);
                let uses: Vec<WireUseId> = self.wire_uses.iter().filter(|(_, u)| u.wire == id).map(|(u, _)| u).collect();
                for u in uses {
                    self.detach_use(u);
                }
                // Its coedges are left free.
                for c in self.wires[id].coedges.clone() {
                    if let Some(coedge) = self.coedges.get_mut(c) {
                        coedge.wire = Id::INVALID;
                        coedge.face = Id::INVALID;
                    }
                }
            }
            NodeKind::Face => {
                let id = FaceId::from_index(index);
                let record = self.faces[id].clone();
                // A surface shared with other faces stays theirs.
                let sharer = self.faces.iter().find(|&(f, r)| f != id && r.surface == record.surface && record.surface.is_valid()).map(|(f, _)| f);
                match sharer {
                    Some(other) => self.face_surfaces.get_mut(record.surface).expect("the surface").face = other,
                    None => self.drop_slot(record.surface),
                }
                self.drop_slot(record.mesh);
                // Its wires are left free.
                for u in record.wires {
                    self.detach_use(u);
                }
                let uses: Vec<FaceUseId> = self.face_uses.iter().filter(|(_, u)| u.face == id).map(|(u, _)| u).collect();
                for u in uses {
                    self.detach_use(u);
                }
            }
            NodeKind::Shell => {
                let id = ShellId::from_index(index);
                for u in self.shells[id].faces.clone() {
                    self.face_uses.remove(u);
                }
                let uses: Vec<ShellUseId> = self.shell_uses.iter().filter(|(_, u)| u.shell == id).map(|(u, _)| u).collect();
                for u in uses {
                    self.detach_use(u);
                }
            }
            NodeKind::Solid => {
                let id = SolidId::from_index(index);
                for u in self.solids[id].shells.clone() {
                    self.shell_uses.remove(u);
                }
                let uses: Vec<SolidUseId> = self.solid_uses.iter().filter(|(_, u)| u.solid == id).map(|(u, _)| u).collect();
                for u in uses {
                    self.detach_use(u);
                }
            }
            NodeKind::CompSolid => {
                let id = CompSolidId::from_index(index);
                for u in self.comp_solids[id].solids.clone() {
                    self.solid_uses.remove(u);
                }
            }
            NodeKind::Compound => {
                let id = CompoundId::from_index(index);
                for u in self.compounds[id].children.clone() {
                    self.compound_children.remove(u);
                }
            }
            NodeKind::Product => {
                let id = ProductId::from_index(index);
                // Its own instances, and every instance of it.
                let own: Vec<InstanceId> = self.products[id].children.iter().filter_map(|&u| self.get(u)).map(|u| u.instance).collect();
                for u in self.products[id].children.clone() {
                    self.product_children.remove(u);
                }
                self.products[id].children.clear();
                for i in own {
                    self.remove_node_impl(i.into());
                }
            }
            NodeKind::Instance => {
                let id = InstanceId::from_index(index);
                let uses: Vec<ProductChildId> = self.product_children.iter().filter(|(_, u)| u.instance == id).map(|(u, _)| u).collect();
                for u in uses {
                    self.detach_use(u);
                }
            }
        }
        // Whatever kind: out of compounds, parts rooted at it emptied, instances of it gone.
        let uses: Vec<CompoundChildId> = self.compound_children.iter().filter(|(_, c)| c.child == node).map(|(u, _)| u).collect();
        for u in uses {
            self.detach_use(u);
        }
        let parts: Vec<ProductId> = self.products.iter().filter(|(_, p)| p.shape.is_some_and(|s| s.node == node)).map(|(p, _)| p).collect();
        for p in parts {
            self.products[p].shape = None;
        }
        self.remove_record(node);
        let instances: Vec<InstanceId> = self.instances.iter().filter(|(_, i)| i.target == node).map(|(i, _)| i).collect();
        for i in instances {
            self.remove_node_impl(i.into());
        }
    }

    /// The direct children of `node` as stored (whatever their placement).
    pub(super) fn child_nodes(&self, node: NodeId) -> Vec<NodeId> {
        self.children_local(Shape::new(node)).map(|s| s.node).collect()
    }

    pub(super) fn remove_subtree_impl(&mut self, node: NodeId) {
        if !self.is_live(node) {
            return;
        }
        // Everything below, top down (each after the nodes above it).
        let mut order: Vec<NodeId> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut queue = std::collections::VecDeque::from([node]);
        while let Some(n) = queue.pop_front() {
            if !seen.insert(n) {
                continue;
            }
            order.push(n);
            queue.extend(self.child_nodes(n));
        }
        self.remove_node_impl(node);
        // Below it, what nothing else uses any more; repeated until nothing more goes.
        loop {
            let mut changed = false;
            for &n in &order[1..] {
                if self.is_live(n) && self.parent_links(n).is_empty() {
                    self.remove_node_impl(n);
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }

    pub(super) fn replace_node_impl(&mut self, old: NodeId, new: NodeId) -> Result<()> {
        if old.kind != new.kind {
            return Err(Error::InvalidInput(format!("cannot replace {old:?} by a node of another kind, {new:?}")));
        }
        if !self.is_live(old) || !self.is_live(new) {
            return Err(Error::InvalidInput(format!("{old:?} or {new:?} is missing")));
        }
        if old == new {
            return Ok(());
        }
        let (o, n) = (old.index, new.index);
        macro_rules! redirect {
            ($table:ident, $field:ident) => {{
                let ids: Vec<_> = self.$table.iter().filter(|(_, r)| r.$field.index() == o).map(|(id, _)| id).collect();
                for id in ids {
                    self.$table[id].$field = Id::from_index(n);
                }
            }};
        }
        match old.kind {
            NodeKind::Vertex => redirect!(vertex_uses, vertex),
            NodeKind::Edge => redirect!(coedges, edge),
            NodeKind::CoEdge => {
                // The coedge's place in its wire.
                let wire = self.coedges[CoEdgeId::from_index(o)].wire;
                if self.get(wire).is_some() {
                    for c in self.wires[wire].coedges.iter_mut() {
                        if c.index() == o {
                            *c = CoEdgeId::from_index(n);
                        }
                    }
                    self.coedges[CoEdgeId::from_index(n)].wire = wire;
                    self.coedges[CoEdgeId::from_index(o)].wire = Id::INVALID;
                }
            }
            NodeKind::Wire => redirect!(wire_uses, wire),
            NodeKind::Face => redirect!(face_uses, face),
            NodeKind::Shell => redirect!(shell_uses, shell),
            NodeKind::Solid => redirect!(solid_uses, solid),
            NodeKind::Instance => redirect!(product_children, instance),
            NodeKind::Compound | NodeKind::CompSolid | NodeKind::Product => {}
        }
        let ids: Vec<CompoundChildId> = self.compound_children.iter().filter(|(_, c)| c.child == old).map(|(id, _)| id).collect();
        for id in ids {
            self.compound_children[id].child = new;
        }
        let ids: Vec<ProductId> = self.products.iter().filter(|(_, p)| p.shape.is_some_and(|s| s.node == old)).map(|(id, _)| id).collect();
        for id in ids {
            if let Some(shape) = &mut self.products[id].shape {
                shape.node = new;
            }
        }
        let ids: Vec<InstanceId> = self.instances.iter().filter(|(_, i)| i.target == old).map(|(id, _)| id).collect();
        for id in ids {
            self.instances[id].target = new;
        }
        self.remove_node_impl(old);
        Ok(())
    }

    /// Redirects every use of each key of `map` to its value (same kind; a value that is
    /// itself a key followed on), in one pass over the uses. The replaced nodes are left
    /// unused, not removed. Coedges cannot be replaced so.
    pub(super) fn replace_nodes_impl(&mut self, map: &HashMap<NodeId, NodeId>) -> Result<()> {
        let mut resolved: HashMap<NodeId, NodeId> = HashMap::default();
        for (&old, &new) in map {
            if old.kind != new.kind || old.kind == NodeKind::CoEdge {
                return Err(Error::InvalidInput(format!("cannot replace {old:?} by {new:?}")));
            }
            let mut target = new;
            let mut steps = 0;
            while let Some(&next) = map.get(&target) {
                target = next;
                steps += 1;
                if steps > map.len() {
                    return Err(Error::InvalidInput("replacements in a cycle".into()));
                }
            }
            if !self.is_live(target) {
                return Err(Error::InvalidInput(format!("{target:?} is missing")));
            }
            resolved.insert(old, target);
        }
        let to = |kind: NodeKind, index: u32| resolved.get(&NodeId { kind, index }).map(|n| n.index);
        macro_rules! redirect {
            ($table:ident, $field:ident, $kind:expr) => {{
                let ids: Vec<_> = self.$table.iter().filter_map(|(id, r)| to($kind, r.$field.index()).map(|n| (id, n))).collect();
                for (id, n) in ids {
                    self.$table[id].$field = Id::from_index(n);
                }
            }};
        }
        redirect!(vertex_uses, vertex, NodeKind::Vertex);
        redirect!(coedges, edge, NodeKind::Edge);
        redirect!(wire_uses, wire, NodeKind::Wire);
        redirect!(face_uses, face, NodeKind::Face);
        redirect!(shell_uses, shell, NodeKind::Shell);
        redirect!(solid_uses, solid, NodeKind::Solid);
        let ids: Vec<_> = self.compound_children.iter().filter_map(|(id, c)| resolved.get(&c.child).map(|&n| (id, n))).collect();
        for (id, n) in ids {
            self.compound_children[id].child = n;
        }
        let ids: Vec<_> = self.instances.iter().filter_map(|(id, i)| resolved.get(&i.target).map(|&n| (id, n))).collect();
        for (id, n) in ids {
            self.instances[id].target = n;
        }
        Ok(())
    }

    /// A copy of the node `node` of `source` (and everything below it) in this model, with
    /// `copied` mapping nodes already copied.
    pub(super) fn import_node(&mut self, source: &Model, node: NodeId, copied: &mut HashMap<NodeId, NodeId>) -> Result<NodeId> {
        if let Some(&n) = copied.get(&node) {
            return Ok(n);
        }
        let missing = || Error::InvalidInput(format!("no {node:?} in the source model"));
        let index = node.index;
        let new: NodeId = match node.kind {
            NodeKind::Vertex => {
                let v = source.get(VertexId::from_index(index)).ok_or_else(missing)?.clone();
                self.insert(v).into()
            }
            NodeKind::Edge => {
                let e = source.get(EdgeId::from_index(index)).ok_or_else(missing)?.clone();
                let edge = self.insert(Edge { curve: Id::INVALID, polyline: Id::INVALID, tolerance: e.tolerance, start: Id::INVALID, end: Id::INVALID });
                copied.insert(node, edge.into());
                let mut ends = [Id::INVALID; 2];
                for (slot, u) in ends.iter_mut().zip([e.start, e.end]) {
                    if let Some(u) = source.get(u) {
                        let vertex = VertexId::from_index(self.import_node(source, u.vertex.into(), copied)?.index);
                        *slot = self.insert(VertexUse { edge, vertex, orientation: u.orientation });
                    }
                }
                let curve = source.get(e.curve).map(|c| self.insert(EdgeCurve { edge, ..c.clone() })).unwrap_or(Id::INVALID);
                let polyline = source.get(e.polyline).map(|p| self.insert(EdgePolyline { edge, ..p.clone() })).unwrap_or(Id::INVALID);
                let record = &mut self.edges[edge];
                (record.start, record.end, record.curve, record.polyline) = (ends[0], ends[1], curve, polyline);
                edge.into()
            }
            NodeKind::CoEdge => {
                let c = source.get(CoEdgeId::from_index(index)).ok_or_else(missing)?.clone();
                let edge = EdgeId::from_index(self.import_node(source, c.edge.into(), copied)?.index);
                let coedge = self.insert(CoEdge {
                    wire: Id::INVALID,
                    edge,
                    face: Id::INVALID,
                    orientation: c.orientation,
                    pcurve: Id::INVALID,
                    pcurve_polyline: Id::INVALID,
                    mesh_polyline: Id::INVALID,
                });
                let pcurve = source.get(c.pcurve).map(|p| self.insert(PCurve { coedge, ..p.clone() })).unwrap_or(Id::INVALID);
                let polyline = source.get(c.pcurve_polyline).map(|p| self.insert(PCurvePolyline { coedge, ..p.clone() })).unwrap_or(Id::INVALID);
                let mesh = source.get(c.mesh_polyline).map(|p| self.insert(EdgeMeshPolyline { coedge, ..p.clone() })).unwrap_or(Id::INVALID);
                let record = &mut self.coedges[coedge];
                (record.pcurve, record.pcurve_polyline, record.mesh_polyline) = (pcurve, polyline, mesh);
                coedge.into()
            }
            NodeKind::Wire => {
                let w = source.get(WireId::from_index(index)).ok_or_else(missing)?.clone();
                let mut coedges = Vec::with_capacity(w.coedges.len());
                for c in w.coedges {
                    coedges.push(CoEdgeId::from_index(self.import_node(source, c.into(), copied)?.index));
                }
                let wire = self.insert(Wire { coedges: coedges.clone() });
                for c in coedges {
                    self.coedges[c].wire = wire;
                }
                wire.into()
            }
            NodeKind::Face => {
                let f = source.get(FaceId::from_index(index)).ok_or_else(missing)?.clone();
                let face = self.insert(Face { surface: Id::INVALID, mesh: Id::INVALID, tolerance: f.tolerance, wires: Vec::new() });
                copied.insert(node, face.into());
                let surface = source.get(f.surface).map(|s| self.insert(FaceSurface { face, ..s.clone() })).unwrap_or(Id::INVALID);
                let mesh = source.get(f.mesh).map(|m| self.insert(FaceMesh { face, ..m.clone() })).unwrap_or(Id::INVALID);
                let mut wires = Vec::new();
                for u in f.wires.iter().filter_map(|&u| source.get(u)) {
                    let wire = WireId::from_index(self.import_node(source, u.wire.into(), copied)?.index);
                    wires.push(self.insert(WireUse { face, wire, orientation: u.orientation }));
                    for c in self.wires[wire].coedges.clone() {
                        self.coedges[c].face = face;
                    }
                }
                let record = &mut self.faces[face];
                (record.surface, record.mesh, record.wires) = (surface, mesh, wires);
                face.into()
            }
            NodeKind::Shell => {
                let s = source.get(ShellId::from_index(index)).ok_or_else(missing)?.clone();
                let shell = self.insert(Shell::default());
                for u in s.faces.iter().filter_map(|&u| source.get(u)) {
                    let face = FaceId::from_index(self.import_node(source, u.face.into(), copied)?.index);
                    let used = self.insert(FaceUse { shell, face, orientation: u.orientation });
                    self.shells[shell].faces.push(used);
                }
                shell.into()
            }
            NodeKind::Solid => {
                let s = source.get(SolidId::from_index(index)).ok_or_else(missing)?.clone();
                let solid = self.insert(Solid::default());
                for u in s.shells.iter().filter_map(|&u| source.get(u)) {
                    let shell = ShellId::from_index(self.import_node(source, u.shell.into(), copied)?.index);
                    let used = self.insert(ShellUse { solid, shell, orientation: u.orientation });
                    self.solids[solid].shells.push(used);
                }
                solid.into()
            }
            NodeKind::CompSolid => {
                let s = source.get(CompSolidId::from_index(index)).ok_or_else(missing)?.clone();
                let comp_solid = self.insert(CompSolid::default());
                for u in s.solids.iter().filter_map(|&u| source.get(u)) {
                    let solid = SolidId::from_index(self.import_node(source, u.solid.into(), copied)?.index);
                    let used = self.insert(SolidUse { comp_solid, solid, orientation: u.orientation });
                    self.comp_solids[comp_solid].solids.push(used);
                }
                comp_solid.into()
            }
            NodeKind::Compound => {
                let c = source.get(CompoundId::from_index(index)).ok_or_else(missing)?.clone();
                let compound = self.insert(Compound::default());
                for u in c.children.iter().filter_map(|&u| source.get(u)).cloned().collect::<Vec<_>>() {
                    let child = self.import_node(source, u.child, copied)?;
                    let used = self.insert(CompoundChild { compound, child, ..u });
                    self.compounds[compound].children.push(used);
                }
                compound.into()
            }
            NodeKind::Product => {
                let p = source.get(ProductId::from_index(index)).ok_or_else(missing)?.clone();
                let product = self.insert(Product::default());
                copied.insert(node, product.into());
                let shape = match p.shape {
                    Some(s) => Some(Shape { node: self.import_node(source, s.node, copied)?, ..s }),
                    None => None,
                };
                self.products[product].shape = shape;
                for u in p.children.iter().filter_map(|&u| source.get(u)).cloned().collect::<Vec<_>>() {
                    let instance = InstanceId::from_index(self.import_node(source, u.instance.into(), copied)?.index);
                    let used = self.insert(ProductChild { product, instance, placement: u.placement });
                    self.products[product].children.push(used);
                }
                product.into()
            }
            NodeKind::Instance => {
                let i = source.get(InstanceId::from_index(index)).ok_or_else(missing)?.clone();
                let target = self.import_node(source, i.target, copied)?;
                self.insert(Instance { target }).into()
            }
        };
        copied.insert(node, new);
        Ok(new)
    }
}
