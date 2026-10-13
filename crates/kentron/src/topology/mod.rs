//! Boundary representation: the [`Model`] arena owning all topology, geometry and meshes
//! of one logical assembly.
//!
//! # Memory model
//!
//! There is no reference counting. A [`Model`] owns one [`Table`] per kind of record:
//!
//! * topology [`entities`] (vertex, edge, coedge, wire, face, shell, solid, compound,
//!   product, instance);
//! * *use* entries placing a child in a parent with an [`Orientation`] (and, for
//!   compounds and products, a placement), so shared children are stored once;
//! * geometry and mesh [`slots`], each owned by one entity and holding its geometry by value.
//!
//! Records refer to each other by typed `u32` ids ([`Id<T>`]). A [`Shape`] is a `Copy`
//! value naming a node, an orientation and a placement.
//!
//! Removal is soft (flagged in place); [`Model::compact`] reclaims space and renumbers ids.
//! [`StableId`]s survive compaction, for history and external references. Recording a
//! modification history is opt-in ([`Model::enable_history`]).
//!
//! Mutation needs `&mut Model`, so the borrow checker rules out concurrent writers while
//! reads of `&Model` can run in parallel. Geometry not attached to a model is an ordinary
//! owned value.
//!
//! Trade-offs: ids are only meaningful in the model that issued them; operands of a binary
//! operation must live in the same model (or be imported, which copies); parallel
//! algorithms read `&Model` and apply their edits afterwards. Equal geometry used by
//! several entities is duplicated until [`Model::deduplicate`] merges it.
//!
//! The layout follows OCCT's `BRepGraph`.

mod collect;
mod edit;
pub use collect::CollectOptions;
pub use edit::UseRecord;
pub mod entities;
mod explore;
mod history;
mod dedup;
mod id;
mod links;
mod navigate;
mod query;
mod derived;
mod revision;
mod shape;
pub mod slots;
mod substitution;
mod table;
mod wire_curve;

pub use entities::*;
pub use explore::{AncestorExplorer, Explorer};
pub use history::{ChangeKind, History, HistoryEvent};
pub use query::Children;
pub use id::{Id, Node, NodeId, NodeKind, StableId};
pub use navigate::Relation;
pub use revision::{ModelUid, Revision, Stamp};
pub use shape::Shape;
pub(crate) use shape::ByOccurrence;
pub use slots::*;
pub use substitution::Substitution;
pub use table::Table;
pub use wire_curve::WireCurve;

/// Storage for every record of one assembly. See the [module docs](self).
#[doc(alias = "BRepGraph")]
#[derive(Debug, Clone, Default)]
pub struct Model {
    vertices: Table<Vertex>,
    edges: Table<Edge>,
    coedges: Table<CoEdge>,
    wires: Table<Wire>,
    faces: Table<Face>,
    shells: Table<Shell>,
    solids: Table<Solid>,
    compounds: Table<Compound>,
    comp_solids: Table<CompSolid>,
    products: Table<Product>,
    instances: Table<Instance>,

    vertex_uses: Table<VertexUse>,
    wire_uses: Table<WireUse>,
    face_uses: Table<FaceUse>,
    shell_uses: Table<ShellUse>,
    solid_uses: Table<SolidUse>,
    compound_children: Table<CompoundChild>,
    product_children: Table<ProductChild>,

    edge_curves: Table<EdgeCurve>,
    edge_polylines: Table<EdgePolyline>,
    pcurves: Table<PCurve>,
    pcurve_polylines: Table<PCurvePolyline>,
    edge_mesh_polylines: Table<EdgeMeshPolyline>,
    face_surfaces: Table<FaceSurface>,
    face_meshes: Table<FaceMesh>,

    /// Next stable id, shared by all tables (never reused).
    next_stable_id: u32,
    uid: ModelUid,
    /// `None` until [`Model::enable_history`].
    history: Option<History>,
    /// Values derived from parts of the model, kept while those parts are unchanged.
    derived: derived::Derived,
}

// Operations read a model from several threads at once.
const _: fn() = || {
    fn shared<T: Send + Sync>() {}
    shared::<Model>();
};

/// A record type stored in one of the model's tables.
pub trait Record: Sized {
    fn table(model: &Model) -> &Table<Self>;
    fn table_mut(model: &mut Model) -> &mut Table<Self>;
}

macro_rules! records {
    ($($ty:ty => $field:ident),* $(,)?) => {
        $(
            impl Record for $ty {
                fn table(model: &Model) -> &Table<Self> {
                    &model.$field
                }
                fn table_mut(model: &mut Model) -> &mut Table<Self> {
                    &mut model.$field
                }
            }
        )*
    };
}

records! {
    Vertex => vertices, Edge => edges, CoEdge => coedges, Wire => wires, Face => faces,
    Shell => shells, Solid => solids, Compound => compounds, CompSolid => comp_solids,
    Product => products, Instance => instances,
    VertexUse => vertex_uses, WireUse => wire_uses, FaceUse => face_uses,
    ShellUse => shell_uses, SolidUse => solid_uses, CompoundChild => compound_children,
    ProductChild => product_children,
    EdgeCurve => edge_curves, EdgePolyline => edge_polylines, PCurve => pcurves,
    PCurvePolyline => pcurve_polylines, EdgeMeshPolyline => edge_mesh_polylines,
    FaceSurface => face_surfaces, FaceMesh => face_meshes,
}

impl Model {
    pub fn new() -> Self {
        Self::default()
    }

    /// True if the model holds no topology node.
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
            && self.edges.is_empty()
            && self.coedges.is_empty()
            && self.wires.is_empty()
            && self.faces.is_empty()
            && self.shells.is_empty()
            && self.solids.is_empty()
            && self.compounds.is_empty()
            && self.comp_solids.is_empty()
            && self.products.is_empty()
            && self.instances.is_empty()
    }

    /// The number of (non-removed) nodes of `kind`.
    pub fn count(&self, kind: NodeKind) -> usize {
        match kind {
            NodeKind::Instance => self.instances.len(),
            NodeKind::Product => self.products.len(),
            NodeKind::Compound => self.compounds.len(),
            NodeKind::CompSolid => self.comp_solids.len(),
            NodeKind::Solid => self.solids.len(),
            NodeKind::Shell => self.shells.len(),
            NodeKind::Face => self.faces.len(),
            NodeKind::Wire => self.wires.len(),
            NodeKind::CoEdge => self.coedges.len(),
            NodeKind::Edge => self.edges.len(),
            NodeKind::Vertex => self.vertices.len(),
        }
    }

    /// The table storing `T`.
    pub fn table<T: Record>(&self) -> &Table<T> {
        T::table(self)
    }

    /// The record for `id`, or `None` if it is invalid or removed.
    pub fn get<T: Record>(&self, id: Id<T>) -> Option<&T> {
        T::table(self).get(id)
    }

    /// The stable id of the record `id` (also for removed records).
    pub fn stable_id<T: Record>(&self, id: Id<T>) -> Option<StableId> {
        T::table(self).stable_id(id).map(StableId)
    }

    /// Low-level insertion that maintains no links between records; building topology
    /// should go through the higher-level builders, which keep parents and children
    /// consistent.
    pub fn insert<T: Record>(&mut self, record: T) -> Id<T> {
        let stable_id = self.next_stable_id;
        self.next_stable_id += 1;
        T::table_mut(self).push(record, stable_id)
    }

    /// The number of slots for nodes of `kind`, removed ones included: the exclusive upper
    /// bound of their indices.
    pub fn slot_count(&self, kind: NodeKind) -> usize {
        match kind {
            NodeKind::Instance => self.instances.capacity_len(),
            NodeKind::Product => self.products.capacity_len(),
            NodeKind::Compound => self.compounds.capacity_len(),
            NodeKind::CompSolid => self.comp_solids.capacity_len(),
            NodeKind::Solid => self.solids.capacity_len(),
            NodeKind::Shell => self.shells.capacity_len(),
            NodeKind::Face => self.faces.capacity_len(),
            NodeKind::Wire => self.wires.capacity_len(),
            NodeKind::CoEdge => self.coedges.capacity_len(),
            NodeKind::Edge => self.edges.capacity_len(),
            NodeKind::Vertex => self.vertices.capacity_len(),
        }
    }

    /// The stable id of a node of any kind (also for removed nodes).
    pub fn node_stable_id(&self, node: NodeId) -> Option<StableId> {
        self.node_slot(node).map(|(stable_id, _, _)| stable_id)
    }

    /// The stable id, revision and removed flag of a node's slot.
    fn node_slot(&self, node: NodeId) -> Option<(StableId, Revision, bool)> {
        macro_rules! slot {
            ($table:expr) => {{
                let table = $table;
                let id = Id::from_index(node.index);
                Some((StableId(table.stable_id(id)?), table.revision(id)?, table.is_removed(id)))
            }};
        }
        match node.kind {
            NodeKind::Instance => slot!(&self.instances),
            NodeKind::Product => slot!(&self.products),
            NodeKind::Compound => slot!(&self.compounds),
            NodeKind::CompSolid => slot!(&self.comp_solids),
            NodeKind::Solid => slot!(&self.solids),
            NodeKind::Shell => slot!(&self.shells),
            NodeKind::Face => slot!(&self.faces),
            NodeKind::Wire => slot!(&self.wires),
            NodeKind::CoEdge => slot!(&self.coedges),
            NodeKind::Edge => slot!(&self.edges),
            NodeKind::Vertex => slot!(&self.vertices),
        }
    }

    /// The live record with stable id `stable_id`, or `None` if there is none (removed,
    /// never issued, or of another type).
    #[doc(alias = "BRepGraph::UIDsView::NodeIdFrom")]
    pub fn find_stable<T: Record>(&self, stable_id: StableId) -> Option<Id<T>> {
        T::table(self).find_stable(stable_id.0)
    }

    /// The live node with stable id `stable_id`, of whatever kind.
    pub fn find_stable_node(&self, stable_id: StableId) -> Option<NodeId> {
        macro_rules! find {
            ($($table:ident => $kind:ident),*) => {{
                $(
                    if let Some(id) = self.$table.find_stable(stable_id.0) {
                        return Some(NodeId { kind: NodeKind::$kind, index: id.index() });
                    }
                )*
                None
            }};
        }
        find!(vertices => Vertex, edges => Edge, coedges => CoEdge, wires => Wire, faces => Face, shells => Shell,
            solids => Solid, compounds => Compound, comp_solids => CompSolid, products => Product, instances => Instance)
    }

    /// The model's modification history, if recording was turned on with
    /// [`Model::enable_history`].
    pub fn history(&self) -> Option<&History> {
        self.history.as_ref()
    }

    pub fn history_mut(&mut self) -> Option<&mut History> {
        self.history.as_mut()
    }

    /// Starts recording history (operations that record do nothing until then) and returns
    /// the history. Calling it again keeps the existing history.
    pub fn enable_history(&mut self) -> &mut History {
        self.history.get_or_insert_with(History::new)
    }

    /// Stops recording and returns the history recorded so far. To pause recording but
    /// keep the history, use [`History::set_enabled`] instead.
    pub fn disable_history(&mut self) -> Option<History> {
        self.history.take()
    }

    /// Runs `edit` on `target` and, if history is enabled, records that `target` became the
    /// nodes it returns (none: `target` was deleted). Nothing is recorded if `edit` panics.
    #[doc(alias = "BRepGraph::EditorView::GenOps::ApplyModification")]
    pub fn apply_modification<F>(&mut self, target: impl Into<NodeId>, operation: &str, edit: F)
    where
        F: FnOnce(&mut Model, NodeId) -> Vec<NodeId>,
    {
        let target = target.into();
        let images = edit(self, target);
        if let Some(history) = &mut self.history {
            history.record(operation, target, &images);
        }
    }

    /// The surface of `face`, or `None` if it has none (or `face` is invalid).
    pub fn face_surface(&self, face: FaceId) -> Option<&crate::geometry::surface::Surface> {
        let face = self.get(face)?;
        self.get(face.surface).map(|slot| &slot.surface)
    }

    /// Merges identical geometry (see [`DeduplicateOptions`] for topology).
    #[doc(alias = "BRepGraph_Deduplicate")]
    pub fn deduplicate(&mut self) -> DeduplicateReport {
        self.deduplicate_with(DeduplicateOptions::default())
    }

    pub fn deduplicate_with(&mut self, options: DeduplicateOptions) -> DeduplicateReport {
        self.deduplicate_impl(options)
    }

    /// Drops removed records from every table and remaps all stored ids (including those
    /// in the history). Invalidates outstanding [`Id`]s; [`StableId`]s and revisions are
    /// preserved.
    #[doc(alias = "BRepGraph_Compact")]
    pub fn compact(&mut self) -> CompactReport {
        let node_slots = |m: &Model| NodeKind::ALL.iter().map(|&k| m.slot_count(k)).sum::<usize>();
        let nodes_before = node_slots(self);
        let slots_before: Vec<usize> = NodeKind::ALL.iter().map(|&k| self.slot_count(k)).collect();
        let vertices = self.vertices.compact();
        let edges = self.edges.compact();
        let coedges = self.coedges.compact();
        let wires = self.wires.compact();
        let faces = self.faces.compact();
        let shells = self.shells.compact();
        let solids = self.solids.compact();
        let compounds = self.compounds.compact();
        let comp_solids = self.comp_solids.compact();
        let products = self.products.compact();
        let instances = self.instances.compact();
        let vertex_uses = self.vertex_uses.compact();
        let wire_uses = self.wire_uses.compact();
        let face_uses = self.face_uses.compact();
        let shell_uses = self.shell_uses.compact();
        let solid_uses = self.solid_uses.compact();
        let compound_children = self.compound_children.compact();
        let product_children = self.product_children.compact();
        let edge_curves = self.edge_curves.compact();
        let edge_polylines = self.edge_polylines.compact();
        let pcurves = self.pcurves.compact();
        let pcurve_polylines = self.pcurve_polylines.compact();
        let edge_mesh_polylines = self.edge_mesh_polylines.compact();
        let face_surfaces = self.face_surfaces.compact();
        let face_meshes = self.face_meshes.compact();
        fn at<T>(map: &[Id<T>], id: Id<T>) -> Id<T> {
            map.get(id.index() as usize).copied().unwrap_or(Id::INVALID)
        }
        let node_map = |node: NodeId| -> NodeId {
            let new = match node.kind {
                NodeKind::Vertex => at(&vertices, Id::from_index(node.index)).index(),
                NodeKind::Edge => at(&edges, Id::from_index(node.index)).index(),
                NodeKind::CoEdge => at(&coedges, Id::from_index(node.index)).index(),
                NodeKind::Wire => at(&wires, Id::from_index(node.index)).index(),
                NodeKind::Face => at(&faces, Id::from_index(node.index)).index(),
                NodeKind::Shell => at(&shells, Id::from_index(node.index)).index(),
                NodeKind::Solid => at(&solids, Id::from_index(node.index)).index(),
                NodeKind::Compound => at(&compounds, Id::from_index(node.index)).index(),
                NodeKind::CompSolid => at(&comp_solids, Id::from_index(node.index)).index(),
                NodeKind::Product => at(&products, Id::from_index(node.index)).index(),
                NodeKind::Instance => at(&instances, Id::from_index(node.index)).index(),
            };
            NodeId { kind: node.kind, index: new }
        };
        for e in self.edges.items_mut() {
            (e.curve, e.polyline, e.start, e.end) = (at(&edge_curves, e.curve), at(&edge_polylines, e.polyline), at(&vertex_uses, e.start), at(&vertex_uses, e.end));
        }
        for c in self.coedges.items_mut() {
            c.wire = at(&wires, c.wire);
            c.edge = at(&edges, c.edge);
            c.face = at(&faces, c.face);
            c.pcurve = at(&pcurves, c.pcurve);
            c.pcurve_polyline = at(&pcurve_polylines, c.pcurve_polyline);
            c.mesh_polyline = at(&edge_mesh_polylines, c.mesh_polyline);
        }
        for w in self.wires.items_mut() {
            w.coedges = w.coedges.iter().map(|&c| at(&coedges, c)).filter(|c| c.is_valid()).collect();
        }
        for f in self.faces.items_mut() {
            (f.surface, f.mesh) = (at(&face_surfaces, f.surface), at(&face_meshes, f.mesh));
            f.wires = f.wires.iter().map(|&u| at(&wire_uses, u)).filter(|u| u.is_valid()).collect();
        }
        for s in self.shells.items_mut() {
            s.faces = s.faces.iter().map(|&u| at(&face_uses, u)).filter(|u| u.is_valid()).collect();
        }
        for s in self.solids.items_mut() {
            s.shells = s.shells.iter().map(|&u| at(&shell_uses, u)).filter(|u| u.is_valid()).collect();
        }
        for c in self.compounds.items_mut() {
            c.children = c.children.iter().map(|&u| at(&compound_children, u)).filter(|u| u.is_valid()).collect();
        }
        for c in self.comp_solids.items_mut() {
            c.solids = c.solids.iter().map(|&u| at(&solid_uses, u)).filter(|u| u.is_valid()).collect();
        }
        for p in self.products.items_mut() {
            if let Some(shape) = &mut p.shape {
                shape.node = node_map(shape.node);
            }
            p.children = p.children.iter().map(|&u| at(&product_children, u)).filter(|u| u.is_valid()).collect();
        }
        for i in self.instances.items_mut() {
            i.target = node_map(i.target);
        }
        for u in self.vertex_uses.items_mut() {
            (u.edge, u.vertex) = (at(&edges, u.edge), at(&vertices, u.vertex));
        }
        for u in self.wire_uses.items_mut() {
            (u.face, u.wire) = (at(&faces, u.face), at(&wires, u.wire));
        }
        for u in self.face_uses.items_mut() {
            (u.shell, u.face) = (at(&shells, u.shell), at(&faces, u.face));
        }
        for u in self.shell_uses.items_mut() {
            (u.solid, u.shell) = (at(&solids, u.solid), at(&shells, u.shell));
        }
        for u in self.solid_uses.items_mut() {
            (u.comp_solid, u.solid) = (at(&comp_solids, u.comp_solid), at(&solids, u.solid));
        }
        for u in self.compound_children.items_mut() {
            (u.compound, u.child) = (at(&compounds, u.compound), node_map(u.child));
        }
        for u in self.product_children.items_mut() {
            (u.product, u.instance) = (at(&products, u.product), at(&instances, u.instance));
        }
        for s in self.edge_curves.items_mut() {
            s.edge = at(&edges, s.edge);
        }
        for s in self.edge_polylines.items_mut() {
            s.edge = at(&edges, s.edge);
        }
        for s in self.pcurves.items_mut() {
            s.coedge = at(&coedges, s.coedge);
        }
        for s in self.pcurve_polylines.items_mut() {
            s.coedge = at(&coedges, s.coedge);
        }
        for s in self.edge_mesh_polylines.items_mut() {
            s.coedge = at(&coedges, s.coedge);
        }
        for s in self.face_surfaces.items_mut() {
            s.face = at(&faces, s.face);
        }
        for s in self.face_meshes.items_mut() {
            s.face = at(&faces, s.face);
        }
        if let Some(history) = &mut self.history {
            history.remap(&node_map);
        }
        self.derived.clear();
        let mut renumbered = 0;
        macro_rules! count_moved {
            ($($map:ident),*) => {$(
                renumbered += $map.iter().enumerate().filter(|(i, id)| id.is_valid() && id.index() as usize != *i).count();
            )*};
        }
        count_moved!(vertices, edges, coedges, wires, faces, shells, solids, compounds, comp_solids, products, instances);
        let mut removed = [0; 11];
        for (i, &kind) in NodeKind::ALL.iter().enumerate() {
            removed[i] = slots_before[i] - self.slot_count(kind);
        }
        CompactReport { nodes_before, nodes_after: node_slots(self), removed, renumbered }
    }
}

/// Options for [`Model::deduplicate_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeduplicateOptions {
    /// Only count what would be merged; leave the model unchanged.
    pub analyze_only: bool,
    /// Record entity merges in the model's history, if it is enabled: one `"deduplicate"`
    /// event per merged node, [`ChangeKind::Replaced`] by the node it was merged into. Sharing geometry is not a
    /// topological change and is never recorded.
    pub record_history: bool,
    /// Also merge topology nodes that become identical once their geometry is shared
    /// (vertices at the same point, edges on the same curve between the same vertices,
    /// wires of the same coedges, faces on the same surface with the same wires), where
    /// this does not change the shape. Merged nodes are removed and their uses redirected.
    #[doc(alias = "MergeEntitiesWhenSafe")]
    pub merge_entities: bool,
    /// Distance under which geometry is considered equal.
    pub linear_tolerance: f64,
    /// Angle under which directions are considered equal.
    pub angular_tolerance: f64,
}

impl Default for DeduplicateOptions {
    fn default() -> Self {
        Self {
            analyze_only: false,
            record_history: true,
            merge_entities: false,
            linear_tolerance: crate::tolerance::LINEAR,
            angular_tolerance: crate::tolerance::ANGULAR,
        }
    }
}

/// What [`Model::deduplicate_with`] found and merged.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DeduplicateReport {
    /// Distinct surfaces among the faces' surfaces.
    #[doc(alias = "NbCanonicalSurfaces")]
    pub distinct_surfaces: usize,
    /// Distinct curves among the edges' curves.
    #[doc(alias = "NbCanonicalCurves")]
    pub distinct_curves: usize,
    /// Faces redirected to an equal surface of another face.
    #[doc(alias = "NbSurfaceRewrites")]
    pub shared_surfaces: usize,
    /// Edges redirected to an equal curve of another edge.
    #[doc(alias = "NbCurveRewrites")]
    pub shared_curves: usize,
    /// Topology nodes merged into an identical one, per kind (indexed as
    /// [`NodeKind::ALL`]); with `analyze_only`, those that would be.
    pub merged: [usize; 11],
    /// True if nodes were actually merged (false with `analyze_only`).
    #[doc(alias = "IsEntityMergeApplied")]
    pub merge_applied: bool,
    /// Events added to the history.
    pub history_events: usize,
}

impl DeduplicateReport {
    pub fn merged(&self, kind: NodeKind) -> usize {
        self.merged[kind as usize]
    }
}

/// What [`Model::compact`] dropped.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CompactReport {
    /// Topology nodes (removed ones included) before, and nodes after.
    pub nodes_before: usize,
    pub nodes_after: usize,
    /// Removed nodes dropped, per kind (indexed as [`NodeKind::ALL`]).
    pub removed: [usize; 11],
    /// Surviving nodes whose id changed.
    pub renumbered: usize,
}

impl CompactReport {
    pub fn removed(&self, kind: NodeKind) -> usize {
        self.removed[kind as usize]
    }
}

impl<T: Record> std::ops::Index<Id<T>> for Model {
    type Output = T;
    /// Panics on an invalid or removed id.
    fn index(&self, id: Id<T>) -> &T {
        &T::table(self)[id]
    }
}

impl<T: Record> std::ops::IndexMut<Id<T>> for Model {
    fn index_mut(&mut self, id: Id<T>) -> &mut T {
        &mut T::table_mut(self)[id]
    }
}
