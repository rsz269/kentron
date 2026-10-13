//! Collecting garbage: removing what is not reachable from the shapes still in use, and
//! freeing what the removed records held.
//!
//! A model that keeps every intermediate result of a long sequence of operations (as a
//! scripting session does) grows without bound; [`Model::collect_garbage`] takes the roots
//! still in use and drops the rest. Nothing is renumbered, so ids of what is kept stay
//! valid; the slots of what was dropped are reused by records made later.

use crate::geometry::curve::{Curve, Line};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{Plane, Surface};
use crate::math::{Direction, Direction2, Point, Point2};

use super::id::{Id, NodeId, NodeKind};
use super::table::Table;
use super::Model;

/// One mark per slot of a table.
struct Marks<T> {
    marked: Vec<bool>,
    _table: std::marker::PhantomData<T>,
}

impl<T> Marks<T> {
    fn new(table: &Table<T>) -> Self {
        Self { marked: vec![false; table.capacity_len()], _table: std::marker::PhantomData }
    }

    /// Marks `id`; false if it was marked already (or is invalid).
    fn mark(&mut self, id: Id<T>) -> bool {
        match self.marked.get_mut(id.index() as usize) {
            Some(m) if !*m => {
                *m = true;
                true
            }
            _ => false,
        }
    }

    fn is_marked(&self, id: Id<T>) -> bool {
        self.marked.get(id.index() as usize).copied().unwrap_or(false)
    }
}

/// Soft-removes every unmarked live entry of `table`, emptying it with `release`; returns
/// how many were removed.
fn sweep<T>(table: &mut Table<T>, marks: &Marks<T>, release: impl Fn(&mut T)) -> usize {
    let dead: Vec<Id<T>> = table.iter().map(|(id, _)| id).filter(|&id| !marks.is_marked(id)).collect();
    for &id in &dead {
        table.remove_releasing(id, &release);
    }
    dead.len()
}

/// The release of a record that holds nothing on the heap.
fn unchanged<T>(_: &mut T) {}

impl Model {
    /// Removes every node that is neither one of `roots` nor reachable from one (as a
    /// child, or through a coedge's wire and face), with the uses, geometry and meshes only
    /// such nodes have, and frees what the removed records held; their slots are reused by
    /// records made later, so ids of removed nodes must not be kept. Ids of what is kept stay
    /// valid. Returns the number of nodes removed.
    ///
    /// Values derived from removed nodes stay in the model's cache (bounded on its own; a
    /// record made in a reused slot does not match them): see
    /// [`CollectOptions::forget_derived`].
    pub fn collect_garbage(&mut self, roots: impl IntoIterator<Item = NodeId>) -> usize {
        self.collect_garbage_with(roots, CollectOptions::default())
    }

    pub fn collect_garbage_with(&mut self, roots: impl IntoIterator<Item = NodeId>, options: CollectOptions) -> usize {
        let mut vertices = Marks::new(&self.vertices);
        let mut edges = Marks::new(&self.edges);
        let mut coedges = Marks::new(&self.coedges);
        let mut wires = Marks::new(&self.wires);
        let mut faces = Marks::new(&self.faces);
        let mut shells = Marks::new(&self.shells);
        let mut solids = Marks::new(&self.solids);
        let mut compounds = Marks::new(&self.compounds);
        let mut comp_solids = Marks::new(&self.comp_solids);
        let mut products = Marks::new(&self.products);
        let mut instances = Marks::new(&self.instances);
        let mut vertex_uses = Marks::new(&self.vertex_uses);
        let mut wire_uses = Marks::new(&self.wire_uses);
        let mut face_uses = Marks::new(&self.face_uses);
        let mut shell_uses = Marks::new(&self.shell_uses);
        let mut solid_uses = Marks::new(&self.solid_uses);
        let mut compound_children = Marks::new(&self.compound_children);
        let mut product_children = Marks::new(&self.product_children);
        let mut edge_curves = Marks::new(&self.edge_curves);
        let mut edge_polylines = Marks::new(&self.edge_polylines);
        let mut pcurves = Marks::new(&self.pcurves);
        let mut pcurve_polylines = Marks::new(&self.pcurve_polylines);
        let mut edge_mesh_polylines = Marks::new(&self.edge_mesh_polylines);
        let mut face_surfaces = Marks::new(&self.face_surfaces);
        let mut face_meshes = Marks::new(&self.face_meshes);

        let mut stack: Vec<NodeId> = roots.into_iter().collect();
        while let Some(node) = stack.pop() {
            let index = node.index;
            match node.kind {
                NodeKind::Vertex => {
                    vertices.mark(Id::from_index(index));
                }
                NodeKind::Edge => {
                    let id = Id::from_index(index);
                    let Some(record) = self.edges.get(id).filter(|_| edges.mark(id)) else { continue };
                    edge_curves.mark(record.curve);
                    edge_polylines.mark(record.polyline);
                    for u in [record.start, record.end] {
                        if let Some(used) = self.vertex_uses.get(u).filter(|_| vertex_uses.mark(u)) {
                            stack.push(used.vertex.into());
                        }
                    }
                }
                NodeKind::CoEdge => {
                    let id = Id::from_index(index);
                    let Some(record) = self.coedges.get(id).filter(|_| coedges.mark(id)) else { continue };
                    pcurves.mark(record.pcurve);
                    pcurve_polylines.mark(record.pcurve_polyline);
                    edge_mesh_polylines.mark(record.mesh_polyline);
                    stack.push(record.edge.into());
                    // (What a coedge is on stays, as its pcurve means nothing without it.)
                    if record.wire.is_valid() {
                        stack.push(record.wire.into());
                    }
                    if record.face.is_valid() {
                        stack.push(record.face.into());
                    }
                }
                NodeKind::Wire => {
                    let id = Id::from_index(index);
                    let Some(record) = self.wires.get(id).filter(|_| wires.mark(id)) else { continue };
                    stack.extend(record.coedges.iter().map(|&c| NodeId::from(c)));
                }
                NodeKind::Face => {
                    let id = Id::from_index(index);
                    let Some(record) = self.faces.get(id).filter(|_| faces.mark(id)) else { continue };
                    face_surfaces.mark(record.surface);
                    face_meshes.mark(record.mesh);
                    for &u in &record.wires {
                        if let Some(used) = self.wire_uses.get(u).filter(|_| wire_uses.mark(u)) {
                            stack.push(used.wire.into());
                        }
                    }
                }
                NodeKind::Shell => {
                    let id = Id::from_index(index);
                    let Some(record) = self.shells.get(id).filter(|_| shells.mark(id)) else { continue };
                    for &u in &record.faces {
                        if let Some(used) = self.face_uses.get(u).filter(|_| face_uses.mark(u)) {
                            stack.push(used.face.into());
                        }
                    }
                }
                NodeKind::Solid => {
                    let id = Id::from_index(index);
                    let Some(record) = self.solids.get(id).filter(|_| solids.mark(id)) else { continue };
                    for &u in &record.shells {
                        if let Some(used) = self.shell_uses.get(u).filter(|_| shell_uses.mark(u)) {
                            stack.push(used.shell.into());
                        }
                    }
                }
                NodeKind::CompSolid => {
                    let id = Id::from_index(index);
                    let Some(record) = self.comp_solids.get(id).filter(|_| comp_solids.mark(id)) else { continue };
                    for &u in &record.solids {
                        if let Some(used) = self.solid_uses.get(u).filter(|_| solid_uses.mark(u)) {
                            stack.push(used.solid.into());
                        }
                    }
                }
                NodeKind::Compound => {
                    let id = Id::from_index(index);
                    let Some(record) = self.compounds.get(id).filter(|_| compounds.mark(id)) else { continue };
                    for &u in &record.children {
                        if let Some(used) = self.compound_children.get(u).filter(|_| compound_children.mark(u)) {
                            stack.push(used.child);
                        }
                    }
                }
                NodeKind::Product => {
                    let id = Id::from_index(index);
                    let Some(record) = self.products.get(id).filter(|_| products.mark(id)) else { continue };
                    stack.extend(record.shape.map(|s| s.node));
                    for &u in &record.children {
                        if let Some(used) = self.product_children.get(u).filter(|_| product_children.mark(u)) {
                            stack.push(used.instance.into());
                        }
                    }
                }
                NodeKind::Instance => {
                    let id = Id::from_index(index);
                    let Some(record) = self.instances.get(id).filter(|_| instances.mark(id)) else { continue };
                    stack.push(record.target);
                }
            }
        }

        // A slot shared by several records (see `Model::deduplicate`) names one of them as
        // its owner: one that is kept, if that one goes.
        for (id, record) in self.edges.iter().filter(|&(id, _)| edges.is_marked(id)) {
            if self.edge_curves.get(record.curve).is_some_and(|s| !edges.is_marked(s.edge)) {
                self.edge_curves.get_mut_keeping_content(record.curve).expect("a live slot").edge = id;
            }
            if self.edge_polylines.get(record.polyline).is_some_and(|s| !edges.is_marked(s.edge)) {
                self.edge_polylines.get_mut_keeping_content(record.polyline).expect("a live slot").edge = id;
            }
        }
        for (id, record) in self.coedges.iter().filter(|&(id, _)| coedges.is_marked(id)) {
            if self.pcurves.get(record.pcurve).is_some_and(|s| !coedges.is_marked(s.coedge)) {
                self.pcurves.get_mut_keeping_content(record.pcurve).expect("a live slot").coedge = id;
            }
            if self.pcurve_polylines.get(record.pcurve_polyline).is_some_and(|s| !coedges.is_marked(s.coedge)) {
                self.pcurve_polylines.get_mut_keeping_content(record.pcurve_polyline).expect("a live slot").coedge = id;
            }
            if self.edge_mesh_polylines.get(record.mesh_polyline).is_some_and(|s| !coedges.is_marked(s.coedge)) {
                self.edge_mesh_polylines.get_mut_keeping_content(record.mesh_polyline).expect("a live slot").coedge = id;
            }
        }
        for (id, record) in self.faces.iter().filter(|&(id, _)| faces.is_marked(id)) {
            if self.face_surfaces.get(record.surface).is_some_and(|s| !faces.is_marked(s.face)) {
                self.face_surfaces.get_mut_keeping_content(record.surface).expect("a live slot").face = id;
            }
            if self.face_meshes.get(record.mesh).is_some_and(|s| !faces.is_marked(s.face)) {
                self.face_meshes.get_mut_keeping_content(record.mesh).expect("a live slot").face = id;
            }
        }

        let mut removed = 0;
        removed += sweep(&mut self.vertices, &vertices, unchanged);
        removed += sweep(&mut self.edges, &edges, unchanged);
        removed += sweep(&mut self.coedges, &coedges, unchanged);
        removed += sweep(&mut self.wires, &wires, |r| r.coedges = Vec::new());
        removed += sweep(&mut self.faces, &faces, |r| r.wires = Vec::new());
        removed += sweep(&mut self.shells, &shells, |r| r.faces = Vec::new());
        removed += sweep(&mut self.solids, &solids, |r| r.shells = Vec::new());
        removed += sweep(&mut self.compounds, &compounds, |r| r.children = Vec::new());
        removed += sweep(&mut self.comp_solids, &comp_solids, |r| r.solids = Vec::new());
        removed += sweep(&mut self.products, &products, |r| *r = Default::default());
        removed += sweep(&mut self.instances, &instances, unchanged);
        sweep(&mut self.vertex_uses, &vertex_uses, unchanged);
        sweep(&mut self.wire_uses, &wire_uses, unchanged);
        sweep(&mut self.face_uses, &face_uses, unchanged);
        sweep(&mut self.shell_uses, &shell_uses, unchanged);
        sweep(&mut self.solid_uses, &solid_uses, unchanged);
        sweep(&mut self.compound_children, &compound_children, unchanged);
        sweep(&mut self.product_children, &product_children, unchanged);
        sweep(&mut self.edge_curves, &edge_curves, |s| s.curve = Curve::Line(Line { origin: Point::ORIGIN, direction: Direction::X }));
        sweep(&mut self.edge_polylines, &edge_polylines, |s| s.polyline = Default::default());
        sweep(&mut self.pcurves, &pcurves, |s| s.curve = Curve2::Line(Line2 { origin: Point2::ORIGIN, direction: Direction2::X }));
        sweep(&mut self.pcurve_polylines, &pcurve_polylines, |s| s.polyline = Default::default());
        sweep(&mut self.edge_mesh_polylines, &edge_mesh_polylines, |s| s.polyline = Default::default());
        sweep(&mut self.face_surfaces, &face_surfaces, |s| s.surface = Surface::Plane(Plane::XY));
        sweep(&mut self.face_meshes, &face_meshes, |s| s.mesh = Default::default());

        if options.forget_derived && removed > 0 {
            let marked = |node: NodeId| {
                let id = node.index;
                match node.kind {
                    NodeKind::Vertex => vertices.is_marked(Id::from_index(id)),
                    NodeKind::Edge => edges.is_marked(Id::from_index(id)),
                    NodeKind::CoEdge => coedges.is_marked(Id::from_index(id)),
                    NodeKind::Wire => wires.is_marked(Id::from_index(id)),
                    NodeKind::Face => faces.is_marked(Id::from_index(id)),
                    NodeKind::Shell => shells.is_marked(Id::from_index(id)),
                    NodeKind::Solid => solids.is_marked(Id::from_index(id)),
                    NodeKind::Compound => compounds.is_marked(Id::from_index(id)),
                    NodeKind::CompSolid => comp_solids.is_marked(Id::from_index(id)),
                    NodeKind::Product => products.is_marked(Id::from_index(id)),
                    NodeKind::Instance => instances.is_marked(Id::from_index(id)),
                }
            };
            self.derived.retain_nodes(marked);
        }
        removed
    }
}

/// Options for [`Model::collect_garbage_with`].
///
/// pykentron's profiles (`PYKENTRON_PROFILE`) are an allocator and these options; a Rust
/// program has the same choices, the allocator by declaring it:
///
/// | profile | global allocator | collection |
/// |---|---|---|
/// | `low-memory` | the system's (none declared) | `forget_derived: true` |
/// | `balanced` | the system's (none declared) | the defaults |
/// | `fast` | `#[global_allocator] static A: mimalloc::MiMalloc = mimalloc::MiMalloc;` | the defaults |
///
/// On a complex five-part CadQuery model (pykentron), against no collection: `low-memory` peaked at
/// 62% of the memory and took ~5–10% longer, `balanced` 72% and ~3% longer, `fast` the same
/// memory and ~17% less time.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CollectOptions {
    /// Also forget the values derived from removed nodes (sampled boundaries, boxes, ...),
    /// which otherwise stay until the cache is full. Frees more memory, but measured on
    /// CadQuery scripts it makes later operations 5–10% slower: the freed memory, scattered
    /// through the heap, is slower for the allocator to reuse.
    pub forget_derived: bool,
}

#[cfg(test)]
mod tests {
    use crate::primitives::Cuboid;
    use super::CollectOptions;
    use crate::topology::{Model, NodeKind};
    use crate::{properties, validity};

    fn counts(model: &Model) -> [usize; 4] {
        [NodeKind::Solid, NodeKind::Face, NodeKind::Edge, NodeKind::Vertex].map(|k| model.count(k))
    }

    #[test]
    fn keeps_what_the_roots_reach_and_removes_the_rest() {
        let mut model = Model::new();
        let kept = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let dropped = Cuboid::new(4.0, 5.0, 6.0).build(&mut model).unwrap();
        let face_of_dropped = model.explore(dropped, NodeKind::Face).next().unwrap();
        assert_eq!(counts(&model), [2, 12, 24, 16]);

        let removed = model.collect_garbage([kept.node, face_of_dropped.node]);
        // The face keeps its wire, 4 edges and 4 vertices; the rest of the second box goes.
        assert_eq!(counts(&model), [1, 7, 16, 12]);
        assert!(removed > 0);
        assert!(!model.is_live(dropped.node));
        assert!(model.is_live(face_of_dropped.node));
        assert!(validity::check(&model, kept).is_valid());
        assert!((properties::volume(&model, kept).mass - 6.0).abs() < 1e-9);
        assert_eq!(model.explore(face_of_dropped, NodeKind::Edge).count(), 4);

        // Collecting again with the same roots removes nothing.
        assert_eq!(model.collect_garbage([kept.node, face_of_dropped.node]), 0);
    }

    #[test]
    fn values_derived_from_removed_nodes_can_be_forgotten() {
        for forget_derived in [false, true] {
            let mut model = Model::new();
            let kept = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
            let dropped = Cuboid::new(4.0, 5.0, 6.0).build(&mut model).unwrap();
            for shape in [kept, dropped] {
                for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
                    model.derived(face.node, 0, || 1u8);
                }
            }
            model.collect_garbage_with([kept.node], CollectOptions { forget_derived });
            assert_eq!(model.derived_count(), if forget_derived { 6 } else { 12 });
        }
    }

    #[test]
    fn records_made_later_reuse_the_slots_of_removed_ones() {
        let mut model = Model::new();
        let kept = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        Cuboid::new(4.0, 5.0, 6.0).build(&mut model).unwrap();
        model.collect_garbage([kept.node]);
        let slots: Vec<usize> = NodeKind::ALL.iter().map(|&k| model.slot_count(k)).collect();

        let new = Cuboid::new(7.0, 8.0, 9.0).build(&mut model).unwrap();
        assert_eq!(NodeKind::ALL.iter().map(|&k| model.slot_count(k)).collect::<Vec<_>>(), slots, "no new slots");
        assert_eq!(counts(&model), [2, 12, 24, 16]);
        for (shape, volume) in [(kept, 6.0), (new, 504.0)] {
            assert!(validity::check(&model, shape).is_valid());
            assert!((properties::volume(&model, shape).mass - volume).abs() < 1e-9);
        }
    }
}
