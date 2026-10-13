//! Hand translations of OCCT's TKBRep GTests (src/ModelingData/TKBRep/GTests).

mod brep_adaptor_comp_curve;
mod brep_graph_assembly;
mod brep_graph_batch_ops;
mod brep_graph_build;
mod brep_graph_child_explorer;
mod brep_graph_children;
mod brep_graph_compact;
mod brep_graph_copy;
mod brep_graph_deduplicate;
mod brep_graph_editing;
mod brep_graph_geometry;
mod brep_graph_history;
mod brep_graph_ids;
mod brep_graph_meshes;
mod brep_graph_misc;
mod brep_graph_iterator;
mod brep_graph_parent_explorer;
mod brep_graph_parents;
mod brep_graph_related;
mod brep_graph_reverse;
mod brep_graph_revisions;
mod brep_graph_seams;
mod brep_graph_sharing;
mod brep_graph_sparse_model;
mod brep_graph_storage;
mod brep_graph_tool;
mod brep_graph_transform;
mod brep_graph_validate;
mod brep_graph_views;
mod brep_tool;
mod brep_tools;
mod brep_tools_gtrsf_modification;
mod brep_tools_reshape;
mod brep_tools_trsf_modification;
mod graph_support;
mod modification_support;
mod top_exp;
mod topo_ds_builder;
mod topo_ds_iterator;
mod topo_ds_tshape;

use kentron::math::Point;
use kentron::primitives::Cuboid;
use kentron::topology::{Model, NodeKind, Shape};
use kentron::make;

/// A box built in `model`.
pub fn make_box(model: &mut Model, dx: f64, dy: f64, dz: f64) -> Shape {
    Cuboid::new(dx, dy, dz).build(model).expect("box creation failed")
}

/// A straight edge between two points.
pub fn make_edge(model: &mut Model, p1: Point, p2: Point) -> Shape {
    make::edge_between(model, p1, p2).expect("edge creation failed")
}

/// The first `kind` occurrence below `shape`.
pub fn first(model: &Model, shape: Shape, kind: NodeKind) -> Shape {
    model.explore(shape, kind).next().unwrap_or_else(|| panic!("no {kind:?} found"))
}
