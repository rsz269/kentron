//! BRepGraph_DefsIterator_Test.cxx and BRepGraph_RefsIterator_Test.cxx
//!
//! Single-level children are read from the entity records (`model[shell].faces`, then
//! `model[face_use].face`). Not ported: `VertexOfEdge_EnumeratesBoundaryVerticesOnly` /
//! `VertexOfEdge_ExposesBoundaryVertexRefsOnly` (kentron edges have no internal vertices)
//! and `ChildOfCompound_SkipsOutOfRangeChildNode` (both files; they corrupt records
//! through OCCT's mutation guards).

use kentron::math::{Point, Transform};
use kentron::topology::{
    Compound, CompSolid, Edge, Face, Model, NodeKind, Shape, Shell, Solid, Wire,
};

use crate::graph_support::{box_model, compound, first, id_of};

#[test]
fn box_hierarchy_traverses_single_level_children() {
    let (model, _) = box_model();
    assert_eq!(model[first::<Solid>(&model)].shells.len(), 1);
    assert_eq!(model[first::<Shell>(&model)].faces.len(), 6);
    assert_eq!(model[first::<Face>(&model)].wires.len(), 1);
    assert_eq!(model[first::<Wire>(&model)].coedges.len(), 4);
    let edge = &model[first::<Edge>(&model)];
    assert!(model.get(edge.start).is_some() && model.get(edge.end).is_some());
}

#[test]
fn edge_of_wire_yields_edge_definitions() {
    let (model, _) = box_model();
    let coedge_id = model[first::<Wire>(&model)].coedges[0];
    let coedge = &model[coedge_id];
    let edge = model.get(coedge.edge).expect("edge");
    assert!(model.get(edge.start).is_some());
    assert!(model.get(edge.end).is_some());
}

#[test]
fn co_edge_of_wire_yields_co_edge_definitions() {
    let (model, _) = box_model();
    let coedge_id = model[first::<Wire>(&model)].coedges[0];
    assert!(model.get(model[coedge_id].edge).is_some());
}

#[test]
fn vertex_of_edge_exposes_boundary_vertex_ref() {
    let (model, _) = box_model();
    let edge_id = first::<Edge>(&model);
    let start_use = &model[model[edge_id].start];
    assert_eq!(start_use.edge, edge_id);
    assert!(model.get(start_use.vertex).is_some());
}

/// The fixture box plus a compound of its solid and a loose vertex.
fn box_and_vertex_compound() -> (Model, Shape, Shape) {
    let (mut model, solid) = box_model();
    let vertex = Shape::from(model.add_vertex(Point::new(1.0, 2.0, 3.0), 0.01));
    let c = compound(&mut model, &[solid, vertex]);
    (model, c, vertex)
}

#[test]
fn child_of_compound_enumerates_heterogeneous_children() {
    let (model, c, _) = box_and_vertex_compound();
    let kinds: Vec<NodeKind> = model[id_of::<Compound>(c)]
        .children
        .iter()
        .map(|&child| model[child].child.kind)
        .collect();
    assert_eq!(kinds, [NodeKind::Solid, NodeKind::Vertex]);
}

#[test]
fn child_of_compound_skips_removed_child_node() {
    let (mut model, c, vertex) = box_and_vertex_compound();
    model.remove_node(vertex.node);
    assert_eq!(model[id_of::<Compound>(c)].children.len(), 1);
    assert_eq!(model.children(c).count(), 1);
}

#[test]
fn solid_of_comp_solid_enumerates_direct_solids() {
    let (mut model, _) = box_model();
    let solid_a = model.new_empty(NodeKind::Solid);
    let solid_b = model.new_empty(NodeKind::Solid);
    let comp_solid = model.new_empty(NodeKind::CompSolid);
    model.add(comp_solid, solid_a).unwrap();
    model.add(comp_solid, solid_b).unwrap();
    assert_eq!(model[id_of::<CompSolid>(comp_solid)].solids.len(), 2);
}

#[test]
fn occurrence_of_product_enumerates_direct_occurrences() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid.node).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    assert_eq!(model[assembly].children.len(), 2);
}

#[test]
fn removed_wire_ref_is_skipped() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let wire_uses = model[face].wires.clone();
    assert_eq!(wire_uses.len(), 1);
    model.remove_use(wire_uses[0]);
    assert!(model[face].wires.is_empty());
}

// RefsIterator-specific checks.

#[test]
fn box_hierarchy_yields_reference_ids() {
    let (model, _) = box_model();
    let solid = &model[first::<Solid>(&model)];
    assert!(solid.shells.iter().all(|&u| model.get(u).is_some()));
    let shell = &model[first::<Shell>(&model)];
    assert_eq!(shell.faces.iter().filter(|&&u| model.get(u).is_some()).count(), 6);
}

#[test]
fn current_id_resolves_to_expected_entry() {
    let (model, _) = box_model();
    let wire_use = model[first::<Face>(&model)].wires[0];
    assert!(model.get(model[wire_use].wire).is_some());
}

#[test]
fn child_of_compound_enumerates_child_refs() {
    let (model, c, _) = box_and_vertex_compound();
    let children = &model[id_of::<Compound>(c)].children;
    assert_eq!(model[children[0]].child.kind, NodeKind::Solid);
    assert_eq!(model[children[1]].child.kind, NodeKind::Vertex);
}
