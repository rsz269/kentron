//! BRepGraph_ShapesViewImport_Test.cxx (programmatic construction through the editor, in
//! spite of its name) and BRepGraph_ScenarioMatrix_Test.cxx (mixed edit/validate
//! scenarios).
//!
//! Not ported: mutation-guard tests (`Mutable*Definition`, `MutInvalidTopologyDefs_*`),
//! `AddWireToFace_*` (the same edges in two wires of one face), `AddWithHistory_*`
//! (TopoDS import with BRepTools_History), and the ScenarioMatrix cases that reconstruct
//! or re-import TopoDS shapes.

use std::collections::HashSet;

use kentron::geometry::curve::{Curve, Line};
use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Point, Transform, Vector};
use kentron::primitives::{Cuboid, Cylinder};
use kentron::topology::{
    CoEdge, Edge, EdgeId, Face, FaceId, Model, NodeId, NodeKind, Orientation, Shape, Shell,
    Solid, Vertex, VertexId, WireId,
};
use kentron::validity::{self, ProblemKind, Severity};

use crate::graph_support::{box_model, compound, first, id_of};

/// A 10 x 10 square on the XY plane: (vertices, edges, wire).
fn square(model: &mut Model) -> (Vec<VertexId>, Vec<EdgeId>, WireId) {
    let p = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    let v: Vec<VertexId> = p.iter().map(|&(x, y)| model.add_vertex(Point::new(x, y, 0.0), 0.001)).collect();
    let e: Vec<EdgeId> = (0..4)
        .map(|i| {
            let (a, b) = (v[i], v[(i + 1) % 4]);
            let edge = model.add_edge(a, b, 0.001).unwrap();
            let (pa, pb) = (model[a].point, model[b].point);
            let line = Curve::Line(Line::new(pa, Direction::from_vector(pb - pa)));
            model.set_edge_curve(edge.into(), line, 0.001).unwrap();
            model.set_edge_range(edge.into(), 0.0, 10.0).unwrap();
            edge
        })
        .collect();
    let coedges: Vec<_> = e.iter().map(|&e| model.add_coedge(e, Orientation::Forward).unwrap()).collect();
    let wire = model.add_wire(&coedges).unwrap();
    (v, e, wire)
}

fn square_face(model: &mut Model) -> (Vec<VertexId>, Vec<EdgeId>, WireId, FaceId) {
    let (v, e, wire) = square(model);
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, Surface::Plane(Plane::new(Point::ORIGIN, Direction::Z)), 0.001).unwrap();
    model.add(face, wire.into()).unwrap();
    (v, e, wire, id_of(face))
}

// ---- construction ----

#[test]
fn add_vertex_and_edge() {
    let mut model = Model::new();
    let (v, e, _) = square(&mut model);
    assert_eq!(model.count(NodeKind::Vertex), 4);
    assert_eq!(model.count(NodeKind::Edge), 4);
    assert_eq!(model.point(v[2]), Point::new(10.0, 10.0, 0.0));
    assert!(model.edge_curve(e[0]).is_some());
    assert_eq!(model.edge_range(e[0]), (0.0, 10.0));
}

#[test]
fn add_wire_closed_rectangle() {
    let mut model = Model::new();
    let (_, _, wire) = square(&mut model);
    assert_eq!(model[wire].coedges.len(), 4);
    assert!(model.is_closed_wire(wire));
}

#[test]
fn add_edge_with_invalid_vertex_fails_without_mutation() {
    let mut model = Model::new();
    let v = model.add_vertex(Point::ORIGIN, 0.001);
    assert!(model.add_edge(v, VertexId::from_index(42), 0.001).is_err());
    assert_eq!(model.slot_count(NodeKind::Edge), 0);
}

#[test]
fn add_coedge_with_invalid_edge_fails_without_mutation() {
    let mut model = Model::new();
    assert!(model.add_coedge(EdgeId::from_index(17), Orientation::Forward).is_err());
    assert_eq!(model.slot_count(NodeKind::CoEdge), 0);
    assert!(model.add_wire(&[kentron::topology::CoEdgeId::from_index(3)]).is_err());
    assert_eq!(model.slot_count(NodeKind::Wire), 0);
}

#[test]
fn adding_an_invalid_child_fails_without_mutation() {
    let mut model = Model::new();
    let face = model.new_empty(NodeKind::Face);
    assert!(model.add(face, WireId::from_index(9).into()).is_err());
    assert!(model[id_of::<Face>(face)].wires.is_empty());
    let shell = model.new_empty(NodeKind::Shell);
    assert!(model.add(shell, FaceId::from_index(4).into()).is_err());
    assert!(model[id_of::<Shell>(shell)].faces.is_empty());
}

#[test]
fn full_solid_programmatic_construction() {
    let mut model = Model::new();
    let (_, _, _, face) = square_face(&mut model);
    let shell = model.new_empty(NodeKind::Shell);
    model.add(shell, face.into()).unwrap();
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell).unwrap();
    assert_eq!(model.count(NodeKind::Solid), 1);
    assert_eq!(model.count(NodeKind::Shell), 1);
    assert_eq!(model.count(NodeKind::Face), 1);
    assert_eq!(model[id_of::<Solid>(solid)].shells.len(), 1);
    assert_eq!(model[id_of::<Shell>(shell)].faces.len(), 1);
}

#[test]
fn compound_and_comp_solid_with_children() {
    let mut model = Model::new();
    let a = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
    let b = Cuboid::from_corners(Point::new(2.0, 0.0, 0.0), Point::new(3.0, 1.0, 1.0)).build(&mut model).unwrap();
    let c = compound(&mut model, &[a, b]);
    assert_eq!(model.children(c).count(), 2);
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add_all(cs, &[a, b]).unwrap();
    assert_eq!(model.children(cs).count(), 2);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn two_box_faces_in_one_model() {
    let (source, solid) = box_model();
    let faces: Vec<Shape> = source.explore(solid, NodeKind::Face).take(2).collect();
    let mut model = Model::new();
    for &f in &faces {
        model.import(&source, f).unwrap();
    }
    assert_eq!(model.count(NodeKind::Face), 2);
}

// ---- removal ----

#[test]
fn remove_vertex_marks_it_removed() {
    let mut model = Model::new();
    let v = model.add_vertex(Point::ORIGIN, 0.001);
    model.remove_node(v);
    assert!(model.get(v).is_none());
    // Removing again, or removing a nonexistent node, does nothing.
    model.remove_node(v);
    model.remove_node(NodeId { kind: NodeKind::Face, index: 999 });
}

#[test]
fn remove_face_from_box_keeps_the_rest() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    model.remove_node(face);
    assert!(model.get(face).is_none());
    assert_eq!(model.count(NodeKind::Face), 5);
    assert_eq!(model.count(NodeKind::Edge), 12);
    let live: Vec<FaceId> = model.table::<Face>().iter().map(|(f, _)| f).collect();
    assert!(!live.contains(&face));
}

#[test]
fn remove_subtree_of_lone_face_removes_everything_below() {
    let mut model = Model::new();
    let (v, e, wire, face) = square_face(&mut model);
    model.remove_subtree(face);
    assert!(model.get(face).is_none());
    assert!(model.get(wire).is_none());
    assert!(e.iter().all(|&e| model.get(e).is_none()));
    assert!(v.iter().all(|&v| model.get(v).is_none()));
    assert_eq!(model.count(NodeKind::CoEdge), 0);
}

#[test]
fn remove_subtree_of_solid_removes_shells_and_faces() {
    let (mut model, solid) = box_model();
    model.remove_subtree(solid.node);
    for kind in [NodeKind::Solid, NodeKind::Shell, NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
        assert_eq!(model.count(kind), 0, "{kind:?}");
    }
}

#[test]
fn remove_subtree_of_box_face_keeps_shared_edges_and_vertices() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    model.remove_subtree(face);
    assert!(model.get(face).is_none());
    assert_eq!(model.count(NodeKind::Face), 5);
    assert_eq!(model.count(NodeKind::Edge), 12);
    assert_eq!(model.count(NodeKind::Vertex), 8);
}

// ---- adjacency ----

fn edges_of(model: &Model, face: FaceId) -> HashSet<EdgeId> {
    model.descendants_of::<Edge>(face).into_iter().collect()
}

#[test]
fn box_faces_sharing_edges() {
    let (model, _) = box_model();
    let faces: Vec<FaceId> = model.table::<Face>().iter().map(|(f, _)| f).collect();
    let mut pairs = 0;
    for (i, &a) in faces.iter().enumerate() {
        for &b in &faces[i + 1..] {
            if !edges_of(&model, a).is_disjoint(&edges_of(&model, b)) {
                pairs += 1;
            }
        }
    }
    assert_eq!(pairs, 12);
}

#[test]
fn box_face_edges_and_edge_vertices() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    assert_eq!(edges_of(&model, face).len(), 4);
    let edge = first::<Edge>(&model);
    let vertices: HashSet<VertexId> = model.descendants_of::<Vertex>(edge).into_iter().collect();
    assert_eq!(vertices.len(), 2);
    assert!(model.is_manifold_edge(edge));
}

#[test]
fn box_vertex_has_three_edges_and_edge_has_four_neighbours() {
    let (model, _) = box_model();
    let vertex = first::<Vertex>(&model);
    let edges: HashSet<EdgeId> = model.ancestors_of::<Edge>(vertex).into_iter().collect();
    assert_eq!(edges.len(), 3);
    let edge = first::<Edge>(&model);
    let adjacent: HashSet<EdgeId> = model
        .descendants_of::<Vertex>(edge)
        .into_iter()
        .flat_map(|v| model.ancestors_of::<Edge>(v))
        .filter(|&e| e != edge)
        .collect();
    assert!(adjacent.len() >= 4);
}

#[test]
fn free_edge_has_no_faces() {
    let mut model = Model::new();
    let a = model.add_vertex(Point::ORIGIN, 0.001);
    let b = model.add_vertex(Point::new(1.0, 0.0, 0.0), 0.001);
    let edge = model.add_edge(a, b, 0.001).unwrap();
    assert!(model.ancestors_of::<Face>(edge).is_empty());
}

#[test]
fn queries_on_invalid_ids_are_empty() {
    let (model, _) = box_model();
    assert!(model.ancestors_of::<Edge>(VertexId::from_index(999)).is_empty());
    assert!(model.related(EdgeId::from_index(999).into()).is_empty());
}

// ---- scenarios ----

#[test]
fn split_cylinder_seam_keeps_pcurves_and_faces() {
    let mut model = Model::new();
    Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
    let seam = model.table::<CoEdge>().iter().find(|&(c, _)| model.is_seam(c)).map(|(_, c)| c.edge).unwrap();
    let (first, last) = model.edge_range(seam);
    let vertex = model.add_vertex(Point::new(5.0, 0.0, 7.5), 1.0e-7);
    let (a, b) = model.split_edge(seam, vertex, 0.5 * (first + last)).unwrap();
    for half in [a, b] {
        let coedges: Vec<_> = model.table::<CoEdge>().iter().filter(|(_, c)| c.edge == half).map(|(id, _)| id).collect();
        assert!(!coedges.is_empty());
        for coedge in coedges {
            assert!(model.pcurve(coedge).is_some());
            assert!(model.get(model[coedge].face).is_some());
            let (f, l) = model.coedge_range(coedge);
            assert!(f < l);
        }
    }
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn split_box_edge_retires_old_vertex_uses() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start_use, end_use) = (model[edge].start, model[edge].end);
    let (first, last) = model.edge_range(edge);
    let vertex = model.add_vertex(Point::new(5.0, 10.0, 15.0), 1.0e-7);
    model.split_edge(edge, vertex, 0.5 * (first + last)).unwrap();
    assert!(model.get(start_use).is_none());
    assert!(model.get(end_use).is_none());
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn product_cycle_is_reported() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let a = model.new_assembly();
    let b = model.new_assembly();
    model.add_instance(a, part, Transform::IDENTITY).unwrap();
    model.add_instance(a, b, Transform::IDENTITY).unwrap();
    let in_b = model.add_instance(b, part, Transform::translation(Vector::new(1.0, 0.0, 0.0))).unwrap();
    assert!(validity::check_model(&model).is_valid());
    // Make the instance in `b` point back at `a`: a -> b -> a.
    model[in_b].target = a.into();
    let report = validity::check_model(&model);
    assert!(!report.is_valid());
    assert!(report.has(ProblemKind::ProductStructure, Severity::Error));
}

#[test]
fn shared_part_between_two_assemblies_does_not_copy_topology() {
    let (mut model, solid) = box_model();
    let solids = model.count(NodeKind::Solid);
    let faces = model.count(NodeKind::Face);
    let part = model.new_part(solid).unwrap();
    let a = model.new_assembly();
    let b = model.new_assembly();
    model.add_instance(a, part, Transform::translation(Vector::new(100.0, 0.0, 0.0))).unwrap();
    model.add_instance(b, part, Transform::translation(Vector::new(0.0, 100.0, 0.0))).unwrap();
    assert_eq!(model.count(NodeKind::Solid), solids);
    assert_eq!(model.count(NodeKind::Face), faces);
    assert!(validity::check_model_structure(&model).is_valid());
    assert_eq!(model.children(Shape::from(a)).count(), 1);
    assert_eq!(model.children(Shape::from(b)).count(), 1);
}

#[test]
fn comp_solid_of_three_boxes_each_solid_has_one_parent() {
    let mut model = Model::new();
    let solids: Vec<Shape> = [0.0, 20.0, 40.0]
        .iter()
        .map(|&x| Cuboid::from_corners(Point::new(x, 0.0, 0.0), Point::new(x + 5.0, 5.0, 5.0)).build(&mut model).unwrap())
        .collect();
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add_all(cs, &solids).unwrap();
    for s in solids {
        assert_eq!(model.parents(s.node), vec![cs.node]);
    }
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn compound_with_free_wire_edge_and_vertex_is_valid() {
    let mut model = Model::new();
    let edge = crate::make_edge(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
    let wire = kentron::make::wire(&mut model, &[edge]).unwrap();
    let free_edge = crate::make_edge(&mut model, Point::new(20.0, 0.0, 0.0), Point::new(30.0, 0.0, 0.0));
    let vertex = kentron::make::vertex(&mut model, Point::new(50.0, 50.0, 50.0));
    compound(&mut model, &[wire, free_edge, vertex]);
    assert_eq!(model.count(NodeKind::Compound), 1);
    assert!(model.count(NodeKind::Wire) >= 1);
    assert!(model.count(NodeKind::Edge) >= 2);
    assert!(validity::check_model_structure(&model).is_valid());
}
