//! BRepGraph_Tool_Test.cxx and BRepGraph_Convenience_Test.cxx: per-entity queries.
//!
//! Not ported: `Face_Bounds_CylinderFaceReturnsSurfaceBounds` (OCCT's face "bounds" are its
//! surface's domain; kentron's come from the face boundary); `Face_CachedTriangulation_*`
//! (runtime mesh caches, to be ported with BRepGraph_CacheMesh_Test);
//! `ShapesView_RemoveShape_RemovesFoundNode` (lookup of imported TopoDS shapes); the
//! node-id factory checks (Rust ids need no factories).

use kentron::mesh::{MeshPolyline, TriangleMesh};
use kentron::primitives;
use kentron::tolerance;
use kentron::topology::{
    CoEdge, Edge, Face, FaceId, Model, NodeKind, Orientation, Shell, Vertex, Wire,
};
use kentron::math::Point;

use crate::graph_support::{box_model, first};

fn cylinder_model() -> Model {
    let mut model = Model::new();
    primitives::Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
    model
}

#[test]
fn face_nb_wires_box_face_has_one_wire() {
    let (model, _) = box_model();
    for (face, record) in model.table::<Face>().iter() {
        assert_eq!(record.wires.len(), 1, "{face:?}");
    }
}

#[test]
fn face_bounds_box_face_has_finite_bounds() {
    let (model, _) = box_model();
    let (u_min, u_max, v_min, v_max) = model.face_bounds(first::<Face>(&model));
    assert!(u_min < u_max && v_min < v_max);
    assert!(u_min < tolerance::INFINITE && v_min < tolerance::INFINITE);
}

#[test]
fn wire_face_of_returns_valid_face() {
    let (model, _) = box_model();
    assert!(model.face_of_wire(first::<Wire>(&model)).is_some());
}

#[test]
fn wire_is_outer_first_wire_of_box_face_is_outer() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let outer = model.outer_wire(face).expect("outer wire");
    let first_use = model[face].wires[0];
    assert_eq!(model[first_use].wire, outer, "the first wire is the outer one");
}

#[test]
fn edge_nb_faces_box_edge_has_exactly_two_faces() {
    let (model, _) = box_model();
    for (edge, _) in model.table::<Edge>().iter() {
        assert_eq!(model.ancestors_of::<Face>(edge).len(), 2, "{edge:?}");
    }
}

/// Manifold: exactly two faces; boundary: exactly one.
#[test]
fn edge_is_manifold_box_edges_are_manifold() {
    let (model, _) = box_model();
    for (edge, _) in model.table::<Edge>().iter() {
        assert_ne!(model.ancestors_of::<Face>(edge).len(), 1, "{edge:?} is a boundary edge");
    }
}

#[test]
fn edge_is_closed_box_edges_are_open() {
    let (model, _) = box_model();
    assert!(model.count(NodeKind::Edge) > 0);
    for (edge, _) in model.table::<Edge>().iter() {
        assert!(!model.is_closed_edge(edge), "{edge:?}");
    }
}

#[test]
fn edge_is_closed_derived_from_vertex_topology() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, end) = model.edge_end_vertices(edge);
    assert_ne!(start, end);
    assert!(!model.is_closed_edge(edge));
}

fn single_triangle() -> TriangleMesh {
    TriangleMesh::new(vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)], vec![[0, 1, 2]])
}

#[test]
fn co_edge_polygon_on_triangulation_round_trip() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let face = model[coedge].face;
    assert!(model.coedge_mesh_polyline(coedge).is_none());

    model.set_face_mesh(face.into(), Some(single_triangle())).unwrap();
    let polyline = MeshPolyline::new(vec![0, 1]);
    model.set_coedge_mesh_polyline(coedge, Some(polyline.clone())).unwrap();
    assert_eq!(model.coedge_mesh_polyline(coedge), Some(&polyline));
}

#[test]
fn edge_polygon_on_triangulation_resolves_via_face() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let (edge, face) = (model[coedge].edge, model[coedge].face);
    let found = model.find_coedge(edge, face).expect("coedge of edge on face");
    assert!(model.coedge_mesh_polyline(found).is_none());

    model.set_face_mesh(face.into(), Some(single_triangle())).unwrap();
    model.set_coedge_mesh_polyline(coedge, Some(MeshPolyline::new(vec![0, 1]))).unwrap();
    let found = model.find_coedge(edge, face).unwrap();
    assert!(model.coedge_mesh_polyline(found).is_some());
}

#[test]
fn vertex_nb_edges_box_vertex_has_three_edges() {
    let (model, _) = box_model();
    for (vertex, _) in model.table::<Vertex>().iter() {
        assert_eq!(model.ancestors_of::<Edge>(vertex).len(), 3, "{vertex:?}");
    }
}

#[test]
fn shell_is_closed_box_shell_is_closed() {
    let (model, _) = box_model();
    assert!(model.is_closed_shell(first::<Shell>(&model)));
}

#[test]
fn shell_nb_faces_box_shell_has_six_faces() {
    let (model, _) = box_model();
    assert_eq!(model[first::<Shell>(&model)].faces.len(), 6);
}

// Convenience_Test

#[test]
fn edge_def_start_end_differ_for_non_closed() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    if !model.is_closed_edge(edge) {
        let record = &model[edge];
        assert_ne!(model[record.start].vertex, model[record.end].vertex);
    }
}

#[test]
fn edge_ops_find_by_vertices_finds_directed_edge() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, end) = model.edge_end_vertices(edge);
    assert_eq!(model.find_edge(start, end, false), Some(edge));
}

#[test]
fn edge_ops_find_by_vertices_reverse_requires_explicit_flag() {
    let (model, _) = box_model();
    let edge = first::<Edge>(&model);
    let (start, end) = model.edge_end_vertices(edge);
    assert_ne!(start, end);
    assert_eq!(model.find_edge(end, start, false), None);
    assert_eq!(model.find_edge(end, start, true), Some(edge));
}

#[test]
fn face_surface_all_box_faces() {
    let (model, _) = box_model();
    for (face, record) in model.table::<Face>().iter() {
        assert!(model.get(record.surface).is_some(), "{face:?}");
    }
}

#[test]
fn find_p_curve_co_edge_id_valid_pair() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let edge = model[model[model[face].wires[0]].wire].coedges[0];
    let edge = model[edge].edge;
    let coedge = model.find_coedge(edge, face).expect("coedge");
    assert!(model.get(model[coedge].pcurve).is_some(), "the coedge has a pcurve");
}

#[test]
fn find_p_curve_co_edge_id_invalid_pair_returns_null() {
    let (model, _) = box_model();
    assert_eq!(model.find_coedge(first::<Edge>(&model), FaceId::from_index(9999)), None);
}

#[test]
fn shell_face_refs_all_valid() {
    let (model, _) = box_model();
    for &face_use in &model[first::<Shell>(&model)].faces {
        assert!(model.get(model[face_use].face).is_some());
    }
}

/// The seam of a cylinder has two coedges on the lateral face, one per orientation.
#[test]
fn find_p_curve_co_edge_id_with_orientation_seam_edge() {
    let model = cylinder_model();
    let seam = model
        .table::<CoEdge>()
        .iter()
        .find(|(id, _)| model.is_seam(*id))
        .map(|(_, c)| (c.edge, c.face))
        .expect("a cylinder has a seam");
    let forward = model.find_coedge_oriented(seam.0, seam.1, Orientation::Forward).expect("forward");
    let reversed = model.find_coedge_oriented(seam.0, seam.1, Orientation::Reversed).expect("reversed");
    assert_ne!(forward, reversed);
}
