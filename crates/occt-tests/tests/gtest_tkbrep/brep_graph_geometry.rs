//! BRepGraph_Geometry_Test.cxx
//!
//! OCCT compares surface *handles* to detect sharing; kentron geometry slots hold values,
//! so "same surface" means equal values. Not ported: `RepId_FactoryMethods`,
//! `RepId_UntypedArithmetic_PreservesKindAndIndex` (id mechanics, see brep_graph_ids) and
//! the derived-state cache part of `DerivedStateCache_LazyAndFreshAfterPCurveRangeMutation`.

use std::collections::HashMap;

use kentron::geometry::curve2::{Curve2, Line2};
use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Direction2, Frame, Point, Point2, Transform, Vector};
use kentron::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};
use kentron::primitives::{self, Cuboid};
use kentron::tolerance;
use kentron::topology::{
    CoEdge, CoEdgeId, Edge, EdgeCurve, Face, FaceId, FaceSurface, Model, NodeId, NodeKind, PCurve,
    Shape, Shell, Solid, Vertex, Wire,
};

use crate::graph_support::{box_model, compound, first, id_of};
use crate::make_box;

fn sphere_model() -> Model {
    let mut model = Model::new();
    primitives::Sphere::new(15.0).build(&mut model).unwrap();
    model
}

fn cylinder_model() -> Model {
    let mut model = Model::new();
    primitives::Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
    model
}

fn coedges(model: &Model) -> Vec<CoEdgeId> {
    model.table::<CoEdge>().iter().map(|(id, _)| id).collect()
}

/// Faces other than `face` whose surface equals its surface.
fn same_domain_faces(model: &Model, face: FaceId) -> Vec<FaceId> {
    let Some(surface) = model.face_surface(face) else { return Vec::new() };
    model
        .table::<Face>()
        .iter()
        .map(|(id, _)| id)
        .filter(|&other| other != face && model.face_surface(other) == Some(surface))
        .collect()
}

#[test]
fn sphere_all_faces_same_surface() {
    let model = sphere_model();
    let first_surface = model.face_surface(first::<Face>(&model)).expect("surface");
    for (face, _) in model.table::<Face>().iter() {
        assert_eq!(model.face_surface(face), Some(first_surface));
    }
}

#[test]
fn sphere_all_faces_share_surface() {
    let model = sphere_model();
    let first_surface = model.face_surface(first::<Face>(&model)).expect("surface");
    let same = model.table::<Face>().iter().filter(|(f, _)| model.face_surface(*f) == Some(first_surface)).count();
    assert_eq!(same, model.count(NodeKind::Face));
}

#[test]
fn box_curve3d_valid_for_all12_edges() {
    let (model, _) = box_model();
    assert_eq!(model.count(NodeKind::Edge), 12);
    for (edge, _) in model.table::<Edge>().iter() {
        assert!(model.edge_curve(edge).is_some(), "{edge:?}");
    }
}

#[test]
fn box_find_p_curve_co_edge_id_all_edge_face_pairs_valid() {
    let (model, _) = box_model();
    let mut count = 0;
    for (_, coedge) in model.table::<CoEdge>().iter() {
        assert!(model.find_coedge(coedge.edge, coedge.face).is_some());
        count += 1;
    }
    assert!(count > 0);
}

#[test]
fn co_edge_face_id_valid() {
    let (model, _) = box_model();
    for (id, coedge) in model.table::<CoEdge>().iter() {
        assert!(model.get(coedge.face).is_some(), "{id:?}");
    }
}

#[test]
fn co_edge_param_range_non_zero() {
    let (model, _) = box_model();
    for coedge in coedges(&model) {
        let (first, last) = model.coedge_range(coedge);
        assert!((last - first).abs() > tolerance::PARAMETRIC, "{coedge:?}");
    }
}

#[test]
fn set_p_curve_two_arg_preserves_existing_range() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let before = model.coedge_range(coedge);
    let replacement = Curve2::Line(Line2::new(Point2::new(100.0, 200.0), Direction2::X));
    model.set_pcurve(coedge, replacement.clone()).unwrap();
    let after = model.coedge_range(coedge);
    assert!((after.0 - before.0).abs() <= tolerance::PARAMETRIC);
    assert!((after.1 - before.1).abs() <= tolerance::PARAMETRIC);
    assert_eq!(model.pcurve(coedge), Some(&replacement));
}

#[test]
fn face_def_surface_is_not_null() {
    let (model, _) = box_model();
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.face_surface(face).is_some(), "{face:?}");
    }
}

#[test]
fn same_domain_faces_simple_box_empty() {
    let (model, _) = box_model();
    for (face, _) in model.table::<Face>().iter() {
        assert!(same_domain_faces(&model, face).is_empty(), "{face:?}");
    }
}

/// Two faces on equal planes have separate surface slots but the same domain.
#[test]
fn same_domain_faces_shared_surface_across_owned_uses() {
    let mut model = Model::new();
    let plane = Surface::Plane(Plane::XY);
    let face_a = model.new_empty(NodeKind::Face);
    let face_b = model.new_empty(NodeKind::Face);
    model.set_face_surface(face_a, plane.clone(), 1.0e-7).unwrap();
    model.set_face_surface(face_b, plane, 1.0e-7).unwrap();
    let (a, b) = (id_of::<Face>(face_a), id_of::<Face>(face_b));

    assert_ne!(model[a].surface, model[b].surface);
    assert_eq!(same_domain_faces(&model, a), [b]);
    assert_eq!(same_domain_faces(&model, b), [a]);
}

#[test]
fn compound_with_moved_child_shared_solid_def() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 10.0, 10.0);
    let moved = solid.moved(&Transform::translation(Vector::new(100.0, 0.0, 0.0)));
    let c = compound(&mut model, &[solid, moved]);

    assert_eq!(model.count(NodeKind::Solid), 1);
    let children = &model[id_of::<kentron::topology::Compound>(c)].children;
    assert_eq!(children.len(), 2);
    assert!(model[children[0]].placement.is_identity());
    assert!(!model[children[1]].placement.is_identity());
}

/// Analytic box faces have no mesh until one is computed.
#[test]
fn face_def_triangulation_null_for_analytic_no_crash() {
    let (model, _) = box_model();
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.face_mesh(face).is_none(), "{face:?}");
    }
}

#[test]
fn def_counts_match_nb() {
    let (model, _) = box_model();
    assert_eq!(model.table::<Solid>().iter().count(), model.count(NodeKind::Solid));
    assert_eq!(model.table::<Shell>().iter().count(), model.count(NodeKind::Shell));
    assert_eq!(model.table::<Face>().iter().count(), model.count(NodeKind::Face));
    assert_eq!(model.table::<Wire>().iter().count(), model.count(NodeKind::Wire));
    assert_eq!(model.table::<Edge>().iter().count(), model.count(NodeKind::Edge));
    assert_eq!(model.table::<Vertex>().iter().count(), model.count(NodeKind::Vertex));
}

#[test]
fn edge_def_all_curves3d_non_null() {
    let (model, _) = box_model();
    for (edge, _) in model.table::<Edge>().iter() {
        if !model.is_degenerate(edge) {
            assert!(model.edge_curve(edge).is_some(), "{edge:?}");
        }
    }
}

#[test]
fn all_co_edges_have_curve2d() {
    let (model, _) = box_model();
    let all = coedges(&model);
    assert!(!all.is_empty());
    for coedge in all {
        assert!(model.pcurve(coedge).is_some(), "{coedge:?}");
    }
}

/// Without a stored pcurve, a coedge on a plane gets the projection of its edge.
#[test]
fn co_edge_p_curve_adaptor_falls_back_on_plane_when_stored_p_curve_removed() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let edge = model[coedge].edge;
    assert!(model.pcurve(coedge).is_some());

    model.clear_pcurve(coedge).unwrap();
    assert!(model.pcurve(coedge).is_none());
    let (_, first_param, last_param) = model.effective_pcurve(coedge).expect("projected pcurve");
    let (edge_first, edge_last) = model.edge_range(edge);
    assert!((first_param - edge_first).abs() <= tolerance::PARAMETRIC);
    assert!((last_param - edge_last).abs() <= tolerance::PARAMETRIC);
}

/// The component of a face: its first solid ancestor, else shell, else itself.
fn component_of(model: &Model, face: FaceId) -> NodeId {
    for kind in [NodeKind::Solid, NodeKind::Shell] {
        if let Some(root) = model.explore_up(face).of_kind(kind).next() {
            return root.node;
        }
    }
    face.into()
}

fn face_counts_by_component(model: &Model) -> HashMap<NodeId, usize> {
    let mut counts = HashMap::new();
    for (face, _) in model.table::<Face>().iter() {
        *counts.entry(component_of(model, face)).or_default() += 1;
    }
    counts
}

fn two_box_compound() -> Model {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 10.0, 10.0);
    let box2 = Cuboid::new(60.0, 60.0, 60.0)
        .with_frame(Frame::new(Point::new(50.0, 50.0, 50.0), Direction::Z, Direction::X))
        .build(&mut model)
        .unwrap();
    compound(&mut model, &[box1, box2]);
    model
}

#[test]
fn connected_components_single_box_one_component() {
    let (model, _) = box_model();
    let counts = face_counts_by_component(&model);
    assert_eq!(counts.len(), 1);
    assert_eq!(counts.values().copied().collect::<Vec<_>>(), [6]);
}

#[test]
fn connected_components_two_box_compound_two_components() {
    assert_eq!(face_counts_by_component(&two_box_compound()).len(), 2);
}

#[test]
fn connected_components_two_box_compound_faces_grouped_per_root() {
    let counts = face_counts_by_component(&two_box_compound());
    assert_eq!(counts.len(), 2);
    assert!(counts.values().all(|&c| c == 6));
}

#[test]
fn connected_components_two_box_compound_cover_all_faces() {
    let model = two_box_compound();
    let total: usize = face_counts_by_component(&model).values().sum();
    assert_eq!(total, model.count(NodeKind::Face));
}

#[test]
fn box_edge_def_same_parameter_is_set() {
    let (model, _) = box_model();
    for coedge in coedges(&model) {
        assert!(model.is_same_parameter(coedge), "{coedge:?}");
    }
}

#[test]
fn box_edge_def_same_range_is_set() {
    let (model, _) = box_model();
    for coedge in coedges(&model) {
        assert!(model.is_same_range(coedge), "{coedge:?}");
    }
}

/// Shifting a coedge's range breaks both the same-range and same-parameter properties.
#[test]
fn derived_state_cache_lazy_and_fresh_after_p_curve_range_mutation() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let edge = model[coedge].edge;
    assert!(model.is_same_range(coedge));

    let (first_param, last_param) = model.edge_range(edge);
    model.set_coedge_range(coedge, first_param + 1.0, last_param).unwrap();
    assert!(!model.is_same_range(coedge));
    assert!(!model.is_same_parameter(coedge));
}

#[test]
fn cylinder_seam_edge_has_two_co_edges() {
    let model = cylinder_model();
    let (coedge, pair) = coedges(&model)
        .into_iter()
        .find_map(|c| model.seam_pair(c).map(|p| (c, p)))
        .expect("a cylinder has a seam");
    assert_ne!(model[coedge].orientation, model[pair].orientation);
    assert_eq!(model[coedge].face, model[pair].face);
}

#[test]
fn cylinder_seam_edge_find_p_curve_co_edge_id_with_orientation() {
    let model = cylinder_model();
    let seam = coedges(&model).into_iter().find(|&c| model.is_seam(c)).expect("seam");
    let (edge, face) = (model[seam].edge, model[seam].face);
    let forward = model.find_coedge_oriented(edge, face, kentron::topology::Orientation::Forward);
    let reversed = model.find_coedge_oriented(edge, face, kentron::topology::Orientation::Reversed);
    assert!(forward.is_some() && reversed.is_some());
    assert_ne!(forward, reversed);
}

#[test]
fn box_find_p_curve_co_edge_id_matches_co_edge() {
    let (model, _) = box_model();
    for (id, coedge) in model.table::<CoEdge>().iter() {
        assert_eq!(model.find_coedge_oriented(coedge.edge, coedge.face, coedge.orientation), Some(id));
    }
}

#[test]
fn box_rep_counts_match_topology() {
    let (model, _) = box_model();
    assert!(model.table::<FaceSurface>().len() > 0);
    assert!(model.table::<EdgeCurve>().len() > 0);
    assert!(model.table::<PCurve>().len() > 0);
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.face_surface(face).is_some());
    }
    for (id, coedge) in model.table::<CoEdge>().iter() {
        if coedge.pcurve.is_valid() {
            assert!(model.pcurve(id).is_some());
        }
    }
}

#[test]
fn cylinder_triangulation_reps_populated() {
    let model = cylinder_model();
    for (face, record) in model.table::<Face>().iter() {
        assert_eq!(model.get(record.mesh).is_some(), model.face_mesh(face).is_some());
    }
}

#[test]
fn compound_two_boxes_surface_dedup() {
    let mut model = Model::new();
    let box1 = make_box(&mut model, 10.0, 10.0, 10.0);
    let box2 = make_box(&mut model, 20.0, 20.0, 20.0);
    compound(&mut model, &[box1, box2]);
    assert_eq!(model.count(NodeKind::Face), 12);
    let surfaces = model.table::<FaceSurface>().len();
    assert!(surfaces > 0 && surfaces <= 12);
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.face_surface(face).is_some());
    }
}

#[test]
fn box_polygon2d_rep_matches_inline() {
    let (model, _) = box_model();
    for (id, coedge) in model.table::<CoEdge>().iter() {
        if model.get(coedge.pcurve_polyline).is_some() {
            assert!(model.pcurve_polyline(id).is_some());
        }
    }
}

/// Clearing a geometry slot keeps the slot id on the entity; setting it again revives the
/// same slot.
#[test]
fn clear_use_owners_marks_records_removed_and_set_reuses_slots() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let face_shape = Shape::from(face);
    let surface_slot = model[face].surface;
    let surface = model.face_surface(face).expect("surface").clone();

    let edge = first::<Edge>(&model);
    let edge_shape = Shape::from(edge);
    let curve_slot = model[edge].curve;
    let (curve, _, _) = model.edge_curve(edge).expect("curve");

    let coedge = first::<CoEdge>(&model);
    let pcurve_slot = model[coedge].pcurve;
    let pcurve = model.pcurve(coedge).expect("pcurve").clone();

    let surfaces = model.table::<FaceSurface>().len();
    model.clear_face_surface(face_shape).unwrap();
    assert_eq!(model[face].surface, surface_slot);
    assert!(model.face_surface(face).is_none());
    assert_eq!(model.table::<FaceSurface>().len(), surfaces - 1);
    model.set_face_surface(face_shape, surface, 1.0e-7).unwrap();
    assert_eq!(model[face].surface, surface_slot);
    assert_eq!(model.table::<FaceSurface>().len(), surfaces);

    let curves = model.table::<EdgeCurve>().len();
    model.clear_edge_curve(edge_shape).unwrap();
    assert_eq!(model[edge].curve, curve_slot);
    assert!(model.edge_curve(edge).is_none());
    assert_eq!(model.table::<EdgeCurve>().len(), curves - 1);
    model.set_edge_curve(edge_shape, curve, 1.0e-7).unwrap();
    assert_eq!(model[edge].curve, curve_slot);
    assert_eq!(model.table::<EdgeCurve>().len(), curves);

    let pcurves = model.table::<PCurve>().len();
    model.clear_pcurve(coedge).unwrap();
    assert_eq!(model[coedge].pcurve, pcurve_slot);
    assert!(model.pcurve(coedge).is_none());
    assert_eq!(model.table::<PCurve>().len(), pcurves - 1);
    model.set_pcurve_with_range(coedge, pcurve, 0.0, 1.0).unwrap();
    assert_eq!(model[coedge].pcurve, pcurve_slot);
    assert_eq!(model.table::<PCurve>().len(), pcurves);

    // meshes and polylines
    let mesh = TriangleMesh::new(
        vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![[0, 1, 2]],
    );
    model.set_face_mesh(face_shape, Some(mesh.clone())).unwrap();
    let mesh_slot = model[face].mesh;
    model.set_face_mesh(face_shape, None).unwrap();
    assert_eq!(model[face].mesh, mesh_slot);
    assert!(model.face_mesh(face).is_none());
    model.set_face_mesh(face_shape, Some(mesh)).unwrap();
    assert_eq!(model[face].mesh, mesh_slot);

    let polyline = Polyline::new(vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)]);
    model.set_edge_polyline(edge_shape, Some(polyline.clone())).unwrap();
    let polyline_slot = model[edge].polyline;
    model.set_edge_polyline(edge_shape, None).unwrap();
    assert_eq!(model[edge].polyline, polyline_slot);
    assert!(model.edge_polyline(edge).is_none());
    model.set_edge_polyline(edge_shape, Some(polyline)).unwrap();
    assert_eq!(model[edge].polyline, polyline_slot);

    let polyline2 = Polyline2::new(vec![Point2::ORIGIN, Point2::new(1.0, 0.0)]);
    model.set_pcurve_polyline(coedge, Some(polyline2.clone())).unwrap();
    let polyline2_slot = model[coedge].pcurve_polyline;
    model.set_pcurve_polyline(coedge, None).unwrap();
    assert_eq!(model[coedge].pcurve_polyline, polyline2_slot);
    model.set_pcurve_polyline(coedge, Some(polyline2)).unwrap();
    assert_eq!(model[coedge].pcurve_polyline, polyline2_slot);

    let mesh_polyline = MeshPolyline::new(vec![0, 1]);
    model.set_coedge_mesh_polyline(coedge, Some(mesh_polyline.clone())).unwrap();
    let mesh_polyline_slot = model[coedge].mesh_polyline;
    model.set_coedge_mesh_polyline(coedge, None).unwrap();
    assert_eq!(model[coedge].mesh_polyline, mesh_polyline_slot);
    assert!(model.coedge_mesh_polyline(coedge).is_none());
    model.set_coedge_mesh_polyline(coedge, Some(mesh_polyline)).unwrap();
    assert_eq!(model[coedge].mesh_polyline, mesh_polyline_slot);
}
