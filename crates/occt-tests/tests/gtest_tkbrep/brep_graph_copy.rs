//! BRepGraph_Copy_Test.cxx: copying whole models and single shapes.
//!
//! A whole-graph copy (`BRepGraph_Copy::Perform`) is `Model::clone` in kentron: geometry is
//! held by value, so a clone never shares it (no `GeomPolicy::Share`), and there are no
//! runtime caches or cache services to drop. Copying one node into another graph
//! (`CopyNode`) is `Model::import`.
//!
//! Not ported: `CopyGraph_{Drops,CanCopyFresh}Runtime*CacheMesh`,
//! `CopyBox_DoesNotCopyTransientCacheServices`, `CopyBox_SharedGeometry`,
//! `CopyGraph_PreservesSupplementAttachmentsAndUids` (internal vertices kept in a side
//! layer), `CopyGraph_NonEmptyTargetUsesSetterPathForPersistentMesh` (a copy into a
//! non-empty graph is an import).

use kentron::boolean;
use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Point, Point2};
use kentron::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};
use kentron::primitives::{Cuboid, Cylinder};
use kentron::properties;
use kentron::topology::{
    CoEdgeId, Face, FaceId, Id, Model, NodeId, NodeKind, Shape,
};

use crate::graph_support::{box_model, compound, first, id_of};

fn total_face_area(model: &Model) -> f64 {
    model
        .table::<Face>()
        .iter()
        .map(|(face, _)| properties::area(model, face.into()).mass.abs())
        .sum()
}

/// The first coedge of the face's first wire.
fn first_coedge_of_face(model: &Model, face: FaceId) -> CoEdgeId {
    let wire = model[model[face].wires[0]].wire;
    model[wire].coedges[0]
}

fn mesh_samples() -> (TriangleMesh, Polyline, Polyline2, MeshPolyline) {
    (
        TriangleMesh::new(
            vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 5.0, 6.0), Point::new(7.0, 8.0, 9.0)],
            vec![[0, 1, 2]],
        ),
        Polyline::new(vec![Point::new(10.0, 11.0, 12.0), Point::new(13.0, 14.0, 15.0)]),
        Polyline2::new(vec![Point2::new(16.0, 17.0), Point2::new(18.0, 19.0)]),
        MeshPolyline::new(vec![1, 2]),
    )
}

/// Attaches the sample meshes to `face`, its first coedge and that coedge's edge.
fn attach_meshes(model: &mut Model, face: FaceId) {
    let (mesh, polyline, polyline2, mesh_polyline) = mesh_samples();
    let coedge = first_coedge_of_face(model, face);
    let edge = model[coedge].edge;
    model.set_face_mesh(face.into(), Some(mesh)).unwrap();
    model.set_edge_polyline(edge.into(), Some(polyline)).unwrap();
    model.set_pcurve_polyline(coedge, Some(polyline2)).unwrap();
    model.set_coedge_mesh_polyline(coedge, Some(mesh_polyline)).unwrap();
}

fn assert_meshes_copied(model: &Model, face: FaceId) {
    let (mesh, polyline, polyline2, mesh_polyline) = mesh_samples();
    let coedge = first_coedge_of_face(model, face);
    let edge = model[coedge].edge;
    assert_eq!(model.face_mesh(face), Some(&mesh));
    assert_eq!(model.edge_polyline(edge), Some(&polyline));
    assert_eq!(model.pcurve_polyline(coedge), Some(&polyline2));
    assert_eq!(model.coedge_mesh_polyline(coedge), Some(&mesh_polyline));
}

#[test]
fn copy_box_face_count() {
    let (model, solid) = box_model();
    let copy = model.clone();
    assert_eq!(copy.count(NodeKind::Face), 6);
    assert_eq!(copy.explore(solid, NodeKind::Face).count(), 6);
}

#[test]
fn copy_keeps_removed_instances_removed() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, kentron::math::Transform::IDENTITY).unwrap();
    model.remove_instance(instance).unwrap();
    let copy = model.clone();
    assert!(copy.get(instance).is_none());
    assert!(copy[assembly].children.is_empty());
    assert_eq!(copy.count(NodeKind::Instance), 0);
}

#[test]
fn copy_box_area_preserved() {
    let (model, solid) = box_model();
    let area = properties::area(&model, solid).mass;
    let copy = model.clone();
    let copied = total_face_area(&copy);
    assert!((copied - area).abs() <= area * 0.01);
}

#[test]
fn copy_box_geometry_is_independent() {
    let (model, _) = box_model();
    let mut copy = model.clone();
    let face = first::<Face>(&copy);
    let moved = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 100.0), Direction::Z));
    copy.set_face_surface(face.into(), moved.clone(), 1.0e-7).unwrap();
    assert_eq!(copy.surface(face), Some(&moved));
    assert_ne!(model.surface(face), Some(&moved));
}

#[test]
fn import_face_copies_meshes() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    attach_meshes(&mut model, face);
    let mut target = Model::new();
    let copy = target.import(&model, face.into()).unwrap();
    assert_meshes_copied(&target, id_of(copy));
}

#[test]
fn clone_copies_meshes() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    attach_meshes(&mut model, face);
    let copy = model.clone();
    assert_meshes_copied(&copy, face);
}

#[test]
fn import_does_not_revive_removed_meshes() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    attach_meshes(&mut model, face);
    let coedge = first_coedge_of_face(&model, face);
    let edge = model[coedge].edge;
    model.set_face_mesh(face.into(), None).unwrap();
    model.set_edge_polyline(edge.into(), None).unwrap();
    let mut target = Model::new();
    let copy: FaceId = id_of(target.import(&model, face.into()).unwrap());
    let copy_edge = target[first_coedge_of_face(&target, copy)].edge;
    assert_eq!(target.face_mesh(copy), None);
    assert_eq!(target.edge_polyline(copy_edge), None);
}

#[test]
fn copy_keeps_history() {
    let (mut model, _) = box_model();
    let original: NodeId = Id::<Face>::from_index(0).into();
    let replacement: NodeId = Id::<Face>::from_index(1).into();
    model.enable_history().record("CopyHistory", original, &[replacement]);
    let copy = model.clone();
    assert_eq!(copy.history().unwrap().modified(original), Some(&[replacement][..]));
}

#[test]
fn copy_keeps_deleted_history() {
    let (mut model, _) = box_model();
    let face: NodeId = first::<Face>(&model).into();
    model.remove_node(face);
    model.enable_history().record_deleted("Delete", &[face]);
    let copy = model.clone();
    assert!(copy.get(Id::<Face>::from_index(0)).is_none());
    assert!(copy.history().unwrap().is_deleted(face));
    assert_eq!(copy.history().unwrap().len(), model.history().unwrap().len());
}

#[test]
fn copy_cylinder_face_count() {
    let mut model = Model::new();
    Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
    assert_eq!(model.clone().count(NodeKind::Face), model.count(NodeKind::Face));
}

#[test]
fn import_single_face() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let mut target = Model::new();
    let copy = target.import(&model, face.into()).unwrap();
    assert_eq!(copy.kind(), NodeKind::Face);
    assert_eq!(target.count(NodeKind::Face), 1);
    assert!(target.count(NodeKind::Vertex) > 0);
    assert!(target.count(NodeKind::Edge) > 0);
    assert!(target.count(NodeKind::Wire) > 0);
    assert_eq!(target.count(NodeKind::Shell), 0);
    assert_eq!(target.count(NodeKind::Solid), 0);
    assert!(properties::area(&target, copy).mass.abs() > 1.0);
}

#[test]
fn copy_faces_only_compound() {
    let (source, solid) = box_model();
    let mut model = Model::new();
    let faces: Vec<Shape> = source
        .explore(solid, NodeKind::Face)
        .map(|f| model.import(&source, f).unwrap())
        .collect();
    compound(&mut model, &faces);
    assert_eq!(model.count(NodeKind::Face), 6);
    assert_eq!(model.count(NodeKind::Solid), 0);
    assert_eq!(model.count(NodeKind::Shell), 0);
    let copy = model.clone();
    assert_eq!(copy.count(NodeKind::Face), 6);
    assert_eq!(copy.count(NodeKind::Solid), 0);
    assert_eq!(copy.count(NodeKind::Shell), 0);
}

#[test]
fn copy_box_same_parameter_preserved() {
    let (model, _) = box_model();
    let copy = model.clone();
    for (coedge, _) in copy.table::<kentron::topology::CoEdge>().iter() {
        assert!(copy.is_same_parameter(coedge), "{coedge:?} lost same-parameter");
        assert!(copy.is_same_range(coedge), "{coedge:?} lost same-range");
    }
}

#[test]
fn fused_boxes_area_preserved() {
    let mut model = Model::new();
    let a = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let b = Cuboid::from_corners(Point::new(10.0, 0.0, 0.0), Point::new(20.0, 10.0, 10.0))
        .build(&mut model)
        .unwrap();
    let fused = boolean::union(&mut model, a, b).unwrap();
    let area = properties::area(&model, fused).mass;
    let mut target = Model::new();
    let copy = target.import(&model, fused).unwrap();
    assert_eq!(target.explore(copy, NodeKind::Face).count(), model.explore(fused, NodeKind::Face).count());
    assert_eq!(target.explore(copy, NodeKind::Edge).count(), model.explore(fused, NodeKind::Edge).count());
    let copied = properties::area(&target, copy).mass;
    assert!((copied - area).abs() <= area * 0.01);
}

#[test]
fn copy_box_stable_ids_preserved() {
    let (model, _) = box_model();
    let copy = model.clone();
    for kind in [
        NodeKind::Vertex,
        NodeKind::Edge,
        NodeKind::Wire,
        NodeKind::Face,
        NodeKind::Shell,
        NodeKind::Solid,
    ] {
        for index in 0..model.slot_count(kind) as u32 {
            let node = NodeId { kind, index };
            assert_eq!(copy.node_stable_id(node), model.node_stable_id(node));
        }
    }
    let face = first::<Face>(&model);
    assert_eq!(copy.stable_id(copy[face].surface), model.stable_id(model[face].surface));
}
