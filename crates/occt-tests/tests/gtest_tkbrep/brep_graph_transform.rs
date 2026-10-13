//! BRepGraph_Transform_Test.cxx: transforming whole models and single shapes.
//!
//! `BRepGraph_Transform::Perform(src, dst, t)` is `clone` + `modify::transform_model`;
//! `TransformNode(src, dst, node, t)` is `import` + `modify::transform_in_place`;
//! `MoveRef` is `Model::move_compound_child` / `move_instance`.
//!
//! Not ported: `LocationOnly_NoCopyGeom` (OCCT's shared-geometry mode puts the transform on
//! the root occurrence; in kentron that is `Shape::moved`), runtime mesh-cache tests
//! (`TransformNode_CopyGeomAndMesh_DropsRuntimeCache`), and
//! `PerformIntoNonEmptyTargetDoesNotMoveExistingVertices` (no transform-into-target).

use kentron::math::{Point, Point2, Transform, Transformable, Vector};
use kentron::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};
use kentron::modify::{self, TransformOptions};
use kentron::primitives::Cuboid;
use kentron::properties;
use kentron::tolerance;
use kentron::topology::{
    CoEdgeId, Face, FaceId, Model, NodeKind, Shape, Shell, Solid, Vertex,
};

use crate::graph_support::{box_model, compound, first, id_of};

fn keep_meshes() -> TransformOptions {
    TransformOptions { copy_meshes: true }
}

fn drop_meshes() -> TransformOptions {
    TransformOptions { copy_meshes: false }
}

fn assert_points_near(a: Point, b: Point) {
    assert!(a.distance(b) <= tolerance::LINEAR, "{a:?} != {b:?}");
}

/// Every vertex of `result` is the same-index vertex of `source` moved by `t`.
fn assert_vertices_transformed(source: &Model, result: &Model, t: &Transform) {
    let src: Vec<Point> = source.table::<Vertex>().iter().map(|(_, v)| v.point).collect();
    let res: Vec<Point> = result.table::<Vertex>().iter().map(|(_, v)| v.point).collect();
    assert_eq!(res.len(), src.len());
    for (s, r) in src.into_iter().zip(res) {
        assert_points_near(r, s.transformed(t));
    }
}

/// `node` of `model` imported into a fresh model: the source, and the source transformed.
fn import_and_transform(model: &Model, node: Shape, t: &Transform) -> (Model, Model, Shape) {
    let mut source = Model::new();
    source.import(model, node).unwrap();
    let mut result = Model::new();
    let copy = result.import(model, node).unwrap();
    modify::transform_in_place(&mut result, copy, t, drop_meshes()).unwrap();
    (source, result, copy)
}

fn first_coedge_of_face(model: &Model, face: FaceId) -> CoEdgeId {
    let wire = model[model[face].wires[0]].wire;
    model[wire].coedges[0]
}

#[test]
fn translate_box_face_count() {
    let (model, _) = box_model();
    let mut result = model.clone();
    let t = Transform::translation(Vector::new(100.0, 200.0, 300.0));
    modify::transform_model(&mut result, &t, drop_meshes()).unwrap();
    assert_eq!(result.count(NodeKind::Face), 6);
}

#[test]
fn translate_box_area_preserved() {
    let (model, solid) = box_model();
    let area = properties::area(&model, solid).mass;
    let mut result = model.clone();
    modify::transform_model(&mut result, &Transform::translation(Vector::new(50.0, 0.0, 0.0)), drop_meshes())
        .unwrap();
    let moved = properties::area(&result, solid).mass;
    assert!((moved - area).abs() <= area * 0.01);
}

#[test]
fn translate_box_vertex_points_shifted() {
    let (model, _) = box_model();
    let t = Transform::translation(Vector::new(100.0, 200.0, 300.0));
    let mut result = model.clone();
    modify::transform_model(&mut result, &t, drop_meshes()).unwrap();
    assert_vertices_transformed(&model, &result, &t);
}

#[test]
fn transform_single_face() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let t = Transform::translation(Vector::new(10.0, 20.0, 30.0));
    let (source, result, copy) = import_and_transform(&model, face.into(), &t);
    assert_eq!(copy.kind(), NodeKind::Face);
    assert_eq!(result.count(NodeKind::Face), 1);
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_keeps_mesh_transformed() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let mesh = TriangleMesh::new(
        vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 5.0, 6.0), Point::new(7.0, 8.0, 9.0)],
        vec![[0, 1, 2]],
    )
    .with_deflection(0.1);
    model.set_face_mesh(face.into(), Some(mesh)).unwrap();
    let t = Transform::translation(Vector::new(5.0, 10.0, 15.0));
    let mut result = model.clone();
    modify::transform_model(&mut result, &t, keep_meshes()).unwrap();
    let moved = result.face_mesh(face).expect("mesh dropped");
    assert_eq!(moved.nodes().len(), 3);
    assert_points_near(moved.nodes()[0], Point::new(6.0, 12.0, 18.0));
    assert_points_near(moved.nodes()[2], Point::new(12.0, 18.0, 24.0));
    assert!((moved.deflection() - 0.1).abs() <= tolerance::LINEAR);
    // The source is untouched.
    assert_points_near(model.face_mesh(face).unwrap().nodes()[0], Point::new(1.0, 2.0, 3.0));
}

#[test]
fn transform_without_meshes_drops_them() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let mesh = TriangleMesh::new(
        vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![[0, 1, 2]],
    );
    model.set_face_mesh(face.into(), Some(mesh)).unwrap();
    let mut result = model.clone();
    modify::transform_model(&mut result, &Transform::translation(Vector::new(1.0, 0.0, 0.0)), drop_meshes())
        .unwrap();
    assert_eq!(result.face_mesh(face), None);
}

#[test]
fn move_compound_child_composes_placement() {
    let (mut model, solid) = box_model();
    let c = compound(&mut model, &[solid]);
    let compound = model[id_of::<kentron::topology::Compound>(c)].children.clone();
    assert_eq!(compound.len(), 1);
    let child = compound[0];
    assert!(model[child].placement.is_identity());
    let t = Transform::translation(Vector::new(42.0, 0.0, 0.0));
    model.move_compound_child(child, &t).unwrap();
    assert!((model[child].placement.translation_part().x - 42.0).abs() <= tolerance::LINEAR);
    model.move_compound_child(child, &t).unwrap();
    assert!((model[child].placement.translation_part().x - 84.0).abs() <= tolerance::LINEAR);
}

#[test]
fn move_instance_composes_placement_and_rejects_scaling() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    let child = model[assembly].children[0];
    model.move_instance(instance, &Transform::translation(Vector::new(1.0, 2.0, 3.0))).unwrap();
    assert!((model[child].placement.translation_part().x - 1.0).abs() <= tolerance::LINEAR);
    assert!(model.move_instance(instance, &Transform::scaling(Point::ORIGIN, 2.0)).is_err());
    assert!((model[child].placement.translation_part().x - 1.0).abs() <= tolerance::LINEAR);
}

#[test]
fn transform_node_face_vertex_points_shifted() {
    let (model, _) = box_model();
    let t = Transform::translation(Vector::new(5.0, 10.0, 15.0));
    let (source, result, _) = import_and_transform(&model, first::<Face>(&model).into(), &t);
    assert_eq!(result.count(NodeKind::Face), 1);
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_node_shell() {
    let (model, _) = box_model();
    let t = Transform::translation(Vector::new(100.0, 0.0, 0.0));
    let (source, result, _) = import_and_transform(&model, first::<Shell>(&model).into(), &t);
    assert_eq!(result.count(NodeKind::Shell), 1);
    assert_eq!(result.count(NodeKind::Face), model.count(NodeKind::Face));
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_node_solid() {
    let (model, _) = box_model();
    let t = Transform::translation(Vector::new(0.0, 50.0, 0.0));
    let (source, result, _) = import_and_transform(&model, first::<Solid>(&model).into(), &t);
    assert_eq!(result.count(NodeKind::Solid), 1);
    assert_eq!(result.count(NodeKind::Shell), model.count(NodeKind::Shell));
    assert_eq!(result.count(NodeKind::Face), model.count(NodeKind::Face));
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_node_vertex() {
    let (model, _) = box_model();
    let vertex = first::<Vertex>(&model);
    let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    let (_, result, copy) = import_and_transform(&model, vertex.into(), &t);
    assert_eq!(result.count(NodeKind::Vertex), 1);
    let expected = model[vertex].point + Vector::new(1.0, 2.0, 3.0);
    assert_points_near(result.point(copy), expected);
}

#[test]
fn transform_node_compound() {
    let mut model = Model::new();
    let a = Cuboid::new(5.0, 5.0, 5.0).build(&mut model).unwrap();
    let b = Cuboid::new(3.0, 3.0, 3.0).build(&mut model).unwrap();
    let c = compound(&mut model, &[a, b]);
    let t = Transform::translation(Vector::new(20.0, 0.0, 0.0));
    let (source, result, _) = import_and_transform(&model, c, &t);
    assert_eq!(result.count(NodeKind::Compound), 1);
    assert_eq!(result.count(NodeKind::Solid), model.count(NodeKind::Solid));
    assert_eq!(result.count(NodeKind::Face), model.count(NodeKind::Face));
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_node_comp_solid() {
    let mut model = Model::new();
    let solid = Cuboid::new(4.0, 4.0, 4.0).build(&mut model).unwrap();
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add(cs, solid).unwrap();
    let t = Transform::translation(Vector::new(0.0, 0.0, 7.0));
    let (source, result, _) = import_and_transform(&model, cs, &t);
    assert_eq!(result.count(NodeKind::CompSolid), 1);
    assert_eq!(result.count(NodeKind::Solid), model.count(NodeKind::Solid));
    assert_vertices_transformed(&source, &result, &t);
}

#[test]
fn transform_node_point_mirror_vertex_points_mirrored() {
    let mut model = Model::new();
    let solid = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let t = Transform::mirror_point(Point::ORIGIN);
    let (source, result, copy) = import_and_transform(&model, solid, &t);
    assert_eq!(result.count(NodeKind::Face), model.count(NodeKind::Face));
    assert_eq!(result.count(NodeKind::Vertex), model.count(NodeKind::Vertex));
    assert_vertices_transformed(&source, &result, &t);
    assert!(kentron::validity::check(&result, copy).is_valid());
}

#[test]
fn transform_node_with_meshes_transforms_them() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let coedge = first_coedge_of_face(&model, face);
    let edge = model[coedge].edge;
    let mesh = TriangleMesh::new(
        vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![[0, 1, 2]],
    );
    model.set_face_mesh(face.into(), Some(mesh)).unwrap();
    model
        .set_edge_polyline(edge.into(), Some(Polyline::new(vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 5.0, 6.0)])))
        .unwrap();
    model.set_pcurve_polyline(coedge, Some(Polyline2::new(vec![Point2::new(7.0, 8.0), Point2::new(9.0, 10.0)])))
        .unwrap();
    model.set_coedge_mesh_polyline(coedge, Some(MeshPolyline::new(vec![0, 1]))).unwrap();

    let mut result = Model::new();
    let solid = result.import(&model, first::<Solid>(&model).into()).unwrap();
    let t = Transform::translation(Vector::new(5.0, 0.0, 0.0));
    modify::transform_in_place(&mut result, solid, &t, keep_meshes()).unwrap();

    let face = id_of::<Face>(crate::first(&result, solid, NodeKind::Face));
    let coedge = first_coedge_of_face(&result, face);
    let edge = result[coedge].edge;
    assert_points_near(result.edge_polyline(edge).unwrap().nodes()[0], Point::new(6.0, 2.0, 3.0));
    // Parameter-space polylines do not move.
    assert_eq!(result.pcurve_polyline(coedge).unwrap().nodes()[0], Point2::new(7.0, 8.0));
    assert_points_near(result.face_mesh(face).unwrap().nodes()[1], Point::new(6.0, 0.0, 0.0));
    assert_eq!(result.coedge_mesh_polyline(coedge).unwrap().nodes(), &[0, 1]);
}
