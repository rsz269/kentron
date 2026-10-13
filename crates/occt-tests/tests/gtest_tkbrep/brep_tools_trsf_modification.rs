//! BRepTools_TrsfModification_Test.cxx
//!
//! Checked through [`modify::transform`]. Not ported: `ParameterlessPolygonOnTriangulationDoesNotThrow`
//! (a mesh polyline bound to a face without a wire). `Polygon3D_DistinctRepresentationLocation`
//! is adapted: kentron has no separate representation locations.

use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Point, Transform, Transformable, Vector};
use kentron::mesh::Polyline;
use kentron::modify::{self, TransformOptions};
use kentron::topology::{Model, NodeKind};

use crate::modification_support::{assert_point_near, mesh_only_face, triangle_mesh};

fn scale(factor: f64) -> Transform {
    Transform::scaling(Point::ORIGIN, factor)
}

const COPY_MESHES: TransformOptions = TransformOptions { copy_meshes: true };

/// A face with exact geometry drops its mesh unless meshes are copied.
#[test]
fn copy_mesh_flag_skips_cached_mesh_on_geometric_face() {
    let mut model = Model::new();
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, Surface::Plane(Plane::XY), 0.1).unwrap();
    model.set_face_mesh(face, Some(triangle_mesh(false))).unwrap();

    let without = modify::transform(&mut model, face, &scale(2.0), TransformOptions::default()).unwrap();
    assert!(model.face_mesh(without).is_none());

    let with = modify::transform(&mut model, face, &scale(2.0), COPY_MESHES).unwrap();
    let mesh = model.face_mesh(with).expect("copied mesh");
    assert_point_near(mesh.nodes()[1], Point::new(2.0, 0.0, 0.0));
    assert_eq!(mesh.deflection(), 0.5);
    assert_point_near(model.face_mesh(face).unwrap().nodes()[1], Point::new(1.0, 0.0, 0.0));
}

#[test]
fn polygon3d_distinct_representation_location() {
    let nodes = vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 5.0, 6.0)];
    let placement = Transform::translation(Vector::new(0.0, 20.0, 0.0));
    let mut model = Model::new();
    let edge = model.new_empty(NodeKind::Edge);
    model.set_edge_polyline(edge, Some(Polyline::new(nodes.clone()))).unwrap();
    let placed = edge.with_placement(placement);

    let result = modify::transform(&mut model, placed, &scale(2.0), COPY_MESHES).unwrap();
    let polyline = model.edge_polyline(result).expect("polyline");
    for (node, source) in polyline.nodes().iter().zip(&nodes) {
        let expected = source.transformed(&placement).transformed(&scale(2.0));
        assert_point_near(node.transformed(&result.placement), expected);
    }
}

#[test]
fn negative_scale_updates_winding_normals_and_deflection() {
    let mut model = Model::new();
    let face = mesh_only_face(&mut model, triangle_mesh(true));

    let result = modify::transform(&mut model, face, &scale(-2.0), COPY_MESHES).unwrap();
    let mesh = model.face_mesh(result).expect("mesh");
    assert_eq!(mesh.triangles()[0], [0, 2, 1]);
    assert!(mesh.normals().expect("normals")[0].angle(Direction::new(0.0, 1.0, -1.0)) < 1.0e-6);
    assert_eq!(mesh.deflection(), 0.5);
}
