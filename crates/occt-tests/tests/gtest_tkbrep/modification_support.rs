//! Helpers shared by the BRepTools_*Modification ports.

use kentron::math::{AffineTransform, Direction, Matrix3, Point, Vector};
use kentron::mesh::TriangleMesh;
use kentron::topology::{Model, NodeKind, Shape};

/// One triangle (0,0,0), (1,0,0), (0,1,1), deflection 0.25, optionally with normals
/// (0,-1,1) at every node.
pub fn triangle_mesh(with_normals: bool) -> TriangleMesh {
    let mesh = TriangleMesh::new(
        vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 1.0)],
        vec![[0, 1, 2]],
    )
    .with_deflection(0.25);
    match with_normals {
        true => mesh.with_normals(vec![Direction::new(0.0, -1.0, 1.0); 3]),
        false => mesh,
    }
}

/// A face carrying only a mesh (no surface).
pub fn mesh_only_face(model: &mut Model, mesh: TriangleMesh) -> Shape {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_mesh(face, Some(mesh)).unwrap();
    face
}

/// The linear map with the given matrix rows.
pub fn linear(rows: [[f64; 3]; 3]) -> AffineTransform {
    AffineTransform::new(Matrix3::from_rows(rows), Vector::ZERO)
}

pub fn diagonal(a: f64, b: f64, c: f64) -> AffineTransform {
    linear([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
}

pub fn assert_point_near(actual: Point, expected: Point) {
    let d = actual.distance(expected);
    assert!(d < 1.0e-12, "{actual:?} is {d} away from {expected:?}");
}
