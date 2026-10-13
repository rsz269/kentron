//! BRepTools_GTrsfModification_Test.cxx
//!
//! OCCT tests the per-representation callbacks of its modifier framework; here the same
//! behavior is checked on [`AffineTransformable`] meshes and on
//! [`modify::transform_affine`]. Not ported (OCCT-only concepts):
//! `Curve2dAndVertexParameterPreserveParametrization` and
//! `PolygonOnTriangulation_PreservesDeflection` (pcurves and mesh polylines bound to a
//! face without a wire, vertex parameters on edges); the cached min/max and
//! representation-location parts of the other tests are dropped.
//! `Tolerance_TracksCurrentTransformation` deviates for the shear case (see
//! GTEST_PORTING.md).

use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::geometry::curve::{BezierCurve, Curve};
use kentron::geometry::surface::{BezierSurface, Surface};
use kentron::math::{
    AffineTransform, AffineTransformable, Direction, Matrix3, Point, Transform, Transformable,
    Vector,
};
use kentron::mesh::Polyline;
use kentron::modify;
use kentron::tolerance;
use kentron::topology::{Model, NodeKind};

use crate::modification_support::{
    assert_point_near, diagonal, linear, mesh_only_face, triangle_mesh,
};

#[test]
fn triangulation_reflection_preserves_winding_and_normals() {
    let source = triangle_mesh(true);
    let result = source.affine_transformed(&diagonal(-2.0, 3.0, 4.0));

    assert_eq!(result.triangles()[0], [0, 2, 1], "a reflection reverses the winding");
    let expected_normal = Direction::new(0.0, -1.0 / 3.0, 1.0 / 4.0);
    for n in result.normals().expect("normals") {
        assert!(n.angle(expected_normal) < 1.0e-6);
    }
    assert_eq!(result.deflection(), 1.0);
    assert_eq!(source.deflection(), 0.25);
    assert_point_near(source.nodes()[1], Point::new(1.0, 0.0, 0.0));
}

#[test]
fn triangulation_negative_uniform_scale_preserves_orientation() {
    let t = AffineTransform::from(Transform::scaling(Point::ORIGIN, -2.0));
    let result = triangle_mesh(true).affine_transformed(&t);

    assert_eq!(result.triangles()[0], [0, 2, 1]);
    let normal = result.normals().expect("normals")[0];
    assert!(normal.angle(Direction::new(0.0, 1.0, -1.0)) < 1.0e-6);
    assert_eq!(result.deflection(), 0.5);
}

#[test]
fn triangulation_normal_transformation_is_scale_invariant() {
    let source = triangle_mesh(true);
    let source_normal = source.normals().unwrap()[0];
    for scale in [tolerance::PARAMETRIC, tolerance::INFINITE] {
        let result = source.affine_transformed(&diagonal(scale, scale, scale));
        assert!(result.normals().unwrap()[0].angle(source_normal) < 1.0e-6, "scale {scale}");
    }
}

#[test]
fn triangulation_singular_affinity_handles_normals() {
    let source = triangle_mesh(true);

    let rank_two = source.affine_transformed(&diagonal(2.0, 3.0, 0.0));
    let normal = rank_two.normals().expect("normals kept")[0];
    assert!(normal.angle(Direction::Z) < 1.0e-6);

    let rank_one = source.affine_transformed(&diagonal(2.0, 0.0, 0.0));
    let normal = rank_one.normals().expect("normals kept")[0];
    assert!(normal.angle(source.normals().unwrap()[0]) < 1.0e-6, "collapsed normals are kept");
}

fn bezier_surface(poles: [[(f64, f64, f64); 2]; 2]) -> BezierSurface {
    BezierSurface::new(
        poles.iter().map(|row| row.iter().map(|&(x, y, z)| Point::new(x, y, z)).collect()).collect(),
    )
}

#[test]
fn surface_located_bezier_transforms_evaluated_geometry() {
    let surface = bezier_surface([
        [(0.0, 0.0, 0.0), (0.0, 3.0, 1.0)],
        [(2.0, 0.0, 1.0), (2.0, 3.0, 2.0)],
    ]);
    let placement = Transform::translation(Vector::new(4.0, -2.0, 5.0));
    let mut model = Model::new();
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, Surface::Bezier(surface.clone()), 0.2).unwrap();
    let placed = face.with_placement(placement);

    let g = AffineTransform::new(
        Matrix3::from_rows([[2.0, 0.5, 0.0], [0.0, 3.0, 0.25], [0.0, 0.0, 4.0]]),
        Vector::new(7.0, 8.0, 9.0),
    );
    let result = modify::transform_affine(&mut model, placed, &g).unwrap();
    assert!(result.placement.is_identity());

    let expected = g.transform_point(surface.point(0.35, 0.65).transformed(&placement));
    let actual = model.surface(result).expect("surface").point(0.35, 0.65);
    assert_point_near(actual, expected);
    assert_point_near(surface.pole(1, 1), Point::new(2.0, 3.0, 2.0));
    // Deviation (GTEST_PORTING.md): the tolerance grows by the largest stretch, not by
    // the largest diagonal entry (OCCT's 0.8).
    assert!((model.tolerance(result) - 0.2 * g.max_stretch()).abs() < 1.0e-14);
}

#[test]
fn curve_located_bezier_preserves_range_and_source() {
    let curve = BezierCurve::new(vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
        Point::new(3.0, 2.0, 1.0),
    ]);
    let placement = Transform::translation(Vector::new(-3.0, 4.0, 2.0));
    let mut model = Model::new();
    let edge = model.new_empty(NodeKind::Edge);
    model.set_edge_curve(edge, Curve::Bezier(curve.clone()), 0.1).unwrap();
    model.set_edge_range(edge, 0.2, 0.8).unwrap();
    let placed = edge.with_placement(placement);

    let g = AffineTransform::new(
        Matrix3::from_rows([[1.0, 0.5, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]),
        Vector::new(2.0, -1.0, 5.0),
    );
    let result = modify::transform_affine(&mut model, placed, &g).unwrap();
    assert!(result.placement.is_identity());

    let (result_curve, first, last) = model.edge_curve(result).expect("curve");
    assert_eq!(first, 0.2);
    assert_eq!(last, 0.8);
    let expected = g.transform_point(curve.point(0.4).transformed(&placement));
    assert_point_near(result_curve.point(0.4), expected);
    assert_point_near(curve.poles()[1], Point::new(1.0, 2.0, 0.0));
}

/// Adapted: OCCT also gives the polyline its own "representation location", which kentron
/// does not have; the edge's placement is kept.
#[test]
fn polygon3d_representation_location_transforms_to_edge_frame() {
    let nodes = vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 5.0, 6.0)];
    let polyline = Polyline::new(nodes.clone()).with_parameters(vec![0.25, 0.75]).with_deflection(0.5);
    let placement = Transform::translation(Vector::new(0.0, 20.0, 0.0));
    let mut model = Model::new();
    let edge = model.new_empty(NodeKind::Edge);
    model.set_edge_polyline(edge, Some(polyline.clone())).unwrap();
    let placed = edge.with_placement(placement);

    let g = diagonal(2.0, 3.0, 4.0);
    let result = modify::transform_affine(&mut model, placed, &g).unwrap();
    let result_polyline = model.edge_polyline(result).expect("polyline");

    for (i, node) in result_polyline.nodes().iter().enumerate() {
        let expected = g.transform_point(nodes[i].transformed(&placement));
        assert_point_near(node.transformed(&result.placement), expected);
    }
    assert_eq!(result_polyline.parameters(), Some(&[0.25, 0.75][..]));
    assert_eq!(result_polyline.deflection(), 2.0);
    assert_eq!(polyline.deflection(), 0.5);
}

#[test]
fn missing_polygon_representations_return_false() {
    let mut model = Model::new();
    let edge = model.new_empty(NodeKind::Edge);
    let result = modify::transform_affine(&mut model, edge, &AffineTransform::IDENTITY).unwrap();
    assert!(model.edge_polyline(result).is_none());
}

#[test]
fn modifier_located_surface_and_mesh_stay_coincident() {
    let surface = bezier_surface([
        [(0.0, 0.0, 0.0), (0.0, 1.0, 1.0)],
        [(1.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
    ]);
    let mesh = triangle_mesh(false);
    let placement = Transform::translation(Vector::new(10.0, 20.0, 30.0));
    let mut model = Model::new();
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, Surface::Bezier(surface.clone()), 0.1).unwrap();
    model.set_face_mesh(face, Some(mesh.clone())).unwrap();
    let placed = face.with_placement(placement);

    let g = AffineTransform::new(
        Matrix3::from_rows([[2.0, 0.5, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]),
        Vector::new(5.0, 6.0, 7.0),
    );
    let result = modify::transform_affine(&mut model, placed, &g).unwrap();

    let expected = g.transform_point(surface.point(0.25, 0.75).transformed(&placement));
    let actual = model.surface(result).expect("surface").point(0.25, 0.75);
    assert_point_near(actual.transformed(&result.placement), expected);

    let result_mesh = model.face_mesh(result).expect("mesh");
    for (node, source) in result_mesh.nodes().iter().zip(mesh.nodes()) {
        let expected = g.transform_point(source.transformed(&placement));
        assert_point_near(node.transformed(&result.placement), expected);
    }
}

#[test]
fn triangulation_located_face_uses_local_coordinates() {
    let mesh = triangle_mesh(false);
    let placement = Transform::translation(Vector::new(10.0, 20.0, 30.0));
    let mut model = Model::new();
    let face = mesh_only_face(&mut model, mesh.clone()).with_placement(placement);

    let g = diagonal(2.0, 3.0, 4.0);
    let result = modify::transform_affine(&mut model, face, &g).unwrap();
    let result_mesh = model.face_mesh(result).expect("mesh");
    for (node, source) in result_mesh.nodes().iter().zip(mesh.nodes()) {
        let expected = g.transform_point(source.transformed(&placement));
        assert_point_near(node.transformed(&result.placement), expected);
    }
}

#[test]
fn tolerance_tracks_current_transformation() {
    let mut model = Model::new();
    let vertex = model.new_empty(NodeKind::Vertex);
    model.set_vertex_point(vertex, Point::new(1.0, 2.0, 3.0), 0.1).unwrap();

    let mut check = |g: &AffineTransform, tolerance: f64, point: Option<Point>| {
        let result = modify::transform_affine(&mut model, vertex, g).unwrap();
        assert!((model.tolerance(result) - tolerance).abs() < 1.0e-14, "tolerance");
        if let Some(p) = point {
            assert_point_near(model.point(result), p);
        }
    };

    check(&diagonal(2.0, 3.0, 4.0), 0.4, Some(Point::new(2.0, 6.0, 12.0)));
    // Deviation: OCCT keeps 0.1 under this shear; the tolerance zone is stretched by the
    // largest singular value (the golden ratio).
    let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
    check(&linear([[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), 0.1 * golden, Some(Point::new(3.0, 2.0, 3.0)));
    let c = 0.5_f64.sqrt();
    check(&linear([[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]]), 0.1, None);
    let scale = AffineTransform::from(Transform::scaling(Point::ORIGIN, 3.0));
    check(&scale, 0.3, Some(Point::new(3.0, 6.0, 9.0)));
}
