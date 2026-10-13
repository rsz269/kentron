//! BRepAlgoAPI_Cut_Test.cxx

use std::f64::consts::PI;

use crate::*;
use kentron::geometry::convert::surface_to_bspline;
use kentron::geometry::surface::{BezierSurface, Surface};
use kentron::topology::ShellId;
use kentron::{make, primitives, validity};

#[test]
fn modalg_bug_24706_cut_result_shell_is_closed() {
    let mut model = Model::new();
    let box_ = create_box(&mut model, pnt(0.0, 0.0, 0.0), 10.0, 10.0, 10.0);
    let sphere = create_sphere(&mut model, pnt(0.0, 0.0, 0.0), 2.0);
    let result = perform_difference(&mut model, box_, sphere);
    assert!(validity::check(&model, result).is_valid());
    assert_eq!(model.explore(result, NodeKind::Solid).count(), 1);
    assert_eq!(model.explore(result, NodeKind::Shell).count(), 1);
    let shell = model.explore(result, NodeKind::Shell).next().expect("a shell");
    assert!(model.is_closed_shell(ShellId::try_from(shell.node).unwrap()));
}

#[test]
fn modalg_bug_25242_cut_fused_disjoint_boxes_keeps_topology() {
    let mut model = Model::new();
    let box_a = create_box(&mut model, pnt(0.0, 0.0, 0.0), 10.0, 10.0, 10.0);
    let box_b = create_box(&mut model, pnt(20.0, 0.0, 0.0), 10.0, 10.0, 10.0);
    let fused = perform_union(&mut model, box_a, box_b);
    let result = perform_difference(&mut model, fused, box_a);
    assert_eq!(result.kind(), NodeKind::Compound);
    let unique = |kind| model.unique_subshapes(result, kind).len();
    assert_eq!(unique(NodeKind::Vertex), 8);
    assert_eq!(unique(NodeKind::Edge), 12);
    assert_eq!(unique(NodeKind::Wire), 6);
    assert_eq!(unique(NodeKind::Face), 6);
    assert_eq!(unique(NodeKind::Shell), 1);
    assert_eq!(unique(NodeKind::Solid), 1);
    assert_near!(get_volume(&model, result), 1000.0, kentron::tolerance::LINEAR);
}

#[test]
fn nurbs_box_minus_rotated_square_box_f3() {
    let mut model = Model::new();
    let mut nurbs_box = create_box(&mut model, pnt(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    nurbs_box = convert_to_nurbs(&mut model, nurbs_box);
    let r = f64::sqrt(2.0) / 2.0;
    let mut square_box = create_box(&mut model, pnt(0.0, 0.0, 0.0), r, r, 1.0);
    square_box = rotate_shape(square_box, Axis::new(pnt(0.0, 0.0, 0.0), Direction::Z), 45.0 * PI / 180.0);
    let result = perform_difference(&mut model, nurbs_box, square_box);
    validate_result(&model, result, 5.91421, -1.0, false);
}

#[test]
fn rotated_square_box_minus_nurbs_box_f4() {
    let mut model = Model::new();
    let mut nurbs_box = create_box(&mut model, pnt(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    nurbs_box = convert_to_nurbs(&mut model, nurbs_box);
    let r = f64::sqrt(2.0) / 2.0;
    let mut square_box = create_box(&mut model, pnt(0.0, 0.0, 0.0), r, r, 1.0);
    square_box = rotate_shape(square_box, Axis::new(pnt(0.0, 0.0, 0.0), Direction::Z), 45.0 * PI / 180.0);
    let result = perform_difference(&mut model, square_box, nurbs_box);
    validate_result(&model, result, 2.91421, -1.0, false);
}

#[test]
fn box_minus_prism_g8() {
    let mut model = Model::new();
    let box_mb = create_box(&mut model, pnt(-0.5, -0.5, -0.5), 1.0, 1.0, 1.0);
    let face = get_face_by_index(&model, box_mb, 1);
    let semi_inf_prism = create_prism(&mut model, face, vec3(1000.0, 0.0, 0.0));
    let box_ab = create_box(&mut model, pnt(0.0, -1.0, -1.0), 2.0, 2.0, 2.0);
    let result = perform_difference(&mut model, box_ab, semi_inf_prism);
    validate_result(&model, result, 30.0, -1.0, false);
}

#[test]
fn complex_polygon_prism_minus_box_h1() {
    let mut model = Model::new();
    let mut points = Vec::new();
    points.push(pnt(0.0, 0.0, 0.0));
    points.push(pnt(1.0, 0.0, 0.0));
    points.push(pnt(1.0, 3.0, 0.0));
    points.push(pnt(2.0, 3.0, 0.0));
    points.push(pnt(2.0, 0.0, 0.0));
    points.push(pnt(3.0, 0.0, 0.0));
    points.push(pnt(3.0, 5.0, 0.0));
    points.push(pnt(0.0, 5.0, 0.0));
    let wire = create_polygon_wire(&mut model, &points, true);
    let face = create_face_from_wire(&mut model, wire);
    let solid = create_prism(&mut model, face, vec3(0.0, 0.0, 2.0));
    let box_ = create_box(&mut model, pnt(-1.0, 2.0, 1.0), 5.0, 1.0, 3.0);
    let result = perform_difference(&mut model, solid, box_);
    validate_result(&model, result, 68.0, -1.0, false);
}

#[test]
fn complex_polygon_prism_minus_box_h2() {
    let mut model = Model::new();
    let mut points = Vec::new();
    points.push(pnt(0.0, 0.0, 0.0));
    points.push(pnt(1.0, 0.0, 0.0));
    points.push(pnt(1.0, 3.0, 0.0));
    points.push(pnt(2.0, 3.0, 0.0));
    points.push(pnt(2.0, 0.0, 0.0));
    points.push(pnt(3.0, 0.0, 0.0));
    points.push(pnt(3.0, 5.0, 0.0));
    points.push(pnt(0.0, 5.0, 0.0));
    let wire = create_polygon_wire(&mut model, &points, true);
    let face = create_face_from_wire(&mut model, wire);
    let solid = create_prism(&mut model, face, vec3(0.0, 0.0, 2.0));
    let box_ = create_box(&mut model, pnt(-1.0, 2.0, 1.0), 5.0, 1.0, 3.0);
    let result = perform_difference(&mut model, solid, box_);
    validate_result(&model, result, 68.0, -1.0, false);
}

#[test]
fn bug825_sphere_half_space_cut() {
    let mut model = Model::new();
    let size = 50.0;
    let poles = vec![vec![pnt(-size, 0.0, -size), pnt(-size, 0.0, size)], vec![pnt(size, 0.0, -size), pnt(size, 0.0, size)]];
    let bezier = Surface::Bezier(BezierSurface::new(poles));
    let bspline = surface_to_bspline(&bezier).expect("B-spline conversion");
    let face = make::face_on_surface(&mut model, Surface::BSpline(bspline)).expect("face");
    let hsp = primitives::half_space(&mut model, face, pnt(0.0, size, 0.0)).expect("half-space");
    let sph1 = create_sphere(&mut model, pnt(0.0, 0.0, 0.0), 25.0);
    let sph2 = create_sphere(&mut model, pnt(0.0, 0.00001, 0.0), 25.0);
    for sphere in [sph1, sph2] {
        let result = perform_difference(&mut model, sphere, hsp);
        assert_near!(get_surface_area(&model, result), 5890.42, 5.89042);
        assert_eq!(model.explore(result, NodeKind::Face).count(), 2);
        assert_eq!(model.explore(result, NodeKind::Shell).count(), 1);
        assert_eq!(model.explore(result, NodeKind::Solid).count(), 1);
        assert_eq!(result.kind(), NodeKind::Compound);
    }
}
