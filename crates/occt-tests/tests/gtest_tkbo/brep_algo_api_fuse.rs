//! BRepAlgoAPI_Fuse_Test.cxx

use std::f64::consts::PI;

use crate::*;

#[test]
fn rotated_sphere_plus_box_d5() {
    let mut model = Model::new();
    let mut sphere = create_sphere(&mut model, pnt(0.0, 0.0, 0.0), 1.0);
    let trsf1 = Transform::rotation(Axis::new(pnt(0.0, 0.0, 0.0), Direction::Z), -PI / 2.0);
    let trsf2 = Transform::rotation(Axis::new(pnt(0.0, 0.0, 0.0), Direction::Y), -PI / 4.0);
    let trsf_combined = trsf2 * trsf1;
    sphere = sphere.moved(&trsf_combined);
    let box_ = create_box(&mut model, pnt(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    let result = perform_union(&mut model, sphere, box_);
    validate_result(&model, result, 14.6393, -1.0, false);
}

#[test]
fn box_plus_rotated_sphere_d6() {
    let mut model = Model::new();
    let mut sphere = create_sphere(&mut model, pnt(0.0, 0.0, 0.0), 1.0);
    let trsf1 = Transform::rotation(Axis::new(pnt(0.0, 0.0, 0.0), Direction::Z), -PI / 2.0);
    let trsf2 = Transform::rotation(Axis::new(pnt(0.0, 0.0, 0.0), Direction::Y), -PI / 4.0);
    let trsf_combined = trsf2 * trsf1;
    sphere = sphere.moved(&trsf_combined);
    let box_ = create_box(&mut model, pnt(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    let result = perform_union(&mut model, box_, sphere);
    validate_result(&model, result, 14.6393, -1.0, false);
}

#[test]
fn profile_based_prisms_d9() {
    let mut model = Model::new();
    let plane1 = Plane::from_frame(Frame::from_z(pnt(0.0, 0.0, 0.0), Direction::Z));
    let profile_steps1 = vec![Step::AlongX(100.0), Step::AlongY(100.0), Step::AlongX(-200.0), Step::AlongY(-200.0), Step::AlongX(100.0)];
    let profile1 = create_profile(&mut model, plane1, &profile_steps1);
    let prism1 = create_prism(&mut model, profile1, vec3(0.0, 0.0, 100.0));
    let plane2 = Plane::from_frame(Frame::from_z(pnt(0.0, 0.0, 100.0), Direction::Z));
    let profile_steps2 = vec![Step::AlongX(-100.0), Step::AlongY(100.0), Step::AlongX(100.0)];
    let profile2 = create_profile(&mut model, plane2, &profile_steps2);
    let prism2 = create_prism(&mut model, profile2, vec3(0.0, 0.0, 100.0));
    let result = perform_union(&mut model, prism1, prism2);
    validate_result(&model, result, 180000.0, -1.0, false);
}

#[test]
fn complex_profile_with_scaling_e1() {
    let mut model = Model::new();
    let scale = 100.0;
    let plane1 = Plane::from_frame(Frame::from_z(pnt(0.0, 0.0, 0.0), Direction::Z));
    let profile_steps1 = vec![Step::Arc { radius: 50.0 * scale, angle: 180.0 }, Step::AlongX(-100.0 * scale), Step::Arc { radius: 50.0 * scale, angle: 180.0 }];
    let profile1 = create_profile(&mut model, plane1, &profile_steps1);
    let prism1 = create_prism(&mut model, profile1, vec3(0.0, 0.0, 30.0 * scale));
    let plane2 = Plane::from_frame(Frame::from_z(pnt(-200.0 * scale, -50.0 * scale, 0.0), Direction::Z));
    let profile_steps2 = vec![Step::AlongX(300.0 * scale), Step::AlongY(200.0 * scale), Step::AlongX(-300.0 * scale), Step::AlongY(-200.0 * scale)];
    let profile2 = create_profile(&mut model, plane2, &profile_steps2);
    let prism2 = create_prism(&mut model, profile2, vec3(0.0, 0.0, -50.0 * scale));
    let result = perform_union(&mut model, prism2, prism1);
    validate_result(&model, result, 1.85425e+09, -1.0, false);
}

#[test]
fn occ277_fuse_and_common_of_overlapping_boxes() {
    // OCCT checks that the results are not null; here a failed operation is an error.
    let mut model = Model::new();
    let shape1 = create_box(&mut model, Point::ORIGIN, 100.0, 100.0, 100.0);
    let shape2 = create_box(&mut model, pnt(50.0, 50.0, 50.0), 200.0, 200.0, 200.0);
    perform_union(&mut model, shape1, shape2);
    perform_intersection(&mut model, shape1, shape2);
}
