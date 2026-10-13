//! BRepPrimAPI_MakeRevol_Test.cxx

use std::f64::consts::{FRAC_PI_2, TAU};

use kentron::geometry::curve::{Circle, Curve};
use kentron::geometry::surface::Surface;
use kentron::math::{Axis, Direction, Frame, Point};
use kentron::primitives::Sphere;
use kentron::topology::{FaceId, Model, NodeKind};
use kentron::{make, sweep, validity};

fn check_toroidal_revolution(x_direction: Direction, radius: f64) {
    let mut model = Model::new();
    let frame = Frame::new(Point::new(0.0, 100.00001, 0.0), Direction::X, x_direction);
    let circle = Circle::new(frame, radius);
    let edge = make::edge_from_curve(&mut model, Curve::Circle(circle)).unwrap();

    let result = sweep::revolve(&mut model, edge, Axis::Z, FRAC_PI_2).unwrap().shape;
    assert!(validity::check(&model, result).is_valid());

    let face = model.explore(result, NodeKind::Face).next().expect("no face");
    let surface = model.face_surface(FaceId::try_from(face.node).unwrap()).expect("no surface");
    assert!(matches!(surface, Surface::Torus(_)), "expected a torus, got {surface:?}");
}

/// From tests/bugs/modalg_5/bug23954: revolving circles in the four orientations used by
/// DRAW must give toroidal (not general revolution) surfaces.
#[test]
fn modalg_bug_23954_circle_revolutions_keep_toroidal_surface() {
    check_toroidal_revolution(Direction::Y, 100.0);
    check_toroidal_revolution(Direction::Z, 100.0);
    check_toroidal_revolution(Direction::NEG_Y, 100.0);
    check_toroidal_revolution(Direction::Z, 10.0);
}

/// From tests/bugs/modalg_5/bug25578: revolving a face of a partial sphere about either
/// axis must not fail at the singularity.
#[test]
fn modalg_bug_25578_partial_sphere_revolutions() {
    let mut model = Model::new();
    let sphere = Sphere::new(5.0)
        .with_latitudes(0.0, FRAC_PI_2)
        .with_angle(FRAC_PI_2)
        .build(&mut model)
        .unwrap();
    let face = model.explore(sphere, NodeKind::Face).next().expect("no face");

    sweep::revolve(&mut model, face, Axis::Z, TAU).expect("revolution about Z failed");
    sweep::revolve(&mut model, face, Axis::X, FRAC_PI_2).expect("revolution about X failed");
}
