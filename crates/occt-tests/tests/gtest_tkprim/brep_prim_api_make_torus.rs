//! BRepPrimAPI_MakeTorus_Test.cxx

use kentron::geometry::ParametricSurface;
use std::f64::consts::FRAC_PI_2;

use occt_tests::assert_near;
use kentron::geometry::surface::Surface;
use kentron::primitives::Torus;
use kentron::topology::{FaceId, Model, NodeKind};

#[test]
fn lateral_face_parameterization() {
    let (major, minor) = (5.0, 1.0);
    let mut model = Model::new();
    let shape = Torus::new(major, minor).build(&mut model).unwrap();

    let face = model.explore(shape, NodeKind::Face).next().expect("no faces found");
    let surface = model.face_surface(FaceId::try_from(face.node).unwrap()).expect("no surface");
    let Surface::Torus(torus) = surface else {
        panic!("lateral face is not a toroidal surface: {surface:?}");
    };

    assert_near!(torus.major_radius(), major, 1e-10);
    assert_near!(torus.minor_radius(), minor, 1e-10);
    assert_near!(torus.point(0.0, FRAC_PI_2).z, minor, 1e-10, "top of the tube at v=PI/2");
    assert_near!(torus.point(0.0, -FRAC_PI_2).z, -minor, 1e-10, "bottom of the tube at v=-PI/2");
}

// Equal end heights. OCCT asserts these partial tori are *invalid*, documenting a known
// limitation of its implementation (a degenerate edge on the axis). The cross-section is
// well defined (the region between the arc and its chord; see `full_torus_equal_heights`,
// which OCCT does build), so kentron expects valid solids.

