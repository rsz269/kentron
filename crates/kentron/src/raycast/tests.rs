use super::*;
use crate::convert::to_nurbs;
use crate::math::Axis;
use crate::primitives::{Cone, Cuboid, Cylinder, Sphere, Torus};

/// The ray's meetings, those within `1e-6` of each other read as one.
fn crossings(model: &Model, shape: Shape, origin: [f64; 3], direction: [f64; 3]) -> Vec<f64> {
    let mut hits = raycast(model, shape, Point::new(origin[0], origin[1], origin[2]), Vector::new(direction[0], direction[1], direction[2])).unwrap();
    hits.dedup_by(|later, kept| *later - *kept <= 1e-6);
    hits
}

fn assert_near(got: &[f64], want: &[f64], tol: f64) {
    assert_eq!(got.len(), want.len(), "{got:?} against {want:?}");
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() <= tol, "{got:?} against {want:?}");
    }
}

#[test]
fn box_faces() {
    let mut model = Model::new();
    let block = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
    assert_near(&crossings(&model, block, [-5.0, 10.0, 15.0], [1.0, 0.0, 0.0]), &[5.0, 15.0], 1e-12);
    // The distances are in millimetres whatever the direction's length.
    assert_near(&crossings(&model, block, [-5.0, 10.0, 15.0], [2.0, 0.0, 0.0]), &[5.0, 15.0], 1e-12);
    // From inside, the one face ahead; a ray past the box, none.
    assert_near(&crossings(&model, block, [5.0, 10.0, 15.0], [0.0, 0.0, -1.0]), &[15.0], 1e-12);
    assert!(crossings(&model, block, [-5.0, 25.0, 15.0], [1.0, 0.0, 0.0]).is_empty());
    assert!(crossings(&model, block, [15.0, 10.0, 15.0], [1.0, 0.0, 0.0]).is_empty());
    // Along a diagonal through two opposite corners.
    let d = (100.0f64 + 400.0 + 900.0).sqrt();
    assert_near(&crossings(&model, block, [-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]), &[0.1 * d, 1.1 * d], 1e-9);
}

#[test]
fn ray_through_an_edge_lists_each_face() {
    let mut model = Model::new();
    let block = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
    // A ray at 45° across the edge x = 10, z = 30 meets the faces x = 10 and z = 30 there.
    let hits = raycast(&model, block, Point::new(1.0, 10.0, 21.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let at = 9.0 * 2f64.sqrt();
    assert!(hits.len() >= 2, "{hits:?}");
    assert!(hits.iter().all(|h| (h - at).abs() <= 1e-9), "{hits:?}");
}

#[test]
fn round_solids() {
    let mut model = Model::new();
    let sphere = Sphere::new(1.5).build(&mut model).unwrap();
    let r = (2.25f64 - 0.25).sqrt();
    assert_near(&crossings(&model, sphere, [-5.0, 0.3, 0.4], [1.0, 0.0, 0.0]), &[5.0 - r, 5.0 + r], 1e-12);
    // Through the poles.
    assert_near(&crossings(&model, sphere, [0.0, 0.0, -4.0], [0.0, 0.0, 1.0]), &[2.5, 5.5], 1e-12);

    let cylinder = Cylinder::new(2.0, 3.0).build(&mut model).unwrap();
    assert_near(&crossings(&model, cylinder, [-5.0, 0.0, 1.5], [1.0, 0.0, 0.0]), &[3.0, 7.0], 1e-12);
    assert_near(&crossings(&model, cylinder, [1.0, 0.0, -1.0], [0.0, 0.0, 1.0]), &[1.0, 4.0], 1e-12);
    // Above the cylinder, nothing.
    assert!(crossings(&model, cylinder, [-5.0, 0.0, 3.5], [1.0, 0.0, 0.0]).is_empty());

    // A cone from radius 2 at z = 0 to 1 at z = 3: radius 1.5 at z = 1.5.
    let cone = Cone::new(2.0, 1.0, 3.0).build(&mut model).unwrap();
    assert_near(&crossings(&model, cone, [-5.0, 0.0, 1.5], [1.0, 0.0, 0.0]), &[3.5, 6.5], 1e-12);

    // A torus of radii 3 and 1 about z: across it in its plane, and at z = 0.5.
    let torus = Torus::new(3.0, 1.0).build(&mut model).unwrap();
    assert_near(&crossings(&model, torus, [-10.0, 0.0, 0.0], [1.0, 0.0, 0.0]), &[6.0, 8.0, 12.0, 14.0], 1e-9);
    let h = 0.75f64.sqrt();
    assert_near(&crossings(&model, torus, [-10.0, 0.0, 0.5], [1.0, 0.0, 0.0]), &[10.0 - 3.0 - h, 10.0 - 3.0 + h, 10.0 + 3.0 - h, 10.0 + 3.0 + h], 1e-9);
    // Down through the hole, nothing; down through the tube, its top and bottom.
    assert!(crossings(&model, torus, [0.0, 0.0, 5.0], [0.0, 0.0, -1.0]).is_empty());
    assert_near(&crossings(&model, torus, [3.0, 0.0, 5.0], [0.0, 0.0, -1.0]), &[4.0, 6.0], 1e-9);
}

#[test]
fn sampled_surfaces() {
    // The same solids as B-spline surfaces: found through the polyhedron and refined.
    let mut model = Model::new();
    let sphere = Sphere::new(1.5).build(&mut model).unwrap();
    let sphere = to_nurbs(&mut model, sphere).unwrap();
    let r = (2.25f64 - 0.25).sqrt();
    assert_near(&crossings(&model, sphere, [-5.0, 0.3, 0.4], [1.0, 0.0, 0.0]), &[5.0 - r, 5.0 + r], 1e-9);
    let torus = Torus::new(3.0, 1.0).build(&mut model).unwrap();
    let torus = to_nurbs(&mut model, torus).unwrap();
    assert_near(&crossings(&model, torus, [-10.0, 0.2, 0.5], [1.0, 0.0, 0.0]), &torus_hits(0.2, 0.5), 1e-9);
    let block = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
    let block = to_nurbs(&mut model, block).unwrap();
    assert_near(&crossings(&model, block, [-5.0, 10.0, 15.0], [1.0, 0.0, 0.0]), &[5.0, 15.0], 1e-9);
    assert_near(&crossings(&model, block, [3.0, 7.0, 40.0], [0.0, 0.0, -1.0]), &[10.0, 40.0], 1e-9);
}

/// Where the line `(-10 + t, y, z)` meets the torus of radii 3 and 1 about z: the roots of
/// `(√(x² + y²) − 3)² + z² = 1`.
fn torus_hits(y: f64, z: f64) -> Vec<f64> {
    let mut out = Vec::new();
    for rho in [3.0 - (1.0 - z * z).sqrt(), 3.0 + (1.0 - z * z).sqrt()] {
        let x = (rho * rho - y * y).sqrt();
        out.extend([10.0 - x, 10.0 + x]);
    }
    out.sort_by(f64::total_cmp);
    out
}

#[test]
fn placed_shapes() {
    // A box turned a quarter about z and moved: the ray meets it where the turned box is.
    let mut model = Model::new();
    let block = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
    let turn = Transform::rotation(Axis::Z, std::f64::consts::FRAC_PI_2);
    let placed = block.moved(&turn).moved(&Transform::translation(Vector::new(100.0, 0.0, 0.0)));
    // Turned, it spans x in [80, 100], y in [0, 10].
    assert_near(&crossings(&model, placed, [0.0, 5.0, 15.0], [1.0, 0.0, 0.0]), &[80.0, 100.0], 1e-9);
    let raycaster = Raycaster::new(&model, placed).unwrap();
    let all = raycaster.cast_all(&[(Point::new(0.0, 5.0, 15.0), Vector::new(1.0, 0.0, 0.0)), (Point::new(90.0, 5.0, -10.0), Vector::new(0.0, 0.0, 1.0))]).unwrap();
    assert_near(&all[0], &[80.0, 100.0], 1e-9);
    assert_near(&all[1], &[10.0, 40.0], 1e-9);
}

#[test]
fn zero_direction_is_refused() {
    let mut model = Model::new();
    let block = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
    assert!(raycast(&model, block, Point::ORIGIN, Vector::new(0.0, 0.0, 0.0)).is_err());
}

#[test]
fn far_from_the_origin() {
    // A cylinder and a sphere some 700 km out: the meetings as exact as near the origin.
    let mut model = Model::new();
    let at = Vector::new(-335418.7, 104584.2, 625130.0);
    let cylinder = Cylinder::new(2.0, 4.0).build(&mut model).unwrap().moved(&Transform::translation(at));
    let o = Point::ORIGIN + at + Vector::new(-10.0, 0.3, 2.0);
    let h = (4.0f64 - 0.09).sqrt();
    assert_near(&crossings(&model, cylinder, [o.x, o.y, o.z], [1.0, 0.0, 0.0]), &[10.0 - h, 10.0 + h], 1e-9);
    let sphere = Sphere::new(2.0).build(&mut model).unwrap().moved(&Transform::translation(at));
    let o = Point::ORIGIN + at + Vector::new(0.3, -10.0, 0.4);
    let h = (4.0f64 - 0.25).sqrt();
    assert_near(&crossings(&model, sphere, [o.x, o.y, o.z], [0.0, 1.0, 0.0]), &[10.0 - h, 10.0 + h], 1e-9);
}
