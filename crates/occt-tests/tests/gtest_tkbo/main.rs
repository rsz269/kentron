//! Hand translations of OCCT's TKBO GTests (src/ModelingAlgorithms/TKBO/GTests).
//!
//! The helpers below mirror `BOPTest_Utilities.pxx` and the `BRepAlgoAPI_TestBase`
//! fixture. Shapes are moved as OCCT's `BRepBuilderAPI_Transform` moves them for rigid
//! transformations: by placement, without copying.

mod brep_algo_api_cut;
mod brep_algo_api_fuse;

use std::f64::consts::PI;

pub use occt_tests::assert_near;
pub use kentron::geometry::surface::Plane;
pub use kentron::math::{Axis, Direction, Frame, Point, Transform, Vector};
pub use kentron::sketch::{Outcome, Step};
pub use kentron::topology::{Model, NodeKind, Shape};

use kentron::primitives::{Cone, Cuboid, Cylinder, Sphere};
use kentron::{blend, boolean, convert, make, properties, sketch, sweep};

/// `gp_Pnt(x, y, z)`.
pub fn pnt(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

/// `gp_Vec(x, y, z)`.
pub fn vec3(x: f64, y: f64, z: f64) -> Vector {
    Vector::new(x, y, z)
}

/// `gp_Dir(x, y, z)`.
pub fn dir(x: f64, y: f64, z: f64) -> Direction {
    Direction::new(x, y, z)
}

/// `BOPTest_Utilities::DefaultTolerance()`.
pub const DEFAULT_TOLERANCE: f64 = 1.0e-6;

pub fn get_surface_area(model: &Model, shape: Shape) -> f64 {
    properties::area(model, shape).mass
}

pub fn get_volume(model: &Model, shape: Shape) -> f64 {
    properties::volume(model, shape).mass
}

/// True if the shape has (next to) no surface area.
pub fn is_empty(model: &Model, shape: Shape, tolerance: f64) -> bool {
    get_surface_area(model, shape) <= tolerance
}

pub fn create_unit_box(model: &mut Model) -> Shape {
    create_box(model, Point::ORIGIN, 1.0, 1.0, 1.0)
}

/// `BRepPrimAPI_MakeBox(corner, dx, dy, dz)`: negative sizes extend backwards from the corner.
pub fn create_box(model: &mut Model, corner: Point, x: f64, y: f64, z: f64) -> Shape {
    Cuboid::from_corners(corner, corner + vec3(x, y, z)).build(model).expect("box")
}

pub fn create_unit_sphere(model: &mut Model) -> Shape {
    Sphere::new(1.0).build(model).expect("sphere")
}

pub fn create_sphere(model: &mut Model, center: Point, radius: f64) -> Shape {
    Sphere::new(radius).with_frame(Frame::new(center, Direction::Z, Direction::X)).build(model).expect("sphere")
}

pub fn create_cylinder(model: &mut Model, radius: f64, height: f64) -> Shape {
    Cylinder::new(radius, height).build(model).expect("cylinder")
}

pub fn create_cone(model: &mut Model, r1: f64, r2: f64, height: f64) -> Shape {
    Cone::new(r1, r2, height).build(model).expect("cone")
}

/// `BRepBuilderAPI_NurbsConvert`.
pub fn convert_to_nurbs(model: &mut Model, shape: Shape) -> Shape {
    convert::to_nurbs(model, shape).expect("NURBS conversion")
}

/// A face on the closed polygon through four points.
pub fn create_polygon_face(model: &mut Model, p1: Point, p2: Point, p3: Point, p4: Point) -> Shape {
    let wire = create_polygon_wire(model, &[p1, p2, p3, p4], true);
    create_face_from_wire(model, wire)
}

pub fn rotate_shape(shape: Shape, axis: Axis, angle: f64) -> Shape {
    shape.moved(&Transform::rotation(axis, angle))
}

pub fn translate_shape(shape: Shape, v: Vector) -> Shape {
    shape.moved(&Transform::translation(v))
}

pub fn create_vertex(model: &mut Model, p: Point) -> Shape {
    make::vertex(model, p)
}

pub fn create_edge(model: &mut Model, p1: Point, p2: Point) -> Shape {
    make::edge_between(model, p1, p2).expect("edge")
}

/// `BRepBuilderAPI_MakePolygon` through `points`, closed if `close`.
pub fn create_polygon_wire(model: &mut Model, points: &[Point], close: bool) -> Shape {
    let n = points.len();
    let count = if close { n } else { n - 1 };
    let edges: Vec<Shape> = (0..count).map(|i| create_edge(model, points[i], points[(i + 1) % n])).collect();
    make::wire(model, &edges).expect("polygon wire")
}

/// `BOPTest_Utilities::CreateProfile`: the face drawn by `steps` in `plane`.
pub fn create_profile(model: &mut Model, plane: Plane, steps: &[Step]) -> Shape {
    sketch::profile(model, *plane.frame(), steps, Outcome::Face).expect("profile")
}

/// `CreateProfile` in the XY plane.
pub fn create_profile_from_steps(model: &mut Model, steps: &[Step]) -> Shape {
    create_profile(model, Plane::new(Point::ORIGIN, Direction::Z), steps)
}

pub fn create_rectangular_profile(model: &mut Model, corner: Point, x: f64, y: f64) -> Shape {
    let points = [corner, corner + vec3(x, 0.0, 0.0), corner + vec3(x, y, 0.0), corner + vec3(0.0, y, 0.0)];
    create_polygon_wire(model, &points, true)
}

/// `BRepBuilderAPI_MakeFace(wire)`.
pub fn create_face_from_wire(model: &mut Model, wire: Shape) -> Shape {
    make::planar_face(model, wire).expect("face from wire")
}

/// `BRepPrimAPI_MakePrism(profile, v)`.
pub fn create_prism(model: &mut Model, profile: Shape, v: Vector) -> Shape {
    sweep::extrude(model, profile, v).expect("prism").shape
}

pub fn create_rectangular_prism(model: &mut Model, corner: Point, x: f64, y: f64, z: f64) -> Shape {
    let wire = create_rectangular_profile(model, corner, x, y);
    let face = create_face_from_wire(model, wire);
    create_prism(model, face, vec3(0.0, 0.0, z))
}

/// `BRepPrimAPI_MakeRevol(profile, axis, angle)`.
pub fn create_revolution(model: &mut Model, profile: Shape, axis: Axis, angle: f64) -> Shape {
    sweep::revolve(model, profile, axis, angle).expect("revolution").shape
}

/// The `index`-th face (from 1) met exploring `shape`, as DRAW's `nexplode`.
pub fn get_face_by_index(model: &Model, shape: Shape, index: usize) -> Shape {
    model.explore(shape, NodeKind::Face).nth(index - 1).expect("face index")
}

fn rotate(shape: Shape, axis: Axis, angle_deg: f64) -> Shape {
    rotate_shape(shape, axis, angle_deg * PI / 180.0)
}

pub fn rotate_z(shape: Shape, angle_deg: f64) -> Shape {
    rotate(shape, Axis::Z, angle_deg)
}

pub fn rotate_y(shape: Shape, angle_deg: f64) -> Shape {
    rotate(shape, Axis::Y, angle_deg)
}

pub fn rotate_x(shape: Shape, angle_deg: f64) -> Shape {
    rotate(shape, Axis::X, angle_deg)
}

/// Z by -90°, then Y by -45°: the rotation many DRAW boolean scripts apply.
pub fn rotate_standard(shape: Shape) -> Shape {
    rotate_y(rotate_z(shape, -90.0), -45.0)
}

pub fn translate(shape: Shape, dx: f64, dy: f64, dz: f64) -> Shape {
    translate_shape(shape, vec3(dx, dy, dz))
}

/// 90° about the Y axis through (0, 0, 1).
pub fn rotate_y90(shape: Shape) -> Shape {
    rotate(shape, Axis::new(pnt(0.0, 0.0, 1.0), Direction::Y), 90.0)
}

pub fn create_identical_boxes(model: &mut Model) -> (Shape, Shape) {
    (create_box(model, Point::ORIGIN, 1.0, 1.0, 1.0), create_box(model, Point::ORIGIN, 1.0, 1.0, 1.0))
}

/// A fillet of radius `radius` on the `edge_index`-th edge (from 1) met exploring `shape`.
pub fn create_blend(model: &mut Model, shape: Shape, edge_index: usize, radius: f64) -> Shape {
    let Some(edge) = model.explore(shape, NodeKind::Edge).nth(edge_index - 1) else { return shape };
    blend::fillet(model, shape, &[edge], radius).expect("fillet")
}

/// `BRepPrimAPI_MakeCylinder` on the plane's frame (DRAW's `pcylinder` on a plane).
pub fn create_cylinder_on_plane(model: &mut Model, plane: Plane, radius: f64, height: f64) -> Shape {
    Cylinder::new(radius, height).with_frame(*plane.frame()).build(model).expect("cylinder")
}

pub fn perform_union(model: &mut Model, a: Shape, b: Shape) -> Shape {
    boolean::union(model, a, b).expect("fuse failed")
}

pub fn perform_difference(model: &mut Model, a: Shape, b: Shape) -> Shape {
    boolean::difference(model, a, b).expect("cut failed")
}

pub fn perform_intersection(model: &mut Model, a: Shape, b: Shape) -> Shape {
    boolean::intersection(model, a, b).expect("common failed")
}

/// `BRepAlgoAPI_TestBase::ValidateResult`: the result's area and volume, where expected
/// (non-negative), or that it is empty. OCCT accepts areas within ±5000 of the expected
/// value, which accepts nearly anything; here they must agree to 1e-4 relative (OCCT's
/// expectations are given to about six digits).
pub fn validate_result(model: &Model, result: Shape, area: f64, volume: f64, empty: bool) {
    if empty {
        assert!(is_empty(model, result, DEFAULT_TOLERANCE), "result should be empty");
        return;
    }
    if area >= 0.0 {
        let actual = get_surface_area(model, result);
        assert!((actual - area).abs() <= 1e-4 * area.max(1.0), "surface area {actual} != {area}");
    }
    if volume >= 0.0 {
        let actual = get_volume(model, result);
        assert!((actual - volume).abs() <= DEFAULT_TOLERANCE, "volume {actual} != {volume}");
    }
}
