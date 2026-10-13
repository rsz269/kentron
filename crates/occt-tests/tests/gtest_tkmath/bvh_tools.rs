//! BVH_Tools_Test.cxx (BVH_Tools: distances and projections of boxes, rays and triangles)
//!
//! `PointBoxSquareDistance` and `BoxBoxSquareDistance` are `BoundingBox::distance_to_point`
//! and `distance_to_box` (squared here, as OCCT's), `PointBoxProjection` is `closest_point`,
//! `RayBoxIntersection` is `ray_interval` (`None` for a miss), `PointTriangleProjection` is
//! `math::closest_point_on_triangle`.
//!
//! Not ported: the `*Float` tests (no single-precision boxes).

use occt_tests::assert_near;
use kentron::math::{closest_point_on_triangle, BoundingBox, BoundingBox2, Point, Point2, Vector, Vector2};
use kentron::tolerance;

const TOL: f64 = tolerance::LINEAR;

fn bx(min: (f64, f64, f64), max: (f64, f64, f64)) -> BoundingBox {
    BoundingBox::new(Point::new(min.0, min.1, min.2), Point::new(max.0, max.1, max.2))
}

fn bx2(min: (f64, f64), max: (f64, f64)) -> BoundingBox2 {
    BoundingBox2::new(Point2::new(min.0, min.1), Point2::new(max.0, max.1))
}

fn unit() -> BoundingBox {
    bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
}

fn pt(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

fn sq_point(b: &BoundingBox, p: Point) -> f64 {
    b.distance_to_point(p).powi(2)
}

fn sq_box(a: &BoundingBox, b: &BoundingBox) -> f64 {
    a.distance_to_box(b).powi(2)
}

fn ray(b: &BoundingBox, o: Point, d: Vector) -> Option<(f64, f64)> {
    b.ray_interval(o, d)
}

fn assert_point(p: Point, x: f64, y: f64, z: f64) {
    assert_near!(p.x, x, TOL);
    assert_near!(p.y, y, TOL);
    assert_near!(p.z, z, TOL);
}

fn diagonal() -> Vector {
    Vector::new(1.0, 1.0, 1.0) * (1.0 / 3f64.sqrt())
}

#[test]
fn point_box_square_distance() {
    let b = unit();
    assert_near!(sq_point(&b, pt(0.5, 0.5, 0.5)), 0.0, TOL);
    assert_near!(sq_point(&b, pt(0.5, 0.5, 1.0)), 0.0, TOL);
    assert_near!(sq_point(&b, pt(2.0, 0.5, 0.5)), 1.0, TOL);
    assert_near!(sq_point(&b, pt(2.0, 2.0, 2.0)), 3.0, TOL);
}

#[test]
fn box_box_square_distance() {
    let b1 = unit();
    assert_near!(sq_box(&b1, &bx((2.0, 0.0, 0.0), (3.0, 1.0, 1.0))), 1.0, TOL);
    assert_near!(sq_box(&b1, &bx((0.5, 0.5, 0.5), (1.5, 1.5, 1.5))), 0.0, TOL);
}

#[test]
fn point_box_projection() {
    let b = unit();
    assert_point(b.closest_point(pt(0.5, 0.5, 0.5)), 0.5, 0.5, 0.5);
    assert_point(b.closest_point(pt(2.0, 0.5, 0.5)), 1.0, 0.5, 0.5);
}

#[test]
fn ray_box_intersection() {
    let b = unit();
    let (enter, leave) = ray(&b, pt(-1.0, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
    assert!(ray(&b, pt(-1.0, 5.0, 0.5), Vector::new(1.0, 0.0, 0.0)).is_none());
    // From inside.
    let (enter, leave) = ray(&b, pt(0.5, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).expect("a hit");
    assert!(enter <= 0.0);
    assert_near!(leave, 0.5, TOL);
}

#[test]
fn ray_box_intersection_parallel_ray() {
    let b = unit();
    assert!(ray(&b, pt(-1.0, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).is_some());
    assert!(ray(&b, pt(-1.0, 2.0, 0.5), Vector::new(1.0, 0.0, 0.0)).is_none());
}

fn triangle() -> (Point, Point, Point) {
    (pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0))
}

fn triangle2() -> (Point, Point, Point) {
    (pt(0.0, 0.0, 0.0), pt(2.0, 0.0, 0.0), pt(0.0, 2.0, 0.0))
}

fn project(p: Point, (a, b, c): (Point, Point, Point)) -> Point {
    closest_point_on_triangle(p, a, b, c)
}

#[test]
fn point_triangle_projection() {
    assert_point(project(pt(-1.0, -1.0, 0.0), triangle()), 0.0, 0.0, 0.0);
    assert_point(project(pt(0.5, -1.0, 0.0), triangle()), 0.5, 0.0, 0.0);
    assert_point(project(pt(0.25, 0.25, 1.0), triangle()), 0.25, 0.25, 0.0);
}

#[test]
fn point_box_square_distance2_d() {
    let b = bx2((0.0, 0.0), (1.0, 1.0));
    assert_near!(b.distance_to_point(Point2::new(0.5, 0.5)).powi(2), 0.0, TOL);
    assert_near!(b.distance_to_point(Point2::new(2.0, 0.5)).powi(2), 1.0, TOL);
}

#[test]
fn ray_box_intersection_diagonal_ray() {
    let (enter, leave) = ray(&unit(), pt(-1.0, -1.0, -1.0), diagonal()).expect("a hit");
    assert!(enter > 0.0);
    assert!(leave > enter);
}

#[test]
fn ray_box_intersection_negative_direction() {
    let (enter, leave) = ray(&unit(), pt(2.0, 0.5, 0.5), Vector::new(-1.0, 0.0, 0.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
}

#[test]
fn ray_box_intersection_touching_edge() {
    assert!(ray(&unit(), pt(-1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).is_some());
}

#[test]
fn ray_box_intersection_touching_corner() {
    assert!(ray(&unit(), pt(-1.0, -1.0, -1.0), diagonal()).is_some());
}

#[test]
fn box_box_square_distance_touching() {
    assert_near!(sq_box(&unit(), &bx((1.0, 0.0, 0.0), (2.0, 1.0, 1.0))), 0.0, TOL);
}

#[test]
fn box_box_square_distance_one_inside_other() {
    assert_near!(sq_box(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((2.0, 2.0, 2.0), (3.0, 3.0, 3.0))), 0.0, TOL);
}

#[test]
fn box_box_square_distance_corner_to_corner() {
    assert_near!(sq_box(&unit(), &bx((2.0, 2.0, 2.0), (3.0, 3.0, 3.0))), 3.0, TOL);
}

#[test]
fn point_box_projection_negative_coords() {
    assert_point(bx((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)).closest_point(pt(-5.0, -5.0, -5.0)), -1.0, -1.0, -1.0);
}

#[test]
fn point_triangle_projection_on_edge01() {
    assert_point(project(pt(1.0, -1.0, 0.0), triangle2()), 1.0, 0.0, 0.0);
}

#[test]
fn point_triangle_projection_on_edge12() {
    assert_point(project(pt(2.0, 2.0, 0.0), triangle2()), 1.0, 1.0, 0.0);
}

#[test]
fn point_triangle_projection_on_edge20() {
    assert_point(project(pt(-1.0, 1.0, 0.0), triangle2()), 0.0, 1.0, 0.0);
}

#[test]
fn point_triangle_projection_vertex1() {
    assert_point(project(pt(3.0, -1.0, 0.0), triangle2()), 2.0, 0.0, 0.0);
}

#[test]
fn point_triangle_projection_vertex2() {
    assert_point(project(pt(-1.0, 3.0, 0.0), triangle2()), 0.0, 2.0, 0.0);
}

#[test]
fn point_triangle_projection3_d() {
    let p = closest_point_on_triangle(pt(0.3, 0.3, 0.8), pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 1.0));
    assert!(p.x >= 0.0);
    assert!(p.y >= 0.0);
    assert!(p.x + p.y <= 1.0 + TOL);
}

#[test]
fn box_box_square_distance2_d() {
    assert_near!(bx2((0.0, 0.0), (1.0, 1.0)).distance_to_box(&bx2((3.0, 0.0), (4.0, 1.0))).powi(2), 4.0, TOL);
}

#[test]
fn point_box_projection2_d() {
    let p = bx2((0.0, 0.0), (1.0, 1.0)).closest_point(Point2::new(2.0, 2.0));
    assert_near!(p.x, 1.0, TOL);
    assert_near!(p.y, 1.0, TOL);
}

#[test]
fn point_box_square_distance_at_vertex() {
    assert_near!(sq_point(&unit(), pt(0.0, 0.0, 0.0)), 0.0, TOL);
}

#[test]
fn point_box_square_distance_at_edge() {
    assert_near!(sq_point(&unit(), pt(0.5, 0.0, 0.0)), 0.0, TOL);
}

#[test]
fn ray_box_intersection_behind_ray() {
    // The box behind the origin: no hit.
    assert!(ray(&unit(), pt(2.0, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).is_none_or(|(_, leave)| leave < 0.0));
}

#[test]
fn ray_box_intersection_y_axis() {
    let (enter, leave) = ray(&unit(), pt(0.5, -1.0, 0.5), Vector::new(0.0, 1.0, 0.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
}

#[test]
fn ray_box_intersection_z_axis() {
    let (enter, leave) = ray(&unit(), pt(0.5, 0.5, -1.0), Vector::new(0.0, 0.0, 1.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
}

#[test]
fn point_box_square_distance_negative_box() {
    assert_near!(sq_point(&bx((-2.0, -2.0, -2.0), (-1.0, -1.0, -1.0)), pt(0.0, 0.0, 0.0)), 3.0, TOL);
}

#[test]
fn box_box_square_distance_edge_to_edge() {
    assert_near!(sq_box(&unit(), &bx((2.0, 2.0, 0.0), (3.0, 3.0, 1.0))), 2.0, TOL);
}

#[test]
fn point_box_square_distance_large_box() {
    let b = bx((0.0, 0.0, 0.0), (1000.0, 1000.0, 1000.0));
    assert_near!(sq_point(&b, pt(500.0, 500.0, 500.0)), 0.0, TOL);
    assert_near!(sq_point(&b, pt(1001.0, 500.0, 500.0)), 1.0, TOL);
}

#[test]
fn ray_box_intersection_large_box() {
    let b = bx((0.0, 0.0, 0.0), (1000.0, 1000.0, 1000.0));
    let (enter, leave) = ray(&b, pt(-100.0, 500.0, 500.0), Vector::new(1.0, 0.0, 0.0)).expect("a hit");
    assert_near!(enter, 100.0, TOL);
    assert_near!(leave, 1100.0, TOL);
}

#[test]
fn point_triangle_projection_centroid() {
    let p = closest_point_on_triangle(pt(1.0, 1.0, 5.0), pt(0.0, 0.0, 0.0), pt(3.0, 0.0, 0.0), pt(0.0, 3.0, 0.0));
    assert_point(p, 1.0, 1.0, 0.0);
}

#[test]
fn box_box_square_distance_same_box() {
    assert_near!(sq_box(&unit(), &unit()), 0.0, TOL);
}

#[test]
fn point_box_projection_all_corners() {
    assert_point(unit().closest_point(pt(-1.0, -1.0, -1.0)), 0.0, 0.0, 0.0);
    assert_point(unit().closest_point(pt(2.0, 2.0, 2.0)), 1.0, 1.0, 1.0);
}

#[test]
fn ray_box_intersection2_d() {
    let b = bx2((0.0, 0.0), (1.0, 1.0));
    let (enter, leave) = b.ray_interval(Point2::new(-1.0, 0.5), Vector2::new(1.0, 0.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
    assert!(b.ray_interval(Point2::new(-1.0, 2.0), Vector2::new(1.0, 0.0)).is_none());
    assert!(b.ray_interval(Point2::new(-1.0, 0.5), Vector2::new(1.0, 0.0)).is_some());
    let (enter, leave) = b.ray_interval(Point2::new(0.5, -1.0), Vector2::new(0.0, 1.0)).expect("a hit");
    assert_near!(enter, 1.0, TOL);
    assert_near!(leave, 2.0, TOL);
}

#[test]
fn ray_box_intersection_early_exit() {
    assert!(ray(&unit(), pt(-2.0, 0.5, 0.5), Vector::new(-1.0, 0.0, 0.0)).is_none());
    assert!(ray(&unit(), pt(0.5, -2.0, 0.5), Vector::new(0.0, -1.0, 0.0)).is_none());
}

#[test]
fn ray_box_intersection_parallel_ray_early_exit() {
    assert!(ray(&unit(), pt(-1.0, 2.0, 0.5), Vector::new(1.0, 0.0, 0.0)).is_none());
    assert!(ray(&unit(), pt(2.0, -1.0, 0.5), Vector::new(0.0, 1.0, 0.0)).is_none());
    assert!(ray(&unit(), pt(-1.0, -1.0, 0.5), Vector::new(0.0, 0.0, 1.0)).is_none());
}

#[test]
fn ray_box_intersection2_d_parallel_both_axes() {
    let b = bx2((0.0, 0.0), (1.0, 1.0));
    // A null direction: a hit if and only if the origin is inside.
    assert!(b.ray_interval(Point2::new(0.5, 0.5), Vector2::new(0.0, 0.0)).is_some());
    assert!(b.ray_interval(Point2::new(2.0, 2.0), Vector2::new(0.0, 0.0)).is_none());
}

#[test]
fn ray_box_intersection_negative_time() {
    assert!(ray(&unit(), pt(2.0, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).is_none());
    let (_, leave) = ray(&unit(), pt(0.5, 0.5, 0.5), Vector::new(1.0, 0.0, 0.0)).expect("a hit");
    assert!(leave >= 0.0);
}
