//! BVH_Box_Test.cxx (BVH_Box: `math::BoundingBox` and `BoundingBox2`)
//!
//! `IsValid` is `!is_empty()`, `Combine` is `add_box`, `IsOut(box)` is `!overlaps`,
//! `Contains(box, hasOverlap)` is `(contains_box, overlaps)`, `Clear` is `EMPTY`, `Center(axis)`
//! is a coordinate of `center()`. A 4×4 matrix transformation is an `AffineTransform`.
//!
//! Not ported: `Box4D`, `CenterByAxis4D` (no 4-dimensional boxes); the `Float*` and
//! `*_Float_*` tests (no single-precision boxes); `Constexpr_*` (compile-time construction is
//! `BoundingBox::new` and `EMPTY`, both `const`, exercised by the other tests); `ModifyCorner*`
//! (the corners are public fields: setting one is plain assignment).
//!
//! Deviation: `Area_DegenerateBox2D` expects OCCT's surface-area-heuristic value for a flat
//! rectangle (width + height, 5); kentron's `area()` is the rectangle's area, 0.

use occt_tests::assert_near;
use kentron::math::{
    AffineTransform, AffineTransformable, BoundingBox, BoundingBox2, Matrix3, Point, Point2, Transform, Transformable,
    Vector,
};
use kentron::tolerance;

const TOL: f64 = tolerance::LINEAR;

fn bx(min: (f64, f64, f64), max: (f64, f64, f64)) -> BoundingBox {
    BoundingBox::new(Point::new(min.0, min.1, min.2), Point::new(max.0, max.1, max.2))
}

fn bx2(min: (f64, f64), max: (f64, f64)) -> BoundingBox2 {
    BoundingBox2::new(Point2::new(min.0, min.1), Point2::new(max.0, max.1))
}

fn assert_point(p: Point, x: f64, y: f64, z: f64) {
    assert_near!(p.x, x, TOL);
    assert_near!(p.y, y, TOL);
    assert_near!(p.z, z, TOL);
}

fn assert_size(v: Vector, x: f64, y: f64, z: f64) {
    assert_near!(v.x, x, TOL);
    assert_near!(v.y, y, TOL);
    assert_near!(v.z, z, TOL);
}

#[test]
fn default_constructor() {
    assert!(BoundingBox::default().is_empty());
}

#[test]
fn constructor_with_corners() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 2.0, 3.0));
    assert!(!b.is_empty());
    assert_point(b.min, 0.0, 0.0, 0.0);
    assert_point(b.max, 1.0, 2.0, 3.0);
}

#[test]
fn add() {
    let mut b = BoundingBox::EMPTY;
    b.add_point(Point::new(1.0, 2.0, 3.0));
    assert!(!b.is_empty());
    b.add_point(Point::new(-1.0, -2.0, -3.0));
    assert_point(b.min, -1.0, -2.0, -3.0);
    assert_point(b.max, 1.0, 2.0, 3.0);
}

#[test]
fn combine() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    b.add_box(&bx((2.0, 2.0, 2.0), (3.0, 3.0, 3.0)));
    assert_point(b.min, 0.0, 0.0, 0.0);
    assert_point(b.max, 3.0, 3.0, 3.0);
}

#[test]
fn size() {
    assert_size(bx((0.0, 0.0, 0.0), (2.0, 3.0, 4.0)).size(), 2.0, 3.0, 4.0);
}

#[test]
fn center() {
    assert_point(bx((0.0, 0.0, 0.0), (2.0, 4.0, 6.0)).center(), 1.0, 2.0, 3.0);
}

#[test]
fn area() {
    // 2 (1·2 + 2·3 + 1·3) = 22
    assert_near!(bx((0.0, 0.0, 0.0), (1.0, 2.0, 3.0)).area(), 22.0, TOL);
}

#[test]
fn is_out() {
    let b1 = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    assert!(!b1.overlaps(&bx((2.0, 2.0, 2.0), (3.0, 3.0, 3.0))));
    assert!(b1.overlaps(&bx((0.5, 0.5, 0.5), (1.5, 1.5, 1.5))));
}

#[test]
fn clear() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    assert!(!b.is_empty());
    b = BoundingBox::EMPTY;
    assert!(b.is_empty());
}

#[test]
fn box2_d() {
    let mut b = BoundingBox2::EMPTY;
    b.add_point(Point2::new(0.0, 0.0));
    b.add_point(Point2::new(1.0, 1.0));
    assert!(!b.is_empty());
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.min.y, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
    assert_near!(b.max.y, 1.0, TOL);
    assert_near!(b.area(), 1.0, TOL);
}

#[test]
fn single_point_box() {
    let mut b = BoundingBox::EMPTY;
    b.add_point(Point::new(5.0, 5.0, 5.0));
    assert!(!b.is_empty());
    assert_near!(b.min.x, 5.0, TOL);
    assert_near!(b.max.x, 5.0, TOL);
    assert_size(b.size(), 0.0, 0.0, 0.0);
    assert_near!(b.area(), 0.0, TOL);
}

#[test]
fn negative_coordinates() {
    let b = bx((-10.0, -20.0, -30.0), (-5.0, -10.0, -15.0));
    assert!(!b.is_empty());
    assert_point(b.center(), -7.5, -15.0, -22.5);
    assert_size(b.size(), 5.0, 10.0, 15.0);
}

#[test]
fn large_values() {
    let b = bx((1e10, 1e10, 1e10), (1e10 + 1.0, 1e10 + 2.0, 1e10 + 3.0));
    assert!(!b.is_empty());
    let s = b.size();
    assert_near!(s.x, 1.0, 1e-5);
    assert_near!(s.y, 2.0, 1e-5);
    assert_near!(s.z, 3.0, 1e-5);
}

#[test]
fn combine_with_invalid() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    b.add_box(&BoundingBox::EMPTY);
    assert!(!b.is_empty());
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
}

#[test]
fn add_to_invalid() {
    let mut b = BoundingBox::EMPTY;
    b.add_box(&bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)));
    assert!(!b.is_empty());
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
}

#[test]
fn is_out_touching_boxes() {
    assert!(bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).overlaps(&bx((1.0, 0.0, 0.0), (2.0, 1.0, 1.0))));
}

#[test]
fn is_out_touching_at_edge() {
    assert!(bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).overlaps(&bx((1.0, 1.0, 0.0), (2.0, 2.0, 1.0))));
}

#[test]
fn is_out_touching_at_corner() {
    assert!(bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).overlaps(&bx((1.0, 1.0, 1.0), (2.0, 2.0, 2.0))));
}

#[test]
fn area_unit_cube() {
    assert_near!(bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).area(), 6.0, TOL);
}

#[test]
fn center_at_origin() {
    assert_point(bx((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)).center(), 0.0, 0.0, 0.0);
}

#[test]
fn multiple_adds() {
    let b = BoundingBox::from_points([
        Point::new(5.0, 3.0, 1.0),
        Point::new(-2.0, 7.0, 4.0),
        Point::new(1.0, -1.0, 8.0),
        Point::new(0.0, 2.0, -3.0),
    ]);
    assert_point(b.min, -2.0, -1.0, -3.0);
    assert_point(b.max, 5.0, 7.0, 8.0);
}

#[test]
fn flat_box2_d() {
    let b = bx((0.0, 0.0, 5.0), (3.0, 4.0, 5.0));
    assert!(!b.is_empty());
    assert_size(b.size(), 3.0, 4.0, 0.0);
    assert_near!(b.area(), 24.0, TOL);
}

#[test]
fn flat_box1_d() {
    let b = bx((0.0, 5.0, 5.0), (10.0, 5.0, 5.0));
    assert!(!b.is_empty());
    assert_size(b.size(), 10.0, 0.0, 0.0);
    assert!(b.area() >= 0.0);
}

#[test]
fn combine_multiple() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    b.add_box(&bx((2.0, 0.0, 0.0), (3.0, 1.0, 1.0)));
    b.add_box(&bx((0.0, 2.0, 0.0), (1.0, 3.0, 1.0)));
    assert_point(b.min, 0.0, 0.0, 0.0);
    assert_point(b.max, 3.0, 3.0, 1.0);
}

#[test]
fn is_out_partial_overlap() {
    let (b1, b2) = (bx((0.0, 0.0, 0.0), (2.0, 2.0, 2.0)), bx((1.0, 1.0, 1.0), (3.0, 3.0, 3.0)));
    assert!(b1.overlaps(&b2));
    assert!(b2.overlaps(&b1));
}

#[test]
fn is_out_contained() {
    let (b1, b2) = (bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), bx((2.0, 2.0, 2.0), (3.0, 3.0, 3.0)));
    assert!(b1.overlaps(&b2));
    assert!(b2.overlaps(&b1));
}

#[test]
fn box2_d_is_out() {
    let b1 = bx2((0.0, 0.0), (1.0, 1.0));
    assert!(!b1.overlaps(&bx2((2.0, 2.0), (3.0, 3.0))));
    assert!(b1.overlaps(&bx2((0.5, 0.5), (1.5, 1.5))));
}

#[test]
fn box2_d_center() {
    let c = bx2((0.0, 0.0), (4.0, 6.0)).center();
    assert_near!(c.x, 2.0, TOL);
    assert_near!(c.y, 3.0, TOL);
}

#[test]
fn box2_d_size() {
    let s = bx2((1.0, 2.0), (4.0, 7.0)).size();
    assert_near!(s.x, 3.0, TOL);
    assert_near!(s.y, 5.0, TOL);
}

#[test]
fn very_small_box() {
    let small = 1e-10;
    let b = bx((0.0, 0.0, 0.0), (small, small, small));
    assert!(!b.is_empty());
    let s = b.size();
    assert_near!(s.x, small, 1e-15);
    assert_near!(s.y, small, 1e-15);
    assert_near!(s.z, small, 1e-15);
}

#[test]
fn symmetric_box() {
    let b = bx((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0));
    assert_point(b.center(), 0.0, 0.0, 0.0);
    assert_size(b.size(), 10.0, 10.0, 10.0);
}

#[test]
fn non_cubic_box() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 10.0, 100.0));
    assert_size(b.size(), 1.0, 10.0, 100.0);
    assert_near!(b.area(), 2220.0, TOL);
}

#[test]
fn clear_and_reuse() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    assert!(!b.is_empty());
    b = BoundingBox::EMPTY;
    assert!(b.is_empty());
    b.add_point(Point::new(5.0, 5.0, 5.0));
    b.add_point(Point::new(10.0, 10.0, 10.0));
    assert!(!b.is_empty());
    assert_near!(b.min.x, 5.0, TOL);
    assert_near!(b.max.x, 10.0, TOL);
}

#[test]
fn is_out_same_box() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    assert!(b.overlaps(&b));
}

#[test]
fn combine_same_box() {
    let mut b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    b.add_box(&bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)));
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
}

#[test]
fn is_out_near_miss() {
    assert!(!bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).overlaps(&bx((1.0001, 0.0, 0.0), (2.0, 1.0, 1.0))));
}

#[test]
fn add_duplicate_point() {
    let mut b = BoundingBox::EMPTY;
    for _ in 0..3 {
        b.add_point(Point::new(1.0, 2.0, 3.0));
    }
    assert_near!(b.min.x, 1.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
    assert_near!(b.size().x, 0.0, TOL);
}

#[test]
fn single_point_constructor() {
    let p = Point::new(5.0, 10.0, 15.0);
    let b = BoundingBox::new(p, p);
    assert!(!b.is_empty());
    assert_point(b.min, 5.0, 10.0, 15.0);
    assert_point(b.max, 5.0, 10.0, 15.0);
    assert_size(b.size(), 0.0, 0.0, 0.0);
}

#[test]
fn single_point_constructor2_d() {
    let p = Point2::new(3.0, 7.0);
    let b = BoundingBox2::new(p, p);
    assert!(!b.is_empty());
    assert_near!(b.min.x, 3.0, TOL);
    assert_near!(b.min.y, 7.0, TOL);
    assert_near!(b.max.x, 3.0, TOL);
    assert_near!(b.max.y, 7.0, TOL);
}

#[test]
fn single_point_constructor_origin() {
    let b = BoundingBox::new(Point::ORIGIN, Point::ORIGIN);
    assert!(!b.is_empty());
    assert_point(b.center(), 0.0, 0.0, 0.0);
}

#[test]
fn center_by_axis() {
    let c = bx((0.0, 0.0, 0.0), (10.0, 20.0, 30.0)).center();
    assert_near!(c.x, 5.0, TOL);
    assert_near!(c.y, 10.0, TOL);
    assert_near!(c.z, 15.0, TOL);
}

#[test]
fn center_by_axis2_d() {
    let c = bx2((2.0, 4.0), (8.0, 12.0)).center();
    assert_near!(c.x, 5.0, TOL);
    assert_near!(c.y, 8.0, TOL);
}

#[test]
fn center_by_axis_negative() {
    assert_point(bx((-10.0, -20.0, -30.0), (10.0, 20.0, 30.0)).center(), 0.0, 0.0, 0.0);
}

#[test]
fn is_out_point_inside() {
    let b = bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    for p in [(5.0, 5.0, 5.0), (1.0, 1.0, 1.0), (9.0, 9.0, 9.0)] {
        assert!(b.contains(Point::new(p.0, p.1, p.2)));
    }
}

#[test]
fn is_out_point_outside() {
    let b = bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    for p in [(11.0, 5.0, 5.0), (5.0, 11.0, 5.0), (5.0, 5.0, 11.0), (-1.0, 5.0, 5.0), (15.0, 15.0, 15.0)] {
        assert!(!b.contains(Point::new(p.0, p.1, p.2)));
    }
}

#[test]
fn is_out_point_on_boundary() {
    let b = bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    for p in [(0.0, 5.0, 5.0), (10.0, 5.0, 5.0), (5.0, 0.0, 5.0), (5.0, 10.0, 5.0), (0.0, 0.0, 0.0), (10.0, 10.0, 10.0)] {
        assert!(b.contains(Point::new(p.0, p.1, p.2)));
    }
}

#[test]
fn is_out_point_invalid_box() {
    assert!(!BoundingBox::EMPTY.contains(Point::ORIGIN));
    assert!(!BoundingBox::EMPTY.contains(Point::new(5.0, 5.0, 5.0)));
}

#[test]
fn is_out_point_2_d() {
    let b = bx2((0.0, 0.0), (5.0, 5.0));
    assert!(b.contains(Point2::new(2.5, 2.5)));
    assert!(!b.contains(Point2::new(6.0, 2.5)));
    assert!(!b.contains(Point2::new(2.5, 6.0)));
    assert!(b.contains(Point2::new(0.0, 0.0)));
}

#[test]
fn is_out_point_negative_coords() {
    let b = bx((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0));
    assert!(b.contains(Point::ORIGIN));
    assert!(b.contains(Point::new(-3.0, -3.0, -3.0)));
    assert!(!b.contains(Point::new(-6.0, 0.0, 0.0)));
}

/// `Contains(other, hasOverlap)`: (contained, overlapping).
fn contains(a: &BoundingBox, b: &BoundingBox) -> (bool, bool) {
    (a.contains_box(b), a.overlaps(b))
}

#[test]
fn contains_fully_contained() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((2.0, 2.0, 2.0), (8.0, 8.0, 8.0))), (true, true));
}

#[test]
fn contains_not_contained() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((8.0, 8.0, 8.0), (12.0, 12.0, 12.0))), (false, true));
}

#[test]
fn contains_no_overlap() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (5.0, 5.0, 5.0)), &bx((10.0, 10.0, 10.0), (15.0, 15.0, 15.0))), (false, false));
}

#[test]
fn contains_same_box() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))), (true, true));
}

#[test]
fn contains_touching_boundary() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((0.0, 0.0, 0.0), (5.0, 5.0, 5.0))), (true, true));
}

#[test]
fn contains_invalid_box() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &BoundingBox::EMPTY), (false, false));
}

#[test]
fn contains_invalid_container() {
    assert_eq!(contains(&BoundingBox::EMPTY, &bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))), (false, false));
}

#[test]
fn contains_partial_overlap_one_axis() {
    assert_eq!(contains(&bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)), &bx((5.0, 5.0, 5.0), (15.0, 8.0, 8.0))), (false, true));
}

#[test]
fn contains_2_d() {
    let (a, b) = (bx2((0.0, 0.0), (10.0, 10.0)), bx2((2.0, 2.0), (8.0, 8.0)));
    assert!(a.contains_box(&b));
    assert!(a.overlaps(&b));
}

#[test]
fn contains_by_corners_fully_contained() {
    let b = bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    assert_eq!(contains(&b, &bx((2.0, 2.0, 2.0), (8.0, 8.0, 8.0))), (true, true));
}

#[test]
fn contains_by_corners_not_contained() {
    let b = bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    assert_eq!(contains(&b, &bx((8.0, 8.0, 8.0), (12.0, 12.0, 12.0))), (false, true));
}

/// A 4×4 matrix's upper rows: a linear part and a translation.
fn affine(diagonal: [f64; 3], translation: [f64; 3]) -> AffineTransform {
    let [a, b, c] = diagonal;
    AffineTransform::new(
        Matrix3::from_rows([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]),
        Vector::new(translation[0], translation[1], translation[2]),
    )
}

#[test]
fn transform_identity() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).transformed(&Transform::IDENTITY);
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
}

#[test]
fn transform_translation() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).transformed(&Transform::translation(Vector::new(5.0, 10.0, 15.0)));
    assert_point(b.min, 5.0, 10.0, 15.0);
    assert_point(b.max, 6.0, 11.0, 16.0);
}

#[test]
fn transform_scale() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).affine_transformed(&affine([2.0, 3.0, 4.0], [0.0; 3]));
    assert_near!(b.min.x, 0.0, TOL);
    assert_point(b.max, 2.0, 3.0, 4.0);
}

#[test]
fn transformed_identity() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let t = b.affine_transformed(&affine([1.0; 3], [0.0; 3]));
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
    assert_near!(t.min.x, 0.0, TOL);
    assert_near!(t.max.x, 1.0, TOL);
}

#[test]
fn transformed_translation() {
    let b = bx((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let t = b.affine_transformed(&affine([1.0; 3], [10.0, 20.0, 30.0]));
    assert_near!(b.min.x, 0.0, TOL);
    assert_near!(b.max.x, 1.0, TOL);
    assert_point(t.min, 10.0, 20.0, 30.0);
    assert_point(t.max, 11.0, 21.0, 31.0);
}

#[test]
fn transform_invalid_box() {
    assert!(BoundingBox::EMPTY.transformed(&Transform::translation(Vector::new(10.0, 20.0, 30.0))).is_empty());
}

#[test]
fn transformed_invalid_box() {
    assert!(BoundingBox::EMPTY.affine_transformed(&affine([1.0; 3], [10.0, 20.0, 30.0])).is_empty());
}

#[test]
fn area_degenerate_box2_d() {
    // Deviation (see the module comment): OCCT gives width + height (5).
    let area = bx2((0.0, 0.0), (5.0, 0.0)).area();
    assert!(area >= 0.0);
    assert_near!(area, 0.0, TOL);
}

#[test]
fn is_out_two_corners_overlapping() {
    assert!(bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)).overlaps(&bx((5.0, 5.0, 5.0), (15.0, 15.0, 15.0))));
}

#[test]
fn is_out_two_corners_disjoint() {
    assert!(!bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)).overlaps(&bx((20.0, 20.0, 20.0), (30.0, 30.0, 30.0))));
}

#[test]
fn is_out_two_corners_contained() {
    assert!(bx((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)).overlaps(&bx((2.0, 2.0, 2.0), (8.0, 8.0, 8.0))));
}
