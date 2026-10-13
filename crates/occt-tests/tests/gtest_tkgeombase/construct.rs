//! gce_MakeCirc2d_Test.cxx, gce_MakeCone_Test.cxx, gce_MakeCylinder_Test.cxx,
//! gce_MakeElips_Test.cxx, gce_MakeHypr_Test.cxx, GC_MakeArcOfCircle_Test.cxx,
//! GC_MakeCircle2d_Test.cxx, GC_MakeConicalSurface_Test.cxx, GC_MakeParabola2d_Test.cxx,
//! GC_MakePlane_Test.cxx, GC_MakeSegment2d_Test.cxx.

use std::f64::consts::PI;

use kentron::geometry::construct::ConstructionError;
use kentron::geometry::curve::{Circle, Ellipse, Hyperbola, Line};
use kentron::geometry::curve2::{Circle2, Hyperbola2, Line2, Parabola2};
use kentron::geometry::surface::{Cone, Cylinder, Plane};
use kentron::geometry::ParametricCurve;
use kentron::math::{Axis, Axis2, Direction, Direction2, Frame, Frame2, Point, Vector};
use kentron::tolerance::{ANGULAR, LINEAR, RESOLUTION};

use crate::{assert_near, assert_point_near, p, p2};

// ---- gce_MakeCirc2d ----

#[test]
fn circ2d_three_points_all_coincident_within_resolution_zero_radius() {
    let delta = 0.25 * RESOLUTION;
    let c = Circle2::through_points(p2(0.0, 0.0), p2(delta, 0.0), p2(0.0, delta)).unwrap();
    assert_near(c.radius(), 0.0, LINEAR);
}

#[test]
fn circ2d_three_points_one_near_coincident_pair() {
    let c = Circle2::through_points(p2(0.0, 0.0), p2(1.0e-150, 0.0), p2(0.0, 1.0)).unwrap();
    assert!(c.radius() > 0.0);
}

#[test]
fn circ2d_three_points_small_distinct() {
    let c = Circle2::through_points(p2(0.0, 0.0), p2(2.0e-150, 0.0), p2(0.0, 1.0)).unwrap();
    assert!(c.radius() > 0.0);
}

// ---- gce_MakeCone ----

fn base_cone() -> Cone {
    Cone::new(Frame::from_z(p(0.0, 0.0, 0.0), Direction::Z), PI / 6.0, 10.0)
}

#[track_caller]
fn assert_x_orthogonal_to_axis(frame: &Frame) {
    assert_near(frame.z().dot(frame.x()), 0.0, ANGULAR);
}

#[test]
fn cone_from_cone_and_distance() {
    let base = base_cone();
    let dist = 2.0;
    let cone = base.offset(dist).unwrap();
    assert_near(cone.radius(), base.radius() + dist / base.semi_angle().cos(), LINEAR);
    assert_near(cone.semi_angle(), base.semi_angle(), LINEAR);
    assert_point_near(cone.frame().origin(), base.frame().origin(), LINEAR);
    assert!(cone.frame().z().is_parallel(base.frame().z(), ANGULAR));
}

#[test]
fn cone_from_cone_and_distance_too_negative() {
    let base = base_cone();
    let dist = -(base.radius() * base.semi_angle().cos() + 1.0);
    assert_eq!(base.offset(dist), Err(ConstructionError::NegativeRadius));
}

#[test]
fn cone_from_cone_and_point() {
    let base = base_cone();
    let v = 4.0;
    let target = 12.0;
    let cone = base.parallel_through(p(target + v * base.semi_angle().tan(), 0.0, v)).unwrap();
    assert_near(cone.radius(), target, LINEAR);
}

#[test]
fn cone_from_cone_and_point_no_valid_parallel_cone() {
    assert_eq!(
        base_cone().parallel_through(p(1.0, 0.0, 20.0)),
        Err(ConstructionError::NegativeRadius)
    );
}

#[test]
fn cone_from_cone_and_point_two_candidates_chooses_nearest() {
    let base = base_cone();
    let point = p(5.0, 0.0, -20.0);
    let cone = base.parallel_through(point).unwrap();

    let radius = Line::from(base.frame().z_axis()).distance_to_point(point);
    let v = (point - base.frame().origin()).dot(Vector::from(base.frame().z()));
    let (tan, cos) = (base.semi_angle().tan(), base.semi_angle().cos());
    let cand1 = radius - v * tan;
    let cand2 = -radius - v * tan;
    let d1 = ((cand1 - base.radius()) * cos).abs();
    let d2 = ((cand2 - base.radius()) * cos).abs();
    assert_near(cone.radius(), if d1 <= d2 { cand1 } else { cand2 }, LINEAR);
}

#[test]
fn cone_from_two_points_and_radii_y_branch() {
    let (p1, p2_) = (p(0.0, 0.0, 0.0), p(0.0, 10.0, 10.0));
    let cone = Cone::from_points_radii(p1, p2_, 3.0, 2.0).unwrap();
    assert!(cone.frame().z().is_parallel(Direction::from_vector(p2_ - p1), ANGULAR));
    assert_x_orthogonal_to_axis(cone.frame());
}

#[test]
fn cone_from_two_points_and_radii_z_branch() {
    let (p1, p2_) = (p(0.0, 0.0, 0.0), p(0.0, 0.0, 10.0));
    let cone = Cone::from_points_radii(p1, p2_, 4.0, 2.0).unwrap();
    assert!(cone.frame().z().is_parallel(Direction::from_vector(p2_ - p1), ANGULAR));
    assert_x_orthogonal_to_axis(cone.frame());
}

#[test]
fn cone_from_four_points_y_branch() {
    let (p1, p2_) = (p(0.0, 0.0, 0.0), p(0.0, 10.0, 10.0));
    let cone = Cone::from_four_points(p1, p2_, p(3.0, 0.0, 0.0), p(5.0, 4.0, 4.0)).unwrap();
    assert!(cone.frame().z().is_parallel(Direction::from_vector(p2_ - p1), ANGULAR));
    assert_x_orthogonal_to_axis(cone.frame());
}

#[test]
fn cone_from_four_points_z_branch() {
    let (p1, p2_) = (p(0.0, 0.0, 0.0), p(0.0, 0.0, 10.0));
    let cone = Cone::from_four_points(p1, p2_, p(4.0, 0.0, 0.0), p(2.0, 0.0, 5.0)).unwrap();
    assert!(cone.frame().z().is_parallel(Direction::from_vector(p2_ - p1), ANGULAR));
    assert_x_orthogonal_to_axis(cone.frame());
}

// ---- gce_MakeCylinder ----

#[test]
fn cylinder_from_axis_y_branch() {
    let axis = Axis::new(p(0.0, 0.0, 0.0), Direction::new(0.0, 1.0, 1.0));
    let cyl = Cylinder::from_axis(axis, 5.0).unwrap();
    assert_near(cyl.radius(), 5.0, LINEAR);
    assert!(cyl.frame().z().is_parallel(axis.direction, ANGULAR));
    assert_x_orthogonal_to_axis(cyl.frame());
}

#[test]
fn cylinder_from_points_z_branch() {
    let (p1, p2_) = (p(0.0, 0.0, 0.0), p(0.0, 0.0, 10.0));
    let cyl = Cylinder::from_three_points(p1, p2_, p(3.0, 4.0, 0.0)).unwrap();
    assert_near(cyl.radius(), 5.0, LINEAR);
    assert!(cyl.frame().z().is_parallel(Direction::from_vector(p2_ - p1), ANGULAR));
    assert_x_orthogonal_to_axis(cyl.frame());
}

#[test]
fn cylinder_from_axis_negative_radius() {
    assert_eq!(
        Cylinder::from_axis(Axis::Z, -1.0),
        Err(ConstructionError::NegativeRadius)
    );
}

// ---- gce_MakeElips ----

#[test]
fn ellipse_from_points() {
    let e = Ellipse::from_points(p(10.0, 0.0, 0.0), p(0.0, 5.0, 0.0), p(0.0, 0.0, 0.0)).unwrap();
    assert_near(e.major_radius(), 10.0, LINEAR);
    assert_near(e.minor_radius(), 5.0, LINEAR);
}

#[test]
fn ellipse_from_points_collinear() {
    assert_eq!(
        Ellipse::from_points(p(10.0, 0.0, 0.0), p(5.0, 0.0, 0.0), p(0.0, 0.0, 0.0)),
        Err(ConstructionError::InvertAxis)
    );
}

#[test]
fn ellipse_from_points_tiny_scale() {
    assert_eq!(
        Ellipse::from_points(
            p(2.0 * RESOLUTION, 0.0, 0.0),
            p(0.0, 1.1 * RESOLUTION, 0.0),
            p(0.0, 0.0, 0.0)
        ),
        Err(ConstructionError::NullAxis)
    );
}

// ---- gce_MakeHypr ----

#[test]
fn hyperbola_from_frame_valid_radii() {
    let h = Hyperbola::try_new(Frame::WORLD, 10.0, 5.0).unwrap();
    assert_near(h.major_radius(), 10.0, LINEAR);
    assert_near(h.minor_radius(), 5.0, LINEAR);
}

#[test]
fn hyperbola_from_frame_major_less_than_minor() {
    let h = Hyperbola::try_new(Frame::WORLD, 3.0, 7.0).unwrap();
    assert_near(h.major_radius(), 3.0, LINEAR);
    assert_near(h.minor_radius(), 7.0, LINEAR);
}

#[test]
fn hyperbola_from_frame_equal_and_zero_radii() {
    assert!(Hyperbola::try_new(Frame::WORLD, 5.0, 5.0).is_ok());
    assert!(Hyperbola::try_new(Frame::WORLD, 0.0, 0.0).is_ok());
}

#[test]
fn hyperbola_from_frame_negative_radii() {
    for (a, b) in [(-1.0, 5.0), (5.0, -1.0), (-3.0, -2.0)] {
        assert_eq!(Hyperbola::try_new(Frame::WORLD, a, b), Err(ConstructionError::NegativeRadius));
    }
}

#[test]
fn hyperbola_from_points() {
    let h = Hyperbola::from_points(p(10.0, 0.0, 0.0), p(0.0, 5.0, 0.0), p(0.0, 0.0, 0.0)).unwrap();
    assert_near(h.major_radius(), 10.0, LINEAR);
    assert_near(h.minor_radius(), 5.0, LINEAR);
}

#[test]
fn hyperbola_from_points_major_less_than_minor() {
    let h = Hyperbola::from_points(p(3.0, 0.0, 0.0), p(0.0, 7.0, 0.0), p(0.0, 0.0, 0.0)).unwrap();
    assert_near(h.major_radius(), 3.0, LINEAR);
    assert_near(h.minor_radius(), 7.0, LINEAR);
}

#[test]
fn hyperbola_from_points_confused() {
    let cases: [(Point, Point, Point); 3] = [
        (p(1.0, 2.0, 3.0), p(4.0, 5.0, 6.0), p(1.0, 2.0, 3.0)),
        (p(4.0, 5.0, 6.0), p(1.0, 2.0, 3.0), p(1.0, 2.0, 3.0)),
        (p(5.0, 5.0, 5.0), p(5.0, 5.0, 5.0), p(0.0, 0.0, 0.0)),
    ];
    for (s1, s2, center) in cases {
        assert_eq!(Hyperbola::from_points(s1, s2, center), Err(ConstructionError::ConfusedPoints));
    }
}

#[test]
fn hyperbola_from_points_collinear() {
    assert_eq!(
        Hyperbola::from_points(p(10.0, 0.0, 0.0), p(5.0, 0.0, 0.0), p(0.0, 0.0, 0.0)),
        Err(ConstructionError::ColinearPoints)
    );
}

#[test]
fn hyperbola_from_points_tiny_scale() {
    assert_eq!(
        Hyperbola::from_points(
            p(2.0 * RESOLUTION, 0.0, 0.0),
            p(0.0, 1.1 * RESOLUTION, 0.0),
            p(0.0, 0.0, 0.0)
        ),
        Err(ConstructionError::ConfusedPoints)
    );
}

#[test]
fn hyperbola2_from_major_axis() {
    let h = Hyperbola2::try_new(Frame2::from_axis(Axis2::X, true), 10.0, 5.0).unwrap();
    assert_near(h.major_radius(), 10.0, LINEAR);
    assert_near(h.minor_radius(), 5.0, LINEAR);
    let h = Hyperbola2::try_new(Frame2::from_axis(Axis2::X, true), 3.0, 7.0).unwrap();
    assert_near(h.major_radius(), 3.0, LINEAR);
    assert_near(h.minor_radius(), 7.0, LINEAR);
}

#[test]
fn hyperbola2_negative_radius() {
    assert_eq!(
        Hyperbola2::try_new(Frame2::WORLD, -1.0, 5.0),
        Err(ConstructionError::NegativeRadius)
    );
    assert_eq!(
        Hyperbola2::try_new(Frame2::WORLD, 5.0, -2.0),
        Err(ConstructionError::NegativeRadius)
    );
    assert!(Hyperbola2::try_new(Frame2::WORLD, 8.0, 4.0).is_ok());
}

#[test]
fn hyperbola2_from_points() {
    let h = Hyperbola2::from_points(p2(10.0, 0.0), p2(0.0, 5.0), p2(0.0, 0.0)).unwrap();
    assert_near(h.major_radius(), 10.0, LINEAR);
    assert_near(h.minor_radius(), 5.0, LINEAR);
    let h = Hyperbola2::from_points(p2(3.0, 0.0), p2(0.0, 7.0), p2(0.0, 0.0)).unwrap();
    assert_near(h.major_radius(), 3.0, LINEAR);
    assert_near(h.minor_radius(), 7.0, LINEAR);
}

#[test]
fn hyperbola2_from_points_degenerate() {
    assert_eq!(
        Hyperbola2::from_points(p2(1.0, 2.0), p2(4.0, 5.0), p2(1.0, 2.0)),
        Err(ConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Hyperbola2::from_points(p2(4.0, 5.0), p2(1.0, 2.0), p2(1.0, 2.0)),
        Err(ConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Hyperbola2::from_points(p2(10.0, 0.0), p2(5.0, 0.0), p2(0.0, 0.0)),
        Err(ConstructionError::ColinearPoints)
    );
}

// ---- GC_MakeArcOfCircle ----

#[test]
fn arc_from_points_and_tangent() {
    assert!(Circle::arc_tangent(p(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), p(0.0, 1.0, 0.0))
        .is_ok());
}

#[test]
fn arc_from_three_points() {
    assert!(Circle::arc_through(p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(-1.0, 0.0, 0.0)).is_ok());
}

#[test]
fn arc_from_points_and_tangent_confusion_scale() {
    let scale = 2.0 * LINEAR;
    let p1 = p(scale, 0.0, 0.0);
    let arc = Circle::arc_tangent(p1, Vector::new(0.0, 1.0, 0.0), p(0.0, scale, 0.0)).unwrap();
    assert_point_near(arc.start_point(), p1, LINEAR);
}

// ---- GC_MakeCircle2d ----

fn base_circle2() -> Circle2 {
    Circle2::new(Frame2::from_axis(Axis2::X, true), 5.0)
}

#[test]
fn circle2_concentric_through_point() {
    assert_near(base_circle2().concentric_through(p2(0.0, 10.0)).radius(), 10.0, LINEAR);
}

#[test]
fn circle2_concentric_offset() {
    assert_near(base_circle2().concentric_offset(-2.0).unwrap().radius(), 3.0, LINEAR);
}

#[test]
fn circle2_concentric_through_center() {
    assert_near(base_circle2().concentric_through(p2(0.0, 0.0)).radius(), 0.0, LINEAR);
}

// ---- GC_MakeConicalSurface ----

#[test]
fn cone_angle_below_resolution() {
    assert_eq!(Cone::try_new(Frame::WORLD, 0.5 * RESOLUTION, 2.0), Err(ConstructionError::BadAngle));
}

#[test]
fn cone_valid_angle() {
    assert!(Cone::try_new(Frame::WORLD, PI / 6.0, 2.0).is_ok());
}

// ---- GC_MakeParabola2d ----

#[track_caller]
fn check_parabola2(
    directrix: Axis2,
    focus: (f64, f64),
    ccw: bool,
    focal: f64,
    apex: (f64, f64),
    parameter: f64,
    coefficients: [f64; 6],
) {
    const TOL: f64 = 1.0e-12;
    let parab = Parabola2::from_directrix_focus(directrix, p2(focus.0, focus.1), ccw);
    assert_near(parab.focal(), focal, TOL);
    assert_near(parab.apex().x, apex.0, TOL);
    assert_near(parab.apex().y, apex.1, TOL);
    assert_near(parab.semi_latus_rectum(), parameter, TOL);
    for (i, (got, expected)) in parab.implicit_coefficients().into_iter().zip(coefficients).enumerate()
    {
        assert!((got - expected).abs() <= TOL, "coefficient {i}: {got} != {expected}");
    }
}

/// OCC26747_1.
#[test]
fn parabola2_opening_right() {
    let directrix = Axis2::new(p2(0.0, 3.0), Direction2::Y);
    check_parabola2(directrix, (1.0, 3.0), true, 0.5, (0.5, 3.0), 1.0, [0.0, 1.0, 0.0, -1.0, -3.0, 10.0]);
}

/// OCC26747_2.
#[test]
fn parabola2_opening_left() {
    let directrix = Axis2::new(p2(0.0, 0.0), Direction2::Y);
    check_parabola2(directrix, (-1.0, 3.0), false, 0.5, (-0.5, 3.0), 1.0, [0.0, 1.0, 0.0, 1.0, -3.0, 10.0]);
}

/// OCC26747_3.
#[test]
fn parabola2_degenerate() {
    let directrix = Axis2::new(p2(0.0, 0.0), Direction2::Y);
    check_parabola2(directrix, (0.0, 3.0), false, 0.0, (0.0, 3.0), 0.0, [0.0, 1.0, 0.0, 0.0, -3.0, 9.0]);
}

// ---- GC_MakePlane ----

#[test]
fn plane_from_plane_and_distance() {
    assert_near(Plane::XY.offset(5.0).origin().z, 5.0, LINEAR);
}

#[test]
fn plane_from_plane_and_point() {
    let point = p(2.0, -3.0, 7.0);
    assert_point_near(Plane::XY.parallel_through(point).origin(), point, LINEAR);
}

#[test]
fn plane_from_axis() {
    let axis = Axis::new(p(1.0, 2.0, 3.0), Direction::Y);
    assert!(Plane::from_axis(axis).normal().is_parallel(axis.direction, ANGULAR));
}

// ---- GC_MakeSegment2d ----

#[test]
fn segment2_from_points() {
    assert_eq!(
        Line2::segment(p2(0.0, 0.0), p2(0.5 * RESOLUTION, 0.0)).err(),
        Some(ConstructionError::ConfusedPoints)
    );
    assert!(Line2::segment(p2(0.0, 0.0), p2(1.0e-150, 0.0)).is_ok());
    assert_eq!(
        Line2::segment(p2(0.0, 0.0), p2(RESOLUTION, 0.0)).err(),
        Some(ConstructionError::ConfusedPoints)
    );
}

#[test]
fn segment2_from_point_direction() {
    let seg = |x: f64, y: f64| Line2::segment_to_projection(p2(0.0, 0.0), Direction2::X, p2(x, y));
    assert_eq!(seg(0.5 * RESOLUTION, 3.0).err(), Some(ConstructionError::ConfusedPoints));
    assert!(seg(2.0 * RESOLUTION, 3.0).is_ok());
    assert_eq!(seg(RESOLUTION, 5.0).err(), Some(ConstructionError::ConfusedPoints));
}
