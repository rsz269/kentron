//! BndLib_Test.cxx and the GeomBndLib_*_Test.cxx files.
//!
//! OCCT has two bounding-box implementations (BndLib and the newer GeomBndLib); kentron has
//! one (`geometry::bounds`). The GeomBndLib tests that compare the two
//! (`*_CompareWithBndLib`, `AdaptorConstructor_*`, `TrimmedCurveHandle_*`,
//! `*_UsesBasisSpecialization`, `SamplingHelpers*`) are not ported; those with concrete
//! expectations are, under a `geom_bnd_lib_` prefix. BndLib_Add3dCurve/Add2dCurve/AddSurface
//! cases identical to a `BndLib::Add` case (full circles, ellipses, spheres, the cylinder and
//! cone patches) are covered by that case, since kentron has one entry point for both.

use std::f64::consts::PI;

use kentron::geometry::bounds;
use kentron::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola,
    TrimmedCurve,
};
use kentron::geometry::curve2::{
    BSplineCurve2, BezierCurve2, Circle2, Curve2, Ellipse2, Hyperbola2, Line2, Parabola2,
    TrimmedCurve2,
};
use kentron::geometry::surface::{
    BezierSurface, Cone, Cylinder, ExtrusionSurface, OffsetSurface, Plane, RevolutionSurface,
    Sphere, Surface, Torus,
};
use kentron::geometry::ParametricSurface;
use kentron::math::{Axis, BoundingBox, BoundingBox2, Direction, Direction2, Frame, Frame2, Point};
use kentron::tolerance::{INFINITE, LINEAR};

use crate::{assert_near, p, p2};

const INF: f64 = f64::INFINITY;

/// Checks `[xmin, xmax, ymin, ymax, zmin, zmax]`.
#[track_caller]
fn assert_box(b: &BoundingBox, expected: [f64; 6], tolerance: f64) {
    let actual = [b.min.x, b.max.x, b.min.y, b.max.y, b.min.z, b.max.z];
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= tolerance, "component {i}: {a} != {e} in {b:?}");
    }
}

/// Checks `[xmin, xmax, ymin, ymax]`.
#[track_caller]
fn assert_box2(b: &BoundingBox2, expected: [f64; 4], tolerance: f64) {
    let actual = [b.min.x, b.max.x, b.min.y, b.max.y];
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= tolerance, "component {i}: {a} != {e} in {b:?}");
    }
}

/// `inner` lies within `outer` (both enlarged by `tolerance`).
#[track_caller]
fn assert_within(inner: &BoundingBox, outer: &BoundingBox, tolerance: f64) {
    assert!(inner.min.x >= outer.min.x - tolerance && inner.max.x <= outer.max.x + tolerance);
    assert!(inner.min.y >= outer.min.y - tolerance && inner.max.y <= outer.max.y + tolerance);
    assert!(inner.min.z >= outer.min.z - tolerance && inner.max.z <= outer.max.z + tolerance);
}

#[track_caller]
fn assert_contains_surface(b: &BoundingBox, s: &Surface, (u1, u2, v1, v2): (f64, f64, f64, f64), nu: usize, nv: usize) {
    for i in 0..=nu {
        let u = u1 + (u2 - u1) * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = v1 + (v2 - v1) * j as f64 / nv as f64;
            let q = s.point(u, v);
            assert!(b.enlarged(LINEAR).contains(q), "{q:?} at ({u}, {v}) outside {b:?}");
        }
    }
}

fn line(origin: Point, direction: Direction) -> Curve {
    Curve::Line(Line::new(origin, direction))
}

fn frame_z(origin: Point) -> Frame {
    Frame::from_z(origin, Direction::Z)
}

fn circle(frame: Frame, radius: f64) -> Curve {
    Curve::Circle(Circle::new(frame, radius))
}

// ---- BndLib: lines ----

#[test]
fn line_finite_segment() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), Direction::X), 0., 10., 0.);
    assert_box(&b, [0., 10., 0., 0., 0., 0.], LINEAR);
}

#[test]
fn line_positive_direction_negative_infinite() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), Direction::X), -INF, 5., 0.);
    assert_eq!(b.min.x, -INF, "xmin should be open");
    assert!(b.max.x.is_finite(), "xmax should not be open");
}

#[test]
fn line_negative_direction_negative_infinite() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), -Direction::X), -INF, 5., 0.);
    assert_eq!(b.max.x, INF, "xmax should be open");
    assert!(b.min.x.is_finite(), "xmin should not be open");
}

#[test]
fn line_negative_direction_positive_infinite() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), -Direction::X), -5., INF, 0.);
    assert_eq!(b.min.x, -INF, "xmin should be open");
    assert!(b.max.x.is_finite(), "xmax should not be open");
}

#[test]
fn line_diagonal_direction_both_infinite() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), Direction::new(-1., -1., -1.)), -INF, INF, 0.);
    assert_eq!(b, BoundingBox::WHOLE);
}

// ---- BndLib: conics ----

#[test]
fn circle_full() {
    let c = circle(Frame::WORLD, 5.0);
    assert_box(&bounds::curve(&c, 0.), [-5., 5., -5., 5., 0., 0.], LINEAR);
    assert_box(&bounds::curve_between(&c, 0., 2. * PI, 0.), [-5., 5., -5., 5., 0., 0.], LINEAR);
}

#[test]
fn circle_arc_first_quadrant() {
    let b = bounds::curve_between(&circle(Frame::WORLD, 5.0), 0., PI / 2., 0.);
    assert_box(&b, [0., 5., 0., 5., 0., 0.], LINEAR);
}

#[test]
fn circle_rotated_axis() {
    let b = bounds::curve(&circle(Frame::from_z(p(0., 0., 0.), Direction::Y), 3.0), 0.);
    assert_box(&b, [-3., 3., 0., 0., -3., 3.], LINEAR);
}

#[test]
fn ellipse_full() {
    let e = Curve::Ellipse(Ellipse::new(Frame::WORLD, 10., 5.));
    assert_box(&bounds::curve(&e, 0.), [-10., 10., -5., 5., 0., 0.], LINEAR);
}

#[test]
fn ellipse_arc() {
    let e = Curve::Ellipse(Ellipse::new(Frame::WORLD, 10., 5.));
    let b = bounds::curve_between(&e, 0., PI / 2., 0.);
    assert_near(b.min.x, 0., LINEAR);
    assert_near(b.max.x, 10., LINEAR);
    assert_near(b.min.y, 0., LINEAR);
    assert_near(b.max.y, 5., LINEAR);
}

#[test]
fn hyperbola_simple_arc() {
    let h = Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 3., 2.));
    let b = bounds::curve_between(&h, -1., 1., 0.);
    assert_near(b.min.x, 3., LINEAR);
    assert_near(b.max.x, 3. * 1f64.cosh(), LINEAR);
    assert_near(b.min.y, 2. * (-1f64).sinh(), LINEAR);
    assert_near(b.max.y, 2. * 1f64.sinh(), LINEAR);
}

#[test]
fn hyperbola_rotated_multiple_extrema() {
    let frame = Frame::new(p(0., 0., 0.), Direction::new(1., 1., 1.), Direction::new(1., -1., 0.));
    let b = bounds::curve_between(&Curve::Hyperbola(Hyperbola::new(frame, 5., 3.)), -2., 2., 0.);
    assert!(!b.is_empty());
    assert!(b.min.x < b.max.x);
    assert!(b.min.y < b.max.y);
}

#[test]
fn hyperbola_infinite_parameter() {
    let h = Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 3., 2.));
    let b = bounds::curve_between(&h, 0., INF, 0.);
    assert_eq!(b.max.x, INF);
    assert_eq!(b.max.y, INF);
}

#[test]
fn parabola_finite_arc() {
    let b = bounds::curve_between(&Curve::Parabola(Parabola::new(Frame::WORLD, 2.)), -4., 4., 0.);
    assert_near(b.min.x, 0., LINEAR);
    assert_near(b.max.x, 2., LINEAR);
    assert_near(b.min.y, -4., LINEAR);
    assert_near(b.max.y, 4., LINEAR);
}

#[test]
fn parabola_infinite_parameter() {
    let b = bounds::curve_between(&Curve::Parabola(Parabola::new(Frame::WORLD, 2.)), 0., INF, 0.);
    assert_eq!(b.max.x, INF);
    assert_eq!(b.max.y, INF);
}

// ---- BndLib: elementary surfaces ----

fn sphere(frame: Frame, radius: f64) -> Surface {
    Surface::Sphere(Sphere::new(frame, radius))
}

#[test]
fn sphere_full() {
    assert_box(&bounds::surface(&sphere(Frame::WORLD, 7.), 0.), [-7., 7., -7., 7., -7., 7.], LINEAR);
}

#[test]
fn sphere_patch() {
    let b = bounds::surface_within(&sphere(Frame::WORLD, 5.), 0., 2. * PI, 0., PI / 2., 0.);
    assert_box(&b, [-5., 5., -5., 5., 0., 5.], LINEAR);
}

#[test]
fn sphere_translated() {
    let b = bounds::surface(&sphere(frame_z(p(10., 20., 30.)), 3.), 0.);
    assert_box(&b, [7., 13., 17., 23., 27., 33.], LINEAR);
}

fn cylinder(frame: Frame, radius: f64) -> Surface {
    Surface::Cylinder(Cylinder::new(frame, radius))
}

#[test]
fn cylinder_finite_patch() {
    let b = bounds::surface_within(&cylinder(Frame::WORLD, 4.), 0., 2. * PI, 0., 10., 0.);
    assert_box(&b, [-4., 4., -4., 4., 0., 10.], LINEAR);
}

#[test]
fn cylinder_infinite_v() {
    let b = bounds::surface_within(&cylinder(Frame::WORLD, 4.), 0., 2. * PI, 0., INF, 0.);
    assert_eq!(b.max.z, INF, "zmax should be open");
    assert!(b.min.z.is_finite(), "zmin should not be open");
}

#[test]
fn cylinder_partial_u() {
    let b = bounds::surface_within(&cylinder(Frame::WORLD, 5.), 0., PI / 2., 0., 10., 0.);
    assert_box(&b, [0., 5., 0., 5., 0., 10.], LINEAR);
}

fn cone(semi_angle: f64, radius: f64) -> Surface {
    Surface::Cone(Cone::new(Frame::WORLD, semi_angle, radius))
}

#[test]
fn cone_finite_patch() {
    let b = bounds::surface_within(&cone(PI / 4., 2.), 0., 2. * PI, 0., 5., 0.);
    let max_radius = 2. + 5. * (PI / 4.).sin();
    assert_near(b.max.x, max_radius, LINEAR);
    assert_near(b.max.y, max_radius, LINEAR);
    assert_near(b.max.z, 5. * (PI / 4.).cos(), LINEAR);
    assert_near(b.min.z, 0., LINEAR);
}

#[test]
fn cone_negative_v() {
    let b = bounds::surface_within(&cone(PI / 6., 3.), 0., 2. * PI, -5., 0., 0.);
    assert_near(b.min.z, -5. * (PI / 6.).cos(), LINEAR);
    assert_near(b.max.z, 0., LINEAR);
}

#[test]
fn cone_infinite_v() {
    let b = bounds::surface_within(&cone(PI / 4., 2.), 0., 2. * PI, 0., INF, 0.);
    assert_eq!(b.max.z, INF);
}

fn torus(frame: Frame, major: f64, minor: f64) -> Surface {
    Surface::Torus(Torus::new(frame, major, minor))
}

#[test]
fn torus_full() {
    let b = bounds::surface(&torus(Frame::WORLD, 10., 3.), 0.);
    assert_box(&b, [-13., 13., -13., 13., -3., 3.], LINEAR);
}

#[test]
fn torus_patch() {
    let b = bounds::surface_within(&torus(Frame::WORLD, 10., 3.), 0., PI / 2., 0., PI, 0.);
    assert!(!b.is_empty());
    assert!(b.min.x >= -LINEAR, "xmin should be >= 0 in the first quadrant");
    assert!(b.min.y >= -LINEAR, "ymin should be >= 0 in the first quadrant");
    assert!(b.max.z <= 3. + LINEAR, "zmax should not exceed r");
}

#[test]
fn torus_negative_v_parameter() {
    let b = bounds::surface_within(&torus(Frame::WORLD, 10., 3.), 0., 2. * PI, -PI, 0., 0.);
    assert!(!b.is_empty());
    assert!(b.min.z <= 0.);
    assert!(b.min.z >= -3. - LINEAR);
}

#[test]
fn torus_large_negative_v_parameter() {
    let b = bounds::surface_within(&torus(Frame::WORLD, 10., 3.), 0., 2. * PI, -3. * PI, -2. * PI, 0.);
    assert!(!b.is_empty());
    assert!(b.min.x <= -13. + LINEAR && b.max.x >= 13. - LINEAR);
    assert!(b.min.y <= -13. + LINEAR && b.max.y >= 13. - LINEAR);
}

#[test]
fn torus_translated() {
    let b = bounds::surface(&torus(frame_z(p(5., 5., 5.)), 8., 2.), 0.);
    assert_box(&b, [-5., 15., -5., 15., 3., 7.], LINEAR);
}

// ---- BndLib: 2D conics ----

#[test]
fn line2_finite_segment() {
    let l = Curve2::Line(Line2::new(p2(0., 0.), Direction2::new(1., 1.)));
    let e = 10. / 2f64.sqrt();
    assert_box2(&bounds::curve2_between(&l, 0., 10., 0.), [0., e, 0., e], LINEAR);
}

#[test]
fn circle2_full() {
    let c = Curve2::Circle(Circle2::new(Frame2::WORLD, 5.));
    assert_box2(&bounds::curve2(&c, 0.), [-5., 5., -5., 5.], LINEAR);
}

#[test]
fn ellipse2_full() {
    let e = Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 8., 4.));
    assert_box2(&bounds::curve2(&e, 0.), [-8., 8., -4., 4.], LINEAR);
}

#[test]
fn hyperbola2_simple_arc() {
    let h = Curve2::Hyperbola(Hyperbola2::new(Frame2::WORLD, 3., 2.));
    let b = bounds::curve2_between(&h, -1., 1., 0.);
    assert_box2(&b, [3., 3. * 1f64.cosh(), 2. * (-1f64).sinh(), 2. * 1f64.sinh()], LINEAR);
}

#[test]
fn parabola2_finite_arc() {
    let pb = Curve2::Parabola(Parabola2::new(Frame2::WORLD, 2.));
    assert_box2(&bounds::curve2_between(&pb, -4., 4., 0.), [0., 2., -4., 4.], LINEAR);
}

// ---- BndLib: tolerance ----

#[test]
fn circle_with_tolerance() {
    let b = bounds::curve(&circle(Frame::WORLD, 5.), 0.5);
    assert_box(&b, [-5.5, 5.5, -5.5, 5.5, -0.5, 0.5], LINEAR);
}

#[test]
fn sphere_with_tolerance() {
    let b = bounds::surface(&sphere(Frame::WORLD, 3.), 1.);
    assert_box(&b, [-4., 4., -4., 4., -4., 4.], LINEAR);
}

// ---- BndLib_Add3dCurve ----

#[test]
fn add3d_circle_arc() {
    let b = bounds::curve_between(&circle(Frame::WORLD, 5.), 0., PI / 2., 0.);
    assert_near(b.min.x, 0., LINEAR);
    assert_near(b.max.x, 5., LINEAR);
    assert_near(b.min.y, 0., LINEAR);
    assert_near(b.max.y, 5., LINEAR);
}

#[test]
fn add3d_line_segment() {
    let segment = Curve::Trimmed(TrimmedCurve::new(line(p(0., 0., 0.), Direction::X), 0., 10.));
    assert_box(&bounds::curve(&segment, 0.), [0., 10., 0., 0., 0., 0.], LINEAR);
}

#[test]
fn add3d_bezier_curve() {
    let bezier = Curve::Bezier(BezierCurve::new(vec![
        p(0., 0., 0.),
        p(1., 2., 0.),
        p(3., 2., 0.),
        p(4., 0., 0.),
    ]));
    let b = bounds::curve(&bezier, 0.);
    assert!(!b.is_empty());
    assert!(b.min.x <= 0. && b.max.x >= 4.);
    assert!(b.min.y <= 0. && b.max.y > 1., "the curve rises above y = 0");
}

#[test]
fn add3d_bspline_curve() {
    let bspline = Curve::BSpline(BSplineCurve::new(
        vec![p(0., 0., 0.), p(1., 1., 0.), p(2., 1., 0.), p(3., 0., 0.)],
        vec![0., 1.],
        vec![4, 4],
        3,
    ));
    let b = bounds::curve(&bspline, 0.);
    assert!(!b.is_empty());
    assert!(b.min.x <= 0. && b.max.x >= 3.);
}

#[test]
fn add3d_tight_circle() {
    let b = bounds::curve_tight(&circle(Frame::WORLD, 5.), 0.);
    assert_near(b.min.x, -5., LINEAR);
    assert_near(b.max.x, 5., LINEAR);
    assert_near(b.min.y, -5., LINEAR);
    assert_near(b.max.y, 5., LINEAR);
}

#[test]
fn add3d_circle_rotated() {
    let b = bounds::curve(&circle(Frame::from_z(p(0., 0., 0.), Direction::X), 4.), 0.);
    assert_box(&b, [0., 0., -4., 4., -4., 4.], LINEAR);
}

// ---- BndLib_Add2dCurve ----

#[test]
fn add2d_circle_arc() {
    let c = Curve2::Circle(Circle2::new(Frame2::WORLD, 5.));
    assert_box2(&bounds::curve2_between(&c, 0., PI / 2., 0.), [0., 5., 0., 5.], LINEAR);
}

#[test]
fn add2d_line_segment() {
    let l = Curve2::Line(Line2::new(p2(0., 0.), Direction2::new(2., 1.)));
    let segment = Curve2::Trimmed(TrimmedCurve2::new(l, 0., 125f64.sqrt()));
    let b = bounds::curve2(&segment, 0.);
    assert_near(b.min.x, 0., LINEAR);
    assert_near(b.max.x, 10., 0.01);
    assert_near(b.min.y, 0., LINEAR);
    assert_near(b.max.y, 5., 0.01);
}

#[test]
fn add2d_bezier_curve() {
    let bezier = Curve2::Bezier(BezierCurve2::new(vec![p2(0., 0.), p2(1., 3.), p2(3., 3.), p2(4., 0.)]));
    let b = bounds::curve2(&bezier, 0.);
    assert!(b.min.x <= 0. && b.max.x >= 4.);
    assert!(b.min.y <= 0. && b.max.y > 1., "the curve rises above y = 0");
}

#[test]
fn add2d_tight_ellipse() {
    let e = Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 6., 3.));
    assert_box2(&bounds::curve2_between_tight(&e, 0., 2. * PI, 0.), [-6., 6., -3., 3.], LINEAR);
}

#[test]
fn add2d_offset_circle_center() {
    let frame = Frame2::new(p2(5., 5.), Direction2::X, true);
    let c = Curve2::Circle(Circle2::new(frame, 3.));
    assert_box2(&bounds::curve2(&c, 0.), [2., 8., 2., 8.], LINEAR);
}

// ---- BndLib_AddSurface ----

#[test]
fn add_surface_plane() {
    let b = bounds::surface_within(&Surface::Plane(Plane::XY), -5., 5., -5., 5., 0.);
    assert_box(&b, [-5., 5., -5., 5., 0., 0.], LINEAR);
}

#[test]
fn add_surface_torus_full() {
    let b = bounds::surface(&torus(Frame::WORLD, 10., 3.), 0.);
    assert!(b.min.x <= -13. + LINEAR && b.max.x >= 13. - LINEAR);
    assert!(b.min.y <= -13. + LINEAR && b.max.y >= 13. - LINEAR);
    assert_near(b.min.z, -3., LINEAR);
    assert_near(b.max.z, 3., LINEAR);
}

#[test]
fn add_surface_bezier() {
    let s = Surface::Bezier(BezierSurface::new(vec![
        vec![p(0., 0., 0.), p(0., 10., 2.)],
        vec![p(10., 0., 2.), p(10., 10., 0.)],
    ]));
    let b = bounds::surface(&s, 0.);
    assert!(b.min.x <= 0. && b.max.x >= 10.);
    assert!(b.min.y <= 0. && b.max.y >= 10.);
    assert!(b.min.z <= 0. && b.max.z >= 2.);
}

#[test]
fn add_surface_tight_sphere() {
    let b = bounds::surface_tight(&sphere(Frame::WORLD, 4.), 0.);
    assert_box(&b, [-4., 4., -4., 4., -4., 4.], LINEAR);
}

#[test]
fn add_surface_cylinder_translated() {
    let b = bounds::surface_within(&cylinder(frame_z(p(5., 5., 0.)), 3.), 0., 2. * PI, 0., 10., 0.);
    assert_box(&b, [2., 8., 2., 8., 0., 10.], LINEAR);
}

#[test]
fn add_surface_plane_tilted() {
    let plane = Surface::Plane(Plane::new(p(0., 0., 0.), Direction::new(0., 1., 1.)));
    assert!(!bounds::surface_within(&plane, 0., 1., 0., 1., 0.).is_empty());
}

#[test]
fn add_surface_large_parameters() {
    let poles = (0..3).map(|i| (0..3).map(|j| p(i as f64 * 10., j as f64 * 10., 0.)).collect()).collect();
    let b = bounds::surface_within(&Surface::Bezier(BezierSurface::new(poles)), 0., 1., 0., 1., 0.);
    assert!(b.min.x <= 1. && b.max.x >= 19.);
    assert!(b.min.y <= 1. && b.max.y >= 19.);
}

#[test]
fn add_surface_large_parameter_values_keep_precision() {
    let (offset, range) = (1.0e10, 100.);
    let b = bounds::surface_within(&Surface::Plane(Plane::XY), offset, offset + range, offset, offset + range, 0.);
    assert_near(b.min.x, offset, 1.);
    assert_near(b.max.x, offset + range, 1.);
    assert_near(b.min.y, offset, 1.);
    assert_near(b.max.y, offset + range, 1.);
    assert!(b.max.x - b.min.x > 0.9 * range);
    assert!(b.max.y - b.min.y > 0.9 * range);
}

// ---- GeomBndLib (concrete expectations) ----

#[test]
fn geom_bnd_lib_line_finite_segment() {
    let b = bounds::curve_between(&line(p(0., 0., 0.), Direction::new(1., 1., 0.)), 0., 10., 0.);
    let len = 10. / 2f64.sqrt();
    assert_box(&b, [0., len, 0., len, 0., 0.], LINEAR);
}

#[test]
fn geom_bnd_lib_circle_arc() {
    let b = bounds::curve_between(&circle(Frame::WORLD, 10.), 0., PI / 2., 0.);
    assert_box(&b, [0., 10., 0., 10., 0., 0.], LINEAR);
}

#[test]
fn geom_bnd_lib_hyperbola_finite_arc() {
    let b = bounds::curve_between(&Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 5., 3.)), -1., 1., 0.);
    assert!(b.min.x <= 5. + LINEAR && b.max.x >= 5. * 1f64.cosh() - LINEAR);
    assert!(b.min.y <= -3. * 1f64.sinh() + LINEAR && b.max.y >= 3. * 1f64.sinh() - LINEAR);
}

#[test]
fn geom_bnd_lib_parabola_finite_arc() {
    let b = bounds::curve_between(&Curve::Parabola(Parabola::new(Frame::WORLD, 2.)), -5., 5., 0.);
    assert!(b.min.x <= LINEAR && b.max.x >= 25. / 8. - LINEAR);
    assert!(b.min.y <= -5. + LINEAR && b.max.y >= 5. - LINEAR);
}

#[test]
fn geom_bnd_lib_bezier_and_bspline() {
    let bezier = Curve::Bezier(BezierCurve::new(vec![p(0., 0., 0.), p(5., 10., 0.), p(10., 0., 0.)]));
    let bspline = Curve::BSpline(BSplineCurve::new(
        vec![p(0., 0., 0.), p(2., 5., 0.), p(8., 5., 0.), p(10., 0., 0.)],
        vec![0., 0.5, 1.],
        vec![3, 1, 3],
        2,
    ));
    for c in [bezier, bspline] {
        let b = bounds::curve(&c, LINEAR);
        assert!(b.min.x <= LINEAR && b.max.x >= 10. - LINEAR);
        assert!(b.min.y <= LINEAR && b.max.y >= 0.);
        assert_near(b.min.z, 0., LINEAR);
        assert_near(b.max.z, 0., LINEAR);
    }
}

#[test]
fn geom_bnd_lib_tight_circle_within_fast() {
    let c = circle(Frame::WORLD, 5.);
    assert_within(&bounds::curve_tight(&c, LINEAR), &bounds::curve(&c, LINEAR), LINEAR);
}

#[test]
fn geom_bnd_lib_2d_line_and_conics() {
    let l = Curve2::Line(Line2::new(p2(0., 0.), Direction2::new(1., 1.)));
    let e = 10. / 2f64.sqrt();
    assert_box2(&bounds::curve2_between(&l, 0., 10., 0.), [0., e, 0., e], LINEAR);
    let c = Curve2::Circle(Circle2::new(Frame2::WORLD, 10.));
    assert_box2(&bounds::curve2_between(&c, 0., PI / 2., 0.), [0., 10., 0., 10.], LINEAR);
    let el = Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, 10., 5.));
    assert_box2(&bounds::curve2(&el, 0.), [-10., 10., -5., 5.], LINEAR);
    let c5 = Curve2::Circle(Circle2::new(Frame2::WORLD, 5.));
    assert_box2(&bounds::curve2_tight(&c5, 0.), [-5., 5., -5., 5.], LINEAR);
}

#[test]
fn geom_bnd_lib_2d_hyperbola_and_parabola() {
    let h = Curve2::Hyperbola(Hyperbola2::new(Frame2::WORLD, 5., 3.));
    let b = bounds::curve2_between(&h, -1., 1., 0.);
    assert!(b.min.x <= 5. + LINEAR && b.max.x >= 5. * 1f64.cosh() - LINEAR);
    assert!(b.min.y <= -3. * 1f64.sinh() + LINEAR && b.max.y >= 3. * 1f64.sinh() - LINEAR);
    let pb = Curve2::Parabola(Parabola2::new(Frame2::WORLD, 2.));
    let b = bounds::curve2_between(&pb, -5., 5., 0.);
    assert!(b.min.x <= LINEAR && b.max.x >= 25. / 8. - LINEAR);
    assert!(b.min.y <= -5. + LINEAR && b.max.y >= 5. - LINEAR);
}

#[test]
fn geom_bnd_lib_2d_bezier_and_bspline() {
    let bezier = Curve2::Bezier(BezierCurve2::new(vec![p2(0., 0.), p2(5., 10.), p2(10., 0.)]));
    let bspline = Curve2::BSpline(BSplineCurve2::new(
        vec![p2(0., 0.), p2(2., 5.), p2(8., 5.), p2(10., 0.)],
        vec![0., 0.5, 1.],
        vec![3, 1, 3],
        2,
    ));
    for c in [bezier, bspline] {
        let b = bounds::curve2(&c, LINEAR);
        assert!(b.min.x <= LINEAR && b.max.x >= 10. - LINEAR);
        assert!(b.min.y <= LINEAR && b.max.y >= 0.);
    }
}

#[test]
fn geom_bnd_lib_offset_circle_along_y_is_not_empty() {
    let c = Curve::Offset(OffsetCurve::new(circle(Frame::WORLD, 3.), 1.5, Direction::Y));
    assert!(!bounds::curve(&c, LINEAR).is_empty());
}

#[test]
fn geom_bnd_lib_offset_plane_patch_is_not_thickened() {
    let s = Surface::Offset(OffsetSurface::new(Surface::Plane(Plane::XY), 3.));
    let b = bounds::surface_within(&s, -4., 6., -2., 7., LINEAR);
    assert_near(b.min.z, 3. - LINEAR, 1e-12);
    assert_near(b.max.z, 3. + LINEAR, 1e-12);
}

#[test]
fn geom_bnd_lib_offset_plane_tight_within_fast() {
    let s = Surface::Offset(OffsetSurface::new(Surface::Plane(Plane::XY), 3.));
    let range = (-5., 5., -5., 5.);
    let fast = bounds::surface_within(&s, -5., 5., -5., 5., LINEAR);
    let tight = bounds::surface_within_tight(&s, -5., 5., -5., 5., LINEAR);
    assert_within(&tight, &fast, LINEAR);
    assert_contains_surface(&tight, &s, range, 10, 10);
}

#[test]
fn geom_bnd_lib_offset_torus_tight_contains_surface() {
    let s = Surface::Offset(OffsetSurface::new(torus(Frame::WORLD, 10., 3.), 1.));
    let fast = bounds::surface_within(&s, 0., 2. * PI, 0., 2. * PI, LINEAR);
    let tight = bounds::surface_within_tight(&s, 0., 2. * PI, 0., 2. * PI, LINEAR);
    assert_within(&tight, &fast, LINEAR);
    assert_contains_surface(&tight, &s, (0., 2. * PI, 0., 2. * PI), 36, 36);
}

#[test]
fn geom_bnd_lib_offset_sphere_tight_contains_surface() {
    let s = Surface::Offset(OffsetSurface::new(sphere(Frame::WORLD, 5.), 2.));
    let fast = bounds::surface(&s, LINEAR);
    let tight = bounds::surface_tight(&s, LINEAR);
    assert_within(&tight, &fast, LINEAR);
    assert_contains_surface(&tight, &s, (0., 2. * PI, -PI / 2., PI / 2.), 36, 18);
}

#[test]
fn geom_bnd_lib_extrusion_infinite_u_is_open_along_basis() {
    let s = Surface::Extrusion(ExtrusionSurface::new(line(p(5., 0., 0.), Direction::Z), Direction::X));
    let b = bounds::surface_within(&s, -INFINITE, INFINITE, 0., 10., LINEAR);
    assert!(!b.is_empty() && b != BoundingBox::WHOLE);
    assert!(b.min.x.is_finite() && b.max.x.is_finite());
    assert!(b.min.y.is_finite() && b.max.y.is_finite());
    assert_eq!(b.min.z, -INF);
    assert_eq!(b.max.z, INF);
    let tol = 2. * LINEAR;
    assert_near(b.min.x, 5., tol);
    assert_near(b.max.x, 15., tol);
    assert_near(b.min.y, 0., tol);
    assert_near(b.max.y, 0., tol);
}

#[test]
fn geom_bnd_lib_extrusion_bspline_along_diagonal_tight() {
    let basis = Curve::BSpline(BSplineCurve::new(
        vec![p(0., 0., 0.), p(3., 5., 0.), p(7., -2., 0.), p(10., 3., 0.)],
        vec![0., 0.5, 1.],
        vec![3, 1, 3],
        2,
    ));
    let s = Surface::Extrusion(ExtrusionSurface::new(basis, Direction::new(1., 1., 1.)));
    let b = bounds::surface_within_tight(&s, 0., 1., -3., 7., LINEAR);
    let (lo, hi) = (-3. / 3f64.sqrt(), 7. / 3f64.sqrt());
    assert_box(&b, [lo, 10. + hi, lo, 3. + hi, lo, hi], 1e-3);
    assert_contains_surface(&b, &s, (0., 1., -3., 7.), 10, 10);
}

#[test]
fn geom_bnd_lib_revolution_infinite_v_is_open_along_axis() {
    let s = Surface::Revolution(RevolutionSurface::new(line(p(5., 0., 0.), Direction::Z), Axis::Z));
    let b = bounds::surface_within(&s, 0., 2. * PI, -INFINITE, INFINITE, LINEAR);
    assert!(!b.is_empty() && b != BoundingBox::WHOLE);
    assert!(b.min.x.is_finite() && b.max.x.is_finite());
    assert!(b.min.y.is_finite() && b.max.y.is_finite());
    assert_eq!(b.min.z, -INF);
    assert_eq!(b.max.z, INF);
    assert_box(&BoundingBox::new(p(b.min.x, b.min.y, 0.), p(b.max.x, b.max.y, 0.)), [-5., 5., -5., 5., 0., 0.], 2. * LINEAR);
}

#[test]
fn geom_bnd_lib_surface_patches() {
    assert_box(&bounds::surface_within(&Surface::Plane(Plane::XY), -5., 5., -3., 3., 0.), [-5., 5., -3., 3., 0., 0.], LINEAR);
    let b = bounds::surface_within(&cylinder(Frame::WORLD, 3.), 0., 2. * PI, 0., 10., LINEAR);
    assert_box(&b, [-3., 3., -3., 3., 0., 10.], LINEAR);
    let b = bounds::surface_within(&cone(PI / 4., 2.), 0., 2. * PI, 0., 5., LINEAR);
    let r = 2. + 5. * (PI / 4.).sin();
    assert!(b.max.x >= r - LINEAR && b.max.y >= r - LINEAR);
    assert!(b.min.x <= -r + LINEAR && b.min.y <= -r + LINEAR);
}

#[test]
fn geom_bnd_lib_sphere_and_torus() {
    assert_box(&bounds::surface(&sphere(Frame::WORLD, 5.), 0.), [-5., 5., -5., 5., -5., 5.], LINEAR);
    let b = bounds::surface_within(&sphere(Frame::WORLD, 5.), 0., PI / 2., 0., PI / 2., LINEAR);
    assert!(b.min.x >= -LINEAR && b.min.y >= -LINEAR && b.min.z >= -LINEAR);
    assert!(b.max.x <= 5. + LINEAR && b.max.y <= 5. + LINEAR && b.max.z <= 5. + LINEAR);
    assert_box(&bounds::surface(&torus(Frame::WORLD, 10., 3.), 0.), [-13., 13., -13., 13., -3., 3.], LINEAR);
}

#[test]
fn geom_bnd_lib_bezier_surface() {
    let s = Surface::Bezier(BezierSurface::new(vec![
        vec![p(0., 0., 0.), p(10., 0., 0.)],
        vec![p(0., 10., 0.), p(10., 10., 5.)],
    ]));
    let b = bounds::surface(&s, LINEAR);
    assert!(b.min.x <= LINEAR && b.max.x >= 10. - LINEAR);
    assert!(b.min.y <= LINEAR && b.max.y >= 10. - LINEAR);
    assert!(b.min.z <= LINEAR && b.max.z >= 5. - LINEAR);
}

#[test]
fn geom_bnd_lib_tight_sphere_within_fast() {
    let s = sphere(Frame::WORLD, 5.);
    assert_within(&bounds::surface_tight(&s, LINEAR), &bounds::surface(&s, LINEAR), LINEAR);
}

#[test]
fn geom_bnd_lib_tight_tilted_torus() {
    let s = torus(Frame::from_z(p(0., 0., 0.), Direction::new(1., 0., 1.)), 8., 3.);
    let fast = bounds::surface(&s, LINEAR);
    let tight = bounds::surface_tight(&s, LINEAR);
    assert_within(&tight, &fast, LINEAR);
    let xz = 8. / 2f64.sqrt() + 3.;
    assert_box(&tight, [-xz, xz, -11., 11., -xz, xz], 10. * LINEAR);
}

#[test]
fn geom_bnd_lib_torus_negative_v_min_no_under_bounding() {
    let s = torus(Frame::WORLD, 10., 3.);
    let range = (0., 2. * PI, -PI / 8., PI / 8.);
    let b = bounds::surface_within(&s, range.0, range.1, range.2, range.3, LINEAR);
    assert_contains_surface(&b, &s, range, 72, 16);
}
