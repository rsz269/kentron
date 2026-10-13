//! Geom2dAPI_InterCurveCurve_Test.cxx: 2D curve intersection (and the two extrema cases
//! filed under it).
//!
//! `OCC29289_EllipseIntersectionNewtonRoot` checks the roots with OCCT's Newton solver;
//! here each intersection point is checked on both ellipses instead. Not ported:
//! `PointRejectsZeroIndex`, `Init_NullHandle_RaisesException` (no 1-based index or null
//! curves).

use kentron::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Line2, TrimmedCurve2};
use kentron::geometry::{extrema, intersect, ParametricCurve};
use kentron::math::{Direction2, Frame2};
use kentron::tolerance;

use crate::{assert_point2_near, p2};

fn ellipse(x: f64, y: f64, dir: Direction2, major: f64, minor: f64) -> Curve2 {
    Curve2::Ellipse(Ellipse2::new(Frame2::new(p2(x, y), dir, true), major, minor))
}

fn occ23706_curve1() -> Curve2 {
    let poles = vec![p2(2.0, 0.0), p2(3.0, -1.0), p2(5.0, 5.0), p2(5.0, 5.0), p2(6.0, 8.0), p2(4.0, 7.0)];
    Curve2::BSpline(
        BSplineCurve2::new(poles, vec![0.0, 0.2, 0.3, 0.4, 0.5], vec![3, 1, 1, 1, 3], 2).with_weights(vec![1.0; 6]),
    )
}

fn occ23706_curve2(perturbed: bool) -> Curve2 {
    let p = if perturbed { p2(5.001, 5.01) } else { p2(5.0, 5.0) };
    let poles = vec![p2(6.0, 3.0), p, p, p2(3.0, 9.0), p2(2.0, 11.0)];
    Curve2::BSpline(BSplineCurve2::new(poles, vec![1.0, 2.0, 3.0, 5.0], vec![3, 1, 1, 3], 2).with_weights(vec![1.0; 5]))
}

#[test]
fn intersecting_ellipses_points_lie_on_both() {
    let a = ellipse(0.0, 0.0, Direction2::X, 2.0, 1.0);
    let b = ellipse(0.5, 0.5, Direction2::new(1.0, 1.0), 2.0, 1.0);
    let points = intersect::curves2(&a, &b, 1.0e-7);
    assert!(!points.is_empty());
    for x in points {
        assert_point2_near(a.point(x.t1), x.point, 1.0e-7);
        assert_point2_near(b.point(x.t2), x.point, 1.0e-7);
    }
}

#[test]
fn bspline_intersection_perturbed() {
    let points = intersect::curves2(&occ23706_curve1(), &occ23706_curve2(true), 0.001);
    assert_eq!(points.len(), 2);
    assert!((points[0].point.x - 5.0024).abs() <= 1.0e-4 && (points[0].point.y - 5.0072).abs() <= 1.0e-4);
    assert!((points[1].point.x - 4.0024).abs() <= 1.0e-4 && (points[1].point.y - 7.0012).abs() <= 1.0e-4);
}

#[test]
fn bspline_intersection_exact() {
    let points = intersect::curves2(&occ23706_curve1(), &occ23706_curve2(false), 0.001);
    assert_eq!(points.len(), 2);
    assert_point2_near(points[0].point, p2(5.0, 5.0), 1.0e-7);
    let q = points[1].point;
    assert!((4.0..5.0).contains(&q.x) && (7.0..8.0).contains(&q.y));
}

#[test]
fn trimmed_ellipse_and_line() {
    let e = Ellipse2::new(
        Frame2::new(p2(-610.348096534595, -710.720096056787), Direction2::new(0.999999902285153, 0.000442074298181498), true),
        15.0066332711999,
        0.291884102212871,
    );
    let arc = Curve2::Trimmed(TrimmedCurve2::new(Curve2::Ellipse(e), 3.09462291909258, 9.37780822627216));
    let line = Line2::new(p2(-625.34430362036, -680.713463264921), Direction2::new(-0.000252749178714602, -0.999999968058926));
    let segment = Curve2::Trimmed(TrimmedCurve2::new(Curve2::Line(line), 0.0, 30.0132665523925));
    let points = intersect::curves2(&arc, &segment, 0.001);
    assert_eq!(points.len(), 1);
    assert!((points[0].point.x - -625.35188801291508).abs() <= 1.0e-7);
    assert!((points[0].point.y - -710.72104765746838).abs() <= 1.0e-7);
}

#[test]
fn intersection_parameters_within_trimmed_ranges() {
    let circle = |y: f64| Curve2::Circle(Circle2::new(Frame2::new(p2(25.0, y), Direction2::X, true), 155.0));
    let a = TrimmedCurve2::new(circle(-25.0), 1.57079632679490, 2.97959469729228);
    let b = TrimmedCurve2::new(circle(25.0), 3.30359060633978, 4.71238898038469);
    let points = intersect::curves2(&Curve2::Trimmed(a.clone()), &Curve2::Trimmed(b.clone()), tolerance::LINEAR);
    let x = points.first().expect("no intersection");
    assert!((a.first_parameter()..=a.last_parameter()).contains(&x.t1));
    assert!((b.first_parameter()..=b.last_parameter()).contains(&x.t2));
    assert!(a.point(x.t1).distance(b.point(x.t2)) < 1.0e-7);
}

#[test]
fn tangent_ellipses_touch() {
    let a = ellipse(0.0, 0.0, Direction2::X, 2.0, 1.0);
    let b = ellipse(4.0, 0.0, Direction2::X, 2.0, 1.0);
    let result = extrema::curve_curve2(&a, &b);
    assert!(result.min_distance().expect("no extrema") <= 7.0e-5);
}

#[test]
fn tangent_circles_touch() {
    let circle = |x: f64| Curve2::Circle(Circle2::new(Frame2::new(p2(x, 0.0), Direction2::X, true), 2.0));
    let result = extrema::curve_curve2(&circle(0.0), &circle(4.0));
    let isolated = result.isolated();
    assert!(isolated.len() >= 4);
    // Sorted by distance: the touching pair comes first.
    assert!(isolated[0].distance <= tolerance::PARAMETRIC);
    assert_point2_near(isolated[0].p1, isolated[0].p2, tolerance::PARAMETRIC);
}
