//! Geom2dGcc_Circ2d2TanOn_Test.cxx and Geom2dGcc_Circ2d2TanRad_Test.cxx: circles tangent
//! to curves.

use kentron::geometry::curve2::{BSplineCurve2, BezierCurve2, Circle2, Curve2, Ellipse2, Line2};
use kentron::geometry::tangency::{self, Qualifier, Tangent};
use kentron::geometry::{extrema, ParametricCurve};
use kentron::math::{Direction2, Frame2};
use kentron::tolerance;

use crate::p2;

fn circle_length(c: &Circle2) -> f64 {
    kentron::geometry::measure::length2(&Curve2::Circle(*c))
}

#[test]
fn tangent_and_centered_on_normal_does_not_panic() {
    let c1 = Curve2::Bezier(BezierCurve2::new(vec![p2(0.0, 0.0), p2(0.0, 1.0), p2(6.0, 0.0)]));
    let c2 = Curve2::Bezier(BezierCurve2::new(vec![p2(0.0, 0.0), p2(0.0, 1.5), p2(6.0, 0.0)]));
    for i in 0..100 {
        let u = i as f64 / 99.0;
        let (p, tangent) = c1.d1(u);
        let normal = Line2::new(p, Direction2::new(-tangent.y, tangent.x));
        let normal = Curve2::Line(normal);
        let g1 = extrema::point_curve2(p, &c1)[0].parameter;
        let g3 = extrema::point_curve2(p, &normal)[0].parameter;
        let _ = tangency::centered_on(
            Tangent::Curve(&c1, Qualifier::Unqualified),
            Tangent::Curve(&c2, Qualifier::Unqualified),
            &normal,
            1.0e-9,
            (g1, g1, g3),
        );
    }
}

#[test]
fn tangent_to_line_and_bezier_with_radius() {
    let line = Curve2::Line(Line2::new(p2(100.0, 0.0), Direction2::new(-1.0, 0.0)));
    let bezier = Curve2::Bezier(BezierCurve2::new(vec![p2(0.0, 0.0), p2(50.0, 50.0), p2(0.0, 100.0)]));
    let solutions = tangency::with_radius(
        Tangent::Curve(&line, Qualifier::Outside),
        Tangent::Curve(&bezier, Qualifier::Outside),
        10.0,
        1.0e-7,
    )
    .unwrap();
    assert!(!solutions.is_empty());
    for s in &solutions {
        let (center, r) = (s.circle.center(), s.circle.radius());
        for t in &s.tangencies {
            assert!((t.point.distance(center) - r).abs() / r * 100.0 <= 1.0);
        }
    }
}

fn slanted_line(x: f64, y: f64, dx: f64, dy: f64) -> Curve2 {
    Curve2::Line(Line2::new(p2(x, y), Direction2::new(dx, dy)))
}

#[test]
fn point_and_line_with_radius() {
    let line = slanted_line(0.0, 500.0, 2000.0, 500.0);
    let solve = |radius| {
        tangency::with_radius(Tangent::Curve(&line, Qualifier::Unqualified), Tangent::Point(p2(1600.0, 300.0)), radius, tolerance::LINEAR)
            .unwrap()
    };
    let solutions = solve(400.0);
    assert_eq!(solutions.len(), 2);
    for s in &solutions {
        let expected = 2513.2741228718346;
        assert!((circle_length(&s.circle) - expected).abs() <= expected * 1.0e-8);
    }
    assert!(solve(200.0).is_empty());
}

#[test]
fn point_and_line_with_radius_no_solution() {
    let line = slanted_line(3000.0, 1700.0, 1000.0, 1550.0);
    for radius in [300.0, 150.0] {
        let solutions = tangency::with_radius(
            Tangent::Curve(&line, Qualifier::Unqualified),
            Tangent::Point(p2(1400.0, 1400.0)),
            radius,
            tolerance::LINEAR,
        )
        .unwrap();
        assert!(solutions.is_empty());
    }
}

#[test]
fn circle_and_bspline_with_radius_does_not_fail() {
    let circle = Curve2::Circle(Circle2::new(
        Frame2::new(p2(1.79356127430613, -11.9148900515605), Direction2::X, true),
        12.0491272292007,
    ));
    let poles = vec![
        p2(0.0, 0.0),
        p2(-2.14533906935467, 4.92056269221697),
        p2(-0.537398638327716, 10.0937503357781),
        p2(10.3564232507092, 8.73929826869582),
        p2(16.6366922939454, 3.11858911810943),
        p2(16.39, -0.262),
    ];
    let bspline = Curve2::BSpline(
        BSplineCurve2::new(poles, vec![0.0, 8.6971225126475, 16.903955884905, 28.0087167558108], vec![4, 1, 1, 4], 3)
            .with_weights(vec![1.0; 6]),
    );
    let result = tangency::with_radius(
        Tangent::Curve(&circle, Qualifier::Unqualified),
        Tangent::Curve(&bspline, Qualifier::Unqualified),
        2.0,
        tolerance::LINEAR,
    );
    assert!(result.is_ok());
}

#[test]
fn tangent_to_segment_with_center() {
    let segment = Curve2::BSpline(BSplineCurve2::new(vec![p2(-10.0, -5.0), p2(10.0, -5.0)], vec![0.0, 1.0], vec![2, 2], 1));
    let solutions = tangency::with_center(Tangent::Curve(&segment, Qualifier::Unqualified), p2(1.0, 1.0), tolerance::LINEAR).unwrap();
    assert_eq!(solutions.len(), 1);
    let c = &solutions[0].circle;
    assert!(c.center().distance(p2(1.0, 1.0)) <= 1.0e-7);
    assert!((c.radius() - 6.0).abs() <= 1.0e-7);
}

#[test]
fn tangent_to_two_ellipses_with_radius() {
    let ellipse = |x: f64| Curve2::Ellipse(Ellipse2::new(Frame2::new(p2(x, 0.0), Direction2::X, true), 2.0, 1.0));
    let (a, b) = (ellipse(0.0), ellipse(4.0));
    let solutions = tangency::with_radius(
        Tangent::Curve(&a, Qualifier::Unqualified),
        Tangent::Curve(&b, Qualifier::Unqualified),
        3.0,
        1.0e-9,
    )
    .unwrap();
    assert!(!solutions.is_empty());
    for s in &solutions {
        assert!((s.circle.radius() - 3.0).abs() <= 1.0e-6);
    }
    assert!(solutions[0].circle.center().distance(p2(5.0, 0.0)) < 10.0);
}
