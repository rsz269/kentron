//! ExtremaPC2d_Hyperbola_Test.cxx, ExtremaPC2d_Parabola_Test.cxx,
//! ExtremaPC2d_BezierCurve_Test.cxx, ExtremaPC2d_BSplineCurve_Test.cxx.
//!
//! The null-handle and handle-ownership checks have no counterpart (a search borrows its
//! curve); "mutation is observed" becomes a search on the edited curve.

use kentron::geometry::curve2::{BSplineCurve2, BezierCurve2, Curve2, Hyperbola2, Parabola2};
use kentron::geometry::extrema::{ExtremaError, PointCurveSearch, SearchMode};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction2, Frame2, Point2};
use kentron::tolerance::PARAMETRIC;

use crate::{assert_near, p2};

const TOL: f64 = 1.0e-8;

/// The frame at `origin` with X along +y and Y along +x (clockwise).
fn swapped_frame(origin: Point2) -> Frame2 {
    Frame2::from_axes(origin, Direction2::Y, Direction2::X).unwrap()
}

/// Checks a singleton domain and a reversed one on `c`.
#[track_caller]
fn check_singleton_and_invalid(c: &Curve2) {
    let s = PointCurveSearch::with_domain(c, 0.4, 0.4);
    assert!(s.find(p2(0., 0.), TOL).unwrap().is_empty());
    assert_eq!(s.find_with_endpoints(p2(0., 0.), TOL).unwrap().len(), 1);
    let s = PointCurveSearch::with_domain(c, 1., 0.);
    assert_eq!(s.find(p2(0., 0.), TOL), Err(ExtremaError::InvalidInput));
}

mod hyperbola {
    use super::*;

    fn hyperbola(major: f64, minor: f64) -> Curve2 {
        Curve2::Hyperbola(Hyperbola2::new(Frame2::WORLD, major, minor))
    }

    #[test]
    fn center_has_vertex_minimum() {
        let r = PointCurveSearch::new(&hyperbola(2., 1.)).find(p2(0., 0.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
        assert_near(r[0].square_distance, 4., TOL);
        assert!(r[0].is_minimum);
    }

    #[test]
    fn major_axis_point_has_three_extrema() {
        let root = 2f64.acosh();
        let r = PointCurveSearch::new(&hyperbola(2., 1.)).find(p2(5., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_eq!(r.iter().filter(|e| e.is_minimum).count(), 2);
        assert_eq!(r.iter().filter(|e| e.is_maximum).count(), 1);
        assert!(r.iter().any(|e| (e.parameter + root).abs() <= TOL));
        assert!(r.iter().any(|e| (e.parameter - root).abs() <= TOL));
        assert_near(r.min_square_distance(), 4., TOL);
        assert_near(r.max_square_distance(), 9., TOL);
    }

    #[test]
    fn root_near_the_domain_end_is_clamped() {
        let c = hyperbola(2., 1.);
        let upper = 2f64.acosh() - 0.5 * PARAMETRIC;
        let r = PointCurveSearch::with_domain(&c, 0.5, upper).find(p2(5., 0.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].parameter, upper);
    }

    #[test]
    fn search_modes() {
        let c = hyperbola(2., 1.);
        let s = PointCurveSearch::new(&c);
        let r = s.find_with(p2(5., 0.), TOL, SearchMode::Min, false).unwrap();
        assert_eq!(r.len(), 2);
        assert!(r.iter().all(|e| e.is_minimum));
        let r = s.find_with(p2(5., 0.), TOL, SearchMode::Max, false).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
        assert!(r[0].is_maximum);
    }

    #[test]
    fn point_on_curve() {
        let c = hyperbola(2., 1.);
        let r = PointCurveSearch::new(&c).find(c.point(0.7), TOL).unwrap();
        assert_near(r.min_square_distance(), 0., TOL);
        assert_near(r.min().unwrap().parameter, 0.7, TOL);
    }

    #[test]
    fn clockwise_frame() {
        let c = Curve2::Hyperbola(Hyperbola2::new(swapped_frame(p2(1., 2.)), 2., 1.));
        let r = PointCurveSearch::new(&c).find(p2(1., 7.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.min_square_distance(), 4., TOL);
        assert_near(r.max_square_distance(), 9., TOL);
    }

    #[test]
    fn bounded_arc_adds_end_minima() {
        let c = hyperbola(2., 1.);
        let r = PointCurveSearch::with_domain(&c, -0.5, 0.5).find_with_endpoints(p2(5., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.max().unwrap().parameter, 0., TOL);
        assert_eq!(r.iter().filter(|e| e.is_minimum).count(), 2);
    }

    #[test]
    fn zero_axes() {
        let r = PointCurveSearch::new(&hyperbola(0., 1.)).find(p2(3., 2.), TOL).unwrap();
        assert!(!r.is_empty());
        assert_near(r.min_square_distance(), 9., TOL);
        let r = PointCurveSearch::new(&hyperbola(0., 0.)).find(p2(3., 4.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
    }

    #[test]
    fn singleton_and_invalid_domains() {
        check_singleton_and_invalid(&hyperbola(2., 1.));
    }
}

mod parabola {
    use super::*;

    fn parabola(focal: f64) -> Curve2 {
        Curve2::Parabola(Parabola2::new(Frame2::WORLD, focal))
    }

    #[test]
    fn vertex_and_behind_vertex_have_a_unique_minimum() {
        let c = parabola(1.);
        let s = PointCurveSearch::new(&c);
        let r = s.find(p2(0., 0.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
        assert_near(r[0].square_distance, 0., TOL);
        let r = s.find(p2(-2., 0.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].square_distance, 4., TOL);
    }

    #[test]
    fn axis_point_ahead_has_three_extrema() {
        let c = parabola(1.);
        let s = PointCurveSearch::new(&c);
        let r = s.find(p2(3., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.min_square_distance(), 8., TOL);
        assert_near(r.max_square_distance(), 9., TOL);
        assert_eq!(r.iter().filter(|e| e.is_minimum).count(), 2);
        assert_eq!(s.find_with(p2(3., 0.), TOL, SearchMode::Min, false).unwrap().len(), 2);
        let r = s.find_with(p2(3., 0.), TOL, SearchMode::Max, false).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
    }

    #[test]
    fn known_cubic_has_three_roots() {
        let r = PointCurveSearch::new(&parabola(1.)).find(p2(3.75, -0.75), TOL).unwrap();
        assert_eq!(r.len(), 3);
        for t in [-3., 1., 2.] {
            assert!(r.iter().any(|e| (e.parameter - t).abs() <= TOL), "root {t}");
        }
    }

    #[test]
    fn root_near_the_domain_end_is_clamped() {
        let c = parabola(1.);
        let upper = 1. - 0.5 * PARAMETRIC;
        let r = PointCurveSearch::with_domain(&c, 0.5, upper).find(p2(3.75, -0.75), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].parameter, upper);
    }

    #[test]
    fn point_on_curve() {
        let c = parabola(1.);
        let r = PointCurveSearch::new(&c).find(c.point(2.), TOL).unwrap();
        assert_near(r.min_square_distance(), 0., TOL);
        assert_near(r.min().unwrap().parameter, 2., TOL);
    }

    #[test]
    fn clockwise_frame() {
        let c = Curve2::Parabola(Parabola2::new(swapped_frame(p2(1., 2.)), 1.));
        let r = PointCurveSearch::new(&c).find(p2(1., 5.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.min_square_distance(), 8., TOL);
    }

    #[test]
    fn bounded_arc_adds_end_minima() {
        let c = parabola(1.);
        let r = PointCurveSearch::with_domain(&c, -1., 1.).find_with_endpoints(p2(3., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.max().unwrap().parameter, 0., TOL);
    }

    #[test]
    fn zero_focal_is_its_axis_line() {
        let c = Curve2::Parabola(Parabola2::new(Frame2::new(p2(1., 2.), Direction2::Y, true), 0.));
        let r = PointCurveSearch::new(&c).find(p2(4., 7.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 5., TOL);
        assert_near(r[0].square_distance, 9., TOL);
    }

    #[test]
    fn singleton_and_invalid_domains() {
        check_singleton_and_invalid(&parabola(1.));
    }
}

mod bezier {
    use super::*;

    fn linear() -> Curve2 {
        Curve2::Bezier(BezierCurve2::new(vec![p2(0., 0.), p2(10., 0.)]))
    }

    fn quadratic() -> BezierCurve2 {
        BezierCurve2::new(vec![p2(-2., 0.), p2(0., 2.), p2(2., 0.)])
    }

    #[test]
    fn linear_projection_is_exact() {
        let r = PointCurveSearch::new(&linear()).find(p2(4., 3.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0.4, TOL);
        assert_near(r[0].square_distance, 9., TOL);
    }

    #[test]
    fn search_modes_and_ends() {
        let c = linear();
        let s = PointCurveSearch::new(&c);
        assert!(s.find_with(p2(4., 3.), TOL, SearchMode::Max, false).unwrap().is_empty());
        assert!(s.find(p2(-2., 0.), TOL).unwrap().is_empty());
        let r = s.find_with(p2(-2., 0.), TOL, SearchMode::Min, true).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0., TOL);
        let r = s.find_with(p2(-2., 0.), TOL, SearchMode::Max, true).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 1., TOL);
        assert_near(r[0].square_distance, 144., TOL);
    }

    #[test]
    fn symmetric_quadratic_has_a_central_minimum() {
        let c = Curve2::Bezier(quadratic());
        let r = PointCurveSearch::new(&c).find_with(p2(0., 0.), TOL, SearchMode::Min, false).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0.5, TOL);
        assert_near(r[0].square_distance, 1., TOL);
    }

    #[test]
    fn rational_quarter_circle_is_infinite_from_its_center() {
        let c = Curve2::Bezier(BezierCurve2::rational(
            vec![p2(1., 0.), p2(1., 1.), p2(0., 1.)],
            vec![1., std::f64::consts::FRAC_1_SQRT_2, 1.],
        ));
        let r = PointCurveSearch::new(&c).find(p2(0., 0.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 1., TOL);
    }

    #[test]
    fn constant_curves_are_infinite() {
        let c = Curve2::Bezier(BezierCurve2::new(vec![p2(2., -1.); 3]));
        let r = PointCurveSearch::new(&c).find(p2(5., 3.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);

        let n = BezierCurve2::MAX_DEGREE + 1;
        let weights = (1..=n).map(|i| 1. + 0.01 * i as f64).collect();
        let c = Curve2::Bezier(BezierCurve2::rational(vec![p2(2., -1.); n], weights));
        let r = PointCurveSearch::new(&c).find(p2(5., 3.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
    }

    #[test]
    fn edited_pole_changes_the_result() {
        let mut bezier = quadratic();
        let q = p2(0., 1.5);
        let before = PointCurveSearch::new(&Curve2::Bezier(bezier.clone())).find(q, TOL).unwrap().min_square_distance();
        bezier.set_pole(1, p2(3., -4.));
        let c = Curve2::Bezier(bezier);
        let r = PointCurveSearch::new(&c).find(q, TOL).unwrap();
        assert!((r.min_square_distance() - before).abs() > 1.0e-4);
        let min = r.min().unwrap();
        assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    }

    #[test]
    fn singleton_and_invalid_domains() {
        check_singleton_and_invalid(&linear());
    }
}

mod bspline {
    use super::*;

    fn linear() -> BSplineCurve2 {
        BSplineCurve2::new(vec![p2(0., 0.), p2(10., 0.)], vec![0., 1.], vec![2, 2], 1)
    }

    fn c1_curve() -> BSplineCurve2 {
        BSplineCurve2::new(
            vec![p2(0., 0.), p2(1., 1.), p2(2., 1.), p2(3., 0.)],
            vec![0., 0.5, 1.],
            vec![3, 1, 3],
            2,
        )
    }

    #[test]
    fn linear_projection_is_exact() {
        let r = PointCurveSearch::new(&Curve2::BSpline(linear())).find(p2(6., 4.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0.6, TOL);
        assert_near(r[0].square_distance, 16., TOL);
    }

    #[test]
    fn partial_domain_uses_the_end_when_asked() {
        let c = Curve2::BSpline(linear());
        let s = PointCurveSearch::with_domain(&c, 0., 0.5);
        assert!(s.find(p2(8., 3.), TOL).unwrap().is_empty());
        let r = s.find_with(p2(8., 3.), TOL, SearchMode::Min, true).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 0.5, TOL);
    }

    #[test]
    fn c1_junction_point_is_found() {
        let c = Curve2::BSpline(c1_curve());
        let r = PointCurveSearch::new(&c).find_with(c.point(0.5), TOL, SearchMode::Min, false).unwrap();
        assert!(!r.is_empty());
        assert_near(r.min_square_distance(), 0., TOL);
        assert_near(r.min().unwrap().parameter, 0.5, TOL);
    }

    #[test]
    fn constant_curve_is_infinite() {
        let mut bspline = linear();
        bspline.set_pole(0, p2(2., -1.));
        bspline.set_pole(1, p2(2., -1.));
        let r = PointCurveSearch::new(&Curve2::BSpline(bspline)).find(p2(5., 3.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
    }

    #[test]
    fn edited_pole_changes_the_result() {
        let mut bspline = c1_curve();
        let q = p2(1.5, 2.);
        let before = PointCurveSearch::new(&Curve2::BSpline(bspline.clone())).find(q, TOL).unwrap().min_square_distance();
        bspline.set_pole(1, p2(1., -3.));
        let c = Curve2::BSpline(bspline);
        let r = PointCurveSearch::new(&c).find(q, TOL).unwrap();
        assert!((r.min_square_distance() - before).abs() > 1.0e-4);
        let min = r.min().unwrap();
        assert_near(min.point.distance(c.point(min.parameter)), 0., TOL);
    }

    #[test]
    fn search_modes_filter_extrema() {
        let c = Curve2::BSpline(c1_curve());
        let s = PointCurveSearch::new(&c);
        assert!(s.find_with(p2(1.5, 2.), TOL, SearchMode::Min, false).unwrap().iter().all(|e| e.is_minimum));
        assert!(s.find_with(p2(1.5, 2.), TOL, SearchMode::Max, false).unwrap().iter().all(|e| e.is_maximum));
    }

    #[test]
    fn singleton_and_invalid_domains() {
        check_singleton_and_invalid(&Curve2::BSpline(linear()));
    }
}
