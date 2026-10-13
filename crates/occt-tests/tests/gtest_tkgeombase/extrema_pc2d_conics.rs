//! ExtremaPC2d_Line_Test.cxx, ExtremaPC2d_Circle_Test.cxx, ExtremaPC2d_Ellipse_Test.cxx.
//!
//! Deviation: over several periods kentron reports the extrema sorted by parameter (OCCT
//! lists them period by period, minima first).

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use kentron::geometry::curve2::{Circle2, Curve2, Ellipse2, Line2};
use kentron::geometry::extrema::{ExtremaError, PointCurveSearch, SearchMode};
use kentron::geometry::ParametricCurve;
use kentron::math::{Direction2, Frame2, Point2, Vector2};
use kentron::tolerance::INFINITE;

use crate::{assert_near, p2};

fn line(origin: Point2, direction: Direction2) -> Curve2 {
    Curve2::Line(Line2::new(origin, direction))
}

fn x_line() -> Curve2 {
    line(p2(0., 0.), Direction2::X)
}

fn circle(radius: f64) -> Curve2 {
    Curve2::Circle(Circle2::new(Frame2::WORLD, radius))
}

fn ellipse(major: f64, minor: f64) -> Curve2 {
    Curve2::Ellipse(Ellipse2::new(Frame2::WORLD, major, minor))
}

/// The frame at `origin` with X along +y and Y along +x (clockwise).
fn swapped_frame(origin: Point2) -> Frame2 {
    Frame2::from_axes(origin, Direction2::Y, Direction2::X).unwrap()
}

// ---- lines ----

mod line {
    use super::*;

    const TOL: f64 = 1.0e-10;

    #[test]
    fn unbounded_projection_is_exact() {
        let c = line(p2(1., 2.), Direction2::X);
        let r = PointCurveSearch::new(&c).find(p2(4., 6.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert_near(r[0].parameter, 3., TOL);
        assert_near(r[0].square_distance, 16., TOL);
        assert!(r[0].is_minimum);
    }

    #[test]
    fn reversed_and_diagonal_directions_keep_the_parameter() {
        let c = line(p2(1., 2.), -Direction2::X);
        assert_near(PointCurveSearch::new(&c).find(p2(-2., 5.), TOL).unwrap()[0].parameter, 3., TOL);
        let c = line(p2(0., 0.), Direction2::new(1., 1.));
        let r = PointCurveSearch::new(&c).find(p2(5., 5.), TOL).unwrap();
        assert_near(r[0].parameter, 5. * 2f64.sqrt(), TOL);
        assert_near(r[0].square_distance, 0., TOL);
    }

    #[test]
    fn search_modes_filter_the_projection() {
        let c = x_line();
        let s = PointCurveSearch::new(&c);
        assert_eq!(s.find_with(p2(4., 3.), TOL, SearchMode::Min, false).unwrap().len(), 1);
        assert!(s.find_with(p2(4., 3.), TOL, SearchMode::Max, false).unwrap().is_empty());
    }

    #[test]
    fn bounded_projection_and_ends() {
        let c = x_line();
        let s = PointCurveSearch::with_domain(&c, 0., 10.);
        let r = s.find_with_endpoints(p2(4., 3.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.min().unwrap().parameter, 4., TOL);
        assert_near(r.min_square_distance(), 9., TOL);
        assert_near(r.max().unwrap().parameter, 10., TOL);

        assert!(s.find(p2(-2., 3.), TOL).unwrap().is_empty());
        let r = s.find_with_endpoints(p2(-2., 3.), TOL).unwrap();
        assert_eq!(r.len(), 2);
        assert_near(r.min().unwrap().parameter, 0., TOL);
        assert_near(r.min_square_distance(), 13., TOL);
        assert_near(r.max().unwrap().parameter, 10., TOL);

        let r = s.find_with_endpoints(p2(0., 3.), TOL).unwrap();
        assert_eq!(r.len(), 2, "a projection on an end is not duplicated");
        assert_near(r.min().unwrap().parameter, 0., TOL);
        assert_near(r.max().unwrap().parameter, 10., TOL);
    }

    #[test]
    fn singleton_domain_honors_modes() {
        let c = x_line();
        let s = PointCurveSearch::with_domain(&c, 2., 2.);
        let r = s.find_with_endpoints(p2(5., 3.), TOL).unwrap();
        assert_eq!(r.len(), 1);
        assert!(r[0].is_minimum && r[0].is_maximum);
        let r = s.find_with(p2(5., 3.), TOL, SearchMode::Max, true).unwrap();
        assert_eq!(r.len(), 1);
        assert!(!r[0].is_minimum && r[0].is_maximum);
    }

    #[test]
    fn invalid_input_and_accessors() {
        let c = line(p2(1., 2.), Direction2::Y);
        let s = PointCurveSearch::with_domain(&c, 3., 1.);
        assert_eq!(s.find(p2(0., 0.), TOL), Err(ExtremaError::InvalidInput));
        assert_eq!(s.find(p2(0., 0.), 0.), Err(ExtremaError::InvalidInput));
        assert_eq!(s.domain(), (3., 1.));
        assert_eq!(s.curve().point(2.).distance(p2(1., 4.)), 0.);
    }
}

// ---- circles ----

mod circle {
    use super::*;

    const TOL: f64 = 1.0e-9;

    #[test]
    fn outside_and_inside_points() {
        let c = circle(5.);
        let s = PointCurveSearch::new(&c);
        let r = s.find(p2(13., 0.), TOL).unwrap();
        assert_eq!(r.len(), 2);
        assert_near(r.min().unwrap().parameter, 0., TOL);
        assert_near(r.min_square_distance(), 64., TOL);
        assert_near(r.max_square_distance(), 324., TOL);
        let r = s.find(p2(2., 0.), TOL).unwrap();
        assert_near(r.min_square_distance(), 9., TOL);
        assert_near(r.max_square_distance(), 49., TOL);
    }

    #[test]
    fn search_modes_return_subsets() {
        let c = circle(5.);
        let s = PointCurveSearch::new(&c);
        let r = s.find_with(p2(13., 0.), TOL, SearchMode::Min, false).unwrap();
        assert_eq!(r.len(), 1);
        assert!(r[0].is_minimum);
        let r = s.find_with(p2(13., 0.), TOL, SearchMode::Max, false).unwrap();
        assert_eq!(r.len(), 1);
        assert!(r[0].is_maximum);
        assert_near(r[0].parameter, PI, TOL);
    }

    #[test]
    fn clockwise_frame_keeps_its_sense() {
        let frame = swapped_frame(p2(2., -1.));
        let c = Curve2::Circle(Circle2::new(frame, 4.));
        let angle: f64 = 0.7;
        let offset = Vector2::from(frame.x()) * angle.cos() + Vector2::from(frame.y()) * angle.sin();
        let q = frame.origin() + offset * 7.;
        let r = PointCurveSearch::new(&c).find(q, TOL).unwrap();
        assert_near(r.min().unwrap().parameter, angle, TOL);
        assert_near(r.min_square_distance(), 9., TOL);
    }

    #[test]
    fn center_and_zero_radius_are_infinite() {
        let r = PointCurveSearch::new(&circle(5.)).find(p2(0., 0.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
        let r = PointCurveSearch::new(&circle(0.)).find(p2(3., 4.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
    }

    #[test]
    fn shifted_full_period_has_unique_representatives() {
        let c = circle(5.);
        let r = PointCurveSearch::with_domain(&c, 7., 7. + 2. * PI).find_with_endpoints(p2(13., 0.), TOL).unwrap();
        assert_eq!(r.len(), 2);
        assert_near(r.min().unwrap().parameter, 4. * PI, TOL);
        assert_near(r.max().unwrap().parameter, 3. * PI, TOL);
    }

    #[test]
    fn multiple_periods_return_every_extremum() {
        let c = circle(5.);
        let r = PointCurveSearch::with_domain(&c, 0., 4. * PI).find(p2(13., 0.), TOL).unwrap();
        assert_eq!(r.len(), 4);
        for (e, t) in r.iter().zip([0., PI, 2. * PI, 3. * PI]) {
            assert_near(e.parameter, t, TOL);
        }
    }

    #[test]
    fn half_circle_ends() {
        let c = circle(5.);
        let r = PointCurveSearch::with_domain(&c, FRAC_PI_2, 1.5 * PI).find_with_endpoints(p2(10., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.max().unwrap().parameter, PI, TOL);
        assert_eq!(r.iter().filter(|e| e.is_minimum).count(), 2);
    }

    #[test]
    fn singleton_and_invalid_domains() {
        let c = circle(5.);
        let s = PointCurveSearch::with_domain(&c, 0.4, 0.4);
        assert!(s.find(p2(8., 0.), TOL).unwrap().is_empty());
        assert!(s.find(p2(0., 0.), TOL).unwrap().is_empty());
        assert_eq!(s.find_with_endpoints(p2(8., 0.), TOL).unwrap().len(), 1);
        let s = PointCurveSearch::with_domain(&c, 1., 0.);
        assert_eq!(s.find(p2(0., 0.), TOL), Err(ExtremaError::InvalidInput));
        assert_eq!(s.find(p2(0., 0.), 0.), Err(ExtremaError::InvalidInput));
    }

    #[test]
    fn large_tolerance_does_not_expand_the_domain() {
        let c = circle(1.);
        assert!(PointCurveSearch::with_domain(&c, 1., 2.).find(p2(2., 0.), 2.).unwrap().is_empty());
    }

    #[test]
    fn non_finite_query_is_rejected() {
        let c = circle(5.);
        assert_eq!(PointCurveSearch::new(&c).find(p2(INFINITE, 0.), TOL), Err(ExtremaError::InvalidInput));
    }

    #[test]
    fn unrepresentable_periodic_step_is_a_numerical_error() {
        let c = circle(5.);
        let first = 1.0e17;
        let s = PointCurveSearch::with_domain(&c, first, first + 100.);
        assert_eq!(s.find(p2(12., 0.), TOL), Err(ExtremaError::NumericalError));
    }

    #[test]
    fn accessors() {
        let frame = Frame2::new(p2(2., 3.), Direction2::Y, false);
        let c = Curve2::Circle(Circle2::new(frame, 4.));
        let s = PointCurveSearch::with_domain(&c, -1., 2.);
        assert_eq!(s.domain(), (-1., 2.));
        assert_eq!(s.curve().point(0.3), c.point(0.3));
    }
}

// ---- ellipses ----

mod ellipse {
    use super::*;

    const TOL: f64 = 1.0e-8;

    #[test]
    fn major_axis_outside_has_two_extrema() {
        let r = PointCurveSearch::new(&ellipse(5., 3.)).find(p2(10., 0.), TOL).unwrap();
        assert_eq!(r.len(), 2);
        assert_near(r.min().unwrap().parameter, 0., TOL);
        assert_near(r.min_square_distance(), 25., TOL);
        assert_near(r.max_square_distance(), 225., TOL);
    }

    #[test]
    fn center_has_four_axis_extrema() {
        let c = ellipse(5., 3.);
        let s = PointCurveSearch::new(&c);
        let r = s.find(p2(0., 0.), TOL).unwrap();
        assert_eq!(r.len(), 4);
        assert_eq!(r.iter().filter(|e| e.is_minimum).count(), 2);
        assert_eq!(r.iter().filter(|e| e.is_maximum).count(), 2);
        assert_near(r.min_square_distance(), 9., TOL);
        assert_near(r.max_square_distance(), 25., TOL);
        assert_eq!(s.find_with(p2(0., 0.), TOL, SearchMode::Min, false).unwrap().len(), 2);
        assert_eq!(s.find_with(p2(0., 0.), TOL, SearchMode::Max, false).unwrap().len(), 2);
    }

    #[test]
    fn point_on_curve() {
        let c = ellipse(5., 3.);
        let r = PointCurveSearch::new(&c).find(c.point(FRAC_PI_4), TOL).unwrap();
        assert!(r.len() >= 2);
        assert_near(r.min().unwrap().parameter, FRAC_PI_4, TOL);
        assert_near(r.min_square_distance(), 0., TOL);
    }

    #[test]
    fn clockwise_frame_maps_local_axes() {
        let c = Curve2::Ellipse(Ellipse2::new(swapped_frame(p2(1., 2.)), 5., 3.));
        let r = PointCurveSearch::new(&c).find(p2(11., 2.), TOL).unwrap();
        assert_near(r.min().unwrap().parameter, FRAC_PI_2, TOL);
        assert_near(r.min_square_distance(), 49., TOL);
    }

    #[test]
    fn circular_ellipse_center_is_infinite() {
        let r = PointCurveSearch::new(&ellipse(5., 5.)).find(p2(0., 0.), TOL).unwrap();
        assert!(r.is_infinite());
        assert_near(r.infinite_square_distance.unwrap(), 25., TOL);
    }

    #[test]
    fn shifted_period_and_half_arc() {
        let c = ellipse(5., 3.);
        let r = PointCurveSearch::with_domain(&c, 7., 7. + 2. * PI).find_with_endpoints(p2(10., 0.), TOL).unwrap();
        assert_eq!(r.len(), 2);
        let r = PointCurveSearch::with_domain(&c, FRAC_PI_2, 1.5 * PI).find_with_endpoints(p2(10., 0.), TOL).unwrap();
        assert_eq!(r.len(), 3);
        assert_near(r.max().unwrap().parameter, PI, TOL);
    }

    #[test]
    fn multiple_periods_return_every_extremum() {
        let c = ellipse(5., 3.);
        let r = PointCurveSearch::with_domain(&c, 0., 4. * PI).find(p2(12., 0.), TOL).unwrap();
        assert_eq!(r.len(), 4);
        assert_eq!(r.iter().filter(|e| e.parameter < 2. * PI).count(), 2);
    }

    #[test]
    fn singleton_and_invalid_domains() {
        let c = ellipse(5., 3.);
        let s = PointCurveSearch::with_domain(&c, 0.4, 0.4);
        assert!(s.find(p2(0., 0.), TOL).unwrap().is_empty());
        assert_eq!(s.find_with_endpoints(p2(0., 0.), TOL).unwrap().len(), 1);
        let s = PointCurveSearch::with_domain(&c, 1., 0.);
        assert_eq!(s.find(p2(0., 0.), TOL), Err(ExtremaError::InvalidInput));
    }
}
