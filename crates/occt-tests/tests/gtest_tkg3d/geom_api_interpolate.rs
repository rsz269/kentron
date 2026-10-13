//! GeomAPI_Interpolate_Test.cxx.

use kentron::geometry::{fit, ParametricCurve};
use kentron::math::Point;
use kentron::tolerance;

#[test]
fn end_tangents_are_preserved() {
    let points: Vec<Point> = (0..5).map(|i| Point::new(i as f64 * 1.57, (i as f64 * 1.57).sin(), 0.0)).collect();
    let curve = fit::interpolate(&points, false, tolerance::LINEAR).unwrap();
    let first = curve.d1(curve.first_parameter()).1;
    let last = curve.d1(curve.last_parameter()).1;
    let again = fit::interpolate_with_end_tangents(&points, first, last, false, tolerance::LINEAR).unwrap();
    let (first2, last2) = (again.d1(again.first_parameter()).1, again.d1(again.last_parameter()).1);
    assert!(first.approx_eq(first2, tolerance::LINEAR) && first.is_parallel(first2, tolerance::ANGULAR));
    assert!(last.approx_eq(last2, tolerance::LINEAR) && last.is_parallel(last2, tolerance::ANGULAR));
}
