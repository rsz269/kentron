//! Hand translations of OCCT's TKG2d GTests (src/ModelingData/TKG2d/GTests).
//!
//! `Copy*` tests (copy constructors, `Copy()`, copy independence) are not ported: geometry
//! values are plain `Clone` data in kentron.

mod geom2d_api_inter_curve_curve;
mod geom2d_bezier_curve;
mod geom2d_bspline_curve;
mod geom2d_conics;
mod geom2d_curve_eval;
mod geom2d_gcc;
mod geom2d_line;
mod geom2d_offset_curve;
mod geom2d_transformation;
mod geom2d_trimmed_curve;

use kentron::math::{Point2, Vector2};

/// `gp_Pnt2d::IsEqual(other, tol)`.
#[track_caller]
pub fn assert_point2_near(actual: Point2, expected: Point2, tolerance: f64) {
    assert!(
        actual.distance(expected) <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

#[track_caller]
pub fn assert_vector2_near(actual: Vector2, expected: Vector2, tolerance: f64) {
    assert!(
        (actual.x - expected.x).abs() <= tolerance && (actual.y - expected.y).abs() <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

pub fn p2(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}
