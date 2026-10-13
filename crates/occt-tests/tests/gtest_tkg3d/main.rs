//! Hand translations of OCCT's TKG3d GTests (src/ModelingData/TKG3d/GTests).
//!
//! `Copy*` tests are not ported (geometry values are plain `Clone` data); see
//! GTEST_PORTING.md for the other skipped files.

mod geom_api_extrema_curve_curve;
mod geom_api_interpolate;
mod geom_bezier_curve;
mod geom_bezier_surface;
mod geom_bspline_curve;
mod geom_bspline_surface;
mod geom_curve_eval;
mod geom_elementary;
mod geom_offset_surface;
mod geom_surface_eval;

use kentron::math::{Point, Vector};

#[track_caller]
pub fn assert_point_near(actual: Point, expected: Point, tolerance: f64) {
    assert!(
        actual.distance(expected) <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

#[track_caller]
pub fn assert_vector_near(actual: Vector, expected: Vector, tolerance: f64) {
    assert!((actual - expected).length() <= tolerance, "{actual:?} is not within {tolerance} of {expected:?}");
}

pub fn p(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}
