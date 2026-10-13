//! Hand translations of OCCT's TKGeomBase GTests (src/ModelingData/TKGeomBase/GTests).
//!
//! See GTEST_PORTING.md for the files that are not ported.

mod bndlib;
mod construct;
mod extrema_pc_bezier;
mod extrema_pc_bspline;
mod extrema_pc_circle;
mod extrema_pc_comparison;
mod extrema_pc_curve;
mod extrema_pc_extended;
mod extrema_pc_ellipse;
mod extrema_pc_hyperbola;
mod extrema_pc_line;
mod extrema_pc_offset;
mod extrema_pc2d_conics;
mod extrema_pc2d_curves;
mod extrema_pc2d_general;
mod extrema_pc_parabola;
mod extrema_pc_search_mode;
mod convert;
mod gcpnts;
mod extrema_ext_pc;
mod geomlib;
mod geomlib_check_curve_on_surface;
mod gprop;
mod intana;
mod lprop;
mod project;

use kentron::math::{Point, Point2, Vector};

#[track_caller]
pub fn assert_point_near(actual: Point, expected: Point, tolerance: f64) {
    assert!(
        actual.distance(expected) <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

#[track_caller]
pub fn assert_point2_near(actual: Point2, expected: Point2, tolerance: f64) {
    assert!(
        actual.distance(expected) <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

#[track_caller]
pub fn assert_vector_near(actual: Vector, expected: Vector, tolerance: f64) {
    assert!(
        (actual - expected).length() <= tolerance,
        "{actual:?} is not within {tolerance} of {expected:?}"
    );
}

#[track_caller]
pub fn assert_near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not within {tolerance} of {expected}"
    );
}

pub fn p(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

pub fn p2(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// Checks the closest extremum's parameter and squared distance (`None` skips a check).
#[track_caller]
pub fn assert_min<P>(
    result: &kentron::geometry::extrema::PointCurveExtrema<P>,
    parameter: Option<f64>,
    square_distance: Option<f64>,
    tolerance: f64,
) {
    let min = result.min().expect("no extrema");
    if let Some(t) = parameter {
        assert_near(min.parameter, t, tolerance);
    }
    if let Some(d) = square_distance {
        assert_near(min.square_distance, d, tolerance);
    }
}

/// A planar 3D curve and its 2D counterpart give the same extrema for a point above the
/// plane (`(8, 3, 6)` vs `(8, 3)`), with squared distances differing by the height².
#[track_caller]
pub fn assert_planar_matches_2d(
    c3: &kentron::geometry::curve::Curve,
    c2: &kentron::geometry::curve2::Curve2,
    tolerance: f64,
) {
    use kentron::geometry::extrema::PointCurveSearch;
    let r3 = PointCurveSearch::new(c3).find(p(8., 3., 6.), tolerance).unwrap();
    let r2 = PointCurveSearch::new(c2).find(p2(8., 3.), tolerance).unwrap();
    assert_eq!(r3.len(), r2.len());
    for (e3, e2) in r3.iter().zip(&r2) {
        assert_near(e3.parameter, e2.parameter, tolerance);
        assert_near(e3.square_distance, e2.square_distance + 36., tolerance);
    }
}

/// A small deterministic pseudo-random generator (xorshift64*) for the randomized OCCT
/// tests; the sequences differ from OCCT's `std::mt19937`, the properties tested do not.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    /// Uniform in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let bits = self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11;
        lo + (hi - lo) * (bits as f64 / (1u64 << 53) as f64)
    }
}

/// The smallest and largest distance from `q` to its orthogonal projections onto `curve`
/// (kentron's other point/curve query, used as the reference where OCCT compares with its
/// older `Extrema_ExtPC`).
pub fn projection_distances(
    q: Point,
    curve: &kentron::geometry::curve::Curve,
) -> Option<(f64, f64)> {
    let projections = kentron::geometry::extrema::point_curve(q, curve);
    Some((projections.first()?.distance, projections.last()?.distance))
}
