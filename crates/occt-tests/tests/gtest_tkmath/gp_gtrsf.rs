//! gp_GTrsf_Test.cxx (gp_GTrsf is `AffineTransform`)

use occt_tests::assert_near;
use kentron::math::{AffineTransform, Direction, Point};
use kentron::tolerance;

/// OCC13963: an affinity with ratio 1 is the identity.
#[test]
fn occ13963_set_affinity_ratio1_is_identity() {
    let t = AffineTransform::plane_affinity(Point::ORIGIN, Direction::new(1.0, -1.0, 0.0), 1.0);
    let result = t.transform_point(Point::new(1.0, 0.0, 0.0));
    assert_near!(result.x, 1.0, tolerance::LINEAR);
    assert_near!(result.y, 0.0, tolerance::LINEAR);
    assert_near!(result.z, 0.0, tolerance::LINEAR);
}

/// OCC13963: (1,0,0) with the plane of normal (1,-1,0)/sqrt(2) and ratio 2 gives
/// (1.5,-0.5,0).
#[test]
fn occ13963_set_affinity_ratio2_scales_along_normal() {
    let t = AffineTransform::plane_affinity(Point::ORIGIN, Direction::new(1.0, -1.0, 0.0), 2.0);
    let result = t.transform_point(Point::new(1.0, 0.0, 0.0));
    assert_near!(result.x, 1.5, tolerance::LINEAR);
    assert_near!(result.y, -0.5, tolerance::LINEAR);
    assert_near!(result.z, 0.0, tolerance::LINEAR);
}
