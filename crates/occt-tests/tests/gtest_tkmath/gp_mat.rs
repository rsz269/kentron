//! gp_Mat_Test.cxx

use occt_tests::assert_double_eq;
use kentron::math::Matrix3;

/// OCC22595: the default matrix is all zeros.
#[test]
fn occ22595_default_constructor() {
    let m = Matrix3::default();
    for row in 0..3 {
        for col in 0..3 {
            assert_double_eq!(0.0, m[(row, col)]);
        }
    }
}
