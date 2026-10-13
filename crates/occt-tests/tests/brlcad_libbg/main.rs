//! Hand translations of BRL-CAD's libbg tests (`src/libbg/tests`): one module per test source,
//! each BRL-CAD test (named in its CMakeLists.txt) with a counterpart using kentron's API.
//! BRL-CAD's return codes become geometric assertions; see `GTEST_PORTING.md` for the tests
//! with no kentron counterpart.

mod lseg_lseg;
mod lseg_pt;
mod plane_dist;
mod plane_isect;
mod plane_pt;
mod tri_area;
mod tri_closest_pt;
mod tri_tri_isect;
mod tri_tri_isect_coplanar;
mod trimesh_area;

use kentron::math::Point;

/// A point from BRL-CAD's `x,y,z` arguments.
pub fn p(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

/// BRL-CAD's default distance tolerance (`BN_TOL_INIT_TOL`), where a test relies on it.
pub const BN_TOL_DIST: f64 = 0.0005;
