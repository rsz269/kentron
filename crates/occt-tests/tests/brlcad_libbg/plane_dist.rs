//! `src/libbg/tests/plane_dist.c` (`bg_test plane_dist`): the distance from a point to a
//! segment, with its closest point, and between points.

use crate::plane_pt::segment_case;
use crate::{p, BN_TOL_DIST};
use kentron::math::Point;

/// `bg_dist_pnt3_lseg3`: where on `ab` the point nearest `pt` lies, the distance to it (for a
/// point on the segment, BRL-CAD's parameter along it instead) and the point.
fn segment(a: Point, b: Point, pt: Point, case: u8, distance: f64, pca: Point) {
    let (found, q) = segment_case(a, b, pt, BN_TOL_DIST);
    assert_eq!(found, case);
    let d = if case == 0 { q.distance(a) / b.distance(a) } else { pt.distance(q) };
    assert!((d - distance).abs() <= 1e-4, "{d} != {distance}");
    assert!(q.distance(pca) <= 1e-4, "{q:?} != {pca:?}");
}

/// The closest point and distance only: BRL-CAD's three `pt3_line3` cases run the segment
/// test with return codes of another function, which its driver accepts because it fails
/// only when every value differs (see GTEST_PORTING.md).
fn segment_point(a: Point, b: Point, pt: Point, distance: f64, pca: Point) {
    let (_, q) = segment_case(a, b, pt, BN_TOL_DIST);
    assert!((pt.distance(q) - distance).abs() <= 1e-4, "{} != {distance}", pt.distance(q));
    assert!(q.distance(pca) <= 1e-4, "{q:?} != {pca:?}");
}

const A: Point = Point::new(-5.0, -5.0, -5.0);
const B: Point = Point::new(5.0, 5.0, 5.0);
const THIRD: f64 = -10.0 / 3.0;

#[test]
fn bg_plane_dist_pt3_line3_case0_1() { segment_point(A, p(1.0, 1.0, 1.0), A, 0.0, A); }
#[test]
fn bg_plane_dist_pt3_line3_case1_1() { segment_point(A, p(1.0, 1.0, 1.0), p(0.0, 0.0, 0.0), 0.0, p(0.0, 0.0, 0.0)); }
#[test]
fn bg_plane_dist_pt3_line3_case2_1() { segment_point(A, B, p(0.0, -5.0, -5.0), 4.08248, p(THIRD, THIRD, THIRD)); }

#[test]
fn bg_plane_dist_pt3_lseg3_case0_1() { segment(A, B, p(0.0, 0.0, 0.0), 0, 0.5, p(0.0, 0.0, 0.0)); }
#[test]
fn bg_plane_dist_pt3_lseg3_case1_1() { segment(A, B, A, 1, 0.0, A); }
#[test]
fn bg_plane_dist_pt3_lseg3_case2_1() { segment(A, B, B, 2, 0.0, B); }
#[test]
fn bg_plane_dist_pt3_lseg3_case3_1() { segment(A, B, p(-10.0, -5.0, -10.0), 3, 7.07107, A); }
#[test]
fn bg_plane_dist_pt3_lseg3_case4_1() { segment(A, B, p(5.0, 10.0, 10.0), 4, 7.07107, B); }
#[test]
fn bg_plane_dist_pt3_lseg3_case5_1() { segment(A, B, p(0.0, -5.0, -5.0), 5, 4.08248, p(THIRD, THIRD, THIRD)); }

/// `bg_dist_pnt3_pnt3`.
#[test]
fn bg_plane_dist_pt3_pt3_same() { assert_eq!(p(0.0, 0.0, 0.0).distance(p(0.0, 0.0, 0.0)), 0.0); }
#[test]
fn bg_plane_dist_pt3_pt3_1() { assert!((p(0.0, 0.0, 0.0).distance(p(1.0, 1.0, 1.0)) - 1.732051).abs() <= 1e-6); }
