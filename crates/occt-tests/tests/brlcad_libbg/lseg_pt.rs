//! `src/libbg/tests/lseg_pt.c` (`bg_lseg_pt_dist`): the point of a segment nearest a point,
//! and their distance.

use crate::p;
use kentron::math::{closest_point_on_segment, Point};

fn check(pt: Point, pca: Point, distance: f64) {
    let (a, b) = (p(-5.0, -5.0, -5.0), p(5.0, 5.0, 5.0));
    let q = closest_point_on_segment(pt, a, b);
    assert!(q.distance(pca) <= 1e-4, "{q:?} != {pca:?}");
    assert!((pt.distance(q) - distance).abs() <= 1e-4, "{} != {distance}", pt.distance(q));
}

const THIRD: f64 = -10.0 / 3.0;

#[test]
fn bg_lseg_pt_dist_case0() { check(p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), 0.0); }
#[test]
fn bg_lseg_pt_dist_case1() { check(p(-5.0, -5.0, -5.0), p(-5.0, -5.0, -5.0), 0.0); }
#[test]
fn bg_lseg_pt_dist_case2() { check(p(5.0, 5.0, 5.0), p(5.0, 5.0, 5.0), 0.0); }
#[test]
fn bg_lseg_pt_dist_case3() { check(p(-10.0, -5.0, -10.0), p(-5.0, -5.0, -5.0), 7.07107); }
#[test]
fn bg_lseg_pt_dist_case4() { check(p(5.0, 10.0, 10.0), p(5.0, 5.0, 5.0), 7.07107); }
#[test]
fn bg_lseg_pt_dist_case5() { check(p(0.0, -5.0, -5.0), p(THIRD, THIRD, THIRD), 4.08248); }
