//! `src/libbg/tests/tri_closest_pt.c` (`bg_tri_closest_pt`): the distance from a point to a
//! triangle.

use crate::p;
use kentron::math::{closest_point_on_triangle, Point};

fn distance(pt: Point) -> f64 {
    let (a, b, c) = (p(-1.0, -1.0, 0.0), p(1.0, -1.0, 0.0), p(0.0, 1.0, 0.0));
    closest_point_on_triangle(pt, a, b, c).distance(pt)
}

#[test]
fn bg_tri_pt_dist_coplanar_center() { assert!(distance(p(0.0, 0.0, 0.0)) <= 1e-12); }
#[test]
fn bg_tri_pt_dist_coplanar_vert_closest() { assert!((distance(p(-2.0, -1.0, 0.0)) - 1.0).abs() <= 1e-12); }
