//! `src/libbg/tests/plane_pt.c` (`bg_test plane_pt`): collinear and distinct points, distances
//! from lines and segments, planes through three points, points on planes, coplanar points.

use crate::{p, BN_TOL_DIST};
use kentron::geometry::construct::ConstructionError;
use kentron::geometry::curve::Line;
use kentron::geometry::surface::Plane;
use kentron::geometry::ParametricSurface;
use kentron::math::{closest_point_on_segment, Direction, Point, Vector};
use kentron::tolerance::LINEAR;

/// `bg_3pts_collinear`: no plane passes through the points alone.
fn collinear(a: Point, b: Point, c: Point) -> bool {
    matches!(Plane::through_points(a, b, c), Err(ConstructionError::ColinearPoints | ConstructionError::ConfusedPoints))
}

#[test]
fn bg_plane_3pts_collinear_false_1() { assert!(!collinear(p(0.0, 0.0, 0.0), p(1.0, 1.0, 1.0), p(6.0, 3.0, 8.0))); }
#[test]
fn bg_plane_3pts_collinear_false_2() { assert!(!collinear(p(0.0, 0.0, 0.0), p(6.0, 3.0, 8.0), p(1.0, 1.0, 1.0))); }
#[test]
fn bg_plane_3pts_collinear_false_3() { assert!(!collinear(p(1.0, 1.0, 1.0), p(0.0, 0.0, 0.0), p(6.0, 3.0, 8.0))); }
#[test]
fn bg_plane_3pts_collinear_false_4() { assert!(!collinear(p(1.0, 1.0, 1.0), p(6.0, 3.0, 8.0), p(0.0, 0.0, 0.0))); }
#[test]
fn bg_plane_3pts_collinear_false_5() { assert!(!collinear(p(6.0, 3.0, 8.0), p(1.0, 1.0, 1.0), p(0.0, 0.0, 0.0))); }
#[test]
fn bg_plane_3pts_collinear_false_6() { assert!(!collinear(p(6.0, 3.0, 8.0), p(0.0, 0.0, 0.0), p(1.0, 1.0, 1.0))); }
#[test]
fn bg_plane_3pts_collinear_same() { assert!(collinear(p(5.3, 8.1, 9.6), p(5.3, 8.1, 9.6), p(5.3, 8.1, 9.6))); }
#[test]
fn bg_plane_3pts_collinear_partial_same_1() { assert!(collinear(p(6.0, 2.0, 4.0), p(6.0, 2.0, 4.0), p(9.3, 8.1, 4.6))); }
#[test]
fn bg_plane_3pts_collinear_partial_same_2() { assert!(collinear(p(9.3, 8.1, 4.6), p(6.0, 2.0, 4.0), p(6.0, 2.0, 4.0))); }
#[test]
fn bg_plane_3pts_collinear_partial_same_3() { assert!(collinear(p(6.0, 2.0, 4.0), p(9.3, 8.1, 4.6), p(6.0, 2.0, 4.0))); }

const L1: [f64; 3] = [1.8, 4.6, 2.2];
const L2: [f64; 3] = [5.4, 13.8, 6.6];
const L3: [f64; 3] = [16.2, 41.4, 19.8];
fn pa(c: [f64; 3]) -> Point { p(c[0], c[1], c[2]) }

#[test]
fn bg_plane_3pts_collinear_true_1() { assert!(collinear(pa(L1), pa(L2), pa(L3))); }
#[test]
fn bg_plane_3pts_collinear_true_2() { assert!(collinear(pa(L1), pa(L3), pa(L2))); }
#[test]
fn bg_plane_3pts_collinear_true_3() { assert!(collinear(pa(L2), pa(L1), pa(L3))); }
#[test]
fn bg_plane_3pts_collinear_true_4() { assert!(collinear(pa(L2), pa(L3), pa(L1))); }
#[test]
fn bg_plane_3pts_collinear_true_5() { assert!(collinear(pa(L3), pa(L1), pa(L2))); }
#[test]
fn bg_plane_3pts_collinear_true_6() { assert!(collinear(pa(L3), pa(L2), pa(L1))); }

/// `bg_3pts_distinct`: no two of the points coincide.
fn distinct(a: Point, b: Point, c: Point) -> bool {
    a.distance(b) > LINEAR && b.distance(c) > LINEAR && a.distance(c) > LINEAR
}

const D1: [f64; 3] = [5.3, 8.6, 9.1];
const D2: [f64; 3] = [8.7, 4.1, 9.3];
const D3: [f64; 3] = [7.5, 9.2, 4.1];

#[test]
fn bg_plane_3pts_distinct_false_all() { assert!(!distinct(pa(D1), pa(D1), pa(D1))); }
#[test]
fn bg_plane_3pts_distinct_false_1() { assert!(!distinct(pa(D1), pa(D1), pa(D2))); }
#[test]
fn bg_plane_3pts_distinct_false_2() { assert!(!distinct(pa(D1), pa(D2), pa(D1))); }
#[test]
fn bg_plane_3pts_distinct_false_3() { assert!(!distinct(pa(D2), pa(D1), pa(D1))); }
#[test]
fn bg_plane_3pts_distinct_true_1() { assert!(distinct(pa(D2), pa(D3), pa(D1))); }
#[test]
fn bg_plane_3pts_distinct_true_2() { assert!(distinct(pa(D2), pa(D1), pa(D3))); }
#[test]
fn bg_plane_3pts_distinct_true_3() { assert!(distinct(pa(D1), pa(D2), pa(D3))); }
#[test]
fn bg_plane_3pts_distinct_true_4() { assert!(distinct(pa(D1), pa(D3), pa(D2))); }
#[test]
fn bg_plane_3pts_distinct_true_5() { assert!(distinct(pa(D3), pa(D2), pa(D1))); }
#[test]
fn bg_plane_3pts_distinct_true_6() { assert!(distinct(pa(D3), pa(D1), pa(D2))); }

/// `bg_distsq_line3_pt3`: the squared distance from the line through `pt` along `dir`.
fn line_distance_squared(pt: Point, dir: Vector, a: Point, expected: f64) {
    let line = Line::new(pt, Direction::from_vector(dir));
    let d = line.distance_to_point(a);
    assert!((d * d - expected).abs() <= 1e-5, "{} != {expected}", d * d);
}

#[test]
fn bg_plane_distsq_line3_pt3_1() { line_distance_squared(p(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0), p(0.0, 0.0, 0.0), 0.0); }
#[test]
fn bg_plane_distsq_line3_pt3_2() { line_distance_squared(p(-1.0, -1.0, -1.0), Vector::new(1.0, 1.0, 1.0), p(0.0, 0.0, 0.0), 0.0); }
#[test]
fn bg_plane_distsq_line3_pt3_3() { line_distance_squared(p(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0), p(0.0, 5.0, 0.0), 16.666666); }
#[test]
fn bg_plane_distsq_line3_pt3_4() { line_distance_squared(p(5.0, 6.0, 4.0), Vector::new(9.0, 3.0, 8.0), p(-5.0, 2.0, 4.0), 48.441558); }

/// Where on segment `ab` the point nearest `pt` lies (BRL-CAD's cases of
/// `bg_distsq_pt3_lseg3_v2`): 0 `pt` on the segment's inside, 1 at `a`, 2 at `b`, 3 nearest
/// `a`, 4 nearest `b`, 5 nearest a point inside.
pub fn segment_case(a: Point, b: Point, pt: Point, tolerance: f64) -> (u8, Point) {
    let q = closest_point_on_segment(pt, a, b);
    let case = if pt.distance(a) <= tolerance {
        1
    } else if pt.distance(b) <= tolerance {
        2
    } else if pt.distance(q) <= tolerance {
        0
    } else if q.distance(a) <= tolerance {
        3
    } else if q.distance(b) <= tolerance {
        4
    } else {
        5
    };
    (case, q)
}

fn segment_distance_squared(pt: Point, expected_case: u8, expected: f64) {
    let (a, b) = (p(-5.0, -5.0, -5.0), p(5.0, 5.0, 5.0));
    let (case, q) = segment_case(a, b, pt, BN_TOL_DIST);
    assert_eq!(case, expected_case);
    let d = pt.distance(q);
    assert!((d * d - expected).abs() <= 1e-5, "{} != {expected}", d * d);
}

#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case0_1() { segment_distance_squared(p(0.0, 0.0, 0.0), 0, 0.0); }
#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case1_1() { segment_distance_squared(p(-5.0, -5.0, -5.0), 1, 0.0); }
#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case2_1() { segment_distance_squared(p(5.0, 5.0, 5.0), 2, 0.0); }
#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case3_1() { segment_distance_squared(p(-10.0, -5.0, -10.0), 3, 50.0); }
#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case4_1() { segment_distance_squared(p(5.0, 10.0, 10.0), 4, 50.0); }
#[test]
fn bg_plane_distsq_pt3_lseg3_v2_case5_1() { segment_distance_squared(p(0.0, -5.0, -5.0), 5, 16.666666); }

/// `bg_make_plane_3pts`: a plane through the points, passing through all three.
fn plane_through(a: Point, b: Point, c: Point, succeeds: bool) {
    match Plane::through_points(a, b, c) {
        Ok(plane) => {
            assert!(succeeds, "a plane through points that should not define one");
            for q in [a, b, c] {
                assert!(plane.distance(q) <= 1e-12 * 30.0, "{q:?} off the plane");
            }
        }
        Err(_) => assert!(!succeeds, "no plane through points that define one"),
    }
}

#[test]
fn bg_plane_mk_plane_3pts_fail_same() { plane_through(p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), false); }
#[test]
fn bg_plane_mk_plane_3pts_fail_collinear_1() { plane_through(p(0.0, 0.0, 0.0), p(5.0, 8.0, 2.0), p(0.0, 0.0, 0.0), false); }
#[test]
fn bg_plane_mk_plane_3pts_fail_collinear_2() { plane_through(p(5.0, 8.0, 2.0), p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), false); }
#[test]
fn bg_plane_mk_plane_3pts_fail_collinear_3() { plane_through(p(-10.0, -10.0, -10.0), p(0.0, 0.0, 0.0), p(5.0, 5.0, 5.0), false); }
#[test]
fn bg_plane_mk_plane_3pts_succeed() { plane_through(p(-10.0, -10.0, -10.0), p(7.0, 3.0, 8.0), p(19.0, 4.0, 2.0), true); }

/// `bg_plane_closest_pt` then `bg_plane_pt_at`: a point of the plane at distance `w` from
/// the origin along the unit normal of `N` (BRL-CAD's plane), found again from its
/// parameters.
fn plane_round_trip(n: [f64; 4], pt: Point) {
    let normal = Vector::new(n[0], n[1], n[2]);
    let origin = Point::ORIGIN + normal * (n[3] / normal.length());
    let plane = Plane::new(origin, Direction::from_vector(normal));
    let (u, v) = plane.parameters(pt);
    assert!(plane.point(u, v).distance(pt) <= 1e-9 * (1.0 + (pt - Point::ORIGIN).length()), "{:?} != {pt:?}", plane.point(u, v));
}

#[test]
fn bg_plane_pt_project_1() { plane_round_trip([0.0, 0.0, 1.0, 1.0], p(2.0, 2.0, 1.0)); }
#[test]
fn bg_plane_pt_project_2() {
    plane_round_trip([0.737397, 0.461098, 0.493593, 236.3], p(77.526016919831704, 185.674010580054585, 189.465011325663767));
}

/// `bg_coplanar_pts`: the fourth point on the plane of the first three, within BRL-CAD's
/// default tolerance.
fn coplanar(q: [Point; 4]) -> bool {
    let plane = Plane::through_points(q[0], q[1], q[2]).expect("a plane");
    plane.distance(q[3]) <= BN_TOL_DIST
}

const C1: [f64; 3] = [-3709.082985807021941582207, 3510.525439930336233373964, 3905.461506563813145476161];
const C2: [f64; 3] = [-3152.697964017485901422333, 2183.289655637726355053019, 5294.303067153545271139592];
const C3: [f64; 3] = [-5625.550060079097420384642, 3027.217726896207295794738, 4211.349734112171063316055];

#[test]
fn bg_plane_pts_coplanar_1() { assert!(coplanar([p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 2.0, 0.0)])); }
#[test]
fn bg_plane_pts_coplanar_2() { assert!(!coplanar([p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 2.0, 1.0)])); }
#[test]
fn bg_plane_pts_coplanar_3() {
    assert!(coplanar([pa(C1), pa(C2), pa(C3), p(-5069.165038289564108708873, 1699.981942603596962726442, 5600.191294701904553221539)]));
}
#[test]
fn bg_plane_pts_coplanar_4() { assert!(coplanar([pa(C1), pa(C2), pa(C3), p(-4019.727, 2883.893, 4507.108)])); }
#[test]
fn bg_plane_pts_coplanar_5() { assert!(!coplanar([pa(C1), pa(C2), pa(C3), p(-4000.0, 2000.0, 4000.0)])); }
