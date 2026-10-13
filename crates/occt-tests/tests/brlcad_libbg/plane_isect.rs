//! `src/libbg/tests/plane_isect.c` (`bg_test plane_isect`): a line against a segment, and two
//! segments, each case of BRL-CAD's return codes asserted geometrically through kentron's
//! closest points of segments.

use crate::{p, BN_TOL_DIST};
use kentron::math::{closest_points_on_segments, Point, Vector};

/// How far lines are extended to find where they meet with segments' closest points.
const REACH: f64 = 1e4;

/// The parameters (along `a + s·u` and `b + t·v`) where two non-parallel lines come nearest,
/// and their distance there.
fn lines_meet(a: Point, u: Vector, b: Point, v: Vector) -> (f64, f64, f64) {
    let (s, t) = closest_points_on_segments(a - u * REACH, a + u * REACH, b - v * REACH, b + v * REACH);
    let (s, t) = (-REACH + s * 2.0 * REACH, -REACH + t * 2.0 * REACH);
    (s, t, (a + u * s).distance(b + v * t))
}

fn parallel(u: Vector, v: Vector) -> bool {
    u.cross(v).length() <= 1e-12 * u.length() * v.length()
}

/// BRL-CAD's `bg_isect_line_lseg` cases for the line `pt + t·d` and the segment `ab`, with
/// `t` in units of length: -4 `a` and `b` coincide, -3/-2 the line crosses the segment's line
/// before `a`/after `b`, -1 no crossing, 0 along the segment's line (`t` at `a`), 1 at `a`,
/// 2 at `b`, 3 between them.
fn line_segment(pt: Point, d: Vector, a: Point, b: Point) -> (i8, f64) {
    let unit = d * (1.0 / d.length());
    let ab = b - a;
    if ab.length() <= BN_TOL_DIST {
        return (-4, 0.0);
    }
    if parallel(unit, ab) {
        let off = (a - pt) - unit * (a - pt).dot(unit);
        return if off.length() <= BN_TOL_DIST { (0, (a - pt).dot(unit)) } else { (-1, 0.0) };
    }
    let (t, s, gap) = lines_meet(pt, unit, a, ab);
    if gap > BN_TOL_DIST {
        return (-1, 0.0);
    }
    let along = s * ab.length();
    let case = if along.abs() <= BN_TOL_DIST {
        1
    } else if (along - ab.length()).abs() <= BN_TOL_DIST {
        2
    } else if s < 0.0 {
        -3
    } else if s > 1.0 {
        -2
    } else {
        3
    };
    (case, t)
}

fn check_line(a: Point, b: Point, case: i8, t: f64) {
    let (found, ft) = line_segment(p(-5.0, -5.0, -5.0), Vector::new(1.0, 1.0, 1.0), a, b);
    assert_eq!(found, case);
    if case >= -3 && case != -1 {
        assert!((ft - t).abs() <= 1e-3, "{ft} != {t}");
    }
}

#[test]
fn bg_plane_isect_line_lseg_case_minus4_1() { check_line(p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), -4, 0.0); }
#[test]
fn bg_plane_isect_line_lseg_case_minus3_1() { check_line(p(0.0, 0.0, 1.0), p(0.0, 0.0, 5.0), -3, 8.660); }
#[test]
fn bg_plane_isect_line_lseg_case_minus2_1() { check_line(p(0.0, 0.0, -5.0), p(0.0, 0.0, -1.0), -2, 8.660); }
#[test]
fn bg_plane_isect_line_lseg_case_minus1_1() { check_line(p(-5.0, -5.0, -6.0), p(0.0, 0.0, -1.0), -1, 0.0); }
#[test]
fn bg_plane_isect_line_lseg_case0_1() { check_line(p(1.0, 1.0, 1.0), p(5.0, 5.0, 5.0), 0, 10.392); }
#[test]
fn bg_plane_isect_line_lseg_case1_1() { check_line(p(0.0, 0.0, 0.0), p(0.0, 0.0, 5.0), 1, 8.660); }
#[test]
fn bg_plane_isect_line_lseg_case2_1() { check_line(p(0.0, 0.0, -5.0), p(0.0, 0.0, 0.0), 2, 8.660); }
#[test]
fn bg_plane_isect_line_lseg_case3_1() { check_line(p(0.0, 0.0, -5.0), p(0.0, 0.0, 5.0), 3, 8.660); }

/// BRL-CAD's `bg_isect_lseg3_lseg3` cases for the segments `p + s·pdir` and `q + t·qdir`
/// (`s`, `t` from 0 to 1): -3 the lines cross outside a segment, -2 parallel, -1 on one line
/// apart, 0 on one line overlapping, 1 crossing; with the parameters where the lines cross
/// (on one line: of `q` and of `q + qdir` along the first).
fn segments(pp: Point, pdir: Vector, q: Point, qdir: Vector) -> (i8, f64, f64) {
    if parallel(pdir, qdir) {
        let off = (q - pp) - pdir * ((q - pp).dot(pdir) / pdir.dot(pdir));
        if off.length() > BN_TOL_DIST {
            return (-2, 0.0, 0.0);
        }
        let along = |x: Point| (x - pp).dot(pdir) / pdir.dot(pdir);
        let (a, b) = (along(q), along(q + qdir));
        let (lo, hi) = (a.min(b), a.max(b));
        let slack = BN_TOL_DIST / pdir.length();
        let overlap = hi >= -slack && lo <= 1.0 + slack;
        return (if overlap { 0 } else { -1 }, a, b);
    }
    let (s, t, gap) = lines_meet(pp, pdir, q, qdir);
    if gap > BN_TOL_DIST {
        return (-3, s, t);
    }
    let inside = |x: f64, len: f64| x >= -BN_TOL_DIST / len && x <= 1.0 + BN_TOL_DIST / len;
    (if inside(s, pdir.length()) && inside(t, qdir.length()) { 1 } else { -3 }, s, t)
}

fn check_segments(pp: Point, pdir: Vector, q: Point, qdir: Vector, case: i8, d0: f64, d1: f64) {
    let (found, s, t) = segments(pp, pdir, q, qdir);
    assert_eq!(found, case);
    if case != -2 {
        assert!((s - d0).abs() <= 1e-5 && (t - d1).abs() <= 1e-5, "({s}, {t}) != ({d0}, {d1})");
    }
}

fn v(x: f64, y: f64, z: f64) -> Vector {
    Vector::new(x, y, z)
}

#[test]
fn bg_plane_isect_lseg3_lseg3_case_minus3_1() { check_segments(p(0.0, 0.0, 0.0), v(0.0, 1.0, 1.0), p(0.0, 0.0, 1.0), v(0.0, -1.0, 1.0), -3, 0.5, -0.5); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case_minus3_2() { check_segments(p(0.0, 0.0, 0.0), v(0.0, 1.0, 1.0), p(0.0, 1.0, 0.0), v(0.0, -0.25, 0.25), -3, 0.5, 2.0); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case_minus2_1() { check_segments(p(0.0, 0.0, 0.0), v(0.0, 1.0, 1.0), p(1.0, 1.0, 1.0), v(0.0, 1.0, 1.0), -2, 0.0, 0.0); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case_minus1_1() { check_segments(p(-5.0, -5.0, -5.0), v(3.0, 3.0, 3.0), p(0.0, 0.0, 0.0), v(-1.0, -1.0, -1.0), -1, 1.666667, 1.333333); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case0_1() { check_segments(p(-5.0, -5.0, -5.0), v(3.0, 3.0, 3.0), p(0.0, 0.0, 0.0), v(-2.0, -2.0, -2.0), 0, 1.666667, 1.0); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case0_2() { check_segments(p(-5.0, -5.0, -5.0), v(3.0, 3.0, 3.0), p(0.0, 0.0, 0.0), v(-3.0, -3.0, -3.0), 0, 1.666667, 0.666667); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case0_3() { check_segments(p(-5.0, -5.0, -5.0), v(5.0, 5.0, 5.0), p(-4.0, -4.0, -4.0), v(2.0, 2.0, 2.0), 0, 0.2, 0.6); }
#[test]
fn bg_plane_isect_lseg3_lseg3_case1_1() { check_segments(p(0.0, 0.0, 0.0), v(0.0, 1.0, 1.0), p(0.0, 1.0, 0.0), v(0.0, -1.0, 1.0), 1, 0.5, 0.5); }
