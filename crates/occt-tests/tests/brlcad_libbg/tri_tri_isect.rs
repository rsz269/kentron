//! `src/libbg/tests/tri_tri_isect.c` (`bg_tri_tri_isect`): whether two triangles meet.

use crate::p;
use kentron::math::{triangles_intersect, Point};
use kentron::tolerance::LINEAR;

const T: [Point; 3] = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)];

#[test]
fn bg_tri_tri_isect_cp_null() { assert!(!triangles_intersect(T, [p(1.00001, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.00001, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_cp_vertex() { assert!(triangles_intersect(T, [p(1.00001, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.0, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_cp_edge() { assert!(triangles_intersect(T, [p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.0, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_cp_full() { assert!(triangles_intersect(T, [p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(0.7, 0.8, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_null() { assert!(!triangles_intersect(T, [p(1.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(1.0, 1.0, 1.0)], LINEAR)); }

/// Not one of BRL-CAD's cases: triangles crossing (each through the other's inside), and one
/// touching another's inside with a corner, for the non-coplanar paths.
#[test]
fn crossing_and_touching_triangles() {
    let pierced = [p(0.2, 0.1, -1.0), p(0.8, 0.1, -1.0), p(0.5, 0.5, 1.0)];
    assert!(triangles_intersect(T, pierced, LINEAR));
    let touching = [p(0.6, 0.3, 0.0), p(0.6, 0.3, 1.0), p(1.6, 0.3, 1.0)];
    assert!(triangles_intersect(T, touching, LINEAR));
    let above = [p(0.6, 0.3, 1e-5), p(0.6, 0.3, 1.0), p(1.6, 0.3, 1.0)];
    assert!(!triangles_intersect(T, above, LINEAR));
}
