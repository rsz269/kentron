//! `src/libbg/tests/tri_tri_isect_coplanar.cpp` (`bg_tri_tri_isect_coplanar`): whether two
//! triangles in one plane meet. BRL-CAD's variants requiring an overlap of non-zero area
//! (`*_area`) have no kentron counterpart (see GTEST_PORTING.md).

use crate::p;
use kentron::math::{triangles_intersect, Point};
use kentron::tolerance::LINEAR;

const T: [Point; 3] = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)];

#[test]
fn bg_tri_tri_isect_coplanar_null_noarea() { assert!(!triangles_intersect(T, [p(1.00001, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.00001, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_coplanar_vertex_noarea() { assert!(triangles_intersect(T, [p(1.00001, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.0, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_coplanar_edge_noarea() { assert!(triangles_intersect(T, [p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(1.0, 1.0, 0.0)], LINEAR)); }
#[test]
fn bg_tri_tri_isect_coplanar_full_noarea() { assert!(triangles_intersect(T, [p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(0.7, 0.8, 0.0)], LINEAR)); }
