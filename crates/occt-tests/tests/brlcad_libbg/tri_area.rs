//! `src/libbg/tests/tri_area.c` (`bg_tri_area`): a triangle's area stays positive far from
//! the origin (up to 1e16) and for tiny sides (down to 1e-15).

use crate::p;
use kentron::math::triangle_area;

#[test]
fn bg_tri_area() {
    let mut base: f64 = 1.0;
    let mut area = 1.0;
    while area > 0.0 {
        base *= 10.0;
        area = triangle_area(p(base, base, 0.0), p(base + 1.0, base, 0.0), p(base + 0.5, base, 1.0));
    }
    assert!(base >= 1e16, "the area vanished at {base}");
    // (From BRL-CAD's VUNITIZE_TOL down.)
    let mut base: f64 = 1e-7;
    area = 1.0;
    while area > 0.0 {
        base *= 0.5;
        area = triangle_area(p(10.0, 10.0, 0.0), p(10.0 + base, 10.0, 0.0), p(10.0 + 0.5 * base, 10.0, 5.0));
    }
    assert!(base <= 1e-15, "the area vanished at {base}");
}
