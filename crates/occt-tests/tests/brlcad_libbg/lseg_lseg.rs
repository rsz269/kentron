//! `src/libbg/tests/lseg_lseg.c` (`bg_lseg_lseg_dist`): the closest points of two segments.
//! BRL-CAD's test prints them; here they are checked against a dense search.

use crate::p;
use kentron::math::closest_points_on_segments;

#[test]
fn bg_lseg_lseg_1() {
    let (p0, p1) = (p(3.41819463201248075, 7.64533159204993584, 22.99994836139265786), p(3.52252929469979303, 7.79661877760347455, 22.99995349026097102));
    let (q0, q1) = (p(3.44524589136147830, 7.67444957869714628, 22.91984784277452647), p(3.56555225936148323, 7.98564760063074353, 23.37334866995422900));
    let (s, t) = closest_points_on_segments(p0, p1, q0, q1);
    assert!((0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t));
    let distance = (p0 + (p1 - p0) * s).distance(q0 + (q1 - q0) * t);
    let n = 2000;
    let mut best = f64::INFINITY;
    for i in 0..=n {
        let a = p0 + (p1 - p0) * (i as f64 / n as f64);
        for j in 0..=n {
            best = best.min(a.distance(q0 + (q1 - q0) * (j as f64 / n as f64)));
        }
    }
    assert!(distance <= best + 1e-12, "{distance} > {best}");
    assert!(best - distance <= 1e-4, "the search found much nearer points: {best} vs {distance}");
}
