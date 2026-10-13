//! gp_Torus_Test.cxx

use kentron::geometry::ParametricSurface;
use std::f64::consts::PI;

use occt_tests::assert_near;
use kentron::geometry::surface::Torus;
use kentron::math::{Direction, Frame, Point};

/// OCC26746: the implicit equation's coefficients vanish on the torus.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn occ26746_coefficients_correct() {
    // from the original DRAW script:
    // torus tr 55.52514413 2.070076585 73.83409062
    //          0.37231784651136368 0.58886674834874120 0.71736697293527607
    //          0.80682335496555135 0.17666016102759910 -0.56376170618524390
    //          87.08479625 23.14682176
    let frame = Frame::new(
        Point::new(55.52514413, 2.070076585, 73.83409062),
        Direction::new(0.37231784651136368, 0.58886674834874120, 0.71736697293527607),
        Direction::new(0.80682335496555135, 0.17666016102759910, -0.56376170618524390),
    );
    let torus = Torus::new(frame, 87.08479625, 23.14682176);

    let tolerance = 3.0e-7;
    let nb_pts_max = 5;
    let step = 2.0 * PI / nb_pts_max as f64;
    let coeffs = torus.implicit_coefficients();

    let mut u = 0.0;
    for u_ind in 0..=nb_pts_max {
        let mut v = 0.0;
        for v_ind in 0..=nb_pts_max {
            let p = torus.point(u, v);
            let (x1, y1, z1) = (p.x, p.y, p.z);
            let (x2, y2, z2) = (x1 * x1, y1 * y1, z1 * z1);
            let (x3, y3, z3) = (x2 * x1, y2 * y1, z2 * z1);
            let (x4, y4, z4) = (x2 * x2, y2 * y2, z2 * z2);

            let monomials = [
                x4, y4, z4,
                x3 * y1, x3 * z1, y3 * x1, y3 * z1, z3 * x1, z3 * y1,
                x2 * y2, x2 * z2, y2 * z2,
                x2 * y1 * z1, x1 * y2 * z1, x1 * y1 * z2,
                x3, y3, z3,
                x2 * y1, x2 * z1, y2 * x1, y2 * z1, z2 * x1, z2 * y1,
                x1 * y1 * z1,
                x2, y2, z2,
                x1 * y1, x1 * z1, y1 * z1,
                x1, y1, z1,
                1.0,
            ];
            let delta: f64 = coeffs.iter().zip(monomials).map(|(c, m)| c * m).sum();

            assert_near!(
                delta,
                0.0,
                tolerance,
                "torus equation not satisfied at (u={u}, v={v}), delta={delta}"
            );

            v = if v_ind == nb_pts_max { 2.0 * PI } else { v + step };
        }
        u = if u_ind == nb_pts_max { 2.0 * PI } else { u + step };
    }
}
