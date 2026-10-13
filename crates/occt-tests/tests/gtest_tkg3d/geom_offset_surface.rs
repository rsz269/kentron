//! Geom_OffsetSurface_Test.cxx: the equivalent non-offset surface and derivative
//! consistency between an offset surface and its equivalent.
//!
//! The `ThreeWay` tests also compare OCCT's adaptor evaluation; kentron has a single
//! evaluation path, so they are merged into the two-way checks. Copy tests are not ported.

use std::f64::consts::{FRAC_PI_3, FRAC_PI_4, FRAC_PI_6, PI};

use kentron::geometry::surface::{Cone, Cylinder, OffsetSurface, Plane, Sphere, Surface, Torus, TrimmedSurface};
use kentron::geometry::ParametricSurface;
use kentron::math::{Direction, Frame};
use kentron::tolerance;

use crate::{assert_point_near, p};

const TOL: f64 = tolerance::LINEAR;

fn cylinder(radius: f64) -> Surface {
    Surface::Cylinder(Cylinder::new(Frame::WORLD, radius))
}

fn torus(major: f64, minor: f64) -> Surface {
    Surface::Torus(Torus::new(Frame::WORLD, major, minor))
}

fn equivalent(basis: Surface, offset: f64) -> Surface {
    OffsetSurface::new(basis, offset).equivalent_surface().expect("no equivalent surface")
}

#[test]
fn plane_becomes_translated_plane() {
    let Surface::Plane(plane) = equivalent(Surface::Plane(Plane::new(p(1.0, 2.0, 3.0), Direction::Z)), 5.0) else {
        panic!("not a plane")
    };
    assert_point_near(plane.origin(), p(1.0, 2.0, 8.0), TOL);
}

#[test]
fn cylinder_radius_changes() {
    for (radius, offset, expected) in [(10.0, 3.0, 13.0), (10.0, -3.0, 7.0), (5.0, -8.0, 3.0)] {
        let Surface::Cylinder(c) = equivalent(cylinder(radius), offset) else { panic!("not a cylinder") };
        assert!((c.radius() - expected).abs() <= TOL);
    }
}

#[test]
fn sphere_radius_changes() {
    let Surface::Sphere(s) = equivalent(Surface::Sphere(Sphere::new(Frame::from_z(p(1.0, 2.0, 3.0), Direction::Z), 10.0)), 5.0)
    else {
        panic!("not a sphere")
    };
    assert!((s.radius() - 15.0).abs() <= TOL);
    assert_point_near(s.frame().origin(), p(1.0, 2.0, 3.0), TOL);
    let Surface::Sphere(s) = equivalent(Surface::Sphere(Sphere::new(Frame::WORLD, 5.0)), -8.0) else { panic!("not a sphere") };
    assert!((s.radius() - 3.0).abs() <= TOL);
}

#[test]
fn cone_reference_radius_changes() {
    let angle = PI / 6.0;
    let Surface::Cone(c) = equivalent(Surface::Cone(Cone::new(Frame::WORLD, angle, 10.0)), 2.0) else { panic!("not a cone") };
    assert!((c.semi_angle() - angle).abs() <= TOL);
    assert!((c.radius() - (10.0 + 2.0 * angle.cos())).abs() <= TOL);
}

#[test]
fn torus_minor_radius_changes() {
    for (minor, offset, expected) in [(5.0, 2.0, 7.0), (5.0, -2.0, 3.0), (3.0, -5.0, 2.0)] {
        let Surface::Torus(t) = equivalent(torus(20.0, minor), offset) else { panic!("not a torus") };
        assert!((t.major_radius() - 20.0).abs() <= TOL);
        assert!((t.minor_radius() - expected).abs() <= TOL);
    }
}

#[test]
fn zero_offset_gives_the_basis() {
    assert_eq!(equivalent(cylinder(10.0), 0.0), cylinder(10.0));
}

#[test]
fn trimmed_plane_stays_trimmed() {
    let trimmed = Surface::Trimmed(TrimmedSurface::new(Surface::Plane(Plane::XY), 0.0, 10.0, 0.0, 5.0));
    let Surface::Trimmed(t) = equivalent(trimmed, 3.0) else { panic!("not trimmed") };
    let (u1, u2, v1, v2) = t.trim();
    assert!(u1.abs() <= TOL && (u2 - 10.0).abs() <= TOL && v1.abs() <= TOL && (v2 - 5.0).abs() <= TOL);
}

#[test]
fn degenerate_cylinder_has_no_equivalent() {
    assert_eq!(OffsetSurface::new(cylinder(1.0), -1.0 + 0.5 * TOL).equivalent_surface(), None);
}

/// Compares the offset surface and its equivalent up to order `order` on a grid.
fn assert_consistent(basis: Surface, offset: f64, us: &[f64], vs: &[f64], order: u32) {
    let off = OffsetSurface::new(basis, offset);
    let eq = off.equivalent_surface().expect("no equivalent surface");
    let tol = 1.0e-10;
    for &u in us {
        for &v in vs {
            assert_point_near(off.point(u, v), eq.point(u, v), tol);
            let ((_, a1, a2), (_, b1, b2)) = (off.d1(u, v), eq.d1(u, v));
            assert!((a1 - b1).length() < tol && (a2 - b2).length() < tol, "d1 at ({u}, {v})");
            if order >= 2 {
                let (a, b) = (off.d2(u, v), eq.d2(u, v));
                for (x, y) in [(a.duu, b.duu), (a.dvv, b.dvv), (a.duv, b.duv)] {
                    assert!((x - y).length() < tol, "d2 at ({u}, {v})");
                }
            }
            if order >= 3 {
                let (a, b) = (off.d3(u, v), eq.d3(u, v));
                for (x, y) in [(a.duuu, b.duuu), (a.dvvv, b.dvvv), (a.duuv, b.duuv), (a.duvv, b.duvv)] {
                    assert!((x - y).length() < tol, "d3 at ({u}, {v})");
                }
            }
        }
    }
}

fn steps(from: f64, to: f64, step: f64) -> Vec<f64> {
    let mut v = Vec::new();
    let mut x = from;
    while x <= to + 1e-12 {
        v.push(x);
        x += step;
    }
    v
}

fn angles(step: f64) -> Vec<f64> {
    steps(0.0, 2.0 * PI - 1e-9, step)
}

#[test]
fn derivatives_plane() {
    let s = Surface::Plane(Plane::new(p(1.0, 2.0, 3.0), Direction::Z));
    assert_consistent(s, 5.0, &steps(-5.0, 5.0, 2.5), &steps(-5.0, 5.0, 2.5), 2);
}

#[test]
fn derivatives_cylinder() {
    assert_consistent(cylinder(10.0), 3.0, &angles(FRAC_PI_4), &steps(-5.0, 5.0, 2.5), 3);
    assert_consistent(cylinder(10.0), -3.0, &angles(FRAC_PI_3), &steps(-3.0, 3.0, 2.0), 2);
    assert_consistent(cylinder(5.0), -8.0, &angles(FRAC_PI_3), &steps(-3.0, 3.0, 2.0), 1);
}

#[test]
fn derivatives_sphere() {
    let s = Surface::Sphere(Sphere::new(Frame::from_z(p(1.0, 2.0, 3.0), Direction::Z), 10.0));
    assert_consistent(s, 5.0, &angles(FRAC_PI_4), &steps(-PI / 3.0, PI / 3.0, FRAC_PI_6), 3);
}

#[test]
fn derivatives_cone() {
    let s = Surface::Cone(Cone::new(Frame::WORLD, PI / 6.0, 10.0));
    assert_consistent(s, 2.0, &angles(FRAC_PI_4), &steps(-5.0, 5.0, 2.5), 2);
}

#[test]
fn derivatives_torus() {
    assert_consistent(torus(20.0, 5.0), 2.0, &angles(FRAC_PI_4), &angles(FRAC_PI_4), 3);
}

#[test]
fn dn_cylinder() {
    let off = OffsetSurface::new(cylinder(10.0), 3.0);
    let eq = off.equivalent_surface().unwrap();
    for nu in 0..=3 {
        for nv in 0..=3 {
            if nu + nv == 0 {
                continue;
            }
            let (a, b) = (off.dn(FRAC_PI_3, 2.5, nu, nv), eq.dn(FRAC_PI_3, 2.5, nu, nv));
            assert!((a - b).length() < 1.0e-10, "dn({nu}, {nv})");
        }
    }
}
