//! The extremum nearest a guess, checked strictly against an independent brute-force list of
//! every candidate in the range.

use std::f64::consts::PI;

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Ellipse, Hyperbola, Line, OffsetCurve, Parabola};
use crate::geometry::ext_pc::tests::{brute_extrema, brute_f, Lcg};
use crate::math::Direction;

fn cubic() -> Curve {
    // Knots 0..4 of multiplicities 4 2 3 1 4: C1 at 1, a C0 corner at 2, C2 at 3.
    Curve::BSpline(BSplineCurve::new(
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(2.0, 2.5, 0.5),
            Point::new(3.0, 1.0, 1.0),
            Point::new(4.0, 0.0, 1.0),
            Point::new(4.5, 1.5, 0.5),
            Point::new(5.0, 3.0, 0.0),
            Point::new(6.0, 2.0, -1.0),
            Point::new(7.0, 0.0, 0.0),
            Point::new(8.0, 1.0, 1.0),
        ],
        vec![0.0, 1.0, 2.0, 3.0, 4.0],
        vec![4, 2, 3, 1, 4],
        3,
    ))
}

fn tilted() -> Frame {
    Frame::new(Point::new(0.5, -0.5, 0.25), Direction::new(0.2, -0.1, 1.0), Direction::new(1.0, 0.4, 0.0))
}

/// The curves of the test and their ranges.
fn curves() -> Vec<(&'static str, GeomAdaptor)> {
    let bez = Curve::Bezier(BezierCurve::new(vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 3.0, 0.0),
        Point::new(2.0, -3.0, 1.0),
        Point::new(3.0, 3.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
    ]));
    let c1 = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.5, 0.0), Point::new(2.0, -0.5, 0.0), Point::new(3.0, 1.0, 0.0), Point::new(4.0, 0.0, 0.0), Point::new(5.0, 1.0, 0.0), Point::new(6.0, -0.5, 0.0)],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 1, 2, 4],
        3,
    ));
    let offset = Curve::Offset(OffsetCurve::new(c1, 0.8, Direction::Z));
    let circle = Curve::Circle(Circle::new(tilted(), 2.0));
    vec![
        ("bspline", GeomAdaptor::new(&cubic(), 0.0, 4.0)),
        ("bspline_sub", GeomAdaptor::new(&cubic(), 0.5, 3.5)),
        ("bezier", GeomAdaptor::new(&bez, 0.0, 1.0)),
        ("circle", GeomAdaptor::new(&circle, 0.0, 2.0 * PI)),
        ("circle_arc", GeomAdaptor::new(&circle, 1.0, 4.0)),
        ("line", GeomAdaptor::new(&Curve::Line(Line::new(Point::new(1.0, 1.0, 0.0), Direction::new(1.0, 0.5, 0.2))), -2.0, 3.0)),
        ("ellipse", GeomAdaptor::new(&Curve::Ellipse(Ellipse::new(tilted(), 3.0, 1.0)), 0.0, 2.0 * PI)),
        ("hyperbola", GeomAdaptor::new(&Curve::Hyperbola(Hyperbola::new(tilted(), 1.0, 0.6)), -2.0, 2.0)),
        ("parabola", GeomAdaptor::new(&Curve::Parabola(Parabola::new(tilted(), 0.5)), -3.0, 3.0)),
        ("offset", GeomAdaptor::new(&offset, 0.0, 3.0)),
    ]
}

/// The candidates of the contract, by brute force: `(u, is_min)`, ascending.
fn candidates(name: &str, a: &GeomAdaptor, p: Point, lo: f64, hi: f64) -> Vec<(f64, bool)> {
    if name.starts_with("circle") {
        // Every point equidistant when the point projects onto the centre.
        let f = tilted();
        let w = p - f.origin();
        let (x, y) = (w.dot(Vector::from(f.x())), w.dot(Vector::from(f.y())));
        if x.hypot(y) < 1e-7 {
            return Vec::new();
        }
    }
    brute_extrema(a, p, lo, hi, 4000, true).into_iter().map(|(u, m, _)| (u, m)).collect()
}

/// The parameters (and kinds) acceptable as the nearest: the nearest one, and the other of a
/// tie whose gaps differ by less than `1e-9`.
fn expected(list: &[(f64, bool)], u0: f64) -> Vec<(f64, bool)> {
    let best = list.iter().map(|c| (c.0 - u0).abs()).fold(f64::INFINITY, f64::min);
    let mut near: Vec<(f64, bool)> = list.iter().copied().filter(|c| (c.0 - u0).abs() - best < 1e-9).collect();
    near.sort_by(|a, b| b.0.total_cmp(&a.0));
    near
}

#[test]
fn locate_ext_pc_is_the_nearest_extremum() {
    let mut rng = Lcg(2026);
    let mut cases = 0;
    let mut found_some = 0;
    for (name, a) in curves() {
        let (f, l) = (a.first_parameter(), a.last_parameter());
        let w = l - f;
        for k in 0..40 {
            let p = Point::new(-2.0 + 9.0 * rng.next(), -3.0 + 6.0 * rng.next(), -2.0 + 4.0 * rng.next());
            let (mut lo, mut hi) = (f, l);
            if name == "bspline_sub" || (k % 5 == 4 && !name.starts_with("circle") && name != "ellipse") {
                let (x, y) = (f + w * rng.next(), f + w * rng.next());
                (lo, hi) = (x.min(y), x.max(y));
                if hi - lo < 1e-2 {
                    continue;
                }
            }
            let list = candidates(name, &a, p, lo, hi);
            let u0 = match k % 8 {
                0 => lo,
                1 => hi,
                2 if !list.is_empty() => list[list.len() / 2].0,
                3 => hi + 0.3 * w,
                4 => lo - 0.1 * w,
                _ => lo + (hi - lo) * rng.next(),
            };
            check_case(name, &a, p, u0, lo, hi, &list);
            cases += 1;
            found_some += usize::from(!list.is_empty());
        }
    }
    // The degenerate circles: the point on the axis.
    for name in ["circle", "circle_arc"] {
        let a = curves().into_iter().find(|c| c.0 == name).unwrap().1;
        let axis = Vector::from(tilted().z());
        let p = tilted().origin() + axis * 3.0000001;
        assert!(locate_ext_pc(p, &a, 1.5, a.first_parameter(), a.last_parameter(), 1e-9).is_none());
    }
    assert!(cases >= 380, "{cases}");
    assert!(found_some > cases / 2);
}

fn check_case(name: &str, a: &GeomAdaptor, p: Point, u0: f64, lo: f64, hi: f64, list: &[(f64, bool)]) {
    let label = format!("{name}: P {p:?}, u0 {u0}, [{lo}, {hi}]");
    let found = locate_ext_pc(p, a, u0, lo, hi, 1e-9);
    assert_eq!(found.is_some(), !list.is_empty(), "{label}: {found:?} vs {list:?}");
    let Some(found) = found else { return };
    let ok = expected(list, u0);
    let matched = ok.iter().find(|e| (e.0 - found.parameter).abs() <= 1e-9);
    let Some(&(u, is_min)) = matched else { panic!("{label}: {found:?}, expected {ok:?} of {list:?}") };
    assert_eq!(found.is_min, is_min, "{label}: kind at {u}");
    let q = a.value(found.parameter);
    assert!(q.distance(found.point) <= 1e-12 * (1.0 + q.distance(Point::ORIGIN)), "{label}");
    assert!((found.square_distance - q.distance_squared(p)).abs() <= 1e-12 * (1.0 + found.square_distance), "{label}");
    if brute_f(a, p, u).abs() < 1e-7 {
        assert!(brute_f(a, p, found.parameter).abs() < 1e-7, "{label}");
    }
}

#[test]
fn closed_forms() {
    // Circle: the projected direction is the minimum, its opposite the maximum.
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 2.0));
    let a = GeomAdaptor::new(&circle, 0.0, 2.0 * PI);
    let p = Point::new(1.0, 1.0, 0.5);
    let near_min = locate_ext_pc(p, &a, 1.0, 0.0, 2.0 * PI, 1e-9).unwrap();
    assert!((near_min.parameter - PI / 4.0).abs() < 1e-14 && near_min.is_min);
    let near_max = locate_ext_pc(p, &a, 4.0, 0.0, 2.0 * PI, 1e-9).unwrap();
    assert!((near_max.parameter - 5.0 * PI / 4.0).abs() < 1e-14 && !near_max.is_min);
    // Line: the foot, or nothing outside the range.
    let line = Curve::Line(Line::new(Point::ORIGIN, Direction::X));
    let l = GeomAdaptor::new(&line, 0.0, 1.0);
    assert!((locate_ext_pc(Point::new(0.25, 1.0, 0.0), &l, 0.9, 0.0, 1.0, 1e-9).unwrap().parameter - 0.25).abs() < 1e-15);
    assert!(locate_ext_pc(Point::new(2.0, 1.0, 0.0), &l, 0.9, 0.0, 1.0, 1e-9).is_none());
    // y = x³ − x as a Bézier, the point on the y axis: three extrema, the nearest chosen.
    let bez = Curve::Bezier(BezierCurve::new(vec![
        Point::new(-1.5, -1.875, 0.0),
        Point::new(-0.5, 3.875, 0.0),
        Point::new(0.5, -3.875, 0.0),
        Point::new(1.5, 1.875, 0.0),
    ]));
    let b = GeomAdaptor::new(&bez, 0.0, 1.0);
    let p = Point::new(-0.577, -0.5, 0.0);
    let list = candidates("bezier", &b, p, 0.0, 1.0);
    assert!(list.len() >= 3, "{list:?}");
    for (u, _) in &list {
        assert!((locate_ext_pc(p, &b, *u + 1e-3, 0.0, 1.0, 1e-9).unwrap().parameter - u).abs() < 1e-9);
        // A guess exactly at an extremum returns it.
        assert!((locate_ext_pc(p, &b, *u, 0.0, 1.0, 1e-9).unwrap().parameter - u).abs() < 1e-12);
    }
    // A degree-1 polyline with the point outside a vertex: the corner.
    let poly = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0)],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        1,
    ));
    let a = GeomAdaptor::new(&poly, 0.0, 2.0);
    let e = locate_ext_pc(Point::new(3.0, -1.0, 0.0), &a, 0.2, 0.0, 2.0, 1e-9).unwrap();
    assert_eq!(e.parameter, 1.0);
    assert!(e.is_min);
}
