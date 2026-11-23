//! Point–curve extrema against closed forms and against an independent brute-force search
//! (dense sampling of `F`, bisection, one-sided values at breakpoints).

use std::f64::consts::PI;

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Curve, OffsetCurve};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Frame};

/// A view hiding the kind of the curve it wraps (so that only the general path applies).
pub(crate) struct Plain<'a>(pub &'a dyn CurveAdaptor);

impl CurveAdaptor for Plain<'_> {
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter()
    }
    fn last_parameter(&self) -> f64 {
        self.0.last_parameter()
    }
    fn intervals(&self, c: Continuity) -> Vec<f64> {
        self.0.intervals(c)
    }
    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }
    fn period(&self) -> f64 {
        self.0.period()
    }
    fn value(&self, u: f64) -> Point {
        self.0.value(u)
    }
    fn d1(&self, u: f64) -> (Point, Vector) {
        self.0.d1(u)
    }
    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        self.0.d2(u)
    }
    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
        self.0.d3(u)
    }
    fn dn(&self, u: f64, n: u32) -> Vector {
        assert!(n <= 3, "dn of order {n}");
        self.0.dn(u, n)
    }
    fn resolution(&self, r3d: f64) -> f64 {
        self.0.resolution(r3d)
    }
}

/// The wave `(t, ε sin kt, 0)`, a view without a kind.
struct Wave {
    eps: f64,
    k: f64,
}

impl CurveAdaptor for Wave {
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        10.0
    }
    fn intervals(&self, _: Continuity) -> Vec<f64> {
        vec![0.0, 10.0]
    }
    fn is_periodic(&self) -> bool {
        false
    }
    fn period(&self) -> f64 {
        0.0
    }
    fn value(&self, t: f64) -> Point {
        Point::new(t, self.eps * (self.k * t).sin(), 0.0)
    }
    fn d1(&self, t: f64) -> (Point, Vector) {
        (self.value(t), Vector::new(1.0, self.eps * self.k * (self.k * t).cos(), 0.0))
    }
    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        let (p, d) = self.d1(t);
        (p, d, Vector::new(0.0, -self.eps * self.k * self.k * (self.k * t).sin(), 0.0))
    }
    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        let (p, d, dd) = self.d2(t);
        (p, d, dd, Vector::new(0.0, -self.eps * self.k.powi(3) * (self.k * t).cos(), 0.0))
    }
    fn dn(&self, t: f64, n: u32) -> Vector {
        match n {
            1 => self.d1(t).1,
            2 => self.d2(t).2,
            3 => self.d3(t).3,
            _ => panic!("dn of order {n}"),
        }
    }
    fn resolution(&self, r3d: f64) -> f64 {
        r3d / (1.0 + self.eps * self.k)
    }
}

/// `F = (C − P)·T` evaluated independently (the tangent from `C″` where `C′` vanishes).
pub(crate) fn brute_f(curve: &dyn CurveAdaptor, p: Point, u: f64) -> f64 {
    let (c, d1, d2) = curve.d2(u);
    let w = c - p;
    let s = d1.length();
    if s > 0.0 {
        w.dot(d1) / s
    } else if d2.length() > 0.0 {
        w.dot(d2) / d2.length()
    } else {
        0.0
    }
}

/// Every extremum of the distance on `[lo, hi]` by brute force: `per_piece` samples of each
/// C2 piece, sign changes bisected to `1e-13`, corners from `F` at `b ∓ 1e-12`, and (with
/// `stationary_ends`) range ends where `|F| < 1e-7`. `(u, is_min, at_breakpoint)`, ascending.
pub(crate) fn brute_extrema(curve: &dyn CurveAdaptor, p: Point, lo: f64, hi: f64, per_piece: usize, stationary_ends: bool) -> Vec<(f64, bool, bool)> {
    let mut breaks: Vec<f64> = curve.intervals(Continuity::C2).into_iter().filter(|&b| b > lo + 1e-12 && b < hi - 1e-12).collect();
    if curve.is_periodic() {
        breaks.clear();
    }
    let f = |u: f64| brute_f(curve, p, u);
    let mut out = Vec::new();
    let mut stops = vec![lo];
    stops.extend_from_slice(&breaks);
    stops.push(hi);
    for w in stops.windows(2) {
        let a = if w[0] == lo { lo } else { w[0] + 1e-12 };
        let b = if w[1] == hi { hi } else { w[1] - 1e-12 };
        let mut prev = (a, f(a));
        for i in 1..=per_piece {
            let u = if i == per_piece { b } else { a + (b - a) * i as f64 / per_piece as f64 };
            let fu = f(u);
            if fu == 0.0 && i < per_piece {
                // (An exact zero: compared across with the next sample.)
                continue;
            }
            if (prev.1 < 0.0 && fu > 0.0) || (prev.1 > 0.0 && fu < 0.0) {
                let (mut x0, mut f0, mut x1) = (prev.0, prev.1, u);
                while x1 - x0 > 1e-13 {
                    let m = 0.5 * (x0 + x1);
                    let fm = f(m);
                    if fm == 0.0 {
                        x0 = m;
                        x1 = m;
                        break;
                    }
                    if (fm < 0.0) == (f0 < 0.0) {
                        x0 = m;
                        f0 = fm;
                    } else {
                        x1 = m;
                    }
                }
                out.push((0.5 * (x0 + x1), prev.1 < 0.0, false));
            }
            prev = (u, fu);
        }
    }
    for &b in &breaks {
        let (l, r) = (f(b - 1e-12), f(b + 1e-12));
        if l < 0.0 && r > 0.0 {
            out.push((b, true, true));
        } else if l > 0.0 && r < 0.0 {
            out.push((b, false, true));
        }
    }
    if stationary_ends {
        for (end, inward) in [(lo, 1.0), (hi, -1.0)] {
            let at = if breaks.is_empty() { end } else { end };
            let fe = f(at);
            if fe.abs() < 1e-7 {
                out.push((end, inward * fe >= 0.0, false));
            }
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.dedup_by(|b, a| (b.0 - a.0).abs() < 1e-9);
    out
}

fn whole(curve: &Curve) -> GeomAdaptor {
    GeomAdaptor::new(curve, curve.first_parameter(), curve.last_parameter())
}

/// A deterministic sequence in `[0, 1)`.
pub(crate) struct Lcg(pub u64);

impl Lcg {
    pub(crate) fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Checks a result against the brute-force list (parameters at `1e-9`, kinds, records).
fn check(curve: &dyn CurveAdaptor, p: Point, lo: f64, hi: f64, found: &[ExtPCPoint], label: &str) {
    let expected = brute_extrema(curve, p, lo, hi, 4000, false);
    let got: Vec<(f64, bool)> = found.iter().map(|e| (e.parameter, e.is_min)).collect();
    assert_eq!(got.len(), expected.len(), "{label}: {got:?} vs {expected:?}");
    for (g, e) in got.iter().zip(&expected) {
        assert!((g.0 - e.0).abs() < 1e-9, "{label}: {got:?} vs {expected:?}");
        assert_eq!(g.1, e.1, "{label}: kind at {}", g.0);
    }
    for e in found {
        let q = curve.value(e.parameter);
        assert!(q.distance(e.point) <= 1e-12 * (1.0 + q.distance(Point::ORIGIN)));
        assert!((e.square_distance - q.distance_squared(p)).abs() <= 1e-12 * (1.0 + e.square_distance));
    }
}

#[test]
fn circle_off_centre_and_at_the_centre() {
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 2.0));
    let a = GeomAdaptor::new(&circle, 0.0, 2.0 * PI);
    let plain = Plain(&a);
    let p = Point::new(0.5, 0.3, 1.0);
    let found = ext_pc_other(&plain, p, 0.0, 2.0 * PI, 1e-10).unwrap();
    let t = 0.3f64.atan2(0.5);
    assert_eq!(found.len(), 2);
    assert!((found[0].parameter - t).abs() < 1e-12 && found[0].is_min);
    assert!((found[1].parameter - (t + PI)).abs() < 1e-12 && !found[1].is_min);
    // At the centre (above it): equidistant everywhere.
    assert!(ext_pc_other(&plain, Point::new(0.0, 0.0, 3.0), 0.0, 2.0 * PI, 1e-10).is_none());
    // Over a range shifted by a period: parameters in [uinf, uinf + 2π).
    let found = ext_pc_other(&plain, p, 2.0 * PI + 1.0, 4.0 * PI + 1.0, 1e-10).unwrap();
    assert_eq!(found.len(), 2);
    assert!(found.iter().all(|e| e.parameter >= 2.0 * PI + 1.0 && e.parameter < 4.0 * PI + 1.0));
}

#[test]
fn through_the_point() {
    let bez = Curve::Bezier(BezierCurve::new(vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 2.0, 0.0), Point::new(2.0, -1.0, 1.0), Point::new(3.0, 0.5, 0.0)]));
    let a = whole(&bez);
    let plain = Plain(&a);
    // Interior: a minimum at zero distance.
    let p = a.value(0.37);
    let found = ext_pc_other(&plain, p, 0.0, 1.0, 1e-10).unwrap();
    let hit = found.iter().find(|e| (e.parameter - 0.37).abs() < 1e-9).expect("the point itself");
    assert!(hit.is_min && hit.square_distance < 1e-24);
    check(&plain, p, 0.0, 1.0, &found, "interior");
    // At the range start: reported first.
    let found = ext_pc_other(&plain, Point::ORIGIN, 0.0, 1.0, 1e-10).unwrap();
    assert_eq!(found[0].parameter, 0.0);
    assert!(found[0].is_min && found[0].square_distance == 0.0);
    // A tiny range at the point.
    let found = ext_pc_other(&plain, Point::ORIGIN, 0.0, 1e-15, 1e-10).unwrap();
    assert_eq!(found.len(), 1);
    assert!(ext_pc_other(&plain, Point::new(5.0, 5.0, 5.0), 0.3, 0.3, 1e-10).is_none());
}

#[test]
fn corners_of_a_polyline() {
    let poly = Curve::BSpline(BSplineCurve::new(
        vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0), Point::new(4.0, 2.0, 0.0)],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![2, 1, 1, 2],
        1,
    ));
    let a = whole(&poly);
    let plain = Plain(&a);
    // Outside the first corner: it is the nearest point, a corner minimum exactly at 1.
    let p = Point::new(3.0, -1.0, 0.0);
    let found = ext_pc_other(&plain, p, 0.0, 3.0, 1e-10).unwrap();
    let corner = found.iter().find(|e| e.parameter == 1.0).expect("the corner");
    assert!(corner.is_min);
    check(&plain, p, 0.0, 3.0, &found, "polyline");
    // Inside the second corner it is a maximum.
    let p = Point::new(1.0, 3.0, 0.0);
    let found = ext_pc_other(&plain, p, 0.0, 3.0, 1e-10).unwrap();
    check(&plain, p, 0.0, 3.0, &found, "polyline inside");
}

#[test]
fn a_fast_wiggle() {
    let w = Wave { eps: 0.2, k: 25.0 };
    let p = Point::new(5.0, -0.2, 0.0);
    let found = ext_pc_other(&w, p, 0.0, 10.0, 1e-10).unwrap();
    assert!(found.len() > 15, "{}", found.len());
    check(&w, p, 0.0, 10.0, &found, "wave");
    for pair in found.windows(2) {
        assert_ne!(pair[0].is_min, pair[1].is_min);
    }
}

/// Random points against B-splines with knots of every continuity, Béziers and offset
/// curves with cusps, over whole ranges and sub-ranges: the same extrema as brute force.
#[test]
fn agrees_with_brute_force() {
    let cubic = Curve::BSpline(BSplineCurve::new(
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
    ));
    let quartic = Curve::Bezier(BezierCurve::new(vec![
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
    let offset = Curve::Offset(OffsetCurve::new(c1.clone(), 0.8, Direction::Z));
    let curves = [("cubic", cubic), ("quartic", quartic), ("c1", c1), ("offset", offset)];
    let mut rng = Lcg(7);
    let mut cases = 0;
    for (name, curve) in &curves {
        let a = whole(curve);
        let plain = Plain(&a);
        let (f, l) = (a.first_parameter(), a.last_parameter());
        for _ in 0..25 {
            let p = Point::new(-1.0 + 10.0 * rng.next(), -2.0 + 5.0 * rng.next(), -1.0 + 2.0 * rng.next() * if *name == "offset" { 0.0 } else { 1.0 });
            let (mut lo, mut hi) = (f, l);
            if rng.next() < 0.4 {
                let (x, y) = (f + (l - f) * rng.next(), f + (l - f) * rng.next());
                (lo, hi) = (x.min(y), x.max(y));
            }
            if hi - lo < 1e-3 {
                continue;
            }
            let found = ext_pc_other(&plain, p, lo, hi, 1e-10).unwrap_or_default();
            check(&plain, p, lo, hi, &found, &format!("{name} {p:?} [{lo}, {hi}]"));
            cases += 1;
        }
    }
    assert!(cases > 80);
}

/// A range of a few hundred units in the last place, sampled to its bottom: the sampling
/// stops where midpoints stop being new (it once recursed without end there).
#[test]
fn a_range_of_a_few_ulps() {
    struct Fine(Wave);
    impl CurveAdaptor for Fine {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            10.0
        }
        fn intervals(&self, c: Continuity) -> Vec<f64> {
            self.0.intervals(c)
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn value(&self, u: f64) -> Point {
            self.0.value(u)
        }
        fn d1(&self, u: f64) -> (Point, Vector) {
            self.0.d1(u)
        }
        fn d2(&self, u: f64) -> (Point, Vector, Vector) {
            self.0.d2(u)
        }
        fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
            self.0.d3(u)
        }
        fn dn(&self, u: f64, n: u32) -> Vector {
            self.0.dn(u, n)
        }
        fn resolution(&self, _: f64) -> f64 {
            1e-300
        }
    }
    let c = Fine(Wave { eps: 0.2, k: 25.0 });
    let lo = 3.0;
    let hi = lo + 300.0 * f64::EPSILON * lo;
    // (The point beside the curve: the direction to it turns fast across the range.)
    let p = c.value(0.5 * (lo + hi)) + Vector::new(0.0, 0.0, 1e-14);
    // (No stretch tolerance: a list, whatever it holds.)
    assert!(ext_pc_other(&c, p, lo, hi, 0.0).is_some());
}
