//! Arc lengths and their inverse (and the arc-length reparametrisation built on them)
//! against closed forms and against references computed to a much finer tolerance.

use std::f64::consts::PI;

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Ellipse, Line};
use crate::geometry::curvlin::CurvlinFunc;
use crate::math::{Direction, Frame, Point, Vector};

/// The helix `(r cos t, r sin t, h t)`, a view without a kind.
struct Helix {
    r: f64,
    h: f64,
    first: f64,
    last: f64,
}

impl CurveAdaptor for Helix {
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn intervals(&self, _: Continuity) -> Vec<f64> {
        vec![self.first, self.last]
    }
    fn is_periodic(&self) -> bool {
        false
    }
    fn period(&self) -> f64 {
        0.0
    }
    fn value(&self, t: f64) -> Point {
        Point::new(self.r * t.cos(), self.r * t.sin(), self.h * t)
    }
    fn d1(&self, t: f64) -> (Point, Vector) {
        (self.value(t), Vector::new(-self.r * t.sin(), self.r * t.cos(), self.h))
    }
    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        let (p, d1) = self.d1(t);
        (p, d1, Vector::new(-self.r * t.cos(), -self.r * t.sin(), 0.0))
    }
    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        let (p, d1, d2) = self.d2(t);
        (p, d1, d2, Vector::new(self.r * t.sin(), -self.r * t.cos(), 0.0))
    }
    fn dn(&self, t: f64, n: u32) -> Vector {
        match n {
            1 => self.d1(t).1,
            2 => self.d2(t).2,
            _ => self.d3(t).3,
        }
    }
    fn resolution(&self, r3d: f64) -> f64 {
        r3d / (self.r * self.r + self.h * self.h).sqrt()
    }
}

fn whole(curve: &Curve) -> GeomAdaptor {
    use crate::geometry::ParametricCurve;
    GeomAdaptor::new(curve, curve.first_parameter(), curve.last_parameter())
}

/// A cubic B-spline over `[0, 4]` with a C1 knot at 1, a C0 corner at 2 and a C2 knot at 3.
fn cornered_cubic() -> Curve {
    let poles = vec![
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
    ];
    Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![4, 2, 3, 1, 4], 3))
}

/// A quintic B-spline with single knots.
fn quintic() -> Curve {
    let poles: Vec<Point> = (0..9)
        .map(|i| {
            let t = i as f64;
            Point::new(t, (t * 1.3).sin() * 2.0, (t * 0.7).cos())
        })
        .collect();
    Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![6, 1, 1, 1, 6], 5))
}

/// A periodic cubic B-spline over `[0, 6]`.
fn periodic() -> Curve {
    let poles: Vec<Point> = (0..6)
        .map(|i| {
            let t = i as f64 * PI / 3.0;
            Point::new(2.0 * t.cos(), t.sin(), 0.4 * (2.0 * t).sin())
        })
        .collect();
    Curve::BSpline(BSplineCurve::new_periodic(poles, (0..=6).map(|i| i as f64).collect(), vec![1; 7], 3))
}

#[test]
fn kronrod_rules_are_exact_for_polynomials() {
    for degree in 0..=22 {
        let f = |x: f64| x.powi(degree);
        let (a, b) = (-0.3f64, 1.7f64);
        let exact = (b.powi(degree + 1) - a.powi(degree + 1)) / (degree + 1) as f64;
        let (k, e) = kronrod(&f, a, b);
        assert!((k - exact).abs() <= 1e-14 * exact.abs().max(1.0), "degree {degree}: {k} vs {exact}");
        if degree <= 13 {
            assert!(e <= 1e-14 * exact.abs().max(1.0), "degree {degree}: Gauss off by {e}");
        }
    }
}

#[test]
fn constant_speed_curves_are_exact() {
    let line = Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::new(1.0, 2.0, 2.0)));
    let a = GeomAdaptor::new(&line, -2.0, 5.0);
    assert!((length(&a, None) - 7.0).abs() < 1e-14);
    assert!((length_between(&a, 4.0, -1.0, Some(1e-3)) - 5.0).abs() < 1e-14);
    assert!((AbscissaPoint::new(&a, 2.5, 0.5).unwrap().parameter().unwrap() - 3.0).abs() < 1e-14);
    assert!((AbscissaPoint::with_tolerance(1e-6, &a, -2.5, 0.5).parameter().unwrap() + 2.0).abs() < 1e-14);
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 3.0));
    let c = GeomAdaptor::new(&circle, 0.5, 2.0);
    assert!((length(&c, None) - 4.5).abs() < 1e-15);
    let u = AbscissaPoint::with_guess(&c, -1.5, 1.5, 1.0).parameter().unwrap();
    assert!((u - 1.0).abs() < 1e-15);
    let segment = Curve::BSpline(BSplineCurve::new(vec![Point::ORIGIN, Point::new(3.0, 4.0, 0.0)], vec![0.0, 2.0], vec![2, 2], 1));
    let s = whole(&segment);
    assert!((length(&s, None) - 5.0).abs() < 1e-14);
    assert!((AbscissaPoint::new(&s, 1.0, 0.0).unwrap().parameter().unwrap() - 0.4).abs() < 1e-15);
    let bez = Curve::Bezier(BezierCurve::new(vec![Point::ORIGIN, Point::new(0.0, 0.0, 2.0)]));
    let b = whole(&bez);
    assert!((length_between(&b, 0.25, 0.75, None) - 1.0).abs() < 1e-15);
}

#[test]
fn helix_length_and_abscissa() {
    let h = Helix { r: 2.0, h: 0.5, first: -1.0, last: 9.0 };
    let speed = (4.0f64 + 0.25).sqrt();
    assert!((length(&h, None) - 10.0 * speed).abs() < 1e-12);
    assert!((length_between(&h, 0.0, 3.0, Some(1e-6)) - 3.0 * speed).abs() < 1e-6);
    for a in [-3.0, -0.5, 0.25, 7.0, 15.0] {
        let u = AbscissaPoint::new(&h, a, 2.0).unwrap().parameter().unwrap();
        assert!((u - (2.0 + a / speed)).abs() < 1e-11, "{a}: {u}");
        let u = AbscissaPoint::with_guess_and_tolerance(&h, a, 2.0, 0.0, 1e-5).parameter().unwrap();
        assert!((u - (2.0 + a / speed)).abs() < 1e-5 / speed);
    }
}

/// The complete elliptic integral of the second kind, by the arithmetic-geometric mean.
fn elliptic_e(k: f64) -> f64 {
    let (mut a, mut g) = (1.0f64, (1.0 - k * k).sqrt());
    let mut sum = 0.5 * k * k;
    let mut power = 0.5;
    for _ in 0..40 {
        let c = 0.5 * (a - g);
        let (na, ng) = (0.5 * (a + g), (a * g).sqrt());
        a = na;
        g = ng;
        power *= 2.0;
        sum += power * c * c;
        if c.abs() < 1e-17 {
            break;
        }
    }
    PI / (2.0 * a) * (1.0 - sum)
}

#[test]
fn ellipse_quarter_and_parabola() {
    let (a, b) = (3.0, 1.2);
    let ellipse = Curve::Ellipse(Ellipse::new(Frame::WORLD, a, b));
    let e = GeomAdaptor::new(&ellipse, 0.0, PI / 2.0);
    let k = (1.0f64 - b * b / (a * a)).sqrt();
    let exact = a * elliptic_e(k);
    assert!((length(&e, None) - exact).abs() < 1e-9 * exact, "{} vs {exact}", length(&e, None));
    for tol in [1e-3, 1e-6, 1e-8, 1e-10] {
        assert!((length(&e, Some(tol)) - exact).abs() <= tol);
    }
    // y = x² from 0 to 1 as a quadratic Bézier.
    let parabola = Curve::Bezier(BezierCurve::new(vec![Point::ORIGIN, Point::new(0.5, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)]));
    let p = whole(&parabola);
    let exact = (2.0 * 5f64.sqrt() + 2f64.asinh()) / 4.0;
    assert!((length(&p, None) - exact).abs() < 1e-9 * exact);
    // The point at half the length: x with (2x√(1+4x²) + asinh 2x)/4 = exact/2.
    let u = AbscissaPoint::new(&p, exact / 2.0, 0.0).unwrap().parameter().unwrap();
    let s = |x: f64| (2.0 * x * (1.0 + 4.0 * x * x).sqrt() + (2.0 * x).asinh()) / 4.0;
    assert!((s(u) - exact / 2.0).abs() < 1e-10);
}

/// A reference length to `1e-13`.
fn reference(curve: &dyn CurveAdaptor, u1: f64, u2: f64) -> f64 {
    length_between(curve, u1, u2, Some(1e-13))
}

#[test]
fn tolerances_are_honoured_across_knots() {
    for curve in [cornered_cubic(), quintic(), periodic()] {
        let a = whole(&curve);
        let (f, l) = (a.first_parameter(), a.last_parameter());
        let exact = reference(&a, f, l);
        let none = length(&a, None);
        assert!((none - exact).abs() <= 1e-9 * exact, "{none} vs {exact}");
        for tol in [1e-3, 1e-5, 1e-7, 1e-9] {
            assert!((length(&a, Some(tol)) - exact).abs() <= tol);
            let part = length_between(&a, l - 0.3 * (l - f), f + 0.1, Some(tol));
            assert!((part - reference(&a, f + 0.1, l - 0.3 * (l - f))).abs() <= 2.0 * tol);
        }
    }
    // A degree-1 polyline: the sum of its sides.
    let poly = Curve::BSpline(BSplineCurve::new(
        vec![Point::ORIGIN, Point::new(3.0, 4.0, 0.0), Point::new(3.0, 4.0, 2.0), Point::new(0.0, 0.0, 2.0)],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![2, 1, 1, 2],
        1,
    ));
    let a = whole(&poly);
    assert!((length(&a, None) - 12.0).abs() < 1e-12);
    // Across the corners in both directions.
    assert!((AbscissaPoint::new(&a, 6.0, 0.0).unwrap().parameter().unwrap() - 1.5).abs() < 1e-12);
    assert!((AbscissaPoint::new(&a, -6.5, 2.5).unwrap().parameter().unwrap() - 0.6).abs() < 1e-12);
}

#[test]
fn abscissa_inverts_length() {
    for curve in [cornered_cubic(), quintic(), periodic()] {
        let a = whole(&curve);
        let (f, l) = (a.first_parameter(), a.last_parameter());
        let total = length(&a, None);
        for k in 0..=10 {
            let target = k as f64 * total / 10.0;
            let u = AbscissaPoint::new(&a, target, f).unwrap().parameter().unwrap();
            assert!((reference(&a, f, u) - target).abs() <= 1e-9 * total, "{k}");
            // Back from the end, and from the middle both ways.
            let v = AbscissaPoint::with_guess(&a, -target, l, f).parameter().unwrap();
            assert!((reference(&a, v, l) - target).abs() <= 1e-9 * total);
            let mid = 0.5 * (f + l);
            let half = 0.45 * target;
            for sign in [1.0, -1.0] {
                let w = AbscissaPoint::with_tolerance(1e-8, &a, sign * half, mid).parameter().unwrap();
                assert!((reference(&a, mid, w) - half).abs() <= 1e-8);
                assert!(sign * (w - mid) >= 0.0);
            }
        }
        // Zero abscissa: the start itself.
        assert_eq!(AbscissaPoint::new(&a, 0.0, f + 0.3).unwrap().parameter(), Some(f + 0.3));
        assert_eq!(AbscissaPoint::new(&a, 5e-8, f + 0.3).unwrap().parameter(), Some(f + 0.3));
    }
}

#[test]
fn errors_on_short_curves_and_bad_abscissae() {
    let line = Curve::Line(Line::new(Point::ORIGIN, Direction::X));
    let tiny = GeomAdaptor::new(&line, 0.0, 1e-9);
    assert!(AbscissaPoint::new(&tiny, 1e-10, 0.0).is_err());
    let curve = quintic();
    let a = whole(&curve);
    let mut f = CurvlinFunc::new(&a, 1e-8);
    assert!(f.get_u_parameter(1.5).is_err());
    assert!(f.get_u_parameter(-0.1).is_err());
    assert!(f.get_u_parameter(f64::NAN).is_err());
    assert!(f.trim(-0.1, 0.5, 1e-9).is_err());
    assert!(f.trim(0.2, 1.1, 1e-9).is_err());
}

#[test]
fn arc_length_round_trip() {
    for curve in [cornered_cubic(), quintic(), periodic()] {
        let a = whole(&curve);
        let f = CurvlinFunc::new(&a, 1e-10);
        assert!((f.get_s_parameter(a.last_parameter()) - 1.0).abs() < 1e-12);
        assert_eq!(f.get_u_parameter(0.0).unwrap(), a.first_parameter());
        assert!((f.get_u_parameter(1.0).unwrap() - a.last_parameter()).abs() < 1e-14);
        let total = reference(&a, a.first_parameter(), a.last_parameter());
        assert!((f.get_length() - total).abs() < 1e-10);
        for k in 0..=20 {
            let s = k as f64 / 20.0;
            let u = f.get_u_parameter(s).unwrap();
            assert!((f.get_s_parameter(u) - s).abs() < 1e-8);
            assert!((reference(&a, a.first_parameter(), u) - s * total).abs() < 1e-9);
        }
    }
}

#[test]
fn curvilinear_circle_and_trim() {
    let r = 2.0;
    let circle = Curve::Circle(Circle::new(Frame::WORLD, r));
    let a = GeomAdaptor::new(&circle, 0.5, 3.0);
    let mut f = CurvlinFunc::new(&a, 1e-8);
    let l = 2.5 * r;
    assert!((f.get_length() - l).abs() < 1e-14);
    assert_eq!(f.nb_intervals(Continuity::C2), 1);
    for s in [0.0, 0.1, 0.5, 0.77, 1.0] {
        let u = f.get_u_parameter(s).unwrap();
        assert!((u - (0.5 + 2.5 * s)).abs() < 1e-14);
        let d1 = f.eval_case1(s, 1).unwrap().unwrap();
        let t = Vector::new(-u.sin(), u.cos(), 0.0) * l;
        assert!((Vector::new(d1[0], d1[1], d1[2]) - t).length() < 1e-12);
        let d2 = f.eval_case1(s, 2).unwrap().unwrap();
        let n = Vector::new(-u.cos(), -u.sin(), 0.0) * (l * l / r);
        assert!((Vector::new(d2[0], d2[1], d2[2]) - n).length() < 1e-12);
        let d0 = f.eval_case1(s, 0).unwrap().unwrap();
        assert!((d0[0] - r * u.cos()).abs() < 1e-14);
        assert!(f.eval_case1(s, 3).unwrap().is_none());
    }
    f.trim(0.2, 0.7, 1e-9).unwrap();
    assert_eq!((f.first_parameter(), f.last_parameter()), (0.2, 0.7));
    assert_eq!(f.intervals(Continuity::C0), vec![0.2, 0.7]);
    assert!((f.get_s_parameter(f.curve().first_parameter()) - 0.2).abs() < 1e-12);
    assert!((f.curve().last_parameter() - (0.5 + 2.5 * 0.7)).abs() < 1e-14);
    assert!((f.get_length() - l).abs() < 1e-14);
    // A trim shorter than its tolerance changes nothing.
    f.trim(0.3, 0.3 + 1e-12, 1e-9).unwrap();
    assert_eq!(f.first_parameter(), 0.2);
}

#[test]
fn curvilinear_intervals_and_derivatives_of_a_spline() {
    let curve = cornered_cubic();
    let a = whole(&curve);
    let tol = 1e-8;
    let f = CurvlinFunc::new(&a, tol);
    let total = reference(&a, 0.0, 4.0);
    let iv = f.intervals(Continuity::C2);
    assert_eq!(f.nb_intervals(Continuity::C2), 3);
    assert_eq!(iv.len(), 4);
    assert_eq!((iv[0], iv[3]), (0.0, 1.0));
    for (k, knot) in [(1, 1.0), (2, 2.0)] {
        assert!((iv[k] - reference(&a, 0.0, knot) / total).abs() < 10.0 * tol / total);
    }
    // First derivative: length L along the tangent; second: by finite differences of the first.
    for s in [0.1, 0.3, 0.62, 0.9] {
        let d1 = f.eval_case1(s, 1).unwrap().unwrap();
        assert!((Vector::new(d1[0], d1[1], d1[2]).length() - f.get_length()).abs() < 1e-9 * total);
        let h = 1e-5;
        let (p, m) = (f.eval_case1(s + h, 1).unwrap().unwrap(), f.eval_case1(s - h, 1).unwrap().unwrap());
        let d2 = f.eval_case1(s, 2).unwrap().unwrap();
        for i in 0..3 {
            let fd = (p[i] - m[i]) / (2.0 * h);
            assert!((fd - d2[i]).abs() < 1e-4 * (1.0 + d2[i].abs()), "{s} {i}: {fd} vs {}", d2[i]);
        }
    }
}

#[test]
fn owning_a_view() {
    let curve = quintic();
    let f = CurvlinFunc::owning(Box::new(whole(&curve)), 1e-7);
    let u = f.get_u_parameter(0.4).unwrap();
    assert!((f.get_s_parameter(u) - 0.4).abs() < 1e-7);
}
