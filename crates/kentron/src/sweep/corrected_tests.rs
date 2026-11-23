//! The rotation-minimizing frame against references computed independently of the law: a
//! fine fourth-order Runge–Kutta integration of the transport equation `N′ = −(N·T′) T`, and,
//! for curves with non-vanishing curvature, the closed form `N = N_F` turned about `T` by
//! `−∫τ ds` (the Frenet normal and the curve's own torsion).

use super::*;
use crate::geometry::curve::{BSplineCurve, Line};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Point, Transformable};

/// Runge–Kutta steps over the whole range of a piece, at least.
const REFERENCE_STEPS: usize = 200_000;

/// The rotation-minimizing frames of `curve` travelled from `first` to `last` at
/// `fractions` (each in [0, 1]), by Runge–Kutta integration of parallel transport from the
/// Frenet normal at the start; on a closed piece (ends within `1e-7` of its length,
/// tangents within `1e-9` rad) the holonomy is taken out by a twist even in arc length.
/// Tangent and normal at each fraction.
pub(crate) fn rmf_reference(curve: &Curve, first: f64, last: f64, fractions: &[f64]) -> Vec<(Vector, Vector)> {
    let sign = if last >= first { 1.0 } else { -1.0 };
    let (lo, hi) = (first.min(last), first.max(last));
    let whole = GeomAdaptor::new(curve, lo, hi);
    let breaks = whole.intervals(Continuity::C2);
    // (Each stop: a parameter, and the fraction it reports, if any; in the order of travel.)
    let mut stops: Vec<(f64, Option<usize>)> = breaks.iter().map(|&b| (b, None)).collect();
    for (i, &f) in fractions.iter().enumerate() {
        stops.push((if f == 1.0 { last } else { first + (last - first) * f }, Some(i)));
    }
    stops.sort_by(|a, b| (sign * a.0).total_cmp(&(sign * b.0)));
    // (The C2 stretch between two stops, for one-sided derivatives at its ends.)
    let stretch = |x0: f64, x1: f64| -> GeomAdaptor {
        let m = 0.5 * (x0 + x1);
        let k = breaks.partition_point(|&b| b <= m).clamp(1, breaks.len() - 1);
        whole.trim(breaks[k - 1], breaks[k])
    };
    // (The unit tangent along the travel, its rate in x, and the speed.)
    let tangent_rate = |view: &GeomAdaptor, x: f64| -> (Vector, Vector, f64) {
        let (_, d1, d2) = view.d2(x);
        let speed = d1.length();
        let d = d1 * (1.0 / speed);
        let rate = (d2 - d * d2.dot(d)) * (sign / speed);
        (d * sign, rate, speed)
    };
    let start_normal = FrenetLaw::new(curve, first, last).at_fraction(0.0).1;
    let mut normal = start_normal;
    let t0 = tangent_rate(&stretch(stops[0].0, stops[1].0), stops[0].0).0;
    let mut out = vec![(Vector::ZERO, Vector::ZERO, 0.0); fractions.len()];
    let mut length = 0.0;
    let range = (last - first).abs();
    let mut previous = stops[0].0;
    let mut tangent = t0;
    if let Some(i) = stops[0].1 {
        out[i] = (t0, normal, 0.0);
    }
    for &(x, report) in &stops[1..] {
        if x != previous {
            let view = stretch(previous, x);
            let steps = ((REFERENCE_STEPS as f64 * (x - previous).abs() / range).ceil() as usize).max(8);
            let h = (x - previous) / steps as f64;
            let rhs = |x: f64, n: Vector| -> (Vector, f64) {
                let (t, rate, speed) = tangent_rate(&view, x);
                (t * (-n.dot(rate)), speed)
            };
            // (Across a break: the least rotation between the tangents either side.)
            let t_in = tangent_rate(&view, previous).0;
            normal = super::super::frenet::least_rotation(normal, tangent, t_in);
            for k in 0..steps {
                let x0 = previous + h * k as f64;
                let x1 = if k + 1 == steps { x } else { x0 + h };
                let hh = x1 - x0;
                let (k1, s0) = rhs(x0, normal);
                let (k2, sm) = rhs(x0 + 0.5 * hh, normal + k1 * (0.5 * hh));
                let (k3, _) = rhs(x0 + 0.5 * hh, normal + k2 * (0.5 * hh));
                let (k4, s1) = rhs(x1, normal + k3 * hh);
                normal = normal + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (hh / 6.0);
                let t = tangent_rate(&view, x1).0;
                normal = (normal - t * t.dot(normal)).normalized();
                length += (s0 + 4.0 * sm + s1) * hh.abs() / 6.0;
            }
            tangent = tangent_rate(&view, x).0;
            previous = x;
        }
        if let Some(i) = report {
            out[i] = (tangent, normal, length);
        }
    }
    // The closure: the holonomy taken out evenly in arc length.
    let gap = curve.point(first).distance(curve.point(last));
    let closed = gap <= 1e-7 * length && t0.cross(tangent).length().atan2(t0.dot(tangent)) <= 1e-9;
    if closed {
        let holonomy = signed_angle(start_normal, normal, tangent);
        return out.into_iter().map(|(t, n, s)| (t, rotate_about(n, t, -holonomy * s / length))).collect();
    }
    out.into_iter().map(|(t, n, _)| (t, n)).collect()
}

/// The largest sine of the angle between the law's frames and the reference's at
/// `fractions`, checking that the law's are orthonormal.
fn worst_against_reference(curve: &Curve, first: f64, last: f64, fractions: &[f64]) -> f64 {
    let law = CorrectedFrenet::new(curve, first, last).unwrap();
    let reference = rmf_reference(curve, first, last, fractions);
    let mut worst = 0.0f64;
    for (&f, (t, n)) in fractions.iter().zip(reference) {
        let (lt, ln) = law.at_fraction(f);
        worst = worst.max(lt.cross(t).length()).max(ln.cross(n).length());
        assert!(ln.dot(n) > 0.0 && lt.dot(t) > 0.0, "at {f}: {ln:?} vs {n:?}");
        assert!(lt.dot(ln).abs() < 1e-12 && (ln.length() - 1.0).abs() < 1e-12 && (lt.length() - 1.0).abs() < 1e-12, "at {f}: not orthonormal");
    }
    worst
}

fn fractions(n: usize) -> Vec<f64> {
    (0..=n).map(|i| i as f64 / n as f64).collect()
}

/// A B-spline through points of the helix `(r cos t, r sin t, c t)`, `t` from 0 over
/// `turns` turns.
pub(crate) fn helix(r: f64, c: f64, turns: f64, points: usize) -> Curve {
    let p: Vec<Point> = (0..=points)
        .map(|i| {
            let t = std::f64::consts::TAU * turns * i as f64 / points as f64;
            Point::new(r * t.cos(), r * t.sin(), c * t)
        })
        .collect();
    Curve::BSpline(crate::geometry::fit::interpolate(&p, false, 1e-9).unwrap())
}

/// Along a helix the rotation-minimizing normal is the Frenet normal turned about the
/// tangent by `−∫τ ds`: the closed form with the curve's own torsion (integrated by Gauss
/// rules per knot span), to 1e-9 rad; and that angle is the helix's `−τ s` to the accuracy
/// of the fit.
#[test]
fn along_a_helix_the_normal_lags_the_frenet_normal_by_the_torsion_integral() {
    let (r, c) = (3.0, 0.8);
    let curve = helix(r, c, 3.0, 240);
    let Curve::BSpline(spline) = &curve else { unreachable!() };
    let (first, last) = (spline.first_parameter(), spline.last_parameter());
    let law = CorrectedFrenet::new(&curve, first, last).unwrap();
    let frenet = FrenetLaw::new(&curve, first, last);
    // (Torsion times speed, and speed.)
    let integrand = |u: f64| {
        let (_, d1, d2, d3) = curve.d3(u);
        let x = d1.cross(d2);
        (x.dot(d3) / x.length_squared() * d1.length(), d1.length())
    };
    let knots: Vec<f64> = spline.knots().to_vec();
    let integral = |to: f64| {
        let (mut angle, mut length) = (0.0, 0.0);
        for w in knots.windows(2) {
            let (a, b) = (w[0], w[1].min(to));
            if b <= a {
                break;
            }
            let pieces = 4;
            for k in 0..pieces {
                let (p, q) = (a + (b - a) * k as f64 / pieces as f64, a + (b - a) * (k + 1) as f64 / pieces as f64);
                let (m, h) = (0.5 * (p + q), 0.5 * (q - p));
                for &(x, wt) in crate::geometry::measure::GAUSS.iter() {
                    let (ts, s) = integrand(m + h * x);
                    angle += wt * h * ts;
                    length += wt * h * s;
                }
            }
        }
        (angle, length)
    };
    let tau = c / (r * r + c * c);
    let mut worst = 0.0f64;
    for f in fractions(97) {
        let u = if f == 1.0 { last } else { first + (last - first) * f };
        let (t, nf) = frenet.at_fraction(f);
        let (angle, s) = integral(u);
        let expected = rotate_about(nf, t, -angle);
        let (lt, n) = law.at_fraction(f);
        assert!(lt.cross(t).length() < 1e-12);
        worst = worst.max(n.cross(expected).length());
        assert!(n.dot(expected) > 0.0);
        // (The fit's torsion is the helix's, near enough.)
        assert!((angle - tau * s).abs() < 1e-3 * (1.0 + tau * s), "{angle} vs {}", tau * s);
    }
    assert!(worst < 1e-9, "{worst:e}");
}

/// The law against the Runge–Kutta transport along non-planar curves: a helix both ways,
/// a twisted cubic with a C1 knot both ways, and a Bézier curve.
#[test]
fn the_law_follows_parallel_transport() {
    use crate::geometry::curve::BezierCurve;
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let c1_knot = Curve::BSpline(BSplineCurve::new(
        vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 2.0, 1.0), p(3.0, 0.0, 2.0), p(4.0, -1.0, 1.0), p(5.0, 1.0, 0.0), p(6.0, 2.0, 2.0)],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 2, 1, 4],
        3,
    ));
    let bezier = Curve::Bezier(BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 1.0, 2.0), p(2.0, -1.0, 3.0), p(0.0, -1.0, 4.0)]));
    let helix = helix(2.0, 0.5, 2.0, 120);
    let Curve::BSpline(h) = &helix else { unreachable!() };
    let (h0, h1) = (h.first_parameter(), h.last_parameter());
    for (k, (curve, first, last)) in [(&c1_knot, 0.0, 3.0), (&c1_knot, 3.0, 0.0), (&bezier, 0.0, 1.0), (&helix, h0, h1), (&helix, h1, h0)].into_iter().enumerate() {
        let worst = worst_against_reference(curve, first, last, &fractions(61));
        assert!(worst < 1e-9, "case {k}: {worst:e}");
    }
}

/// A closed space curve (a periodic trefoil, both ways): the frame comes back to its start,
/// and follows the transport with the holonomy taken out evenly in arc length.
#[test]
fn a_closed_spine_closes_its_frame() {
    let mut model = crate::topology::Model::new();
    let shape = crate::io::brep::read(&mut model, &mut include_str!("tests/occt_periodic_trefoil.brep").as_bytes()).unwrap();
    let edge = model.unique_subshapes(shape, crate::topology::NodeKind::Edge)[0];
    let (curve, t0, t1) = model.edge_curve(edge).unwrap();
    for (first, last) in [(t0, t1), (t1, t0)] {
        let law = CorrectedFrenet::new(&curve, first, last).unwrap();
        let (a, b) = (law.at_fraction(0.0), law.at_fraction(1.0));
        assert!(a.0.cross(b.0).length() < 1e-9 && a.1.cross(b.1).length() < 1e-9 && a.1.dot(b.1) > 0.0, "{a:?} vs {b:?}");
        assert!(law.closure.is_some_and(|c| c.angle.abs() > 1e-3), "a twisted closed curve has holonomy");
        let worst = worst_against_reference(&curve, first, last, &fractions(53));
        assert!(worst < 1e-9, "{worst:e}");
    }
}

/// Along planar curves the normal stays in the plane: the Frenet normal up to sign (here,
/// an S with an inflection, in a tilted plane, both ways).
#[test]
fn along_a_planar_curve_the_normal_stays_in_the_plane() {
    let p = |x: f64, y: f64| Point::new(x, y, 0.0);
    let s = Curve::BSpline(BSplineCurve::new(vec![p(0.0, 0.0), p(1.0, 2.0), p(2.0, 2.0), p(3.0, -2.0), p(4.0, -2.0), p(5.0, 0.0)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], 3));
    let tilt = crate::math::Transform::rotation(crate::math::Axis::new(Point::ORIGIN, Direction::from_vector(Vector::new(1.0, 1.0, 0.3))), 0.7);
    let s = s.transformed(&tilt);
    let plane_normal = Vector::new(0.0, 0.0, 1.0).transformed(&tilt);
    for (first, last) in [(0.0, 3.0), (3.0, 0.0)] {
        let law = CorrectedFrenet::new(&s, first, last).unwrap();
        let frenet = FrenetLaw::new(&s, first, last);
        for f in fractions(40) {
            let (t, n) = law.at_fraction(f);
            assert!(n.dot(plane_normal).abs() < 1e-12, "at {f}: {n:?}");
            let (ft, fnormal) = frenet.at_fraction(f);
            assert!(t.cross(ft).length() < 1e-12 && n.cross(fnormal).length() < 1e-12, "at {f}");
        }
    }
}

/// Along a line the frame is constant: the tangent and its default normal.
#[test]
fn along_a_line_the_frame_is_constant() {
    let line = Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::NEG_Y));
    let law = CorrectedFrenet::new(&line, 0.0, 5.0).unwrap();
    for f in fractions(7) {
        let (t, n) = law.at_fraction(f);
        assert_eq!((t, n), (Vector::new(0.0, -1.0, 0.0), Vector::new(0.0, 0.0, -1.0)));
    }
    assert!(CorrectedFrenet::new(&line, 2.0, 2.0).is_err());
}
