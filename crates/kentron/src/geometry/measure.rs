//! Lengths along curves.

use super::curve::Curve;
use super::curve2::Curve2;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// Gauss-Legendre nodes and weights on [-1, 1] (order 10).
const GAUSS: [(f64, f64); 10] = [
    (-0.9739065285171717, 0.0666713443086881),
    (-0.8650633666889845, 0.1494513491505806),
    (-0.6794095682990244, 0.2190863625159820),
    (-0.4333953941292472, 0.2692667193099963),
    (-0.1488743389816312, 0.2955242247147529),
    (0.1488743389816312, 0.2955242247147529),
    (0.4333953941292472, 0.2692667193099963),
    (0.6794095682990244, 0.2190863625159820),
    (0.8650633666889845, 0.1494513491505806),
    (0.9739065285171717, 0.0666713443086881),
];

fn gauss(speed: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let (m, h) = (0.5 * (a + b), 0.5 * (b - a));
    GAUSS.iter().map(|&(x, w)| w * speed(m + h * x)).sum::<f64>() * h
}

/// The integral of `speed` over `[a, b]`, adaptively: halved until two halves agree with
/// the whole to a relative 1e-13, or to `floor` (where the piece adds next to nothing to the
/// whole integral, its own relative error being rounding noise).
fn adaptive(speed: &dyn Fn(f64) -> f64, a: f64, b: f64, whole: f64, depth: usize, floor: f64) -> f64 {
    let m = 0.5 * (a + b);
    let (left, right) = (gauss(speed, a, m), gauss(speed, m, b));
    if depth == 0 || (left + right - whole).abs() <= (1e-13 * (left + right).abs()).max(floor).max(1e-300) {
        return left + right;
    }
    adaptive(speed, a, m, left, depth - 1, floor) + adaptive(speed, m, b, right, depth - 1, floor)
}

/// The pieces `[t1, t2]` splits into at the breaks (each period's, for a periodic curve).
fn pieces(breaks: &[f64], period: Option<f64>, t1: f64, t2: f64) -> Vec<f64> {
    let mut cuts = vec![t1, t2];
    for &k in breaks {
        match period {
            Some(p) => {
                let mut t = k + ((t1 - k) / p).ceil() * p;
                while t < t2 {
                    cuts.push(t);
                    t += p;
                }
            }
            None if t1 < k && k < t2 => cuts.push(k),
            None => {}
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    cuts
}

/// The arc length of a curve with speed `speed` from `t1` to `t2` (negative if `t2 < t1`).
fn arc_length(speed: &dyn Fn(f64) -> f64, breaks: &[f64], period: Option<f64>, t1: f64, t2: f64) -> f64 {
    let (lo, hi, sign) = if t1 <= t2 { (t1, t2, 1.0) } else { (t2, t1, -1.0) };
    let cuts = pieces(breaks, period, lo, hi);
    // Several pieces per span keep the quadrature's error small on long ranges.
    let n = 8;
    let parts: Vec<(f64, f64, f64)> = cuts
        .windows(2)
        .flat_map(|w| (0..n).map(move |i| (w[0] + (w[1] - w[0]) * i as f64 / n as f64, w[0] + (w[1] - w[0]) * (i + 1) as f64 / n as f64)))
        .map(|(a, b)| (a, b, gauss(speed, a, b)))
        .collect();
    // (A floor for each piece's error: a part in 1e-15 of the whole.)
    let floor = 1e-15 * parts.iter().map(|p| p.2.abs()).sum::<f64>();
    let total: f64 = parts.iter().map(|&(a, b, whole)| adaptive(speed, a, b, whole, 20, floor)).sum();
    sign * total
}

/// The parameter at arc length `distance` from `from`: Newton steps on the length,
/// kept inside a bracket.
fn parameter_at(speed: &dyn Fn(f64) -> f64, breaks: &[f64], period: Option<f64>, from: f64, distance: f64) -> f64 {
    if distance == 0.0 {
        return from;
    }
    let sign = distance.signum();
    // A bracket: step out until the length passes the distance.
    let mut step = (distance.abs() / speed(from).max(1e-12)).max(1e-9);
    let mut far = from + sign * step;
    for _ in 0..200 {
        if arc_length(speed, breaks, period, from, far).abs() >= distance.abs() {
            break;
        }
        step *= 2.0;
        far = from + sign * step;
    }
    let (mut lo, mut hi) = if sign > 0.0 { (from, far) } else { (far, from) };
    let mut t = from + sign * (distance.abs() / speed(from).max(1e-12)).min(step);
    for _ in 0..100 {
        let f = arc_length(speed, breaks, period, from, t) - distance;
        if f.abs() <= 1e-12 * distance.abs().max(1.0) {
            return t;
        }
        // f grows with t.
        if f > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
        let v = speed(t);
        let next = if v > 0.0 { t - f / v } else { 0.5 * (lo + hi) };
        t = if next > lo && next < hi { next } else { 0.5 * (lo + hi) };
        if hi - lo <= 1e-15 * (1.0 + t.abs()) {
            break;
        }
    }
    t
}

/// The length of a bounded curve over its whole domain.
#[doc(alias = "GCPnts_AbscissaPoint::Length")]
pub fn length(curve: &Curve) -> f64 {
    length_between(curve, curve.first_parameter(), curve.last_parameter())
}

/// The length of the part of `curve` between parameters `t1` and `t2`.
pub fn length_between(curve: &Curve, t1: f64, t2: f64) -> f64 {
    arc_length(&|t| curve.d1(t).1.length(), &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), t1, t2)
}

/// The parameter of the point at arc length `distance` from parameter `from` (backwards
/// for a negative distance).
#[doc(alias = "GCPnts_AbscissaPoint")]
pub fn parameter_at_length(curve: &Curve, from: f64, distance: f64) -> f64 {
    let speed = |t: f64| curve.d1(t).1.length();
    parameter_at(&speed, &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), from, distance)
}

/// [`length`] for a 2D curve.
pub fn length2(curve: &Curve2) -> f64 {
    length_between2(curve, curve.first_parameter(), curve.last_parameter())
}

pub fn length_between2(curve: &Curve2, t1: f64, t2: f64) -> f64 {
    arc_length(&|t| curve.d1(t).1.length(), &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), t1, t2)
}

pub fn parameter_at_length2(curve: &Curve2, from: f64, distance: f64) -> f64 {
    let speed = |t: f64| curve.d1(t).1.length();
    parameter_at(&speed, &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), from, distance)
}

/// The largest distance between a 3D curve and a curve on a surface over a common
/// parameter range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveDeviation {
    pub max_distance: f64,
    /// A parameter where the maximum is reached.
    pub parameter: f64,
}

/// How far `curve` strays from `surface(pcurve(t))` for `t` in `[t1, t2]` (comparing
/// points at equal parameters). The maximum is searched globally, so it is found even
/// when the deviation oscillates quickly.
#[doc(alias = "GeomLib_CheckCurveOnSurface")]
pub fn curve_on_surface_deviation(
    curve: &Curve,
    pcurve: &Curve2,
    surface: &Surface,
    t1: f64,
    t2: f64,
) -> CurveDeviation {
    let distance = |t: f64| {
        let q = pcurve.point(t);
        surface.point(q.x, q.y).distance(curve.point(t))
    };
    let (lo, hi) = (t1.min(t2), t1.max(t2));
    // Samples, then golden-section refinement around each local maximum.
    const N: usize = 256;
    let ts: Vec<f64> = (0..=N).map(|i| lo + (hi - lo) * i as f64 / N as f64).collect();
    let ds: Vec<f64> = ts.iter().map(|&t| distance(t)).collect();
    let mut best = CurveDeviation { max_distance: ds[0], parameter: ts[0] };
    for i in 0..=N {
        let left = if i > 0 { ds[i - 1] } else { f64::NEG_INFINITY };
        let right = if i < N { ds[i + 1] } else { f64::NEG_INFINITY };
        if ds[i] < left || ds[i] < right {
            continue;
        }
        let (mut a, mut b) = (ts[i.saturating_sub(1)], ts[(i + 1).min(N)]);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..60 {
            let (x1, x2) = (b - g * (b - a), a + g * (b - a));
            if distance(x1) > distance(x2) {
                b = x2;
            } else {
                a = x1;
            }
        }
        let mut t = 0.5 * (a + b);
        // Then exactly, where the distance's slope (S - C)·(dS/dt - C') changes sign.
        let slope = |t: f64| {
            let (q, dq) = pcurve.d1(t);
            let (sp, su, sv) = surface.d1(q.x, q.y);
            let (cp, dc) = curve.d1(t);
            (sp - cp).dot(su * dq.x + sv * dq.y - dc)
        };
        let (mut lo, mut hi) = (ts[i.saturating_sub(1)], ts[(i + 1).min(N)]);
        if slope(lo) > 0.0 && slope(hi) < 0.0 {
            for _ in 0..200 {
                let m = 0.5 * (lo + hi);
                if slope(m) > 0.0 {
                    lo = m;
                } else {
                    hi = m;
                }
                if hi - lo <= 1e-15 * (1.0 + m.abs()) {
                    break;
                }
            }
            t = 0.5 * (lo + hi);
        }
        for (t, d) in [(t, distance(t)), (ts[i], ds[i])] {
            if d > best.max_distance {
                best = CurveDeviation { max_distance: d, parameter: t };
            }
        }
    }
    best
}
