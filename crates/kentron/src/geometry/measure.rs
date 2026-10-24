//! Lengths along curves.

use super::curve::Curve;
use super::curve2::Curve2;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// Gauss-Legendre nodes and weights on [-1, 1] (order 10).
pub(crate) const GAUSS: [(f64, f64); 10] = [
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
/// kept inside a bracket, each length measured on from the last point it was known at; no
/// further than `bounds` (a bounded curve's ends, where it stops).
fn parameter_at(speed: &dyn Fn(f64) -> f64, breaks: &[f64], period: Option<f64>, bounds: Option<(f64, f64)>, from: f64, distance: f64) -> f64 {
    if distance == 0.0 {
        return from;
    }
    let sign = distance.signum();
    let target = distance.abs();
    let along = |a: f64, b: f64| arc_length(speed, breaks, period, a, b).abs();
    let bound = bounds.map(|(lo, hi)| if sign > 0.0 { hi } else { lo });
    // A bracket: step out until the length passes the distance (or the curve ends short of
    // it, or as good as: its end).
    let mut step = (target / speed(from).max(1e-12)).max(1e-9);
    let (mut near, mut near_length) = (from, 0.0);
    let (mut far, mut far_length);
    loop {
        far = from + sign * step;
        if let Some(b) = bound.filter(|&b| (far - b) * sign >= 0.0) {
            far = b;
            far_length = near_length + along(near, far);
            if far_length <= target * (1.0 + 1e-12) {
                return far;
            }
            break;
        }
        far_length = near_length + along(near, far);
        if far_length >= target || step > 1e300 {
            break;
        }
        (near, near_length) = (far, far_length);
        step *= 2.0;
    }
    let ((mut lo, mut hi), mut known) = (if sign > 0.0 { (near, far) } else { (far, near) }, (near, near_length));
    let mut t = (from + sign * (target / speed(from).max(1e-12)).min(step)).clamp(lo, hi);
    for _ in 0..100 {
        // (Away from `from` past the point last known, the length grows; back towards it, it
        // shrinks.)
        let (kt, kl) = known;
        let length = kl + if (t - kt) * sign >= 0.0 { along(kt, t) } else { -along(t, kt) };
        known = (t, length);
        let f = length - target;
        if f.abs() <= 1e-12 * target.max(1.0) {
            return t;
        }
        if (f > 0.0) == (sign > 0.0) {
            hi = t;
        } else {
            lo = t;
        }
        let v = speed(t) * sign;
        let next = if v != 0.0 { t - f / v } else { 0.5 * (lo + hi) };
        t = if next > lo && next < hi { next } else { 0.5 * (lo + hi) };
        if hi - lo <= 1e-15 * (1.0 + t.abs()) {
            break;
        }
    }
    t
}

/// The ends where `curve` stops, if bounded and not periodic.
fn stops(first: f64, last: f64, periodic: bool) -> Option<(f64, f64)> {
    (!periodic && first.is_finite() && last.is_finite()).then_some((first, last))
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
    let bounds = stops(curve.first_parameter(), curve.last_parameter(), curve.is_periodic());
    parameter_at(&speed, &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), bounds, from, distance)
}

/// Lengths along a curve from the start of a range up to each of the curve's breaks in it
/// (and its end), for points at given lengths found without measuring from the start each
/// time.
#[derive(Debug, Clone)]
pub struct LengthTable {
    /// The parameters: the range's ends and the breaks between.
    pub parameters: Vec<f64>,
    /// The length from the range's start to each.
    pub lengths: Vec<f64>,
}

impl LengthTable {
    /// The table of `curve` over `[first, last]`.
    pub fn new(curve: &Curve, first: f64, last: f64) -> LengthTable {
        let mut parameters = vec![first];
        parameters.extend(curve.breaks().into_iter().filter(|&t| t > first && t < last));
        parameters.push(last);
        let mut lengths = vec![0.0];
        for w in parameters.windows(2) {
            let l = lengths.last().expect("a length") + length_between(curve, w[0], w[1]);
            lengths.push(l);
        }
        LengthTable { parameters, lengths }
    }

    /// The whole length.
    pub fn total(&self) -> f64 {
        *self.lengths.last().expect("a length")
    }

    /// The parameter of `curve` (the table's) at `distance` along from the start: from the
    /// break before it.
    pub fn parameter_at(&self, curve: &Curve, distance: f64) -> f64 {
        let n = self.parameters.len();
        if distance <= 0.0 {
            return if distance == 0.0 { self.parameters[0] } else { parameter_at_length(curve, self.parameters[0], distance) };
        }
        if distance >= self.total() * (1.0 - 1e-14) {
            return if distance <= self.total() { self.parameters[n - 1] } else { parameter_at_length(curve, self.parameters[n - 1], distance - self.total()) };
        }
        let i = self.lengths.partition_point(|&l| l <= distance).clamp(1, n - 1) - 1;
        parameter_at_length(curve, self.parameters[i], distance - self.lengths[i])
    }
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
    let bounds = stops(curve.first_parameter(), curve.last_parameter(), curve.is_periodic());
    parameter_at(&speed, &curve.breaks(), curve.period().filter(|_| curve.is_periodic()), bounds, from, distance)
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

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::f64::consts::PI;

    use super::*;
    use crate::math::Point;

    #[test]
    fn parameters_at_lengths_are_found_from_where_the_length_was_last_known() {
        // Twenty turns of a helix (radius 2, pitch 1) as a B-spline: the parameter a length
        // along, forwards and backwards, out to the very end, found in about as many
        // evaluations of the speed as the whole length takes, not anew over the whole way at
        // every step (nor, at the end, searched for past it).
        let points: Vec<Point> = (0..=800).map(|i| {
            let a = 40.0 * PI * i as f64 / 800.0;
            Point::new(2.0 * a.cos(), 2.0 * a.sin(), a / (2.0 * PI))
        }).collect();
        let curve = Curve::BSpline(crate::geometry::fit::interpolate(&points, false, 1e-9).unwrap());
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let whole = length(&curve);
        let calls = Cell::new(0usize);
        let speed = |t: f64| {
            calls.set(calls.get() + 1);
            curve.d1(t).1.length()
        };
        for (from, distance) in [(first, whole), (first, 0.3 * whole), (last, -whole), (last, -0.5 * whole)] {
            calls.set(0);
            let t = parameter_at(&speed, &curve.breaks(), None, Some((first, last)), from, distance);
            let along = length_between(&curve, from, t);
            assert!((along - distance).abs() < 1e-9 * whole, "{along} not {distance}");
            assert!(calls.get() < 400_000, "{} evaluations", calls.get());
        }
    }

    #[test]
    fn a_table_of_lengths_finds_points_from_the_break_before() {
        // Five turns of a helix as a B-spline of many spans: the parameter a length along,
        // from the table, is where that length from the start ends.
        let points: Vec<Point> = (0..=200).map(|i| {
            let a = 10.0 * PI * i as f64 / 200.0;
            Point::new(2.0 * a.cos(), 2.0 * a.sin(), a / (2.0 * PI))
        }).collect();
        let curve = Curve::BSpline(crate::geometry::fit::interpolate(&points, false, 1e-9).unwrap());
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let table = LengthTable::new(&curve, first, last);
        assert!((table.total() - length(&curve)).abs() < 1e-9 * table.total());
        for f in [0.0, 1e-6, 0.137, 0.5, 0.731, 0.999999, 1.0] {
            let distance = f * table.total();
            let t = table.parameter_at(&curve, distance);
            assert!((length_between(&curve, first, t) - distance).abs() < 1e-9 * table.total(), "at {f}");
        }
    }
}
