//! A pcurve made to follow its 3D curve parameter for parameter (what OCCT's
//! `BRepLib::SameParameter` sees to, by `Approx_SameParameter`), kentron's own way.
//!
//! A pcurve projected point by point and interpolated runs along the right path but not
//! always at the right parameters: where the 3D curve's speed jumps (a rational circle's
//! double knots), a smooth interpolant of the projected points cannot follow it. Here the
//! curve's points at sample parameters are projected onto the surface (each from the old
//! pcurve's point, so in the same period and away from singular points as it is), and a
//! cubic B-spline fitted to them by least squares on knots that keep the curve's own breaks
//! (as many derivatives continuous as there), its spans halved where it strays, until it
//! follows the curve within the target or stops gaining.

use crate::math::{Point, Point2};

use super::curve::Curve;
use super::curve2::Curve2;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// Sample points fitted in each span, between its ends.
const SAMPLES: usize = 5;
/// Rounds of halving the spans that stray.
const ROUNDS: usize = 12;
/// The most spans the fit has.
const MAX_SPANS: usize = 400;

/// The largest distance of the curve on the surface `pcurve` from `curve` at the same
/// parameters, sampled at `params`.
fn deviation(curve: &Curve, surface: &Surface, pcurve: &Curve2, params: &[f64]) -> f64 {
    params.iter().map(|&t| {
        let q = pcurve.point(t);
        surface.point(q.x, q.y).distance(curve.point(t))
    }).fold(0.0, f64::max)
}

/// The breaks of `curve` strictly inside `(first, last)` and how many derivatives are
/// continuous there (0: only the point).
fn breaks(curve: &Curve, first: f64, last: f64) -> Vec<(f64, usize)> {
    let margin = 1e-9 * (last - first).abs().max(1.0);
    match curve {
        Curve::BSpline(b) => b
            .knots()
            .iter()
            .zip(b.multiplicities())
            .filter(|(k, _)| **k > first + margin && **k < last - margin)
            .map(|(&k, &m)| (k, b.degree().saturating_sub(m)))
            .collect(),
        Curve::Trimmed(t) => breaks(t.basis(), first, last),
        _ => Vec::new(),
    }
}

/// The knot lines of `surface` (a B-spline's) where it is less smooth than cubic: `(along
/// v?, knot, continuous derivatives)`.
fn surface_knots(surface: &Surface) -> Vec<(bool, f64, usize)> {
    let mut surface = surface;
    while let Surface::Trimmed(t) = surface {
        surface = t.basis();
    }
    let Surface::BSpline(b) = surface else { return Vec::new() };
    let mut out = Vec::new();
    for (along_v, knots, mults, degree) in [(false, b.u_knots(), b.u_multiplicities(), b.u_degree()), (true, b.v_knots(), b.v_multiplicities(), b.v_degree())] {
        let n = knots.len();
        for (i, (&k, &m)) in knots.iter().zip(mults).enumerate() {
            let smooth = degree.saturating_sub(m);
            if i > 0 && i + 1 < n && smooth < 2 {
                out.push((along_v, k, smooth));
            }
        }
    }
    out
}

/// Where the curve's feet on `surface` cross the surface's less smooth knot lines, between
/// `params` (with their feet `uv`): the parameters, found by halving, and how smooth there.
fn crossings(curve: &Curve, surface: &Surface, pcurve: &Curve2, params: &[f64], uv: &[Point2]) -> Vec<(f64, usize)> {
    let knots = surface_knots(surface);
    let mut out = Vec::new();
    if knots.is_empty() {
        return out;
    }
    let at = |q: Point2, along_v: bool| if along_v { q.y } else { q.x };
    for i in 1..params.len() {
        for &(along_v, k, smooth) in &knots {
            let (a, b) = (at(uv[i - 1], along_v) - k, at(uv[i], along_v) - k);
            if a == 0.0 || a * b >= 0.0 {
                continue;
            }
            let (mut lo, mut hi, mut slo) = (params[i - 1], params[i], a);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    break;
                }
                let q = foot(surface, curve.point(mid), pcurve.point(mid)).unwrap_or_else(|| pcurve.point(mid));
                let f = at(q, along_v) - k;
                if (f < 0.0) == (slo < 0.0) {
                    (lo, slo) = (mid, f);
                } else {
                    hi = mid;
                }
            }
            out.push((0.5 * (lo + hi), smooth));
        }
    }
    out
}

/// The point of `surface` nearest `p`, searched from `seed` by Gauss-Newton steps (`None`
/// where the surface is singular there or the search runs off).
fn foot(surface: &Surface, p: Point, seed: Point2) -> Option<Point2> {
    let (mut u, mut v) = (seed.x, seed.y);
    let mut best = surface.point(u, v).distance(p);
    for _ in 0..30 {
        let (s, su, sv) = surface.d1(u, v);
        let r = p - s;
        let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a * c - b * b;
        if !(det > 1e-12 * a * c) || !(a > 0.0 && c > 0.0) {
            return None;
        }
        let (gu, gv) = (su.dot(r), sv.dot(r));
        let (du, dv) = ((c * gu - b * gv) / det, (a * gv - b * gu) / det);
        // (Steps halved while they make it worse.)
        let mut step = 1.0;
        let mut moved = false;
        for _ in 0..8 {
            let (nu, nv) = (u + step * du, v + step * dv);
            let d = surface.point(nu, nv).distance(p);
            if d <= best {
                (u, v, best, moved) = (nu, nv, d, true);
                break;
            }
            step *= 0.5;
        }
        let size = (du * du * a + dv * dv * c).sqrt();
        if !moved || step * size < 1e-13 * (1.0 + best) {
            break;
        }
    }
    let (_, su, sv) = surface.d1(u, v);
    // (Not where the step from the seed is long next to how far it was: it ran off.)
    let travel = (su * (u - seed.x) + sv * (v - seed.y)).length();
    let start = surface.point(seed.x, seed.y).distance(p);
    (travel <= 10.0 * start + 1e-9).then_some(Point2::new(u, v))
}

/// `pcurve` (of `curve` over `[first, last]`, on `surface`) remade where it strays from the
/// curve at the same parameters by more than `target`: the pcurve and its largest distance
/// from the curve. Its ends stay as they are.
pub(crate) fn same_parameter(curve: &Curve, first: f64, last: f64, surface: &Surface, pcurve: &Curve2, target: f64) -> (Curve2, f64) {
    let mut breaks = breaks(curve, first, last);
    // The pieces the knots bound: those of the curve's breaks to begin with.
    let mut ends: Vec<f64> = std::iter::once(first).chain(breaks.iter().map(|b| b.0)).chain(std::iter::once(last)).collect();
    let checks = |ends: &[f64]| -> Vec<f64> {
        ends.windows(2).flat_map(|w| (0..=2 * SAMPLES + 2).map(move |i| w[0] + (w[1] - w[0]) * i as f64 / (2 * SAMPLES + 2) as f64)).collect()
    };
    let before = deviation(curve, surface, pcurve, &checks(&ends));
    if before <= target {
        return (pcurve.clone(), before);
    }
    // (And where it crosses the surface's knots of less continuity: there too the pcurve's
    // speed jumps.)
    if !surface_knots(surface).is_empty() {
        let n = 64 * (breaks.len() + 1);
        let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
        let uv: Vec<Point2> = params.iter().map(|&t| foot(surface, curve.point(t), pcurve.point(t)).unwrap_or_else(|| pcurve.point(t))).collect();
        let margin = 1e-9 * (last - first).abs().max(1.0);
        for (t, smooth) in crossings(curve, surface, pcurve, &params, &uv) {
            if t <= first + margin || t >= last - margin {
                continue;
            }
            match breaks.iter_mut().find(|b| (b.0 - t).abs() <= margin) {
                Some(b) => b.1 = b.1.min(smooth),
                None => breaks.push((t, smooth)),
            }
        }
        breaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    ends = std::iter::once(first).chain(breaks.iter().map(|b| b.0)).chain(std::iter::once(last)).collect();
    let (start, end) = (pcurve.point(first), pcurve.point(last));
    let mut best = (pcurve.clone(), before);
    let (mut previous, mut stalled) = (f64::INFINITY, 0);
    for _ in 0..ROUNDS {
        // The samples, the curve's points' feet on the surface.
        let mut params = Vec::new();
        let mut feet: Vec<[f64; 3]> = Vec::new();
        for (k, w) in ends.windows(2).enumerate() {
            for i in usize::from(k > 0)..=SAMPLES + 1 {
                let t = w[0] + (w[1] - w[0]) * i as f64 / (SAMPLES + 1) as f64;
                let seed = if t == first { start } else if t == last { end } else { pcurve.point(t) };
                let uv = if t == first || t == last { seed } else { foot(surface, curve.point(t), seed).unwrap_or(seed) };
                params.push(t);
                feet.push([uv.x, uv.y, 0.0]);
            }
        }
        // Cubic, of the pieces' ends as knots: the curve's breaks as continuous as it is
        // there, those of halving fully.
        let mut flat = vec![first; 4];
        for &e in &ends[1..ends.len() - 1] {
            let smooth = breaks.iter().find(|b| b.0 == e).map_or(2, |b| b.1.min(2));
            flat.extend(std::iter::repeat_n(e, 3 - smooth));
        }
        flat.extend([last; 4]);
        let Some(fit) = super::fit::least_squares_with_knots::<Point2>(&feet, &params, 3, flat) else { break };
        let fit = Curve2::BSpline(fit);
        // The pieces that stray, halved.
        let mut worst: f64 = 0.0;
        let mut next = vec![first];
        for w in ends.windows(2) {
            let d = deviation(curve, surface, &fit, &checks(w));
            worst = worst.max(d);
            if d > target && next.len() + ends.len() < MAX_SPANS {
                next.push(0.5 * (w[0] + w[1]));
            }
            next.push(w[1]);
        }
        if worst < 0.9 * previous {
            (previous, stalled) = (worst, 0);
        } else {
            stalled += 1;
        }
        if worst < best.1 {
            best = (fit, worst);
        }
        // (Done when within the target, or when halving has stopped gaining.)
        if worst <= target || next.len() == ends.len() || stalled == 3 {
            break;
        }
        ends = next;
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::{BSplineCurve, BSplineSurface, KnotVector};
    use crate::math::Point;

    /// A rational quadratic circle (double knots: its speed jumps there) on a bent B-spline
    /// surface, with a pcurve through its points at the wrong parameters: remade to follow it.
    #[test]
    fn rational_circle_followed() {
        // (Flat, but unevenly parametrized.)
        let xs = [-3.0, -2.0, 0.5, 3.0];
        let poles: Vec<Vec<Point>> = (0..4).map(|i| (0..4).map(|j| Point::new(xs[i], xs[3 - j] * -1.0, 0.0)).collect()).collect();
        let kv = KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3);
        let surface = Surface::BSpline(BSplineSurface::new(poles, kv.clone(), kv));
        let h = std::f64::consts::FRAC_1_SQRT_2;
        let pts = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0), (-1.0, -1.0), (0.0, -1.0), (1.0, -1.0), (1.0, 0.0)];
        let circle = BSplineCurve::new(pts.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect(), vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![3, 2, 2, 2, 3], 2)
            .with_weights((0..9).map(|i| if i % 2 == 1 { h } else { 1.0 }).collect());
        let curve = Curve::BSpline(circle);
        // The feet at evenly spread parameters, interpolated smoothly (as the projection does).
        let n = 40;
        let mut uv = Vec::new();
        let mut params = Vec::new();
        for i in 0..=n {
            let t = 4.0 * i as f64 / n as f64;
            let p = curve.point(t);
            let guess = Point2::new(0.5, 0.5);
            uv.push(foot(&surface, p, guess).unwrap());
            params.push(t);
        }
        let smooth = Curve2::BSpline(crate::geometry::fit::interpolate_with_parameters(&uv, &params).unwrap());
        let (fixed, dev) = same_parameter(&curve, 0.0, 4.0, &surface, &smooth, 1e-7);
        let check: Vec<f64> = (0..=997).map(|i| 4.0 * i as f64 / 997.0).collect();
        let was = deviation(&curve, &surface, &smooth, &check);
        let now = deviation(&curve, &surface, &fixed, &check);
        assert!(was > 1e-4, "the smooth pcurve strays {was}");
        assert!(now < 1e-6 && dev < 1e-6, "remade: {now} ({dev})");
    }
}
