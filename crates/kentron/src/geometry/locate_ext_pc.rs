//! The local extremum of the distance from a point to a curve whose parameter is nearest a
//! guess.
//!
//! Lines and conics list the stationary points of their whole (infinite or periodic) curve
//! in closed form — the foot of the perpendicular on a line; the projected direction and its
//! opposite on a circle; the roots of `(b² − a²) sin t cos t + a x sin t − b y cos t = 0` on
//! an ellipse; a cubic in `t` for a parabola and a quartic in `eᵗ` for a hyperbola — and
//! keep those in the range.
//!
//! Other curves are searched outward from the guess in windows that double in width, with
//! the sampling and root refinement of [`ext_pc`](super::ext_pc): sign changes of
//! `F = (C − P)·T`, corners at breakpoints, stretches where `F` vanishes, and range ends where
//! `F` vanishes. As soon as a window holds a candidate, the nearest one is certain (everything
//! nearer has been searched), so tracking a point along a curve from its previous projection
//! usually costs one small window.
//!
//! The answer is the candidate nearest the guess in parameter, the larger one of two equally
//! near; `None` when the range holds none.

use std::f64::consts::PI;

use crate::geometry::adaptor::{in_period, CurveAdaptor, CurveKind};
use crate::geometry::curve::Curve;
use crate::geometry::ext_pc::{parametric_tolerance, sort_and_merge, Candidate, CandidateKind, ExtPCPoint, Searcher};
use crate::geometry::polynomial_roots::{trigonometric_roots, PolynomialRoots};
use crate::math::{Frame, Point, Vector};

/// The extremum found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LocatedExtremum {
    pub parameter: f64,
    pub point: Point,
    pub square_distance: f64,
    pub is_min: bool,
}

impl From<ExtPCPoint> for LocatedExtremum {
    fn from(e: ExtPCPoint) -> LocatedExtremum {
        LocatedExtremum { parameter: e.parameter, point: e.point, square_distance: e.square_distance, is_min: e.is_min }
    }
}

/// `|F|` below this is zero: a stationary range end, or a stretch.
const STATIONARY: f64 = 1e-7;

/// Closed-form solutions closer than this are one.
const SAME: f64 = 1e-7;

/// Windows of a search before giving up on growing them (each doubles the last).
const MAX_WINDOWS: usize = 80;

/// The local extremum of `|C(u) − p|²` on `[umin, usup]` whose parameter is nearest `u0` (the
/// larger of two equally near), or `None` if the range has none. `tolf` is the caller's
/// accuracy on the distance function; roots are refined to full precision whatever it is.
pub(crate) fn locate_ext_pc(p: Point, curve: &dyn CurveAdaptor, u0: f64, umin: f64, usup: f64, tolf: f64) -> Option<LocatedExtremum> {
    let _ = tolf;
    if !(umin <= usup) || u0.is_nan() {
        return None;
    }
    let tolu = parametric_tolerance(curve, umin, usup);
    let best = match closed_form(p, curve, umin, usup, tolu) {
        Some(list) => nearest(list?, u0)?,
        None => general(p, curve, u0, umin, usup, tolu)?,
    };
    let point = curve.value(best.0);
    Some(LocatedExtremum { parameter: best.0, point, square_distance: point.distance_squared(p), is_min: best.1 })
}

/// The candidate `(u, is_min)` nearest `u0`, the larger of two equally near.
fn nearest(list: Vec<(f64, bool)>, u0: f64) -> Option<(f64, bool)> {
    let mut best: Option<(f64, bool)> = None;
    for c in list {
        best = match best {
            None => Some(c),
            Some(b) => {
                let (gb, gc) = ((b.0 - u0).abs(), (c.0 - u0).abs());
                let tie = (gb - gc).abs() <= 4.0 * f64::EPSILON * (u0.abs() + b.0.abs().max(c.0.abs()));
                if (tie && c.0 > b.0) || (!tie && gc < gb) { Some(c) } else { Some(b) }
            }
        };
    }
    best
}

/// The general search: windows around the guess (clamped into the range) doubling until one
/// holds a candidate.
fn general(p: Point, curve: &dyn CurveAdaptor, u0: f64, umin: f64, usup: f64, tolu: f64) -> Option<(f64, bool)> {
    let centre = u0.clamp(umin, usup);
    let mut searcher = Searcher::new(curve, p, umin, usup, STATIONARY);
    let mut found: Vec<Candidate> = Vec::new();
    let width = usup - umin;
    let mut radius = (width / 64.0).max(tolu);
    let (mut lo, mut hi) = (centre, centre);
    let mut ends_checked = [false, false];
    for _ in 0..MAX_WINDOWS {
        let (new_lo, new_hi) = ((centre - radius).max(umin), (centre + radius).min(usup));
        if lo == hi && lo == centre {
            searcher.search(new_lo, new_hi, &mut found);
        } else {
            if new_lo < lo {
                searcher.search(new_lo, lo, &mut found);
            }
            if new_hi > hi {
                searcher.search(hi, new_hi, &mut found);
            }
        }
        lo = new_lo;
        hi = new_hi;
        if searcher.failed {
            return None;
        }
        for (k, end) in [umin, usup].into_iter().enumerate() {
            if !ends_checked[k] && (end - lo).abs().min((end - hi).abs()) == 0.0 {
                ends_checked[k] = true;
                if searcher.end_value(end).abs() < STATIONARY {
                    found.push(Candidate { u: end, is_min: end_is_min(&mut searcher, end, k == 0), kind: CandidateKind::Root });
                }
            }
        }
        if !found.is_empty() {
            sort_and_merge(&mut found, 0.0);
            let list = found
                .iter()
                .map(|c| match c.kind {
                    // (A stretch: its point nearest the guess.)
                    CandidateKind::Stretch(a, b) => (u0.clamp(a, b), c.is_min),
                    _ => (c.u, c.is_min),
                })
                .collect();
            return nearest(list, u0);
        }
        if lo <= umin && hi >= usup {
            return None;
        }
        radius *= 2.0;
    }
    None
}

/// Whether a stationary range end is a minimum: `F` there and just inside, signed by the
/// direction into the range (the distance grows into the range).
fn end_is_min(searcher: &mut Searcher, end: f64, at_start: bool) -> bool {
    let f = searcher.end_value(end);
    if at_start { f >= 0.0 } else { f <= 0.0 }
}

/// The candidates of lines and conics, `Some(None)` where no isolated extremum exists, `None`
/// for other curves.
#[allow(clippy::option_option)]
fn closed_form(p: Point, curve: &dyn CurveAdaptor, umin: f64, usup: f64, tolu: f64) -> Option<Option<Vec<(f64, bool)>>> {
    let kind = curve.curve_type()?;
    let basis = curve.basis_curve()?;
    // Stationary points of the whole curve: (parameter, is_min), and whether periodic.
    let (solutions, periodic, slack): (Vec<(f64, bool)>, bool, f64) = match (kind, basis) {
        (CurveKind::Line, Curve::Line(l)) => {
            let t = Vector::from(l.direction).dot(p - l.origin);
            (vec![(t, true)], false, 1e-7)
        }
        (CurveKind::Circle, Curve::Circle(c)) => {
            let (x, y) = local(c.frame(), p);
            if x.hypot(y) < 1e-7 {
                return Some(None);
            }
            let t = y.atan2(x);
            (vec![(t, true), (t + PI, false)], true, tolu)
        }
        (CurveKind::Ellipse, Curve::Ellipse(e)) => {
            let (x, y) = local(e.frame(), p);
            let (a, b) = (e.major_radius(), e.minor_radius());
            if (a - b).abs() < 1e-7 && x.hypot(y) < 1e-7 {
                return Some(None);
            }
            // g(t) = (b² − a²) sin t cos t + a x sin t − b y cos t.
            let g = |t: f64| {
                let (s, c) = t.sin_cos();
                (b * b - a * a) * s * c + a * x * s - b * y * c
            };
            let dg = |t: f64| {
                let (s, c) = t.sin_cos();
                (b * b - a * a) * (c * c - s * s) + a * x * c + b * y * s
            };
            let roots = trigonometric_roots(0.0, 0.5 * (b * b - a * a), -b * y, a * x, 0.0, 0.0, 2.0 * PI);
            if !roots.done || roots.infinite {
                return Some(None);
            }
            let ts: Vec<f64> = roots.roots.iter().map(|&t| polish(&g, &dg, t)).collect();
            let list = classify(ts, &dg, Some(2.0 * PI));
            (list, true, tolu)
        }
        (CurveKind::Parabola, Curve::Parabola(pa)) => {
            let (x, y) = local(pa.frame(), p);
            let f = pa.focal();
            // g(t) = t³/(4f) + (2f − x) t − 2f y.
            let g = |t: f64| t * t * t / (4.0 * f) + (2.0 * f - x) * t - 2.0 * f * y;
            let dg = |t: f64| 3.0 * t * t / (4.0 * f) + 2.0 * f - x;
            let roots = PolynomialRoots::cubic(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
            if !roots.done || roots.infinite {
                return Some(None);
            }
            let ts: Vec<f64> = roots.roots.iter().map(|&t| polish(&g, &dg, t)).collect();
            (classify(ts, &dg, None), false, tolu)
        }
        (CurveKind::Hyperbola, Curve::Hyperbola(h)) => {
            let (x, y) = local(h.frame(), p);
            let (a, b) = (h.major_radius(), h.minor_radius());
            let s2 = a * a + b * b;
            // g(t) = (a² + b²) sinh t cosh t − a x sinh t − b y cosh t; times 4e²ᵗ, a quartic in
            // v = eᵗ: (a² + b²)v⁴ − 2(ax + by)v³ + 2(ax − by)v − (a² + b²).
            let g = |t: f64| s2 * t.sinh() * t.cosh() - a * x * t.sinh() - b * y * t.cosh();
            let dg = |t: f64| s2 * (2.0 * t).cosh() - a * x * t.cosh() - b * y * t.sinh();
            let roots = PolynomialRoots::quartic(s2, -2.0 * (a * x + b * y), 0.0, 2.0 * (a * x - b * y), -s2);
            if !roots.done || roots.infinite {
                return Some(None);
            }
            let ts: Vec<f64> = roots.roots.iter().filter(|&&v| v > 0.0 && v.is_finite()).map(|&v| polish(&g, &dg, v.ln())).collect();
            (classify(ts, &dg, None), false, tolu)
        }
        _ => return None,
    };
    let mut out = Vec::new();
    for (t, is_min) in solutions {
        let t = if periodic || curve.is_periodic() {
            let period = 2.0 * PI;
            let mut u = in_period(t, umin, umin + period);
            if u >= umin + period - 4.0 * f64::EPSILON * (umin.abs() + period) {
                u = umin;
            }
            u
        } else {
            t
        };
        if t >= umin - slack && t <= usup + slack {
            out.push((t, is_min));
        }
    }
    Some(Some(out))
}

/// The coordinates of `p` projected into the plane of `frame` (along its X and Y axes).
fn local(frame: &Frame, p: Point) -> (f64, f64) {
    let w = p - frame.origin();
    (w.dot(Vector::from(frame.x())), w.dot(Vector::from(frame.y())))
}

/// Sorted stationary parameters, those within `SAME` of one kept (modulo `period`) merged,
/// each a minimum where the derivative `dg` of the half-distance's derivative is positive.
fn classify(mut ts: Vec<f64>, dg: &dyn Fn(f64) -> f64, period: Option<f64>) -> Vec<(f64, bool)> {
    ts.sort_by(f64::total_cmp);
    let mut kept: Vec<f64> = Vec::new();
    for t in ts {
        let close = |k: &f64| {
            let d = (t - k).abs();
            d < SAME || period.is_some_and(|p| (d - p).abs() < SAME)
        };
        if !kept.iter().any(close) {
            kept.push(t);
        }
    }
    kept.into_iter().map(|t| (t, dg(t) > 0.0)).collect()
}

/// Newton steps on `g` from `t`, kept while they lower `|g|`.
fn polish(g: &dyn Fn(f64) -> f64, dg: &dyn Fn(f64) -> f64, mut t: f64) -> f64 {
    let mut gt = g(t);
    for _ in 0..8 {
        let d = dg(t);
        if gt == 0.0 || d == 0.0 || !d.is_finite() {
            break;
        }
        let next = t - gt / d;
        let gn = g(next);
        if !(gn.abs() < gt.abs()) {
            break;
        }
        t = next;
        gt = gn;
    }
    t
}

#[cfg(test)]
mod tests;
