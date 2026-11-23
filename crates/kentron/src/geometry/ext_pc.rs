//! The local extrema of the distance from a point to a curve (any curve view), and the
//! search machinery shared with [`locate_ext_pc`](super::locate_ext_pc).
//!
//! The squared distance `d(u) = |C(u) − P|²` has `d′ = 2|C′|·F` with
//! `F(u) = (C(u) − P)·T(u)`, `T` the unit tangent. Its extrema are
//! - the sign changes of `F` inside a smooth span (a curve passing through `P`, and a cusp
//!   where the tangent reverses, are sign changes too);
//! - corners at breakpoints of the curve's continuity where the one-sided limits of `F` have
//!   opposite signs;
//! - stretches where `F` vanishes (the curve keeps its distance from `P`).
//!
//! Each smooth span is sampled adaptively: a sample interval is halved while the tangent or
//! the direction from `P` turns by more than about twenty degrees across it, so that `F`
//! cannot hide many oscillations between samples. Sign changes between samples are refined
//! by Brent's method to full precision. Between samples of one sign, `F` may still dip
//! across zero and back (two extrema close together); dips are recognised from the samples'
//! derivatives (both ends heading towards zero) and from sampled local minima of `|F|`, the
//! turning point of `F` is found as a root of `F′`, and if `F` has changed sign there, both
//! roots are refined.
//!
//! Only `value`, `d1`, `d2` and `d3` of the curve are evaluated; breakpoints are approached
//! from inside each span (at a relative `1e-12` from them), since a view evaluates a knot on
//! one side only.

use crate::geometry::adaptor::{in_period, CurveAdaptor};
use crate::geometry::brent_minimum::BrentMinimum;
use crate::geometry::function_roots::{bracketed_root, Function};
use crate::geometry::Continuity;
use crate::math::{Point, Vector};

/// One local extremum of the distance from a point to a curve.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ExtPCPoint {
    pub parameter: f64,
    pub point: Point,
    pub square_distance: f64,
    pub is_min: bool,
}

/// The largest turn (radians) of the tangent, or of the direction from the point, across a
/// sample interval.
const MAX_TURN: f64 = 0.35;

/// Sample intervals of a span at the start.
const INITIAL_SAMPLES: usize = 8;

/// Halvings of an initial sample interval at most.
const MAX_DEPTH: u32 = 16;

/// Samples of one span at most.
const MAX_SAMPLES: usize = 4096;

/// What kind of extremum a candidate is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum CandidateKind {
    /// A sign change of `F` inside a span.
    Root,
    /// A sign change of `F` across a breakpoint.
    Corner,
    /// `F` vanishes over `[lo, hi]`: every point there is stationary.
    Stretch(f64, f64),
}

/// A candidate extremum of the distance.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Candidate {
    pub u: f64,
    pub is_min: bool,
    pub kind: CandidateKind,
}

/// `F` and what the sampler needs at one parameter.
#[derive(Clone, Copy, Debug)]
struct Sample {
    u: f64,
    f: f64,
    /// `F′`, NaN where the tangent is not defined by the first derivative.
    df: f64,
    tangent: Vector,
    /// The unit direction from the point (zero at the point).
    towards: Vector,
    /// How fast (radians per unit parameter) the tangent and the direction from the point
    /// turn.
    rate: f64,
}

/// The search for extrema of the distance from `p` to `curve` over `[lo, hi]`.
pub(crate) struct Searcher<'a> {
    curve: &'a dyn CurveAdaptor,
    p: Point,
    lo: f64,
    hi: f64,
    /// `|F|` up to this over a whole span is a stretch.
    stretch: f64,
    /// Breakpoints strictly inside `(lo, hi)`, increasing.
    breaks: Vec<f64>,
    /// An evaluation failed (non-finite values).
    pub failed: bool,
}

/// A closure as a [`Function`].
struct Closure<F: FnMut(f64) -> f64>(F);

impl<F: FnMut(f64) -> f64> Function for Closure<F> {
    fn value(&mut self, x: f64, value: &mut f64) -> bool {
        *value = (self.0)(x);
        value.is_finite()
    }
}

/// The offset from a breakpoint at which its sides are evaluated.
pub(crate) fn side_offset(b: f64) -> f64 {
    (1e-12 * b.abs().max(1.0)).max(16.0 * f64::EPSILON * b.abs())
}

impl<'a> Searcher<'a> {
    /// A search over `[lo, hi]`, with `|F| ≤ stretch` over a whole span counting as a stretch.
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, p: Point, lo: f64, hi: f64, stretch: f64) -> Searcher<'a> {
        let list = curve.intervals(Continuity::C2);
        let mut breaks: Vec<f64> = Vec::new();
        if list.len() > 2 {
            let inner = &list[1..list.len() - 1];
            if curve.is_periodic() && curve.period() > 0.0 {
                // Every translate of the curve's breakpoints by whole periods.
                let period = curve.period();
                let (f0, _) = (list[0], list[list.len() - 1]);
                for &k in inner {
                    let base = in_period(k, f0, f0 + period);
                    let j0 = ((lo - base) / period).floor() as i64 - 1;
                    let j1 = ((hi - base) / period).ceil() as i64 + 1;
                    for j in j0..=j1.min(j0 + 100_000) {
                        breaks.push(base + j as f64 * period);
                    }
                }
            } else {
                breaks.extend_from_slice(inner);
            }
        }
        breaks.retain(|&b| b > lo + side_offset(b) && b < hi - side_offset(b));
        breaks.sort_by(f64::total_cmp);
        breaks.dedup_by(|b, a| (*b - *a).abs() <= side_offset(*a));
        Searcher { curve, p, lo, hi, stretch, breaks, failed: false }
    }

    /// Whether `u` is one of the breakpoints inside the range.
    fn is_break(&self, u: f64) -> bool {
        let i = self.breaks.partition_point(|&b| b < u - side_offset(u));
        i < self.breaks.len() && (self.breaks[i] - u).abs() <= side_offset(u)
    }

    /// `F` and its derivative at `u`; `None` on non-finite values.
    fn sample(&mut self, u: f64) -> Option<Sample> {
        let (c, d1, d2) = self.curve.d2(u);
        let w = c - self.p;
        let s = d1.length();
        let wl = w.length();
        let towards = if wl > 0.0 { w / wl } else { Vector::ZERO };
        let rate = if s > 0.0 { d1.cross(d2).length() / (s * s) } else { 0.0 } + if wl > 0.0 { towards.cross(d1).length() / wl } else { 0.0 };
        let out = if s > 0.0 {
            let t = d1 / s;
            let f = w.dot(t);
            let df = s + w.dot(d2) / s - f * d1.dot(d2) / (s * s);
            Sample { u, f, df, tangent: t, towards, rate }
        } else {
            // (The tangent from the first derivative that does not vanish.)
            let mut t = d2;
            if t.length() == 0.0 {
                t = self.curve.d3(u).3;
            }
            let tl = t.length();
            let t = if tl > 0.0 { t / tl } else { Vector::ZERO };
            Sample { u, f: w.dot(t), df: f64::NAN, tangent: t, towards, rate }
        };
        let finite = out.f.is_finite() && c.x.is_finite() && c.y.is_finite() && c.z.is_finite();
        if !finite {
            self.failed = true;
            return None;
        }
        Some(out)
    }

    fn f(&mut self, u: f64) -> f64 {
        self.sample(u).map_or(f64::NAN, |s| s.f)
    }

    /// `F` at the parameter where a span's end `u` is evaluated: just inside the span when `u`
    /// is a breakpoint (`inward` the direction into the span).
    fn inside(&self, u: f64, inward: f64) -> f64 {
        if self.is_break(u) { u + inward * side_offset(u) } else { u }
    }

    /// The candidates in `[a, b] ⊂ [lo, hi]`: roots inside the spans, stretches, and corners
    /// at breakpoints in `[a, b]`. Unordered; duplicates are possible near span ends.
    pub(crate) fn search(&mut self, a: f64, b: f64, out: &mut Vec<Candidate>) {
        if !(b > a) {
            return;
        }
        let i0 = self.breaks.partition_point(|&k| k < a);
        let i1 = self.breaks.partition_point(|&k| k <= b);
        let inner: Vec<f64> = self.breaks[i0..i1].iter().copied().filter(|&k| k > a && k < b).collect();
        let mut stops = vec![a];
        stops.extend_from_slice(&inner);
        stops.push(b);
        for w in stops.windows(2) {
            let (lo, hi) = (self.inside(w[0], 1.0), self.inside(w[1], -1.0));
            if hi > lo {
                self.span(lo, hi, out);
            }
        }
        for &k in &self.breaks[i0..i1].to_vec() {
            self.corner(k, out);
        }
    }

    /// A corner at the breakpoint `k`, if the sides' `F` have opposite signs.
    fn corner(&mut self, k: f64, out: &mut Vec<Candidate>) {
        let delta = side_offset(k);
        let (fl, fr) = (self.f(k - delta), self.f(k + delta));
        if !fl.is_finite() || !fr.is_finite() {
            return;
        }
        let is_min = if fl < 0.0 && fr > 0.0 || fl == 0.0 && fr > 0.0 || fl < 0.0 && fr == 0.0 {
            true
        } else if fl > 0.0 && fr < 0.0 || fl == 0.0 && fr < 0.0 || fl > 0.0 && fr == 0.0 {
            false
        } else {
            return;
        };
        out.push(Candidate { u: k, is_min, kind: CandidateKind::Corner });
    }

    /// Samples of `[lo, hi]` (inside one smooth span), refined where the tangent or the
    /// direction from the point turns too much.
    fn samples(&mut self, lo: f64, hi: f64) -> Option<Vec<Sample>> {
        let n = INITIAL_SAMPLES;
        let mut first = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let u = if i == n { hi } else { lo + (hi - lo) * (i as f64 / n as f64) };
            first.push(self.sample(u)?);
        }
        let min_width = (hi - lo) / (n as f64 * f64::from(1u32 << MAX_DEPTH));
        let mut out = Vec::with_capacity(4 * n);
        out.push(first[0]);
        for w in first.windows(2) {
            self.refine(w[0], w[1], min_width, &mut out)?;
            out.push(w[1]);
        }
        Some(out)
    }

    /// Appends the samples strictly between `a` and `b` that the turning criterion wants: the
    /// midpoint always, and the halves refined further while the tangent or the direction from
    /// the point turns by more than `MAX_TURN` across them (from the three points, or from the
    /// rates of turning at them).
    fn refine(&mut self, a: Sample, b: Sample, min_width: f64, out: &mut Vec<Sample>) -> Option<()> {
        let turn = |x: Vector, y: Vector| -> f64 {
            if x == Vector::ZERO || y == Vector::ZERO {
                0.0
            } else {
                x.cross(y).length().atan2(x.dot(y))
            }
        };
        let width = b.u - a.u;
        if width <= min_width || out.len() >= MAX_SAMPLES {
            return Some(());
        }
        let m = self.sample(0.5 * (a.u + b.u))?;
        let smooth = |p: &Sample, q: &Sample| turn(p.tangent, q.tangent) <= MAX_TURN && turn(p.towards, q.towards) <= MAX_TURN;
        let slow = a.rate.max(m.rate).max(b.rate) * width <= 2.0 * MAX_TURN;
        if smooth(&a, &m) && smooth(&m, &b) && smooth(&a, &b) && slow {
            out.push(m);
            return Some(());
        }
        self.refine(a, m, min_width, out)?;
        out.push(m);
        self.refine(m, b, min_width, out)
    }

    /// A root of `F` bracketed by `[a, b]`, to full precision.
    fn root(&mut self, a: f64, b: f64) -> Option<f64> {
        let mut g = Closure(|x: f64| self.f(x));
        let r = bracketed_root(&mut g, a, b, 0.0, 400, f64::EPSILON);
        if r.value.is_finite() { Some(r.root) } else { None }
    }

    /// The candidates of the span `[lo, hi]` (already moved inside off breakpoints).
    fn span(&mut self, lo: f64, hi: f64, out: &mut Vec<Candidate>) {
        let Some(s) = self.samples(lo, hi) else { return };
        let n = s.len() - 1;
        if s.iter().all(|x| x.f.abs() <= self.stretch) {
            out.push(Candidate { u: 0.5 * (lo + hi), is_min: true, kind: CandidateKind::Stretch(lo, hi) });
            return;
        }
        let root = |is_min: bool, u: f64| Candidate { u, is_min, kind: CandidateKind::Root };
        // Exact zeros at samples between samples of opposite signs.
        for i in 1..n {
            if s[i].f == 0.0 {
                let (l, r) = (s[i - 1].f, s[i + 1].f);
                if l < 0.0 && r > 0.0 {
                    out.push(root(true, s[i].u));
                } else if l > 0.0 && r < 0.0 {
                    out.push(root(false, s[i].u));
                }
            }
        }
        let mut dips: Vec<(usize, usize)> = Vec::new();
        for i in 0..n {
            let (f0, f1) = (s[i].f, s[i + 1].f);
            if f0 == 0.0 || f1 == 0.0 {
                continue;
            }
            if (f0 < 0.0) != (f1 < 0.0) {
                if let Some(u) = self.root(s[i].u, s[i + 1].u) {
                    out.push(root(f0 < 0.0, u));
                }
                continue;
            }
            let sign = f0.signum();
            if sign * s[i].df < 0.0 && sign * s[i + 1].df > 0.0 {
                dips.push((i, i + 1));
            }
        }
        // Samples nearest zero among their neighbours (all of one sign).
        for i in 1..n {
            let f = s[i].f;
            if f == 0.0 {
                continue;
            }
            let sign = f.signum();
            let higher = |j: usize| s[j].f != 0.0 && s[j].f.signum() == sign && sign * s[j].f >= sign * f;
            if higher(i - 1) && higher(i + 1) {
                dips.retain(|&(a, b)| !(a >= i - 1 && b <= i + 1));
                dips.push((i - 1, i + 1));
            }
        }
        for (i, j) in dips {
            let (a, b) = (s[i].u, s[j].u);
            let sign = s[i].f.signum();
            let start = if j == i + 2 { s[i + 1].u } else { 0.5 * (a + b) };
            // The least `sign·F` in the bracket, stopping at the first value across zero.
            let mut across: Option<f64> = None;
            let mut minimiser = BrentMinimum::new(1e-10, 100, 1e-15);
            let mut objective = |x: f64| -> Option<f64> {
                let v = sign * self.f(x);
                if !v.is_finite() {
                    return None;
                }
                if v < 0.0 {
                    across = Some(x);
                    return None;
                }
                Some(v)
            };
            minimiser.perform(&mut objective, a, start, b);
            if let Some(m) = across {
                // `F` crosses zero and comes back: two extrema.
                if let Some(u) = self.root(a, m) {
                    out.push(root(sign < 0.0, u));
                }
                if let Some(u) = self.root(m, b) {
                    out.push(root(sign > 0.0, u));
                }
            }
        }
    }

    /// `F` at a range end, evaluated inside the range if it is a breakpoint of the curve.
    pub(crate) fn end_value(&mut self, u: f64) -> f64 {
        let inward = if u <= self.lo { 1.0 } else { -1.0 };
        let at = if self.is_break(u) { u + inward * side_offset(u) } else { u };
        self.f(at)
    }
}

/// Sorts candidates by parameter and drops any within `tolu` of the one kept before it.
pub(crate) fn sort_and_merge(candidates: &mut Vec<Candidate>, tolu: f64) {
    candidates.sort_by(|a, b| a.u.total_cmp(&b.u));
    let mut kept: Vec<Candidate> = Vec::with_capacity(candidates.len());
    for c in candidates.drain(..) {
        if let Some(last) = kept.last() {
            if c.u - last.u <= tolu {
                continue;
            }
        }
        kept.push(c);
    }
    *candidates = kept;
}

/// The parametric resolution the extrema are located to.
pub(crate) fn parametric_tolerance(curve: &dyn CurveAdaptor, lo: f64, hi: f64) -> f64 {
    let r = curve.resolution(1e-7);
    if r.is_finite() && r > 0.0 { r } else { 1e-9 * (hi - lo).abs().max(1.0) }
}

/// Every local extremum of the squared distance from `p` to `curve` on `[uinf, usup]`,
/// ascending: stationary points inside the smooth spans, corners at breakpoints, one point of
/// each stretch where `|F| ≤ tolf`, and range ends where the curve passes through `p`.
/// `None` when no search can be made: a range within the parametric resolution (and no end at
/// `p`), a curve equidistant from `p` over the whole range, or non-finite evaluations.
/// Parameters of periodic curves are brought into `[uinf, uinf + period)`.
pub(crate) fn ext_pc_other(curve: &dyn CurveAdaptor, p: Point, uinf: f64, usup: f64, tolf: f64) -> Option<Vec<ExtPCPoint>> {
    let tolu = parametric_tolerance(curve, uinf, usup);
    let record = |u: f64, is_min: bool| {
        let q = curve.value(u);
        ExtPCPoint { parameter: u, point: q, square_distance: q.distance_squared(p), is_min }
    };
    let infinite = |u: f64| u.abs() >= 1e100;
    let zero_end = |u: f64| !infinite(u) && curve.value(u).distance_squared(p) < 1e-14;
    if !(usup - uinf > tolu) {
        let mut out = Vec::new();
        if zero_end(uinf) {
            out.push(record(uinf, true));
        } else if zero_end(usup) {
            out.push(record(usup, true));
        }
        return if out.is_empty() { None } else { Some(out) };
    }
    let mut searcher = Searcher::new(curve, p, uinf, usup, tolf);
    let mut candidates = Vec::new();
    searcher.search(uinf, usup, &mut candidates);
    if searcher.failed {
        return None;
    }
    if !candidates.is_empty() && candidates.iter().all(|c| matches!(c.kind, CandidateKind::Stretch(..))) {
        let covered: f64 = candidates.iter().map(|c| if let CandidateKind::Stretch(a, b) = c.kind { b - a } else { 0.0 }).sum();
        // (Equidistant over the whole range, up to the sides of breakpoints.)
        if covered >= (usup - uinf) * (1.0 - 1e-9) - 2.0 * side_offset(usup) * (searcher.breaks.len() as f64 + 1.0) {
            return None;
        }
    }
    let periodic = curve.is_periodic() && curve.period() > 0.0;
    let mut points: Vec<(f64, bool)> = Vec::new();
    for c in &candidates {
        let u = if periodic { in_period(c.u, uinf, uinf + curve.period()) } else { c.u };
        if u >= uinf - tolu && u <= usup + tolu {
            points.push((u, c.is_min));
        }
    }
    for end in [uinf, usup] {
        if zero_end(end) && !points.iter().any(|&(u, _)| (u - end).abs() <= tolu) {
            points.push((end, true));
        }
    }
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out: Vec<ExtPCPoint> = Vec::with_capacity(points.len());
    for (u, is_min) in points {
        if out.last().is_some_and(|l: &ExtPCPoint| u - l.parameter <= tolu) {
            continue;
        }
        out.push(record(u, is_min));
    }
    Some(out)
}

#[cfg(test)]
pub(crate) mod tests;
