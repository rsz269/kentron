//! A face's parameter-space curve (pcurve) from a 3D curve lying on or near its surface.
//!
//! The curve's points are taken to their feet on the surface (the nearest surface points,
//! [`surface_analysis::SurfaceInverse`]), the first by a search of the whole surface and
//! each next one by Newton steps from where the last one predicts it. The feet's parameters
//! are unwrapped across seams (neighbouring feet never more than half a period apart),
//! placed so that the start lies in the bounds' first period, and made definite where the
//! surface collapses (a pole's longitude follows the curve's neighbouring points, or at an
//! end the direction it arrives from).
//!
//! The result is a line where the feet move affinely in the curve's parameter (isos, circles
//! about an axis, helices), with a boundary iso's fixed coordinate exactly the bound;
//! otherwise a piecewise cubic Hermite interpolant of the feet and of their exact
//! derivatives (from the linearised foot condition), its spans halved where it strays from
//! the feet in 3D, written as a C1 B-spline (C0 at the curve's own corners).

use crate::math::{Direction2, Point, Point2, Vector};
use crate::tolerance::LINEAR;

use super::adaptor::in_period;
use super::bspline::BSplineCurve;
use super::curve::Curve;
use super::curve2::{Curve2, Line2};
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;
use super::surface_analysis::{period_shift, SingularPlace, SurfaceInverse};

/// The most feet a pcurve is made of.
const MAX_NODES: usize = 400;

/// Intervals the curve is cut into at the start, at least.
const MIN_INTERVALS: usize = 8;

/// The pcurve builder for one surface.
pub(crate) struct PCurveProjector<'a> {
    surface: &'a Surface,
    inverse: SurfaceInverse<'a>,
    preci: f64,
}

/// A foot of the curve: its parameter, surface parameters and the derivatives of the
/// surface parameters along the curve on either side (`None` where unknown).
#[derive(Clone, Copy, Debug)]
struct Node {
    t: f64,
    uv: [f64; 2],
    left: Option<[f64; 2]>,
    right: Option<[f64; 2]>,
    /// On a singular place of the surface (its free coordinate made up).
    on_place: bool,
}

/// How a direction of the parameter space repeats: its period, and whether it is only a
/// geometric closure (the surface not defined beyond its bounds).
#[derive(Clone, Copy, Debug)]
struct Repeat {
    period: f64,
    pseudo: bool,
}

/// What is known of the surface while one curve is projected.
struct Setting {
    tol: f64,
    shape_tol: f64,
    bounds: (f64, f64, f64, f64),
    repeat: [Option<Repeat>; 2],
    /// The parametric resolutions of `tol`.
    res: [f64; 2],
    places: Vec<SingularPlace>,
    /// Coordinates fixed exactly on a boundary iso.
    snapped: [Option<f64>; 2],
}

impl Setting {
    /// `uv` with the coordinates along a closure that is no period brought into the bounds
    /// (where the surface is defined).
    fn wrap(&self, uv: [f64; 2]) -> [f64; 2] {
        let (u0, u1, v0, v1) = self.bounds;
        let mut out = uv;
        for (k, (lo, hi)) in [(u0, u1), (v0, v1)].into_iter().enumerate() {
            if let Some(r) = self.repeat[k] {
                if r.pseudo && (out[k] < lo || out[k] > hi) {
                    out[k] = in_period(out[k], lo, hi);
                }
            }
        }
        out
    }
}

impl<'a> PCurveProjector<'a> {
    pub(crate) fn new(surface: &'a Surface, preci: f64) -> Self {
        PCurveProjector { surface, inverse: SurfaceInverse::new(surface), preci }
    }

    /// The pcurve of `curve` over `[first, last]` on the surface, parametrised like the
    /// curve; `tol_first` and `tol_last` are the end vertices' tolerances (negative: the
    /// linear tolerance), used to tell an end on a pole. `None` only for a degenerate range or
    /// a curve or surface that does not evaluate.
    pub(crate) fn perform(&mut self, curve: &Curve, first: f64, last: f64, tol_first: f64, tol_last: f64) -> Option<Curve2> {
        let scale = 1.0 + first.abs().max(last.abs());
        if !(first.is_finite() && last.is_finite()) || !(last - first > 4.0 * f64::EPSILON * scale) {
            return None;
        }
        let (start, end) = (curve.point(first), curve.point(last));
        if !finite(start) || !finite(end) {
            return None;
        }
        if let Some(plane) = plane_basis(self.surface) {
            if let Some(pc) = super::project::analytic_pcurve(curve, &Surface::Plane(plane)) {
                return Some(pc);
            }
        }
        let tol = self.preci.max(LINEAR);
        let shape_tol = if self.preci > 1.0 { LINEAR } else { tol };
        let bounds = self.inverse.bounds();
        let s = self.surface;
        let repeat_of = |periodic: Option<f64>, closed: bool, lo: f64, hi: f64| match periodic {
            Some(p) if p > 0.0 => Some(Repeat { period: p, pseudo: false }),
            _ if closed && hi - lo > 0.0 && (hi - lo).is_finite() => Some(Repeat { period: hi - lo, pseudo: true }),
            _ => None,
        };
        let u_periodic = s.u_period().filter(|_| s.is_u_periodic());
        let v_periodic = s.v_period().filter(|_| s.is_v_periodic());
        let u_closed = u_periodic.is_none() && self.inverse.closed_in_u(tol);
        let v_closed = v_periodic.is_none() && self.inverse.closed_in_v(tol);
        let repeat = [repeat_of(u_periodic, u_closed, bounds.0, bounds.1), repeat_of(v_periodic, v_closed, bounds.2, bounds.3)];
        let resolution = |r: f64| if r.is_finite() && r > 0.0 { r } else { 1e-12 };
        let res = [
            resolution(super::resolution::u_resolution(s, bounds, tol)),
            resolution(super::resolution::v_resolution(s, bounds, tol)),
        ];
        let places: Vec<SingularPlace> = self.inverse.singular_places().iter().copied().filter(|p| p.spread <= tol).collect();
        let mut setting = Setting { tol, shape_tol, bounds, repeat, res, places, snapped: [None, None] };
        let tol_first = if tol_first < 0.0 { LINEAR } else { tol_first };
        let tol_last = if tol_last < 0.0 { LINEAR } else { tol_last };

        let mut nodes = self.first_nodes(curve, first, last, &setting)?;
        self.settle_places(curve, &mut nodes, &setting, tol_first, tol_last);
        place_periods(&mut nodes, &setting);
        self.snap_to_bounds(curve, &mut nodes, &mut setting);
        if let Some(line) = self.straight(curve, &nodes, &setting) {
            return Some(line);
        }
        self.refine(curve, &mut nodes, &setting);
        Some(Curve2::BSpline(assemble(&nodes, self.surface)))
    }

    /// The feet at the starting parameters: the curve's breaks and an even spread between,
    /// each foot found from the one before.
    fn first_nodes(&mut self, curve: &Curve, first: f64, last: f64, setting: &Setting) -> Option<Vec<Node>> {
        let breaks = curve_breaks(curve, first, last);
        let mut ends: Vec<(f64, bool)> = vec![(first, false)];
        ends.extend(breaks.iter().copied());
        ends.push((last, false));
        let spans = ends.len() - 1;
        let per_span = MIN_INTERVALS.div_ceil(spans).max(2).min((MAX_NODES / 2 / spans).max(1));
        let mut params: Vec<(f64, bool)> = Vec::with_capacity(spans * per_span + 1);
        for (k, w) in ends.windows(2).enumerate() {
            if k == 0 {
                params.push(w[0]);
            }
            for i in 1..per_span {
                params.push((w[0].0 + (w[1].0 - w[0].0) * i as f64 / per_span as f64, false));
            }
            params.push(w[1]);
        }
        let mut nodes: Vec<Node> = Vec::with_capacity(params.len());
        for &(t, corner) in &params {
            let p = curve.point(t);
            if !finite(p) {
                return None;
            }
            let (uv, _gap) = match nodes.last() {
                None => {
                    let (uv, gap) = self.inverse.nearest(p, self.preci);
                    ([uv.x, uv.y], gap)
                }
                Some(prev) => {
                    let predicted = if self.inverse.is_elementary() { prev.uv } else { predict(prev, t, setting) };
                    self.foot_near(p, predicted, setting)
                }
            };
            if !(uv[0].is_finite() && uv[1].is_finite()) {
                return None;
            }
            let mut node = Node { t, uv, left: None, right: None, on_place: false };
            if self.place_of(setting, p, uv).is_some() {
                // (Its free coordinate is made up later; no derivative there.)
                node.on_place = true;
            } else {
                self.derivatives(curve, &mut node, corner, first, last, setting);
            }
            nodes.push(node);
        }
        Some(nodes)
    }

    /// The foot of `p` searched near `hint` (unwrapped parameters), unwrapped next to it.
    fn foot_near(&mut self, p: Point, hint: [f64; 2], setting: &Setting) -> ([f64; 2], f64) {
        let inside = setting.wrap(hint);
        let max_gap = 10.0 * setting.tol;
        let (uv, gap) = self.inverse.nearest_from(Point2::new(inside[0], inside[1]), p, self.preci, max_gap);
        let mut uv = [uv.x, uv.y];
        for k in 0..2 {
            if let Some(r) = setting.repeat[k] {
                uv[k] += period_shift(uv[k], hint[k], r.period);
            }
            if let Some(b) = setting.snapped[k] {
                uv[k] = b;
            }
        }
        (uv, gap)
    }

    /// The derivatives of the node's surface parameters along the curve (both sides at a
    /// corner of the curve).
    fn derivatives(&self, curve: &Curve, node: &mut Node, corner: bool, first: f64, last: f64, setting: &Setting) {
        let surface = self.surface;
        let uv = setting.wrap(node.uv);
        let at = |t: f64| {
            let (p, d) = curve.d1(t);
            foot_derivative(surface, uv, p, d)
        };
        if corner {
            let delta = 1e-9 * (last - first).max(node.t.abs().max(1.0) * 1e-3);
            node.left = at(node.t - delta);
            node.right = at(node.t + delta);
        } else {
            let d = at(node.t);
            node.left = d;
            node.right = d;
        }
    }

    /// Feet on the surface's singular places: the free coordinate taken from the nearest
    /// foot off the place, or at an end on a true pole from the direction the curve arrives
    /// from; a curve wholly on a place gets its fixed coordinate exactly and the free one
    /// running linearly between its ends'.
    fn settle_places(&mut self, curve: &Curve, nodes: &mut [Node], setting: &Setting, tol_first: f64, tol_last: f64) {
        if setting.places.is_empty() {
            return;
        }
        let n = nodes.len();
        let mut which: Vec<Option<usize>> = vec![None; n];
        for (i, node) in nodes.iter_mut().enumerate() {
            let p = curve.point(node.t);
            which[i] = self.place_of(setting, p, node.uv);
            if which[i].is_some() {
                node.on_place = true;
                node.left = None;
                node.right = None;
            }
        }
        if which.iter().all(Option::is_some) {
            // The whole curve on one place (a degenerate edge).
            let place = setting.places[which[0].expect("place")];
            let (fixed, free) = if place.fixed_is_u { (0, 1) } else { (1, 0) };
            let (a, b) = (nodes[0].uv[free], nodes[n - 1].uv[free]);
            let (t0, t1) = (nodes[0].t, nodes[n - 1].t);
            for node in nodes.iter_mut() {
                node.uv[fixed] = place.fixed;
                node.uv[free] = a + (b - a) * (node.t - t0) / (t1 - t0);
                node.left = Some(derivative_of(fixed, free, (b - a) / (t1 - t0)));
                node.right = node.left;
            }
            return;
        }
        for i in 0..n {
            let Some(k) = which[i] else { continue };
            let place = setting.places[k];
            let free = if place.fixed_is_u { 1 } else { 0 };
            // The nearest node off the place.
            let before = (0..i).rev().find(|&j| which[j].is_none());
            let after = (i + 1..n).find(|&j| which[j].is_none());
            let source = match (before, after) {
                (Some(b), Some(a)) => {
                    if nodes[i].t - nodes[b].t <= nodes[a].t - nodes[i].t { b } else { a }
                }
                (Some(b), None) => b,
                (None, Some(a)) => a,
                (None, None) => continue,
            };
            let mut value = nodes[source].uv[free];
            let at_end = i == 0 || i == n - 1;
            let end_tol = if i == 0 { tol_first } else { tol_last };
            if at_end && place.spread <= LINEAR && curve.point(nodes[i].t).distance(place.point) <= end_tol.max(LINEAR) {
                if let Some(w) = self.arrival(curve, &place, nodes[i].t, i == 0, nodes[source]) {
                    value = w;
                }
            }
            if let Some(r) = setting.repeat[free] {
                value += period_shift(value, nodes[source].uv[free], r.period);
            }
            nodes[i].uv[free] = value;
        }
    }

    /// The free coordinate on a pole that the curve reaches at its end `t` (its start if
    /// `at_start`): the one whose iso leaves the pole in the direction the curve comes from;
    /// `None` where that direction is not clear.
    fn arrival(&self, curve: &Curve, place: &SingularPlace, t: f64, at_start: bool, neighbour: Node) -> Option<f64> {
        let s = self.surface;
        let (_, d) = curve.d1(t);
        let mut dir = if at_start { d } else { -d };
        if !(dir.length() > 1e-12 * (1.0 + (curve.point(t) - Point::ORIGIN).length())) {
            dir = curve.point(neighbour.t) - place.point;
        }
        let len = dir.length();
        if !(len > 0.0) || !len.is_finite() {
            return None;
        }
        let dir = dir * (1.0 / len);
        let (fixed_k, free_k) = if place.fixed_is_u { (0, 1) } else { (1, 0) };
        let inward = if neighbour.uv[fixed_k] >= place.fixed { 1.0 } else { -1.0 };
        // The cosine of the angle between the iso leaving the pole at `w` and the direction,
        // and its derivative in `w`.
        let leaving = |w: f64| -> (f64, f64) {
            let (u, v) = if place.fixed_is_u { (place.fixed, w) } else { (w, place.fixed) };
            let d2 = s.d2(u, v);
            let e = if place.fixed_is_u { d2.du } else { d2.dv } * inward;
            let de = d2.duv * inward;
            let l = e.length();
            if !(l > 0.0) {
                return (-2.0, 0.0);
            }
            let c = e.dot(dir) / l;
            (c, de.dot(dir) / l - c * e.dot(de) / (l * l))
        };
        // The free range: a period about the neighbour's value where it repeats.
        let period = match free_k {
            0 => s.u_period().filter(|_| s.is_u_periodic()),
            _ => s.v_period().filter(|_| s.is_v_periodic()),
        };
        let (lo, hi) = match period {
            Some(p) => (neighbour.uv[free_k] - 0.5 * p, neighbour.uv[free_k] + 0.5 * p),
            None => place.free_range,
        };
        if !(hi > lo) || !(hi - lo).is_finite() {
            return None;
        }
        const STEPS: usize = 72;
        let step = (hi - lo) / STEPS as f64;
        let mut best = (f64::NEG_INFINITY, lo);
        for i in 0..=STEPS {
            let w = lo + step * i as f64;
            let (c, _) = leaving(w);
            if c > best.0 {
                best = (c, w);
            }
        }
        // The maximum, where the derivative changes sign, by bisection about the best sample.
        let (mut a, mut b) = (best.1 - step, best.1 + step);
        if period.is_none() {
            a = a.max(lo);
            b = b.min(hi);
        }
        let mut w = best.1;
        if leaving(a).1 > 0.0 && leaving(b).1 < 0.0 {
            for _ in 0..100 {
                let m = 0.5 * (a + b);
                if m <= a || m >= b {
                    break;
                }
                if leaving(m).1 > 0.0 {
                    a = m;
                } else {
                    b = m;
                }
            }
            w = 0.5 * (a + b);
        }
        (leaving(w).0 > 0.9).then_some(w)
    }

    /// Fixes a coordinate at a bound exactly where every foot lies on that boundary iso
    /// (and it is no iso the surface collapses).
    fn snap_to_bounds(&self, curve: &Curve, nodes: &mut [Node], setting: &mut Setting) {
        let (u0, u1, v0, v1) = setting.bounds;
        let s = self.surface;
        let near = [
            resolution_or(super::resolution::u_resolution(s, setting.bounds, LINEAR)),
            resolution_or(super::resolution::v_resolution(s, setting.bounds, LINEAR)),
        ];
        for (k, bound) in [(0, u0), (0, u1), (1, v0), (1, v1)] {
            if setting.snapped[k].is_some() || !bound.is_finite() || bound.abs() >= 1e100 {
                continue;
            }
            let collapsed = setting.places.iter().any(|p| p.fixed_is_u == (k == 0) && (p.fixed - bound).abs() <= near[k]);
            if collapsed {
                continue;
            }
            let on = nodes.iter().all(|n| {
                if (n.uv[k] - bound).abs() > near[k] {
                    return false;
                }
                // (And in 3D: the curve's point within the linear tolerance of the iso's.)
                let mut uv = n.uv;
                uv[k] = bound;
                s.point(uv[0], uv[1]).distance(curve.point(n.t)) <= s.point(n.uv[0], n.uv[1]).distance(curve.point(n.t)) + LINEAR
            });
            if on {
                for n in nodes.iter_mut() {
                    n.uv[k] = bound;
                    for d in [&mut n.left, &mut n.right] {
                        if let Some(d) = d.as_mut() {
                            d[k] = 0.0;
                        }
                    }
                }
                setting.snapped[k] = Some(bound);
            }
        }
    }

    /// The pcurve as a line (or a degree-1 B-spline) where the feet move affinely in the
    /// curve's parameter: checked at the nodes and between them, in 3D, by the distance of
    /// the curve's point from the surface normal line through the straight pcurve's point.
    fn straight(&self, curve: &Curve, nodes: &[Node], setting: &Setting) -> Option<Curve2> {
        let s = self.surface;
        let n = nodes.len();
        let (a, z) = (nodes[0], nodes[n - 1]);
        let (t0, t1) = (a.t, z.t);
        let b = [(z.uv[0] - a.uv[0]) / (t1 - t0), (z.uv[1] - a.uv[1]) / (t1 - t0)];
        let at = |t: f64| [a.uv[0] + b[0] * (t - t0), a.uv[1] + b[1] * (t - t0)];
        let fits = |t: f64| -> bool {
            let uv = at(t);
            let uv = setting.wrap(uv);
            let (q, du, dv) = s.d1(uv[0], uv[1]);
            let c = curve.point(t);
            let w = c - q;
            let normal = du.cross(dv);
            let nl = normal.length();
            if nl > 1e-9 * du.length() * dv.length() && nl > 1e-300 {
                let n = normal * (1.0 / nl);
                (w - n * w.dot(n)).length() <= setting.shape_tol
            } else {
                // (A degenerate normal, as at a pole: the 3D distance against the foot's.)
                let foot = s.point(uv[0], uv[1]);
                c.distance(foot) <= setting.shape_tol + nodes.iter().map(|n| curve.point(n.t).distance(s.point(n.uv[0], n.uv[1]))).fold(0.0, f64::max)
            }
        };
        if !nodes.iter().all(|node| fits(node.t)) || !nodes.windows(2).all(|w| fits(0.5 * (w[0].t + w[1].t))) {
            return None;
        }
        let speed = b[0].hypot(b[1]);
        if (speed - 1.0).abs() <= 1e-9 {
            let d = Direction2::new(b[0], b[1]);
            let origin = Point2::new(a.uv[0] - t0 * d.x(), a.uv[1] - t0 * d.y());
            let origin = Point2::new(
                setting.snapped[0].unwrap_or(origin.x),
                setting.snapped[1].unwrap_or(origin.y),
            );
            return Some(Curve2::Line(Line2::new(origin, d)));
        }
        let start = Point2::new(a.uv[0], a.uv[1]);
        let mut end = Point2::new(z.uv[0], z.uv[1]);
        if start.distance(end) == 0.0 {
            // (All feet one point: the end moved by a negligible amount, for a non-zero length.)
            let k = if setting.snapped[0].is_some() { 1 } else { 0 };
            let nudge = 1e-3 * setting.res[k].min(1e-6);
            if k == 0 { end.x += nudge } else { end.y += nudge }
        }
        Some(Curve2::BSpline(BSplineCurve::new(vec![start, end], vec![t0, t1], vec![2, 2], 1)))
    }

    /// Cuts in three the spans whose Hermite interpolant strays from the feet by more than
    /// half the tolerance in 3D (checked at their thirds, which become nodes), until all are
    /// within it or the node budget is spent.
    fn refine(&mut self, curve: &Curve, nodes: &mut Vec<Node>, setting: &Setting) {
        let target = 0.5 * setting.tol;
        let (first, last) = (nodes[0].t, nodes[nodes.len() - 1].t);
        let shortest = 1e-10 * (last - first).max(1e-300);
        let mut done: Vec<bool> = vec![false; nodes.len() - 1];
        loop {
            let tangents = tangents(nodes);
            let mut next: Vec<Node> = Vec::with_capacity(nodes.len() * 3);
            let mut next_done: Vec<bool> = Vec::with_capacity(nodes.len() * 3);
            let mut added = 0;
            for i in 0..nodes.len() - 1 {
                next.push(nodes[i]);
                if done[i] {
                    next_done.push(true);
                    continue;
                }
                let (a, b) = (nodes[i], nodes[i + 1]);
                let budget = nodes.len() + added + 2 <= MAX_NODES;
                if b.t - a.t <= shortest || !budget {
                    next_done.push(true);
                    continue;
                }
                let probe = |this: &mut Self, s: f64| {
                    let t = a.t + (b.t - a.t) * s;
                    let guess = hermite(a.t, a.uv, tangents[i].1, b.t, b.uv, tangents[i + 1].0, t);
                    let node = this.node_at(curve, t, guess, setting, first, last);
                    let (g, f) = (setting.wrap(guess), setting.wrap(node.uv));
                    let error = this.surface.point(g[0], g[1]).distance(this.surface.point(f[0], f[1]));
                    (node, error)
                };
                let (n1, e1) = probe(self, 1.0 / 3.0);
                let (n2, e2) = probe(self, 2.0 / 3.0);
                let fine = |e: f64| e.is_nan() || e <= target;
                if fine(e1) && fine(e2) {
                    next_done.push(true);
                    continue;
                }
                next.push(n1);
                next.push(n2);
                next_done.extend([false, false, false]);
                added += 2;
            }
            next.push(nodes[nodes.len() - 1]);
            *nodes = next;
            done = next_done;
            if added == 0 {
                break;
            }
        }
    }

    /// The singular place the curve's point `p` (or its foot `uv`) lies on.
    fn place_of(&self, setting: &Setting, p: Point, uv: [f64; 2]) -> Option<usize> {
        if setting.places.is_empty() {
            return None;
        }
        let w = setting.wrap(uv);
        let foot = self.surface.point(w[0], w[1]);
        setting.places.iter().position(|place| p.distance(place.point) <= setting.tol + place.spread || foot.distance(place.point) <= setting.tol + place.spread)
    }

    /// The node at `t`, its foot searched from `guess` (on a singular place, `guess` itself
    /// with any snapped coordinate).
    fn node_at(&mut self, curve: &Curve, t: f64, guess: [f64; 2], setting: &Setting, first: f64, last: f64) -> Node {
        let p = curve.point(t);
        let (uv, _) = self.foot_near(p, guess, setting);
        if let Some(k) = self.place_of(setting, p, uv) {
            // (The free coordinate made up: the interpolant's.)
            let free = if setting.places[k].fixed_is_u { 1 } else { 0 };
            let mut uv = uv;
            uv[free] = guess[free];
            return Node { t, uv, left: None, right: None, on_place: true };
        }
        let mut node = Node { t, uv, left: None, right: None, on_place: false };
        self.derivatives(curve, &mut node, false, first, last, setting);
        for (k, b) in setting.snapped.iter().enumerate() {
            if b.is_some() {
                for d in [&mut node.left, &mut node.right] {
                    if let Some(d) = d.as_mut() {
                        d[k] = 0.0;
                    }
                }
            }
        }
        node
    }
}

/// The plane a trimmed or offset plane is made of.
fn plane_basis(surface: &Surface) -> Option<super::surface::Plane> {
    match surface {
        Surface::Plane(p) => Some(*p),
        Surface::Trimmed(t) => plane_basis(t.basis()),
        Surface::Offset(o) => {
            // (An offset plane is the plane moved along its normal, parametrised alike.)
            let Surface::Plane(p) = o.basis() else { return None };
            let moved = surface.point(0.0, 0.0) - o.basis().point(0.0, 0.0);
            Some(crate::math::Transformable::transformed(p, &crate::math::Transform::translation(moved)))
        }
        _ => None,
    }
}

fn finite(p: Point) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

fn resolution_or(r: f64) -> f64 {
    if r.is_finite() && r > 0.0 { r } else { 1e-12 }
}

/// The breaks of `curve` strictly inside `(first, last)`, and whether it is only C0 there
/// (its speed may jump).
fn curve_breaks(curve: &Curve, first: f64, last: f64) -> Vec<(f64, bool)> {
    let margin = 1e-9 * (last - first).abs().max(1e-300);
    match curve {
        Curve::BSpline(b) => b
            .knots()
            .iter()
            .zip(b.multiplicities())
            .filter(|(k, _)| **k > first + margin && **k < last - margin)
            .map(|(&k, &m)| (k, m >= b.degree()))
            .collect(),
        Curve::Trimmed(t) => curve_breaks(t.basis(), first, last),
        Curve::Offset(o) => {
            let mut out = curve_breaks(o.basis(), first, last);
            // (An offset loses a derivative: the basis's C1 knots are its corners.)
            if let Curve::BSpline(b) = o.basis() {
                for (t, c0) in out.iter_mut() {
                    let m = b.knots().iter().position(|k| k == t).map_or(0, |i| b.multiplicities()[i]);
                    *c0 = *c0 || m + 1 >= b.degree();
                }
            }
            out
        }
        _ => Vec::new(),
    }
}

/// The derivative of the foot `uv` of the curve's point `c` (moving at `dc`) on the surface:
/// from differentiating `(S − C)·S_u = (S − C)·S_v = 0` along the curve. `None` where that
/// system is singular (a degenerate surface point, or a point at a centre of curvature).
fn foot_derivative(surface: &Surface, uv: [f64; 2], c: Point, dc: Vector) -> Option<[f64; 2]> {
    let d = surface.d2(uv[0], uv[1]);
    let r = d.point - c;
    let (a11, a12, a22) = (d.du.dot(d.du) + r.dot(d.duu), d.du.dot(d.dv) + r.dot(d.duv), d.dv.dot(d.dv) + r.dot(d.dvv));
    let (b1, b2) = (dc.dot(d.du), dc.dot(d.dv));
    let det = a11 * a22 - a12 * a12;
    let (lu, lv) = (d.du.length(), d.dv.length());
    let scale = lu * lu * lv * lv;
    // (Not where the surface is nearly degenerate: a side collapsing to a point.)
    if !(lu.min(lv) > 1e-6 * lu.max(lv)) || !(d.du.cross(d.dv).length() > 1e-6 * lu * lv) {
        return None;
    }
    if !(det.abs() > 1e-10 * scale) || !det.is_finite() {
        return None;
    }
    let out = [(b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det];
    (out[0].is_finite() && out[1].is_finite()).then_some(out)
}

/// A derivative whose fixed coordinate does not move.
fn derivative_of(fixed: usize, free: usize, rate: f64) -> [f64; 2] {
    let mut d = [0.0; 2];
    d[fixed] = 0.0;
    d[free] = rate;
    d
}

/// Where the node's foot is expected at `t`, from its derivative.
fn predict(node: &Node, t: f64, setting: &Setting) -> [f64; 2] {
    let Some(d) = node.right.filter(|_| !node.on_place) else { return node.uv };
    let step = [d[0] * (t - node.t), d[1] * (t - node.t)];
    // (A step across a large part of the bounds is no prediction to trust.)
    let (u0, u1, v0, v1) = setting.bounds;
    let reach = |k: usize, lo: f64, hi: f64| {
        let width = setting.repeat[k].map_or(hi - lo, |r| r.period);
        step[k].abs() <= 0.25 * width || !width.is_finite()
    };
    if reach(0, u0, u1) && reach(1, v0, v1) && step[0].is_finite() && step[1].is_finite() {
        [node.uv[0] + step[0], node.uv[1] + step[1]]
    } else {
        node.uv
    }
}

/// The derivatives on either side of each node, the unknown ones from the neighbouring
/// feet (central differences inside, secants at the ends and on either side of a corner).
fn tangents(nodes: &[Node]) -> Vec<([f64; 2], [f64; 2])> {
    let n = nodes.len();
    let secant = |i: usize, j: usize| {
        let h = nodes[j].t - nodes[i].t;
        [(nodes[j].uv[0] - nodes[i].uv[0]) / h, (nodes[j].uv[1] - nodes[i].uv[1]) / h]
    };
    (0..n)
        .map(|i| {
            let central = || {
                if i == 0 {
                    secant(0, 1)
                } else if i == n - 1 {
                    secant(n - 2, n - 1)
                } else {
                    secant(i - 1, i + 1)
                }
            };
            let left = nodes[i].left.unwrap_or_else(|| if nodes[i].right.is_some() && i > 0 { secant(i - 1, i) } else { central() });
            let right = nodes[i].right.unwrap_or_else(|| if nodes[i].left.is_some() && i + 1 < n { secant(i, i + 1) } else { central() });
            (left, right)
        })
        .collect()
}

/// The cubic Hermite interpolant of `(t0, p0, m0)` and `(t1, p1, m1)` at `t`.
fn hermite(t0: f64, p0: [f64; 2], m0: [f64; 2], t1: f64, p1: [f64; 2], m1: [f64; 2], t: f64) -> [f64; 2] {
    let h = t1 - t0;
    let s = (t - t0) / h;
    let (s2, s3) = (s * s, s * s * s);
    let (h00, h10, h01, h11) = (2.0 * s3 - 3.0 * s2 + 1.0, s3 - 2.0 * s2 + s, -2.0 * s3 + 3.0 * s2, s3 - s2);
    [0, 1].map(|k| h00 * p0[k] + h10 * h * m0[k] + h01 * p1[k] + h11 * h * m1[k])
}

/// Shifts the feet by whole periods so that the first lies in the bounds' first period (on
/// the seam: at its start when the pcurve moves on to greater values, else at its end), and
/// brings feet that a closure (no period) took past the bounds back within them.
fn place_periods(nodes: &mut [Node], setting: &Setting) {
    let (u0, u1, v0, v1) = setting.bounds;
    for k in 0..2 {
        let Some(r) = setting.repeat[k] else { continue };
        let (lo, hi) = if k == 0 { (u0, u1) } else { (v0, v1) };
        if !lo.is_finite() || lo.abs() >= 1e100 {
            continue;
        }
        let res = setting.res[k];
        let x0 = nodes[0].uv[k];
        let mut placed = in_period(x0, lo, lo + r.period);
        let on_seam = placed - lo <= res || lo + r.period - placed <= res;
        if on_seam {
            let decreasing = match nodes.iter().find(|n| (n.uv[k] - x0).abs() > res) {
                Some(n) => n.uv[k] < x0,
                // (Along the seam: its start.)
                None => false,
            };
            let start = if placed - lo <= res { placed } else { placed - r.period };
            placed = if decreasing { start + r.period } else { start };
        }
        let shift = placed - x0;
        if shift != 0.0 {
            for n in nodes.iter_mut() {
                n.uv[k] += shift;
            }
            nodes[0].uv[k] = placed;
        }
        if r.pseudo {
            for n in nodes.iter_mut() {
                if n.uv[k] > hi + res {
                    n.uv[k] -= r.period;
                } else if n.uv[k] < lo - res {
                    n.uv[k] += r.period;
                }
            }
        }
    }
}

/// The nodes as a cubic B-spline: double knots (C1) at the nodes, triple ones (C0) where the
/// derivatives on either side differ; on a sphere, poles kept within its latitudes.
fn assemble(nodes: &[Node], surface: &Surface) -> super::curve2::BSplineCurve2 {
    let tangents = tangents(nodes);
    let n = nodes.len();
    let p = |k: usize| Point2::new(nodes[k].uv[0], nodes[k].uv[1]);
    let mut poles: Vec<Point2> = vec![p(0)];
    let mut knots = vec![nodes[0].t];
    let mut mults = vec![4];
    for i in 0..n - 1 {
        let h = nodes[i + 1].t - nodes[i].t;
        let (m0, m1) = (tangents[i].1, tangents[i + 1].0);
        poles.push(Point2::new(nodes[i].uv[0] + h / 3.0 * m0[0], nodes[i].uv[1] + h / 3.0 * m0[1]));
        poles.push(Point2::new(nodes[i + 1].uv[0] - h / 3.0 * m1[0], nodes[i + 1].uv[1] - h / 3.0 * m1[1]));
        if i + 1 < n - 1 {
            let (l, r) = tangents[i + 1];
            let smooth = (l[0] - r[0]).abs() <= 1e-12 * (1.0 + l[0].abs()) && (l[1] - r[1]).abs() <= 1e-12 * (1.0 + l[1].abs());
            knots.push(nodes[i + 1].t);
            if smooth {
                mults.push(2);
            } else {
                poles.push(p(i + 1));
                mults.push(3);
            }
        }
    }
    poles.push(p(n - 1));
    knots.push(nodes[n - 1].t);
    mults.push(4);
    if let Surface::Sphere(_) = untrimmed(surface) {
        let (_, _, v0, v1) = surface.bounds();
        for q in poles.iter_mut() {
            q.y = q.y.clamp(v0, v1);
        }
    }
    BSplineCurve::new(poles, knots, mults, 3)
}

fn untrimmed(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => untrimmed(t.basis()),
        s => s,
    }
}

#[cfg(test)]
mod tests;
