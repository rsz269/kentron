//! Sections of two surfaces traced by marching along them.
//!
//! Points where the surfaces meet are found from seeds and followed along the curve: each
//! step goes along the direction both surfaces share (the cross product of their normals)
//! and comes back onto both by Newton steps in both surfaces' parameters at once, the step's
//! length held; the step shortens where the curve turns. The seeds come from the base
//! surface's parameter range, split into cells while a cell may hold a point of the other
//! surface (the distance to it at the cell's corners and middle no larger than the cell),
//! at the changes of sign of that distance along the smallest cells' sides; and at the points
//! where edges of either face meet the other.
//!
//! A branch ends where it leaves the base's range or the other surface's bounds (at the
//! point where it crosses them), where it comes back to where it started (a closed curve),
//! or where the surfaces turn tangent (where the direction is lost). Both are followed a
//! little past their bounds, so that a section leaving a face crosses its boundary first.
//!
//! Where the surfaces meet nearly tangent along a branch (a fillet and the face it rolls
//! on), or a branch cannot be followed, the sections are not traced here: the caller traces
//! them on its grid.

use crate::geometry::ParametricSurface;
use crate::geometry::surface::Surface;
use crate::math::{Point, Vector};

use super::geom;

/// A range of a surface's parameters, with the periods of those that wrap round in it
/// (`u_period`, `v_period`) and the surface's own periods (`u_wrap`, `v_wrap`: a parameter
/// a whole number of them off is the same point).
#[derive(Clone, Copy, Debug)]
pub(super) struct Range {
    pub u: (f64, f64),
    pub v: (f64, f64),
    pub u_period: Option<f64>,
    pub v_period: Option<f64>,
    pub u_wrap: Option<f64>,
    pub v_wrap: Option<f64>,
}

impl Range {
    fn holds(&self, u: f64, v: f64) -> bool {
        let within = |x: f64, (a, b): (f64, f64), period: Option<f64>, wrap: Option<f64>| {
            if period.is_some() {
                return true;
            }
            let (lo, hi) = (a.min(b), a.max(b));
            // (Moved by whole periods of the surface to the middle of the range.)
            let x = wrap.map_or(x, |w| x + ((0.5 * (lo + hi) - x) / w).round() * w);
            x >= lo && x <= hi
        };
        within(u, self.u, self.u_period, self.u_wrap) && within(v, self.v, self.v_period, self.v_wrap)
    }
}

/// A traced section: the base's parameters along it and the other's, closed or not.
pub(super) struct Chain {
    pub uv: Vec<(f64, f64)>,
    pub st: Vec<(f64, f64)>,
    pub closed: bool,
}

/// The sine of the angle between the surfaces' normals below which they count as meeting
/// tangent.
const TANGENT: f64 = 0.05;

/// A point on both surfaces: its parameters on each, and where it is.
#[derive(Clone, Copy, Debug)]
struct Both {
    u: f64,
    v: f64,
    s: f64,
    t: f64,
    p: Point,
}

/// The chains of the base's parameters `(u, v)` along the sections of `base` (over `range`)
/// and `other` (within its bounds, or where it is unbounded near `other_range`, the other
/// face's), each closed or not; `None` where that is not done (the caller traces them
/// otherwise). The seeds are sought on an `nu` by `nv` lattice of cells, and at `crossings`
/// (points where edges of either face meet the other); `leaf` is the longest step (in
/// space); `slack` how far a crossing may lie off the surfaces.
pub(super) fn trace(base: &Surface, other: &Surface, range: Range, (nu, nv): (usize, usize), leaf: f64, crossings: &[Point], slack: f64, other_range: ((f64, f64), (f64, f64))) -> Option<Vec<Chain>> {
    if !(leaf > 0.0) {
        return None;
    }
    let (s0, s1, t0, t1) = other.bounds();
    let infinite = |x: f64| !x.is_finite() || x.abs() >= crate::tolerance::INFINITE;
    // The range followed on the other surface: a little past its bounds, as it is carried on
    // across its edges (so that a section leaving its face crosses the face's boundary
    // before it ends); where it is unbounded, the face's range, widened likewise.
    let past = |a: f64, b: f64| (a - 0.05 * (b - a), b + 0.05 * (b - a));
    let along = |lo: f64, hi: f64, (a, b): (f64, f64)| -> Option<(f64, f64)> {
        if !infinite(lo) && !infinite(hi) {
            return Some(past(lo, hi));
        }
        if infinite(a) || infinite(b) || !(b > a) {
            return None;
        }
        let (a, b) = (a - 0.1 * (b - a), b + 0.1 * (b - a));
        Some((a.max(lo), b.min(hi)))
    };
    let (Some(u), Some(v)) = (along(s0, s1, other_range.0), along(t0, t1, other_range.1)) else {
        return None;
    };
    let bounds = Range { u, v, u_period: other.u_period(), v_period: other.v_period(), u_wrap: other.u_period(), v_wrap: other.v_period() };
    let own = (s0, s1, t0, t1);
    let mut seeds = seeds(base, other, range, nu, nv);
    // (And the points where edges of either face meet the other, which the sections pass:
    // a section shorter than the lattice's cells, across a narrow face, is found there.)
    let lattice = seeds.len();
    for &p in crossings {
        let [u, v] = geom::uv_of(base, p);
        let [s, t] = geom::uv_of(other, p);
        if range.holds(u, v) {
            seeds.push(Both { u, v, s, t, p });
        }
    }
    // Followed a little past the base's bounds too, where the range reaches them (as far as
    // points on both surfaces may lie off them, by `slack`): a section there (across a
    // narrow face, at a corner of the base where edges of the other meet it only within
    // their tolerance) is followed past them rather than missed.
    let range = {
        let (bu0, bu1, bv0, bv1) = base.bounds();
        let (um, vm) = (0.5 * (range.u.0 + range.u.1), 0.5 * (range.v.0 + range.v.1));
        let widened = |(lo, hi): (f64, f64), (b0, b1): (f64, f64), periodic: bool, speed: &dyn Fn(f64) -> f64| -> (f64, f64) {
            if periodic {
                return (lo, hi);
            }
            let margin = |at: f64| (2.0 * slack / speed(at).max(1e-300)).min(0.05 * (hi - lo).abs());
            let at = |x: f64, b: f64| (x - b).abs() <= 1e-12 * (b1 - b0).abs().max(1.0);
            (if at(lo, b0) { lo - margin(lo) } else { lo }, if at(hi, b1) { hi + margin(hi) } else { hi })
        };
        let u = widened(range.u, (bu0, bu1), range.u_period.is_some(), &|u| base.d1(u, vm).1.length());
        let v = widened(range.v, (bv0, bv1), range.v_period.is_some(), &|v| base.d1(um, v).2.length());
        Range { u, v, ..range }
    };
    let tracer = Tracer { base, other, range, bounds, own, step: leaf, stopped: std::cell::Cell::new(false) };
    let mut chains: Vec<(Vec<Both>, bool)> = Vec::new();
    for (k, seed) in seeds.into_iter().enumerate() {
        // (Seeds on a branch traced already are passed over.)
        if chains.iter().any(|(points, closed)| near_polyline(points, *closed, seed.p, 0.75 * leaf)) {
            continue;
        }
        // (A crossing's point on both surfaces square to the curve from it, so as not to
        // slide along the curve away from it.)
        let start = if k < lattice {
            tracer.onto_both(seed)
        } else {
            tracer.direction(&seed).and_then(|dir| tracer.stepped(&seed, dir, 0.0, seed)).and_then(|b| tracer.held(b)).or_else(|| tracer.onto_both(seed))
        };
        let start = start.filter(|start| k < lattice || start.p.distance(seed.p) <= leaf);
        let Some(start) = start else {
            // (A crossing on both surfaces with no point of the sections near it in the
            // range: the section there runs past the base's bounds, which the point lies
            // beyond by its tolerance. Not followed here.)
            let on = |surface: &Surface| {
                let [x, y] = geom::uv_of(surface, seed.p);
                surface.point(x, y).distance(seed.p) <= slack
            };
            // (Only a crossing beyond the base's bounds: its nearest point is on them.)
            let beyond = || {
                let [x, y] = geom::uv_of(base, seed.p);
                let (u0, u1, v0, v1) = base.bounds();
                let at = |x: f64, a: f64, b: f64| (x - a).abs() <= 1e-9 * (b - a).abs().max(1.0) || (x - b).abs() <= 1e-9 * (b - a).abs().max(1.0);
                (base.u_period().is_none() && at(x, u0, u1)) || (base.v_period().is_none() && at(y, v0, v1))
            };
            // (Nor where the surfaces touch there, meeting tangent: no section is followed
            // from such a point, as along the edge a fillet rolls on a face by.)
            let touching = || tracer.sine(&seed) < TANGENT;
            if k >= lattice && on(base) && on(other) && beyond() && !touching() {
                return None;
            }
            continue;
        };
        if chains.iter().any(|(points, closed)| near_polyline(points, *closed, start.p, 0.75 * leaf)) {
            continue;
        }
        let (ahead, closed) = tracer.march(start, 1.0);
        let chain = if closed {
            let mut c = vec![start];
            c.extend(ahead);
            (c, true)
        } else {
            let (behind, _) = tracer.march(start, -1.0);
            let mut c: Vec<Both> = behind.into_iter().rev().collect();
            c.push(start);
            c.extend(ahead);
            (c, false)
        };
        // (Points bunched together, where the steps shortened as the branch ended, thinned:
        // curves fitted through them would overshoot. The ends kept.)
        let chain = (thinned(chain.0, 0.05 * leaf), chain.1);
        // (A branch stopped short: the section is not followed here.)
        if tracer.stopped.get() {
            return None;
        }
        // (Nor where the surfaces meet nearly tangent anywhere along it, as a fillet the face
        // it rolls on: where they only touch, or cross at a shallow angle, the sections the
        // grid finds are the ones the faces are split by.)
        if chain.0.iter().any(|b| tracer.sine(b) < TANGENT) {
            return None;
        }
        if chain.0.len() >= 2 {
            chains.push(chain);
        }
    }
    Some(
        chains
            .into_iter()
            .map(|(points, closed)| Chain {
                uv: unwrapped(points.iter().map(|b| (b.u, b.v)), range),
                st: unwrapped(points.iter().map(|b| (b.s, b.t)), tracer.bounds),
                closed,
            })
            .collect(),
    )
}

/// `points` without those nearer than `spacing` to the one kept before them (the last kept
/// in place of the one before it, where that is too near).
fn thinned(points: Vec<Both>, spacing: f64) -> Vec<Both> {
    let n = points.len();
    let mut out: Vec<Both> = Vec::with_capacity(n);
    for (i, b) in points.into_iter().enumerate() {
        match out.last() {
            Some(last) if last.p.distance(b.p) < spacing => {
                // (The last kept, in place of the one before it unless that is the first:
                // a branch shorter than the spacing keeps both its ends.)
                if i + 1 == n && out.len() > 1 {
                    *out.last_mut().expect("a point") = b;
                } else if i + 1 == n && last.p.distance(b.p) > 0.0 {
                    out.push(b);
                }
            }
            _ => out.push(b),
        }
    }
    out
}

/// Parameters along a chain, unwrapped across the periods of `range` (the first brought
/// into it).
fn unwrapped(points: impl Iterator<Item = (f64, f64)>, range: Range) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    for (mut u, mut v) in points {
        match out.last() {
            Some(&(pu, pv)) => {
                if let Some(period) = range.u_wrap {
                    u += ((pu - u) / period).round() * period;
                }
                if let Some(period) = range.v_wrap {
                    v += ((pv - v) / period).round() * period;
                }
            }
            None => {
                let middle = |x: f64, (a, b): (f64, f64), wrap: Option<f64>| wrap.map_or(x, |w| x + ((0.5 * (a + b) - x) / w).round() * w);
                (u, v) = (middle(u, range.u, range.u_wrap), middle(v, range.v, range.v_wrap));
            }
        }
        out.push((u, v));
    }
    out
}

/// Whether `p` lies within `reach` of the polyline through `points` (closed or not).
fn near_polyline(points: &[Both], closed: bool, p: Point, reach: f64) -> bool {
    if points.len() == 1 {
        return points[0].p.distance(p) <= reach;
    }
    let n = points.len();
    let segments = if closed { n } else { n - 1 };
    (0..segments).any(|i| {
        let (a, b) = (points[i].p, points[(i + 1) % n].p);
        let d = b - a;
        let len2 = d.dot(d);
        let s = if len2 > 0.0 { ((p - a).dot(d) / len2).clamp(0.0, 1.0) } else { 0.0 };
        (a + d * s).distance(p) <= reach
    })
}

/// Seeds of the sections: points of the base (with a guess on the other surface) where the
/// distance to the other changes sign along the sides of the cells of an `nu` by `nv`
/// lattice over the range, sought only in blocks of cells that may hold a point of the
/// other (each node's distance found once).
fn seeds(base: &Surface, other: &Surface, range: Range, nu: usize, nv: usize) -> Vec<Both> {
    let (u0, u1, v0, v1) = (range.u.0, range.u.1, range.v.0, range.v.1);
    let at = |i: usize, j: usize| (u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
    // (Each distance with the foot found for it, which starts the next search.)
    let last = std::cell::Cell::new(None);
    let mut nodes: crate::hashing::FastMap<(usize, usize), (f64, Point)> = Default::default();
    let mut node = |i: usize, j: usize| -> (f64, Point) {
        *nodes.entry((i, j)).or_insert_with(|| {
            let (u, v) = at(i, j);
            let p = base.point(u, v);
            let (d, foot) = geom::signed_distance_from(other, p, last.get());
            last.set(foot);
            (d, p)
        })
    };
    let mut out = Vec::new();
    // Blocks of cells (by node indices), halved while they may hold a point of the other.
    let mut stack: Vec<(usize, usize, usize, usize)> = vec![(0, nu, 0, nv)];
    while let Some((ia, ib, ja, jb)) = stack.pop() {
        let corners = [(ia, ja), (ib, ja), (ib, jb), (ia, jb)];
        let values: Vec<(f64, Point)> = corners.iter().map(|&(i, j)| node(i, j)).collect();
        let (ic, jc) = ((ia + ib) / 2, (ja + jb) / 2);
        let (dc, pc) = node(ic, jc);
        let radius = values.iter().map(|(_, p)| p.distance(pc)).fold(0.0, f64::max);
        let near = values.iter().map(|(d, _)| d.abs()).fold(dc.abs(), f64::min);
        let sign = values[0].0 >= 0.0;
        let changes = values.iter().any(|(d, _)| (*d >= 0.0) != sign) || (dc >= 0.0) != sign;
        // (No point of the other within the block: farther from it than the block is wide.)
        if !changes && near > 1.5 * radius + geom::on_tol() {
            continue;
        }
        if ib - ia > 1 || jb - ja > 1 {
            // Halved along the longer side (in cells).
            if ib - ia >= jb - ja {
                stack.push((ia, ic, ja, jb));
                stack.push((ic, ib, ja, jb));
            } else {
                stack.push((ia, ib, ja, jc));
                stack.push((ia, ib, jc, jb));
            }
            continue;
        }
        // A single cell: along each side changing sign, where the distance is 0 (linearly).
        for k in 0..4 {
            let (a, b) = (k, (k + 1) % 4);
            let ((da, pa), (db, pb)) = (values[a], values[b]);
            if (da >= 0.0) == (db >= 0.0) || da == db {
                continue;
            }
            let w = da / (da - db);
            let ((ua, va), (ub, vb)) = (at(corners[a].0, corners[a].1), at(corners[b].0, corners[b].1));
            let (u, v) = (ua + (ub - ua) * w, va + (vb - va) * w);
            let p = pa + (pb - pa) * w;
            let [s, t] = geom::uv_of(other, base.point(u, v));
            out.push(Both { u, v, s, t, p });
        }
    }
    out
}

/// What marching needs: both surfaces, the ranges it stays in, the longest step.
struct Tracer<'a> {
    base: &'a Surface,
    other: &'a Surface,
    range: Range,
    bounds: Range,
    /// The other surface's own bounds (`bounds` reach past them).
    own: (f64, f64, f64, f64),
    step: f64,
    /// A branch stopped short: where the surfaces turned tangent, or its steps failed.
    stopped: std::cell::Cell<bool>,
}

impl Tracer<'_> {
    /// The other surface's point and derivatives at `(s, t)`; past its bounds, carried on
    /// across its edge in its tangent plane, square to the edge.
    fn other_d1(&self, s: f64, t: f64) -> (Point, Vector, Vector) {
        let (s0, s1, t0, t1) = self.own;
        let clamp = |x: f64, lo: f64, hi: f64, periodic: bool| if periodic { x } else { x.clamp(lo, hi) };
        let (sc, tc) = (clamp(s, s0, s1, self.bounds.u_period.is_some()), clamp(t, t0, t1, self.bounds.v_period.is_some()));
        let (q, qs, qt) = self.other.d1(sc, tc);
        let (ds, dt) = (s - sc, t - tc);
        if ds == 0.0 && dt == 0.0 {
            return (q, qs, qt);
        }
        let across = |a: Vector, b: Vector| {
            let bb = b.dot(b);
            if bb > 0.0 { a - b * (a.dot(b) / bb) } else { a }
        };
        let (es, et) = (if ds != 0.0 { across(qs, qt) } else { qs }, if dt != 0.0 { across(qt, qs) } else { qt });
        (q + es * ds + et * dt, es, et)
    }

    /// The point both surfaces share near `guess` (Newton steps moving all four parameters
    /// the least); `None` if they do not converge within the ranges.
    fn onto_both(&self, guess: Both) -> Option<Both> {
        let (mut u, mut v, mut s, mut t) = (guess.u, guess.v, guess.s, guess.t);
        for _ in 0..30 {
            let (p, pu, pv) = self.base.d1(u, v);
            let (q, qs, qt) = self.other_d1(s, t);
            let f = p - q;
            if f.length() <= 1e-12 * (1.0 + (p - Point::ORIGIN).length()) {
                return (self.range.holds(u, v) && self.bounds.holds(s, t)).then_some(Both { u, v, s, t, p });
            }
            let step = min_norm([pu, pv, -qs, -qt], f)?;
            (u, v, s, t) = (u + step[0], v + step[1], s + step[2], t + step[3]);
            if ![u, v, s, t].iter().all(|x| x.is_finite()) {
                return None;
            }
        }
        None
    }

    /// `b` if within the ranges, or a hair outside the base's (a crossing on its bounds),
    /// moved onto them.
    fn held(&self, b: Both) -> Option<Both> {
        if !self.bounds.holds(b.s, b.t) {
            return None;
        }
        if self.range.holds(b.u, b.v) {
            return Some(b);
        }
        let snap = |x: f64, (lo, hi): (f64, f64)| {
            let hair = 1e-8 * (hi - lo).abs();
            if x < lo && x >= lo - hair { lo } else if x > hi && x <= hi + hair { hi } else { x }
        };
        let (u, v) = (snap(b.u, self.range.u), snap(b.v, self.range.v));
        if !self.range.holds(u, v) {
            return None;
        }
        self.on_bound(&Both { u, v, ..b }, if u != b.u { 0 } else { 1 }, if u != b.u { u } else { v })
    }

    /// The sine of the angle between the surfaces' normals at a point.
    fn sine(&self, b: &Both) -> f64 {
        let (_, pu, pv) = self.base.d1(b.u, b.v);
        let (_, qs, qt) = self.other_d1(b.s, b.t);
        let (n1, n2) = (pu.cross(pv), qs.cross(qt));
        let scale = n1.length() * n2.length();
        if scale > 0.0 { n1.cross(n2).length() / scale } else { 0.0 }
    }

    /// The direction both surfaces share at a point, unit; `None` where they are tangent.
    fn direction(&self, b: &Both) -> Option<Vector> {
        let (_, pu, pv) = self.base.d1(b.u, b.v);
        let (_, qs, qt) = self.other_d1(b.s, b.t);
        let (n1, n2) = (pu.cross(pv), qs.cross(qt));
        let d = n1.cross(n2);
        let len = d.length();
        (len > 1e-6 * n1.length() * n2.length() && len > 0.0).then(|| d * (1.0 / len))
    }

    /// The point on both surfaces `h` along `dir` from `from` (on the plane square to `dir`
    /// there), from a guess.
    fn stepped(&self, from: &Both, dir: Vector, h: f64, guess: Both) -> Option<Both> {
        let (mut u, mut v, mut s, mut t) = (guess.u, guess.v, guess.s, guess.t);
        for _ in 0..20 {
            let (p, pu, pv) = self.base.d1(u, v);
            let (q, qs, qt) = self.other_d1(s, t);
            let f = p - q;
            let g = (p - from.p).dot(dir) - h;
            let scale = 1e-12 * (1.0 + (p - Point::ORIGIN).length());
            if f.length() <= scale && g.abs() <= scale {
                return Some(Both { u, v, s, t, p });
            }
            // Rows: the three of base - other, and the plane's; columns u, v, s, t.
            let m = [
                [pu.x, pv.x, -qs.x, -qt.x],
                [pu.y, pv.y, -qs.y, -qt.y],
                [pu.z, pv.z, -qs.z, -qt.z],
                [pu.dot(dir), pv.dot(dir), 0.0, 0.0],
            ];
            let step = solve4(m, [-f.x, -f.y, -f.z, -g])?;
            (u, v, s, t) = (u + step[0], v + step[1], s + step[2], t + step[3]);
            if ![u, v, s, t].iter().all(|x| x.is_finite()) {
                return None;
            }
        }
        None
    }

    /// The branch from `start` one way (`sense` +1 or -1 along the shared direction): its
    /// points after `start`, and whether it came back to `start` (closed).
    fn march(&self, start: Both, sense: f64) -> (Vec<Both>, bool) {
        let mut out: Vec<Both> = Vec::new();
        // (From a point where the surfaces touch, no branch: they meet there only.)
        let Some(mut dir) = self.direction(&start).map(|d| d * sense) else {
            return (out, false);
        };
        // Where the steps fail: the branch ends there where the surfaces turn tangent
        // (meeting there, or touching along a curve, as a fillet the face it rolls on), else
        // it is not followed.
        let failed = |at: &Both| {
            if self.sine(at) > TANGENT {
                self.stopped.set(true);
            }
        };
        let mut at = start;
        let mut h = self.step;
        let smallest = 1e-6 * self.step;
        let mut travelled = 0.0;
        // (A bound on the steps: a curve winding without end, as along a tangency.)
        for _ in 0..20_000 {
            let guess_point = at.p + dir * h;
            let guess = self.guess_near(&at, guess_point);
            // (Not where the steps led to another crossing of the plane, far along the
            // curve or on another branch: a step reaches about as far as it is long.)
            let next = self.stepped(&at, dir, h, guess).filter(|n| n.p.distance(at.p) <= 1.5 * h);
            let Some(next) = next else {
                if h <= smallest {
                    failed(&at);
                    break;
                }
                h *= 0.5;
                continue;
            };
            let Some(next_dir) = self.direction(&next).map(|d| if d.dot(dir) < 0.0 { d * -1.0 } else { d }) else {
                // Tangent there: the branch stops (at the last point found).
                if h <= smallest {
                    failed(&at);
                    break;
                }
                h *= 0.5;
                continue;
            };
            // Shorter where the curve turns faster.
            let turn = dir.dot(next_dir).clamp(-1.0, 1.0).acos();
            if turn > 0.15 && h > smallest {
                h *= 0.5;
                continue;
            }
            // And where its middle is not found near the chord's (as the curve is evaluated
            // between its samples, square to the chord): where the surfaces turn tangent, as
            // the point found in a plane slides along them.
            let chord = next.p - at.p;
            let middle = at.p + chord * 0.5;
            let mid_guess = Both { u: 0.5 * (at.u + next.u), v: 0.5 * (at.v + next.v), s: 0.5 * (at.s + next.s), t: 0.5 * (at.t + next.t), p: middle };
            let across = chord.normalized();
            let true_middle = self.stepped(&Both { p: middle, ..at }, across, 0.0, mid_guess);
            if h > smallest && !true_middle.is_some_and(|m| m.p.distance(middle) <= 0.1 * chord.length()) {
                h *= 0.5;
                continue;
            }
            // Out of the ranges: the point where it leaves them ends the branch.
            if !self.range.holds(next.u, next.v) || !self.bounds.holds(next.s, next.t) {
                if let Some(end) = self.exit(&at, &next, dir, h) {
                    out.push(end);
                }
                return (out, false);
            }
            // Back at the start: a closed curve.
            travelled += h;
            // (In the parameters too, up to whole periods: a curve reaching across a closed
            // but not periodic surface's ends comes back to where it started only in space.)
            let back = |x: f64, x0: f64, step: f64, wrap: Option<f64>| {
                let d = wrap.map_or(x0 - x, |w| (x0 - x) - ((x0 - x) / w).round() * w);
                d.abs() <= 3.0 * step + 1e-12
            };
            let returned = back(next.u, start.u, (next.u - at.u).abs(), self.range.u_wrap) && back(next.v, start.v, (next.v - at.v).abs(), self.range.v_wrap);
            if returned && travelled > 2.5 * self.step && next.p.distance(start.p) <= 1.01 * h && (start.p - at.p).dot(dir) > 0.0 {
                // (The step back kept, unless at the start: the closing chord no longer than
                // a step.)
                if next.p.distance(start.p) > 0.05 * h {
                    out.push(next);
                }
                return (out, true);
            }
            out.push(next);
            at = next;
            dir = next_dir;
            if turn < 0.05 {
                h = (h * 1.5).min(self.step);
            }
        }
        (out, false)
    }

    /// A guess at the parameters of the point on both surfaces near `p`, from those at
    /// `at` moved by the derivatives.
    fn guess_near(&self, at: &Both, p: Point) -> Both {
        let d = p - at.p;
        let (_, pu, pv) = self.base.d1(at.u, at.v);
        let (_, qs, qt) = self.other_d1(at.s, at.t);
        let solve2 = |a: Vector, b: Vector| {
            let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
            let det = aa * bb - ab * ab;
            if det.abs() <= 1e-300 {
                return (0.0, 0.0);
            }
            let (ra, rb) = (d.dot(a), d.dot(b));
            ((ra * bb - rb * ab) / det, (rb * aa - ra * ab) / det)
        };
        let (du, dv) = solve2(pu, pv);
        let (ds, dt) = solve2(qs, qt);
        Both { u: at.u + du, v: at.v + dv, s: at.s + ds, t: at.t + dt, p }
    }

    /// Where the branch leaves the ranges between `at` (within them) and `next` (out of
    /// them): on the bound it crosses, the parameter held there and the others found by
    /// Newton steps; else by halving the step.
    fn exit(&self, at: &Both, next: &Both, dir: Vector, h: f64) -> Option<Both> {
        // The bound crossed: (which parameter, its value there).
        let crossed = |x: f64, (a, b): (f64, f64), period: Option<f64>, wrap: Option<f64>| -> Option<f64> {
            if period.is_some() {
                return None;
            }
            let (lo, hi) = (a.min(b), a.max(b));
            let x = wrap.map_or(x, |w| x + ((0.5 * (lo + hi) - x) / w).round() * w);
            if x < lo { Some(lo) } else if x > hi { Some(hi) } else { None }
        };
        let r = &self.range;
        let b = &self.bounds;
        let candidates = [
            (0, crossed(next.u, r.u, r.u_period, r.u_wrap)),
            (1, crossed(next.v, r.v, r.v_period, r.v_wrap)),
            (2, crossed(next.s, b.u, b.u_period, b.u_wrap)),
            (3, crossed(next.t, b.v, b.v_period, b.v_wrap)),
        ];
        // (A point a hair from a corner of the base's range: at the corner.)
        let cornered = |end: Both| -> Both {
            let near = |x: f64, (a, b): (f64, f64)| {
                let hair = 1e-4 * (b - a).abs();
                if (x - a).abs() <= hair { Some(a) } else if (x - b).abs() <= hair { Some(b) } else { None }
            };
            let (u, v) = (near(end.u, r.u).filter(|_| r.u_period.is_none()), near(end.v, r.v).filter(|_| r.v_period.is_none()));
            let (Some(u), Some(v)) = (u, v) else { return end };
            let p = self.base.point(u, v);
            // (The other's parameters there by Newton steps on the distance.)
            let (mut s, mut t) = (end.s, end.t);
            for _ in 0..20 {
                let (q, qs, qt) = self.other_d1(s, t);
                let f = q - p;
                let (a, b, c) = (qs.dot(qs), qs.dot(qt), qt.dot(qt));
                let det = a * c - b * b;
                if det.abs() <= 1e-300 {
                    break;
                }
                let (rs, rt) = (f.dot(qs), f.dot(qt));
                let (ds, dt) = ((rs * c - rt * b) / det, (rt * a - rs * b) / det);
                (s, t) = (s - ds, t - dt);
                if ds.abs() + dt.abs() <= 1e-15 {
                    break;
                }
            }
            if self.other_d1(s, t).0.distance(p) <= geom::on_tol() { Both { u, v, s, t, p } } else { end }
        };
        for (which, bound) in candidates {
            let Some(bound) = bound else { continue };
            if let Some(end) = self.on_bound(at, which, bound).filter(|e| e.p.distance(at.p) <= 2.0 * h + 1e-9) {
                return Some(cornered(end));
            }
        }
        // (Else by halving the step.)
        let (mut lo, mut hi) = (0.0, h);
        let mut best: Option<Both> = None;
        let mut from = *at;
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            let guess = self.guess_near(&from, at.p + dir * mid);
            let Some(p) = self.stepped(at, dir, mid, guess) else {
                hi = mid;
                continue;
            };
            if self.range.holds(p.u, p.v) && self.bounds.holds(p.s, p.t) {
                lo = mid;
                best = Some(p);
                from = p;
            } else {
                hi = mid;
            }
            if hi - lo <= 1e-12 * self.step.max(1.0) {
                break;
            }
        }
        best
    }

    /// The point on both surfaces with parameter `which` (0 u, 1 v, 2 s, 3 t) at `bound`,
    /// near `at`.
    fn on_bound(&self, at: &Both, which: usize, bound: f64) -> Option<Both> {
        let mut x = [at.u, at.v, at.s, at.t];
        // (Moved by whole periods next to the bound first.)
        let wrap = [self.range.u_wrap, self.range.v_wrap, self.bounds.u_wrap, self.bounds.v_wrap][which];
        if let Some(w) = wrap {
            x[which] += ((bound - x[which]) / w).round() * w;
        }
        x[which] = bound;
        for _ in 0..30 {
            let (p, pu, pv) = self.base.d1(x[0], x[1]);
            let (q, qs, qt) = self.other_d1(x[2], x[3]);
            let f = p - q;
            if f.length() <= 1e-12 * (1.0 + (p - Point::ORIGIN).length()) {
                let b = Both { u: x[0], v: x[1], s: x[2], t: x[3], p };
                return (self.range.holds(b.u, b.v) && self.bounds.holds(b.s, b.t)).then_some(b);
            }
            let columns = [pu, pv, -qs, -qt];
            let free: Vec<usize> = (0..4).filter(|&k| k != which).collect();
            let m = [
                [columns[free[0]].x, columns[free[1]].x, columns[free[2]].x],
                [columns[free[0]].y, columns[free[1]].y, columns[free[2]].y],
                [columns[free[0]].z, columns[free[1]].z, columns[free[2]].z],
            ];
            let step = solve3(m, [-f.x, -f.y, -f.z])?;
            for (k, &i) in free.iter().enumerate() {
                x[i] += step[k];
            }
            if !x.iter().all(|v| v.is_finite()) {
                return None;
            }
        }
        None
    }
}

/// The least-norm `x` with `sum columns[k] x[k] = -f` (four unknowns, three equations).
fn min_norm(columns: [Vector; 4], f: Vector) -> Option<[f64; 4]> {
    // x = -J^T (J J^T)^-1 f.
    let mut m = [[0.0; 3]; 3];
    for c in &columns {
        let a = [c.x, c.y, c.z];
        for r in 0..3 {
            for k in 0..3 {
                m[r][k] += a[r] * a[k];
            }
        }
    }
    let y = solve3(m, [f.x, f.y, f.z])?;
    let y = Vector::new(y[0], y[1], y[2]);
    Some(columns.map(|c| -c.dot(y)))
}

/// `m x = b` for 3x3 `m` (Cramer's rule); `None` if `m` is (nearly) singular.
fn solve3(m: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det3 = |a: [[f64; 3]; 3]| a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let det = det3(m);
    let scale = m.iter().flatten().fold(0.0f64, |s, x| s.max(x.abs()));
    if !(det.abs() > 1e-24 * scale * scale * scale) {
        return None;
    }
    let mut out = [0.0; 3];
    for (k, o) in out.iter_mut().enumerate() {
        let mut a = m;
        for r in 0..3 {
            a[r][k] = b[r];
        }
        *o = det3(a) / det;
    }
    Some(out)
}

/// `m x = b` for 4x4 `m` (Gaussian elimination with partial pivoting); `None` if singular.
pub(super) fn solve4(mut m: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    let scale = m.iter().flatten().fold(0.0f64, |s, x| s.max(x.abs()));
    if !(scale > 0.0) {
        return None;
    }
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))?;
        if m[pivot][col].abs() <= 1e-14 * scale {
            return None;
        }
        m.swap(col, pivot);
        b.swap(col, pivot);
        for r in col + 1..4 {
            let factor = m[r][col] / m[col][col];
            for k in col..4 {
                m[r][k] -= factor * m[col][k];
            }
            b[r] -= factor * b[col];
        }
    }
    let mut x = [0.0; 4];
    for r in (0..4).rev() {
        let sum: f64 = (r + 1..4).map(|k| m[r][k] * x[k]).sum();
        x[r] = (b[r] - sum) / m[r][r];
    }
    Some(x)
}
