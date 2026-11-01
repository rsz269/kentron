//! Drafted extrusions of straight-edged outlines whose edges the draft closes: the outline's
//! offsets (each edge moving in at unit speed, as the sides lean in) followed through the
//! events where an edge shrinks to nothing (its neighbours then meet) and where a reflex
//! corner runs into an edge across (the outline splits), up to the top's offset: the
//! straight skeleton's roof, cut off at the top. Each side is the part of its edge's
//! drafted plane between the edge, the skeleton's arcs from its ends and the edge's pieces
//! at the top (or the ridges where the outline closes below it).

use crate::error::{Error, Result};
use crate::make;
use crate::math::{Point, Vector};
use crate::properties;
use crate::sew::{self, SewOptions};
use crate::topology::{Model, NodeKind, Shape};

type P2 = [f64; 2];

fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn scale(a: P2, k: f64) -> P2 {
    [a[0] * k, a[1] * k]
}
fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn length(a: P2) -> f64 {
    dot(a, a).sqrt()
}

/// An edge of the outline: where it starts, its direction, and the normal it moves along
/// (into the outline for a draft leaning in, out of it for one leaning out).
struct Line {
    start: P2,
    direction: P2,
    normal: P2,
}

/// A corner of the moving outline: where it was at time `t`, its velocity since, and the
/// edges before and after it.
#[derive(Clone, Copy)]
struct Corner {
    at: P2,
    t: f64,
    velocity: P2,
    before: usize,
    after: usize,
}

impl Corner {
    fn position(&self, t: f64) -> P2 {
        add(self.at, scale(self.velocity, t - self.t))
    }
}

/// Where a side face's boundary runs: from one point at a time (offset) to another.
type Segment = ((P2, f64), (P2, f64));

struct Wavefront {
    lines: Vec<Line>,
    corners: Vec<Corner>,
    loops: Vec<Vec<usize>>,
    /// Per edge, the pieces of its side face's boundary.
    sides: Vec<Vec<Segment>>,
    /// The outlines at the top.
    tops: Vec<Vec<P2>>,
    size: f64,
}

impl Wavefront {
    /// The corner between edges `before` and `after` at `at`, time `t`: moving so as to stay on
    /// both edges as they move.
    fn corner(&self, at: P2, t: f64, before: usize, after: usize) -> Corner {
        let (a, b) = (self.lines[before].normal, self.lines[after].normal);
        let k = 1.0 + dot(a, b);
        // (Edges turned back on one another: the corner stays where they meet.)
        let velocity = if k > 1e-12 { scale(add(a, b), 1.0 / k) } else { [0.0, 0.0] };
        Corner { at, t, velocity, before, after }
    }

    fn reflex(&self, c: &Corner) -> bool {
        // (Turning against the way the edges move: the corner runs into the outline.)
        let (a, b) = (&self.lines[c.before], &self.lines[c.after]);
        cross(a.direction, b.direction) * dot(a.normal, [-a.direction[1], a.direction[0]]) < -1e-12
    }

    /// The pieces of the outline's edges between the loops' corners at time `t`, each added
    /// to its edge's side.
    fn close_loop(&mut self, lp: &[usize], t: f64) {
        for (i, &c) in lp.iter().enumerate() {
            let d = lp[(i + 1) % lp.len()];
            let (a, b) = (self.corners[c], self.corners[d]);
            self.sides[a.after].push(((a.position(t), t), (b.position(t), t)));
        }
    }

    /// Loop `lp` at time `t` with its spikes taken out: a corner between edges turned back on
    /// one another (lying along one line, as where a split meets an edge along which the
    /// corner's own edge runs) is the tip of a spike of no width. Along the stretch both
    /// edges share, each edge's side ends in a ridge; the edge reaching less far from the
    /// tip ends there, the other running on to the next corner, made anew between it and
    /// that corner's other edge. Returns the loop left (fewer than 3 corners: closed).
    fn without_spikes(&mut self, mut lp: Vec<usize>, t: f64) -> Vec<usize> {
        'again: loop {
            let m = lp.len();
            if m < 3 {
                return lp;
            }
            for k in 0..m {
                let r = self.corners[lp[k]];
                let (a, b) = (&self.lines[r.before], &self.lines[r.after]);
                if 1.0 + dot(a.normal, b.normal) > 1e-9 || dot(a.direction, b.direction) > -1.0 + 1e-9 {
                    continue;
                }
                let (pk, qk) = ((k + m - 1) % m, (k + 1) % m);
                let (p, q) = (self.corners[lp[pk]], self.corners[lp[qk]]);
                let tip = r.position(t);
                let (pp, qp) = (p.position(t), q.position(t));
                // (How far each neighbour lies back from the tip along the first edge.)
                let back = |x: P2| dot(sub(tip, x), a.direction);
                let (bp, bq) = (back(pp), back(qp));
                self.end(lp[k], t);
                if bq < bp - 1e-12 * self.size {
                    // The first edge runs on from p to q; the second ends at q.
                    self.sides[r.before].push(((qp, t), (tip, t)));
                    self.sides[r.after].push(((tip, t), (qp, t)));
                    self.end(lp[qk], t);
                    let c = self.corner(qp, t, r.before, q.after);
                    self.corners.push(c);
                    lp[qk] = self.corners.len() - 1;
                } else if bp < bq - 1e-12 * self.size {
                    // The second edge runs on from p to q; the first ends at p.
                    self.sides[r.before].push(((pp, t), (tip, t)));
                    self.sides[r.after].push(((tip, t), (pp, t)));
                    self.end(lp[pk], t);
                    let c = self.corner(pp, t, p.before, r.after);
                    self.corners.push(c);
                    lp[pk] = self.corners.len() - 1;
                } else {
                    // Both end at the one point: p and q one corner.
                    self.sides[r.before].push(((pp, t), (tip, t)));
                    self.sides[r.after].push(((tip, t), (qp, t)));
                    self.end(lp[pk], t);
                    self.end(lp[qk], t);
                    let c = self.corner(pp, t, p.before, q.after);
                    self.corners.push(c);
                    let (tip_id, q_id) = (lp[k], lp[qk]);
                    lp[pk] = self.corners.len() - 1;
                    lp.retain(|&v| v != tip_id && v != q_id);
                    continue 'again;
                }
                lp.remove(k);
                continue 'again;
            }
            return lp;
        }
    }

    /// The corner `c` ending at time `t`: the arc it ran along added to both its edges' sides.
    fn end(&mut self, c: usize, t: f64) {
        let k = self.corners[c];
        let arc = ((k.at, k.t), (k.position(t), t));
        if length(sub(arc.0 .0, arc.1 .0)) > 1e-12 * self.size || (arc.0 .1 - arc.1 .1).abs() > 1e-12 * self.size {
            self.sides[k.before].push(arc);
            self.sides[k.after].push(arc);
        }
    }
}

/// The next event after `now`: (time, loop, kind) with kind 0 an edge shrinking to nothing
/// (between positions `i` and `i + 1` of the loop) and 1 a reflex corner (position `i`) running
/// into an edge (between positions `j` and `j + 1` of loop `other`).
#[derive(Clone, Copy, Debug)]
enum Event {
    Edge { t: f64, lp: usize, i: usize },
    Split { t: f64, lp: usize, i: usize, other: usize, j: usize },
}

impl Event {
    fn time(&self) -> f64 {
        match *self {
            Event::Edge { t, .. } | Event::Split { t, .. } => t,
        }
    }
}

fn next_event(w: &Wavefront, now: f64, until: f64) -> Option<Event> {
    let eps = 1e-12 * w.size;
    let mut best: Option<Event> = None;
    let better = |e: &Event, best: &Option<Event>| best.as_ref().is_none_or(|b| e.time() < b.time() - eps || (e.time() <= b.time() + eps && matches!(e, Event::Edge { .. }) && matches!(b, Event::Split { .. })));
    for (l, lp) in w.loops.iter().enumerate() {
        let n = lp.len();
        for i in 0..n {
            let (a, b) = (w.corners[lp[i]], w.corners[lp[(i + 1) % n]]);
            let line = &w.lines[a.after];
            // (When the two corners' positions along the edge meet.)
            let along = |c: &Corner| (dot(sub(c.at, line.start), line.direction), dot(c.velocity, line.direction), c.t);
            let ((pa, va, ta), (pb, vb, tb)) = (along(&a), along(&b));
            let dv = vb - va;
            if dv < -1e-14 {
                let t = (pa - va * ta - pb + vb * tb) / dv;
                if t >= now - eps && t <= until + eps {
                    let e = Event::Edge { t: t.max(now), lp: l, i };
                    if better(&e, &best) {
                        best = Some(e);
                    }
                }
            }
        }
        for i in 0..n {
            let r = w.corners[lp[i]];
            if !w.reflex(&r) {
                continue;
            }
            for (o, other) in w.loops.iter().enumerate() {
                let m = other.len();
                for j in 0..m {
                    let (ca, cb) = (other[j], other[(j + 1) % m]);
                    if ca == lp[i] || cb == lp[i] {
                        continue;
                    }
                    let (a, b) = (w.corners[ca], w.corners[cb]);
                    let line = &w.lines[a.after];
                    // (When the corner reaches the edge's line, as both move.)
                    let speed = 1.0 - dot(r.velocity, line.normal);
                    if speed <= 1e-12 {
                        continue;
                    }
                    let t = dot(sub(sub(r.at, scale(r.velocity, r.t)), line.start), line.normal) / speed;
                    if !(t > now + eps && t <= until + eps) {
                        continue;
                    }
                    let x = r.position(t);
                    let (pa, pb) = (a.position(t), b.position(t));
                    let span = dot(sub(pb, pa), line.direction);
                    let s = dot(sub(x, pa), line.direction);
                    if span > 0.0 && s > eps && s < span - eps {
                        let e = Event::Split { t, lp: l, i, other: o, j };
                        if better(&e, &best) {
                            best = Some(e);
                        }
                    }
                }
            }
        }
    }
    best
}

/// The outline `points` (2D, in order, closed, counter-clockwise) offset by its straight
/// skeleton up to `depth` (in for a positive one, out for a negative): per edge, its side's
/// boundary (pieces at times 0 to |depth|), and the outlines at the top.
fn skeleton(points: &[P2], depth: f64) -> Result<(Vec<Vec<Segment>>, Vec<Vec<P2>>)> {
    let n = points.len();
    let inward = depth > 0.0;
    let until = depth.abs();
    let size = points.iter().map(|p| length(sub(*p, points[0]))).fold(0.0, f64::max).max(1e-9);
    let lines: Vec<Line> = (0..n)
        .map(|i| {
            let (a, b) = (points[i], points[(i + 1) % n]);
            let d = sub(b, a);
            let l = length(d);
            let direction = scale(d, 1.0 / l);
            let left = [-direction[1], direction[0]];
            Line { start: a, direction, normal: if inward { left } else { scale(left, -1.0) } }
        })
        .collect();
    let mut w = Wavefront { lines, corners: Vec::new(), loops: Vec::new(), sides: vec![Vec::new(); n], tops: Vec::new(), size };
    for i in 0..n {
        let c = w.corner(points[i], 0.0, (i + n - 1) % n, i);
        w.corners.push(c);
        w.sides[i].push(((points[i], 0.0), (points[(i + 1) % n], 0.0)));
    }
    w.loops.push((0..n).collect());
    let mut now = 0.0;
    let mut steps = 0;
    while let Some(event) = next_event(&w, now, until) {
        steps += 1;
        if steps > 10 * n + 100 {
            return Err(Error::Failed("a drafted outline's offsets that do not settle".into()));
        }
        now = event.time();
        match event {
            Event::Edge { t, lp, i } => {
                let lp_corners = w.loops[lp].clone();
                let m = lp_corners.len();
                let (a, b) = (lp_corners[i], lp_corners[(i + 1) % m]);
                let x = scale(add(w.corners[a].position(t), w.corners[b].position(t)), 0.5);
                w.end(a, t);
                w.end(b, t);
                // (The vanished edge's side ends at the point.)
                let (before, after) = (w.corners[a].before, w.corners[b].after);
                let c = w.corner(x, t, before, after);
                w.corners.push(c);
                let id = w.corners.len() - 1;
                let mut rest: Vec<usize> = Vec::with_capacity(m - 1);
                for k in 0..m {
                    let v = lp_corners[(i + 2 + k) % m];
                    if v != a && v != b {
                        rest.push(v);
                    }
                }
                rest.push(id);
                let rest = w.without_spikes(rest, t);
                if rest.len() <= 2 {
                    // (Closed: the rest of the outline is a ridge at this time.)
                    w.close_loop(&rest, t);
                    for &v in &rest {
                        w.end(v, t);
                    }
                    w.loops.remove(lp);
                } else {
                    w.loops[lp] = rest;
                }
            }
            Event::Split { t, lp, i, other, j } => {
                let r = w.loops[lp][i];
                let (ca, cb) = (w.loops[other][j], w.loops[other][(j + 1) % w.loops[other].len()]);
                let edge = w.corners[ca].after;
                let x = w.corners[r].position(t);
                w.end(r, t);
                let (before, after) = (w.corners[r].before, w.corners[r].after);
                let r1 = w.corner(x, t, before, edge);
                let r2 = w.corner(x, t, edge, after);
                w.corners.push(r1);
                w.corners.push(r2);
                let (id1, id2) = (w.corners.len() - 2, w.corners.len() - 1);
                if lp == other {
                    // One loop in two: r1 then on from cb round to r's corner before; r2 then on
                    // from r's corner after round to ca.
                    let lpv = w.loops[lp].clone();
                    let m = lpv.len();
                    let pos = |v: usize| lpv.iter().position(|&x| x == v).expect("in the loop");
                    let (pr, pb) = (pos(r), pos(cb));
                    let mut first = vec![id1];
                    let mut k = pb;
                    while k != pr {
                        first.push(lpv[k]);
                        k = (k + 1) % m;
                    }
                    let mut second = vec![id2];
                    let mut k = (pr + 1) % m;
                    while k != pb {
                        second.push(lpv[k]);
                        k = (k + 1) % m;
                    }
                    w.loops.remove(lp);
                    for part in [first, second] {
                        let part = w.without_spikes(part, t);
                        if part.len() <= 2 {
                            w.close_loop(&part, t);
                            for &v in &part {
                                w.end(v, t);
                            }
                        } else {
                            w.loops.push(part);
                        }
                    }
                } else {
                    // Two loops in one (a hole's: not met with outlines of one face).
                    return Err(Error::Failed("a drafted outline running into another".into()));
                }
            }
        }
    }
    // At the top: what is left.
    for lp in std::mem::take(&mut w.loops) {
        w.close_loop(&lp, until);
        let top: Vec<P2> = lp.iter().map(|&c| w.corners[c].position(until)).collect();
        for &c in &lp {
            w.end(c, until);
        }
        w.tops.push(top);
    }
    Ok((w.sides, w.tops))
}

/// The pieces of a side's boundary chained into one closed polygon (pieces of no length
/// left out); `None` where they do not close.
fn chained(pieces: &[Segment], tol: f64) -> Option<Vec<(P2, f64)>> {
    let close = |a: (P2, f64), b: (P2, f64)| length(sub(a.0, b.0)) <= tol && (a.1 - b.1).abs() <= tol;
    let mut left: Vec<Segment> = pieces.iter().copied().filter(|s| !close(s.0, s.1)).collect();
    if left.is_empty() {
        return None;
    }
    let first = left.remove(0);
    let mut out = vec![first.0, first.1];
    while !left.is_empty() {
        let end = *out.last().expect("points");
        let k = left.iter().position(|s| close(s.0, end) || close(s.1, end))?;
        let s = left.remove(k);
        let next = if close(s.0, end) { s.1 } else { s.0 };
        out.push(next);
    }
    if !close(out[0], *out.last().expect("points")) {
        return None;
    }
    out.pop();
    Some(out)
}

/// The drafted extrusion of the straight-edged outline `corners` (3D, in order, closed)
/// along `v`, its sides leaning in by the offset `d` at the top (out for a negative one);
/// `frame` the outline's plane's origin and axes, the outline counter-clockwise in it.
pub(super) fn drafted_by_skeleton(model: &mut Model, profile: Shape, corners: &[Point], origin: Point, x: Vector, y: Vector, v: Vector, d: f64) -> Result<Shape> {
    let flat: Vec<P2> = corners.iter().map(|p| [(*p - origin).dot(x), (*p - origin).dot(y)]).collect();
    let (sides, tops) = skeleton(&flat, d)?;
    let until = d.abs();
    let lift = |(p, t): (P2, f64)| origin + x * p[0] + y * p[1] + v * (t / until);
    let size = flat.iter().map(|p| length(sub(*p, flat[0]))).fold(0.0, f64::max).max(1e-9);
    let tol = 1e-9 * size;
    let polygon = |model: &mut Model, points: &[Point]| -> Result<Shape> {
        let vertices: Vec<Shape> = points.iter().map(|&p| make::vertex(model, p)).collect();
        let edges = (0..vertices.len()).map(|k| make::edge(model, vertices[k], vertices[(k + 1) % vertices.len()])).collect::<Result<Vec<_>>>()?;
        let wire = make::wire(model, &edges)?;
        make::planar_face(model, wire)
    };
    let mut faces = vec![profile];
    for pieces in &sides {
        let ring = chained(pieces, tol).ok_or_else(|| Error::Failed("a drafted side whose boundary does not close".into()))?;
        let points: Vec<Point> = ring.into_iter().map(lift).collect();
        faces.push(polygon(model, &points)?);
    }
    for top in &tops {
        let points: Vec<Point> = top.iter().map(|&p| lift((p, until))).collect();
        // (A top closed to a line or a point has no face.)
        let area: f64 = (0..top.len()).map(|k| cross(top[k], top[(k + 1) % top.len()])).sum::<f64>().abs();
        if area > tol * size {
            faces.push(polygon(model, &points)?);
        }
    }
    let sewn = sew::sew(model, &faces, &SewOptions { tolerance: 1e-7 * size + crate::tolerance::LINEAR, ..Default::default() })?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let [shell] = shells[..] else {
        return Err(Error::Failed(format!("a drafted extrusion's faces sewn into {} shells, not one", shells.len())));
    };
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    if properties::volume(model, solid).mass < 0.0 {
        let turned = model.new_empty(NodeKind::Solid);
        model.add(turned, shell.reversed())?;
        return Ok(turned);
    }
    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 2 x 1 rectangle drafted in by 0.5 over a height of 1: a hip roof with a ridge 1 long,
    /// V = h · w · (2 l + l') / 6 with l' the ridge.
    #[test]
    fn a_rectangle_drafted_to_its_ridge() {
        let mut model = Model::new();
        let pts = [Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(2.0, 1.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let vs: Vec<Shape> = pts.iter().map(|&p| make::vertex(&mut model, p)).collect();
        let es: Vec<Shape> = (0..4).map(|k| make::edge(&mut model, vs[k], vs[(k + 1) % 4]).unwrap()).collect();
        let wire = make::wire(&mut model, &es).unwrap();
        let face = make::planar_face(&mut model, wire).unwrap();
        let h = 1.0;
        let solid = drafted_by_skeleton(&mut model, face, &pts, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, h), 0.5).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid(), "{:?}", crate::validity::check(&model, solid));
        let volume = properties::volume(&model, solid).mass;
        let expected = h * 1.0 * (2.0 * 2.0 + 1.0) / 6.0;
        assert!((volume - expected).abs() < 1e-9, "{volume} vs {expected}");
    }

    /// Two 4 x 4 squares joined by a neck 1 wide (its two sides a little apart along it),
    /// drafted in by 0.8: the neck closes at 0.5, a reflex corner running into the edge
    /// across it (the outline splits in two), and the top is two squares, each 2.4 wide.
    #[test]
    fn a_dumbbell_drafted_past_its_neck_closing() {
        let mut model = Model::new();
        let pts = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(4.0, 1.5, 0.0),
            Point::new(6.0, 1.5, 0.0),
            Point::new(6.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 4.0, 0.0),
            Point::new(6.2, 4.0, 0.0),
            Point::new(6.2, 2.5, 0.0),
            Point::new(4.3, 2.5, 0.0),
            Point::new(4.3, 4.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
        ];
        let n = pts.len();
        let vs: Vec<Shape> = pts.iter().map(|&p| make::vertex(&mut model, p)).collect();
        let es: Vec<Shape> = (0..n).map(|k| make::edge(&mut model, vs[k], vs[(k + 1) % n]).unwrap()).collect();
        let wire = make::wire(&mut model, &es).unwrap();
        let face = make::planar_face(&mut model, wire).unwrap();
        let solid = drafted_by_skeleton(&mut model, face, &pts, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, 1.0), 0.8).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid(), "{:?}", crate::validity::check(&model, solid));
        let tops: Vec<f64> = model.explore(solid, NodeKind::Face).filter_map(|f| {
            let a = properties::area(&model, f);
            ((a.center_of_mass.z - 1.0).abs() < 1e-9).then_some(a.mass)
        }).collect();
        assert_eq!(tops.len(), 2, "{tops:?}");
        // (Each blob's top: the square offset by 0.8, its side toward the neck notched.)
        assert!(tops.iter().all(|&a| a > 2.3 * 2.3 && a < 2.5 * 2.5), "{tops:?}");
        // Between the offsets' areas at the top and the bottom.
        let volume = properties::volume(&model, solid).mass;
        assert!(volume > 2.0 * 2.4 * 2.4 && volume < 32.0 + 2.0, "{volume}");
    }

    /// An L of two arms 1 wide drafted in until both close (its reflex corner running in
    /// along the valley between them): valid, its volume the integral of the offsets' areas.
    #[test]
    fn an_l_drafted_past_an_arm_closing() {
        let mut model = Model::new();
        let pts = [Point::new(0.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0), Point::new(4.0, 1.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(1.0, 3.0, 0.0), Point::new(0.0, 3.0, 0.0)];
        let vs: Vec<Shape> = pts.iter().map(|&p| make::vertex(&mut model, p)).collect();
        let es: Vec<Shape> = (0..6).map(|k| make::edge(&mut model, vs[k], vs[(k + 1) % 6]).unwrap()).collect();
        let wire = make::wire(&mut model, &es).unwrap();
        let face = make::planar_face(&mut model, wire).unwrap();
        let solid = drafted_by_skeleton(&mut model, face, &pts, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, 1.0), 0.5).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid(), "{:?}", crate::validity::check(&model, solid));
        // Both arms 1 wide: at offset t the L is (1 - 2t) wide, (6 - 2t) long in all, until both
        // close at t = 0.5, the top: V = 2 ∫ (1 - 2t)(6 - 2t) dt over [0, 0.5] = 17/6.
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - 17.0 / 6.0).abs() < 1e-9, "{volume}");
    }
}
