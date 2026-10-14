//! Splitting edges at their split points into pieces, and faces along the pieces lying
//! on them into regions, traced in each face's parameter space.

use crate::hashing::{FastMap as HashMap, FastSet as HashSet};

use crate::geometry::curve2::Curve2;
use crate::geometry::ParametricCurve;

use super::geom::{self, Uv, on_tol};
use super::graph::{Graph, Location};

/// A stretch of an edge between consecutive split points.
pub(super) struct Piece {
    pub edge: usize,
    pub first: f64,
    pub last: f64,
    pub start: usize,
    pub end: usize,
    /// The piece standing for this one in the result (pieces of the two operands that
    /// coincide are one edge there).
    pub rep: usize,
}

/// A piece as it lies on a face: its polyline in the face's (u, v), from its start to its
/// end vertex, and its pcurve there when known from the operand (same parametrization).
pub(super) struct ArcOnFace {
    pub piece: usize,
    pub uv: Vec<Uv>,
    pub pcurve: Option<Curve2>,
    /// The directions in which the arc leaves its start and its end (backwards), measured
    /// a short fixed distance along it: arcs leaving a node tangentially (touching circles)
    /// are told apart by how they curve away.
    pub leave_start: Uv,
    pub leave_end: Uv,
}

/// A boundary arc of a face over its whole edge, kept with the face (see
/// [`Graph::whole_boundary_arcs`]).
pub(super) struct WholeArc {
    pub range: (f64, f64),
    pub uv: Vec<Uv>,
    pub leave_start: Uv,
    pub leave_end: Uv,
}

/// The direction from `uv_at(t_end)` to the point a distance `step` along the arc (towards
/// `t_other`), from the exact parametrization.
fn departure(uv_at: impl Fn(f64) -> Uv, t_end: f64, t_other: f64, step: f64) -> Uv {
    let p0 = uv_at(t_end);
    let span = t_other - t_end;
    // The parameter speed near the end, then the parameter giving the wanted distance.
    let h = 1e-6 * span;
    let q = uv_at(t_end + h);
    let speed = ((q[0] - p0[0]).hypot(q[1] - p0[1]) / h.abs()).max(1e-300);
    let dt = (step / speed).min(0.25 * span.abs()) * span.signum();
    let p = uv_at(t_end + dt);
    [p[0] - p0[0], p[1] - p0[1]]
}

/// A region of a split face: its loops as (arc, forward) sequences, the first one outer
/// (counter-clockwise in (u, v)), and their polygons.
#[derive(Clone)]
pub(super) struct Region {
    pub loops: Vec<Vec<(usize, bool)>>,
    pub polygons: Vec<Vec<Uv>>,
    pub interior: Uv,
}

impl Region {
    /// The region of a face kept whole without being split (its loops and polygons are the
    /// face's own; see `Kept::far`).
    pub fn whole() -> Region {
        Region { loops: Vec::new(), polygons: Vec::new(), interior: [f64::NAN; 2] }
    }
}

impl Graph {
    /// Splits every edge at its split points; `pieces[i]` lists edge `i`'s pieces in order.
    pub fn pieces(&mut self) -> (Vec<Piece>, Vec<Vec<usize>>) {
        let mut pieces: Vec<Piece> = Vec::new();
        let mut of_edge = Vec::with_capacity(self.edges.len());
        for (ei, e) in self.edges.iter_mut().enumerate() {
            let mut stops = vec![(e.first, e.start)];
            let mut inner: Vec<(f64, usize)> = e.paves.clone();
            inner.sort_by(|a, b| a.0.total_cmp(&b.0));
            stops.extend(inner);
            stops.push((e.last, e.end));
            let mut list = Vec::new();
            for w in stops.windows(2) {
                let ((t0, v0), (t1, v1)) = (w[0], w[1]);
                if v0 == v1 && (t1 - t0).abs() <= 1e-9 * (e.last - e.first).abs().max(1e-12) {
                    continue;
                }
                list.push(pieces.len());
                let rep = pieces.len();
                pieces.push(Piece { edge: ei, first: t0, last: t1, start: v0, end: v1, rep });
            }
            of_edge.push(list);
        }
        // Coincident pieces (same ends, same middle point) become one.
        let mut by_ends: HashMap<(usize, usize), Vec<usize>> = HashMap::default();
        for i in 0..pieces.len() {
            let p = &pieces[i];
            if self.edges[p.edge].curve.is_none() {
                continue;
            }
            let key = (p.start.min(p.end), p.start.max(p.end));
            let candidates = by_ends.entry(key).or_default();
            let probe = |q: &Piece, s: f64| self.edges[q.edge].point(q.first + (q.last - q.first) * s, &self.vertices);
            // (Points of one on the other: the curves may be parametrized differently, as an
            // arc and a B-spline along it.)
            let found = candidates.iter().copied().find(|&j| {
                let q = &pieces[j];
                let Some(curve) = &self.edges[q.edge].curve else { return false };
                [0.5, 0.25, 0.75].iter().all(|&s| geom::project(curve, q.first.min(q.last), q.first.max(q.last), probe(p, s)).1 <= on_tol())
            });
            match found {
                Some(j) => pieces[i].rep = pieces[j].rep,
                None => candidates.push(i),
            }
        }
        (pieces, of_edge)
    }

    /// A piece's polyline on face `fi`, from 3D samples mapped to the face's (u, v), the
    /// periodic u unwrapped and brought into the face's range at the middle.
    fn piece_uv(&self, fi: usize, piece: &Piece) -> Vec<Uv> {
        let face = &self.faces[fi];
        let e = &self.edges[piece.edge];
        let ts = geom::sample_parameters(e.curve.as_ref(), piece.first, piece.last);
        let points: Vec<crate::math::Point> = match &e.curve {
            Some(c) => c.points_along(&ts),
            None => ts.iter().map(|&t| e.point(t, &self.vertices)).collect(),
        };
        // At a pole u is undefined: it is taken from the neighbouring sample.
        let at_pole: Vec<bool> = points.iter().map(|&p| geom::is_pole(&face.surface, p)).collect();
        let mut uv: Vec<Uv> = points.iter().map(|&p| geom::uv_of(&face.surface, p)).collect();
        let n = uv.len();
        for i in 0..n {
            if at_pole[i] {
                let neighbour = if i + 1 < n && !at_pole[i + 1] { i + 1 } else { i.saturating_sub(1) };
                uv[i][0] = uv[neighbour][0];
            }
        }
        for i in 1..n {
            uv[i] = face.unwrap(uv[i], uv[i - 1]);
        }
        let mid = uv[uv.len() / 2];
        let shift = [face.shift_u(mid[0]) - mid[0], face.shift_v(mid[1]) - mid[1]];
        for p in &mut uv {
            p[0] += shift[0];
            p[1] += shift[1];
        }
        uv
    }

    /// The pieces lying on face `fi`: its boundary's, the section edges', and those of the
    /// other operand's edges on its surface that are inside it.
    ///
    /// `whole`, when given, holds [`Self::whole_boundary_arcs`] of the face, used for the
    /// boundary's edges left whole.
    pub fn arcs_on_face(&self, fi: usize, pieces: &mut Vec<Piece>, of_edge: &[Vec<usize>], whole: Option<&[WholeArc]>) -> Vec<ArcOnFace> {
        let face = &self.faces[fi];
        let mut arcs = Vec::new();
        for (k, b) in face.boundary.iter().enumerate() {
            for &pi in &of_edge[b.edge] {
                let p = &pieces[pi];
                if let Some(w) = whole.and_then(|w| w.get(k)).filter(|w| w.range == (p.first, p.last)) {
                    arcs.push(ArcOnFace { piece: pi, uv: w.uv.clone(), pcurve: Some(b.pcurve.clone()), leave_start: w.leave_start, leave_end: w.leave_end });
                    continue;
                }
                arcs.push(self.boundary_arc(fi, b, pi, p.first, p.last));
            }
        }
        // Each piece once (coincident pieces of the two operands are one): those of the
        // boundary first.
        let mut reps: crate::hashing::FastSet<usize> = arcs.iter().map(|a: &ArcOnFace| pieces[a.piece].rep).collect();
        for &e in &face.sections {
            for &pi in &of_edge[e] {
                if reps.insert(pieces[pi].rep) {
                    arcs.push(self.arc_from_3d(fi, &pieces[pi], pi));
                }
            }
        }
        for &e in &face.on_edges {
            for &pi in &of_edge[e] {
                let p = &pieces[pi];
                if reps.contains(&p.rep) {
                    continue;
                }
                let mid = self.edges[e].point(0.5 * (p.first + p.last), &self.vertices);
                if face.locate(mid, &self.edges, &self.vertices) == Location::Inside {
                    reps.insert(p.rep);
                    arcs.push(self.arc_from_3d(fi, p, pi));
                }
            }
        }
        self.split_degenerate_arcs(&mut arcs, pieces);
        arcs
    }

    /// An arc of the boundary use `b` of face `fi` from `first` to `last`.
    fn boundary_arc(&self, fi: usize, b: &super::graph::BoundaryUse, piece: usize, first: f64, last: f64) -> ArcOnFace {
        let uv = geom::sample_pcurve(&b.pcurve, first, last);
        let at = |t: f64| {
            let q = b.pcurve.point(t);
            [q.x, q.y]
        };
        let step = self.departure_step(fi);
        let (leave_start, leave_end) = (departure(at, first, last, step), departure(at, last, first, step));
        ArcOnFace { piece, uv, pcurve: Some(b.pcurve.clone()), leave_start, leave_end }
    }

    /// Each boundary use of face `fi` as an arc over its whole edge (which the face keeps
    /// while the other operand does not cut it), with the edge's range.
    pub fn whole_boundary_arcs(&self, fi: usize) -> Vec<WholeArc> {
        let face = &self.faces[fi];
        face.boundary
            .iter()
            .map(|b| {
                let e = &self.edges[b.edge];
                let arc = self.boundary_arc(fi, b, 0, e.first, e.last);
                WholeArc { range: (e.first, e.last), uv: arc.uv, leave_start: arc.leave_start, leave_end: arc.leave_end }
            })
            .collect()
    }

    /// The departure step of face `fi`'s arcs (which [`Self::whole_boundary_arcs`] go by).
    pub fn arcs_step(&self, fi: usize) -> f64 {
        self.departure_step(fi)
    }

    /// The distance along arcs at which their departure directions are measured.
    fn departure_step(&self, fi: usize) -> f64 {
        let face = &self.faces[fi];
        let (lo, hi) = face.u_range;
        1e-4 * (hi - lo).abs().max(face.bbox.0.distance(face.bbox.1)).max(1e-3)
    }

    /// An arc of a piece lying on face `fi`, from its 3D geometry.
    fn arc_from_3d(&self, fi: usize, piece: &Piece, index: usize) -> ArcOnFace {
        let uv = self.piece_uv(fi, piece);
        let face = &self.faces[fi];
        let e = &self.edges[piece.edge];
        let near = |reference: Uv| {
            move |t: f64| {
                face.unwrap(geom::uv_of(&face.surface, e.point(t, &self.vertices)), reference)
            }
        };
        let step = self.departure_step(fi);
        let (start, end) = (uv[0], *uv.last().expect("polyline"));
        let leave_start = departure(near(start), piece.first, piece.last, step);
        let leave_end = departure(near(end), piece.last, piece.first, step);
        ArcOnFace { piece: index, uv, pcurve: None, leave_start, leave_end }
    }

    /// A degenerate edge (a pole) is a stretch of u in parameter space: arcs ending at its
    /// vertex split it there, into new pieces.
    fn split_degenerate_arcs(&self, arcs: &mut Vec<ArcOnFace>, pieces: &mut Vec<Piece>) {
        let degenerate: Vec<usize> = (0..arcs.len()).filter(|&i| self.edges[pieces[arcs[i].piece].edge].curve.is_none()).collect();
        for i in degenerate {
            let piece = &pieces[arcs[i].piece];
            let (vertex, edge, first, last) = (piece.start, piece.edge, piece.first, piece.last);
            let (a, b) = (arcs[i].uv[0], *arcs[i].uv.last().expect("polyline"));
            let (lo, hi) = (a[0].min(b[0]), a[0].max(b[0]));
            let eps = 1e-9 * (hi - lo).max(1.0);
            let mut cuts: Vec<f64> = Vec::new();
            for (j, arc) in arcs.iter().enumerate() {
                if j == i || self.edges[pieces[arc.piece].edge].curve.is_none() {
                    continue;
                }
                let p = &pieces[arc.piece];
                for (end, uv) in [(p.start, arc.uv[0]), (p.end, *arc.uv.last().expect("polyline"))] {
                    if end == vertex && (uv[1] - a[1]).abs() <= 1e-6 && uv[0] > lo + eps && uv[0] < hi - eps {
                        cuts.push(uv[0]);
                    }
                }
            }
            if cuts.is_empty() {
                continue;
            }
            // From the arc's start to its end, in u.
            cuts.sort_by(|x, y| if a[0] <= b[0] { x.total_cmp(y) } else { y.total_cmp(x) });
            cuts.dedup_by(|x, y| (*x - *y).abs() <= eps);
            let pcurve = arcs[i].pcurve.clone();
            let param = |u: f64| first + (last - first) * (u - a[0]) / (b[0] - a[0]);
            let mut stops = vec![a[0]];
            stops.extend(cuts);
            stops.push(b[0]);
            let mut replacement = Vec::new();
            for w in stops.windows(2) {
                let index = pieces.len();
                pieces.push(Piece { edge, first: param(w[0]), last: param(w[1]), start: vertex, end: vertex, rep: index });
                let (from, to) = ([w[0], a[1]], [w[1], a[1]]);
                replacement.push(ArcOnFace {
                    piece: index,
                    uv: vec![from, to],
                    pcurve: pcurve.clone(),
                    leave_start: [to[0] - from[0], 0.0],
                    leave_end: [from[0] - to[0], 0.0],
                });
            }
            arcs[i] = replacement.remove(0);
            arcs.extend(replacement);
        }
    }
}

/// Splits a face into the regions bounded by `arcs`, keeping those inside `face_polylines`
/// (the original boundary).
///
/// With `keep_slits`, arcs that bound nothing (dangling, or with the same region on both
/// sides) are kept as slits: loops along the arc and back, within the region holding them,
/// so that a later operation meeting them splits there (as OCCT's internal edges).
pub(super) fn split_face(arcs: &[ArcOnFace], pieces: &[Piece], face_polylines: &[Vec<Uv>], uv_scale: f64, keep_slits: bool) -> Vec<Region> {
    // Nodes: arc ends, merged when they are the same vertex at the same (u, v).
    let eps = 1e-7 * uv_scale.max(1.0);
    let mut nodes: Vec<(usize, Uv)> = Vec::new();
    let node = |vertex: usize, uv: Uv, nodes: &mut Vec<(usize, Uv)>| {
        match nodes.iter().position(|&(v, q)| v == vertex && ((q[0] - uv[0]).powi(2) + (q[1] - uv[1]).powi(2)).sqrt() <= eps.max(1e-6)) {
            Some(i) => i,
            None => {
                nodes.push((vertex, uv));
                nodes.len() - 1
            }
        }
    };
    let mut ends: Vec<(usize, usize)> = Vec::with_capacity(arcs.len());
    for a in arcs {
        let p = &pieces[a.piece];
        let start = node(p.start, a.uv[0], &mut nodes);
        let end = node(p.end, *a.uv.last().expect("polyline"), &mut nodes);
        ends.push((start, end));
    }
    // The same result edge between the same nodes twice (a piece of each operand) is one arc.
    let mut alive: Vec<bool> = vec![true; arcs.len()];
    let mut seen: HashSet<(usize, usize, usize)> = HashSet::default();
    for (i, a) in arcs.iter().enumerate() {
        let (s, e) = ends[i];
        let key = (pieces[a.piece].rep, s.min(e), s.max(e));
        if s != e && !seen.insert(key) {
            alive[i] = false;
        }
    }
    // Dangling arcs bound nothing.
    let duplicate: Vec<bool> = alive.iter().map(|a| !a).collect();
    loop {
        let mut degree = vec![0usize; nodes.len()];
        for (i, &(s, e)) in ends.iter().enumerate() {
            if alive[i] {
                degree[s] += 1;
                degree[e] += 1;
            }
        }
        let mut changed = false;
        for (i, &(s, e)) in ends.iter().enumerate() {
            if alive[i] && (degree[s] < 2 || degree[e] < 2) {
                alive[i] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let cycles = loop {
        let cycles = trace(arcs, &ends, &alive, &nodes);
        // An arc with the same cycle on both sides separates nothing.
        let mut bridges = Vec::new();
        for cycle in &cycles {
            let mut sides: HashMap<usize, usize> = HashMap::default();
            for &(a, _) in cycle {
                *sides.entry(a).or_default() += 1;
            }
            bridges.extend(sides.into_iter().filter(|&(a, n)| n == 2 && ends[a].0 != ends[a].1).map(|(a, _)| a));
        }
        if bridges.is_empty() {
            break cycles;
        }
        for a in bridges {
            alive[a] = false;
        }
    };
    let polygon = |cycle: &Vec<(usize, bool)>| -> Vec<Uv> {
        let mut points: Vec<Uv> = Vec::new();
        for &(a, forward) in cycle {
            let uv = &arcs[a].uv;
            let ordered: Vec<Uv> = if forward { uv.clone() } else { uv.iter().rev().copied().collect() };
            points.extend(&ordered[..ordered.len() - 1]);
        }
        points
    };
    let closed = |mut p: Vec<Uv>| {
        if let Some(&first) = p.first() {
            p.push(first);
        }
        p
    };
    let mut outers: Vec<(Vec<(usize, bool)>, Vec<Uv>)> = Vec::new();
    let mut holes: Vec<(Vec<(usize, bool)>, Vec<Uv>, Uv)> = Vec::new();
    for cycle in cycles {
        let poly = polygon(&cycle);
        if geom::signed_area(&poly) > 0.0 {
            outers.push((cycle, closed(poly)));
        } else {
            let probe = poly[0];
            holes.push((cycle, closed(poly), probe));
        }
    }
    if keep_slits {
        for (i, a) in arcs.iter().enumerate() {
            if alive[i] || duplicate[i] || ends[i].0 == ends[i].1 {
                continue;
            }
            let mut poly = a.uv.clone();
            poly.extend(a.uv.iter().rev().skip(1));
            // (Probed between its first two points, clear of where it meets other arcs.)
            let probe = [0.5 * (a.uv[0][0] + a.uv[1][0]), 0.5 * (a.uv[0][1] + a.uv[1][1])];
            holes.push((vec![(i, true), (i, false)], poly, probe));
        }
    }
    let mut regions: Vec<Region> = outers
        .into_iter()
        .map(|(cycle, poly)| Region { loops: vec![cycle], polygons: vec![poly], interior: [0.0, 0.0] })
        .collect();
    // The boxes of the outer loops, to rule out most regions cheaply.
    let boxes: Vec<(Uv, Uv)> = regions
        .iter()
        .map(|r| {
            r.polygons[0].iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| {
                ([lo[0].min(p[0]), lo[1].min(p[1])], [hi[0].max(p[0]), hi[1].max(p[1])])
            })
        })
        .collect();
    for (cycle, poly, probe) in holes {
        // A point of the hole's boundary strictly inside the smallest enclosing region.
        let in_box = |i: usize| {
            let (lo, hi) = boxes[i];
            probe[0] >= lo[0] - eps && probe[0] <= hi[0] + eps && probe[1] >= lo[1] - eps && probe[1] <= hi[1] + eps
        };
        let container = regions
            .iter()
            .enumerate()
            .filter(|&(i, _)| in_box(i))
            .filter(|(_, r)| geom::inside_polylines(&r.polygons[..1], probe) && !on_polyline(&r.polygons[0], probe, eps))
            .min_by(|a, b| geom::signed_area(&a.1.polygons[0]).total_cmp(&geom::signed_area(&b.1.polygons[0])))
            .map(|(i, _)| i);
        if let Some(i) = container {
            regions[i].loops.push(cycle);
            regions[i].polygons.push(poly);
        }
    }
    // Interior points keep clear of every arc on the face, including those bounding
    // nothing (they may lie on the other operand's boundary).
    let obstacles: Vec<&[Uv]> = arcs.iter().map(|a| a.uv.as_slice()).collect();
    let obstacles = geom::Obstacles::new(&obstacles);
    for r in &mut regions {
        r.interior = geom::interior_point_among(&r.polygons, &obstacles);
    }
    regions.retain(|r| geom::inside_polylines(face_polylines, r.interior));
    regions
}

fn on_polyline(poly: &[Uv], p: Uv, eps: f64) -> bool {
    poly.windows(2).any(|w| {
        let (a, b) = (w[0], w[1]);
        let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
        let len2 = ab[0] * ab[0] + ab[1] * ab[1];
        let t = if len2 > 0.0 { ((ap[0] * ab[0] + ap[1] * ab[1]) / len2).clamp(0.0, 1.0) } else { 0.0 };
        ((ap[0] - ab[0] * t).powi(2) + (ap[1] - ab[1] * t).powi(2)).sqrt() <= eps
    })
}

/// The cycles of the arrangement's faces: each half-arc is followed by the next one
/// clockwise around its end, so that bounded faces run counter-clockwise.
fn trace(arcs: &[ArcOnFace], ends: &[(usize, usize)], alive: &[bool], nodes: &[(usize, Uv)]) -> Vec<Vec<(usize, bool)>> {
    // Outgoing half-arcs at each node, by angle.
    let mut around: Vec<Vec<(f64, usize, bool)>> = vec![Vec::new(); nodes.len()];
    for (i, a) in arcs.iter().enumerate() {
        if !alive[i] {
            continue;
        }
        let (s, e) = ends[i];
        let (out, back) = (a.leave_start, a.leave_end);
        around[s].push((out[1].atan2(out[0]), i, true));
        around[e].push((back[1].atan2(back[0]), i, false));
    }
    for list in &mut around {
        list.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    let end_of = |(a, forward): (usize, bool)| if forward { ends[a].1 } else { ends[a].0 };
    let mut visited: HashSet<(usize, bool)> = HashSet::default();
    let mut cycles = Vec::new();
    for i in 0..arcs.len() {
        if !alive[i] {
            continue;
        }
        for start in [(i, true), (i, false)] {
            if visited.contains(&start) {
                continue;
            }
            let mut cycle = Vec::new();
            let mut h = start;
            loop {
                visited.insert(h);
                cycle.push(h);
                let node = end_of(h);
                let list = &around[node];
                let twin = (h.0, !h.1);
                let k = list.iter().position(|&(_, a, f)| (a, f) == twin).expect("twin half-arc");
                let (_, a, f) = list[(k + list.len() - 1) % list.len()];
                h = (a, f);
                if h == start || visited.contains(&h) {
                    break;
                }
            }
            cycles.push(cycle);
        }
    }
    cycles
}
