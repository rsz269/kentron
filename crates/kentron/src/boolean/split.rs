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
    /// For an arc of the face's own boundary that bounds it, whether the face runs along the
    /// arc (its region on the arc's left as it goes) or against it.
    pub along: Option<bool>,
}

/// The parameters of the point of `surface` nearest `p` by Gauss-Newton steps from `guess`
/// (within the surface's bounds).
fn foot_from(surface: &crate::geometry::surface::Surface, p: crate::math::Point, guess: Uv) -> Uv {
    use crate::geometry::ParametricSurface;
    let (u0, u1, v0, v1) = surface.bounds();
    let [mut u, mut v] = guess;
    for _ in 0..30 {
        let (q, su, sv) = surface.d1(u, v);
        let d = p - q;
        let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a * c - b * b;
        if !(det.abs() > 1e-300) {
            break;
        }
        let (ru, rv) = (d.dot(su), d.dot(sv));
        let (du, dv) = ((ru * c - rv * b) / det, (rv * a - ru * b) / det);
        (u, v) = ((u + du).clamp(u0.min(u1), u0.max(u1)), (v + dv).clamp(v0.min(v1), v0.max(v1)));
        if du.abs() + dv.abs() <= 1e-15 * (1.0 + u.abs() + v.abs()) {
            break;
        }
    }
    [u, v]
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
            // (Two sections, traced where different faces meet along one curve, within a
            // little of their length: one curve.)
            let chord = self.vertices.points[p.start].distance(self.vertices.points[p.end]);
            let found = candidates.iter().copied().find(|&j| {
                let q = &pieces[j];
                let Some(curve) = &self.edges[q.edge].curve else { return false };
                let sections = self.edges[p.edge].operand.is_none() && self.edges[q.edge].operand.is_none();
                let within = if sections { on_tol().max(2e-4 * chord) } else { on_tol() };
                [0.5, 0.25, 0.75].iter().all(|&s| geom::project(curve, q.first.min(q.last), q.first.max(q.last), probe(p, s)).1 <= within)
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
        // (The first point too, on a seam of a surface closed without a period: by the next.)
        if n >= 2 {
            uv[0] = face.unwrap(uv[0], uv[1]);
        }
        // (Near such a seam a point's foot may be found on the wrong side of it, on the bound
        // there: found again from its neighbour's parameters, where they are well away from it.)
        if n >= 2 && face.closed_bounds.iter().any(Option::is_some) {
            for pass in 0..2 {
                let order: Vec<usize> = if pass == 0 { (1..n).collect() } else { (0..n - 1).rev().collect() };
                for i in order {
                    let j = if pass == 0 { i - 1 } else { i + 1 };
                    let stuck = face.closed_bounds.iter().enumerate().any(|(k, b)| {
                        b.is_some_and(|(lo, hi)| {
                            let at = |x: f64| (uv[i][k] - x).abs() <= 1e-9 * (hi - lo).abs();
                            (at(lo) || at(hi)) && (uv[j][k] - uv[i][k]).abs() > 0.25 * (hi - lo).abs()
                        })
                    });
                    if stuck {
                        uv[i] = foot_from(&face.surface, points[i], uv[j]);
                    }
                }
            }
        }
        let mid = uv[uv.len() / 2];
        let mut shift = [face.shift_u(mid[0]) - mid[0], face.shift_v(mid[1]) - mid[1]];
        // (A face spanning more than a period, as a band winding round a cylinder, holds
        // copies of a point in several: the arc goes where its middle is inside the face.)
        let wider = |period: Option<f64>, (lo, hi): (f64, f64)| period.is_some_and(|p| hi - lo > p * (1.0 + 1e-9));
        if wider(face.u_period, face.u_range) || wider(face.v_period, face.v_range) {
            // (Between its middle samples: of an arc of two, its ends are on the face's
            // boundary, in none of the copies.)
            let k = uv.len() / 2;
            let mid = if uv.len() % 2 == 0 { [0.5 * (uv[k - 1][0] + uv[k][0]), 0.5 * (uv[k - 1][1] + uv[k][1])] } else { mid };
            if let Some(q) = face.representatives(mid).find(|&q| face.inside_polylines(q)) {
                shift = [q[0] - mid[0], q[1] - mid[1]];
            }
        }
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
                    arcs.push(ArcOnFace { piece: pi, uv: w.uv.clone(), pcurve: Some(b.pcurve.clone()), leave_start: w.leave_start, leave_end: w.leave_end, along: b.bounds.then_some(b.along) });
                    continue;
                }
                arcs.push(self.boundary_arc(fi, b, pi, p.first, p.last));
            }
        }
        // Each piece once (coincident pieces of the two operands are one): those of the
        // boundary first.
        let mut reps: crate::hashing::FastSet<usize> = arcs.iter().map(|a: &ArcOnFace| pieces[a.piece].rep).collect();
        // (Not a piece of a section followed past a free-form face's surface's bounds, on the
        // surface carried on beyond them, its points (those within it) off the surface itself
        // (their feet held on its bounds) by more than the face's edges may stray, and by a tenth of its
        // length (a fitted section strays a little between its fit points): not on the face.)
        let reach = (2.0 * face.drift).max(4.0 * geom::on_tol());
        let freeform = geom::freeform(&face.surface);
        let off_surface = |p: &Piece| {
            let edge = &self.edges[p.edge];
            let length = edge.point(p.first, &self.vertices).distance(edge.point(p.last, &self.vertices));
            freeform && [0.25, 0.5, 0.75].iter().all(|&s| {
                let x = edge.point(p.first + (p.last - p.first) * s, &self.vertices);
                geom::surface_point(&face.surface, geom::uv_of(&face.surface, x)).distance(x) > reach.max(0.1 * length)
            })
        };
        for &e in &face.sections {
            for &pi in &of_edge[e] {
                if !off_surface(&pieces[pi]) && reps.insert(pieces[pi].rep) {
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
                // (Inside; or, a piece nearly along the boundary (as an arc nearly tangent to
                // it, cutting a sliver off the face), on the face's side of the boundary edge
                // nearest it, off the edge itself.)
                let at = face.locate(mid, &self.edges, &self.vertices);
                let sliver = || {
                    let reach = 1e-3 * face.extent.0.distance(face.extent.1);
                    !face.near_boundary(mid, geom::TOL, &self.edges, &self.vertices) && face.inside_by_nearest_edge(mid, reach, &self.edges) == Some(true)
                };
                if at == Location::Inside || sliver() {
                    reps.insert(p.rep);
                    arcs.push(self.arc_from_3d(fi, p, pi));
                }
            }
        }
        self.close_open_ends(fi, &mut arcs, pieces, &mut reps);
        self.split_degenerate_arcs(&mut arcs, pieces);
        arcs
    }

    /// Adds to the arcs of face `fi` the pieces of the other operand's edges from an end
    /// left open (that of a section stopping inside the face, met by no other arc) to a vertex of its other arcs
    /// (its boundary's, sections', or those of the other operand's edges on it),
    /// lying on its surface: where an edge crosses the band in which a surface the face meets
    /// nearly tangent is within tolerance of it (as a rim crossing the line where a fillet
    /// meets the face it rolls on), the sections of both stop where the band begins, and the
    /// stretch of the edge across it, found on neither, is where they meet.
    fn close_open_ends(&self, fi: usize, arcs: &mut Vec<ArcOnFace>, pieces: &[Piece], reps: &mut HashSet<usize>) {
        let face = &self.faces[fi];
        let mut ends: HashMap<usize, usize> = HashMap::default();
        for a in arcs.iter() {
            let p = &pieces[a.piece];
            *ends.entry(p.start).or_default() += 1;
            *ends.entry(p.end).or_default() += 1;
        }
        // Ends of sections met by no other arc, and the vertices of the face's boundary.
        let open: HashSet<usize> = arcs
            .iter()
            .map(|a| &pieces[a.piece])
            .filter(|p| self.edges[p.edge].operand.is_none())
            .flat_map(|p| [p.start, p.end])
            .filter(|v| ends[v] == 1 && face.vertices_on.contains(v))
            .collect();
        if open.is_empty() {
            return;
        }
        let mut found = Vec::new();
        for (pi, p) in pieces.iter().enumerate() {
            let e = &self.edges[p.edge];
            if e.operand != Some(1 - face.operand) || e.curve.is_none() || p.start == p.end || reps.contains(&p.rep) {
                continue;
            }
            if !(open.contains(&p.start) && ends.contains_key(&p.end) || open.contains(&p.end) && ends.contains_key(&p.start)) {
                continue;
            }
            // (Not a piece a few tolerances long: points the merging of vertices missed, not a
            // stretch across such a band, as wide as the root of the tolerance times the
            // surface's radius.)
            if self.vertices.points[p.start].distance(self.vertices.points[p.end]) < 100.0 * on_tol() {
                continue;
            }
            let on = (1..4).all(|k| {
                let q = e.point(p.first + (p.last - p.first) * k as f64 / 4.0, &self.vertices);
                geom::signed_distance(&face.surface, q).abs() <= on_tol() && face.locate(q, &self.edges, &self.vertices) != Location::Outside
            });
            if on {
                found.push(pi);
            }
        }
        // (Only where one such piece goes on from an end: of two, neither is known to be it.)
        let mut from: HashMap<usize, usize> = HashMap::default();
        for &pi in &found {
            let p = &pieces[pi];
            *from.entry(if open.contains(&p.start) { p.start } else { p.end }).or_default() += 1;
        }
        found.retain(|&pi| {
            let p = &pieces[pi];
            from[&if open.contains(&p.start) { p.start } else { p.end }] == 1
        });
        for pi in found {
            if reps.insert(pieces[pi].rep) {
                arcs.push(self.arc_from_3d(fi, &pieces[pi], pi));
            }
        }
    }

    /// The interior point of a region of face `fi` too thin for its arcs' polylines (a strip
    /// narrower than their sag, as between two near-parallel cuts), taken again from the arcs
    /// sampled `DENSER` times as densely from their exact geometry; none where an arc's
    /// samples do not end where its polyline does (at a pole, say).
    pub fn thin_interior(&self, fi: usize, arcs: &[ArcOnFace], pieces: &[Piece], region: &Region) -> Option<Uv> {
        const DENSER: usize = 16;
        let face = &self.faces[fi];
        let mut polygons = Vec::with_capacity(region.loops.len());
        for lp in &region.loops {
            let mut points: Vec<Uv> = Vec::new();
            for &(a, forward) in lp {
                let arc = &arcs[a];
                let piece = &pieces[arc.piece];
                let (start, end) = (arc.uv[0], *arc.uv.last()?);
                let at = |t: f64| -> Uv {
                    match &arc.pcurve {
                        Some(c) => {
                            let q = c.point(t);
                            [q.x, q.y]
                        }
                        None => geom::uv_of(&face.surface, self.edges[piece.edge].point(t, &self.vertices)),
                    }
                };
                let n = (arc.uv.len() - 1).max(1) * DENSER;
                let mut uv: Vec<Uv> = Vec::with_capacity(n + 1);
                for i in 0..=n {
                    let p = at(piece.first + (piece.last - piece.first) * i as f64 / n as f64);
                    uv.push(face.unwrap(p, uv.last().copied().unwrap_or(start)));
                }
                let reach = 1e-6 * (end[0] - start[0]).hypot(end[1] - start[1]).max(1.0);
                let gap = |p: Uv, q: Uv| (p[0] - q[0]).hypot(p[1] - q[1]);
                if gap(uv[0], start) > reach || gap(uv[n], end) > reach {
                    return None;
                }
                (uv[0], uv[n]) = (start, end);
                if !forward {
                    uv.reverse();
                }
                points.extend(&uv[..n]);
            }
            points.push(*points.first()?);
            polygons.push(points);
        }
        Some(geom::interior_point(&polygons))
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
        ArcOnFace { piece, uv, pcurve: Some(b.pcurve.clone()), leave_start, leave_end, along: b.bounds.then_some(b.along) }
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
                let p = e.point(t, &self.vertices);
                let uv = face.unwrap(geom::uv_of(&face.surface, p), reference);
                // (On a seam of a surface closed without a period, maybe the foot's other side:
                // found again from the reference.)
                let on_seam = face.closed_bounds.iter().enumerate().any(|(k, b)| b.is_some_and(|(lo, hi)| (uv[k] - lo).abs().min((uv[k] - hi).abs()) <= 1e-9 * (hi - lo).abs()));
                if on_seam { foot_from(&face.surface, p, reference) } else { uv }
            }
        };
        let step = self.departure_step(fi);
        let (start, end) = (uv[0], *uv.last().expect("polyline"));
        let leave_start = departure(near(start), piece.first, piece.last, step);
        let leave_end = departure(near(end), piece.last, piece.first, step);
        ArcOnFace { piece: index, uv, pcurve: None, leave_start, leave_end, along: None }
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
                    along: None,
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
    // (The same vertex at two nodes, one of them ending one arc: the arcs' ends a little
    // apart in (u, v), as edges meeting within their tolerances may be (or a section's end
    // and the boundary's pcurves, straying past the surface's bounds): one node.)
    let mut degree = vec![0usize; nodes.len()];
    for &(s, e) in &ends {
        degree[s] += 1;
        degree[e] += 1;
    }
    let near = 1e-3 * uv_scale.max(1.0);
    let mut merged: Vec<usize> = (0..nodes.len()).collect();
    for i in 0..nodes.len() {
        if degree[i] != 1 {
            continue;
        }
        let partners: Vec<usize> = (0..nodes.len()).filter(|&j| j != i && nodes[j].0 == nodes[i].0).collect();
        if let [j] = partners[..] {
            let (p, q) = (nodes[i].1, nodes[j].1);
            // (Of two such ends, the later into the earlier.)
            if (degree[j] > 1 || j < i) && ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt() <= near {
                merged[i] = j;
            }
        }
    }
    for e in &mut ends {
        *e = (merged[e.0], merged[e.1]);
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
        // (A node left ending one arc, by others dropped from it, is the other node of its
        // vertex near it, as above.)
        let lone = (0..nodes.len()).find_map(|i| {
            if degree[i] != 1 {
                return None;
            }
            let partners: Vec<usize> = (0..nodes.len()).filter(|&j| j != i && degree[j] > 0 && nodes[j].0 == nodes[i].0).collect();
            let (p, q) = (nodes[i].1, nodes[*partners.first()?].1);
            (partners.len() == 1 && ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt() <= near).then(|| (i, partners[0]))
        });
        if let Some((i, j)) = lone {
            for e in &mut ends {
                *e = (if e.0 == i { j } else { e.0 }, if e.1 == i { j } else { e.1 });
            }
            continue;
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
        // (Not a region whose outer loop runs along the same arcs, the other way round: the
        // cycle round the outside of that loop, its point taken off it where the arcs' ends
        // are a little apart.)
        let along = |r: &Region| r.loops[0].iter().any(|&(a, _)| cycle.iter().any(|&(b, _)| a == b));
        let container = regions
            .iter()
            .enumerate()
            .filter(|&(i, _)| in_box(i))
            .filter(|(_, r)| !along(r))
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
    // Within the face: its point inside the face's boundary, or, a sliver along the boundary
    // (thinner than its sampling), running along an arc of that boundary the way the face
    // does, the face on its left.
    let sliver = |r: &Region| {
        let polygon = &r.polygons[0];
        let perimeter: f64 = polygon.windows(2).map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt()).sum();
        2.0 * geom::signed_area(polygon).abs() <= 1e-5 * uv_scale.max(1.0) * perimeter
    };
    // (The face's boundary indexed where many regions are tested against a long one.)
    let segments: usize = face_polylines.iter().map(Vec::len).sum();
    let index = (regions.len() > 2 && segments >= 64).then(|| geom::InsidePolylines::new(face_polylines));
    let inside = |p: Uv| index.as_ref().map_or_else(|| geom::inside_polylines(face_polylines, p), |index| index.contains(p));
    regions.retain(|r| inside(r.interior) || sliver(r) && r.loops.iter().flatten().any(|&(a, forward)| arcs[a].along == Some(forward)));
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
        untangle_tangent(list, arcs);
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

/// Orders again the half-arcs leaving a node (sorted by their departures' angles) whose
/// departures are nearly the same, by their directions farther along: arcs leaving nearly
/// tangent to one another (as a section meeting a boundary where its surfaces are tangent
/// along it) curve apart by less, a step from the node, than the noise in their points.
fn untangle_tangent(list: &mut [(f64, usize, bool)], arcs: &[ArcOnFace]) {
    use std::f64::consts::{PI, TAU};
    const TANGENT: f64 = 1e-3;
    let n = list.len();
    // (Two half-arcs are in the one order round a node either way.)
    if n < 3 {
        return;
    }
    let gap = |k: usize| if k + 1 < n { list[k + 1].0 - list[k].0 } else { list[0].0 + TAU - list[n - 1].0 };
    let Some(wide) = (0..n).find(|&k| gap(k) >= TANGENT) else { return };
    if (0..n).all(|k| gap(k) >= TANGENT) {
        return;
    }
    // From past a wide gap round, the angles increasing.
    list.rotate_left((wide + 1) % n);
    let mut unwrapped: Vec<f64> = list.iter().map(|h| h.0).collect();
    for k in 1..n {
        while unwrapped[k] < unwrapped[k - 1] {
            unwrapped[k] += TAU;
        }
    }
    let polyline = |&(_, i, forward): &(f64, usize, bool)| -> Vec<Uv> {
        let mut uv = arcs[i].uv.clone();
        if !forward {
            uv.reverse();
        }
        uv
    };
    let extent = |uv: &[Uv]| uv.iter().map(|p| (p[0] - uv[0][0]).hypot(p[1] - uv[0][1])).fold(0.0, f64::max);
    let mut k = 0;
    while k < n {
        let mut end = k + 1;
        while end < n && unwrapped[end] - unwrapped[end - 1] < TANGENT {
            end += 1;
        }
        if end - k >= 2 {
            let base = unwrapped[k];
            let lines: Vec<Vec<Uv>> = list[k..end].iter().map(polyline).collect();
            let reach = 0.25 * lines.iter().map(|uv| extent(uv)).fold(f64::INFINITY, f64::min);
            // Each one's turn from the first's departure, at the first of its points that far
            // from the node.
            let turn = |uv: &[Uv]| -> f64 {
                let o = uv[0];
                let p = uv.iter().find(|p| (p[0] - o[0]).hypot(p[1] - o[1]) >= reach).unwrap_or(&uv[uv.len() - 1]);
                let a = (p[1] - o[1]).atan2(p[0] - o[0]) - base;
                a - TAU * ((a + PI) / TAU).floor()
            };
            if reach > 0.0 {
                let turns: Vec<f64> = lines.iter().map(|uv| turn(uv)).collect();
                let mut order: Vec<usize> = (0..end - k).collect();
                order.sort_by(|&a, &b| turns[a].total_cmp(&turns[b]));
                let group: Vec<(f64, usize, bool)> = order.iter().map(|&j| list[k + j]).collect();
                list[k..end].copy_from_slice(&group);
            }
        }
        k = end;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Arcs between vertices (by index), straight through the given (u, v) points.
    fn arcs(paths: &[(usize, usize, Vec<Uv>)]) -> (Vec<ArcOnFace>, Vec<Piece>) {
        let towards = |a: Uv, b: Uv| -> Uv {
            let (du, dv) = (b[0] - a[0], b[1] - a[1]);
            let l = (du * du + dv * dv).sqrt();
            [du / l, dv / l]
        };
        let mut out = (Vec::new(), Vec::new());
        for (i, (start, end, uv)) in paths.iter().enumerate() {
            let n = uv.len();
            out.1.push(Piece { edge: i, first: 0.0, last: 1.0, start: *start, end: *end, rep: i });
            out.0.push(ArcOnFace { piece: i, uv: uv.clone(), pcurve: None, leave_start: towards(uv[0], uv[1]), leave_end: towards(uv[n - 1], uv[n - 2]), along: None });
        }
        out
    }

    /// A square face with a square hole, its outer edges' ends a little apart in (u, v) at
    /// the vertices they share (as a fillet's edges may end, within their tolerances): one
    /// region with the hole, not the hole's inside, nor the outer loop again as a hole.
    #[test]
    fn a_face_whose_edges_end_a_little_apart() {
        let gap = 3e-5;
        let corners: [Uv; 4] = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]];
        let mut paths = Vec::new();
        for i in 0..4 {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            // (Each edge ends a little inside the corner where the next starts.)
            let b = [b[0] + gap * (5.0 - b[0]).signum(), b[1] + gap * (5.0 - b[1]).signum()];
            paths.push((i, (i + 1) % 4, vec![a, [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])], b]));
        }
        let hole = [[4.0, 4.0], [4.0, 6.0], [6.0, 6.0], [6.0, 4.0]];
        for i in 0..4 {
            paths.push((4 + i, 4 + (i + 1) % 4, vec![hole[i], hole[(i + 1) % 4]]));
        }
        let (arcs, pieces) = arcs(&paths);
        let mut outer: Vec<Uv> = corners.to_vec();
        outer.push(corners[0]);
        let mut inner: Vec<Uv> = hole.to_vec();
        inner.push(hole[0]);
        let regions = split_face(&arcs, &pieces, &[outer, inner], 10.0, false);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].loops.len(), 2);
    }


    /// The boundary of the unit square, its left and right sides split at mid-height
    /// (vertices 4 and 5), as paths for [`arcs`].
    fn square_split_at_mid_height() -> (Vec<(usize, usize, Vec<Uv>)>, Vec<Vec<Uv>>) {
        let paths = vec![
            (0, 1, vec![[0.0, 0.0], [1.0, 0.0]]),
            (1, 4, vec![[1.0, 0.0], [1.0, 0.5]]),
            (4, 2, vec![[1.0, 0.5], [1.0, 1.0]]),
            (2, 3, vec![[1.0, 1.0], [0.0, 1.0]]),
            (3, 5, vec![[0.0, 1.0], [0.0, 0.5]]),
            (5, 0, vec![[0.0, 0.5], [0.0, 0.0]]),
        ];
        (paths, vec![vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]])
    }

    /// A cut across a face ending a little off the boundary's vertices in (u, v) (a
    /// section's ends on the surface's bounds, the boundary's pcurves straying past them):
    /// its ends are those vertices, and it splits the face in two.
    #[test]
    fn a_cut_ending_a_little_off_the_boundary() {
        let (mut paths, boundary) = square_split_at_mid_height();
        paths.push((4, 5, vec![[1.0 - 3e-6, 0.5], [0.5, 0.52], [3e-6, 0.5]]));
        let (arcs, pieces) = arcs(&paths);
        assert_eq!(split_face(&arcs, &pieces, &boundary, 1.0, false).len(), 2);
    }

    /// The same, the cut going on from one end along a stretch that ends nowhere (dropped,
    /// it leaves that end the cut's alone): still two regions.
    #[test]
    fn a_cut_ending_off_the_boundary_where_a_stray_piece_goes_on() {
        let (mut paths, boundary) = square_split_at_mid_height();
        paths.push((4, 5, vec![[1.0 - 3e-6, 0.5], [0.5, 0.52], [3e-6, 0.5]]));
        paths.push((5, 6, vec![[3e-6, 0.5], [3e-3, 0.4]]));
        let (arcs, pieces) = arcs(&paths);
        assert_eq!(split_face(&arcs, &pieces, &boundary, 1.0, false).len(), 2);
    }

    /// A cut leaving a vertex of the boundary tangent to it (a section where the surfaces
    /// are tangent along the boundary), its departure a hair the wrong side of it (by less
    /// than the noise in its points): ordered round the vertex by where it goes, it splits
    /// the face in two.
    #[test]
    fn a_cut_leaving_the_boundary_tangent() {
        let paths = vec![
            (0, 4, vec![[0.0, 0.0], [2.0, 0.0]]),
            (4, 1, vec![[2.0, 0.0], [4.0, 0.0]]),
            (1, 2, vec![[4.0, 0.0], [4.0, 1.0]]),
            (2, 3, vec![[4.0, 1.0], [0.0, 1.0]]),
            (3, 5, vec![[0.0, 1.0], [0.0, 0.5]]),
            (5, 0, vec![[0.0, 0.5], [0.0, 0.0]]),
            (4, 5, vec![[2.0, 0.0], [1.5, 0.03], [1.0, 0.12], [0.5, 0.28], [0.0, 0.5]]),
        ];
        let (mut arcs, pieces) = arcs(&paths);
        arcs[6].leave_start = [-1.0, -1e-5];
        let boundary = vec![vec![[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [0.0, 1.0], [0.0, 0.0]]];
        let regions = split_face(&arcs, &pieces, &boundary, 4.0, false);
        assert_eq!(regions.len(), 2);
    }
}
