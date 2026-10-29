//! Sewing: joining faces whose boundaries meet within a tolerance into shells.
//!
//! Faces given separately (or shells with gaps) share no edges where they meet. Sewing finds
//! the free edges (used by one face only) that run along each other within the tolerance and
//! replaces each such pair by one shared edge:
//!
//! 1. Free edges are cut where an end of another free edge lies along them, so that a long
//!    edge facing several short ones is matched piece by piece.
//! 2. Pieces whose ends and curves lie within the tolerance of each other are paired, the
//!    closest first; each piece joins at most one other (the result is manifold).
//! 3. The faces are rebuilt: a pair becomes one edge, with the curve of one piece and, on
//!    the other piece's face, that face's own boundary reparametrized along the curve. The
//!    vertices met at the pair's ends become one, and tolerances grow to cover the gaps.
//! 4. Faces connected by shared edges make shells, their orientations made to agree.
//!
//! The result is a shell, a face (a face sewn to nothing), or a compound of those.

use std::collections::{HashMap, VecDeque};

use crate::error::{Error, Result};
use crate::geometry::curve2::BSplineCurve2;
use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::extrema::{nearest_on_curve, NearestOnCurve};
use crate::geometry::fit;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{BoundingBox, Bvh, Point, Point2, Transformable};
use crate::topology::{ByOccurrence, CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

/// How [`sew`] joins faces.
#[doc(alias = "BRepBuilderAPI_Sewing")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SewOptions {
    /// The largest gap closed: boundaries this close to each other are joined.
    pub tolerance: f64,
    /// Cut free edges where the ends of other free edges lie along them.
    #[doc(alias = "SetCutting")]
    pub cutting: bool,
    /// Let more than two faces share an edge. Not implemented.
    #[doc(alias = "SetNonManifoldMode")]
    pub non_manifold: bool,
}

impl Default for SewOptions {
    fn default() -> Self {
        SewOptions { tolerance: 1.0e-6, cutting: true, non_manifold: false }
    }
}

/// The faces of `shapes` joined along boundaries that meet within `options.tolerance`. The
/// shapes are left as they are; the result is made of new faces.
#[doc(alias = "BRepBuilderAPI_Sewing::Perform")]
#[doc(alias = "sewing")]
pub fn sew(model: &mut Model, shapes: &[Shape], options: &SewOptions) -> Result<Shape> {
    if options.non_manifold {
        todo!("non-manifold sewing");
    }
    let mut input = Input::collect(model, shapes)?;
    let tolerance = options.tolerance;
    if options.cutting {
        input.cut(tolerance);
    }
    input.split_edges();
    input.pair(tolerance);
    input.build(model, tolerance)
}

struct Vertex {
    point: Point,
    tolerance: f64,
}

struct Edge {
    /// `None` for a degenerate edge.
    curve: Option<Curve>,
    first: f64,
    last: f64,
    start: usize,
    end: usize,
    tolerance: f64,
    /// The number of face boundaries using the edge.
    uses: usize,
    /// Where the edge is cut: (parameter, vertex), in increasing parameter.
    cuts: Vec<(f64, usize)>,
    /// Its pieces, in increasing parameter.
    sections: Vec<usize>,
    /// The input edge it is.
    source: crate::topology::NodeId,
}

impl Edge {
    fn free(&self) -> bool {
        self.uses == 1 && self.curve.is_some()
    }
}

/// An edge as used by a face's boundary.
struct EdgeUse {
    edge: usize,
    /// The boundary runs along the edge's curve.
    forward: bool,
    pcurve: Curve2,
}

struct Face {
    /// The input face it is (for the history of what the sewing made of it).
    source: crate::topology::NodeId,
    surface: Surface,
    tolerance: f64,
    reversed: bool,
    wires: Vec<Vec<EdgeUse>>,
}

/// A piece of an edge between two vertices (the whole edge if it is not cut).
struct Section {
    edge: usize,
    first: f64,
    last: f64,
    start: usize,
    end: usize,
    /// The section whose edge stands for this one in the result (itself if not paired),
    /// and whether the two run the same way.
    rep: usize,
    same: bool,
}

struct Input {
    vertices: Vec<Vertex>,
    /// Union-find over `vertices`: vertices merged into one.
    parent: Vec<usize>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    sections: Vec<Section>,
}

impl Input {
    fn collect(model: &Model, shapes: &[Shape]) -> Result<Input> {
        let mut input = Input { vertices: Vec::new(), parent: Vec::new(), edges: Vec::new(), faces: Vec::new(), sections: Vec::new() };
        let mut vertex_of: ByOccurrence<usize> = ByOccurrence::default();
        let mut edge_of: ByOccurrence<usize> = ByOccurrence::default();
        let mut faces_seen: ByOccurrence<()> = ByOccurrence::default();
        let mut vertex = |input: &mut Input, v: Shape| -> usize {
            if let Some(&i) = vertex_of.get(&v) {
                return i;
            }
            input.vertices.push(Vertex { point: model.point(v), tolerance: model.tolerance(v) });
            input.parent.push(input.parent.len());
            vertex_of.insert(&v, input.parent.len() - 1);
            input.parent.len() - 1
        };
        for &shape in shapes {
            for face in model.explore(shape, NodeKind::Face) {
                if faces_seen.get(&face).is_some() {
                    continue;
                }
                faces_seen.insert(&face, ());
                let surface = model
                    .surface(face)
                    .ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?
                    .transformed(&face.placement);
                let mut wires = Vec::new();
                for wire in model.children(Shape { orientation: Orientation::Forward, ..face }) {
                    let mut uses = Vec::new();
                    for coedge in model.children(wire) {
                        let id = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                        let edge_id: EdgeId = model[id].edge;
                        let edge_shape = Shape { placement: coedge.placement, ..Shape::from(edge_id) };
                        let index = match edge_of.get(&edge_shape) {
                            Some(&i) => i,
                            None => {
                                let (a, b) = model.edge_vertices(edge_shape);
                                let (start, end) = (vertex(&mut input, a), vertex(&mut input, b));
                                let (curve, first, last) = match model.edge_curve(edge_shape) {
                                    Some((c, first, last)) => (Some(c), first, last),
                                    None => (None, 0.0, 0.0),
                                };
                                let tolerance = model.tolerance(edge_shape);
                                input.edges.push(Edge { curve, first, last, start, end, tolerance, uses: 0, cuts: Vec::new(), sections: Vec::new(), source: edge_id.into() });
                                edge_of.insert(&edge_shape, input.edges.len() - 1);
                                input.edges.len() - 1
                            }
                        };
                        let (pcurve, first, last) = match model.effective_pcurve(id) {
                            Some(found) => found,
                            None => projected_pcurve(&input.edges[index], &surface)?,
                        };
                        let edge = &mut input.edges[index];
                        if edge.curve.is_none() {
                            (edge.first, edge.last) = (first, last);
                        }
                        edge.uses += 1;
                        let orientation = model[id].orientation.compose(coedge.orientation);
                        uses.push(EdgeUse { edge: index, forward: orientation != Orientation::Reversed, pcurve });
                    }
                    // A wire used reversed runs through its coedges backwards.
                    if wire.orientation == Orientation::Reversed {
                        uses.reverse();
                    }
                    wires.push(uses);
                }
                let reversed = face.orientation == Orientation::Reversed;
                input.faces.push(Face { source: face.node, surface, tolerance: model.tolerance(face), reversed, wires });
            }
        }
        if input.faces.is_empty() {
            return Err(Error::InvalidInput("nothing to sew: no faces".into()));
        }
        Ok(input)
    }

    fn root(&mut self, mut v: usize) -> usize {
        while self.parent[v] != v {
            self.parent[v] = self.parent[self.parent[v]];
            v = self.parent[v];
        }
        v
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.root(a), self.root(b));
        if a != b {
            self.parent[a.max(b)] = a.min(b);
        }
    }

    /// Cuts each free edge where an end of another free edge lies along it (within
    /// `tolerance`, away from its own ends).
    fn cut(&mut self, tolerance: f64) {
        let free: Vec<usize> = (0..self.edges.len()).filter(|&e| self.edges[e].free()).collect();
        let boxes: Vec<BoundingBox> = free.iter().map(|&e| edge_box(&self.edges[e], tolerance)).collect();
        let tree = Bvh::new(&boxes);
        let mut ends: Vec<usize> = free.iter().flat_map(|&e| [self.edges[e].start, self.edges[e].end]).collect();
        ends.sort_unstable();
        ends.dedup();
        for v in ends {
            let p = self.vertices[v].point;
            for k in tree.containing(p) {
                let edge = &self.edges[free[k]];
                let curve = edge.curve.as_ref().expect("a free edge has a curve");
                let (a, b) = (curve.point(edge.first), curve.point(edge.last));
                if v == edge.start || v == edge.end || p.distance(a) <= tolerance || p.distance(b) <= tolerance {
                    continue;
                }
                let found = nearest_on_curve(p, curve, edge.first, edge.last);
                let inside = found.parameter > edge.first.min(edge.last) && found.parameter < edge.first.max(edge.last);
                if found.distance <= tolerance && inside {
                    self.edges[free[k]].cuts.push((found.parameter, v));
                }
            }
        }
        // Cuts closer than the tolerance to each other are one, their vertices merged.
        for e in free {
            let mut cuts = std::mem::take(&mut self.edges[e].cuts);
            if cuts.is_empty() {
                continue;
            }
            cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
            let curve = self.edges[e].curve.clone().expect("a free edge has a curve");
            let mut kept: Vec<(f64, usize)> = Vec::new();
            for (t, v) in cuts {
                match kept.last() {
                    Some(&(s, w)) if curve.point(s).distance(curve.point(t)) <= tolerance => self.join(w, v),
                    _ => kept.push((t, v)),
                }
            }
            self.edges[e].cuts = kept;
        }
    }

    /// The edges' sections, between their ends and cuts.
    fn split_edges(&mut self) {
        for e in 0..self.edges.len() {
            let edge = &self.edges[e];
            let mut stops = vec![(edge.first, edge.start)];
            stops.extend(edge.cuts.iter().copied());
            stops.push((edge.last, edge.end));
            let mut sections = Vec::with_capacity(stops.len() - 1);
            for w in stops.windows(2) {
                let rep = self.sections.len();
                self.sections.push(Section { edge: e, first: w[0].0, last: w[1].0, start: w[0].1, end: w[1].1, rep, same: true });
                sections.push(rep);
            }
            self.edges[e].sections = sections;
        }
    }

    fn section_curve(&self, s: usize) -> Option<&Curve> {
        self.edges[self.sections[s].edge].curve.as_ref()
    }

    /// Pairs free sections lying along each other within `tolerance`, the closest pairs
    /// first, merging the vertices at their ends.
    fn pair(&mut self, tolerance: f64) {
        let free: Vec<usize> = (0..self.sections.len()).filter(|&s| self.edges[self.sections[s].edge].free()).collect();
        let boxes: Vec<BoundingBox> = free.iter().map(|&s| self.section_box(s, tolerance)).collect();
        let tree = Bvh::new(&boxes);
        let mut candidates: Vec<(f64, usize, usize, bool)> = Vec::new();
        for (i, j) in tree.overlapping_pairs(&tree) {
            let (s, t) = (free[i], free[j]);
            if s >= t || self.sections[s].edge == self.sections[t].edge {
                continue;
            }
            if let Some((deviation, same)) = self.match_sections(s, t, tolerance) {
                candidates.push((deviation, s, t, same));
            }
        }
        candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut paired = vec![false; self.sections.len()];
        for (_, s, t, same) in candidates {
            if paired[s] || paired[t] {
                continue;
            }
            let (s0, s1) = (self.sections[s].start, self.sections[s].end);
            let (t0, t1) = if same { (self.sections[t].start, self.sections[t].end) } else { (self.sections[t].end, self.sections[t].start) };
            // Two edges from one vertex (consecutive edges of a boundary, folding back within
            // the tolerance) are not a seam to close. (A vertex where an edge is cut is not
            // the edge's own.)
            let own = |section: usize, v: usize| {
                let edge = &self.edges[self.sections[section].edge];
                v == edge.start || v == edge.end
            };
            if [(s0, t0), (s1, t1)].iter().any(|&(a, b)| a == b && own(s, a) && own(t, b)) {
                continue;
            }
            let [a0, a1, b0, b1] = [s0, s1, t0, t1].map(|v| self.root(v));
            // Joining must not bring a section's two ends together (unless it is closed).
            let closed = a0 == a1 && b0 == b1;
            if !closed && (a0 == a1 || a0 == b1 || b0 == a1 || b0 == b1) {
                continue;
            }
            self.join(a0, b0);
            self.join(a1, b1);
            paired[s] = true;
            paired[t] = true;
            self.sections[t].rep = s;
            self.sections[t].same = same;
        }
    }

    fn section_box(&self, s: usize, margin: f64) -> BoundingBox {
        let section = &self.sections[s];
        let curve = self.section_curve(s).expect("a free section has a curve");
        curve_box(curve, section.first, section.last, margin)
    }

    /// The largest distance between two sections whose ends meet within `tolerance` and whose
    /// curves stay within it, and whether they run the same way.
    fn match_sections(&self, s: usize, t: usize, tolerance: f64) -> Option<(f64, bool)> {
        let (a, b) = (&self.sections[s], &self.sections[t]);
        let (ca, cb) = (self.section_curve(s)?, self.section_curve(t)?);
        let (a0, a1) = (ca.point(a.first), ca.point(a.last));
        let (b0, b1) = (cb.point(b.first), cb.point(b.last));
        let same = a0.distance(b0) <= tolerance && a1.distance(b1) <= tolerance;
        let opposite = a0.distance(b1) <= tolerance && a1.distance(b0) <= tolerance;
        let same = match (same, opposite) {
            (false, false) => return None,
            (true, false) => true,
            (false, true) => false,
            // Both ends together (closed or tiny sections): by the points a quarter along.
            (true, true) => {
                let q = ca.point(lerp(a.first, a.last, 0.25));
                q.distance(cb.point(lerp(b.first, b.last, 0.25))) <= q.distance(cb.point(lerp(b.first, b.last, 0.75)))
            }
        };
        let ends = if same { a0.distance(b0).max(a1.distance(b1)) } else { a0.distance(b1).max(a1.distance(b0)) };
        let mut deviation = ends;
        for i in 1..8 {
            let p = ca.point(lerp(a.first, a.last, i as f64 / 8.0));
            deviation = deviation.max(nearest_on_curve(p, cb, b.first, b.last).distance);
            if deviation > tolerance {
                return None;
            }
        }
        let q = cb.point(lerp(b.first, b.last, 0.5));
        deviation = deviation.max(nearest_on_curve(q, ca, a.first, a.last).distance);
        (deviation <= tolerance).then_some((deviation, same))
    }

    /// The result's faces, grouped into shells.
    fn build(mut self, model: &mut Model, tolerance: f64) -> Result<Shape> {
        // Vertices: one per group of merged ones, at their middle.
        let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
        for s in &self.sections {
            for v in [s.start, s.end] {
                members.entry(v).or_default();
            }
        }
        let used: Vec<usize> = members.into_keys().collect();
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for v in used {
            let r = self.root(v);
            groups.entry(r).or_default().push(v);
        }
        let mut out_vertex: HashMap<usize, (VertexId, Point)> = HashMap::new();
        let mut vertex_tolerance: HashMap<usize, f64> = HashMap::new();
        let mut roots: Vec<usize> = groups.keys().copied().collect();
        roots.sort_unstable();
        for r in roots {
            let group = &groups[&r];
            let mut sum = crate::math::Vector::new(0.0, 0.0, 0.0);
            for &v in group {
                sum = sum + (self.vertices[v].point - Point::ORIGIN);
            }
            let point = Point::ORIGIN + sum / group.len() as f64;
            let tolerance = group.iter().map(|&v| point.distance(self.vertices[v].point) + self.vertices[v].tolerance).fold(0.0, f64::max);
            out_vertex.insert(r, (model.add_vertex(point, tolerance), point));
            vertex_tolerance.insert(r, tolerance);
        }
        // (Sections no longer than the tolerance whose ends were merged into one vertex
        // are left out, as OCCT's sewing drops the edges it makes degenerate.)
        let collapsed: Vec<bool> = (0..self.sections.len())
            .map(|s| {
                let section = &self.sections[s];
                let Some(curve) = &self.edges[section.edge].curve else { return false };
                if section.start == section.end || self.parent_root(section.start) != self.parent_root(section.end) {
                    return false;
                }
                const SAMPLES: usize = 8;
                let points: Vec<Point> = (0..=SAMPLES).map(|i| curve.point(lerp(section.first, section.last, i as f64 / SAMPLES as f64))).collect();
                points.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>() <= tolerance
            })
            .collect();
        // Edges: one per representative section.
        let mut out_edge: HashMap<usize, (EdgeId, f64)> = HashMap::new();
        for s in 0..self.sections.len() {
            if self.sections[s].rep != s || collapsed[s] {
                continue;
            }
            let section = &self.sections[s];
            let edge = &self.edges[section.edge];
            let (start, end) = (self.parent_root(section.start), self.parent_root(section.end));
            let id = model.add_edge(out_vertex[&start].0, out_vertex[&end].0, edge.tolerance)?;
            let mut tolerance = edge.tolerance;
            if let Some(curve) = &edge.curve {
                model.set_edge_curve(id.into(), curve.clone(), edge.tolerance)?;
                model.set_edge_range(id.into(), section.first, section.last)?;
                for (v, t) in [(start, section.first), (end, section.last)] {
                    let gap = out_vertex[&v].1.distance(curve.point(t));
                    let vt = vertex_tolerance.get_mut(&v).expect("vertex");
                    *vt = vt.max(gap);
                }
            }
            for other in 0..self.sections.len() {
                if self.sections[other].rep == s {
                    if other != s {
                        tolerance = tolerance.max(self.edges[self.sections[other].edge].tolerance);
                    }
                    // (Each input edge sewn into it, or a piece of one, made it.)
                    crate::boolean::record(model.uid(), self.edges[self.sections[other].edge].source, id.into());
                }
            }
            out_edge.insert(s, (id, tolerance));
        }
        // Faces, and how they use the edges: (face, forward) per edge.
        let mut uses: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();
        let mut faces = Vec::with_capacity(self.faces.len());
        for (fi, face) in self.faces.iter().enumerate() {
            let shape = model.new_empty(NodeKind::Face);
            crate::boolean::record(model.uid(), face.source, shape.node);
            model.set_face_surface(shape, face.surface.clone(), face.tolerance)?;
            for wire in &face.wires {
                let mut coedges = Vec::new();
                for u in wire {
                    let edge = &self.edges[u.edge];
                    let mut sections = edge.sections.clone();
                    if !u.forward {
                        sections.reverse();
                    }
                    for s in sections {
                        let section = &self.sections[s];
                        if collapsed[section.rep] {
                            // (Kept, as a degenerate edge of this face's own, where its pcurve
                            // still runs across the face's parameters: a ruled surface's end
                            // between curves that meet there, say. Without it the face's
                            // boundary would not close.)
                            let (a, b) = (u.pcurve.point(section.first), u.pcurve.point(section.last));
                            if a.distance(b) > 1e-9 * (1.0 + (a - Point2::ORIGIN).length()) {
                                let vertex = out_vertex[&self.parent_root(section.start)].0;
                                let degenerate = model.add_edge(vertex, vertex, edge.tolerance)?;
                                let coedge = model.add_coedge(degenerate, if u.forward { Orientation::Forward } else { Orientation::Reversed })?;
                                model.set_pcurve_with_range(coedge, u.pcurve.clone(), section.first, section.last)?;
                                coedges.push(coedge);
                            }
                            continue;
                        }
                        let (id, tolerance) = out_edge.get_mut(&section.rep).expect("edge");
                        let forward = u.forward == section.same;
                        let coedge = model.add_coedge(*id, if forward { Orientation::Forward } else { Orientation::Reversed })?;
                        let rep = &self.sections[section.rep];
                        let pcurve = if section.rep == s {
                            u.pcurve.clone()
                        } else {
                            let curve = self.edges[rep.edge].curve.as_ref().expect("a paired section has a curve");
                            let other = edge.curve.as_ref().expect("a paired section has a curve");
                            let (pcurve, deviation) = transfer(curve, (rep.first, rep.last), other, section, &u.pcurve, &face.surface)?;
                            *tolerance = tolerance.max(deviation);
                            pcurve
                        };
                        model.set_pcurve_with_range(coedge, pcurve, rep.first, rep.last)?;
                        coedges.push(coedge);
                        uses.entry(section.rep).or_default().push((fi, forward));
                    }
                }
                let wire = model.add_wire(&coedges)?;
                model.add(shape, wire.into())?;
            }
            faces.push(shape);
        }
        // Tolerances covering the gaps.
        for (&s, &(id, tolerance)) in &out_edge {
            model.set_edge_tolerance(id, tolerance)?;
            let section = &self.sections[s];
            for v in [section.start, section.end] {
                let r = self.parent_root(v);
                let vt = vertex_tolerance.get_mut(&r).expect("vertex");
                *vt = vt.max(tolerance);
            }
        }
        for (r, tolerance) in vertex_tolerance {
            let (id, point) = out_vertex[&r];
            model.set_vertex_point(id.into(), point, tolerance)?;
        }
        let shells = shells(&self.faces, &uses);
        assemble(model, &faces, shells)
    }

    /// The root of `v`'s group, without path compression (for use through `&self`).
    fn parent_root(&self, mut v: usize) -> usize {
        while self.parent[v] != v {
            v = self.parent[v];
        }
        v
    }
}

/// The pcurve, on `surface`, of the edge `curve` over `range`, from the boundary of another
/// face running along it: that face's `pcurve` over `section`, of the curve `other`,
/// reparametrized to follow `curve`. Returns the pcurve and its largest distance from `curve`.
fn transfer(curve: &Curve, range: (f64, f64), other: &Curve, section: &Section, pcurve: &Curve2, surface: &Surface) -> Result<(Curve2, f64)> {
    // (Finer where the interpolation strays measurably from the curve: an ellipse's
    // pcurve on a cylinder takes more than a few points to follow closely.)
    let size = curve.point(range.0).distance(curve.point(0.5 * (range.0 + range.1))) + curve.point(0.5 * (range.0 + range.1)).distance(curve.point(range.1));
    let mut found = transfer_with(curve, range, other, section, pcurve, surface, 16)?;
    for n in [64, 256] {
        if found.1 <= 1e-10 * (1.0 + size) {
            break;
        }
        let finer = transfer_with(curve, range, other, section, pcurve, surface, n)?;
        if finer.1 < found.1 {
            found = finer;
        }
    }
    Ok(found)
}

/// [`transfer`] by interpolating `n + 1` points.
fn transfer_with(curve: &Curve, range: (f64, f64), other: &Curve, section: &Section, pcurve: &Curve2, surface: &Surface, n: usize) -> Result<(Curve2, f64)> {
    #[allow(non_snake_case)]
    let N = n;
    let ts: Vec<f64> = (0..=N).map(|i| lerp(range.0, range.1, i as f64 / N as f64)).collect();
    let (lo, hi) = (section.first.min(section.last), section.first.max(section.last));
    let onto = NearestOnCurve::new(other, lo, hi);
    let mut ss: Vec<f64> = ts.iter().map(|&t| onto.nearest(curve.point(t)).parameter).collect();
    (ss[0], ss[N]) = if section.same { (section.first, section.last) } else { (section.last, section.first) };
    // Keep the parameters monotonic (projections may wobble where curves are close to
    // tangent to each other's normals).
    let increasing = ss[N] >= ss[0];
    for i in 1..=N {
        if (ss[i] < ss[i - 1]) == increasing && ss[i] != ss[i - 1] {
            ss[i] = ss[i - 1];
        }
    }
    let uv: Vec<Point2> = ss.iter().map(|&s| pcurve.point(s)).collect();
    let affine = (0..=N).all(|i| (ss[i] - lerp(ss[0], ss[N], i as f64 / N as f64)).abs() <= 1e-9 * (1.0 + (ss[N] - ss[0]).abs()));
    let result = if affine && matches!(pcurve, Curve2::Line(_)) {
        Curve2::BSpline(BSplineCurve2::new(vec![uv[0], uv[N]], vec![range.0, range.1], vec![2, 2], 1))
    } else {
        Curve2::BSpline(fit::interpolate_with_parameters(&uv, &ts)?)
    };
    let mut deviation: f64 = 0.0;
    for i in 0..=2 * N {
        let t = lerp(range.0, range.1, i as f64 / (2 * N) as f64);
        let uv = result.point(t);
        deviation = deviation.max(surface.point(uv.x, uv.y).distance(curve.point(t)));
    }
    Ok((result, deviation))
}

/// The pcurve of an edge on a face that has none stored: the edge's curve projected onto
/// the face's surface.
fn projected_pcurve(edge: &Edge, surface: &Surface) -> Result<(Curve2, f64, f64)> {
    let curve = edge.curve.as_ref().ok_or_else(|| Error::InvalidInput("a degenerate edge without a pcurve".into()))?;
    let (first, last) = (edge.first, edge.last);
    let pcurve = crate::geometry::project::analytic_pcurve(curve, surface)
        .or_else(|| crate::geometry::project::pcurve_on_range(curve, surface, first, last, 1e-7).map(Curve2::BSpline))
        .or_else(|| crate::geometry::project::polyline_pcurve_on_range(curve, surface, first, last, 1e-7).map(Curve2::BSpline))
        .ok_or_else(|| Error::InvalidInput("no pcurve for a face's edge".into()))?;
    Ok((pcurve, first, last))
}

/// Faces grouped by shared edges, each group's orientations made to agree (each shared edge
/// used once forward and once reversed). Returns per shell its faces with their orientation.
fn shells(faces: &[Face], uses: &HashMap<usize, Vec<(usize, bool)>>) -> Vec<Vec<(usize, bool)>> {
    // Per face: (neighbour, whether the two must have opposite flips).
    let mut links: Vec<Vec<(usize, bool)>> = vec![Vec::new(); faces.len()];
    let mut keys: Vec<usize> = uses.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        let around = &uses[&key];
        let distinct = |a: &(usize, bool), b: &(usize, bool)| a.0 != b.0;
        match &around[..] {
            [a, b] if distinct(a, b) => {
                // Same direction on both faces: one must be flipped.
                let flip = a.1 == b.1;
                links[a.0].push((b.0, flip));
                links[b.0].push((a.0, flip));
            }
            [a, rest @ ..] => {
                for b in rest.iter().filter(|b| distinct(a, b)) {
                    links[a.0].push((b.0, false));
                    links[b.0].push((a.0, false));
                }
            }
            [] => {}
        }
    }
    let mut flipped: Vec<Option<bool>> = vec![None; faces.len()];
    let mut shells = Vec::new();
    for start in 0..faces.len() {
        if flipped[start].is_some() {
            continue;
        }
        flipped[start] = Some(faces[start].reversed);
        let mut shell = Vec::new();
        let mut queue = VecDeque::from([start]);
        while let Some(f) = queue.pop_front() {
            let reversed = flipped[f].expect("visited");
            shell.push((f, reversed));
            for &(g, flip) in &links[f] {
                if flipped[g].is_none() {
                    flipped[g] = Some(reversed != flip);
                    queue.push_back(g);
                }
            }
        }
        shell.sort_unstable();
        shells.push(shell);
    }
    shells
}

/// One face per lone face, a shell per group; a compound if there are several.
fn assemble(model: &mut Model, faces: &[Shape], shells: Vec<Vec<(usize, bool)>>) -> Result<Shape> {
    let oriented = |f: usize, reversed: bool| if reversed { faces[f].reversed() } else { faces[f] };
    let mut parts = Vec::with_capacity(shells.len());
    for members in shells {
        if let [(f, reversed)] = members[..] {
            parts.push(oriented(f, reversed));
            continue;
        }
        let shell = model.new_empty(NodeKind::Shell);
        for (f, reversed) in members {
            model.add(shell, oriented(f, reversed))?;
        }
        parts.push(shell);
    }
    if parts.len() == 1 {
        return Ok(parts[0]);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for part in parts {
        model.add(compound, part)?;
    }
    Ok(compound)
}

fn lerp(a: f64, b: f64, s: f64) -> f64 {
    a + (b - a) * s
}

/// The box of a curve over `[first, last]`, enlarged by `margin`.
fn curve_box(curve: &Curve, first: f64, last: f64, margin: f64) -> BoundingBox {
    let points: Vec<Point> = (0..=32).map(|i| curve.point(lerp(first, last, i as f64 / 32.0))).collect();
    // Chords may miss a little of the curve.
    let chord = points.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max);
    BoundingBox::from_points(points).enlarged(margin + 0.1 * chord)
}

fn edge_box(edge: &Edge, margin: f64) -> BoundingBox {
    curve_box(edge.curve.as_ref().expect("a free edge has a curve"), edge.first, edge.last, margin)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::make;

    /// A face through the corners, all moved by `shift` (so that faces leave gaps).
    fn square(model: &mut Model, corners: [[f64; 3]; 4], shift: f64) -> Shape {
        let points: Vec<Point> = corners.iter().map(|c| Point::new(c[0] + shift, c[1] - 0.5 * shift, c[2] + 0.7 * shift)).collect();
        let edges: Vec<Shape> = (0..4).map(|i| make::edge_between(model, points[i], points[(i + 1) % 4]).unwrap()).collect();
        let wire = make::wire(model, &edges).unwrap();
        make::planar_face(model, wire).unwrap()
    }

    #[test]
    fn six_loose_squares_make_a_closed_cube() {
        let mut model = Model::new();
        let c = |x: f64, y: f64, z: f64| [x, y, z];
        let faces = [
            [c(0., 0., 0.), c(0., 1., 0.), c(1., 1., 0.), c(1., 0., 0.)],
            [c(0., 0., 1.), c(1., 0., 1.), c(1., 1., 1.), c(0., 1., 1.)],
            [c(0., 0., 0.), c(1., 0., 0.), c(1., 0., 1.), c(0., 0., 1.)],
            [c(0., 1., 0.), c(0., 1., 1.), c(1., 1., 1.), c(1., 1., 0.)],
            [c(0., 0., 0.), c(0., 0., 1.), c(0., 1., 1.), c(0., 1., 0.)],
            [c(1., 0., 0.), c(1., 1., 0.), c(1., 1., 1.), c(1., 0., 1.)],
        ];
        let shapes: Vec<Shape> = faces.iter().enumerate().map(|(i, &f)| square(&mut model, f, 1e-4 * i as f64)).collect();
        // Too small a tolerance sews nothing: six lone faces.
        let loose = sew(&mut model, &shapes, &SewOptions { tolerance: 1e-6, ..Default::default() }).unwrap();
        assert_eq!(loose.kind(), NodeKind::Compound);
        assert_eq!(model.unique_subshapes(loose, NodeKind::Edge).len(), 24);
        let shell = sew(&mut model, &shapes, &SewOptions { tolerance: 1e-3, ..Default::default() }).unwrap();
        assert_eq!(shell.kind(), NodeKind::Shell);
        assert_eq!(model.unique_subshapes(shell, NodeKind::Face).len(), 6);
        assert_eq!(model.unique_subshapes(shell, NodeKind::Edge).len(), 12);
        assert_eq!(model.unique_subshapes(shell, NodeKind::Vertex).len(), 8);
        // Consistently oriented: every edge used once each way.
        let mut uses: HashMap<crate::topology::NodeId, i32> = HashMap::new();
        for face in model.explore(shell, NodeKind::Face) {
            for edge in model.explore(face, NodeKind::Edge) {
                *uses.entry(edge.node).or_default() += if edge.orientation == Orientation::Forward { 1 } else { -1 };
            }
        }
        assert!(uses.values().all(|&n| n == 0));
        let volume = crate::properties::volume(&model, shell).mass;
        assert!((volume.abs() - 1.0).abs() < 1e-2, "volume {volume}");
        assert!(crate::validity::check(&model, shell).errors().next().is_none());
    }

    #[test]
    fn edges_sewn_to_a_point_stay_where_they_cross_the_face() {
        // A ruled face between two half circles meeting at both ends (one in the XY plane,
        // one in the XZ plane), sewn to the half discs they bound (a closed shell): the ruled
        // face's two rulings at the ends become points in space, but cross its parameters. They stay (as degenerate
        // edges), and the face keeps its area.
        use crate::geometry::curve::Circle;
        use crate::math::{Direction, Frame};
        let mut model = Model::new();
        let half = |model: &mut Model, frame: Frame| make::edge_from_curve_range(model, Curve::Circle(Circle::new(frame, 1.0)), 0.0, std::f64::consts::PI).unwrap();
        let flat = half(&mut model, Frame::WORLD);
        let upright = half(&mut model, Frame::new(Point::ORIGIN, -Direction::Y, Direction::X));
        let sections = [make::wire(&mut model, &[flat]).unwrap(), make::wire(&mut model, &[upright]).unwrap()];
        let options = crate::loft::LoftOptions { ruled: true, ..Default::default() };
        let ruled = crate::loft::loft(&mut model, &sections, &options).unwrap().shape;
        let disc = |model: &mut Model, frame: Frame| {
            let rim = half(model, frame);
            let across = make::edge_between(model, Point::new(-1.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)).unwrap();
            let wire = make::wire(model, &[rim, across]).unwrap();
            make::planar_face(model, wire).unwrap()
        };
        let discs = [disc(&mut model, Frame::WORLD), disc(&mut model, Frame::new(Point::ORIGIN, -Direction::Y, Direction::X))];
        let area = crate::properties::area(&model, ruled).mass + std::f64::consts::PI;
        let sewn = sew(&mut model, &[ruled, discs[0], discs[1]], &SewOptions::default()).unwrap();
        assert_eq!(sewn.kind(), NodeKind::Shell);
        let after = crate::properties::area(&model, sewn).mass;
        assert!((after - area).abs() < 1e-6 * area, "{after} not {area}");
        assert!(crate::validity::check(&model, sewn).errors().next().is_none());
    }

    #[test]
    fn a_long_edge_is_cut_to_meet_two_short_ones() {
        let mut model = Model::new();
        let c = |x: f64, y: f64, z: f64| [x, y, z];
        // One square below two half-width ones: the long top edge meets two short edges.
        let big = square(&mut model, [c(0., 0., 0.), c(2., 0., 0.), c(2., 1., 0.), c(0., 1., 0.)], 0.0);
        let left = square(&mut model, [c(0., 1., 0.), c(1., 1., 0.), c(1., 2., 0.), c(0., 2., 0.)], 0.0);
        let right = square(&mut model, [c(1., 1., 0.), c(2., 1., 0.), c(2., 2., 0.), c(1., 2., 0.)], 0.0);
        let shell = sew(&mut model, &[big, left, right], &SewOptions::default()).unwrap();
        assert_eq!(shell.kind(), NodeKind::Shell);
        assert_eq!(model.unique_subshapes(shell, NodeKind::Face).len(), 3);
        // 4 + 4 + 4 edges, less the 3 shared, the big face's top edge being cut in two.
        assert_eq!(model.unique_subshapes(shell, NodeKind::Edge).len(), 10);
        assert!(crate::validity::check(&model, shell).errors().next().is_none());
    }

    /// Sewing records where its edges come from: the long edge both its pieces, each short
    /// edge the one it is sewn into with a piece of it.
    #[test]
    fn sewing_records_its_edges() {
        let mut model = Model::new();
        let c = |x: f64, y: f64, z: f64| [x, y, z];
        let big = square(&mut model, [c(0., 0., 0.), c(2., 0., 0.), c(2., 1., 0.), c(0., 1., 0.)], 0.0);
        let left = square(&mut model, [c(0., 1., 0.), c(1., 1., 0.), c(1., 2., 0.), c(0., 2., 0.)], 0.0);
        let right = square(&mut model, [c(1., 1., 0.), c(2., 1., 0.), c(2., 2., 0.), c(1., 2., 0.)], 0.0);
        let on_y1 = |model: &Model, f: Shape| {
            model.explore(f, NodeKind::Edge).find(|&e| { let (a, b) = model.edge_vertices(e); model.point(a).y == 1.0 && model.point(b).y == 1.0 }).unwrap()
        };
        let (long, short) = (on_y1(&model, big), on_y1(&model, left));
        crate::boolean::record_begin(&model);
        let shell = sew(&mut model, &[big, left, right], &SewOptions::default()).unwrap();
        let history = crate::boolean::record_end(&model, shell);
        let images = |e: Shape| history.iter().filter(|(from, _)| *from == e.node).map(|(_, to)| *to).collect::<Vec<_>>();
        assert_eq!(images(long).len(), 2);
        assert_eq!(images(short).len(), 1);
        assert!(images(long).contains(&images(short)[0]));
    }
}
