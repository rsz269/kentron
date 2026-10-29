//! Sewing the kept regions into the result: shared vertices and edges, faces with pcurves,
//! shells (kept apart where solids only touch along an edge) and solids with their voids.

use crate::hashing::FastMap as HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::project;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Transform2, Transformable2, Vector, Vector2};
use crate::properties;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

use super::geom::{self, Uv, TOL};
use super::graph::Graph;
use super::split::{ArcOnFace, Piece, Region};

/// A region kept in the result, turned inside out if `flip`.
pub(super) struct Kept {
    pub face: usize,
    pub region: Region,
    pub flip: bool,
    /// `face` is a far face (in `Graph::far`), kept whole as it is (`region` unused).
    pub far: bool,
}

struct OutEdge {
    id: EdgeId,
    curve: Option<Curve>,
    first: f64,
    last: f64,
    /// The largest distance found between the edge's geometry and its vertices or pcurves.
    deviation: f64,
    /// An operand's edge, kept whole.
    reused: bool,
}

/// One use of a result edge by a face (in `kept`), from which the direction from the edge
/// into the face and the face's outward normal are found where needed.
struct EdgeUse {
    face: usize,
    /// Whether the face runs along the edge's curve.
    along: bool,
    /// +1 if the face's normal points out of the result.
    sign: f64,
}

/// The result, its edges made from section curves (where the operands' faces meet), and
/// whether its shells are closed (each edge used by an even number of faces).
/// `invalid_operands` is asked (only then) where the result has a void inside no solid: a
/// bug with valid operands, which is reported; with invalid ones (as OCCT carries on) the
/// void is dropped.
pub(super) fn build(
    model: &mut Model,
    graph: &Graph,
    pieces: &[Piece],
    arcs_of: &[Vec<ArcOnFace>],
    kept: &[Kept],
    invalid_operands: &dyn Fn(&Model) -> bool,
) -> Result<Option<(Shape, Vec<Shape>, bool)>> {
    let mut out = Output::default();
    let faces = out.faces(model, graph, pieces, arcs_of, kept)?;
    if out.far_conflict {
        return Ok(None);
    }
    out.set_tolerances(model, graph, pieces)?;
    let vectors = |edge: EdgeId, u: &EdgeUse| -> (Vector, Vector) {
        let (curve, first, last) = model.edge_curve(Shape::from(edge)).expect("uses are of edges with curves");
        let (p, tangent) = curve.d1(0.5 * (first + last));
        let tangent = if u.along { tangent } else { -tangent };
        into_face(kept_output(graph, &kept[u.face]), p, tangent, u.sign, tangent.length() * (last - first))
    };
    let shells = group_shells(&out.uses, faces.len(), vectors);
    let closed = out.uses.values().all(|u| u.len() % 2 == 0);
    let result = assemble(model, graph, kept, &faces, shells, invalid_operands)?;
    Ok(Some((result, out.sections(graph, pieces), closed)))
}

/// Several results from one graph, sharing their vertices, edges and faces: every region
/// in `all` is made into a face once (as it is), and each selection of them (by index in
/// `all`, turned inside out where marked) is assembled into solids.
pub(super) fn build_partition(
    model: &mut Model,
    graph: &Graph,
    pieces: &[Piece],
    arcs_of: &[Vec<ArcOnFace>],
    all: &[Kept],
    selections: &[Vec<(usize, bool)>],
    invalid_operands: &dyn Fn(&Model) -> bool,
) -> Result<Vec<Shape>> {
    let mut out = Output::default();
    let faces = out.faces(model, graph, pieces, arcs_of, all)?;
    out.set_tolerances(model, graph, pieces)?;
    let mut results = Vec::with_capacity(selections.len());
    for selection in selections {
        if selection.is_empty() {
            results.push(model.new_empty(NodeKind::Compound));
            continue;
        }
        let local: HashMap<usize, (usize, bool)> = selection.iter().enumerate().map(|(li, &(i, flip))| (i, (li, flip))).collect();
        let kept: Vec<Kept> = selection.iter().map(|&(i, flip)| Kept { face: all[i].face, region: all[i].region.clone(), flip, far: false }).collect();
        let cell_faces: Vec<Shape> = selection.iter().map(|&(i, flip)| if flip { faces[i].reversed() } else { faces[i] }).collect();
        let mut uses: HashMap<EdgeId, Vec<EdgeUse>> = HashMap::default();
        for (&edge, list) in &out.uses {
            for u in list {
                if let Some(&(li, flip)) = local.get(&u.face) {
                    uses.entry(edge).or_default().push(EdgeUse { face: li, along: u.along, sign: if flip { -u.sign } else { u.sign } });
                }
            }
        }
        let vectors = |edge: EdgeId, u: &EdgeUse| -> (Vector, Vector) {
            let (curve, first, last) = model.edge_curve(Shape::from(edge)).expect("uses are of edges with curves");
            let (p, tangent) = curve.d1(0.5 * (first + last));
            let tangent = if u.along { tangent } else { -tangent };
            into_face(kept_output(graph, &kept[u.face]), p, tangent, u.sign, tangent.length() * (last - first))
        };
        let shells = group_shells(&uses, cell_faces.len(), vectors);
        results.push(assemble(model, graph, &kept, &cell_faces, shells, invalid_operands)?);
    }
    Ok(results)
}

/// The face's own surface of a kept region.
fn kept_output<'a>(graph: &'a Graph, k: &Kept) -> &'a crate::geometry::surface::Surface {
    if k.far { &graph.far[k.face].data.output } else { &graph.faces[k.face].output }
}

/// The result's vertices and edges, made once and shared (keyed by the graph's vertices and
/// representative pieces), and how faces use the edges.
#[derive(Default)]
pub(super) struct Output {
    vertices: HashMap<usize, VertexId>,
    edges: HashMap<usize, OutEdge>,
    /// Per result edge: its uses by the result's faces (by index in `kept`).
    uses: HashMap<EdgeId, Vec<EdgeUse>>,
    /// An edge of a far face (kept whole, as it is) that the result could not keep.
    pub far_conflict: bool,
}

impl Output {
    /// The result's edges made from section curves.
    pub fn sections(&self, graph: &Graph, pieces: &[Piece]) -> Vec<Shape> {
        let mut sections: Vec<(usize, EdgeId)> = self
            .edges
            .iter()
            .filter(|(&rep, _)| graph.edges[pieces[rep].edge].operand.is_none())
            .map(|(&rep, e)| (rep, e.id))
            .collect();
        sections.sort_unstable();
        sections.into_iter().map(|(_, id)| Shape::from(id)).collect()
    }

    /// The edge standing for piece `index` (made if needed), oriented like the piece.
    pub fn edge(&mut self, model: &mut Model, graph: &Graph, pieces: &[Piece], index: usize) -> Result<Shape> {
        let piece = &pieces[index];
        let rep = &pieces[piece.rep];
        let same = piece.rep == index || same_direction(graph, piece, rep);
        let out = edge_for(model, graph, rep, piece.rep, &mut self.vertices, &mut self.edges)?;
        let shape = Shape::from(out.id);
        Ok(if same { shape } else { shape.reversed() })
    }

    /// Each made face's and edge's operand face or edge, for [`super::record_end`]: a kept
    /// region's face from its operand's face, a result edge from each operand edge with a
    /// piece it stands for.
    pub fn record(&self, model: crate::topology::ModelUid, graph: &Graph, pieces: &[Piece], kept: &[Kept], faces: &[Shape]) {
        for (k, face) in kept.iter().zip(faces) {
            let source = if k.far { Some(graph.far[k.face].face) } else { graph.faces[k.face].regions_key };
            if let Some(source) = source {
                super::record(model, source.into(), face.node);
            }
        }
        for piece in pieces {
            if let (Some(out), Some(source)) = (self.edges.get(&piece.rep), graph.edges[piece.edge].source) {
                super::record(model, source.into(), out.id.into());
            }
        }
    }

    /// The vertex standing for the graph's vertex `v`.
    pub fn vertex(&mut self, model: &mut Model, graph: &Graph, v: usize) -> VertexId {
        *self.vertices.entry(v).or_insert_with(|| result_vertex(model, graph, v))
    }

    /// Pairs of result faces (indices in `kept`) sharing an edge.
    pub fn neighbours(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for around in self.uses.values() {
            for w in around.windows(2) {
                pairs.push((w[0].face, w[1].face));
            }
        }
        pairs
    }

    /// The faces of the kept regions, with their outward orientation.
    pub fn faces(&mut self, model: &mut Model, graph: &Graph, pieces: &[Piece], arcs_of: &[Vec<ArcOnFace>], kept: &[Kept]) -> Result<Vec<Shape>> {
        let Output { vertices, edges, uses, far_conflict } = self;
        // Far faces keep their edges: one built anew in place of such an edge would leave
        // them bounded by edges the result does not have.
        let far_edges: crate::hashing::FastSet<EdgeId> = kept.iter().filter(|k| k.far).flat_map(|k| graph.far[k.face].data.edges.iter().map(|e| e.edge)).collect();
        let mut faces: Vec<Shape> = Vec::new();
        for (ki, k) in kept.iter().enumerate() {
            if k.far {
                let far = &graph.far[k.face];
                let sign = far.outward * if k.flip { -1.0 } else { 1.0 };
                for e in far.data.edges.iter().filter(|e| e.curve.is_some()) {
                    uses.entry(e.edge).or_default().push(EdgeUse { face: ki, along: e.along, sign });
                }
                let face = Shape::from(far.face);
                faces.push(if sign > 0.0 { face } else { face.reversed() });
                continue;
            }
            let face_data = &graph.faces[k.face];
            let surface = face_data.output.clone();
            let sign = face_data.outward * if k.flip { -1.0 } else { 1.0 };
            // A face kept whole (bounded by its own edges, each kept whole) is the operand's
            // face itself: only how it uses its edges is recorded.
            if let Some(origin) = face_data.origin {
                // The pieces it is bounded by (its arcs), each with whether the face runs
                // along it.
                let arcs: Vec<(usize, bool)> = k.region.loops.iter().flatten().copied().collect();
                let own = |ai: usize| {
                    let arc = &arcs_of[k.face][ai];
                    let piece = &pieces[arc.piece];
                    piece.rep == arc.piece && face_data.boundary.iter().any(|b| b.edge == piece.edge)
                };
                let along_pieces: Option<Vec<(usize, bool)>> = (arcs.len() == face_data.boundary.len() && arcs.iter().all(|&(ai, _)| own(ai)))
                    .then(|| arcs.iter().map(|&(ai, forward)| (arcs_of[k.face][ai].piece, forward)).collect());
                if let Some(along_pieces) = along_pieces {
                    let mut whole = true;
                    for &(pi, _) in &along_pieces {
                        let out = edge_for(model, graph, &pieces[pi], pi, vertices, edges)?;
                        *far_conflict |= !out.reused && graph.edges[pieces[pi].edge].origin.is_some_and(|id| far_edges.contains(&id));
                        whole &= out.reused;
                    }
                    if whole {
                        // The face's own pcurves may stray from its edges: the edges'
                        // tolerances are to cover that too.
                        struct Strays(HashMap<EdgeId, f64>);
                        let strays = model.derived(origin.into(), model.face_fingerprint(origin), || {
                            let mut strays: HashMap<EdgeId, f64> = HashMap::default();
                            for wire in model.children(Shape::from(origin)).collect::<Vec<_>>() {
                                for coedge in model.children(wire).collect::<Vec<_>>() {
                                    let Ok(coedge) = crate::topology::CoEdgeId::try_from(coedge.node) else { continue };
                                    let edge = model[coedge].edge;
                                    let (Some((pcurve, first, last)), Some((curve, ..))) = (model.effective_pcurve(coedge), model.edge_curve(Shape::from(edge))) else { continue };
                                    let gap = (0..=8)
                                        .map(|i| {
                                            let t = first + (last - first) * i as f64 / 8.0;
                                            let uv = pcurve.point(t);
                                            surface.point(uv.x, uv.y).distance(curve.point(t))
                                        })
                                        .fold(0.0, f64::max);
                                    let entry = strays.entry(edge).or_insert(0.0);
                                    *entry = entry.max(gap);
                                }
                            }
                            Strays(strays)
                        });
                        for &(pi, along) in &along_pieces {
                            let out = edges.get_mut(&pi).expect("edge");
                            out.deviation = out.deviation.max(strays.0.get(&out.id).copied().unwrap_or(0.0));
                            if out.curve.is_some() {
                                uses.entry(out.id).or_default().push(EdgeUse { face: ki, along, sign });
                            }
                        }
                        let face = Shape::from(origin);
                        faces.push(if sign > 0.0 { face } else { face.reversed() });
                        continue;
                    }
                }
            }
            let face = model.new_empty(NodeKind::Face);
            model.set_face_surface(face, surface.clone(), TOL)?;
            for cycle in &k.region.loops {
                let mut coedges = Vec::with_capacity(cycle.len());
                for &(ai, forward) in cycle {
                    let arc = &arcs_of[k.face][ai];
                    let piece = &pieces[arc.piece];
                    let rep = &pieces[piece.rep];
                    let same_direction = piece.rep == arc.piece || same_direction(graph, piece, rep);
                    let out = edge_for(model, graph, rep, piece.rep, vertices, edges)?;
                    let conflict = !out.reused && graph.edges[rep.edge].origin.is_some_and(|id| far_edges.contains(&id));
                    *far_conflict |= conflict;
                    let orientation = if forward == same_direction { Orientation::Forward } else { Orientation::Reversed };
                    let coedge = model.add_coedge(out.id, orientation)?;
                    let pcurve = match (&arc.pcurve, &out.curve) {
                        (Some(p), _) if piece.rep == arc.piece => p.clone(),
                        (Some(p), None) => p.clone(),
                        (_, Some(curve)) => {
                            // The arc's middle in (u, v): on its own pcurve if it has one (its
                            // polyline may be just its two ends).
                            let middle = match &arc.pcurve {
                                Some(p) => {
                                    let q = p.point(0.5 * (piece.first + piece.last));
                                    [q.x, q.y]
                                }
                                None => arc.uv[arc.uv.len() / 2],
                            };
                            // (On the surface computed with where it is the output one cut
                            // down near the face (the same parameters): carried on past the face,
                            // the output one may fold back over it, a foot found on the fold.)
                            let near = match (&face_data.surface, &surface) {
                                (crate::geometry::surface::Surface::BSpline(_), crate::geometry::surface::Surface::BSpline(_)) => &face_data.surface,
                                _ => &surface,
                            };
                            // (Across a seam where the surface closes without a period, its two
                            // uses taken to their sides by the bounds' width, as by a period.)
                            let closed = |k: usize| face_data.closed_bounds[k].map(|(lo, hi)| hi - lo);
                            let periods = [face_data.u_period.or_else(|| closed(0)), face_data.v_period.or_else(|| closed(1))];
                            // (From the arc's end at the curve's start, where the start is on a
                            // seam: its foot alone may be found on either side.)
                            let start = [arc.uv[0], arc.uv[arc.uv.len() - 1]]
                                .into_iter()
                                .find(|uv| geom::surface_point(near, *uv).distance(curve.point(out.first)) <= 1e3 * TOL);
                            aligned_pcurve(curve, near, out.first, out.last, middle, periods, face_data.drift, start)?
                        }
                        (None, None) => return Err(Error::InvalidInput("a degenerate edge without a pcurve".into())),
                    };
                    if let Some(curve) = &out.curve {
                        for i in 0..=8 {
                            let t = out.first + (out.last - out.first) * i as f64 / 8.0;
                            let uv = pcurve.point(t);
                            out.deviation = out.deviation.max(surface.point(uv.x, uv.y).distance(curve.point(t)));
                        }
                        uses.entry(out.id).or_default().push(EdgeUse { face: ki, along: orientation == Orientation::Forward, sign });
                    }
                    model.set_pcurve_with_range(coedge, pcurve, out.first, out.last)?;
                    coedges.push(coedge);
                }
                let wire = model.add_wire(&coedges)?;
                model.add(face, wire.into())?;
            }
            faces.push(if sign > 0.0 { face } else { face.reversed() });
        }
        self.record(model.uid(), graph, pieces, kept, &faces);
        Ok(faces)
    }

    /// Tolerances covering the gaps found between the edges' curves, pcurves and vertices.
    pub fn set_tolerances(&self, model: &mut Model, graph: &Graph, pieces: &[Piece]) -> Result<()> {
        let (vertices, edges) = (&self.vertices, &self.edges);
        let mut vertex_tolerance: HashMap<usize, f64> = HashMap::default();
        for (rep, out) in edges {
            let piece = &pieces[*rep];
            if let Some(curve) = &out.curve {
                let gaps = [
                    curve.point(out.first).distance(graph.vertices.points[piece.start]),
                    curve.point(out.last).distance(graph.vertices.points[piece.end]),
                ];
                let tolerance = (out.deviation.max(gaps[0]).max(gaps[1]) * 1.5).max(TOL);
                // An operand's edge is only ever made more tolerant (and left alone, unedited,
                // when tolerant enough).
                if !out.reused || tolerance > model.tolerance(Shape::from(out.id)) {
                    model.set_edge_tolerance(out.id, tolerance)?;
                }
                for (v, gap) in [(piece.start, gaps[0]), (piece.end, gaps[1])] {
                    let t = vertex_tolerance.entry(v).or_insert(TOL);
                    *t = t.max(gap * 1.5).max(tolerance);
                }
            }
        }
        for (v, tolerance) in vertex_tolerance {
            let id = vertices[&v];
            let reused = graph.vertex_origin.get(&v) == Some(&id);
            if reused {
                // An operand's vertex keeps its point and is only made more tolerant.
                if tolerance > model.tolerance(Shape::from(id)) {
                    let point = model.point(Shape::from(id));
                    model.set_vertex_point(id.into(), point, tolerance)?;
                }
            } else if tolerance > TOL {
                model.set_vertex_point(id.into(), graph.vertices.points[v], tolerance)?;
            }
        }
        Ok(())
    }
}

/// True if two coincident pieces run the same way.
fn same_direction(graph: &Graph, piece: &Piece, rep: &Piece) -> bool {
    if piece.start != piece.end {
        return piece.start == rep.start;
    }
    let at = |p: &Piece, s: f64| graph.edges[p.edge].point(p.first + (p.last - p.first) * s, &graph.vertices);
    at(piece, 0.25).distance(at(rep, 0.25)) < at(piece, 0.25).distance(at(rep, 0.75))
}

fn edge_for<'a>(
    model: &mut Model,
    graph: &Graph,
    rep: &Piece,
    rep_index: usize,
    vertices: &mut HashMap<usize, VertexId>,
    edges: &'a mut HashMap<usize, OutEdge>,
) -> Result<&'a mut OutEdge> {
    if !edges.contains_key(&rep_index) {
        let mut vertex = |v: usize, model: &mut Model| *vertices.entry(v).or_insert_with(|| result_vertex(model, graph, v));
        let (start, end) = (vertex(rep.start, model), vertex(rep.end, model));
        let edge = &graph.edges[rep.edge];
        let curve = edge.curve.clone();
        // An operand's edge kept whole, between the result's vertices, is reused.
        let whole = rep.first == edge.first && rep.last == edge.last && rep.start == edge.start && rep.end == edge.end;
        let original = edge.origin.filter(|&id| {
            let (a, b) = model.edge_vertices(Shape::from(id));
            whole && (VertexId::try_from(a.node).ok(), VertexId::try_from(b.node).ok()) == (Some(start), Some(end))
        });
        let (id, reused) = match original {
            Some(id) => (id, true),
            None => {
                let id = model.add_edge(start, end, TOL)?;
                if let Some(c) = &curve {
                    model.set_edge_curve(id.into(), c.clone(), TOL)?;
                    model.set_edge_range(id.into(), rep.first, rep.last)?;
                }
                (id, false)
            }
        };
        edges.insert(rep_index, OutEdge { id, curve, first: rep.first, last: rep.last, deviation: 0.0, reused });
    }
    Ok(edges.get_mut(&rep_index).expect("edge"))
}

/// The result's vertex for the graph's vertex `v`: the operand's own where it is one.
fn result_vertex(model: &mut Model, graph: &Graph, v: usize) -> VertexId {
    match graph.vertex_origin.get(&v) {
        Some(&id) => id,
        None => model.add_vertex(graph.vertices.points[v], TOL),
    }
}

/// The pcurve of `curve` over `[first, last]` on `surface`, shifted by whole periods of
/// `u` to have its middle at the arc's, `target`; straying from the surface by no more than
/// the face's edges may (`drift`) where nothing nearer holds. Fitted from `start`, the
/// parameters of the curve's first point, where given.
fn aligned_pcurve(curve: &Curve, surface: &crate::geometry::surface::Surface, first: f64, last: f64, target: Uv, periods: [Option<f64>; 2], drift: f64, start: Option<Uv>) -> Result<Curve2> {
    // An analytic pcurve may hold on part of the curve only (a meridian up to the poles,
    // beyond which it would leave the surface's parameter domain): it must match over the
    // edge's range, within the domain.
    let (_, _, v_min, v_max) = surface.bounds();
    let analytic = project::analytic_pcurve(curve, surface).filter(|p| {
        (0..=8).all(|i| {
            let t = first + (last - first) * i as f64 / 8.0;
            let uv = p.point(t);
            uv.y >= v_min - 1e-9 && uv.y <= v_max + 1e-9 && surface.point(uv.x, uv.y).distance(curve.point(t)) <= TOL
        })
    });
    let pcurve = match analytic {
        Some(p) => p,
        None => {
            // A traced section (interpolated through traced points) may stray from the surface
            // by more than the fine tolerances: a looser fit then, the edge's tolerance
            // covering the gap.
            let size = {
                let (a, b) = (curve.point(first), curve.point(last));
                a.distance(curve.point(0.5 * (first + last))) + curve.point(0.5 * (first + last)).distance(b)
            };
            let loose = vec![1e-3 * size.max(1e-3), 1e-2 * size.max(1e-3)];
            // (Then a polyline a little looser, where the curve turns too sharply in (u, v) for
            // a smooth fit, as where a face narrows to nothing, and strays a little.)
            // (From the curve's first point's foot as the operation finds it, where that is on
            // it: the projection's own may be on another part of the surface passing nearly as
            // near, as the next turn of one winding round.)
            let start = start.map(|uv| (uv[0], uv[1])).or_else(|| {
                let p = curve.point(first);
                let uv = geom::uv_of(surface, p);
                (surface.point(uv[0], uv[1]).distance(p) <= 1e3 * TOL).then_some((uv[0], uv[1]))
            });
            project::pcurve_on_range_from(curve, surface, first, last, &[1e-8, 1e-6, 1e-4], start)
                .or_else(|| project::polyline_pcurve_on_range_from(curve, surface, first, last, 1e-7, start))
                .or_else(|| project::polyline_pcurve_on_range_from(curve, surface, first, last, 1e-5, start))
                .or_else(|| project::pcurve_on_range_from(curve, surface, first, last, &loose, start))
                // (A section in the tolerance of the face's edges, as where faces meet within
                // it past the surface's domain: as far off as those edges may be.)
                .or_else(|| (drift > TOL).then(|| project::polyline_pcurve_on_range_from(curve, surface, first, last, drift, start)).flatten())
                // (Or twice as far, as an edge of the other shape taken as lying on a face
                // whose edges are that loose (a filling's), but leaving it near a corner: the
                // edge's tolerance then covers how far it strays.)
                .or_else(|| (drift > TOL).then(|| project::polyline_pcurve_on_range_from(curve, surface, first, last, 2.0 * drift, start)).flatten())
                .map(Curve2::BSpline)
                .ok_or_else(|| Error::InvalidInput("no pcurve for a section edge".into()))?
        }
    };
    // The arc's middle against the pcurve's, in each periodic direction.
    let mid = pcurve.point(0.5 * (first + last));
    let shift = |period: Option<f64>, target: f64, at: f64| period.map_or(0.0, |p| ((target - at) / p).round() * p);
    let (du, dv) = (shift(periods[0], target[0], mid.x), shift(periods[1], target[1], mid.y));
    if du == 0.0 && dv == 0.0 {
        return Ok(pcurve);
    }
    Ok(pcurve.transformed(&Transform2::translation(Vector2::new(du, dv))))
}

/// Where a face leaves an edge it bounds: the direction into it, across the edge (bent the
/// way the face curves, a little way in, so that faces leaving tangent to each other, as
/// round solids touching along a line do, still come in their order around the edge), and
/// its normal out of the result.
fn into_face(surface: &crate::geometry::surface::Surface, p: Point, tangent: Vector, sign: f64, length: f64) -> (Vector, Vector) {
    let normal = geom::normal_at_point(surface, p);
    // (Where the edge has no direction, or the normal runs along it, as at a face's degenerate
    // point, the way into the face is unknown: some way square to the normal, the shells'
    // checks judging after.)
    let across = normal.cross(tangent);
    if !(across.length() > 1e-12 * tangent.length()) || !(across.length() > 0.0) {
        let basis = [normal, tangent].into_iter().find(|v| v.length() > 0.0).unwrap_or(Vector::new(0.0, 0.0, 1.0)).normalized();
        let other = if basis.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
        return (basis.cross(other).normalized(), normal * sign);
    }
    let into = across.normalized();
    if matches!(surface, crate::geometry::surface::Surface::Plane(_)) {
        return (into, normal * sign);
    }
    // The point `h` along `into` on the tangent plane, taken onto the surface by its (u, v).
    let uv = geom::uv_of(surface, p);
    let (at, su, sv) = surface.d1(uv[0], uv[1]);
    let h = 1e-4 * length.max(TOL);
    let (uu, uv_, vv) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let det = uu * vv - uv_ * uv_;
    let (a, b) = (h * into.dot(su), h * into.dot(sv));
    let along = tangent.normalized();
    let chord = surface.point(uv[0] + (a * vv - b * uv_) / det, uv[1] + (b * uu - a * uv_) / det) - at;
    let chord = chord - along * chord.dot(along);
    let bent = if det > 0.0 && chord.length() > 0.0 { chord.normalized() } else { into };
    (bent, normal * sign)
}

/// Groups faces into shells by shared edges. Where more than two faces meet at an edge
/// (solids touching along it), each is joined to its neighbour around the edge on the side
/// of its material, so that each solid keeps a shell of its own.
fn group_shells(uses: &HashMap<EdgeId, Vec<EdgeUse>>, count: usize, vectors: impl Fn(EdgeId, &EdgeUse) -> (Vector, Vector)) -> Vec<Vec<usize>> {
    let mut group: Vec<usize> = (0..count).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    fn join(group: &mut [usize], i: usize, j: usize) {
        let (ri, rj) = (root(group, i), root(group, j));
        group[ri] = rj;
    }
    let mut keys: Vec<EdgeId> = uses.keys().copied().collect();
    keys.sort_unstable_by_key(|k| k.index());
    for key in keys {
        let around = &uses[&key];
        if around.len() <= 2 {
            if let [a, b] = &around[..] {
                join(&mut group, a.face, b.face);
            }
            continue;
        }
        // The edge's direction: across a use's normal and into-direction (the first whose
        // are not degenerate; where none is, the faces joined in turn).
        let found: Vec<(Vector, Vector)> = around.iter().map(|u| vectors(key, u)).collect();
        let Some(first) = found.iter().position(|(into, outward)| into.cross(*outward).length() > 0.0) else {
            for w in around.windows(2) {
                join(&mut group, w[0].face, w[1].face);
            }
            continue;
        };
        let d = found[first].0.cross(found[first].1).normalized();
        let (a, b) = (found[first].0, d.cross(found[first].0));
        let mut sorted: Vec<(f64, usize, bool)> = around
            .iter()
            .zip(&found)
            .map(|(u, &(into, outward))| (into.dot(b).atan2(into.dot(a)), u.face, d.cross(into).dot(outward) < 0.0))
            .collect();
        sorted.sort_by(|x, y| x.0.total_cmp(&y.0));
        for k in 0..sorted.len() {
            if sorted[k].2 {
                join(&mut group, sorted[k].1, sorted[(k + 1) % sorted.len()].1);
            }
        }
    }
    let mut members: HashMap<usize, Vec<usize>> = HashMap::default();
    for i in 0..count {
        members.entry(root(&mut group, i)).or_default().push(i);
    }
    let mut shells: Vec<Vec<usize>> = members.into_values().collect();
    shells.sort();
    shells
}

/// Shells into solids: those enclosing a positive volume are outer shells, the others
/// voids, put in the smallest solid around them. The result is a compound of solids.
fn assemble(
    model: &mut Model,
    graph: &Graph,
    kept: &[Kept],
    faces: &[Shape],
    shells: Vec<Vec<usize>>,
    invalid_operands: &dyn Fn(&Model) -> bool,
) -> Result<Shape> {
    let mut built: Vec<(Shape, f64, Vec<usize>)> = Vec::new();
    // A lone shell of several faces bounds the result from outside: its volume (slow to
    // integrate on free-form faces) is not needed to tell.
    let lone = shells.len() == 1 && shells[0].len() > 1;
    // The faces' terms, found for all faces of shells of several at once.
    let in_many: Vec<usize> = if lone { Vec::new() } else { shells.iter().filter(|m| m.len() > 1).flatten().copied().collect() };
    let model_now: &Model = model;
    let terms: HashMap<usize, f64> = in_many.iter().copied().zip(crate::parallel::map(&in_many, |&i| properties::face_volume_term(model_now, faces[i]))).collect();
    for members in shells {
        let shell = model.new_empty(NodeKind::Shell);
        for &i in &members {
            model.add(shell, faces[i])?;
        }
        // From the faces' terms, which the model keeps for faces met before; the terms
        // about the origin add up to the volume only for a closed shell, and a lone face
        // (a sliver of a face, or a closed one such as a whole sphere) is integrated as
        // before.
        let volume = if lone {
            1.0
        } else if members.len() == 1 {
            properties::volume(model, shell).mass
        } else {
            members.iter().map(|&i| terms[&i]).sum()
        };
        built.push((shell, volume, members));
    }
    // Slivers enclosing no volume (stray fragments of faces along nearly coincident
    // surfaces) bound nothing.
    let largest = built.iter().map(|(_, v, _)| v.abs()).fold(0.0, f64::max);
    built.retain(|(_, v, members)| members.len() > 1 || v.abs() > 1e-9 * largest);
    let (outer, voids): (Vec<_>, Vec<_>) = built.into_iter().partition(|(_, v, _)| *v > 0.0);
    let mut solids: Vec<(Vec<Shape>, f64, Vec<usize>)> = outer.into_iter().map(|(s, v, m)| (vec![s], v, m)).collect();
    for (shell, _, members) in voids {
        let k = &kept[members[0]];
        let face = region_face(graph, k);
        let interior = if k.far {
            let own: Vec<&[Uv]> = face.polygons.iter().map(|p| p.as_slice()).collect();
            geom::interior_point_among(face.polygons, &geom::Obstacles::new(&own))
        } else {
            k.region.interior
        };
        let probe = face.surface.point(interior[0], interior[1]);
        let around = solids
            .iter()
            .enumerate()
            .filter(|(_, (_, _, m))| encloses(graph, kept, m, probe))
            .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .map(|(i, _)| i);
        match around {
            Some(i) => solids[i].0.push(shell),
            None if invalid_operands(model) => {}
            None => return Err(Error::InvalidInput("the result has a void outside any solid".into())),
        }
    }
    let compound = model.new_empty(NodeKind::Compound);
    for (shells, _, _) in solids {
        let solid = model.new_empty(NodeKind::Solid);
        for shell in shells {
            model.add(solid, shell)?;
        }
        model.add(compound, solid)?;
    }
    Ok(compound)
}

/// True if the regions `members` enclose `p` (ray casting).
/// A kept region as rays and probes see it: its face's surface and box, and its polygons
/// (a far face's own boundary) in the face's (u, v), with the face's periods and ranges.
struct RegionFace<'a> {
    surface: &'a crate::geometry::surface::Surface,
    bbox: (Point, Point),
    polygons: &'a [Vec<Uv>],
    periods: [Option<f64>; 2],
    ranges: ((f64, f64), (f64, f64)),
}

fn region_face<'a>(graph: &'a Graph, k: &'a Kept) -> RegionFace<'a> {
    if k.far {
        let data = &graph.far[k.face].data;
        let sample = &data.sample;
        RegionFace { surface: &data.surface, bbox: sample.bbox, polygons: &sample.polylines, periods: data.periods, ranges: (sample.u_range, sample.v_range) }
    } else {
        let face = &graph.faces[k.face];
        RegionFace { surface: &face.surface, bbox: face.bbox, polygons: &k.region.polygons, periods: [face.u_period, face.v_period], ranges: (face.u_range, face.v_range) }
    }
}

fn encloses(graph: &Graph, kept: &[Kept], members: &[usize], p: Point) -> bool {
    const DIRECTIONS: [(f64, f64, f64); 3] = [(0.5773, 0.3141, 0.7536), (-0.2718, 0.8414, 0.4142), (0.6931, -0.5403, 0.1732)];
    'rays: for (x, y, z) in DIRECTIONS {
        let d = Vector::new(x, y, z).normalized();
        let mut crossings = 0;
        for &i in members {
            let face = region_face(graph, &kept[i]);
            let Some(hits) = geom::ray_hits(face.surface, p, d, &face.bbox) else { continue 'rays };
            for t in hits {
                let uv = geom::uv_of(face.surface, p + d * t);
                let copies = super::graph::representatives(uv, face.periods, face.ranges.0, face.ranges.1);
                if copies.into_iter().any(|uv| geom::inside_polylines(face.polygons, uv)) {
                    crossings += 1;
                }
            }
        }
        return crossings % 2 == 1;
    }
    false
}
