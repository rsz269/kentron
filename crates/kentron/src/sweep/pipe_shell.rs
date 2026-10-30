//! Sweeps through several sections along a spine: each section placed on the spine and
//! seen in the spine's frame there, the sections blended along it by how far along it they
//! stand (carried on past the first and the last), the blend carried by the frame.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction2, Point, Point2, Transform, Transformable};
use crate::properties;
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

use super::pipe::{motions, section_placement, Framing, PieceMotion, SpinePiece};
use super::PipeOptions;

/// Samples per spine piece for the frames and the swept surfaces.
const SAMPLES: usize = 64;

/// A section as poles (in the spine's frame at its start) of each of its edges, all fitted
/// alike; a vertex has every pole at its point.
enum SectionPoles {
    Wire(Vec<Vec<Point>>),
    Vertex(Point),
}

/// The section's edges (in the wire's order, as it runs), its edge curves fitted through
/// `n + 1` points evenly spread in their parameter, in the coordinates `to_reference`
/// takes them to.
fn section_poles(model: &Model, section: Shape, to_reference: &Transform, n: usize) -> Result<SectionPoles> {
    let edges: Vec<Shape> = model.explore(section, NodeKind::Edge).collect();
    if edges.is_empty() {
        let v = model.explore(section, NodeKind::Vertex).next().ok_or_else(|| Error::InvalidInput("an empty sweep section".into()))?;
        return Ok(SectionPoles::Vertex(model.point(v).transformed(to_reference)));
    }
    let mut out = Vec::with_capacity(edges.len());
    for edge in edges {
        let (curve, a, b) = model.edge_curve(edge).ok_or_else(|| Error::InvalidInput("a sweep section edge without a curve".into()))?;
        let (a, b) = if edge.orientation == Orientation::Reversed { (b, a) } else { (a, b) };
        let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let points: Vec<Point> = params.iter().map(|&u| curve.point(a + (b - a) * u).transformed(to_reference)).collect();
        let fit = crate::geometry::fit::interpolate_with_parameters(&points, &params)?;
        out.push(fit.poles().to_vec());
    }
    Ok(SectionPoles::Wire(out))
}

/// The value at `g` of the law through `(g_s, value_s)` (sorted): the polynomial through
/// them all where there are four or fewer (linear between two, as OCCT's; a parabola through
/// three), else through the four nearest `g`; carried on past the ends.
fn blend(at: &[(f64, Point)], g: f64) -> Point {
    match at {
        [] => Point::ORIGIN,
        [(_, p)] => *p,
        _ => {
            // The window of (up to) four sections around `g`.
            let n = at.len();
            let k = at.windows(2).position(|w| g <= w[1].0).unwrap_or(n - 2);
            let lo = if n <= 4 { 0 } else { k.saturating_sub(1).min(n - 4) };
            let window = &at[lo..(lo + 4).min(n)];
            // Lagrange's form (Neville's would do as well for so few).
            let mut out = crate::math::Vector::new(0.0, 0.0, 0.0);
            for (j, &(gj, pj)) in window.iter().enumerate() {
                let mut w = 1.0;
                for (m, &(gm, _)) in window.iter().enumerate() {
                    if m != j && (gj - gm).abs() > 1e-300 {
                        w *= (g - gm) / (gj - gm);
                    }
                }
                out = out + (pj - Point::ORIGIN) * w;
            }
            Point::ORIGIN + out
        }
    }
}

/// Sweeps `sections` (wires of as many edges, or vertices) along the spine's `pieces`:
/// a shell, or with `solid` a solid closed by planar faces at the ends.
pub(super) fn pipe_shell(model: &mut Model, sections: &[(Shape, PipeOptions)], pieces: &[SpinePiece], solid: bool, framing: &Framing) -> Result<super::Sweep> {
    if sections.is_empty() || pieces.is_empty() {
        return Err(Error::InvalidInput("a sweep needs a spine and sections".into()));
    }
    let (reference, _) = motions(pieces, SAMPLES, (0, 0.0), framing)?;
    // How far along the spine a fraction of a piece is, each piece as long as it is (sections
    // are blended by it, as OCCT's by their abscissa on the path; within a piece, by its
    // parameter): from the spine's start, and at each junction.
    let mut junctions = vec![0.0];
    for (k, piece) in pieces.iter().enumerate() {
        let length = crate::geometry::measure::length_between(&piece.curve, piece.first.min(piece.last), piece.first.max(piece.last));
        junctions.push(junctions[k] + length);
    }
    let along = |k: usize, s: f64| junctions[k] + (junctions[k + 1] - junctions[k]) * s;
    // (The samples of the frames along each piece, that far along.)
    let sampled: Vec<Vec<f64>> = reference.iter().enumerate().map(|(k, motion)| motion.params.iter().map(|&s| along(k, s)).collect()).collect();
    // Each section: where it stands (how far along the spine) and its motion into the frame
    // at the spine's start.
    let mut placed: Vec<(f64, Shape, Transform)> = Vec::with_capacity(sections.len());
    for &(section, options) in sections {
        // (Where it stands found before it is placed, as OCCT finds it.)
        let (k, s) = section_placement(model, section, pieces)?;
        let placing = super::placing(model, section, pieces, &options)?;
        let section = section.moved(&placing);
        let (at_anchor, _) = motions(pieces, SAMPLES, (k, s), framing)?;
        placed.push((along(k, s), section, at_anchor[0].transforms[0]));
    }
    placed.sort_by(|a, b| a.0.total_cmp(&b.0));
    let first_section = placed[0].1.node;
    // The sections seen from the spine's start (all in one plane there), then made alike as a
    // loft's: turned the same way, started at matching vertices, split to match.
    // (Each section's edges, by node, with the input's edge they are or are part of.)
    let mut sources: Vec<HashMap<NodeId, NodeId>> = Vec::with_capacity(placed.len());
    for p in placed.iter_mut() {
        let original: Vec<NodeId> = model.explore(p.1, NodeKind::Edge).map(|e| e.node).collect();
        p.1 = crate::modify::transform(model, p.1, &p.2, crate::modify::TransformOptions::default())?;
        p.2 = Transform::IDENTITY;
        sources.push(model.explore(p.1, NodeKind::Edge).map(|e| e.node).zip(original).collect());
    }
    if placed.iter().filter(|(_, s, _)| s.kind() == NodeKind::Wire).count() > 1 {
        let shapes: Vec<Shape> = placed.iter().map(|(_, s, _)| *s).collect();
        let matched = crate::loft::compatible_sections(model, &shapes, crate::loft::Matching::ByRays)?;
        for ((p, (m, pieces)), sources) in placed.iter_mut().zip(matched).zip(sources.iter_mut()) {
            p.1 = m;
            *sources = pieces.into_iter().filter_map(|(edge, copy)| sources.get(&copy).map(|&source| (edge, source))).collect();
        }
    }
    let edge_sources: Vec<Vec<Option<NodeId>>> =
        placed.iter().zip(&sources).map(|((_, s, _), sources)| model.explore(*s, NodeKind::Edge).map(|e| sources.get(&e.node).copied()).collect()).collect();
    let mut generated: Vec<(NodeId, Vec<Shape>)> = Vec::new();
    // The sections' edges fitted alike (enough points for the closest fit).
    let mut n = 32;
    let poles = loop {
        let poles: Vec<SectionPoles> = placed.iter().map(|(_, s, t)| section_poles(model, *s, t, n)).collect::<Result<_>>()?;
        let fine: Vec<SectionPoles> = placed.iter().map(|(_, s, t)| section_poles(model, *s, t, 2 * n)).collect::<Result<_>>()?;
        let mut error = 0.0f64;
        for (a, b) in poles.iter().zip(&fine) {
            if let (SectionPoles::Wire(a), SectionPoles::Wire(b)) = (a, b) {
                for (pa, pb) in a.iter().zip(b) {
                    let (ca, cb) = (fit_curve(pa, n)?, fit_curve(pb, 2 * n)?);
                    for i in 0..=4 * n {
                        let u = i as f64 / (4 * n) as f64;
                        error = error.max(ca.point(u).distance(cb.point(u)));
                    }
                }
            }
        }
        if error <= 1e-7 || n >= 512 {
            break poles;
        }
        n *= 2;
    };
    let edge_count = poles.iter().filter_map(|p| match p { SectionPoles::Wire(e) => Some(e.len()), _ => None }).max().unwrap_or(0);
    if edge_count == 0 {
        return Err(Error::InvalidInput("a sweep needs a section with edges".into()));
    }
    if poles.iter().any(|p| matches!(p, SectionPoles::Wire(e) if e.len() != edge_count)) {
        return Err(Error::InvalidInput("sweep sections with different numbers of edges are not supported yet".into()));
    }
    let closed = poles.iter().all(|p| match p {
        SectionPoles::Wire(e) => e[0][0].distance(*e[edge_count - 1].last().expect("poles")) <= 1e3 * tolerance::LINEAR,
        SectionPoles::Vertex(_) => true,
    });
    let pole_count = n + 1;
    // The law of pole i of edge e: (global position, reference point) per section.
    let law = |e: usize, i: usize| -> Vec<(f64, Point)> {
        placed
            .iter()
            .zip(&poles)
            .map(|((g, _, _), p)| (*g, match p { SectionPoles::Wire(edges) => edges[e][i], SectionPoles::Vertex(v) => *v }))
            .collect()
    };
    let junction = |j: usize| -> Transform {
        if j < reference.len() { reference[j].transforms[0] } else { *reference[j - 1].transforms.last().expect("samples") }
    };
    let vertex_count = if closed { edge_count } else { edge_count + 1 };
    let vertex_law = |v: usize| -> Vec<(f64, Point)> { if v < edge_count { law(v, 0) } else { law(edge_count - 1, pole_count - 1) } };
    let m = pieces.len();
    let mut builder = Builder { model, junctions: &junctions, vertices: HashMap::new(), edges: HashMap::new(), paths: HashMap::new() };
    let shell = builder.model.new_empty(NodeKind::Shell);
    let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
    let u_knots = fit_curve(&vec![Point::ORIGIN; pole_count], n)?;
    for k in 0..m {
        let motion = &reference[k];
        for e in 0..edge_count {
            let (v0, v1) = (e, if closed { (e + 1) % vertex_count } else { e + 1 });
            // The surface: each pole's law carried by the frames along the piece.
            // (Fitted alike, as few poles along as hold the samples: few where the frames and
            // the sections' blend vary smoothly, as along a straight spine.)
            let samples: Vec<Vec<Point>> = (0..pole_count).map(|i| carried_points(motion, &sampled[k], &law(e, i))).collect();
            let paths = crate::geometry::fit::approximate_alike(&samples, &motion.params, 3, 1e-8)?;
            let v = KnotVector::new(paths[0].knots().to_vec(), paths[0].multiplicities().to_vec(), paths[0].degree());
            let u = KnotVector::new(u_knots.knots().to_vec(), u_knots.multiplicities().to_vec(), u_knots.degree());
            let swept = BSplineSurface::try_new(paths.iter().map(|p| p.poles().to_vec()).collect(), u, v)?;
            // (Planar where it lies in a plane, its edges' pcurves found from it.)
            let planar = super::pipe::planar_surface(&swept);
            let surface = planar.clone().map(Surface::Plane).unwrap_or(Surface::BSpline(swept));
            let bottom = builder.junction_edge(e, k, v0, v1, &junction(k), &|i| law(e, i), pole_count, n, &vertex_law)?;
            let top = builder.junction_edge(e, k + 1, v0, v1, &junction(k + 1), &|i| law(e, i), pole_count, n, &vertex_law)?;
            let right = builder.vertex_path(v1, k, motion, &sampled[k], &vertex_law(v1), &junction, &vertex_law)?;
            let left = builder.vertex_path(v0, k, motion, &sampled[k], &vertex_law(v0), &junction, &vertex_law)?;
            // Where the section is a point, the face closes at it with no edge, unless both
            // its sides are one edge (a seam): then by a degenerate edge of its own there.
            let point = |builder: &mut Builder, j: usize| -> Result<Option<EdgeId>> {
                if v0 != v1 {
                    return Ok(None);
                }
                let v = builder.vertex(v0, j, &junction(j), &vertex_law);
                Ok(Some(builder.model.add_edge(v, v, tolerance::LINEAR)?))
            };
            let bottom = match bottom { Some(b) => Some(b), None => point(&mut builder, k)? };
            let top = match top { Some(t) => Some(t), None => point(&mut builder, k + 1)? };
            let uses: Vec<(EdgeId, Orientation, Curve2)> = [
                bottom.map(|b| (b, Orientation::Forward, line(0.0, 0.0, 1.0, 0.0))),
                Some((right, Orientation::Forward, line(1.0, 0.0, 0.0, 1.0))),
                top.map(|t| (t, Orientation::Reversed, line(0.0, 1.0, 1.0, 0.0))),
                Some((left, Orientation::Reversed, line(0.0, 0.0, 0.0, 1.0))),
            ]
            .into_iter()
            .flatten()
            .collect();
            let face = builder.model.new_empty(NodeKind::Face);
            builder.model.set_face_surface(face, surface, tolerance::LINEAR)?;
            let mut coedges = Vec::with_capacity(4);
            for (edge, orientation, pcurve) in uses {
                let coedge = builder.model.add_coedge(edge, orientation)?;
                if planar.is_none() {
                    builder.model.set_pcurve_with_range(coedge, pcurve, 0.0, 1.0)?;
                }
                coedges.push(coedge);
            }
            let wire = builder.model.add_wire(&coedges)?;
            builder.model.add(face, wire.into())?;
            builder.model.add(shell, face)?;
            for source in edge_sources.iter().filter_map(|s| s.get(e).copied().flatten()) {
                match generated.iter_mut().find(|(n, _)| *n == source) {
                    Some((_, faces)) => faces.push(face),
                    None => generated.push((source, vec![face])),
                }
            }
        }
    }
    if !solid {
        return Ok(super::Sweep { shape: shell, generated, ends: Vec::new(), copies: Vec::new() });
    }
    let mut caps = Vec::new();
    // Closed by planar faces at the ends (none where the sweep ends in a point).
    for j in [0, m] {
        let edges: Vec<Shape> = (0..edge_count).filter_map(|e| builder.edges.get(&(e, j)).copied().flatten()).map(Shape::from).collect();
        if edges.len() != edge_count {
            continue;
        }
        let wire = crate::make::wire(builder.model, &edges)?;
        let cap = crate::make::planar_face(builder.model, wire)?;
        // (Using the section's edges the other way from the sides: those run along the
        // first section's edges and against the last's.)
        let first = edges[0].node;
        let along = builder.model.explore(cap, NodeKind::Edge).find(|e| e.node == first).is_some_and(|e| e.orientation == Orientation::Forward);
        let cap = if along == (j == 0) { cap.reversed() } else { cap };
        builder.model.add(shell, cap)?;
        caps.push((j, cap));
    }
    let volume = properties::closed_shell_volume(builder.model, shell);
    let result = builder.model.new_empty(NodeKind::Solid);
    builder.model.add(result, if volume < 0.0 { shell.reversed() } else { shell })?;
    // (The first section's caps: its first and last faces.)
    let ends = match caps[..] {
        [(0, first), (_, last)] => vec![(first_section, first, last)],
        _ => Vec::new(),
    };
    Ok(super::Sweep { shape: result, generated, ends, copies: Vec::new() })
}

/// The curve through `poles` as fitted (the knots of a fit through `n + 1` points).
fn fit_curve(poles: &[Point], n: usize) -> Result<BSplineCurve> {
    let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let knots = crate::geometry::fit::interpolate_with_parameters(&vec![Point::ORIGIN; n + 1], &params)?;
    Ok(BSplineCurve::new(poles.to_vec(), knots.knots().to_vec(), knots.multiplicities().to_vec(), knots.degree()))
}

/// A point following `law` (in the reference frame) carried by the frames along a piece, at
/// its samples (`along`, how far along the spine each is).
fn carried_points(motion: &PieceMotion, along: &[f64], law: &[(f64, Point)]) -> Vec<Point> {
    along.iter().zip(&motion.transforms).map(|(&g, t)| blend(law, g).transformed(t)).collect()
}

/// The cubic through a point following `law` carried by the frames along a piece,
/// parametrized by the fraction along it.
fn carried(motion: &PieceMotion, along: &[f64], law: &[(f64, Point)]) -> Result<BSplineCurve> {
    crate::geometry::fit::interpolate_with_parameters(&carried_points(motion, along, law), &motion.params)
}

/// The sweep's vertices and edges, made once and shared by the faces meeting along them.
struct Builder<'m> {
    model: &'m mut Model,
    /// How far along the spine each junction is.
    junctions: &'m [f64],
    vertices: HashMap<(usize, usize), VertexId>,
    /// The edges of the sections at the junctions; `None` where the section is a point.
    edges: HashMap<(usize, usize), Option<EdgeId>>,
    paths: HashMap<(usize, usize), EdgeId>,
}

impl Builder<'_> {
    fn vertex(&mut self, v: usize, j: usize, at: &Transform, vertex_law: &dyn Fn(usize) -> Vec<(f64, Point)>) -> VertexId {
        if let Some(&id) = self.vertices.get(&(v, j)) {
            return id;
        }
        let p = blend(&vertex_law(v), self.junctions[j]).transformed(at);
        // (Where the section is a point, its vertices are one.)
        let same = self.vertices.iter().find(|(&(_, jj), &id)| jj == j && self.model[id].point.distance(p) <= tolerance::LINEAR).map(|(_, &id)| id);
        let id = same.unwrap_or_else(|| self.model.add_vertex(p, tolerance::LINEAR));
        self.vertices.insert((v, j), id);
        id
    }

    /// The section's edge `e` at junction `j` (from vertex `v0` to `v1`); a degenerate
    /// edge of its own where the section is a point there (each face its own).
    #[allow(clippy::too_many_arguments)]
    fn junction_edge(
        &mut self,
        e: usize,
        j: usize,
        v0: usize,
        v1: usize,
        at: &Transform,
        law: &dyn Fn(usize) -> Vec<(f64, Point)>,
        pole_count: usize,
        n: usize,
        vertex_law: &dyn Fn(usize) -> Vec<(f64, Point)>,
    ) -> Result<Option<EdgeId>> {
        if let Some(found) = self.edges.get(&(e, j)) {
            return Ok(*found);
        }
        let poles: Vec<Point> = (0..pole_count).map(|i| blend(&law(i), self.junctions[j]).transformed(at)).collect();
        let (a, b) = (self.vertex(v0, j, at, vertex_law), self.vertex(v1, j, at, vertex_law));
        if poles.iter().all(|p| p.distance(poles[0]) <= tolerance::LINEAR) {
            // A point: the faces close at it, with no edge.
            self.edges.insert((e, j), None);
            return Ok(None);
        }
        let edge = self.model.add_edge(a, b, tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), Curve::BSpline(fit_curve(&poles, n)?), tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), 0.0, 1.0)?;
        self.edges.insert((e, j), Some(edge));
        Ok(Some(edge))
    }

    fn vertex_path(
        &mut self,
        v: usize,
        k: usize,
        motion: &PieceMotion,
        along: &[f64],
        law: &[(f64, Point)],
        junction: &dyn Fn(usize) -> Transform,
        vertex_law: &dyn Fn(usize) -> Vec<(f64, Point)>,
    ) -> Result<EdgeId> {
        if let Some(&id) = self.paths.get(&(v, k)) {
            return Ok(id);
        }
        let (a, b) = (self.vertex(v, k, &junction(k), vertex_law), self.vertex(v, k + 1, &junction(k + 1), vertex_law));
        let curve = carried(motion, along, law)?;
        let edge = self.model.add_edge(a, b, tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), Curve::BSpline(curve), tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), 0.0, 1.0)?;
        self.paths.insert((v, k), edge);
        Ok(edge)
    }
}
