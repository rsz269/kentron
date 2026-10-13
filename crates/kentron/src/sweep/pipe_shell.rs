//! Sweeps through several sections along a spine: each section placed on the spine and
//! seen in the spine's frame there, the sections blended along it (linearly, carried on
//! past the first and the last), the blend carried by the frame.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction2, Point, Point2, Transform, Transformable};
use crate::properties;
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

use super::pipe::{motions, section_placement, PieceMotion, SpinePiece};
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

/// The value at `g` of the law through `(g_s, value_s)` (sorted): linear between
/// neighbours, carried on linearly past the ends.
fn blend(at: &[(f64, Point)], g: f64) -> Point {
    match at {
        [] => Point::ORIGIN,
        [(_, p)] => *p,
        _ => {
            let k = at.windows(2).position(|w| g <= w[1].0).unwrap_or(at.len() - 2);
            let ((g0, p0), (g1, p1)) = (at[k], at[k + 1]);
            let l = if (g1 - g0).abs() > 1e-300 { (g - g0) / (g1 - g0) } else { 0.0 };
            p0 + (p1 - p0) * l
        }
    }
}

/// Sweeps `sections` (wires of as many edges, or vertices) along the spine's `pieces`:
/// a shell, or with `solid` a solid closed by planar faces at the ends.
pub(super) fn pipe_shell(model: &mut Model, sections: &[(Shape, PipeOptions)], pieces: &[SpinePiece], solid: bool) -> Result<Shape> {
    if sections.is_empty() || pieces.is_empty() {
        return Err(Error::InvalidInput("a sweep needs a spine and sections".into()));
    }
    let (reference, _) = motions(pieces, SAMPLES, (0, 0.0))?;
    // Each section: where it stands (as a global spine position) and its motion into the
    // frame at the spine's start.
    let mut placed: Vec<(f64, Shape, Transform)> = Vec::with_capacity(sections.len());
    for &(section, options) in sections {
        let placing = super::placing(model, section, pieces, &options)?;
        let section = section.moved(&placing);
        let (k, s) = section_placement(model, section, pieces)?;
        let (at_anchor, _) = motions(pieces, SAMPLES, (k, s))?;
        placed.push((k as f64 + s, section, at_anchor[0].transforms[0]));
    }
    placed.sort_by(|a, b| a.0.total_cmp(&b.0));
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
    let mut builder = Builder { model, vertices: HashMap::new(), edges: HashMap::new(), paths: HashMap::new() };
    let shell = builder.model.new_empty(NodeKind::Shell);
    let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
    let u_knots = fit_curve(&vec![Point::ORIGIN; pole_count], n)?;
    for k in 0..m {
        let motion = &reference[k];
        for e in 0..edge_count {
            let (v0, v1) = (e, if closed { (e + 1) % vertex_count } else { e + 1 });
            // The surface: each pole's law carried by the frames along the piece.
            let paths: Vec<BSplineCurve> = (0..pole_count).map(|i| carried(motion, k, &law(e, i))).collect::<Result<_>>()?;
            let v = KnotVector::new(paths[0].knots().to_vec(), paths[0].multiplicities().to_vec(), paths[0].degree());
            let u = KnotVector::new(u_knots.knots().to_vec(), u_knots.multiplicities().to_vec(), u_knots.degree());
            let surface = Surface::BSpline(BSplineSurface::try_new(paths.iter().map(|p| p.poles().to_vec()).collect(), u, v)?);
            let bottom = builder.junction_edge(e, k, v0, v1, &junction(k), &|i| law(e, i), pole_count, n, &vertex_law)?;
            let top = builder.junction_edge(e, k + 1, v0, v1, &junction(k + 1), &|i| law(e, i), pole_count, n, &vertex_law)?;
            let right = builder.vertex_path(v1, k, motion, &vertex_law(v1), &junction, &vertex_law)?;
            let left = builder.vertex_path(v0, k, motion, &vertex_law(v0), &junction, &vertex_law)?;
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
                builder.model.set_pcurve_with_range(coedge, pcurve, 0.0, 1.0)?;
                coedges.push(coedge);
            }
            let wire = builder.model.add_wire(&coedges)?;
            builder.model.add(face, wire.into())?;
            builder.model.add(shell, face)?;
        }
    }
    if !solid {
        return Ok(shell);
    }
    // Closed by planar faces at the ends (none where the sweep ends in a point).
    for j in [0, m] {
        let edges: Vec<Shape> = (0..edge_count).filter_map(|e| builder.edges.get(&(e, j)).copied().flatten()).map(Shape::from).collect();
        if edges.len() != edge_count {
            continue;
        }
        let wire = crate::make::wire(builder.model, &edges)?;
        let cap = crate::make::planar_face(builder.model, wire)?;
        builder.model.add(shell, cap)?;
    }
    let volume = properties::closed_shell_volume(builder.model, shell);
    let result = builder.model.new_empty(NodeKind::Solid);
    builder.model.add(result, if volume < 0.0 { shell.reversed() } else { shell })?;
    Ok(result)
}

/// The curve through `poles` as fitted (the knots of a fit through `n + 1` points).
fn fit_curve(poles: &[Point], n: usize) -> Result<BSplineCurve> {
    let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let knots = crate::geometry::fit::interpolate_with_parameters(&vec![Point::ORIGIN; n + 1], &params)?;
    Ok(BSplineCurve::new(poles.to_vec(), knots.knots().to_vec(), knots.multiplicities().to_vec(), knots.degree()))
}

/// The cubic through a point following `law` (in the reference frame) carried by the
/// frames along piece `k`, parametrized by the fraction along it.
fn carried(motion: &PieceMotion, k: usize, law: &[(f64, Point)]) -> Result<BSplineCurve> {
    let points: Vec<Point> = motion.params.iter().zip(&motion.transforms).map(|(&s, t)| blend(law, k as f64 + s).transformed(t)).collect();
    crate::geometry::fit::interpolate_with_parameters(&points, &motion.params)
}

/// The sweep's vertices and edges, made once and shared by the faces meeting along them.
struct Builder<'m> {
    model: &'m mut Model,
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
        let p = blend(&vertex_law(v), j as f64).transformed(at);
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
        let poles: Vec<Point> = (0..pole_count).map(|i| blend(&law(i), j as f64).transformed(at)).collect();
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
        law: &[(f64, Point)],
        junction: &dyn Fn(usize) -> Transform,
        vertex_law: &dyn Fn(usize) -> Vec<(f64, Point)>,
    ) -> Result<EdgeId> {
        if let Some(&id) = self.paths.get(&(v, k)) {
            return Ok(id);
        }
        let (a, b) = (self.vertex(v, k, &junction(k), vertex_law), self.vertex(v, k + 1, &junction(k + 1), vertex_law));
        let curve = carried(motion, k, law)?;
        let edge = self.model.add_edge(a, b, tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), Curve::BSpline(curve), tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), 0.0, 1.0)?;
        self.paths.insert((v, k), edge);
        Ok(edge)
    }
}
