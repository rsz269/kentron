//! Pipes along general spines: the profile moved rigidly along the spine, turning with its
//! tangent and no more (a rotation-minimizing frame, which is OCCT's corrected Frenet
//! trihedron: the Frenet frame with its torsion taken out).

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve, TrimmedCurve};
use crate::geometry::surface::{BSplineSurface, KnotVector};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::Surface;
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Direction2, Frame, Point, Point2, Transform, Transformable, Vector};
use crate::properties;
use crate::tolerance;
use crate::topology::{EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape, VertexId, WireId};

use super::Sweep;

/// One edge of the spine, run from `first` to `last` (which may be the smaller).
pub(super) struct SpinePiece {
    pub curve: Curve,
    pub first: f64,
    pub last: f64,
}

impl SpinePiece {
    fn at(&self, s: f64) -> f64 {
        self.first + (self.last - self.first) * s
    }

    /// The point and unit tangent (in the direction of travel) at fraction `s`.
    pub(super) fn point_tangent(&self, s: f64) -> Result<(Point, Vector)> {
        let (p, d) = self.curve.d1(self.at(s));
        let d = if self.last < self.first { -d } else { d };
        let length = d.length();
        if !(length > 1e-300) {
            return Err(Error::InvalidInput("a pipe's spine with no tangent".into()));
        }
        Ok((p, d * (1.0 / length)))
    }
}

/// The rigid motions along one spine piece: at fractions `params` (from 0 to 1), the
/// transforms taking the start of the whole spine to there.
pub(super) struct PieceMotion {
    pub params: Vec<f64>,
    pub transforms: Vec<Transform>,
}

/// Rotation-minimizing frames along the pieces (double reflection, Wang et al. 2008),
/// carried on from piece to piece and sampled finely: the motions taking the frame at
/// `anchor` (piece, fraction) to those at `samples` points per piece, and to the points
/// halfway between them (to check fits).
pub(super) fn motions(pieces: &[SpinePiece], samples: usize, anchor: (usize, f64)) -> Result<(Vec<PieceMotion>, Vec<Vec<Transform>>)> {
    const FINE: usize = 16;
    let (p0, t0) = pieces[0].point_tangent(0.0)?;
    let r0 = {
        let a = if t0.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
        let r = a - t0 * a.dot(t0);
        r * (1.0 / r.length())
    };
    let frame = |p: Point, t: Vector, r: Vector| -> Frame { Frame::new(p, Direction::from_vector(t), Direction::from_vector(r)) };
    // One double reflection from (x, t, r) to the point x1 with tangent t1: the new r.
    let reflect = |x: Point, t: Vector, r: Vector, x1: Point, t1: Vector| -> Vector {
        let v1 = x1 - x;
        let c1 = v1.dot(v1);
        let (r_l, t_l) = if c1 > 1e-300 { (r - v1 * (2.0 / c1 * v1.dot(r)), t - v1 * (2.0 / c1 * v1.dot(t))) } else { (r, t) };
        let v2 = t1 - t_l;
        let c2 = v2.dot(v2);
        let r1 = if c2 > 1e-300 { r_l - v2 * (2.0 / c2 * v2.dot(r_l)) } else { r_l };
        let r1 = r1 - t1 * r1.dot(t1);
        r1 * (1.0 / r1.length())
    };
    // The raw frames (t, r) carried along; the recorded ones turned by `correction`, which
    // at a corner of the spine keeps the frame as it was at the end of the last piece
    // (OCCT's "modified" transition: the profile does not turn there).
    let (mut x, mut t, mut r) = (p0, t0, r0);
    let mut correction = Correction::identity();
    let mut anchor_frame = None;
    let mut sampled: Vec<(Vec<f64>, Vec<Frame>, Vec<Frame>)> = Vec::with_capacity(pieces.len());
    for (k, piece) in pieces.iter().enumerate() {
        let (p, tangent) = piece.point_tangent(0.0)?;
        if k > 0 {
            if p.distance(x) > 1e3 * tolerance::LINEAR {
                return Err(Error::InvalidInput("a pipe's spine whose edges do not meet".into()));
            }
            // The raw frame turned the least onto the new tangent; the correction then
            // takes it to the corrected frame at the end of the last piece.
            let (ct, cr) = (correction.apply(t), correction.apply(r));
            let axis = t.cross(tangent);
            r = if axis.length() > 1e-12 {
                let (sin, cos) = (axis.length(), t.dot(tangent));
                let a = axis * (1.0 / axis.length());
                r * cos + a.cross(r) * sin + a * (a.dot(r) * (1.0 - cos))
            } else {
                r
            };
            t = tangent;
            let flat = r - t * r.dot(t);
            r = flat * (1.0 / flat.length());
            correction = Correction::between((t, r), (ct, cr));
        }
        x = p;
        let framed = |correction: &Correction, x: Point, t: Vector, r: Vector| frame(x, correction.apply(t), correction.apply(r));
        if anchor.0 == k && anchor.1 <= 0.0 {
            anchor_frame = Some(framed(&correction, x, t, r));
        }
        let n = samples * 2 * FINE;
        let (mut params, mut frames, mut half) = (vec![0.0], vec![framed(&correction, x, t, r)], Vec::with_capacity(samples));
        for i in 1..=n {
            let (s0, s) = ((i - 1) as f64 / n as f64, i as f64 / n as f64);
            if anchor.0 == k && anchor_frame.is_none() && anchor.1 > s0 && anchor.1 <= s {
                let (xa, ta) = piece.point_tangent(anchor.1)?;
                anchor_frame = Some(framed(&correction, xa, ta, reflect(x, t, r, xa, ta)));
            }
            let (x1, t1) = piece.point_tangent(s)?;
            r = reflect(x, t, r, x1, t1);
            (x, t) = (x1, t1);
            if i % (2 * FINE) == 0 {
                params.push(s);
                frames.push(framed(&correction, x, t, r));
            } else if i % (2 * FINE) == FINE {
                half.push(framed(&correction, x, t, r));
            }
        }
        sampled.push((params, frames, half));
    }
    let anchor_frame = anchor_frame.ok_or_else(|| Error::Failed("a pipe's profile placed off its spine".into()))?;
    let to = |f: &Frame| Transform::frame_to_frame(&anchor_frame, f);
    let mut out = Vec::with_capacity(pieces.len());
    let mut halves = Vec::with_capacity(pieces.len());
    for (params, frames, half) in sampled {
        out.push(PieceMotion { params, transforms: frames.iter().map(to).collect() });
        halves.push(half.iter().map(to).collect());
    }
    Ok((out, halves))
}

/// A rotation taking one orthonormal frame `(t, r, t × r)` onto another.
struct Correction {
    from: [Vector; 3],
    to: [Vector; 3],
}

impl Correction {
    fn identity() -> Correction {
        let basis = [Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, 1.0)];
        Correction { from: basis, to: basis }
    }

    fn between((t, r): (Vector, Vector), (u, s): (Vector, Vector)) -> Correction {
        Correction { from: [t, r, t.cross(r)], to: [u, s, u.cross(s)] }
    }

    fn apply(&self, v: Vector) -> Vector {
        self.to[0] * v.dot(self.from[0]) + self.to[1] * v.dot(self.from[1]) + self.to[2] * v.dot(self.from[2])
    }
}

/// The cubic through the moved copies of `p` at the samples, parametrized by the fraction
/// along the piece.
fn path(motion: &PieceMotion, p: Point) -> Result<BSplineCurve> {
    let points: Vec<Point> = motion.transforms.iter().map(|t| p.transformed(t)).collect();
    crate::geometry::fit::interpolate_with_parameters(&points, &motion.params)
}

/// The state of one pipe: the copies of the profile at the spine's junctions (the start,
/// between pieces, the end), what its vertices and edges sweep along each piece.
struct Piper<'m> {
    model: &'m mut Model,
    motions: Vec<PieceMotion>,
    /// The spine closes on itself, the profile coming back to where it started: the last
    /// junction is the first.
    closed: bool,
    vertex_copies: HashMap<(VertexId, usize), VertexId>,
    edge_copies: HashMap<(EdgeId, usize), EdgeId>,
    paths: HashMap<(VertexId, usize), EdgeId>,
    faces: HashMap<(EdgeId, usize), Option<FaceId>>,
    generated: Vec<(NodeId, Vec<Shape>)>,
}

impl Piper<'_> {
    /// The motion to junction `j` (0 the start).
    fn junction(&self, j: usize) -> Transform {
        let j = self.junction_index(j);
        if j < self.motions.len() {
            self.motions[j].transforms[0]
        } else {
            *self.motions[j - 1].transforms.last().expect("samples")
        }
    }

    /// Junction `j`, the last being the first on a closed spine.
    fn junction_index(&self, j: usize) -> usize {
        if self.closed && j == self.motions.len() { 0 } else { j }
    }

    /// True where the profile stands at junction `j` as it is.
    fn unmoved(&self, j: usize) -> bool {
        let t = self.junction(j);
        let rows = t.rows();
        let identity = Transform::IDENTITY.rows();
        rows.iter().flatten().zip(identity.iter().flatten()).all(|(a, b)| (a - b).abs() <= 1e-12)
    }

    fn record(&mut self, input: impl Into<NodeId>, output: Shape) {
        let input = input.into();
        match self.generated.iter_mut().find(|(n, _)| *n == input) {
            Some((_, list)) => list.push(output),
            None => self.generated.push((input, vec![output])),
        }
    }

    fn edge_vertices(&self, e: EdgeId) -> (VertexId, VertexId) {
        let edge = &self.model[e];
        (self.model[edge.start].vertex, self.model[edge.end].vertex)
    }

    fn vertex_copy(&mut self, v: VertexId, j: usize) -> VertexId {
        let j = self.junction_index(j);
        if self.unmoved(j) {
            return v;
        }
        if let Some(&c) = self.vertex_copies.get(&(v, j)) {
            return c;
        }
        let p = self.model[v].point.transformed(&self.junction(j));
        let c = self.model.add_vertex(p, self.model[v].tolerance);
        self.vertex_copies.insert((v, j), c);
        c
    }

    fn edge_copy(&mut self, e: EdgeId, j: usize) -> Result<EdgeId> {
        let j = self.junction_index(j);
        if self.unmoved(j) {
            return Ok(e);
        }
        if let Some(&c) = self.edge_copies.get(&(e, j)) {
            return Ok(c);
        }
        let (a, b) = self.edge_vertices(e);
        let (a, b) = (self.vertex_copy(a, j), self.vertex_copy(b, j));
        let tol = self.model[e].tolerance;
        let c = self.model.add_edge(a, b, tol)?;
        if let Some(slot) = self.model.get(self.model[e].curve) {
            let (curve, first, last) = (slot.curve.transformed(&self.junction(j)), slot.first, slot.last);
            self.model.set_edge_curve(c.into(), curve, tol)?;
            self.model.set_edge_range(c.into(), first, last)?;
        }
        self.edge_copies.insert((e, j), c);
        Ok(c)
    }

    /// The edge swept by vertex `v` along piece `k`, parametrized by the fraction along it.
    fn vertex_path(&mut self, v: VertexId, k: usize) -> Result<EdgeId> {
        if let Some(&e) = self.paths.get(&(v, k)) {
            return Ok(e);
        }
        let (a, b) = (self.vertex_copy(v, k), self.vertex_copy(v, k + 1));
        let curve = path(&self.motions[k], self.model[v].point)?;
        let edge = self.model.add_edge(a, b, tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), Curve::BSpline(curve), tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), 0.0, 1.0)?;
        self.paths.insert((v, k), edge);
        self.record(v, edge.into());
        Ok(edge)
    }

    /// The face swept by edge `e` along piece `k`: on the surface `(u, s)`, u the edge's
    /// parameter and s the fraction along the piece, its normal `∂u × ∂s`; `None` for an
    /// edge without a curve.
    fn edge_face(&mut self, e: EdgeId, k: usize) -> Result<Option<FaceId>> {
        if let Some(&f) = self.faces.get(&(e, k)) {
            return Ok(f);
        }
        let Some(slot) = self.model.get(self.model[e].curve).cloned() else {
            self.faces.insert((e, k), None);
            return Ok(None);
        };
        let (curve, a, b) = (slot.curve, slot.first, slot.last);
        // The profile edge as a B-spline with its own parameter: exact where the conversion
        // keeps it, otherwise a cubic through points along it (closer and closer).
        let exact = crate::geometry::convert::curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(curve.clone(), a, b)));
        let size = curve.point(a).distance(curve.point(0.5 * (a + b))) + curve.point(0.5 * (a + b)).distance(curve.point(b));
        let keeps = |bs: &BSplineCurve| {
            (bs.first_parameter() - a).abs() <= 1e-12 * (1.0 + a.abs())
                && (bs.last_parameter() - b).abs() <= 1e-12 * (1.0 + b.abs())
                && (0..=8).all(|i| {
                    let t = a + (b - a) * i as f64 / 8.0;
                    curve.point(t).distance(bs.point(t)) <= 1e-9 * (1.0 + size)
                })
        };
        // A closed edge, smooth across its ends, sweeps a surface periodic round it.
        let smooth_loop = {
            let ((pa, da), (pb, db)) = (curve.d1(a), curve.d1(b));
            pa.distance(pb) <= tolerance::LINEAR && da.length() > 1e-12 && db.length() > 1e-12 && da.normalized().dot(db.normalized()) > 1.0 - 1e-9
        };
        let profile = match exact.filter(|bs| keeps(bs)) {
            Some(bs) => bs,
            None => {
                let mut n = 32;
                loop {
                    let fit = if smooth_loop {
                        let params: Vec<f64> = (0..n).map(|i| a + (b - a) * i as f64 / n as f64).collect();
                        let points: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
                        crate::geometry::fit::interpolate_periodic_with_parameters(&points, &params, b - a)?
                    } else {
                        let params: Vec<f64> = (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect();
                        let points: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
                        crate::geometry::fit::interpolate_with_parameters(&points, &params)?
                    };
                    let error = (0..n).map(|i| {
                        let t = a + (b - a) * (i as f64 + 0.5) / n as f64;
                        curve.point(t).distance(fit.point(t))
                    }).fold(0.0, f64::max);
                    if error <= 1e-7 * size.max(1.0) || n >= 1024 {
                        break fit;
                    }
                    n *= 2;
                }
            }
        };
        // Each pole moved along the piece (the same knots for all).
        let motion = &self.motions[k];
        let paths: Vec<BSplineCurve> = profile.poles().iter().map(|&q| path(motion, q)).collect::<Result<_>>()?;
        let v = KnotVector::new(paths[0].knots().to_vec(), paths[0].multiplicities().to_vec(), paths[0].degree());
        let u = if profile.is_periodic() {
            KnotVector::periodic(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
        } else {
            KnotVector::new(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
        };
        let poles: Vec<Vec<Point>> = paths.iter().map(|p: &BSplineCurve| p.poles().to_vec()).collect();
        let mut swept = BSplineSurface::try_new(poles, u, v)?;
        if let Some(w) = profile.weights() {
            let rows = paths[0].poles().len();
            swept = swept.with_weights(w.iter().map(|&wi| vec![wi; rows]).collect());
        }
        let surface = Surface::BSpline(swept);
        let (u0, u1) = (a, b);
        let (start, end) = self.edge_vertices(e);
        let (start_path, end_path) = (self.vertex_path(start, k)?, self.vertex_path(end, k)?);
        let (bottom, top) = (self.edge_copy(e, k)?, self.edge_copy(e, k + 1)?);
        let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
        let across = |_: &Model, _: EdgeId, s: f64| -> Result<Curve2> { Ok(line(0.0, s, 1.0, 0.0)) };
        // Counter-clockwise in (u, s): the edge at the start, its end vertex's path, the
        // edge at the end backwards, its start vertex's path backwards.
        let uses = [
            (bottom, Orientation::Forward, across(self.model, bottom, 0.0)?, a, b),
            (end_path, Orientation::Forward, line(u1, 0.0, 0.0, 1.0), 0.0, 1.0),
            (top, Orientation::Reversed, across(self.model, top, 1.0)?, a, b),
            (start_path, Orientation::Reversed, line(u0, 0.0, 0.0, 1.0), 0.0, 1.0),
        ];
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(4);
        for (edge, orientation, pcurve, first, last) in uses {
            let coedge = self.model.add_coedge(edge, orientation)?;
            self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = self.model.add_wire(&coedges)?;
        self.model.add(face, wire.into())?;
        let id = FaceId::from_index(face.node.index);
        self.faces.insert((e, k), Some(id));
        self.record(e, face);
        Ok(Some(id))
    }

    /// The faces swept by edge `e` (used with `orientation`) along every piece.
    fn edge_faces(&mut self, e: EdgeId, orientation: Orientation) -> Result<Vec<Shape>> {
        let mut out = Vec::new();
        for k in 0..self.motions.len() {
            if let Some(face) = self.edge_face(e, k)? {
                out.push(Shape { orientation, ..Shape::from(face) });
            }
        }
        Ok(out)
    }

    /// The face at junction `j`: the profile face, or a moved copy of it (its edges the
    /// copies at the junction, its pcurves the same on the moved surface).
    fn face_copy(&mut self, face: FaceId, j: usize) -> Result<Shape> {
        if self.unmoved(self.junction_index(j)) {
            return Ok(Shape::from(face));
        }
        let tol = self.model[face].tolerance;
        let surface = self.model.face_surface(face).ok_or_else(|| Error::InvalidInput("a profile face without a surface".into()))?.transformed(&self.junction(j));
        let copy = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(copy, surface, tol)?;
        for wire_use in self.model[face].wires.clone() {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            let mut coedges = Vec::new();
            for coedge in self.model[wire].coedges.clone() {
                let (edge, coedge_orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                let moved = self.edge_copy(edge, j)?;
                let c = self.model.add_coedge(moved, coedge_orientation)?;
                if let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) {
                    self.model.set_pcurve_with_range(c, pcurve, first, last)?;
                }
                coedges.push(c);
            }
            let w = self.model.add_wire(&coedges)?;
            self.model.add(copy, Shape { orientation: wire_orientation, ..Shape::from(w) })?;
        }
        Ok(copy)
    }

    /// The solid swept by a face used with `orientation`: closed by the face's copies at
    /// the start and the end (none along a closed spine).
    fn face_solid(&mut self, face: FaceId, orientation: Orientation) -> Result<Shape> {
        let m = self.motions.len();
        let shell = self.model.new_empty(NodeKind::Shell);
        if !self.closed {
            let start = self.face_copy(face, 0)?;
            self.model.add(shell, Shape { orientation: start.orientation.compose(orientation), ..start })?;
        }
        for wire_use in self.model[face].wires.clone() {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            for coedge in self.model[wire].coedges.clone() {
                let (edge, coedge_orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                // (A side runs along the edge where the face it closes runs the other way.)
                let used = coedge_orientation.compose(wire_orientation).compose(orientation);
                for side in self.edge_faces(edge, used.reversed())? {
                    self.model.add(shell, side)?;
                }
            }
        }
        if !self.closed {
            // The copy at the end closes it, facing the other way.
            let end = self.face_copy(face, m)?;
            self.model.add(shell, Shape { orientation: end.orientation.compose(orientation).reversed(), ..end })?;
        }
        let volume = properties::closed_shell_volume(self.model, shell);
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
        self.record(face, solid);
        Ok(solid)
    }
}

/// Sweeps `profile` (unplaced) along the spine's `pieces`, in order.
pub(super) fn pipe_along(model: &mut Model, profile: Shape, pieces: &[SpinePiece]) -> Result<Sweep> {
    if pieces.is_empty() {
        return Err(Error::InvalidInput("a pipe needs a spine".into()));
    }
    // Enough samples for the profile's farthest points to follow their paths closely.
    let size = {
        let b = properties::bounding_box(model, profile);
        let spine: f64 = pieces.iter().map(|p| {
            let (a, b) = (p.curve.point(p.first), p.curve.point(p.last));
            a.distance(p.curve.point(0.5 * (p.first + p.last))) + p.curve.point(0.5 * (p.first + p.last)).distance(b)
        }).sum();
        (b.max - b.min).length() + spine
    };
    let probes: Vec<Point> = {
        let b = properties::bounding_box(model, profile);
        vec![b.min, b.max, Point::new(b.min.x, b.max.y, b.min.z), Point::new(b.max.x, b.min.y, b.max.z)]
    };
    let anchor = section_placement(model, profile, pieces)?;
    let mut samples = 16;
    let motions = loop {
        let (motions, halves) = motions(pieces, samples, anchor)?;
        let mut error = 0.0f64;
        for (motion, half) in motions.iter().zip(&halves) {
            for &p in &probes {
                let fitted = path(motion, p)?;
                for (i, t) in half.iter().enumerate() {
                    let s = 0.5 * (motion.params[i] + motion.params[i + 1]);
                    error = error.max(fitted.point(s).distance(p.transformed(t)));
                }
            }
        }
        if error <= 1e-7 * size.max(1.0) || samples >= 512 {
            break motions;
        }
        samples *= 2;
    };
    // A closed spine (smoothly), the profile coming back to where it started: no ends.
    let closed = {
        let (start, t0) = pieces[0].point_tangent(0.0)?;
        let (end, t1) = pieces[pieces.len() - 1].point_tangent(1.0)?;
        let (first, last) = (motions[0].transforms[0], *motions[motions.len() - 1].transforms.last().expect("samples"));
        let same = first.rows().iter().flatten().zip(last.rows().iter().flatten()).all(|(a, b)| (a - b).abs() <= 1e-7 * (1.0 + a.abs()));
        start.distance(end) <= tolerance::LINEAR && t0.dot(t1) >= 1.0 - 1e-9 && same
    };
    let mut piper = Piper {
        model,
        motions,
        closed,
        vertex_copies: HashMap::new(),
        edge_copies: HashMap::new(),
        paths: HashMap::new(),
        faces: HashMap::new(),
        generated: Vec::new(),
    };
    let shape = match profile.kind() {
        NodeKind::Vertex => {
            let v = VertexId::from_index(profile.node.index);
            let edges: Vec<EdgeId> = (0..piper.motions.len()).map(|k| piper.vertex_path(v, k)).collect::<Result<_>>()?;
            if let [edge] = edges[..] {
                Shape::from(edge)
            } else {
                let coedges = edges.iter().map(|&e| piper.model.add_coedge(e, Orientation::Forward)).collect::<Result<Vec<_>>>()?;
                Shape::from(piper.model.add_wire(&coedges)?)
            }
        }
        NodeKind::Edge => {
            let e = EdgeId::from_index(profile.node.index);
            let faces = piper.edge_faces(e, profile.orientation)?;
            match faces[..] {
                [face] => face,
                _ => {
                    let shell = piper.model.new_empty(NodeKind::Shell);
                    for f in faces {
                        piper.model.add(shell, f)?;
                    }
                    shell
                }
            }
        }
        NodeKind::Wire => {
            let wire = WireId::from_index(profile.node.index);
            let shell = piper.model.new_empty(NodeKind::Shell);
            for coedge in piper.model[wire].coedges.clone() {
                let c = &piper.model[coedge];
                let (edge, orientation) = (c.edge, c.orientation.compose(profile.orientation));
                for f in piper.edge_faces(edge, orientation)? {
                    piper.model.add(shell, f)?;
                }
            }
            shell
        }
        NodeKind::Face => {
            let face = FaceId::from_index(profile.node.index);
            piper.face_solid(face, profile.orientation)?
        }
        NodeKind::Shell | NodeKind::Compound => {
            // Each face (or sub-shape) swept on its own, together.
            let compound = piper.model.new_empty(NodeKind::Compound);
            let faces: Vec<Shape> = piper.model.explore(profile, NodeKind::Face).collect();
            for f in faces {
                let solid = piper.face_solid(FaceId::from_index(f.node.index), f.orientation)?;
                piper.model.add(compound, solid)?;
            }
            compound
        }
        kind => return Err(Error::InvalidInput(format!("a pipe of a {kind:?}"))),
    };
    let generated = std::mem::take(&mut piper.generated);
    Ok(Sweep { shape, generated })
}

/// Where along the spine the profile stands, as OCCT places it (GeomFill_SectionPlacement):
/// a planar profile where the spine crosses its plane nearest its centre (an end of the
/// spine lying in it counting), or else the end nearer the plane; another profile at the
/// spine's point nearest its centre. The profile moves with the spine's frame there.
pub(super) fn section_placement(model: &Model, profile: Shape, pieces: &[SpinePiece]) -> Result<(usize, f64)> {
    let (center, axis) = profile_axis(model, profile);
    let axis = axis.map(|n| (center, n));
    const STEPS: usize = 512;
    let Some((origin, normal)) = axis else {
        // The spine's point nearest the profile.
        let mut best = (f64::INFINITY, (0, 0.0));
        for (k, piece) in pieces.iter().enumerate() {
            for i in 0..=STEPS {
                let s = i as f64 / STEPS as f64;
                let d = piece.curve.point(piece.at(s)).distance(center);
                if d < best.0 {
                    best = (d, (k, s));
                }
            }
        }
        return Ok(best.1);
    };
    let side = |p: Point| (p - origin).dot(normal);
    let last = pieces.len() - 1;
    let mut best: (f64, (usize, f64)) = (f64::INFINITY, (0, 0.0));
    let (first_point, last_point) = (pieces[0].curve.point(pieces[0].at(0.0)), pieces[last].curve.point(pieces[last].at(1.0)));
    if side(first_point).abs() <= 1e-5 {
        best = (first_point.distance(origin), (0, 0.0));
    }
    if side(last_point).abs() <= 1e-5 && last_point.distance(origin) < best.0 {
        best = (last_point.distance(origin), (last, 1.0));
    }
    // Where the spine crosses the plane.
    for (k, piece) in pieces.iter().enumerate() {
        let f = |s: f64| side(piece.curve.point(piece.at(s)));
        let mut previous = f(0.0);
        for i in 1..=STEPS {
            let (s0, s1) = ((i - 1) as f64 / STEPS as f64, i as f64 / STEPS as f64);
            let value = f(s1);
            // A sample on the plane where the spine crosses it (as at a vertex of the spine;
            // not along a stretch of it lying in the plane).
            let crossing = || piece.point_tangent(s1).is_ok_and(|(_, t)| t.dot(normal).abs() > 1e-3 * normal.length());
            if value.abs() <= 1e-12 * (1.0 + (piece.curve.point(piece.at(s1)) - origin).length()) && crossing() {
                let d = piece.curve.point(piece.at(s1)).distance(origin);
                if d < best.0 {
                    best = (d, (k, s1));
                }
            } else if (previous < 0.0) != (value < 0.0) {
                let (mut lo, mut hi) = (s0, s1);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) < 0.0) == (previous < 0.0) { lo = mid } else { hi = mid }
                }
                let s = 0.5 * (lo + hi);
                let d = piece.curve.point(piece.at(s)).distance(origin);
                if d < best.0 {
                    best = (d, (k, s));
                }
            }
            previous = value;
        }
    }
    if best.0.is_finite() {
        return Ok(best.1);
    }
    Ok(if side(first_point).abs() <= side(last_point).abs() { (0, 0.0) } else { (last, 1.0) })
}

/// A profile's centre and, if it is planar, its normal: a conic's own (its centre and
/// axis), or the plane of inertia of points along its edges (their centroid); its
/// vertices' centroid otherwise.
pub(super) fn profile_axis(model: &Model, profile: Shape) -> (Point, Option<Vector>) {
    use crate::geometry::curve::Curve as C;
    let edges: Vec<(Curve, f64, f64)> = model.explore(profile, NodeKind::Edge).filter_map(|e| model.edge_curve(e)).collect();
    match &edges[..] {
        [(C::Circle(c), ..)] => return (c.frame().origin(), Some(Vector::from(c.frame().z()))),
        [(C::Ellipse(c), ..)] => return (c.frame().origin(), Some(Vector::from(c.frame().z()))),
        _ => {}
    }
    let mut points = Vec::new();
    for (c, a, b) in &edges {
        for i in 0..=10 {
            points.push(c.point(a + (b - a) * i as f64 / 10.0));
        }
    }
    if let Some((c, n)) = plane_of_inertia(&points) {
        return (c, Some(n));
    }
    let vertices: Vec<Point> = model.explore(profile, NodeKind::Vertex).map(|v| model.point(v)).collect();
    let n = vertices.len().max(1) as f64;
    (vertices.iter().fold(Point::ORIGIN, |acc, p| acc + (*p - Point::ORIGIN) * (1.0 / n)), None)
}

/// The plane through the points' centroid with the least spread across it (the normal
/// the direction of least inertia), `None` if they lie on a line.
fn plane_of_inertia(points: &[Point]) -> Option<(Point, Vector)> {
    if points.len() < 3 {
        return None;
    }
    let n = points.len() as f64;
    let c = points.iter().fold(Point::ORIGIN, |acc, p| acc + (*p - Point::ORIGIN) * (1.0 / n));
    let mut m = [[0.0f64; 3]; 3];
    for p in points {
        let d = [p.x - c.x, p.y - c.y, p.z - c.z];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] * d[j];
            }
        }
    }
    // Jacobi rotations to diagonalize the symmetric matrix.
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (mut p, mut q, mut largest) = (0, 1, 0.0f64);
        for i in 0..3 {
            for j in i + 1..3 {
                if m[i][j].abs() > largest {
                    (p, q, largest) = (i, j, m[i][j].abs());
                }
            }
        }
        if largest <= 1e-300 {
            break;
        }
        let theta = 0.5 * (2.0 * m[p][q]).atan2(m[q][q] - m[p][p]);
        let (sn, cs) = theta.sin_cos();
        for k in 0..3 {
            let (a, b) = (m[k][p], m[k][q]);
            m[k][p] = cs * a - sn * b;
            m[k][q] = sn * a + cs * b;
        }
        for k in 0..3 {
            let (a, b) = (m[p][k], m[q][k]);
            m[p][k] = cs * a - sn * b;
            m[q][k] = sn * a + cs * b;
        }
        for k in 0..3 {
            let (a, b) = (v[k][p], v[k][q]);
            v[k][p] = cs * a - sn * b;
            v[k][q] = sn * a + cs * b;
        }
    }
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| m[a][a].total_cmp(&m[b][b]));
    let scale = m[order[2]][order[2]].max(1e-300);
    if m[order[1]][order[1]] <= 1e-12 * scale {
        return None;
    }
    let k = order[0];
    Some((c, Vector::new(v[0][k], v[1][k], v[2][k])))
}
