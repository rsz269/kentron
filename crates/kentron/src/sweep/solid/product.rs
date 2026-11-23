//! The swept solid under a steady motion, built as the product of the boundary of `W` and the
//! motion, capped by the rest of the body at the start and `W` at the end (see the module
//! above).

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::BSplineCurve;
use crate::geometry::surface::{BSplineSurface, KnotVector};
use crate::geometry::curve::{Circle, Curve, Line, TrimmedCurve};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{Cylinder, Plane, RevolutionSurface, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Direction2, Frame, Point, Point2, Transform, Transformable, Vector};
use crate::sweep::pipe::{path, PieceMotion};
use crate::tolerance;
use crate::topology::{EdgeId, FaceId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{contour, Steady};

/// Which way the motion carries a face (or a piece of one) relative to its outward normal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Outwards somewhere (`V·n > 0`): in `W`.
    Front,
    /// Inwards, or along the face everywhere.
    Back,
}

/// The swept solid of `solid` (unplaced; changed: its faces are divided where the motion
/// grazes them) under `steady`.
pub(super) fn sweep(model: &mut Model, solid: Shape, steady: &Steady) -> Result<Shape> {
    let bounds = crate::properties::bounding_box(model, solid);
    let size = (bounds.max - bounds.min).length();
    // The characteristic curves within a billionth of the body's size.
    let tolerance = 1e-9 * size;
    // Faces closed round a seam are split, so that each curve is followed within one period.
    for _ in 0..2 {
        crate::heal::split_closed_faces(model, solid)?;
    }
    for face in model.unique_subshapes(solid, NodeKind::Face) {
        let cuts = contour::cuts(model, face, steady, tolerance)?;
        if !cuts.is_empty() {
            crate::divide::divide_face(model, face, &cuts)?;
        }
    }
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let speed = speed_bound(model, solid, steady);
    let sides: Vec<Side> = faces.iter().map(|&f| side(model, f, steady, speed)).collect::<Result<_>>()?;
    let folds = split_folds(model, &faces, &sides, steady)?;
    // (How far the body reaches from a turn's axis: the farthest corner of its box.)
    let reach = (0..8)
        .map(|i| {
            let p = Point::new(if i & 1 == 0 { bounds.min.x } else { bounds.max.x }, if i & 2 == 0 { bounds.min.y } else { bounds.max.y }, if i & 4 == 0 { bounds.min.z } else { bounds.max.z });
            match steady {
                Steady::Turn { axis, .. } | Steady::Screw { axis, .. } => {
                    let r = p - axis.origin;
                    (r - Vector::from(axis.direction) * r.dot(Vector::from(axis.direction))).length()
                }
                Steady::Slide { .. } => 0.0,
            }
        })
        .fold(0.0, f64::max);
    let mut builder = Builder::new(model, *steady, size, reach)?;
    let product = builder.product(&faces, &sides)?;
    if folds == 0 {
        return Ok(product);
    }
    resolved(model, product)
}

/// The product where it runs into itself (see the module above): its faces split where they
/// cross, the pieces of space it winds round at least once, joined. Itself where it does not.
fn resolved(model: &mut Model, product: Shape) -> Result<Shape> {
    let faces = model.unique_subshapes(product, NodeKind::Face);
    let cells = crate::boolean::enclosed_solids(model, &faces)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    if cells.is_empty() {
        return Err(Error::Failed("a sweep running into itself, not split into pieces".into()));
    }
    if cells.len() == 1 {
        return Ok(product);
    }
    let classifier = crate::boolean::PointClassifier::new(model, product)?;
    let mut kept = Vec::new();
    for &cell in &cells {
        let p = crate::boolean::point_inside(model, cell)?.ok_or_else(|| Error::Failed("a piece of a sweep with no point inside".into()))?;
        let winding = classifier.winding(p).ok_or_else(|| Error::Failed("how often a sweep winds round a point of it".into()))?;
        if winding >= 1 {
            kept.push(cell);
        }
    }
    super::super::joined(model, &kept)
}

/// The edges between a face in `W` and one outside it (each with whether the face in `W` runs
/// along the edge's curve), but degenerate ones: the boundary of `W`.
fn boundary_edges(model: &Model, faces: &[Shape], sides: &[Side]) -> Result<Vec<(EdgeId, bool)>> {
    // Each edge's uses: the face (by index) and whether it runs along the edge's curve.
    let mut uses: HashMap<EdgeId, Vec<(usize, bool)>> = HashMap::new();
    for (k, &face) in faces.iter().enumerate() {
        let loops = crate::heal::face_loops(model, face).ok_or_else(|| Error::Failed("a face without pcurves".into()))?;
        let reversed_face = face.orientation == Orientation::Reversed;
        for u in loops.iter().flatten() {
            uses.entry(u.edge).or_default().push((k, u.reversed == reversed_face));
        }
    }
    let mut edges: Vec<EdgeId> = uses.keys().copied().collect();
    edges.sort_by_key(|e| e.index());
    let mut out = Vec::new();
    for e in edges {
        let list = &uses[&e];
        let fronts: Vec<&(usize, bool)> = list.iter().filter(|(k, _)| sides[*k] == Side::Front).collect();
        if list.len() != 2 || fronts.len() != 1 || list[0].0 == list[1].0 || model.is_degenerate(Shape::from(e)) {
            continue;
        }
        out.push((e, fronts[0].1));
    }
    Ok(out)
}

/// Splits the edges bounding `W` where the surface they sweep folds: where the motion runs
/// along the edge (its tangent crossed with the velocity, the swept surface's normal,
/// vanishes and turns over), as a characteristic curve does where it turns past the direction
/// of a slide. Each piece then sweeps a surface without a fold, meeting the next along the
/// fold's path. Returns how many folds were found.
fn split_folds(model: &mut Model, faces: &[Shape], sides: &[Side], steady: &Steady) -> Result<usize> {
    let mut count = 0;
    for (e, _) in boundary_edges(model, faces, sides)? {
        let Some((curve, a, b)) = model.edge_curve(Shape::from(e)) else { continue };
        let normal = |s: f64| {
            let (p, d) = curve.d1(s);
            d.cross(steady.velocity(p))
        };
        let samples = match &curve {
            Curve::BSpline(bs) => 64.max(4 * bs.poles().len()),
            _ => 64,
        };
        let inside = |at: f64| at - a > 1e-9 * (b - a) && b - at > 1e-9 * (b - a);
        let mut roots = Vec::new();
        let mut previous = (a, normal(a));
        for i in 1..=samples {
            let s = a + (b - a) * i as f64 / samples as f64;
            let n = normal(s);
            if n.dot(previous.1) < 0.0 {
                let reference = previous.1;
                let f = |x: f64| normal(previous.0 + (s - previous.0) * x).dot(reference);
                let x = super::contour::root(&f, f(0.0), f(1.0));
                let at = previous.0 + (s - previous.0) * x;
                if inside(at) {
                    roots.push(at);
                    count += 1;
                }
            }
            previous = (s, n);
        }
        roots.sort_by(f64::total_cmp);
        let mut piece = e;
        for at in roots {
            let tolerance = model[piece].tolerance;
            let v = model.add_vertex(curve.point(at), tolerance);
            piece = model.split_edge(piece, v, at)?.1;
        }
    }
    Ok(count)
}

/// The largest speed over the body's box.
fn speed_bound(model: &Model, solid: Shape, steady: &Steady) -> f64 {
    let b = crate::properties::bounding_box(model, solid);
    let mut speed = 0.0f64;
    for i in 0..8 {
        let p = Point::new(if i & 1 == 0 { b.min.x } else { b.max.x }, if i & 2 == 0 { b.min.y } else { b.max.y }, if i & 4 == 0 { b.min.z } else { b.max.z });
        speed = speed.max(steady.velocity(p).length());
    }
    speed
}

/// Which way the motion carries `face` (as its shell uses it): from `V·n` at points inside
/// it, the largest either way.
fn side(model: &Model, face: Shape, steady: &Steady, speed: f64) -> Result<Side> {
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
    let polygons = contour::loop_polygons(model, face).ok_or_else(|| Error::Failed("a face without pcurves".into()))?;
    let gauge = contour::Gauge { surface, steady };
    let sign = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    let (u0, u1, v0, v1) = model.face_uv_bounds(face);
    const N: usize = 24;
    let (mut most, mut least) = (0.0f64, 0.0f64);
    for i in 0..N {
        for j in 0..N {
            let q = Point2::new(u0 + (u1 - u0) * (i as f64 + 0.5) / N as f64, v0 + (v1 - v0) * (j as f64 + 0.5) / N as f64);
            if !contour::within(&polygons, q) {
                continue;
            }
            let g = sign * gauge.g(q.x, q.y);
            most = most.max(g);
            least = least.min(g);
        }
    }
    // (V·n of a face the motion runs along is zero but for rounding. Points near a piece's
    // boundary may fall just across the curve bounding it, as the loops are sampled: the
    // piece goes by the larger.)
    let still = 1e-9 * speed;
    Ok(if most > still && most >= -least { Side::Front } else { Side::Back })
}

/// Builds the product's faces: the paths of vertices, the copies at the end, the faces edges
/// sweep.
///
/// A screw is built in stages of at most a quarter turn (each the same screw, carrying the
/// copies the stage before ended at), so that no face winds round its axis: each is a strip
/// short enough that a point's place on it is found from any start on it.
struct Builder<'m> {
    model: &'m mut Model,
    /// The motion of one stage.
    steady: Steady,
    /// Where the stages start and the last ends: the motion's places.
    junctions: Vec<Transform>,
    /// A stage of the screw as the pipes' motions sample it (its points' paths quintics
    /// through them).
    screw: Option<PieceMotion>,
    size: f64,
    vertex_copies: HashMap<(VertexId, usize), VertexId>,
    edge_copies: HashMap<(EdgeId, usize), EdgeId>,
    paths: HashMap<(VertexId, usize), EdgeId>,
}

impl<'m> Builder<'m> {
    fn new(model: &'m mut Model, steady: Steady, size: f64, reach: f64) -> Result<Builder<'m>> {
        let (stage, stages) = match steady {
            Steady::Screw { axis, angle, rise } => {
                let n = (angle / std::f64::consts::FRAC_PI_2).ceil().max(1.0);
                (Steady::Screw { axis, angle: angle / n, rise: rise / n }, n as usize)
            }
            other => (other, 1),
        };
        let screw = match stage {
            Steady::Screw { axis, angle, rise } => Some(screw_motion(axis, angle, rise, size, reach)?),
            _ => None,
        };
        let junctions = (0..=stages).map(|k| steady.at(k as f64 / stages as f64)).collect();
        Ok(Builder { model, steady: stage, junctions, screw, size, vertex_copies: HashMap::new(), edge_copies: HashMap::new(), paths: HashMap::new() })
    }

    fn product(&mut self, faces: &[Shape], sides: &[Side]) -> Result<Shape> {
        let shell = self.model.new_empty(NodeKind::Shell);
        // The faces outside W where the motion starts.
        for (k, &face) in faces.iter().enumerate() {
            if sides[k] == Side::Back {
                self.model.add(shell, face)?;
            }
        }
        // What the boundary of W sweeps, and W at the end.
        let stages = self.junctions.len() - 1;
        for (e, along) in boundary_edges(self.model, faces, sides)? {
            for k in 0..stages {
                let face = self.side_face(e, along, k)?;
                self.model.add(shell, face)?;
            }
        }
        for (k, &face) in faces.iter().enumerate() {
            if sides[k] == Side::Front {
                let copy = self.face_copy(face)?;
                self.model.add(shell, copy)?;
            }
        }
        let volume = crate::properties::closed_shell_volume(self.model, shell);
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
        Ok(solid)
    }

    /// A vertex's copy at junction `k` (itself at the start).
    fn vertex_copy(&mut self, v: VertexId, k: usize) -> VertexId {
        if k == 0 {
            return v;
        }
        if let Some(&c) = self.vertex_copies.get(&(v, k)) {
            return c;
        }
        let (p, tol) = (self.model[v].point, self.model[v].tolerance);
        let c = self.model.add_vertex(p.transformed(&self.junctions[k]), tol);
        self.vertex_copies.insert((v, k), c);
        c
    }

    /// An edge's copy at junction `k` (itself at the start).
    fn edge_copy(&mut self, e: EdgeId, k: usize) -> Result<EdgeId> {
        if k == 0 {
            return Ok(e);
        }
        if let Some(&c) = self.edge_copies.get(&(e, k)) {
            return Ok(c);
        }
        let (a, b) = self.model.edge_end_vertices(e);
        let (a, b) = (self.vertex_copy(a, k), self.vertex_copy(b, k));
        let tol = self.model[e].tolerance;
        let c = self.model.add_edge(a, b, tol)?;
        if let Some((curve, first, last)) = self.model.edge_curve(Shape::from(e)) {
            self.model.set_edge_curve(c.into(), curve.transformed(&self.junctions[k]), tol)?;
            self.model.set_edge_range(c.into(), first, last)?;
        }
        self.edge_copies.insert((e, k), c);
        Ok(c)
    }

    /// The copy at the end of a face, used as the body uses it (its edges the copies of its
    /// edges, its pcurves the same on the moved surface).
    fn face_copy(&mut self, face: Shape) -> Result<Shape> {
        let id = FaceId::from_index(face.node.index);
        let tol = self.model[id].tolerance;
        let end = *self.junctions.last().expect("junctions");
        let stages = self.junctions.len() - 1;
        let surface = self.model.face_surface(id).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.transformed(&end);
        let copy = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(copy, surface, tol)?;
        for wire_use in self.model[id].wires.clone() {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            let mut coedges = Vec::new();
            for coedge in self.model[wire].coedges.clone() {
                let (edge, orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                let moved = self.edge_copy(edge, stages)?;
                let c = self.model.add_coedge(moved, orientation)?;
                if let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) {
                    self.model.set_pcurve_with_range(c, pcurve, first, last)?;
                }
                coedges.push(c);
            }
            let w = self.model.add_wire(&coedges)?;
            self.model.add(copy, Shape { orientation: wire_orientation, ..Shape::from(w) })?;
        }
        Ok(Shape { orientation: face.orientation, ..copy })
    }

    /// The edge a vertex's path makes over stage `k`, from its copy at junction `k` to that at
    /// the next.
    fn vertex_path(&mut self, v: VertexId, k: usize) -> Result<EdgeId> {
        if let Some(&e) = self.paths.get(&(v, k)) {
            return Ok(e);
        }
        let (start, end) = (self.vertex_copy(v, k), self.vertex_copy(v, k + 1));
        let p = self.model[start].point;
        let tol = self.model[v].tolerance.max(tolerance::LINEAR);
        let edge = self.model.add_edge(start, end, tol)?;
        let (curve, first, last) = match self.steady {
            Steady::Slide { by } => (Curve::Line(Line::new(p, Direction::from_vector(by))), 0.0, by.length()),
            Steady::Turn { axis, angle } => (Curve::Circle(turn_circle(axis, p)?), 0.0, angle),
            Steady::Screw { .. } => (Curve::BSpline(path(self.screw.as_ref().expect("a screw's motion"), p)?), 0.0, 1.0),
        };
        self.model.set_edge_curve(edge.into(), curve, tol)?;
        self.model.set_edge_range(edge.into(), first, last)?;
        self.paths.insert((v, k), edge);
        Ok(edge)
    }

    /// The face edge `e` sweeps over stage `k`, used so that it bounds the sweep of the face
    /// in `W` beside it: that face runs along the edge's curve if `along`.
    fn side_face(&mut self, e: EdgeId, along: bool, k: usize) -> Result<Shape> {
        let (start, finish) = self.model.edge_end_vertices(e);
        let (start_path, end_path) = (self.vertex_path(start, k)?, self.vertex_path(finish, k)?);
        let top = self.edge_copy(e, k + 1)?;
        let e = self.edge_copy(e, k)?;
        let (curve, a, b) = self.model.edge_curve(Shape::from(e)).ok_or_else(|| Error::Failed("an edge without a curve".into()))?;
        let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
        // The surface; the (u, v) of the curve's point at s and moved to t; and the uses
        // round the face, counter-clockwise in (u, v).
        type Uses = Vec<(EdgeId, Orientation, Curve2, f64, f64)>;
        let (surface, uv, uses): (Surface, Box<dyn Fn(f64, f64) -> Point2>, Uses) = match self.steady {
            Steady::Slide { by } => {
                let length = by.length();
                match &curve {
                    Curve::Line(l) => {
                        let x = Vector::from(l.direction);
                        let across = by - x * by.dot(x);
                        if across.length() <= tolerance::LINEAR {
                            return Err(Error::Failed("a straight edge slid along itself bounding W".into()));
                        }
                        let plane = Plane::from_frame(Frame::new(l.origin, Direction::from_vector(x.cross(across)), l.direction));
                        let (du, dv) = (by.dot(x), across.length());
                        let path = |s: f64| (line(s, 0.0, du / length, dv / length), 0.0, length);
                        let (p0, p1) = (path(a), path(b));
                        (
                            Surface::Plane(plane),
                            Box::new(move |s, t| Point2::new(s + du * t, dv * t)),
                            vec![
                                (e, Orientation::Forward, line(0.0, 0.0, 1.0, 0.0), a, b),
                                (end_path, Orientation::Forward, p1.0, p1.1, p1.2),
                                (top, Orientation::Reversed, line(du, dv, 1.0, 0.0), a, b),
                                (start_path, Orientation::Reversed, p0.0, p0.1, p0.2),
                            ],
                        )
                    }
                    Curve::Circle(c) if Vector::from(c.frame().z()).cross(by).length() <= 1e-12 * length => {
                        let height = by.dot(Vector::from(c.frame().z()));
                        let uses = rectangle_uses(e, top, start_path, end_path, a, b, height, true);
                        (Surface::Cylinder(Cylinder::new(*c.frame(), c.radius())), Box::new(move |s, t| Point2::new(s, height * t)), uses)
                    }
                    other => {
                        // The extrusion as the B-spline surface it is: the curve's poles and
                        // their copies at the end, joined straight.
                        let profile = Profile::of(other, a, b, self.size)?;
                        let surface = profiled_surface(&profile.spline, &|q| Ok(vec![q, q + by]), KnotVector::new(vec![0.0, length], vec![2, 2], 1))?;
                        let uses = profile.uses(e, top, start_path, end_path, a, b, length);
                        (surface, Box::new(move |s, t| Point2::new(profile.at(s), length * t)), uses)
                    }
                }
            }
            Steady::Turn { axis, angle } => {
                let revolution = RevolutionSurface::new(curve.clone(), axis);
                let uses = rectangle_uses(e, top, start_path, end_path, a, b, angle, false);
                (Surface::Revolution(revolution), Box::new(move |s, t| Point2::new(angle * t, s)), uses)
            }
            Steady::Screw { .. } => {
                let profile = Profile::of(&curve, a, b, self.size)?;
                let motion = self.screw.as_ref().expect("a screw's motion");
                let knots = path(motion, Point::ORIGIN)?;
                let knots = KnotVector::new(knots.knots().to_vec(), knots.multiplicities().to_vec(), knots.degree());
                let swept = profiled_surface(&profile.spline, &|q| Ok(path(motion, q)?.poles().to_vec()), knots)?;
                let uses = profile.uses(e, top, start_path, end_path, a, b, 1.0);
                (swept, Box::new(move |s, t| Point2::new(profile.at(s), t)), uses)
            }
        };
        // Facing as the boundary of the sweep of the face beside it: along (the way that face
        // runs the edge) crossed with the motion.
        let (sm, tm) = (0.5 * (a + b), 0.5);
        let ds = 1e-6 * (b - a);
        let moved = |s: f64, t: f64| curve.point(s).transformed(&self.steady.at(t));
        let along_s = (moved(sm + ds, tm) - moved(sm - ds, tm)) * if along { 1.0 } else { -1.0 };
        let along_t = moved(sm, tm + 1e-6) - moved(sm, tm - 1e-6);
        let q = uv(sm, tm);
        let (_, su, sv) = surface.d1(q.x, q.y);
        let orientation = if su.cross(sv).dot(along_s.cross(along_t)) >= 0.0 { Orientation::Forward } else { Orientation::Reversed };
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, o, pcurve, first, last) in uses {
            let coedge = self.model.add_coedge(edge, o)?;
            self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = self.model.add_wire(&coedges)?;
        self.model.add(face, wire.into())?;
        Ok(Shape { orientation, ..face })
    }
}

/// The surface through the paths of `profile`'s poles (each pole's path's poles over `knots`,
/// the same for all): u the profile's parameter, v the motion's. The motion is rigid, so the
/// surface's iso-curves are the profile moved, exactly where the paths are.
fn profiled_surface(profile: &BSplineCurve, path_poles: &dyn Fn(Point) -> Result<Vec<Point>>, knots: KnotVector) -> Result<Surface> {
    let poles: Vec<Vec<Point>> = profile.poles().iter().map(|&q| path_poles(q)).collect::<Result<_>>()?;
    let u = if profile.is_periodic() {
        KnotVector::periodic(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
    } else {
        KnotVector::new(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
    };
    let rows = poles[0].len();
    let mut swept = BSplineSurface::try_new(poles, u, knots)?;
    if let Some(w) = profile.weights() {
        swept = swept.with_weights(w.iter().map(|&wi| vec![wi; rows]).collect());
    }
    Ok(Surface::BSpline(swept))
}

/// The uses round a side face whose surface has the curve's parameter `s` along one direction
/// and the motion (`0..extent` over the motion) along the other: `s_first` says whether `s` is
/// u. Counter-clockwise in (u, v).
#[allow(clippy::too_many_arguments)]
fn rectangle_uses(e: EdgeId, top: EdgeId, start_path: EdgeId, end_path: EdgeId, a: f64, b: f64, extent: f64, s_first: bool) -> Vec<(EdgeId, Orientation, Curve2, f64, f64)> {
    let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
    let (f, r) = (Orientation::Forward, Orientation::Reversed);
    // In (s, m): the edge at m = 0, the end path at s = b, the copy at m = extent, the start
    // path at s = a; counter-clockwise when s is u and the extent positive.
    let mut uses = if s_first {
        vec![
            (e, f, line(0.0, 0.0, 1.0, 0.0), a, b),
            (end_path, f, line(b, 0.0, 0.0, extent.signum()), 0.0, path_range(extent)),
            (top, r, line(0.0, extent, 1.0, 0.0), a, b),
            (start_path, r, line(a, 0.0, 0.0, extent.signum()), 0.0, path_range(extent)),
        ]
    } else {
        // (u the motion, v the curve's parameter.)
        vec![
            (e, f, line(0.0, 0.0, 0.0, 1.0), a, b),
            (end_path, f, line(0.0, b, extent.signum(), 0.0), 0.0, path_range(extent)),
            (top, r, line(extent, 0.0, 0.0, 1.0), a, b),
            (start_path, r, line(0.0, a, extent.signum(), 0.0), 0.0, path_range(extent)),
        ]
    };
    // Counter-clockwise: reversed where the rectangle is mirrored.
    let mirrored = (extent < 0.0) != !s_first;
    if mirrored {
        uses.reverse();
        for u in &mut uses {
            u.1 = u.1.reversed();
        }
    }
    uses
}

/// The parameter range of a vertex's path: its length along a slide (signed extents taken as
/// their size), the angle of a turn, the screw's fraction.
fn path_range(extent: f64) -> f64 {
    extent.abs()
}

/// The circle a turn about `axis` carries `p` round, its parameter the angle turned.
fn turn_circle(axis: crate::math::Axis, p: Point) -> Result<Circle> {
    let d = Vector::from(axis.direction);
    let foot = axis.origin + d * (p - axis.origin).dot(d);
    let radial = p - foot;
    let x = Direction::try_from_vector(radial).filter(|_| radial.length() > tolerance::LINEAR).ok_or_else(|| Error::InvalidInput("a point on the turn's axis".into()))?;
    Ok(Circle::new(Frame::new(foot, axis.direction, x), radial.length()))
}

/// The screw sampled as the pipes' motions are, finely enough that the path of each point
/// within `reach` of the axis (the quintics through its places, velocities and accelerations)
/// lies within 1e-10 of `size` of its helix.
fn screw_motion(axis: crate::math::Axis, angle: f64, rise: f64, size: f64, reach: f64) -> Result<PieceMotion> {
    let along = Vector::from(axis.direction);
    let at = move |f: f64| Transform::translation(along * (rise * f)) * Transform::rotation(axis, angle * f);
    let rates = move |f: f64, p: Point| -> (Vector, Vector) {
        let q = p.transformed(&at(f));
        let r = (q - axis.origin) - along * (q - axis.origin).dot(along);
        let turning = along.cross(r);
        (turning * angle + along * rise, along.cross(turning) * (angle * angle))
    };
    // (Probed at the farthest the body reaches from the axis.)
    let far = axis.origin + Vector::from(perpendicular(axis.direction)) * reach;
    let mut samples = 4usize.max((angle / (std::f64::consts::PI / 8.0)).ceil() as usize);
    loop {
        let params: Vec<f64> = (0..=samples).map(|i| i as f64 / samples as f64).collect();
        let transforms: Vec<Transform> = params.iter().map(|&f| at(f)).collect();
        let motion = PieceMotion { params, transforms, rates: Some(Box::new(rates)) };
        let curve = path(&motion, far)?;
        let error = (0..4 * samples).map(|i| (i as f64 + 0.5) / (4 * samples) as f64).map(|f| curve.point(f).distance(far.transformed(&at(f)))).fold(0.0, f64::max);
        if error <= 1e-10 * size || samples >= 4096 {
            return Ok(motion);
        }
        samples *= 2;
    }
}

fn perpendicular(d: Direction) -> Direction {
    let v = Vector::from(d);
    let other = if v.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    Direction::from_vector(v.cross(other))
}

/// An edge's curve over `[a, b]` as a B-spline: exactly where it converts so (a circle's
/// arcs as rational quadratics), else cubics through its points within 1e-10 of the size. Where
/// the B-spline's parameter is not the edge's, `sigma` maps the edge's to it.
struct Profile {
    spline: BSplineCurve,
    /// The B-spline's parameter (as the first coordinate) against the edge's, within 1e-12 of
    /// the size in space; `None` where they are one.
    sigma: Option<crate::geometry::curve2::BSplineCurve2>,
}

impl Profile {
    fn of(curve: &Curve, a: f64, b: f64, size: f64) -> Result<Profile> {
        let exact = crate::geometry::convert::curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(curve.clone(), a, b)));
        if let Some(bs) = exact {
            let same = (bs.first_parameter() - a).abs() <= 1e-12 * (1.0 + a.abs()) && (bs.last_parameter() - b).abs() <= 1e-12 * (1.0 + b.abs());
            let keeps = (0..=16).all(|i| {
                let t = a + (b - a) * i as f64 / 16.0;
                curve.point(t).distance(bs.point(t)) <= 1e-12 * (1.0 + size)
            });
            if same && keeps {
                return Ok(Profile { spline: bs, sigma: None });
            }
            if let Some(sigma) = reparametrization(curve, a, b, &bs, size) {
                return Ok(Profile { spline: bs, sigma: Some(sigma) });
            }
        }
        let mut n = 32;
        loop {
            let params: Vec<f64> = (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect();
            let points: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
            let fit = crate::geometry::fit::interpolate_with_parameters(&points, &params)?;
            let error = (0..n).map(|i| a + (b - a) * (i as f64 + 0.5) / n as f64).map(|t| curve.point(t).distance(fit.point(t))).fold(0.0, f64::max);
            if error <= 1e-10 * size || n >= 8192 {
                return Ok(Profile { spline: fit, sigma: None });
            }
            n *= 2;
        }
    }

    /// The B-spline's parameter at the edge's `s`.
    fn at(&self, s: f64) -> f64 {
        match &self.sigma {
            Some(sigma) => sigma.point(s).x,
            None => s,
        }
    }

    /// The uses round the side face of a surface with the B-spline's parameter along u and the
    /// motion (`0..extent`) along v ([`rectangle_uses`], the edges' pcurves mapped).
    #[allow(clippy::too_many_arguments)]
    fn uses(&self, e: EdgeId, top: EdgeId, start_path: EdgeId, end_path: EdgeId, a: f64, b: f64, extent: f64) -> Vec<(EdgeId, Orientation, Curve2, f64, f64)> {
        let mut uses = rectangle_uses(e, top, start_path, end_path, a, b, extent, true);
        if let Some(sigma) = &self.sigma {
            let (ua, ub) = (self.at(a), self.at(b));
            let raised = {
                let poles = sigma.poles().iter().map(|q| Point2::new(q.x, q.y + extent)).collect();
                crate::geometry::curve2::BSplineCurve2::new(poles, sigma.knots().to_vec(), sigma.multiplicities().to_vec(), sigma.degree())
            };
            let line = |x: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, 0.0), Direction2::new(0.0, dy)));
            for u in &mut uses {
                if u.0 == e {
                    u.2 = Curve2::BSpline(sigma.clone());
                } else if u.0 == top {
                    u.2 = Curve2::BSpline(raised.clone());
                } else if u.0 == start_path {
                    u.2 = line(ua, extent.signum());
                } else if u.0 == end_path {
                    u.2 = line(ub, extent.signum());
                }
            }
        }
        uses
    }
}

/// The map from `curve`'s parameter over `[a, b]` to that of `spline` (the same points), as
/// cubics through samples, within 1e-12 of `size` in space, in pieces between the points of
/// the spline's knots (the map is smooth only between them); `None` if it does not settle.
fn reparametrization(curve: &Curve, a: f64, b: f64, spline: &BSplineCurve, size: f64) -> Option<crate::geometry::curve2::BSplineCurve2> {
    let (u0, u1) = (spline.first_parameter(), spline.last_parameter());
    let hair = 1e-12 * (1.0 + size);
    // Newton on the distance: the parameter on one curve of the other's point.
    let newton = |target: Point, guess: f64, lo: f64, hi: f64, on: &dyn Fn(f64) -> (Point, Vector)| -> Option<f64> {
        let mut x = guess.clamp(lo, hi);
        for _ in 0..50 {
            let (p, d) = on(x);
            let dd = d.dot(d);
            if !(dd > 0.0) {
                return None;
            }
            let step = (target - p).dot(d) / dd;
            x = (x + step).clamp(lo, hi);
            if step.abs() <= 1e-15 * (hi - lo) {
                break;
            }
        }
        (on(x).0.distance(target) <= hair).then_some(x)
    };
    let on_spline = |u: f64| spline.d1(u);
    let on_curve = |s: f64| curve.d1(s);
    // The pieces: between the edge's parameters of the spline's knots.
    let mut breaks = vec![(a, u0)];
    for &k in spline.knots() {
        if k > u0 + 1e-12 * (u1 - u0) && k < u1 - 1e-12 * (u1 - u0) {
            let guess = a + (b - a) * (k - u0) / (u1 - u0);
            breaks.push((newton(spline.point(k), guess, a, b, &on_curve)?, k));
        }
    }
    breaks.push((b, u1));
    let mut pieces: Vec<crate::geometry::curve2::BSplineCurve2> = Vec::new();
    for w in breaks.windows(2) {
        let ((sa, ua), (sb, ub)) = (w[0], w[1]);
        let mut n = 8;
        let piece = loop {
            if n > 4096 {
                return None;
            }
            let params: Vec<f64> = (0..=n).map(|i| sa + (sb - sa) * i as f64 / n as f64).collect();
            let mut us: Vec<f64> = Vec::with_capacity(params.len());
            for (i, &x) in params.iter().enumerate() {
                let u = if i == 0 {
                    ua
                } else if i == n {
                    ub
                } else {
                    newton(curve.point(x), us[i - 1] + (ub - ua) / n as f64, ua, ub, &on_spline)?
                };
                us.push(u);
            }
            let points: Vec<Point2> = us.iter().map(|&u| Point2::new(u, 0.0)).collect();
            let sigma = crate::geometry::fit::interpolate_with_parameters(&points, &params).ok()?;
            let error = (0..n).map(|i| sa + (sb - sa) * (i as f64 + 0.5) / n as f64).map(|x| curve.point(x).distance(spline.point(sigma.point(x).x.clamp(u0, u1)))).fold(0.0, f64::max);
            if error <= hair {
                break sigma;
            }
            n *= 2;
        };
        pieces.push(piece);
    }
    joined_pieces(&pieces)
}

/// Cubic B-splines end to end (each clamped, each starting where the one before ends) as one,
/// continuous across the joins.
fn joined_pieces(pieces: &[crate::geometry::curve2::BSplineCurve2]) -> Option<crate::geometry::curve2::BSplineCurve2> {
    let degree = pieces.first()?.degree();
    if pieces.iter().any(|p| p.degree() != degree) {
        return None;
    }
    let (mut poles, mut knots, mut multiplicities): (Vec<Point2>, Vec<f64>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
    for (k, piece) in pieces.iter().enumerate() {
        let (pk, pm) = (piece.knots(), piece.multiplicities());
        let skip = if k == 0 { 0 } else { 1 };
        poles.extend_from_slice(&piece.poles()[skip..]);
        for (i, (&x, &m)) in pk.iter().zip(pm).enumerate() {
            if i == 0 && k > 0 {
                // (The join: as often as the degree.)
                *multiplicities.last_mut()? = degree;
                continue;
            }
            knots.push(x);
            multiplicities.push(m);
        }
    }
    crate::geometry::curve2::BSplineCurve2::try_new(poles, knots, multiplicities, degree).ok()
}
