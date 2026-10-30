//! Sweeping a profile through space: extrusion, revolution and pipes.
//!
//! The dimension goes up by one: a vertex sweeps an edge, an edge a face, a face a solid.

mod corner;
mod draft;
pub use draft::extrude_drafted;
mod middle;
mod pipe;
pub use middle::middle_path;
mod pipe_shell;

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, TAU};

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::project;
use crate::geometry::surface::{Cone, Cylinder, ExtrusionSurface, Plane, RevolutionSurface, Sphere, Surface, Torus};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Axis, Direction, Direction2, Frame, Point, Point2, Transform, Transformable, Vector};
use crate::properties;
use crate::tolerance::{self, ANGULAR};
use crate::topology::{EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape, VertexId, WireId};

/// The result of a sweep, with the history linking profile sub-shapes to what they swept.
#[derive(Debug, Clone, PartialEq)]
pub struct Sweep {
    pub shape: Shape,
    generated: Vec<(NodeId, Vec<Shape>)>,
    /// Each profile face's first and last faces in the result (the face itself, closing the
    /// start, and its copy closing the end), where the sweep made them.
    ends: Vec<(NodeId, Shape, Shape)>,
    /// Each profile edge and vertex with its copy at the end (an extrusion's).
    copies: Vec<(NodeId, Shape)>,
}

impl Sweep {
    /// The same sweep with every resulting shape moved by `t`.
    fn moved(self, t: &Transform) -> Sweep {
        let generated = self.generated.into_iter().map(|(n, s)| (n, s.iter().map(|s| s.moved(t)).collect())).collect();
        let ends = self.ends.into_iter().map(|(n, a, b)| (n, a.moved(t), b.moved(t))).collect();
        let copies = self.copies.into_iter().map(|(n, c)| (n, c.moved(t))).collect();
        Sweep { shape: self.shape.moved(t), generated, ends, copies }
    }

    /// The shapes swept by `input` (a sub-shape of the profile); empty if `input` is not
    /// part of the profile.
    pub fn generated(&self, input: Shape) -> Vec<Shape> {
        self.generated.iter().filter(|(node, _)| *node == input.node).flat_map(|(_, s)| s.clone()).collect()
    }

    /// Every profile sub-shape (by node) with the shapes it swept.
    pub fn history(&self) -> &[(NodeId, Vec<Shape>)] {
        &self.generated
    }

    /// Each profile face (by node) with the result's first and last faces made from it.
    pub fn ends(&self) -> &[(NodeId, Shape, Shape)] {
        &self.ends
    }

    /// Each profile edge and vertex (by node) with its copy at the end of an extrusion.
    pub fn copies(&self) -> &[(NodeId, Shape)] {
        &self.copies
    }
}

/// Sweeps `profile` along the vector `v`.
///
/// A vertex sweeps a straight edge, an edge a face, a wire a shell and a face a solid,
/// closed by the profile at the start and a translated copy at the end; a compound a
/// compound of its parts' sweeps. Straight edges
/// sweep planes and circles swept along their axis cylinders; other edges are not
/// supported yet.
#[doc(alias = "BRepPrimAPI_MakePrism")]
pub fn extrude(model: &mut Model, profile: Shape, v: Vector) -> Result<Sweep> {
    if !(v.length() > tolerance::LINEAR) || !v.length().is_finite() {
        return Err(Error::InvalidInput(format!("extrusion vector must be non-zero, got {v:?}")));
    }
    if profile.placement != Transform::IDENTITY {
        // Sweep the shape where it is defined, then place the result like the profile.
        let local = v.transformed(&profile.placement.inverse());
        let sweep = extrude(model, Shape { placement: Transform::IDENTITY, ..profile }, local)?;
        return Ok(sweep.moved(&profile.placement));
    }
    let mut x = Extruder {
        model,
        v,
        translation: Transform::translation(v),
        vertex_paths: HashMap::new(),
        end_vertices: HashMap::new(),
        end_edges: HashMap::new(),
        faces: HashMap::new(),
        generated: Vec::new(),
        ends: Vec::new(),
        copies: Vec::new(),
    };
    let shape = match profile.kind() {
        NodeKind::Vertex => {
            let vertex = VertexId::from_index(profile.node.index);
            x.check(vertex)?;
            Shape::from(x.vertex_path(vertex)?)
        }
        NodeKind::Edge => {
            let edge = EdgeId::from_index(profile.node.index);
            x.check(edge)?;
            let (face, orientation) =
                x.edge_face(edge)?.ok_or_else(|| Error::InvalidInput("a degenerate edge sweeps nothing".into()))?;
            Shape { orientation: orientation.compose(profile.orientation), ..Shape::from(face) }
        }
        NodeKind::Wire => {
            let wire = WireId::from_index(profile.node.index);
            x.check(wire)?;
            let shell = x.model.new_empty(NodeKind::Shell);
            let coedges = x.model[wire].coedges.clone();
            for coedge in coedges {
                let c = &x.model[coedge];
                let (edge, orientation) = (c.edge, c.orientation.compose(profile.orientation));
                if let Some((face, canonical)) = x.edge_face(edge)? {
                    x.model.add(shell, Shape { orientation: canonical.compose(orientation), ..Shape::from(face) })?;
                }
            }
            shell
        }
        NodeKind::Face => {
            let face = FaceId::from_index(profile.node.index);
            x.check(face)?;
            x.face_solid(face, profile.orientation)?
        }
        // A compound: each of its parts swept, as a compound.
        NodeKind::Compound => {
            let parts: Vec<Shape> = x.model.children(profile).collect();
            let compound = x.model.new_empty(NodeKind::Compound);
            let (mut generated, mut ends, mut copies) = (Vec::new(), Vec::new(), Vec::new());
            for part in parts {
                let sweep = extrude(x.model, part, v)?;
                generated.extend(sweep.generated);
                ends.extend(sweep.ends);
                copies.extend(sweep.copies);
                x.model.add(compound, sweep.shape)?;
            }
            return Ok(Sweep { shape: compound, generated, ends, copies });
        }
        kind => todo!("extruding a {kind:?}"),
    };
    let generated = std::mem::take(&mut x.generated);
    let ends = std::mem::take(&mut x.ends);
    let copies = std::mem::take(&mut x.copies);
    Ok(Sweep { shape, generated, ends, copies })
}

/// Sweeps `profile` along the `spine` (an edge, or a wire of edges meeting smoothly),
/// keeping its relation to the spine's tangent (OCCT's corrected Frenet trihedron: the
/// profile turns with the tangent and no more): along a straight spine the profile is
/// translated, along a circular arc it turns about the arc's axis, along other curves its
/// faces sweep B-spline surfaces.
#[doc(alias = "BRepOffsetAPI_MakePipe")]
pub fn pipe(model: &mut Model, profile: Shape, spine: Shape) -> Result<Sweep> {
    pipe_framed(model, profile, spine, &PipeOptions::default())
}

/// [`pipe`], the profile carried in the spine's Frenet frame where `options.frenet` (turning
/// with its principal normal, as about a helix's axis; where the spine runs straight,
/// turning no more than its tangent), or as `options.guide` turns it.
fn pipe_framed(model: &mut Model, profile: Shape, spine: Shape, options: &PipeOptions) -> Result<Sweep> {
    use crate::geometry::ParametricCurve;
    let pieces = spine_pieces(model, spine)?;
    let framing = framing(model, &pieces, options)?;
    if let ([piece], pipe::Framing::Minimal | pipe::Framing::Frenet) = (&pieces[..], &framing) {
        let (curve, first, last) = (&piece.curve, piece.first, piece.last);
        // (From the profile carried back to the spine's start from where it stands on it,
        // as the pipe carries it: by OCCT's placement an end, or another point along it.)
        let at = || -> Result<f64> {
            let (_, s) = pipe::section_placement(model, profile, &pieces)?;
            Ok(first + (last - first) * s)
        };
        match curve {
            Curve::Line(_) => {
                let back = curve.point(first) - curve.point(at()?);
                return extrude(model, profile.moved(&Transform::translation(back)), curve.point(last) - curve.point(first));
            }
            Curve::Circle(c) => {
                let axis = c.frame().z_axis();
                let axis = if last >= first { axis } else { Axis::new(axis.origin, -axis.direction) };
                let back = Transform::rotation(axis, -(at()? - first).abs());
                return revolve(model, profile.moved(&back), axis, (last - first).abs());
            }
            _ => {}
        }
    }
    if profile.placement != Transform::IDENTITY {
        // Sweep the shape where it is defined, then place the result like the profile.
        let inverse = profile.placement.inverse();
        let local: Vec<pipe::SpinePiece> = pieces.iter().map(|p| pipe::SpinePiece { curve: p.curve.transformed(&inverse), first: p.first, last: p.last }).collect();
        let local_options = PipeOptions { guide: options.guide.map(|g| g.moved(&inverse)), ..*options };
        let framing = self::framing(model, &local, &local_options)?;
        let sweep = pipe::pipe_along(model, Shape { placement: Transform::IDENTITY, ..profile }, &local, &framing)?;
        return Ok(sweep.moved(&profile.placement));
    }
    pipe::pipe_along(model, profile, &pieces, &framing)
}

/// The frames `options` carry a profile in along the spine's `pieces`.
fn framing(model: &Model, pieces: &[pipe::SpinePiece], options: &PipeOptions) -> Result<pipe::Framing> {
    Ok(match options.guide {
        Some(Guide::Wire { wire, curvilinear }) => pipe::Framing::Guide(pipe::GuideLaw::new(pieces, &spine_pieces(model, wire)?, curvilinear)?),
        Some(Guide::Fixed(d)) => pipe::Framing::Fixed(Vector::from(d)),
        None if options.frenet => pipe::Framing::Frenet,
        None => pipe::Framing::Minimal,
    })
}

/// The edges of a pipe's spine, in order, each run in the direction the spine takes it.
fn spine_pieces(model: &Model, spine: Shape) -> Result<Vec<pipe::SpinePiece>> {
    let edges: Vec<Shape> = match spine.kind() {
        NodeKind::Edge => vec![spine],
        NodeKind::Wire => model.explore(spine, NodeKind::Edge).collect(),
        kind => return Err(Error::InvalidInput(format!("a pipe's spine is an edge or a wire, not a {kind:?}"))),
    };
    let mut pieces = Vec::with_capacity(edges.len());
    for edge in edges {
        let (curve, first, last) = model.edge_curve(edge).ok_or_else(|| Error::InvalidInput("a pipe's spine has no curve".into()))?;
        let (first, last) = if edge.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
        pieces.push(pipe::SpinePiece { curve, first, last });
    }
    Ok(pieces)
}

/// What turns a pipe's profile about its spine, in place of the spine alone.
#[doc(alias = "BRepOffsetAPI_MakePipeShell::SetMode")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Guide {
    /// A wire (or edge) the profile keeps facing: the frame's second axis points from the
    /// spine's point to the guide's point dividing its length as the spine's point divides
    /// the spine's (`curvilinear`), else where the plane square to the spine there meets
    /// the guide.
    #[doc(alias = "GeomFill_GuideTrihedronAC")]
    Wire { wire: Shape, curvilinear: bool },
    /// No turning at all: the profile carried along the spine by translations alone, the
    /// frame's first axis this direction (what `correction` turns the profile square to).
    #[doc(alias = "GeomFill_Fixed")]
    Fixed(Direction),
}

impl Guide {
    fn moved(self, t: &Transform) -> Guide {
        match self {
            Guide::Wire { wire, curvilinear } => Guide::Wire { wire: wire.moved(t), curvilinear },
            Guide::Fixed(d) => Guide::Fixed(d.transformed(t)),
        }
    }
}

/// How a pipe places its profile on the spine, beyond keeping it where it is.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PipeOptions {
    /// Moves the profile onto the spine: its point nearest the spine's start to where the
    /// profile stands on the spine.
    #[doc(alias = "WithContact")]
    pub contact: bool,
    /// Turns the profile about its centre to be square to the spine there.
    #[doc(alias = "WithCorrection")]
    pub correction: bool,
    /// Carries the profile in the spine's Frenet frame (see [`pipe`]'s frame otherwise).
    #[doc(alias = "BRepOffsetAPI_MakePipeShell::SetMode")]
    pub frenet: bool,
    /// What else turns the profile (in place of `frenet`).
    pub guide: Option<Guide>,
}

/// [`pipe`], the profile first placed on the spine as `options` say.
#[doc(alias = "BRepOffsetAPI_MakePipeShell")]
pub fn pipe_with(model: &mut Model, profile: Shape, spine: Shape, options: &PipeOptions) -> Result<Sweep> {
    if !options.contact && !options.correction {
        return pipe_framed(model, profile, spine, options);
    }
    let pieces = spine_pieces(model, spine)?;
    let placing = placing(model, profile, &pieces, options)?;
    pipe_framed(model, profile.moved(&placing), spine, options)
}

/// How a pipe goes round a corner of its spine (where the tangent jumps).
#[doc(alias = "BRepBuilderAPI_TransitionMode")]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Corner {
    /// The profile is carried on unturned (and the pipe's surfaces are not cut or filled):
    /// what [`pipe`] does.
    #[default]
    #[doc(alias = "BRepBuilderAPI_Transformed")]
    Unturned,
    /// The profile turns with the spine, and the sweeps on either side meet in the plane
    /// bisecting the corner.
    #[doc(alias = "BRepBuilderAPI_RightCorner")]
    Mitred,
    /// The profile turns with the spine, turned about the corner to fill its outside.
    #[doc(alias = "BRepBuilderAPI_RoundCorner")]
    Rounded,
}

/// [`pipe_with`], going round the spine's corners as `corner` says: the profile turned at
/// each by the least rotation taking the tangent before it onto the tangent after it. A
/// face sweeps a solid, a closed planar wire the shell bounding the solid its face would
/// sweep, less the ends. Along lines and arcs for mitred corners.
#[doc(alias = "BRepOffsetAPI_MakePipeShell::SetTransitionMode")]
pub fn pipe_cornered(model: &mut Model, profile: Shape, spine: Shape, options: &PipeOptions, corner: Corner) -> Result<Shape> {
    if corner == Corner::Unturned {
        return Ok(pipe_with(model, profile, spine, options)?.shape);
    }
    let pieces = spine_pieces(model, spine)?;
    let placing = placing(model, profile, &pieces, options)?;
    corner::pipe_cornered(model, profile.moved(&placing), &pieces, corner)
}

/// Sweeps several `sections` (wires of as many edges, or vertices, each placed as its
/// options say) along `spine`: each section seen in the spine's frame where it stands,
/// the sections blended along the spine (linearly between them, carried on past the first
/// and the last) and carried by its frame. A shell, or with `solid` a solid closed by
/// planar faces at the ends.
#[doc(alias = "BRepOffsetAPI_MakePipeShell")]
pub fn pipe_through(model: &mut Model, sections: &[(Shape, PipeOptions)], spine: Shape, solid: bool) -> Result<Shape> {
    pipe_through_with(model, sections, spine, solid, None).map(|s| s.shape)
}

/// [`pipe_through`], the frames turned by `guide` (if any), with its history: the faces
/// each section edge swept (or part of one, where sections were split to match), and the
/// first section's (as a profile face's) first and last faces, the caps of a solid.
pub fn pipe_through_with(model: &mut Model, sections: &[(Shape, PipeOptions)], spine: Shape, solid: bool, guide: Option<Guide>) -> Result<Sweep> {
    let pieces = spine_pieces(model, spine)?;
    let framing = framing(model, &pieces, &PipeOptions { guide, ..Default::default() })?;
    // (Each section corrected square to the guide's fixed direction, where it has one.)
    let sections: Vec<(Shape, PipeOptions)> = sections.iter().map(|&(s, o)| (s, PipeOptions { guide: guide.or(o.guide), ..o })).collect();
    pipe_shell::pipe_shell(model, &sections, &pieces, solid, &framing)
}

/// The motion placing `profile` on the spine as `options` say (contact, correction).
fn placing(model: &Model, profile: Shape, pieces: &[pipe::SpinePiece], options: &PipeOptions) -> Result<Transform> {
    use crate::geometry::ParametricCurve;
    let mut placing = Transform::IDENTITY;
    if !options.contact && !options.correction {
        return Ok(placing);
    }
    let (k, s) = pipe::section_placement(model, profile, pieces)?;
    let piece = &pieces[k];
    let t = piece.first + (piece.last - piece.first) * s;
    let (at, d) = piece.curve.d1(t);
    // (Square to a fixed frame's direction where the frame is fixed.)
    let tangent = match options.guide {
        Some(Guide::Fixed(d)) => Vector::from(d),
        _ if piece.last < piece.first => -d,
        _ => d,
    };
    let tangent = tangent * (1.0 / tangent.length());
    let (center, normal) = pipe::profile_axis(model, profile);
    if options.correction {
        if let Some(n) = normal {
            let axis = n.cross(tangent);
            if axis.length() > ANGULAR {
                let angle = n.dot(tangent).clamp(-1.0, 1.0).acos();
                placing = Transform::rotation(Axis::new(center, Direction::from_vector(axis)), angle);
            }
        }
    }
    if options.contact {
        // The profile's point nearest the spine's start, moved to where it stands.
        let start = pieces[0].curve.point(pieces[0].first);
        let mut nearest = (f64::INFINITY, center);
        for e in model.explore(profile, NodeKind::Edge).collect::<Vec<_>>() {
            if let Some((c, a, b)) = model.edge_curve(e) {
                for i in 0..=200 {
                    let p = c.point(a + (b - a) * i as f64 / 200.0);
                    // (Of points equally near, the last: where it is a tie, OCCT's choice
                    // comes down to rounding; this one matches its results.)
                    if p.distance(start) <= nearest.0 {
                        nearest = (p.distance(start), p);
                    }
                }
            }
        }
        if model.explore(profile, NodeKind::Edge).next().is_none() {
            nearest = (0.0, center);
        }
        let p = nearest.1.transformed(&placing);
        placing = Transform::translation(at - p) * placing;
    }
    Ok(placing)
}

/// Sweeps `profile` around `axis` by `angle` (radians). Where possible the swept surfaces
/// are recognized as elementary ones (a circle revolved about a coplanar axis gives a
/// torus, not a general surface of revolution).
#[doc(alias = "BRepPrimAPI_MakeRevol")]
///
/// A vertex sweeps an edge (a circular arc, or a degenerate edge if it lies on the axis),
/// an edge a face, a wire a shell and a face a solid. A partial revolution of a face is
/// closed by the profile itself at the start and by a rotated copy at the end; a full one
/// has no caps (the swept faces close on themselves along seams). Edges lying on
/// the axis sweep nothing. A negative angle revolves the other way.
pub fn revolve(model: &mut Model, profile: Shape, axis: Axis, angle: f64) -> Result<Sweep> {
    if !(angle.is_finite() && angle != 0.0 && angle.abs() <= TAU + tolerance::PARAMETRIC) {
        return Err(Error::InvalidInput(format!("revolution angle must be non-zero and at most 2π in size, got {angle}")));
    }
    let (axis, angle) =
        if angle < 0.0 { (Axis::new(axis.origin, -axis.direction), -angle) } else { (axis, angle.min(TAU)) };
    if profile.placement != Transform::IDENTITY {
        // Sweep the shape where it is defined, then place the result like the profile.
        let local = axis.transformed(&profile.placement.inverse());
        let sweep = revolve(model, Shape { placement: Transform::IDENTITY, ..profile }, local, angle)?;
        return Ok(sweep.moved(&profile.placement));
    }
    let mut r = Revolver {
        model,
        axis,
        angle,
        full: angle >= TAU - tolerance::PARAMETRIC,
        rotation: Transform::rotation(axis, angle),
        vertex_paths: HashMap::new(),
        end_vertices: HashMap::new(),
        end_edges: HashMap::new(),
        faces: HashMap::new(),
        generated: Vec::new(),
        ends: Vec::new(),
        copies: Vec::new(),
    };
    let shape = match profile.kind() {
        NodeKind::Vertex => {
            let vertex = VertexId::from_index(profile.node.index);
            r.check(vertex)?;
            Shape::from(r.vertex_path(vertex)?)
        }
        NodeKind::Edge => {
            let edge = EdgeId::from_index(profile.node.index);
            r.check(edge)?;
            let (face, orientation) = r
                .edge_face(edge)?
                .ok_or_else(|| Error::InvalidInput("the edge lies on the axis".into()))?;
            Shape { orientation: orientation.compose(profile.orientation), ..Shape::from(face) }
        }
        NodeKind::Wire => {
            let wire = WireId::from_index(profile.node.index);
            r.check(wire)?;
            let shell = r.model.new_empty(NodeKind::Shell);
            let coedges = r.model[wire].coedges.clone();
            for coedge in coedges {
                let c = &r.model[coedge];
                let (edge, orientation) = (c.edge, c.orientation.compose(profile.orientation));
                if let Some((face, canonical)) = r.edge_face(edge)? {
                    r.model.add(shell, Shape { orientation: canonical.compose(orientation), ..Shape::from(face) })?;
                }
            }
            shell
        }
        NodeKind::Face => {
            let face = FaceId::from_index(profile.node.index);
            r.check(face)?;
            r.face_solid(face, profile.orientation)?
        }
        kind => todo!("revolving a {kind:?}"),
    };
    let generated = std::mem::take(&mut r.generated);
    let (ends, copies) = (std::mem::take(&mut r.ends), std::mem::take(&mut r.copies));
    Ok(Sweep { shape, generated, ends, copies })
}

/// The state of one revolution: what each profile vertex and edge has swept, and their
/// rotated copies at the end of a partial revolution.
struct Revolver<'m> {
    model: &'m mut Model,
    axis: Axis,
    angle: f64,
    full: bool,
    rotation: Transform,
    vertex_paths: HashMap<VertexId, EdgeId>,
    end_vertices: HashMap<VertexId, VertexId>,
    end_edges: HashMap<EdgeId, EdgeId>,
    /// The face swept by each edge (none for edges on the axis or degenerate edges), with
    /// the orientation that makes its normal `∂u × ∂v` of the sweep.
    faces: HashMap<EdgeId, Option<(FaceId, Orientation)>>,
    generated: Vec<(NodeId, Vec<Shape>)>,
    ends: Vec<(NodeId, Shape, Shape)>,
    copies: Vec<(NodeId, Shape)>,
}

/// How a swept edge's face is parametrized: an elementary surface with `u` the rotation
/// angle and `v = v0 + sigma·t` along the edge's curve, or a plane (perpendicular to the
/// axis) whose pcurves are computed from the 3D curves.
enum SweptSurface {
    Ruled { surface: Surface, v0: f64, sigma: f64 },
    Planar(Surface),
}

impl Revolver<'_> {
    fn check<T: crate::topology::Node>(&self, id: crate::topology::Id<T>) -> Result<()> {
        self.model.get(id).map(|_| ()).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))
    }

    fn tolerance(&self) -> f64 {
        tolerance::LINEAR
    }

    /// The foot of `p` on the axis and the vector from it to `p`.
    fn radial(&self, p: Point) -> (Point, Vector) {
        let z = Vector::from(self.axis.direction);
        let foot = self.axis.origin + z * (p - self.axis.origin).dot(z);
        (foot, p - foot)
    }

    fn on_axis(&self, p: Point) -> bool {
        self.radial(p).1.length() <= self.tolerance()
    }

    fn point(&self, v: VertexId) -> Point {
        self.model[v].point
    }

    fn edge_vertices(&self, e: EdgeId) -> (VertexId, VertexId) {
        let edge = &self.model[e];
        (self.model[edge.start].vertex, self.model[edge.end].vertex)
    }

    fn record(&mut self, input: impl Into<NodeId>, output: Shape) {
        self.generated.push((input.into(), vec![output]));
    }

    /// The vertex at the end of the revolution (the same one for a full revolution or a
    /// vertex on the axis).
    fn end_vertex(&mut self, v: VertexId) -> VertexId {
        if self.full || self.on_axis(self.point(v)) {
            return v;
        }
        if let Some(&end) = self.end_vertices.get(&v) {
            return end;
        }
        let p = self.point(v).transformed(&self.rotation);
        let end = self.model.add_vertex(p, self.model[v].tolerance);
        self.end_vertices.insert(v, end);
        self.copies.push((v.into(), Shape::from(end)));
        end
    }

    /// The arc swept by a vertex (degenerate on the axis), from it to its end copy,
    /// parametrized by the rotation angle. Arcs are shared by the faces meeting along them;
    /// a degenerate edge is not: like OCCT's `BRepPrimAPI_MakeRevol`, each call (each face)
    /// gets its own, as a degenerate edge shared by several faces is not valid.
    fn vertex_path(&mut self, v: VertexId) -> Result<EdgeId> {
        if let Some(&e) = self.vertex_paths.get(&v) {
            return Ok(e);
        }
        let p = self.point(v);
        let end = self.end_vertex(v);
        let tol = self.tolerance();
        // (At least the vertex's tolerance, as OCCT's sweeps give the edges they generate.)
        let edge_tol = self.model[v].tolerance.max(tol);
        let edge = self.model.add_edge(v, end, edge_tol)?;
        let (center, radial) = self.radial(p);
        if let Some(x) = Direction::try_from_vector(radial).filter(|_| radial.length() > tol) {
            let circle = Circle::new(Frame::new(center, self.axis.direction, x), radial.length());
            self.model.set_edge_curve(edge.into(), Curve::Circle(circle), edge_tol)?;
            self.model.set_edge_range(edge.into(), 0.0, self.angle)?;
            self.vertex_paths.insert(v, edge);
        }
        self.record(v, edge.into());
        Ok(edge)
    }

    /// True for a straight edge lying on the axis (it sweeps nothing).
    fn is_on_axis(&self, e: EdgeId) -> bool {
        let (a, b) = self.edge_vertices(e);
        let straight = matches!(self.model.get(self.model[e].curve).map(|s| &s.curve), Some(Curve::Line(_)));
        straight && self.on_axis(self.point(a)) && self.on_axis(self.point(b))
    }

    /// The edge at the end of the revolution: a rotated copy, or the edge itself for a
    /// full revolution or an edge on the axis.
    fn end_edge(&mut self, e: EdgeId) -> Result<EdgeId> {
        if self.full || self.is_on_axis(e) {
            return Ok(e);
        }
        if let Some(&end) = self.end_edges.get(&e) {
            return Ok(end);
        }
        let (a, b) = self.edge_vertices(e);
        let (a, b) = (self.end_vertex(a), self.end_vertex(b));
        let tol = self.model[e].tolerance;
        let end = self.model.add_edge(a, b, tol)?;
        if let Some(slot) = self.model.get(self.model[e].curve) {
            let (curve, first, last) = (slot.curve.transformed(&self.rotation), slot.first, slot.last);
            self.model.set_edge_curve(end.into(), curve, tol)?;
            self.model.set_edge_range(end.into(), first, last)?;
        }
        self.end_edges.insert(e, end);
        self.copies.push((e.into(), Shape::from(end)));
        Ok(end)
    }

    /// The surface swept by `curve`, parametrized as described at [`SweptSurface`].
    fn swept_surface(&self, curve: &Curve, first: f64, last: f64) -> SweptSurface {
        let z = Vector::from(self.axis.direction);
        let tol = self.tolerance();
        let general = || SweptSurface::Ruled {
            surface: Surface::Revolution(RevolutionSurface::new(curve.clone(), self.axis)),
            v0: 0.0,
            sigma: 1.0,
        };
        match curve {
            Curve::Line(line) => {
                let d = Vector::from(line.direction);
                let (foot, radial) = self.radial(line.origin);
                if d.cross(z).length() <= ANGULAR {
                    // Parallel to the axis: a cylinder, with v the height above `foot`.
                    let Some(x) = Direction::try_from_vector(radial) else { return general() };
                    let cylinder = Cylinder::new(Frame::new(foot, self.axis.direction, x), radial.length());
                    return SweptSurface::Ruled { surface: Surface::Cylinder(cylinder), v0: 0.0, sigma: d.dot(z) };
                }
                if d.dot(z).abs() <= ANGULAR {
                    return SweptSurface::Planar(Surface::Plane(Plane::from_frame(Frame::from_z(foot, self.axis.direction))));
                }
                if (line.origin - self.axis.origin).dot(d.cross(z)).abs() > tol {
                    return general(); // skew: a hyperboloid
                }
                // Coplanar with the axis: a cone around it, with v along the line (pointing
                // away from the axis direction's negative side).
                let x = match Direction::try_from_vector(radial).filter(|_| radial.length() > tol) {
                    Some(x) => x,
                    None => match Direction::try_from_vector(self.radial(line.origin + d).1) {
                        Some(x) => x,
                        None => return general(),
                    },
                };
                let sigma = d.dot(z).signum();
                let g = d * sigma;
                let semi_angle = g.dot(Vector::from(x)).atan2(g.dot(z));
                let cone = Cone::new(Frame::new(foot, self.axis.direction, x), semi_angle, radial.length());
                SweptSurface::Ruled { surface: Surface::Cone(cone), v0: 0.0, sigma }
            }
            Curve::Circle(circle) => {
                let f = circle.frame();
                let (c, n) = (f.origin(), Vector::from(f.z()));
                let coplanar = n.dot(z).abs() <= ANGULAR && (c - self.axis.origin).dot(n).abs() <= tol;
                if !coplanar {
                    return general();
                }
                // The circle lies in a plane through the axis: a sphere if centered on it,
                // otherwise a torus. `v` is the angle in the meridian plane from the radial
                // direction `e`.
                let (center_foot, center_radial) = self.radial(c);
                let on_axis = center_radial.length() <= tol;
                let mid = circle.point(0.5 * (first + last));
                let e = if on_axis { self.radial(mid).1 } else { center_radial };
                let Some(e) = Direction::try_from_vector(e) else { return general() };
                let (ev, x, y) = (Vector::from(e), Vector::from(f.x()), Vector::from(f.y()));
                let sigma = (x.dot(ev) * y.dot(z) - x.dot(z) * y.dot(ev)).signum();
                let phase = x.dot(z).atan2(x.dot(ev));
                let frame = Frame::new(center_foot, self.axis.direction, e);
                if on_axis {
                    // Latitudes must stay within [-π/2, π/2] (one side of the axis).
                    let latitude = |t: f64| phase + sigma * t;
                    let shift = (latitude(0.5 * (first + last)) / TAU).round() * TAU;
                    let (lo, hi) = {
                        let (a, b) = (latitude(first) - shift, latitude(last) - shift);
                        (a.min(b), a.max(b))
                    };
                    if lo < -FRAC_PI_2 - tolerance::PARAMETRIC || hi > FRAC_PI_2 + tolerance::PARAMETRIC {
                        return general();
                    }
                    let sphere = Sphere::new(frame, circle.radius());
                    SweptSurface::Ruled { surface: Surface::Sphere(sphere), v0: phase - shift, sigma }
                } else if circle.radius() <= center_radial.length() {
                    let torus = Torus::new(frame, center_radial.length(), circle.radius());
                    SweptSurface::Ruled { surface: Surface::Torus(torus), v0: phase, sigma }
                } else {
                    general()
                }
            }
            _ => general(),
        }
    }

    /// The face swept by an edge, with the orientation that makes its normal the sweep's
    /// `∂u × ∂v` (u the angle, v the edge's parameter); `None` for an edge on the axis or
    /// a degenerate edge.
    fn edge_face(&mut self, e: EdgeId) -> Result<Option<(FaceId, Orientation)>> {
        if let Some(&face) = self.faces.get(&e) {
            return Ok(face);
        }
        let Some(slot) = self.model.get(self.model[e].curve).cloned() else {
            self.faces.insert(e, None);
            return Ok(None);
        };
        if self.is_on_axis(e) {
            self.faces.insert(e, None);
            return Ok(None);
        }
        let (curve, a, b) = (slot.curve, slot.first, slot.last);
        let (start, end) = self.edge_vertices(e);
        let (start_path, end_path) = (self.vertex_path(start)?, self.vertex_path(end)?);
        let end_edge = self.end_edge(e)?;
        let swept = self.swept_surface(&curve, a, b);
        let angle = self.angle;
        let line = |origin: (f64, f64), direction: (f64, f64)| {
            Curve2::Line(Line2::new(Point2::new(origin.0, origin.1), Direction2::new(direction.0, direction.1)))
        };
        // The loop, counter-clockwise in (angle, edge parameter): the start vertex's arc, the
        // end copy, the end vertex's arc backwards, the edge backwards.
        let mut uses: Vec<(EdgeId, Orientation, Option<(Curve2, f64, f64)>)> = Vec::new();
        let surface = match &swept {
            SweptSurface::Ruled { surface, v0, sigma } => {
                let v = |t: f64| v0 + sigma * t;
                uses.push((start_path, Orientation::Forward, Some((line((0.0, v(a)), (1.0, 0.0)), 0.0, angle))));
                uses.push((end_edge, Orientation::Forward, Some((line((angle, *v0), (0.0, *sigma)), a, b))));
                uses.push((end_path, Orientation::Reversed, Some((line((0.0, v(b)), (1.0, 0.0)), 0.0, angle))));
                uses.push((e, Orientation::Reversed, Some((line((0.0, *v0), (0.0, *sigma)), a, b))));
                surface.clone()
            }
            SweptSurface::Planar(plane) => {
                // On a plane, pcurves follow from the 3D curves; arcs of vertices on the axis
                // collapse to points and are left out. A full revolution gives an annulus
                // (or a disc) bounded by the two circles alone, as separate wires: like OCCT,
                // no seam on a plane.
                let mut boundary = vec![(start_path, Orientation::Forward)];
                if !self.full {
                    boundary.push((end_edge, Orientation::Forward));
                }
                boundary.push((end_path, Orientation::Reversed));
                if !self.full {
                    boundary.push((e, Orientation::Reversed));
                }
                for (edge, orientation) in boundary {
                    if let Some(slot) = self.model.get(self.model[edge].curve) {
                        let pcurve = project::analytic_pcurve(&slot.curve, plane).expect("edge in the plane");
                        uses.push((edge, orientation, Some((pcurve, slot.first, slot.last))));
                    }
                }
                plane.clone()
            }
        };
        // Orient the loop counter-clockwise about the surface's normal.
        let t = 0.5 * (a + b);
        let (p, dp) = curve.d1(t);
        let sweep_normal = Vector::from(self.axis.direction).cross(p - self.axis.origin).cross(dp);
        let surface_normal = match &swept {
            SweptSurface::Ruled { surface, v0, sigma } => {
                let (_, su, sv) = surface.d1(0.0, v0 + sigma * t);
                su.cross(sv)
            }
            SweptSurface::Planar(Surface::Plane(plane)) => Vector::from(plane.normal()),
            SweptSurface::Planar(_) => unreachable!("planar sweeps are planes"),
        };
        let canonical = if sweep_normal.dot(surface_normal) < 0.0 {
            uses.reverse();
            for u in &mut uses {
                u.1 = u.1.reversed();
            }
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        let tol = self.tolerance();
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tol)?;
        // One loop, or one wire per circle for a full planar sweep.
        let separate = self.full && matches!(swept, SweptSurface::Planar(_));
        let loops: Vec<Vec<_>> = if separate { uses.into_iter().map(|u| vec![u]).collect() } else { vec![uses] };
        for uses in loops {
            let mut coedges = Vec::with_capacity(uses.len());
            for (edge, orientation, pcurve) in uses {
                let coedge = self.model.add_coedge(edge, orientation)?;
                if let Some((curve, first, last)) = pcurve {
                    self.model.set_pcurve_with_range(coedge, curve, first, last)?;
                }
                coedges.push(coedge);
            }
            let wire = self.model.add_wire(&coedges)?;
            self.model.add(face, wire.into())?;
        }
        let id = FaceId::from_index(face.node.index);
        self.faces.insert(e, Some((id, canonical)));
        self.record(e, Shape { orientation: canonical, ..face });
        Ok(Some((id, canonical)))
    }

    /// The solid swept by a face used with `orientation`.
    fn face_solid(&mut self, face: FaceId, orientation: Orientation) -> Result<Shape> {
        let shell = self.model.new_empty(NodeKind::Shell);
        // A partial revolution is closed by the profile at its start and a rotated copy at
        // its end; a full one needs no caps (the profile's edges become seams).
        if !self.full {
            self.model.add(shell, Shape { orientation, ..Shape::from(face) })?;
        }
        let wire_uses = self.model[face].wires.clone();
        let mut end_wires = Vec::new();
        for wire_use in wire_uses {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            let coedges = self.model[wire].coedges.clone();
            let mut end_coedges = Vec::new();
            for coedge in coedges {
                let (edge, coedge_orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                let used = coedge_orientation.compose(wire_orientation).compose(orientation);
                if let Some((side, canonical)) = self.edge_face(edge)? {
                    self.model.add(shell, Shape { orientation: canonical.compose(used), ..Shape::from(side) })?;
                }
                if !self.full {
                    let end = self.end_edge(edge)?;
                    let copy = self.model.add_coedge(end, coedge_orientation)?;
                    if let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) {
                        self.model.set_pcurve_with_range(copy, pcurve, first, last)?;
                    }
                    end_coedges.push(copy);
                }
            }
            end_wires.push((end_coedges, wire_orientation));
        }
        if !self.full {
            // The rotated copy of the face closes the end, facing the other way.
            let tol = self.model[face].tolerance;
            let surface = self.model.face_surface(face).expect("face with a surface").transformed(&self.rotation);
            let end_face = self.model.new_empty(NodeKind::Face);
            self.model.set_face_surface(end_face, surface, tol)?;
            for (coedges, wire_orientation) in end_wires {
                let wire = self.model.add_wire(&coedges)?;
                self.model.add(end_face, Shape { orientation: wire_orientation, ..Shape::from(wire) })?;
            }
            self.model.add(shell, Shape { orientation: orientation.reversed(), ..end_face })?;
            self.ends.push((face.into(), Shape::from(face), end_face));
        }
        // Orient the shell so that the solid has a positive volume.
        let volume = properties::closed_shell_volume(self.model, shell);
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
        self.record(face, solid);
        Ok(solid)
    }
}

/// The state of one extrusion: what each profile vertex and edge has swept, and their
/// translated copies at the end.
struct Extruder<'m> {
    model: &'m mut Model,
    v: Vector,
    translation: Transform,
    vertex_paths: HashMap<VertexId, EdgeId>,
    end_vertices: HashMap<VertexId, VertexId>,
    end_edges: HashMap<EdgeId, EdgeId>,
    /// The face swept by each edge (none for degenerate edges), with the orientation that
    /// makes its normal `v × (edge tangent)`.
    faces: HashMap<EdgeId, Option<(FaceId, Orientation)>>,
    generated: Vec<(NodeId, Vec<Shape>)>,
    ends: Vec<(NodeId, Shape, Shape)>,
    copies: Vec<(NodeId, Shape)>,
}

impl Extruder<'_> {
    fn check<T: crate::topology::Node>(&self, id: crate::topology::Id<T>) -> Result<()> {
        self.model.get(id).map(|_| ()).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))
    }

    fn point(&self, v: VertexId) -> Point {
        self.model[v].point
    }

    fn edge_vertices(&self, e: EdgeId) -> (VertexId, VertexId) {
        let edge = &self.model[e];
        (self.model[edge.start].vertex, self.model[edge.end].vertex)
    }

    fn record(&mut self, input: impl Into<NodeId>, output: Shape) {
        self.generated.push((input.into(), vec![output]));
    }

    /// The translated copy of a vertex.
    fn end_vertex(&mut self, v: VertexId) -> VertexId {
        if let Some(&end) = self.end_vertices.get(&v) {
            return end;
        }
        let tol = self.model[v].tolerance;
        let end = self.model.add_vertex(self.point(v) + self.v, tol);
        self.end_vertices.insert(v, end);
        self.copies.push((v.into(), Shape::from(end)));
        end
    }

    /// The straight edge swept by a vertex, from it to its copy, parametrized by length.
    fn vertex_path(&mut self, v: VertexId) -> Result<EdgeId> {
        if let Some(&e) = self.vertex_paths.get(&v) {
            return Ok(e);
        }
        let end = self.end_vertex(v);
        // (At least the vertex's tolerance, as OCCT's sweeps give the edges they generate.)
        let tol = self.model[v].tolerance.max(tolerance::LINEAR);
        let edge = self.model.add_edge(v, end, tol)?;
        let line = Line::new(self.point(v), Direction::from_vector(self.v));
        self.model.set_edge_curve(edge.into(), Curve::Line(line), tol)?;
        self.model.set_edge_range(edge.into(), 0.0, self.v.length())?;
        self.vertex_paths.insert(v, edge);
        self.record(v, edge.into());
        Ok(edge)
    }

    /// The translated copy of an edge.
    fn end_edge(&mut self, e: EdgeId) -> Result<EdgeId> {
        if let Some(&end) = self.end_edges.get(&e) {
            return Ok(end);
        }
        let (a, b) = self.edge_vertices(e);
        let (a, b) = (self.end_vertex(a), self.end_vertex(b));
        let tol = self.model[e].tolerance;
        let end = self.model.add_edge(a, b, tol)?;
        if let Some(slot) = self.model.get(self.model[e].curve) {
            let (curve, first, last) = (slot.curve.transformed(&self.translation), slot.first, slot.last);
            self.model.set_edge_curve(end.into(), curve, tol)?;
            self.model.set_edge_range(end.into(), first, last)?;
        }
        self.end_edges.insert(e, end);
        self.copies.push((e.into(), Shape::from(end)));
        Ok(end)
    }

    /// The face swept by an edge, with the orientation that makes its normal
    /// `v × (edge tangent)`; `None` for a degenerate edge. A straight edge sweeps a plane
    /// (u along the edge, v across), a circle swept along its axis a cylinder (u the
    /// circle's angle, v the height).
    fn edge_face(&mut self, e: EdgeId) -> Result<Option<(FaceId, Orientation)>> {
        if let Some(&face) = self.faces.get(&e) {
            return Ok(face);
        }
        let Some(slot) = self.model.get(self.model[e].curve).cloned() else {
            self.faces.insert(e, None);
            return Ok(None);
        };
        let (curve, a, b) = (slot.curve, slot.first, slot.last);
        let (start, end) = self.edge_vertices(e);
        let (start_path, end_path) = (self.vertex_path(start)?, self.vertex_path(end)?);
        let end_edge = self.end_edge(e)?;
        let len = self.v.length();
        let line = |origin: (f64, f64), direction: (f64, f64)| {
            Curve2::Line(Line2::new(Point2::new(origin.0, origin.1), Direction2::new(direction.0, direction.1)))
        };
        // The surface and, in its (u, v), the origin of the edge's pcurve `(t + u0, 0)` and
        // the offset `(du, dv)` of the end copy; the paths run from `(t + u0, 0)` to there.
        let (surface, u0, du, dv) = match &curve {
            Curve::Line(l) => {
                let x = Vector::from(l.direction);
                let across = self.v - x * self.v.dot(x);
                let y = Direction::try_from_vector(across)
                    .filter(|_| across.length() > tolerance::LINEAR)
                    .ok_or_else(|| Error::InvalidInput("a straight edge extruded along itself".into()))?;
                let plane = Plane::from_frame(Frame::new(l.origin, Direction::from_vector(x.cross(across)), l.direction));
                (Surface::Plane(plane), 0.0, self.v.dot(x), self.v.dot(Vector::from(y)))
            }
            Curve::Circle(c) if Vector::from(c.frame().z()).cross(self.v).length() <= ANGULAR * len.max(1.0) => {
                let height = self.v.dot(Vector::from(c.frame().z()));
                (Surface::Cylinder(Cylinder::new(*c.frame(), c.radius())), 0.0, 0.0, height)
            }
            // Any other curve sweeps a surface of linear extrusion: u the curve's parameter,
            // v the distance along the vector.
            other => {
                let extrusion = ExtrusionSurface::new(other.clone(), Direction::from_vector(self.v));
                (Surface::Extrusion(extrusion), 0.0, 0.0, len)
            }
        };
        let path = |t: f64| (line((t + u0, 0.0), (du / len, dv / len)), 0.0, len);
        let mut uses: Vec<(EdgeId, Orientation, (Curve2, f64, f64))> = vec![
            (e, Orientation::Forward, (line((u0, 0.0), (1.0, 0.0)), a, b)),
            (end_path, Orientation::Forward, path(b)),
            (end_edge, Orientation::Reversed, (line((u0 + du, dv), (1.0, 0.0)), a, b)),
            (start_path, Orientation::Reversed, path(a)),
        ];
        // Counter-clockwise in (u, v) when the copy is on the +v side.
        if dv < 0.0 {
            uses.reverse();
            for u in &mut uses {
                u.1 = u.1.reversed();
            }
        }
        let (_, dp) = curve.d1(0.5 * (a + b));
        let (_, su, sv) = surface.d1(0.5 * (a + b) + u0, 0.0);
        let canonical = if self.v.cross(dp).dot(su.cross(sv)) < 0.0 { Orientation::Reversed } else { Orientation::Forward };
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, orientation, (pcurve, first, last)) in uses {
            let coedge = self.model.add_coedge(edge, orientation)?;
            self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = self.model.add_wire(&coedges)?;
        self.model.add(face, wire.into())?;
        let id = FaceId::from_index(face.node.index);
        self.faces.insert(e, Some((id, canonical)));
        self.record(e, Shape { orientation: canonical, ..face });
        Ok(Some((id, canonical)))
    }

    /// The solid swept by a face used with `orientation`: the face, the faces its edges
    /// sweep, and its translated copy.
    fn face_solid(&mut self, face: FaceId, orientation: Orientation) -> Result<Shape> {
        let shell = self.model.new_empty(NodeKind::Shell);
        self.model.add(shell, Shape { orientation, ..Shape::from(face) })?;
        let wire_uses = self.model[face].wires.clone();
        let mut end_wires = Vec::new();
        for wire_use in wire_uses {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            let coedges = self.model[wire].coedges.clone();
            let mut end_coedges = Vec::new();
            for coedge in coedges {
                let (edge, coedge_orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                let used = coedge_orientation.compose(wire_orientation).compose(orientation);
                if let Some((side, canonical)) = self.edge_face(edge)? {
                    self.model.add(shell, Shape { orientation: canonical.compose(used), ..Shape::from(side) })?;
                }
                let end = self.end_edge(edge)?;
                let copy = self.model.add_coedge(end, coedge_orientation)?;
                if let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) {
                    self.model.set_pcurve_with_range(copy, pcurve, first, last)?;
                }
                end_coedges.push(copy);
            }
            end_wires.push((end_coedges, wire_orientation));
        }
        // The translated copy of the face closes the end, facing the other way.
        let tol = self.model[face].tolerance;
        let surface = self.model.face_surface(face).expect("face with a surface").transformed(&self.translation);
        let end_face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(end_face, surface, tol)?;
        for (coedges, wire_orientation) in end_wires {
            let wire = self.model.add_wire(&coedges)?;
            self.model.add(end_face, Shape { orientation: wire_orientation, ..Shape::from(wire) })?;
        }
        self.model.add(shell, Shape { orientation: orientation.reversed(), ..end_face })?;
        // Orient the shell so that the solid has a positive volume.
        let volume = properties::closed_shell_volume(self.model, shell);
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
        self.record(face, solid);
        self.ends.push((face.into(), Shape::from(face), end_face));
        Ok(solid)
    }
}

/// A curve's kind, for messages.
fn curve_kind(curve: &Curve) -> String {
    format!("{curve:?}").split('(').next().unwrap_or_default().to_lowercase()
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    use super::*;
    use crate::make;
    use crate::validity;

    /// A planar face through the points (x, 0, z), in order.
    fn polygon(model: &mut Model, points: &[(f64, f64)]) -> Shape {
        let vertices: Vec<Shape> = points.iter().map(|&(x, z)| make::vertex(model, Point::new(x, 0.0, z))).collect();
        let edges: Vec<Shape> = (0..vertices.len())
            .map(|i| make::edge(model, vertices[i], vertices[(i + 1) % vertices.len()]).unwrap())
            .collect();
        let wire = make::wire(model, &edges).unwrap();
        let p: Vec<Point> = points.iter().map(|&(x, z)| Point::new(x, 0.0, z)).collect();
        let normal = (p[1] - p[0]).cross(p[2] - p[1]);
        face_on_plane(model, wire, Plane::new(p[0], Direction::from_vector(normal)))
    }

    fn face_on_plane(model: &mut Model, wire: Shape, plane: Plane) -> Shape {
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, Surface::Plane(plane), tolerance::LINEAR).unwrap();
        model.add(face, wire).unwrap();
        face
    }

    #[test]
    fn a_square_piped_along_a_planar_s_curve_keeps_its_area_across() {
        use crate::geometry::curve::BezierCurve;
        // The square across the spine at its start (the spine starts along +y), its centre
        // on it: the volume is the square's area times the spine's length (Pappus).
        let mut model = Model::new();
        let half = 0.5;
        let vertices: Vec<Shape> = [(-half, -half), (half, -half), (half, half), (-half, half)]
            .iter()
            .map(|&(x, z)| make::vertex(&mut model, Point::new(x, 0.0, z)))
            .collect();
        let edges: Vec<Shape> = (0..4).map(|i| make::edge(&mut model, vertices[i], vertices[(i + 1) % 4]).unwrap()).collect();
        let wire = make::wire(&mut model, &edges).unwrap();
        let square = face_on_plane(&mut model, wire, Plane::new(Point::ORIGIN, Direction::new(0.0, -1.0, 0.0)));
        let s_curve = Curve::Bezier(BezierCurve::new(vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
        ]));
        let spine = make::edge_from_curve(&mut model, s_curve).unwrap();
        let solid = pipe(&mut model, square, spine).unwrap().shape;
        let report = validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        let length = properties::length(&model, spine).mass;
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - length).abs() < 1e-4 * length, "{volume} vs {length}");
    }

    fn surfaces(model: &Model, shape: Shape) -> Vec<&'static str> {
        let mut kinds: Vec<&'static str> = model
            .explore(shape, NodeKind::Face)
            .map(|f| match model.face_surface(FaceId::from_index(f.node.index)).unwrap() {
                Surface::Plane(_) => "plane",
                Surface::Cylinder(_) => "cylinder",
                Surface::Cone(_) => "cone",
                Surface::Sphere(_) => "sphere",
                Surface::Torus(_) => "torus",
                Surface::Revolution(_) => "revolution",
                _ => "other",
            })
            .collect();
        kinds.sort();
        kinds
    }

    #[track_caller]
    fn assert_valid_solid(model: &Model, shape: Shape, volume: f64) {
        let report = validity::check(model, shape);
        assert!(report.problems.is_empty(), "{report:?}");
        let props = properties::volume(model, shape);
        assert!((props.mass - volume).abs() <= 1e-9 * volume.max(1.0), "{} != {volume}", props.mass);
    }

    #[test]
    fn sections_turned_or_started_apart_are_matched() {
        // Equal circles across a straight spine, the middle one started a quarter turn on
        // and the last one running the other way round: still a cylinder.
        let mut model = Model::new();
        let spine = make::edge_between(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0)).unwrap();
        let mut sections = Vec::new();
        for (x, normal, start) in [(0.0, Direction::X, Direction::Y), (5.0, Direction::X, Direction::Z), (10.0, -Direction::X, Direction::Y)] {
            let circle = Curve::Circle(Circle::new(Frame::new(Point::new(x, 0.0, 0.0), normal, start), 1.0));
            let edge = make::edge_from_curve(&mut model, circle).unwrap();
            sections.push((make::wire(&mut model, &[edge]).unwrap(), PipeOptions::default()));
        }
        let solid = pipe_through(&mut model, &sections, spine, true).unwrap();
        let report = validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - 10.0 * PI).abs() < 1e-6 * volume, "{volume}");
    }

    /// Circles of radius 1.5, 1.5, 1 and 1 at the ends of a spine of a line (5 long), a half
    /// circle (radius 4) and a line (CadQuery's multi-section sweep example): the radius is
    /// the cubic through them by how far along the spine each stands (as OCCT's, by their
    /// abscissa), not by their places among its edges; the volume π ∫ r² along it.
    #[test]
    fn sections_blend_by_how_far_along_the_spine_they_stand() {
        let mut model = Model::new();
        let first = make::edge_between(&mut model, Point::new(-5.0, 0.0, 4.0), Point::new(0.0, 0.0, 4.0)).unwrap();
        let half = Curve::Circle(Circle::new(Frame::new(Point::ORIGIN, Direction::Y, Direction::Z), 4.0));
        let half = make::edge_from_curve_range(&mut model, half, 0.0, PI).unwrap();
        let last = make::edge_between(&mut model, Point::new(0.0, 0.0, -4.0), Point::new(-5.0, 0.0, -4.0)).unwrap();
        let spine = make::wire(&mut model, &[first, half, last]).unwrap();
        let mut sections = Vec::new();
        for (x, z, r) in [(-5.0, 4.0, 1.5), (0.0, 4.0, 1.5), (0.0, -4.0, 1.0), (-5.0, -4.0, 1.0)] {
            let circle = Curve::Circle(Circle::new(Frame::new(Point::new(x, 0.0, z), Direction::X, Direction::Y), r));
            let edge = make::edge_from_curve(&mut model, circle).unwrap();
            sections.push((make::wire(&mut model, &[edge]).unwrap(), PipeOptions::default()));
        }
        let solid = pipe_through(&mut model, &sections, spine, true).unwrap();
        let report = validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        let at = [0.0, 5.0, 5.0 + 4.0 * PI, 10.0 + 4.0 * PI];
        let radius = |s: f64| {
            (0..4)
                .map(|i| [1.5, 1.5, 1.0, 1.0][i] * (0..4).filter(|&j| j != i).map(|j| (s - at[j]) / (at[i] - at[j])).product::<f64>())
                .sum::<f64>()
        };
        // (Simpson's rule, far finer than the cubic needs.)
        let n = 2000;
        let h = at[3] / n as f64;
        let expected = PI * h / 3.0 * (0..=n).map(|i| radius(i as f64 * h).powi(2) * if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 }).sum::<f64>();
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - expected).abs() < 1e-5 * expected, "{volume} vs {expected}");
    }

    /// The face of the polygon through `points` (in the XY plane).
    fn polygon_xy(model: &mut Model, points: &[(f64, f64)]) -> Shape {
        let p: Vec<Point> = points.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
        let edges: Vec<Shape> = (0..p.len()).map(|i| make::edge_between(model, p[i], p[(i + 1) % p.len()]).unwrap()).collect();
        let wire = make::wire(model, &edges).unwrap();
        make::planar_face(model, wire).unwrap()
    }

    /// A square along a straight spine, its frames pointing at a guide that turns a quarter
    /// round the spine: the square turns with it (a twisted prism, of the square's volume
    /// all the same), its corner at (1, 0) ending at (0, 1).
    #[test]
    fn a_guide_turns_the_profile() {
        let mut model = Model::new();
        let square = polygon_xy(&mut model, &[(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]);
        let spine = make::edge_between(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 10.0)).unwrap();
        // (A quarter turn of a helix about the spine, as a B-spline.)
        let points: Vec<Point> = (0..=32).map(|i| { let a = FRAC_PI_2 * i as f64 / 32.0; Point::new(2.0 * a.cos(), 2.0 * a.sin(), 10.0 * i as f64 / 32.0) }).collect();
        let helix = crate::geometry::fit::interpolate(&points, false, 1e-9).unwrap();
        let guide = make::edge_from_curve(&mut model, Curve::BSpline(helix)).unwrap();
        let options = PipeOptions { guide: Some(Guide::Wire { wire: guide, curvilinear: true }), ..Default::default() };
        let swept = pipe_with(&mut model, square, spine, &options).unwrap();
        assert!(validity::check(&model, swept.shape).is_valid());
        let volume = properties::volume(&model, swept.shape).mass;
        assert!((volume - 20.0).abs() < 1e-3 * 20.0, "{volume}");
        let &(_, _, last) = swept.ends().first().unwrap();
        let corner = model.explore(last, NodeKind::Vertex).map(|v| model.point(v)).min_by(|a, b| a.distance(Point::new(0.0, 1.0, 10.0)).total_cmp(&b.distance(Point::new(0.0, 1.0, 10.0)))).unwrap();
        assert!(corner.distance(Point::new(0.0, 1.0, 10.0)) < 1e-4, "{corner:?}");
    }

    /// A fixed frame carries a rectangle by translations alone, corrected square to the
    /// frame's direction first: along an arc, the rectangle turned to face X sweeps its area
    /// times how far the arc goes along X.
    #[test]
    fn a_fixed_frame_translates_the_profile() {
        let mut model = Model::new();
        let rect = polygon_xy(&mut model, &[(-0.5, -0.25), (0.5, -0.25), (0.5, 0.25), (-0.5, 0.25)]);
        let arc = Curve::Circle(Circle::new(Frame::new(Point::new(2.0, 0.0, 0.0), -Direction::Y, Direction::X), 2.0));
        let spine = make::edge_from_curve_range(&mut model, arc, PI / 2.0, PI).unwrap();
        let spine = make::wire(&mut model, &[spine]).unwrap();
        let options = PipeOptions { correction: true, guide: Some(Guide::Fixed(Direction::X)), ..Default::default() };
        let swept = pipe_with(&mut model, rect, spine, &options).unwrap();
        assert!(validity::check(&model, swept.shape).is_valid());
        let (a, b) = {
            let (a, b) = model.edge_vertices(model.explore(spine, NodeKind::Edge).next().unwrap());
            (model.point(a), model.point(b))
        };
        let volume = properties::volume(&model, swept.shape).mass;
        let expected = 0.5 * (b.x - a.x).abs();
        assert!((volume - expected).abs() < 1e-4 * expected, "{volume} {expected}");
        for (_, first, last) in swept.ends() {
            for f in [first, last] {
                let n = match model.surface(*f).unwrap().transformed(&f.placement) { Surface::Plane(p) => p.normal(), other => panic!("{other:?}") };
                assert!(n.dot(Direction::X).abs() > 1.0 - 1e-9, "{n:?}");
            }
        }
        // (The sides swept by the edges along Z lie in the planes y = ±0.25: planes.)
        assert_eq!(surfaces(&model, swept.shape).iter().filter(|s| **s == "plane").count(), 4);
    }

    /// A profile off a straight spine, its plane meeting the spine's end: it stands there,
    /// and sweeps back from it (the spine's start carrying it below its plane).
    #[test]
    fn a_profile_standing_at_the_spine_end_sweeps_back() {
        let mut model = Model::new();
        let square = polygon_xy(&mut model, &[(4.0, 4.0), (5.0, 4.0), (5.0, 5.0), (4.0, 5.0)]);
        let spine = make::edge_between(&mut model, Point::new(0.0, -3.0, -3.0), Point::ORIGIN).unwrap();
        let swept = pipe(&mut model, square, spine).unwrap();
        let b = properties::bounding_box(&model, swept.shape);
        assert!((b.min.z + 3.0).abs() < 1e-6 && b.max.z.abs() < 1e-6, "{b:?}");
        assert!((b.min.y - 1.0).abs() < 1e-6 && (b.max.y - 5.0).abs() < 1e-6, "{b:?}");
    }

    /// A square swept through to a larger one: each edge of either sweeps one side, the
    /// first section's caps its first and last faces.
    #[test]
    fn sweeps_through_sections_record_their_sides_and_caps() {
        let mut model = Model::new();
        let spine = make::edge_between(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 1.0)).unwrap();
        let square = |model: &mut Model, h: f64, z: f64| {
            let p = [(-h, -h), (h, -h), (h, h), (-h, h)].map(|(x, y)| Point::new(x, y, z));
            let edges: Vec<Shape> = (0..4).map(|i| make::edge_between(model, p[i], p[(i + 1) % 4]).unwrap()).collect();
            make::wire(model, &edges).unwrap()
        };
        let (small, large) = (square(&mut model, 0.5, 0.0), square(&mut model, 1.0, 1.0));
        let swept = pipe_through_with(&mut model, &[(small, PipeOptions::default()), (large, PipeOptions::default())], spine, true, None).unwrap();
        let faces = model.unique_subshapes(swept.shape, NodeKind::Face);
        assert_eq!(faces.len(), 6);
        for (a, b) in model.unique_subshapes(small, NodeKind::Edge).into_iter().zip(model.unique_subshapes(large, NodeKind::Edge)) {
            let (sa, sb) = (swept.generated(a), swept.generated(b));
            assert_eq!(sa.len(), 1);
            assert_eq!(sa, sb);
            assert!(faces.iter().any(|f| f.node == sa[0].node));
        }
        let &(source, first, last) = swept.ends().first().unwrap();
        assert_eq!(source, small.node);
        assert!((properties::area(&model, first).mass - 1.0).abs() < 1e-9);
        assert!((properties::area(&model, last).mass - 4.0).abs() < 1e-9);
    }

    #[test]
    fn extruded_polygon_is_a_prism() {
        // An L-shaped profile (area 3) in the XZ plane, extruded straight and slanted:
        // volume = area x height.
        for v in [Vector::new(0.0, 2.0, 0.0), Vector::new(1.0, 2.0, 0.5), Vector::new(0.0, -2.0, 0.0)] {
            let mut model = Model::new();
            let face = polygon(&mut model, &[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]);
            let solid = extrude(&mut model, face, v).unwrap().shape;
            assert_valid_solid(&model, solid, 6.0);
            assert_eq!(surfaces(&model, solid), vec!["plane"; 8]);
        }
    }

    /// An extruded face's history: it closes the start itself, its copy moved along the
    /// extrusion closes the end, and each of its edges sweeps one of the side faces.
    #[test]
    fn extrusions_record_their_ends_and_sides() {
        let mut model = Model::new();
        let face = polygon(&mut model, &[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]);
        let v = Vector::new(0.0, 3.0, 0.0);
        let swept = extrude(&mut model, face, v).unwrap();
        let &(source, first, last) = swept.ends().first().unwrap();
        assert_eq!(source, face.node);
        assert_eq!(first.node, face.node);
        let (a, b) = (properties::area(&model, first).center_of_mass, properties::area(&model, last).center_of_mass);
        assert!((b - a - v).length() < 1e-9, "{a:?} {b:?}");
        let faces = model.unique_subshapes(swept.shape, NodeKind::Face);
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            let sides = swept.generated(edge);
            assert_eq!(sides.len(), 1);
            assert!(faces.iter().any(|f| f.node == sides[0].node));
            // (Its copy at the end: the edge moved along the extrusion.)
            let &(_, copy) = swept.copies().iter().find(|(n, _)| *n == edge.node).unwrap();
            let (a, b) = (properties::length(&model, edge).center_of_mass, properties::length(&model, copy).center_of_mass);
            assert!((b - a - v).length() < 1e-9);
        }
    }

    #[test]
    fn placed_profiles() {
        // A profile moved by a placement sweeps the same solid as the profile built there.
        let t = Transform::rotation(Axis::new(Point::new(1.0, 2.0, 0.0), Direction::Z), 0.7)
            * Transform::translation(Vector::new(3.0, 0.0, 1.0));
        let mut model = Model::new();
        let face = polygon(&mut model, &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)]);
        let placed = face.moved(&t);
        // The unit square (in the XZ plane, 1 to 2 from the Z axis) swept along its normal
        // makes a 1 x 1 x 2 box; revolved about the Z axis, a ring of volume 2π·1.5.
        let prism = extrude(&mut model, placed, Vector::new(0.0, 2.0, 0.0).transformed(&t)).unwrap().shape;
        assert_valid_solid(&model, prism, 2.0);
        let ring = revolve(&mut model, placed, Axis::Z.transformed(&t), TAU).unwrap().shape;
        assert_valid_solid(&model, ring, 3.0 * PI);
        assert_eq!(ring.placement, t);
    }

    #[test]
    fn extruded_disc_is_a_cylinder() {
        let mut model = Model::new();
        let circle = make::edge_from_curve(
            &mut model,
            Curve::Circle(Circle::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X), 2.0)),
        )
        .unwrap();
        let wire = make::wire(&mut model, &[circle]).unwrap();
        let disc = make::planar_face(&mut model, wire).unwrap();
        let solid = extrude(&mut model, disc, Vector::new(0.0, 0.0, 3.0)).unwrap().shape;
        assert_valid_solid(&model, solid, 12.0 * PI);
        assert_eq!(surfaces(&model, solid), vec!["cylinder", "plane", "plane"]);
        let area = properties::area(&model, solid).mass;
        assert!((area - (8.0 * PI + 12.0 * PI)).abs() < 1e-9, "{area}");
    }

    #[test]
    fn extruded_vertex_and_wire() {
        let mut model = Model::new();
        let p = make::vertex(&mut model, Point::new(1.0, 2.0, 3.0));
        let edge = extrude(&mut model, p, Vector::new(0.0, 0.0, 4.0)).unwrap().shape;
        assert_eq!(model.edge_range(edge), (0.0, 4.0));
        let (a, b) = model.edge_vertices(edge);
        assert_eq!(model.point(b) - model.point(a), Vector::new(0.0, 0.0, 4.0));
        let face = polygon(&mut model, &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let wire = model.children(face).next().unwrap();
        let shell = extrude(&mut model, wire, Vector::new(0.0, 3.0, 0.0)).unwrap().shape;
        // The square lies in the XZ plane: swept along Y its four sides make a tube.
        assert!((properties::area(&model, shell).mass - 12.0).abs() < 1e-9);
        assert!(extrude(&mut model, face, Vector::ZERO).is_err());
    }

    #[test]
    fn rectangle_sweeps_a_ring() {
        for angle in [TAU, FRAC_PI_2, 4.0] {
            let mut model = Model::new();
            let face = polygon(&mut model, &[(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)]);
            let solid = revolve(&mut model, face, Axis::Z, angle).unwrap().shape;
            assert_valid_solid(&model, solid, 5.0 * PI * angle / TAU);
            let mut expected = vec!["cylinder", "cylinder", "plane", "plane"];
            if angle < TAU {
                expected.extend(["plane", "plane"]);
            }
            expected.sort();
            assert_eq!(surfaces(&model, solid), expected, "angle {angle}");
        }
    }

    #[test]
    fn triangle_on_the_axis_sweeps_a_cone() {
        for angle in [TAU, 1.0] {
            let mut model = Model::new();
            let face = polygon(&mut model, &[(0.0, 0.0), (2.0, 0.0), (0.0, 3.0)]);
            let solid = revolve(&mut model, face, Axis::Z, angle).unwrap().shape;
            assert_valid_solid(&model, solid, 4.0 * PI * angle / TAU);
            assert!(surfaces(&model, solid).contains(&"cone"));
        }
    }

    #[test]
    fn faces_meeting_on_the_axis_get_their_own_degenerate_edges() {
        // Two cones meet at the origin: each has its own degenerate edge there (OCCT's
        // bugs/modalg_1/bug15036 counts them).
        let mut model = Model::new();
        let face = polygon(&mut model, &[(0.0, 0.0), (2.0, 1.0), (1.0, 2.0)]);
        let solid = revolve(&mut model, face, Axis::Z, TAU).unwrap().shape;
        assert_valid_solid(&model, solid, 3.0 * PI);
        let degenerate: Vec<_> = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .filter(|&e| model.is_degenerate(e))
            .collect();
        assert_eq!(degenerate.len(), 2);
        let faces = model.unique_subshapes(solid, NodeKind::Face);
        for &e in &degenerate {
            let users = faces.iter().filter(|&&f| model.unique_subshapes(f, NodeKind::Edge).contains(&e));
            assert_eq!(users.count(), 1);
        }
    }

    #[test]
    fn half_disc_sweeps_a_sphere() {
        let mut model = Model::new();
        let r = 2.0;
        let arc = make::edge_from_curve_range(
            &mut model,
            Curve::Circle(Circle::new(Frame::new(Point::ORIGIN, Direction::NEG_Y, Direction::X), r)),
            -FRAC_PI_2,
            FRAC_PI_2,
        )
        .unwrap();
        let (bottom, top) = {
            let e = &model[EdgeId::from_index(arc.node.index)];
            (Shape::from(model[e.start].vertex), Shape::from(model[e.end].vertex))
        };
        let axis_edge = make::edge(&mut model, top, bottom).unwrap();
        let wire = make::wire(&mut model, &[arc, axis_edge]).unwrap();
        let face = face_on_plane(&mut model, wire, Plane::new(Point::ORIGIN, Direction::NEG_Y));
        let solid = revolve(&mut model, face, Axis::Z, TAU).unwrap().shape;
        assert_valid_solid(&model, solid, 4.0 / 3.0 * PI * r * r * r);
        assert_eq!(surfaces(&model, solid), ["sphere"], "a full revolution has no caps");
    }

    #[test]
    fn circle_sweeps_a_torus() {
        for angle in [TAU, FRAC_PI_2] {
            let mut model = Model::new();
            let frame = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Y, Direction::X);
            let circle = make::edge_from_curve(&mut model, Curve::Circle(Circle::new(frame, 1.0))).unwrap();
            let wire = make::wire(&mut model, &[circle]).unwrap();
            let face = face_on_plane(&mut model, wire, Plane::from_frame(frame));
            let solid = revolve(&mut model, face, Axis::Z, angle).unwrap().shape;
            assert_valid_solid(&model, solid, 10.0 * PI * PI * angle / TAU);
            assert!(surfaces(&model, solid).contains(&"torus"));
        }
    }

    #[test]
    fn negative_angle_revolves_the_other_way() {
        let mut model = Model::new();
        let face = polygon(&mut model, &[(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)]);
        let solid = revolve(&mut model, face, Axis::Z, -1.0).unwrap().shape;
        assert_valid_solid(&model, solid, 5.0 * PI / TAU);
        let center = properties::volume(&model, solid).center_of_mass;
        assert!(center.y < 0.0, "{center:?}");
    }

    #[test]
    fn vertex_and_edge_profiles() {
        let mut model = Model::new();
        let v = make::vertex(&mut model, Point::new(3.0, 0.0, 1.0));
        let sweep = revolve(&mut model, v, Axis::Z, FRAC_PI_2).unwrap();
        let (curve, first, last) = model.edge_curve(sweep.shape).unwrap();
        assert!(curve.point(last).distance(Point::new(0.0, 3.0, 1.0)) <= 1e-12);
        assert!(curve.point(first).distance(Point::new(3.0, 0.0, 1.0)) <= 1e-12);
        assert_eq!(sweep.generated(v), [sweep.shape]);

        // A slanted segment sweeps a conical band: Pappus' area 2π·d·L.
        let e = make::edge_between(&mut model, Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 1.0)).unwrap();
        let face = revolve(&mut model, e, Axis::Z, TAU).unwrap().shape;
        assert!(validity::check(&model, face).problems.is_empty());
        assert_eq!(surfaces(&model, face), ["cone"]);
        let area = properties::area(&model, face).mass;
        assert!((area - TAU * 1.5 * 2f64.sqrt()).abs() <= 1e-9, "{area}");
    }

    #[test]
    fn invalid_angles_are_rejected() {
        let mut model = Model::new();
        let v = make::vertex(&mut model, Point::new(1.0, 0.0, 0.0));
        assert!(revolve(&mut model, v, Axis::Z, 0.0).is_err());
        assert!(revolve(&mut model, v, Axis::Z, 7.0).is_err());
        assert!(revolve(&mut model, v, Axis::Z, f64::NAN).is_err());
    }
}
