//! Changing the geometry of whole shapes.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::convert::ArcParameters;
use crate::geometry::curve2::{BSplineCurve2, Curve2, Line2};
use crate::geometry::surface::Surface;
use crate::geometry::{project, ParametricCurve, ParametricSurface};
use crate::math::{AffineTransform, AffineTransformable, Direction2, Point2, Transform, Transform2, Transformable, Transformable2, Vector2};
use crate::topology::{CoEdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId, EdgeId};

/// Options for [`transform`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TransformOptions {
    /// Transform and keep the meshes of faces that have exact geometry too. When false
    /// they are dropped (they can be recomputed); meshes of mesh-only faces are always
    /// transformed.
    pub copy_meshes: bool,
}

/// A transformed copy of `shape` (geometry rewritten, not just re-placed: see
/// [`Shape::moved`] for that). The result has an identity placement.
///
/// [`Shape::moved`]: crate::topology::Shape::moved
#[doc(alias = "BRepBuilderAPI_Transform")]
pub fn transform(model: &mut Model, shape: Shape, t: &Transform, options: TransformOptions) -> Result<Shape> {
    // A reflection turns the surfaces' normals (∂u × ∂v) inwards: faces are used the other
    // way round, so that solids keep their material inside (as OCCT does).
    let mut copier = Copier { reverse_faces: t.is_negative(), meshes: options.copy_meshes, ..Copier::default() };
    let copy = copier.copy(model, shape.moved(t))?;
    let reversed = copier.reverse_faces && shape.kind() == NodeKind::Face;
    let orientation = if reversed { shape.orientation.reversed() } else { shape.orientation };
    Ok(Shape { orientation, ..copy })
}

/// Transforms the geometry of `shape` and everything below it in place (points, curves,
/// surfaces; meshes too if `options.copy_meshes`, otherwise they are dropped). Nodes shared
/// with other shapes move with it. Pcurves are unchanged: they live in surface parameters.
pub fn transform_in_place(model: &mut Model, shape: Shape, t: &Transform, options: TransformOptions) -> Result<()> {
    let nodes: Vec<NodeId> = std::iter::once(shape).chain(model.explore_all(shape)).map(|s| s.node).collect();
    transform_nodes(model, &nodes, t, options)
}

/// Transforms the geometry of `nodes` (each once) in place by `t`, and the placements of
/// the uses within them so that placed parts move as a whole.
fn transform_nodes(model: &mut Model, nodes: &[NodeId], t: &Transform, options: TransformOptions) -> Result<()> {
    use crate::topology::{CoEdge, FaceId, Product, ProductChild};
    let mut seen = std::collections::HashSet::new();
    let nodes: Vec<NodeId> = nodes.iter().copied().filter(|n| seen.insert(*n)).collect();
    let s = t.scale_factor();
    let tolerance_scale = s.abs().max(1.0);
    let affine = AffineTransform::from(*t);
    let inverse = t.inverse();
    let conjugate = |p: Transform| *t * p * inverse;
    for &node in &nodes {
        match node.kind {
            NodeKind::Vertex => {
                let v = VertexId::from_index(node.index);
                let Some(record) = model.get(v) else { continue };
                let (point, tolerance) = (record.point.transformed(t), record.tolerance * tolerance_scale);
                model.set_vertex_point(v.into(), point, tolerance)?;
            }
            NodeKind::Edge => {
                let e = EdgeId::from_index(node.index);
                let Some(record) = model.get(e) else { continue };
                let tolerance = record.tolerance * tolerance_scale;
                if let Some((curve, first, last)) = model.edge_curve(e) {
                    let scale = parameter_scale(&curve, s);
                    model.set_edge_curve(e.into(), curve.transformed(t), tolerance)?;
                    model.set_edge_range(e.into(), first * scale, last * scale)?;
                    // The coedges' ranges follow the edge's parameters.
                    let coedges: Vec<CoEdgeId> = model.table::<CoEdge>().iter().filter(|(_, c)| c.edge == e).map(|(c, _)| c).collect();
                    for c in coedges {
                        if model.pcurve(c).is_some() {
                            let (f, l) = model.coedge_range(c);
                            model.set_coedge_range(c, f * scale, l * scale)?;
                        }
                    }
                }
                match options.copy_meshes {
                    true => {
                        if let Some(p) = model.edge_polyline(e).cloned() {
                            model.set_edge_polyline(e.into(), Some(p.affine_transformed(&affine)))?;
                        }
                    }
                    false => model.set_edge_polyline(e.into(), None)?,
                }
            }
            NodeKind::Face => {
                let f = FaceId::from_index(node.index);
                let Some(record) = model.get(f) else { continue };
                let tolerance = record.tolerance * tolerance_scale;
                let has_surface = model.face_surface(f).is_some();
                if let Some(surface) = model.face_surface(f).cloned() {
                    let scales = parameter_scales(&surface, s);
                    model.set_face_surface(f.into(), surface.transformed(t), tolerance)?;
                    let coedges: Vec<CoEdgeId> = model.descendants_of::<CoEdge>(f);
                    for c in coedges {
                        if let Some(pcurve) = model.pcurve(c).cloned() {
                            let scaled = scaled_pcurve(&pcurve, scales).ok_or_else(|| Error::InvalidInput("a pcurve that cannot follow the scaling".into()))?;
                            model.set_pcurve(c, scaled)?;
                        }
                        if let Some(p) = model.pcurve_polyline(c).cloned() {
                            let mapped = p.mapped(|q| Point2::new(q.x * scales.0, q.y * scales.1));
                            model.set_pcurve_polyline(c, Some(mapped))?;
                        }
                        if !options.copy_meshes && has_surface {
                            model.set_coedge_mesh_polyline(c, None)?;
                        }
                    }
                }
                // Meshes of faces with exact geometry follow it if asked; a mesh-only face
                // always keeps (and moves) its mesh.
                match (options.copy_meshes || !has_surface, model.face_mesh(f).cloned()) {
                    (true, Some(mesh)) => model.set_face_mesh(f.into(), Some(mesh.affine_transformed(&affine)))?,
                    (false, Some(_)) => model.set_face_mesh(f.into(), None)?,
                    _ => {}
                }
            }
            NodeKind::Compound => {
                let c = crate::topology::CompoundId::from_index(node.index);
                let Some(record) = model.get(c) else { continue };
                for u in record.children.clone() {
                    let p = model[u].placement;
                    model.move_compound_child(u, &(conjugate(p) * p.inverse()))?;
                }
            }
            NodeKind::Product => {
                let id = crate::topology::ProductId::from_index(node.index);
                let Some(record) = model.get::<Product>(id).cloned() else { continue };
                for u in record.children {
                    if let Some(child) = model.get::<ProductChild>(u) {
                        let instance = child.instance;
                        let p = child.placement;
                        model.move_instance(instance, &(conjugate(p) * p.inverse()))?;
                    }
                }
            }
            _ => {}
        }
    }
    // Under a reflection, faces are used the other way round in their shells.
    if t.is_negative() {
        let faces: std::collections::HashSet<NodeId> = nodes.iter().copied().filter(|n| n.kind == NodeKind::Face).collect();
        let shells: Vec<Shape> = nodes.iter().filter(|n| n.kind == NodeKind::Shell).map(|&n| Shape::from(n)).collect();
        for shell in shells {
            let children: Vec<Shape> = model.children_local(shell).collect();
            for child in children.into_iter().filter(|c| faces.contains(&c.node)) {
                model.replace_child(shell, child, child.reversed())?;
            }
        }
    }
    Ok(())
}

/// Transforms all geometry of the model in place, as [`transform_in_place`] does for one
/// shape.
#[doc(alias = "BRepGraph_Transform::Perform")]
pub fn transform_model(model: &mut Model, t: &Transform, options: TransformOptions) -> Result<()> {
    let mut nodes: Vec<NodeId> = Vec::new();
    for kind in NodeKind::ALL {
        for index in 0..model.slot_count(kind) {
            let node = NodeId { kind, index: index as u32 };
            if model.node_stable_id(node).is_some() {
                nodes.push(node);
            }
        }
    }
    transform_nodes(model, &nodes, t, options)
}

/// A copy of `shape` deformed by a general affine map. Curves and surfaces become
/// B-splines where their type is not preserved; meshes are transformed (reflections keep
/// their triangles' winding consistent); tolerances and mesh deflections grow by the
/// largest stretch factor of the map. The result has an identity placement.
#[doc(alias = "BRepBuilderAPI_GTransform")]
pub fn transform_affine(model: &mut Model, shape: Shape, t: &AffineTransform) -> Result<Shape> {
    let mut copier = AffineCopier { t: *t, stretch: t.max_stretch(), copied: HashMap::new(), edge_parameters: HashMap::new() };
    let copy = copier.copy(model, shape)?;
    let reversed = t.matrix().determinant() < 0.0 && shape.kind() == NodeKind::Face;
    let orientation = if reversed { shape.orientation.reversed() } else { shape.orientation };
    Ok(Shape { orientation, ..copy })
}

/// Copies a shape deformed by an affine map (its placement applied first).
struct AffineCopier {
    t: AffineTransform,
    stretch: f64,
    /// Copies of nodes reached at a placement.
    copied: HashMap<(NodeId, [[u64; 4]; 3]), NodeId>,
    /// How each copied edge's curve is parametrized against its original's.
    edge_parameters: HashMap<NodeId, CurveParameters>,
}

impl AffineCopier {
    /// The map for geometry stored under `placement`.
    fn map(&self, placement: &Transform) -> AffineTransform {
        let p = AffineTransform::from(*placement);
        AffineTransform::new(self.t.matrix() * p.matrix(), self.t.matrix() * p.translation() + self.t.translation())
    }

    fn copy(&mut self, model: &mut Model, shape: Shape) -> Result<Shape> {
        let key = (shape.node, shape.placement.rows().map(|r| r.map(f64::to_bits)));
        if let Some(&n) = self.copied.get(&key) {
            return Ok(Shape::new(n));
        }
        let map = self.map(&shape.placement);
        let local = Shape { placement: Transform::IDENTITY, orientation: Orientation::Forward, ..shape };
        let node: NodeId = match shape.kind() {
            NodeKind::Vertex => {
                let v = VertexId::try_from(shape.node).expect("a vertex");
                let record = model.get(v).ok_or_else(|| Error::InvalidInput(format!("no {v:?}")))?.clone();
                model.add_vertex(map.transform_point(record.point), record.tolerance * self.stretch).into()
            }
            NodeKind::Edge => {
                let e = EdgeId::try_from(shape.node).expect("an edge");
                let record = model.get(e).ok_or_else(|| Error::InvalidInput(format!("no {e:?}")))?.clone();
                let ends: Vec<Shape> = model.children_local(local).collect();
                let edge = if ends.len() == 2 {
                    let a = self.copy(model, Shape { placement: shape.placement * ends[0].placement, ..ends[0] })?;
                    let b = self.copy(model, Shape { placement: shape.placement * ends[1].placement, ..ends[1] })?;
                    let (a, b) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
                    let (a, b) = if ends[0].orientation == Orientation::Reversed { (b, a) } else { (a, b) };
                    model.add_edge(a, b, record.tolerance * self.stretch)?
                } else {
                    let edge = EdgeId::try_from(model.new_empty(NodeKind::Edge).node).expect("an edge");
                    model.set_edge_tolerance(edge, record.tolerance * self.stretch)?;
                    edge
                };
                if let Some((curve, first, last)) = model.edge_curve(local) {
                    let (curve, first, last, parameters) = affine_curve(&curve, first, last, &map)?;
                    model.set_edge_curve(edge.into(), curve, record.tolerance * self.stretch)?;
                    model.set_edge_range(edge.into(), first, last)?;
                    self.edge_parameters.insert(edge.into(), parameters);
                }
                if let Some(p) = model.edge_polyline(local).cloned() {
                    model.set_edge_polyline(edge.into(), Some(p.affine_transformed(&map)))?;
                }
                edge.into()
            }
            NodeKind::Face => {
                let f = crate::topology::FaceId::try_from(shape.node).expect("a face");
                let record = model.get(f).ok_or_else(|| Error::InvalidInput(format!("no {f:?}")))?.clone();
                let face = model.new_empty(NodeKind::Face);
                let tolerance = record.tolerance * self.stretch;
                model.set_face_tolerance(crate::topology::FaceId::try_from(face.node).expect("a face"), tolerance)?;
                let mut parameters = SurfaceParameters::Axes(None, None);
                if let Some(surface) = model.face_surface(f).cloned() {
                    let bounds = model.face_bounds(f);
                    let (mapped, p) = affine_surface(&surface, bounds, &map)?;
                    model.set_face_surface(face, mapped, tolerance)?;
                    parameters = p;
                }
                if let Some(mesh) = model.face_mesh(f).cloned() {
                    model.set_face_mesh(face, Some(mesh.affine_transformed(&map)))?;
                }
                let wires: Vec<Shape> = model.children_local(local).collect();
                for wire in wires {
                    let mut coedges = Vec::new();
                    for c in model.children_local(wire).collect::<Vec<_>>() {
                        let id = CoEdgeId::try_from(c.node).expect("a coedge");
                        let (edge, orientation) = (model[id].edge, model[id].orientation);
                        let placed = Shape { placement: shape.placement * wire.placement, ..Shape::from(edge) };
                        let new_edge = EdgeId::try_from(self.copy(model, placed)?.node).expect("an edge");
                        let coedge = model.add_coedge(new_edge, orientation)?;
                        if let Some((pcurve, first, last)) = model.effective_pcurve(id) {
                            let edge_parameters = self.edge_parameters.get(&new_edge.into()).copied().unwrap_or(CurveParameters::Same);
                            let (pcurve, first, last) = affine_pcurve(&pcurve, first, last, edge_parameters, &parameters)?;
                            model.set_pcurve_with_range(coedge, pcurve, first, last)?;
                        }
                        coedges.push(coedge);
                    }
                    let w = model.add_wire(&coedges)?;
                    model.add(face, Shape { orientation: wire.orientation, ..Shape::from(w) })?;
                }
                face.node
            }
            NodeKind::Shell | NodeKind::Solid | NodeKind::CompSolid | NodeKind::Compound => {
                let copy = model.new_empty(shape.kind());
                let reverse = self.t.matrix().determinant() < 0.0;
                for child in model.children_local(local).collect::<Vec<_>>() {
                    let placed = Shape { orientation: Orientation::Forward, placement: shape.placement * child.placement, ..child };
                    let c = self.copy(model, placed)?;
                    let flip = reverse && child.kind() == NodeKind::Face;
                    model.add(copy, Shape { orientation: if flip { child.orientation.reversed() } else { child.orientation }, ..c })?;
                }
                copy.node
            }
            other => return Err(Error::InvalidInput(format!("cannot transform a {other:?}"))),
        };
        self.copied.insert(key, node);
        Ok(Shape::new(node))
    }
}

/// How the curve of an edge copied by [`transform_affine`] is parametrized against the
/// original's.
#[derive(Debug, Clone, Copy)]
enum CurveParameters {
    Same,
    /// The copy's parameter is the original's times this (a line's, which is a distance).
    Scaled(f64),
    /// The B-spline form of an arc of a conic, whose own parameter is the original's.
    Arc(ArcParameters),
}

impl CurveParameters {
    /// The original's parameter at the copy's `t`.
    fn original(&self, t: f64) -> f64 {
        match self {
            CurveParameters::Same => t,
            CurveParameters::Scaled(k) => t / k,
            CurveParameters::Arc(arc) => arc.angle(t),
        }
    }

    /// The copy's parameter at the original's `t`.
    fn copied(&self, t: f64) -> f64 {
        match self {
            CurveParameters::Same => t,
            CurveParameters::Scaled(k) => t * k,
            CurveParameters::Arc(arc) => arc.parameter(t),
        }
    }
}

/// How the surface of a face copied by [`transform_affine`] is parametrized against the
/// original's.
#[derive(Debug, Clone, Copy)]
enum SurfaceParameters {
    /// `(u, v) ↦ A·(u, v)`: a plane mapped to a plane, its frame no longer orthonormal.
    Linear([[f64; 2]; 2]),
    /// u and v each on their own: kept (`None`) or the B-spline form of an arc in it.
    Axes(Option<ArcParameters>, Option<ArcParameters>),
}

impl SurfaceParameters {
    /// The copy's (u, v) at the original's `p`.
    fn copied(&self, p: Point2) -> Point2 {
        match self {
            SurfaceParameters::Linear(a) => Point2::new(a[0][0] * p.x + a[0][1] * p.y, a[1][0] * p.x + a[1][1] * p.y),
            SurfaceParameters::Axes(u, v) => Point2::new(u.map_or(p.x, |a| a.parameter(p.x)), v.map_or(p.y, |a| a.parameter(p.y))),
        }
    }
}

/// How the B-spline form of `curve` over `[a, b]` is parametrized, where it is not by the
/// curve's own parameter: an arc of a circle or ellipse, a piece of hyperbola.
fn arc_parameters(curve: &Curve, a: f64, b: f64) -> Option<ArcParameters> {
    match curve {
        Curve::Circle(_) | Curve::Ellipse(_) => Some(ArcParameters::elliptic(a, b)),
        Curve::Hyperbola(_) => Some(ArcParameters::hyperbolic(a, b)),
        Curve::Trimmed(t) => arc_parameters(t.basis(), a, b),
        _ => None,
    }
}

/// A curve deformed by an affine map, with its new range and how its parameter changed:
/// lines stay lines (their parameter, a distance, scaled), Bézier and B-spline curves keep
/// their parameters, other curves become B-splines first.
fn affine_curve(curve: &Curve, first: f64, last: f64, t: &AffineTransform) -> Result<(Curve, f64, f64, CurveParameters)> {
    use crate::geometry::curve::{Line, TrimmedCurve};
    use crate::math::Direction;
    Ok(match curve {
        Curve::Line(l) => {
            let d = t.transform_vector(crate::math::Vector::from(l.direction));
            let k = d.length();
            let direction = Direction::try_from_vector(d).ok_or_else(|| Error::InvalidInput("the map collapses a line".into()))?;
            (Curve::Line(Line::new(t.transform_point(l.origin), direction)), first * k, last * k, CurveParameters::Scaled(k))
        }
        Curve::Bezier(b) => (Curve::Bezier(b.affine_transformed(t)), first, last, CurveParameters::Same),
        Curve::BSpline(b) => (Curve::BSpline(b.affine_transformed(t)), first, last, CurveParameters::Same),
        Curve::Trimmed(c) if matches!(c.basis(), Curve::Bezier(_) | Curve::BSpline(_)) => {
            let (basis, _, _, _) = affine_curve(c.basis(), first, last, t)?;
            (Curve::Trimmed(TrimmedCurve::new(basis, c.first_parameter(), c.last_parameter())), first, last, CurveParameters::Same)
        }
        Curve::Trimmed(c) => affine_curve(c.basis(), first, last, t)?,
        other => {
            let piece = Curve::Trimmed(TrimmedCurve::new(other.clone(), first, last));
            let b = crate::geometry::convert::curve_to_bspline(&piece)
                .ok_or_else(|| Error::InvalidInput(format!("no B-spline form of {other:?}")))?;
            let (f, l) = (b.first_parameter(), b.last_parameter());
            let parameters = arc_parameters(other, first, last).map_or(CurveParameters::Same, CurveParameters::Arc);
            (Curve::BSpline(b.affine_transformed(t)), f, l, parameters)
        }
    })
}

/// A surface deformed by an affine map, and how its parameters changed: a plane stays a
/// plane, Bézier and B-spline surfaces map their poles, other surfaces become B-splines
/// over the face's bounds first.
fn affine_surface(surface: &Surface, (u1, u2, v1, v2): (f64, f64, f64, f64), t: &AffineTransform) -> Result<(Surface, SurfaceParameters)> {
    use crate::geometry::surface::Plane;
    use crate::math::{Direction, Frame, Vector};
    let same = SurfaceParameters::Axes(None, None);
    Ok(match surface {
        Surface::Bezier(b) => (Surface::Bezier(b.affine_transformed(t)), same),
        Surface::BSpline(b) => (Surface::BSpline(b.affine_transformed(t)), same),
        Surface::Plane(plane) => {
            // The frame's axes map to the new plane's u and v directions, no longer
            // orthonormal: a new frame along the first, and (u, v) mapped linearly.
            let frame = plane.frame();
            let (x, y) = (t.transform_vector(Vector::from(frame.x())), t.transform_vector(Vector::from(frame.y())));
            let (z, x_hint) = Direction::try_from_vector(x.cross(y)).zip(Direction::try_from_vector(x)).ok_or_else(|| Error::InvalidInput("the map collapses a plane".into()))?;
            let image = Frame::new(t.transform_point(frame.origin()), z, x_hint);
            let (ex, ey) = (Vector::from(image.x()), Vector::from(image.y()));
            let a = [[x.dot(ex), y.dot(ex)], [x.dot(ey), y.dot(ey)]];
            (Surface::Plane(Plane::from_frame(image)), SurfaceParameters::Linear(a))
        }
        other => {
            let (axes, basis) = surface_arcs(other, (u1, u2), (v1, v2));
            let b = crate::geometry::convert::surface_patch_to_bspline(basis, (u1, u2), (v1, v2))
                .ok_or_else(|| Error::InvalidInput("no B-spline form of the surface".into()))?;
            (Surface::BSpline(b.affine_transformed(t)), axes)
        }
    })
}

/// How [`surface_patch_to_bspline`](crate::geometry::convert::surface_patch_to_bspline)
/// parametrizes `surface` over the rectangle: in which of u and v its B-spline form runs
/// round an arc. Also the surface the conversion works on (a trimmed one's basis).
fn surface_arcs<'a>(surface: &'a Surface, (u1, u2): (f64, f64), (v1, v2): (f64, f64)) -> (SurfaceParameters, &'a Surface) {
    let arc = |a, b| Some(ArcParameters::elliptic(a, b));
    let axes = match surface {
        Surface::Cylinder(_) | Surface::Cone(_) => SurfaceParameters::Axes(arc(u1, u2), None),
        Surface::Sphere(_) | Surface::Torus(_) => SurfaceParameters::Axes(arc(u1, u2), arc(v1, v2)),
        Surface::Extrusion(e) => SurfaceParameters::Axes(arc_parameters(e.basis(), u1, u2), None),
        Surface::Revolution(r) => SurfaceParameters::Axes(arc(u1, u2), arc_parameters(r.basis(), v1, v2)),
        Surface::Trimmed(t) => return surface_arcs(t.basis(), (u1, u2), (v1, v2)),
        _ => SurfaceParameters::Axes(None, None),
    };
    (axes, surface)
}

/// The pcurve, over `[first, last]`, of an edge copied by [`transform_affine`] on a face
/// copied with it, with its new range: the original pcurve, run at the original edge's
/// parameter and carried into the new surface's parameters. Exact where that keeps it a
/// line, a B-spline, or the B-spline form of an arc; otherwise fitted.
fn affine_pcurve(pcurve: &Curve2, first: f64, last: f64, edge: CurveParameters, surface: &SurfaceParameters) -> Result<(Curve2, f64, f64)> {
    use crate::geometry::curve2::TrimmedCurve2;
    let (a, b) = (edge.copied(first), edge.copied(last));
    let map_poles = |c: &BSplineCurve2| {
        let mut c = c.clone();
        for i in 0..c.poles().len() {
            let p = surface.copied(c.pole(i));
            c.set_pole(i, p);
        }
        c
    };
    match (surface, edge, pcurve) {
        (SurfaceParameters::Axes(None, None), CurveParameters::Same, _) => return Ok((pcurve.clone(), first, last)),
        (SurfaceParameters::Linear(m), CurveParameters::Same | CurveParameters::Scaled(_), Curve2::Line(line)) => {
            // A line's parameter, a distance, scales as the edge's does (both are the length
            // the map gives a unit step along it).
            let d = Vector2::new(m[0][0] * line.direction.x() + m[0][1] * line.direction.y(), m[1][0] * line.direction.x() + m[1][1] * line.direction.y());
            let k = if let CurveParameters::Scaled(k) = edge { k } else { 1.0 };
            if (d.length() - k).abs() <= 1e-12 * k {
                return Ok((Curve2::Line(Line2::new(surface.copied(line.origin), Direction2::from_vector(d))), a, b));
            }
        }
        (SurfaceParameters::Linear(_), CurveParameters::Same, Curve2::BSpline(c)) => return Ok((Curve2::BSpline(map_poles(c)), a, b)),
        (SurfaceParameters::Linear(_), CurveParameters::Arc(arc), _) if arc.first == first => {
            // The pcurve of an arc is an arc over the same range, whose B-spline form runs
            // with the edge's.
            let mut basis = pcurve;
            while let Curve2::Trimmed(t) = basis {
                basis = t.basis();
            }
            if matches!(basis, Curve2::Circle(_) | Curve2::Ellipse(_) | Curve2::Hyperbola(_)) {
                let piece = Curve2::Trimmed(TrimmedCurve2::new(basis.clone(), first, last));
                if let Some(c) = crate::geometry::convert::curve2_to_bspline(&piece) {
                    return Ok((Curve2::BSpline(map_poles(&c)), a, b));
                }
            }
        }
        _ => {}
    }
    // Fitted, broken where the edge's B-spline form or the surface's has a knot.
    let along = |t: f64| pcurve.point(edge.original(t));
    let mut breaks: Vec<f64> = match edge {
        CurveParameters::Arc(arc) => arc.knots(),
        _ => Vec::new(),
    };
    if let SurfaceParameters::Axes(u, v) = surface {
        for (axis, arc) in [(0, u), (1, v)] {
            let Some(arc) = arc else { continue };
            let coordinate = |t: f64| if axis == 0 { along(t).x } else { along(t).y };
            for knot in arc.knots() {
                breaks.extend(crossings(&|t| coordinate(t) - knot, a, b));
            }
        }
    }
    let margin = 1e-9 * (b - a);
    breaks.retain(|&k| k > a + margin && k < b - margin);
    breaks.sort_by(f64::total_cmp);
    breaks.dedup_by(|x, y| (*x - *y).abs() <= margin);
    let options = crate::geometry::convert::CurveApproximationOptions { tolerance: 1e-10, max_segments: 256, ..Default::default() };
    let point = |t: f64| surface.copied(along(t));
    let fitted = crate::geometry::convert::approximate_pieces(&point, a, b, &breaks, 5, &options)
        .ok_or_else(|| Error::Failed("no pcurve on the mapped surface within 1e-10".into()))?;
    Ok((Curve2::BSpline(fitted.curve), a, b))
}

/// Where `f` changes sign over `[a, b]` (from one side of 0 to the other, not where it
/// only touches or runs along 0), to the last bit: the samples' sign changes, bisected.
fn crossings(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> Vec<f64> {
    const SAMPLES: usize = 256;
    let mut found = Vec::new();
    // The last sample off 0.
    let mut last: Option<(f64, bool)> = None;
    for i in 0..=SAMPLES {
        let t = a + (b - a) * i as f64 / SAMPLES as f64;
        let value = f(t);
        if value == 0.0 {
            continue;
        }
        if let Some((from, negative)) = last.filter(|&(_, negative)| negative != (value < 0.0)) {
            let (mut lo, mut hi) = (from, t);
            for _ in 0..100 {
                let middle = 0.5 * (lo + hi);
                if middle <= lo || middle >= hi {
                    break;
                }
                let v = f(middle);
                if v != 0.0 && (v < 0.0) == negative { lo = middle } else { hi = middle }
            }
            found.push(0.5 * (lo + hi));
        }
        last = Some((t, value < 0.0));
    }
    found
}

/// A deep copy of `shape`: new nodes all the way down, with geometry (and, if
/// `copy_meshes`, meshes) duplicated, so that editing the copy leaves `shape` untouched.
#[doc(alias = "BRepBuilderAPI_Copy")]
pub fn copy(model: &mut Model, shape: Shape, copy_meshes: bool) -> Result<Shape> {
    let copy = Copier { meshes: copy_meshes, ..Copier::default() }.copy(model, Shape { orientation: Orientation::Forward, ..shape })?;
    Ok(Shape { orientation: shape.orientation, ..copy })
}

/// [`copy`] of `shape` with all its placements applied to the geometry (the copy placed
/// as it is, nowhere), and the copies of `edges` (edges of `shape`, as placed in it).
pub(crate) fn copy_unplaced(model: &mut Model, shape: Shape, edges: &[Shape]) -> Result<(Shape, Vec<Shape>)> {
    let mut copier = Copier::default();
    let copy = copier.copy(model, Shape { orientation: Orientation::Forward, ..shape })?;
    let copies = edges
        .iter()
        .map(|e| {
            let found = copier.edges.get(&e.node).and_then(|l| l.iter().find(|(t, _)| *t == e.placement));
            let &(_, id) = found.ok_or_else(|| Error::InvalidInput("an edge not of the shape".into()))?;
            Ok(Shape { orientation: e.orientation, ..Shape::from(id) })
        })
        .collect::<Result<Vec<Shape>>>()?;
    Ok((Shape { orientation: shape.orientation, ..copy }, copies))
}

/// How a curve's parameter changes when the curve is scaled by `s` (OCCT's
/// `Geom_Curve::TransformedParameter`): a distance along a line scales, an angle does not.
fn parameter_scale(curve: &Curve, s: f64) -> f64 {
    match curve {
        Curve::Line(_) | Curve::Parabola(_) => s.abs(),
        Curve::Trimmed(c) => parameter_scale(c.basis(), s),
        Curve::Offset(c) => parameter_scale(c.basis(), s),
        Curve::Circle(_) | Curve::Ellipse(_) | Curve::Hyperbola(_) | Curve::Bezier(_) | Curve::BSpline(_) => 1.0,
    }
}

/// How a surface's (u, v) change when the surface is scaled by `s` (OCCT's
/// `Geom_Surface::ParametricTransformation`), as factors on u and v.
fn parameter_scales(surface: &Surface, s: f64) -> (f64, f64) {
    let s = s.abs();
    match surface {
        Surface::Plane(_) => (s, s),
        Surface::Cylinder(_) | Surface::Cone(_) => (1.0, s),
        Surface::Extrusion(e) => (parameter_scale(e.basis(), s), s),
        Surface::Revolution(r) => (1.0, parameter_scale(r.basis(), s)),
        Surface::Trimmed(t) => parameter_scales(t.basis(), s),
        Surface::Offset(o) => parameter_scales(o.basis(), s),
        Surface::Sphere(_) | Surface::Torus(_) | Surface::Bezier(_) | Surface::BSpline(_) => (1.0, 1.0),
    }
}

/// A pcurve with its (u, v) scaled by `(su, sv)`, where that keeps it the same kind of
/// curve (`None` otherwise).
fn scaled_pcurve(pcurve: &Curve2, (su, sv): (f64, f64)) -> Option<Curve2> {
    if su == 1.0 && sv == 1.0 {
        return Some(pcurve.clone());
    }
    if su == sv {
        return Some(pcurve.transformed(&Transform2::scaling(Point2::ORIGIN, su)));
    }
    match pcurve {
        Curve2::Line(line) => {
            // Through the images of two of its points (its parameter stays a distance, which
            // matches the edge's only along u or v: checked by the caller).
            let scale = |p: Point2| Point2::new(p.x * su, p.y * sv);
            let (a, b) = (scale(line.point(0.0)), scale(line.point(1.0)));
            Some(Curve2::Line(Line2::new(a, Direction2::from_vector(b - a))))
        }
        // An affine map of the poles maps the curve, its parameter unchanged.
        Curve2::BSpline(b) => {
            let poles = b.poles().iter().map(|p| Point2::new(p.x * su, p.y * sv)).collect();
            let (knots, multiplicities) = (b.knots().to_vec(), b.multiplicities().to_vec());
            let mut scaled = if b.is_periodic() {
                crate::geometry::curve2::BSplineCurve2::new_periodic(poles, knots, multiplicities, b.degree())
            } else {
                crate::geometry::curve2::BSplineCurve2::new(poles, knots, multiplicities, b.degree())
            };
            if let Some(w) = b.weights() {
                scaled = scaled.with_weights(w.to_vec());
            }
            Some(Curve2::BSpline(scaled))
        }
        _ => None,
    }
}

/// Copies topology with its geometry transformed by each shape's placement, sharing what
/// the original shares.
#[derive(Default)]
struct Copier {
    /// Use the copied faces reversed (under a reflection).
    reverse_faces: bool,
    /// Copy the meshes of faces with exact geometry too (mesh-only faces always keep theirs).
    meshes: bool,
    vertices: HashMap<NodeId, Vec<(Transform, VertexId)>>,
    edges: HashMap<NodeId, Vec<(Transform, EdgeId)>>,
}

impl Copier {
    /// A copy of the node of `shape` (forward, its placement applied to the geometry).
    fn copy(&mut self, model: &mut Model, shape: Shape) -> Result<Shape> {
        match shape.kind() {
            NodeKind::Vertex => Ok(self.vertex(model, shape).into()),
            NodeKind::Edge => Ok(self.edge(model, shape)?.0.into()),
            NodeKind::Face => self.face(model, shape),
            NodeKind::Wire => {
                // A wire outside a face: its edges copied, in order (no pcurves to carry).
                let coedges: Vec<Shape> = model.children_local(shape).collect();
                let mut ids = Vec::new();
                for c in coedges {
                    let id = CoEdgeId::try_from(c.node).expect("a wire holds coedges");
                    let (edge, orientation) = (model[id].edge, model[id].orientation);
                    let placed = Shape { placement: shape.placement * c.placement, ..Shape::from(edge) };
                    let (new_edge, _) = self.edge(model, placed)?;
                    ids.push(model.add_coedge(new_edge, orientation)?);
                }
                Ok(model.add_wire(&ids)?.into())
            }
            NodeKind::Shell | NodeKind::Solid | NodeKind::CompSolid | NodeKind::Compound => {
                let copy = model.new_empty(shape.kind());
                let children: Vec<Shape> = model.children_local(shape).collect();
                for child in children {
                    let placed = Shape { orientation: Orientation::Forward, placement: shape.placement * child.placement, ..child };
                    let c = self.copy(model, placed)?;
                    let reversed = self.reverse_faces && child.kind() == NodeKind::Face;
                    let orientation = if reversed { child.orientation.reversed() } else { child.orientation };
                    model.add(copy, Shape { orientation, ..c })?;
                }
                Ok(copy)
            }
            other => Err(Error::InvalidInput(format!("cannot transform a {other:?}"))),
        }
    }

    fn vertex(&mut self, model: &mut Model, shape: Shape) -> VertexId {
        if let Some(&(_, v)) = self.vertices.get(&shape.node).and_then(|l| l.iter().find(|(t, _)| *t == shape.placement)) {
            return v;
        }
        let tolerance = model.tolerance(shape) * shape.placement.scale_factor().abs().max(1.0);
        let v = model.add_vertex(model.point(shape), tolerance);
        self.vertices.entry(shape.node).or_default().push((shape.placement, v));
        v
    }

    /// The copy of an edge, and the factor its parameters were scaled by.
    fn edge(&mut self, model: &mut Model, shape: Shape) -> Result<(EdgeId, f64)> {
        let forward = Shape { orientation: Orientation::Forward, ..shape };
        let original = model.edge_curve(Shape { placement: Transform::IDENTITY, ..forward });
        let scale = original.as_ref().map_or(1.0, |(c, _, _)| parameter_scale(c, shape.placement.scale_factor()));
        if let Some(&(_, e)) = self.edges.get(&shape.node).and_then(|l| l.iter().find(|(t, _)| *t == shape.placement)) {
            return Ok((e, scale));
        }
        let tolerance = model.tolerance(shape) * shape.placement.scale_factor().abs().max(1.0);
        let e = match model.children_local(forward).count() {
            2 => {
                let (a, b) = model.edge_vertices(forward);
                let (a, b) = (self.vertex(model, a), self.vertex(model, b));
                model.add_edge(a, b, tolerance)?
            }
            // An edge without vertices (a polyline only).
            _ => {
                let e = EdgeId::try_from(model.new_empty(NodeKind::Edge).node).expect("an edge");
                model.set_edge_tolerance(e, tolerance)?;
                e
            }
        };
        let had_curve = original.is_some();
        if let Some((curve, first, last)) = original {
            model.set_edge_curve(e.into(), curve.transformed(&shape.placement), tolerance)?;
            model.set_edge_range(e.into(), first * scale, last * scale)?;
        }
        if self.meshes || !had_curve {
            let local = Shape { placement: Transform::IDENTITY, ..forward };
            if let Some(p) = model.edge_polyline(local).cloned() {
                model.set_edge_polyline(e.into(), Some(p.affine_transformed(&AffineTransform::from(shape.placement))))?;
            }
        }
        self.edges.entry(shape.node).or_default().push((shape.placement, e));
        Ok((e, scale))
    }

    fn face(&mut self, model: &mut Model, shape: Shape) -> Result<Shape> {
        let local = Shape { placement: Transform::IDENTITY, ..shape };
        let mesh = model.face_mesh(local).cloned().map(|m| m.affine_transformed(&AffineTransform::from(shape.placement)));
        let Some(original) = model.surface(local).cloned() else {
            // A face carrying only a mesh: the mesh moves.
            let copy = model.new_empty(NodeKind::Face);
            model.set_face_mesh(copy, mesh)?;
            return Ok(copy);
        };
        let surface = original.transformed(&shape.placement);
        let uv_scales = parameter_scales(&original, shape.placement.scale_factor());
        let tolerance = model.tolerance(shape) * shape.placement.scale_factor().abs().max(1.0);
        let copy = model.new_empty(NodeKind::Face);
        model.set_face_surface(copy, surface.clone(), tolerance)?;
        if self.meshes {
            model.set_face_mesh(copy, mesh)?;
        }
        let wires: Vec<Shape> = model.children_local(shape).collect();
        for wire in wires {
            // Each coedge's new edge, orientation, range, pcurve and (where it was refitted)
            // the shape's own pcurve in the new parameters.
            let mut uses: Vec<(EdgeId, Orientation, f64, f64, Curve2, Option<Curve2>)> = Vec::new();
            let coedge_shapes: Vec<Shape> = model.children_local(wire).collect();
            for c in coedge_shapes {
                let id = CoEdgeId::try_from(c.node).expect("a wire holds coedges");
                let (edge, orientation) = (model[id].edge, model[id].orientation);
                let placed = Shape { placement: shape.placement * wire.placement * c.placement, ..Shape::from(edge) };
                let (new_edge, scale) = self.edge(model, placed)?;
                let (pcurve, first, last) = model
                    .effective_pcurve(id)
                    .ok_or_else(|| Error::InvalidInput("transforming a face without pcurves".into()))?;
                let (first, last) = (first * scale, last * scale);
                // The pcurve in the new parameters, checked against the edge; refitted where
                // its own parametrization changes differently.
                let curve = model.edge_curve(new_edge);
                let fits = |p: &Curve2| match &curve {
                    None => true,
                    Some((c, _, _)) => (0..=4).all(|i| {
                        let t = first + (last - first) * i as f64 / 4.0;
                        let q = p.point(t);
                        // A wrong parameter mapping is off by far more than the fit of an
                        // approximate pcurve.
                        surface.point(q.x, q.y).distance(c.point(t)) <= (10.0 * tolerance).max(1e-6)
                    }),
                };
                let scaled = scaled_pcurve(&pcurve, uv_scales);
                let (pcurve, own) = match scaled.clone().filter(|p| fits(p)) {
                    Some(p) => (p, None),
                    None => {
                        let (c, _, _) = curve.as_ref().ok_or_else(|| Error::InvalidInput("no pcurve for a degenerate edge".into()))?;
                        // (Refitted as closely as it goes; else, where the shape's own pcurve
                        // strays from its edge already, that one as it is.)
                        match project::pcurve_on_range_within(c, &surface, first, last, &[1e-9, 1e-7, 1e-5]) {
                            // (In the same period as the one it replaces: on a periodic surface
                            // the fit may come out a period off, on the seam's other side.)
                            Some(p) => {
                                let p = Curve2::BSpline(p);
                                let by = scaled.as_ref().map_or(Vector2::new(0.0, 0.0), |old| {
                                    let (was, now) = (old.point(first), p.point(first));
                                    let shift = |d: f64, period: Option<f64>| period.map_or(0.0, |period| (d / period).round() * period);
                                    Vector2::new(shift(was.x - now.x, surface.u_period()), shift(was.y - now.y, surface.v_period()))
                                });
                                (if by.length() > 0.0 { p.transformed(&Transform2::translation(by)) } else { p }, scaled)
                            }
                            None => (scaled.ok_or_else(|| Error::InvalidInput("no pcurve for a transformed edge".into()))?, None),
                        }
                    }
                };
                uses.push((new_edge, orientation, first, last, pcurve, own));
            }
            // (Refitted pcurves that leave the loop open where the shape's own close it (an
            // edge whose curve strays from them within its loose tolerance, half a turn round
            // a cylinder): the shape's own for the loop.)
            let closes = |pcurves: &[&Curve2]| -> bool {
                let n = uses.len();
                (0..n).all(|i| {
                    let j = (i + 1) % n;
                    let end = |k: usize, start: bool| {
                        let (_, orientation, first, last, ..) = &uses[k];
                        let at_start = start != (*orientation == Orientation::Reversed);
                        let uv = pcurves[k].point(if at_start { *first } else { *last });
                        surface.point(uv.x, uv.y)
                    };
                    end(i, false).distance(end(j, true)) <= (10.0 * tolerance).max(crate::tolerance::LINEAR)
                })
            };
            let chosen: Vec<&Curve2> = uses.iter().map(|u| &u.4).collect();
            let own: Vec<&Curve2> = uses.iter().map(|u| u.5.as_ref().unwrap_or(&u.4)).collect();
            let keep_own = uses.iter().any(|u| u.5.is_some()) && !closes(&chosen) && closes(&own);
            let mut ids = Vec::new();
            for (k, (new_edge, orientation, first, last, ..)) in uses.iter().enumerate() {
                let pcurve = if keep_own { own[k].clone() } else { chosen[k].clone() };
                let coedge = model.add_coedge(*new_edge, *orientation)?;
                model.set_pcurve_with_range(coedge, pcurve, *first, *last)?;
                ids.push(coedge);
            }
            let w = model.add_wire(&ids)?;
            model.add(copy, Shape { orientation: wire.orientation, ..Shape::from(w) })?;
        }
        Ok(copy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Point, Vector};
    use crate::primitives::{Cone, Cuboid, Cylinder, Sphere, Torus};
    use crate::{boolean, properties, validity};

    #[test]
    fn scaled_primitives() {
        let mut model = Model::new();
        let t = Transform::scaling(Point::new(1.0, 2.0, 3.0), 10.0) * Transform::translation(Vector::new(0.5, 0.0, 0.0));
        let shapes = [
            Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap(),
            Cylinder::new(1.0, 2.0).with_angle(1.0).build(&mut model).unwrap(),
            Cone::new(2.0, 1.0, 3.0).build(&mut model).unwrap(),
            Sphere::new(1.5).with_latitudes(-1.0, 1.2).with_angle(2.0).build(&mut model).unwrap(),
            Torus::new(3.0, 1.0).with_angle(4.0).build(&mut model).unwrap(),
        ];
        for shape in shapes {
            let (v, a) = (properties::volume(&model, shape).mass, properties::area(&model, shape).mass);
            let scaled = transform(&mut model, shape, &t, TransformOptions::default()).unwrap();
            let report = validity::check(&model, scaled);
            assert!(report.problems.is_empty(), "{report:?}");
            let (v2, a2) = (properties::volume(&model, scaled).mass, properties::area(&model, scaled).mass);
            assert!((v2 - 1000.0 * v).abs() < 1e-9 * v2, "{v2} != 1000 x {v}");
            assert!((a2 - 100.0 * a).abs() < 1e-9 * a2, "{a2} != 100 x {a}");
            // Booleans see the new geometry.
            let corner = model.unique_subshapes(scaled, NodeKind::Vertex).iter().map(|&v| model.point(v)).next().unwrap();
            let block = Cuboid::from_corners(Point::new(-100.0, -100.0, -100.0), Point::new(100.0, 100.0, corner.z + 0.37)).build(&mut model).unwrap();
            let common = boolean::intersection(&mut model, scaled, block).unwrap();
            let report = validity::check(&model, common);
            assert!(report.problems.is_empty(), "{report:?}");
        }
    }

    /// Stretches, shears and mirrors (a negative determinant), and one shear that also
    /// turns the parts.
    fn stretches() -> Vec<AffineTransform> {
        use crate::math::Matrix3;
        let rows = [
            [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]],
            [[1.3, 0.4, 0.0], [0.0, 0.8, 0.3], [0.2, 0.0, 1.1]],
            [[-1.2, 0.3, 0.1], [0.2, 0.7, -0.4], [0.0, 0.5, 1.6]],
            [[0.34, 0.05, -1.56], [-0.9, 0.6, 0.2], [0.4, 1.2, 0.3]],
        ];
        rows.iter().map(|r| AffineTransform::new(Matrix3::from_rows(*r), Vector::new(3.0, -1.0, 2.0))).collect()
    }

    /// Primitives with the closed form of how far inside each point is (negative
    /// outside): a box, a cylinder, a frustum, a cone to its apex, a ball, a ring.
    fn primitives_with_depth(model: &mut Model) -> Vec<(Shape, Box<dyn Fn(Point) -> f64>)> {
        let radial = |p: Point| (p.x * p.x + p.y * p.y).sqrt();
        vec![
            (Cuboid::new(1.0, 2.0, 3.0).build(model).unwrap(), Box::new(|p: Point| [p.x, 1.0 - p.x, p.y, 2.0 - p.y, p.z, 3.0 - p.z].into_iter().fold(f64::MAX, f64::min))),
            (Cylinder::new(1.5, 2.0).build(model).unwrap(), Box::new(move |p: Point| (1.5 - radial(p)).min(p.z).min(2.0 - p.z))),
            // Within the side's slope the depth below the side is scaled, which keeps its sign.
            (Cone::new(2.0, 1.0, 3.0).build(model).unwrap(), Box::new(move |p: Point| (2.0 - p.z / 3.0 - radial(p)).min(p.z).min(3.0 - p.z))),
            (Cone::new(1.5, 0.0, 2.5).build(model).unwrap(), Box::new(move |p: Point| (1.5 * (1.0 - p.z / 2.5) - radial(p)).min(p.z).min(2.5 - p.z))),
            (Sphere::new(1.5).build(model).unwrap(), Box::new(|p: Point| 1.5 - (p - Point::ORIGIN).length())),
            (Torus::new(3.0, 1.0).build(model).unwrap(), Box::new(move |p: Point| 1.0 - ((radial(p) - 3.0).powi(2) + p.z * p.z).sqrt())),
        ]
    }

    #[test]
    fn affine_maps_of_primitives_are_valid_with_the_volumes_and_centroids_they_map_to() {
        let mut model = Model::new();
        let solids = primitives_with_depth(&mut model);
        for t in stretches() {
            let det = t.matrix().determinant();
            for (shape, _) in &solids {
                let mapped = transform_affine(&mut model, *shape, &t).unwrap();
                let report = validity::check(&model, mapped);
                assert!(report.problems.is_empty(), "{report:?}");
                let (v, v2) = (properties::volume(&model, *shape), properties::volume(&model, mapped));
                assert!((v2.mass - det.abs() * v.mass).abs() < 1e-7 * v2.mass, "volume {} against {}", v2.mass, det.abs() * v.mass);
                assert!(v2.center_of_mass.distance(t.transform_point(v.center_of_mass)) < 1e-6, "centroid");
            }
        }
    }

    /// Points of a primitive's box, mapped, lie in its image exactly where they lie in it
    /// (the images' faces B-splines with seams and poles: the classifier's feet and rays
    /// there).
    #[test]
    fn affine_maps_of_primitives_place_points_as_their_preimages() {
        use crate::boolean::{PointClassifier, PointLocation};
        let mut model = Model::new();
        let solids = primitives_with_depth(&mut model);
        for t in stretches() {
            for (shape, depth) in &solids {
                let mapped = transform_affine(&mut model, *shape, &t).unwrap();
                let bounds = crate::properties::bounding_box(&model, *shape);
                let (lo, hi) = (bounds.min, bounds.max);
                let classifier = PointClassifier::new(&model, mapped).unwrap();
                let mut seed = 0x2545_f491_4f6c_dd1du64;
                let mut next = || {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    (seed >> 11) as f64 / (1u64 << 53) as f64
                };
                for _ in 0..600 {
                    let p = Point::new(lo.x + (hi.x - lo.x) * next(), lo.y + (hi.y - lo.y) * next(), lo.z + (hi.z - lo.z) * next());
                    let d = depth(p);
                    if d.abs() < 1e-4 {
                        continue;
                    }
                    let q = t.transform_point(p);
                    let want = if d > 0.0 { PointLocation::Inside } else { PointLocation::Outside };
                    assert_eq!(classifier.locate(q), want, "{p:?} at depth {d} mapped to {q:?}");
                }
            }
        }
    }
}

/// Which sub-shapes [`set_tolerances`] and [`limit_tolerances`] change.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ToleranceTarget {
    #[default]
    All,
    Vertices,
    Edges,
    Faces,
}

/// The vertices, edges and faces of `shape` (each once) that `target` names.
fn tolerance_targets(model: &Model, shape: Shape, target: ToleranceTarget) -> Vec<Shape> {
    let kinds: &[NodeKind] = match target {
        ToleranceTarget::All => &[NodeKind::Vertex, NodeKind::Edge, NodeKind::Face],
        ToleranceTarget::Vertices => &[NodeKind::Vertex],
        ToleranceTarget::Edges => &[NodeKind::Edge],
        ToleranceTarget::Faces => &[NodeKind::Face],
    };
    kinds.iter().flat_map(|&kind| model.unique_subshapes(shape, kind)).collect()
}

fn set_tolerance(model: &mut Model, s: Shape, tolerance: f64) -> Result<()> {
    match s.kind() {
        NodeKind::Vertex => model.set_vertex_tolerance(VertexId::try_from(s.node).expect("a vertex"), tolerance),
        NodeKind::Edge => model.set_edge_tolerance(EdgeId::try_from(s.node).expect("an edge"), tolerance),
        NodeKind::Face => model.set_face_tolerance(crate::topology::FaceId::try_from(s.node).expect("a face"), tolerance),
        _ => Ok(()),
    }
}

/// Sets the tolerance of every vertex, edge and face of `shape` (or of one kind of them)
/// to `tolerance`.
#[doc(alias = "ShapeFix_ShapeTolerance::SetTolerance")]
pub fn set_tolerances(model: &mut Model, shape: Shape, tolerance: f64, target: ToleranceTarget) -> Result<()> {
    for s in tolerance_targets(model, shape, target) {
        set_tolerance(model, s, tolerance)?;
    }
    Ok(())
}

/// Brings the tolerances of the vertices, edges and faces of `shape` (or of one kind of
/// them) into `[min, max]`. Returns whether any changed.
#[doc(alias = "ShapeFix_ShapeTolerance::LimitTolerance")]
pub fn limit_tolerances(model: &mut Model, shape: Shape, min: f64, max: f64, target: ToleranceTarget) -> Result<bool> {
    let mut changed = false;
    for s in tolerance_targets(model, shape, target) {
        let tolerance = model.tolerance(s);
        let limited = tolerance.clamp(min, max.max(min));
        if limited != tolerance {
            set_tolerance(model, s, limited)?;
            changed = true;
        }
    }
    Ok(changed)
}

/// Makes the tolerances of `shape` consistent: an edge's at least those of its faces, a
/// vertex's at least those of its edges and its distance from their curves' ends.
#[doc(alias = "BRepLib::UpdateTolerances")]
pub fn update_tolerances(model: &mut Model, shape: Shape) -> Result<()> {
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let face_tolerance = model.tolerance(face);
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            if model.tolerance(edge) < face_tolerance {
                set_tolerance(model, edge, face_tolerance)?;
            }
        }
    }
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        let edge = Shape { orientation: Orientation::Forward, ..edge };
        let edge_tolerance = model.tolerance(edge);
        let (start, end) = model.edge_vertices(edge);
        let ends = model.edge_curve(edge).map(|(curve, first, last)| (curve.point(first), curve.point(last)));
        for (vertex, end_point) in [(start, ends.map(|e| e.0)), (end, ends.map(|e| e.1))] {
            let gap = end_point.map_or(0.0, |p| model.point(vertex).distance(p));
            let needed = edge_tolerance.max(gap);
            if model.tolerance(vertex) < needed {
                set_tolerance(model, vertex, needed)?;
            }
        }
    }
    Ok(())
}
