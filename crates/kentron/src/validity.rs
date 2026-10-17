//! Checking that a shape is a valid boundary representation.

use std::collections::{HashMap, HashSet};

use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Vector};
use crate::tolerance::{LINEAR, PARAMETRIC};
use crate::properties;
use crate::topology::{
    CoEdgeId, EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape, ShellId, SolidId, VertexId,
    VertexUseId, WireId,
};

/// One defect found by a check.
#[doc(alias = "BRepGraph_Validate::Issue")]
#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    /// The node the problem was found on (the owner, for a use or slot record).
    pub node: NodeId,
    pub kind: ProblemKind,
    pub severity: Severity,
    pub description: String,
}

/// Whether a problem makes the shape or model invalid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Severity {
    Error,
    /// Suspicious but allowed (e.g. a face with a surface but no boundary).
    Warning,
}

/// What kind of defect a [`Problem`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProblemKind {
    /// An id refers to a removed or nonexistent record, or to a record of the wrong kind.
    DanglingReference,
    /// A use entry or slot disagrees with the record it belongs to.
    InconsistentLink,
    /// A live geometry or mesh slot that no record references. (A slot referenced by
    /// several records is fine: that is shared geometry, see `Model::deduplicate`.)
    OrphanSlot,
    /// Consecutive coedges of a wire do not share a vertex.
    DisconnectedWire,
    /// A face has a surface but no wires.
    FaceWithoutBoundary,
    /// A shell bounding a solid is open: an edge is not shared by exactly two of its faces,
    /// in opposite directions (OCCT's `BRepCheck_NotClosed`, `BadOrientationOfSubshape`).
    NotClosed,
    /// A product graph cycle, or a root product that is also instantiated.
    ProductStructure,
    /// A face's boundary crosses itself: two of its edges (or one) cross away from a vertex
    /// they share (OCCT's `BRepCheck_SelfIntersectingWire`).
    SelfIntersection,
    /// Any geometric defect found by [`check`].
    Geometry,
}

/// The problems found; [`Report::is_valid`] when there is no error.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Report {
    pub problems: Vec<Problem>,
}

impl Report {
    /// True if no problem is an error (warnings are allowed).
    pub fn is_valid(&self) -> bool {
        self.errors().next().is_none()
    }

    pub fn errors(&self) -> impl Iterator<Item = &Problem> {
        self.problems.iter().filter(|p| p.severity == Severity::Error)
    }

    pub fn warnings(&self) -> impl Iterator<Item = &Problem> {
        self.problems.iter().filter(|p| p.severity == Severity::Warning)
    }

    /// True if some problem of `kind` has `severity`.
    pub fn has(&self, kind: ProblemKind, severity: Severity) -> bool {
        self.problems.iter().any(|p| p.kind == kind && p.severity == severity)
    }
}

/// Two sub-shapes of a shape that meet although they share no vertex (see
/// [`self_interferences`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interference {
    pub first: Shape,
    pub second: Shape,
}

/// The pairs of `shape`'s vertices, edges and faces that touch or cross although they
/// share no vertex: within the sum of their tolerances of each other. Vertices first, then
/// edges, then faces; each pair in that order. (Sub-shapes sharing a vertex are not
/// checked, where OCCT's checker also finds them crossing away from it.)
#[doc(alias = "BOPAlgo_CheckerSI")]
#[doc(alias = "bopcheck")]
pub fn self_interferences(model: &mut Model, shape: Shape) -> crate::error::Result<Vec<Interference>> {
    struct Item {
        shape: Shape,
        vertices: Vec<NodeId>,
        bbox: crate::math::BoundingBox,
        tolerance: f64,
    }
    let mut items = Vec::new();
    for kind in [NodeKind::Vertex, NodeKind::Edge, NodeKind::Face] {
        for s in model.unique_subshapes(shape, kind) {
            if kind == NodeKind::Edge && (model.is_degenerate(s) || model.edge_curve(s).is_none()) {
                continue;
            }
            let tolerance = model.tolerance(s);
            let vertices = model.unique_subshapes(s, NodeKind::Vertex).into_iter().map(|v| v.node).collect();
            let bbox = properties::bounding_box(model, s).enlarged(tolerance);
            items.push(Item { shape: s, vertices, bbox, tolerance });
        }
    }
    let mut found = Vec::new();
    for i in 0..items.len() {
        for j in i + 1..items.len() {
            let (a, b) = (&items[i], &items[j]);
            if !a.bbox.overlaps(&b.bbox) || a.vertices.iter().any(|v| b.vertices.contains(v)) {
                continue;
            }
            let (first, second, reach) = (a.shape, b.shape, a.tolerance + b.tolerance);
            if crate::distance::distance_within(model, first, second, reach)?.value <= reach {
                found.push(Interference { first, second });
            }
        }
    }
    Ok(found)
}

/// Checks the topology and geometry of `shape` and everything below it: closed wires and
/// shells, consistent orientations, vertices and pcurves within tolerance of the curves and
/// surfaces they lie on, ...
#[doc(alias = "BRepCheck_Analyzer")]
pub fn check(model: &Model, shape: Shape) -> Report {
    let mut checker = Checker { model, report: Report::default() };
    let mut seen = HashSet::new();
    for occurrence in std::iter::once(shape).chain(model.explore_all(shape)) {
        if !seen.insert(occurrence.node) {
            continue;
        }
        let node = occurrence.node;
        match node.kind {
            NodeKind::Vertex => checker.vertex(VertexId::from_index(node.index)),
            NodeKind::Edge => checker.edge(EdgeId::from_index(node.index)),
            NodeKind::CoEdge => checker.coedge(CoEdgeId::from_index(node.index)),
            NodeKind::Wire => checker.wire(WireId::from_index(node.index)),
            NodeKind::Face => checker.face(FaceId::from_index(node.index)),
            NodeKind::Shell => checker.shell(ShellId::from_index(node.index)),
            NodeKind::Solid => checker.solid(occurrence),
            // A compound only groups shapes, which are checked in turn.
            NodeKind::Compound => {}
            kind => todo!("checking a {kind:?}"),
        }
    }
    checker.report
}

/// The distance from `p` to a surface, where there is a closed form. (Elsewhere a
/// vertex is still checked through its edges: each lies on its edge's curve, which the
/// same-parameter check ties to the surface.)
fn surface_distance(surface: &Surface, p: Point) -> Option<f64> {
    Some(match surface {
        Surface::Plane(plane) => plane.distance(p),
        Surface::Cylinder(cylinder) => {
            let (_, v) = cylinder.parameters(p);
            let axis_point = cylinder.frame().origin() + Vector::from(cylinder.frame().z()) * v;
            ((p - axis_point).length() - cylinder.radius()).abs()
        }
        Surface::Sphere(sphere) => (p.distance(sphere.frame().origin()) - sphere.radius()).abs(),
        Surface::Cone(cone) => {
            // In the half-plane through the axis and `p`, the distance to the generating
            // line {(R + v·sin α, v·cos α)}, or past the apex (where R + v·sin α < 0) to its
            // mirror image {(-(R + v·sin α), v·cos α)}.
            let d = p - cone.frame().origin();
            let height = d.dot(Vector::from(cone.frame().z()));
            let radial = (d - Vector::from(cone.frame().z()) * height).length();
            let (sin, cos) = cone.semi_angle_sin_cos();
            ((radial - cone.radius()) * cos - height * sin).abs().min(((radial + cone.radius()) * cos + height * sin).abs())
        }
        Surface::Torus(torus) => {
            let d = p - torus.frame().origin();
            let height = d.dot(Vector::from(torus.frame().z()));
            let radial = (d - Vector::from(torus.frame().z()) * height).length();
            ((radial - torus.major_radius()).hypot(height) - torus.minor_radius()).abs()
        }
        _ => return None,
    })
}

struct Checker<'a> {
    model: &'a Model,
    report: Report,
}

impl Checker<'_> {
    fn problem(&mut self, node: impl Into<NodeId>, kind: ProblemKind, severity: Severity, description: String) {
        self.report.problems.push(Problem { node: node.into(), kind, severity, description });
    }

    fn error(&mut self, node: impl Into<NodeId>, kind: ProblemKind, description: String) {
        self.problem(node, kind, Severity::Error, description);
    }

    fn vertex(&mut self, id: VertexId) {
        let Some(v) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        let p = v.point;
        if ![p.x, p.y, p.z, v.tolerance].iter().all(|c| c.is_finite()) || v.tolerance < 0.0 {
            self.error(id, ProblemKind::Geometry, format!("{id:?} has a non-finite point or bad tolerance"));
        }
    }

    /// The vertex at the start or end of an edge, if its use is consistent.
    fn edge_vertex(&mut self, edge: EdgeId, use_id: VertexUseId) -> Option<VertexId> {
        let Some(vertex_use) = self.model.get(use_id) else {
            self.error(edge, ProblemKind::DanglingReference, format!("{edge:?} has a missing vertex use"));
            return None;
        };
        if vertex_use.edge != edge {
            self.error(edge, ProblemKind::InconsistentLink, format!("{use_id:?} names another edge"));
        }
        if self.model.get(vertex_use.vertex).is_none() {
            self.error(edge, ProblemKind::DanglingReference, format!("{edge:?} uses a missing vertex"));
            return None;
        }
        Some(vertex_use.vertex)
    }

    fn edge(&mut self, id: EdgeId) {
        let Some(edge) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        let (start, end, tolerance) = (edge.start, edge.end, edge.tolerance);
        let vertices = [self.edge_vertex(id, start), self.edge_vertex(id, end)];
        let Some(slot) = self.model.get(edge.curve) else { return };
        if slot.edge != id {
            self.error(id, ProblemKind::InconsistentLink, format!("the curve of {id:?} names another edge"));
        }
        if !(slot.first <= slot.last) {
            self.error(id, ProblemKind::Geometry, format!("{id:?} has a reversed range"));
            return;
        }
        // Each end vertex lies on the curve at its end of the range.
        for (vertex, t) in vertices.into_iter().zip([slot.first, slot.last]) {
            let Some(vertex) = vertex else { continue };
            let v = self.model.get(vertex).expect("checked vertex");
            let gap = slot.curve.point(t).distance(v.point);
            if gap > tolerance.max(v.tolerance) {
                self.error(id, ProblemKind::Geometry, format!("{vertex:?} is {gap} away from the curve of {id:?} at {t}"));
            }
        }
    }

    /// The vertices a coedge runs from and to, following its orientation.
    fn coedge_ends(&mut self, id: CoEdgeId) -> Option<(VertexId, VertexId)> {
        let coedge = self.model.get(id)?;
        let edge_id = coedge.edge;
        let reversed = coedge.orientation == Orientation::Reversed;
        let edge = self.model.get(edge_id)?;
        let (start, end) = (edge.start, edge.end);
        let start = self.edge_vertex(edge_id, start)?;
        let end = self.edge_vertex(edge_id, end)?;
        Some(if reversed { (end, start) } else { (start, end) })
    }

    fn coedge(&mut self, id: CoEdgeId) {
        let Some(coedge) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        let (edge, wire, face) = (coedge.edge, coedge.wire, coedge.face);
        if self.model.get(edge).is_none() {
            self.error(id, ProblemKind::DanglingReference, format!("{id:?} uses a missing edge"));
        }
        match self.model.get(wire) {
            Some(w) if w.coedges.contains(&id) => {}
            _ => self.error(id, ProblemKind::InconsistentLink, format!("{id:?} is not in its wire {wire:?}")),
        }
        if face.is_valid() && self.model.get(face).is_none() {
            self.error(id, ProblemKind::DanglingReference, format!("{id:?} names a missing face"));
        }
    }

    fn wire(&mut self, id: WireId) {
        let Some(wire) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        let coedges = wire.coedges.clone();
        let mut ends = Vec::with_capacity(coedges.len());
        for &c in &coedges {
            match self.model.get(c) {
                Some(coedge) if coedge.wire == id => ends.push(self.coedge_ends(c)),
                Some(_) => self.error(id, ProblemKind::InconsistentLink, format!("{c:?} names another wire")),
                None => self.error(id, ProblemKind::DanglingReference, format!("{id:?} has a missing coedge")),
            }
        }
        // On a face the boundary must be closed: each coedge ends where the next starts.
        let on_face = coedges.first().and_then(|&c| self.model.get(c)).is_some_and(|c| c.face.is_valid());
        let n = ends.len();
        let links = if on_face { n } else { n.saturating_sub(1) };
        for i in 0..links {
            if let (Some((_, end)), Some((start, _))) = (ends[i], ends[(i + 1) % n]) {
                if end != start {
                    self.error(
                        id,
                        ProblemKind::DisconnectedWire,
                        format!("coedges {i} and {} of {id:?} do not share a vertex", (i + 1) % n),
                    );
                }
            }
        }
    }

    fn face(&mut self, id: FaceId) {
        let Some(face) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        let wire_uses = face.wires.clone();
        let tolerance = face.tolerance;
        let surface = self.model.get(face.surface).map(|slot| slot.surface.clone());
        if surface.is_none() && self.model.get(face.mesh).is_none() {
            self.error(id, ProblemKind::Geometry, format!("{id:?} has neither a surface nor a mesh"));
        }
        if wire_uses.is_empty() && surface.is_some() {
            self.problem(id, ProblemKind::FaceWithoutBoundary, Severity::Warning, format!("{id:?} has no wires"));
        }
        let mut all = Vec::new();
        for use_id in wire_uses {
            let Some(wire_use) = self.model.get(use_id) else {
                self.error(id, ProblemKind::DanglingReference, format!("{id:?} has a missing wire use"));
                continue;
            };
            if wire_use.face != id {
                self.error(id, ProblemKind::InconsistentLink, format!("{use_id:?} names another face"));
            }
            let Some(wire) = self.model.get(wire_use.wire) else {
                self.error(id, ProblemKind::DanglingReference, format!("{id:?} uses a missing wire"));
                continue;
            };
            let coedges = wire.coedges.clone();
            for &c in &coedges {
                let Some(coedge) = self.model.get(c) else { continue };
                if coedge.face != id {
                    self.error(id, ProblemKind::InconsistentLink, format!("{c:?} does not name its face"));
                }
                let Some(surface) = &surface else { continue };
                let Some((start, _)) = self.coedge_ends(c) else { continue };
                let vertex = self.model.get(start).expect("checked vertex");
                let distance = surface_distance(surface, vertex.point).unwrap_or(0.0);
                if distance > tolerance + vertex.tolerance {
                    self.error(id, ProblemKind::Geometry, format!("{start:?} is {distance} away from the surface of {id:?}"));
                }
                self.same_parameter(id, c, surface, tolerance);
            }
            if let Some(surface) = &surface {
                self.closed_in_parameters(id, &coedges, surface, tolerance);
            }
            all.extend(coedges);
        }
        if let Some(surface) = &surface {
            self.boundary_crossings(id, &all, surface);
        }
    }

    /// The surface at a coedge's pcurve point follows the edge's curve (or, for a
    /// degenerate edge, stays at its vertex) at equal parameters.
    fn same_parameter(&mut self, face: FaceId, coedge: CoEdgeId, surface: &Surface, face_tolerance: f64) {
        let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) else {
            if !matches!(surface, Surface::Plane(_)) {
                self.error(face, ProblemKind::Geometry, format!("{coedge:?} has no pcurve on {face:?}"));
            }
            return;
        };
        let edge_id = self.model.get(coedge).expect("checked coedge").edge;
        let Some(edge) = self.model.get(edge_id) else { return };
        let tolerance = edge.tolerance + face_tolerance;
        let curve = self.model.get(edge.curve);
        let pole = match self.model.get(edge.start) {
            Some(u) => self.model.get(u.vertex).map(|v| v.point),
            None => None,
        };
        if let Some(slot) = curve {
            if (slot.first - first).abs() > PARAMETRIC || (slot.last - last).abs() > PARAMETRIC {
                self.error(face, ProblemKind::Geometry, format!("{coedge:?} and {edge_id:?} have different ranges"));
                return;
            }
        }
        const SAMPLES: usize = 8;
        for i in 0..=SAMPLES {
            let t = first + (last - first) * i as f64 / SAMPLES as f64;
            let uv = pcurve.point(t);
            let on_surface = surface.point(uv.x, uv.y);
            let expected = match (curve, pole) {
                (Some(slot), _) => slot.curve.point(t),
                (None, Some(pole)) => pole,
                (None, None) => return,
            };
            let gap = on_surface.distance(expected);
            if gap > tolerance {
                self.error(
                    face,
                    ProblemKind::Geometry,
                    format!("the pcurve of {coedge:?} is {gap} away from {edge_id:?} at {t}"),
                );
                return;
            }
        }
    }

    /// Consecutive pcurves of a loop meet in the surface's parameter space (so that seams
    /// are used on the correct sides).
    fn closed_in_parameters(&mut self, face: FaceId, coedges: &[CoEdgeId], surface: &Surface, tolerance: f64) {
        let ends: Vec<_> = coedges
            .iter()
            .map(|&c| {
                let (pcurve, first, last) = self.model.effective_pcurve(c)?;
                let reversed = self.model.get(c)?.orientation == Orientation::Reversed;
                let (a, b) = (pcurve.point(first), pcurve.point(last));
                Some(if reversed { (b, a) } else { (a, b) })
            })
            .collect();
        let periods = [surface.u_period(), surface.v_period()];
        for i in 0..ends.len() {
            let (Some((_, end)), Some((start, _))) = (ends[i], ends[(i + 1) % ends.len()]) else { continue };
            // They meet on the surface within the face's tolerance or the tolerance of the
            // vertex where they meet (as OCCT's BRepCheck_Wire), and on the same side of any
            // seam (not a period apart).
            let at = self.coedge_ends(coedges[i]).and_then(|(_, v)| self.model.get(v)).map_or(0.0, |v| v.tolerance);
            let gap = surface.point(end.x, end.y).distance(surface.point(start.x, start.y));
            let across_seam = [(end.x - start.x).abs(), (end.y - start.y).abs()]
                .iter()
                .zip(periods)
                .any(|(d, period)| period.is_some_and(|p| *d > 0.5 * p));
            if gap > 2.0 * tolerance.max(at).max(LINEAR) || across_seam {
                self.error(
                    face,
                    ProblemKind::DisconnectedWire,
                    format!("the pcurves of {:?} and {:?} do not meet on {face:?}", coedges[i], coedges[(i + 1) % ends.len()]),
                );
            }
        }
    }

    /// The face's boundary does not cross itself: no two of its coedges' pcurves (in one
    /// loop or two) cross, nor one itself, other than at a vertex they share (OCCT's
    /// `BRepCheck_SelfIntersectingWire`, `BRepCheck_IntersectingWires`).
    fn boundary_crossings(&mut self, face: FaceId, coedges: &[CoEdgeId], surface: &Surface) {
        if let Some((e1, e2, p)) = boundary_crossing(self.model, coedges, surface) {
            self.error(face, ProblemKind::SelfIntersection, format!("the boundary of {face:?} crosses itself ({e1:?} and {e2:?}) at {p:?}"));
        }
    }

    fn shell(&mut self, id: ShellId) {
        let Some(shell) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        for &use_id in &shell.faces.clone() {
            match self.model.get(use_id) {
                Some(face_use) if face_use.shell != id => {
                    self.error(id, ProblemKind::InconsistentLink, format!("{use_id:?} names another shell"))
                }
                Some(face_use) if self.model.get(face_use.face).is_none() => {
                    self.error(id, ProblemKind::DanglingReference, format!("{id:?} uses a missing face"))
                }
                Some(_) => {}
                None => self.error(id, ProblemKind::DanglingReference, format!("{id:?} has a missing face use")),
            }
        }
    }

    fn solid(&mut self, solid: Shape) {
        let id = SolidId::from_index(solid.node.index);
        let Some(record) = self.model.get(id) else {
            return self.error(id, ProblemKind::DanglingReference, format!("{id:?} is missing"));
        };
        for &use_id in &record.shells.clone() {
            let Some(shell_use) = self.model.get(use_id) else {
                self.error(id, ProblemKind::DanglingReference, format!("{id:?} has a missing shell use"));
                continue;
            };
            if shell_use.solid != id {
                self.error(id, ProblemKind::InconsistentLink, format!("{use_id:?} names another solid"));
            }
            self.closed_shell(id, Shape::from(shell_use.shell));
        }
        if self.report.is_valid() {
            // (Its shells are closed: from its faces' terms, which the model keeps.)
            let forward = Shape { orientation: Orientation::Forward, ..solid };
            let volume = properties::solid_volume(self.model, forward).map_or_else(|| properties::volume(self.model, forward).mass, |(v, _)| v);
            if volume <= 0.0 {
                self.error(id, ProblemKind::Geometry, format!("{id:?} encloses no volume (inside out?): {volume}"));
            }
        }
    }

    /// Every edge of a shell bounding a solid must be used as many times in one direction as
    /// in the other: usually twice, by two faces (or by one face along a seam); four times
    /// where the solid touches itself along the edge (as OCCT's `BRepCheck_Shell` allows).
    fn closed_shell(&mut self, solid: SolidId, shell: Shape) {
        for (edge, (forward, reversed)) in edge_uses(self.model, shell) {
            if forward != reversed {
                self.error(
                    solid,
                    ProblemKind::NotClosed,
                    format!("{edge:?} is used {forward} times forward and {reversed} times reversed in a shell of {solid:?}"),
                );
            }
        }
    }
}

/// The faces of `shape` whose boundaries cross themselves (see [`check`]), each with two
/// edges crossing there (or one edge twice).
pub fn boundary_crossings(model: &Model, shape: Shape) -> Vec<(Shape, [Shape; 2])> {
    let mut out = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face) else { continue };
        let Ok(id) = FaceId::try_from(face.node) else { continue };
        let Some(record) = model.get(id) else { continue };
        let coedges: Vec<CoEdgeId> = record.wires.iter().filter_map(|&u| model.get(u)).filter_map(|u| model.get(u.wire)).flat_map(|w| w.coedges.iter().copied()).collect();
        if let Some((e1, e2, _)) = boundary_crossing(model, &coedges, surface) {
            out.push((face, [Shape::from(NodeId::from(e1)), Shape::from(NodeId::from(e2))]));
        }
    }
    out
}

/// Two edges whose pcurves cross on `surface` (or one crossing itself), and where, among the
/// coedges `coedges` of a face: away from a vertex they share, and not only touching.
fn boundary_crossing(model: &Model, coedges: &[CoEdgeId], surface: &Surface) -> Option<(EdgeId, EdgeId, Point)> {
    // Each pcurve as a polyline; the crossings of their sides refined on the pcurves.
    struct Run {
        edge: EdgeId,
        curve: Curve2,
        params: Vec<f64>,
        points: Vec<Point2>,
    }
    let mut runs = Vec::new();
    for &c in coedges {
        let Some(coedge) = model.get(c) else { continue };
        let edge = coedge.edge;
        let Some((curve, first, last)) = model.effective_pcurve(c) else { continue };
        let n = match &curve {
            Curve2::Line(_) => 1,
            other => (8 * (other.breaks().iter().filter(|&&t| t > first && t < last).count() + 1)).clamp(16, 512),
        };
        let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
        let points = params.iter().map(|&t| curve.point(t)).collect();
        runs.push(Run { edge, curve, params, points });
    }
    let mut sides: Vec<(usize, usize, f64, f64)> = Vec::new();
    let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
    for (r, run) in runs.iter().enumerate() {
        for (i, w) in run.points.windows(2).enumerate() {
            sides.push((r, i, w[0].x.min(w[1].x), w[0].x.max(w[1].x)));
        }
        for p in &run.points {
            lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
        }
    }
    if !(lo.x <= hi.x) {
        return None;
    }
    let size = (hi - lo).length().max(1e-12);
    sides.sort_by(|a, b| a.2.total_cmp(&b.2));
    let cross = |a: Point2, b: Point2, c: Point2| (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    for k in 0..sides.len() {
        let (r1, i1, _, x1) = sides[k];
        for &(r2, i2, x0, _) in &sides[k + 1..] {
            if x0 > x1 {
                break;
            }
            let (run1, run2) = (&runs[r1], &runs[r2]);
            if r1 == r2 && i1.abs_diff(i2) <= 1 || r1 != r2 && run1.edge == run2.edge {
                continue;
            }
            let (a, b) = (run1.points[i1], run1.points[i1 + 1]);
            let (c, d) = (run2.points[i2], run2.points[i2 + 1]);
            let (d1, d2, d3, d4) = (cross(c, d, a), cross(c, d, b), cross(a, b, c), cross(a, b, d));
            if !(d1 * d2 < 0.0 && d3 * d4 < 0.0) {
                continue;
            }
            // Where the pcurves themselves cross, from where the sides do (none where
            // they only come close).
            let (alpha, beta) = (d1 / (d1 - d2), d3 / (d3 - d4));
            let span = |run: &Run, i: usize, f: f64| run.params[i] + f * (run.params[i + 1] - run.params[i]);
            let range = |r: &Run| (r.params[0], *r.params.last().unwrap());
            let Some((s, t)) = pcurves_crossing((&run1.curve, range(run1)), (&run2.curve, range(run2)), span(run1, i1, alpha), span(run2, i2, beta), size) else {
                continue;
            };
            if r1 == r2 && (s - t).abs() <= 1e-6 * (run1.params.last().unwrap() - run1.params[0]).abs() {
                continue;
            }
            // (Not at a vertex the two edges share.)
            let uv = run1.curve.point(s);
            let p = surface.point(uv.x, uv.y);
            let reach = model.get(run1.edge).map_or(0.0, |e| e.tolerance) + model.get(run2.edge).map_or(0.0, |e| e.tolerance);
            let ends = |e: EdgeId| -> Vec<VertexId> {
                model.get(e).map_or(Vec::new(), |e| [e.start, e.end].iter().filter_map(|&u| model.get(u).map(|u| u.vertex)).collect())
            };
            let (ends1, ends2) = (ends(run1.edge), ends(run2.edge));
            let at_vertex = ends1.iter().filter(|v| ends2.contains(v)).filter_map(|&v| model.get(v)).any(|v| v.point.distance(p) <= v.tolerance + reach + LINEAR);
            if at_vertex {
                continue;
            }
            return Some((run1.edge, run2.edge, p));
        }
    }
None
}


/// Where two pcurves cross near the parameters `s` and `t` (Newton's method, within their
/// ranges): none if it does not settle on a crossing (to `size`'s precision), or the
/// pcurves only touch there.
fn pcurves_crossing((c1, (a1, b1)): (&Curve2, (f64, f64)), (c2, (a2, b2)): (&Curve2, (f64, f64)), mut s: f64, mut t: f64, size: f64) -> Option<(f64, f64)> {
    for _ in 0..30 {
        let (p, dp) = c1.d1(s);
        let (q, dq) = c2.d1(t);
        let (fx, fy) = (p.x - q.x, p.y - q.y);
        if fx.hypot(fy) <= 1e-10 * size {
            // (Crossing, not touching: the pcurves running apart there.)
            let sine = (dp.x * dq.y - dp.y * dq.x) / (dp.x.hypot(dp.y) * dq.x.hypot(dq.y));
            return (sine.abs() > 1e-3).then_some((s, t));
        }
        // [dp, -dq] (ds, dt) = -(p - q)
        let det = -dp.x * dq.y + dp.y * dq.x;
        if det.abs() <= 1e-300 {
            return None;
        }
        let ds = (-fx * -dq.y + dq.x * -fy) / det;
        let dt = (dp.x * -fy - dp.y * -fx) / det;
        s = (s + ds).clamp(a1.min(b1), a1.max(b1));
        t = (t + dt).clamp(a2.min(b2), a2.max(b2));
    }
    None
}

/// How many times each (not degenerate) edge of a shell is used forward and reversed.
fn edge_uses(model: &Model, shell: Shape) -> HashMap<NodeId, (usize, usize)> {
    let mut uses: HashMap<NodeId, (usize, usize)> = HashMap::new();
    for edge in model.explore(shell, NodeKind::Edge) {
        if model.is_degenerate(edge) {
            continue;
        }
        let entry = uses.entry(edge.node).or_default();
        match edge.orientation {
            Orientation::Forward => entry.0 += 1,
            Orientation::Reversed => entry.1 += 1,
            _ => {}
        }
    }
    uses
}

/// Checks the consistency of the model's records: everything [`check_model_structure`]
/// checks, plus topological sanity such as coedges whose edges connect along their wire
/// and seam coedge pairs with opposite orientations. Geometry is not checked (see
/// [`check`]).
#[doc(alias = "BRepGraph_Validate")]
pub fn check_model(model: &Model) -> Report {
    let mut report = check_model_structure(model);
    let mut problem = |node: NodeId, kind: ProblemKind, severity: Severity, description: String| {
        report.problems.push(Problem { node, kind, severity, description });
    };
    // Coedges of a wire connect where they follow each other.
    let ends = |c: CoEdgeId| -> Option<(VertexId, VertexId)> {
        let coedge = model.get(c)?;
        let edge = model.get(coedge.edge)?;
        let (a, b) = (model.get(edge.start)?.vertex, model.get(edge.end)?.vertex);
        Some(if coedge.orientation == Orientation::Reversed { (b, a) } else { (a, b) })
    };
    let meets = |u: VertexId, v: VertexId| match (model.get(u), model.get(v)) {
        _ if u == v => true,
        (Some(p), Some(q)) => p.point.distance(q.point) <= p.tolerance + q.tolerance,
        _ => false,
    };
    for (wire, w) in model.table::<crate::topology::Wire>().iter() {
        for pair in w.coedges.windows(2) {
            if let (Some(a), Some(b)) = (ends(pair[0]), ends(pair[1])) {
                if !meets(a.1, b.0) {
                    problem(wire.into(), ProblemKind::DisconnectedWire, Severity::Warning, format!("{:?} does not start where {:?} ends", pair[1], pair[0]));
                }
            }
        }
    }
    // Seam coedges on one face run opposite ways.
    for (face, _) in model.table::<crate::topology::Face>().iter() {
        let mut by_edge: HashMap<EdgeId, Vec<Orientation>> = HashMap::new();
        for c in model.descendants_of::<crate::topology::CoEdge>(face) {
            if let Some(coedge) = model.get(c) {
                by_edge.entry(coedge.edge).or_default().push(coedge.orientation);
            }
        }
        for (edge, uses) in by_edge {
            if uses.len() == 2 && uses[0] == uses[1] {
                problem(face.into(), ProblemKind::InconsistentLink, Severity::Error, format!("the seam {edge:?} is used twice the same way"));
            }
        }
    }
    for (face, f) in model.table::<crate::topology::Face>().iter() {
        if f.wires.is_empty() && model.get(f.surface).is_some() {
            problem(face.into(), ProblemKind::FaceWithoutBoundary, Severity::Warning, format!("{face:?} has a surface but no boundary"));
        }
    }
    // An instance places a product or a topology root, not an instance.
    for (instance, i) in model.table::<crate::topology::Instance>().iter() {
        if i.target.kind == NodeKind::Instance {
            problem(instance.into(), ProblemKind::ProductStructure, Severity::Error, format!("{instance:?} places an instance"));
        }
    }
    // The product graph has no cycle.
    for (product, _) in model.table::<crate::topology::Product>().iter() {
        let mut stack: Vec<crate::topology::ProductId> = vec![product];
        let mut seen = HashSet::new();
        let mut steps = 0;
        while let Some(p) = stack.pop() {
            steps += 1;
            let Some(record) = model.get(p) else { continue };
            for u in record.children.iter().filter_map(|&u| model.get(u)) {
                let Some(instance) = model.get(u.instance) else { continue };
                let Ok(target) = crate::topology::ProductId::try_from(instance.target) else { continue };
                if target == product {
                    problem(product.into(), ProblemKind::ProductStructure, Severity::Error, format!("{product:?} instances itself"));
                    stack.clear();
                    break;
                }
                if seen.insert(target) && steps < 100_000 {
                    stack.push(target);
                }
            }
        }
    }
    report
}

/// The structural part of [`check_model`]: ids referring to removed or missing records,
/// use entries disagreeing with their parent or child, slots disagreeing with their owner.
pub fn check_model_structure(model: &Model) -> Report {
    use crate::topology::{
        CoEdge, CompSolid, Compound, CompoundChild, Edge, Face, FaceUse, Instance, Product, ProductChild, Shell,
        ShellUse, Solid, SolidUse, VertexUse, Wire, WireUse,
    };
    let mut report = Report::default();
    let mut problem = |node: NodeId, kind: ProblemKind, description: String| {
        report.problems.push(Problem { node, kind, severity: Severity::Error, description });
    };
    // A valid id must name a live record.
    macro_rules! live {
        ($owner:expr, $id:expr) => {{
            let id = $id;
            if id.is_valid() && model.get(id).is_none() {
                problem($owner, ProblemKind::DanglingReference, format!("{:?} refers to the missing {:?}", $owner, id));
            }
        }};
    }
    let live_node = |n: NodeId| match n.kind {
        NodeKind::Vertex => model.get(VertexId::from_index(n.index)).is_some(),
        NodeKind::Edge => model.get(EdgeId::from_index(n.index)).is_some(),
        NodeKind::CoEdge => model.get(CoEdgeId::from_index(n.index)).is_some(),
        NodeKind::Wire => model.get(WireId::from_index(n.index)).is_some(),
        NodeKind::Face => model.get(FaceId::from_index(n.index)).is_some(),
        NodeKind::Shell => model.get(ShellId::from_index(n.index)).is_some(),
        NodeKind::Solid => model.get(SolidId::from_index(n.index)).is_some(),
        NodeKind::Compound => model.get(crate::topology::CompoundId::from_index(n.index)).is_some(),
        NodeKind::CompSolid => model.get(crate::topology::CompSolidId::from_index(n.index)).is_some(),
        NodeKind::Product => model.get(crate::topology::ProductId::from_index(n.index)).is_some(),
        NodeKind::Instance => model.get(crate::topology::InstanceId::from_index(n.index)).is_some(),
    };
    for (id, e) in model.table::<Edge>().iter() {
        let owner: NodeId = id.into();
        for u in [e.start, e.end] {
            live!(owner, u);
            if let Some(used) = model.get(u) {
                live!(owner, used.vertex);
                if used.edge != id {
                    problem(owner, ProblemKind::InconsistentLink, format!("{u:?} of {id:?} belongs to {:?}", used.edge));
                }
            }
        }
        if let Some(slot) = model.get(e.curve) {
            if model.get(slot.edge).is_none_or(|o| o.curve != e.curve) {
                problem(owner, ProblemKind::InconsistentLink, format!("the curve of {id:?} names {:?} as its owner", slot.edge));
            }
        }
    }
    for (id, c) in model.table::<CoEdge>().iter() {
        let owner: NodeId = id.into();
        live!(owner, c.edge);
        if !c.edge.is_valid() {
            problem(owner, ProblemKind::DanglingReference, format!("{id:?} has no edge"));
        }
        live!(owner, c.wire);
        live!(owner, c.face);
        if let Some(w) = model.get(c.wire) {
            if !w.coedges.contains(&id) {
                problem(owner, ProblemKind::InconsistentLink, format!("{id:?} is not in its {:?}", c.wire));
            }
        }
        if let Some(slot) = model.get(c.pcurve) {
            if slot.coedge != id {
                problem(owner, ProblemKind::InconsistentLink, format!("the pcurve of {id:?} names {:?}", slot.coedge));
            }
        }
    }
    for (id, w) in model.table::<Wire>().iter() {
        for &c in &w.coedges {
            live!(NodeId::from(id), c);
            if model.get(c).is_some_and(|c| c.wire != id) {
                problem(id.into(), ProblemKind::InconsistentLink, format!("{c:?} of {id:?} names another wire"));
            }
        }
    }
    for (id, f) in model.table::<Face>().iter() {
        for &u in &f.wires {
            live!(NodeId::from(id), u);
            if let Some(used) = model.get(u) {
                live!(NodeId::from(id), used.wire);
                if used.face != id {
                    problem(id.into(), ProblemKind::InconsistentLink, format!("{u:?} of {id:?} names {:?}", used.face));
                }
            }
        }
        if let Some(slot) = model.get(f.surface) {
            if model.get(slot.face).is_none_or(|o| o.surface != f.surface) {
                problem(id.into(), ProblemKind::InconsistentLink, format!("the surface of {id:?} names {:?} as its owner", slot.face));
            }
        }
    }
    macro_rules! container {
        ($ty:ty, $list:ident, $parent:ident, $child:ident) => {
            for (id, record) in model.table::<$ty>().iter() {
                for &u in &record.$list {
                    live!(NodeId::from(id), u);
                    if let Some(used) = model.get(u) {
                        live!(NodeId::from(id), used.$child);
                        if used.$parent != id {
                            problem(id.into(), ProblemKind::InconsistentLink, format!("{u:?} of {id:?} names another parent"));
                        }
                    }
                }
            }
        };
    }
    container!(Shell, faces, shell, face);
    container!(Solid, shells, solid, shell);
    container!(CompSolid, solids, comp_solid, solid);
    container!(Product, children, product, instance);
    for (id, c) in model.table::<Compound>().iter() {
        for &u in &c.children {
            live!(NodeId::from(id), u);
            if let Some(used) = model.get(u) {
                if used.compound != id {
                    problem(id.into(), ProblemKind::InconsistentLink, format!("{u:?} of {id:?} names another compound"));
                }
                if !live_node(used.child) {
                    problem(id.into(), ProblemKind::DanglingReference, format!("{id:?} holds the missing {:?}", used.child));
                }
            }
        }
    }
    for (id, p) in model.table::<Product>().iter() {
        if let Some(shape) = p.shape {
            if !live_node(shape.node) {
                problem(id.into(), ProblemKind::DanglingReference, format!("the shape of {id:?} is missing"));
            }
        }
    }
    for (id, i) in model.table::<Instance>().iter() {
        if !live_node(i.target) {
            problem(id.into(), ProblemKind::DanglingReference, format!("{id:?} places the missing {:?}", i.target));
        }
    }
    // Use records their parents do not list.
    macro_rules! listed {
        ($ty:ty, $parent_ty:ty, $parent:ident, $list:ident) => {
            for (u, used) in model.table::<$ty>().iter() {
                if model.get(used.$parent).is_some_and(|p| !p.$list.contains(&u)) {
                    problem(used.$parent.into(), ProblemKind::InconsistentLink, format!("{u:?} is not listed by its parent"));
                }
            }
        };
    }
    listed!(WireUse, Face, face, wires);
    listed!(FaceUse, Shell, shell, faces);
    listed!(ShellUse, Solid, solid, shells);
    listed!(SolidUse, CompSolid, comp_solid, solids);
    listed!(CompoundChild, Compound, compound, children);
    listed!(ProductChild, Product, product, children);
    for (u, used) in model.table::<VertexUse>().iter() {
        if model.get(used.edge).is_some_and(|e| e.start != u && e.end != u) {
            problem(used.edge.into(), ProblemKind::InconsistentLink, format!("{u:?} is not an end of {:?}", used.edge));
        }
    }
    // Slots no record uses.
    let curves: HashSet<_> = model.table::<Edge>().iter().map(|(_, e)| e.curve).collect();
    for (slot, s) in model.table::<crate::topology::EdgeCurve>().iter() {
        if !curves.contains(&slot) {
            problem(s.edge.into(), ProblemKind::OrphanSlot, format!("{slot:?} is used by no edge"));
        }
    }
    let surfaces: HashSet<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    for (slot, s) in model.table::<crate::topology::FaceSurface>().iter() {
        if !surfaces.contains(&slot) {
            problem(s.face.into(), ProblemKind::OrphanSlot, format!("{slot:?} is used by no face"));
        }
    }
    let pcurves: HashSet<_> = model.table::<CoEdge>().iter().map(|(_, c)| c.pcurve).collect();
    for (slot, s) in model.table::<crate::topology::PCurve>().iter() {
        if !pcurves.contains(&slot) {
            problem(s.coedge.into(), ProblemKind::OrphanSlot, format!("{slot:?} is used by no coedge"));
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::Plane;
    use crate::math::Point;
    use crate::primitives::Cuboid;
    use crate::tolerance::LINEAR;

    /// A face on the XY plane bounded by the polygon through `points` (straight edges).
    fn polygon_face(model: &mut Model, points: &[Point]) -> Shape {
        let vertices: Vec<_> = points.iter().map(|&p| model.add_vertex(p, LINEAR)).collect();
        let coedges: Vec<_> = (0..vertices.len())
            .map(|i| {
                let edge = model.add_edge(vertices[i], vertices[(i + 1) % vertices.len()], LINEAR).unwrap();
                model.add_coedge(edge, Orientation::Forward).unwrap()
            })
            .collect();
        let wire = model.add_wire(&coedges).unwrap();
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, Surface::Plane(Plane::XY), LINEAR).unwrap();
        model.add(face, wire.into()).unwrap();
        face
    }

    #[test]
    fn box_is_valid() {
        let mut model = Model::new();
        let solid = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let report = check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        assert!(report.problems.is_empty(), "{report:?}");
    }

    #[test]
    fn curved_primitives_are_valid() {
        use crate::primitives::{Cone, Cylinder, Sphere};
        use std::f64::consts::PI;
        let mut model = Model::new();
        for shape in [
            Cylinder::new(2.0, 5.0).build(&mut model).unwrap(),
            Cylinder::new(2.0, 5.0).with_angle(1.0).build(&mut model).unwrap(),
            Cylinder::new(2.0, 5.0).with_angle(1.5 * PI).build(&mut model).unwrap(),
            Sphere::new(3.0).build(&mut model).unwrap(),
            Sphere::new(3.0).with_angle(2.0).build(&mut model).unwrap(),
            Sphere::new(3.0).with_latitudes(-0.3, 1.1).build(&mut model).unwrap(),
            Sphere::new(3.0).with_latitudes(0.0, PI / 2.0).with_angle(PI / 2.0).build(&mut model).unwrap(),
            Sphere::new(3.0).with_latitudes(-1.0, 0.5).with_angle(4.0).build(&mut model).unwrap(),
            Cone::new(3.0, 1.0, 4.0).build(&mut model).unwrap(),
            Cone::new(1.0, 3.0, 4.0).build(&mut model).unwrap(),
            Cone::new(3.0, 0.0, 4.0).build(&mut model).unwrap(),
            Cone::new(0.0, 3.0, 4.0).build(&mut model).unwrap(),
            Cone::new(3.0, 1.0, 4.0).with_angle(2.5).build(&mut model).unwrap(),
            Cone::new(3.0, 0.0, 4.0).with_angle(PI).build(&mut model).unwrap(),
            Cone::new(0.0, 3.0, 4.0).with_angle(5.0).build(&mut model).unwrap(),
        ] {
            let report = check(&model, shape);
            assert!(report.problems.is_empty(), "{report:?}");
        }
    }

    /// Swapping the pcurves of a cylinder's seam puts each half on the wrong side.
    #[test]
    fn swapped_seam_pcurves_are_reported() {
        use crate::primitives::Cylinder;
        let mut model = Model::new();
        let solid = Cylinder::new(2.0, 5.0).build(&mut model).unwrap();
        let seam = model.table::<crate::topology::CoEdge>().iter().map(|(id, _)| id).find(|&c| model.is_seam(c)).unwrap();
        let mate = model.seam_pair(seam).unwrap();
        let (a, b) = (model.pcurve(seam).unwrap().clone(), model.pcurve(mate).unwrap().clone());
        model.set_pcurve(seam, b).unwrap();
        model.set_pcurve(mate, a).unwrap();
        let report = check(&model, solid);
        assert!(report.has(ProblemKind::DisconnectedWire, Severity::Error), "{report:?}");
    }

    /// An edge curve that misses the edge's vertices and the face's pcurve.
    #[test]
    fn edge_curve_off_its_vertices_is_reported() {
        use crate::geometry::curve::{Circle, Curve};
        use crate::math::{Direction, Frame};
        use crate::primitives::Cylinder;
        let mut model = Model::new();
        let solid = Cylinder::new(2.0, 5.0).build(&mut model).unwrap();
        let circle = model
            .explore(solid, NodeKind::Edge)
            .find(|&e| matches!(model.get(model.get(crate::topology::EdgeId::from_index(e.node.index)).unwrap().curve).map(|s| &s.curve), Some(Curve::Circle(_))))
            .unwrap();
        let bigger = Curve::Circle(Circle::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X), 2.5));
        model.set_edge_curve(circle, bigger, LINEAR).unwrap();
        model.set_edge_range(circle, 0.0, std::f64::consts::TAU).unwrap();
        let report = check(&model, solid);
        assert!(report.has(ProblemKind::Geometry, Severity::Error), "{report:?}");
    }

    #[test]
    fn invalid_cones_and_zones_are_rejected() {
        use crate::primitives::{Cone, Sphere};
        let mut model = Model::new();
        assert!(Cone::new(2.0, 2.0, 3.0).build(&mut model).is_err(), "equal radii");
        assert!(Cone::new(0.0, 0.0, 3.0).build(&mut model).is_err(), "no radius");
        assert!(Cone::new(2.0, 1.0, 0.0).build(&mut model).is_err(), "no height");
        assert!(Sphere::new(1.0).with_latitudes(0.5, 0.2).build(&mut model).is_err(), "reversed latitudes");
    }

    #[test]
    fn open_shell_is_reported() {
        let mut model = Model::new();
        let solid = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
        let shell = model.new_empty(NodeKind::Shell);
        for &face in &faces[..5] {
            model.add(shell, face).unwrap();
        }
        let open = model.new_empty(NodeKind::Solid);
        model.add(open, shell).unwrap();
        let report = check(&model, open);
        assert!(report.has(ProblemKind::NotClosed, Severity::Error), "{report:?}");
    }

    #[test]
    fn inside_out_solid_is_reported() {
        let mut model = Model::new();
        let solid = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
        let shell = model.new_empty(NodeKind::Shell);
        for face in faces {
            model.add(shell, face.reversed()).unwrap();
        }
        let inverted = model.new_empty(NodeKind::Solid);
        model.add(inverted, shell).unwrap();
        let report = check(&model, inverted);
        assert!(report.has(ProblemKind::Geometry, Severity::Error), "{report:?}");
    }

    #[test]
    fn disconnected_wire_is_reported() {
        let mut model = Model::new();
        let p = |x, y| Point::new(x, y, 0.0);
        let v: Vec<_> = [p(0., 0.), p(1., 0.), p(5., 5.), p(6., 5.)].iter().map(|&q| model.add_vertex(q, LINEAR)).collect();
        let e1 = model.add_edge(v[0], v[1], LINEAR).unwrap();
        let e2 = model.add_edge(v[2], v[3], LINEAR).unwrap();
        let c1 = model.add_coedge(e1, Orientation::Forward).unwrap();
        let c2 = model.add_coedge(e2, Orientation::Forward).unwrap();
        let wire = model.add_wire(&[c1, c2]).unwrap();
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, Surface::Plane(Plane::XY), LINEAR).unwrap();
        model.add(face, wire.into()).unwrap();
        let report = check(&model, face);
        assert!(report.has(ProblemKind::DisconnectedWire, Severity::Error), "{report:?}");
    }

    /// A bow tie: its boundary crosses itself in the middle.
    #[test]
    fn self_crossing_boundary_is_reported() {
        let mut model = Model::new();
        let p = |x, y| Point::new(x, y, 0.0);
        let mut polygon = |points: &[Point]| {
            let vertices: Vec<Shape> = points.iter().map(|&q| crate::make::vertex(&mut model, q)).collect();
            let edges: Vec<Shape> = (0..vertices.len()).map(|i| crate::make::edge(&mut model, vertices[i], vertices[(i + 1) % vertices.len()]).unwrap()).collect();
            let wire = crate::make::wire(&mut model, &edges).unwrap();
            crate::make::face_on(&mut model, Surface::Plane(Plane::XY), wire).unwrap()
        };
        let square = polygon(&[p(0., 0.), p(1., 0.), p(1., 1.), p(0., 1.)]);
        let bow_tie = polygon(&[p(0., 0.), p(1., 1.), p(1., 0.), p(0., 1.)]);
        assert!(check(&model, square).problems.is_empty());
        let report = check(&model, bow_tie);
        assert!(report.has(ProblemKind::SelfIntersection, Severity::Error), "{report:?}");
    }

    #[test]
    fn vertex_off_its_plane_is_reported() {
        let mut model = Model::new();
        let p = Point::new;
        let flat = polygon_face(&mut model, &[p(0., 0., 0.), p(1., 0., 0.), p(1., 1., 0.)]);
        assert!(check(&model, flat).is_valid());
        let lifted = polygon_face(&mut model, &[p(0., 0., 0.), p(1., 0., 0.), p(1., 1., 0.5)]);
        let report = check(&model, lifted);
        assert!(report.has(ProblemKind::Geometry, Severity::Error), "{report:?}");
    }
}
