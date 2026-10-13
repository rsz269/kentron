//! Faces far from the other operand: their boxes clear of where the operands' boxes
//! overlap, so that nothing of the other operand meets them and they lie outside it. They
//! stay out of the boolean graph: they stand only for counting ray crossings when the other
//! operand's points are classified, and go into the result whole, as they are, or not at
//! all. What they need is found once per face and kept in the model with it.

use crate::hashing::FastSet as HashSet;
use std::sync::Arc;

use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::geometry::ParametricSurface;
use crate::math::{BoundingBox, Point};
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape};

use super::geom::{self, on_tol};
use super::graph::{self, FaceSample, Location};

/// One of a face's boundary edges, for telling points on it.
pub(super) struct BoundaryEdge {
    pub edge: EdgeId,
    /// `None` for a degenerate edge.
    pub curve: Option<Curve>,
    pub first: f64,
    pub last: f64,
    /// The edge's start (all of a degenerate edge).
    pub start: Point,
    pub bbox: (Point, Point),
    /// Whether the face runs along the edge's curve.
    pub along: bool,
}

/// What a face far from the other operand is, from the model alone (unplaced).
pub(super) struct FaceData {
    /// The surface computed with (a torus's as a surface of revolution).
    pub surface: Surface,
    /// The face's own surface.
    pub output: Surface,
    pub sample: FaceSample,
    pub periods: [Option<f64>; 2],
    pub edges: Vec<BoundaryEdge>,
    /// A box sure to hold the face; `None` on surfaces without one.
    pub hull: Option<BoundingBox>,
}

/// A face of an operand left out of the graph, far from the other operand.
pub(super) struct FarFace {
    pub operand: usize,
    pub face: FaceId,
    /// +1 if `∂u × ∂v` points out of the operand, -1 if into it.
    pub outward: f64,
    pub data: Arc<FaceData>,
    /// Its place among the operands' faces.
    pub seq: usize,
}

impl FaceData {
    /// Where `p` (on the face's surface) lies with respect to the face: on its boundary
    /// within [`on_tol`], else inside or outside.
    pub fn locate(&self, p: Point) -> Location {
        for e in &self.edges {
            if !geom::box_contains(&e.bbox, p) {
                continue;
            }
            let d = match &e.curve {
                Some(c) => geom::project(c, e.first, e.last, p).1,
                None => e.start.distance(p),
            };
            if d <= on_tol() {
                return Location::Boundary;
            }
        }
        let sample = &self.sample;
        // (Beyond the surface's edge, as in Face::locate.)
        let reach = 1e-3 * sample.extent.0.distance(sample.extent.1) + on_tol();
        let Some(uv) = geom::uv_within(&self.surface, p, reach) else { return Location::Outside };
        let inside = graph::representatives(uv, self.periods, sample.u_range, sample.v_range)
            .into_iter()
            .any(|uv| geom::inside_polylines(&sample.polylines, uv));
        if inside { Location::Inside } else { Location::Outside }
    }
}

/// The face's data, kept in the model with it; `None` for a face the graph must handle
/// itself (placed, without a surface or pcurves, or a planar B-spline it recomputes).
pub(super) fn face_data(model: &Model, face: Shape) -> Option<Arc<FaceData>> {
    if face.placement != crate::math::Transform::IDENTITY {
        return None;
    }
    let id = FaceId::try_from(face.node).ok()?;
    let original = model.surface(face)?;
    if graph::planar_bspline(original).is_some() {
        return None;
    }
    let kept = model.derived(face.node, model.face_fingerprint(id), || Kept(compute(model, Shape::from(id), original.clone())));
    kept.0.clone()
}

/// A face's data as the model keeps it.
struct Kept(Option<Arc<FaceData>>);

fn compute(model: &Model, face: Shape, output: Surface) -> Option<Arc<FaceData>> {
    let surface = match &output {
        Surface::Torus(t) => graph::torus_as_revolution(t),
        s => s.clone(),
    };
    let mut pcurves: Vec<(Curve2, f64, f64)> = Vec::new();
    let mut edges = Vec::new();
    for wire in model.children(face) {
        for coedge in model.children(wire) {
            let coedge_id = CoEdgeId::try_from(coedge.node).ok()?;
            let edge: EdgeId = model[coedge_id].edge;
            let (pcurve, pfirst, plast) = model.effective_pcurve(coedge_id)?;
            let edge_shape = Shape { placement: coedge.placement, ..Shape::from(edge) };
            let (curve, first, last) = match model.edge_curve(edge_shape) {
                Some((c, a, b)) => (Some(c), a, b),
                None => (None, pfirst, plast),
            };
            let (start, _) = model.edge_vertices(edge_shape);
            let start = model.point(start);
            let bbox = graph::Edge::sampled_box(curve.as_ref(), first, last, start);
            let along = model[coedge_id].orientation.compose(wire.orientation) == Orientation::Forward;
            pcurves.push((pcurve, first, last));
            edges.push(BoundaryEdge { edge, curve, first, last, start, bbox, along });
        }
    }
    let ranged: Vec<(&Curve2, f64, f64)> = pcurves.iter().map(|(p, a, b)| (p, *a, *b)).collect();
    let sample = graph::sample_face(&surface, &ranged);
    let hull = match &surface {
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_) => {
            Some(BoundingBox::new(sample.extent.0, sample.extent.1))
        }
        // A torus computed as the revolution of its meridian circle.
        Surface::Revolution(r) if matches!(r.basis(), Curve::Circle(_)) => Some(BoundingBox::new(sample.extent.0, sample.extent.1)),
        Surface::BSpline(_) | Surface::Bezier(_) => crate::properties::control_box(&surface),
        _ => None,
    };
    let periods = [surface.u_period(), surface.v_period()];
    Some(Arc::new(FaceData { surface, output, sample, periods, edges, hull }))
}

/// Finds, for all of the operands' unplaced faces and edges at once, what the graph keeps in
/// the model for them (faces' data, edges' boxes), so that building the graph finds it
/// there: for operands of many faces (for a few, going through them twice costs more than
/// it saves).
/// Each face's data found by [`prepare`] (unplaced faces, by node): what [`choose`] and the
/// graph look up again.
pub(super) type Prepared = crate::hashing::FastMap<NodeId, Option<Arc<FaceData>>>;

/// [`face_data`], from what [`prepare`] found where it has the face.
fn face_data_in(model: &Model, face: Shape, prepared: &Prepared) -> Option<Arc<FaceData>> {
    match prepared.get(&face.node) {
        Some(data) if face.placement == crate::math::Transform::IDENTITY => data.clone(),
        _ => face_data(model, face),
    }
}

pub(super) fn prepare(model: &Model, shapes: [Shape; 2], shape_edges: [&[Shape]; 2]) -> Prepared {
    let mut faces: Vec<FaceId> = Vec::new();
    let mut edges: Vec<EdgeId> = Vec::new();
    for (shape, shape_edges) in shapes.into_iter().zip(shape_edges) {
        for face in model.unique_subshapes(shape, NodeKind::Face) {
            if let (Ok(id), true) = (FaceId::try_from(face.node), face.placement == crate::math::Transform::IDENTITY) {
                faces.push(id);
            }
        }
        for &edge in shape_edges {
            if let (Ok(id), true) = (EdgeId::try_from(edge.node), edge.placement == crate::math::Transform::IDENTITY) {
                edges.push(id);
            }
        }
    }
    if faces.len() < 64 {
        return Prepared::default();
    }
    faces.sort_unstable_by_key(|f| f.index());
    faces.dedup();
    edges.sort_unstable_by_key(|e| e.index());
    edges.dedup();
    // The data, and its samples kept as the graph keeps a face's (the same: the same
    // surface, boundary and ranges).
    let found = crate::parallel::map(&faces, |&id| {
        let data = face_data(model, Shape::from(id));
        if let Some(data) = &data {
            model.derived(id.into(), model.face_fingerprint(id), || data.sample.clone());
        }
        data
    });
    crate::parallel::map(&edges, |&id| graph::warm_edge_box(model, id));
    faces.into_iter().map(NodeId::from).zip(found).collect()
}

/// Faces far from the other operand (unplaced, by node), with their data.
pub(super) type FarData = crate::hashing::FastMap<NodeId, Arc<FaceData>>;

/// The faces of `a` and of `b` (unplaced, by node) far from the other operand: only between
/// operands made of solids whose every face has a sure box.
pub(super) fn choose(model: &Model, a: Shape, b: Shape, prepared: &Prepared) -> [FarData; 2] {
    let none = || [FarData::default(), FarData::default()];
    if !graph::made_of_solids(model, a) || !graph::made_of_solids(model, b) {
        return none();
    }
    let mut hulls: Vec<Vec<(NodeId, BoundingBox)>> = Vec::new();
    let mut edges: Vec<HashSet<EdgeId>> = Vec::new();
    let mut boundaries: Vec<Vec<Arc<FaceData>>> = Vec::new();
    for shape in [a, b] {
        let (mut list, mut own, mut data) = (Vec::new(), HashSet::default(), Vec::new());
        for face in model.explore(shape, NodeKind::Face) {
            let Some(d) = face_data_in(model, face, prepared) else { return none() };
            let Some(hull) = d.hull else { return none() };
            own.extend(d.edges.iter().map(|e| e.edge));
            list.push((face.node, hull));
            data.push(d);
        }
        hulls.push(list);
        edges.push(own);
        boundaries.push(data);
    }
    let joined = |list: &[(NodeId, BoundingBox)]| list.iter().fold(BoundingBox::EMPTY, |mut b, (_, h)| {
        b.add_box(h);
        b
    });
    let (ba, bb) = (joined(&hulls[0]), joined(&hulls[1]));
    if ba.is_empty() || bb.is_empty() {
        return none();
    }
    // Where the boxes overlap (none: every face is far), widened by as much as the
    // operation lets things meet across (free-form sections are fitted to a thousandth of
    // their size).
    let margin = 1e-3 * ba.min.distance(ba.max).max(bb.min.distance(bb.max)) + on_tol();
    let (lo, hi) = (
        Point::new(ba.min.x.max(bb.min.x) - margin, ba.min.y.max(bb.min.y) - margin, ba.min.z.max(bb.min.z) - margin),
        Point::new(ba.max.x.min(bb.max.x) + margin, ba.max.y.min(bb.max.y) + margin, ba.max.z.min(bb.max.z) + margin),
    );
    let common = (lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z).then(|| BoundingBox::new(lo, hi));
    // (Operands sharing nodes, as results of earlier operations do, meet where they share an
    // edge.)
    let far = |k: usize| {
        hulls[k]
            .iter()
            .zip(&boundaries[k])
            .filter(|((_, h), d)| common.is_none_or(|c| !c.overlaps(h)) && !d.edges.iter().any(|e| edges[1 - k].contains(&e.edge)))
            .map(|((n, _), d)| (*n, d.clone()))
            .collect()
    };
    [far(0), far(1)]
}
