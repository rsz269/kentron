//! Faces far from the other operand: their boxes clear of where the operands' boxes
//! overlap, so that nothing of the other operand meets them and they lie outside it; or,
//! within that overlap, clear of every box of the other operand's faces, so that nothing of
//! it meets them either and they lie wholly inside or outside it (found once the graph is
//! built, alike for faces joined by edges). They stay out of the boolean graph: they stand
//! only for counting ray crossings when the other operand's points are classified, and go
//! into the result whole, as they are, or not at all. What they need is found once per face
//! and kept in the model with it.

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
    /// Where the face is placed (its data is unplaced).
    pub placement: crate::math::Transform,
    /// +1 if `∂u × ∂v` points out of the operand, -1 if into it.
    pub outward: f64,
    pub data: Arc<FaceData>,
    /// Its place among the operands' faces.
    pub seq: usize,
    /// True for a face within the operands' overlap, apart from the other operand's faces:
    /// inside or outside it, to be found ([`states`]); else outside it.
    pub apart: bool,
}

impl FarFace {
    /// `b` (a box of the unplaced face) placed with the face.
    pub fn placed(&self, b: (Point, Point)) -> (Point, Point) {
        if self.placement == crate::math::Transform::IDENTITY {
            return b;
        }
        use crate::math::Transformable;
        let placed = BoundingBox::from_points((0..8).map(|i| {
            Point::new(if i & 1 == 0 { b.0.x } else { b.1.x }, if i & 2 == 0 { b.0.y } else { b.1.y }, if i & 4 == 0 { b.0.z } else { b.1.z }).transformed(&self.placement)
        }));
        (placed.min, placed.max)
    }
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
    let fingerprint = model.face_fingerprint(id);
    let kept = model.derived(face.node, fingerprint, || Kept(compute(model, Shape::from(id), original.clone())));
    // (Its samples kept as the graph keeps a face's, the same: the same surface, boundary
    // and ranges.)
    if let Some(data) = &kept.0 {
        model.derived(face.node, fingerprint, || data.sample.clone());
    }
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
    let mut wires = Vec::new();
    for wire in model.children(face) {
        let wire_start = pcurves.len();
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
        wires.push(wire_start..pcurves.len());
    }
    // (Aligned by periods as the graph aligns a face's pcurves.)
    let mut aligned: Vec<Vec<(Curve2, f64, f64, bool)>> = wires.iter().map(|r| r.clone().map(|i| (pcurves[i].0.clone(), pcurves[i].1, pcurves[i].2, edges[i].along)).collect()).collect();
    graph::align_wires(&mut aligned, [surface.u_period(), surface.v_period()]);
    for (r, wire) in wires.iter().zip(aligned) {
        for (i, (pcurve, _, _, _)) in r.clone().zip(wire) {
            pcurves[i].0 = pcurve;
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
    // The data (and its samples, kept for the graph).
    let found = crate::parallel::map(&faces, |&id| face_data(model, Shape::from(id)));
    crate::parallel::map(&edges, |&id| graph::warm_edge_box(model, id));
    faces.into_iter().map(NodeId::from).zip(found).collect()
}

/// Faces far from the other operand (unplaced, by node), with their data and whether they
/// lie within the operands' overlap, only apart from the other's faces.
pub(super) type FarData = crate::hashing::FastMap<NodeId, (Arc<FaceData>, bool)>;

/// The faces of `a` and of `b` (unplaced, by node) far from the other operand: only between
/// operands made of solids whose every face has a sure box.
pub(super) fn choose(model: &Model, a: Shape, b: Shape, prepared: &Prepared, placed_ok: [bool; 2]) -> [FarData; 2] {
    let none = || [FarData::default(), FarData::default()];
    if !graph::made_of_solids(model, a) || !graph::made_of_solids(model, b) {
        return none();
    }
    // Each face's box (a placed face's, its unplaced one's placed: placed faces are left out
    // only where `placed_ok`, but all bound where the operands meet).
    let mut hulls: Vec<Vec<(NodeId, BoundingBox)>> = Vec::new();
    let mut placed: Vec<HashSet<NodeId>> = Vec::new();
    // (A face used more than once is not left out: far faces are known by node.)
    let mut repeated: Vec<HashSet<NodeId>> = Vec::new();
    let mut edges: Vec<HashSet<EdgeId>> = Vec::new();
    let mut boundaries: Vec<Vec<Arc<FaceData>>> = Vec::new();
    for shape in [a, b] {
        let (mut list, mut moved, mut own, mut data) = (Vec::new(), HashSet::default(), HashSet::default(), Vec::new());
        let (mut seen, mut again) = (HashSet::default(), HashSet::default());
        for face in model.explore(shape, NodeKind::Face) {
            if !seen.insert(face.node) {
                again.insert(face.node);
            }
            let unplaced = Shape { placement: crate::math::Transform::IDENTITY, ..face };
            let Some(d) = face_data_in(model, unplaced, prepared) else { return none() };
            let Some(hull) = d.hull else { return none() };
            let hull = if face.placement == crate::math::Transform::IDENTITY {
                hull
            } else {
                moved.insert(face.node);
                BoundingBox::from_points((0..8).map(|i| {
                    use crate::math::Transformable;
                    Point::new(
                        if i & 1 == 0 { hull.min.x } else { hull.max.x },
                        if i & 2 == 0 { hull.min.y } else { hull.max.y },
                        if i & 4 == 0 { hull.min.z } else { hull.max.z },
                    )
                    .transformed(&face.placement)
                }))
            };
            own.extend(d.edges.iter().map(|e| e.edge));
            list.push((face.node, hull));
            data.push(d);
        }
        hulls.push(list);
        placed.push(moved);
        repeated.push(again);
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
    // (Within the overlap, the boxes of the other operand's faces there: a face of one
    // overlapping a face of the other overlaps it within the overlap.)
    let widened = |h: &BoundingBox| BoundingBox::new(
        Point::new(h.min.x - margin, h.min.y - margin, h.min.z - margin),
        Point::new(h.max.x + margin, h.max.y + margin, h.max.z + margin),
    );
    let near: Vec<Vec<BoundingBox>> = hulls.iter().map(|list| list.iter().filter(|(_, h)| common.is_some_and(|c| c.overlaps(h))).map(|(_, h)| *h).collect()).collect();
    let trees: Vec<crate::math::Bvh> = near.iter().map(|boxes| crate::math::Bvh::new(boxes)).collect();
    let far = |k: usize| {
        hulls[k]
            .iter()
            .zip(&boundaries[k])
            .filter(|((n, _), d)| !repeated[k].contains(n) && !d.edges.iter().any(|e| edges[1 - k].contains(&e.edge)))
            .filter_map(|((n, h), d)| {
                if common.is_none_or(|c| !c.overlaps(h)) {
                    return (placed_ok[k] || !placed[k].contains(n)).then(|| (*n, (d.clone(), false)));
                }
                // (Unplaced only: it may be kept, inside the other operand as well as outside.)
                (!placed[k].contains(n) && trees[1 - k].overlapping(&widened(h)).is_empty()).then(|| (*n, (d.clone(), true)))
            })
            .collect()
    };
    [far(0), far(1)]
}

/// Where each far face of `graph` lies with respect to the other operand: outside it, but
/// for those apart from its faces within the operands' overlap, found by a few of each set
/// of them joined by edges (they lie alike). `None` if those disagree or one is found on the
/// other operand (it should not be: the graph is then built without far faces).
pub(super) fn states(graph: &graph::Graph) -> Option<Vec<super::classify::State>> {
    use super::classify::State;
    let mut states = vec![State::Outside; graph.far.len()];
    let apart: Vec<usize> = (0..graph.far.len()).filter(|&i| graph.far[i].apart).collect();
    if apart.is_empty() {
        return Some(states);
    }
    // Sets joined by edges (of one operand).
    let mut parent: Vec<usize> = (0..graph.far.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut by_edge: crate::hashing::FastMap<(usize, EdgeId), usize> = Default::default();
    for &i in &apart {
        let f = &graph.far[i];
        for e in &f.data.edges {
            if let Some(&j) = by_edge.get(&(f.operand, e.edge)) {
                let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                parent[a] = b;
            } else {
                by_edge.insert((f.operand, e.edge), i);
            }
        }
    }
    let mut sets: crate::hashing::FastMap<usize, Vec<usize>> = Default::default();
    for &i in &apart {
        let r = root(&mut parent, i);
        sets.entry(r).or_default().push(i);
    }
    let sets: Vec<Vec<usize>> = sets.into_values().collect();
    let classify = |i: usize| {
        let f = &graph.far[i];
        let polylines = &f.data.sample.polylines;
        let own: Vec<&[geom::Uv]> = polylines.iter().map(|p| p.as_slice()).collect();
        let uv = geom::interior_point_among(polylines, &geom::Obstacles::new(&own));
        let p = geom::surface_point(&f.data.surface, uv);
        let outward = geom::normal_at(&f.data.surface, uv) * f.outward;
        graph.classify(p, outward, 1 - f.operand, 0.0)
    };
    const SAMPLES: usize = 3;
    let found = crate::parallel::map(&sets, |set| {
        let sampled: Vec<usize> = if set.len() <= SAMPLES { set.clone() } else { (0..SAMPLES).map(|k| set[k * set.len() / SAMPLES]).collect() };
        let first = classify(sampled[0]);
        (matches!(first, State::Inside | State::Outside) && sampled[1..].iter().all(|&i| classify(i) == first)).then_some(first)
    });
    for (set, state) in sets.iter().zip(found) {
        let state = state?;
        for &i in set {
            states[i] = state;
        }
    }
    Some(states)
}
