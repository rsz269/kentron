//! Boolean operations and sections.
//!
//! Implemented for shapes bounded by planes, quadrics and surfaces of revolution and
//! extrusion; other surfaces are `todo!()`. Operations on solids are described below; those on
//! other shapes (faces, edges, vertices) in [`mixed`]. The algorithm follows OCCT's general
//! fuse:
//!
//! 1. Every edge of one operand is intersected with every face of the other, splitting the
//!    edge where it meets the face (or noting that it lies on the face's surface). Edges
//!    are also split at any vertex lying on them, so both operands agree on their meeting
//!    points.
//! 2. Every face of one operand is intersected with every face of the other: the curves
//!    where their surfaces meet ([`intersect::quadrics`](crate::geometry::intersect::quadrics),
//!    or branches of [`intersect::cylinder_quadric`](crate::geometry::intersect::cylinder_quadric)
//!    for crossing cylinders) are cut at the points found in step 1, and the stretches
//!    inside both faces become section edges.
//! 3. Each face is split into regions along the pieces of edges lying on it, traced in
//!    the face's parameter space.
//! 4. Each region is classified against the other operand: inside, outside, or on its
//!    boundary with the same or the opposite outward normal.
//! 5. The operation keeps some regions (turning those of the tool inside out for a difference),
//!    which are sewn along shared edges into shells and solids.
//!
//! The result is a compound of solids, as in OCCT; it is empty when nothing remains.

mod arguments;
mod build;
mod cells;
mod classify;
mod far;
mod fuse;
mod fuzzy;
mod geom;

/// A point well inside the region closed polylines (each ending where it starts, in a
/// surface's parameters) bound: of the middles of its spans along a few lines across it, the
/// one farthest from the boundary.
pub(crate) fn interior_uv(polylines: &[Vec<[f64; 2]>]) -> [f64; 2] {
    geom::interior_point(polylines)
}
mod graph;
mod interference;
mod intersect;
mod march;
mod mixed;
mod split;
mod volume;

use crate::error::{Error, Result};
use crate::math::BoundingBox;
use crate::topology::{Model, Shape};

pub use arguments::{combine, Combination};
pub use cells::{cells, Cell};
pub use fuzzy::snapped;
pub use interference::{interference, Interference};
pub use volume::{enclosed_solids, enclosed_solids_unsplit};
pub(crate) use volume::enclosed_solids_at_once;
pub(crate) use cells::point_inside;
pub(crate) use geom::on_tol;

use build::Kept;
use classify::State;
use graph::Graph;

/// The union of `a` and `b`.
#[doc(alias = "BRepAlgoAPI_Fuse")]
#[doc(alias = "fuse")]
pub fn union(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    // An operand may be a compound of solids that overlap one another: like OCCT's general
    // fuse, every solid is fused with those it meets (boxes overlapping or touching), one
    // after another; groups apart from each other are only gathered. (Two solids apart are
    // gathered too, and any shapes apart: each as it was, placed where it was, as OCCT's
    // fuse keeps what it does not change.)
    let apart = {
        let boxed = |s: Shape| crate::properties::tight_bounding_box(model, s).enlarged(geom::on_tol());
        let (ba, bb) = (boxed(a), boxed(b));
        !ba.is_empty() && !bb.is_empty() && !ba.overlaps(&bb)
    };
    // (Other shapes apart, faces as solids: gathered as they are.)
    if apart && (only_solids(model, a).is_none() || only_solids(model, b).is_none()) {
        let result = model.new_empty(crate::topology::NodeKind::Compound);
        for operand in [a, b] {
            if operand.kind() == crate::topology::NodeKind::Compound {
                for child in model.children(operand).collect::<Vec<_>>() {
                    model.add(result, child)?;
                }
            } else {
                model.add(result, operand)?;
            }
        }
        return Ok(result);
    }
    if let (Some(sa), Some(sb)) = (only_solids(model, a), only_solids(model, b)) {
        if sa.len() > 1 || sb.len() > 1 || apart {
            let solids: Vec<Shape> = sa.into_iter().chain(sb).collect();
            let boxes: Vec<BoundingBox> = solids.iter().map(|&s| crate::properties::tight_bounding_box(model, s).enlarged(geom::on_tol())).collect();
            let mut group: Vec<usize> = (0..solids.len()).collect();
            fn root(group: &mut [usize], mut i: usize) -> usize {
                while group[i] != i {
                    group[i] = group[group[i]];
                    i = group[i];
                }
                i
            }
            for i in 0..solids.len() {
                for j in i + 1..solids.len() {
                    if boxes[i].overlaps(&boxes[j]) {
                        let (ri, rj) = (root(&mut group, i), root(&mut group, j));
                        group[ri.max(rj)] = ri.min(rj);
                    }
                }
            }
            let result = model.new_empty(crate::topology::NodeKind::Compound);
            for r in 0..solids.len() {
                let members: Vec<usize> = (0..solids.len()).filter(|&i| root(&mut group, i) == r).collect();
                if members.is_empty() {
                    continue;
                }
                // (In layers of solids whose boxes are apart, each layer fused at once with
                // what the ones before make: many solids each meeting only one other, as
                // tiles laid on a floor, in a few operations, not one each on the whole
                // growing.)
                let mut layers: Vec<Vec<usize>> = Vec::new();
                for &i in &members {
                    match layers.iter_mut().find(|layer| layer.iter().all(|&j| !boxes[i].overlaps(&boxes[j]))) {
                        Some(layer) => layer.push(i),
                        None => layers.push(vec![i]),
                    }
                }
                let gathered = |model: &mut Model, layer: &[usize]| -> Result<Shape> {
                    if let [one] = layer {
                        return Ok(solids[*one]);
                    }
                    let c = model.new_empty(crate::topology::NodeKind::Compound);
                    for &i in layer {
                        model.add(c, solids[i])?;
                    }
                    Ok(c)
                };
                let mut fused = gathered(model, &layers[0])?;
                for layer in &layers[1..] {
                    let next = gathered(model, layer)?;
                    fused = run(model, fused, next, Operation::Union)?;
                }
                for solid in model.explore(fused, crate::topology::NodeKind::Solid).collect::<Vec<_>>() {
                    model.add(result, solid)?;
                }
            }
            return Ok(result);
        }
    }
    run(model, a, b, Operation::Union)
}

/// The solids of `shape` when every face of it belongs to one of them.
fn only_solids(model: &Model, shape: Shape) -> Option<Vec<Shape>> {
    use crate::topology::NodeKind;
    let solids: Vec<Shape> = model.explore(shape, NodeKind::Solid).collect();
    let faces_in_solids: usize = solids.iter().map(|&s| model.explore(s, NodeKind::Face).count()).sum();
    (!solids.is_empty() && model.explore(shape, NodeKind::Face).count() == faces_in_solids).then_some(solids)
}

/// `a` minus `b`.
#[doc(alias = "BRepAlgoAPI_Cut")]
#[doc(alias = "cut")]
pub fn difference(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    if let Some(each) = each_solid_apart(model, a, b, Operation::Difference)? {
        return Ok(each);
    }
    run(model, a, b, Operation::Difference)
}

thread_local! {
    /// Set while [`cells`] splits shapes (which itself runs these operations): there an
    /// operand's solids are taken one by one.
    static SPLITTING_CELLS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// For operands made of solids meeting one another: the operation by [`combine`], their
/// solids split by one another (`None` otherwise). (Within [`cells`]: the operation on each
/// solid of `a`, gathered.)
fn each_solid_apart(model: &mut Model, a: Shape, b: Shape, operation: Operation) -> Result<Option<Shape>> {
    let meets = |model: &Model, s: Shape| only_solids(model, s).is_some_and(|solids| solids.len() > 1 && arguments::meeting(model, &solids));
    if SPLITTING_CELLS.with(|f| f.get()) {
        if !meets(model, a) {
            return Ok(None);
        }
        let gathered = model.new_empty(crate::topology::NodeKind::Compound);
        for solid in only_solids(model, a).expect("checked") {
            let part = run(model, solid, b, operation)?;
            for piece in model.explore(part, crate::topology::NodeKind::Solid).collect::<Vec<_>>() {
                model.add(gathered, piece)?;
            }
        }
        return Ok(Some(gathered));
    }
    if !(meets(model, a) || meets(model, b)) || only_solids(model, a).is_none() || only_solids(model, b).is_none() {
        return Ok(None);
    }
    let combination = match operation {
        Operation::Intersection => Combination::Intersection,
        Operation::Difference => Combination::Difference,
        Operation::Union => return Ok(None),
    };
    Ok(Some(combine(model, &[a], &[b], combination)?))
}

/// [`difference`] built once, at the current tolerances, its shells closed or not: for
/// callers keeping only some of its faces (and checking what they build from them).
pub(crate) fn difference_once(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    Ok(run_tolerant(model, a, b, Operation::Difference)?.0)
}

/// [`difference`], given `first`, [`difference_once`]'s result for the same operands: what
/// it would have built again first is not.
pub(crate) fn difference_after(model: &mut Model, a: Shape, b: Shape, first: Shape) -> Result<Shape> {
    // (Where `difference` takes another way, that way.)
    let meets = |model: &Model, s: Shape| only_solids(model, s).is_some_and(|solids| solids.len() > 1 && arguments::meeting(model, &solids));
    if SPLITTING_CELLS.with(|f| f.get()) || meets(model, a) || meets(model, b) || is_empty(model, a) || is_empty(model, b) {
        return difference(model, a, b);
    }
    if crate::tolerance::loosened() || closed_shells(model, first) {
        return Ok(first);
    }
    // (As [`run_with_sections`] goes on: again with the looser tolerances, taken if its only
    // problem is a face's boundary crossing itself; else the first.)
    let crossing = |p: &crate::validity::Problem| p.kind == crate::validity::ProblemKind::SelfIntersection;
    match crate::tolerance::with_loosened(true, || run_tolerant(model, a, b, Operation::Difference)) {
        Ok((r, _, _)) if crate::validity::check(model, r).errors().all(crossing) => Ok(r),
        _ => Ok(first),
    }
}

/// The intersection of `a` and `b`.
#[doc(alias = "BRepAlgoAPI_Common")]
#[doc(alias = "common")]
pub fn intersection(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    if let Some(each) = each_solid_apart(model, a, b, Operation::Intersection)? {
        return Ok(each);
    }
    run(model, a, b, Operation::Intersection)
}

/// The recordings under way ([`record_begin`] to [`record_end`]), one per model: each operand
/// face and edge, with a result face or edge made from it (by every build, retries and
/// pieces of a split operation included). (Shared by threads: builds may run on any.)
static RECORDING: std::sync::Mutex<Vec<(crate::topology::ModelUid, Vec<(crate::topology::NodeId, crate::topology::NodeId)>)>> = std::sync::Mutex::new(Vec::new());

/// How many recordings are under way (checked first: with none, recording costs nothing).
static RECORDING_ON: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Records a made result face or edge's source in `model` (by its uid), while it is being
/// recorded (by the booleans, sewing and fillets).
pub(crate) fn record(model: crate::topology::ModelUid, source: crate::topology::NodeId, made: crate::topology::NodeId) {
    if RECORDING_ON.load(std::sync::atomic::Ordering::Relaxed) == 0 {
        return;
    }
    let mut all = RECORDING.lock().unwrap_or_else(|p| p.into_inner());
    if let Some((_, list)) = all.iter_mut().find(|(m, _)| *m == model) {
        list.push((source, made));
    }
}

/// Starts recording where the result faces and edges of `model`'s booleans (and sewing,
/// fillets and chamfers) come from (for [`record_end`]).
pub fn record_begin(model: &Model) {
    let mut all = RECORDING.lock().unwrap_or_else(|p| p.into_inner());
    all.retain(|(m, _)| *m != model.uid());
    all.push((model.uid(), Vec::new()));
    RECORDING_ON.store(all.len(), std::sync::atomic::Ordering::Relaxed);
}

/// Stops recording: each operand face or edge (unplaced, by node) and a face or edge of
/// `result` made from it since [`record_begin`] (an operand's face or edge the result keeps
/// as it is maps to itself). Only what `result` has: what an attempt that was not kept, or
/// a step after the boolean, made is left out.
#[doc(alias = "BRepAlgoAPI_BuilderAlgo::Modified")]
pub fn record_end(model: &Model, result: Shape) -> Vec<(crate::topology::NodeId, crate::topology::NodeId)> {
    let list = {
        let mut all = RECORDING.lock().unwrap_or_else(|p| p.into_inner());
        let at = all.iter().position(|(m, _)| *m == model.uid());
        let list = at.map(|at| all.remove(at).1).unwrap_or_default();
        RECORDING_ON.store(all.len(), std::sync::atomic::Ordering::Relaxed);
        list
    };
    let present: crate::hashing::FastSet<crate::topology::NodeId> = [crate::topology::NodeKind::Face, crate::topology::NodeKind::Edge].into_iter().flat_map(|k| model.explore(result, k).map(|s| s.node).collect::<Vec<_>>()).collect();
    // (Through the steps of an operation made of several, as a fillet's tools cut one after
    // another: each made face or edge back to the inputs it came from, those no step made.)
    let mut sources: crate::hashing::FastMap<crate::topology::NodeId, Vec<crate::topology::NodeId>> = Default::default();
    for &(source, made) in &list {
        if source != made {
            sources.entry(made).or_default().push(source);
        }
    }
    let mut out = Vec::new();
    let mut seen = crate::hashing::FastSet::default();
    for &(source, made) in &list {
        if !present.contains(&made) {
            continue;
        }
        let mut stack = vec![source];
        let mut visited = crate::hashing::FastSet::default();
        while let Some(n) = stack.pop() {
            if !visited.insert(n) {
                continue;
            }
            match sources.get(&n) {
                Some(earlier) if n != made => stack.extend(earlier.iter().copied()),
                _ => {
                    if seen.insert((n, made)) {
                        out.push((n, made));
                    }
                }
            }
        }
    }
    out
}

/// The edges where `a` and `b` meet: where their faces cross, and pieces of edges of either
/// lying on the other; and the points where they touch elsewhere.
#[doc(alias = "BRepAlgoAPI_Section")]
pub fn section(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    if [a, b].iter().any(|&s| is_empty(model, s)) {
        return Ok(model.new_empty(crate::topology::NodeKind::Compound));
    }
    // Operands of solids meeting one another: the edges where any two of their solids meet.
    let meets = |model: &Model, s: Shape| only_solids(model, s).is_some_and(|solids| solids.len() > 1 && arguments::meeting(model, &solids));
    if !SPLITTING_CELLS.with(|f| f.get()) && only_solids(model, a).is_some() && only_solids(model, b).is_some() && (meets(model, a) || meets(model, b)) {
        return combine(model, &[a], &[b], Combination::Section);
    }
    let graph = operands(model, a, b, false)?;
    mixed::section(model, graph)
}

/// [`union`], and the result's edges where the operands' faces meet.
#[doc(alias = "BRepAlgoAPI_BooleanOperation::SectionEdges")]
pub fn union_with_sections(model: &mut Model, a: Shape, b: Shape) -> Result<(Shape, Vec<Shape>)> {
    run_with_sections(model, a, b, Operation::Union)
}

/// [`difference`], and the result's edges where the operands' faces meet.
pub fn difference_with_sections(model: &mut Model, a: Shape, b: Shape) -> Result<(Shape, Vec<Shape>)> {
    run_with_sections(model, a, b, Operation::Difference)
}

/// Where a point lies with respect to a solid.
#[doc(alias = "TopAbs_State")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointLocation {
    Inside,
    Outside,
    /// On its boundary (within 1e-6).
    Boundary,
}

/// Where `p` lies with respect to the solids of `shape` (inside any of them counts).
#[doc(alias = "BRepClass3d_SolidClassifier")]
#[doc(alias = "isInside")]
pub fn classify_point(model: &Model, shape: Shape, p: crate::math::Point) -> Result<PointLocation> {
    Ok(PointClassifier::new(model, shape)?.locate(p))
}

/// [`classify_point`] for several points of one shape: the shape is prepared once.
pub struct PointClassifier {
    graph: Graph,
}

impl PointClassifier {
    pub fn new(model: &Model, shape: Shape) -> Result<PointClassifier> {
        // (A large shape's faces' samples and edges' boxes found for all at once first.)
        let edges = model.unique_subshapes(shape, crate::topology::NodeKind::Edge);
        far::prepare(model, [shape, shape], [&edges, &edges]);
        let mut graph = Graph::default();
        graph.add_operand(model, shape, 0, &Default::default())?;
        Ok(PointClassifier { graph })
    }

    /// How many times the shape's faces wind round `p` (off them): 1 inside a solid, 0
    /// outside, more where a shell passing through itself covers it more than once; `None`
    /// where that cannot be told.
    pub(crate) fn winding(&self, p: crate::math::Point) -> Option<i32> {
        self.graph.winding(0, p)
    }

    /// Where `p`, a point of the surface of a face the shape has, lies with respect to that
    /// face's bounds: inside, on its boundary (within 1e-6) or outside (as OCCT's
    /// BRepClass_FaceClassifier tells IN, ON and OUT); inside or on any of the shape's faces
    /// counts.
    pub(crate) fn locate_on_faces(&self, p: crate::math::Point) -> PointLocation {
        let mut found = PointLocation::Outside;
        for face in &self.graph.faces {
            match face.locate(p, &self.graph.edges, &self.graph.vertices) {
                graph::Location::Inside => return PointLocation::Inside,
                graph::Location::Boundary => found = PointLocation::Boundary,
                graph::Location::Outside => {}
            }
        }
        found
    }

    /// Where `p` lies with respect to the shape's solids.
    pub fn locate(&self, p: crate::math::Point) -> PointLocation {
        match self.graph.classify(p, crate::math::Vector::new(0.0, 0.0, 1.0), 0, geom::ON_TOL) {
            State::Inside => PointLocation::Inside,
            State::Outside => PointLocation::Outside,
            State::OnSame | State::OnOpposite => PointLocation::Boundary,
        }
    }
}

/// True for a shape with nothing in it: no faces, edges or vertices (a face without edges,
/// as of a whole sphere, is something).
fn is_empty(model: &Model, s: Shape) -> bool {
    use crate::topology::NodeKind;
    [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex].iter().all(|&k| model.explore(s, k).next().is_none())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Operation {
    Union,
    Difference,
    Intersection,
}

impl Operation {
    /// Whether a region of operand `k` (0 for the object, 1 for the tool) in `state`
    /// belongs to the result, and if so whether it is turned inside out. Regions shared by
    /// both operands are taken from the object.
    fn keeps(self, k: usize, state: State) -> Option<bool> {
        use State::*;
        match (self, k, state) {
            (Operation::Union, _, Outside) | (Operation::Union, 0, OnSame) => Some(false),
            (Operation::Intersection, _, Inside) | (Operation::Intersection, 0, OnSame) => Some(false),
            (Operation::Difference, 0, Outside) | (Operation::Difference, 0, OnOpposite) => Some(false),
            (Operation::Difference, 1, Inside) => Some(true),
            _ => None,
        }
    }
}

pub(super) fn run(model: &mut Model, a: Shape, b: Shape, operation: Operation) -> Result<Shape> {
    // An intersection of solids, one placed (the larger, moved about, as when checks slide
    // one shape along another): done where it is unplaced, the other moved to it, the
    // result placed back. Its faces' kept data (far faces, regions, boxes) then serve each
    // time.
    let rigid = |t: &crate::math::Transform| !t.is_identity() && (t.scale_factor() - 1.0).abs() <= 1e-12 && !t.is_negative();
    let faces = |s: Shape| model.explore(s, crate::topology::NodeKind::Face).count();
    let (fa, fb) = (faces(a), faces(b));
    let solids = |s: Shape| only_solids(model, s).is_some();
    let local = if operation != Operation::Intersection || !solids(a) || !solids(b) {
        None
    } else if rigid(&a.placement) && fa >= fb {
        Some((a.placement, Shape { placement: crate::math::Transform::IDENTITY, ..a }, b.moved(&a.placement.inverse())))
    } else if rigid(&b.placement) && fb > fa {
        Some((b.placement, a.moved(&b.placement.inverse()), Shape { placement: crate::math::Transform::IDENTITY, ..b }))
    } else {
        None
    };
    if let Some((t, a, b)) = local {
        let result = run_with_sections(model, a, b, operation)?.0;
        return Ok(result.moved(&t));
    }
    Ok(run_with_sections(model, a, b, operation)?.0)
}

fn operands(model: &Model, a: Shape, b: Shape, use_far: bool) -> Result<Graph> {
    operands_placed(model, a, b, use_far, [false, false])
}

/// [`operands`]; with `placed_far[k]`, faces of operand `k` placed (moved) may be left out
/// as far ones too (where the operation keeps nothing outside the other operand of `k`:
/// they are only counted, never built into the result).
fn operands_placed(model: &Model, a: Shape, b: Shape, use_far: bool, placed_far: [bool; 2]) -> Result<Graph> {
    let mut graph = Graph::default();
    // The operands' edges, found once for the steps that go through them all.
    let edges = [a, b].map(|s| model.unique_subshapes(s, crate::topology::NodeKind::Edge));
    let prepared = far::prepare(model, [a, b], [&edges[0], &edges[1]]);
    let far = if use_far { far::choose(model, a, b, &prepared, placed_far) } else { Default::default() };
    graph.add_operand(model, a, 0, &far[0])?;
    graph.add_operand(model, b, 1, &far[1])?;
    for k in 0..2 {
        // (Not on free-form faces, where it is not needed.)
        let freeform = graph.faces.iter().filter(|f| f.operand == k).any(|f| geom::freeform(&f.surface))
            || graph.far.iter().filter(|f| f.operand == k).any(|f| geom::freeform(&f.data.surface));
        if freeform {
            // (A box holding its faces: enough to tell points far outside it.)
            let mut held = crate::math::BoundingBox::EMPTY;
            for b in &graph.face_bounds[k] {
                held.add_box(b);
            }
            if !held.is_empty() && held.is_bounded() {
                graph.tight_boxes[k] = Some((held.min, held.max));
            }
            continue;
        }
        // The faces' grid boxes (holding them) spare finding the extremes of faces that
        // lie within the box of the others.
        let tight = crate::properties::tight_bounding_box_within(model, [a, b][k], &edges[k], &graph.face_bounds[k]);
        if !tight.is_empty() {
            graph.tight_boxes[k] = Some((tight.min, tight.max));
        }
    }
    Ok(graph)
}

fn run_with_sections(model: &mut Model, a: Shape, b: Shape, operation: Operation) -> Result<(Shape, Vec<Shape>)> {
    // An empty operand (as OCCT): nothing added, nothing taken away, nothing in common.
    let empty = |s: Shape, model: &Model| is_empty(model, s);
    match (empty(a, model), empty(b, model), operation) {
        (_, true, Operation::Union | Operation::Difference) => return Ok((a, Vec::new())),
        (true, false, Operation::Union) => return Ok((b, Vec::new())),
        (true, _, _) | (_, true, _) => return Ok((model.new_empty(crate::topology::NodeKind::Compound), Vec::new())),
        _ => {}
    }
    // Where the result's shells are not closed (or none could be built), again with the
    // looser tolerances: shapes read from files may put points meant to coincide (a
    // cylinder through a box's corners) only nearly together. That result if it is valid
    // (or its only problem is a face's boundary crossing itself: still closed, better than
    // the first).
    let first = run_tolerant(model, a, b, operation);
    if crate::tolerance::loosened() || first.as_ref().is_ok_and(|r| r.2) {
        return first.map(|(r, s, _)| (r, s));
    }
    let crossing = |p: &crate::validity::Problem| p.kind == crate::validity::ProblemKind::SelfIntersection;
    match crate::tolerance::with_loosened(true, || run_tolerant(model, a, b, operation)) {
        Ok((r, s, _)) if crate::validity::check(model, r).errors().all(crossing) => return Ok((r, s)),
        _ => {}
    }
    // (A union or an intersection the other way round, the same operation: where the
    // operands meet in pieces about as small as the tolerance (walls of plates nearly
    // touching at a corner), which of them is traced first may decide whether the pieces
    // close up. Taken where its shells close. Only where the first was built, its shells
    // left open: one that failed outright may take far longer the other way round.)
    if matches!(operation, Operation::Union | Operation::Intersection) && first.is_ok() {
        if let Ok((r, s, true)) = run_tolerant(model, b, a, operation) {
            return Ok((r, s));
        }
    }
    first.map(|(r, s, _)| (r, s))
}

/// Whether every shell of `shape` is closed (each edge used by its faces once each way).
pub(crate) fn closed_shells(model: &Model, shape: Shape) -> bool {
    crate::validity::shells_closed(model, shape)
}

/// [`run_with_sections`] at the current tolerances, with whether the result's shells are
/// closed.
fn run_tolerant(model: &mut Model, a: Shape, b: Shape, operation: Operation) -> Result<(Shape, Vec<Shape>, bool)> {
    // With faces far from the other operand left out of the graph, unless that turns out
    // not to hold (an open operand, a shared edge the result cannot keep): then with all.
    if let Some(result) = run_graph(model, a, b, operation, true)? {
        return Ok(result);
    }
    run_graph(model, a, b, operation, false)?.ok_or_else(|| Error::InvalidInput("a boolean operation that could not be built".into()))
}

/// The operation on the operands' graph; `None` if far faces (with `use_far`) cannot be
/// left out after all.
fn run_graph(model: &mut Model, a: Shape, b: Shape, operation: Operation, use_far: bool) -> Result<Option<(Shape, Vec<Shape>, bool)>> {
    let placed_far = [0, 1].map(|k| operation.keeps(k, State::Outside).is_none());
    let mut graph = operands_placed(model, a, b, use_far, placed_far)?;
    if !graph.far.is_empty() && (graph.open != [false, false] || graph.dimension != [3, 3] || !graph.only_solids()) {
        return Ok(None);
    }
    if !graph.only_solids() {
        return Ok(Some((mixed::run(model, graph, operation)?, Vec::new(), true)));
    }
    let (pieces, arcs_of, regions, far_states) = classified(model, &mut graph);
    // Faces far from the other operand lie outside it, whole (or inside it, those apart from
    // its faces where the operands overlap).
    let Some(far_states) = far_states else { return Ok(None) };
    let mut kept = Vec::new();
    for (fi, region, state) in regions {
        if let Some(flip) = operation.keeps(graph.faces[fi].operand, state) {
            kept.push(Kept { face: fi, region, flip, far: false });
        }
    }
    for (i, f) in graph.far.iter().enumerate() {
        if let Some(flip) = operation.keeps(f.operand, far_states[i]) {
            kept.push(Kept { face: i, region: split::Region::whole(), flip, far: true });
        }
    }
    // In the operands' order of faces.
    kept.sort_by_key(|k| if k.far { graph.far[k.face].seq } else { graph.faces[k.face].seq });
    let invalid_operands = |model: &Model| [a, b].iter().any(|&s| !crate::validity::check(model, s).is_valid());
    build::build(model, &graph, &pieces, &arcs_of, &kept, &invalid_operands)
}

/// The operands' faces split along where they meet, and each region with where it lies
/// with respect to the other operand; with where the far faces lie ([`far::states`]).
#[allow(clippy::type_complexity)]
fn classified(model: &Model, graph: &mut Graph) -> (Vec<split::Piece>, Vec<Vec<split::ArcOnFace>>, Vec<(usize, split::Region, State)>, Option<Vec<State>>) {
    graph.intersect_edges_with_faces();
    graph.intersect_edges_on_faces();
    graph.share_vertices();
    graph.intersect_faces();
    graph.drop_spurious_paves();
    let (mut pieces, of_edge) = graph.pieces();
    // (A face's boundary arcs over whole edges are the same each time: kept in the model
    // with the face, while their departure step stays the same, as when the face is moved.)
    let whole: Vec<Option<std::sync::Arc<WholeArcs>>> = crate::parallel::map(&(0..graph.faces.len()).collect::<Vec<_>>(), |&fi| {
        let id = graph.faces[fi].regions_key?;
        let step = graph.arcs_step(fi);
        let kept = model.derived(id.into(), model.face_fingerprint(id), || WholeArcs(step, graph.whole_boundary_arcs(fi)));
        ((kept.0 - step).abs() <= 1e-9 * step).then_some(kept)
    });
    let arcs_of: Vec<_> = (0..graph.faces.len()).map(|fi| graph.arcs_on_face(fi, &mut pieces, &of_edge, whole[fi].as_ref().map(|w| &w.1[..]))).collect();
    // Each face's regions and where they lie, found for all faces at once.
    let faces: Vec<usize> = (0..arcs_of.len()).collect();
    let graph: &Graph = graph;
    // (Narrower, as twice its area over its perimeter, than a thousandth of its size: its
    // holes' areas and perimeters counted, as a ring's between circles close together.)
    let thin = |polygons: &[Vec<geom::Uv>]| {
        let (lo, hi) = geom::polyline_box(&polygons[0]);
        let perimeter: f64 = polygons.iter().flat_map(|p| p.windows(2)).map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1])).sum();
        let area: f64 = polygons.iter().map(|p| geom::signed_area(p)).sum();
        2.0 * area.abs() < 1e-3 * (hi[0] - lo[0]).hypot(hi[1] - lo[1]) * perimeter
    };
    let found = crate::parallel::map(&faces, |&fi| {
        let arcs = &arcs_of[fi];
        let face = &graph.faces[fi];
        let scale = face.bbox.0.distance(face.bbox.1);
        // A face met by nothing (its arcs its own edges, each whole, in order) is one region,
        // the same each time: kept in the model with the face.
        let untouched = face.regions_key.filter(|_| {
            arcs.len() == face.boundary.len()
                && arcs.iter().zip(&face.boundary).all(|(a, b)| {
                    let (p, e) = (&pieces[a.piece], &graph.edges[b.edge]);
                    p.edge == b.edge && p.first == e.first && p.last == e.last
                })
        });
        let regions_of_face = match untouched {
            Some(id) => {
                let fingerprint = model.face_fingerprint(id);
                (*model.derived(id.into(), fingerprint, || FaceRegions(split::split_face(arcs, &pieces, &face.polylines, scale, false)))).0.clone()
            }
            None => {
                let mut regions = split::split_face(arcs, &pieces, &face.polylines, scale, false);
                // (Regions thinner than their polylines' sampling: their interior points from
                // denser ones, the coarse ones' possibly off the region.)
                for r in &mut regions {
                    if thin(&r.polygons) {
                        if let Some(p) = graph.thin_interior(fi, arcs, &pieces, r) {
                            r.interior = p;
                        }
                    }
                }
                regions
            }
        };
        regions_of_face.into_iter().map(|region| (fi, region, State::Outside)).collect::<Vec<_>>()
    });
    let mut regions: Vec<(usize, split::Region, State)> = found.into_iter().flatten().collect();
    // Where each region lies: regions joined across edges the other operand does not touch
    // lie alike, so of a large group of them only some are classified (a spread of them),
    // unless they disagree or one lies on the other operand (then all are). (Far faces joined
    // to them count as classified regions of their groups, as they would in the graph.)
    let far_states = far::states(graph);
    let roots = clean_components(graph, &pieces, &arcs_of, &regions, far_states.as_deref());
    let mut groups: crate::hashing::FastMap<usize, Vec<usize>> = Default::default();
    for (ri, &root) in roots[..regions.len()].iter().enumerate() {
        groups.entry(root).or_default().push(ri);
    }
    let mut far_votes: crate::hashing::FastMap<usize, Vec<State>> = Default::default();
    for (&root, &state) in roots[regions.len()..].iter().zip(far_states.iter().flatten()) {
        far_votes.entry(root).or_default().push(state);
    }
    let votes = |g: &[usize]| far_votes.get(&roots[g[0]]).map_or(&[][..], |v| &v[..]);
    let classify = |ri: usize| {
        let (fi, region, _) = &regions[ri];
        let face = &graph.faces[*fi];
        let at = |uv: geom::Uv| {
            let p = geom::surface_point(&face.surface, uv);
            let outward = geom::normal_at(&face.surface, uv) * face.outward;
            graph.classify_on(p, outward, 1 - face.operand, face.drift, Some(*fi))
        };
        let state = at(region.interior);
        if !matches!(state, State::OnSame | State::OnOpposite) {
            return state;
        }
        // (On the other operand at the point taken: lying on it, or only touching it there,
        // as a face tangent to the other's along a line through that point. Other points of
        // the region then tell: lying on it, all are on it.)
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in region.polygons.iter().flatten() {
            lo = [lo[0].min(p[0]), lo[1].min(p[1])];
            hi = [hi[0].max(p[0]), hi[1].max(p[1])];
        }
        const SPREAD: [f64; 5] = [0.2763, 0.4127, 0.5891, 0.7177, 0.1531];
        let mut off = Vec::new();
        'samples: for &a in &SPREAD {
            for &b in SPREAD.iter().rev() {
                let uv = [lo[0] + (hi[0] - lo[0]) * a, lo[1] + (hi[1] - lo[1]) * b];
                if !geom::inside_polylines(&region.polygons, uv) {
                    continue;
                }
                let s = at(uv);
                if matches!(s, State::OnSame | State::OnOpposite) {
                    return state;
                }
                off.push(s);
                if off.len() >= 3 {
                    break 'samples;
                }
            }
        }
        // (Off it at all others: only touching it. They agree, or the point taken stands.)
        match off.first() {
            Some(&s) if off.iter().all(|&t| t == s) => s,
            _ => state,
        }
    };
    const SAMPLES: usize = 7;
    let sample = |g: &[usize]| -> Vec<usize> {
        if g.len() <= SAMPLES { g.to_vec() } else { (0..SAMPLES).map(|k| g[k * g.len() / SAMPLES]).collect() }
    };
    let sampled: Vec<usize> = groups.values().flat_map(|g| sample(g)).collect();
    let states = crate::parallel::map(&sampled, |&ri| classify(ri));
    let mut known: crate::hashing::FastMap<usize, State> = sampled.iter().copied().zip(states).collect();
    let mut rest: Vec<usize> = Vec::new();
    for g in groups.values() {
        // (Unless the sample agrees, all: a region lying on the other operand, or a group
        // not all on one side after all.)
        let s = sample(g);
        let unanimous = s.iter().all(|ri| known[ri] == known[&s[0]]) && votes(g).iter().all(|&v| v == known[&s[0]]) && matches!(known[&s[0]], State::Inside | State::Outside);
        if g.len() > SAMPLES && !unanimous {
            rest.extend(g.iter().filter(|ri| !known.contains_key(ri)));
        }
    }
    let states = crate::parallel::map(&rest, |&ri| classify(ri));
    known.extend(rest.iter().copied().zip(states));
    // Each group's regions as most of those classified (those lying on the other operand
    // apart).
    for g in groups.values() {
        let (mut inside, mut outside) = (0, 0);
        for state in g.iter().filter_map(|ri| known.get(ri)).chain(votes(g)) {
            match state {
                State::Inside => inside += 1,
                State::Outside => outside += 1,
                _ => {}
            }
        }
        let majority = (inside != outside).then(|| if inside > outside { State::Inside } else { State::Outside });
        for &ri in g {
            regions[ri].2 = match (known.get(&ri), majority) {
                (Some(&s), _) if !matches!(s, State::Inside | State::Outside) => s,
                (_, Some(m)) => m,
                (Some(&s), None) => s,
                (None, None) => State::Outside,
            };
        }
    }
    (pieces, arcs_of, regions, far_states)
}

/// `a` and `b` split by one another: the solids inside both, inside `a` only and inside
/// `b` only, sharing the faces, edges and vertices where they meet (each made once).
#[doc(alias = "BOPAlgo_Builder")]
pub(crate) fn partition(model: &mut Model, a: Shape, b: Shape) -> Result<[Shape; 3]> {
    let mut graph = operands(model, a, b, false)?;
    if !graph.only_solids() {
        return Err(Error::InvalidInput("a partition of shapes other than solids".into()));
    }
    let (pieces, arcs_of, regions, _) = classified(model, &mut graph);
    // Each region once (coplanar ones of `a` only, the other operand's being the same),
    // and the cells it bounds: as it is, or turned inside out.
    use State::*;
    let mut all = Vec::new();
    let mut cells: [Vec<(usize, bool)>; 3] = Default::default();
    let mut regions = regions;
    regions.sort_by_key(|r| graph.faces[r.0].seq);
    for (fi, region, state) in regions {
        let operand = graph.faces[fi].operand;
        let uses: &[(usize, bool)] = match (operand, state) {
            (0, Inside) => &[(0, false), (2, true)],
            (0, Outside) => &[(1, false)],
            (0, OnSame) => &[(0, false)],
            (0, OnOpposite) => &[(1, false), (2, true)],
            (_, Inside) => &[(0, false), (1, true)],
            (_, Outside) => &[(2, false)],
            _ => &[],
        };
        if uses.is_empty() {
            continue;
        }
        let index = all.len();
        all.push(Kept { face: fi, region, flip: false, far: false });
        for &(cell, flip) in uses {
            cells[cell].push((index, flip));
        }
    }
    let invalid_operands = |model: &Model| [a, b].iter().any(|&s| !crate::validity::check(model, s).is_valid());
    let built = build::build_partition(model, &graph, &pieces, &arcs_of, &all, &cells, &invalid_operands)?;
    Ok([built[0], built[1], built[2]])
}

/// The regions of a face met by nothing, kept in the model.
struct FaceRegions(Vec<split::Region>);

/// A face's [`split::WholeArc`]s, kept with it, and their departure step.
struct WholeArcs(f64, Vec<split::WholeArc>);

/// For each region, its group: regions of one operand joined across a piece of its own
/// edge that the other operand leaves alone (no section, not shared with it, not on its
/// faces) lie on the same side of it (each group then takes its majority, so that one
/// region classified wrongly, a ray miscounting near a thin feature, follows its
/// neighbours; and a large group needs only some of its regions classified). Then, given
/// where they lie (`far`), for each far face (after the regions), its group: joined to
/// regions across its edges alike, and to the far faces it shares an edge with (an unplaced
/// one, whose edges are the graph's).
fn clean_components(graph: &Graph, pieces: &[split::Piece], arcs_of: &[Vec<split::ArcOnFace>], regions: &[(usize, split::Region, State)], far: Option<&[State]>) -> Vec<usize> {
    use crate::hashing::{FastMap as HashMap, FastSet as HashSet};
    // Pieces touched by the other operand: shared with its pieces, or of edges on its faces.
    let mut by_rep: HashMap<usize, HashSet<Option<usize>>> = HashMap::default();
    for p in pieces {
        by_rep.entry(p.rep).or_default().insert(graph.edges[p.edge].operand);
    }
    let on_other: HashSet<usize> = graph.faces.iter().flat_map(|f| f.on_edges.iter().chain(&f.sections).copied()).collect();
    // (And pieces of edges only partly on the other's faces, taken as arcs on them.)
    let arcs_on_other: HashSet<usize> = arcs_of
        .iter()
        .enumerate()
        .flat_map(|(fi, arcs)| arcs.iter().map(move |a| (fi, a.piece)))
        .filter(|&(fi, pi)| graph.edges[pieces[pi].edge].operand == Some(1 - graph.faces[fi].operand))
        .map(|(_, pi)| pieces[pi].rep)
        .collect();
    let clean = |piece: usize, operand: usize| {
        let p = &pieces[piece];
        let e = &graph.edges[p.edge];
        e.operand == Some(operand) && e.curve.is_some() && !on_other.contains(&p.edge) && !arcs_on_other.contains(&p.rep) && by_rep.get(&p.rep).is_some_and(|ops| ops.len() == 1)
    };
    // The face of the other operand a piece lies on (where its edge does only in part, so
    // not found as on it): the regions on its two sides may lie on the two sides of that
    // face.
    let lies_on = |piece: &split::Piece, operand: usize| -> Option<usize> {
        use crate::geometry::ParametricCurve as _;
        let curve = graph.edges[piece.edge].curve.as_ref()?;
        let p = curve.point(0.5 * (piece.first + piece.last));
        graph.index().faces_containing(1 - operand, p).into_iter().find(|&g| {
            let face = &graph.faces[g];
            geom::signed_distance(&face.surface, p).abs() <= geom::on_tol() && face.locate(p, &graph.edges, &graph.vertices) != graph::Location::Outside
        })
    };
    // Which side of face `g` region `ri` leaves the piece of its arc `ai` (run `forward` or
    // not) towards: by a point a little into the region (to the left of its loop) from the
    // arc's middle; 0 where the region runs along the face.
    let side = |ri: usize, ai: usize, forward: bool, g: usize| -> i8 {
        let fi = regions[ri].0;
        let uv = &arcs_of[fi][ai].uv;
        if uv.len() < 2 {
            return 0;
        }
        let length: f64 = uv.windows(2).map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1])).sum();
        let k = (uv.len() / 2).max(1);
        let (a, b) = (uv[k - 1], uv[k]);
        let sign = if forward { 1.0 } else { -1.0 };
        let (dx, dy) = ((b[0] - a[0]) * sign, (b[1] - a[1]) * sign);
        let norm = dx.hypot(dy);
        if !(norm > 0.0) {
            return 0;
        }
        let step = 0.01 * length / norm;
        let q = [0.5 * (a[0] + b[0]) - dy * step, 0.5 * (a[1] + b[1]) + dx * step];
        let d = geom::signed_distance(&graph.faces[g].surface, geom::surface_point(&graph.faces[fi].surface, q));
        if d.abs() <= 10.0 * geom::on_tol() { 0 } else if d > 0.0 { 1 } else { -1 }
    };
    let joined = if far.is_some() { graph.far.len() } else { 0 };
    let mut parent: Vec<usize> = (0..regions.len() + joined).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while parent[r] != r {
            r = parent[r];
        }
        parent[i] = r;
        r
    }
    let mut first_by_rep: HashMap<usize, (usize, i8)> = HashMap::default();
    for (ri, (fi, region, _)) in regions.iter().enumerate() {
        let operand = graph.faces[*fi].operand;
        for &(ai, forward) in region.loops.iter().flatten() {
            let piece = arcs_of[*fi][ai].piece;
            if !clean(piece, operand) {
                continue;
            }
            let here = lies_on(&pieces[piece], operand).map_or(0, |g| side(ri, ai, forward, g));
            match first_by_rep.get(&pieces[piece].rep) {
                // (Not across a face of the other operand passing between them.)
                Some(&(_, there)) if here * there < 0 => {}
                Some(&(rj, _)) => {
                    let (a, b) = (find(&mut parent, ri), find(&mut parent, rj));
                    parent[a] = b;
                }
                None => {
                    first_by_rep.insert(pieces[piece].rep, (ri, here));
                }
            }
        }
    }
    // (A far face's edges are clear of the other operand's faces, as the face is: nothing of
    // it passes between the face and the regions beyond them.)
    let mut far_edges: HashMap<(usize, crate::topology::EdgeId), Vec<usize>> = HashMap::default();
    for (i, f) in graph.far.iter().enumerate().take(joined) {
        if f.placement == crate::math::Transform::IDENTITY {
            for e in f.data.edges.iter().filter(|e| e.curve.is_some()) {
                far_edges.entry((f.operand, e.edge)).or_default().push(i);
            }
        }
    }
    // (Far faces sharing an edge, which is not in the graph, alike too.)
    for faces in far_edges.values() {
        for w in faces.windows(2) {
            let (a, b) = (find(&mut parent, regions.len() + w[0]), find(&mut parent, regions.len() + w[1]));
            parent[a] = b;
        }
    }
    for (piece, p) in pieces.iter().enumerate() {
        let e = &graph.edges[p.edge];
        let (Some(operand), Some(origin)) = (e.operand, e.origin) else { continue };
        let Some(faces) = far_edges.get(&(operand, origin)) else { continue };
        if !clean(piece, operand) {
            continue;
        }
        for &i in faces {
            let node = regions.len() + i;
            match first_by_rep.get(&p.rep) {
                Some(&(rj, _)) => {
                    let (a, b) = (find(&mut parent, node), find(&mut parent, rj));
                    parent[a] = b;
                }
                None => {
                    first_by_rep.insert(p.rep, (node, 0));
                }
            }
        }
    }
    (0..regions.len() + joined).map(|ri| find(&mut parent, ri)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    use crate::math::{Direction, Frame, Point, Vector};
    use crate::primitives::{Cone, Cuboid, Cylinder, Sphere};
    use crate::topology::NodeKind;
    use crate::{properties, validity};

    /// Two unit boxes, one moved 5 along x, the other -15, fused: apart, each is kept as it
    /// was, placed where it was (as OCCT's fuse keeps them: a selection by nearness then
    /// finds the operand itself among the result's solids).
    /// Two balls turned alike, overlapping: what they share is the lens of closed form,
    /// `π (a + b − d)² (d² + 2d(a + b) − 3(a − b)²) / 12d`, whichever way they are turned
    /// (each ball's box was once cut short for some turns, and points of the other inside the
    /// ball taken to be outside it).
    #[test]
    fn turned_balls_share_their_lens() {
        let (a, b, d) = (1.5, 2.0, 0.8636146524836521);
        let r = [[0.9788398224233601, 0.18094247081611328, -0.09556371849799977], [0.1493341997978843, -0.9509430121093625, -0.2709366798554569], [-0.13989960260461975, 0.25093268016608106, -0.9578417881966399]];
        let column = |j: usize| Direction::new(r[0][j], r[1][j], r[2][j]);
        let mut model = Model::new();
        let first = Sphere::new(a).with_frame(Frame::new(Point::ORIGIN, column(2), column(0))).build(&mut model).unwrap();
        let second = Sphere::new(b).with_frame(Frame::new(Point::ORIGIN + Vector::from(column(0)) * d, column(2), column(0))).build(&mut model).unwrap();
        let common = intersection(&mut model, first, second).unwrap();
        let lens = PI * (a + b - d).powi(2) * (d * d + 2.0 * d * (a + b) - 3.0 * (a - b).powi(2)) / (12.0 * d);
        let volume = crate::properties::volume(&model, common).mass;
        assert!((volume - lens).abs() < 1e-6 * lens, "volume {volume}, lens {lens}");
    }

    #[test]
    fn solids_apart_are_kept_as_they_are() {
        let mut model = Model::new();
        let unit = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let near = unit.moved(&crate::math::Transform::translation(Vector::new(5.0, 0.0, 0.0)));
        let far = unit.moved(&crate::math::Transform::translation(Vector::new(-15.0, 0.0, 0.0)));
        let fused = union(&mut model, near, far).unwrap();
        let solids: Vec<Shape> = model.explore(fused, NodeKind::Solid).collect();
        assert_eq!(solids.len(), 2);
        assert!(solids.contains(&near) && solids.contains(&far));
        assert!((properties::volume(&model, fused).mass - 2.0).abs() < 1e-12);
    }

    /// A door (25 x 4.5 x 55) with its lower corners cut away by a block (5 x 4.5 x 20, as
    /// thick as the door) whose upper edge is chamfered back all but 1e-6 of its width,
    /// leaving a sliver of its top face 1e-6 wide standing on the door's side: the sliver is
    /// not taken as lying on the side (it stands across it, apart from it by more than the
    /// faces' tolerances), and stays in the result, closing it (OCCT's volume).
    #[test]
    fn a_sliver_standing_on_a_face_is_not_on_it() {
        let mut model = Model::new();
        let door = Cuboid::from_corners(Point::new(-12.5, -2.25, -27.5), Point::new(12.5, 2.25, 27.5)).build(&mut model).unwrap();
        let block = Cuboid::from_corners(Point::new(-2.5, -2.25, -10.0), Point::new(2.5, 2.25, 10.0)).build(&mut model).unwrap();
        let top_right: Vec<Shape> = model
            .unique_subshapes(block, NodeKind::Edge)
            .into_iter()
            .filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                (p.x - 2.5).abs() < 1e-9 && (q.x - 2.5).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9 && (q.z - 10.0).abs() < 1e-9
            })
            .collect();
        assert_eq!(top_right.len(), 1);
        let block = crate::blend::chamfer(&mut model, block, &top_right, 4.999999).unwrap();
        let cut = |model: &mut Model, body: Shape, at: crate::math::Transform| {
            let tool = crate::modify::transform(model, block, &at, Default::default()).unwrap();
            difference(model, body, tool).unwrap()
        };
        let lower_left = crate::math::Transform::translation(Vector::new(-10.0, 0.0, -17.5));
        let once = cut(&mut model, door, lower_left);
        assert!(validity::check(&model, once).is_valid(), "{:?}", validity::check(&model, once));
        let twice = cut(&mut model, once, crate::math::Transform::rotation(crate::math::Axis::Z, PI) * lower_left);
        assert!(validity::check(&model, twice).is_valid(), "{:?}", validity::check(&model, twice));
        let v = properties::volume(&model, twice).mass;
        assert!((v - 5400.0).abs() < 1e-3, "{v}");
    }

    fn cuboid(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
        let p = Point::new(corner.0, corner.1, corner.2);
        Cuboid::from_corners(p, p + Vector::new(size.0, size.1, size.2)).build(model).unwrap()
    }

    /// Checks the result is valid and has the expected volume, area and number of solids.
    fn check(model: &Model, shape: Shape, volume: f64, area: f64, solids: usize) {
        let report = validity::check(model, shape);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(model, shape).mass;
        let a = properties::area(model, shape).mass;
        assert!((v - volume).abs() <= 1e-9 * volume.max(1.0), "volume {v} != {volume}");
        assert!((a - area).abs() <= 1e-9 * area.max(1.0), "area {a} != {area}");
        assert_eq!(model.unique_subshapes(shape, NodeKind::Solid).len(), solids);
    }

    #[test]
    fn overlapping_boxes() {
        // Two unit cubes overlapping in a 0.5 x 1 x 1 slab.
        type Op = fn(&mut Model, Shape, Shape) -> Result<Shape>;
        let cases: [(Op, f64, f64); 3] = [(union, 1.5, 8.0), (difference, 0.5, 4.0), (intersection, 0.5, 4.0)];
        for (op, volume, area) in cases {
            let mut model = Model::new();
            let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
            let b = cuboid(&mut model, (0.5, 0.0, 0.0), (1.0, 1.0, 1.0));
            let result = op(&mut model, a, b).unwrap();
            check(&model, result, volume, area, 1);
        }
    }

    #[test]
    fn corner_overlap() {
        // Cubes overlapping in a 0.5 cube at a corner.
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let b = cuboid(&mut model, (0.5, 0.5, 0.5), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, b).unwrap();
        check(&model, fused, 2.0 - 0.125, 12.0 - 6.0 * 0.25, 1);
        let cut_result = difference(&mut model, a, b).unwrap();
        check(&model, cut_result, 1.0 - 0.125, 6.0, 1);
        let common_result = intersection(&mut model, a, b).unwrap();
        check(&model, common_result, 0.125, 1.5, 1);
    }

    #[test]
    fn mirrored_operands() {
        // A box reflected across x = 0 (left-handed plane frames) fused with an upright one.
        use crate::math::{Direction, Point, Transform};
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.5, 0.0, 0.0), (1.0, 1.0, 1.0));
        let mirror = Transform::mirror_plane(Point::ORIGIN, Direction::new(1.0, 0.0, 0.0));
        let a = crate::modify::transform(&mut model, a, &mirror, Default::default()).unwrap();
        let b = cuboid(&mut model, (-1.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, b).unwrap();
        check(&model, fused, 1.5, 8.0, 1);
        let cut_result = difference(&mut model, a, b).unwrap();
        check(&model, cut_result, 0.5, 4.0, 1);
    }

    #[test]
    fn a_twisted_ruled_loft_and_a_slab() {
        // A ruled loft between squares turned 30° apart (bilinear B-spline sides) and a slab
        // across it: the slab cuts it in two (OCCT gives these values).
        use crate::loft::{loft, LoftOptions};
        let mut model = Model::new();
        let square = |model: &mut Model, z: f64, turn: f64, r: f64| {
            let corners: Vec<Point> = (0..4)
                .map(|k| {
                    let a = turn + k as f64 * PI / 2.0;
                    Point::new(r * a.cos(), r * a.sin(), z)
                })
                .collect();
            let edges: Vec<Shape> = (0..4).map(|k| crate::make::edge_between(model, corners[k], corners[(k + 1) % 4]).unwrap()).collect();
            crate::make::wire(model, &edges).unwrap()
        };
        let (bottom, top) = (square(&mut model, 0.0, 0.0, 1.5), square(&mut model, 2.0, PI / 6.0, 1.0));
        let ruled = loft(&mut model, &[bottom, top], &LoftOptions { solid: true, ruled: true, ..Default::default() }).unwrap().shape;
        let slab = cuboid(&mut model, (-2.0, -2.0, 0.5), (4.0, 4.0, 1.0));
        let cut_result = difference(&mut model, ruled, slab).unwrap();
        check(&model, cut_result, 3.114182544031942, 19.65470359808339, 2);
        let fused = union(&mut model, ruled, slab).unwrap();
        check(&model, fused, 19.114182544031934, 55.632589281053406, 1);
        let common_result = intersection(&mut model, ruled, slab).unwrap();
        check(&model, common_result, 2.95120159687027, 13.035647114615353, 1);
    }

    #[test]
    fn empty_operands() {
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let nothing = model.new_empty(NodeKind::Compound);
        assert_eq!(difference(&mut model, a, nothing).unwrap(), a);
        assert_eq!(union(&mut model, nothing, a).unwrap(), a);
        let common = intersection(&mut model, a, nothing).unwrap();
        assert_eq!(model.unique_subshapes(common, NodeKind::Solid).len(), 0);
        let cut = difference(&mut model, nothing, a).unwrap();
        assert_eq!(model.unique_subshapes(cut, NodeKind::Solid).len(), 0);
    }

    #[test]
    fn hole_through_a_box() {
        // A bar through a block: the difference is a block with a square hole.
        let mut model = Model::new();
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (3.0, 3.0, 1.0));
        let bar = cuboid(&mut model, (1.0, 1.0, -1.0), (1.0, 1.0, 3.0));
        let result = difference(&mut model, block, bar).unwrap();
        check(&model, result, 8.0, 2.0 * 8.0 + 12.0 + 4.0, 1);
    }

    /// A block whose top is a face of a surface carried on far past it and folding back over
    /// it (a hair's breadth above), fused with a box through its top: points on the top, and
    /// crossings of it, are found on the top itself, not on the fold passing just above.
    #[test]
    fn a_face_of_a_surface_folding_back_over_it() {
        use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
        let mut model = Model::new();
        let (gap, w) = (0.002, 4.0);
        // Along u: the top from x 0 to 10 at z 1, then back over it at z 1 + gap; across v,
        // y 0 to w.
        let row = |y: f64| vec![Point::new(0.0, y, 1.0), Point::new(10.0, y, 1.0), Point::new(10.0, y, 1.0 + gap), Point::new(0.0, y, 1.0 + gap)];
        let (r0, r1) = (row(0.0), row(w));
        let poles: Vec<Vec<Point>> = r0.iter().zip(&r1).map(|(&a, &b)| vec![a, b]).collect();
        let u = KnotVector::new(vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], 1);
        let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let surface = Surface::BSpline(BSplineSurface::new(poles, u, v));
        let top = crate::make::face_on_domain(&mut model, surface, 0.0, 1.0, 0.0, 1.0).unwrap();
        // The block under it: the top's edges extruded down to z 0, closed by the bottom.
        let sides = crate::sweep::extrude(&mut model, top, Vector::new(0.0, 0.0, -1.0)).unwrap().shape;
        let block = sides;
        let volume = properties::volume(&model, block).mass.abs();
        assert!((volume - 40.0).abs() < 1e-9, "{volume}");
        let tower = cuboid(&mut model, (4.0, 1.0, 0.5), (2.0, 2.0, 1.0));
        let fused = union(&mut model, block, tower).unwrap();
        let report = validity::check(&model, fused);
        assert!(report.is_valid(), "{report:?}");
        let v = properties::volume(&model, fused).mass;
        assert!((v - 42.0).abs() < 1e-6, "{v}");
    }

    /// A rod cut back to a torus's side and fused to it again: the rod's end face lies on
    /// the torus, bounded by traced sections that stray from it between their fit points by
    /// more than their tolerance; they lie on the torus's face as far as they stray from
    /// their own (the same surface), which they are split along (else the rod is lost, or
    /// the shell left open).
    #[test]
    fn a_piece_cut_to_a_surface_fuses_back_onto_it() {
        let mut model = Model::new();
        let torus = crate::primitives::Torus::new(10.0, 4.0).build(&mut model).unwrap();
        let frame = Frame::new(Point::new(-8.0, 0.0, 0.0), -Direction::X, Direction::Y);
        let rod = Cylinder::new(2.0, 10.0).with_frame(frame).build(&mut model).unwrap();
        let piece = difference(&mut model, rod, torus).unwrap();
        let fused = union(&mut model, torus, piece).unwrap();
        let report = validity::check(&model, fused);
        assert!(report.is_valid(), "{report:?}");
        let expected = properties::volume(&model, torus).mass + properties::volume(&model, piece).mass;
        let v = properties::volume(&model, fused).mass;
        assert!((v - expected).abs() < 1e-6 * expected, "{v} vs {expected}");
    }

    #[test]
    fn faces_far_from_the_other_operand_pass_through() {
        let mut model = Model::new();
        // An L of two blocks, cut at one end: the faces at the other end are far from the
        // cutter and stay the same faces.
        let long = cuboid(&mut model, (0.0, 0.0, 0.0), (10.0, 1.0, 1.0));
        let leg = cuboid(&mut model, (0.0, 1.0, 0.0), (1.0, 5.0, 1.0));
        let l = union(&mut model, long, leg).unwrap();
        let frame = Frame::new(Point::new(9.5, 0.5, -1.0), Direction::Z, Direction::X);
        let cutter = Cylinder::new(0.25, 3.0).with_frame(frame).build(&mut model).unwrap();
        let cut = difference(&mut model, l, cutter).unwrap();
        assert!(validity::check(&model, cut).is_valid());
        let expected = 10.0 + 5.0 - PI * 0.25 * 0.25;
        assert!((properties::volume(&model, cut).mass - expected).abs() < 1e-6);
        let nodes = |s: Shape, model: &Model| model.explore(s, NodeKind::Face).map(|f| f.node).collect::<crate::hashing::FastSet<_>>();
        let (before, after) = (nodes(l, &model), nodes(cut, &model));
        assert!(before.intersection(&after).count() >= 5, "the far faces are the same faces");
        // Solids apart: every face is far.
        let apart = cuboid(&mut model, (20.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let both = union(&mut model, cut, apart).unwrap();
        assert!(validity::check(&model, both).is_valid());
        assert!((properties::volume(&model, both).mass - expected - 1.0).abs() < 1e-6);
        assert_eq!(model.explore(both, NodeKind::Solid).count(), 2);
        let common = intersection(&mut model, cut, apart).unwrap();
        assert_eq!(properties::volume(&model, common).mass, 0.0);
    }

    #[test]
    fn faces_apart_from_the_other_operand_where_they_overlap_pass_through() {
        let mut model = Model::new();
        // A post through a block's top: the post's bottom, inside the block, is apart from
        // the block's faces, and the same face in the results (as its top).
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
        let post = cuboid(&mut model, (4.0, 4.0, 8.0), (2.0, 2.0, 4.0));
        let at = |model: &Model, s: Shape, z: f64| {
            model.explore(s, NodeKind::Face).map(|f| f.node).find(|&n| {
                let b = properties::bounding_box(model, Shape::from(crate::topology::FaceId::try_from(n).unwrap()));
                (b.min.z - z).abs() < 1e-5 && (b.max.z - z).abs() < 1e-5 && b.max.x - b.min.x < 3.0
            })
        };
        let (bottom, top) = (at(&model, post, 8.0).unwrap(), at(&model, post, 12.0).unwrap());
        // (Left out of the graph: the bottom found inside the block, the top beyond where
        // the operands' boxes overlap.)
        let graph = operands(&model, block, post, true).unwrap();
        let states = far::states(&graph).unwrap();
        let left_out = |n: crate::topology::NodeId| graph.far.iter().zip(&states).find(|(f, _)| crate::topology::NodeId::from(f.face) == n).map(|(f, s)| (f.apart, *s));
        assert_eq!(left_out(bottom), Some((true, State::Inside)));
        assert_eq!(left_out(top), Some((false, State::Outside)));
        let fused = union(&mut model, block, post).unwrap();
        check(&model, fused, 1008.0, 616.0, 1);
        assert_eq!(at(&model, fused, 12.0), Some(top));
        let common = intersection(&mut model, block, post).unwrap();
        check(&model, common, 8.0, 24.0, 1);
        assert_eq!(at(&model, common, 8.0), Some(bottom));
        let cut = difference(&mut model, block, post).unwrap();
        check(&model, cut, 992.0, 616.0, 1);
        let rest = difference(&mut model, post, block).unwrap();
        check(&model, rest, 8.0, 24.0, 1);
        assert_eq!(at(&model, rest, 12.0), Some(top));
    }

    /// Two struts of one radius ending at a node, square to each other, their axes meeting
    /// there: their sides meet in halves of ellipses crossing where the sides touch, at the
    /// ends of the struts. On the first strut's side the section runs across its seam from
    /// one of those points to the other, traced across the second onto the other branch
    /// (beyond the second strut's end): it ends at both, each passed only near.
    #[test]
    fn struts_meeting_square_at_their_ends_cross_at_the_points_their_sides_touch() {
        use crate::geometry::ParametricCurve as _;
        use crate::math::Point;
        let t = 0.5 * 100.0 * 2f64.sqrt() / (1.0 + 2.0 * 2f64.sqrt());
        let s = 100.0 / (1.0 + 2.0 * 2f64.sqrt());
        let node = Point::new(2.0 * t, 0.0, t);
        let mut model = Model::new();
        let a = strut(&mut model, Point::new(2.0 * t + s, 0.0, t), node, 5.0);
        let b = strut(&mut model, Point::new(2.0 * t, t, 0.0), node, 5.0);
        let mut graph = operands(&model, a, b, true).unwrap();
        classified(&model, &mut graph);
        let side = graph.faces.iter().position(|f| f.operand == 0 && matches!(f.surface, crate::geometry::surface::Surface::BSpline(_))).unwrap();
        let ends: Vec<Point> = graph.faces[side].sections.iter().flat_map(|&e| {
            let (c, e) = (graph.edges[e].curve.as_ref().unwrap(), &graph.edges[e]);
            [c.point(e.first), c.point(e.last)]
        }).collect();
        let h = 5.0 * std::f64::consts::FRAC_1_SQRT_2;
        for touching in [Point::new(2.0 * t, -h, t - h), Point::new(2.0 * t, h, t + h)] {
            assert!(ends.iter().any(|p| p.distance(touching) <= 1e-3), "no section ends at {touching:?}: {ends:?}");
        }
    }

    #[test]
    fn far_faces_vote_with_the_regions_they_are_joined_to() {
        let mut model = Model::new();
        // A post through a block's top: the block's top is split, its sides and bottom left
        // out of the graph, and joined to the top's regions across the top's edges.
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
        let post = cuboid(&mut model, (4.0, 4.0, 8.0), (2.0, 2.0, 4.0));
        let mut graph = operands(&model, block, post, true).unwrap();
        let (pieces, arcs_of, regions, far_states) = classified(&model, &mut graph);
        let far_states = far_states.unwrap();
        let roots = clean_components(&graph, &pieces, &arcs_of, &regions, Some(&far_states));
        assert_eq!(roots.len(), regions.len() + graph.far.len());
        let block_far: Vec<usize> = (0..graph.far.len()).filter(|&i| graph.far[i].operand == 0).collect();
        assert_eq!(block_far.len(), 5);
        // (All in the group of the top's region outside the post.)
        let top = (0..regions.len()).find(|&ri| graph.faces[regions[ri].0].operand == 0 && regions[ri].2 == State::Outside).unwrap();
        for i in block_far {
            assert_eq!(far_states[i], State::Outside);
            assert_eq!(roots[regions.len() + i], roots[top]);
        }
    }

    #[test]
    fn void_inside_a_box() {
        let mut model = Model::new();
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (3.0, 3.0, 3.0));
        let inner = cuboid(&mut model, (1.0, 1.0, 1.0), (1.0, 1.0, 1.0));
        let result = difference(&mut model, block, inner).unwrap();
        check(&model, result, 26.0, 54.0 + 6.0, 1);
        let solid = model.unique_subshapes(result, NodeKind::Solid)[0];
        assert_eq!(model.unique_subshapes(solid, NodeKind::Shell).len(), 2);
    }

    #[test]
    fn disjoint_and_touching_boxes() {
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let far = cuboid(&mut model, (2.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, far).unwrap();
        check(&model, fused, 2.0, 12.0, 2);
        let common_result = intersection(&mut model, a, far).unwrap();
        assert!(model.unique_subshapes(common_result, NodeKind::Solid).is_empty());
        // Sharing a whole face: the union is one 2 x 1 x 1 box.
        let next = cuboid(&mut model, (1.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, next).unwrap();
        check(&model, fused, 2.0, 10.0, 1);
        let cut_result = difference(&mut model, a, next).unwrap();
        check(&model, cut_result, 1.0, 6.0, 1);
    }

    #[test]
    fn identical_boxes() {
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let b = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, b).unwrap();
        check(&model, fused, 1.0, 6.0, 1);
        let common_result = intersection(&mut model, a, b).unwrap();
        check(&model, common_result, 1.0, 6.0, 1);
        let cut_result = difference(&mut model, a, b).unwrap();
        assert!(model.unique_subshapes(cut_result, NodeKind::Solid).is_empty());
    }

    #[test]
    fn cut_splits_a_bar_in_two() {
        let mut model = Model::new();
        let bar = cuboid(&mut model, (0.0, 0.0, 0.0), (3.0, 1.0, 1.0));
        let slab = cuboid(&mut model, (1.0, -1.0, -1.0), (1.0, 3.0, 3.0));
        let result = difference(&mut model, bar, slab).unwrap();
        check(&model, result, 2.0, 12.0, 2);
    }

    #[test]
    fn solids_touching_along_an_edge() {
        // Cubes meeting along an edge stay two solids with a shell each (sharing the edge);
        // filling in the other two corners of the 2 x 2 square makes one solid again.
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (2.0, 1.0, 1.0));
        let b = cuboid(&mut model, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0));
        let fused = union(&mut model, a, b).unwrap();
        check(&model, fused, 3.0, 14.0, 1);
        let c = cuboid(&mut model, (0.0, 1.0, 0.0), (1.0, 1.0, 1.0));
        let d = cuboid(&mut model, (1.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let corner_a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let left = union(&mut model, corner_a, b).unwrap();
        let both = union(&mut model, left, c).unwrap();
        let both = union(&mut model, both, d).unwrap();
        check(&model, both, 4.0, 16.0, 1);
        let touching = union(&mut model, corner_a, b).unwrap();
        check(&model, touching, 2.0, 12.0, 2);
    }

    #[test]
    fn tiles_laid_on_a_floor() {
        // A floor and a compound of 64 tiles sunk a little into it, apart from each other:
        // fused in two operations (the floor, then the tiles at once), one solid.
        let mut model = Model::new();
        let floor = cuboid(&mut model, (0.0, 0.0, 0.0), (16.0, 16.0, 1.0));
        let tiles = model.new_empty(NodeKind::Compound);
        for i in 0..8 {
            for j in 0..8 {
                let tile = cuboid(&mut model, (i as f64 * 2.0 + 0.25, j as f64 * 2.0 + 0.25, 0.5), (1.5, 1.5, 1.0));
                model.add(tiles, tile).unwrap();
            }
        }
        let fused = union(&mut model, floor, tiles).unwrap();
        check(&model, fused, 256.0 + 64.0 * 1.5 * 1.5 * 0.5, 2.0 * 256.0 + 4.0 * 16.0 + 64.0 * 4.0 * 1.5 * 0.5, 1);
    }

    #[test]
    fn cylinders_touching_along_a_line() {
        // Side by side, touching where their sides are tangent: the faces leave the line
        // they share in the same direction (told apart by the way they curve), and the
        // cylinders stay two solids, as on OCCT.
        let mut model = Model::new();
        let (r, h) = (13.25, 65.0);
        let a = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::X, r, h);
        let b = cylinder(&mut model, (0.0, 2.0 * r, 0.0), Direction::X, r, h);
        let fused = union(&mut model, a, b).unwrap();
        let pi = std::f64::consts::PI;
        check(&model, fused, 2.0 * pi * r * r * h, 4.0 * pi * r * (h + r), 2);
    }

    #[test]
    fn cylinder_grazing_a_box_edge() {
        // A cylinder poking 0.01 past a box's side, its end flush with the box's front: the
        // box's edges along the side dip into it between samples (OCCT gives these values).
        let mut model = Model::new();
        let bar = cuboid(&mut model, (-19.8, -54.5, -37.73), (12.4, 3.5, 47.73));
        let boss = cylinder(&mut model, (-12.0, -54.5, -24.7), Direction::Y, 4.61, 40.0);
        let fused = union(&mut model, bar, boss).unwrap();
        check(&model, fused, 4508.434770483955, 2661.863807650023, 1);
    }

    #[test]
    fn hole_crossing_an_edge_at_a_sample() {
        // A hole through a block's side: the side's edges cross it twice, once exactly where
        // they are sampled (the interval beside that point holds the other crossing).
        let mut model = Model::new();
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (32.0, 10.0, 5.0));
        let drill = cylinder(&mut model, (10.35, 1.2, -1.0), Direction::Z, 1.25, 7.0);
        let result = difference(&mut model, block, drill).unwrap();
        let (r, d) = (1.25f64, 1.2f64);
        let cut_off = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
        let removed = 5.0 * (PI * r * r - cut_off);
        let arc = 2.0 * r * (PI - (d / r).acos());
        let area = 2.0 * (320.0 + 160.0 + 50.0) - 2.0 * (PI * r * r - cut_off) + 5.0 * arc - 5.0 * 2.0 * (r * r - d * d).sqrt();
        check(&model, result, 1600.0 - removed, area, 1);
    }

    fn cylinder(model: &mut Model, origin: (f64, f64, f64), axis: Direction, radius: f64, height: f64) -> Shape {
        let frame = Frame::from_z(Point::new(origin.0, origin.1, origin.2), axis);
        Cylinder::new(radius, height).with_frame(frame).build(model).unwrap()
    }

    fn sphere(model: &mut Model, center: (f64, f64, f64), radius: f64) -> Shape {
        let frame = Frame::new(Point::new(center.0, center.1, center.2), Direction::Z, Direction::X);
        Sphere::new(radius).with_frame(frame).build(model).unwrap()
    }

    /// Like `check`, to a relative accuracy (for results with approximated geometry).
    fn check_near(model: &Model, shape: Shape, volume: f64, area: f64, relative: f64) {
        let report = validity::check(model, shape);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(model, shape).mass;
        let a = properties::area(model, shape).mass;
        assert!((v - volume).abs() <= relative * volume.abs().max(1.0), "volume {v} != {volume}");
        assert!((a - area).abs() <= relative * area.abs().max(1.0), "area {a} != {area}");
    }

    /// Union, intersection and difference of `a` and `b` are valid and add up: |A ∪ B| +
    /// |A ∩ B| = |A| + |B| and |A − B| = |A| − |A ∩ B|, for volumes and (union and
    /// intersection) areas.
    fn check_consistent(model: &mut Model, a: Shape, b: Shape, relative: f64) -> (f64, f64) {
        let volume = |m: &Model, s: Shape| properties::volume(m, s).mass;
        let area = |m: &Model, s: Shape| properties::area(m, s).mass;
        let (va, vb, aa, ab) = (volume(model, a), volume(model, b), area(model, a), area(model, b));
        let merged = union(model, a, b).unwrap();
        let both = intersection(model, a, b).unwrap();
        let rest = difference(model, a, b).unwrap();
        for s in [merged, both, rest] {
            let report = validity::check(model, s);
            assert!(report.problems.is_empty(), "{report:?}");
        }
        let (vu, vc, vr) = (volume(model, merged), volume(model, both), volume(model, rest));
        let close = |x: f64, y: f64| (x - y).abs() <= relative * y.abs().max(1.0);
        assert!(close(vu + vc, va + vb), "volumes: {vu} + {vc} != {va} + {vb}");
        assert!(close(vr, va - vc), "difference volume {vr} != {va} - {vc}");
        assert!(vc > 0.0);
        let (au, ac) = (area(model, merged), area(model, both));
        assert!(close(au + ac, aa + ab), "areas: {au} + {ac} != {aa} + {ab}");
        (vc, ac)
    }

    /// A strip swept from a short curve (its first pole doubled: no tangent there, as font
    /// outlines' curves have), crossing a sphere beside its seam: the section is where
    /// the strip crosses the sphere (two curves), nothing along the seam (points beyond the
    /// spline's end are not on the strip's surface, though along its normal at that end,
    /// which vanishes, they are no distance from it).
    #[test]
    fn section_of_a_swept_strip_and_a_sphere() {
        let mut model = Model::new();
        let ball = sphere(&mut model, (0.0, 0.0, 0.0), 50.0);
        // (A glyph's piece of outline, met projecting text onto a sphere.)
        let (a, b) = (Point::new(86.07779878153892, -437.6899883531731, -59.34253200971563), Point::new(86.08473793833254, -437.66262713896106, -59.53699334792927));
        let curve = crate::geometry::curve::Curve::Bezier(crate::geometry::curve::BezierCurve::new(vec![a, a, b]));
        let edge = crate::make::edge_from_curve(&mut model, curve).unwrap();
        let along = Vector::new(-0.18758591913397418, 0.9735687357550348, 0.13028982962229022) * 800.0;
        let strip = crate::sweep::extrude(&mut model, edge, along).unwrap().shape;
        let cut = section(&mut model, strip, ball).unwrap();
        let edges = model.unique_subshapes(cut, NodeKind::Edge);
        assert_eq!(edges.len(), 2, "{} section edges", edges.len());
    }

    /// A box and a slab under a wave touching the box's top along a line: the regions of
    /// the top and of the wave either side of the line, on the other operand only along it,
    /// are outside and inside it.
    #[test]
    fn box_and_a_wave_touching_its_top() {
        let mut model = Model::new();
        let block = cuboid(&mut model, (-5.0, -1.0, -2.5), (10.0, 2.0, 5.0));
        let wave = |x: f64| 0.25 + 0.75 * (PI * x / 6.0).cos();
        let points: Vec<Point> = (0..=48).map(|i| -12.0 + 0.5 * i as f64).map(|x| Point::new(x, wave(x), -2.5)).collect();
        let curve = crate::geometry::curve::Curve::BSpline(crate::geometry::fit::interpolate(&points, false, 1e-9).unwrap());
        let (first, last) = (points[0], points[48]);
        let edges = [
            crate::make::edge_from_curve(&mut model, curve).unwrap(),
            crate::make::edge_between(&mut model, last, Point::new(12.0, -2.0, -2.5)).unwrap(),
            crate::make::edge_between(&mut model, Point::new(12.0, -2.0, -2.5), Point::new(-12.0, -2.0, -2.5)).unwrap(),
            crate::make::edge_between(&mut model, Point::new(-12.0, -2.0, -2.5), first).unwrap(),
        ];
        let profile = crate::make::wire_within(&mut model, &edges, 1e-7).unwrap();
        let face = crate::make::planar_face(&mut model, profile).unwrap();
        let slab = crate::sweep::extrude(&mut model, face, Vector::new(0.0, 0.0, 5.0)).unwrap().shape;
        let result = intersection(&mut model, block, slab).unwrap();
        assert!(validity::check(&model, result).is_valid());
        // Over -5 <= x <= 5, from y = -1 up to the wave, 5 deep.
        let volume = 5.0 * (12.5 + 0.75 * 12.0 / PI * (5.0 * PI / 6.0).sin());
        let v = properties::volume(&model, result).mass;
        assert!((v - volume).abs() <= 1e-4 * volume, "volume {v} != {volume}");
    }

    /// Two triangular bars twisted 40° either way from one triangle (a knurl's): where they
    /// part, each one's faces meet the other's from a corner of the triangle to where the
    /// bars' edges cross, along curves leaving and rejoining the edges they end on (found
    /// past them a hair at both ends, inside between): the overlap they bound closes, and
    /// the overlap and the rest of each make it up.
    #[test]
    fn twisted_bars_from_one_start() {
        twisted_bars(include_str!("tests/twisted_bars_from_one_start_a.brep"), include_str!("tests/twisted_bars_from_one_start_b.brep"), 2.5821);
    }

    /// [`twisted_bars_from_one_start`], the bars 14 long (the knurl's, twisting more
    /// slowly): the curve from a corner to where the edges cross is a step long, and is
    /// followed up to the faces' bounds in shorter ones.
    #[test]
    fn long_twisted_bars_from_one_start() {
        twisted_bars(include_str!("tests/twisted_bars_long_a.brep"), include_str!("tests/twisted_bars_long_b.brep"), 12.0498);
    }

    /// Two bars meeting where they start: the overlap, the rest of the first and the union
    /// sound, adding up, the union as OCCT's (`occt`).
    fn twisted_bars(a: &str, b: &str, occt: f64) {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let a = read(a);
        let b = read(b);
        let whole = properties::volume(&model, a).mass;
        let common = intersection(&mut model, a, b).unwrap();
        let cut = difference(&mut model, a, b).unwrap();
        let fused = union(&mut model, a, b).unwrap();
        for r in [common, cut, fused] {
            assert!(validity::check(&model, r).is_valid(), "{:?}", validity::check(&model, r));
        }
        let (c, d, u) = (properties::volume(&model, common).mass, properties::volume(&model, cut).mass, properties::volume(&model, fused).mass);
        assert!((c + d - whole).abs() <= 1e-6 * whole, "{c} + {d} != {whole}");
        assert!((u - (2.0 * whole - c)).abs() <= 1e-6 * u, "{u}");
        assert!((u - occt).abs() <= 2e-4 * occt, "{u}");
    }

    /// The top of a threaded container (a cylinder face whose band winds round four times
    /// up its thread, in one face) halved by a block through its axis, the block's face
    /// meeting the cylinder along its seam: the lines on the band where the block's face
    /// meets it go each in its own turn, so the halves part (half the part either side).
    #[test]
    fn thread_band_halved_through_its_seam() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let part = read(include_str!("tests/thread_band_halved_a.brep"));
        let block = read(include_str!("tests/thread_band_halved_b.brep"));
        let whole = properties::volume(&model, part).mass;
        let cut = difference(&mut model, part, block).unwrap();
        let common = intersection(&mut model, part, block).unwrap();
        for half in [cut, common] {
            assert!(validity::check(&model, half).is_valid(), "{:?}", validity::check(&model, half));
        }
        let (a, b) = (properties::volume(&model, cut).mass, properties::volume(&model, common).mass);
        assert!((a + b - whole).abs() <= 1e-6 * whole, "{a} + {b} != {whole}");
        assert!((a - b).abs() <= 1e-3 * whole, "{a} against {b}");
    }

    /// A slot's end, its rounded corners (cylinders of radius 0.5) chamfered at the rim, cut
    /// by a longer slot of the same width: the second's side is tangent to the first's
    /// corners along their lines, and the chamfer's edge on a corner leaves that line within
    /// tolerance of the side for a while, met where it leaves; that stretch of it closes the
    /// line the side is split along, and the regions of the first's faces either side of it
    /// lie on either side of the second (OCCT's volume).
    #[test]
    fn slot_cut_past_a_rounded_slot_end() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let part = read(include_str!("tests/slot_past_a_rounded_end_a.brep"));
        let slot = read(include_str!("tests/slot_past_a_rounded_end_b.brep"));
        let result = difference(&mut model, part, slot).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        let v = properties::volume(&model, result).mass;
        assert!((v - 12.481824).abs() <= 1e-5 * 12.481824, "volume {v}");
    }

    /// Two uneven stucco blocks (lofts, from cqfantasy's round block generator, Apache-2.0)
    /// side by side, touching along a face, their top B-spline faces meeting nearly tangent
    /// along an edge each block has: a section between those faces ends on one block's edge,
    /// and the other's, lying along it, is split there too (left whole, it did not match the
    /// pieces, and the shell stayed open). The union is both blocks.
    #[test]
    fn stucco_blocks_side_by_side() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let a = read(include_str!("tests/stucco_blocks_side_by_side_a.brep"));
        let b = read(include_str!("tests/stucco_blocks_side_by_side_b.brep"));
        let both = properties::volume(&model, a).mass + properties::volume(&model, b).mass;
        let result = union(&mut model, a, b).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        assert_eq!(model.explore(result, NodeKind::Solid).count(), 1);
        let v = properties::volume(&model, result).mass;
        assert!((v - both).abs() <= 1e-6 * both, "volume {v} against {both}");
    }

    /// A strip of rays from the origin through a glyph's quadratic edge (the sweep a conical
    /// projection sections a sphere with), 1.4 wide and 90 long, and a sphere of radius 50
    /// about the origin: one arc across it, 1.081 long. (Sampled by spans alone, the strip's
    /// grid has cells 22 long and 0.25 wide; the foot of a point on its side edge, found from
    /// the nearest sample 8.8 away, settled 1.0 off, and the arc was lost.)
    #[test]
    fn a_long_ray_strip_sectioned_by_a_sphere() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let strip = read(include_str!("tests/ray_strip_on_sphere_a.brep"));
        let sphere = read(include_str!("tests/ray_strip_on_sphere_b.brep"));
        let result = section(&mut model, strip, sphere).unwrap();
        let edges = model.unique_subshapes(result, NodeKind::Edge);
        assert_eq!(edges.len(), 1);
        let length = properties::length(&model, result).mass;
        assert!((length - 1.081).abs() < 1e-3, "{length}");
    }

    /// A box and a slab under a wave (a B-spline) rising just through the box's top, by 2e-5
    /// over 0.013, and again by 5e-4: the top meets the wave in pairs of lines close
    /// together, and is kept between them.
    #[test]
    fn box_and_a_wave_just_through_its_top() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let block = read(include_str!("tests/box_under_a_wave_a.brep"));
        let slab = read(include_str!("tests/box_under_a_wave_b.brep"));
        let result = intersection(&mut model, block, slab).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        // (OCCT's.)
        let v = properties::volume(&model, result).mass;
        assert!((v - 187.5).abs() <= 1e-4 * 187.5, "volume {v}");
        let top = model.explore(result, NodeKind::Face).filter(|&f| (properties::area(&model, f).center_of_mass.y - 1.0).abs() < 1e-9).count();
        assert_eq!(top, 2);
    }

    /// A cylinder (an extruded circle) cut through by a slab 0.001 thick tilted 20°: the strip
    /// of its side within the slab is narrower than the sag of its boundaries' polylines, its
    /// interior point taken from denser ones; it goes, leaving two valid solids (OCCT's
    /// volumes).
    #[test]
    fn a_cylinder_cut_by_a_thin_tilted_slab() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let cylinder = read(include_str!("tests/cylinder_cut_by_thin_slab_a.brep"));
        let slab = read(include_str!("tests/cylinder_cut_by_thin_slab_b.brep"));
        let result = difference(&mut model, cylinder, slab).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        let mut volumes: Vec<f64> = model.explore(result, NodeKind::Solid).map(|s| properties::volume(&model, s).mass).collect();
        volumes.sort_by(f64::total_cmp);
        assert_eq!(volumes.len(), 2, "{volumes:?}");
        for (v, occt) in volumes.iter().zip([45363.99, 55230.79]) {
            assert!((v - occt).abs() <= 1e-5 * occt, "{volumes:?}");
        }
    }

    /// A ring with a thread fused into its bore (the bore one face, a band winding twice
    /// round, its parameters spanning two periods), cut by a disc with a chamfered top at
    /// the bottom of the bore and six screw holes from the top: the circle where the disc's
    /// top meets the bore goes in the period the bore is there, splitting it, and rays from
    /// the holes' bottoms meet the bore in whichever period (else the bore's bottom stays,
    /// or a hole's bottom goes, the shell open).
    #[test]
    fn a_band_winding_round_a_cylinder_cut_across() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let ring = read(include_str!("tests/threaded_ring_a.brep"));
        let tools = read(include_str!("tests/threaded_ring_b.brep"));
        // (The disc takes the ring from r 75 to its side, r 76.5 at its foot (0.0015 high)
        // narrowing to 75.0015 at its top, z 1.5: π ∫ (R² - 75²) dz; the holes are r 1.55,
        // 15 deep.)
        let cone = 75.0 * (1.5f64.powi(2) - 0.0015f64.powi(2)) + (1.5f64.powi(3) - 0.0015f64.powi(3)) / 3.0;
        let expected = std::f64::consts::PI * (cone + (76.5f64.powi(2) - 75.0f64.powi(2)) * 0.0015 + 6.0 * 1.55f64.powi(2) * 15.0);
        let common = intersection(&mut model, ring, tools).unwrap();
        let v = properties::volume(&model, common).mass;
        assert!((v - expected).abs() < 1e-9 * expected, "{v} vs {expected}");
        let volume = properties::volume(&model, ring).mass;
        let result = difference(&mut model, ring, tools).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        // (The thread's twisted faces integrated to a few parts in 10⁵.)
        let removed = volume - properties::volume(&model, result).mass;
        assert!((removed - expected).abs() < 1e-4 * volume, "{removed} vs {expected}");
    }

    /// A thread whose flank is one pitch long, so that its turns touch along a helix, cut
    /// across by a box: the box's top and bottom meet the flank in circles closing where the
    /// turns touch, there ending on the edges of two turns. The same away from the origin: a
    /// point there lies on one turn, the nearest of the surface's samples on the next (once
    /// that turn's foot, 5e-5 off, was taken: the circle lost its second end and the thread
    /// was all outside the box).
    #[test]
    fn a_thread_whose_turns_touch_cut_across_away_from_the_origin() {
        let common = |z0: f64| {
            let mut model = Model::new();
            // (As pykentron's `makeHelix`: a spline through points and tangents, 32 a turn,
            // by the length along it; the ISO thread's internal profile, swept in the
            // helix's Frenet frame.)
            let (r, turns) = (3.0, 5.0);
            let speed = (2.0 * PI).hypot(1.0);
            let steps = (turns * 32.0) as usize;
            let params: Vec<f64> = (0..=steps).map(|i| turns * speed * i as f64 / steps as f64).collect();
            let (points, tangents): (Vec<Point>, Vec<Option<Vector>>) = params
                .iter()
                .map(|&s| {
                    let (u, du, dv) = (2.0 * PI * s / speed, 2.0 * PI / speed, 1.0 / speed);
                    (Point::new(r * u.cos(), r * u.sin(), z0 + s / speed), Some(Vector::new(-r * u.sin() * du, r * u.cos() * du, dv)))
                })
                .collect();
            let helix = crate::geometry::fit::interpolate_with_tangents(&points, Some(&params), &tangents, false, false, 1e-9).unwrap();
            let helix = crate::make::edge_from_curve(&mut model, crate::geometry::curve::Curve::BSpline(helix)).unwrap();
            let spine = crate::make::wire(&mut model, &[helix]).unwrap();
            let depth = 5.0 / 8.0 * 3f64.sqrt() / 2.0;
            let corners = [Point::new(r, 0.0, z0 - 0.5), Point::new(r, 0.0, z0 + 0.5), Point::new(r - depth, 0.0, z0)];
            let sides: Vec<Shape> = (0..3).map(|i| crate::make::edge_between(&mut model, corners[i], corners[(i + 1) % 3]).unwrap()).collect();
            let outline = crate::make::wire(&mut model, &sides).unwrap();
            let profile = crate::make::planar_face(&mut model, outline).unwrap();
            let options = crate::sweep::PipeOptions { frenet: true, ..Default::default() };
            let thread = crate::sweep::pipe_with(&mut model, profile, spine, &options).unwrap().shape;
            let tool = Cuboid::from_corners(Point::new(-6.5, -6.5, z0 + 1.0), Point::new(6.5, 6.5, z0 + turns - 1.0)).build(&mut model).unwrap();
            let common = intersection(&mut model, thread, tool).unwrap();
            assert!(validity::check(&model, common).is_valid(), "{:?}", validity::check(&model, common));
            properties::volume(&model, common).mass
        };
        // (Three of the five turns' worth: the thread's cross-section, 0.2706, about the
        // helix, 2π 2.9 round at its centroid.)
        let (near, far) = (common(-1.0), common(50.0));
        let expected = 3.0 * 0.5 * 5.0 / 8.0 * 3f64.sqrt() / 2.0 * (2.0 * PI * (3.0 - 5.0 / 8.0 * 3f64.sqrt() / 2.0 / 3.0)).hypot(1.0);
        assert!((near - expected).abs() < 2e-3 * expected, "{near} vs {expected}");
        assert!((far - near).abs() < 1e-6 * near, "{far} vs {near}");
    }

    /// A socket with a thread fused into its bore (the bore's ceiling at z 5.37 meeting the
    /// thread's swept faces in fitted sections, straying 2e-6 from them, more than their
    /// tolerance), and a box 20 x 20 from z 5.27, just under the ceiling, its corners in the
    /// thread: where the box's sides meet those sections the sections along the sides pass,
    /// so the sides are split there, and the box's bottom, in the bore, is outside the
    /// socket (else a slab 0.1 thick is taken as inside it: 1280).
    #[test]
    fn a_box_reaching_just_into_a_threaded_bore() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let socket = read(include_str!("tests/threaded_socket_a.brep"));
        let block = read(include_str!("tests/threaded_socket_b.brep"));
        let volume = properties::volume(&model, socket).mass;
        // (The socket above the ceiling within the box, 20 x 20 x 3.1, and slivers of the
        // thread below it.)
        let common = intersection(&mut model, socket, block).unwrap();
        assert!(validity::check(&model, common).is_valid(), "{:?}", validity::check(&model, common));
        let v = properties::volume(&model, common).mass;
        assert!(v > 1240.0 && v < 1241.0, "{v}");
        let result = difference(&mut model, socket, block).unwrap();
        assert!(validity::check(&model, result).is_valid(), "{:?}", validity::check(&model, result));
        let left = properties::volume(&model, result).mass;
        assert!((left + v - volume).abs() < 1e-6 * volume, "{left} + {v} vs {volume}");
    }

    /// A tube from `a` to `b` of radius `r`, as a lattice's strut: two circles square to the
    /// line swept through along it (its side a B-spline closed round without a period).
    fn strut(model: &mut Model, a: crate::math::Point, b: crate::math::Point, r: f64) -> Shape {
        use crate::geometry::curve::{Circle, Curve};
        use crate::math::{Direction, Frame};
        let axis = Direction::from_vector(b - a);
        let spine = crate::make::edge_between(model, a, b).unwrap();
        let circles: Vec<(Shape, crate::sweep::PipeOptions)> = [a, b]
            .iter()
            .map(|&c| {
                let edge = crate::make::edge_from_curve(model, Curve::Circle(Circle::new(Frame::from_z(c, axis), r))).unwrap();
                (crate::make::wire(model, &[edge]).unwrap(), Default::default())
            })
            .collect();
        crate::sweep::pipe_through(model, &circles, spine, true).unwrap()
    }

    /// A shelled lid, its top's edges rounded, cut by two countersunk holes where the rounds
    /// begin: the countersinks' rims cross the lines where the rounds meet the top tangent
    /// (one on a straight round, one where they turn a corner), and the sections of the
    /// rounds with them meet those lines tangent or stop short of them. The result is closed,
    /// with OCCT's volume (0.169687).
    #[test]
    fn countersunk_holes_over_rounded_edges() {
        use crate::offset::Join;
        let mut model = Model::new();
        let at = |x: f64, y: f64| Point::new(x, y, 0.25);
        let corners = [at(1.1, 0.75), at(1.65, -0.75), at(-1.65, -0.75)];
        let edges: Vec<Shape> = (0..3).map(|k| crate::make::edge_between(&mut model, corners[k], corners[(k + 1) % 3]).unwrap()).collect();
        let wire = crate::make::wire(&mut model, &edges).unwrap();
        let profile = crate::make::planar_face(&mut model, wire).unwrap();
        let prism = crate::sweep::extrude(&mut model, profile, Vector::new(0.0, 0.0, 0.5)).unwrap().shape;
        // Edges whose vertices all satisfy `test`.
        let edges_where = |model: &Model, shape: Shape, test: &dyn Fn(&[Point]) -> bool| -> Vec<Shape> {
            model
                .unique_subshapes(shape, NodeKind::Edge)
                .into_iter()
                .filter(|&e| test(&model.unique_subshapes(e, NodeKind::Vertex).into_iter().map(|v| model.point(v)).collect::<Vec<_>>()))
                .collect()
        };
        let upright = edges_where(&model, prism, &|p| p.len() == 2 && (p[0].x - p[1].x).abs() < 1e-9 && (p[0].y - p[1].y).abs() < 1e-9);
        let rounded = crate::blend::fillet(&mut model, prism, &upright, 0.5).unwrap();
        let top = edges_where(&model, rounded, &|p| p.iter().all(|p| (p.z - 0.75).abs() < 1e-9));
        let rounded = crate::blend::fillet(&mut model, rounded, &top, 0.25).unwrap();
        let underside = model
            .unique_subshapes(rounded, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(crate::geometry::surface::Surface::Plane(p)) if (p.origin().z - 0.25).abs() < 1e-9))
            .unwrap();
        let mut lid = crate::offset::thick_solid(&mut model, rounded, &[underside], -0.06, Join::Arc).unwrap();
        for y in [-0.666854, 0.333146] {
            let bore = cylinder(&mut model, (0.655871, y, -1.378836), Direction::Z, 0.15, 1.378836 + 0.687068);
            let sink = cone(&mut model, (0.655871, y, 0.687068), Direction::Z, 0.15, 0.225, 0.75 - 0.687068);
            let tool = union(&mut model, bore, sink).unwrap();
            lid = difference(&mut model, lid, tool).unwrap();
        }
        let report = validity::check(&model, lid);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(&model, lid).mass;
        assert!((v - 0.169687).abs() < 1e-5, "volume {v}");
    }

    /// Struts meeting end to end at a point of a lattice (an octagon's and a square's of a
    /// truncated cuboctahedron), their sides touching at the points their ends' circles
    /// share and crossing each other's seams: the sections end on those points and go on
    /// across the seams, and the overlap is as OCCT's (a common of 69.04, 166.67, 96.22).
    #[test]
    fn struts_meeting_at_a_lattice_node() {
        use crate::math::Point;
        let t = 0.5 * 100.0 * 2f64.sqrt() / (1.0 + 2.0 * 2f64.sqrt());
        let s = 100.0 / (1.0 + 2.0 * 2f64.sqrt());
        let node = Point::new(2.0 * t, 0.0, t);
        let cases = [
            ((Point::new(t, 0.0, 2.0 * t), Point::new(t, 0.0, 2.0 * t + s)), (Point::new(t, 0.0, 2.0 * t + s), Point::new(2.0 * t, 0.0, 100.0 - t)), 69.0356),
            ((Point::new(2.0 * t + s, 0.0, t), node), (Point::new(2.0 * t, t, 0.0), node), 166.6667),
            ((node, Point::new(t, 0.0, 2.0 * t)), (Point::new(2.0 * t, t, 0.0), node), 96.2250),
        ];
        for ((a0, a1), (b0, b1), occt) in cases {
            let mut model = Model::new();
            let a = strut(&mut model, a0, a1, 5.0);
            let b = strut(&mut model, b0, b1, 5.0);
            let common = intersection(&mut model, a, b).unwrap();
            assert!(validity::check(&model, common).is_valid(), "{:?}", validity::check(&model, common));
            let v = properties::volume(&model, common).mass;
            assert!((v - occt).abs() <= 1e-3 * occt, "{v} vs {occt}");
        }
    }

    /// A cylinder and a strip winding twenty times round it (a B-spline ruled between two
    /// helices in 63 spans, as a thread's flank, from inside it to outside): they meet along
    /// a curve winding with it, found all along (the strip's points found on the turn they
    /// are on) and split where it crosses the cylinder's seam, each turn.
    #[test]
    fn a_cylinder_and_a_strip_winding_round_it() {
        cylinder_and_winding_strip(20.0, 2.2);
    }

    /// As [`a_cylinder_and_a_strip_winding_round_it`], winding five times (a short thread's
    /// flank: its turns as close, though fewer) from just inside the cylinder.
    #[test]
    fn a_cylinder_and_a_strip_winding_round_it_a_few_times() {
        cylinder_and_winding_strip(5.0, 2.499);
    }

    /// A cylinder and a strip winding twelve times round it, from the side away from the
    /// cylinder's seam: points a whole turn apart along the strip are alike, so its box is
    /// found from points in each of its spans (from twelve along it it would be a sliver
    /// on that side, clear of the seam, whose crossings would then be missed). The curve
    /// where they meet is split at each crossing of the seam and follows the strip.
    #[test]
    fn a_cylinder_and_a_strip_winding_round_it_whole_turns_from_its_far_side() {
        use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
        use crate::geometry::ParametricCurve as _;
        let (turns, n, r) = (12.0, 241, 2.5);
        let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let helix = |radius: f64| {
            let at = |t: f64| 2.0 * PI * turns * t + PI;
            let points: Vec<Point> = params.iter().map(|&t| Point::new(radius * at(t).cos(), radius * at(t).sin(), turns * t)).collect();
            crate::geometry::fit::interpolate_with_parameters(&points, &params).unwrap()
        };
        let (inner, outer) = (helix(2.2), helix(3.0));
        let poles: Vec<Vec<Point>> = inner.poles().iter().zip(outer.poles()).map(|(&a, &b)| vec![a, b]).collect();
        let u = KnotVector::new(inner.knots().to_vec(), inner.multiplicities().to_vec(), inner.degree());
        let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let mut model = Model::new();
        let strip = crate::make::face_on_surface(&mut model, Surface::BSpline(BSplineSurface::new(poles, u, v))).unwrap();
        let bar = cylinder(&mut model, (0.0, 0.0, -1.0), Direction::Z, r, turns + 2.0);
        let cut = section(&mut model, bar, strip).unwrap();
        let edges = model.unique_subshapes(cut, NodeKind::Edge);
        assert!(edges.len() > turns as usize && edges.len() <= 2 * turns as usize + 2, "{} edges", edges.len());
        // (On the cylinder all along, between the points fitted through as well.)
        for e in edges {
            let (curve, first, last) = model.edge_curve(e).unwrap();
            let far = (0..=400).map(|k| { let p = curve.point(first + (last - first) * k as f64 / 400.0); (p.x.hypot(p.y) - r).abs() }).fold(0.0, f64::max);
            assert!(far <= 1e-4, "a section edge {far:.2e} off the cylinder");
        }
    }

    /// A cylinder of radius 2.5 and a strip winding `turns` times round it from radius `inner`
    /// to 3.
    fn cylinder_and_winding_strip(turns: f64, inner: f64) {
        use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
        let (n, r) = (64, 2.5);
        let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let helix = |radius: f64| {
            let points: Vec<Point> = params.iter().map(|&t| Point::new(radius * (2.0 * PI * turns * t).cos(), radius * (2.0 * PI * turns * t).sin(), turns * t)).collect();
            crate::geometry::fit::interpolate_with_parameters(&points, &params).unwrap()
        };
        let (inner, outer) = (helix(inner), helix(3.0));
        let poles: Vec<Vec<Point>> = inner.poles().iter().zip(outer.poles()).map(|(&a, &b)| vec![a, b]).collect();
        let u = KnotVector::new(inner.knots().to_vec(), inner.multiplicities().to_vec(), inner.degree());
        let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let surface = BSplineSurface::new(poles, u, v);
        // The curve where the strip is at radius r, found finely: at each of many u where the
        // strip reaches across it (not all: its interpolated edges wobble), by bisection along
        // v (outwards along the strip's rulings).
        let fine: f64 = {
            use crate::geometry::ParametricSurface;
            let at = |u: f64| {
                let radius = |v: f64| {
                    let p = surface.point(u, v);
                    p.x.hypot(p.y) - r
                };
                if radius(0.0) >= 0.0 || radius(1.0) <= 0.0 {
                    return None;
                }
                let (mut lo, mut hi) = (0.0, 1.0);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if radius(mid) < 0.0 { lo = mid } else { hi = mid }
                }
                Some(surface.point(u, 0.5 * (lo + hi)))
            };
            let points: Vec<Option<Point>> = (0..=20000).map(|i| at(i as f64 / 20000.0)).collect();
            points.windows(2).filter_map(|w| Some(w[0]?.distance(w[1]?))).sum()
        };
        let mut model = Model::new();
        let strip = crate::make::face_on_surface(&mut model, Surface::BSpline(surface)).unwrap();
        let bar = cylinder(&mut model, (0.0, 0.0, -1.0), Direction::Z, r, turns + 2.0);
        let cut = section(&mut model, bar, strip).unwrap();
        let edges = model.unique_subshapes(cut, NodeKind::Edge);
        let length = properties::length(&model, cut).mass;
        assert!((length - fine).abs() <= 5e-4 * fine, "length {length} != {fine}");
        // (Split at each of its crossings of the cylinder's seam, once a turn.)
        assert!(edges.len() > turns as usize && edges.len() <= 2 * turns as usize + 2, "{} edges", edges.len());
    }

    /// A tube (radii 3 and 4, 4.875 high) and a thread on it (as CadQuery's contrib example
    /// makes one: four faces ruled between helices, two turns, fading out to its root at both
    /// ends), sunk 0.001 into the tube's outer wall. Where the thread fades out its flanks
    /// meet the wall along curves ending on its crest, which there crosses the wall in a
    /// section 0.0002 long, between points where the crest's edges meet it: the union is
    /// whole, nearly the sum of the two.
    #[test]
    fn a_tube_and_a_thread_fading_out_into_its_wall() {
        use crate::geometry::curve::Curve;
        let (radius, pitch, height, d, depth, aspect, frac) = (4.0 - 1e-3, 2.0, 4.0, 0.5, 0.5, 10.0, 0.1);
        let helix = |r_eps: f64, d: f64| {
            move |t: f64| {
                let (z, r) = if t > frac && t < 1.0 - frac {
                    (height * t + d, radius + r_eps)
                } else if t <= frac {
                    (height * t + d * (PI / 2.0 * t / frac).sin(), radius + r_eps * (PI / 2.0 * t / frac).sin())
                } else {
                    (height * t - d * (2.0 * PI - PI / 2.0 * (1.0 - t) / frac).sin(), radius - r_eps * (2.0 * PI - PI / 2.0 * (1.0 - t) / frac).sin())
                };
                let a = 2.0 * PI * height / pitch * t;
                Point::new(-r * a.sin(), r * a.cos(), z)
            }
        };
        let mut model = Model::new();
        let n = 200;
        let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let mut wire = |model: &mut Model, f: &dyn Fn(f64) -> Point| {
            let points: Vec<Point> = params.iter().map(|&t| f(t)).collect();
            let curve = crate::geometry::fit::interpolate_with_parameters(&points, &params).unwrap();
            let edge = crate::make::edge_from_curve(model, Curve::BSpline(curve)).unwrap();
            crate::make::wire(model, &[edge]).unwrap()
        };
        let (e1_bottom, e1_top) = (wire(&mut model, &helix(0.0, -d)), wire(&mut model, &helix(0.0, d)));
        let (e2_bottom, e2_top) = (wire(&mut model, &helix(depth, -d / aspect)), wire(&mut model, &helix(depth, d / aspect)));
        let ruled = crate::loft::LoftOptions { ruled: true, ..Default::default() };
        let mut faces = Vec::new();
        for pair in [[e1_bottom, e1_top], [e2_bottom, e2_top], [e1_bottom, e2_bottom], [e1_top, e2_top]] {
            let shell = crate::loft::loft(&mut model, &pair, &ruled).unwrap().shape;
            faces.extend(model.unique_subshapes(shell, NodeKind::Face));
        }
        let shell = crate::sew::sew(&mut model, &faces, &Default::default()).unwrap();
        let shell = model.unique_subshapes(shell, NodeKind::Shell)[0];
        let mut thread = model.new_empty(NodeKind::Solid);
        model.add(thread, shell).unwrap();
        if properties::volume(&model, thread).mass < 0.0 {
            thread = model.new_empty(NodeKind::Solid);
            model.add(thread, shell.reversed()).unwrap();
        }
        let outer = cylinder(&mut model, (0.0, 0.0, -d), Direction::Z, 4.0, height + 1.75 * d);
        let inner = cylinder(&mut model, (0.0, 0.0, -d), Direction::Z, 3.0, height + 1.75 * d);
        let tube = difference(&mut model, outer, inner).unwrap();
        let result = union(&mut model, tube, thread).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let sum = properties::volume(&model, tube).mass + properties::volume(&model, thread).mass;
        let volume = properties::volume(&model, result).mass;
        assert!(volume < sum && volume > sum - 0.1, "{volume} against {sum}");
    }

    #[test]
    fn box_with_a_cylindrical_hole() {
        let mut model = Model::new();
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (4.0, 4.0, 2.0));
        let drill = cylinder(&mut model, (2.0, 2.0, -1.0), Direction::Z, 1.0, 4.0);
        let result = difference(&mut model, block, drill).unwrap();
        check(&model, result, 32.0 - 2.0 * PI, 64.0 + 2.0 * PI, 1);
        let plug = intersection(&mut model, block, drill).unwrap();
        check(&model, plug, 2.0 * PI, 6.0 * PI, 1);
    }

    #[test]
    fn half_a_cylinder() {
        // A box cuts a cylinder along its axis (the plane meets it in two lines).
        let mut model = Model::new();
        let rod = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 1.0, 2.0);
        let half = cuboid(&mut model, (0.0, -2.0, -1.0), (2.0, 4.0, 4.0));
        let result = difference(&mut model, rod, half).unwrap();
        check(&model, result, PI, 3.0 * PI + 4.0, 1);
    }

    #[test]
    fn oblique_cut_through_a_cylinder() {
        // A tilted half-space keeps the cylinder below z = 2 + x·tan α: the cap is an
        // ellipse; by symmetry the volume is that of height 2.
        let mut model = Model::new();
        let rod = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 1.0, 4.0);
        let alpha: f64 = 0.4;
        let slab = cuboid(&mut model, (-5.0, -5.0, -10.0), (10.0, 10.0, 10.0));
        let tilt = crate::math::Transform::translation(Vector::new(0.0, 0.0, 2.0))
            * crate::math::Transform::rotation(crate::math::Axis::Y, -alpha);
        let below = intersection(&mut model, rod, slab.moved(&tilt)).unwrap();
        // The cap's pcurve on the cylinder (a sine wave) is approximated.
        check_near(&model, below, 2.0 * PI, 4.0 * PI + PI + PI / alpha.cos(), 1e-8);
    }

    #[test]
    fn hemisphere() {
        let mut model = Model::new();
        let ball = sphere(&mut model, (0.0, 0.0, 0.0), 1.0);
        let upper = cuboid(&mut model, (-2.0, -2.0, 0.0), (4.0, 4.0, 2.0));
        let dome = intersection(&mut model, ball, upper).unwrap();
        check_near(&model, dome, 2.0 * PI / 3.0, 3.0 * PI, 1e-9);
        let bowl = difference(&mut model, ball, upper).unwrap();
        check_near(&model, bowl, 2.0 * PI / 3.0, 3.0 * PI, 1e-9);
    }

    #[test]
    fn sphere_and_box_anywhere() {
        let mut model = Model::new();
        let ball = sphere(&mut model, (0.3, 0.2, 0.1), 1.0);
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        check_consistent(&mut model, ball, block, 1e-9);
    }

    #[test]
    fn parallel_cylinders() {
        // Two unit cylinders whose axes are 1 apart overlap in a lens.
        let mut model = Model::new();
        let a = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 1.0, 2.0);
        let b = cylinder(&mut model, (1.0, 0.0, -1.0), Direction::Z, 1.0, 4.0);
        let lens = 2.0 * (0.5f64).acos() - 0.5 * 3.0f64.sqrt();
        let (volume, _) = check_consistent(&mut model, a, b, 1e-9);
        assert!((volume - 2.0 * lens).abs() < 1e-9, "{volume}");
    }

    #[test]
    fn crossing_cylinders() {
        // A thin cylinder through a thick one, axes crossing at right angles.
        let mut model = Model::new();
        let thick = cylinder(&mut model, (0.0, 0.0, -3.0), Direction::Z, 2.0, 6.0);
        let thin = cylinder(&mut model, (-4.0, 0.0, 0.0), Direction::X, 1.0, 8.0);
        let (volume, _) = check_consistent(&mut model, thick, thin, 1e-6);
        // The intersection's volume: ∫ over the thin disc of the chord across the thick one.
        let n = 20000;
        let exact: f64 = (0..n)
            .map(|i| {
                let y = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
                4.0 * (1.0 - y * y).sqrt() * (4.0 - y * y).sqrt() * 2.0 / n as f64
            })
            .sum();
        assert!((volume - exact).abs() < 1e-5 * exact, "{volume} vs {exact}");
        // Offset so that the thin cylinder only bites into the thick one.
        let biting = cylinder(&mut model, (-4.0, 1.5, 0.0), Direction::X, 1.0, 8.0);
        check_consistent(&mut model, thick, biting, 1e-6);
    }

    fn cone(model: &mut Model, origin: (f64, f64, f64), axis: Direction, r1: f64, r2: f64, height: f64) -> Shape {
        let frame = Frame::from_z(Point::new(origin.0, origin.1, origin.2), axis);
        Cone::new(r1, r2, height).with_frame(frame).build(model).unwrap()
    }

    #[test]
    fn frustum_of_a_cone() {
        // A box cuts a cone (radius 2 down to an apex 3 up) at half its height.
        let mut model = Model::new();
        let spike = cone(&mut model, (0.0, 0.0, 0.0), Direction::Z, 2.0, 0.0, 3.0);
        let lower = cuboid(&mut model, (-3.0, -3.0, -1.0), (6.0, 6.0, 2.5));
        let frustum = intersection(&mut model, spike, lower).unwrap();
        let slant = (1.5f64 * 1.5 + 1.0).sqrt();
        check(&model, frustum, PI * 0.5 * 7.0, PI * 3.0 * slant + 5.0 * PI, 1);
        let tip = difference(&mut model, spike, lower).unwrap();
        check(&model, tip, PI * 1.5 / 3.0, PI * slant + PI, 1);
    }

    #[test]
    fn half_a_cone() {
        // A plane through the axis cuts the cone along two generators.
        let mut model = Model::new();
        let spike = cone(&mut model, (0.0, 0.0, 0.0), Direction::Z, 2.0, 0.0, 3.0);
        let half = cuboid(&mut model, (0.0, -3.0, -1.0), (3.0, 6.0, 5.0));
        let result = difference(&mut model, spike, half).unwrap();
        let slant = 13.0f64.sqrt();
        check(&model, result, 2.0 * PI, PI * slant + 2.0 * PI + 6.0, 1);
    }

    #[test]
    fn cones_against_other_solids() {
        let mut model = Model::new();
        let spike = cone(&mut model, (0.0, 0.0, 0.0), Direction::Z, 2.0, 0.5, 3.0);
        // A plane parallel to a generator meets the cone in a parabola.
        let alpha = (1.5f64 / 3.0).atan();
        let slab = cuboid(&mut model, (-5.0, -5.0, -5.0), (10.0, 10.0, 5.0));
        let tilt = crate::math::Transform::translation(Vector::new(0.5, 0.0, 1.0))
            * crate::math::Transform::rotation(crate::math::Axis::Y, -alpha);
        check_consistent(&mut model, spike, slab.moved(&tilt), 1e-6);
        let pipe = cylinder(&mut model, (-4.0, 0.3, 1.2), Direction::X, 0.6, 8.0);
        check_consistent(&mut model, spike, pipe, 1e-6);
        let ball = sphere(&mut model, (0.4, 0.2, 2.6), 0.9);
        check_consistent(&mut model, spike, ball, 1e-6);
        let other = cone(&mut model, (2.5, 0.0, 1.0), Direction::new(-1.0, 0.0, 0.1), 1.0, 0.3, 4.0);
        check_consistent(&mut model, spike, other, 1e-6);
    }

    #[test]
    fn torus_as_a_surface_of_revolution() {
        let frame = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::from_vector(Vector::new(1.0, 1.0, 2.0)), Direction::X);
        use crate::geometry::ParametricSurface;
        let torus = crate::geometry::surface::Torus::new(frame, 10.0, 3.0);
        let revolution = graph::torus_as_revolution(&torus);
        for (u, v) in [(0.0, 0.0), (0.3, 1.2), (2.5, 4.0), (5.0, 6.0)] {
            assert!(torus.point(u, v).distance(revolution.point(u, v)) < 1e-12, "({u}, {v})");
        }
    }

    #[test]
    fn tubes_turning_a_corner() {
        // Tubes (radius 0.5) along arcs of two circles, one ending where they cross and the
        // other starting there: they touch above and below the corner, where the branches of
        // their section cross. Volumes of the unions as OCCT's.
        use crate::geometry::curve::{Circle, Curve};
        use crate::math::Axis;
        for (r2, expected) in [(9.1, 8.065398), (12.0, 10.09317)] {
            let mut model = Model::new();
            let (c1, c2) = (Point::new(0.0, 10.0, 0.0), Point::new(0.0, -9.0, 0.0));
            let along = (100.0 - r2 * r2 + 19.0 * 19.0) / 38.0;
            let corner = Point::new((100.0f64 - along * along).sqrt(), 10.0 - along, 0.0);
            let disc = |model: &mut Model, centre: Point, up: f64| {
                let normal = Vector::new(0.0, 0.0, up).cross(corner - centre);
                let circle = Circle::new(Frame::new(corner, Direction::from_vector(normal), Direction::X), 0.5);
                let edge = crate::make::edge_from_curve(model, Curve::Circle(circle)).unwrap();
                let wire = crate::make::wire(model, &[edge]).unwrap();
                crate::make::planar_face(model, wire).unwrap()
            };
            let f1 = disc(&mut model, c1, 1.0);
            let start = f1.moved(&crate::math::Transform::rotation(Axis::new(c1, Direction::Z), 0.6));
            let a = crate::sweep::revolve(&mut model, start, Axis::new(c1, -Direction::Z), 0.6).unwrap().shape;
            let f2 = disc(&mut model, c2, -1.0);
            let b = crate::sweep::revolve(&mut model, f2, Axis::new(c2, -Direction::Z), 0.6).unwrap().shape;
            let r = union(&mut model, a, b).unwrap();
            assert!(crate::validity::check(&model, r).is_valid());
            let v = crate::properties::volume(&model, r).mass;
            assert!((v - expected).abs() < 1e-5 * expected, "r2 {r2}: {v} (OCCT {expected})");
        }
    }

    /// A quarter of a small cylinder (a fillet's tool) cut by a torus far larger than it
    /// (R 216, r 60), crossing it from one straight edge to the other: on the torus's grid,
    /// fine for the whole torus, the section was a cell or two wide and lost, and the quarter
    /// kept whole. Volume from OCCT.
    #[test]
    fn a_quarter_cylinder_cut_by_a_torus_far_larger_than_it() {
        use crate::math::{Direction, Frame};
        let mut model = Model::new();
        let block = cuboid(&mut model, (-11.0, 52.5, 24.0), (1.5, 1.5, 8.0));
        let cylinder = crate::primitives::Cylinder::new(1.5, 10.0).with_frame(Frame::new(Point::new(-9.5, 52.5, 23.0), Direction::Z, Direction::X)).build(&mut model).unwrap();
        let quarter = intersection(&mut model, block, cylinder).unwrap();
        let torus = crate::primitives::Torus::new(216.0, 60.0).with_frame(Frame::new(Point::new(0.0, 0.0, 216.0), Direction::Y, Direction::X)).build(&mut model).unwrap();
        let cut = intersection(&mut model, quarter, torus).unwrap();
        assert!(crate::validity::check(&model, cut).is_valid());
        let v = crate::properties::volume(&model, cut).mass;
        assert!((v - 7.291729).abs() < 1e-3, "{v}");
    }

    /// A concave fillet's tool, cut back to the torus its edge runs into, fused onto the
    /// corner of a wall and a rib under that torus (a crop of a beer mat holder). The wall's
    /// edge along the torus (read from a file, straying from it by 2e-6) and the tool's own
    /// along it are one edge; the tool's cylinder, tangent to the rib's side, meets the torus
    /// in a curve tangent to the rib's edge there, met at the corner's vertex. Taken a hair
    /// apart each, they left an edge 2.5e-4 long and faces of no area. Volume and faces from
    /// OCCT.
    #[test]
    fn a_torus_trimmed_fillet_tool_fused_onto_its_corner() {
        let mut model = Model::new();
        let read = |model: &mut Model, text: &str| crate::io::brep::read(model, &mut text.as_bytes()).unwrap();
        let corner = read(&mut model, include_str!("tests/torus_trimmed_fillet_tool_a.brep"));
        let tool = read(&mut model, include_str!("tests/torus_trimmed_fillet_tool_b.brep"));
        let fused = union(&mut model, corner, tool).unwrap();
        assert!(crate::validity::check(&model, fused).is_valid());
        assert_eq!(model.explore(fused, crate::topology::NodeKind::Face).count(), 14);
        for e in model.unique_subshapes(fused, crate::topology::NodeKind::Edge) {
            let length = crate::properties::length(&model, e).mass;
            assert!(length > 1e-3, "an edge {length} long at {:?}", model.point(model.edge_vertices(e).0));
        }
        let v = crate::properties::volume(&model, fused).mass;
        assert!((v - 425.656548).abs() < 1e-3, "{v}");
    }

    #[test]
    fn half_a_torus() {
        // A torus (R 10, r 3) cut by the half-space x > 0: half its volume 2π²Rr² and area
        // 4π²Rr, plus the two disks.
        let mut model = Model::new();
        let torus = crate::primitives::Torus::new(10.0, 3.0).build(&mut model).unwrap();
        check_near(&model, torus, 2.0 * PI * PI * 90.0, 4.0 * PI * PI * 30.0, 1e-9);
        let block = cuboid(&mut model, (0.0, -20.0, -20.0), (20.0, 40.0, 40.0));
        let half = intersection(&mut model, torus, block).unwrap();
        check_near(&model, half, PI * PI * 90.0, 2.0 * PI * PI * 30.0 + 2.0 * PI * 9.0, 1e-6);
    }

    /// A prism whose side is a straight B-spline extruded (a surface of extrusion lying in a
    /// plane), joined to a box whose side is that plane: the extruded curve lies in the plane
    /// along all of it, the surfaces one there, met along no line (sampled, every sample was a
    /// root, hundreds of lines whose crossings took seconds).
    #[test]
    fn prism_on_a_box_along_a_straight_extruded_spline() {
        use crate::geometry::curve::{BSplineCurve, Curve};
        let mut model = Model::new();
        let corners = [(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)].map(|(x, y)| Point::new(x, y, 0.0));
        let straight = |model: &mut Model, a: Point, b: Point| {
            let line = Curve::BSpline(BSplineCurve::new(vec![a, a + (b - a) * 0.5, b], vec![0.0, 0.5, 1.0], vec![2, 1, 2], 1));
            crate::make::edge_from_curve(model, line).unwrap()
        };
        let edges = [
            straight(&mut model, corners[0], corners[1]),
            crate::make::edge_between(&mut model, corners[1], corners[2]).unwrap(),
            straight(&mut model, corners[2], corners[3]),
            crate::make::edge_between(&mut model, corners[3], corners[0]).unwrap(),
        ];
        let wire = crate::make::wire(&mut model, &edges).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let prism = crate::sweep::extrude(&mut model, face, Vector::new(0.0, 0.0, 5.0)).unwrap().shape;
        assert!(model.explore(prism, NodeKind::Face).any(|f| matches!(model.surface(f), Some(crate::geometry::surface::Surface::Extrusion(_)))));
        // (The box beside the prism's side at y = 6, along part of it.)
        let block = cuboid(&mut model, (2.0, 6.0, 1.0), (5.0, 4.0, 3.0));
        let joined = union(&mut model, prism, block).unwrap();
        let report = validity::check(&model, joined);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(&model, joined).mass;
        assert!((v - (300.0 + 60.0)).abs() < 1e-8 * 360.0, "{v}");
    }

    /// A prism whose side is a parabolic arc (a quadratic B-spline, y = 6t(1 - t) over
    /// x = 10t; the arch's area 10) extruded along z, and its prism.
    fn arched_prism(model: &mut Model, z: f64, height: f64) -> Shape {
        use crate::geometry::curve::{BSplineCurve, Curve};
        let at = |x: f64, y: f64| Point::new(x, y, z);
        let arc = Curve::BSpline(BSplineCurve::new(vec![at(0.0, 0.0), at(5.0, 3.0), at(10.0, 0.0)], vec![0.0, 1.0], vec![3, 3], 2));
        let edges = [crate::make::edge_from_curve(model, arc).unwrap(), crate::make::edge_between(model, at(10.0, 0.0), at(0.0, 0.0)).unwrap()];
        let wire = crate::make::wire(model, &edges).unwrap();
        let face = crate::make::planar_face(model, wire).unwrap();
        crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, height)).unwrap().shape
    }

    /// A prism on a closed periodic cubic B-spline (a wavy ring of ten poles), its curve's
    /// seam at +x, raised by `z` and moved across by `(dx, dy)`.
    fn wavy_prism(model: &mut Model, (dx, dy, z): (f64, f64, f64), height: f64) -> Shape {
        use crate::geometry::curve::{BSplineCurve, Curve};
        let poles = (0..10)
            .map(|i| {
                let (angle, r) = (i as f64 * PI / 5.0, if i % 2 == 0 { 5.0 } else { 4.0 });
                Point::new(dx + r * angle.cos(), dy + r * angle.sin(), z)
            })
            .collect();
        let knots: Vec<f64> = (0..=10).map(|k| 3.0 + k as f64).collect();
        let curve = BSplineCurve::try_new_periodic(poles, knots, vec![1; 11], 3).unwrap();
        let edge = crate::make::edge_from_curve(model, Curve::BSpline(curve)).unwrap();
        let wire = crate::make::wire(model, &[edge]).unwrap();
        let face = crate::make::planar_face(model, wire).unwrap();
        crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, height)).unwrap().shape
    }

    /// Two wavy prisms, the second moved across (so its side crosses the first's) and up (so
    /// its bottom cuts the first's side along a closed spline, its seam inside the second): the
    /// stretch of that curve inside the second's bottom runs across the seam. Each result's
    /// shells close; union and intersection add up to the operands, difference and
    /// intersection to the first; and the intersection is the cross-sections' overlap (that of
    /// the prisms moved across only, over their height) over the height they share.
    #[test]
    fn prisms_on_a_periodic_spline_crossing_its_seam() {
        let (h, d) = (5.0, (0.8, 0.4));
        let mut model = Model::new();
        let a = wavy_prism(&mut model, (0.0, 0.0, 0.0), h);
        let b = wavy_prism(&mut model, (d.0, d.1, 1.0), h);
        let across = wavy_prism(&mut model, (d.0, d.1, 0.0), h);
        let volume = |model: &Model, s: Shape| properties::volume(model, s).mass;
        let va = volume(&model, a);
        let mut results = Vec::new();
        for op in [union, intersection, difference] {
            let r = op(&mut model, a, b).unwrap();
            let report = validity::check(&model, r);
            assert!(report.problems.is_empty(), "{report:?}");
            results.push(volume(&model, r));
        }
        let [joined, common, cut] = results[..] else { unreachable!() };
        assert!((joined + common - 2.0 * va).abs() < 1e-8 * va, "{joined} + {common} against {}", 2.0 * va);
        assert!((cut + common - va).abs() < 1e-8 * va, "{cut} + {common} against {va}");
        let overlap = intersection(&mut model, a, across).unwrap();
        let section = volume(&model, overlap) / h;
        assert!((common - section * (h - 1.0)).abs() < 1e-8 * va, "{common} against {}", section * (h - 1.0));
    }

    /// A prism along z, `height` tall, on a conic (an ellipse or circle in a plane square to
    /// z) made a periodic rational B-spline, exactly, of three spans, its origin (where the
    /// prism's side has its seam) moved to knot `origin`, a third of the way round from the
    /// conic's for each.
    fn conic_prism(model: &mut Model, conic: crate::geometry::curve::Curve, height: f64, origin: usize) -> Shape {
        use crate::geometry::curve::Curve;
        let mut spline = crate::geometry::convert::curve_to_bspline(&conic).unwrap();
        assert!(crate::geometry::ParametricCurve::is_periodic(&spline));
        assert_eq!(spline.knots().len(), 4);
        spline.set_origin(origin);
        let edge = crate::make::edge_from_curve(model, Curve::BSpline(spline)).unwrap();
        let wire = crate::make::wire(model, &[edge]).unwrap();
        let face = crate::make::planar_face(model, wire).unwrap();
        crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, height)).unwrap().shape
    }

    /// The ellipse with semi-axes `a` (along x turned by `turn` about z) and `b`, centred on
    /// the z axis in the plane at height `z`.
    fn turned_ellipse(a: f64, b: f64, turn: f64, z: f64) -> crate::geometry::curve::Curve {
        let x = Direction::new(turn.cos(), turn.sin(), 0.0);
        crate::geometry::curve::Curve::Ellipse(crate::geometry::curve::Ellipse::new(Frame::new(Point::new(0.0, 0.0, z), Direction::Z, x), a, b))
    }

    /// The area common to the ellipse of semi-axes `a > b` about the origin and the same
    /// turned by a small `turn`: in polar form the two cross where the angle is `turn / 2`
    /// plus quarter turns, each inside the other on alternate quarters, and the ellipse's
    /// sector from angle 0 to φ has the area `ab/2 · ψ(φ)`, `ψ(φ) = atan2(a sin φ, b cos φ)`.
    fn turned_ellipses_overlap(a: f64, b: f64, turn: f64) -> f64 {
        let psi = |phi: f64| (a * phi.sin()).atan2(b * phi.cos());
        let h = 0.5 * turn;
        a * b * (psi(0.5 * PI + h) - psi(h) + psi(PI - h) - psi(0.5 * PI - h))
    }

    /// Prisms on an ellipse and on the same turned a little about its centre (a thousandth
    /// to a ten-thousandth of a turn), the second raised: their sides cross along four lines
    /// at shallow angles, between which the outlines lie a hair apart (2e-3 to 2e-4 at most).
    /// Every result is valid, and its volume the outlines' overlap (in closed form) times the
    /// height they share, or what is left of an operand: the slivers between the outlines
    /// kept apart, each side cut to its own outline (two outlines between the same points,
    /// a little apart, are not one edge).
    #[test]
    fn prisms_on_an_ellipse_and_the_same_turned_a_little() {
        let (a, b, h, raised) = (5.0, 3.0, 4.0, 1.5);
        for turn in [1e-3, 1e-4] {
            let mut model = Model::new();
            let first = conic_prism(&mut model, turned_ellipse(a, b, 0.0, 0.0), h, 1);
            let second = conic_prism(&mut model, turned_ellipse(a, b, turn, raised), h, 1);
            let whole = PI * a * b * h;
            let common = turned_ellipses_overlap(a, b, turn) * (h - raised);
            let expected = [(intersection as fn(&mut Model, Shape, Shape) -> Result<Shape>, first, second, common), (intersection, second, first, common), (union, first, second, 2.0 * whole - common), (difference, first, second, whole - common), (difference, second, first, whole - common)];
            for (k, (op, x, y, want)) in expected.into_iter().enumerate() {
                let r = op(&mut model, x, y).unwrap();
                let report = validity::check(&model, r);
                assert!(report.problems.is_empty(), "turned {turn}, result {k}: {report:?}");
                let v = properties::volume(&model, r).mass;
                assert!((v - want).abs() < 1e-10 * whole, "turned {turn}, result {k}: {v} against {want}");
            }
        }
    }

    /// Prisms on a circle and on the same moved across by 3e-6 of its radius 3: their sides
    /// cross at an angle of 1e-6, the outlines within the tolerance of each other for a
    /// stretch either side of each crossing (crossed there, not run along), and 3e-6 apart at
    /// most between (the sides' seams where they are farthest apart). The intersection is the
    /// lens the circles share (in closed form) times the height they share, and the results
    /// add up (to the part in 1e-9 the faces' integrals on these rational outlines reach).
    #[test]
    fn prisms_on_a_circle_and_the_same_moved_a_hair() {
        let (r, d, h, raised) = (3.0, 3e-6, 4.0, 1.5);
        let mut model = Model::new();
        let circle = |x: f64, z: f64| crate::geometry::curve::Curve::Circle(crate::geometry::curve::Circle::new(Frame::new(Point::new(x, 0.0, z), Direction::Z, Direction::X), r));
        let first = conic_prism(&mut model, circle(0.0, 0.0), h, 0);
        let second = conic_prism(&mut model, circle(d, raised), h, 0);
        let whole = properties::volume(&model, first).mass;
        assert!((whole - PI * r * r * h).abs() < 1e-11 * whole, "{whole}");
        let lens = 2.0 * r * r * (0.5 * d / r).acos() - 0.5 * d * (4.0 * r * r - d * d).sqrt();
        let common = lens * (h - raised);
        for (op, x, y, want) in [(intersection as fn(&mut Model, Shape, Shape) -> Result<Shape>, first, second, common), (intersection, second, first, common), (union, first, second, 2.0 * whole - common), (difference, first, second, whole - common)] {
            let r = op(&mut model, x, y).unwrap();
            let report = validity::check(&model, r);
            assert!(report.problems.is_empty(), "{report:?}");
            let v = properties::volume(&model, r).mass;
            assert!((v - want).abs() < 1e-9 * whole, "{v} against {want}");
        }
    }

    /// Prisms on a closed periodic B-spline outline (eight poles about a circle, its knots from
    /// 3) and on the same moved across by 2.9e-6, the second lowered by a fifteenth of their
    /// height: the outlines cross twice at angles of about 1e-6, within the classifying
    /// distance of each other for a long stretch about each crossing. The bottom
    /// outline of each meets the other's side where the outlines cross, found as the
    /// sides' lines of meeting are (in closed form, the classifying distance not taken for
    /// meeting along a stretch). Every result is valid and they add up; the intersection
    /// lies between the prism's cross-section and that less a strip as wide as the move
    /// along the outline, over the height they share.
    #[test]
    fn prisms_on_an_outline_and_the_same_moved_a_hair() {
        use crate::geometry::curve::{BSplineCurve, Curve};
        let poles = [
            (2.4775224573227153, 0.0),
            (1.9866654191666364, 1.9866654191666362),
            (1.776982384600636e-16, 2.9020324649324856),
            (-2.4924892437560295, 2.49248924375603),
            (-3.201553642958206, 3.920772421147315e-16),
            (-1.7113199793120286, -1.7113199793120282),
            (-6.092011759060008e-16, -3.3163367393229923),
            (1.625845109015556, -1.6258451090155568),
        ];
        let (h, lowered, (dx, dy)) = (8.366476807755433, 0.5641892695822941, (2.8569289529656317e-06, -6.410179343273588e-07));
        let prism = |model: &mut Model, (dx, dy, z): (f64, f64, f64)| {
            let curve = BSplineCurve::try_new_periodic(poles.iter().map(|&(x, y)| Point::new(x + dx, y + dy, z)).collect(), (3..=11).map(f64::from).collect(), vec![1; 9], 3).unwrap();
            let edge = crate::make::edge_from_curve(model, Curve::BSpline(curve)).unwrap();
            let wire = crate::make::wire(model, &[edge]).unwrap();
            let face = crate::make::planar_face(model, wire).unwrap();
            crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, h)).unwrap().shape
        };
        let mut model = Model::new();
        let a = prism(&mut model, (0.0, 0.0, 0.0));
        let b = prism(&mut model, (dx, dy, -lowered));
        let whole = properties::volume(&model, a).mass;
        let section = whole / h;
        let outline = model.explore(a, NodeKind::Edge).filter(|&e| matches!(model.edge_curve(e), Some((Curve::BSpline(_), _, _)))).map(|e| properties::length(&model, e).mass).next().unwrap();
        let mut volumes = Vec::new();
        for (op, x, y) in [(intersection as fn(&mut Model, Shape, Shape) -> Result<Shape>, a, b), (intersection, b, a), (union, a, b), (difference, a, b), (difference, b, a)] {
            let r = op(&mut model, x, y).unwrap();
            let report = validity::check(&model, r);
            assert!(report.problems.is_empty(), "{report:?}");
            volumes.push(properties::volume(&model, r).mass);
        }
        let [common, other_way, joined, a_less, b_less] = volumes[..] else { unreachable!() };
        let shared = h - lowered;
        let strip = (dx * dx + dy * dy).sqrt() * outline * shared;
        assert!(common <= section * shared + 1e-9 * whole && common >= section * shared - strip, "{common} against {}", section * shared);
        for (sum, want) in [(other_way, common), (joined + common, 2.0 * whole), (a_less + common, whole), (b_less + common, whole)] {
            assert!((sum - want).abs() < 1e-9 * whole, "{sum} against {want}");
        }
    }

    /// Prisms on an ellipse and on the same turned by 3e-7 about its centre: the outlines
    /// cross four times, at most 6e-7 apart in between, within the classifying distance all
    /// round: one outline as far as the operations can tell, the sides lying on one another.
    /// Every result is valid, the union the prism made taller, and the volumes are those of
    /// the exact solids within the slivers' (the outlines' areas apart, times the height).
    #[test]
    fn prisms_on_an_ellipse_and_the_same_turned_a_hair() {
        let (a, b, h, raised, turn) = (5.0, 3.0, 4.0, 1.5, 3e-7);
        let mut model = Model::new();
        let first = conic_prism(&mut model, turned_ellipse(a, b, 0.0, 0.0), h, 1);
        let second = conic_prism(&mut model, turned_ellipse(a, b, turn, raised), h, 1);
        let whole = PI * a * b * h;
        let common = turned_ellipses_overlap(a, b, turn) * (h - raised);
        let slivers = (PI * a * b - turned_ellipses_overlap(a, b, turn)) * h;
        assert!(slivers > 0.0 && slivers < 1e-4, "{slivers}");
        for (op, x, y, want) in [(intersection as fn(&mut Model, Shape, Shape) -> Result<Shape>, first, second, common), (intersection, second, first, common), (union, first, second, 2.0 * whole - common), (difference, first, second, whole - common)] {
            let r = op(&mut model, x, y).unwrap();
            let report = validity::check(&model, r);
            assert!(report.problems.is_empty(), "{report:?}");
            let v = properties::volume(&model, r).mass;
            assert!((v - want).abs() < slivers + 1e-10 * whole, "{v} against {want}");
        }
    }

    /// Two arched prisms, one on the other's side surface along all of it (moved along the
    /// direction of extrusion): their curved sides are one surface, met along no line, and
    /// the union is the prism made taller.
    #[test]
    fn prisms_on_one_extruded_surface() {
        let mut model = Model::new();
        let (lower, upper) = (arched_prism(&mut model, 0.0, 5.0), arched_prism(&mut model, 2.0, 5.0));
        let joined = union(&mut model, lower, upper).unwrap();
        let report = validity::check(&model, joined);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(&model, joined).mass;
        assert!((v - 10.0 * 7.0).abs() < 1e-8 * 70.0, "{v}");
    }

    /// An arched prism and a cylinder along the same direction, crossing its arc: they meet
    /// along two lines, and the union's volume is the cross-sections' union (the overlap
    /// integrated across x) times the height.
    #[test]
    fn an_extruded_arc_across_a_parallel_cylinder() {
        let mut model = Model::new();
        let prism = arched_prism(&mut model, 0.0, 5.0);
        let (cx, cy, r) = (5.0, 1.5, 2.0);
        let cylinder = cylinder(&mut model, (cx, cy, -1.0), Direction::Z, r, 7.0);
        let joined = union(&mut model, prism, cylinder).unwrap();
        let report = validity::check(&model, joined);
        assert!(report.problems.is_empty(), "{report:?}");
        // (Vertical lines through points on both the circle and the arch: the two where
        // they cross, each perhaps in pieces.)
        let top = |x: f64| 6.0 * (x / 10.0) * (1.0 - x / 10.0);
        let mut crossings: Vec<(f64, f64)> = Vec::new();
        for e in model.explore(joined, NodeKind::Edge) {
            let Some((crate::geometry::curve::Curve::Line(_), _, _)) = model.edge_curve(e) else { continue };
            let (a, b) = model.edge_vertices(e);
            let (p, q) = (model.point(a), model.point(b));
            let on_circle = ((p.x - cx).hypot(p.y - cy) - r).abs() < 1e-6;
            if (p - q).cross(Vector::new(0.0, 0.0, 1.0)).length() < 1e-9 && on_circle && (p.y - top(p.x)).abs() < 1e-6 && !crossings.iter().any(|c| (c.0 - p.x).hypot(c.1 - p.y) < 1e-6) {
                crossings.push((p.x, p.y));
            }
        }
        assert_eq!(crossings.len(), 2, "{crossings:?}");
        // (The overlap of the arch and the disc, by Simpson's rule in x.)
        let overlap_at = |x: f64| {
            let h = (r * r - (x - cx) * (x - cx)).max(0.0).sqrt();
            let top = 6.0 * (x / 10.0) * (1.0 - x / 10.0);
            (top.min(cy + h) - (cy - h).max(0.0)).max(0.0)
        };
        let n = 200_000;
        let (a, b) = (cx - r, cx + r);
        let step = (b - a) / n as f64;
        let overlap = (0..=n).map(|i| overlap_at(a + step * i as f64) * if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 }).sum::<f64>() * step / 3.0;
        let v = properties::volume(&model, joined).mass;
        let expected = 5.0 * (10.0 - overlap) + PI * r * r * 7.0;
        assert!((v - expected).abs() < 1e-6 * expected, "{v} against {expected}");
    }

    /// A unit cube cut by a slab across the middle of its bottom (CadQuery's history test):
    /// recorded, the cube's bottom face makes the two strips the cut leaves of it, an edge
    /// across the slab's two pieces, and the slab's bottom face (on the cube's, taken away)
    /// nothing; an untouched face is kept as itself.
    #[test]
    fn a_cut_records_where_its_faces_come_from() {
        let mut model = Model::new();
        let cube = cuboid(&mut model, (-0.5, -0.5, 0.0), (1.0, 1.0, 1.0));
        let slab = cuboid(&mut model, (-0.5, -0.25, 0.0), (1.0, 0.5, 0.1));
        let at = |model: &Model, shape: Shape, z: f64| model.explore(shape, NodeKind::Face).find(|&f| (properties::area(model, f).center_of_mass.z - z).abs() < 1e-9).unwrap();
        let (bottom, top, slab_bottom) = (at(&model, cube, 0.0), at(&model, cube, 1.0), at(&model, slab, 0.0));
        let side_edge = model.explore(bottom, NodeKind::Edge).find(|&e| { let (a, b) = model.edge_vertices(e); model.point(a).x > 0.49 && model.point(b).x > 0.49 }).unwrap();
        record_begin(&model);
        let result = difference(&mut model, cube, slab).unwrap();
        let history = record_end(&model, result);
        let images = |s: Shape| history.iter().filter(|(from, _)| *from == s.node).count();
        assert_eq!(images(bottom), 2);
        assert_eq!(images(side_edge), 2);
        assert_eq!(images(slab_bottom), 0);
        assert!(history.contains(&(top.node, top.node)) || images(top) == 1);
    }

    /// Faces a cut keeps whole are recorded too, those of surfaces it takes as simpler ones
    /// (a sweep's planar B-spline sides, seen as planes) as well.
    #[test]
    fn a_cut_records_the_faces_it_keeps_whole() {
        let mut model = Model::new();
        let spine = crate::make::edge_between(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 1.0)).unwrap();
        let p = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)].map(|(x, y)| Point::new(x, y, 0.0));
        let edges: Vec<Shape> = (0..4).map(|i| crate::make::edge_between(&mut model, p[i], p[(i + 1) % 4]).unwrap()).collect();
        let square = crate::make::wire(&mut model, &edges).unwrap();
        let prism = crate::sweep::pipe_through(&mut model, &[(square, Default::default())], spine, true).unwrap();
        let hole = cuboid(&mut model, (-0.1, -0.1, -1.0), (0.2, 0.2, 3.0));
        let sides: Vec<Shape> = model.explore(prism, NodeKind::Face).filter(|&f| properties::area(&model, f).center_of_mass.z.abs() > 0.4 && properties::area(&model, f).center_of_mass.z < 0.6).collect();
        assert_eq!(sides.len(), 4);
        record_begin(&model);
        let result = difference(&mut model, prism, hole).unwrap();
        let history = record_end(&model, result);
        for side in sides {
            assert_eq!(history.iter().filter(|(from, _)| *from == side.node).count(), 1);
        }
    }

    #[test]
    fn extruded_ellipse() {
        // An elliptic prism (semi-axes 10 and 5, height 30): volume π a b h, and half of it
        // where a box covers x > 0.
        let mut model = Model::new();
        let ellipse = crate::geometry::curve::Ellipse::new(Frame::WORLD, 10.0, 5.0);
        let edge = crate::make::edge_from_curve(&mut model, crate::geometry::curve::Curve::Ellipse(ellipse)).unwrap();
        let wire = crate::make::wire(&mut model, &[edge]).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let prism = crate::sweep::extrude(&mut model, face, Vector::new(0.0, 0.0, 30.0)).unwrap().shape;
        let volume = PI * 10.0 * 5.0 * 30.0;
        assert!((properties::volume(&model, prism).mass - volume).abs() < 1e-8 * volume);
        let block = cuboid(&mut model, (0.0, -20.0, -10.0), (20.0, 40.0, 50.0));
        let half = intersection(&mut model, prism, block).unwrap();
        let report = validity::check(&model, half);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(&model, half).mass;
        assert!((v - volume / 2.0).abs() < 1e-6 * volume, "{v} != {}", volume / 2.0);
    }

    #[test]
    fn surface_of_revolution() {
        // A stadium in a plane parallel to the axis, revolved: its arcs sweep surfaces of
        // revolution that are neither tori nor quadrics (bcut_simple/G6).
        let mut model = Model::new();
        let block = cuboid(&mut model, (0.0, 0.0, 0.0), (100.0, 100.0, 40.0));
        let z = 40.0;
        let p = |x: f64, y: f64| Point::new(x, y, z);
        let arc = |model: &mut Model, center: Point, start: Point, end: Point| {
            let circle = crate::geometry::curve::Circle::new(Frame::new(center, Direction::Z, Direction::from_vector(start - center)), 10.0);
            let angle = {
                let (a, b) = (start - center, end - center);
                let t = a.cross(b).z.atan2(a.dot(b));
                if t <= 0.0 { t + 2.0 * PI } else { t }
            };
            crate::make::edge_from_curve_range(model, crate::geometry::curve::Curve::Circle(circle), 0.0, angle).unwrap()
        };
        let e1 = crate::make::edge_between(&mut model, p(50.0, 20.0), p(50.0, 70.0)).unwrap();
        let e2 = arc(&mut model, p(40.0, 70.0), p(50.0, 70.0), p(30.0, 70.0));
        let e3 = crate::make::edge_between(&mut model, p(30.0, 70.0), p(30.0, 20.0)).unwrap();
        let e4 = arc(&mut model, p(40.0, 20.0), p(30.0, 20.0), p(50.0, 20.0));
        let wire = crate::make::wire(&mut model, &[e1, e2, e3, e4]).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let ring = crate::sweep::revolve(&mut model, face, crate::math::Axis::new(Point::new(0.0, 0.0, 50.0), Direction::Y), 2.0 * PI)
            .unwrap()
            .shape;
        // The ring itself: area 2π ∮ ρ ds_m (ds_m the arc length in the meridian half-plane,
        // shorter than the profile's own for these arcs off the axis's planes) and volume
        // 2π x̄ A (each point of the stadium turns at a distance whose projection is x).
        let rho = |x: f64| (x * x + 100.0f64).sqrt();
        let n = 20000;
        let mut area = 50.0 * (rho(30.0) + rho(50.0));
        for i in 0..n {
            let a = PI * (i as f64 + 0.5) / n as f64;
            let x = 40.0 + 10.0 * a.cos();
            let (dx, dy) = (-10.0 * a.sin(), 10.0 * a.cos());
            area += 2.0 * rho(x) * ((x * dx / rho(x)).powi(2) + dy * dy).sqrt() * PI / n as f64;
        }
        let area = 2.0 * PI * area;
        let volume = 2.0 * PI * 40.0 * (20.0 * 50.0 + PI * 100.0);
        let a = properties::area(&model, ring).mass;
        let v = properties::volume(&model, ring).mass;
        assert!((a - area).abs() < 1e-6 * area, "ring area {a} != {area}");
        assert!((v - volume).abs() < 1e-6 * volume, "ring volume {v} != {volume}");
        // OCCT's bcut_simple/G6 expects 41187.4 for the difference's area (checked independently by
        // integrating over the box's faces outside the ring and the ring's inside the box).
        let (_, common_area) = check_consistent(&mut model, block, ring, 1e-6);
        let _ = common_area;
        let r = difference(&mut model, block, ring).unwrap();
        let a = properties::area(&model, r).mass;
        assert!((a - 41187.4).abs() < 1e-5 * 41187.4, "difference area {a}");
    }

    #[test]
    fn cones_with_parallel_axes() {
        // bopcommon_simple/ZN2: OCCT expects the intersection to have area 340.979.
        let mut model = Model::new();
        let a = cone(&mut model, (0.0, 0.0, 0.0), Direction::Z, 8.0, 4.0, 8.0);
        let b = cone(&mut model, (4.0, 0.0, 0.0), Direction::Z, 8.0, 4.0, 4.0);
        let (_, area) = check_consistent(&mut model, a, b, 1e-6);
        assert!((area - 340.979).abs() < 1e-5 * 340.979, "intersection area {area}");
    }

    /// A cylinder (radius 54, 11 long, along X from x = 0 at z = 6) and a cone about a crossing
    /// axis far off (along Y at z = 216, 45°, apex at y = -216): their section's branches, by
    /// angle round the cylinder, meet where its generator only grazes the cone, at the
    /// cylinder's end and with the linear term of their quadratic zero there by symmetry.
    /// Joined there, at the point the end's edge meets the cone, the section bounds the faces
    /// (OCCT's common: 15166.010).
    #[test]
    fn a_cylinder_and_a_cone_about_crossing_axes_far_apart() {
        let mut model = Model::new();
        let a = cylinder(&mut model, (0.0, 0.0, 6.0), Direction::X, 54.0, 11.0);
        let b = cone(&mut model, (0.0, -216.0, 216.0), Direction::Y, 0.0, 216.0, 216.0);
        let both = intersection(&mut model, a, b).unwrap();
        let report = validity::check(&model, both);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, both).mass;
        assert!((volume - 15166.010).abs() < 1e-5 * 15166.010, "{volume}");
    }

    /// A torus cut by a box turned about its axis, one face through the axis: that face meets
    /// the torus in two of its meridian circles, exactly.
    #[test]
    fn a_plane_through_a_torus_axis_meets_it_in_meridian_circles() {
        use crate::geometry::curve::Curve;
        use crate::primitives::Torus;
        let mut model = Model::new();
        let torus = Torus::new(10.0, 2.0).build(&mut model).unwrap();
        let turn = 0.3f64;
        // (Its face at local y = 0 through the axis, from 15 one way along it to 15 the other.)
        let frame = Frame::new(Point::new(-15.0 * turn.cos(), -15.0 * turn.sin(), -5.0), Direction::Z, Direction::new(turn.cos(), turn.sin(), 0.0));
        let slab = Cuboid::new(30.0, 20.0, 10.0).with_frame(frame).build(&mut model).unwrap();
        let curves = section(&mut model, torus, slab).unwrap();
        let circles: Vec<f64> = model.unique_subshapes(curves, NodeKind::Edge).into_iter().filter_map(|e| match model.edge_curve(e) {
            Some((Curve::Circle(c), _, _)) if (c.radius() - 2.0).abs() < 1e-9 => Some(c.radius()),
            _ => None,
        }).collect();
        assert!(circles.len() >= 2, "{circles:?}");
    }

    /// A loft between two ellipses (closed in u without a period) cut by a box's top: the
    /// cut, a closed section starting on the loft's seam, runs against u there; its pcurve is
    /// fitted from the arc's own start on the seam, not a polyline whose chord the wrong way
    /// round matched the curve at its middle alone (OCCT's volume, 182.3066).
    #[test]
    fn a_loft_cut_across_its_seam() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let loft = read(include_str!("tests/ellipse_loft_a.brep"));
        let slab = read(include_str!("tests/ellipse_loft_b.brep"));
        let result = intersection(&mut model, loft, slab).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 182.3066).abs() < 1e-3, "{volume}");
    }

    /// A square bar lofted round at its end (to a circle through the square's middles) and
    /// cut down to a cylinder through that circle: the loft's sides meet the cylinder along
    /// its top circle and in curves from the square's crossings up to it, where the surfaces
    /// touch. The curves end on the circle, the corners outside the cylinder cut away (as
    /// OCCT: 15 faces, 1387.617 by its loft).
    #[test]
    fn a_square_lofted_round_cut_to_a_cylinder() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let bar = read(include_str!("tests/square_lofted_round_a.brep"));
        let cylinder = read(include_str!("tests/square_lofted_round_b.brep"));
        let result = intersection(&mut model, bar, cylinder).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        assert_eq!(model.unique_subshapes(result, crate::topology::NodeKind::Face).len(), 15);
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 1387.559).abs() < 1e-2, "{volume}");
    }

    /// A block with chamfered curved edges fused to a frame whose free-form face comes down
    /// onto the block's flat side, touching it: the section along that touch and the one
    /// across the chamfer meet the side's corner at a shallow angle, their pcurves crossing
    /// a few millionths from the vertex they share. The vertex covers it: a valid solid.
    #[test]
    fn a_chamfered_block_fused_where_a_curved_face_touches_it() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let block = read(include_str!("tests/chamfered_block_on_curved_frame_a.brep"));
        let frame = read(include_str!("tests/chamfered_block_on_curved_frame_b.brep"));
        let result = union(&mut model, block, frame).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 31.72242).abs() < 1e-4, "{volume}");
    }

    /// Two lattice nodes fused, each cut flat at the same height into a small square cap,
    /// one cap turned 45° to the other, and struts whose cylinders touch the caps' plane
    /// along the squares' sides: the corners of each cap beyond the other lie near those
    /// cylinders within the tolerance, but a plane only touches a cylinder, so they are not
    /// on it. The caps' union is kept whole (an octagon and eight corners, as OCCT).
    #[test]
    fn lattice_node_caps_turned_to_one_another_fused() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let a = read(include_str!("tests/lattice_node_caps_a.brep"));
        let b = read(include_str!("tests/lattice_node_caps_b.brep"));
        let result = union(&mut model, a, b).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 0.181326).abs() < 1e-5, "{volume}");
    }

    /// A few degrees of a revolved ring cut by a thick disc square to it: the disc's flat end
    /// meets the ring's torus in a closed curve, traced on a grid; through its first few dozen
    /// samples it strayed from the torus by thousandths, the disc's end face its area off by
    /// a tenth. Refined until it follows both surfaces (OCCT's volume).
    #[test]
    fn a_ring_segment_cut_by_a_disc() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let ring = read(include_str!("tests/ring_segment_cut_by_a_disc_a.brep"));
        let disc = read(include_str!("tests/ring_segment_cut_by_a_disc_b.brep"));
        let result = difference(&mut model, ring, disc).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 56644.4472).abs() < 0.02, "{volume}");
    }

    /// A cylinder with a window through its side, the window's pcurves kept a period off the
    /// rest of the face (as a model may keep them): its window is still a hole in the face,
    /// not a face of its own, when the solid meets another (the same volume as without the
    /// shift).
    #[test]
    fn a_hole_kept_a_period_off_its_face() {
        use crate::geometry::curve2::Curve2;
        use crate::math::{Transform2, Transformable2, Vector2};
        use crate::topology::CoEdgeId;
        let mut model = Model::new();
        let tube = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 10.0, 20.0);
        let window = cuboid(&mut model, (-12.0, -2.0, 5.0), (4.0, 4.0, 6.0));
        let windowed = difference(&mut model, tube, window).unwrap();
        let cut = cuboid(&mut model, (-20.0, -20.0, -1.0), (40.0, 40.0, 16.0));
        let expected = {
            let copy = crate::modify::copy(&mut model, windowed, false).unwrap();
            let r = intersection(&mut model, copy, cut).unwrap();
            properties::volume(&model, r).mass
        };
        // The side face's second wire (the window), a period on.
        let side = model.unique_subshapes(windowed, NodeKind::Face).into_iter().find(|&f| matches!(model.surface(f), Some(crate::geometry::surface::Surface::Cylinder(_))) && model.children(f).count() == 2).unwrap();
        let hole = model.children(side).nth(1).unwrap();
        let coedges: Vec<CoEdgeId> = model.children(hole).map(|c| CoEdgeId::try_from(c.node).unwrap()).collect();
        for c in coedges {
            let (p, _, _) = model.effective_pcurve(c).unwrap();
            let shifted: Curve2 = p.transformed(&Transform2::translation(Vector2::new(2.0 * PI, 0.0)));
            model.set_pcurve(c, shifted).unwrap();
        }
        let r = intersection(&mut model, windowed, cut).unwrap();
        let report = validity::check(&model, r);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, r).mass;
        assert!((volume - expected).abs() < 1e-6 * expected, "{volume} {expected}");
    }

    #[test]
    fn rotated_box() {
        // A cube turned 45° about Z, centred on a face of another: the difference removes a prism
        // with a square-diagonal cross-section.
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
        let h = std::f64::consts::FRAC_1_SQRT_2;
        let frame = Frame::new(Point::new(1.0, 1.0 - h, -1.0), Direction::Z, Direction::new(h, h, 0.0));
        let b = Cuboid::new(1.0, 1.0, 3.0).with_frame(frame).build(&mut model).unwrap();
        let result = difference(&mut model, a, b).unwrap();
        check(&model, result, 4.0 - 1.0, 2.0 * 3.0 + 8.0 + 4.0 * 1.0, 1);
    }

    /// A disc with a groove 0.002 deep round its side (radius 16.448 between z 1.05 and 1.95,
    /// 16.45 below and above), a box cut into its rim across the groove's upper step: that
    /// step, a ring 0.002 wide, is split into a C and the box's square, the C (thinner than
    /// its polylines' sampling, its point taken from them off the ring's own) kept.
    #[test]
    fn a_ring_a_thousandth_of_its_radius_wide_cut_across_keeps_its_rest() {
        let mut model = Model::new();
        let (r, groove) = (16.45, 16.448);
        let below = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, r, 1.05);
        let middle = cylinder(&mut model, (0.0, 0.0, 1.05), Direction::Z, groove, 0.9);
        let above = cylinder(&mut model, (0.0, 0.0, 1.95), Direction::Z, r, 1.55);
        let disc = union(&mut model, below, middle).unwrap();
        let disc = union(&mut model, disc, above).unwrap();
        let block = cuboid(&mut model, (-0.7, 15.75, 1.2), (1.4, 1.4, 3.5));
        let result = difference(&mut model, disc, block).unwrap();
        // The block's part inside: |x| <= 0.7, y >= 15.75, within the circle at each height.
        let section = |r: f64| {
            let f = |x: f64| 0.5 * x * (r * r - x * x).sqrt() + 0.5 * r * r * (x / r).asin();
            f(0.7) - f(-0.7) - 1.4 * 15.75
        };
        let volume = PI * (r * r * 2.6 + groove * groove * 0.9) - 0.75 * section(groove) - 1.55 * section(r);
        let report = validity::check(&model, result);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(&model, result).mass;
        assert!((v - volume).abs() < 1e-6 * volume, "{v} != {volume}");
        assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), 12);
    }

    /// A prism on a regular polygon of 32 and of 128 sides (corners at radius 62.5, 1 high),
    /// a cylinder 1 high standing on it, its radius halfway into the sides' sag: the two
    /// faces on the plane between them cross at two points by each side. The union's faces
    /// there are the polygon's corners past the circle and the disc's rims past the sides,
    /// each a few thousandths across, thinner than the polylines' sag (the disc far from the
    /// polygon's side faces, its circle's chords cutting inside the corners' points; its
    /// rims' points off the disc's own polylines). OCCT's volumes and faces.
    #[test]
    fn a_cylinder_standing_on_a_polygon_just_past_its_sides() {
        for (n, faces) in [(32usize, 99usize), (128, 387)] {
            let mut model = Model::new();
            let big = 62.5f64;
            let r = big * (1.0 - 0.5 * (1.0 - (PI / n as f64).cos()));
            let corners: Vec<Point> = (0..n).map(|i| {
                let a = 2.0 * PI * i as f64 / n as f64;
                Point::new(big * a.cos(), big * a.sin(), 0.0)
            }).collect();
            let edges: Vec<Shape> = (0..n).map(|i| crate::make::edge_between(&mut model, corners[i], corners[(i + 1) % n]).unwrap()).collect();
            let wire = crate::make::wire(&mut model, &edges).unwrap();
            let face = crate::make::planar_face(&mut model, wire).unwrap();
            let polygon = crate::sweep::extrude(&mut model, face, Vector::new(0.0, 0.0, 1.0)).unwrap().shape;
            let disc = cylinder(&mut model, (0.0, 0.0, 1.0), Direction::Z, r, 1.0);
            let result = union(&mut model, polygon, disc).unwrap();
            let report = validity::check(&model, result);
            assert!(report.problems.is_empty(), "{n}: {report:?}");
            let volume = 0.5 * n as f64 * big * big * (2.0 * PI / n as f64).sin() + PI * r * r;
            let v = properties::volume(&model, result).mass;
            assert!((v - volume).abs() < 1e-6 * volume, "{n}: {v} != {volume}");
            assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), faces, "{n}");
        }
    }

    /// A disc 2 thick with a recess of radius 62.495 from below to its middle, a hole of
    /// radius 62.5 from above to the same plane (its seam turned 1 degree): left between them
    /// on that plane, a ring 0.005 wide, narrower than the circles' polylines' sag. The hole's
    /// bottom is split by the recess's circle into the ring and a disc on the recess's
    /// ceiling; the ring's inner circle, its points off the bottom's polylines, is still its
    /// hole, and its point is taken where it is that thin. OCCT's volume and faces.
    #[test]
    fn a_ring_left_between_a_hole_and_a_recess_a_hair_narrower() {
        let mut model = Model::new();
        let disc = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 70.0, 2.0);
        let recess = cylinder(&mut model, (0.0, 0.0, 0.0), Direction::Z, 62.495, 1.0);
        let recessed = difference(&mut model, disc, recess).unwrap();
        let turned = 1f64.to_radians();
        let frame = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::new(turned.cos(), turned.sin(), 0.0));
        let hole = Cylinder::new(62.5, 1.0).with_frame(frame).build(&mut model).unwrap();
        let result = difference(&mut model, recessed, hole).unwrap();
        let report = validity::check(&model, result);
        assert!(report.problems.is_empty(), "{report:?}");
        let volume = PI * (70.0f64.powi(2) * 2.0 - 62.495f64.powi(2) - 62.5f64.powi(2));
        let v = properties::volume(&model, result).mass;
        assert!((v - volume).abs() < 1e-6 * volume, "{v} != {volume}");
        assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), 6);
    }

    /// The section of `a` and `b` read from the fixtures: its edges' lengths and ends.
    fn fixture_section(a: &str, b: &str) -> Vec<(f64, Point, Point)> {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let (a, b) = (read(a), read(b));
        let s = section(&mut model, a, b).unwrap();
        model
            .unique_subshapes(s, NodeKind::Edge)
            .into_iter()
            .map(|e| {
                let (v0, v1) = model.edge_vertices(e);
                (properties::length(&model, e).mass, model.point(v0), model.point(v1))
            })
            .collect()
    }

    /// The end walls of two thickened plates mirrored onto one another (cells of a Schwarz
    /// lattice), each from the arc their middles share outwards, twisting across the mirror
    /// plane: their section is the arc and two branches running off it where the walls cross
    /// the plane, the surfaces tangent there. The trace of the second branch stops just past
    /// the arc, its points bunching as it ends, and was met with the arc only within the
    /// last tenth of its parameters, too short a way: the arc is split at both branches, as
    /// OCCT's.
    #[test]
    fn mirrored_walls_crossing_split_their_shared_arc_at_both_branches() {
        let edges = fixture_section(include_str!("tests/mirrored_walls_crossing_a.brep"), include_str!("tests/mirrored_walls_crossing_b.brep"));
        for junction in [Point::new(3.98, 0.2580, 1.9732), Point::new(3.98, 1.9669, 0.3022)] {
            let meeting = edges.iter().filter(|(_, p, q)| p.distance(junction) < 1e-3 || q.distance(junction) < 1e-3).count();
            assert_eq!(meeting, 3, "{junction:?}: {edges:?}");
        }
        let total: f64 = edges.iter().map(|e| e.0).sum();
        assert!((total - 4.6262).abs() < 1e-3, "{total}");
    }

    /// The side walls of the same plates, both lying nearly flat and continuing one another
    /// across the mirror plane, from a corner they share: their section runs from that corner
    /// along the plane (0.7502 long, as OCCT's). Met by the corner's vertex only within the
    /// tolerance it was dropped, as a vertex of the edges already.
    #[test]
    fn mirrored_walls_continuing_meet_along_the_mirror_plane() {
        let edges = fixture_section(include_str!("tests/mirrored_walls_continuing_a.brep"), include_str!("tests/mirrored_walls_continuing_b.brep"));
        let corner = Point::new(3.98, 1.99, 0.0);
        assert!(
            edges.iter().any(|(l, p, q)| (l - 0.7502).abs() < 1e-3 && (p.distance(corner) < 1e-6 || q.distance(corner) < 1e-6)),
            "{edges:?}"
        );
    }

    /// An offset face of a thickened plate (a cell of a Schwarz lattice) and its mirror image
    /// across x = 1.99, meeting end to end there near tangent: each one's edge meets the
    /// other's face a little apart (4e-5), and branches of their section from either side stop
    /// short of those points. Carried on a whole step to them (as two fillets meeting at a
    /// corner), each to a different one, they made a small loop of two edges between the two
    /// points, and the union of the two plates was left open.
    #[test]
    fn mirrored_offset_faces_leave_no_loops_where_their_edges_meet() {
        use crate::math::{Direction, Transform};
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut include_str!("tests/mirrored_offset_faces_a.brep").as_bytes()).unwrap();
        let mirror = Transform::mirror_plane(Point::new(1.99, 0.0, 0.0), Direction::X);
        let b = crate::modify::transform(&mut model, a, &mirror, Default::default()).unwrap();
        let s = section(&mut model, a, b).unwrap();
        for e in model.unique_subshapes(s, NodeKind::Edge) {
            let (v0, v1) = model.edge_vertices(e);
            let length = properties::length(&model, e).mass;
            assert!(model.point(v0).distance(model.point(v1)) > 1e-4 || length < 1e-4, "a loop of {length} at {:?}", model.point(v0));
        }
    }

    /// Panics where two edges of the section of `a` and `b` (BRep text) cover one stretch:
    /// their ends a hair apart, and their middles.
    fn section_covering_nothing_twice(a: &str, b: &str) {
        use crate::geometry::ParametricCurve;
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut a.as_bytes()).unwrap();
        let b = crate::io::brep::read(&mut model, &mut b.as_bytes()).unwrap();
        let s = section(&mut model, a, b).unwrap();
        let edges: Vec<(Point, Point, Point)> = model
            .unique_subshapes(s, NodeKind::Edge)
            .into_iter()
            .map(|e| {
                let (c, t0, t1) = model.edge_curve(e).unwrap();
                let (v0, v1) = model.edge_vertices(e);
                (model.point(v0), model.point(v1), c.point(0.5 * (t0 + t1)))
            })
            .collect();
        for (i, &(a0, a1, am)) in edges.iter().enumerate() {
            for &(b0, b1, bm) in &edges[i + 1..] {
                let near = |p: Point, q: Point| p.distance(q) <= 2e-3;
                let ends = (near(a0, b0) && near(a1, b1)) || (near(a0, b1) && near(a1, b0));
                assert!(!(ends && near(am, bm) && a0.distance(a1) > 1e-2), "the stretch {a0:?} to {a1:?} twice");
            }
        }
    }

    /// The two end walls of a thickened plate (a cell of a gyroid lattice) along the cell's
    /// face x = 10, and those of its image mirrored across x = 10 and y = 5 (the next cell's):
    /// both hold the plate's edge in that face, one B-spline, the image's the same curve run
    /// the other way. Vertices found on either edge where sections end on it split the other
    /// too; taken for different curves, each was split at its own, and the stretch between
    /// two such vertices a hair apart was covered by a piece of each, leaving a loop.
    #[test]
    fn mirrored_walls_share_their_edge_run_the_other_way() {
        section_covering_nothing_twice(include_str!("tests/mirrored_walls_reversed_edge_a.brep"), include_str!("tests/mirrored_walls_reversed_edge_b.brep"));
    }

    /// As above, the walls of two other cells of the lattice along their shared face: a
    /// section ending where the walls meet nearly tangent met the two edges there (one of
    /// each, coincident) a little apart, and each took its own vertex; the other edge, with a
    /// vertex of its own that near, did not take it, and the edges were split a hair apart.
    /// One vertex for both.
    #[test]
    fn mirrored_walls_meet_their_coincident_edges_at_one_vertex() {
        section_covering_nothing_twice(include_str!("tests/mirrored_walls_section_ends_a.brep"), include_str!("tests/mirrored_walls_section_ends_b.brep"));
    }

    /// An elliptic base, its sloped side ending on the ellipse of its top, and two planks
    /// laid on the top: one cut to the ellipse (its upright wall rising from an arc of it),
    /// one running past the base. The section of the sloped side with the plank's wall,
    /// traced along the arc, is the base's own edge there, and is taken for it (the two
    /// between the same vertices, within a little of their length): covering it twice, the
    /// side made no region, and the union lost it.
    #[test]
    fn planks_on_an_elliptic_base_cut_to_its_edge() {
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut include_str!("tests/elliptic_base_a.brep").as_bytes()).unwrap();
        let b = crate::io::brep::read(&mut model, &mut include_str!("tests/planks_cut_to_the_ellipse_b.brep").as_bytes()).unwrap();
        let common = intersection(&mut model, a, b).unwrap();
        let common = properties::volume(&model, common).mass;
        let expected = properties::volume(&model, a).mass + properties::volume(&model, b).mass - common;
        let u = union(&mut model, a, b).unwrap();
        assert!(validity::check(&model, u).is_valid());
        assert!((properties::volume(&model, u).mass - expected).abs() <= 1e-6 * expected, "{} vs {expected}", properties::volume(&model, u).mass);
    }

    /// A block's corner hollowed round (a quarter circle, R 2, about the z axis) and a slab
    /// set on it whose underside ends on the same quarter circle, made a rational B-spline:
    /// both curves the one circle, the edges are split only where they end, not wherever
    /// the search for crossings met them all along. The union is the two together.
    #[test]
    fn a_slab_on_a_hollowed_corner_along_its_arc_made_a_b_spline() {
        use crate::geometry::curve::{BSplineCurve, Curve};
        let mut model = Model::new();
        let block = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 2.0, 1.0)).build(&mut model).unwrap();
        let round = Cylinder::new(2.0, 3.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, -1.0), Direction::Z)).build(&mut model).unwrap();
        let a = difference(&mut model, block, round).unwrap();
        let w = 0.5f64.sqrt();
        let arc = BSplineCurve::new(vec![Point::new(2.0, 0.0, 1.0), Point::new(2.0, 2.0, 1.0), Point::new(0.0, 2.0, 1.0)], vec![0.0, 1.0], vec![3, 3], 2).with_weights(vec![1.0, w, 1.0]);
        let corners = [(0.0, 2.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0), (2.0, 0.0)].map(|(x, y)| Point::new(x, y, 1.0));
        let mut sides = vec![crate::make::edge_from_curve(&mut model, Curve::BSpline(arc)).unwrap()];
        for i in 0..4 {
            sides.push(crate::make::edge_between(&mut model, corners[i], corners[i + 1]).unwrap());
        }
        let wire = crate::make::wire(&mut model, &sides).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let b = crate::sweep::extrude(&mut model, face, Vector::new(0.0, 0.0, 1.0)).unwrap().shape;
        let expected = properties::volume(&model, a).mass + properties::volume(&model, b).mass;
        let u = union(&mut model, a, b).unwrap();
        assert!(validity::check(&model, u).is_valid());
        assert!((properties::volume(&model, u).mass - expected).abs() <= 1e-6 * expected, "{} vs {expected}", properties::volume(&model, u).mass);
        let most = model.unique_subshapes(u, NodeKind::Face).iter().map(|&f| model.unique_subshapes(f, NodeKind::Edge).len()).max().unwrap();
        assert!(most <= 8, "a face of {most} edges");
    }

    /// A blade lofted through ellipses and a cap on its end: a half ellipse turned half a turn
    /// about the axis it ends on, its flat face the blade's end, both bounded by the one
    /// ellipse. A vertex found on one operand's ellipse did not split the other's (not taken
    /// for one curve), and the cap's edges along it, made of the blade's, had their pcurves
    /// fitted on the other branch of the surface (past the axis, the half ellipse turned half
    /// a turn passes through the same points) near the poles: the cap lost, the union
    /// invalid. The union is the two together.
    #[test]
    fn a_revolved_cap_on_a_lofted_blade_end() {
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut include_str!("tests/lofted_blade_a.brep").as_bytes()).unwrap();
        let b = crate::io::brep::read(&mut model, &mut include_str!("tests/revolved_cap_on_blade_b.brep").as_bytes()).unwrap();
        let expected = properties::volume(&model, a).mass + properties::volume(&model, b).mass;
        let u = union(&mut model, a, b).unwrap();
        assert!(validity::check(&model, u).is_valid());
        assert!((properties::volume(&model, u).mass - expected).abs() <= 1e-6 * expected, "{} vs {expected}", properties::volume(&model, u).mass);
        assert_eq!(model.unique_subshapes(u, NodeKind::Face).len(), 8);
        // (Along the ellipse the operands' edges split each other only where they end: a few
        // edges there, not hundreds.)
        assert!(model.unique_subshapes(u, NodeKind::Edge).len() < 40, "{} edges", model.unique_subshapes(u, NodeKind::Edge).len());
        // The blade again with the capped one: within it, its end face inside, though found
        // past the cap's axis where the ray meets the cap near its pole (the cap's
        // parameters there on the other branch, off the face, the hit missed).
        let again = union(&mut model, a, u).unwrap();
        assert!(validity::check(&model, again).is_valid());
        assert!((properties::volume(&model, again).mass - expected).abs() <= 1e-6 * expected, "{} vs {expected}", properties::volume(&model, again).mass);
        let common = intersection(&mut model, a, u).unwrap();
        let blade = properties::volume(&model, a).mass;
        assert!((properties::volume(&model, common).mass - blade).abs() <= 1e-6 * blade, "{} vs {blade}", properties::volume(&model, common).mass);
    }

    /// At a corner of two cells of a gyroid lattice: one plate's offset face and wall along
    /// the cell's top face, and the walls of the neighbouring cell's plate (its mirror image)
    /// along the face below and their shared face. The latter's short edge across its walls
    /// at the corner crosses the offset face 1.8e-4 off the face's edge (both 0.05 long
    /// there), within the face's loose tolerance: taken for crossing that edge, it was split
    /// 2.1e-4 off its own curve, and the sections through the corner ended at that point
    /// twice, a lens between them.
    #[test]
    fn mirrored_plates_corner_edges_cross_the_face_not_its_edge() {
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut include_str!("tests/mirrored_plates_corner_a.brep").as_bytes()).unwrap();
        let b = crate::io::brep::read(&mut model, &mut include_str!("tests/mirrored_plates_corner_b.brep").as_bytes()).unwrap();
        let s = section(&mut model, a, b).unwrap();
        let mut ends: Vec<(u32, u32)> = Vec::new();
        for e in model.unique_subshapes(s, NodeKind::Edge) {
            if properties::length(&model, e).mass <= 1e-2 {
                continue;
            }
            let (v0, v1) = model.edge_vertices(e);
            let key = (v0.node.index.min(v1.node.index), v0.node.index.max(v1.node.index));
            assert!(!ends.contains(&key), "two section edges from {:?} to {:?}", model.point(v0), model.point(v1));
            ends.push(key);
        }
    }

    /// An end wall of a thickened plate (a cell of a Schwarz lattice) and its mirror image
    /// across y = 1.99, hinged on the arc they share there and twisting through the plane:
    /// their section is the arc and two branches, where the walls cross the plane, from the
    /// arc out to the offset edge. Towards the arc the walls cross ever more slightly (both
    /// lie in the plane along it), and the second branch's last stretch was taken for the
    /// walls only touching and dropped; and where it ends on the offset edge the point there
    /// lies off the traced curve by more than the tolerance (the walls crossing at a sine of
    /// 0.03): the branch was lost, and the plate's union with its image left open.
    #[test]
    fn mirrored_hinged_walls_crossing_the_plane_at_two_branches() {
        use crate::math::{Direction, Transform};
        let mut model = Model::new();
        let a = crate::io::brep::read(&mut model, &mut include_str!("tests/mirrored_walls_hinged_a.brep").as_bytes()).unwrap();
        let mirror = Transform::mirror_plane(Point::new(0.0, 1.99, 0.0), Direction::Y);
        let b = crate::modify::transform(&mut model, a, &mirror, Default::default()).unwrap();
        let s = section(&mut model, a, b).unwrap();
        let branches: Vec<(Point, Point)> = model
            .unique_subshapes(s, NodeKind::Edge)
            .into_iter()
            .filter(|&e| (properties::length(&model, e).mass - 0.45).abs() < 1e-3)
            .map(|e| {
                let (v0, v1) = model.edge_vertices(e);
                (model.point(v0), model.point(v1))
            })
            .collect();
        assert_eq!(branches.len(), 2, "{branches:?}");
        for end in [Point::new(0.5405, 1.99, 0.0702), Point::new(0.1042, 1.99, 0.5349)] {
            assert!(branches.iter().any(|(p, q)| p.distance(end) < 1e-3 || q.distance(end) < 1e-3), "{end:?}: {branches:?}");
        }
    }

    /// Two lines crossing at 0.03 radians, from points on them 0.5 short of the crossing:
    /// their closest points are the crossing itself (alternate projections from there close
    /// in by only cos² of the angle each round, leaving them well apart).
    #[test]
    fn edges_crossing_at_a_slight_angle_met_where_they_cross() {
        use crate::geometry::curve::{Curve, Line};
        let line = |angle: f64| Curve::Line(Line { origin: Point::new(1.0, 2.0, 3.0), direction: Direction::new(angle.cos(), angle.sin(), 0.0) });
        let (a, b) = (line(0.0), line(0.03));
        let (u, s) = intersect::closest_points(&a, (-10.0, 10.0), &b, (-10.0, 10.0), (-0.5, -0.5));
        assert!(u.abs() < 1e-9 && s.abs() < 1e-9, "{u} {s}");
    }

    /// The section of `a` and `b` mirrored across the plane `normal` · (p - `origin`) = 0:
    /// its edges' lengths and ends.
    fn mirrored_section(a: &str, b: &str, origin: Point, normal: Direction) -> Vec<(f64, Point, Point)> {
        use crate::math::Transform;
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let (a, b) = (read(a), read(b));
        let b = crate::modify::transform(&mut model, b, &Transform::mirror_plane(origin, normal), Default::default()).unwrap();
        let s = section(&mut model, a, b).unwrap();
        model
            .unique_subshapes(s, NodeKind::Edge)
            .into_iter()
            .map(|e| {
                let (v0, v1) = model.edge_vertices(e);
                (properties::length(&model, e).mass, model.point(v0), model.point(v1))
            })
            .collect()
    }

    /// A countersunk screw's head and its threaded shank (cq_warehouse's iso14581 M6, made by
    /// CadQuery): the thread's crest, a strip of the cylinder of radius 3, reaches into the
    /// head across its base, whose circle (radius 3 too) lies on the crest for 2.4 of its
    /// length, crossing the strip at the thread's lead. The crest is cut along that stretch
    /// (the circle lying along the crest, off its edges beyond): its piece inside the head
    /// was kept whole, the shell left open (area 780.13). OCCT's union: 778.9465.
    #[test]
    fn a_thread_crest_reaching_into_a_head_across_its_base_circle() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let head = read(include_str!("tests/countersunk_head_a.brep"));
        let shank = read(include_str!("tests/countersunk_threaded_shank_b.brep"));
        let result = union(&mut model, head, shank).unwrap();
        assert!(crate::validity::check(&model, result).is_valid());
        let area = crate::properties::area(&model, result).mass;
        assert!((area - 778.9465).abs() < 1e-2, "{area}");
    }

    /// A thickened Schwarz plate's wall lying flat in a face of its cell (z = 0), and its
    /// mirror image across x = 1.99, continuing it from the corner edge they share: their
    /// section is that corner, in the plane. Their corner edges leave it a hair apart, each
    /// within the classifying distance of the other's wall for a while, weaving across it:
    /// each was met by the other's wall at every sample there, 0.0035 apart, splitting both
    /// walls into a staircase of slivers. Lying along it there and turning away, never
    /// crossing it, they meet it nowhere within that stretch.
    #[test]
    fn mirrored_flat_walls_meet_along_their_corner_alone() {
        let wall = include_str!("tests/thickened_plate_flat_wall.brep");
        let edges = mirrored_section(wall, wall, Point::new(1.99, 0.0, 0.0), Direction::X);
        assert_eq!(edges.len(), 1, "{edges:?}");
        assert!((edges[0].0 - 0.44994).abs() < 1e-4, "{edges:?}");
    }

    /// The same flat wall and the plate's offset face mirrored across y = 1.99: the offset
    /// faces continue each other across that plane, and the wall's offset edge passes along
    /// the mirrored one within 1.5e-5 for a little way near the plane, crossing it there.
    /// That piece of their section, 0.0011 long, runs along the wall's own edge through both
    /// its ends (within the faces' tolerance): the edge itself, not another cut a hair beside
    /// it, which left a sliver face of no area, inside out.
    #[test]
    fn a_section_along_a_walls_own_edge_is_that_edge() {
        let edges = mirrored_section(include_str!("tests/thickened_plate_flat_wall.brep"), include_str!("tests/mirrored_offset_faces_a.brep"), Point::new(0.0, 1.99, 0.0), Direction::Y);
        assert!(edges.is_empty(), "{edges:?}");
    }

    /// A disc with screw holes round its rim, their axes on its side, each blind hole's drill
    /// tip a cone whose apex lies on the side: cut by the next hole's tool, away from them.
    /// The cones' halves end at the apex on a degenerate edge whose pcurve stops a little
    /// apart in u from where the section beside it reaches the apex: their ends are one
    /// point there, and the halves were lost, the shell left open (OCCT's volume).
    #[test]
    fn a_disc_with_drill_tips_on_its_side_cut_again() {
        let mut model = Model::new();
        let mut read = |text: &str| crate::io::brep::read(&mut model, &mut text.as_bytes()).unwrap();
        let disc = read(include_str!("tests/disk_with_rim_holes_a.brep"));
        let tool = read(include_str!("tests/disk_with_rim_holes_b.brep"));
        let result = difference(&mut model, disc, tool).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let volume = properties::volume(&model, result).mass;
        assert!((volume - 355415.434).abs() < 1e-2, "{volume}");
    }

    /// The solid paraboloid x² + y² ≤ 8z, z ≤ `height`: the parabola of the xz plane from its
    /// vertex at the origin to the top's rim, closed by the top's radius and the axis, turned
    /// a whole turn about the z axis.
    fn paraboloid(model: &mut Model, height: f64) -> Shape {
        use crate::geometry::curve::{Curve, Parabola};
        let frame = Frame::new(Point::new(0.0, 0.0, 0.0), Direction::Y, Direction::Z);
        let rim = (8.0 * height).sqrt();
        let side = crate::make::edge_from_curve_range(model, Curve::Parabola(Parabola::new(frame, 2.0)), 0.0, rim).unwrap();
        let top = crate::make::edge_between(model, Point::new(rim, 0.0, height), Point::new(0.0, 0.0, height)).unwrap();
        let axis = crate::make::edge_between(model, Point::new(0.0, 0.0, height), Point::new(0.0, 0.0, 0.0)).unwrap();
        let wire = crate::make::wire(model, &[side, top, axis]).unwrap();
        let section = crate::make::planar_face(model, wire).unwrap();
        crate::sweep::revolve(model, section, crate::math::Axis::Z, 2.0 * PI).unwrap().shape
    }

    /// Checks the result is valid, with the expected volume and area to `relative`.
    fn measures(model: &Model, shape: Shape, volume: f64, area: f64, relative: f64) {
        let report = validity::check(model, shape);
        assert!(report.is_valid(), "{report:?}");
        let v = properties::volume(model, shape).mass;
        let a = properties::area(model, shape).mass;
        assert!((v - volume).abs() <= relative * volume, "volume {v} != {volume}");
        assert!((a - area).abs() <= relative * area, "area {a} != {area}");
    }

    /// A block (−10, −10, 0) … (10, 10, 10) united with the paraboloid x² + y² ≤ 8z, z ≤ 12.5,
    /// its vertex on the block's bottom face from inside, standing out of its top. The
    /// parabola runs on across the axis past its vertex: turned, the surface's parameters of
    /// a point of the paraboloid were found on either half, half a turn apart, and the side
    /// was lost; with the basis trimmed to the edge (as OCCT's) the union is whole:
    /// 4000 + 225π, its area 1600 + 20π + (4π/3)(29√29 − 48√6).
    #[test]
    fn a_paraboloid_tip_on_a_face_united() {
        let mut model = Model::new();
        let block = Cuboid::from_corners(Point::new(-10.0, -10.0, 0.0), Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
        let solid = paraboloid(&mut model, 12.5);
        let result = union(&mut model, block, solid).unwrap();
        let area = 1600.0 + 20.0 * PI + 4.0 * PI / 3.0 * (29.0 * 29f64.sqrt() - 48.0 * 6f64.sqrt());
        measures(&model, result, 4000.0 + 225.0 * PI, area, 1e-7);
    }

    /// The paraboloid x² + y² ≤ 8z, z ≤ 8, common with a block below the plane z = h = 10⁻⁴:
    /// its tip, 4πh², its area 8πh + (32π/3)((1 + h/2)^{3/2} − 1).
    #[test]
    fn a_paraboloid_tip_a_ten_thousandth_high() {
        let mut model = Model::new();
        let solid = paraboloid(&mut model, 8.0);
        let h = 1e-4;
        let below = Cuboid::from_corners(Point::new(-10.0, -10.0, -1.0), Point::new(10.0, 10.0, h)).build(&mut model).unwrap();
        let result = intersection(&mut model, solid, below).unwrap();
        let area = 8.0 * PI * h + 32.0 * PI / 3.0 * ((1.0 + 0.5 * h).powf(1.5) - 1.0);
        measures(&model, result, 4.0 * PI * h * h, area, 1e-7);
    }

    /// A cube (−5, −5, −10) … (5, 5, 0) common with the half-space above z = s·y, s = 10⁻⁶:
    /// the planes meet in the x axis, a millionth of a radian apart (parallel only within
    /// OCCT's 10⁻⁸), and what is common is a wedge 10 · 5 · 5s / 2 under the half of the top
    /// face where y < 0, its area 50 + 50√(1 + s²) + 50s + 25s.
    #[test]
    fn a_wedge_a_millionth_steep() {
        let mut model = Model::new();
        let s = 1e-6;
        let cube = Cuboid::from_corners(Point::new(-5.0, -5.0, -10.0), Point::new(5.0, 5.0, 0.0)).build(&mut model).unwrap();
        let plane = crate::geometry::surface::Surface::Plane(crate::geometry::surface::Plane::new(Point::new(0.0, 0.0, 0.0), Direction::new(0.0, -s, 1.0)));
        let face = crate::make::face_on_domain(&mut model, plane, -1000.0, 1000.0, -1000.0, 1000.0).unwrap();
        let above = crate::primitives::half_space(&mut model, face, Point::new(0.0, 0.0, 1.0)).unwrap();
        let result = intersection(&mut model, cube, above).unwrap();
        measures(&model, result, 125.0 * s, 50.0 + 50.0 * (1.0 + s * s).sqrt() + 75.0 * s, 1e-7);
    }

    /// Two balls of radius 5 whose centres are d = 10 − 10⁻⁶ apart: their lens, a millionth
    /// thick, lies everywhere within the classifying distance of the other ball, but two
    /// spheres not one meet only along a circle: π(4r + d)(2r − d)²/12, its area
    /// 4πr(r − d/2).
    #[test]
    fn a_lens_a_millionth_thick() {
        let mut model = Model::new();
        let (r, d) = (5.0, 9.999999);
        let first = Sphere::new(r).build(&mut model).unwrap();
        let second = Sphere::new(r).with_frame(Frame::from_z(Point::new(d, 0.0, 0.0), Direction::Z)).build(&mut model).unwrap();
        let result = intersection(&mut model, first, second).unwrap();
        let thickness: f64 = 2.0 * r - d;
        measures(&model, result, PI * (4.0 * r + d) * thickness * thickness / 12.0, 4.0 * PI * r * (r - 0.5 * d), 1e-6);
    }

    /// A ball-end mill's raster: a plate (−10, −10, −20) … (20, 10, 0) less a cylinder of
    /// radius 3 ending at c = (0, 0, −1) at slope s1, a ball of radius 3 about c, and a
    /// cylinder of radius 3 starting at c at slope s2. Each cylinder is tangent to the ball
    /// along a great circle; the two circles cross where the three surfaces touch, and the
    /// cylinders meet past them within a millionth of each other: there a sliver of the
    /// first, inside the second, was taken as lying on it (the result had a void outside
    /// any solid), and the lune of the ball left between the circles, its boundary's
    /// pcurve swinging round near the pole, was found inside out. The volume against the
    /// tools' union integrated along x (exact on each line, from the cylinders' and the
    /// ball's equations) over an 800 × 400 grid in z and θ, y = 3 sin θ: within 10⁻⁴ of a
    /// grid ten times finer.
    #[test]
    fn a_ball_between_cylinders_tangent_to_it_along_crossing_circles() {
        // The x-interval of a cylinder (start o, unit axis (ax, 0, az), length l, radius 3)
        // on the line along x through (y, z): |ax ζ − az wx| ≤ √(9 − wy²), 0 ≤ ax wx + az ζ ≤ l.
        let cylinder = |o: [f64; 3], ax: f64, az: f64, l: f64, y: f64, z: f64| -> Option<(f64, f64)> {
            let (wy, zeta) = (y - o[1], z - o[2]);
            let h = (9.0 - wy * wy).max(0.0).sqrt();
            let (a, b) = ((zeta * ax - h) / az, (zeta * ax + h) / az);
            let (a, b) = if az > 0.0 { (a, b) } else { (b, a) };
            let (lo, hi) = (a.max(-zeta * az / ax), b.min((l - zeta * az) / ax));
            (lo < hi).then_some((lo + o[0], hi + o[0]))
        };
        let reference = |s1: f64, s2: f64| {
            let (n1, n2) = ((1.0 + s1 * s1).sqrt(), (1.0 + s2 * s2).sqrt());
            let (nz, nt, z0) = (800, 400, -4.5);
            let mut removed = 0.0;
            for i in 0..nz {
                let z = z0 * (1.0 - (i as f64 + 0.5) / nz as f64);
                for j in 0..nt {
                    let theta = PI * ((j as f64 + 0.5) / nt as f64 - 0.5);
                    let y = 3.0 * theta.sin();
                    let mut runs: Vec<(f64, f64)> = Vec::new();
                    let r = 9.0 - y * y - (z + 1.0) * (z + 1.0);
                    if r > 0.0 {
                        runs.push((-r.sqrt(), r.sqrt()));
                    }
                    runs.extend(cylinder([-5.0, 0.0, -1.0 - 5.0 * s1], 1.0 / n1, s1 / n1, 5.0 * n1, y, z));
                    runs.extend(cylinder([0.0, 0.0, -1.0], 1.0 / n2, s2 / n2, 5.0, y, z));
                    runs.sort_by(|a, b| a.0.total_cmp(&b.0));
                    let (mut length, mut end) = (0.0, f64::NEG_INFINITY);
                    for (a, b) in runs {
                        if b > a.max(end) {
                            length += b - a.max(end);
                            end = b;
                        }
                    }
                    removed += length * 3.0 * theta.cos() * (PI / nt as f64) * (-z0 / nz as f64);
                }
            }
            30.0 * 20.0 * 20.0 - removed
        };
        for (s1, s2) in [(-0.03, -0.025), (-0.01, -0.005), (0.01, 0.015), (-0.02, -0.005)] {
            let mut model = Model::new();
            let plate = Cuboid::from_corners(Point::new(-10.0, -10.0, -20.0), Point::new(20.0, 10.0, 0.0)).build(&mut model).unwrap();
            let c = Point::new(0.0, 0.0, -1.0);
            let from = Point::new(-5.0, 0.0, -1.0 - 5.0 * s1);
            let first = Cylinder::new(3.0, from.distance(c)).with_frame(Frame::from_z(from, Direction::new(5.0, 0.0, 5.0 * s1))).build(&mut model).unwrap();
            let ball = Sphere::new(3.0).with_frame(Frame::from_z(c, Direction::Z)).build(&mut model).unwrap();
            let second = Cylinder::new(3.0, 5.0).with_frame(Frame::from_z(c, Direction::new(5.0, 0.0, 5.0 * s2))).build(&mut model).unwrap();
            let mut body = plate;
            for tool in [first, ball, second] {
                body = difference(&mut model, body, tool).unwrap();
            }
            let report = validity::check(&model, body);
            assert!(report.is_valid(), "{s1} {s2}: {report:?}");
            let (v, expected) = (properties::volume(&model, body).mass, reference(s1, s2));
            assert!((v - expected).abs() <= 1e-3, "{s1} {s2}: volume {v} != {expected}");
        }
    }

    /// A ball-end mill's move from (75, 20, −6) to (80, 20, −5.8) (its ball's centre), as a
    /// raster builds it: a plate (70, y0, −15) … (90, 30, 0) less the ball and the shank (the
    /// cylinder up from its centre) where it ends, and the front half of the shank's side slid
    /// along the move. That slid face's bottom, the half circle at the start slid along the
    /// move, ends on the ball's equator: the two meet along it and along a curve through the
    /// points where the slid circle's ends touch the ball, the equator crossing it there.
    /// That curve was dropped as only touching for its one point nearly tangent, and on one
    /// side cut where it passes over the equator, taken for it traced again; the plate's side
    /// y = y0 met no line of the slid face, its half circle from −π/2 to π/2 searched only
    /// from 0. With the slid face alone the removed volume is 5 (y-width) ∫ min(1, (z + 6)/0.2)
    /// dz over [−6, 0], 5.9 × 5 × the width; with all three, against the tools' union
    /// integrated along x (exact on each line) over a 1600 × 800 grid in z (split at the
    /// shank's bottom) and θ, y = 20 + 3 sin θ: within 2 × 10⁻⁴ of a grid four times finer.
    #[test]
    fn a_ball_end_mill_slid_shank_ending_on_its_ball() {
        use crate::geometry::surface::{Cylinder as CylinderSurface, Surface};
        use std::f64::consts::FRAC_PI_2;
        let h = 0.2;
        let zc = -6.0 + h;
        let reference = |y0: f64, all: bool| {
            let (za, zb) = (zc - 3.0, 0.0);
            let (t0, t1) = (((y0 - 20.0) / 3.0).max(-1.0).asin(), FRAC_PI_2);
            let (nz, nt) = (1600, 800);
            let mut removed = 0.0;
            for (lo, hi) in [(za, zc), (zc, zb)] {
                let n = ((hi - lo) / (zb - za) * nz as f64).ceil() as usize;
                for i in 0..n {
                    let z = lo + (hi - lo) * (i as f64 + 0.5) / n as f64;
                    for j in 0..nt {
                        let theta = t0 + (t1 - t0) * (j as f64 + 0.5) / nt as f64;
                        let (y, w) = (20.0 + 3.0 * theta.sin(), 3.0 * theta.cos());
                        let mut runs: Vec<(f64, f64)> = Vec::new();
                        let r = 9.0 - (y - 20.0).powi(2) - (z - zc).powi(2);
                        if all && r > 0.0 {
                            runs.push((80.0 - r.sqrt(), 80.0 + r.sqrt()));
                        }
                        if all && z >= zc {
                            runs.push((80.0 - w, 80.0 + w));
                        }
                        // (The slid face's point at height z − hs, s along the move.)
                        let (s0, s1) = (((z - 21.0) / h).max(0.0), ((z + 6.0) / h).min(1.0));
                        if s0 < s1 {
                            runs.push((75.0 + w + 5.0 * s0, 75.0 + w + 5.0 * s1));
                        }
                        runs.sort_by(|a, b| a.0.total_cmp(&b.0));
                        let (mut length, mut end) = (0.0, f64::NEG_INFINITY);
                        for (a, b) in runs {
                            if b > a.max(end) {
                                length += b - a.max(end);
                                end = b;
                            }
                        }
                        removed += length * w * (t1 - t0) / nt as f64 * (hi - lo) / n as f64;
                    }
                }
            }
            20.0 * (30.0 - y0) * 15.0 - removed
        };
        for (y0, all) in [(10.0, true), (18.0, true), (10.0, false), (18.0, false)] {
            let mut model = Model::new();
            let mut body = Cuboid::from_corners(Point::new(70.0, y0, -15.0), Point::new(90.0, 30.0, 0.0)).build(&mut model).unwrap();
            let (from, to) = (Point::new(75.0, 20.0, -6.0), Point::new(80.0, 20.0, zc));
            let side = Surface::Cylinder(CylinderSurface::new(Frame::new(from, Direction::Z, Direction::X), 3.0));
            let half = crate::make::face_on_domain(&mut model, side, -FRAC_PI_2, FRAC_PI_2, 0.0, 27.0).unwrap();
            let slid = crate::sweep::extrude(&mut model, half, to - from).unwrap().shape;
            let slid = model.explore(slid, NodeKind::Solid).next().unwrap();
            let mut tools = vec![slid];
            if all {
                let ball = Sphere::new(3.0).with_frame(Frame::from_z(to, Direction::Z)).build(&mut model).unwrap();
                let shank = Cylinder::new(3.0, 30.0).with_frame(Frame::from_z(to, Direction::Z)).build(&mut model).unwrap();
                tools = vec![ball, shank, slid];
            }
            for tool in tools {
                body = difference(&mut model, body, tool).unwrap();
            }
            let report = validity::check(&model, body);
            assert!(report.is_valid(), "{y0} {all}: {report:?}");
            let v = properties::volume(&model, body).mass;
            let expected = if all { reference(y0, true) } else { 20.0 * (30.0 - y0) * 15.0 - 5.0 * (23.0 - y0.max(17.0)) * 5.9 };
            assert!((expected - reference(y0, all)).abs() <= 1e-3, "{y0} {all}: the reference {} != {expected}", reference(y0, all));
            assert!((v - expected).abs() <= 1e-3, "{y0} {all}: volume {v} != {expected}");
        }
    }

    /// A piece of a ball-end mill's raster over the height field z = −5 + 4 sin(πx/50)
    /// sin(πy/40), its ball's centre 3 above it.
    #[derive(Clone, Copy)]
    enum Mill {
        /// A move's cylinder, from the ball's centre over one point of the field to over
        /// the other.
        Along(f64, f64, f64, f64),
        /// The ball over a point.
        Ball(f64, f64),
        /// The shank (the cylinder up from the ball's centre) over a point.
        Shank(f64, f64),
    }

    /// The ball's centre over (x, y).
    fn mill_centre(x: f64, y: f64) -> Point {
        Point::new(x, y, -2.0 + 4.0 * (PI * x / 50.0).sin() * (PI * y / 40.0).sin())
    }

    impl Mill {
        fn build(self, model: &mut Model) -> Shape {
            match self {
                Mill::Along(x0, y0, x1, y1) => {
                    let (a, b) = (mill_centre(x0, y0), mill_centre(x1, y1));
                    let d = b - a;
                    Cylinder::new(3.0, d.length()).with_frame(Frame::from_z(a, Direction::new(d.x, d.y, d.z))).build(model).unwrap()
                }
                Mill::Ball(x, y) => Sphere::new(3.0).with_frame(Frame::from_z(mill_centre(x, y), Direction::Z)).build(model).unwrap(),
                Mill::Shank(x, y) => Cylinder::new(3.0, 30.0).with_frame(Frame::from_z(mill_centre(x, y), Direction::Z)).build(model).unwrap(),
            }
        }

        /// Its y about which it lies.
        fn y(self) -> f64 {
            match self {
                Mill::Along(_, y, ..) | Mill::Ball(_, y) | Mill::Shank(_, y) => y,
            }
        }

        /// Its x-interval on the line along x through (y, z), exactly. A move's (along x,
        /// start o, axis to e, radius 3): |ax ζ − az wx| ≤ √(9 − wy²), 0 ≤ ax wx + az ζ ≤ l,
        /// taken from its end further along x.
        fn on_line(self, y: f64, z: f64) -> Option<(f64, f64)> {
            match self {
                Mill::Along(x0, y0, x1, y1) => {
                    let (o, e) = if x0 < x1 { (mill_centre(x0, y0), mill_centre(x1, y1)) } else { (mill_centre(x1, y1), mill_centre(x0, y0)) };
                    let d = e - o;
                    let l = d.length();
                    let (ax, az) = (d.x / l, d.z / l);
                    let (wy, zeta) = (y - o.y, z - o.z);
                    if wy * wy >= 9.0 {
                        return None;
                    }
                    let h = (9.0 - wy * wy).sqrt();
                    let (a, b) = ((zeta * ax - h) / az, (zeta * ax + h) / az);
                    let (a, b) = if az > 0.0 { (a, b) } else { (b, a) };
                    let (lo, hi) = (a.max(-zeta * az / ax), b.min((l - zeta * az) / ax));
                    (lo < hi).then_some((lo + o.x, hi + o.x))
                }
                Mill::Ball(x, y0) => {
                    let o = mill_centre(x, y0);
                    let r = 9.0 - (y - o.y).powi(2) - (z - o.z).powi(2);
                    (r > 0.0).then(|| (o.x - r.sqrt(), o.x + r.sqrt()))
                }
                Mill::Shank(x, y0) => {
                    let o = mill_centre(x, y0);
                    let r = 9.0 - (y - o.y).powi(2);
                    (r > 0.0 && z >= o.z).then(|| (o.x - r.sqrt(), o.x + r.sqrt()))
                }
            }
        }
    }

    /// A plate (x0, y0, −20) … (x1, y1, 0) less the `mills` in turn, the last a move's
    /// cylinder: the result valid, and the volume the last removes against the move, below
    /// the plate's top, less the others integrated along x (exact on each line) over a grid
    /// of 800 in z (split at the shanks' bottoms) by 2400 in y (split where the others' sides
    /// end, y = m + w sin θ on each part): within 5 × 10⁻⁵ of a grid twice as fine for the
    /// raster pieces below.
    fn mill_raster(plate: [f64; 4], mills: &[Mill]) {
        let (&last, others) = mills.split_last().unwrap();
        let Mill::Along(x0, y0, x1, y1) = last else { panic!("the last a move") };
        let reference = || {
            let (a, b) = (mill_centre(x0, y0), mill_centre(x1, y1));
            let (z0, z1) = (a.z.min(b.z) - 3.1, (a.z.max(b.z) + 3.1).min(0.0));
            let (nz, ny) = (800, 2400);
            let mut zs = vec![z0, z1];
            zs.extend(others.iter().filter_map(|m| match *m {
                Mill::Shank(x, y) => Some(mill_centre(x, y).z).filter(|&z| z > z0 && z < z1),
                _ => None,
            }));
            zs.sort_by(f64::total_cmp);
            let mut ys = vec![y0 - 3.0, y0 + 3.0];
            ys.extend(others.iter().flat_map(|m| [m.y() - 3.0, m.y() + 3.0]).filter(|&y| y > y0 - 3.0 && y < y0 + 3.0));
            ys.sort_by(f64::total_cmp);
            ys.dedup();
            let mut removed = 0.0;
            for (zlo, zhi) in zs.iter().zip(&zs[1..]) {
                let n = ((zhi - zlo) / (z1 - z0) * nz as f64).ceil() as usize;
                for i in 0..n {
                    let z = zlo + (zhi - zlo) * (i as f64 + 0.5) / n as f64;
                    for (ylo, yhi) in ys.iter().zip(&ys[1..]) {
                        let (m, w) = (0.5 * (ylo + yhi), 0.5 * (yhi - ylo));
                        let k = (ny as f64 * (yhi - ylo) / 6.0).round() as usize;
                        for j in 0..k {
                            let theta = PI * ((j as f64 + 0.5) / k as f64 - 0.5);
                            let y = m + w * theta.sin();
                            let Some((t0, t1)) = last.on_line(y, z) else { continue };
                            let mut runs: Vec<(f64, f64)> = others.iter().filter_map(|m| m.on_line(y, z)).collect();
                            runs.sort_by(|a, b| a.0.total_cmp(&b.0));
                            let (mut length, mut end) = (t1 - t0, t0);
                            for (a, b) in runs {
                                let (a, b) = (a.max(end), b.min(t1));
                                if b > a {
                                    length -= b - a;
                                    end = b;
                                }
                            }
                            removed += length * w * theta.cos() * (PI / k as f64) * (zhi - zlo) / n as f64;
                        }
                    }
                }
            }
            removed
        };
        let mut model = Model::new();
        let mut body = Cuboid::from_corners(Point::new(plate[0], plate[1], -20.0), Point::new(plate[2], plate[3], 0.0)).build(&mut model).unwrap();
        for m in others {
            let tool = m.build(&mut model);
            body = difference(&mut model, body, tool).unwrap();
        }
        let before = properties::volume(&model, body).mass;
        let tool = last.build(&mut model);
        let result = difference(&mut model, body, tool).unwrap();
        let report = validity::check(&model, result);
        assert!(report.is_valid(), "{report:?}");
        let (removed, expected) = (before - properties::volume(&model, result).mass, reference());
        assert!((removed - expected).abs() <= 1e-3, "removed {removed} != {expected}");
    }

    /// A ball-end mill's raster (see [`Mill`]): the shank at (85, 30), the move from there
    /// to (80, 30), the ball and the shank at (80, 32), then the move on from there to
    /// (85, 32), its cylinder tangent to the ball along the circle it starts from. The first
    /// move's end cap, cut by the ball, is left a sliver between the arc the ball cuts and
    /// its circle, the two meeting nearly tangent; the last move's side crosses the sliver,
    /// from the one to the other, within the classifying distance of the arc for half the
    /// way. That piece of their section was dropped, located by its middle on the cap's
    /// boundary, and with it the side's region inside the plate: the move removed nothing.
    #[test]
    fn a_ball_end_mill_move_across_a_sliver_of_the_cap_before() {
        use Mill::*;
        mill_raster([70.0, 24.0, 92.0, 38.0], &[Shank(85.0, 30.0), Along(85.0, 30.0, 80.0, 30.0), Ball(80.0, 32.0), Shank(80.0, 32.0), Along(80.0, 32.0, 85.0, 32.0)]);
    }

    /// As [`a_ball_end_mill_move_across_a_sliver_of_the_cap_before`], the moves from (90, 14)
    /// to (85, 14) and from (80, 16) to (85, 16) and the balls at their ends, the shank at
    /// (85, 16), then the move on to (90, 16): the sliver it crosses, of the first ball
    /// between the second and the first move's side, nowhere wider than the classifying
    /// distance, no point of the piece across it clear of the boundary to locate it by. The
    /// side of the arc its middle is on tells, by the arc's own curve.
    #[test]
    fn a_ball_end_mill_move_across_a_sliver_narrower_than_the_classifying_distance() {
        use Mill::*;
        mill_raster([76.0, 8.0, 98.0, 22.0], &[Along(90.0, 14.0, 85.0, 14.0), Ball(85.0, 14.0), Along(80.0, 16.0, 85.0, 16.0), Ball(85.0, 16.0), Shank(85.0, 16.0), Along(85.0, 16.0, 90.0, 16.0)]);
    }

    /// A ball-end mill's raster (see [`Mill`]): the move from (80, 18) to (75, 18), the move
    /// from (70, 20) to (75, 20), the ball and the shank there, then the move on to
    /// (80, 20). Its end cap is square to the first move's cylinder but for 4.8 × 10⁻⁴
    /// radians: they meet in an ellipse, not the circle about the cylinder's axis taken
    /// for it within a thousandth of a radian (straying from the cap's plane by 1.4 × 10⁻³),
    /// which missed the point where the move's side meets them both, ended no section on
    /// the cap and left its side one region outside the plate.
    #[test]
    fn a_ball_end_mill_move_ending_on_a_cap_nearly_square_to_a_cylinder() {
        use Mill::*;
        mill_raster([66.0, 12.0, 88.0, 26.0], &[Along(80.0, 18.0, 75.0, 18.0), Along(70.0, 20.0, 75.0, 20.0), Ball(75.0, 20.0), Shank(75.0, 20.0), Along(75.0, 20.0, 80.0, 20.0)]);
    }

    /// A ball-end mill's raster (see [`Mill`]) by the plate's top, the ball's centre above it
    /// at (40, 20): the moves from (45, 24) to (50, 24) and from (40, 24) to (45, 24), the
    /// move from (45, 20) to (50, 20), then the move before it, from (40, 20) to (45, 20),
    /// ending where that one starts, their axes 1.3 degrees apart. Its side crosses the
    /// other's (equal cylinders, along an ellipse nearly square to both) at a sine of 0.017,
    /// and so the edge where the other's side meets the first row's: one of that edge's
    /// samples lay 6.7 × 10⁻⁸ from the side, within the tolerance, and was taken for where
    /// the edge meets it, 3.9 × 10⁻⁶ along the edge from where it does. The ellipse passed
    /// the vertex put there by more than the classifying distance and was not ended on it, so
    /// its piece across the other move's side was not made; the sections on the move's side
    /// did not close, it stayed one region outside the plate, and the move removed nothing.
    #[test]
    fn a_ball_end_mill_move_crossing_the_next_ones_side_beside_a_sample_of_its_edge() {
        use Mill::*;
        mill_raster([36.0, 14.0, 54.0, 28.0], &[Along(45.0, 24.0, 50.0, 24.0), Along(40.0, 24.0, 45.0, 24.0), Along(45.0, 20.0, 50.0, 20.0), Along(40.0, 20.0, 45.0, 20.0)]);
    }

    /// The volume inside both a ring about z (centre circle of radius `big` at height `z0`;
    /// the torus of tube radius `tori.0` less that of `tori.1`, above `z_min`) and a tube along
    /// y through (`big`, ·, `z0`) (radius `tube.0` less `tube.1`, from `ys.0` to `ys.1`): along
    /// each line in x the exact length inside both, by the midpoint rule across y and z (in
    /// stretches between the heights where the lengths jump).
    fn ring_tube_overlap(big: f64, z0: f64, tori: (f64, f64), z_min: f64, tube: (f64, f64), ys: (f64, f64)) -> f64 {
        // (The x > 0 half of a torus's section along the line: the tube lies there.)
        let torus = |r: f64, y: f64, dz: f64| -> Option<(f64, f64)> {
            let w = (r * r - dz * dz).sqrt();
            let hi = (big + w) * (big + w) - y * y;
            (w.is_finite() && hi > 0.0).then(|| (((big - w) * (big - w) - y * y).max(0.0).sqrt(), hi.sqrt()))
        };
        let disc = |r: f64, dz: f64| -> Option<(f64, f64)> {
            let w = (r * r - dz * dz).sqrt();
            w.is_finite().then(|| (big - w, big + w))
        };
        let overlap = |a: Option<(f64, f64)>, b: Option<(f64, f64)>| match (a, b) {
            (Some(a), Some(b)) => (a.1.min(b.1) - a.0.max(b.0)).max(0.0),
            _ => 0.0,
        };
        let length = |y: f64, z: f64| {
            let dz = z - z0;
            if z < z_min {
                return 0.0;
            }
            let (t1, t0) = (torus(tori.0, y, dz), torus(tori.1, y, dz));
            let (c1, c0) = (disc(tube.0, dz), disc(tube.1, dz));
            overlap(t1, c1) - overlap(t1, c0) - overlap(t0, c1) + overlap(t0, c0)
        };
        let nodes = |breaks: &[f64], lo: f64, hi: f64, n: usize| -> Vec<(f64, f64)> {
            let mut b: Vec<f64> = breaks.iter().copied().filter(|&x| x > lo && x < hi).chain([lo, hi]).collect();
            b.sort_by(f64::total_cmp);
            b.windows(2)
                .flat_map(|w| {
                    let k = ((n as f64 * (w[1] - w[0]) / (hi - lo)).round() as usize).max(4);
                    let h = (w[1] - w[0]) / k as f64;
                    (0..k).map(move |i| (w[0] + h * (i as f64 + 0.5), h))
                })
                .collect()
        };
        let reach = tori.0.min(tube.0);
        let zs = nodes(&[z0 - tori.1, z0 + tori.1, z0 - tube.1, z0 + tube.1, z_min], z0 - reach, z0 + reach, 1200);
        let ys = nodes(&[], ys.0.min(ys.1), ys.0.max(ys.1), 2400);
        ys.iter().map(|&(y, hy)| zs.iter().map(|&(z, hz)| length(y, z) * hz).sum::<f64>() * hy).sum()
    }

    /// `a` op `b` for each operation, valid, against the volumes from those of `a`, `b` and
    /// their overlap `common` (within the quadrature's error, a two-thousandth of it).
    fn check_against_overlap(model: &mut Model, a: Shape, b: Shape, (va, vb, common): (f64, f64, f64), what: &str) {
        let tolerance = 5e-4 * common + 1e-6 * (va + vb);
        for (name, expected) in [("union", va + vb - common), ("intersection", common), ("difference", va - common), ("reverse difference", vb - common)] {
            let r = match name {
                "union" => union(model, a, b),
                "intersection" => intersection(model, a, b),
                "difference" => difference(model, a, b),
                _ => difference(model, b, a),
            }
            .unwrap();
            let report = validity::check(model, r);
            assert!(report.is_valid(), "{what}, {name}: {report:?}");
            let v = properties::volume(model, r).mass;
            assert!((v - expected).abs() <= tolerance, "{what}, {name}: volume {v} against {expected}");
        }
    }

    /// A torus and a cylinder of its tube's radius whose axis is tangent to its centre
    /// circle, the cylinder ending in the meridian plane through the point of tangency (as
    /// in a compressor's housing): the two touch along that meridian circle, the cylinder's
    /// end edge, and the branch of their section where the cylinder runs into the torus
    /// crosses that circle at its top and bottom. The surfaces part there only to fourth
    /// order, and the traced branch passed the crossings 1e-4 off: they were not found, the
    /// branch was dropped and neither face was split (the results invalid, the difference
    /// the whole torus). Now the crossings are where the surfaces' curvatures across the
    /// circle change order. With the torus's seam on that meridian and off it, the cylinder
    /// either way along its axis.
    #[test]
    fn a_torus_and_a_cylinder_touching_along_a_meridian() {
        let (big, r, z0, length) = (12.8, 3.315, 3.7227, 16.885);
        for (turn, way) in [(0.0f64, -1.0), (37f64.to_radians(), -1.0), (0.0, 1.0), (-71f64.to_radians(), 1.0)] {
            let mut model = Model::new();
            let frame = Frame::new(Point::new(0.0, 0.0, z0), Direction::Z, Direction::new(turn.cos(), turn.sin(), 0.0));
            let torus = crate::primitives::Torus::new(big, r).with_frame(frame).build(&mut model).unwrap();
            let (start, axis) = if way < 0.0 { (-length, Direction::Y) } else { (length, -Direction::Y) };
            let pipe = cylinder(&mut model, (big, start, z0), axis, r, length);
            let common = ring_tube_overlap(big, z0, (r, 0.0), f64::NEG_INFINITY, (r, 0.0), (0.0, start));
            let volumes = (2.0 * PI * PI * big * r * r, PI * r * r * length, common);
            check_against_overlap(&mut model, torus, pipe, volumes, &format!("seam turned {turn}, cylinder {way}"));
        }
    }

    /// The same contact where a tube (a cylinder bored through) meets tori: one of its bore's
    /// radius (touching the bore along the meridian circle, the bore's end edge), one of its
    /// outside's, and the shell between them cut down to its top (the half of a compressor's
    /// volute). Also a smaller ring, and the tube's end a hair (within the tolerance) either
    /// side of the meridian plane, its end edge lying on the torus nearly tangent.
    #[test]
    fn a_tube_against_tori_touching_it_along_a_meridian() {
        let cases = [(12.8, 3.315, 0.3377, 3.7227, 16.885, 0.0, 0.0), (12.8, 3.315, 0.3377, 3.7227, 16.885, 0.4, 1e-3), (12.8, 3.315, 0.3377, 3.7227, 16.885, -0.4, -1e-3), (8.0, 2.5, 0.4, 2.6, 10.0, 0.0, 0.0)];
        for (big, r, wall, z0, length, turn, shift) in cases {
            let mut model = Model::new();
            let bore = r - wall;
            let frame = Frame::new(Point::new(0.0, 0.0, z0), Direction::Z, Direction::new(f64::cos(turn), f64::sin(turn), 0.0));
            let torus = |model: &mut Model, tube: f64| crate::primitives::Torus::new(big, tube).with_frame(frame).build(model).unwrap();
            let outside = cylinder(&mut model, (big, shift - length, z0), Direction::Y, r, length);
            let inside = cylinder(&mut model, (big, shift - length, z0), Direction::Y, bore, length);
            let pipe = difference(&mut model, outside, inside).unwrap();
            let ys = (shift - length, shift);
            let v_pipe = PI * (r * r - bore * bore) * length;
            let what = |name: &str| format!("ring {big}, turned {turn}, shifted {shift}: {name}");
            for tube in [bore, r] {
                let t = torus(&mut model, tube);
                let common = ring_tube_overlap(big, z0, (tube, 0.0), f64::NEG_INFINITY, (r, bore), ys);
                check_against_overlap(&mut model, pipe, t, (v_pipe, 2.0 * PI * PI * big * tube * tube, common), &what(&format!("torus {tube}")));
            }
            // The shell between the tori above z = r / 2, against the cylinder whole.
            let (outer, inner) = (torus(&mut model, r), torus(&mut model, bore));
            let shell = difference(&mut model, outer, inner).unwrap();
            let z_min = 0.5 * r;
            let slab = cuboid(&mut model, (-2.0 * big, -2.0 * big, z_min), (4.0 * big, 4.0 * big, 2.0 * z0));
            let top = intersection(&mut model, shell, slab).unwrap();
            // (Pappus: the area of each disc above the cut turned about the axis at `big`.)
            let above = |radius: f64| {
                let d = z0 - z_min;
                PI * radius * radius - (radius * radius * (d / radius).acos() - d * (radius * radius - d * d).sqrt())
            };
            let v_top = 2.0 * PI * big * (above(r) - above(bore));
            let common = ring_tube_overlap(big, z0, (r, bore), z_min, (r, 0.0), ys);
            check_against_overlap(&mut model, top, outside, (v_top, PI * r * r * length, common), &what("shell's top"));
        }
    }

    /// A tube (a cylinder bored through, along y from `-length` to 0) cut where a torus of its
    /// bore's radius runs into it, moved by `moved`; and the top of the shell between that
    /// torus and one of the tube's outside radius, with the tube's cylinder cut from it: a
    /// compressor's outlet and the half of its volute it joins. Unmoved, the two touch along
    /// pieces of that cylinder and of the torus, and across the tube's end, and each operand
    /// has its own copy of the section of the cylinder with the torus.
    fn outlet_and_volute(model: &mut Model, (big, r, wall, z0, length): (f64, f64, f64, f64, f64), turn: f64, moved: Vector) -> (Shape, Shape) {
        let bore = r - wall;
        let frame = Frame::new(Point::new(0.0, 0.0, z0), Direction::Z, Direction::new(f64::cos(turn), f64::sin(turn), 0.0));
        let torus = |model: &mut Model, tube: f64| crate::primitives::Torus::new(big, tube).with_frame(frame).build(model).unwrap();
        let outside = cylinder(model, (big, -length, z0), Direction::Y, r, length);
        let inside = cylinder(model, (big, -length, z0), Direction::Y, bore, length);
        let pipe = difference(model, outside, inside).unwrap();
        let t = torus(model, bore);
        let outlet = difference(model, pipe, t).unwrap();
        let outlet = crate::modify::transform(model, outlet, &crate::math::Transform::translation(moved), Default::default()).unwrap();
        let (outer, inner) = (torus(model, r), torus(model, bore));
        let shell = difference(model, outer, inner).unwrap();
        let slab = cuboid(model, (-2.0 * big, -2.0 * big, 0.5 * r), (4.0 * big, 4.0 * big, 2.0 * z0));
        let top = intersection(model, shell, slab).unwrap();
        let cut = cylinder(model, (big, -length, z0), Direction::Y, r, length);
        let volute = difference(model, top, cut).unwrap();
        (outlet, volute)
    }

    /// The outlet and volute of [`outlet_and_volute`]: each operand's section of the cylinder
    /// with the torus was fitted only to within a ten-thousandth, the two passed each other
    /// back and forth, crossing at points found along both, and the faces between them were
    /// slivers (the union invalid and 33.5 short, the intersection failing, the difference
    /// invalid). The operands against quadrature of exact lengths along lines (within its
    /// error, 1e-4 of them); the operations against their volumes. Touching, the torus's seam
    /// turned, a smaller ring, and the outlet moved 1e-3 and 0.1 away along the tube and
    /// across it.
    #[test]
    fn an_outlet_against_the_volute_it_touches_along_a_torus_and_a_cylinder() {
        let housing = (12.8, 3.315, 0.3377, 3.7227, 16.885);
        let still = Vector::new(0.0, 0.0, 0.0);
        let cases = [
            (housing, 0.0, still),
            (housing, 0.4, still),
            ((8.0, 2.5, 0.4, 2.6, 10.0), 0.0, still),
            (housing, 0.0, Vector::new(0.0, -1e-3, 0.0)),
            (housing, 0.0, Vector::new(0.0, -0.1, 0.0)),
            (housing, 0.0, Vector::new(1e-3, 0.0, 0.0)),
            (housing, 0.0, Vector::new(0.1, 0.0, 0.0)),
        ];
        for (sizes, turn, moved) in cases {
            let (big, r, wall, z0, length) = sizes;
            let bore = r - wall;
            let mut model = Model::new();
            let (outlet, volute) = outlet_and_volute(&mut model, sizes, turn, moved);
            let what = format!("ring {big}, turned {turn}, moved {moved:?}");
            let z_min = 0.5 * r;
            // (Pappus: the area of each disc above the cut turned about the axis at `big`.)
            let above = |radius: f64| {
                let d = z0 - z_min;
                PI * radius * radius - (radius * radius * (d / radius).acos() - d * (radius * radius - d * d).sqrt())
            };
            let expected = [
                PI * (r * r - bore * bore) * length - ring_tube_overlap(big, z0, (bore, 0.0), f64::NEG_INFINITY, (r, bore), (-length, 0.0)),
                2.0 * PI * big * (above(r) - above(bore)) - ring_tube_overlap(big, z0, (r, bore), z_min, (r, 0.0), (-length, 0.0)),
            ];
            let mut volumes = [0.0; 2];
            for ((part, expected), volume) in [outlet, volute].into_iter().zip(expected).zip(&mut volumes) {
                let report = validity::check(&model, part);
                assert!(report.is_valid(), "{what}: {report:?}");
                *volume = properties::volume(&model, part).mass;
                assert!((*volume - expected).abs() <= 1e-4 * expected, "{what}: {volume} against {expected}");
            }
            check_against_overlap(&mut model, outlet, volute, (volumes[0], volumes[1], 0.0), &what);
        }
    }
}
