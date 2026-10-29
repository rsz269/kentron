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
mod geom;
mod graph;
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
pub use volume::{enclosed_solids, enclosed_solids_unsplit};
pub(crate) use cells::point_inside;

use build::Kept;
use classify::State;
use graph::Graph;

/// The union of `a` and `b`.
#[doc(alias = "BRepAlgoAPI_Fuse")]
#[doc(alias = "fuse")]
pub fn union(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    // An operand may be a compound of solids that overlap one another: like OCCT's general
    // fuse, every solid is fused with those it meets (boxes overlapping or touching), one
    // after another; groups apart from each other are only gathered.
    if let (Some(sa), Some(sb)) = (only_solids(model, a), only_solids(model, b)) {
        if sa.len() > 1 || sb.len() > 1 {
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

/// The intersection of `a` and `b`.
#[doc(alias = "BRepAlgoAPI_Common")]
#[doc(alias = "common")]
pub fn intersection(model: &mut Model, a: Shape, b: Shape) -> Result<Shape> {
    if let Some(each) = each_solid_apart(model, a, b, Operation::Intersection)? {
        return Ok(each);
    }
    run(model, a, b, Operation::Intersection)
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
    !crate::validity::check(model, shape).problems.iter().any(|p| p.kind == crate::validity::ProblemKind::NotClosed)
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
    let (pieces, arcs_of, regions) = classified(model, &mut graph);
    let mut kept = Vec::new();
    for (fi, region, state) in regions {
        if let Some(flip) = operation.keeps(graph.faces[fi].operand, state) {
            kept.push(Kept { face: fi, region, flip, far: false });
        }
    }
    // Faces far from the other operand lie outside it, whole.
    for (i, f) in graph.far.iter().enumerate() {
        if let Some(flip) = operation.keeps(f.operand, State::Outside) {
            kept.push(Kept { face: i, region: split::Region::whole(), flip, far: true });
        }
    }
    // In the operands' order of faces.
    kept.sort_by_key(|k| if k.far { graph.far[k.face].seq } else { graph.faces[k.face].seq });
    let invalid_operands = |model: &Model| [a, b].iter().any(|&s| !crate::validity::check(model, s).is_valid());
    build::build(model, &graph, &pieces, &arcs_of, &kept, &invalid_operands)
}

/// The operands' faces split along where they meet, and each region with where it lies
/// with respect to the other operand.
#[allow(clippy::type_complexity)]
fn classified(model: &Model, graph: &mut Graph) -> (Vec<split::Piece>, Vec<Vec<split::ArcOnFace>>, Vec<(usize, split::Region, State)>) {
    graph.intersect_edges_with_faces();
    graph.intersect_edges_on_faces();
    graph.share_vertices();
    graph.intersect_faces();
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
    // (Narrower, as twice its area over its perimeter, than a thousandth of its size.)
    let thin = |polygon: &[geom::Uv]| {
        let (lo, hi) = geom::polyline_box(polygon);
        let perimeter: f64 = polygon.windows(2).map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1])).sum();
        2.0 * geom::signed_area(polygon).abs() < 1e-3 * (hi[0] - lo[0]).hypot(hi[1] - lo[1]) * perimeter
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
                    if thin(&r.polygons[0]) {
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
    // unless they disagree or one lies on the other operand (then all are).
    let roots = clean_components(graph, &pieces, &arcs_of, &regions);
    let mut groups: crate::hashing::FastMap<usize, Vec<usize>> = Default::default();
    for (ri, &root) in roots.iter().enumerate() {
        groups.entry(root).or_default().push(ri);
    }
    let classify = |ri: usize| {
        let (fi, region, _) = &regions[ri];
        let face = &graph.faces[*fi];
        let at = |uv: geom::Uv| {
            let p = geom::surface_point(&face.surface, uv);
            let outward = geom::normal_at(&face.surface, uv) * face.outward;
            graph.classify(p, outward, 1 - face.operand, face.drift)
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
        let unanimous = s.iter().all(|ri| known[ri] == known[&s[0]]) && matches!(known[&s[0]], State::Inside | State::Outside);
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
        for ri in g {
            match known.get(ri) {
                Some(State::Inside) => inside += 1,
                Some(State::Outside) => outside += 1,
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
    (pieces, arcs_of, regions)
}

/// `a` and `b` split by one another: the solids inside both, inside `a` only and inside
/// `b` only, sharing the faces, edges and vertices where they meet (each made once).
#[doc(alias = "BOPAlgo_Builder")]
pub(crate) fn partition(model: &mut Model, a: Shape, b: Shape) -> Result<[Shape; 3]> {
    let mut graph = operands(model, a, b, false)?;
    if !graph.only_solids() {
        return Err(Error::InvalidInput("a partition of shapes other than solids".into()));
    }
    let (pieces, arcs_of, regions) = classified(model, &mut graph);
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
/// neighbours; and a large group needs only some of its regions classified).
fn clean_components(graph: &Graph, pieces: &[split::Piece], arcs_of: &[Vec<split::ArcOnFace>], regions: &[(usize, split::Region, State)]) -> Vec<usize> {
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
    let mut parent: Vec<usize> = (0..regions.len()).collect();
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
    (0..regions.len()).map(|ri| find(&mut parent, ri)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    use crate::math::{Direction, Frame, Point, Vector};
    use crate::primitives::{Cone, Cuboid, Cylinder, Sphere};
    use crate::topology::NodeKind;
    use crate::{properties, validity};

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
}
