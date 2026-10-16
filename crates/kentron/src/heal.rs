//! Shape healing: repairs for shapes made elsewhere (read from files) that do not quite
//! hold together.
//!
//! [`fix_shape`] makes a shape pass kentron's checks where tolerances or orientations are
//! all that is wrong: tolerances grow to cover the gaps between edges' curves and their
//! curves on faces, and between vertices and the curves and surfaces they lie on (within
//! a limit), and solids turned inside out are turned back. [`convert_to_bspline`],
//! [`split_closed_faces`] and [`split_by_angle`] reshape faces (as OCCT's ShapeCustom and
//! ShapeUpgrade).

use crate::error::Result;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

/// Settings for [`fix_shape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FixOptions {
    /// The precision of the fixes: tolerances are at least this.
    pub precision: f64,
    /// Tolerances do not grow beyond this.
    pub max_tolerance: f64,
}

impl Default for FixOptions {
    fn default() -> Self {
        Self { precision: crate::tolerance::LINEAR, max_tolerance: 1.0 }
    }
}

/// Gives the edges of `shape`'s planar faces the pcurves they lack (their curves projected
/// onto the planes), in place; returns how many were added.
#[doc(alias = "BRepLib::BuildPCurveForEdgesOnPlane")]
#[doc(alias = "bopaddpcs")]
pub fn add_plane_pcurves(model: &mut Model, shape: Shape) -> Result<usize> {
    let mut added = 0;
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface @ crate::geometry::surface::Surface::Plane(_)) = model.surface(face).cloned() else { continue };
        for coedge in model.unique_subshapes(face, NodeKind::CoEdge) {
            let coedge = CoEdgeId::try_from(coedge.node).expect("a coedge");
            if model.pcurve(coedge).is_some() {
                continue;
            }
            let Some((curve, first, last)) = model.edge_curve(Shape::from(model[coedge].edge)) else { continue };
            let Some(pcurve) = crate::geometry::project::analytic_pcurve(&curve, &surface) else { continue };
            model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            added += 1;
        }
    }
    Ok(added)
}

/// Repairs `shape` in place, returning it: the tolerances of its edges cover the distance between
/// their curves and their pcurves on their faces' surfaces (as OCCT's
/// `BRepLib::SameParameter` updates them), those of vertices their distances from the
/// curves and surfaces they bound, all within `[precision, max_tolerance]`; solids
/// enclosing negative volume are reversed (as `ShapeFix_Solid`).
#[doc(alias = "ShapeFix_Shape")]
#[doc(alias = "fixshape")]
pub fn fix_shape(model: &mut Model, shape: Shape, options: &FixOptions) -> Result<Shape> {
    let clamp = |t: f64| t.clamp(options.precision, options.max_tolerance.max(options.precision));
    // Edges: their curves against their pcurves (and their faces' tolerances). The
    // coedges of each edge, kept for the vertices.
    let mut coedges_of_edge: std::collections::HashMap<EdgeId, Vec<CoEdgeId>> = Default::default();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face).cloned() else { continue };
        let face_tolerance = model.tolerance(face);
        for coedge in model.unique_subshapes(face, NodeKind::CoEdge) {
            let coedge = CoEdgeId::try_from(coedge.node).expect("a coedge");
            let edge = model[coedge].edge;
            coedges_of_edge.entry(edge).or_default().push(coedge);
            let gap = pcurve_gap(model, coedge, &surface);
            let needed = clamp(gap * (1.0 + 1e-6) + 1e-12).max(face_tolerance);
            if model.tolerance(Shape::from(edge)) < needed && gap <= options.max_tolerance {
                model.set_edge_tolerance(edge, needed)?;
            }
        }
    }
    // Vertices: the ends of their edges' curves and pcurves.
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        let edge = Shape { orientation: Orientation::Forward, ..edge };
        let edge_id = EdgeId::try_from(edge.node).expect("an edge");
        let edge_tolerance = model.tolerance(edge);
        let (start, end) = model.edge_vertices(edge);
        let ends = model.edge_curve(edge).map(|(c, first, last)| (c.point(first), c.point(last)));
        for (vertex, at) in [(start, ends.map(|e| e.0)), (end, ends.map(|e| e.1))] {
            let p = model.point(vertex);
            let mut needed = edge_tolerance;
            if let Some(q) = at {
                needed = needed.max(p.distance(q) * (1.0 + 1e-6));
            }
            // The pcurves' ends on their surfaces.
            for &coedge in coedges_of_edge.get(&edge_id).into_iter().flatten() {
                let face = model[coedge].face;
                let (Some(surface), Some((pcurve, first, last))) = (model.surface(Shape::from(face)), model.effective_pcurve(coedge)) else { continue };
                let near = |t: f64| {
                    let uv = pcurve.point(t);
                    surface.point(uv.x, uv.y).distance(p)
                };
                needed = needed.max(near(first).min(near(last)) * (1.0 + 1e-6));
            }
            let needed = clamp(needed);
            if model.tolerance(vertex) < needed {
                model.set_vertex_tolerance(VertexId::try_from(vertex.node).expect("a vertex"), needed)?;
            }
        }
    }
    // Solids inside out: their shells turned round.
    for solid in model.unique_subshapes(shape, NodeKind::Solid) {
        let solid = Shape { orientation: Orientation::Forward, ..solid };
        if crate::properties::volume(model, solid).mass < 0.0 {
            for shell in model.children(solid).collect::<Vec<_>>() {
                model.replace_child(solid, shell, shell.reversed())?;
            }
        }
    }
    Ok(shape)
}

/// The largest distance between a coedge's edge curve and its pcurve on `surface` (at
/// equal parameters, sampled).
fn pcurve_gap(model: &Model, coedge: CoEdgeId, surface: &crate::geometry::surface::Surface) -> f64 {
    let Some((pcurve, first, last)) = model.effective_pcurve(coedge) else { return 0.0 };
    let edge = Shape::from(model[coedge].edge);
    let Some((curve, _, _)) = model.edge_curve(edge) else { return 0.0 };
    const SAMPLES: usize = 32;
    (0..=SAMPLES)
        .map(|i| {
            let t = first + (last - first) * i as f64 / SAMPLES as f64;
            let uv = pcurve.point(t);
            surface.point(uv.x, uv.y).distance(curve.point(t))
        })
        .fold(0.0, f64::max)
}

/// Which faces [`convert_to_bspline`] converts, by their surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BSplineConversion {
    /// Surfaces of linear extrusion.
    pub extrusion: bool,
    /// Surfaces of revolution, and cylinders, cones, spheres and tori (as OCCT's DT_ToBspl
    /// converts them to revolutions first).
    pub revolution: bool,
    pub offset: bool,
    pub plane: bool,
}

impl Default for BSplineConversion {
    /// All but planes (OCCT's default `ero`).
    fn default() -> Self {
        Self { extrusion: true, revolution: true, offset: true, plane: false }
    }
}

/// Converts the surfaces of `shape`'s faces that `kinds` names to B-splines over each
/// face's parameter range, in place; the edges keep their curves and get new pcurves.
/// Returns `shape`.
#[doc(alias = "ShapeCustom::ConvertToBSpline")]
#[doc(alias = "DT_ToBspl")]
pub fn convert_to_bspline(model: &mut Model, shape: Shape, kinds: &BSplineConversion) -> Result<Shape> {
    use crate::geometry::surface::{Surface, TrimmedSurface};
    let selected = |s: &Surface| match s {
        Surface::Extrusion(_) => kinds.extrusion,
        Surface::Revolution(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_) => kinds.revolution,
        Surface::Offset(_) => kinds.offset,
        Surface::Plane(_) => kinds.plane,
        _ => false,
    };
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let face = Shape { orientation: Orientation::Forward, ..face };
        let Some(old) = model.surface(face).cloned() else { continue };
        let basis = match &old {
            Surface::Trimmed(t) => t.basis().clone(),
            s => s.clone(),
        };
        if !selected(&basis) {
            continue;
        }
        let (u0, u1, v0, v1) = model.face_uv_bounds(face);
        let Some(bspline) = crate::geometry::convert::surface_to_bspline(&Surface::Trimmed(TrimmedSurface::new(basis, u0, u1, v0, v1))) else {
            continue;
        };
        let new = Surface::BSpline(bspline);
        // The pcurves, projected anew (old ones first: a seam's two uses keep their sides).
        let coedges: Vec<CoEdgeId> = model.unique_subshapes(face, NodeKind::CoEdge).iter().map(|c| CoEdgeId::try_from(c.node).expect("a coedge")).collect();
        let old_pcurves: Vec<Option<(crate::geometry::curve2::Curve2, f64, f64)>> = coedges.iter().map(|&c| model.effective_pcurve(c)).collect();
        let (nu0, nu1, nv0, nv1) = new.bounds();
        // The old parameter rectangle onto the new one (exact on its sides).
        let rescale = |p: crate::math::Point2| {
            let f = |x: f64, a: f64, b: f64, na: f64, nb: f64| if b > a { na + (x - a) / (b - a) * (nb - na) } else { na };
            crate::math::Point2::new(f(p.x, u0, u1, nu0, nu1), f(p.y, v0, v1, nv0, nv1))
        };
        let mut fresh = Vec::with_capacity(coedges.len());
        for (k, &c) in coedges.iter().enumerate() {
            let edge = Shape::from(model[c].edge);
            let tolerance = model.tolerance(edge).max(crate::tolerance::LINEAR);
            let pcurve = match model.edge_curve(edge) {
                Some((curve, first, last)) => {
                    let piece = crate::geometry::curve::Curve::Trimmed(crate::geometry::curve::TrimmedCurve::new(curve, first, last));
                    crate::geometry::project::pcurve(&piece, &new, tolerance.max(1e-7))
                        .map(|p| (crate::geometry::curve2::Curve2::BSpline(p), first, last))
                }
                // A degenerate edge (a pole): the old pcurve's ends, carried over.
                None => old_pcurves[k].as_ref().map(|(p, first, last)| {
                    let (a, b) = (rescale(p.point(*first)), rescale(p.point(*last)));
                    let segment = crate::geometry::curve2::BSplineCurve2::new(vec![a, b], vec![*first, *last], vec![2, 2], 1);
                    (crate::geometry::curve2::Curve2::BSpline(segment), *first, *last)
                }),
            };
            fresh.push(pcurve);
        }
        // Seams: both uses of an edge projected onto one side; the use that was on the low
        // side of the closed direction goes to the low side.
        for i in 0..coedges.len() {
            for j in i + 1..coedges.len() {
                if model[coedges[i]].edge != model[coedges[j]].edge {
                    continue;
                }
                let (Some(oi), Some(oj)) = (&old_pcurves[i], &old_pcurves[j]) else { continue };
                let mid = |p: &(crate::geometry::curve2::Curve2, f64, f64)| p.0.point(0.5 * (p.1 + p.2));
                let (pi, pj) = (mid(oi), mid(oj));
                // The closed direction: where the old uses differ most.
                let across_u = (pi.x - pj.x).abs() >= (pi.y - pj.y).abs();
                let (lo, hi) = if across_u { (nu0, nu1) } else { (nv0, nv1) };
                let low_first = if across_u { pi.x < pj.x } else { pi.y < pj.y };
                for (k, low) in [(i, low_first), (j, !low_first)] {
                    if let Some((p, a, b)) = &fresh[k] {
                        let at = if low { lo } else { hi };
                        fresh[k] = Some((pin(p, across_u, at), *a, *b));
                    }
                }
            }
        }
        if fresh.iter().any(|p| p.is_none()) {
            continue; // a pcurve that could not be found: the face stays as it was
        }
        let tolerance = model.tolerance(face);
        model.set_face_surface(face, new, tolerance)?;
        for (c, p) in coedges.into_iter().zip(fresh) {
            let (pcurve, first, last) = p.expect("checked");
            model.set_pcurve_with_range(c, pcurve, first, last)?;
        }
    }
    fix_shape(model, shape, &FixOptions::default())
}

/// `pcurve` with its u (or v, if `!along_u`) coordinate set to `at` everywhere (a seam's
/// side).
fn pin(pcurve: &crate::geometry::curve2::Curve2, along_u: bool, at: f64) -> crate::geometry::curve2::Curve2 {
    use crate::geometry::curve2::Curve2;
    match pcurve {
        Curve2::BSpline(b) => {
            let poles = b.poles().iter().map(|p| if along_u { crate::math::Point2::new(at, p.y) } else { crate::math::Point2::new(p.x, at) }).collect();
            let mut c = b.clone();
            for (i, p) in Vec::<crate::math::Point2>::into_iter(poles).enumerate() {
                c.set_pole(i, p);
            }
            Curve2::BSpline(c)
        }
        other => other.clone(),
    }
}

/// One edge of a wire as the wire runs it: its curve over `[lo, hi]`, run from `hi` to
/// `lo` if `reversed`.
struct Run {
    coedge: CoEdgeId,
    edge: EdgeId,
    curve: crate::geometry::curve::Curve,
    lo: f64,
    hi: f64,
    reversed: bool,
}

impl Run {
    fn start(&self) -> f64 {
        if self.reversed { self.hi } else { self.lo }
    }

    fn end(&self) -> f64 {
        if self.reversed { self.lo } else { self.hi }
    }

    /// The distance from `p` to the edge, and the parameter of its nearest point.
    fn distance(&self, p: crate::math::Point) -> (f64, f64) {
        let found = crate::geometry::extrema::nearest_on_curve(p, &self.curve, self.lo, self.hi);
        (found.distance, found.parameter)
    }
}

/// What a tail cuts off an edge: nothing, the part next to the tail's tip up to a
/// parameter, or the whole edge.
#[derive(Clone, Copy, PartialEq)]
enum Cut {
    Nothing,
    UpTo(f64),
    Whole,
}

/// Removes the tails of the wires of `shape`'s faces (in place; returns `shape`): where two
/// edges meet at an angle below `max_angle` (radians) and run back within `max_width` of
/// each other, the stretch where they do is cut off and the two edges joined where it
/// ends (OCCT's ShapeFix_Wire::FixTails, as `fixshape -maxtaila -maxtailw` asks for).
#[doc(alias = "ShapeFix_Wire::FixTails")]
#[doc(alias = "ShapeAnalysis_Wire::CheckTail")]
pub fn fix_wire_tails(model: &mut Model, shape: Shape, max_angle: f64, max_width: f64) -> Result<Shape> {
    let max_sine = max_angle.sin().max(0.0);
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let forward = Shape { orientation: Orientation::Forward, ..face };
        for wire in model.children(forward).collect::<Vec<_>>() {
            // Where a tail was just cut off, the next pair there is checked without the angle
            // (as OCCT's).
            let mut joined: Option<VertexId> = None;
            for _ in 0..10_000 {
                let runs = wire_runs(model, wire);
                if runs.len() < 2 {
                    break;
                }
                let n = runs.len();
                let mut fixed = false;
                for i in 0..n {
                    let (a, b) = (&runs[i], &runs[(i + 1) % n]);
                    if a.edge == b.edge {
                        continue;
                    }
                    let at = edge_vertex_at(model, b, true);
                    let check_angle = joined != Some(at);
                    let Some(cuts) = tail(a, b, if check_angle { max_sine } else { -1.0 }, max_width) else { continue };
                    let removed = cuts.iter().filter(|c| **c != Cut::Nothing).count();
                    let whole = cuts.iter().filter(|c| **c == Cut::Whole).count();
                    // (A wire keeps an edge at least.)
                    if n + 2 - whole < 1 + removed {
                        continue;
                    }
                    joined = Some(cut_tail(model, a, b, cuts)?);
                    fixed = true;
                    break;
                }
                if !fixed {
                    break;
                }
            }
        }
    }
    Ok(shape)
}

/// The edges of `wire` in its order, as it runs them.
fn wire_runs(model: &Model, wire: Shape) -> Vec<Run> {
    let mut order: Vec<Shape> = model.children(wire).collect();
    if wire.orientation == Orientation::Reversed {
        order.reverse();
    }
    order
        .into_iter()
        .filter_map(|c| {
            let coedge = CoEdgeId::try_from(c.node).ok()?;
            let edge = model[coedge].edge;
            let (curve, lo, hi) = model.edge_curve(Shape::from(edge))?;
            let reversed = model[coedge].orientation.compose(c.orientation) == Orientation::Reversed;
            Some(Run { coedge, edge, curve, lo, hi, reversed })
        })
        .collect()
}

/// The vertex where `run` starts (or ends) as the wire runs it.
fn edge_vertex_at(model: &Model, run: &Run, start: bool) -> VertexId {
    let (first, last) = model.edge_end_vertices(run.edge);
    if start != run.reversed { first } else { last }
}

/// The tail where `a` ends and `b` begins (OCCT's ShapeAnalysis_Wire::CheckTail): what to
/// cut off each, if they meet at a sine below `max_sine` (unchecked if negative) and run
/// back within `max_width` of each other there.
fn tail(a: &Run, b: &Run, max_sine: f64, max_width: f64) -> Option<[Cut; 2]> {
    const CONFUSION: f64 = 1e-7;
    let (tol2, tol3, tol4) = (max_width + 0.5 * CONFUSION, max_width + CONFUSION, max_width + 1.5 * CONFUSION);
    let runs = [a, b];
    // The tip: where `a` ends and `b` begins.
    let near = [a.end(), b.start()];
    let far = [a.start(), b.end()];
    if a.curve.point(near[0]).distance(b.curve.point(near[1])) > tol2 {
        return None;
    }
    if max_sine >= 0.0 {
        // Both leaving the tip: back along `a`, on along `b`.
        let leaving = |r: &Run, t: f64, forward: bool| {
            let d = r.curve.d1(t).1;
            let along = if r.reversed != forward { d } else { -d };
            along.normalized()
        };
        let (da, db) = (leaving(a, near[0], false), leaving(b, near[1], true));
        if da.dot(db) < 0.0 || da.cross(db).length_squared() > max_sine * max_sine {
            return None;
        }
    }
    // How far each runs within the width of the other.
    let mut cut1 = far;
    let mut cut2 = [0.0; 2];
    let mut dist = [0.0; 2];
    let mut whole = [true, true];
    for k in 0..2 {
        let (me, other) = (runs[k], runs[1 - k]);
        let (d, t) = other.distance(me.curve.point(far[k]));
        dist[k] = d;
        cut2[k] = t;
        if d <= tol2 {
            continue;
        }
        whole[k] = false;
        let (mut inside, mut outside) = (near[k], far[k]);
        for _ in 0..200 {
            let mid = 0.5 * (inside + outside);
            let (d, t) = other.distance(me.curve.point(mid));
            cut2[k] = t;
            if d <= tol2 {
                inside = mid;
            } else {
                outside = mid;
                if d <= tol3 {
                    break;
                }
            }
        }
        cut1[k] = outside;
    }
    // The stretch must stay within the width throughout.
    for k in 0..2 {
        let (me, other) = (runs[k], runs[1 - k]);
        for step in 1..23 {
            let t = near[k] + (cut1[k] - near[k]) * step as f64 / 23.0;
            if other.distance(me.curve.point(t)).0 > tol4 {
                return None;
            }
        }
    }
    // Both wholly within: both go.
    if whole[0] && whole[1] && a.curve.point(far[0]).distance(b.curve.point(far[1])) <= tol3 {
        return Some([Cut::Whole, Cut::Whole]);
    }
    // One edge cut where it leaves the other, the other where that point projects onto it.
    let mut first = 0;
    if whole[0] || whole[1] {
        first = if whole[0] { 0 } else { 1 };
        if dist[1 - first] < dist[first] && whole[1 - first] {
            first = 1 - first;
        }
    }
    let mut params = [0.0; 2];
    params[first] = cut1[first];
    params[1 - first] = cut2[first];
    let mut cuts = [Cut::Nothing; 2];
    for k in 0..2 {
        let r = runs[k];
        let p_conf = 1e-9 * (1.0 + (r.hi - r.lo).abs());
        cuts[k] = if (params[k] - far[k]).abs() <= p_conf {
            Cut::Whole
        } else if (params[k] - near[k]).abs() <= p_conf {
            Cut::Nothing
        } else {
            let (lo, hi) = (params[k].min(near[k]), params[k].max(near[k]));
            let (flo, fhi) = (params[k].min(far[k]), params[k].max(far[k]));
            if crate::geometry::measure::length_between(&r.curve, flo, fhi) <= CONFUSION {
                Cut::Whole
            } else if crate::geometry::measure::length_between(&r.curve, lo, hi) <= CONFUSION {
                Cut::Nothing
            } else {
                Cut::UpTo(params[k])
            }
        };
    }
    (cuts != [Cut::Nothing; 2]).then_some(cuts)
}

/// Cuts the tail where `a` ends and `b` begins as `cuts` say, joining what is left of them
/// at one vertex (returned).
fn cut_tail(model: &mut Model, a: &Run, b: &Run, cuts: [Cut; 2]) -> Result<VertexId> {
    let wire = model[a.coedge].wire;
    // Each edge's last kept vertex towards the tip.
    let mut ends = [edge_vertex_at(model, a, false), edge_vertex_at(model, b, true)];
    for (k, run) in [a, b].into_iter().enumerate() {
        match cuts[k] {
            Cut::Nothing => {}
            Cut::Whole => {
                model.remove_coedge(wire, run.coedge)?;
                // (Its far end: where `a` starts, where `b` ends.)
                ends[k] = edge_vertex_at(model, run, k == 0);
            }
            Cut::UpTo(t) => {
                let point = run.curve.point(t);
                let vertex = model.add_vertex(point, model.tolerance(Shape::from(run.edge)));
                let (low, high) = model.split_edge(run.edge, vertex, t)?;
                // The half next to the tip goes.
                let tip_param = if k == 0 { run.end() } else { run.start() };
                let tip_half = if tip_param >= t { high } else { low };
                let coedge = model
                    .parents(tip_half.into())
                    .into_iter()
                    .filter_map(|n| CoEdgeId::try_from(n).ok())
                    .find(|&c| model[c].wire == wire)
                    .ok_or_else(|| crate::Error::Failed("a split edge lost its wire".into()))?;
                model.remove_coedge(wire, coedge)?;
                ends[k] = vertex;
            }
        }
    }
    let [keep, other] = ends;
    if keep != other {
        let gap = model.point(Shape::from(keep)).distance(model.point(Shape::from(other)));
        let tolerance = model.tolerance(Shape::from(keep)).max(model.tolerance(Shape::from(other))).max(0.5 * gap + 1e-7);
        model.replace_node(other, keep)?;
        model.set_vertex_tolerance(keep, tolerance)?;
    }
    Ok(keep)
}

/// Splits the faces of `shape` closed by a seam in two, midway between the seam's two sides
/// (in place; returns `shape`).
#[doc(alias = "ShapeUpgrade_ShapeDivideClosed")]
#[doc(alias = "DT_ClosedSplit")]
pub fn split_closed_faces(model: &mut Model, shape: Shape) -> Result<Shape> {
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(loops) = face_loops(model, face) else { continue };
        // A seam: an edge the face uses twice; split across the gap between its two uses.
        let mut split: Option<(bool, f64)> = None;
        'search: for lp in &loops {
            for (i, a) in lp.iter().enumerate() {
                for b in &lp[i + 1..] {
                    if a.edge != b.edge {
                        continue;
                    }
                    let (ba, bb) = (pcurve_box(a), pcurve_box(b));
                    // The gap between the two uses' boxes, in u and in v.
                    let gap = |lo1: f64, hi1: f64, lo2: f64, hi2: f64| if lo1 < lo2 { (hi1, lo2) } else { (hi2, lo1) };
                    let (uf, ul) = gap(ba.0, ba.1, bb.0, bb.1);
                    let (vf, vl) = gap(ba.2, ba.3, bb.2, bb.3);
                    split = Some(if ul - uf > vl - vf { (true, 0.5 * (uf + ul)) } else { (false, 0.5 * (vf + vl)) });
                    break 'search;
                }
            }
        }
        if let Some((along_u, at)) = split {
            split_face(model, shape, face, along_u, at)?;
        }
    }
    fix_shape(model, shape, &FixOptions::default())
}

/// Splits the faces of `shape` on surfaces of revolution (and cylinders, cones, spheres and
/// tori) into pieces spanning at most `max_angle` (radians) around the axis, in equal parts
/// (in place; returns `shape`).
#[doc(alias = "ShapeUpgrade_ShapeDivideAngle")]
#[doc(alias = "DT_SplitAngle")]
pub fn split_by_angle(model: &mut Model, shape: Shape, max_angle: f64) -> Result<Shape> {
    use crate::geometry::surface::Surface;
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face) else { continue };
        let basis = match surface {
            Surface::Trimmed(t) => t.basis(),
            Surface::Offset(o) => o.basis(),
            s => s,
        };
        if !matches!(basis, Surface::Revolution(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) {
            continue;
        }
        let (u0, u1, _, _) = model.face_uv_bounds(face);
        let length = u1 - u0;
        let pieces = ((length - crate::tolerance::ANGULAR) / max_angle).floor() as usize + 1;
        // Split off one piece after another (the rest keeps the higher parameters).
        let mut rest = face;
        for k in 1..pieces {
            let at = u0 + length * k as f64 / pieces as f64;
            let parts = split_face(model, shape, rest, true, at)?;
            let Some(&last) = parts.last() else { break };
            rest = last;
        }
    }
    fix_shape(model, shape, &FixOptions::default())
}

/// The parameters strictly inside `[first, last]` where `curve` is not C1 by its
/// representation: a B-spline's knots of multiplicity at least its degree (a trimmed or
/// offset one's basis'); none for other curves.
pub fn curve_breaks(curve: &crate::geometry::curve::Curve, first: f64, last: f64) -> Vec<f64> {
    use crate::geometry::curve::Curve;
    match curve {
        Curve::Trimmed(t) => curve_breaks(t.basis(), first, last),
        Curve::Offset(o) => curve_breaks(o.basis(), first, last),
        Curve::BSpline(b) => knot_breaks(b.knots(), b.multiplicities(), b.degree(), first, last),
        _ => Vec::new(),
    }
}

/// [`curve_breaks`] of a curve in a surface's parameters.
pub fn curve2_breaks(curve: &crate::geometry::curve2::Curve2, first: f64, last: f64) -> Vec<f64> {
    use crate::geometry::curve2::Curve2;
    match curve {
        Curve2::Trimmed(t) => curve2_breaks(t.basis(), first, last),
        Curve2::Offset(o) => curve2_breaks(o.basis(), first, last),
        Curve2::BSpline(b) => knot_breaks(b.knots(), b.multiplicities(), b.degree(), first, last),
        _ => Vec::new(),
    }
}

/// The u and v parameters strictly inside `[u0, u1] x [v0, v1]` where `surface` is not C1
/// across its iso-lines by its representation: a B-spline's knots of multiplicity at least
/// its degree in that direction (a trimmed or offset one's basis'); none for others.
pub fn surface_breaks(surface: &crate::geometry::surface::Surface, u0: f64, u1: f64, v0: f64, v1: f64) -> (Vec<f64>, Vec<f64>) {
    use crate::geometry::surface::Surface;
    match surface {
        Surface::Trimmed(t) => surface_breaks(t.basis(), u0, u1, v0, v1),
        Surface::Offset(o) => surface_breaks(o.basis(), u0, u1, v0, v1),
        Surface::BSpline(b) => (
            knot_breaks(b.u_knots(), b.u_multiplicities(), b.u_degree(), u0, u1),
            knot_breaks(b.v_knots(), b.v_multiplicities(), b.v_degree(), v0, v1),
        ),
        _ => (Vec::new(), Vec::new()),
    }
}

/// The knots strictly inside `[first, last]` (by more than a relative 1e-9) of multiplicity
/// at least `degree`.
fn knot_breaks(knots: &[f64], multiplicities: &[usize], degree: usize, first: f64, last: f64) -> Vec<f64> {
    let (lo, hi) = (first.min(last), first.max(last));
    let margin = 1e-9 * (hi - lo).abs().max(1e-300);
    knots.iter().zip(multiplicities).filter(|&(&k, &m)| m >= degree.max(1) && k > lo + margin && k < hi - margin).map(|(&k, _)| k).collect()
}

/// Divides the faces and edges of `shape` where their surfaces and curves are not C1 (see
/// [`surface_breaks`], [`curve_breaks`]): faces along those iso-lines, edges at those
/// parameters (in place; returns `shape`, or for a face divided a shell of its pieces).
#[doc(alias = "ShapeUpgrade_ShapeDivideContinuity")]
#[doc(alias = "DT_ShapeDivide")]
pub fn split_at_breaks(model: &mut Model, shape: Shape) -> Result<Shape> {
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face).cloned() else { continue };
        let (u0, u1, v0, v1) = model.face_uv_bounds(face);
        let (us, vs) = surface_breaks(&surface, u0, u1, v0, v1);
        if us.is_empty() && vs.is_empty() {
            continue;
        }
        let cuts: Vec<crate::divide::Cut> = us
            .iter()
            .map(|&u| crate::divide::Cut::iso(&surface, true, u, v0, v1))
            .chain(vs.iter().map(|&v| crate::divide::Cut::iso(&surface, false, v, u0, u1)))
            .collect();
        let pieces = crate::divide::divide_face(model, face, &cuts)?;
        // (A face on its own, divided: a shell of its pieces, as OCCT gives.)
        if shape.kind() == NodeKind::Face && pieces.len() > 1 {
            let shell = model.new_empty(NodeKind::Shell);
            for piece in pieces {
                model.add(shell, piece)?;
            }
            return split_at_breaks(model, shell);
        }
    }
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        let Some((curve, first, last)) = model.edge_curve(edge) else { continue };
        let mut breaks = curve_breaks(&curve, first, last);
        // (From the last: the edge split keeps its first part as the edge to split next.)
        breaks.sort_by(|a, b| b.total_cmp(a));
        let mut rest = EdgeId::try_from(edge.node).expect("an edge");
        for t in breaks {
            let v = model.add_vertex(curve.point(t), model.tolerance(edge));
            let (before, _) = model.split_edge(rest, v, t)?;
            rest = before;
        }
    }
    Ok(shape)
}

/// A use of an edge in a face's loop, in the loop's direction of travel.
#[derive(Clone)]
pub(crate) struct LoopUse {
    pub coedge: CoEdgeId,
    pub edge: EdgeId,
    /// The pcurve and its range; the loop runs from `first` to `last` (decreasing if the
    /// use is reversed).
    pub pcurve: crate::geometry::curve2::Curve2,
    pub first: f64,
    pub last: f64,
    pub reversed: bool,
}

impl LoopUse {
    pub(crate) fn at(&self, s: f64) -> crate::math::Point2 {
        self.pcurve.point(self.first + (self.last - self.first) * s)
    }
}

/// The loops of a face as it is used forward (outer first), each a list of edge uses along
/// the loop; `None` if a use has no pcurve or the face is placed.
pub(crate) fn face_loops(model: &Model, face: Shape) -> Option<Vec<Vec<LoopUse>>> {
    if face.placement != crate::math::Transform::IDENTITY {
        return None;
    }
    let forward = Shape { orientation: Orientation::Forward, ..face };
    let mut loops = Vec::new();
    for wire in model.children(forward).collect::<Vec<_>>() {
        let mut order: Vec<Shape> = model.children(wire).collect();
        if wire.orientation == Orientation::Reversed {
            order.reverse();
        }
        let mut uses = Vec::with_capacity(order.len());
        for c in order {
            let coedge = CoEdgeId::try_from(c.node).ok()?;
            // (The coedge as the wire gives it already carries the wire's orientation.)
            let reversed = model[coedge].orientation.compose(c.orientation) == Orientation::Reversed;
            let (pcurve, a, b) = model.effective_pcurve(coedge)?;
            let (first, last) = if reversed { (b, a) } else { (a, b) };
            uses.push(LoopUse { coedge, edge: model[coedge].edge, pcurve, first, last, reversed });
        }
        loops.push(uses);
    }
    Some(loops)
}

/// The parameter box of a use's pcurve (u0, u1, v0, v1), from samples.
fn pcurve_box(u: &LoopUse) -> (f64, f64, f64, f64) {
    let mut b = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..=20 {
        let p = u.at(i as f64 / 20.0);
        b = (b.0.min(p.x), b.1.max(p.x), b.2.min(p.y), b.3.max(p.y));
    }
    b
}

/// Splits `face` of `shape` along its iso-line `u = at` (or `v = at`): when the line crosses
/// the outer loop twice and no hole, two faces replace it in its shells. Returns the faces
/// (the one below `at` first), or the face alone if it is not split.
fn split_face(model: &mut Model, shape: Shape, face: Shape, along_u: bool, at: f64) -> Result<Vec<Shape>> {
    use crate::geometry::curve2::{BSplineCurve2, Curve2};
    use crate::math::Point2;
    let unsplit = Ok(vec![face]);
    let Some(loops) = face_loops(model, face) else { return unsplit };
    let Some(surface) = model.surface(face).cloned() else { return unsplit };
    let coord = |p: Point2| if along_u { p.x } else { p.y };
    let scale = {
        let (u0, u1, v0, v1) = model.face_uv_bounds(face);
        if along_u { u1 - u0 } else { v1 - v0 }.abs().max(1e-300)
    };
    let close = 1e-9 * scale.max(1.0);
    // Crossings of a loop: (use index, fraction along it), where the line passes into the
    // use's interior, or at its end vertex.
    let crossings = |lp: &[LoopUse]| -> Vec<(usize, f64)> {
        let mut out = Vec::new();
        for (i, u) in lp.iter().enumerate() {
            let f = |s: f64| coord(u.at(s)) - at;
            // At its end (the vertex shared with the next use).
            if f(1.0).abs() <= close {
                out.push((i, 1.0));
                continue;
            }
            const N: usize = 32;
            let mut prev = f(0.0);
            for k in 1..=N {
                let s = k as f64 / N as f64;
                let cur = f(s);
                // A sample on the line (not the use's end, taken above).
                if cur.abs() <= close && k < N && prev.abs() > close {
                    out.push((i, s));
                }
                if prev.abs() > close && cur.abs() > close && (prev < 0.0) != (cur < 0.0) {
                    let (mut lo, mut hi, flo) = (s - 1.0 / N as f64, s, prev);
                    for _ in 0..80 {
                        let m = 0.5 * (lo + hi);
                        if (f(m) < 0.0) == (flo < 0.0) { lo = m } else { hi = m }
                    }
                    out.push((i, 0.5 * (lo + hi)));
                }
                prev = cur;
            }
        }
        out
    };
    let outer = crossings(&loops[0]);
    if outer.len() != 2 || loops[1..].iter().any(|lp| !crossings(lp).is_empty()) {
        return unsplit;
    }
    let tolerance = model.tolerance(face).max(crate::tolerance::LINEAR);
    // The crossings made vertices one at a time (two may be on one edge): an edge crossed
    // inside is split there, the loop read again. A degenerate edge (a pole) is not split in
    // the model: its use is cut in two below.
    let mut pieces: Vec<Piece> = Vec::new();
    let mut cross_at: Vec<(usize, VertexId, Point2)> = Vec::new(); // (piece the crossing ends)
    for _ in 0..4 {
        let Some(loops) = face_loops(model, face) else { return unsplit };
        let outer = crossings(&loops[0]);
        if outer.len() != 2 {
            return unsplit;
        }
        let inside_real_edge = outer.iter().find(|&&(i, s)| s < 1.0 && model.edge_curve(Shape::from(loops[0][i].edge)).is_some());
        if let Some(&(i, s)) = inside_real_edge {
            use crate::geometry::ParametricCurve;
            let u = &loops[0][i];
            let edge_shape = Shape::from(u.edge);
            let t = u.first + (u.last - u.first) * s;
            let (curve, _, _) = model.edge_curve(edge_shape).expect("checked");
            let v = model.add_vertex(curve.point(t), model.tolerance(edge_shape));
            model.split_edge(u.edge, v, t)?;
            continue;
        }
        // Every crossing at a vertex or inside a degenerate use: the pieces of the loop.
        for (k, u) in loops[0].iter().enumerate() {
            let (start, end) = model.edge_vertices(Shape::from(u.edge));
            let end_vertex = VertexId::try_from(if u.reversed { start } else { end }.node).expect("a vertex");
            let piece = |first: f64, last: f64, edge: EdgeId| Piece { edge, reversed: u.reversed, pcurve: u.pcurve.clone(), first, last };
            match outer.iter().find(|&&(i, _)| i == k) {
                Some(&(_, s)) if s < 1.0 => {
                    // Inside a degenerate use: two new degenerate edges at its vertex.
                    let t = u.first + (u.last - u.first) * s;
                    let edge_tolerance = model.tolerance(Shape::from(u.edge));
                    let e1 = model.add_edge(end_vertex, end_vertex, edge_tolerance)?;
                    let e2 = model.add_edge(end_vertex, end_vertex, edge_tolerance)?;
                    pieces.push(piece(u.first, t, e1));
                    cross_at.push((pieces.len() - 1, end_vertex, u.pcurve.point(t)));
                    pieces.push(piece(t, u.last, e2));
                }
                Some(_) => {
                    pieces.push(piece(u.first, u.last, u.edge));
                    cross_at.push((pieces.len() - 1, end_vertex, u.pcurve.point(u.last)));
                }
                None => pieces.push(piece(u.first, u.last, u.edge)),
            }
        }
        break;
    }
    let Some(loops) = face_loops(model, face) else { return unsplit };
    if cross_at.len() != 2 {
        return unsplit;
    }
    // The new edge along the iso-line, from the crossing lower along it to the higher.
    let other = |p: Point2| if along_u { p.y } else { p.x };
    let (low, high) = if other(cross_at[0].2) <= other(cross_at[1].2) { (0, 1) } else { (1, 0) };
    let (a, b) = (other(cross_at[low].2), other(cross_at[high].2));
    if !(b - a > close) {
        return unsplit;
    }
    use crate::geometry::ParametricSurface;
    let iso = if along_u { surface.u_iso(at) } else { surface.v_iso(at) };
    let edge = model.add_edge(cross_at[low].1, cross_at[high].1, tolerance)?;
    model.set_edge_curve(edge.into(), iso, tolerance)?;
    model.set_edge_range(edge.into(), a, b)?;
    let point = |s: f64| if along_u { Point2::new(at, s) } else { Point2::new(s, at) };
    let segment = Curve2::BSpline(BSplineCurve2::new(vec![point(a), point(b)], vec![a, b], vec![2, 2], 1));
    // The two loops: from after one crossing round to the other, closed along the new edge
    // (from the crossing the run ends at back to the one it started after).
    let n = pieces.len();
    let (i, j) = (cross_at[0].0, cross_at[1].0);
    let run = |from: usize, to: usize| -> Vec<Piece> {
        let mut out = Vec::new();
        let mut k = (from + 1) % n;
        loop {
            out.push(pieces[k].clone());
            if k == to {
                break;
            }
            k = (k + 1) % n;
        }
        out
    };
    let shells: Vec<Shape> = model
        .unique_subshapes(shape, NodeKind::Shell)
        .into_iter()
        .filter(|s| model.children(Shape { orientation: Orientation::Forward, ..*s }).any(|f| f.node == face.node))
        .collect();
    let face_tolerance = model.tolerance(face);
    let add_loop = |model: &mut Model, face: Shape, pieces: &[Piece], closing: Option<CoEdgeId>| -> Result<()> {
        let mut coedges = Vec::with_capacity(pieces.len() + 1);
        for p in pieces {
            let c = model.add_coedge(p.edge, if p.reversed { Orientation::Reversed } else { Orientation::Forward })?;
            model.set_pcurve_with_range(c, p.pcurve.clone(), p.first.min(p.last), p.first.max(p.last))?;
            coedges.push(c);
        }
        coedges.extend(closing);
        let wire = model.add_wire(&coedges)?;
        model.add(face, wire.into())
    };
    let mut new_faces = Vec::new();
    for (from, to) in [(i, j), (j, i)] {
        let to_is_high = to == cross_at[high].0;
        let closing = model.add_coedge(edge, if to_is_high { Orientation::Reversed } else { Orientation::Forward })?;
        model.set_pcurve_with_range(closing, segment.clone(), a, b)?;
        let new_face = model.new_empty(NodeKind::Face);
        model.set_face_surface(new_face, surface.clone(), face_tolerance)?;
        add_loop(model, new_face, &run(from, to), Some(closing))?;
        new_faces.push(new_face);
    }
    // Holes: with the face on their side.
    let first_low = coord(pieces[(i + 1) % n].pcurve.point(0.5 * (pieces[(i + 1) % n].first + pieces[(i + 1) % n].last))) < at;
    for hole in &loops[1..] {
        let low_side = coord(hole[0].at(0.5)) < at;
        let target = if low_side == first_low { new_faces[0] } else { new_faces[1] };
        let hole_pieces: Vec<Piece> = hole.iter().map(|u| Piece { edge: u.edge, reversed: u.reversed, pcurve: u.pcurve.clone(), first: u.first, last: u.last }).collect();
        add_loop(model, target, &hole_pieces, None)?;
    }
    for shell in shells {
        let forward = Shape { orientation: Orientation::Forward, ..shell };
        let used = model.children(forward).find(|f| f.node == face.node).expect("the face");
        model.remove(forward, used)?;
        for &f in &new_faces {
            model.add(forward, Shape { orientation: used.orientation, ..f })?;
        }
    }
    if !first_low {
        new_faces.reverse();
    }
    Ok(new_faces)
}

/// A piece of a loop for a new face: an edge used along (or against) it over part of its
/// pcurve, from `first` to `last` in the direction of travel.
#[derive(Clone)]
struct Piece {
    edge: EdgeId,
    reversed: bool,
    pcurve: crate::geometry::curve2::Curve2,
    first: f64,
    last: f64,
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_wire_tail_is_cut_off() {
        // OCCT's heal/wire_tails_composed A1: a triangle with a long thin spike to (9, 0);
        // cut where its sides part by 0.7, OCCT's wire is 9.0547 long, of three edges.
        let mut model = Model::new();
        let p = [crate::math::Point::new(-1.0, 0.0, 0.0), crate::math::Point::new(9.0, 0.0, 0.0), crate::math::Point::new(1.0, 1.0, 0.0)];
        let edges: Vec<Shape> = (0..3).map(|i| crate::make::edge_between(&mut model, p[i], p[(i + 1) % 3]).unwrap()).collect();
        let wire = crate::make::wire(&mut model, &edges).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let fixed = fix_wire_tails(&mut model, face, 10f64.to_radians(), 0.7).unwrap();
        assert_eq!(model.unique_subshapes(fixed, NodeKind::Edge).len(), 3);
        assert_eq!(model.unique_subshapes(fixed, NodeKind::Vertex).len(), 3);
        let length = crate::properties::length(&model, fixed).mass;
        assert!((length - 9.054743851305652).abs() < 1e-4, "{length}");
    }
    use crate::geometry::surface::Surface;
    use crate::primitives::{Cuboid, Cylinder, Sphere};

    #[test]
    fn a_face_is_divided_where_its_surface_bends() {
        // A sheet of degree 2 across u, bent at u = 1 (a double knot): two faces, meeting
        // there, and no edge running round the bend.
        use crate::geometry::surface::{BSplineSurface, KnotVector};
        use crate::math::Point;
        let row = |x: f64, z: f64| vec![Point::new(x, 0.0, z), Point::new(x, 1.0, z)];
        let poles = vec![row(0.0, 0.0), row(0.5, 0.2), row(1.0, 1.0), row(1.5, 0.2), row(2.0, 0.0)];
        let surface = Surface::BSpline(BSplineSurface::new(poles, KnotVector::new(vec![0.0, 1.0, 2.0], vec![3, 2, 3], 2), KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1)));
        assert_eq!(surface_breaks(&surface, 0.0, 2.0, 0.0, 1.0), (vec![1.0], vec![]));
        let mut model = Model::new();
        let face = crate::make::face_on_domain(&mut model, surface, 0.0, 2.0, 0.0, 1.0).unwrap();
        let area = crate::properties::area(&model, face).mass;
        let divided = split_at_breaks(&mut model, face).unwrap();
        assert_eq!(divided.kind(), NodeKind::Shell);
        assert_eq!(model.unique_subshapes(divided, NodeKind::Face).len(), 2);
        assert!((crate::properties::area(&model, divided).mass - area).abs() < 1e-6 * area);
        for edge in model.unique_subshapes(divided, NodeKind::Edge) {
            let (curve, first, last) = model.edge_curve(edge).unwrap();
            assert!(curve_breaks(&curve, first, last).is_empty());
        }
        assert!(problems(&model, divided).is_empty(), "{:?}", problems(&model, divided));
    }


    use crate::properties::volume;

    fn problems(model: &Model, s: Shape) -> Vec<String> {
        crate::validity::check(model, s).errors().map(|p| p.description.clone()).collect()
    }

    #[test]
    fn surfaces_become_bsplines_keeping_the_solid() {
        for (name, build) in [
            ("cylinder", Box::new(|m: &mut Model| Cylinder::new(5.0, 10.0).build(m).unwrap()) as Box<dyn Fn(&mut Model) -> Shape>),
            ("sphere", Box::new(|m: &mut Model| Sphere::new(4.0).build(m).unwrap())),
        ] {
            let mut model = Model::new();
            let solid = build(&mut model);
            let before = volume(&model, solid).mass;
            let converted = convert_to_bspline(&mut model, solid, &BSplineConversion::default()).unwrap();
            let surfaces: Vec<Surface> = model.unique_subshapes(converted, NodeKind::Face).iter().filter_map(|&f| model.surface(f).cloned()).collect();
            assert!(surfaces.iter().all(|s| matches!(s, Surface::BSpline(_) | Surface::Plane(_))), "{name}");
            assert!(surfaces.iter().any(|s| matches!(s, Surface::BSpline(_))), "{name}");
            // (Volumes of rational B-spline faces are integrated to about 1e-4.)
            assert!((volume(&model, converted).mass - before).abs() < 1e-3 * before, "{name}: {} vs {before}", volume(&model, converted).mass);
            assert!(problems(&model, converted).is_empty(), "{name}: {:?}", problems(&model, converted));
        }
    }

    #[test]
    fn closed_faces_are_split_in_two() {
        let mut model = Model::new();
        let c = Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
        let before = volume(&model, c).mass;
        let split = split_closed_faces(&mut model, c).unwrap();
        let faces = model.unique_subshapes(split, NodeKind::Face);
        let sides = faces.iter().filter(|&&f| matches!(model.surface(f), Some(Surface::Cylinder(_)))).count();
        assert_eq!(sides, 2);
        assert!((volume(&model, split).mass - before).abs() < 1e-9 * before);
        assert!(problems(&model, split).is_empty(), "{:?}", problems(&model, split));
    }

    #[test]
    fn a_sphere_is_split_into_quarters() {
        let mut model = Model::new();
        let s = Sphere::new(4.0).build(&mut model).unwrap();
        let before = volume(&model, s).mass;
        let split = split_by_angle(&mut model, s, 95f64.to_radians()).unwrap();
        assert_eq!(model.unique_subshapes(split, NodeKind::Face).len(), 4);
        assert!((volume(&model, split).mass - before).abs() < 1e-9 * before);
        assert!(problems(&model, split).is_empty(), "{:?}", problems(&model, split));
    }

    #[test]
    fn an_inside_out_solid_is_turned_back() {
        let mut model = Model::new();
        let b = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let shell = model.children(b).next().unwrap();
        model.replace_child(b, shell, shell.reversed()).unwrap();
        assert!(volume(&model, b).mass < 0.0);
        let fixed = fix_shape(&mut model, b, &FixOptions::default()).unwrap();
        assert!((volume(&model, fixed).mass - 6.0).abs() < 1e-9);
        assert!(problems(&model, fixed).is_empty());
    }
}
