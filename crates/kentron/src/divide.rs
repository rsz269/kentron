//! Dividing faces: along curves in their parameters (the iso-lines of a grid, edges lying on
//! them), each piece a face of its own.
//!
//! The face's loops and the cuts make a planar arrangement in the face's (u, v): the loops'
//! edges are split where the cuts cross them, the cuts' pieces within the face become edges,
//! and walking round the arrangement with the face on the left gives the pieces' loops
//! (holes no cut reaches going to the piece around them).

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::heal::{face_loops, LoopUse};
use crate::math::{Direction2, Point2, Transformable2, Vector2};
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

/// A curve to divide a face along: in the face's parameters over `[first, last]`, and in
/// space over the same parameters.
#[derive(Clone, Debug)]
pub struct Cut {
    pub pcurve: Curve2,
    pub curve: Curve,
    pub first: f64,
    pub last: f64,
}

impl Cut {
    /// The iso-line `u = at` (with `along_u`) or `v = at` of `surface`, over `[from, to]` of
    /// the other parameter.
    pub fn iso(surface: &crate::geometry::surface::Surface, along_u: bool, at: f64, from: f64, to: f64) -> Cut {
        let (pcurve, curve) = if along_u {
            (Curve2::Line(Line2::new(Point2::new(at, 0.0), Direction2::new(0.0, 1.0))), surface.u_iso(at))
        } else {
            (Curve2::Line(Line2::new(Point2::new(0.0, at), Direction2::new(1.0, 0.0))), surface.v_iso(at))
        };
        Cut { pcurve, curve, first: from, last: to }
    }
}

/// How many pieces to divide a face into ([`divide_by_number`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Parts {
    /// About this many, of equal area, as near squares as the face's extents allow.
    Total(usize),
    /// This many along u by this many along v.
    Grid(usize, usize),
}

/// `face` divided by a grid of its iso-lines into pieces of about equal area (as OCCT's
/// `ShapeUpgrade_ShapeDivideArea` splitting by number): a compound of the pieces, or the
/// face itself where it is not divided.
///
/// The grid spans the face's parameter box (widened by a hundredth across a direction that
/// is not periodic, within the surface's bounds). For a total of `n`, a piece is a square
/// of side `√(area / n)`: across the box as many as the lengths of its middle iso-lines
/// hold (at least one), and at least two across the shorter direction.
#[doc(alias = "ShapeUpgrade_ShapeDivideArea")]
#[doc(alias = "DT_SplitByNumber")]
pub fn divide_by_number(model: &mut Model, face: Shape, parts: Parts) -> Result<Shape> {
    let surface = model.surface(face).cloned().ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
    let area = crate::properties::area(model, face).mass;
    let (mut uf, mut ul, mut vf, mut vl) = model.face_uv_bounds(face);
    let (su0, su1, sv0, sv1) = surface.bounds();
    // (A little wider, so that the pieces' pcurves lie within the grid.)
    if surface.u_period().is_none() {
        let du = 0.01 * (ul - uf);
        uf -= du.min(uf - su0);
        ul += du.min(su1 - ul);
    }
    if surface.v_period().is_none() {
        let dv = 0.01 * (vl - vf);
        vf -= dv.min(vf - sv0);
        vl += dv.min(sv1 - vl);
    }
    let (nu, nv) = match parts {
        Parts::Grid(nu, nv) => (nu.max(1), nv.max(1)),
        Parts::Total(n) => {
            if n <= 1 || area - area / n as f64 <= crate::tolerance::LINEAR {
                return Ok(face);
            }
            // The face's extents: its middle iso-lines' lengths.
            let (um, vm) = (0.5 * (uf + ul), 0.5 * (vf + vl));
            let u_size = crate::geometry::measure::length_between(&surface.v_iso(vm), uf, ul);
            let v_size = crate::geometry::measure::length_between(&surface.u_iso(um), vf, vl);
            let side = (area / n as f64).sqrt();
            let (mut nu, mut nv) = (((u_size / side) as usize).max(1), ((v_size / side) as usize).max(1));
            // The second direction (the shorter, measured by the surface's greatest rates of
            // change) is divided at least in two.
            let (mut su, mut sv) = (0.0f64, 0.0f64);
            for i in 0..=4 {
                for j in 0..=4 {
                    let (_, du, dv) = surface.d1(uf + (ul - uf) * i as f64 / 4.0, vf + (vl - vf) * j as f64 / 4.0);
                    (su, sv) = (su.max(du.length()), sv.max(dv.length()));
                }
            }
            let u_first = (ul - uf) * su > (vl - vf) * sv;
            if u_first && nv == 1 {
                nv = 2;
            } else if !u_first && nu == 1 {
                nu = 2;
            }
            (nu, nv)
        }
    };
    let mut cuts = Vec::with_capacity(nu + nv);
    // (The lines run a little past the box, so that they cross a seam on its side.)
    let (pu, pv) = (0.01 * (ul - uf), 0.01 * (vl - vf));
    for i in 1..nu {
        cuts.push(Cut::iso(&surface, true, uf + (ul - uf) * i as f64 / nu as f64, vf - pv, vl + pv));
    }
    for j in 1..nv {
        cuts.push(Cut::iso(&surface, false, vf + (vl - vf) * j as f64 / nv as f64, uf - pu, ul + pu));
    }
    if cuts.is_empty() {
        return Ok(face);
    }
    let pieces = divide_face(model, face, &cuts)?;
    let compound = model.new_empty(NodeKind::Compound);
    for p in pieces {
        model.add(compound, p)?;
    }
    Ok(compound)
}

/// `shape` with its faces divided along `edges` lying on them (as OCCT's
/// `BRepFeat_SplitShape`): on `face` only, if given, else on whichever faces of `shape` each
/// edge lies on. The faces are replaced in place; returns `shape`.
#[doc(alias = "BRepFeat_SplitShape")]
#[doc(alias = "splitshape")]
pub fn divide_by_edges(model: &mut Model, shape: Shape, face: Option<Shape>, edges: &[Shape]) -> Result<Shape> {
    let faces = match face {
        Some(f) => vec![f],
        None => model.unique_subshapes(shape, NodeKind::Face),
    };
    let mut cuts_of: Vec<(Shape, Vec<Cut>)> = Vec::new();
    for &e in edges {
        let Some((curve, first, last)) = model.edge_curve(e) else { continue };
        let tolerance = model.tolerance(e).max(crate::tolerance::LINEAR);
        for &f in &faces {
            let Some(cut) = edge_on_face(model, f, &curve, first, last, tolerance) else { continue };
            match cuts_of.iter_mut().find(|(g, _)| g.node == f.node) {
                Some((_, cuts)) => cuts.extend(cut),
                None => cuts_of.push((f, cut)),
            }
        }
    }
    for (f, cuts) in cuts_of {
        let pieces = divide_face(model, f, &cuts)?;
        // (The shape a face itself: its pieces together.)
        if f.node == shape.node && pieces.len() > 1 {
            let compound = model.new_empty(NodeKind::Compound);
            for p in pieces {
                model.add(compound, p)?;
            }
            return Ok(compound);
        }
    }
    Ok(shape)
}

/// The cuts on `face` along `curve` over `[first, last]` where it lies on the face's surface
/// (within `tolerance`): its pcurve, and on a periodic surface the same shifted by a
/// period each way (the parts off the face's domain are left out); `None` off the face.
fn edge_on_face(model: &Model, face: Shape, curve: &Curve, first: f64, last: f64, tolerance: f64) -> Option<Vec<Cut>> {
    use crate::geometry::extrema::nearest_on_surface;
    let surface = model.surface(face)?;
    let (u0, u1, v0, v1) = model.face_uv_bounds(face);
    let (du, dv) = (surface.u_period(), surface.v_period());
    const N: usize = 64;
    let ts: Vec<f64> = (0..=N).map(|i| first + (last - first) * i as f64 / N as f64).collect();
    let mut uv: Vec<Point2> = Vec::with_capacity(ts.len());
    for &t in &ts {
        let (u, v, d) = nearest_on_surface(curve.point(t), surface, uv.last().map(|q: &Point2| (q.x, q.y)));
        if d > 10.0 * tolerance {
            return None;
        }
        let mut q = Point2::new(u, v);
        if let Some(&last) = uv.last() {
            let unwrap = |x: f64, previous: f64, period: Option<f64>| period.map_or(x, |p| x + ((previous - x) / p).round() * p);
            q = Point2::new(unwrap(q.x, last.x, du), unwrap(q.y, last.y, dv));
        }
        uv.push(q);
    }
    // (Some of it within the face's parameter box, in some copy.) In closed form where it
    // has one, in the period of the points found.
    let pcurve = match crate::geometry::project::analytic_pcurve(curve, surface) {
        Some(exact) if exact.point(first).distance(uv[0]) <= 1e-6 * (1.0 + uv[0].x.abs() + uv[0].y.abs()) || du.is_some() || dv.is_some() => {
            let at = exact.point(first);
            let shift = |d: f64, period: Option<f64>| period.map_or(0.0, |p| (d / p).round() * p);
            let (su, sv) = (shift(uv[0].x - at.x, du), shift(uv[0].y - at.y, dv));
            let moved = exact.transformed(&crate::math::Transform2::translation(Vector2::new(su, sv)));
            if moved.point(first).distance(uv[0]) <= 1e-6 * (1.0 + uv[0].x.abs() + uv[0].y.abs()) {
                moved
            } else {
                Curve2::BSpline(crate::geometry::fit::interpolate_with_parameters(&uv, &ts).ok()?)
            }
        }
        _ => Curve2::BSpline(crate::geometry::fit::interpolate_with_parameters(&uv, &ts).ok()?),
    };
    let shifts = |period: Option<f64>| period.map_or(vec![0.0], |p| vec![-p, 0.0, p]);
    let mut cuts = Vec::new();
    for su in shifts(du) {
        for sv in shifts(dv) {
            let moved = pcurve.transformed(&crate::math::Transform2::translation(Vector2::new(su, sv)));
            let near = uv.iter().any(|q| (u0..=u1).contains(&(q.x + su)) && (v0..=v1).contains(&(q.y + sv)));
            if near {
                cuts.push(Cut { pcurve: moved, curve: curve.clone(), first, last });
            }
        }
    }
    (!cuts.is_empty()).then_some(cuts)
}

/// A crossing of a cut with a loop use or another cut: where along the cut, and the point.
#[derive(Clone, Copy, Debug)]
struct Crossing {
    at: f64,
    uv: Point2,
}

/// Where two curves in the plane cross over their ranges: the parameters on each.
fn crossings(a: &Curve2, a0: f64, a1: f64, b: &Curve2, b0: f64, b1: f64, eps: f64) -> Vec<(f64, f64)> {
    let samples = |c: &Curve2, t0: f64, t1: f64| -> Vec<(f64, Point2)> {
        let n = if matches!(c, Curve2::Line(_)) { 1 } else { 64 };
        (0..=n).map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            (t, c.point(t))
        }).collect()
    };
    let (pa, pb) = (samples(a, a0, a1), samples(b, b0, b1));
    let mut out: Vec<(f64, f64)> = Vec::new();
    for wa in pa.windows(2) {
        for wb in pb.windows(2) {
            // Segment against segment, then Newton on the curves from there.
            let (p, r) = (wa[0].1, wa[1].1 - wa[0].1);
            let (q, s) = (wb[0].1, wb[1].1 - wb[0].1);
            let den = r.x * s.y - r.y * s.x;
            if den.abs() <= 1e-300 {
                continue;
            }
            let w = q - p;
            let (x, y) = ((w.x * s.y - w.y * s.x) / den, (w.x * r.y - w.y * r.x) / den);
            let slack = 1e-6;
            if !(-slack..=1.0 + slack).contains(&x) || !(-slack..=1.0 + slack).contains(&y) {
                continue;
            }
            let (mut t, mut u) = (wa[0].0 + (wa[1].0 - wa[0].0) * x, wb[0].0 + (wb[1].0 - wb[0].0) * y);
            // (Newton on a(t) - b(u) = 0: [a', -b'] [dt, du]ᵀ = -(a - b).)
            for _ in 0..30 {
                let ((pa, da), (pb, db)) = (a.d1(t), b.d1(u));
                let f = pa - pb;
                let det = -da.x * db.y + db.x * da.y;
                if det.abs() <= 1e-300 {
                    break;
                }
                let dt = (f.x * db.y - db.x * f.y) / det;
                let du = (da.y * f.x - da.x * f.y) / det;
                t += dt;
                u += du;
                if dt.abs() + du.abs() <= 1e-15 * (1.0 + t.abs() + u.abs()) {
                    break;
                }
            }
            // (Within the ranges: where both end at the crossing, Newton may step a little past.)
            let (lo_a, hi_a, lo_b, hi_b) = (a0.min(a1), a0.max(a1), b0.min(b1), b0.max(b1));
            let span_a = (hi_a - lo_a).max(1e-300);
            let (t, u) = (t.clamp(lo_a, hi_a), u.clamp(lo_b, hi_b));
            if a.point(t).distance(b.point(u)) > eps {
                continue;
            }
            if !out.iter().any(|&(x, _)| (x - t).abs() <= 1e-9 * span_a) {
                out.push((t, u));
            }
        }
    }
    out
}

/// One way along a piece of an edge in the arrangement: the edge, whether along its own
/// direction, its pcurve from `from` to `to` (in the direction of travel), and its end nodes.
#[derive(Clone)]
struct Arc {
    edge: EdgeId,
    along: bool,
    pcurve: Curve2,
    from: f64,
    to: f64,
    start: usize,
    end: usize,
}

impl Arc {
    fn at(&self, s: f64) -> Point2 {
        self.pcurve.point(self.from + (self.to - self.from) * s)
    }
}

/// Whether `edge` is one of the loops' edges.
fn boundary_of(loops: &[Vec<LoopUse>], edge: EdgeId) -> bool {
    loops.iter().flatten().any(|u| u.edge == edge)
}

/// Even-odd: whether `q` lies within the closed polylines.
fn within(polygons: &[Vec<Point2>], q: Point2) -> bool {
    let mut inside = false;
    for poly in polygons {
        for i in 0..poly.len() {
            let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
            if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                inside = !inside;
            }
        }
    }
    inside
}

fn polygon_of(arcs: &[&Arc]) -> Vec<Point2> {
    arcs.iter().flat_map(|a| (0..16).map(move |i| a.at(i as f64 / 16.0))).collect()
}

fn signed_area(points: &[Point2]) -> f64 {
    let n = points.len();
    0.5 * (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum::<f64>()
}

/// `face` divided along `cuts` (in its parameters; parts of them off the face are left
/// out): the pieces, which replace it in the shells holding it (its edges the cuts cross
/// split there, in its neighbours too). The face alone if no cut divides it.
pub fn divide_face(model: &mut Model, face: Shape, cuts: &[Cut]) -> Result<Vec<Shape>> {
    let failed = |what: &str| Error::Failed(format!("dividing a face: {what}"));
    let surface = model.surface(face).cloned().ok_or_else(|| failed("a face without a surface"))?;
    let loops = face_loops(model, face).ok_or_else(|| failed("a face without pcurves"))?;
    let (u0, u1, v0, v1) = model.face_uv_bounds(face);
    let scale = (u1 - u0).hypot(v1 - v0).max(1e-300);
    let eps = 1e-9 * scale.max(1.0);
    let tolerance = model.tolerance(face).max(crate::tolerance::LINEAR);
    let loop_polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect();

    // Where the cuts cross the loops: on each use, at its fraction along it; on each cut.
    let mut on_uses: Vec<Vec<Vec<(f64, Point2)>>> = loops.iter().map(|l| vec![Vec::new(); l.len()]).collect();
    let mut on_cuts: Vec<Vec<Crossing>> = vec![Vec::new(); cuts.len()];
    // (Near a use's ends: at its vertex. Only where the loop passes from one side of the cut
    // to the other: a cut touching a loop does not divide the face there.)
    let snap = (1e-6 * scale).max(eps);
    for (li, l) in loops.iter().enumerate() {
        for (ui, u) in l.iter().enumerate() {
            let next = &l[(ui + 1) % l.len()];
            for (ci, c) in cuts.iter().enumerate() {
                for (t, tau) in crossings(&u.pcurve, u.first, u.last, &c.pcurve, c.first, c.last, eps.max(1e-7 * scale)) {
                    let s = (t - u.first) / (u.last - u.first);
                    let uv = u.pcurve.point(t);
                    // (At the use's start: that is the previous use's end, taken there; a closed
                    // use's start and end are one, taken at its end.)
                    let at_start = s <= 1e-9 || uv.distance(u.at(0.0)) <= snap;
                    let at_end = s >= 1.0 - 1e-9 || uv.distance(u.at(1.0)) <= snap;
                    if at_start && (!at_end || s < 0.5) {
                        continue;
                    }
                    let (pc, tangent) = c.pcurve.d1(tau);
                    let side = |q: Point2| {
                        let w = q - pc;
                        tangent.x * w.y - tangent.y * w.x
                    };
                    let step = 1e-3;
                    let (before, after) = if at_end { (u.at(1.0 - step), next.at(step)) } else { (u.at((s - step).max(0.0)), u.at((s + step).min(1.0))) };
                    if (side(before) > 0.0) == (side(after) > 0.0) {
                        continue;
                    }
                    let (s, uv) = if at_end { (1.0, u.at(1.0)) } else { (s, uv) };
                    on_uses[li][ui].push((s, uv));
                    on_cuts[ci].push(Crossing { at: tau, uv });
                }
            }
        }
    }
    // Where cuts cross each other (in the face or not: the pieces off it are dropped).
    let mut cut_crossings: Vec<(usize, f64, usize, f64, Point2)> = Vec::new();
    for i in 0..cuts.len() {
        for j in i + 1..cuts.len() {
            for (t, u) in crossings(&cuts[i].pcurve, cuts[i].first, cuts[i].last, &cuts[j].pcurve, cuts[j].first, cuts[j].last, eps.max(1e-7 * scale)) {
                cut_crossings.push((i, t, j, u, cuts[i].pcurve.point(t)));
            }
        }
    }
    // (Cuts dividing nothing: none crossing the loops or another within the face, none closed
    // within it.)
    let closed_inside = |c: &Cut| c.pcurve.point(c.first).distance(c.pcurve.point(c.last)) <= (1e-6 * scale).max(eps) && within(&loop_polygons, c.pcurve.point(0.5 * (c.first + c.last)));
    let crossing_inside = cut_crossings.iter().any(|x| within(&loop_polygons, x.4));
    if on_cuts.iter().all(|c| c.is_empty()) && !cuts.iter().any(closed_inside) && !crossing_inside {
        return Ok(vec![face]);
    }

    // The loops' edges split where crossed inside (at one vertex where two uses of a seam
    // are crossed at one point of it); degenerate ones are cut below.
    // (Each with the face's point there, which the cuts' edges end at.)
    let mut splits: Vec<(EdgeId, f64, Point2)> = Vec::new();
    for (li, l) in loops.iter().enumerate() {
        for (ui, u) in l.iter().enumerate() {
            if model.edge_curve(Shape::from(u.edge)).is_none() {
                continue;
            }
            for &(s, uv) in &on_uses[li][ui] {
                if s < 1.0 {
                    let t = u.first + (u.last - u.first) * s;
                    if !splits.iter().any(|&(e, x, _)| e == u.edge && (x - t).abs() <= 1e-9 * (u.last - u.first).abs()) {
                        splits.push((u.edge, t, uv));
                    }
                }
            }
        }
    }
    splits.sort_by(|a, b| (a.0.index(), a.1).partial_cmp(&(b.0.index(), b.1)).expect("parameters"));
    let mut current: Vec<(EdgeId, EdgeId)> = Vec::new(); // (original, piece holding the rest)
    for (edge, t, uv) in splits {
        let piece = current.iter().find(|c| c.0 == edge).map_or(edge, |c| c.1);
        let (curve, ..) = model.edge_curve(Shape::from(piece)).ok_or_else(|| failed("an edge without a curve"))?;
        // (Its tolerance covering the face's point there, where the cuts end.)
        let p = curve.point(t);
        let v = model.add_vertex(p, model.tolerance(Shape::from(piece)).max(2.0 * p.distance(surface.point(uv.x, uv.y))));
        let (_, rest) = model.split_edge(piece, v, t)?;
        match current.iter_mut().find(|c| c.0 == edge) {
            Some(c) => c.1 = rest,
            None => current.push((edge, rest)),
        }
    }
    let loops = face_loops(model, face).ok_or_else(|| failed("a face without pcurves"))?;

    // The nodes: a vertex at a point of the face's (u, v).
    let mut nodes: Vec<(VertexId, Point2)> = Vec::new();
    let node_eps = (1e-6 * scale).max(eps);
    let node = |nodes: &mut Vec<(VertexId, Point2)>, v: VertexId, uv: Point2| -> usize {
        match nodes.iter().position(|&(w, q)| w == v && q.distance(uv) <= node_eps) {
            Some(i) => i,
            None => {
                nodes.push((v, uv));
                nodes.len() - 1
            }
        }
    };
    let ends = |model: &Model, u: &LoopUse| -> (VertexId, VertexId) {
        let (a, b) = model.edge_end_vertices(u.edge);
        if u.reversed { (b, a) } else { (a, b) }
    };
    let mut arcs: Vec<Arc> = Vec::new();
    for l in &loops {
        for u in l {
            let (a, b) = ends(model, u);
            let (sa, sb) = (node(&mut nodes, a, u.at(0.0)), node(&mut nodes, b, u.at(1.0)));
            if model.edge_curve(Shape::from(u.edge)).is_some() {
                arcs.push(Arc { edge: u.edge, along: !u.reversed, pcurve: u.pcurve.clone(), from: u.first, to: u.last, start: sa, end: sb });
                continue;
            }
            // A degenerate edge (a pole) crossed inside: cut in pieces, a new degenerate edge
            // at its vertex each.
            let mut stops: Vec<f64> = loops_crossings_on(&on_uses, &loops, u).into_iter().filter(|&s| s > 0.0 && s < 1.0).collect();
            stops.sort_by(f64::total_cmp);
            if stops.is_empty() {
                arcs.push(Arc { edge: u.edge, along: !u.reversed, pcurve: u.pcurve.clone(), from: u.first, to: u.last, start: sa, end: sb });
                continue;
            }
            let mut previous = (0.0, sa);
            for s in stops.into_iter().chain([1.0]) {
                let to_node = if s >= 1.0 { sb } else { node(&mut nodes, a, u.at(s)) };
                let (t0, t1) = (u.first + (u.last - u.first) * previous.0, u.first + (u.last - u.first) * s);
                let e = model.add_edge(a, a, model.tolerance(Shape::from(u.edge)))?;
                // (Its pcurve runs along the travel: the edge forward.)
                arcs.push(Arc { edge: e, along: true, pcurve: u.pcurve.clone(), from: t0, to: t1, start: previous.1, end: to_node });
                previous = (s, to_node);
            }
        }
    }

    // The cuts' pieces within the face, between their crossings with the loops and each
    // other.
    let mut cut_vertices: Vec<Option<VertexId>> = vec![None; cut_crossings.len()];
    for (ci, c) in cuts.iter().enumerate() {
        let mut stops: Vec<(f64, usize)> = Vec::new();
        for x in &on_cuts[ci] {
            // (The node already made where the loop is crossed.)
            let n = nodes.iter().position(|&(_, q)| q.distance(x.uv) <= node_eps.max(1e-6 * scale)).ok_or_else(|| failed("a crossing off the loops' nodes"))?;
            stops.push((x.at, n));
        }
        for (k, &(i, ti, j, tj, uv)) in cut_crossings.iter().enumerate() {
            let at = if i == ci { ti } else if j == ci { tj } else { continue };
            // (Where the cuts meet on the loops, at the node there.)
            if let Some(n) = nodes.iter().position(|&(_, q)| q.distance(uv) <= node_eps.max(1e-6 * scale)) {
                stops.push((at, n));
                continue;
            }
            if !within(&loop_polygons, uv) {
                continue;
            }
            let v = match cut_vertices[k] {
                Some(v) => v,
                None => {
                    let v = model.add_vertex(surface.point(uv.x, uv.y), tolerance);
                    cut_vertices[k] = Some(v);
                    v
                }
            };
            stops.push((at, node(&mut nodes, v, uv)));
        }
        stops.sort_by(|a, b| a.0.total_cmp(&b.0));
        stops.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-12 * (c.last - c.first).abs() && a.1 == b.1);
        // A closed cut crossing nothing, within the face: one closed edge.
        let (start, end) = (c.pcurve.point(c.first), c.pcurve.point(c.last));
        if stops.is_empty() && start.distance(end) <= node_eps && within(&loop_polygons, c.pcurve.point(0.5 * (c.first + c.last))) {
            let v = model.add_vertex(c.curve.point(c.first), tolerance);
            let n = node(&mut nodes, v, start);
            stops = vec![(c.first, n), (c.last, n)];
        }
        for w in stops.windows(2) {
            let ((ta, na), (tb, nb)) = (w[0], w[1]);
            if tb - ta <= 1e-12 * (c.last - c.first).abs() || (na == nb && stops.len() > 2) {
                continue;
            }
            if !within(&loop_polygons, c.pcurve.point(0.5 * (ta + tb))) {
                continue;
            }
            let e = model.add_edge(nodes[na].0, nodes[nb].0, tolerance)?;
            model.set_edge_curve(e.into(), c.curve.clone(), tolerance)?;
            model.set_edge_range(e.into(), ta, tb)?;
            arcs.push(Arc { edge: e, along: true, pcurve: c.pcurve.clone(), from: ta, to: tb, start: na, end: nb });
            arcs.push(Arc { edge: e, along: false, pcurve: c.pcurve.clone(), from: tb, to: ta, start: nb, end: na });
        }
    }

    // Pieces of cuts lying along another piece (cuts on one curve twice) are taken once.
    let mut kept: Vec<usize> = Vec::new();
    let mut dropped: Vec<EdgeId> = Vec::new();
    for (i, a) in arcs.iter().enumerate() {
        let doubled = !boundary_of(&loops, a.edge)
            && a.along
            && kept.iter().any(|&k| {
                let b = &arcs[k];
                ((b.start, b.end) == (a.start, a.end) || (b.start, b.end) == (a.end, a.start)) && b.at(0.5).distance(a.at(0.5)) <= node_eps.max(1e-6 * scale)
            });
        if doubled {
            dropped.push(a.edge);
        } else {
            kept.push(i);
        }
    }
    arcs.retain(|a| !dropped.contains(&a.edge));

    // Pieces of cuts ending inside the face, at a node no other edge reaches, divide nothing:
    // left out (repeatedly, as leaving one out may leave another so).
    let boundary: std::collections::HashSet<EdgeId> = loops.iter().flatten().map(|u| u.edge).collect();
    loop {
        let mut degree = vec![0usize; nodes.len()];
        for a in arcs.iter().filter(|a| a.along || boundary.contains(&a.edge)) {
            degree[a.start] += 1;
            degree[a.end] += 1;
        }
        let before = arcs.len();
        arcs.retain(|a| boundary.contains(&a.edge) || (degree[a.start] > 1 && degree[a.end] > 1));
        if arcs.len() == before {
            break;
        }
    }

    // Round the pieces: from each arc on, at each node the next arc turning the most to the
    // left (the smallest clockwise angle from the way back).
    let leave = |a: &Arc| -> Vector2 {
        let d = a.at(1e-4) - a.at(0.0);
        let len = d.length();
        if len > 0.0 { d * (1.0 / len) } else { d }
    };
    let back = |a: &Arc| -> Vector2 {
        let d = a.at(1.0 - 1e-4) - a.at(1.0);
        let len = d.length();
        if len > 0.0 { d * (1.0 / len) } else { d }
    };
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for (i, a) in arcs.iter().enumerate() {
        outgoing[a.start].push(i);
    }
    let mut used = vec![false; arcs.len()];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for first in 0..arcs.len() {
        if used[first] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut a = first;
        loop {
            used[a] = true;
            cycle.push(a);
            let r = back(&arcs[a]);
            let base = r.y.atan2(r.x);
            // (Not straight back along the same edge unless nothing else leaves.)
            let next = outgoing[arcs[a].end]
                .iter()
                .copied()
                .filter(|&b| !(arcs[b].edge == arcs[a].edge && arcs[b].along != arcs[a].along) || outgoing[arcs[a].end].len() == 1)
                .min_by(|&x, &y| {
                    let cw = |b: usize| {
                        let d = leave(&arcs[b]);
                        let angle = base - d.y.atan2(d.x);
                        let angle = angle.rem_euclid(std::f64::consts::TAU);
                        if angle <= 1e-12 { std::f64::consts::TAU } else { angle }
                    };
                    cw(x).total_cmp(&cw(y))
                })
                .ok_or_else(|| failed("a node with no way on"))?;
            if next == first {
                break;
            }
            if used[next] || cycle.len() > arcs.len() {
                return Err(failed("a walk round a piece not closing"));
            }
            a = next;
        }
        cycles.push(cycle);
    }

    // Outer loops (counter-clockwise) make the pieces; the others go in the piece around them.
    let polygons: Vec<Vec<Point2>> = cycles.iter().map(|c| polygon_of(&c.iter().map(|&i| &arcs[i]).collect::<Vec<_>>())).collect();
    let areas: Vec<f64> = polygons.iter().map(|p| signed_area(p)).collect();
    let outer: Vec<usize> = (0..cycles.len()).filter(|&i| areas[i] > 0.0).collect();
    if outer.is_empty() {
        return Err(failed("no piece"));
    }
    let mut groups: Vec<Vec<usize>> = outer.iter().map(|&i| vec![i]).collect();
    for i in (0..cycles.len()).filter(|&i| areas[i] <= 0.0) {
        // (Just to the left of it, where the piece around it lies: not on the loop, which
        // may also bound the piece it holes out of the face.)
        let a = &arcs[cycles[i][0]];
        let ahead = a.at(0.5 + 1e-4) - a.at(0.5 - 1e-4);
        let ahead = ahead * (1.0 / ahead.length().max(1e-300));
        let q = a.at(0.5) + Vector2::new(-ahead.y, ahead.x) * (1e-6 * scale).max(eps);
        let around = outer.iter().enumerate().filter(|&(_, &o)| within(&polygons[o..=o], q)).min_by(|a, b| areas[*a.1].total_cmp(&areas[*b.1])).map(|(k, _)| k);
        groups[around.unwrap_or(0)].push(i);
    }
    let face_tolerance = model.tolerance(face);
    let mut pieces = Vec::with_capacity(groups.len());
    for group in groups {
        let f = model.new_empty(NodeKind::Face);
        model.set_face_surface(f, surface.clone(), face_tolerance)?;
        for c in group {
            let mut coedges = Vec::with_capacity(cycles[c].len());
            for &i in &cycles[c] {
                let a = &arcs[i];
                let co = model.add_coedge(a.edge, if a.along { Orientation::Forward } else { Orientation::Reversed })?;
                model.set_pcurve_with_range(co, a.pcurve.clone(), a.from.min(a.to), a.from.max(a.to))?;
                coedges.push(co);
            }
            let wire = model.add_wire(&coedges)?;
            model.add(f, wire.into())?;
        }
        pieces.push(f);
    }
    // In the shells (and compounds) holding the face, the pieces instead.
    let parents: Vec<crate::topology::NodeId> = model.parents(face.node).into_iter().filter(|n| matches!(n.kind, NodeKind::Shell | NodeKind::Compound)).collect();
    for parent in parents {
        let holder = Shape::new(parent);
        let uses: Vec<Shape> = model.children_local(holder).filter(|f| f.node == face.node).collect();
        for used in uses {
            model.remove(holder, used)?;
            for &f in &pieces {
                model.add(holder, Shape { orientation: used.orientation, placement: used.placement, ..f })?;
            }
        }
    }
    Ok(pieces.into_iter().map(|f| Shape { orientation: face.orientation, ..f }).collect())
}

/// The fractions along `u` (a use of the loops read before the edges were split, matched by
/// its edge and pcurve range) where cuts cross it.
fn loops_crossings_on(on_uses: &[Vec<Vec<(f64, Point2)>>], _loops: &[Vec<LoopUse>], u: &LoopUse) -> Vec<f64> {
    // (A degenerate use is not split in the model: its crossings are where its points lie.)
    let mut out = Vec::new();
    for l in on_uses {
        for crossings in l {
            for &(_, uv) in crossings {
                // On this use: the point at some fraction of it.
                let n = 64;
                if let Some(s) = (0..=n).map(|i| i as f64 / n as f64).min_by(|&a, &b| u.at(a).distance(uv).total_cmp(&u.at(b).distance(uv))) {
                    let span = u.at(0.0).distance(u.at(1.0)).max(1e-300);
                    if u.at(s).distance(uv) <= 0.5 * span / n as f64 {
                        // Refined on the straight-enough piece.
                        let (a, b) = (u.at(0.0), u.at(1.0));
                        let ab = b - a;
                        let f = ((uv - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
                        if !out.iter().any(|&x: &f64| (x - f).abs() <= 1e-9) {
                            out.push(f);
                        }
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cylinder, Sphere};

    fn faces(model: &Model, s: Shape) -> Vec<Shape> {
        model.unique_subshapes(s, NodeKind::Face)
    }

    /// OCCT's heal/split_number A1 and A2: a hemisphere's face in 10 and in 102 pieces (5 by
    /// 2 and 17 by 6), each valid, together the face.
    #[test]
    fn a_hemisphere_in_squares() {
        for (n, expected) in [(10, 10), (100, 102)] {
            let mut model = Model::new();
            let solid = Sphere::new(10.0).with_latitudes(0.0, std::f64::consts::FRAC_PI_2).build(&mut model).unwrap();
            let face = faces(&model, solid).into_iter().find(|&f| matches!(model.surface(f), Some(crate::geometry::surface::Surface::Sphere(_)))).unwrap();
            let area = crate::properties::area(&model, face).mass;
            let divided = divide_by_number(&mut model, face, Parts::Total(n)).unwrap();
            let pieces = faces(&model, divided);
            assert_eq!(pieces.len(), expected);
            let total: f64 = pieces.iter().map(|&f| crate::properties::area(&model, f).mass).sum();
            assert!((total - area).abs() <= 1e-6 * area, "{total} {area}");
            for f in pieces {
                assert!(crate::validity::check(&model, f).is_valid(), "{:?}", crate::validity::check(&model, f));
            }
            // The solid holding the face: the pieces instead, still closed.
            assert!(crate::validity::check(&model, solid).is_valid(), "{:?}", crate::validity::check(&model, solid));
        }
    }

    /// A cut touching a hole of a face (a cylinder's side less a hole reaching to u = 2π/3
    /// exactly, as in OCCT's heal/split_number D1) does not divide it there.
    #[test]
    fn a_cut_touching_a_hole_passes_by() {
        let mut model = Model::new();
        let cylinder = Cylinder::new(10.0, 50.0).build(&mut model).unwrap();
        let frame = crate::math::Frame::from_z(crate::math::Point::new(0.0, -20.0, 30.0), crate::math::Direction::Y);
        let drill = Cylinder::new(5.0, 40.0).with_frame(frame).build(&mut model).unwrap();
        let cut = crate::boolean::difference(&mut model, cylinder, drill).unwrap();
        let side = faces(&model, cut).into_iter().find(|&f| matches!(model.surface(f), Some(crate::geometry::surface::Surface::Cylinder(_)))).unwrap();
        let divided = divide_by_number(&mut model, side, Parts::Grid(3, 2)).unwrap();
        assert_eq!(faces(&model, divided).len(), 6);
        assert!(crate::validity::check(&model, cut).is_valid(), "{:?}", crate::validity::check(&model, cut));
    }
}
