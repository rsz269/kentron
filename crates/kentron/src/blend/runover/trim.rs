//! A run-over fillet around a closed edge, built directly: the fillet's faces between the
//! sections where the contacts switch, and the solid's faces the contacts cross trimmed
//! along them, without boolean operations.
//!
//! 1. The sections are taken from a switch on (from the edge's own start where nothing
//!    switches), and between consecutive switches the ball's arcs make a surface (exact
//!    arcs across, interpolated along), bounded by the arcs at its ends and the contact
//!    curves along its sides.
//! 2. A contact curve lies on a face (where the ball rolls on it, a new edge, the side of
//!    the fillet surface) or on a face's edge (where the ball rolls over it, a piece of that
//!    edge). The edges the contacts switch on are split there, and the seams of periodic
//!    faces where the contact curves cross them.
//! 3. Each face the contact curves run on is rebuilt by walking around its kept part: along
//!    a run of contact curves on it (the part away from the edge on its left), then along its
//!    boundary from where the run ends to where the next one starts, and so on round.
//!    Boundary loops no run meets are kept or dropped whole, as is a face no run crosses
//!    among those the ball reaches: dropped where they lie between the edge and the
//!    contacts.

use crate::hashing::FastMap as HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::extrema::{nearest_on_curve, nearest_on_surface};
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction2, Point, Point2, Transform2, Transformable2, Vector, Vector2};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{Section, Side, Support};

const TOL: f64 = 1e-7;

/// An edge as a loop runs along it (with its curve or against it), with its pcurve on the
/// face over the edge's range.
#[derive(Clone)]
struct Use {
    edge: EdgeId,
    forward: bool,
    pcurve: Curve2,
    first: f64,
    last: f64,
    /// Whether it is of the face's own boundary (its pcurve placed in the face's domain).
    own: bool,
}

impl Use {
    fn start(&self, model: &Model) -> VertexId {
        let (a, b) = model.edge_end_vertices(self.edge);
        if self.forward { a } else { b }
    }

    fn end(&self, model: &Model) -> VertexId {
        let (a, b) = model.edge_end_vertices(self.edge);
        if self.forward { b } else { a }
    }

    /// The point of the pcurve at fraction `s` of the way along the loop.
    fn uv(&self, s: f64) -> Point2 {
        let t = if self.forward { self.first + (self.last - self.first) * s } else { self.last - (self.last - self.first) * s };
        self.pcurve.point(t)
    }

    fn shifted(&self, by: Vector2) -> Use {
        if by.x == 0.0 && by.y == 0.0 {
            return self.clone();
        }
        Use { pcurve: self.pcurve.transformed(&Transform2::translation(by)), ..self.clone() }
    }
}

/// The fillet surface between two switches, with its parameter ranges (`u` across the
/// fillet from the first contact to the second, `v` along it).
struct Patch {
    surface: Surface,
    u: (f64, f64),
    v: (f64, f64),
}

/// What a contact curve along a stretch is: new edges on a face (in order along the
/// fillet, run forward), or a piece of an edge the ball rolls over (with whether the edge
/// runs along the fillet).
#[derive(Clone)]
enum Contact {
    OnFace(usize, Vec<EdgeId>),
    OnEdge(EdgeId, bool),
}

/// Which side of a run of contact curves on a face is kept: away from the fillet's edge
/// (a run of side `k`'s contacts), or toward or away from the ball's centre (where the
/// fillet crosses a face at an end of its edge).
#[derive(Clone, Copy)]
enum Kept {
    AwayFromEdge(usize),
    Ball { center: Point, toward: bool },
}

/// Where the fillet's faces end, between stretches or at an end of the edge: the vertices
/// on each side, the edge across, and on each side the `v` where it meets the contact curve.
struct Boundary {
    vertices: [VertexId; 2],
    edge: EdgeId,
    v: [f64; 2],
    /// Across the fillet at an end of the edge (its pcurve on the end's patch, as `u` runs);
    /// otherwise an arc between patches, along `v = const` on both.
    across: Option<Curve2>,
}

/// `shape` (which this edits: a copy) with its `edge` rounded by the ball of `radius` in
/// `sections` (around a closed edge, the last the first again), touching what `supports`
/// say, switching at `breaks`, on the faces of `sides`; along an open edge, `tangents` are
/// its directions at its ends.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(model: &mut Model, shape: Shape, edge: Shape, sections: &[Section], supports: &[[Support; 2]], breaks: &[usize], sides: &[Side; 2], radius: f64, convex: bool, closed: bool, tangents: [Vector; 2]) -> Result<Shape> {
    let failed = |what: &str| Error::Failed(format!("a run-over fillet built directly: {what}"));
    let size = 1.0 + radius;
    let near = 1e-6 * size;
    let m = sections.len() - 1;
    let (secs, supp, mut starts): (Vec<Section>, Vec<[Support; 2]>, Vec<usize>) = if closed {
        // From a switch on, round to it again (the edge's parameter carried on past its end,
        // increasing).
        let s0 = breaks.first().copied().unwrap_or(0) % m;
        let period = sections[m].t - sections[0].t;
        let secs = (0..=m).map(|j| Section { t: sections[(s0 + j) % m].t + if s0 + j >= m { period } else { 0.0 }, ..sections[(s0 + j) % m] }).collect();
        let supp = (0..=m).map(|j| supports[(s0 + j) % m]).collect();
        (secs, supp, breaks.iter().map(|&b| (b + m - s0) % m).chain([0]).collect())
    } else {
        // The end sections carried on past the ends along the edge (brought into the plane
        // square to it there), for the fillet to reach across the faces at its ends.
        let carry = |s: &Section, direction: Vector, step: f64| {
            let by = direction * (2.0 * radius);
            let flat = |p: Point| p - direction * (p - s.at).dot(direction);
            Section { t: s.t + step, at: s.at + by, center: flat(s.center) + by, contacts: [flat(s.contacts[0]) + by, flat(s.contacts[1]) + by], uv: [None, None] }
        };
        // (Its parameter on as fast as it runs at the end.)
        let rate = |a: &Section, b: &Section| (b.t - a.t) / a.at.distance(b.at).max(1e-300);
        let first = carry(&sections[0], -tangents[0], -2.0 * radius * rate(&sections[0], &sections[1]));
        let last = carry(&sections[m], tangents[1], 2.0 * radius * rate(&sections[m - 1], &sections[m]));
        let secs = std::iter::once(first).chain(sections.iter().copied()).chain([last]).collect();
        let supp = std::iter::once(supports[0]).chain(supports.iter().copied()).chain([supports[m]]).collect();
        (secs, supp, breaks.iter().map(|&b| b + 1).chain([0]).collect())
    };
    let m = secs.len() - 1;
    starts.sort_unstable();
    starts.dedup();
    let stretches: Vec<(usize, usize)> = starts.iter().enumerate().map(|(i, &a)| (a, starts.get(i + 1).copied().unwrap_or(m))).collect();
    let count = stretches.len();
    // The boundaries: before each stretch, and (along an open edge) after the last.
    let bounds = if closed { count } else { count + 1 };
    let next = |i: usize| (i + 1) % bounds;
    // (What the contacts lie on along a stretch: at its sections past the switch.)
    let support_of = |i: usize| -> [Support; 2] {
        let (a, b) = stretches[i];
        supp[(a + 1).min(b)]
    };

    // The fillet surfaces.
    let patches = stretches.iter().map(|&(a, b)| fillet_surface(&secs[a..=b])).collect::<Result<Vec<_>>>()?;
    // The contacts' pcurves on the faces they run on, along each stretch (over the edge's
    // parameter, as the fillet's sides are).
    let pcurve_along = |i: usize, k: usize| -> Option<Curve2> {
        let (a, b) = stretches[i];
        let (f, _) = support_of(i)[k];
        let face = &sides[k].0[f];
        let mut points: Vec<Point2> = Vec::with_capacity(b - a + 1);
        for s in &secs[a..=b] {
            let mut q = match s.uv[k] {
                Some(q) => q,
                // (At a switch, where the contact is on the face's edge, and past the ends.)
                None => {
                    let (u, v, _) = nearest_on_surface(s.contacts[k], &face.surface, points.last().map(|p| (p.x, p.y)));
                    Point2::new(u, v)
                }
            };
            if let Some(&last) = points.last() {
                let unwrap = |x: f64, previous: f64, period: Option<f64>| period.map_or(x, |p| x + ((previous - x) / p).round() * p);
                q = Point2::new(unwrap(q.x, last.x, face.u_period), unwrap(q.y, last.y, face.v_period));
            }
            points.push(q);
        }
        let ts: Vec<f64> = secs[a..=b].iter().map(|s| s.t).collect();
        crate::geometry::fit::interpolate_with_parameters(&points, &ts).ok().map(Curve2::BSpline)
    };
    let edge_id = |s: Shape| EdgeId::try_from(s.node).map_err(|_| failed("an edge that is not one"));
    let mut pcurves: HashMap<EdgeId, Curve2> = HashMap::default();
    let mut runs: HashMap<crate::topology::NodeId, Vec<(Vec<EdgeId>, Kept)>> = HashMap::default();
    let mut pieces: HashMap<EdgeId, Vec<EdgeId>> = HashMap::default();

    // The vertices at the switches; on the edges the contacts switch on, the edges split
    // there.
    let mut boundaries: Vec<Option<Boundary>> = (0..bounds).map(|_| None).collect();
    let mut switch_vertices: Vec<Option<[VertexId; 2]>> = vec![None; bounds];
    let mut splits: Vec<(EdgeId, f64, VertexId)> = Vec::new();
    for i in 0..bounds {
        if !closed && (i == 0 || i == count) {
            continue;
        }
        let mut pair = Vec::with_capacity(2);
        for k in 0..2 {
            let p = secs[stretches[i].0].contacts[k];
            let (before, after) = (support_of((i + count - 1) % count)[k], support_of(i)[k]);
            let on_edge = if before != after || before.1.is_some() || after.1.is_some() {
                sides[k].0.iter().flat_map(|f| f.edges.iter().zip(&f.edge_shapes)).find_map(|((curve, a, b), &e)| {
                    let foot = nearest_on_curve(p, curve, *a, *b);
                    (foot.distance <= near).then_some((e, foot.parameter))
                })
            } else {
                None
            };
            pair.push(match on_edge {
                Some((e, t)) => {
                    let e = edge_id(e)?;
                    let (a, b) = model.edge_end_vertices(e);
                    if let Some(&v) = [a, b].iter().find(|&&v| model.point(Shape::from(v)).distance(p) <= near) {
                        v
                    } else {
                        match splits.iter().find(|s| s.0 == e && model.point(Shape::from(s.2)).distance(p) <= near) {
                            Some(s) => s.2,
                            None => {
                                // (On the edge, as near the contact as it is.)
                                let (curve, ..) = model.edge_curve(Shape::from(e)).ok_or_else(|| failed("an edge without a curve"))?;
                                let q = curve.point(t);
                                let v = model.add_vertex(q, (2.0 * q.distance(p)).max(TOL));
                                splits.push((e, t, v));
                                v
                            }
                        }
                    }
                }
                None => model.add_vertex(p, TOL),
            });
        }
        switch_vertices[i] = Some([pair[0], pair[1]]);
    }
    splits.sort_by(|a, b| (a.0.index(), a.1).partial_cmp(&(b.0.index(), b.1)).expect("parameters"));
    for (e, t, v) in splits {
        split_piece(model, &mut pieces, e, t, v)?;
    }
    // The arcs between stretches.
    for i in 0..bounds {
        let Some(vertices) = switch_vertices[i] else { continue };
        let patch = &patches[i];
        let e = model.add_edge(vertices[0], vertices[1], TOL)?;
        model.set_edge_curve(e.into(), patch.surface.v_iso(patch.v.0), TOL)?;
        model.set_edge_range(e.into(), patch.u.0, patch.u.1)?;
        let before = &patches[(i + count - 1) % count];
        boundaries[i] = Some(Boundary { vertices, edge: e, v: [patch.v.0, patch.v.0], across: None });
        let _ = before;
    }
    // At the ends of an open edge, where the fillet crosses the face there.
    if !closed {
        let (va, vb) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
        for (end, vertex) in [(0, va), (count, vb)] {
            let stretch = if end == 0 { 0 } else { count - 1 };
            if support_of(stretch).iter().any(|s| *s != (0, None)) {
                return Err(failed("an end where the ball touches more than the edge's faces"));
            }
            let cut = end_cut(model, shape, edge, vertex, sides, &patches[stretch], end == 0, near)?;
            // The edges at the end split where the contacts meet them.
            let mut vertices = Vec::with_capacity(2);
            for (e, p) in [cut.a, cut.b] {
                let (curve, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| failed("an edge without a curve"))?;
                let foot = nearest_on_curve(p, &curve, first, last);
                if foot.distance > near.max(1e3 * TOL) {
                    return Err(failed("the fillet meeting an edge at its end off it"));
                }
                let v = model.add_vertex(foot.point, (2.0 * foot.distance).max(TOL));
                split_piece(model, &mut pieces, e, foot.parameter, v)?;
                vertices.push(v);
            }
            let vertices = [vertices[0], vertices[1]];
            let e = model.add_edge(vertices[0], vertices[1], TOL)?;
            model.set_edge_curve(e.into(), Curve::BSpline(cut.curve), TOL)?;
            model.set_edge_range(e.into(), 0.0, 1.0)?;
            raise_tolerance(model, e, cut.deviation)?;
            pcurves.insert(e, Curve2::BSpline(cut.on_face));
            let s = &secs[if end == 0 { 1 } else { m - 1 }];
            runs.entry(cut.face.node).or_default().push((vec![e], Kept::Ball { center: s.center, toward: convex }));
            boundaries[end] = Some(Boundary { vertices, edge: e, v: cut.v, across: Some(Curve2::BSpline(cut.on_fillet)) });
        }
    }
    let boundaries: Vec<Boundary> = boundaries.into_iter().map(|b| b.expect("a boundary")).collect();

    // The contact curves along each stretch.
    let mut contacts: Vec<[Contact; 2]> = Vec::with_capacity(count);
    for i in 0..count {
        let patch = &patches[i];
        let (from, to) = (&boundaries[i], &boundaries[next(i)]);
        let mut pair = Vec::with_capacity(2);
        for k in 0..2 {
            let (a, b) = (from.vertices[k], to.vertices[k]);
            pair.push(match support_of(i)[k] {
                (f, None) => {
                    let curve = patch.surface.u_iso(if k == 0 { patch.u.0 } else { patch.u.1 });
                    let e = model.add_edge(a, b, TOL)?;
                    model.set_edge_curve(e.into(), curve, TOL)?;
                    // (From where the fillet's faces end on this side: at the arcs, the
                    // patch's ends.)
                    let range = (if from.across.is_some() { from.v[k] } else { patch.v.0 }, if to.across.is_some() { to.v[k] } else { patch.v.1 });
                    model.set_edge_range(e.into(), range.0, range.1)?;
                    let pcurve = pcurve_along(i, k).ok_or_else(|| failed("no pcurve of a contact curve"))?;
                    pcurves.insert(e, pcurve);
                    Contact::OnFace(f, vec![e])
                }
                (f, Some(r)) => {
                    let original = edge_id(sides[k].0[f].edge_shapes[r])?;
                    let candidates = pieces.get(&original).cloned().unwrap_or_else(|| vec![original]);
                    let piece = candidates
                        .into_iter()
                        .find(|&p| {
                            let (x, y) = model.edge_end_vertices(p);
                            (x == a && y == b) || (x == b && y == a)
                        })
                        .ok_or_else(|| failed("no piece of an edge the ball rolls over between its switches"))?;
                    Contact::OnEdge(piece, model.edge_end_vertices(piece).0 == a)
                }
            });
        }
        contacts.push([pair[0].clone(), pair[1].clone()]);
    }

    // Where contact curves on a periodic face cross its seam (in its parameters): both split
    // there.
    for i in 0..count {
        for k in 0..2 {
            let Contact::OnFace(f, ref mut along) = contacts[i][k] else { continue };
            let face = &sides[k].0[f];
            for seam in seams(model, face.face) {
                // (The seam's u, or v, where the face is periodic that way.)
                let Some(across) = seam_parameter(model, face.face, seam, (face.u_period, face.v_period)) else { continue };
                let mut j = 0;
                while j < along.len() {
                    let pcurve = pcurves[&along[j]].clone();
                    let (_, first, last) = model.edge_curve(Shape::from(along[j])).ok_or_else(|| failed("a contact edge without a curve"))?;
                    let Some(t) = crossing(&pcurve, first, last, across) else {
                        j += 1;
                        continue;
                    };
                    // (On the seam: the face's point at the pcurve's crossing.)
                    let (curve, ..) = model.edge_curve(Shape::from(along[j])).expect("a curve");
                    let uv = pcurve.point(t);
                    let p = face.surface.point(uv.x, uv.y);
                    let off = p.distance(curve.point(t));
                    let seam_now = pieces.get(&seam).cloned().unwrap_or_else(|| vec![seam]);
                    let (gap, ts) = seam_now
                        .iter()
                        .find_map(|&s| {
                            let (c, a, b) = model.edge_curve(Shape::from(s))?;
                            let foot = nearest_on_curve(p, &c, a, b);
                            (foot.distance <= near).then_some((foot.distance, foot.parameter))
                        })
                        .ok_or_else(|| failed("a seam crossing off the seam"))?;
                    // (Its tolerance covering the seam, as near as it passes.)
                    let v = model.add_vertex(p, (2.0 * gap.max(off)).max(TOL));
                    let (h1, h2) = model.split_edge(along[j], v, t)?;
                    pcurves.insert(h1, pcurve.clone());
                    pcurves.insert(h2, pcurve);
                    along.splice(j..=j, [h1, h2]);
                    split_piece(model, &mut pieces, seam, ts, v)?;
                    j += 2;
                }
            }
        }
    }

    // The runs of contact curves on each face, each with the part away from the edge kept.
    for k in 0..2 {
        let on = |i: usize| match &contacts[i][k] {
            Contact::OnFace(f, _) => Some(*f),
            Contact::OnEdge(..) => None,
        };
        let all_same = (0..count).all(|i| on(i).is_some() && on(i) == on(0));
        let begin = if !closed || all_same { 0 } else { (0..count).find(|&i| on(i) != on((i + count - 1) % count)).unwrap_or(0) };
        let mut i = 0;
        while i < count {
            let at = (begin + i) % count;
            let Some(f) = on(at) else {
                i += 1;
                continue;
            };
            let mut edges = Vec::new();
            while i < count && on((begin + i) % count) == Some(f) {
                if let Contact::OnFace(_, along) = &contacts[(begin + i) % count][k] {
                    edges.extend(along.iter().copied());
                }
                i += 1;
            }
            runs.entry(sides[k].0[f].face.node).or_default().push((edges, Kept::AwayFromEdge(k)));
        }
    }

    // The sections, for which side of the contacts a point of a face lies on.
    let removed = |q: Point| -> bool {
        let s = secs[..m].iter().min_by(|a, b| a.at.distance(q).total_cmp(&b.at.distance(q))).expect("sections");
        let (a, b, w, e) = (s.contacts[0] - s.center, s.contacts[1] - s.center, q - s.center, s.at - s.center);
        let n = a.cross(b);
        let slack = -1e-9 * radius * radius * radius;
        n.dot(a.cross(w)) >= slack && n.dot(w.cross(b)) >= slack && w.length() <= e.length() + near
    };

    // The faces: those the runs cross rebuilt, those the ball reaches without dropped where
    // they lie between the edge and the contacts, the rest as they are.
    let reached: Vec<crate::topology::NodeId> = sides.iter().flat_map(|s| s.0.iter().map(|f| f.face.node)).collect();
    let mut faces: Vec<Shape> = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        if let Some(face_runs) = runs.get(&face.node) {
            faces.extend(rebuilt(model, face, face_runs, &pcurves, &secs, &removed, near)?);
        } else if reached.contains(&face.node) {
            let surface = model.surface(face).cloned().ok_or_else(|| failed("a face without a surface"))?;
            let (u0, u1, v0, v1) = model.face_uv_bounds(face);
            if !removed(surface.point(0.5 * (u0 + u1), 0.5 * (v0 + v1))) {
                faces.push(face);
            }
        } else {
            faces.push(face);
        }
    }

    // The fillet faces, each bounded by the arcs (or crossings) at its ends and the contact
    // curves along it.
    for i in 0..count {
        let patch = &patches[i];
        let (u0, u1, v0, v1) = (patch.u.0, patch.u.1, patch.v.0, patch.v.1);
        let line = |origin: Point2, along_u: bool| Curve2::Line(Line2::new(origin, if along_u { Direction2::new(1.0, 0.0) } else { Direction2::new(0.0, 1.0) }));
        let (from, to) = (&boundaries[i], &boundaries[next(i)]);
        let mut uses: Vec<(EdgeId, bool, Curve2)> = vec![(from.edge, true, from.across.clone().unwrap_or_else(|| line(Point2::new(0.0, v0), true)))];
        let side = |model: &mut Model, k: usize| -> Result<Vec<(EdgeId, bool, Curve2)>> {
            Ok(match &contacts[i][k] {
                Contact::OnFace(_, along) => along.iter().map(|&e| (e, true, line(Point2::new(if k == 0 { u0 } else { u1 }, 0.0), false))).collect(),
                Contact::OnEdge(e, forward) => {
                    let (curve, first, last) = model.edge_curve(Shape::from(*e)).ok_or_else(|| failed("an edge without a curve"))?;
                    let (pcurve, deviation) = pcurve_along_side(&curve, first, last, patch, if k == 0 { u0 } else { u1 }).ok_or_else(|| failed("no pcurve on the fillet of an edge it rolls over"))?;
                    raise_tolerance(model, *e, deviation)?;
                    vec![(*e, *forward, pcurve)]
                }
            })
        };
        uses.extend(side(model, 1)?);
        uses.push((to.edge, false, to.across.clone().unwrap_or_else(|| line(Point2::new(0.0, v1), true))));
        uses.extend(side(model, 0)?.into_iter().rev().map(|(e, forward, p)| (e, !forward, p)));
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, patch.surface.clone(), TOL)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (e, forward, pcurve) in uses {
            let c = model.add_coedge(e, if forward { Orientation::Forward } else { Orientation::Reversed })?;
            let (_, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| failed("an edge without a curve"))?;
            model.set_pcurve_with_range(c, pcurve, first, last)?;
            coedges.push(c);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, wire.into())?;
        // Facing out of the solid: toward the ball's centre at a concave edge, away at a
        // convex one.
        let (p, du, dv) = patch.surface.d1(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let s = &secs[(stretches[i].0 + stretches[i].1) / 2];
        let toward = du.cross(dv).dot(s.center - p) > 0.0;
        faces.push(if toward != convex { face } else { face.reversed() });
    }

    let shell = model.new_empty(NodeKind::Shell);
    for f in faces {
        model.add(shell, f)?;
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    Ok(solid)
}

/// Where the fillet crosses the face at an end of its edge: that face, the crossing (as a
/// curve over `u` in [0, 1], on the fillet's patch and on the face), how far those stray,
/// the `v` where it meets the fillet's sides, and the edges at the end of the faces along
/// the edge with where it meets them.
struct EndCut {
    face: Shape,
    curve: crate::geometry::curve::BSplineCurve,
    on_fillet: crate::geometry::curve2::BSplineCurve2,
    on_face: crate::geometry::curve2::BSplineCurve2,
    deviation: f64,
    v: [f64; 2],
    a: (EdgeId, Point),
    b: (EdgeId, Point),
}

/// The crossing of `patch` (the first stretch's `at_start`, else the last's)
/// with the face at `vertex`, the end of `edge`: the face whose boundary holds the other
/// edges at the vertex of the two faces along the edge.
#[allow(clippy::too_many_arguments)]
fn end_cut(model: &Model, shape: Shape, edge: Shape, vertex: Shape, sides: &[Side; 2], patch: &Patch, at_start: bool, near: f64) -> Result<EndCut> {
    let failed = |what: &str| Error::Failed(format!("a run-over fillet's end: {what}"));
    let at_vertex = |e: Shape| {
        let (a, b) = model.edge_vertices(e);
        a.node == vertex.node || b.node == vertex.node
    };
    let others = |face: Shape| -> Vec<Shape> { model.unique_subshapes(face, NodeKind::Edge).into_iter().filter(|&e| e.node != edge.node && at_vertex(e)).collect() };
    let (ea, eb) = (others(sides[0].0[0].face), others(sides[1].0[0].face));
    let ([a], [b]) = (&ea[..], &eb[..]) else { return Err(failed("not three edges at the end")) };
    let ends: Vec<Shape> = model
        .unique_subshapes(shape, NodeKind::Face)
        .into_iter()
        .filter(|&f| f.node != sides[0].0[0].face.node && f.node != sides[1].0[0].face.node)
        .filter(|&f| {
            let edges = model.unique_subshapes(f, NodeKind::Edge);
            edges.iter().any(|e| e.node == a.node) && edges.iter().any(|e| e.node == b.node)
        })
        .collect();
    let [face] = ends[..] else { return Err(failed("not one face at the end")) };
    let surface = model.surface(face).cloned().ok_or_else(|| failed("a face without a surface"))?;
    let hint = std::cell::Cell::new(None);
    let distance = |p: Point| {
        let (u, v, _) = nearest_on_surface(p, &surface, hint.get());
        hint.set(Some((u, v)));
        let (foot, du, dv) = surface.d1(u, v);
        (p - foot).dot(du.cross(dv))
    };
    let (v0, v1) = patch.v;
    // (From the patch's end past the edge's end inward: the first crossing.)
    let vs: Vec<f64> = {
        let (from, to) = if at_start { (v0, v1) } else { (v1, v0) };
        (0..=128).map(|j| from + (to - from) * j as f64 / 128.0).collect()
    };
    let cross = |u: f64| -> Option<f64> {
        let d = |v: f64| distance(patch.surface.point(u, v));
        let mut previous = (vs[0], d(vs[0]));
        for &v in &vs[1..] {
            let here = (v, d(v));
            if (previous.1 <= 0.0) != (here.1 <= 0.0) {
                let (mut lo, mut hi) = (previous.0, here.0);
                let low = previous.1 <= 0.0;
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if (d(mid) <= 0.0) == low { lo = mid } else { hi = mid }
                }
                return Some(0.5 * (lo + hi));
            }
            previous = here;
        }
        None
    };
    let mut n = 16;
    loop {
        let us: Vec<f64> = (0..=n).map(|j| j as f64 / n as f64).collect();
        let vs: Vec<f64> = us.iter().map(|&u| cross(u)).collect::<Option<_>>().ok_or_else(|| failed("the fillet not reaching across the face at the end"))?;
        let points: Vec<Point> = us.iter().zip(&vs).map(|(&u, &v)| patch.surface.point(u, v)).collect();
        let curve = crate::geometry::fit::interpolate_with_parameters(&points, &us)?;
        let on_fillet = crate::geometry::fit::interpolate_with_parameters(&us.iter().zip(&vs).map(|(&u, &v)| Point2::new(u, v)).collect::<Vec<_>>(), &us)?;
        let mut face_uv: Vec<Point2> = Vec::with_capacity(points.len());
        for &p in &points {
            let (u, v, _) = nearest_on_surface(p, &surface, face_uv.last().map(|q: &Point2| (q.x, q.y)));
            let mut q = Point2::new(u, v);
            if let Some(&last) = face_uv.last() {
                let unwrap = |x: f64, previous: f64, period: Option<f64>| period.map_or(x, |p| x + ((previous - x) / p).round() * p);
                q = Point2::new(unwrap(q.x, last.x, surface.u_period()), unwrap(q.y, last.y, surface.v_period()));
            }
            face_uv.push(q);
        }
        let on_face = crate::geometry::fit::interpolate_with_parameters(&face_uv, &us)?;
        let mut deviation: f64 = 0.0;
        for w in us.windows(2) {
            let u = 0.5 * (w[0] + w[1]);
            let p = curve.point(u);
            let (q, r) = (on_fillet.point(u), on_face.point(u));
            deviation = deviation.max(patch.surface.point(q.x, q.y).distance(p)).max(surface.point(r.x, r.y).distance(p));
        }
        if deviation <= 0.1 * near || n >= 256 {
            let (first, last) = (points[0], points[n]);
            let id = |e: &Shape| EdgeId::try_from(e.node).map_err(|_| failed("an edge that is not one"));
            return Ok(EndCut { face, curve, on_fillet, on_face, deviation, v: [vs[0], vs[n]], a: (id(a)?, first), b: (id(b)?, last) });
        }
        n *= 2;
    }
}

/// The fillet surface through `secs`: across (u in [0, 1]) each section's arc from the first
/// contact to the second, a rational quadratic; along the edge (v, its parameter) a cubic
/// interpolation of the arcs' homogeneous control points.
fn fillet_surface(secs: &[Section]) -> Result<Patch> {
    use crate::geometry::fit;
    use crate::geometry::surface::{BSplineSurface, KnotVector};
    let ts: Vec<f64> = secs.iter().map(|s| s.t).collect();
    let mut middle = Vec::with_capacity(secs.len());
    let mut weights = Vec::with_capacity(secs.len());
    for s in secs {
        let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        let w = half.cos();
        // The corner of the arc's control polygon, where its end tangents meet.
        let corner = s.center + (a + b).normalized() * (a.length() / w);
        middle.push(Point::ORIGIN + (corner - Point::ORIGIN) * w);
        weights.push(Point::new(w, 0.0, 0.0));
    }
    let rows = [
        fit::interpolate_with_parameters(&secs.iter().map(|s| s.contacts[0]).collect::<Vec<_>>(), &ts)?,
        fit::interpolate_with_parameters(&middle, &ts)?,
        fit::interpolate_with_parameters(&secs.iter().map(|s| s.contacts[1]).collect::<Vec<_>>(), &ts)?,
    ];
    let w = fit::interpolate_with_parameters(&weights, &ts)?;
    let w_poles: Vec<f64> = w.poles().iter().map(|p| p.x).collect();
    if w_poles.iter().any(|&w| w <= 0.0) {
        return Err(Error::InvalidInput("the fillet turns too sharply".into()));
    }
    let poles = vec![
        rows[0].poles().to_vec(),
        rows[1].poles().iter().zip(&w_poles).map(|(p, w)| Point::ORIGIN + (*p - Point::ORIGIN) * (1.0 / w)).collect(),
        rows[2].poles().to_vec(),
    ];
    let n = w_poles.len();
    let u = KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2);
    let v = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
    let surface = BSplineSurface::try_new(poles, u, v)?.with_weights(vec![vec![1.0; n], w_poles, vec![1.0; n]]);
    Ok(Patch { surface: Surface::BSpline(surface), u: (0.0, 1.0), v: (ts[0], ts[ts.len() - 1]) })
}

/// Splits whichever piece of `edge` (split before, as recorded in `pieces`) holds curve
/// parameter `t`, at `vertex`.
fn split_piece(model: &mut Model, pieces: &mut HashMap<EdgeId, Vec<EdgeId>>, edge: EdgeId, t: f64, vertex: VertexId) -> Result<()> {
    let list = pieces.entry(edge).or_insert_with(|| vec![edge]);
    let at = list
        .iter()
        .position(|&p| model.edge_curve(Shape::from(p)).is_some_and(|(_, a, b)| a.min(b) < t && t < a.max(b)))
        .ok_or_else(|| Error::Failed("a split off the edge's pieces".into()))?;
    let (h1, h2) = model.split_edge(list[at], vertex, t)?;
    list.splice(at..=at, [h1, h2]);
    Ok(())
}

/// The seam edges of a face: those its boundary runs along twice.
fn seams(model: &Model, face: Shape) -> Vec<EdgeId> {
    let mut count: HashMap<EdgeId, usize> = HashMap::default();
    for wire in model.children(Shape { orientation: Orientation::Forward, ..face }) {
        for coedge in model.children(wire) {
            if let Ok(id) = CoEdgeId::try_from(coedge.node) {
                *count.entry(model[id].edge).or_default() += 1;
            }
        }
    }
    let mut seams: Vec<EdgeId> = count.into_iter().filter(|&(_, n)| n >= 2).map(|(e, _)| e).collect();
    seams.sort_unstable_by_key(|e| e.index());
    seams
}

/// Where `pcurve` over `[first, last]` crosses a seam at `across` (`(u, period)` or
/// `(v, period)` as `Seam` says), away from its ends.
fn crossing(pcurve: &Curve2, first: f64, last: f64, across: Seam) -> Option<f64> {
    const N: usize = 64;
    let (at, period, along_u) = across;
    let x = |t: f64| {
        let p = pcurve.point(t);
        if along_u { p.x } else { p.y }
    };
    // (Which copy of the seam a value lies past.)
    let side = |t: f64| ((x(t) - at) / period).floor();
    let span = last - first;
    for i in 0..N {
        let (mut lo, mut hi) = (first + span * i as f64 / N as f64, first + span * (i + 1) as f64 / N as f64);
        if side(lo) == side(hi) {
            continue;
        }
        let low = side(lo);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if side(mid) == low { lo = mid } else { hi = mid }
        }
        let t = 0.5 * (lo + hi);
        if (t - first).abs() > 1e-9 * span.abs() && (last - t).abs() > 1e-9 * span.abs() {
            return Some(t);
        }
    }
    None
}

/// A seam's place across a periodic face: its u (or v), the period, and whether it is u.
type Seam = (f64, f64, bool);

/// Where `seam` lies across `face`, periodic as `periods` say.
fn seam_parameter(model: &Model, face: Shape, seam: EdgeId, periods: (Option<f64>, Option<f64>)) -> Option<Seam> {
    for wire in model.children(Shape { orientation: Orientation::Forward, ..face }) {
        for coedge in model.children(wire) {
            let Ok(id) = CoEdgeId::try_from(coedge.node) else { continue };
            if model[id].edge != seam {
                continue;
            }
            let (pcurve, first, last) = model.effective_pcurve(id)?;
            let (a, b) = (pcurve.point(first), pcurve.point(last));
            return match periods {
                (Some(p), _) if (a.x - b.x).abs() <= 1e-9 * (1.0 + p) => Some((a.x, p, true)),
                (_, Some(p)) if (a.y - b.y).abs() <= 1e-9 * (1.0 + p) => Some((a.y, p, false)),
                _ => None,
            };
        }
    }
    None
}

/// The pcurve on `patch` of `curve` over `[first, last]` lying along its side `u = side`
/// (an edge the ball rolls over), and how far it strays.
fn pcurve_along_side(curve: &Curve, first: f64, last: f64, patch: &Patch, side: f64) -> Option<(Curve2, f64)> {
    const N: usize = 32;
    let iso = patch.surface.u_iso(side);
    let ts: Vec<f64> = (0..=N).map(|i| first + (last - first) * i as f64 / N as f64).collect();
    let points: Vec<Point2> = ts.iter().map(|&t| Point2::new(side, nearest_on_curve(curve.point(t), &iso, patch.v.0, patch.v.1).parameter)).collect();
    let pcurve = Curve2::BSpline(crate::geometry::fit::interpolate_with_parameters(&points, &ts).ok()?);
    let deviation = (0..2 * N).map(|i| first + (last - first) * (i as f64 + 0.5) / (2 * N) as f64).map(|t| {
        let uv = pcurve.point(t);
        patch.surface.point(uv.x, uv.y).distance(curve.point(t))
    }).fold(0.0, f64::max);
    Some((pcurve, deviation))
}

fn raise_tolerance(model: &mut Model, edge: EdgeId, deviation: f64) -> Result<()> {
    let tolerance = (2.0 * deviation).max(TOL);
    if tolerance > model.tolerance(Shape::from(edge)) {
        model.set_edge_tolerance(edge, tolerance)?;
    }
    Ok(())
}

/// `face` rebuilt around `runs` of contact curves on it (each from side `k` of the
/// fillet): the faces of its kept part.
fn rebuilt(model: &mut Model, face: Shape, runs: &[(Vec<EdgeId>, Kept)], pcurves: &HashMap<EdgeId, Curve2>, secs: &[Section], removed: &dyn Fn(Point) -> bool, near: f64) -> Result<Vec<Shape>> {
    let failed = |what: &str| Error::Failed(format!("rebuilding a face along a run-over fillet: {what}"));
    let node = Shape { orientation: Orientation::Forward, ..face };
    let surface = model.surface(node).cloned().ok_or_else(|| failed("a face without a surface"))?;
    let (fu0, fu1, fv0, fv1) = model.face_uv_bounds(node);
    let periods = (surface.u_period(), surface.v_period());
    // A shift bringing `p` into the face's domain, in each periodic direction.
    let into_domain = |p: Point2| -> Vector2 {
        let shift = |x: f64, lo: f64, period: Option<f64>| period.map_or(0.0, |period| -((x - lo) / period).floor() * period);
        Vector2::new(shift(p.x, fu0, periods.0), shift(p.y, fv0, periods.1))
    };
    let uv_near = 1e-6 * (1.0 + (fu1 - fu0).hypot(fv1 - fv0));

    // The face's own loops, in its own frame (its material on their left).
    let mut loops: Vec<Vec<Use>> = Vec::new();
    for wire in model.children(node) {
        let mut uses = Vec::new();
        for coedge in model.children(wire) {
            let Ok(id) = CoEdgeId::try_from(coedge.node) else { continue };
            let (pcurve, first, last) = model.effective_pcurve(id).ok_or_else(|| failed("an edge without a pcurve"))?;
            let orientation = model[id].orientation.compose(coedge.orientation);
            uses.push(Use { edge: model[id].edge, forward: orientation != Orientation::Reversed, pcurve, first, last, own: true });
        }
        if wire.orientation == Orientation::Reversed {
            uses.reverse();
        }
        loops.push(uses);
    }

    // The runs as chains of uses, each with the kept part on its left, its pcurves on the
    // face continuous along it.
    let mut chains: Vec<Vec<Use>> = Vec::new();
    for (edges, kept) in runs {
        let mut uses = Vec::with_capacity(edges.len());
        for &e in edges {
            let (curve, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| failed("a contact edge without a curve"))?;
            let pcurve = pcurves.get(&e).cloned().ok_or_else(|| failed("a contact curve without its pcurve"))?;
            let deviation = (0..=16).map(|i| first + (last - first) * i as f64 / 16.0).map(|t| {
                let uv = pcurve.point(t);
                surface.point(uv.x, uv.y).distance(curve.point(t))
            }).fold(0.0, f64::max);
            raise_tolerance(model, e, deviation)?;
            uses.push(Use { edge: e, forward: true, pcurve, first, last, own: false });
        }
        // (Kept on the left: away from the edge, where from a section the contact lies away
        // from the edge; or toward the ball's centre or away.)
        let middle = &uses[uses.len() / 2];
        let (curve, first, last) = model.edge_curve(Shape::from(middle.edge)).expect("a curve");
        let (c, tangent) = curve.d1(0.5 * (first + last));
        let kept_way = match *kept {
            Kept::AwayFromEdge(k) => {
                let s = secs.iter().min_by(|a, b| a.contacts[k].distance(c).total_cmp(&b.contacts[k].distance(c))).expect("sections");
                c - s.at
            }
            Kept::Ball { center, toward } => if toward { center - c } else { c - center },
        };
        let (u, v, _) = nearest_on_surface(c, &surface, None);
        let (_, du, dv) = surface.d1(u, v);
        if du.cross(dv).cross(tangent).dot(kept_way) < 0.0 {
            uses.reverse();
            for u in &mut uses {
                u.forward = !u.forward;
            }
        }
        // (A closed run starting where it crosses a seam, if it does.)
        let closed = uses[0].start(model) == uses[uses.len() - 1].end(model);
        if closed {
            let seam_vertices: Vec<VertexId> = seams(model, node).into_iter().flat_map(|s| {
                let (a, b) = model.edge_end_vertices(s);
                [a, b]
            }).collect();
            if let Some(at) = uses.iter().position(|u| seam_vertices.contains(&u.start(model))) {
                uses.rotate_left(at);
            }
        }
        for j in 1..uses.len() {
            let end = uses[j - 1].uv(1.0);
            let start = uses[j].uv(0.0);
            let shift = |d: f64, period: Option<f64>| period.map_or(0.0, |p| (d / p).round() * p);
            uses[j] = uses[j].shifted(Vector2::new(shift(end.x - start.x, periods.0), shift(end.y - start.y, periods.1)));
        }
        chains.push(uses);
    }

    // Round the kept part: along a chain, then along the face's loop from where it ends to
    // where the next starts, and so on until back at the first.
    let ends: Vec<(VertexId, Point2)> = chains
        .iter()
        .map(|c| {
            let last = &c[c.len() - 1];
            (last.end(model), last.uv(1.0) + into_domain(last.uv(0.9)))
        })
        .collect();
    let starts: Vec<(VertexId, Point2)> = chains.iter().map(|c| (c[0].start(model), c[0].uv(0.0) + into_domain(c[0].uv(0.1)))).collect();
    let open: Vec<bool> = chains.iter().map(|c| c[0].start(model) != c[c.len() - 1].end(model) || loops.iter().flatten().any(|u| u.start(model) == c[0].start(model))).collect();
    let mut used = vec![false; chains.len()];
    let mut visited: Vec<Vec<bool>> = loops.iter().map(|l| vec![false; l.len()]).collect();
    let mut result: Vec<Vec<Use>> = Vec::new();
    for first in 0..chains.len() {
        if used[first] || !open[first] {
            continue;
        }
        let mut path: Vec<Use> = Vec::new();
        let mut current = first;
        'round: for _ in 0..=chains.len() {
            used[current] = true;
            path.extend(chains[current].iter().cloned());
            let (vertex, uv) = ends[current];
            let (li, mut j) = loops
                .iter()
                .enumerate()
                .flat_map(|(li, l)| l.iter().enumerate().map(move |(j, u)| (li, j, u)))
                .filter(|(_, _, u)| u.start(model) == vertex)
                .min_by(|a, b| a.2.uv(0.0).distance(uv).total_cmp(&b.2.uv(0.0).distance(uv)))
                .map(|(li, j, _)| (li, j))
                .ok_or_else(|| failed("a contact curve ending off the face's boundary"))?;
            for _ in 0..loops[li].len() {
                let u = &loops[li][j];
                visited[li][j] = true;
                path.push(u.clone());
                let (w, at) = (u.end(model), u.uv(1.0));
                let next = (0..chains.len()).filter(|&d| open[d] && starts[d].0 == w).min_by(|&a, &b| starts[a].1.distance(at).total_cmp(&starts[b].1.distance(at)));
                if let Some(d) = next.filter(|&d| starts[d].1.distance(at) <= uv_near.max(near)) {
                    if d == first {
                        break 'round;
                    }
                    current = d;
                    continue 'round;
                }
                j = (j + 1) % loops[li].len();
            }
            return Err(failed("a face's boundary not leading to a contact curve"));
        }
        result.push(path);
    }
    // Closed runs off the boundary are loops of their own; loops no run meets are kept
    // unless they lie between the edge and the contacts.
    for (c, chain) in chains.iter().enumerate() {
        if !open[c] {
            result.push(chain.clone());
        }
    }
    let run_vertices: Vec<VertexId> = chains.iter().flat_map(|c| [c[0].start(model), c[c.len() - 1].end(model)]).collect();
    for (li, l) in loops.iter().enumerate() {
        if visited[li].iter().any(|&v| v) || l.iter().any(|u| run_vertices.contains(&u.start(model))) {
            continue;
        }
        let (curve, first, last) = model.edge_curve(Shape::from(l[0].edge)).ok_or_else(|| failed("an edge without a curve"))?;
        if !removed(curve.point(0.5 * (first + last))) {
            result.push(l.clone());
        }
    }

    // Each loop's pcurves continuous round it, from one of the face's own (in its domain);
    // outer loops (counter-clockwise) make faces, holes go in the one around them.
    let mut polygons: Vec<(Vec<Use>, Vec<Point2>, f64)> = Vec::new();
    for mut l in result {
        if let Some(at) = l.iter().position(|u| u.own) {
            l.rotate_left(at);
        } else {
            let shift = into_domain(l[0].uv(0.5));
            l[0] = l[0].shifted(shift);
        }
        for j in 1..l.len() {
            let (end, start) = (l[j - 1].uv(1.0), l[j].uv(0.0));
            let shift = |d: f64, period: Option<f64>| period.map_or(0.0, |p| (d / p).round() * p);
            l[j] = l[j].shifted(Vector2::new(shift(end.x - start.x, periods.0), shift(end.y - start.y, periods.1)));
        }
        let points: Vec<Point2> = l.iter().flat_map(|u| (0..16).map(move |i| u.uv(i as f64 / 16.0))).collect();
        let area = 0.5 * (0..points.len()).map(|i| {
            let (a, b) = (points[i], points[(i + 1) % points.len()]);
            a.x * b.y - b.x * a.y
        }).sum::<f64>();
        polygons.push((l, points, area));
    }
    let (outer, holes): (Vec<_>, Vec<_>) = polygons.into_iter().partition(|p| p.2 > 0.0);
    if outer.is_empty() {
        return Err(failed("no part of a face kept"));
    }
    let mut groups: Vec<Vec<Vec<Use>>> = outer.iter().map(|o| vec![o.0.clone()]).collect();
    for (l, points, _) in holes {
        let within = |polygon: &[Point2], q: Point2| {
            let mut inside = false;
            for i in 0..polygon.len() {
                let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
                if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                    inside = !inside;
                }
            }
            inside
        };
        let at = (0..outer.len()).find(|&o| within(&outer[o].1, points[0])).unwrap_or(0);
        groups[at].push(l);
    }
    let tolerance = model.tolerance(node);
    let mut faces = Vec::with_capacity(groups.len());
    for group in groups {
        let f = model.new_empty(NodeKind::Face);
        model.set_face_surface(f, surface.clone(), tolerance)?;
        for l in group {
            let mut coedges = Vec::with_capacity(l.len());
            for u in l {
                let c = model.add_coedge(u.edge, if u.forward { Orientation::Forward } else { Orientation::Reversed })?;
                model.set_pcurve_with_range(c, u.pcurve, u.first, u.last)?;
                coedges.push(c);
            }
            let wire = model.add_wire(&coedges)?;
            model.add(f, wire.into())?;
        }
        faces.push(Shape { orientation: face.orientation, ..f });
    }
    Ok(faces)
}
