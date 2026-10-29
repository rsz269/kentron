//! Rolling-ball fillets of several edges built together: each ends on the face at its end
//! (as a single fillet does) or, where three filleted edges meet at a vertex, at a corner
//! patch spanning the three fillets' ends.
//!
//! At such a corner each of the three faces holds the contact curves of the two fillets
//! along its edges; they cross at a point, where both stop on that face. Each fillet stops
//! along the curve across it joining its two points (straight in its parameters), and the
//! patch (a Coons surface through the three such curves) closes the corner, as OCCT's
//! three-fillet corner does with its filling surface.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{BSplineCurve2, Curve2, Line2};
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{fit, project, ParametricCurve, ParametricSurface};
use crate::math::{Direction2, Point, Point2};
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{contact_vertex, cover_pcurve, joined_extensions, off_patch, reach_curve, reach_snapped, end_cut, extension, fillet_surface, new_face, rebuild_face, rolling_with, unwrapped, EndCut, Profile, Rolling, Use, NOT_REACHED};

/// Samples along a fillet's stop curve at a corner, and across a corner patch.
const ACROSS: usize = 32;
const PATCH: usize = 16;

/// One fillet: its rolling data, surface and spine parameters.
struct Stripe {
    rolling: Rolling,
    surface: Surface,
    ts: Vec<f64>,
}

/// How a fillet ends at one of its edge's vertices.
enum End {
    /// On the face there (the faces' other edges at the vertex shortened).
    Face(EndCut),
    /// At an open shell's free boundary: each contact curve stops where it crosses its
    /// face's free edge at the vertex (those edges shortened there).
    Free { vertex: VertexId, v: [f64; 2], edges: [(EdgeId, f64); 2] },
    /// At a mitre: two fillets meeting along their surfaces' intersection, from the point
    /// their contacts share on their common face to where one reaches its far face; the
    /// other then crosses that face, to the edge between their far faces.
    /// `pieces`: the mitre's curves (0 the intersection, 1 the cut beyond it) in order from
    /// the first side to the second, with the coedge's orientation that way and the pcurve
    /// on this fillet; `at` the mitre's points at the two sides.
    Mitre { mitre: usize, v: [f64; 2], at: [usize; 2], pieces: Vec<(usize, Orientation, BSplineCurve2)> },
    /// At the corner patch `corner`: stopping on each side at the point the two fillets
    /// on that face share (by face: the point's index in the corner).
    Corner { corner: usize, at: [usize; 2], v: [f64; 2] },
    /// Narrowed to nothing at its edge's vertex, where the faces are tangent: both contacts
    /// end there (at `v`), the fillet's face closing on a degenerate edge.
    Point { vertex: VertexId, v: [f64; 2] },
}

/// Three fillets meeting at a vertex: the points where they stop on each face (with the
/// face) and the fillets (by index) with which end of each.
struct Corner {
    points: Vec<(Shape, Point)>,
    stripes: Vec<(usize, usize)>,
}

/// Two fillets meeting at a vertex where three faces meet: their intersection curve (over
/// [0, 1]) from the point their contacts share on their common face to the point where the
/// first to get there reaches its far face, and the other's crossing of that face from
/// there to the edge between the two far faces (none where both get there together).
struct Mitre {
    /// The common face's point, the intersection's end and (with a cut) the cut's end.
    points: Vec<Point>,
    curve: BSplineCurve,
    cut: Option<MitreCut>,
    /// The edges the mitre shortens: each with its vertex at the corner, the parameter of
    /// the mitre's point it now ends at, and that point's index. The edge between the far
    /// faces, to the last point; none where the fillets share both faces, joined across a
    /// section; where their faces are tangent pairwise instead, the edges between each pair
    /// to the points at their sides.
    edges: Vec<(EdgeId, VertexId, f64, usize)>,
    deviation: f64,
}

/// The longer fillet's crossing of the shorter one's far face.
struct MitreCut {
    curve: BSplineCurve,
    face: Shape,
    on_face: Curve2,
}

/// `shape` with `edges` (meeting each other only three at a time, at corners where three
/// faces meet) rounded by rolling-ball fillets of the given radii (as laws over the fraction
/// of each edge's length), built together.
pub(in crate::blend) fn fillet_network(model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)]) -> Result<Shape> {
    let profiles: Vec<(Shape, Profile)> = edges.iter().map(|&(edge, radius)| (edge, Profile::Ball(radius))).collect();
    blend_network(model, shape, &profiles)
}

/// [`fillet_network`] for fillets or chamfers (each edge's profile).
pub(in crate::blend) fn blend_network(model: &mut Model, shape: Shape, edges: &[(Shape, Profile)]) -> Result<Shape> {
    // (Further past the ends where the end faces are not reached.)
    match blend_network_past(model, shape, edges, 1) {
        Err(Error::InvalidInput(m)) if m == NOT_REACHED => blend_network_past(model, shape, edges, 3),
        other => other,
    }
}

/// [`blend_network`] with the sections carried `far` times as far past the edges' ends.
fn blend_network_past(model: &mut Model, shape: Shape, edges: &[(Shape, Profile)], far: usize) -> Result<Shape> {
    let mut stripes = Vec::with_capacity(edges.len());
    for (edge, profile) in edges {
        let edge = *edge;
        let rolling = rolling_with(model, shape, edge, profile, far)?;
        if rolling.closed {
            return Err(Error::InvalidInput("closed edges in a network of fillets are not supported yet".into()));
        }
        let surface = fillet_surface(&rolling.sections)?;
        let ts = rolling.sections.iter().map(|s| s.t).collect();
        stripes.push(Stripe { rolling, surface, ts });
    }
    // The fillets' ends at each vertex.
    let mut at_vertex: HashMap<VertexId, Vec<(usize, usize)>> = HashMap::new();
    for (i, s) in stripes.iter().enumerate() {
        for (k, v) in [s.rolling.va, s.rolling.vb].into_iter().enumerate() {
            at_vertex.entry(VertexId::try_from(v.node).expect("a vertex")).or_default().push((i, k));
        }
    }
    let mut ends: Vec<[Option<End>; 2]> = stripes.iter().map(|_| [None, None]).collect();
    let mut corners: Vec<Corner> = Vec::new();
    let mut mitres: Vec<Mitre> = Vec::new();
    let mut vertices: Vec<(VertexId, Vec<(usize, usize)>)> = at_vertex.into_iter().collect();
    vertices.sort_by_key(|(v, _)| v.index());
    for (vertex, meeting) in vertices {
        match meeting[..] {
            [(i, k)] => {
                let s = &stripes[i];
                if s.rolling.vanishing[k] {
                    let t = if k == 0 { s.rolling.t0 } else { s.rolling.t1 };
                    ends[i][k] = Some(End::Point { vertex, v: [t, t] });
                    continue;
                }
                if let Some(end) = free_end(model, shape, s, k)? {
                    ends[i][k] = Some(end);
                    continue;
                }
                let cut = end_cut(model, shape, s.rolling.id, if k == 0 { s.rolling.va } else { s.rolling.vb }, &s.rolling.sides, &s.surface, &s.rolling.scan, k == 1, s.rolling.target)?;
                ends[i][k] = Some(End::Face(cut));
            }
            [a, b] => {
                // (Edges going on into one another: the fillets join across, else they meet
                // in a mitre.)
                let (m, ea, eb) = match joint(model, shape, vertex, &stripes, a, b) {
                    Some(joined) => joined?,
                    None => mitre(model, shape, vertex, &stripes, a, b)?,
                };
                ends[a.0][a.1] = Some(ea.into_end(mitres.len()));
                ends[b.0][b.1] = Some(eb.into_end(mitres.len()));
                mitres.push(m);
            }
            [_, _, _] => {
                let corner = corner(model, shape, vertex, &stripes, &meeting)?;
                for (n, &(i, k)) in meeting.iter().enumerate() {
                    ends[i][k] = Some(corner.1[n].take_end(corners.len()));
                }
                corners.push(corner.0);
            }
            _ => return Err(Error::InvalidInput("fillets meeting other than three at a corner are not supported here".into())),
        }
    }
    let ends: Vec<[End; 2]> = ends.into_iter().map(|[a, b]| [a.expect("an end"), b.expect("an end")]).collect();
    build(model, shape, &stripes, &ends, &corners, &mitres)
}

/// A fillet's end at a mitre, before the mitre is numbered.
struct PendingMitre {
    v: [f64; 2],
    at: [usize; 2],
    pieces: Vec<(usize, Orientation, BSplineCurve2)>,
}

impl PendingMitre {
    fn into_end(self, mitre: usize) -> End {
        End::Mitre { mitre, v: self.v, at: self.at, pieces: self.pieces }
    }
}

/// u on a fillet `w` of the way across from side `c` (0 at u = 0).
fn across(c: usize, w: f64) -> f64 {
    if c == 0 { w } else { 1.0 - w }
}

/// Solves `s(u, x) = l(y, z)` for `(x, y, z)` from `guess` (Newton, steps halved while they
/// do not bring the points closer).
fn meet(s: &Stripe, u: f64, l: &Stripe, guess: [f64; 3]) -> Option<[f64; 3]> {
    let gap = |p: [f64; 3]| -> [f64; 3] {
        let d = s.surface.point(u, p[0]) - l.surface.point(p[1], p[2]);
        [d.x, d.y, d.z]
    };
    let norm = |f: [f64; 3]| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
    let scale = 1.0 + (s.surface.point(u, guess[0]) - Point::ORIGIN).length();
    let mut p = guess;
    let mut f = gap(p);
    for _ in 0..100 {
        if norm(f) <= 1e-11 * scale {
            return Some(p);
        }
        let h = 1e-7;
        let col = |i: usize| {
            let mut q = p;
            q[i] += h;
            let g = gap(q);
            [(g[0] - f[0]) / h, (g[1] - f[1]) / h, (g[2] - f[2]) / h]
        };
        let step = super::super::solve3([col(0), col(1), col(2)], f)?;
        let mut t = 1.0;
        loop {
            let q = [p[0] - t * step[0], p[1] - t * step[1], p[2] - t * step[2]];
            let g = gap(q);
            if norm(g) < norm(f) || t < 1e-6 {
                (p, f) = (q, g);
                break;
            }
            t *= 0.5;
        }
    }
    (norm(f) <= 1e-8 * scale).then_some(p)
}

/// The intersection of fillets `s` and `l` (common sides `cs` and `cl`) from `far`, where
/// it meets `s`'s far side (`s`'s v, `l`'s u and v), back across `s` to `near`, where their
/// contacts meet on the common face: per step `s`'s fraction across, its v, and `l`'s u and
/// v, from the common face.
fn march(s: &Stripe, cs: usize, l: &Stripe, far: [f64; 3], near: [f64; 3]) -> Option<Vec<[f64; 4]>> {
    let n = 2 * ACROSS;
    let mut out = vec![[1.0, far[0], far[1], far[2]]];
    let mut p = far;
    let mut delta = [0.0; 3];
    for k in (1..n).rev() {
        let w = k as f64 / n as f64;
        let guess = [p[0] + delta[0], p[1] + delta[1], p[2] + delta[2]];
        let q = meet(s, across(cs, w), l, guess)?;
        delta = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        p = q;
        out.push([w, p[0], p[1], p[2]]);
    }
    out.push([0.0, near[0], near[1], near[2]]);
    out.reverse();
    Some(out)
}

/// Where fillet `s` (common side `cs`, ending at its edge's end `end`) reaches its far side
/// on fillet `l` (common side `cl`): `s`'s v, `l`'s u and v, if within `l`'s width.
fn far_meeting(s: &Stripe, cs: usize, end: usize, l: &Stripe, cl: usize, v_near: f64) -> Option<[f64; 3]> {
    let vs = if end == 0 { s.rolling.t0 } else { s.rolling.t1 };
    let p = s.surface.point(across(cs, 1.0), vs);
    let (u, v, _) = crate::geometry::extrema::newton_on_surface(p, &l.surface, (across(cl, 0.5), v_near));
    let found = meet(s, across(cs, 1.0), l, [vs, u, v])?;
    let fraction = if cl == 0 { found[1] } else { 1.0 - found[1] };
    (-1e-9..=1.0 + 1e-9).contains(&fraction).then_some(found)
}

/// The mitre where fillets `a` and `b` (fillet, end) meet at `vertex`: they share one face
/// (where their contacts cross), and their other two faces meet along an edge there.
fn mitre(model: &Model, shape: Shape, vertex: VertexId, stripes: &[Stripe], a: (usize, usize), b: (usize, usize)) -> Result<(Mitre, PendingMitre, PendingMitre)> {
    let (sa, sb) = (&stripes[a.0], &stripes[b.0]);
    if sa.rolling.convex != sb.rolling.convex {
        return Err(Error::InvalidInput("two fillets of opposite kinds meeting at a vertex are not supported yet".into()));
    }
    // The common face: side `ia` of `a`, `ib` of `b`.
    let common = (0..2).flat_map(|i| (0..2).map(move |j| (i, j))).filter(|&(i, j)| sa.rolling.sides[i].face.node == sb.rolling.sides[j].face.node).collect::<Vec<_>>();
    let [(ia, ib)] = common[..] else {
        return Err(Error::InvalidInput("two fillets meeting at a vertex without one common face".into()));
    };
    // The edge between their far faces at the vertex.
    let v_node: crate::topology::NodeId = vertex.into();
    let at_vertex = |e: Shape| {
        let (p, q) = model.edge_vertices(e);
        p.node == v_node || q.node == v_node
    };
    let edges_of = |face: Shape| -> Vec<EdgeId> { model.explore(face, NodeKind::Edge).filter(|&e| at_vertex(e)).filter_map(|e| EdgeId::try_from(e.node).ok()).collect() };
    let (on_a, on_b) = (edges_of(sa.rolling.sides[1 - ia].face), edges_of(sb.rolling.sides[1 - ib].face));
    let shared: Vec<EdgeId> = on_a.iter().copied().filter(|e| on_b.contains(e) && *e != sa.rolling.id && *e != sb.rolling.id).collect();
    let [third] = shared[..] else {
        return Err(Error::InvalidInput("two fillets meeting at a vertex whose other faces do not meet along an edge there".into()));
    };
    let faces_at = model.explore(shape, NodeKind::Face).filter(|&f| model.explore(f, NodeKind::Vertex).any(|v| v.node == v_node)).map(|f| f.node).collect::<std::collections::HashSet<_>>();
    if faces_at.len() != 3 {
        return Err(Error::InvalidInput("two fillets meeting at a vertex of other than three faces are not supported yet".into()));
    }
    let (va1, vb1, p1) = crossing(&sa.rolling.contacts[ia], &sa.ts, a.1, &sb.rolling.contacts[ib], &sb.ts, b.1, model.point(Shape::from(vertex)))?;
    // Across the one that gets to its far face first (`s`), `l` crossing that face after.
    let by_a = far_meeting(sa, ia, a.1, sb, ib, vb1).and_then(|far| march(sa, ia, sb, far, [va1, across(ib, 0.0), vb1]));
    let (s_is_a, path) = match by_a {
        Some(path) => (true, path),
        None => (
            false,
            far_meeting(sb, ib, b.1, sa, ia, va1)
                .and_then(|far| march(sb, ib, sa, far, [vb1, across(ia, 0.0), va1]))
                .ok_or_else(|| Error::InvalidInput("no intersection of two fillets meeting at a vertex".into()))?,
        ),
    };
    let ((s, cs), (l, cl)) = if s_is_a { ((sa, ia), (sb, ib)) } else { ((sb, ib), (sa, ia)) };
    let params: Vec<f64> = path.iter().map(|p| p[0]).collect();
    let points: Vec<Point> = path.iter().map(|p| {
        let (ps, pl) = (s.surface.point(across(cs, p[0]), p[1]), l.surface.point(p[2], p[3]));
        ps + (pl - ps) * 0.5
    }).collect();
    let curve = fit::interpolate_with_parameters(&points, &params)?;
    let on_s = fit::interpolate_with_parameters(&path.iter().map(|p| Point2::new(across(cs, p[0]), p[1])).collect::<Vec<_>>(), &params)?;
    let on_l = fit::interpolate_with_parameters(&path.iter().map(|p| Point2::new(p[2], p[3])).collect::<Vec<_>>(), &params)?;
    let mut deviation: f64 = 0.0;
    for w in params.windows(2) {
        let m = 0.5 * (w[0] + w[1]);
        let p = curve.point(m);
        let (qs, ql) = (on_s.point(m), on_l.point(m));
        deviation = deviation.max(s.surface.point(qs.x, qs.y).distance(p)).max(l.surface.point(ql.x, ql.y).distance(p));
    }
    let end = *path.last().expect("a path");
    let q = *points.last().expect("points");
    let (vs_far, l_end_fraction) = (end[1], if cl == 0 { end[2] } else { 1.0 - end[2] });
    // `l` on across `s`'s far face, to its own far side.
    let far_face = &s.rolling.sides[1 - cs];
    let size = 1.0 + (q - Point::ORIGIN).length();
    let mut cut = None;
    let mut cut_on_l = None;
    let mut last = (q, end[3]);
    if l_end_fraction < 1.0 - 1e-7 {
        let n = ACROSS;
        let mut v = end[3];
        let mut ws = Vec::with_capacity(n + 1);
        let mut pts = Vec::with_capacity(n + 1);
        let mut uv_l = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let frac = l_end_fraction + (1.0 - l_end_fraction) * k as f64 / n as f64;
            let u = across(cl, frac);
            if k > 0 {
                for _ in 0..60 {
                    let g = far_face.distance(l.surface.point(u, v))?.0;
                    let h = 1e-7;
                    let dg = (far_face.distance(l.surface.point(u, v + h))?.0 - g) / h;
                    if dg.abs() <= 1e-300 {
                        break;
                    }
                    v -= g / dg;
                    if g.abs() <= 1e-12 * size {
                        break;
                    }
                }
                if far_face.distance(l.surface.point(u, v))?.0.abs() > 1e-7 * size {
                    return Err(Error::InvalidInput("a fillet at a mitre does not cross the other's far face".into()));
                }
            }
            ws.push(k as f64 / n as f64);
            pts.push(if k == 0 { q } else { l.surface.point(u, v) });
            uv_l.push(Point2::new(u, v));
        }
        let curve = fit::interpolate_with_parameters(&pts, &ws)?;
        let on_l = fit::interpolate_with_parameters(&uv_l, &ws)?;
        let face_uv: Vec<Point2> = pts.iter().map(|&p| far_face.uv(p)).collect::<Result<_>>()?;
        let on_face = match project::analytic_pcurve(&Curve::BSpline(curve.clone()), &far_face.surface) {
            Some(p) => p,
            None => Curve2::BSpline(fit::interpolate_with_parameters(&unwrapped(face_uv, &far_face.surface), &ws)?),
        };
        for w in ws.windows(2) {
            let m = 0.5 * (w[0] + w[1]);
            let p = curve.point(m);
            let ql = on_l.point(m);
            deviation = deviation.max(l.surface.point(ql.x, ql.y).distance(p)).max(far_face.distance(p)?.0.abs());
        }
        last = (*pts.last().expect("points"), v);
        cut = Some(MitreCut { curve, face: far_face.face, on_face });
        cut_on_l = Some(on_l);
    }
    // The last point on the edge between the far faces.
    let (curve3, first, last_t) = model.edge_curve(Shape::from(third)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a mitre".into()))?;
    let on_edge = crate::geometry::extrema::nearest_on_curve(last.0, &curve3, first, last_t);
    if on_edge.distance > (1e-6 * size).max(1e3 * s.rolling.target) {
        return Err(Error::InvalidInput("two fillets at a vertex whose mitre misses the edge between their far faces".into()));
    }
    let mut mitre_points = vec![p1, q];
    if cut.is_some() {
        mitre_points.push(last.0);
    }
    let (vs1, vl1) = if s_is_a { (va1, vb1) } else { (vb1, va1) };
    // Each fillet's stop, from its first side to its second.
    let stop = |c: usize, v_common: f64, v_far: f64, far_at: usize, mut pieces: Vec<(usize, Orientation, BSplineCurve2)>| {
        let mut v = [0.0; 2];
        let mut at = [0; 2];
        (v[c], at[c]) = (v_common, 0);
        (v[1 - c], at[1 - c]) = (v_far, far_at);
        // (The pieces run from the common side; the other way round where that is the second.)
        if c == 1 {
            pieces.reverse();
            for piece in &mut pieces {
                piece.1 = piece.1.reversed();
            }
        }
        PendingMitre { v, at, pieces }
    };
    let ps = stop(cs, vs1, vs_far, 1, vec![(0, Orientation::Forward, on_s)]);
    let mut l_pieces = vec![(0, Orientation::Forward, on_l)];
    if let Some(on) = cut_on_l {
        l_pieces.push((1, Orientation::Forward, on));
    }
    let pl = stop(cl, vl1, last.1, mitre_points.len() - 1, l_pieces);
    let last = mitre_points.len() - 1;
    let m = Mitre { points: mitre_points, curve, cut, edges: vec![(third, vertex, on_edge.parameter, last)], deviation };
    Ok(if s_is_a { (m, ps, pl) } else { (m, pl, ps) })
}

/// Where fillets `a` and `b` (fillet, end) meet at `vertex` along edges going on into one
/// another (tangent there), sharing one face whose contact curves join there, and their
/// other two faces meeting along an edge at the vertex: the fillets joined across, from the
/// common face's contact at the vertex to where the far contacts cross that edge, straight
/// in each fillet's parameters (a mitre whose curve is that). `None` where the edges are not
/// tangent there.
fn joint(model: &Model, shape: Shape, vertex: VertexId, stripes: &[Stripe], a: (usize, usize), b: (usize, usize)) -> Option<Result<(Mitre, PendingMitre, PendingMitre)>> {
    let (sa, sb) = (&stripes[a.0], &stripes[b.0]);
    // The edges' directions away from the vertex: opposite where one goes on into the other.
    let away = |s: &Stripe, k: usize| -> Option<crate::math::Vector> {
        let (curve, first, last) = model.edge_curve(Shape::from(s.rolling.id))?;
        let (_, d) = curve.d1(if k == 0 { first } else { last });
        Some(if k == 0 { d } else { -d }.normalized())
    };
    let (da, db) = (away(sa, a.1)?, away(sb, b.1)?);
    if da.dot(db) > -1.0 + 1e-4 {
        return None;
    }
    Some((|| {
        if sa.rolling.convex != sb.rolling.convex {
            return Err(Error::InvalidInput("two fillets of opposite kinds meeting at a vertex are not supported yet".into()));
        }
        let common = (0..2).flat_map(|i| (0..2).map(move |j| (i, j))).filter(|&(i, j)| sa.rolling.sides[i].face.node == sb.rolling.sides[j].face.node).collect::<Vec<_>>();
        let at_end = |s: &Stripe, k: usize| if k == 0 { s.rolling.t0 } else { s.rolling.t1 };
        // (Both faces in common, or each tangent there to one of the other's: one edge going
        // on between the same faces, or faces going on smoothly, as a box's edge into the
        // mitre where fillets of its two sides meet; the fillets joined across their
        // sections at the vertex.)
        let (va1, vb1) = (at_end(sa, a.1), at_end(sb, b.1));
        let tangent_sides = || {
            let p = model.point(Shape::from(vertex));
            // (From a search of the whole surface, not from the last foot found, wherever the
            // sections' solving left it: a free-form face's foot found from afar may be another.)
            let normal = |s: &Stripe, i: usize| {
                let side = &s.rolling.sides[i];
                let saved = side.last.replace(None);
                let found = side.distance(p).ok().map(|d| d.2);
                side.last.set(saved);
                found
            };
            let tangent = |i: usize, j: usize| matches!((normal(sa, i), normal(sb, j)), (Some(x), Some(y)) if x.dot(y) >= 1.0 - crate::tolerance::ANGULAR.sqrt());
            if tangent(0, 0) && tangent(1, 1) {
                Some(0)
            } else if tangent(0, 1) && tangent(1, 0) {
                Some(1)
            } else {
                None
            }
        };
        let paired = match common.len() {
            2 => Some(common.iter().find(|&&(i, _)| i == 0).expect("side 0 in common").1),
            0 => tangent_sides(),
            _ => None,
        };
        if let Some(jb) = paired {
            // (Both faces in common, with another edge of one of them through the vertex, a
            // sphere's seam: the fillets joined where their contacts on that face cross it,
            // for the contacts' ends to lie on it, not across their sections through the
            // vertex, whose contacts there miss it a little.)
            let (va1, vb1) = if common.len() == 2 { seam_joint(model, vertex, sa, sb, jb, va1, vb1).unwrap_or((va1, vb1)) } else { (va1, vb1) };
            let ends: Vec<Point> = (0..2).map(|i| {
                let (pa, pb) = (sa.rolling.contacts[i].point(va1), sb.rolling.contacts[if i == 0 { jb } else { 1 - jb }].point(vb1));
                pa + (pb - pa) * 0.5
            }).collect();
            let n = ACROSS;
            let ws: Vec<f64> = (0..=n).map(|k| k as f64 / n as f64).collect();
            let (uv_a, uv_b): (Vec<Point2>, Vec<Point2>) = ws.iter().map(|&w| (Point2::new(w, va1), Point2::new(across(jb, w), vb1))).unzip();
            let points: Vec<Point> = uv_a.iter().zip(&uv_b).map(|(p, q)| {
                let (x, y) = (sa.surface.point(p.x, p.y), sb.surface.point(q.x, q.y));
                x + (y - x) * 0.5
            }).collect();
            let curve = fit::interpolate_with_parameters(&points, &ws)?;
            let mut deviation: f64 = 0.0;
            for (k, &w) in ws.iter().enumerate() {
                let p = curve.point(w);
                deviation = deviation.max(sa.surface.point(uv_a[k].x, uv_a[k].y).distance(p)).max(sb.surface.point(uv_b[k].x, uv_b[k].y).distance(p));
            }
            for (i, e) in ends.iter().enumerate() {
                deviation = deviation.max(sa.rolling.contacts[i].point(va1).distance(*e));
            }
            let on_a = BSplineCurve2::new(vec![uv_a[0], uv_a[n]], vec![0.0, 1.0], vec![2, 2], 1);
            let on_b = BSplineCurve2::new(vec![uv_b[0], uv_b[n]], vec![0.0, 1.0], vec![2, 2], 1);
            let stop = |c: usize, v_at: f64, on: BSplineCurve2| {
                let mut at = [0; 2];
                (at[c], at[1 - c]) = (0, 1);
                let mut pieces = vec![(0, Orientation::Forward, on)];
                if c == 1 {
                    pieces[0].1 = pieces[0].1.reversed();
                }
                PendingMitre { v: [v_at, v_at], at, pieces }
            };
            // (Faces tangent pairwise: the edge between each pair, through the vertex, ends
            // where the fillets' contacts meet it.)
            let mut edges = Vec::new();
            if common.is_empty() {
                let v_node: crate::topology::NodeId = vertex.into();
                for (i, end) in ends.iter().enumerate() {
                    let (fa, fb) = (sa.rolling.sides[i].face, sb.rolling.sides[if i == 0 { jb } else { 1 - jb }].face);
                    let at = |face: Shape| -> Vec<EdgeId> {
                        model.explore(face, NodeKind::Edge).filter(|&e| { let (p, q) = model.edge_vertices(e); p.node == v_node || q.node == v_node }).filter_map(|e| EdgeId::try_from(e.node).ok()).collect()
                    };
                    let between: Vec<EdgeId> = at(fa).into_iter().filter(|e| at(fb).contains(e)).collect();
                    let [edge] = between[..] else {
                        return Err(Error::InvalidInput("two fillets joining at a vertex whose tangent faces do not meet along an edge there".into()));
                    };
                    let (curve, first, last) = model.edge_curve(Shape::from(edge)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a joint".into()))?;
                    let on = crate::geometry::extrema::nearest_on_curve(*end, &curve, first, last);
                    deviation = deviation.max(on.distance);
                    edges.push((edge, vertex, on.parameter, i));
                }
            } else {
                // (Both faces in common: any other edge of theirs through the vertex, as a
                // sphere's seam reaching the edge, ends where the fillet's contact crosses it.)
                let v_node: crate::topology::NodeId = vertex.into();
                for (i, end) in ends.iter().enumerate() {
                    let others: Vec<EdgeId> = model
                        .explore(sa.rolling.sides[i].face, NodeKind::Edge)
                        .filter(|&e| { let (p, q) = model.edge_vertices(e); p.node == v_node || q.node == v_node })
                        .filter_map(|e| EdgeId::try_from(e.node).ok())
                        .filter(|&e| e != sa.rolling.id && e != sb.rolling.id)
                        .collect();
                    for edge in others {
                        if edges.iter().any(|x: &(EdgeId, VertexId, f64, usize)| x.0 == edge) {
                            continue;
                        }
                        let (curve, first, last) = model.edge_curve(Shape::from(edge)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a joint".into()))?;
                        let on = crate::geometry::extrema::nearest_on_curve(*end, &curve, first, last);
                        deviation = deviation.max(on.distance);
                        edges.push((edge, vertex, on.parameter, i));
                    }
                }
            }
            let m = Mitre { points: ends, curve, cut: None, edges, deviation };
            return Ok((m, stop(0, va1, on_a), stop(jb, vb1, on_b)));
        }
        let [(ia, ib)] = common[..] else {
            return Err(Error::InvalidInput("two fillets joining at a vertex without one common face".into()));
        };
        let v_node: crate::topology::NodeId = vertex.into();
        let at_vertex = |e: Shape| {
            let (p, q) = model.edge_vertices(e);
            p.node == v_node || q.node == v_node
        };
        let edges_of = |face: Shape| -> Vec<EdgeId> { model.explore(face, NodeKind::Edge).filter(|&e| at_vertex(e)).filter_map(|e| EdgeId::try_from(e.node).ok()).collect() };
        let (on_a, on_b) = (edges_of(sa.rolling.sides[1 - ia].face), edges_of(sb.rolling.sides[1 - ib].face));
        let shared: Vec<EdgeId> = on_a.iter().copied().filter(|e| on_b.contains(e) && *e != sa.rolling.id && *e != sb.rolling.id).collect();
        let [third] = shared[..] else {
            return Err(Error::InvalidInput("two fillets joining at a vertex whose other faces do not meet along an edge there".into()));
        };
        let (curve3, first3, last3) = model.edge_curve(Shape::from(third)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a joint".into()))?;
        // Each fillet's spine parameter at the vertex, and where its far contact comes
        // nearest the edge between the far faces (from the middle of the fillet to that end).
        let far_crossing = |s: &Stripe, c: usize, k: usize| -> (f64, f64, Point) {
            let contact = &s.rolling.contacts[1 - c];
            let (mid, end) = (0.5 * (s.rolling.t0 + s.rolling.t1), at_end(s, k));
            let (lo, hi) = (s.ts[0], s.ts[s.ts.len() - 1]);
            // (Some way past the vertex too: the far contact may cross the edge beyond it.)
            let far = (end + (end - mid)).clamp(lo.min(hi), lo.max(hi));
            let gap = |w: f64| crate::geometry::extrema::nearest_on_curve(contact.point(w), &curve3, first3, last3).distance;
            let n = 256;
            let (mut best, mut at) = (f64::INFINITY, mid);
            for j in 0..=n {
                let w = mid + (far - mid) * j as f64 / n as f64;
                let d = gap(w);
                if d < best {
                    (best, at) = (d, w);
                }
            }
            let step = (far - mid).abs() / n as f64;
            let (mut x, mut y) = (at - step, at + step);
            for _ in 0..100 {
                let (p, q) = (x + (y - x) / 3.0, y - (y - x) / 3.0);
                if gap(p) < gap(q) { y = q } else { x = p }
            }
            let w = 0.5 * (x + y);
            (w, gap(w), contact.point(w))
        };
        let (wa, gap_a, qa) = far_crossing(sa, ia, a.1);
        let (wb, gap_b, qb) = far_crossing(sb, ib, b.1);
        // (Where the edge slants across the joint, one far contact crosses it before its
        // fillet's end, the other only past its own, on its face carried on past its patch,
        // which strays from the first's face there: that one's nearest approach to the first
        // crossing then.)
        let own = |s: &Stripe, w: f64| w >= s.rolling.t0.min(s.rolling.t1) - 1e-9 && w <= s.rolling.t0.max(s.rolling.t1) + 1e-9;
        let nearest_to = |s: &Stripe, c: usize, k: usize, q: Point| -> (f64, f64, Point) {
            let contact = &s.rolling.contacts[1 - c];
            let (lo, hi) = (s.ts[0].min(s.ts[s.ts.len() - 1]), s.ts[0].max(s.ts[s.ts.len() - 1]));
            let end = at_end(s, k);
            let span = (s.rolling.t1 - s.rolling.t0).abs();
            let (x0, x1) = ((end - 0.25 * span).max(lo), (end + 0.25 * span).min(hi));
            let d = |w: f64| contact.point(w).distance(q);
            let n = 256;
            let at = (0..=n).map(|j| x0 + (x1 - x0) * j as f64 / n as f64).min_by(|x, y| d(*x).total_cmp(&d(*y))).unwrap_or(end);
            let step = (x1 - x0) / n as f64;
            let (mut x, mut y) = (at - step, at + step);
            for _ in 0..100 {
                let (p, r) = (x + (y - x) / 3.0, y - (y - x) / 3.0);
                if d(p) < d(r) { y = r } else { x = p }
            }
            let w = 0.5 * (x + y);
            let p = contact.point(w);
            (w, crate::geometry::extrema::nearest_on_curve(p, &curve3, first3, last3).distance, p)
        };
        // (That one need only come near the first's crossing, which the joint goes through.)
        let near_at = |q: Point| (1e-6 * (1.0 + (q - Point::ORIGIN).length())).max(1e3 * sa.rolling.target.max(sb.rolling.target)).max(2.0 * model.tolerance(Shape::from(third)));
        // (The far contacts, each square to its own spine, cross the edge a little apart:
        // within a twentieth of the fillets' width, the joint's tolerance covers that.)
        let width = sa.rolling.contacts[0].point(at_end(sa, a.1)).distance(sa.rolling.contacts[1].point(at_end(sa, a.1)));
        let meet = |gap_a: f64, gap_b: f64, qa: Point, qb: Point, q: Point| gap_a.max(gap_b) <= near_at(q) && qa.distance(qb) <= near_at(q).max(0.05 * width);
        let (((wa, _, qa), (wb, _, qb)), q) = if meet(gap_a, gap_b, qa, qb, qa + (qb - qa) * 0.5) {
            (((wa, gap_a, qa), (wb, gap_b, qb)), qa + (qb - qa) * 0.5)
        } else {
            let slant = match (own(sa, wa), own(sb, wb)) {
                (true, false) => Some((((wa, gap_a, qa), nearest_to(sb, ib, b.1, qa)), qa)),
                (false, true) => Some(((nearest_to(sa, ia, a.1, qb), (wb, gap_b, qb)), qb)),
                _ => None,
            };
            match slant {
                Some((((wa, ga, qa), (wb, gb, qb)), q)) if meet(ga.min(gb), ga.min(gb), qa, qb, q) => (((wa, ga, qa), (wb, gb, qb)), q),
                _ => return Err(Error::InvalidInput("two fillets joining at a vertex whose far contacts miss the edge between their far faces".into())),
            }
        };
        let on_edge = crate::geometry::extrema::nearest_on_curve(q, &curve3, first3, last3);
        // (Where the far faces are tangent along their edge (a side going on into a rounded
        // corner), fillets' balls (not chamfers') where their far contacts cross it are one, touching
        // the common face and that point: the fillets meet across its section, each running
        // on to it.)
        let balls = !(sa.rolling.sections.iter().any(|x| x.straight) || sb.rolling.sections.iter().any(|x| x.straight));
        let smooth = balls && matches!(
            (sa.rolling.sides[1 - ia].distance(q).map(|d| d.2), sb.rolling.sides[1 - ib].distance(q).map(|d| d.2)),
            (Ok(x), Ok(y)) if x.dot(y) >= 1.0 - crate::tolerance::ANGULAR.sqrt().max(crate::blend::KINK)
        );
        let (va1, vb1) = if smooth { (wa, wb) } else { (va1, vb1) };
        // The common face's point: both contacts at the vertex.
        let p1 = {
            let (pa, pb) = (sa.rolling.contacts[ia].point(va1), sb.rolling.contacts[ib].point(vb1));
            pa + (pb - pa) * 0.5
        };
        // Across each, straight in its parameters from the common side to the far one.
        let n = ACROSS;
        let ws: Vec<f64> = (0..=n).map(|k| k as f64 / n as f64).collect();
        let line = |c: usize, v0: f64, v1: f64| -> Vec<Point2> { ws.iter().map(|&w| Point2::new(across(c, w), v0 + (v1 - v0) * w)).collect() };
        let (uv_a, uv_b) = (line(ia, va1, wa), line(ib, vb1, wb));
        let points: Vec<Point> = uv_a.iter().zip(&uv_b).map(|(p, q)| {
            let (x, y) = (sa.surface.point(p.x, p.y), sb.surface.point(q.x, q.y));
            x + (y - x) * 0.5
        }).collect();
        let curve = fit::interpolate_with_parameters(&points, &ws)?;
        let on_a_fillet = BSplineCurve2::new(vec![uv_a[0], uv_a[n]], vec![0.0, 1.0], vec![2, 2], 1);
        let on_b_fillet = BSplineCurve2::new(vec![uv_b[0], uv_b[n]], vec![0.0, 1.0], vec![2, 2], 1);
        let mut deviation: f64 = p1.distance(sa.rolling.contacts[ia].point(va1)).max(q.distance(qa)).max(q.distance(qb));
        for k in 0..=4 * n {
            let w = k as f64 / (4 * n) as f64;
            let p = curve.point(w);
            let (x, y) = (on_a_fillet.point(w), on_b_fillet.point(w));
            deviation = deviation.max(sa.surface.point(x.x, x.y).distance(p)).max(sb.surface.point(y.x, y.y).distance(p));
        }
        let stop = |c: usize, v_common: f64, v_far: f64, on: BSplineCurve2| {
            let mut v = [0.0; 2];
            let mut at = [0; 2];
            (v[c], at[c]) = (v_common, 0);
            (v[1 - c], at[1 - c]) = (v_far, 1);
            let mut pieces = vec![(0, Orientation::Forward, on)];
            if c == 1 {
                pieces.reverse();
                for piece in &mut pieces {
                    piece.1 = piece.1.reversed();
                }
            }
            PendingMitre { v, at, pieces }
        };
        let pa = stop(ia, va1, wa, on_a_fillet);
        let pb = stop(ib, vb1, wb, on_b_fillet);
        let m = Mitre { points: vec![p1, q], curve, cut: None, edges: vec![(third, vertex, on_edge.parameter, 1)], deviation };
        Ok((m, pa, pb))
    })())
}

/// Where fillets `sa` and `sb` (side 0 of `sa` on side `jb` of `sb`), joined at `vertex`
/// between the same two faces, cross another edge of one of those faces through the vertex
/// (from their ends there, `va1` and `vb1`): each one's spine parameter where its contact
/// on that face meets the edge, if both do within a tenth of their lengths.
fn seam_joint(model: &Model, vertex: VertexId, sa: &Stripe, sb: &Stripe, jb: usize, va1: f64, vb1: f64) -> Option<(f64, f64)> {
    let v_node: crate::topology::NodeId = vertex.into();
    for i in 0..2 {
        let Some(edge) = model
            .explore(sa.rolling.sides[i].face, NodeKind::Edge)
            .filter(|&e| { let (p, q) = model.edge_vertices(e); p.node == v_node || q.node == v_node })
            .filter_map(|e| EdgeId::try_from(e.node).ok())
            .find(|&e| e != sa.rolling.id && e != sb.rolling.id)
        else {
            continue;
        };
        let (curve, first, last) = model.edge_curve(Shape::from(edge))?;
        // Along the contact, from `t` a tenth of the fillet's length either way: where it
        // comes nearest the edge (sampled, then by thirds).
        let crossing = |s: &Stripe, c: usize, t: f64| -> (f64, f64) {
            let contact = &s.rolling.contacts[c];
            let (lo, hi) = (s.ts[0].min(s.ts[s.ts.len() - 1]), s.ts[0].max(s.ts[s.ts.len() - 1]));
            let reach = 0.1 * (s.rolling.t1 - s.rolling.t0).abs();
            let (x0, x1) = ((t - reach).max(lo), (t + reach).min(hi));
            let gap = |w: f64| crate::geometry::extrema::nearest_on_curve(contact.point(w), &curve, first, last).distance;
            let n = 64;
            let (mut best, mut at) = (f64::INFINITY, t);
            for j in 0..=n {
                let w = x0 + (x1 - x0) * j as f64 / n as f64;
                let d = gap(w);
                if d < best {
                    (best, at) = (d, w);
                }
            }
            let step = (x1 - x0) / n as f64;
            let (mut x, mut y) = ((at - step).max(x0), (at + step).min(x1));
            for _ in 0..100 {
                let (p, q) = (x + (y - x) / 3.0, y - (y - x) / 3.0);
                if gap(p) < gap(q) { y = q } else { x = p }
            }
            let w = 0.5 * (x + y);
            (w, gap(w))
        };
        let (wa, gap_a) = crossing(sa, i, va1);
        let (wb, gap_b) = crossing(sb, if i == 0 { jb } else { 1 - jb }, vb1);
        let size = 1.0 + (model.point(Shape::from(vertex)) - Point::ORIGIN).length();
        let near = (1e-9 * size).max(10.0 * sa.rolling.target.max(sb.rolling.target));
        return (gap_a.max(gap_b) <= near).then_some((wa, wb));
    }
    None
}

/// Fillet `s`'s end `k` (0 at its edge's start), where from its vertex each face goes on
/// only along a free edge (one no other face uses) and no face spans the two: each contact
/// curve stops where it crosses that edge, as OCCT's fillets end at a free boundary.
fn free_end(model: &Model, shape: Shape, s: &Stripe, k: usize) -> Result<Option<End>> {
    let vertex = if k == 0 { s.rolling.va } else { s.rolling.vb };
    let at_vertex = |e: EdgeId| {
        let (a, b) = model.edge_vertices(Shape::from(e));
        a.node == vertex.node || b.node == vertex.node
    };
    let uses = |e: EdgeId| model.explore(shape, NodeKind::Face).filter(|&f| model.explore(f, NodeKind::Edge).any(|x| x.node == crate::topology::NodeId::from(e))).count();
    let mut edges = [(EdgeId::INVALID, 0.0); 2];
    let mut v = [0.0; 2];
    for (i, side) in s.rolling.sides.iter().enumerate() {
        let others: Vec<EdgeId> = model.explore(side.face, NodeKind::Edge).filter_map(|e| EdgeId::try_from(e.node).ok()).filter(|&e| e != s.rolling.id && at_vertex(e)).collect();
        let [other] = others[..] else { return Ok(None) };
        if uses(other) != 1 {
            return Ok(None);
        }
        let (curve, first, last) = model.edge_curve(Shape::from(other)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        // Along the half of the contact curve towards this end (past the edge included).
        let (lo, hi) = (s.ts[0], s.ts[s.ts.len() - 1]);
        let mid = 0.5 * (lo + hi);
        let (a, b) = if k == 0 { (lo, mid) } else { (mid, hi) };
        let gap = |w: f64| crate::geometry::extrema::nearest_on_curve(s.rolling.contacts[i].point(w), &curve, first, last).distance;
        let n = 128;
        let (mut best, mut at) = (f64::INFINITY, a);
        for j in 0..=n {
            let w = a + (b - a) * j as f64 / n as f64;
            let d = gap(w);
            if d < best {
                (best, at) = (d, w);
            }
        }
        let step = (b - a) / n as f64;
        let (mut x, mut y) = ((at - step).max(a), (at + step).min(b));
        for _ in 0..100 {
            let (p, q) = (x + (y - x) / 3.0, y - (y - x) / 3.0);
            if gap(p) < gap(q) { y = q } else { x = p }
        }
        let w = 0.5 * (x + y);
        let point = s.rolling.contacts[i].point(w);
        let on = crate::geometry::extrema::nearest_on_curve(point, &curve, first, last);
        // (Within the free edge's tolerance, too: its curve strays from its faces that much.)
        let edge_tolerance = model.tolerance(Shape::from(other));
        if on.distance > (1e3 * s.rolling.target).max(1e-6 * (1.0 + (point - Point::ORIGIN).length())).max(2.0 * edge_tolerance) {
            return Err(Error::InvalidInput("a fillet's contact curve does not reach the free edge at its end".into()));
        }
        v[i] = w;
        edges[i] = (other, on.parameter);
    }
    // (No face joins the two edges: else it is an end face.)
    let spanned = model.explore(shape, NodeKind::Face).any(|f| {
        let es: Vec<_> = model.explore(f, NodeKind::Edge).map(|e| e.node).collect();
        es.contains(&crate::topology::NodeId::from(edges[0].0)) && es.contains(&crate::topology::NodeId::from(edges[1].0))
    });
    if spanned {
        return Ok(None);
    }
    Ok(Some(End::Free { vertex: VertexId::try_from(vertex.node).expect("a vertex"), v, edges }))
}

/// A fillet's end at a corner, before the corner is numbered.
struct PendingEnd {
    at: [usize; 2],
    v: [f64; 2],
}

impl PendingEnd {
    fn take_end(&self, corner: usize) -> End {
        End::Corner { corner, at: self.at, v: self.v }
    }
}

/// The corner where three fillets (`meeting`: fillet and end) meet at `vertex`.
fn corner(model: &Model, shape: Shape, vertex: VertexId, stripes: &[Stripe], meeting: &[(usize, usize)]) -> Result<(Corner, Vec<PendingEnd>)> {
    let _ = shape;
    // The faces: each shared by two of the fillets.
    let mut faces: Vec<Shape> = Vec::new();
    for &(i, _) in meeting {
        for side in &stripes[i].rolling.sides {
            if !faces.iter().any(|f| f.node == side.face.node) {
                faces.push(side.face);
            }
        }
    }
    if faces.len() != 3 {
        return Err(Error::InvalidInput("a corner of fillets where other than three faces meet".into()));
    }
    let side_on = |i: usize, face: Shape| stripes[i].rolling.sides.iter().position(|s| s.face.node == face.node);
    let mut points = Vec::with_capacity(3);
    let mut v: Vec<[f64; 2]> = vec![[0.0; 2]; meeting.len()];
    let mut at: Vec<[usize; 2]> = vec![[0; 2]; meeting.len()];
    for (f, &face) in faces.iter().enumerate() {
        let on: Vec<(usize, usize)> = meeting.iter().enumerate().filter_map(|(n, &(i, _))| side_on(i, face).map(|s| (n, s))).collect();
        let [(n1, s1), (n2, s2)] = on[..] else {
            return Err(Error::InvalidInput("a corner face not shared by two of its fillets".into()));
        };
        let (a, b) = (&stripes[meeting[n1].0], &stripes[meeting[n2].0]);
        let (t1, t2, p) = crossing(&a.rolling.contacts[s1], &a.ts, meeting[n1].1, &b.rolling.contacts[s2], &b.ts, meeting[n2].1, model.point(Shape::from(vertex)))?;
        let _ = model;
        points.push((face, p));
        (v[n1][s1], at[n1][s1]) = (t1, f);
        (v[n2][s2], at[n2][s2]) = (t2, f);
    }
    let pending = (0..meeting.len()).map(|n| PendingEnd { at: at[n], v: v[n] }).collect();
    Ok((Corner { points, stripes: meeting.to_vec() }, pending))
}

/// Where two contact curves on one face cross near the ends `end_a`, `end_b` of their
/// fillets (0 the start of the edge, 1 its end): their parameters and the point. Of several
/// crossings, the one nearest `corner` (the fillets' vertex).
fn crossing(a: &BSplineCurve, ts_a: &[f64], end_a: usize, b: &BSplineCurve, ts_b: &[f64], end_b: usize, corner: Point) -> Result<(f64, f64, Point)> {
    // Over the half of each curve towards the corner (its ends past the edge included).
    let half = |ts: &[f64], end: usize| {
        let (lo, hi) = (ts[0], ts[ts.len() - 1]);
        let mid = 0.5 * (lo + hi);
        if end == 0 { (lo, mid) } else { (mid, hi) }
    };
    let ((a0, a1), (b0, b1)) = (half(ts_a, end_a), half(ts_b, end_b));
    let n = 64;
    let at = |i: usize, j: usize| (a0 + (a1 - a0) * i as f64 / n as f64, b0 + (b1 - b0) * j as f64 / n as f64);
    let grid: Vec<Vec<f64>> = (0..=n).map(|i| (0..=n).map(|j| { let (ta, tb) = at(i, j); a.point(ta).distance(b.point(tb)) }).collect()).collect();
    // Newton on the gap's square, from a sample.
    let refine = |(mut ta, mut tb): (f64, f64)| -> (f64, f64) {
        for _ in 0..50 {
            let (pa, da, dda) = a.d2(ta);
            let (pb, db, ddb) = b.d2(tb);
            let r = pa - pb;
            let g = [r.dot(da), -r.dot(db)];
            let h = [[da.dot(da) + r.dot(dda), -da.dot(db)], [-da.dot(db), db.dot(db) - r.dot(ddb)]];
            let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
            if det.abs() <= 1e-300 {
                break;
            }
            let (dx, dy) = ((g[0] * h[1][1] - g[1] * h[0][1]) / det, (h[0][0] * g[1] - h[1][0] * g[0]) / det);
            ta = (ta - dx).clamp(a0.min(a1), a0.max(a1));
            tb = (tb - dy).clamp(b0.min(b1), b0.max(b1));
            if dx.abs() + dy.abs() <= 1e-14 * (1.0 + ta.abs() + tb.abs()) {
                break;
            }
        }
        (ta, tb)
    };
    // (Fitted through sections, two contact curves on one face may pass a few millionths
    // apart where they cross.) From each of the samples' local minima, the crossings found.
    let close = |ta: f64, tb: f64| {
        let (pa, pb) = (a.point(ta), b.point(tb));
        pa.distance(pb) <= 1e-5 * (1.0 + (pa - Point::ORIGIN).length())
    };
    let mut found: Vec<(f64, f64)> = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            let d = grid[i][j];
            let local = (i.saturating_sub(1)..=(i + 1).min(n)).all(|x| (j.saturating_sub(1)..=(j + 1).min(n)).all(|y| grid[x][y] >= d));
            if local {
                let (ta, tb) = refine(at(i, j));
                if close(ta, tb) {
                    found.push((ta, tb));
                }
            }
        }
    }
    let Some(&(ta, tb)) = found.iter().min_by(|x, y| a.point(x.0).distance(corner).total_cmp(&a.point(y.0).distance(corner))) else {
        return Err(Error::InvalidInput("fillets at a corner whose contact curves do not cross".into()));
    };
    let (pa, pb) = (a.point(ta), b.point(tb));
    Ok((ta, tb, pa + (pb - pa) * 0.5))
}

/// The curve across fillet `stripe` from its first side at `v0` to its second at `v1`,
/// straight in its parameters, over [0, 1]: in 3D, and on the fillet.
fn stop_curve(stripe: &Stripe, v0: f64, v1: f64) -> Result<(BSplineCurve, BSplineCurve2)> {
    let on_fillet = BSplineCurve2::new(vec![Point2::new(0.0, v0), Point2::new(1.0, v1)], vec![0.0, 1.0], vec![2, 2], 1);
    let params: Vec<f64> = (0..=ACROSS).map(|k| k as f64 / ACROSS as f64).collect();
    let points: Vec<Point> = params.iter().map(|&w| stripe.surface.point(w, v0 + (v1 - v0) * w)).collect();
    Ok((fit::interpolate_with_parameters(&points, &params)?, on_fillet))
}

/// Builds the rounded solid.
fn build(model: &mut Model, shape: Shape, stripes: &[Stripe], ends: &[[End; 2]], corners: &[Corner], mitres: &[Mitre]) -> Result<Shape> {
    // Within twice the contact curves' distance from their faces (between the fillet's
    // ends: past them the contacts run on through sections past the edge's ends, which may
    // stray from the faces' surfaces).
    let mut deviation: f64 = 0.0;
    let mut size: f64 = 1.0;
    let v_of = |end: &End| match end {
        End::Face(cut) => cut.v,
        End::Free { v, .. } | End::Mitre { v, .. } | End::Corner { v, .. } | End::Point { v, .. } => *v,
    };
    for (s, ends) in stripes.iter().zip(ends) {
        let (from, to) = (v_of(&ends[0]), v_of(&ends[1]));
        for i in 0..2 {
            let (a, b) = (from[i], to[i]);
            for k in 0..=100 {
                let p = s.rolling.contacts[i].point(a + (b - a) * k as f64 / 100.0);
                deviation = deviation.max(s.rolling.sides[i].distance(p)?.0.abs());
                size = size.max((p - Point::ORIGIN).length());
            }
        }
    }
    let deviation = mitres.iter().fold(deviation, |d, m| d.max(m.deviation));
    let tolerance = (deviation * 2.0).max(crate::tolerance::LINEAR).max(1e-8 * size);
    // The mitres' points as vertices, their curves as edges (each shared by two fillets).
    let mitre_vertices: Vec<Vec<VertexId>> = mitres.iter().map(|m| m.points.iter().map(|&p| model.add_vertex(p, tolerance)).collect()).collect();
    let mut mitre_edges: Vec<Vec<EdgeId>> = Vec::with_capacity(mitres.len());
    // (face, edge, pcurve, the edge's vertices) inserted into end faces
    let mut inserts: Vec<(Shape, (EdgeId, Curve2, [VertexId; 2]))> = Vec::new();
    for (m, vs) in mitres.iter().zip(&mitre_vertices) {
        let e = model.add_edge(vs[0], vs[1], tolerance)?;
        model.set_edge_curve(e.into(), Curve::BSpline(m.curve.clone()), tolerance)?;
        model.set_edge_range(e.into(), 0.0, 1.0)?;
        let mut edges = vec![e];
        if let Some(cut) = &m.cut {
            let e = model.add_edge(vs[1], vs[2], tolerance)?;
            model.set_edge_curve(e.into(), Curve::BSpline(cut.curve.clone()), tolerance)?;
            model.set_edge_range(e.into(), 0.0, 1.0)?;
            inserts.push((cut.face, (e, cut.on_face.clone(), [vs[1], vs[2]])));
            edges.push(e);
        }
        mitre_edges.push(edges);
    }
    // The corners' points as vertices.
    let corner_vertices: Vec<Vec<VertexId>> = corners.iter().map(|c| c.points.iter().map(|&(_, p)| model.add_vertex(p, tolerance)).collect()).collect();
    // Per fillet and end: the contacts' parameters there and their vertices; the curve
    // across (3D and on the fillet) and, for an end on a face, its pcurve there.
    struct Stop {
        v: [f64; 2],
        vertices: [VertexId; 2],
        curve: BSplineCurve,
        on_fillet: Curve2,
        /// At a mitre: its edges from the first side to the second, each with the coedge's
        /// orientation that way and its pcurve on the fillet.
        mitre: Option<Vec<(EdgeId, Orientation, Curve2)>>,
    }
    let mut stops: Vec<[Stop; 2]> = Vec::with_capacity(stripes.len());
    // (edge, its vertex at a fillet's end, the vertex replacing it, the parameter there)
    let mut trims: Vec<(EdgeId, VertexId, VertexId, f64)> = Vec::new();
    // (edges a fillet takes away whole, its contact ending at their far vertex)
    let mut removed: Vec<EdgeId> = Vec::new();
    let mut across_edges: Vec<[EdgeId; 2]> = Vec::with_capacity(stripes.len());
    for (i, s) in stripes.iter().enumerate() {
        let mut pair = Vec::with_capacity(2);
        for end in &ends[i] {
            pair.push(match end {
                End::Face(cut) => {
                    if cut.beyond.is_some() {
                        return Err(Error::InvalidInput("a fillet running over onto the next face is not supported yet in a network".into()));
                    }
                    if cut.fan.is_some() {
                        return Err(Error::InvalidInput("a fillet ending across several faces is not supported yet in a network".into()));
                    }
                    let vertices = [contact_vertex(model, Some(cut), 0, s.rolling.contacts[0].point(cut.v[0]), tolerance)?, contact_vertex(model, Some(cut), 1, s.rolling.contacts[1].point(cut.v[1]), tolerance)?];
                    for (j, &(old, t)) in cut.edges.iter().enumerate() {
                        if cut.onto[j].is_some() {
                            removed.push(old);
                        } else if cut.extended[j] {
                            let (e, on_side, on_end, vs) = extension(model, cut, j, &s.rolling.sides[j], vertices[j], tolerance)?;
                            inserts.push((s.rolling.sides[j].face, (e, on_side, vs)));
                            inserts.push((cut.face.face, (e, on_end, vs)));
                        } else {
                            trims.push((old, cut.vertex, vertices[j], t));
                        }
                    }
                    Stop { v: cut.v, vertices, curve: cut.curve.clone(), on_fillet: Curve2::BSpline(cut.on_fillet.clone()), mitre: None }
                }
                End::Free { vertex, v, edges } => {
                    let vertices = [0, 1].map(|j| model.add_vertex(s.rolling.contacts[j].point(v[j]), tolerance));
                    for (j, &(old, t)) in edges.iter().enumerate() {
                        trims.push((old, *vertex, vertices[j], t));
                    }
                    let (curve, on_fillet) = stop_curve(s, v[0], v[1])?;
                    Stop { v: *v, vertices, curve, on_fillet: Curve2::BSpline(on_fillet), mitre: None }
                }
                End::Mitre { mitre, v, at, pieces } => {
                    let vertices = [mitre_vertices[*mitre][at[0]], mitre_vertices[*mitre][at[1]]];
                    let pieces = pieces.iter().map(|(k, o, p)| (mitre_edges[*mitre][*k], *o, Curve2::BSpline(p.clone()))).collect();
                    Stop { v: *v, vertices, curve: mitres[*mitre].curve.clone(), on_fillet: Curve2::BSpline(BSplineCurve2::new(vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)], vec![0.0, 1.0], vec![2, 2], 1)), mitre: Some(pieces) }
                }
                End::Corner { corner, at, v } => {
                    let (curve, on_fillet) = stop_curve(s, v[0], v[1])?;
                    Stop { v: *v, vertices: [corner_vertices[*corner][at[0]], corner_vertices[*corner][at[1]]], curve, on_fillet: Curve2::BSpline(on_fillet), mitre: None }
                }
                End::Point { vertex, v } => {
                    let (curve, on_fillet) = stop_curve(s, v[0], v[1])?;
                    Stop { v: *v, vertices: [*vertex, *vertex], curve, on_fillet: Curve2::BSpline(on_fillet), mitre: None }
                }
            });
        }
        let [a, b]: [Stop; 2] = pair.try_into().map_err(|_| Error::Failed("a fillet's ends".into()))?;
        // The edges across at both ends.
        let mut edges = [EdgeId::INVALID; 2];
        for (k, stop) in [&a, &b].into_iter().enumerate() {
            if let Some(pieces) = &stop.mitre {
                edges[k] = pieces[0].0;
                continue;
            }
            let e = model.add_edge(stop.vertices[0], stop.vertices[1], tolerance)?;
            // (Where the fillet narrows to nothing: degenerate, without a curve.)
            if !matches!(ends[i][k], End::Point { .. }) {
                model.set_edge_curve(e.into(), Curve::BSpline(stop.curve.clone()), tolerance)?;
                model.set_edge_range(e.into(), 0.0, 1.0)?;
            }
            edges[k] = e;
            if let End::Face(cut) = &ends[i][k] {
                let pcurve = match project::analytic_pcurve(&Curve::BSpline(cut.curve.clone()), &cut.face.surface) {
                    Some(p) => p,
                    None => Curve2::BSpline(cut.on_face.clone()),
                };
                inserts.push((cut.face.face, (e, pcurve, stop.vertices)));
            }
        }
        across_edges.push(edges);
        stops.push([a, b]);
    }
    for (m, vs) in mitres.iter().zip(&mitre_vertices) {
        for &(edge, vertex, t, k) in &m.edges {
            trims.push((edge, vertex, vs[k], t));
        }
    }
    // The contact curves as edges, between the fillets' stops.
    let mut contact_edges: Vec<[EdgeId; 2]> = Vec::with_capacity(stripes.len());
    for (i, s) in stripes.iter().enumerate() {
        let mut edges = [EdgeId::INVALID; 2];
        for j in 0..2 {
            let (a, b) = (&stops[i][0], &stops[i][1]);
            let e = model.add_edge(a.vertices[j], b.vertices[j], tolerance)?;
            model.set_edge_curve(e.into(), Curve::BSpline(s.rolling.contacts[j].clone()), tolerance)?;
            model.set_edge_range(e.into(), a.v[j], b.v[j])?;
            edges[j] = e;
        }
        contact_edges.push(edges);
    }
    // The edges shortened at fillets' ends, each once.
    let mut shortened: Vec<(EdgeId, EdgeId, f64, f64)> = Vec::new();
    let mut done: Vec<EdgeId> = Vec::new();
    for &(old, ..) in &trims {
        if done.contains(&old) || removed.contains(&old) {
            continue;
        }
        done.push(old);
        let (a, b) = model.edge_vertices(Shape::from(old));
        let (curve, first, last) = model.edge_curve(Shape::from(old)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        let (a_id, b_id) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
        let (mut from, mut to, mut lo, mut hi) = (a_id, b_id, first, last);
        for &(_, at, new_vertex, t) in trims.iter().filter(|c| c.0 == old) {
            if at == a_id {
                (from, lo) = (new_vertex, t);
            } else if at == b_id {
                (to, hi) = (new_vertex, t);
            }
            reach_curve(model, new_vertex, &curve, t)?;
        }
        // (No closer than the edge it replaces holds to its faces.)
        let kept = tolerance.max(model.tolerance(Shape::from(old)));
        let e = model.add_edge(from, to, kept)?;
        model.set_edge_curve(e.into(), curve, kept)?;
        model.set_edge_range(e.into(), lo, hi)?;
        crate::boolean::record(model.uid(), old.into(), e.into());
        shortened.push((old, e, lo, hi));
    }
    // The faces: each touched once, with every change to it.
    let mut touched: Vec<Shape> = Vec::new();
    for s in stripes {
        for side in &s.rolling.sides {
            if !touched.iter().any(|f| f.node == side.face.node) {
                touched.push(side.face);
            }
        }
    }
    for (face, _) in &inserts {
        if !touched.iter().any(|f| f.node == face.node) {
            touched.push(*face);
        }
    }
    let mut rebuilt: Vec<(Shape, Shape)> = Vec::new();
    for &face in &touched {
        let mut replaces: Vec<(EdgeId, Vec<(EdgeId, f64, f64)>, Curve2)> = Vec::new();
        // (A side face whose contact leaves its patch takes its surface carried on.)
        let mut side_carried: Option<Surface> = None;
        for (i, s) in stripes.iter().enumerate() {
            for (j, side) in s.rolling.sides.iter().enumerate() {
                if side.face.node != face.node {
                    continue;
                }
                let pcurve = match project::analytic_pcurve(&Curve::BSpline(s.rolling.contacts[j].clone()), &side.surface) {
                    Some(p) => p,
                    None => Curve2::BSpline(fit::interpolate_with_parameters(
                        &unwrapped(s.rolling.sections.iter().map(|x| x.contact_uv[j]).collect(), &side.surface),
                        &s.ts,
                    )?),
                };
                cover_pcurve(model, contact_edges[i][j], &Curve::BSpline(s.rolling.contacts[j].clone()), &pcurve, &side.surface, stops[i][0].v[j], stops[i][1].v[j])?;
                let (a, b) = (stops[i][0].v[j], stops[i][1].v[j]);
                if off_patch(model, side, &pcurve, a.min(b), a.max(b)) {
                    side_carried = Some(side.surface.clone());
                }
                replaces.push((s.rolling.id, vec![(contact_edges[i][j], stops[i][0].v[j], stops[i][1].v[j])], pcurve));
            }
        }
        let replace_refs: Vec<(EdgeId, &[(EdgeId, f64, f64)], Curve2)> = replaces.iter().map(|(old, pieces, p)| (*old, pieces.as_slice(), p.clone())).collect();
        let face_inserts: Vec<(EdgeId, Curve2, [VertexId; 2])> = inserts.iter().filter(|(f, _)| f.node == face.node).map(|(_, x)| x.clone()).collect();
        // (An end face whose surface a cut carried on past its patch takes it.)
        let carried = ends.iter().flatten().find_map(|end| match end {
            End::Face(cut) if cut.carried && cut.face.face.node == face.node => Some(&cut.face.surface),
            _ => None,
        }).or(side_carried.as_ref());
        let new = rebuild_face(model, face, &replace_refs, &shortened, &removed, &face_inserts, carried)?;
        crate::boolean::record(model.uid(), face.node, new.node);
        rebuilt.push((face, new));
    }
    // The fillet faces: across at the start (u from 0 to 1), along the second contact, back
    // across at the end, back along the first: counter-clockwise in (u, v).
    let line = |o: (f64, f64), d: (f64, f64)| Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1)));
    let mut added: Vec<Shape> = Vec::new();
    for (i, s) in stripes.iter().enumerate() {
        let [a, b] = &stops[i];
        let mut uses: Vec<Use> = Vec::new();
        match &a.mitre {
            Some(pieces) => uses.extend(pieces.iter().map(|(e, o, p)| (*e, *o, p.clone(), 0.0, 1.0))),
            None => uses.push((across_edges[i][0], Orientation::Forward, a.on_fillet.clone(), 0.0, 1.0)),
        }
        uses.push((contact_edges[i][1], Orientation::Forward, line((1.0, 0.0), (0.0, 1.0)), a.v[1], b.v[1]));
        match &b.mitre {
            Some(pieces) => uses.extend(pieces.iter().rev().map(|(e, o, p)| (*e, o.reversed(), p.clone(), 0.0, 1.0))),
            None => uses.push((across_edges[i][1], Orientation::Reversed, b.on_fillet.clone(), 0.0, 1.0)),
        }
        uses.push((contact_edges[i][0], Orientation::Reversed, line((0.0, 0.0), (0.0, 1.0)), a.v[0], b.v[0]));
        let fillet = new_face(model, s.surface.clone(), vec![(uses, Orientation::Forward)], tolerance)?;
        crate::boolean::record(model.uid(), s.rolling.id.into(), fillet.node);
        // Outward: away from the ball's centre for a convex edge (the ball is in the solid).
        // (At the section nearest the middle of the edge: not one past its ends, where the
        // ball may lie the other way.)
        let middle = 0.5 * (s.rolling.t0 + s.rolling.t1);
        let mid = s.rolling.sections.iter().min_by(|a, b| (a.t - middle).abs().total_cmp(&(b.t - middle).abs())).expect("sections");
        let (p, du, dv) = s.surface.d1(0.5, mid.t);
        let away = (p - mid.center).dot(du.cross(dv)) > 0.0;
        added.push(if away == s.rolling.convex { fillet } else { fillet.reversed() });
    }
    // The corner patches, each using the fillets' stop curves the other way from them.
    let fillet_faces = added.clone();
    for (c, corner) in corners.iter().enumerate() {
        let face = patch(model, corner, &corner_vertices[c], &stops.iter().map(|s| [(s[0].vertices, s[0].curve.clone()), (s[1].vertices, s[1].curve.clone())]).collect::<Vec<_>>(), &across_edges, tolerance)?;
        let (i, k) = corner.stripes[0];
        let edge = across_edges[i][k];
        let used = |model: &Model, face: Shape| -> Option<Orientation> {
            model.explore(face, NodeKind::CoEdge).find(|c| model[crate::topology::CoEdgeId::try_from(c.node).expect("a coedge")].edge == edge).map(|c| {
                let id = crate::topology::CoEdgeId::try_from(c.node).expect("a coedge");
                model[id].orientation.compose(c.orientation)
            })
        };
        let (a, b) = (used(model, fillet_faces[i]), used(model, face));
        added.push(if a.is_some() && a == b { face.reversed() } else { face });
    }
    // (A contact that stays put, as where the ball stands still on the axis of a corner's
    // cylinder as wide as it: its edge of no length collapsed to a point, the side face's
    // boundary passing through it and the fillet face closing there on a degenerate edge,
    // as at a pole.)
    for (i, edges) in contact_edges.iter().enumerate() {
        for &e in edges {
            let (a, b) = model.edge_vertices(Shape::from(e));
            // (Of no length at all: the network's tolerance grows with the contacts' deviation,
            // far past what a contact that moves can be short of.)
            if a.node == b.node || !(crate::properties::length(model, Shape::from(e)).mass <= crate::tolerance::LINEAR) {
                continue;
            }
            let (a, b) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
            let gap = model.point(Shape::from(a)).distance(model.point(Shape::from(b)));
            if gap > crate::tolerance::LINEAR {
                continue;
            }
            model.replace_node(b, a)?;
            for parent in model.parents(e.into()) {
                let Ok(coedge) = CoEdgeId::try_from(parent) else { continue };
                let wire = model[coedge].wire;
                if model.face_of_wire(wire).is_some_and(|f| f != FaceId::try_from(added[i].node).expect("a face")) {
                    model.remove_coedge(wire, coedge)?;
                }
            }
            model.clear_edge_curve(Shape::from(e))?;
        }
    }
    // The solid with the faces replaced and the new ones added (or, for a shape of shells
    // only, the shell).
    let mut solids = Vec::new();
    let mut holders: Vec<Shape> = model.explore(shape, NodeKind::Solid).collect();
    let shells_only = holders.is_empty();
    if shells_only {
        holders = model.explore(shape, NodeKind::Shell).collect();
    }
    for solid in holders {
        let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
        if shells_only && !faces.iter().any(|f| rebuilt.iter().any(|(old, _)| old.node == f.node)) {
            solids.push(solid);
            continue;
        }
        let shell = model.new_empty(NodeKind::Shell);
        for f in faces {
            match rebuilt.iter().find(|(old, _)| old.node == f.node) {
                Some(&(_, new)) => model.add(shell, Shape { orientation: f.orientation, ..new })?,
                None => model.add(shell, f)?,
            }
        }
        for &f in &added {
            model.add(shell, f)?;
        }
        if shells_only {
            solids.push(Shape { orientation: solid.orientation, ..shell });
            continue;
        }
        let new_solid = model.new_empty(NodeKind::Solid);
        model.add(new_solid, shell)?;
        solids.push(new_solid);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for s in solids {
        model.add(compound, s)?;
    }
    let snapped: Vec<VertexId> = ends.iter().flatten().filter_map(|end| if let End::Face(cut) = end { Some(cut.onto) } else { None }).flat_map(|o| o.into_iter().flatten()).collect();
    reach_snapped(model, shape, compound, &snapped)?;
    let extended = ends.iter().flatten().any(|end| matches!(end, End::Face(cut) if cut.extended.iter().any(|&x| x)));
    joined_extensions(model, compound, extended)
}

/// The patch closing a corner: a Coons surface through the three fillets' stop curves
/// (the third side collapsed to a point), outward like the fillets.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn patch(
    model: &mut Model,
    corner: &Corner,
    vertices: &[VertexId],
    stops: &[[([VertexId; 2], BSplineCurve); 2]],
    across: &[[EdgeId; 2]],
    tolerance: f64,
) -> Result<Shape> {
    // The three stop curves, each joining two of the corner's points.
    let sides: Vec<(EdgeId, [usize; 2], BSplineCurve)> = corner
        .stripes
        .iter()
        .map(|&(i, k)| {
            let (vs, curve) = &stops[i][k];
            let index = |v: VertexId| vertices.iter().position(|&w| w == v).expect("a corner vertex");
            (across[i][k], [index(vs[0]), index(vs[1])], curve.clone())
        })
        .collect();
    // Corner points 0, 1, 2: the bottom from 0 to 1, the right from 1 to 2, the left from
    // 0 to 2 (each curve run from one point to the other, reversed where needed), the top
    // collapsed at 2.
    let find = |a: usize, b: usize| -> (usize, bool) {
        let k = sides.iter().position(|s| (s.1 == [a, b]) || (s.1 == [b, a])).expect("a side");
        (k, sides[k].1 == [a, b])
    };
    let side_point = |(k, forward): (usize, bool), w: f64| sides[k].2.point(if forward { w } else { 1.0 - w });
    let (bottom, right, left) = (find(0, 1), find(1, 2), find(0, 2));
    let p = [0, 1, 2].map(|n| model[vertices[n]].point);
    let coons = |u: f64, v: f64| -> Point {
        let b = side_point(bottom, u);
        let l = side_point(left, v);
        let r = side_point(right, v);
        let sum = (b - Point::ORIGIN) * (1.0 - v) + (p[2] - Point::ORIGIN) * v + (l - Point::ORIGIN) * (1.0 - u) + (r - Point::ORIGIN) * u
            - ((p[0] - Point::ORIGIN) * ((1.0 - u) * (1.0 - v)) + (p[1] - Point::ORIGIN) * (u * (1.0 - v)) + (p[2] - Point::ORIGIN) * v);
        Point::ORIGIN + sum
    };
    // Interpolated on a grid: rows along u, then their poles along v.
    let params: Vec<f64> = (0..=PATCH).map(|k| k as f64 / PATCH as f64).collect();
    let rows: Vec<BSplineCurve> = params
        .iter()
        .map(|&v| fit::interpolate_with_parameters(&params.iter().map(|&u| coons(u, v)).collect::<Vec<_>>(), &params))
        .collect::<Result<_>>()?;
    let count = rows[0].poles().len();
    let columns: Vec<BSplineCurve> = (0..count)
        .map(|n| fit::interpolate_with_parameters(&rows.iter().map(|r| r.poles()[n]).collect::<Vec<_>>(), &params))
        .collect::<Result<_>>()?;
    let poles: Vec<Vec<Point>> = columns.iter().map(|c| c.poles().to_vec()).collect();
    let u = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
    let v = KnotVector::new(columns[0].knots().to_vec(), columns[0].multiplicities().to_vec(), columns[0].degree());
    let surface = Surface::BSpline(BSplineSurface::try_new(poles, u, v)?);
    // The boundary: counter-clockwise in (u, v).
    let segment = |a: (f64, f64), b: (f64, f64)| Curve2::BSpline(BSplineCurve2::new(vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)], vec![0.0, 1.0], vec![2, 2], 1));
    let along = |(k, forward): (usize, bool), from: (f64, f64), to: (f64, f64)| -> Use {
        // The edge runs from its first vertex to its second over [0, 1]: its pcurve runs
        // the same way, the coedge forward where the boundary does.
        if forward {
            (sides[k].0, Orientation::Forward, segment(from, to), 0.0, 1.0)
        } else {
            (sides[k].0, Orientation::Reversed, segment(to, from), 0.0, 1.0)
        }
    };
    let apex = model.add_edge(vertices[2], vertices[2], tolerance)?;
    let uses: Vec<Use> = vec![
        along(bottom, (0.0, 0.0), (1.0, 0.0)),
        along(right, (1.0, 0.0), (1.0, 1.0)),
        (apex, Orientation::Forward, segment((1.0, 1.0), (0.0, 1.0)), 0.0, 1.0),
        {
            let (k, forward) = left;
            // The left side runs from point 0 up to 2; the boundary comes down it.
            if forward {
                (sides[k].0, Orientation::Reversed, segment((0.0, 0.0), (0.0, 1.0)), 0.0, 1.0)
            } else {
                (sides[k].0, Orientation::Forward, segment((0.0, 1.0), (0.0, 0.0)), 0.0, 1.0)
            }
        },
    ];
    let tolerance = tolerance.max(deviation_from_sides(&surface, &sides, bottom, right, left));
    new_face(model, surface, vec![(uses, Orientation::Forward)], tolerance)
}

/// How far the patch's boundary strays from the stop curves it is built through.
fn deviation_from_sides(surface: &Surface, sides: &[(EdgeId, [usize; 2], BSplineCurve)], bottom: (usize, bool), right: (usize, bool), left: (usize, bool)) -> f64 {
    let at = |(k, forward): (usize, bool), w: f64| sides[k].2.point(if forward { w } else { 1.0 - w });
    let mut worst: f64 = 0.0;
    for n in 0..=64 {
        let w = n as f64 / 64.0;
        worst = worst
            .max(surface.point(w, 0.0).distance(at(bottom, w)))
            .max(surface.point(1.0, w).distance(at(right, w)))
            .max(surface.point(0.0, w).distance(at(left, w)));
    }
    2.0 * worst
}
