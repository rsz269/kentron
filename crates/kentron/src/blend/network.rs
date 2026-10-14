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
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{cover_pcurve, end_cut, extension, fillet_surface, new_face, rebuild_face, rolling_past, unwrapped, EndCut, Rolling, Use, NOT_REACHED};

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
    /// The edge between the far faces: its vertex at the corner, and the parameter of the
    /// last point.
    edge: (EdgeId, VertexId, f64),
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
    // (Further past the ends where the end faces are not reached.)
    match fillet_network_past(model, shape, edges, 1) {
        Err(Error::InvalidInput(m)) if m == NOT_REACHED => fillet_network_past(model, shape, edges, 3),
        other => other,
    }
}

/// [`fillet_network`] with the sections carried `far` times as far past the edges' ends.
fn fillet_network_past(model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)], far: usize) -> Result<Shape> {
    let mut stripes = Vec::with_capacity(edges.len());
    for &(edge, radius) in edges {
        let rolling = rolling_past(model, shape, edge, radius, far)?;
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
                if let Some(end) = free_end(model, shape, s, k)? {
                    ends[i][k] = Some(end);
                    continue;
                }
                let cut = end_cut(model, shape, s.rolling.id, if k == 0 { s.rolling.va } else { s.rolling.vb }, &s.rolling.sides, &s.surface, &s.ts, k == 1, s.rolling.target)?;
                ends[i][k] = Some(End::Face(cut));
            }
            [a, b] => {
                let (m, ea, eb) = mitre(model, shape, vertex, &stripes, a, b)?;
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
    let (va1, vb1, p1) = crossing(&sa.rolling.contacts[ia], &sa.ts, a.1, &sb.rolling.contacts[ib], &sb.ts, b.1)?;
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
    let m = Mitre { points: mitre_points, curve, cut, edge: (third, vertex, on_edge.parameter), deviation };
    Ok(if s_is_a { (m, ps, pl) } else { (m, pl, ps) })
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
        if on.distance > (1e3 * s.rolling.target).max(1e-6 * (1.0 + (point - Point::ORIGIN).length())) {
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
    let _ = (shape, vertex);
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
        let (t1, t2, p) = crossing(&a.rolling.contacts[s1], &a.ts, meeting[n1].1, &b.rolling.contacts[s2], &b.ts, meeting[n2].1)?;
        let _ = model;
        points.push((face, p));
        (v[n1][s1], at[n1][s1]) = (t1, f);
        (v[n2][s2], at[n2][s2]) = (t2, f);
    }
    let pending = (0..meeting.len()).map(|n| PendingEnd { at: at[n], v: v[n] }).collect();
    Ok((Corner { points, stripes: meeting.to_vec() }, pending))
}

/// Where two contact curves on one face cross near the ends `end_a`, `end_b` of their
/// fillets (0 the start of the edge, 1 its end): their parameters and the point.
fn crossing(a: &BSplineCurve, ts_a: &[f64], end_a: usize, b: &BSplineCurve, ts_b: &[f64], end_b: usize) -> Result<(f64, f64, Point)> {
    // Over the half of each curve towards the corner (its ends past the edge included).
    let half = |ts: &[f64], end: usize| {
        let (lo, hi) = (ts[0], ts[ts.len() - 1]);
        let mid = 0.5 * (lo + hi);
        if end == 0 { (lo, mid) } else { (mid, hi) }
    };
    let ((a0, a1), (b0, b1)) = (half(ts_a, end_a), half(ts_b, end_b));
    let n = 64;
    let mut best = (f64::INFINITY, a0, b0);
    for i in 0..=n {
        let ta = a0 + (a1 - a0) * i as f64 / n as f64;
        let pa = a.point(ta);
        for j in 0..=n {
            let tb = b0 + (b1 - b0) * j as f64 / n as f64;
            let d = pa.distance(b.point(tb));
            if d < best.0 {
                best = (d, ta, tb);
            }
        }
    }
    // Newton on the gap's square.
    let (mut ta, mut tb) = (best.1, best.2);
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
    let (pa, pb) = (a.point(ta), b.point(tb));
    let size = 1.0 + (pa - Point::ORIGIN).length();
    if pa.distance(pb) > 1e-6 * size {
        return Err(Error::InvalidInput("fillets at a corner whose contact curves do not cross".into()));
    }
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
    // Within twice the contact curves' distance from their faces.
    let mut deviation: f64 = 0.0;
    let mut size: f64 = 1.0;
    for s in stripes {
        let (a, b) = (s.ts[0], s.ts[s.ts.len() - 1]);
        for i in 0..2 {
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
    let mut across_edges: Vec<[EdgeId; 2]> = Vec::with_capacity(stripes.len());
    for (i, s) in stripes.iter().enumerate() {
        let mut pair = Vec::with_capacity(2);
        for end in &ends[i] {
            pair.push(match end {
                End::Face(cut) => {
                    let vertices = [0, 1].map(|j| model.add_vertex(s.rolling.contacts[j].point(cut.v[j]), tolerance));
                    for (j, &(old, t)) in cut.edges.iter().enumerate() {
                        if cut.extended[j] {
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
            model.set_edge_curve(e.into(), Curve::BSpline(stop.curve.clone()), tolerance)?;
            model.set_edge_range(e.into(), 0.0, 1.0)?;
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
        trims.push((m.edge.0, m.edge.1, *vs.last().expect("mitre points"), m.edge.2));
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
        if done.contains(&old) {
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
        }
        let e = model.add_edge(from, to, tolerance.max(model.tolerance(Shape::from(old))))?;
        model.set_edge_curve(e.into(), curve, tolerance)?;
        model.set_edge_range(e.into(), lo, hi)?;
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
                replaces.push((s.rolling.id, vec![(contact_edges[i][j], stops[i][0].v[j], stops[i][1].v[j])], pcurve));
            }
        }
        let replace_refs: Vec<(EdgeId, &[(EdgeId, f64, f64)], Curve2)> = replaces.iter().map(|(old, pieces, p)| (*old, pieces.as_slice(), p.clone())).collect();
        let face_inserts: Vec<(EdgeId, Curve2, [VertexId; 2])> = inserts.iter().filter(|(f, _)| f.node == face.node).map(|(_, x)| x.clone()).collect();
        let new = rebuild_face(model, face, &replace_refs, &shortened, &face_inserts)?;
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
        // Outward: away from the ball's centre for a convex edge (the ball is in the solid).
        let mid = &s.rolling.sections[s.rolling.sections.len() / 2];
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
    Ok(compound)
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
