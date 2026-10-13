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

use super::{end_cut, fillet_surface, new_face, rebuild_face, rolling, unwrapped, EndCut, Rolling, Use};

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

/// `shape` with `edges` (meeting each other only three at a time, at corners where three
/// faces meet) rounded by rolling-ball fillets of the given radii (as laws over the fraction
/// of each edge's length), built together.
pub(in crate::blend) fn fillet_network(model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)]) -> Result<Shape> {
    let mut stripes = Vec::with_capacity(edges.len());
    for &(edge, radius) in edges {
        let rolling = rolling(model, shape, edge, radius)?;
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
    let mut vertices: Vec<(VertexId, Vec<(usize, usize)>)> = at_vertex.into_iter().collect();
    vertices.sort_by_key(|(v, _)| v.index());
    for (vertex, meeting) in vertices {
        match meeting[..] {
            [(i, k)] => {
                let s = &stripes[i];
                let cut = end_cut(model, shape, s.rolling.id, if k == 0 { s.rolling.va } else { s.rolling.vb }, &s.rolling.sides, &s.surface, &s.ts, k == 1, s.rolling.target)?;
                ends[i][k] = Some(End::Face(cut));
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
    build(model, shape, &stripes, &ends, &corners)
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
fn build(model: &mut Model, shape: Shape, stripes: &[Stripe], ends: &[[End; 2]], corners: &[Corner]) -> Result<Shape> {
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
    let tolerance = (deviation * 2.0).max(crate::tolerance::LINEAR).max(1e-8 * size);
    // The corners' points as vertices.
    let corner_vertices: Vec<Vec<VertexId>> = corners.iter().map(|c| c.points.iter().map(|&(_, p)| model.add_vertex(p, tolerance)).collect()).collect();
    // Per fillet and end: the contacts' parameters there and their vertices; the curve
    // across (3D and on the fillet) and, for an end on a face, its pcurve there.
    struct Stop {
        v: [f64; 2],
        vertices: [VertexId; 2],
        curve: BSplineCurve,
        on_fillet: Curve2,
    }
    let mut stops: Vec<[Stop; 2]> = Vec::with_capacity(stripes.len());
    // (edge, its vertex at a fillet's end, the vertex replacing it, the parameter there)
    let mut trims: Vec<(EdgeId, VertexId, VertexId, f64)> = Vec::new();
    // (face, edge, pcurve, the edge's vertices) inserted into end faces
    let mut inserts: Vec<(Shape, (EdgeId, Curve2, [VertexId; 2]))> = Vec::new();
    let mut across_edges: Vec<[EdgeId; 2]> = Vec::with_capacity(stripes.len());
    for (i, s) in stripes.iter().enumerate() {
        let mut pair = Vec::with_capacity(2);
        for end in &ends[i] {
            pair.push(match end {
                End::Face(cut) => {
                    let vertices = [0, 1].map(|j| model.add_vertex(s.rolling.contacts[j].point(cut.v[j]), tolerance));
                    for (j, &(old, t)) in cut.edges.iter().enumerate() {
                        trims.push((old, cut.vertex, vertices[j], t));
                    }
                    Stop { v: cut.v, vertices, curve: cut.curve.clone(), on_fillet: Curve2::BSpline(cut.on_fillet.clone()) }
                }
                End::Corner { corner, at, v } => {
                    let (curve, on_fillet) = stop_curve(s, v[0], v[1])?;
                    Stop { v: *v, vertices: [corner_vertices[*corner][at[0]], corner_vertices[*corner][at[1]]], curve, on_fillet: Curve2::BSpline(on_fillet) }
                }
            });
        }
        let [a, b]: [Stop; 2] = pair.try_into().map_err(|_| Error::Failed("a fillet's ends".into()))?;
        // The edges across at both ends.
        let mut edges = [EdgeId::INVALID; 2];
        for (k, stop) in [&a, &b].into_iter().enumerate() {
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
        let uses: Vec<Use> = vec![
            (across_edges[i][0], Orientation::Forward, a.on_fillet.clone(), 0.0, 1.0),
            (contact_edges[i][1], Orientation::Forward, line((1.0, 0.0), (0.0, 1.0)), a.v[1], b.v[1]),
            (across_edges[i][1], Orientation::Reversed, b.on_fillet.clone(), 0.0, 1.0),
            (contact_edges[i][0], Orientation::Reversed, line((0.0, 0.0), (0.0, 1.0)), a.v[0], b.v[0]),
        ];
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
    // The solid with the faces replaced and the new ones added.
    let mut solids = Vec::new();
    for solid in model.explore(shape, NodeKind::Solid).collect::<Vec<_>>() {
        let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
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
