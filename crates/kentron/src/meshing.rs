//! Triangulating shapes.
//!
//! 1. Each edge is cut into a polyline, adaptively, until each chord lies within the linear
//!    deflection of the curve and turns by at most the angular deflection; faces sharing an
//!    edge share its nodes.
//! 2. Each face's boundary is mapped into its surface's parameters through the pcurves; on
//!    curved surfaces a grid of interior points is added, spaced by the surface's curvature.
//! 3. The points are triangulated (a Delaunay triangulation in parameter space scaled to
//!    the surface's metric, with the boundary segments then forced in by edge flips), and
//!    the triangles outside the face dropped.

use crate::hashing::{FastMap as HashMap, FastSet};

use crate::error::{Error, Result};
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Transform, Transformable};
use crate::mesh::{MeshPolyline, Polyline, TriangleMesh};
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeKind, Orientation, Shape, WireId};

/// How finely [`mesh_with`] triangulates.
#[doc(alias = "IMeshTools_Parameters")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshOptions {
    /// The largest distance between the mesh and the geometry.
    #[doc(alias = "Deflection")]
    pub deflection: f64,
    /// The largest turn, in radians, between neighbouring segments along a curve.
    #[doc(alias = "Angle")]
    pub angle: f64,
}

impl Default for MeshOptions {
    fn default() -> Self {
        Self { deflection: 1e-3, angle: 0.5 }
    }
}

/// Stores a triangle mesh on every face of `shape`, and on every coedge the chain of that
/// mesh's nodes along its edge (so adjacent faces' meshes match along shared edges; a seam
/// gets one chain per coedge). The meshes deviate from the exact geometry by at most
/// `deflection`. Existing meshes that are fine enough are kept.
#[doc(alias = "BRepMesh_IncrementalMesh")]
pub fn mesh(model: &mut Model, shape: Shape, deflection: f64) -> Result<()> {
    mesh_with(model, shape, MeshOptions { deflection, ..MeshOptions::default() })
}

/// [`mesh`] with an angular deflection too.
pub fn mesh_with(model: &mut Model, shape: Shape, options: MeshOptions) -> Result<()> {
    if !(options.deflection > 0.0) || !(options.angle > 0.0) {
        return Err(Error::InvalidInput("mesh deflections must be positive".into()));
    }
    let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..shape };
    let unplaced = |s: Shape| Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..s };
    // The edges' polylines, then the faces' meshes (on them), each found for all at once and
    // then stored in order.
    let edges: Vec<Shape> = model
        .unique_subshapes(local, NodeKind::Edge)
        .into_iter()
        .map(unplaced)
        .filter(|&edge| !model.edge_polyline(edge).is_some_and(|p| p.deflection() > 0.0 && p.deflection() <= options.deflection))
        .collect();
    let edges = once_each(edges);
    let model_now: &Model = model;
    let polylines = crate::parallel::map(&edges, |&edge| edge_polyline(model_now, edge, options));
    for (&edge, polyline) in edges.iter().zip(polylines) {
        if let Some(polyline) = polyline {
            model.set_edge_polyline(edge, Some(polyline))?;
        }
    }
    // As OCCT does, a mesh recording no deviation (of a plane, say) is exact.
    let faces: Vec<Shape> = model
        .unique_subshapes(local, NodeKind::Face)
        .into_iter()
        .map(unplaced)
        .filter(|&face| !model.face_mesh(face).is_some_and(|m| m.deflection() <= options.deflection))
        .collect();
    let faces = once_each(faces);
    let model_now: &Model = model;
    let meshes = crate::parallel::map(&faces, |&face| face_mesh(model_now, face, options));
    for (&face, found) in faces.iter().zip(meshes) {
        if let Some((mesh, chains)) = found? {
            model.set_face_mesh(face, Some(mesh))?;
            for (coedge, nodes) in chains {
                model.set_coedge_mesh_polyline(coedge, Some(MeshPolyline::new(nodes).with_deflection(options.deflection)))?;
            }
        }
    }
    Ok(())
}

/// The shapes in order, each node once (unplaced, a shape placed several times is one).
fn once_each(shapes: Vec<Shape>) -> Vec<Shape> {
    let mut seen = crate::hashing::FastSet::default();
    shapes.into_iter().filter(|s| seen.insert(s.node)).collect()
}

/// The triangles of `shape`'s faces, meshed as [`mesh_with`] does (and stored there),
/// placed and turned to face outwards, with the nodes shared between faces merged.
#[doc(alias = "tessellate")]
pub fn triangulate(model: &mut Model, shape: Shape, options: MeshOptions) -> Result<(Vec<Point>, Vec<[u32; 3]>)> {
    mesh_with(model, shape, options)?;
    let mut nodes: Vec<Point> = Vec::new();
    let mut index: HashMap<[u64; 3], u32> = HashMap::default();
    let mut triangles = Vec::new();
    for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
        let Some(mesh) = model.face_mesh(Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face }) else { continue };
        let ids: Vec<u32> = mesh
            .nodes()
            .iter()
            .map(|p| {
                let q = p.transformed(&face.placement);
                *index.entry([q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]).or_insert_with(|| {
                    nodes.push(q);
                    (nodes.len() - 1) as u32
                })
            })
            .collect();
        for t in mesh.triangles() {
            let [a, b, c] = t.map(|i| ids[i as usize]);
            if a == b || b == c || a == c {
                continue;
            }
            triangles.push(if face.orientation == Orientation::Reversed { [a, c, b] } else { [a, b, c] });
        }
    }
    Ok((nodes, triangles))
}

/// An edge's polyline: parameters and points within the deflections.
fn edge_polyline(model: &Model, edge: Shape, options: MeshOptions) -> Option<Polyline> {
    let (curve, first, last) = model.edge_curve(edge)?;
    let params = curve_parameters(&curve, first, last, options);
    let mut points: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
    // The ends at the vertices exactly (shared with the other edges there).
    let (a, b) = model.edge_vertices(edge);
    points[0] = model.point(a);
    let n = points.len();
    points[n - 1] = model.point(b);
    Some(Polyline::new(points).with_parameters(params).with_deflection(options.deflection))
}

/// Parameters along a curve: subdivided until each chord is within the deflection and
/// turns by at most the angle.
fn curve_parameters(curve: &crate::geometry::curve::Curve, first: f64, last: f64, options: MeshOptions) -> Vec<f64> {
    use crate::geometry::curve::Curve;
    if let Curve::Line(_) = curve {
        return vec![first, last];
    }
    let start = 4;
    let mut params = vec![first];
    let mut stack: Vec<(f64, f64, u32)> = (0..start)
        .rev()
        .map(|i| {
            let (a, b) = (first + (last - first) * i as f64 / start as f64, first + (last - first) * (i + 1) as f64 / start as f64);
            (a, b, 0)
        })
        .collect();
    while let Some((a, b, depth)) = stack.pop() {
        let (pa, ta) = curve.d1(a);
        let (pb, tb) = curve.d1(b);
        let m = 0.5 * (a + b);
        let pm = curve.point(m);
        let chord = pa + (pb - pa) * 0.5;
        let sag = pm.distance(chord);
        let (la, lb) = (ta.length(), tb.length());
        let turn = if la > 0.0 && lb > 0.0 { (ta.dot(tb) / (la * lb)).clamp(-1.0, 1.0).acos() } else { 0.0 };
        if depth < 20 && (sag > options.deflection || turn > options.angle) {
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
        } else {
            params.push(b);
        }
    }
    params
}

/// A boundary node of a face: its (u, v) and the key of its 3D node.
#[derive(Clone, Copy)]
struct BoundaryNode {
    uv: Point2,
    key: NodeKey,
}

/// What a mesh node stands for: a node of an edge's polyline (or a vertex), or a point
/// inside the face.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum NodeKey {
    Vertex(crate::topology::VertexId),
    OnEdge(EdgeId, usize),
    Inside(usize),
}

/// A face's mesh and, per coedge, its chain of mesh nodes; `None` for a face without a
/// surface or pcurves.
#[allow(clippy::type_complexity)]
fn face_mesh(model: &Model, face: Shape, options: MeshOptions) -> Result<Option<(TriangleMesh, Vec<(CoEdgeId, Vec<u32>)>)>> {
    let Some(surface) = model.surface(face).cloned() else { return Ok(None) };
    let id = FaceId::try_from(face.node).expect("a face");
    // The boundary loops in (u, v), with their 3D nodes.
    let mut loops: Vec<Vec<BoundaryNode>> = Vec::new();
    let mut chains: Vec<(CoEdgeId, Vec<NodeKey>)> = Vec::new();
    let mut points: HashMap<NodeKey, Point> = HashMap::default();
    for wire in model.children_local(face).collect::<Vec<_>>() {
        let wid = WireId::try_from(wire.node).expect("a wire");
        let mut coedges: Vec<(CoEdgeId, Orientation)> =
            model.get(wid).expect("a wire").coedges.iter().map(|&c| (c, model[c].orientation)).collect();
        if wire.orientation == Orientation::Reversed {
            coedges.reverse();
            for c in &mut coedges {
                c.1 = c.1.reversed();
            }
        }
        let mut nodes: Vec<BoundaryNode> = Vec::new();
        for (c, orientation) in coedges {
            let edge = model[c].edge;
            let Some((pcurve, pf, pl)) = model.effective_pcurve(c) else { return Ok(None) };
            let (a, b) = model.edge_vertices(Shape::from(edge));
            let (va, vb) = (crate::topology::VertexId::try_from(a.node).expect("a vertex"), crate::topology::VertexId::try_from(b.node).expect("a vertex"));
            // (parameter on the pcurve, key, point) along the edge's direction.
            let mut along: Vec<(f64, NodeKey, Point)> = match model.edge_polyline(Shape::from(edge)) {
                Some(polyline) if model.edge_curve(Shape::from(edge)).is_some() => {
                    let (ef, el) = model.edge_range(Shape::from(edge));
                    let params = polyline.parameters().expect("meshed edges have parameters");
                    let n = params.len();
                    params
                        .iter()
                        .zip(polyline.nodes())
                        .enumerate()
                        .map(|(i, (&t, &p))| {
                            // The edge's parameter mapped onto the pcurve's range.
                            let s = if el != ef { (t - ef) / (el - ef) } else { 0.0 };
                            let key = if i == 0 { NodeKey::Vertex(va) } else if i == n - 1 { NodeKey::Vertex(vb) } else { NodeKey::OnEdge(edge, i) };
                            (pf + (pl - pf) * s, key, p)
                        })
                        .collect()
                }
                _ => {
                    // A degenerate edge: its vertex all along a stretch of (u, v).
                    let p = model.point(a);
                    let k = 8;
                    (0..=k).map(|i| (pf + (pl - pf) * i as f64 / k as f64, NodeKey::Vertex(va), p)).collect()
                }
            };
            if orientation == Orientation::Reversed {
                along.reverse();
            }
            let mut chain = Vec::with_capacity(along.len());
            for (j, (t, key, p)) in along.iter().enumerate() {
                points.insert(*key, *p);
                chain.push(*key);
                // The last node starts the next coedge.
                if j + 1 < along.len() {
                    nodes.push(BoundaryNode { uv: uv_of(&surface, *p, pcurve.point(*t)), key: *key });
                }
            }
            chains.push((c, chain));
        }
        if nodes.len() >= 3 {
            loops.push(nodes);
        }
    }
    if loops.is_empty() {
        return Ok(None);
    }
    // Scale (u, v) to lengths on the surface.
    let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
    for n in loops.iter().flatten() {
        lo = Point2::new(lo.x.min(n.uv.x), lo.y.min(n.uv.y));
        hi = Point2::new(hi.x.max(n.uv.x), hi.y.max(n.uv.y));
    }
    let (mut scale, steps) = metric(&surface, lo, hi, options);
    // Triangulated in units of the grid's cells, so that the triangles keep within a cell or
    // so along the way the surface bends, however long the cells are the other way.
    if (steps.0 > 1 || steps.1 > 1) && hi.x > lo.x && hi.y > lo.y {
        scale = (steps.0 as f64 / (hi.x - lo.x), steps.1 as f64 / (hi.y - lo.y));
    }
    let scaled = |p: Point2| Point2::new(p.x * scale.0, p.y * scale.1);
    // Interior points on a grid, away from the boundary.
    let polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().map(|n| scaled(n.uv)).collect()).collect();
    let mut inside: Vec<Point2> = Vec::new();
    if steps.0 > 1 || steps.1 > 1 {
        let (su, sv) = (scale.0 * (hi.x - lo.x) / steps.0 as f64, scale.1 * (hi.y - lo.y) / steps.1 as f64);
        let clearance = 0.3 * su.min(sv).max(1e-12);
        for i in 0..=steps.0 {
            for j in 0..=steps.1 {
                let uv = Point2::new(lo.x + (hi.x - lo.x) * i as f64 / steps.0 as f64, lo.y + (hi.y - lo.y) * j as f64 / steps.1 as f64);
                let q = scaled(uv);
                if inside_polygons(&polygons, q) && distance_to_polygons(&polygons, q) > clearance {
                    inside.push(uv);
                }
            }
        }
    }
    // The triangulation, refined where the triangles stray from a curved surface, then the
    // 3D mesh.
    let boundary_points: Vec<Point2> = polygons.iter().flatten().copied().collect();
    let boundary_3d: Vec<Point> = loops.iter().flatten().map(|n| points[&n.key]).collect();
    let mut segments = Vec::new();
    let mut offset = 0;
    for l in &loops {
        for i in 0..l.len() {
            segments.push((offset + i, offset + (i + 1) % l.len()));
        }
        offset += l.len();
    }
    let boundary_edges: FastSet<(usize, usize)> = segments.iter().flat_map(|&(a, b)| [(a, b), (b, a)]).collect();
    let curved = !matches!(surface, Surface::Plane(_));
    let unscale = |q: Point2| Point2::new(q.x / scale.0, q.y / scale.1);
    let mut triangles;
    let mut rounds = 0;
    loop {
        let mut all = boundary_points.clone();
        all.extend(inside.iter().map(|&p| scaled(p)));
        triangles = constrained_delaunay(&all, &segments);
        triangles.retain(|t| {
            let c = Point2::new((all[t[0]].x + all[t[1]].x + all[t[2]].x) / 3.0, (all[t[0]].y + all[t[1]].y + all[t[2]].y) / 3.0);
            inside_polygons(&polygons, c)
        });
        rounds += 1;
        if !curved || rounds > 8 || all.len() > 200_000 {
            break;
        }
        // Inner edges whose chords sag from the surface by more than the deflection get a
        // point at their middle.
        let point3 = |i: usize| if i < boundary_3d.len() { boundary_3d[i] } else { let uv = inside[i - boundary_3d.len()]; surface.point(uv.x, uv.y) };
        let mut added: FastSet<(usize, usize)> = FastSet::default();
        let mut new_points = Vec::new();
        for t in &triangles {
            for k in 0..3 {
                let (a, b) = (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]));
                if boundary_edges.contains(&(a, b)) || !added.insert((a, b)) {
                    continue;
                }
                let mid = unscale(Point2::new(0.5 * (all[a].x + all[b].x), 0.5 * (all[a].y + all[b].y)));
                let on_surface = surface.point(mid.x, mid.y);
                let (pa, pb) = (point3(a), point3(b));
                let chord = pa + (pb - pa) * 0.5;
                // Not on (or right by) the boundary, whose nodes are shared with the
                // neighbouring faces.
                let length = all[a].distance(all[b]);
                if on_surface.distance(chord) > options.deflection && distance_to_polygons(&polygons, scaled(mid)) > 0.01 * length {
                    new_points.push(mid);
                }
            }
        }
        if new_points.is_empty() {
            break;
        }
        inside.extend(new_points);
    }
    let deflection = if curved { achieved_deflection(&surface, &triangles, &boundary_3d, &inside, &polygons, unscale) } else { 0.0 };
    let mut keys: Vec<NodeKey> = loops.iter().flatten().map(|n| n.key).collect();
    for (k, &uv) in inside.iter().enumerate() {
        let key = NodeKey::Inside(k);
        points.insert(key, surface.point(uv.x, uv.y));
        keys.push(key);
    }
    let mut node_of: HashMap<NodeKey, u32> = HashMap::default();
    let mut nodes: Vec<Point> = Vec::new();
    let mut index = |key: NodeKey| *node_of.entry(key).or_insert_with(|| {
        nodes.push(points[&key]);
        (nodes.len() - 1) as u32
    });
    let mut kept = Vec::new();
    for [a, b, c] in triangles {
        let t = [index(keys[a]), index(keys[b]), index(keys[c])];
        if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
            kept.push(t);
        }
    }
    let chains = chains.into_iter().map(|(c, ks)| (c, ks.into_iter().map(&mut index).collect())).collect();
    let _ = id;
    Ok(Some((TriangleMesh::new(nodes, kept).with_deflection(deflection), chains)))
}

/// How far a face's triangles stray from its surface: the largest distance between the
/// middle of a triangle, or of one of its sides, and the surface point it stands for.
fn achieved_deflection(
    surface: &Surface,
    triangles: &[[usize; 3]],
    boundary_3d: &[Point],
    inside: &[Point2],
    polygons: &[Vec<Point2>],
    unscale: impl Fn(Point2) -> Point2,
) -> f64 {
    let boundary_uv: Vec<Point2> = polygons.iter().flatten().map(|&q| unscale(q)).collect();
    let uv = |i: usize| if i < boundary_uv.len() { boundary_uv[i] } else { inside[i - boundary_uv.len()] };
    let point3 = |i: usize| if i < boundary_3d.len() { boundary_3d[i] } else { let p = inside[i - boundary_3d.len()]; surface.point(p.x, p.y) };
    let gap = |ids: &[usize]| {
        let n = ids.len() as f64;
        let mid_uv = ids.iter().fold(Point2::new(0.0, 0.0), |m, &i| Point2::new(m.x + uv(i).x / n, m.y + uv(i).y / n));
        let mid = ids.iter().fold(Point::new(0.0, 0.0, 0.0), |m, &i| { let p = point3(i); Point::new(m.x + p.x / n, m.y + p.y / n, m.z + p.z / n) });
        surface.point(mid_uv.x, mid_uv.y).distance(mid)
    };
    triangles.iter().map(|t| gap(t).max(gap(&[t[0], t[1]])).max(gap(&[t[1], t[2]])).max(gap(&[t[2], t[0]]))).fold(0.0, f64::max)
}

/// The parameters of a boundary point: the surface's own (in the period of `estimate`, the
/// pcurve's), so that the boundary in (u, v) follows the edges exactly even where a pcurve
/// is parametrized otherwise than its edge; the estimate where the point has no parameters
/// of its own (a pole, an apex).
fn uv_of(surface: &Surface, p: Point, estimate: Point2) -> Point2 {
    if surface.point(estimate.x, estimate.y).distance(p) <= 1e-7 {
        return estimate;
    }
    let (u, v) = match surface.parameters(p) {
        Some(uv) => uv,
        None => {
            let (u, v, d) = crate::geometry::extrema::nearest_on_surface(p, surface, Some((estimate.x, estimate.y)));
            if d > 1e-3 {
                return estimate;
            }
            (u, v)
        }
    };
    let near = |x: f64, target: f64, period: Option<f64>| period.map_or(x, |q| x + ((target - x) / q).round() * q);
    let uv = Point2::new(near(u, estimate.x, surface.u_period()), near(v, estimate.y, surface.v_period()));
    // At a singular point u is free: the pcurve's.
    let (_, du, dv) = surface.d1(uv.x, uv.y);
    if du.cross(dv).length() <= 1e-12 * (du.length() + dv.length()).powi(2).max(1e-300) {
        return estimate;
    }
    uv
}

/// The scale of u and v to lengths on the surface (their mean speeds), and the numbers of
/// grid steps across `[lo, hi]` keeping within the deflections.
fn metric(surface: &Surface, lo: Point2, hi: Point2, options: MeshOptions) -> ((f64, f64), (usize, usize)) {
    if let Surface::Plane(_) = surface {
        return ((1.0, 1.0), (1, 1));
    }
    let n = 6;
    let (mut speed_u, mut speed_v, mut bend_u, mut bend_v, mut count) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0);
    let (mut turn_u, mut turn_v) = (0.0f64, 0.0f64);
    for i in 0..=n {
        for j in 0..=n {
            let (u, v) = (lo.x + (hi.x - lo.x) * i as f64 / n as f64, lo.y + (hi.y - lo.y) * j as f64 / n as f64);
            let d = surface.d2(u, v);
            let (lu, lv) = (d.du.length(), d.dv.length());
            speed_u += lu;
            speed_v += lv;
            count += 1.0;
            // Bending: the second derivatives across the tangent plane.
            let normal = d.du.cross(d.dv);
            let nl = normal.length();
            if nl > 0.0 {
                let nn = normal * (1.0 / nl);
                bend_u = bend_u.max(d.duu.dot(nn).abs());
                bend_v = bend_v.max(d.dvv.dot(nn).abs());
                if lu > 0.0 {
                    turn_u = turn_u.max(d.duu.dot(nn).abs() / lu);
                }
                if lv > 0.0 {
                    turn_v = turn_v.max(d.dvv.dot(nn).abs() / lv);
                }
            }
        }
    }
    let (su, sv) = ((speed_u / count).max(1e-12), (speed_v / count).max(1e-12));
    // A step h along u deviates by about bend h² / 8 and turns by turn·h.
    let step = |bend: f64, turn: f64, span: f64| -> usize {
        if span <= 0.0 || (bend <= 1e-12 && turn <= 1e-12) {
            return 1;
        }
        let by_sag = if bend > 0.0 { (8.0 * options.deflection / bend).sqrt() } else { f64::INFINITY };
        let by_turn = if turn > 0.0 { options.angle / turn } else { f64::INFINITY };
        ((span / by_sag.min(by_turn)).ceil() as usize).clamp(1, 400)
    };
    let (mut nu, mut nv) = (step(bend_u, turn_u, hi.x - lo.x), step(bend_v, turn_v, hi.y - lo.y));
    // A grid of a single step one way has no inner points (its ends lie on the boundary),
    // and the triangles would then span the face from side to side: at least one middle row.
    if nu > 1 {
        nv = nv.max(2);
    }
    if nv > 1 {
        nu = nu.max(2);
    }
    ((su, sv), (nu, nv))
}

/// Even-odd rule over all loops.
fn inside_polygons(polygons: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for poly in polygons {
        let n = poly.len();
        for i in 0..n {
            let (a, b) = (poly[i], poly[(i + 1) % n]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > p.x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

fn distance_to_polygons(polygons: &[Vec<Point2>], p: Point2) -> f64 {
    let mut best = f64::INFINITY;
    for poly in polygons {
        let n = poly.len();
        for i in 0..n {
            let (a, b) = (poly[i], poly[(i + 1) % n]);
            let ab = b - a;
            let len2 = ab.dot(ab);
            let t = if len2 > 0.0 { ((p - a).dot(ab) / len2).clamp(0.0, 1.0) } else { 0.0 };
            best = best.min(p.distance(a + ab * t));
        }
    }
    best
}

// --- constrained Delaunay triangulation -------------------------------------------------

/// Twice the signed area of `abc` (positive counter-clockwise).
fn orient(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// True if `d` lies inside the circumcircle of the counter-clockwise triangle `abc`.
fn in_circle(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let (ax, ay) = (a.x - d.x, a.y - d.y);
    let (bx, by) = (b.x - d.x, b.y - d.y);
    let (cx, cy) = (c.x - d.x, c.y - d.y);
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay) + (cx * cx + cy * cy) * (ax * by - bx * ay);
    det > 0.0
}

/// Triangles (counter-clockwise) over `points` containing every segment (as an edge, or a
/// chain of edges through points lying on it). Coincident points are triangulated once.
fn constrained_delaunay(points: &[Point2], segments: &[(usize, usize)]) -> Vec<[usize; 3]> {
    let n = points.len();
    if n < 3 {
        return Vec::new();
    }
    let (mut lo, mut hi) = (points[0], points[0]);
    for p in points {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let size = (hi.x - lo.x).max(hi.y - lo.y).max(1e-12);
    let eps = 1e-12 * size;
    let mid = Point2::new(0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y));
    let mut pts: Vec<Point2> = points.to_vec();
    // A super triangle.
    pts.push(Point2::new(mid.x - 20.0 * size, mid.y - 10.0 * size));
    pts.push(Point2::new(mid.x + 20.0 * size, mid.y - 10.0 * size));
    pts.push(Point2::new(mid.x, mid.y + 20.0 * size));
    let mut tris = Triangulation::new(&pts, [n, n + 1, n + 2]);
    // Coincident points are one.
    let mut same: Vec<usize> = (0..n).collect();
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::default();
    let cell = (size * 1e-9).max(1e-300);
    for i in 0..n {
        let key = ((pts[i].x / cell).round() as i64, (pts[i].y / cell).round() as i64);
        let mut found = None;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(list) = grid.get(&(key.0 + dx, key.1 + dy)) {
                    for &j in list {
                        if pts[j].distance(pts[i]) <= eps {
                            found = Some(j);
                            break 'search;
                        }
                    }
                }
            }
        }
        match found {
            Some(j) => same[i] = j,
            None => {
                grid.entry(key).or_default().push(i);
                tris.insert(i);
            }
        }
    }
    for &(a, b) in segments {
        let (a, b) = (same[a], same[b]);
        if a != b {
            tris.force_edge(a, b, eps);
        }
    }
    tris.live()
        .filter(|t| t.iter().all(|&v| v < n))
        .collect()
}

/// A triangulation kept as triangles and a map from directed edges to triangles.
struct Triangulation<'a> {
    pts: &'a [Point2],
    tris: Vec<Option<[usize; 3]>>,
    edges: HashMap<(usize, usize), usize>,
    last: usize,
}

impl<'a> Triangulation<'a> {
    fn new(pts: &'a [Point2], first: [usize; 3]) -> Self {
        let mut t = Triangulation { pts, tris: Vec::new(), edges: HashMap::default(), last: 0 };
        t.add(first);
        t
    }

    fn add(&mut self, t: [usize; 3]) -> usize {
        let i = self.tris.len();
        self.tris.push(Some(t));
        for k in 0..3 {
            self.edges.insert((t[k], t[(k + 1) % 3]), i);
        }
        self.last = i;
        i
    }

    fn remove(&mut self, i: usize) -> [usize; 3] {
        let t = self.tris[i].take().expect("a live triangle");
        for k in 0..3 {
            if self.edges.get(&(t[k], t[(k + 1) % 3])) == Some(&i) {
                self.edges.remove(&(t[k], t[(k + 1) % 3]));
            }
        }
        t
    }

    fn live(&self) -> impl Iterator<Item = [usize; 3]> + '_ {
        self.tris.iter().flatten().copied()
    }

    /// The triangle containing `p`, by walking from the last one made.
    fn locate(&self, p: Point2) -> Option<usize> {
        let mut i = if self.tris[self.last].is_some() { self.last } else { self.tris.iter().position(Option::is_some)? };
        for _ in 0..self.tris.len() + 10 {
            let t = self.tris[i]?;
            let mut moved = false;
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if orient(self.pts[a], self.pts[b], p) < 0.0 {
                    if let Some(&j) = self.edges.get(&(b, a)) {
                        i = j;
                        moved = true;
                        break;
                    }
                }
            }
            if !moved {
                return Some(i);
            }
        }
        // Walking failed (degenerate configurations): search them all.
        self.tris.iter().position(|t| {
            t.is_some_and(|t| (0..3).all(|k| orient(self.pts[t[k]], self.pts[t[(k + 1) % 3]], p) >= 0.0))
        })
    }

    /// Bowyer-Watson insertion of point `p`.
    fn insert(&mut self, p: usize) {
        let Some(start) = self.locate(self.pts[p]) else { return };
        // The cavity: triangles whose circumcircle holds the point, grown from the one
        // containing it.
        let mut bad = vec![start];
        let mut seen = FastSet::from_iter([start]);
        let mut k = 0;
        while k < bad.len() {
            let t = self.tris[bad[k]].expect("live");
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                if let Some(&j) = self.edges.get(&(b, a)) {
                    if seen.insert(j) {
                        let u = self.tris[j].expect("live");
                        if in_circle(self.pts[u[0]], self.pts[u[1]], self.pts[u[2]], self.pts[p]) {
                            bad.push(j);
                        }
                    }
                }
            }
            k += 1;
        }
        let removed: Vec<[usize; 3]> = bad.iter().map(|&i| self.remove(i)).collect();
        let inner: FastSet<(usize, usize)> =
            removed.iter().flat_map(|t| (0..3).map(move |k| (t[k], t[(k + 1) % 3]))).collect();
        for t in &removed {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if !inner.contains(&(b, a)) {
                    self.add([a, b, p]);
                }
            }
        }
    }

    /// Makes `a`-`b` an edge, flipping the edges crossing it; through points lying on the
    /// segment it is forced piecewise.
    fn force_edge(&mut self, a: usize, b: usize, eps: f64) {
        let mut pending = vec![(a, b)];
        let mut guard = 0;
        while let Some((a, b)) = pending.pop() {
            guard += 1;
            if guard > 10_000 {
                return;
            }
            if self.edges.contains_key(&(a, b)) || self.edges.contains_key(&(b, a)) {
                continue;
            }
            let (pa, pb) = (self.pts[a], self.pts[b]);
            let length = pa.distance(pb);
            // A point on the segment splits it.
            let on = self.tris.iter().flatten().flatten().copied().filter(|&v| v != a && v != b).find(|&v| {
                let pv = self.pts[v];
                let t = (pv - pa).dot(pb - pa) / (length * length);
                t > 0.0 && t < 1.0 && orient(pa, pb, pv).abs() / length <= eps.max(1e-9 * length)
            });
            if let Some(v) = on {
                pending.push((v, b));
                pending.push((a, v));
                continue;
            }
            // Sloan's method: the edges crossing the segment are flipped in turn (those whose
            // quadrilateral is not convex wait for their neighbours), new edges still
            // crossing it going back in the queue.
            let mut queue: std::collections::VecDeque<(usize, usize)> = self
                .tris
                .iter()
                .flatten()
                .flat_map(|t| (0..3).map(move |k| (t[k], t[(k + 1) % 3])))
                .filter(|&(u, v)| u < v && u != a && u != b && v != a && v != b)
                .filter(|&(u, v)| self.crosses(u, v, a, b))
                .collect();
            let mut budget = 50 * (queue.len() + 10);
            while let Some((u, v)) = queue.pop_front() {
                budget -= 1;
                if budget == 0 {
                    break;
                }
                match self.flip(u, v) {
                    Some((w, x)) => {
                        if w != a && w != b && x != a && x != b && self.crosses(w, x, a, b) {
                            queue.push_back((w, x));
                        }
                    }
                    None if self.edges.contains_key(&(u, v)) => queue.push_back((u, v)),
                    None => {}
                }
            }
        }
    }

    fn crosses(&self, u: usize, v: usize, a: usize, b: usize) -> bool {
        let (pu, pv, pa, pb) = (self.pts[u], self.pts[v], self.pts[a], self.pts[b]);
        orient(pa, pb, pu) * orient(pa, pb, pv) < 0.0 && orient(pu, pv, pa) * orient(pu, pv, pb) < 0.0
    }

    /// Flips edge `u`-`v` (either way round) if the two triangles beside it make a convex
    /// quadrilateral; the new edge.
    fn flip(&mut self, u: usize, v: usize) -> Option<(usize, usize)> {
        let (i, j) = (*self.edges.get(&(u, v))?, *self.edges.get(&(v, u))?);
        let (t1, t2) = (self.tris[i]?, self.tris[j]?);
        let other = |t: [usize; 3]| t.into_iter().find(|&x| x != u && x != v).expect("a third vertex");
        let (w, x) = (other(t1), other(t2));
        // t1 = (u, v, w) counter-clockwise, t2 = (v, u, x).
        let (pu, pv, pw, px) = (self.pts[u], self.pts[v], self.pts[w], self.pts[x]);
        if orient(pw, px, pv) <= 0.0 || orient(px, pw, pu) <= 0.0 {
            return None;
        }
        self.remove(i);
        self.remove(j);
        self.add([w, x, v]);
        self.add([x, w, u]);
        Some((w.min(x), w.max(x)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Cylinder, Sphere};
    use crate::properties;

    fn area_and_volume(nodes: &[Point], triangles: &[[u32; 3]]) -> (f64, f64) {
        let (mut area, mut volume) = (0.0, 0.0);
        for t in triangles {
            let [a, b, c] = t.map(|i| nodes[i as usize]);
            let (ab, ac) = (b - a, c - a);
            area += 0.5 * ab.cross(ac).length();
            volume += (a - Point::ORIGIN).dot((b - Point::ORIGIN).cross(c - Point::ORIGIN)) / 6.0;
        }
        (area, volume)
    }

    /// Closed: every undirected edge used twice, once each way.
    fn closed(triangles: &[[u32; 3]]) -> bool {
        let mut count: HashMap<(u32, u32), i32> = HashMap::default();
        for t in triangles {
            for k in 0..3 {
                *count.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
            }
        }
        count.iter().all(|(&(a, b), &n)| n == 1 && count.get(&(b, a)) == Some(&1))
    }

    #[test]
    fn a_box_is_twelve_triangles() {
        let mut model = Model::new();
        let b = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let (nodes, triangles) = triangulate(&mut model, b, MeshOptions::default()).unwrap();
        assert_eq!((nodes.len(), triangles.len()), (8, 12));
        assert!(closed(&triangles));
        let (area, volume) = area_and_volume(&nodes, &triangles);
        assert!((area - 22.0).abs() < 1e-12 && (volume - 6.0).abs() < 1e-12);
    }

    #[test]
    fn curved_solids_mesh_closed_and_close() {
        let mut model = Model::new();
        for (shape, deflection) in [
            (Cylinder::new(2.0, 3.0).build(&mut model).unwrap(), 1e-3),
            (Sphere::new(1.5).build(&mut model).unwrap(), 1e-3),
        ] {
            let (nodes, triangles) = triangulate(&mut model, shape, MeshOptions { deflection, angle: 0.3 }).unwrap();
            assert!(closed(&triangles), "not closed");
            let (_, volume) = area_and_volume(&nodes, &triangles);
            let exact = properties::volume(&model, shape).mass;
            assert!((volume - exact).abs() < 0.01 * exact, "volume {volume} against {exact}");
            assert!(volume < exact + 1e-9, "an inscribed mesh is smaller");
        }
    }

    #[test]
    fn a_tall_cylinder_keeps_within_the_deflection_and_records_it() {
        let mut model = Model::new();
        let cylinder = Cylinder::new(4.0, 6.0).build(&mut model).unwrap();
        mesh_with(&mut model, cylinder, MeshOptions { deflection: 0.02, angle: 0.2 }).unwrap();
        for face in model.explore(cylinder, NodeKind::Face).collect::<Vec<_>>() {
            let mesh = model.face_mesh(face).unwrap();
            let Some(Surface::Cylinder(_)) = model.surface(face) else {
                assert_eq!(mesh.deflection(), 0.0, "a plane's mesh is exact");
                continue;
            };
            // Every side's middle lies within the deflection of the radius.
            let nodes = mesh.nodes();
            let mut worst = 0.0f64;
            for t in mesh.triangles() {
                for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    let (p, q) = (nodes[a as usize], nodes[b as usize]);
                    worst = worst.max(4.0 - (0.5 * (p.x + q.x)).hypot(0.5 * (p.y + q.y)));
                }
            }
            assert!(worst <= 0.02, "sags by {worst}");
            assert!(mesh.deflection() > 0.0 && mesh.deflection() <= 0.02 && mesh.deflection() >= worst - 1e-12);
        }
    }
}
