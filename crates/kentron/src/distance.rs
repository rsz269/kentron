//! The distance between two shapes.
//!
//! Both shapes are meshed (finely, for their size); the closest pair of triangles is found
//! through bounding volume hierarchies, and the two points refined on the faces' exact
//! surfaces. A shape lying inside a solid of the other is at distance zero.

use crate::boolean::{PointClassifier, PointLocation};
use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::geometry::surface::Surface;
use std::collections::{HashMap, HashSet};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{BoundingBox, Bvh, Point, Transform, Transformable, Vector};
use crate::meshing::{self, MeshOptions};
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// The smallest distance between two shapes, and a point on each where it is reached.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShapeDistance {
    pub value: f64,
    pub on_first: Point,
    pub on_second: Point,
}

/// The smallest distance between `a` and `b` (zero where they touch, cross, or one lies
/// inside a solid of the other).
#[doc(alias = "BRepExtrema_DistShapeShape")]
#[doc(alias = "distmini")]
pub fn distance(model: &mut Model, a: Shape, b: Shape) -> Result<ShapeDistance> {
    distance_within(model, a, b, f64::INFINITY)
}

/// [`distance`] where it is at most `limit`; beyond it, some value above `limit` (the
/// shapes' meshes showing them that far apart, not refined on their exact geometry).
pub fn distance_within(model: &mut Model, a: Shape, b: Shape, limit: f64) -> Result<ShapeDistance> {

    // What depends on each shape alone (its size, its points' classification, its
    // triangles) is kept in the model with it: shapes are often measured against several.
    let (print_a, print_b) = (model.shape_fingerprint(a), model.shape_fingerprint(b));
    let print = |s: Shape| if s == a { print_a } else { print_b };
    // Wholly inside a solid of the other: zero. (A shape crossing the other's boundary is
    // caught by the triangles below; one wholly inside has all its vertices inside, so a
    // sample of them must agree, not one point a ray may misjudge.)
    for (outer, inner) in [(a, b), (b, a)] {
        if model.explore(outer, NodeKind::Solid).next().is_none() {
            continue;
        }
        let vertices: Vec<Shape> = model.explore(inner, NodeKind::Vertex).collect();
        if vertices.is_empty() {
            continue;
        }
        let step = (vertices.len() / 8).max(1);
        let sample: Vec<Point> = vertices.iter().step_by(step).take(8).map(|&v| model.point(v)).collect();
        let kept = model.derived(outer.node, print(outer), || KeptClassifier(PointClassifier::new(model, outer)));
        let classifier = kept.0.as_ref().map_err(Clone::clone)?;
        let all_inside = sample.iter().all(|&p| classifier.locate(p) != PointLocation::Outside);
        if all_inside {
            return Ok(ShapeDistance { value: 0.0, on_first: sample[0], on_second: sample[0] });
        }
    }
    // A vertex and faces on elementary surfaces, or edges: exactly, without meshes.
    let exactly = |model: &Model, v: Shape, s: Shape| vertex_to_faces(model, v, s).or_else(|| vertex_to_edges(model, v, s));
    if let Some((value, near)) = exactly(model, a, b) {
        return Ok(ShapeDistance { value, on_first: model.point(a), on_second: near });
    }
    if let Some((value, near)) = exactly(model, b, a) {
        return Ok(ShapeDistance { value, on_first: near, on_second: model.point(b) });
    }
    let size = |s: Shape, model: &Model| {
        model
            .derived(s.node, print(s), || {
                let b = crate::properties::bounding_box(model, s);
                KeptSize(if b.is_empty() { 0.0 } else { b.min.distance(b.max) })
            })
            .0
    };
    // Each shape is meshed finely for its own size (a point has none).
    let (size_a, size_b) = (size(a, model), size(b, model));
    let scale = size_a.max(size_b).max(1e-9);
    let options = |size: f64| MeshOptions { deflection: 1e-3 * if size > 0.0 { size } else { scale }, angle: 0.3, relative: false };
    let (options_a, options_b) = (options(size_a), options(size_b));
    let kept_a = Triangles::kept(model, a, options_a, print_a)?;
    let kept_b = Triangles::kept(model, b, options_b, print_b)?;
    let (ta, tb) = (&kept_a.0, &kept_b.0);

    // Points (vertices, edges) without faces.
    let (Some(bva), Some(bvb)) = (kept_a.1.as_ref(), kept_b.1.as_ref()) else {
        return Ok(ShapeDistance { value: f64::INFINITY, on_first: Point::ORIGIN, on_second: Point::ORIGIN });
    };
    let Some((_, _, nearest)) = bva.closest_pair(bvb, |i, j| triangle_triangle(ta.corners(i), tb.corners(j)).0) else {
        return Ok(ShapeDistance { value: f64::INFINITY, on_first: Point::ORIGIN, on_second: Point::ORIGIN });
    };
    // The meshes stray from the exact geometry (a concave face's towards the other shape, so
    // that meshes can even cross across a small gap): every pair of pieces the meshes bring
    // within their deviations of the nearest is refined on the exact geometry.
    let slack = ta.deviation.iter().fold(0.0f64, |m, &d| m.max(d)) + tb.deviation.iter().fold(0.0f64, |m, &d| m.max(d));
    let below = nearest - 2.0 * slack - 1e-12 * scale;
    if below > limit {
        return Ok(ShapeDistance { value: below, on_first: Point::ORIGIN, on_second: Point::ORIGIN });
    }
    // Meshes within their deviations of each other (touching, or crossing): where a point of
    // one's boundary lies inside a solid of the other (clear of its boundary), the solids
    // overlap, at distance zero. The nearest points on the exact geometry may not tell so:
    // where their faces cross at a shallow angle, the meshes cross over a long stretch about
    // the crossing, and the exact points found near where they do may still lie apart by
    // the angle times the way left to it.
    if nearest <= 2.0 * slack {
        if let Some(p) = overlap_point(model, [a, b], [print_a, print_b], [ta, tb])? {
            return Ok(ShapeDistance { value: 0.0, on_first: p, on_second: p });
        }
    }
    let threshold = nearest + 2.0 * slack + 1e-12 * scale;
    // (Per pair of pieces, its nearest pair of triangles: the first by their indices among
    // equally near ones, so that the points do not depend on the order pairs are met in.)
    let mut candidates: crate::hashing::FastMap<(usize, usize), (f64, (usize, usize), Point, Point)> = Default::default();
    for (i, j) in bva.pairs_within(bvb, threshold) {
        let (d, p, q) = triangle_triangle(ta.corners(i), tb.corners(j));
        if d <= threshold {
            let entry = candidates.entry((ta.face[i], tb.face[j])).or_insert((d, (i, j), p, q));
            if d < entry.0 || (d == entry.0 && (i, j) < entry.1) {
                *entry = (d, (i, j), p, q);
            }
        }
    }
    let mut candidates: Vec<((usize, usize), (f64, Point, Point))> = candidates.into_iter().map(|(k, (d, _, p, q))| (k, (d, p, q))).collect();
    // (Equally near pairs by their pieces, so that the points given are the same each time.)
    candidates.sort_by(|x, y| x.1 .0.total_cmp(&y.1 .0).then(x.0.cmp(&y.0)));
    let reach = 20.0 * options_a.deflection.max(options_b.deflection);
    let mut best: Option<(f64, Point, Point)> = None;
    for ((ka, kb), (d, p, q)) in candidates {
        if best.is_some_and(|b| d - 2.0 * slack > b.0) {
            break;
        }
        // On the exact geometry: the pieces, and the faces' boundary edges near the other
        // piece, each pair kept when its points are surely on the pieces (not out on a
        // surface beyond a face's boundary, nor by it: there the edges' answer counts).
        // Without any, the meshes' estimate.
        let within = d + 2.0 * slack + reach;
        let (near_a, near_b) = (ta.edges_near(ka, q, within), tb.edges_near(kb, p, within));
        let model: &Model = model;
        let mut pairs: Vec<(&Exact, &Exact, bool, bool)> = vec![(ta.exact[ka].get(model), tb.exact[kb].get(model), true, true)];
        pairs.extend(near_a.iter().map(|&e| (ta.edges[e].0.get(model), tb.exact[kb].get(model), false, true)));
        pairs.extend(near_b.iter().map(|&e| (ta.exact[ka].get(model), tb.edges[e].0.get(model), true, false)));
        for &ea in &near_a {
            pairs.extend(near_b.iter().map(|&eb| (ta.edges[ea].0.get(model), tb.edges[eb].0.get(model), false, false)));
        }
        let mut found: Option<(f64, Point, Point)> = None;
        for (ga, gb, check_a, check_b) in pairs {
            let Some(r) = refine(ga, gb, p, q, within.max(reach), scale) else { continue };
            if (check_a && !ta.surely_on_piece(model, ka, r.1, scale)) || (check_b && !tb.surely_on_piece(model, kb, r.2, scale)) {
                continue;
            }
            if found.is_none_or(|f| r.0 < f.0) {
                found = Some(r);
            }
        }
        let found = found.unwrap_or((d, p, q));
        if best.is_none_or(|b| found.0 < b.0) {
            best = Some(found);
        }
    }
    let (value, on_first, on_second) = best.expect("the nearest pair is a candidate");
    Ok(ShapeDistance { value, on_first, on_second })
}

/// A point of either shape's boundary lying inside a solid of the other (not on its boundary):
/// among the nodes of its triangles in the box the two shapes' triangles share, a spread of at
/// most [`OVERLAP_SAMPLES`]. `None` where none of those does (or neither has a solid).
fn overlap_point(model: &Model, shapes: [Shape; 2], prints: [u64; 2], triangles: [&Triangles; 2]) -> Result<Option<Point>> {
    let boxes = triangles.map(|t| {
        let mut b = BoundingBox::EMPTY;
        for &p in &t.nodes {
            b.add_point(p);
        }
        b
    });
    let shared = |p: Point| (0..2).all(|k| boxes[k].contains(p));
    for (outer, inner) in [(0, 1), (1, 0)] {
        if model.explore(shapes[outer], NodeKind::Solid).next().is_none() {
            continue;
        }
        let nodes: Vec<Point> = triangles[inner].nodes.iter().copied().filter(|&p| shared(p)).collect();
        if nodes.is_empty() {
            continue;
        }
        let kept = model.derived(shapes[outer].node, prints[outer], || KeptClassifier(PointClassifier::new(model, shapes[outer])));
        let classifier = kept.0.as_ref().map_err(Clone::clone)?;
        let step = nodes.len().div_ceil(OVERLAP_SAMPLES).max(1);
        if let Some(&p) = nodes.iter().step_by(step).find(|&&p| classifier.locate(p) == PointLocation::Inside) {
            return Ok(Some(p));
        }
    }
    Ok(None)
}

/// How many points of a shape's boundary [`overlap_point`] tries.
const OVERLAP_SAMPLES: usize = 64;

/// The distance from `vertex` to `faces` (a face, or a shell or compound holding faces only,
/// none of them placed) on planes, cylinders, cones, spheres and tori, and the point of them
/// nearest: as `BRepExtrema_DistShapeShape` finds it between a vertex and a face, the
/// nearest of the point's extrema on each face's surface that lie inside the face
/// (`BRepExtrema_ExtPF`, classified as `BRepClass_FaceClassifier` does) and of its nearest
/// points on the faces' boundary edges (`BRepExtrema_ExtPC`, their ends included). `None`
/// where it is not told so: other shapes or surfaces, edges on other curves than lines and
/// circles (see [`edge_nearest`]), a point on a surface's axis or at its centre (its extrema
/// not isolated there), or an extremum on a face's boundary (left to the meshes and the
/// refinement on the exact geometry).
#[doc(alias = "BRepExtrema_ExtPF")]
fn vertex_to_faces(model: &Model, vertex: Shape, faces: Shape) -> Option<(f64, Point)> {
    if vertex.kind() != NodeKind::Vertex {
        return None;
    }
    // (Faces only, under shells and compounds.)
    fn faces_only(model: &Model, s: Shape) -> bool {
        match s.kind() {
            NodeKind::Face => true,
            NodeKind::Shell | NodeKind::Compound => model.children(s).all(|c| faces_only(model, c)),
            _ => false,
        }
    }
    if !faces_only(model, faces) {
        return None;
    }
    let p = model.point(vertex);
    let mut best: Option<(f64, Point)> = None;
    let mut found = false;
    for face in model.explore(faces, NodeKind::Face) {
        found = true;
        if face.placement != Transform::IDENTITY {
            return None;
        }
        let extrema = surface_extrema(model.surface(face)?, p)?;
        // Its boundary edges' nearest points (their ends included).
        let mut nearest: Option<(f64, Point)> = None;
        for edge in model.explore(face, NodeKind::Edge) {
            let q = edge_nearest(model, edge, p)?;
            let d = p.distance(q);
            if nearest.is_none_or(|n| d < n.0) {
                nearest = Some((d, q));
            }
        }
        // The extrema nearer than those, inside the face (on its boundary: not told here).
        // (The face's classifier kept in the model with it, as a solid's is.)
        let mut classifier = None;
        for q in extrema {
            let d = p.distance(q);
            if nearest.is_some_and(|n| d >= n.0) {
                continue;
            }
            let kept = classifier.get_or_insert_with(|| model.derived(face.node, model.shape_fingerprint(face), || KeptClassifier(PointClassifier::new(model, face))));
            let Ok(c) = &kept.0 else { return None };
            match c.locate_on_faces(q) {
                PointLocation::Inside => nearest = Some((d, q)),
                // On the boundary: as near as the boundary's nearest point, within the
                // tolerance, when that is (a sphere's face whose rim lies where the point's
                // foot does); else not told here.
                PointLocation::Boundary => {
                    if nearest.is_some_and(|n| n.0 - d <= crate::tolerance::LINEAR) {
                        continue;
                    }
                    return None;
                }
                PointLocation::Outside => {}
            }
        }
        let Some(nearest) = nearest else { return None };
        if best.is_none_or(|b| nearest.0 < b.0) {
            best = Some(nearest);
        }
    }
    best.filter(|_| found)
}

/// The distance from `vertex` to `edges` (an edge, or a wire or compound holding edges only,
/// on lines and circles) and the point of them nearest: the nearest of each edge's nearest
/// point to it, its ends included (`BRepExtrema_ExtPC`), or of its vertex where it has no
/// curve. `None` for other shapes and curves.
fn vertex_to_edges(model: &Model, vertex: Shape, edges: Shape) -> Option<(f64, Point)> {
    if vertex.kind() != NodeKind::Vertex {
        return None;
    }
    fn edges_only(model: &Model, s: Shape) -> bool {
        match s.kind() {
            NodeKind::Edge => true,
            NodeKind::Wire | NodeKind::Compound => model.children(s).all(|c| edges_only(model, c)),
            _ => false,
        }
    }
    if !edges_only(model, edges) {
        return None;
    }
    let p = model.point(vertex);
    let mut best: Option<(f64, Point)> = None;
    for edge in model.explore(edges, NodeKind::Edge) {
        let q = edge_nearest(model, edge, p)?;
        let d = p.distance(q);
        if best.is_none_or(|b| d < b.0) {
            best = Some((d, q));
        }
    }
    best
}

/// The point of `edge` nearest `p`, its ends included (its vertex where it has no curve):
/// on a line or a circle, where the distance along it has one minimum its samples lead to
/// surely; `None` on other curves (a free-form curve's samples may miss its nearest point,
/// which its mesh finds).
fn edge_nearest(model: &Model, edge: Shape, p: Point) -> Option<Point> {
    match model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }) {
        Some((curve @ (Curve::Line(_) | Curve::Circle(_)), first, last)) => Some(crate::geometry::extrema::nearest_on_curve(p, &curve, first, last).point),
        // Other curves: the nearest of the distance's minima over the range, its ends
        // included, all found by the search over the curve's samples (each span's).
        Some((curve, first, last)) => nearest_by_search(&curve, first, last, p),
        None => Some(model.point(model.edge_vertices(edge).0)),
    }
}

/// The point of `curve` over `[first, last]` nearest to `p`: the nearest of the distance's
/// minima there, the range's ends included. `None` where the search fails.
fn nearest_by_search(curve: &Curve, first: f64, last: f64, p: Point) -> Option<Point> {
    use crate::geometry::extrema::{PointCurveSearch, SearchMode};
    let (lo, hi) = (first.min(last), first.max(last));
    let found = PointCurveSearch::with_domain(curve, lo, hi).find_with(p, crate::tolerance::LINEAR, SearchMode::Min, true).ok()?;
    let ends = [lo, hi].map(|t| curve.point(t));
    found.iter().map(|e| e.point).chain(ends).min_by(|a, b| p.distance(*a).total_cmp(&p.distance(*b)))
}

/// The points of `surface` (a plane, cylinder, cone, sphere or torus, whole) where the
/// distance from `p` is stationary: its foot on a plane; on the others, in the plane through
/// the axis (or centre) and `p`, the feet on the meridians there (both sides of the axis).
/// Where `p` lies on the axis or at the centre, the extrema are not isolated: each is a
/// circle round the axis (the whole sphere about its centre), and one point of it stands for
/// it, in the plane through the axis and the frame's x direction. (A face holds such a
/// connected set wholly, or not at all, or meets it on its boundary, whose nearest point is
/// then as near: so one point of it tells, for [`vertex_to_faces`].) `None` for other
/// surfaces, and for a point on a torus's core circle.
fn surface_extrema(surface: &Surface, p: Point) -> Option<Vec<Point>> {
    let tiny = |d: Vector| d.length() <= 1e-12 * (1.0 + (p - Point::ORIGIN).length());
    // `p` from the axis through `origin` along `z`: the unit direction away from it in the
    // plane through the axis and `p`; for `p` on the axis, the frame's x direction (`x`).
    let away = |origin: Point, z: Vector, x: Vector| -> Option<Vector> {
        let d = p - origin;
        let radial = d - z * d.dot(z);
        Some(if tiny(radial) { x } else { radial * (1.0 / radial.length()) })
    };
    Some(match surface {
        Surface::Plane(plane) => {
            let f = plane.frame();
            let n = Vector::from(f.x()).cross(Vector::from(f.y()));
            vec![p - n * (p - f.origin()).dot(n)]
        }
        Surface::Cylinder(c) => {
            let (o, z) = (c.frame().origin(), Vector::from(c.frame().z()));
            let e = away(o, z, Vector::from(c.frame().x()))?;
            let at = o + z * (p - o).dot(z);
            vec![at + e * c.radius(), at - e * c.radius()]
        }
        Surface::Sphere(s) => {
            let d = p - s.frame().origin();
            let e = if tiny(d) { Vector::from(s.frame().x()) } else { d * (1.0 / d.length()) };
            vec![s.frame().origin() + e * s.radius(), s.frame().origin() - e * s.radius()]
        }
        Surface::Cone(c) => {
            // Each meridian line through the plane: from the circle at the frame, along the
            // generator (both nappes).
            let (o, z) = (c.frame().origin(), Vector::from(c.frame().z()));
            let e = away(o, z, Vector::from(c.frame().x()))?;
            let (sin, cos) = c.semi_angle_sin_cos();
            [e, -e]
                .iter()
                .map(|&e| {
                    let (start, along) = (o + e * c.radius(), e * sin + z * cos);
                    start + along * (p - start).dot(along)
                })
                .collect()
        }
        Surface::Torus(t) => {
            let (o, z) = (t.frame().origin(), Vector::from(t.frame().z()));
            let e = away(o, z, Vector::from(t.frame().x()))?;
            let mut found = Vec::with_capacity(4);
            for centre in [o + e * t.major_radius(), o - e * t.major_radius()] {
                let d = p - centre;
                if tiny(d) {
                    return None;
                }
                let d = d * (1.0 / d.length());
                found.extend([centre + d * t.minor_radius(), centre - d * t.minor_radius()]);
            }
            found
        }
        _ => return None,
    })
}

/// A face of a shape near a point, with its point nearest that point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NearFace {
    pub face: Shape,
    pub point: Point,
    pub distance: f64,
}

/// The faces of `shape` within `slack` of being the nearest to `p`, each with its point
/// nearest `p`, nearest first (equally near ones in the order the shape lists its faces):
/// as `BRepExtrema_DistShapeShape` between a vertex and each face finds it, the nearest of
/// the point's projection onto the face's surface where the face holds it
/// (`BRepExtrema_ExtPF`, with `BRepClass_FaceClassifier`) and its projections onto the
/// face's boundary edges (`BRepExtrema_ExtPC`, their ends included). The faces' meshes only
/// pick which faces and edges to project onto, and tell whether a face holds a projection by
/// its boundary; on faces on planes, cylinders, cones, spheres and tori, they are made only
/// where needed for that (see [`nearest_elementary_faces`]). Faces without a surface are
/// left out.
#[doc(alias = "BRepExtrema_ExtPF")]
pub(crate) fn nearest_faces(model: &mut Model, shape: Shape, p: Point, slack: f64) -> Result<Vec<NearFace>> {
    nearest_faces_by(model, shape, p, slack, true)
}

/// [`nearest_faces`]; without `exact_first`, always through the meshes.
fn nearest_faces_by(model: &mut Model, shape: Shape, p: Point, slack: f64, exact_first: bool) -> Result<Vec<NearFace>> {
    let print = model.shape_fingerprint(shape);
    let size = model
        .derived(shape.node, print, || {
            let b = crate::properties::bounding_box(model, shape);
            KeptSize(if b.is_empty() { 0.0 } else { b.min.distance(b.max) })
        })
        .0;
    let scale = size.max(1e-9);
    // (Meshed as finely as for `distance`, so that the meshes serve both.)
    let options = MeshOptions { deflection: 1e-3 * scale, angle: 0.3, relative: false };
    if exact_first {
        if let Some(found) = nearest_elementary_faces(model, shape, print, p, slack, scale, options.deflection)? {
            return Ok(found);
        }
    }
    let kept = Triangles::kept(model, shape, options, print)?;
    let t = &kept.0;
    let Some(bvh) = kept.1.as_ref() else { return Ok(Vec::new()) };
    let Some((_, nearest)) = bvh.nearest(p, |i| point_triangle(p, t.corners(i)).distance(p)) else { return Ok(Vec::new()) };
    // A face within `slack` of the nearest has its mesh within the meshes' deviations more.
    let deviation = t.deviation.iter().fold(0.0f64, |m, &d| m.max(d));
    let reach = nearest + 2.0 * deviation + slack + 1e-12 * scale;
    let query = BoundingBox::new(p - Vector::new(reach, reach, reach), p + Vector::new(reach, reach, reach));
    let mut pieces: Vec<usize> = bvh.overlapping(&query).into_iter().filter(|&i| point_triangle(p, t.corners(i)).distance(p) <= reach).map(|i| t.face[i]).collect();
    pieces.sort_unstable();
    pieces.dedup();
    let mut found = Vec::new();
    for k in pieces {
        let Source::Face(face) = t.exact[k].source else { continue };
        let model: &Model = model;
        let mut best: Option<(f64, Point)> = None;
        let mut keep = |d: f64, q: Point| {
            if best.is_none_or(|b| d < b.0) {
                best = Some((d, q));
            }
        };
        // The projection onto the surface, where the face holds it: surely, on its mesh and
        // clear of its boundary; else as the face's classifier finds it.
        if let Exact::Surface(s) = t.exact[k].get(model) {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, s, None);
            let q = s.point(u, v);
            let on = t.surely_on_piece(model, k, q, scale)
                || (t.near_piece(k, q, scale) && {
                    let kept = model.derived(face.node, model.shape_fingerprint(face), || KeptClassifier(PointClassifier::new(model, face)));
                    let classifier = kept.0.as_ref().map_err(Clone::clone)?;
                    classifier.locate(q) == PointLocation::Boundary
                });
            if on {
                keep(q.distance(p), q);
            }
        }
        // The boundary edges near.
        for &e in &t.piece_edges[k] {
            let (lazy, points, deviation) = &t.edges[e];
            if polyline_distance(points, p) > reach + 2.0 * deviation {
                continue;
            }
            if let Exact::Curve(c, first, last) = lazy.get(model) {
                let q = crate::geometry::extrema::nearest_on_curve(p, c, *first, *last).point;
                keep(q.distance(p), q);
            }
        }
        if let Some((distance, point)) = best {
            found.push(NearFace { face, point, distance });
        }
    }
    // (Stable: equally near faces stay in the shape's order.)
    found.sort_by(|a, b| a.distance.total_cmp(&b.distance));
    if let Some(least) = found.first().map(|f| f.distance) {
        found.retain(|f| f.distance <= least + slack);
    }
    Ok(found)
}

/// A shape's faces on planes, cylinders, cones, spheres and tori, with their boundary edges
/// and boxes, kept in the model with it (`None` where a face lies on another surface).
struct KeptElementaryFaces(Option<(Vec<ElementaryFace>, Vec<BoundaryEdge>)>);

struct ElementaryFace {
    face: Shape,
    /// Its surface, placed.
    surface: Surface,
    /// A box around it (its vertices, edges and surface).
    bounds: BoundingBox,
    /// Its boundary edges with curves (indices into the shape's edges), each once.
    edges: Vec<usize>,
    /// Those of them it bounds on one side only (not its seams).
    rims: Vec<usize>,
}

struct BoundaryEdge {
    curve: Curve,
    first: f64,
    last: f64,
    bounds: BoundingBox,
}

impl KeptElementaryFaces {
    fn gather(model: &Model, shape: Shape) -> KeptElementaryFaces {
        let mut faces = Vec::new();
        let mut edges: Vec<BoundaryEdge> = Vec::new();
        let mut known: HashMap<(crate::topology::NodeId, [[u64; 4]; 3]), usize> = HashMap::new();
        for face in model.explore(shape, NodeKind::Face) {
            let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
            let Some(surface) = model.surface(local) else { continue };
            if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) {
                return KeptElementaryFaces(None);
            }
            let surface = surface.transformed(&face.placement);
            // (Each edge with how many times the face uses it: twice, a seam.)
            let mut boundary: Vec<(usize, usize)> = Vec::new();
            for edge in model.explore(face, NodeKind::Edge) {
                let rows = edge.placement.rows().map(|r| r.map(f64::to_bits));
                let i = match known.get(&(edge.node, rows)) {
                    Some(&i) => i,
                    None => {
                        let Some((curve, first, last)) = model.edge_curve(edge) else { continue };
                        known.insert((edge.node, rows), edges.len());
                        edges.push(BoundaryEdge { curve, first, last, bounds: crate::properties::bounding_box(model, edge) });
                        edges.len() - 1
                    }
                };
                match boundary.iter_mut().find(|(e, _)| *e == i) {
                    Some((_, uses)) => *uses += 1,
                    None => boundary.push((i, 1)),
                }
            }
            faces.push(ElementaryFace {
                face,
                surface,
                bounds: crate::properties::bounding_box(model, face),
                edges: boundary.iter().map(|&(e, _)| e).collect(),
                rims: boundary.iter().filter(|&&(_, uses)| uses == 1).map(|&(e, _)| e).collect(),
            });
        }
        KeptElementaryFaces(Some((faces, edges)))
    }
}

/// [`nearest_faces`] on a shape whose faces all lie on planes, cylinders, cones, spheres and
/// tori (`None` for other shapes, and where a face's classifier cannot be made), on the
/// faces' exact geometry alone: the point's projections onto those surfaces are found in
/// closed form, so no mesh is needed to pick the faces and edges to project onto; the
/// faces' and edges' boxes (and the distance to each face's whole surface) rule out those
/// that cannot be near enough, and the face's classifier tells whether it holds the
/// projection onto its surface.
fn nearest_elementary_faces(model: &mut Model, shape: Shape, print: u64, p: Point, slack: f64, scale: f64, deflection: f64) -> Result<Option<Vec<NearFace>>> {
    let model: &Model = model;
    let kept = model.derived(shape.node, print, || KeptElementaryFaces::gather(model, shape));
    let Some((faces, edges)) = &kept.0 else { return Ok(None) };
    // (The boxes sample curved edges and surfaces, with margins: by a little more still,
    // as the meshes' deviations are allowed for.)
    let margin = 1e-3 * scale;
    // How far from a face's rims its mesh may tell otherwise than its classifier whether
    // it holds a point: the rims' polylines and its mesh each deviate by up to
    // `deflection`, and a point within twice its mesh's deviation of it counts as on it.
    let band = 4.0 * deflection + 1e-9 * scale;
    let below = |b: &BoundingBox| if b.is_empty() { 0.0 } else { (b.distance_to_point(p) - margin).max(0.0) };
    let mut order: Vec<(f64, usize)> = faces.iter().enumerate().map(|(i, f)| (below(&f.bounds), i)).collect();
    order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut least = f64::INFINITY;
    let mut found: Vec<(usize, NearFace)> = Vec::new();
    for (lower, i) in order {
        // No face further on comes within `slack` of the nearest.
        if lower > least + slack {
            break;
        }
        let f = &faces[i];
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &f.surface, None);
        let q = f.surface.point(u, v);
        // None of the face (on its surface) is nearer than its surface.
        if q.distance(p) > least + slack + 1e-9 * scale {
            continue;
        }
        // Whether the face holds the projection onto its surface: inside it, as its
        // classifier finds it. By its boundary (within the classifier's tolerance), and by
        // its rims within what the meshes may deviate from them, it is as its mesh tells
        // (`None`, left open); but a plane's mesh is the polygon of its rims' polylines,
        // whose sides on straight rims are those rims themselves.
        let classifier = model.derived(f.face.node, model.shape_fingerprint(f.face), || KeptClassifier(PointClassifier::new(model, f.face)));
        let Ok(classifier) = &classifier.0 else { return Ok(None) };
        let planar = matches!(f.surface, Surface::Plane(_));
        let by_rim = || {
            f.rims.iter().any(|&e| {
                let edge = &edges[e];
                !(planar && matches!(edge.curve, Curve::Line(_)))
                    && (edge.bounds.is_empty() || edge.bounds.distance_to_point(q) <= band)
                    && crate::geometry::extrema::nearest_on_curve(q, &edge.curve, edge.first, edge.last).point.distance(q) <= band
            })
        };
        let holds = match classifier.locate_on_faces(q) {
            PointLocation::Boundary => None,
            _ if by_rim() => None,
            PointLocation::Inside => Some(true),
            PointLocation::Outside => Some(false),
        };
        if holds == Some(false) && interior_minimum_elsewhere(&f.surface, p, q, scale, |x| classifier.locate_on_faces(x) != PointLocation::Outside) {
            // The face's nearest point may lie inside it away from the surface's nearest:
            // left to the meshes.
            return Ok(None);
        }
        // Its boundary edges' nearest point (the first of the nearest), those that cannot
        // come nearer passed over.
        let mut on_edges: Option<(f64, Point)> = None;
        for &e in &f.edges {
            let edge = &edges[e];
            if on_edges.is_some_and(|b| below(&edge.bounds) >= b.0) {
                continue;
            }
            let x = crate::geometry::extrema::nearest_on_curve(p, &edge.curve, edge.first, edge.last).point;
            let d = x.distance(p);
            if on_edges.is_none_or(|b| d < b.0) {
                on_edges = Some((d, x));
            }
        }
        // The projection where the face holds it, unless an edge comes nearer still. Left
        // open, where the edges' nearest point is the projection but for rounding (as on a
        // crease), as though it held it (the two the same but for rounding); otherwise
        // left to the meshes.
        let holds = match holds {
            Some(holds) => holds,
            None if on_edges.is_some_and(|(_, x)| x.distance(q) <= 1e-12 * scale) => true,
            None => return Ok(None),
        };
        let d = q.distance(p);
        let best = match on_edges {
            Some(e) if !holds || e.0 < d => Some(e),
            _ if holds => Some((d, q)),
            _ => None,
        };
        if let Some((distance, point)) = best {
            least = least.min(distance);
            found.push((i, NearFace { face: f.face, point, distance }));
        }
    }
    // In the shape's order, then (stably) nearest first.
    found.sort_by_key(|&(i, _)| i);
    let mut found: Vec<NearFace> = found.into_iter().map(|(_, f)| f).collect();
    found.sort_by(|a, b| a.distance.total_cmp(&b.distance));
    if let Some(least) = found.first().map(|f| f.distance) {
        found.retain(|f| f.distance <= least + slack);
    }
    Ok(Some(found))
}

/// Whether a face on `surface` (an elementary one) whose `holds` tells its points may have
/// its point nearest `p` inside it other than `q`, the surface's point nearest `p` (which
/// the face does not hold). Off the surface's nearest point, its distance from `p` has a
/// local minimum only on a cone, at the foot on the other generator on `p`'s side of the
/// axis (the stationary points on the far side, and a sphere's, cylinder's or torus's
/// other ones, are passes or maxima, going round the axis or centre towards `p`); and
/// where `p` lies on the axis or at the centre, the nearest points are not isolated.
fn interior_minimum_elsewhere(surface: &Surface, p: Point, q: Point, scale: f64, holds: impl Fn(Point) -> bool) -> bool {
    let Some(extrema) = surface_extrema(surface, p) else { return true };
    let Surface::Cone(c) = surface else { return false };
    let (o, z) = (c.frame().origin(), Vector::from(c.frame().z()));
    let d = p - o;
    let away = d - z * d.dot(z);
    extrema.into_iter().any(|x| (x - o).dot(away) > 0.0 && x.distance(q) > 1e-9 * scale && holds(x))
}

/// The nearest points of two pieces' exact geometry, from the meshes' nearest points `p` and
/// `q`, staying within `reach` of them (none when a piece has no exact geometry).
fn refine(ea: &Exact, eb: &Exact, p: Point, q: Point, reach: f64, scale: f64) -> Option<(f64, Point, Point)> {
    let mut best: Option<(f64, Point, Point)> = None;
    let (mut wa, mut wb) = (None, None);
    // Alternating projections: slow where the pieces touch tangentially (each round gains
    // little), so stopping when they no longer gain.
    let mut x = p;
    for _ in 0..2000 {
        let (y, b) = eb.nearest(x, wb)?;
        let (next, a) = ea.nearest(y, wa)?;
        x = next;
        (wa, wb) = (Some(a), Some(b));
        if x.distance(p) > reach || y.distance(q) > reach {
            break;
        }
        let dist = x.distance(y);
        if best.is_some_and(|b| dist >= b.0) {
            break;
        }
        let stalled = best.is_some_and(|b| dist >= b.0 * (1.0 - 1e-6));
        best = Some((dist, x, y));
        if stalled || dist <= 1e-12 {
            break;
        }
    }
    // Then Levenberg-Marquardt on both pieces' parameters at once: fast where they touch
    // (the residual goes to zero).
    let (Some(a0), Some(b0)) = (wa, wb) else { return best };
    let (na, nb) = (ea.dimension(), eb.dimension());
    let n = na + nb;
    if n == 0 {
        return best;
    }
    let split = |w: &[f64; 4]| ([w[0], if na > 1 { w[1] } else { 0.0 }], [w[na], if nb > 1 { w[na + 1] } else { 0.0 }]);
    let mut w = [0.0; 4];
    w[..na].copy_from_slice(&a0[..na]);
    w[na..n].copy_from_slice(&b0[..nb]);
    let residual = |w: &[f64; 4]| {
        let (u, v) = split(w);
        ea.evaluate(u).0 - eb.evaluate(v).0
    };
    let mut r = residual(&w);
    let mut lambda = 1e-6;
    for _ in 0..100 {
        let (u, v) = split(&w);
        let (_, da) = ea.evaluate(u);
        let (_, db) = eb.evaluate(v);
        let columns: Vec<Vector> = da[..na].iter().copied().chain(db[..nb].iter().map(|&d| -d)).collect();
        let mut m = [[0.0; 4]; 4];
        let mut g = [0.0; 4];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] = if i < n && j < n { columns[i].dot(columns[j]) } else if i == j { 1.0 } else { 0.0 };
            }
            if i < n {
                m[i][i] = m[i][i] * (1.0 + lambda) + 1e-300;
                g[i] = -columns[i].dot(r);
            }
        }
        let Some(step) = solve4(m, g) else { break };
        let mut trial = [w[0] + step[0], w[1] + step[1], w[2] + step[2], w[3] + step[3]];
        let (u, v) = split(&trial);
        let (u, v) = (ea.clamp(u), eb.clamp(v));
        trial[..na].copy_from_slice(&u[..na]);
        trial[na..n].copy_from_slice(&v[..nb]);
        let rt = residual(&trial);
        if rt.length() < r.length() {
            (w, r) = (trial, rt);
            lambda *= 0.3;
            let (u, v) = split(&w);
            let (x, y) = (ea.evaluate(u).0, eb.evaluate(v).0);
            if x.distance(p) > reach || y.distance(q) > reach {
                break;
            }
            if best.is_none_or(|b| r.length() < b.0) {
                best = Some((r.length(), x, y));
            }
            if r.length() <= 1e-13 * (1.0 + scale) {
                break;
            }
        } else {
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
        }
    }
    best
}

/// The solution of a 4 x 4 linear system (Gaussian elimination with partial pivoting).
fn solve4(mut m: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))?;
        if m[pivot][col].abs() <= 1e-300 {
            return None;
        }
        m.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..4 {
            let f = m[row][col] / m[col][col];
            for k in col..4 {
                m[row][k] -= f * m[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0; 4];
    for row in (0..4).rev() {
        let mut sum = b[row];
        for k in row + 1..4 {
            sum -= m[row][k] * x[k];
        }
        x[row] = sum / m[row][row];
    }
    Some(x)
}

/// A shape's faces' triangles, placed, with the face each comes from.
struct Triangles {
    nodes: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    face: Vec<usize>,
    /// Per piece (face, free edge or free vertex), its exact geometry, placed.
    exact: Vec<Lazy>,
    /// Per piece, how far its mesh strays from the exact geometry.
    deviation: Vec<f64>,
    /// Per piece, its triangles (consecutive).
    ranges: Vec<std::ops::Range<usize>>,
    /// The faces' boundary edges: exact curve, polyline and its deviation.
    edges: Vec<(Lazy, Vec<Point>, f64)>,
    /// Per piece, its boundary edges (none but for faces).
    piece_edges: Vec<Vec<usize>>,
}

/// The exact geometry a piece of [`Triangles`] approximates.
enum Exact {
    Surface(Surface),
    Curve(Curve, f64, f64),
    Point(Point),
    /// A face with only a mesh.
    Unknown,
}

/// A shape's size (its box's diagonal), kept in the model with it.
struct KeptSize(f64);

/// A shape's points' classifier (or why it could not be made), kept in the model with it.
struct KeptClassifier(Result<PointClassifier>);

/// A shape's triangles at a mesh deviation, with their box tree, kept in the model with it.
struct KeptTriangles(Triangles, Option<Bvh>);

/// A piece's or edge's exact geometry, made (placed) the first time it is needed: only the
/// few pieces near the nearest pair are refined on theirs.
struct Lazy {
    source: Source,
    made: std::sync::OnceLock<Exact>,
    /// A face's classifier against its bounds (made the first time a point by its boundary
    /// is asked about; none where it could not be made).
    bounds: std::sync::OnceLock<Option<PointClassifier>>,
}

/// What a piece's exact geometry is made from.
enum Source {
    /// A face with a surface, placed.
    Face(Shape),
    /// An edge with a curve, placed.
    Edge(Shape),
    Point(Point),
    Unknown,
}

impl Lazy {
    fn new(source: Source) -> Lazy {
        Lazy { source, made: std::sync::OnceLock::new(), bounds: std::sync::OnceLock::new() }
    }

    /// Whether `p`, a point of a face's surface, lies inside the face's bounds and off its
    /// boundary (by the face's own edges, not its mesh).
    fn inside_bounds(&self, model: &Model, p: Point) -> bool {
        let Source::Face(face) = self.source else { return false };
        let classifier = self.bounds.get_or_init(|| PointClassifier::new(model, face).ok());
        classifier.as_ref().is_some_and(|c| c.locate_on_faces(p) == PointLocation::Inside)
    }

    fn get(&self, model: &Model) -> &Exact {
        self.made.get_or_init(|| match self.source {
            Source::Face(face) => {
                let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
                model.surface(local).map_or(Exact::Unknown, |s| Exact::Surface(s.transformed(&face.placement)))
            }
            Source::Edge(edge) => model.edge_curve(edge).map_or(Exact::Unknown, |(c, first, last)| Exact::Curve(c, first, last)),
            Source::Point(p) => Exact::Point(p),
            Source::Unknown => Exact::Unknown,
        })
    }

    fn is_surface(&self) -> bool {
        matches!(self.source, Source::Face(_))
    }
}

impl Exact {
    /// The point of the geometry nearest `x`, with its parameters (starting from `hint` on a
    /// surface).
    fn nearest(&self, x: Point, hint: Option<[f64; 2]>) -> Option<(Point, [f64; 2])> {
        match self {
            Exact::Surface(s) => {
                let (u, v, _) = crate::geometry::extrema::nearest_on_surface(x, s, hint.map(|h| (h[0], h[1])));
                Some((s.point(u, v), [u, v]))
            }
            Exact::Curve(c, first, last) => {
                let found = crate::geometry::extrema::nearest_on_curve(x, c, *first, *last);
                Some((found.point, [found.parameter, 0.0]))
            }
            Exact::Point(p) => Some((*p, [0.0; 2])),
            Exact::Unknown => None,
        }
    }

    /// How many parameters the geometry has.
    fn dimension(&self) -> usize {
        match self {
            Exact::Surface(_) => 2,
            Exact::Curve(..) => 1,
            _ => 0,
        }
    }

    /// The point at `w` and the derivatives along each parameter.
    fn evaluate(&self, w: [f64; 2]) -> (Point, [Vector; 2]) {
        match self {
            Exact::Surface(s) => {
                let (p, du, dv) = s.d1(w[0], w[1]);
                (p, [du, dv])
            }
            Exact::Curve(c, ..) => {
                let (p, d) = c.d1(w[0]);
                (p, [d, Vector::new(0.0, 0.0, 0.0)])
            }
            Exact::Point(p) => (*p, [Vector::new(0.0, 0.0, 0.0); 2]),
            Exact::Unknown => unreachable!("no geometry"),
        }
    }

    /// `w` kept within an edge's range.
    fn clamp(&self, w: [f64; 2]) -> [f64; 2] {
        match self {
            Exact::Curve(_, first, last) => [w[0].clamp(first.min(*last), first.max(*last)), 0.0],
            _ => w,
        }
    }
}

impl Triangles {
    /// The shape's triangles (meshed first, if not yet finely enough) and their box tree,
    /// kept in the model with it.
    /// `print` is the shape's fingerprint (which its meshes do not change).
    fn kept(model: &mut Model, shape: Shape, options: MeshOptions, print: u64) -> Result<std::sync::Arc<KeptTriangles>> {
        meshing::mesh_with(model, shape, options)?;
        // (The meshes as they are after meshing: they are what is gathered.)
        let key = print ^ model.mesh_fingerprint(shape).rotate_left(13) ^ options.deflection.to_bits().rotate_left(29) ^ options.angle.to_bits().rotate_left(47);
        let model: &Model = model;
        Ok(model.derived(shape.node, key, || {
            let triangles = Triangles::gather(model, shape);
            let bvh = triangles.bvh();
            KeptTriangles(triangles, bvh)
        }))
    }

    /// The triangles of the shape's faces' meshes (and its free edges and vertices).
    fn gather(model: &Model, shape: Shape) -> Triangles {
        let mut out = Triangles { nodes: Vec::new(), triangles: Vec::new(), face: Vec::new(), exact: Vec::new(), deviation: Vec::new(), ranges: Vec::new(), edges: Vec::new(), piece_edges: Vec::new() };
        let mut known: HashMap<crate::topology::NodeId, Vec<(Transform, usize)>> = HashMap::new();
        for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
            let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
            let Some(mesh) = model.face_mesh(local) else { continue };
            let base = out.nodes.len() as u32;
            out.nodes.extend(mesh.nodes().iter().map(|p| p.transformed(&face.placement)));
            let k = out.exact.len();
            for t in mesh.triangles() {
                out.triangles.push(t.map(|i| i + base));
                out.face.push(k);
            }
            out.exact.push(Lazy::new(if model.surface(local).is_some() { Source::Face(face) } else { Source::Unknown }));
            out.deviation.push(mesh.deflection());
            out.ranges.push(out.triangles.len() - mesh.triangles().len()..out.triangles.len());
            let mut boundary = Vec::new();
            for edge in model.explore(face, NodeKind::Edge).collect::<Vec<_>>() {
                let seen = known.entry(edge.node).or_default();
                if let Some(&(_, i)) = seen.iter().find(|(t, _)| *t == edge.placement) {
                    boundary.push(i);
                    continue;
                }
                let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..edge };
                let (true, Some(polyline)) = (model.has_edge_curve(edge), model.edge_polyline(local)) else { continue };
                let points: Vec<Point> = polyline.nodes().iter().map(|p| p.transformed(&edge.placement)).collect();
                let deviation = sag(model, edge, &points).min(polyline.deflection());
                out.edges.push((Lazy::new(Source::Edge(edge)), points, deviation));
                seen.push((edge.placement, out.edges.len() - 1));
                boundary.push(out.edges.len() - 1);
            }
            out.piece_edges.push(boundary);
        }
        // Edges and vertices outside the faces, as triangles with coincident corners.
        let in_faces: HashSet<_> = model.explore(shape, NodeKind::Face).flat_map(|f| model.explore(f, NodeKind::Edge)).map(|e| e.node).collect();
        let mut on_edges: HashSet<_> = HashSet::new();
        for edge in model.explore(shape, NodeKind::Edge).collect::<Vec<_>>() {
            on_edges.extend(model.explore(edge, NodeKind::Vertex).map(|v| v.node));
            if in_faces.contains(&edge.node) {
                continue;
            }
            let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..edge };
            let Some(polyline) = model.edge_polyline(local) else { continue };
            let base = out.nodes.len() as u32;
            out.nodes.extend(polyline.nodes().iter().map(|p| p.transformed(&edge.placement)));
            let k = out.exact.len();
            let start = out.triangles.len();
            for i in 1..polyline.nodes().len() as u32 {
                out.triangles.push([base + i - 1, base + i, base + i]);
                out.face.push(k);
            }
            out.exact.push(Lazy::new(if model.has_edge_curve(edge) { Source::Edge(edge) } else { Source::Unknown }));
            out.deviation.push(polyline.deflection());
            out.ranges.push(start..out.triangles.len());
            out.piece_edges.push(Vec::new());
        }
        for vertex in model.explore(shape, NodeKind::Vertex).collect::<Vec<_>>() {
            if on_edges.insert(vertex.node) {
                let p = model.point(vertex);
                out.nodes.push(p);
                let n = out.nodes.len() as u32 - 1;
                out.ranges.push(out.triangles.len()..out.triangles.len() + 1);
                out.triangles.push([n, n, n]);
                out.face.push(out.exact.len());
                out.exact.push(Lazy::new(Source::Point(p)));
                out.deviation.push(0.0);
                out.piece_edges.push(Vec::new());
            }
        }
        out
    }

    fn triangle_box(&self, i: usize) -> BoundingBox {
        let [a, b, c] = self.corners(i);
        BoundingBox::new(
            Point::new(a.x.min(b.x).min(c.x), a.y.min(b.y).min(c.y), a.z.min(b.z).min(c.z)),
            Point::new(a.x.max(b.x).max(c.x), a.y.max(b.y).max(c.y), a.z.max(b.z).max(c.z)),
        )
    }

    /// Whether `p` lies within the mesh's deviation of piece `k`'s mesh (so on the piece, not
    /// out on its surface's extension).
    fn near_piece(&self, k: usize, p: Point, scale: f64) -> bool {
        let reach = 2.0 * self.deviation[k] + 1e-9 * scale;
        self.ranges[k].clone().any(|i| point_triangle(p, self.corners(i)).distance(p) <= reach)
    }

    /// Whether `p` lies on piece `k` for sure: on a face's mesh and inside its boundary, as
    /// the boundary's polylines show where it is clear of them and as the face's exact bounds
    /// do by them (a point on the boundary, or beyond it on the surface, is left to the
    /// boundary's edges; as OCCT's BRepExtrema_ExtFF keeps a face's extremum where its
    /// BRepClass_FaceClassifier finds it in the face). Edges and points: always.
    fn surely_on_piece(&self, model: &Model, k: usize, p: Point, scale: f64) -> bool {
        if !self.exact[k].is_surface() {
            return true;
        }
        if !self.near_piece(k, p, scale) {
            return false;
        }
        let clear = self.piece_edges[k].iter().all(|&e| {
            let (_, points, deviation) = &self.edges[e];
            polyline_distance(points, p) > 2.0 * (deviation + self.deviation[k]) + 1e-9 * scale
        });
        clear || self.exact[k].inside_bounds(model, p)
    }

    /// Piece `k`'s boundary edges coming within `within` of `p`.
    fn edges_near(&self, k: usize, p: Point, within: f64) -> Vec<usize> {
        self.piece_edges[k].iter().copied().filter(|&e| polyline_distance(&self.edges[e].1, p) <= within).collect()
    }

    fn corners(&self, i: usize) -> [Point; 3] {
        self.triangles[i].map(|k| self.nodes[k as usize])
    }

    fn bvh(&self) -> Option<Bvh> {
        if self.triangles.is_empty() {
            return None;
        }
        let boxes: Vec<BoundingBox> = (0..self.triangles.len()).map(|i| self.triangle_box(i)).collect();
        Some(Bvh::new(&boxes))
    }
}

/// The distance from `p` to a polyline.
/// How far an edge's curve strays from its polyline `points` (placed as the edge), as the
/// middles of its segments are from the curve, with a margin: the deflection a polyline
/// records is the one asked for, which a line's (exact) or a gentle curve's is well within,
/// and a point inside a face by its boundary is left to the boundary's edges only within
/// what they stray.
fn sag(model: &Model, edge: Shape, points: &[Point]) -> f64 {
    let Some((curve, first, last)) = model.edge_curve(edge) else { return f64::INFINITY };
    if matches!(curve, Curve::Line(_)) {
        return 0.0;
    }
    let most = points
        .windows(2)
        .map(|w| {
            let middle = w[0] + (w[1] - w[0]) * 0.5;
            crate::geometry::extrema::nearest_on_curve(middle, &curve, first, last).distance
        })
        .fold(0.0, f64::max);
    1.5 * most
}

fn polyline_distance(points: &[Point], p: Point) -> f64 {
    if points.len() == 1 {
        return points[0].distance(p);
    }
    points
        .windows(2)
        .map(|w| {
            let (a, b) = (w[0], w[1]);
            let ab = b - a;
            let t = if ab.dot(ab) > 0.0 { ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0) } else { 0.0 };
            (a + ab * t).distance(p)
        })
        .fold(f64::INFINITY, f64::min)
}

/// The closest point of triangle `t` to `p` (Ericson, Real-Time Collision Detection 5.1.5).
fn point_triangle(p: Point, [a, b, c]: [Point; 3]) -> Point {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    a + ab * (vb * denom) + ac * (vc * denom)
}

/// The closest points of segments `p1 q1` and `p2 q2`.
fn segment_segment(p1: Point, q1: Point, p2: Point, q2: Point) -> (Point, Point) {
    let (d1, d2, r) = (q1 - p1, q2 - p2, p1 - p2);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (s, t) = if a <= 1e-300 && e <= 1e-300 {
        (0.0, 0.0)
    } else if a <= 1e-300 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(r);
        if e <= 1e-300 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            let mut s = if denom > 0.0 { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
            let mut t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
            (s, t)
        }
    };
    (p1 + d1 * s, p2 + d2 * t)
}

/// Where segment `p q` crosses triangle `t`, if it does.
fn segment_triangle(p: Point, q: Point, [a, b, c]: [Point; 3]) -> Option<Point> {
    let n = (b - a).cross(c - a);
    let (dp, dq) = ((p - a).dot(n), (q - a).dot(n));
    if dp * dq > 0.0 || (dp == 0.0 && dq == 0.0) {
        return None;
    }
    let x = p + (q - p) * (dp / (dp - dq));
    let inside = |u: Point, v: Point| (v - u).cross(x - u).dot(n) >= 0.0;
    (inside(a, b) && inside(b, c) && inside(c, a)).then_some(x)
}

/// The distance between two triangles and the closest points.
fn triangle_triangle(s: [Point; 3], t: [Point; 3]) -> (f64, Point, Point) {
    for k in 0..3 {
        if let Some(x) = segment_triangle(s[k], s[(k + 1) % 3], t) {
            return (0.0, x, x);
        }
        if let Some(x) = segment_triangle(t[k], t[(k + 1) % 3], s) {
            return (0.0, x, x);
        }
    }
    let mut best = (f64::INFINITY, s[0], t[0]);
    let mut consider = |p: Point, q: Point| {
        let d = p.distance(q);
        if d < best.0 {
            best = (d, p, q);
        }
    };
    for k in 0..3 {
        consider(s[k], point_triangle(s[k], t));
        consider(point_triangle(t[k], s), t[k]);
        for l in 0..3 {
            let (p, q) = segment_segment(s[k], s[(k + 1) % 3], t[l], t[(l + 1) % 3]);
            consider(p, q);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Direction, Frame};
    use crate::primitives::{Cuboid, Cylinder, Sphere};

    /// Prisms on an ellipse (as a periodic rational B-spline) and on the same turned by
    /// 1e-4 about its centre, the second raised: they overlap (the ellipses' centre lies in
    /// both over the height they share), at distance zero either way round. Their sides
    /// cross at angles of about 1e-4, the meshes crossing over long stretches about each
    /// crossing; the nearest points refined on the exact sides near where the meshes cross
    /// were found apart.
    #[test]
    fn prisms_crossing_at_shallow_angles_are_at_no_distance() {
        use crate::geometry::curve::{Curve, Ellipse};
        let prism = |model: &mut Model, turn: f64, z: f64| {
            let x = Direction::new(turn.cos(), turn.sin(), 0.0);
            let ellipse = Curve::Ellipse(Ellipse::new(Frame::new(Point::new(0.0, 0.0, z), Direction::Z, x), 5.0, 3.0));
            let mut spline = crate::geometry::convert::curve_to_bspline(&ellipse).unwrap();
            spline.set_origin(1);
            let edge = crate::make::edge_from_curve(model, Curve::BSpline(spline)).unwrap();
            let wire = crate::make::wire(model, &[edge]).unwrap();
            let face = crate::make::planar_face(model, wire).unwrap();
            crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, 4.0)).unwrap().shape
        };
        let mut model = Model::new();
        let a = prism(&mut model, 0.0, 0.0);
        let b = prism(&mut model, 1e-4, 1.5);
        for (x, y) in [(a, b), (b, a)] {
            let d = distance(&mut model, x, y).unwrap();
            assert_eq!(d.value, 0.0, "{d:?}");
        }
    }

    /// The nearest faces of shapes on elementary surfaces, found without meshes where that
    /// can be told, are those the meshes find: from points all about them, by their edges
    /// and corners (within the classifiers' tolerance and the meshes' deviations), and
    /// across a cone's apex (where the surface's nearest point is on its other nappe).
    #[test]
    fn nearest_faces_without_meshes_are_those_the_meshes_find() {
        use crate::primitives::{Cone, Torus};
        let tilted = Frame::new(Point::new(1.5, -2.0, 0.7), Direction::new(0.3, 0.4, 0.86), Direction::X);
        let mut model = Model::new();
        let mut shapes = vec![
            Cuboid::new(5.0, 7.0, 3.0).build(&mut model).unwrap(),
            Cylinder::new(2.5, 6.0).with_frame(tilted).build(&mut model).unwrap(),
            Cone::new(3.0, 0.0, 5.0).build(&mut model).unwrap(),
            Sphere::new(2.75).with_frame(tilted).build(&mut model).unwrap(),
            Torus::new(6.0, 2.0).with_frame(tilted).build(&mut model).unwrap(),
        ];
        let block = Cuboid::from_corners(Point::new(0.5, -10.0, -10.0), Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
        let ball = Sphere::new(3.0).build(&mut model).unwrap();
        shapes.push(crate::boolean::difference(&mut model, ball, block).unwrap());
        let mut state = 0x9e3779b97f4a7c15u64;
        let mut random = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let (mut exact, mut total) = (0, 0);
        for shape in shapes {
            let b = crate::properties::bounding_box(&model, shape);
            let mut points: Vec<Point> = (0..60)
                .map(|_| {
                    let t = [random() * 1.4 - 0.2, random() * 1.4 - 0.2, random() * 1.4 - 0.2];
                    Point::new(b.min.x + (b.max.x - b.min.x) * t[0], b.min.y + (b.max.y - b.min.y) * t[1], b.min.z + (b.max.z - b.min.z) * t[2])
                })
                .collect();
            for edge in model.unique_subshapes(shape, NodeKind::Edge) {
                let Some((c, first, last)) = model.edge_curve(edge) else { continue };
                for k in 0..=4 {
                    let x = c.point(first + (last - first) * k as f64 / 4.0);
                    for h in [0.0, 1e-7, 1e-4, 0.1] {
                        points.push(x + Vector::new(random() - 0.5, random() - 0.5, random() - 0.5) * h);
                    }
                }
            }
            points.push(Point::new(0.9232894240000435, 3.0863540094379074, 5.211782089823887));
            for p in points {
                let meshed = nearest_faces_by(&mut model, shape, p, crate::tolerance::LINEAR, false).unwrap();
                assert_eq!(nearest_faces(&mut model, shape, p, crate::tolerance::LINEAR).unwrap(), meshed, "{p:?}");
                let print = model.shape_fingerprint(shape);
                let scale = model.derived(shape.node, print, || KeptSize(0.0)).0.max(1e-9);
                total += 1;
                exact += nearest_elementary_faces(&mut model, shape, print, p, crate::tolerance::LINEAR, scale, 1e-3 * scale).unwrap().is_some() as usize;
            }
        }
        // (Most without meshes.)
        assert!(2 * exact > total, "{exact} of {total}");
    }

    /// Points inside a face 10.4 x 40.1, by its edges (as a hole's corner just inside the
    /// outline it is cut from): on the face, at no distance; by its corner outside, as far
    /// as the corner.
    #[test]
    fn points_inside_a_face_by_its_edges() {
        let mut model = Model::new();
        let corners = [(-5.2, -20.05), (5.2, -20.05), (5.2, 20.05), (-5.2, 20.05)];
        let edges: Vec<Shape> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                crate::make::edge_between(&mut model, Point::new(a.0, a.1, 0.0), Point::new(b.0, b.1, 0.0)).unwrap()
            })
            .collect();
        let wire = crate::make::wire(&mut model, &edges).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        for (x, y, expected) in [(4.0, 20.0, 0.0), (0.0, 20.04, 0.0), (5.19, -20.0, 0.0), (5.3, 20.15, 0.02f64.sqrt())] {
            let v = crate::make::vertex(&mut model, Point::new(x, y, 0.0));
            let d = distance(&mut model, face, v).unwrap().value;
            assert!((d - expected).abs() <= 1e-9, "({x}, {y}): {d} != {expected}");
        }
    }

    #[test]
    fn distances_between_simple_solids() {
        let mut model = Model::new();
        let a = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let at = |x: f64, y: f64, z: f64| Frame::new(Point::new(x, y, z), Direction::Z, Direction::X);
        let b = Cuboid::new(1.0, 1.0, 1.0).with_frame(at(3.0, 0.0, 0.0)).build(&mut model).unwrap();
        assert!((distance(&mut model, a, b).unwrap().value - 2.0).abs() < 1e-12);
        // Touching along a face.
        let c = Cuboid::new(1.0, 1.0, 1.0).with_frame(at(1.0, 0.5, 0.0)).build(&mut model).unwrap();
        assert_eq!(distance(&mut model, a, c).unwrap().value, 0.0);
        // A sphere of radius 1 centred 3 above the box's top: 2 away.
        let s = Sphere::new(1.0).with_frame(at(0.5, 0.5, 4.0)).build(&mut model).unwrap();
        assert!((distance(&mut model, a, s).unwrap().value - 2.0).abs() < 1e-6);
        // A cylinder lying on the box's top (touching along a line).
        let lying = Frame::new(Point::new(-1.0, 0.5, 1.5), Direction::X, Direction::Y);
        let cyl = Cylinder::new(0.5, 3.0).with_frame(lying).build(&mut model).unwrap();
        let dc = distance(&mut model, a, cyl).unwrap().value;
        assert!(dc < 1e-7, "{dc}");
        // A torus resting on the box's top (touching along a circle).
        let ring = crate::primitives::Torus::new(0.3, 0.1).with_frame(at(0.5, 0.5, 1.1)).build(&mut model).unwrap();
        assert!(distance(&mut model, a, ring).unwrap().value < 1e-9);
        // A small box inside the big one.
        let inner = Cuboid::new(0.2, 0.2, 0.2).with_frame(at(0.4, 0.4, 0.4)).build(&mut model).unwrap();
        assert_eq!(distance(&mut model, a, inner).unwrap().value, 0.0);
    }

    #[test]
    fn a_rod_in_a_hole_is_not_touching_though_the_meshes_cross() {
        let mut model = Model::new();
        let at = |z: f64| Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X);
        let block = Cuboid::new(4.0, 4.0, 2.0).with_frame(Frame::new(Point::new(-2.0, -2.0, 0.0), Direction::Z, Direction::X)).build(&mut model).unwrap();
        let cutter = Cylinder::new(1.0, 4.0).with_frame(at(-1.0)).build(&mut model).unwrap();
        let holed = crate::boolean::difference(&mut model, block, cutter).unwrap();
        // 0.002 of clearance, less than the meshes' deviation: the hole's chords cut into the rod.
        let rod = Cylinder::new(0.998, 4.0).with_frame(at(-1.0)).build(&mut model).unwrap();
        let d = distance(&mut model, holed, rod).unwrap().value;
        assert!((d - 0.002).abs() < 1e-9, "{d}");
    }

    /// A shaft pressed into a tube's bore by less than the meshes stray: its side crosses the
    /// tube's ends by their inner edges, where the meshes alone cannot tell it on them.
    #[test]
    fn a_shaft_pressed_into_a_bore_by_a_hair_touches_it() {
        let mut model = Model::new();
        let at = |z: f64| Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X);
        for press in [1.005, 1.01] {
            let outer = Cylinder::new(2.0, 4.0).build(&mut model).unwrap();
            let bore = Cylinder::new(1.0, 4.0).build(&mut model).unwrap();
            let tube = crate::boolean::difference(&mut model, outer, bore).unwrap();
            let shaft = Cylinder::new(press, 6.0).with_frame(at(-1.0)).build(&mut model).unwrap();
            let d = distance(&mut model, tube, shaft).unwrap().value;
            assert!(d < 1e-9, "{press}: {d}");
        }
    }

    #[test]
    fn distances_from_points_and_edges() {
        let mut model = Model::new();
        let cylinder = Cylinder::new(4.0, 6.0).build(&mut model).unwrap();
        let side = model.explore(cylinder, NodeKind::Face).find(|&f| matches!(model.surface(f), Some(Surface::Cylinder(_)))).unwrap();
        // Inside the side face, off its axis (the meshes' chords would say less).
        let inside = crate::make::vertex(&mut model, Point::new(3.0, 1.0, 3.0));
        let d = distance(&mut model, inside, side).unwrap().value;
        assert!((d - (4.0 - 10f64.sqrt())).abs() < 1e-9, "{d}");
        let outside = crate::make::vertex(&mut model, Point::new(0.0, 10.0, 3.0));
        assert!((distance(&mut model, side, outside).unwrap().value - 6.0).abs() < 1e-9);
        // A point and a straight edge.
        let (p, q) = (crate::make::vertex(&mut model, Point::new(0.0, 0.0, 10.0)), crate::make::vertex(&mut model, Point::new(1.0, 0.0, 10.0)));
        let edge = crate::make::edge(&mut model, p, q).unwrap();
        assert!((distance(&mut model, edge, outside).unwrap().value - (100.0f64 + 49.0).sqrt()).abs() < 1e-9);
    }

    /// A point and faces whose nearest points are not isolated (on a cylinder's axis, at a
    /// sphere's centre), or lie on the face's rim; and a point and a B-spline edge: found
    /// exactly, without meshes, as the meshes find them.
    #[test]
    fn distances_found_exactly_where_extrema_are_not_isolated_or_on_a_rim() {
        let mut model = Model::new();
        let cylinder = Cylinder::new(4.0, 6.0).build(&mut model).unwrap();
        let side = model.explore(cylinder, NodeKind::Face).find(|&f| matches!(model.surface(f), Some(Surface::Cylinder(_)))).unwrap();
        // On the axis, within the side's height: its radius; above it, to its rim.
        for (z, want) in [(3.0, 4.0), (8.0, (16.0f64 + 4.0).sqrt())] {
            let v = crate::make::vertex(&mut model, Point::new(0.0, 0.0, z));
            let exact = vertex_to_faces(&model, v, side).expect("found exactly").0;
            assert!((exact - want).abs() < 1e-9, "{z}: {exact}");
            let meshed = distance_within(&mut model, side, v, f64::INFINITY).unwrap().value;
            assert!((exact - meshed).abs() < 1e-6, "{z}: {exact} {meshed}");
        }
        // A hemisphere from its centre (every point of it as near), and from a point whose
        // foot lies on its rim.
        let ball = Sphere::new(2.0).build(&mut model).unwrap();
        let below = Cuboid::from_corners(Point::new(-5.0, -5.0, -5.0), Point::new(5.0, 5.0, 0.0)).build(&mut model).unwrap();
        let half = crate::boolean::difference(&mut model, ball, below).unwrap();
        let dome = model.explore(half, NodeKind::Face).find(|&f| matches!(model.surface(f), Some(Surface::Sphere(_)))).unwrap();
        for (p, want) in [(Point::new(0.0, 0.0, 0.0), 2.0), (Point::new(3.0, 0.0, 0.0), 1.0), (Point::new(0.0, 4.0, -3.0), (4.0f64 - 2.0).hypot(3.0))] {
            let v = crate::make::vertex(&mut model, p);
            let exact = vertex_to_faces(&model, v, dome).expect("found exactly").0;
            assert!((exact - want).abs() < 1e-9, "{p:?}: {exact}");
        }
        // A B-spline edge (a circle as one).
        let arc = crate::geometry::curve::Curve::Circle(crate::geometry::curve::Circle::new(Frame::from_z(Point::ORIGIN, Direction::Z), 5.0));
        let spline = crate::geometry::convert::curve_to_bspline(&arc).unwrap();
        let edge = crate::make::edge_from_curve(&mut model, crate::geometry::curve::Curve::BSpline(spline)).unwrap();
        for (p, want) in [(Point::new(1.0, 1.0, 1.0), (5.0 - 2f64.sqrt()).hypot(1.0)), (Point::new(9.0, 3.0, 0.0), 90f64.sqrt() - 5.0)] {
            let v = crate::make::vertex(&mut model, p);
            let exact = vertex_to_edges(&model, v, edge).expect("found exactly").0;
            assert!((exact - want).abs() < 1e-9, "{p:?}: {exact}");
        }
    }
}
