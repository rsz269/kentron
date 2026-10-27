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
            if (check_a && !ta.surely_on_piece(ka, r.1, scale)) || (check_b && !tb.surely_on_piece(kb, r.2, scale)) {
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
        Lazy { source, made: std::sync::OnceLock::new() }
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

    /// Whether `p` lies on piece `k` for sure: on a face's mesh and clear of its boundary
    /// (a point by the boundary is left to the boundary's edges). Edges and points: always.
    fn surely_on_piece(&self, k: usize, p: Point, scale: f64) -> bool {
        if !self.exact[k].is_surface() {
            return true;
        }
        self.near_piece(k, p, scale)
            && self.piece_edges[k].iter().all(|&e| {
                let (_, points, deviation) = &self.edges[e];
                polyline_distance(points, p) > 2.0 * (deviation + self.deviation[k]) + 1e-9 * scale
            })
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
}

