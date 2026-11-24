//! Polyhedra grown with rounded joins, built directly.
//!
//! A solid grown by `d > 0` with its joins rounded (by arcs) is the set of points within `d`
//! of it: its Minkowski sum with a ball of radius `d`. Its boundary is made of the points at
//! distance `d` whose nearest point on the solid lies on one face, edge or vertex of it. For a
//! solid of planar faces those parts are, where nothing far apart on the solid comes within
//! `2d` of itself across the outside:
//!
//! * each face moved `d` along its outward normal;
//! * along each convex edge, the part of the cylinder of radius `d` about it between its
//!   faces' normals (the normals of the points round the edge turning from one face's to the
//!   other's);
//! * at each vertex whose edges are all convex, the part of the sphere of radius `d` about it
//!   bounded by arcs of great circles from each face's normal to the next one's round the
//!   vertex (the directions whose nearest point is the vertex: the cone the normals span).
//!
//! Along a concave edge the two faces' offsets meet in a line (no point has its nearest point
//! on the edge alone). Where a concave edge ends at a vertex with two convex edges square to
//! it (an inner corner of a prism's outline, at its cap), the cylinders of those two edges
//! meet in a mitre: the ellipse where the plane through the vertex equally far from both
//! edges' lines cuts them, from the cap's corner moved to the end of the concave edge's line.
//!
//! The parts meet edge to edge, so the result is built from them exactly (planes, cylinders and
//! spheres sharing edges by construction) without moving the faces and rounding where they
//! part by fillets. It is built so for
//!
//! * a convex polyhedron (every edge convex): it is its own Minkowski sum's boundary's
//!   description for every `d`;
//! * a right prism over a simple polygon (two caps, the other faces rectangles square to
//!   them): its Minkowski sum's section at each height is the outline grown by `r` in its
//!   plane, `r` from 0 (over the caps) to `d` (along the sides), so the parts above are right
//!   when the outline grown by every `r` up to `d` has the parts its corners and sides give
//!   ([`outline_grows_regularly`]: no side shrinking away between inner corners, nothing
//!   coming within `r` of a side it is not next to).
//!
//! Other solids (curved faces, faces with holes, other concave edges) go the general way
//! ([`super::shape`]).

use std::collections::{BTreeMap, HashMap};

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Ellipse, Line};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{Cylinder, Plane, Sphere, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Direction2, Frame, Point, Point2, Transformable, Vector, Vector2};
use crate::topology::{EdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

/// How far apart (as the sine of the angle between them) directions taken as square to one
/// another or parallel may be.
const SQUARE: f64 = 1e-9;

/// A face of the polyhedron: its outward normal and its corners (by index), anticlockwise
/// about the normal.
struct Face {
    shape: Shape,
    plane: Plane,
    normal: Vector,
    corners: Vec<usize>,
    /// The edges from each corner to the next.
    edges: Vec<NodeId>,
}

/// An edge of the polyhedron: from `start` to `end`, as `faces[0]` runs (`faces[1]` the
/// other way), its direction, length, whether it is convex, and the angle its faces'
/// normals turn through.
struct Edge {
    start: usize,
    end: usize,
    faces: [usize; 2],
    original: NodeId,
    axis: Vector,
    length: f64,
    convex: bool,
    angle: f64,
}

/// A vertex where a concave edge ends between two convex edges square to it: the face the
/// convex edges bound (square to the concave edge), and where the offsets of the concave
/// edge's faces meet in the plane through the vertex square to it.
struct Mitre {
    cap: usize,
    convex: [usize; 2],
    meet: Point,
}

/// `shape` grown by `distance` (> 0) with rounded joins, built directly where `shape` is a
/// convex polyhedron or a right prism over a simple polygon growing regularly (see the
/// module's documentation): a solid of one closed shell whose faces are planar, each bounded
/// by one wire of straight edges, no two faces meeting tangentially. `None` otherwise.
pub(super) fn polyhedron_grown(model: &mut Model, shape: Shape, distance: f64) -> Result<Option<Shape>> {
    if !(distance > 0.0) {
        return Ok(None);
    }
    let Some((points, faces)) = polyhedron(model, shape) else { return Ok(None) };
    // Each directed edge (from corner to corner) to the face it bounds that way round and
    // the original edge.
    let mut half: HashMap<(usize, usize), (usize, NodeId)> = HashMap::new();
    for (f, face) in faces.iter().enumerate() {
        let n = face.corners.len();
        for k in 0..n {
            if half.insert((face.corners[k], face.corners[(k + 1) % n]), (f, face.edges[k])).is_some() {
                return Ok(None);
            }
        }
    }
    // Each edge used once each way (a closed manifold shell), convex or concave.
    let mut edges: Vec<Edge> = Vec::new();
    for (&(a, b), &(f, original)) in &half {
        let Some(&(g, other)) = half.get(&(b, a)) else { return Ok(None) };
        if other != original {
            return Ok(None);
        }
        if a > b {
            continue;
        }
        // (`f` runs from a to b, anticlockwise about its normal.)
        let along = points[b] - points[a];
        let length = along.length();
        if !(length > crate::tolerance::LINEAR) {
            return Ok(None);
        }
        let axis = along * (1.0 / length);
        let (nf, ng) = (faces[f].normal, faces[g].normal);
        let turn = nf.cross(ng);
        let (sin, cos) = (turn.dot(axis), nf.dot(ng));
        // (Convex where the normals turn about the edge as it runs along `f`; not tangent.)
        if !(sin.abs() > 1e3 * crate::tolerance::ANGULAR.sqrt()) || (turn - axis * sin).length() > SQUARE {
            return Ok(None);
        }
        edges.push(Edge { start: a, end: b, faces: [f, g], original, axis, length, convex: sin > 0.0, angle: sin.abs().atan2(cos) });
    }
    // Euler's characteristic of a sphere.
    if points.len() + faces.len() != edges.len() + 2 {
        return Ok(None);
    }
    // (In a fixed order, whatever the map's.)
    edges.sort_by_key(|e| (e.start, e.end));
    let mitres = if edges.iter().all(|e| e.convex) {
        // Each face's corners convex (so it is for a valid convex polyhedron; checked, cheaply).
        for face in &faces {
            let n = face.corners.len();
            for k in 0..n {
                let (p, q, r) = (points[face.corners[(k + n - 1) % n]], points[face.corners[k]], points[face.corners[(k + 1) % n]]);
                if !((q - p).cross(r - q).dot(face.normal) > 0.0) {
                    return Ok(None);
                }
            }
        }
        BTreeMap::new()
    } else {
        if !prism_grows_regularly(&points, &faces, &edges, distance) {
            return Ok(None);
        }
        let Some(mitres) = mitres(&points, &faces, &edges, distance) else { return Ok(None) };
        mitres
    };
    build(model, &points, &faces, &edges, &mitres, distance).map(Some)
}

/// The corners and the faces of `shape` where it is one solid of planar faces, each with
/// one wire of straight edges.
fn polyhedron(model: &Model, shape: Shape) -> Option<(Vec<Point>, Vec<Face>)> {
    if shape.kind() != NodeKind::Solid || model.unique_subshapes(shape, NodeKind::Shell).len() != 1 {
        return None;
    }
    let mut points: Vec<Point> = Vec::new();
    let mut nodes: Vec<NodeId> = Vec::new();
    let mut index = |model: &Model, v: Shape| -> usize {
        match nodes.iter().position(|&n| n == v.node) {
            Some(i) => i,
            None => {
                nodes.push(v.node);
                points.push(model.point(v));
                nodes.len() - 1
            }
        }
    };
    let mut faces = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(Surface::Plane(plane)) = model.surface(face).map(|s| s.transformed(&face.placement)) else { return None };
        if model.children(face).count() != 1 {
            return None;
        }
        let normal = Vector::from(plane.normal()) * if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
        // The face's edges as it uses them, chained into one loop.
        let mut runs: Vec<(usize, usize, NodeId)> = Vec::new();
        for edge in model.explore(face, NodeKind::Edge) {
            let Some((Curve::Line(_), _, _)) = model.edge_curve(edge) else { return None };
            let (a, b) = model.edge_vertices(edge);
            let (a, b) = (index(model, a), index(model, b));
            if a == b {
                return None;
            }
            runs.push((a, b, edge.node));
        }
        let first = runs.first()?;
        let (mut corners, mut edges) = (vec![first.0], vec![first.2]);
        let mut at = first.1;
        while at != corners[0] {
            if corners.len() > runs.len() || corners.contains(&at) {
                return None;
            }
            let next = runs.iter().find(|r| r.0 == at)?;
            corners.push(at);
            edges.push(next.2);
            at = next.1;
        }
        if corners.len() != runs.len() || corners.len() < 3 {
            return None;
        }
        faces.push(Face { shape: face, plane, normal, corners, edges });
    }
    // Every corner on its faces' planes.
    for face in &faces {
        if face.corners.iter().any(|&c| face.plane.distance(points[c]).abs() > crate::tolerance::LINEAR.max(model.tolerance(face.shape))) {
            return None;
        }
    }
    Some((points, faces))
}

/// Whether the polyhedron is a right prism (two faces facing apart, its caps, holding every
/// corner; every other face square to them; every edge not in a cap square to the caps)
/// whose outline grows regularly by `d`.
fn prism_grows_regularly(points: &[Point], faces: &[Face], edges: &[Edge], d: f64) -> bool {
    let caps = |top: usize, bottom: usize| -> bool {
        let up = faces[top].normal;
        if (up + faces[bottom].normal).length() > SQUARE {
            return false;
        }
        if faces.iter().enumerate().any(|(f, face)| f != top && f != bottom && face.normal.dot(up).abs() > SQUARE) {
            return false;
        }
        let (n, m) = (faces[top].corners.len(), faces[bottom].corners.len());
        if n != m || 2 * n != points.len() {
            return false;
        }
        let in_cap = |f: usize, c: usize| faces[f].corners.contains(&c);
        edges.iter().all(|e| {
            let both = |f: usize| in_cap(f, e.start) && in_cap(f, e.end);
            both(top) || both(bottom) || e.axis.cross(up).length() <= SQUARE
        })
    };
    for top in 0..faces.len() {
        for bottom in 0..faces.len() {
            if top == bottom || !caps(top, bottom) {
                continue;
            }
            // The outline: the top's corners in its plane, anticlockwise.
            let up = faces[top].normal;
            let x = Vector::from(faces[top].plane.frame().x());
            let x = (x - up * x.dot(up)).normalized();
            let y = up.cross(x);
            let o = faces[top].plane.origin();
            let outline: Vec<Point2> = faces[top].corners.iter().map(|&c| Point2::new((points[c] - o).dot(x), (points[c] - o).dot(y))).collect();
            return outline_grows_regularly(&outline, d);
        }
    }
    false
}

/// Whether the simple polygon `q` (anticlockwise) grown by every `r` in (0, `d`] is bounded by
/// the parts its sides and corners give: each side moved `r` out, cut back at inner (reflex)
/// corners to where it meets the next side's (along the corner's bisector), and an arc of
/// radius `r` about each outer corner. So it is where
///
/// * every side moved by `d` keeps a length, cut back at its inner corners;
/// * every point of those parts at `d` (each side moved, each arc), and of each inner
///   corner's bisector at every `r` up to `d`, is at least `r` from every side not at that
///   side's or corner's ends (sides at its ends it is `r` from by construction).
///
/// (A point of a side's part at `r` is no nearer another side than the same point moved on
/// to `d` along its normal, by `d - r` at most; one cut back by `d` is on the corner's
/// bisector at the `r` it is cut back by.)
fn outline_grows_regularly(q: &[Point2], d: f64) -> bool {
    let n = q.len();
    let margin = crate::tolerance::LINEAR;
    let side = |i: usize| (q[i % n], q[(i + 1) % n]);
    let dir = |i: usize| {
        let (a, b) = side(i);
        (b - a).normalized()
    };
    // Outward normals (right of the sides, the outline anticlockwise).
    let normal = |i: usize| {
        let t = dir(i);
        Vector2::new(t.y, -t.x)
    };
    let prev = |i: usize| (i + n - 1) % n;
    // Corner i (at q[i], between sides i - 1 and i): inner, and how far it cuts its sides back
    // per unit grown.
    let mut inner = vec![false; n];
    let mut cut = vec![0.0; n];
    for i in 0..n {
        let (s, c) = (dir(prev(i)).cross(dir(i)), dir(prev(i)).dot(dir(i)));
        if s.abs() <= SQUARE {
            return false;
        }
        inner[i] = s < 0.0;
        if inner[i] {
            // tan of half the angle the normals turn through.
            cut[i] = s.abs() / (1.0 + c);
        }
    }
    let sides: Vec<(Point2, Point2)> = (0..n)
        .map(|i| {
            let (a, b) = side(i);
            let (t, nu) = (dir(i), normal(i));
            (a + nu * d + t * (d * cut[i]), b + nu * d + t * (-d * cut[(i + 1) % n]))
        })
        .collect();
    for i in 0..n {
        let (a, b) = side(i);
        if !((b - a).length() - d * (cut[i] + cut[(i + 1) % n]) > margin) {
            return false;
        }
    }
    // Sides j at corner i's ends: i - 1 and i.
    let at_corner = |i: usize, j: usize| j == prev(i) || j == i;
    let next_to = |i: usize, j: usize| j == i || j == (i + 1) % n || j == prev(i);
    for i in 0..n {
        for j in 0..n {
            let (a, b) = side(j);
            if !next_to(i, j) && segment_segment(sides[i].0, sides[i].1, a, b) < d + margin {
                return false;
            }
            if at_corner(i, j) {
                continue;
            }
            if inner[i] {
                let (n0, n1) = (normal(prev(i)), normal(i));
                let w = (n0 + n1) * (1.0 / (1.0 + n0.dot(n1)));
                if bisector_gap(q[i], w, d, a, b) < margin {
                    return false;
                }
            } else if arc_segment(q[i], d, normal(prev(i)), normal(i), a, b) < d + margin {
                return false;
            }
        }
    }
    true
}

/// The distance from `p` to the segment from `a` to `b`.
fn point_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    p.distance(a + ab * t)
}

/// The distance between the segments `a`-`b` and `c`-`d`.
fn segment_segment(a: Point2, b: Point2, c: Point2, d: Point2) -> f64 {
    let side = |p: Point2, q: Point2, r: Point2| (q - p).cross(r - p);
    let (s1, s2, s3, s4) = (side(a, b, c), side(a, b, d), side(c, d, a), side(c, d, b));
    if s1 * s2 < 0.0 && s3 * s4 < 0.0 {
        return 0.0;
    }
    point_segment(a, c, d).min(point_segment(b, c, d)).min(point_segment(c, a, b)).min(point_segment(d, a, b))
}

/// The distance from the arc of radius `r` about `c`, its directions turning anticlockwise
/// from `n0` to `n1` (less than half a turn), to the segment `a`-`b`: least where the arc
/// ends, faces an end of the segment, or faces square to the segment's line; nothing where
/// they cross.
fn arc_segment(c: Point2, r: f64, n0: Vector2, n1: Vector2, a: Point2, b: Point2) -> f64 {
    let within = |v: Vector2| n0.cross(v) >= 0.0 && v.cross(n1) >= 0.0;
    let ab = b - a;
    let t = ab.normalized();
    // Where the circle crosses the segment's line, on the segment and the arc.
    let foot = a + t * (c - a).dot(t);
    let h = foot.distance(c);
    if h <= r {
        let half = (r * r - h * h).sqrt();
        for p in [foot + t * half, foot + t * -half] {
            let s = (p - a).dot(t);
            if s >= 0.0 && s <= ab.length() && within(p - c) {
                return 0.0;
            }
        }
    }
    let mut candidates = vec![n0, n1];
    let nu = Vector2::new(t.y, -t.x);
    for v in [a - c, b - c, nu, -nu] {
        if v.length() > 0.0 && within(v) {
            candidates.push(v.normalized());
        }
    }
    candidates.into_iter().map(|v| point_segment(c + v * r, a, b)).fold(f64::INFINITY, f64::min)
}

/// The least over `r` in [0, `d`] of how much further than `r` the point `c + r w` is from
/// the segment `a`-`b` (convex in `r`: a distance to a convex set along a line, less `r`).
fn bisector_gap(c: Point2, w: Vector2, d: f64, a: Point2, b: Point2) -> f64 {
    let gap = |r: f64| point_segment(c + w * r, a, b) - r;
    let (mut lo, mut hi) = (0.0, d);
    let golden = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..100 {
        let (m1, m2) = (hi - golden * (hi - lo), lo + golden * (hi - lo));
        if gap(m1) <= gap(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    gap(0.0).min(gap(d)).min(gap(0.5 * (lo + hi)))
}

/// Each vertex where a concave edge ends (by corner), as a [`Mitre`]; `None` where one does
/// not end between two convex edges square to it, or the offsets of its faces do not meet
/// beyond them along the convex edges.
fn mitres(points: &[Point], faces: &[Face], edges: &[Edge], d: f64) -> Option<BTreeMap<usize, Mitre>> {
    let mut mitres = BTreeMap::new();
    for e in edges.iter().filter(|e| !e.convex) {
        for corner in [e.start, e.end] {
            let at: Vec<usize> = (0..edges.len()).filter(|&i| edges[i].start == corner || edges[i].end == corner).collect();
            let convex: Vec<usize> = at.iter().copied().filter(|&i| edges[i].convex).collect();
            let [a, b] = convex[..] else { return None };
            if at.len() != 3 || edges[a].axis.dot(e.axis).abs() > SQUARE || edges[b].axis.dot(e.axis).abs() > SQUARE {
                return None;
            }
            let cap = edges[a].faces.into_iter().find(|f| edges[b].faces.contains(f))?;
            let [f, g] = e.faces;
            if e.faces.contains(&cap) {
                return None;
            }
            let (nf, ng) = (faces[f].normal, faces[g].normal);
            let meet = points[corner] + (nf + ng) * (d / (1.0 + nf.dot(ng)));
            // (Beyond the corner along both convex edges, short of their other ends.)
            for i in [a, b] {
                let from = if edges[i].start == corner { edges[i].axis } else { -edges[i].axis };
                let along = (meet - points[corner]).dot(from);
                if !(along > 0.0 && along < edges[i].length) {
                    return None;
                }
            }
            mitres.insert(corner, Mitre { cap, convex: [a, b], meet });
        }
    }
    Some(mitres)
}

/// The grown solid, from its parts (see the module's documentation).
fn build(model: &mut Model, points: &[Point], faces: &[Face], edges: &[Edge], mitres: &BTreeMap<usize, Mitre>, d: f64) -> Result<Shape> {
    let tolerance = crate::tolerance::LINEAR;
    let line2 = |o: (f64, f64), dir: (f64, f64)| Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(dir.0, dir.1)));
    let point_of = |model: &Model, v: VertexId| model.get(v).expect("vertex").point;
    // Each face's corners moved; at a mitre, the faces of its concave edge share where their
    // offsets meet.
    let mut meets: HashMap<usize, VertexId> = HashMap::new();
    for (&corner, mitre) in mitres {
        meets.insert(corner, model.add_vertex(mitre.meet, tolerance));
    }
    let mut moved: HashMap<(usize, usize), VertexId> = HashMap::new();
    for (f, face) in faces.iter().enumerate() {
        for &c in &face.corners {
            let vertex = match mitres.get(&c) {
                Some(mitre) if mitre.cap != f => meets[&c],
                _ => model.add_vertex(points[c] + face.normal * d, tolerance),
            };
            moved.insert((f, c), vertex);
        }
    }
    // Each face's edges moved: by (face, start corner), with its orientation there and
    // length. A concave edge's faces share one.
    let mut lines: HashMap<(usize, usize), (EdgeId, Orientation, f64)> = HashMap::new();
    let line = |model: &mut Model, a: VertexId, b: VertexId| -> Result<(EdgeId, f64)> {
        let (pa, pb) = (point_of(model, a), point_of(model, b));
        let length = pa.distance(pb);
        let edge = model.add_edge(a, b, tolerance)?;
        model.set_edge_curve(edge.into(), Curve::Line(Line::new(pa, Direction::from_vector(pb - pa))), tolerance)?;
        model.set_edge_range(edge.into(), 0.0, length)?;
        Ok((edge, length))
    };
    for e in edges {
        let [f, g] = e.faces;
        if e.convex {
            let (edge, length) = line(model, moved[&(f, e.start)], moved[&(f, e.end)])?;
            lines.insert((f, e.start), (edge, Orientation::Forward, length));
            let (edge, length) = line(model, moved[&(g, e.end)], moved[&(g, e.start)])?;
            lines.insert((g, e.end), (edge, Orientation::Forward, length));
        } else {
            let (edge, length) = line(model, meets[&e.start], meets[&e.end])?;
            lines.insert((f, e.start), (edge, Orientation::Forward, length));
            lines.insert((g, e.end), (edge, Orientation::Reversed, length));
        }
    }
    // Each convex edge's arcs about its ends where they are not mitres, from its first face's
    // normal to its second's: by (edge, end).
    let mut arcs: HashMap<(usize, usize), (EdgeId, Curve)> = HashMap::new();
    for (i, e) in edges.iter().enumerate().filter(|(_, e)| e.convex) {
        let [f, g] = e.faces;
        for corner in [e.start, e.end] {
            if mitres.contains_key(&corner) {
                continue;
            }
            let frame = Frame::new(points[corner], Direction::from_vector(e.axis), Direction::from_vector(faces[f].normal));
            let circle = Curve::Circle(Circle::new(frame, d));
            let edge = model.add_edge(moved[&(f, corner)], moved[&(g, corner)], tolerance)?;
            model.set_edge_curve(edge.into(), circle.clone(), tolerance)?;
            model.set_edge_range(edge.into(), 0.0, e.angle)?;
            arcs.insert((i, corner), (edge, circle));
        }
    }
    // Each mitre's ellipse, from the cap's corner moved to where the concave edge's faces'
    // offsets meet: the cut of the convex edges' cylinders by the plane through the corner
    // equally far from their lines.
    let mut ellipses: HashMap<usize, (EdgeId, Curve, f64, f64)> = HashMap::new();
    for (&corner, mitre) in mitres {
        let v = points[corner];
        let from = |i: usize| if edges[i].start == corner { edges[i].axis } else { -edges[i].axis };
        let (a, b) = (from(mitre.convex[0]), from(mitre.convex[1]));
        let m = (a - b).normalized();
        let minor = a.cross(m).normalized();
        let major = m.cross(minor);
        let (start, end) = (moved[&(mitre.cap, corner)], meets[&corner]);
        let (ps, pe) = (point_of(model, start), point_of(model, end));
        let make = |flip: bool| {
            let z = if flip { -m } else { m };
            let frame = Frame::new(v, Direction::from_vector(z), Direction::from_vector(major));
            let ellipse = Ellipse::new(frame, d / m.dot(a).abs(), d);
            let param = |p: Point| {
                let r = p - v;
                (r.dot(Vector::from(frame.y())) / d).atan2(r.dot(Vector::from(frame.x())) * m.dot(a).abs() / d)
            };
            let (t0, mut t1) = (param(ps), param(pe));
            while t1 < t0 {
                t1 += std::f64::consts::TAU;
            }
            (ellipse, t0, t1)
        };
        let (mut ellipse, mut t0, mut t1) = make(false);
        if t1 - t0 > std::f64::consts::PI {
            (ellipse, t0, t1) = make(true);
        }
        let curve = Curve::Ellipse(ellipse);
        if curve.point(t0).distance(ps) > 1e-9 * (1.0 + d) || curve.point(t1).distance(pe) > 1e-9 * (1.0 + d) {
            return Err(Error::Failed("a mitre's ellipse missing its ends".into()));
        }
        let edge = model.add_edge(start, end, tolerance)?;
        model.set_edge_curve(edge.into(), curve.clone(), tolerance)?;
        model.set_edge_range(edge.into(), t0, t1)?;
        ellipses.insert(corner, (edge, curve, t0, t1));
    }
    let shell = model.new_empty(NodeKind::Shell);
    let uid = model.uid();
    // A face on `surface` bounded by one loop of (edge, orientation, pcurve, range).
    let add_face = |model: &mut Model, surface: Surface, uses: Vec<(EdgeId, Orientation, Curve2, f64, f64)>| -> Result<Shape> {
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, surface, tolerance)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, orientation, pcurve, first, last) in uses {
            let coedge = model.add_coedge(edge, orientation)?;
            model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, wire.into())?;
        model.add(shell, face)?;
        Ok(face)
    };
    let projected = |curve: &Curve, surface: &Surface, first: f64, last: f64| -> Result<Curve2> {
        crate::geometry::project::pcurve_on_range_within(curve, surface, first, last, &[1e-7])
            .map(Curve2::BSpline)
            .ok_or_else(|| Error::Failed("an arc's pcurve on a grown solid's face".into()))
    };
    // The faces moved.
    for (f, face) in faces.iter().enumerate() {
        let x = Direction::from_vector(Vector::from(face.plane.frame().x()));
        let plane = Plane::from_frame(Frame::new(face.plane.origin() + face.normal * d, Direction::from_vector(face.normal), x));
        let surface = Surface::Plane(plane);
        let mut uses = Vec::with_capacity(face.corners.len());
        for &c in &face.corners {
            let (edge, orientation, length) = lines[&(f, c)];
            let (curve, _, _) = model.edge_curve(Shape::from(edge)).expect("a line");
            let pcurve = crate::geometry::project::analytic_pcurve(&curve, &surface).expect("a line on its plane");
            uses.push((edge, orientation, pcurve, 0.0, length));
        }
        let new = add_face(model, surface, uses)?;
        crate::boolean::record(uid, face.shape.node, new.node);
        for (k, &c) in face.corners.iter().enumerate() {
            crate::boolean::record(uid, face.edges[k], lines[&(f, c)].0.into());
        }
    }
    // The convex edges' cylinders: (θ, v) from the first face's normal at the edge's start,
    // anticlockwise (0, 0), (angle, 0), (angle, length), (0, length), short of the ends of
    // a face's side at mitres.
    for (i, e) in edges.iter().enumerate().filter(|(_, e)| e.convex) {
        let [f, g] = e.faces;
        let frame = Frame::new(points[e.start], Direction::from_vector(e.axis), Direction::from_vector(faces[f].normal));
        let surface = Surface::Cylinder(Cylinder::new(frame, d));
        // How far along the edge face `h`'s side starts at `corner`.
        let short = |h: usize, corner: usize| match mitres.get(&corner) {
            Some(mitre) if mitre.cap != h => (mitre.meet - points[corner]).dot(e.axis).abs(),
            _ => 0.0,
        };
        let (along_f, _, _) = lines[&(f, e.start)];
        let (along_g, _, _) = lines[&(g, e.end)];
        let (a, l) = (e.angle, e.length);
        // An end at `corner` (at `v` along the edge), run from face `from`'s side to the
        // other's.
        let end = |model: &Model, corner: usize, from: usize, v: f64| -> Result<(EdgeId, Orientation, Curve2, f64, f64)> {
            match ellipses.get(&corner) {
                Some((edge, curve, t0, t1)) => {
                    let mut pcurve = projected(curve, &surface, *t0, *t1)?;
                    // (Within the cylinder's turn about the edge, not a period on.)
                    if pcurve.point(0.5 * (t0 + t1)).x > a + 1.0 {
                        pcurve = crate::math::Transformable2::transformed(&pcurve, &crate::math::Transform2::translation(Vector2::new(-std::f64::consts::TAU, 0.0)));
                    }
                    let starts_here = model.edge_vertices(Shape::from(*edge)).0.node == NodeId::from(moved[&(from, corner)]);
                    Ok((*edge, if starts_here { Orientation::Forward } else { Orientation::Reversed }, pcurve, *t0, *t1))
                }
                None => {
                    let (edge, _) = &arcs[&(i, corner)];
                    Ok((*edge, if from == f { Orientation::Forward } else { Orientation::Reversed }, line2((0.0, v), (1.0, 0.0)), 0.0, a))
                }
            }
        };
        let (_, _, len_g) = lines[&(g, e.end)];
        let (_, _, len_f) = lines[&(f, e.start)];
        let uses = vec![
            end(model, e.start, f, 0.0)?,
            (along_g, Orientation::Reversed, line2((a, l - short(g, e.end)), (0.0, -1.0)), 0.0, len_g),
            end(model, e.end, g, l)?,
            (along_f, Orientation::Reversed, line2((0.0, short(f, e.start)), (0.0, 1.0)), 0.0, len_f),
        ];
        let new = add_face(model, surface, uses)?;
        crate::boolean::record(uid, e.original, new.node);
    }
    // The spheres of the vertices whose edges are all convex: each arc as its cylinder does
    // not use it (the cylinder runs its start's arc forward, its end's backward), chained
    // round the vertex.
    for (corner, &point) in points.iter().enumerate() {
        if mitres.contains_key(&corner) {
            continue;
        }
        let runs: Vec<(usize, bool)> = (0..edges.len()).filter(|&i| edges[i].start == corner || edges[i].end == corner).map(|i| (i, edges[i].end == corner)).collect();
        if runs.iter().any(|&(i, _)| !edges[i].convex) {
            return Err(Error::Failed("a vertex of a grown solid neither convex nor a mitre".into()));
        }
        // Each run as (from vertex, to vertex) on the sphere.
        let ends = |(i, forward): (usize, bool)| {
            let (a, b) = (moved[&(edges[i].faces[0], corner)], moved[&(edges[i].faces[1], corner)]);
            if forward { (a, b) } else { (b, a) }
        };
        let mut chain = vec![runs[0]];
        while chain.len() < runs.len() {
            let to = ends(*chain.last().expect("a run")).1;
            let Some(&next) = runs.iter().find(|&&r| ends(r).0 == to) else { break };
            chain.push(next);
        }
        if chain.len() < runs.len() || ends(chain[chain.len() - 1]).1 != ends(chain[0]).0 {
            return Err(Error::Failed("a vertex's arcs not in one loop".into()));
        }
        // The sphere's poles and seam away from the patch: its frame's X axis opposite the
        // direction the normals round the vertex spread about (the patch about u = π).
        let mut out = Vector::new(0.0, 0.0, 0.0);
        for &(i, _) in &runs {
            let e = &edges[i];
            out = out - (points[if e.start == corner { e.end } else { e.start }] - point).normalized();
        }
        let out = out.normalized();
        let z = perpendicular(out);
        let surface = Surface::Sphere(Sphere::new(Frame::new(point, Direction::from_vector(z), Direction::from_vector(-out)), d));
        let mut uses = Vec::with_capacity(chain.len());
        for &(i, forward) in &chain {
            let (edge, circle) = &arcs[&(i, corner)];
            let pcurve = projected(circle, &surface, 0.0, edges[i].angle)?;
            uses.push((*edge, if forward { Orientation::Forward } else { Orientation::Reversed }, pcurve, 0.0, edges[i].angle));
        }
        let new = add_face(model, surface, uses)?;
        for &(i, _) in &runs {
            crate::boolean::record(uid, edges[i].original, new.node);
        }
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    Ok(solid)
}

/// A unit vector square to the unit vector `v`.
fn perpendicular(v: Vector) -> Vector {
    let other = if v.x.abs() < 0.6 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    v.cross(other).normalized()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cuboid;
    use std::f64::consts::PI;

    fn valid(model: &Model, s: Shape) {
        let problems: Vec<String> = crate::validity::check(model, s).errors().map(|p| p.description.clone()).collect();
        assert!(problems.is_empty(), "{problems:?}");
    }

    /// Volume, area and face count.
    fn measures(model: &Model, s: Shape) -> (f64, f64, usize) {
        (crate::properties::volume(model, s).mass, crate::properties::area(model, s).mass, model.unique_subshapes(s, NodeKind::Face).len())
    }

    /// `solid` grown by `d` directly and the general way (its faces moved, the edges where
    /// they part filleted): the same solid.
    fn both_ways(model: &mut Model, solid: Shape, d: f64) -> (f64, f64, usize) {
        let direct = polyhedron_grown(model, solid, d).unwrap().expect("built directly");
        valid(model, direct);
        let general = super::super::shape::grown_the_general_way(model, solid, d).unwrap();
        valid(model, general);
        let (a, b) = (measures(model, direct), measures(model, general));
        assert!((a.0 - b.0).abs() <= 1e-9 * a.0, "volumes {} and {}", a.0, b.0);
        assert!((a.1 - b.1).abs() <= 1e-7 * a.1, "areas {} and {}", a.1, b.1);
        assert_eq!(a.2, b.2);
        a
    }

    /// A box a by b by c grown by d: its volume abc + 2d(ab + bc + ca) + πd²(a + b + c) +
    /// 4/3πd³, its area 2(ab + bc + ca) + 2πd(a + b + c) + 4πd², in 6 + 12 + 8 faces.
    #[test]
    fn a_box_grown_is_its_minkowski_sum_with_a_ball() {
        for (a, b, c, d) in [(10.0, 10.0, 10.0, 1.0), (2.34, 6.42, 5.01, 1.02), (1.0, 7.0, 0.3, 4.5)] {
            let mut model = Model::new();
            let frame = Frame::new(Point::new(1.0, -2.0, 0.5), Direction::from_vector(Vector::new(1.0, 2.0, 2.0)), Direction::X);
            let solid = Cuboid::new(a, b, c).with_frame(frame).build(&mut model).unwrap();
            let (v, s, n) = both_ways(&mut model, solid, d);
            let volume = a * b * c + 2.0 * d * (a * b + b * c + c * a) + PI * d * d * (a + b + c) + 4.0 / 3.0 * PI * d.powi(3);
            let area = 2.0 * (a * b + b * c + c * a) + 2.0 * PI * d * (a + b + c) + 4.0 * PI * d * d;
            assert!((v - volume).abs() < 1e-9 * volume, "{v} != {volume}");
            assert!((s - area).abs() < 1e-7 * area, "{s} != {area}");
            assert_eq!(n, 26);
        }
    }

    /// Steiner's formula for a convex polyhedron grown by d: V + Sd + Md² + 4/3πd³, M half
    /// the sum over the edges of their lengths times the angles their faces' normals turn
    /// through.
    fn steiner(model: &Model, solid: Shape, d: f64) -> f64 {
        let (points, faces) = polyhedron(model, solid).unwrap();
        let mut mean = 0.0;
        for (f, face) in faces.iter().enumerate() {
            for other in faces.iter().skip(f + 1) {
                let shared: Vec<usize> = face.corners.iter().copied().filter(|c| other.corners.contains(c)).collect();
                if let [p, q] = shared[..] {
                    mean += 0.5 * points[p].distance(points[q]) * face.normal.angle(other.normal);
                }
            }
        }
        let (v0, s0) = (crate::properties::volume(model, solid).mass, crate::properties::area(model, solid).mass);
        v0 + s0 * d + mean * d * d + 4.0 / 3.0 * PI * d.powi(3)
    }

    #[test]
    fn a_tetrahedron_grown() {
        let mut model = Model::new();
        let corners = [Point::new(0.0, 0.0, 0.0), Point::new(4.0, 0.5, 0.0), Point::new(1.0, 3.0, 0.2), Point::new(1.5, 1.0, 2.5)];
        let solid = crate::blend::tetrahedron(&mut model, corners).unwrap();
        let d = 0.7;
        let (v, _, n) = both_ways(&mut model, solid, d);
        assert_eq!(n, 4 + 6 + 4);
        let want = steiner(&model, solid, d);
        assert!((v - want).abs() < 1e-9 * want, "{v} != {want}");
    }

    /// A pyramid on a rectangle (a box's top collapsed to a point): four faces at its apex.
    #[test]
    fn a_pyramid_grown() {
        let mut model = Model::new();
        let solid = crate::primitives::Wedge::with_top(4.0, 3.0, 5.0, 1.5, 2.0, 1.5, 2.0).build(&mut model).unwrap();
        let d = 0.8;
        let (v, _, n) = both_ways(&mut model, solid, d);
        assert_eq!(n, 5 + 8 + 5);
        let want = steiner(&model, solid, d);
        assert!((v - want).abs() < 1e-9 * want, "{v} != {want}");
    }

    /// The prism of height `h` over the outline `corners` (anticlockwise, in z = 0).
    fn prism(model: &mut Model, corners: &[(f64, f64)], h: f64) -> Shape {
        let n = corners.len();
        let edges: Vec<Shape> = (0..n)
            .map(|i| {
                let ((x0, y0), (x1, y1)) = (corners[i], corners[(i + 1) % n]);
                crate::make::edge_between(model, Point::new(x0, y0, 0.0), Point::new(x1, y1, 0.0)).unwrap()
            })
            .collect();
        let wire = crate::make::wire(model, &edges).unwrap();
        let face = crate::make::planar_face(model, wire).unwrap();
        crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, h)).unwrap().shape
    }

    /// A prism t high over an outline of area A and perimeter P grown by d, its outline grown
    /// by r of area A + Pr + kr² (k half the outer corners' turns less the inner corners'
    /// tan(φ/2)): t(A + Pd + kd²) + 2(Ad + πPd²/4 + 2kd³/3) of volume, t(P + 2kd) +
    /// 2(A + πPd/2 + 2kd²) of area.
    fn prism_grown(a: f64, p: f64, k: f64, t: f64, d: f64) -> (f64, f64) {
        (t * (a + p * d + k * d * d) + 2.0 * (a * d + PI * p * d * d / 4.0 + 2.0 * k * d.powi(3) / 3.0), t * (p + 2.0 * k * d) + 2.0 * (a + PI * p * d / 2.0 + 2.0 * k * d * d))
    }

    /// An L-shaped plate grown: its inner corner's edge stays sharp, the rounds of the edges
    /// beside it at the caps meeting in mitres.
    #[test]
    fn an_l_plate_grown_mitres_its_rounds_over_its_inner_corner() {
        let k = 5.0 * PI / 4.0 - 1.0;
        for (corners, a, p, t, d) in [
            (vec![(0.0, 0.0), (20.0, 0.0), (20.0, 5.0), (5.0, 5.0), (5.0, 15.0), (0.0, 15.0)], 150.0, 70.0, 8.0, 1.0),
            (vec![(0.0, 0.0), (6.474, 0.0), (6.474, 2.209), (2.209, 2.209), (2.209, 5.820), (0.0, 5.820)], 6.474 * 2.209 + 2.209 * (5.820 - 2.209), 2.0 * (6.474 + 5.820), 2.089, 1.593),
        ] {
            let mut model = Model::new();
            let solid = prism(&mut model, &corners, t);
            let (v, s, n) = both_ways(&mut model, solid, d);
            let (volume, area) = prism_grown(a, p, k, t, d);
            assert!((v - volume).abs() < 1e-9 * volume, "{v} != {volume}");
            assert!((s - area).abs() < 1e-7 * area, "{s} != {area}");
            // 8 faces moved, 17 rounds (every edge but the inner one), 10 balls.
            assert_eq!(n, 8 + 17 + 10);
        }
    }

    /// A chevron (inner corner not square) and a cross (four inner corners, sides between
    /// them) grown, against the general way.
    #[test]
    fn prisms_with_inner_corners_grown() {
        let mut model = Model::new();
        let chevron = prism(&mut model, &[(0.0, 0.0), (4.0, 3.0), (8.0, 0.0), (8.0, 2.0), (4.0, 5.0), (0.0, 2.0)], 1.5);
        both_ways(&mut model, chevron, 0.4);
        let cross = prism(&mut model, &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (3.0, 1.0), (3.0, 2.0), (2.0, 2.0), (2.0, 3.0), (1.0, 3.0), (1.0, 2.0), (0.0, 2.0), (0.0, 1.0), (1.0, 1.0)], 0.5);
        both_ways(&mut model, cross, 0.3);
    }

    /// Outlines that do not grow regularly are left to the general way: a slot narrower than
    /// twice the distance closes; an arm shorter than the distance past an inner corner loses
    /// its side.
    #[test]
    fn outlines_growing_irregularly_are_left_to_the_general_way() {
        let mut model = Model::new();
        let u = prism(&mut model, &[(0.0, 0.0), (5.0, 0.0), (5.0, 4.0), (3.0, 4.0), (3.0, 1.0), (2.0, 1.0), (2.0, 4.0), (0.0, 4.0)], 1.0);
        assert!(polyhedron_grown(&mut model, u, 0.6).unwrap().is_none());
        assert!(polyhedron_grown(&mut model, u, 0.4).unwrap().is_some());
        let l = prism(&mut model, &[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 5.0), (0.0, 5.0)], 1.0);
        assert!(polyhedron_grown(&mut model, l, 2.5).unwrap().is_none());
        assert!(polyhedron_grown(&mut model, l, 1.5).unwrap().is_some());
    }

    /// Solids with concave edges not a prism's, or curved faces, are left to the general way.
    #[test]
    fn other_solids_are_left_to_the_general_way() {
        let mut model = Model::new();
        let lower = Cuboid::new(4.0, 4.0, 1.0).build(&mut model).unwrap();
        let upper = Cuboid::new(2.0, 2.0, 2.0).build(&mut model).unwrap();
        let stepped = crate::boolean::union(&mut model, lower, upper).unwrap();
        assert!(polyhedron_grown(&mut model, stepped, 0.5).unwrap().is_none());
        let cylinder = crate::primitives::Cylinder::new(1.0, 2.0).build(&mut model).unwrap();
        assert!(polyhedron_grown(&mut model, cylinder, 0.5).unwrap().is_none());
    }
}
