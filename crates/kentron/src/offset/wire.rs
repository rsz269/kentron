//! Offsets of planar wires.
//!
//! [`planar_wire`] offsets a wire within its plane. Each edge is offset (lines and circles
//! exactly, other curves by interpolation of offset points), arcs fill the gaps at convex
//! corners, and the raw contour so made is cut where it crosses itself: of its pieces, those
//! lying at the full distance from the wire are kept (the others pass closer to the wire,
//! where concave corners or small features overlap) and chained into closed wires. An open
//! wire is offset on both sides, around its ends, into a closed contour.

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Line2};
use crate::geometry::fit;
use crate::geometry::surface::{Plane, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Direction2, Frame2, Point, Point2, Vector, Vector2};
use crate::topology::{Model, NodeKind, Orientation, Shape};

use super::Join;

/// The offset of a planar wire within its plane at `distance`: outwards for a positive
/// distance (inwards for a negative one) from a closed wire; around an open wire at the
/// distance's magnitude. The result is a wire, a compound of wires, or an empty compound if
/// nothing remains (a wire offset inwards by more than its half width).
#[doc(alias = "BRepOffsetAPI_MakeOffset")]
#[doc(alias = "mkoffset")]
pub fn planar_wire(model: &mut Model, wire: Shape, distance: f64, join: Join) -> Result<Shape> {
    let (edges, closed) = read_edges(model, wire)?;
    let flatness = model.explore(wire, NodeKind::Edge).map(|e| model.tolerance(e)).fold(0.0, f64::max);
    let plane = plane_of(&edges, flatness)?;
    let path: Vec<Piece> = edges.iter().map(|e| e.on_plane(&plane)).collect::<Result<_>>()?;
    let mut path = cut_spine(path);
    let scale = path.iter().map(|p| p.start().distance(p.end())).fold(0.0, f64::max).max(1e-9);
    // Within the gaps the wire leaves between its edges' curves.
    let n = path.len();
    let links = if closed { n } else { n - 1 };
    let gaps = (0..links).map(|i| path[i].end().distance(path[(i + 1) % n].start())).fold(0.0, f64::max);
    let tolerance = (1e-7 * scale.max(1.0)).max(gaps);
    // The side: to the right of the path's direction of travel.
    let (path, side_distance) = if closed {
        let ccw = signed_area(&path) > 0.0;
        let right = if ccw == (distance > 0.0) { distance.abs() } else { -distance.abs() };
        (path, right)
    } else {
        // There and back: the right side of the way back is the left side of the wire.
        let back: Vec<Piece> = path.iter().rev().map(Piece::reversed).collect();
        path.extend(back);
        (path, distance.abs())
    };
    if distance == 0.0 {
        return Err(Error::InvalidInput("an offset by zero".into()));
    }
    let raw = raw_contour(&path, side_distance, tolerance, join);
    let loops = trim(&raw, &path, side_distance.abs(), tolerance);
    let loops: Vec<Vec<Piece>> = if join == Join::Intersection { loops.into_iter().map(join_collinear).collect() } else { loops };
    build(model, &plane, &loops)
}

/// The offset of an open planar wire on one side only, at `distance`: to the right of its
/// direction of travel (seen along the plane's normal turned towards +z, else +y, else +x)
/// for a positive distance, to the left for a negative one: an open wire (or several, where
/// the offset breaks), without the ends going round the wire's ends. A closed wire is offset
/// as by [`planar_wire`].
#[doc(alias = "BRepOffsetAPI_MakeOffset::Init")]
#[doc(alias = "openoffset")]
pub fn open_planar_wire(model: &mut Model, wire: Shape, distance: f64, join: Join) -> Result<Shape> {
    let (edges, closed) = read_edges(model, wire)?;
    if closed {
        return planar_wire(model, wire, distance, join);
    }
    if distance == 0.0 {
        return Err(Error::InvalidInput("an offset by zero".into()));
    }
    let flatness = model.explore(wire, NodeKind::Edge).map(|e| model.tolerance(e)).fold(0.0, f64::max);
    let plane = plane_of(&edges, flatness)?;
    let forward: Vec<Piece> = edges.iter().map(|e| e.on_plane(&plane)).collect::<Result<_>>()?;
    let forward = cut_spine(forward);
    let scale = forward.iter().map(|p| p.start().distance(p.end())).fold(0.0, f64::max).max(1e-9);
    let n = forward.len();
    let gaps = (0..n - 1).map(|i| forward[i].end().distance(forward[i + 1].start())).fold(0.0, f64::max);
    let tolerance = (1e-7 * scale.max(1.0)).max(gaps);
    // There and back, as around the wire; then the side wanted, its ends cut off.
    let mut path = forward.clone();
    path.extend(forward.iter().rev().map(Piece::reversed));
    let raw = raw_contour(&path, distance.abs(), tolerance, join);
    let loops = trim(&raw, &path, distance.abs(), tolerance);
    let loops: Vec<Vec<Piece>> = if join == Join::Intersection { loops.into_iter().map(join_collinear).collect() } else { loops };
    let z = Vector::from(plane.normal());
    let upward = if z.z.abs() > 1e-12 { z.z > 0.0 } else if z.y.abs() > 1e-12 { z.y > 0.0 } else { z.x > 0.0 };
    let want_right = (distance > 0.0) == upward;
    let (start, end) = (forward[0].start(), forward[n - 1].end());
    // Where a point lies: beyond an end of the wire (a cap), or on its right.
    let place = |p: Point2| -> (bool, bool) {
        let mut best = (f64::INFINITY, false, false);
        for (k, piece) in forward.iter().enumerate() {
            let samples = 64;
            for i in 0..=samples {
                let s = i as f64 / samples as f64;
                let q = piece.at(s);
                let d = q.distance(p);
                if d < best.0 {
                    let t = piece.tangent(piece.first + (piece.last - piece.first) * s);
                    let at_end = (k == 0 && i == 0) || (k == n - 1 && i == samples);
                    let beyond = at_end && {
                        let away = p - q;
                        let outwards = if k == 0 && i == 0 { t * -1.0 } else { t };
                        away.dot(outwards) > 1e-6 * away.length()
                    };
                    best = (d, beyond, right(t).dot(p - q) > 0.0);
                }
            }
        }
        let _ = (start, end);
        (best.1, best.2)
    };
    let mut chains: Vec<Vec<Piece>> = Vec::new();
    for pieces in loops {
        let keep: Vec<bool> = pieces.iter().map(|p| { let (cap, on_right) = place(p.at(0.5)); !cap && on_right == want_right }).collect();
        // Runs of kept pieces, round the loop from one that is not.
        let Some(first_dropped) = keep.iter().position(|k| !k) else {
            chains.push(pieces);
            continue;
        };
        let m = pieces.len();
        let mut run: Vec<Piece> = Vec::new();
        for i in 1..=m {
            let k = (first_dropped + i) % m;
            if keep[k] {
                run.push(pieces[k].clone());
            } else if !run.is_empty() {
                chains.push(std::mem::take(&mut run));
            }
        }
    }
    build_open(model, &plane, &chains)
}

/// Open wires through `chains` of pieces (each following on from the one before).
fn build_open(model: &mut Model, plane: &Plane, chains: &[Vec<Piece>]) -> Result<Shape> {
    let tolerance = crate::tolerance::LINEAR;
    let mut wires = Vec::new();
    for pieces in chains {
        let curves: Vec<Curve> = pieces.iter().map(|p| p.curve.to_3d(plane.frame())).collect();
        let mut points: Vec<Point> = curves.iter().zip(pieces).map(|(c, p)| c.point(p.first)).collect();
        points.push(curves[curves.len() - 1].point(pieces[pieces.len() - 1].last));
        let vertices: Vec<_> = points.iter().map(|&p| model.add_vertex(p, tolerance)).collect();
        let mut coedges = Vec::with_capacity(pieces.len());
        for (i, (p, curve)) in pieces.iter().zip(&curves).enumerate() {
            let (a, b) = (vertices[i], vertices[i + 1]);
            let reversed = p.last < p.first;
            let (start, end, first, last) = if reversed { (b, a, p.last, p.first) } else { (a, b, p.first, p.last) };
            let edge = model.add_edge(start, end, tolerance)?;
            model.set_edge_curve(edge.into(), curve.clone(), tolerance)?;
            model.set_edge_range(edge.into(), first, last)?;
            let gap = curve.point(first).distance(model.point(Shape::from(start))).max(curve.point(last).distance(model.point(Shape::from(end))));
            if gap > tolerance {
                model.set_edge_tolerance(edge, gap * 1.5)?;
            }
            coedges.push(model.add_coedge(edge, if reversed { Orientation::Reversed } else { Orientation::Forward })?);
        }
        wires.push(Shape::from(model.add_wire(&coedges)?));
    }
    match wires[..] {
        [one] => Ok(one),
        _ => {
            let compound = model.new_empty(NodeKind::Compound);
            for w in wires {
                model.add(compound, w)?;
            }
            Ok(compound)
        }
    }
}

/// The wire's pieces cut at their inflections and curvature extrema, as OCCT's
/// BRepFill_OffsetWire prepares its spine (so the offset has an edge per cut piece); a
/// wire of one edge is cut even where it has none.
fn cut_spine(path: Vec<Piece>) -> Vec<Piece> {
    let forced = path.len() == 1;
    let mut out = Vec::with_capacity(path.len());
    for piece in path {
        let (lo, hi) = (piece.first.min(piece.last), piece.first.max(piece.last));
        let mut cuts = cut_points(&piece.curve, lo, hi);
        // A circle of half a turn or more stays whole.
        if matches!(piece.curve, Curve2::Circle(_)) && hi - lo >= std::f64::consts::PI {
            cuts.clear();
        } else if forced {
            cuts = match cuts.len() {
                // OCCT cuts into three equal parts, but moves its start on after the first
                // cut, which leaves two: a third and two thirds.
                0 => equal_cuts(&piece.curve, lo, hi, 3),
                // Two parts: the longer (by parameter) is halved.
                1 if cuts[0] > 0.5 * (lo + hi) => {
                    let mut c = equal_cuts(&piece.curve, lo, cuts[0], 2);
                    c.push(cuts[0]);
                    c
                }
                1 => {
                    let mut c = vec![cuts[0]];
                    c.extend(equal_cuts(&piece.curve, cuts[0], hi, 2));
                    c
                }
                _ => cuts,
            };
        }
        let mut bounds = vec![lo];
        bounds.extend(cuts);
        bounds.push(hi);
        let parts = bounds.windows(2).map(|w| Piece { first: w[0], last: w[1], ..piece.clone() });
        if piece.last < piece.first {
            let mut parts: Vec<Piece> = parts.map(|p| p.reversed()).collect();
            parts.reverse();
            out.extend(parts);
        } else {
            out.extend(parts);
        }
    }
    out
}

/// Where OCCT's MAT2d_CutCurve cuts a curve over `[lo, hi]`: at its inflections and
/// curvature extrema (as GeomLProp_CurAndInf2d finds them: conics in closed form, B-splines
/// within each span of C3 continuity), skipping those within the confusion tolerance of a
/// cut or the end.
fn cut_points(curve: &Curve2, lo: f64, hi: f64) -> Vec<f64> {
    let basis = match curve {
        Curve2::Trimmed(t) => t.basis(),
        c => c,
    };
    let mut found: Vec<f64> = match basis {
        Curve2::Line(_) | Curve2::Circle(_) => Vec::new(),
        Curve2::Ellipse(_) => {
            let quarters = (0..4).map(|k| k as f64 * std::f64::consts::FRAC_PI_2);
            quarters.map(|q| lo + (q - lo).rem_euclid(std::f64::consts::TAU)).filter(|&t| t <= hi).collect()
        }
        Curve2::Parabola(_) | Curve2::Hyperbola(_) => [0.0].into_iter().filter(|&t| lo <= t && t <= hi).collect(),
        Curve2::BSpline(b) => {
            // Spans of C3 continuity: knots of multiplicity above degree - 3 bound them.
            let p = b.degree();
            let mut bounds = vec![lo];
            for (&k, &m) in b.knots().iter().zip(b.multiplicities()) {
                if lo < k && k < hi && m + 3 > p {
                    bounds.push(k);
                }
            }
            bounds.push(hi);
            bounds.windows(2).flat_map(|w| numeric_features(curve, w[0], w[1])).collect()
        }
        _ => numeric_features(curve, lo, hi),
    };
    found.sort_by(f64::total_cmp);
    found.dedup();
    keep_apart(curve, lo, hi, found.into_iter())
}

/// Curvature extrema and inflections strictly inside `[lo, hi]`, as OCCT's
/// GeomLProp_NumericCurInf2d: sign changes of the signed curvature's derivative (100
/// samples, values below 1e-4 of the interval counting as zero) and of the sine of the
/// angle between the first two derivatives (30 samples, below 1e-6 zero; none where the
/// second derivative is below 1e-4).
fn numeric_features(curve: &Curve2, lo: f64, hi: f64) -> Vec<f64> {
    let slope = |t: f64| {
        let (_, d1, d2, d3) = curve.d3(t);
        let s2 = d1.dot(d1);
        let s3 = s2 * s2.sqrt();
        if s3 * s2 < 1e-290 {
            return 0.0;
        }
        d1.cross(d3) / s3 - 3.0 * d1.cross(d2) * d1.dot(d2) / (s3 * s2)
    };
    let bend = |t: f64| {
        let (_, d1, d2) = curve.d2(t);
        // A second derivative below 1e-4 counts as none (as OCCT's FuncCurNul2d).
        let n = d1.length() * d2.length();
        if d2.length() < 1e-4 || n <= 0.0 { 0.0 } else { d1.cross(d2) / n }
    };
    let mut roots = sign_changes(&slope, lo, hi, 100, 1e-4 * (hi - lo));
    roots.extend(sign_changes(&bend, lo, hi, 30, 1e-6));
    roots
}

/// The roots of `f` where it changes sign between `samples` equal steps over `[lo, hi]`,
/// by bisection; values within `null` of zero have no sign.
fn sign_changes(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64, samples: usize, null: f64) -> Vec<f64> {
    // The ends evaluated as any parameter: at a knot, on the next span (so a curvature
    // extremum at a knot shows as a sign change at the end of the span before it, as in
    // OCCT).
    let ts: Vec<f64> = (0..=samples).map(|i| lo + (hi - lo) * i as f64 / samples as f64).collect();
    let fs: Vec<f64> = ts.iter().map(|&t| f(t)).collect();
    let mut roots = Vec::new();
    for i in 0..samples {
        if fs[i].abs() <= null || fs[i + 1].abs() <= null || (fs[i] < 0.0) == (fs[i + 1] < 0.0) {
            continue;
        }
        let (mut a, mut b, fa) = (ts[i], ts[i + 1], fs[i]);
        for _ in 0..80 {
            let m = 0.5 * (a + b);
            if (f(m) < 0.0) == (fa < 0.0) {
                a = m;
            } else {
                b = m;
            }
        }
        roots.push(0.5 * (a + b));
    }
    roots
}

/// `parts` equal parameter steps from `lo`, cut as OCCT's CutCurve does (each step taken
/// from the last cut).
fn equal_cuts(curve: &Curve2, lo: f64, hi: f64, parts: usize) -> Vec<f64> {
    let step = (hi - lo) / parts as f64;
    let mut cuts = Vec::new();
    let mut from = lo;
    for i in 1..parts {
        let t = from + i as f64 * step;
        if keep_apart(curve, from, hi, std::iter::once(t)).len() == 1 {
            cuts.push(t);
            from = t;
        }
    }
    cuts
}

/// The candidate cuts (increasing) that are clear of the previous cut and of the end.
fn keep_apart(curve: &Curve2, lo: f64, hi: f64, candidates: impl Iterator<Item = f64>) -> Vec<f64> {
    const PARAMETER: f64 = 1e-8;
    const DISTANCE: f64 = 1e-6;
    let (mut from, end) = (lo, curve.point(hi));
    let mut start = curve.point(lo);
    let mut cuts = Vec::new();
    for t in candidates {
        let p = curve.point(t);
        if t - from > PARAMETER && p.distance(start) > DISTANCE {
            if hi - t < PARAMETER || end.distance(p) < DISTANCE {
                break;
            }
            cuts.push(t);
            (from, start) = (t, p);
        }
    }
    cuts
}

/// An edge of the wire, along the wire.
struct WireEdge {
    curve: Curve,
    first: f64,
    last: f64,
}

impl WireEdge {
    fn at(&self, s: f64) -> Point {
        self.curve.point(self.first + (self.last - self.first) * s)
    }

    /// The edge in the plane's parameters.
    fn on_plane(&self, plane: &Plane) -> Result<Piece> {
        let surface = Surface::Plane(*plane);
        // A trimmed curve has its basis's parameters (and the edge its own range).
        let basis = match &self.curve {
            Curve::Trimmed(t) => t.basis(),
            c => c,
        };
        // A curve slightly off the plane (B-spline poles, say) projected onto it first.
        let exact = crate::geometry::project::analytic_pcurve(basis, &surface).or_else(|| {
            let flat = crate::geometry::project::curve_on_plane(basis, plane, plane.normal(), true)?;
            crate::geometry::project::analytic_pcurve(&flat, &surface)
        });
        let curve = match exact {
            Some(c) => c,
            None => {
                // Interpolated through points of the edge, over its own parameters.
                let (lo, hi) = (self.first.min(self.last), self.first.max(self.last));
                let size = self.at(0.0).distance(self.at(0.5)) + self.at(0.5).distance(self.at(1.0));
                let curve = interpolate_adaptively(|t| uv(plane, self.curve.point(t)), lo, hi, 1e-9 * size.max(1.0))
                    .ok_or_else(|| Error::InvalidInput("an edge that cannot be interpolated".into()))?;
                Curve2::BSpline(curve)
            }
        };
        Ok(Piece { curve, first: self.first, last: self.last, origin: 0 })
    }
}

fn uv(plane: &Plane, p: Point) -> Point2 {
    let (u, v) = plane.parameters(p);
    Point2::new(u, v)
}

/// The wire's edges in order along it, and whether it is closed (ends at its first vertex).
fn read_edges(model: &Model, wire: Shape) -> Result<(Vec<WireEdge>, bool)> {
    // (edge, start vertex, end vertex) along each edge's direction in the wire.
    let mut edges = Vec::new();
    for edge in model.explore(wire, NodeKind::Edge) {
        let Some((curve, first, last)) = model.edge_curve(edge) else { continue };
        let (a, b) = model.edge_vertices(edge);
        let (first, last) = if edge.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
        edges.push((WireEdge { curve, first, last }, a, b));
    }
    if edges.is_empty() {
        return Err(Error::InvalidInput("an offset of a wire without edges".into()));
    }
    // In the order they connect, by shared vertices (a wire used reversed lists its edges
    // backwards).
    let same = |a: &Shape, b: &Shape| a.node == b.node && a.placement == b.placement;
    let mut chain = vec![edges.remove(0)];
    while !edges.is_empty() {
        let (head, tail) = (chain[0].1, chain[chain.len() - 1].2);
        // Shared vertices, else vertices within their tolerances of each other.
        let meet = |a: &Shape, b: &Shape| {
            same(a, b) || model.point(*a).distance(model.point(*b)) <= model.tolerance(*a) + model.tolerance(*b)
        };
        // An edge continuing the chain at either end, turned round if it runs the other way
        // (edges used internally have no direction along the wire).
        let fits = |e: &(WireEdge, Shape, Shape), m: &dyn Fn(&Shape, &Shape) -> bool| {
            m(&e.1, &tail) || m(&e.2, &head) || m(&e.2, &tail) || m(&e.1, &head)
        };
        // (Placed by the test that found it: vertices' tolerances may be large enough for
        // all of them to meet.)
        let found = edges.iter().position(|e| fits(e, &same)).map(|i| (i, true)).or_else(|| edges.iter().position(|e| fits(e, &meet)).map(|i| (i, false)));
        let Some((i, shared)) = found else {
            return Err(Error::InvalidInput("a wire whose edges do not connect".into()));
        };
        let joins = |a: &Shape, b: &Shape| if shared { same(a, b) } else { meet(a, b) };
        let (mut e, mut a, mut b) = edges.remove(i);
        if !(joins(&a, &tail) || joins(&b, &head)) {
            std::mem::swap(&mut e.first, &mut e.last);
            std::mem::swap(&mut a, &mut b);
        }
        if joins(&a, &tail) {
            chain.push((e, a, b));
        } else {
            chain.insert(0, (e, a, b));
        }
    }
    let meet = |a: &Shape, b: &Shape| same(a, b) || model.point(*a).distance(model.point(*b)) <= model.tolerance(*a) + model.tolerance(*b);
    let closed = meet(&chain[0].1, &chain[chain.len() - 1].2);
    // Each curve trimmed to where its vertices are (curves may stop short of or run past
    // them, within the vertices' tolerances).
    let mut out = Vec::with_capacity(chain.len());
    for (mut e, a, b) in chain {
        let (pa, pb) = (model.point(a), model.point(b));
        let size = e.at(0.0).distance(e.at(1.0)).max(e.at(0.0).distance(e.at(0.5)));
        let near = 1e-9 * size.max(1.0);
        let (lo, hi) = (e.first.min(e.last), e.first.max(e.last));
        let project = |p: Point| crate::geometry::extrema::nearest_on_curve(p, &e.curve, lo, hi).parameter;
        let first = if e.at(0.0).distance(pa) > near { project(pa) } else { e.first };
        let last = if e.at(1.0).distance(pb) > near { project(pb) } else { e.last };
        if (first - last).abs() > 0.0 {
            (e.first, e.last) = (first, last);
        }
        out.push(e);
    }
    Ok((out, closed))
}

/// The wire's plane: normal along the Newell normal of its points (for a straight wire, any
/// plane through it). The wire must lie in it within `tolerance` (its edges').
fn plane_of(edges: &[WireEdge], tolerance: f64) -> Result<Plane> {
    let points: Vec<Point> = edges.iter().flat_map(|e| (0..16).map(move |i| e.at(i as f64 / 16.0))).chain([edges[edges.len() - 1].at(1.0)]).collect();
    let n = points.len();
    let newell = (0..n).fold(Vector::ZERO, |acc, i| acc + (points[i] - Point::ORIGIN).cross(points[(i + 1) % n] - Point::ORIGIN));
    let size = points.iter().map(|p| p.distance(points[0])).fold(0.0, f64::max).max(1e-12);
    let normal = match Direction::try_from_vector(newell).filter(|_| newell.length() > 1e-9 * size * size) {
        Some(d) => d,
        None => {
            // Straight: a plane containing the line.
            let along = points.iter().map(|&p| p - points[0]).max_by(|a, b| a.length().total_cmp(&b.length())).unwrap_or(Vector::ZERO);
            let along = Direction::try_from_vector(along).ok_or_else(|| Error::InvalidInput("an offset of a wire of no extent".into()))?;
            let helper = if along.z().abs() < 0.9 { Vector::new(0.0, 0.0, 1.0) } else { Vector::new(1.0, 0.0, 0.0) };
            Direction::from_vector(Vector::from(along).cross(helper).cross(Vector::from(along)))
        }
    };
    let plane = Plane::new(points[0], normal);
    if points.iter().any(|&p| plane.distance(p) > (1e-6 * size.max(1.0)).max(2.0 * tolerance)) {
        return Err(Error::InvalidInput("an offset of a wire that is not planar".into()));
    }
    Ok(plane)
}

/// A curve in the plane traversed from `first` to `last` (decreasing if `last < first`).
#[derive(Clone)]
struct Piece {
    curve: Curve2,
    first: f64,
    last: f64,
    /// For raw offset pieces, which one (pieces cut from it keep the number).
    origin: usize,
}

impl Piece {
    fn at(&self, s: f64) -> Point2 {
        self.curve.point(self.first + (self.last - self.first) * s)
    }

    fn start(&self) -> Point2 {
        self.curve.point(self.first)
    }

    fn end(&self) -> Point2 {
        self.curve.point(self.last)
    }

    /// The unit tangent in the direction of travel at parameter `t`.
    fn tangent(&self, t: f64) -> Vector2 {
        let (_, d) = self.curve.d1(t);
        let d = if self.last < self.first { -d } else { d };
        let length = d.length();
        if length > 0.0 { d * (1.0 / length) } else { d }
    }

    fn reversed(&self) -> Piece {
        Piece { first: self.last, last: self.first, ..self.clone() }
    }
}

/// The right-hand normal of a direction of travel.
fn right(t: Vector2) -> Vector2 {
    Vector2::new(t.y, -t.x)
}

/// The signed area enclosed by a closed path (positive counter-clockwise).
fn signed_area(path: &[Piece]) -> f64 {
    let points: Vec<Point2> = path.iter().flat_map(|p| (0..16).map(move |i| p.at(i as f64 / 16.0))).collect();
    let n = points.len();
    (0..n).map(|i| (points[i] - Point2::ORIGIN).cross(points[(i + 1) % n] - Point2::ORIGIN)).sum::<f64>() * 0.5
}

/// The offset of a piece to the right of its direction of travel by `d` (to the left if
/// negative). `None` where it collapses (a circle offset to its centre).
fn offset_piece(piece: &Piece, d: f64) -> Option<Piece> {
    // Relative to the curve's own direction.
    let d_curve = if piece.last < piece.first { -d } else { d };
    let curve = match &piece.curve {
        Curve2::Line(l) => {
            let n = right(Vector2::from(l.direction));
            Curve2::Line(Line2::new(l.origin + n * d_curve, l.direction))
        }
        Curve2::Circle(c) => {
            let f = c.frame();
            let ccw = f.x().x() * f.y().y() - f.x().y() * f.y().x() > 0.0;
            // The right of a counter-clockwise circle is outwards.
            let radius = c.radius() + if ccw { d_curve } else { -d_curve };
            if radius.abs() <= 1e-12 * c.radius().max(1.0) {
                return None;
            }
            if radius > 0.0 {
                Curve2::Circle(Circle2::new(*f, radius))
            } else {
                // Through the centre: the opposite point at each angle.
                let turned = Frame2::from_axes(f.origin(), -f.x(), -f.y()).expect("perpendicular axes");
                Curve2::Circle(Circle2::new(turned, -radius))
            }
        }
        other => {
            let (lo, hi) = (piece.first.min(piece.last), piece.first.max(piece.last));
            let exact = |t: f64| {
                let (p, dt) = other.d1(t);
                let length = dt.length();
                let n = if length > 0.0 { right(dt * (1.0 / length)) } else { Vector2::new(0.0, 0.0) };
                p + n * d_curve
            };
            let fitted = interpolate_adaptively(&exact, lo, hi, d.abs() * 1e-6)?;
            Curve2::BSpline(fitted)
        }
    };
    Some(Piece { curve, first: piece.first, last: piece.last, origin: 0 })
}

/// A B-spline through `f` over `[lo, hi]` (at its own parameters), with samples doubled
/// until the curve stays within `tolerance` of `f` between them.
fn interpolate_adaptively(f: impl Fn(f64) -> Point2, lo: f64, hi: f64, tolerance: f64) -> Option<BSplineCurve2> {
    let mut n = 64;
    loop {
        let params: Vec<f64> = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).collect();
        let points: Vec<Point2> = params.iter().map(|&t| f(t)).collect();
        let curve = fit::interpolate_with_parameters(&points, &params).ok()?;
        let error = params.windows(2).map(|w| {
            let t = 0.5 * (w[0] + w[1]);
            curve.point(t).distance(f(t))
        }).fold(0.0, f64::max);
        if error <= tolerance || n >= 4096 {
            return Some(curve);
        }
        n *= 2;
    }
}

/// Arcs at corners shorter than this are left out (the offset edges meet across the gap
/// within tolerance). Chosen to match OCCT's results, which have no arcs at corners where
/// they would be shorter than about this.
const SHORTEST_ARC: f64 = 1.0e-5;

/// The offsets of the path's pieces, with arcs filling the corners turning away from the
/// offset side; numbered in order along the contour.
fn raw_contour(path: &[Piece], d: f64, tolerance: f64, join: Join) -> Vec<Piece> {
    let mut raw = Vec::new();
    let n = path.len();
    for i in 0..n {
        if let Some(p) = offset_piece(&path[i], d) {
            raw.push(p);
        }
        // The corner at the end of piece i.
        let next = &path[(i + 1) % n];
        let corner = path[i].end();
        let (t1, t2) = (path[i].tangent(path[i].last), next.tangent(next.first));
        let cross = t1.cross(t2);
        let dot = t1.dot(t2);
        let turn = cross.atan2(dot);
        // Turning left opens a gap on the right, and the other way round.
        let opens = if d > 0.0 { turn > 1e-9 } else { turn < -1e-9 };
        // A full reversal (the end of an open wire) opens on both sides: go round the end.
        let reversal = dot < 0.0 && cross.abs() <= 1e-9;
        if !(opens || reversal) {
            continue;
        }
        let a = corner + right(t1) * d;
        // Joined by intersection: the offsets carried on along their tangents to where they
        // meet (exact for lines), unless the corner nearly turns back.
        if join == Join::Intersection && !reversal && turn.abs() < 0.95 * std::f64::consts::PI {
            let reach = d.abs() * (0.5 * turn.abs()).tan();
            let meet = a + t1 * reach;
            for (from, direction) in [(a, t1), (meet, t2)] {
                if reach > SHORTEST_ARC {
                    let line = Line2::new(from, Direction2::from_vector(direction));
                    raw.push(Piece { curve: Curve2::Line(line), first: 0.0, last: reach, origin: 0 });
                }
            }
            continue;
        }
        let x = (a - corner) * (1.0 / d.abs());
        let ccw = if reversal { d > 0.0 } else { turn > 0.0 };
        let sweep = if reversal { std::f64::consts::PI } else { turn.abs() };
        if sweep * d.abs() <= SHORTEST_ARC {
            continue;
        }
        let frame = Frame2::new(corner, Direction2::from_vector(x), ccw);
        raw.push(Piece { curve: Curve2::Circle(Circle2::new(frame, d.abs())), first: 0.0, last: sweep, origin: 0 });
    }
    for (k, p) in raw.iter_mut().enumerate() {
        p.origin = k;
    }
    raw
}

/// A polyline of a piece: parameters and points.
fn sample(piece: &Piece) -> Vec<(f64, Point2)> {
    let (lo, hi) = (piece.first.min(piece.last), piece.first.max(piece.last));
    let n = match &piece.curve {
        Curve2::Line(_) => 1,
        // A few points per knot span within the range.
        Curve2::BSpline(b) => (2 * b.knots().iter().filter(|&&k| k > lo && k < hi).count()).clamp(64, 256),
        _ => 64,
    };
    (0..=n)
        .map(|i| {
            let t = piece.first + (piece.last - piece.first) * i as f64 / n as f64;
            (t, piece.curve.point(t))
        })
        .collect()
}

/// Where two pieces cross: parameters on each, refined by Newton steps from the crossing
/// of their polylines.
fn crossings(a: &Piece, b: &Piece, pa: &[(f64, Point2)], pb: &[(f64, Point2)]) -> Vec<(f64, f64)> {
    let mut found = Vec::new();
    for i in 0..pa.len() - 1 {
        let (p0, p1) = (pa[i].1, pa[i + 1].1);
        let (lo, hi) = (Point2::new(p0.x.min(p1.x), p0.y.min(p1.y)), Point2::new(p0.x.max(p1.x), p0.y.max(p1.y)));
        for j in 0..pb.len() - 1 {
            let (q0, q1) = (pb[j].1, pb[j + 1].1);
            if q0.x.max(q1.x) < lo.x || q0.x.min(q1.x) > hi.x || q0.y.max(q1.y) < lo.y || q0.y.min(q1.y) > hi.y {
                continue;
            }
            let (r, s) = (p1 - p0, q1 - q0);
            let denom = r.cross(s);
            if denom.abs() <= 1e-300 {
                continue;
            }
            let qp = q0 - p0;
            let (u, v) = (qp.cross(s) / denom, qp.cross(r) / denom);
            // Slightly beyond the segments too, so that crossings at sample points are found.
            if !(-1e-9..=1.0 + 1e-9).contains(&u) || !(-1e-9..=1.0 + 1e-9).contains(&v) {
                continue;
            }
            let mut ta = pa[i].0 + (pa[i + 1].0 - pa[i].0) * u;
            let mut tb = pb[j].0 + (pb[j + 1].0 - pb[j].0) * v;
            // Newton on A(ta) - B(tb) = 0.
            for _ in 0..20 {
                let (ca, da) = a.curve.d1(ta);
                let (cb, db) = b.curve.d1(tb);
                // Solve da·x - db·y = -(A - B) by Cramer's rule.
                let (rhs, c2) = (-(ca - cb), -db);
                let det = da.cross(c2);
                if det.abs() <= 1e-300 {
                    break;
                }
                let (x, y) = (rhs.cross(c2) / det, da.cross(rhs) / det);
                ta += x;
                tb += y;
                if x.abs() + y.abs() <= 1e-14 * (1.0 + ta.abs() + tb.abs()) {
                    break;
                }
            }
            let inside = |t: f64, p: &Piece| t >= p.first.min(p.last) - 1e-12 && t <= p.first.max(p.last) + 1e-12;
            if inside(ta, a) && inside(tb, b) && a.curve.point(ta).distance(b.curve.point(tb)) <= 1e-6 {
                found.push((ta.clamp(a.first.min(a.last), a.first.max(a.last)), tb.clamp(b.first.min(b.last), b.first.max(b.last))));
            }
        }
    }
    found
}

/// The path's pieces with their polylines, in a tree for nearest-distance queries.
struct PathIndex<'a> {
    path: &'a [Piece],
    samples: Vec<Vec<(f64, Point2)>>,
    tree: crate::math::Bvh,
}

impl<'a> PathIndex<'a> {
    fn new(path: &'a [Piece]) -> Self {
        let samples: Vec<Vec<(f64, Point2)>> = path.iter().map(sample).collect();
        let boxes: Vec<crate::math::BoundingBox> = samples
            .iter()
            .map(|pts| {
                let b = crate::math::BoundingBox::from_points(pts.iter().map(|&(_, p)| Point::new(p.x, p.y, 0.0)));
                // Chords may miss a little of the curve.
                let span = pts[0].1.distance(pts[pts.len() - 1].1).max(b.size().length());
                b.enlarged(0.05 * span)
            })
            .collect();
        PathIndex { path, samples, tree: crate::math::Bvh::new(&boxes) }
    }

    /// The distance from `p` to the path.
    fn distance(&self, p: Point2) -> f64 {
        let q = Point::new(p.x, p.y, 0.0);
        self.tree.nearest(q, |i| distance_to_piece(&self.path[i], &self.samples[i], p)).map_or(f64::INFINITY, |(_, d)| d)
    }
}

/// The distance from `p` to a piece: from the nearest point of its polyline, refined by
/// Newton steps on (C(t) - p)·C'(t) = 0.
fn distance_to_piece(piece: &Piece, points: &[(f64, Point2)], p: Point2) -> f64 {
    let (mut t, mut d) = (points[0].0, f64::INFINITY);
    for w in points.windows(2) {
        let (a, b) = (w[0].1, w[1].1);
        let ab = b - a;
        let s = if ab.length_squared() > 0.0 { ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0) } else { 0.0 };
        let dq = (a + ab * s).distance(p);
        if dq < d {
            d = dq;
            t = w[0].0 + (w[1].0 - w[0].0) * s;
        }
    }
    let (lo, hi) = (piece.first.min(piece.last), piece.first.max(piece.last));
    for _ in 0..20 {
        let (c, d1, d2) = piece.curve.d2(t);
        let f = (c - p).dot(d1);
        let df = d1.dot(d1) + (c - p).dot(d2);
        if df.abs() <= 1e-300 {
            break;
        }
        let next = (t - f / df).clamp(lo, hi);
        let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
        t = next;
        if done {
            break;
        }
    }
    // The curve's distance (the polyline's only seeds the search: its chords pass inside the
    // curve's bends, closer to points there).
    piece.curve.point(t).distance(p)
}

/// The raw contour cut where it crosses itself, keeping the pieces at the full distance from
/// the path, chained into closed loops.
fn trim(raw: &[Piece], path: &[Piece], d: f64, tolerance: f64) -> Vec<Vec<Piece>> {
    // Slivers cut off by crossings are short beside the chaining tolerance, and beside the
    // distance (a wire's gaps may make the tolerance as long as offsets of its short edges).
    let shortest = (100.0 * tolerance.max(1e-7 * d)).min(1e-2 * d);
    let samples: Vec<Vec<(f64, Point2)>> = raw.iter().map(sample).collect();
    let mut cuts: Vec<Vec<f64>> = raw.iter().map(|p| vec![p.first, p.last]).collect();
    let boxes: Vec<(Point2, Point2)> = samples
        .iter()
        .map(|pts| {
            let (mut lo, mut hi) = (pts[0].1, pts[0].1);
            for &(_, p) in pts {
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
            // Chords may miss a little of the curve.
            let margin = 0.05 * lo.distance(hi) + tolerance;
            (Point2::new(lo.x - margin, lo.y - margin), Point2::new(hi.x + margin, hi.y + margin))
        })
        .collect();
    let overlap = |a: &(Point2, Point2), b: &(Point2, Point2)| a.0.x <= b.1.x && b.0.x <= a.1.x && a.0.y <= b.1.y && b.0.y <= a.1.y;
    // Pieces by the left side of their boxes, for a sweep.
    let mut order: Vec<usize> = (0..raw.len()).collect();
    order.sort_by(|&a, &b| boxes[a].0.x.total_cmp(&boxes[b].0.x));
    let mut pairs = Vec::new();
    for (k, &i) in order.iter().enumerate() {
        pairs.push((i, i));
        for &j in &order[k + 1..] {
            if boxes[j].0.x > boxes[i].1.x {
                break;
            }
            if overlap(&boxes[i], &boxes[j]) {
                pairs.push((i.min(j), i.max(j)));
            }
        }
    }
    pairs.sort_unstable();
    for (i, j) in pairs {
        {
            for (ta, tb) in crossings(&raw[i], &raw[j], &samples[i], &samples[j]) {
                if i == j && (ta - tb).abs() <= 1e-9 * (1.0 + ta.abs()) {
                    continue;
                }
                cuts[i].push(ta);
                cuts[j].push(tb);
            }
        }
    }
    let index = PathIndex::new(path);
    // Pieces at the full distance from the path (else they pass where the offset is not).
    let mut kept: Vec<Piece> = Vec::new();
    for (piece, cuts) in raw.iter().zip(cuts) {
        let mut ts = cuts;
        let piece_length = (0..8).map(|k| piece.at(k as f64 / 8.0).distance(piece.at((k + 1) as f64 / 8.0))).sum::<f64>();
        let increasing = piece.last >= piece.first;
        ts.sort_by(|a, b| if increasing { a.total_cmp(b) } else { b.total_cmp(a) });
        ts.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * (1.0 + a.abs()));
        for w in ts.windows(2) {
            let part = Piece { first: w[0], last: w[1], ..piece.clone() };
            // Slivers cut off by crossings closer than the chaining tolerance are dropped
            // (their neighbours meet across them); whole pieces, however short, are kept,
            // and so are most of one (a short piece clipped at a corner).
            let length = (0..4).map(|k| part.at(k as f64 / 4.0).distance(part.at((k + 1) as f64 / 4.0))).sum::<f64>();
            let sliver = length < 0.5 * piece_length;
            if sliver && length <= shortest {
                continue;
            }
            let probes = [0.25, 0.5, 0.75].map(|s| index.distance(part.at(s)));
            if probes.iter().all(|&dist| dist >= d * (1.0 - 1e-4) - tolerance) {
                kept.push(part);
            }
        }
    }
    let ends: Vec<f64> = raw.iter().map(|p| p.last).collect();
    chain_loops(kept, &ends, 100.0 * tolerance.max(1e-7 * d), 1e-3 * d)
}

/// Kept pieces chained end to start into closed loops; consecutive parts of one raw piece
/// are joined back.
/// `ends[k]`: where raw piece `k` ends (elsewhere its parts end at crossings); ends within
/// `near` meet.
fn chain_loops(mut pieces: Vec<Piece>, ends: &[f64], near: f64, gap: f64) -> Vec<Vec<Piece>> {
    let mut loops = Vec::new();
    // Across wider gaps (where nearly tangent edges of the wire do not quite meet) only when
    // nothing is nearer.
    let far = near.max(gap);
    while !pieces.is_empty() {
        let mut chain = vec![pieces.remove(0)];
        let closed = loop {
            let tail = chain[chain.len() - 1].end();
            let closing = tail.distance(chain[0].start());
            // (Unless another piece goes on from nearer: the chain is not round yet.)
            if closing <= near && !pieces.iter().any(|p| p.start().distance(tail) < 0.5 * closing) {
                break true;
            }
            let last = &chain[chain.len() - 1];
            let continuing = |p: &Piece| p.origin == last.origin && p.first == last.last;
            // At a crossing the contour turns onto the other curve; where a raw piece ends it
            // goes on to the next one.
            if last.last != ends[last.origin] {
                let other = pieces.iter().position(|p| !continuing(p) && p.start().distance(tail) <= near);
                if let Some(i) = other {
                    chain.push(pieces.remove(i));
                    continue;
                }
            }
            if let Some(i) = pieces.iter().position(continuing) {
                chain.push(pieces.remove(i));
                continue;
            }
            let nearest = pieces
                .iter()
                .enumerate()
                .map(|(i, p)| (i, p.start().distance(tail)))
                .min_by(|a, b| a.1.total_cmp(&b.1));
            // A lone piece closes only on itself (a full circle), not across a gap.
            let may_close = chain.len() > 1 && closing <= far;
            match nearest {
                Some((i, d)) if d <= near || (d <= far && (d < closing || chain.len() == 1)) => chain.push(pieces.remove(i)),
                _ if may_close => break true,
                _ => break false,
            }
        };
        if !closed {
            continue; // an open remnant: not part of a closed offset
        }
        // Slivers closing on themselves where curves cross nearly tangentially.
        let length: f64 = chain.iter().map(|p| (0..4).map(|k| p.at(k as f64 / 4.0).distance(p.at((k + 1) as f64 / 4.0))).sum::<f64>()).sum();
        if length <= 10.0 * far {
            continue;
        }
        // Join consecutive parts of one raw piece (also across the loop's start).
        let mut merged: Vec<Piece> = Vec::new();
        for p in chain {
            match merged.last_mut() {
                Some(m) if m.origin == p.origin && m.last == p.first => m.last = p.last,
                _ => merged.push(p),
            }
        }
        if merged.len() > 1 {
            let (first, last) = (&merged[0], &merged[merged.len() - 1]);
            if first.origin == last.origin && last.last == first.first {
                let tail = merged.pop().expect("piece");
                merged[0].first = tail.first;
            }
        }
        loops.push(merged);
    }
    loops
}

/// Wires of the loops' pieces placed back in 3D.
/// A loop with consecutive straight pieces along one line joined (an offset carried on to
/// a corner, and the corner's extension, make one edge).
fn join_collinear(mut pieces: Vec<Piece>) -> Vec<Piece> {
    let line_of = |p: &Piece| -> Option<(Point2, Vector2)> {
        let Curve2::Line(_) = p.curve else { return None };
        let (a, b) = (p.start(), p.end());
        let d = b - a;
        let length = d.length();
        (length > 0.0).then(|| (a, d * (1.0 / length)))
    };
    let mut i = 0;
    while pieces.len() > 2 && i < pieces.len() {
        let j = (i + 1) % pieces.len();
        let merged = match (line_of(&pieces[i]), line_of(&pieces[j])) {
            (Some((a, da)), Some((_, db))) if da.cross(db).abs() <= 1e-12 && da.dot(db) > 0.0 => {
                let end = pieces[j].end();
                let length = a.distance(end);
                Some(Piece { curve: Curve2::Line(Line2::new(a, Direction2::from_vector(da))), first: 0.0, last: length, origin: pieces[i].origin })
            }
            _ => None,
        };
        match merged {
            Some(p) => {
                pieces[i] = p;
                pieces.remove(j);
                if j < i {
                    i -= 1;
                }
            }
            None => i += 1,
        }
    }
    pieces
}

fn build(model: &mut Model, plane: &Plane, loops: &[Vec<Piece>]) -> Result<Shape> {
    let tolerance = crate::tolerance::LINEAR;
    let mut wires = Vec::new();
    for pieces in loops {
        let curves: Vec<Curve> = pieces.iter().map(|p| p.curve.to_3d(plane.frame())).collect();
        let starts: Vec<_> = pieces
            .iter()
            .zip(&curves)
            .map(|(p, c)| model.add_vertex(c.point(p.first), tolerance))
            .collect();
        let mut coedges = Vec::with_capacity(pieces.len());
        for (i, (p, curve)) in pieces.iter().zip(&curves).enumerate() {
            let (a, b) = (starts[i], starts[(i + 1) % starts.len()]);
            let reversed = p.last < p.first;
            let (start, end, first, last) = if reversed { (b, a, p.last, p.first) } else { (a, b, p.first, p.last) };
            let edge = model.add_edge(start, end, tolerance)?;
            model.set_edge_curve(edge.into(), curve.clone(), tolerance)?;
            model.set_edge_range(edge.into(), first, last)?;
            // The gap between the curve's ends and the shared vertices.
            let gap = curve.point(first).distance(model.point(Shape::from(start))).max(curve.point(last).distance(model.point(Shape::from(end))));
            if gap > tolerance {
                model.set_edge_tolerance(edge, gap * 1.5)?;
                for v in [start, end] {
                    let point = model.point(Shape::from(v));
                    let t = model.tolerance(Shape::from(v)).max(gap * 1.5);
                    model.set_vertex_point(v.into(), point, t)?;
                }
            }
            coedges.push(model.add_coedge(edge, if reversed { Orientation::Reversed } else { Orientation::Forward })?);
        }
        wires.push(Shape::from(model.add_wire(&coedges)?));
    }
    if wires.len() == 1 {
        return Ok(wires[0]);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for w in wires {
        model.add(compound, w)?;
    }
    Ok(compound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::make;
    use crate::properties;

    fn polyline(model: &mut Model, points: &[[f64; 2]], close: bool) -> Shape {
        let p: Vec<Point> = points.iter().map(|c| Point::new(c[0], c[1], 0.0)).collect();
        let n = if close { p.len() } else { p.len() - 1 };
        let edges: Vec<Shape> = (0..n).map(|i| make::edge_between(model, p[i], p[(i + 1) % p.len()]).unwrap()).collect();
        make::wire(model, &edges).unwrap()
    }

    fn counts(model: &Model, s: Shape) -> (usize, usize) {
        (model.unique_subshapes(s, NodeKind::Edge).len(), model.unique_subshapes(s, NodeKind::Wire).len())
    }

    #[test]
    fn a_square_offset_outwards_gets_rounded_corners() {
        let mut model = Model::new();
        let square = polyline(&mut model, &[[0., 0.], [10., 0.], [10., 10.], [0., 10.]], true);
        let out = planar_wire(&mut model, square, 1.0, Join::Arc).unwrap();
        let length = properties::length(&model, out).mass;
        assert!((length - (40.0 + 2.0 * std::f64::consts::PI)).abs() < 1e-6, "length {length}");
        assert_eq!(counts(&model, out), (8, 1));
        assert!(crate::validity::check(&model, out).errors().next().is_none());
    }

    #[test]
    fn corners_joined_by_intersection_stay_sharp() {
        // An L shape offset outwards with sharp corners: 6 edges, as OCCT gives.
        let mut model = Model::new();
        let l = polyline(&mut model, &[[0., 0.], [10., 0.], [10., 3.], [3., 3.], [3., 8.], [0., 8.]], true);
        let out = planar_wire(&mut model, l, 1.5, Join::Intersection).unwrap();
        assert_eq!(counts(&model, out), (6, 1));
        // 13 x 11 less the 8.5 x 6.5 notch.
        assert!((properties::length(&model, out).mass - 2.0 * (13.0 + 11.0)).abs() < 1e-9);
        assert!(crate::validity::check(&model, out).errors().next().is_none());
    }

    #[test]
    fn a_square_offset_inwards_shrinks() {
        let mut model = Model::new();
        let square = polyline(&mut model, &[[0., 0.], [10., 0.], [10., 10.], [0., 10.]], true);
        let inside = planar_wire(&mut model, square, -1.0, Join::Arc).unwrap();
        assert!((properties::length(&model, inside).mass - 32.0).abs() < 1e-6);
        assert_eq!(counts(&model, inside), (4, 1));
        // Inwards by more than half the width, nothing remains.
        let gone = planar_wire(&mut model, square, -6.0, Join::Arc).unwrap();
        assert_eq!(counts(&model, gone), (0, 0));
    }

    #[test]
    fn a_concave_notch_is_trimmed() {
        let mut model = Model::new();
        // An L shape: offsetting outwards trims the offsets meeting in its inner corner.
        let l = polyline(&mut model, &[[0., 0.], [10., 0.], [10., 4.], [4., 4.], [4., 10.], [0., 10.]], true);
        let out = planar_wire(&mut model, l, 1.0, Join::Arc).unwrap();
        // The outline grows by quarter circles at the five convex corners; the concave corner
        // cuts 1 off each of its two sides.
        let expected = 40.0 + 5.0 * std::f64::consts::FRAC_PI_2 - 2.0;
        let length = properties::length(&model, out).mass;
        assert!((length - expected).abs() < 1e-6, "length {length}, expected {expected}");
        assert_eq!(counts(&model, out), (11, 1));
    }

    #[test]
    fn an_open_wire_is_offset_around() {
        let mut model = Model::new();
        let segment = polyline(&mut model, &[[0., 0.], [10., 0.]], false);
        let around = planar_wire(&mut model, segment, 1.0, Join::Arc).unwrap();
        let length = properties::length(&model, around).mass;
        assert!((length - (20.0 + 2.0 * std::f64::consts::PI)).abs() < 1e-6, "length {length}");
        // A wire of one edge is cut in two first (as OCCT does): two pieces on each side
        // and an arc round each end.
        assert_eq!(counts(&model, around), (6, 1));
    }
}
