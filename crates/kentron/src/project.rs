//! Projecting edges onto faces.
//!
//! [`normal`] projects along the faces' normals: each point of an edge goes to the foot of
//! the perpendicular dropped from it to a face's surface. Following the feet along the edge
//! gives a curve on the surface, kept where the feet lie inside the face (and within a
//! distance, if one is given):
//!
//! 1. The edge is sampled and each sample's foot found by Newton steps from the last one's
//!    (searching the whole surface where that fails), a foot counting where the line to it
//!    is normal to the surface: not where the steps stopped at the surface's bounds.
//! 2. Runs of samples whose feet count make the projected pieces, their ends found by
//!    bisection between the last sample that counts and the first that does not.
//! 3. Each piece is fitted, through feet at more and more samples until the fit passes
//!    within the tolerance of the feet halfway between them, and becomes an edge of its own
//!    (sharing no vertices with the others).

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::extrema::{nearest_on_surface, newton_foot};
use crate::geometry::fit;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Transform, Transformable};
use crate::topology::{CoEdgeId, Model, NodeKind, Shape};
use crate::tolerance;

/// How [`normal`] projects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalProjection {
    /// How far the projected curves may stray from the feet of the perpendiculars.
    pub tolerance: f64,
    /// An edge is projected onto a face only if it comes this near the face's surface and
    /// each piece of its projection starts this near the edge; a piece then runs on,
    /// however far.
    pub max_distance: Option<f64>,
    /// Whether projections end at the faces' boundaries; otherwise they run on over the
    /// faces' whole surfaces.
    pub within_faces: bool,
}

impl Default for NormalProjection {
    fn default() -> Self {
        NormalProjection { tolerance: 1e-4, max_distance: None, within_faces: true }
    }
}

/// The edges of `shapes` projected onto the faces of `target` along their normals: a
/// compound of edges, one per piece of an edge's projection onto a face (see the
/// [module](self) documentation).
#[doc(alias = "BRepOffsetAPI_NormalProjection")]
#[doc(alias = "BRepAlgo_NormalProjection")]
pub fn normal(model: &mut Model, shapes: &[Shape], target: Shape, options: &NormalProjection) -> Result<Shape> {
    let faces: Vec<Target> = model.unique_subshapes(target, NodeKind::Face).into_iter().filter_map(|f| Target::of(model, f)).collect();
    if faces.is_empty() {
        return Err(Error::InvalidInput("nothing to project onto: no faces".into()));
    }
    let mut curves: Vec<(Curve, f64, f64)> = Vec::new();
    for &shape in shapes {
        for e in model.unique_subshapes(shape, NodeKind::Edge) {
            if let Some(curve) = model.edge_curve(e) {
                curves.push(curve);
            }
        }
    }
    let tol = options.tolerance.max(tolerance::LINEAR);
    let mut pieces: Vec<(Curve, f64, f64)> = Vec::new();
    for (curve, first, last) in &curves {
        for face in &faces {
            for (t0, t1, seed) in face.runs(curve, *first, *last, options) {
                if let Some(piece) = face.fit(curve, t0, t1, seed, tol) {
                    pieces.push(piece);
                }
            }
        }
    }
    let result = model.new_empty(NodeKind::Compound);
    for (curve, first, last) in pieces {
        let (a, b) = (curve.point(first), curve.point(last));
        if a.distance(b) <= tol && curve_length(&curve, first, last) <= tol {
            continue;
        }
        let start = model.add_vertex(a, tol);
        // A closed piece (an edge's whole loop) starts and ends at one vertex.
        let end = if a.distance(b) <= tol { start } else { model.add_vertex(b, tol) };
        let edge = model.add_edge(start, end, tol)?;
        model.set_edge_curve(edge.into(), curve, tol)?;
        model.set_edge_range(edge.into(), first, last)?;
        model.add(result, edge.into())?;
    }
    Ok(result)
}

/// A rough length of `curve` over `[first, last]` (of a polyline through samples).
fn curve_length(curve: &Curve, first: f64, last: f64) -> f64 {
    const N: usize = 16;
    let points: Vec<Point> = (0..=N).map(|i| curve.point(first + (last - first) * i as f64 / N as f64)).collect();
    points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// A face projected onto: its surface where it is, and its boundary in its parameters.
struct Target {
    surface: Surface,
    /// The face's loops, sampled in `(u, v)`, for locating feet.
    loops: Vec<Vec<Point2>>,
    u_period: Option<f64>,
    v_period: Option<f64>,
    /// The size of the surface's bounds, scaling what counts as on it.
    size: f64,
    /// How near the loops (in parameters) counts as on them.
    uv_tolerance: f64,
}

/// The foot of a perpendicular: its parameters, the point and its distance.
#[derive(Clone, Copy)]
struct Foot {
    uv: (f64, f64),
    distance: f64,
    /// Whether the line to it is normal to the surface (rather than ending at its bounds).
    normal: bool,
}

impl Target {
    fn of(model: &Model, face: Shape) -> Option<Target> {
        let surface = model.surface(face)?;
        let surface = if face.placement == Transform::IDENTITY { surface.clone() } else { surface.transformed(&face.placement) };
        let mut loops: Vec<Vec<Point2>> = Vec::new();
        for coedge in model.unique_subshapes(face, NodeKind::CoEdge) {
            let Ok(id) = CoEdgeId::try_from(coedge.node) else { continue };
            let Some((pcurve, first, last)) = model.effective_pcurve(id) else { continue };
            const N: usize = 256;
            loops.push((0..=N).map(|i| pcurve.point(first + (last - first) * i as f64 / N as f64)).collect());
        }
        let (u0, u1, v0, v1) = surface.bounds();
        let finite = |x: f64| if x.is_finite() { x } else { 0.0 };
        let corners = [(u0, v0), (u1, v1), (u0, v1), (u1, v0)].map(|(u, v)| surface.point(finite(u), finite(v)));
        let size = corners[0].distance(corners[1]) + corners[2].distance(corners[3]);
        let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
        for q in loops.iter().flatten() {
            (lo.x, lo.y, hi.x, hi.y) = (lo.x.min(q.x), lo.y.min(q.y), hi.x.max(q.x), hi.y.max(q.y));
        }
        let uv_tolerance = if loops.is_empty() { 0.0 } else { 1e-9 * lo.distance(hi) };
        Some(Target { u_period: surface.u_period(), v_period: surface.v_period(), surface, loops, size, uv_tolerance })
    }

    /// The foot from `p` found afresh: the nearest point of the surface if the line to it is
    /// normal there, else the nearest of the feet found from a grid over the surface. Feet
    /// are wherever the line from `p` is normal to the surface: nearest points, but also
    /// saddles and farthest points (the other side of a cylinder), as a foot followed along a
    /// curve may become. `near` hints where the nearest point is.
    fn foot(&self, p: Point, near: Option<(f64, f64)>) -> Option<Foot> {
        let (u, v, _) = nearest_on_surface(p, &self.surface, near);
        if !(u.is_finite() && v.is_finite()) {
            return None;
        }
        let (point, du, dv) = self.surface.d1(u, v);
        let nearest = self.judge(p, (u, v), point, du, dv);
        if nearest.normal {
            return Some(nearest);
        }
        // Else the nearest of the feet reached by steps from a grid over the surface.
        let (u0, u1, v0, v1) = self.finite_bounds();
        const N: usize = 6;
        let mut best: Option<Foot> = None;
        for i in 0..=N {
            for j in 0..=N {
                let start = (u0 + (u1 - u0) * i as f64 / N as f64, v0 + (v1 - v0) * j as f64 / N as f64);
                let f = self.stationary(p, start);
                if f.normal && best.is_none_or(|b| f.distance < b.distance) {
                    best = Some(f);
                }
            }
        }
        Some(best.unwrap_or(nearest))
    }

    /// The foot from `p` followed from `near`, the foot from `from`: by Newton steps to where
    /// the distance is stationary, else (from the bounds, where the last foot stopped, toward
    /// the surface's inside) to the nearest point. `None` where the steps jump onto another
    /// branch (as across a cylinder): along one the distance changes no faster than the
    /// point moves (its derivative is that of the point toward the foot).
    fn follow(&self, p: Point, near: (f64, f64), from: Point) -> Option<Foot> {
        let before = self.surface.point(near.0, near.1).distance(from);
        let reach = 1.5 * p.distance(from) + 1e-7 * (1.0 + self.size);
        let near_enough = |f: &Foot| (f.distance - before).abs() <= reach;
        let f = self.stationary(p, near);
        if f.normal && near_enough(&f) {
            return Some(f);
        }
        let f = newton_foot(p, &self.surface, near);
        Some(self.judge(p, (f.u, f.v), f.point, f.du, f.dv)).filter(near_enough)
    }

    /// The foot at `curve(to)` followed from `uv`, the foot at `curve(from)`: in one step, or
    /// where that jumps, in halves (down to `depth` times), as where the foot runs fast.
    fn trace(&self, curve: &Curve, from: f64, uv: (f64, f64), to: f64, depth: u32) -> Option<Foot> {
        if let Some(f) = self.follow(curve.point(to), uv, curve.point(from)) {
            return Some(f);
        }
        if depth == 0 {
            return None;
        }
        let middle = 0.5 * (from + to);
        let m = self.trace(curve, from, uv, middle, depth - 1)?;
        self.trace(curve, middle, m.uv, to, depth - 1)
    }

    /// The surface's bounds, its infinite directions cut near the origin.
    fn finite_bounds(&self) -> (f64, f64, f64, f64) {
        let (u0, u1, v0, v1) = self.surface.bounds();
        let finite = |a: f64, b: f64| if a.abs() < 1e9 && b.abs() < 1e9 { (a, b) } else { (-1e3, 1e3) };
        let ((u0, u1), (v0, v1)) = (finite(u0, u1), finite(v0, v1));
        (u0, u1, v0, v1)
    }

    /// Where Newton steps from `start` toward a point of the surface whose normal passes
    /// through `p` end (any such point: the steps solve for where the distance is
    /// stationary), kept within the bounds.
    fn stationary(&self, p: Point, start: (f64, f64)) -> Foot {
        let (u0, u1, v0, v1) = self.finite_bounds();
        let (u_periodic, v_periodic) = (self.surface.is_u_periodic(), self.surface.is_v_periodic());
        let clamp = |u: f64, v: f64| (if u_periodic { u } else { u.clamp(u0, u1) }, if v_periodic { v } else { v.clamp(v0, v1) });
        let (mut u, mut v) = clamp(start.0, start.1);
        let (u_most, v_most) = (0.25 * (u1 - u0), 0.25 * (v1 - v0));
        for _ in 0..40 {
            let d = self.surface.d2(u, v);
            let r = d.point - p;
            let (gu, gv) = (r.dot(d.du), r.dot(d.dv));
            let (a, b, c) = (d.du.dot(d.du) + r.dot(d.duu), d.du.dot(d.dv) + r.dot(d.duv), d.dv.dot(d.dv) + r.dot(d.dvv));
            let det = a * c - b * b;
            if !(det.abs() > 1e-300) {
                break;
            }
            let (su, sv) = ((-gu * c + gv * b) / det, (-gv * a + gu * b) / det);
            // Steps of at most a quarter of the range, keeping their direction.
            let scale = (su.abs() / u_most).max(sv.abs() / v_most).max(1.0);
            let (nu, nv) = clamp(u + su / scale, v + sv / scale);
            let done = (nu - u).abs() <= 1e-13 * (u1 - u0) && (nv - v).abs() <= 1e-13 * (v1 - v0);
            (u, v) = (nu, nv);
            if done {
                break;
            }
        }
        let (point, du, dv) = self.surface.d1(u, v);
        self.judge(p, (u, v), point, du, dv)
    }

    fn judge(&self, p: Point, uv: (f64, f64), point: Point, du: crate::math::Vector, dv: crate::math::Vector) -> Foot {
        let r = p - point;
        let distance = r.length();
        let square = |d: crate::math::Vector| d.length() <= 1e-12 || r.dot(d).abs() <= 1e-6 * distance * d.length();
        let normal = distance <= 1e-9 * (1.0 + self.size) || (square(du) && square(dv));
        Foot { uv, distance, normal }
    }

    /// Whether the parameters lie inside the face's loops (even-odd).
    fn inside(&self, (u, v): (f64, f64)) -> bool {
        if self.loops.is_empty() {
            return true;
        }
        // On the boundary (within a billionth of the face's size in parameters) counts as in.
        let on = |a: Point2, b: Point2, q: Point2| {
            let (ab, aq) = (b - a, q - a);
            let s = (aq.dot(ab) / ab.dot(ab).max(1e-300)).clamp(0.0, 1.0);
            (aq - ab * s).length() <= self.uv_tolerance
        };
        if self.loops.iter().any(|l| l.windows(2).any(|w| on(w[0], w[1], Point2::new(u, v)))) {
            return true;
        }
        let shifts = |period: Option<f64>| match period {
            Some(p) => vec![0.0, -p, p],
            None => vec![0.0],
        };
        for su in shifts(self.u_period) {
            for sv in shifts(self.v_period) {
                let (u, v) = (u + su, v + sv);
                let mut crossings = 0usize;
                for polygon in &self.loops {
                    for w in polygon.windows(2) {
                        let (a, b) = (w[0], w[1]);
                        if (a.y > v) != (b.y > v) && a.x + (v - a.y) / (b.y - a.y) * (b.x - a.x) > u {
                            crossings += 1;
                        }
                    }
                }
                if crossings % 2 == 1 {
                    return true;
                }
            }
        }
        false
    }

    /// Whether a foot counts: normal to the surface and (if asked) in the face.
    fn counts(&self, foot: &Foot, options: &NormalProjection) -> bool {
        foot.normal && (!options.within_faces || self.inside(foot.uv))
    }

    /// The parameter ranges of `curve` over `[first, last]` whose feet count, each with the
    /// foot at its start: runs of samples whose feet count, each followed from the last
    /// (none if the curve stays farther than the maximum distance, or a run starts so).
    fn runs(&self, curve: &Curve, first: f64, last: f64, options: &NormalProjection) -> Vec<(f64, f64, (f64, f64))> {
        const N: usize = 128;
        let ts: Vec<f64> = (0..=N).map(|i| first + (last - first) * i as f64 / N as f64).collect();
        const DEPTH: u32 = 10;
        // The last sample's foot.
        let mut previous: Option<(f64, (f64, f64))> = None;
        // The feet that count, and whether each was followed from the sample before.
        let mut feet: Vec<Option<Foot>> = Vec::with_capacity(ts.len());
        let mut followed: Vec<bool> = Vec::with_capacity(ts.len());
        let mut feet_any: Vec<Option<Foot>> = Vec::with_capacity(ts.len());
        // Feet are searched for afresh where one followed is lost, and every few samples
        // (feet entering over the bounds are followed from where the last stopped there).
        const SEARCH_EVERY: usize = 8;
        let mut searched = 0;
        // Whether the last foot was followed (and is lost if the next is not).
        let mut following = false;
        for (i, &t) in ts.iter().enumerate() {
            let traced = previous.and_then(|(from, uv)| self.trace(curve, from, uv, t, DEPTH));
            let (foot, joined) = match traced {
                Some(f) if f.normal => (Some(f), true),
                traced if following || previous.is_none() || i - searched >= SEARCH_EVERY => {
                    searched = i;
                    (self.foot(curve.point(t), traced.or(previous.map(|(_, uv)| Foot { uv, distance: 0.0, normal: false })).map(|f| f.uv)), false)
                }
                traced => (traced, false),
            };
            following = joined;
            previous = foot.map(|f| (t, f.uv));
            followed.push(joined);
            feet_any.push(foot);
            feet.push(foot.filter(|f| self.counts(f, options)));
        }
        // Following a foot that counts at `good` toward `bad`, where it stops counting.
        let edge_of = |good: f64, mut uv: (f64, f64), bad: f64| -> (f64, (f64, f64)) {
            let (mut good, mut bad) = (good, bad);
            for _ in 0..60 {
                if (good - bad).abs() <= 1e-13 * (1.0 + (last - first).abs()) {
                    break;
                }
                let mid = 0.5 * (good + bad);
                match self.trace(curve, good, uv, mid, DEPTH).filter(|f| self.counts(f, options)) {
                    Some(f) => {
                        good = mid;
                        uv = f.uv;
                    }
                    None => bad = mid,
                }
            }
            (good, uv)
        };
        if let Some(d) = options.max_distance {
            let nearest = feet_any.iter().flatten().map(|f| f.distance).fold(f64::INFINITY, f64::min);
            if nearest > d {
                return Vec::new();
            }
        }
        let mut runs = Vec::new();
        let mut i = 0;
        while i < ts.len() {
            let Some(f) = feet[i] else {
                i += 1;
                continue;
            };
            let mut j = i;
            while j + 1 < ts.len() && feet[j + 1].is_some() && followed[j + 1] {
                j += 1;
            }
            // A foot found afresh may have been there samples before: followed back to them.
            let (mut k, mut uv) = (i, f.uv);
            while k > 0 && !followed[k] {
                match self.trace(curve, ts[k], uv, ts[k - 1], DEPTH).filter(|g| g.normal && self.counts(g, options)) {
                    Some(g) => {
                        k -= 1;
                        uv = g.uv;
                    }
                    None => break,
                }
            }
            let (start, uv) = if k > 0 { edge_of(ts[k], uv, ts[k - 1]) } else { (ts[0], uv) };
            let end = if j + 1 < ts.len() { edge_of(ts[j], feet[j].unwrap().uv, ts[j + 1]).0 } else { ts[j] };
            if end > start {
                // A piece starting too far makes none count.
                let from = curve.point(start);
                if options.max_distance.is_some_and(|d| self.surface.point(uv.0, uv.1).distance(from) > d) {
                    return Vec::new();
                }
                runs.push((start, end, uv));
            }
            i = j + 1;
        }
        runs
    }

    /// The projection of `curve` over `[t0, t1]` as a B-spline in the same parameter, its
    /// feet followed from `seed`.
    fn fit(&self, curve: &Curve, t0: f64, t1: f64, seed: (f64, f64), tolerance: f64) -> Option<(Curve, f64, f64)> {
        // The foot at `t` followed from `near`, the foot at `from` (by Newton steps from it
        // even where they jump: the run's feet do not, but the steps may stray).
        let follow = |t: f64, near: (f64, f64), from: f64| {
            let p = curve.point(t);
            let f = self.trace(curve, from, near, t, 10).unwrap_or_else(|| {
                let f = newton_foot(p, &self.surface, near);
                Foot { uv: (f.u, f.v), distance: 0.0, normal: false }
            });
            (f.uv, self.surface.point(f.uv.0, f.uv.1))
        };
        let mut intervals = 4;
        let mut best = None;
        while intervals <= 2048 {
            let params: Vec<f64> = (0..=intervals).map(|i| t0 + (t1 - t0) * i as f64 / intervals as f64).collect();
            let (mut near, mut from) = (seed, t0);
            let mut uvs = Vec::with_capacity(params.len());
            let mut points = Vec::with_capacity(params.len());
            for &t in &params {
                let (uv, p) = follow(t, near, from);
                (near, from) = (uv, t);
                uvs.push(uv);
                points.push(p);
            }
            let fitted = fit::interpolate_with_parameters(&points, &params).ok()?;
            let error = params
                .windows(2)
                .zip(&uvs)
                .map(|(w, &uv)| {
                    let t = 0.5 * (w[0] + w[1]);
                    fitted.point(t).distance(follow(t, uv, w[0]).1)
                })
                .fold(0.0, f64::max);
            best = Some(Curve::BSpline(fitted));
            if error <= tolerance {
                break;
            }
            intervals *= 2;
        }
        best.map(|c| (c, t0, t1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cylinder;

    #[test]
    fn a_line_onto_a_cylinder() {
        // A line above a cylinder of radius 5 and height 10, along its axis, projects onto
        // its side (and nowhere onto its caps, whose normals do not reach it).
        let mut model = Model::new();
        let cylinder = Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
        let line = crate::make::edge_between(&mut model, Point::new(8.0, 0.0, -3.0), Point::new(8.0, 0.0, 13.0)).unwrap();
        let r = normal(&mut model, &[line], cylinder, &NormalProjection::default()).unwrap();
        let edges = model.unique_subshapes(r, NodeKind::Edge);
        assert_eq!(edges.len(), 1);
        let length = crate::properties::length(&model, r).mass;
        assert!((length - 10.0).abs() < 1e-6, "{length}");
        let (curve, a, b) = model.edge_curve(edges[0]).unwrap();
        for t in [a, 0.5 * (a + b), b] {
            let p = curve.point(t);
            assert!(((p.x * p.x + p.y * p.y).sqrt() - 5.0).abs() < 1e-4 && p.y.abs() < 1e-6, "{p:?}");
        }
    }

    #[test]
    fn a_ring_onto_a_cylinder_and_a_maximum_distance() {
        // A circle of radius 8 around a cylinder of radius 5 projects onto a closed circle of
        // radius 5 (one edge on one vertex), but not within a distance 2 of it (it keeps 3
        // away).
        use crate::geometry::curve::Circle;
        use crate::math::{Direction, Frame};
        let mut model = Model::new();
        let cylinder = Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
        let circle = Circle::new(Frame::new(Point::new(0.0, 0.0, 4.0), Direction::Z, Direction::X), 8.0);
        let ring = crate::make::edge_from_curve(&mut model, Curve::Circle(circle)).unwrap();
        let r = normal(&mut model, &[ring], cylinder, &NormalProjection::default()).unwrap();
        assert_eq!(model.unique_subshapes(r, NodeKind::Edge).len(), 1);
        assert_eq!(model.unique_subshapes(r, NodeKind::Vertex).len(), 1);
        let length = crate::properties::length(&model, r).mass;
        assert!((length - std::f64::consts::TAU * 5.0).abs() < 1e-3, "{length}");
        let near = NormalProjection { max_distance: Some(2.0), ..Default::default() };
        let r = normal(&mut model, &[ring], cylinder, &near).unwrap();
        assert!(model.unique_subshapes(r, NodeKind::Edge).is_empty());
    }
}

/// An edge's projection along a direction: the edge made (run as the edge projected runs),
/// and its pcurve on the surface projected onto over the edge's range.
struct AlongEdge {
    edge: Shape,
    pcurve: crate::geometry::curve2::BSplineCurve2,
    first: f64,
    last: f64,
}

/// The point of `surface` met from `p` along `direction` (moving either way), from `start`
/// (its parameters and how far along): Newton's method on `S(u, v) - (p + s d)`.
fn meet_along(surface: &Surface, p: Point, direction: crate::math::Vector, start: (f64, f64, f64)) -> Option<(f64, f64, f64)> {
    let (mut u, mut v, mut s) = start;
    let scale = 1.0 + (p - Point::ORIGIN).length();
    for _ in 0..50 {
        let (q, du, dv) = surface.d1(u, v);
        let f = q - (p + direction * s);
        if f.length() <= 1e-12 * scale {
            return Some((u, v, s));
        }
        // [du dv -d] (Δu, Δv, Δs) = -f
        let columns = [[du.x, du.y, du.z], [dv.x, dv.y, dv.z], [-direction.x, -direction.y, -direction.z]];
        let step = solve3(columns, [f.x, f.y, f.z])?;
        u -= step[0];
        v -= step[1];
        s -= step[2];
    }
    let (q, _, _) = surface.d1(u, v);
    ((q - (p + direction * s)).length() <= 1e-8 * scale).then_some((u, v, s))
}

/// Edge `edge`'s projection onto `surface` along `direction`, between the vertices `ends`
/// (its own start and end projected), the nearest projection where there are several.
fn edge_along(model: &mut Model, edge: Shape, surface: &Surface, direction: crate::math::Vector, ends: [crate::topology::VertexId; 2], tolerance: f64) -> Result<AlongEdge> {
    use crate::geometry::curve2::BSplineCurve2;
    let forward = Shape { orientation: crate::topology::Orientation::Forward, ..edge };
    let (curve, first, last) = model.edge_curve(forward).ok_or_else(|| Error::InvalidInput("a degenerate edge projects to nothing".into()))?;
    let (u_period, v_period) = (surface.u_period(), surface.v_period());
    let solve = |n: usize| -> Result<(Vec<f64>, Vec<Point>, Vec<Point2>)> {
        let ts: Vec<f64> = (0..=n).map(|k| first + (last - first) * k as f64 / n as f64).collect();
        let mut points = Vec::with_capacity(n + 1);
        let mut uvs: Vec<Point2> = Vec::with_capacity(n + 1);
        let mut previous: Option<(f64, f64, f64)> = None;
        for &t in &ts {
            let p = curve.point(t);
            let start = match previous {
                Some(s) => s,
                None => {
                    let (u, v, _) = nearest_on_surface(p, surface, None);
                    (u, v, (surface.point(u, v) - p).dot(direction))
                }
            };
            let (u, v, s) = meet_along(surface, p, direction, start).ok_or_else(|| Error::Failed("an edge does not project onto the surface along the direction".into()))?;
            previous = Some((u, v, s));
            points.push(p + direction * s);
            // (Continuous across the surface's periods.)
            let (mut u, mut v) = (u, v);
            if let Some(last) = uvs.last() {
                if let Some(pu) = u_period {
                    u -= ((u - last.x) / pu).round() * pu;
                }
                if let Some(pv) = v_period {
                    v -= ((v - last.y) / pv).round() * pv;
                }
            }
            uvs.push(Point2::new(u, v));
        }
        Ok((ts, points, uvs))
    };
    // More samples until the fitted curve holds to fresh projections halfway between.
    let mut n = 16;
    let (curve3, pcurve) = loop {
        let (ts, points, uvs) = solve(n)?;
        let curve3 = fit::interpolate_with_parameters(&points, &ts)?;
        let pcurve: BSplineCurve2 = fit::interpolate_with_parameters(&uvs, &ts)?;
        let (fine_ts, fine_points, _) = solve(2 * n)?;
        let worst = fine_ts.iter().zip(&fine_points).map(|(&t, &p)| curve3.point(t).distance(p)).fold(0.0, f64::max);
        if worst <= tolerance || n >= 512 {
            break (curve3, pcurve);
        }
        n *= 2;
    };
    let e = model.add_edge(ends[0], ends[1], tolerance)?;
    model.set_edge_curve(e.into(), Curve::BSpline(curve3), tolerance)?;
    model.set_edge_range(e.into(), first, last)?;
    Ok(AlongEdge { edge: Shape { orientation: edge.orientation, ..Shape::from(e) }, pcurve, first, last })
}

/// The edges of `wire` projected onto `surface` along `direction`, in order along it, with
/// their vertices shared (each of `vertices` projected once).
fn wire_along_surface(
    model: &mut Model,
    wire: Shape,
    surface: &Surface,
    direction: crate::math::Vector,
    vertices: &mut std::collections::HashMap<crate::topology::NodeId, crate::topology::VertexId>,
    tolerance: f64,
) -> Result<Vec<AlongEdge>> {
    let mut out = Vec::new();
    for edge in model.ordered_edges(wire) {
        let forward = Shape { orientation: crate::topology::Orientation::Forward, ..edge };
        let (a, b) = model.edge_vertices(forward);
        let mut ends = [crate::topology::VertexId::INVALID; 2];
        for (k, v) in [a, b].into_iter().enumerate() {
            ends[k] = match vertices.get(&v.node) {
                Some(&w) => w,
                None => {
                    let p = model.point(v);
                    let (u, w, _) = nearest_on_surface(p, surface, None);
                    let (u, w, s) = meet_along(surface, p, direction, (u, w, (surface.point(u, w) - p).dot(direction)))
                        .ok_or_else(|| Error::Failed("a vertex does not project onto the surface along the direction".into()))?;
                    let id = model.add_vertex({ let (a, b) = (surface.point(u, w), p + direction * s); a + (b - a) * 0.5 }, tolerance);
                    vertices.insert(v.node, id);
                    id
                }
            };
        }
        out.push(edge_along(model, edge, surface, direction, ends, tolerance)?);
    }
    Ok(out)
}

/// `wire` projected onto the surface of the face `target` along `direction` (the projection
/// nearest the wire): a wire.
#[doc(alias = "BRepProj_Projection")]
pub fn wire_along(model: &mut Model, wire: Shape, target: Shape, direction: crate::math::Vector) -> Result<Shape> {
    let surface = target_surface(model, target)?;
    let direction = unit(direction)?;
    let mut vertices = std::collections::HashMap::new();
    let edges = wire_along_surface(model, wire, &surface, direction, &mut vertices, 1e-6)?;
    let coedges: Vec<CoEdgeId> = edges.iter().map(|e| model.add_coedge(crate::topology::EdgeId::try_from(e.edge.node).expect("an edge"), e.edge.orientation)).collect::<Result<_>>()?;
    Ok(Shape::from(model.add_wire(&coedges)?))
}

/// `face` projected onto the surface of the face `target` along `direction`: a face on that
/// surface bounded by the projections of `face`'s wires (CadQuery's `Face.project`).
pub fn face_along(model: &mut Model, face: Shape, target: Shape, direction: crate::math::Vector) -> Result<Shape> {
    let surface = target_surface(model, target)?;
    let direction = unit(direction)?;
    let tolerance = 1e-6;
    let wires: Vec<Shape> = model.children(face).collect();
    let mut vertices = std::collections::HashMap::new();
    let mut loops = Vec::with_capacity(wires.len());
    for &w in &wires {
        loops.push(wire_along_surface(model, w, &surface, direction, &mut vertices, tolerance)?);
    }
    let new = model.new_empty(NodeKind::Face);
    model.set_face_surface(new, surface.clone(), tolerance)?;
    let (u_period, v_period) = (surface.u_period(), surface.v_period());
    for (k, edges) in loops.iter_mut().enumerate() {
        // Each pcurve in the periods where the one before it ends.
        let mut previous_end: Option<Point2> = None;
        let mut uv = Vec::new();
        for e in edges.iter_mut() {
            let reversed = e.edge.orientation == crate::topology::Orientation::Reversed;
            let (start, end) = if reversed { (e.last, e.first) } else { (e.first, e.last) };
            if let Some(target) = previous_end {
                let at = e.pcurve.point(start);
                let shift = Point2::new(
                    u_period.map_or(0.0, |p| ((target.x - at.x) / p).round() * p),
                    v_period.map_or(0.0, |p| ((target.y - at.y) / p).round() * p),
                );
                if shift.x != 0.0 || shift.y != 0.0 {
                    let poles = e.pcurve.poles().iter().map(|q| Point2::new(q.x + shift.x, q.y + shift.y)).collect();
                    e.pcurve = crate::geometry::curve2::BSplineCurve2::new(poles, e.pcurve.knots().to_vec(), e.pcurve.multiplicities().to_vec(), e.pcurve.degree());
                }
            }
            previous_end = Some(e.pcurve.point(end));
            uv.extend((0..16).map(|i| e.pcurve.point(start + (end - start) * i as f64 / 16.0)));
        }
        // The outer loop counter-clockwise in the surface's parameters, holes clockwise.
        let n = uv.len();
        let area: f64 = (0..n).map(|i| uv[i].x * uv[(i + 1) % n].y - uv[(i + 1) % n].x * uv[i].y).sum();
        let flip = (k == 0) == (area < 0.0);
        let mut coedges = Vec::with_capacity(edges.len());
        let order: Vec<&AlongEdge> = if flip { edges.iter().rev().collect() } else { edges.iter().collect() };
        for e in order {
            let o = if flip { e.edge.orientation.reversed() } else { e.edge.orientation };
            let c = model.add_coedge(crate::topology::EdgeId::try_from(e.edge.node).expect("an edge"), o)?;
            model.set_pcurve_with_range(c, crate::geometry::curve2::Curve2::BSpline(e.pcurve.clone()), e.first, e.last)?;
            coedges.push(c);
        }
        let w = model.add_wire(&coedges)?;
        model.add(new, Shape::from(w))?;
    }
    Ok(new)
}

/// The surface of the face `target`, placed as it is.
fn target_surface(model: &Model, target: Shape) -> Result<Surface> {
    let face = model.unique_subshapes(target, NodeKind::Face).into_iter().next().ok_or_else(|| Error::InvalidInput("nothing to project onto: no face".into()))?;
    Ok(model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.transformed(&face.placement))
}

/// `v` of unit length (an error if it has none).
fn unit(v: crate::math::Vector) -> Result<crate::math::Vector> {
    if v.length() > 1e-300 { Ok(v.normalized()) } else { Err(Error::InvalidInput("a projection along a zero direction".into())) }
}

/// The solution of the 3×3 system with these columns (Cramer's rule), if not singular.
fn solve3(c: [[f64; 3]; 3], r: [f64; 3]) -> Option<[f64; 3]> {
    let det = |a: [[f64; 3]; 3]| a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2]) + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]);
    let d = det(c);
    if d.abs() <= 1e-300 {
        return None;
    }
    let mut out = [0.0; 3];
    for (i, o) in out.iter_mut().enumerate() {
        let mut m = c;
        m[i] = r;
        *o = det(m) / d;
    }
    Some(out)
}
