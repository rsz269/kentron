//! Fillets whose ball rolls off the faces at the edge onto the faces beyond them (where a
//! face is narrower than the fillet, along the edge or at its ends), as OCCT's do.
//!
//! In each section square to the edge the ball's centre lies at its radius from both sides
//! of the edge, each side a face at the edge with its neighbours across edges away from
//! the edge's ends. Distances are to the faces as bounded (to a face's edges where its
//! surface's nearest point lies off it), so a ball touching a face beyond, or rolling over
//! the edge between two, comes out of the same equations.
//!
//! The solid is then rebuilt without tangent meetings for a boolean to find:
//!
//! 1. A wedge is cut away: in each section, from one contact out and round beyond the edge
//!    to the other, closed by the chord between the contacts. Its faces cross the solid's
//!    faces, along the contact curves and away from them.
//! 2. The lens between that chord and the ball's arc is lofted the same way and kept
//!    within the solid (trimmed by the faces at the edge's ends, and where the ball hangs
//!    over a hole, by the faces there).
//! 3. The chord faces of both are dropped, and the rest sewn along the contact curves and
//!    the chords' traces on the end faces. Faces split along the way are merged again.
//!
//! At a concave edge the fillet adds what lies between the faces and the ball. The wedge
//! there runs from within the ball out through each contact, along a line crossing the face
//! (so that cutting it away splits the faces exactly along the contacts); the solid's faces
//! it leaves and the ball's arcs between the contacts are sewn. Around a closed edge the
//! sections make rings, with no ends.

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::extrema::{nearest_on_curve, nearest_on_surface};
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Frame, Point, Point2, Transform, Transformable, Vector};
use crate::topology::{CoEdgeId, Model, NodeKind, Orientation, Shape};
use crate::{boolean, loft, make, sew, unify};

mod trim;

/// A face measured from: its surface where it is, its boundary in its parameters (to tell
/// points on it) and its edges' curves.
struct Bounded {
    face: Shape,
    /// A box around it.
    bbox: crate::math::BoundingBox,
    surface: Surface,
    loops: Vec<Vec<Point2>>,
    /// The pcurves `loops` samples, each with its range.
    loop_curves: Vec<(crate::geometry::curve2::Curve2, f64, f64)>,
    edges: Vec<(Curve, f64, f64)>,
    /// The edges themselves.
    edge_shapes: Vec<Shape>,
    /// The faces across each of `edges`.
    beyond: Vec<Vec<Surface>>,
    u_period: Option<f64>,
    v_period: Option<f64>,
    uv_tolerance: f64,
    near: std::cell::Cell<Option<(f64, f64)>>,
}

impl Bounded {
    fn of(model: &Model, face: Shape) -> Option<Bounded> {
        Self::of_with(model, face, true)
    }

    /// [`Bounded::of`], with the faces beyond its edges and its box (for measuring distances
    /// to it among others) found or not: only telling points on it from points off it, they
    /// are not needed, and finding them is the slow part.
    fn of_with(model: &Model, face: Shape, whole: bool) -> Option<Bounded> {
        let surface = model.surface(face)?;
        let surface = if face.placement == Transform::IDENTITY { surface.clone() } else { surface.transformed(&face.placement) };
        let mut loops = Vec::new();
        let mut loop_curves = Vec::new();
        for coedge in model.unique_subshapes(face, NodeKind::CoEdge) {
            let Ok(id) = CoEdgeId::try_from(coedge.node) else { continue };
            let Some((pcurve, first, last)) = model.effective_pcurve(id) else { continue };
            const N: usize = 128;
            loops.push((0..=N).map(|i| pcurve.point(first + (last - first) * i as f64 / N as f64)).collect::<Vec<Point2>>());
            loop_curves.push((pcurve, first, last));
        }
        let mut edges = Vec::new();
        let mut edge_shapes = Vec::new();
        for e in model.unique_subshapes(face, NodeKind::Edge) {
            let e = Shape { orientation: Orientation::Forward, ..e };
            let Some(curve) = model.edge_curve(e) else { continue };
            edges.push(curve);
            edge_shapes.push(e);
        }
        let ids: Vec<crate::topology::EdgeId> = edge_shapes.iter().filter_map(|e| crate::topology::EdgeId::try_from(e.node).ok()).collect();
        let beyond: Vec<Vec<Surface>> = if !whole {
            vec![Vec::new(); ids.len()]
        } else {
            model
            .faces_around(&ids)
            .into_iter()
            .map(|faces| faces.into_iter().filter(|&f| crate::topology::NodeId::from(f) != face.node).filter_map(|f| model.surface(Shape::from(f)).cloned()).collect())
            .collect()
        };
        let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
        for q in loops.iter().flatten() {
            (lo.x, lo.y, hi.x, hi.y) = (lo.x.min(q.x), lo.y.min(q.y), hi.x.max(q.x), hi.y.max(q.y));
        }
        let uv_tolerance = 1e-7 * lo.distance(hi);
        let bbox = if whole { crate::properties::bounding_box(model, face) } else { crate::math::BoundingBox::default() };
        Some(Bounded { face, bbox, edge_shapes, u_period: surface.u_period(), v_period: surface.v_period(), surface, loops, loop_curves, edges, beyond, uv_tolerance, near: std::cell::Cell::new(None) })
    }

    /// Whether `(u, v)` lies on the face (inside its loops, even-odd, or on them).
    fn contains(&self, uv: (f64, f64)) -> bool {
        self.locate(uv) != Some(false)
    }

    /// Where `(u, v)` lies: inside the face (`Some(true)`), on its loops (`None`) or outside
    /// (`Some(false)`).
    fn locate(&self, (u, v): (f64, f64)) -> Option<bool> {
        // How far `q` lies from the segment `a b`.
        let gap = |a: Point2, b: Point2, q: Point2| {
            let (ab, aq) = (b - a, q - a);
            let s = (aq.dot(ab) / ab.dot(ab).max(1e-300)).clamp(0.0, 1.0);
            (aq - ab * s).length()
        };
        let shifts = |period: Option<f64>| match period {
            Some(p) => vec![0.0, -p, p],
            None => vec![0.0],
        };
        for su in shifts(self.u_period) {
            for sv in shifts(self.v_period) {
                let q = Point2::new(u + su, v + sv);
                // (The ray a little off `q`, by less than the loops are told from: through a
                // point where a loop closes, its two ends a rounding error apart, as a circle's
                // start level with `q`, it would cross it once too few.)
                let y = q.y + 0.37 * self.uv_tolerance;
                let crosses = |a: Point2, b: Point2| (a.y > y) != (b.y > y) && a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x) > q.x;
                let mut crossings = 0;
                for (l, (pcurve, first, last)) in self.loops.iter().zip(&self.loop_curves) {
                    let step = (last - first) / (l.len() - 1) as f64;
                    for (i, w) in l.windows(2).enumerate() {
                        // (Near a curve its samples stray from, its own points: within their
                        // sag a point is told from the loop by the loop itself.)
                        if gap(w[0], w[1], q) <= 0.25 * w[0].distance(w[1]) + self.uv_tolerance {
                            const FINE: usize = 128;
                            let at = |j: usize| if j == 0 { w[0] } else if j == FINE { w[1] } else { pcurve.point(first + step * (i as f64 + j as f64 / FINE as f64)) };
                            let mut a = w[0];
                            for j in 1..=FINE {
                                let b = at(j);
                                if gap(a, b, q) <= self.uv_tolerance {
                                    return None;
                                }
                                crossings += usize::from(crosses(a, b));
                                a = b;
                            }
                        } else {
                            crossings += usize::from(crosses(w[0], w[1]));
                        }
                    }
                }
                if crossings % 2 == 1 {
                    return Some(true);
                }
            }
        }
        Some(false)
    }

    /// The distance from `p` to the face, the nearest point of it, and where that lies: on
    /// the face (`None`) or on its edge of that index.
    fn nearest(&self, p: Point) -> (f64, Point, Option<usize>) {
        let (u, v, d) = nearest_on_surface(p, &self.surface, self.near.get());
        // (A foot on the surface's bounds and on the face's loops, the point lying beyond
        // them: the nearest point is on one of the face's edges, found there with which.)
        let on_bounds = {
            use crate::geometry::ParametricSurface;
            let (u0, u1, v0, v1) = self.surface.bounds();
            (self.u_period.is_none() && (u <= u0 || u >= u1)) || (self.v_period.is_none() && (v <= v0 || v >= v1))
        };
        if u.is_finite() && v.is_finite() && self.locate((u, v)).unwrap_or(!on_bounds) {
            self.near.set(Some((u, v)));
            return (d, self.surface.point(u, v), None);
        }
        self.edges
            .iter()
            .enumerate()
            .map(|(i, (c, a, b))| {
                let foot = nearest_on_curve(p, c, *a, *b);
                let q = c.point(foot.parameter);
                (q.distance(p), q, Some(i))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .unwrap_or((f64::INFINITY, p, None))
    }
}

/// One side of the edge: the faces to keep the ball's radius from.
struct Side(Vec<Bounded>);

/// What a contact lies on: a face of the side, and an edge of it or the face itself.
type Support = (usize, Option<usize>);

impl Side {
    /// The distance from `p` to the side, the nearest point and what that lies on: of the
    /// points nearest on its faces, those toward `at` (the ball's side facing the edge, the
    /// part of it the fillet is: a face its far side reaches is no matter).
    fn nearest(&self, p: Point, at: Point) -> (f64, Point, Support) {
        // (Nearest boxes first; a face whose box lies further than a point found is passed
        // over.)
        let mut order: Vec<(f64, usize)> = self.0.iter().enumerate().map(|(i, f)| (f.bbox.distance_to_point(p), i)).collect();
        order.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut best = (f64::INFINITY, p, (0, None));
        for (reach, i) in order {
            if reach > best.0 {
                break;
            }
            let (d, q, edge) = self.0[i].nearest(p);
            // (Of two as near, the first: an edge two faces share is as near on both.)
            let nearer = d < best.0 || (d == best.0 && i < best.2 .0);
            if nearer && (q - p).dot(at - p) > 0.0 {
                best = (d, q, (i, edge));
            }
        }
        best
    }
}

/// The faces of each side of `edge` (between `faces`): each face and its neighbours across
/// edges meeting neither end of `edge` (the faces beyond it, not those at the edge's ends).
fn sides(model: &Model, shape: Shape, edge: Shape, faces: [Shape; 2]) -> [Side; 2] {
    let (a, b) = model.edge_vertices(edge);
    let ends = [a.node, b.node];
    let all = model.unique_subshapes(shape, NodeKind::Face);
    let side = |face: Shape, other: Shape| -> Side {
        let mut members = vec![face];
        for e in model.unique_subshapes(face, NodeKind::Edge) {
            if e.node == edge.node {
                continue;
            }
            let (p, q) = model.edge_vertices(e);
            if ends.contains(&p.node) || ends.contains(&q.node) {
                continue;
            }
            for &f in &all {
                let fresh = f.node != face.node && f.node != other.node && !members.iter().any(|m| m.node == f.node);
                if fresh && model.unique_subshapes(f, NodeKind::Edge).iter().any(|x| x.node == e.node) {
                    members.push(f);
                }
            }
        }
        Side(members.into_iter().filter_map(|f| Bounded::of(model, f)).collect())
    };
    [side(faces[0], faces[1]), side(faces[1], faces[0])]
}

/// A section of the fillet: the point of the edge (or its continuation), the ball's centre
/// and where it touches the two sides.
#[derive(Clone, Copy)]
struct Section {
    /// The edge's parameter.
    t: f64,
    at: Point,
    center: Point,
    contacts: [Point; 2],
    /// The contacts' parameters on the faces they lie on (not on edges).
    uv: [Option<Point2>; 2],
    /// The ball's radius there.
    radius: f64,
}

/// Where [`roll`] puts a section at each switch of what the contacts lie on: nowhere, or
/// one found by bisection (where the contact first lies on what it switches to), or found
/// faster where the ball rolls from a face onto its edge or back.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Switches {
    Ignored,
    Bisected,
    Found,
}

/// The ball along an edge, section by section (see [`roll`]).
struct Rolled {
    sections: Vec<Section>,
    /// What the contacts of each section lie on.
    supports: Vec<[Support; 2]>,
    /// The sections where the contacts switch.
    breaks: Vec<usize>,
    side_set: [Side; 2],
    convex: bool,
    closed: bool,
    /// The edge's tangents at its start and end.
    tangents: [Vector; 2],
}

/// The curve a fillet runs along: one edge's, or a chain of edges going on into one another
/// smoothly. Piece `i` runs from `a` to `b` of its curve (the chain's way along it, which
/// may be against the curve's), at spine parameters `i` to `i + 1`.
struct Spine {
    pieces: Vec<(Curve, f64, f64)>,
    /// Each piece's length up to each of [`Spine::SAMPLES`] steps along it, from the
    /// spine's start.
    lengths: Vec<Vec<f64>>,
    closed: bool,
}

impl Spine {
    const SAMPLES: usize = 64;
    /// The fewest sections along one piece.
    const LEAST: usize = 4;

    fn new(pieces: Vec<(Curve, f64, f64)>, closed: bool) -> Spine {
        let mut total = 0.0;
        let lengths = pieces
            .iter()
            .map(|(c, a, b)| {
                let points: Vec<Point> = (0..=Self::SAMPLES).map(|i| c.point(a + (b - a) * i as f64 / Self::SAMPLES as f64)).collect();
                let start = total;
                let along: Vec<f64> = std::iter::once(start).chain(points.windows(2).scan(start, |l, w| { *l += w[0].distance(w[1]); Some(*l) })).collect();
                total = *along.last().expect("samples");
                along
            })
            .collect();
        Spine { pieces, lengths, closed }
    }

    /// The parameter range: 0 to the number of pieces.
    fn end(&self) -> f64 {
        self.pieces.len() as f64
    }

    /// The piece at spine parameter `s`, and the fraction along it.
    fn piece(&self, s: f64) -> (usize, f64) {
        let i = (s.floor().max(0.0) as usize).min(self.pieces.len() - 1);
        (i, (s - i as f64).clamp(0.0, 1.0))
    }

    /// The point at `s` and the unit direction along the spine there.
    fn d1(&self, s: f64) -> (Point, Vector) {
        let (i, f) = self.piece(s);
        let (c, a, b) = &self.pieces[i];
        let (p, d) = c.d1(a + (b - a) * f);
        (p, d.normalized() * (b - a).signum())
    }

    /// About `count` parameters along the spine, from its start to its end: each piece's
    /// share by its length, at least [`Spine::LEAST`] (evenly in its own parameter within
    /// it). Sections spaced so are about as far apart on long pieces as on short ones, and
    /// a short piece is not passed over.
    fn parameters(&self, count: usize) -> Vec<f64> {
        let total = *self.lengths.last().and_then(|l| l.last()).expect("samples");
        let mut ts = Vec::with_capacity(count + 1);
        for (i, along) in self.lengths.iter().enumerate() {
            let length = along.last().expect("samples") - along[0];
            let n = ((count as f64 * length / total.max(1e-300)).round() as usize).max(Self::LEAST);
            ts.extend((0..n).map(|k| i as f64 + k as f64 / n as f64));
        }
        ts.push(self.end());
        ts
    }

    /// The fraction of the spine's length up to `s`.
    fn fraction(&self, s: f64) -> f64 {
        let (i, f) = self.piece(s);
        let along = &self.lengths[i];
        let x = f * Self::SAMPLES as f64;
        let k = (x.floor() as usize).min(Self::SAMPLES - 1);
        let l = along[k] + (along[k + 1] - along[k]) * (x - k as f64);
        let total = *self.lengths.last().and_then(|l| l.last()).expect("samples");
        l / total.max(1e-300)
    }
}

/// The faces of `shape` bounded by `edge`, each with the direction of the edge in its
/// boundary (+1 along the edge's own).
fn faces_of(model: &Model, shape: Shape, edge: Shape) -> Vec<(Shape, f64)> {
    let mut faces: Vec<(Shape, f64)> = Vec::new();
    for face in model.explore(shape, NodeKind::Face) {
        for wire in model.children(face) {
            for coedge in model.children(wire) {
                let Ok(c) = CoEdgeId::try_from(coedge.node) else { continue };
                if crate::topology::NodeId::from(model[c].edge) == edge.node {
                    faces.push((face, if model[c].orientation.compose(wire.orientation) == Orientation::Reversed { -1.0 } else { 1.0 }));
                }
            }
        }
    }
    faces
}

/// The ball of radius `radius(s)` along `edge` of `shape` (`s` the fraction of the edge's
/// length from its start), in `count` sections and one more where the contacts change what
/// they lie on (the ball rolling onto another face, or over an edge; found by bisection,
/// where the fillet's faces end). Around a closed edge the last section is the first.
fn roll(model: &Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64, count: usize, switches: Switches) -> Result<Rolled> {
    let faces = faces_of(model, shape, edge);
    let [(f0, along0), (f1, along1)] = faces[..] else {
        return Err(Error::InvalidInput(format!("a fillet edge must bound two faces, it bounds {}", faces.len())));
    };
    let (curve, t0, t1) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
    let (va, vb) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
    let side_set = sides(model, shape, edge, [f0, f1]);
    let spine = Spine::new(vec![(curve, t0, t1)], va.node == vb.node);
    roll_on(model, &spine, side_set, [(f0, along0), (f1, along1)], radius, count, switches)
}

/// [`roll`] along `spine`, between the sides `side_set`; `faces` the faces at the spine's
/// start, each with the spine's direction in its boundary (+1 along it); `count` sections
/// for each piece of the spine.
fn roll_on(model: &Model, spine: &Spine, side_set: [Side; 2], faces: [(Shape, f64); 2], radius: &dyn Fn(f64) -> f64, count: usize, switches: Switches) -> Result<Rolled> {
    let [(f0, along0), (f1, along1)] = faces;
    let (t0, t1) = (0.0, spine.end());
    let closed = spine.closed;
    let count = count * spine.pieces.len();
    // The radius at a parameter of the spine (by the fraction of its length up to there),
    // and the largest, for scales.
    let radius_at = |t: f64| radius(spine.fraction(t));
    let largest = largest(radius);
    if !(largest > 0.0) {
        return Err(Error::InvalidInput(format!("fillet radius must be positive, got {largest}")));
    }
    // The faces' outward normals at the start of the edge, for the ball's first place.
    let outward = |model: &Model, face: Shape, p: Point| -> Option<Vector> {
        let surface = model.surface(face)?.transformed(&face.placement);
        let (u, v, _) = nearest_on_surface(p, &surface, None);
        let (_, du, dv) = surface.d1(u, v);
        let n = du.cross(dv);
        // (None at a singular point, as a cone's apex.)
        if n.length() <= 1e-300 {
            return None;
        }
        let n = n.normalized();
        Some(if face.orientation == Orientation::Reversed { -n } else { n })
    };
    let path = spine;
    let spine = |t: f64| -> (Point, Vector) { path.d1(t) };
    // The first place: between the faces' tangent planes at the start of the edge (or just
    // along it, where a face's surface is singular there).
    let ((p0, d0), (n0, n1)) = [0.0, 1e-4, 1e-2]
        .iter()
        .find_map(|&f| {
            let (p, d) = spine(t0 + (t1 - t0) * f);
            Some(((p, d), (outward(model, f0, p)?, outward(model, f1, p)?)))
        })
        .ok_or_else(|| Error::InvalidInput("a fillet's faces have no normals at its edge".into()))?;
    // (Into each face, square to the edge.)
    let (w0, w1) = (n0.cross(d0 * along0).normalized(), n1.cross(d0 * along1).normalized());
    // (Convex where the faces' normals point apart from where the material lies between;
    // the ball lies between the faces either way, in the material or out of it.)
    let convex = n0.dot(w1) < 0.0;
    let half = 0.5 * w0.dot(w1).clamp(-1.0, 1.0).acos();
    let mut offset = (w0 + w1).normalized() * (radius_at(t0) / half.sin().max(1e-3));
    // In the section through a point of the edge, Newton steps from a nearby centre (as an
    // offset from the edge): the section, and what the contacts lie on.
    let solve = |t: f64, offset: Vector| -> Result<(Section, [Support; 2])> {
        let (at, tangent) = spine(t);
        let radius = radius_at(t);
        let ex = (offset - tangent * offset.dot(tangent)).normalized();
        let ey = tangent.cross(ex);
        // The distances to the sides less the radius, the nearest points and what they lie
        // on; the distances' gradients are the directions from those points.
        let f = |x: f64, y: f64| {
            let p = at + ex * x + ey * y;
            let (a, b) = (side_set[0].nearest(p, at), side_set[1].nearest(p, at));
            (p, [a.0 - radius, b.0 - radius], [a.1, b.1], [a.2, b.2])
        };
        let (mut x, mut y) = (offset.dot(ex), offset.dot(ey));
        for _ in 0..80 {
            let (p, v, feet, on) = f(x, y);
            if v[0].abs().max(v[1].abs()) <= 1e-10 * radius.max(1.0) {
                let uv = |k: usize, (face, rim): Support| if rim.is_none() { side_set[k].0[face].near.get().map(|(u, v)| Point2::new(u, v)) } else { None };
                return Ok((Section { t, at, center: p, contacts: feet, uv: [uv(0, on[0]), uv(1, on[1])], radius }, on));
            }
            // (Nothing on a side facing the edge, or the centre on it: no section here.)
            let (d0, d1) = (p - feet[0], p - feet[1]);
            if !(v[0].is_finite() && v[1].is_finite()) || d0.length() <= 1e-12 * radius || d1.length() <= 1e-12 * radius {
                break;
            }
            let n = [d0.normalized(), d1.normalized()];
            let j = [[n[0].dot(ex), n[0].dot(ey)], [n[1].dot(ex), n[1].dot(ey)]];
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            if !(det.abs() > 1e-300) {
                break;
            }
            let (dx, dy) = ((v[0] * j[1][1] - v[1] * j[0][1]) / det, (v[1] * j[0][0] - v[0] * j[1][0]) / det);
            let scale = (dx.hypot(dy) / (0.25 * radius)).max(1.0);
            x -= dx / scale;
            y -= dy / scale;
        }
        Err(Error::InvalidInput("no fillet section here (the ball rolling over)".into()))
    };
    // Where one side switches between a face and an edge of it: from the edge to its
    // surface's point nearest the ball's centre (into the face before the switch, out of it
    // after), and whether the face is the first support.
    // (The side switching, the face and its edge, and whether the face is the first support.)
    let switching = |previous: [Support; 2], support: [Support; 2]| -> Option<(usize, usize, usize, bool)> {
        let changed: Vec<usize> = (0..2).filter(|&k| previous[k] != support[k]).collect();
        let [k] = changed[..] else { return None };
        match (previous[k], support[k]) {
            ((f, None), (g, Some(e))) if f == g => Some((k, f, e, true)),
            ((g, Some(e)), (f, None)) if f == g => Some((k, f, e, false)),
            _ => None,
        }
    };
    let onto_edge = |s: &Section, previous: [Support; 2], support: [Support; 2]| -> Option<(Vector, bool)> {
        let (k, face, rim, first) = switching(previous, support)?;
        let bounded = &side_set[k].0[face];
        let (u, v, _) = nearest_on_surface(s.center, &bounded.surface, bounded.near.get());
        let q = bounded.surface.point(u, v);
        let (curve, a, b) = &bounded.edges[rim];
        Some((q - nearest_on_curve(q, curve, *a, *b).point, first))
    };
    let mut sections: Vec<Section> = Vec::with_capacity(count + 1);
    let mut supports: Vec<[Support; 2]> = Vec::new();
    let mut breaks: Vec<usize> = Vec::new();
    let mut last_t = t0;
    let ts: Vec<f64> = if path.pieces.len() == 1 { (0..=count).map(|i| t0 + (t1 - t0) * i as f64 / count as f64).collect() } else { path.parameters(count) };
    for t in ts {
        let (section, support) = solve(t, offset)?;
        if let (Some(&previous), Some(&before)) = (supports.last(), sections.last()) {
            let (mut previous, mut before) = (previous, before);
            // (Each switch within the step in turn: the two sides may switch close together.)
            for _ in 0..4 {
                if previous == support || switches == Switches::Ignored {
                    break;
                }
                let (mut lo, mut hi, mut start) = (last_t, t, before.center - before.at);
                let (v_lo, v_hi) = (onto_edge(&before, previous, support), onto_edge(&section, previous, support));
                let (s, reached) = match (v_lo, v_hi) {
                    // (A section already at the switch, as a symmetrical solid's may be: no way
                    // into the face to measure along, and that section is the switch.)
                    (Some((v_lo, first)), Some((v_hi, _))) if switches == Switches::Found && (if first { v_lo } else { v_hi }).length() <= 1e-9 * (1.0 + largest) => {
                        (if first { before } else { section }, support)
                    }
                    // Where the ball rolls from a face onto its edge (or back): how far into the
                    // face the surface's point nearest the centre lies, through nothing at the
                    // switch (along the way into the face where it lies on it), by secants
                    // through points on the face, where it runs smoothly (past the edge it runs
                    // on otherwise).
                    (Some((v_lo, first)), Some((v_hi, _))) if switches == Switches::Found => {
                        let into = if first { v_lo } else { v_hi };
                        let into = if into.length() > 0.0 { into.normalized() } else { into };
                        let measure = |s: &Section| onto_edge(s, previous, support).map_or(0.0, |(v, _)| v.dot(into));
                        let (g_lo, g_hi) = (v_lo.dot(into), v_hi.dot(into));
                        let (mut inside, mut past) = if first { ((lo, g_lo), (hi, g_hi)) } else { ((hi, g_hi), (lo, g_lo)) };
                        let mut before_inside: Option<(f64, f64)> = None;
                        let mut found = None;
                        for _ in 0..60 {
                            let (a, ga) = inside;
                            let (b, gb) = past;
                            let secant = before_inside.filter(|&(_, gp)| gp != ga).map(|(p, gp)| a - ga * (a - p) / (ga - gp));
                            let falsi = a - ga * (b - a) / (gb - ga);
                            let between = |x: f64| (x - a) * (x - b) < 0.0;
                            let mid = secant.filter(|&x| between(x)).or(Some(falsi).filter(|&x| between(x))).unwrap_or(0.5 * (a + b));
                            let (s, _) = solve(mid, start)?;
                            let g = measure(&s);
                            if g.abs() <= 1e-9 * (1.0 + largest) || (b - a).abs() <= 1e-9 * (t1 - t0).abs() {
                                found = Some(s);
                                break;
                            }
                            if g > 0.0 {
                                before_inside = Some(inside);
                                inside = (mid, g);
                            } else {
                                past = (mid, g);
                            }
                            start = s.center - s.at;
                        }
                        match found {
                            Some(s) => (s, support),
                            None => (solve(past.0, start)?.0, support),
                        }
                    }
                    _ => {
                        // Otherwise (and for booleans, which want the switch where the
                        // contact first lies on the edge) by bisection.
                        while (hi - lo).abs() > 1e-12 * (t1 - t0).abs() {
                            let mid = 0.5 * (lo + hi);
                            let (s, k) = solve(mid, start)?;
                            if k == previous {
                                lo = mid;
                                start = s.center - s.at;
                            } else {
                                hi = mid;
                            }
                        }
                        solve(hi, start)?
                    }
                };
                // (A switch found a hair from the section before or this one (as where the
                // contacts switch just where two of a chain's edges meet): that section is the
                // break, not one beside it, which would make the faces slivers.)
                let near = 1e-6 * (t1 - t0).abs();
                if (s.t - last_t).abs() <= near {
                    let k = sections.len() - 1;
                    if breaks.last() != Some(&k) {
                        breaks.push(k);
                    }
                    supports[k] = reached;
                } else if (t - s.t).abs() <= near {
                    breaks.push(sections.len());
                    break;
                } else {
                    breaks.push(sections.len());
                    sections.push(s);
                    supports.push(reached);
                }
                (previous, before, last_t) = (reached, s, s.t);
            }
        }
        offset = section.center - section.at;
        sections.push(section);
        supports.push(support);
        last_t = t;
    }
    if closed {
        // (Exactly the first, for the loft to close on itself.)
        let n = sections.len();
        sections[n - 1] = Section { t: sections[n - 1].t, ..sections[0] };
        if supports[n - 1] != supports[0] {
            return Err(Error::InvalidInput("a run-over fillet whose contacts switch at the closed edge's start".into()));
        }
    }
    Ok(Rolled { sections, supports, breaks, side_set, convex, closed, tangents: [spine(t0).1, spine(t1).1] })
}

/// Whether the ball of `radius` along `edge` of `shape` touches anything but the edge's two
/// faces somewhere (where the fillet must roll on onto the faces beyond).
pub(super) fn needed(model: &Model, shape: Shape, edge: Shape, radius: f64) -> bool {
    needed_along(model, shape, edge, &|_| radius)
}

/// [`needed`] for the ball of radius `radius(s)` (`s` the fraction of the edge's length).
pub(super) fn needed_along(model: &Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> bool {
    let Ok(rolled) = roll(model, shape, edge, radius, 16, Switches::Ignored) else { return false };
    // (Not at the ends of an open edge, where the contacts meet the faces there.)
    let n = rolled.supports.len();
    let inner = if rolled.closed { &rolled.supports[..] } else { &rolled.supports[1..n - 1] };
    inner.iter().any(|s| s.iter().any(|&support| support != (0, None)))
}

/// The largest of `radius(s)` along an edge (`s` from 0 to 1), sampled.
fn largest(radius: &dyn Fn(f64) -> f64) -> f64 {
    (0..=32).map(|i| radius(i as f64 / 32.0)).fold(0.0, f64::max)
}

/// `shape` with `edge` rounded by a ball of `radius` rolling on the faces at the edge and on
/// beyond them (see the [module](self) documentation).
pub(super) fn fillet(model: &mut Model, shape: Shape, edge: Shape, radius: f64) -> Result<Shape> {
    fillet_along(model, shape, edge, &|_| radius)
}

/// [`fillet`] by a ball of radius `radius(s)` (`s` the fraction of the edge's length from
/// its start).
pub(super) fn fillet_along(model: &mut Model, shape: Shape, edge: Shape, law: &dyn Fn(f64) -> f64) -> Result<Shape> {
    // (Its largest radius, for scales.)
    let radius = largest(law);
    // Built directly (on a copy, whose edges it splits), where it can be; otherwise by
    // boolean operations.
    if let Ok(result) = directly(model, shape, edge, law) {
        return Ok(result);
    }
    // (Finer along curved edges, where the faces through the sections are smooth.)
    let straight = matches!(model.edge_curve(edge), Some((Curve::Line(_), _, _)));
    let rolled = roll(model, shape, edge, law, if straight { 48 } else { 192 }, Switches::Bisected)?;
    by_booleans(model, shape, rolled, radius)
}

/// The fillet of `shape` the ball `rolled` makes (its largest radius `radius`), built by
/// boolean operations: the wedge cut away, the lens within the solid kept, and the two
/// sewn (see the [module](self) documentation).
fn by_booleans(model: &mut Model, shape: Shape, rolled: Rolled, radius: f64) -> Result<Shape> {
    let Rolled { mut sections, supports, breaks, side_set, convex, closed, tangents } = rolled;
    if !convex && !closed {
        return Err(Error::InvalidInput("run-over fillets of open concave edges are not supported yet".into()));
    }
    // (Faces through the sections lie on what the ball touches where that is a plane or a
    // rim, ruled ones exactly. Along a curved face, smooth ones, on it within the
    // interpolation's error.)
    let curved = supports.iter().any(|s| (0..2).any(|k| s[k].1.is_none() && !matches!(side_set[k].0[s[k].0].surface, Surface::Plane(_))));
    // (Not over a rim with more of the solid beyond, where the face there is not a quadric:
    // the faces through the sections would run along it too near to cut cheaply.)
    for (section, support) in sections.iter().zip(&supports) {
        for k in 0..2 {
            let (face, Some(rim)) = support[k] else { continue };
            let quadric = |s: &Surface| matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_));
            if side_set[k].0[face].beyond[rim].iter().all(quadric) {
                continue;
            }
            let past = section.contacts[k] + (section.contacts[k] - section.center).normalized() * (1e-3 * section.radius);
            if !matches!(boolean::classify_point(model, shape, past), Ok(boolean::PointLocation::Outside)) {
                return Err(Error::InvalidInput("run-over fillets past bosses that are not quadrics are not supported yet".into()));
            }
        }
    }
    let (sections, breaks) = if closed {
        (sections, breaks)
    } else {
        // Past the ends, the end sections carried on along the edge's tangent there.
        let (first, last) = (sections[0], sections[sections.len() - 1]);
        // (Brought into the plane square to the edge there: a loft's end sections are planar.)
        let carry = |s: Section, by: Vector| {
            let t = by.normalized();
            let flat = |p: Point| p - t * (p - s.at).dot(t);
            Section { at: s.at + by, center: flat(s.center) + by, contacts: [flat(s.contacts[0]) + by, flat(s.contacts[1]) + by], ..s }
        };
        // (Where the fillet narrows to nothing toward an end, only its last section there:
        // the lens closes to its point once.)
        let mut breaks = breaks;
        while sections.len() > 2 && narrows(&sections[sections.len() - 1]) && narrows(&sections[sections.len() - 2]) {
            let gone = sections.len() - 2;
            sections.remove(gone);
            breaks.retain(|&b| b != gone);
            for b in breaks.iter_mut().filter(|b| **b > gone) {
                *b -= 1;
            }
        }
        while sections.len() > 2 && narrows(&sections[0]) && narrows(&sections[1]) {
            sections.remove(1);
            breaks.retain(|&b| b != 1);
            for b in breaks.iter_mut().filter(|b| **b > 1) {
                *b -= 1;
            }
        }
        let (first, last) = (sections[0], sections[sections.len() - 1]);
        // (Not where the fillet narrows to nothing at the end, its contacts meeting there (a
        // face narrowing to the edge's end): the lens ends at their point, and nothing lies
        // past it to cut.)
        // (Nor where the solid goes on past the end (the edge ending on a face across it, as
        // a pocket's rim at the pocket's inner corner): the wedge carried on would cut into it.
        // Told by the middle of the lens's chord a quarter of the radius on.)
        let goes_on = |model: &Model, s: &Section, by: Vector| {
            let p = s.contacts[0] + (s.contacts[1] - s.contacts[0]) * 0.5 + by.normalized() * (0.25 * s.radius);
            matches!(boolean::classify_point(model, shape, p), Ok(boolean::PointLocation::Inside))
        };
        let (narrow_first, narrow_last) = (narrows(&first) || goes_on(model, &first, -tangents[0]), narrows(&last) || goes_on(model, &last, tangents[1]));
        if !narrow_first {
            sections.insert(0, carry(first, tangents[0] * -first.radius));
        }
        if !narrow_last {
            sections.push(carry(last, tangents[1] * last.radius));
        }
        let shift = usize::from(!narrow_first);
        (sections, breaks.into_iter().map(|b| b + shift).collect())
    };
    // Only the sections the faces need: where the fillet runs straight (the ball on planes),
    // its ends; and where the contacts switch.
    let (sections, breaks) = if curved { (sections, breaks) } else { thinned(&sections, &breaks, 1e-7 * radius.max(1.0), !closed) };

    // The wedge: from one contact out, round beyond the edge and in to the other, closed by
    // the chord between them; and the lens between that chord and the ball's arc. (Out
    // from the contacts the way from the ball's centre to the edge: across whatever faces
    // meet there, and smoothly along the edge, unlike the faces' normals where the ball
    // rolls over a rim.) Around a closed edge, a ring of them.
    let mut wedge_wires = Vec::new();
    let mut lens_wires = Vec::new();
    for s in &sections[..if closed { sections.len() - 1 } else { sections.len() }] {
        let toward = (s.at - s.center).normalized();
        // (Further out where a contact is the edge's point itself, as where the ball reaches
        // the end of a face narrowing to nothing: the point out from it is not this one.)
        let at_edge = s.contacts.iter().any(|c| c.distance(s.at) <= 1e-3 * s.radius);
        let beyond = s.at + toward * (if at_edge { 0.75 } else { 0.5 } * s.radius);
        if convex {
            let out = |c: Point| c + toward * (0.5 * s.radius);
            wedge_wires.push(polygon(model, &[s.contacts[0], out(s.contacts[0]), beyond, out(s.contacts[1]), s.contacts[1]])?);
            lens_wires.push(if narrows(s) { make::vertex(model, s.contacts[1]) } else { lens(model, s)? });
        } else {
            // (At a concave edge, the wedge along a line from within the ball out through
            // each contact, crossing the face there: the solid is cut exactly along the
            // contacts. The line leans from the ball's radius toward the edge, into the
            // solid even where the ball stands over a rim. And the ball's arc alone.)
            let across = |c: Point| ((c - s.center).normalized() + toward).normalized() * s.radius;
            let (inner, deep) = (|c: Point| c - across(c) * 0.25, |c: Point| c + across(c) * 0.5);
            wedge_wires.push(polygon(model, &[inner(s.contacts[0]), deep(s.contacts[0]), beyond, deep(s.contacts[1]), inner(s.contacts[1])])?);
            lens_wires.push(arc(model, s)?);
        }
    }
    if closed {
        wedge_wires.push(wedge_wires[0]);
        lens_wires.push(lens_wires[0]);
    }
    // (Both the same way: then the chords are the same faces in both, and their edges along
    // the contacts lie on the faces there.)
    let options = loft::LoftOptions { solid: !closed, ruled: !curved, spacing: loft::Spacing::ByMeanDistance, ..Default::default() };
    let solid = |model: &mut Model, wires: &[Shape], what: &str| -> Result<Shape> {
        let lofted = loft::loft_breaking(model, wires, &options, &breaks).map_err(|e| Error::Failed(format!("lofting the {what}: {e}")))?;
        if !closed {
            return Ok(lofted);
        }
        let shells = model.unique_subshapes(lofted, NodeKind::Shell);
        let [shell] = shells[..] else {
            return Err(Error::Failed(format!("a {what} ring of {} shells", shells.len())));
        };
        solid_of(model, shell)
    };
    let wedge = solid(model, &wedge_wires, "wedge")?;
    // (A face on one of `walls`: its corners and middle on its surface, the boolean
    // operations rebuilding flat ones as planes.)
    let on_wall = |model: &Model, f: Shape, walls: &[Surface]| {
        let Some(surface) = model.surface(f).map(|s| s.transformed(&f.placement)) else { return false };
        let (u0, u1, v0, v1) = model.face_uv_bounds(f);
        let mut points: Vec<Point> = model.unique_subshapes(f, NodeKind::Vertex).iter().map(|&v| model.point(v)).collect();
        points.push(surface.point(0.5 * (u0 + u1), 0.5 * (v0 + v1)));
        // (And its edges' middles: the middle of a lens's end, its parameters' box's, may lie
        // on the chord, its arc's does not.)
        for e in model.unique_subshapes(f, NodeKind::Edge) {
            if let Some((c, a, b)) = model.edge_curve(Shape { orientation: Orientation::Forward, ..e }) {
                points.push(c.point(0.5 * (a + b)));
            }
        }
        let tol = 1e-6 * (1.0 + radius);
        walls.iter().any(|w| points.iter().all(|&p| nearest_on_surface(p, w, None).2 <= tol))
    };
    let keep = |model: &Model, s: Shape, walls: &[Surface]| -> Vec<Shape> {
        model.unique_subshapes(s, NodeKind::Face).into_iter().filter(|&f| !on_wall(model, f, walls)).collect()
    };
    // (The cut built once: only some of its faces are kept, and what is sewn from them is
    // checked below.)
    let cut = boolean::difference_once(model, shape, wedge).map_err(|e| Error::Failed(format!("cutting the wedge: {e}")))?;
    // (Its shells left open, as where the wedge's faces meet the solid's along curves
    // only nearly the same: with the looser tolerances too.)
    let cut = if closed_shells(model, cut) { cut } else { boolean::difference(model, shape, wedge).map_err(|e| Error::Failed(format!("cutting the wedge: {e}")))? };
    // (Or nothing cut at all, the wedge's edges along the contacts not found to split the
    // faces they lie on (as where a contact ends on an edge between tangent faces): with the
    // looser tolerances.)
    let untouched = (super::volume_of(model, cut) - super::volume_of(model, shape)).abs() <= 1e-9 * super::volume_of(model, shape).abs().max(1.0);
    let cut = if untouched { crate::tolerance::with_loosened(true, || boolean::difference_once(model, shape, wedge)).map_err(|e| Error::Failed(format!("cutting the wedge: {e}")))? } else { cut };
    // (At a convex edge, its shells still open along a face kept, not only along the
    // wedge's chords that are left out: the lens, sewn on along the contacts, does not
    // close that, and the pieces sewn would not close (on the public CadQuery scripts, never
    // once), found only slowly.)
    if convex && !closed_shells(model, cut) {
        let kept = keep(model, cut, &chords(model, wedge, &sections));
        let mut uses: crate::hashing::FastMap<crate::topology::NodeId, usize> = Default::default();
        for f in model.unique_subshapes(cut, NodeKind::Face) {
            for e in model.explore(f, NodeKind::Edge) {
                *uses.entry(e.node).or_default() += 1;
            }
        }
        if kept.iter().any(|&f| model.explore(f, NodeKind::Edge).any(|e| uses.get(&e.node) == Some(&1))) {
            return Err(Error::Failed("a run-over fillet whose wedge's cut leaves the shell open".into()));
        }
    }
    let faces = if convex {
        let lens_solid = solid(model, &lens_wires, "lens")?;
        let (wedge_walls, lens_walls) = (chords(model, wedge, &sections), chords(model, lens_solid, &sections));
        // The lens within the solid (trimmed by the faces at the edge's ends, and where the
        // ball hangs over a hole or past a rim, by the faces there). Both less their chords,
        // sewn.
        let trimmed = boolean::intersection(model, lens_solid, shape).map_err(|e| Error::Failed(format!("trimming the lens: {e}")))?;
        let mut faces = keep(model, cut, &wedge_walls);
        faces.extend(keep(model, trimmed, &lens_walls));
        faces
    } else {
        // At a concave edge, what the wedge takes away is for the ball's arc to bound (the
        // wedge lies within the solid but for its part in the ball): the solid's faces as
        // the wedge trims them, and the arc, sewn.
        let arcs = loft::loft_breaking(model, &lens_wires, &loft::LoftOptions { solid: false, ruled: !curved, spacing: loft::Spacing::ByMeanDistance, ..Default::default() }, &breaks).map_err(|e| Error::Failed(format!("lofting the fillet: {e}")))?;
        let walls: Vec<Surface> = model.unique_subshapes(wedge, NodeKind::Face).iter().filter_map(|&f| model.surface(f).map(|s| s.transformed(&f.placement))).collect();
        let mut faces = keep(model, cut, &walls);
        faces.extend(model.unique_subshapes(arcs, NodeKind::Face));
        faces
    };
    let sewn = sew::sew(model, &faces, &sew::SewOptions { tolerance: 1e-5 * radius.max(1.0), ..Default::default() })?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let [shell] = shells[..] else {
        return Err(Error::Failed(format!("a run-over fillet whose pieces do not close ({} shells)", shells.len())));
    };
    // (Open, its pieces not meeting edge to edge, it cannot come out right: told before
    // the faces are merged and the result measured and checked.)
    let mut uses: std::collections::HashMap<crate::topology::NodeId, (usize, usize)> = std::collections::HashMap::new();
    for e in model.explore(shell, NodeKind::Edge) {
        if !model.is_degenerate(e) {
            let used = uses.entry(e.node).or_default();
            if e.orientation == Orientation::Reversed { used.1 += 1 } else { used.0 += 1 }
        }
    }
    if uses.values().any(|&used| used != (1, 1)) {
        return Err(Error::Failed("a run-over fillet whose pieces leave the shell open".into()));
    }
    let solid = solid_of(model, shell)?;
    let result = unify::unify_same_domain(model, solid)?;
    // (Edges sewn from lofted pieces and the solid's faces stray from their faces' pcurves
    // by a little more than the classifying distance: made as tolerant as that, up to a
    // hundred thousandth of the radius.)
    tolerant_of_pcurves(model, result, 1e-5 * (1.0 + radius))?;
    // (Checked: valid, and what it takes away (at a concave edge, adds) no more than the
    // wedge, nor less than nothing.)
    let (before, after, wedge_volume) = (super::volume_of(model, shape), super::volume_of(model, result), super::volume_of(model, wedge));
    let slack = 1e-6 * before.abs().max(1.0);
    let change = if convex { before - after } else { after - before };
    if !crate::validity::check(model, result).is_valid() || change < -slack || change > wedge_volume + slack {
        return Err(Error::Failed("a run-over fillet that does not come out right".into()));
    }
    Ok(result)
}

/// `shape` with the chain `edges` (going on into one another smoothly at shared vertices)
/// rounded by one ball of `radius` rolling along them all, on the faces at the edges and on
/// beyond them, built by boolean operations (see the [module](self) documentation).
pub(super) fn fillet_chain(model: &mut Model, shape: Shape, edges: &[Shape], radius: f64) -> Result<Shape> {
    let curved = edges.iter().any(|&e| !matches!(model.edge_curve(e), Some((Curve::Line(_), _, _))));
    let (rolled, _) = roll_chain(model, shape, edges, &|_| radius, if curved { 96 } else { 48 }, Switches::Bisected)?;
    by_booleans(model, shape, rolled, radius)
}

/// Whether the ball of `radius` along the chain `edges` of `shape` touches anything but the
/// faces at the edges somewhere (away from an open chain's ends).
pub(super) fn chain_needed(model: &Model, shape: Shape, edges: &[Shape], radius: f64) -> bool {
    let rolled = roll_chain(model, shape, edges, &|_| radius, 8, Switches::Ignored);
    let Ok((rolled, own)) = rolled else { return false };
    // (Each section's contacts on the faces at the edge it lies along; at a vertex between
    // two, on either's.)
    let n = rolled.supports.len();
    let range = if rolled.closed { 0..n } else { 1..n - 1 };
    range.into_iter().any(|i| {
        let (s, t) = (rolled.supports[i], rolled.sections[i].t);
        let last = own.len() - 1;
        let pieces = [(t.floor().max(0.0) as usize).min(last), ((t.ceil().max(1.0) as usize) - 1).min(last)];
        (0..2).any(|k| s[k].1.is_some() || !pieces.iter().any(|&p| own[p][k] == s[k].0))
    })
}

/// Whether `edges` make one chain (see [`chain_order`]).
pub(super) fn is_chain(model: &Model, edges: &[Shape]) -> bool {
    chain_order(model, edges).is_some()
}

/// `edges` in order along the chain they make (each next one from the vertex the last ends
/// at), each with whether the chain runs against it; and whether the chain closes. `None`
/// where they are not one chain (three meeting at a vertex, or apart).
fn chain_order(model: &Model, edges: &[Shape]) -> Option<(Vec<(Shape, bool)>, bool)> {
    use crate::topology::NodeId;
    let ends: Vec<(NodeId, NodeId)> = edges
        .iter()
        .map(|&e| {
            let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..e });
            (a.node, b.node)
        })
        .collect();
    let n = edges.len();
    let mut uses: crate::hashing::FastMap<NodeId, usize> = Default::default();
    for &(a, b) in &ends {
        *uses.entry(a).or_default() += 1;
        *uses.entry(b).or_default() += 1;
    }
    if n == 0 || uses.values().any(|&u| u > 2) {
        return None;
    }
    let start = (0..n).find_map(|i| if uses[&ends[i].0] == 1 { Some((i, false)) } else if uses[&ends[i].1] == 1 { Some((i, true)) } else { None });
    let closed = start.is_none();
    let (mut i, mut against) = start.unwrap_or((0, false));
    let mut used = vec![false; n];
    let mut order = Vec::with_capacity(n);
    loop {
        used[i] = true;
        order.push((edges[i], against));
        let tip = if against { ends[i].0 } else { ends[i].1 };
        match (0..n).find(|&j| !used[j] && (ends[j].0 == tip || ends[j].1 == tip)) {
            Some(j) => {
                against = ends[j].0 != tip;
                i = j;
            }
            None => break,
        }
    }
    (order.len() == n).then_some((order, closed))
}

/// A face's outward normal at its surface's point nearest `p`.
fn outward_normal(model: &Model, face: Shape, p: Point) -> Option<Vector> {
    let surface = model.surface(face)?.transformed(&face.placement);
    let (u, v, _) = nearest_on_surface(p, &surface, None);
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    if n.length() <= 1e-300 {
        return None;
    }
    let n = n.normalized();
    Some(if face.orientation == Orientation::Reversed { -n } else { n })
}

/// [`roll`] along the chain `edges`: its sides each the faces at the edges on that side (one
/// edge's face going on into the next edge's on the same side, the same face or tangent to
/// it), with their neighbours. Also, for each edge, its faces' places among its sides'.
fn roll_chain(model: &Model, shape: Shape, edges: &[Shape], radius: &dyn Fn(f64) -> f64, count: usize, switches: Switches) -> Result<(Rolled, Vec<[usize; 2]>)> {
    let (order, closed) = chain_order(model, edges).ok_or_else(|| Error::InvalidInput("fillet edges that are not one chain".into()))?;
    let mut pieces = Vec::with_capacity(order.len());
    let mut sided: Vec<[(Shape, f64); 2]> = Vec::with_capacity(order.len());
    for &(edge, against) in &order {
        let faces = faces_of(model, shape, edge);
        let [(g0, a0), (g1, a1)] = faces[..] else {
            return Err(Error::InvalidInput(format!("a fillet edge must bound two faces, it bounds {}", faces.len())));
        };
        let (curve, t0, t1) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
        let (a, b, sign) = if against { (t1, t0, -1.0) } else { (t0, t1, 1.0) };
        let mut pair = [(g0, a0 * sign), (g1, a1 * sign)];
        if let Some(previous) = sided.last() {
            let same = |x: Shape, y: Shape| x.node == y.node;
            let keep = if same(pair[0].0, previous[0].0) || same(pair[1].0, previous[1].0) {
                true
            } else if same(pair[0].0, previous[1].0) || same(pair[1].0, previous[0].0) {
                false
            } else {
                let p = curve.point(a);
                let normal = |f: Shape| outward_normal(model, f, p).unwrap_or(Vector::new(0.0, 0.0, 0.0));
                let n0 = normal(pair[0].0);
                n0.dot(normal(previous[0].0)) >= n0.dot(normal(previous[1].0))
            };
            if !keep {
                pair.swap(0, 1);
            }
        }
        pieces.push((curve, a, b));
        sided.push(pair);
    }
    let chain: Vec<crate::topology::NodeId> = order.iter().map(|(e, _)| e.node).collect();
    let ends: Vec<crate::topology::NodeId> = if closed {
        Vec::new()
    } else {
        let vertex = |(e, against): (Shape, bool), last: bool| {
            let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..e });
            if against == last { a.node } else { b.node }
        };
        vec![vertex(order[0], false), vertex(order[order.len() - 1], true)]
    };
    let all = model.unique_subshapes(shape, NodeKind::Face);
    let primary = |k: usize| {
        let mut faces: Vec<Shape> = Vec::new();
        for pair in &sided {
            if !faces.iter().any(|f| f.node == pair[k].0.node) {
                faces.push(pair[k].0);
            }
        }
        faces
    };
    let side = |k: usize| -> (Side, Vec<usize>) {
        let (mine, other) = (primary(k), primary(1 - k));
        let mut members = mine.clone();
        for face in &mine {
            for e in model.unique_subshapes(*face, NodeKind::Edge) {
                if chain.contains(&e.node) {
                    continue;
                }
                let (p, q) = model.edge_vertices(e);
                if ends.contains(&p.node) || ends.contains(&q.node) {
                    continue;
                }
                for &f in &all {
                    let fresh = !members.iter().any(|m| m.node == f.node) && !other.iter().any(|o| o.node == f.node);
                    if fresh && model.unique_subshapes(f, NodeKind::Edge).iter().any(|x| x.node == e.node) {
                        members.push(f);
                    }
                }
            }
        }
        // (Each edge's face on this side, by its place among the side's faces measured.)
        let kept: Vec<Shape> = members.iter().copied().filter(|&f| Bounded::of(model, f).is_some()).collect();
        let own: Vec<usize> = sided.iter().map(|pair| kept.iter().position(|f| f.node == pair[k].0.node).unwrap_or(usize::MAX)).collect();
        let bounded: Vec<Bounded> = kept.iter().filter_map(|&f| Bounded::of(model, f)).collect();
        (Side(bounded), own)
    };
    let ((s0, own0), (s1, own1)) = (side(0), side(1));
    let spine = Spine::new(pieces, closed);
    let rolled = roll_on(model, &spine, [s0, s1], sided[0], radius, count, switches)?;
    Ok((rolled, own0.into_iter().zip(own1).map(|(a, b)| [a, b]).collect()))
}

/// The edges of `shape` made as tolerant as their faces' pcurves stray from them, where that
/// is more than their tolerance and no more than `limit`; the vertices then as tolerant as
/// their edges.
fn tolerant_of_pcurves(model: &mut Model, shape: Shape, limit: f64) -> Result<()> {
    let mut wider: Vec<(crate::topology::EdgeId, f64)> = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face).cloned() else { continue };
        let coedges: Vec<Shape> = model.children(face).flat_map(|w| model.children(w).collect::<Vec<_>>()).collect();
        for coedge in coedges {
            let Ok(c) = CoEdgeId::try_from(coedge.node) else { continue };
            let edge = model[c].edge;
            let Some((pcurve, p0, p1)) = model.effective_pcurve(c) else { continue };
            // (At the same parameters on both, as the edge's own curve has them; a degenerate
            // edge's point, its vertex.)
            let curve = model.get(edge).and_then(|e| model.get(e.curve)).map(|slot| slot.curve.clone());
            let (start, _) = model.edge_vertices(Shape::from(edge));
            let pole = model.point(start);
            let stray = (0..=16)
                .map(|k| {
                    let t = p0 + (p1 - p0) * k as f64 / 16.0;
                    let uv = pcurve.point(t);
                    surface.point(uv.x, uv.y).distance(curve.as_ref().map_or(pole, |c| c.point(t)))
                })
                .fold(0.0, f64::max);
            if stray > model.tolerance(Shape::from(edge)) && stray <= limit {
                wider.push((edge, 1.5 * stray));
            }
        }
    }
    for (edge, tolerance) in wider {
        if tolerance > model.tolerance(Shape::from(edge)) {
            model.set_edge_tolerance(edge, tolerance)?;
        }
    }
    crate::modify::update_tolerances(model, shape)
}

/// [`fillet_along`] built directly (see [`trim`]).
fn directly(model: &mut Model, shape: Shape, edge: Shape, law: &dyn Fn(f64) -> f64) -> Result<Shape> {
    let (copy, edge) = copied(model, shape, edge)?;
    let straight = matches!(model.edge_curve(edge), Some((Curve::Line(_), _, _)));
    let Rolled { sections, supports, breaks, side_set, convex, closed, tangents } = roll(model, copy, edge, law, if straight { 32 } else { 64 }, Switches::Found)?;
    let result = trim::build(model, copy, edge, &sections, &supports, &breaks, &side_set, largest(law), convex, closed, tangents)?;
    // (Checked: what is kept or dropped amiss leaves the shell open. Built this way it is
    // turned out, its faces as they were and the fillet's facing out.)
    if !closed_shells(model, result) {
        return Err(Error::Failed("a run-over fillet that does not come out right".into()));
    }
    // (The contact curves' pcurves pass through their sections only, straying between by what
    // the fit leaves.)
    tolerant_of_pcurves(model, result, 1e-5 * (1.0 + largest(law)))?;
    Ok(result)
}

/// Whether every shell of `shape` is closed: each of its edges run along once each way.
fn closed_shells(model: &Model, shape: Shape) -> bool {
    model.unique_subshapes(shape, NodeKind::Shell).into_iter().all(|shell| {
        let mut uses: crate::hashing::FastMap<crate::topology::EdgeId, (usize, usize)> = Default::default();
        for face in model.children(shell) {
            for wire in model.children(face) {
                for coedge in model.children(wire) {
                    let Ok(id) = CoEdgeId::try_from(coedge.node) else { return false };
                    let entry = uses.entry(model[id].edge).or_default();
                    if model[id].orientation.compose(coedge.orientation) == Orientation::Reversed { entry.1 += 1 } else { entry.0 += 1 }
                }
            }
        }
        !uses.is_empty() && uses.values().all(|&n| n == (1, 1))
    })
}

/// A copy of `shape`, and its edge where `edge` of `shape` is.
fn copied(model: &mut Model, shape: Shape, edge: Shape) -> Result<(Shape, Shape)> {
    let at = |model: &Model, e: Shape| model.edge_curve(Shape { orientation: Orientation::Forward, ..e }).map(|(c, a, b)| [c.point(a), c.point(0.5 * (a + b))]);
    let [start, middle] = at(model, edge).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
    let copy = crate::modify::copy(model, shape, false)?;
    let near = 1e-9 * (1.0 + start.distance(middle));
    let found = model
        .unique_subshapes(copy, NodeKind::Edge)
        .into_iter()
        .find(|&e| at(model, e).is_some_and(|[s, m]| s.distance(start) <= near && m.distance(middle) <= near))
        .ok_or_else(|| Error::Failed("the edge not found in a copy of its shape".into()))?;
    Ok((copy, Shape { orientation: Orientation::Forward, ..found }))
}

/// The solid bounded by the closed `shell` (its faces turned out).
pub(super) fn solid_of(model: &mut Model, shell: Shape) -> Result<Shape> {
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    if super::volume_of(model, solid) >= 0.0 {
        return Ok(solid);
    }
    let turned = model.new_empty(NodeKind::Solid);
    model.add(turned, shell.reversed())?;
    Ok(turned)
}

/// The surfaces of a loft's faces along the chords between the contacts: those all of
/// whose vertices are contacts of the sections, of two sections or more (not an end cap).
fn chords(model: &Model, solid: Shape, sections: &[Section]) -> Vec<Surface> {
    let size = sections[0].at.distance(sections[sections.len() - 1].at) + sections[0].at.distance(sections[0].center);
    let section_of = |p: Point| sections.iter().position(|s| s.contacts.iter().any(|&c| c.distance(p) <= 1e-7 * (1.0 + size)));
    model
        .unique_subshapes(solid, NodeKind::Face)
        .into_iter()
        .filter(|&f| {
            let found: Vec<Option<usize>> = model.unique_subshapes(f, NodeKind::Vertex).iter().map(|&v| section_of(model.point(v))).collect();
            let mut at: Vec<usize> = found.iter().flatten().copied().collect();
            at.sort_unstable();
            at.dedup();
            if !(found.iter().all(Option::is_some) && at.len() >= 2) {
                return false;
            }
            // (Through the middle of a chord it spans, unlike the arc's faces.)
            let Some(surface) = model.surface(f).map(|s| s.transformed(&f.placement)) else { return false };
            let s = &sections[at[0]];
            let middle = s.contacts[0] + (s.contacts[1] - s.contacts[0]) * 0.5;
            nearest_on_surface(middle, &surface, None).2 <= 1e-6 * (1.0 + size)
        })
        .filter_map(|f| model.surface(f).cloned())
        .collect()
}

/// Whether the fillet narrows to nothing at section `s`: its contacts (nearly) one point.
fn narrows(s: &Section) -> bool {
    s.contacts[0].distance(s.contacts[1]) <= 1e-3 * s.radius
}

/// A closed wire of straight edges through `points`.
/// The sections that stay (with the breaks among them renumbered): the ends and breaks,
/// and between them those off the straight line from one kept section to the next by more
/// than `tolerance`.
fn thinned(sections: &[Section], breaks: &[usize], tolerance: f64, ends: bool) -> (Vec<Section>, Vec<usize>) {
    let points = |s: &Section| [s.at, s.center, s.contacts[0], s.contacts[1]];
    // Whether those strictly between `k` and `j` lie on the lines from one to the other.
    let straight = |k: usize, j: usize| {
        let (a, b) = (points(&sections[k]), points(&sections[j]));
        let along = (sections[j].at - sections[k].at).length();
        (k + 1..j).all(|i| {
            let fraction = if along > 0.0 { (sections[i].at - sections[k].at).length() / along } else { 0.0 };
            points(&sections[i]).iter().zip(a.iter().zip(&b)).all(|(&p, (&from, &to))| (p - (from + (to - from) * fraction)).length() <= tolerance)
        })
    };
    // (The breaks, and the sections at the edge's ends (with `ends`, past which they are
    // carried), kept: the faces there are cut by the faces at the ends.)
    let ends = if ends { vec![1, sections.len().saturating_sub(2)] } else { Vec::new() };
    let mut stops: Vec<usize> = breaks.iter().copied().chain(ends).collect();
    stops.sort_unstable();
    let mut kept = vec![0];
    let mut kept_breaks = Vec::new();
    let mut k = 0;
    while k + 1 < sections.len() {
        // As far on as stays straight, but not past a break.
        let limit = stops.iter().copied().find(|&b| b > k).unwrap_or(sections.len() - 1).min(sections.len() - 1);
        let mut j = k + 1;
        while j < limit && straight(k, j + 1) {
            j += 1;
        }
        if breaks.contains(&j) {
            kept_breaks.push(kept.len());
        }
        kept.push(j);
        k = j;
    }
    (kept.iter().map(|&i| sections[i]).collect(), kept_breaks)
}

fn polygon(model: &mut Model, points: &[Point]) -> Result<Shape> {
    let vertices: Vec<Shape> = points.iter().map(|&p| make::vertex(model, p)).collect();
    let n = vertices.len();
    let edges: Vec<Shape> = (0..n).map(|i| make::edge(model, vertices[i], vertices[(i + 1) % n])).collect::<Result<_>>()?;
    make::wire(model, &edges)
}

/// The ball's arc in a section, from the first contact to the second on the edge's side.
fn arc(model: &mut Model, s: &Section) -> Result<Shape> {
    let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
    let angle = a.cross(b).length().atan2(a.dot(b));
    let circle = Circle::new(Frame::new(s.center, Direction::from_vector(a.cross(b)), Direction::from_vector(a)), s.radius);
    let edge = make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, angle)?;
    make::wire(model, &[edge])
}

/// The lens in a section between the ball's arc (from the first contact to the second,
/// on the edge's side) and the chord back.
fn lens(model: &mut Model, s: &Section) -> Result<Shape> {
    let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
    let angle = a.cross(b).length().atan2(a.dot(b));
    let normal = a.cross(b);
    // (A section whose contacts meet, the faces tangent there: no arc across.)
    let (Some(z), Some(x)) = (Direction::try_from_vector(normal), Direction::try_from_vector(a)) else {
        return Err(Error::InvalidInput("a fillet's section has no arc across (its faces are tangent there)".into()));
    };
    let circle = Circle::new(Frame::new(s.center, z, x), s.radius);
    let arc = make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, angle)?;
    let (start, end) = model.edge_vertices(arc);
    let chord = make::edge(model, end, start)?;
    make::wire(model, &[arc, chord])
}

/// Faces as [`rolls_off`] measures them, each made once (for several edges of one face).
#[derive(Default)]
pub(super) struct KnownFaces(std::cell::RefCell<crate::hashing::FastMap<crate::topology::NodeId, Option<std::rc::Rc<Bounded>>>>);

impl KnownFaces {
    fn get(&self, model: &Model, face: Shape) -> Option<std::rc::Rc<Bounded>> {
        self.0.borrow_mut().entry(face.node).or_insert_with(|| Bounded::of_with(model, face, false).map(std::rc::Rc::new)).clone()
    }
}

/// Whether the ball of the fillet along a straight edge between `faces`, touching them at
/// `contacts` in the section at the edge's start (the edge running along `tangent` for
/// `length`, from `ends[0]` to `ends[1]`), touches them off their bounds somewhere along it:
/// where it must roll onto the faces beyond. (Not where a contact leaves its face across an
/// edge from one of the edge's ends, as a face narrowing to its end: the fillet's end meets
/// the faces there as at any end. Unless `narrowing`: for an edge of a chain, whose ends are
/// not the fillet's.)
pub(super) fn rolls_off(model: &Model, faces: [Shape; 2], tangent: Vector, length: f64, contacts: [Point; 2], ends: [Point; 2], known: &KnownFaces, narrowing: bool) -> bool {
    const SAMPLES: usize = 32;
    faces.iter().zip(contacts).any(|(&face, contact)| {
        let Some(bounded) = known.get(model, face) else { return false };
        // (Whether the face's boundary edge nearest `p` runs from one of the edge's ends.)
        let from_an_end = |p: Point| {
            let gap = |&(ref c, first, last): &(Curve, f64, f64)| (0..=64).map(|i| c.point(first + (last - first) * i as f64 / 64.0).distance(p)).fold(f64::INFINITY, f64::min);
            let Some(nearest) = (0..bounded.edges.len()).min_by(|&i, &j| gap(&bounded.edges[i]).total_cmp(&gap(&bounded.edges[j]))) else { return false };
            let (c, first, last) = &bounded.edges[nearest];
            [c.point(*first), c.point(*last)].iter().any(|q| ends.iter().any(|e| q.distance(*e) <= crate::tolerance::LINEAR))
        };
        (0..=SAMPLES).any(|i| {
            // (Not at the very ends, where the contact meets the faces there.)
            let s = length * (0.02 + 0.96 * i as f64 / SAMPLES as f64);
            let p = contact + tangent.normalized() * s;
            let (u, v, d) = nearest_on_surface(p, &bounded.surface, None);
            d <= 1e-6 * (1.0 + length) && !bounded.contains((u, v)) && (narrowing || !from_an_end(p))
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cylinder;

    /// A disc's middle level with where its circle starts and ends: on the face, though the
    /// circle's two ends lie a rounding error apart there and a ray from the point along
    /// that level passes between them.
    #[test]
    fn a_point_level_with_where_a_loop_closes_lies_within_it() {
        let mut model = Model::new();
        let solid = Cylinder::new(5.0, 12.0).build(&mut model).unwrap();
        let discs: Vec<Shape> = model.unique_subshapes(solid, NodeKind::Face).into_iter().filter(|&f| matches!(model.surface(f), Some(Surface::Plane(_)))).collect();
        assert_eq!(discs.len(), 2);
        for disc in discs {
            let bounded = Bounded::of_with(&model, disc, false).unwrap();
            for p in [Point::new(2.0, 0.0, 0.0), Point::new(-2.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0)] {
                let (u, v, _) = nearest_on_surface(p, &bounded.surface, None);
                assert_eq!(bounded.locate((u, v)), Some(true), "{p:?}");
            }
            let (u, v, _) = nearest_on_surface(Point::new(7.0, 0.0, 0.0), &bounded.surface, None);
            assert_eq!(bounded.locate((u, v)), Some(false));
        }
    }
}
