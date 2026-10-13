//! Fillets made by a ball rolling along any edge between two faces.
//!
//! The ball's centre is found section by section along the edge (continued a little past
//! its ends): in the plane square to the edge, at the ball's radius from both faces. The
//! fillet surface interpolates the arcs between the contact points (exact circles across,
//! a cubic B-spline along the edge, rational). The solid is then rebuilt around it: the
//! two faces lose the strips between the edge and the contact curves, the faces at the
//! edge's ends lose (or gain) the corner cut off by the fillet surface, and the edges
//! there are shortened to where the contact curves meet them.

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{BSplineCurve2, Curve2, Line2};
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{fit, project, ParametricCurve, ParametricSurface};
use crate::math::{Direction2, Point, Point2, Vector};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

#[path = "network.rs"]
mod network;
pub(super) use network::fillet_network;

/// Sections along the edge, and past each of its ends.
const SECTIONS: usize = 48;
const PAST_END: usize = 8;
/// Samples across the fillet where it meets an end face.
const ACROSS: usize = 16;

/// A face along the edge: its surface, +1 where the surface's normal points out of the
/// solid (-1 otherwise), and the direction of the edge in its boundary (+1 if the face's
/// boundary runs along the edge's own direction).
#[derive(Clone)]
struct Side {
    face: Shape,
    surface: Surface,
    sign: f64,
    along: f64,
    /// The parameters of the last point projected (where the surface has no closed-form
    /// inverse, the next projection starts there).
    last: std::cell::Cell<Option<(f64, f64)>>,
    /// The face's range of `v`: a surface of revolution's inverse may take a point to its
    /// meridian's half across the axis (the same point, with the normal the other way).
    v_range: (f64, f64),
}

impl Side {
    /// The signed distance from the surface (positive outside the solid), the foot point
    /// and the outward unit normal there.
    /// The parameters of `p`'s foot on the surface.
    fn uv(&self, p: Point) -> Result<Point2> {
        let (u, v) = self.parameters(p);
        Ok(Point2::new(u, v))
    }

    /// The parameters of `p`'s foot on the surface: in closed form where there is one,
    /// otherwise the nearest point (from the last one found).
    fn parameters(&self, p: Point) -> (f64, f64) {
        if let Some(uv) = self.surface.parameters(p) {
            let (v0, v1) = self.v_range;
            let margin = 0.25 * (v1 - v0);
            let outside = !(uv.1 >= v0 - margin && uv.1 <= v1 + margin);
            if !(matches!(self.surface, Surface::Revolution(_)) && outside && (v1 - v0).is_finite()) {
                return uv;
            }
            // (Found again from within the face's range.)
            let start = self.last.get().filter(|l| l.1 >= v0 - margin && l.1 <= v1 + margin).unwrap_or((uv.0 + std::f64::consts::PI, 0.5 * (v0 + v1)));
            let (u, v, _) = crate::geometry::extrema::newton_on_surface(p, &self.surface, start);
            self.last.set(Some((u, v)));
            return (u, v);
        }
        // (Locally from the last foot, the points following one another; the whole
        // surface searched the first time.)
        let (u, v, _) = match self.last.get() {
            Some(start) => crate::geometry::extrema::newton_on_surface(p, &self.surface, start),
            None => crate::geometry::extrema::nearest_on_surface(p, &self.surface, None),
        };
        self.last.set(Some((u, v)));
        (u, v)
    }

    fn distance(&self, p: Point) -> Result<(f64, Point, Vector)> {
        let (u, v) = self.parameters(p);
        let (foot, du, dv) = self.surface.d1(u, v);
        let n = du.cross(dv);
        if n.length() <= 1e-300 {
            return Err(Error::InvalidInput("a fillet section meets a singular point of a face".into()));
        }
        let n = n.normalized() * self.sign;
        Ok(((p - foot).dot(n), foot, n))
    }
}

/// A section of the fillet: the spine point on the edge (or its continuation), the ball's
/// centre and its contact points with the two faces.
struct Section {
    t: f64,
    center: Point,
    contacts: [Point; 2],
    /// The contacts' parameters on the faces.
    contact_uv: [Point2; 2],
    radius: f64,
}

/// The ball's sections along `edge` of `shape` (and past its ends), with what the fillet
/// is built from: the edge, its faces, the contact curves and whether it is convex.
fn rolling(model: &Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<Rolling> {
    let id = EdgeId::try_from(edge.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not an edge")))?;
    let forward = Shape { orientation: Orientation::Forward, ..edge };
    let (curve, t0, t1) = model.edge_curve(forward).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
    let (va, vb) = model.edge_vertices(forward);
    let closed = va.node == vb.node;
    let sides = sides_of(model, shape, id)?;
    // The fraction of the edge's length from its start, for the radius law.
    let lengths: Vec<f64> = {
        let n = 256;
        let mut acc = vec![0.0];
        for k in 1..=n {
            let (a, b) = (t0 + (t1 - t0) * (k - 1) as f64 / n as f64, t0 + (t1 - t0) * k as f64 / n as f64);
            let last = *acc.last().expect("lengths");
            acc.push(last + curve.point(a).distance(curve.point(0.5 * (a + b))) + curve.point(0.5 * (a + b)).distance(curve.point(b)));
        }
        acc
    };
    let fraction = |t: f64| {
        let x = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0) * (lengths.len() - 1) as f64;
        let k = (x.floor() as usize).min(lengths.len() - 2);
        let l = lengths[k] + (lengths[k + 1] - lengths[k]) * (x - k as f64);
        l / lengths[lengths.len() - 1]
    };

    // Convex where the material lies between the faces: the ball rolls inside it (at
    // distance -r from both faces), else outside (+r).
    let mid = 0.5 * (t0 + t1);
    let (o, d) = curve.d1(mid);
    let tangent = d.normalized();
    let normals = [sides[0].distance(o)?.2, sides[1].distance(o)?.2];
    let into = |i: usize| normals[i].cross(tangent * sides[i].along).normalized();
    if into(0).cross(into(1)).length() <= crate::tolerance::ANGULAR.sqrt() {
        return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
    }
    let convex = normals[0].dot(into(1)) < 0.0;
    let side_of_ball = if convex { -1.0 } else { 1.0 };

    // The spine: points of the edge, continued past its ends (by a few radii) along the
    // faces' intersection.
    let reach = if closed { 0.0 } else { 4.0 * (0..=8).map(|i| radius(i as f64 / 8.0)).fold(0.0, f64::max) };
    let past = |t_end: f64| -> f64 { reach / curve.d1(t_end).1.length() };
    let (t_min, t_max) = (t0 - past(t0), t1 + past(t1));
    let spine_at = |t: f64| -> Result<(Point, Vector)> {
        if (t0..=t1).contains(&t) {
            let (p, d) = curve.d1(t);
            return Ok((p, d.normalized()));
        }
        // Past an end: the edge's curve continued (keeping its parametrization smooth),
        // brought onto both faces' surfaces.
        let (p, d) = curve.d1(t);
        let q = onto_both(&sides, p, d.normalized())?;
        let n = sides[0].distance(q)?.2.cross(sides[1].distance(q)?.2);
        Ok((q, if n.dot(d) < 0.0 { -n.normalized() } else { n.normalized() }))
    };
    // (Each section from the centre of the last found, relative to its point of the spine.)
    let last_offset: std::cell::Cell<Option<Vector>> = std::cell::Cell::new(None);
    let section_at = |t: f64| -> Result<Section> {
        let r = radius(fraction(t));
        if !(r > crate::tolerance::LINEAR) {
            return Err(Error::InvalidInput(format!("fillet radius must be positive, got {r}")));
        }
        let (o, tangent) = spine_at(t)?;
        let found = match section(&sides, o, tangent, r, side_of_ball, last_offset.get()) {
            Ok(found) => Ok((o, Section { t, ..found })),
            // On the edge, where the faces' surfaces meet a seam or an end of theirs: just
            // inside it.
            Err(e) if (t0.min(t1)..=t0.max(t1)).contains(&t) => {
                let inside = t + (0.5 * (t0 + t1) - t).signum() * 1e-4 * (t1 - t0).abs();
                let (o, tangent) = spine_at(inside)?;
                section(&sides, o, tangent, r, side_of_ball, last_offset.get()).map(|found| (o, Section { t: inside, ..found })).map_err(|_| e)
            }
            Err(e) => Err(e),
        };
        let (o, s) = found?;
        last_offset.set(Some(s.center - o));
        Ok(s)
    };
    // More sections until the contact curves interpolated through them hold to fresh
    // sections midway.
    let size = curve.point(t0).distance(curve.point(0.5 * (t0 + t1))) + curve.point(t1).distance(curve.point(0.5 * (t0 + t1))) + reach;
    let target = 1e-9 * size.max(1.0);
    let mut n = SECTIONS;
    let mut previous_worst = f64::INFINITY;
    let mut previous_poor = false;
    let (sections, contacts) = loop {
        let past_n = if closed { 0 } else { (PAST_END * n / SECTIONS).max(PAST_END) };
        // Past each end, outwards as far as the ball still fits between the faces (their
        // surfaces may end or degenerate further on).
        let past = |end: f64, far: f64| -> Vec<Section> {
            (1..=past_n).map(|j| end + (far - end) * j as f64 / past_n as f64).map_while(|t| section_at(t).ok()).collect()
        };
        let before = past(t0, t_min);
        let after = past(t1, t_max);
        let mut sections: Vec<Section> = before.into_iter().rev().collect();
        for k in 0..=n {
            sections.push(section_at(t0 + (t1 - t0) * k as f64 / n as f64)?);
        }
        sections.extend(after);
        let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();
        let contact = |i: usize| -> Result<BSplineCurve> {
            let points: Vec<Point> = sections.iter().map(|s| s.contacts[i]).collect();
            fit::interpolate_with_parameters(&points, &ts)
        };
        let contacts = [contact(0)?, contact(1)?];
        let mut worst: f64 = 0.0;
        for w in ts.windows(2) {
            let mid = section_at(0.5 * (w[0] + w[1]))?;
            for i in 0..2 {
                worst = worst.max(contacts[i].point(mid.t).distance(mid.contacts[i]));
            }
        }
        // (Not further where doubling the sections has twice running hardly helped, and the
        // gap is still far from closed: the contact curves are not smooth there (the faces'
        // own seams), and twice as many would not close it. Errors jumping about instead
        // (sections finding other balls) are left to the full count.)
        let ratio = previous_worst / worst;
        let poor = (0.25..4.0).contains(&ratio);
        let stalled = n >= 4 * SECTIONS && worst > 100.0 * target && poor && previous_poor;
        if worst <= target || stalled || n >= 32 * SECTIONS {
            break (sections, contacts);
        }
        previous_worst = worst;
        previous_poor = poor;
        n *= 2;
    };
    Ok(Rolling { id, va, vb, closed, sides, sections, contacts, convex, target, t0, t1 })
}

/// What a rolling-ball fillet is made from (see [`rolling`]).
struct Rolling {
    id: EdgeId,
    va: Shape,
    vb: Shape,
    closed: bool,
    sides: [Side; 2],
    sections: Vec<Section>,
    contacts: [BSplineCurve; 2],
    convex: bool,
    target: f64,
    t0: f64,
    t1: f64,
}

/// `shape` with `edge` rounded by a ball of radius `radius(s)`, `s` the fraction of the
/// edge's length from its start.
pub(super) fn rolling_ball(model: &mut Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<Shape> {
    let Rolling { id, va, vb, closed, sides, sections, contacts, convex, target, .. } = rolling(model, shape, edge, radius)?;
    let surface = fillet_surface(&sections)?;
    let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();

    // Where the fillet meets the faces at the edge's ends.
    let mut ends = Vec::new();
    if !closed {
        for (vertex, at_end) in [(va, false), (vb, true)] {
            ends.push(end_cut(model, shape, id, vertex, &sides, &surface, &ts, at_end, target)?);
        }
    }
    let mut builder = Rebuild { model, id, sides: &sides, surface: &surface, contacts: &contacts, ends: &ends, ts: &ts };
    builder.build(shape, closed, &sections, convex)
}

/// The solid a rolling-ball fillet of `edge` takes away (at a convex edge) or adds (at a
/// concave one), for a boolean operation, with whether it is convex: in each section the
/// region between the ball's arc and the edge, bounded by the rays from the ball's centre
/// through the contacts (carried well past the edge) and closed across beyond it. A convex
/// tool runs on past the edge's ends, for the faces there to trim it; a concave one stops
/// at them.
pub(super) fn rolling_tool(model: &mut Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<(Shape, bool)> {
    use crate::geometry::curve::Circle;
    use crate::math::{Direction, Frame};
    let Rolling { sections, convex, closed, t0, t1, .. } = rolling(model, shape, edge, radius)?;
    if closed {
        return Err(Error::InvalidInput("rolling-ball tools along closed edges are not supported yet".into()));
    }
    let span = (t1 - t0).abs() * 1e-9;
    let kept: Vec<&Section> = sections.iter().filter(|s| convex || (s.t >= t0.min(t1) - span && s.t <= t0.max(t1) + span)).collect();
    let mut wires = Vec::with_capacity(kept.len());
    for s in kept {
        let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        // (Not where the faces turn tangent: the arc shrinks to nothing there.)
        if half <= 1e-3 || s.contacts[0].distance(s.contacts[1]) <= 1e3 * crate::tolerance::LINEAR {
            continue;
        }
        // Far enough along the rays for the closing segment to pass beyond the edge's
        // corner (which lies within the wedge, at r / sin(half) from the centre).
        let reach = (s.radius / half.sin().max(1e-3) + 0.5 * s.radius) / half.cos().max(1e-3);
        let far = [s.center + a.normalized() * reach, s.center + b.normalized() * reach];
        let normal = a.cross(b);
        let frame = Frame::new(s.center, Direction::from_vector(normal), Direction::from_vector(a));
        // (The lines start where the arc ends, the contacts being only nearly on it.)
        let circle = Circle::new(frame, a.length());
        let ends = [circle.point(0.0), circle.point(2.0 * half)];
        let arc = crate::make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, 2.0 * half)?;
        let edges = [
            arc,
            crate::make::edge_between(model, ends[1], far[1])?,
            crate::make::edge_between(model, far[1], far[0])?,
            crate::make::edge_between(model, far[0], ends[0])?,
        ];
        wires.push(crate::make::wire(model, &edges)?);
    }
    let tool = crate::loft::loft(model, &wires, &crate::loft::LoftOptions { solid: true, ruled: false, ..Default::default() })?.shape;
    Ok((tool, convex))
}

/// The two faces along the edge.
fn sides_of(model: &Model, shape: Shape, id: EdgeId) -> Result<[Side; 2]> {
    let mut sides = Vec::new();
    for face in model.explore(shape, NodeKind::Face) {
        for wire in model.children(face) {
            for coedge in model.children(wire) {
                let c = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                if model[c].edge != id {
                    continue;
                }
                let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.clone();
                let sign = |o: Orientation| if o == Orientation::Reversed { -1.0 } else { 1.0 };
                sides.push(Side {
                    face,
                    surface,
                    sign: sign(face.orientation),
                    along: sign(model[c].orientation.compose(wire.orientation)),
                    last: Default::default(),
                    v_range: { let b = model.face_uv_bounds(face); (b.2, b.3) },
                });
            }
        }
    }
    match <[Side; 2]>::try_from(sides) {
        Ok(sides) => Ok(sides),
        Err(sides) => Err(Error::InvalidInput(format!("a fillet edge must bound two faces, it bounds {}", sides.len()))),
    }
}

/// The point on both faces' surfaces nearest to `guess` in the plane through it square to
/// `tangent` (Newton).
fn onto_both(sides: &[Side; 2], guess: Point, tangent: Vector) -> Result<Point> {
    newton(guess, |x| Ok([sides[0].distance(x)?.0, sides[1].distance(x)?.0, (x - guess).dot(tangent)]))
}

/// The ball of radius `r` touching both faces, centred in the plane through `o` square to
/// `tangent`, on the given side of the faces (-1 inside the solid): from the centre at
/// offset `near` from `o` (a neighbouring section's), else from between the faces' tangent
/// planes.
fn section(sides: &[Side; 2], o: Point, tangent: Vector, r: f64, side: f64, near: Option<Vector>) -> Result<Section> {
    if let Some(offset) = near {
        let start = o + (offset - tangent * offset.dot(tangent));
        let solved = newton(start, |q| Ok([sides[0].distance(q)?.0 - side * r, sides[1].distance(q)?.0 - side * r, (q - o).dot(tangent)]));
        if let Ok(center) = solved {
            let contacts = [sides[0].distance(center)?.1, sides[1].distance(center)?.1];
            let contact_uv = [sides[0].uv(center)?, sides[1].uv(center)?];
            return Ok(Section { t: 0.0, center, contacts, contact_uv, radius: r });
        }
    }
    let normals = [sides[0].distance(o)?.2, sides[1].distance(o)?.2];
    let into: Vec<Vector> = (0..2).map(|i| normals[i].cross(tangent * sides[i].along)).collect();
    let tiny = 1e-9;
    let unit = |v: Vector| if v.length() > tiny { v * (1.0 / v.length()) } else { v };
    let sum = unit(into[0]) + unit(into[1]);
    if into.iter().any(|v| v.length() <= tiny) || sum.length() <= tiny {
        return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
    }
    let into: Vec<Vector> = into.into_iter().map(|v| v.normalized()).collect();
    let half = 0.5 * into[0].dot(into[1]).clamp(-1.0, 1.0).acos();
    let start = o + sum.normalized() * (r / half.sin());
    let center = newton(start, |q| {
        Ok([sides[0].distance(q)?.0 - side * r, sides[1].distance(q)?.0 - side * r, (q - o).dot(tangent)])
    })?;
    let contacts = [sides[0].distance(center)?.1, sides[1].distance(center)?.1];
    let contact_uv = [sides[0].uv(center)?, sides[1].uv(center)?];
    Ok(Section { t: 0.0, center, contacts, contact_uv, radius: r })
}

/// Newton's method on a function of a point with three components, from `start`.
fn newton(start: Point, f: impl Fn(Point) -> Result<[f64; 3]>) -> Result<Point> {
    let scale = 1.0 + (start - Point::ORIGIN).length();
    let h = 1e-7 * scale;
    let mut x = start;
    for _ in 0..60 {
        let fx = f(x)?;
        if fx.iter().all(|v| v.abs() <= 1e-12 * scale) {
            return Ok(x);
        }
        let columns: Vec<[f64; 3]> = [Vector::new(h, 0.0, 0.0), Vector::new(0.0, h, 0.0), Vector::new(0.0, 0.0, h)]
            .iter()
            .map(|&d| {
                let g = f(x + d)?;
                Ok([(g[0] - fx[0]) / h, (g[1] - fx[1]) / h, (g[2] - fx[2]) / h])
            })
            .collect::<Result<_>>()?;
        let step = super::solve3([columns[0], columns[1], columns[2]], fx)
            .ok_or_else(|| Error::InvalidInput("no fillet section here (singular system)".into()))?;
        x = x - Vector::new(step[0], step[1], step[2]);
    }
    if f(x)?.iter().all(|v| v.abs() <= 1e-8 * scale) {
        Ok(x)
    } else {
        Err(Error::InvalidInput("no fillet section here (no convergence)".into()))
    }
}

/// The fillet surface: across (u in [0, 1]) the arc of each section from the contact on the
/// first face to the one on the second, as a rational quadratic; along the edge (v, the
/// spine's parameter) a cubic interpolation of the arcs' homogeneous control points.
fn fillet_surface(sections: &[Section]) -> Result<Surface> {
    let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();
    let mut middle = Vec::with_capacity(sections.len());
    let mut weights = Vec::with_capacity(sections.len());
    for s in sections {
        let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        let w = half.cos();
        // The corner of the arc's control polygon, where its end tangents meet.
        let corner = s.center + (a + b).normalized() * (s.radius / w);
        middle.push(Point::ORIGIN + (corner - Point::ORIGIN) * w);
        weights.push(Point::new(w, 0.0, 0.0));
    }
    let rows = [
        fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[0]).collect::<Vec<_>>(), &ts)?,
        fit::interpolate_with_parameters(&middle, &ts)?,
        fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[1]).collect::<Vec<_>>(), &ts)?,
    ];
    let w = fit::interpolate_with_parameters(&weights, &ts)?;
    let w_poles: Vec<f64> = w.poles().iter().map(|p| p.x).collect();
    if w_poles.iter().any(|&w| w <= 0.0) {
        return Err(Error::InvalidInput("the fillet turns too sharply".into()));
    }
    let poles = vec![
        rows[0].poles().to_vec(),
        rows[1].poles().iter().zip(&w_poles).map(|(p, w)| Point::ORIGIN + (*p - Point::ORIGIN) * (1.0 / w)).collect(),
        rows[2].poles().to_vec(),
    ];
    let n = w_poles.len();
    let u = KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2);
    let v = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
    let surface = BSplineSurface::try_new(poles, u, v)?.with_weights(vec![vec![1.0; n], w_poles, vec![1.0; n]]);
    Ok(Surface::BSpline(surface))
}

/// Where the fillet ends at a vertex of the edge: the face there, the curve along which the
/// fillet surface crosses it (from the contact curve on the first face to the one on the
/// second) and the two edges at the vertex it shortens.
struct EndCut {
    vertex: VertexId,
    face: Side,
    /// The crossing as a 3D curve over [0, 1] and as a pcurve on the fillet surface.
    curve: BSplineCurve,
    on_fillet: BSplineCurve2,
    /// The crossing on the end face, over [0, 1].
    on_face: BSplineCurve2,
    /// The fillet's v where it crosses at u = 0 and u = 1.
    v: [f64; 2],
    /// Per face along the edge: the edge at the vertex it shares with the end face, and the
    /// parameter where the contact curve meets it.
    edges: [(EdgeId, f64); 2],
}

#[allow(clippy::too_many_arguments)]
fn end_cut(
    model: &Model,
    shape: Shape,
    id: EdgeId,
    vertex: Shape,
    sides: &[Side; 2],
    fillet: &Surface,
    ts: &[f64],
    at_end: bool,
    target: f64,
) -> Result<EndCut> {
    let v_id = VertexId::try_from(vertex.node).expect("a vertex");
    // The other edges at the vertex: one on each face along the edge.
    let at_vertex = |e: EdgeId| {
        let (a, b) = model.edge_vertices(Shape::from(e));
        a.node == vertex.node || b.node == vertex.node
    };
    let face_edges = |face: Shape| -> Vec<EdgeId> {
        model.explore(face, NodeKind::Edge).filter_map(|e| EdgeId::try_from(e.node).ok()).filter(|&e| e != id && at_vertex(e)).collect()
    };
    let neighbours: Vec<Vec<EdgeId>> = sides.iter().map(|s| face_edges(s.face)).collect();
    let (Some(&a), Some(&b)) = (neighbours[0].first(), neighbours[1].first()) else {
        return Err(Error::InvalidInput("rolling-ball fillets ending where no other edge continues each face are not supported yet".into()));
    };
    if neighbours[0].len() != 1 || neighbours[1].len() != 1 || a == b {
        return Err(Error::InvalidInput("rolling-ball fillets ending at vertices with more than three edges are not supported yet".into()));
    }
    // The end face: the one using both edges.
    let end_faces: Vec<Shape> = model
        .explore(shape, NodeKind::Face)
        .filter(|&f| f.node != sides[0].face.node && f.node != sides[1].face.node)
        .filter(|&f| {
            let edges: Vec<EdgeId> = model.explore(f, NodeKind::Edge).filter_map(|e| EdgeId::try_from(e.node).ok()).collect();
            edges.contains(&a) && edges.contains(&b)
        })
        .collect();
    let [face] = end_faces[..] else {
        return Err(Error::InvalidInput("rolling-ball fillets ending at vertices with more than three faces are not supported yet".into()));
    };
    let sign = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    let end = Side { face, surface: model.surface(face).expect("a face's surface").clone(), sign, along: 1.0, last: Default::default(), v_range: { let b = model.face_uv_bounds(face); (b.2, b.3) } };

    // Across the fillet, where it leaves the end face's surface: going from the middle of
    // the edge towards this end, the outward distance turns positive.
    let mid = 0.5 * (ts[0] + ts[ts.len() - 1]);
    let order: Vec<f64> = if at_end {
        ts.iter().copied().filter(|&t| t >= mid).collect()
    } else {
        ts.iter().rev().copied().filter(|&t| t <= mid).collect()
    };
    // (Where the fillet crosses the end face's surface on its way to this end, either way
    // across it: the end face may lie on either side of the fillet. Of several, the one
    // nearest the edge's vertex.)
    let vertex_point = model.point(vertex);
    let crossing = |u: f64| -> Result<(f64, f64)> {
        let d = |v: f64| -> Result<f64> { Ok(end.distance(fillet.point(u, v))?.0) };
        // (Into the end face's outside first, as the fillet usually leaves the solid there.)
        let mut best: Option<(f64, f64)> = None;
        for outwards in [true, false] {
        if best.is_some() {
            break;
        }
        for w in order.windows(2) {
            let (da, db) = (d(w[0])?, d(w[1])?);
            if (da <= 0.0) != (db <= 0.0) && (da <= 0.0) == outwards {
                let (mut lo, mut hi) = (w[0], w[1]);
                for _ in 0..80 {
                    let m = 0.5 * (lo + hi);
                    if (d(m)? <= 0.0) == (da <= 0.0) {
                        lo = m;
                    } else {
                        hi = m;
                    }
                }
                let v = 0.5 * (lo + hi);
                // (Outwards the first from the middle; the other way, the one nearest the
                // vertex.)
                if outwards {
                    best = Some((0.0, v));
                    break;
                }
                let gap = fillet.point(u, v).distance(vertex_point);
                if best.is_none_or(|(g, _)| gap < g) {
                    best = Some((gap, v));
                }
            }
        }
        }
        best.map(|(_, v)| (u, v)).ok_or_else(|| Error::InvalidInput("the fillet does not reach the face at the end of its edge".into()))
    };
    // Where the contact curves meet the edges at the vertex.
    // (Where a contact curve would meet the edge past its far end, the fillet is wider than
    // the face there: it would run over onto the next face, which is not built here. Told
    // from the crossing's ends before it is refined.)
    let meet = |e: EdgeId, p: Point| -> Result<f64> {
        let (c, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        let t = parameter_on(&c, first, last, p);
        if c.point(t).distance(p) > (1e3 * target).max(1e-6 * (1.0 + (p - Point::ORIGIN).length())) {
            return Err(Error::InvalidInput("a fillet running over onto the next face is not supported yet".into()));
        }
        Ok(t)
    };
    let ends = [crossing(0.0)?, crossing(1.0)?];
    let edges = [(a, meet(a, fillet.point(ends[0].0, ends[0].1))?), (b, meet(b, fillet.point(ends[1].0, ends[1].1))?)];
    // More samples across until the interpolated crossing stays on the end face (or doubling
    // them hardly helps).
    let mut across = ACROSS;
    let mut previous_worst = f64::INFINITY;
    let (uv, curve, on_fillet) = loop {
        let uv: Vec<(f64, f64)> = (0..=across).map(|j| crossing(j as f64 / across as f64)).collect::<Result<_>>()?;
        let params: Vec<f64> = uv.iter().map(|&(u, _)| u).collect();
        let points: Vec<Point> = uv.iter().map(|&(u, v)| fillet.point(u, v)).collect();
        let curve = fit::interpolate_with_parameters(&points, &params)?;
        let on_fillet = fit::interpolate_with_parameters(&uv.iter().map(|&(u, v)| Point2::new(u, v)).collect::<Vec<_>>(), &params)?;
        let mut worst: f64 = 0.0;
        for w in params.windows(2) {
            let u = 0.5 * (w[0] + w[1]);
            let p = curve.point(u);
            let q = on_fillet.point(u);
            worst = worst.max(end.distance(p)?.0.abs()).max(fillet.point(q.x, q.y).distance(p));
        }
        let stalled = across >= 4 * ACROSS && worst > 0.25 * previous_worst;
        if worst <= target || stalled || across >= 64 * ACROSS {
            break (uv, curve, on_fillet);
        }
        previous_worst = worst;
        across *= 2;
    };
    let params: Vec<f64> = uv.iter().map(|&(u, _)| u).collect();
    let face_uv: Vec<Point2> = uv.iter().map(|&(u, v)| end.uv(fillet.point(u, v))).collect::<Result<_>>()?;
    let on_face = fit::interpolate_with_parameters(&unwrapped(face_uv, &end.surface), &params)?;
    Ok(EndCut { vertex: v_id, face: end, curve, on_fillet, on_face, v: [uv[0].1, uv[across].1], edges })
}

/// The parameter of the point of `curve` over `[first, last]` nearest to `p`.
fn parameter_on(curve: &Curve, first: f64, last: f64, p: Point) -> f64 {
    let n = 200;
    let gap = |t: f64| curve.point(t).distance(p);
    let mut best = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).map(|t| (t, gap(t))).min_by(|a, b| a.1.total_cmp(&b.1)).expect("samples").0;
    let step = (last - first) / n as f64;
    let (mut lo, mut hi) = ((best - step).max(first.min(last)), (best + step).min(first.max(last)));
    let ratio = 0.618_033_988_749_895;
    for _ in 0..100 {
        let (x1, x2) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
        if gap(x1) < gap(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    best = 0.5 * (lo + hi);
    best
}

/// Rebuilds the solid with the fillet.
struct Rebuild<'a> {
    model: &'a mut Model,
    id: EdgeId,
    sides: &'a [Side; 2],
    surface: &'a Surface,
    contacts: &'a [BSplineCurve; 2],
    ends: &'a [EndCut],
    ts: &'a [f64],
}

/// A coedge to build: edge, orientation, pcurve and range.
type Use = (EdgeId, Orientation, Curve2, f64, f64);

impl Rebuild<'_> {
    fn build(&mut self, shape: Shape, closed: bool, sections: &[Section], convex: bool) -> Result<Shape> {
        // Within twice the contact curves' distance from their faces, and no finer than
        // the model's size allows pcurves to be fitted (1e-8 of it).
        let size = {
            let (a, b) = (self.ts[0], *self.ts.last().expect("sections"));
            let (p, q) = (self.contacts[0].point(a), self.contacts[0].point(b));
            p.distance(q) + (p - Point::ORIGIN).length() + 1.0
        };
        let tolerance = (self.deviation()? * 2.0).max(crate::tolerance::LINEAR).max(1e-8 * size);
        // The fillet's v range along each contact curve.
        let (v_start, v_end) = match self.ends {
            [start, end] => ([start.v[0], start.v[1]], [end.v[0], end.v[1]]),
            _ => ([self.ts[0]; 2], [*self.ts.last().expect("sections"); 2]),
        };
        // Vertices at the ends of the contact curves.
        let corner = |model: &mut Model, i: usize, v: f64| model.add_vertex(self.contacts[i].point(v), tolerance);
        let starts = [corner(self.model, 0, v_start[0]), corner(self.model, 1, v_start[1])];
        let finishes = if closed { starts } else { [corner(self.model, 0, v_end[0]), corner(self.model, 1, v_end[1])] };
        // Along a closed edge, where each contact curve crosses the other edges of its face
        // at the edge's vertex (a seam): its vertex there.
        let vertex_id = VertexId::try_from(self.model.edge_vertices(Shape::from(self.id)).0.node).expect("a vertex");
        let at_vertex: [Vec<EdgeId>; 2] = [0, 1].map(|i| {
            self.model
                .explore(self.sides[i].face, NodeKind::Edge)
                .filter_map(|e| EdgeId::try_from(e.node).ok())
                .filter(|&e| e != self.id)
                .filter(|&e| {
                    let (a, b) = self.model.edge_vertices(Shape::from(e));
                    let v: crate::topology::NodeId = vertex_id.into();
                    a.node == v || b.node == v
                })
                .collect()
        });
        let mut splits: [Option<(f64, VertexId)>; 2] = [None, None];
        if closed {
            for i in 0..2 {
                let Some(&other) = at_vertex[i].first() else { continue };
                let Some((curve, first, last)) = self.model.edge_curve(Shape::from(other)) else { continue };
                let gap = |v: f64| crate::geometry::extrema::nearest_on_curve(self.contacts[i].point(v), &curve, first, last).distance;
                let n = 400;
                let (mut best, mut at) = (f64::INFINITY, v_start[i]);
                for k in 0..=n {
                    let v = v_start[i] + (v_end[i] - v_start[i]) * k as f64 / n as f64;
                    let d = gap(v);
                    if d < best {
                        (best, at) = (d, v);
                    }
                }
                let step = (v_end[i] - v_start[i]) / n as f64;
                let (mut lo, mut hi) = (at - step, at + step);
                for _ in 0..80 {
                    let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
                    if gap(a) < gap(b) { hi = b } else { lo = a }
                }
                let v = 0.5 * (lo + hi);
                let span = (v_end[i] - v_start[i]).abs();
                if v > v_start[i] + 1e-6 * span && v < v_end[i] - 1e-6 * span {
                    splits[i] = Some((v, self.model.add_vertex(self.contacts[i].point(v), tolerance)));
                }
            }
        }
        // The contact curves as edges (in pieces where split), in order along them.
        let mut contact_edges: [Vec<(EdgeId, f64, f64)>; 2] = [Vec::new(), Vec::new()];
        for i in 0..2 {
            let pieces = match splits[i] {
                Some((v, x)) => vec![(starts[i], x, v_start[i], v), (x, finishes[i], v, v_end[i])],
                None => vec![(starts[i], finishes[i], v_start[i], v_end[i])],
            };
            for (a, b, lo, hi) in pieces {
                let e = self.model.add_edge(a, b, tolerance)?;
                self.model.set_edge_curve(e.into(), Curve::BSpline(self.contacts[i].clone()), tolerance)?;
                self.model.set_edge_range(e.into(), lo, hi)?;
                contact_edges[i].push((e, lo, hi));
            }
        }
        // The fillet's boundary across at each end: the end cut, or for a closed edge the arc
        // of the first section (a seam).
        let mut across = Vec::new();
        let mut shortened: Vec<(EdgeId, EdgeId, f64, f64)> = Vec::new();
        // (edge, its vertex at a fillet's end, the vertex replacing it, the parameter there)
        let mut trims: Vec<(EdgeId, VertexId, VertexId, f64)> = Vec::new();
        if closed {
            let s = &sections[0];
            let arc = arc_curve(s);
            let e = self.model.add_edge(starts[0], starts[1], tolerance)?;
            self.model.set_edge_curve(e.into(), arc, tolerance)?;
            self.model.set_edge_range(e.into(), 0.0, 1.0)?;
            across.push(e);
            // Other edges of the faces at the edge's vertex (a seam) now end on the contact
            // curves.
            let vertex = self.model.edge_vertices(Shape::from(self.id)).0;
            let v_id = VertexId::try_from(vertex.node).expect("a vertex");
            for i in 0..2 {
                let others: Vec<EdgeId> = self
                    .model
                    .explore(self.sides[i].face, NodeKind::Edge)
                    .filter_map(|e| EdgeId::try_from(e.node).ok())
                    .filter(|&e| e != self.id && !shortened.iter().any(|s| s.0 == e))
                    .collect();
                for old in others {
                    let (a, b) = self.model.edge_vertices(Shape::from(old));
                    let (a_id, b_id) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
                    if a_id != v_id && b_id != v_id {
                        continue;
                    }
                    let (curve, first, last) = self.model.edge_curve(Shape::from(old)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet".into()))?;
                    // Ending where the contact curve crosses it.
                    let (meet, at) = match splits[i] {
                        Some((v, x)) => (x, self.contacts[i].point(v)),
                        None => (starts[i], self.contacts[i].point(v_start[i])),
                    };
                    let t = parameter_on(&curve, first, last, at);
                    let (from, to, lo, hi) = if a_id == v_id { (meet, b_id, t, last) } else { (a_id, meet, first, t) };
                    let e = self.model.add_edge(from, to, tolerance.max(self.model.tolerance(Shape::from(old))))?;
                    self.model.set_edge_curve(e.into(), curve, tolerance)?;
                    self.model.set_edge_range(e.into(), lo, hi)?;
                    shortened.push((old, e, lo, hi));
                }
            }
        } else {
            for (k, cut) in self.ends.iter().enumerate() {
                let vs = if k == 0 { starts } else { finishes };
                let e = self.model.add_edge(vs[0], vs[1], tolerance)?;
                self.model.set_edge_curve(e.into(), Curve::BSpline(cut.curve.clone()), tolerance)?;
                self.model.set_edge_range(e.into(), 0.0, 1.0)?;
                across.push(e);
                // The edges at the vertex now end where the contact curves meet them.
                for (i, &(old, t)) in cut.edges.iter().enumerate() {
                    trims.push((old, cut.vertex, vs[i], t));
                }
            }
            // Each trimmed edge once, at one end or both.
            let mut done: Vec<EdgeId> = Vec::new();
            for &(old, ..) in &trims {
                if done.contains(&old) {
                    continue;
                }
                done.push(old);
                let (a, b) = self.model.edge_vertices(Shape::from(old));
                let (curve, first, last) = self.model.edge_curve(Shape::from(old)).expect("a curve");
                let (mut from, mut to) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
                let (mut lo, mut hi) = (first, last);
                let (a_id, b_id) = (from, to);
                for &(_, at, new_vertex, t) in trims.iter().filter(|c| c.0 == old) {
                    if at == a_id {
                        (from, lo) = (new_vertex, t);
                    } else if at == b_id {
                        (to, hi) = (new_vertex, t);
                    }
                }
                let e = self.model.add_edge(from, to, tolerance.max(self.model.tolerance(Shape::from(old))))?;
                self.model.set_edge_curve(e.into(), curve, tolerance)?;
                self.model.set_edge_range(e.into(), lo, hi)?;
                shortened.push((old, e, lo, hi));
            }
        }

        // The faces along the edge, with the contact curves for the edge, and the end faces,
        // with the end cuts.
        let mut rebuilt: Vec<(Shape, Shape)> = Vec::new();
        for (i, side) in self.sides.iter().enumerate() {
            let pcurve = match project::analytic_pcurve(&Curve::BSpline(self.contacts[i].clone()), &side.surface) {
                Some(p) => p,
                None => Curve2::BSpline(fit::interpolate_with_parameters(
                    &unwrapped(sections.iter().map(|s| s.contact_uv[i]).collect(), &side.surface),
                    self.ts,
                )?),
            };
            let face = rebuild_face(self.model, side.face, &[(self.id, &contact_edges[i], pcurve)], &shortened, &[])?;
            rebuilt.push((side.face, face));
        }
        // Each end face once, with the cuts of both ends where it is at both.
        let mut end_faces: Vec<(Shape, Vec<(EdgeId, Curve2, [VertexId; 2])>)> = Vec::new();
        for (k, cut) in self.ends.iter().enumerate() {
            let pcurve = match project::analytic_pcurve(&Curve::BSpline(cut.curve.clone()), &cut.face.surface) {
                Some(p) => p,
                None => Curve2::BSpline(cut.on_face.clone()),
            };
            let vs = if k == 0 { starts } else { finishes };
            let insert = (across[k], pcurve, vs);
            match end_faces.iter_mut().find(|(f, _)| f.node == cut.face.face.node) {
                Some((_, inserts)) => inserts.push(insert),
                None => end_faces.push((cut.face.face, vec![insert])),
            }
        }
        for (face, inserts) in end_faces {
            if self.sides.iter().any(|s| s.face.node == face.node) {
                todo!("rolling-ball fillets ending on one of their own faces");
            }
            let new = rebuild_face(self.model, face, &[], &shortened, &inserts)?;
            rebuilt.push((face, new));
        }

        // The fillet face: across at the start (u from 0 to 1), along the second contact,
        // back across at the end, back along the first contact: counter-clockwise in (u, v).
        let line = |o: (f64, f64), d: (f64, f64)| Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1)));
        let across_pcurve = |k: usize| -> Curve2 {
            if closed {
                line((0.0, if k == 0 { v_start[0] } else { v_end[0] }), (1.0, 0.0))
            } else {
                Curve2::BSpline(self.ends[k].on_fillet.clone())
            }
        };
        let last_across = if closed { across[0] } else { across[1] };
        let mut uses: Vec<Use> = vec![(across[0], Orientation::Forward, across_pcurve(0), 0.0, 1.0)];
        for &(e, lo, hi) in &contact_edges[1] {
            uses.push((e, Orientation::Forward, line((1.0, 0.0), (0.0, 1.0)), lo, hi));
        }
        uses.push((last_across, Orientation::Reversed, across_pcurve(1), 0.0, 1.0));
        for &(e, lo, hi) in contact_edges[0].iter().rev() {
            uses.push((e, Orientation::Reversed, line((0.0, 0.0), (0.0, 1.0)), lo, hi));
        }
        let fillet = new_face(self.model, self.surface.clone(), vec![(uses, Orientation::Forward)], tolerance)?;
        // Outward: away from the ball's centre for a convex edge (the ball is in the solid).
        let mid = &sections[sections.len() / 2];
        let (p, du, dv) = self.surface.d1(0.5, mid.t);
        let away = (p - mid.center).dot(du.cross(dv)) > 0.0;
        let fillet = if away == convex { fillet } else { fillet.reversed() };

        // The solid with the faces replaced and the fillet added.
        let mut solids = Vec::new();
        for solid in self.model.explore(shape, NodeKind::Solid).collect::<Vec<_>>() {
            let faces: Vec<Shape> = self.model.explore(solid, NodeKind::Face).collect();
            let touched = faces.iter().any(|f| rebuilt.iter().any(|(old, _)| old.node == f.node));
            if !touched {
                solids.push(solid);
                continue;
            }
            let shell = self.model.new_empty(NodeKind::Shell);
            for f in faces {
                match rebuilt.iter().find(|(old, _)| old.node == f.node) {
                    Some(&(_, new)) => {
                        self.model.add(shell, Shape { orientation: f.orientation, ..new })?;
                        if new.node == rebuilt[1].1.node {
                            self.model.add(shell, fillet)?;
                        }
                    }
                    None => self.model.add(shell, f)?,
                }
            }
            let new_solid = self.model.new_empty(NodeKind::Solid);
            self.model.add(new_solid, shell)?;
            solids.push(new_solid);
        }
        let compound = self.model.new_empty(NodeKind::Compound);
        for s in solids {
            self.model.add(compound, s)?;
        }
        Ok(compound)
    }

    /// The largest distance of the contact curves from their faces' surfaces.
    fn deviation(&self) -> Result<f64> {
        let mut worst: f64 = 0.0;
        for i in 0..2 {
            let (a, b) = (self.ts[0], *self.ts.last().expect("sections"));
            for k in 0..=200 {
                let p = self.contacts[i].point(a + (b - a) * k as f64 / 200.0);
                worst = worst.max(self.sides[i].distance(p)?.0.abs());
            }
        }
        Ok(worst)
    }

}

/// A copy of `face` (same surface, forward) whose boundary has `replace`'s edge for the
/// filleted one, the shortened edges for the old ones, and `insert`'s edge between the
/// two shortened edges (which end at `insert`'s vertices).
fn rebuild_face(
    model: &mut Model,
    face: Shape,
    replaces: &[(EdgeId, &[(EdgeId, f64, f64)], Curve2)],
    shortened: &[(EdgeId, EdgeId, f64, f64)],
    inserts: &[(EdgeId, Curve2, [VertexId; 2])],
) -> Result<Shape> {
    let forward = Shape { orientation: Orientation::Forward, ..face };
    let surface = model.surface(forward).expect("a face's surface").clone();
    let tolerance = model.tolerance(forward);
    let mut wires = Vec::new();
    for wire in model.children(forward).collect::<Vec<_>>() {
        let mut uses: Vec<Use> = Vec::new();
        let mut coedges: Vec<Shape> = model.children_local(wire).collect();
        // (Starting at a coedge kept, for the replacing pieces to follow one.)
        let replaced = |model: &Model, c: &Shape| replaces.iter().any(|r| r.0 == model[CoEdgeId::try_from(c.node).expect("a wire holds coedges")].edge);
        if let Some(k) = coedges.iter().position(|c| !replaced(model, c)) {
            coedges.rotate_left(k);
        }
        for coedge in coedges {
            let c = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
            let (edge, orientation) = (model[c].edge, model[c].orientation);
            let (pcurve, first, last) = model.effective_pcurve(c).ok_or_else(|| Error::InvalidInput("a face without pcurves".into()))?;
            if let Some((_, pieces, p)) = replaces.iter().find(|r| r.0 == edge) {
                {
                    // In the same periods as the edge it replaces; its pieces in the order
                    // the wire runs.
                    let (a, b) = (pieces[0].1, pieces[pieces.len() - 1].2);
                    let target = pcurve.point(0.5 * (first + last));
                    let shifted = align(p, 0.5 * (a + b), target, &surface);
                    let mut ordered: Vec<&(EdgeId, f64, f64)> = if orientation == Orientation::Reversed { pieces.iter().rev().collect() } else { pieces.iter().collect() };
                    // (A closed edge's pieces go round from where the coedge before them
                    // ends: a seam may reach them at any of their vertices.)
                    if let Some(previous) = uses.last() {
                        let (s, e) = model.edge_vertices(Shape::from(previous.0));
                        let arrived = if previous.1 == Orientation::Reversed { s.node } else { e.node };
                        let start_of = |model: &Model, edge: EdgeId| {
                            let (s, e) = model.edge_vertices(Shape::from(edge));
                            if orientation == Orientation::Reversed { e.node } else { s.node }
                        };
                        if let Some(k) = ordered.iter().position(|p| start_of(model, p.0) == arrived) {
                            ordered.rotate_left(k);
                        }
                    }
                    for &&(new, lo, hi) in &ordered {
                        // Each piece in the periods where the coedge before it ends (a
                        // closed edge's pieces lie on either side of a seam).
                        let curve = match uses.last() {
                            Some(previous) if pieces.len() > 1 => {
                                let target = previous.2.point(if previous.1 == Orientation::Reversed { previous.3 } else { previous.4 });
                                let from = if orientation == Orientation::Reversed { hi } else { lo };
                                align(p, from, target, &surface)
                            }
                            _ => shifted.clone(),
                        };
                        uses.push((new, orientation, curve, lo, hi));
                    }
                    continue;
                }
            }
            if let Some(&(_, new, a, b)) = shortened.iter().find(|s| s.0 == edge) {
                uses.push((new, orientation, pcurve, a, b));
                continue;
            }
            uses.push((edge, orientation, pcurve, first, last));
        }
        for &(edge, ref pcurve, [a, b]) in inserts {
            // Between the coedge ending at one of the edge's vertices and the one starting
            // at the other.
            let n = uses.len();
            let ends = |model: &Model, u: &Use| {
                let (s, e) = model.edge_vertices(Shape::from(u.0));
                let (s, e) = (VertexId::try_from(s.node).expect("vertex"), VertexId::try_from(e.node).expect("vertex"));
                if u.1 == Orientation::Reversed { (e, s) } else { (s, e) }
            };
            if let Some(k) = (0..n).find(|&k| {
                let (_, end) = ends(model, &uses[k]);
                let (start, _) = ends(model, &uses[(k + 1) % n]);
                (end == a && start == b) || (end == b && start == a)
            }) {
                let (_, end) = ends(model, &uses[k]);
                let orientation = if end == a { Orientation::Forward } else { Orientation::Reversed };
                // Starting where the previous coedge ends, in the surface's periods.
                let previous = &uses[k];
                let target = previous.2.point(if previous.1 == Orientation::Reversed { previous.3 } else { previous.4 });
                let from = if orientation == Orientation::Reversed { 1.0 } else { 0.0 };
                let shifted = align(pcurve, from, target, &surface);
                uses.insert(k + 1, (edge, orientation, shifted, 0.0, 1.0));
            }
        }
        wires.push((uses, wire.orientation));
    }
    new_face(model, surface, wires, tolerance)
}

fn new_face(model: &mut Model, surface: Surface, wires: Vec<(Vec<Use>, Orientation)>, tolerance: f64) -> Result<Shape> {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, tolerance)?;
    for (uses, orientation) in wires {
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, o, pcurve, first, last) in uses {
            let c = model.add_coedge(edge, o)?;
            model.set_pcurve_with_range(c, pcurve, first, last)?;
            coedges.push(c);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, Shape { orientation, ..Shape::from(wire) })?;
    }
    Ok(face)
}

/// Parameters along a curve made continuous across the surface's periods.
fn unwrapped(mut uv: Vec<Point2>, surface: &Surface) -> Vec<Point2> {
    for i in 1..uv.len() {
        if let Some(period) = surface.u_period() {
            uv[i].x += ((uv[i - 1].x - uv[i].x) / period).round() * period;
        }
        if let Some(period) = surface.v_period() {
            uv[i].y += ((uv[i - 1].y - uv[i].y) / period).round() * period;
        }
    }
    uv
}

/// `pcurve` moved by whole periods of the surface's parameters so that its point at `t`
/// is nearest to `target`.
fn align(pcurve: &Curve2, t: f64, target: Point2, surface: &Surface) -> Curve2 {
    use crate::math::{Transform2, Transformable2, Vector2};
    let p = pcurve.point(t);
    let shift = |period: Option<f64>, d: f64| period.map_or(0.0, |period| (d / period).round() * period);
    let (du, dv) = (shift(surface.u_period(), target.x - p.x), shift(surface.v_period(), target.y - p.y));
    if du == 0.0 && dv == 0.0 {
        return pcurve.clone();
    }
    pcurve.transformed(&Transform2::translation(Vector2::new(du, dv)))
}

/// The arc of a section, from the first contact to the second, over [0, 1] (a rational
/// quadratic B-spline).
fn arc_curve(s: &Section) -> Curve {
    let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
    let w = (0.5 * a.cross(b).length().atan2(a.dot(b))).cos();
    let corner = s.center + (a + b).normalized() * (s.radius / w);
    let curve = BSplineCurve::new(vec![s.contacts[0], corner, s.contacts[1]], vec![0.0, 1.0], vec![3, 3], 2).with_weights(vec![1.0, w, 1.0]);
    Curve::BSpline(curve)
}
