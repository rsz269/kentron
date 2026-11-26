//! Sweeps known in closed form: tools whose swept solid is built face by face, without
//! following where the motion grazes them or resolving a product that runs into itself.
//!
//! # A ball
//!
//! A ball of radius `r` covers what lies within `r` of its centre, whichever way it is turned
//! about the centre; so under any rigid motion it sweeps the points within `r` of its centre's
//! path: a tube about the path (where the motion grazes the ball, the great circle square to
//! its centre's velocity, swept) and a half ball at each end (behind the motion at the start,
//! ahead of it at the end). Where no two pieces of that boundary meet (the tube does not run
//! into itself, nor into the end balls), that is the swept solid's whole boundary:
//!
//! * slid: a cylinder between two half spheres (a capsule);
//! * turned about an axis farther from the centre than `r`: a piece of the torus about the axis
//!   through the centre, its ends the half spheres (taken where the turn and the reach of the
//!   two end balls round the axis, `2·asin(r / R)` at a distance `R`, leave a gap; a whole turn
//!   or more, the whole torus);
//! * screwed: the tube about the centre's helix, which the screw carries along itself: the
//!   circle square to the helix at the start, screwed (B-spline surfaces through its poles'
//!   helices, in quarter turns, within about 1e-10 of the size, as [`super::product`] sweeps a
//!   curve), taken where the tube is embedded: the helix's curvature radius `(R² + c²) / R`
//!   (`c` the rise per radian) beyond `r`, and centres on it that are `2r` apart staying at
//!   least that far apart as they part further (no turn of the tube running into the next).
//!   A screw about an axis through the centre is a slide along it; one that does not rise, a
//!   turn.
//!
//! # A ball-end mill
//!
//! A half ball of radius `r` under a cylinder of the same radius about the axis through its
//! centre, flat at a height `h` above the centre, is convex, so slid it fills the convex hull of
//! its two places. Seen along the slide `d`, its outline is the half great circle square to `d`
//! below its equator, the two lines of the cylinder square to `d`'s part across the axis, and
//! the half of the top's rim the slide carries in under the top (none where `d` lies in the
//! top's plane: the whole top is then swept along itself). The boundary is that outline slid
//! (a half cylinder about the centre's path, two planes, and a cylinder of lines along `d`
//! through the rim's half circle), the mill behind the slide at the start and ahead of it at
//! the end (each a lune of the half ball, between its equator and the outline's half circle,
//! and a half of the cylinder), and the top where the slide leaves it. Slid along its axis, the
//! mill is the mill lengthened.
//!
//! Anything else (other tools, other motions, the conditions failing) is left to the general
//! sweep.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Line2};
use crate::geometry::surface::{Cylinder, ExtrusionSurface, KnotVector, Plane, Sphere, Surface, Torus};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Transformable, Vector};
use crate::sweep::pipe::path;
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{product, RigidMotion, Steady};

/// A tool whose sweep is known in closed form.
#[derive(Debug, Clone, Copy)]
pub(super) enum Tool {
    /// The ball about `centre`.
    Ball { centre: Point, radius: f64 },
    /// A ball-end mill: the half of the ball about `centre` below its equator (square to
    /// `axis`), under the cylinder about the axis from the centre up to `height` above it,
    /// flat there.
    BallEnd { centre: Point, axis: Direction, radius: f64, height: f64 },
}

/// The tool `solid` is, if it is one of [`Tool`]'s: told from its faces' surfaces (all on one
/// sphere; or on one sphere, one cylinder of its radius about an axis through its centre, and
/// one plane square to that axis), its vertices' reach along the axis and its volume.
pub(super) fn recognised(model: &Model, solid: Shape) -> Option<Tool> {
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let (mut sphere, mut cylinder, mut plane): (Option<(Point, f64)>, Option<(Axis, f64)>, Option<Plane>) = (None, None, None);
    let mut size = 0.0f64;
    let surfaces: Vec<Surface> = faces
        .iter()
        .map(|&f| {
            let s = model.surface(f)?;
            Some(if f.placement == Transform::IDENTITY { s.clone() } else { s.transformed(&f.placement) })
        })
        .collect::<Option<_>>()?;
    for s in &surfaces {
        match s {
            Surface::Sphere(b) => size = size.max(b.radius()),
            Surface::Cylinder(c) => size = size.max(c.radius()),
            _ => {}
        }
    }
    if !(size > 0.0) {
        return None;
    }
    let tol = 1e-9 * size;
    for s in &surfaces {
        match s {
            Surface::Sphere(b) => {
                let (c, r) = (b.frame().origin(), b.radius());
                match sphere {
                    None => sphere = Some((c, r)),
                    Some((c0, r0)) if c0.distance(c) <= tol && (r0 - r).abs() <= tol => {}
                    _ => return None,
                }
            }
            Surface::Cylinder(cy) => {
                let (axis, r) = (Axis::new(cy.frame().origin(), cy.frame().z()), cy.radius());
                match cylinder {
                    None => cylinder = Some((axis, r)),
                    Some((a0, r0)) if (r0 - r).abs() <= tol && Vector::from(a0.direction).cross(Vector::from(axis.direction)).length() <= 1e-12 && off_axis(a0, axis.origin) <= tol => {}
                    _ => return None,
                }
            }
            Surface::Plane(p) => match plane {
                None => plane = Some(*p),
                Some(p0) if Vector::from(p0.normal()).cross(Vector::from(p.normal())).length() <= 1e-12 && (p.origin() - p0.origin()).dot(Vector::from(p0.normal())).abs() <= tol => {}
                _ => return None,
            },
            _ => return None,
        }
    }
    let (centre, radius) = sphere?;
    let volume = crate::properties::volume(model, solid).mass;
    let close = |v: f64, expected: f64| (v - expected).abs() <= 1e-6 * expected;
    match (cylinder, plane) {
        (None, None) => close(volume, 4.0 / 3.0 * PI * radius.powi(3)).then_some(Tool::Ball { centre, radius }),
        (Some((axis, r)), Some(top)) => {
            if (r - radius).abs() > tol || off_axis(axis, centre) > tol || Vector::from(top.normal()).cross(Vector::from(axis.direction)).length() > 1e-12 {
                return None;
            }
            let along = (top.origin() - centre).dot(Vector::from(axis.direction));
            let up = if along > 0.0 { axis.direction } else { -axis.direction };
            let height = along.abs();
            if !(height > tol) {
                return None;
            }
            let u = Vector::from(up);
            let reach_ok = model.unique_subshapes(solid, NodeKind::Vertex).iter().all(|v| {
                let p = model[VertexId::from_index(v.node.index)].point.transformed(&v.placement);
                let h = (p - centre).dot(u);
                h >= -radius - tol && h <= height + tol
            });
            (reach_ok && close(volume, 2.0 / 3.0 * PI * radius.powi(3) + PI * radius * radius * height)).then_some(Tool::BallEnd { centre, axis: up, radius, height })
        }
        _ => None,
    }
}

/// How far `p` lies from `axis`.
fn off_axis(axis: Axis, p: Point) -> f64 {
    (p - axis.origin).cross(Vector::from(axis.direction)).length()
}

/// The solid `tool` sweeps under `motion`, built face by face; `None` where the motion is not
/// one taken for the tool (see the module's notes), for the general sweep.
pub(super) fn sweep(model: &mut Model, tool: Tool, motion: &RigidMotion) -> Option<Result<Shape>> {
    match tool {
        Tool::Ball { centre, radius } => ball(model, centre, radius, motion),
        Tool::BallEnd { centre, axis, radius, height } => match *motion {
            RigidMotion::Slide { by } if by.length() > 0.0 && by.length().is_finite() => Some(ball_end_slid(model, centre, Vector::from(axis), radius, height, by)),
            _ => None,
        },
    }
}

/// A ball's sweep (see the module's notes).
fn ball(model: &mut Model, centre: Point, r: f64, motion: &RigidMotion) -> Option<Result<Shape>> {
    if let RigidMotion::Carried(_) = motion {
        return None;
    }
    let (start, steady) = motion.steady().ok()?;
    let c0 = centre.transformed(&start);
    match steady {
        Steady::Slide { by } => Some(capsule(model, c0, by, r)),
        Steady::Turn { axis, angle } => turned_ball(model, c0, r, axis, angle),
        Steady::Screw { axis, angle, rise } => {
            if rise == 0.0 {
                return turned_ball(model, c0, r, axis, angle);
            }
            let big = off_axis(axis, c0);
            if big <= 1e-12 * (r + rise.abs()) {
                // (Its centre on the axis: slid along it.)
                return Some(capsule(model, c0, Vector::from(axis.direction) * rise, r));
            }
            helical_embedded(big, rise / angle, angle, r).then(|| screwed_ball(model, c0, r, axis, angle, rise))
        }
        Steady::Carried(_) => None,
    }
}

/// A ball at `c` turned by `angle` (positive) about `axis`: the ball itself where its centre
/// is on the axis, else a piece of torus with half balls at its ends (see the module's notes);
/// `None` where the axis meets the ball or the ends' half balls would meet.
fn turned_ball(model: &mut Model, c: Point, r: f64, axis: Axis, angle: f64) -> Option<Result<Shape>> {
    let big = off_axis(axis, c);
    if big <= 1e-12 * r {
        return Some(whole_ball(model, c, r));
    }
    if !(big > r * (1.0 + 1e-9)) {
        return None;
    }
    if angle < TAU && angle + 2.0 * (r / big).asin() > TAU * (1.0 - 1e-12) {
        return None;
    }
    Some(torus_piece(model, c, r, axis, angle.min(TAU)))
}

/// True where the tube of radius `r` about the helix of radius `big` rising `c` a radian, over
/// `angle`, is embedded (see the module's notes).
fn helical_embedded(big: f64, c: f64, angle: f64, r: f64) -> bool {
    let curvature = big / (big * big + c * c);
    if !(r * curvature < 1.0 - 1e-9) {
        return false;
    }
    // The squared distance between centres `delta` apart in angle.
    let apart = |delta: f64| 2.0 * big * big * (1.0 - delta.cos()) + c * c * delta * delta;
    let reach = 4.0 * r * r;
    // Where centres first part by 2r (the distance grows up to half a turn): none within the
    // motion or half a turn, a tube short enough not to reach itself, or one whose turns
    // overlap.
    let near = angle.min(PI);
    if apart(near) < reach {
        return angle <= PI;
    }
    let (mut lo, mut hi) = (0.0, near);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if apart(mid) < reach { lo = mid } else { hi = mid }
    }
    let first = hi;
    // Beyond it, never back within 2r: sampled finely, each sampled minimum refined.
    let n = ((angle - first) / (PI / 64.0)).ceil().max(1.0) as usize;
    let at = |i: usize| first + (angle - first) * i as f64 / n as f64;
    for i in 0..=n {
        let d = at(i);
        if apart(d) <= reach * (1.0 + 1e-9) && i > 0 {
            return false;
        }
        if i > 0 && i < n && apart(d) <= apart(at(i - 1)) && apart(d) <= apart(at(i + 1)) {
            let (mut a, mut b) = (at(i - 1), at(i + 1));
            for _ in 0..100 {
                let (m1, m2) = (a + (b - a) / 3.0, b - (b - a) / 3.0);
                if apart(m1) < apart(m2) { b = m2 } else { a = m1 }
            }
            if apart(0.5 * (a + b)) <= reach * (1.0 + 1e-9) {
                return false;
            }
        }
    }
    true
}

/// A side of a face on the rectangle `[u1, u2] × [v1, v2]` of its surface, counter-clockwise
/// from the side at `v1`.
enum Side {
    /// An edge running along the side (an iso-curve), its pcurve the side's line.
    Edge(EdgeId),
    /// A side collapsed to a point (a pole): a degenerate edge of the face's own.
    Pole,
    /// An edge used as given: its orientation and pcurve over its range.
    Given(EdgeId, Orientation, Curve2, f64, f64),
}

/// The faces of a swept solid made one by one, sharing vertices and edges.
struct Build<'m> {
    model: &'m mut Model,
    tol: f64,
    vertices: Vec<VertexId>,
    shell: Shape,
}

impl<'m> Build<'m> {
    fn new(model: &'m mut Model, size: f64) -> Build<'m> {
        let shell = model.new_empty(NodeKind::Shell);
        Build { model, tol: tolerance::LINEAR.max(1e-9 * size), vertices: Vec::new(), shell }
    }

    /// The vertex at `p`, made if none is there yet.
    fn vertex(&mut self, p: Point) -> VertexId {
        if let Some(&v) = self.vertices.iter().find(|&&v| self.model[v].point.distance(p) <= self.tol) {
            return v;
        }
        let v = self.model.add_vertex(p, tolerance::LINEAR);
        self.vertices.push(v);
        v
    }

    /// An edge along `curve` over `[a, b]`.
    fn edge(&mut self, curve: Curve, a: f64, b: f64) -> Result<EdgeId> {
        let (start, end) = (self.vertex(curve.point(a)), self.vertex(curve.point(b)));
        let e = self.model.add_edge(start, end, tolerance::LINEAR)?;
        self.model.set_edge_curve(e.into(), curve, tolerance::LINEAR)?;
        self.model.set_edge_range(e.into(), a, b)?;
        Ok(e)
    }

    fn line(&mut self, from: Point, along: Vector) -> Result<EdgeId> {
        self.edge(Curve::Line(Line::new(from, Direction::from_vector(along))), 0.0, along.length())
    }

    fn circle(&mut self, frame: Frame, r: f64, a: f64, b: f64) -> Result<EdgeId> {
        self.edge(Curve::Circle(Circle::new(frame, r)), a, b)
    }

    /// A face on `[u1, u2] × [v1, v2]` of `surface` bounded by `sides`, facing away from
    /// `inside`.
    fn rectangle(&mut self, surface: Surface, (u1, u2, v1, v2): (f64, f64, f64, f64), sides: [Side; 4], inside: Point) -> Result<()> {
        let mut uses: Vec<(EdgeId, Orientation, Curve2, f64, f64)> = Vec::with_capacity(4);
        for (i, side) in sides.into_iter().enumerate() {
            // The side's own parameter s over [sa, sb], increasing, its point in (u, v), and
            // whether the loop runs along it with s.
            let (sa, sb) = if i % 2 == 0 { (u1, u2) } else { (v1, v2) };
            let uv = move |s: f64| match i {
                0 => Point2::new(s, v1),
                1 => Point2::new(u2, s),
                2 => Point2::new(s, v2),
                _ => Point2::new(u1, s),
            };
            let with_s = i < 2;
            let unit = if i % 2 == 0 { Direction2::new(1.0, 0.0) } else { Direction2::new(0.0, 1.0) };
            match side {
                Side::Given(e, o, pc, a, b) => uses.push((e, o, pc, a, b)),
                Side::Pole => {
                    let mid = uv(0.5 * (sa + sb));
                    let v = self.vertex(surface.point(mid.x, mid.y));
                    let e = self.model.add_edge(v, v, tolerance::LINEAR)?;
                    let origin = if i % 2 == 0 { Point2::new(0.0, uv(0.0).y) } else { Point2::new(uv(0.0).x, 0.0) };
                    uses.push((e, if with_s { Orientation::Forward } else { Orientation::Reversed }, Curve2::Line(Line2::new(origin, unit)), sa, sb));
                }
                Side::Edge(e) => {
                    let (curve, t0, t1) = self.model.edge_curve(Shape::from(e)).ok_or_else(|| Error::Failed("an edge without a curve".into()))?;
                    let on = |s: f64| {
                        let q = uv(s);
                        surface.point(q.x, q.y)
                    };
                    let k = (sb - sa) / (t1 - t0);
                    // s at the edge's t, the side's way or the other.
                    let maps = [(sa - k * t0, k), (sb + k * t0, -k)];
                    let fits = |&(s0, ds): &(f64, f64)| [0.0, 0.25, 0.5, 0.75, 1.0].iter().all(|&f| {
                        let t = t0 + (t1 - t0) * f;
                        curve.point(t).distance(on(s0 + ds * t)) <= 10.0 * self.tol
                    });
                    let &(s0, ds) = maps.iter().find(|m| fits(m)).ok_or_else(|| Error::Failed("an edge not along its side".into()))?;
                    let start = uv(s0);
                    let pcurve = if (ds.abs() - 1.0).abs() <= 1e-12 {
                        let d = if ds > 0.0 { unit } else { -unit };
                        Curve2::Line(Line2::new(start, d))
                    } else {
                        let (p0, p1) = (uv(s0 + ds * t0), uv(s0 + ds * t1));
                        Curve2::BSpline(BSplineCurve2::new(vec![p0, p1], vec![t0, t1], vec![2, 2], 1))
                    };
                    let orientation = if with_s == (ds > 0.0) { Orientation::Forward } else { Orientation::Reversed };
                    uses.push((e, orientation, pcurve, t0, t1));
                }
            }
        }
        let (um, vm) = (0.5 * (u1 + u2), 0.5 * (v1 + v2));
        let (p, su, sv) = surface.d1(um, vm);
        let outward = su.cross(sv).dot(p - inside) > 0.0;
        self.face(surface, vec![uses], outward)
    }

    /// A face on a plane facing along `normal`, bounded by the loop of `edges` (end to end, in
    /// order, either way round).
    fn planar(&mut self, origin: Point, normal: Vector, edges: &[EdgeId]) -> Result<()> {
        let plane = Plane::new(origin, Direction::from_vector(normal));
        let (x, y) = (Vector::from(plane.frame().x()), Vector::from(plane.frame().y()));
        let flat = |p: Point| Point2::new((p - origin).dot(x), (p - origin).dot(y));
        let flat_v = |v: Vector| (v.dot(x), v.dot(y));
        // Each edge's way round: from where the one before ends.
        let ends = |m: &Model, e: EdgeId| {
            let (a, b) = m.edge_end_vertices(e);
            (a, b)
        };
        let mut uses = Vec::new();
        let (a0, b0) = ends(self.model, edges[0]);
        let (a1, b1) = ends(self.model, edges[edges.len().min(2) - 1]);
        let mut at = if b0 == a1 || b0 == b1 { b0 } else { a0 };
        if at == a0 && !(a0 == a1 || a0 == b1) {
            return Err(Error::Failed("a planar face's edges not end to end".into()));
        }
        let first_forward = at == b0;
        for (i, &e) in edges.iter().enumerate() {
            let (a, b) = ends(self.model, e);
            let forward = if i == 0 {
                first_forward
            } else if a == at {
                true
            } else if b == at {
                false
            } else {
                return Err(Error::Failed("a planar face's edges not end to end".into()));
            };
            if i > 0 {
                at = if forward { b } else { a };
            }
            let (curve, t0, t1) = self.model.edge_curve(Shape::from(e)).ok_or_else(|| Error::Failed("an edge without a curve".into()))?;
            let pcurve = match &curve {
                Curve::Line(l) => {
                    let (dx, dy) = flat_v(Vector::from(l.direction));
                    Curve2::Line(Line2::new(flat(l.origin), Direction2::new(dx, dy)))
                }
                Curve::Circle(c) => {
                    let (fx, fy) = (flat_v(Vector::from(c.frame().x())), flat_v(Vector::from(c.frame().y())));
                    let ccw = fx.0 * fy.1 - fx.1 * fy.0 > 0.0;
                    Curve2::Circle(Circle2::new(Frame2::new(flat(c.frame().origin()), Direction2::new(fx.0, fx.1), ccw), c.radius()))
                }
                _ => return Err(Error::Failed("a planar face's edge neither a line nor a circle".into())),
            };
            uses.push((e, if forward { Orientation::Forward } else { Orientation::Reversed }, pcurve, t0, t1));
        }
        // Counter-clockwise about the normal.
        let mut points = Vec::new();
        for (_, o, pc, t0, t1) in &uses {
            for j in 0..16 {
                let f = j as f64 / 16.0;
                let t = if *o == Orientation::Forward { t0 + (t1 - t0) * f } else { t1 - (t1 - t0) * f };
                points.push(pc.point(t));
            }
        }
        let n = points.len();
        let area: f64 = (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum();
        if area < 0.0 {
            uses.reverse();
            for u in &mut uses {
                u.1 = u.1.reversed();
            }
        }
        self.face(Surface::Plane(plane), vec![uses], true)
    }

    fn face(&mut self, surface: Surface, loops: Vec<Vec<(EdgeId, Orientation, Curve2, f64, f64)>>, forward: bool) -> Result<()> {
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        for uses in loops {
            let mut coedges = Vec::with_capacity(uses.len());
            for (edge, o, pcurve, first, last) in uses {
                let coedge = self.model.add_coedge(edge, o)?;
                self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
                coedges.push(coedge);
            }
            let wire = self.model.add_wire(&coedges)?;
            self.model.add(face, wire.into())?;
        }
        let orientation = if forward { Orientation::Forward } else { Orientation::Reversed };
        self.model.add(self.shell, Shape { orientation, ..face })
    }

    /// The solid the faces bound.
    fn solid(self) -> Result<Shape> {
        let volume = crate::properties::closed_shell_volume(self.model, self.shell);
        if !(volume > 0.0) {
            return Err(Error::Failed("an exact sweep's faces bound no volume".into()));
        }
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, self.shell)?;
        Ok(solid)
    }
}

/// A direction square to `d`.
fn square_to(d: Vector) -> Vector {
    let other = if d.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    let s = d.cross(other);
    s * (1.0 / s.length())
}

fn unit(v: Vector) -> Vector {
    v * (1.0 / v.length())
}

/// The ball about `c`.
fn whole_ball(model: &mut Model, c: Point, r: f64) -> Result<Shape> {
    let mut b = Build::new(model, r + (c - Point::ORIGIN).length());
    let sphere = Sphere::new(Frame::from_z(c, Direction::Z), r);
    let seam = b.edge(sphere.u_iso(0.0), -FRAC_PI_2, FRAC_PI_2)?;
    b.rectangle(Surface::Sphere(sphere), (0.0, TAU, -FRAC_PI_2, FRAC_PI_2), [Side::Pole, Side::Edge(seam), Side::Pole, Side::Edge(seam)], c)?;
    b.solid()
}

/// The ball about `c` slid by `by`: a cylinder between two half spheres.
fn capsule(model: &mut Model, c: Point, by: Vector, r: f64) -> Result<Shape> {
    let length = by.length();
    let mut b = Build::new(model, r + length + (c - Point::ORIGIN).length());
    let along = unit(by);
    let x = square_to(along);
    let end = c + by;
    let frame = |o: Point| Frame::new(o, Direction::from_vector(along), Direction::from_vector(x));
    let start_circle = b.circle(frame(c), r, 0.0, TAU)?;
    let end_circle = b.circle(frame(end), r, 0.0, TAU)?;
    let seam = b.line(c + x * r, by)?;
    // The half spheres' seams: their meridians through c + r·x, from the pole to the equator
    // behind and from the equator to the pole ahead.
    let start_sphere = Sphere::new(frame(c), r);
    let end_sphere = Sphere::new(frame(end), r);
    let start_meridian = b.edge(start_sphere.u_iso(0.0), -FRAC_PI_2, 0.0)?;
    let end_meridian = b.edge(end_sphere.u_iso(0.0), 0.0, FRAC_PI_2)?;
    let mid = c + by * 0.5;
    b.rectangle(Surface::Cylinder(Cylinder::new(frame(c), r)), (0.0, TAU, 0.0, length), [Side::Edge(start_circle), Side::Edge(seam), Side::Edge(end_circle), Side::Edge(seam)], mid)?;
    b.rectangle(Surface::Sphere(start_sphere), (0.0, TAU, -FRAC_PI_2, 0.0), [Side::Pole, Side::Edge(start_meridian), Side::Edge(start_circle), Side::Edge(start_meridian)], c)?;
    b.rectangle(Surface::Sphere(end_sphere), (0.0, TAU, 0.0, FRAC_PI_2), [Side::Edge(end_circle), Side::Edge(end_meridian), Side::Pole, Side::Edge(end_meridian)], end)?;
    b.solid()
}

/// The ball about `c` turned by `angle` (positive, at most a whole turn) about `axis`, which
/// lies farther from `c` than `r`: the torus piece and, short of a whole turn, the half balls
/// at its ends.
fn torus_piece(model: &mut Model, c: Point, r: f64, axis: Axis, angle: f64) -> Result<Shape> {
    let z = Vector::from(axis.direction);
    let foot = axis.origin + z * (c - axis.origin).dot(z);
    let big = c.distance(foot);
    let mut b = Build::new(model, big + r + (foot - Point::ORIGIN).length());
    let x = unit(c - foot);
    let torus = Torus::new(Frame::new(foot, axis.direction, Direction::from_vector(x)), big, r);
    let whole = angle >= TAU;
    // The torus's u runs round the axis, v round the tube: the meridian circles at the ends
    // (u-isos), the seam round the outside (the v-iso at 0).
    let start_circle = b.edge(torus.u_iso(0.0), 0.0, TAU)?;
    let outside = b.edge(torus.v_iso(0.0), 0.0, angle)?;
    let core = |u: f64| foot + (x * u.cos() + z.cross(x) * u.sin()) * big;
    if whole {
        b.rectangle(Surface::Torus(torus), (0.0, TAU, 0.0, TAU), [Side::Edge(outside), Side::Edge(start_circle), Side::Edge(outside), Side::Edge(start_circle)], core(PI))?;
        return b.solid();
    }
    let end_circle = b.edge(torus.u_iso(angle), 0.0, TAU)?;
    let mid = 0.5 * angle;
    b.rectangle(Surface::Torus(torus), (0.0, angle, 0.0, TAU), [Side::Edge(outside), Side::Edge(end_circle), Side::Edge(outside), Side::Edge(start_circle)], core(mid))?;
    // The half balls: about each end's centre, the pole behind the motion at the start and
    // ahead of it at the end, the equator the meridian circle there.
    for (u, ahead) in [(0.0, false), (angle, true)] {
        let centre = core(u);
        let radial = unit(centre - foot);
        let tangent = z.cross(radial);
        let pole = if ahead { tangent } else { -tangent };
        let sphere = Sphere::new(Frame::new(centre, Direction::from_vector(pole), Direction::from_vector(radial)), r);
        let meridian = b.edge(sphere.u_iso(0.0), 0.0, FRAC_PI_2)?;
        let circle = if ahead { end_circle } else { start_circle };
        b.rectangle(Surface::Sphere(sphere), (0.0, TAU, 0.0, FRAC_PI_2), [Side::Edge(circle), Side::Edge(meridian), Side::Pole, Side::Edge(meridian)], centre)?;
    }
    b.solid()
}

/// The ball about `c` screwed by `angle` (positive) about `axis`, rising `rise`: the tube about
/// its centre's helix, in quarter turns, and the half balls at its ends.
fn screwed_ball(model: &mut Model, c: Point, r: f64, axis: Axis, angle: f64, rise: f64) -> Result<Shape> {
    let z = Vector::from(axis.direction);
    let foot = axis.origin + z * (c - axis.origin).dot(z);
    let big = c.distance(foot);
    let size = big + r + rise.abs() + (foot - Point::ORIGIN).length();
    let x = unit(c - foot);
    let tangent = unit(z.cross(x) * (angle * big) + z * rise);
    let stages = ((angle / FRAC_PI_2).ceil() as usize).max(1);
    let (step_angle, step_rise) = (angle / stages as f64, rise / stages as f64);
    let step = Transform::translation(z * step_rise) * Transform::rotation(axis, step_angle);
    let motion = product::screw_motion(axis, step_angle, step_rise, size, big + r)?;
    // The circle square to the helix at the start, as a B-spline, screwed.
    let circle = Curve::Circle(Circle::new(Frame::new(c, Direction::from_vector(tangent), Direction::from_vector(x)), r));
    let profile = product::Profile::of(&circle, 0.0, TAU, size)?;
    let knots = path(&motion, Point::ORIGIN)?;
    let knots = KnotVector::new(knots.knots().to_vec(), knots.multiplicities().to_vec(), knots.degree());
    let surface = product::profiled_surface(&profile.spline, &|q| Ok(path(&motion, q)?.poles().to_vec()), knots)?;
    let (u0, u1) = (profile.spline.first_parameter(), profile.spline.last_parameter());
    // The circle's pcurve on a stage's surface at v: its angle mapped to the B-spline's
    // parameter.
    let on_tube = |v: f64| -> Curve2 {
        match &profile.sigma {
            Some(sigma) => {
                let poles = sigma.poles().iter().map(|q| Point2::new(q.x, q.y + v)).collect();
                Curve2::BSpline(BSplineCurve2::new(poles, sigma.knots().to_vec(), sigma.multiplicities().to_vec(), sigma.degree()))
            }
            None => Curve2::Line(Line2::new(Point2::new(0.0, v), Direction2::new(1.0, 0.0))),
        }
    };
    let mut b = Build::new(model, size);
    let mut places = vec![Transform::IDENTITY];
    for k in 0..stages {
        places.push(step * places[k]);
    }
    let circles: Vec<EdgeId> = places.iter().map(|m| b.edge(circle.transformed(m), 0.0, TAU)).collect::<Result<_>>()?;
    let seam_curve = surface.u_iso(u0);
    for k in 0..stages {
        let stage = surface.transformed(&places[k]);
        let seam = b.edge(seam_curve.transformed(&places[k]), 0.0, 1.0)?;
        let inside = c.transformed(&(Transform::translation(z * (0.5 * step_rise)) * Transform::rotation(axis, 0.5 * step_angle) * places[k]));
        let sides = [
            Side::Given(circles[k], Orientation::Forward, on_tube(0.0), 0.0, TAU),
            Side::Edge(seam),
            Side::Given(circles[k + 1], Orientation::Reversed, on_tube(1.0), 0.0, TAU),
            Side::Edge(seam),
        ];
        b.rectangle(stage, (u0, u1, 0.0, 1.0), sides, inside)?;
    }
    // The half balls: behind the motion at the start, ahead at the end.
    for (k, ahead) in [(0, false), (stages, true)] {
        let m = places[k];
        let centre = c.transformed(&m);
        let (t, radial) = (tangent.transformed(&m), x.transformed(&m));
        let pole = if ahead { t } else { -t };
        let sphere = Sphere::new(Frame::new(centre, Direction::from_vector(pole), Direction::from_vector(radial)), r);
        let meridian = b.edge(sphere.u_iso(0.0), 0.0, FRAC_PI_2)?;
        b.rectangle(Surface::Sphere(sphere), (0.0, TAU, 0.0, FRAC_PI_2), [Side::Edge(circles[k]), Side::Edge(meridian), Side::Pole, Side::Edge(meridian)], centre)?;
    }
    b.solid()
}

/// A ball-end mill (the ball about `c` of radius `r`, the cylinder up `a` to `h` above `c`)
/// slid by `d` (see the module's notes).
fn ball_end_slid(model: &mut Model, c: Point, a: Vector, r: f64, h: f64, d: Vector) -> Result<Shape> {
    let length = d.length();
    let size = r + h + length + (c - Point::ORIGIN).length();
    let rise = d.dot(a);
    let across = d - a * rise;
    if across.length() <= 1e-9 * length {
        // (Along its axis: the mill lengthened.)
        let low = if rise > 0.0 { c } else { c + d };
        return mill(model, low, a, r, h + rise.abs(), size);
    }
    let mut b = Build::new(model, size);
    let e1 = unit(across);
    let s = a.cross(e1);
    let ahead = unit(d);
    let beta = rise.atan2(across.length());
    // (The up direction square to the slide in the plane of the axis and the slide.)
    let w = ahead.cross(s);
    let dir = Direction::from_vector;
    let end = c + d;
    let top = a * h;
    // The outline's half circle below the equator, its lines along the slide, the equators'
    // halves, the cylinder's lines up from them, the rims' halves, the top's lines.
    let below = |o: Point| Frame::new(o, dir(ahead), dir(s));
    let level = |o: Point| Frame::new(o, dir(a), dir(s));
    let low_start = b.circle(below(c), r, PI, TAU)?;
    let low_end = b.circle(below(end), r, PI, TAU)?;
    let side_plus = b.line(c + s * r, d)?;
    let side_minus = b.line(c - s * r, d)?;
    let equator_behind = b.circle(level(c), r, 0.0, PI)?;
    let equator_ahead = b.circle(level(end), r, PI, TAU)?;
    let up_plus_start = b.line(c + s * r, top)?;
    let up_minus_start = b.line(c - s * r, top)?;
    let up_plus_end = b.line(end + s * r, top)?;
    let up_minus_end = b.line(end - s * r, top)?;
    let rim_behind_start = b.circle(level(c + top), r, 0.0, PI)?;
    let rim_ahead_end = b.circle(level(end + top), r, PI, TAU)?;
    let top_plus = b.line(c + top + s * r, d)?;
    let top_minus = b.line(c + top - s * r, d)?;
    let mid = c + d * 0.5;
    // The half cylinder about the centre's path.
    b.rectangle(Surface::Cylinder(Cylinder::new(below(c), r)), (PI, TAU, 0.0, length), [Side::Edge(low_start), Side::Edge(side_plus), Side::Edge(low_end), Side::Edge(side_minus)], mid)?;
    // The cylinder's halves behind at the start and ahead at the end.
    b.rectangle(Surface::Cylinder(Cylinder::new(level(c), r)), (0.0, PI, 0.0, h), [Side::Edge(equator_behind), Side::Edge(up_minus_start), Side::Edge(rim_behind_start), Side::Edge(up_plus_start)], c + top * 0.5)?;
    b.rectangle(Surface::Cylinder(Cylinder::new(level(end), r)), (PI, TAU, 0.0, h), [Side::Edge(equator_ahead), Side::Edge(up_plus_end), Side::Edge(rim_ahead_end), Side::Edge(up_minus_end)], end + top * 0.5)?;
    // The half ball's lunes, on spheres whose poles are the outline's ends (±s), between the
    // meridian of the outline's half circle (u = 0) and the equator's half.
    let lune = |o: Point| Sphere::new(Frame::new(o, dir(s), dir(-w)), r);
    b.rectangle(Surface::Sphere(lune(c)), (0.0, beta + FRAC_PI_2, -FRAC_PI_2, FRAC_PI_2), [Side::Pole, Side::Edge(equator_behind), Side::Pole, Side::Edge(low_start)], c)?;
    b.rectangle(Surface::Sphere(lune(end)), (beta + 1.5 * PI, TAU, -FRAC_PI_2, FRAC_PI_2), [Side::Pole, Side::Edge(low_end), Side::Pole, Side::Edge(equator_ahead)], end)?;
    // The rim's half that the slide carries in under the top, and the top where the slide
    // leaves it; level, the top swept along itself.
    let (wall_plus_top, wall_minus_top) = (top_plus, top_minus);
    if rise.abs() <= 1e-12 * length {
        b.planar(c + top, a, &[rim_behind_start, top_minus, rim_ahead_end, top_plus])?;
    } else {
        let up = rise > 0.0;
        let rim = Curve::Circle(Circle::new(level(c + top), r));
        let (u1, u2) = if up { (0.0, PI) } else { (PI, TAU) };
        let carried = ExtrusionSurface::new(rim, dir(ahead));
        if up {
            let rim_behind_end = b.circle(level(end + top), r, 0.0, PI)?;
            b.rectangle(Surface::Extrusion(carried), (u1, u2, 0.0, length), [Side::Edge(rim_behind_start), Side::Edge(top_minus), Side::Edge(rim_behind_end), Side::Edge(top_plus)], mid)?;
            b.planar(end + top, a, &[rim_behind_end, rim_ahead_end])?;
        } else {
            let rim_ahead_start = b.circle(level(c + top), r, PI, TAU)?;
            b.rectangle(Surface::Extrusion(carried), (u1, u2, 0.0, length), [Side::Edge(rim_ahead_start), Side::Edge(top_plus), Side::Edge(rim_ahead_end), Side::Edge(top_minus)], mid)?;
            b.planar(c + top, a, &[rim_behind_start, rim_ahead_start])?;
        }
    }
    // The two planes the cylinder's lines square to the slide sweep.
    b.planar(c + s * r, s, &[side_plus, up_plus_end, wall_plus_top, up_plus_start])?;
    b.planar(c - s * r, -s, &[side_minus, up_minus_end, wall_minus_top, up_minus_start])?;
    b.solid()
}

/// The ball-end mill about `c`: the half ball below, the cylinder up `a` to `h` above `c`.
fn mill(model: &mut Model, c: Point, a: Vector, r: f64, h: f64, size: f64) -> Result<Shape> {
    let mut b = Build::new(model, size);
    let x = square_to(a);
    let frame = |o: Point| Frame::new(o, Direction::from_vector(a), Direction::from_vector(x));
    let equator = b.circle(frame(c), r, 0.0, TAU)?;
    let rim = b.circle(frame(c + a * h), r, 0.0, TAU)?;
    let seam = b.line(c + x * r, a * h)?;
    let sphere = Sphere::new(frame(c), r);
    let meridian = b.edge(sphere.u_iso(0.0), -FRAC_PI_2, 0.0)?;
    b.rectangle(Surface::Sphere(sphere), (0.0, TAU, -FRAC_PI_2, 0.0), [Side::Pole, Side::Edge(meridian), Side::Edge(equator), Side::Edge(meridian)], c)?;
    b.rectangle(Surface::Cylinder(Cylinder::new(frame(c), r)), (0.0, TAU, 0.0, h), [Side::Edge(equator), Side::Edge(seam), Side::Edge(rim), Side::Edge(seam)], c + a * (0.5 * h))?;
    b.planar(c + a * h, a, &[rim])?;
    b.solid()
}
