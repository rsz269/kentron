//! Rolled fillet faces made exact where they are: a ball of one radius rolling with its
//! centre along a line sweeps a cylinder, round a circle a torus, about a point a sphere.
//! The rolled fillet is a B-spline fitted to the ball's arcs; the exact surface is what OCCT
//! makes of such fillets, and what later fillets and booleans on them take as one (their
//! edges with planes lines and circles, their fillets in turn elementary).

use std::collections::HashSet;

use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::curve2::Curve2;
use crate::geometry::project;
use crate::geometry::surface::{Cylinder, Sphere, Surface, Torus};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Frame, Point, Transform2, Transformable2, Vector, Vector2};
use crate::topology::{CoEdgeId, Model, NodeId, NodeKind, Orientation, Shape};

/// What made a face exact, to undo it.
pub(super) struct Undo {
    face: Shape,
    surface: Surface,
    tolerance: f64,
    pcurves: Vec<(CoEdgeId, Curve2, f64, f64)>,
    /// Its wires and its uses in shells, turned round (the exact surface facing the other way).
    turned: Vec<(Shape, Shape)>,
}

/// What made an edge exact, to undo it.
pub(super) struct EdgeUndo {
    edge: Shape,
    curve: Curve,
    range: (f64, f64),
    tolerance: f64,
    pcurves: Vec<(CoEdgeId, Curve2, f64, f64)>,
}

/// What [`exact_faces`] did.
pub(super) struct Done {
    faces: Vec<Undo>,
    edges: Vec<EdgeUndo>,
}

impl Done {
    pub(super) fn is_empty(&self) -> bool {
        self.faces.is_empty() && self.edges.is_empty()
    }
}

/// The B-spline faces of `shape` not among `before` (faces a fillet made) that are a
/// cylinder, torus or sphere to within their tolerance, made so: their surfaces replaced,
/// their pcurves found again on them; and the B-spline edges of those faces that are lines
/// or circular arcs, made so (their pcurves found again on their faces). Returns what was
/// done, to [`undo`] it.
pub(super) fn exact_faces(model: &mut Model, shape: Shape, before: &HashSet<NodeId>, before_edges: &HashSet<NodeId>) -> Done {
    let mut faces = Vec::new();
    let mut new_faces = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        if before.contains(&face.node) {
            continue;
        }
        let face = Shape { orientation: Orientation::Forward, ..face };
        new_faces.push(face);
        if let Some(undo) = make_exact(model, shape, face) {
            faces.push(undo);
        }
    }
    let mut edges = Vec::new();
    let mut seen = HashSet::new();
    for face in new_faces {
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            // (Only the fillet's own edges: the shape's, which other shapes may share, stay.)
            if !before_edges.contains(&edge.node) && seen.insert(edge.node) {
                if let Some(undo) = exact_edge(model, shape, Shape { orientation: Orientation::Forward, ..edge }) {
                    edges.push(undo);
                }
            }
        }
    }
    Done { faces, edges }
}

/// Undoes [`exact_faces`].
pub(super) fn undo(model: &mut Model, done: Done) {
    // (The records are the fillet's own: put back as they were.)
    for undo in done.edges.into_iter().rev() {
        let _ = model.set_edge_curve(undo.edge, undo.curve, undo.tolerance);
        let _ = model.set_edge_range(undo.edge, undo.range.0, undo.range.1);
        for (c, p, a, b) in undo.pcurves {
            let _ = model.set_pcurve_with_range(c, p, a, b);
        }
    }
    for undo in done.faces.into_iter().rev() {
        for (parent, turned) in undo.turned.iter().rev() {
            let _ = model.replace_child(*parent, *turned, turned.reversed());
        }
        let _ = model.set_face_surface(undo.face, undo.surface, undo.tolerance);
        for (c, p, a, b) in undo.pcurves {
            let _ = model.set_pcurve_with_range(c, p, a, b);
        }
    }
}

/// A circle through three points: its centre, radius and the normal of its plane.
fn circle_through(a: Point, b: Point, c: Point) -> Option<(Point, f64, Vector)> {
    let (ab, ac) = (b - a, c - a);
    let n = ab.cross(ac);
    let n2 = n.dot(n);
    if !(n2 > 1e-24 * ab.dot(ab) * ac.dot(ac)) {
        return None;
    }
    let to_centre = (n.cross(ab) * ac.dot(ac) + ac.cross(n) * ab.dot(ab)) * (0.5 / n2);
    let centre = a + to_centre;
    Some((centre, to_centre.length(), n * (1.0 / n2.sqrt())))
}

/// A unit vector square to `v`.
fn square_to(v: Vector) -> Vector {
    let t = if v.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    v.cross(t).normalized()
}

/// The cylinder, torus or sphere the B-spline `surface` is over `bounds` (to within `tol`),
/// its seam away from the middle of the face (at `middle`).
fn exact_surface(surface: &Surface, bounds: (f64, f64, f64, f64), tol: f64) -> Option<Surface> {
    let (u0, u1, v0, v1) = bounds;
    // The arcs across (u) at several places along (v): their centres and radii.
    const ALONG: usize = 9;
    let mut arcs = Vec::with_capacity(ALONG);
    for k in 0..ALONG {
        let v = v0 + (v1 - v0) * k as f64 / (ALONG - 1) as f64;
        let (a, b, c) = (surface.point(u0, v), surface.point(0.5 * (u0 + u1), v), surface.point(u1, v));
        arcs.push(circle_through(a, b, c)?);
    }
    let radius = arcs[0].1;
    if arcs.iter().any(|a| (a.1 - radius).abs() > tol) {
        return None;
    }
    let middle = surface.point(0.5 * (u0 + u1), 0.5 * (v0 + v1));
    let centres: Vec<Point> = arcs.iter().map(|a| a.0).collect();
    let (first, last) = (centres[0], centres[ALONG - 1]);
    let candidate = if centres.iter().all(|c| c.distance(first) <= tol) {
        // About a point: a sphere, its poles square to the face's middle.
        let out = Direction::try_from_vector(middle - first)?;
        let z = square_to(Vector::from(out));
        Surface::Sphere(Sphere::new(Frame::new(first, Direction::from_vector(z), Direction::from_vector(-Vector::from(out))), radius))
    } else {
        let along = last - first;
        let d = Vector::from(Direction::try_from_vector(along)?);
        let off_line = |c: Point| ((c - first) - d * (c - first).dot(d)).length();
        if centres.iter().all(|&c| off_line(c) <= tol) && arcs.iter().all(|a| a.2.cross(d).length() <= 1e-6) {
            // Along a line: a cylinder, its seam opposite the face's middle.
            let radial = Direction::try_from_vector((middle - first) - d * (middle - first).dot(d))?;
            Surface::Cylinder(Cylinder::new(Frame::new(first, Direction::from_vector(d), Direction::from_vector(-Vector::from(radial))), radius))
        } else {
            // Round a circle: a torus.
            let (centre, major, axis) = circle_through(first, centres[ALONG / 2], last)?;
            if centres.iter().any(|&c| {
                let r = c - centre;
                r.dot(axis).abs() > tol || ((r - axis * r.dot(axis)).length() - major).abs() > tol
            }) {
                return None;
            }
            // (Each arc in a plane through the axis: square to the circle's tangent there.)
            if arcs.iter().any(|a| a.2.dot(axis).abs() > 1e-6 || a.2.dot(a.0 - centre).abs() > 1e-6 * major.max(1.0)) {
                return None;
            }
            if !(major > tol) {
                return None;
            }
            let radial = Direction::try_from_vector((middle - centre) - axis * (middle - centre).dot(axis))?;
            Surface::Torus(Torus::new(Frame::new(centre, Direction::from_vector(axis), Direction::from_vector(-Vector::from(radial))), major, radius))
        }
    };
    // The face's surface on it everywhere.
    const GRID: usize = 11;
    for i in 0..GRID {
        for j in 0..GRID {
            let (u, v) = (u0 + (u1 - u0) * i as f64 / (GRID - 1) as f64, v0 + (v1 - v0) * j as f64 / (GRID - 1) as f64);
            let p = surface.point(u, v);
            let (s, t) = candidate.parameters(p)?;
            if candidate.point(s, t).distance(p) > tol {
                return None;
            }
        }
    }
    Some(candidate)
}

/// `face` (a face of `shape`) made exact, or `None` where it is not one of those surfaces or
/// its pcurves cannot be found on it.
fn make_exact(model: &mut Model, shape: Shape, face: Shape) -> Option<Undo> {
    let surface = model.surface(face)?.clone();
    if !matches!(surface, Surface::BSpline(_)) {
        return None;
    }
    let tolerance = model.tolerance(face);
    let bounds = model.face_uv_bounds(face);
    let size = {
        let (u0, u1, v0, v1) = bounds;
        surface.point(u0, v0).distance(surface.point(u1, v1)) + surface.point(u0, v1).distance(surface.point(u1, v0))
    };
    let tol = tolerance.max(crate::tolerance::LINEAR).max(1e-9 * size).min(1e-6 * size.max(1.0));
    let exact = exact_surface(&surface, bounds, tol)?;
    // Its pcurves on the exact surface: each edge's curve, near where the face is (the seam
    // away from it), the edge used once.
    let coedges: Vec<CoEdgeId> = model.unique_subshapes(face, NodeKind::CoEdge).into_iter().filter_map(|c| CoEdgeId::try_from(c.node).ok()).collect();
    let mut edges_seen = HashSet::new();
    let mut new_pcurves = Vec::with_capacity(coedges.len());
    let mut old_pcurves = Vec::with_capacity(coedges.len());
    let periods = [exact.u_period(), exact.v_period()];
    let (mu, mv) = exact.parameters(surface.point(0.5 * (bounds.0 + bounds.1), 0.5 * (bounds.2 + bounds.3)))?;
    for &c in &coedges {
        let edge = model[c].edge;
        if !edges_seen.insert(edge) {
            return None;
        }
        let (curve, first, last) = model.edge_curve(Shape::from(edge))?;
        let pcurve = match project::analytic_pcurve(&curve, &exact) {
            Some(p) => p,
            None => Curve2::BSpline(project::pcurve_on_range_within(&curve, &exact, first, last, &[1e-9, 1e-7, 1e-5])?),
        };
        // (Each in the periods about the face's middle.)
        let mid = pcurve.point(0.5 * (first + last));
        let shift = |period: Option<f64>, target: f64, at: f64| period.map_or(0.0, |p| ((target - at) / p).round() * p);
        let (du, dv) = (shift(periods[0], mu, mid.x), shift(periods[1], mv, mid.y));
        let pcurve = if du != 0.0 || dv != 0.0 { pcurve.transformed(&Transform2::translation(Vector2::new(du, dv))) } else { pcurve };
        // (On the edge, at the edge's own parameters.)
        let gap = (0..=8).map(|i| {
            let t = first + (last - first) * i as f64 / 8.0;
            let uv = pcurve.point(t);
            exact.point(uv.x, uv.y).distance(curve.point(t))
        });
        if gap.fold(0.0, f64::max) > tol.max(model.tolerance(Shape::from(edge))) {
            return None;
        }
        let (old, a, b) = model.effective_pcurve(c)?;
        old_pcurves.push((c, old, a, b));
        new_pcurves.push((c, pcurve, first, last));
    }
    // Facing as the face did: else its wires turned round, and its uses in shells.
    let (p, du, dv) = surface.d1(0.5 * (bounds.0 + bounds.1), 0.5 * (bounds.2 + bounds.3));
    let old_normal = du.cross(dv);
    let (s, t) = exact.parameters(p)?;
    let (_, eu, ev) = exact.d1(s, t);
    let turn = old_normal.dot(eu.cross(ev)) < 0.0;
    let mut turned = Vec::new();
    if turn {
        for wire in model.children(face).collect::<Vec<_>>() {
            model.replace_child(face, wire, wire.reversed()).ok()?;
            turned.push((face, wire.reversed()));
        }
        for shell in model.unique_subshapes(shape, NodeKind::Shell) {
            let shell = Shape { orientation: Orientation::Forward, ..shell };
            if let Some(used) = model.children(shell).find(|f| f.node == face.node) {
                model.replace_child(shell, used, used.reversed()).ok()?;
                turned.push((shell, used.reversed()));
            }
        }
    }
    model.set_face_surface(face, exact, tolerance).ok()?;
    for (c, p, a, b) in new_pcurves {
        model.set_pcurve_with_range(c, p, a, b).ok()?;
    }
    Some(Undo { face, surface, tolerance, pcurves: old_pcurves, turned })
}

/// `edge` (of `shape`), a B-spline on a straight line or a circle to within its tolerance,
/// made a line or a circular arc: its pcurves on its faces found again for the new curve's
/// parameters. `None` where it is neither, or a pcurve cannot be found.
fn exact_edge(model: &mut Model, shape: Shape, edge: Shape) -> Option<EdgeUndo> {
    let (curve, first, last) = model.edge_curve(edge)?;
    if !matches!(curve, Curve::BSpline(_)) {
        return None;
    }
    let tolerance = model.tolerance(edge);
    const N: usize = 16;
    let samples: Vec<Point> = (0..=N).map(|i| curve.point(first + (last - first) * i as f64 / N as f64)).collect();
    let (start, end) = (samples[0], samples[N]);
    let length: f64 = samples.windows(2).map(|w| w[0].distance(w[1])).sum();
    let tol = tolerance.max(crate::tolerance::LINEAR).min(1e-6 * length.max(1.0));
    // A straight line from start to end, or the circle through start, middle and end, with
    // the new parameter at each sample (increasing along the edge).
    let (exact, params): (Curve, Vec<f64>) = {
        let chord = end - start;
        let line = Direction::try_from_vector(chord).filter(|_| chord.length() > tol).and_then(|d| {
            let dv = Vector::from(d);
            let ts: Vec<f64> = samples.iter().map(|&p| (p - start).dot(dv)).collect();
            let on = samples.iter().zip(&ts).all(|(&p, &t)| (start + dv * t).distance(p) <= tol);
            let rising = ts.windows(2).all(|w| w[1] > w[0]);
            (on && rising).then(|| (Curve::Line(Line::new(start, d)), ts))
        });
        match line {
            Some(l) => l,
            None => {
                let (centre, radius, normal) = circle_through(start, samples[N / 2], end)?;
                // (The normal the way the edge turns, so that the angle rises along it.)
                let turn = (samples[1] - centre).cross(samples[2] - samples[1]);
                let normal = if turn.dot(normal) < 0.0 { -normal } else { normal };
                let x = Direction::try_from_vector(start - centre)?;
                let frame = Frame::new(centre, Direction::try_from_vector(normal)?, x);
                let circle = Circle::new(frame, radius);
                let mut ts = Vec::with_capacity(samples.len());
                let mut previous = -1.0;
                for &p in &samples {
                    let r = p - centre;
                    let mut t = r.dot(Vector::from(frame.y())).atan2(r.dot(Vector::from(frame.x())));
                    if t < -1e-12 {
                        t += std::f64::consts::TAU;
                    }
                    while t < previous - 1e-9 {
                        t += std::f64::consts::TAU;
                    }
                    if circle.point(t).distance(p) > tol {
                        return None;
                    }
                    ts.push(t);
                    previous = t;
                }
                if !ts.windows(2).all(|w| w[1] > w[0]) {
                    return None;
                }
                (Curve::Circle(circle), ts)
            }
        }
    };
    let (a, b) = (params[0], params[N]);
    // Its pcurves on its faces, at the new parameters, each where the old one was (in the
    // surface's periods).
    let edge_id = crate::topology::EdgeId::try_from(edge.node).ok()?;
    let mut new_pcurves = Vec::new();
    let mut old_pcurves = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let uses: Vec<CoEdgeId> = model.unique_subshapes(face, NodeKind::CoEdge).into_iter().filter_map(|c| CoEdgeId::try_from(c.node).ok()).filter(|&c| model[c].edge == edge_id).collect();
        // (A seam, used twice by one face, a pcurve on either side of it: left as it is.)
        if uses.len() > 1 {
            return None;
        }
        for c in uses {
            let surface = model.surface(Shape { orientation: Orientation::Forward, ..face })?.clone();
            let (old, oa, ob) = model.effective_pcurve(c)?;
            let pcurve = match project::analytic_pcurve(&exact, &surface) {
                Some(p) => p,
                None => Curve2::BSpline(project::pcurve_on_range_within(&exact, &surface, a, b, &[1e-9, 1e-7, 1e-5])?),
            };
            let target = old.point(0.5 * (oa + ob));
            let mid = pcurve.point(0.5 * (a + b));
            let shift = |period: Option<f64>, target: f64, at: f64| period.map_or(0.0, |p| ((target - at) / p).round() * p);
            let (du, dv) = (shift(surface.u_period(), target.x, mid.x), shift(surface.v_period(), target.y, mid.y));
            let pcurve = if du != 0.0 || dv != 0.0 { pcurve.transformed(&Transform2::translation(Vector2::new(du, dv))) } else { pcurve };
            // (On the new curve, and where the old pcurve ran: its ends at the old ones'.)
            let off = (0..=8).map(|i| {
                let t = a + (b - a) * i as f64 / 8.0;
                let uv = pcurve.point(t);
                surface.point(uv.x, uv.y).distance(exact.point(t))
            });
            if off.fold(0.0, f64::max) > tol.max(tolerance) {
                return None;
            }
            for (t_new, t_old) in [(a, oa), (b, ob)] {
                let (p, q) = (pcurve.point(t_new), old.point(t_old));
                if surface.point(p.x, p.y).distance(surface.point(q.x, q.y)) > 10.0 * tol.max(tolerance) || (p - q).length() > 0.5 * surface.u_period().unwrap_or(f64::INFINITY) {
                    return None;
                }
            }
            old_pcurves.push((c, old, oa, ob));
            new_pcurves.push((c, pcurve));
        }
    }
    if new_pcurves.is_empty() {
        return None;
    }
    model.set_edge_curve(edge, exact, tolerance).ok()?;
    model.set_edge_range(edge, a, b).ok()?;
    for (c, p) in new_pcurves {
        model.set_pcurve_with_range(c, p, a, b).ok()?;
    }
    Some(EdgeUndo { edge, curve, range: (first, last), tolerance, pcurves: old_pcurves })
}
