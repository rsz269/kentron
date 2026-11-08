//! Offsets joined by arcs, found as what they are: the points within the distance of the
//! solid's faces (or beyond it, inside, for an offset inwards). Moving the faces one by one
//! and rounding the edges where they part fails where an offset closes a notch or passes a
//! feature narrower than twice it (the moved faces pass through one another, and the
//! rounding's tools reach past faces they do not touch).
//!
//! The boundary of the space within a distance of the faces lies on the faces moved by the
//! distance, on tubes about the edges whose faces' offsets part and on balls about their
//! vertices: those (with the solid's own faces, for a hollowed solid) divide space into
//! cells, each wholly in the offset or out of it, kept or not by the distance of a point
//! inside it. More surfaces only divide the cells further: a curved face is moved over the
//! whole of its parameters' box. Faces on planes, cylinders, cones, spheres and tori (moved
//! through their axis or centre, turned round it: only dividing the cells), tubes about
//! straight and circular edges and pipes about others.

use crate::boolean::{PointClassifier, PointLocation};
use crate::error::Result;
use crate::geometry::surface::{Cone, Cylinder, Sphere, Surface, Torus};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::heal::face_loops;
use crate::math::{Direction, Frame, Point, Point2, Transform, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::regularize::winding;

/// What is made of the solid.
pub(super) enum Made<'a> {
    /// Its faces moved `distance` outwards (inwards if negative).
    Offset,
    /// Hollowed to walls `distance` thick (outwards if positive), open where `removed`
    /// were.
    Thick { removed: &'a [NodeId] },
}

/// A face of the solid: its surface, which way it faces, and its loops in the surface's
/// parameters (as polygons) and edges.
struct Bounded {
    shape: Shape,
    surface: Surface,
    reversed: bool,
    loops: Vec<Vec<Point2>>,
    edges: Vec<(crate::geometry::curve::Curve, f64, f64)>,
    /// The edges' nodes, and whether each is measured from (not where the face meets a removed
    /// face tangentially: the wall ends square to the edge there, see `strips`).
    edge_nodes: Vec<NodeId>,
    measured: Vec<bool>,
    u_period: Option<f64>,
    v_period: Option<f64>,
}

impl Bounded {
    /// The parameters of the foot of the perpendicular from `p` to the surface, in the
    /// period of the face's loops where `p` falls in them.
    fn foot(&self, p: Point) -> (Point2, bool) {
        let (u, v) = self.surface.parameters(p).unwrap_or_else(|| {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &self.surface, None);
            (u, v)
        });
        let shifts = |period: Option<f64>| -> &[f64] { if period.is_some() { &[0.0, -1.0, 1.0, -2.0, 2.0] } else { &[0.0] } };
        for &k in shifts(self.u_period) {
            for &l in shifts(self.v_period) {
                let q = Point2::new(u + k * self.u_period.unwrap_or(0.0), v + l * self.v_period.unwrap_or(0.0));
                if winding(&self.loops, q) != 0 {
                    return (q, true);
                }
            }
        }
        (Point2::new(u, v), false)
    }

    /// The distance from `p` to the face.
    fn distance(&self, p: Point) -> f64 {
        let (q, on) = self.foot(p);
        if on {
            return self.surface.point(q.x, q.y).distance(p);
        }
        self.edges.iter().zip(&self.measured).filter(|(_, &m)| m).map(|((c, a, b), _)| crate::geometry::extrema::nearest_on_curve(p, c, *a, *b).distance).fold(f64::INFINITY, f64::min)
    }

    /// The face's outward normal at `p` (on it).
    fn normal(&self, p: Point) -> Vector {
        let (q, _) = self.foot(p);
        let (_, du, dv) = self.surface.d1(q.x, q.y);
        let n = du.cross(dv);
        let n = n * (1.0 / n.length().max(1e-300));
        if self.reversed { -n } else { n }
    }
}

/// `solid` offset by `distance` as `made` says, joined by arcs; `None` where it is not a
/// solid of faces on planes, cylinders, cones, spheres and tori, or its cells cannot be told
/// apart.
pub(super) fn offset(model: &mut Model, solid: Shape, distance: f64, made: Made) -> Result<Option<Shape>> {
    // (A compound of one solid: that solid.)
    let solid = match (solid.kind(), &model.unique_subshapes(solid, NodeKind::Solid)[..]) {
        (NodeKind::Solid, _) => solid,
        (NodeKind::Compound, [one]) => *one,
        _ => return Ok(None),
    };
    if distance == 0.0 {
        return Ok(None);
    }
    let removed: &[NodeId] = match made {
        Made::Offset => &[],
        Made::Thick { removed } => removed,
    };
    // (Faces placed in the solid (as read from some files): measured on a copy with the
    // placements applied; not where faces are removed, which are named by node.)
    let placed = model.unique_subshapes(solid, NodeKind::Face).iter().any(|f| f.placement != Transform::IDENTITY);
    let solid = match (placed, removed.is_empty()) {
        (false, _) => solid,
        (true, true) => crate::modify::copy_unplaced(model, solid, &[])?.0,
        (true, false) => return Ok(None),
    };
    let mut faces = Vec::new();
    for face in model.unique_subshapes(solid, NodeKind::Face) {
        let Some(surface @ (Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_))) = model.surface(face).cloned() else { return Ok(None) };
        let Some(uv) = face_loops(model, face) else { return Ok(None) };
        let loops = uv.iter().map(|l| l.iter().flat_map(|u| (0..16).map(move |i| u.at(i as f64 / 16.0))).collect()).collect();
        let (mut edges, mut edge_nodes) = (Vec::new(), Vec::new());
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            // (A degenerate edge, at a pole, is a point its neighbours' ends are at already.)
            let Some(curve) = model.edge_curve(edge) else { continue };
            edges.push(curve);
            edge_nodes.push(edge.node);
        }
        let measured = vec![true; edges.len()];
        let (u_period, v_period) = (surface.u_period(), surface.v_period());
        faces.push(Bounded { shape: face, surface, reversed: face.orientation == Orientation::Reversed, loops, edges, edge_nodes, measured, u_period, v_period });
    }
    let kept: Vec<usize> = (0..faces.len()).filter(|&i| !removed.contains(&faces[i].shape.node)).collect();
    let (strips, all_tangent, openings) = strips(model, solid, &mut faces, removed, distance.abs())?;
    // (Walls grown outwards with faces removed end on the removed faces' planes, not where
    // the distance takes them: but where every kept face meets the removed ones
    // tangentially, where they end on the strips.)
    if distance > 0.0 && !removed.is_empty() && !all_tangent {
        return Ok(None);
    }
    let inside = PointClassifier::new(model, solid)?;
    let size = crate::properties::bounding_box(model, solid).size().length().max(1e-9);
    let (sign, d) = (distance.signum(), distance.abs());
    let mut sheets: Vec<Shape> = Vec::new();
    // The kept faces moved: planar ones as they are, others over their parameters' box, and
    // those beside a removed face over a box grown past it (their edges moved there no
    // longer lie on it; the pieces then left loose outside the solid are dropped).
    let removed_edges: Vec<NodeId> = faces
        .iter()
        .filter(|f| removed.contains(&f.shape.node))
        .flat_map(|f| model.unique_subshapes(f.shape, NodeKind::Edge).into_iter().map(|e| e.node))
        .collect();
    for &i in &kept {
        let face = &faces[i];
        let along = if face.reversed { -distance } else { distance };
        let beside_removed = model.unique_subshapes(face.shape, NodeKind::Edge).iter().any(|e| removed_edges.contains(&e.node));
        let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
        for q in face.loops.iter().flatten() {
            lo = Point2::new(lo.x.min(q.x), lo.y.min(q.y));
            hi = Point2::new(hi.x.max(q.x), hi.y.max(q.y));
        }
        // (A plane beside a cone grown too: their offsets need not meet within the plane's
        // own bounds, a cone's moving across its axis as it narrows.)
        let beside_cone = matches!(face.surface, Surface::Plane(_))
            && model.unique_subshapes(face.shape, NodeKind::Edge).iter().any(|e| faces.iter().any(|g| matches!(g.surface, Surface::Cone(_)) && model.unique_subshapes(g.shape, NodeKind::Edge).iter().any(|x| x.node == e.node)));
        let grow = if beside_removed { (4.0 * d).max(0.25 * size) } else if beside_cone { 4.0 * d } else { 0.0 };
        let moved = match &face.surface {
            Surface::Plane(plane) if !beside_removed && !beside_cone => {
                let n = Vector::from(plane.normal());
                crate::modify::transform(model, face.shape, &Transform::translation(n * along), Default::default())?
            }
            Surface::Plane(plane) => {
                let n = Vector::from(plane.normal());
                let f = plane.frame();
                let moved = crate::geometry::surface::Plane::from_frame(Frame::new(f.origin() + n * along, f.z(), f.x()));
                crate::make::face_on_domain(model, Surface::Plane(moved), lo.x - grow, hi.x + grow, lo.y - grow, hi.y + grow)?
            }
            Surface::Cylinder(c) => {
                // (Moved through its axis: the cylinder on the other side, turned round.)
                let r = c.radius() + along;
                let f = c.frame();
                let frame = if r >= 0.0 { *f } else { Frame::new(f.origin(), f.z(), -f.x()) };
                if r.abs() <= 1e-9 * size {
                    continue;
                }
                let turn = (grow / r.abs()).min(std::f64::consts::PI);
                let (u0, u1) = if hi.x - lo.x + 2.0 * turn >= std::f64::consts::TAU { (lo.x, lo.x + std::f64::consts::TAU) } else { (lo.x - turn, hi.x + turn) };
                crate::make::face_on_domain(model, Surface::Cylinder(Cylinder::new(frame, r.abs())), u0, u1, lo.y - grow, hi.y + grow)?
            }
            Surface::Cone(c) => {
                // (Moved along its normals: the cone of the same angle about the same axis,
                // its radius grown by the distance's share across the axis, shifted along it
                // by the rest; where that radius would be negative, the same cone from its
                // apex, its parameter along it shifted to match.)
                let (sin, cos) = c.semi_angle_sin_cos();
                let f = c.frame();
                let z = Vector::from(f.z());
                let r = c.radius() + along * cos;
                let origin = f.origin() - z * (along * sin);
                let (cone, shift) = if r >= 0.0 {
                    (Cone::new(Frame::new(origin, f.z(), f.x()), c.semi_angle(), r), 0.0)
                } else {
                    let apex = -r / sin;
                    (Cone::new(Frame::new(origin + z * (apex * cos), f.z(), f.x()), c.semi_angle(), 0.0), -apex)
                };
                // (Grown past its bounds each way: trimmed by another face, its offset need not
                // meet that face's within them.)
                // (All the way round: only dividing the cells further.)
                let grow = grow.max(4.0 * d);
                let (u0, u1) = (lo.x, lo.x + std::f64::consts::TAU);
                // (Not past its apex, on the side of it the face is: the moved cone's radius
                // there nothing. To the apex itself, where the sheet closes on a point: stopped
                // short of it, it leaves a hole there, and the cells inside it run on out.)
                let (mut v0, mut v1) = (lo.y - grow + shift, hi.y + grow + shift);
                let apex = -cone.radius() / cone.semi_angle_sin_cos().0;
                if 0.5 * (lo.y + hi.y) + shift >= apex {
                    v0 = v0.max(apex);
                } else {
                    v1 = v1.min(apex);
                }
                crate::make::face_on_domain(model, Surface::Cone(cone), u0, u1, v0, v1)?
            }
            Surface::Sphere(sphere) => {
                // (The whole sphere about the same centre, its radius grown by the distance;
                // through the centre, the same sphere turned round: only dividing the cells
                // further.)
                let r = sphere.radius() + along;
                if r.abs() <= 1e-9 * size {
                    continue;
                }
                let moved = Surface::Sphere(Sphere::new(*sphere.frame(), r.abs()));
                let half = 0.5 * std::f64::consts::PI;
                crate::make::face_on_domain(model, moved, 0.0, std::f64::consts::TAU, -half, half)?
            }
            Surface::Torus(torus) => {
                // (The whole torus, its tube's radius grown by the distance, as the sphere.)
                let r = torus.minor_radius() + along;
                if r.abs() <= 1e-9 * size {
                    continue;
                }
                let moved = Surface::Torus(Torus::new(*torus.frame(), torus.major_radius(), r.abs()));
                crate::make::face_on_domain(model, moved, 0.0, std::f64::consts::TAU, 0.0, std::f64::consts::TAU)?
            }
            _ => return Ok(None),
        };
        sheets.push(moved);
    }
    // Tubes about the edges whose kept faces' offsets part (convex edges going out, concave
    // ones going in), balls about their ends.
    let mut ends: Vec<(Point, usize)> = Vec::new();
    for edge in model.unique_subshapes(solid, NodeKind::Edge) {
        let using: Vec<usize> = kept.iter().copied().filter(|&i| model.unique_subshapes(faces[i].shape, NodeKind::Edge).iter().any(|e| e.node == edge.node)).collect();
        let [f, g] = using[..] else { continue };
        let Some((curve, c0, c1)) = model.edge_curve(edge) else { return Ok(None) };
        let (middle, along) = {
            let (p, d1) = curve.d1(0.5 * (c0 + c1));
            (p, d1 * (1.0 / d1.length().max(1e-300)))
        };
        let (nf, ng) = (faces[f].normal(middle), faces[g].normal(middle));
        if nf.cross(ng).length() < 1e-9 {
            continue;
        }
        // (Convex where the other face turns back from the way across the edge into one.)
        let across = nf.cross(along);
        let across = across * (1.0 / across.length());
        let step = 1e-4 * size;
        let into = if faces[f].distance(middle + across * step) < 1e-3 * step { across } else { -across };
        let convex = into.dot(ng) < 0.0;
        if convex != (sign > 0.0) {
            continue;
        }
        let (a, b) = (curve.point(c0), curve.point(c1));
        let length = a.distance(b);
        // (About a circular edge, a torus: the faces there surfaces of revolution about its
        // axis, the tube's cross-section the same all round.)
        if let crate::geometry::curve::Curve::Circle(circle) = &curve {
            let f = circle.frame();
            let z = Vector::from(f.z());
            let radial = |p: Point| {
                let w = p - f.origin();
                let w = w - z * w.dot(z);
                w * (1.0 / w.length().max(1e-300))
            };
            let e = radial(middle);
            let angle_of = |n: Vector| n.dot(z).atan2(n.dot(e));
            let (x, y) = (nf * sign, ng * sign);
            let (vx, vy) = (angle_of(x), angle_of(y));
            let mut sweep = vy - vx;
            if sweep > std::f64::consts::PI {
                sweep -= std::f64::consts::TAU;
            } else if sweep < -std::f64::consts::PI {
                sweep += std::f64::consts::TAU;
            }
            let (v0, v1) = if sweep >= 0.0 { (vx, vx + sweep) } else { (vx + sweep, vx) };
            let torus = Surface::Torus(Torus::new(*f, circle.radius(), d));
            sheets.push(crate::make::face_on_domain(model, torus, c0.min(c1), c0.max(c1), v0, v1)?);
            if length > 1e-9 * size {
                for p in [a, b] {
                    match ends.iter_mut().find(|(q, _)| q.distance(p) <= 1e-9 * size) {
                        Some((_, n)) => *n += 1,
                        None => ends.push((p, 1)),
                    }
                }
            }
            continue;
        }
        // (About any other edge, the whole pipe round it: only dividing the cells further
        // where the faces' offsets do not part.)
        if !matches!(curve, crate::geometry::curve::Curve::Line(_)) && length > 1e-9 * size {
            let (start, tangent) = curve.d1(c0);
            let Some(t) = Direction::try_from_vector(tangent) else { return Ok(None) };
            let ring = crate::make::edge_from_curve(model, crate::geometry::curve::Curve::Circle(crate::geometry::curve::Circle::new(Frame::from_z(start, t), d)))?;
            let along = crate::make::edge_from_curve_range(model, curve.clone(), c0, c1)?;
            let spine = crate::make::wire(model, &[along])?;
            let Ok(pipe) = crate::sweep::pipe(model, ring, spine) else { return Ok(None) };
            sheets.extend(model.unique_subshapes(pipe.shape, NodeKind::Face));
            for p in [a, b] {
                match ends.iter_mut().find(|(q, _)| q.distance(p) <= 1e-9 * size) {
                    Some((_, n)) => *n += 1,
                    None => ends.push((p, 1)),
                }
            }
            continue;
        }
        if length <= 1e-12 * size {
            return Ok(None);
        }
        // (The tube's angle running from one face's offset to the other's.)
        let (x, y) = (nf * sign, ng * sign);
        let mut t = (b - a) * (1.0 / length);
        let mut origin = a;
        if t.cross(x).dot(y) < 0.0 {
            t = -t;
            origin = b;
        }
        let angle = x.dot(y).clamp(-1.0, 1.0).acos();
        let (Some(z), Some(xd)) = (Direction::try_from_vector(t), Direction::try_from_vector(x)) else { continue };
        let cylinder = Surface::Cylinder(Cylinder::new(Frame::new(origin, z, xd), d));
        sheets.push(crate::make::face_on_domain(model, cylinder, 0.0, angle, 0.0, length)?);
        for p in [a, b] {
            match ends.iter_mut().find(|(q, _)| q.distance(p) <= 1e-9 * size) {
                Some((_, n)) => *n += 1,
                None => ends.push((p, 1)),
            }
        }
    }
    // (A ball all of whose points lie nearer other faces bounds nothing: left out, as the
    // cell it would be, touching the rest only along the tubes and faces it meets, is not
    // always found as a cavity of the cell around it.)
    let bounds = |p: Point| {
        (0..64).any(|i| {
            let z = 1.0 - (2.0 * i as f64 + 1.0) / 64.0;
            let phi = i as f64 * std::f64::consts::PI * (3.0 - 5f64.sqrt());
            let r = (1.0 - z * z).sqrt();
            let q = p + Vector::new(r * phi.cos(), r * phi.sin(), z) * d;
            kept.iter().all(|&i| faces[i].distance(q) >= d * (1.0 - 1e-6))
        })
    };
    for p in ends.into_iter().filter(|&(p, n)| n > 1 && bounds(p)).map(|(p, _)| p) {
        let ball = crate::primitives::Sphere::new(d).with_frame(Frame::from_z(p, Direction::Z)).build(model)?;
        sheets.extend(model.unique_subshapes(ball, NodeKind::Face));
    }
    sheets.extend(strips);
    // A hollowed solid's walls end on its faces.
    if matches!(made, Made::Thick { .. }) {
        sheets.push(solid);
    }
    let cells = crate::boolean::enclosed_solids(model, &sheets)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    let mut chosen = Vec::new();
    for &cell in &cells {
        let Some(p) = crate::boolean::point_inside(model, cell)? else { return Ok(None) };
        // (Over a removed face met tangentially all round, within the distance on the side the
        // walls grow: the opening, as far as the strips round it (points there may be within
        // the distance of the edges ending at its corners).)
        let opening = openings.iter().any(|&r| {
            let (q, on) = faces[r].foot(p);
            let foot = faces[r].surface.point(q.x, q.y);
            on && foot.distance(p) < d && (p - foot).dot(faces[r].normal(foot)) * sign > 0.0
        });
        let within = !opening && kept.iter().map(|&i| faces[i].distance(p)).fold(f64::INFINITY, f64::min) < d;
        // (A point on the solid's boundary is within the distance: it lies in an offset
        // outwards and out of one inwards; a hollowed solid's cells are bounded by the
        // solid's faces, and lie on one side of them.)
        let keep = match (inside.locate(p), &made, sign > 0.0) {
            (PointLocation::Inside, Made::Offset, grown) => grown || !within,
            (PointLocation::Outside, Made::Offset, grown) => grown && within,
            (PointLocation::Inside, Made::Thick { .. }, grown) => !grown && within,
            (PointLocation::Outside, Made::Thick { .. }, grown) => grown && within,
            (_, Made::Offset, grown) => grown,
            _ => return Ok(None),
        };
        if keep {
            chosen.push(cell);
        }
    }
    if chosen.is_empty() {
        return Ok(None);
    }
    // (Pieces of one moved face may come from its copy and the box grown past it.)
    super::regularize::joined(model, &chosen, false).map(Some)
}

/// Where a kept face meets a removed one tangentially (a fillet running into a removed top),
/// the wall ends square to the edge, as OCCT's do: on the strip the edge sweeps along the
/// removed face's normal (the kept face's offset reaches no further than that strip, where it
/// is the removed face's offset, which is not there). The points within the distance of the
/// edge itself, past the strip, are not the wall's: the edge is no longer measured from. The
/// strips, a little longer and wider than the walls, as sheets dividing the cells; whether
/// every edge between a kept face and a removed one is such an edge; and the removed faces
/// (by index) all of whose edges with kept faces are.
fn strips(model: &mut Model, solid: Shape, faces: &mut [Bounded], removed: &[NodeId], d: f64) -> Result<(Vec<Shape>, bool, Vec<usize>)> {
    use crate::geometry::curve::Curve;
    let mut out = Vec::new();
    let mut all_tangent = true;
    // (Removed faces with an edge met at an angle.)
    let mut angled: Vec<usize> = Vec::new();
    if removed.is_empty() {
        return Ok((out, all_tangent, Vec::new()));
    }
    for edge in model.unique_subshapes(solid, NodeKind::Edge) {
        let using: Vec<usize> = (0..faces.len()).filter(|&i| faces[i].edge_nodes.contains(&edge.node)).collect();
        let [f, g] = using[..] else { continue };
        let (k, r) = match (removed.contains(&faces[f].shape.node), removed.contains(&faces[g].shape.node)) {
            (false, true) => (f, g),
            (true, false) => (g, f),
            _ => continue,
        };
        let Some((curve, c0, c1)) = model.edge_curve(edge) else { continue };
        let middle = curve.point(0.5 * (c0 + c1));
        let (nk, nr) = (faces[k].normal(middle), faces[r].normal(middle));
        if nk.cross(nr).length() > 1e-6 {
            all_tangent = false;
            angled.push(r);
            continue;
        }
        let Some(n) = Direction::try_from_vector(nr) else { continue };
        if let Some(i) = faces[k].edge_nodes.iter().position(|&e| e == edge.node) {
            faces[k].measured[i] = false;
        }
        let reach = 3.0 * d;
        let sheet = match &curve {
            Curve::Line(_) => {
                let (a, b) = (curve.point(c0), curve.point(c1));
                let length = a.distance(b);
                let Some(t) = Direction::try_from_vector(b - a) else { continue };
                let Some(normal) = Direction::try_from_vector(Vector::from(t).cross(Vector::from(n))) else { continue };
                let plane = crate::geometry::surface::Plane::from_frame(Frame::new(a, normal, t));
                crate::make::face_on_domain(model, Surface::Plane(plane), -reach, length + reach, -reach, reach)?
            }
            Curve::Circle(circle) if Vector::from(circle.frame().z()).cross(Vector::from(n)).length() < 1e-9 => {
                let cylinder = Cylinder::new(*circle.frame(), circle.radius());
                crate::make::face_on_domain(model, Surface::Cylinder(cylinder), 0.0, std::f64::consts::TAU, -reach, reach)?
            }
            _ => {
                let swept = crate::geometry::surface::ExtrusionSurface::new(curve.clone(), n);
                let moved = crate::make::face_on_domain(model, Surface::Extrusion(swept), c0.min(c1), c0.max(c1), 0.0, 2.0 * reach)?;
                crate::modify::transform(model, moved, &Transform::translation(Vector::from(n) * -reach), Default::default())?
            }
        };
        out.push(sheet);
    }
    let openings = (0..faces.len()).filter(|&i| removed.contains(&faces[i].shape.node) && !angled.contains(&i)).collect();
    Ok((out, all_tangent, openings))
}
