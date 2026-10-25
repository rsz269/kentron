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
//! whole of its parameters' box. Faces on planes and cylinders (a cylinder moved through its
//! axis turned round it), tubes about straight edges.

use crate::boolean::{PointClassifier, PointLocation};
use crate::error::Result;
use crate::geometry::surface::{Cylinder, Surface};
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
    u_period: Option<f64>,
}

impl Bounded {
    /// The parameters of the foot of the perpendicular from `p` to the surface, in the
    /// period of the face's loops where `p` falls in them.
    fn foot(&self, p: Point) -> (Point2, bool) {
        let (u, v) = self.surface.parameters(p).unwrap_or_else(|| {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &self.surface, None);
            (u, v)
        });
        let shifts: &[f64] = if self.u_period.is_some() { &[0.0, -1.0, 1.0, -2.0, 2.0] } else { &[0.0] };
        for &k in shifts {
            let q = Point2::new(u + k * self.u_period.unwrap_or(0.0), v);
            if winding(&self.loops, q) != 0 {
                return (q, true);
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
        self.edges.iter().map(|(c, a, b)| crate::geometry::extrema::nearest_on_curve(p, c, *a, *b).distance).fold(f64::INFINITY, f64::min)
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
/// solid of planar and cylindrical faces whose offsets part along straight edges only, or
/// its cells cannot be told apart.
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
    // (Walls grown outwards with faces removed end on the removed faces' planes, not
    // where the distance takes them.)
    if distance > 0.0 && !removed.is_empty() {
        return Ok(None);
    }
    let mut faces = Vec::new();
    for face in model.unique_subshapes(solid, NodeKind::Face) {
        let Some(surface @ (Surface::Plane(_) | Surface::Cylinder(_))) = model.surface(face).cloned() else { return Ok(None) };
        let Some(uv) = face_loops(model, face) else { return Ok(None) };
        let loops = uv.iter().map(|l| l.iter().flat_map(|u| (0..16).map(move |i| u.at(i as f64 / 16.0))).collect()).collect();
        let mut edges = Vec::new();
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            let Some(curve) = model.edge_curve(edge) else { return Ok(None) };
            edges.push(curve);
        }
        let u_period = surface.u_period();
        faces.push(Bounded { shape: face, surface, reversed: face.orientation == Orientation::Reversed, loops, edges, u_period });
    }
    let kept: Vec<usize> = (0..faces.len()).filter(|&i| !removed.contains(&faces[i].shape.node)).collect();
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
        let grow = if beside_removed { (4.0 * d).max(0.25 * size) } else { 0.0 };
        let moved = match &face.surface {
            Surface::Plane(plane) if !beside_removed => {
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
        // (Tubes about straight edges only.)
        let (a, b) = (curve.point(c0), curve.point(c1));
        let length = a.distance(b);
        if !matches!(curve, crate::geometry::curve::Curve::Line(_)) || length <= 1e-12 * size {
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
    // A hollowed solid's walls end on its faces.
    if matches!(made, Made::Thick { .. }) {
        sheets.push(solid);
    }
    let cells = crate::boolean::enclosed_solids(model, &sheets)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    let mut chosen = Vec::new();
    for &cell in &cells {
        let Some(p) = crate::boolean::point_inside(model, cell)? else { return Ok(None) };
        let within = kept.iter().map(|&i| faces[i].distance(p)).fold(f64::INFINITY, f64::min) < d;
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
