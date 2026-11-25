//! Draft angles: faces tilted against a pull direction, as moulded parts need to leave
//! their mould.
//!
//! A drafted face keeps its place on a neutral plane and tilts by the draft angle relative
//! to the pull direction: a plane turns about its line on the neutral plane, a cylinder or
//! cone about the pull direction becomes the cone through its circle on the neutral plane.
//! Faces meeting a drafted one tangentially are drafted with it (those that can be). The shape keeps its
//! topology: its edges and vertices are found again where the new surfaces meet.
//!
//! Faces moved rigidly ([`move_faces`]) are the same direct edit: the moved faces' surfaces
//! carried by the motion, the faces round them kept on theirs and met again.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::surface::{Cone, Cylinder, Plane, Surface};
use crate::geometry::ParametricSurface;
use crate::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use crate::tolerance::{ANGULAR, LINEAR};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

/// A face to draft: by `angle` (radians) to the pull direction, staying in place where it
/// meets `neutral`.
#[doc(alias = "BRepOffsetAPI_DraftAngle::Add")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DraftedFace {
    pub face: Shape,
    pub angle: f64,
    pub neutral: Plane,
}

/// `shape` with `faces` drafted against the pull `direction` (and the faces tangent to them
/// drafted alike). Faces with an angle below 1e-4 stay as they are.
#[doc(alias = "BRepOffsetAPI_DraftAngle")]
#[doc(alias = "depouille")]
pub fn draft(model: &mut Model, shape: Shape, direction: Direction, faces: &[DraftedFace]) -> Result<Shape> {
    match drafted(model, shape, direction, faces, true) {
        Ok(drafted) => Ok(drafted),
        Err(e) => cut_by_drafted_planes(model, shape, direction, faces)?.ok_or(e),
    }
}

/// Planar faces drafted so far that they close over the part, their tilted planes meeting
/// one another or another face before the part ends (a ridge, an apex), which the draft
/// keeping the part's topology cannot follow: where each face's plane bounds the whole part
/// (the part lies on its inner side) and the draft moves the face inwards all over it, the
/// draft takes away just what lies beyond the tilted planes, the faces round them trimmed
/// where they meet them: the part cut by the tilted planes. `None` where that does not hold
/// (a face not a plane, or drafted with faces tangent to it; a part with curved faces, which
/// may bulge past a plane between its vertices; a plane not bounding the part; a face moved
/// outwards somewhere, as below its neutral plane).
fn cut_by_drafted_planes(model: &mut Model, shape: Shape, direction: Direction, faces: &[DraftedFace]) -> Result<Option<Shape>> {
    let faces: Vec<&DraftedFace> = faces.iter().filter(|f| f.angle.abs() > 1e-4).collect();
    let neighbours = tangent_neighbours(model, shape);
    if faces.is_empty() || faces.iter().any(|f| neighbours.get(&f.face.node).is_some_and(|n| !n.is_empty())) {
        return Ok(None);
    }
    let placed = |model: &Model, f: Shape| model.surface(f).map(|s| s.transformed(&f.placement));
    let planar = |s: &Surface| {
        let mut s = s;
        while let Surface::Trimmed(t) = s {
            s = t.basis();
        }
        matches!(s, Surface::Plane(_))
    };
    if !model.unique_subshapes(shape, NodeKind::Face).into_iter().all(|f| placed(model, f).is_some_and(|s| planar(&s))) {
        return Ok(None);
    }
    let corners: Vec<Point> = model.unique_subshapes(shape, NodeKind::Vertex).into_iter().map(|v| model.point(v)).collect();
    // A face's outward normal on `surface` (a plane, the face's or its drafted one) and a point
    // of it.
    let outward = |surface: &Surface, face: Shape| {
        let (o, du, dv) = surface.d1(0.0, 0.0);
        let n = du.cross(dv).normalized();
        (o, if face.orientation == Orientation::Reversed { -n } else { n })
    };
    let mut tools = Vec::new();
    for f in faces {
        let Some(surface) = placed(model, f.face) else { return Ok(None) };
        let (o, n) = outward(&surface, f.face);
        if corners.iter().any(|&p| (p - o).dot(n) > LINEAR) {
            return Ok(None);
        }
        let Some(tilted) = drafted_surface(&surface, f.face.orientation == Orientation::Reversed, direction, f.angle, &f.neutral) else { return Ok(None) };
        let (o, n) = outward(&tilted, f.face);
        let own: Vec<Point> = model.unique_subshapes(f.face, NodeKind::Vertex).into_iter().map(|v| model.point(v)).collect();
        // (Moved inwards: the face as it was lies on or beyond its tilted plane.)
        if own.iter().any(|&p| (p - o).dot(n) < -LINEAR) {
            return Ok(None);
        }
        tools.push((tilted, o - n));
    }
    let mut cut = shape;
    for (plane, inside) in tools {
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, plane, LINEAR)?;
        let half = crate::primitives::half_space(model, face, inside)?;
        cut = crate::boolean::intersection(model, cut, half)?;
    }
    Ok(Some(cut))
}

/// `draft`; `ordered`: a straight edge the draft turns round (its faces' new planes meeting
/// past its ends, as where drafted sides close over a face) fails it, as it fails OCCT's
/// Draft_Modification (Draft_EdgeRecomputation). A drafted prism, longer than the part of it
/// a feature keeps, lets its far end turn round.
pub(crate) fn drafted(model: &mut Model, shape: Shape, direction: Direction, faces: &[DraftedFace], ordered: bool) -> Result<Shape> {
    // Each face to draft (by node), with its settings: the given faces, then those tangent
    // to them.
    let mut drafted: HashMap<NodeId, (f64, Plane)> = HashMap::new();
    // A face with no angle to speak of is left as it is (as OCCT does).
    let mut queue: Vec<(Shape, f64, Plane)> = faces.iter().filter(|f| f.angle.abs() > 1e-4).map(|f| (f.face, f.angle, f.neutral)).collect();
    let neighbours = tangent_neighbours(model, shape);
    while let Some((face, angle, neutral)) = queue.pop() {
        if drafted.contains_key(&face.node) {
            continue;
        }
        drafted.insert(face.node, (angle, neutral));
        for &next in neighbours.get(&face.node).into_iter().flatten() {
            // Only neighbours that can be drafted: OCCT follows edges marked tangent, which
            // kentron does not record; a face it could not draft would make it fail.
            let draftable = model.surface(next).is_some_and(|s| {
                drafted_surface(&s.transformed(&next.placement), next.orientation == Orientation::Reversed, direction, angle, &neutral).is_some()
            });
            if draftable {
                queue.push((next, angle, neutral));
            }
        }
    }
    let new_surface = |face: Shape, surface: &Surface| -> Result<Surface> {
        match drafted.get(&face.node) {
            None => Ok(surface.clone()),
            Some(&(angle, neutral)) => {
                drafted_surface(surface, face.orientation == Orientation::Reversed, direction, angle, &neutral)
                    .ok_or_else(|| Error::Failed("a face that cannot be drafted (not a plane, or a cylinder or cone about the pull direction)".into()))
            }
        }
    };
    let drafted = crate::offset::shape::with_surfaces(model, shape, &new_surface, ordered)?;
    // (A loop the draft moves over a seam: the seam moved off it.)
    crate::heal::move_crossed_seams(model, drafted)?;
    Ok(drafted)
}

/// `shape` with `faces` moved rigidly by `motion` (`p ↦ R·p + t`, `R` a rotation), as a
/// direct edit: the moved faces' surfaces carried by the motion, the faces round them kept
/// on their surfaces, extended or trimmed to meet them, the edges and vertices found again
/// where the surfaces meet (as Draft_Modification finds them for drafted faces). A hole's
/// face moved across a plate moves the hole; a plane moved along itself stays the plane it
/// was.
///
/// `faces` must name one face of `shape` or more, and `motion` be rigid (InvalidInput
/// otherwise); where the faces round the moved ones cannot meet them, it fails.
#[doc(alias = "BRepTools_Modifier")]
#[doc(alias = "BRepTools_Modification")]
pub fn move_faces(model: &mut Model, shape: Shape, faces: &[Shape], motion: &Transform) -> Result<Shape> {
    if faces.is_empty() {
        return Err(Error::InvalidInput("no faces to move".into()));
    }
    let own: Vec<NodeId> = model.unique_subshapes(shape, NodeKind::Face).into_iter().map(|f| f.node).collect();
    if let Some(f) = faces.iter().find(|f| f.kind() != NodeKind::Face || !own.contains(&f.node)) {
        return Err(Error::InvalidInput(format!("{:?} is not a face of the shape", f.node)));
    }
    if !is_rigid(motion) {
        return Err(Error::InvalidInput("a face move is rigid: a rotation and a translation".into()));
    }
    let moved: Vec<NodeId> = faces.iter().map(|f| f.node).collect();
    let new_surface = |face: Shape, surface: &Surface| -> Result<Surface> {
        if !moved.contains(&face.node) {
            return Ok(surface.clone());
        }
        let carried = surface.transformed(motion);
        // A plane moved along itself is the same plane (as Draft_Modification keeps a face
        // whose draft turns it by no angle).
        if let (Surface::Plane(old), Surface::Plane(new)) = (surface, &carried) {
            let along = Vector::from(old.normal()).dot(Vector::from(new.normal())) >= 1.0 - ANGULAR;
            if along && (new.origin() - old.origin()).dot(Vector::from(old.normal())).abs() <= LINEAR {
                return Ok(surface.clone());
            }
        }
        Ok(carried)
    };
    // Each edge between a moved face and one kept is the section of their surfaces: where
    // they have none (a plane turned parallel to the plane it met), the edit fails, as
    // Draft_Modification fails an edge it cannot recompute (Draft_EdgeRecomputation).
    let faces_of_edge = faces_of_edges(model, shape);
    let placed = |f: Shape| model.surface(f).map(|s| s.transformed(&f.placement));
    for pair in faces_of_edge.values() {
        let [a, b] = pair[..] else { continue };
        if moved.contains(&a.node) == moved.contains(&b.node) {
            continue;
        }
        let (Some(sa), Some(sb)) = (placed(a), placed(b)) else { continue };
        let (na, nb) = (new_surface(a, &sa)?, new_surface(b, &sb)?);
        use crate::geometry::intersect::{quadrics, QuadricIntersection as Q};
        if matches!(quadrics(&na, &nb, LINEAR), Q::Empty | Q::Points(_)) {
            return Err(Error::Failed("a face round the moved ones cannot reach them (their surfaces do not meet)".into()));
        }
    }
    // (Where the surfaces round a moved face do not reach it, no vertex or edge is found:
    // the edit fails, as Draft_Modification fails it.)
    let moved = crate::offset::shape::with_surfaces(model, shape, &new_surface, true).map_err(|e| match e {
        Error::InvalidInput(m) | Error::Failed(m) => Error::Failed(format!("the faces round the moved ones do not meet them: {m}")),
    })?;
    crate::heal::move_crossed_seams(model, moved)?;
    Ok(moved)
}

/// Whether `t` moves rigidly: no scale, no reflection, and its matrix a rotation (`RᵀR = I`
/// to the angular tolerance).
fn is_rigid(t: &Transform) -> bool {
    let rows = t.rows();
    let column = |c: usize| Vector::new(rows[0][c], rows[1][c], rows[2][c]);
    let orthonormal = (0..3).all(|i| (0..3).all(|j| (column(i).dot(column(j)) - if i == j { 1.0 } else { 0.0 }).abs() <= ANGULAR));
    orthonormal && !t.is_negative()
}

/// For each edge of `shape`, the faces using it (each once).
fn faces_of_edges(model: &Model, shape: Shape) -> HashMap<NodeId, Vec<Shape>> {
    let mut faces_of_edge: HashMap<NodeId, Vec<Shape>> = HashMap::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            let list = faces_of_edge.entry(edge.node).or_default();
            if !list.iter().any(|f| f.node == face.node) {
                list.push(face);
            }
        }
    }
    faces_of_edge
}

/// For each face of `shape`, the faces meeting it tangentially along an edge.
fn tangent_neighbours(model: &Model, shape: Shape) -> HashMap<NodeId, Vec<Shape>> {
    let faces_of_edge = faces_of_edges(model, shape);
    let mut out: HashMap<NodeId, Vec<Shape>> = HashMap::new();
    for (edge, faces) in &faces_of_edge {
        let [a, b] = faces[..] else { continue };
        let edge_shape = Shape::from(*edge);
        let Some((curve, first, last)) = model.edge_curve(edge_shape) else { continue };
        use crate::geometry::ParametricCurve;
        let p = curve.point(0.5 * (first + last));
        let normal = |face: Shape| -> Option<Vector> {
            let surface = model.surface(face)?.transformed(&face.placement);
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
            let (_, du, dv) = surface.d1(u, v);
            let n = du.cross(dv);
            let n = n * (1.0 / n.length());
            Some(if face.orientation == Orientation::Reversed { -n } else { n })
        };
        let (Some(na), Some(nb)) = (normal(a), normal(b)) else { continue };
        // Tangent (G1): the same outward normal.
        if na.dot(nb) >= 1.0 - 1e-6 {
            out.entry(a.node).or_default().push(b);
            out.entry(b.node).or_default().push(a);
        }
    }
    out
}

/// The surface of a face drafted by `angle` against the pull `direction`, staying in place
/// where it meets `neutral` (`reversed`: the face uses its surface reversed): a plane turned
/// about its line on the neutral plane so that its outward normal makes `sin angle` with the
/// pull; a cylinder or cone about the pull direction the cone through its circle on the
/// neutral plane, leaning by `angle` (a cone back to a cylinder at no angle). The frame of
/// the result is the input's (a plane's turned with it). `None` where the face cannot be
/// drafted: another kind of surface, a plane parallel to the neutral plane or hinged along
/// the pull, an angle no plane through the hinge makes, a cylinder or cone not about the
/// pull direction or the neutral plane's normal.
fn drafted_surface(surface: &Surface, reversed: bool, direction: Direction, angle: f64, neutral: &Plane) -> Option<Surface> {
    let mut surface = surface;
    while let Surface::Trimmed(t) = surface {
        surface = t.basis();
    }
    let d = Vector::from(direction);
    let m = Vector::from(neutral.normal());
    // The face's outward side: +1 along the frame's normal (plane) or away from the axis.
    let outward = |right_handed: bool| if right_handed != reversed { 1.0 } else { -1.0 };
    // The signed height along the axis `z` from `origin` to the neutral plane.
    let height = |origin: Point, z: Vector| {
        let zm = z.dot(m);
        (zm.abs() > ANGULAR).then(|| (neutral.origin() - origin).dot(m) / zm)
    };
    // A cylinder's or cone's drafted generator: its semi-angle and radius at the frame's
    // origin, given the radius on the neutral plane.
    let leaning = |sigma: f64, z: Vector, h: f64, on_neutral: f64| {
        let s = if d.dot(z) < 0.0 { -sigma * angle } else { sigma * angle };
        let r = on_neutral + h * s.tan();
        if r >= 0.0 { (-s, r) } else { (s, -r) }
    };
    let about_pull = |z: Vector| d.dot(m).abs() > 1.0 - ANGULAR && d.dot(z).abs() > 1.0 - ANGULAR;
    match surface {
        Surface::Plane(plane) => {
            let frame = plane.frame();
            let n = Vector::from(frame.z());
            let o = frame.origin();
            let hinge = n.cross(m);
            if !(hinge.length() > 0.0) {
                return None;
            }
            let l = hinge.normalized();
            // A point of the hinge: from the origin, within the plane, square to the hinge.
            let across = l.cross(n);
            let at = o + across * ((neutral.origin() - o).dot(m) / across.dot(m));
            let a = d.dot(l);
            if a.abs() > 1.0 - ANGULAR {
                return None;
            }
            let rho = (1.0 - a * a).sqrt();
            if rho <= angle.sin().abs() {
                return None;
            }
            // Turning σn by θ about l, its component along d is c cos θ − b sin θ =
            // ρ cos(θ + φ); set to sin α.
            let sigma = outward(frame.is_right_handed());
            let b = sigma * d.dot(n.cross(l));
            let c = sigma * d.dot(n);
            let phi = (b / rho).atan2(c / rho);
            let theta0 = (angle.sin() / rho).clamp(-1.0, 1.0).acos();
            let mut theta = theta0 - phi;
            if theta.cos() < 0.0 {
                theta = -theta0 - phi;
            }
            // (Brought within a half turn by half turns: the same plane, its frame's normal
            // possibly flipped.)
            while theta.abs() > std::f64::consts::PI {
                theta -= std::f64::consts::PI * theta.signum();
            }
            if theta.abs() <= ANGULAR {
                return Some(surface.clone());
            }
            let turn = Transform::rotation(Axis { origin: at, direction: Direction::from_vector(l) }, theta);
            Some(Surface::Plane(plane.transformed(&turn)))
        }
        Surface::Cylinder(cylinder) => {
            let frame = *cylinder.frame();
            let z = Vector::from(frame.z());
            if !about_pull(z) {
                return None;
            }
            if angle.abs() <= ANGULAR {
                return Some(surface.clone());
            }
            let h = height(frame.origin(), z)?;
            let (semi, r) = leaning(outward(frame.is_right_handed()), z, h, cylinder.radius());
            Some(Surface::Cone(Cone::new(frame, semi, r)))
        }
        Surface::Cone(cone) => {
            let frame = *cone.frame();
            let z = Vector::from(frame.z());
            if !about_pull(z) {
                return None;
            }
            let h = height(frame.origin(), z)?;
            let on_neutral = cone.radius() + h * cone.semi_angle().tan();
            if on_neutral <= 0.0 {
                return None;
            }
            if angle.abs() <= ANGULAR {
                return Some(Surface::Cylinder(Cylinder::new(frame, on_neutral)));
            }
            let (semi, r) = leaning(outward(frame.is_right_handed()), z, h, on_neutral);
            if (semi - cone.semi_angle()).abs() <= ANGULAR {
                return Some(surface.clone());
            }
            Some(Surface::Cone(Cone::new(frame, semi, r)))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Cylinder as CylinderPrimitive};
    use crate::properties::volume;

    fn valid(model: &Model, s: Shape) {
        let problems: Vec<String> = crate::validity::check(model, s).errors().map(|p| p.description.clone()).collect();
        assert!(problems.is_empty(), "{problems:?}");
    }

    /// The face of `shape` on a plane through `x`-coordinate `x` (parallel to YZ).
    fn face_at_x(model: &Model, shape: Shape, x: f64) -> Shape {
        model
            .unique_subshapes(shape, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Plane(p)) if (p.origin().x - x).abs() < 1e-9 && p.normal().x().abs() > 0.5))
            .expect("a face at x")
    }

    /// The face of `shape` on a plane through `y`-coordinate `y` (parallel to XZ).
    fn face_at_y(model: &Model, shape: Shape, y: f64) -> Shape {
        model
            .unique_subshapes(shape, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Plane(p)) if (p.origin().y - y).abs() < 1e-9 && p.normal().y().abs() > 0.5))
            .expect("a face at y")
    }

    /// Sides drafted so far that they meet before the top: the part closes along a ridge
    /// (two opposite sides) or at a point (all four), the top gone. Its volume, the
    /// integral of the section `(w − nx·τz)(d − ny·τz)` up to where it closes.
    #[test]
    fn sides_drafted_past_meeting_close_the_part() {
        let (w, d, h) = (4.0, 6.0, 5.0);
        let angle = 0.6f64;
        let tau = angle.tan();
        let integral = |nx: f64, ny: f64| {
            let top = [(nx, w), (ny, d)].iter().filter(|(n, _)| *n > 0.0).map(|(n, l)| l / (n * tau)).fold(h, f64::min);
            w * d * top - (w * ny + d * nx) * tau * top * top / 2.0 + nx * ny * tau * tau * top.powi(3) / 3.0
        };
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        for (sides, nx, ny) in [(vec![0, 1], 2.0, 0.0), (vec![0, 1, 2, 3], 2.0, 2.0), (vec![0, 2], 1.0, 1.0)] {
            let mut model = Model::new();
            let b = Cuboid::new(w, d, h).build(&mut model).unwrap();
            let faces: Vec<DraftedFace> = sides
                .iter()
                .map(|&i| {
                    let face = match i {
                        0 => face_at_x(&model, b, 0.0),
                        1 => face_at_x(&model, b, w),
                        2 => face_at_y(&model, b, 0.0),
                        _ => face_at_y(&model, b, d),
                    };
                    DraftedFace { face, angle, neutral: bottom }
                })
                .collect();
            let drafted = draft(&mut model, b, Direction::Z, &faces).unwrap();
            let expected = integral(nx, ny);
            let got = volume(&model, drafted).mass;
            assert!((got - expected).abs() < 1e-9 * expected, "{sides:?}: {got} against {expected}");
            valid(&model, drafted);
        }
    }

    #[test]
    fn a_box_side_leans_in() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let side = face_at_x(&model, b, 10.0);
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        let drafted = draft(&mut model, b, Direction::Z, &[DraftedFace { face: side, angle: 5f64.to_radians(), neutral: bottom }]).unwrap();
        // The side stays at x = 10 on the bottom and leans in going up (OCCT: 5212.6).
        let expected = 20.0 * (300.0 - 0.5 * 30.0 * 30.0 * 5f64.to_radians().tan());
        assert!((volume(&model, drafted).mass - expected).abs() < 1e-6 * expected, "{}", volume(&model, drafted).mass);
        valid(&model, drafted);
    }

    /// A drafted box recorded: each of its faces makes one face of the result, the drafted
    /// side the leaning one.
    #[test]
    fn drafts_record_their_faces() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let side = face_at_x(&model, b, 10.0);
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        crate::boolean::record_begin(&model);
        let drafted = draft(&mut model, b, Direction::Z, &[DraftedFace { face: side, angle: 5f64.to_radians(), neutral: bottom }]).unwrap();
        let history = crate::boolean::record_end(&model, drafted);
        for face in model.unique_subshapes(b, NodeKind::Face) {
            assert_eq!(history.iter().filter(|(from, _)| *from == face.node).count(), 1);
        }
        let &(_, leaning) = history.iter().find(|(from, _)| *from == side.node).unwrap();
        let area = crate::properties::area(&model, Shape::from(leaning)).mass;
        assert!((area - 20.0 * 30.0 / 5f64.to_radians().cos()).abs() < 1e-6 * area, "{area}");
    }

    #[test]
    fn a_cylinder_becomes_a_cone() {
        let mut model = Model::new();
        let c = CylinderPrimitive::new(5.0, 10.0).build(&mut model).unwrap();
        let side = model
            .unique_subshapes(c, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Cylinder(_))))
            .unwrap();
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        let drafted = draft(&mut model, c, Direction::Z, &[DraftedFace { face: side, angle: 5f64.to_radians(), neutral: bottom }]).unwrap();
        // A frustum from radius 5 at the bottom (OCCT: 655.987).
        let top = 5.0 - 10.0 * 5f64.to_radians().tan();
        let expected = std::f64::consts::PI * 10.0 / 3.0 * (25.0 + 5.0 * top + top * top);
        assert!((volume(&model, drafted).mass - expected).abs() < 1e-6 * expected, "{}", volume(&model, drafted).mass);
        valid(&model, drafted);
    }

    /// OCCT's draft/angle J5: a cylinder lying along X (its seam turned toward -Y, +Z)
    /// with a small upright one fused on; drafting the small one widens it below its top,
    /// and its loop on the large one runs over that one's seam, which is moved off it.
    #[test]
    fn a_loop_drafted_over_a_seam_moves_it() {
        use crate::math::Frame;
        let mut model = Model::new();
        let lying = Frame::new(Point::ORIGIN, Direction::X, Direction::from_vector(Vector::new(0.0, -0.5, 1.0)));
        let large = CylinderPrimitive::new(10.0, 20.0).with_frame(lying).build(&mut model).unwrap();
        let upright = Frame::from_z(Point::new(10.0, -0.2, 5.0), Direction::Z);
        let small = CylinderPrimitive::new(4.0, 10.0).with_frame(upright).build(&mut model).unwrap();
        let fused = crate::boolean::union(&mut model, large, small).unwrap();
        let side = model
            .unique_subshapes(fused, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Cylinder(c)) if (c.radius() - 4.0).abs() < 1e-9))
            .unwrap();
        let top = Plane::new(Point::new(0.0, 0.0, 15.0), Direction::Z);
        let drafted = draft(&mut model, fused, Direction::Z, &[DraftedFace { face: side, angle: 4f64.to_radians(), neutral: top }]).unwrap();
        valid(&model, drafted);
        // (OCCT: 2018.64.)
        let area = crate::properties::area(&model, drafted).mass;
        assert!((area - 2018.64).abs() < 1e-3 * 2018.64, "{area}");
    }

    /// The outline through `corners` in the plane z = 0 extruded `h` up, its sides drafted
    /// `degrees` about its floor, inward as they rise.
    fn sides_drafted(model: &mut Model, corners: &[(f64, f64)], h: f64, degrees: f64) -> Result<Shape> {
        let n = corners.len();
        let edges: Vec<Shape> = (0..n)
            .map(|i| {
                let ((x0, y0), (x1, y1)) = (corners[i], corners[(i + 1) % n]);
                crate::make::edge_between(model, Point::new(x0, y0, 0.0), Point::new(x1, y1, 0.0)).unwrap()
            })
            .collect();
        let wire = crate::make::wire(model, &edges).unwrap();
        let floor = crate::make::planar_face(model, wire).unwrap();
        let block = crate::sweep::extrude(model, floor, Vector::new(0.0, 0.0, h)).unwrap().shape;
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        let sides: Vec<DraftedFace> = model
            .unique_subshapes(block, NodeKind::Face)
            .into_iter()
            .filter(|&f| matches!(model.surface(f), Some(Surface::Plane(p)) if p.normal().z().abs() < 0.5))
            .map(|face| DraftedFace { face, angle: degrees.to_radians(), neutral: bottom })
            .collect();
        assert_eq!(sides.len(), n);
        draft(model, block, Direction::Z, &sides)
    }

    /// A box's sides drafted until its top is nearly gone: the long sides meet at 5√3, just
    /// above the top at 8; the section at height z is (a − 2zt)(b − 2zt), t = tan 30°.
    #[test]
    fn a_box_drafted_to_a_narrow_top() {
        let mut model = Model::new();
        let (a, b, h) = (30.0, 10.0, 8.0);
        let drafted = sides_drafted(&mut model, &[(0.0, 0.0), (a, 0.0), (a, b), (0.0, b)], h, 30.0).unwrap();
        valid(&model, drafted);
        let t = 30f64.to_radians().tan();
        let expected = a * b * h - (a + b) * t * h * h + 4.0 * t * t * h * h * h / 3.0;
        assert!((volume(&model, drafted).mass - expected).abs() < 1e-6 * expected, "{}", volume(&model, drafted).mass);
    }

    /// The same box 10 high: its long sides meet at 5√3, below the top, whose short edges the
    /// draft would turn round. Keeping the topology would make a solid inside out at the top:
    /// the box is cut by the tilted planes instead, closing along a ridge where the long sides
    /// meet (the volume, the integral of the section `(a − 2τz)(b − 2τz)` up to there).
    #[test]
    fn a_box_drafted_until_its_top_is_gone_closes_along_a_ridge() {
        let mut model = Model::new();
        let (a, b) = (30.0, 10.0);
        let drafted = sides_drafted(&mut model, &[(0.0, 0.0), (a, 0.0), (a, b), (0.0, b)], 10.0, 30.0).unwrap();
        valid(&model, drafted);
        let t = 30f64.to_radians().tan();
        let top = b / (2.0 * t);
        let expected = a * b * top - (a + b) * t * top * top + 4.0 * t * t * top.powi(3) / 3.0;
        assert!((volume(&model, drafted).mass - expected).abs() < 1e-9 * expected, "{}", volume(&model, drafted).mass);
    }

    /// An L-block 8 high, its arms 5 wide, its sides drafted 30°: they close over each arm at
    /// 2.5/tan 30° ≈ 4.33, below the top. OCCT fails it; kentron made a wrong solid.
    #[test]
    fn an_l_block_drafted_until_its_arms_close_is_refused() {
        let mut model = Model::new();
        let outline = [(0.0, 0.0), (20.0, 0.0), (20.0, 5.0), (5.0, 5.0), (5.0, 15.0), (0.0, 15.0)];
        assert!(sides_drafted(&mut model, &outline, 8.0, 30.0).is_err());
    }

    /// The face of `shape` whose surface passes through `p` (on the shape's boundary).
    fn face_through(model: &Model, shape: Shape, p: Point) -> Shape {
        model
            .unique_subshapes(shape, NodeKind::Face)
            .into_iter()
            .find(|&f| {
                let s = model.surface(f).unwrap().transformed(&f.placement);
                let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &s, None);
                s.point(u, v).distance(p) < 1e-9 && crate::boolean::classify_point(model, f, p).unwrap() == crate::boolean::PointLocation::Boundary
            })
            .expect("a face through the point")
    }

    fn locate(model: &Model, s: Shape, x: f64, y: f64, z: f64) -> crate::boolean::PointLocation {
        crate::boolean::classify_point(model, s, Point::new(x, y, z)).unwrap()
    }

    /// A box's top raised 3 and slid sideways along itself: the box 3 taller.
    #[test]
    fn a_box_top_moved_up() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let top = face_through(&model, b, Point::new(5.0, 10.0, 30.0));
        let moved = move_faces(&mut model, b, &[top], &Transform::translation(Vector::new(1.0, -2.0, 3.0))).unwrap();
        valid(&model, moved);
        let v = volume(&model, moved).mass;
        assert!((v - 6000.0 - 200.0 * 3.0).abs() < 1e-6, "{v}");
    }

    /// A box's side x = 10 turned by θ about its bottom edge (x = 10, z = 0): the side leans,
    /// the box a prism under a trapezoid of top 10 + h·tan θ.
    #[test]
    fn a_box_side_turned_about_its_bottom_edge_leans() {
        let (w, d, h) = (10.0, 20.0, 30.0);
        for theta in [0.2f64, -0.25] {
            let mut model = Model::new();
            let b = Cuboid::new(w, d, h).build(&mut model).unwrap();
            let side = face_at_x(&model, b, w);
            // (A turn by θ about +y carries +z toward +x: the top edge out to w + h·tan θ.)
            let edge = Axis::new(Point::new(w, 0.0, 0.0), Direction::Y);
            let moved = move_faces(&mut model, b, &[side], &Transform::rotation(edge, theta)).unwrap();
            valid(&model, moved);
            let expected = d * h * (w + 0.5 * h * theta.tan());
            let v = volume(&model, moved).mass;
            assert!((v - expected).abs() < 1e-6 * expected, "{theta}: {v} vs {expected}");
        }
    }

    /// A plate's through hole moved across it: the plate's volume unchanged, the hole at its
    /// new place.
    #[test]
    fn a_hole_moved_across_a_plate() {
        let mut model = Model::new();
        let (w, d, h, r) = (12.0, 8.0, 2.0, 1.5);
        let plate = Cuboid::new(w, d, h).build(&mut model).unwrap();
        let drill = CylinderPrimitive::new(r, 4.0).with_frame(crate::math::Frame::from_z(Point::new(4.0, 4.0, -1.0), Direction::Z)).build(&mut model).unwrap();
        let holed = crate::boolean::difference(&mut model, plate, drill).unwrap();
        let before = volume(&model, holed).mass;
        let wall = face_through(&model, holed, Point::new(4.0 + r, 4.0, 1.0));
        let moved = move_faces(&mut model, holed, &[wall], &Transform::translation(Vector::new(4.5, -1.5, 0.0))).unwrap();
        valid(&model, moved);
        let v = volume(&model, moved).mass;
        assert!((v - before).abs() < 1e-6 * before, "{v} vs {before}");
        use crate::boolean::PointLocation::*;
        assert_eq!(locate(&model, moved, 4.0, 4.0, 1.0), Inside);
        assert_eq!(locate(&model, moved, 8.5, 2.5, 1.0), Outside);
        assert_eq!(locate(&model, moved, 8.5 + r, 2.5, 1.0), Boundary);
        assert_eq!(locate(&model, moved, 8.5, 2.5 - r - 0.1, 1.0), Inside);
    }

    /// A notch through a block's side y = d (x from 2 to 4, 3 deep), its three faces moved 5
    /// along the side: the notch from 7 to 9.
    #[test]
    fn a_notch_moved_along_a_side() {
        let mut model = Model::new();
        let (w, d, h) = (12.0, 8.0, 3.0);
        let block = Cuboid::new(w, d, h).build(&mut model).unwrap();
        let cut = Cuboid::from_corners(Point::new(2.0, d - 3.0, -1.0), Point::new(4.0, d + 1.0, h + 1.0)).build(&mut model).unwrap();
        let notched = crate::boolean::difference(&mut model, block, cut).unwrap();
        let faces: Vec<Shape> = [Point::new(2.0, d - 1.5, 1.5), Point::new(4.0, d - 1.5, 1.5), Point::new(3.0, d - 3.0, 1.5)]
            .into_iter()
            .map(|p| face_through(&model, notched, p))
            .collect();
        let moved = move_faces(&mut model, notched, &faces, &Transform::translation(Vector::new(5.0, 0.0, 0.0))).unwrap();
        valid(&model, moved);
        let v = volume(&model, moved).mass;
        let expected = w * d * h - 2.0 * 3.0 * h;
        assert!((v - expected).abs() < 1e-6 * expected, "{v} vs {expected}");
        use crate::boolean::PointLocation::*;
        assert_eq!(locate(&model, moved, 3.0, d - 1.0, 1.5), Inside);
        assert_eq!(locate(&model, moved, 8.0, d - 1.0, 1.5), Outside);
        assert_eq!(locate(&model, moved, 8.0, d - 3.5, 1.5), Inside);
    }

    /// A box's top slid and turned along itself: the same plane, the same box.
    #[test]
    fn a_plane_moved_along_itself_is_unchanged() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let top = face_through(&model, b, Point::new(5.0, 10.0, 30.0));
        let turn = Transform::rotation(Axis::new(Point::new(1.0, 2.0, 0.0), Direction::Z), 0.3) * Transform::translation(Vector::new(2.0, 1.0, 0.0));
        let moved = move_faces(&mut model, b, &[top], &turn).unwrap();
        valid(&model, moved);
        let new_top = face_through(&model, moved, Point::new(5.0, 10.0, 30.0));
        assert_eq!(model.surface(new_top), model.surface(top));
        assert!((volume(&model, moved).mass - 6000.0).abs() < 1e-6);
    }

    #[test]
    fn face_moves_are_refused_or_fail() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let top = face_through(&model, b, Point::new(5.0, 10.0, 30.0));
        let shift = Transform::translation(Vector::new(0.0, 0.0, 1.0));
        assert!(matches!(move_faces(&mut model, b, &[], &shift), Err(Error::InvalidInput(_))));
        let scaled = Transform::scaling(Point::ORIGIN, 2.0);
        assert!(matches!(move_faces(&mut model, b, &[top], &scaled), Err(Error::InvalidInput(_))));
        let mirrored = Transform::mirror_plane(Point::new(0.0, 0.0, 15.0), Direction::Z);
        assert!(matches!(move_faces(&mut model, b, &[top], &mirrored), Err(Error::InvalidInput(_))));
        let sheared = Transform::from_rows([[1.0, 0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]).unwrap();
        assert!(matches!(move_faces(&mut model, b, &[top], &sheared), Err(Error::InvalidInput(_))));
        // The top stood upright, parallel to the sides x = 0 and x = 10: they cannot meet it.
        let upright = Transform::rotation(Axis::new(Point::new(5.0, 10.0, 30.0), Direction::Y), std::f64::consts::FRAC_PI_2);
        assert!(matches!(move_faces(&mut model, b, &[top], &upright), Err(Error::Failed(_))));
    }

    #[test]
    fn a_tiny_angle_changes_nothing() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let side = face_at_x(&model, b, 10.0);
        let bottom = Plane::new(Point::ORIGIN, Direction::Z);
        let drafted = draft(&mut model, b, Direction::Z, &[DraftedFace { face: side, angle: 1e-5, neutral: bottom }]).unwrap();
        assert!((volume(&model, drafted).mass - 6000.0).abs() < 1e-6);
    }
}
