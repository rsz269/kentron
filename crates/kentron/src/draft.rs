//! Draft angles: faces tilted against a pull direction, as moulded parts need to leave
//! their mould.
//!
//! A drafted face keeps its place on a neutral plane and tilts by the draft angle relative
//! to the pull direction: a plane turns about its line on the neutral plane, a cylinder or
//! cone about the pull direction becomes the cone through its circle on the neutral plane.
//! Faces meeting a drafted one tangentially are drafted with it (those that can be). The shape keeps its
//! topology: its edges and vertices are found again where the new surfaces meet.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::surface::{Cone, Cylinder, Plane, Surface};
use crate::geometry::ParametricSurface;
use crate::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use crate::tolerance::ANGULAR;
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
    let drafted = crate::offset::shape::with_surfaces(model, shape, &new_surface)?;
    // (A loop the draft moves over a seam: the seam moved off it.)
    crate::heal::move_crossed_seams(model, drafted)?;
    Ok(drafted)
}

/// For each face of `shape`, the faces meeting it tangentially along an edge.
fn tangent_neighbours(model: &Model, shape: Shape) -> HashMap<NodeId, Vec<Shape>> {
    let mut faces_of_edge: HashMap<NodeId, Vec<Shape>> = HashMap::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        for edge in model.unique_subshapes(face, NodeKind::Edge) {
            let list = faces_of_edge.entry(edge.node).or_default();
            if !list.iter().any(|f| f.node == face.node) {
                list.push(face);
            }
        }
    }
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

/// The surface of a drafted face, as OCCT's Draft_Modification::NewSurface makes it.
/// `reversed`: the face uses its surface reversed. `None` if the face cannot be drafted.
fn drafted_surface(surface: &Surface, reversed: bool, direction: Direction, angle: f64, neutral: &Plane) -> Option<Surface> {
    let surface = match surface {
        Surface::Trimmed(t) => t.basis(),
        s => s,
    };
    let d = Vector::from(direction);
    // The frame's handedness and the face's orientation together give the outward side.
    let flip = |right_handed: bool| right_handed == reversed;
    match surface {
        Surface::Plane(plane) => {
            let (axis, theta) = plane_rotation(plane, flip(plane.frame().is_right_handed()), direction, angle, neutral)?;
            if theta.abs() <= ANGULAR {
                return Some(surface.clone());
            }
            Some(Surface::Plane(plane.transformed(&Transform::rotation(axis, theta))))
        }
        Surface::Cylinder(c) => {
            let z = Vector::from(c.frame().z());
            let n = Vector::from(neutral.normal());
            if d.dot(n).abs() <= 1.0 - ANGULAR || d.dot(z).abs() <= 1.0 - ANGULAR {
                return None;
            }
            if angle.abs() <= ANGULAR {
                return Some(surface.clone());
            }
            let mut alpha = if flip(c.frame().is_right_handed()) { -angle } else { angle };
            if d.dot(z) < 0.0 {
                alpha = -alpha;
            }
            // Where the axis crosses the neutral plane, as a height along it.
            let height = axis_height(c.frame().origin(), z, neutral)?;
            let mut radius = c.radius() + height * alpha.tan();
            if radius < 0.0 {
                radius = -radius;
            } else {
                alpha = -alpha;
            }
            Some(Surface::Cone(Cone::new(*c.frame(), alpha, radius)))
        }
        Surface::Cone(c) => {
            let z = Vector::from(c.frame().z());
            let n = Vector::from(neutral.normal());
            if d.dot(n).abs() <= 1.0 - ANGULAR || d.dot(z).abs() <= 1.0 - ANGULAR {
                return None;
            }
            let height = axis_height(c.frame().origin(), z, neutral)?;
            // The cone's circle on the neutral plane.
            let (sin, cos) = c.semi_angle_sin_cos();
            let circle_radius = c.radius() + height / cos * sin;
            if circle_radius <= 0.0 {
                return None;
            }
            if angle.abs() <= ANGULAR {
                return Some(Surface::Cylinder(Cylinder::new(*c.frame(), circle_radius)));
            }
            let mut alpha = if flip(c.frame().is_right_handed()) { -angle } else { angle };
            if d.dot(z) < 0.0 {
                alpha = -alpha;
            }
            let mut radius = circle_radius + height * alpha.tan();
            if radius < 0.0 {
                radius = -radius;
            } else {
                alpha = -alpha;
            }
            if (alpha - c.semi_angle()).abs() < ANGULAR {
                return Some(surface.clone());
            }
            Some(Surface::Cone(Cone::new(*c.frame(), alpha, radius)))
        }
        _ => None,
    }
}

/// The height along an axis (from `origin`, along `z`) where it crosses `plane`.
fn axis_height(origin: Point, z: Vector, plane: &Plane) -> Option<f64> {
    let n = Vector::from(plane.normal());
    let along = z.dot(n);
    if along.abs() <= ANGULAR {
        return None;
    }
    Some((plane.origin() - origin).dot(n) / along)
}

/// The axis (the plane's line on the neutral plane) and angle turning a plane to its draft
/// (OCCT's FindRotation). `flip`: the outward normal is opposite the frame's z.
fn plane_rotation(plane: &Plane, flip: bool, direction: Direction, angle: f64, neutral: &Plane) -> Option<(Axis, f64)> {
    let n = Vector::from(plane.frame().z());
    let m = Vector::from(neutral.normal());
    let line = Direction::try_from_vector(n.cross(m))?;
    let nx = Vector::from(line);
    // A point of both planes: the one nearest the plane's origin.
    let origin = {
        let p0 = plane.frame().origin();
        let gap = (neutral.origin() - p0).dot(m);
        // Move within the plane, across the line, until on the neutral plane.
        let across = n.cross(nx);
        let rate = across.dot(m);
        if rate.abs() <= 1e-300 {
            return None;
        }
        p0 + across * (gap / rate)
    };
    let ny = n.cross(nx);
    let d = Vector::from(direction);
    let a = d.dot(nx);
    if a.abs() > 1.0 - ANGULAR {
        return None;
    }
    let (mut b, mut c) = (d.dot(ny), d.dot(n));
    if flip {
        (b, c) = (-b, -c);
    }
    let denom = (1.0 - a * a).sqrt();
    let sin = angle.sin();
    if denom <= sin.abs() {
        return None;
    }
    let phi = (b / denom).atan2(c / denom);
    let theta0 = (sin / denom).acos();
    let mut theta = theta0 - phi;
    if theta.cos() < 0.0 {
        theta = -theta0 - phi;
    }
    while theta.abs() > std::f64::consts::PI {
        theta += if theta < 0.0 { std::f64::consts::PI } else { -std::f64::consts::PI };
    }
    Some((Axis::new(origin, line), theta))
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
