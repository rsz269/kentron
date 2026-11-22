//! The local properties of a shape's boundary at points: its outward normal and principal
//! curvatures, on the face (or faces) the boundary point lies on, as OCCT's
//! `BRepLProp_SLProps` gives them on a face (normals as `GeomLib::NormEstim` finds them for
//! a mesh's nodes, also at a sphere's pole).
//!
//! A boundary point on an edge or a vertex lies on several faces. Where they meet smoothly
//! (as a fillet meets the faces it rounds) the boundary still has a tangent plane there:
//! their outward normals agree within [`tolerance::ANGULAR`](crate::tolerance::ANGULAR), as
//! `BRepLib::ContinuityOfFaces` tells G1 faces (the sine between them within it, and not
//! opposite). Else the point lies on a crease and has neither normal nor curvatures.

use crate::distance::{self, NearFace};
use crate::error::Result;
use crate::geometry::props::{estimated_normal, surface_point};
use crate::geometry::surface::Surface;
use crate::math::{Direction, Point, Transformable, Vector};
use crate::tolerance::{ANGULAR, LINEAR};
use crate::topology::{Model, Orientation, Shape};

/// How far off a shape's boundary a point may lie and still have a normal from
/// [`normals`]: ten times [`tolerance::LINEAR`](crate::tolerance::LINEAR).
pub const ON_BOUNDARY: f64 = 10.0 * LINEAR;

/// The principal curvatures of a shape's boundary at a point, `max ≥ min`, in reciprocal
/// length: positive where the boundary bends away from its outward normal (as a ball's
/// does), negative where it bends towards it (as round a hole).
#[doc(alias = "BRepLProp_SLProps::MaxCurvature")]
#[doc(alias = "BRepLProp_SLProps::MinCurvature")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrincipalCurvatures {
    pub max: f64,
    pub min: f64,
}

/// The outward unit normal of `shape`'s boundary (its faces, each facing as it is used) at
/// each of `points`: `None` for a point off the boundary by more than [`ON_BOUNDARY`], on a
/// crease (an edge or vertex where faces meet at an angle), or where the surface has no
/// normal (a cone's apex).
#[doc(alias = "BRepLProp_SLProps::Normal")]
#[doc(alias = "GeomLib::NormEstim")]
pub fn normals(model: &mut Model, shape: Shape, points: &[Point]) -> Result<Vec<Option<Direction>>> {
    points
        .iter()
        .map(|&p| {
            let Some(at) = boundary_point(model, shape, p)? else { return Ok(None) };
            if at.distance > ON_BOUNDARY {
                return Ok(None);
            }
            Ok(at.normal())
        })
        .collect()
}

/// The principal curvatures of `shape`'s boundary at the boundary point nearest each of
/// `points` (which need not lie on the boundary): `None` where that point lies on a crease,
/// or where the surface's curvatures are undefined there (where its parameters are
/// singular, at a sphere's pole or a cone's apex). Where faces meet smoothly but bend
/// differently (a fillet and a face it rounds, at the edge between them), the curvatures of
/// the face nearest the given point (the first the shape lists, where equally near).
#[doc(alias = "BRepLProp_SLProps")]
#[doc(alias = "GeomLProp_SLProps")]
pub fn curvatures(model: &mut Model, shape: Shape, points: &[Point]) -> Result<Vec<Option<PrincipalCurvatures>>> {
    points
        .iter()
        .map(|&p| {
            let Some(at) = boundary_point(model, shape, p)? else { return Ok(None) };
            if at.normal().is_none() {
                return Ok(None);
            }
            let face = &at.faces[0];
            let Some(k) = surface_point(&face.surface, face.uv.0, face.uv.1, LINEAR).curvatures else { return Ok(None) };
            // (Signed along the surface's normal `D1U ^ D1V`, positive bending towards it:
            // along the outward normal, the other way round.)
            let sign = if face.reversed { 1.0 } else { -1.0 };
            let (a, b) = (sign * k.max, sign * k.min);
            Ok(Some(PrincipalCurvatures { max: a.max(b), min: a.min(b) }))
        })
        .collect()
}

/// The boundary point nearest a point, and the faces it lies on.
struct BoundaryPoint {
    /// How far the given point lies from it.
    distance: f64,
    /// The faces, the one nearest the given point first.
    faces: Vec<OnFace>,
}

/// A face a boundary point lies on, and the point's parameters on its surface.
struct OnFace {
    surface: Surface,
    uv: (f64, f64),
    reversed: bool,
}

impl OnFace {
    /// The face's outward normal at the point.
    fn normal(&self) -> Option<Direction> {
        let n = estimated_normal(&self.surface, self.uv.0, self.uv.1, LINEAR)?;
        Some(if self.reversed { -n } else { n })
    }
}

impl BoundaryPoint {
    /// The boundary's outward normal there: the faces', where they all agree.
    fn normal(&self) -> Option<Direction> {
        let first = self.faces[0].normal()?;
        for face in &self.faces[1..] {
            let n = face.normal()?;
            if !(Vector::from(first).cross(Vector::from(n)).length() <= ANGULAR && first.dot(n) > 0.0) {
                return None;
            }
        }
        Some(first)
    }
}

/// The point of `shape`'s boundary nearest `p` and the faces it lies on: those whose own
/// nearest points lie within [`LINEAR`] of it. `None` for a shape without faces.
fn boundary_point(model: &mut Model, shape: Shape, p: Point) -> Result<Option<BoundaryPoint>> {
    let near = distance::nearest_faces(model, shape, p, LINEAR)?;
    let Some(&NearFace { point: nearest, distance, .. }) = near.first() else { return Ok(None) };
    let mut faces = Vec::new();
    for f in near.iter().filter(|f| f.point.distance(nearest) <= LINEAR) {
        let local = Shape { orientation: Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..f.face };
        let Some(surface) = model.surface(local).map(|s| s.transformed(&f.face.placement)) else { continue };
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(f.point, &surface, None);
        faces.push(OnFace { surface, uv: (u, v), reversed: f.face.orientation == Orientation::Reversed });
    }
    if faces.is_empty() {
        return Ok(None);
    }
    Ok(Some(BoundaryPoint { distance, faces }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Frame;
    use crate::primitives::{Cone, Cuboid, Cylinder, Sphere, Torus};
    use crate::topology::NodeKind;

    fn normal(model: &mut Model, shape: Shape, p: [f64; 3]) -> Option<Direction> {
        normals(model, shape, &[Point::new(p[0], p[1], p[2])]).unwrap()[0]
    }

    fn bend(model: &mut Model, shape: Shape, p: [f64; 3]) -> Option<[f64; 2]> {
        curvatures(model, shape, &[Point::new(p[0], p[1], p[2])]).unwrap()[0].map(|k| [k.max, k.min])
    }

    fn assert_normal(got: Option<Direction>, want: [f64; 3]) {
        let got = got.expect("a normal");
        let want = Vector::new(want[0], want[1], want[2]).normalized();
        assert!((Vector::from(got) - want).length() < 1e-9, "{got:?} against {want:?}");
    }

    fn assert_bend(got: Option<[f64; 2]>, want: [f64; 2]) {
        let got = got.expect("curvatures");
        assert!((got[0] - want[0]).abs() < 1e-9 && (got[1] - want[1]).abs() < 1e-9, "{got:?} against {want:?}");
    }

    #[test]
    fn a_box_has_its_faces_normals_and_none_on_its_edges_and_corners() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        assert_normal(normal(&mut model, b, [5.0, 10.0, 30.0]), [0.0, 0.0, 1.0]);
        assert_normal(normal(&mut model, b, [0.0, 3.0, 4.0]), [-1.0, 0.0, 0.0]);
        // Right by an edge, on a face.
        assert_normal(normal(&mut model, b, [10.0, 1e-4, 15.0]), [1.0, 0.0, 0.0]);
        // On an edge, a crease, and on a corner.
        assert_eq!(normal(&mut model, b, [0.0, 0.0, 15.0]), None);
        assert_eq!(normal(&mut model, b, [10.0, 20.0, 30.0]), None);
        // Off the boundary, inside and out; within its tolerance, on it.
        assert_eq!(normal(&mut model, b, [5.0, 10.0, 15.0]), None);
        assert_eq!(normal(&mut model, b, [5.0, 10.0, 30.0 + 1e-5]), None);
        assert_normal(normal(&mut model, b, [5.0, 10.0, 30.0 + 5e-7]), [0.0, 0.0, 1.0]);
        // Flat everywhere off its creases, from points inside and out.
        assert_bend(bend(&mut model, b, [5.0, 10.0, 29.0]), [0.0, 0.0]);
        assert_bend(bend(&mut model, b, [5.0, -2.0, 15.0]), [0.0, 0.0]);
        // Nearest an edge or a corner: a crease.
        assert_eq!(bend(&mut model, b, [-1.0, -1.0, 15.0]), None);
        assert_eq!(bend(&mut model, b, [11.0, 21.0, 31.0]), None);
    }

    #[test]
    fn a_cylinder_bends_round_its_axis() {
        let mut model = Model::new();
        let c = Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
        let a = 0.7f64;
        assert_normal(normal(&mut model, c, [5.0 * a.cos(), 5.0 * a.sin(), 10.0]), [a.cos(), a.sin(), 0.0]);
        assert_bend(bend(&mut model, c, [5.0, 0.0, 10.0]), [0.2, 0.0]);
        assert_bend(bend(&mut model, c, [4.0 * a.cos(), 4.0 * a.sin(), 7.0]), [0.2, 0.0]);
        assert_bend(bend(&mut model, c, [7.0 * a.cos(), 7.0 * a.sin(), 13.0]), [0.2, 0.0]);
        // Its seam is smooth; its rims are creases.
        assert_normal(normal(&mut model, c, [5.0, 0.0, 3.0]), [1.0, 0.0, 0.0]);
        assert_eq!(normal(&mut model, c, [5.0 * a.cos(), 5.0 * a.sin(), 20.0]), None);
        assert_normal(normal(&mut model, c, [1.0, 2.0, 20.0]), [0.0, 0.0, 1.0]);
        assert_bend(bend(&mut model, c, [1.0, 2.0, 19.0]), [0.0, 0.0]);
    }

    #[test]
    fn a_ball_bends_alike_every_way() {
        let mut model = Model::new();
        let s = Sphere::new(4.0).build(&mut model).unwrap();
        let d = Vector::new(1.0, -2.0, 2.0).normalized();
        let on = |r: f64| [r * d.x, r * d.y, r * d.z];
        assert_normal(normal(&mut model, s, on(4.0)), [d.x, d.y, d.z]);
        assert_bend(bend(&mut model, s, on(4.0)), [0.25, 0.25]);
        assert_bend(bend(&mut model, s, on(1.5)), [0.25, 0.25]);
        assert_bend(bend(&mut model, s, on(9.0)), [0.25, 0.25]);
        // At its poles too, though its parameters are singular there (a normal only).
        assert_normal(normal(&mut model, s, [0.0, 0.0, 4.0]), [0.0, 0.0, 1.0]);
        assert_normal(normal(&mut model, s, [0.0, 0.0, -4.0]), [0.0, 0.0, -1.0]);
    }

    #[test]
    fn a_hole_bends_towards_its_normal() {
        let mut model = Model::new();
        let b = Cuboid::from_corners(Point::new(-10.0, -10.0, 0.0), Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
        let c = Cylinder::new(3.0, 20.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, -5.0), Direction::Z)).build(&mut model).unwrap();
        let holed = crate::boolean::difference(&mut model, b, c).unwrap();
        let a = 2.0f64;
        assert_normal(normal(&mut model, holed, [3.0 * a.cos(), 3.0 * a.sin(), 5.0]), [-a.cos(), -a.sin(), 0.0]);
        assert_bend(bend(&mut model, holed, [3.0 * a.cos(), 3.0 * a.sin(), 5.0]), [0.0, -1.0 / 3.0]);
        assert_bend(bend(&mut model, holed, [2.0 * a.cos(), 2.0 * a.sin(), 4.0]), [0.0, -1.0 / 3.0]);
        assert_bend(bend(&mut model, holed, [4.0 * a.cos(), 4.0 * a.sin(), 6.0]), [0.0, -1.0 / 3.0]);
        // Where the hole meets the top: a crease.
        assert_eq!(normal(&mut model, holed, [3.0 * a.cos(), 3.0 * a.sin(), 10.0]), None);
    }

    #[test]
    fn a_cone_bends_as_its_lean_and_radius_say_and_has_no_normal_at_its_apex() {
        let mut model = Model::new();
        let (bottom, height) = (4.0, 6.0);
        let c = Cone::new(bottom, 0.0, height).build(&mut model).unwrap();
        // The side leans α from the axis: tan α = 4 / 6.
        let alpha = (bottom / height).atan();
        let (a, z) = (1.1f64, 2.0);
        let rho = bottom * (1.0 - z / height);
        let on = [rho * a.cos(), rho * a.sin(), z];
        assert_normal(normal(&mut model, c, on), [a.cos() * alpha.cos(), a.sin() * alpha.cos(), alpha.sin()]);
        assert_bend(bend(&mut model, c, on), [alpha.cos() / rho, 0.0]);
        // At its apex, neither.
        assert_eq!(normal(&mut model, c, [0.0, 0.0, height]), None);
        assert_eq!(bend(&mut model, c, [0.0, 0.0, height + 1.0]), None);
    }

    #[test]
    fn a_ring_bends_round_its_tube_and_round_its_axis() {
        let mut model = Model::new();
        let (big, small) = (6.0, 2.0);
        let t = Torus::new(big, small).build(&mut model).unwrap();
        for phi in [0.3f64, 1.2, 2.5, -2.0] {
            let rho = big + small * phi.cos();
            let around = phi.cos() / rho;
            let want = [around.max(1.0 / small), around.min(1.0 / small)];
            let a = 0.4f64;
            let on = [rho * a.cos(), rho * a.sin(), small * phi.sin()];
            assert_normal(normal(&mut model, t, on), [phi.cos() * a.cos(), phi.cos() * a.sin(), phi.sin()]);
            assert_bend(bend(&mut model, t, on), want);
        }
    }

    #[test]
    fn a_filleted_edge_meets_its_faces_smoothly() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        // The edge along y at x = 10, z = 10.
        let edge = model
            .unique_subshapes(b, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (c, f, l) = model.edge_curve(e).unwrap();
                let (p, q) = (crate::geometry::ParametricCurve::point(&c, f), crate::geometry::ParametricCurve::point(&c, l));
                [p, q].iter().all(|x| (x.x - 10.0).abs() < 1e-9 && (x.z - 10.0).abs() < 1e-9)
            })
            .unwrap();
        let r = 2.0;
        let f = crate::blend::fillet(&mut model, b, &[edge], r).unwrap();
        // Where the fillet meets the top and the side: their normals.
        assert_normal(normal(&mut model, f, [10.0 - r, 5.0, 10.0]), [0.0, 0.0, 1.0]);
        assert_normal(normal(&mut model, f, [10.0, 5.0, 10.0 - r]), [1.0, 0.0, 0.0]);
        // On the fillet: round its axis, at x = z = 10 - r.
        let a = 0.6f64;
        let on = [10.0 - r + r * a.cos(), 5.0, 10.0 - r + r * a.sin()];
        assert_normal(normal(&mut model, f, on), [a.cos(), 0.0, a.sin()]);
        assert_bend(bend(&mut model, f, on), [1.0 / r, 0.0]);
        // The fillet's ends, where it meets the box's front and back, are creases.
        assert_eq!(normal(&mut model, f, [10.0 - r + r * a.cos(), 0.0, 10.0 - r + r * a.sin()]), None);
        // An edge left sharp.
        assert_eq!(normal(&mut model, f, [0.0, 5.0, 10.0]), None);
    }

    #[test]
    fn a_shape_without_faces_has_no_boundary() {
        let mut model = Model::new();
        let v = crate::make::vertex(&mut model, Point::new(1.0, 2.0, 3.0));
        assert_eq!(normal(&mut model, v, [1.0, 2.0, 3.0]), None);
        assert_eq!(bend(&mut model, v, [1.0, 2.0, 3.0]), None);
    }
}

