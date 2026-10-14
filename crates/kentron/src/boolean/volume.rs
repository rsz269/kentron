//! Solids from the faces that enclose them.
//!
//! The faces of the arguments are split by one another (their union), then gathered into
//! closed shells: around each edge, the faces meeting there are ordered by the direction in
//! which they leave it, and the space on one side of a face continues across the edge to the
//! next face in that direction. Each piece of space so bounded is a solid, or the space
//! outside some faces altogether (a shell turned inwards): a cavity of the smallest solid
//! around it, or nothing if no solid holds it. Faces with a free edge enclose nothing and are
//! left out first.

use crate::error::{Error, Result};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Vector};
use crate::properties;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape};

use super::PointLocation;

/// The solids enclosed by the faces of `shapes` (solids, shells, faces), split by one
/// another: each closed piece of space they bound, sharing the faces where they meet. Faces
/// bounding no such piece are left out. One solid, or a compound of several (empty if the
/// faces enclose nothing).
#[doc(alias = "BOPAlgo_MakerVolume")]
#[doc(alias = "mkvolume")]
pub fn enclosed_solids(model: &mut Model, shapes: &[Shape]) -> Result<Shape> {
    enclose(model, shapes, true)
}

/// [`enclosed_solids`] of faces taken as they are, not split by one another: they must meet
/// along shared edges already.
pub fn enclosed_solids_unsplit(model: &mut Model, shapes: &[Shape]) -> Result<Shape> {
    enclose(model, shapes, false)
}

fn enclose(model: &mut Model, shapes: &[Shape], split: bool) -> Result<Shape> {
    // The arguments as faces (a solid's shells), split by one another.
    let mut sheets = Vec::new();
    for &shape in shapes {
        let faces = model.unique_subshapes(shape, NodeKind::Face);
        if faces.is_empty() {
            continue;
        }
        sheets.push(match shape.kind() {
            NodeKind::Face | NodeKind::Shell => shape,
            _ => compound(model, &faces)?,
        });
    }
    let Some(&first) = sheets.first() else {
        return Err(Error::InvalidInput("no faces to enclose solids with".into()));
    };
    let mut fused = first;
    if split {
        // (Keeping the arcs that bound nothing yet: later sheets may close regions with
        // them, as if all were fused at once.)
        for &sheet in &sheets[1..] {
            fused = super::mixed::with_slits(|| super::union(model, fused, sheet))?;
        }
    } else {
        fused = compound(model, &sheets)?;
    }
    let mut faces: Vec<Shape> = Vec::new();
    for f in model.unique_subshapes(fused, NodeKind::Face) {
        faces.push(without_slits(model, Shape { orientation: Orientation::Forward, ..f })?);
    }
    // Shells of half-faces (a face and the side of it the space lies on), dropping faces
    // that bound nothing until none is left.
    // (Faces of different sheets lying on one another, between the same edges, are one.)
    let faces = distinct_faces(model, faces);
    let mut faces = faces;
    let shells = loop {
        let sides = Sides::around_edges(model, &faces)?;
        let loose = sides.loose_faces();
        if loose.is_empty() {
            break sides.shells();
        }
        faces = faces.into_iter().enumerate().filter(|(i, _)| !loose.contains(i)).map(|(_, f)| f).collect();
    };
    // Each shell as a shell of the model, its faces facing out of its piece of space.
    let mut built = Vec::with_capacity(shells.len());
    for half_faces in &shells {
        let shell = model.new_empty(NodeKind::Shell);
        for &(i, side) in half_faces {
            // The space on the side the face's normal points to: the face turned away.
            model.add(shell, if side { faces[i].reversed() } else { faces[i] })?;
        }
        built.push((shell, properties::volume(model, shell).mass));
    }
    // (Shells enclosing next to nothing, between faces lying on one another, are left out.)
    let mut size = crate::math::BoundingBox::EMPTY;
    for &f in &faces {
        size.add_box(&properties::bounding_box(model, f));
    }
    let negligible = 1e-9 * (size.max - size.min).length().powi(3);
    let built: Vec<(Shape, f64)> = built.into_iter().map(|(s, v)| (s, if v.abs() <= negligible { 0.0 } else { v })).collect();
    let mut solids: Vec<(Shape, f64)> = Vec::new();
    for &(shell, volume) in built.iter().filter(|(_, v)| *v > 0.0) {
        let solid = model.new_empty(NodeKind::Solid);
        model.add(solid, shell)?;
        solids.push((solid, volume));
    }
    // A shell turned inwards bounds the space outside its faces: a cavity of the smallest
    // solid around it (not one it touches), if any.
    for (k, &(shell, volume)) in built.iter().enumerate() {
        if volume >= 0.0 {
            continue;
        }
        let Some(probe) = point_on(model, &faces[shells[k][0].0]) else { continue };
        let around = solids
            .iter()
            .filter(|(solid, _)| super::classify_point(model, *solid, probe).is_ok_and(|l| l == PointLocation::Inside))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|s| s.0);
        if let Some(solid) = around {
            model.add(solid, shell)?;
        }
    }
    // Tolerances grown to cover the gaps where pieces were joined.
    for &(solid, _) in &solids {
        crate::modify::update_tolerances(model, solid)?;
    }
    match solids.as_slice() {
        [(solid, _)] => Ok(*solid),
        _ => compound(model, &solids.iter().map(|s| s.0).collect::<Vec<_>>()),
    }
}

/// `faces` without those lying on another of them (on its surface, between the same
/// edges).
fn distinct_faces(model: &Model, faces: Vec<Shape>) -> Vec<Shape> {
    use crate::geometry::ParametricSurface;
    let mut kept: Vec<(Vec<crate::topology::NodeId>, Shape)> = Vec::new();
    for f in faces {
        let mut key: Vec<_> = model.unique_subshapes(f, NodeKind::Edge).iter().map(|e| e.node).collect();
        key.sort_unstable();
        let same = |r: Shape| -> bool {
            let (Some(s), Some(rs), Ok(id)) = (model.surface(f), model.surface(r), crate::topology::FaceId::try_from(f.node)) else { return false };
            let (u0, u1, v0, v1) = model.face_bounds(id);
            let p = s.point(0.5 * (u0 + u1), 0.5 * (v0 + v1));
            super::geom::signed_distance(rs, p).abs() <= super::geom::on_tol()
        };
        if !kept.iter().any(|(k, r)| *k == key && same(*r)) {
            kept.push((key, f));
        }
    }
    kept.into_iter().map(|(_, f)| f).collect()
}

/// `face` without its slits (loops running along edges and back), as a new face if it had
/// any.
fn without_slits(model: &mut Model, face: Shape) -> Result<Shape> {
    let wires: Vec<Shape> = model.children(face).collect();
    let slit = |model: &Model, wire: Shape| {
        let mut uses: std::collections::HashMap<crate::topology::NodeId, i32> = std::collections::HashMap::new();
        for c in model.children_local(wire) {
            let id = crate::topology::CoEdgeId::try_from(c.node).expect("a wire holds coedges");
            let o = if model[id].orientation == Orientation::Reversed { -1 } else { 1 };
            *uses.entry(model[id].edge.into()).or_default() += o;
        }
        !uses.is_empty() && uses.values().all(|&n| n == 0)
    };
    if !wires.iter().any(|&w| slit(model, w)) {
        return Ok(face);
    }
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.clone();
    let tolerance = model.tolerance(face);
    let new = model.new_empty(NodeKind::Face);
    model.set_face_surface(new, surface, tolerance)?;
    for w in wires {
        if !slit(model, w) {
            model.add(new, w)?;
        }
    }
    Ok(new)
}

/// How the faces meet around their edges: at each edge, the faces leaving it (a face twice
/// along a seam), in turn about the edge.
struct Sides {
    count: usize,
    /// Per edge: (face, direction about the edge's tangent in radians, whether turning
    /// positively about the tangent from the face goes to the side its normal points to).
    around: Vec<Vec<(usize, f64, bool)>>,
}

impl Sides {
    fn around_edges(model: &Model, faces: &[Shape]) -> Result<Sides> {
        let mut by_edge: Vec<(EdgeId, Vec<(usize, Point, Vector, Vector, Vector)>)> = Vec::new();
        for (i, &face) in faces.iter().enumerate() {
            let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
            let loops = crate::heal::face_loops(model, face).ok_or_else(|| Error::InvalidInput("a face without pcurves".into()))?;
            for u in loops.iter().flatten() {
                if model.is_degenerate(Shape::from(u.edge)) {
                    continue;
                }
                let Some((curve, first, last)) = model.edge_curve(Shape::from(u.edge)) else { continue };
                let (point, tangent) = curve.d1(0.5 * (first + last));
                // Into the face from the edge: to the left of its loop in the parameters.
                let uv = u.at(0.5);
                let ahead = u.at(0.5 + 1e-4) - u.at(0.5 - 1e-4);
                let (_, du, dv) = surface.d1(uv.x, uv.y);
                let inward = du * -ahead.y + dv * ahead.x;
                let normal = du.cross(dv);
                match by_edge.iter_mut().find(|(e, _)| *e == u.edge) {
                    Some((_, list)) => list.push((i, point, tangent, inward, normal)),
                    None => by_edge.push((u.edge, vec![(i, point, tangent, inward, normal)])),
                }
            }
        }
        let mut around = Vec::with_capacity(by_edge.len());
        for (_, list) in by_edge {
            let tangent = list[0].2.normalized();
            let flat = |v: Vector| v - tangent * v.dot(tangent);
            let x = flat(list[0].3).normalized();
            let y = tangent.cross(x);
            let mut sheets: Vec<(usize, f64, bool)> = list
                .iter()
                .map(|&(i, _, _, inward, normal)| {
                    let d = flat(inward);
                    (i, d.dot(y).atan2(d.dot(x)), normal.dot(tangent.cross(d)) > 0.0)
                })
                .collect();
            sheets.sort_by(|a, b| a.1.total_cmp(&b.1));
            around.push(sheets);
        }
        Ok(Sides { count: faces.len(), around })
    }

    /// The faces with an edge no other face meets (bounding no closed space).
    fn loose_faces(&self) -> Vec<usize> {
        let mut loose: Vec<usize> = self.around.iter().filter(|s| s.len() == 1).map(|s| s[0].0).collect();
        loose.sort_unstable();
        loose.dedup();
        loose
    }

    /// The closed shells: half-faces `(face, side)` (side true for the side the face's
    /// normal points to) joined across each edge to the next face about it on that side.
    fn shells(&self) -> Vec<Vec<(usize, bool)>> {
        let half = |face: usize, side: bool| 2 * face + side as usize;
        let mut parent: Vec<usize> = (0..2 * self.count).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        for sheets in &self.around {
            let n = sheets.len();
            for (k, &(face, _, positive_is_normal)) in sheets.iter().enumerate() {
                for side in [false, true] {
                    // Turning from the face towards that side, up to the next face.
                    let positive = side == positive_is_normal;
                    let next = if positive { sheets[(k + 1) % n] } else { sheets[(k + n - 1) % n] };
                    // The space lies on the next face's side facing back.
                    let next_side = positive != next.2;
                    let (a, b) = (root(&mut parent, half(face, side)), root(&mut parent, half(next.0, next_side)));
                    parent[a] = b;
                }
            }
        }
        let mut shells: Vec<(usize, Vec<(usize, bool)>)> = Vec::new();
        for face in 0..self.count {
            for side in [false, true] {
                let r = root(&mut parent, half(face, side));
                match shells.iter_mut().find(|(key, _)| *key == r) {
                    Some((_, list)) => list.push((face, side)),
                    None => shells.push((r, vec![(face, side)])),
                }
            }
        }
        shells.into_iter().map(|(_, list)| list).collect()
    }
}

/// A point inside a face (just off the middle of its first edge).
fn point_on(model: &Model, face: &Shape) -> Option<Point> {
    let surface = model.surface(*face)?;
    let loops = crate::heal::face_loops(model, *face)?;
    let u = loops.first()?.first()?;
    let ahead = u.at(0.5 + 1e-3) - u.at(0.5 - 1e-3);
    let uv = u.at(0.5);
    let inward = crate::math::Vector2::new(-ahead.y, ahead.x) * 0.5;
    Some(surface.point(uv.x + inward.x, uv.y + inward.y))
}

fn compound(model: &mut Model, parts: &[Shape]) -> Result<Shape> {
    let c = model.new_empty(NodeKind::Compound);
    for &p in parts {
        model.add(c, p)?;
    }
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::{Plane, Sphere, Surface};
    use crate::math::{Direction, Frame};

    #[test]
    fn two_concentric_spheres_enclose_a_ball_and_a_shell() {
        let mut model = Model::new();
        let sphere = |model: &mut Model, r: f64| {
            let s = Surface::Sphere(Sphere::new(Frame::WORLD, r));
            crate::make::face_on_domain(model, s, 0.0, std::f64::consts::TAU, -std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2).unwrap()
        };
        let (outer, inner) = (sphere(&mut model, 2.0), sphere(&mut model, 1.0));
        let result = enclosed_solids(&mut model, &[outer, inner]).unwrap();
        let solids = model.unique_subshapes(result, NodeKind::Solid);
        assert_eq!(solids.len(), 2);
        let mut volumes: Vec<f64> = solids.iter().map(|&s| properties::volume(&model, s).mass).collect();
        volumes.sort_by(f64::total_cmp);
        let ball = |r: f64| 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
        assert!((volumes[0] - ball(1.0)).abs() < 1e-3 * ball(1.0), "{volumes:?}");
        assert!((volumes[1] - (ball(2.0) - ball(1.0))).abs() < 1e-3 * ball(2.0), "{volumes:?}");
    }

    #[test]
    fn a_plane_across_a_box_encloses_two_halves() {
        let mut model = Model::new();
        let b = crate::primitives::Cuboid::new(2.0, 2.0, 2.0).build(&mut model).unwrap();
        let plane = Surface::Plane(Plane::new(Point::new(0.0, 0.0, 1.0), Direction::Z));
        let cutter = crate::make::face_on_domain(&mut model, plane, -1.0, 3.0, -1.0, 3.0).unwrap();
        let result = enclosed_solids(&mut model, &[b, cutter]).unwrap();
        let solids = model.unique_subshapes(result, NodeKind::Solid);
        assert_eq!(solids.len(), 2);
        for s in solids {
            assert!((properties::volume(&model, s).mass - 4.0).abs() < 1e-9);
        }
        // The halves share the face between them.
        assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), 11);
    }
}
