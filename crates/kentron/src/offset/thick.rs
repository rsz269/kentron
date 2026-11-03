//! Shells thickened into solids: the solid between a face or shell and its offset.

use crate::error::{Error, Result};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::Join;

/// The solid between `shape` (a face, a shell or faces) and its offset by `distance` along
/// its faces' normals (against them for a negative distance), joined as `join` says: the
/// faces, their offsets and, along each free edge, a ruled wall between the edge and its
/// offset, sewn. A closed shell thickens into a solid of two shells, the one inside the
/// other a cavity.
#[doc(alias = "thickshell")]
#[doc(alias = "BRepOffset_MakeOffset::Initialize(thickening)")]
pub fn thick_shell(model: &mut Model, shape: Shape, distance: f64, join: Join) -> Result<Shape> {
    let faces = model.unique_subshapes(shape, NodeKind::Face);
    if faces.is_empty() {
        return Err(Error::InvalidInput("a thick shell of no faces".into()));
    }
    if distance == 0.0 {
        return Err(Error::InvalidInput("a thick shell of no thickness".into()));
    }
    let shell = if shape.kind() == NodeKind::Shell {
        shape
    } else {
        let shell = model.new_empty(NodeKind::Shell);
        for &f in &faces {
            model.add(shell, f)?;
        }
        shell
    };
    let moved = super::shape::shape(model, shell, distance, join)?;
    // The free edges of each (used by one face only, and not degenerate), with their faces.
    let free = |model: &Model, s: Shape| -> Vec<(Shape, Shape)> {
        let mut uses: crate::hashing::FastMap<NodeId, (usize, Shape, Shape)> = Default::default();
        for f in model.unique_subshapes(s, NodeKind::Face) {
            for e in model.explore(f, NodeKind::Edge) {
                let entry = uses.entry(e.node).or_insert((0, e, f));
                entry.0 += 1;
            }
        }
        let mut out: Vec<(Shape, Shape)> = uses.into_values().filter(|&(n, e, _)| n == 1 && !model.is_degenerate(e)).map(|(_, e, f)| (e, f)).collect();
        out.sort_by_key(|(e, _)| e.node);
        out
    };
    let middle = |model: &Model, e: Shape| -> Option<Point> {
        let (c, a, b) = model.edge_curve(Shape { orientation: Orientation::Forward, ..e })?;
        Some(c.point(0.5 * (a + b)))
    };
    let normal = |model: &Model, f: Shape, p: Point| -> Option<Vector> {
        let surface = model.surface(f)?.transformed(&f.placement);
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
        let (_, du, dv) = surface.d1(u, v);
        let n = du.cross(dv);
        (n.length() > 1e-300).then(|| n.normalized() * if f.orientation == Orientation::Reversed { -1.0 } else { 1.0 })
    };
    let mine = free(model, shell);
    let theirs: Vec<(Shape, Point)> = free(model, moved).into_iter().filter_map(|(e, _)| Some((e, middle(model, e)?))).collect();
    let mut pieces: Vec<Shape> = faces.clone();
    pieces.extend(model.unique_subshapes(moved, NodeKind::Face));
    // A ruled wall between each free edge and its offset (the offset edge nearest where the
    // edge's middle moves, run the same way).
    for &(e, f) in &mine {
        let (Some(m), Some(n)) = (middle(model, e), middle(model, e).and_then(|m| normal(model, f, m))) else { continue };
        let target = m + n * distance;
        let Some(&(g, _)) = theirs.iter().min_by(|a, b| a.1.distance(target).total_cmp(&b.1.distance(target))) else { continue };
        let ends = |model: &Model, x: Shape| {
            let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..x });
            (model.point(a), model.point(b))
        };
        let ((a0, a1), (b0, b1)) = (ends(model, e), ends(model, g));
        let g = if a0.distance(b0) + a1.distance(b1) > a0.distance(b1) + a1.distance(b0) { Shape { orientation: Orientation::Reversed, ..g } } else { Shape { orientation: Orientation::Forward, ..g } };
        let we = crate::make::wire(model, &[Shape { orientation: Orientation::Forward, ..e }])?;
        let wg = crate::make::wire(model, &[g])?;
        let wall = crate::loft::loft(model, &[we, wg], &crate::loft::LoftOptions { solid: false, ruled: true, ..Default::default() })?.shape;
        pieces.extend(model.unique_subshapes(wall, NodeKind::Face));
    }
    // (A closed shell, no free edges: it and its offset bound the solid as they are.)
    let shells = if mine.is_empty() {
        let mut shells = vec![shell];
        shells.extend(model.unique_subshapes(moved, NodeKind::Shell));
        shells
    } else {
        let size = crate::properties::bounding_box(model, moved).size().length().max(crate::properties::bounding_box(model, shell).size().length());
        let sewn = crate::sew::sew(model, &pieces, &crate::sew::SewOptions { tolerance: 1e-6 * size.max(1.0), ..Default::default() })?;
        model.unique_subshapes(sewn, NodeKind::Shell)
    };
    if shells.is_empty() {
        return Err(Error::Failed("a thick shell whose pieces do not sew".into()));
    }
    // The outermost shell (enclosing the most) facing out, the others (cavities) in.
    let volume_of = |model: &mut Model, s: Shape| -> Result<f64> {
        let solid = model.new_empty(NodeKind::Solid);
        model.add(solid, s)?;
        Ok(crate::properties::volume(model, solid).mass)
    };
    let mut volumes = Vec::with_capacity(shells.len());
    for &s in &shells {
        volumes.push(volume_of(model, s)?);
    }
    let outer = (0..shells.len()).max_by(|&a, &b| volumes[a].abs().total_cmp(&volumes[b].abs())).expect("a shell");
    let solid = model.new_empty(NodeKind::Solid);
    for (i, &s) in shells.iter().enumerate() {
        let outwards = i == outer;
        model.add(solid, if (volumes[i] > 0.0) == outwards { s } else { s.reversed() })?;
    }
    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Sphere};

    /// A sphere's face (closed) thickened 10 outwards: a solid of two shells, its volume the
    /// spheres' difference (OCCT's bugs offset bug26233).
    #[test]
    fn a_sphere_thickened_into_a_hollow_ball() {
        let mut model = Model::new();
        let ball = Sphere::new(100.0).build(&mut model).unwrap();
        let face = model.unique_subshapes(ball, NodeKind::Face)[0];
        let thick = thick_shell(&mut model, face, 10.0, Join::Intersection).unwrap();
        assert!(crate::validity::check(&model, thick).is_valid());
        assert_eq!(model.unique_subshapes(thick, NodeKind::Shell).len(), 2);
        let v = crate::properties::volume(&model, thick).mass;
        let expected = 4.0 / 3.0 * std::f64::consts::PI * (110f64.powi(3) - 100f64.powi(3));
        assert!((v - expected).abs() < 1e-6 * expected, "{v}");
    }

    /// A box's top face thickened 2 downwards: the slab under it, walled at its edges.
    #[test]
    fn a_square_thickened_into_a_slab() {
        let mut model = Model::new();
        let block = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 6.0, 4.0)).build(&mut model).unwrap();
        let top = model
            .unique_subshapes(block, NodeKind::Face)
            .into_iter()
            .find(|&f| crate::properties::area(&model, f).center_of_mass.z > 3.9)
            .unwrap();
        let thick = thick_shell(&mut model, top, -2.0, Join::Intersection).unwrap();
        assert!(crate::validity::check(&model, thick).is_valid());
        assert_eq!(model.unique_subshapes(thick, NodeKind::Face).len(), 6);
        let v = crate::properties::volume(&model, thick).mass;
        assert!((v - 120.0).abs() < 1e-9 * 120.0, "{v}");
    }
}
