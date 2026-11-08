//! Offsets inwards of prisms, by their profile's offset: a solid whose faces are two parallel
//! planar caps and sides ruled along the caps' normal (planes, cylinders, extrusions), as an
//! extruded sketch is. Moved inwards, its side faces' offsets can cross themselves (a spline
//! profile turning tighter than the distance) or pass one another (a narrow part closing),
//! and offset surfaces cannot be untangled; its profile's offset in its plane can (offsets of
//! planar wires trim their own crossings and split). The inside of the offset is that
//! profile's offset swept between the caps moved in, the joins at its concave corners the
//! profile's (arcs round them, as the arc join's rolling balls along the edges, or the sides
//! extended): exactly.

use crate::error::{Error, Result};
use crate::geometry::surface::Surface;
use crate::math::{Transform, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Shape};

use super::Join;

/// A prism: its caps (the one the sides start from, and the other) and the vector between
/// their planes.
struct Prism {
    bottom: Shape,
    top: Shape,
    height: Vector,
}

/// `solid` as a prism, if it is one with a side that is not planar (prisms of planar faces
/// alone are offset face by face, their crossings trimmed in their planes).
fn prism(model: &Model, solid: Shape) -> Option<Prism> {
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    // The direction of the sides' rulings.
    let along = faces.iter().find_map(|&f| match model.surface(f)?.transformed(&f.placement) {
        Surface::Cylinder(c) => Some(Vector::from(c.frame().z())),
        Surface::Extrusion(e) => Some(Vector::from(e.direction())),
        _ => None,
    })?;
    let parallel = |v: Vector| v.cross(along).length() <= 1e-9 * v.length().max(1e-300);
    let mut caps = Vec::new();
    for &f in &faces {
        match model.surface(f)?.transformed(&f.placement) {
            Surface::Plane(p) if parallel(Vector::from(p.normal())) => caps.push((f, Vector::from(p.frame().origin() - crate::math::Point::ORIGIN).dot(along))),
            Surface::Plane(p) if Vector::from(p.normal()).dot(along).abs() <= 1e-9 => {}
            Surface::Cylinder(c) if parallel(Vector::from(c.frame().z())) => {}
            Surface::Extrusion(e) if parallel(Vector::from(e.direction())) => {}
            _ => return None,
        }
    }
    let [(a, za), (b, zb)] = caps[..] else { return None };
    let length = along.length();
    if !((zb - za).abs() > 1e-9 * length) {
        return None;
    }
    let (bottom, top, rise) = if za < zb { (a, b, zb - za) } else { (b, a, za - zb) };
    Some(Prism { bottom, top, height: along * (rise / (length * length)) })
}

/// The inside of `solid` (a prism) at `distance` (positive) from its faces, but past the caps
/// in `open` (which the inside then runs out through): `None` where `solid` is not a prism or
/// nothing is left.
fn inside(model: &mut Model, solid: Shape, distance: f64, open: &[NodeId], join: Join) -> Result<Option<Shape>> {
    let Some(p) = prism(model, solid) else { return Ok(None) };
    let h = p.height.length();
    let up = p.height * (1.0 / h);
    let lo = if open.contains(&p.bottom.node) { -distance } else { distance };
    let hi = h + if open.contains(&p.top.node) { distance } else { -distance };
    if !(hi - lo > crate::tolerance::LINEAR) {
        return Ok(None);
    }
    // The profile: the bottom cap's outer boundary offset inwards, its holes outwards.
    let wires = model.unique_subshapes(p.bottom, NodeKind::Wire);
    let face = crate::topology::FaceId::try_from(p.bottom.node).map_err(|_| Error::InvalidInput("a cap that is not a face".into()))?;
    let outer_id = model.outer_wire(face).ok_or_else(|| Error::InvalidInput("a cap without a boundary".into()))?;
    let outer = wires.iter().position(|w| w.node == NodeId::from(outer_id)).ok_or_else(|| Error::InvalidInput("a cap without a boundary".into()))?;
    let swept = |model: &mut Model, wire: Shape, by: f64, from: f64, to: f64| -> Result<Vec<Shape>> {
        let offset = super::planar_wire(model, wire, by, join)?;
        let loops = if offset.kind() == NodeKind::Wire { vec![offset] } else { model.unique_subshapes(offset, NodeKind::Wire) };
        let mut out = Vec::new();
        for l in loops {
            let face = crate::make::planar_face(model, l)?;
            let placed = crate::modify::transform(model, face, &Transform::translation(up * from), Default::default())?;
            out.push(crate::sweep::extrude(model, placed, up * (to - from))?.shape);
        }
        Ok(out)
    };
    let mut pieces = swept(model, wires[outer], -distance, lo, hi)?;
    if pieces.is_empty() {
        return Ok(None);
    }
    let mut inner = pieces.remove(0);
    for piece in pieces {
        inner = crate::boolean::union(model, inner, piece)?;
    }
    for (i, &w) in wires.iter().enumerate() {
        if i == outer {
            continue;
        }
        for hole in swept(model, w, distance, lo - distance, hi + distance)? {
            inner = crate::boolean::difference(model, inner, hole)?;
        }
    }
    Ok(Some(inner))
}

/// `solid` (a prism) offset inwards by `distance` (negative): the inside at that distance.
pub(super) fn offset(model: &mut Model, solid: Shape, distance: f64, join: Join) -> Result<Option<Shape>> {
    if !(distance < 0.0) {
        return Ok(None);
    }
    let (copy, _) = crate::modify::copy_unplaced(model, solid, &[])?;
    inside(model, copy, -distance, &[], join)
}

/// `solid` (a prism) hollowed to walls `thickness` thick inwards (negative), open at the
/// caps in `removed`: the solid less its inside at that distance, run out through them.
pub(super) fn thick(model: &mut Model, solid: Shape, removed: &[NodeId], thickness: f64, join: Join) -> Result<Option<Shape>> {
    if !(thickness < 0.0) {
        return Ok(None);
    }
    let Some(p) = prism(model, solid) else { return Ok(None) };
    if removed.iter().any(|r| *r != p.bottom.node && *r != p.top.node) {
        return Ok(None);
    }
    // (On a copy: an offset tried face by face before may have taken the solid's faces into
    // its walls. Its caps, the copy's at the same ends.)
    let (copy, _) = crate::modify::copy_unplaced(model, solid, &[])?;
    let Some(q) = prism(model, copy) else { return Ok(None) };
    let open: Vec<NodeId> = [(p.bottom, q.bottom), (p.top, q.top)].into_iter().filter(|(c, _)| removed.contains(&c.node)).map(|(_, d)| d.node).collect();
    let Some(cavity) = inside(model, copy, -thickness, &open, join)? else { return Ok(None) };
    crate::boolean::difference(model, copy, cavity).map(Some)
}
