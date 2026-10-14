//! Shapes split by one another into cells.

use crate::error::Result;
use crate::topology::{Model, NodeKind, Shape};

/// A piece of a set of shapes split by one another: a solid (a face, an edge, for shapes
/// of lower dimension) lying inside some of them and outside the others.
#[doc(alias = "BOPAlgo_CellsBuilder")]
#[derive(Debug, Clone, PartialEq)]
pub struct Cell {
    /// Whether the cell lies inside each of the shapes, in their order.
    pub inside: Vec<bool>,
    pub shape: Shape,
}

/// `shapes` split by one another into cells: every piece of space (of their faces, edges,
/// for shapes of lower dimension) that lies inside one or more of them, with which. The
/// pieces of one shape together make it up; the cells of a union, an intersection or a
/// difference of them are those inside the right ones.
///
/// Solid cells share the faces, edges and vertices where they meet, as do the halves of
/// a solid split by a planar face (with the face's part between them); a face of a solid
/// split by later shapes is no longer shared with the cells that met it whole.
#[doc(alias = "BOPAlgo_Builder")]
#[doc(alias = "BOPAlgo_CellsBuilder")]
#[doc(alias = "general fuse")]
pub fn cells(model: &mut Model, shapes: &[Shape]) -> Result<Vec<Cell>> {
    let before = super::SPLITTING_CELLS.with(|f| f.replace(true));
    let found = split_all(model, shapes);
    super::SPLITTING_CELLS.with(|f| f.set(before));
    let found = found?;
    // (Split one after another, cells met by one shape have each their copies of its parts:
    // made one.)
    super::arguments::glue(model, &found.iter().map(|c| c.shape).collect::<Vec<_>>())?;
    Ok(found)
}

/// [`cells`], the shapes split one after another.
fn split_all(model: &mut Model, shapes: &[Shape]) -> Result<Vec<Cell>> {
    // Pieces inside the shapes added so far as marked, each split by every shape added
    // after it: by the next, into the part inside it and the part outside; and the part of
    // the next shape outside all of them. A shape of lower dimension than a piece (a face
    // across a solid) splits it without taking anything away, its own part inside the
    // piece a cell inside both.
    // (Each piece with its box once found: a piece meets many shapes.)
    let mut pieces: Vec<(Vec<bool>, Shape, Option<crate::math::BoundingBox>)> = Vec::new();
    for (k, &shape) in shapes.iter().enumerate() {
        let shape_box = crate::properties::bounding_box(model, shape);
        let shape_dimension = dimension(model, shape);
        let mut rest = Some(shape);
        let mut next = Vec::with_capacity(2 * pieces.len() + 1);
        for (mut inside, piece, piece_box) in pieces {
            let piece_box = piece_box.unwrap_or_else(|| crate::properties::bounding_box(model, piece));
            // (A part of the piece lies within its box; a part inside the shape too, within
            // both boxes.)
            let both = meet(&piece_box, &shape_box).unwrap_or(piece_box);
            if !piece_box.overlaps(&shape_box) {
                inside.push(false);
                next.push((inside, piece, Some(piece_box)));
                continue;
            }
            let piece_dimension = dimension(model, piece);
            if piece_dimension > shape_dimension {
                let (parts, cut) = split_by(model, piece, shape)?;
                if let Some(r) = rest {
                    // The shape's part inside the piece: where the piece's parts meet, when
                    // they were split by it (one face shared by both).
                    let within = match cut {
                        Some(cut) => cut,
                        None => super::intersection(model, r, piece)?,
                    };
                    if !super::is_empty(model, within) {
                        let mut marks = inside.clone();
                        marks.push(true);
                        next.push((marks, within, Some(both)));
                    }
                    rest = Some(super::difference(model, r, piece)?).filter(|&r| !super::is_empty(model, r));
                }
                inside.push(false);
                for part in parts {
                    next.push((inside.clone(), part, Some(piece_box)));
                }
                continue;
            }
            // Solids split together, sharing the faces where they meet.
            if piece_dimension == 3 && shape_dimension == 3 {
                if let Some(r) = rest {
                    let [within, without, beyond] = super::partition(model, piece, r)?;
                    rest = Some(beyond).filter(|&r| !super::is_empty(model, r));
                    if !super::is_empty(model, within) {
                        let mut marks = inside.clone();
                        marks.push(true);
                        next.push((marks, within, Some(both)));
                    }
                    if !super::is_empty(model, without) {
                        inside.push(false);
                        next.push((inside, without, Some(piece_box)));
                    }
                    continue;
                }
            }
            let within = super::intersection(model, piece, shape)?;
            let without = super::difference(model, piece, shape)?;
            if let Some(r) = rest {
                rest = Some(if shape_dimension > piece_dimension {
                    // (A lower-dimensional piece takes nothing away: it splits.)
                    let (parts, _) = split_by(model, r, piece)?;
                    compound(model, &parts)?
                } else {
                    super::difference(model, r, piece)?
                })
                .filter(|&r| !super::is_empty(model, r));
            }
            if !super::is_empty(model, within) {
                let mut marks = inside.clone();
                marks.push(true);
                next.push((marks, within, Some(both)));
            }
            if !super::is_empty(model, without) {
                inside.push(false);
                next.push((inside, without, Some(piece_box)));
            }
        }
        if let Some(r) = rest {
            let mut marks = vec![false; k];
            marks.push(true);
            next.push((marks, r, Some(shape_box)));
        }
        pieces = next;
    }
    // Each piece's solids (faces, edges) are cells of their own.
    let mut out = Vec::new();
    for (inside, piece, _) in pieces {
        for part in parts(model, piece) {
            out.push(Cell { inside: inside.clone(), shape: part });
        }
    }
    Ok(out)
}

/// The box where two boxes overlap; `None` if they do not.
fn meet(a: &crate::math::BoundingBox, b: &crate::math::BoundingBox) -> Option<crate::math::BoundingBox> {
    use crate::math::Point;
    let lo = Point::new(a.min.x.max(b.min.x), a.min.y.max(b.min.y), a.min.z.max(b.min.z));
    let hi = Point::new(a.max.x.min(b.max.x), a.max.y.min(b.max.y), a.max.z.min(b.max.z));
    let m = crate::math::BoundingBox::new(lo, hi);
    (!m.is_empty()).then_some(m)
}

/// 3 for a shape with solids, 2 with faces, 1 with edges, 0 otherwise.
fn dimension(model: &Model, shape: Shape) -> u8 {
    if !model.unique_subshapes(shape, NodeKind::Solid).is_empty() {
        3
    } else if !model.unique_subshapes(shape, NodeKind::Face).is_empty() {
        2
    } else if !model.unique_subshapes(shape, NodeKind::Edge).is_empty() {
        1
    } else {
        0
    }
}

/// The parts of `shape` either side of `cutter` of lower dimension: a planar face across a
/// solid splits it in two, sharing the face where they meet (also given); anything else
/// leaves it whole.
fn split_by(model: &mut Model, shape: Shape, cutter: Shape) -> Result<(Vec<Shape>, Option<Shape>)> {
    use crate::geometry::surface::Surface;
    use crate::math::Transformable;
    let faces = model.unique_subshapes(cutter, NodeKind::Face);
    let whole = Ok((parts(model, shape), None));
    let ([face], 3) = (faces.as_slice(), dimension(model, shape)) else { return whole };
    let Some(Surface::Plane(plane)) = model.surface(*face).map(|s| s.transformed(&face.placement)) else { return whole };
    // Across the solid: the face covers its shadow on the plane.
    let b = crate::properties::bounding_box(model, shape);
    let reach = (b.max - b.min).length();
    let corners = [b.min, b.max, crate::math::Point::new(b.min.x, b.min.y, b.max.z), crate::math::Point::new(b.max.x, b.max.y, b.min.z), crate::math::Point::new(b.min.x, b.max.y, b.min.z), crate::math::Point::new(b.max.x, b.min.y, b.max.z), crate::math::Point::new(b.min.x, b.max.y, b.max.z), crate::math::Point::new(b.max.x, b.min.y, b.min.z)];
    let (u0, u1, v0, v1) = model.face_uv_bounds(*face);
    let covered = corners.iter().all(|&c| {
        let (u, v) = plane.parameters(c);
        u >= u0 && u <= u1 && v >= v0 && v <= v1
    });
    let side = |c: crate::math::Point| (c - plane.origin()).dot(crate::math::Vector::from(plane.normal()));
    let crosses = corners.iter().any(|&c| side(c) > 0.0) && corners.iter().any(|&c| side(c) < 0.0);
    if !covered || !crosses {
        return whole;
    }
    let (u, v) = plane.parameters(b.min + (b.max - b.min) * 0.5);
    let bound = crate::make::face_on_domain(model, Surface::Plane(plane), u - reach, u + reach, v - reach, v + reach)?;
    let reference = plane.origin() + crate::math::Vector::from(plane.normal()) * reach.max(1e-9);
    let half = crate::primitives::half_space(model, bound, reference)?;
    let [within, without, _] = super::partition(model, shape, half)?;
    let mut out = parts(model, within);
    out.extend(parts(model, without));
    // The faces the two sides share: those on the plane.
    let on_plane: Vec<Shape> = model
        .unique_subshapes(within, NodeKind::Face)
        .into_iter()
        .filter(|&f| matches!(model.surface(f).map(|s| s.transformed(&f.placement)), Some(Surface::Plane(q)) if q.distance(plane.origin()) <= crate::tolerance::LINEAR && crate::math::Vector::from(q.normal()).cross(crate::math::Vector::from(plane.normal())).length() <= crate::tolerance::ANGULAR))
        .collect();
    let cut = if on_plane.is_empty() { None } else { Some(compound(model, &on_plane)?) };
    Ok((out, cut))
}

/// A compound of `parts`.
fn compound(model: &mut Model, parts: &[Shape]) -> Result<Shape> {
    let c = model.new_empty(NodeKind::Compound);
    for &p in parts {
        model.add(c, p)?;
    }
    Ok(c)
}

/// The separate parts of `shape`: its solids, or else its faces, or else its edges.
fn parts(model: &Model, shape: Shape) -> Vec<Shape> {
    for kind in [NodeKind::Solid, NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
        let found = model.unique_subshapes(shape, kind);
        if !found.is_empty() {
            return found;
        }
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point;
    use crate::primitives::Cuboid;

    #[test]
    fn two_overlapping_boxes_make_three_cells() {
        let mut model = Model::new();
        let a = Cuboid::new(2.0, 2.0, 2.0).build(&mut model).unwrap();
        let b = Cuboid::from_corners(Point::new(1.0, 0.0, 0.0), Point::new(3.0, 2.0, 2.0)).build(&mut model).unwrap();
        let found = cells(&mut model, &[a, b]).unwrap();
        // The cells meet on shared faces: three of six faces, two shared.
        let all = compound(&mut model, &found.iter().map(|c| c.shape).collect::<Vec<_>>()).unwrap();
        assert_eq!(model.unique_subshapes(all, NodeKind::Face).len(), 16);
        let volume = |inside: [bool; 2]| -> f64 {
            found.iter().filter(|c| c.inside == inside).map(|c| crate::properties::volume(&model, c.shape).mass).sum()
        };
        assert_eq!(found.len(), 3);
        assert!((volume([true, true]) - 4.0).abs() < 1e-9);
        assert!((volume([true, false]) - 4.0).abs() < 1e-9);
        assert!((volume([false, true]) - 4.0).abs() < 1e-9);
    }
}
