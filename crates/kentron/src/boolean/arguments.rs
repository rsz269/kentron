//! Boolean operations of several objects with several tools, the arguments split by one
//! another.

use crate::error::Result;
use crate::math::BoundingBox;
use crate::topology::{Model, NodeKind, Shape};

use super::{cells, geom, only_solids, PointClassifier, PointLocation};

/// What [`combine`] makes of the objects and the tools.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Combination {
    /// Everything inside any argument, as one.
    Union,
    /// What lies inside an object and a tool.
    Intersection,
    /// What lies inside an object and no tool.
    Difference,
    /// What lies inside a tool and no object.
    ReversedDifference,
    /// The edges where the arguments' boundaries meet.
    Section,
}

/// `combination` of several `objects` with several `tools`, every argument split by all the
/// others (those of its own side too: objects overlapping or touching one another split one
/// another, as OCCT's operations of several arguments do). An intersection or difference
/// is the pieces of space between the arguments' boundaries that it keeps, each a solid of
/// its own (sharing the faces where they meet); a section, the pieces of the boundaries'
/// edges lying on two arguments' boundaries, split where they meet. Where neither side's
/// arguments meet one another, or the arguments are not all solids, it is the operation of
/// each side's arguments together.
#[doc(alias = "BRepAlgoAPI_BooleanOperation")]
#[doc(alias = "BOPAlgo_BOP")]
pub fn combine(model: &mut Model, objects: &[Shape], tools: &[Shape], combination: Combination) -> Result<Shape> {
    let solids = |model: &Model, list: &[Shape]| -> Option<Vec<Shape>> {
        let mut out = Vec::new();
        for &s in list {
            out.extend(only_solids(model, s)?);
        }
        Some(out)
    };
    if let (Some(os), Some(ts)) = (solids(model, objects), solids(model, tools)) {
        if combination != Combination::Union && !os.is_empty() && !ts.is_empty() && (meeting(model, &os) || meeting(model, &ts)) {
            return split_together(model, &os, &ts, combination);
        }
    }
    // Each side's arguments together, then the operation of the two.
    let together = |model: &mut Model, list: &[Shape]| -> Result<Option<Shape>> {
        let Some((&first, rest)) = list.split_first() else { return Ok(None) };
        let mut all = first;
        for &s in rest {
            all = super::union(model, all, s)?;
        }
        Ok(Some(all))
    };
    if combination == Combination::Union {
        let all: Vec<Shape> = objects.iter().chain(tools).copied().collect();
        return Ok(together(model, &all)?.unwrap_or_else(|| model.new_empty(NodeKind::Compound)));
    }
    let (Some(o), Some(t)) = (together(model, objects)?, together(model, tools)?) else {
        return Ok(model.new_empty(NodeKind::Compound));
    };
    match combination {
        Combination::Intersection => super::intersection(model, o, t),
        Combination::Difference => super::difference(model, o, t),
        Combination::ReversedDifference => super::difference(model, t, o),
        Combination::Section => super::section(model, o, t),
        Combination::Union => unreachable!("handled above"),
    }
}

/// Whether two of `solids` meet (their boxes, touching included, overlap).
pub(super) fn meeting(model: &Model, solids: &[Shape]) -> bool {
    let boxes: Vec<BoundingBox> = solids.iter().map(|&s| crate::properties::bounding_box(model, s).enlarged(geom::on_tol())).collect();
    (0..boxes.len()).any(|i| (i + 1..boxes.len()).any(|j| boxes[i].overlaps(&boxes[j])))
}

/// [`combine`] of solids by their cells.
fn split_together(model: &mut Model, objects: &[Shape], tools: &[Shape], combination: Combination) -> Result<Shape> {
    let all: Vec<Shape> = objects.iter().chain(tools).copied().collect();
    let pieces = cells(model, &all)?;
    let result = model.new_empty(NodeKind::Compound);
    let n = objects.len();
    if combination == Combination::Section {
        // The cells' edges on the boundaries of two arguments (cells share the edges where
        // they meet).
        let mut edges: Vec<Shape> = Vec::new();
        let mut seen = crate::hashing::FastSet::default();
        for cell in &pieces {
            for e in model.unique_subshapes(cell.shape, NodeKind::Edge) {
                if seen.insert(e.node) {
                    edges.push(e);
                }
            }
        }
        let boxes: Vec<BoundingBox> = all.iter().map(|&s| crate::properties::bounding_box(model, s).enlarged(1e-6)).collect();
        let classifiers: Vec<PointClassifier> = all.iter().map(|&s| PointClassifier::new(model, s)).collect::<Result<_>>()?;
        for e in edges {
            let Some(p) = middle(model, e) else { continue };
            let on = (0..all.len()).filter(|&k| boxes[k].contains(p) && classifiers[k].locate(p) == PointLocation::Boundary).count();
            if on >= 2 {
                model.add(result, e)?;
            }
        }
        return Ok(result);
    }
    for cell in pieces {
        let (in_object, in_tool) = (cell.inside[..n].iter().any(|&x| x), cell.inside[n..].iter().any(|&x| x));
        let keep = match combination {
            Combination::Intersection => in_object && in_tool,
            Combination::Difference => in_object && !in_tool,
            Combination::ReversedDifference => in_tool && !in_object,
            Combination::Union | Combination::Section => unreachable!("handled by the callers"),
        };
        if keep {
            for solid in model.explore(cell.shape, NodeKind::Solid).collect::<Vec<_>>() {
                model.add(result, solid)?;
            }
        }
    }
    Ok(result)
}

/// The point half way along an edge's range.
fn middle(model: &Model, edge: Shape) -> Option<crate::math::Point> {
    use crate::geometry::ParametricCurve;
    let (c, first, last) = model.edge_curve(edge)?;
    Some(c.point(0.5 * (first + last)))
}

/// Makes `shapes` (cells, built one after another) share the vertices, edges and faces they
/// have in common: those at one point; edges between the same vertices lying on one another
/// and parametrized alike; faces bounded by the same edges on one surface (used the other
/// way round where their normals are opposite).
pub(super) fn glue(model: &mut Model, shapes: &[Shape]) -> Result<()> {
    use crate::geometry::{ParametricCurve, ParametricSurface};
    use crate::hashing::FastMap;
    use crate::topology::{FaceId, Orientation};
    let tol = geom::on_tol();
    let gather = |model: &Model, kind: NodeKind| -> Vec<Shape> {
        let mut seen = crate::hashing::FastSet::default();
        shapes.iter().flat_map(|&s| model.unique_subshapes(s, kind)).filter(|s| s.placement.is_identity() && seen.insert(s.node)).collect()
    };
    // Vertices at one point.
    let vertices = gather(model, NodeKind::Vertex);
    let cell = |p: crate::math::Point| [(p.x / (10.0 * tol)).round() as i64, (p.y / (10.0 * tol)).round() as i64, (p.z / (10.0 * tol)).round() as i64];
    let mut grid: FastMap<[i64; 3], Vec<Shape>> = FastMap::default();
    let mut merged: std::collections::HashMap<crate::topology::NodeId, crate::topology::NodeId> = Default::default();
    for v in vertices {
        let p = model.point(v);
        let c = cell(p);
        let near = (-1..=1).flat_map(|i| (-1..=1).flat_map(move |j| (-1..=1).map(move |k| [c[0] + i, c[1] + j, c[2] + k])));
        let found = near.filter_map(|k| grid.get(&k)).flatten().copied().find(|&w| model.point(w).distance(p) <= tol);
        match found {
            Some(w) => {
                merged.insert(v.node, w.node);
            }
            None => grid.entry(c).or_default().push(v),
        }
    }
    model.replace_nodes(&merged)?;
    merged.clear();
    // Edges between the same vertices, alike.
    let edges = gather(model, NodeKind::Edge);
    let mut by_ends: FastMap<(crate::topology::NodeId, crate::topology::NodeId), Vec<Shape>> = FastMap::default();
    for e in edges {
        let (a, b) = model.edge_vertices(e);
        let key = if a.node <= b.node { (a.node, b.node) } else { (b.node, a.node) };
        let Some((c, first, last)) = model.edge_curve(e) else { continue };
        let candidates = by_ends.entry(key).or_default();
        let alike = candidates.iter().copied().find(|&r| {
            let Some((rc, rf, rl)) = model.edge_curve(r) else { return false };
            (rf - first).abs() <= 1e-9 * (1.0 + first.abs())
                && (rl - last).abs() <= 1e-9 * (1.0 + last.abs())
                && [0.25, 0.5, 0.75].iter().all(|&s| {
                    let t = first + (last - first) * s;
                    c.point(t).distance(rc.point(t)) <= tol
                })
        });
        match alike {
            Some(r) => {
                merged.insert(e.node, r.node);
            }
            None => candidates.push(e),
        }
    }
    model.replace_nodes(&merged)?;
    // Faces bounded by the same edges, on one surface.
    let faces = gather(model, NodeKind::Face);
    // (The shells using each face, found once.)
    let mut shells_of: FastMap<crate::topology::NodeId, Vec<Shape>> = FastMap::default();
    for shell in gather(model, NodeKind::Shell) {
        for f in model.children_local(shell).map(|c| c.node).collect::<crate::hashing::FastSet<_>>() {
            shells_of.entry(f).or_default().push(shell);
        }
    }
    let mut by_edges: FastMap<Vec<crate::topology::NodeId>, Vec<Shape>> = FastMap::default();
    for f in faces {
        let mut key: Vec<_> = model.unique_subshapes(f, NodeKind::Edge).iter().map(|e| e.node).collect();
        key.sort_unstable();
        let Ok(id) = FaceId::try_from(f.node) else { continue };
        let Some(surface) = model.surface(f).cloned() else { continue };
        let (u0, u1, v0, v1) = model.face_bounds(id);
        let p = surface.point(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let candidates = by_edges.entry(key).or_default();
        let same = candidates.iter().copied().find(|&r| model.surface(r).is_some_and(|rs| geom::signed_distance(rs, p).abs() <= tol));
        let Some(r) = same else {
            candidates.push(f);
            continue;
        };
        let rs = model.surface(r).cloned().expect("checked");
        let opposite = geom::normal_at_point(&surface, p).dot(geom::normal_at_point(&rs, p)) < 0.0;
        // Each shell using the face uses the kept one instead, the same way round.
        for &shell in shells_of.get(&f.node).map(Vec::as_slice).unwrap_or_default() {
            let uses: Vec<Shape> = model.children_local(shell).filter(|c| c.node == f.node).collect();
            for old in uses {
                let flip = |o: Orientation| if opposite { o.reversed() } else { o };
                let new = Shape { node: r.node, orientation: flip(old.orientation), placement: old.placement };
                model.replace_child(shell, old, new)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point;
    use crate::primitives::Cuboid;

    /// Counts of distinct vertices, edges, faces and solids.
    fn counts(model: &Model, s: Shape) -> [usize; 4] {
        [NodeKind::Vertex, NodeKind::Edge, NodeKind::Face, NodeKind::Solid].map(|k| model.unique_subshapes(s, k).len())
    }

    #[test]
    fn overlapping_objects_split_one_another() {
        // Two overlapping boxes as the objects, a third as the tool (OCCT's counts).
        let mut model = Model::new();
        let block = |model: &mut Model, x: f64, y: f64, z: f64| Cuboid::from_corners(Point::new(x, y, z), Point::new(x + 10.0, y + 10.0, z + 10.0)).build(model).unwrap();
        let (o1, o2, t) = (block(&mut model, 0.0, 0.0, 0.0), block(&mut model, 5.0, 5.0, 0.0), block(&mut model, 2.0, 2.0, 5.0));
        let cases = [
            (Combination::Intersection, [20, 34, 18, 3], 440.0),
            (Combination::Difference, [32, 52, 24, 3], 1310.0),
            (Combination::ReversedDifference, [22, 34, 14, 1], 560.0),
            (Combination::Union, [34, 52, 20, 1], 2310.0),
        ];
        for (combination, expected, volume) in cases {
            let r = combine(&mut model, &[o1, o2], &[t], combination).unwrap();
            assert_eq!(counts(&model, r), expected, "{combination:?}");
            let v = crate::properties::volume(&model, r).mass;
            assert!((v - volume).abs() < 1e-9 * volume, "{combination:?} {v}");
            let report = crate::validity::check(&model, r);
            assert!(report.problems.is_empty(), "{combination:?} {report:?}");
        }
        // The section: where any two of the boxes meet.
        let r = combine(&mut model, &[o1, o2], &[t], Combination::Section).unwrap();
        assert_eq!(counts(&model, r)[..2], [22, 28]);
    }
}
