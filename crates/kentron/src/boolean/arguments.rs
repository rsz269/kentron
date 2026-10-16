//! Boolean operations of several objects with several tools, the arguments split by one
//! another.

use crate::error::Result;
use crate::math::BoundingBox;
use crate::topology::{Model, NodeKind, Shape};

use super::{geom, only_solids, PointClassifier, PointLocation};

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
    // Each argument's solids apart (one may hold several), else the argument whole.
    // (A compound's members are arguments of their own: all are split by one another.)
    fn members(model: &Model, s: Shape, out: &mut Vec<Shape>) {
        if s.kind() == NodeKind::Compound {
            for c in model.children(s).collect::<Vec<_>>() {
                members(model, c, out);
            }
        } else {
            out.push(s);
        }
    }
    let flatten = |model: &Model, list: &[Shape]| -> Vec<Shape> {
        let mut out = Vec::new();
        for &s in list {
            members(model, s, &mut out);
        }
        out.into_iter().flat_map(|s| only_solids(model, s).unwrap_or_else(|| vec![s])).filter(|&s| !super::is_empty(model, s)).collect()
    };
    let (os, ts) = (flatten(model, objects), flatten(model, tools));
    let dims = |model: &Model, list: &[Shape]| -> Vec<u8> { list.iter().map(|&s| super::cells::dimension(model, s)).collect() };
    let (od, td) = (dims(model, &os), dims(model, &ts));
    // (As OCCT's: a union of arguments of one dimension; a difference of objects of no
    // higher dimension than the tools.)
    let lowest = |d: &[u8]| d.iter().copied().min().unwrap_or(0);
    let highest = |d: &[u8]| d.iter().copied().max().unwrap_or(0);
    let refuse = |what: &str| Err(crate::error::Error::InvalidInput(format!("the {what} of arguments of these dimensions is not defined")));
    match combination {
        Combination::Union if od.iter().chain(&td).any(|&d| d != od.first().or(td.first()).copied().unwrap_or(0)) => return refuse("union"),
        Combination::Difference if !ts.is_empty() && highest(&od) > lowest(&td) => return refuse("difference"),
        Combination::ReversedDifference if !os.is_empty() && highest(&td) > lowest(&od) => return refuse("difference"),
        _ => {}
    }
    let mixed = od.iter().chain(&td).any(|&d| Some(d) != od.first().copied());
    let all_solids = od.iter().chain(&td).all(|&d| d == 3);
    if combination != Combination::Union && !os.is_empty() && !ts.is_empty() && (mixed || (all_solids && (meeting(model, &os) || meeting(model, &ts)))) {
        // An intersection whose pieces (each object's with each tool) lie apart: they are its
        // cells (another argument's boundary splits none: it would have to pass through the
        // inside of both).
        if combination == Combination::Intersection && all_solids {
            if let Some(apart) = pieces_apart(model, &os, &ts)? {
                return Ok(apart);
            }
        }
        return split_together(model, &os, &ts, combination);
    }
    // A section of objects alone: where they meet one another.
    if combination == Combination::Section && ts.is_empty() && os.len() > 1 {
        return split_together(model, &os, &ts, combination);
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

/// The intersections of each of `objects` with each of `tools` (solids), gathered, when
/// they lie apart from one another; `None` if two of them meet.
fn pieces_apart(model: &mut Model, objects: &[Shape], tools: &[Shape]) -> Result<Option<Shape>> {
    let boxes = |model: &Model, list: &[Shape]| -> Vec<BoundingBox> { list.iter().map(|&s| crate::properties::bounding_box(model, s).enlarged(geom::on_tol())).collect() };
    let (ob, tb) = (boxes(model, objects), boxes(model, tools));
    let mut pieces: Vec<(Shape, BoundingBox)> = Vec::new();
    for (i, &o) in objects.iter().enumerate() {
        for (j, &t) in tools.iter().enumerate() {
            if !ob[i].overlaps(&tb[j]) {
                continue;
            }
            let common = super::run(model, o, t, super::Operation::Intersection)?;
            for solid in model.explore(common, NodeKind::Solid).collect::<Vec<_>>() {
                let b = crate::properties::bounding_box(model, solid).enlarged(geom::on_tol());
                if pieces.iter().any(|(_, other)| other.overlaps(&b)) {
                    return Ok(None);
                }
                pieces.push((solid, b));
            }
        }
    }
    let result = model.new_empty(NodeKind::Compound);
    for (solid, _) in pieces {
        model.add(result, solid)?;
    }
    Ok(Some(result))
}

/// Whether two of `solids` meet (their boxes, touching included, overlap).
pub(super) fn meeting(model: &Model, solids: &[Shape]) -> bool {
    let boxes: Vec<BoundingBox> = solids.iter().map(|&s| crate::properties::bounding_box(model, s).enlarged(geom::on_tol())).collect();
    (0..boxes.len()).any(|i| (i + 1..boxes.len()).any(|j| boxes[i].overlaps(&boxes[j])))
}

/// [`combine`] of solids by their cells.
fn split_together(model: &mut Model, objects: &[Shape], tools: &[Shape], combination: Combination) -> Result<Shape> {
    let all: Vec<Shape> = objects.iter().chain(tools).copied().collect();
    // (A section needs no solid split by faces: only where the arguments meet.)
    let pieces = super::cells::cells_with(model, &all, combination != Combination::Section)?;
    let result = model.new_empty(NodeKind::Compound);
    let n = objects.len();
    if combination == Combination::Section {
        return section_of_cells(model, &all, &pieces);
    }
    // Kept faces and edges gathered in a shell or wire where their argument is one.
    let containers: Vec<Option<Shape>> = all
        .iter()
        .map(|&s| {
            if model.explore(s, NodeKind::Solid).next().is_none() && model.explore(s, NodeKind::Shell).next().is_some() {
                Some(model.new_empty(NodeKind::Shell))
            } else if model.explore(s, NodeKind::Face).next().is_none() && model.explore(s, NodeKind::Wire).next().is_some() {
                Some(model.new_empty(NodeKind::Wire))
            } else {
                None
            }
        })
        .collect();
    let dims: Vec<u8> = all.iter().map(|&s| super::cells::dimension(model, s)).collect();
    let mut used = vec![false; all.len()];
    // (Container first or not, the argument, the dimension, the part.)
    let mut items: Vec<(u8, usize, u8, Shape)> = Vec::new();
    for cell in pieces {
        let (in_object, in_tool) = (cell.inside[..n].iter().any(|&x| x), cell.inside[n..].iter().any(|&x| x));
        let keep = match combination {
            Combination::Intersection => in_object && in_tool,
            Combination::Difference => in_object && !in_tool,
            Combination::ReversedDifference => in_tool && !in_object,
            Combination::Union | Combination::Section => unreachable!("handled by the callers"),
        };
        if !keep {
            continue;
        }
        let dim = super::cells::dimension(model, cell.shape);
        // The argument it is a part of (of its dimension).
        let owner = (0..all.len()).find(|&k| cell.inside[k] && dims[k] == dim);
        match owner.and_then(|k| containers[k].map(|c| (k, c))) {
            Some((k, container)) if dim == 2 => {
                for f in model.unique_subshapes(cell.shape, NodeKind::Face) {
                    model.add(container, f)?;
                }
                used[k] = true;
            }
            Some((k, container)) if dim == 1 => {
                for e in model.unique_subshapes(cell.shape, NodeKind::Edge) {
                    let id = crate::topology::EdgeId::try_from(e.node).expect("an edge");
                    let c = model.add_coedge(id, e.orientation)?;
                    let wire = crate::topology::WireId::try_from(container.node).expect("a wire");
                    model.add(wire.into(), c.into())?;
                }
                used[k] = true;
            }
            _ => {
                for part in super::cells::parts_of(model, cell.shape) {
                    items.push((1, owner.unwrap_or(all.len()), dim, part));
                }
            }
        }
    }
    for (k, c) in containers.into_iter().enumerate() {
        if let (Some(c), true) = (c, used[k]) {
            items.push((0, k, dims[k], c));
        }
    }
    // An edge (or point) lying on a kept face or solid is part of it (a face on a solid's
    // boundary is kept apart, as OCCT does).
    let mut within: crate::hashing::FastSet<crate::topology::NodeId> = Default::default();
    for &(_, _, d, item) in &items {
        for kind in [NodeKind::Edge, NodeKind::Vertex] {
            if super::cells::dimension_of_kind(kind) < d {
                within.extend(model.unique_subshapes(item, kind).iter().map(|s| s.node));
            }
        }
    }
    // As OCCT's results list them: the shells and wires first, then the rest, each in the
    // order of the arguments they are parts of.
    items.sort_by_key(|&(container, k, _, _)| (container, k));
    for (container, _, _, item) in items {
        if container == 1 && within.contains(&item.node) {
            continue;
        }
        model.add(result, item)?;
    }
    Ok(result)
}

/// The section of `all` from their cells: the cells' edges lying on two arguments, each
/// stretch once (those of cells of lower dimension first: a solid left whole keeps edges
/// its crossing faces' pieces split), split at the section's vertices lying on them; and
/// the points where the arguments touch elsewhere.
fn section_of_cells(model: &mut Model, all: &[Shape], pieces: &[super::Cell]) -> Result<Shape> {
    use crate::geometry::ParametricCurve;
    use crate::hashing::{FastMap, FastSet};
    use crate::topology::NodeId;
    let result = model.new_empty(NodeKind::Compound);
    let mut edges: Vec<(u8, Shape)> = Vec::new();
    let mut seen = FastSet::default();
    for cell in pieces {
        let dim = super::cells::dimension(model, cell.shape);
        for e in model.unique_subshapes(cell.shape, NodeKind::Edge) {
            if seen.insert(e.node) {
                edges.push((dim, e));
            }
        }
    }
    edges.sort_by_key(|&(d, _)| d);
    let boxes: Vec<BoundingBox> = all.iter().map(|&s| crate::properties::bounding_box(model, s).enlarged(1e-6)).collect();
    let dims: Vec<u8> = all.iter().map(|&s| super::cells::dimension(model, s)).collect();
    let classifiers: Vec<Option<PointClassifier>> = all.iter().zip(&dims).map(|(&s, &d)| if d == 3 { PointClassifier::new(model, s).map(Some) } else { Ok(None) }).collect::<Result<_>>()?;
    // (Within an argument's own tolerance, for tolerant shapes read from files.)
    let tolerances: Vec<f64> = all.iter().map(|&s| max_tolerance(model, s).max(1e-6)).collect();
    let on_two = |model: &mut Model, p: crate::math::Point| -> Result<bool> {
        let mut on = 0;
        for k in 0..all.len() {
            if boxes[k].contains(p) && on_argument(model, &classifiers[k], all[k], p, tolerances[k])? {
                on += 1;
                if on >= 2 {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    };
    // A grid of cells over the arguments, to find kept edges and vertices near a point.
    let mut scene = BoundingBox::EMPTY;
    for b in &boxes {
        scene.add_box(b);
    }
    let h = ((scene.max - scene.min).length() / 128.0).max(1e-9);
    let cell_of = |p: crate::math::Point| [(p.x / h).floor() as i64, (p.y / h).floor() as i64, (p.z / h).floor() as i64];
    let cells_of_box = |b: &BoundingBox| {
        let (lo, hi) = (cell_of(b.min), cell_of(b.max));
        (lo[0]..=hi[0]).flat_map(move |i| (lo[1]..=hi[1]).flat_map(move |j| (lo[2]..=hi[2]).map(move |k| [i, j, k])))
    };
    let tol = 1e-6;
    let mut kept: Vec<(Shape, crate::geometry::curve::Curve, f64, f64)> = Vec::new();
    let mut grid: FastMap<[i64; 3], Vec<usize>> = FastMap::default();
    let mut ends: FastSet<NodeId> = FastSet::default();
    for (_, e) in edges {
        let Some((c, first, last)) = model.edge_curve(e) else { continue };
        let p = c.point(0.5 * (first + last));
        if !on_two(model, p)? {
            continue;
        }
        // Once: a stretch already kept (another cell's copy, or a whole solid's edge its faces'
        // pieces split) is passed over.
        let covered = grid.get(&cell_of(p)).is_some_and(|near| {
            near.iter().any(|&i| {
                let (_, kc, kf, kl) = &kept[i];
                geom::project(kc, kf.min(*kl), kf.max(*kl), p).1 <= tol
            })
        });
        if covered {
            continue;
        }
        let b = crate::properties::bounding_box(model, e).enlarged(tol);
        for g in cells_of_box(&b) {
            grid.entry(g).or_default().push(kept.len());
        }
        for v in model.unique_subshapes(e, NodeKind::Vertex) {
            ends.insert(v.node);
        }
        kept.push((e, c, first, last));
    }
    // The points where they touch elsewhere (an edge piercing a face): vertices on two
    // arguments ending no edge of the section.
    let mut points: Vec<Shape> = Vec::new();
    let mut seen_v = FastSet::default();
    for cell in pieces {
        for v in model.unique_subshapes(cell.shape, NodeKind::Vertex) {
            if !ends.contains(&v.node) && seen_v.insert(v.node) {
                let p = model.point(v);
                if on_two(model, p)? {
                    points.push(v);
                }
            }
        }
    }
    // Each edge split at the section's vertices lying inside it.
    let mut vertex_grid: FastMap<[i64; 3], Vec<Shape>> = FastMap::default();
    let mut all_vertices: Vec<Shape> = Vec::new();
    let mut seen_end = FastSet::default();
    for (e, ..) in &kept {
        for v in model.unique_subshapes(*e, NodeKind::Vertex) {
            if seen_end.insert(v.node) {
                all_vertices.push(v);
            }
        }
    }
    all_vertices.extend(points.iter().copied());
    for &v in &all_vertices {
        vertex_grid.entry(cell_of(model.point(v))).or_default().push(v);
    }
    for (e, c, first, last) in kept {
        let b = crate::properties::bounding_box(model, e).enlarged(tol);
        let (a, z) = model.edge_vertices(e);
        let mut inner: Vec<(f64, Shape)> = Vec::new();
        for g in cells_of_box(&b) {
            for &v in vertex_grid.get(&g).map(Vec::as_slice).unwrap_or_default() {
                if v.node == a.node || v.node == z.node || inner.iter().any(|(_, w)| w.node == v.node) {
                    continue;
                }
                let found = geom::project(&c, first.min(last), first.max(last), model.point(v));
                let margin = 1e-9 * (last - first).abs();
                if found.1 <= tol && found.0 > first.min(last) + margin && found.0 < first.max(last) - margin {
                    inner.push((found.0, v));
                }
            }
        }
        if inner.is_empty() {
            model.add(result, e)?;
            continue;
        }
        inner.sort_by(|x, y| x.0.total_cmp(&y.0));
        let (lo, hi) = (first.min(last), first.max(last));
        let (start, end) = if first <= last { (a, z) } else { (z, a) };
        let mut stops = vec![(lo, start)];
        stops.extend(inner);
        stops.push((hi, end));
        let tolerance = model.tolerance(e);
        for w in stops.windows(2) {
            let ((t0, v0), (t1, v1)) = (w[0], w[1]);
            let (v0, v1) = (crate::topology::VertexId::try_from(v0.node).expect("a vertex"), crate::topology::VertexId::try_from(v1.node).expect("a vertex"));
            let piece = model.add_edge(v0, v1, tolerance)?;
            model.set_edge_curve(piece.into(), c.clone(), tolerance)?;
            model.set_edge_range(piece.into(), t0, t1)?;
            model.add(result, piece.into())?;
        }
    }
    for v in points {
        model.add(result, v)?;
    }
    Ok(result)
}

/// Whether `p` lies on argument `shape`: on a solid's boundary (by its classifier), on a
/// shell's, a face's or a wire's self (within `tolerance`).
fn on_argument(model: &mut Model, classifier: &Option<PointClassifier>, shape: Shape, p: crate::math::Point, tolerance: f64) -> Result<bool> {
    Ok(match classifier {
        Some(c) => c.locate(p) == PointLocation::Boundary,
        None => {
            // (A face on an elementary surface farther from the point than the tolerance,
            // as its surface is: off it, without measuring.)
            if shape.kind() == NodeKind::Face {
                use crate::geometry::surface::Surface;
                if let Some(surface @ (Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_))) = model.surface(shape) {
                    if geom::signed_distance(surface, p).abs() > 2.0 * tolerance {
                        return Ok(false);
                    }
                }
            }
            let v = crate::make::vertex(model, p);
            crate::distance::distance(model, v, shape)?.value <= tolerance
        }
    })
}

/// The largest tolerance of `shape`'s faces, edges and vertices.
fn max_tolerance(model: &Model, shape: Shape) -> f64 {
    [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex].iter().flat_map(|&k| model.unique_subshapes(shape, k)).map(|s| model.tolerance(s)).fold(0.0, f64::max)
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

    #[test]
    fn arguments_of_mixed_dimensions() {
        let mut model = Model::new();
        let block = |model: &mut Model, a: [f64; 3], b: [f64; 3]| Cuboid::from_corners(Point::new(a[0], a[1], a[2]), Point::new(b[0], b[1], b[2])).build(model).unwrap();
        // A solid and the shell of an overlapping box as objects, a box as the tool: the
        // cut keeps the shell's parts outside the tool, and the solid's, split by the shell
        // (OCCT's counts, the shell first).
        let s1 = block(&mut model, [0.0; 3], [10.0; 3]);
        let s2 = block(&mut model, [5.0, 0.0, 0.0], [15.0, 10.0, 10.0]);
        let shell = model.explore(s2, NodeKind::Shell).next().unwrap();
        let tool = block(&mut model, [5.0, 0.0, 5.0], [15.0, 10.0, 15.0]);
        let cut = combine(&mut model, &[s1, shell], &[tool], Combination::Difference).unwrap();
        let kinds: Vec<NodeKind> = model.children(cut).map(|c| c.kind()).collect();
        assert_eq!(kinds, [NodeKind::Shell, NodeKind::Solid, NodeKind::Solid]);
        let v = crate::properties::volume(&model, model.explore(cut, NodeKind::Solid).collect::<Vec<_>>()[0]).mass + crate::properties::volume(&model, model.explore(cut, NodeKind::Solid).collect::<Vec<_>>()[1]).mass;
        assert!((v - 750.0).abs() < 1e-9, "{v}");
        // Not defined: a union of different dimensions, a solid cut by a face.
        let face = model.explore(s2, NodeKind::Face).next().unwrap();
        assert!(combine(&mut model, &[s1], &[face], Combination::Union).is_err());
        assert!(combine(&mut model, &[s1], &[face], Combination::Difference).is_err());
        // A face across a solid: their common part is the face's inside the solid; an edge
        // through the solid meets its boundary at two points.
        let plane = crate::geometry::surface::Surface::Plane(crate::geometry::surface::Plane::new(Point::new(5.0, 5.0, 5.0), crate::math::Direction::X));
        let across = crate::make::face_on_domain(&mut model, plane, -20.0, 20.0, -20.0, 20.0).unwrap();
        let common = combine(&mut model, &[s1], &[across], Combination::Intersection).unwrap();
        assert!((crate::properties::area(&model, common).mass - 100.0).abs() < 1e-9);
        let through = crate::make::edge_between(&mut model, Point::new(-5.0, 5.0, 5.0), Point::new(15.0, 5.0, 5.0)).unwrap();
        let section = combine(&mut model, &[s1], &[through], Combination::Section).unwrap();
        assert_eq!(counts(&model, section)[..2], [2, 0]);
    }
}
