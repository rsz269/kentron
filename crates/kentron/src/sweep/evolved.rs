//! Evolved solids: a profile swept along a planar spine, keeping its place beside the spine
//! (its offset from the spine and its height), turning about the spine's corners, the
//! self-intersections removed.
//!
//! As OCCT's evolved shapes in volume mode (BRepFill_AdvancedEvolved): the profile is swept
//! along the spine from where it stands (its nearest point on the spine), the sweep turned
//! about the spine's corners (a rounded transition); the sweep's free edges level with the
//! spine bound planar lids; the space the sweep and the lids enclose, split where they
//! cross, is the solid (each closed piece of it, joined).

use crate::error::{Error, Result};
use crate::geometry::extrema::nearest_on_curve;
use crate::geometry::ParametricCurve;
use crate::math::{Point, Point2, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::{Corner, PipeOptions};

/// How [`evolved`] finishes its result.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EvolvedOptions {
    /// A solid (the sweep closed by its lids); else the sweep's faces alone, as a shell.
    pub solid: bool,
}

/// `profile` (a wire) swept along `spine` (a planar wire or face), keeping its place beside
/// the spine, turning about the spine's corners: the solid the sweep and its lids enclose,
/// or (unless `options.solid`) that solid's faces but the lids.
#[doc(alias = "BRepOffsetAPI_MakeEvolved")]
#[doc(alias = "BRepFill_AdvancedEvolved")]
pub fn evolved(model: &mut Model, spine: Shape, profile: Shape, options: EvolvedOptions) -> Result<Shape> {
    let spine = match spine.kind() {
        NodeKind::Wire => spine,
        NodeKind::Face => *model.unique_subshapes(spine, NodeKind::Wire).first().ok_or_else(|| Error::InvalidInput("an evolved shape's spine face has no wire".into()))?,
        kind => return Err(Error::InvalidInput(format!("an evolved shape's spine is a wire or a face, not a {kind:?}"))),
    };
    if profile.kind() != NodeKind::Wire {
        return Err(Error::InvalidInput("an evolved shape's profile is a wire".into()));
    }
    let plane = crate::make::plane_of_wire(model, spine)?;
    let normal = Vector::from(plane.normal());
    let spine = starting_by(model, spine, profile)?;
    let swept = super::pipe_cornered(model, profile, spine, &PipeOptions::default(), Corner::Rounded)?;
    // The lids: faces the sweep's free edges level with the spine bound.
    let size = crate::properties::bounding_box(model, swept).size().length().max(1e-9);
    let height = |p: Point| normal.dot(p - Point::ORIGIN);
    let level = |model: &Model, e: Shape| -> bool {
        let Some((curve, a, b)) = model.edge_curve(e) else { return false };
        let heights: Vec<f64> = (0..=4).map(|k| height(curve.point(a + (b - a) * k as f64 / 4.0))).collect();
        heights.iter().all(|h| (h - heights[0]).abs() <= 1e-7 * size)
    };
    // (Each edge's use by the sweep's faces: a free edge's one use, as its face runs it.)
    let mut uses: std::collections::HashMap<NodeId, (Shape, usize)> = std::collections::HashMap::new();
    for face in model.explore(swept, NodeKind::Face).collect::<Vec<_>>() {
        for e in model.explore(face, NodeKind::Edge).collect::<Vec<_>>() {
            uses.entry(e.node).or_insert((e, 0)).1 += 1;
        }
    }
    let free: Vec<Shape> = uses.values().filter(|&&(e, count)| count == 1 && !model.is_degenerate(e) && level(model, e)).map(|&(e, _)| e).collect();
    // The free edges in the spine's plane, as they run: polylines in its coordinates, by height.
    let frame = *plane.frame();
    let (x, y) = (Vector::from(frame.x()), Vector::from(frame.y()));
    let flat = |p: Point| Point2::new(x.dot(p - frame.origin()), y.dot(p - frame.origin()));
    let mut runs: Vec<(f64, Vec<Point2>)> = Vec::new();
    for &e in &free {
        let (curve, a, b) = model.edge_curve(e).expect("a level edge has a curve");
        let (a, b) = if e.orientation == Orientation::Reversed { (b, a) } else { (a, b) };
        runs.push((height(curve.point(a)), (0..=64).map(|k| flat(curve.point(a + (b - a) * k as f64 / 64.0))).collect()));
    }
    // Each lid piece as often as the free edges at its height wind round it, turned to close
    // the sweep (its edges run the other way to theirs).
    let counted = model.new_empty(NodeKind::Compound);
    model.add(counted, swept)?;
    if !free.is_empty() {
        let edges = model.new_empty(NodeKind::Compound);
        for &e in &free {
            model.add(edges, Shape { orientation: Orientation::Forward, ..e })?;
        }
        let lids = crate::make::faces_from_edges(model, edges, false)?;
        for lid in model.unique_subshapes(lids, NodeKind::Face) {
            let Some(q) = interior_point(model, lid) else { continue };
            let at = height(q);
            let around: Vec<&Vec<Point2>> = runs.iter().filter(|(h, _)| (h - at).abs() <= 1e-6 * size).map(|(_, r)| r).collect();
            let w = winding(&around, flat(q));
            if w == 0 {
                continue;
            }
            // (Its plane's normal along the spine's, or against it.)
            let up = matches!(model.surface(lid), Some(crate::geometry::surface::Surface::Plane(p)) if Vector::from(p.normal()).dot(normal) > 0.0) == (lid.orientation == Orientation::Forward);
            let lid = if (w < 0) == up { lid } else { lid.reversed() };
            model.add(counted, lid)?;
            for _ in 1..w.unsigned_abs() {
                let copy = crate::modify::copy(model, lid, false)?;
                model.add(counted, copy)?;
            }
        }
    }
    let faces: Vec<Shape> = model.unique_subshapes(counted, NodeKind::Face);
    let cells = crate::boolean::enclosed_solids(model, &faces)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    if cells.is_empty() {
        return Err(Error::Failed("an evolved shape enclosing nothing".into()));
    }
    // The cells the sweep and its lids wind round (taken the way the largest is wound).
    let classifier = crate::boolean::PointClassifier::new(model, counted)?;
    let mut wound = Vec::with_capacity(cells.len());
    for &cell in &cells {
        let w = match crate::boolean::point_inside(model, cell)? {
            Some(p) => classifier.winding(p),
            None => None,
        };
        wound.push((cell, w, crate::properties::volume(model, cell).mass.abs()));
    }
    let sign = wound.iter().filter(|(_, w, _)| w.is_some_and(|w| w != 0)).max_by(|a, b| a.2.total_cmp(&b.2)).map_or(1, |(_, w, _)| w.unwrap().signum());
    let cells: Vec<Shape> = wound.iter().filter(|(_, w, _)| w.is_some_and(|w| sign * w >= 1)).map(|&(c, _, _)| c).collect();
    if cells.is_empty() {
        return Err(Error::Failed("an evolved shape enclosing nothing".into()));
    }
    let solid = joined(model, &cells)?;
    if options.solid {
        return Ok(solid);
    }
    // The faces but the lids (planar, level with the spine).
    let shell = model.new_empty(NodeKind::Shell);
    for face in model.unique_subshapes(solid, NodeKind::Face) {
        let lid = matches!(model.surface(face), Some(crate::geometry::surface::Surface::Plane(p)) if Vector::from(p.normal()).cross(normal).length() <= 1e-9);
        if !lid {
            model.add(shell, face)?;
        }
    }
    Ok(shell)
}

/// `spine` starting where `profile` stands on it (its point nearest the profile), where it is
/// closed: the sweep runs on round from there. An open spine as it is.
fn starting_by(model: &mut Model, spine: Shape, profile: Shape) -> Result<Shape> {
    let edges = model.ordered_edges(spine);
    let (first_vertex, _) = model.edge_vertices(edges[0]);
    let (_, last_vertex) = model.edge_vertices(edges[edges.len() - 1]);
    let closed = first_vertex.node == last_vertex.node || model.point(first_vertex).distance(model.point(last_vertex)) <= model.tolerance(first_vertex).max(model.tolerance(last_vertex));
    if !closed {
        return Ok(spine);
    }
    // Each spine edge's curve, as the spine runs (`from` to `to`).
    let mut runs = Vec::with_capacity(edges.len());
    for &e in &edges {
        let (curve, a, b) = model.edge_curve(e).ok_or_else(|| Error::InvalidInput("an evolved shape's spine edge has no curve".into()))?;
        runs.push(if e.orientation == Orientation::Reversed { (curve, b, a) } else { (curve, a, b) });
    }
    // The profile's points (its edges sampled), and the spine's point nearest them.
    let mut points = Vec::new();
    for e in model.unique_subshapes(profile, NodeKind::Edge) {
        if let Some((curve, a, b)) = model.edge_curve(e) {
            points.extend((0..=32).map(|k| curve.point(a + (b - a) * k as f64 / 32.0)));
        }
    }
    let mut best = (f64::INFINITY, 0, 0.0);
    for (k, (curve, from, to)) in runs.iter().enumerate() {
        for &p in &points {
            let near = nearest_on_curve(p, curve, from.min(*to), from.max(*to));
            if near.distance < best.0 {
                best = (near.distance, k, near.parameter);
            }
        }
    }
    let (_, k, t) = best;
    let scale = crate::properties::bounding_box(model, spine).size().length().max(1e-9);
    let close = |a: f64, b: f64, (curve, _, _): &(crate::geometry::curve::Curve, f64, f64)| curve.point(a).distance(curve.point(b)) <= 1e-9 * scale;
    // The edges from there round: the rest of edge k, the others, edge k up to there.
    let mut order: Vec<(usize, f64, f64)> = Vec::new();
    let (from, to) = (runs[k].1, runs[k].2);
    if !close(t, to, &runs[k]) {
        order.push((k, t, to));
    }
    for i in 1..runs.len() {
        let j = (k + i) % runs.len();
        order.push((j, runs[j].1, runs[j].2));
    }
    if !close(from, t, &runs[k]) {
        order.push((k, from, t));
    }
    if order.len() == edges.len() && order[0].0 == 0 {
        return Ok(spine);
    }
    let mut made = Vec::with_capacity(order.len());
    for (j, a, b) in order {
        let curve = runs[j].0.clone();
        let edge = crate::make::edge_from_curve_range(model, curve, a.min(b), a.max(b))?;
        made.push(if a > b { edge.reversed() } else { edge });
    }
    crate::make::wire(model, &made)
}

/// The cells joined: the faces two of them share are inside; faces left in pieces on one
/// surface merged.
fn joined(model: &mut Model, cells: &[Shape]) -> Result<Shape> {
    let mut uses: std::collections::HashMap<NodeId, (Shape, usize)> = std::collections::HashMap::new();
    let mut order = Vec::new();
    for &cell in cells {
        for face in model.explore(cell, NodeKind::Face).collect::<Vec<_>>() {
            let entry = uses.entry(face.node).or_insert_with(|| {
                order.push(face.node);
                (face, 0)
            });
            entry.1 += 1;
        }
    }
    let shell = model.new_empty(NodeKind::Shell);
    for node in order {
        let (face, count) = uses[&node];
        if count == 1 {
            model.add(shell, face)?;
        }
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    crate::unify::unify_same_domain_with(model, solid, crate::unify::UnifyOptions { same_surface_only: true, ..Default::default() })
}

/// How many times the polylines (closed loops, cut anywhere into runs) wind round `q`.
fn winding(runs: &[&Vec<Point2>], q: Point2) -> i32 {
    let mut w = 0;
    for run in runs {
        for pair in run.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let side = (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
            if a.y <= q.y && b.y > q.y && side > 0.0 {
                w += 1;
            } else if a.y > q.y && b.y <= q.y && side < 0.0 {
                w -= 1;
            }
        }
    }
    w
}

/// A point inside the planar `face`: just inside the middle of a side of its outer loop.
fn interior_point(model: &Model, face: Shape) -> Option<Point> {
    use crate::geometry::ParametricSurface;
    let loops = crate::heal::face_loops(model, Shape { orientation: Orientation::Forward, ..face })?;
    let polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect();
    let outer = polygons.first()?;
    let size = outer.iter().fold(0.0f64, |m, p| m.max((*p - outer[0]).length()));
    // (Closed: each loop's first point again at its end, for `winding`'s runs.)
    let closed: Vec<Vec<Point2>> = polygons.iter().map(|p| p.iter().copied().chain(std::iter::once(p[0])).collect()).collect();
    let closed_refs: Vec<&Vec<Point2>> = closed.iter().collect();
    for i in 0..outer.len() {
        let (a, b) = (outer[i], outer[(i + 1) % outer.len()]);
        let along = b - a;
        if along.length() <= 1e-12 * size {
            continue;
        }
        let k = 1e-4 * size / along.length();
        for sign in [1.0, -1.0] {
            let q = Point2::new(0.5 * (a.x + b.x) - sign * along.y * k, 0.5 * (a.y + b.y) + sign * along.x * k);
            if winding(&closed_refs, q) != 0 {
                return Some(model.surface(face)?.point(q.x, q.y));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn polyline(model: &mut Model, points: &[[f64; 3]]) -> Shape {
        let edges: Vec<Shape> = points.windows(2).map(|w| crate::make::edge_between(model, Point::new(w[0][0], w[0][1], w[0][2]), Point::new(w[1][0], w[1][1], w[1][2])).unwrap()).collect();
        crate::make::wire(model, &edges).unwrap()
    }

    /// OCCT's evolved/voluved AGT002: a closed polygonal spine (a notched rectangle) and a
    /// closed profile standing on it, leaning in: its volume and area as OCCT's.
    #[test]
    fn a_profile_swept_round_a_notched_rectangle() {
        let mut model = Model::new();
        let spine = polyline(&mut model, &[[0.0, 0.0, 0.0], [0.0, 10.0, 0.0], [20.0, 10.0, 0.0], [20.0, 0.0, 0.0], [13.0, 0.0, 0.0], [13.0, 5.0, 0.0], [6.0, 5.0, 0.0], [6.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let profile = polyline(&mut model, &[[-2.0, 2.0, 0.0], [0.0, 2.0, 5.0], [0.0, 2.0, 10.0], [1.0, 2.0, 10.0], [1.0, 2.0, 0.0], [-2.0, 2.0, 0.0]]);
        let solid = evolved(&mut model, spine, profile, EvolvedOptions { solid: true }).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid());
        let volume = crate::properties::volume(&model, solid).mass;
        let area = crate::properties::area(&model, solid).mass;
        assert!((volume - 1023.79).abs() < 1e-2, "{volume}");
        assert!((area - 1649.58).abs() < 1e-2, "{area}");
    }
}
