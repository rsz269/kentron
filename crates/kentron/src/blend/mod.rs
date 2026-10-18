//! Rounding and bevelling edges of solids: fillets and chamfers.

mod ball;
mod corner2d;
mod runover;

pub use corner2d::{chamfer_2d, fillet_2d, fillet_wire_2d, Corner2d};

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::surface::Surface;
use std::f64::consts::TAU;

use crate::geometry::ParametricSurface;
use crate::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeId, NodeKind, Orientation, Shape};
use crate::{boolean, make, sweep, tolerance};

/// `shape` with `edges` rounded by a rolling ball of `radius`.
///
/// Supported where the fillet is a cylinder or a torus: straight edges between planes, or
/// between a plane and a cylinder parallel to the edge; circular edges between faces
/// coaxial with them (planes, cylinders, cones, spheres, tori). No two of the edges may
/// share a vertex (no corner blends yet). The fillet is made by cutting away (at a convex
/// edge) or adding (at a concave one) the solid swept by the region between the edge and
/// the ball's arc in the edge's cross-section. A cut runs past the edge's ends so that the
/// faces there trim it; added material stops at them (exact where the end faces are square
/// to the edge). The result is a compound of solids, like a boolean operation's.
///
/// As in OCCT, the rounding runs on along edges continuing an edge tangentially (a
/// rounded outline's lines and arcs): those are filleted too, their tools meeting exactly
/// at the joints.
#[doc(alias = "BRepFilletAPI_MakeFillet")]
pub fn fillet(model: &mut Model, shape: Shape, edges: &[Shape], radius: f64) -> Result<Shape> {
    let pairs: Vec<(Shape, f64)> = edges.iter().map(|&e| (e, radius)).collect();
    fillet_edges(model, shape, &pairs)
}

/// [`fillet`] with a radius per edge.
#[doc(alias = "BRepFilletAPI_MakeFillet::Add")]
pub fn fillet_edges(model: &mut Model, shape: Shape, edges: &[(Shape, f64)]) -> Result<Shape> {
    let edges = &sharp(model, shape, edges, |e| e.0)?;
    on_solids(model, shape, |model, shape| each_solid(model, shape, edges, |e| e.0, fillet_solid_edges))
}

/// [`fillet`] of all `edges`, those where the faces meet smoothly too (as an offset's
/// parting edges, which may meet so but for rounding).
pub(crate) fn fillet_all(model: &mut Model, shape: Shape, edges: &[Shape], radius: f64) -> Result<Shape> {
    let edges: Vec<(Shape, f64)> = edges.iter().map(|&e| (e, radius)).collect();
    on_solids(model, shape, |model, shape| each_solid(model, shape, &edges, |e| e.0, fillet_solid_edges))
}

/// `op` on each solid of `shape` with the edges of `edges` it has (by node and placement:
/// one solid placed twice has the same edges twice, each pair apart), the solids without
/// any left as they are; `op` on `shape` itself where it has one solid, or more than solids.
fn each_solid<T: Copy>(model: &mut Model, shape: Shape, edges: &[T], edge: impl Fn(&T) -> Shape, op: impl Fn(&mut Model, Shape, &[T]) -> Result<Shape>) -> Result<Shape> {
    let solids: Vec<Shape> = model.explore(shape, NodeKind::Solid).collect();
    let shells_in = |model: &Model, s: Shape| model.explore(s, NodeKind::Shell).count();
    let only_solids = shells_in(model, shape) == solids.iter().map(|&s| shells_in(model, s)).sum::<usize>();
    if solids.len() < 2 || !only_solids {
        return op(model, shape, edges);
    }
    let mut pieces = Vec::new();
    for &solid in &solids {
        let own: Vec<Shape> = model.explore(solid, NodeKind::Edge).collect();
        let mine: Vec<T> = edges.iter().copied().filter(|e| own.iter().any(|o| o.is_same(&edge(e)))).collect();
        if mine.is_empty() {
            pieces.push(solid);
        } else {
            let done = op(model, solid, &mine)?;
            pieces.extend(model.explore(done, NodeKind::Solid));
        }
    }
    let out = model.new_empty(NodeKind::Compound);
    for piece in pieces {
        model.add(out, piece)?;
    }
    Ok(out)
}

/// Whether the faces at `edge` (of `shape`) meet smoothly there: a seam (a face meeting
/// itself), or two faces tangent along it. As OCCT, such edges are left out of fillets and
/// chamfers (there is no corner to round).
fn smooth(model: &Model, shape: Shape, edge: Shape) -> bool {
    use crate::geometry::ParametricCurve;
    let mut uses: Vec<NodeId> = Vec::new();
    for face in model.explore(shape, NodeKind::Face) {
        for e in model.explore(face, NodeKind::Edge) {
            if e.node == edge.node {
                uses.push(face.node);
            }
        }
    }
    match uses[..] {
        [a, b] if a == b => true,
        [a, b] => {
            let Some((curve, first, last)) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }) else { return false };
            [0.2, 0.5, 0.8].iter().all(|&s| {
                let p = curve.point(first + (last - first) * s);
                match (normal_at(model, shape, a, p), normal_at(model, shape, b, p)) {
                    (Some(n), Some(m)) => n.dot(m).abs() >= 1.0 - tolerance::ANGULAR.sqrt(),
                    _ => false,
                }
            })
        }
        _ => false,
    }
}

/// `edges` less those where the faces meet smoothly; an error if that leaves none.
fn sharp<T: Copy>(model: &Model, shape: Shape, edges: &[T], edge: impl Fn(&T) -> Shape) -> Result<Vec<T>> {
    let kept: Vec<T> = edges.iter().copied().filter(|e| !smooth(model, shape, edge(e))).collect();
    if kept.is_empty() && !edges.is_empty() {
        return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
    }
    Ok(kept)
}

/// [`fillet_edges`] of a shape with solids.
fn fillet_solid_edges(model: &mut Model, shape: Shape, edges: &[(Shape, f64)]) -> Result<Shape> {
    let rounded = fillet_solid_edges_as(model, shape, edges, true)?;
    // (Where the cylinder and torus fillets' tools leave a face's boundary crossing itself,
    // as where a fillet runs over another feature of its faces: all rolled instead, which
    // goes round such features.)
    // (Either taken only where it leaves the solid as sound as it was.)
    let before = fault_kinds(model, shape);
    let crossing = crate::validity::ProblemKind::SelfIntersection;
    if !before.contains(&crossing) && crosses_itself(model, rounded) {
        if let Ok(rolled) = fillet_solid_edges_as(model, shape, edges, false).and_then(|r| sound(model, &before, r)) {
            return Ok(rolled);
        }
    }
    Ok(rounded)
}

/// [`fillet_solid_edges`], the cylinder and torus fillets made as tools for boolean
/// operations where `elementary`, else rolled as the others.
fn fillet_solid_edges_as(model: &mut Model, shape: Shape, edges: &[(Shape, f64)], elementary: bool) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    let edges = &tangent_chains(model, shape, edges);
    // The edges' geometry and adjacent faces, all taken from the original shape. Edges
    // whose fillet is a cylinder or a torus get a tool for a boolean operation; the others
    // a rolling-ball fillet, built first (it keeps the other edges).
    let mut classified: Vec<(Shape, f64, Option<Corner>)> = Vec::new();
    for &(edge, radius) in edges {
        if !(radius > tolerance::LINEAR) {
            return Err(Error::InvalidInput(format!("fillet radius must be positive, got {radius}")));
        }
        let corner = Corner::of(model, shape, edge, false)?.filter(|_| elementary && elementary_around(model, shape, edge));
        classified.push((edge, radius, corner));
    }
    // The edges meeting at their ends, in groups. A group with a rolling-ball edge is
    // rounded edge by edge, each fillet ending on the faces (and fillets) at its ends, where
    // that succeeds; otherwise, and for groups of cylinder and torus fillets, the fillets
    // are tools for boolean operations meeting as they cut each other.
    let ends_of = |model: &Model, edge: Shape| -> Vec<Point> {
        let (a, b) = model.edge_vertices(edge);
        if a.node == b.node { Vec::new() } else { vec![model.point(a), model.point(b)] }
    };
    let n = classified.len();
    let mut group: Vec<usize> = (0..n).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    for i in 0..n {
        for j in i + 1..n {
            let (pi, pj) = (ends_of(model, classified[i].0), ends_of(model, classified[j].0));
            if pi.iter().any(|p| pj.iter().any(|q| p.distance(*q) <= tolerance::LINEAR)) {
                let (a, b) = (root(&mut group, i), root(&mut group, j));
                group[a.max(b)] = a.min(b);
            }
        }
    }
    let placement = shape.placement;
    if classified.iter().any(|(e, _, c)| c.is_none() && e.placement != placement) {
        todo!("rolling-ball fillets of edges placed apart from their shape");
    }
    // (Radii in the shape's own units, where its placement scales it.)
    let unscale = 1.0 / placement.scale_factor().abs();
    let mut local = Shape { placement: Transform::IDENTITY, ..shape };
    let mut corners = Vec::new();
    let mut rolling_tools = Vec::new();
    let roots: Vec<usize> = { let mut r: Vec<usize> = (0..n).map(|i| root(&mut group, i)).collect(); r.sort_unstable(); r.dedup(); r };
    for r in roots {
        let members: Vec<usize> = (0..n).filter(|&i| root(&mut group, i) == r).collect();
        if members.iter().all(|&i| classified[i].2.is_some()) {
            // An edge alone whose ball would touch its faces off their bounds rolls on onto
            // the faces beyond (where that can be built).
            if let [i] = members[..] {
                let (edge, radius, corner) = &classified[i];
                let corner = corner.as_ref().expect("an elementary fillet");
                let local_edge = Shape { placement: Transform::IDENTITY, ..*edge };
                // (Along a straight convex edge, from the corner's own ball.)
                let rolls = match corner.path {
                    Path::Line { .. } if corner.convex => rolls_off(model, local, corner, *radius),
                    _ => runover::needed(model, local, local_edge, radius * unscale),
                };
                if rolls {
                    let before = fault_kinds(model, local);
                    if let Ok(rounded) = runover::fillet(model, local, local_edge, radius * unscale).and_then(|r| sound(model, &before, r)) {
                        local = rounded;
                        continue;
                    }
                }
            }
            for &i in &members {
                let (_, radius, corner) = &classified[i];
                corners.push((corner.clone().expect("an elementary fillet"), *radius));
            }
            continue;
        }
        // Edge by edge on the solid as it is rebuilt (the edges found again by their
        // geometry: a fillet shortens the edges at its ends).
        let radii: Vec<f64> = members.iter().map(|&i| classified[i].1 * unscale).collect();
        let laws: Vec<Box<dyn Fn(f64) -> f64>> = radii.iter().map(|&r| Box::new(move |_| r) as Box<dyn Fn(f64) -> f64>).collect();
        let edges: Vec<(Shape, &dyn Fn(f64) -> f64)> = members.iter().zip(&laws).map(|(&i, law)| (Shape { placement: Transform::IDENTITY, ..classified[i].0 }, law.as_ref())).collect();
        // (Each way taken only where it leaves the solid as sound as it was.)
        let before = fault_kinds(model, local);
        // (One edge whose ball touches the faces beyond: rolling on onto them first.)
        let mut rolled_over = false;
        if let [i] = members[..] {
            let edge = Shape { placement: Transform::IDENTITY, ..classified[i].0 };
            if runover::needed(model, local, edge, radii[0]) {
                rolled_over = true;
                if let Ok(rounded) = runover::fillet(model, local, edge, radii[0]).and_then(|r| sound(model, &before, r)) {
                    local = rounded;
                    continue;
                }
            }
        }
        let roll = |model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)]| {
            ball::fillet_network(model, shape, edges).or_else(|_| fillet_in_turn(model, shape, edges, rolled_over))
        };
        let rounded = roll(model, local, &edges).or_else(|e| match &members[..] {
            // One edge whose ball rolls off its faces: on onto the faces beyond (unless that
            // was tried first, on the same shape).
            [i] if !rolled_over => runover::fillet(model, local, Shape { placement: Transform::IDENTITY, ..classified[*i].0 }, radii[0]).map_err(|_| e),
            _ => Err(e),
        });
        // (Where a face's boundary now crosses itself, the fillets run over another feature
        // on their faces: rounded without it, the feature then put back; or one runs off its
        // faces: on onto the faces beyond.)
        let crossing = crate::validity::ProblemKind::SelfIntersection;
        let rounded = match rounded {
            Ok(r) if !before.contains(&crossing) && fault_kinds(model, r).contains(&crossing) => over_features(model, local, r, &edges, &roll, 3).or_else(|e| match &members[..] {
                [i] if !rolled_over => runover::fillet(model, local, Shape { placement: Transform::IDENTITY, ..classified[*i].0 }, radii[0]),
                _ => Err(e),
            }),
            other => other,
        };
        let rounded = rounded.and_then(|r| sound(model, &before, r));
        match rounded {
            Ok(rounded) => local = rounded,
            Err(_) => {
                // As tools instead (made on the shape as it was).
                for &i in &members {
                    let (edge, radius, corner) = &classified[i];
                    match corner {
                        Some(c) => corners.push((c.clone(), *radius)),
                        None => {
                            let (tool, convex) = ball::rolling_tool(model, Shape { placement: Transform::IDENTITY, ..shape }, Shape { placement: Transform::IDENTITY, ..*edge }, &|_| radius * unscale)?;
                            rolling_tools.push((tool.moved(&placement), convex));
                        }
                    }
                }
            }
        }
    }
    let shape = local.moved(&placement);
    let corners = merge_arcs(corners);
    // Where filleted edges meet: two make a mitre (the fillets cut one after the other meet
    // along their intersection); three need a corner patch.
    let mut meetings: Vec<(Point, Vec<usize>)> = Vec::new();
    for (i, (corner, _)) in corners.iter().enumerate() {
        for p in corner.ends() {
            match meetings.iter_mut().find(|(q, _)| q.distance(p) <= tolerance::LINEAR) {
                Some((_, at)) => at.push(i),
                None => meetings.push((p, vec![i])),
            }
        }
    }
    let mut patches = Vec::new();
    for (p, at) in &meetings {
        match at.len() {
            1 => {}
            // (Concave fillets, or more than three, are left meeting as they cut each other.)
            3 => patches.extend(CornerPatch::of(&corners, at, *p)?),
            _ => {}
        }
    }
    let mut result = shape;
    // (A step leaving the solid's shells open, where they were closed, and the solid in
    // pieces: the steps after it split those pieces by one another, and by the open edges
    // ever further, and the fillet fails there. A step only leaving them open may be undone
    // by the next, which takes the open part away.)
    let closed = boolean::closed_shells(model, shape);
    let keep = |model: &Model, result: Shape| -> Result<Shape> {
        if closed && model.explore(result, NodeKind::Solid).nth(1).is_some() && !boolean::closed_shells(model, result) {
            return Err(Error::Failed("a fillet's tool left the solid's shells open".into()));
        }
        Ok(result)
    };
    let mitred = mitres(&corners);
    for (((corner, radius), extend), mitre) in corners.iter().zip(extensions(model, shape, &corners)).zip(&mitred) {
        let tool = corner.tool(model, *radius, mitre_runs(corner, *radius, extend, mitre))?;
        let tool = mitred_tool(model, tool, mitre)?;
        result = if corner.convex { boolean::difference(model, result, tool)? } else { boolean::union(model, result, tool)? };
        result = keep(model, result)?;
    }
    for (tool, convex) in rolling_tools {
        result = if convex { boolean::difference(model, result, tool)? } else { boolean::union(model, result, tool)? };
        result = keep(model, result)?;
    }
    for patch in patches {
        let tool = patch.tool(model)?;
        result = boolean::difference(model, result, tool)?;
        result = keep(model, result)?;
    }
    Ok(result)
}

/// Whether `corner`'s fillet of `radius` (along a straight edge, convex) would touch its
/// faces off their bounds: where the ball must roll on onto the faces beyond.
fn rolls_off(model: &Model, shape: Shape, corner: &Corner, radius: f64) -> bool {
    let Path::Line { length } = corner.path else { return false };
    if !corner.convex {
        return false;
    }
    let Ok((_, t1, t2)) = corner.ball(radius) else { return false };
    let faces = model.unique_subshapes(shape, NodeKind::Face);
    let find = |node: NodeId| faces.iter().copied().find(|f| f.node == node);
    let (Some(a), Some(b)) = (find(corner.faces[0].0), find(corner.faces[1].0)) else { return false };
    runover::rolls_off(model, [a, b], corner.tangent, length, [t1, t2])
}

/// `op` on `shape`, or where `shape` has no solids but closed shells (as revolving a closed
/// wire makes), on those shells made solids, the result's shells given back (one shell, or
/// a compound of them): OCCT rounds shells as well.
fn on_solids(model: &mut Model, shape: Shape, op: impl FnOnce(&mut Model, Shape) -> Result<Shape>) -> Result<Shape> {
    if !model.unique_subshapes(shape, NodeKind::Solid).is_empty() {
        return op(model, shape);
    }
    let shells = model.unique_subshapes(shape, NodeKind::Shell);
    if shells.is_empty() || !shells.iter().all(|&s| crate::topology::ShellId::try_from(s.node).is_ok_and(|id| model.is_closed_shell(id))) {
        return op(model, shape);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for shell in shells {
        let solid = model.new_empty(NodeKind::Solid);
        model.add(solid, shell)?;
        // (Its faces facing out.)
        let solid = if crate::properties::volume(model, solid).mass < 0.0 {
            let turned = model.new_empty(NodeKind::Solid);
            model.add(turned, shell.reversed())?;
            turned
        } else {
            solid
        };
        model.add(compound, solid)?;
    }
    let result = op(model, compound)?;
    let shells = model.unique_subshapes(result, NodeKind::Shell);
    match &shells[..] {
        [one] => Ok(*one),
        many => {
            let out = model.new_empty(NodeKind::Compound);
            for &s in many {
                model.add(out, s)?;
            }
            Ok(out)
        }
    }
}

/// Where two straight corners' tools meet at an end of their edges, at an angle, on a plane
/// they share, at a corner concave on that plane (a pocket's rim), and no third there: per
/// corner and end, the vertex and the normal of the mitre plane there (through the vertex,
/// along the bisector of the two edges and the plane's normal), pointing to the corner's own
/// edge. The tools then run past the vertex and are cut back to the mitre, as OCCT mitres
/// them.
fn mitres(corners: &[(Corner, f64)]) -> Vec<[Option<(Point, Vector)>; 2]> {
    // The direction from an end of the edge into it.
    let inward = |c: &Corner, at_end: bool| -> Vector {
        let p = if at_end { c.end } else { c.origin };
        c.direction_at(p) * if at_end { -1.0 } else { 1.0 }
    };
    corners
        .iter()
        .enumerate()
        .map(|(i, (corner, _))| {
            [false, true].map(|at_end| -> Option<(Point, Vector)> {
                // (Straight edges only: where an arc's blend meets another's, their surfaces
                // meet otherwise.)
                if !matches!(corner.path, Path::Line { .. }) || !corner.convex {
                    return None;
                }
                let p = if at_end { corner.end } else { corner.origin };
                let others: Vec<(usize, bool)> = corners
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .flat_map(|(j, (o, _))| {
                        let mut ends = Vec::new();
                        if !o.ends().is_empty() {
                            if o.origin.distance(p) <= tolerance::LINEAR {
                                ends.push((j, false));
                            }
                            if o.end.distance(p) <= tolerance::LINEAR {
                                ends.push((j, true));
                            }
                        }
                        ends
                    })
                    .collect();
                let [(j, other_end)] = others[..] else { return None };
                let other = &corners[j].0;
                if !other.convex || !matches!(other.path, Path::Line { .. }) {
                    return None;
                }
                let (di, dj) = (inward(corner, at_end), inward(other, other_end));
                if di.dot(dj).abs() >= 1.0 - tolerance::ANGULAR {
                    return None;
                }
                // The plane both edges bound.
                let k = corner.faces.iter().position(|f| other.faces.iter().any(|g| g.0 == f.0))?;
                let Surface::Plane(plane) = untrimmed(&corner.faces[k].1) else { return None };
                // (Only where the corner is concave on it: the other edge leaving the vertex
                // away from the side of this one the plane's face lies on. At a convex corner
                // the tools, square at their ends, meet as they should cutting each other.)
                if dj.dot(corner.into[k]) >= 0.0 {
                    return None;
                }
                let n = Vector::from(plane.normal());
                let b = di + dj;
                let m = b.cross(n);
                if !(m.length() > 1e-9) {
                    return None;
                }
                let m = m.normalized();
                Some((p, if m.dot(di) > 0.0 { m } else { -m }))
            })
        })
        .collect()
}

/// The tool's runs past its edge's ends, far enough past a mitred end for the mitre to cut
/// it back across its whole section.
fn mitre_runs(corner: &Corner, size: f64, run: [f64; 2], mitre: &[Option<(Point, Vector)>; 2]) -> [f64; 2] {
    let cap = corner.run_cap(size);
    [0, 1].map(|k| if mitre[k].is_some() { run[k].max((4.0 * size).min(cap)) } else { run[k] })
}

/// The tool cut back to its mitres: what lies on its own edge's side of each.
fn mitred_tool(model: &mut Model, tool: Shape, mitre: &[Option<(Point, Vector)>; 2]) -> Result<Shape> {
    let mut tool = tool;
    for &(p, m) in mitre.iter().flatten() {
        let b = crate::properties::bounding_box(model, tool);
        let reach = 4.0 * (b.max - b.min).length().max(1.0);
        let frame = Frame::from_z(p, Direction::from_vector(m));
        let corner = p - Vector::from(frame.x()) * reach - Vector::from(frame.y()) * reach;
        let half = crate::primitives::Cuboid::new(2.0 * reach, 2.0 * reach, reach)
            .with_frame(Frame::new(corner, frame.z(), frame.x()))
            .build(model)?;
        tool = boolean::intersection(model, tool, half)?;
    }
    Ok(tool)
}

/// For each corner (with its tool's size), how far its tool runs past the edge's start and
/// end (0 for not at all): not where another of the edges continues it tangentially (the
/// tools meet there), nor where there is material just past the end, within the strip the
/// tool takes away (a cut would carve a groove into the material carrying on there). Else
/// as far as it may, unless its path comes back into the solid (as an arc's circle may):
/// then only until its cross-section is clear of the solid, and as far again.
fn extensions(model: &Model, shape: Shape, corners: &[(Corner, f64)]) -> Vec<[f64; 2]> {
    let joints: Vec<[bool; 2]> = corners
        .iter()
        .enumerate()
        .map(|(i, (corner, _))| {
            [corner.origin, corner.end].map(|p| {
                let d = corner.direction_at(p);
                corners.iter().enumerate().any(|(j, (other, _))| {
                    j != i
                        && other.ends().iter().any(|q| q.distance(p) <= tolerance::LINEAR)
                        && d.dot(other.direction_at(p)).abs() >= 1.0 - tolerance::ANGULAR
                })
            })
        })
        .collect();
    let inside = |p: Point| matches!(boolean::classify_point(model, shape, p), Ok(boolean::PointLocation::Inside));
    let run = |corner: &Corner, size: f64, at_end: bool| -> f64 {
        let cap = corner.run_cap(size);
        if cap <= 0.0 {
            return 0.0;
        }
        let probes = |s: f64| {
            let (p, into) = corner.past(at_end, s);
            let inward = (into[0] + into[1]).normalized();
            (p + inward * (0.3 * size), [p + into[0] * (0.8 * size) + inward * (0.1 * size), p + into[1] * (0.8 * size) + inward * (0.1 * size)])
        };
        // (Material along the edge's line carrying on: no run at all.)
        if inside(probes(size).0) {
            return 0.0;
        }
        // All the way, unless the path comes back into material (as an arc's circle may):
        // then until the tool's cross-section is clear of the solid, and as far again (short
        // of where it comes back).
        let step = (0.5 * size).max(cap / 64.0);
        let mut s = size;
        let mut clear = None;
        while s <= cap {
            let (along, section) = probes(s);
            if inside(along) {
                return clear.map_or(0.0, |c: f64| (2.0 * c).min(s - size).max(c));
            }
            if clear.is_none() && !section.iter().any(|&q| inside(q)) {
                clear = Some(s);
            }
            s += step;
        }
        cap
    };
    corners
        .iter()
        .zip(&joints)
        .map(|((corner, size), joint)| [if joint[0] { 0.0 } else { run(corner, *size, false) }, if joint[1] { 0.0 } else { run(corner, *size, true) }])
        .collect()
}

/// An edge to chamfer, with how far: measured from `face` (one of the two faces at the edge;
/// either if `None`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChamferedEdge {
    pub edge: Shape,
    pub size: ChamferSize,
    pub face: Option<Shape>,
}

/// How far a chamfer cuts back along the faces at its edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChamferSize {
    /// The distance along the face, and along the other face.
    Distances(f64, f64),
    /// The distance along the face, and the angle (radians) the chamfer makes with it.
    #[doc(alias = "SetDistAngle")]
    DistanceAngle(f64, f64),
}

/// `shape` with `edges` bevelled: a flat face (a cone along circular edges) from the
/// faces' points `distance` away from each edge (measured along the faces, square to the
/// edge).
///
/// Supported where the edge is straight or circular and the faces' cross-sections there
/// are lines or circles, as for [`fillet`]; edges continuing an edge tangentially are
/// chamfered too. Two chamfered edges meeting at a vertex meet in a mitre; three meet at
/// the point where their faces do. The result is a compound of solids, like a boolean
/// operation's.
#[doc(alias = "BRepFilletAPI_MakeChamfer")]
pub fn chamfer(model: &mut Model, shape: Shape, edges: &[Shape], distance: f64) -> Result<Shape> {
    let edges: Vec<ChamferedEdge> = edges.iter().map(|&edge| ChamferedEdge { edge, size: ChamferSize::Distances(distance, distance), face: None }).collect();
    chamfer_edges(model, shape, &edges)
}

/// [`chamfer`] with two distances per edge, the first along a given face.
#[doc(alias = "BRepFilletAPI_MakeChamfer::Add")]
pub fn chamfer_edges(model: &mut Model, shape: Shape, edges: &[ChamferedEdge]) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to chamfer".into()));
    }
    let edges = &sharp(model, shape, edges, |e| e.edge)?;
    on_solids(model, shape, |model, shape| each_solid(model, shape, edges, |e| e.edge, chamfer_solid_edges))
}

/// [`chamfer_edges`] of a shape with solids.
fn chamfer_solid_edges(model: &mut Model, shape: Shape, edges: &[ChamferedEdge]) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to chamfer".into()));
    }
    for e in edges {
        let valid = match e.size {
            ChamferSize::Distances(a, b) => a > tolerance::LINEAR && b > tolerance::LINEAR,
            ChamferSize::DistanceAngle(d, angle) => d > tolerance::LINEAR && angle > 0.0 && angle < std::f64::consts::PI,
        };
        if !valid {
            return Err(Error::InvalidInput(format!("invalid chamfer size {:?}", e.size)));
        }
    }
    // The edges continuing them tangentially (indices into `edges` as their "radius").
    let indexed: Vec<(Shape, f64)> = edges.iter().enumerate().map(|(i, e)| (e.edge, i as f64)).collect();
    let chains = tangent_chains(model, shape, &indexed);
    // Where an edge's faces have other cross-sections than lines and circles: every edge
    // chamfered in turn, section by section along it.
    let mut elementary = true;
    for &(edge, _) in &chains {
        if Corner::of(model, shape, edge, true)?.is_none() {
            elementary = false;
            break;
        }
    }
    let firsts = chain_faces(model, shape, &chains, edges);
    if !elementary {
        return chamfer_in_turn(model, shape, &chains, edges, &firsts);
    }
    let mut corners: Vec<(Corner, [f64; 2])> = Vec::new();
    for (k, (edge, index)) in chains.into_iter().enumerate() {
        let wanted = edges[index as usize];
        let Some(corner) = Corner::of(model, shape, edge, true)? else {
            return Err(Error::Failed("a chamfered edge's faces changed kind".into()));
        };
        // The first distance along the chosen face (for an edge carrying the chosen one on,
        // the face on the same side).
        let chosen = firsts[k];
        let first_on_second = chosen.is_some_and(|f| corner.faces[1].0 == f);
        let (d1, d2) = match wanted.size {
            ChamferSize::Distances(a, b) => (a, b),
            ChamferSize::DistanceAngle(d, angle) => {
                // In the cross-section's triangle (the faces' tangents there): the angle at
                // the edge between the faces, `angle` at the point along the face, the sines'
                // law gives the distance along the other.
                let between = corner.into[0].dot(corner.into[1]).clamp(-1.0, 1.0).acos();
                let opposite = (between + angle).sin();
                if !(opposite > tolerance::ANGULAR) {
                    return Err(Error::InvalidInput(format!("a chamfer at {angle} rad does not meet the other face")));
                }
                (d, d * angle.sin() / opposite)
            }
        };
        corners.push((corner, if first_on_second { [d2, d1] } else { [d1, d2] }));
    }
    let sized: Vec<(Corner, f64)> = corners.iter().map(|(c, d)| (c.clone(), d[0].max(d[1]))).collect();
    let mut result = shape;
    let mitred = mitres(&sized);
    for (((corner, distances), extend), mitre) in corners.iter().zip(extensions(model, shape, &sized)).zip(&mitred) {
        let section = corner.bevel_section(model, *distances)?;
        let tool = corner.sweep_section(model, section, mitre_runs(corner, distances[0].max(distances[1]), extend, mitre))?;
        let tool = mitred_tool(model, tool, mitre)?;
        result = if corner.convex { boolean::difference(model, result, tool)? } else { boolean::union(model, result, tool)? };
    }
    for tool in chamfer_corner_tools(model, &corners)? {
        result = boolean::difference(model, result, tool)?;
    }
    Ok(result)
}

/// Each of `chains`' edges' face for its chamfer's first distance: the chosen face of the
/// edge it carries on (`edges`, by index), carried along the chain from vertex to vertex:
/// the same face where the next edge has it, else the face across from the other one, else
/// the face tangent to it there.
fn chain_faces(model: &Model, shape: Shape, chains: &[(Shape, f64)], edges: &[ChamferedEdge]) -> Vec<Option<NodeId>> {
    let faces_of = |edge: Shape| -> Vec<NodeId> {
        model.explore(shape, NodeKind::Face).filter(|&f| model.explore(f, NodeKind::Edge).any(|e| e.node == edge.node)).map(|f| f.node).collect()
    };
    let faces: Vec<Vec<NodeId>> = chains.iter().map(|&(e, _)| faces_of(e)).collect();
    let mut firsts: Vec<Option<NodeId>> = chains
        .iter()
        .map(|&(e, index)| {
            let chosen = edges[index as usize].face.map(|f| f.node)?;
            faces_of(e).contains(&chosen).then_some(chosen)
        })
        .collect();
    // (Never set where none was chosen.)
    let chosen_any = |k: usize| edges[chains[k].1 as usize].face.is_some();
    loop {
        let mut changed = false;
        for k in 0..chains.len() {
            if firsts[k].is_some() || !chosen_any(k) || faces[k].len() != 2 {
                continue;
            }
            let (a, b) = model.edge_vertices(chains[k].0);
            for j in 0..chains.len() {
                let Some(f) = firsts[j] else { continue };
                if chains[j].1 != chains[k].1 || faces[j].len() != 2 {
                    continue;
                }
                let (c, d) = model.edge_vertices(chains[j].0);
                let shared = [a, b].into_iter().find(|v| v.node == c.node || v.node == d.node);
                let Some(v) = shared else { continue };
                let other = if faces[j][0] == f { faces[j][1] } else { faces[j][0] };
                let face = if faces[k].contains(&f) {
                    Some(f)
                } else if faces[k].contains(&other) {
                    faces[k].iter().copied().find(|&g| g != other)
                } else {
                    // (The one whose normal there is nearer the chosen face's.)
                    let p = model.point(v);
                    let n = normal_at(model, shape, f, p);
                    faces[k].iter().copied().max_by(|&x, &y| {
                        let m = |g| n.zip(normal_at(model, shape, g, p)).map_or(0.0, |(n, q)| n.dot(q).abs());
                        m(x).total_cmp(&m(y))
                    })
                };
                if face.is_some() {
                    firsts[k] = face;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            return firsts;
        }
    }
}

/// [`chamfer_solid_edges`] edge by edge on the solid as it is rebuilt (each found again by
/// its curve), section by section along it; in the order given, else in reverse.
fn chamfer_in_turn(model: &mut Model, shape: Shape, chains: &[(Shape, f64)], edges: &[ChamferedEdge], firsts: &[Option<NodeId>]) -> Result<Shape> {
    // (On the shape unplaced, in its own units, the result placed as the shape is.)
    let placement = shape.placement;
    if chains.iter().any(|(e, _)| e.placement != placement) {
        todo!("chamfers of edges placed apart from their shape, between faces of this kind");
    }
    let unscale = 1.0 / placement.scale_factor().abs();
    let shape = Shape { placement: Transform::IDENTITY, ..shape };
    let mut wanted = Vec::with_capacity(chains.len());
    for (k, &(edge, index)) in chains.iter().enumerate() {
        let edge = Shape { placement: Transform::IDENTITY, orientation: Orientation::Forward, ..edge };
        let (curve, first, last) = model.edge_curve(edge).ok_or_else(|| Error::InvalidInput("a degenerate edge has no chamfer".into()))?;
        let e = edges[index as usize];
        let (distance, second) = match e.size {
            ChamferSize::Distances(a, b) => (a * unscale, ball::ChamferSecond::Distance(b * unscale)),
            ChamferSize::DistanceAngle(d, angle) => (d * unscale, ball::ChamferSecond::Angle(angle)),
        };
        wanted.push((curve, first, last, firsts[k], distance, second));
    }
    // All together first (meeting in mitres and corners, ending on free boundaries).
    let profiles: Vec<(Shape, ball::Profile)> = chains
        .iter()
        .zip(&wanted)
        .map(|(&(edge, _), (_, _, _, face, distance, second))| (Shape { placement: Transform::IDENTITY, ..edge }, ball::Profile::Chamfer { first: *face, distance: *distance, second: *second }))
        .collect();
    // (Each way taken only where it leaves the shape as sound as it was.)
    let before = fault_kinds(model, shape);
    let mut failure = match ball::blend_network(model, shape, &profiles).and_then(|done| sound(model, &before, done)) {
        Ok(done) => return Ok(done.moved(&placement)),
        Err(e) => Some(e),
    };
    let orders: &[bool] = if wanted.len() > 1 { &[false, true] } else { &[false] };
    'orders: for &reverse in orders {
        let mut sequence: Vec<_> = wanted.iter().collect();
        if reverse {
            sequence.reverse();
        }
        let mut result = shape;
        for (curve, first, last, face, distance, second) in sequence {
            let Some(edge) = find_edge(model, result, curve, *first, *last) else {
                failure = Some(Error::Failed("a chamfered edge was lost by another chamfer".into()));
                continue 'orders;
            };
            let profile = ball::Profile::Chamfer { first: *face, distance: *distance, second: *second };
            match ball::rolling_chamfer(model, result, edge, &profile) {
                Ok(done) => result = done,
                Err(e) => {
                    failure = Some(e);
                    continue 'orders;
                }
            }
        }
        match sound(model, &before, result) {
            Ok(result) => return Ok(result.moved(&placement)),
            Err(e) => failure = Some(e),
        }
    }
    Err(failure.expect("a failure"))
}

/// Where three convex chamfers of straight edges between planes meet at a vertex, the
/// corner OCCT closes with a triangle: through the points where each two chamfers meet
/// their common face. The tools cutting it: pyramids from the triangle to a point just
/// outside the vertex (their sides within what the chamfers took away).
fn chamfer_corner_tools(model: &mut Model, corners: &[(Corner, [f64; 2])]) -> Result<Vec<Shape>> {
    let planar = |c: &Corner| c.convex && matches!(c.path, Path::Line { .. }) && c.traces.iter().all(|t| matches!(t, Trace::Line));
    let mut vertices: Vec<(Point, Vec<usize>)> = Vec::new();
    for (i, (corner, _)) in corners.iter().enumerate() {
        for p in corner.ends() {
            match vertices.iter_mut().find(|(q, _)| q.distance(p) <= tolerance::LINEAR) {
                Some((_, at)) => at.push(i),
                None => vertices.push((p, vec![i])),
            }
        }
    }
    let mut tools = Vec::new();
    for (v, at) in vertices {
        if at.len() != 3 || !at.iter().all(|&i| planar(&corners[i].0)) {
            continue;
        }
        // The point where chamfers `a` and `b` meet their common face.
        let meeting = |a: usize, b: usize| -> Option<Point> {
            let ((ca, da), (cb, db)) = (&corners[a], &corners[b]);
            let ia = (0..2).find(|&k| (0..2).any(|l| ca.faces[k].0 == cb.faces[l].0))?;
            let ib = (0..2).find(|&l| cb.faces[l].0 == ca.faces[ia].0)?;
            // v + into_a d_a + s t_a = v + into_b d_b + t t_b, in the face's plane.
            let (pa, ta) = (v + ca.into[ia] * da[ia], ca.tangent);
            let (pb, tb) = (v + cb.into[ib] * db[ib], cb.tangent);
            let (aa, ab, bb) = (ta.dot(ta), ta.dot(tb), tb.dot(tb));
            let r = pb - pa;
            let det = aa * bb - ab * ab;
            if det.abs() <= 1e-12 {
                return None;
            }
            let s = (r.dot(ta) * bb - r.dot(tb) * ab) / det;
            Some(pa + ta * s)
        };
        let [a, b, c] = [at[0], at[1], at[2]];
        let (Some(p), Some(q), Some(r)) = (meeting(a, b), meeting(b, c), meeting(c, a)) else { continue };
        let size = corners[a].1.iter().chain(&corners[b].1).chain(&corners[c].1).fold(0.0f64, |m, &d| m.max(d));
        let inward = Point::barycenter(&[(p, 1.0), (q, 1.0), (r, 1.0)]) - v;
        if inward.length() <= tolerance::LINEAR {
            continue;
        }
        let apex = v - inward.normalized() * size;
        tools.push(tetrahedron(model, [p, q, r, apex])?);
    }
    Ok(tools)
}

/// The solid tetrahedron with corners `points`.
fn tetrahedron(model: &mut Model, points: [Point; 4]) -> Result<Shape> {
    let mut faces = Vec::with_capacity(4);
    for [a, b, c] in [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]] {
        let (a, b, c) = (points[a], points[b], points[c]);
        let sides = [make::edge_between(model, a, b)?, make::edge_between(model, b, c)?, make::edge_between(model, c, a)?];
        let wire = make::wire(model, &sides)?;
        faces.push(make::planar_face(model, wire)?);
    }
    let sewn = crate::sew::sew(model, &faces, &crate::sew::SewOptions::default())?;
    let shell = model.explore(sewn, NodeKind::Shell).next().ok_or_else(|| Error::Failed("a tetrahedron's faces did not sew".into()))?;
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    // (Outward: a positive volume.)
    if crate::properties::volume(model, solid).mass < 0.0 {
        let solid = model.new_empty(NodeKind::Solid);
        model.add(solid, shell.reversed())?;
        return Ok(solid);
    }
    Ok(solid)
}

/// `shape` with `edges` rounded one after another by rolling-ball fillets, each found again
/// in the solid as the ones before left it. With `rolled_over`, a single edge's ball
/// rolling on over the faces beyond has been tried on `shape` already (and failed).
fn fillet_in_turn(model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)], rolled_over: bool) -> Result<Shape> {
    // Each edge's curve and range, to find it by.
    let mut wanted = Vec::with_capacity(edges.len());
    for &(edge, radius) in edges {
        let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
        wanted.push((curve, first, last, radius));
    }
    // In the order given, else in reverse (the later fillet may be the one that has to roll
    // over onto the earlier; one edge has one order).
    let mut failure = None;
    let orders: &[bool] = if wanted.len() > 1 { &[false, true] } else { &[false] };
    for &order in orders {
        let mut sequence: Vec<_> = wanted.iter().collect();
        if order {
            sequence.reverse();
        }
        match in_order(model, shape, &sequence, rolled_over && wanted.len() == 1) {
            Ok(rounded) => return Ok(rounded),
            Err(e) => failure = Some(e),
        }
    }
    Err(failure.expect("a failure"))
}

/// [`fillet_in_turn`] in the given order: each edge's ball rolling on its faces, or rolling
/// on over the faces beyond them (a constant radius; unless `rolled_over`, for the only
/// edge).
fn in_order(model: &mut Model, shape: Shape, sequence: &[&(Curve, f64, f64, &dyn Fn(f64) -> f64)], rolled_over: bool) -> Result<Shape> {
    let mut shape = shape;
    for &(curve, first, last, radius) in sequence {
        let edge = find_edge(model, shape, curve, *first, *last).ok_or_else(|| Error::Failed("a filleted edge was lost by another fillet".into()))?;
        let r0 = radius(0.0);
        let constant = [0.25, 0.5, 0.75, 1.0].iter().all(|&s| (radius(s) - r0).abs() <= tolerance::LINEAR);
        shape = match ball::rolling_ball(model, shape, edge, *radius) {
            Ok(rounded) => rounded,
            Err(e) if constant && !rolled_over => runover::fillet(model, shape, edge, r0).map_err(|_| e)?,
            Err(e) => return Err(e),
        };
    }
    Ok(shape)
}

/// The edge of `shape` on `curve` within `[first, last]` (the longest such: what is left of
/// an edge of that curve and range).
fn find_edge(model: &Model, shape: Shape, curve: &Curve, first: f64, last: f64) -> Option<Shape> {
    use crate::geometry::ParametricCurve;
    let size = curve.point(first).distance(curve.point(0.5 * (first + last))) + curve.point(0.5 * (first + last)).distance(curve.point(last));
    let tol = 1e-6 * size.max(1.0);
    model
        .explore(shape, NodeKind::Edge)
        .filter_map(|e| {
            let (c, a, b) = model.edge_curve(Shape { orientation: Orientation::Forward, ..e })?;
            let mid = c.point(0.5 * (a + b));
            let on = crate::geometry::extrema::nearest_on_curve(mid, curve, first, last);
            let ends_on = [c.point(a), c.point(b)].iter().all(|p| crate::geometry::extrema::nearest_on_curve(*p, curve, first, last).distance <= tol);
            (on.distance <= tol && ends_on).then(|| (e, c.point(a).distance(c.point(b)) + c.point(a).distance(mid)))
        })
        .max_by(|x, y| x.1.total_cmp(&y.1))
        .map(|(e, _)| Shape { orientation: Orientation::Forward, ..e })
}

/// `edges` with the edges continuing them tangentially (sharing a face, their directions
/// parallel at the common vertex), with the same radius, repeatedly.
/// The unit normal of face `face` (a node of `shape`) at its point nearest `p`.
fn normal_at(model: &Model, shape: Shape, face: NodeId, p: Point) -> Option<Vector> {
    let face = model.explore(shape, NodeKind::Face).find(|f| f.node == face)?;
    let surface = model.surface(face)?.transformed(&face.placement);
    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    (n.length() > 1e-300).then(|| n.normalized())
}

fn tangent_chains(model: &Model, shape: Shape, edges: &[(Shape, f64)]) -> Vec<(Shape, f64)> {
    use crate::geometry::ParametricCurve;
    let all: Vec<Shape> = model.explore(shape, NodeKind::Edge).collect();
    // Each edge's end vertices with its direction there, and its faces.
    let ends = |edge: Shape| -> Vec<(NodeId, Vector)> {
        let Some((curve, first, last)) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }) else { return Vec::new() };
        let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
        [(a.node, first), (b.node, last)].into_iter().map(|(v, t)| (v, curve.d1(t).1.normalized())).collect()
    };
    let faces_of = |edge: Shape| -> Vec<NodeId> {
        model
            .explore(shape, NodeKind::Face)
            .filter(|&f| model.explore(f, NodeKind::Edge).any(|e| e.node == edge.node))
            .map(|f| f.node)
            .collect()
    };
    let mut out: Vec<(Shape, f64)> = edges.to_vec();
    let mut i = 0;
    while i < out.len() {
        let (edge, radius) = out[i];
        let (mine, my_faces) = (ends(edge), faces_of(edge));
        for &other in &all {
            if out.iter().any(|(e, _)| e.node == other.node) {
                continue;
            }
            let theirs = ends(other);
            let their_faces = faces_of(other);
            // (Between the same two faces, as nearly tangent as the faces' normals below: an
            // edge carried on by a fitted curve leaves its vertex very nearly, not exactly, as
            // the edge arrives.)
            let same_faces = their_faces.len() == 2 && their_faces.iter().all(|f| my_faces.contains(f));
            // (Elsewhere within about 1/2000 rad: edges of models made elsewhere meet that nearly
            // tangent, as continuous as the faces either side.)
            let tangent = if same_faces { tolerance::ANGULAR.sqrt() } else { 1e5 * tolerance::ANGULAR };
            let Some(at) = mine.iter().find_map(|(v, d)| {
                theirs.iter().any(|(w, e)| w == v && d.dot(*e).abs() >= 1.0 - tangent).then_some(*v)
            }) else {
                continue;
            };
            // Between the same faces, or one of them and one tangent there to the other (a
            // rounded outline's line and arc): the fillet carries on smoothly.
            let point = model.point(Shape::from(crate::topology::VertexId::try_from(at).expect("a vertex")));
            let carries_on = their_faces.iter().any(|f| my_faces.contains(f))
                && their_faces.iter().filter(|f| !my_faces.contains(f)).all(|f| {
                    my_faces.iter().filter(|g| !their_faces.contains(g)).any(|g| match (normal_at(model, shape, *f, point), normal_at(model, shape, *g, point)) {
                        (Some(a), Some(b)) => a.dot(b).abs() >= 1.0 - tolerance::ANGULAR.sqrt(),
                        _ => false,
                    })
                });
            if carries_on {
                out.push((other, radius));
            }
        }
        i += 1;
    }
    out
}

/// A trimmed surface's basis (a face's boundary trims it already), other surfaces as they
/// are.
fn untrimmed(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => untrimmed(t.basis()),
        s => s,
    }
}

/// True if every face at the edge's vertices is elementary (a plane, quadric or torus), so
/// that booleans meet a fillet tool there in closed form.
fn elementary_around(model: &Model, shape: Shape, edge: Shape) -> bool {
    let (a, b) = model.edge_vertices(edge);
    model.explore(shape, NodeKind::Face).all(|face| {
        let touches = model.explore(face, NodeKind::Vertex).any(|v| v.node == a.node || v.node == b.node);
        !touches
            || matches!(
                model.surface(face).map(untrimmed),
                Some(Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_))
            )
    })
}

/// Joins arcs of one circle that follow each other with the same faces' traces and radius
/// (a section circle split at seams) into one edge, swept by one revolution.
fn merge_arcs(mut corners: Vec<(Corner, f64)>) -> Vec<(Corner, f64)> {
    loop {
        let mut joined = false;
        'pairs: for i in 0..corners.len() {
            for j in 0..corners.len() {
                if i == j {
                    continue;
                }
                let ((a, ra), (b, rb)) = (&corners[i], &corners[j]);
                let (Path::Arc { axis: xa, radius: ka, angle: aa }, Path::Arc { axis: xb, radius: kb, angle: ab }) = (&a.path, &b.path) else { continue };
                let same_circle = Vector::from(xa.direction).dot(Vector::from(xb.direction)) > 1.0 - tolerance::ANGULAR
                    && coaxial(*xa, *xb)
                    && (ka - kb).abs() <= tolerance::LINEAR;
                let same_traces = a.traces.iter().zip(&b.traces).all(|(s, t)| std::mem::discriminant(s) == std::mem::discriminant(t));
                if !same_circle || !same_traces || a.convex != b.convex || (ra - rb).abs() > tolerance::LINEAR || a.end.distance(b.origin) > tolerance::LINEAR {
                    continue;
                }
                let angle = (aa + ab).min(TAU);
                let (axis, radius, end) = (*xa, *ka, b.end);
                corners[i].0.path = Path::Arc { axis, radius, angle };
                corners[i].0.end = end;
                corners.remove(j);
                joined = true;
                break 'pairs;
            }
        }
        if !joined {
            return corners;
        }
    }
}

/// `shape` with edges rounded by a ball whose radius varies along each: through the given
/// `(u, r)` pairs, `u` the position along the edge (rescaled so that the first and last pairs
/// are at its ends, by length), with a cubic spline flat at both ends (one pair gives a
/// constant radius). On a closed edge the spline is periodic and, as in OCCT, the last pair
/// takes the first's radius (it lies where the edge closes).
#[doc(alias = "BRepFilletAPI_MakeFillet::Add")]
#[doc(alias = "BRepFilletAPI_MakeFillet::SetRadius")]
pub fn fillet_evolving(model: &mut Model, shape: Shape, edges: &[(Shape, Vec<(f64, f64)>)]) -> Result<Shape> {
    on_solids(model, shape, |model, shape| fillet_solid_evolving(model, shape, edges))
}

/// [`fillet_evolving`] of a shape with solids.
fn fillet_solid_evolving(model: &mut Model, shape: Shape, edges: &[(Shape, Vec<(f64, f64)>)]) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    // A placed shape is filleted where it is defined, then placed back.
    let placement = shape.placement;
    if edges.iter().any(|(e, _)| e.placement != placement) {
        todo!("rolling-ball fillets of edges placed apart from their shape");
    }
    let unscale = 1.0 / placement.scale_factor().abs();
    let laws: Vec<RadiusLaw> = edges
        .iter()
        .map(|(edge, law)| {
            let (a, b) = model.edge_vertices(*edge);
            RadiusLaw::new(law, a.node == b.node)
        })
        .collect::<Result<_>>()?;
    let radii: Vec<Box<dyn Fn(f64) -> f64 + '_>> = laws.iter().map(|law| Box::new(move |s: f64| law.at(s) * unscale) as Box<dyn Fn(f64) -> f64>).collect();
    // Edges meeting at their ends, in groups, rounded together (as `fillet_edges` does).
    let ends_of = |model: &Model, edge: Shape| -> Vec<Point> {
        let (a, b) = model.edge_vertices(edge);
        if a.node == b.node { Vec::new() } else { vec![model.point(a), model.point(b)] }
    };
    let n = edges.len();
    let mut group: Vec<usize> = (0..n).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    for i in 0..n {
        for j in i + 1..n {
            let (pi, pj) = (ends_of(model, edges[i].0), ends_of(model, edges[j].0));
            if pi.iter().any(|p| pj.iter().any(|q| p.distance(*q) <= tolerance::LINEAR)) {
                let (a, b) = (root(&mut group, i), root(&mut group, j));
                group[a.max(b)] = a.min(b);
            }
        }
    }
    let mut result = Shape { placement: Transform::IDENTITY, ..shape };
    let roots: Vec<usize> = { let mut r: Vec<usize> = (0..n).map(|i| root(&mut group, i)).collect(); r.sort_unstable(); r.dedup(); r };
    for r in roots {
        let members: Vec<usize> = (0..n).filter(|&i| root(&mut group, i) == r).collect();
        let group_edges: Vec<(Shape, &dyn Fn(f64) -> f64)> = members.iter().map(|&i| (Shape { placement: Transform::IDENTITY, ..edges[i].0 }, radii[i].as_ref())).collect();
        let roll = |model: &mut Model, shape: Shape, edges: &[(Shape, &dyn Fn(f64) -> f64)]| match edges {
            [(edge, law)] => ball::rolling_ball(model, shape, *edge, *law),
            _ => ball::fillet_network(model, shape, edges).or_else(|_| fillet_in_turn(model, shape, edges, false)),
        };
        let built = roll(model, result, &group_edges);
        let before = fault_kinds(model, result);
        let crossing = crate::validity::ProblemKind::SelfIntersection;
        // (Where a face's boundary now crosses itself, the fillet runs over another
        // feature on its faces: rounded without it, the feature then put back; or off its
        // faces: rolling on onto the faces beyond.)
        let built = match built {
            Ok(rounded) if !before.contains(&crossing) && fault_kinds(model, rounded).contains(&crossing) => {
                over_features(model, result, rounded, &group_edges, &roll, 3).or_else(|e| match &group_edges[..] {
                    [(edge, law)] => runover::fillet_along(model, result, *edge, *law),
                    _ => Err(e),
                })
            }
            other => other,
        };
        // (Taken only where it leaves the solid as sound as it was.)
        let built = built.and_then(|r| sound(model, &before, r));
        result = match built {
            Ok(rounded) => rounded,
            // Else as tools, cut away or added (the booleans meeting them).
            Err(_) => {
                let original = result;
                let mut rounded = result;
                for &(edge, law) in &group_edges {
                    let (tool, convex) = ball::rolling_tool(model, original, edge, law)?;
                    rounded = if convex { boolean::difference(model, rounded, tool)? } else { boolean::union(model, rounded, tool)? };
                }
                rounded
            }
        };
    }
    Ok(result.moved(&placement))
}

/// `solid` with `edges` rounded by `fillet` where the fillets run over other features on
/// their faces (a boss or a pocket whose loop the ball's contacts cross, which `fillet`
/// leaves whole, as in `crossed`): the features beyond the loops crossed removed, the
/// edges rounded on what remains, and the features put back (the pockets cut again, the
/// bosses joined again and carried on down to the fillets).
fn over_features(
    model: &mut Model,
    solid: Shape,
    crossed: Shape,
    edges: &[(Shape, &dyn Fn(f64) -> f64)],
    fillet: &dyn Fn(&mut Model, Shape, &[(Shape, &dyn Fn(f64) -> f64)]) -> Result<Shape>,
    depth: usize,
) -> Result<Shape> {
    let failed = || Error::Failed("a fillet running over other features of its faces".into());
    let is_edge = |e: &Shape| edges.iter().any(|(f, _)| f.node == e.node);
    // The faces along the edges, the inner loops of theirs the fillets cross, and the faces
    // beyond those.
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let along: Vec<Shape> = faces.iter().copied().filter(|&f| model.unique_subshapes(f, NodeKind::Edge).iter().any(is_edge)).collect();
    let crossing: std::collections::HashSet<NodeId> =
        crate::validity::boundary_crossings(model, crossed).into_iter().flat_map(|(_, pair)| pair.map(|e| e.node)).collect();
    let mut loops = std::collections::HashSet::new();
    for &face in &along {
        for wire in model.explore(face, NodeKind::Wire).skip(1).collect::<Vec<_>>() {
            let wire_edges = model.unique_subshapes(wire, NodeKind::Edge);
            if !wire_edges.iter().any(is_edge) && wire_edges.iter().any(|e| crossing.contains(&e.node)) {
                loops.extend(wire_edges.into_iter().map(|e| e.node));
            }
        }
    }
    let features: Vec<Shape> = faces
        .iter()
        .copied()
        .filter(|f| !along.iter().any(|g| g.node == f.node))
        .filter(|&f| model.unique_subshapes(f, NodeKind::Edge).iter().any(|e| loops.contains(&e.node)))
        .collect();
    if features.is_empty() {
        return Err(failed());
    }
    let plain = crate::defeature::remove_faces(model, solid, &features)?.shape;
    let [plain] = model.unique_subshapes(plain, NodeKind::Solid)[..] else { return Err(failed()) };
    // The edges again on the plain solid (by their ends and middles), each radius law
    // turned round where its edge is.
    let ends = |model: &Model, e: Shape| -> Option<[Point; 3]> {
        let (a, b) = model.edge_vertices(e);
        let (curve, t0, t1) = model.edge_curve(e)?;
        Some([model.point(a), crate::geometry::ParametricCurve::point(&curve, 0.5 * (t0 + t1)), model.point(b)])
    };
    let candidates = model.unique_subshapes(plain, NodeKind::Edge);
    let mut found: Vec<(Shape, bool)> = Vec::new();
    for &(edge, _) in edges {
        let [a, m, b] = ends(model, edge).ok_or_else(failed)?;
        let reach = 10.0 * model.tolerance(edge).max(tolerance::LINEAR);
        let same = candidates.iter().find_map(|&c| {
            let [p, n, q] = ends(model, c)?;
            (n.distance(m) <= reach).then_some(())?;
            match (p.distance(a) <= reach && q.distance(b) <= reach, p.distance(b) <= reach && q.distance(a) <= reach) {
                (true, _) => Some((c, false)),
                (_, true) => Some((c, true)),
                _ => None,
            }
        });
        found.push(same.ok_or_else(failed)?);
    }
    let laws: Vec<Box<dyn Fn(f64) -> f64 + '_>> = found
        .iter()
        .zip(edges)
        .map(|(&(_, turned), &(_, law))| if turned { Box::new(move |s: f64| law(1.0 - s)) as Box<dyn Fn(f64) -> f64> } else { Box::new(law) })
        .collect();
    let plain_edges: Vec<(Shape, &dyn Fn(f64) -> f64)> = found.iter().zip(&laws).map(|(&(e, _), law)| (e, law.as_ref())).collect();
    let mut rounded = fillet(model, plain, &plain_edges)?;
    // (Running over more features: those taken away in turn.)
    if crosses_itself(model, rounded) {
        if depth <= 1 {
            return Err(failed());
        }
        rounded = over_features(model, plain, rounded, &plain_edges, fillet, depth - 1)?;
    }
    // The pockets (in the plain solid, not in `solid`) and the bosses (the other way).
    let pockets = boolean::difference(model, plain, solid)?;
    let bosses = boolean::difference(model, solid, plain)?;
    // (Only their solids: what else a difference leaves, as faces where the solids touch,
    // is nothing to cut or join.)
    let some = |model: &Model, s: Shape| {
        !model.unique_subshapes(s, NodeKind::Solid).is_empty() && crate::properties::volume(model, s).mass.abs() > tolerance::LINEAR.powi(3)
    };
    let solids_of = |model: &mut Model, s: Shape| -> Result<Shape> {
        let compound = model.new_empty(NodeKind::Compound);
        for solid in model.unique_subshapes(s, NodeKind::Solid) {
            model.add(compound, solid)?;
        }
        Ok(compound)
    };
    let (pockets, bosses) = (solids_of(model, pockets)?, solids_of(model, bosses)?);
    if some(model, pockets) {
        rounded = boolean::difference(model, rounded, pockets)?;
    }
    if some(model, bosses) {
        // (A boss standing on a face along the edges is carried on down to the fillet, as
        // its footprint on the face, pushed into the plain solid past the fillet's reach.)
        let reach = 4.0 * edges.iter().flat_map(|(_, law)| (0..=8).map(move |i| law(i as f64 / 8.0))).fold(0.0, f64::max);
        let planes: Vec<(crate::geometry::surface::Plane, Vector)> = along
            .iter()
            .filter_map(|&f| match model.surface(f)? {
                Surface::Plane(plane) => {
                    let outward = Vector::from(plane.normal()) * if f.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
                    Some((plane.clone(), outward))
                }
                _ => None,
            })
            .collect();
        for face in model.unique_subshapes(bosses, NodeKind::Face) {
            let Some(Surface::Plane(own)) = model.surface(face) else { continue };
            let Some(&(_, outward)) = planes.iter().find(|(plane, _)| {
                Vector::from(plane.normal()).cross(Vector::from(own.normal())).length() <= tolerance::ANGULAR && plane.distance(own.frame().origin()) <= tolerance::LINEAR
            }) else {
                continue;
            };
            let footprint = crate::modify::copy(model, face, false)?;
            let prism = crate::sweep::extrude(model, footprint, outward * -reach)?.shape;
            let below = boolean::intersection(model, prism, plain)?;
            let below = solids_of(model, below)?;
            if some(model, below) {
                rounded = boolean::union(model, rounded, below)?;
            }
        }
        rounded = boolean::union(model, rounded, bosses)?;
    }
    match crate::validity::check(model, rounded).is_valid() {
        true => Ok(rounded),
        false => Err(failed()),
    }
}

/// The kinds of fault `shape` has.
fn fault_kinds(model: &Model, shape: Shape) -> std::collections::HashSet<crate::validity::ProblemKind> {
    crate::validity::check(model, shape).errors().map(|p| p.kind).collect()
}

/// `result` if it has no kind of fault the shape it was made from had not (`before`, that
/// shape's): an operation leaves a shape as sound as it found it, else fails.
fn sound(model: &Model, before: &std::collections::HashSet<crate::validity::ProblemKind>, result: Shape) -> Result<Shape> {
    let new: Vec<_> = fault_kinds(model, result).into_iter().filter(|k| !before.contains(k)).collect();
    match new[..] {
        [] => Ok(result),
        _ => Err(Error::Failed(format!("the fillet left the shape with faults it had not: {new:?}"))),
    }
}

/// Whether a face of `shape` has a boundary crossing itself.
fn crosses_itself(model: &Model, shape: Shape) -> bool {
    crate::validity::check(model, shape).errors().any(|p| p.kind == crate::validity::ProblemKind::SelfIntersection)
}

/// A radius along an edge: the cubic spline through `(s, r)` pairs (`s` from 0 to 1) with
/// zero slope at both ends (OCCT's `Law_Interpol` with null end derivatives), or, along a
/// closed edge, the periodic one (OCCT's periodic `Law_Interpol`).
#[doc(alias = "Law_Interpol")]
struct RadiusLaw {
    s: Vec<f64>,
    r: Vec<f64>,
    /// Second derivatives at the pairs.
    m: Vec<f64>,
    periodic: bool,
}

impl RadiusLaw {
    fn new(pairs: &[(f64, f64)], closed: bool) -> Result<RadiusLaw> {
        let (first, last) = match pairs {
            [] => return Err(Error::InvalidInput("a radius law needs at least one value".into())),
            [(_, r)] => return Ok(RadiusLaw { s: vec![0.0, 1.0], r: vec![*r, *r], m: vec![0.0, 0.0], periodic: false }),
            [f, .., l] => (f.0, l.0),
        };
        if !(last > first) {
            return Err(Error::InvalidInput("radius law parameters must increase".into()));
        }
        let s: Vec<f64> = pairs.iter().map(|p| (p.0 - first) / (last - first)).collect();
        let mut r: Vec<f64> = pairs.iter().map(|p| p.1).collect();
        if s.windows(2).any(|w| !(w[0] < w[1])) || r.iter().any(|&r| !(r > 0.0)) {
            return Err(Error::InvalidInput("radius law parameters must increase and radii be positive".into()));
        }
        if closed {
            // The last pair is where the edge closes: OCCT gives it the first's radius
            // (ChFiDS_FilSpine::ComputeLaw), and the spline goes on smoothly round the seam.
            let n = s.len() - 1;
            r[n] = r[0];
            let h: Vec<f64> = s.windows(2).map(|w| w[1] - w[0]).collect();
            // Second derivatives m[0..n] (m[n] = m[0]): a cyclic tridiagonal system, small
            // enough to solve densely.
            let mut a = vec![vec![0.0; n + 1]; n];
            for i in 0..n {
                let prev = (i + n - 1) % n;
                let hp = h[prev];
                let rp = if i == 0 { r[n - 1] } else { r[i - 1] };
                a[i][prev] += hp / 6.0;
                a[i][i] += (hp + h[i]) / 3.0;
                a[i][(i + 1) % n] += h[i] / 6.0;
                a[i][n] = (r[i + 1] - r[i]) / h[i] - (r[i] - rp) / hp;
            }
            // Gaussian elimination with partial pivoting.
            for c in 0..n {
                let p = (c..n).max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs())).expect("a row");
                a.swap(c, p);
                for row in c + 1..n {
                    let f = a[row][c] / a[c][c];
                    for k in c..=n {
                        a[row][k] -= f * a[c][k];
                    }
                }
            }
            let mut m = vec![0.0; n + 1];
            for i in (0..n).rev() {
                let tail: f64 = (i + 1..n).map(|k| a[i][k] * m[k]).sum();
                m[i] = (a[i][n] - tail) / a[i][i];
            }
            m[n] = m[0];
            return Ok(RadiusLaw { s, r, m, periodic: true });
        }
        // The clamped cubic spline (slopes 0 at both ends): a tridiagonal system for the
        // second derivatives.
        let n = s.len();
        let h: Vec<f64> = s.windows(2).map(|w| w[1] - w[0]).collect();
        let mut diag = vec![0.0; n];
        let mut upper = vec![0.0; n];
        let mut lower = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        diag[0] = h[0] / 3.0;
        upper[0] = h[0] / 6.0;
        rhs[0] = (r[1] - r[0]) / h[0];
        for i in 1..n - 1 {
            lower[i] = h[i - 1] / 6.0;
            diag[i] = (h[i - 1] + h[i]) / 3.0;
            upper[i] = h[i] / 6.0;
            rhs[i] = (r[i + 1] - r[i]) / h[i] - (r[i] - r[i - 1]) / h[i - 1];
        }
        lower[n - 1] = h[n - 2] / 6.0;
        diag[n - 1] = h[n - 2] / 3.0;
        rhs[n - 1] = -(r[n - 1] - r[n - 2]) / h[n - 2];
        // Thomas algorithm.
        for i in 1..n {
            let f = lower[i] / diag[i - 1];
            diag[i] -= f * upper[i - 1];
            rhs[i] -= f * rhs[i - 1];
        }
        let mut m = vec![0.0; n];
        m[n - 1] = rhs[n - 1] / diag[n - 1];
        for i in (0..n - 1).rev() {
            m[i] = (rhs[i] - upper[i] * m[i + 1]) / diag[i];
        }
        Ok(RadiusLaw { s, r, m, periodic: false })
    }

    fn at(&self, x: f64) -> f64 {
        let x = if self.periodic { x.rem_euclid(1.0) } else { x.clamp(0.0, 1.0) };
        let k = self.s.windows(2).position(|w| x <= w[1]).unwrap_or(self.s.len() - 2);
        let (a, b) = (self.s[k], self.s[k + 1]);
        let h = b - a;
        let (p, q) = ((b - x) / h, (x - a) / h);
        p * self.r[k] + q * self.r[k + 1] + ((p * p * p - p) * self.m[k] + (q * q * q - q) * self.m[k + 1]) * h * h / 6.0
    }
}

/// Where three convex edges filleted with the same radius meet: the fillets end on the
/// ball touching the three faces. Each fillet stops at its cross-section through the
/// ball's centre; between those planes, on the vertex's side, the solid keeps only the
/// ball. The edge tools leave the region where the fillets overlap there (outside the ball,
/// which lies within each fillet); the tool cuts what remains between the planes outside
/// the ball, up to past the vertex.
struct CornerPatch {
    vertex: Point,
    center: Point,
    radius: f64,
    /// Through the centre, each normal towards the vertex.
    planes: Vec<Vector>,
}

impl CornerPatch {
    /// `None` where the corner is not one such (concave edges, unequal radii, other than
    /// three faces): the fillets are left meeting as they cut each other.
    fn of(corners: &[(Corner, f64)], at: &[usize], vertex: Point) -> Result<Option<CornerPatch>> {
        let radius = corners[at[0]].1;
        let mut faces: Vec<(NodeId, Surface, f64)> = Vec::new();
        for &i in at {
            let (corner, r) = &corners[i];
            if !corner.convex || (r - radius).abs() > tolerance::LINEAR {
                return Ok(None);
            }
            for face in &corner.faces {
                if !faces.iter().any(|f| f.0 == face.0) {
                    faces.push(face.clone());
                }
            }
        }
        let [a, b, c] = &faces[..] else {
            return Ok(None);
        };
        if faces.iter().any(|f| f.1.parameters(vertex).is_none()) {
            return Ok(None);
        }
        // Outward distances from the three surfaces; the centre lies `radius` inside each.
        let distance = |(_, surface, sign): &(NodeId, Surface, f64), p: Point| -> f64 {
            let (u, v) = surface.parameters(p).unwrap_or_else(|| todo!("corner blends on faces of this kind"));
            let (foot, du, dv) = surface.d1(u, v);
            (p - foot).dot(du.cross(dv).normalized()) * sign
        };
        let residual = |p: Point| [distance(a, p) + radius, distance(b, p) + radius, distance(c, p) + radius];
        // Newton from the centre for three planes square to each other.
        let normal_at = |face: &(NodeId, Surface, f64)| {
            let (u, v) = face.1.parameters(vertex).expect("an elementary surface");
            let (_, du, dv) = face.1.d1(u, v);
            du.cross(dv).normalized() * face.2
        };
        let mut center = vertex - (normal_at(a) + normal_at(b) + normal_at(c)) * radius;
        let h = 1e-7 * radius.max(1.0);
        for _ in 0..50 {
            let f = residual(center);
            if f.iter().all(|x| x.abs() <= 1e-12 * radius.max(1.0)) {
                break;
            }
            let columns: Vec<[f64; 3]> = [Vector::new(h, 0.0, 0.0), Vector::new(0.0, h, 0.0), Vector::new(0.0, 0.0, h)]
                .iter()
                .map(|&d| {
                    let g = residual(center + d);
                    [(g[0] - f[0]) / h, (g[1] - f[1]) / h, (g[2] - f[2]) / h]
                })
                .collect();
            let step = solve3([columns[0], columns[1], columns[2]], f)
                .ok_or_else(|| Error::InvalidInput("no ball touches the three faces at a corner".into()))?;
            center = center - Vector::new(step[0], step[1], step[2]);
        }
        if residual(center).iter().any(|x| x.abs() > 1e-9 * radius.max(1.0)) {
            return Err(Error::InvalidInput("no ball touches the three faces at a corner".into()));
        }
        // Each fillet's cross-section through the centre: square to a straight edge, through
        // the axis of a circular one.
        let planes = at
            .iter()
            .map(|&i| {
                let n = match corners[i].0.path {
                    Path::Line { .. } => corners[i].0.tangent,
                    Path::Arc { axis, .. } => Vector::from(axis.direction).cross(center - foot(axis, center)).normalized(),
                };
                if (vertex - center).dot(n) < 0.0 { -n } else { n }
            })
            .collect();
        Ok(Some(CornerPatch { vertex, center, radius, planes }))
    }

    fn tool(&self, model: &mut Model) -> Result<Shape> {
        let reach = self.vertex.distance(self.center) * 1.05 + self.radius * 0.05;
        // The spheres' poles away from the corner.
        let away = Direction::from_vector(self.center - self.vertex);
        let frame = Frame::new(self.center, away, Direction::from_vector(perpendicular(Vector::from(away))));
        let outer = crate::primitives::Sphere::new(reach).with_frame(frame).build(model)?;
        let ball = crate::primitives::Sphere::new(self.radius).with_frame(frame).build(model)?;
        let mut tool = boolean::difference(model, outer, ball)?;
        for &n in &self.planes {
            // The half-space on the vertex's side of the plane, bounded by a square larger
            // than the sphere.
            let x = perpendicular(n);
            let y = n.cross(x).normalized();
            let s = 2.0 * reach;
            let corners = [(-s, -s), (s, -s), (s, s), (-s, s)].map(|(i, j)| self.center + x * i + y * j);
            let sides: Vec<Shape> = (0..4).map(|k| make::edge_between(model, corners[k], corners[(k + 1) % 4])).collect::<Result<_>>()?;
            let wire = make::wire(model, &sides)?;
            let face = make::planar_face(model, wire)?;
            let half = crate::primitives::half_space(model, face, self.center + n * reach)?;
            tool = boolean::intersection(model, tool, half)?;
        }
        Ok(tool)
    }
}

/// A unit vector perpendicular to `v`.
fn perpendicular(v: Vector) -> Vector {
    let helper = if v.x.abs() < 0.6 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    v.cross(helper).normalized()
}

/// The solution of the 3 x 3 system with the given columns and right-hand side, if the
/// matrix is not singular (Cramer's rule).
fn solve3(columns: [[f64; 3]; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    let det = |c: [[f64; 3]; 3]| {
        c[0][0] * (c[1][1] * c[2][2] - c[2][1] * c[1][2]) - c[1][0] * (c[0][1] * c[2][2] - c[2][1] * c[0][2])
            + c[2][0] * (c[0][1] * c[1][2] - c[1][1] * c[0][2])
    };
    let d = det(columns);
    if d.abs() <= 1e-300 {
        return None;
    }
    Some(std::array::from_fn(|k| {
        let mut c = columns;
        c[k] = rhs;
        det(c) / d
    }))
}

/// How the cross-section sweeps along the edge.
#[derive(Clone)]
enum Path {
    Line { length: f64 },
    /// About the edge circle's axis, from the cross-section at the edge's start.
    Arc { axis: Axis, radius: f64, angle: f64 },
}

/// Where a face meets the edge's cross-section plane, near the edge.
#[derive(Clone, Copy)]
enum Trace {
    Line,
    Circle { center: Point, radius: f64 },
}

/// An edge where two faces meet, seen in its cross-section at `origin` (its start).
#[derive(Clone)]
struct Corner {
    origin: Point,
    /// The edge's direction at `origin`.
    tangent: Vector,
    path: Path,
    end: Point,
    /// From the edge into each face, perpendicular to the edge.
    into: [Vector; 2],
    traces: [Trace; 2],
    /// The two faces: their nodes and surfaces, with +1 where the surface's normal points
    /// out of the solid, -1 where it points in.
    faces: [(NodeId, Surface, f64); 2],
    /// True if the material lies between the faces (on the side of the smaller angle).
    convex: bool,
}

impl Corner {
    /// `None` where the fillet is neither a cylinder nor a torus.
    fn of(model: &Model, shape: Shape, edge: Shape, unwrap_trimmed: bool) -> Result<Option<Corner>> {
        let id = EdgeId::try_from(edge.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not an edge")))?;
        let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
        let (origin, end) = (model.point(a), model.point(b));
        // (A chamfer takes a trimmed curve as its basis, whose parameters are the trimmed
        // one's; fillets leave such edges to the rolling ball, which blends them where they
        // meet others as OCCT does.)
        let mut curve = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge });
        while let (true, Some((Curve::Trimmed(t), first, last))) = (unwrap_trimmed, &curve) {
            curve = Some((t.basis().clone(), *first, *last));
        }
        let (tangent, path) = match curve {
            Some((Curve::Line(_), _, _)) => ((end - origin).normalized(), Path::Line { length: origin.distance(end) }),
            Some((Curve::Circle(c), first, last)) => {
                let axis = c.frame().z_axis();
                let tangent = Vector::from(axis.direction).cross(origin - c.center()).normalized();
                (tangent, Path::Arc { axis, radius: c.radius(), angle: (last - first).min(TAU) })
            }
            _ => return Ok(None),
        };
        // A point of the edge a fraction of its length from `origin`.
        let along_edge = |s: f64| match &path {
            Path::Line { length } => origin + tangent * (s * length),
            Path::Arc { axis, angle, .. } => origin.transformed(&Transform::rotation(*axis, s * angle)),
        };
        // Each face using the edge: its outward normal at `origin`, the edge's direction as
        // the face's boundary runs (counter-clockwise about that normal), and its trace.
        let mut sides: Vec<(Vector, Vector, Trace, (NodeId, Surface, f64))> = Vec::new();
        for face in model.explore(shape, NodeKind::Face) {
            for wire in model.children(face) {
                for coedge in model.children(wire) {
                    let c = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                    if model[c].edge != id {
                        continue;
                    }
                    let surface = model
                        .surface(face)
                        .map(|s| untrimmed(s).transformed(&face.placement))
                        .ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
                    let sign = |o: Orientation| if o == Orientation::Reversed { -1.0 } else { 1.0 };
                    // At a singular point of the surface (a pole, an apex) the normal is its
                    // limit along the edge.
                    // (A surface with no closed-form inverse: not a cylinder or torus fillet.)
                    if surface.parameters(origin).is_none() {
                        return Ok(None);
                    }
                    let normal_at = |p: Point| {
                        let (u, v) = surface.parameters(p)?;
                        let (_, du, dv) = surface.d1(u, v);
                        let n = du.cross(dv);
                        (n.length() > 1e-9 * (du.length() + dv.length()).powi(2)).then(|| n.normalized())
                    };
                    let normal = normal_at(origin)
                        .or_else(|| normal_at(along_edge(1e-6)))
                        .ok_or_else(|| Error::InvalidInput("no normal at a fillet edge".into()))?
                        * sign(face.orientation);
                    let along = tangent * sign(model[c].orientation.compose(wire.orientation));
                    let Some(face_trace) = trace(&surface, &path, origin) else { return Ok(None) };
                    sides.push((normal, along, face_trace, (face.node, surface, sign(face.orientation))));
                }
            }
        }
        let [(n1, t1, trace1, _), (n2, t2, trace2, _)] = sides.iter().map(|s| (s.0, s.1, s.2, ())).collect::<Vec<_>>()[..] else {
            return Err(Error::InvalidInput(format!("a fillet edge must bound two faces, it bounds {}", sides.len())));
        };
        // The face lies to the left of its boundary seen from outside.
        let into = [n1.cross(t1).normalized(), n2.cross(t2).normalized()];
        if into[0].cross(into[1]).length() <= tolerance::ANGULAR.sqrt() {
            return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
        }
        // The material is between the faces if each face's outward normal points away from
        // the other face.
        let convex = n1.dot(into[1]) < 0.0;
        let faces = [sides[0].3.clone(), sides[1].3.clone()];
        Ok(Some(Corner { origin, tangent, path, end, into, traces: [trace1, trace2], faces, convex }))
    }

    /// The ends of the edge (none for a closed one).
    fn ends(&self) -> Vec<Point> {
        match self.path {
            Path::Arc { angle, .. } if angle >= TAU - tolerance::ANGULAR => Vec::new(),
            _ => vec![self.origin, self.end],
        }
    }

    /// The edge's direction at one of its ends (up to its sign).
    fn direction_at(&self, p: Point) -> Vector {
        match &self.path {
            Path::Line { .. } => self.tangent,
            Path::Arc { axis, .. } => {
                let z = Vector::from(axis.direction);
                let d = p - axis.origin;
                z.cross(d - z * d.dot(z)).normalized()
            }
        }
    }

    /// The solid between the edge and the fillet, run past the edge's start and end by
    /// `run` for a cut.
    fn tool(&self, model: &mut Model, radius: f64, run: [f64; 2]) -> Result<Shape> {
        let section = self.section(model, radius)?;
        self.sweep_section(model, section, run)
    }

    /// The furthest a tool of `size` across may run past either end of the edge: a cut may
    /// run past (the faces there trim it), added material must stop at them; and an arc's
    /// not round to meet itself.
    fn run_cap(&self, size: f64) -> f64 {
        if !self.convex {
            return 0.0;
        }
        match self.path {
            Path::Line { length } => size + length,
            Path::Arc { radius, angle, .. } if angle < TAU - tolerance::ANGULAR => (size + radius * angle).min(0.5 * (TAU - angle) * radius),
            Path::Arc { .. } => 0.0,
        }
    }

    /// The point `s` along the edge's path past its start (or its end, `at_end`), and the
    /// directions into the faces there.
    fn past(&self, at_end: bool, s: f64) -> (Point, [Vector; 2]) {
        match self.path {
            Path::Line { .. } => {
                let p = if at_end { self.end + self.tangent * s } else { self.origin - self.tangent * s };
                (p, self.into)
            }
            Path::Arc { axis, radius, angle } => {
                let turn = Transform::rotation(axis, if at_end { angle + s / radius } else { -s / radius });
                (self.origin.transformed(&turn), self.into.map(|v| v.transformed(&turn)))
            }
        }
    }

    /// The solid swept by a cross-section at `origin` along the edge, run past the edge's
    /// start and end by `run`.
    fn sweep_section(&self, model: &mut Model, section: Shape, run: [f64; 2]) -> Result<Shape> {
        match self.path {
            Path::Line { length } => {
                let [before, after] = run;
                let section = section.moved(&Transform::translation(self.tangent * -before));
                Ok(sweep::extrude(model, section, self.tangent * (length + before + after))?.shape)
            }
            Path::Arc { axis, radius: arc_radius, angle } => {
                let [before, after] = run.map(|r| r / arc_radius);
                let section = section.moved(&Transform::rotation(axis, -before));
                Ok(sweep::revolve(model, section, axis, (angle + before + after).min(TAU))?.shape)
            }
        }
    }

    /// The cross-section of the tool at `origin`: bounded by the faces' traces from the edge
    /// to where the ball touches them, and the ball's arc between those points.
    fn section(&self, model: &mut Model, radius: f64) -> Result<Shape> {
        let o = self.origin;
        let (center, t1, t2) = self.ball(radius)?;
        let leg = |model: &mut Model, i: usize, from: Point, to: Point| match self.traces[i] {
            Trace::Line => make::edge_between(model, from, to),
            Trace::Circle { center, .. } => arc(model, center, from, to),
        };
        let leg1 = leg(model, 0, o, t1)?;
        let ball = arc(model, center, t1, t2)?;
        let leg2 = leg(model, 1, t2, o)?;
        let wire = make::wire(model, &[leg1, ball, leg2])?;
        make::planar_face(model, wire)
    }

    /// The ball of `radius` in the cross-section at `origin`: its centre and where it
    /// touches the two faces (their surfaces).
    fn ball(&self, radius: f64) -> Result<(Point, Point, Point)> {
        let o = self.origin;
        let [w1, w2] = self.into;
        // In-plane coordinates: x along w1, y towards w2.
        let ex = w1;
        let ey = (w2 - w1 * w1.dot(w2)).normalized();
        let to2 = |p: Point| [(p - o).dot(ex), (p - o).dot(ey)];
        let to3 = |q: [f64; 2]| o + ex * q[0] + ey * q[1];
        // Each trace offset by `radius` towards the other face.
        let offsets: Vec<Offset> = (0..2)
            .map(|i| {
                let (w, other) = if i == 0 { (w1, w2) } else { (w2, w1) };
                let m = (other - w * w.dot(other)).normalized();
                match self.traces[i] {
                    Trace::Line => Ok(Offset::Line { point: to2(o + m * radius), direction: [w.dot(ex), w.dot(ey)] }),
                    Trace::Circle { center, radius: rho } => {
                        let outward = (o - center).dot(m) > 0.0;
                        let r = if outward { rho + radius } else { rho - radius };
                        if !(r > tolerance::LINEAR) {
                            return Err(Error::InvalidInput(format!("fillet radius {radius} too large for a face of radius {rho}")));
                        }
                        Ok(Offset::Circle { center: to2(center), radius: r })
                    }
                }
            })
            .collect::<Result<_>>()?;
        // The ball's centre: where the offsets meet, nearest to where it would be between
        // the faces' tangent planes.
        let half = 0.5 * w1.dot(w2).clamp(-1.0, 1.0).acos();
        let estimate = to2(o + (w1 + w2).normalized() * (radius / half.sin()));
        let center = intersections(&offsets[0], &offsets[1])
            .into_iter()
            .min_by(|a, b| dist2(*a, estimate).total_cmp(&dist2(*b, estimate)))
            .ok_or_else(|| Error::InvalidInput("no fillet of this radius fits the edge".into()))?;
        let center = to3(center);
        let touch = |i: usize| -> Point {
            match (self.traces[i], &offsets[i]) {
                (Trace::Line, Offset::Line { direction, .. }) => {
                    // The foot of the centre on the trace (a line through `o`).
                    let d = ex * direction[0] + ey * direction[1];
                    o + d * (center - o).dot(d)
                }
                (Trace::Circle { center: c, radius: rho }, _) => c + (center - c).normalized() * rho,
                _ => unreachable!(),
            }
        };
        Ok((center, touch(0), touch(1)))
    }
}

impl Corner {
    /// The cross-section of a chamfer's tool at `origin`: bounded by the faces' traces from
    /// the edge to the points `distances` along them, and the segment between those points.
    fn bevel_section(&self, model: &mut Model, distances: [f64; 2]) -> Result<Shape> {
        let o = self.origin;
        let foot = |i: usize| -> Result<Point> {
            let (w, d) = (self.into[i], distances[i]);
            Ok(match self.traces[i] {
                Trace::Line => o + w * d,
                Trace::Circle { center, radius: rho } => {
                    // `d` along the circle from `o`, setting off along `w`.
                    if d >= std::f64::consts::PI * rho {
                        return Err(Error::InvalidInput(format!("chamfer distance {d} too large for a face of radius {rho}")));
                    }
                    let (r, t) = (o - center, w.normalized() * rho);
                    let angle = d / rho;
                    center + r * angle.cos() + t * angle.sin()
                }
            })
        };
        let (t1, t2) = (foot(0)?, foot(1)?);
        let leg = |model: &mut Model, i: usize, from: Point, to: Point| match self.traces[i] {
            Trace::Line => make::edge_between(model, from, to),
            Trace::Circle { center, .. } => arc(model, center, from, to),
        };
        let leg1 = leg(model, 0, o, t1)?;
        let bevel = make::edge_between(model, t1, t2)?;
        let leg2 = leg(model, 1, t2, o)?;
        let wire = make::wire(model, &[leg1, bevel, leg2])?;
        make::planar_face(model, wire)
    }
}

/// A face's trace in the cross-section of the edge at `origin`, if the fillet along the
/// edge is a cylinder or a torus.
fn trace(surface: &Surface, path: &Path, origin: Point) -> Option<Trace> {
    Some(match (surface, path) {
        (Surface::Plane(_), _) => Trace::Line,
        (Surface::Cylinder(c), Path::Line { .. }) => {
            let axis = c.frame().z_axis();
            Trace::Circle { center: foot(axis, origin), radius: c.radius() }
        }
        (Surface::Cylinder(_) | Surface::Cone(_), Path::Arc { axis, .. }) if coaxial(surface_axis(surface), *axis) => Trace::Line,
        (Surface::Sphere(s), Path::Arc { axis, .. }) if foot(*axis, s.frame().origin()).distance(s.frame().origin()) <= tolerance::LINEAR => {
            Trace::Circle { center: s.frame().origin(), radius: s.radius() }
        }
        (Surface::Torus(t), Path::Arc { axis, .. }) if coaxial(t.frame().z_axis(), *axis) => {
            // The tube's circle in the meridian half-plane through `origin`.
            let f = t.frame();
            let radial = (origin - foot(f.z_axis(), origin)).normalized();
            Trace::Circle { center: f.origin() + radial * t.major_radius(), radius: t.minor_radius() }
        }
        _ => return None,
    })
}

fn surface_axis(surface: &Surface) -> Axis {
    match surface {
        Surface::Cylinder(c) => c.frame().z_axis(),
        Surface::Cone(c) => c.frame().z_axis(),
        _ => unreachable!(),
    }
}

/// True if two axes lie on the same line.
fn coaxial(a: Axis, b: Axis) -> bool {
    Vector::from(a.direction).cross(Vector::from(b.direction)).length() <= tolerance::ANGULAR
        && foot(a, b.origin).distance(b.origin) <= tolerance::LINEAR
}

/// The point of the axis's line nearest to `p`.
fn foot(axis: Axis, p: Point) -> Point {
    let d = Vector::from(axis.direction);
    axis.origin + d * (p - axis.origin).dot(d)
}

/// The shorter arc of the circle about `center` from `from` to `to` (in the plane of the
/// three points).
fn arc(model: &mut Model, center: Point, from: Point, to: Point) -> Result<Shape> {
    let (a, b) = (from - center, to - center);
    let normal = Direction::from_vector(a.cross(b));
    let angle = a.cross(b).length().atan2(a.dot(b));
    let circle = Circle::new(Frame::new(center, normal, Direction::from_vector(a)), a.length());
    make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, angle)
}

/// A trace offset towards the fillet's centre, in the cross-section's 2D coordinates.
enum Offset {
    Line { point: [f64; 2], direction: [f64; 2] },
    Circle { center: [f64; 2], radius: f64 },
}

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

/// Where two offsets meet.
fn intersections(a: &Offset, b: &Offset) -> Vec<[f64; 2]> {
    match (a, b) {
        (Offset::Line { point: p, direction: d }, Offset::Line { point: q, direction: e }) => {
            // p + s d = q + t e
            let det = d[0] * -e[1] + e[0] * d[1];
            if det.abs() <= 1e-15 {
                return Vec::new();
            }
            let r = [q[0] - p[0], q[1] - p[1]];
            let s = (r[0] * -e[1] + e[0] * r[1]) / det;
            vec![[p[0] + s * d[0], p[1] + s * d[1]]]
        }
        (Offset::Line { point: p, direction: d }, Offset::Circle { center: c, radius: r })
        | (Offset::Circle { center: c, radius: r }, Offset::Line { point: p, direction: d }) => {
            // |p + s d - c|² = r², with d a unit vector.
            let f = [p[0] - c[0], p[1] - c[1]];
            let b = f[0] * d[0] + f[1] * d[1];
            let disc = b * b - (f[0] * f[0] + f[1] * f[1] - r * r);
            if disc < 0.0 {
                return Vec::new();
            }
            let root = disc.sqrt();
            [-b - root, -b + root].iter().map(|s| [p[0] + s * d[0], p[1] + s * d[1]]).collect()
        }
        (Offset::Circle { center: c1, radius: r1 }, Offset::Circle { center: c2, radius: r2 }) => {
            let d = dist2(*c1, *c2);
            if d <= 1e-15 || d > r1 + r2 || d < (r1 - r2).abs() {
                return Vec::new();
            }
            let along = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
            let h = (r1 * r1 - along * along).max(0.0).sqrt();
            let u = [(c2[0] - c1[0]) / d, (c2[1] - c1[1]) / d];
            let m = [c1[0] + along * u[0], c1[1] + along * u[1]];
            vec![[m[0] - h * u[1], m[1] + h * u[0]], [m[0] + h * u[1], m[1] - h * u[0]]]
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::primitives::Cuboid;
    use crate::{properties, validity};

    /// A 20 x 20 x 20 box with a boss (a cylinder of radius 3 and height 5) standing on its
    /// top face 3.5 from its front edge: rounding that edge (R 4) runs under the boss, which
    /// stands on the fillet, carried on down to it.
    #[test]
    fn rounding_runs_under_a_boss() {
        use crate::primitives::Cylinder;
        let mut model = Model::new();
        let block = Cuboid::new(20.0, 20.0, 20.0).build(&mut model).unwrap();
        let frame = Frame::from_z(Point::new(10.0, 3.5, 20.0), Direction::Z);
        let boss = Cylinder::new(3.0, 5.0).with_frame(frame).build(&mut model).unwrap();
        let shape = boolean::union(&mut model, block, boss).unwrap();
        let front = model.unique_subshapes(shape, NodeKind::Edge).into_iter().find(|&e| {
            let (a, b) = model.edge_vertices(e);
            let (p, q) = (model.point(a), model.point(b));
            p.y.abs() < 1e-9 && q.y.abs() < 1e-9 && (p.z - 20.0).abs() < 1e-9 && (q.z - 20.0).abs() < 1e-9
        }).unwrap();
        let rounded = fillet_evolving(&mut model, shape, &[(front, vec![(0.0, 4.0)])]).unwrap();
        let report = validity::check(&model, rounded);
        assert!(report.is_valid(), "{report:?}");
        // The box less the fillet's corner, the boss, and the boss's footprint carried down
        // between the top and the fillet (less than its area times the radius).
        let volume = properties::volume(&model, rounded).mass;
        let without = 8000.0 - 20.0 * 16.0 * (1.0 - PI / 4.0) + PI * 9.0 * 5.0;
        assert!(volume > without + 1.0 && volume < without + PI * 9.0 * 4.0, "{volume} {without}");
    }

    /// A plate (100 x 76.2, 10 thick) whose back edge bulges into an arc (25.4 high, of a
    /// circle 61.9 across that comes back into the plate near its front corners): the arc's
    /// top edge bevelled (3) or rounded (R 3) is cut only along the arc, run past its ends
    /// clear of the plate but not round the circle into it (OCCT gives these values).
    #[test]
    fn tools_along_an_arc_stop_short_of_its_circle_coming_back() {
        use crate::geometry::curve::{Circle, Curve};
        let plate = |model: &mut Model| {
            let radius = (50.0f64.powi(2) + 25.4f64.powi(2)) / (2.0 * 25.4);
            let centre = Point::new(50.0, 76.2 + 25.4 - radius, 0.0);
            let start = (76.2 - centre.y).atan2(50.0);
            let arc = make::edge_from_curve_range(model, Curve::Circle(Circle::new(Frame::from_z(centre, Direction::Z), radius)), start, PI - start).unwrap();
            let corners = [Point::new(0.0, 76.2, 0.0), Point::new(0.0, 0.0, 0.0), Point::new(100.0, 0.0, 0.0), Point::new(100.0, 76.2, 0.0)];
            let mut edges = vec![arc];
            for w in corners.windows(2) {
                edges.push(make::edge_between(model, w[0], w[1]).unwrap());
            }
            let wire = make::wire(model, &edges).unwrap();
            let face = make::planar_face(model, wire).unwrap();
            crate::sweep::extrude(model, face, Vector::new(0.0, 0.0, 10.0)).unwrap().shape
        };
        let top_arc = |model: &Model, shape: Shape| {
            model.unique_subshapes(shape, NodeKind::Edge).into_iter().find(|&e| {
                matches!(model.edge_curve(e), Some((Curve::Circle(_), _, _))) && (model.point(model.edge_vertices(e).0).z - 10.0).abs() < 1e-9
            }).unwrap()
        };
        for (round, expected) in [(false, 93449.591), (true, 93751.569)] {
            let mut model = Model::new();
            let shape = plate(&mut model);
            let edge = top_arc(&model, shape);
            let result = if round { fillet(&mut model, shape, &[edge], 3.0) } else { chamfer(&mut model, shape, &[edge], 3.0) }.unwrap();
            assert!(validity::check(&model, result).is_valid());
            assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), 7);
            let volume = properties::volume(&model, result).mass;
            assert!((volume - expected).abs() < 1e-2, "{volume} not {expected}");
        }
    }

    /// A box with a hole: the hole's vertical edges (where its cylinder meets itself) left out
    /// of the edges rounded, as OCCT leaves them; only those, an error.
    #[test]
    fn seams_are_left_out_of_fillets() {
        use crate::primitives::Cylinder;
        let mut model = Model::new();
        let block = Cuboid::new(20.0, 20.0, 10.0).build(&mut model).unwrap();
        let bore = Cylinder::new(4.0, 30.0).with_frame(Frame::from_z(Point::new(10.0, 10.0, -10.0), Direction::Z)).build(&mut model).unwrap();
        let shape = boolean::difference(&mut model, block, bore).unwrap();
        let vertical = |model: &Model, cylindrical: bool| -> Vec<Shape> {
            model.unique_subshapes(shape, NodeKind::Edge).into_iter().filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9 && (p.z - q.z).abs() > 1.0 && ((p.x - 10.0).hypot(p.y - 10.0) < 5.0) == cylindrical
            }).collect()
        };
        let (corners, seams) = (vertical(&model, false), vertical(&model, true));
        // (The bore's vertical edges: its seam, and any edges between pieces of its cylinder.)
        assert!(corners.len() == 4 && !seams.is_empty());
        let alone = fillet(&mut model, shape, &corners, 2.0).unwrap();
        let all: Vec<Shape> = corners.iter().chain(&seams).copied().collect();
        let with_seams = fillet(&mut model, shape, &all, 2.0).unwrap();
        let (a, b) = (properties::volume(&model, alone).mass, properties::volume(&model, with_seams).mass);
        assert!((a - b).abs() < 1e-9, "{a} {b}");
        assert!(fillet(&mut model, shape, &seams, 2.0).is_err());
    }

    /// Two copies of one box (one node placed twice): each copy's edges bevelled on that copy
    /// alone, though the copies share their edges' nodes.
    #[test]
    fn copies_of_one_solid_are_bevelled_each_on_its_own() {
        let mut model = Model::new();
        let block = Cuboid::new(2.0, 2.0, 5.0).build(&mut model).unwrap();
        let pair = model.new_empty(NodeKind::Compound);
        for dy in [0.0, 4.0] {
            model.add(pair, block.moved(&Transform::translation(Vector::new(0.0, dy, 0.0)))).unwrap();
        }
        let mut vertical: Vec<Shape> = Vec::new();
        for e in model.explore(pair, NodeKind::Edge) {
            let (a, b) = model.edge_vertices(e);
            if (model.point(a).z - model.point(b).z).abs() > 1.0 && !vertical.iter().any(|v| v.is_same(&e)) {
                vertical.push(e);
            }
        }
        assert_eq!(vertical.len(), 8);
        let bevelled = chamfer(&mut model, pair, &vertical, 0.5).unwrap();
        assert_eq!(model.unique_subshapes(bevelled, NodeKind::Solid).len(), 2);
        let volume = properties::volume(&model, bevelled).mass;
        let expected = 2.0 * (20.0 - 4.0 * 0.125 * 5.0);
        assert!((volume - expected).abs() < 1e-9, "{volume} not {expected}");
    }

    /// A block (20 x 20 x 6) with a pocket (10 x 8, 3 deep) in its top: the pocket's rim
    /// bevelled (1) or rounded (R 1), the blends mitred at its corners, which are concave
    /// on the top face (OCCT gives these values).
    #[test]
    fn blends_round_a_pocket_rim_meet_in_mitres() {
        for (round, expected) in [(false, 2140.6667), (true, 2151.8909)] {
            let mut model = Model::new();
            let block = Cuboid::new(20.0, 20.0, 6.0).build(&mut model).unwrap();
            let pocket = Cuboid::new(10.0, 8.0, 4.0).with_frame(Frame::from_z(Point::new(5.0, 6.0, 3.0), Direction::Z)).build(&mut model).unwrap();
            let shape = boolean::difference(&mut model, block, pocket).unwrap();
            let rim: Vec<Shape> = model.unique_subshapes(shape, NodeKind::Edge).into_iter().filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                (p.z - 6.0).abs() < 1e-9 && (q.z - 6.0).abs() < 1e-9 && p.x > 1.0 && p.x < 19.0 && p.y > 1.0 && p.y < 19.0
            }).collect();
            assert_eq!(rim.len(), 4);
            let result = if round { fillet(&mut model, shape, &rim, 1.0) } else { chamfer(&mut model, shape, &rim, 1.0) }.unwrap();
            assert!(validity::check(&model, result).is_valid());
            assert_eq!(model.unique_subshapes(result, NodeKind::Face).len(), 15);
            let volume = properties::volume(&model, result).mass;
            assert!((volume - expected).abs() < 1e-3, "{volume} not {expected}");
        }
    }

    #[test]
    fn rounding_runs_along_tangent_edges() {
        // A 40 x 20 x 10 box with its vertical edges rounded (R 5): rounding one straight top
        // edge (R 2) rounds the whole top outline, lines and arcs (OCCT gives these values).
        let mut model = Model::new();
        let block = Cuboid::new(40.0, 20.0, 10.0).build(&mut model).unwrap();
        let vertical: Vec<Shape> = model
            .explore(block, NodeKind::Edge)
            .filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                (model.point(a).z - model.point(b).z).abs() > 1.0
            })
            .collect();
        let rounded = fillet(&mut model, block, &vertical, 5.0).unwrap();
        let top_line = model
            .explore(rounded, NodeKind::Edge)
            .find(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                p.z > 9.99 && q.z > 9.99 && (p.y - q.y).abs() < 1e-9 && p.distance(q) > 1.0
            })
            .unwrap();
        let result = fillet(&mut model, rounded, &[top_line], 2.0).unwrap();
        let report = validity::check(&model, result);
        assert!(report.problems.is_empty(), "{report:?}");
        assert!((properties::volume(&model, result).mass - 7692.167397995441).abs() < 1e-6);
        assert!((properties::area(&model, result).mass - 2573.8193424316737).abs() < 1e-6);
        assert_eq!(model.explore(result, NodeKind::Face).count(), 18);
    }

    #[test]
    fn rounded_box_edge() {
        // The vertical edge of a 300 x 300 x 100 box at the origin, radius 100: the side
        // faces lose 100 x 100 each, top and bottom a corner of 100² - π 100² / 4, and a
        // quarter cylinder is added.
        let mut model = Model::new();
        let block = Cuboid::new(300.0, 300.0, 100.0).build(&mut model).unwrap();
        let edge = model.explore(block, NodeKind::Edge).next().unwrap();
        let (a, b) = model.edge_vertices(edge);
        assert_eq!((model.point(a).x, model.point(a).y, model.point(b).x, model.point(b).y), (0.0, 0.0, 0.0, 0.0));
        let rounded = fillet(&mut model, block, &[edge], 100.0).unwrap();
        let report = validity::check(&model, rounded);
        assert!(report.problems.is_empty(), "{report:?}");
        let corner = 1e4 - PI * 1e4 / 4.0;
        let area = 300_000.0 - 2.0 * 1e4 - 2.0 * corner + PI / 2.0 * 1e4;
        let volume = 300.0 * 300.0 * 100.0 - corner * 100.0;
        assert!((properties::area(&model, rounded).mass - area).abs() < 1e-6 * area);
        assert!((properties::volume(&model, rounded).mass - volume).abs() < 1e-6 * volume);
    }

    /// The edges of `shape` with both ends at height `z`.
    fn edges_at(model: &Model, shape: Shape, z: f64) -> Vec<Shape> {
        model
            .unique_subshapes(shape, NodeKind::Edge)
            .into_iter()
            .filter(|&e| {
                let (p, q) = model.edge_vertices(e);
                (model.point(p).z - z).abs() < 1e-12 && (model.point(q).z - z).abs() < 1e-12
            })
            .collect()
    }

    fn check_volume(model: &Model, shape: Shape, volume: f64) {
        let report = validity::check(model, shape);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(model, shape).mass;
        assert!((v - volume).abs() < 1e-9 * volume.max(1.0), "{v} != {volume}");
    }

    #[test]
    fn chamfered_box_edges() {
        // One edge of a unit cube: a prism of d²/2 cut away.
        let mut model = Model::new();
        let cube = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let top = edges_at(&model, cube, 1.0);
        assert_eq!(top.len(), 4);
        let d = 0.1;
        let one = chamfer(&mut model, cube, &top[..1], d).unwrap();
        check_volume(&model, one, 1.0 - d * d / 2.0);
        // The top's four edges, meeting in mitres: the top slab less a frustum is cut away.
        let four = chamfer(&mut model, cube, &top, d).unwrap();
        let frustum = d / 3.0 * (1.0 + (1.0 - 2.0 * d).powi(2) + (1.0 - 2.0 * d));
        check_volume(&model, four, 1.0 - d + frustum);
        assert_eq!(model.unique_subshapes(four, NodeKind::Face).len(), 10);
        // Asymmetric, 0.2 along the top: the top face 0.6 across.
        let top_face = model
            .unique_subshapes(cube, NodeKind::Face)
            .into_iter()
            .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| model.point(v).z == 1.0))
            .unwrap();
        let edges: Vec<ChamferedEdge> = top.iter().map(|&edge| ChamferedEdge { edge, size: ChamferSize::Distances(0.2, 0.1), face: Some(top_face) }).collect();
        let lopsided = chamfer_edges(&mut model, cube, &edges).unwrap();
        let frustum = 0.1 / 3.0 * (1.0 + 0.36 + 0.6);
        check_volume(&model, lopsided, 1.0 - 0.1 + frustum);
    }

    #[test]
    fn chamfered_edge_between_bspline_faces() {
        // A unit cube made of B-spline faces: chamfered section by section, the same prisms
        // cut away as from the cube of planes.
        let mut model = Model::new();
        let planes = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let cube = crate::convert::to_nurbs(&mut model, planes).unwrap();
        let top = edges_at(&model, cube, 1.0);
        let d = 0.1;
        let one = chamfer(&mut model, cube, &top[..1], d).unwrap();
        check_volume(&model, one, 1.0 - d * d / 2.0);
        assert_eq!(model.unique_subshapes(one, NodeKind::Face).len(), 7);
        // 0.2 along the top face, 0.1 down the side; and 0.1 along the top at 45 degrees.
        let top_face = model
            .unique_subshapes(cube, NodeKind::Face)
            .into_iter()
            .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| (model.point(v).z - 1.0).abs() < 1e-9))
            .unwrap();
        let lopsided = chamfer_edges(&mut model, cube, &[ChamferedEdge { edge: top[0], size: ChamferSize::Distances(0.2, 0.1), face: Some(top_face) }]).unwrap();
        check_volume(&model, lopsided, 1.0 - 0.2 * 0.1 / 2.0);
        let angled = chamfer_edges(&mut model, cube, &[ChamferedEdge { edge: top[0], size: ChamferSize::DistanceAngle(0.1, PI / 4.0), face: Some(top_face) }]).unwrap();
        check_volume(&model, angled, 1.0 - d * d / 2.0);
    }

    #[test]
    fn chamfer_joined_across_a_split_face() {
        // A cube of B-spline faces whose top is split in two across the edge along y = 0:
        // that edge in two halves between the front and each half of the top, chamfered as
        // one (the halves' chamfers joined where the top's halves meet).
        let mut model = Model::new();
        let planes = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let cube = crate::convert::to_nurbs(&mut model, planes).unwrap();
        let top = model
            .unique_subshapes(cube, NodeKind::Face)
            .into_iter()
            .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| (model.point(v).z - 1.0).abs() < 1e-9))
            .unwrap();
        let surface = model.surface(top).unwrap().clone();
        let (u0, u1, v0, v1) = model.face_uv_bounds(top);
        // (Across x = 0.5, whichever parameter runs along x.)
        let along_u = (surface.point(u1, v0).x - surface.point(u0, v0).x).abs() > 0.5;
        let cut = if along_u { crate::divide::Cut::iso(&surface, true, 0.5 * (u0 + u1), v0, v1) } else { crate::divide::Cut::iso(&surface, false, 0.5 * (v0 + v1), u0, u1) };
        assert_eq!(crate::divide::divide_face(&mut model, top, &[cut]).unwrap().len(), 2);
        let front: Vec<Shape> = edges_at(&model, cube, 1.0)
            .into_iter()
            .filter(|&e| {
                let (p, q) = model.edge_vertices(e);
                model.point(p).y.abs() < 1e-9 && model.point(q).y.abs() < 1e-9
            })
            .collect();
        assert_eq!(front.len(), 2);
        let d = 0.1;
        let chamfered = chamfer(&mut model, cube, &front[..1], d).unwrap();
        check_volume(&model, chamfered, 1.0 - d * d / 2.0);
    }

    #[test]
    fn chamfered_box_corner() {
        // Three edges at a corner of a unit cube: three chamfer faces, and the corner closed
        // by a triangle on x + y + z = 2d (as OCCT closes it).
        let mut model = Model::new();
        let cube = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let at_origin: Vec<Shape> = model
            .unique_subshapes(cube, NodeKind::Edge)
            .into_iter()
            .filter(|&e| {
                let (p, q) = model.edge_vertices(e);
                model.point(p).distance(Point::ORIGIN) < 1e-12 || model.point(q).distance(Point::ORIGIN) < 1e-12
            })
            .collect();
        assert_eq!(at_origin.len(), 3);
        let d = 0.1;
        let bevelled = chamfer(&mut model, cube, &at_origin, d).unwrap();
        assert_eq!(model.unique_subshapes(bevelled, NodeKind::Face).len(), 6 + 3 + 1);
        let report = validity::check(&model, bevelled);
        assert!(report.problems.is_empty(), "{report:?}");
        let corner_face = model
            .unique_subshapes(bevelled, NodeKind::Face)
            .into_iter()
            .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| {
                let p = model.point(v);
                (p.x + p.y + p.z - 2.0 * d).abs() < 1e-9
            }))
            .expect("the corner triangle");
        let area = properties::area(&model, corner_face).mass;
        // An equilateral triangle of side d√2.
        assert!((area - 3f64.sqrt() / 4.0 * 2.0 * d * d).abs() < 1e-12, "{area}");
    }

    #[test]
    fn chamfered_cylinder_rim() {
        // The top rim of a unit cylinder: a cone; the revolved triangle (Pappus) cut away.
        let mut model = Model::new();
        let cylinder = crate::primitives::Cylinder::new(1.0, 1.0).build(&mut model).unwrap();
        let rim = edges_at(&model, cylinder, 1.0);
        let d = 0.1;
        let bevelled = chamfer(&mut model, cylinder, &rim, d).unwrap();
        check_volume(&model, bevelled, PI - TAU * (d * d / 2.0) * (1.0 - d / 3.0));
        assert_eq!(model.unique_subshapes(bevelled, NodeKind::Face).len(), 4);
    }

    #[test]
    fn chamfered_concave_edge() {
        // The inner edge of an L-shaped solid: a prism of d²/2 added.
        let mut model = Model::new();
        let a = Cuboid::new(2.0, 1.0, 1.0).build(&mut model).unwrap();
        let b = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 2.0)).build(&mut model).unwrap();
        let l = boolean::union(&mut model, a, b).unwrap();
        let edge = edges_at(&model, l, 1.0)
            .into_iter()
            .find(|&e| {
                let (p, q) = model.edge_vertices(e);
                model.point(p).x == 1.0 && model.point(q).x == 1.0
            })
            .unwrap();
        let bevelled = chamfer(&mut model, l, &[edge], 0.25).unwrap();
        check_volume(&model, bevelled, 3.0 + 0.25 * 0.25 / 2.0);
    }

    #[test]
    fn concave_edge() {
        // The inner edge of an L-shaped solid gets material added.
        let mut model = Model::new();
        let a = Cuboid::new(2.0, 1.0, 1.0).build(&mut model).unwrap();
        let b = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 2.0)).build(&mut model).unwrap();
        let l = boolean::union(&mut model, a, b).unwrap();
        // The inner edge runs along y at x = 1, z = 1.
        let edge = model
            .unique_subshapes(l, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (p, q) = model.edge_vertices(e);
                let (p, q) = (model.point(p), model.point(q));
                p.x == 1.0 && q.x == 1.0 && p.z == 1.0 && q.z == 1.0
            })
            .unwrap();
        let rounded = fillet(&mut model, l, &[edge], 0.25).unwrap();
        let report = validity::check(&model, rounded);
        assert!(report.problems.is_empty(), "{report:?}");
        let added = 0.25 * 0.25 * (1.0 - PI / 4.0);
        let v = properties::volume(&model, rounded).mass;
        assert!((v - (3.0 + added)).abs() < 1e-9, "{v} != {}", 3.0 + added);
    }

    /// The volume a convex fillet of radius `r` removes along a circular edge of radius `big`
    /// between a plane and a cylinder, over `angle` (Pappus: the corner's area times the
    /// path of its centroid).
    fn torus_corner_volume(big: f64, r: f64, angle: f64) -> f64 {
        let area = r * r * (1.0 - PI / 4.0);
        let centroid = r * (5.0 / 6.0 - PI / 4.0) / (1.0 - PI / 4.0);
        area * (big - centroid) * angle
    }

    fn circular_edges_at(model: &Model, shape: Shape, z: f64) -> Vec<Shape> {
        model
            .unique_subshapes(shape, NodeKind::Edge)
            .into_iter()
            .filter(|&e| matches!(model.edge_curve(e), Some((Curve::Circle(_), _, _))))
            .filter(|&e| (model.point(model.edge_vertices(e).0).z - z).abs() < 1e-9)
            .collect()
    }

    #[test]
    fn rounded_cylinder_rims() {
        // Each rim of a cylinder (R 50, h 100), and of a quarter of one, radius 10.
        for angle in [2.0 * PI, PI / 2.0] {
            for z in [0.0, 100.0] {
                let mut model = Model::new();
                let cylinder = crate::primitives::Cylinder::new(50.0, 100.0).with_angle(angle).build(&mut model).unwrap();
                let rim = circular_edges_at(&model, cylinder, z);
                assert_eq!(rim.len(), 1);
                let rounded = fillet(&mut model, cylinder, &rim, 10.0).unwrap();
                let report = validity::check(&model, rounded);
                assert!(report.problems.is_empty(), "angle {angle}, z {z}: {report:?}");
                let volume = angle / 2.0 * 2500.0 * 100.0 - torus_corner_volume(50.0, 10.0, angle);
                let v = properties::volume(&model, rounded).mass;
                assert!((v - volume).abs() < 1e-6 * volume, "angle {angle}, z {z}: {v} != {volume}");
            }
        }
    }

    #[test]
    fn box_corner() {
        // The three edges at a corner of a cube (side 10), radius 2: each edge's fillet
        // removes a corner prism along its length less 2r, and the corner keeps the ball's
        // octant out of the cube [0, r]³.
        let mut model = Model::new();
        let block = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        let at_origin: Vec<Shape> = model
            .unique_subshapes(block, NodeKind::Edge)
            .into_iter()
            .filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                model.point(a) == Point::ORIGIN || model.point(b) == Point::ORIGIN
            })
            .collect();
        assert_eq!(at_origin.len(), 3);
        let rounded = fillet(&mut model, block, &at_origin, 2.0).unwrap();
        let report = validity::check(&model, rounded);
        assert!(report.problems.is_empty(), "{report:?}");
        let r: f64 = 2.0;
        let prism = r * r * (1.0 - PI / 4.0) * (10.0 - r);
        let corner = r * r * r - PI / 6.0 * r * r * r;
        let volume = 1000.0 - 3.0 * prism - corner;
        let v = properties::volume(&model, rounded).mass;
        assert!((v - volume).abs() < 1e-6 * volume, "{v} != {volume}");
    }

    /// Rounds `edge` of `shape` with the rolling ball and checks the result against the
    /// expected volume.
    fn check_ball(model: &mut Model, shape: Shape, edge: Shape, radius: f64, volume: f64) {
        let rounded = ball::rolling_ball(model, shape, edge, &|_| radius).unwrap();
        let report = validity::check(model, rounded);
        assert!(report.problems.is_empty(), "{report:?}");
        let v = properties::volume(model, rounded).mass;
        assert!((v - volume).abs() < 1e-6 * volume, "{v} != {volume}");
    }

    #[test]
    fn rolling_ball_on_a_box_edge() {
        // As rounded_box_edge, with the general construction.
        let mut model = Model::new();
        let block = Cuboid::new(300.0, 300.0, 100.0).build(&mut model).unwrap();
        let edge = model
            .unique_subshapes(block, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                p.x == 0.0 && p.y == 0.0 && q.x == 0.0 && q.y == 0.0
            })
            .unwrap();
        let corner = 1e4 - PI * 1e4 / 4.0;
        check_ball(&mut model, block, edge, 100.0, 300.0 * 300.0 * 100.0 - corner * 100.0);
    }

    #[test]
    fn rolling_ball_on_a_cylinder_rim() {
        // As rounded_cylinder_rims (full cylinder, top rim), with the general construction.
        let mut model = Model::new();
        let cylinder = crate::primitives::Cylinder::new(50.0, 100.0).build(&mut model).unwrap();
        let rim = circular_edges_at(&model, cylinder, 100.0);
        let volume = PI * 2500.0 * 100.0 - torus_corner_volume(50.0, 10.0, 2.0 * PI);
        check_ball(&mut model, cylinder, rim[0], 10.0, volume);
    }

    #[test]
    fn rolling_ball_on_a_concave_edge() {
        // An L-shaped prism (one face per side), its inner edge along y.
        let mut model = Model::new();
        let corners = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)];
        let points: Vec<Point> = corners.iter().map(|&(x, z)| Point::new(x, 0.0, z)).collect();
        let sides: Vec<Shape> = (0..6).map(|i| make::edge_between(&mut model, points[i], points[(i + 1) % 6]).unwrap()).collect();
        let wire = make::wire(&mut model, &sides).unwrap();
        let face = make::planar_face(&mut model, wire).unwrap();
        let l = sweep::extrude(&mut model, face, Vector::new(0.0, 1.0, 0.0)).unwrap().shape;
        assert!((properties::volume(&model, l).mass - 3.0).abs() < 1e-12);
        let edge = model
            .unique_subshapes(l, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (p, q) = model.edge_vertices(e);
                let (p, q) = (model.point(p), model.point(q));
                p.x == 1.0 && q.x == 1.0 && p.z == 1.0 && q.z == 1.0
            })
            .unwrap();
        check_ball(&mut model, l, edge, 0.25, 3.0 + 0.25 * 0.25 * (1.0 - PI / 4.0));
    }
}

#[cfg(test)]
mod radius_law_tests {
    use super::RadiusLaw;

    #[test]
    fn periodic_law_on_a_closed_edge_closes_smoothly_with_the_first_radius() {
        let law = RadiusLaw::new(&[(0.0, 1.0), (1.0, 4.0), (2.0, 3.0)], true).unwrap();
        // Through the pairs, the last one taking the first's radius (OCCT).
        assert!((law.at(0.0) - 1.0).abs() < 1e-12 && (law.at(0.5) - 4.0).abs() < 1e-12 && (law.at(1.0) - 1.0).abs() < 1e-12);
        // Smooth across the seam: the slope just before 1 is the slope just after 0.
        let h = 1e-6;
        let (before, after) = ((law.at(1.0) - law.at(1.0 - h)) / h, (law.at(h) - law.at(0.0)) / h);
        assert!((before - after).abs() < 1e-3, "slopes {before} and {after}");
        // Symmetric here (two equal spans): r(x) = r(1 - x).
        assert!((law.at(0.2) - law.at(0.8)).abs() < 1e-12);
        // A constant law stays constant.
        let flat = RadiusLaw::new(&[(0.0, 2.0)], true).unwrap();
        assert_eq!(flat.at(0.3), 2.0);
    }

    #[test]
    fn open_law_is_flat_at_both_ends() {
        let law = RadiusLaw::new(&[(0.0, 1.0), (1.0, 4.0), (2.0, 3.0)], false).unwrap();
        assert!((law.at(1.0) - 3.0).abs() < 1e-12);
        let h = 1e-6;
        assert!(((law.at(h) - law.at(0.0)) / h).abs() < 1e-3 && ((law.at(1.0) - law.at(1.0 - h)) / h).abs() < 1e-3);
    }
}
