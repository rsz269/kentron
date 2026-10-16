//! Pipes along spines with corners (where the tangent jumps), the profile turned at each
//! corner with the spine: by the least rotation taking the tangent before the corner onto
//! the tangent after it, about their common normal.
//!
//! Each edge of the spine is swept on its own (along a line an extrusion, along an arc a
//! revolution, else a pipe), the profile carried from each to the next. At a mitred corner
//! ([`Corner::Mitred`]) the two edges' sweeps are carried on past it along their curves, far
//! enough to reach the plane bisecting it, and joined: split by one another, where what is
//! left of either beyond the plane (its end carried on past the corner) is dropped. Their
//! surfaces then meet where they cross, as OCCT trims sweeps meeting at corners against each
//! other (in the plane, where the pieces are straight: the profile turned at the corner is
//! its own mirror image across it).
//!
//! A closed profile is swept as the face it bounds; a shell wanted is the boundary of the
//! solid swept, less its ends. At a round corner ([`Corner::Rounded`]) the profile is turned
//! about the axis through the corner square to both tangents, joining the pieces' sweeps.

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::ParametricCurve;
use crate::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use crate::topology::{Model, NodeKind, Shape};
use crate::{boolean, make, properties, tolerance};

use super::pipe::SpinePiece;
use super::Corner;

/// The profile swept along `pieces`, going round their corners as `corner` says.
pub(super) fn pipe_cornered(model: &mut Model, profile: Shape, pieces: &[SpinePiece], corner: Corner) -> Result<Shape> {
    // Edges of the spine too short to matter (a millionth of it, where files leave a step
    // back and forth) are left out: they would make corners turning right back.
    let length = |p: &SpinePiece| {
        let middle = p.curve.point(0.5 * (p.first + p.last));
        p.curve.point(p.first).distance(middle) + middle.distance(p.curve.point(p.last))
    };
    let total: f64 = pieces.iter().map(length).sum();
    let pieces: Vec<SpinePiece> = pieces
        .iter()
        .filter(|p| length(p) > 1e-6 * total)
        .map(|p| SpinePiece { curve: p.curve.clone(), first: p.first, last: p.last })
        .collect();
    let pieces = &pieces[..];
    let n = pieces.len();
    // Each piece's start and end, with the tangent there.
    let ends: Vec<(Point, Vector, Point, Vector)> = pieces
        .iter()
        .map(|p| {
            let (a, ta) = p.point_tangent(0.0)?;
            let (b, tb) = p.point_tangent(1.0)?;
            Ok((a, ta, b, tb))
        })
        .collect::<Result<_>>()?;
    // A spine closed on itself goes round a corner (or on smoothly) where it closes too.
    let closed = n > 1 && ends[n - 1].2.distance(ends[0].0) <= tolerance::LINEAR;
    let wanted_shell = profile.kind() == NodeKind::Wire;
    let face = match profile.kind() {
        NodeKind::Face => profile,
        // (An open wire is swept as it is round round corners.)
        NodeKind::Wire => match make::planar_face(model, profile) {
            Ok(face) => face,
            Err(_) if corner == Corner::Rounded => profile,
            Err(_) => return Err(Error::InvalidInput("a pipe round corners needs a closed planar profile".into())),
        },
        kind => return Err(Error::InvalidInput(format!("a pipe round corners sweeps a wire or a face, not a {kind:?}"))),
    };
    // Where the profile stands at the start of each piece.
    let (k0, s0) = super::pipe::section_placement(model, face, pieces)?;
    let mut starts = vec![Transform::IDENTITY; n];
    starts[k0] = carried(&pieces[k0], s0, 0.0)?;
    for k in k0 + 1..n {
        starts[k] = turn(ends[k - 1].2, ends[k - 1].3, ends[k].1) * carried(&pieces[k - 1], 0.0, 1.0)? * starts[k - 1];
    }
    for k in (0..k0).rev() {
        let end = turn(ends[k].2, ends[k].3, ends[k + 1].1).inverse() * starts[k + 1];
        starts[k] = carried(&pieces[k], 1.0, 0.0)? * end;
    }
    // The corners: after piece k (and before the next).
    let corners: Vec<usize> = (0..n)
        .filter(|&k| k + 1 < n || closed)
        .filter(|&k| ends[k].3.dot(ends[(k + 1) % n].1) < 1.0 - 1e-12)
        .collect();
    if corner == Corner::Rounded {
        return rounded(model, face, pieces, &ends, &starts, &corners, closed, wanted_shell);
    }
    for &k in &corners {
        if ends[k].3.dot(ends[(k + 1) % n].1) <= (std::f64::consts::PI - 0.025).cos() {
            return Err(Error::InvalidInput("a pipe's spine turns back on itself at a corner".into()));
        }
    }
    // How far the profile reaches from the spine (where it stands).
    let reach = {
        let b = properties::bounding_box(model, face);
        let p = pieces[k0].curve.point(pieces[k0].first + (pieces[k0].last - pieces[k0].first) * s0);
        let corners = [b.min.x, b.max.x].into_iter().flat_map(|x| [b.min.y, b.max.y].into_iter().flat_map(move |y| [b.min.z, b.max.z].into_iter().map(move |z| Point::new(x, y, z))));
        corners.map(|q| q.distance(p)).fold(0.0, f64::max)
    };
    let half_angle = |k: usize| 0.5 * ends[k].3.dot(ends[(k + 1) % n].1).clamp(-1.0, 1.0).acos();
    let bisector = |k: usize| {
        let normal = ends[k].3 + ends[(k + 1) % n].1;
        (ends[k].2, normal * (1.0 / normal.length()))
    };
    // How far past corner `c` the piece before it reaches the bisecting plane (and before
    // it, back along its tangent, the piece after: the profile oblique to the spine, they
    // are not each other's mirror images): the farthest of the profile's points there, along
    // the tangent, and half as far again (the piece may curve away).
    let profile_points: Vec<Point> = {
        let mut points = Vec::new();
        for e in model.unique_subshapes(face, NodeKind::Edge) {
            if let Some((curve, a, b)) = model.edge_curve(e) {
                points.extend((0..=16).map(|i| curve.point(a + (b - a) * i as f64 / 16.0)));
            }
        }
        points
    };
    let extension = |c: usize| -> Result<f64> {
        let at_end = carried(&pieces[c], 0.0, 1.0)? * starts[c];
        let (p, normal) = bisector(c);
        let cos = ends[c].3.dot(normal);
        let farthest = profile_points.iter().map(|&q| -(q.transformed(&at_end) - p).dot(normal) / cos).fold(0.0, f64::max);
        Ok(1.5 * farthest + 1e-3 * reach)
    };
    let extension_back = |c: usize| -> Result<f64> {
        let next = (c + 1) % n;
        let (p, normal) = bisector(c);
        let cos = ends[next].1.dot(normal);
        let farthest = profile_points.iter().map(|&q| (q.transformed(&starts[next]) - p).dot(normal) / cos).fold(0.0, f64::max);
        Ok(1.5 * farthest + 1e-3 * reach)
    };
    // Each piece swept, carried on past its corners and cut there by the plane bisecting
    // the corner: the profile turned at the corner, its section there is the next piece's.
    // The pieces' faces are sewn (not fused: faces of two pieces lying in one plane stay
    // two faces, meeting at the corner), less the sections where the pieces meet.
    let spine_length: f64 = pieces.iter().map(length).sum();
    let plane_size = 4.0 * (reach + spine_length);
    let cut = |model: &mut Model, part: Shape, (p, normal): (Point, Vector), keep_before: bool| -> Result<Shape> {
        let plane = crate::geometry::surface::Plane::new(p, Direction::try_from_vector(normal).ok_or_else(|| Error::Failed("a corner of no direction".into()))?);
        let (u, v) = plane.parameters(p);
        let bound = make::face_on_domain(model, crate::geometry::surface::Surface::Plane(plane), u - plane_size, u + plane_size, v - plane_size, v + plane_size)?;
        // (The side dropped: past the plane for the piece before the corner.)
        let away = if keep_before { normal } else { -normal };
        let half = crate::primitives::half_space(model, bound, p + away * (1e-3 * plane_size))?;
        boolean::difference(model, part, half)
    };
    // Where consecutive pieces meet (after piece k): the plane bisecting a corner, else
    // the section of the profile carried there (the same for both).
    let (profile_centre, profile_normal) = super::pipe::profile_axis(model, face);
    let section_plane = |k: usize| -> Result<Option<(Point, Vector)>> {
        if corners.contains(&k) {
            return Ok(Some(bisector(k)));
        }
        let at_end = carried(&pieces[k], 0.0, 1.0)? * starts[k];
        Ok(profile_normal.map(|n| (profile_centre.transformed(&at_end), n.transformed(&at_end))))
    };
    let lies_in = |model: &Model, f: Shape, (p, normal): (Point, Vector)| {
        let unit = normal * (1.0 / normal.length());
        matches!(model.surface(f), Some(crate::geometry::surface::Surface::Plane(_)))
            && model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| (model.point(v) - p).dot(unit).abs() <= 1e3 * tolerance::LINEAR)
    };
    let mut faces: Vec<Shape> = Vec::new();
    // (Each piece also as swept, carried on past its corners, uncut: for the other way of
    // joining them, below.)
    let mut extended: Vec<(Shape, (f64, f64))> = Vec::new();
    for (k, piece) in pieces.iter().enumerate() {
        let previous = (k + n - 1) % n;
        let joined_before = k > 0 || closed;
        let joined_after = k + 1 < n || closed;
        let before = if joined_before && corners.contains(&previous) { extension_back(previous)? } else { 0.0 };
        let after = if corners.contains(&k) { extension(k)? } else { 0.0 };
        let spine = carried_on_edge(model, piece, before, after)?;
        let start = if before > 0.0 { carried_back(piece, before)? * starts[k] } else { starts[k] };
        let swept = super::pipe(model, face.moved(&start), spine)?.shape;
        extended.push((swept, (before, after)));
        let mut part = swept;
        if before > 0.0 {
            part = cut(model, part, bisector(previous), false)?;
        }
        if after > 0.0 {
            part = cut(model, part, bisector(k), true)?;
        }
        // Its faces, less those where it meets the pieces either side.
        let mut dropped: Vec<(Point, Vector)> = Vec::new();
        if joined_before {
            dropped.extend(section_plane(previous)?);
        }
        if joined_after {
            dropped.extend(section_plane(k)?);
        }
        for f in model.unique_subshapes(part, NodeKind::Face) {
            if !dropped.iter().any(|&plane| lies_in(model, f, plane)) {
                faces.push(f);
            }
        }
    }
    // (Within the cuts' edges' tolerances: on curved pieces the two cuts' curves of a joint,
    // each fitted on its own, agree only that closely.)
    let gap = faces
        .iter()
        .flat_map(|&f| model.unique_subshapes(f, NodeKind::Edge))
        .map(|e| model.tolerance(e))
        .fold(0.0, f64::max);
    let options = crate::sew::SewOptions { tolerance: (2.0 * gap).clamp(1e-6, 1e-3 * reach.max(1e-9)), ..Default::default() };
    let sewn = crate::sew::sew(model, &faces, &options)?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let closed_up = |model: &Model, shell: Shape| {
        let mut uses: std::collections::HashMap<crate::topology::NodeId, usize> = std::collections::HashMap::new();
        for e in model.explore(shell, NodeKind::Edge) {
            *uses.entry(e.node).or_default() += 1;
        }
        uses.values().all(|&u| u == 2)
    };
    let solid = match shells[..] {
        [shell] if closed_up(model, shell) => {
            let solid = model.new_empty(NodeKind::Solid);
            model.add(solid, shell)?;
            if properties::volume(model, solid).mass < 0.0 {
                let turned = model.new_empty(NodeKind::Solid);
                model.add(turned, shell.reversed())?;
                turned
            } else {
                solid
            }
        }
        // Where the pieces' sections at a corner differ (curved pieces, each cut by the
        // plane on its own): the pieces carried on past the corners fused, what is left of
        // either beyond a plane dropped.
        _ => {
            let joint = |k: usize, toward_next: bool| {
                let (p, normal) = bisector(k);
                let length = extended[k].1 .1.max(extended[(k + 1) % n].1 .0);
                (p, if toward_next { normal } else { -normal }, 2.0 * (reach + length))
            };
            let mut result = extended[0].0;
            for k in 1..n {
                let mut planes = vec![joint(k - 1, true)];
                if closed && k == n - 1 {
                    planes.push(joint(k, false));
                }
                result = joined(model, result, extended[k].0, &planes)?;
            }
            let options = crate::unify::UnifyOptions { same_surface_only: true, ..Default::default() };
            let result = crate::unify::unify_same_domain_with(model, result, options)?;
            let solids = model.unique_subshapes(result, NodeKind::Solid);
            let [solid] = solids[..] else {
                return Err(Error::Failed("a pipe round corners whose pieces do not join".into()));
            };
            solid
        }
    };
    if !wanted_shell {
        return Ok(solid);
    }
    // The boundary less the ends (where the profile stands at the spine's ends).
    let caps: Vec<(Point, Vector)> = match (closed, profile_normal) {
        (false, Some(normal)) => {
            let at_end = carried(&pieces[n - 1], 0.0, 1.0)? * starts[n - 1];
            vec![(profile_centre.transformed(&starts[0]), normal.transformed(&starts[0])), (profile_centre.transformed(&at_end), normal.transformed(&at_end))]
        }
        _ => Vec::new(),
    };
    let shell = model.new_empty(NodeKind::Shell);
    for f in model.unique_subshapes(solid, NodeKind::Face) {
        if !caps.iter().any(|&plane| lies_in(model, f, plane)) {
            model.add(shell, f)?;
        }
    }
    Ok(shell)
}

/// The pipe round its corners (after the pieces `corners`): each piece swept on its own,
/// and at each corner the profile turned about the axis through it square to both
/// tangents, as OCCT's round corners are; their faces sewn, with the ends closed for a
/// solid.
#[allow(clippy::too_many_arguments)]
fn rounded(
    model: &mut Model,
    face: Shape,
    pieces: &[SpinePiece],
    ends: &[(Point, Vector, Point, Vector)],
    starts: &[Transform],
    corners: &[usize],
    closed: bool,
    wanted_shell: bool,
) -> Result<Shape> {
    let n = pieces.len();
    // A face swept: the union of the pieces' solids (overlapping inside each corner) and
    // the knees (filling it outside), less the ends for a shell.
    if face.kind() == NodeKind::Face {
        let mut solid: Option<Shape> = None;
        for (k, piece) in pieces.iter().enumerate() {
            let spine = carried_on_edge(model, piece, 0.0, 0.0)?;
            let mut parts = vec![super::pipe(model, face.moved(&starts[k]), spine)?.shape];
            if corners.contains(&k) {
                let (p, t1, t2) = (ends[k].2, ends[k].3, ends[(k + 1) % n].1);
                let axis = t1.cross(t2);
                let angle = axis.length().atan2(t1.dot(t2));
                let at_end = carried(piece, 0.0, 1.0)? * starts[k];
                // (The part of the section outside the corner, past the axis: the rest turns
                // inside the pieces, which overlap there.)
                let outside = (t1 * t1.dot(t2) - t2) * (1.0 / (t1 * t1.dot(t2) - t2).length().max(1e-300));
                let section = face.moved(&at_end);
                let b = properties::bounding_box(model, section);
                let size = (b.max - b.min).length().max(1e-9);
                let plane = crate::geometry::surface::Plane::new(p, Direction::from_vector(outside));
                let (u, v) = plane.parameters(p);
                let bound = make::face_on_domain(model, crate::geometry::surface::Surface::Plane(plane), u - 4.0 * size, u + 4.0 * size, v - 4.0 * size, v + 4.0 * size)?;
                let inside = crate::primitives::half_space(model, bound, p - outside * (1e-3 * size))?;
                let outer = boolean::difference(model, section, inside)?;
                for f in model.unique_subshapes(outer, NodeKind::Face) {
                    parts.push(super::revolve(model, f, Axis::new(p, Direction::from_vector(axis)), angle)?.shape);
                }
            }
            for part in parts {
                solid = Some(match solid {
                    None => part,
                    Some(s) => boolean::union(model, s, part)?,
                });
            }
        }
        let solids = model.unique_subshapes(solid.expect("a piece"), NodeKind::Solid);
        let [solid] = solids[..] else {
            return Err(Error::Failed("a pipe round corners whose pieces do not join".into()));
        };
        if !wanted_shell {
            return Ok(solid);
        }
        let (centre, normal) = super::pipe::profile_axis(model, face);
        let caps: Vec<(Point, Vector)> = match (closed, normal) {
            (false, Some(normal)) => {
                let at_end = carried(&pieces[n - 1], 0.0, 1.0)? * starts[n - 1];
                vec![(centre.transformed(&starts[0]), normal.transformed(&starts[0])), (centre.transformed(&at_end), normal.transformed(&at_end))]
            }
            _ => Vec::new(),
        };
        let in_cap = |model: &Model, f: Shape| {
            matches!(model.surface(f), Some(crate::geometry::surface::Surface::Plane(_)))
                && caps.iter().any(|&(p, t)| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| (model.point(v) - p).dot(t * (1.0 / t.length())).abs() <= 1e3 * tolerance::LINEAR))
        };
        let shell = model.new_empty(NodeKind::Shell);
        for f in model.unique_subshapes(solid, NodeKind::Face) {
            if !in_cap(model, f) {
                model.add(shell, f)?;
            }
        }
        return Ok(shell);
    }
    // (`face` an open wire swept as it is.)
    let wires: Vec<Shape> = if face.kind() == NodeKind::Wire { vec![face] } else { model.children(face).collect() };
    let [wire] = wires[..] else {
        return Err(Error::InvalidInput("round corners of pipes with profiles with holes are not supported yet".into()));
    };
    let mut faces: Vec<Shape> = Vec::new();
    for (k, piece) in pieces.iter().enumerate() {
        let spine = carried_on_edge(model, piece, 0.0, 0.0)?;
        let swept = super::pipe(model, wire.moved(&starts[k]), spine)?.shape;
        faces.extend(model.unique_subshapes(swept, NodeKind::Face));
        if corners.contains(&k) {
            let (p, t1, t2) = (ends[k].2, ends[k].3, ends[(k + 1) % n].1);
            let axis = t1.cross(t2);
            let angle = axis.length().atan2(t1.dot(t2));
            let at_end = carried(piece, 0.0, 1.0)? * starts[k];
            let knee = super::revolve(model, wire.moved(&at_end), Axis::new(p, Direction::from_vector(axis)), angle)?.shape;
            faces.extend(model.unique_subshapes(knee, NodeKind::Face));
        }
    }
    if !closed && !wanted_shell && face.kind() == NodeKind::Face {
        let at_end = carried(&pieces[n - 1], 0.0, 1.0)? * starts[n - 1];
        faces.push(face.moved(&starts[0]));
        faces.push(face.moved(&at_end));
    }
    let faces: Vec<Shape> = faces.into_iter().map(|f| if f.placement == Transform::IDENTITY { f } else { crate::modify::copy(model, f, false).unwrap_or(f) }).collect();
    let sewn = crate::sew::sew(model, &faces, &Default::default())?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let [shell] = shells[..] else {
        return Err(Error::Failed("a pipe round corners whose pieces do not sew".into()));
    };
    if wanted_shell || face.kind() == NodeKind::Wire {
        return Ok(shell);
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    if properties::volume(model, solid).mass < 0.0 {
        let turned = model.new_empty(NodeKind::Solid);
        model.add(turned, shell.reversed())?;
        return Ok(turned);
    }
    Ok(solid)
}

/// The least rotation at a corner at `p` taking the tangent `t1` onto `t2`.
fn turn(p: Point, t1: Vector, t2: Vector) -> Transform {
    let axis = t1.cross(t2);
    if axis.length() <= 1e-12 {
        return Transform::IDENTITY;
    }
    Transform::rotation(Axis::new(p, Direction::from_vector(axis)), axis.length().atan2(t1.dot(t2)))
}

/// The motion carrying the profile along `piece` from fraction `from` to `to` (0 or 1),
/// turning with the tangent and no more.
fn carried(piece: &SpinePiece, from: f64, to: f64) -> Result<Transform> {
    let at = |s: f64| piece.first + (piece.last - piece.first) * s;
    match &piece.curve {
        Curve::Line(_) => Ok(Transform::translation(piece.curve.point(at(to)) - piece.curve.point(at(from)))),
        Curve::Circle(c) => Ok(Transform::rotation(c.frame().z_axis(), at(to) - at(from))),
        _ => {
            let single = [SpinePiece { curve: piece.curve.clone(), first: piece.first, last: piece.last }];
            let (motions, _) = super::pipe::motions(&single, 16, (0, from))?;
            let transforms = &motions[0].transforms;
            Ok(if to <= 0.0 { transforms[0] } else { *transforms.last().expect("samples") })
        }
    }
}

/// How fast `piece`'s curve runs (a line or an arc, the curves pieces are carried on along).
fn speed(piece: &SpinePiece) -> Result<f64> {
    match &piece.curve {
        Curve::Line(_) => Ok(1.0),
        Curve::Circle(c) => Ok(c.radius()),
        _ => Err(Error::Failed("mitred corners of pipes along lines and arcs only".into())),
    }
}

/// The motion carrying the profile from the start of `piece` back along its curve carried on
/// before it, by `length`.
fn carried_back(piece: &SpinePiece, length: f64) -> Result<Transform> {
    let sign = if piece.last >= piece.first { 1.0 } else { -1.0 };
    match &piece.curve {
        Curve::Line(l) => Ok(Transform::translation(Vector::from(l.axis().direction) * (-sign * length))),
        Curve::Circle(c) => Ok(Transform::rotation(c.frame().z_axis(), -sign * length / c.radius())),
        _ => Err(Error::Failed("mitred corners of pipes along lines and arcs only".into())),
    }
}

/// An edge along `piece` carried on along its curve by `before` before its start and
/// `after` past its end, run the way the piece runs.
fn carried_on_edge(model: &mut Model, piece: &SpinePiece, before: f64, after: f64) -> Result<Shape> {
    let speed = if before == 0.0 && after == 0.0 { 1.0 } else { speed(piece)? };
    let forward = piece.last >= piece.first;
    let (lo, hi) = if forward { (piece.first - before / speed, piece.last + after / speed) } else { (piece.last - after / speed, piece.first + before / speed) };
    // (An arc's range shifted by whole turns to start in [0, 2π).)
    let (lo, hi) = match &piece.curve {
        Curve::Circle(_) => {
            let shift = lo.div_euclid(std::f64::consts::TAU) * std::f64::consts::TAU;
            (lo - shift, hi - shift)
        }
        _ => (lo, hi),
    };
    let edge = make::edge_from_curve_range(model, piece.curve.clone(), lo, hi)?;
    Ok(if forward { edge } else { edge.reversed() })
}

/// The union of `a` and `b` (overlapping where they are carried on past corners between
/// them) joined across `planes` (each a point, the normal pointing from `a` to `b`, and how
/// far from the point it divides them): split by one another, what is left of either beyond
/// a plane, near it, is that one carried on past the corner and is dropped. Their surfaces
/// then meet where they cross (as OCCT trims sweeps meeting at corners against each other).
fn joined(model: &mut Model, a: Shape, b: Shape, planes: &[(Point, Vector, f64)]) -> Result<Shape> {
    // A piece's end beyond a corner: wholly near it, and beyond the plane.
    let beyond = |model: &Model, solid: Shape, from_a: bool| {
        let centre = properties::volume(model, solid).center_of_mass;
        let b = properties::bounding_box(model, solid);
        planes.iter().any(|&(o, normal, reach)| {
            let side = (centre - o).dot(normal);
            let near = [b.min.x, b.max.x].iter().all(|&x| [b.min.y, b.max.y].iter().all(|&y| [b.min.z, b.max.z].iter().all(|&z| Point::new(x, y, z).distance(o) <= reach)));
            near && if from_a { side > 0.0 } else { side < 0.0 }
        })
    };
    // Split together, the pieces share the faces where they meet: those of two pieces kept
    // are inside the result, the others bound it.
    let [both, only_a, only_b] = boolean::partition(model, a, b)?;
    let mut kept: Vec<Shape> = model.unique_subshapes(both, NodeKind::Solid);
    for (from_a, pieces) in [(true, only_a), (false, only_b)] {
        for solid in model.unique_subshapes(pieces, NodeKind::Solid) {
            if !beyond(model, solid, from_a) {
                kept.push(solid);
            }
        }
    }
    let mut uses: std::collections::HashMap<crate::topology::NodeId, (Shape, usize)> = std::collections::HashMap::new();
    let mut order = Vec::new();
    for &solid in &kept {
        for face in model.explore(solid, NodeKind::Face).collect::<Vec<_>>() {
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
    Ok(solid)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Circle;
    use crate::math::Frame;
    use crate::sweep::{pipe_cornered, PipeOptions};

    /// A disc of radius `r` at `centre`, square to `normal` (its boundary wire, or the face).
    fn disc(model: &mut Model, centre: Point, normal: Direction, r: f64, face: bool) -> Shape {
        let circle = Circle::new(Frame::new(centre, normal, Direction::Z), r);
        let edge = make::edge_from_curve(model, Curve::Circle(circle)).unwrap();
        let wire = make::wire(model, &[edge]).unwrap();
        if face { make::planar_face(model, wire).unwrap() } else { wire }
    }

    /// The wire through `points` (closed if the last is the first).
    fn polyline(model: &mut Model, points: &[Point]) -> Shape {
        let vertices: Vec<Shape> = points.iter().map(|&p| make::vertex(model, p)).collect();
        let edges: Vec<Shape> = vertices.windows(2).map(|w| make::edge(model, w[0], w[1]).unwrap()).collect();
        make::wire(model, &edges).unwrap()
    }

    #[test]
    fn a_tube_mitred_round_a_corner() {
        // A tube of radius 1 along an L of legs 10: mitred, it is as big as its spine is long
        // (what the joint gains outside the corner it loses inside), two cylinders and two
        // ends; swept from the circle, the cylinders alone.
        let mut model = Model::new();
        let spine = polyline(&mut model, &[Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0), Point::new(10.0, 10.0, 0.0)]);
        let face = disc(&mut model, Point::ORIGIN, Direction::X, 1.0, true);
        let solid = pipe_cornered(&mut model, face, spine, &PipeOptions::default(), Corner::Mitred).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid());
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - std::f64::consts::PI * 20.0).abs() < 1e-6, "{volume}");
        assert_eq!(model.unique_subshapes(solid, NodeKind::Face).len(), 4);
        let wire = disc(&mut model, Point::ORIGIN, Direction::X, 1.0, false);
        let shell = pipe_cornered(&mut model, wire, spine, &PipeOptions::default(), Corner::Mitred).unwrap();
        assert_eq!(shell.kind(), NodeKind::Shell);
        assert_eq!(model.unique_subshapes(shell, NodeKind::Face).len(), 2);
    }

    #[test]
    fn a_ring_mitred_round_a_square() {
        // Round a closed square of sides 10: four cylinders meeting at four mitres, no ends.
        let mut model = Model::new();
        let corners = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)].map(|(x, y)| Point::new(x, y, 0.0));
        let spine = polyline(&mut model, &corners);
        let face = disc(&mut model, Point::new(5.0, 0.0, 0.0), Direction::X, 1.0, true);
        let solid = pipe_cornered(&mut model, face, spine, &PipeOptions::default(), Corner::Mitred).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid());
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - std::f64::consts::PI * 40.0).abs() < 1e-6, "{volume}");
        assert_eq!(model.unique_subshapes(solid, NodeKind::Face).len(), 4);
    }
}
