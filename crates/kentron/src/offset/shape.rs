//! Offsets of solids and shells, keeping their topology.
//!
//! Each face moves along its outward normal: its surface is replaced by the offset surface
//! (a plane shifted, a cylinder, sphere, cone or torus with its radius changed). Where faces
//! meet, the new edge is where their offset surfaces meet (joined by intersection), traced
//! point by point in the faces' parameters from the old edge and extended or cut back to the
//! new vertices, each found where its faces' offset surfaces meet. Offsets small enough to
//! leave every face in place are supported; faces vanishing or new ones appearing (as offsets
//! grow past features) are not.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{Curve, Line};
use crate::geometry::curve2::{BSplineCurve2, Curve2};
use crate::geometry::extrema::nearest_on_curve;
use crate::geometry::fit;
use crate::geometry::surface::{OffsetSurface, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Point, Point2, Transformable, Vector};
use crate::topology::{ByOccurrence, CoEdgeId, EdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

use super::Join;

/// `shape` (a solid or shell) with each face moved `distance` along its outward normal
/// (inwards for a negative distance), its faces meeting where their offsets intersect.
#[doc(alias = "BRepOffsetAPI_MakeOffsetShape")]
#[doc(alias = "BRepOffset_MakeOffset")]
pub fn shape(model: &mut Model, shape: Shape, distance: f64, join: Join) -> Result<Shape> {
    shape_with_faces(model, shape, distance, &[], join)
}

/// [`shape`] with some faces (of `shape`) moved by distances of their own.
#[doc(alias = "BRepOffset_MakeOffset::SetOffsetOnFace")]
#[doc(alias = "offsetonface")]
pub fn shape_with_faces(model: &mut Model, shape: Shape, distance: f64, faces: &[(Shape, f64)], join: Join) -> Result<Shape> {
    let own: Vec<(NodeId, f64)> = faces.iter().map(|(f, d)| (f.node, *d)).collect();
    let moved = offset_faces(model, shape, distance, &own, &[], join)?;
    let result = rebuild(model, shape, &moved.faces, &moved.new_faces)?;
    let radius = moved.parting.first().map_or(0.0, |&(a, _)| own.iter().find(|(n, _)| *n == moved.faces[a].shape.node).map_or(distance, |(_, d)| *d).abs());
    let result = round_parting(model, result, &moved, radius)?;
    super::regularize::regularized(model, result)
}

/// `shape` (a solid) hollowed out to walls of `thickness`, open where `removed` (faces of
/// it) were: the walls grow outwards for a positive thickness, inwards for a negative one.
///
/// The walls are bounded by the solid's faces and their offsets (the removed faces staying
/// where they are), and closed at each removed face by a rim: of the face and its moved
/// version, the larger with the other cut out.
#[doc(alias = "BRepOffsetAPI_MakeThickSolid")]
#[doc(alias = "thickshell")]
pub fn thick_solid(model: &mut Model, shape: Shape, removed: &[Shape], thickness: f64, join: Join) -> Result<Shape> {
    // A compound of one solid: that solid.
    let solids = model.unique_subshapes(shape, NodeKind::Solid);
    let shape = match (shape.kind(), &solids[..]) {
        (NodeKind::Solid, _) => shape,
        (_, [solid]) => *solid,
        _ => return Err(Error::InvalidInput("a thick solid of a shape that is not one solid".into())),
    };
    // The walls reuse the solid's faces: placed faces are first moved where they are.
    if shape.placement != crate::math::Transform::IDENTITY {
        let unplaced = Shape { placement: crate::math::Transform::IDENTITY, ..shape };
        let baked = crate::modify::transform(model, unplaced, &shape.placement, Default::default())?;
        let removed: Vec<Shape> = removed
            .iter()
            .map(|f| corresponding_face(model, shape, baked, *f))
            .collect::<Result<_>>()?;
        return thick_solid(model, baked, &removed, thickness, join);
    }
    let kept: Vec<NodeId> = removed.iter().map(|f| f.node).collect();
    let moved = offset_faces(model, shape, thickness, &[], &kept, join)?;
    let outwards = thickness > 0.0;
    let shell = model.new_empty(NodeKind::Shell);
    for (i, face) in moved.faces.iter().enumerate() {
        let old = face.shape;
        let new = Shape { orientation: old.orientation, ..moved.new_faces[i] };
        if !kept.contains(&old.node) {
            // The wall's outside and inside.
            let (outside, inside) = if outwards { (new, old) } else { (old, new) };
            model.add(shell, outside)?;
            model.add(shell, inside.reversed())?;
            continue;
        }
        // The rim: the larger face, the smaller one cut out of it.
        let (larger, smaller) = if outwards { (moved.new_faces[i], Shape { orientation: Orientation::Forward, ..old }) } else {
            let copy = model.new_empty(NodeKind::Face);
            model.set_face_surface(copy, face.surface.clone(), model.tolerance(old))?;
            for wire in copy_wires(model, Shape { orientation: Orientation::Forward, ..old }, false)? {
                model.add(copy, wire.into())?;
            }
            (copy, moved.new_faces[i])
        };
        for hole in copy_wires(model, smaller, true)? {
            model.add(larger, hole.into())?;
        }
        model.add(shell, Shape { orientation: old.orientation, ..larger })?;
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    let solid = round_parting(model, solid, &moved, thickness.abs())?;
    super::regularize::regularized(model, solid)
}

/// The face of `copy` (`original` moved where it is placed) standing for `face` of
/// `original`: the one at the same position in the order of their faces.
fn corresponding_face(model: &Model, original: Shape, copy: Shape, face: Shape) -> Result<Shape> {
    let faces = model.unique_subshapes(original, NodeKind::Face);
    let i = faces
        .iter()
        .position(|f| f.node == face.node)
        .ok_or_else(|| Error::InvalidInput("a removed face not of the solid".into()))?;
    model.unique_subshapes(copy, NodeKind::Face).get(i).copied().ok_or_else(|| Error::InvalidInput("a face lost moving the solid".into()))
}

/// New wires through the edges of `face`'s wires (as it is used: forward), run backwards if
/// `reverse` (to cut the face out of another one on the same surface).
fn copy_wires(model: &mut Model, face: Shape, reverse: bool) -> Result<Vec<crate::topology::WireId>> {
    let mut wires = Vec::new();
    for wire in model.children(face).collect::<Vec<_>>() {
        let mut order: Vec<Shape> = model.children(wire).collect();
        if wire.orientation == Orientation::Reversed {
            order.reverse();
        }
        let mut coedges = Vec::with_capacity(order.len());
        for coedge in order {
            let id = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
            let orientation = model[id].orientation.compose(coedge.orientation);
            let edge = model[id].edge;
            let (pcurve, first, last) = model.effective_pcurve(id).ok_or_else(|| Error::InvalidInput("a face without pcurves".into()))?;
            let new = model.add_coedge(edge, if reverse { orientation.reversed() } else { orientation })?;
            model.set_pcurve_with_range(new, pcurve, first, last)?;
            coedges.push(new);
        }
        if reverse {
            coedges.reverse();
        }
        wires.push(model.add_wire(&coedges)?);
    }
    Ok(wires)
}

/// The faces of a shape and their offsets (forward), in the same order.
struct Moved {
    faces: Vec<OldFace>,
    new_faces: Vec<Shape>,
    /// The faces (by index) of the edges whose faces' offsets part, joined by arcs: their
    /// offsets meet sharply first, the edge rounded after.
    parting: Vec<(usize, usize)>,
}

/// A face's use of an edge: the face (by index), the coedge and its curve in the face.
struct Use {
    face: usize,
    pcurve: Curve2,
    orientation: Orientation,
}

struct OldEdge {
    shape: Shape,
    curve: Option<Curve>,
    first: f64,
    last: f64,
    start: usize,
    end: usize,
    uses: Vec<Use>,
}

struct OldFace {
    shape: Shape,
    surface: Surface,
    offset: Surface,
    /// For a cone, the side of its apex the face is on (the sign of the radius over it).
    nappe: f64,
}

/// The offset of every face of `shape` (by `distance`, or by its own in `own`) but those in
/// `kept` (which keep their surfaces), with the same topology.
fn offset_faces(model: &mut Model, shape: Shape, distance: f64, own: &[(NodeId, f64)], kept: &[NodeId], join: Join) -> Result<Moved> {
    if distance == 0.0 && own.is_empty() {
        return Err(Error::InvalidInput("an offset by zero".into()));
    }
    let distance_of = |node: NodeId| -> f64 {
        if kept.contains(&node) {
            return 0.0;
        }
        own.iter().find(|(n, _)| *n == node).map_or(distance, |(_, d)| *d)
    };
    let new_surface = |face: Shape, surface: &Surface| -> Result<Surface> {
        let d = distance_of(face.node);
        let along = if face.orientation == Orientation::Reversed { -d } else { d };
        if d == 0.0 { Ok(surface.clone()) } else { offset_surface(surface, along) }
    };
    // Joined by arcs, faces whose offsets part (at convex edges going out, concave ones going
    // in) are joined by rolling-ball faces; elsewhere arcs and intersections agree. The
    // offsets meet sharply there first, and the edges are rounded after (by the ball's
    // radius: its centre then runs along the old edge, as the arc join's does).
    let parting = std::cell::RefCell::new(Vec::new());
    let check = |faces: &[OldFace], edges: &[OldEdge]| -> Result<()> {
        if join == Join::Arc {
            for e in edges.iter().filter(|e| parts(faces, e, &distance_of)) {
                let (a, b) = (e.uses[0].face, e.uses[1].face);
                if distance_of(faces[a].shape.node) != distance_of(faces[b].shape.node) {
                    todo!("offsets joined by arcs between faces moved by different distances");
                }
                parting.borrow_mut().push((a, b));
            }
        }
        Ok(())
    };
    let scale = own.iter().map(|(_, d)| d.abs()).fold(distance.abs(), f64::max);
    let mut moved = replace_surfaces(model, shape, &new_surface, &check, scale)?;
    moved.parting = parting.into_inner();
    Ok(moved)
}

/// `shape` (holding `moved`'s new faces) with the edges where the offsets of parting faces
/// meet rounded by `radius`.
fn round_parting(model: &mut Model, shape: Shape, moved: &Moved, radius: f64) -> Result<Shape> {
    if moved.parting.is_empty() {
        return Ok(shape);
    }
    let edges_of = |model: &Model, face: Shape| -> Vec<NodeId> { model.unique_subshapes(face, NodeKind::Edge).into_iter().map(|e| e.node).collect() };
    let mut wanted: Vec<NodeId> = Vec::new();
    for &(a, b) in &moved.parting {
        let (ea, eb) = (edges_of(model, moved.new_faces[a]), edges_of(model, moved.new_faces[b]));
        for n in ea {
            if eb.contains(&n) && !wanted.contains(&n) {
                wanted.push(n);
            }
        }
    }
    // (The edges as the result holds them.)
    let edges: Vec<Shape> = model.unique_subshapes(shape, NodeKind::Edge).into_iter().filter(|e| wanted.contains(&e.node)).collect();
    if edges.is_empty() {
        return Ok(shape);
    }
    // Closed shells outside solids are rounded as the solids they bound (the fillets cut
    // and add material), and given back as shells.
    if model.unique_subshapes(shape, NodeKind::Solid).is_empty() {
        let shells = model.unique_subshapes(shape, NodeKind::Shell);
        let compound = model.new_empty(NodeKind::Compound);
        for shell in shells {
            let solid = model.new_empty(NodeKind::Solid);
            model.add(solid, shell)?;
            model.add(compound, solid)?;
        }
        let rounded = crate::blend::fillet(model, compound, &edges, radius)?;
        let shells = model.unique_subshapes(rounded, NodeKind::Shell);
        return match (shape.kind(), &shells[..]) {
            (NodeKind::Shell, [one]) => Ok(*one),
            _ => {
                let out = model.new_empty(NodeKind::Compound);
                for s in shells {
                    model.add(out, s)?;
                }
                Ok(out)
            }
        };
    }
    crate::blend::fillet(model, shape, &edges, radius)
}

/// `shape` with its faces' surfaces replaced by `new_surface` (given each face as the shape
/// uses it and its surface, placed), keeping the topology: the vertices and edges found
/// again where the new surfaces meet.
pub(crate) fn with_surfaces(model: &mut Model, shape: Shape, new_surface: &dyn Fn(Shape, &Surface) -> Result<Surface>) -> Result<Shape> {
    let moved = replace_surfaces(model, shape, new_surface, &|_, _| Ok(()), 0.0)?;
    rebuild(model, shape, &moved.faces, &moved.new_faces)
}

/// The faces of `shape` on the surfaces `new_surface` gives them, with the same topology:
/// `check` sees the faces and edges before anything is built; `scale` is a length the
/// tolerances grow with (beside the shape's size).
fn replace_surfaces(
    model: &mut Model,
    shape: Shape,
    new_surface: &dyn Fn(Shape, &Surface) -> Result<Surface>,
    check: &dyn Fn(&[OldFace], &[OldEdge]) -> Result<()>,
    scale: f64,
) -> Result<Moved> {
    // Faces and their new surfaces.
    let mut faces: Vec<OldFace> = Vec::new();
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.transformed(&face.placement);
        let offset = new_surface(face, &surface)?;
        let nappe = match &offset {
            Surface::Cone(c) => {
                let (_, _, v0, v1) = model.face_uv_bounds(face);
                (c.radius() + 0.5 * (v0 + v1) * c.semi_angle_sin_cos().0).signum()
            }
            _ => 0.0,
        };
        faces.push(OldFace { shape: face, surface, offset, nappe });
    }
    // Edges, with their uses by the faces; vertices.
    let mut vertices: Vec<(Shape, Point)> = Vec::new();
    let mut vertex_of: ByOccurrence<usize> = ByOccurrence::default();
    let mut edges: Vec<OldEdge> = Vec::new();
    let mut edge_of: ByOccurrence<usize> = ByOccurrence::default();
    for (fi, face) in faces.iter().enumerate() {
        let forward = Shape { orientation: Orientation::Forward, ..face.shape };
        for wire in model.children(forward) {
            let mut order: Vec<Shape> = model.children(wire).collect();
            if wire.orientation == Orientation::Reversed {
                order.reverse();
            }
            for coedge in order {
                let id = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                let edge_shape = Shape { placement: coedge.placement, ..Shape::from(model[id].edge) };
                let index = match edge_of.get(&edge_shape) {
                    Some(&i) => i,
                    None => {
                        let (a, b) = model.edge_vertices(edge_shape);
                        let mut vertex = |v: Shape| match vertex_of.get(&v) {
                            Some(&i) => i,
                            None => {
                                vertices.push((v, model.point(v)));
                                vertex_of.insert(&v, vertices.len() - 1);
                                vertices.len() - 1
                            }
                        };
                        let (start, end) = (vertex(a), vertex(b));
                        let (curve, first, last) = match model.edge_curve(edge_shape) {
                            Some((c, f, l)) => (Some(c), f, l),
                            None => (None, 0.0, 0.0),
                        };
                        edges.push(OldEdge { shape: edge_shape, curve, first, last, start, end, uses: Vec::new() });
                        edge_of.insert(&edge_shape, edges.len() - 1);
                        edges.len() - 1
                    }
                };
                let (pcurve, first, last) = model
                    .effective_pcurve(id)
                    .ok_or_else(|| Error::InvalidInput("an offset of a face without pcurves".into()))?;
                if edges[index].curve.is_none() {
                    (edges[index].first, edges[index].last) = (first, last);
                }
                let orientation = model[id].orientation.compose(coedge.orientation);
                edges[index].uses.push(Use { face: fi, pcurve, orientation });
            }
        }
    }
    check(&faces, &edges)?;
    // The new vertices: where the new surfaces of their faces meet.
    let mut around: Vec<Vec<(usize, Point2)>> = vec![Vec::new(); vertices.len()];
    for edge in &edges {
        for u in &edge.uses {
            for (v, t) in [(edge.start, edge.first), (edge.end, edge.last)] {
                let p = u.pcurve.point(t);
                if !around[v].iter().any(|&(f, q)| f == u.face && q.distance(p) <= 1e-9) {
                    around[v].push((u.face, p));
                }
            }
        }
    }
    let scale = vertices.iter().map(|(_, p)| p.distance(vertices[0].1)).fold(scale, f64::max);
    let new_points: Vec<Point> = (0..vertices.len())
        .map(|v| {
            // A vertex of a cone face alone (its apex) goes to the offset cone's apex.
            let mut on: Vec<usize> = around[v].iter().map(|&(f, _)| f).collect();
            on.dedup();
            if let ([f], Surface::Cone(c)) = (&on[..], &faces[around[v][0].0].offset) {
                let _ = f;
                let (frame, angle) = (c.frame(), c.semi_angle());
                return Ok(frame.origin() + Vector::from(frame.z()) * (-c.radius() / angle.tan()));
            }
            meet(&faces, &around[v], scale)
        })
        .collect::<Result<_>>()?;
    let tolerance = 1e-7 * scale.max(1.0);
    let new_vertices: Vec<VertexId> = new_points.iter().map(|&p| model.add_vertex(p, tolerance)).collect();
    // The new edges, with their pcurves on the faces using them.
    let mut new_edges: Vec<(EdgeId, Vec<(Curve2, f64, f64)>)> = Vec::with_capacity(edges.len());
    for edge in &edges {
        new_edges.push(trace_edge(model, &faces, edge, &new_points, &new_vertices, tolerance)?);
    }
    // The faces, wired as before.
    let mut use_index: HashMap<usize, usize> = HashMap::new();
    let mut new_faces: Vec<Shape> = Vec::with_capacity(faces.len());
    for (fi, face) in faces.iter().enumerate() {
        let new = model.new_empty(NodeKind::Face);
        model.set_face_surface(new, face.offset.clone(), tolerance)?;
        let forward = Shape { orientation: Orientation::Forward, ..face.shape };
        for wire in model.children(forward).collect::<Vec<_>>() {
            let mut coedges = Vec::new();
            // A wire used reversed runs through its coedges backwards.
            let mut order: Vec<Shape> = model.children(wire).collect();
            if wire.orientation == Orientation::Reversed {
                order.reverse();
            }
            for coedge in order {
                let id = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                let ei = *edge_of.get(&Shape { placement: coedge.placement, ..Shape::from(model[id].edge) }).expect("edge");
                // This face's next use of the edge (a seam is used twice).
                let k = use_index.entry(ei).or_insert(0);
                let uses: Vec<usize> = (0..edges[ei].uses.len()).filter(|&i| edges[ei].uses[i].face == fi).collect();
                let which = uses[*k % uses.len()];
                if uses.len() > 1 {
                    *k += 1;
                }
                let (edge_id, pcurves) = &new_edges[ei];
                if pcurves.is_empty() {
                    continue; // a degenerate edge gone (its face no longer reaches the pole)
                }
                let (pcurve, first, last) = pcurves[which].clone();
                let new_coedge = model.add_coedge(*edge_id, edges[ei].uses[which].orientation)?;
                model.set_pcurve_with_range(new_coedge, pcurve, first, last)?;
                coedges.push(new_coedge);
            }
            use_index.retain(|&e, _| !edges[e].uses.iter().all(|u| u.face == fi));
            let new_wire = model.add_wire(&coedges)?;
            model.add(new, new_wire.into())?;
        }
        new_faces.push(new);
    }
    Ok(Moved { faces, new_faces, parting: Vec::new() })
}

/// Whether the offsets of an edge's two (moved) faces part, leaving a gap there: at a convex
/// edge for an outward offset, a concave one for an inward one.
fn parts(faces: &[OldFace], edge: &OldEdge, distance_of: &dyn Fn(NodeId) -> f64) -> bool {
    let Some(curve) = &edge.curve else { return false };
    let [u1, u2] = &edge.uses[..] else { return false };
    let (d1, d2) = (distance_of(faces[u1.face].shape.node), distance_of(faces[u2.face].shape.node));
    if u1.face == u2.face || d1 == 0.0 || d2 == 0.0 || (d1 > 0.0) != (d2 > 0.0) {
        return false;
    }
    let distance = d1;
    let t = 0.5 * (edge.first + edge.last);
    // Outward normals, and the edge's direction along the first face's boundary.
    let outward = |u: &Use| {
        let f = &faces[u.face];
        let uv = u.pcurve.point(t);
        let (_, du, dv) = f.surface.d1(uv.x, uv.y);
        let n = du.cross(dv);
        if f.shape.orientation == Orientation::Reversed { -n } else { n }
    };
    let (n1, n2) = (outward(u1).normalized(), outward(u2).normalized());
    let along = curve.d1(t).1;
    let reversed = (u1.orientation == Orientation::Reversed) != (faces[u1.face].shape.orientation == Orientation::Reversed);
    let along = if reversed { -along } else { along };
    let turn = n1.cross(n2);
    if turn.length() <= 1e-6 {
        return false; // tangent faces
    }
    let convex = turn.dot(along) > 0.0;
    convex == (distance > 0.0)
}

/// The containers of the new faces, nested like the old ones.
fn rebuild(model: &mut Model, shape: Shape, faces: &[OldFace], new_faces: &[Shape]) -> Result<Shape> {
    let face_index = |s: &Shape| faces.iter().position(|f| f.shape.node == s.node && f.shape.placement == s.placement);
    let mut shells = Vec::new();
    for shell in model.unique_subshapes(shape, NodeKind::Shell) {
        let new_shell = model.new_empty(NodeKind::Shell);
        // Each face as the shell uses it.
        let forward = Shape { orientation: Orientation::Forward, ..shell };
        for face in model.children(forward).collect::<Vec<_>>() {
            let i = face_index(&face).ok_or_else(|| Error::InvalidInput("a shell's face not offset".into()))?;
            model.add(new_shell, Shape { orientation: face.orientation, ..new_faces[i] })?;
        }
        shells.push((shell, new_shell));
    }
    match shape.kind() {
        NodeKind::Solid => {
            let solid = model.new_empty(NodeKind::Solid);
            for shell in model.children(shape).collect::<Vec<_>>() {
                let (_, new_shell) = shells.iter().find(|(s, _)| s.node == shell.node).ok_or_else(|| Error::InvalidInput("a solid without shells".into()))?;
                model.add(solid, Shape { orientation: shell.orientation, ..*new_shell })?;
            }
            Ok(solid)
        }
        NodeKind::Shell => Ok(shells[0].1),
        NodeKind::Face => Ok(Shape { orientation: faces[0].shape.orientation, ..new_faces[0] }),
        _ => {
            // A compound: of its solids (or shells).
            let compound = model.new_empty(NodeKind::Compound);
            for solid in model.unique_subshapes(shape, NodeKind::Solid) {
                let new_solid = model.new_empty(NodeKind::Solid);
                for shell in model.children(solid).collect::<Vec<_>>() {
                    if let Some((_, new_shell)) = shells.iter().find(|(s, _)| s.node == shell.node) {
                        model.add(new_solid, Shape { orientation: shell.orientation, ..*new_shell })?;
                    }
                }
                model.add(compound, new_solid)?;
            }
            // The shells in no solid.
            let in_solids: Vec<NodeId> = model.unique_subshapes(shape, NodeKind::Solid).into_iter().flat_map(|s| model.children(s).map(|c| c.node).collect::<Vec<_>>()).collect();
            for (shell, new_shell) in &shells {
                if !in_solids.contains(&shell.node) {
                    model.add(compound, Shape { orientation: shell.orientation, ..*new_shell })?;
                }
            }
            Ok(compound)
        }
    }
}

/// The surface `distance` away along its normal `∂u × ∂v`, parametrized like it (each point
/// of the offset at the parameters of the point it comes from): a plane, quadric or torus
/// where the offset is one, a general offset surface otherwise.
fn offset_surface(surface: &Surface, distance: f64) -> Result<Surface> {
    let offset = OffsetSurface::new(surface.clone(), distance);
    let analytic = |s: &Surface| matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_));
    let elementary = match surface {
        Surface::Trimmed(t) => analytic(t.basis()),
        s => analytic(s),
    };
    match offset.equivalent_surface() {
        Some(s) => Ok(s),
        None if elementary => Err(Error::InvalidInput("an offset collapsing a face to its axis".into())),
        None => Ok(Surface::Offset(offset)),
    }
}

/// The point of `surface` nearest `p`, by Gauss-Newton steps from `uv`: its parameters, the
/// point and the unit normal there.
fn foot(surface: &Surface, p: Point, uv: Point2) -> (Point2, Point, Vector) {
    foot_on(surface, p, uv, 0.0)
}

/// [`foot`], staying on one side of a cone's apex (`nappe`, the sign of the radius; 0 for
/// either side).
fn foot_on(surface: &Surface, p: Point, mut uv: Point2, nappe: f64) -> (Point2, Point, Vector) {
    let keep = |uv: Point2| match surface {
        Surface::Cone(c) if nappe != 0.0 => {
            let apex = -c.radius() / c.semi_angle_sin_cos().0;
            let on = (c.radius() + uv.y * c.semi_angle_sin_cos().0) * nappe;
            if on < 0.0 { Point2::new(uv.x, apex + nappe.signum() * c.semi_angle_sin_cos().0.signum() * 1e-6 * (1.0 + apex.abs())) } else { uv }
        }
        _ => uv,
    };
    // Closed-form parameters where the surface has them (exact at poles, where steps on
    // the parameters stall), moved by whole periods next to the previous ones.
    if let Some((u, v)) = surface.parameters(p) {
        let near = |x: f64, previous: f64, period: Option<f64>| match period {
            Some(t) => x + ((previous - x) / t).round() * t,
            None => x,
        };
        uv = Point2::new(near(u, uv.x, surface.u_period()), near(v, uv.y, surface.v_period()));
    }
    uv = keep(uv);
    for _ in 0..50 {
        let (q, du, dv) = surface.d1(uv.x, uv.y);
        let r = q - p;
        let (a, b, c) = (du.dot(du), du.dot(dv), dv.dot(dv));
        let (g1, g2) = (r.dot(du), r.dot(dv));
        let det = a * c - b * b;
        // At a pole (one derivative vanishing) only the other parameter moves the point.
        let (x, y) = if det.abs() > 1e-12 * a.max(c).powi(2) {
            ((-g1 * c + g2 * b) / det, (-g2 * a + g1 * b) / det)
        } else if c > a && c > 0.0 {
            (0.0, -g2 / c)
        } else if a > 0.0 {
            (-g1 / a, 0.0)
        } else {
            break;
        };
        uv = keep(Point2::new(uv.x + x, uv.y + y));
        if x.abs() + y.abs() <= 1e-15 * (1.0 + uv.x.abs() + uv.y.abs()) {
            break;
        }
    }
    let (q, du, dv) = surface.d1(uv.x, uv.y);
    let n = du.cross(dv);
    let n = if n.length() > 0.0 { n.normalized() } else { n };
    (uv, q, n)
}

/// Where the offset surfaces of a vertex's faces meet, from its parameters on each: the point
/// on all of them (least squares on their tangent planes, nearest the start where they leave
/// it free).
fn meet(faces: &[OldFace], around: &[(usize, Point2)], scale: f64) -> Result<Point> {
    let mut uvs: Vec<(usize, Point2)> = around.to_vec();
    let starts: Vec<Point> = uvs.iter().map(|&(f, uv)| faces[f].offset.point(uv.x, uv.y)).collect();
    let mut x = Point::ORIGIN + starts.iter().fold(Vector::ZERO, |acc, p| acc + (*p - Point::ORIGIN)) / starts.len() as f64;
    let x0 = x;
    let weight = 1e-12;
    for _ in 0..60 {
        // Normal equations of Σ (n·(x + δ - q))² + w |x + δ - x0|².
        let mut m = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        for (f, uv) in uvs.iter_mut() {
            let (new_uv, q, n) = foot_on(&faces[*f].offset, x, *uv, faces[*f].nappe);
            *uv = new_uv;
            let n = [n.x, n.y, n.z];
            let r = (q - x).dot(Vector::new(n[0], n[1], n[2]));
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] += n[i] * n[j];
                }
                rhs[i] += n[i] * r;
            }
        }
        let back = x0 - x;
        for i in 0..3 {
            m[i][i] += weight;
        }
        rhs[0] += weight * back.x;
        rhs[1] += weight * back.y;
        rhs[2] += weight * back.z;
        let delta = solve3(m, rhs).ok_or_else(|| Error::InvalidInput("an offset vertex not found".into()))?;
        x = x + delta;
        if delta.length() <= 1e-13 * scale.max(1.0) {
            break;
        }
    }
    Ok(x)
}

fn solve3(m: [[f64; 3]; 3], r: [f64; 3]) -> Option<Vector> {
    let det = |a: [[f64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d = det(m);
    if d.abs() <= 1e-300 {
        return None;
    }
    let column = |k: usize| {
        let mut a = m;
        for i in 0..3 {
            a[i][k] = r[i];
        }
        det(a) / d
    };
    Some(Vector::new(column(0), column(1), column(2)))
}

/// The new edge standing for `edge`: between its new vertices, on the offset surfaces of its
/// faces (where they meet, or along the face for a seam or between tangent faces), with a
/// pcurve per use.
fn trace_edge(
    model: &mut Model,
    faces: &[OldFace],
    edge: &OldEdge,
    points: &[Point],
    vertices: &[VertexId],
    tolerance: f64,
) -> Result<(EdgeId, Vec<(Curve2, f64, f64)>)> {
    let (a, b) = (points[edge.start], points[edge.end]);
    let id = model.add_edge(vertices[edge.start], vertices[edge.end], tolerance)?;
    let Some(curve) = &edge.curve else {
        // Degenerate (a pole or apex): where the new vertex is on each face, across the same
        // range of the other parameter; gone (no pcurves) where the face no longer reaches
        // its surface's pole (the offsets of the neighbouring faces meet away from it).
        let singular = edge.uses.iter().all(|u| {
            let (uv, _, _) = foot(&faces[u.face].offset, a, u.pcurve.point(edge.first));
            let (_, du, dv) = faces[u.face].offset.d1(uv.x, uv.y);
            du.length().min(dv.length()) <= 1e-6 * du.length().max(dv.length())
        });
        if !singular {
            return Ok((id, Vec::new()));
        }
        let pcurves = edge
            .uses
            .iter()
            .map(|u| {
                let (p0, p1) = (u.pcurve.point(edge.first), u.pcurve.point(edge.last));
                let (uv, _, _) = foot(&faces[u.face].offset, a, p0);
                let (s, t) = if (p1.x - p0.x).abs() >= (p1.y - p0.y).abs() {
                    (Point2::new(p0.x, uv.y), Point2::new(p1.x, uv.y))
                } else {
                    (Point2::new(uv.x, p0.y), Point2::new(uv.x, p1.y))
                };
                (segment2(s, t, edge.first, edge.last), edge.first, edge.last)
            })
            .collect();
        return Ok((id, pcurves));
    };
    // The parameters of the old curve where the new vertices are (beyond its ends where the
    // edge grows); a closed edge keeps its range.
    let (t0, t1) = if edge.start == edge.end {
        (edge.first, edge.last)
    } else {
        // Near the old end: within a few times the distance moved, in parameter (a closed
        // curve comes back to the point a period away).
        let near = |p: Point, t: f64| {
            let speed = curve.d1(t).1.length().max(1e-300);
            let reach = 4.0 * p.distance(curve.point(t)) / speed + 1e-9 * (1.0 + t.abs());
            nearest_on_curve(p, curve, t - reach, t + reach).parameter
        };
        (near(a, edge.first), near(b, edge.last))
    };
    let distinct: Vec<usize> = {
        let mut f: Vec<usize> = edge.uses.iter().map(|u| u.face).collect();
        f.sort_unstable();
        f.dedup();
        f
    };
    // (A trimmed plane is a plane: the face bounds it.)
    let plane = |s: &Surface| match s {
        Surface::Trimmed(t) => matches!(t.basis(), Surface::Plane(_)),
        s => matches!(s, Surface::Plane(_)),
    };
    let straight = matches!(curve, Curve::Line(_)) && distinct.iter().all(|&f| plane(&faces[f].offset));
    let n = if straight { 1 } else { 32 };
    let ts: Vec<f64> = (0..=n).map(|i| t0 + (t1 - t0) * i as f64 / n as f64).collect();
    // Per sample: the point and the parameters on each use's face.
    let mut xs: Vec<Point> = Vec::with_capacity(ts.len());
    let mut uvs: Vec<Vec<Point2>> = vec![Vec::with_capacity(ts.len()); edge.uses.len()];
    for &t in &ts {
        // Each use's parameters: those of the old pcurve (offsets keep the parametrization),
        // corrected where two faces' offsets cross.
        let mut on: Vec<Point2> = edge.uses.iter().map(|u| u.pcurve.point(t)).collect();
        let x = if distinct.len() == 2 {
            let i = edge.uses.iter().position(|u| u.face == distinct[0]).expect("use");
            let j = edge.uses.iter().position(|u| u.face == distinct[1]).expect("use");
            let (fi, fj) = (&faces[distinct[0]].offset, &faces[distinct[1]].offset);
            let (p, d) = curve.d1(t);
            match intersect(fi, fj, on[i], on[j], p, d) {
                Some((x, ui, uj)) => {
                    (on[i], on[j]) = (ui, uj);
                    x
                }
                // Tangent faces: the nearest point of one new surface (for an offset, the
                // edge moved along the common normal), at its parameters on both.
                None => {
                    let (ui, x, _) = foot(fi, p, on[i]);
                    on[i] = ui;
                    on[j] = foot(fj, x, on[j]).0;
                    x
                }
            }
        } else {
            // A seam (both uses on one face) or an edge of one face: along the face.
            faces[distinct[0]].offset.point(on[0].x, on[0].y)
        };
        xs.push(x);
        for (k, uv) in on.into_iter().enumerate() {
            uvs[k].push(uv);
        }
    }
    // The ends exactly at the new vertices, also in each face's parameters (the traced
    // parameters are poor at a pole).
    xs[0] = a;
    xs[n] = b;
    for (k, u) in edge.uses.iter().enumerate() {
        let f = &faces[u.face];
        uvs[k][0] = foot_on(&f.offset, a, uvs[k][0], f.nappe).0;
        uvs[k][n] = foot_on(&f.offset, b, uvs[k][n], f.nappe).0;
    }
    // Between two elementary surfaces, the conic they meet in, where it runs through both new
    // vertices: exact, where the samples taken along the old edge may miss part of the new
    // one (a hyperbola turning back past the old edge's end).
    let mut conic = if n > 1 && distinct.len() == 2 && coaxial_circle(curve, &xs).is_none() {
        exact_section(&faces[distinct[0]].offset, &faces[distinct[1]].offset, a, b, xs[n / 2])
    } else {
        None
    };
    if let Some((c, first, last)) = &conic {
        // The samples again, along the conic, each one's parameters found from the last; their
        // ends where the traced ones are (whole periods away, the parameters moved by them;
        // else the traced edge kept).
        let points: Vec<Point> = (0..=n).map(|i| if i == 0 { a } else if i == n { b } else { c.point(first + (last - first) * i as f64 / n as f64) }).collect();
        let mut marched: Vec<Vec<Point2>> = Vec::with_capacity(edge.uses.len());
        let mut agree = true;
        for (k, u) in edge.uses.iter().enumerate() {
            let f = &faces[u.face];
            let mut uv = vec![uvs[k][0]];
            for i in 1..=n {
                uv.push(foot_on(&f.offset, points[i], uv[i - 1], f.nappe).0);
            }
            let whole = |d: f64, period: Option<f64>| match period {
                Some(t) => (d / t).round() * t,
                None => 0.0,
            };
            let (dx, dy) = (uvs[k][n].x - uv[n].x, uvs[k][n].y - uv[n].y);
            let shift = Point2::new(whole(dx, f.offset.u_period()), whole(dy, f.offset.v_period())) - Point2::new(0.0, 0.0);
            if (dx - shift.x).abs() + (dy - shift.y).abs() > 1e-6 * (1.0 + uv[n].x.abs() + uv[n].y.abs()) || shift.x != 0.0 || shift.y != 0.0 {
                // (A pcurve ending a period away runs round the face: not this edge.)
                agree = false;
                break;
            }
            marched.push(uv);
        }
        if agree {
            xs = points;
            uvs = marched;
        } else {
            conic = None;
        }
    }
    let curve3 = if n == 1 {
        Curve::Line(Line { origin: a, direction: Direction::try_from_vector(b - a).ok_or_else(|| Error::InvalidInput("an offset edge of no length".into()))? })
    } else if let Some(circle) = coaxial_circle(curve, &xs) {
        // (A circle between surfaces about its axis moves to a circle about it, at the same
        // angles.)
        Curve::Circle(circle)
    } else if let Some((c, _, _)) = &conic {
        c.clone()
    } else {
        Curve::BSpline(fit::interpolate_with_parameters(&xs, &ts)?)
    };
    // Lines by length, conics over their own parameters; other curves over the old curve's
    // parameters, which keep increasing from the start vertex to the end one unless the
    // offset turned the edge round.
    let (first, last, params) = match &conic {
        _ if n == 1 => (0.0, a.distance(b), vec![0.0, a.distance(b)]),
        Some((_, first, last)) => (*first, *last, (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect()),
        None => (t0, t1, ts.clone()),
    };
    if first >= last {
        return Err(Error::InvalidInput("an offset edge turned round (the offset is too large)".into()));
    }
    model.set_edge_curve(id.into(), curve3, tolerance)?;
    model.set_edge_range(id.into(), first, last)?;
    let pcurves = uvs
        .into_iter()
        .map(|uv| {
            let pcurve = if uv.len() == 2 {
                segment2(uv[0], uv[1], params[0], params[1])
            } else {
                Curve2::BSpline(fit::interpolate_with_parameters(&uv, &params).expect("increasing parameters"))
            };
            (pcurve, first, last)
        })
        .collect();
    Ok((id, pcurves))
}

/// Where elementary surfaces `s1` and `s2` meet in a conic through `a` and `b`: the conic,
/// run from `a` to `b` (over the arc nearest `middle` on a closed one), and its parameters
/// there.
fn exact_section(s1: &Surface, s2: &Surface, a: Point, b: Point, middle: Point) -> Option<(Curve, f64, f64)> {
    use crate::geometry::intersect::{quadrics, QuadricIntersection as Q};
    let elementary = |s: &Surface| matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_));
    if !elementary(s1) || !elementary(s2) {
        return None;
    }
    let scale = 1.0 + (a - Point::ORIGIN).length().max((b - Point::ORIGIN).length());
    let tolerance = 1e-9 * scale;
    let curves: Vec<Curve> = match quadrics(s1, s2, tolerance) {
        Q::Circles(v) => v.into_iter().map(Curve::Circle).collect(),
        Q::Ellipses(v) => v.into_iter().map(Curve::Ellipse).collect(),
        Q::Parabolas(v) => v.into_iter().map(Curve::Parabola).collect(),
        Q::Hyperbolas(v) => v.into_iter().map(Curve::Hyperbola).collect(),
        _ => return None,
    };
    // A point's parameter on a conic, in closed form (in the conic's frame).
    let parameter = |c: &Curve, p: Point| -> Option<f64> {
        let local = |frame: &crate::math::Frame| {
            let d = p - frame.origin();
            (d.dot(Vector::from(frame.x())), d.dot(Vector::from(frame.y())))
        };
        Some(match c {
            Curve::Circle(c) => {
                let (x, y) = local(c.frame());
                y.atan2(x).rem_euclid(std::f64::consts::TAU)
            }
            Curve::Ellipse(e) => {
                let (x, y) = local(e.frame());
                (y / e.minor_radius()).atan2(x / e.major_radius()).rem_euclid(std::f64::consts::TAU)
            }
            Curve::Hyperbola(h) => (local(h.frame()).1 / h.minor_radius()).asinh(),
            Curve::Parabola(q) => local(q.frame()).1,
            _ => return None,
        })
    };
    for mut c in curves {
        let (Some(ta), Some(tb)) = (parameter(&c, a), parameter(&c, b)) else { continue };
        if c.point(ta).distance(a) > 1e3 * tolerance || c.point(tb).distance(b) > 1e3 * tolerance {
            continue;
        }
        let range = match c.period() {
            // (Of the two arcs, the one through the middle of the traced edge.)
            Some(period) => {
                let forward = (ta, if tb > ta { tb } else { tb + period });
                let mid = |(x, y): (f64, f64)| c.point(0.5 * (x + y)).distance(middle);
                let mut turned = c.clone();
                turned.reverse();
                let (ra, rb) = (turned.reversed_parameter(ta).rem_euclid(period), turned.reversed_parameter(tb).rem_euclid(period));
                let backward = (ra, if rb > ra { rb } else { rb + period });
                if mid(forward) <= turned.point(0.5 * (backward.0 + backward.1)).distance(middle) {
                    forward
                } else {
                    c = turned;
                    backward
                }
            }
            None if ta < tb => (ta, tb),
            None => {
                let (ra, rb) = (c.reversed_parameter(ta), c.reversed_parameter(tb));
                c.reverse();
                (ra, rb)
            }
        };
        // (Where it is the edge: through the traced middle, not another stretch of it.)
        let (lo, hi) = range;
        let near_middle = c.point(nearest_on_curve(middle, &c, lo, hi).parameter).distance(middle);
        if near_middle <= 1e-6 * scale {
            return Some((c, lo, hi));
        }
    }
    None
}

/// The circle about the axis of `curve` (a circle) through `points`, if they all lie at one
/// height along it and one distance from it; its X axis is the old circle's, so that
/// points keep their angles.
fn coaxial_circle(curve: &Curve, points: &[Point]) -> Option<crate::geometry::curve::Circle> {
    let Curve::Circle(c) = curve else { return None };
    let frame = c.frame();
    let (o, z) = (frame.origin(), Vector::from(frame.z()));
    let (heights, radii): (Vec<f64>, Vec<f64>) = points
        .iter()
        .map(|&p| {
            let v = p - o;
            let h = v.dot(z);
            (h, (v - z * h).length())
        })
        .unzip();
    let size = c.radius().max(radii[0]).max(1e-300);
    let tol = 1e-9 * (size + heights[0].abs());
    let same = |xs: &[f64]| xs.iter().all(|x| (x - xs[0]).abs() <= tol);
    if !(radii[0] > tol) || !same(&heights) || !same(&radii) {
        return None;
    }
    let new = crate::math::Frame::new(o + z * heights[0], frame.z(), frame.x());
    Some(crate::geometry::curve::Circle::new(new, radii[0]))
}

fn segment2(a: Point2, b: Point2, t0: f64, t1: f64) -> Curve2 {
    Curve2::BSpline(BSplineCurve2::new(vec![a, b], vec![t0, t1], vec![2, 2], 1))
}

/// The point where two surfaces meet in the plane through `p` across `d`, by Newton steps on
/// their parameters from `a` and `b`; `None` if they do not cross there (tangent surfaces).
fn intersect(s1: &Surface, s2: &Surface, a: Point2, b: Point2, p: Point, d: Vector) -> Option<(Point, Point2, Point2)> {
    let (mut z1, mut z2) = (a, b);
    // Tangent along the edge: no crossing to find.
    let (_, du1, dv1) = s1.d1(a.x, a.y);
    let (_, du2, dv2) = s2.d1(b.x, b.y);
    let (n1, n2) = (du1.cross(dv1), du2.cross(dv2));
    if n1.cross(n2).length() <= 1e-9 * n1.length() * n2.length() {
        return None;
    }
    for _ in 0..50 {
        let (x1, u1, v1) = s1.d1(z1.x, z1.y);
        let (x2, u2, v2) = s2.d1(z2.x, z2.y);
        let f = x1 - x2;
        let g = (x1 - p).dot(d);
        // Jacobian columns: ∂/∂u1, ∂/∂v1, ∂/∂u2, ∂/∂v2 of (f, g).
        let jac = [
            [u1.x, v1.x, -u2.x, -v2.x],
            [u1.y, v1.y, -u2.y, -v2.y],
            [u1.z, v1.z, -u2.z, -v2.z],
            [u1.dot(d), v1.dot(d), 0.0, 0.0],
        ];
        let step = solve4(jac, [-f.x, -f.y, -f.z, -g])?;
        z1 = Point2::new(z1.x + step[0], z1.y + step[1]);
        z2 = Point2::new(z2.x + step[2], z2.y + step[3]);
        if step.iter().map(|s| s.abs()).sum::<f64>() <= 1e-14 * (1.0 + z1.x.abs() + z1.y.abs() + z2.x.abs() + z2.y.abs()) {
            break;
        }
    }
    let x = s1.point(z1.x, z1.y);
    (x.distance(s2.point(z2.x, z2.y)) <= 1e-7 * (1.0 + x.distance(p))).then_some((x, z1, z2))
}

/// Gaussian elimination with partial pivoting.
fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() <= 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..4 {
            let f = a[row][col] / a[col][col];
            for k in col..4 {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0; 4];
    for row in (0..4).rev() {
        let mut s = b[row];
        for k in row + 1..4 {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cone, Cuboid, Cylinder};
    use crate::properties::volume;

    fn valid(model: &Model, s: Shape) {
        let problems: Vec<String> = crate::validity::check(model, s).errors().map(|p| p.description.clone()).collect();
        assert!(problems.is_empty(), "{problems:?}");
    }

    #[test]
    fn a_box_grows_with_rounded_edges() {
        // Joined by arcs, a box grown by r is its Minkowski sum with a ball: faces moved,
        // quarter cylinders along the edges, eighth spheres at the corners.
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        let r = 1.0;
        let grown = shape(&mut model, b, r, Join::Arc).unwrap();
        valid(&model, grown);
        let expected = 1000.0 + 600.0 * r + 120.0 * std::f64::consts::PI * r * r / 4.0 + 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
        let v = volume(&model, grown).mass;
        assert!((v - expected).abs() < 1e-6 * expected, "{v} != {expected}");
        assert_eq!(model.unique_subshapes(grown, NodeKind::Face).len(), 6 + 12 + 8);
    }

    #[test]
    fn a_face_moves_by_its_own_distance() {
        // A 10 x 10 x 10 box grown by 1, its top by 2: 12 x 12 x 13.
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        let top = model
            .unique_subshapes(b, NodeKind::Face)
            .into_iter()
            .find(|&f| model.unique_subshapes(f, NodeKind::Vertex).iter().all(|&v| model.point(v).z == 10.0))
            .unwrap();
        let grown = shape_with_faces(&mut model, b, 1.0, &[(top, 2.0)], Join::Intersection).unwrap();
        valid(&model, grown);
        assert!((volume(&model, grown).mass - 12.0 * 12.0 * 13.0).abs() < 1e-9);
    }

    #[test]
    fn a_box_grows_and_shrinks_by_its_faces() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let out = shape(&mut model, b, 1.0, Join::Intersection).unwrap();
        assert!((volume(&model, out).mass - 12.0 * 22.0 * 32.0).abs() < 1e-6);
        valid(&model, out);
        let inside = shape(&mut model, b, -1.0, Join::Intersection).unwrap();
        assert!((volume(&model, inside).mass - 8.0 * 18.0 * 28.0).abs() < 1e-6);
        // Arc joins agree with intersections where no faces part (a box shrinking).
        let arcs = shape(&mut model, b, -1.0, Join::Arc).unwrap();
        assert!((volume(&model, arcs).mass - 8.0 * 18.0 * 28.0).abs() < 1e-6);
    }

    #[test]
    fn a_cylinder_and_a_cone_offset_outwards() {
        let mut model = Model::new();
        let c = Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
        let out = shape(&mut model, c, 1.0, Join::Intersection).unwrap();
        let expected = std::f64::consts::PI * 36.0 * 17.0;
        assert!((volume(&model, out).mass - expected).abs() < 1e-6 * expected);
        valid(&model, out);
        // A cone grows to a new apex: 2.6 higher (1 / sin of its half angle), its base 1 lower.
        let cone = Cone::new(5.0, 0.0, 12.0).build(&mut model).unwrap();
        let out = shape(&mut model, cone, 1.0, Join::Intersection).unwrap();
        let expected = std::f64::consts::PI / 3.0 * 6.5 * 6.5 * 15.6;
        assert!((volume(&model, out).mass - expected).abs() < 1e-4 * expected, "{}", volume(&model, out).mass);
    }

    #[test]
    fn a_cup_from_a_cylinder() {
        let mut model = Model::new();
        let c = Cylinder::new(5.0, 15.0).build(&mut model).unwrap();
        // The top face (at z = 15) removed, walls 1 thick outside.
        let top = model
            .unique_subshapes(c, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Plane(p)) if p.origin().z > 10.0))
            .unwrap();
        let cup = thick_solid(&mut model, c, &[top], 1.0, Join::Intersection).unwrap();
        let expected = std::f64::consts::PI * (36.0 * 16.0 - 25.0 * 15.0);
        // Within the interpolation of the traced edges' pcurves.
        assert!((volume(&model, cup).mass - expected).abs() < 1e-5 * expected, "{}", volume(&model, cup).mass);
        valid(&model, cup);
        // Inside: walls 1 thick within the cylinder.
        let inner = thick_solid(&mut model, c, &[top], -1.0, Join::Intersection).unwrap();
        let expected = std::f64::consts::PI * (25.0 * 15.0 - 16.0 * 14.0);
        assert!((volume(&model, inner).mass - expected).abs() < 1e-5 * expected, "{}", volume(&model, inner).mass);
    }
}
