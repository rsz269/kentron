//! The middle path of a pipe-like shape, between two of its sections.
//!
//! As OCCT's `BRepOffsetAPI_MiddlePath`: from each vertex of the start section a path of the
//! shape's edges runs (one step at a time, each to a vertex not met before) to the end
//! section; at each step the edges between the paths' vertices make a section, and the path
//! wanted runs through the sections' centres (of area where a section is planar, of length
//! otherwise). Between two sections it is a line where every path's edge is one, an arc
//! where every path's edge is an arc about the same axis through the same angle, and
//! otherwise a curve interpolated through the centres, along the paths' mean direction.

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::{fit, ParametricCurve};
use crate::math::{Direction, Frame, Point, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};
use crate::{make, properties};

/// The middle path of `shape` from section `start` to section `end` (each a wire, a face
/// (its outer wire) or an edge of `shape`): a wire.
#[doc(alias = "BRepOffsetAPI_MiddlePath")]
#[doc(alias = "middlepath")]
pub fn middle_path(model: &mut Model, shape: Shape, start: Shape, end: Shape) -> Result<Shape> {
    let invalid = |m: &str| Error::InvalidInput(format!("middle path: {m}"));
    let start_edges = section_edges(model, start)?;
    let end_edges = section_edges(model, end)?;
    let ends_of = |model: &Model, e: Shape| {
        let (a, b) = model.edge_vertices(e);
        (a.node, b.node)
    };
    let closed = ends_of(model, start_edges[0]).0 == ends_of(model, *start_edges.last().expect("edges")).1;
    let mut start_vertices: Vec<NodeId> = start_edges.iter().map(|&e| ends_of(model, e).0).collect();
    if !closed {
        start_vertices.push(ends_of(model, *start_edges.last().expect("edges")).1);
    }
    let end_vertices: Vec<NodeId> = end_edges.iter().flat_map(|&e| { let (a, b) = ends_of(model, e); [a, b] }).collect();
    let in_start = |e: NodeId| start_edges.iter().any(|s| s.node == e);
    // The shape's edges at each vertex.
    let edges: Vec<Shape> = model.unique_subshapes(shape, NodeKind::Edge).into_iter().map(|e| Shape { orientation: Orientation::Forward, ..e }).collect();
    let at = |v: NodeId| -> Vec<Shape> {
        edges.iter().copied().filter(|&e| { let (a, b) = ends_of(model, e); a == v || b == v }).collect()
    };
    // An edge oriented away from `v`, with its far vertex.
    let away = |model: &Model, e: Shape, v: NodeId| -> (Shape, NodeId) {
        let (a, b) = ends_of(model, e);
        if a == v { (e, b) } else { (e.reversed(), a) }
    };
    // The paths: from each start vertex, edge by edge (staying at a vertex for a step where
    // more than one edge leads on to a vertex not met yet, as OCCT's paths do).
    let mut seen: Vec<NodeId> = start_vertices.clone();
    let mut paths: Vec<Vec<Step>> = Vec::new();
    for &v in &start_vertices {
        let first = at(v).into_iter().find(|e| !in_start(e.node)).ok_or_else(|| invalid("a start vertex with no edge leaving the section"))?;
        let (e, next) = away(model, first, v);
        seen.push(next);
        paths.push(vec![Step::Edge(e, next)]);
    }
    loop {
        let mut grown = Vec::new();
        for path in &mut paths {
            let v = path.last().expect("a path").end();
            if end_vertices.contains(&v) {
                continue;
            }
            let last_edge = path.iter().rev().find_map(|s| match s { Step::Edge(e, _) => Some(e.node), Step::Stay(_) => None });
            let candidates: Vec<(Shape, NodeId)> = at(v).into_iter().filter(|e| Some(e.node) != last_edge).map(|e| away(model, e, v)).filter(|(_, w)| !seen.contains(w)).collect();
            match candidates[..] {
                [] => {}
                [(e, w)] => {
                    path.push(Step::Edge(e, w));
                    grown.push(w);
                }
                _ => path.push(Step::Stay(v)),
            }
        }
        if grown.is_empty() {
            break;
        }
        seen.extend(grown);
    }
    // (Steps where every path stayed are none; paths that got to the end first stay there.)
    let steps = paths.iter().map(|p| p.len()).max().unwrap_or(0);
    for path in &mut paths {
        while path.len() < steps {
            let v = path.last().expect("a path").end();
            path.push(Step::Stay(v));
        }
    }
    let mut keep: Vec<bool> = (0..steps).map(|i| paths.iter().any(|p| matches!(p[i], Step::Edge(..)))).collect();
    if keep.is_empty() {
        keep.push(true);
    }
    let paths: Vec<Vec<Step>> = paths.into_iter().map(|p| p.into_iter().zip(&keep).filter(|(_, k)| **k).map(|(s, _)| s).collect()).collect();
    let steps = paths[0].len();
    // The sections: the start, the edges between the paths' vertices after each step, the end.
    let mut centres = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let section: Vec<Shape> = if i == 0 {
            start_edges.clone()
        } else if i == steps {
            end_edges.clone()
        } else {
            let points: Vec<NodeId> = paths.iter().map(|p| p[i - 1].end()).collect();
            let lateral: Vec<NodeId> = paths.iter().flat_map(|p| p.iter().filter_map(|s| s.edge().map(|e| e.node))).collect();
            let n = points.len();
            let links = if closed { n } else { n - 1 };
            let mut found = Vec::with_capacity(links);
            for j in 0..links {
                let (a, b) = (points[j], points[(j + 1) % n]);
                if a == b {
                    continue;
                }
                // (Where no edge joins them, a line: OCCT draws one across their face.)
                let edge = edges.iter().copied().find(|e| !lateral.contains(&e.node) && { let (p, q) = ends_of(model, *e); (p == a && q == b) || (p == b && q == a) });
                let edge = match edge {
                    Some(e) => e,
                    None => {
                        let point = |v: NodeId| model.point(Shape::from(crate::topology::VertexId::try_from(v).expect("a vertex")));
                        let (p, q) = (point(a), point(b));
                        make::edge_between(model, p, q)?
                    }
                };
                found.push(edge);
            }
            found
        };
        centres.push(centre(model, &section, closed)?);
    }
    // The path's edges.
    let mut pieces: Vec<Option<Shape>> = vec![None; steps];
    for (i, piece) in pieces.iter_mut().enumerate() {
        let kinds: Vec<Option<Curve>> = paths.iter().map(|p| p[i].edge().and_then(|e| model.edge_curve(e)).map(|(c, _, _)| basis(c))).collect();
        if kinds.iter().all(|k| matches!(k, Some(Curve::Line(_)))) {
            if centres[i].distance(centres[i + 1]) > crate::tolerance::LINEAR {
                *piece = Some(make::edge_between(model, centres[i], centres[i + 1])?);
            }
        } else if kinds.iter().all(|k| matches!(k, Some(Curve::Circle(_)))) {
            *piece = arc(model, &paths, i, centres[i], centres[i + 1])?;
        }
    }
    // The rest interpolated through the centres, along the paths' mean direction.
    let mut i = 0;
    while i < steps {
        if pieces[i].is_some() {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < steps && pieces[j].is_none() {
            j += 1;
        }
        let points: Vec<Point> = centres[i..=j].to_vec();
        let tangents: Vec<Option<Vector>> = (i..=j)
            .map(|k| {
                let mut sum = Vector::ZERO;
                for p in &paths {
                    if k > i || k == j {
                        if let Some(e) = (k > 0).then(|| p[k - 1].edge()).flatten() { sum = sum + direction(model, e, false); }
                    }
                    if k < j || k == i {
                        if let Some(e) = (k < steps).then(|| p[k].edge()).flatten() { sum = sum + direction(model, e, true); }
                    }
                }
                (sum.length() > 1e-12).then(|| sum.normalized())
            })
            .collect();
        let curve = fit::interpolate_with_tangents(&points, None, &tangents, false, true, 1e-9)?;
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        pieces[i] = Some(make::edge_from_curve_range(model, Curve::BSpline(curve), a, b)?);
        i = j;
    }
    let edges: Vec<Shape> = pieces.into_iter().flatten().collect();
    if edges.is_empty() {
        return Err(invalid("the sections' centres coincide"));
    }
    make::wire(model, &edges)
}

/// A step of a path: along an edge (oriented along the path) to a vertex, or staying at one.
#[derive(Clone, Copy)]
enum Step {
    Edge(Shape, NodeId),
    Stay(NodeId),
}

impl Step {
    fn end(&self) -> NodeId {
        match *self {
            Step::Edge(_, v) | Step::Stay(v) => v,
        }
    }

    fn edge(&self) -> Option<Shape> {
        match *self {
            Step::Edge(e, _) => Some(e),
            Step::Stay(_) => None,
        }
    }
}

/// A section's edges, in order along it.
fn section_edges(model: &mut Model, section: Shape) -> Result<Vec<Shape>> {
    let wire = match section.kind() {
        NodeKind::Wire => section,
        NodeKind::Face => model.children(section).next().ok_or_else(|| Error::InvalidInput("middle path: a face with no wire".into()))?,
        NodeKind::Edge => make::wire(model, &[section])?,
        _ => return Err(Error::InvalidInput("middle path: a section is a wire, a face or an edge".into())),
    };
    let edges = model.ordered_edges(wire);
    if edges.is_empty() {
        return Err(Error::InvalidInput("middle path: an empty section".into()));
    }
    Ok(edges)
}

/// A section's centre: of the area it bounds where it is planar (an open one closed by a
/// line), of its length otherwise.
fn centre(model: &mut Model, section: &[Shape], closed: bool) -> Result<Point> {
    let mut edges = section.to_vec();
    if !closed {
        let first = model.point(model.edge_vertices(section[0]).0);
        let last = model.point(model.edge_vertices(*section.last().expect("edges")).1);
        if first.distance(last) > crate::tolerance::LINEAR {
            edges.push(make::edge_between(model, last, first)?);
        }
    }
    let wire = make::wire(model, &edges)?;
    if let Ok(face) = make::planar_face(model, wire) {
        let p = properties::area(model, face);
        if p.mass.abs() > 1e-12 {
            return Ok(p.center_of_mass);
        }
    }
    Ok(properties::length(model, wire).center_of_mass)
}

/// A trimmed curve's basis, other curves as they are.
fn basis(curve: Curve) -> Curve {
    match curve {
        Curve::Trimmed(t) => t.basis().clone(),
        c => c,
    }
}

/// The unit tangent of `edge` (as oriented) at its start (`at_start`) or end.
fn direction(model: &Model, edge: Shape, at_start: bool) -> Vector {
    let Some((curve, first, last)) = model.edge_curve(edge) else { return Vector::ZERO };
    let reversed = edge.orientation == Orientation::Reversed;
    let t = if at_start != reversed { first } else { last };
    let d = curve.d1(t).1;
    let d = if reversed { -d } else { d };
    if d.length() > 1e-300 { d.normalized() } else { Vector::ZERO }
}

/// The arc between two centres where every path's edge at `step` is an arc about the same
/// axis through the same angle (`None` otherwise).
fn arc(model: &mut Model, paths: &[Vec<Step>], step: usize, from: Point, to: Point) -> Result<Option<Shape>> {
    let mut axis: Option<(Point, Vector, f64)> = None;
    for p in paths {
        let Some(edge) = p[step].edge() else { return Ok(None) };
        let Some((Curve::Circle(c), first, last)) = model.edge_curve(edge).map(|(c, a, b)| (basis(c), a, b)) else { return Ok(None) };
        let (o, z) = (c.center(), Vector::from(c.frame().z()));
        let angle = (last - first).abs();
        match axis {
            None => axis = Some((o, z, angle)),
            Some((o0, z0, a0)) => {
                let off = (o - o0) - z0 * (o - o0).dot(z0);
                if z.cross(z0).length() > 1e-7 || off.length() > 1e-5 || (angle - a0).abs() > 1e-7 {
                    return Ok(None);
                }
            }
        }
    }
    let Some((o, z, angle)) = axis else { return Ok(None) };
    // About the axis through the start's centre, turning the way the paths do.
    let centre = o + z * (from - o).dot(z);
    let v1 = from - centre;
    if v1.length() <= crate::tolerance::LINEAR {
        return Ok(None);
    }
    let tangent = paths[0][step].edge().map_or(Vector::ZERO, |e| direction(model, e, true));
    let z = if z.cross(v1).dot(tangent) < 0.0 { -z } else { z };
    let circle = Circle::new(Frame::new(centre, Direction::from_vector(z), Direction::from_vector(v1)), v1.length());
    let _ = to;
    Ok(Some(make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, angle)?))
}
