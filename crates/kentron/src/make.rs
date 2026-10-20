//! Elementary topology construction: vertices, edges, wires.

use crate::error::{Error, Result};
use crate::geometry::curve::{Curve, Line};
use crate::geometry::surface::{Plane, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Point, Vector};
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId, WireId};

/// The vertex named by `shape`.
fn vertex_id(model: &Model, shape: Shape) -> Result<VertexId> {
    let id = VertexId::try_from(shape.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a vertex")))?;
    model.get(id).ok_or_else(|| Error::InvalidInput(format!("no {id:?}")))?;
    Ok(id)
}

/// A new vertex at `p`. Two calls with the same point make two distinct vertices.
#[doc(alias = "BRepBuilderAPI_MakeVertex")]
pub fn vertex(model: &mut Model, p: Point) -> Shape {
    model.add_vertex(p, tolerance::LINEAR).into()
}

/// A straight edge between two existing vertices. Fails if they coincide.
#[doc(alias = "BRepBuilderAPI_MakeEdge")]
pub fn edge(model: &mut Model, start: Shape, end: Shape) -> Result<Shape> {
    let (start, end) = (vertex_id(model, start)?, vertex_id(model, end)?);
    let (a, b) = (model[start].point, model[end].point);
    let direction = Direction::try_from_vector(b - a)
        .filter(|_| a.distance(b) > tolerance::LINEAR)
        .ok_or_else(|| Error::InvalidInput("edge between coincident points".into()))?;
    let edge = model.add_edge(start, end, tolerance::LINEAR)?;
    model.set_edge_curve(edge.into(), Curve::Line(Line::new(a, direction)), tolerance::LINEAR)?;
    model.set_edge_range(edge.into(), 0.0, a.distance(b))?;
    Ok(edge.into())
}

/// A straight edge between two points, with new vertices. Fails if they coincide.
pub fn edge_between(model: &mut Model, start: Point, end: Point) -> Result<Shape> {
    if start.distance(end) <= tolerance::LINEAR {
        return Err(Error::InvalidInput("edge between coincident points".into()));
    }
    let (a, b) = (vertex(model, start), vertex(model, end));
    edge(model, a, b)
}

/// An edge along the whole of `curve` (closed if the curve is closed, e.g. a circle).
/// Fails for unbounded curves such as lines.
pub fn edge_from_curve(model: &mut Model, curve: Curve) -> Result<Shape> {
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    if first <= -tolerance::INFINITE || last >= tolerance::INFINITE {
        return Err(Error::InvalidInput("an edge needs a bounded curve".into()));
    }
    edge_from_curve_range(model, curve, first, last)
}

/// An edge along `curve` restricted to `[first, last]`. For a periodic curve, equal
/// parameters mean one full period starting there.
pub fn edge_from_curve_range(model: &mut Model, curve: Curve, first: f64, last: f64) -> Result<Shape> {
    let (first, last) = match curve.period() {
        Some(period) if (last - first).abs() <= tolerance::PARAMETRIC => (first, first + period),
        _ => (first, last),
    };
    let slack = tolerance::PARAMETRIC;
    let in_domain = curve.period().is_some()
        || (first >= curve.first_parameter() - slack && last <= curve.last_parameter() + slack);
    if !(first < last) || !in_domain || first <= -tolerance::INFINITE || last >= tolerance::INFINITE {
        return Err(Error::InvalidInput(format!("invalid edge range [{first}, {last}]")));
    }
    let (a, b) = (curve.point(first), curve.point(last));
    let start = model.add_vertex(a, tolerance::LINEAR);
    // A closed curve (a full circle) starts and ends at the same vertex.
    let end = if a.distance(b) <= tolerance::LINEAR { start } else { model.add_vertex(b, tolerance::LINEAR) };
    let edge = model.add_edge(start, end, tolerance::LINEAR)?;
    model.set_edge_curve(edge.into(), curve, tolerance::LINEAR)?;
    model.set_edge_range(edge.into(), first, last)?;
    Ok(edge.into())
}

/// A wire through `edges`, in order. Consecutive edges must meet: share a vertex, or have
/// end vertices within their tolerances, which are then merged (the later edge is made to
/// use the earlier one's vertex), as are the wire's two ends if they meet.
#[doc(alias = "BRepBuilderAPI_MakeWire")]
pub fn wire(model: &mut Model, edges: &[Shape]) -> Result<Shape> {
    wire_within(model, edges, 0.0)
}

/// [`wire`], ends also counting as meeting where they lie within `tolerance` of each other
/// (as OCCT's `ShapeAnalysis_FreeBounds` joins them): the vertex kept where two are merged
/// across a gap has its tolerance raised to cover it.
pub fn wire_within(model: &mut Model, edges: &[Shape], tolerance: f64) -> Result<Shape> {
    let gap_tolerance = tolerance;
    let ends = |model: &Model, e: EdgeId| -> Result<(VertexId, VertexId)> {
        let edge = model.get(e).ok_or_else(|| Error::InvalidInput(format!("no {e:?}")))?;
        let vertex = |u| model.get(u).map(|u: &crate::topology::VertexUse| u.vertex);
        match (vertex(edge.start), vertex(edge.end)) {
            (Some(a), Some(b)) => Ok((a, b)),
            _ => Err(Error::InvalidInput(format!("{e:?} has no vertices"))),
        }
    };
    // Vertices that coincide within their tolerances count as the same (and get merged).
    let meet = |model: &Model, a: VertexId, b: VertexId| {
        let (va, vb) = (&model[a], &model[b]);
        a == b || va.point.distance(vb.point) <= va.tolerance.max(vb.tolerance).max(tolerance::LINEAR).max(gap_tolerance)
    };
    // Merging `gone` into `kept` across a gap, for edge `edge`: `kept`'s tolerance and the
    // edge's cover the gap.
    let cover = |model: &mut Model, kept: VertexId, gone: VertexId, edge: EdgeId| -> Result<()> {
        let (k, g) = (&model[kept], &model[gone]);
        let gap = k.point.distance(g.point) + g.tolerance;
        if gap > k.tolerance {
            model.set_vertex_tolerance(kept, gap)?;
        }
        if gap > model[edge].tolerance {
            model.set_edge_tolerance(edge, gap)?;
        }
        Ok(())
    };
    // A wire holds its edges where they are: a placed edge is copied with its placement
    // applied.
    let mut placed = Vec::with_capacity(edges.len());
    for &e in edges {
        placed.push(if e.placement == crate::math::Transform::IDENTITY { e } else { crate::modify::copy(model, e, false)? });
    }
    let ids: Vec<EdgeId> = placed
        .iter()
        .map(|e| EdgeId::try_from(e.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not an edge"))))
        .collect::<Result<_>>()?;
    if ids.is_empty() {
        return Err(Error::InvalidInput("a wire needs at least one edge".into()));
    }
    // Orient each edge to continue from the previous one; the first edge runs as given,
    // unless only the other way leads to the second.
    let mut orientations = Vec::with_capacity(ids.len());
    let (mut start, mut end) = ends(model, ids[0])?;
    if placed[0].orientation == Orientation::Reversed {
        (start, end) = (end, start);
    }
    if let Some(&next) = ids.get(1) {
        let (a, b) = ends(model, next)?;
        if !meet(model, end, a) && !meet(model, end, b) && (meet(model, start, a) || meet(model, start, b)) {
            (start, end) = (end, start);
        }
    }
    let first_start = start;
    let natural = ends(model, ids[0])?;
    // (A closed edge keeps its given direction.)
    let forward = if natural.0 == natural.1 { placed[0].orientation != Orientation::Reversed } else { (start, end) == natural };
    orientations.push(if forward { Orientation::Forward } else { Orientation::Reversed });
    for &id in &ids[1..] {
        let (a, b) = ends(model, id)?;
        let (orientation, joined, next_end) = if meet(model, end, a) {
            (Orientation::Forward, true, b)
        } else if meet(model, end, b) {
            (Orientation::Reversed, false, a)
        } else {
            return Err(Error::InvalidInput(format!("{id:?} does not continue the wire")));
        };
        // Make the edge use the wire's vertex where it only coincides with it.
        let record = &model[id];
        let joining_use = if joined { record.start } else { record.end };
        cover(model, end, if joined { a } else { b }, id)?;
        model.replace_edge_vertex(id, joining_use, end)?;
        orientations.push(orientation);
        end = next_end;
    }
    // Close the wire on its first vertex if its ends meet.
    if ids.len() > 1 && end != first_start && meet(model, end, first_start) {
        let last = *ids.last().expect("edges");
        let record = &model[last];
        let closing_use = if *orientations.last().expect("edges") == Orientation::Forward { record.end } else { record.start };
        cover(model, first_start, end, last)?;
        model.replace_edge_vertex(last, closing_use, first_start)?;
    }
    let coedges = ids
        .iter()
        .zip(orientations)
        .map(|(&e, o)| model.add_coedge(e, o))
        .collect::<Result<Vec<_>>>()?;
    Ok(model.add_wire(&coedges)?.into())
}

/// A face bounded by `wire`, on the plane containing it. Fails unless the wire is closed
/// and planar.
#[doc(alias = "BRepBuilderAPI_MakeFace")]
pub fn planar_face(model: &mut Model, wire: Shape) -> Result<Shape> {
    // A placed wire: the face of the wire where it is defined, placed like it.
    if wire.placement != crate::math::Transform::IDENTITY {
        let face = planar_face(model, Shape { placement: crate::math::Transform::IDENTITY, ..wire })?;
        return Ok(face.moved(&wire.placement));
    }
    let id = WireId::try_from(wire.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a wire")))?;
    if model.get(id).is_none() {
        return Err(Error::InvalidInput(format!("no {id:?}")));
    }
    if !model.is_closed_wire(id) {
        return Err(Error::InvalidInput("a planar face needs a closed wire".into()));
    }
    let (points, normal, tolerance) = wire_polygon(model, wire)?;
    let normal = Direction::try_from_vector(normal)
        .filter(|_| normal.length() > tolerance)
        .ok_or_else(|| Error::InvalidInput("the wire encloses no area".into()))?;
    let plane = Plane::new(points[0], normal);
    if points.iter().any(|&p| plane.distance(p) > tolerance) {
        return Err(Error::InvalidInput("the wire is not planar".into()));
    }
    face_with_boundary(model, Surface::Plane(plane), wire)
}

/// A face on the plane of the open `wire`, bounded by it as it is: an invalid face (its
/// boundary does not close), as OCCT's `BRepBuilderAPI_MakeFace` makes from an open wire.
/// For callers that report a face's validity rather than fail on it (CadQuery's
/// `Face.makeFromWires`); [`planar_face`] fails instead.
pub fn open_face(model: &mut Model, wire: Shape) -> Result<Shape> {
    let plane = plane_of_wire(model, wire)?;
    face_with_boundary(model, Surface::Plane(plane), wire)
}

/// The plane containing `wire`, open or closed (an open one as if closed by a segment from
/// its end to its start), its normal such that the wire runs counter-clockwise about it.
/// Fails unless the wire is planar and not straight.
#[doc(alias = "BRepLib_FindSurface")]
pub fn plane_of_wire(model: &Model, wire: Shape) -> Result<Plane> {
    let (points, normal, tolerance) = wire_polygon(model, wire)?;
    let normal = Direction::try_from_vector(normal)
        .filter(|_| normal.length() > tolerance * tolerance)
        .ok_or_else(|| Error::InvalidInput("the wire lies on no one plane".into()))?;
    let plane = Plane::new(points[0], normal);
    if points.iter().any(|&p| plane.distance(p) > tolerance) {
        return Err(Error::InvalidInput("the wire is not planar".into()));
    }
    Ok(plane)
}

/// Points along a closed wire in its direction of travel (a few per curved edge), its
/// Newell normal (pointing so that the wire runs counter-clockwise about it, with the
/// enclosed area as length) and the largest edge tolerance.
fn wire_polygon(model: &Model, wire: Shape) -> Result<(Vec<Point>, Vector, f64)> {
    let mut points = Vec::new();
    let mut tolerance = tolerance::LINEAR;
    for edge in model.explore(wire, NodeKind::Edge) {
        let (curve, first, last) = model
            .edge_curve(edge)
            .ok_or_else(|| Error::InvalidInput("a planar face needs edges with curves".into()))?;
        tolerance = tolerance.max(model.tolerance(edge));
        let samples = if matches!(curve, Curve::Line(_)) { 1 } else { 8 };
        for i in 0..samples {
            let s = i as f64 / samples as f64;
            let t = if edge.orientation == Orientation::Reversed { last - (last - first) * s } else { first + (last - first) * s };
            points.push(curve.point(t));
        }
    }
    let n = points.len();
    let normal = (0..n).fold(Vector::ZERO, |acc, i| {
        let (a, b) = (points[i] - Point::ORIGIN, points[(i + 1) % n] - Point::ORIGIN);
        acc + a.cross(b) * 0.5
    });
    Ok((points, normal, tolerance))
}

/// A face on `surface` bounded by `wire`, with pcurves derived where the surface allows.
fn face_with_boundary(model: &mut Model, surface: Surface, wire: Shape) -> Result<Shape> {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, tolerance::LINEAR)?;
    model.add(face, wire)?;
    Ok(face)
}

/// A face on `surface` bounded by the closed `wire`. On a plane the wire is used as the
/// outer boundary: reversed if it runs clockwise about the plane's normal (like OCCT's
/// `BRepBuilderAPI_MakeFace` with `Inside`); on other surfaces it is used as given.
#[doc(alias = "BRepBuilderAPI_MakeFace")]
pub fn face_on(model: &mut Model, surface: Surface, wire: Shape) -> Result<Shape> {
    let id = WireId::try_from(wire.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a wire")))?;
    if model.get(id).is_none() || !model.is_closed_wire(id) {
        return Err(Error::InvalidInput("a face needs a closed wire".into()));
    }
    let wire = match &surface {
        Surface::Plane(plane) => {
            let (points, normal, tolerance) = wire_polygon(model, wire)?;
            if let Some(&p) = points.iter().find(|&&p| plane.distance(p) > tolerance) {
                return Err(Error::InvalidInput(format!("the wire leaves the plane at {p:?}")));
            }
            if normal.dot(Vector::from(plane.normal())) < 0.0 { wire.reversed() } else { wire }
        }
        _ => wire,
    };
    face_with_boundary(model, surface, wire)
}

/// `face` with `holes` (closed wires lying on its surface, within it) cut out: a new face on
/// the same surface, bounded by the face's wires and the holes, each hole run clockwise in
/// the surface's parameters (as OCCT's `BRepBuilderAPI_MakeFace::Add`).
#[doc(alias = "BRepBuilderAPI_MakeFace::Add")]
pub fn face_with_holes(model: &mut Model, face: Shape, holes: &[Shape]) -> Result<Shape> {
    use crate::geometry::curve2::{BSplineCurve2, Curve2};
    use crate::geometry::ParametricSurface;
    use crate::math::Point2;
    use crate::topology::CoEdgeId;
    // (A placed face or hole: a copy with its placement applied.)
    let face = if face.placement == crate::math::Transform::IDENTITY { face } else { crate::modify::copy(model, face, false)? };
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.clone();
    let new = model.new_empty(NodeKind::Face);
    let tolerance = model.tolerance(face).max(tolerance::LINEAR);
    model.set_face_surface(new, surface.clone(), tolerance)?;
    // The face's own wires, their coedges again with the same curves on the surface.
    for w in model.children(face).collect::<Vec<_>>() {
        let mut coedges = Vec::new();
        // (Its coedges as the wire holds them; the wire's own orientation kept on the new one.)
        for c in model.explore(Shape { orientation: Orientation::Forward, ..w }, NodeKind::CoEdge).collect::<Vec<_>>() {
            let id = CoEdgeId::try_from(c.node).expect("a coedge");
            let (edge, orientation) = (model[id].edge, model[id].orientation.compose(c.orientation));
            let (pcurve, first, last) = model.effective_pcurve(id).ok_or_else(|| Error::InvalidInput("a face's edge without a curve on it".into()))?;
            let new_c = model.add_coedge(edge, orientation)?;
            model.set_pcurve_with_range(new_c, pcurve, first, last)?;
            coedges.push(new_c);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(new, Shape { orientation: w.orientation, ..Shape::from(wire) })?;
    }
    // The holes: each edge's curve on the surface, carried on across the surface's periods
    // from the one before; turned to run clockwise where they run the other way.
    let flip = face.orientation == Orientation::Reversed;
    let (u_period, v_period) = (surface.u_period(), surface.v_period());
    for &hole in holes {
        let hole = if hole.placement == crate::math::Transform::IDENTITY { hole } else { crate::modify::copy(model, hole, false)? };
        let id = WireId::try_from(hole.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a wire")))?;
        if !model.is_closed_wire(id) {
            return Err(Error::InvalidInput("a hole needs a closed wire".into()));
        }
        let mut uses: Vec<(EdgeId, Orientation, BSplineCurve2, f64, f64)> = Vec::new();
        let mut previous: Option<Point2> = None;
        let mut uv: Vec<Point2> = Vec::new();
        // (In the order the wire holds its edges: which way round it runs is decided below.)
        for c in model.explore(Shape { orientation: Orientation::Forward, ..hole }, NodeKind::CoEdge).collect::<Vec<_>>() {
            let cid = CoEdgeId::try_from(c.node).expect("a coedge");
            let (edge, orientation) = (model[cid].edge, model[cid].orientation.compose(c.orientation));
            let (curve, first, last) = model.edge_curve(Shape::from(edge)).ok_or_else(|| Error::InvalidInput("a degenerate edge in a hole".into()))?;
            // (Within the edge's tolerance; else, a wire only near the surface (made on another
            // it approximates), the nearest points of the surface, each found from the last,
            // interpolated: adjacent edges' ends at the same points. The edge's tolerance then
            // raised to how far it strays.)
            let tol = model[edge].tolerance.max(tolerance::LINEAR);
            let mut pcurve = crate::geometry::project::pcurve_on_range(&curve, &surface, first, last, tol)
                .or_else(|| {
                    let n = 64;
                    let ts: Vec<f64> = (0..=n).map(|k| first + (last - first) * k as f64 / n as f64).collect();
                    let mut from: Option<(f64, f64)> = None;
                    let uv: Vec<Point2> = ts
                        .iter()
                        .map(|&t| {
                            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(curve.point(t), &surface, from);
                            from = Some((u, v));
                            Point2::new(u, v)
                        })
                        .collect();
                    crate::geometry::fit::interpolate_with_parameters(&uv, &ts).ok()
                })
                .ok_or_else(|| Error::InvalidInput("a hole's edge does not lie on the face's surface".into()))?;
            let strays = (0..=32)
                .map(|k| {
                    let t = first + (last - first) * k as f64 / 32.0;
                    let q = pcurve.point(t);
                    surface.point(q.x, q.y).distance(curve.point(t))
                })
                .fold(0.0, f64::max);
            if strays > model[edge].tolerance {
                model.set_edge_tolerance(edge, 1.5 * strays)?;
            }
            let (a, b) = model.edge_vertices(Shape::from(edge));
            for v in [a, b] {
                let v = VertexId::try_from(v.node).expect("a vertex");
                if 1.5 * strays > model[v].tolerance {
                    model.set_vertex_tolerance(v, 1.5 * strays)?;
                }
            }
            let (start, end) = if orientation == Orientation::Reversed { (last, first) } else { (first, last) };
            if let Some(target) = previous {
                let at = pcurve.point(start);
                let shift = Point2::new(
                    u_period.map_or(0.0, |p| ((target.x - at.x) / p).round() * p),
                    v_period.map_or(0.0, |p| ((target.y - at.y) / p).round() * p),
                );
                if shift.x != 0.0 || shift.y != 0.0 {
                    let poles = pcurve.poles().iter().map(|q| Point2::new(q.x + shift.x, q.y + shift.y)).collect();
                    pcurve = BSplineCurve2::new(poles, pcurve.knots().to_vec(), pcurve.multiplicities().to_vec(), pcurve.degree());
                }
            }
            previous = Some(pcurve.point(end));
            uv.extend((0..16).map(|k| pcurve.point(start + (end - start) * k as f64 / 16.0)));
            uses.push((edge, orientation, pcurve, first, last));
        }
        let n = uv.len();
        let area: f64 = (0..n).map(|k| uv[k].x * uv[(k + 1) % n].y - uv[(k + 1) % n].x * uv[k].y).sum();
        if (area > 0.0) != flip {
            uses.reverse();
            for u in &mut uses {
                u.1 = u.1.reversed();
            }
        }
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, orientation, pcurve, first, last) in uses {
            let c = model.add_coedge(edge, orientation)?;
            model.set_pcurve_with_range(c, Curve2::BSpline(pcurve), first, last)?;
            coedges.push(c);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(new, Shape::from(wire))?;
    }
    Ok(Shape { orientation: face.orientation, ..new })
}

/// A face covering the whole natural domain of a bounded surface (a sphere, torus, trimmed
/// surface, ...), with its boundary built from the domain's edges (seams and degenerate
/// edges at poles included). Fails for unbounded surfaces such as planes and cylinders.
pub fn face_on_surface(model: &mut Model, surface: Surface) -> Result<Shape> {
    use crate::geometry::ParametricSurface;
    let (u1, u2, v1, v2) = surface.bounds();
    face_on_domain(model, surface, u1, u2, v1, v2)
}

/// A face on the rectangle `[u1, u2] × [v1, v2]` of `surface`'s parameters, bounded by the
/// iso-curves along its sides. A side spanning a whole period of a periodic direction meets
/// the opposite side in a seam, and a side collapsed to a point (a pole) is a degenerate
/// edge.
#[doc(alias = "BRepBuilderAPI_MakeFace")]
pub fn face_on_domain(model: &mut Model, surface: Surface, u1: f64, u2: f64, v1: f64, v2: f64) -> Result<Shape> {
    use crate::geometry::curve2::{Curve2, Line2};
    use crate::geometry::ParametricSurface;
    use crate::math::{Direction2, Point2};
    if ![u1, u2, v1, v2].iter().all(|x| x.is_finite()) || !(u1 < u2 && v1 < v2) {
        return Err(Error::InvalidInput(format!("a face needs a finite parameter domain, got [{u1}, {u2}] × [{v1}, {v2}]")));
    }
    let closes = |period: Option<f64>, a: f64, b: f64| period.is_some_and(|p| ((b - a) - p).abs() <= tolerance::PARAMETRIC);
    let (u_closed, v_closed) = (closes(surface.u_period(), u1, u2), closes(surface.v_period(), v1, v2));
    let corners = [(u1, v1), (u2, v1), (u2, v2), (u1, v2)];
    let mut vertices: Vec<VertexId> = Vec::new();
    for &(u, v) in &corners {
        let p = surface.point(u, v);
        let existing = vertices.iter().copied().find(|&id| model[id].point.distance(p) <= tolerance::LINEAR);
        vertices.push(existing.unwrap_or_else(|| model.add_vertex(p, tolerance::LINEAR)));
    }
    // The four sides, counter-clockwise in (u, v): (iso curve, its range, 2D line, whether
    // it runs backwards, the side whose edge it shares in a seam).
    let sides = [
        (surface.v_iso(v1), (u1, u2), Line2::new(Point2::new(0.0, v1), Direction2::new(1.0, 0.0)), false, None),
        (surface.u_iso(u2), (v1, v2), Line2::new(Point2::new(u2, 0.0), Direction2::new(0.0, 1.0)), false, None),
        (surface.v_iso(v2), (u1, u2), Line2::new(Point2::new(0.0, v2), Direction2::new(1.0, 0.0)), true, v_closed.then_some(0)),
        (surface.u_iso(u1), (v1, v2), Line2::new(Point2::new(u1, 0.0), Direction2::new(0.0, 1.0)), true, u_closed.then_some(1)),
    ];
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface.clone(), tolerance::LINEAR)?;
    let mut coedges = Vec::new();
    let mut edges = Vec::new();
    for (i, (curve, (first, last), line, backwards, seam_of)) in sides.into_iter().enumerate() {
        let edge = match seam_of {
            Some(j) => edges[j],
            None => {
                let (a, b) = (vertices[i], vertices[(i + 1) % 4]);
                let (start, end) = if backwards { (b, a) } else { (a, b) };
                let edge = model.add_edge(start, end, tolerance::LINEAR)?;
                // A side collapsed to a point (a pole) is a degenerate edge.
                let collapsed = curve.point(first).distance(curve.point(last)) <= tolerance::LINEAR
                    && curve.point(0.5 * (first + last)).distance(curve.point(first)) <= tolerance::LINEAR;
                if !collapsed {
                    model.set_edge_curve(edge.into(), curve, tolerance::LINEAR)?;
                    model.set_edge_range(edge.into(), first, last)?;
                }
                edge
            }
        };
        edges.push(edge);
        let coedge = model.add_coedge(edge, if backwards { Orientation::Reversed } else { Orientation::Forward })?;
        model.set_pcurve_with_range(coedge, Curve2::Line(line), first, last)?;
        coedges.push(coedge);
    }
    let wire = model.add_wire(&coedges)?;
    model.add(face, wire.into())?;
    Ok(face)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::NodeKind;

    /// The edges of a box bound its six faces; two circles crossing each other, inside a
    /// third, its regions: the lens, the two crescents, and the rest with both as a hole.
    #[test]
    fn faces_from_box_edges_and_crossing_circles() {
        use crate::geometry::curve::Circle;
        use crate::math::Frame;
        let mut model = Model::new();
        let solid = crate::primitives::Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let edges = model.new_empty(NodeKind::Compound);
        for e in model.unique_subshapes(solid, NodeKind::Edge) {
            model.add(edges, e).unwrap();
        }
        let faces = super::faces_from_edges(&mut model, edges, false).unwrap();
        assert_eq!(model.unique_subshapes(faces, NodeKind::Face).len(), 6);
        let area = crate::properties::area(&model, faces).mass;
        assert!((area - 22.0).abs() < 1e-9, "{area}");
        let circles = model.new_empty(NodeKind::Compound);
        for (x, r) in [(0.0, 10.0), (-2.0, 5.0), (2.0, 5.0)] {
            let c = Circle::new(Frame::new(Point::new(x, 0.0, 0.0), Direction::Z, Direction::X), r);
            let e = super::edge_from_curve(&mut model, Curve::Circle(c)).unwrap();
            model.add(circles, e).unwrap();
        }
        let regions = super::faces_from_edges(&mut model, circles, false).unwrap();
        let faces = model.unique_subshapes(regions, NodeKind::Face);
        assert_eq!(faces.len(), 4);
        assert_eq!(model.unique_subshapes(regions, NodeKind::Wire).len(), 5);
        let total = crate::properties::area(&model, regions).mass;
        assert!((total - std::f64::consts::PI * 100.0).abs() < 1e-6, "{total}");
    }

    /// Separately built edges whose ends coincide become one closed loop.
    #[test]
    fn wire_merges_coincident_vertices() {
        let mut model = Model::new();
        let p = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let edges: Vec<Shape> = (0..3).map(|i| edge_between(&mut model, p[i], p[(i + 1) % 3]).unwrap()).collect();
        let wire = wire(&mut model, &edges).unwrap();
        assert_eq!(model.unique_subshapes(wire, NodeKind::Vertex).len(), 3);
        assert!(model.is_closed_wire(crate::topology::WireId::try_from(wire.node).unwrap()));
    }

    #[test]
    fn wire_rejects_a_gap() {
        let mut model = Model::new();
        let a = edge_between(&mut model, Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)).unwrap();
        let b = edge_between(&mut model, Point::new(2.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0)).unwrap();
        assert!(wire(&mut model, &[a, b]).is_err());
    }

    /// A hole (a patch between two rulings and two circles, run either way) cut in a
    /// cylinder's side: the face loses the patch's area.
    #[test]
    fn holes_are_cut_in_curved_faces() {
        use crate::geometry::curve::Circle;
        use crate::math::Frame;
        let mut model = Model::new();
        let solid = crate::primitives::Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
        let side = model.unique_subshapes(solid, NodeKind::Face).into_iter().find(|&f| matches!(model.surface(f), Some(Surface::Cylinder(_)))).unwrap();
        let (a, b) = (0.5f64, 1.5f64);
        let at = |t: f64, z: f64| Point::new(5.0 * t.cos(), 5.0 * t.sin(), z);
        let arc = |model: &mut Model, z: f64| edge_from_curve_range(model, Curve::Circle(Circle::new(Frame::from_z(Point::new(0.0, 0.0, z), Direction::Z), 5.0)), a, b).unwrap();
        let low = arc(&mut model, 3.0);
        let right = edge_between(&mut model, at(b, 3.0), at(b, 6.0)).unwrap();
        let high = arc(&mut model, 6.0);
        let left = edge_between(&mut model, at(a, 6.0), at(a, 3.0)).unwrap();
        let hole = wire(&mut model, &[low, right, high, left]).unwrap();
        for hole in [hole, hole.reversed()] {
            let face = face_with_holes(&mut model, side, &[hole]).unwrap();
            assert!(crate::validity::check(&model, face).is_valid());
            let area = crate::properties::area(&model, face).mass;
            let expected = 2.0 * std::f64::consts::PI * 5.0 * 10.0 - 5.0 * (b - a) * 3.0;
            assert!((area - expected).abs() < 1e-6, "{area} != {expected}");
        }
    }
}

/// The planar faces `edges` bound (OCCT's `BOPAlgo_Tools::EdgesToWires` then `WiresToFaces`):
/// in each plane the edges span, every region they enclose, holes and all, a compound of
/// them. Unless `shared` (the edges already meet at shared vertices), the edges are first
/// split where they cross or touch. Edges bounding nothing are left out.
#[doc(alias = "BOPAlgo_Tools::EdgesToWires")]
#[doc(alias = "BOPAlgo_Tools::WiresToFaces")]
#[doc(alias = "edgestofaces")]
pub fn faces_from_edges(model: &mut Model, edges: Shape, shared: bool) -> Result<Shape> {
    let input: Vec<Shape> = model.unique_subshapes(edges, NodeKind::Edge).into_iter().filter(|&e| !model.is_degenerate(e)).collect();
    if input.is_empty() {
        return Err(Error::InvalidInput("no edges to make faces of".into()));
    }
    // The edges meeting where they cross: split by one another.
    let pieces: Vec<Shape> = if shared || input.len() == 1 {
        input
    } else {
        let cells = crate::boolean::cells(model, &input)?;
        cells.into_iter().map(|c| c.shape).collect()
    };
    let samples = |model: &Model, e: Shape| -> Vec<Point> {
        model.edge_curve(Shape { orientation: Orientation::Forward, ..e }).map_or_else(Vec::new, |(c, a, b)| (0..=8).map(|i| c.point(a + (b - a) * i as f64 / 8.0)).collect())
    };
    let points: Vec<Vec<Point>> = pieces.iter().map(|&e| samples(model, e)).collect();
    let size = {
        let all: Vec<Point> = points.iter().flatten().copied().collect();
        let b = crate::math::BoundingBox::from_points(all);
        b.min.distance(b.max)
    };
    let tol = pieces.iter().map(|&e| model.tolerance(e)).fold(tolerance::LINEAR, f64::max).max(1e-9 * size);
    // The planes the edges span: of each curved edge lying in one, and of each two edges
    // meeting at an angle at a vertex.
    let mut planes: Vec<Plane> = Vec::new();
    let mut add = |plane: Plane| {
        let n = Vector::from(plane.normal());
        let known = planes.iter().any(|q| Vector::from(q.normal()).cross(n).length() <= 1e-9 && (plane.origin() - q.origin()).dot(Vector::from(q.normal())).abs() <= tol);
        if !known {
            planes.push(plane);
        }
    };
    let plane_through = |ps: &[Point]| -> Option<Plane> {
        let o = *ps.first()?;
        let far = ps.iter().copied().max_by(|a, b| a.distance(o).total_cmp(&b.distance(o)))?;
        let normal = ps.iter().map(|&p| (far - o).cross(p - o)).max_by(|a, b| a.length().total_cmp(&b.length()))?;
        let plane = Plane::new(o, Direction::try_from_vector(normal)?);
        let n = Vector::from(plane.normal());
        ps.iter().all(|&p| (p - o).dot(n).abs() <= tol).then_some(plane)
    };
    for (k, &e) in pieces.iter().enumerate() {
        let straight = matches!(model.edge_curve(e), Some((Curve::Line(_), _, _)));
        if !straight {
            if let Some(p) = plane_through(&points[k]) {
                add(p);
            }
        }
    }
    // (Each edge's direction away from each of its ends.)
    let mut ends: Vec<(Point, Vector)> = Vec::new();
    for &e in &pieces {
        let Some((curve, a, b)) = model.edge_curve(Shape { orientation: Orientation::Forward, ..e }) else { continue };
        let (pa, da) = curve.d1(a);
        let (pb, db) = curve.d1(b);
        ends.push((pa, da));
        ends.push((pb, -db));
    }
    for i in 0..ends.len() {
        for j in i + 1..ends.len() {
            let ((p, d), (q, f)) = (ends[i], ends[j]);
            if p.distance(q) > tol {
                continue;
            }
            if let Some(n) = Direction::try_from_vector(d.cross(f)).filter(|_| d.normalized().cross(f.normalized()).length() > 1e-6) {
                add(Plane::new(p, n));
            }
        }
    }
    // In each plane, the regions its edges enclose.
    let out = model.new_empty(NodeKind::Compound);
    for plane in planes {
        let n = Vector::from(plane.normal());
        let on: Vec<Shape> = (0..pieces.len()).filter(|&k| !points[k].is_empty() && points[k].iter().all(|&p| (p - plane.origin()).dot(n).abs() <= tol)).map(|k| pieces[k]).collect();
        for face in plane_regions(model, &plane, &on, tol)? {
            model.add(out, face)?;
        }
    }
    Ok(out)
}

/// The faces of the regions `edges` (lying in `plane`, meeting at shared vertices) enclose
/// there: the edges' arrangement traced loop by loop, each turning at every vertex to the
/// next edge clockwise; the loops running anticlockwise (in the plane's parameters) bound
/// faces, the clockwise ones are holes of the least face holding them (or the outside).
/// Edges bounding nothing (ending at a vertex of no other) are left out.
fn plane_regions(model: &mut Model, plane: &Plane, edges: &[Shape], tol: f64) -> Result<Vec<Shape>> {
    use crate::math::Point2;
    // The edges by their end vertices, with their curves sampled in the plane's parameters.
    struct Arc {
        edge: Shape,
        ends: [VertexId; 2],
        uv: Vec<Point2>,
    }
    let mut arcs: Vec<Arc> = Vec::new();
    for &e in edges {
        let e = Shape { orientation: Orientation::Forward, ..e };
        let Some((curve, a, b)) = model.edge_curve(e) else { continue };
        let (va, vb) = model.edge_vertices(e);
        let (Ok(va), Ok(vb)) = (VertexId::try_from(va.node), VertexId::try_from(vb.node)) else { continue };
        let n = 64;
        let uv = (0..=n).map(|i| { let (u, v) = plane.parameters(curve.point(a + (b - a) * i as f64 / n as f64)); Point2::new(u, v) }).collect();
        if !arcs.iter().any(|x| x.edge.node == e.node) {
            arcs.push(Arc { edge: e, ends: [va, vb], uv });
        }
    }
    // (Edges ending where no other does bound nothing: taken away, repeatedly.)
    loop {
        let degree = |v: VertexId, arcs: &[Arc]| arcs.iter().map(|a| a.ends.iter().filter(|&&x| x == v).count()).sum::<usize>();
        let before = arcs.len();
        let kept: Vec<bool> = arcs.iter().map(|a| a.ends.iter().all(|&v| degree(v, &arcs) >= 2)).collect();
        let mut k = 0;
        arcs.retain(|_| { k += 1; kept[k - 1] });
        if arcs.len() == before {
            break;
        }
    }
    // Half-edges: (arc, forward); each leaves a vertex at an angle (of its chord to a point
    // a little along it, which also tells apart arcs leaving tangent to each other).
    let leaving = |a: &Arc, forward: bool| -> (VertexId, f64) {
        let m = a.uv.len() - 1;
        let (from, to) = if forward { (a.uv[0], a.uv[1]) } else { (a.uv[m], a.uv[m - 1]) };
        (if forward { a.ends[0] } else { a.ends[1] }, (to.y - from.y).atan2(to.x - from.x))
    };
    let halves: Vec<(usize, bool)> = (0..arcs.len()).flat_map(|i| [(i, true), (i, false)]).collect();
    let mut out_of: crate::hashing::FastMap<VertexId, Vec<(f64, usize)>> = Default::default();
    for (h, &(i, forward)) in halves.iter().enumerate() {
        let (v, angle) = leaving(&arcs[i], forward);
        out_of.entry(v).or_default().push((angle, h));
    }
    for list in out_of.values_mut() {
        list.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    let twin = |h: usize| h ^ 1;
    // After arriving along h at its end, the next half-edge: the one leaving that vertex
    // clockwise next from the way back.
    let next = |h: usize| -> Option<usize> {
        let back = twin(h);
        let (i, forward) = halves[back];
        let (v, _) = leaving(&arcs[i], forward);
        let list = out_of.get(&v)?;
        let k = list.iter().position(|&(_, x)| x == back)?;
        Some(list[(k + list.len() - 1) % list.len()].1)
    };
    let mut used = vec![false; halves.len()];
    let mut loops: Vec<(Vec<usize>, f64, Vec<Point2>)> = Vec::new();
    for start in 0..halves.len() {
        if used[start] {
            continue;
        }
        let mut lp = Vec::new();
        let mut h = start;
        while !used[h] {
            used[h] = true;
            lp.push(h);
            match next(h) {
                Some(x) => h = x,
                None => break,
            }
        }
        if h != start {
            continue;
        }
        let polygon: Vec<Point2> = lp
            .iter()
            .flat_map(|&h| {
                let (i, forward) = halves[h];
                let mut pts = arcs[i].uv.clone();
                if !forward {
                    pts.reverse();
                }
                pts.pop();
                pts
            })
            .collect();
        let n = polygon.len();
        let area = 0.5 * (0..n).map(|k| polygon[k].x * polygon[(k + 1) % n].y - polygon[(k + 1) % n].x * polygon[k].y).sum::<f64>();
        loops.push((lp, area, polygon));
    }
    let inside = |p: Point2, polygon: &[Point2]| {
        let n = polygon.len();
        let mut c = false;
        for k in 0..n {
            let (a, b) = (polygon[k], polygon[(k + 1) % n]);
            if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y) {
                c = !c;
            }
        }
        c
    };
    let small = tol * tol;
    let outer: Vec<usize> = (0..loops.len()).filter(|&k| loops[k].1 > small).collect();
    let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); loops.len()];
    for k in (0..loops.len()).filter(|&k| loops[k].1 < -small) {
        // (A point of the hole's own boundary, off the faces' in general: its middle sample.)
        let polygon = &loops[k].2;
        let p = polygon[polygon.len() / 2];
        // (Not a face's own boundary run the other way: the outside of it.)
        let twins: Vec<usize> = loops[k].0.iter().map(|&h| twin(h)).collect();
        let holder = outer
            .iter()
            .copied()
            .filter(|&o| !loops[o].0.iter().all(|h| twins.contains(h)) && loops[o].1 > -loops[k].1 + small && inside(p, &loops[o].2))
            .min_by(|&a, &b| loops[a].1.total_cmp(&loops[b].1));
        if let Some(o) = holder {
            holes_of[o].push(k);
        }
    }
    let wire_of = |model: &mut Model, lp: &[usize]| -> Result<Shape> {
        let uses: Vec<Shape> = lp
            .iter()
            .map(|&h| {
                let (i, forward) = halves[h];
                Shape { orientation: if forward { Orientation::Forward } else { Orientation::Reversed }, ..arcs[i].edge }
            })
            .collect();
        wire(model, &uses)
    };
    let mut faces = Vec::new();
    for &o in &outer {
        let w = wire_of(model, &loops[o].0)?;
        let mut face = face_on(model, Surface::Plane(*plane), w)?;
        if !holes_of[o].is_empty() {
            let holes: Vec<Shape> = holes_of[o].iter().map(|&k| wire_of(model, &loops[k].0)).collect::<Result<_>>()?;
            face = face_with_holes(model, face, &holes)?;
        }
        faces.push(face);
    }
    Ok(faces)
}
