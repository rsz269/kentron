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
