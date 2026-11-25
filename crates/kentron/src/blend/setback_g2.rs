//! Curvature-continuous fillets set back from the corners where three of them meet.
//!
//! Each fillet is the G2 blend of tangent offset `radius` at the equal-area thumbweight (the
//! order-2 blend section of [`super::fillet_g2`] that takes off what the round of `radius`
//! does), stopped on its section `setback` from the corner; the corner is closed by the convex
//! patch of [`super::g2_corner`], which meets the three fillets and the three faces with their
//! curvature. Built in three steps: each corner's cube of side `setback` is cut away, the
//! fillets are cut along the edges left (their tools run on into the cubes' room, cutting
//! nothing there), and in each cube's room the three walls the cut left are replaced by the
//! patch, the faces' flats grown to its boundary.
//!
//! The patch is fitted as a B-spline (a graph over the plane square to the corner's diagonal,
//! see [`super::g2_corner::Graph`]): its edges are the fillets' end sections and the flats'
//! boundaries exactly, the patch within its edges' tolerance of them.

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{fit, ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Vector};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId, WireId};
use crate::{boolean, tolerance};

use super::g2_corner::{diagonal, onto_plane, CornerSupport, Graph, UnitSection, PLANE};
use super::sections::{fillet_g2_weight, Section};
use super::setback::set_back_corners;
use super::{convex_solid, find_edge, volume_of};

/// `shape` with `edges` rounded by curvature-continuous fillets, as [`super::fillet_setback`]
/// rounds them with circular ones: each the G2 blend of tangent offset `radius` (meeting each
/// face where a ball of `radius` would) at the thumbweight whose section takes off what the
/// ball's round does, except at every corner where three of the edges meet: there each fillet
/// stops `setback` short of the corner, on its section square to its edge, and one convex
/// patch closes the corner, meeting the three fillets and the three faces with their tangent
/// planes and their curvature (see [`super::g2_corner`]). On each face the flat ends along the
/// fillets' section scaled to span from one fillet's contact line to the other's, its apex
/// where those lines would cross.
///
/// `radius` must be positive and `setback` greater than `radius` (at `setback = radius` the
/// patch closes to points on the faces). Supported at corners of three planes meeting square
/// to one another along straight convex edges, each edge long enough to hold its set-back
/// stops, and the fillets otherwise where [`super::fillet_g2`] supports them.
pub fn fillet_setback_g2(model: &mut Model, shape: Shape, edges: &[Shape], radius: f64, setback: f64) -> Result<Shape> {
    if !(radius > tolerance::LINEAR) || !radius.is_finite() {
        return Err(Error::InvalidInput(format!("a fillet radius of {radius}")));
    }
    if !(setback > radius) || !setback.is_finite() {
        return Err(Error::InvalidInput(format!("a setback of {setback} not greater than the radius {radius}")));
    }
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    let section = UnitSection::equal_area_square();
    let corners = set_back_corners(model, shape, edges, radius, setback)?;
    if corners.is_empty() {
        return fillet_g2_weight(model, shape, edges, radius, section.w);
    }
    let support = CornerSupport::new(section.clone(), radius, setback)?;
    // Each corner's cube cut away: a box from the far corner of the cube, `setback` in along
    // each edge, out past the faces.
    let mut cut = shape;
    let before = volume_of(model, shape);
    for &(vertex, along) in &corners {
        let points: Vec<Point> = (0..8).map(|k| (0..3).fold(vertex, |p, i| p + along[i] * if k >> i & 1 == 1 { setback } else { -setback })).collect();
        let tool = convex_solid(model, &points)?.ok_or_else(|| Error::Failed("a set-back corner's cube is flat".into()))?;
        cut = boolean::difference(model, cut, tool)?;
    }
    let removed = before - volume_of(model, cut);
    let expected = corners.len() as f64 * setback.powi(3);
    if (removed - expected).abs() > 1e-6 * expected.max(before) {
        return Err(Error::Failed("a set-back corner's cube cut more than the corner".into()));
    }
    // The named edges as the cubes left them.
    let mut kept = Vec::with_capacity(edges.len());
    for &edge in edges {
        let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("an edge without a curve".into()))?;
        let found = find_edge(model, cut, &curve, first, last).ok_or_else(|| Error::Failed("an edge lost to a set-back corner's cube".into()))?;
        kept.push(found);
    }
    let rounded = fillet_g2_weight(model, cut, &kept, radius, section.w)?;
    // The fillets' faces and the cubes' walls are new; the faces each corner grows are cut by
    // the cubes, new too: none is shared with `shape`, which is left as it was.
    let given: std::collections::HashSet<crate::topology::NodeId> = model.explore(shape, NodeKind::Face).map(|f| f.node).collect();
    let solids: Vec<Shape> = model.explore(rounded, NodeKind::Solid).collect();
    for &(vertex, along) in &corners {
        // (The solid holding the corner's room: the end of its first fillet's section.)
        let mark = vertex + along[0] * setback + along[1] * radius;
        let solid = solids
            .iter()
            .copied()
            .find(|&s| model.explore(s, NodeKind::Vertex).any(|v| model.point(v).distance(mark) <= 1e-6 * setback))
            .ok_or_else(|| Error::Failed("a set-back corner's solid not found".into()))?;
        close_corner(model, solid, vertex, along, &support, edges, &given)?;
    }
    let report = crate::validity::check(model, rounded);
    if !report.is_valid() {
        return Err(Error::Failed(format!("a G2 set-back corner left the solid invalid: {:?}", report.problems.first())));
    }
    Ok(rounded)
}

/// A corner's frame: its vertex and its edges' unit directions.
struct Frame {
    vertex: Point,
    along: [Vector; 3],
}

impl Frame {
    fn world(&self, u: [f64; 3]) -> Point {
        (0..3).fold(self.vertex, |p, k| p + self.along[k] * u[k])
    }

    fn local(&self, p: Point) -> [f64; 3] {
        let d = p - self.vertex;
        self.along.map(|e| d.dot(e))
    }

    fn direction(&self, u: [f64; 3]) -> Vector {
        (0..3).fold(Vector::new(0.0, 0.0, 0.0), |v, k| v + self.along[k] * u[k])
    }
}

/// Replaces, in `solid`, the three walls of the room the cube at `vertex` left (each square to
/// an edge, `setback` along it, between the fillet's end section and the other walls) by the
/// corner's patch, growing each face's flat to the patch's boundary on it.
///
/// The walls leave a hole bounded by the edges each shares with a face that stays: the
/// fillets' end sections (where a section meets its faces so nearly tangentially that the
/// cut found its end a little short, with the short pieces of the faces' traces closing it),
/// and on each face two edges from the flat's corner a, `setback` along one edge and `radius`
/// in from the other, through the cube's corner to b, the other way about. On each face those
/// two give way to the flat's boundary from a to b; the rest of the hole, with those, bounds
/// the patch.
fn close_corner(model: &mut Model, solid: Shape, vertex: Point, along: [Vector; 3], support: &CornerSupport, named: &[Shape], given: &std::collections::HashSet<crate::topology::NodeId>) -> Result<()> {
    let (r, s) = (support.radius, support.setback);
    let frame = Frame { vertex, along };
    let tol = 1e-6 * s;
    let failed = |what: &str| Error::Failed(format!("a G2 set-back corner's {what} not found"));
    let shell = model.explore(solid, NodeKind::Shell).next().ok_or_else(|| failed("shell"))?;
    let faces: Vec<Shape> = model.explore(shell, NodeKind::Face).collect();
    // The plane a face lies on, as (axis, offset along it), where it is square to an edge.
    let plane_of = |model: &Model, f: Shape| -> Option<(usize, f64)> {
        let Some(Surface::Plane(plane)) = model.surface(f) else { return None };
        let n = Vector::from(plane.normal());
        let k = (0..3).find(|&k| n.cross(along[k]).length() <= 1e-9)?;
        Some((k, frame.local(plane.origin())[k]))
    };
    let in_room = |model: &Model, f: Shape| model.explore(f, NodeKind::Vertex).all(|v| frame.local(model.point(v)).iter().all(|&x| x > -tol && x < s + tol));
    let walls: Vec<Shape> = faces.iter().copied().filter(|&f| plane_of(model, f).is_some_and(|(_, d)| (d - s).abs() <= tol) && in_room(model, f)).collect();
    if walls.len() != 3 {
        return Err(failed("room's walls"));
    }
    // The hole's edges, each with the face beside it that stays.
    let mut uses: Vec<(Shape, usize)> = Vec::new();
    for &wall in &walls {
        for e in model.explore(wall, NodeKind::Edge) {
            match uses.iter_mut().find(|(x, _)| x.node == e.node) {
                Some(u) => u.1 += 1,
                None => uses.push((Shape { orientation: Orientation::Forward, ..e }, 1)),
            }
        }
    }
    let hole: Vec<(Shape, Shape)> = uses
        .into_iter()
        .filter(|&(_, n)| n == 1)
        .map(|(e, _)| {
            let beside = faces.iter().copied().find(|f| !walls.iter().any(|w| w.node == f.node) && model.explore(*f, NodeKind::Edge).any(|x| x.node == e.node));
            beside.map(|f| (e, f)).ok_or_else(|| failed("face beside the room"))
        })
        .collect::<Result<_>>()?;
    let ends = |model: &Model, e: Shape| {
        let (p, q) = model.edge_vertices(e);
        (p, q)
    };
    // On face i, the hole's edges beside it run from a_i through the cube's corner to b_i.
    let point = |i: usize, a: f64, b: f64| {
        let mut u = [0.0; 3];
        u[(i + 1) % 3] = a;
        u[(i + 2) % 3] = b;
        u
    };
    let mut flats: Vec<(EdgeId, Shape, Shape, Vec<Shape>)> = Vec::new();
    for i in 0..3 {
        let chain: Vec<Shape> = hole.iter().filter(|(_, f)| plane_of(model, *f).is_some_and(|(k, d)| k == i && d.abs() <= tol)).map(|(e, _)| *e).collect();
        let (a, b) = (point(i, s, r), point(i, r, s));
        let find = |u: [f64; 3]| -> Result<Shape> {
            chain.iter().flat_map(|&e| {
                let (p, q) = ends(model, e);
                [p, q]
            }).find(|&v| model.point(v).distance(frame.world(u)) <= tol).ok_or_else(|| failed("flat's corner"))
        };
        let (va, vb) = (find(a)?, find(b)?);
        let face = hole.iter().find(|(e, _)| e.node == chain[0].node).map(|(_, f)| *f).expect("the chain's face");
        if given.contains(&face.node) {
            return Err(failed("face of its own (one shared with the solid given)"));
        }
        let (j, l) = ((i + 1) % 3, (i + 2) % 3);
        let mut apex = [0.0; 3];
        apex[j] = r;
        apex[l] = r;
        let curve = Section::g2(frame.world(a), frame.world(apex), frame.world(b), support.section.w).curve()?;
        let edge = model.add_edge(VertexId::try_from(va.node).expect("a vertex"), VertexId::try_from(vb.node).expect("a vertex"), tolerance::LINEAR)?;
        model.set_edge_curve(Shape::from(edge), Curve::BSpline(curve), tolerance::LINEAR)?;
        model.set_edge_range(Shape::from(edge), 0.0, 1.0)?;
        // In the face's wire the chain's coedges, one after another (round the wire's end
        // perhaps), give way to the new edge's, run the same way.
        let wire = model.explore(face, NodeKind::Wire).find(|w| model.explore(*w, NodeKind::Edge).any(|e| e.node == chain[0].node)).ok_or_else(|| failed("face's wire"))?;
        let wire_id = WireId::try_from(wire.node).expect("a wire");
        let coedges = model[wire_id].coedges.clone();
        let member: Vec<bool> = coedges.iter().map(|&c| chain.iter().any(|e| e.node == Shape::from(model[c].edge).node)).collect();
        let n = coedges.len();
        let first = (0..n).find(|&k| member[k] && !member[(k + n - 1) % n]).ok_or_else(|| failed("chain in the face's wire"))?;
        let starts = {
            let c = coedges[first];
            let (p, q) = ends(model, Shape::from(model[c].edge));
            if model[c].orientation == Orientation::Forward { p } else { q }
        };
        let forward = starts.node == va.node;
        let new = model.add_coedge(edge, if forward { Orientation::Forward } else { Orientation::Reversed })?;
        let mut order: Vec<CoEdgeId> = Vec::with_capacity(n);
        for k in 0..n {
            let at = (first + k) % n;
            if k == 0 {
                order.push(new);
            } else if !member[at] {
                order.push(coedges[at]);
            }
        }
        if order.len() + chain.len() != n + 1 || (0..chain.len()).any(|k| !member[(first + k) % n]) {
            return Err(failed("chain run together in the face's wire"));
        }
        for (k, &c) in coedges.iter().enumerate() {
            if member[k] {
                model.remove_coedge(wire_id, c)?;
            }
        }
        model.append_coedge(wire_id, new)?;
        model.set_coedge_order(wire_id, &order)?;
        flats.push((edge, va, vb, chain));
    }
    // The walls go (and with them the edges and vertices only they held).
    for &wall in &walls {
        model.remove(shell, wall)?;
        model.remove_subtree(wall.node);
    }
    // The patch's loop: the hole's edges round from a_0, each face's chain given way to its
    // flat's edge.
    let mut rest: Vec<Shape> = hole.iter().map(|(e, _)| *e).filter(|e| !flats.iter().any(|f| f.3.iter().any(|c| c.node == e.node))).collect();
    let mut loop_: Vec<(Shape, Orientation)> = Vec::new();
    let start = flats[0].1;
    let mut at = start;
    for _ in 0..rest.len() + 3 {
        if let Some(f) = flats.iter().find(|f| f.1.node == at.node || f.2.node == at.node) {
            if !loop_.iter().any(|(e, _)| e.node == Shape::from(f.0).node) {
                let forward = f.1.node == at.node;
                loop_.push((Shape::from(f.0), if forward { Orientation::Forward } else { Orientation::Reversed }));
                at = if forward { f.2 } else { f.1 };
                if at.node == start.node {
                    break;
                }
                continue;
            }
        }
        let k = rest.iter().position(|&e| {
            let (p, q) = ends(model, e);
            p.node == at.node || q.node == at.node
        }).ok_or_else(|| failed("patch's loop"))?;
        let e = rest.remove(k);
        let (p, q) = ends(model, e);
        let forward = p.node == at.node;
        loop_.push((e, if forward { Orientation::Forward } else { Orientation::Reversed }));
        at = if forward { q } else { p };
        if at.node == start.node {
            break;
        }
    }
    if at.node != start.node || !rest.is_empty() {
        return Err(failed("patch's closed loop"));
    }
    // The patch: the graph of its height, fitted finer the larger the corner is against the
    // radius (24 spans to the radius), each way over its boundary's box.
    let width = {
        let b = support.boundary(64);
        let xs: Vec<[f64; 2]> = b.iter().flatten().map(|u| onto_plane(*u)).collect();
        let span = |k: usize| xs.iter().map(|x| x[k]).fold(f64::NEG_INFINITY, f64::max) - xs.iter().map(|x| x[k]).fold(f64::INFINITY, f64::min);
        span(0).max(span(1))
    };
    let spans = ((24.0 * width / r).ceil() as usize).clamp(24, 128);
    let graph = support.graph(spans)?;
    let surface = graph_surface(&frame, &graph)?;
    // Each edge's pcurve its projection onto the plane.
    let mut uses: Vec<(EdgeId, Orientation, Curve2, f64, f64)> = Vec::new();
    let mut worst: f64 = 0.0;
    for (e, o) in loop_ {
        let (curve, first, last) = model.edge_curve(e).ok_or_else(|| failed("patch edge's curve"))?;
        let (pcurve, deviation) = on_graph(&frame, &graph, &curve, first, last)?;
        worst = worst.max(deviation);
        uses.push((EdgeId::try_from(e.node).expect("an edge"), o, pcurve, first, last));
    }
    // The loop as it runs, about the plane's normal: counter-clockwise in the surface's
    // parameters (its own), else reversed; the face turned to face out, away from the solid,
    // which lies along the diagonal beyond the patch.
    let polygon: Vec<Point2> = uses
        .iter()
        .flat_map(|(_, o, pcurve, first, last)| {
            (0..16).map(move |k| {
                let t = k as f64 / 16.0;
                let t = if *o == Orientation::Forward { first + (last - first) * t } else { last - (last - first) * t };
                pcurve.point(t)
            })
        })
        .collect();
    let area: f64 = (0..polygon.len()).map(|k| {
        let (p, q) = (polygon[k], polygon[(k + 1) % polygon.len()]);
        p.x * q.y - q.x * p.y
    }).sum();
    if area < 0.0 {
        uses.reverse();
        for u in &mut uses {
            u.1 = u.1.reversed();
        }
    }
    let (_, du, dv) = surface.d1(0.5 * (graph.lo[0] + graph.hi[0]), 0.5 * (graph.lo[1] + graph.hi[1]));
    let inward = frame.direction([1.0, 1.0, 1.0]);
    let outward_forward = du.cross(dv).dot(inward) < 0.0;
    let tolerance = (1.5 * worst.max(graph.deviation)).max(tolerance::LINEAR);
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, tolerance)?;
    let mut coedges = Vec::with_capacity(uses.len());
    for (edge, o, pcurve, first, last) in uses {
        if model[edge].tolerance < tolerance {
            model.set_edge_tolerance(edge, tolerance)?;
        }
        let c = model.add_coedge(edge, o)?;
        model.set_pcurve_with_range(c, pcurve, first, last)?;
        coedges.push(c);
    }
    let w = model.add_wire(&coedges)?;
    model.add(face, Shape::from(w))?;
    for v in model.unique_subshapes(face, NodeKind::Vertex) {
        let id = VertexId::try_from(v.node).expect("a vertex");
        if model[id].tolerance < tolerance {
            model.set_vertex_tolerance(id, tolerance)?;
        }
    }
    model.add(shell, if outward_forward { face } else { face.reversed() })?;
    // (The patch comes of the corner's named edges.)
    for &edge in named {
        let (p, q) = model.edge_vertices(edge);
        if model.point(p).distance(vertex) <= tol || model.point(q).distance(vertex) <= tol {
            boolean::record(model.uid(), edge.node, face.node);
        }
    }
    Ok(())
}

/// The graph as a B-spline surface of `frame`'s space: its poles over the knots' Greville
/// points, `corner + x·PLANE + h·diagonal`.
fn graph_surface(frame: &Frame, graph: &Graph) -> Result<Surface> {
    let (gu, gv) = (graph.greville(0), graph.greville(1));
    let e = [frame.direction(PLANE[0]), frame.direction(PLANE[1])];
    let up = frame.direction([1.0, 1.0, 1.0]) * (1.0 / 3f64.sqrt());
    let poles: Vec<Vec<Point>> = gu.iter().enumerate().map(|(i, &x)| gv.iter().enumerate().map(|(j, &y)| frame.vertex + e[0] * x + e[1] * y + up * graph.heights[i][j]).collect()).collect();
    let knots = |k: usize| {
        let (knots, mults) = graph.knots(k);
        KnotVector::new(knots, mults, Graph::DEGREE)
    };
    Ok(Surface::BSpline(BSplineSurface::try_new(poles, knots(0), knots(1))?))
}

/// The pcurve on the graph's surface of an edge's `curve` from `first` to `last` (its
/// projection onto the plane, exactly so for a B-spline: the projection is linear), and how
/// far the curve strays from the surface.
fn on_graph(frame: &Frame, graph: &Graph, curve: &Curve, first: f64, last: f64) -> Result<(Curve2, f64)> {
    let flat = |p: Point| {
        let x = onto_plane(frame.local(p));
        Point2::new(x[0], x[1])
    };
    let pcurve = match curve {
        Curve::BSpline(b) => {
            let poles: Vec<Point2> = b.poles().iter().map(|&p| flat(p)).collect();
            let mut c = crate::geometry::curve2::BSplineCurve2::try_new(poles, b.knots().to_vec(), b.multiplicities().to_vec(), b.degree())?;
            if let Some(w) = b.weights() {
                c = c.with_weights(w.to_vec());
            }
            Curve2::BSpline(c)
        }
        _ => {
            let k = 64;
            let ts: Vec<f64> = (0..=k).map(|i| first + (last - first) * i as f64 / k as f64).collect();
            let points: Vec<Point2> = ts.iter().map(|&t| flat(curve.point(t))).collect();
            Curve2::BSpline(fit::interpolate_with_parameters(&points, &ts)?)
        }
    };
    let k = 128;
    let deviation = (0..=k)
        .map(|i| {
            let t = first + (last - first) * i as f64 / k as f64;
            let u = frame.local(curve.point(t));
            (graph.height(onto_plane(u)) - diagonal(u)).abs()
        })
        .fold(0.0, f64::max);
    Ok((pcurve, deviation))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::{classify_point, PointLocation};
    use crate::properties::{curvatures, normals};

    fn block(model: &mut Model, w: [f64; 3]) -> Shape {
        crate::primitives::Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(w[0], w[1], w[2])).build(model).unwrap()
    }

    fn corner_edges(model: &Model, block: Shape, at: Point) -> Vec<Shape> {
        model
            .unique_subshapes(block, NodeKind::Edge)
            .into_iter()
            .filter(|&e| {
                let (a, b) = model.edge_vertices(e);
                model.point(a).distance(at) < 1e-9 || model.point(b).distance(at) < 1e-9
            })
            .collect()
    }

    /// The box `w` with its corner at `w` set back: the model, the result and the corner's
    /// support function.
    fn one_corner(w: [f64; 3], r: f64, s: f64) -> (Model, Shape, CornerSupport) {
        let mut model = Model::new();
        let b = block(&mut model, w);
        let edges = corner_edges(&model, b, Point::new(w[0], w[1], w[2]));
        assert_eq!(edges.len(), 3);
        let out = fillet_setback_g2(&mut model, b, &edges, r, s).unwrap();
        (model, out, CornerSupport::new(UnitSection::equal_area_square(), r, s).unwrap())
    }

    /// The corner's frame at `w`: its edges run back along the axes (x, then y, then z).
    fn frame_at(w: [f64; 3]) -> Frame {
        Frame { vertex: Point::new(w[0], w[1], w[2]), along: [Vector::new(-1.0, 0.0, 0.0), Vector::new(0.0, -1.0, 0.0), Vector::new(0.0, 0.0, -1.0)] }
    }

    /// The volume the support function gives: the box less the corner's cube, plus what of
    /// the cube the solid holds ([`super::super::g2_corner::volume_in_cube`]), less each
    /// fillet's section beyond the cube along its edge.
    fn reference(w: [f64; 3], c: &CornerSupport) -> f64 {
        let (r, s) = (c.radius, c.setback);
        let section = 0.5 * c.section.area_share() * r * r;
        w[0] * w[1] * w[2] - s.powi(3) + super::super::g2_corner::volume_in_cube(c, 32) - section * w.iter().map(|x| x - s).sum::<f64>()
    }

    /// Volume, validity and the point tests at one corner, for set-backs from just past the
    /// radius to six radii: the solid's volume is the support function's within the fitted
    /// patch's tolerance; where the fillets' contact lines would cross on each face the
    /// solid holds nothing (the flat stops short), the flat holds points just beyond its
    /// boundary, and past the set-back each edge is its fillet.
    #[test]
    fn corners_hold_their_support_function_s_volume() {
        let w = [12.0, 14.0, 16.0];
        for (r, s) in [(1.0, 1.08), (1.0, 1.6), (1.3, 3.6), (0.8, 4.8)] {
            let (model, out, c) = one_corner(w, r, s);
            assert!(crate::validity::check(&model, out).is_valid());
            let got = volume_of(&model, out);
            let want = reference(w, &c);
            assert!((got - want).abs() < 2e-4 * r.powi(3), "{r} {s}: volume {got} against {want}");
            let frame = frame_at(w);
            let inside = |u: [f64; 3]| matches!(classify_point(&model, out, frame.world(u)).unwrap(), PointLocation::Inside);
            for i in 0..3 {
                let (j, l) = ((i + 1) % 3, (i + 2) % 3);
                let mut away = [0.0; 3];
                away[i] = 1e-5 * r;
                away[j] = r;
                away[l] = r;
                assert!(!inside(away), "{r} {s}: {away:?}");
                // Just past the flat's boundary's middle, on the face's side away from the
                // corner, the flat.
                let m = c.section.point(0.5);
                let mut flat = away;
                flat[j] = r + (s - r) * m[0] + 1e-3 * r;
                flat[l] = r + (s - r) * m[1] + 1e-3 * r;
                assert!(inside(flat), "{r} {s}: {flat:?}");
                // A little past the set-back, across the fillet's section's middle.
                let mut fillet = [0.0; 3];
                fillet[i] = s + 0.05 * r;
                fillet[j] = r * m[0];
                fillet[l] = r * m[1];
                let (mut a, mut b) = (fillet, fillet);
                a[j] += 1e-3 * r;
                a[l] += 1e-3 * r;
                b[j] -= 1e-3 * r;
                b[l] -= 1e-3 * r;
                assert!(inside(a) && !inside(b), "{r} {s}: {fillet:?}");
            }
        }
    }

    /// The built patch bends as the patch it is fitted to: its normals and curvatures on the
    /// solid against the support function's, near each fillet's seam and each flat's boundary
    /// and across the middle; and just past the seam the fillet bends as its section does.
    /// (The support function's own curvature tends to the fillets' and the faces', see
    /// `the_patch_s_curvature_tends_to_the_fillets_and_the_faces`.)
    #[test]
    fn the_built_patch_bends_as_its_support_function_does() {
        let w = [12.0, 14.0, 16.0];
        let (r, s) = (1.0, 2.9);
        let (mut model, out, c) = one_corner(w, r, s);
        let frame = frame_at(w);
        let sec = &c.section;
        let mut normals_at = Vec::new();
        for t in [0.2, 0.35, 0.5, 0.7] {
            let (p, q) = (sec.point(t - 1e-5), sec.point(t + 1e-5));
            let tangent = [(q[0] - p[0]) / 2e-5, (q[1] - p[1]) / 2e-5];
            let normal = [-tangent[1], tangent[0]];
            let normal = if normal[0] < 0.0 { [-normal[0], -normal[1]] } else { normal };
            for i in 0..3 {
                let (j, l) = ((i + 1) % 3, (i + 2) % 3);
                for e in [1e-3, 1e-2, 0.3] {
                    let mut n = [0.0; 3];
                    n[i] = e;
                    n[j] = normal[0];
                    n[l] = normal[1];
                    normals_at.push(n);
                    let mut m = [0.0; 3];
                    m[i] = 1.0;
                    m[j] = e * normal[0];
                    m[l] = e * normal[1];
                    normals_at.push(m);
                }
                // Just past fillet i's seam, across its section at t.
                let p = sec.point(t);
                let mut u = [0.0; 3];
                u[i] = s + 1e-3 * s;
                u[j] = r * p[0];
                u[l] = r * p[1];
                let k = curvatures(&mut model, out, &[frame.world(u)]).unwrap()[0].unwrap();
                let want = c.curvatures([1e-12, normal[0], normal[1]]);
                assert!((k.max - want[0]).abs() < 1e-3 && k.min.abs() < 1e-9, "{t} {i}: {k:?} {want:?}");
            }
        }
        let mut worst = (0.0f64, 0.0f64);
        for n in normals_at {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let outward = frame.direction(n.map(|v| -v / len));
            // (On the solid where a ray along the normal meets it: the fitted patch lies within
            // its tolerance of the point, farther than a normal's query reaches.)
            let near = frame.world(c.at(n).point);
            let from = near + outward * (0.01 * r);
            let hit = crate::raycast::raycast(&model, out, from, -outward).unwrap().into_iter().filter(|&d| d > 0.0).fold(f64::INFINITY, f64::min);
            assert!((hit - 0.01 * r).abs() < 1e-4, "{hit}");
            let point = from - outward * hit;
            let k = curvatures(&mut model, out, &[point]).unwrap()[0].unwrap();
            let want = c.curvatures(n);
            worst.0 = worst.0.max((k.max - want[0]).abs()).max((k.min - want[1]).abs());
            let got = normals(&mut model, out, &[point]).unwrap()[0].unwrap();
            worst.1 = worst.1.max(Vector::from(got).cross(outward).length());
        }
        assert!(worst.0 < 0.01 / r && worst.1 < 1e-3, "{worst:?}");
    }

    /// The corner of a cube's three equal edges is the same however they are named: points
    /// about the patch classify as their images under each swap of the axes do (those
    /// within the fitted patch's tolerance of it left out).
    #[test]
    fn a_square_corner_is_the_same_under_every_swap_of_its_edges() {
        let w = [10.0, 10.0, 10.0];
        let (r, s) = (1.2, 3.3);
        let (model, out, c) = one_corner(w, r, s);
        let frame = frame_at(w);
        let inside = |u: [f64; 3]| matches!(classify_point(&model, out, frame.world(u)).unwrap(), PointLocation::Inside);
        let mut seed = 7u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut compared = 0;
        for _ in 0..400 {
            let u = [s * next(), s * next(), s * next()];
            // (Not near the patch: its height over the point's foot more than 1e-4 away.)
            if let Some((h, _)) = c.height(onto_plane(u)) {
                if (h - diagonal(u)).abs() < 1e-4 {
                    continue;
                }
            }
            let at = inside(u);
            for perm in [[1, 2, 0], [2, 0, 1], [1, 0, 2], [0, 2, 1], [2, 1, 0]] {
                assert_eq!(inside([u[perm[0]], u[perm[1]], u[perm[2]]]), at, "{u:?} {perm:?}");
            }
            compared += 1;
        }
        assert!(compared > 300);
    }

    /// All twelve edges of a box set back at all eight corners: every corner as one alone,
    /// the fillets between corners cut at both ends.
    #[test]
    fn a_whole_box_sets_back_every_corner() {
        let w = [9.0, 10.0, 11.0];
        let (r, s) = (1.0, 2.5);
        let mut model = Model::new();
        let b = block(&mut model, w);
        let edges = model.unique_subshapes(b, NodeKind::Edge);
        assert_eq!(edges.len(), 12);
        let out = fillet_setback_g2(&mut model, b, &edges, r, s).unwrap();
        assert!(crate::validity::check(&model, out).is_valid());
        let c = CornerSupport::new(UnitSection::equal_area_square(), r, s).unwrap();
        let section = 0.5 * c.section.area_share() * r * r;
        let corner = s.powi(3) - super::super::g2_corner::volume_in_cube(&c, 32);
        let want = w[0] * w[1] * w[2] - 8.0 * corner - 4.0 * section * w.iter().map(|x| x - 2.0 * s).sum::<f64>();
        let got = volume_of(&model, out);
        assert!((got - want).abs() < 1e-3, "{got} against {want}");
    }

    /// Refused: a set-back not past the radius, a radius not positive, edges too short for
    /// their stops at both ends, and a corner whose faces are not square to one another.
    #[test]
    fn set_backs_it_cannot_make_are_refused() {
        let mut model = Model::new();
        let b = block(&mut model, [10.0, 10.0, 10.0]);
        let edges = model.unique_subshapes(b, NodeKind::Edge);
        assert!(matches!(fillet_setback_g2(&mut model, b, &edges, 1.0, 1.0), Err(Error::InvalidInput(_))));
        assert!(matches!(fillet_setback_g2(&mut model, b, &edges, 0.0, 2.0), Err(Error::InvalidInput(_))));
        assert!(matches!(fillet_setback_g2(&mut model, b, &edges, 1.0, 5.5), Err(Error::InvalidInput(_))));
        // A wedge's corner: one face slanted.
        let points = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(0.0, 0.0, 10.0),
            Point::new(10.0, 0.0, 6.0),
            Point::new(0.0, 10.0, 10.0),
            Point::new(10.0, 10.0, 6.0),
        ];
        let wedge = convex_solid(&mut model, &points).unwrap().unwrap();
        let at = Point::new(10.0, 10.0, 6.0);
        let edges = corner_edges(&model, wedge, at);
        assert_eq!(edges.len(), 3);
        assert!(fillet_setback_g2(&mut model, wedge, &edges, 1.0, 2.5).is_err());
    }
}
