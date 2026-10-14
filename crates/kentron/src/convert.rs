//! Converting the geometry of shapes.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point2, Transform, Transformable};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

/// A copy of `shape` whose curves and surfaces are all B-splines. Placements are applied to
/// the geometry (the copy is unplaced). Curves and surfaces convert exactly (conics and
/// surfaces of revolution as rational B-splines; surfaces over their faces' parameter
/// ranges), lines, planes and B-splines keeping their parametrization. Where a conversion
/// changes it (along a circle), pcurves are fitted anew on the converted surface and edge
/// tolerances raised to cover their deviation. Offset curves and surfaces are not
/// supported.
#[doc(alias = "BRepBuilderAPI_NurbsConvert")]
pub fn to_nurbs(model: &mut Model, shape: Shape) -> Result<Shape> {
    let mut copier = Copier { vertices: HashMap::new(), edges: HashMap::new() };
    let copy = copier.copy(model, Shape { orientation: Orientation::Forward, ..shape })?;
    Ok(Shape { orientation: shape.orientation, ..copy })
}

/// A degree-1 B-spline through `a` and `b` over `[first, last]`.
fn segment<P: crate::geometry::ControlPoint>(a: P, b: P, first: f64, last: f64) -> Result<BSplineCurve<P>> {
    BSplineCurve::try_new(vec![a, b], vec![first, last], vec![2, 2], 1)
}

fn curve_to_nurbs(curve: &Curve, first: f64, last: f64) -> Result<Curve> {
    let trimmed = Curve::Trimmed(crate::geometry::curve::TrimmedCurve::new(curve.clone(), first, last));
    Ok(Curve::BSpline(match curve {
        Curve::Line(_) => segment(curve.point(first), curve.point(last), first, last)?,
        Curve::BSpline(b) => b.clone(),
        other => crate::geometry::convert::curve_to_bspline(&trimmed)
            .ok_or_else(|| Error::InvalidInput(format!("NURBS conversion of a {} is not supported", kind(other))))?,
    }))
}

fn pcurve_to_nurbs(curve: &Curve2, first: f64, last: f64) -> Result<Curve2> {
    let trimmed = Curve2::Trimmed(crate::geometry::curve2::TrimmedCurve2::new(curve.clone(), first, last));
    Ok(Curve2::BSpline(match curve {
        Curve2::Line(_) => {
            let (a, b) = (curve.point(first), curve.point(last));
            segment(Point2::new(a.x, a.y), Point2::new(b.x, b.y), first, last)?
        }
        Curve2::BSpline(b) => b.clone(),
        other => crate::geometry::convert::curve2_to_bspline(&trimmed)
            .ok_or_else(|| Error::InvalidInput(format!("NURBS conversion of a 2D {} is not supported", kind(other))))?,
    }))
}

/// Whether converting keeps the curve's parametrization (else it changes between knots, as
/// a circle's).
fn keeps_parameters(curve: &Curve) -> bool {
    match curve {
        Curve::Line(_) | Curve::Bezier(_) | Curve::BSpline(_) => true,
        Curve::Trimmed(t) => keeps_parameters(t.basis()),
        _ => false,
    }
}

fn keeps_parameters2(curve: &Curve2) -> bool {
    match curve {
        Curve2::Line(_) | Curve2::Bezier(_) | Curve2::BSpline(_) => true,
        Curve2::Trimmed(t) => keeps_parameters2(t.basis()),
        _ => false,
    }
}

/// Whether a circle or ellipse (converted alike: the same parametrization over the same
/// range).
fn elliptic(curve: &Curve) -> bool {
    matches!(curve, Curve::Circle(_) | Curve::Ellipse(_))
}

fn elliptic2(curve: &Curve2) -> bool {
    matches!(curve, Curve2::Circle(_) | Curve2::Ellipse(_))
}

/// The pcurve of an edge on the converted surface `new_surface` (of `surface`): the old
/// one converted where its parametrization follows the edge's, else fitted through the
/// points of the converted edge curve `new_curve` on the new surface (or, on a degenerate
/// edge, a segment between its ends mapped onto it). Returns it with the largest distance
/// between the edge curve and the surface along it.
#[allow(clippy::too_many_arguments)]
fn new_pcurve(
    curve: Option<&Curve>,
    new_curve: Option<&Curve>,
    pcurve: &Curve2,
    first: f64,
    last: f64,
    surface: &Surface,
    new_surface: &Surface,
    surface_kept: bool,
) -> Result<(Curve2, f64)> {
    let follows = match curve {
        None => keeps_parameters2(pcurve),
        Some(c) => (keeps_parameters(c) && keeps_parameters2(pcurve)) || (elliptic(c) && elliptic2(pcurve)),
    };
    if surface_kept && follows {
        return Ok((pcurve_to_nurbs(pcurve, first, last)?, 0.0));
    }
    let (u0, u1, v0, v1) = new_surface.bounds();
    let inside = |p: Point2| (p.x.clamp(u0, u1), p.y.clamp(v0, v1));
    let onto = |q: crate::math::Point, guess: Point2| {
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(q, new_surface, Some(inside(guess)));
        Point2::new(u, v)
    };
    let Some(new_curve) = new_curve else {
        // A degenerate edge (at a pole): its ends, each parameter found along a regular
        // line of the other.
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let map = |p: Point2| Point2::new(onto(surface.point(p.x, vm), Point2::new(p.x, vm)).x, onto(surface.point(um, p.y), Point2::new(um, p.y)).y);
        let (a, b) = (map(pcurve.point(first)), map(pcurve.point(last)));
        return Ok((Curve2::BSpline(segment(a, b, first, last)?), 0.0));
    };
    let spans = match new_curve {
        Curve::BSpline(b) => b.knots().len(),
        _ => 1,
    };
    let n = (64 + 8 * spans).min(1024);
    let at = |i: usize| first + (last - first) * i as f64 / n as f64;
    let params: Vec<f64> = (0..=n).map(at).collect();
    let mut points = Vec::with_capacity(n + 1);
    for &t in &params {
        // (Converting moves points along the curves and surfaces only a little: the old
        // pcurve's point is close, and on the right side of a seam.)
        points.push(onto(new_curve.point(t), pcurve.point(t)));
    }
    let fitted = Curve2::BSpline(crate::geometry::fit::interpolate_with_parameters(&points, &params)?);
    let mut deviation = 0.0f64;
    for i in 0..n {
        let t = 0.5 * (at(i) + at(i + 1));
        let p = fitted.point(t);
        deviation = deviation.max(new_surface.point(p.x, p.y).distance(new_curve.point(t)));
    }
    Ok((fitted, deviation))
}

fn kind(x: &impl std::fmt::Debug) -> String {
    format!("{x:?}").split('(').next().unwrap_or_default().to_lowercase()
}

/// The surface as a B-spline over the parameter rectangle `[u1, u2] x [v1, v2]`.
fn surface_to_nurbs(surface: &Surface, (u1, u2, v1, v2): (f64, f64, f64, f64)) -> Result<Surface> {
    Ok(Surface::BSpline(match surface {
        Surface::Plane(_) => {
            let poles = vec![
                vec![surface.point(u1, v1), surface.point(u1, v2)],
                vec![surface.point(u2, v1), surface.point(u2, v2)],
            ];
            BSplineSurface::try_new(poles, KnotVector::new(vec![u1, u2], vec![2, 2], 1), KnotVector::new(vec![v1, v2], vec![2, 2], 1))?
        }
        Surface::BSpline(b) => b.clone(),
        other => {
            let trimmed = Surface::Trimmed(crate::geometry::surface::TrimmedSurface::new(other.clone(), u1, u2, v1, v2));
            crate::geometry::convert::surface_to_bspline(&trimmed)
                .ok_or_else(|| Error::InvalidInput(format!("NURBS conversion of a {} surface is not supported", kind(other))))?
        }
    }))
}

/// Copies topology with converted geometry, sharing what the original shares.
struct Copier {
    vertices: HashMap<NodeId, Vec<(Transform, VertexId)>>,
    edges: HashMap<NodeId, Vec<(Transform, EdgeId)>>,
}

impl Copier {
    /// A copy of the node of `shape` (forward, its placement applied to the geometry).
    fn copy(&mut self, model: &mut Model, shape: Shape) -> Result<Shape> {
        match shape.kind() {
            NodeKind::Vertex => Ok(self.vertex(model, shape).into()),
            NodeKind::Edge => Ok(self.edge(model, shape)?.into()),
            NodeKind::Face => self.face(model, shape),
            NodeKind::Wire => Err(Error::InvalidInput("converting a wire outside a face is not supported".into())),
            NodeKind::Shell | NodeKind::Solid | NodeKind::CompSolid | NodeKind::Compound => {
                let copy = model.new_empty(shape.kind());
                let children: Vec<Shape> = model.children_local(shape).collect();
                for child in children {
                    let placed = Shape { orientation: Orientation::Forward, placement: shape.placement * child.placement, ..child };
                    let c = self.copy(model, placed)?;
                    model.add(copy, Shape { orientation: child.orientation, ..c })?;
                }
                Ok(copy)
            }
            other => Err(Error::InvalidInput(format!("cannot convert a {other:?}"))),
        }
    }

    fn vertex(&mut self, model: &mut Model, shape: Shape) -> VertexId {
        if let Some(&(_, v)) = self.vertices.get(&shape.node).and_then(|l| l.iter().find(|(t, _)| *t == shape.placement)) {
            return v;
        }
        let tolerance = model.tolerance(shape);
        let v = model.add_vertex(model.point(shape), tolerance);
        self.vertices.entry(shape.node).or_default().push((shape.placement, v));
        v
    }

    fn edge(&mut self, model: &mut Model, shape: Shape) -> Result<EdgeId> {
        if let Some(&(_, e)) = self.edges.get(&shape.node).and_then(|l| l.iter().find(|(t, _)| *t == shape.placement)) {
            return Ok(e);
        }
        let forward = Shape { orientation: Orientation::Forward, ..shape };
        let (a, b) = model.edge_vertices(forward);
        let (a, b) = (self.vertex(model, a), self.vertex(model, b));
        let tolerance = model.tolerance(shape);
        let e = model.add_edge(a, b, tolerance)?;
        if let Some((curve, first, last)) = model.edge_curve(forward) {
            model.set_edge_curve(e.into(), curve_to_nurbs(&curve, first, last)?, tolerance)?;
            model.set_edge_range(e.into(), first, last)?;
        }
        self.edges.entry(shape.node).or_default().push((shape.placement, e));
        Ok(e)
    }

    fn face(&mut self, model: &mut Model, shape: Shape) -> Result<Shape> {
        let surface = model
            .surface(shape)
            .ok_or_else(|| Error::InvalidInput("converting a face without a surface".into()))?
            .transformed(&shape.placement);
        let tolerance = model.tolerance(shape);
        let copy = model.new_empty(NodeKind::Face);
        let wires: Vec<Shape> = model.children_local(shape).collect();
        let mut new_wires = Vec::new();
        let (mut u1, mut u2, mut v1, mut v2) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for wire in wires {
            let mut coedges = Vec::new();
            let coedge_shapes: Vec<Shape> = model.children_local(wire).collect();
            for c in coedge_shapes {
                let id = CoEdgeId::try_from(c.node).expect("a wire holds coedges");
                let (edge, orientation) = (model[id].edge, model[id].orientation);
                let placed = Shape { placement: shape.placement * wire.placement * c.placement, ..Shape::from(edge) };
                let old_curve = model.edge_curve(placed).map(|(c, _, _)| c);
                let new_edge = self.edge(model, placed)?;
                let (pcurve, first, last) = model
                    .effective_pcurve(id)
                    .ok_or_else(|| Error::InvalidInput("converting a face without pcurves".into()))?;
                for i in 0..=64 {
                    let p = pcurve.point(first + (last - first) * i as f64 / 64.0);
                    (u1, u2, v1, v2) = (u1.min(p.x), u2.max(p.x), v1.min(p.y), v2.max(p.y));
                }
                coedges.push((new_edge, orientation, old_curve, pcurve, first, last));
            }
            new_wires.push((coedges, wire.orientation));
        }
        // (Sampled pcurves may bulge a little beyond their samples: a margin, within the
        // surface's own bounds where not periodic.)
        let (bu1, bu2, bv1, bv2) = surface.bounds();
        let (du, dv) = (0.01 * (u2 - u1), 0.01 * (v2 - v1));
        let widen = |a: f64, b: f64, d: f64, lo: f64, hi: f64, period: Option<f64>| {
            if matches!(surface, Surface::Plane(_) | Surface::BSpline(_)) {
                (a, b)
            } else if let Some(period) = period {
                // (Not beyond a period, which trimming would wrap round.)
                let d = d.min(0.5 * (period - (b - a))).max(0.0);
                (a - d, b + d)
            } else {
                ((a - d).max(lo), (b + d).min(hi))
            }
        };
        let (u1, u2) = widen(u1, u2, du, bu1, bu2, surface.u_period());
        let (v1, v2) = widen(v1, v2, dv, bv1, bv2, surface.v_period());
        let new_surface = surface_to_nurbs(&surface, (u1, u2, v1, v2))?;
        let surface_kept = matches!(surface, Surface::Plane(_) | Surface::BSpline(_) | Surface::Bezier(_));
        model.set_face_surface(copy, new_surface.clone(), tolerance)?;
        for (coedges, orientation) in new_wires {
            let mut ids = Vec::new();
            for (edge, o, old_curve, pcurve, first, last) in coedges {
                let new_curve = model.edge_curve(Shape::from(edge)).map(|(c, _, _)| c);
                let (pcurve, deviation) = new_pcurve(old_curve.as_ref(), new_curve.as_ref(), &pcurve, first, last, &surface, &new_surface, surface_kept)?;
                if deviation > model.tolerance(Shape::from(edge)) {
                    model.set_edge_tolerance(edge, 1.01 * deviation)?;
                    // (A vertex's tolerance covers its edges'.)
                    let (a, b) = model.edge_vertices(Shape::from(edge));
                    for v in [a, b] {
                        let id = VertexId::try_from(v.node).expect("an edge's vertex");
                        if model.tolerance(v) < 1.01 * deviation {
                            model.set_vertex_tolerance(id, 1.01 * deviation)?;
                        }
                    }
                }
                let c = model.add_coedge(edge, o)?;
                model.set_pcurve_with_range(c, pcurve, first, last)?;
                ids.push(c);
            }
            let wire = model.add_wire(&ids)?;
            model.add(copy, Shape { orientation, ..Shape::from(wire) })?;
        }
        Ok(copy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Point, Vector};
    use crate::primitives::Cuboid;
    use crate::{boolean, properties, validity};

    #[test]
    fn nurbs_box() {
        let mut model = Model::new();
        let block = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let moved = block.moved(&Transform::translation(Vector::new(1.0, 0.0, 0.0)));
        let nurbs = to_nurbs(&mut model, moved).unwrap();
        let report = validity::check(&model, nurbs);
        assert!(report.problems.is_empty(), "{report:?}");
        assert!((properties::volume(&model, nurbs).mass - 6.0).abs() < 1e-12);
        assert!((properties::area(&model, nurbs).mass - 22.0).abs() < 1e-12);
        assert_eq!(model.unique_subshapes(nurbs, NodeKind::Edge).len(), 12);
        for face in model.explore(nurbs, NodeKind::Face) {
            assert!(matches!(model.surface(face), Some(Surface::BSpline(_))));
        }
        let bbox_corner = model.unique_subshapes(nurbs, NodeKind::Vertex).iter().map(|&v| model.point(v)).fold(Point::ORIGIN, |m, p| Point::new(m.x.max(p.x), m.y.max(p.y), m.z.max(p.z)));
        assert_eq!(bbox_corner, Point::new(2.0, 2.0, 3.0));
        // Booleans see through the conversion.
        let other = Cuboid::from_corners(Point::new(0.0, 1.0, 0.0), Point::new(3.0, 3.0, 1.0)).build(&mut model).unwrap();
        let fused = boolean::union(&mut model, nurbs, other).unwrap();
        let report = validity::check(&model, fused);
        assert!(report.problems.is_empty(), "{report:?}");
        assert!((properties::volume(&model, fused).mass - (6.0 + 6.0 - 1.0)).abs() < 1e-9);
    }

    #[test]
    fn nurbs_round_solids() {
        use crate::primitives::{Cone, Cylinder, Sphere, Torus};
        use std::f64::consts::PI;
        let mut model = Model::new();
        let solids = [
            (Cylinder::new(2.0, 3.0).build(&mut model).unwrap(), PI * 4.0 * 3.0),
            (Sphere::new(1.5).build(&mut model).unwrap(), 4.0 / 3.0 * PI * 1.5f64.powi(3)),
            (Cone::new(2.0, 1.0, 3.0).build(&mut model).unwrap(), PI * 3.0 / 3.0 * (4.0 + 2.0 + 1.0)),
            (Torus::new(3.0, 1.0).build(&mut model).unwrap(), 2.0 * PI * PI * 3.0),
        ];
        for (solid, volume) in solids {
            let nurbs = to_nurbs(&mut model, solid).unwrap();
            let report = validity::check(&model, nurbs);
            assert!(report.problems.is_empty(), "{report:?}");
            for face in model.explore(nurbs, NodeKind::Face) {
                assert!(matches!(model.surface(face), Some(Surface::BSpline(_))));
            }
            for edge in model.explore(nurbs, NodeKind::Edge) {
                assert!(model.edge_curve(edge).is_none_or(|(c, _, _)| matches!(c, Curve::BSpline(_))));
                assert!(model.tolerance(edge) < 1e-6, "{}", model.tolerance(edge));
            }
            let v = properties::volume(&model, nurbs).mass;
            assert!((v - volume).abs() < 1e-6 * volume, "{v} for {volume}");
        }
    }
}
