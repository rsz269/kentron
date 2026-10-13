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
/// the geometry (the copy is unplaced). Lines and planes convert exactly, keeping their
/// parametrization (a plane becomes a bilinear patch over its face's parameter range);
/// other geometry is not supported yet.
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
    Ok(Curve::BSpline(match curve {
        Curve::Line(_) => segment(curve.point(first), curve.point(last), first, last)?,
        Curve::BSpline(b) => b.clone(),
        other => todo!("NURBS conversion of a {}", kind(other)),
    }))
}

fn pcurve_to_nurbs(curve: &Curve2, first: f64, last: f64) -> Result<Curve2> {
    Ok(Curve2::BSpline(match curve {
        Curve2::Line(_) => {
            let (a, b) = (curve.point(first), curve.point(last));
            segment(Point2::new(a.x, a.y), Point2::new(b.x, b.y), first, last)?
        }
        Curve2::BSpline(b) => b.clone(),
        other => todo!("NURBS conversion of a 2D {}", kind(other)),
    }))
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
        other => todo!("NURBS conversion of a {} surface", kind(other)),
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
                let new_edge = self.edge(model, placed)?;
                let (pcurve, first, last) = model
                    .effective_pcurve(id)
                    .ok_or_else(|| Error::InvalidInput("converting a face without pcurves".into()))?;
                for i in 0..=16 {
                    let p = pcurve.point(first + (last - first) * i as f64 / 16.0);
                    (u1, u2, v1, v2) = (u1.min(p.x), u2.max(p.x), v1.min(p.y), v2.max(p.y));
                }
                coedges.push((new_edge, orientation, pcurve_to_nurbs(&pcurve, first, last)?, first, last));
            }
            new_wires.push((coedges, wire.orientation));
        }
        model.set_face_surface(copy, surface_to_nurbs(&surface, (u1, u2, v1, v2))?, tolerance)?;
        for (coedges, orientation) in new_wires {
            let mut ids = Vec::new();
            for (edge, o, pcurve, first, last) in coedges {
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
}
