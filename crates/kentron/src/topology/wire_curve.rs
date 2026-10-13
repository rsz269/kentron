//! A wire viewed as one curve.

use crate::math::{Point, Vector};

use super::shape::Shape;
use super::Model;

/// The curve running along a wire's edges, in wire order and orientation, parametrized by
/// the concatenation of the edges' parameter ranges.
#[doc(alias = "BRepAdaptor_CompCurve")]
#[derive(Debug, Clone)]
pub struct WireCurve<'a> {
    model: &'a Model,
    wire: Shape,
}

impl<'a> WireCurve<'a> {
    pub fn new(model: &'a Model, wire: Shape) -> Self {
        Self { model, wire }
    }

    /// The wire's edges in order, each with its range and where it starts on the wire's
    /// parameter.
    fn pieces(&self) -> Vec<(Shape, f64, f64, f64)> {
        let mut out = Vec::new();
        let mut at = 0.0;
        for edge in self.model.ordered_edges(self.wire) {
            if self.model.edge_curve(edge).is_none() {
                continue;
            }
            let (first, last) = self.model.edge_range(edge);
            out.push((edge, first, last, at));
            at += last - first;
        }
        out
    }

    pub fn first_parameter(&self) -> f64 {
        0.0
    }

    pub fn last_parameter(&self) -> f64 {
        self.pieces().last().map_or(0.0, |&(_, first, last, at)| at + (last - first))
    }

    pub fn point(&self, t: f64) -> Point {
        self.d1(t).0
    }

    /// Point and first derivative at `t`.
    pub fn d1(&self, t: f64) -> (Point, Vector) {
        use crate::geometry::ParametricCurve;
        let (edge, u) = self.edge_at(t);
        let (curve, _, _) = self.model.edge_curve(edge).expect("an edge with a curve");
        let (p, d) = curve.d1(u);
        (p, if edge.orientation == super::Orientation::Reversed { -d } else { d })
    }

    /// The edge containing parameter `t`, and the corresponding parameter on that edge.
    pub fn edge_at(&self, t: f64) -> (Shape, f64) {
        let pieces = self.pieces();
        let piece = pieces
            .iter()
            .find(|&&(_, first, last, at)| t <= at + (last - first))
            .or(pieces.last())
            .copied()
            .expect("a wire with edges");
        let (edge, first, last, at) = piece;
        let s = t - at;
        let u = if edge.orientation == super::Orientation::Reversed { last - s } else { first + s };
        (edge, u)
    }
}
