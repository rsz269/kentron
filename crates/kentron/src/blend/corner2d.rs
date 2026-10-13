//! Rounding and bevelling the corner where two edges of a plane meet (OCCT's ChFi2d).

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::{measure, ParametricCurve};
use crate::math::{Point, Vector};
use crate::topology::{Model, NodeKind, Orientation, Shape};
use crate::{make, tolerance};

/// A corner rounded or bevelled: the two edges cut back, and the edge between them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Corner2d {
    /// The first edge, from its far end to where the corner edge begins.
    pub first: Shape,
    /// The fillet's arc, or the chamfer's segment, from the first edge to the second.
    pub corner: Shape,
    /// The second edge, from where the corner edge ends to its far end.
    pub second: Shape,
}

/// An edge seen from the corner: its curve, the parameter at the corner and at the far
/// end.
struct Arm {
    curve: Curve,
    at_corner: f64,
    far: f64,
}

impl Arm {
    /// The edge's curve, oriented from the vertex it shares with `other` (within
    /// tolerance).
    fn of(model: &Model, edge: Shape, other: Shape) -> Result<(Arm, Point)> {
        let (curve, first, last) = model
            .edge_curve(Shape { orientation: Orientation::Forward, ..edge })
            .ok_or_else(|| Error::InvalidInput("a degenerate edge has no corner".into()))?;
        let (oc, of, ol) = model
            .edge_curve(Shape { orientation: Orientation::Forward, ..other })
            .ok_or_else(|| Error::InvalidInput("a degenerate edge has no corner".into()))?;
        let others = [oc.point(of), oc.point(ol)];
        let tol = 1e-6 * (curve.point(first).distance(curve.point(last)) + 1.0);
        for (at, far) in [(first, last), (last, first)] {
            let p = curve.point(at);
            if others.iter().any(|q| q.distance(p) <= tol) {
                return Ok((Arm { curve, at_corner: at, far }, p));
            }
        }
        Err(Error::InvalidInput("the edges do not meet".into()))
    }

    /// The parameter `length` along the edge from the corner.
    fn at_length(&self, length: f64) -> f64 {
        let sign = if self.far > self.at_corner { 1.0 } else { -1.0 };
        measure::parameter_at_length(&self.curve, self.at_corner, sign * length)
    }

    /// The unit direction leaving the corner, at parameter `t`.
    fn direction(&self, t: f64) -> Vector {
        let d = self.curve.d1(t).1;
        let d = if self.far > self.at_corner { d } else { -d };
        d.normalized()
    }

    fn within(&self, t: f64) -> bool {
        let (lo, hi) = (self.at_corner.min(self.far), self.at_corner.max(self.far));
        let slack = 1e-9 * (hi - lo).abs().max(1.0);
        t >= lo - slack && t <= hi + slack
    }

    /// The edge from its far end to `t` (as the first edge of a corner) or from `t` to its
    /// far end (as the second).
    fn trimmed(&self, model: &mut Model, t: f64, leading: bool) -> Result<Shape> {
        let (a, b) = (t.min(self.far), t.max(self.far));
        let edge = make::edge_from_curve_range(model, self.curve.clone(), a, b)?;
        // Running away from the corner if `t < far`; as a first edge it runs towards it.
        let away = self.far > t;
        Ok(if away == leading { edge.reversed() } else { edge })
    }
}

/// The corner where edges `first` and `second` (in one plane) meet, rounded by an arc of
/// `radius` tangent to both: the edges cut back to where the arc touches them.
///
/// The arc's centre lies inside the corner (on the side of the smaller angle between the
/// edges), where the curves offset by `radius` towards each other meet; the edges may be
/// lines, arcs or other planar curves.
#[doc(alias = "ChFi2d_FilletAPI")]
#[doc(alias = "BRepFilletAPI_MakeFillet2d::AddFillet")]
pub fn fillet_2d(model: &mut Model, first: Shape, second: Shape, radius: f64) -> Result<Corner2d> {
    if !(radius > tolerance::LINEAR) {
        return Err(Error::InvalidInput(format!("fillet radius must be positive, got {radius}")));
    }
    let (a, corner) = Arm::of(model, first, second)?;
    let (b, _) = Arm::of(model, second, first)?;
    let (ta, tb) = (a.direction(a.at_corner), b.direction(b.at_corner));
    if ta.dot(tb) <= -1.0 + tolerance::ANGULAR {
        return Err(Error::InvalidInput("the edges meet tangentially: no corner to round".into()));
    }
    // Points a little along each edge, for the plane and the sides where the tangents
    // leave the corner together (a cusp).
    let probe = |arm: &Arm| arm.curve.point(arm.at_length(0.5 * radius.min(measure::length_between(&arm.curve, arm.at_corner.min(arm.far), arm.at_corner.max(arm.far)))));
    let (pa, pb) = (probe(&a), probe(&b));
    let normal = [ta.cross(tb), ta.cross(pb - corner), tb.cross(pa - corner)]
        .into_iter()
        .find(|n| n.length() > tolerance::ANGULAR * (1.0 + (pb - corner).length() + (pa - corner).length()))
        .ok_or_else(|| Error::InvalidInput("the edges meet along a line: no corner to round".into()))?
        .normalized();
    // Each curve offset towards the other (inside the corner).
    let side = |dir: Vector, towards: Point| if (towards - corner).dot(normal.cross(dir)) >= 0.0 { 1.0 } else { -1.0 };
    let (sa, sb) = (side(ta, pb), side(tb, pa));
    let offset_a = |s: f64| a.curve.point(s) + normal.cross(a.direction(s)) * (sa * radius);
    let offset_b = |t: f64| b.curve.point(t) + normal.cross(b.direction(t)) * (sb * radius);
    // From the corner of the edges' tangent lines (a cusp: about a radius along).
    let half = 0.5 * ta.dot(tb).clamp(-1.0, 1.0).acos();
    let reach = if half > 1e-3 { radius / half.tan() } else { radius };
    let (mut s, mut t) = (a.at_length(reach), b.at_length(reach));
    for _ in 0..60 {
        let residual = offset_a(s) - offset_b(t);
        if residual.length() <= 1e-12 * (radius + 1.0) {
            break;
        }
        // Gauss-Newton on the two parameters (derivatives by central differences).
        let (hs, ht) = (1e-7 * (1.0 + s.abs()), 1e-7 * (1.0 + t.abs()));
        let ds = (offset_a(s + hs) - offset_a(s - hs)) * (0.5 / hs);
        let dt = (offset_b(t + ht) - offset_b(t - ht)) * (-0.5 / ht);
        let (m11, m12, m22) = (ds.dot(ds), ds.dot(dt), dt.dot(dt));
        let det = m11 * m22 - m12 * m12;
        if det.abs() <= 1e-300 {
            break;
        }
        let (r1, r2) = (ds.dot(residual), dt.dot(residual));
        s -= (m22 * r1 - m12 * r2) / det;
        t -= (m11 * r2 - m12 * r1) / det;
    }
    let center = offset_a(s);
    if !(center.distance(offset_b(t)) <= 1e-7 * (radius + 1.0)) || !a.within(s) || !b.within(t) {
        return Err(Error::InvalidInput(format!("no fillet of radius {radius} fits the corner at {corner:?}")));
    }
    let (pa, pb) = (a.curve.point(s), b.curve.point(t));
    let arc = super::arc(model, center, pa, pb)?;
    Ok(Corner2d { first: a.trimmed(model, s, true)?, corner: arc, second: b.trimmed(model, t, false)? })
}

/// The corner where edges `first` and `second` meet, bevelled by a segment from the point
/// `d1` along the first edge from the corner to the point `d2` along the second.
#[doc(alias = "ChFi2d_ChamferAPI")]
#[doc(alias = "BRepFilletAPI_MakeFillet2d::AddChamfer")]
pub fn chamfer_2d(model: &mut Model, first: Shape, second: Shape, d1: f64, d2: f64) -> Result<Corner2d> {
    if !(d1 > tolerance::LINEAR && d2 > tolerance::LINEAR) {
        return Err(Error::InvalidInput(format!("chamfer distances must be positive, got {d1}, {d2}")));
    }
    let (a, _) = Arm::of(model, first, second)?;
    let (b, _) = Arm::of(model, second, first)?;
    let (s, t) = (a.at_length(d1), b.at_length(d2));
    if !a.within(s) || !b.within(t) {
        return Err(Error::InvalidInput("a chamfer distance is longer than its edge".into()));
    }
    let segment = make::edge_between(model, a.curve.point(s), b.curve.point(t))?;
    Ok(Corner2d { first: a.trimmed(model, s, true)?, corner: segment, second: b.trimmed(model, t, false)? })
}

/// `wire` (in one plane) with its corners rounded by `radius`: at `vertices` (each where
/// two of its edges meet), or at every vertex where two edges meet at an angle if none are
/// given.
#[doc(alias = "BRepFilletAPI_MakeFillet2d")]
pub fn fillet_wire_2d(model: &mut Model, wire: Shape, radius: f64, vertices: Option<&[Shape]>) -> Result<Shape> {
    let mut edges = model.ordered_edges(wire);
    if edges.is_empty() {
        return Err(Error::InvalidInput("a wire without edges".into()));
    }
    let closed = model.explore(wire, NodeKind::Edge).count() > 1 && {
        let (a, _) = model.edge_vertices(edges[0]);
        let (_, b) = model.edge_vertices(edges[edges.len() - 1]);
        model.point(a).distance(model.point(b)) <= tolerance::LINEAR * 10.0
    };
    let wanted = |model: &Model, p: Point| match vertices {
        None => true,
        Some(vs) => vs.iter().any(|&v| model.point(v).distance(p) <= tolerance::LINEAR * 10.0),
    };
    let n = edges.len();
    let joints = if closed { n } else { n - 1 };
    let mut corners: Vec<Option<Shape>> = vec![None; n];
    for i in 0..joints {
        let j = (i + 1) % n;
        let (_, end) = model.edge_vertices(edges[i]);
        if !wanted(model, model.point(end)) {
            continue;
        }
        match fillet_2d(model, edges[i], edges[j], radius) {
            Ok(c) => {
                edges[i] = c.first;
                edges[j] = c.second;
                corners[i] = Some(c.corner);
            }
            // (Where the edges meet tangentially there is no corner; elsewhere it is an
            // error only if asked for.)
            Err(e) if vertices.is_some() => return Err(e),
            Err(_) => {}
        }
    }
    let mut out = Vec::with_capacity(2 * n);
    for i in 0..n {
        out.push(edges[i]);
        if let Some(c) = corners[i] {
            out.push(c);
        }
    }
    make::wire(model, &out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Circle;
    use crate::math::{Direction, Frame};
    use crate::properties;
    use std::f64::consts::PI;

    fn segment(model: &mut Model, a: Point, b: Point) -> Shape {
        make::edge_between(model, a, b).unwrap()
    }

    #[test]
    fn rounded_right_angle() {
        // OCCT's fillet2d A1: the arc is a quarter circle.
        let mut model = Model::new();
        let e1 = segment(&mut model, Point::ORIGIN, Point::new(100.0, 0.0, 0.0));
        let e2 = segment(&mut model, Point::ORIGIN, Point::new(0.0, 100.0, 0.0));
        let c = fillet_2d(&mut model, e1, e2, 10.0).unwrap();
        assert!((properties::length(&model, c.corner).mass - 5.0 * PI).abs() < 1e-9);
        assert!((properties::length(&model, c.first).mass - 90.0).abs() < 1e-9);
        assert!((properties::length(&model, c.second).mass - 90.0).abs() < 1e-9);
    }

    #[test]
    fn rounded_line_and_arc() {
        // OCCT's fillet2d A3: a line ending where a half circle begins, the wire 253.58 long
        // once rounded.
        let mut model = Model::new();
        let e1 = segment(&mut model, Point::new(200.0, 0.0, 0.0), Point::new(100.0, 0.0, 0.0));
        let circle = Circle::new(Frame::from_z(Point::new(50.0, 0.0, 0.0), Direction::Z), 50.0);
        let e2 = make::edge_from_curve_range(&mut model, Curve::Circle(circle), 0.0, PI).unwrap();
        let wire = make::wire(&mut model, &[e1, e2]).unwrap();
        let rounded = fillet_wire_2d(&mut model, wire, 10.0, None).unwrap();
        let length = properties::length(&model, rounded).mass;
        assert!((length - 253.58).abs() < 0.01, "{length}");
    }

    #[test]
    fn rounded_cusp() {
        // OCCT's fillet2d A5: two half circles leaving their common point the same way (an
        // S-bend's cusp); OCCT's wire is 275.293 long once rounded.
        let mut model = Model::new();
        let half = |model: &mut Model, x: f64| {
            let circle = Circle::new(Frame::from_z(Point::new(x, 0.0, 0.0), Direction::Z), 50.0);
            make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, PI).unwrap()
        };
        let (e1, e2) = (half(&mut model, 150.0), half(&mut model, 50.0));
        let c = fillet_2d(&mut model, e1, e2, 10.0).unwrap();
        let wire = make::wire(&mut model, &[c.first, c.corner, c.second]).unwrap();
        let length = properties::length(&model, wire).mass;
        assert!((length - 275.293).abs() < 1e-3, "{length}");
    }

    #[test]
    fn bevelled_right_angle() {
        // OCCT's chamfer2d A1.
        let mut model = Model::new();
        let e1 = segment(&mut model, Point::ORIGIN, Point::new(100.0, 0.0, 0.0));
        let e2 = segment(&mut model, Point::ORIGIN, Point::new(0.0, 100.0, 0.0));
        let c = chamfer_2d(&mut model, e1, e2, 10.0, 20.0).unwrap();
        assert!((properties::length(&model, c.corner).mass - 500f64.sqrt()).abs() < 1e-9);
    }
}
