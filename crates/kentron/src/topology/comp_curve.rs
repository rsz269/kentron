//! A wire read as one curve, for the guided pipe frames.
//!
//! [`EdgeWire`] is a light, owned copy of a wire: its edges in connection order, each with
//! its curve, range, the direction the wire runs it and its end vertices. Guide preparation
//! turns, opens and cuts such copies without touching the [`Model`].
//!
//! [`CompCurve`] lays the wire's non-degenerate edges end to end as one parametric curve:
//! edge `i` (from 1) takes the parameters `[k_{i−1}, k_i]`, with knots `k_i = i`, mapped
//! affinely onto its own range in the direction the wire runs it. Derivatives scale by the
//! map's slope. An interior knot belongs to the edge after it, unless a trim tolerance moves
//! the decision (see [`CompCurve::trim`]).

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;
use crate::math::{Point, Transform, Vector};

use super::{Model, NodeId, NodeKind, Orientation, Shape, WireId};

/// Identities for wire vertices: model vertices (node and placement) and vertices made by
/// cutting edges.
#[derive(Clone, Debug, Default)]
pub(crate) struct ShapeKeys {
    known: Vec<(NodeId, Transform)>,
    fresh: usize,
}

impl ShapeKeys {
    pub(crate) fn new() -> ShapeKeys {
        ShapeKeys::default()
    }

    /// The identity of a model vertex.
    pub(crate) fn of(&mut self, shape: Shape) -> usize {
        if let Some(i) = self.known.iter().position(|&(n, t)| n == shape.node && t == shape.placement) {
            return 2 * i;
        }
        self.known.push((shape.node, shape.placement));
        2 * (self.known.len() - 1)
    }

    /// A new identity.
    pub(crate) fn fresh(&mut self) -> usize {
        self.fresh += 1;
        2 * self.fresh - 1
    }
}

/// An end of a wire edge.
#[derive(Clone, Copy, Debug)]
pub(crate) struct WireVertex {
    pub key: usize,
    pub point: Point,
    pub tolerance: f64,
}

/// Whether two optional vertices are the same one.
pub(crate) fn same_vertex(a: Option<WireVertex>, b: Option<WireVertex>) -> bool {
    matches!((a, b), (Some(a), Some(b)) if a.key == b.key)
}

/// An edge of a wire as the wire runs it.
#[derive(Clone, Debug)]
pub(crate) struct WireEdge {
    /// The edge's curve over its range; `None` for a degenerate edge.
    pub curve: Option<GeomAdaptor>,
    pub first: f64,
    pub last: f64,
    /// The wire runs the edge from `last` to `first`.
    pub reversed: bool,
    /// The vertices the wire enters and leaves the edge at.
    pub start: Option<WireVertex>,
    pub end: Option<WireVertex>,
    pub degenerate: bool,
}

impl WireEdge {
    /// The edge run the other way.
    pub(crate) fn reversed(&self) -> WireEdge {
        WireEdge { reversed: !self.reversed, start: self.end, end: self.start, ..self.clone() }
    }

    /// The parameter the wire enters the edge at, and the one it leaves it at.
    pub(crate) fn run(&self) -> (f64, f64) {
        if self.reversed { (self.last, self.first) } else { (self.first, self.last) }
    }

    /// The point at the curve parameter `t`.
    pub(crate) fn point(&self, t: f64) -> Option<Point> {
        self.curve.as_ref().map(|c| c.value(t))
    }

    /// The derivative along the run at the curve parameter `t`.
    pub(crate) fn run_derivative(&self, t: f64) -> Option<Vector> {
        let (_, d) = self.curve.as_ref()?.d1(t);
        Some(if self.reversed { -d } else { d })
    }

    /// The edge cut at the curve parameter `t` (strictly inside its range): the part the
    /// wire runs first and the part after, joined at a new vertex `at`.
    pub(crate) fn split(&self, t: f64, at: WireVertex) -> (WireEdge, WireEdge) {
        let low = WireEdge { last: t, ..self.clone() };
        let high = WireEdge { first: t, ..self.clone() };
        let (mut before, mut after) = if self.reversed { (high, low) } else { (low, high) };
        before.end = Some(at);
        after.start = Some(at);
        (before, after)
    }
}

/// A wire's edges in connection order, each as the wire runs it.
#[derive(Clone, Debug)]
pub(crate) struct EdgeWire {
    pub edges: Vec<WireEdge>,
    pub closed: bool,
}

impl EdgeWire {
    /// The wire (or lone edge) `wire` of `model`, its curves placed.
    pub(crate) fn from_model(model: &Model, wire: Shape, keys: &mut ShapeKeys) -> Result<EdgeWire> {
        let (edges, closed) = match wire.kind() {
            NodeKind::Edge => {
                let (a, b) = model.edge_vertices(wire);
                (vec![wire], a.node == b.node)
            }
            NodeKind::Wire => {
                let id = WireId::try_from(wire.node).map_err(|_| Error::InvalidInput("not a wire".into()))?;
                (model.ordered_edges(wire), model.is_closed_wire(id))
            }
            kind => return Err(Error::InvalidInput(format!("a guide or spine that is a {kind:?}"))),
        };
        let mut out = Vec::with_capacity(edges.len());
        for edge in edges {
            let degenerate = model.is_degenerate(edge);
            let (curve, first, last) = match model.edge_curve(edge) {
                Some((c, f, l)) => (Some(GeomAdaptor::new(&c, f, l)), f, l),
                None if degenerate => (None, 0.0, 0.0),
                None => return Err(Error::InvalidInput("an edge without a curve".into())),
            };
            let (a, b) = model.edge_vertices(edge);
            let mut vertex = |v: Shape| WireVertex { key: keys.of(Shape { orientation: Orientation::Forward, ..v }), point: model.point(v), tolerance: model.tolerance(v) };
            let (start, end) = (vertex(a), vertex(b));
            let reversed = edge.orientation == Orientation::Reversed;
            out.push(WireEdge { curve, first, last, reversed, start: Some(start), end: Some(end), degenerate });
        }
        if out.is_empty() {
            return Err(Error::InvalidInput("a wire without edges".into()));
        }
        Ok(EdgeWire { edges: out, closed })
    }

    /// The vertices the wire starts and ends at (for a closed wire the last edge's end,
    /// twice).
    pub(crate) fn ends(&self) -> (Option<WireVertex>, Option<WireVertex>) {
        let last = self.edges.last().and_then(|e| e.end);
        if self.closed {
            return (last, last);
        }
        (self.edges.first().and_then(|e| e.start), last)
    }

    /// The wire run the other way.
    pub(crate) fn reversed(&self) -> EdgeWire {
        EdgeWire { edges: self.edges.iter().rev().map(WireEdge::reversed).collect(), closed: self.closed }
    }
}

/// One edge of a composite curve: its curve (`view`, possibly restricted by a trim) and the
/// range `[first, last]` the affine map runs over.
#[derive(Clone, Debug)]
struct Piece {
    view: GeomAdaptor,
    first: f64,
    last: f64,
    reversed: bool,
}

/// A wire as one curve.
#[derive(Clone, Debug)]
pub(crate) struct CompCurve {
    pieces: Vec<Piece>,
    /// `k_0 … k_n`.
    knots: Vec<f64>,
    first: f64,
    last: f64,
    /// Moves the choice of edge at a knot towards the inside of the range.
    tolerance: f64,
    closed: bool,
}

impl CompCurve {
    /// The wire `wire` of `model` as one curve.
    pub(crate) fn new(model: &Model, wire: Shape) -> Result<CompCurve> {
        CompCurve::from_wire(EdgeWire::from_model(model, wire, &mut ShapeKeys::new())?)
    }

    /// The prepared wire `wire` as one curve, its knots at the edge indices.
    pub(crate) fn from_wire(wire: EdgeWire) -> Result<CompCurve> {
        let pieces: Vec<Piece> = wire
            .edges
            .iter()
            .filter(|e| !e.degenerate)
            .filter_map(|e| Some(Piece { view: e.curve.clone()?, first: e.first, last: e.last, reversed: e.reversed }))
            .collect();
        if pieces.is_empty() {
            return Err(Error::InvalidInput("a wire without a curve".into()));
        }
        let knots: Vec<f64> = (0..=pieces.len()).map(|i| i as f64).collect();
        let last = pieces.len() as f64;
        Ok(CompCurve { pieces, knots, first: 0.0, last, tolerance: 0.0, closed: wire.closed })
    }

    pub(crate) fn is_closed(&self) -> bool {
        self.closed
    }

    /// The slope of piece `i`'s map (negative for a reversed piece).
    fn slope(&self, i: usize) -> f64 {
        let p = &self.pieces[i];
        let d = (p.last - p.first) / (self.knots[i + 1] - self.knots[i]);
        if p.reversed { -d } else { d }
    }

    /// The edge parameter of `w` on piece `i`.
    fn map(&self, i: usize, w: f64) -> f64 {
        let p = &self.pieces[i];
        let along = (w - self.knots[i]) * ((p.last - p.first) / (self.knots[i + 1] - self.knots[i]));
        if p.reversed { p.last - along } else { p.first + along }
    }

    /// The composite parameter of the edge parameter `x` on piece `i`.
    fn unmap(&self, i: usize, x: f64) -> f64 {
        let p = &self.pieces[i];
        let along = if p.reversed { p.last - x } else { x - p.first };
        self.knots[i] + along * ((self.knots[i + 1] - self.knots[i]) / (p.last - p.first))
    }

    /// The piece holding `w`.
    fn holding(&self, w: f64) -> usize {
        let shifted = if w - self.first < self.last - w { w + self.tolerance } else { w - self.tolerance };
        let n = self.pieces.len();
        self.knots[1..n].partition_point(|&k| k <= shifted)
    }

    /// The edge holding `w` (as a piece's curve and range) and the edge parameter there.
    pub(crate) fn edge(&self, w: f64) -> (&GeomAdaptor, f64) {
        let i = self.holding(w);
        (&self.pieces[i].view, self.map(i, w))
    }

    pub(crate) fn continuity(&self) -> Continuity {
        if self.pieces.len() > 1 { Continuity::C0 } else { self.pieces[0].view.continuity() }
    }

    /// The same composite over `[first, last]`, the edge holding `first` restricted to start
    /// there and the one holding `last` to end there; `tolerance` moves the choice of edge
    /// at a knot inwards. Always from the untrimmed edges.
    pub(crate) fn trim(&self, first: f64, last: f64, tolerance: f64) -> CompCurve {
        let mut out = CompCurve { first, last, tolerance, ..self.clone() };
        // (Each piece's view back to its edge's whole range.)
        for p in &mut out.pieces {
            p.view = p.view.trim(p.first, p.last);
        }
        let (i, j) = (out.holding(first), out.holding(last));
        let (x0, x1) = (out.map(i, first), out.map(j, last));
        if i == j {
            let p = &mut out.pieces[i];
            p.view = p.view.trim(x0.min(x1), x0.max(x1));
        } else {
            for (k, x, starts) in [(i, x0, true), (j, x1, false)] {
                let p = &mut out.pieces[k];
                // (The run starts at the low end of a forward edge, the high end of a reversed one.)
                let (lo, hi) = (p.view.first_parameter(), p.view.last_parameter());
                p.view = if starts != p.reversed { p.view.trim(x, hi) } else { p.view.trim(lo, x) };
            }
        }
        out
    }

    /// The derivatives up to `order` at `w`, scaled by the map's slope.
    fn derivative(&self, w: f64, n: u32) -> Vector {
        let i = self.holding(w);
        self.pieces[i].view.dn(self.map(i, w), n) * self.slope(i).powi(n as i32)
    }
}

impl CurveAdaptor for CompCurve {
    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let n = self.pieces.len();
        let mut out = vec![self.first];
        for i in 0..n {
            let lo = if i == 0 { self.first } else { self.knots[i].max(self.first) };
            let hi = if i + 1 == n { self.last } else { self.knots[i + 1].min(self.last) };
            if !(hi > lo) {
                continue;
            }
            let (x0, x1) = (self.map(i, lo), self.map(i, hi));
            let view = self.pieces[i].view.trim(x0.min(x1), x0.max(x1));
            let mut breaks: Vec<f64> = view.intervals(continuity).into_iter().map(|x| self.unmap(i, x)).collect();
            if self.pieces[i].reversed {
                breaks.reverse();
            }
            for w in breaks.into_iter().skip(1) {
                let w = if w > hi { hi } else { w };
                if w > *out.last().expect("a bound") {
                    out.push(w);
                }
            }
            // (The piece's end exactly: its knot, or the range's end.)
            *out.last_mut().expect("a bound") = hi;
        }
        if *out.last().expect("a bound") != self.last {
            out.push(self.last);
        }
        out
    }

    fn is_periodic(&self) -> bool {
        false
    }

    fn period(&self) -> f64 {
        self.last - self.first
    }

    fn value(&self, w: f64) -> Point {
        let i = self.holding(w);
        self.pieces[i].view.value(self.map(i, w))
    }

    fn d1(&self, w: f64) -> (Point, Vector) {
        let i = self.holding(w);
        let (p, d1) = self.pieces[i].view.d1(self.map(i, w));
        (p, d1 * self.slope(i))
    }

    fn d2(&self, w: f64) -> (Point, Vector, Vector) {
        let i = self.holding(w);
        let s = self.slope(i);
        let (p, d1, d2) = self.pieces[i].view.d2(self.map(i, w));
        (p, d1 * s, d2 * (s * s))
    }

    fn d3(&self, w: f64) -> (Point, Vector, Vector, Vector) {
        let i = self.holding(w);
        let s = self.slope(i);
        let (p, d1, d2, d3) = self.pieces[i].view.d3(self.map(i, w));
        (p, d1 * s, d2 * (s * s), d3 * (s * s * s))
    }

    fn dn(&self, w: f64, n: u32) -> Vector {
        self.derivative(w, n)
    }

    fn resolution(&self, r3d: f64) -> f64 {
        self.pieces.iter().map(|p| p.view.resolution(r3d)).fold(f64::INFINITY, f64::min)
    }

    fn basis_curve(&self) -> Option<&Curve> {
        Some(self.pieces[0].view.curve())
    }

    fn trimmed(&self, first: f64, last: f64, tolerance: f64) -> Box<dyn CurveAdaptor> {
        Box::new(self.trim(first, last, tolerance))
    }
}

#[cfg(test)]
#[path = "comp_curve_tests.rs"]
pub(crate) mod tests;
