//! OCCT's `BRepAdaptor_CompCurve`, ported: a wire's edges read as one curve.
//!
//! The ported algorithms that take a wire (the guide of a pipe, the sections that
//! `BRepFill_CompatibleWires` makes compatible) cut it into new edges and build new wires
//! from them; those are kept here as owned [`EdgeWire`]s rather than in a [`Model`]: their
//! edges in `TopoDS_Iterator` order, each with its curve, range, orientation and vertices,
//! and identities standing for OCCT's `IsSame`. [`EdgeWire::explore`] orders them as
//! `BRepTools_WireExplorer` (without a face) does.

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;
use crate::math::{Point, Transform, Vector};

use super::entities::Orientation;
use super::id::{NodeId, NodeKind};
use super::shape::Shape;
use super::{Model, WireId};

/// Identities for vertices, edges and wires, standing for OCCT's `IsSame`: one per node and
/// placement met in a model, and a new one for each shape the ported algorithms make.
#[derive(Debug, Default)]
pub(crate) struct ShapeKeys {
    seen: Vec<(NodeId, Transform)>,
    made: usize,
}

impl ShapeKeys {
    pub(crate) fn new() -> ShapeKeys {
        ShapeKeys::default()
    }

    /// The identity of a shape of the model (its node at its placement).
    pub(crate) fn of(&mut self, shape: Shape) -> usize {
        let index = match self.seen.iter().position(|&(n, t)| n == shape.node && t == shape.placement) {
            Some(i) => i,
            None => {
                self.seen.push((shape.node, shape.placement));
                self.seen.len() - 1
            }
        };
        // (Even identities for the model's shapes, odd ones for those made.)
        2 * index
    }

    /// A new identity.
    pub(crate) fn fresh(&mut self) -> usize {
        self.made += 1;
        2 * self.made - 1
    }
}

/// A vertex of an [`EdgeWire`] (`TopoDS_Vertex` with its `BRep_TVertex`), with the
/// orientation of the use it was read from.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WireVertex {
    pub key: usize,
    pub point: Point,
    pub tolerance: f64,
    pub orientation: Orientation,
}

impl WireVertex {
    #[doc(alias = "TopoDS_Shape::IsSame")]
    pub(crate) fn is_same(&self, other: &WireVertex) -> bool {
        self.key == other.key
    }

    /// The same vertex oriented `orientation` (`TopoDS_Shape::Oriented`).
    pub(crate) fn oriented(&self, orientation: Orientation) -> WireVertex {
        WireVertex { orientation, ..*self }
    }
}

/// `IsSame` of two vertices that may be null (two null ones are the same).
pub(crate) fn same_vertex(a: Option<WireVertex>, b: Option<WireVertex>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.is_same(&b),
        (None, None) => true,
        _ => false,
    }
}

/// An edge of an [`EdgeWire`] (`TopoDS_Edge` with its `BRep_TEdge`): its 3D curve (with the
/// placement applied, as `BRep_Tool::Curve` gives it) over its range, its orientation in the
/// wire, and its vertices as stored, the first at the start of the range (oriented
/// forward), the second at its end (reversed).
#[derive(Debug, Clone)]
pub(crate) struct WireEdge {
    pub key: usize,
    /// `BRepAdaptor_Curve(E)`: the curve over the edge's range; `None` for a degenerate edge.
    pub curve: Option<GeomAdaptor>,
    /// `BRep_Tool::Range`.
    pub first: f64,
    pub last: f64,
    /// `BRep_Tool::Tolerance(E)`.
    pub tolerance: f64,
    pub orientation: Orientation,
    pub vertices: [Option<WireVertex>; 2],
    /// `BRep_Tool::Degenerated`.
    pub degenerated: bool,
}

impl WireEdge {
    #[doc(alias = "TopoDS_Shape::IsSame")]
    pub(crate) fn is_same(&self, other: &WireEdge) -> bool {
        self.key == other.key
    }

    /// `Reversed()`.
    pub(crate) fn reversed(&self) -> WireEdge {
        WireEdge { orientation: self.orientation.reversed(), ..self.clone() }
    }

    /// The edge's vertices as `TopoDS_Iterator(E, cumOri)` meets them: as stored, their
    /// orientations composed with the edge's when `cumulate`.
    pub(crate) fn children(&self, cumulate: bool) -> Vec<WireVertex> {
        let orientations = [Orientation::Forward, Orientation::Reversed];
        self.vertices
            .iter()
            .zip(orientations)
            .filter_map(|(v, o)| v.map(|v| v.oriented(if cumulate { o.compose(self.orientation) } else { o })))
            .collect()
    }

    /// `TopExp::Vertices(E, VFirst, VLast, CumOri)`: the vertex oriented forward and the
    /// one oriented reversed among the edge's children.
    #[doc(alias = "TopExp::Vertices")]
    pub(crate) fn ends(&self, cumulate: bool) -> (Option<WireVertex>, Option<WireVertex>) {
        let (mut first, mut last) = (None, None);
        for v in self.children(cumulate) {
            match v.orientation {
                Orientation::Forward => first = Some(v),
                Orientation::Reversed => last = Some(v),
                _ => {}
            }
        }
        (first, last)
    }

    /// `TopExp::FirstVertex(E, CumOri)`.
    pub(crate) fn first_vertex(&self, cumulate: bool) -> Option<WireVertex> {
        self.ends(cumulate).0
    }

    /// `TopExp::LastVertex(E, CumOri)`.
    pub(crate) fn last_vertex(&self, cumulate: bool) -> Option<WireVertex> {
        self.ends(cumulate).1
    }

    /// `BRep_Tool::Parameter(V, E)` for a vertex of the edge (its parameter by the
    /// orientation it has in the edge; on a closed edge, by the orientation `v` has). `None`
    /// where OCCT looks for a point representation (a vertex not of the edge), where it
    /// throws `Standard_NoSuchObject` without one.
    #[doc(alias = "BRep_Tool::Parameter")]
    pub(crate) fn vertex_parameter(&self, v: &WireVertex) -> Option<f64> {
        let mut found: Option<WireVertex> = None;
        let mut rev = false;
        for current in self.children(false) {
            if v.is_same(&current) {
                if found.is_none() {
                    found = Some(current);
                } else {
                    rev = self.orientation == Orientation::Reversed;
                    if current.orientation == v.orientation {
                        found = Some(current);
                    }
                }
            }
        }
        match found.map(|f| f.orientation) {
            Some(Orientation::Forward) => Some(if rev { self.last } else { self.first }),
            Some(Orientation::Reversed) => Some(if rev { self.first } else { self.last }),
            _ => None,
        }
    }

    /// `TopExp::CommonVertex(E1, E2, V)`.
    #[doc(alias = "TopExp::CommonVertex")]
    pub(crate) fn common_vertex(&self, other: &WireEdge) -> Option<WireVertex> {
        let (first1, last1) = self.ends(false);
        let (first2, last2) = other.ends(false);
        let same = |a: Option<WireVertex>, b: Option<WireVertex>| matches!((a, b), (Some(a), Some(b)) if a.is_same(&b)) || (a.is_none() && b.is_none());
        if same(first1, first2) || same(first1, last2) {
            return first1;
        }
        if same(last1, first2) || same(last1, last2) {
            return last1;
        }
        None
    }
}

/// A wire as the ported algorithms read and make one: its edges in `TopoDS_Iterator` order
/// (orientations composed with the wire's) and its `Closed()` flag.
#[derive(Debug, Clone)]
pub(crate) struct EdgeWire {
    pub edges: Vec<WireEdge>,
    pub closed: bool,
}

/// One step of [`EdgeWire::explore`]: the index of the edge in [`EdgeWire::edges`], and
/// `CurrentVertex()`, the vertex it starts from (oriented forward).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Explored {
    pub index: usize,
    pub vertex: Option<WireVertex>,
}

impl EdgeWire {
    /// The wire `wire` of `model` as OCCT reads it: its edges in order with their 3D curves
    /// (with placements applied), ranges, orientations and vertices; `Closed()` as kentron
    /// writes the flag in a BRep file. Fails for an edge with no 3D curve that is not
    /// degenerate (OCCT's `BRepAdaptor_Curve` then reads a curve on a surface, not ported).
    pub(crate) fn from_model(model: &Model, wire: Shape, keys: &mut ShapeKeys) -> Result<EdgeWire> {
        let id = WireId::try_from(wire.node).map_err(|_| Error::InvalidInput("not a wire".into()))?;
        let mut edges = Vec::new();
        // (A wire's children are its coedges, each over its edge.)
        let uses: Vec<Shape> = model.children(wire).flat_map(|c| if c.kind() == NodeKind::Edge { vec![c] } else { model.children(c).collect() }).collect();
        for edge in uses {
            if edge.kind() != NodeKind::Edge {
                continue;
            }
            let degenerated = model.is_degenerate(edge);
            let (curve, first, last) = match model.edge_curve(edge) {
                Some((curve, first, last)) => (Some(GeomAdaptor::new(&curve, first, last)), first, last),
                None if degenerated => (None, 0.0, 0.0),
                None => return Err(Error::InvalidInput("an edge of the wire has no 3D curve".into())),
            };
            let (start, end) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
            let vertex = |v: Shape, orientation, keys: &mut ShapeKeys| WireVertex {
                key: keys.of(v),
                point: model.point(v),
                tolerance: model.tolerance(v),
                orientation,
            };
            let vertices = [Some(vertex(start, Orientation::Forward, keys)), Some(vertex(end, Orientation::Reversed, keys))];
            edges.push(WireEdge { key: keys.of(edge), curve, first, last, tolerance: model.tolerance(edge), orientation: edge.orientation, vertices, degenerated });
        }
        Ok(EdgeWire { edges, closed: model.is_closed_wire(id) })
    }

    /// `BRepTools_WireExplorer(W)` (without a face): the edges in the order they connect,
    /// from the free vertex oriented forward of an open wire, or the first edge's start on a
    /// closed one.
    #[doc(alias = "BRepTools_WireExplorer")]
    pub(crate) fn explore(&self) -> Vec<Explored> {
        // `myMap`: per vertex (by identity, in order of first use), the edges starting at
        // it; `vmap`: the vertices met once (an indexed map: a removal moves the last key
        // into the gap).
        let mut map: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut vmap: Vec<WireVertex> = Vec::new();
        let toggle = |vmap: &mut Vec<WireVertex>, v: WireVertex| match vmap.iter().position(|w| w.is_same(&v)) {
            Some(i) => {
                vmap.swap_remove(i);
            }
            None => vmap.push(v),
        };
        for (i, e) in self.edges.iter().enumerate() {
            if matches!(e.orientation, Orientation::Internal | Orientation::External) {
                continue;
            }
            let (v1, v2) = e.ends(true);
            if let Some(v1) = v1 {
                match map.iter_mut().find(|(k, _)| *k == v1.key) {
                    Some((_, l)) => l.push(i),
                    None => map.push((v1.key, vec![i])),
                }
                toggle(&mut vmap, v1.oriented(Orientation::Forward));
            }
            if let Some(v2) = v2 {
                toggle(&mut vmap, v2.oriented(Orientation::Reversed));
            }
            // (Edges with an infinite end: not met in the wires read here.)
        }
        // `myDoubles`: the edges met twice.
        let mut doubles: Vec<usize> = Vec::new();
        for (i, e) in self.edges.iter().enumerate() {
            if self.edges[..i].iter().any(|f| f.is_same(e)) && !doubles.contains(&e.key) {
                doubles.push(e.key);
            }
        }
        let start = if !vmap.is_empty() {
            vmap.iter().find(|v| v.orientation == Orientation::Forward).copied()
        } else {
            self.edges
                .iter()
                .find(|e| !matches!(e.orientation, Orientation::Internal | Orientation::External))
                .and_then(|e| e.ends(true).0)
        };
        let mut out = Vec::new();
        let Some(start) = start else { return out };
        let Some((_, list)) = map.iter_mut().find(|(k, _)| *k == start.key) else { return out };
        if list.is_empty() {
            return out;
        }
        let mut edge = list.remove(0);
        out.push(Explored { index: edge, vertex: self.edges[edge].first_vertex(true) });
        // (`Next()`.)
        while let Some(vertex) = self.edges[edge].last_vertex(true) {
            let Some((_, list)) = map.iter_mut().find(|(k, _)| *k == vertex.key) else { break };
            if list.is_empty() {
                break;
            }
            if list.len() == 1 {
                let next = list[0];
                if !same_vertex(self.edges[next].first_vertex(true), Some(vertex)) {
                    break;
                }
                edge = next;
                list.clear();
            } else {
                // Degenerate edges first, then edges met twice, then the first other one.
                let current = &self.edges[edge];
                if let Some(j) = list.iter().position(|&k| !self.edges[k].is_same(current) && self.edges[k].degenerated) {
                    edge = list.remove(j);
                } else if let Some(j) = list.iter().position(|&k| doubles.contains(&self.edges[k].key) && !self.edges[k].is_same(current)) {
                    edge = list.remove(j);
                } else if let Some(j) = list.iter().position(|&k| !self.edges[k].is_same(current)) {
                    edge = list.remove(j);
                } else {
                    break;
                }
            }
            out.push(Explored { index: edge, vertex: Some(vertex) });
        }
        out
    }

    /// The edges in [`Self::explore`] order, each with `CurrentVertex()`.
    pub(crate) fn explored_edges(&self) -> Vec<(&WireEdge, Option<WireVertex>)> {
        self.explore().into_iter().map(|x| (&self.edges[x.index], x.vertex)).collect()
    }

    /// `TopExp::Vertices(W, Vfirst, Vlast)`: the free vertices of an open wire (the one
    /// oriented forward first), the last edge's end twice on a closed one, none otherwise.
    #[doc(alias = "TopExp::Vertices")]
    pub(crate) fn ends(&self) -> (Option<WireVertex>, Option<WireVertex>) {
        let mut vmap: Vec<WireVertex> = Vec::new();
        let mut last_v2 = None;
        for e in &self.edges {
            let (v1, v2) = if e.orientation == Orientation::Reversed {
                let (a, b) = e.ends(false);
                (b, a)
            } else {
                e.ends(false)
            };
            for v in [v1.map(|v| v.oriented(Orientation::Forward)), v2.map(|v| v.oriented(Orientation::Reversed))].into_iter().flatten() {
                match vmap.iter().position(|w| w.is_same(&v)) {
                    Some(i) => {
                        vmap.remove(i);
                    }
                    None => vmap.push(v),
                }
            }
            last_v2 = v2;
        }
        if vmap.is_empty() {
            (last_v2.map(|v| v.oriented(Orientation::Forward)), last_v2.map(|v| v.oriented(Orientation::Reversed)))
        } else if vmap.len() == 2 {
            let first = vmap.iter().find(|v| v.orientation == Orientation::Forward).copied();
            let last = vmap.iter().find(|v| v.orientation == Orientation::Reversed).copied();
            (first, last)
        } else {
            (None, None)
        }
    }
}

/// An edge's curve as the composite curve keeps it (`BRepAdaptor_Curve`): the edge, and its
/// curve over the edge's range or a part of it.
#[derive(Debug, Clone)]
struct EdgeCurve {
    edge: WireEdge,
    curve: GeomAdaptor,
}

/// OCCT's `BRepAdaptor_CompCurve`: the wire's non-degenerate edges, in
/// `BRepTools_WireExplorer` order, as one curve: the `i`-th edge (from 0) over `[i, i + 1]`
/// (or over its stretch of the wire's length, by curvilinear abscissa), each mapped
/// linearly onto its range, run backwards where it goes against the wire's direction of
/// parsing.
#[doc(alias = "BRepAdaptor_CompCurve")]
#[derive(Debug, Clone)]
pub(crate) struct CompCurve {
    /// `myWire`.
    wire: EdgeWire,
    /// `IsbyAC`: the edges' lengths (in explorer order, the degenerate ones left out), where
    /// the knots are curvilinear abscissae.
    lengths: Option<Vec<f64>>,
    /// `myCurves`.
    curves: Vec<EdgeCurve>,
    /// `myKnots`.
    knots: Vec<f64>,
    /// `TFirst`, `TLast`, `PTol`.
    t_first: f64,
    t_last: f64,
    p_tol: f64,
    /// `CurIndex` (1-based).
    cur_index: usize,
    /// `Forward`.
    forward: bool,
}

impl CompCurve {
    /// `BRepAdaptor_CompCurve(W)` for the wire `wire` of `model`.
    pub(crate) fn new(model: &Model, wire: Shape) -> Result<CompCurve> {
        let wire = EdgeWire::from_model(model, wire, &mut ShapeKeys::new())?;
        CompCurve::from_wire(wire)
    }

    /// `BRepAdaptor_CompCurve(W, false)`: knots at the edges' indices. Fails for a wire
    /// without non-degenerate edges (where OCCT's adaptor is left without curves).
    pub(crate) fn from_wire(wire: EdgeWire) -> Result<CompCurve> {
        CompCurve::initialize(wire, None)
    }

    /// `BRepAdaptor_CompCurve(W, true)`: knots at the curvilinear abscissae of the edges'
    /// ends, with `length` measuring an edge (`GCPnts_AbscissaPoint::Length`).
    pub(crate) fn by_abscissa(wire: EdgeWire, length: &dyn Fn(&GeomAdaptor) -> f64) -> Result<CompCurve> {
        let lengths = wire
            .explored_edges()
            .into_iter()
            .filter(|(e, _)| !e.degenerated)
            .map(|(e, _)| length(e.curve.as_ref().expect("a non-degenerate edge has a curve")))
            .collect();
        CompCurve::initialize(wire, Some(lengths))
    }

    /// `BRepAdaptor_CompCurve(W, AC, First, Last, Tol)`.
    pub(crate) fn trimmed(wire: EdgeWire, lengths: Option<Vec<f64>>, first: f64, last: f64, tolerance: f64) -> Result<CompCurve> {
        let mut c = CompCurve::initialize(wire, lengths)?;
        c.trim_ends(first, last, tolerance);
        Ok(c)
    }

    /// `Initialize(W, AC)`.
    fn initialize(wire: EdgeWire, lengths: Option<Vec<f64>>) -> Result<CompCurve> {
        let explored: Vec<WireEdge> = wire.explored_edges().into_iter().filter(|(e, _)| !e.degenerated).map(|(e, _)| e.clone()).collect();
        let nb_edge = explored.len();
        if nb_edge == 0 {
            return Err(Error::InvalidInput("a wire without edges".into()));
        }
        let cur_index = nb_edge.div_ceil(2);
        let mut knots = vec![0.0; nb_edge + 1];
        let mut curves = Vec::with_capacity(nb_edge);
        for (ii, edge) in explored.into_iter().enumerate() {
            let curve = edge.curve.clone().expect("a non-degenerate edge has a curve");
            knots[ii + 1] = match &lengths {
                Some(l) => knots[ii] + l[ii],
                None => (ii + 1) as f64,
            };
            curves.push(EdgeCurve { edge, curve });
        }
        // The direction of parsing is kept unless the first edge, as it is oriented, goes
        // away from the second.
        let mut forward = true;
        if nb_edge > 2 || (nb_edge == 2 && !wire.closed) {
            let or = curves[0].edge.orientation;
            let vi = curves[0].edge.common_vertex(&curves[1].edge);
            let vl = curves[0].edge.last_vertex(false);
            if same_vertex(vi, vl) {
                if or == Orientation::Reversed {
                    forward = false;
                }
            } else if or != Orientation::Reversed {
                forward = false;
            }
        }
        let t_last = knots[nb_edge];
        Ok(CompCurve { wire, lengths, curves, knots, t_first: 0.0, t_last, p_tol: 0.0, cur_index, forward })
    }

    /// The trimming part of `Initialize(W, AC, First, Last, Tol)`: the extreme edges cut
    /// at `first` and `last`.
    fn trim_ends(&mut self, first: f64, last: f64, tolerance: f64) {
        self.t_first = first;
        self.t_last = last;
        self.p_tol = tolerance;
        let (mut i1, mut i2) = (self.cur_index, self.cur_index);
        let (f, _) = self.prepare(self.t_first, &mut i1);
        let (l, _) = self.prepare(self.t_last, &mut i2);
        self.cur_index = (i1 + i2) / 2;
        if i1 == i2 {
            let c = &mut self.curves[i1 - 1].curve;
            *c = if l > f { c.trim(f, l) } else { c.trim(l, f) };
        } else {
            let c1 = &mut self.curves[i1 - 1].curve;
            let k = c1.last_parameter();
            *c1 = if k > f { c1.trim(f, k) } else { c1.trim(k, f) };
            let c2 = &mut self.curves[i2 - 1].curve;
            let k = c2.first_parameter();
            *c2 = if k <= l { c2.trim(k, l) } else { c2.trim(l, k) };
        }
    }

    /// `Wire()`.
    pub(crate) fn wire(&self) -> &EdgeWire {
        &self.wire
    }

    /// `Edge(U, E, UonE)`: the edge carrying parameter `u` and the parameter on it.
    pub(crate) fn edge(&self, u: f64) -> (&WireEdge, f64) {
        let mut index = self.cur_index;
        let (u_on_edge, _) = self.prepare(u, &mut index);
        (&self.curves[index - 1].edge, u_on_edge)
    }

    /// `Continuity()`: C0 for several edges, the edge's own for one.
    pub(crate) fn continuity(&self) -> Continuity {
        if self.curves.len() > 1 {
            return Continuity::C0;
        }
        self.curves[0].curve.continuity()
    }

    /// `Trim(First, Last, Tol)`: the same wire between `first` and `last`.
    pub(crate) fn trim(&self, first: f64, last: f64, tolerance: f64) -> CompCurve {
        CompCurve::trimmed(self.wire.clone(), self.lengths.clone(), first, last, tolerance).expect("the wire has edges")
    }

    /// `IsClosed()`: the wire's `Closed()` flag.
    pub(crate) fn is_closed(&self) -> bool {
        self.wire.closed
    }

    /// `GetType()`: OCCT 8.1 gives `GeomAbs_OtherCurve` whatever the wire (its one-edge
    /// case commented out as temporary); `None` here.
    #[doc(alias = "GetType")]
    pub(crate) fn curve_type(&self) -> Option<CurveKind> {
        None
    }

    /// The curve `Line()`, `Circle()`, `BSpline()`, ... read their answer from: the first
    /// edge's (with its placement applied).
    pub(crate) fn basis_curve(&self) -> &Curve {
        self.curves[0].curve.curve()
    }

    /// `Prepare(W, Delta, CurIndex)`: the edge (1-based index, searched from `index`)
    /// carrying the wire's parameter `w`, the parameter on it and the ratio of the edge's
    /// parameter to the wire's (negative where the edge runs backwards). Near a knot, the
    /// edge after it for a parameter in the first half of the wire, the one before in the
    /// second, by `PTol`.
    #[doc(alias = "Prepare")]
    fn prepare(&self, w: f64, index: &mut usize) -> (f64, f64) {
        let eps = if w - self.t_first < self.t_last - w { self.p_tol } else { -self.p_tol };
        let w_test = w + eps;
        let knot = |i: usize| self.knots[i - 1];
        let mut found = false;
        if knot(*index) > w_test {
            let mut ii = *index - 1;
            while ii > 0 && !found {
                if knot(ii) <= w_test {
                    *index = ii;
                    found = true;
                }
                ii -= 1;
            }
            if !found {
                *index = 1;
            }
        } else if knot(*index + 1) <= w_test {
            let mut ii = *index + 1;
            while ii <= self.curves.len() && !found {
                if knot(ii + 1) > w_test {
                    *index = ii;
                    found = true;
                }
                ii += 1;
            }
            if !found {
                *index = self.curves.len();
            }
        }
        let edge = &self.curves[*index - 1].edge;
        let or = edge.orientation;
        let reverse = (self.forward && or == Orientation::Reversed) || (!self.forward && or != Orientation::Reversed);
        let (f, l) = (edge.first, edge.last);
        let mut delta = knot(*index + 1) - knot(*index);
        if delta > self.p_tol * 1.0e-9 {
            delta = (l - f) / delta;
        }
        if reverse {
            delta *= -1.0;
            (l + (w - knot(*index)) * delta, delta)
        } else {
            (f + (w - knot(*index)) * delta, delta)
        }
    }

    /// `InvPrepare(index, First, Delta)`: for the edge `index` (1-based), the parameters
    /// of `T = Ti + (t - First) * Delta`.
    #[doc(alias = "InvPrepare")]
    fn inv_prepare(&self, index: usize) -> (f64, f64) {
        let edge = &self.curves[index - 1].edge;
        let or = edge.orientation;
        let reverse = (self.forward && or == Orientation::Reversed) || (!self.forward && or != Orientation::Reversed);
        let (f, l) = (edge.first, edge.last);
        let mut delta = self.knots[index] - self.knots[index - 1];
        if l - f > self.p_tol * 1.0e-9 {
            delta /= l - f;
        }
        if reverse { (l, -delta) } else { (f, delta) }
    }

    /// The edge (1-based) and the parameter on it, with the ratio, for `u` (`Prepare` from
    /// `CurIndex`).
    fn locate(&self, u: f64) -> (&GeomAdaptor, f64, f64) {
        let mut index = self.cur_index;
        let (u, d) = self.prepare(u, &mut index);
        (&self.curves[index - 1].curve, u, d)
    }
}

impl CurveAdaptor for CompCurve {
    fn first_parameter(&self) -> f64 {
        self.t_first
    }

    fn last_parameter(&self) -> f64 {
        self.t_last
    }

    /// `NbIntervals(S)`, `Intervals(T, S)`: each edge's intervals, mapped onto the wire's
    /// parameter in the direction of parsing.
    fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let mut t = Vec::new();
        for ii in 1..=self.curves.len() {
            let ti = self.curves[ii - 1].curve.intervals(continuity);
            let (f, delta) = self.inv_prepare(ii);
            let big_f = self.knots[ii - 1];
            // (The first edge's first bound, then each edge's others.)
            let skip = if ii == 1 { 0 } else { 1 };
            if delta < 0.0 {
                for &x in ti.iter().rev().skip(skip) {
                    t.push(big_f + (x - f) * delta);
                }
            } else {
                for &x in ti.iter().skip(skip) {
                    t.push(big_f + (x - f) * delta);
                }
            }
        }
        t
    }

    /// `IsPeriodic()`: never.
    fn is_periodic(&self) -> bool {
        false
    }

    /// `Period()`: `TLast - TFirst`.
    fn period(&self) -> f64 {
        self.t_last - self.t_first
    }

    fn value(&self, u: f64) -> Point {
        let (c, u, _) = self.locate(u);
        c.value(u)
    }

    fn d1(&self, u: f64) -> (Point, Vector) {
        let (c, u, d) = self.locate(u);
        let (p, d1) = c.d1(u);
        (p, d1 * d)
    }

    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        let (c, u, d) = self.locate(u);
        let (p, d1, d2) = c.d2(u);
        (p, d1 * d, d2 * (d * d))
    }

    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
        let (c, u, d) = self.locate(u);
        let (p, d1, d2, d3) = c.d3(u);
        (p, d1 * d, d2 * (d * d), d3 * (d * d * d))
    }

    fn dn(&self, u: f64, n: u32) -> Vector {
        let (c, u, d) = self.locate(u);
        c.dn(u, n) * d.powi(n as i32)
    }

    /// `Resolution(R3d)`: the smallest of the edges' (on their own parameters).
    fn resolution(&self, r3d: f64) -> f64 {
        let mut res = 1.0e200;
        for c in &self.curves {
            let r = c.curve.resolution(r3d);
            if r < res {
                res = r;
            }
        }
        res
    }
}

#[cfg(test)]
#[path = "comp_curve_tests.rs"]
pub(crate) mod tests;
