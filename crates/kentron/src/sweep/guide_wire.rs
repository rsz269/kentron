//! The guide wire of a pipe swept with a guide, as OCCT's `BRepFill_PipeShell::Set(Guide,
//! ...)` prepares it before the guide trihedra receive it as a `BRepAdaptor_CompCurve`
//! ([`CompCurve`](crate::topology::comp_curve::CompCurve)), ported:
//!
//! * an open guide along an open spine goes through `BRepFill_CompatibleWires` (its path for
//!   open sections: `SearchOrigin`, then `SameNumberByACR`), which turns it to run the
//!   spine's way and cuts its edges so that both wires have as many;
//! * a closed guide is cut and turned at the point nearest the spine's start
//!   (`BRepFill::SearchOrigin`);
//! * with `BRepFill::ComputeACR` and `BRepFill::InsertACR`, and the pieces of `BRepLib`
//!   (`BRepLib_MakeEdge`, `BRepLib_MakeWire`), `BRepLProp` and `BRepExtrema` these call.
//!
//! The wires are [`EdgeWire`]s: the new edges are not put in the model.
//!
//! Left out: `BRepFill_CompatibleWires`' paths for closed sections (`ComputeOrigin`,
//! `SameNumberByPolarMethod`) and its history (`Generated`), which the guide does not use.

use std::f64::consts::{FRAC_PI_2, PI};

use crate::error::{Error, Result};
use crate::geometry::adaptor::{epsilon, in_period, CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::ext_pc::{ext_pc_other, ExtPCPoint};
use crate::geometry::{Continuity, ParametricCurve};
use crate::math::{Point, Vector};
use crate::topology::comp_curve::{same_vertex, EdgeWire, ShapeKeys, WireEdge, WireVertex};
use crate::topology::{Model, Orientation, Shape};

use super::frenet::{angle, direction};

/// `Precision::Confusion()`, also `BRepLib::Precision()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1.0e-9;
/// `Precision::Angular()`.
const ANGULAR: f64 = 1.0e-12;
/// `gp::Resolution()`.
const RESOLUTION: f64 = f64::MIN_POSITIVE;
/// `Precision::Infinite()`.
const INFINITE: f64 = 2.0e100;

/// `Precision::IsInfinite`.
fn is_infinite(x: f64) -> bool {
    x.abs() >= 0.5 * INFINITE
}

/// `gp_Vec::Angle`, which throws `gp_VectorWithNullMagnitude` for a null vector.
fn vector_angle(a: Vector, b: Vector) -> Result<f64> {
    if a.length() <= RESOLUTION || b.length() <= RESOLUTION {
        return Err(Error::Failed("gp_Vec::Angle of a null vector".into()));
    }
    Ok(angle(a, b))
}

/// `gp_Vec::IsEqual(Other, LinearTolerance, AngularTolerance)`.
fn vectors_equal(a: Vector, b: Vector, linear: f64, angular: f64) -> bool {
    let (ma, mb) = (a.length(), b.length());
    let val = (ma - mb).abs();
    if ma <= linear || mb <= linear { val <= linear } else { val <= linear && angle(a, b) <= angular }
}

fn curve_of(edge: &WireEdge) -> Result<&GeomAdaptor> {
    edge.curve.as_ref().ok_or_else(|| Error::Failed("an edge without a 3D curve".into()))
}

fn point_of(v: Option<WireVertex>) -> Result<Point> {
    v.map(|v| v.point).ok_or_else(|| Error::Failed("BRep_Tool::Pnt of a null vertex".into()))
}

/// The reduced curvilinear abscissae of a wire's edges (`BRepFill::ComputeACR`): at `i`
/// (from 1) the length of the wire up to the end of its `i`-th edge (in
/// `BRepTools_WireExplorer` order) over the whole length, which is at `0`; for a wire of
/// length below `Precision::Confusion()`, the unreduced lengths with `1` at the end.
///
/// `length` measures an edge's curve over its range: OCCT's `GCPnts_AbscissaPoint::Length`.
#[doc(alias = "BRepFill::ComputeACR")]
pub(crate) fn compute_acr(wire: &EdgeWire, length: &dyn Fn(&GeomAdaptor) -> f64) -> Vec<f64> {
    let edges = wire.explored_edges();
    let nb_edges = edges.len();
    let mut acr = vec![0.0; nb_edges + 1];
    for (i, (edge, _)) in edges.iter().enumerate() {
        acr[i + 1] = acr[i];
        if !edge.degenerated {
            if let Some(curve) = &edge.curve {
                acr[i + 1] += length(curve);
            }
        }
    }
    acr[0] = acr[nb_edges];
    if acr[0] > CONFUSION {
        for i in 1..=nb_edges {
            acr[i] /= acr[0];
        }
    } else {
        // A wire reduced to a point.
        acr[nb_edges] = 1.0;
    }
    acr
}

/// `ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)`: `u1` into `[first, last)`, `u2`
/// into `(u1, u1 + period]`.
fn adjust_periodic(first: f64, last: f64, preci: f64, u1: f64, u2: f64) -> (f64, f64) {
    if is_infinite(first) || is_infinite(last) {
        return (first, last);
    }
    let period = last - first;
    if period < epsilon(last) {
        return (first, last);
    }
    let mut u1 = u1 - ((u1 - first) / period).floor() * period;
    if last - u1 < preci {
        u1 -= period;
    }
    let mut u2 = u2 - ((u2 - u1) / period).floor() * period;
    if u2 - u1 < preci {
        u2 += period;
    }
    (u1, u2)
}

/// A new vertex (`BRep_Builder::MakeVertex(V, P, Tol)`).
fn new_vertex(point: Point, tolerance: f64, keys: &mut ShapeKeys) -> WireVertex {
    WireVertex { key: keys.fresh(), point, tolerance, orientation: Orientation::Forward }
}

/// `BRepLib_MakeEdge(C, V1, V2, p1, p2)`: an edge on `curve` (a trimmed curve's basis)
/// between the parameters, through the vertices given or new ones; the errors where it is
/// not done.
#[doc(alias = "BRepLib_MakeEdge")]
fn make_edge(curve: &Curve, vv1: Option<WireVertex>, vv2: Option<WireVertex>, pp1: f64, pp2: f64, keys: &mut ShapeKeys) -> Result<WireEdge> {
    let mut curve = curve;
    while let Curve::Trimmed(t) = curve {
        curve = t.basis();
    }
    let (mut p1, mut p2) = (pp1, pp2);
    let (cf, cl) = (curve.first_parameter(), curve.last_parameter());
    let (mut v1, mut v2);
    if curve.is_periodic() {
        (p1, p2) = adjust_periodic(cf, cl, P_CONFUSION, p1, p2);
        v1 = vv1;
        v2 = vv2;
    } else {
        if p1 < p2 {
            v1 = vv1;
            v2 = vv2;
        } else {
            v2 = vv1;
            v1 = vv2;
            std::mem::swap(&mut p1, &mut p2);
        }
        if cf - p1 > P_CONFUSION || p2 - cl > P_CONFUSION {
            return Err(Error::Failed("BRepLib_MakeEdge: parameter out of range".into()));
        }
        if p2 - p1 <= RESOLUTION {
            return Err(Error::Failed("BRepLib_MakeEdge: line through identical points".into()));
        }
    }
    let ca = GeomAdaptor::new(curve, cf, cl);
    let p1inf = p1 <= -0.5 * INFINITE;
    let p2inf = p2 >= 0.5 * INFINITE;
    let point1 = if p1inf { Point::ORIGIN } else { ca.value(p1) };
    let point2 = if p2inf { Point::ORIGIN } else { ca.value(p2) };
    let preci = CONFUSION;
    let closed = !p1inf && !p2inf && point1.distance(point2) <= preci;
    let mut degenerated = false;
    if closed {
        match (v1, v2) {
            (None, None) => {
                v1 = Some(new_vertex(point1, preci, keys));
                v2 = v1;
            }
            (None, Some(_)) => v1 = v2,
            (Some(_), None) => v2 = v1,
            (Some(a), Some(b)) => {
                if !a.is_same(&b) || point1.distance(a.point) > preci.max(a.tolerance) {
                    return Err(Error::Failed("BRepLib_MakeEdge: different points on a closed curve".into()));
                }
                let middle = ca.value((p1 + p2) / 2.0);
                if point1.distance(middle) < preci {
                    degenerated = true;
                }
            }
        }
    } else {
        for (v, inf, p) in [(&mut v1, p1inf, point1), (&mut v2, p2inf, point2)] {
            if inf {
                if v.is_some() {
                    return Err(Error::Failed("BRepLib_MakeEdge: a point at an infinite parameter".into()));
                }
            } else {
                match v {
                    None => *v = Some(new_vertex(p, preci, keys)),
                    Some(w) if p.distance(w.point) > preci.max(w.tolerance) => {
                        return Err(Error::Failed("BRepLib_MakeEdge: the point and the parameter differ".into()));
                    }
                    Some(_) => {}
                }
            }
        }
    }
    Ok(WireEdge {
        key: keys.fresh(),
        curve: Some(GeomAdaptor::new(curve, p1, p2)),
        first: p1,
        last: p2,
        tolerance: preci,
        orientation: Orientation::Forward,
        vertices: [v1.map(|v| v.oriented(Orientation::Forward)), v2.map(|v| v.oriented(Orientation::Reversed))],
        degenerated,
    })
}

/// The curve of an edge as `BRep_Tool::Curve(E, first, last)` gives it.
fn edge_curve(edge: &WireEdge) -> Result<&Curve> {
    Ok(curve_of(edge)?.curve())
}

/// `CutEdge` (BRepFill.cxx): `edge` cut at `param`, the piece from `v_ref` first.
fn cut_edge(edge: &WireEdge, param: f64, v_ref: Option<WireVertex>, keys: &mut ShapeKeys) -> Result<(WireEdge, WireEdge)> {
    let c = edge_curve(edge)?;
    let (first, last) = (edge.first, edge.last);
    let vi = Some(new_vertex(curve_of(edge)?.value(param), CONFUSION, keys));
    let (vf, vl) = edge.ends(false);
    if same_vertex(v_ref, vf) {
        let e1 = make_edge(c, vf, vi, first, param, keys)?;
        let e2 = make_edge(c, vi, vl, param, last, keys)?;
        Ok((e1, e2))
    } else {
        let e2 = make_edge(c, vf, vi, first, param, keys)?;
        let e1 = make_edge(c, vi, vl, param, last, keys)?;
        Ok((e1, e2))
    }
}

/// `TrimEdge` (BRepFill.cxx): `edge` cut at the reduced abscissae `cuts` (of `[t0, t1]`,
/// mapped linearly onto its range), the pieces from its first vertex on if `seq_order`,
/// from its last otherwise; each oriented as `edge`.
fn trim_edge(edge: &WireEdge, cuts: &[f64], t0: f64, t1: f64, seq_order: bool, keys: &mut ShapeKeys) -> Result<Vec<WireEdge>> {
    let ndec = cuts.len();
    let c = edge_curve(edge)?;
    let (first, last) = (edge.first, edge.last);
    let orientation = edge.orientation;
    let (vf, vl) = edge.ends(false);
    let mut s = Vec::with_capacity(ndec + 1);
    let oriented = |mut e: WireEdge| {
        e.orientation = orientation;
        e
    };
    if seq_order {
        let mut m0 = first;
        let mut v0 = vf;
        for (j, &cut) in cuts.iter().enumerate() {
            let m1 = (cut - t0) * (last - first) / (t1 - t0) + first;
            let piece = oriented(make_edge(c, v0, None, m0, m1, keys)?);
            m0 = m1;
            v0 = piece.last_vertex(false);
            s.push(piece);
            if j == ndec - 1 {
                s.push(oriented(make_edge(c, v0, vl, m0, last, keys)?));
            }
        }
    } else {
        let mut m1 = last;
        let mut v1 = vl;
        for j in (0..ndec).rev() {
            let m0 = (cuts[j] - t0) * (last - first) / (t1 - t0) + first;
            let piece = oriented(make_edge(c, None, v1, m0, m1, keys)?);
            m1 = m0;
            v1 = piece.first_vertex(false);
            s.push(piece);
            if j == 0 {
                s.push(oriented(make_edge(c, vf, v1, first, m1, keys)?));
            }
        }
    }
    Ok(s)
}

/// `BRepLib_MakeWire`, edge by edge (`Add(E)`): each edge connected to the wire by a vertex
/// it shares (or one within tolerance, which it then shares, moved between both), and
/// turned to run on from the wire's end or up to its start.
#[doc(alias = "BRepLib_MakeWire")]
#[derive(Default)]
struct MakeWire {
    edges: Vec<WireEdge>,
    /// `myVertices`.
    vertices: Vec<WireVertex>,
    /// `VF`, `VL`: the wire's ends.
    vf: Option<WireVertex>,
    vl: Option<WireVertex>,
    closed: bool,
}

impl MakeWire {
    fn add_vertex(&mut self, v: WireVertex) {
        if !self.vertices.iter().any(|w| w.is_same(&v)) {
            self.vertices.push(v);
        }
    }

    /// `B.UpdateVertex(V, P, Tol)`: the vertex moved, wherever it is used.
    fn update_vertex(&mut self, key: usize, point: Point, tolerance: f64) {
        let update = |v: &mut WireVertex| {
            if v.key == key {
                v.point = point;
                v.tolerance = tolerance;
            }
        };
        self.vertices.iter_mut().for_each(update);
        for e in &mut self.edges {
            e.vertices.iter_mut().flatten().for_each(update);
        }
        for v in [&mut self.vf, &mut self.vl].into_iter().flatten() {
            update(v);
        }
    }

    /// `Add(E)` (checking the proximity of vertices).
    fn add(&mut self, e: &WireEdge, keys: &mut ShapeKeys) -> Result<()> {
        let init;
        let mut my_edge;
        // (`myError` set to `BRepLib_NonManifoldWire`.)
        let mut non_manifold = false;
        let mut my_vertex: Option<WireVertex> = None;
        if self.edges.is_empty() {
            init = true;
            my_edge = e.clone();
            for v in my_edge.children(true) {
                self.add_vertex(v);
            }
        } else {
            init = self.closed;
            let ee = WireEdge { orientation: Orientation::Forward, ..e.clone() };
            let (mut forward, mut reverse) = (false, false);
            let mut connected = false;
            let mut copy_edge = false;
            if self.vf.is_none() || self.vl.is_none() {
                non_manifold = true;
            }
            let mut decide = |vw: &WireVertex, ve: &WireVertex, non_manifold: &mut bool, vf: Option<WireVertex>, vl: Option<WireVertex>| {
                if *non_manifold {
                    return;
                }
                if same_vertex(vf, vl) {
                    if !same_vertex(vf, Some(*vw)) {
                        *non_manifold = true;
                    }
                } else if same_vertex(vf, Some(*vw)) {
                    if ve.orientation == Orientation::Forward { reverse = true } else { forward = true }
                } else if same_vertex(vl, Some(*vw)) {
                    if ve.orientation == Orientation::Reversed { reverse = true } else { forward = true }
                } else {
                    *non_manifold = true;
                }
            };
            for ve in ee.children(true) {
                if self.vertices.iter().any(|w| w.is_same(&ve)) {
                    connected = true;
                    my_vertex = Some(ve);
                    decide(&ve, &ve, &mut non_manifold, self.vf, self.vl);
                } else {
                    // A vertex of the wire within tolerance of it.
                    if let Some(vw) = self.vertices.iter().find(|vw| {
                        let l = ve.point.distance(vw.point);
                        l < ve.tolerance || l < vw.tolerance
                    }) {
                        copy_edge = true;
                        let vw = *vw;
                        decide(&vw, &ve, &mut non_manifold, self.vf, self.vl);
                    }
                    if copy_edge {
                        connected = true;
                    }
                }
            }
            if !connected {
                return Err(Error::Failed("BRepLib_MakeWire: disconnected wire".into()));
            }
            if !copy_edge {
                my_edge = ee;
                for v in my_edge.children(true) {
                    self.add_vertex(v);
                }
            } else {
                // A copy of the edge (`EmptyCopied`) sharing the near vertices.
                my_edge = WireEdge { key: keys.fresh(), ..ee.clone() };
                for (slot, ve) in ee.vertices.iter().enumerate() {
                    let Some(ve) = *ve else { continue };
                    let near = self.vertices.iter().find(|vw| {
                        let l = ve.point.distance(vw.point);
                        l < ve.tolerance || l < vw.tolerance
                    });
                    match near.copied() {
                        Some(vw) => {
                            let l = ve.point.distance(vw.point);
                            let (tol_w, tol_e) = (vw.tolerance, ve.tolerance);
                            let mut maxtol = 0.5 * (tol_w + tol_e + l);
                            let (cw, ce);
                            if maxtol > tol_w && maxtol > tol_e {
                                cw = (maxtol - tol_e) / l;
                                ce = 1.0 - cw;
                            } else if maxtol > tol_w {
                                maxtol = tol_e;
                                cw = 0.0;
                                ce = 1.0;
                            } else {
                                maxtol = tol_w;
                                cw = 1.0;
                                ce = 0.0;
                            }
                            let pc = Point::new(
                                cw * vw.point.x + ce * ve.point.x,
                                cw * vw.point.y + ce * ve.point.y,
                                cw * vw.point.z + ce * ve.point.z,
                            );
                            self.update_vertex(vw.key, pc, maxtol);
                            let vw = WireVertex { point: pc, tolerance: maxtol, ..vw };
                            my_vertex = Some(vw.oriented(ve.orientation));
                            my_edge.vertices[slot] = Some(vw.oriented(ve.orientation));
                        }
                        None => {
                            self.add_vertex(ve);
                        }
                    }
                }
            }
            // The orientation: as given where undecided.
            if (forward == reverse && e.orientation == Orientation::Reversed) || (reverse && !forward) {
                my_edge = my_edge.reversed();
            }
        }
        self.edges.push(my_edge.clone());
        self.closed = false;
        if init {
            let wire = EdgeWire { edges: self.edges.clone(), closed: false };
            (self.vf, self.vl) = wire.ends();
        } else {
            if !non_manifold {
                let (v1, v2) = my_edge.ends(false);
                let v_ref = if same_vertex(v1, my_vertex) {
                    v2
                } else if same_vertex(v2, my_vertex) {
                    v1
                } else {
                    non_manifold = true;
                    None
                };
                if !same_vertex(self.vf, self.vl) {
                    if same_vertex(self.vf, my_vertex) {
                        self.vf = v_ref;
                    } else if same_vertex(self.vl, my_vertex) {
                        self.vl = v_ref;
                    } else {
                        non_manifold = true;
                    }
                }
            }
            if non_manifold {
                self.vf = None;
                self.vl = None;
            }
        }
        if let (Some(a), Some(b)) = (self.vf, self.vl) {
            if a.is_same(&b) {
                self.closed = true;
            }
        }
        Ok(())
    }

    fn wire(self) -> EdgeWire {
        EdgeWire { edges: self.edges, closed: self.closed }
    }
}

/// `BRepFill::InsertACR(wire, ACRcuts, prec)`: the wire with its edges cut at the reduced
/// curvilinear abscissae `cuts` (those within `prec` of an edge's ends left out).
#[doc(alias = "BRepFill::InsertACR")]
pub(crate) fn insert_acr(wire: &EdgeWire, cuts: &[f64], prec: f64, length: &dyn Fn(&GeomAdaptor) -> f64, keys: &mut ShapeKeys) -> Result<EdgeWire> {
    let acr = compute_acr(wire, length);
    let mut made = MakeWire::default();
    let mut t1 = 0.0;
    for (n, (edge, vertex)) in wire.explored_edges().into_iter().enumerate() {
        let t0 = t1;
        t1 = acr[n + 1];
        let paradec: Vec<f64> = cuts.iter().copied().filter(|&c| t0 + prec < c && c < t1 - prec).collect();
        if paradec.is_empty() || edge.degenerated {
            made.add(edge, keys)?;
        } else {
            // Cut in the wire's direction of parsing (the wire is forward; the cuts of a
            // reversed edge counted from its other end).
            let so = same_vertex(vertex, edge.first_vertex(false));
            let ndec = paradec.len();
            let sr: Vec<f64> = if edge.orientation == Orientation::Forward {
                paradec
            } else {
                (0..ndec).map(|j| t0 + t1 - paradec[ndec - 1 - j]).collect()
            };
            for piece in trim_edge(edge, &sr, t0, t1, so, keys)? {
                made.add(&piece, keys)?;
            }
        }
    }
    Ok(made.wire())
}

/// `BRepLProp_CLProps` (`GeomLProp_CLPropsBase` over a `BRepAdaptor_Curve`): a curve's
/// point and derivatives at a parameter, evaluated to the order asked for.
#[doc(alias = "BRepLProp_CLProps")]
struct CurveProps<'a> {
    curve: &'a GeomAdaptor,
    u: f64,
    order: usize,
    lin_tol: f64,
    point: Point,
    derivatives: [Vector; 3],
    /// `myTangentStatus`: `None` undecided; `Some(None)` undefined; `Some(Some(order))`
    /// defined by the first derivative not null.
    tangent: Option<Option<usize>>,
}

impl<'a> CurveProps<'a> {
    fn new(curve: &'a GeomAdaptor, u: f64, order: usize, lin_tol: f64) -> CurveProps<'a> {
        let mut props = CurveProps { curve, u, order, lin_tol, point: Point::ORIGIN, derivatives: [Vector::ZERO; 3], tangent: None };
        props.evaluate();
        props
    }

    /// `EvalDerivatives`.
    fn evaluate(&mut self) {
        let c = self.curve;
        match self.order {
            0 => self.point = c.value(self.u),
            1 => (self.point, self.derivatives[0]) = c.d1(self.u),
            2 => (self.point, self.derivatives[0], self.derivatives[1]) = c.d2(self.u),
            _ => (self.point, self.derivatives[0], self.derivatives[1], self.derivatives[2]) = c.d3(self.u),
        }
    }

    /// `D1()`, `D2()`, `D3()` (`EnsureDeriv`).
    fn derivative(&mut self, k: usize) -> Vector {
        if self.order < k {
            self.order = k;
            self.evaluate();
        }
        self.derivatives[k - 1]
    }

    /// `IsTangentDefined()`.
    fn is_tangent_defined(&mut self) -> bool {
        if let Some(status) = self.tangent {
            return status.is_some();
        }
        let tol_sq = self.lin_tol * self.lin_tol;
        // (`myCN` is 4: the first three derivatives are tried.)
        for k in 1..=3 {
            let v = self.derivative(k);
            if v.length_squared() > tol_sq {
                self.tangent = Some(Some(k));
                return true;
            }
        }
        self.tangent = Some(None);
        false
    }

    /// `Tangent(D)` once the tangent is defined (`ComputeTangent`): the first derivative not
    /// null, turned along the curve.
    fn tangent_direction(&mut self) -> Vector {
        let order = self.tangent.flatten().unwrap_or(1);
        if order == 1 {
            return direction(self.derivatives[0]);
        }
        const DIVISION_FACTOR: f64 = 1.0e-3;
        const MIN_STEP: f64 = 1.0e-7;
        let (sup, inf) = (self.curve.last_parameter(), self.curve.first_parameter());
        let du = if sup >= f64::MAX || inf <= -f64::MAX { 0.0 } else { sup - inf };
        let delta = (du * DIVISION_FACTOR).max(MIN_STEP);
        let mut v = self.derivatives[order - 1];
        let other = if self.u - inf < delta { self.u + delta } else { self.u - delta };
        let p1 = self.curve.value(self.u.min(other));
        let p2 = self.curve.value(self.u.max(other));
        if v.dot(p2 - p1) < 0.0 {
            v = -v;
        }
        direction(v)
    }
}

/// The continuity order `BRepLProp::Continuity` evaluates a curve to: 3 from C3, 2 at C2, 1
/// at C1, none otherwise.
fn props_order(continuity: Continuity) -> usize {
    match continuity {
        Continuity::C3 | Continuity::CN => 3,
        Continuity::C2 => 2,
        Continuity::C1 => 1,
        _ => 0,
    }
}

/// `BRepLProp::Continuity(C1, C2, u1, u2, tl, ta)`: the continuity where edge `e1` at
/// `u1` meets `e2` at `u2` (their derivatives along the edges' orientations). Fails where
/// the points are farther apart than `tl` ("Courbes non jointives").
#[doc(alias = "BRepLProp::Continuity")]
fn edges_continuity(e1: &WireEdge, e2: &WireEdge, u1: f64, u2: f64, tl: f64, ta: f64) -> Result<Continuity> {
    let (c1, c2) = (curve_of(e1)?, curve_of(e2)?);
    let (n1, n2) = (props_order(c1.continuity()), props_order(c2.continuity()));
    let mut clp1 = CurveProps::new(c1, u1, n1, tl);
    let mut clp2 = CurveProps::new(c2, u2, n2, tl);
    if clp1.point.distance(clp2.point) > tl {
        return Err(Error::Failed("BRepLProp::Continuity: the curves do not join".into()));
    }
    let mut cont = Continuity::C0;
    let mut done = false;
    let min = n1.min(n2);
    if min >= 1 {
        let mut d1 = clp1.derivative(1);
        let mut d2 = clp2.derivative(1);
        if e1.orientation == Orientation::Reversed {
            d1 = -d1;
        }
        if e2.orientation == Orientation::Reversed {
            d2 = -d2;
        }
        if vectors_equal(d1, d2, tl, ta) {
            cont = Continuity::C1;
        } else if clp1.is_tangent_defined() && clp2.is_tangent_defined() {
            let mut dir1 = clp1.tangent_direction();
            let mut dir2 = clp2.tangent_direction();
            if e1.orientation == Orientation::Reversed {
                dir1 = -dir1;
            }
            if e2.orientation == Orientation::Reversed {
                dir2 = -dir2;
            }
            if angle(dir1, dir2) <= ta {
                cont = Continuity::G1;
            }
            done = true;
        } else {
            done = true;
        }
    }
    if min >= 2 && !done {
        let d1 = clp1.derivative(2);
        let d2 = clp2.derivative(2);
        if vectors_equal(d1, d2, tl, ta) {
            cont = Continuity::C2;
        }
    }
    if e1.is_same(e2) && c1.is_periodic() && cont >= Continuity::G1 {
        cont = Continuity::CN;
    }
    Ok(cont)
}

/// `WireContinuity` (BRepFill_CompatibleWires.cxx): the least continuity where the wire's
/// edges meet (and where its last edge meets its first, if they touch); `CN` for a wire
/// with a degenerate edge.
fn wire_continuity(wire: &EdgeWire) -> Result<Continuity> {
    let mut cont_w = Continuity::CN;
    let edges: Vec<&WireEdge> = wire.explored_edges().into_iter().map(|(e, _)| e).collect();
    let nb_edges = edges.len();
    if edges.iter().any(|e| e.degenerated) {
        return Ok(cont_w);
    }
    let mut test = true;
    for j in 1..=nb_edges {
        let (e1, e2) = if j == nb_edges { (edges[nb_edges - 1], edges[0]) } else { (edges[j - 1], edges[j]) };
        let v1 = e1.ends(true).1.ok_or_else(|| Error::Failed("an edge without a last vertex".into()))?;
        let v2 = e2.ends(true).0.ok_or_else(|| Error::Failed("an edge without a first vertex".into()))?;
        let missing = || Error::Failed("BRep_Tool::Parameter: the vertex is not on the edge".into());
        let u1 = e1.vertex_parameter(&v1).ok_or_else(missing)?;
        let u2 = e2.vertex_parameter(&v2).ok_or_else(missing)?;
        let eps = v2.tolerance + v1.tolerance;
        if j == nb_edges {
            test = curve_of(e1)?.value(u1).distance(curve_of(e2)?.value(u2)) <= eps;
        }
        if test {
            let cont = edges_continuity(e1, e2, u1, u2, eps, ANGULAR)?;
            if cont <= cont_w {
                cont_w = cont;
            }
        }
    }
    Ok(cont_w)
}

/// Whether `PlaneOfWire` (BRepFill_CompatibleWires.cxx) finds a plane for the wire: one it
/// lies on (`BRepLib_FindSurface`), or one across the axes of inertia of its edges unless
/// two of its radii of gyration are equal (a straight wire). Kentron's own plane finder
/// and inertia (`make::plane_of_wire`, `properties::length`) stand for OCCT's here: what
/// is asked is only whether the wire is straight.
fn plane_of_wire_found(model: &Model, wire: Shape) -> bool {
    if crate::make::plane_of_wire(model, wire).is_ok() {
        return true;
    }
    let props = crate::properties::length(model, wire);
    if props.mass <= 0.0 {
        return false;
    }
    let [m1, m2, m3] = props.principal_moments();
    let radius = |m: f64| (m.max(0.0) / props.mass).sqrt();
    let (r1, r2, r3) = (radius(m1), radius(m2), radius(m3));
    let r_max = r1.max(r2).max(r3);
    let tol = CONFUSION;
    let equal = |a: f64, b: f64| (r_max - a).abs() < tol && (r_max - b).abs() < tol;
    !(equal(r1, r2) || equal(r1, r3) || equal(r2, r3))
}

/// A section of `BRepFill_CompatibleWires` (`myWork(i)`): the wire, its identity (for
/// `IsSame` between sections) and the model's wire it was read from, while unchanged.
struct Section {
    wire: EdgeWire,
    identity: usize,
    shape: Option<Shape>,
}

/// `BRepFill_CompatibleWires::SearchOrigin()`: each open wire after the first turned (its
/// edges reversed, in reverse order) where it runs against the one before; by the
/// distances between their ends where either is straight, else by the angle between their
/// chords.
fn search_origin_open(model: &Model, sections: &mut [Section], degen1: bool, degen2: bool, keys: &mut ShapeKeys) -> Result<()> {
    let (ideb, ifin) = (if degen1 { 1 } else { 0 }, if degen2 { sections.len() - 1 } else { sections.len() });
    let v_closed = !degen1 && !degen2 && sections[ideb].identity == sections[ifin - 1].identity;
    if sections[ideb..ifin].iter().any(|s| s.wire.closed) {
        return Err(Error::Failed("BRepFill_CompatibleWires: open and closed sections (NotSameTopology)".into()));
    }
    let line = |s: &Section| s.shape.is_some_and(|w| !plane_of_wire_found(model, w));
    let (vdeb, vfin) = sections[ideb].wire.ends();
    let (mut pdeb, mut pfin) = (point_of(vdeb)?, point_of(vfin)?);
    let mut isline0 = line(&sections[ideb]);
    let first_edge = |w: &EdgeWire| w.explored_edges().first().map(|(e, _)| (*e).clone());
    let mut e0 = first_edge(&sections[ideb].wire);
    for section in sections.iter_mut().take(ifin).skip(ideb + 1) {
        let wire = &section.wire;
        let seq: Vec<WireEdge> = wire.explored_edges().into_iter().map(|(e, _)| e.clone()).collect();
        let e = seq.first().cloned();
        let (vdeb, vfin) = wire.ends();
        let isline = line(section);
        let parcours;
        if isline0 || isline {
            // Straight segments: by the distances between the ends.
            let (p1, p2) = (point_of(vdeb)?, point_of(vfin)?);
            let dist1 = pdeb.distance(p1) + pfin.distance(p2);
            let dist2 = pdeb.distance(p2) + pfin.distance(p1);
            parcours = dist2 >= dist1;
        } else {
            let (p1, mut p2, p1o, mut p2o) = (point_of(vdeb)?, point_of(vfin)?, pdeb, pfin);
            if p1.distance(p2) <= CONFUSION || p1o.distance(p2o) <= CONFUSION {
                let missing = || Error::Failed("a wire without edges".into());
                let c0 = curve_of(e0.as_ref().ok_or_else(missing)?)?;
                let c = curve_of(e.as_ref().ok_or_else(missing)?)?;
                p2o = c0.value(c0.first_parameter() + CONFUSION);
                p2 = c.value(c.first_parameter() + CONFUSION);
            }
            let straight = vector_angle(p2o - p1o, p2 - p1)?;
            parcours = straight < FRAC_PI_2;
        }
        let edges: Vec<WireEdge> = if parcours { seq } else { seq.iter().rev().map(|e| e.reversed()).collect() };
        *section = Section { wire: EdgeWire { edges, closed: false }, identity: keys.fresh(), shape: None };
        if parcours {
            pdeb = point_of(vdeb)?;
            pfin = point_of(vfin)?;
        } else {
            pfin = point_of(vdeb)?;
            pdeb = point_of(vfin)?;
        }
        isline0 = isline;
        e0 = e;
    }
    if v_closed {
        let last = sections.len() - 1;
        sections[last] = Section { wire: sections[0].wire.clone(), identity: sections[0].identity, shape: sections[0].shape };
    }
    Ok(())
}

/// `BRepFill_CompatibleWires::SameNumberByACR(report)`: every wire's edges cut at the
/// reduced curvilinear abscissae of all the wires' vertices, so that they have as many.
/// Fails ("Failed" status) where they still do not.
fn same_number_by_acr(sections: &mut [Section], report: bool, degen1: bool, degen2: bool, length: &dyn Fn(&GeomAdaptor) -> f64, keys: &mut ShapeKeys) -> Result<()> {
    let (ideb, ifin) = (if degen1 { 1 } else { 0 }, if degen2 { sections.len() - 1 } else { sections.len() });
    let v_closed = !degen1 && !degen2 && sections[ideb].identity == sections[ifin - 1].identity;
    let nb_sects = sections.len();
    let counts: Vec<usize> = sections.iter().map(|s| s.wire.explore().len()).collect();
    let (nbmin, nbmax) = (*counts.iter().min().unwrap_or(&0), *counts.iter().max().unwrap_or(&0));
    if nbmax > 1 && (report || nbmin < nbmax) {
        // The table of cuts (1-based, as OCCT's).
        let mut nbdec = (nbmax - 1) * nb_sects + 1;
        let mut dec = vec![0.0; nbdec + 2];
        dec[2] = 1.0;
        let mut wire_len = vec![0.0; nb_sects];
        for (i, section) in sections.iter().enumerate() {
            let acr = compute_acr(&section.wire, length);
            wire_len[i] = acr[0];
            for &a in &acr[1..acr.len().saturating_sub(1)] {
                let mut k = 1;
                while dec[k] < a {
                    k += 1;
                    if k > nbdec {
                        break;
                    }
                }
                if k > 1 && k <= nbdec && dec[k - 1] < a && a < dec[k] {
                    for l in (k..nbdec).rev() {
                        dec[l + 1] = dec[l];
                    }
                    dec[k] = a;
                }
            }
        }
        let mut k = 1;
        while dec[k] < 1.0 {
            k += 1;
            if k > nbdec {
                break;
            }
        }
        nbdec = k - 1;
        let dec2: Vec<f64> = dec[1..=nbdec].to_vec();
        // Cuts leaving pieces too short on every wire are dropped.
        let dec3: Vec<f64> = (0..nbdec)
            .filter(|&k| {
                let knot1 = dec2[k];
                let knot2 = if k == nbdec - 1 { 1.0 } else { dec2[k + 1] };
                (0..nb_sects).any(|i| (knot2 - knot1) * wire_len[i] > CONFUSION)
            })
            .map(|k| dec2[k])
            .collect();
        for (i, section) in sections.iter_mut().enumerate() {
            let mut tol = CONFUSION;
            if wire_len[i] > RESOLUTION {
                tol /= wire_len[i];
            }
            let wire = insert_acr(&section.wire, &dec3, tol, length, keys)?;
            *section = Section { wire, identity: keys.fresh(), shape: None };
        }
    }
    if v_closed {
        let last = sections.len() - 1;
        sections[last] = Section { wire: sections[0].wire.clone(), identity: sections[0].identity, shape: sections[0].shape };
    }
    let counts: Vec<usize> = sections[ideb..ifin].iter().map(|s| s.wire.explore().len()).collect();
    if counts.iter().min() != counts.iter().max() {
        return Err(Error::Failed("BRepFill_CompatibleWires: the wires keep different numbers of edges".into()));
    }
    Ok(())
}

/// The guide of a pipe along the open `spine` as `BRepFill_PipeShell::Set` makes it from
/// the open `guide` (`BRepFill_CompatibleWires` on both, with `SetPercent(0.1)`, then
/// `Shape().Value(2)`): turned to run the spine's way and, where their numbers of edges
/// differ or both are tangent-continuous, cut so that it has as many edges as the spine.
///
/// `length` measures an edge's curve over its range: OCCT's `GCPnts_AbscissaPoint::Length`.
/// Fails where OCCT's status is not `Done` ("Uncompatible wires"), and for sections whose
/// ends meet (the paths for closed sections are not ported).
#[doc(alias = "BRepFill_CompatibleWires")]
pub(crate) fn compatible_guide(model: &Model, spine: Shape, guide: Shape, length: &dyn Fn(&GeomAdaptor) -> f64) -> Result<EdgeWire> {
    Ok(compatible_wires(model, spine, guide, length)?.1)
}

/// `BRepFill_CompatibleWires` on two open wires (`SetPercent(0.1)`, `Perform()`): both as
/// `Shape()` gives them (see [`compatible_guide`]).
pub(crate) fn compatible_wires(model: &Model, spine: Shape, guide: Shape, length: &dyn Fn(&GeomAdaptor) -> f64) -> Result<(EdgeWire, EdgeWire)> {
    let mut keys = ShapeKeys::new();
    let mut sections = Vec::new();
    for shape in [spine, guide] {
        let shape = Shape { orientation: Orientation::Forward, ..shape };
        let wire = EdgeWire::from_model(model, shape, &mut keys)?;
        sections.push(Section { wire, identity: keys.of(shape), shape: Some(shape) });
    }
    // The report by curvilinear abscissa: for different numbers of edges, or wires all at
    // least C1.
    let mut cont_s = Continuity::CN;
    let mut counts = Vec::new();
    for s in &sections {
        let cont = wire_continuity(&s.wire)?;
        if cont < cont_s {
            cont_s = cont;
        }
        counts.push(s.wire.explore().len());
    }
    let (nbmin, nbmax) = (*counts.iter().min().unwrap_or(&0), *counts.iter().max().unwrap_or(&0));
    let report = nbmax != nbmin || cont_s >= Continuity::C1;
    // Wires reduced to points at either end.
    let all_degenerated = |s: &Section| s.wire.explored_edges().iter().all(|(e, _)| e.degenerated);
    let degen1 = all_degenerated(&sections[0]);
    let degen2 = all_degenerated(&sections[sections.len() - 1]);
    let (ideb, ifin) = (if degen1 { 1 } else { 0 }, if degen2 { sections.len() - 1 } else { sections.len() });
    let (mut all_closed, mut all_open) = (true, true);
    for s in &sections[ideb.min(ifin)..ifin] {
        let (v1, v2) = s.wire.ends();
        let closed = s.wire.closed || same_vertex(v1, v2);
        all_closed &= closed;
        all_open &= !closed;
    }
    if all_closed {
        return Err(Error::Failed("BRepFill_CompatibleWires: closed sections (not ported)".into()));
    } else if all_open {
        search_origin_open(model, &mut sections, degen1, degen2, &mut keys)?;
        if report {
            same_number_by_acr(&mut sections, report, degen1, degen2, length, &mut keys)?;
        }
    } else {
        return Err(Error::Failed("BRepFill_CompatibleWires: open and closed sections (NotSameTopology)".into()));
    }
    let guide = sections.swap_remove(1).wire;
    Ok((sections.swap_remove(0).wire, guide))
}

/// Where the point is nearest a wire (`BRepExtrema_SolutionElem` on the wire's side): at a
/// vertex, or on an edge (its index in [`EdgeWire::edges`]) at a parameter.
#[derive(Debug, Clone, Copy)]
enum Support {
    Vertex(WireVertex),
    Edge(usize, f64),
}

#[derive(Debug, Clone, Copy)]
struct Solution {
    distance: f64,
    point: Point,
    support: Support,
}

/// `Extrema_ExtPC` over an edge's curve (`BRepExtrema_ExtPC(V, E)`): `None` where not done.
/// Lines and circles take `Extrema_ExtPElC`'s solutions; other curves the path for any curve
/// (OCCT's own paths for B-splines, Bézier curves and the other conics are not ported).
fn edge_extrema(p: Point, edge: &WireEdge) -> Option<Vec<ExtPCPoint>> {
    let curve = edge.curve.as_ref()?;
    let tol = edge.tolerance.min(CONFUSION);
    let tolf = curve.resolution(tol).max(P_CONFUSION);
    let (uinf, usup) = (edge.first, edge.last);
    // `mytolu`.
    let tolu = curve.resolution(CONFUSION);
    let t3d = CONFUSION;
    // `Extrema_ExtPElC`: (parameter, point, square distance, minimum), or not done.
    let elementary: Option<Vec<(f64, Point, f64, bool)>> = match (curve.kind(), curve.curve()) {
        (CurveKind::Line, Curve::Line(l)) => {
            let v1 = Vector::from(l.direction);
            let dist = v1.dot(p - l.origin);
            if dist >= uinf - t3d && dist <= usup + t3d {
                let q = l.origin + v1 * dist;
                Some(vec![(dist, q, p.distance_squared(q), true)])
            } else {
                None
            }
        }
        (CurveKind::Circle, Curve::Circle(c)) => {
            let o = c.center();
            let axe = Vector::from(c.frame().z());
            let pp = p + axe * -((p - o).dot(axe));
            let opp = pp - o;
            if opp.length() < t3d {
                None
            } else {
                // `gp_Dir::AngleWithRef` of the X axis to `OPp` about the axis.
                let x = Vector::from(c.frame().x());
                let d = direction(opp);
                let cross = x.cross(d);
                let cos = x.dot(d);
                let sin = cross.length();
                #[allow(clippy::approx_constant)]
                let mut ang = if cos > -0.70710678118655 && cos < 0.70710678118655 {
                    cos.acos()
                } else if cos < 0.0 {
                    PI - sin.asin()
                } else {
                    sin.asin()
                };
                if cross.dot(direction(axe)) < 0.0 {
                    ang = -ang;
                }
                let mut usol = [ang, 0.0];
                if usol[0] + PI < ANGULAR {
                    usol[0] = -PI;
                } else if usol[0] - PI > -ANGULAR {
                    usol[0] = PI;
                }
                usol[1] = usol[0] + PI;
                let r = c.radius();
                let tolu_c = if r > RESOLUTION { t3d / r } else { INFINITE };
                let (_, u0) = adjust_periodic(uinf, uinf + 2.0 * PI, tolu_c, uinf, usol[0]);
                let (_, u1) = adjust_periodic(uinf, uinf + 2.0 * PI, tolu_c, uinf, usol[1]);
                usol = [u0, u1];
                for u in &mut usol {
                    let gap = *u - 2.0 * PI - uinf;
                    if gap < tolu_c && gap > -tolu_c {
                        *u = uinf;
                    }
                }
                let mut out = Vec::new();
                for (k, &us) in usol.iter().enumerate() {
                    if uinf - us < tolu_c && us - usup < tolu_c {
                        let cu = curve.curve().point(us);
                        out.push((us, cu, cu.distance_squared(p), k == 0));
                    }
                }
                Some(out)
            }
        }
        _ => return ext_pc_other(curve, p, uinf, usup, tolf),
    };
    let elementary = elementary?;
    let mut out = Vec::new();
    for (u, q, sq, is_min) in elementary {
        let u = if curve.is_periodic() { in_period(u, uinf, uinf + curve.period()) } else { u };
        if u >= uinf - tolu && u <= usup + tolu {
            out.push(ExtPCPoint { parameter: u, point: q, square_distance: sq, is_min });
        }
    }
    Some(out)
}

/// `isOnBoundary` (BRepExtrema_DistanceSS.cxx): whether the point at `t` on the edge is at
/// one of its vertices.
fn is_on_boundary(edge: &WireEdge, point: Point, t: f64, ptol: f64) -> bool {
    edge.children(true).iter().any(|v| edge.vertex_parameter(v).is_some_and(|vp| (vp - t).abs() < ptol) && v.point.distance(point) < v.tolerance)
}

/// `BRepExtrema_DistShapeShape(V, W)` for a vertex at `p` and a wire, as far as it goes
/// for those (`DistanceVertVert`, then `DistanceMapMap` of the vertex and the edges by
/// `BRepExtrema_DistanceSS`): the solutions, none where it is not done. The edges are
/// tried in their order in the wire, without the filter by bounding boxes (which leaves
/// out none that could be nearest) and on one thread, so solutions equally near come in
/// that order.
fn vertex_wire_distance(p: Point, wire: &EdgeWire) -> Vec<Solution> {
    let eps = CONFUSION;
    // `TopExp::MapShapes(W, TopAbs_VERTEX / TopAbs_EDGE)`.
    let mut vertices: Vec<WireVertex> = Vec::new();
    let mut edges: Vec<usize> = Vec::new();
    for (i, e) in wire.edges.iter().enumerate() {
        if !edges.iter().any(|&j| wire.edges[j].is_same(e)) {
            edges.push(i);
        }
        for v in e.children(true) {
            if !vertices.iter().any(|w| w.is_same(&v)) {
                vertices.push(v);
            }
        }
    }
    let mut dist_ref = match vertices.first() {
        Some(v) => p.distance(v.point),
        None => 1.0e30,
    };
    let mut solutions: Vec<Solution> = Vec::new();
    let merge = |task: Vec<Solution>, task_dist: f64, solutions: &mut Vec<Solution>, dist_ref: &mut f64| {
        if task_dist < *dist_ref - eps {
            *solutions = task;
            *dist_ref = task_dist;
        } else if (task_dist - *dist_ref).abs() < eps {
            solutions.extend(task);
            if *dist_ref > task_dist {
                *dist_ref = task_dist;
            }
        }
    };
    // `DistanceVertVert`.
    let mut task = Vec::new();
    let mut task_dist = dist_ref;
    for v in &vertices {
        let d = p.distance(v.point);
        let sol = Solution { distance: d, point: v.point, support: Support::Vertex(*v) };
        if d < task_dist - eps {
            task.clear();
            task.push(sol);
            task_dist = d;
        } else if (d - task_dist).abs() < eps {
            task.push(sol);
            if task_dist > d {
                task_dist = d;
            }
        }
    }
    merge(task, task_dist, &mut solutions, &mut dist_ref);
    // `DistanceMapMap(V, edges)`.
    let mut task: Vec<Solution> = Vec::new();
    let mut task_dist = dist_ref;
    for &i in &edges {
        let edge = &wire.edges[i];
        if edge.degenerated {
            continue;
        }
        // `BRepExtrema_DistanceSS::Perform(V, E)`.
        let mut dst_ref = task_dist;
        let mut found: Vec<Solution> = Vec::new();
        let extrema = edge_extrema(p, edge).unwrap_or_default();
        if !extrema.is_empty() {
            let dstmin = extrema.iter().map(|x| x.square_distance).fold(extrema[0].square_distance, f64::min).sqrt();
            if dstmin < dst_ref - eps || (dstmin - dst_ref).abs() < eps {
                for x in &extrema {
                    if (dstmin - x.square_distance.sqrt()).abs() < eps {
                        let pt = x.point;
                        if found.iter().all(|s| s.point.distance(pt) > CONFUSION) && !is_on_boundary(edge, pt, x.parameter, P_CONFUSION) {
                            if dst_ref > dstmin {
                                dst_ref = dstmin;
                            }
                            found.push(Solution { distance: dstmin, point: pt, support: Support::Edge(i, x.parameter) });
                        }
                    }
                }
            }
        }
        if !found.is_empty() {
            if dst_ref < task_dist - eps {
                task = found;
                task_dist = dst_ref;
            } else if (dst_ref - task_dist).abs() < eps {
                task.extend(found);
                if task_dist > dst_ref {
                    task_dist = dst_ref;
                }
            }
        }
    }
    merge(task, task_dist, &mut solutions, &mut dist_ref);
    // The solutions farther than the distance dropped (as OCCT's loop does it, passing
    // over the one after each dropped).
    let mut i = 0;
    while i < solutions.len() {
        if solutions[i].distance > dist_ref + eps {
            solutions.remove(i);
        }
        i += 1;
    }
    solutions
}

/// `MakeWire` (BRepFill.cxx): the edges from `rangdeb` (1-based) on, round to it, or the
/// other way round reversed; a closed wire.
fn make_closed_wire(edges: &[WireEdge], rangdeb: usize, forward: bool) -> EdgeWire {
    let n = edges.len();
    let order: Vec<WireEdge> = if forward {
        (rangdeb..=n).chain(1..rangdeb).map(|r| edges[r - 1].clone()).collect()
    } else {
        (1..=rangdeb).rev().chain((rangdeb + 1..=n).rev()).map(|r| edges[r - 1].reversed()).collect()
    };
    EdgeWire { edges: order, closed: true }
}

/// `BRepFill::SearchOrigin(W, P, Dir, Tol)`: the closed wire `wire` (its shape taken
/// forward) started at the point nearest `p` (at a vertex within `tolerance`, else at a new
/// one cutting the edge there) and run the way `dir` points there.
#[doc(alias = "BRepFill::SearchOrigin")]
pub(crate) fn search_origin(model: &Model, wire: Shape, p: Point, dir: Vector, tolerance: f64) -> Result<EdgeWire> {
    let mut keys = ShapeKeys::new();
    let wire = EdgeWire::from_model(model, Shape { orientation: Orientation::Forward, ..wire }, &mut keys)?;
    search_origin_of(&wire, p, dir, tolerance, &mut keys)
}

/// [`search_origin`] of a wire read already (forward).
pub(crate) fn search_origin_of(wire: &EdgeWire, p: Point, dir: Vector, tolerance: f64, keys: &mut ShapeKeys) -> Result<EdgeWire> {
    if !wire.closed {
        return Err(Error::InvalidInput("BRepFill::SearchOrigin: the wire must be closed".into()));
    }
    let mut v = new_vertex(p, tolerance, keys);
    let mut theparam = 1.0e101;
    let mut cut: Option<usize> = None;
    let solutions = vertex_wire_distance(p, wire);
    if !solutions.is_empty() {
        let mut isol = 0;
        let mut dss = p.distance(solutions[0].point);
        for (i, s) in solutions.iter().enumerate().skip(1) {
            if dss > p.distance(s.point) {
                dss = p.distance(s.point);
                isol = i;
            }
        }
        match solutions[isol].support {
            Support::Vertex(w) => v = w,
            Support::Edge(index, t) => {
                let e = &wire.edges[index];
                let (vf, vl) = e.ends(false);
                let (vf, vl) = (vf.ok_or_else(|| Error::Failed("an edge without vertices".into()))?, vl.ok_or_else(|| Error::Failed("an edge without vertices".into()))?);
                let mut dist = p.distance(vf.point);
                if dist < tolerance {
                    v = vl;
                }
                let d = p.distance(vl.point);
                if d < tolerance && d < dist {
                    v = vf;
                    dist = d;
                }
                if dist > tolerance {
                    cut = Some(index);
                    theparam = t;
                }
            }
        }
    }
    // The edges, the one cut in two, and the index (1-based) of the one to start from.
    let mut edges: Vec<WireEdge> = Vec::new();
    let mut rangdeb = 0;
    for (e, current) in wire.explored_edges() {
        if cut.is_some_and(|i| e.is_same(&wire.edges[i])) {
            let (e1, e2) = cut_edge(e, theparam, current, keys)?;
            edges.push(e1);
            edges.push(e2);
            rangdeb = edges.len();
        } else {
            edges.push(e.clone());
        }
        if cut.is_none() && same_vertex(Some(v), current) {
            rangdeb = edges.len();
        }
    }
    if rangdeb == 0 {
        rangdeb = edges.len();
    }
    // The direction of parsing.
    let e = edges.get(rangdeb - 1).ok_or_else(|| Error::InvalidInput("a wire without edges".into()))?;
    if cut.is_none() {
        theparam = e.vertex_parameter(&v).ok_or_else(|| Error::Failed("BRep_Tool::Parameter: the vertex is not on the edge".into()))?;
    }
    let (_, mut ve) = curve_of(e)?.d1(theparam);
    if e.orientation == Orientation::Reversed {
        ve = -ve;
    }
    let mut ang = vector_angle(ve, dir)?;
    if ang > PI {
        ang = 2.0 * PI - ang;
    }
    let forward = ang <= FRAC_PI_2;
    Ok(make_closed_wire(&edges, rangdeb, forward))
}

#[cfg(test)]
#[path = "guide_wire_tests.rs"]
mod tests;
