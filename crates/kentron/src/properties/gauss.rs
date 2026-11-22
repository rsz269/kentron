//! A face's area and its moments integrated to a requested precision: the adaptive Gauss
//! integration of BRepGProp::SurfaceProperties called with an Eps (BRepGProp_Sinert,
//! BRepGProp_Gauss), with the integration orders, knots and subdivisions of BRepGProp_Face
//! and the edges of BRepGProp_Domain.

use std::f64::consts::PI;

use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Transform, Transformable, Vector, Vector2};
use crate::topology::{FaceId, Model, Orientation, Shape};

/// The least range of an interval integrated, and the least value of the measure.
const EPS_PARAM: f64 = 1e-12;
const EPS_DIM: f64 = 1e-30;
/// The points of the rule that estimates a rule's error, as a part of its points.
const ERROR_ALGEBR_RATIO: f64 = 2.0 / 3.0;
/// The most Gauss points on an interval (math::GaussPointsMax), and the most intervals.
const GAUSS_POINTS_MAX: usize = 61;
const SUBS_POWER: usize = 32;
const SM: usize = SUBS_POWER * GAUSS_POINTS_MAX + 1;
/// Precision::Infinite.
const INFINITE: f64 = 2e100;

/// The measure, first moments and second moments (`Ixx = ∫ y² + z²`, `Ixy = ∫ x y`, ...)
/// of an area about a point.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Inertia {
    pub mass: f64,
    pub ix: f64,
    pub iy: f64,
    pub iz: f64,
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

/// Sums and products, made safe for infinite values on an unbounded face (Add, AddInf,
/// Mult, MultInf).
#[derive(Clone, Copy)]
struct Arithmetic {
    infinite: bool,
}

impl Arithmetic {
    fn add(self, a: f64, b: f64) -> f64 {
        if !self.infinite {
            return a + b;
        }
        let (pos, neg) = (|x: f64| x >= 0.5 * INFINITE, |x: f64| x <= -0.5 * INFINITE);
        if pos(a) {
            return if neg(b) { 0.0 } else { INFINITE };
        }
        if pos(b) {
            return if neg(a) { 0.0 } else { INFINITE };
        }
        if neg(a) {
            return if pos(b) { 0.0 } else { -INFINITE };
        }
        if neg(b) {
            return if pos(a) { 0.0 } else { -INFINITE };
        }
        a + b
    }

    fn mult(self, a: f64, b: f64) -> f64 {
        if !self.infinite {
            return a * b;
        }
        if a == 0.0 || b == 0.0 {
            return 0.0;
        }
        let (pos, neg) = (|x: f64| x >= 0.5 * INFINITE, |x: f64| x <= -0.5 * INFINITE);
        if pos(a) {
            return if b < 0.0 { -INFINITE } else { INFINITE };
        }
        if pos(b) {
            return if a < 0.0 { -INFINITE } else { INFINITE };
        }
        if neg(a) {
            return if b < 0.0 { INFINITE } else { -INFINITE };
        }
        if neg(b) {
            return if a < 0.0 { INFINITE } else { -INFINITE };
        }
        a * b
    }

    fn add_inertia(self, from: &Inertia, to: &mut Inertia) {
        to.mass = self.add(to.mass, from.mass);
        to.ix = self.add(to.ix, from.ix);
        to.iy = self.add(to.iy, from.iy);
        to.iz = self.add(to.iz, from.iz);
        to.ixx = self.add(to.ixx, from.ixx);
        to.iyy = self.add(to.iyy, from.iyy);
        to.izz = self.add(to.izz, from.izz);
        to.ixy = self.add(to.ixy, from.ixy);
        to.ixz = self.add(to.ixz, from.ixz);
        to.iyz = self.add(to.iyz, from.iyz);
    }

    /// computeSInertiaOfElementaryPart: the element `|N| weight` at `p`, about `location`.
    fn add_element(self, p: Point, normal: Vector, location: Point, weight: f64, out: &mut Inertia) {
        let ds = self.mult(normal.length(), weight);
        let x = self.add(p.x, -location.x);
        let y = self.add(p.y, -location.y);
        let z = self.add(p.z, -location.z);
        out.mass = self.add(out.mass, ds);
        let (xds, yds, zds) = (self.mult(x, ds), self.mult(y, ds), self.mult(z, ds));
        out.ix = self.add(out.ix, xds);
        out.iy = self.add(out.iy, yds);
        out.iz = self.add(out.iz, zds);
        out.ixy = self.add(out.ixy, self.mult(x, yds));
        out.iyz = self.add(out.iyz, self.mult(y, zds));
        out.ixz = self.add(out.ixz, self.mult(x, zds));
        let (xxds, yyds, zzds) = (self.mult(x, xds), self.mult(y, yds), self.mult(z, zds));
        out.ixx = self.add(out.ixx, self.add(yyds, zzds));
        out.iyy = self.add(out.iyy, self.add(xxds, zzds));
        out.izz = self.add(out.izz, self.add(xxds, yyds));
    }
}

/// The coefficient of an integration order for a precision (BRepGProp_Face's SCoeff and
/// LCoeff, found by experiment in OCCT).
fn s_coeff(eps: f64) -> f64 {
    if eps < 0.1 { -0.15 * (1.0 + eps.log10()) + 0.75 } else { 0.75 }
}

fn l_coeff(eps: f64) -> f64 {
    if eps < 0.1 { -0.50 * (1.0 + eps.log10()) + 0.25 } else { 0.25 }
}

/// A surface by its kind, as an adaptor sees it (a trimmed surface as its basis).
enum SurfaceKind<'a> {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Bezier(usize, usize),
    BSpline(&'a crate::geometry::surface::BSplineSurface),
    Other,
}

fn surface_kind(surface: &Surface) -> SurfaceKind<'_> {
    match surface {
        Surface::Plane(_) => SurfaceKind::Plane,
        Surface::Cylinder(_) => SurfaceKind::Cylinder,
        Surface::Cone(_) => SurfaceKind::Cone,
        Surface::Sphere(_) => SurfaceKind::Sphere,
        Surface::Torus(_) => SurfaceKind::Torus,
        Surface::Bezier(b) => SurfaceKind::Bezier(b.u_degree(), b.v_degree()),
        Surface::BSpline(b) => SurfaceKind::BSpline(b),
        Surface::Trimmed(t) => surface_kind(t.basis()),
        _ => SurfaceKind::Other,
    }
}

/// A pcurve by its kind, as an adaptor sees it (a trimmed curve as its basis).
enum CurveKind<'a> {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    Bezier(usize),
    BSpline(&'a crate::geometry::curve2::BSplineCurve2),
    Other,
}

fn curve_kind(curve: &Curve2) -> CurveKind<'_> {
    match curve {
        Curve2::Line(_) => CurveKind::Line,
        Curve2::Circle(_) => CurveKind::Circle,
        Curve2::Ellipse(_) => CurveKind::Ellipse,
        Curve2::Parabola(_) => CurveKind::Parabola,
        Curve2::Hyperbola(_) => CurveKind::Hyperbola,
        Curve2::Bezier(b) => CurveKind::Bezier(b.degree()),
        Curve2::BSpline(b) => CurveKind::BSpline(b),
        Curve2::Trimmed(t) => curve_kind(t.basis()),
        Curve2::Offset(_) => CurveKind::Other,
    }
}

/// The parameter of the reversed curve at the point of parameter `t` (Geom2d_Curve's
/// ReversedParameter).
fn reversed_parameter(curve: &Curve2, t: f64) -> f64 {
    match curve {
        Curve2::Line(_) | Curve2::Parabola(_) | Curve2::Hyperbola(_) => -t,
        Curve2::Circle(_) | Curve2::Ellipse(_) => 2.0 * PI - t,
        Curve2::Bezier(_) => 1.0 - t,
        Curve2::BSpline(b) => b.first_parameter() + b.last_parameter() - t,
        Curve2::Trimmed(c) => reversed_parameter(c.basis(), t),
        Curve2::Offset(c) => reversed_parameter(c.basis(), t),
    }
}

/// A face as BRepGProp_Face: its surface, placed, over the box of its parameters.
struct Face<'a> {
    surface: &'a Surface,
    placement: Transform,
    reversed: bool,
    bounds: (f64, f64, f64, f64),
}

impl Face<'_> {
    /// The point and `Su × Sv` (reversed on a reversed face).
    fn normal(&self, u: f64, v: f64) -> (Point, Vector) {
        let (p, su, sv) = self.surface.d1(u, v);
        let t = &self.placement;
        let n = su.transformed(t).cross(sv.transformed(t));
        (p.transformed(t), if self.reversed { -n } else { n })
    }

    fn s_int_order(&self, eps: f64) -> usize {
        let (nu, nv) = match surface_kind(self.surface) {
            SurfaceKind::Plane => (1, 1),
            SurfaceKind::Cylinder | SurfaceKind::Cone => (2, 1),
            SurfaceKind::Sphere | SurfaceKind::Torus => (2, 2),
            SurfaceKind::Bezier(u, v) => (u, v),
            SurfaceKind::BSpline(b) => (b.u_degree(), b.v_degree()),
            SurfaceKind::Other => (2, 2),
        };
        ((s_coeff(eps) * ((nu + 1).max(nv + 1)) as f64).ceil() as usize).min(GAUSS_POINTS_MAX)
    }

    fn su_int_subs(&self) -> usize {
        match surface_kind(self.surface) {
            SurfaceKind::Plane | SurfaceKind::Bezier(..) | SurfaceKind::Other => 1,
            SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => 3,
            SurfaceKind::BSpline(b) => b.u_knots().len() - 1,
        }
    }

    fn sv_int_subs(&self) -> usize {
        match surface_kind(self.surface) {
            SurfaceKind::Plane | SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Bezier(..) | SurfaceKind::Other => 1,
            SurfaceKind::Sphere => 2,
            SurfaceKind::Torus => 3,
            SurfaceKind::BSpline(b) => b.v_knots().len() - 1,
        }
    }

    fn u_knots(&self) -> Vec<f64> {
        let (u1, u2, _, _) = self.bounds;
        match surface_kind(self.surface) {
            SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => vec![0.0, PI * 2.0 / 3.0, PI * 4.0 / 3.0, PI * 2.0],
            SurfaceKind::BSpline(b) => b.u_knots().to_vec(),
            _ => vec![u1, u2],
        }
    }

    /// (As OCCT, the box's bounds in u for planes, cylinders, cones and others.)
    fn v_knots(&self) -> Vec<f64> {
        let (u1, u2, _, _) = self.bounds;
        match surface_kind(self.surface) {
            SurfaceKind::Sphere => vec![-PI / 2.0, 0.0, PI / 2.0],
            SurfaceKind::Torus => vec![0.0, PI * 2.0 / 3.0, PI * 4.0 / 3.0, PI * 2.0],
            SurfaceKind::BSpline(b) => b.v_knots().to_vec(),
            _ => vec![u1, u2],
        }
    }
}

/// A pcurve loaded as BRepGProp_Face::Load(edge): reversed for a reversed edge.
struct Loaded {
    curve: Curve2,
    reversed: bool,
    first: f64,
    last: f64,
}

impl Loaded {
    fn new(curve: Curve2, a: f64, b: f64, reversed: bool) -> Self {
        let (first, last) = if reversed { (reversed_parameter(&curve, b), reversed_parameter(&curve, a)) } else { (a, b) };
        Loaded { curve, reversed, first, last }
    }

    fn d1(&self, l: f64) -> (Point2, Vector2) {
        if self.reversed {
            let (p, d) = self.curve.d1(reversed_parameter(&self.curve, l));
            (p, -d)
        } else {
            self.curve.d1(l)
        }
    }

    fn l_int_subs(&self) -> usize {
        match curve_kind(&self.curve) {
            CurveKind::Circle | CurveKind::Ellipse => 3,
            CurveKind::BSpline(b) => b.knots().len() - 1,
            _ => 1,
        }
    }

    fn l_knots(&self) -> Vec<f64> {
        match curve_kind(&self.curve) {
            CurveKind::Circle | CurveKind::Ellipse => vec![0.0, PI * 2.0 / 3.0, PI * 4.0 / 3.0, PI * 2.0],
            CurveKind::BSpline(b) if self.reversed => b.knots().iter().rev().map(|&k| reversed_parameter(&self.curve, k)).collect(),
            CurveKind::BSpline(b) => b.knots().to_vec(),
            _ => vec![self.first, self.last],
        }
    }

    fn l_int_order(&self, face: &Face, eps: f64) -> usize {
        let (a, b) = if self.reversed { (reversed_parameter(&self.curve, self.last), reversed_parameter(&self.curve, self.first)) } else { (self.first, self.last) };
        let bounds = crate::geometry::bounds::curve2_between(&self.curve, a, b, 1e-7);
        let (_, _, v1, v2) = face.bounds;
        let dv = v2 - v1;
        let r = if dv > f64::EPSILON { ((bounds.max.y - bounds.min.y) / dv).min(1.0) } else { 1.0 };
        let r_int = (face.sv_int_subs() as f64 * r).ceil() as usize;
        let l_subs = self.l_int_subs();
        let ns = (face.s_int_order(1.0) * r_int / l_subs).max(1) as f64;
        let nl: f64 = match curve_kind(&self.curve) {
            CurveKind::Line => 1.0,
            CurveKind::Circle | CurveKind::Ellipse | CurveKind::Parabola => 6.0,
            CurveKind::Hyperbola => 9.0,
            CurveKind::Bezier(d) => d as f64,
            CurveKind::BSpline(b) => b.degree() as f64,
            CurveKind::Other => 9.0,
        };
        let nl = nl.max(ns);
        let nn = if l_subs <= 4 { (l_coeff(eps) * (nl + 1.0)).ceil() } else { nl + 1.0 };
        (nn as usize).min(GAUSS_POINTS_MAX)
    }
}

/// The `n` Gauss–Legendre points and weights on [-1, 1] (math::GaussPoints, GaussWeights),
/// the roots of the Legendre polynomial by Newton's method.
fn gauss(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut points, mut weights) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n.div_ceil(2) {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                (p0, p1) = (p1, p2);
            }
            derivative = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let step = p1 / derivative;
            x -= step;
            if step.abs() <= 1e-16 {
                break;
            }
        }
        points[i] = -x;
        points[n - 1 - i] = x;
        let w = 2.0 / ((1.0 - x * x) * derivative * derivative);
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    (points, weights)
}

/// fillIntervalBounds: `a` to `b` cut at the knots between them, into `first` and `last`
/// (grown as needed, as for the knots' bisections); the number of intervals.
fn fill_interval_bounds(a: f64, b: f64, knots: &[f64], first: &mut Vec<f64>, last: &mut Vec<f64>, errors: &mut Vec<f64>, inertias: &mut Vec<Inertia>) -> usize {
    let size = knots.len().max((knots.len() - 1) * SUBS_POWER + 1);
    if size > first.len() {
        first.resize(size, 0.0);
        last.resize(size, 0.0);
        errors.clear();
        errors.resize(size, 0.0);
        inertias.resize(size, Inertia::default());
    }
    let from = knots.partition_point(|&k| k <= a);
    let to = from + knots[from..].partition_point(|&k| k < b);
    let inside = &knots[from..to];
    first[0] = a;
    first[1..=inside.len()].copy_from_slice(inside);
    last[..inside.len()].copy_from_slice(inside);
    last[inside.len()] = b;
    inside.len() + 1
}

/// maxSubs: `n` intervals bisected up to SUBS_POWER times as many.
fn max_subs(n: usize) -> usize {
    n * SUBS_POWER + 1
}

/// The first index of the largest value (std::max_element).
fn argmax(values: &[f64]) -> usize {
    let mut best = 0;
    for (i, &v) in values.iter().enumerate() {
        if v > values[best] {
            best = i;
        }
    }
    best
}

/// The area of `face` and its moments about `location`, integrated to the relative error
/// `eps` (BRepGProp_Gauss::Compute with an Eps, for BRepGProp_Sinert): an error estimated
/// and the intervals along the boundary bisected where it is largest only for `eps` in
/// (0, 0.001); the error, relative to the area, estimated for a negative `eps` or one
/// below 0.001. Returns the moments and that error (`eps`'s size where not estimated).
pub(super) fn face_area(model: &Model, face: Shape, location: Point, eps: f64) -> (Inertia, f64) {
    let is_error_calculation = 0.0 > eps || eps < 0.001;
    let is_verify_computation = 0.0 < eps && eps < 0.001;
    let epsilon = eps.abs();

    let id = FaceId::from_index(face.node.index);
    let record = model.get(id).expect("face in model");
    let surface = &model.get(record.surface).expect("face with a surface").surface;
    let the_face = Face { surface, placement: face.placement, reversed: face.orientation == Orientation::Reversed, bounds: model.face_uv_bounds(face) };
    let (bu1, bu2, bv1, bv2) = the_face.bounds;
    let arithmetic = Arithmetic { infinite: [bu1, bu2, bv1, bv2].iter().any(|x| x.abs() >= 0.5 * INFINITE) };
    let (add, mult) = (|a, b| arithmetic.add(a, b), |a, b| arithmetic.mult(a, b));

    // The edges of the face used forward, as BRepGProp_Domain explores them: those inside
    // or outside it left out but for the first.
    let mut edges: Vec<Option<Loaded>> = Vec::new();
    for &wire_use_id in &record.wires {
        let wire_use = model.get(wire_use_id).expect("wire use");
        let wire = model.get(wire_use.wire).expect("wire in model");
        for &coedge_id in &wire.coedges {
            let coedge = model.get(coedge_id).expect("coedge in model");
            let orientation = coedge.orientation.compose(wire_use.orientation);
            if !edges.is_empty() && !matches!(orientation, Orientation::Forward | Orientation::Reversed) {
                continue;
            }
            edges.push(model.effective_pcurve(coedge_id).map(|(curve, a, b)| Loaded::new(curve, a, b, orientation == Orientation::Reversed)));
        }
    }
    let natural_restriction = record.wires.is_empty();

    let mut inertia = Inertia::default();
    let (mut l1s, mut l2s, mut err_l, mut inertia_l) = (vec![0.0; SM], vec![0.0; SM], vec![0.0; SM], vec![Inertia::default(); SM]);
    let (mut u1s, mut u2s, mut err_u, mut inertia_u) = (vec![0.0; SM], vec![0.0; SM], vec![0.0; SM], vec![Inertia::default(); SM]);
    let mut err_ul = vec![0.0; SM];

    let u1 = bu1;
    let mut error_l_max = 0.0f64;
    let (mut eps_total, mut eps_l, mut eps_u) = (0.0f64, 0.0f64, 0.0f64);
    let gl_end = if is_error_calculation { 2 } else { 1 };

    let nb_u = [the_face.s_int_order(epsilon), 0];
    let nb_u = [nb_u[0], (ERROR_ALGEBR_RATIO * nb_u[0] as f64).ceil() as usize];
    let u_gauss = [gauss(nb_u[0]), gauss(nb_u[1])];
    let u_knots = the_face.u_knots();

    let mut next_edge = 0;
    loop {
        let loaded = if natural_restriction {
            None
        } else {
            if next_edge >= edges.len() {
                break;
            }
            next_edge += 1;
            match &edges[next_edge - 1] {
                Some(loaded) => Some(loaded),
                None => return (Inertia { mass: INFINITE, ..Inertia::default() }, INFINITE),
            }
        };
        let nb_l0 = match loaded {
            None => (2 * nb_u[0]).min(GAUSS_POINTS_MAX),
            Some(c) => c.l_int_order(&the_face, epsilon),
        };
        let nb_l = [nb_l0, (ERROR_ALGEBR_RATIO * nb_l0 as f64).ceil() as usize];
        let l_gauss = [gauss(nb_l[0]), gauss(nb_l[1])];
        let (l_knots, l1, l2) = match loaded {
            None => (the_face.v_knots(), bv1, bv2),
            Some(c) => (c.l_knots(), c.first, c.last),
        };
        let mut error_l = 0.0;
        let mut kl_end = 1;
        let mut jl = 0usize;

        if (l2 - l1).abs() > EPS_PARAM {
            let il_sub_end = fill_interval_bounds(l1, l2, &l_knots, &mut l1s, &mut l2s, &mut err_l, &mut inertia_l);
            if err_ul.len() < err_l.len() {
                err_ul.resize(err_l.len(), 0.0);
            }
            let l_max_subs = max_subs(il_sub_end).min(SM);
            for i in inertia_l.iter_mut() {
                i.mass = 0.0;
            }
            err_l.iter_mut().for_each(|e| *e = 0.0);
            err_ul.iter_mut().for_each(|e| *e = 0.0);
            let mut l_range = [0usize; 2];

            loop {
                jl += 1;
                if jl > il_sub_end {
                    let il = argmax(&err_l);
                    l_range = [il, jl - 1];
                    l1s[jl - 1] = (l1s[il] + l2s[il]) * 0.5;
                    l2s[jl - 1] = l2s[il];
                    l2s[il] = l1s[jl - 1];
                } else {
                    l_range[0] = jl - 1;
                }

                if jl == l_max_subs || (l2s[jl - 1] - l1s[jl - 1]).abs() < EPS_PARAM {
                    if kl_end == 1 {
                        inertia_l[jl - 1] = Inertia::default();
                        err_l[jl - 1] = 0.0;
                    } else {
                        jl -= 1;
                        eps_l = error_l;
                        eps_total = eps_l / 0.9;
                        break;
                    }
                } else {
                    for &ils in &l_range[..kl_end] {
                        let lm = 0.5 * (l2s[ils] + l1s[ils]);
                        let lr = 0.5 * (l2s[ils] - l1s[ils]);
                        let (mut cix, mut ciy, mut ciz, mut cixy, mut cixz, mut ciyz) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                        let (mut cdim, mut cixx, mut ciyy, mut cizz) = ([0.0; 2], [0.0; 2], [0.0; 2], [0.0; 2]);

                        for igl in 0..gl_end {
                            cdim[igl] = 0.0;
                            cixx[igl] = 0.0;
                            ciyy[igl] = 0.0;
                            cizz[igl] = 0.0;
                            for il in 0..nb_l[igl] {
                                let l = lm + lr * l_gauss[igl].0[il];
                                let (v, u2, dul) = match loaded {
                                    None => (l, bu2, l_gauss[igl].1[il]),
                                    Some(c) => {
                                        let (puv, vuv) = c.d1(l);
                                        let dul = vuv.y * l_gauss[igl].1[il];
                                        if dul.abs() < EPS_PARAM {
                                            continue;
                                        }
                                        (puv.y.clamp(bv1, bv2), puv.x.clamp(bu1, bu2), dul)
                                    }
                                };
                                err_ul[ils] = 0.0;
                                let mut ku_end = 1;
                                let mut ju = 0usize;
                                if (u2 - u1).abs() < EPS_PARAM {
                                    continue;
                                }
                                let iu_sub_end = fill_interval_bounds(u1, u2, &u_knots, &mut u1s, &mut u2s, &mut err_u, &mut inertia_u);
                                let u_max_subs = max_subs(iu_sub_end).min(SM);
                                for i in inertia_u.iter_mut() {
                                    i.mass = 0.0;
                                }
                                err_u.iter_mut().for_each(|e| *e = 0.0);
                                let mut error_u = 0.0;
                                let mut u_range = [0usize; 2];

                                loop {
                                    ju += 1;
                                    if ju > iu_sub_end {
                                        let iu = argmax(&err_u);
                                        u_range = [iu, ju - 1];
                                        u1s[ju - 1] = (u1s[iu] + u2s[iu]) * 0.5;
                                        u2s[ju - 1] = u2s[iu];
                                        u2s[iu] = u1s[ju - 1];
                                    } else {
                                        u_range[0] = ju - 1;
                                    }

                                    if ju == u_max_subs || (u2s[ju - 1] - u1s[ju - 1]).abs() < EPS_PARAM {
                                        if ku_end == 1 {
                                            err_u[ju - 1] = 0.0;
                                            inertia_u[ju - 1] = Inertia::default();
                                        } else {
                                            ju -= 1;
                                            eps_u = error_u;
                                            eps_total = 10.0 * eps_u * ((u2 - u1) * dul).abs();
                                            eps_l = 0.9 * eps_total;
                                            break;
                                        }
                                    } else {
                                        for &ius in &u_range[..ku_end] {
                                            let mut local = [Inertia::default(); 2];
                                            let length = gl_end - igl;
                                            let um = 0.5 * (u2s[ius] + u1s[ius]);
                                            let ur = 0.5 * (u2s[ius] - u1s[ius]);
                                            for igu in 0..length {
                                                for iu in 0..nb_u[igu] {
                                                    let w = u_gauss[igu].1[iu];
                                                    let u = um + ur * u_gauss[igu].0[iu];
                                                    let (p, n) = the_face.normal(u, v);
                                                    if igu > 0 {
                                                        local[igu].mass += w * n.length();
                                                    } else {
                                                        arithmetic.add_element(p, n, location, w, &mut local[igu]);
                                                    }
                                                }
                                            }
                                            let ui = &mut inertia_u[ius];
                                            ui.mass = mult(local[0].mass, ur);
                                            if igl > 0 {
                                                continue;
                                            }
                                            let d_mass = (local[1].mass - local[0].mass).abs();
                                            ui.ix = mult(local[0].ix, ur);
                                            ui.iy = mult(local[0].iy, ur);
                                            ui.iz = mult(local[0].iz, ur);
                                            ui.ixx = mult(local[0].ixx, ur);
                                            ui.iyy = mult(local[0].iyy, ur);
                                            ui.izz = mult(local[0].izz, ur);
                                            ui.ixy = mult(local[0].ixy, ur);
                                            ui.ixz = mult(local[0].ixz, ur);
                                            ui.iyz = mult(local[0].iyz, ur);
                                            err_u[ius] = mult(d_mass, ur);
                                        }
                                    }
                                    if ju == iu_sub_end {
                                        ku_end = 2;
                                        error_u = err_u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                                    }
                                    if !((error_u - eps_u > 0.0 && eps_u != 0.0) || ku_end == 1) {
                                        break;
                                    }
                                }

                                for iu in &inertia_u[..ju] {
                                    cdim[igl] = add(cdim[igl], mult(iu.mass, dul));
                                    cixx[igl] = add(cixx[igl], mult(iu.ixx, dul));
                                    ciyy[igl] = add(ciyy[igl], mult(iu.iyy, dul));
                                    cizz[igl] = add(cizz[igl], mult(iu.izz, dul));
                                }
                                if igl > 0 {
                                    continue;
                                }
                                err_ul[ils] = error_u * ((u2 - u1) * dul).abs();
                                for iu in &inertia_u[..ju] {
                                    cix = add(cix, mult(iu.ix, dul));
                                    ciy = add(ciy, mult(iu.iy, dul));
                                    ciz = add(ciz, mult(iu.iz, dul));
                                    cixy = add(cixy, mult(iu.ixy, dul));
                                    cixz = add(cixz, mult(iu.ixz, dul));
                                    ciyz = add(ciyz, mult(iu.iyz, dul));
                                }
                            }
                        }

                        let li = &mut inertia_l[ils];
                        li.mass = mult(cdim[0], lr);
                        li.ixx = mult(cixx[0], lr);
                        li.iyy = mult(ciyy[0], lr);
                        li.izz = mult(cizz[0], lr);
                        if gl_end == 2 {
                            let sub_dim = (cdim[1] - cdim[0]).abs();
                            err_l[ils] = add(mult(sub_dim, lr), err_ul[ils]);
                        }
                        li.ix = mult(cix, lr);
                        li.iy = mult(ciy, lr);
                        li.iz = mult(ciz, lr);
                        li.ixy = mult(cixy, lr);
                        li.ixz = mult(cixz, lr);
                        li.iyz = mult(ciyz, lr);
                    }
                }

                // The precision corrected by the measure found so far.
                if jl == il_sub_end {
                    kl_end = 2;
                    let d_dim: f64 = inertia_l[..jl].iter().map(|i| i.mass).sum();
                    let d_dim = (d_dim * epsilon).abs();
                    if d_dim > eps_total {
                        eps_total = d_dim;
                        eps_l = 0.9 * eps_total;
                    }
                }
                if kl_end == 2 {
                    error_l = err_l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                }
                if !((error_l - eps_l > 0.0 && is_verify_computation) || kl_end == 1) {
                    break;
                }
            }

            for li in &inertia_l[..jl] {
                arithmetic.add_inertia(li, &mut inertia);
            }
            error_l_max = error_l_max.max(error_l);
        }

        if natural_restriction {
            break;
        }
    }

    let mass = if inertia.mass.abs() >= EPS_DIM { inertia.mass } else { 0.0 };
    let error = if gl_end == 2 { if mass != 0.0 { error_l_max / mass.abs() } else { 0.0 } } else { epsilon };
    (inertia, error)
}
