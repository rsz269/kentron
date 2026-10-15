//! Bezier curves (2D and 3D) and surfaces, polynomial or rational.
//!
//! Control points are indexed from 0. Weights are optional: a curve or surface without
//! weights is polynomial ("non-rational"); setting any weight makes it rational.

use crate::math::{
    AffineTransform, AffineTransformable, Point, Point2, Transform, Transform2, Transformable,
    Transformable2,
};

use super::control::ControlPoint;
use super::{Continuity, ParametricCurve};

/// A Bezier curve on the parameter range [0, 1]: degree = number of poles - 1. The type
/// parameter is the point type: [`Point`] (3D, the default) or [`Point2`]
/// ([`BezierCurve2`](super::curve2::BezierCurve2)).
#[doc(alias = "Geom_BezierCurve")]
#[derive(Debug, Clone, PartialEq)]
pub struct BezierCurve<P = Point> {
    poles: Vec<P>,
    weights: Option<Vec<f64>>,
    /// The curve as a B-spline, what it is evaluated as (made on first use).
    bspline: super::patches::Shared<super::bspline::BSplineCurve<P>>,
}

impl<P: ControlPoint> BezierCurve<P> {
    /// The largest supported degree.
    pub const MAX_DEGREE: usize = 25;

    /// A polynomial curve. Panics unless there are 2 to `MAX_DEGREE + 1` poles.
    pub fn new(poles: Vec<P>) -> Self {
        Self::try_new(poles).expect("a Bezier curve needs 2 to MAX_DEGREE + 1 poles")
    }

    pub fn try_new(poles: Vec<P>) -> Option<Self> {
        (2..=Self::MAX_DEGREE + 1).contains(&poles.len()).then_some(Self { poles, weights: None, bspline: Default::default() })
    }

    /// A rational curve. Panics unless the pole count is valid and there is one positive
    /// weight per pole.
    pub fn rational(poles: Vec<P>, weights: Vec<f64>) -> Self {
        Self::try_rational(poles, weights).expect("invalid rational Bezier curve")
    }

    pub fn try_rational(poles: Vec<P>, weights: Vec<f64>) -> Option<Self> {
        let valid = weights.len() == poles.len() && weights.iter().all(|&w| w > 0.0);
        let mut curve = Self::try_new(poles)?;
        valid.then(|| {
            curve.weights = Some(weights);
            curve
        })
    }

    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    pub fn poles(&self) -> &[P] {
        &self.poles
    }

    pub fn pole(&self, i: usize) -> P {
        self.poles[i]
    }

    /// The weights of a rational curve; `None` for a polynomial one.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// The weight of pole `i` (1 for a polynomial curve).
    pub fn weight(&self, i: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |w| w[i])
    }

    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    pub fn set_pole(&mut self, i: usize, p: P) {
        self.bspline = Default::default();
        self.poles[i] = p;
    }

    /// Sets a pole and its weight (making the curve rational). Panics if `weight <= 0`.
    pub fn set_pole_with_weight(&mut self, i: usize, p: P, weight: f64) {
        self.bspline = Default::default();
        self.set_pole(i, p);
        self.set_weight(i, weight);
    }

    /// Sets the weight of pole `i`, making the curve rational. Panics if `weight <= 0`.
    pub fn set_weight(&mut self, i: usize, weight: f64) {
        self.bspline = Default::default();
        assert!(weight > 0.0, "non-positive weight {weight}");
        let n = self.poles.len();
        self.weights.get_or_insert_with(|| vec![1.0; n])[i] = weight;
    }

    /// Inserts a pole after pole `i` (raising the degree; the shape changes).
    pub fn insert_pole_after(&mut self, i: usize, p: P) {
        self.bspline = Default::default();
        self.insert_pole_after_with_weight(i, p, 1.0)
    }

    pub fn insert_pole_after_with_weight(&mut self, i: usize, p: P, weight: f64) {
        self.bspline = Default::default();
        self.insert_pole_at(i + 1, p, weight)
    }

    /// Inserts a pole before pole `i`.
    pub fn insert_pole_before(&mut self, i: usize, p: P) {
        self.bspline = Default::default();
        self.insert_pole_before_with_weight(i, p, 1.0)
    }

    pub fn insert_pole_before_with_weight(&mut self, i: usize, p: P, weight: f64) {
        self.bspline = Default::default();
        self.insert_pole_at(i, p, weight)
    }

    fn insert_pole_at(&mut self, i: usize, p: P, weight: f64) {
        self.bspline = Default::default();
        assert!(self.poles.len() <= Self::MAX_DEGREE, "a Bezier curve has at most MAX_DEGREE + 1 poles");
        assert!(weight > 0.0, "non-positive weight {weight}");
        self.poles.insert(i, p);
        if weight != 1.0 || self.weights.is_some() {
            let n = self.poles.len() - 1;
            self.weights.get_or_insert_with(|| vec![1.0; n]).insert(i, weight);
        }
    }

    /// Removes pole `i` (lowering the degree; the shape changes). Panics if only 2 poles
    /// are left.
    pub fn remove_pole(&mut self, i: usize) {
        self.bspline = Default::default();
        assert!(self.poles.len() > 2, "a Bezier curve keeps at least 2 poles");
        self.poles.remove(i);
        if let Some(w) = &mut self.weights {
            w.remove(i);
        }
    }

    /// Raises the degree to `degree` without changing the shape. Does nothing if the
    /// degree is already at least `degree`; panics above `MAX_DEGREE`.
    #[doc(alias = "Increase")]
    pub fn increase_degree(&mut self, degree: usize) {
        self.bspline = Default::default();
        assert!(degree <= Self::MAX_DEGREE, "degree {degree} above MAX_DEGREE");
        while self.degree() < degree {
            // Q_i = i/(n+1) P_{i-1} + (1 - i/(n+1)) P_i, in homogeneous coordinates.
            let h = self.homogeneous();
            let n = h.len();
            let elevated: Vec<[f64; 4]> = (0..=n)
                .map(|i| {
                    let a = i as f64 / n as f64;
                    let prev = if i > 0 { h[i - 1] } else { [0.0; 4] };
                    let this = if i < n { h[i] } else { [0.0; 4] };
                    std::array::from_fn(|k| a * prev[k] + (1.0 - a) * this[k])
                })
                .collect();
            self.set_homogeneous(&elevated);
        }
    }

    /// Restricts the curve to `[u1, u2]` and reparametrizes it to [0, 1].
    pub fn segment(&mut self, u1: f64, u2: f64) {
        self.bspline = Default::default();
        // Pole j of the segment is the blossom at (u1 n-j times, u2 j times).
        let h = self.homogeneous();
        let n = h.len() - 1;
        let segment: Vec<[f64; 4]> = (0..=n)
            .map(|j| {
                let mut level = h.clone();
                for k in 0..n {
                    let t = if k < j { u2 } else { u1 };
                    for i in 0..n - k {
                        level[i] = std::array::from_fn(|c| (1.0 - t) * level[i][c] + t * level[i + 1][c]);
                    }
                }
                level[0]
            })
            .collect();
        self.set_homogeneous(&segment);
    }

    /// The poles as (w x, w y, w z, w).
    fn homogeneous(&self) -> Vec<[f64; 4]> {
        self.poles
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let w = self.weights.as_ref().map_or(1.0, |w| w[i]);
                let c = p.coords();
                [c[0] * w, c[1] * w, c[2] * w, w]
            })
            .collect()
    }

    fn set_homogeneous(&mut self, h: &[[f64; 4]]) {
        self.bspline = Default::default();
        self.poles = h.iter().map(|q| P::from_coords([q[0] / q[3], q[1] / q[3], q[2] / q[3]])).collect();
        self.weights = self.weights.is_some().then(|| h.iter().map(|q| q[3]).collect());
    }

    /// The B-spline form the curve is evaluated as (made once, shared by clones).
    fn evaluated(&self) -> &super::bspline::BSplineCurve<P> {
        self.bspline.get(|| Some(self.to_bspline())).expect("a B-spline form")
    }

    /// The same curve as a B-spline (one span over [0, 1]).
    pub fn to_bspline(&self) -> super::bspline::BSplineCurve<P> {
        let d = self.degree();
        let curve = super::bspline::BSplineCurve::new(self.poles.clone(), vec![0.0, 1.0], vec![d + 1, d + 1], d);
        match &self.weights {
            Some(w) => curve.with_weights(w.clone()),
            None => curve,
        }
    }
}

impl<P: ControlPoint> ParametricCurve for BezierCurve<P> {
    type Point = P;
    type Vector = P::Vector;
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        1.0
    }
    fn is_closed(&self) -> bool {
        self.poles[0] == self.poles[self.poles.len() - 1]
    }
    fn is_periodic(&self) -> bool {
        false
    }
    fn period(&self) -> Option<f64> {
        None
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, t: f64) -> P {
        self.evaluated().point(t)
    }
    fn d1(&self, t: f64) -> (P, P::Vector) {
        self.evaluated().d1(t)
    }
    fn d2(&self, t: f64) -> (P, P::Vector, P::Vector) {
        self.evaluated().d2(t)
    }
    fn d3(&self, t: f64) -> (P, P::Vector, P::Vector, P::Vector) {
        self.evaluated().d3(t)
    }
    fn dn(&self, t: f64, n: u32) -> P::Vector {
        self.evaluated().dn(t, n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        1.0 - t
    }
    fn reverse(&mut self) {
        self.bspline = Default::default();
        self.poles.reverse();
        if let Some(w) = &mut self.weights {
            w.reverse();
        }
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        self.evaluated().resolution(tolerance)
    }
}

impl Transformable for BezierCurve<Point> {
    fn transformed(&self, t: &Transform) -> Self {
        Self { poles: self.poles.iter().map(|p| p.transformed(t)).collect(), bspline: Default::default(), ..self.clone() }
    }
}

impl Transformable2 for BezierCurve<Point2> {
    fn transformed(&self, t: &Transform2) -> Self {
        Self { poles: self.poles.iter().map(|p| p.transformed(t)).collect(), bspline: Default::default(), ..self.clone() }
    }
}

impl AffineTransformable for BezierCurve<Point> {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        Self { poles: self.poles.iter().map(|p| p.affine_transformed(t)).collect(), bspline: Default::default(), ..self.clone() }
    }
}

/// A Bezier surface on [0, 1] x [0, 1], defined by a grid of control points `poles[i][j]`
/// (`i` along u — a "row" —, `j` along v — a "column").
#[doc(alias = "Geom_BezierSurface")]
#[derive(Debug, Clone, PartialEq)]
pub struct BezierSurface {
    poles: Vec<Vec<Point>>,
    weights: Option<Vec<Vec<f64>>>,
    /// The surface as a B-spline, what it is evaluated as (made on first use).
    bspline: super::patches::Shared<super::bspline::BSplineSurface>,
}

impl BezierSurface {
    /// The largest supported degree in each direction.
    pub const MAX_DEGREE: usize = 25;

    /// True if the boundary curves `u = 0` and `u = 1` coincide within `tolerance` for `v`
    /// in `[v1, v2]` (taking weights into account).
    #[doc(alias = "GeomLib::IsBzUClosed")]
    pub fn is_u_closed_between(&self, v1: f64, v2: f64, tolerance: f64) -> bool {
        use super::parametric::ParametricSurface;
        (0..=20).all(|i| {
            let v = v1 + (v2 - v1) * i as f64 / 20.0;
            self.point(0.0, v).distance(self.point(1.0, v)) <= tolerance
        })
    }

    #[doc(alias = "GeomLib::IsBzVClosed")]
    pub fn is_v_closed_between(&self, u1: f64, u2: f64, tolerance: f64) -> bool {
        use super::parametric::ParametricSurface;
        (0..=20).all(|i| {
            let u = u1 + (u2 - u1) * i as f64 / 20.0;
            self.point(u, 0.0).distance(self.point(u, 1.0)) <= tolerance
        })
    }

    /// A polynomial surface. Panics unless the grid is rectangular with 2 to
    /// `MAX_DEGREE + 1` poles in each direction.
    pub fn new(poles: Vec<Vec<Point>>) -> Self {
        Self::try_new(poles).expect("invalid Bezier surface pole grid")
    }

    pub fn try_new(poles: Vec<Vec<Point>>) -> Option<Self> {
        let range = 2..=Self::MAX_DEGREE + 1;
        let valid = range.contains(&poles.len())
            && range.contains(&poles[0].len())
            && poles.iter().all(|row| row.len() == poles[0].len());
        valid.then_some(Self { poles, weights: None, bspline: Default::default() })
    }

    /// A rational surface. Panics unless the weights match the poles and are positive.
    pub fn rational(poles: Vec<Vec<Point>>, weights: Vec<Vec<f64>>) -> Self {
        let mut surface = Self::new(poles);
        assert!(
            weights.len() == surface.poles.len()
                && weights.iter().all(|row| row.len() == surface.poles[0].len())
                && weights.iter().flatten().all(|&w| w > 0.0),
            "invalid Bezier surface weights"
        );
        surface.weights = Some(weights);
        surface
    }

    pub fn u_degree(&self) -> usize {
        self.poles.len() - 1
    }

    pub fn v_degree(&self) -> usize {
        self.poles[0].len() - 1
    }

    pub fn u_pole_count(&self) -> usize {
        self.poles.len()
    }

    pub fn v_pole_count(&self) -> usize {
        self.poles[0].len()
    }

    pub fn poles(&self) -> &[Vec<Point>] {
        &self.poles
    }

    pub fn pole(&self, i: usize, j: usize) -> Point {
        self.poles[i][j]
    }

    pub fn weights(&self) -> Option<&[Vec<f64>]> {
        self.weights.as_deref()
    }

    /// The weight of pole `(i, j)` (1 for a polynomial surface).
    pub fn weight(&self, i: usize, j: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |w| w[i][j])
    }

    /// True if the weights vary along u.
    pub fn is_u_rational(&self) -> bool {
        self.weights.as_ref().is_some_and(|w| w.iter().any(|row| row.iter().zip(&w[0]).any(|(a, b)| (a - b).abs() > f64::EPSILON * a.abs().max(1.0))))
    }

    /// True if the weights vary along v.
    pub fn is_v_rational(&self) -> bool {
        self.weights.as_ref().is_some_and(|w| w.iter().any(|row| row.iter().any(|&x| (x - row[0]).abs() > f64::EPSILON * x.abs().max(1.0))))
    }

    pub fn set_pole(&mut self, i: usize, j: usize, p: Point) {
        self.bspline = Default::default();
        self.poles[i][j] = p;
    }

    pub fn set_pole_with_weight(&mut self, i: usize, j: usize, p: Point, weight: f64) {
        self.bspline = Default::default();
        self.set_pole(i, j, p);
        self.set_weight(i, j, weight);
    }

    /// Panics if `weight <= 0`.
    pub fn set_weight(&mut self, i: usize, j: usize, weight: f64) {
        self.bspline = Default::default();
        assert!(weight > 0.0, "non-positive weight {weight}");
        let (nu, nv) = (self.u_pole_count(), self.v_pole_count());
        self.weights.get_or_insert_with(|| vec![vec![1.0; nv]; nu])[i][j] = weight;
    }

    /// Replaces row `i` (the poles with u index `i`).
    pub fn set_pole_row(&mut self, i: usize, poles: &[Point]) {
        self.bspline = Default::default();
        assert_eq!(poles.len(), self.v_pole_count(), "a row of poles");
        self.poles[i] = poles.to_vec();
    }

    /// Replaces column `j` (the poles with v index `j`).
    pub fn set_pole_col(&mut self, j: usize, poles: &[Point]) {
        self.bspline = Default::default();
        assert_eq!(poles.len(), self.u_pole_count(), "a column of poles");
        for (row, &p) in self.poles.iter_mut().zip(poles) {
            row[j] = p;
        }
    }

    pub fn set_pole_col_with_weights(&mut self, j: usize, poles: &[Point], weights: &[f64]) {
        self.bspline = Default::default();
        self.set_pole_col(j, poles);
        self.set_weight_col(j, weights);
    }

    pub fn set_weight_row(&mut self, i: usize, weights: &[f64]) {
        self.bspline = Default::default();
        for (j, &w) in weights.iter().enumerate() {
            self.set_weight(i, j, w);
        }
    }

    pub fn set_weight_col(&mut self, j: usize, weights: &[f64]) {
        self.bspline = Default::default();
        for (i, &w) in weights.iter().enumerate() {
            self.set_weight(i, j, w);
        }
    }

    pub fn insert_pole_row_after(&mut self, i: usize, poles: &[Point]) {
        self.bspline = Default::default();
        self.insert_row(i + 1, poles, None);
    }

    pub fn insert_pole_row_after_with_weights(&mut self, i: usize, poles: &[Point], weights: &[f64]) {
        self.bspline = Default::default();
        self.insert_row(i + 1, poles, Some(weights));
    }

    pub fn insert_pole_row_before(&mut self, i: usize, poles: &[Point]) {
        self.bspline = Default::default();
        self.insert_row(i, poles, None);
    }

    pub fn insert_pole_col_after(&mut self, j: usize, poles: &[Point]) {
        self.bspline = Default::default();
        self.exchange_uv();
        self.insert_row(j + 1, poles, None);
        self.exchange_uv();
    }

    pub fn insert_pole_col_before(&mut self, j: usize, poles: &[Point]) {
        self.bspline = Default::default();
        self.exchange_uv();
        self.insert_row(j, poles, None);
        self.exchange_uv();
    }

    pub fn remove_pole_row(&mut self, i: usize) {
        self.bspline = Default::default();
        assert!(self.u_pole_count() > 2, "a Bezier surface keeps two rows of poles");
        self.poles.remove(i);
        if let Some(w) = &mut self.weights {
            w.remove(i);
        }
    }

    pub fn remove_pole_col(&mut self, j: usize) {
        self.bspline = Default::default();
        self.exchange_uv();
        self.remove_pole_row(j);
        self.exchange_uv();
    }

    /// Raises the degrees to at least `(u_degree, v_degree)` without changing the shape.
    #[doc(alias = "Increase")]
    pub fn increase_degree(&mut self, u_degree: usize, v_degree: usize) {
        self.bspline = Default::default();
        let mut b = super::convert::surface_to_bspline(&super::surface::Surface::Bezier(self.clone())).expect("a B-spline form");
        b.increase_degree(u_degree.max(self.u_degree()), v_degree.max(self.v_degree()));
        *self = Self::from_bspline_patch(&b);
    }

    /// Restricts the surface to `[u1, u2] x [v1, v2]`, reparametrized to [0, 1]².
    pub fn segment(&mut self, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.bspline = Default::default();
        let mut b = super::convert::surface_to_bspline(&super::surface::Surface::Bezier(self.clone())).expect("a B-spline form");
        b.segment(u1, u2, v1, v2);
        *self = Self::from_bspline_patch(&b);
    }

    /// Swaps the u and v directions (transposes the pole grid).
    pub fn exchange_uv(&mut self) {
        self.bspline = Default::default();
        let (nu, nv) = (self.u_pole_count(), self.v_pole_count());
        self.poles = (0..nv).map(|j| (0..nu).map(|i| self.poles[i][j]).collect()).collect();
        if let Some(w) = &self.weights {
            self.weights = Some((0..nv).map(|j| (0..nu).map(|i| w[i][j]).collect()).collect());
        }
    }

    /// Inserts a row of poles (with weights, or weight 1) at index `i`.
    fn insert_row(&mut self, i: usize, poles: &[Point], weights: Option<&[f64]>) {
        self.bspline = Default::default();
        assert_eq!(poles.len(), self.v_pole_count(), "a row of poles");
        assert!(self.u_pole_count() <= Self::MAX_DEGREE, "too many rows of poles");
        self.poles.insert(i, poles.to_vec());
        match (&mut self.weights, weights) {
            (Some(w), given) => w.insert(i, given.map_or_else(|| vec![1.0; poles.len()], |g| g.to_vec())),
            (None, Some(given)) => {
                let nv = poles.len();
                let mut w = vec![vec![1.0; nv]; self.poles.len()];
                w[i] = given.to_vec();
                self.weights = Some(w);
            }
            (None, None) => {}
        }
    }

    /// The Bézier surface of a one-patch B-spline surface.
    fn from_bspline_patch(b: &super::bspline::BSplineSurface) -> Self {
        let mut s = Self::new(b.poles().to_vec());
        s.weights = b.weights().map(|w| w.to_vec());
        s
    }

    /// What identifies this surface and its unchanged clones (for caches keyed by surface):
    /// that of the B-spline form they share.
    pub(crate) fn id(&self) -> super::SurfaceId {
        self.evaluated().id()
    }

    /// The B-spline form the surface is evaluated as (made once, shared by clones).
    fn evaluated(&self) -> &super::bspline::BSplineSurface {
        self.bspline.get(|| Some(self.to_bspline())).expect("a B-spline form")
    }

    /// The same surface as a B-spline (one patch over [0, 1]²).
    pub fn to_bspline(&self) -> super::bspline::BSplineSurface {
        use super::bspline::{BSplineSurface, KnotVector};
        let (du, dv) = (self.u_degree(), self.v_degree());
        let knots = |d: usize| KnotVector::new(vec![0.0, 1.0], vec![d + 1, d + 1], d);
        let surface = BSplineSurface::new(self.poles.clone(), knots(du), knots(dv));
        match &self.weights {
            Some(w) => surface.with_weights(w.clone()),
            None => surface,
        }
    }
}

impl crate::geometry::ParametricSurface for BezierSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, 1.0, 0.0, 1.0)
    }
    fn is_u_closed(&self) -> bool {
        let last = self.poles.len() - 1;
        (0..self.poles[0].len()).all(|j| self.poles[0][j] == self.poles[last][j])
    }
    fn is_v_closed(&self) -> bool {
        self.poles.iter().all(|row| row[0] == row[row.len() - 1])
    }
    fn is_u_periodic(&self) -> bool {
        false
    }
    fn is_v_periodic(&self) -> bool {
        false
    }
    fn u_period(&self) -> Option<f64> {
        None
    }
    fn v_period(&self) -> Option<f64> {
        None
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn_u(&self, _n: u32) -> bool {
        true
    }
    fn is_cn_v(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, u: f64, v: f64) -> Point {
        self.evaluated().point(u, v)
    }
    fn d1(&self, u: f64, v: f64) -> (Point, crate::math::Vector, crate::math::Vector) {
        self.evaluated().d1(u, v)
    }
    fn d2(&self, u: f64, v: f64) -> crate::geometry::SurfaceD2 {
        self.evaluated().d2(u, v)
    }
    fn d3(&self, u: f64, v: f64) -> crate::geometry::SurfaceD3 {
        self.evaluated().d3(u, v)
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> crate::math::Vector {
        self.evaluated().dn(u, v, nu, nv)
    }
    fn u_iso(&self, u: f64) -> crate::geometry::curve::Curve {
        self.evaluated().u_iso(u)
    }
    fn v_iso(&self, v: f64) -> crate::geometry::curve::Curve {
        self.evaluated().v_iso(v)
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        1.0 - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        1.0 - v
    }
    fn u_reverse(&mut self) {
        self.bspline = Default::default();
        self.poles.reverse();
        if let Some(w) = &mut self.weights {
            w.reverse();
        }
    }
    fn v_reverse(&mut self) {
        self.bspline = Default::default();
        for row in &mut self.poles {
            row.reverse();
        }
        if let Some(w) = &mut self.weights {
            for row in w {
                row.reverse();
            }
        }
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        self.evaluated().resolution(tolerance)
    }
}

impl Transformable for BezierSurface {
    fn transformed(&self, t: &Transform) -> Self {
        let poles = self.poles.iter().map(|row| row.iter().map(|p| p.transformed(t)).collect()).collect();
        BezierSurface { poles, weights: self.weights.clone(), bspline: Default::default() }
    }
}

impl AffineTransformable for BezierSurface {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        // Affine maps commute with (rational) Bezier evaluation: map the poles.
        let poles = self.poles.iter().map(|row| row.iter().map(|p| p.affine_transformed(t)).collect()).collect();
        BezierSurface { poles, weights: self.weights.clone(), bspline: Default::default() }
    }
}

#[cfg(test)]
mod curve_tests {
    use super::*;

    #[test]
    fn elevation_and_segments_keep_the_shape() {
        let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 2.0, 0.0), Point::new(3.0, 1.0, 1.0), Point::new(4.0, 0.0, 0.0)];
        let curve = BezierCurve::rational(poles, vec![1.0, 2.0, 0.5, 1.0]);
        let mut raised = curve.clone();
        raised.increase_degree(6);
        assert_eq!(raised.degree(), 6);
        let mut part = curve.clone();
        part.segment(0.2, 0.7);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!(raised.point(t).distance(curve.point(t)) < 1e-12);
            assert!(part.point(t).distance(curve.point(0.2 + 0.5 * t)) < 1e-12);
        }
    }
}
