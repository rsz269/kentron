//! What every curve and surface can do: evaluation, derivatives, parameter domain.

use crate::math::{Point, Vector};

use super::curve::Curve;

/// Smoothness of a curve or surface (`G` = geometric, tangent/curvature continuity only).
#[doc(alias = "GeomAbs_Shape")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Continuity {
    C0,
    G1,
    C1,
    G2,
    C2,
    C3,
    /// Infinitely differentiable.
    CN,
}

/// A curve `t -> point`, in 3D (`Point`, `Vector`) or 2D (`Point2`, `Vector2`).
#[doc(alias = "Geom_Curve")]
#[doc(alias = "Geom2d_Curve")]
#[doc(alias = "Adaptor3d_Curve")]
pub trait ParametricCurve {
    type Point: Copy;
    type Vector: Copy;

    /// Start of the parameter domain (`-tolerance::INFINITE` for unbounded curves).
    fn first_parameter(&self) -> f64;

    /// End of the parameter domain (`tolerance::INFINITE` for unbounded curves).
    fn last_parameter(&self) -> f64;

    /// True if the start and end points coincide (within `tolerance::LINEAR`).
    fn is_closed(&self) -> bool;

    fn is_periodic(&self) -> bool;

    /// The parameter period of a periodic curve.
    fn period(&self) -> Option<f64>;

    /// The smoothness over the whole domain.
    fn continuity(&self) -> Continuity;

    /// True if the curve is `n` times continuously differentiable.
    fn is_cn(&self, n: u32) -> bool;

    #[doc(alias = "D0")]
    #[doc(alias = "Value")]
    fn point(&self, t: f64) -> Self::Point;

    /// Point and first derivative.
    fn d1(&self, t: f64) -> (Self::Point, Self::Vector);

    /// Point, first and second derivatives.
    fn d2(&self, t: f64) -> (Self::Point, Self::Vector, Self::Vector);

    /// Point, first, second and third derivatives.
    fn d3(&self, t: f64) -> (Self::Point, Self::Vector, Self::Vector, Self::Vector);

    /// The `n`-th derivative. Panics if `n` is 0.
    fn dn(&self, t: f64, n: u32) -> Self::Vector;

    /// The parameter of the point at `t` on the reversed curve.
    fn reversed_parameter(&self, t: f64) -> f64;

    /// Reverses the direction of travel.
    fn reverse(&mut self);

    /// A parameter step that moves a point along the curve by at most `tolerance` (an
    /// estimate, for converting distances into parameter tolerances).
    fn resolution(&self, tolerance: f64) -> f64;

    /// Parameters inside the domain where the curve may be less smooth than elsewhere (a
    /// B-spline's knots), for algorithms that sample it. One period's worth for a
    /// periodic curve.
    fn breaks(&self) -> Vec<f64> {
        Vec::new()
    }

    fn start_point(&self) -> Self::Point {
        self.point(self.first_parameter())
    }

    fn end_point(&self) -> Self::Point {
        self.point(self.last_parameter())
    }
}

/// Point and derivatives up to order 2 of a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceD2 {
    pub point: Point,
    pub du: Vector,
    pub dv: Vector,
    pub duu: Vector,
    pub dvv: Vector,
    pub duv: Vector,
}

/// Point and derivatives up to order 3 of a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceD3 {
    pub point: Point,
    pub du: Vector,
    pub dv: Vector,
    pub duu: Vector,
    pub dvv: Vector,
    pub duv: Vector,
    pub duuu: Vector,
    pub dvvv: Vector,
    pub duuv: Vector,
    pub duvv: Vector,
}

/// A surface `(u, v) -> point`.
#[doc(alias = "Geom_Surface")]
#[doc(alias = "Adaptor3d_Surface")]
pub trait ParametricSurface {
    /// `(u1, u2, v1, v2)`; infinite bounds are `±tolerance::INFINITE`.
    fn bounds(&self) -> (f64, f64, f64, f64);

    fn is_u_closed(&self) -> bool;
    fn is_v_closed(&self) -> bool;
    fn is_u_periodic(&self) -> bool;
    fn is_v_periodic(&self) -> bool;
    fn u_period(&self) -> Option<f64>;
    fn v_period(&self) -> Option<f64>;

    fn continuity(&self) -> Continuity;

    /// True if the surface is `n` times continuously differentiable in u.
    #[doc(alias = "IsCNu")]
    fn is_cn_u(&self, n: u32) -> bool;

    #[doc(alias = "IsCNv")]
    fn is_cn_v(&self, n: u32) -> bool;

    #[doc(alias = "D0")]
    #[doc(alias = "Value")]
    fn point(&self, u: f64, v: f64) -> Point;

    /// Point and first partial derivatives `(point, du, dv)`.
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector);

    fn d2(&self, u: f64, v: f64) -> SurfaceD2;

    fn d3(&self, u: f64, v: f64) -> SurfaceD3;

    /// The partial derivative of order `nu` in u and `nv` in v. Panics if both are 0.
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector;

    /// The curve of constant `u` (parametrized by v).
    fn u_iso(&self, u: f64) -> Curve;

    /// The curve of constant `v` (parametrized by u).
    fn v_iso(&self, v: f64) -> Curve;

    fn u_reversed_parameter(&self, u: f64) -> f64;
    fn v_reversed_parameter(&self, v: f64) -> f64;

    /// Reverses the u direction (and so the normal).
    fn u_reverse(&mut self);

    /// Reverses the v direction (and so the normal).
    fn v_reverse(&mut self);

    /// Parameter steps `(du, dv)` that move a point by at most `tolerance`.
    fn resolution(&self, tolerance: f64) -> (f64, f64);
}


