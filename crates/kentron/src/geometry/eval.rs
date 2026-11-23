//! Evaluation of the elementary curves and surfaces, and dispatch for the `Curve`,
//! `Curve2` and `Surface` enums.

use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::{Axis, Direction, Frame, Point, Point2, Transform, Transformable, Vector, Vector2};
use crate::tolerance::{INFINITE, LINEAR, RESOLUTION};

use super::curve::{Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola, TrimmedCurve};
use super::curve2::{Circle2, Curve2, Ellipse2, Hyperbola2, Line2, OffsetCurve2, Parabola2};
use super::parametric::{Continuity, ParametricCurve, ParametricSurface, SurfaceD2, SurfaceD3};
use super::surface::{Cone, Cylinder, ExtrusionSurface, Plane, RevolutionSurface, Sphere, Surface, Torus};

/// `cos(t + n·π/2)` and `sin(t + n·π/2)`: the `n`-th derivatives of `cos t` and `sin t`.
fn cos_sin_derivative(t: f64, n: u32) -> (f64, f64) {
    cos_sin_derivative_of(t.sin_cos(), n)
}

/// [`cos_sin_derivative`] from `t.sin_cos()` (worked out once for a point and its derivatives).
fn cos_sin_derivative_of((sin, cos): (f64, f64), n: u32) -> (f64, f64) {
    match n % 4 {
        0 => (cos, sin),
        1 => (-sin, cos),
        2 => (-cos, -sin),
        _ => (sin, -cos),
    }
}

// --- enum dispatch --------------------------------------------------------------------------

macro_rules! dispatch_curve {
    ($enum:ident, $point:ty, $vector:ty, [$($variant:ident),*]) => {
        impl ParametricCurve for $enum {
            type Point = $point;
            type Vector = $vector;
            fn first_parameter(&self) -> f64 {
                match self { $($enum::$variant(c) => c.first_parameter(),)* }
            }
            fn last_parameter(&self) -> f64 {
                match self { $($enum::$variant(c) => c.last_parameter(),)* }
            }
            fn is_closed(&self) -> bool {
                match self { $($enum::$variant(c) => c.is_closed(),)* }
            }
            fn is_periodic(&self) -> bool {
                match self { $($enum::$variant(c) => c.is_periodic(),)* }
            }
            fn period(&self) -> Option<f64> {
                match self { $($enum::$variant(c) => c.period(),)* }
            }
            fn continuity(&self) -> Continuity {
                match self { $($enum::$variant(c) => c.continuity(),)* }
            }
            fn is_cn(&self, n: u32) -> bool {
                match self { $($enum::$variant(c) => c.is_cn(n),)* }
            }
            fn point(&self, t: f64) -> $point {
                match self { $($enum::$variant(c) => c.point(t),)* }
            }
            fn d1(&self, t: f64) -> ($point, $vector) {
                match self { $($enum::$variant(c) => c.d1(t),)* }
            }
            fn d2(&self, t: f64) -> ($point, $vector, $vector) {
                match self { $($enum::$variant(c) => c.d2(t),)* }
            }
            fn d3(&self, t: f64) -> ($point, $vector, $vector, $vector) {
                match self { $($enum::$variant(c) => c.d3(t),)* }
            }
            fn dn(&self, t: f64, n: u32) -> $vector {
                match self { $($enum::$variant(c) => c.dn(t, n),)* }
            }
            fn reversed_parameter(&self, t: f64) -> f64 {
                match self { $($enum::$variant(c) => c.reversed_parameter(t),)* }
            }
            fn reverse(&mut self) {
                match self { $($enum::$variant(c) => c.reverse(),)* }
            }
            fn resolution(&self, tolerance: f64) -> f64 {
                match self { $($enum::$variant(c) => c.resolution(tolerance),)* }
            }
            fn breaks(&self) -> Vec<f64> {
                match self { $($enum::$variant(c) => c.breaks(),)* }
            }
        }
    };
}

dispatch_curve!(Curve, Point, Vector, [Line, Circle, Ellipse, Hyperbola, Parabola, Bezier, BSpline, Trimmed, Offset]);
dispatch_curve!(Curve2, Point2, Vector2, [Line, Circle, Ellipse, Hyperbola, Parabola, Bezier, BSpline, Trimmed, Offset]);

macro_rules! dispatch_surface {
    ([$($variant:ident),*]) => {
        impl ParametricSurface for Surface {
            fn bounds(&self) -> (f64, f64, f64, f64) {
                match self { $(Surface::$variant(s) => s.bounds(),)* }
            }
            fn is_u_closed(&self) -> bool {
                match self { $(Surface::$variant(s) => s.is_u_closed(),)* }
            }
            fn is_v_closed(&self) -> bool {
                match self { $(Surface::$variant(s) => s.is_v_closed(),)* }
            }
            fn is_u_periodic(&self) -> bool {
                match self { $(Surface::$variant(s) => s.is_u_periodic(),)* }
            }
            fn is_v_periodic(&self) -> bool {
                match self { $(Surface::$variant(s) => s.is_v_periodic(),)* }
            }
            fn u_period(&self) -> Option<f64> {
                match self { $(Surface::$variant(s) => s.u_period(),)* }
            }
            fn v_period(&self) -> Option<f64> {
                match self { $(Surface::$variant(s) => s.v_period(),)* }
            }
            fn continuity(&self) -> Continuity {
                match self { $(Surface::$variant(s) => s.continuity(),)* }
            }
            fn is_cn_u(&self, n: u32) -> bool {
                match self { $(Surface::$variant(s) => s.is_cn_u(n),)* }
            }
            fn is_cn_v(&self, n: u32) -> bool {
                match self { $(Surface::$variant(s) => s.is_cn_v(n),)* }
            }
            fn point(&self, u: f64, v: f64) -> Point {
                match self { $(Surface::$variant(s) => s.point(u, v),)* }
            }
            fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
                match self { $(Surface::$variant(s) => s.d1(u, v),)* }
            }
            fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
                match self { $(Surface::$variant(s) => s.d2(u, v),)* }
            }
            fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
                match self { $(Surface::$variant(s) => s.d3(u, v),)* }
            }
            fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
                match self { $(Surface::$variant(s) => s.dn(u, v, nu, nv),)* }
            }
            fn u_iso(&self, u: f64) -> Curve {
                match self { $(Surface::$variant(s) => s.u_iso(u),)* }
            }
            fn v_iso(&self, v: f64) -> Curve {
                match self { $(Surface::$variant(s) => s.v_iso(v),)* }
            }
            fn u_reversed_parameter(&self, u: f64) -> f64 {
                match self { $(Surface::$variant(s) => s.u_reversed_parameter(u),)* }
            }
            fn v_reversed_parameter(&self, v: f64) -> f64 {
                match self { $(Surface::$variant(s) => s.v_reversed_parameter(v),)* }
            }
            fn u_reverse(&mut self) {
                match self { $(Surface::$variant(s) => s.u_reverse(),)* }
            }
            fn v_reverse(&mut self) {
                match self { $(Surface::$variant(s) => s.v_reverse(),)* }
            }
            fn resolution(&self, tolerance: f64) -> (f64, f64) {
                match self { $(Surface::$variant(s) => s.resolution(tolerance),)* }
            }
        }
    };
}

dispatch_surface!([Plane, Cylinder, Cone, Sphere, Torus, Bezier, BSpline, Extrusion, Revolution, Offset, Trimmed]);

// --- trimmed curves -------------------------------------------------------------------------

impl<C: ParametricCurve> TrimmedCurve<C> {
    /// Applies the trimming rules of [`TrimmedCurve::new`], returning the normalized range.
    pub(super) fn checked_range(basis: &C, first: f64, mut last: f64) -> (f64, f64) {
        match basis.period() {
            Some(period) => {
                // Shift `last` into (first, first + period].
                last -= ((last - first) / period).floor() * period;
                if last <= first + crate::tolerance::PARAMETRIC {
                    last += period;
                }
            }
            None => {
                let slack = crate::tolerance::PARAMETRIC;
                // A range reversed by at most the confusion tolerance collapses to a point.
                assert!(first <= last + LINEAR, "reversed trimming range [{first}, {last}]");
                last = last.max(first);
                assert!(
                    first >= basis.first_parameter() - slack && last <= basis.last_parameter() + slack,
                    "trimming range [{first}, {last}] outside the curve's domain"
                );
            }
        }
        (first, last)
    }
}

impl<C: ParametricCurve> ParametricCurve for TrimmedCurve<C>
where
    C::Point: PointLike,
{
    type Point = C::Point;
    type Vector = C::Vector;
    fn first_parameter(&self) -> f64 {
        self.range().0
    }
    fn last_parameter(&self) -> f64 {
        self.range().1
    }
    fn is_closed(&self) -> bool {
        self.start_point().gap(self.end_point()) <= LINEAR
    }
    /// True only if the basis is periodic and the range covers one full period.
    fn is_periodic(&self) -> bool {
        let (first, last) = self.range();
        self.basis().period().is_some_and(|p| (last - first - p).abs() <= crate::tolerance::PARAMETRIC)
    }
    fn period(&self) -> Option<f64> {
        if self.is_periodic() { self.basis().period() } else { None }
    }
    fn continuity(&self) -> Continuity {
        self.basis().continuity()
    }
    fn is_cn(&self, n: u32) -> bool {
        self.basis().is_cn(n)
    }
    fn point(&self, t: f64) -> C::Point {
        self.basis().point(t)
    }
    fn d1(&self, t: f64) -> (C::Point, C::Vector) {
        self.basis().d1(t)
    }
    fn d2(&self, t: f64) -> (C::Point, C::Vector, C::Vector) {
        self.basis().d2(t)
    }
    fn d3(&self, t: f64) -> (C::Point, C::Vector, C::Vector, C::Vector) {
        self.basis().d3(t)
    }
    fn dn(&self, t: f64, n: u32) -> C::Vector {
        self.basis().dn(t, n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        self.basis().reversed_parameter(t)
    }
    fn reverse(&mut self) {
        let (first, last) = self.range();
        let basis = self.basis_mut();
        basis.reverse();
        let (first, last) = (basis.reversed_parameter(last), basis.reversed_parameter(first));
        self.set_range_unchecked(first, last);
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        self.basis().resolution(tolerance)
    }
    fn breaks(&self) -> Vec<f64> {
        self.basis().breaks()
    }
}

/// Distance between two points of the same dimension (for closedness checks).
pub trait PointLike: Copy {
    fn gap(self, other: Self) -> f64;
}

impl PointLike for Point {
    fn gap(self, other: Self) -> f64 {
        self.distance(other)
    }
}

impl PointLike for Point2 {
    fn gap(self, other: Self) -> f64 {
        self.distance(other)
    }
}

// --- lines ----------------------------------------------------------------------------------

impl ParametricCurve for Line {
    type Point = Point;
    type Vector = Vector;
    fn first_parameter(&self) -> f64 {
        -INFINITE
    }
    fn last_parameter(&self) -> f64 {
        INFINITE
    }
    fn is_closed(&self) -> bool {
        false
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
    fn point(&self, t: f64) -> Point {
        self.origin + Vector::from(self.direction) * t
    }
    fn d1(&self, t: f64) -> (Point, Vector) {
        (self.point(t), Vector::from(self.direction))
    }
    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        (self.point(t), Vector::from(self.direction), Vector::ZERO)
    }
    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        (self.point(t), Vector::from(self.direction), Vector::ZERO, Vector::ZERO)
    }
    fn dn(&self, _t: f64, n: u32) -> Vector {
        assert!(n > 0, "derivative of order 0");
        if n == 1 { Vector::from(self.direction) } else { Vector::ZERO }
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }
    fn reverse(&mut self) {
        self.direction = -self.direction;
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        tolerance
    }
}

impl ParametricCurve for Line2 {
    type Point = Point2;
    type Vector = Vector2;
    fn first_parameter(&self) -> f64 {
        -INFINITE
    }
    fn last_parameter(&self) -> f64 {
        INFINITE
    }
    fn is_closed(&self) -> bool {
        false
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
    fn point(&self, t: f64) -> Point2 {
        self.origin + Vector2::from(self.direction) * t
    }
    fn d1(&self, t: f64) -> (Point2, Vector2) {
        (self.point(t), Vector2::from(self.direction))
    }
    fn d2(&self, t: f64) -> (Point2, Vector2, Vector2) {
        (self.point(t), Vector2::from(self.direction), Vector2::ZERO)
    }
    fn d3(&self, t: f64) -> (Point2, Vector2, Vector2, Vector2) {
        (self.point(t), Vector2::from(self.direction), Vector2::ZERO, Vector2::ZERO)
    }
    fn dn(&self, _t: f64, n: u32) -> Vector2 {
        assert!(n > 0, "derivative of order 0");
        if n == 1 { Vector2::from(self.direction) } else { Vector2::ZERO }
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }
    fn reverse(&mut self) {
        self.direction = -self.direction;
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        tolerance
    }
}

// --- circles --------------------------------------------------------------------------------

/// `t.sin_cos()` at each of `params`: every 8th afresh, and any not a mean step on from the
/// one before; those between from the last by the rotation of one step (agreeing to a few
/// units in the last place).
fn sin_cos_along(params: &[f64]) -> Vec<(f64, f64)> {
    let n = params.len();
    if n < 3 {
        return params.iter().map(|t| t.sin_cos()).collect();
    }
    let h = (params[n - 1] - params[0]) / (n - 1) as f64;
    let (sin_h, cos_h) = h.sin_cos();
    let even = |i: usize| ((params[i] - params[i - 1]) - h).abs() <= 1e-12 * (1.0 + params[i].abs());
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(n);
    for (i, t) in params.iter().enumerate() {
        out.push(if i % 8 == 0 || !even(i) {
            t.sin_cos()
        } else {
            let (s, c) = out[i - 1];
            (s * cos_h + c * sin_h, c * cos_h - s * sin_h)
        });
    }
    out
}

impl Curve {
    /// The points at `params` (as `point` at each, a circle's and an ellipse's to a few
    /// units in the last place, for fewer sines and cosines where they are evenly spaced).
    pub(crate) fn points_along(&self, params: &[f64]) -> Vec<Point> {
        match self {
            Curve::Circle(c) => sin_cos_along(params).into_iter().map(|sc| c.point_of(sc)).collect(),
            Curve::Ellipse(e) => sin_cos_along(params).into_iter().map(|sc| e.point_of(sc)).collect(),
            _ => params.iter().map(|&t| self.point(t)).collect(),
        }
    }
}

impl Curve2 {
    /// [`Curve::points_along`] in the plane.
    pub(crate) fn points_along(&self, params: &[f64]) -> Vec<Point2> {
        match self {
            Curve2::Circle(c) => sin_cos_along(params).into_iter().map(|sc| c.point_of(sc)).collect(),
            Curve2::Ellipse(e) => sin_cos_along(params).into_iter().map(|sc| e.point_of(sc)).collect(),
            _ => params.iter().map(|&t| self.point(t)).collect(),
        }
    }
}

impl Circle {
    /// The point at the parameter whose `sin_cos()` is given.
    fn point_of(&self, (sin, cos): (f64, f64)) -> Point {
        let f = self.frame();
        f.origin() + (Vector::from(f.x()) * cos + Vector::from(f.y()) * sin) * self.radius()
    }

    /// The `n`th derivative (n > 0) at the parameter whose `sin_cos()` is given.
    fn derivative_of(&self, sin_cos: (f64, f64), n: u32) -> Vector {
        let (cos, sin) = cos_sin_derivative_of(sin_cos, n);
        let f = self.frame();
        (Vector::from(f.x()) * cos + Vector::from(f.y()) * sin) * self.radius()
    }
}

impl Circle2 {
    /// The point at the parameter whose `sin_cos()` is given.
    fn point_of(&self, (sin, cos): (f64, f64)) -> Point2 {
        let f = self.frame();
        f.origin() + (Vector2::from(f.x()) * cos + Vector2::from(f.y()) * sin) * self.radius()
    }

    /// The `n`th derivative (n > 0) at the parameter whose `sin_cos()` is given.
    fn derivative_of(&self, sin_cos: (f64, f64), n: u32) -> Vector2 {
        let (cos, sin) = cos_sin_derivative_of(sin_cos, n);
        let f = self.frame();
        (Vector2::from(f.x()) * cos + Vector2::from(f.y()) * sin) * self.radius()
    }
}

impl Ellipse {
    /// The point at the parameter whose `sin_cos()` is given.
    fn point_of(&self, (sin, cos): (f64, f64)) -> Point {
        let f = self.frame();
        f.origin() + Vector::from(f.x()) * (cos * self.major_radius()) + Vector::from(f.y()) * (sin * self.minor_radius())
    }

    /// The `n`th derivative (n > 0) at the parameter whose `sin_cos()` is given.
    fn derivative_of(&self, sin_cos: (f64, f64), n: u32) -> Vector {
        let (cos, sin) = cos_sin_derivative_of(sin_cos, n);
        let f = self.frame();
        Vector::from(f.x()) * (cos * self.major_radius()) + Vector::from(f.y()) * (sin * self.minor_radius())
    }
}

impl Ellipse2 {
    /// The point at the parameter whose `sin_cos()` is given.
    fn point_of(&self, (sin, cos): (f64, f64)) -> Point2 {
        let f = self.frame();
        f.origin() + Vector2::from(f.x()) * (cos * self.major_radius()) + Vector2::from(f.y()) * (sin * self.minor_radius())
    }

    /// The `n`th derivative (n > 0) at the parameter whose `sin_cos()` is given.
    fn derivative_of(&self, sin_cos: (f64, f64), n: u32) -> Vector2 {
        let (cos, sin) = cos_sin_derivative_of(sin_cos, n);
        let f = self.frame();
        Vector2::from(f.x()) * (cos * self.major_radius()) + Vector2::from(f.y()) * (sin * self.minor_radius())
    }
}

impl ParametricCurve for Circle {
    type Point = Point;
    type Vector = Vector;
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        TAU
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, t: f64) -> Point {
        self.point_of(t.sin_cos())
    }
    // The derivatives from the same sine and cosine (one call to sin_cos for all).
    fn d1(&self, t: f64) -> (Point, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1))
    }
    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2))
    }
    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2), self.derivative_of(sin_cos, 3))
    }
    fn dn(&self, t: f64, n: u32) -> Vector {
        assert!(n > 0, "derivative of order 0");
        self.derivative_of(t.sin_cos(), n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        TAU - t
    }
    /// Flips the frame's Y axis (as OCCT does), so `t` maps to `2π - t`.
    fn reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Circle::new(frame, self.radius());
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        if self.radius() > tolerance { 2.0 * (tolerance / (2.0 * self.radius())).asin() } else { TAU }
    }
}

impl ParametricCurve for Circle2 {
    type Point = Point2;
    type Vector = Vector2;
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        TAU
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, t: f64) -> Point2 {
        self.point_of(t.sin_cos())
    }
    // The derivatives from the same sine and cosine (one call to sin_cos for all).
    fn d1(&self, t: f64) -> (Point2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1))
    }
    fn d2(&self, t: f64) -> (Point2, Vector2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2))
    }
    fn d3(&self, t: f64) -> (Point2, Vector2, Vector2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2), self.derivative_of(sin_cos, 3))
    }
    fn dn(&self, t: f64, n: u32) -> Vector2 {
        assert!(n > 0, "derivative of order 0");
        self.derivative_of(t.sin_cos(), n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        TAU - t
    }
    /// Flips the frame's Y axis, so `t` maps to `2π - t`.
    fn reverse(&mut self) {
        let f = self.frame();
        let frame = crate::math::Frame2::from_axes(f.origin(), f.x(), -f.y()).expect("orthonormal frame");
        *self = Circle2::new(frame, self.radius());
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        if self.radius() > tolerance { 2.0 * (tolerance / (2.0 * self.radius())).asin() } else { TAU }
    }
}

impl ParametricCurve for Ellipse {
    type Point = Point;
    type Vector = Vector;
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        TAU
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, t: f64) -> Point {
        self.point_of(t.sin_cos())
    }
    // The derivatives from the same sine and cosine (one call to sin_cos for all).
    fn d1(&self, t: f64) -> (Point, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1))
    }
    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2))
    }
    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2), self.derivative_of(sin_cos, 3))
    }
    fn dn(&self, t: f64, n: u32) -> Vector {
        assert!(n > 0, "derivative of order 0");
        self.derivative_of(t.sin_cos(), n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        TAU - t
    }
    /// Flips the frame's Y axis, so `t` maps to `2π - t`.
    fn reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Ellipse::new(frame, self.major_radius(), self.minor_radius());
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        let r = self.major_radius();
        if r > tolerance { 2.0 * (tolerance / (2.0 * r)).asin() } else { TAU }
    }
}

impl ParametricCurve for Ellipse2 {
    type Point = Point2;
    type Vector = Vector2;
    fn first_parameter(&self) -> f64 {
        0.0
    }
    fn last_parameter(&self) -> f64 {
        TAU
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn continuity(&self) -> Continuity {
        Continuity::CN
    }
    fn is_cn(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, t: f64) -> Point2 {
        self.point_of(t.sin_cos())
    }
    // The derivatives from the same sine and cosine (one call to sin_cos for all).
    fn d1(&self, t: f64) -> (Point2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1))
    }
    fn d2(&self, t: f64) -> (Point2, Vector2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2))
    }
    fn d3(&self, t: f64) -> (Point2, Vector2, Vector2, Vector2) {
        let sin_cos = t.sin_cos();
        (self.point_of(sin_cos), self.derivative_of(sin_cos, 1), self.derivative_of(sin_cos, 2), self.derivative_of(sin_cos, 3))
    }
    fn dn(&self, t: f64, n: u32) -> Vector2 {
        assert!(n > 0, "derivative of order 0");
        self.derivative_of(t.sin_cos(), n)
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        TAU - t
    }
    /// Flips the frame's Y axis, so `t` maps to `2π - t`.
    fn reverse(&mut self) {
        let f = self.frame();
        let frame = crate::math::Frame2::from_axes(f.origin(), f.x(), -f.y()).expect("orthonormal frame");
        *self = Ellipse2::new(frame, self.major_radius(), self.minor_radius());
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        let r = self.major_radius();
        if r > tolerance { 2.0 * (tolerance / (2.0 * r)).asin() } else { TAU }
    }
}

// --- surfaces -------------------------------------------------------------------------------

/// Derivatives up to order 3 assembled from `dn`.
fn d2_from<S: ParametricSurface + ?Sized>(s: &S, u: f64, v: f64) -> SurfaceD2 {
    SurfaceD2 {
        point: s.point(u, v),
        du: s.dn(u, v, 1, 0),
        dv: s.dn(u, v, 0, 1),
        duu: s.dn(u, v, 2, 0),
        dvv: s.dn(u, v, 0, 2),
        duv: s.dn(u, v, 1, 1),
    }
}

/// The point and derivatives up to the second from the point and `dn` (each derivative
/// order as `dn(nu, nv)`), for surfaces that work out what the derivatives share once.
fn d2_of(point: Point, dn: impl Fn(u32, u32) -> Vector) -> SurfaceD2 {
    SurfaceD2 { point, du: dn(1, 0), dv: dn(0, 1), duu: dn(2, 0), dvv: dn(0, 2), duv: dn(1, 1) }
}

/// [`d2_of`] up to the third derivatives.
fn d3_of(point: Point, dn: impl Fn(u32, u32) -> Vector) -> SurfaceD3 {
    let d2 = d2_of(point, &dn);
    SurfaceD3 { point: d2.point, du: d2.du, dv: d2.dv, duu: d2.duu, dvv: d2.dvv, duv: d2.duv, duuu: dn(3, 0), dvvv: dn(0, 3), duuv: dn(2, 1), duvv: dn(1, 2) }
}

fn d3_from<S: ParametricSurface + ?Sized>(s: &S, u: f64, v: f64) -> SurfaceD3 {
    let d2 = d2_from(s, u, v);
    SurfaceD3 {
        point: d2.point,
        du: d2.du,
        dv: d2.dv,
        duu: d2.duu,
        dvv: d2.dvv,
        duv: d2.duv,
        duuu: s.dn(u, v, 3, 0),
        dvvv: s.dn(u, v, 0, 3),
        duuv: s.dn(u, v, 2, 1),
        duvv: s.dn(u, v, 1, 2),
    }
}

/// The frame's axes as vectors.
fn axes(f: &Frame) -> (Vector, Vector, Vector) {
    (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()))
}

impl RevolutionSurface {
    /// The `nu`th u-derivative of `w` (the basis point from the axis origin, or a basis
    /// derivative) rotated by `u` (given as `u.sin_cos()`).
    fn turned(&self, sin_cos: (f64, f64), w: Vector, nu: u32) -> Vector {
        let z = Vector::from(self.axis().direction);
        // The rotated vector is `w∥ + cos u·w⊥ + sin u·(z × w⊥)`; its u-derivatives drop
        // `w∥` and shift the trigonometric factors.
        let parallel = z * w.dot(z);
        let perpendicular = w - parallel;
        let (cu, su) = cos_sin_derivative_of(sin_cos, nu);
        let turned = perpendicular * cu + z.cross(perpendicular) * su;
        if nu == 0 { parallel + turned } else { turned }
    }
}

impl Surface {
    /// The points at `u` and each of `vs`, as [`ParametricSurface::point`] gives them (with
    /// one sine and cosine of `u` for a cylinder or cone).
    pub(crate) fn points_along_v<const N: usize>(&self, u: f64, vs: [f64; N]) -> [Point; N] {
        match self {
            Surface::Cylinder(c) => {
                let su = u.sin_cos();
                vs.map(|v| c.point_of(su, v))
            }
            Surface::Cone(c) => {
                let su = u.sin_cos();
                vs.map(|v| c.point_of(su, v))
            }
            _ => vs.map(|v| self.point(u, v)),
        }
    }
}

impl Cylinder {
    /// The point at `u` (given as `u.sin_cos()`) and `v`.
    fn point_of(&self, (sin, cos): (f64, f64), v: f64) -> Point {
        let (x, y, z) = axes(self.frame());
        self.frame().origin() + (x * cos + y * sin) * self.radius() + z * v
    }

    /// The derivative `(nu, nv)` (not both 0) at `u` (given as `u.sin_cos()`).
    fn dn_of(&self, sin_cos: (f64, f64), nu: u32, nv: u32) -> Vector {
        let (x, y, z) = axes(self.frame());
        match (nu, nv) {
            (0, 1) => z,
            (0, _) => Vector::ZERO,
            (_, 0) => {
                let (cos, sin) = cos_sin_derivative_of(sin_cos, nu);
                (x * cos + y * sin) * self.radius()
            }
            _ => Vector::ZERO,
        }
    }
}

impl Cone {
    /// The point at `u` (given as `u.sin_cos()`) and `v`.
    fn point_of(&self, (su, cu): (f64, f64), v: f64) -> Point {
        let (x, y, z) = axes(self.frame());
        let (sa, ca) = self.semi_angle_sin_cos();
        self.frame().origin() + (x * cu + y * su) * (self.radius() + v * sa) + z * (v * ca)
    }

    /// The derivative `(nu, nv)` (not both 0) at `u` (given as `u.sin_cos()`) and `v`.
    fn dn_of(&self, sin_cos: (f64, f64), v: f64, nu: u32, nv: u32) -> Vector {
        let (x, y, z) = axes(self.frame());
        let (sa, ca) = self.semi_angle_sin_cos();
        let (cu, su) = cos_sin_derivative_of(sin_cos, nu);
        let radial = x * cu + y * su;
        match (nu, nv) {
            (_, 0) => radial * (self.radius() + v * sa),
            (0, 1) => radial * sa + z * ca,
            (_, 1) => radial * sa,
            _ => Vector::ZERO,
        }
    }
}

impl Sphere {
    /// The point at `u` and `v` (given as their `sin_cos()`).
    fn point_of(&self, (su, cu): (f64, f64), (sv, cv): (f64, f64)) -> Point {
        let (x, y, z) = axes(self.frame());
        self.frame().origin() + ((x * cu + y * su) * cv + z * sv) * self.radius()
    }

    /// The derivative `(nu, nv)` (not both 0) at `u` and `v` (given as their `sin_cos()`).
    fn dn_of(&self, u: (f64, f64), v: (f64, f64), nu: u32, nv: u32) -> Vector {
        let (x, y, z) = axes(self.frame());
        // r·(cos v·cos u, cos v·sin u, sin v): derivatives factor into u and v parts.
        let (cu, su) = cos_sin_derivative_of(u, nu);
        let (cv, sv) = cos_sin_derivative_of(v, nv);
        let height = if nu == 0 { sv } else { 0.0 };
        ((x * cu + y * su) * cv + z * height) * self.radius()
    }
}

impl Torus {
    /// The point at `u` and `v` (given as their `sin_cos()`).
    fn point_of(&self, (su, cu): (f64, f64), (sv, cv): (f64, f64)) -> Point {
        let (x, y, z) = axes(self.frame());
        let (major, minor) = (self.major_radius(), self.minor_radius());
        self.frame().origin() + (x * cu + y * su) * (major + minor * cv) + z * (minor * sv)
    }

    /// The derivative `(nu, nv)` (not both 0) at `u` and `v` (given as their `sin_cos()`).
    fn dn_of(&self, u: (f64, f64), v: (f64, f64), nu: u32, nv: u32) -> Vector {
        let (x, y, z) = axes(self.frame());
        let (cu, su) = cos_sin_derivative_of(u, nu);
        let (cv, sv) = cos_sin_derivative_of(v, nv);
        let (major, minor) = (self.major_radius(), self.minor_radius());
        let radial = if nv == 0 { major + minor * cv } else { minor * cv };
        let height = if nu == 0 { minor * sv } else { 0.0 };
        (x * cu + y * su) * radial + z * height
    }
}

impl ParametricSurface for Plane {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (-INFINITE, INFINITE, -INFINITE, INFINITE)
    }
    fn is_u_closed(&self) -> bool {
        false
    }
    fn is_v_closed(&self) -> bool {
        false
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
        let (x, y, _) = axes(self.frame());
        self.origin() + x * u + y * v
    }
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let (x, y, _) = axes(self.frame());
        (self.point(u, v), x, y)
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        d2_from(self, u, v)
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        d3_from(self, u, v)
    }
    fn dn(&self, _u: f64, _v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        let (x, y, _) = axes(self.frame());
        match (nu, nv) {
            (1, 0) => x,
            (0, 1) => y,
            _ => Vector::ZERO,
        }
    }
    fn u_iso(&self, u: f64) -> Curve {
        Curve::Line(Line::new(self.point(u, 0.0), self.frame().y()))
    }
    fn v_iso(&self, v: f64) -> Curve {
        Curve::Line(Line::new(self.point(0.0, v), self.frame().x()))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        -u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
    fn u_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_x();
        *self = Plane::from_frame(frame);
    }
    fn v_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Plane::from_frame(frame);
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        (tolerance, tolerance)
    }
}

impl ParametricSurface for Cylinder {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, TAU, -INFINITE, INFINITE)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_v_closed(&self) -> bool {
        false
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        false
    }
    fn u_period(&self) -> Option<f64> {
        Some(TAU)
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
        self.point_of(u.sin_cos(), v)
    }
    // The derivatives from the same sine and cosine of u.
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let su = u.sin_cos();
        (self.point_of(su, v), self.dn_of(su, 1, 0), self.dn_of(su, 0, 1))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        let su = u.sin_cos();
        d2_of(self.point_of(su, v), |nu, nv| self.dn_of(su, nu, nv))
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        let su = u.sin_cos();
        d3_of(self.point_of(su, v), |nu, nv| self.dn_of(su, nu, nv))
    }
    fn dn(&self, u: f64, _v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        self.dn_of(u.sin_cos(), nu, nv)
    }
    fn u_iso(&self, u: f64) -> Curve {
        Curve::Line(Line::new(self.point(u, 0.0), self.frame().z()))
    }
    fn v_iso(&self, v: f64) -> Curve {
        let f = self.frame();
        let origin = f.origin() + Vector::from(f.z()) * v;
        Curve::Circle(Circle::new(Frame::from_axes_unchecked(origin, f.x(), f.y(), f.z()), self.radius()))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
    fn u_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Cylinder::new(frame, self.radius());
    }
    fn v_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_z();
        *self = Cylinder::new(frame, self.radius());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let du = if self.radius() > tolerance { tolerance / self.radius() } else { TAU };
        (du, tolerance)
    }
}

/// A cone: `(u, v) -> O + (R + v·sin α)(cos u·X + sin u·Y) + v·cos α·Z`, with `v` the
/// distance along a generating line from the reference circle of radius `R`.
impl ParametricSurface for Cone {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, TAU, -INFINITE, INFINITE)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_v_closed(&self) -> bool {
        false
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        false
    }
    fn u_period(&self) -> Option<f64> {
        Some(TAU)
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
        self.point_of(u.sin_cos(), v)
    }
    // The derivatives from the same sine and cosine of u.
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let su = u.sin_cos();
        (self.point_of(su, v), self.dn_of(su, v, 1, 0), self.dn_of(su, v, 0, 1))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        let su = u.sin_cos();
        d2_of(self.point_of(su, v), |nu, nv| self.dn_of(su, v, nu, nv))
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        let su = u.sin_cos();
        d3_of(self.point_of(su, v), |nu, nv| self.dn_of(su, v, nu, nv))
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        self.dn_of(u.sin_cos(), v, nu, nv)
    }
    /// The generating line at `u`, parametrized by `v`.
    fn u_iso(&self, u: f64) -> Curve {
        let direction = Direction::from_vector(self.dn(u, 0.0, 0, 1));
        Curve::Line(Line::new(self.point(u, 0.0), direction))
    }
    /// The circle at `v` (with its frame flipped where the radius would be negative,
    /// beyond the apex).
    fn v_iso(&self, v: f64) -> Curve {
        let f = self.frame();
        let (sa, ca) = self.semi_angle_sin_cos();
        let radius = self.radius() + v * sa;
        let origin = f.origin() + Vector::from(f.z()) * (v * ca);
        let (x, y) = if radius < 0.0 { (-f.x(), -f.y()) } else { (f.x(), f.y()) };
        Curve::Circle(Circle::new(Frame::from_axes_unchecked(origin, x, y, f.z()), radius.abs()))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
    fn u_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Cone::new(frame, self.semi_angle(), self.radius());
    }
    /// Flips the axis and the semi-angle, so `v` maps to `-v`.
    fn v_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_z();
        *self = Cone::new(frame, -self.semi_angle(), self.radius());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let du = if self.radius() > tolerance { tolerance / self.radius() } else { TAU };
        (du, tolerance)
    }
}

impl ParametricSurface for Sphere {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, TAU, -FRAC_PI_2, FRAC_PI_2)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_v_closed(&self) -> bool {
        false
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        false
    }
    fn u_period(&self) -> Option<f64> {
        Some(TAU)
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
        self.point_of(u.sin_cos(), v.sin_cos())
    }
    // The derivatives from the same sines and cosines.
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        (self.point_of(su, sv), self.dn_of(su, sv, 1, 0), self.dn_of(su, sv, 0, 1))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        d2_of(self.point_of(su, sv), |nu, nv| self.dn_of(su, sv, nu, nv))
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        d3_of(self.point_of(su, sv), |nu, nv| self.dn_of(su, sv, nu, nv))
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        self.dn_of(u.sin_cos(), v.sin_cos(), nu, nv)
    }
    /// The meridian at `u`, parametrized by `v` (a full circle; the sphere uses
    /// `[-π/2, π/2]` of it).
    fn u_iso(&self, u: f64) -> Curve {
        let f = self.frame();
        let (x, y, z) = axes(f);
        let (su, cu) = u.sin_cos();
        let radial = Direction::from_vector(x * cu + y * su);
        let normal = Direction::from_vector(Vector::from(radial).cross(z));
        Curve::Circle(Circle::new(Frame::from_axes_unchecked(f.origin(), radial, f.z(), normal), self.radius()))
    }
    /// The parallel at latitude `v`.
    fn v_iso(&self, v: f64) -> Curve {
        let f = self.frame();
        let (sv, cv) = v.sin_cos();
        let origin = f.origin() + Vector::from(f.z()) * (self.radius() * sv);
        Curve::Circle(Circle::new(Frame::from_axes_unchecked(origin, f.x(), f.y(), f.z()), self.radius() * cv))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
    fn u_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Sphere::new(frame, self.radius());
    }
    fn v_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_z();
        *self = Sphere::new(frame, self.radius());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let r = if self.radius() > tolerance { tolerance / self.radius() } else { PI };
        (r, r)
    }
}

/// A torus: `(u, v) -> O + (R + r·cos v)(cos u·X + sin u·Y) + r·sin v·Z`.
impl ParametricSurface for Torus {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, TAU, 0.0, TAU)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_v_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        true
    }
    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn v_period(&self) -> Option<f64> {
        Some(TAU)
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
        self.point_of(u.sin_cos(), v.sin_cos())
    }
    // The derivatives from the same sines and cosines.
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        (self.point_of(su, sv), self.dn_of(su, sv, 1, 0), self.dn_of(su, sv, 0, 1))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        d2_of(self.point_of(su, sv), |nu, nv| self.dn_of(su, sv, nu, nv))
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        let (su, sv) = (u.sin_cos(), v.sin_cos());
        d3_of(self.point_of(su, sv), |nu, nv| self.dn_of(su, sv, nu, nv))
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        self.dn_of(u.sin_cos(), v.sin_cos(), nu, nv)
    }
    /// The meridian circle at `u`, parametrized by `v`.
    fn u_iso(&self, u: f64) -> Curve {
        let f = self.frame();
        let (x, y, z) = axes(f);
        let (su, cu) = u.sin_cos();
        let radial = x * cu + y * su;
        let normal = Direction::from_vector(radial.cross(z));
        let center = f.origin() + radial * self.major_radius();
        let frame = Frame::from_axes_unchecked(center, Direction::from_vector(radial), f.z(), normal);
        Curve::Circle(Circle::new(frame, self.minor_radius()))
    }
    /// The parallel circle at `v`.
    fn v_iso(&self, v: f64) -> Curve {
        let f = self.frame();
        let (sv, cv) = v.sin_cos();
        let radius = self.major_radius() + self.minor_radius() * cv;
        let origin = f.origin() + Vector::from(f.z()) * (self.minor_radius() * sv);
        // (On a spindle torus's inner part, where the radius would be negative, the circle of
        // its size with x and y reversed (OCCT's `ElSLib::TorusVIso`): the same points.)
        if radius < 0.0 {
            return Curve::Circle(Circle::new(Frame::from_axes_unchecked(origin, -f.x(), -f.y(), f.z()), -radius));
        }
        Curve::Circle(Circle::new(Frame::from_axes_unchecked(origin, f.x(), f.y(), f.z()), radius))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        TAU - v
    }
    fn u_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_y();
        *self = Torus::new(frame, self.major_radius(), self.minor_radius());
    }
    fn v_reverse(&mut self) {
        let mut frame = *self.frame();
        frame.flip_z();
        *self = Torus::new(frame, self.major_radius(), self.minor_radius());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let outer = self.major_radius() + self.minor_radius();
        let du = if outer > tolerance { tolerance / outer } else { TAU };
        let dv = if self.minor_radius() > tolerance { tolerance / self.minor_radius() } else { TAU };
        (du, dv)
    }
}

/// `(u, v) -> basis(v)` rotated by `u` about the axis.
impl ParametricSurface for RevolutionSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, TAU, self.basis().first_parameter(), self.basis().last_parameter())
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_v_closed(&self) -> bool {
        self.basis().is_closed()
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        self.basis().is_periodic()
    }
    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }
    fn v_period(&self) -> Option<f64> {
        self.basis().period()
    }
    fn continuity(&self) -> Continuity {
        self.basis().continuity()
    }
    fn is_cn_u(&self, _n: u32) -> bool {
        true
    }
    fn is_cn_v(&self, n: u32) -> bool {
        self.basis().is_cn(n)
    }
    fn point(&self, u: f64, v: f64) -> Point {
        self.basis().point(v).transformed(&Transform::rotation(self.axis(), u))
    }
    // The derivatives from the same sine and cosine of u and basis point.
    // (The point as [`point`](Self::point) gives it, from the one basis point.)
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let su = u.sin_cos();
        let basis = self.basis().point(v);
        let at = basis - self.axis().origin;
        (basis.transformed(&Transform::rotation(self.axis(), u)), self.turned(su, at, 1), self.turned(su, self.basis().dn(v, 1), 0))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        let su = u.sin_cos();
        let basis = self.basis().point(v);
        let at = basis - self.axis().origin;
        let (d1, d2) = (self.basis().dn(v, 1), self.basis().dn(v, 2));
        SurfaceD2 {
            point: basis.transformed(&Transform::rotation(self.axis(), u)),
            du: self.turned(su, at, 1),
            dv: self.turned(su, d1, 0),
            duu: self.turned(su, at, 2),
            dvv: self.turned(su, d2, 0),
            duv: self.turned(su, d1, 1),
        }
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        d3_from(self, u, v)
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        let w = if nv == 0 { self.basis().point(v) - self.axis().origin } else { self.basis().dn(v, nv) };
        self.turned(u.sin_cos(), w, nu)
    }
    /// The basis rotated by `u`.
    fn u_iso(&self, u: f64) -> Curve {
        self.basis().transformed(&Transform::rotation(self.axis(), u))
    }
    /// The circle described by the basis point at `v`.
    fn v_iso(&self, v: f64) -> Curve {
        let axis = self.axis();
        let z = Vector::from(axis.direction);
        let w = self.basis().point(v) - axis.origin;
        let center = axis.origin + z * w.dot(z);
        let radial = w - z * w.dot(z);
        let frame = match Direction::try_from_vector(radial) {
            Some(x) => Frame::new(center, axis.direction, x),
            None => Frame::from_z(center, axis.direction),
        };
        Curve::Circle(Circle::new(frame, radial.length()))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.basis().reversed_parameter(v)
    }
    /// Reverses the axis, so `u` maps to `2π - u`.
    fn u_reverse(&mut self) {
        let axis = Axis::new(self.axis().origin, -self.axis().direction);
        *self = RevolutionSurface::new(self.basis().clone(), axis);
    }
    fn v_reverse(&mut self) {
        let mut basis = self.basis().clone();
        basis.reverse();
        *self = RevolutionSurface::new(basis, self.axis());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let mid = 0.5 * (self.basis().first_parameter() + self.basis().last_parameter());
        let w = self.basis().point(mid) - self.axis().origin;
        let z = Vector::from(self.axis().direction);
        let radius = (w - z * w.dot(z)).length();
        let du = if radius > tolerance { tolerance / radius } else { TAU };
        (du, self.basis().resolution(tolerance))
    }
}

impl ParametricSurface for ExtrusionSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.basis().first_parameter(), self.basis().last_parameter(), -INFINITE, INFINITE)
    }
    fn is_u_closed(&self) -> bool {
        self.basis().is_closed()
    }
    fn is_v_closed(&self) -> bool {
        false
    }
    fn is_u_periodic(&self) -> bool {
        self.basis().is_periodic()
    }
    fn is_v_periodic(&self) -> bool {
        false
    }
    fn u_period(&self) -> Option<f64> {
        self.basis().period()
    }
    fn v_period(&self) -> Option<f64> {
        None
    }
    fn continuity(&self) -> Continuity {
        self.basis().continuity()
    }
    fn is_cn_u(&self, n: u32) -> bool {
        self.basis().is_cn(n)
    }
    fn is_cn_v(&self, _n: u32) -> bool {
        true
    }
    fn point(&self, u: f64, v: f64) -> Point {
        self.basis().point(u) + Vector::from(self.direction()) * v
    }
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        let (p, du) = self.basis().d1(u);
        (p + Vector::from(self.direction()) * v, du, Vector::from(self.direction()))
    }
    fn d2(&self, u: f64, v: f64) -> SurfaceD2 {
        d2_from(self, u, v)
    }
    fn d3(&self, u: f64, v: f64) -> SurfaceD3 {
        d3_from(self, u, v)
    }
    fn dn(&self, u: f64, _v: f64, nu: u32, nv: u32) -> Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        match (nu, nv) {
            (0, 1) => Vector::from(self.direction()),
            (_, 0) => self.basis().dn(u, nu),
            _ => Vector::new(0.0, 0.0, 0.0),
        }
    }
    /// The straight line through the basis point at `u`.
    fn u_iso(&self, u: f64) -> Curve {
        Curve::Line(Line::new(self.basis().point(u), self.direction()))
    }
    /// The basis moved by `v` along the direction.
    fn v_iso(&self, v: f64) -> Curve {
        self.basis().transformed(&Transform::translation(Vector::from(self.direction()) * v))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.basis().reversed_parameter(u)
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
    fn u_reverse(&mut self) {
        let mut basis = self.basis().clone();
        basis.reverse();
        *self = ExtrusionSurface::new(basis, self.direction());
    }
    fn v_reverse(&mut self) {
        *self = ExtrusionSurface::new(self.basis().clone(), -self.direction());
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        (self.basis().resolution(tolerance), tolerance)
    }
}

// --- parabolas and hyperbolas ------------------------------------------------------------------

/// `ParametricCurve` for an unbounded conic `t -> origin + fx(t)·x + fy(t)·y`, given the
/// coefficient functions' n-th derivatives and how to reverse it (flipping the frame's Y
/// axis maps `t` to `-t`).
macro_rules! open_conic {
    ($ty:ty, $point:ty, $vector:ty, $coefficients:expr, $flip:expr, $speed:expr) => {
        impl ParametricCurve for $ty {
            type Point = $point;
            type Vector = $vector;
            fn first_parameter(&self) -> f64 {
                -INFINITE
            }
            fn last_parameter(&self) -> f64 {
                INFINITE
            }
            fn is_closed(&self) -> bool {
                false
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
            fn point(&self, t: f64) -> $point {
                let (a, b) = ($coefficients)(self, t, 0);
                let f = self.frame();
                f.origin() + <$vector>::from(f.x()) * a + <$vector>::from(f.y()) * b
            }
            fn d1(&self, t: f64) -> ($point, $vector) {
                (self.point(t), self.dn(t, 1))
            }
            fn d2(&self, t: f64) -> ($point, $vector, $vector) {
                (self.point(t), self.dn(t, 1), self.dn(t, 2))
            }
            fn d3(&self, t: f64) -> ($point, $vector, $vector, $vector) {
                (self.point(t), self.dn(t, 1), self.dn(t, 2), self.dn(t, 3))
            }
            fn dn(&self, t: f64, n: u32) -> $vector {
                assert!(n > 0, "derivative of order 0");
                let (a, b) = ($coefficients)(self, t, n);
                let f = self.frame();
                <$vector>::from(f.x()) * a + <$vector>::from(f.y()) * b
            }
            fn reversed_parameter(&self, t: f64) -> f64 {
                -t
            }
            fn reverse(&mut self) {
                ($flip)(self)
            }
            fn resolution(&self, tolerance: f64) -> f64 {
                tolerance / ($speed)(self)
            }
        }
    };
}

/// The n-th derivatives of a parabola's coefficients `(t²/(4f), t)` (`(t, 0)` for `f = 0`).
fn parabola_coefficients(focal: f64, t: f64, n: u32) -> (f64, f64) {
    if focal == 0.0 {
        return match n {
            0 => (t, 0.0),
            1 => (1.0, 0.0),
            _ => (0.0, 0.0),
        };
    }
    match n {
        0 => (t * t / (4.0 * focal), t),
        1 => (t / (2.0 * focal), 1.0),
        2 => (1.0 / (2.0 * focal), 0.0),
        _ => (0.0, 0.0),
    }
}

/// The n-th derivatives of a hyperbola's coefficients `(a cosh t, b sinh t)`.
fn hyperbola_coefficients(major: f64, minor: f64, t: f64, n: u32) -> (f64, f64) {
    let (c, s) = (t.cosh(), t.sinh());
    if n % 2 == 0 { (major * c, minor * s) } else { (major * s, minor * c) }
}

open_conic!(
    Parabola,
    Point,
    Vector,
    |c: &Parabola, t, n| parabola_coefficients(c.focal(), t, n),
    |c: &mut Parabola| {
        let mut frame = *c.frame();
        frame.flip_y();
        *c = Parabola::new(frame, c.focal());
    },
    |c: &Parabola| if c.focal() > 0.0 { 1.0f64.max(1e3 / (2.0 * c.focal())) } else { 1.0 }
);
open_conic!(
    Parabola2,
    Point2,
    Vector2,
    |c: &Parabola2, t, n| parabola_coefficients(c.focal(), t, n),
    |c: &mut Parabola2| {
        let f = c.frame();
        let frame = crate::math::Frame2::from_axes(f.origin(), f.x(), -f.y()).expect("orthonormal frame");
        *c = Parabola2::new(frame, c.focal());
    },
    |c: &Parabola2| if c.focal() > 0.0 { 1.0f64.max(1e3 / (2.0 * c.focal())) } else { 1.0 }
);
open_conic!(
    Hyperbola,
    Point,
    Vector,
    |c: &Hyperbola, t, n| hyperbola_coefficients(c.major_radius(), c.minor_radius(), t, n),
    |c: &mut Hyperbola| {
        let mut frame = *c.frame();
        frame.flip_y();
        *c = Hyperbola::new(frame, c.major_radius(), c.minor_radius());
    },
    |c: &Hyperbola| c.major_radius().max(c.minor_radius()).max(1e-300)
);
open_conic!(
    Hyperbola2,
    Point2,
    Vector2,
    |c: &Hyperbola2, t, n| hyperbola_coefficients(c.major_radius(), c.minor_radius(), t, n),
    |c: &mut Hyperbola2| {
        let f = c.frame();
        let frame = crate::math::Frame2::from_axes(f.origin(), f.x(), -f.y()).expect("orthonormal frame");
        *c = Hyperbola2::new(frame, c.major_radius(), c.minor_radius());
    },
    |c: &Hyperbola2| c.major_radius().max(c.minor_radius()).max(1e-300)
);

// --- offset curves -----------------------------------------------------------------------------

/// The angle within which a C0 B-spline's pieces meet tangentially, for an offset curve's
/// continuity (OCCT's `Precision::Angular`).
const G1_ANGLE: f64 = 1e-12;

/// The arithmetic of OCCT's `gp_XYZ` / `gp_XY` that the offset-curve formulas use, for
/// [`Vector`] and [`Vector2`] alike.
trait OffsetAlgebra:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<f64, Output = Self> + Div<f64, Output = Self> + Neg<Output = Self>
{
    fn dot_with(self, other: Self) -> f64;
}

impl OffsetAlgebra for Vector {
    fn dot_with(self, other: Self) -> f64 {
        self.dot(other)
    }
}

impl OffsetAlgebra for Vector2 {
    fn dot_with(self, other: Self) -> f64 {
        self.dot(other)
    }
}

/// OCCT's `Geom_OffsetCurveUtils::CalculateD0`: the point `p + offset·N/|N|` of a 3D offset
/// curve, `N = normal(d1)` (`d1 ^ direction`); `None` where `|N|` is not above
/// `gp::Resolution`.
fn offset_curve_d0(p: Point, d1: Vector, normal: impl Fn(Vector) -> Vector, offset: f64) -> Option<Point> {
    let n = normal(d1);
    let r = n.dot_with(n).sqrt();
    if r <= RESOLUTION {
        return None;
    }
    Some(p + n * (offset / r))
}

/// OCCT's `Geom2d_OffsetCurveUtils::CalculateD0`: the point `p + offset·N/|N|` of a 2D offset
/// curve, `N = normal(d1)` (`(d1.y, -d1.x)`, normalized as `gp_Dir2d` does); `None` where
/// `|d1|²` is not above `gp::Resolution`.
fn offset_curve2_d0(p: Point2, d1: Vector2, normal: impl Fn(Vector2) -> Vector2, offset: f64) -> Option<Point2> {
    if d1.dot_with(d1) <= RESOLUTION {
        return None;
    }
    let n = normal(d1);
    let square = n.dot_with(n);
    // gp_Dir2d's constructor: an already unit vector kept as it is.
    let n = if (1.0 - RESOLUTION..=1.0 + RESOLUTION).contains(&square) { n } else { n / square.sqrt() };
    Some(p + n * offset)
}

/// OCCT's `Geom_OffsetCurveUtils::CalculateD1` / `Geom2d_OffsetCurveUtils::CalculateD1` (the
/// same steps): the point and first derivative of an offset curve from the basis point `p`
/// and derivatives `d1`, `d2`; `None` where `|N|²` is not above `gp::Resolution`.
fn offset_curve_d1<P: Add<V, Output = P>, V: OffsetAlgebra>(
    p: P,
    d1: V,
    d2: V,
    normal: impl Fn(V) -> V,
    offset: f64,
) -> Option<(P, V)> {
    // P(u) = p(u) + Offset * Ndir / R
    // with R = || p' ^ V|| and Ndir = P' ^ direction (local normal direction)
    // P'(u) = p'(u) + (Offset / R**2) * (DNdir/DU * R -  Ndir * (DR/R))
    let n = normal(d1);
    let mut dn = normal(d2);
    let r2 = n.dot_with(n);
    let r = r2.sqrt();
    let r3 = r * r2;
    let dr = n.dot_with(dn);
    if r3 <= RESOLUTION {
        if r2 <= RESOLUTION {
            return None;
        }
        // OCCT: "We try another computation but the stability is not very good."
        dn = dn * r;
        dn = dn - n * (dr / r);
        dn = dn * (offset / r2);
    } else {
        // OCCT: "Same computation as IICURV in EUCLID-IS because the stability is better".
        dn = dn * (offset / r);
        dn = dn - n * (offset * dr / r3);
    }
    Some((p + n * (offset / r), d1 + dn))
}

/// OCCT's `Geom_OffsetCurveUtils::CalculateD2` / `Geom2d_OffsetCurveUtils::CalculateD2` (the
/// same steps): the point and derivatives to the second of an offset curve from the basis
/// point and derivatives `d1`..`d3`, the basis second derivative turned over where
/// `direction_change` (set by [`adjust_offset_derivatives`]); `None` where `|N|⁴` is not
/// above `gp::Resolution`.
fn offset_curve_d2<P: Add<V, Output = P>, V: OffsetAlgebra>(
    p: P,
    d1: V,
    d2: V,
    d3: V,
    normal: impl Fn(V) -> V,
    offset: f64,
    direction_change: bool,
) -> Option<(P, V, V)> {
    // P"(u) = p"(u) + (Offset / R) * (D2Ndir/DU - DNdir * (2.0 * Dr/ R**2) +
    //         Ndir * ( (3.0 * Dr**2 / R**4) - (D2r / R**2)))
    let n = normal(d1);
    let mut dn = normal(d2);
    let mut d2n = normal(d3);
    let r2 = n.dot_with(n);
    let r = r2.sqrt();
    let r3 = r2 * r;
    let r4 = r2 * r2;
    let r5 = r3 * r2;
    let dr = n.dot_with(dn);
    let d2r = n.dot_with(d2n) + dn.dot_with(dn);
    if r5 <= RESOLUTION {
        if r4 <= RESOLUTION {
            return None;
        }
        // OCCT: "We try another computation but the stability is not very good dixit ISG."
        d2n = d2n - dn * (2.0 * dr / r2);
        d2n = d2n + n * (((3.0 * dr * dr) / r4) - (d2r / r2));
        d2n = d2n * (offset / r);
        dn = dn * r;
        dn = dn - n * (dr / r);
        dn = dn * (offset / r2);
    } else {
        // OCCT: "Same computation as IICURV in EUCLID-IS because the stability is better."
        d2n = d2n * (offset / r);
        d2n = d2n - dn * (2.0 * offset * dr / r3);
        d2n = d2n + n * (offset * (((3.0 * dr * dr) / r5) - (d2r / r3)));
        dn = dn * (offset / r);
        dn = dn - n * (offset * dr / r3);
    }
    let d2 = if direction_change { -d2 } else { d2 };
    Some((p + n * (offset / r), d1 + dn, d2 + d2n))
}

/// The low-`|N|` branch shared by OCCT's `Geom_OffsetCurveUtils::CalculateD3` and
/// `Geom2d_OffsetCurveUtils::CalculateD3`: the offset terms `(DNdir, D2Ndir, D3Ndir)`.
/// Kept as OCCT has it: its second-derivative term subtracts the `Ndir` term that
/// `CalculateD2` adds, and its third-derivative coefficients are not those of the
/// formula's derivative.
#[allow(clippy::too_many_arguments)]
fn offset_d3_terms_low<V: OffsetAlgebra>(
    n: V,
    mut dn: V,
    mut d2n: V,
    mut d3n: V,
    [r, r2, r4, r6]: [f64; 4],
    [dr, d2r, d3r]: [f64; 3],
    offset: f64,
) -> (V, V, V) {
    // V3 = P"' (U) :
    d3n = d3n - d2n * (3.0 * dr / r2);
    d3n = d3n - dn * (3.0 * ((d2r / r2) + (dr * dr / r4)));
    d3n = d3n + n * (6.0 * dr * dr / r4 + 6.0 * dr * d2r / r4 - 15.0 * dr * dr * dr / r6 - d3r);
    d3n = d3n * (offset / r);
    // V2 = P" (U) :
    d2n = d2n - dn * (2.0 * dr / r2);
    d2n = d2n - n * ((3.0 * dr * dr / r4) - (d2r / r2));
    d2n = d2n * (offset / r);
    // V1 = P' (U) :
    dn = dn * r;
    dn = dn - n * (dr / r);
    dn = dn * (offset / r2);
    (dn, d2n, d3n)
}

/// OCCT's `Geom_OffsetCurveUtils::CalculateD3` (3D, `three_d`) or
/// `Geom2d_OffsetCurveUtils::CalculateD3` (2D, whose regular branch applies the offset in
/// another order): the point and derivatives to the third of an offset curve from the basis
/// point and derivatives `d1`..`d4`, the basis third derivative turned over where
/// `direction_change`; `None` where `|N|⁶` is not above `gp::Resolution`. Kept as OCCT has
/// it: the second derivative subtracts the `Ndir` term that `CalculateD2` adds (so it differs
/// from [`offset_curve_d2`]'s), and the third derivative's `DNdir` and `Ndir` coefficients
/// are not those of the formula's derivative.
#[allow(clippy::too_many_arguments)]
fn offset_curve_d3<P: Add<V, Output = P>, V: OffsetAlgebra>(
    p: P,
    d1: V,
    d2: V,
    d3: V,
    d4: V,
    normal: impl Fn(V) -> V,
    offset: f64,
    direction_change: bool,
    three_d: bool,
) -> Option<(P, V, V, V)> {
    // P"'(u) = p"'(u) + (Offset / R) * (D3Ndir - (3.0 * Dr/R**2) * D2Ndir -
    //          (3.0 * D2r / R2) * DNdir + (3.0 * Dr * Dr / R4) * DNdir -
    //          (D3r/R2) * Ndir + (6.0 * Dr * Dr / R4) * Ndir +
    //          (6.0 * Dr * D2r / R4) * Ndir - (15.0 * Dr* Dr* Dr /R6) * Ndir
    let n = normal(d1);
    let dn = normal(d2);
    let d2n = normal(d3);
    let d3n = normal(d4);
    let r2 = n.dot_with(n);
    let r = r2.sqrt();
    let r3 = r2 * r;
    let r4 = r2 * r2;
    let r5 = r3 * r2;
    let r6 = r3 * r3;
    let r7 = r5 * r2;
    let dr = n.dot_with(dn);
    let d2r = n.dot_with(d2n) + dn.dot_with(dn);
    let d3r = n.dot_with(d3n) + 3.0 * dn.dot_with(d2n);
    let (dn, d2n, d3n) = if r7 <= RESOLUTION {
        if r6 <= RESOLUTION {
            return None;
        }
        offset_d3_terms_low(n, dn, d2n, d3n, [r, r2, r4, r6], [dr, d2r, d3r], offset)
    } else if three_d {
        // V3 = P"' (U) :
        let mut d3n = d3n / r;
        d3n = d3n - d2n * (3.0 * dr / r3);
        d3n = d3n - dn * (3.0 * ((d2r / r3) + (dr * dr) / r5));
        d3n = d3n + n * (6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5 - 15.0 * dr * dr * dr / r7 - d3r);
        d3n = d3n * offset;
        // V2 = P" (U) :
        let mut d2n = d2n / r;
        d2n = d2n - dn * (2.0 * dr / r3);
        d2n = d2n - n * ((3.0 * dr * dr / r5) - (d2r / r3));
        d2n = d2n * offset;
        // V1 = P' (U) :
        let dn = dn * (offset / r) - n * (offset * dr / r3);
        (dn, d2n, d3n)
    } else {
        // V3 = P"' (U) :
        let mut d3n = d3n * (offset / r);
        d3n = d3n - d2n * (3.0 * offset * dr / r3);
        d3n = d3n - dn * ((3.0 * offset) * ((d2r / r3) + (dr * dr) / r5));
        d3n = d3n + n * (offset * (6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5 - 15.0 * dr * dr * dr / r7 - d3r));
        // V2 = P" (U) :
        let mut d2n = d2n * (offset / r);
        d2n = d2n - dn * (2.0 * offset * dr / r3);
        d2n = d2n - n * (offset * (((3.0 * dr * dr) / r5) - (d2r / r3)));
        // V1 = P' (U) :
        let dn = dn * (offset / r) - n * (offset * dr / r3);
        (dn, d2n, d3n)
    };
    let d3 = if direction_change { -d3 } else { d3 };
    Some((p + n * (offset / r), d1 + dn, d2 + d2n, d3 + d3n))
}

/// OCCT's `Geom_OffsetCurveUtils::AdjustDerivative` / `Geom2d_OffsetCurveUtils::
/// AdjustDerivative`: at a basis point whose first derivative vanishes, the first non-null
/// derivative of order 2 or 3 (3 at most) stands for it and the following ones for the next
/// `max_derivative - 1`, all turned over when the curve runs against it (a chord to a point
/// `max(1e-3·range, 1e-7)` away tells); the turn returned.
fn adjust_offset_derivatives<C, V>(curve: &C, max_derivative: u32, u: f64, derivatives: [&mut V; 4]) -> bool
where
    C: ParametricCurve<Vector = V>,
    C::Point: Sub<Output = V>,
    V: OffsetAlgebra,
{
    const MIN_STEP: f64 = 1e-7;
    const MAX_DERIVATIVE_ORDER: u32 = 3;
    const DIVISION_FACTOR: f64 = 1e-3;
    let infimum = curve.first_parameter();
    let supremum = curve.last_parameter();
    // OCCT tests against RealLast/RealFirst, which kentron's infinite bounds (as OCCT's,
    // Precision::Infinite) never reach.
    let du = if supremum >= f64::MAX || infimum <= f64::MIN { 0.0 } else { supremum - infimum };
    let delta = (du * DIVISION_FACTOR).max(MIN_STEP);

    // Derivative is approximated by Taylor-series.
    let mut index = 1;
    let mut v;
    loop {
        index += 1;
        v = curve.dn(u, index);
        if !(v.dot_with(v) <= RESOLUTION && index < MAX_DERIVATIVE_ORDER) {
            break;
        }
    }
    let other = if u - infimum < delta { u + delta } else { u - delta };
    let p1 = curve.point(if other < u { other } else { u });
    let p2 = curve.point(if u < other { other } else { u });
    let direction_change = v.dot_with(p2 - p1) < 0.0;
    let sign = if direction_change { -1.0 } else { 1.0 };

    let [d1, d2, d3, d4] = derivatives;
    *d1 = v * sign;
    for (i, d) in [d2, d3, d4].into_iter().enumerate().take(max_derivative as usize - 1) {
        *d = curve.dn(u, index + i as u32 + 1) * sign;
    }
    direction_change
}

/// `ParametricCurve` for an offset curve `B(t) + offset·N/|N|`, with `N` a linear function
/// of the basis tangent (`normal`), evaluated as OCCT's `Geom_OffsetCurve` /
/// `Geom2d_OffsetCurve` (`EvalD0`..`EvalD3`, `EvalDN`) do: `d0` is the dimension's
/// `CalculateD0`, `three_d` picks the dimension's `CalculateD3`. Where OCCT throws
/// (`Geom_UndefinedValue` / `Geom_UndefinedDerivative`: a normal too short), the
/// evaluation is kentron's earlier one (the direct formula, central differences of the first
/// derivative above it). `rebuild` makes the reversed curve from the reversed basis.
macro_rules! offset_curve {
    ($ty:ty, $point:ty, $vector:ty, $normal:expr, $d0:ident, $three_d:expr, $rebuild:expr, $g1:expr) => {
        impl $ty {
            /// Kentron's earlier derivatives of order `n` (central differences of the first),
            /// for where OCCT's evaluation throws.
            fn central_difference(&self, t: f64, n: u32) -> $vector {
                if n == 1 {
                    return self.d1(t).1;
                }
                let h = 1e-4 * (1.0 + t.abs());
                (self.central_difference(t + h, n - 1) - self.central_difference(t - h, n - 1)) * (0.5 / h)
            }
        }

        impl ParametricCurve for $ty {
            type Point = $point;
            type Vector = $vector;
            fn first_parameter(&self) -> f64 {
                self.basis().first_parameter()
            }
            fn last_parameter(&self) -> f64 {
                self.basis().last_parameter()
            }
            fn is_closed(&self) -> bool {
                self.basis().is_closed()
            }
            fn is_periodic(&self) -> bool {
                self.basis().is_periodic()
            }
            fn period(&self) -> Option<f64> {
                self.basis().period()
            }
            fn continuity(&self) -> Continuity {
                // As OCCT's Geom_OffsetCurve: one order less than the basis, its tangent
                // continuity kept; a C0 B-spline basis tangent-continuous counted as G1.
                match self.basis().continuity() {
                    Continuity::CN => Continuity::CN,
                    Continuity::C3 => Continuity::C2,
                    Continuity::C2 => Continuity::C1,
                    Continuity::G2 => Continuity::G2,
                    Continuity::G1 => Continuity::G1,
                    Continuity::C0 if ($g1)(self.basis()) => Continuity::G1,
                    _ => Continuity::C0,
                }
            }
            fn is_cn(&self, n: u32) -> bool {
                self.basis().is_cn(n + 1)
            }
            /// OCCT's `EvalD0`.
            fn point(&self, t: f64) -> $point {
                let (p, d) = self.basis().d1(t);
                let normal = |v: $vector| ($normal)(self, v);
                $d0(p, d, normal, self.offset()).unwrap_or_else(|| {
                    // OCCT throws Geom_UndefinedValue.
                    let n = normal(d);
                    p + n * (self.offset() / n.length())
                })
            }
            /// OCCT's `EvalD1`.
            fn d1(&self, t: f64) -> ($point, $vector) {
                let (p, d, dd) = self.basis().d2(t);
                let normal = |v: $vector| ($normal)(self, v);
                offset_curve_d1(p, d, dd, normal, self.offset()).unwrap_or_else(|| {
                    // OCCT throws Geom_UndefinedDerivative. With n = N/|N|:
                    // n' = (N' - n (n·N'))/|N|, N' = normal(B'').
                    let big = normal(d);
                    let len = big.length();
                    let n = big * (1.0 / len);
                    let big_prime = normal(dd);
                    let n_prime = (big_prime - n * n.dot(big_prime)) * (1.0 / len);
                    (p + n * self.offset(), d + n_prime * self.offset())
                })
            }
            /// OCCT's `EvalD2`.
            fn d2(&self, t: f64) -> ($point, $vector, $vector) {
                let (p, mut d1, mut d2, mut d3) = self.basis().d3(t);
                let mut direction_change = false;
                if d1.length_squared() <= RESOLUTION {
                    let mut unused_d4 = <$vector>::default();
                    direction_change =
                        adjust_offset_derivatives(self.basis(), 3, t, [&mut d1, &mut d2, &mut d3, &mut unused_d4]);
                }
                let normal = |v: $vector| ($normal)(self, v);
                offset_curve_d2(p, d1, d2, d3, normal, self.offset(), direction_change).unwrap_or_else(|| {
                    // OCCT throws Geom_UndefinedDerivative.
                    (self.point(t), self.d1(t).1, self.central_difference(t, 2))
                })
            }
            /// OCCT's `EvalD3`.
            fn d3(&self, t: f64) -> ($point, $vector, $vector, $vector) {
                let (p, mut d1, mut d2, mut d3) = self.basis().d3(t);
                let mut d4 = self.basis().dn(t, 4);
                let mut direction_change = false;
                if d1.length_squared() <= RESOLUTION {
                    direction_change =
                        adjust_offset_derivatives(self.basis(), 4, t, [&mut d1, &mut d2, &mut d3, &mut d4]);
                }
                let normal = |v: $vector| ($normal)(self, v);
                offset_curve_d3(p, d1, d2, d3, d4, normal, self.offset(), direction_change, $three_d).unwrap_or_else(
                    || {
                        // OCCT throws Geom_UndefinedDerivative.
                        (self.point(t), self.d1(t).1, self.central_difference(t, 2), self.central_difference(t, 3))
                    },
                )
            }
            /// OCCT's `EvalDN`: from `d1`..`d3` up to the third order, the basis curve's
            /// derivative above it (OCCT's choice, kept).
            fn dn(&self, t: f64, n: u32) -> $vector {
                assert!(n > 0, "derivative of order 0");
                match n {
                    1 => self.d1(t).1,
                    2 => self.d2(t).2,
                    3 => self.d3(t).3,
                    _ => self.basis().dn(t, n),
                }
            }
            fn reversed_parameter(&self, t: f64) -> f64 {
                self.basis().reversed_parameter(t)
            }
            fn reverse(&mut self) {
                // The reversed basis turns the offset side over: the offset changes sign.
                let mut basis = self.basis().clone();
                basis.reverse();
                *self = ($rebuild)(self, basis);
            }
            fn resolution(&self, tolerance: f64) -> f64 {
                self.basis().resolution(tolerance)
            }
            fn breaks(&self) -> Vec<f64> {
                self.basis().breaks()
            }
        }
    };
}

offset_curve!(
    OffsetCurve,
    Point,
    Vector,
    |c: &OffsetCurve, d: Vector| d.cross(Vector::from(c.direction())),
    offset_curve_d0,
    true,
    |c: &OffsetCurve, basis: Curve| OffsetCurve::new(basis, -c.offset(), c.direction()),
    |b: &Curve| match b {
        Curve::BSpline(b) => b.is_g1(b.first_parameter(), b.last_parameter(), G1_ANGLE),
        Curve::Trimmed(t) => matches!(t.basis(), Curve::BSpline(b) if b.is_g1(t.first_parameter(), t.last_parameter(), G1_ANGLE)),
        _ => false,
    }
);
/// `Geom_OffsetCurveUtils::EvaluateDN(U, …, 3, DN)`, as `GeomAdaptor_Curve::DN(U, 3)`
/// evaluates an offset curve: `EvaluateD3` given one vector for both its first and second
/// derivatives, so that the basis second derivative (written last) stands for the first
/// too, in the test for a vanishing tangent, after `AdjustDerivative` and in `CalculateD3`
/// (which reads them all before it writes). `None` where OCCT throws.
pub(crate) fn offset_adaptor_third_derivative(c: &OffsetCurve, t: f64) -> Option<Vector> {
    let (p, _, d2, d3) = c.basis().d3(t);
    let (mut d1, mut d2, mut d3) = (d2, d2, d3);
    let mut d4 = c.basis().dn(t, 4);
    let mut direction_change = false;
    if d1.length_squared() <= RESOLUTION {
        direction_change = adjust_offset_derivatives(c.basis(), 4, t, [&mut d1, &mut d2, &mut d3, &mut d4]);
        d1 = d2;
    }
    let normal = |v: Vector| v.cross(Vector::from(c.direction()));
    offset_curve_d3(p, d1, d2, d3, d4, normal, c.offset(), direction_change, true).map(|r| r.3)
}

offset_curve!(
    OffsetCurve2,
    Point2,
    Vector2,
    |_c: &OffsetCurve2, d: Vector2| Vector2::new(d.y, -d.x),
    offset_curve2_d0,
    false,
    |c: &OffsetCurve2, basis: Curve2| OffsetCurve2::new(basis, -c.offset()),
    |b: &Curve2| match b {
        Curve2::BSpline(b) => b.is_g1(b.first_parameter(), b.last_parameter(), G1_ANGLE),
        Curve2::Trimmed(t) => matches!(t.basis(), Curve2::BSpline(b) if b.is_g1(t.first_parameter(), t.last_parameter(), G1_ANGLE)),
        _ => false,
    }
);

#[cfg(test)]
mod tests {
    use super::*;

    /// Points stepped along a circle by rotations agree with those worked out one by one.
    #[test]
    fn stepped_conic_samples_agree() {
        let frame = Frame::from_z(Point::new(3.0, -2.0, 1.0), Direction::from_vector(Vector::new(0.3, 0.4, 1.0)));
        let circle = Curve::Circle(Circle::new(frame, 7.5));
        for (first, last, n) in [(0.0, TAU, 180), (-40.0, -39.0, 64), (2.0, 1.5, 9)] {
            let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
            for (t, p) in params.iter().zip(circle.points_along(&params)) {
                assert!(p.distance(circle.point(*t)) <= 1e-13, "at {t}");
            }
        }
    }

    fn xyz(v: Vector) -> [f64; 3] {
        [v.x, v.y, v.z]
    }

    fn xy(v: Vector2) -> [f64; 2] {
        [v.x, v.y]
    }

    /// A degree-5 B-spline over three spans (not planar in 3D), and a cubic whose first and
    /// last two poles coincide (null first derivative at both ends).
    fn offset_test_poles() -> (Vec<Point>, Vec<Point>) {
        let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let regular = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 2.0, 0.5),
            p(2.0, -1.0, 1.5),
            p(3.0, 1.0, -0.5),
            p(4.0, 3.0, 1.0),
            p(5.0, 0.0, 2.0),
            p(6.0, 2.0, 0.5),
            p(7.0, 1.0, 1.0),
        ];
        let singular = vec![p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.5), p(3.0, 1.0, 2.0), p(3.0, 1.0, 2.0)];
        (regular, singular)
    }

    /// A parameter and seven OCCT vectors at it.
    type OffsetReference<const N: usize> = (f64, [f64; N], [f64; N], [f64; N], [f64; N], [f64; N], [f64; N], [f64; N]);

    #[track_caller]
    fn assert_close<const N: usize>(got: [f64; N], want: [f64; N], what: &str) {
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() <= 1e-12 * w.abs().max(1.0), "{what}: {got:?} against OCCT's {want:?}");
        }
    }

    /// Offset curves evaluated as OCCT 8.1's Geom_OffsetCurve (reference values from it; OCCT
    /// 7.9.3 through OCP agrees but for D3's third derivative, which it got wrong another way,
    /// and DN above the third order, which it refused), the quirks of its D3 (a second
    /// derivative unlike D2's) and of the null-tangent ends (AdjustDerivative, turned over at
    /// the far end) included.
    #[test]
    fn offset_curve_matches_occt() {
        use super::super::curve::BSplineCurve;
        let (regular, singular) = offset_test_poles();
        let basis = BSplineCurve::new(regular, vec![0.0, 1.0, 2.0, 3.0], vec![6, 1, 1, 6], 5);
        let curve = OffsetCurve::new(Curve::BSpline(basis), 0.3, Direction::from_vector(Vector::new(0.0, 0.0, 1.0)));
        let regular_reference: [OffsetReference<3>; 4] = [
            (0.3, [1.036512296071307, 0.6143789417766173, 0.632443125], [1.7893588619915195, -0.7175599594789732, 1.44634375], [8.138715986400722, -11.823607633941705, -5.473750000000001], [6.611080814950514, -15.633028011730305, -5.473750000000001], [-1482.5014437103234, -3561.6426319977277, -3.9875000000000056], [-32.24999999999999, -424.08333333333337, 78.41666666666664], [42.49999999999999, 636.3888888888889, -166.38888888888886]),
            (1.0, [2.840948058501177, 0.28369199048109606, 0.6273148148148151], [2.062019116263083, 1.7655196354932279, -0.5439814814814812], [-2.933943110764039, 2.0543401776254697, 1.435185185185185], [-5.087350512490628, 4.569388517046599, 1.435185185185185], [104.4988749086049, -125.45826232271351, 10.138888888888888], [-2.5000000000000004, 21.38888888888889, -38.05555555555555], [5.000000000000001, 36.38888888888888, 32.361111111111114]),
            (1.7, [3.881647381532218, 1.1680625316261035, 0.8423649421296299], [0.93889480624479, 0.14200294541478495, 1.0929126157407398], [-0.2819183073930448, -2.4826138503159685, 1.0587731481481466], [0.334300228464488, -6.556926401639766, 1.0587731481481466], [26.430372788986606, -0.8946041408649976, -8.571527777777776], [0.9999999999999998, 46.8611111111111, -15.402777777777775], [5.000000000000001, 36.38888888888888, 32.361111111111114]),
            (2.6, [5.676073656932933, 0.9966444833107286, 1.0924562962962954], [2.574507846498388, 0.8935546406476844, -1.2379259259259263], [-1.1322765021680574, -0.2996084299722732, -1.1318518518518494], [-1.9628089030237286, 2.0933198650721234, -1.1318518518518494], [48.111415892866155, -179.6727413282667, 29.822222222222226], [28.000000000000004, -248.88888888888897, 144.2222222222222], [42.49999999999999, -511.11111111111114, 249.8611111111111]),
        ];
        for (t, p, d1, d2, d3_d2, d3, d4, d5) in regular_reference {
            let at = format!("t = {t}");
            assert_close(xyz(curve.point(t) - Point::ORIGIN), p, &at);
            let (q, e1) = curve.d1(t);
            assert_close(xyz(q - Point::ORIGIN), p, &at);
            assert_close(xyz(e1), d1, &at);
            let (q, e1, e2) = curve.d2(t);
            assert_close(xyz(q - Point::ORIGIN), p, &at);
            assert_close(xyz(e1), d1, &at);
            assert_close(xyz(e2), d2, &at);
            let (q, e1, e2, e3) = curve.d3(t);
            assert_close(xyz(q - Point::ORIGIN), p, &at);
            assert_close(xyz(e1), d1, &at);
            assert_close(xyz(e2), d3_d2, &at);
            assert_close(xyz(e3), d3, &at);
            for (n, want) in [(1, d1), (2, d2), (3, d3), (4, d4), (5, d5)] {
                assert_close(xyz(curve.dn(t, n)), want, &format!("{at}, DN({n})"));
            }
        }

        let basis = BSplineCurve::new(singular, vec![0.0, 1.0, 2.0], vec![4, 1, 4], 3);
        let curve = OffsetCurve::new(Curve::BSpline(basis), 0.3, Direction::from_vector(Vector::new(0.0, 0.0, 1.0)));
        let singular_reference: [OffsetReference<3>; 2] = [
            (0.0, [0.2683281572999747, -0.13416407864998736, 0.0], [2.9329179606750064, 5.865835921350013, 1.5], [-1.768328157299975, -10.868951216287465, 0.0], [0.2683281572999747, -0.13416407864998736, 0.0], [2.9329179606750064, 5.865835921350013, 1.5], [-4.04911749434976, -9.728556547762572, 0.0], [17.860592970279573, -6.540498834186885, 0.0]),
            (2.0, [2.865835921350013, 0.7316718427000253, 2.0], [5.865835921350013, -2.932917960675006, 4.5], [-7.063966744387541, 7.365835921350013, -6.0], [2.865835921350013, 0.7316718427000253, 2.0], [5.865835921350013, -2.932917960675006, 4.5], [9.076427924137352, -5.353374741600202, 6.0], [6.741744952161866, 8.70389460241793, 0.0]),
        ];
        for (t, p, d1, d2, p3, d3_d1, d3_d2, d3) in singular_reference {
            let at = format!("t = {t}");
            let (q, e1, e2) = curve.d2(t);
            assert_close(xyz(q - Point::ORIGIN), p, &at);
            assert_close(xyz(e1), d1, &at);
            assert_close(xyz(e2), d2, &at);
            let (q, e1, e2, e3) = curve.d3(t);
            assert_close(xyz(q - Point::ORIGIN), p3, &at);
            assert_close(xyz(e1), d3_d1, &at);
            assert_close(xyz(e2), d3_d2, &at);
            assert_close(xyz(e3), d3, &at);
        }
    }

    /// `GeomAdaptor_Curve::DN(U, 3)` of an offset curve, OCCT 8.1's (its second derivative
    /// standing for the first: [`offset_adaptor_third_derivative`]).
    #[test]
    fn offset_curve_adaptor_third_derivative_matches_occt() {
        use super::super::curve::BSplineCurve;
        let (regular, _) = offset_test_poles();
        let basis = BSplineCurve::new(regular, vec![0.0, 1.0, 2.0, 3.0], vec![6, 1, 1, 6], 5);
        let curve = OffsetCurve::new(Curve::BSpline(basis), 0.3, Direction::from_vector(Vector::new(0.0, 0.0, 1.0)));
        let occt = [
            (0.3, [2909.273253710383, -904.90727782838303, -3.9874999999999936]),
            (1.0, [168.80088003012858, -2.5473493712513715, 10.138888888888886]),
            (1.7, [-298.26987879065734, -6.9122246075663973, -8.5715277777777779]),
            (2.6, [86.086461937931446, -271.06461408439219, 29.82222222222223]),
        ];
        for (t, want) in occt {
            assert_close(xyz(offset_adaptor_third_derivative(&curve, t).unwrap()), want, &format!("t = {t}"));
        }
    }

    /// As [`offset_curve_matches_occt`], for Geom2d_OffsetCurve.
    #[test]
    fn offset_curve2_matches_occt() {
        use super::super::curve::BSplineCurve;
        let (regular, singular) = offset_test_poles();
        let flat = |poles: Vec<Point>| poles.into_iter().map(|p| Point2::new(p.x, p.y)).collect::<Vec<_>>();
        let basis = BSplineCurve::new(flat(regular), vec![0.0, 1.0, 2.0, 3.0], vec![6, 1, 1, 6], 5);
        let curve = OffsetCurve2::new(Curve2::BSpline(basis), 0.3);
        let regular_reference: [OffsetReference<2>; 4] = [
            (0.3, [1.036512296071307, 0.6143789417766173], [1.7893588619915195, -0.7175599594789732], [8.138715986400722, -11.823607633941705], [6.611080814950515, -15.633028011730305], [-1482.5014437103234, -3561.6426319977277], [-32.24999999999999, -424.08333333333337], [42.49999999999999, 636.3888888888889]),
            (1.0, [2.840948058501177, 0.28369199048109606], [2.062019116263083, 1.7655196354932279], [-2.933943110764039, 2.0543401776254697], [-5.0873505124906275, 4.569388517046599], [104.49887490860489, -125.45826232271352], [-2.5000000000000004, 21.38888888888889], [5.000000000000001, 36.38888888888888]),
            (1.7, [3.881647381532218, 1.1680625316261035], [0.93889480624479, 0.14200294541478495], [-0.2819183073930448, -2.4826138503159685], [0.33430022846448804, -6.556926401639766], [26.4303727889866, -0.8946041408649976], [0.9999999999999998, 46.8611111111111], [5.000000000000001, 36.38888888888888]),
            (2.6, [5.676073656932933, 0.9966444833107286], [2.574507846498388, 0.8935546406476844], [-1.1322765021680574, -0.2996084299722732], [-1.9628089030237295, 2.0933198650721234], [48.11141589286615, -179.67274132826668], [28.000000000000004, -248.88888888888897], [42.49999999999999, -511.11111111111114]),
        ];
        for (t, p, d1, d2, d3_d2, d3, d4, d5) in regular_reference {
            let at = format!("t = {t}");
            assert_close(xy(curve.point(t) - Point2::ORIGIN), p, &at);
            let (q, e1) = curve.d1(t);
            assert_close(xy(q - Point2::ORIGIN), p, &at);
            assert_close(xy(e1), d1, &at);
            let (q, e1, e2) = curve.d2(t);
            assert_close(xy(q - Point2::ORIGIN), p, &at);
            assert_close(xy(e1), d1, &at);
            assert_close(xy(e2), d2, &at);
            let (q, e1, e2, e3) = curve.d3(t);
            assert_close(xy(q - Point2::ORIGIN), p, &at);
            assert_close(xy(e1), d1, &at);
            assert_close(xy(e2), d3_d2, &at);
            assert_close(xy(e3), d3, &at);
            for (n, want) in [(1, d1), (2, d2), (3, d3), (4, d4), (5, d5)] {
                assert_close(xy(curve.dn(t, n)), want, &format!("{at}, DN({n})"));
            }
        }

        let basis = BSplineCurve::new(flat(singular), vec![0.0, 1.0, 2.0], vec![4, 1, 4], 3);
        let curve = OffsetCurve2::new(Curve2::BSpline(basis), 0.3);
        let singular_reference: [OffsetReference<2>; 2] = [
            (0.0, [0.2683281572999747, -0.13416407864998736], [2.9329179606750064, 5.865835921350013], [-1.768328157299975, -10.868951216287465], [0.2683281572999747, -0.13416407864998736], [2.9329179606750064, 5.865835921350013], [-4.04911749434976, -9.728556547762572], [17.86059297027957, -6.540498834186884]),
            (2.0, [2.865835921350013, 0.7316718427000253], [5.865835921350013, -2.932917960675006], [-7.063966744387541, 7.365835921350013], [2.865835921350013, 0.7316718427000253], [5.865835921350013, -2.932917960675006], [9.076427924137352, -5.353374741600202], [6.7417449521618655, 8.70389460241793]),
        ];
        for (t, p, d1, d2, p3, d3_d1, d3_d2, d3) in singular_reference {
            let at = format!("t = {t}");
            let (q, e1, e2) = curve.d2(t);
            assert_close(xy(q - Point2::ORIGIN), p, &at);
            assert_close(xy(e1), d1, &at);
            assert_close(xy(e2), d2, &at);
            let (q, e1, e2, e3) = curve.d3(t);
            assert_close(xy(q - Point2::ORIGIN), p3, &at);
            assert_close(xy(e1), d3_d1, &at);
            assert_close(xy(e2), d3_d2, &at);
            assert_close(xy(e3), d3, &at);
        }
    }
}
