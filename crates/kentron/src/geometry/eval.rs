//! Evaluation of the elementary curves and surfaces, and dispatch for the `Curve`,
//! `Curve2` and `Surface` enums.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

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
        self.points_along_v_of(u, u.sin_cos(), vs)
    }

    /// [`Self::points_along_v`], given `u.sin_cos()` as `su`.
    pub(crate) fn points_along_v_of<const N: usize>(&self, u: f64, su: (f64, f64), vs: [f64; N]) -> [Point; N] {
        match self {
            Surface::Cylinder(c) => vs.map(|v| c.point_of(su, v)),
            Surface::Cone(c) => vs.map(|v| c.point_of(su, v)),
            _ => vs.map(|v| self.point(u, v)),
        }
    }

    /// The point at `u` and `v`, given `u.sin_cos()` as `su` (the same as [`ParametricSurface::point`]).
    pub(crate) fn point_of_u(&self, u: f64, su: (f64, f64), v: f64) -> Point {
        match self {
            Surface::Cylinder(c) => c.point_of(su, v),
            Surface::Cone(c) => c.point_of(su, v),
            _ => self.point(u, v),
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

// --- offset curves --------------------------------------------------------------------------

/// The offset direction's unnormalised vector `N` from a basis derivative: `B′ × D` for a 3D
/// offset along `direction`; in 2D the derivative turned clockwise, `(B′_y, −B′_x)`. Linear,
/// so `N⁽ʲ⁾` is the map of `B⁽ʲ⁺¹⁾`.
fn offset_normal_3d(direction: [f64; 3]) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |b| [b[1] * direction[2] - b[2] * direction[1], b[2] * direction[0] - b[0] * direction[2], b[0] * direction[1] - b[1] * direction[0]]
}

fn offset_normal_2d(b: [f64; 3]) -> [f64; 3] {
    [b[1], -b[0], 0.0]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `C(t) = B(t) + o·n(t)` and its derivatives of orders `0..=order` (as coordinates), for
/// the basis `basis`, offset `o` and normal map `normal` (see [`offset_normal_3d`]), with
/// `n = N/|N|`. Where the basis tangent vanishes, the limit (see
/// [`offset_basis_derivatives`]). Non-finite where the offset is undefined.
fn offset_derivatives<C>(basis: &C, offset: f64, normal: impl Fn([f64; 3]) -> [f64; 3], t: f64, order: usize) -> Vec<[f64; 3]>
where
    C: ParametricCurve,
    C::Point: super::control::ControlPoint<Vector = C::Vector>,
{
    let d = offset_basis_derivatives(basis, t, order);
    // N⁽ʲ⁾, j = 0..=order.
    let nd: Vec<[f64; 3]> = d[1..].iter().map(|&b| normal(b)).collect();
    // g = N·N and its derivatives by Leibniz' rule.
    let binomial = |m: usize, i: usize| -> f64 { (0..i).fold(1.0, |acc, j| acc * (m - j) as f64 / (j + 1) as f64) };
    let g: Vec<f64> = (0..=order).map(|m| (0..=m).map(|i| binomial(m, i) * dot3(nd[i], nd[m - i])).sum()).collect();
    // s = g^(−1/2): from g s′ = −½ g′ s, differentiated m times.
    let mut s = vec![1.0 / g[0].sqrt()];
    for m in 0..order {
        let mut rhs: f64 = -0.5 * (0..=m).map(|i| binomial(m, i) * g[i + 1] * s[m - i]).sum::<f64>();
        rhs -= (1..=m).map(|i| binomial(m, i) * g[i] * s[m + 1 - i]).sum::<f64>();
        s.push(rhs / g[0]);
    }
    (0..=order)
        .map(|k| {
            let mut n = [0.0; 3];
            for i in 0..=k {
                let c = binomial(k, i) * s[k - i];
                for x in 0..3 {
                    n[x] += c * nd[i][x];
                }
            }
            std::array::from_fn(|x| d[k][x] + offset * n[x])
        })
        .collect()
}

/// The basis point and the derivatives standing for `B′ … B⁽ᵒʳᵈᵉʳ⁺¹⁾` in the offset formulas
/// at `t`, as coordinates. Where the basis tangent vanishes (`|B′|²` not above
/// `RESOLUTION`) the first basis derivative of order `k = 2` or `3` that does not vanish
/// and those after it, `±B⁽ᵏ⁾, ±B⁽ᵏ⁺¹⁾, …`, signed by the direction of travel on the side
/// of larger parameters at the range's start and of smaller parameters elsewhere; NaN if
/// all vanish.
fn offset_basis_derivatives<C>(basis: &C, t: f64, order: usize) -> Vec<[f64; 3]>
where
    C: ParametricCurve,
    C::Point: super::control::ControlPoint<Vector = C::Vector>,
{
    use super::control::ControlPoint;
    let v = |x: C::Vector| C::Point::vector_coords(x);
    let p = |x: C::Point| x.coords();
    let mut d: Vec<[f64; 3]> = match order {
        0 => {
            let (b, d1) = basis.d1(t);
            vec![p(b), v(d1)]
        }
        1 => {
            let (b, d1, d2) = basis.d2(t);
            vec![p(b), v(d1), v(d2)]
        }
        _ => {
            let (b, d1, d2, d3) = basis.d3(t);
            vec![p(b), v(d1), v(d2), v(d3)]
        }
    };
    for n in d.len()..order + 2 {
        d.push(v(basis.dn(t, n as u32)));
    }
    if dot3(d[1], d[1]) > RESOLUTION {
        return d;
    }
    // A vanishing tangent: the same basis values whatever the order asked for.
    let (b, d1, d2, d3) = basis.d3(t);
    let mut all = vec![p(b), v(d1), v(d2), v(d3)];
    let Some(k) = (2..=3).find(|&k| dot3(all[k], all[k]) > RESOLUTION) else {
        return vec![[f64::NAN; 3]; order + 2];
    };
    for n in all.len()..=k + order {
        all.push(v(basis.dn(t, n as u32)));
    }
    let (a, z) = (basis.first_parameter(), basis.last_parameter());
    let at_start = (t - a) <= 1e-9 * (z - a).abs().max(f64::MIN_POSITIVE) && (z - a).is_finite();
    // Near t₀, B′(t) ≈ B⁽ᵏ⁾(t₀)·(t − t₀)^(k−1)/(k−1)!: its sign on the chosen side.
    let sign = if at_start || k % 2 == 1 { 1.0 } else { -1.0 };
    let mut out = vec![all[0]];
    out.extend((0..=order).map(|j| all[k + j].map(|x| sign * x)));
    out
}

impl OffsetCurve {
    fn derivatives(&self, t: f64, order: usize) -> Vec<[f64; 3]> {
        let direction = Vector::from(self.direction());
        offset_derivatives(self.basis(), self.offset(), offset_normal_3d([direction.x, direction.y, direction.z]), t, order)
    }
}

impl OffsetCurve2 {
    fn derivatives(&self, t: f64, order: usize) -> Vec<[f64; 3]> {
        offset_derivatives(self.basis(), self.offset(), offset_normal_2d, t, order)
    }
}

/// An offset curve's continuity, one order below its basis's (a C0 B-spline basis whose
/// pieces meet tangentially over `first..last` gives G1).
fn offset_continuity(basis: Continuity, spline_g1: impl FnOnce() -> bool) -> Continuity {
    match basis {
        Continuity::CN => Continuity::CN,
        Continuity::C3 => Continuity::C2,
        Continuity::C2 => Continuity::C1,
        Continuity::G2 => Continuity::G2,
        Continuity::G1 => Continuity::G1,
        Continuity::C1 => Continuity::C0,
        Continuity::C0 => if spline_g1() { Continuity::G1 } else { Continuity::C0 },
    }
}

macro_rules! offset_curve {
    ($ty:ident, $curve:ident, $point:ty, $vector:ty, $make:expr, $negated:expr) => {
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
                let basis = self.basis();
                offset_continuity(basis.continuity(), || {
                    let angular = crate::tolerance::ANGULAR;
                    match basis {
                        $curve::BSpline(s) => s.is_g1(s.first_parameter(), s.last_parameter(), angular),
                        $curve::Trimmed(t) => match t.basis() {
                            $curve::BSpline(s) => s.is_g1(t.first_parameter(), t.last_parameter(), angular),
                            _ => false,
                        },
                        _ => false,
                    }
                })
            }
            fn is_cn(&self, n: u32) -> bool {
                self.basis().is_cn(n + 1)
            }
            fn point(&self, t: f64) -> $point {
                let d = self.derivatives(t, 0);
                <$point as super::control::ControlPoint>::from_coords(d[0])
            }
            fn d1(&self, t: f64) -> ($point, $vector) {
                let d = self.derivatives(t, 1);
                let v = <$point as super::control::ControlPoint>::vector_from_coords;
                (<$point as super::control::ControlPoint>::from_coords(d[0]), v(d[1]))
            }
            fn d2(&self, t: f64) -> ($point, $vector, $vector) {
                let d = self.derivatives(t, 2);
                let v = <$point as super::control::ControlPoint>::vector_from_coords;
                (<$point as super::control::ControlPoint>::from_coords(d[0]), v(d[1]), v(d[2]))
            }
            fn d3(&self, t: f64) -> ($point, $vector, $vector, $vector) {
                let d = self.derivatives(t, 3);
                let v = <$point as super::control::ControlPoint>::vector_from_coords;
                (<$point as super::control::ControlPoint>::from_coords(d[0]), v(d[1]), v(d[2]), v(d[3]))
            }
            fn dn(&self, t: f64, n: u32) -> $vector {
                assert!(n > 0, "derivative of order 0");
                // (Orders 1 to 3 evaluated as by d1 to d3.)
                <$point as super::control::ControlPoint>::vector_from_coords(self.derivatives(t, n as usize)[n as usize])
            }
            fn reversed_parameter(&self, t: f64) -> f64 {
                self.basis().reversed_parameter(t)
            }
            fn reverse(&mut self) {
                let mut basis = self.basis().clone();
                basis.reverse();
                *self = $make(self, basis, $negated(self));
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

offset_curve!(OffsetCurve, Curve, Point, Vector, |c: &OffsetCurve, basis, offset| OffsetCurve::new(basis, offset, c.direction()), |c: &OffsetCurve| -c.offset());
offset_curve!(OffsetCurve2, Curve2, Point2, Vector2, |_: &OffsetCurve2, basis, offset| OffsetCurve2::new(basis, offset), |c: &OffsetCurve2| -c.offset());

/// The third derivative of an offset curve for the curve adaptor: the exact one, `None`
/// where the offset is undefined (the adaptor then asks the curve itself).
pub(crate) fn offset_adaptor_third_derivative(c: &OffsetCurve, t: f64) -> Option<Vector> {
    let d3 = c.d3(t).3;
    (d3.x.is_finite() && d3.y.is_finite() && d3.z.is_finite()).then_some(d3)
}

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

    use super::super::bspline::BSplineCurve;
    use super::super::curve2::BSplineCurve2;

    /// A cubic whose first two and last two poles coincide: its tangent vanishes at both
    /// ends.
    fn stalled_cubic() -> BSplineCurve {
        let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new(1.0, 2.0, 0.0), Point::new(3.0, 1.0, 0.0), Point::new(4.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0)];
        BSplineCurve::new(poles, vec![0.0, 1.0, 2.0], vec![4, 2, 4], 3)
    }

    fn quintic() -> BSplineCurve {
        let poles: Vec<Point> = (0..8).map(|i| Point::new(i as f64, (0.9 * i as f64).sin(), 0.2 * (0.5 * i as f64).cos())).collect();
        BSplineCurve::new(poles, vec![0.0, 1.0, 1.8, 3.0], vec![6, 1, 1, 6], 5)
    }

    fn close(a: Vector, b: Vector, relative: f64) -> bool {
        (a - b).length() <= relative * (1.0 + b.length())
    }

    /// Every derivative is the central difference of the one below it, and the orders agree
    /// across d1, d2, d3 and dn.
    fn check_derivatives(c: &Curve, ts: &[f64]) {
        let h = 1e-5;
        for &t in ts {
            let (p0, a1) = c.d1(t);
            let (p1, b1, b2) = c.d2(t);
            let (p2, c1, c2, c3) = c.d3(t);
            assert!(p0.distance(c.point(t)) < 1e-13 && p1.distance(p0) < 1e-13 && p2.distance(p0) < 1e-13);
            assert!(close(a1, b1, 1e-13) && close(c1, b1, 1e-13) && close(c.dn(t, 1), a1, 1e-13));
            assert_eq!(c2, b2, "the second derivative of d3 is d2's");
            assert!(close(c.dn(t, 2), b2, 1e-13) && close(c.dn(t, 3), c3, 1e-13));
            let difference = |k: u32| {
                let f = |s: f64| if k == 1 { c.point(s) - Point::ORIGIN } else { c.dn(s, k - 1) };
                (f(t + h) - f(t - h)) * (0.5 / h)
            };
            for k in 1..=5 {
                assert!(close(c.dn(t, k), difference(k), 1e-6), "order {k} at {t}: {:?} {:?}", c.dn(t, k), difference(k));
            }
        }
    }

    #[test]
    fn offset_curves_have_exact_derivatives() {
        let offset = Curve::Offset(OffsetCurve::new(Curve::BSpline(quintic()), 0.3, Direction::Z));
        check_derivatives(&offset, &[0.3, 1.4, 2.6]);
        let stalled = Curve::Offset(OffsetCurve::new(Curve::BSpline(stalled_cubic()), 0.3, Direction::Z));
        check_derivatives(&stalled, &[0.3, 0.6, 1.7]);
        // The adaptor's third derivative is the curve's.
        let Curve::Offset(o) = &offset else { unreachable!() };
        assert_eq!(offset_adaptor_third_derivative(o, 0.7), Some(o.d3(0.7).3));
    }

    #[test]
    fn offsets_of_a_circle_are_circles() {
        let r = 2.0;
        for (o, sign) in [(0.5, 1.0), (-0.5, -1.0)] {
            // B′ × Z points away from the centre of a counterclockwise circle.
            let c = OffsetCurve::new(Curve::Circle(Circle::new(Frame::WORLD, r)), o, Direction::Z);
            let radius = r + sign * 0.5;
            for t in [0.0, 0.7, 2.0, 4.5] {
                assert!((c.point(t) - Point::new(radius * t.cos(), radius * t.sin(), 0.0)).length() < 1e-14);
                for n in 1..=5u32 {
                    let (dc, ds) = cos_sin_derivative(t, n);
                    assert!(close(c.dn(t, n), Vector::new(radius * dc, radius * ds, 0.0), 1e-13), "order {n}");
                }
            }
        }
        // 2D: to the right of a counterclockwise circle is outwards.
        let c2 = OffsetCurve2::new(Curve2::Circle(Circle2::new(crate::math::Frame2::WORLD, r)), 0.5);
        for t in [0.0, 1.3, 3.0] {
            assert!((c2.point(t) - Point2::new(2.5 * t.cos(), 2.5 * t.sin())).length() < 1e-14);
            let (dc, ds) = cos_sin_derivative(t, 3);
            assert!((c2.dn(t, 3) - Vector2::new(2.5 * dc, 2.5 * ds)).length() < 1e-13);
        }
    }

    #[test]
    fn offsets_of_lines_and_parabolas() {
        let line = OffsetCurve::new(Curve::Line(Line::new(Point::ORIGIN, Direction::X)), 2.0, Direction::Z);
        assert!((line.point(3.0) - Point::new(3.0, -2.0, 0.0)).length() < 1e-15);
        for n in 2..=5 {
            assert!(line.dn(1.0, n).length() < 1e-15);
        }
        // The parabola (t, t²) offset by o to its right: closed form by the unit normal
        // (2t, −1)/√(1 + 4t²).
        let poles = vec![Point2::new(-1.0, 1.0), Point2::new(0.0, -1.0), Point2::new(1.0, 1.0)];
        let parabola = Curve2::BSpline(BSplineCurve2::new(poles, vec![-1.0, 1.0], vec![3, 3], 2));
        let o = 0.25;
        let c = OffsetCurve2::new(parabola, o);
        let exact = |t: f64| {
            let q = (1.0 + 4.0 * t * t).sqrt();
            Point2::new(t + o * 2.0 * t / q, t * t - o / q)
        };
        for t in [-0.6, 0.1, 0.8] {
            assert!((c.point(t) - exact(t)).length() < 1e-14);
            let q = (1.0 + 4.0 * t * t).sqrt();
            let d1 = Vector2::new(1.0 + o * 2.0 / q.powi(3), 2.0 * t + o * 4.0 * t / q.powi(3));
            assert!((c.d1(t).1 - d1).length() < 1e-13, "{:?} {d1:?}", c.d1(t).1);
            let h = 1e-4;
            for k in 2..=5u32 {
                let fd = (c.dn(t + h, k - 1) - c.dn(t - h, k - 1)) * (0.5 / h);
                assert!((c.dn(t, k) - fd).length() <= 1e-5 * (1.0 + fd.length()), "order {k}");
            }
        }
    }

    #[test]
    fn offsets_where_the_basis_tangent_vanishes_are_limits() {
        let c = OffsetCurve::new(Curve::BSpline(stalled_cubic()), 0.3, Direction::Z);
        for (t, inward) in [(0.0, 1.0), (2.0, -1.0)] {
            let p = c.point(t);
            assert!(p.x.is_finite() && p.y.is_finite());
            assert_eq!(p, c.d2(t).0);
            assert_eq!(c.d1(t).1, c.d2(t).1);
            assert_eq!(c.d2(t).0, c.d3(t).0);
            // Approached from inside the range, the point tends to the limit.
            let near = c.point(t + inward * 1e-7);
            assert!(near.distance(p) < 1e-5, "{near:?} {p:?}");
            // The first derivative lies along the curve's limit tangent.
            let chord = c.point(t + inward * 1e-6) - p;
            let d1 = c.d1(t).1;
            assert!(chord.cross(d1).length() < 1e-3 * chord.length() * d1.length(), "{chord:?} {d1:?}");
        }
        // A 3D offset along the line's own direction is undefined: no panic.
        let along = OffsetCurve::new(Curve::Line(Line::new(Point::ORIGIN, Direction::Z)), 1.0, Direction::Z);
        assert!(!along.point(1.0).x.is_finite());
        let Curve::Line(_) = along.basis() else { unreachable!() };
        assert!(offset_adaptor_third_derivative(&along, 1.0).is_none());
    }

    #[test]
    fn reversed_offsets_are_the_same_points_backwards() {
        let c = OffsetCurve::new(Curve::BSpline(quintic()), 0.3, Direction::Z);
        let mut r = c.clone();
        r.reverse();
        for t in [0.0, 0.4, 1.5, 3.0] {
            assert!(r.point(r.reversed_parameter(t)).distance(c.point(t)) < 1e-13);
        }
        let poles = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 2.0), Point2::new(3.0, 1.0)];
        let c2 = OffsetCurve2::new(Curve2::BSpline(BSplineCurve2::new(poles, vec![0.0, 1.0], vec![3, 3], 2)), -0.4);
        let mut r2 = c2.clone();
        r2.reverse();
        for t in [0.0, 0.3, 1.0] {
            assert!((r2.point(r2.reversed_parameter(t)) - c2.point(t)).length() < 1e-13);
        }
    }

    #[test]
    fn offset_continuity_is_one_order_below() {
        let c = OffsetCurve::new(Curve::BSpline(quintic()), 0.3, Direction::Z);
        // (A B-spline counts as C3 at most.)
        assert_eq!(c.continuity(), Continuity::C2);
        assert_eq!(OffsetCurve::new(Curve::Circle(Circle::new(Frame::WORLD, 1.0)), 0.1, Direction::Z).continuity(), Continuity::CN);
        // A polyline with its legs in line: C0 but G1.
        let straight = BSplineCurve::new(vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0)], vec![0.0, 1.0, 2.0], vec![2, 1, 2], 1);
        assert_eq!(OffsetCurve::new(Curve::BSpline(straight), 0.1, Direction::Z).continuity(), Continuity::G1);
        assert!(c.is_cn(3) && !c.is_cn(4));
    }
}
