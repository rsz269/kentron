//! Curves and surfaces.
//!
//! Each kind of curve or surface is a plain value type ([`curve::Circle`],
//! [`surface::Torus`], ...) that knows both its shape and its parametrization. The
//! [`curve::Curve`], [`curve2::Curve2`] and [`surface::Surface`] enums hold any of them,
//! e.g. inside a [`crate::topology::Model`].
//!
//! Evaluation goes through the [`ParametricCurve`] and [`ParametricSurface`] traits.
//! Bezier and B-spline curves are generic over the point type, so one implementation
//! serves 3D (`BezierCurve`, the default) and 2D (`BezierCurve2 = BezierCurve<Point2>`)
//! curves. Control points and knots are indexed from 0.

mod bezier;
pub mod bounds;
mod bspline;
pub mod construct;
pub mod control;
pub mod convert;
pub mod curve;
pub mod curve2;
mod eval;
pub mod extend;
pub mod extrema;
pub mod fit;
pub mod intersect;
mod knots;
pub mod measure;
mod parametric;
mod patches;
pub mod project;
pub mod props;
pub mod sample;
pub mod surface;
pub mod tangency;

pub use bspline::KnotDistribution;
pub use control::ControlPoint;
pub use parametric::{Continuity, ParametricCurve, ParametricSurface, SurfaceD2, SurfaceD3};
pub(crate) use patches::SurfaceId;
