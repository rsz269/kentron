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
pub(crate) mod abscissa;
pub(crate) mod adaptor;
pub(crate) mod curvlin;
pub(crate) mod function_root;
pub(crate) mod ext_pc;
pub(crate) mod locate_ext_pc;
// (OCCT's shape analysis and pcurve projection, as the readers add pcurves.)
#[allow(dead_code)]
pub(crate) mod surface_analysis;
#[allow(dead_code)]
pub(crate) mod pcurve_projection;
pub(crate) mod same_parameter;
pub(crate) mod function_roots;
// (The ports below serve the guide-plane trihedron, not wired in yet.)
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) mod int_curve_plane;
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) mod polynomial_roots;
// (For the shape-analysis ports, not wired in yet.)
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) mod resolution;
pub(crate) mod tangential_deflection;
pub(crate) mod brent_minimum;
pub(crate) mod particle_swarm;
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
pub(crate) use bspline::basis_and_first;
pub use control::ControlPoint;
pub use parametric::{Continuity, ParametricCurve, ParametricSurface, SurfaceD2, SurfaceD3};
pub(crate) use patches::SurfaceId;
