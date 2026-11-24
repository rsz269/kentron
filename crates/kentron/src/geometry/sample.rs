//! Sampling points along curves.

use crate::math::{Point, Vector};

use super::curve::Curve;
use super::ParametricCurve;

/// Points sampled along a curve.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSampling {
    /// Increasing parameters of the samples, including both ends.
    pub parameters: Vec<f64>,
    pub points: Vec<Point>,
    /// The largest distance between the curve and the polyline through the points.
    pub deflection: f64,
}

// CLEAN-ROOM GAP: 165 lines removed; rewrite from the cleanroom spec `quasi_uniform_sampling.md`.
