//! The point types Bezier and B-spline geometry can be built on, so that one
//! implementation serves 2D and 3D curves.

use std::fmt::Debug;

use crate::math::{Point, Point2, Vector, Vector2};

/// A point type usable as a control point: [`Point`] or [`Point2`].
pub trait ControlPoint: Copy + PartialEq + Debug {
    /// The matching vector type.
    type Vector: Copy + PartialEq + Debug;

    /// The coordinates (a 2D point has `z = 0`).
    fn coords(self) -> [f64; 3];

    fn from_coords(c: [f64; 3]) -> Self;

    fn vector_coords(v: Self::Vector) -> [f64; 3];

    fn vector_from_coords(c: [f64; 3]) -> Self::Vector;
}

impl ControlPoint for Point {
    type Vector = Vector;

    fn coords(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    fn from_coords(c: [f64; 3]) -> Self {
        Point::new(c[0], c[1], c[2])
    }

    fn vector_coords(v: Vector) -> [f64; 3] {
        [v.x, v.y, v.z]
    }

    fn vector_from_coords(c: [f64; 3]) -> Vector {
        Vector::new(c[0], c[1], c[2])
    }
}

impl ControlPoint for Point2 {
    type Vector = Vector2;

    fn coords(self) -> [f64; 3] {
        [self.x, self.y, 0.0]
    }

    fn from_coords(c: [f64; 3]) -> Self {
        Point2::new(c[0], c[1])
    }

    fn vector_coords(v: Vector2) -> [f64; 3] {
        [v.x, v.y, 0.0]
    }

    fn vector_from_coords(c: [f64; 3]) -> Vector2 {
        Vector2::new(c[0], c[1])
    }
}
