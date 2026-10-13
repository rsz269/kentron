use std::ops::{Add, Sub};

use super::{Vector, Vector2};

/// A position in 3D space.
#[doc(alias = "gp_Pnt")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub const ORIGIN: Point = Point::new(0.0, 0.0, 0.0);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn distance(&self, other: Point) -> f64 {
        self.distance_squared(other).sqrt()
    }

    pub fn distance_squared(&self, other: Point) -> f64 {
        (*self - other).length_squared()
    }

    /// True if the distance to `other` is at most `tolerance`.
    pub fn approx_eq(&self, other: Point, tolerance: f64) -> bool {
        self.distance(other) <= tolerance
    }

    /// The weighted average of `points`; panics if the weights sum to zero.
    pub fn barycenter(points: &[(Point, f64)]) -> Point {
        let total: f64 = points.iter().map(|(_, w)| w).sum();
        assert!(total != 0.0, "barycenter of zero total weight");
        let sum = points
            .iter()
            .fold(Vector::ZERO, |acc, &(p, w)| acc + (p - Point::ORIGIN) * w);
        Point::ORIGIN + sum / total
    }
}

impl Add<Vector> for Point {
    type Output = Point;
    fn add(self, v: Vector) -> Point {
        Point::new(self.x + v.x, self.y + v.y, self.z + v.z)
    }
}

impl Sub<Vector> for Point {
    type Output = Point;
    fn sub(self, v: Vector) -> Point {
        Point::new(self.x - v.x, self.y - v.y, self.z - v.z)
    }
}

/// The vector from `other` to `self`.
impl Sub for Point {
    type Output = Vector;
    fn sub(self, other: Point) -> Vector {
        Vector::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

/// A position in 2D space (also used for (u, v) surface parameters).
#[doc(alias = "gp_Pnt2d")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    pub const ORIGIN: Point2 = Point2::new(0.0, 0.0);

    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn distance(&self, other: Point2) -> f64 {
        self.distance_squared(other).sqrt()
    }

    pub fn distance_squared(&self, other: Point2) -> f64 {
        (*self - other).length_squared()
    }

    pub fn approx_eq(&self, other: Point2, tolerance: f64) -> bool {
        self.distance(other) <= tolerance
    }
}

impl Add<Vector2> for Point2 {
    type Output = Point2;
    fn add(self, v: Vector2) -> Point2 {
        Point2::new(self.x + v.x, self.y + v.y)
    }
}

/// The vector from `other` to `self`.
impl Sub for Point2 {
    type Output = Vector2;
    fn sub(self, other: Point2) -> Vector2 {
        Vector2::new(self.x - other.x, self.y - other.y)
    }
}
