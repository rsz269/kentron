use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

use super::{Direction, Direction2};

/// A displacement in 3D space (also used as a plain coordinate triple).
#[doc(alias = "gp_Vec")]
#[doc(alias = "gp_XYZ")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub const ZERO: Vector = Vector::new(0.0, 0.0, 0.0);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn dot(&self, other: Vector) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: Vector) -> Vector {
        Vector::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    pub fn length_squared(&self) -> f64 {
        self.dot(*self)
    }

    /// Unit vector in the same direction; panics if the length is (near) zero.
    pub fn normalized(&self) -> Vector {
        let length = self.length();
        assert!(length > crate::tolerance::RESOLUTION, "normalizing a zero-length vector");
        *self / length
    }

    /// Angle in [0, PI] between the two vectors.
    pub fn angle(&self, other: Vector) -> f64 {
        self.cross(other).length().atan2(self.dot(other))
    }

    pub fn is_parallel(&self, other: Vector, angular_tolerance: f64) -> bool {
        let angle = self.angle(other);
        angle <= angular_tolerance || std::f64::consts::PI - angle <= angular_tolerance
    }

    pub fn is_opposite(&self, other: Vector, angular_tolerance: f64) -> bool {
        std::f64::consts::PI - self.angle(other) <= angular_tolerance
    }

    pub fn is_perpendicular(&self, other: Vector, angular_tolerance: f64) -> bool {
        (std::f64::consts::FRAC_PI_2 - self.angle(other)).abs() <= angular_tolerance
    }

    /// True if every coordinate differs by at most `tolerance`.
    pub fn approx_eq(&self, other: Vector, tolerance: f64) -> bool {
        (self.x - other.x).abs() <= tolerance
            && (self.y - other.y).abs() <= tolerance
            && (self.z - other.z).abs() <= tolerance
    }
}

impl From<Direction> for Vector {
    fn from(d: Direction) -> Self {
        Vector::new(d.x(), d.y(), d.z())
    }
}

/// Coordinate by index: 0 = x, 1 = y, 2 = z. Panics on other indices.
impl Index<usize> for Vector {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        match i {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("vector index {i} out of range"),
        }
    }
}

impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        match i {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("vector index {i} out of range"),
        }
    }
}

impl Add for Vector {
    type Output = Vector;
    fn add(self, rhs: Vector) -> Vector {
        Vector::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vector {
    type Output = Vector;
    fn sub(self, rhs: Vector) -> Vector {
        Vector::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Vector {
    type Output = Vector;
    fn mul(self, rhs: f64) -> Vector {
        Vector::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Vector {
        rhs * self
    }
}

impl Div<f64> for Vector {
    type Output = Vector;
    fn div(self, rhs: f64) -> Vector {
        Vector::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Neg for Vector {
    type Output = Vector;
    fn neg(self) -> Vector {
        Vector::new(-self.x, -self.y, -self.z)
    }
}

/// A displacement in 2D space.
#[doc(alias = "gp_Vec2d")]
#[doc(alias = "gp_XY")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

impl Vector2 {
    pub const ZERO: Vector2 = Vector2::new(0.0, 0.0);

    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn dot(&self, other: Vector2) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The scalar cross product `x1*y2 - y1*x2`.
    pub fn cross(&self, other: Vector2) -> f64 {
        self.x * other.y - self.y * other.x
    }

    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    pub fn length_squared(&self) -> f64 {
        self.dot(*self)
    }

    /// Panics if the length is (near) zero.
    pub fn normalized(&self) -> Vector2 {
        let length = self.length();
        assert!(length > crate::tolerance::RESOLUTION, "normalizing a zero-length vector");
        *self / length
    }

    /// Signed angle in [-PI, PI] from `self` to `other`.
    pub fn angle(&self, other: Vector2) -> f64 {
        self.cross(other).atan2(self.dot(other))
    }

    pub fn is_parallel(&self, other: Vector2, angular_tolerance: f64) -> bool {
        let angle = self.angle(other).abs();
        angle <= angular_tolerance || std::f64::consts::PI - angle <= angular_tolerance
    }

    pub fn is_opposite(&self, other: Vector2, angular_tolerance: f64) -> bool {
        std::f64::consts::PI - self.angle(other).abs() <= angular_tolerance
    }

    pub fn is_perpendicular(&self, other: Vector2, angular_tolerance: f64) -> bool {
        (std::f64::consts::FRAC_PI_2 - self.angle(other).abs()).abs() <= angular_tolerance
    }
}

impl From<Direction2> for Vector2 {
    fn from(d: Direction2) -> Self {
        Vector2::new(d.x(), d.y())
    }
}

impl Add for Vector2 {
    type Output = Vector2;
    fn add(self, rhs: Vector2) -> Vector2 {
        Vector2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Vector2 {
    type Output = Vector2;
    fn sub(self, rhs: Vector2) -> Vector2 {
        Vector2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Vector2 {
    type Output = Vector2;
    fn mul(self, rhs: f64) -> Vector2 {
        Vector2::new(self.x * rhs, self.y * rhs)
    }
}

impl Mul<Vector2> for f64 {
    type Output = Vector2;
    fn mul(self, rhs: Vector2) -> Vector2 {
        rhs * self
    }
}

impl Div<f64> for Vector2 {
    type Output = Vector2;
    fn div(self, rhs: f64) -> Vector2 {
        Vector2::new(self.x / rhs, self.y / rhs)
    }
}

impl Neg for Vector2 {
    type Output = Vector2;
    fn neg(self) -> Vector2 {
        Vector2::new(-self.x, -self.y)
    }
}
