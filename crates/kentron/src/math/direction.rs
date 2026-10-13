use std::ops::Neg;

use super::{Vector, Vector2};

/// A unit vector in 3D space. Always normalized, so its coordinates are read-only.
#[doc(alias = "gp_Dir")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction {
    x: f64,
    y: f64,
    z: f64,
}

impl Direction {
    pub const X: Direction = Direction { x: 1.0, y: 0.0, z: 0.0 };
    pub const Y: Direction = Direction { x: 0.0, y: 1.0, z: 0.0 };
    pub const Z: Direction = Direction { x: 0.0, y: 0.0, z: 1.0 };
    pub const NEG_X: Direction = Direction { x: -1.0, y: 0.0, z: 0.0 };
    pub const NEG_Y: Direction = Direction { x: 0.0, y: -1.0, z: 0.0 };
    pub const NEG_Z: Direction = Direction { x: 0.0, y: 0.0, z: -1.0 };

    /// Normalizes (x, y, z); `None` if its length is (near) zero.
    pub fn try_new(x: f64, y: f64, z: f64) -> Option<Self> {
        let length = (x * x + y * y + z * z).sqrt();
        (length > crate::tolerance::RESOLUTION).then(|| Self { x: x / length, y: y / length, z: z / length })
    }

    /// Normalizes (x, y, z); panics if its length is (near) zero.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self::try_new(x, y, z).expect("direction from a zero-length vector")
    }

    /// Normalizes `v`; panics if it is (near) zero.
    pub fn from_vector(v: Vector) -> Self {
        Self::new(v.x, v.y, v.z)
    }

    /// Normalizes `v`; `None` if it is (near) zero.
    pub fn try_from_vector(v: Vector) -> Option<Self> {
        Self::try_new(v.x, v.y, v.z)
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn dot(&self, other: Direction) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// The (not normalized) cross product.
    pub fn cross(&self, other: Direction) -> Vector {
        Vector::from(*self).cross(Vector::from(other))
    }

    /// Angle in [0, PI].
    pub fn angle(&self, other: Direction) -> f64 {
        Vector::from(*self).angle(Vector::from(other))
    }

    /// Signed angle in [-PI, PI], positive when `self x other` points the same way as
    /// `reference`.
    pub fn signed_angle(&self, other: Direction, reference: Direction) -> f64 {
        let angle = self.angle(other);
        if self.cross(other).dot(Vector::from(reference)) < 0.0 { -angle } else { angle }
    }

    /// True if the angle between the directions is at most `angular_tolerance`.
    pub fn approx_eq(&self, other: Direction, angular_tolerance: f64) -> bool {
        self.angle(other) <= angular_tolerance
    }

    pub fn is_parallel(&self, other: Direction, angular_tolerance: f64) -> bool {
        Vector::from(*self).is_parallel(Vector::from(other), angular_tolerance)
    }

    pub fn is_opposite(&self, other: Direction, angular_tolerance: f64) -> bool {
        Vector::from(*self).is_opposite(Vector::from(other), angular_tolerance)
    }

    pub fn is_perpendicular(&self, other: Direction, angular_tolerance: f64) -> bool {
        Vector::from(*self).is_perpendicular(Vector::from(other), angular_tolerance)
    }
}

impl Neg for Direction {
    type Output = Direction;
    fn neg(self) -> Direction {
        Direction { x: -self.x, y: -self.y, z: -self.z }
    }
}

/// A unit vector in 2D space.
#[doc(alias = "gp_Dir2d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction2 {
    x: f64,
    y: f64,
}

impl Direction2 {
    pub const X: Direction2 = Direction2 { x: 1.0, y: 0.0 };
    pub const Y: Direction2 = Direction2 { x: 0.0, y: 1.0 };
    pub const NEG_X: Direction2 = Direction2 { x: -1.0, y: 0.0 };
    pub const NEG_Y: Direction2 = Direction2 { x: 0.0, y: -1.0 };

    /// Normalizes (x, y); `None` if its length is (near) zero.
    pub fn try_new(x: f64, y: f64) -> Option<Self> {
        let length = (x * x + y * y).sqrt();
        (length > crate::tolerance::RESOLUTION).then(|| Self { x: x / length, y: y / length })
    }

    /// Normalizes (x, y); panics if its length is (near) zero.
    pub fn new(x: f64, y: f64) -> Self {
        Self::try_new(x, y).expect("direction from a zero-length vector")
    }

    /// Normalizes `v`; panics if it is (near) zero.
    pub fn from_vector(v: Vector2) -> Self {
        Self::new(v.x, v.y)
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn dot(&self, other: Direction2) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The scalar cross product.
    pub fn cross(&self, other: Direction2) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Signed angle in [-PI, PI] from `self` to `other`.
    pub fn angle(&self, other: Direction2) -> f64 {
        self.cross(other).atan2(self.dot(other))
    }

    pub fn approx_eq(&self, other: Direction2, angular_tolerance: f64) -> bool {
        self.angle(other).abs() <= angular_tolerance
    }

    pub fn is_parallel(&self, other: Direction2, angular_tolerance: f64) -> bool {
        Vector2::from(*self).is_parallel(Vector2::from(other), angular_tolerance)
    }

    pub fn is_opposite(&self, other: Direction2, angular_tolerance: f64) -> bool {
        Vector2::from(*self).is_opposite(Vector2::from(other), angular_tolerance)
    }

    pub fn is_perpendicular(&self, other: Direction2, angular_tolerance: f64) -> bool {
        Vector2::from(*self).is_perpendicular(Vector2::from(other), angular_tolerance)
    }
}

impl Neg for Direction2 {
    type Output = Direction2;
    fn neg(self) -> Direction2 {
        Direction2 { x: -self.x, y: -self.y }
    }
}
