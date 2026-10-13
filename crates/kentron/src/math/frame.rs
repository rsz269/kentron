use super::{Axis, Axis2, Direction, Direction2, Point, Point2, Transform2, Transformable2, Vector};

/// A local coordinate system: an origin and three orthonormal axes. Usually right-handed
/// (`z = x × y`); flipping one axis makes it left-handed, which for a surface positioned by
/// the frame reverses its normal.
#[doc(alias = "gp_Ax2")]
#[doc(alias = "gp_Ax3")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    origin: Point,
    x: Direction,
    y: Direction,
    z: Direction,
}

impl Default for Frame {
    fn default() -> Self {
        Self::WORLD
    }
}

impl Frame {
    /// The global coordinate system.
    pub const WORLD: Frame =
        Frame { origin: Point::ORIGIN, x: Direction::X, y: Direction::Y, z: Direction::Z };

    /// Right-handed frame with main axis `z` and X axis from `x_hint` (projected to be
    /// perpendicular to `z`). Panics if `x_hint` is parallel to `z`.
    pub fn new(origin: Point, z: Direction, x_hint: Direction) -> Self {
        let zv = Vector::from(z);
        let hint = Vector::from(x_hint);
        let x = Direction::try_from_vector(hint - zv * hint.dot(zv))
            .expect("frame X hint parallel to its Z axis");
        let y = Direction::from_vector(zv.cross(Vector::from(x)));
        Self { origin, x, y, z }
    }

    /// A frame from axes known to be orthonormal (e.g. transformed ones).
    pub(crate) fn from_axes_unchecked(origin: Point, x: Direction, y: Direction, z: Direction) -> Self {
        Self { origin, x, y, z }
    }

    /// Right-handed frame with main axis `z`; the X axis is chosen automatically.
    pub fn from_z(origin: Point, z: Direction) -> Self {
        // OCCT's gp_Ax2 choice, so that parametrizations match OCCT's: zero the smallest
        // component of `z` and swap the other two.
        let (a, b, c) = (z.x().abs(), z.y().abs(), z.z().abs());
        let x_hint = if b <= a && b <= c {
            if a > c { Direction::new(-z.z(), 0.0, z.x()) } else { Direction::new(z.z(), 0.0, -z.x()) }
        } else if a <= b && a <= c {
            if b > c { Direction::new(0.0, -z.z(), z.y()) } else { Direction::new(0.0, z.z(), -z.y()) }
        } else if a > b {
            Direction::new(-z.y(), z.x(), 0.0)
        } else {
            Direction::new(z.y(), -z.x(), 0.0)
        };
        Self::new(origin, z, x_hint)
    }

    pub fn origin(&self) -> Point {
        self.origin
    }

    pub fn x(&self) -> Direction {
        self.x
    }

    pub fn y(&self) -> Direction {
        self.y
    }

    /// The main direction (normal of the frame's XY plane).
    pub fn z(&self) -> Direction {
        self.z
    }

    /// The main axis: origin and Z direction.
    pub fn z_axis(&self) -> Axis {
        Axis::new(self.origin, self.z)
    }

    pub fn is_right_handed(&self) -> bool {
        self.x.cross(self.y).dot(Vector::from(self.z)) > 0.0
    }

    pub fn flip_x(&mut self) {
        self.x = -self.x;
    }

    pub fn flip_y(&mut self) {
        self.y = -self.y;
    }

    pub fn flip_z(&mut self) {
        self.z = -self.z;
    }

    /// Re-orients the frame so that its main direction is `z`, keeping its handedness.
    pub fn set_z(&mut self, z: Direction) {
        // As OCCT's gp_Ax3::SetDirection.
        let dot = Vector::from(z).dot(Vector::from(self.x));
        if 1.0 - dot.abs() <= crate::tolerance::ANGULAR {
            if dot > 0.0 {
                (self.x, self.y) = (self.y, self.z);
            } else {
                self.x = self.z;
            }
            self.z = z;
        } else {
            let right_handed = self.is_right_handed();
            let zv = Vector::from(z);
            self.z = z;
            self.x = Direction::from_vector(zv.cross(Vector::from(self.x).cross(zv)));
            self.y = Direction::from_vector(if right_handed { zv.cross(Vector::from(self.x)) } else { Vector::from(self.x).cross(zv) });
        }
    }

    /// Moves the origin to `axis.origin` and sets the main direction, keeping handedness.
    pub fn set_z_axis(&mut self, axis: Axis) {
        self.origin = axis.origin;
        self.set_z(axis.direction);
    }

    /// Sets the X direction (projected to be perpendicular to Z), keeping handedness.
    pub fn set_x(&mut self, x: Direction) {
        // As OCCT's gp_Ax3::SetXDirection.
        let dot = Vector::from(x).dot(Vector::from(self.z));
        if 1.0 - dot.abs() <= crate::tolerance::ANGULAR {
            if dot > 0.0 {
                self.z = self.x;
                self.y = -self.y;
            } else {
                self.z = self.x;
            }
            self.x = x;
        } else {
            let right_handed = self.is_right_handed();
            let zv = Vector::from(self.z);
            self.x = Direction::from_vector(zv.cross(Vector::from(x).cross(zv)));
            self.y = Direction::from_vector(if right_handed { zv.cross(Vector::from(self.x)) } else { Vector::from(self.x).cross(zv) });
        }
    }

    /// Sets the Y direction (projected to be perpendicular to Z), keeping handedness.
    pub fn set_y(&mut self, y: Direction) {
        // As OCCT's gp_Ax3::SetYDirection.
        let dot = Vector::from(y).dot(Vector::from(self.z));
        if 1.0 - dot.abs() <= crate::tolerance::ANGULAR {
            if dot > 0.0 {
                self.z = self.y;
                self.x = -self.x;
            } else {
                self.z = self.y;
            }
            self.y = y;
        } else {
            let right_handed = self.is_right_handed();
            let zv = Vector::from(self.z);
            let x = Direction::from_vector(Vector::from(y).cross(zv));
            self.y = Direction::from_vector(zv.cross(Vector::from(x)));
            self.x = if right_handed { x } else { -x };
        }
    }
}

/// A 2D coordinate system: an origin and two orthonormal axes, counter-clockwise
/// (`y` is `x` turned by +90°) or clockwise. Conics in 2D are placed by one; a clockwise
/// frame reverses their direction of travel.
#[doc(alias = "gp_Ax22d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame2 {
    origin: Point2,
    x: Direction2,
    y: Direction2,
}

impl Default for Frame2 {
    fn default() -> Self {
        Self::WORLD
    }
}

impl Frame2 {
    /// The global 2D coordinate system.
    pub const WORLD: Frame2 = Frame2 { origin: Point2::ORIGIN, x: Direction2::X, y: Direction2::Y };

    /// The frame at `origin` with X axis `x`, counter-clockwise if `ccw`.
    pub fn new(origin: Point2, x: Direction2, ccw: bool) -> Self {
        let turned = Direction2::new(-x.y(), x.x());
        Self { origin, x, y: if ccw { turned } else { -turned } }
    }

    /// The frame with X axis `x_axis`, counter-clockwise if `ccw`.
    pub fn from_axis(x_axis: Axis2, ccw: bool) -> Self {
        Self::new(x_axis.origin, x_axis.direction, ccw)
    }

    /// The frame with the given axes; `None` unless they are perpendicular.
    pub fn from_axes(origin: Point2, x: Direction2, y: Direction2) -> Option<Self> {
        x.is_perpendicular(y, crate::tolerance::ANGULAR).then_some(Self { origin, x, y })
    }

    pub fn origin(&self) -> Point2 {
        self.origin
    }

    pub fn x(&self) -> Direction2 {
        self.x
    }

    pub fn y(&self) -> Direction2 {
        self.y
    }

    pub fn x_axis(&self) -> Axis2 {
        Axis2::new(self.origin, self.x)
    }

    pub fn y_axis(&self) -> Axis2 {
        Axis2::new(self.origin, self.y)
    }

    /// True if `y` is `x` turned counter-clockwise.
    pub fn is_ccw(&self) -> bool {
        self.x.cross(self.y) > 0.0
    }
}

impl Transformable2 for Frame2 {
    /// The frame's origin and axes moved by `t` (a reflection makes a counter-clockwise
    /// frame clockwise).
    fn transformed(&self, t: &Transform2) -> Self {
        Frame2 { origin: self.origin.transformed(t), x: self.x.transformed(t), y: self.y.transformed(t) }
    }
}
