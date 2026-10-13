use super::{Direction, Direction2, Point, Point2};

/// An oriented line used for positioning: rotation and symmetry axes, the axis of a
/// cylinder or a revolution, ...
#[doc(alias = "gp_Ax1")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis {
    pub origin: Point,
    pub direction: Direction,
}

impl Axis {
    pub const X: Axis = Axis::new(Point::ORIGIN, Direction::X);
    pub const Y: Axis = Axis::new(Point::ORIGIN, Direction::Y);
    pub const Z: Axis = Axis::new(Point::ORIGIN, Direction::Z);

    pub const fn new(origin: Point, direction: Direction) -> Self {
        Self { origin, direction }
    }
}

/// An oriented line in 2D.
#[doc(alias = "gp_Ax2d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis2 {
    pub origin: Point2,
    pub direction: Direction2,
}

impl Axis2 {
    pub const X: Axis2 = Axis2::new(Point2::ORIGIN, Direction2::X);
    pub const Y: Axis2 = Axis2::new(Point2::ORIGIN, Direction2::Y);

    pub const fn new(origin: Point2, direction: Direction2) -> Self {
        Self { origin, direction }
    }
}
