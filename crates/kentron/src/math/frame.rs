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

    /// The frame with main axis `z` at `origin` and the conventional default X axis.
    ///
    /// The X axis is found by dropping the component of `z` with the smallest magnitude (y
    /// first, then x, then z when magnitudes tie) and swapping the other two with one sign
    /// change, which gives a vector perpendicular to `z`. For the world axes this yields
    /// X = +x for z = +z, X = +z for z = +x and X = +z for z = +y. Circles, conics and
    /// elementary surfaces built from an axis alone take their angular origin from it.
    pub fn from_z(origin: Point, z: Direction) -> Self {
        let (p, q, r) = (z.x(), z.y(), z.z());
        let (a, b, c) = (p.abs(), q.abs(), r.abs());
        let hint = if b <= a && b <= c {
            // y is the smallest component: X lies in the world XZ plane.
            if a > c { Direction::new(-r, 0.0, p) } else { Direction::new(r, 0.0, -p) }
        } else if a <= b && a <= c {
            // x is the smallest component: X lies in the world YZ plane.
            if b > c { Direction::new(0.0, -r, q) } else { Direction::new(0.0, r, -q) }
        } else if a > b {
            // z is the smallest component: X lies in the world XY plane.
            Direction::new(-q, p, 0.0)
        } else {
            Direction::new(q, -p, 0.0)
        };
        Self::new(origin, z, hint)
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

    /// Turns the frame so its main axis becomes `z`, keeping the origin and the handedness.
    ///
    /// Generally the X axis becomes the old X made perpendicular to `z` and Y completes the
    /// frame. When `z` is (nearly) along the old X axis that cannot work, so the axes are
    /// rotated instead: along +X the old Y and Z become X and Y, along -X the old Z becomes
    /// X and Y stays. In those cases `z` is stored as given.
    pub fn set_z(&mut self, z: Direction) {
        let right_handed = self.is_right_handed();
        let t = z.dot(self.x);
        if 1.0 - t.abs() <= crate::tolerance::ANGULAR {
            if t > 0.0 {
                self.x = self.y;
                self.y = self.z;
            } else {
                self.x = self.z;
            }
            self.z = z;
        } else {
            self.z = z;
            self.x = Direction::from_vector(Vector::from(z).cross(self.x.cross(z)));
            self.y = Direction::from_vector(if right_handed { z.cross(self.x) } else { self.x.cross(z) });
        }
    }

    /// Moves the frame's origin to the axis' origin and turns its main axis to the axis'
    /// direction as [`Frame::set_z`] does.
    pub fn set_z_axis(&mut self, axis: Axis) {
        self.origin = axis.origin;
        self.set_z(axis.direction);
    }

    /// Turns the frame about its main axis so its X axis points along `x` (made
    /// perpendicular to Z), keeping the origin and the handedness. When `x` is (nearly)
    /// along Z the old X becomes the main axis instead and `x` is stored as given; along +Z
    /// the Y axis is reversed so the handedness survives.
    pub fn set_x(&mut self, x: Direction) {
        let right_handed = self.is_right_handed();
        let t = x.dot(self.z);
        if 1.0 - t.abs() <= crate::tolerance::ANGULAR {
            if t > 0.0 {
                self.y = -self.y;
            }
            self.z = self.x;
            self.x = x;
        } else {
            self.x = Direction::from_vector(Vector::from(self.z).cross(x.cross(self.z)));
            self.y = Direction::from_vector(if right_handed { self.z.cross(self.x) } else { self.x.cross(self.z) });
        }
    }

    /// Turns the frame about its main axis so its Y axis points along `y` (made
    /// perpendicular to Z), keeping the origin and the handedness. When `y` is (nearly)
    /// along Z the old Y becomes the main axis instead and `y` is stored as given; along +Z
    /// the X axis is reversed so the handedness survives.
    pub fn set_y(&mut self, y: Direction) {
        let right_handed = self.is_right_handed();
        let t = y.dot(self.z);
        if 1.0 - t.abs() <= crate::tolerance::ANGULAR {
            if t > 0.0 {
                self.x = -self.x;
            }
            self.z = self.y;
            self.y = y;
        } else {
            let w = Direction::from_vector(y.cross(self.z));
            self.y = Direction::from_vector(Vector::from(self.z).cross(Vector::from(w)));
            self.x = if right_handed { w } else { -w };
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

#[cfg(test)]
mod tests {
    use super::*;

    fn xyz(d: Direction) -> (f64, f64, f64) {
        // (Signed zeros compare equal.)
        (d.x(), d.y(), d.z())
    }

    #[test]
    fn default_x_axes_of_the_world_axes() {
        let cases = [
            ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            ((0.0, 0.0, -1.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, -1.0, 0.0)),
            ((-1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0)),
            ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            ((0.0, -1.0, 0.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0)),
        ];
        for (z, x, y) in cases {
            let f = Frame::from_z(Point::ORIGIN, Direction::new(z.0, z.1, z.2));
            assert_eq!(xyz(f.x()), x, "X for z = {z:?}");
            assert_eq!(xyz(f.y()), y, "Y for z = {z:?}");
        }
    }

    #[test]
    fn default_x_axes_of_oblique_directions() {
        // (z, the expected X before normalization)
        let cases = [
            // y smallest (ties with x or z go to y)
            ((1.0, 0.0, 1.0), (1.0, 0.0, -1.0)),
            ((1.0, 1.0, 1.0), (1.0, 0.0, -1.0)),
            ((3.0, 1.0, 2.0), (-2.0, 0.0, 3.0)),
            ((-3.0, 1.0, 2.0), (-2.0, 0.0, -3.0)),
            ((2.0, -1.0, 3.0), (3.0, 0.0, -2.0)),
            // x smallest
            ((0.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
            ((1.0, 3.0, 2.0), (0.0, -2.0, 3.0)),
            ((1.0, -2.0, -3.0), (0.0, -3.0, 2.0)),
            // z smallest (x and y tie: the second sign)
            ((1.0, 1.0, 0.0), (1.0, -1.0, 0.0)),
            ((-1.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
            ((3.0, 2.0, 1.0), (-2.0, 3.0, 0.0)),
            ((2.0, -3.0, 1.0), (-3.0, -2.0, 0.0)),
            ((-2.0, 3.0, -1.0), (3.0, 2.0, 0.0)),
        ];
        for (z, x) in cases {
            let zd = Direction::new(z.0, z.1, z.2);
            let f = Frame::from_z(Point::ORIGIN, zd);
            let expected = Direction::new(x.0, x.1, x.2);
            assert!(f.x().approx_eq(expected, 1e-12), "X for z = {z:?}: {:?}", f.x());
            assert!(f.is_right_handed());
            assert!(f.x().dot(zd).abs() < 1e-15);
        }
    }

    #[test]
    fn default_frames_are_orthonormal_and_right_handed() {
        // A deterministic spread of directions over the sphere.
        for i in 0..500 {
            let t = i as f64 * 0.618_033_988_749_895;
            let z = 1.0 - 2.0 * ((i as f64 + 0.5) / 500.0);
            let r = (1.0 - z * z).sqrt();
            let d = Direction::new(r * (t * std::f64::consts::TAU).cos(), r * (t * std::f64::consts::TAU).sin(), z);
            let f = Frame::from_z(Point::new(1.0, 2.0, 3.0), d);
            assert!(f.x().dot(d).abs() < 1e-15);
            assert!(f.y().dot(d).abs() < 1e-15);
            assert!(f.x().dot(f.y()).abs() < 1e-15);
            assert!(f.is_right_handed());
            assert_eq!(f.z(), d);
        }
    }

    fn assert_orthonormal(f: &Frame) {
        for (a, b) in [(f.x(), f.y()), (f.y(), f.z()), (f.z(), f.x())] {
            assert!(a.dot(b).abs() < 1e-12, "{f:?}");
        }
    }

    fn frames() -> [Frame; 2] {
        let right = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(1.0, 2.0, 2.0), Direction::new(0.0, 1.0, -1.0));
        let mut left = right;
        left.flip_y();
        [right, left]
    }

    #[test]
    fn set_z_keeps_handedness_in_every_branch() {
        for f in frames() {
            let handed = f.is_right_handed();
            for (d, generic) in [(f.x(), false), (-f.x(), false), (Direction::new(1.0, -1.0, 0.5), true)] {
                let mut g = f;
                g.set_z(d);
                assert_eq!(g.is_right_handed(), handed);
                assert_eq!(g.z(), d);
                assert_eq!(g.origin(), f.origin());
                if generic {
                    assert_orthonormal(&g);
                }
            }
            let mut g = f;
            g.set_z(f.x());
            assert_eq!((g.x(), g.y()), (f.y(), f.z()));
            let mut g = f;
            g.set_z(-f.x());
            assert_eq!((g.x(), g.y()), (f.z(), f.y()));
            let mut g = f;
            let axis = Axis::new(Point::new(-1.0, 0.0, 5.0), Direction::new(0.0, 0.0, 1.0));
            g.set_z_axis(axis);
            assert_eq!((g.origin(), g.z()), (axis.origin, axis.direction));
            assert_eq!(g.is_right_handed(), handed);
        }
    }

    #[test]
    fn set_x_and_set_y_keep_handedness_in_every_branch() {
        for f in frames() {
            let handed = f.is_right_handed();
            let generic = Direction::new(0.3, -1.0, 0.5);
            for d in [f.z(), -f.z(), generic] {
                let mut g = f;
                g.set_x(d);
                assert_eq!(g.is_right_handed(), handed);
                let mut h = f;
                h.set_y(d);
                assert_eq!(h.is_right_handed(), handed);
                if d == generic {
                    assert_orthonormal(&g);
                    assert_orthonormal(&h);
                    assert_eq!(g.z(), f.z());
                    assert_eq!(h.z(), f.z());
                    let along = |a: Direction| {
                        let v = Vector::from(d) - Vector::from(f.z()) * d.dot(f.z());
                        a.approx_eq(Direction::from_vector(v), 1e-12)
                    };
                    assert!(along(g.x()));
                    assert!(along(h.y()));
                } else {
                    assert_eq!((g.x(), g.z()), (d, f.x()));
                    assert_eq!((h.y(), h.z()), (d, f.y()));
                }
            }
        }
    }
}
