use super::{Point, Point2};

/// An axis-aligned box. Sides may be infinite (an "open" side, e.g. for the box of an
/// infinite line); the empty box has `min > max`.
#[doc(alias = "Bnd_Box")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min: Point,
    pub max: Point,
}

impl Default for BoundingBox {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl BoundingBox {
    /// The box containing nothing; adding a point to it gives the box of that point.
    #[doc(alias = "SetVoid")]
    pub const EMPTY: BoundingBox = BoundingBox {
        min: Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        max: Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    };

    /// The box containing everything.
    #[doc(alias = "SetWhole")]
    pub const WHOLE: BoundingBox = BoundingBox {
        min: Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        max: Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
    };

    pub const fn new(min: Point, max: Point) -> Self {
        Self { min, max }
    }

    /// The smallest box containing `points` (empty if there are none).
    pub fn from_points(points: impl IntoIterator<Item = Point>) -> Self {
        let mut b = Self::EMPTY;
        for p in points {
            b.add_point(p);
        }
        b
    }

    #[doc(alias = "IsVoid")]
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y || self.min.z > self.max.z
    }

    /// True if no side is infinite.
    pub fn is_bounded(&self) -> bool {
        [self.min.x, self.min.y, self.min.z, self.max.x, self.max.y, self.max.z]
            .iter()
            .all(|c| c.is_finite())
    }

    #[doc(alias = "Add")]
    pub fn add_point(&mut self, p: Point) {
        self.min = Point::new(self.min.x.min(p.x), self.min.y.min(p.y), self.min.z.min(p.z));
        self.max = Point::new(self.max.x.max(p.x), self.max.y.max(p.y), self.max.z.max(p.z));
    }

    /// Grows the box to contain `other`.
    #[doc(alias = "Add")]
    pub fn add_box(&mut self, other: &BoundingBox) {
        if !other.is_empty() {
            self.add_point(other.min);
            self.add_point(other.max);
        }
    }

    /// The box grown by `gap` on every side (an empty box stays empty).
    #[doc(alias = "Enlarge")]
    pub fn enlarged(&self, gap: f64) -> BoundingBox {
        if self.is_empty() {
            return *self;
        }
        Self {
            min: Point::new(self.min.x - gap, self.min.y - gap, self.min.z - gap),
            max: Point::new(self.max.x + gap, self.max.y + gap, self.max.z + gap),
        }
    }

    /// True if `p` is inside or on the box.
    #[doc(alias = "IsOut")]
    pub fn contains(&self, p: Point) -> bool {
        (self.min.x..=self.max.x).contains(&p.x)
            && (self.min.y..=self.max.y).contains(&p.y)
            && (self.min.z..=self.max.z).contains(&p.z)
    }
}

/// An axis-aligned rectangle; see [`BoundingBox`].
#[doc(alias = "Bnd_Box2d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox2 {
    pub min: Point2,
    pub max: Point2,
}

impl Default for BoundingBox2 {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl BoundingBox2 {
    pub const EMPTY: BoundingBox2 = BoundingBox2 {
        min: Point2::new(f64::INFINITY, f64::INFINITY),
        max: Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    };

    pub const fn new(min: Point2, max: Point2) -> Self {
        Self { min, max }
    }

    pub fn from_points(points: impl IntoIterator<Item = Point2>) -> Self {
        let mut b = Self::EMPTY;
        for p in points {
            b.add_point(p);
        }
        b
    }

    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y
    }

    pub fn is_bounded(&self) -> bool {
        [self.min.x, self.min.y, self.max.x, self.max.y].iter().all(|c| c.is_finite())
    }

    pub fn add_point(&mut self, p: Point2) {
        self.min = Point2::new(self.min.x.min(p.x), self.min.y.min(p.y));
        self.max = Point2::new(self.max.x.max(p.x), self.max.y.max(p.y));
    }

    pub fn add_box(&mut self, other: &BoundingBox2) {
        if !other.is_empty() {
            self.add_point(other.min);
            self.add_point(other.max);
        }
    }

    pub fn enlarged(&self, gap: f64) -> BoundingBox2 {
        if self.is_empty() {
            return *self;
        }
        Self {
            min: Point2::new(self.min.x - gap, self.min.y - gap),
            max: Point2::new(self.max.x + gap, self.max.y + gap),
        }
    }

    pub fn contains(&self, p: Point2) -> bool {
        (self.min.x..=self.max.x).contains(&p.x) && (self.min.y..=self.max.y).contains(&p.y)
    }
}

// --- measures, relations and distances ----------------------------------------------------

use super::{AffineTransform, AffineTransformable, Transform, Transformable, Vector, Vector2};

impl BoundingBox {
    /// The extent along each axis (zero for an empty box).
    pub fn size(&self) -> Vector {
        if self.is_empty() {
            return Vector::ZERO;
        }
        self.max - self.min
    }

    pub fn center(&self) -> Point {
        Point::new(0.5 * (self.min.x + self.max.x), 0.5 * (self.min.y + self.max.y), 0.5 * (self.min.z + self.max.z))
    }

    /// The surface area of the box (zero for an empty or flat one).
    #[doc(alias = "BVH_Box::Area")]
    pub fn area(&self) -> f64 {
        let s = self.size();
        2.0 * (s.x * s.y + s.y * s.z + s.x * s.z)
    }

    /// True if the boxes share a point (touching counts); an empty box overlaps nothing.
    #[doc(alias = "IsOut")]
    pub fn overlaps(&self, other: &BoundingBox) -> bool {
        !self.is_empty()
            && !other.is_empty()
            && self.min.x <= other.max.x
            && other.min.x <= self.max.x
            && self.min.y <= other.max.y
            && other.min.y <= self.max.y
            && self.min.z <= other.max.z
            && other.min.z <= self.max.z
    }

    /// True if `other` lies within the box (on its boundary included); never for empty
    /// boxes.
    #[doc(alias = "BVH_Box::Contains")]
    pub fn contains_box(&self, other: &BoundingBox) -> bool {
        self.overlaps(other) && self.contains(other.min) && self.contains(other.max)
    }

    /// The point of the box closest to `p` (`p` itself inside).
    #[doc(alias = "BVH_Tools::PointBoxProjection")]
    pub fn closest_point(&self, p: Point) -> Point {
        Point::new(p.x.clamp(self.min.x, self.max.x), p.y.clamp(self.min.y, self.max.y), p.z.clamp(self.min.z, self.max.z))
    }

    /// The distance from `p` to the box (zero inside).
    #[doc(alias = "BVH_Tools::PointBoxSquareDistance")]
    pub fn distance_to_point(&self, p: Point) -> f64 {
        self.closest_point(p).distance(p)
    }

    /// The distance between the boxes (zero if they overlap).
    #[doc(alias = "BVH_Tools::BoxBoxSquareDistance")]
    pub fn distance_to_box(&self, other: &BoundingBox) -> f64 {
        let gap = |lo1: f64, hi1: f64, lo2: f64, hi2: f64| (lo2 - hi1).max(lo1 - hi2).max(0.0);
        Vector::new(
            gap(self.min.x, self.max.x, other.min.x, other.max.x),
            gap(self.min.y, self.max.y, other.min.y, other.max.y),
            gap(self.min.z, self.max.z, other.min.z, other.max.z),
        )
        .length()
    }

    /// Where the line `origin + t·direction` enters and leaves the box, `(t_enter, t_leave)`
    /// (`t_enter` negative when `origin` is inside); `None` if it misses the box or the box is
    /// behind the origin.
    #[doc(alias = "BVH_Tools::RayBoxIntersection")]
    pub fn ray_interval(&self, origin: Point, direction: Vector) -> Option<(f64, f64)> {
        if self.is_empty() {
            return None;
        }
        let (mut enter, mut leave) = (f64::NEG_INFINITY, f64::INFINITY);
        for (o, d, lo, hi) in [
            (origin.x, direction.x, self.min.x, self.max.x),
            (origin.y, direction.y, self.min.y, self.max.y),
            (origin.z, direction.z, self.min.z, self.max.z),
        ] {
            if d == 0.0 {
                if o < lo || o > hi {
                    return None;
                }
                continue;
            }
            let (a, b) = ((lo - o) / d, (hi - o) / d);
            enter = enter.max(a.min(b));
            leave = leave.min(a.max(b));
        }
        (enter <= leave && leave >= 0.0).then_some((enter, leave))
    }
}

impl Transformable for BoundingBox {
    /// The box of the transformed box (of its eight corners).
    fn transformed(&self, t: &Transform) -> Self {
        if self.is_empty() {
            return *self;
        }
        BoundingBox::from_points(corners(self).into_iter().map(|p| p.transformed(t)))
    }
}

impl AffineTransformable for BoundingBox {
    /// The box of the transformed box (of its eight corners).
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        if self.is_empty() {
            return *self;
        }
        BoundingBox::from_points(corners(self).into_iter().map(|p| t.transform_point(p)))
    }
}

fn corners(b: &BoundingBox) -> [Point; 8] {
    std::array::from_fn(|i| {
        Point::new(
            if i & 1 == 0 { b.min.x } else { b.max.x },
            if i & 2 == 0 { b.min.y } else { b.max.y },
            if i & 4 == 0 { b.min.z } else { b.max.z },
        )
    })
}

impl BoundingBox2 {
    /// The extent along each axis (zero for an empty rectangle).
    pub fn size(&self) -> Vector2 {
        if self.is_empty() {
            return Vector2::ZERO;
        }
        self.max - self.min
    }

    pub fn center(&self) -> Point2 {
        Point2::new(0.5 * (self.min.x + self.max.x), 0.5 * (self.min.y + self.max.y))
    }

    pub fn area(&self) -> f64 {
        let s = self.size();
        s.x * s.y
    }

    /// True if the rectangles share a point (touching counts).
    pub fn overlaps(&self, other: &BoundingBox2) -> bool {
        !self.is_empty()
            && !other.is_empty()
            && self.min.x <= other.max.x
            && other.min.x <= self.max.x
            && self.min.y <= other.max.y
            && other.min.y <= self.max.y
    }

    /// True if `other` lies within the rectangle (on its boundary included).
    pub fn contains_box(&self, other: &BoundingBox2) -> bool {
        self.overlaps(other) && self.contains(other.min) && self.contains(other.max)
    }

    pub fn closest_point(&self, p: Point2) -> Point2 {
        Point2::new(p.x.clamp(self.min.x, self.max.x), p.y.clamp(self.min.y, self.max.y))
    }

    pub fn distance_to_point(&self, p: Point2) -> f64 {
        self.closest_point(p).distance(p)
    }

    pub fn distance_to_box(&self, other: &BoundingBox2) -> f64 {
        let gap = |lo1: f64, hi1: f64, lo2: f64, hi2: f64| (lo2 - hi1).max(lo1 - hi2).max(0.0);
        gap(self.min.x, self.max.x, other.min.x, other.max.x).hypot(gap(self.min.y, self.max.y, other.min.y, other.max.y))
    }
}

impl BoundingBox2 {
    /// Where the line `origin + t·direction` enters and leaves the rectangle; see
    /// [`BoundingBox::ray_interval`].
    pub fn ray_interval(&self, origin: Point2, direction: Vector2) -> Option<(f64, f64)> {
        if self.is_empty() {
            return None;
        }
        let (mut enter, mut leave) = (f64::NEG_INFINITY, f64::INFINITY);
        for (o, d, lo, hi) in [(origin.x, direction.x, self.min.x, self.max.x), (origin.y, direction.y, self.min.y, self.max.y)] {
            if d == 0.0 {
                if o < lo || o > hi {
                    return None;
                }
                continue;
            }
            let (a, b) = ((lo - o) / d, (hi - o) / d);
            enter = enter.max(a.min(b));
            leave = leave.min(a.max(b));
        }
        (enter <= leave && leave >= 0.0).then_some((enter, leave))
    }
}
