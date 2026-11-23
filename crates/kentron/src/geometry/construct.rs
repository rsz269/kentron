//! Constructing curves and surfaces from points, distances and other geometry.
//!
//! The constructors live on the constructed types; they report invalid input through
//! [`ConstructionError`] instead of panicking.

use std::fmt;

use crate::math::{Axis, Axis2, Direction, Direction2, Frame, Frame2, Point, Point2, Vector};

use super::curve::{Circle, Ellipse, Hyperbola, TrimmedCurve};
use super::curve2::{Circle2, Hyperbola2, Line2, Parabola2, TrimmedCurve2};
use super::surface::{Cone, Cylinder, Plane};

/// Why a construction failed.
#[doc(alias = "gce_ErrorType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConstructionError {
    /// Points that must be distinct coincide.
    ConfusedPoints,
    /// Points that must span a plane are collinear.
    ColinearPoints,
    /// A radius, or the radius implied by the input, is negative.
    NegativeRadius,
    /// An axis has (near) zero length.
    NullAxis,
    /// The major radius would be smaller than the minor one.
    InvertAxis,
    /// An angle is out of range (e.g. a cone half-angle of zero or ±π/2).
    BadAngle,
}

impl fmt::Display for ConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::ConfusedPoints => "coincident points",
            Self::ColinearPoints => "collinear points",
            Self::NegativeRadius => "negative radius",
            Self::NullAxis => "null axis",
            Self::InvertAxis => "major radius smaller than minor radius",
            Self::BadAngle => "angle out of range",
        })
    }
}

impl std::error::Error for ConstructionError {}

pub type Result<T> = std::result::Result<T, ConstructionError>;

fn check_radius(radius: f64) -> Result<()> {
    if radius < 0.0 { Err(ConstructionError::NegativeRadius) } else { Ok(()) }
}

impl Circle {
    /// The arc from `p1` through `p2` to `p3`.
    #[doc(alias = "GC_MakeArcOfCircle")]
    pub fn arc_through(p1: Point, p2: Point, p3: Point) -> Result<TrimmedCurve> {
        let (b, c) = (p2 - p1, p3 - p1);
        if b.length() <= crate::tolerance::RESOLUTION || c.length() <= crate::tolerance::RESOLUTION || (p3 - p2).length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let normal = b.cross(c);
        if normal.length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        // The circumcenter, and the circle starting at p1 round the normal (so that p2 comes
        // before p3).
        let center = p1 + (c.cross(normal) * b.dot(b) + normal.cross(b) * c.dot(c)) * (0.5 / normal.dot(normal));
        let frame = Frame::new(center, Direction::from_vector(normal), Direction::from_vector(p1 - center));
        let radius = (center.distance(p1) + center.distance(p2) + center.distance(p3)) / 3.0;
        let circle = Circle::new(frame, radius);
        let end = circle_parameter(&circle, p3);
        Ok(TrimmedCurve::new(super::curve::Curve::Circle(circle), 0.0, end))
    }

    /// The arc starting at `p1` with tangent `tangent` there and ending at `p2`.
    #[doc(alias = "GC_MakeArcOfCircle")]
    pub fn arc_tangent(p1: Point, tangent: Vector, p2: Point) -> Result<TrimmedCurve> {
        let chord = p2 - p1;
        if chord.length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let dir = chord * (1.0 / chord.length());
        let tangent = Direction::try_from_vector(tangent).ok_or(ConstructionError::NullAxis)?;
        let axis = Vector::from(tangent).cross(dir);
        if axis.length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        let axis = axis * (1.0 / axis.length());
        // The center: on the chord's bisector and on the normal to the tangent at p1.
        let to_center = Vector::from(tangent).cross(axis);
        let bisector = axis.cross(dir);
        let middle = p1 + chord * 0.5;
        // p1 + s·to_center = middle + t·bisector, in the plane of the arc.
        let (a, b, c) = (to_center.dot(to_center), to_center.dot(bisector), bisector.dot(bisector));
        let w = p1 - middle;
        let (d, e) = (to_center.dot(w), bisector.dot(w));
        let det = a * c - b * b;
        if det.abs() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        let s = (b * e - c * d) / det;
        let center = p1 + to_center * s;
        let radius = 0.5 * (center.distance(p1) + center.distance(p2));
        let frame = Frame::new(center, Direction::from_vector(axis), Direction::from_vector(to_center));
        let circle = Circle::new(frame, radius);
        let (start, end) = (circle_parameter(&circle, p1), circle_parameter(&circle, p2));
        Ok(TrimmedCurve::new(super::curve::Curve::Circle(circle), start, end))
    }
}

impl Circle2 {
    /// The circle through three points. Points coinciding within
    /// [`RESOLUTION`](crate::tolerance::RESOLUTION) give a zero-radius circle.
    #[doc(alias = "gce_MakeCirc2d")]
    pub fn through_points(p1: Point2, p2: Point2, p3: Point2) -> Result<Circle2> {
        let (d12, d13, d23) = (p1.distance(p2), p1.distance(p3), p2.distance(p3));
        let res = crate::tolerance::RESOLUTION;
        if d12 < res && d13 < res && d23 < res {
            return Ok(Circle2::new(Frame2::new(p1, Direction2::X, true), 0.0));
        }
        // Two coincide: the circle on the other two as a diameter.
        let pair = if d12 < res { Some((p1, p3)) } else if d13 < res || d23 < res { Some((p1, p2)) } else { None };
        if let Some((a, b)) = pair {
            let center = Point2::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
            let x = Direction2::try_new(a.x - center.x, a.y - center.y).unwrap_or(Direction2::X);
            return Ok(Circle2::new(Frame2::new(center, x, true), 0.5 * a.distance(b)));
        }
        // The circumcircle, counter-clockwise through the points in order when they turn
        // that way.
        let (bx, by, cx, cy) = (p2.x - p1.x, p2.y - p1.y, p3.x - p1.x, p3.y - p1.y);
        let d = 2.0 * (bx * cy - by * cx);
        if d == 0.0 {
            return Err(ConstructionError::ColinearPoints);
        }
        let (b2, c2) = (bx * bx + by * by, cx * cx + cy * cy);
        let center = Point2::new(p1.x + (cy * b2 - by * c2) / d, p1.y + (bx * c2 - cx * b2) / d);
        let radius = (center.distance(p1) + center.distance(p2) + center.distance(p3)) / 3.0;
        let x = Direction2::try_new(p1.x - center.x, p1.y - center.y).unwrap_or(Direction2::X);
        Ok(Circle2::new(Frame2::new(center, x, d > 0.0), radius))
    }

    /// The concentric circle through `p` (zero radius if `p` is the center).
    #[doc(alias = "GC_MakeCircle2d")]
    pub fn concentric_through(&self, p: Point2) -> Circle2 {
        Circle2::new(*self.frame(), self.center().distance(p))
    }

    /// The concentric circle whose radius differs by `distance`.
    #[doc(alias = "GC_MakeCircle2d")]
    pub fn concentric_offset(&self, distance: f64) -> Result<Circle2> {
        let radius = self.radius() + distance;
        check_radius(radius)?;
        Ok(Circle2::new(*self.frame(), radius))
    }
}

impl Ellipse {
    /// The ellipse centered at `center` with its major apex at `s1` and passing through
    /// `s2` (which fixes the plane and the minor radius).
    #[doc(alias = "gce_MakeElips")]
    pub fn from_points(s1: Point, s2: Point, center: Point) -> Result<Ellipse> {
        let d1 = s1.distance(center);
        if d1 < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::NullAxis);
        }
        let x = Direction::from_vector(s1 - center);
        let d2 = super::curve::Line::new(center, x).distance_to_point(s2);
        if d1 < d2 || d2 < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::InvertAxis);
        }
        let normal = Vector::from(x).cross(s2 - center);
        if normal.length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::InvertAxis);
        }
        Ok(Ellipse::new(Frame::new(center, Direction::from_vector(normal), x), d1, d2))
    }
}

impl Hyperbola {
    /// Like [`new`](Self::new), without panicking. The major radius may be smaller than
    /// the minor one.
    #[doc(alias = "gce_MakeHypr")]
    pub fn try_new(frame: Frame, major_radius: f64, minor_radius: f64) -> Result<Hyperbola> {
        check_radius(major_radius)?;
        check_radius(minor_radius)?;
        Ok(Self::new(frame, major_radius, minor_radius))
    }

    /// The hyperbola centered at `center` with its apex at `s1`; `s2` fixes the plane and
    /// is at distance `minor_radius` from the major axis.
    #[doc(alias = "gce_MakeHypr")]
    pub fn from_points(s1: Point, s2: Point, center: Point) -> Result<Hyperbola> {
        let (a, b) = (s1 - center, s2 - center);
        if a.length() <= crate::tolerance::RESOLUTION || b.length() <= crate::tolerance::RESOLUTION || s1.distance(s2) <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let x = Direction::from_vector(a);
        let minor = super::curve::Line::new(center, x).distance_to_point(s2);
        if minor <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        let normal = Vector::from(x).cross(b);
        if normal.length() <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        Ok(Hyperbola::new(Frame::new(center, Direction::from_vector(normal), x), a.length(), minor))
    }
}

impl Hyperbola2 {
    /// Like [`new`](Self::new), without panicking.
    #[doc(alias = "gce_MakeHypr2d")]
    pub fn try_new(frame: Frame2, major_radius: f64, minor_radius: f64) -> Result<Hyperbola2> {
        check_radius(major_radius)?;
        check_radius(minor_radius)?;
        Ok(Self::new(frame, major_radius, minor_radius))
    }

    /// 2D version of [`Hyperbola::from_points`].
    #[doc(alias = "gce_MakeHypr2d")]
    pub fn from_points(s1: Point2, s2: Point2, center: Point2) -> Result<Hyperbola2> {
        let (ax, ay, bx, by) = (s1.x - center.x, s1.y - center.y, s2.x - center.x, s2.y - center.y);
        let (da, db) = (ax.hypot(ay), bx.hypot(by));
        if da <= crate::tolerance::RESOLUTION || db <= crate::tolerance::RESOLUTION || s1.distance(s2) <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let x = Direction2::try_new(ax, ay).ok_or(ConstructionError::ConfusedPoints)?;
        let minor = (x.x() * by - x.y() * bx).abs();
        if minor <= crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ColinearPoints);
        }
        Ok(Hyperbola2::new(Frame2::new(center, x, x.x() * by - x.y() * bx > 0.0), da, minor))
    }
}

impl Parabola2 {
    /// The parabola with the given directrix and focus, its frame counter-clockwise if
    /// `ccw`. A focus on the directrix gives a degenerate parabola (zero focal length).
    #[doc(alias = "GC_MakeParabola2d")]
    #[doc(alias = "gce_MakeParab2d")]
    pub fn from_directrix_focus(directrix: Axis2, focus: Point2, ccw: bool) -> Parabola2 {
        let (o, d) = (directrix.origin, directrix.direction);
        let along = (focus.x - o.x) * d.x() + (focus.y - o.y) * d.y();
        let foot = Point2::new(o.x + d.x() * along, o.y + d.y() * along);
        let apex = Point2::new(0.5 * (foot.x + focus.x), 0.5 * (foot.y + focus.y));
        let focal = 0.5 * foot.distance(focus);
        // The axis points from the directrix to the focus; with the focus on the directrix,
        // a quarter turn from it (clockwise if counter-clockwise).
        let x = match Direction2::try_new(focus.x - foot.x, focus.y - foot.y).filter(|_| focal > 0.0) {
            Some(x) => x,
            None if ccw => Direction2::new(d.y(), -d.x()),
            None => Direction2::new(-d.y(), d.x()),
        };
        let frame = Frame2::from_axes(apex, x, d).unwrap_or_else(|| Frame2::new(apex, x, ccw));
        Parabola2::new(frame, focal)
    }
}

impl Line2 {
    /// The segment from `p1` to `p2`.
    #[doc(alias = "GC_MakeSegment2d")]
    pub fn segment(p1: Point2, p2: Point2) -> Result<TrimmedCurve2> {
        let length = p1.distance(p2);
        if !(length > crate::tolerance::RESOLUTION) {
            return Err(ConstructionError::ConfusedPoints);
        }
        let direction = Direction2::try_new(p2.x - p1.x, p2.y - p1.y).ok_or(ConstructionError::ConfusedPoints)?;
        Ok(TrimmedCurve2::new(super::curve2::Curve2::Line(Line2::new(p1, direction)), 0.0, length))
    }

    /// The segment of the line through `p1` along `direction`, from `p1` to the projection
    /// of `p2`.
    #[doc(alias = "GC_MakeSegment2d")]
    pub fn segment_to_projection(
        p1: Point2,
        direction: Direction2,
        p2: Point2,
    ) -> Result<TrimmedCurve2> {
        let last = (p2.x - p1.x) * direction.x() + (p2.y - p1.y) * direction.y();
        if !(last.abs() > crate::tolerance::RESOLUTION) {
            return Err(ConstructionError::ConfusedPoints);
        }
        let line = super::curve2::Curve2::Line(Line2::new(p1, direction));
        Ok(TrimmedCurve2::new(line, last.min(0.0), last.max(0.0)))
    }
}

impl Plane {
    /// The plane through three points, with its origin at `a` and its X axis towards `b`, its
    /// normal along `(b - a) × (c - a)`. Fails where the points coincide
    /// ([`ConstructionError::ConfusedPoints`]) or lie on a line
    /// ([`ConstructionError::ColinearPoints`]), within the linear tolerance.
    #[doc(alias = "GC_MakePlane")]
    #[doc(alias = "bg_make_plane_3pts")]
    pub fn through_points(a: Point, b: Point, c: Point) -> Result<Plane> {
        let (ab, ac) = (b - a, c - a);
        let normal = ab.cross(ac);
        let longest = ab.length().max(ac.length()).max((c - b).length());
        // (The normal's length is twice the triangle's area: over the longest side, its
        // height.)
        if longest <= crate::tolerance::LINEAR {
            return Err(ConstructionError::ConfusedPoints);
        }
        if normal.length() / longest <= crate::tolerance::LINEAR {
            return Err(ConstructionError::ColinearPoints);
        }
        Ok(Plane::from_frame(Frame::new(a, Direction::from_vector(normal), Direction::from_vector(ab))))
    }

    /// The plane through `axis.origin` normal to the axis.
    #[doc(alias = "GC_MakePlane")]
    pub fn from_axis(axis: Axis) -> Plane {
        Plane::new(axis.origin, axis.direction)
    }

    /// The parallel plane at signed `distance` along the normal.
    #[doc(alias = "GC_MakePlane")]
    pub fn offset(&self, distance: f64) -> Plane {
        let f = self.frame();
        Plane::from_frame(Frame::new(f.origin() + Vector::from(f.z()) * distance, f.z(), f.x()))
    }

    /// The parallel plane through `p`, with `p` as its origin.
    #[doc(alias = "GC_MakePlane")]
    pub fn parallel_through(&self, p: Point) -> Plane {
        let f = self.frame();
        Plane::from_frame(Frame::new(p, f.z(), f.x()))
    }
}

impl Cylinder {
    /// Like [`new`](Self::new), without panicking.
    #[doc(alias = "gce_MakeCylinder")]
    pub fn try_new(frame: Frame, radius: f64) -> Result<Cylinder> {
        check_radius(radius)?;
        Ok(Self::new(frame, radius))
    }

    /// The cylinder of `radius` around `axis`, its frame at the axis' origin. Fails with
    /// [`ConstructionError::NegativeRadius`] for a negative radius.
    #[doc(alias = "gce_MakeCylinder")]
    pub fn from_axis(axis: Axis, radius: f64) -> Result<Cylinder> {
        check_radius(radius)?;
        Ok(Cylinder::new(constructed_frame(axis.origin, axis.direction), radius))
    }

    /// The cylinder whose axis runs through `p1` and `p2` and whose surface passes through
    /// `p3`, its frame at `p1` with Z towards `p2`. Fails with
    /// [`ConstructionError::ConfusedPoints`] when `p1` and `p2` coincide.
    #[doc(alias = "gce_MakeCylinder")]
    pub fn from_three_points(p1: Point, p2: Point, p3: Point) -> Result<Cylinder> {
        let axis = p2 - p1;
        if axis.length() < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let z = Direction::from_vector(axis);
        let radius = super::curve::Line::new(p1, z).distance_to_point(p3);
        Ok(Cylinder::new(constructed_frame(p1, z), radius))
    }
}

impl Cone {
    /// Like [`new`](Self::new), without panicking.
    #[doc(alias = "gce_MakeCone")]
    #[doc(alias = "GC_MakeConicalSurface")]
    pub fn try_new(frame: Frame, semi_angle: f64, radius: f64) -> Result<Cone> {
        let max = std::f64::consts::FRAC_PI_2 - crate::tolerance::RESOLUTION;
        if semi_angle.abs() <= crate::tolerance::RESOLUTION || semi_angle.abs() >= max {
            return Err(ConstructionError::BadAngle);
        }
        check_radius(radius)?;
        Ok(Self::new(frame, semi_angle, radius))
    }

    /// The parallel cone at distance `distance` (measured perpendicular to the generating
    /// lines; positive outwards).
    #[doc(alias = "gce_MakeCone")]
    pub fn offset(&self, distance: f64) -> Result<Cone> {
        let radius = self.radius() + distance / self.semi_angle().cos();
        check_radius(radius)?;
        Ok(Cone::new(*self.frame(), self.semi_angle(), radius))
    }

    /// The parallel cone through `p`. When two exist, the one nearer to `self` is chosen.
    #[doc(alias = "gce_MakeCone")]
    pub fn parallel_through(&self, p: Point) -> Result<Cone> {
        let f = self.frame();
        let v = (p - f.origin()).dot(Vector::from(f.z()));
        let radial = super::curve::Line::from(f.z_axis()).distance_to_point(p);
        let (tan, cos) = (self.semi_angle().tan(), self.semi_angle().cos());
        // The two cones through p (on either side of the axis), the nearer if both exist.
        let candidates = [radial - v * tan, -radial - v * tan];
        let usable: Vec<f64> = candidates.iter().filter(|&&r| r >= -crate::tolerance::RESOLUTION).map(|&r| r.max(0.0)).collect();
        let radius = usable
            .into_iter()
            .min_by(|a, b| ((a - self.radius()) * cos).abs().total_cmp(&((b - self.radius()) * cos).abs()))
            .ok_or(ConstructionError::NegativeRadius)?;
        Ok(Cone::new(*f, self.semi_angle(), radius))
    }

    /// The cone whose axis runs from `p1` to `p2`, with radius `r1` at `p1` and `r2` at
    /// `p2`. Its frame is at `p1` with Z towards `p2`, and its reference radius is `r1`.
    /// Fails with [`ConstructionError::NullAxis`] when the points coincide,
    /// [`ConstructionError::NegativeRadius`] for a negative radius and
    /// [`ConstructionError::BadAngle`] when the radii are equal (a cylinder) or the cone
    /// would be flat.
    #[doc(alias = "gce_MakeCone")]
    pub fn from_points_radii(p1: Point, p2: Point, r1: f64, r2: f64) -> Result<Cone> {
        let axis = p2 - p1;
        let length = axis.length();
        if length < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::NullAxis);
        }
        if r1 < 0.0 || r2 < 0.0 {
            return Err(ConstructionError::NegativeRadius);
        }
        let angle = check_cone_angle(((r1 - r2) / length).atan().abs())?;
        let z = Direction::from_vector(axis);
        let semi_angle = if r1 > r2 { -angle } else { angle };
        Ok(Cone::new(constructed_frame(p1, z), semi_angle, r1))
    }

    /// The cone with its axis on the line through `p1` and `p2` whose surface passes
    /// through `p3` and `p4`. Its frame is at the foot of `p3` on the axis with Z towards
    /// the foot of `p4`, and its reference radius is the distance of `p3` from the axis.
    ///
    /// The heights of the feet are measured as unsigned distances from `p1`, so the
    /// result only passes through both points when `p1` does not lie between the feet.
    /// Fails with [`ConstructionError::ConfusedPoints`] when `p1` and `p2`, or `p3` and
    /// `p4`, coincide, [`ConstructionError::BadAngle`] when the feet are equally far from
    /// `p1` or the cone would be a cylinder or flat, and [`ConstructionError::NullAxis`]
    /// when the feet coincide.
    #[doc(alias = "gce_MakeCone")]
    pub fn from_four_points(p1: Point, p2: Point, p3: Point, p4: Point) -> Result<Cone> {
        if (p2 - p1).length() < crate::tolerance::RESOLUTION || (p4 - p3).length() < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::ConfusedPoints);
        }
        let a = Direction::from_vector(p2 - p1);
        let av = Vector::from(a);
        let foot = |p: Point| p1 + av * (p - p1).dot(av);
        let (f3, f4) = (foot(p3), foot(p4));
        let (h3, h4) = ((f3 - p1).length(), (f4 - p1).length());
        if (h3 - h4).abs() < crate::tolerance::RESOLUTION {
            return Err(ConstructionError::BadAngle);
        }
        let line = super::curve::Line::new(p1, a);
        let (d3, d4) = (line.distance_to_point(p3), line.distance_to_point(p4));
        let angle = check_cone_angle(((d3 - d4) / (h3 - h4)).atan().abs())?;
        let (r1, r2) = ((f3 - p3).length(), (f4 - p4).length());
        let z = Direction::try_from_vector(f4 - f3).ok_or(ConstructionError::NullAxis)?;
        let semi_angle = if r1 > r2 { -angle } else { angle };
        Ok(Cone::new(constructed_frame(f3, z), semi_angle, r1))
    }
}

/// `angle` if it is usable as a cone's (unsigned) semi-angle: neither (near) zero nor
/// (near) a right angle.
fn check_cone_angle(angle: f64) -> Result<f64> {
    let res = crate::tolerance::RESOLUTION;
    if (std::f64::consts::FRAC_PI_2 - angle).abs() < res || angle < res {
        Err(ConstructionError::BadAngle)
    } else {
        Ok(angle)
    }
}

/// The frame of a cylinder or cone built from points: origin `origin`, main axis `z` and
/// the X axis these constructors use by convention. It lies in the world XY plane whenever
/// `z` has an x component, otherwise in the world YZ plane if `z` has a y component, and is
/// world X when `z` is along world Z.
fn constructed_frame(origin: Point, z: Direction) -> Frame {
    let res = crate::tolerance::RESOLUTION;
    let (p, q, r) = (z.x(), z.y(), z.z());
    let hint = if p.abs() > res {
        Direction::new(-q, p, 0.0)
    } else if q.abs() > res {
        Direction::new(0.0, -r, q)
    } else {
        Direction::new(r, 0.0, -p)
    };
    Frame::new(origin, z, hint)
}

/// The parameter (angle) of a point of a circle.
fn circle_parameter(circle: &Circle, p: Point) -> f64 {
    let f = circle.frame();
    let d = p - f.origin();
    d.dot(Vector::from(f.y())).atan2(d.dot(Vector::from(f.x()))).rem_euclid(std::f64::consts::TAU)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::ParametricSurface;

    fn xyz(d: Direction) -> (f64, f64, f64) {
        (d.x(), d.y(), d.z())
    }

    #[test]
    fn constructed_x_axes() {
        let s2 = std::f64::consts::FRAC_1_SQRT_2;
        let cases = [
            ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            ((0.0, 1.0, 1.0), (0.0, -s2, s2)),
            ((1.0, 1.0, 0.0), (-s2, s2, 0.0)),
        ];
        for (z, x) in cases {
            let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Direction::new(z.0, z.1, z.2));
            let c = Cylinder::from_axis(axis, 2.0).unwrap();
            let got = xyz(c.frame().x());
            for (g, e) in [(got.0, x.0), (got.1, x.1), (got.2, x.2)] {
                assert!((g - e).abs() <= 1e-16, "X for z = {z:?}: {got:?}");
            }
            assert_eq!(c.frame().origin(), axis.origin);
            assert!(c.frame().is_right_handed());
        }
        assert_eq!(Cylinder::from_axis(Axis::new(Point::ORIGIN, Direction::Z), -1.0), Err(ConstructionError::NegativeRadius));
    }

    #[test]
    fn cylinder_through_three_points() {
        let c = Cylinder::from_three_points(Point::ORIGIN, Point::new(0.0, 0.0, 7.0), Point::new(3.0, 4.0, 0.0)).unwrap();
        assert_eq!(c.radius(), 5.0);
        assert_eq!(xyz(c.frame().z()), (0.0, 0.0, 1.0));
        assert_eq!(xyz(c.frame().x()), (1.0, 0.0, 0.0));
        let p = Point::new(1.0, 1.0, 1.0);
        assert_eq!(Cylinder::from_three_points(p, p, Point::ORIGIN), Err(ConstructionError::ConfusedPoints));
    }

    #[test]
    fn cone_from_points_and_radii() {
        let (p1, p2) = (Point::ORIGIN, Point::new(0.0, 0.0, 10.0));
        let cone = Cone::from_points_radii(p1, p2, 4.0, 2.0).unwrap();
        assert!((cone.semi_angle() + 0.2f64.atan()).abs() < 1e-15);
        assert_eq!(cone.radius(), 4.0);
        assert_eq!(xyz(cone.frame().x()), (1.0, 0.0, 0.0));
        assert_eq!(cone.frame().origin(), p1);
        let (u, v) = cone.parameters(Point::new(2.0, 0.0, 10.0));
        assert!(cone.point(u, v).distance(Point::new(2.0, 0.0, 10.0)) < 1e-12);
        let growing = Cone::from_points_radii(p1, p2, 2.0, 4.0).unwrap();
        assert!(growing.semi_angle() > 0.0);
        assert_eq!(Cone::from_points_radii(p1, p2, 3.0, 3.0), Err(ConstructionError::BadAngle));
        assert_eq!(Cone::from_points_radii(p1, p1, 3.0, 2.0), Err(ConstructionError::NullAxis));
        assert_eq!(Cone::from_points_radii(p1, p2, 3.0, -2.0), Err(ConstructionError::NegativeRadius));
    }

    #[test]
    fn cone_from_four_points() {
        let (p1, p2) = (Point::new(1.0, 1.0, 0.0), Point::new(1.0, 1.0, 1.0));
        let (p3, p4) = (Point::new(4.0, 1.0, 2.0), Point::new(1.0, 2.0, 5.0));
        let cone = Cone::from_four_points(p1, p2, p3, p4).unwrap();
        for p in [p3, p4] {
            let (u, v) = cone.parameters(p);
            assert!(cone.point(u, v).distance(p) < 1e-12, "{p:?}");
        }
        assert_eq!(cone.frame().origin(), Point::new(1.0, 1.0, 2.0));
        assert_eq!(cone.radius(), 3.0);
        assert!(cone.semi_angle() < 0.0);
        // Feet equally far from p1, on either side: no angle.
        let p5 = Point::new(2.0, 1.0, -2.0);
        assert_eq!(Cone::from_four_points(p1, p2, p3, p5), Err(ConstructionError::BadAngle));
        assert_eq!(Cone::from_four_points(p1, p2, p3, p3), Err(ConstructionError::ConfusedPoints));
        assert_eq!(Cone::from_four_points(p1, p1, p3, p4), Err(ConstructionError::ConfusedPoints));
    }
}
