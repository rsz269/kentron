use std::ops::Mul;

use super::{
    Axis, Axis2, Direction, Direction2, Frame, Matrix3, Point, Point2, Quaternion, Vector,
    Vector2,
};

/// A similarity transformation of 3D space: rotation and/or reflection, uniform scaling,
/// then translation. `a * b` applies `b` first.
#[doc(alias = "gp_Trsf")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    matrix: Matrix3,
    scale: f64,
    translation: Vector,
}

impl Default for Transform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform {
    pub const IDENTITY: Transform =
        Transform { matrix: Matrix3::IDENTITY, scale: 1.0, translation: Vector::ZERO };

    pub fn translation(v: Vector) -> Self {
        Transform { translation: v, ..Self::IDENTITY }
    }

    /// Rotation by `angle` (radians, counter-clockwise looking against the axis direction)
    /// about `axis`.
    pub fn rotation(axis: Axis, angle: f64) -> Self {
        let k = Vector::from(axis.direction);
        let (sin, cos) = angle.sin_cos();
        let mut m = [[0.0; 3]; 3];
        let cross = [[0.0, -k.z, k.y], [k.z, 0.0, -k.x], [-k.y, k.x, 0.0]];
        for r in 0..3 {
            for c in 0..3 {
                let identity = if r == c { 1.0 } else { 0.0 };
                m[r][c] = cos * identity + sin * cross[r][c] + (1.0 - cos) * k[r] * k[c];
            }
        }
        let matrix = Matrix3::from_rows(m);
        let o = axis.origin - Point::ORIGIN;
        Transform { matrix, scale: 1.0, translation: o - matrix * o }
    }

    /// Rotation about the origin.
    pub fn from_quaternion(q: Quaternion) -> Self {
        Transform { matrix: q.to_matrix(), scale: 1.0, translation: Vector::ZERO }
    }

    /// Uniform scaling about `center`.
    pub fn scaling(center: Point, factor: f64) -> Self {
        let c = center - Point::ORIGIN;
        Transform { matrix: Matrix3::IDENTITY, scale: factor, translation: c - c * factor }
    }

    /// Point reflection (central symmetry) through `center`.
    pub fn mirror_point(center: Point) -> Self {
        Self::scaling(center, -1.0)
    }

    /// Reflection about a line (a half-turn about `axis`).
    pub fn mirror_axis(axis: Axis) -> Self {
        Self::rotation(axis, std::f64::consts::PI)
    }

    /// Reflection across the plane through `origin` with normal `normal`.
    pub fn mirror_plane(origin: Point, normal: Direction) -> Self {
        // A half-turn about the normal combined with a point reflection (scale -1), as OCCT
        // represents it: the matrix stays a rotation.
        let half_turn = Self::rotation(Axis::new(origin, normal), std::f64::consts::PI);
        Self::mirror_point(origin) * half_turn
    }

    /// Maps frame `from` onto frame `to`: geometry expressed relative to `from` ends up at
    /// the same relative position with respect to `to`.
    pub fn frame_to_frame(from: &Frame, to: &Frame) -> Self {
        let axes = |f: &Frame| {
            let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
            Matrix3::from_rows([[x.x, y.x, z.x], [x.y, y.y, z.y], [x.z, y.z, z.z]])
        };
        let mut matrix = axes(to) * axes(from).transpose();
        let mut scale = 1.0;
        if matrix.determinant() < 0.0 {
            // Frames of opposite handedness: keep the matrix a rotation.
            matrix = Matrix3::ZERO - matrix;
            scale = -1.0;
        }
        let from_origin = from.origin() - Point::ORIGIN;
        let translation = (to.origin() - Point::ORIGIN) - matrix * from_origin * scale;
        Transform { matrix, scale, translation }
    }

    /// The transformation `p ↦ M p + t` given by the rows `[m0, m1, m2, t]` of a 3×4
    /// matrix whose 3×3 part `M` is a rotation times a uniform scale (negative for a
    /// reflection). `None` if `M` is singular.
    #[doc(alias = "gp_Trsf::SetValues")]
    pub fn from_rows(rows: [[f64; 4]; 3]) -> Option<Self> {
        let m = |r: usize| [rows[r][0], rows[r][1], rows[r][2]];
        let full = Matrix3::from_rows([m(0), m(1), m(2)]);
        let det = full.determinant();
        if det.abs() < 1e-30 {
            return None;
        }
        let scale = det.signum() * det.abs().cbrt();
        let unit = |r: usize| [rows[r][0] / scale, rows[r][1] / scale, rows[r][2] / scale];
        let matrix = Matrix3::from_rows([unit(0), unit(1), unit(2)]);
        Some(Transform { matrix, scale, translation: Vector::new(rows[0][3], rows[1][3], rows[2][3]) })
    }

    /// The rows `[m0, m1, m2, t]` of the 3×4 matrix of `p ↦ M p + t` (the inverse of
    /// [`from_rows`](Self::from_rows)).
    #[doc(alias = "gp_Trsf::Value")]
    pub fn rows(&self) -> [[f64; 4]; 3] {
        let m = |r: usize, c: usize| self.scale * self.matrix[(r, c)];
        let t = self.translation;
        [[m(0, 0), m(0, 1), m(0, 2), t.x], [m(1, 0), m(1, 1), m(1, 2), t.y], [m(2, 0), m(2, 1), m(2, 2), t.z]]
    }

    /// The transformation applied `n` times (its inverse `-n` times for a negative `n`).
    pub fn powi(&self, n: i32) -> Transform {
        let base = if n < 0 { self.inverse() } else { *self };
        (0..n.unsigned_abs()).fold(Transform::IDENTITY, |t, _| t * base)
    }

    /// The uniform scale factor (negative for a point reflection or a reflection across a
    /// plane, which OCCT represents as a rotation combined with scaling by -1).
    pub fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// True if the transformation reverses orientation (a reflection, or a scaling by a
    /// negative factor).
    pub fn is_negative(&self) -> bool {
        self.matrix.determinant() * self.scale < 0.0
    }

    /// Where the origin goes.
    pub fn translation_part(&self) -> Vector {
        self.translation
    }

    /// The rotational part.
    pub fn rotation_part(&self) -> Quaternion {
        Quaternion::from_matrix(&self.matrix)
    }

    /// Panics if the transformation is not invertible (zero scale).
    pub fn inverse(&self) -> Transform {
        assert!(self.scale != 0.0, "inverse of a singular transformation");
        let matrix = self.matrix.transpose();
        let scale = 1.0 / self.scale;
        Transform { matrix, scale, translation: -(matrix * self.translation) * scale }
    }

    pub fn is_identity(&self) -> bool {
        const EPS: f64 = 1e-12;
        (self.scale - 1.0).abs() <= EPS
            && self.translation.approx_eq(Vector::ZERO, EPS)
            && (0..3).all(|r| (0..3).all(|c| {
                (self.matrix[(r, c)] - Matrix3::IDENTITY[(r, c)]).abs() <= EPS
            }))
    }
}

/// `a * b` applies `b` first, then `a` (as matrices compose).
impl Mul for Transform {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Transform {
        Transform {
            matrix: self.matrix * rhs.matrix,
            scale: self.scale * rhs.scale,
            translation: self.matrix * rhs.translation * self.scale + self.translation,
        }
    }
}

/// Something that can be moved by a [`Transform`]. Only [`transformed`](Self::transformed)
/// is required; the rest are conveniences.
pub trait Transformable: Sized {
    fn transformed(&self, t: &Transform) -> Self;

    fn translated(&self, v: Vector) -> Self {
        self.transformed(&Transform::translation(v))
    }

    fn rotated(&self, axis: Axis, angle: f64) -> Self {
        self.transformed(&Transform::rotation(axis, angle))
    }

    fn scaled(&self, center: Point, factor: f64) -> Self {
        self.transformed(&Transform::scaling(center, factor))
    }

    fn mirrored_point(&self, center: Point) -> Self {
        self.transformed(&Transform::mirror_point(center))
    }

    fn mirrored_axis(&self, axis: Axis) -> Self {
        self.transformed(&Transform::mirror_axis(axis))
    }

    fn mirrored_plane(&self, origin: Point, normal: Direction) -> Self {
        self.transformed(&Transform::mirror_plane(origin, normal))
    }
}

impl Transformable for Point {
    fn transformed(&self, t: &Transform) -> Self {
        Point::ORIGIN + t.matrix * (*self - Point::ORIGIN) * t.scale + t.translation
    }
}

/// Vectors ignore the translation part.
impl Transformable for Vector {
    fn transformed(&self, t: &Transform) -> Self {
        t.matrix * *self * t.scale
    }
}

/// Directions ignore translation and scaling.
impl Transformable for Direction {
    fn transformed(&self, t: &Transform) -> Self {
        let v = t.matrix * Vector::from(*self);
        Direction::from_vector(if t.scale < 0.0 { -v } else { v })
    }
}

impl Transformable for Axis {
    fn transformed(&self, t: &Transform) -> Self {
        Axis::new(self.origin.transformed(t), self.direction.transformed(t))
    }
}

impl Transformable for Frame {
    fn transformed(&self, t: &Transform) -> Self {
        Frame::from_axes_unchecked(
            self.origin().transformed(t),
            self.x().transformed(t),
            self.y().transformed(t),
            self.z().transformed(t),
        )
    }
}

/// A similarity transformation of the plane. `a * b` applies `b` first.
#[doc(alias = "gp_Trsf2d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform2 {
    matrix: [[f64; 2]; 2],
    scale: f64,
    translation: Vector2,
}

impl Default for Transform2 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform2 {
    pub const IDENTITY: Transform2 =
        Transform2 { matrix: [[1.0, 0.0], [0.0, 1.0]], scale: 1.0, translation: Vector2::ZERO };

    pub fn translation(v: Vector2) -> Self {
        Transform2 { translation: v, ..Self::IDENTITY }
    }

    /// The transformation with linear part `matrix` (orthogonal) and `scale`, keeping
    /// `center` fixed.
    fn about(center: Point2, matrix: [[f64; 2]; 2], scale: f64) -> Self {
        let c = Vector2::new(center.x, center.y);
        let moved = mat2(matrix, c) * scale;
        Transform2 { matrix, scale, translation: c - moved }
    }

    /// Counter-clockwise rotation by `angle` (radians) about `center`.
    pub fn rotation(center: Point2, angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self::about(center, [[c, -s], [s, c]], 1.0)
    }

    pub fn scaling(center: Point2, factor: f64) -> Self {
        Self::about(center, [[1.0, 0.0], [0.0, 1.0]], factor)
    }

    pub fn mirror_point(center: Point2) -> Self {
        Self::scaling(center, -1.0)
    }

    /// Reflection across the line `axis`.
    pub fn mirror_axis(axis: Axis2) -> Self {
        let (x, y) = (axis.direction.x(), axis.direction.y());
        let matrix = [[x * x - y * y, 2.0 * x * y], [2.0 * x * y, y * y - x * x]];
        Self::about(axis.origin, matrix, 1.0)
    }

    /// The uniform scale factor (negative for a point reflection, which is a scaling by
    /// -1).
    pub fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// True if the transformation reverses orientation (a reflection across an axis).
    pub fn is_negative(&self) -> bool {
        self.matrix[0][0] * self.matrix[1][1] - self.matrix[0][1] * self.matrix[1][0] < 0.0
    }

    pub fn is_identity(&self) -> bool {
        *self == Self::IDENTITY
    }

    /// The linear part (rotation times scale, possibly a reflection), row-major.
    pub fn matrix(&self) -> [[f64; 2]; 2] {
        let m = self.matrix;
        [[m[0][0] * self.scale, m[0][1] * self.scale], [m[1][0] * self.scale, m[1][1] * self.scale]]
    }

    pub fn translation_part(&self) -> Vector2 {
        self.translation
    }

    /// Panics if the transformation is singular (zero scale).
    pub fn inverse(&self) -> Transform2 {
        assert!(self.scale != 0.0, "inverse of a singular transformation");
        let m = self.matrix;
        let matrix = [[m[0][0], m[1][0]], [m[0][1], m[1][1]]];
        let scale = 1.0 / self.scale;
        Transform2 { matrix, scale, translation: mat2(matrix, self.translation) * -scale }
    }

    /// `self` applied `n` times (the inverse `-n` times for negative `n`; identity for 0).
    pub fn powi(&self, n: i32) -> Transform2 {
        let base = if n < 0 { self.inverse() } else { *self };
        (0..n.unsigned_abs()).fold(Self::IDENTITY, |acc, _| acc * base)
    }
}

/// `matrix · v`.
fn mat2(m: [[f64; 2]; 2], v: Vector2) -> Vector2 {
    Vector2::new(m[0][0] * v.x + m[0][1] * v.y, m[1][0] * v.x + m[1][1] * v.y)
}

/// `a * b` applies `b` first, then `a`.
impl Mul for Transform2 {
    type Output = Transform2;
    fn mul(self, rhs: Transform2) -> Transform2 {
        let (a, b) = (self.matrix, rhs.matrix);
        let matrix = [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ];
        Transform2 {
            matrix,
            scale: self.scale * rhs.scale,
            translation: mat2(a, rhs.translation) * self.scale + self.translation,
        }
    }
}

/// Something that can be moved by a [`Transform2`].
pub trait Transformable2: Sized {
    fn transformed(&self, t: &Transform2) -> Self;

    fn translated(&self, v: Vector2) -> Self {
        self.transformed(&Transform2::translation(v))
    }

    fn rotated(&self, center: Point2, angle: f64) -> Self {
        self.transformed(&Transform2::rotation(center, angle))
    }

    fn scaled(&self, center: Point2, factor: f64) -> Self {
        self.transformed(&Transform2::scaling(center, factor))
    }

    fn mirrored_point(&self, center: Point2) -> Self {
        self.transformed(&Transform2::mirror_point(center))
    }

    fn mirrored_axis(&self, axis: Axis2) -> Self {
        self.transformed(&Transform2::mirror_axis(axis))
    }
}

impl Transformable2 for Point2 {
    fn transformed(&self, t: &Transform2) -> Self {
        let v = mat2(t.matrix, Vector2::new(self.x, self.y)) * t.scale + t.translation;
        Point2::new(v.x, v.y)
    }
}

impl Transformable2 for Vector2 {
    fn transformed(&self, t: &Transform2) -> Self {
        mat2(t.matrix, *self) * t.scale
    }
}

/// Directions ignore translation and scaling (but turn round for a negative scale).
impl Transformable2 for Direction2 {
    fn transformed(&self, t: &Transform2) -> Self {
        let v = mat2(t.matrix, Vector2::from(*self));
        Direction2::from_vector(if t.scale < 0.0 { v * -1.0 } else { v })
    }
}

impl Transformable2 for Axis2 {
    fn transformed(&self, t: &Transform2) -> Self {
        Axis2::new(self.origin.transformed(t), self.direction.transformed(t))
    }
}

/// A general affine transformation of 3D space (may be non-uniform, e.g. scaling along one
/// direction only).
#[doc(alias = "gp_GTrsf")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineTransform {
    matrix: Matrix3,
    translation: Vector,
}

impl Default for AffineTransform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl AffineTransform {
    pub const IDENTITY: AffineTransform =
        AffineTransform { matrix: Matrix3::IDENTITY, translation: Vector::ZERO };

    /// `p -> matrix * p + translation`.
    pub const fn new(matrix: Matrix3, translation: Vector) -> Self {
        Self { matrix, translation }
    }

    /// The largest factor by which the map stretches any direction (the largest singular
    /// value of its matrix).
    pub fn max_stretch(&self) -> f64 {
        // The square root of the largest eigenvalue of MᵀM, by power iteration.
        let m = self.matrix.transpose() * self.matrix;
        let mut v = Vector::new(0.5773, 0.5774, 0.5775);
        let mut eigenvalue = 0.0;
        for _ in 0..200 {
            let w = m * v;
            let length = w.length();
            if length == 0.0 {
                return 0.0;
            }
            eigenvalue = length / v.length();
            v = w * (1.0 / length);
        }
        eigenvalue.sqrt()
    }

    /// Scales distances to the plane through `origin` with normal `normal` by `ratio` (0
    /// projects onto the plane).
    pub fn plane_affinity(origin: Point, normal: Direction, ratio: f64) -> Self {
        // p -> p - (1 - ratio) n (n·(p - origin)).
        let n = Vector::from(normal);
        let k = 1.0 - ratio;
        let row = |i: usize| {
            let ni = [n.x, n.y, n.z][i];
            let id = |j: usize| if i == j { 1.0 } else { 0.0 };
            [id(0) - k * ni * n.x, id(1) - k * ni * n.y, id(2) - k * ni * n.z]
        };
        let matrix = Matrix3::from_rows([row(0), row(1), row(2)]);
        let translation = n * (k * n.dot(origin - Point::ORIGIN));
        Self { matrix, translation }
    }

    pub fn transform_point(&self, p: Point) -> Point {
        Point::ORIGIN + self.matrix * (p - Point::ORIGIN) + self.translation
    }

    /// The linear part.
    pub fn matrix(&self) -> Matrix3 {
        self.matrix
    }

    pub fn translation(&self) -> Vector {
        self.translation
    }

    /// The direction a surface normal `n` takes: by the inverse transpose of the linear
    /// part (the cofactors, where the part is singular); `None` where it collapses.
    pub fn transform_normal(&self, n: Direction) -> Option<Direction> {
        let det = self.matrix.determinant();
        let v = self.matrix.cofactors() * Vector::from(n);
        let v = if det < 0.0 { v * -1.0 } else { v };
        Direction::try_from_vector(v).filter(|_| v.length() > 1e-300)
    }

    /// Transforms a displacement (ignoring the translation part).
    pub fn transform_vector(&self, v: Vector) -> Vector {
        self.matrix * v
    }
}

impl From<Transform> for AffineTransform {
    fn from(t: Transform) -> Self {
        let rows = t.rows();
        let matrix = Matrix3::from_rows(std::array::from_fn(|i| [rows[i][0], rows[i][1], rows[i][2]]));
        Self { matrix, translation: Vector::new(rows[0][3], rows[1][3], rows[2][3]) }
    }
}

/// Something that can be deformed by an [`AffineTransform`].
pub trait AffineTransformable: Sized {
    fn affine_transformed(&self, t: &AffineTransform) -> Self;
}

impl AffineTransformable for Point {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        t.transform_point(*self)
    }
}
