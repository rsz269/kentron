use super::{Matrix3, Vector};

/// Conventions for Euler angles: which axes the three rotations are about, and whether
/// they are about the fixed (extrinsic) or the rotating (intrinsic) axes.
#[doc(alias = "gp_EulerSequence")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EulerSequence {
    /// Classic Euler angles: same as `IntrinsicZXZ`.
    Classic,
    /// Yaw, pitch, roll: same as `IntrinsicZYX`.
    YawPitchRoll,
    ExtrinsicXYZ,
    ExtrinsicXZY,
    ExtrinsicYZX,
    ExtrinsicYXZ,
    ExtrinsicZXY,
    ExtrinsicZYX,
    IntrinsicXYZ,
    IntrinsicXZY,
    IntrinsicYZX,
    IntrinsicYXZ,
    IntrinsicZXY,
    IntrinsicZYX,
    ExtrinsicXYX,
    ExtrinsicXZX,
    ExtrinsicYZY,
    ExtrinsicYXY,
    ExtrinsicZYZ,
    ExtrinsicZXZ,
    IntrinsicXYX,
    IntrinsicXZX,
    IntrinsicYZY,
    IntrinsicYXY,
    IntrinsicZXZ,
    IntrinsicZYZ,
}

impl EulerSequence {
    /// The 24 sequences that name their axes explicitly (all but `Classic` and
    /// `YawPitchRoll`).
    pub const EXPLICIT: [EulerSequence; 24] = {
        use EulerSequence::*;
        [
            ExtrinsicXYZ, ExtrinsicXZY, ExtrinsicYZX, ExtrinsicYXZ, ExtrinsicZXY, ExtrinsicZYX,
            IntrinsicXYZ, IntrinsicXZY, IntrinsicYZX, IntrinsicYXZ, IntrinsicZXY, IntrinsicZYX,
            ExtrinsicXYX, ExtrinsicXZX, ExtrinsicYZY, ExtrinsicYXY, ExtrinsicZYZ, ExtrinsicZXZ,
            IntrinsicXYX, IntrinsicXZX, IntrinsicYZY, IntrinsicYXY, IntrinsicZXZ, IntrinsicZYZ,
        ]
    };
}

/// A rotation quaternion `w + xi + yj + zk`.
#[doc(alias = "gp_Quaternion")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Default for Quaternion {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Quaternion {
    pub const IDENTITY: Quaternion = Quaternion::new(0.0, 0.0, 0.0, 1.0);

    pub const fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    /// Rotation by `angle` (radians) about `axis`.
    pub fn from_axis_angle(axis: Vector, angle: f64) -> Self {
        let n = axis.length();
        let a = if n > 0.0 { axis * (1.0 / n) } else { Vector::new(0.0, 0.0, 1.0) };
        let (sin, cos) = (0.5 * angle).sin_cos();
        Self::new(a.x * sin, a.y * sin, a.z * sin, cos)
    }

    // CLEAN-ROOM GAP: 67 lines removed; rewrite from the cleanroom spec `quaternion_euler.md`.

    pub fn to_matrix(&self) -> Matrix3 {
        let Self { x, y, z, w } = *self;
        let s = 2.0 / (x * x + y * y + z * z + w * w);
        let (x2, y2, z2) = (x * s, y * s, z * s);
        let (xx, xy, xz, yy, yz, zz) = (x * x2, x * y2, x * z2, y * y2, y * z2, z * z2);
        let (wx, wy, wz) = (w * x2, w * y2, w * z2);
        Matrix3::from_rows([
            [1.0 - (yy + zz), xy - wz, xz + wy],
            [xy + wz, 1.0 - (xx + zz), yz - wx],
            [xz - wy, yz + wx, 1.0 - (xx + yy)],
        ])
    }

    /// Rotates `v`.
    pub fn rotate(&self, v: Vector) -> Vector {
        let Self { x, y, z, w } = *self;
        let q = Self::new(v.x * w + v.z * y - v.y * z, v.y * w + v.x * z - v.z * x, v.z * w + v.y * x - v.x * y, v.x * x + v.y * y + v.z * z);
        let n = x * x + y * y + z * z + w * w;
        Vector::new(
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y + y * q.w + z * q.x - x * q.z,
            w * q.z + z * q.w + x * q.y - y * q.x,
        ) * (1.0 / n)
    }
}

impl Quaternion {
    /// The rotation of a rotation matrix.
    pub(crate) fn from_matrix(m: &Matrix3) -> Self {
        let at = |r: usize, c: usize| m[(r - 1, c - 1)];
        let trace = at(1, 1) + at(2, 2) + at(3, 3);
        let (q, big) = if trace > 0.0 {
            let q = Self::new(at(3, 2) - at(2, 3), at(1, 3) - at(3, 1), at(2, 1) - at(1, 2), trace + 1.0);
            (q, q.w)
        } else if at(1, 1) > at(2, 2) && at(1, 1) > at(3, 3) {
            let q = Self::new(1.0 + at(1, 1) - at(2, 2) - at(3, 3), at(1, 2) + at(2, 1), at(1, 3) + at(3, 1), at(3, 2) - at(2, 3));
            (q, q.x)
        } else if at(2, 2) > at(3, 3) {
            let q = Self::new(at(1, 2) + at(2, 1), 1.0 + at(2, 2) - at(1, 1) - at(3, 3), at(2, 3) + at(3, 2), at(1, 3) - at(3, 1));
            (q, q.y)
        } else {
            let q = Self::new(at(1, 3) + at(3, 1), at(2, 3) + at(3, 2), 1.0 + at(3, 3) - at(1, 1) - at(2, 2), at(2, 1) - at(1, 2));
            (q, q.z)
        };
        let s = 0.5 / big.sqrt();
        Self::new(q.x * s, q.y * s, q.z * s, q.w * s)
    }
}

// CLEAN-ROOM GAP: 53 lines removed; rewrite from the cleanroom spec `quaternion_euler.md`.
