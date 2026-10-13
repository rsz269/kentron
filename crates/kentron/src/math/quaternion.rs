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

    pub fn from_euler_angles(sequence: EulerSequence, alpha: f64, beta: f64, gamma: f64) -> Self {
        // Ken Shoemake's conversion (Graphics Gems IV), as OCCT's gp_Quaternion.
        let o = sequence.parameters();
        let (mut a, mut b, mut c) = (alpha, beta, gamma);
        if !o.extrinsic {
            (a, c) = (gamma, alpha);
        }
        if o.odd {
            b = -b;
        }
        let (si, ci) = (0.5 * a).sin_cos();
        let (sj, cj) = (0.5 * b).sin_cos();
        let (sh, ch) = (0.5 * c).sin_cos();
        let (cc, cs, sc, ss) = (ci * ch, ci * sh, si * ch, si * sh);
        let mut v = [0.0; 4]; // w, x, y, z
        if o.two_axes {
            v[o.i] = cj * (cs + sc);
            v[o.j] = sj * (cc + ss);
            v[o.k] = sj * (cs - sc);
            v[0] = cj * (cc - ss);
        } else {
            v[o.i] = cj * sc - sj * cs;
            v[o.j] = cj * ss + sj * cc;
            v[o.k] = cj * cs - sj * sc;
            v[0] = cj * cc + sj * ss;
        }
        if o.odd {
            v[o.j] = -v[o.j];
        }
        Self::new(v[1], v[2], v[3], v[0])
    }

    /// The angles `(alpha, beta, gamma)` of this rotation in `sequence`.
    pub fn euler_angles(&self, sequence: EulerSequence) -> (f64, f64, f64) {
        let m = self.to_matrix();
        let at = |r: usize, c: usize| m[(r - 1, c - 1)];
        let o = sequence.parameters();
        let (mut alpha, mut beta, mut gamma);
        if o.two_axes {
            let sy = at(o.i, o.j).hypot(at(o.i, o.k));
            if sy > 16.0 * f64::EPSILON {
                alpha = at(o.i, o.j).atan2(at(o.i, o.k));
                gamma = at(o.j, o.i).atan2(-at(o.k, o.i));
            } else {
                alpha = (-at(o.j, o.k)).atan2(at(o.j, o.j));
                gamma = 0.0;
            }
            beta = sy.atan2(at(o.i, o.i));
        } else {
            let cy = at(o.i, o.i).hypot(at(o.j, o.i));
            if cy > 16.0 * f64::EPSILON {
                alpha = at(o.k, o.j).atan2(at(o.k, o.k));
                gamma = at(o.j, o.i).atan2(at(o.i, o.i));
            } else {
                alpha = (-at(o.j, o.k)).atan2(at(o.j, o.j));
                gamma = 0.0;
            }
            beta = (-at(o.k, o.i)).atan2(cy);
        }
        if o.odd {
            (alpha, beta, gamma) = (-alpha, -beta, -gamma);
        }
        if !o.extrinsic {
            (alpha, gamma) = (gamma, alpha);
        }
        (alpha, beta, gamma)
    }

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

/// How OCCT's conversion reads an Euler sequence: the axes (1-based, in order), whether
/// the first two make an odd permutation, whether the third axis is the first again, and
/// whether the rotations are about fixed axes.
struct SequenceParameters {
    i: usize,
    j: usize,
    k: usize,
    odd: bool,
    two_axes: bool,
    extrinsic: bool,
}

impl EulerSequence {
    fn parameters(self) -> SequenceParameters {
        use EulerSequence::*;
        let p = |first: usize, odd: bool, two_axes: bool, extrinsic: bool| SequenceParameters {
            i: first,
            j: 1 + (first + usize::from(odd)) % 3,
            k: 1 + (first + usize::from(!odd)) % 3,
            odd,
            two_axes,
            extrinsic,
        };
        // Intrinsic rotations are the extrinsic ones in reverse order (with alpha and gamma
        // swapped by the conversion).
        match self {
            ExtrinsicXYZ => p(1, false, false, true),
            ExtrinsicXZY => p(1, true, false, true),
            ExtrinsicYZX => p(2, false, false, true),
            ExtrinsicYXZ => p(2, true, false, true),
            ExtrinsicZXY => p(3, false, false, true),
            ExtrinsicZYX => p(3, true, false, true),
            IntrinsicXYZ => p(3, true, false, false),
            IntrinsicXZY => p(2, false, false, false),
            IntrinsicYZX => p(1, true, false, false),
            IntrinsicYXZ => p(3, false, false, false),
            IntrinsicZXY => p(2, true, false, false),
            IntrinsicZYX | YawPitchRoll => p(1, false, false, false),
            ExtrinsicXYX => p(1, false, true, true),
            ExtrinsicXZX => p(1, true, true, true),
            ExtrinsicYZY => p(2, false, true, true),
            ExtrinsicYXY => p(2, true, true, true),
            ExtrinsicZXZ => p(3, false, true, true),
            ExtrinsicZYZ => p(3, true, true, true),
            IntrinsicXYX => p(1, false, true, false),
            IntrinsicXZX => p(1, true, true, false),
            IntrinsicYZY => p(2, false, true, false),
            IntrinsicYXY => p(2, true, true, false),
            IntrinsicZXZ | Classic => p(3, false, true, false),
            IntrinsicZYZ => p(3, true, true, false),
        }
    }
}
