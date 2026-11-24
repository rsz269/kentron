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

    /// The rotation by the Euler angles `alpha`, `beta`, `gamma` in `sequence`, the angles
    /// taken about the sequence's axes in the order named. Extrinsic `PQS`: about the fixed
    /// axes, `alpha` about `P` first, so `q_S(γ)·q_Q(β)·q_P(α)`; intrinsic `PQS`: about the
    /// moving axes, `q_P(α)·q_Q(β)·q_S(γ)` (`q_A(t)` the turn by `t` about world axis `A`).
    /// The product as it comes, not sign-normalised.
    #[doc(alias = "SetEulerAngles")]
    pub fn from_euler_angles(sequence: EulerSequence, alpha: f64, beta: f64, gamma: f64) -> Self {
        let (extrinsic, [p, q, s]) = sequence.axes();
        let turn = |axis: usize, angle: f64| {
            let (sin, cos) = (0.5 * angle).sin_cos();
            let mut v = [0.0; 3];
            v[axis] = sin;
            Self::new(v[0], v[1], v[2], cos)
        };
        match extrinsic {
            true => turn(s, gamma).product(&turn(q, beta)).product(&turn(p, alpha)),
            false => turn(p, alpha).product(&turn(q, beta)).product(&turn(s, gamma)),
        }
    }

    /// The Euler angles `(alpha, beta, gamma)` of this rotation in `sequence` (the inverse of
    /// [`Self::from_euler_angles`]; the quaternion need not be of unit length). `alpha` and
    /// `gamma` are in `[−π, π]`; `beta` in `[−π/2, π/2]` for sequences of three different
    /// axes, in `[0, π]` for two-axis sequences whose middle axis follows the first cyclically
    /// (`XYX`, `YZY`, `ZXZ`) and in `[−π, 0]` for the others (`XZX`, `YXY`, `ZYZ`). At gimbal
    /// lock the angle of the leftmost factor of the product (extrinsic `gamma`, intrinsic
    /// `alpha`) is zero and the other carries the whole turn about the shared axis.
    #[doc(alias = "GetEulerAngles")]
    pub fn euler_angles(&self, sequence: EulerSequence) -> (f64, f64, f64) {
        let (extrinsic, [p, q, s]) = sequence.axes();
        // As a product R_A(a)·R_B(b)·R_C(c) from left to right.
        let (a_axis, b_axis, c_axis) = if extrinsic { (s, q, p) } else { (p, q, s) };
        let m = self.to_matrix();
        // Coordinates permuted so that A is x and B is y; an odd permutation (a reflection)
        // turns every rotation the other way.
        let k = 3 - a_axis - b_axis;
        let perm = [a_axis, b_axis, if a_axis == c_axis { k } else { c_axis }];
        let sign = if b_axis == (a_axis + 1) % 3 { 1.0 } else { -1.0 };
        let at = |r: usize, c: usize| m[(perm[r], perm[c])];
        let lock = 16.0 * f64::EPSILON;
        let (a, b, c) = if a_axis == c_axis {
            // R_x(a)·R_y(b)·R_x(c).
            let sin_b = at(0, 1).hypot(at(0, 2));
            let b = sin_b.atan2(at(0, 0));
            if sin_b <= lock { (0.0, b, (-at(1, 2)).atan2(at(1, 1))) } else { (at(1, 0).atan2(-at(2, 0)), b, at(0, 1).atan2(at(0, 2))) }
        } else {
            // R_x(a)·R_y(b)·R_z(c).
            let cos_b = at(0, 0).hypot(at(0, 1));
            let b = at(0, 2).atan2(cos_b);
            if cos_b <= lock { (0.0, b, at(1, 0).atan2(at(1, 1))) } else { ((-at(1, 2)).atan2(at(2, 2)), b, (-at(0, 1)).atan2(at(0, 0))) }
        };
        let (a, b, c) = (sign * a, sign * b, sign * c);
        if extrinsic { (c, b, a) } else { (a, b, c) }
    }

    /// The Hamilton product `self · other` (`other` applied first).
    fn product(&self, other: &Self) -> Self {
        let (a, b) = (self, other);
        Self::new(
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        )
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

impl EulerSequence {
    /// Whether the sequence is extrinsic, and its three axes in the order named
    /// (0, 1, 2 for X, Y, Z).
    fn axes(self) -> (bool, [usize; 3]) {
        use EulerSequence::*;
        const X: usize = 0;
        const Y: usize = 1;
        const Z: usize = 2;
        match self {
            Classic | IntrinsicZXZ => (false, [Z, X, Z]),
            YawPitchRoll | IntrinsicZYX => (false, [Z, Y, X]),
            ExtrinsicXYZ => (true, [X, Y, Z]),
            ExtrinsicXZY => (true, [X, Z, Y]),
            ExtrinsicYZX => (true, [Y, Z, X]),
            ExtrinsicYXZ => (true, [Y, X, Z]),
            ExtrinsicZXY => (true, [Z, X, Y]),
            ExtrinsicZYX => (true, [Z, Y, X]),
            IntrinsicXYZ => (false, [X, Y, Z]),
            IntrinsicXZY => (false, [X, Z, Y]),
            IntrinsicYZX => (false, [Y, Z, X]),
            IntrinsicYXZ => (false, [Y, X, Z]),
            IntrinsicZXY => (false, [Z, X, Y]),
            ExtrinsicXYX => (true, [X, Y, X]),
            ExtrinsicXZX => (true, [X, Z, X]),
            ExtrinsicYZY => (true, [Y, Z, Y]),
            ExtrinsicYXY => (true, [Y, X, Y]),
            ExtrinsicZYZ => (true, [Z, Y, Z]),
            ExtrinsicZXZ => (true, [Z, X, Z]),
            IntrinsicXYX => (false, [X, Y, X]),
            IntrinsicXZX => (false, [X, Z, X]),
            IntrinsicYZY => (false, [Y, Z, Y]),
            IntrinsicYXY => (false, [Y, X, Y]),
            IntrinsicZYZ => (false, [Z, Y, Z]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn all() -> Vec<EulerSequence> {
        let mut v = EulerSequence::EXPLICIT.to_vec();
        v.extend([EulerSequence::Classic, EulerSequence::YawPitchRoll]);
        v
    }

    #[test]
    fn euler_angles_are_the_axis_products() {
        let turn = |axis: usize, t: f64| {
            let mut v = Vector::new(0.0, 0.0, 0.0);
            match axis {
                0 => v.x = 1.0,
                1 => v.y = 1.0,
                _ => v.z = 1.0,
            }
            Quaternion::from_axis_angle(v, t)
        };
        for s in all() {
            let (extrinsic, [p, q, r]) = s.axes();
            for (a, b, c) in [(0.3, -1.1, 2.5), (-2.9, 0.7, 0.2), (7.0, 4.0, -5.0)] {
                let got = Quaternion::from_euler_angles(s, a, b, c);
                let want = if extrinsic { turn(r, c).product(&turn(q, b)).product(&turn(p, a)) } else { turn(p, a).product(&turn(q, b)).product(&turn(r, c)) };
                for (x, y) in [(got.x, want.x), (got.y, want.y), (got.z, want.z), (got.w, want.w)] {
                    assert!((x - y).abs() < 1e-15, "{s:?}");
                }
                // Back to angles in range, the same rotation.
                let (a2, b2, c2) = got.euler_angles(s);
                assert!(a2.abs() <= PI && c2.abs() <= PI);
                if p == r {
                    let forward = q == (p + 1) % 3;
                    assert!(if forward { (0.0..=PI).contains(&b2) } else { (-PI..=0.0).contains(&b2) }, "{s:?} {b2}");
                } else {
                    assert!(b2.abs() <= FRAC_PI_2);
                }
                let (m1, m2) = (got.to_matrix(), Quaternion::from_euler_angles(s, a2, b2, c2).to_matrix());
                for i in 0..3 {
                    for j in 0..3 {
                        assert!((m1[(i, j)] - m2[(i, j)]).abs() < 1e-14, "{s:?}");
                    }
                }
                // A quaternion not of unit length gives the same angles.
                let scaled = Quaternion::new(got.x * 3.0, got.y * 3.0, got.z * 3.0, got.w * 3.0);
                let (a3, b3, c3) = scaled.euler_angles(s);
                assert!((a3 - a2).abs() < 1e-12 && (b3 - b2).abs() < 1e-12 && (c3 - c2).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn gimbal_lock_puts_the_turn_in_the_rightmost_factor() {
        use EulerSequence::*;
        let cases = [
            (ExtrinsicXYZ, FRAC_PI_2, (-0.5, FRAC_PI_2, 0.0)),
            (ExtrinsicXYZ, -FRAC_PI_2, (1.3, -FRAC_PI_2, 0.0)),
            (IntrinsicZYX, FRAC_PI_2, (0.0, FRAC_PI_2, 0.5)),
            (ExtrinsicZXZ, 0.0, (1.3, 0.0, 0.0)),
            (ExtrinsicZXZ, PI, (-0.5, PI, 0.0)),
            (IntrinsicZYZ, PI, (0.0, -PI, 0.5)),
        ];
        for (s, beta, (a, b, c)) in cases {
            let (x, y, z) = Quaternion::from_euler_angles(s, 0.4, beta, 0.9).euler_angles(s);
            assert!((x - a).abs() < 1e-12 && (y - b).abs() < 1e-12 && (z - c).abs() < 1e-12, "{s:?} {beta}: {x} {y} {z}");
        }
    }
}
