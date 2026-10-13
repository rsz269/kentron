//! gp_Quaternion_Test.cxx

use std::f64::consts::PI;

use occt_tests::assert_near;
use kentron::math::{
    Direction, EulerSequence, Frame, Matrix3, Point, Quaternion, Transform, Transformable, Vector,
};
use kentron::tolerance;

/// The axis letters of an explicit sequence, e.g. b"XYZ" for `ExtrinsicXYZ`.
fn axes(seq: EulerSequence) -> [u8; 3] {
    let name = format!("{seq:?}");
    let letters = &name.as_bytes()[name.len() - 3..];
    [letters[0], letters[1], letters[2]]
}

fn sample_quaternion() -> Quaternion {
    Quaternion::new(0.06766916507860499, 0.21848101129786085, 0.11994599260380681, 0.9660744746954637)
}

/// OCC25574: Euler angles round-trip for every sequence.
#[test]
fn occ25574_euler_angles_consistency() {
    let quat = sample_quaternion();
    let r_inv = quat.to_matrix().inverse().unwrap();

    for seq in EulerSequence::EXPLICIT {
        let (alpha, beta, gamma) = quat.euler_angles(seq);
        let quat2 = Quaternion::from_euler_angles(seq, alpha, beta, gamma);

        let diff = quat2.to_matrix() * r_inv - Matrix3::IDENTITY;
        assert!(diff.determinant() <= 1e-5, "round-trip failed for Euler sequence {seq:?}");
    }
}

/// OCC25574: a rotation about one of a sequence's axes leaves that axis fixed.
#[test]
fn occ25574_euler_axis_rotation_preserves_axis() {
    for seq in EulerSequence::EXPLICIT {
        for j in 0..3 {
            let axis = (axes(seq)[j] - b'X') as usize;
            assert!(axis <= 2);

            // 90-degree rotation about the j-th axis of the sequence
            let mut angles = [0.0; 3];
            angles[j] = 0.5 * PI;
            let q2 = Quaternion::from_euler_angles(seq, angles[0], angles[1], angles[2]);

            // unit vector on the rotation axis
            let mut v = Vector::ZERO;
            v[axis] = 1.0;
            let v2 = v.transformed(&Transform::from_quaternion(q2));

            assert!(
                (v - v2).length_squared() <= tolerance::LINEAR_SQUARED,
                "rotation about an axis must not move a point on it: {seq:?}, angle index {j}"
            );
        }
    }
}

/// OCC25574: extrinsic and intrinsic sequences with reversed axes give reversed angles.
#[test]
fn occ25574_extrinsic_intrinsic_correspondence() {
    use EulerSequence::*;
    let alpha = 0.1517461713131;
    let beta = 1.5162198410141;
    let gamma = 1.9313156236541;

    let pairs = [
        (ExtrinsicXYZ, IntrinsicZYX),
        (ExtrinsicXZY, IntrinsicYZX),
        (ExtrinsicYZX, IntrinsicXZY),
        (ExtrinsicYXZ, IntrinsicZXY),
        (ExtrinsicZXY, IntrinsicYXZ),
        (ExtrinsicZYX, IntrinsicXYZ),
    ];

    for (extrinsic, intrinsic) in pairs {
        let quat = Quaternion::from_euler_angles(extrinsic, alpha, beta, gamma);
        let (gamma2, beta2, alpha2) = quat.euler_angles(intrinsic);

        assert_near!(alpha, alpha2, 1e-5, "alpha mismatch for {extrinsic:?}/{intrinsic:?}");
        assert_near!(beta, beta2, 1e-5, "beta mismatch for {extrinsic:?}/{intrinsic:?}");
        assert_near!(gamma, gamma2, 1e-5, "gamma mismatch for {extrinsic:?}/{intrinsic:?}");
    }
}

/// OCC25574: angles recovered as yaw-pitch-roll from a frame built by three rotations.
#[test]
fn occ25574_yaw_pitch_roll_round_trip() {
    let world = Frame::WORLD;
    let alpha = 0.0;
    let beta = -35.0 / 180.0 * PI;
    let gamma = 90.0 / 180.0 * PI;

    // build the rotated frame step by step (yaw-pitch-roll)
    let rot_z = Quaternion::from_axis_angle(world.z().into(), alpha);
    let rot_y = rot_z.rotate(world.y().into());
    let rot_x = rot_z.rotate(world.x().into());

    let rot_yaw = Quaternion::from_axis_angle(rot_y, beta);
    let rot_z2 = rot_yaw.rotate(world.z().into());
    let rot_x2 = rot_yaw.rotate(rot_x);

    let rot_roll = Quaternion::from_axis_angle(rot_x2, gamma);
    let rot_z3 = rot_roll.rotate(rot_z2);

    let result =
        Frame::new(Point::ORIGIN, Direction::from_vector(rot_z3), Direction::from_vector(rot_x2));

    let t = Transform::frame_to_frame(&Frame::WORLD, &result);
    let (a, b, g) = t.rotation_part().euler_angles(EulerSequence::YawPitchRoll);
    assert_near!(alpha, a, 1e-5, "yaw-pitch-roll alpha mismatch");
    assert_near!(beta, b, 1e-5, "yaw-pitch-roll beta mismatch");
    assert_near!(gamma, g, 1e-5, "yaw-pitch-roll gamma mismatch");
}

/// OCC25574 (issue 25946): IntrinsicZYX and ExtrinsicXYZ give the same angles reversed.
#[test]
fn occ25574_intrinsic_zyx_vs_extrinsic_xyz() {
    let quat = sample_quaternion();
    let (alpha, beta, gamma) = quat.euler_angles(EulerSequence::IntrinsicZYX);
    let (alpha2, beta2, gamma2) = quat.euler_angles(EulerSequence::ExtrinsicXYZ);

    assert_near!(alpha, gamma2, 1e-5, "intrinsic ZYX alpha should equal extrinsic XYZ gamma");
    assert_near!(beta, beta2, 1e-5, "beta should match");
    assert_near!(gamma, alpha2, 1e-5, "intrinsic ZYX gamma should equal extrinsic XYZ alpha");
}
