//! gp_Ax3_Test.cxx (gp_Ax3 is `Frame`)

use kentron::math::{Axis, Direction, Frame, Point};
use kentron::tolerance;

/// OCC29406: re-orienting a frame's axes keeps its handedness.
#[test]
fn occ29406_set_direction_preserves_orientation() {
    let tol = tolerance::ANGULAR;

    // main (Z) direction
    {
        let mut f1 = Frame::WORLD;
        let mut f2 = Frame::WORLD;
        let mut f3 = Frame::WORLD;
        let mut f4 = Frame::WORLD;
        f3.flip_z();
        f4.flip_z();

        let right1 = f1.is_right_handed();
        f1.set_z(Direction::X);
        assert_eq!(right1, f1.is_right_handed(), "handedness should be preserved");
        assert!(Direction::X.approx_eq(f1.z(), tol));

        let right2 = f2.is_right_handed();
        f2.set_z(Direction::NEG_X);
        assert_eq!(right2, f2.is_right_handed());
        assert!(Direction::NEG_X.approx_eq(f2.z(), tol));

        let right3 = f3.is_right_handed();
        f3.set_z(Direction::X);
        assert_eq!(right3, f3.is_right_handed());
        assert!(Direction::X.approx_eq(f3.z(), tol));

        let right4 = f4.is_right_handed();
        f4.set_z(Direction::NEG_X);
        assert_eq!(right4, f4.is_right_handed());
        assert!(Direction::NEG_X.approx_eq(f4.z(), tol));

        // set_z_axis
        let mut f5 = Frame::WORLD;
        let mut f6 = Frame::WORLD;
        let axis1 = Axis::new(Point::ORIGIN, Direction::X);
        let axis2 = Axis::new(Point::ORIGIN, Direction::NEG_X);

        let right5 = f5.is_right_handed();
        f5.set_z_axis(axis1);
        assert_eq!(right5, f5.is_right_handed());
        assert!(axis1.direction.approx_eq(f5.z(), tol));

        let right6 = f6.is_right_handed();
        f6.set_z_axis(axis2);
        assert_eq!(right6, f6.is_right_handed());
        assert!(axis2.direction.approx_eq(f6.z(), tol));
    }

    // X direction
    {
        let mut f1 = Frame::WORLD;
        let mut f2 = Frame::WORLD;
        let mut f3 = Frame::WORLD;
        let mut f4 = Frame::WORLD;
        f3.flip_x();
        f4.flip_x();

        let right1 = f1.is_right_handed();
        f1.set_x(Direction::Z);
        assert_eq!(right1, f1.is_right_handed());
        let good_y1 = Direction::from_vector(f1.z().cross(Direction::Z));
        if f1.is_right_handed() {
            assert!(good_y1.approx_eq(f1.y(), tol));
        } else {
            assert!(good_y1.is_opposite(f1.y(), tol));
        }

        let right2 = f2.is_right_handed();
        f2.set_x(Direction::NEG_Z);
        assert_eq!(right2, f2.is_right_handed());

        let right3 = f3.is_right_handed();
        f3.set_x(Direction::Z);
        assert_eq!(right3, f3.is_right_handed());

        let right4 = f4.is_right_handed();
        f4.set_x(Direction::NEG_Z);
        assert_eq!(right4, f4.is_right_handed());
    }

    // Y direction
    {
        let mut f1 = Frame::WORLD;
        let mut f2 = Frame::WORLD;
        let mut f3 = Frame::WORLD;
        let mut f4 = Frame::WORLD;
        f3.flip_y();
        f4.flip_y();

        let right1 = f1.is_right_handed();
        f1.set_y(Direction::Z);
        assert_eq!(right1, f1.is_right_handed());
        let good_x1 = Direction::from_vector(f1.z().cross(Direction::Z));
        if f1.is_right_handed() {
            assert!(good_x1.is_opposite(f1.x(), tol));
        } else {
            assert!(good_x1.approx_eq(f1.x(), tol));
        }

        let right2 = f2.is_right_handed();
        f2.set_y(Direction::NEG_Z);
        assert_eq!(right2, f2.is_right_handed());

        let right3 = f3.is_right_handed();
        f3.set_y(Direction::Z);
        assert_eq!(right3, f3.is_right_handed());

        let right4 = f4.is_right_handed();
        f4.set_y(Direction::NEG_Z);
        assert_eq!(right4, f4.is_right_handed());
    }
}
