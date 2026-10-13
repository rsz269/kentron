//! gp_Trsf_Test.cxx (gp_Trsf is `Transform`)

use kentron::math::{Axis, Direction, Point, Transform, Transformable};
use kentron::tolerance;

/// OCC23361: applying a composed transformation equals applying its factors in turn.
#[test]
fn occ23361_transformation_composition() {
    let p = Point::new(0.0, 0.0, 2.0);
    let t1 = Transform::rotation(Axis::new(p, Direction::Y), -0.49328285294022267);
    let t2 = Transform::rotation(Axis::new(p, Direction::Z), 0.87538474718473880);

    let p1 = Point::new(10.0, 3.0, 4.0);
    let p2 = p1.transformed(&(t2 * t1));
    let p3 = p1.transformed(&t1).transformed(&t2);

    assert!(p2.approx_eq(p3, tolerance::LINEAR));
}
