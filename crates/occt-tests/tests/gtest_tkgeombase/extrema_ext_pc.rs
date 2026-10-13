//! Extrema_ExtPC_Test.cxx.

use kentron::geometry::curve::{Circle, Curve};
use kentron::geometry::extrema;
use kentron::geometry::surface::Cylinder;
use kentron::math::{Direction, Frame};

use crate::{assert_near, p};

/// Bug24945: the projection onto a circle, then its cylinder parameters.
#[test]
fn cylinder_parameter_normalization() {
    let point = p(-1725.97, 843.257, -4.22741e-13);
    let circle = Curve::Circle(Circle::new(
        Frame::new(p(0.0, 843.257, 0.0), -Direction::Y, Direction::X),
        1725.9708621929999,
    ));
    let projections = extrema::point_curve(point, &circle);
    assert!(!projections.is_empty());
    let projected = projections[0].point;
    assert_near(projected.x, -1725.97, 1.0e-2);
    assert_near(projected.y, 843.257, 1.0e-2);
    assert_near(projected.z, 0.0, 1.0e-10);

    let cylinder = Cylinder::new(Frame::new(p(0.0, 2103.87, 0.0), -Direction::Y, -Direction::X), 1890.0);
    let (u, v) = cylinder.parameters(projected);
    assert_near(u, 0.0, 1.0e-4);
    assert_near(v, 1260.613, 1.0e-2);
}
