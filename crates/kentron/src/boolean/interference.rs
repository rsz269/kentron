//! How two solids meet: the volume they share and the gap between them.

use crate::error::Result;
use crate::topology::{Model, NodeKind, Shape};

/// How two solids meet, as they stand: the volume they share, and the gap between them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interference {
    /// The volume of their intersection (a solid wholly inside the other shares its own).
    pub volume: f64,
    /// The least distance from a point of one to a point of the other, insides included:
    /// zero where they overlap or touch.
    pub gap: f64,
}

/// How solids `a` and `b` meet: overlapping, a volume above zero and no gap; touching, no
/// volume and no gap; apart, no volume and the gap between them.
pub fn interference(model: &mut Model, a: Shape, b: Shape) -> Result<Interference> {
    let common = super::intersection(model, a, b)?;
    if model.explore(common, NodeKind::Solid).next().is_some() {
        let volume = crate::properties::volume(model, common).mass;
        return Ok(Interference { volume, gap: 0.0 });
    }
    let gap = crate::distance::distance(model, a, b)?.value;
    Ok(Interference { volume: 0.0, gap })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Axis, Direction, Frame, Point, Transform, Vector};
    use crate::primitives::{Cuboid, Cylinder, Sphere};
    use std::f64::consts::PI;

    /// A turn about a skew axis, so that nothing lies along the coordinate axes.
    fn turn() -> Transform {
        let axis = Axis::new(Point::new(1.0, 2.0, -0.5), Direction::try_from_vector(Vector::new(0.3, -0.7, 0.5)).unwrap());
        Transform::rotation(axis, 0.83)
    }

    fn at(x: f64, y: f64, z: f64) -> Frame {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X)
    }

    fn close(found: f64, expected: f64, relative: f64) {
        assert!((found - expected).abs() <= relative * expected.abs().max(1e-300), "{found} != {expected}");
    }

    /// A 4 x 3 x 2 box and a cube of 3 overlapping it from several places, both turned: they
    /// share the box their spans have in common.
    #[test]
    fn overlapping_boxes_share_the_box_their_spans_have_in_common() {
        let mut model = Model::new();
        let t = turn();
        for o in [[1.3, 0.7, -2.1], [3.2, 2.2, 1.2], [2.0, 1.0, 0.5]] {
            let a = Cuboid::new(4.0, 3.0, 2.0).build(&mut model).unwrap().moved(&t);
            let b = Cuboid::new(3.0, 3.0, 3.0).with_frame(at(o[0], o[1], o[2])).build(&mut model).unwrap().moved(&t);
            let span = |lo: f64, length: f64| length.min(lo + 3.0) - lo.max(0.0);
            let i = interference(&mut model, a, b).unwrap();
            close(i.volume, span(o[0], 4.0) * span(o[1], 3.0) * span(o[2], 2.0), 1e-12);
            assert_eq!(i.gap, 0.0);
        }
    }

    /// A 10 x 20 x 30 box with a corner at the origin and a cylinder of radius 5 and height 20
    /// on the z axis share a quarter of the cylinder: 125π.
    #[test]
    fn a_box_and_a_cylinder_on_its_edge_share_a_quarter_of_the_cylinder() {
        let mut model = Model::new();
        let a = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let b = Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
        let i = interference(&mut model, a, b).unwrap();
        close(i.volume, 125.0 * PI, 1e-12);
        assert!((i.volume - 392.699081699).abs() < 1e-9);
        assert_eq!(i.gap, 0.0);
    }

    /// Cubes face to face share nothing and are no distance apart.
    #[test]
    fn boxes_face_to_face_only_touch() {
        let mut model = Model::new();
        let t = turn();
        let a = Cuboid::new(3.0, 3.0, 3.0).build(&mut model).unwrap().moved(&t);
        let b = Cuboid::new(3.0, 3.0, 3.0).with_frame(at(3.0, 0.0, 0.0)).build(&mut model).unwrap().moved(&t);
        assert_eq!(interference(&mut model, a, b).unwrap(), Interference { volume: 0.0, gap: 0.0 });
    }

    /// Balls of radii 1.5 and 2 with their centres 3.5 + g apart, cubes g apart face to face
    /// and a ball g off a cube's face: the gap g, from a hair to a millimetre.
    #[test]
    fn solids_apart_are_their_gap_apart() {
        let mut model = Model::new();
        let t = turn();
        for g in [1e-6, 1e-4, 1e-3, 1.0] {
            let a = Sphere::new(1.5).build(&mut model).unwrap().moved(&t);
            let b = Sphere::new(2.0).with_frame(at(3.5 + g, 0.0, 0.0)).build(&mut model).unwrap().moved(&t);
            let i = interference(&mut model, a, b).unwrap();
            assert_eq!(i.volume, 0.0);
            close(i.gap, g, 1e-8);
            let a = Cuboid::new(3.0, 3.0, 3.0).build(&mut model).unwrap().moved(&t);
            let b = Cuboid::new(3.0, 3.0, 3.0).with_frame(at(3.0 + g, 0.0, 0.0)).build(&mut model).unwrap().moved(&t);
            let i = interference(&mut model, a, b).unwrap();
            assert_eq!(i.volume, 0.0);
            close(i.gap, g, 1e-8);
            let b = Sphere::new(1.5).with_frame(at(-1.5 - g, 1.5, 1.5)).build(&mut model).unwrap().moved(&t);
            let i = interference(&mut model, a, b).unwrap();
            assert_eq!(i.volume, 0.0);
            close(i.gap, g, 1e-8);
        }
    }

    /// A ball of radius 1 inside a cube of 6 shares its own volume, either way round.
    #[test]
    fn a_solid_inside_the_other_shares_its_own_volume() {
        let mut model = Model::new();
        let t = turn();
        let a = Cuboid::new(6.0, 6.0, 6.0).build(&mut model).unwrap().moved(&t);
        let b = Sphere::new(1.0).with_frame(at(3.0, 3.0, 3.0)).build(&mut model).unwrap().moved(&t);
        for (x, y) in [(a, b), (b, a)] {
            let i = interference(&mut model, x, y).unwrap();
            close(i.volume, 4.0 / 3.0 * PI, 1e-12);
            assert_eq!(i.gap, 0.0);
        }
    }

    /// Balls of radii 1.5 and 2 overlapping share a lens.
    #[test]
    fn overlapping_balls_share_a_lens() {
        let mut model = Model::new();
        let t = turn();
        let (ra, rb) = (1.5f64, 2.0f64);
        for d in [0.8, 1.7, 2.6, 3.3] {
            let a = Sphere::new(ra).build(&mut model).unwrap().moved(&t);
            let b = Sphere::new(rb).with_frame(at(d, 0.0, 0.0)).build(&mut model).unwrap().moved(&t);
            let lens = PI * (ra + rb - d).powi(2) * (d * d + 2.0 * d * (ra + rb) - 3.0 * (ra - rb).powi(2)) / (12.0 * d);
            let i = interference(&mut model, a, b).unwrap();
            close(i.volume, lens, 1e-9);
            assert_eq!(i.gap, 0.0);
        }
    }

    /// Cylinders of radius r on axes crossing square share Steinmetz's solid: 16r³/3.
    #[test]
    fn crossing_cylinders_share_steinmetz_solid() {
        let mut model = Model::new();
        let t = turn();
        let r = 1.25;
        let a = Cylinder::new(r, 6.0).with_frame(at(0.0, 0.0, -3.0)).build(&mut model).unwrap().moved(&t);
        let across = Frame::new(Point::new(-3.0, 0.0, 0.0), Direction::X, Direction::Y);
        let b = Cylinder::new(r, 6.0).with_frame(across).build(&mut model).unwrap().moved(&t);
        let i = interference(&mut model, a, b).unwrap();
        close(i.volume, 16.0 * r.powi(3) / 3.0, 1e-9);
        assert_eq!(i.gap, 0.0);
    }

    /// A shaft pressed into a tube of bore 1, from a hair over the bore up: they share a
    /// thin ring and no gap (their meshes, coarser than the press, need not cross).
    #[test]
    fn a_shaft_pressed_into_a_tube_shares_a_thin_ring() {
        let mut model = Model::new();
        let t = turn();
        for press in [1.005, 1.01, 1.05] {
            let outer = Cylinder::new(2.0, 4.0).build(&mut model).unwrap();
            let bore = Cylinder::new(1.0, 4.0).build(&mut model).unwrap();
            let tube = crate::boolean::difference(&mut model, outer, bore).unwrap().moved(&t);
            let shaft = Cylinder::new(press, 6.0).with_frame(at(0.0, 0.0, -1.0)).build(&mut model).unwrap().moved(&t);
            let i = interference(&mut model, tube, shaft).unwrap();
            close(i.volume, PI * (press * press - 1.0) * 4.0, 1e-12);
            assert_eq!(i.gap, 0.0);
        }
    }
}
