//! Times boolean operations on growing models: a plate cut by a k × k grid of cylinders
//! (holes), and its section by a plane.
//!
//! `cargo run --release -p kentron --example boolean_bench [k...]` (default 4 8 16 24).

use std::time::Instant;

use kentron::boolean;
use kentron::math::{Direction, Frame, Point};
use kentron::primitives::{Cuboid, Cylinder};
use kentron::topology::{Model, NodeKind, Shape};
use kentron::{properties, validity};

/// A 100 × 100 × 2 plate and a compound of k × k cylinders through it.
fn plate_and_pins(model: &mut Model, k: usize) -> (Shape, Shape) {
    let plate = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(100.0, 100.0, 2.0)).build(model).unwrap();
    let pitch = 100.0 / k as f64;
    let pins = model.new_empty(NodeKind::Compound);
    for i in 0..k {
        for j in 0..k {
            let center = Point::new(pitch * (i as f64 + 0.5), pitch * (j as f64 + 0.5), -1.0);
            let pin = Cylinder::new(0.3 * pitch, 4.0).with_frame(Frame::from_z(center, Direction::Z)).build(model).unwrap();
            model.add(pins, pin).unwrap();
        }
    }
    (plate, pins)
}

fn main() {
    let ks: Vec<usize> = std::env::args().skip(1).map(|a| a.parse().expect("a grid size")).collect();
    let ks = if ks.is_empty() { vec![4, 8, 16, 24] } else { ks };
    println!("{:>4} {:>6} {:>10} {:>10} {:>12}", "k", "holes", "cut (s)", "faces", "volume");
    for k in ks {
        let mut model = Model::new();
        let (plate, pins) = plate_and_pins(&mut model, k);
        let start = Instant::now();
        let result = boolean::difference(&mut model, plate, pins).unwrap();
        let seconds = start.elapsed().as_secs_f64();
        assert!(validity::check(&model, result).is_valid());
        let faces = model.explore(result, NodeKind::Face).count();
        println!("{k:>4} {:>6} {seconds:>10.3} {faces:>10} {:>12.4}", k * k, properties::volume(&model, result).mass);
    }
}
