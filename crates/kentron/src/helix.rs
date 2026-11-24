//! Helices: wires winding round an axis, on cylinders or cones, in parts of their own
//! pitch, height and diameters.

use std::f64::consts::{FRAC_PI_2, TAU};

use crate::error::{Error, Result};
use crate::geometry::curve::BSplineCurve;
use crate::geometry::convert::CurveApproximationOptions;
use crate::geometry::curve::Curve;
use crate::geometry::{Continuity, ParametricCurve};
use crate::math::{Direction, Frame, Point, Transform, Transformable, Vector};
use crate::topology::{Model, Shape};

/// A helix about the z axis: at angle `t` (radians) it is `r(t)` from the axis, where
/// `r(t) = start_radius + rise·tan(taper)·t`, at height `rise·t` with `rise = pitch / 2π`
/// (turning counterclockwise seen from above, clockwise if `left_handed`).
#[doc(alias = "HelixGeom_HelixCurve")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelixCurve {
    pub first: f64,
    pub last: f64,
    pub pitch: f64,
    pub start_radius: f64,
    /// The angle of the cone it lies on to the axis (0 on a cylinder).
    pub taper: f64,
    pub left_handed: bool,
}

impl HelixCurve {
    fn rise(&self) -> f64 {
        self.pitch / TAU
    }

    /// How fast the radius grows with the angle.
    fn widening(&self) -> f64 {
        if self.taper.abs() > 1e-10 { self.rise() * self.taper.tan() } else { 0.0 }
    }

    /// The `n`-th derivative (`n` ≥ 0) at `t`.
    fn derivative(&self, t: f64, n: u32) -> Vector {
        // r·(cos, sin)(t): by Leibniz, r⁽ⁿ⁾ terms vanish beyond the first.
        let (k, r) = (self.widening(), self.start_radius + self.widening() * t);
        let shifted = |m: u32| t + m as f64 * FRAC_PI_2;
        let (cos_n, sin_n) = (shifted(n).cos(), shifted(n).sin());
        let (mut x, mut y) = (r * cos_n, r * sin_n);
        if n >= 1 {
            let (cos_m, sin_m) = (shifted(n - 1).cos(), shifted(n - 1).sin());
            x += n as f64 * k * cos_m;
            y += n as f64 * k * sin_m;
        }
        if self.left_handed {
            y = -y;
        }
        let z = match n {
            0 => self.rise() * t,
            1 => self.rise(),
            _ => 0.0,
        };
        Vector::new(x, y, z)
    }
}

impl ParametricCurve for HelixCurve {
    type Point = Point;
    type Vector = Vector;

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn is_closed(&self) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn continuity(&self) -> Continuity {
        Continuity::CN
    }

    fn is_cn(&self, _n: u32) -> bool {
        true
    }

    fn point(&self, t: f64) -> Point {
        Point::ORIGIN + self.derivative(t, 0)
    }

    fn d1(&self, t: f64) -> (Point, Vector) {
        (self.point(t), self.derivative(t, 1))
    }

    fn d2(&self, t: f64) -> (Point, Vector, Vector) {
        (self.point(t), self.derivative(t, 1), self.derivative(t, 2))
    }

    fn d3(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        (self.point(t), self.derivative(t, 1), self.derivative(t, 2), self.derivative(t, 3))
    }

    fn dn(&self, t: f64, n: u32) -> Vector {
        assert!(n > 0, "dn needs n > 0");
        self.derivative(t, n)
    }

    fn reversed_parameter(&self, t: f64) -> f64 {
        self.first + self.last - t
    }

    fn reverse(&mut self) {
        unimplemented!("reversing a helix curve")
    }

    fn resolution(&self, tolerance: f64) -> f64 {
        let widest = self.start_radius.max(self.start_radius + self.widening() * self.last).max(self.start_radius + self.widening() * self.first);
        let speed = (widest * widest + self.widening().powi(2) + self.rise().powi(2)).sqrt();
        tolerance / speed.max(1e-300)
    }
}

/// One part of a helical wire: rising `height` along the axis with `pitch` a turn, its
/// diameter changing linearly to `end_diameter`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelixPart {
    pub height: f64,
    pub pitch: f64,
    pub end_diameter: f64,
}

/// How a helix is approximated by B-splines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelixOptions {
    /// The largest allowed distance from the exact helix.
    pub tolerance: f64,
    pub max_degree: usize,
    pub max_segments: usize,
    pub left_handed: bool,
}

impl Default for HelixOptions {
    fn default() -> Self {
        Self { tolerance: 1e-4, max_degree: 8, max_segments: 1000, left_handed: false }
    }
}

/// A helical wire and how closely its B-splines follow the exact helix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Helix {
    pub wire: Shape,
    pub tolerance_reached: f64,
}

/// A helical wire about the z axis of `frame`, starting `start_diameter / 2` along its x
/// axis, in `parts` one above the other. Each part's turns are B-spline edges approximating
/// the exact helix (see [`HelixCurve`]); where parts meet, their edges are made to share
/// a tangent. Fails for a part with a non-positive height or pitch, a taper of a right
/// angle or more, or a start on the axis.
#[doc(alias = "HelixBRep_BuilderHelix")]
pub fn helix(model: &mut Model, frame: &Frame, start_diameter: f64, parts: &[HelixPart], options: &HelixOptions) -> Result<Helix> {
    if parts.is_empty() {
        return Err(Error::InvalidInput("a helix needs at least one part".into()));
    }
    let tolerance = options.tolerance;
    let approximation = CurveApproximationOptions { tolerance, continuity: Continuity::C1, max_segments: options.max_segments, max_degree: options.max_degree };
    let (origin, z) = (frame.origin(), Vector::from(frame.z()));
    let mut start = origin + Vector::from(frame.x()) * (start_diameter / 2.0);
    let mut diameter = start_diameter;
    let mut reached: f64 = 0.0;
    let mut edges: Vec<Vec<BSplineCurve>> = Vec::with_capacity(parts.len());
    for part in parts {
        let taper = (0.5 * (part.end_diameter - diameter) / part.height).atan();
        if !(taper.abs() < FRAC_PI_2 - 1e-12) {
            return Err(Error::InvalidInput(format!("a helix part tapers at {taper} rad, a right angle or more")));
        }
        if !(part.height >= tolerance) || !(part.pitch >= tolerance) {
            return Err(Error::InvalidInput(format!("a helix part of height {} and pitch {}", part.height, part.pitch)));
        }
        // The part's own frame: on the axis level with its start, x towards the start.
        let foot = origin + z * (start - origin).dot(z);
        let radial = start - foot;
        let radius = radial.length();
        if !(radius >= tolerance) {
            return Err(Error::InvalidInput("a helix part starts on its axis".into()));
        }
        let local = Frame::new(foot, frame.z(), Direction::from_vector(radial));
        let to_world = Transform::frame_to_frame(&Frame::WORLD, &local);
        let turns = part.height / part.pitch;
        let whole = turns.floor();
        let exact = |first: f64, last: f64| HelixCurve { first, last, pitch: part.pitch, start_radius: radius, taper, left_handed: options.left_handed };
        let approximate = |first: f64, last: f64| -> Result<(BSplineCurve, f64)> {
            let made = crate::geometry::convert::approximate_parametric(&exact(first, last), &approximation)
                .ok_or_else(|| Error::Failed("a helix turn could not be approximated".into()))?;
            Ok((made.curve.transformed(&to_world), made.max_error))
        };
        let mut part_edges = Vec::with_capacity(whole as usize + 1);
        let cylinder = taper == 0.0;
        for i in 0..whole as usize {
            let (a, b) = (TAU * i as f64, TAU * (i + 1) as f64);
            if cylinder && i > 0 {
                // On a cylinder every turn is the first one moved up by whole pitches.
                let first: &BSplineCurve = &part_edges[0];
                part_edges.push(first.transformed(&Transform::translation(z * (part.pitch * i as f64))));
                continue;
            }
            let (curve, error) = approximate(a, b)?;
            reached = reached.max(error);
            part_edges.push(curve);
        }
        if turns - whole > 1e-7 * turns.max(1.0) || whole == 0.0 {
            let (curve, error) = approximate(TAU * whole, TAU * turns)?;
            reached = reached.max(error);
            part_edges.push(curve);
        }
        let last = part_edges.last().expect("an edge");
        start = last.point(last.last_parameter());
        diameter = part.end_diameter;
        edges.push(part_edges);
    }
    for k in 1..edges.len() {
        let (before, after) = edges.split_at_mut(k);
        smooth(before[k - 1].last_mut().expect("an edge"), after[0].first_mut().expect("an edge"));
    }
    let curves: Vec<BSplineCurve> = edges.into_iter().flatten().collect();
    let mut vertex = model.add_vertex(curves[0].point(curves[0].first_parameter()), reached);
    let mut shapes = Vec::with_capacity(curves.len());
    for curve in curves {
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        let end = model.add_vertex(curve.point(b), reached);
        let edge = model.add_edge(vertex, end, reached)?;
        model.set_edge_curve(edge.into(), Curve::BSpline(curve), reached)?;
        model.set_edge_range(edge.into(), a, b)?;
        shapes.push(Shape::from(edge));
        vertex = end;
    }
    let wire = crate::make::wire(model, &shapes)?;
    Ok(Helix { wire, tolerance_reached: reached })
}

/// Gives the end of `before` and the start of `after` the mean of their two derivatives
/// there (each keeping its own junction point, `before` its first two poles and `after` its
/// last two), raising a curve's degree one step at a time while it has too little freedom;
/// a curve that cannot be moved is left as it was. Nothing is done where the tangents
/// already agree within `1e-7` rad.
fn smooth(before: &mut BSplineCurve, after: &mut BSplineCurve) {
    let (b, a) = (before.last_parameter(), after.first_parameter());
    let ((p, d1), (q, d2)) = (before.d1(b), after.d1(a));
    if d1.angle(d2) < 1e-7 {
        return;
    }
    let mean = (d1 + d2) * 0.5;
    let moved = |curve: &BSplineCurve, u: f64, point: Point, starting: i32, ending: i32| {
        (curve.degree()..=BSplineCurve::<Point>::MAX_DEGREE).find_map(|degree| {
            let mut trial = curve.clone();
            trial.increase_degree(degree);
            trial.move_point_and_tangent(u, point, mean, crate::tolerance::PARAMETRIC, starting, ending).map(|()| trial)
        })
    };
    if let Some(c) = moved(before, b, p, 1, -1) {
        *before = c;
    }
    if let Some(c) = moved(after, a, q, -1, 1) {
        *after = c;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::NodeKind;

    /// The edges' curves of a wire, in order.
    fn curves(model: &Model, wire: Shape) -> Vec<(Curve, f64, f64)> {
        model.explore(wire, NodeKind::Edge).filter_map(|e| model.edge_curve(e)).collect()
    }

    #[test]
    fn a_helix_follows_the_exact_one_a_turn_an_edge() {
        let mut model = Model::new();
        let part = HelixPart { height: 17.5, pitch: 5.0, end_diameter: 20.0 };
        let made = helix(&mut model, &Frame::WORLD, 20.0, &[part], &HelixOptions::default()).unwrap();
        assert!(made.tolerance_reached <= 1e-4);
        let edges = curves(&model, made.wire);
        assert_eq!(edges.len(), 4, "three turns and a half");
        let exact = HelixCurve { first: 0.0, last: TAU * 3.5, pitch: 5.0, start_radius: 10.0, taper: 0.0, left_handed: false };
        for (c, a, b) in &edges {
            for i in 0..=20 {
                let p = c.point(a + (b - a) * i as f64 / 20.0);
                // On the cylinder, at the height its angle gives (within the tolerance).
                let angle = p.y.atan2(p.x).rem_euclid(TAU);
                let turns = ((p.z / 5.0) - angle / TAU).round();
                let q = exact.point(angle + turns * TAU);
                assert!(p.distance(q) < 2e-4, "{p:?} {q:?}");
            }
        }
        // Ending 3.5 turns on: half way round, at the full height.
        let (c, _, b) = edges.last().unwrap();
        assert!(c.point(*b).distance(Point::new(-10.0, 0.0, 17.5)) < 1e-6);
    }

    #[test]
    fn spiral_parts_meet_smoothly_and_widen() {
        let mut model = Model::new();
        let parts = [HelixPart { height: 10.0, pitch: 4.0, end_diameter: 30.0 }, HelixPart { height: 10.0, pitch: 2.0, end_diameter: 40.0 }];
        let made = helix(&mut model, &Frame::WORLD, 20.0, &parts, &HelixOptions::default()).unwrap();
        let edges = curves(&model, made.wire);
        // Where each edge meets the next: the same point, nearly the same tangent direction
        // (that of the approximations), exactly where the parts meet (made so).
        for w in edges.windows(2) {
            let ((c1, _, l1), (c2, f2, _)) = (&w[0], &w[1]);
            let ((p1, v1), (p2, v2)) = (c1.d1(*l1), c2.d1(*f2));
            let kink = if (p1.z - 10.0).abs() < 1e-6 { 1e-9 } else { 1e-4 };
            assert!(p1.distance(p2) < 1e-9 && v1.angle(v2) < kink, "{p1:?} {p2:?} {}", v1.angle(v2));
        }
        let (c, _, b) = edges.last().unwrap();
        let end = c.point(*b);
        assert!((end.z - 20.0).abs() < 1e-6 && ((end.x * end.x + end.y * end.y).sqrt() - 20.0).abs() < 1e-3, "{end:?}");
        // A part with no height is refused.
        let flat = [HelixPart { height: 0.0, pitch: 4.0, end_diameter: 20.0 }];
        assert!(helix(&mut model, &Frame::WORLD, 20.0, &flat, &HelixOptions::default()).is_err());
    }

    #[test]
    fn edges_by_whole_turns_and_handedness() {
        let mut model = Model::new();
        let options = HelixOptions::default();
        // Three whole turns: three edges, no remainder; less than a turn: one.
        let made = helix(&mut model, &Frame::WORLD, 10.0, &[HelixPart { height: 6.0, pitch: 2.0, end_diameter: 10.0 }], &options).unwrap();
        assert_eq!(curves(&model, made.wire).len(), 3);
        let made = helix(&mut model, &Frame::WORLD, 10.0, &[HelixPart { height: 1.0, pitch: 2.0, end_diameter: 10.0 }], &options).unwrap();
        let edges = curves(&model, made.wire);
        assert_eq!(edges.len(), 1);
        let (c, _, b) = &edges[0];
        assert!(c.point(*b).distance(Point::new(-5.0, 0.0, 1.0)) < 1e-6);
        // Left-handed: clockwise seen from above.
        let left = HelixOptions { left_handed: true, ..options };
        let made = helix(&mut model, &Frame::WORLD, 10.0, &[HelixPart { height: 1.0, pitch: 4.0, end_diameter: 10.0 }], &left).unwrap();
        let edges = curves(&model, made.wire);
        let (c, a, b) = &edges[0];
        assert!(c.point(*b).distance(Point::new(0.0, -5.0, 1.0)) < 1e-6, "{:?}", c.point(*b));
        assert!(c.d1(*a).1.y < 0.0);
        // A cone narrowing, three parts with two junctions, both smooth.
        let parts = [
            HelixPart { height: 4.0, pitch: 2.0, end_diameter: 6.0 },
            HelixPart { height: 3.0, pitch: 1.5, end_diameter: 6.0 },
            HelixPart { height: 2.5, pitch: 1.0, end_diameter: 8.0 },
        ];
        let made = helix(&mut model, &Frame::WORLD, 10.0, &parts, &options).unwrap();
        let edges = curves(&model, made.wire);
        assert_eq!(edges.len(), 2 + 2 + 3);
        for w in edges.windows(2) {
            let ((c1, _, l1), (c2, f2, _)) = (&w[0], &w[1]);
            let ((p1, v1), (p2, v2)) = (c1.d1(*l1), c2.d1(*f2));
            // (Exactly where the parts meet; within the approximation between turns.)
            let kink = if (p1.z - 4.0).abs() < 1e-6 || (p1.z - 7.0).abs() < 1e-6 { 1e-9 } else { 1e-3 };
            assert!(p1.distance(p2) < 1e-9 && v1.angle(v2) < kink, "{p1:?} {}", v1.angle(v2));
        }
        let (c, _, b) = edges.last().unwrap();
        let end = c.point(*b);
        assert!((end.z - 9.5).abs() < 1e-6 && ((end.x * end.x + end.y * end.y).sqrt() - 4.0).abs() < 1e-3, "{end:?}");
    }

    #[test]
    fn invalid_helices_are_refused() {
        let mut model = Model::new();
        let options = HelixOptions::default();
        let part = |height: f64, pitch: f64, end_diameter: f64| [HelixPart { height, pitch, end_diameter }];
        assert!(helix(&mut model, &Frame::WORLD, 10.0, &[], &options).is_err());
        assert!(helix(&mut model, &Frame::WORLD, 10.0, &part(5.0, -1.0, 10.0), &options).is_err());
        assert!(helix(&mut model, &Frame::WORLD, 0.0, &part(5.0, 1.0, 0.0), &options).is_err());
        assert!(helix(&mut model, &Frame::WORLD, 10.0, &part(1e-30, 1.0, 20.0), &options).is_err());
    }
}
