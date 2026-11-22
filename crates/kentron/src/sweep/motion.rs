//! What a solid passes through as it moves rigidly, cut from a blank.
//!
//! The set a moving tool fills is bounded by its surface where the motion runs along it (its
//! characteristic curves, swept), by its sharp edges swept, and by the tool where the motion
//! starts and ends. Three kinds of motion make that a plane problem, and those are the ones
//! [`cut_by_motion`] takes:
//!
//! * **Planar**, a prism moved in the plane square to its edges (turns about axes square to
//!   the plane, slides along it): the cut is the same at every height, so its section is the
//!   plane region the prism's outline sweeps, bounded by the paths of its convex corners, the
//!   envelopes of its sides (where a side's point moves along the side: the foot on it of the
//!   motion's instantaneous centre) and the outline where the motion starts and ends.
//! * **A screw** (a turn with a slide along its axis) carrying a solid of revolution clean
//!   past a rod about the same axis: what it fills there is screwed into itself, so its
//!   section square to the axis is the same at every height, turned. A point of the tool
//!   projects onto the section along the helix through it; the section's boundary is where
//!   the tool's characteristic curve (where the screw's velocity is tangent to the tool) and
//!   its rims project.
//! * **A turn** carrying a solid of revolution clean past a ring sector about the turn's
//!   axis: likewise, in the meridian half-plane (distance from the axis, height), each point
//!   projecting along the circle through it.
//!
//! On a solid of revolution the characteristic curve comes in closed form: at each place
//! round the tool's axis, the meridian's line pieces have one point where the velocity runs
//! along the surface (it is linear along a line), its circles two, opposite (the circle's
//! normal square to the velocity at its centre). Where the velocity there vanishes along the
//! meridian (a wheel set along the helix it grinds), the whole meridian there is
//! characteristic and is a candidate too.
//!
//! The section is then the blank's section less what those candidate curves bound, each
//! piece kept where the swept set does not cover it (how deep a point lies in the tool at
//! its deepest over the motion: sampled, then searched), and is swept back into a solid:
//! extruded, screwed or revolved. The curves are cubics through their points, as close as
//! about 1e-10 of the blank's size.

use std::f64::consts::{PI, TAU};

use super::section::{self, cross, norm, sub, Candidate, Cover, Loop, Outline, Piece, P2};
use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Circle, Curve};
use crate::geometry::surface::{Plane, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Axis, Direction, Frame, Point, Point2, Transform, Vector};
use crate::tolerance;
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// One rigid motion driven by a parameter t: a turn about `axis` by `turn`·t (radians,
/// right-handed about its direction) together with a slide along it by `advance`·t. A screw;
/// a turn alone where `advance` is zero; a slide alone where `turn` is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MotionStep {
    pub axis: Axis,
    pub turn: f64,
    pub advance: f64,
}

/// A rigid motion as t runs from `from` to `to`: its `steps` applied in order, each driven by
/// t, the first to the tool as placed and each next to where the one before put it (so a
/// tool turning in mesh with a wheel, seen from the wheel, is its own turn and then the
/// wheel's taken back).
#[derive(Debug, Clone, PartialEq)]
pub struct Motion {
    pub steps: Vec<MotionStep>,
    pub from: f64,
    pub to: f64,
}

impl Motion {
    /// Where the motion puts the tool at `t`.
    pub fn transform(&self, t: f64) -> Transform {
        let mut out = Transform::IDENTITY;
        for step in &self.steps {
            let s = Transform::translation(Vector::from(step.axis.direction) * (step.advance * t)) * Transform::rotation(step.axis, step.turn * t);
            out = s * out;
        }
        out
    }

    /// The motion's place at `t` as a rotation and translation.
    fn rigid(&self, t: f64) -> Rigid {
        let mut out = Rigid::IDENTITY;
        for step in &self.steps {
            out = Rigid::step(step, t).then(&out);
        }
        out
    }

    /// Where the tool's point `p` (as placed) is at `t`, and its velocity there.
    fn moving(&self, t: f64, p: V3) -> (V3, V3) {
        let (mut q, mut v) = (p, [0.0; 3]);
        for step in &self.steps {
            let s = Rigid::step(step, t);
            q = s.apply(q);
            let a = dir(step.axis.direction);
            let o = pt(step.axis.origin);
            let w = cross3(a, sub3(q, o));
            v = add3(s.rotate(v), add3(scale3(w, step.turn), scale3(a, step.advance)));
        }
        (q, v)
    }
}

/// A moving solid.
#[derive(Debug, Clone, PartialEq)]
pub enum Tool {
    /// A prism: the polygon `outline` in `frame`'s xy plane, from height `from` to `to`
    /// along its z axis.
    Prism { frame: Frame, outline: Vec<Point2>, from: f64, to: f64 },
    /// A solid of revolution about `frame`'s z axis: the region `meridian` bounds (closed, its
    /// pieces end to end, in (distance from the axis, height along it), at no negative
    /// distance) turned whole. A torus is one circle, a wheel with a V rim a triangle.
    Revolved { frame: Frame, meridian: Vec<MeridianPiece> },
}

/// A piece of a meridian, in (distance from the axis, height along it).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeridianPiece {
    Line { from: Point2, to: Point2 },
    /// The arc of the circle about `centre` of `radius` from angle `start` to `end` (radians
    /// from the distance axis towards the height axis; decreasing runs clockwise; a whole
    /// turn is a whole circle).
    Arc { centre: Point2, radius: f64, start: f64, end: f64 },
}

/// What a moving tool is cut from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Blank {
    /// The disc of `radius` about `frame`'s origin in its xy plane, from height 0 to
    /// `height` along its z axis: a wheel, a rod.
    Disc { frame: Frame, radius: f64, height: f64 },
    /// Of the ring about `frame`'s z axis between `radii`, between `heights` along it, what
    /// lies between polar `angles` (radians from its x axis towards y, the second the larger).
    RingSector { frame: Frame, radii: [f64; 2], heights: [f64; 2], angles: [f64; 2] },
}

/// `blank` less everything `tool` passes through as `motion` carries it: what at every
/// instant of the motion lies outside the tool.
///
/// Taken are the motions that make it a plane problem (see the module's notes): a prism
/// moved in the plane square to it, through all of a disc blank's height; a solid of
/// revolution screwed about the axis of a disc blank (a rod), starting wholly below it and
/// ending wholly above it (or the other way); a solid of revolution turned about the axis of
/// a ring sector, starting and ending clear of it. Others are refused, as are blanks and
/// tools that do not match their motion.
pub fn cut_by_motion(model: &mut Model, blank: &Blank, tool: &Tool, motion: &Motion) -> Result<Shape> {
    if !(motion.from < motion.to) || !motion.from.is_finite() || !motion.to.is_finite() {
        return Err(Error::InvalidInput(format!("a motion runs over an interval, not from {} to {}", motion.from, motion.to)));
    }
    match tool {
        Tool::Prism { frame, outline, from, to } => {
            let Blank::Disc { frame: blank_frame, radius, height } = *blank else {
                return Err(Error::InvalidInput("a prism's planar cut is taken from a disc blank".into()));
            };
            planar_cut(model, &blank_frame, radius, height, frame, outline, (*from, *to), motion)
        }
        Tool::Revolved { frame, meridian } => {
            let meridian = Meridian::new(meridian)?;
            let revolved = Revolved { centre: pt(frame.origin()), axis: dir(frame.z()), e1: dir(frame.x()), e2: dir(frame.y()), meridian };
            let [step] = motion.steps[..] else {
                return Err(Error::InvalidInput(
                    "a solid of revolution is cut along a single screw or turn: its swept set is then the same along it".into(),
                ));
            };
            if step.turn == 0.0 {
                return Err(Error::InvalidInput("a slide alone is not taken: only screws and turns".into()));
            }
            match *blank {
                Blank::Disc { frame: blank_frame, radius, height } if step.advance != 0.0 => {
                    screw_cut(model, &blank_frame, radius, height, &revolved, &step, motion)
                }
                Blank::RingSector { frame: blank_frame, radii, heights, angles } if step.advance == 0.0 => {
                    turn_cut(model, &blank_frame, radii, heights, angles, &revolved, &step, motion)
                }
                _ => Err(Error::InvalidInput("a screw cuts a disc blank (a rod), a turn a ring sector".into())),
            }
        }
    }
}

/// Sweeps `profile` along the screw about `axis` rising `lead` a turn (right-handed about the
/// axis's direction; a negative lead turns the other way), up `height` along the axis: each
/// point of it along its helix, as the quintic pieces matching the helix's point, tangent and
/// curvature at samples along it, within about 1e-10 of the profile's size. An edge sweeps the
/// surface through its poles' paths, so a B-spline edge's surface is as close; an edge whose
/// curve converts to a B-spline only with another parameter (a circle's rational arcs) is first
/// fitted with a cubic, as a pipe fits it, within 1e-7 of its size. A face square to the axis
/// sweeps a twisted prism, closed by the face and its screwed copy.
pub fn screw_extrude(model: &mut Model, profile: Shape, axis: Axis, lead: f64, height: f64) -> Result<super::Sweep> {
    use crate::math::Transformable;
    if !(lead != 0.0 && lead.is_finite() && height != 0.0 && height.is_finite()) {
        return Err(Error::InvalidInput(format!("a screw extrusion needs a lead and a height, got {lead} and {height}")));
    }
    {
        // (Told on the profile where it is defined, the axis brought there.)
        let local = Shape { placement: Transform::IDENTITY, ..profile };
        let local_axis = axis.transformed(&profile.placement.inverse());
        if let Some(fault) = super::profile_fault(model, local) {
            return Err(Error::InvalidInput(fault));
        }
        if super::crosses_axis(model, local, local_axis) {
            return Err(Error::InvalidInput("a screw's profile crossing its axis".into()));
        }
        if height.abs() > lead.abs() + crate::tolerance::LINEAR && turns_meet(model, local, local_axis, lead) {
            return Err(Error::InvalidInput(format!("a screw whose turns meet: its profile overlaps itself a lead of {lead} along the axis")));
        }
    }
    let k = lead / TAU;
    let angle = height / k;
    let along = Vector::from(axis.direction);
    let at = move |f: f64| Transform::translation(along * (height * f)) * Transform::rotation(axis, angle * f);
    // A point's velocity and acceleration (in the fraction) where the screw has carried it:
    // turning by `angle` about the axis and rising `height`.
    let rates = move |f: f64, p: Point| -> (Vector, Vector) {
        let q = p.transformed(&at(f));
        let r = (q - axis.origin) - along * (q - axis.origin).dot(along);
        let turning = along.cross(r);
        (turning * angle + along * height, along.cross(turning) * (angle * angle))
    };
    let bounds = crate::properties::bounding_box(model, profile);
    let probes = [bounds.min, bounds.max, Point::new(bounds.min.x, bounds.max.y, bounds.min.z), Point::new(bounds.max.x, bounds.min.y, bounds.max.z)];
    let size = (bounds.max - bounds.min).length() + height.abs();
    let mut samples = 4;
    let motion = loop {
        let params: Vec<f64> = (0..=samples).map(|i| i as f64 / samples as f64).collect();
        let transforms: Vec<Transform> = params.iter().map(|&f| at(f)).collect();
        let motion = super::pipe::PieceMotion { params, transforms, rates: Some(Box::new(rates)) };
        let mut error = 0.0f64;
        for &p in &probes {
            let path = super::pipe::path(&motion, p)?;
            for i in 0..4 * samples {
                let f = (i as f64 + 0.5) / (4 * samples) as f64;
                error = error.max(path.point(f).distance(p.transformed(&at(f))));
            }
        }
        if error <= 1e-10 * size || samples >= 1024 {
            break motion;
        }
        samples *= 2;
    };
    super::pipe::sweep_by_motions(model, profile, vec![motion], false)
}

/// Whether the turns of the screw of the planar face `profile` about `axis`, rising `lead` a
/// turn, meet: the axis's direction lying in the profile's plane, a point of the profile a
/// turn on lies where the profile moved `lead` along the axis does, so the turns meet where
/// the profile's region and its copy moved so overlap (more than along their boundaries).
/// They stay apart wherever the lead is at least the profile's extent along the axis. OCCT has
/// no such check (its sweeps make a solid passing through itself); this is the definition,
/// on the profile's loops sampled as polygons.
fn turns_meet(model: &Model, profile: Shape, axis: Axis, lead: f64) -> bool {
    if profile.kind() != NodeKind::Face {
        return false;
    }
    let Some(Surface::Plane(plane)) = model.surface(profile).cloned() else { return false };
    let (n, z) = (Vector::from(plane.normal()), Vector::from(axis.direction));
    if n.dot(z).abs() > tolerance::ANGULAR {
        return false;
    }
    let across = n.cross(z).normalized();
    let Some(loops) = crate::heal::face_loops(model, profile) else { return false };
    use crate::geometry::ParametricSurface;
    const SAMPLES: usize = 64;
    let polygons: Vec<Vec<Point2>> = loops
        .iter()
        .map(|l| {
            l.iter()
                .flat_map(|u| {
                    let n = if matches!(u.pcurve, crate::geometry::curve2::Curve2::Line(_)) { 1 } else { SAMPLES };
                    (0..n).map(move |k| u.at(k as f64 / n as f64))
                })
                .map(|uv| {
                    let p = plane.point(uv.x, uv.y) - axis.origin;
                    Point2::new(p.dot(across), p.dot(z))
                })
                .collect()
        })
        .collect();
    let (lo, hi) = polygons.iter().flatten().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p.y), hi.max(p.y)));
    if !(hi - lo > lead.abs() + tolerance::LINEAR) {
        return false;
    }
    let moved: Vec<Vec<Point2>> = polygons.iter().map(|p| p.iter().map(|q| Point2::new(q.x, q.y + lead)).collect()).collect();
    let near = tolerance::LINEAR;
    // Whether `q` lies inside the region the loops `region` bound, even-odd (`None` on them).
    let inside = |region: &[Vec<Point2>], q: Point2| -> Option<bool> {
        let mut crossings = 0;
        for p in region {
            for i in 0..p.len() {
                let (a, b) = (p[i], p[(i + 1) % p.len()]);
                let ab = b - a;
                let s = ((q - a).dot(ab) / ab.dot(ab).max(1e-300)).clamp(0.0, 1.0);
                if (q - (a + ab * s)).length() <= near {
                    return None;
                }
                if (a.y > q.y) != (b.y > q.y) && a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) > q.x {
                    crossings += 1;
                }
            }
        }
        Some(crossings % 2 == 1)
    };
    let sides = |region: &[Vec<Point2>]| -> Vec<(Point2, Point2)> { region.iter().flat_map(|p| (0..p.len()).map(move |i| (p[i], p[(i + 1) % p.len()]))).collect() };
    let (ours, theirs) = (sides(&polygons), sides(&moved));
    // (Two sides crossing each other properly, each's ends clear of the other's line.)
    let side_of = |a: Point2, b: Point2, c: Point2| {
        let ab = b - a;
        (ab.x * (c.y - a.y) - ab.y * (c.x - a.x)) / ab.length().max(1e-300)
    };
    let crossing = ours.iter().any(|&(a, b)| {
        theirs.iter().any(|&(c, d)| {
            let (d1, d2, d3, d4) = (side_of(c, d, a), side_of(c, d, b), side_of(a, b, c), side_of(a, b, d));
            ((d1 > near && d2 < -near) || (d1 < -near && d2 > near)) && ((d3 > near && d4 < -near) || (d3 < -near && d4 > near))
        })
    });
    // (Else one region within the other in part: a corner or a side's middle of one strictly
    // inside the other, as where a side of one lies along a side of the other.)
    let within = |sides: &[(Point2, Point2)], region: &[Vec<Point2>]| sides.iter().any(|&(a, b)| [a, a + (b - a) * 0.5].into_iter().any(|q| inside(region, q) == Some(true)));
    crossing || within(&theirs, &polygons) || within(&ours, &moved)
}

// ---- Small vector arithmetic, for the many evaluations the sections take.

type V3 = [f64; 3];

fn pt(p: Point) -> V3 {
    [p.x, p.y, p.z]
}

fn dir(d: Direction) -> V3 {
    [d.x(), d.y(), d.z()]
}

fn add3(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot3(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn length3(a: V3) -> f64 {
    dot3(a, a).sqrt()
}

/// A rotation `m` (rows) and translation `d`: x ↦ m x + d.
#[derive(Clone, Copy)]
struct Rigid {
    m: [[f64; 3]; 3],
    d: V3,
}

impl Rigid {
    const IDENTITY: Rigid = Rigid { m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], d: [0.0; 3] };

    /// A turn by `angle` about the axis through `o` along unit `a`.
    fn turn(o: V3, a: V3, angle: f64) -> Rigid {
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        let [x, y, z] = a;
        let m = [
            [c + x * x * t, x * y * t - z * s, x * z * t + y * s],
            [y * x * t + z * s, c + y * y * t, y * z * t - x * s],
            [z * x * t - y * s, z * y * t + x * s, c + z * z * t],
        ];
        let r = Rigid { m, d: [0.0; 3] };
        Rigid { m, d: sub3(o, r.rotate(o)) }
    }

    fn step(step: &MotionStep, t: f64) -> Rigid {
        let a = dir(step.axis.direction);
        let mut r = Rigid::turn(pt(step.axis.origin), a, step.turn * t);
        r.d = add3(r.d, scale3(a, step.advance * t));
        r
    }

    fn rotate(&self, v: V3) -> V3 {
        let m = &self.m;
        [dot3(m[0], v), dot3(m[1], v), dot3(m[2], v)]
    }

    fn apply(&self, p: V3) -> V3 {
        add3(self.rotate(p), self.d)
    }

    /// The inverse applied: m⁻¹ (p − d).
    fn unapply(&self, p: V3) -> V3 {
        let q = sub3(p, self.d);
        let m = &self.m;
        [m[0][0] * q[0] + m[1][0] * q[1] + m[2][0] * q[2], m[0][1] * q[0] + m[1][1] * q[1] + m[2][1] * q[2], m[0][2] * q[0] + m[1][2] * q[1] + m[2][2] * q[2]]
    }

    /// `self` after `first`.
    fn then(&self, first: &Rigid) -> Rigid {
        let mut m = [[0.0; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = (0..3).map(|k| self.m[i][k] * first.m[k][j]).sum();
            }
        }
        Rigid { m, d: self.apply(first.d) }
    }
}

/// A frame's coordinates: to and from the world's.
#[derive(Clone, Copy)]
struct Coordinates {
    origin: V3,
    x: V3,
    y: V3,
    z: V3,
}

impl Coordinates {
    fn of(frame: &Frame) -> Coordinates {
        Coordinates { origin: pt(frame.origin()), x: dir(frame.x()), y: dir(frame.y()), z: dir(frame.z()) }
    }

    fn local(&self, p: V3) -> V3 {
        let q = sub3(p, self.origin);
        [dot3(q, self.x), dot3(q, self.y), dot3(q, self.z)]
    }

    fn world(&self, l: V3) -> V3 {
        add3(self.origin, add3(scale3(self.x, l[0]), add3(scale3(self.y, l[1]), scale3(self.z, l[2]))))
    }
}

/// The signed distance from `p` to the closed polygon `outline` (positive inside).
fn polygon_depth(outline: &[P2], p: P2) -> f64 {
    let n = outline.len();
    let mut nearest = f64::INFINITY;
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (outline[i], outline[(i + 1) % n]);
        let ab = sub(b, a);
        let f = ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / (ab[0] * ab[0] + ab[1] * ab[1]);
        let f = f.clamp(0.0, 1.0);
        nearest = nearest.min(norm(sub(p, [a[0] + ab[0] * f, a[1] + ab[1] * f])));
        if (a[1] > p[1]) != (b[1] > p[1]) && a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) > p[0] {
            inside = !inside;
        }
    }
    if inside { nearest } else { -nearest }
}

// ---- Planar motions of a prism.

/// The planar cut: the disc's section less the region the prism's outline sweeps,
/// extruded through the disc's height.
#[allow(clippy::too_many_arguments)]
fn planar_cut(model: &mut Model, blank_frame: &Frame, radius: f64, height: f64, tool_frame: &Frame, outline: &[Point2], depth: (f64, f64), motion: &Motion) -> Result<Shape> {
    let plane = Coordinates::of(blank_frame);
    let normal = plane.z;
    let size = 2.0 * radius;
    if !(radius > 0.0 && height > 0.0) {
        return Err(Error::InvalidInput("a disc blank needs a radius and a height".into()));
    }
    if length3(cross3(dir(tool_frame.z()), normal)) > 1e-12 {
        return Err(Error::InvalidInput("a prism's planar cut needs the prism square to the blank's plane".into()));
    }
    // (The prism's heights in the blank's: they must hold the blank's.)
    let tool_z = dot3(sub3(pt(tool_frame.origin()), plane.origin), normal);
    let sense = dot3(dir(tool_frame.z()), normal);
    let (low, high) = ((tool_z + sense * depth.0).min(tool_z + sense * depth.1), (tool_z + sense * depth.0).max(tool_z + sense * depth.1));
    if low > 1e-12 * size || high < height - 1e-12 * size {
        return Err(Error::InvalidInput("a prism's planar cut needs the prism through all of the blank's height".into()));
    }
    for step in &motion.steps {
        let a = dir(step.axis.direction);
        let planar = if step.turn != 0.0 {
            length3(cross3(a, normal)) <= 1e-12 && step.advance == 0.0
        } else {
            dot3(a, normal).abs() <= 1e-12
        };
        if !planar {
            return Err(Error::InvalidInput("a prism's planar cut needs a motion in the blank's plane: turns square to it, slides along it".into()));
        }
    }
    // The outline in the blank's plane, counter-clockwise.
    let tool = Coordinates::of(tool_frame);
    let mut body: Vec<P2> = outline
        .iter()
        .map(|p| {
            let l = plane.local(tool.world([p.x, p.y, 0.0]));
            [l[0], l[1]]
        })
        .collect();
    if body.len() < 3 {
        return Err(Error::InvalidInput("a prism's outline needs three corners or more".into()));
    }
    if section::signed_area(&body) < 0.0 {
        body.reverse();
    }
    let planar = Planar { motion, plane, body };
    let (candidates, slope) = planar.candidates(size)?;
    let tool_size = planar.body.iter().map(|&p| norm(sub(p, planar.body[0]))).fold(0.0, f64::max);
    let cover = PlanarCover { planar: &planar, step: 0.02 * tool_size / slope, slope, enough: 1e-11 * size };
    let loops = section::region(&candidates, &cover, Outline::Disc { centre: [0.0, 0.0], radius })?;
    let map = PlaneMap { origin: plane.origin, x: plane.x, y: plane.y, splines: false };
    let faces = section_faces(model, &candidates, &loops, &map, size)?;
    Ok(super::extrude(model, faces, Vector::new(normal[0], normal[1], normal[2]) * height)?.shape)
}

struct Planar<'a> {
    motion: &'a Motion,
    plane: Coordinates,
    /// The outline as placed, in the plane's coordinates, counter-clockwise.
    body: Vec<P2>,
}

impl Planar<'_> {
    fn lift(&self, p: P2) -> V3 {
        self.plane.world([p[0], p[1], 0.0])
    }

    fn flat(&self, q: V3) -> P2 {
        let l = self.plane.local(q);
        [l[0], l[1]]
    }

    /// Where the outline's point `p` is at `t`, and its velocity there (in the plane).
    fn moving(&self, t: f64, p: P2) -> (P2, P2) {
        let (q, v) = self.motion.moving(t, self.lift(p));
        let w = [dot3(v, self.plane.x), dot3(v, self.plane.y)];
        (self.flat(q), w)
    }

    /// The candidates: the paths of the convex corners, the envelopes of the sides, the
    /// sides where the motion starts and ends. With the fastest any corner moves.
    fn candidates(&self, size: f64) -> Result<(Vec<Candidate<'_>>, f64)> {
        let (t0, t1) = (self.motion.from, self.motion.to);
        let n = self.body.len();
        let mut slope = 0.0f64;
        for i in 0..=1000 {
            let t = t0 + (t1 - t0) * i as f64 / 1000.0;
            for &p in &self.body {
                slope = slope.max(norm(self.moving(t, p).1));
            }
        }
        // (Points just off the outline move a little faster than its corners.)
        let slope = 1.25 * slope + 1e-9;
        let convex: Vec<bool> = (0..n)
            .map(|i| cross(sub(self.body[i], self.body[(i + n - 1) % n]), sub(self.body[(i + 1) % n], self.body[i])) > 0.0)
            .collect();
        let mut corner_breaks: Vec<Vec<f64>> = vec![Vec::new(); n];
        let mut end_breaks: [Vec<Vec<f64>>; 2] = [vec![Vec::new(); n], vec![Vec::new(); n]];
        let mut candidates: Vec<Candidate> = Vec::new();
        // The sides' envelopes: on side i (a to b), the fraction μ where the velocity runs
        // along it, the velocity being affine along the side.
        for i in 0..n {
            let (a, b) = (self.body[i], self.body[(i + 1) % n]);
            let fraction = move |t: f64| -> Option<(f64, P2)> {
                let ((pa, va), (pb, vb)) = (self.moving(t, a), self.moving(t, b));
                let d = sub(pb, pa);
                let denominator = cross(sub(va, vb), d);
                if denominator == 0.0 {
                    return None;
                }
                let mu = cross(va, d) / denominator;
                Some((mu, [pa[0] + d[0] * mu, pa[1] + d[1] * mu]))
            };
            let valid = |t: f64| fraction(t).is_some_and(|(mu, _)| (0.0..=1.0).contains(&mu));
            for (from, to) in runs(&valid, t0, t1, 4096) {
                for s in [from, to] {
                    let Some((mu, _)) = fraction(s) else { continue };
                    if s == t0 || s == t1 {
                        // (Starting or ending with the motion: on the side where it starts or ends.)
                        let which = if s == t0 { 0 } else { 1 };
                        end_breaks[which][i].push(mu);
                    } else {
                        let corner = if mu < 0.5 { i } else { (i + 1) % n };
                        corner_breaks[corner].push(s);
                    }
                }
                candidates.push(Candidate {
                    point: Box::new(move |t: f64| fraction(t).map(|(_, p)| p).unwrap_or_else(|| self.moving(t, a).0)),
                    range: (from, to),
                    closed: false,
                    breaks: Vec::new(),
                });
            }
        }
        // The convex corners' paths.
        let closed = |p: P2| norm(sub(self.moving(t0, p).0, self.moving(t1, p).0)) <= 1e-12 * size;
        for i in 0..n {
            if !convex[i] {
                continue;
            }
            let p = self.body[i];
            candidates.push(Candidate {
                point: Box::new(move |t: f64| self.moving(t, p).0),
                range: (t0, t1),
                closed: closed(p),
                breaks: std::mem::take(&mut corner_breaks[i]),
            });
        }
        // The sides where the motion starts and where it ends.
        for (which, t) in [(0, t0), (1, t1)] {
            for i in 0..n {
                let (a, b) = (self.moving(t, self.body[i]).0, self.moving(t, self.body[(i + 1) % n]).0);
                candidates.push(Candidate {
                    point: Box::new(move |f: f64| [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]),
                    range: (0.0, 1.0),
                    closed: false,
                    breaks: std::mem::take(&mut end_breaks[which][i]),
                });
            }
        }
        Ok((candidates, slope))
    }
}

struct PlanarCover<'a> {
    planar: &'a Planar<'a>,
    step: f64,
    slope: f64,
    enough: f64,
}

impl Cover for PlanarCover<'_> {
    fn depth(&self, p: P2) -> f64 {
        let q = self.planar.lift(p);
        let motion = self.planar.motion;
        let at = |t: f64| polygon_depth(&self.planar.body, self.planar.flat(motion.rigid(t).unapply(q)));
        section::largest(&at, motion.from, motion.to, self.step, self.slope, false, self.enough)
    }
}

/// The maximal intervals of `[a, b]` where `valid` holds: sampled `n` times, each change
/// found by bisection.
fn runs(valid: &dyn Fn(f64) -> bool, a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    let at = |i: usize| a + (b - a) * i as f64 / n as f64;
    let states: Vec<bool> = (0..=n).map(|i| valid(at(i))).collect();
    let change = |lo: f64, hi: f64, at_lo: bool| {
        let (mut lo, mut hi) = (lo, hi);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if valid(mid) == at_lo {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        // (The side where it holds.)
        if at_lo { lo } else { hi }
    };
    let mut out = Vec::new();
    let mut open = if states[0] { Some(a) } else { None };
    for i in 0..n {
        if states[i] != states[i + 1] {
            let s = change(at(i), at(i + 1), states[i]);
            if states[i + 1] {
                open = Some(s);
            } else if let Some(from) = open.take() {
                out.push((from, s));
            }
        }
    }
    if let Some(from) = open {
        out.push((from, b));
    }
    out.retain(|&(x, y)| y > x);
    out
}

// ---- Solids of revolution, screwed or turned.

/// A meridian's pieces, counter-clockwise.
#[derive(Clone)]
struct Meridian {
    pieces: Vec<MeridianPiece>,
}

impl Meridian {
    fn new(pieces: &[MeridianPiece]) -> Result<Meridian> {
        if pieces.is_empty() {
            return Err(Error::InvalidInput("a meridian needs a piece".into()));
        }
        let mut m = Meridian { pieces: pieces.to_vec() };
        let n = m.pieces.len();
        let scale = m.points().iter().map(|p| norm(*p)).fold(1.0, f64::max);
        for i in 0..n {
            let (end, start) = (m.point(i, 1.0), m.point((i + 1) % n, 0.0));
            if norm(sub(end, start)) > 1e-9 * scale {
                return Err(Error::InvalidInput("a meridian's pieces must run end to end, closed".into()));
            }
        }
        if m.points().iter().any(|p| p[0] < -1e-12 * scale) {
            return Err(Error::InvalidInput("a meridian lies at no negative distance from its axis".into()));
        }
        if let Some(MeridianPiece::Arc { radius, .. }) = m.pieces.iter().find(|p| matches!(p, MeridianPiece::Arc { radius, .. } if !(*radius > 0.0))) {
            return Err(Error::InvalidInput(format!("a meridian's arc needs a radius, got {radius}")));
        }
        if section::signed_area(&m.points()) < 0.0 {
            m.pieces.reverse();
            for p in &mut m.pieces {
                *p = match *p {
                    MeridianPiece::Line { from, to } => MeridianPiece::Line { from: to, to: from },
                    MeridianPiece::Arc { centre, radius, start, end } => MeridianPiece::Arc { centre, radius, start: end, end: start },
                };
            }
        }
        Ok(m)
    }

    fn points(&self) -> Vec<P2> {
        (0..self.pieces.len()).flat_map(|i| (0..32).map(move |k| (i, k as f64 / 32.0))).map(|(i, f)| self.point(i, f)).collect()
    }

    /// Piece `i`'s point at fraction `f`.
    fn point(&self, i: usize, f: f64) -> P2 {
        match self.pieces[i] {
            MeridianPiece::Line { from, to } => [from.x + (to.x - from.x) * f, from.y + (to.y - from.y) * f],
            MeridianPiece::Arc { centre, radius, start, end } => {
                let a = start + (end - start) * f;
                [centre.x + radius * a.cos(), centre.y + radius * a.sin()]
            }
        }
    }

    /// Piece `i`'s unit tangent at fraction `f`, as it runs.
    fn tangent(&self, i: usize, f: f64) -> P2 {
        match self.pieces[i] {
            MeridianPiece::Line { from, to } => {
                let d = [to.x - from.x, to.y - from.y];
                let l = norm(d);
                [d[0] / l, d[1] / l]
            }
            MeridianPiece::Arc { start, end, .. } => {
                let a = start + (end - start) * f;
                let s = (end - start).signum();
                [-a.sin() * s, a.cos() * s]
            }
        }
    }

    /// The signed distance from (distance, height) `p` to the meridian (positive inside).
    fn depth(&self, p: P2) -> f64 {
        let mut nearest = f64::INFINITY;
        let mut inside = false;
        for piece in &self.pieces {
            match *piece {
                MeridianPiece::Line { from, to } => {
                    let (a, b) = ([from.x, from.y], [to.x, to.y]);
                    let ab = sub(b, a);
                    let f = (((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / (ab[0] * ab[0] + ab[1] * ab[1])).clamp(0.0, 1.0);
                    nearest = nearest.min(norm(sub(p, [a[0] + ab[0] * f, a[1] + ab[1] * f])));
                    if (a[1] > p[1]) != (b[1] > p[1]) && a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) > p[0] {
                        inside = !inside;
                    }
                }
                MeridianPiece::Arc { centre, radius, start, end } => {
                    let c = [centre.x, centre.y];
                    let within = |angle: f64| {
                        let (lo, hi) = (start.min(end), start.max(end));
                        if hi - lo >= TAU - 1e-12 {
                            return true;
                        }
                        let mut a = angle;
                        while a < lo {
                            a += TAU;
                        }
                        while a > lo + TAU {
                            a -= TAU;
                        }
                        a <= hi
                    };
                    let d = sub(p, c);
                    let angle = d[1].atan2(d[0]);
                    let to_arc = if within(angle) {
                        (norm(d) - radius).abs()
                    } else {
                        let (a, b) = ([c[0] + radius * start.cos(), c[1] + radius * start.sin()], [c[0] + radius * end.cos(), c[1] + radius * end.sin()]);
                        norm(sub(p, a)).min(norm(sub(p, b)))
                    };
                    nearest = nearest.min(to_arc);
                    // (The ray to the right crosses the circle where x = cx ± √(r² − dy²), the
                    // crossing counted where the arc holds it and the ray's point does not lie
                    // on one of the arc's ends' heights' sides unevenly: half-open in the height.)
                    let dy = p[1] - c[1];
                    if dy.abs() < radius {
                        let w = (radius * radius - dy * dy).sqrt();
                        for x in [c[0] - w, c[0] + w] {
                            if x > p[0] {
                                let angle = dy.atan2(x - c[0]);
                                if within(angle) {
                                    // (Half-open: the circle's top and bottom count once.)
                                    inside = !inside;
                                }
                            }
                        }
                    }
                }
            }
        }
        if inside { nearest } else { -nearest }
    }

    /// The largest of ζ·along + ρ·across over the meridian (`across` ≥ 0): how far the solid
    /// reaches along a direction making those parts with its axis and square to it.
    fn reach(&self, along: f64, across: f64) -> f64 {
        let mut best = f64::NEG_INFINITY;
        for (i, piece) in self.pieces.iter().enumerate() {
            for f in [0.0, 1.0] {
                let p = self.point(i, f);
                best = best.max(p[1] * along + p[0] * across);
            }
            if let MeridianPiece::Arc { centre, radius, start, end } = *piece {
                let angle = along.atan2(across);
                let (lo, hi) = (start.min(end), start.max(end));
                let mut a = angle;
                while a < lo {
                    a += TAU;
                }
                if a <= hi || hi - lo >= TAU - 1e-12 {
                    best = best.max(centre.y * along + centre.x * across + radius * along.hypot(across));
                }
            }
        }
        best
    }
}

/// A solid of revolution placed: its centre, axis and the frame's other two axes.
struct Revolved {
    centre: V3,
    axis: V3,
    e1: V3,
    e2: V3,
    meridian: Meridian,
}

impl Revolved {
    fn radial(&self, u: f64) -> V3 {
        let (s, c) = u.sin_cos();
        add3(scale3(self.e1, c), scale3(self.e2, s))
    }

    /// The point at angle `u` round the axis and (distance, height) `m`.
    fn at(&self, u: f64, m: P2) -> V3 {
        add3(self.centre, add3(scale3(self.radial(u), m[0]), scale3(self.axis, m[1])))
    }

    /// The normal at angle `u` of the meridian's outward normal `n`.
    fn normal(&self, u: f64, n: P2) -> V3 {
        add3(scale3(self.radial(u), n[0]), scale3(self.axis, n[1]))
    }

    /// (Distance, height) of world point `q`.
    fn meridian_of(&self, q: V3) -> P2 {
        let d = sub3(q, self.centre);
        let h = dot3(d, self.axis);
        [length3(sub3(d, scale3(self.axis, h))), h]
    }

    fn depth(&self, q: V3) -> f64 {
        self.meridian.depth(self.meridian_of(q))
    }

    /// How far the solid reaches along unit `d`, from its centre.
    fn reach(&self, d: V3) -> f64 {
        let along = dot3(d, self.axis);
        let across = length3(sub3(d, scale3(self.axis, along)));
        self.meridian.reach(along, across)
    }
}

/// A screw's or turn's velocity field, per radian of the turn: ω × (x − o) + k ω.
#[derive(Clone, Copy)]
struct Field {
    origin: V3,
    axis: V3,
    k: f64,
}

impl Field {
    fn at(&self, x: V3) -> V3 {
        add3(cross3(self.axis, sub3(x, self.origin)), scale3(self.axis, self.k))
    }
}

/// How a screw or a turn projects a point onto the section: along the helix through it onto
/// the blank's plane (`Screw`), or along the circle through it onto the meridian half-plane
/// (`Turn`), in the blank's coordinates.
#[derive(Clone, Copy)]
enum Quotient {
    Screw { blank: Coordinates, k: f64 },
    Turn { blank: Coordinates },
}

impl Quotient {
    fn project(&self, q: V3) -> P2 {
        match *self {
            Quotient::Screw { blank, k } => {
                let [x, y, z] = blank.local(q);
                let (s, c) = (-z / k).sin_cos();
                [x * c - y * s, x * s + y * c]
            }
            Quotient::Turn { blank } => {
                let [x, y, z] = blank.local(q);
                [x.hypot(y), z]
            }
        }
    }

    /// The section's point `p` moved by the screw or turn by `s`.
    fn lift(&self, p: P2, s: f64) -> V3 {
        match *self {
            Quotient::Screw { blank, k } => {
                let (sn, c) = s.sin_cos();
                blank.world([p[0] * c - p[1] * sn, p[0] * sn + p[1] * c, k * s])
            }
            Quotient::Turn { blank } => {
                let (sn, c) = s.sin_cos();
                blank.world([p[0] * c, p[0] * sn, p[1]])
            }
        }
    }

    /// How fast the section's point `p` moves as it is screwed or turned.
    fn speed(&self, p: P2) -> f64 {
        match *self {
            Quotient::Screw { k, .. } => norm(p).hypot(k),
            Quotient::Turn { .. } => p[0].abs(),
        }
    }
}

struct QuotientCover<'a> {
    tool: &'a Revolved,
    quotient: Quotient,
    /// The turns to search (for a screw, those carrying the section's plane through the tool).
    range: (f64, f64),
    periodic: bool,
    step_length: f64,
    enough: f64,
}

impl Cover for QuotientCover<'_> {
    fn depth(&self, p: P2) -> f64 {
        let speed = self.quotient.speed(p).max(1e-300);
        let at = |s: f64| self.tool.depth(self.quotient.lift(p, s));
        section::largest(&at, self.range.0, self.range.1, self.step_length / speed, speed, self.periodic, self.enough)
    }
}

/// A first-order trigonometric polynomial in the angle u round the tool's axis:
/// `a + c cos u + s sin u`, what the characteristic condition's terms are (the velocity
/// field being affine and the tool's points and normals turning with u).
#[derive(Clone, Copy)]
struct Trig {
    a: f64,
    c: f64,
    s: f64,
}

impl Trig {
    /// The polynomial `f` is, from its values at 0, π/2 and π.
    fn of(f: &dyn Fn(f64) -> f64) -> Trig {
        let (f0, f1, f2) = (f(0.0), f(0.5 * PI), f(PI));
        let a = 0.5 * (f0 + f2);
        Trig { a, c: f0 - a, s: f1 - a }
    }

    fn at(&self, u: f64) -> f64 {
        self.a + self.c * u.cos() + self.s * u.sin()
    }

    /// Its value at u over 2 sin((u − d)/2), where it vanishes at d: without the cancellation
    /// near d.
    fn reduced(&self, u: f64, d: f64) -> f64 {
        let m = 0.5 * (u + d);
        -self.c * m.sin() + self.s * m.cos()
    }

    fn magnitude(&self) -> f64 {
        self.a.abs() + self.c.hypot(self.s)
    }
}

/// How near a degenerate angle d (in sin((u − d)/2)) the reduced forms take over.
const NEAR: f64 = 1e-4;

/// The degenerate angle nearest `u` if `u` lies near one: with sin((u − d)/2).
fn near(flips: &[f64], u: f64) -> Option<(f64, f64)> {
    flips.iter().map(|&d| (d, (0.5 * (u - d)).sin())).find(|&(_, h)| h.abs() < NEAR)
}

/// The characteristic curves of a solid of revolution under `field`, projected by
/// `quotient`, with its rims: the section's candidates.
fn revolved_candidates<'a>(tool: &'a Revolved, field: Field, quotient: Quotient) -> Vec<Candidate<'a>> {
    let meridian = &tool.meridian;
    let n = meridian.pieces.len();
    // The vertices that are rims: convex corners off the axis.
    let scale = meridian.points().iter().map(|p| norm(*p)).fold(1.0, f64::max);
    let rim: Vec<bool> = (0..n)
        .map(|i| {
            let (t_in, t_out) = (meridian.tangent((i + n - 1) % n, 1.0), meridian.tangent(i, 0.0));
            cross(t_in, t_out) > 1e-9 && meridian.point(i, 0.0)[0] > 1e-9 * scale
        })
        .collect();
    let mut rim_breaks: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut candidates: Vec<Candidate<'a>> = Vec::new();
    let wrap = |u: f64| u.rem_euclid(TAU);
    for i in 0..n {
        match meridian.pieces[i] {
            MeridianPiece::Line { from, to } => {
                let (p0, p1) = ([from.x, from.y], [to.x, to.y]);
                // (A piece on the axis sweeps nothing.)
                if p0[0] <= 1e-9 * scale && p1[0] <= 1e-9 * scale {
                    continue;
                }
                let d = sub(p1, p0);
                let l = norm(d);
                let nu = [d[1] / l, -d[0] / l];
                // The condition at the piece's two ends, at angle u: linear between.
                let condition = |p: P2| Trig::of(&|u: f64| dot3(field.at(tool.at(u, p)), tool.normal(u, nu)));
                let (f0, f1) = (condition(p0), condition(p1));
                let flips = degenerate_angles(&|u: f64| f0.at(u).hypot(f1.at(u)), f0.magnitude() + f1.magnitude());
                // σ along the piece where the condition vanishes (its limit at a degenerate
                // angle, where all of the line does).
                let fraction = {
                    let flips = flips.clone();
                    move |u: f64| -> f64 {
                        let (g0, g1) = match near(&flips, u) {
                            Some((dd, _)) => (f0.reduced(u, dd), f1.reduced(u, dd)),
                            None => (f0.at(u), f1.at(u)),
                        };
                        g0 / (g0 - g1)
                    }
                };
                let point = {
                    let fraction = fraction.clone();
                    move |u: f64| {
                        let s = fraction(u);
                        quotient.project(tool.at(u, [p0[0] + d[0] * s, p0[1] + d[1] * s]))
                    }
                };
                let valid = |u: f64| (0.0..=1.0).contains(&fraction(u));
                for (from_u, to_u, closed) in periodic_runs(&valid, TAU) {
                    if !closed {
                        for u in [from_u, to_u] {
                            let vertex = if fraction(u) < 0.5 { i } else { (i + 1) % n };
                            rim_breaks[vertex].push(wrap(u));
                        }
                    }
                    candidates.push(Candidate { point: Box::new(point.clone()), range: (from_u, to_u), closed, breaks: Vec::new() });
                }
                // Where the whole line is characteristic.
                for &dd in &flips {
                    let s = fraction(dd);
                    rim_breaks[i].push(dd);
                    rim_breaks[(i + 1) % n].push(dd);
                    candidates.push(Candidate {
                        point: Box::new(move |f: f64| quotient.project(tool.at(dd, [p0[0] + d[0] * f, p0[1] + d[1] * f]))),
                        range: (0.0, 1.0),
                        closed: false,
                        breaks: if s > 0.0 && s < 1.0 { vec![s] } else { Vec::new() },
                    });
                }
            }
            MeridianPiece::Arc { centre, radius, start, end } => {
                let m = [centre.x, centre.y];
                let full = (end - start).abs() >= TAU - 1e-12;
                // (A clockwise arc's outward normal points at its centre.)
                let r = if end < start { -radius } else { radius };
                // At angle u, the circle's normal (cos v, sin v) square to the velocity at its
                // centre: A cos v + B sin v = 0, the normal ±(B, −A).
                let a_term = Trig::of(&|u: f64| dot3(field.at(tool.at(u, m)), tool.radial(u)));
                let b_term = Trig::of(&|u: f64| dot3(field.at(tool.at(u, m)), tool.axis));
                let flips = degenerate_angles(&|u: f64| a_term.at(u).hypot(b_term.at(u)), a_term.magnitude() + b_term.magnitude());
                // The normal of branch `sign` at u (u counted on round a second turn where the
                // branches join), its sign flipping at each degenerate angle passed so that it
                // runs on smoothly through it.
                let branch = {
                    let flips = flips.clone();
                    move |u: f64, sign: f64| -> P2 {
                        let w = u.rem_euclid(2.0 * TAU);
                        let (turns, at) = if w >= TAU { (1, w - TAU) } else { (0, w) };
                        let count = flips.iter().filter(|&&dd| dd < at).count() + turns * flips.len();
                        let sign = if count % 2 == 0 { sign } else { -sign };
                        let raw = match near(&flips, at) {
                            // (Near d, (B, −A) is 2 sin((u − d)/2) times the reduced pair; on d's
                            // left the flip is still to come.)
                            Some((dd, h)) => {
                                let s = if h < 0.0 { -1.0 } else { 1.0 };
                                [s * b_term.reduced(at, dd), -s * a_term.reduced(at, dd)]
                            }
                            None => [b_term.at(at), -a_term.at(at)],
                        };
                        let l = norm(raw);
                        [sign * raw[0] / l, sign * raw[1] / l]
                    }
                };
                // The meridian angle of the circle's point with outward normal `nrm`.
                let angle_of = move |nrm: P2| {
                    let a = (r.signum() * nrm[1]).atan2(r.signum() * nrm[0]);
                    let lo = start.min(end);
                    (a - lo).rem_euclid(TAU) + lo
                };
                let within = move |nrm: P2| full || angle_of(nrm) <= start.max(end);
                let point = {
                    let branch = branch.clone();
                    move |u: f64, sign: f64| {
                        let nrm = branch(u, sign);
                        quotient.project(tool.at(u.rem_euclid(TAU), [m[0] + r * nrm[0], m[1] + r * nrm[1]]))
                    }
                };
                // An odd number of flips joins the two branches into one curve twice round.
                let (signs, period) = if flips.len() % 2 == 1 { (vec![1.0], 2.0 * TAU) } else { (vec![1.0, -1.0], TAU) };
                for sign in signs {
                    let valid = |u: f64| within(branch(u, sign));
                    for (from_u, to_u, closed) in periodic_runs(&valid, period) {
                        if !closed {
                            // (Ending at the arc's start or end: a vertex of the meridian.)
                            for u in [from_u, to_u] {
                                let a = angle_of(branch(u, sign));
                                let at_start = (a - start).abs() < (a - end).abs();
                                let vertex = if at_start == (end > start) { i } else { (i + 1) % n };
                                rim_breaks[vertex].push(wrap(u));
                            }
                        }
                        let point = point.clone();
                        candidates.push(Candidate { point: Box::new(move |u: f64| point(u, sign)), range: (from_u, to_u), closed, breaks: Vec::new() });
                    }
                }
                // Where the whole circle is characteristic: there the branches pass through it.
                for &dd in &flips {
                    let limit = [b_term.reduced(dd, dd), -a_term.reduced(dd, dd)];
                    let l = norm(limit);
                    let breaks: Vec<f64> = [1.0, -1.0].iter().map(|&s| angle_of([s * limit[0] / l, s * limit[1] / l])).collect();
                    candidates.push(Candidate {
                        point: Box::new(move |a: f64| quotient.project(tool.at(dd, [m[0] + radius * a.cos(), m[1] + radius * a.sin()]))),
                        range: (start.min(end), start.max(end)),
                        closed: full,
                        breaks,
                    });
                    if !full {
                        rim_breaks[i].push(dd);
                        rim_breaks[(i + 1) % n].push(dd);
                    }
                }
            }
        }
    }
    // The rims: each convex vertex's circle.
    for i in 0..n {
        if !rim[i] {
            continue;
        }
        let p = meridian.point(i, 0.0);
        candidates.push(Candidate {
            point: Box::new(move |u: f64| quotient.project(tool.at(u, p))),
            range: (0.0, TAU),
            closed: true,
            breaks: std::mem::take(&mut rim_breaks[i]),
        });
    }
    candidates
}

/// The angles in [0, 2π) where `magnitude` (a characteristic condition's size, at most
/// about `scale`) vanishes: sampled, each local least searched, those within rounding of
/// zero kept.
fn degenerate_angles(magnitude: &dyn Fn(f64) -> f64, scale: f64) -> Vec<f64> {
    let n = 4096;
    let values: Vec<f64> = (0..n).map(|i| magnitude(TAU * i as f64 / n as f64)).collect();
    let mut out = Vec::new();
    for i in 0..n {
        let (before, v, after) = (values[(i + n - 1) % n], values[i], values[(i + 1) % n]);
        if v > before || v > after || v > 1e-2 * scale {
            continue;
        }
        let h = TAU / n as f64;
        // (Golden section for the least.)
        let (mut lo, mut hi) = (h * i as f64 - h, h * i as f64 + h);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        let (mut x1, mut x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        let (mut f1, mut f2) = (magnitude(x1), magnitude(x2));
        for _ in 0..100 {
            if f1 > f2 {
                lo = x1;
                x1 = x2;
                f1 = f2;
                x2 = lo + g * (hi - lo);
                f2 = magnitude(x2);
            } else {
                hi = x2;
                x2 = x1;
                f2 = f1;
                x1 = hi - g * (hi - lo);
                f1 = magnitude(x1);
            }
        }
        let u = 0.5 * (lo + hi);
        if magnitude(u) <= 1e-12 * scale {
            let u = u.rem_euclid(TAU);
            if !out.iter().any(|&w: &f64| (w - u).abs() < 1e-9) {
                out.push(u);
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

/// The runs of [0, `period`) where `valid` holds, as (from, to, closed): a run round the
/// whole period is closed; the others start after a place where it fails, so none is split
/// at the period's ends.
fn periodic_runs(valid: &dyn Fn(f64) -> bool, period: f64) -> Vec<(f64, f64, bool)> {
    let n = 4096;
    let states: Vec<bool> = (0..n).map(|i| valid(period * i as f64 / n as f64)).collect();
    let Some(start) = states.iter().position(|&s| !s) else {
        return vec![(0.0, period, true)];
    };
    let a = period * start as f64 / n as f64;
    runs(valid, a, a + period, n).into_iter().map(|(x, y)| (x, y, false)).collect()
}

/// The screw cut: the rod's section less the screwed tool's, screwed up the rod.
fn screw_cut(model: &mut Model, blank_frame: &Frame, radius: f64, height: f64, tool: &Revolved, step: &MotionStep, motion: &Motion) -> Result<Shape> {
    let blank = Coordinates::of(blank_frame);
    let size = 2.0 * radius;
    let a = dir(step.axis.direction);
    if length3(cross3(a, blank.z)) > 1e-12 || length3(cross3(sub3(pt(step.axis.origin), blank.origin), blank.z)) > 1e-9 * size {
        return Err(Error::InvalidInput("a screw cuts a rod about its own axis".into()));
    }
    let k = step.advance / step.turn;
    // The tool's extent along the axis, in the blank's heights; at the motion's ends it must
    // lie wholly beyond the blank's ends, so that every point of the blank the screw carries
    // the tool through is passed through within the motion.
    let up = blank.z;
    let centre_height = dot3(sub3(tool.centre, blank.origin), up);
    let (low, high) = (centre_height - tool.reach(scale3(up, -1.0)), centre_height + tool.reach(up));
    let rise = |t: f64| step.advance * t * dot3(a, up);
    let (first, last) = (rise(motion.from), rise(motion.to));
    let clear = if last > first { high + first <= 0.0 && low + last >= height } else { low + first >= height && high + last <= 0.0 };
    if !clear {
        return Err(Error::InvalidInput(
            "a screw cut needs the tool wholly beyond the rod's two ends where the motion starts and ends: its ends would bound the cut otherwise, which is not taken".into(),
        ));
    }
    let quotient = Quotient::Screw { blank, k };
    let field = Field { origin: pt(step.axis.origin), axis: a, k };
    let candidates = revolved_candidates(tool, field, quotient);
    // (A point of the section's plane, screwed by s about the blank's z axis, rises k·s: it is
    // in the tool as placed only for the turns carrying it into the tool's heights.)
    let (s0, s1) = ((low / k).min(high / k), (low / k).max(high / k));
    let margin = 1e-3 * (s1 - s0) + 1e-9;
    let tool_size = 2.0 * tool.reach(tool.axis).max(tool.reach(tool.e1));
    let cover = QuotientCover { tool, quotient, range: (s0 - margin, s1 + margin), periodic: false, step_length: 0.02 * tool_size, enough: 1e-11 * size };
    let loops = section::region(&candidates, &cover, Outline::Disc { centre: [0.0, 0.0], radius })?;
    let map = PlaneMap { origin: blank.origin, x: blank.x, y: blank.y, splines: true };
    let faces = section_faces(model, &candidates, &loops, &map, size)?;
    let axis = Axis::new(blank_frame.origin(), blank_frame.z());
    Ok(screw_extrude(model, faces, axis, TAU * k, height)?.shape)
}

/// The turn cut: the ring sector's meridian section less the turned tool's, revolved
/// through the sector.
#[allow(clippy::too_many_arguments)]
fn turn_cut(model: &mut Model, blank_frame: &Frame, radii: [f64; 2], heights: [f64; 2], angles: [f64; 2], tool: &Revolved, step: &MotionStep, motion: &Motion) -> Result<Shape> {
    let blank = Coordinates::of(blank_frame);
    if !(0.0 <= radii[0] && radii[0] < radii[1] && heights[0] < heights[1] && angles[0] < angles[1] && angles[1] - angles[0] <= TAU) {
        return Err(Error::InvalidInput("a ring sector needs increasing radii (from zero or more), heights and angles".into()));
    }
    let size = (radii[1] - radii[0]).hypot(heights[1] - heights[0]);
    let a = dir(step.axis.direction);
    if length3(cross3(a, blank.z)) > 1e-12 || length3(cross3(sub3(pt(step.axis.origin), blank.origin), blank.z)) > 1e-9 * size {
        return Err(Error::InvalidInput("a turn cuts a ring sector about its own axis".into()));
    }
    // The turns about the blank's z axis the motion makes.
    let sense = dot3(a, blank.z).signum();
    let (first, last) = (sense * step.turn * motion.from, sense * step.turn * motion.to);
    let (first, last) = (first.min(last), first.max(last));
    // The polar angles the tool's part beyond the sector's inner radius spans as placed
    // (sampled on its surface, its extremes lying there, with a margin for the samples'
    // spacing): every turn carrying it over the sector must lie within the motion's.
    if last - first < TAU {
        let meridian = &tool.meridian;
        let (us, fs) = (1440, 64);
        let reference = {
            let l = blank.local(tool.centre);
            l[1].atan2(l[0])
        };
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in 0..meridian.pieces.len() {
            for j in 0..=fs {
                let m = meridian.point(i, j as f64 / fs as f64);
                for k in 0..us {
                    let l = blank.local(tool.at(TAU * k as f64 / us as f64, m));
                    if l[0].hypot(l[1]) >= radii[0] {
                        let angle = (l[1].atan2(l[0]) - reference + PI).rem_euclid(TAU) - PI + reference;
                        lo = lo.min(angle);
                        hi = hi.max(angle);
                    }
                }
            }
        }
        let spacing = TAU / us as f64 + 0.01;
        let (lo, hi) = (lo - spacing, hi + spacing);
        let needed = (angles[0] - hi, angles[1] - lo);
        // (Some whole number of turns may shift what is needed into the motion's.)
        let fits = lo.is_finite() && (-2..=2).any(|w: i32| needed.0 + TAU * w as f64 >= first && needed.1 + TAU * w as f64 <= last);
        if lo.is_finite() && (hi - lo >= TAU || !fits) {
            return Err(Error::InvalidInput(
                "a turn cut needs the tool clear of the ring sector where the motion starts and ends: its ends would bound the cut otherwise, which is not taken".into(),
            ));
        }
    }
    let quotient = Quotient::Turn { blank };
    let field = Field { origin: pt(step.axis.origin), axis: a, k: 0.0 };
    let candidates = revolved_candidates(tool, field, quotient);
    let tool_size = 2.0 * tool.reach(tool.axis).max(tool.reach(tool.e1));
    let cover = QuotientCover { tool, quotient, range: (-PI, PI), periodic: true, step_length: 0.02 * tool_size, enough: 1e-11 * size };
    let loops = section::region(&candidates, &cover, Outline::Rectangle { min: [radii[0], heights[0]], max: [radii[1], heights[1]] })?;
    let (s, c) = angles[0].sin_cos();
    let map = PlaneMap { origin: blank.origin, x: add3(scale3(blank.x, c), scale3(blank.y, s)), y: blank.z, splines: false };
    let faces = section_faces(model, &candidates, &loops, &map, size)?;
    let axis = Axis::new(blank_frame.origin(), blank_frame.z());
    Ok(super::revolve(model, faces, axis, angles[1] - angles[0])?.shape)
}

// ---- The section made faces.

/// A plane's coordinates: (u, v) ↦ origin + u x + v y. With `splines`, the blank's arcs are
/// made rational B-splines (exactly), for a sweep that carries each pole along.
struct PlaneMap {
    origin: V3,
    x: V3,
    y: V3,
    splines: bool,
}

impl PlaneMap {
    fn point(&self, p: P2) -> Point {
        let q = add3(self.origin, add3(scale3(self.x, p[0]), scale3(self.y, p[1])));
        Point::new(q[0], q[1], q[2])
    }
}

/// The faces of the section's loops in the plane `map` (a face for each outer loop, with the
/// holes inside it): a compound if more than one.
fn section_faces(model: &mut Model, candidates: &[Candidate], loops: &[Loop], map: &PlaneMap, size: f64) -> Result<Shape> {
    if loops.is_empty() {
        return Err(Error::Failed("the tool passes through all of the blank: nothing is left".into()));
    }
    let polygons: Vec<Vec<P2>> = loops.iter().map(|l| l.polygon(candidates)).collect();
    let areas: Vec<f64> = polygons.iter().map(|p| section::signed_area(p)).collect();
    let normal = cross3(map.x, map.y);
    let frame = Frame::new(map.point([0.0, 0.0]), Direction::new(normal[0], normal[1], normal[2]), Direction::new(map.x[0], map.x[1], map.x[2]));
    let mut faces = Vec::new();
    for i in 0..polygons.len() {
        if areas[i] <= 0.0 {
            continue;
        }
        let mut wires = vec![loop_wire(model, candidates, &loops[i], map, size)?];
        for (j, hole) in polygons.iter().enumerate() {
            if areas[j] >= 0.0 {
                continue;
            }
            // (A hole goes in the least outer loop holding it.)
            let holder = (0..polygons.len())
                .filter(|&k| areas[k] > 0.0 && section::contains(&polygons[k], hole[0]))
                .min_by(|&x, &y| areas[x].total_cmp(&areas[y]));
            if holder == Some(i) {
                wires.push(loop_wire(model, candidates, &loops[j], map, size)?);
            }
        }
        faces.push(crate::make::face_bounded_on(model, Surface::Plane(Plane::from_frame(frame)), &wires)?);
    }
    match faces[..] {
        [face] => Ok(face),
        _ => {
            let compound = model.new_empty(NodeKind::Compound);
            for f in faces {
                model.add(compound, f)?;
            }
            Ok(compound)
        }
    }
}

/// The wire of a loop: its junctions vertices, its pieces edges (the candidates' cubics
/// through their points, the blank's own circle and lines exactly).
fn loop_wire(model: &mut Model, candidates: &[Candidate], l: &Loop, map: &PlaneMap, size: f64) -> Result<Shape> {
    let n = l.pieces.len();
    let vertices: Vec<_> = l.junctions.iter().map(|&j| model.add_vertex(map.point(j), tolerance::LINEAR)).collect();
    let normal = cross3(map.x, map.y);
    let mut coedges = Vec::with_capacity(n);
    for (i, piece) in l.pieces.iter().enumerate() {
        let (start, end) = (l.junctions[i], l.junctions[(i + 1) % n]);
        let (va, vb) = (vertices[i], vertices[(i + 1) % n]);
        let (curve, first, last) = match *piece {
            Piece::Candidate { .. } => {
                let fit = fit_piece(candidates, piece, start, end, size)?;
                let poles: Vec<Point> = fit.poles().iter().map(|p| map.point([p.x, p.y])).collect();
                let curve = BSplineCurve::try_new(poles, fit.knots().to_vec(), fit.multiplicities().to_vec(), fit.degree())?;
                (Curve::BSpline(curve), 0.0, 1.0)
            }
            Piece::Arc { centre, radius, from, to } => {
                let z = if to > from { normal } else { scale3(normal, -1.0) };
                let frame = Frame::new(map.point(centre), Direction::new(z[0], z[1], z[2]), Direction::new(map.x[0], map.x[1], map.x[2]));
                let (first, last) = if to > from { (from, to) } else { (-from, -to) };
                let circle = Curve::Circle(Circle::new(frame, radius));
                if map.splines {
                    let arc = Curve::Trimmed(crate::geometry::curve::TrimmedCurve::new(circle, first, last));
                    let spline = crate::geometry::convert::curve_to_bspline(&arc).ok_or_else(|| Error::Failed("an arc not made a B-spline".into()))?;
                    let (a, b) = (spline.first_parameter(), spline.last_parameter());
                    (Curve::BSpline(spline), a, b)
                } else {
                    (circle, first, last)
                }
            }
            Piece::Segment { .. } => {
                let (p, q) = (map.point(start), map.point(end));
                let length = p.distance(q);
                (Curve::Line(crate::geometry::curve::Line::new(p, Direction::from_vector(q - p))), 0.0, length)
            }
        };
        let edge = model.add_edge(va, vb, tolerance::LINEAR)?;
        model.set_edge_curve(edge.into(), curve, tolerance::LINEAR)?;
        model.set_edge_range(edge.into(), first, last)?;
        coedges.push(model.add_coedge(edge, Orientation::Forward)?);
    }
    Ok(Shape::from(model.add_wire(&coedges)?))
}

/// The cubic through a candidate piece's points, even in its parameter, from junction
/// `start` to junction `end`: points doubled until it lies within 1e-10 of the blank's size
/// of the candidate between them.
fn fit_piece(candidates: &[Candidate], piece: &Piece, start: P2, end: P2, size: f64) -> Result<crate::geometry::curve2::BSplineCurve2> {
    let mut n = 16;
    loop {
        let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let points: Vec<Point2> = params
            .iter()
            .enumerate()
            .map(|(i, &f)| {
                let p = if i == 0 { start } else if i == n { end } else { piece.point(candidates, f) };
                Point2::new(p[0], p[1])
            })
            .collect();
        let fit = crate::geometry::fit::interpolate_with_parameters(&points, &params)?;
        let error = (0..n)
            .flat_map(|i| [0.25, 0.5, 0.75].map(|h| (i as f64 + h) / n as f64))
            .map(|f| {
                let (p, q) = (fit.point(f), piece.point(candidates, f));
                norm(sub([p.x, p.y], q))
            })
            .fold(0.0, f64::max);
        if error <= 1e-10 * size {
            return Ok(fit);
        }
        if n >= 1 << 16 {
            return Err(Error::Failed(format!("a boundary curve of the section fitted no closer than {error:e}")));
        }
        n *= 2;
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    use super::*;
    use crate::{properties, validity};

    fn measure(model: &Model, shape: Shape) -> (f64, f64) {
        (properties::volume(model, shape).mass, properties::area(model, shape).mass)
    }

    fn close(value: f64, expected: f64, relative: f64) -> bool {
        (value - expected).abs() <= relative * expected.abs()
    }

    /// A square turned whole about a point outside it leaves of a disc the disc inside its
    /// nearest side's circle and the ring outside its corners' circle.
    #[test]
    fn square_turned_whole_leaves_a_disc_and_a_ring() {
        let mut model = Model::new();
        let blank = Blank::Disc { frame: Frame::WORLD, radius: 10.0, height: 1.0 };
        let outline = vec![Point2::new(4.0, -1.0), Point2::new(6.0, -1.0), Point2::new(6.0, 1.0), Point2::new(4.0, 1.0)];
        let tool = Tool::Prism { frame: Frame::WORLD, outline, from: -1.0, to: 2.0 };
        let motion = Motion { steps: vec![MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: TAU, advance: 0.0 }], from: 0.0, to: 1.0 };
        let cut = cut_by_motion(&mut model, &blank, &tool, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        // The disc of radius 4 and the ring from √37 to 10, a unit thick.
        let expected_volume = PI * (16.0 + 100.0 - 37.0);
        let expected_area = 2.0 * expected_volume + TAU * (4.0 + 37f64.sqrt() + 10.0);
        assert!(close(volume, expected_volume, 1e-9), "volume {volume} against {expected_volume}");
        assert!(close(area, expected_area, 1e-9), "area {area} against {expected_area}");
        assert!(validity::check(&model, cut).is_valid());
    }

    /// A torus about the rod's own axis, screwed along it, cuts the ring its tube sweeps
    /// round: the rod falls into a core and a sleeve.
    #[test]
    fn coaxial_torus_screwed_cuts_a_ring() {
        let mut model = Model::new();
        let blank = Blank::Disc { frame: Frame::WORLD, radius: 10.0, height: 5.0 };
        let meridian = vec![MeridianPiece::Arc { centre: Point2::new(5.0, 0.0), radius: 1.0, start: 0.0, end: TAU }];
        let tool = Tool::Revolved { frame: Frame::WORLD, meridian };
        let motion = Motion { steps: vec![MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: TAU, advance: 4.0 }], from: -3.0, to: 3.0 };
        let cut = cut_by_motion(&mut model, &blank, &tool, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        let expected_volume = PI * (16.0 + 100.0 - 36.0) * 5.0;
        let expected_area = 2.0 * PI * (16.0 + 100.0 - 36.0) + TAU * (4.0 + 6.0 + 10.0) * 5.0;
        assert!(close(volume, expected_volume, 1e-8), "volume {volume} against {expected_volume}");
        assert!(close(area, expected_area, 1e-8), "area {area} against {expected_area}");
    }

    /// A ball turned past a ring sector grooves it with the torus its turn sweeps: the
    /// sector's meridian rectangle less the ball's disc, turned a quarter (Pappus).
    #[test]
    fn ball_turned_across_a_ring_sector() {
        let mut model = Model::new();
        let blank = Blank::RingSector { frame: Frame::WORLD, radii: [1.0, 5.0], heights: [-2.0, 2.0], angles: [0.0, FRAC_PI_2] };
        let meridian = vec![
            MeridianPiece::Arc { centre: Point2::new(0.0, 0.0), radius: 1.0, start: -FRAC_PI_2, end: FRAC_PI_2 },
            MeridianPiece::Line { from: Point2::new(0.0, 1.0), to: Point2::new(0.0, -1.0) },
        ];
        let tool = Tool::Revolved { frame: Frame::from_z(Point::new(3.0, 0.0, 0.0), Direction::Z), meridian };
        let motion = Motion { steps: vec![MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 1.0, advance: 0.0 }], from: -1.5, to: 3.2 };
        let cut = cut_by_motion(&mut model, &blank, &tool, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        let expected_volume = FRAC_PI_2 * 3.0 * (16.0 - PI);
        // The two ends (the rectangle less the disc), the top and bottom, the inner and
        // outer walls, and the groove: a quarter of the torus.
        let expected_area = 2.0 * (16.0 - PI) + FRAC_PI_2 * 24.0 + FRAC_PI_2 * 4.0 * (1.0 + 5.0) + 3.0 * PI * PI;
        assert!(close(volume, expected_volume, 1e-8), "volume {volume} against {expected_volume}");
        assert!(close(area, expected_area, 1e-8), "area {area} against {expected_area}");
    }

    /// A face screwed: a disc about its own axis sweeps the cylinder (its rim's helicoid lying
    /// on it), as much as a disc extruded straight (its circle fitted first, within 1e-7).
    #[test]
    fn disc_screwed_is_a_cylinder() {
        let mut model = Model::new();
        let rim = crate::make::edge_from_curve(&mut model, Curve::Circle(Circle::new(Frame::WORLD, 10.0))).unwrap();
        let wire = crate::make::wire(&mut model, &[rim]).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let solid = screw_extrude(&mut model, face, Axis::new(Point::ORIGIN, Direction::Z), 4.0, 5.0).unwrap().shape;
        let (volume, area) = measure(&model, solid);
        assert!(close(volume, 500.0 * PI, 1e-6), "volume {volume}");
        assert!(close(area, 300.0 * PI, 1e-6), "area {area}");
        assert!(validity::check(&model, solid).is_valid());
    }

    /// brep-bench's torus screwed through a rod (a grinding wheel set along the helix it
    /// grinds, its whole tube circle on the rod's side characteristic): against the values
    /// other kernels agree on to 1e-9.
    #[test]
    fn torus_screwed_through_a_rod() {
        let mut model = Model::new();
        let rod = Blank::Disc { frame: Frame::WORLD, radius: 5.0, height: 20.0 };
        let h = (TAU * 5.5 / 30.0).atan();
        let frame = Frame::new(Point::new(11.5, 0.0, 0.0), Direction::new(0.0, h.cos(), -h.sin()), Direction::X);
        let wheel = Tool::Revolved { frame, meridian: vec![MeridianPiece::Arc { centre: Point2::new(6.0, 0.0), radius: 2.0, start: 0.0, end: TAU }] };
        let screw = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 1.0, advance: 30.0 / TAU };
        let motion = Motion { steps: vec![screw], from: (-108f64).to_radians(), to: 348f64.to_radians() };
        let cut = cut_by_motion(&mut model, &rod, &wheel, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        assert!(close(volume, 1460.630002527, 1e-7), "volume {volume}");
        assert!(close(area, 791.5632225563, 1e-7), "area {area}");
        assert!(validity::check(&model, cut).is_valid());
    }

    /// brep-bench's V wheel screwed through a rod: its flanks' characteristic curves and its
    /// sharp rim.
    #[test]
    fn v_wheel_screwed_through_a_rod() {
        let mut model = Model::new();
        let rod = Blank::Disc { frame: Frame::WORLD, radius: 5.0, height: 20.0 };
        let h = (TAU * 4.5 / 30.0).atan();
        let frame = Frame::new(Point::new(14.0, 0.0, 0.0), Direction::new(0.0, h.cos(), -h.sin()), Direction::X);
        let half = 2.5 * (PI / 6.0).tan();
        let (tip, near, far) = (Point2::new(10.0, 0.0), Point2::new(7.5, half), Point2::new(7.5, -half));
        let meridian = vec![MeridianPiece::Line { from: far, to: tip }, MeridianPiece::Line { from: tip, to: near }, MeridianPiece::Line { from: near, to: far }];
        let screw = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 1.0, advance: 30.0 / TAU };
        let motion = Motion { steps: vec![screw], from: (-144f64).to_radians(), to: 384f64.to_radians() };
        let cut = cut_by_motion(&mut model, &rod, &Tool::Revolved { frame, meridian }, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        assert!(close(volume, 1554.38096569, 1e-7), "volume {volume}");
        assert!(close(area, 812.521505271, 1e-7), "area {area}");
    }

    /// brep-bench's leaning torus turned across a ring sector.
    #[test]
    fn torus_turned_across_a_block() {
        let mut model = Model::new();
        let block = Blank::RingSector { frame: Frame::WORLD, radii: [8.0, 12.0], heights: [-1.0, 1.0], angles: [0.0, PI / 3.0] };
        let lean = 25f64.to_radians();
        let frame = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::new(0.0, lean.sin(), lean.cos()), Direction::X);
        let cutter = Tool::Revolved { frame, meridian: vec![MeridianPiece::Arc { centre: Point2::new(2.5, 0.0), radius: 2.0, start: 0.0, end: TAU }] };
        let motion = Motion { steps: vec![MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 1.0, advance: 0.0 }], from: (-140f64).to_radians(), to: 140f64.to_radians() };
        let cut = cut_by_motion(&mut model, &block, &cutter, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        assert!(close(volume, 57.5718996751, 1e-7), "volume {volume}");
        assert!(close(area, 115.83662531, 1e-7), "area {area}");
        assert!(validity::check(&model, cut).is_valid());
    }

    /// brep-bench's rack tooth rolled round a gear blank: involute flanks, trochoid root.
    #[test]
    fn rack_cuts_a_tooth_space() {
        let mut model = Model::new();
        let plane = Frame::new(Point::ORIGIN, Direction::Y, Direction::X);
        let blank = Blank::Disc { frame: plane, radius: 22.0, height: 6.0 };
        let corners = [[18.0, -0.842856], [24.0, -3.026666], [24.0, 3.026666], [18.0, 0.842856]];
        let tooth = Tool::Prism { frame: plane, outline: corners.iter().map(|&[x, y]| Point2::new(x, y)).collect(), from: -2.0, to: 8.0 };
        let slide = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 0.0, advance: 20.0 };
        let turn = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Y), turn: 1.0, advance: 0.0 };
        let motion = Motion { steps: vec![slide, turn], from: -PI / 3.0, to: PI / 3.0 };
        let cut = cut_by_motion(&mut model, &blank, &tooth, &motion).unwrap();
        let (volume, area) = measure(&model, cut);
        assert!(close(volume, 9042.502788955, 1e-7), "volume {volume}");
        assert!(close(area, 3872.68353226, 1e-7), "area {area}");
        assert!(validity::check(&model, cut).is_valid());
    }

    /// A motion that is no plane problem is refused, as is a screw whose ends leave the tool
    /// in the rod.
    #[test]
    fn unsupported_motions_are_refused() {
        let mut model = Model::new();
        let rod = Blank::Disc { frame: Frame::WORLD, radius: 5.0, height: 20.0 };
        let wheel = Tool::Revolved { frame: Frame::from_z(Point::new(11.5, 0.0, 0.0), Direction::Y), meridian: vec![MeridianPiece::Arc { centre: Point2::new(6.0, 0.0), radius: 2.0, start: 0.0, end: TAU }] };
        let screw = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::Z), turn: 1.0, advance: 30.0 / TAU };
        let short = Motion { steps: vec![screw], from: 0.0, to: 1.0 };
        assert!(cut_by_motion(&mut model, &rod, &wheel, &short).is_err());
        let tumble = MotionStep { axis: Axis::new(Point::ORIGIN, Direction::X), turn: 1.0, advance: 0.0 };
        let two = Motion { steps: vec![screw, tumble], from: -10.0, to: 10.0 };
        assert!(cut_by_motion(&mut model, &rod, &wheel, &two).is_err());
    }

    /// The face of the polygon through `points` (x, z) in the xz plane.
    fn meridian_polygon(model: &mut Model, points: &[(f64, f64)]) -> Shape {
        let p: Vec<Point> = points.iter().map(|&(x, z)| Point::new(x, 0.0, z)).collect();
        let edges: Vec<Shape> = (0..p.len()).map(|i| crate::make::edge_between(model, p[i], p[(i + 1) % p.len()]).unwrap()).collect();
        let wire = crate::make::wire(model, &edges).unwrap();
        crate::make::planar_face(model, wire).unwrap()
    }

    /// Two turns of a screw of the rectangle 1 wide (from 3 out) and 2 high: with a lead of 2.5
    /// the turns stay apart (volume M·angle, M = 2·(4² − 3²)/2 = 7, two turns 4π); with a lead
    /// of 1 a turn lies half in the last, refused. A slanted band 2 high but only 0.5 thick
    /// along the axis at each distance from it stays clear of itself a lead of 1 on, as the
    /// screw does (M = 0.5·3.5, two turns).
    #[test]
    fn a_screw_whose_turns_meet_is_refused() {
        let axis = Axis::new(Point::ORIGIN, Direction::Z);
        let rectangle = [(3.0, 0.0), (4.0, 0.0), (4.0, 2.0), (3.0, 2.0)];
        let band = [(3.0, 0.0), (4.0, 1.5), (4.0, 2.0), (3.0, 0.5)];
        for (profile, lead, fits, m) in [(&rectangle, 2.5, true, 7.0), (&rectangle, 1.0, false, 7.0), (&band, 1.0, true, 1.75)] {
            let mut model = Model::new();
            let face = meridian_polygon(&mut model, profile);
            let result = screw_extrude(&mut model, face, axis, lead, 2.0 * lead);
            assert_eq!(result.is_ok(), fits, "lead {lead}: {result:?}");
            if let Ok(s) = result {
                let (volume, _) = measure(&model, s.shape);
                assert!(close(volume, m * 2.0 * TAU, 1e-6), "volume {volume}");
            }
        }
    }
}
