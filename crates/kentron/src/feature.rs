//! Form features: a profile swept into a boss (material added) or a pocket (removed) on a
//! base shape, extending a length, up to a face, between two faces, to the end of the shape
//! or right through it.
//!
//! The tool is the profile (a face, or a planar wire) swept over the span its extent needs
//! (along a direction, possibly with a draft angle, or about an axis), cut off beyond the
//! limiting faces (each the boundary of a half-space, planes extended so that the cut
//! reaches across the tool), and, for a boss that ends with the shape, kept within the
//! shape's extent; the base is then fused with it or cut by it.

use crate::error::{Error, Result};
use crate::geometry::curve::{BezierCurve, Curve};
use crate::geometry::extrema::{curve_surface, CurveSurfaceExtrema};
use crate::geometry::surface::{Plane, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// How the profile is swept into the feature's tool.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureSweep {
    /// Along a direction.
    #[doc(alias = "BRepFeat_MakePrism")]
    Prism(Direction),
    /// Along the profile plane's normal, the sides drafted by an angle (radians; positive
    /// narrows the feature as it goes).
    #[doc(alias = "BRepFeat_MakeDPrism")]
    DraftedPrism(f64),
    /// About an axis.
    #[doc(alias = "BRepFeat_MakeRevol")]
    Revolution(Axis),
}

/// Whether the feature adds material, removes it, or stands alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureMode {
    /// A pocket: the tool cut out of the base.
    Remove,
    /// A boss: the tool fused to the base.
    Add,
    /// The tool alone (as a boss would add it).
    ToolOnly,
}

/// How far the feature extends.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureExtent {
    /// A length along the sweep (an angle in radians about the axis of a revolution);
    /// negative goes the other way.
    Length(f64),
    /// From the profile up to a face (on whichever side of the profile it is).
    Until(Shape),
    /// Between two faces.
    FromUntil(Shape, Shape),
    /// From the profile to the end of the base shape.
    UntilEnd,
    /// Right through the base shape (both ways).
    ThruAll,
    /// A length (an angle about the axis of a revolution) towards a face, stopping at it if
    /// it comes first.
    UntilLength(Shape, f64),
}

/// A form feature on `base`: `profile` (a face, or a planar wire) swept as `sweep` says
/// over `extent`, then fused to the base or cut out of it (`mode`).
#[doc(alias = "BRepFeat_Form")]
#[doc(alias = "featprism")]
#[doc(alias = "featrevol")]
#[doc(alias = "featdprism")]
pub fn feature(model: &mut Model, base: Shape, profile: Shape, sweep: FeatureSweep, mode: FeatureMode, extent: FeatureExtent) -> Result<Shape> {
    let profile = match profile.kind() {
        NodeKind::Face => profile,
        NodeKind::Wire => crate::make::planar_face(model, profile)?,
        _ => {
            // A shell or compound of faces: its first face.
            *model.unique_subshapes(profile, NodeKind::Face).first().ok_or_else(|| Error::InvalidInput("a feature profile without faces".into()))?
        }
    };
    let center = crate::properties::area(model, profile).center_of_mass;
    // A length the whole construction fits in.
    let mut bbox = crate::properties::bounding_box(model, base);
    bbox.add_box(&crate::properties::bounding_box(model, profile));
    // (Limiting faces may be unbounded, made from bare surfaces: they are reached below.)
    let height = (bbox.max - bbox.min).length().max(1e-6);
    let tool = match sweep {
        FeatureSweep::Prism(direction) => prism_tool(model, base, profile, center, Vector::from(direction), None, height, mode, extent)?,
        FeatureSweep::DraftedPrism(angle) => {
            let normal = profile_normal(model, profile)?;
            prism_tool(model, base, profile, center, normal, Some(angle), height, mode, extent)?
        }
        FeatureSweep::Revolution(axis) => revolution_tool(model, base, profile, center, axis, height, extent)?,
    };
    // As OCCT's BRepFeat_Form: of the tool split by the base (the pieces outside it for a
    // boss, inside for a pocket), only those reaching into the stretch the extent selects
    // along the sweep.
    let along = match sweep {
        FeatureSweep::Prism(direction) => Some(Vector::from(direction)),
        FeatureSweep::DraftedPrism(_) => Some(profile_normal(model, profile)?),
        FeatureSweep::Revolution(_) => None,
    };
    let tool = match along.and_then(|d| selection(model, base, center, d * (1.0 / d.length()), mode, extent).map(|r| (d, r))) {
        Some((d, range)) => select_parts(model, base, tool, center, d * (1.0 / d.length()), mode, range)?,
        None => tool,
    };
    match mode {
        FeatureMode::Remove => crate::boolean::difference(model, base, tool),
        FeatureMode::Add => crate::boolean::union(model, base, tool),
        FeatureMode::ToolOnly => Ok(tool),
    }
}

/// The stretch along the sweep (from `center` along unit `d`) a feature's tool pieces must
/// reach into to be kept, as OCCT selects them: up to (or from) a limiting face of the base,
/// between limiting faces, or within the base for a boss to its end or through it. `None`
/// keeps every piece.
fn selection(model: &Model, base: Shape, center: Point, d: Vector, mode: FeatureMode, extent: FeatureExtent) -> Option<(f64, f64)> {
    let base_faces: Vec<crate::topology::NodeId> = model.unique_subshapes(base, NodeKind::Face).iter().map(|f| f.node).collect();
    let in_base = |face: Shape| model.unique_subshapes(face, NodeKind::Face).iter().any(|f| base_faces.contains(&f.node));
    let reach = |face: Shape| crossing_extent(model, face, center, d);
    match extent {
        FeatureExtent::Length(_) => None,
        FeatureExtent::Until(until) | FeatureExtent::UntilLength(until, _) if in_base(until) => {
            let (lo, hi) = reach(until);
            // The side the centre line meets it on (as OCCT's sense of the prism).
            let size = (hi - lo).abs() + extent_along(model, base, center, d).1.abs() + 1.0;
            let ahead = line_hit(model, until, center, d, 4.0 * size).map_or(0.5 * (lo + hi) >= 0.0, |t| t >= 0.0);
            Some(if ahead { (f64::NEG_INFINITY, hi) } else { (lo, f64::INFINITY) })
        }
        FeatureExtent::FromUntil(from, until) if in_base(from) || in_base(until) => {
            let (a, b) = (reach(from), reach(until));
            Some((a.0.min(b.0), a.1.max(b.1)))
        }
        FeatureExtent::UntilEnd | FeatureExtent::ThruAll if mode != FeatureMode::Remove => Some(extent_along(model, base, center, d)),
        _ => None,
    }
}

/// Where the line `center + t d` crosses `shape`'s faces (inside them), as OCCT's
/// BRepFeat::ParametricMinMax measures along the feature's centre line: the least and
/// greatest `t`; the shape's extent along `d` if the line misses it.
fn crossing_extent(model: &Model, shape: Shape, center: Point, d: Vector) -> (f64, f64) {
    let (lo, hi) = extent_along(model, shape, center, d);
    let reach = 2.0 * (lo.abs().max(hi.abs()) + 1.0);
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    for face in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(face).map(|s| s.transformed(&face.placement)) else { continue };
        let mut hits: Vec<f64> = Vec::new();
        if let Surface::Plane(plane) = &surface {
            let n = Vector::from(plane.normal());
            let along = d.dot(n);
            if along.abs() > 1e-12 {
                hits.push((plane.origin() - center).dot(n) / along);
            }
        } else {
            let segment = Curve::Bezier(BezierCurve::new(vec![center + d * (-reach), center + d * reach]));
            if let CurveSurfaceExtrema::Isolated(found) = curve_surface(&segment, &surface) {
                hits.extend(found.iter().filter(|e| e.distance <= 1e-6 * reach.max(1.0)).map(|e| (e.t * 2.0 - 1.0) * reach));
            }
        }
        for t in hits {
            if inside_face(model, face, &surface, center + d * t) {
                range = (range.0.min(t), range.1.max(t));
            }
        }
    }
    if range.0 <= range.1 { range } else { (lo, hi) }
}

/// The pieces of `tool` split by `base` (outside it for a boss, inside for a pocket) that
/// reach into `range` along unit `d` from `center`, as one shape.
fn select_parts(model: &mut Model, base: Shape, tool: Shape, center: Point, d: Vector, mode: FeatureMode, range: (f64, f64)) -> Result<Shape> {
    let pieces = match mode {
        FeatureMode::Remove => crate::boolean::intersection(model, base, tool)?,
        _ => crate::boolean::difference(model, tool, base)?,
    };
    let delta = crate::tolerance::LINEAR * 10.0;
    let kept: Vec<Shape> = model
        .unique_subshapes(pieces, NodeKind::Solid)
        .into_iter()
        .filter(|&solid| {
            let (lo, hi) = crossing_extent(model, solid, center, d);
            lo <= range.1 - delta && hi >= range.0 + delta
        })
        .collect();
    if kept.is_empty() {
        // No pieces at all (the split found none): the tool as it is.
        if model.unique_subshapes(pieces, NodeKind::Solid).is_empty() {
            return Ok(tool);
        }
        return Err(Error::Failed("a feature with no part of its tool where it is to be".into()));
    }
    // For a pocket the kept pieces are cut out, for a boss added.
    compound_of(model, &kept)
}

fn compound_of(model: &mut Model, shapes: &[Shape]) -> Result<Shape> {
    if let [one] = shapes {
        return Ok(*one);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for &s in shapes {
        model.add(compound, s)?;
    }
    Ok(compound)
}

/// How a rib's (or slot's) profile is thickened.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RibSweep {
    /// To `side` on one side of the profile's plane and `other_side` on the other (each a
    /// vector, usually across the plane).
    #[doc(alias = "BRepFeat_MakeLinearForm")]
    Linear { side: Vector, other_side: Vector },
    /// Turned about `axis` (in the profile's plane), as far as `side` from the plane on one
    /// side and `other_side` on the other.
    #[doc(alias = "BRepFeat_MakeRevolutionForm")]
    Revolution { axis: Axis, side: f64, other_side: f64 },
}

/// A rib (material added) or slot (removed) on `base`: the open `wire` in `plane`, closed
/// along the base's section by the plane into a region (on the right of the wire's first
/// edge, as OCCT's check point), thickened as `sweep` says.
#[doc(alias = "BRepFeat_RibSlot")]
#[doc(alias = "featlf")]
#[doc(alias = "featrf")]
pub fn rib(model: &mut Model, base: Shape, wire: Shape, plane: &Plane, sweep: RibSweep, mode: FeatureMode) -> Result<Shape> {
    let profile = rib_profile(model, base, wire, plane, None)?;
    let region = profile.region;
    let tool = match sweep {
        RibSweep::Linear { side, other_side } => {
            let prism = |model: &mut Model, region: Shape| -> Result<Shape> {
                let start = if other_side.length() > 0.0 {
                    crate::modify::transform(model, region, &Transform::translation(other_side), Default::default())?
                } else {
                    region
                };
                Ok(crate::sweep::extrude(model, start, side - other_side)?.shape)
            };
            // As OCCT: a rib slides when, from where its wire meets the base, it runs along
            // the base's boundary both ways; one that does not reaches down through the base
            // and is cut back to it, to sit on the faces it meets.
            let slides = profile.ends.is_none_or(|(a, b)| {
                [a, b].into_iter().all(|p| {
                    [side, other_side].into_iter().filter(|v| v.length() > 0.0).all(|v| {
                        [0.5, 1.0].into_iter().all(|f| matches!(crate::boolean::classify_point(model, base, p + v * f), Ok(crate::boolean::PointLocation::Boundary)))
                    })
                })
            });
            if slides || mode == FeatureMode::Remove {
                prism(model, region)?
            } else {
                let scale = {
                    let b = crate::properties::bounding_box(model, base);
                    (b.max - b.min).length().max(1e-9)
                };
                let deep = rib_profile(model, base, wire, plane, Some(scale))?;
                let tool = prism(model, deep.region)?;
                let outside = crate::boolean::difference(model, tool, base)?;
                // The piece holding the rib (its check point, half way across).
                let inner = profile.check + (side + other_side) * 0.5;
                let pieces = model.unique_subshapes(outside, NodeKind::Solid);
                let kept: Vec<Shape> = pieces
                    .iter()
                    .copied()
                    .filter(|&piece| matches!(crate::boolean::classify_point(model, piece, inner), Ok(crate::boolean::PointLocation::Inside)))
                    .collect();
                if kept.is_empty() { prism(model, region)? } else { compound_of(model, &kept)? }
            }
        }
        RibSweep::Revolution { axis, side, other_side } => {
            // As OCCT: turned each way by asin(thickness / radius) + π/10 (at most π/2), then
            // cut to the slab between the two sides.
            let radius = extent_from_axis(model, region, axis);
            let turn = |h: f64| ((h / radius).clamp(-1.0, 1.0).asin() + std::f64::consts::PI / 10.0).min(std::f64::consts::FRAC_PI_2);
            let (a1, a2) = (turn(side), turn(other_side));
            // Always reaching into the base (on its side of the axis) and cut back to it, as a
            // linear rib that does not slide: closed along the base, the tool's faces would
            // lie on the base's own as the section turned round, not matched by the boolean.
            let deep = match profile.ends {
                Some((a, b)) if mode != FeatureMode::Remove => {
                    let from_axis = |p: Point| {
                        let d = Vector::from(axis.direction);
                        let r = p - axis.origin;
                        (r - d * r.dot(d)).length()
                    };
                    let depth = 0.5 * from_axis(a).min(from_axis(b));
                    if depth > 0.0 { Some(rib_profile(model, base, wire, plane, Some(depth))?.region) } else { None }
                }
                _ => None,
            };
            let region = deep.unwrap_or(region);
            let start = crate::modify::transform(model, region, &Transform::rotation(axis, -a2), Default::default())?;
            let mut tool = crate::sweep::revolve(model, start, axis, a1 + a2)?.shape;
            let normal = profile_normal(model, region)?;
            let center = crate::properties::area(model, region).center_of_mass;
            let size = (crate::properties::bounding_box(model, tool).max - crate::properties::bounding_box(model, tool).min).length();
            for (h, sign) in [(side, 1.0), (other_side, -1.0)] {
                let at = center + normal * (sign * h);
                let cut = Plane::new(at, Direction::try_from_vector(normal).ok_or_else(|| Error::InvalidInput("a degenerate rib plane".into()))?);
                let (u, v) = cut.parameters(at);
                let bound = crate::make::face_on_domain(model, Surface::Plane(cut), u - 2.0 * size, u + 2.0 * size, v - 2.0 * size, v + 2.0 * size)?;
                let half = crate::primitives::half_space(model, bound, at + normal * (sign * 1e-3 * size.max(1e-9)))?;
                tool = crate::boolean::difference(model, tool, half)?;
            }
            if deep.is_some() {
                let outside = crate::boolean::difference(model, tool, base)?;
                let pieces = model.unique_subshapes(outside, NodeKind::Solid);
                let kept: Vec<Shape> = pieces
                    .iter()
                    .copied()
                    .filter(|&piece| matches!(crate::boolean::classify_point(model, piece, profile.check), Ok(crate::boolean::PointLocation::Inside)))
                    .collect();
                if !kept.is_empty() {
                    tool = compound_of(model, &kept)?;
                }
            }
            tool
        }
    };
    match mode {
        FeatureMode::Remove => crate::boolean::difference(model, base, tool),
        FeatureMode::Add => crate::boolean::union(model, base, tool),
        FeatureMode::ToolOnly => Ok(tool),
    }
}

/// A piece of a curve, from `first` to `last` (decreasing if it runs backwards).
#[derive(Clone)]
struct CurvePiece {
    curve: Curve,
    first: f64,
    last: f64,
}

impl CurvePiece {
    fn at(&self, s: f64) -> Point {
        self.curve.point(self.first + (self.last - self.first) * s)
    }

    /// The part from fraction `a` to fraction `b` of the way along.
    fn part(&self, a: f64, b: f64) -> CurvePiece {
        let t = |s: f64| self.first + (self.last - self.first) * s;
        CurvePiece { curve: self.curve.clone(), first: t(a), last: t(b) }
    }
}

/// A chain of pieces, sampled: points with the piece and fraction they are at.
struct Chain {
    pieces: Vec<CurvePiece>,
    samples: Vec<(Point, usize, f64)>,
}

impl Chain {
    fn new(pieces: Vec<CurvePiece>) -> Chain {
        const N: usize = 32;
        let mut samples = Vec::new();
        for (k, p) in pieces.iter().enumerate() {
            for i in 0..=N {
                if i == 0 && k > 0 {
                    continue;
                }
                let s = i as f64 / N as f64;
                samples.push((p.at(s), k, s));
            }
        }
        Chain { pieces, samples }
    }
}

/// The planar region a rib or slot is made of: the open `wire` in `plane` closed along the
/// boundary of `base`'s section by the plane, on the side of OCCT's check point (right of
/// the wire's first edge as it runs, seen along the plane's normal). A closed wire bounds
/// the region itself.
///
/// With `deep`, the wire's ends are joined instead through the base, `deep` into it (straight
/// across, away from the wire), as OCCT's profile of a rib that does not slide on the base.
/// Also gives where the wire meets the base (first and last) and the check point.
fn rib_profile(model: &mut Model, base: Shape, wire: Shape, plane: &Plane, deep: Option<f64>) -> Result<RibProfile> {
    let normal = Vector::from(plane.normal());
    // The wire's edges in order, as it runs.
    let mut wire_pieces = Vec::new();
    for coedge in model.children(wire).collect::<Vec<_>>() {
        let id = crate::topology::CoEdgeId::try_from(coedge.node).map_err(|_| Error::InvalidInput("a rib profile that is not a wire".into()))?;
        let reversed = model[id].orientation.compose(coedge.orientation) == Orientation::Reversed;
        let (curve, first, last) = model.edge_curve(Shape::from(model[id].edge)).ok_or_else(|| Error::InvalidInput("a rib profile edge without a curve".into()))?;
        let curve = curve.transformed(&wire.placement);
        wire_pieces.push(if reversed { CurvePiece { curve, first: last, last: first } } else { CurvePiece { curve, first, last } });
    }
    if wire_pieces.is_empty() {
        return Err(Error::InvalidInput("an empty rib profile".into()));
    }
    let mut wire_pieces = in_order(wire_pieces);
    let wire_start = wire_pieces[0].at(0.0);
    let wire_end = wire_pieces[wire_pieces.len() - 1].at(1.0);
    let scale = {
        let b = crate::properties::bounding_box(model, base);
        (b.max - b.min).length().max(1e-9)
    };
    let tolerance = 1e-7 * scale.max(1.0);
    if wire_start.distance(wire_end) <= tolerance {
        let w = crate::make::wire(model, &model.unique_subshapes(wire, NodeKind::Edge))?;
        let region = crate::make::planar_face(model, w)?;
        return Ok(RibProfile { region, ends: None, check: wire_start });
    }
    // The check point: off the middle of the first edge, to the right of it in the plane.
    let check = {
        let p = &wire_pieces[0];
        let mid = 0.5 * (p.first + p.last);
        let (q, mut t) = p.curve.d1(mid);
        if p.last < p.first {
            t = -t;
        }
        let side = t.cross(normal) * -1.0;
        q + side * (0.01 * scale / side.length().max(1e-300))
    };
    // The base's section by the plane, chained into loops.
    let (u0, v0) = plane.parameters(wire_start);
    let big = crate::make::face_on_domain(model, Surface::Plane(*plane), u0 - 4.0 * scale, u0 + 4.0 * scale, v0 - 4.0 * scale, v0 + 4.0 * scale)?;
    let section = crate::boolean::section(model, base, big)?;
    let mut free: Vec<CurvePiece> = model
        .unique_subshapes(section, NodeKind::Edge)
        .into_iter()
        .filter_map(|e| model.edge_curve(e).map(|(c, a, b)| CurvePiece { curve: c.transformed(&e.placement), first: a, last: b }))
        .collect();
    let mut loops: Vec<Vec<CurvePiece>> = Vec::new();
    while let Some(first) = free.pop() {
        let mut chain = vec![first];
        loop {
            let end = chain[chain.len() - 1].at(1.0);
            if end.distance(chain[0].at(0.0)) <= 1e3 * tolerance && chain.len() > 1 {
                break;
            }
            let Some(k) = free.iter().position(|p| p.at(0.0).distance(end) <= 1e3 * tolerance || p.at(1.0).distance(end) <= 1e3 * tolerance) else { break };
            let next = free.remove(k);
            chain.push(if next.at(0.0).distance(end) <= 1e3 * tolerance { next } else { CurvePiece { first: next.last, last: next.first, ..next } });
        }
        loops.push(chain);
    }
    // As OCCT: ends of the wire off the base are carried on along their tangents (as far
    // as the base's size), to meet it.
    let loop_samples: Vec<Vec<Point>> = loops.iter().map(|lp| Chain::new(lp.clone()).samples.iter().map(|s| s.0).collect()).collect();
    let off_base = |p: Point| {
        !loop_samples.iter().any(|ps| {
            ps.windows(2).any(|w| {
                let seg = w[1] - w[0];
                let t = if seg.dot(seg) > 0.0 { ((p - w[0]).dot(seg) / seg.dot(seg)).clamp(0.0, 1.0) } else { 0.0 };
                p.distance(w[0] + seg * t) <= 1e3 * tolerance
            })
        })
    };
    let (start_off, end_off) = (off_base(wire_start), off_base(wire_end));
    let line = |a: Point, b: Point| CurvePiece { curve: Curve::Bezier(BezierCurve::new(vec![a, b])), first: 0.0, last: 1.0 };
    let unit = |v: Vector| v * (1.0 / v.length().max(1e-300));
    if start_off {
        let p = &wire_pieces[0];
        let (_, t) = p.curve.d1(p.first);
        let t = unit(if p.last < p.first { -t } else { t });
        wire_pieces.insert(0, line(wire_start - t * scale, wire_start));
    }
    if end_off {
        let p = &wire_pieces[wire_pieces.len() - 1];
        let (_, t) = p.curve.d1(p.last);
        let t = unit(if p.last < p.first { -t } else { t });
        wire_pieces.push(line(wire_end, wire_end + t * scale));
    }
    let (start_extension, end_extension) = (start_off.then_some(0), end_off.then_some(wire_pieces.len() - 1));
    // Where the wire crosses a loop (or ends on it): positions along the wire's samples and
    // along the loop's.
    let uv = |p: Point| {
        let (u, v) = plane.parameters(p);
        crate::math::Point2::new(u, v)
    };
    let wire_chain = Chain::new(wire_pieces.clone());
    #[derive(Clone, Copy)]
    struct Crossing {
        along_wire: f64,
        wire_at: (usize, f64),
        lp: usize,
        along_loop: f64,
        loop_at: (usize, f64),
    }
    let mut crossings: Vec<Crossing> = Vec::new();
    for (li, lp) in loops.iter().enumerate() {
        let loop_chain = Chain::new(lp.clone());
        let ls = &loop_chain.samples;
        let ws = &wire_chain.samples;
        for i in 0..ws.len() - 1 {
            let (a, b) = (uv(ws[i].0), uv(ws[i + 1].0));
            for j in 0..ls.len() - 1 {
                let (c, d) = (uv(ls[j].0), uv(ls[j + 1].0));
                let (r, e) = (b - a, d - c);
                let denom = r.cross(e);
                // Ends on the loop count too.
                let near = |p: crate::math::Point2, q0: crate::math::Point2, q1: crate::math::Point2| {
                    let seg = q1 - q0;
                    let t = if seg.dot(seg) > 0.0 { ((p - q0).dot(seg) / seg.dot(seg)).clamp(0.0, 1.0) } else { 0.0 };
                    (p.distance(q0 + seg * t), t)
                };
                let mut hit: Option<(f64, f64)> = None;
                if denom.abs() > 1e-300 {
                    let t = (c - a).cross(e) / denom;
                    let w = (c - a).cross(r) / denom;
                    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&w) {
                        hit = Some((t, w));
                    }
                }
                for (end_index, p) in [(0usize, a), (ws.len() - 1, b)] {
                    let (dist, w) = near(p, c, d);
                    if dist <= 1e3 * tolerance && (i == end_index || i + 1 == end_index) {
                        hit = Some((if end_index == 0 { 0.0 } else { 1.0 }, w));
                    }
                }
                if let Some((t, w)) = hit {
                    // (A step from one piece to the next is on the next: each piece's
                    // samples after the first start past its start.)
                    let interpolate = |s: &[(Point, usize, f64)], i: usize, f: f64| {
                        let (k0, s0) = (s[i].1, s[i].2);
                        let (k1, s1) = (s[i + 1].1, s[i + 1].2);
                        if k0 == k1 { (k0, s0 + (s1 - s0) * f) } else { (k1, s1 * f) }
                    };
                    crossings.push(Crossing {
                        along_wire: i as f64 + t,
                        wire_at: interpolate(ws, i, t),
                        lp: li,
                        along_loop: j as f64 + w,
                        loop_at: interpolate(ls, j, w),
                    });
                }
            }
        }
    }
    crossings.sort_by(|x, y| x.along_wire.total_cmp(&y.along_wire));
    // The crossings nearest the wire's ends: on the wire itself, or else (an end outside the
    // base) on its carrying on.
    let on_wire = |c: &&Crossing| Some(c.wire_at.0) != start_extension && Some(c.wire_at.0) != end_extension;
    let start_outside = start_off && matches!(crate::boolean::classify_point(model, base, wire_start), Ok(crate::boolean::PointLocation::Outside));
    let end_outside = end_off && matches!(crate::boolean::classify_point(model, base, wire_end), Ok(crate::boolean::PointLocation::Outside));
    let first = crossings
        .iter()
        .filter(|c| start_outside && Some(c.wire_at.0) == start_extension)
        .next_back()
        .or_else(|| crossings.iter().find(on_wire))
        .or_else(|| crossings.iter().find(|c| Some(c.wire_at.0) != start_extension));
    let last = crossings
        .iter()
        .find(|c| end_outside && Some(c.wire_at.0) == end_extension)
        .or_else(|| crossings.iter().rev().find(on_wire))
        .or_else(|| crossings.iter().rev().find(|c| Some(c.wire_at.0) != end_extension));
    let (Some(first), Some(last)) = (first, last) else {
        return Err(Error::Failed("a rib profile that does not meet the base".into()));
    };
    if first.lp != last.lp || (last.along_wire - first.along_wire).abs() < 1e-9 {
        return Err(Error::Failed("a rib profile whose ends are not on one boundary of the base".into()));
    }
    let lp = &loops[first.lp];
    let loop_len = Chain::new(lp.clone()).samples.len() as f64 - 1.0;
    // The crossings exactly: Newton steps on the two curves in the plane's parameters.
    let refine = |w: (usize, f64), l: (usize, f64)| -> ((usize, f64), (usize, f64)) {
        let (wp, lpc) = (&wire_pieces[w.0], &lp[l.0]);
        let (mut a, mut b) = (w.1, l.1);
        let f = |a: f64, b: f64| uv(wp.at(a)) - uv(lpc.at(b));
        for _ in 0..30 {
            let r = f(a, b);
            if r.length() <= 1e-14 * scale.max(1.0) {
                break;
            }
            let h = 1e-7;
            let da = (f(a + h, b) - r) * (1.0 / h);
            let db = (f(a, b + h) - r) * (1.0 / h);
            let det = da.cross(db);
            if det.abs() <= 1e-300 {
                break;
            }
            // Solve da·x + db·y = -r.
            let x = (-r).cross(db) / det;
            let y = da.cross(-r) / det;
            a = (a + x).clamp(0.0, 1.0);
            b = (b + y).clamp(0.0, 1.0);
        }
        ((w.0, a), (l.0, b))
    };
    let (first_wire, first_loop) = refine(first.wire_at, first.loop_at);
    let (last_wire, last_loop) = refine(last.wire_at, last.loop_at);
    let first = Crossing { wire_at: first_wire, loop_at: first_loop, ..*first };
    let last = Crossing { wire_at: last_wire, loop_at: last_loop, ..*last };
    let (first, last) = (&first, &last);
    // The wire's part between the crossings.
    let wire_part = pieces_between(&wire_pieces, first.wire_at, last.wire_at);
    // The loop from the last crossing back round to the first: one way or the other.
    let forward = pieces_between_cyclic(lp, last.loop_at, first.loop_at, last.along_loop, first.along_loop, loop_len, true);
    let backward = pieces_between_cyclic(lp, last.loop_at, first.loop_at, last.along_loop, first.along_loop, loop_len, false);
    let encloses = |back: &[CurvePiece]| {
        let mut polygon: Vec<crate::math::Point2> = Vec::new();
        for p in wire_part.iter().chain(back.iter()) {
            for i in 0..16 {
                polygon.push(uv(p.at(i as f64 / 16.0)));
            }
        }
        let c = uv(check);
        let mut inside = false;
        for i in 0..polygon.len() {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            if (a.y > c.y) != (b.y > c.y) && c.x < a.x + (c.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                inside = !inside;
            }
        }
        inside
    };
    let back = if encloses(&forward) { forward } else { backward };
    let (first_point, last_point) = (wire_part[0].at(0.0), wire_part[wire_part.len() - 1].at(1.0));
    // Deep: straight down from the ends, away from the wire's middle, and across.
    let back = match deep {
        Some(depth) => {
            let chord = last_point - first_point;
            let mut down = normal.cross(chord);
            down = down * (1.0 / down.length().max(1e-300));
            // (Away from the wire's middle; a straight wire's region, and the base, are on
            // the check point's side.)
            let middle = wire_part[wire_part.len() / 2].at(0.5);
            let bulge = (middle - first_point).dot(down);
            if bulge > 1e-6 * scale || (bulge.abs() <= 1e-6 * scale && (check - first_point).dot(down) < 0.0) {
                down = -down;
            }
            let line = |a: Point, b: Point| CurvePiece { curve: Curve::Bezier(BezierCurve::new(vec![a, b])), first: 0.0, last: 1.0 };
            let (a, b) = (last_point + down * depth, first_point + down * depth);
            vec![line(last_point, a), line(a, b), line(b, first_point)]
        }
        None => back,
    };
    // The closed loop of new edges, sharing their vertices.
    let all: Vec<CurvePiece> = wire_part.into_iter().chain(back).filter(|p| p.at(0.0).distance(p.at(1.0)) > tolerance || (p.first - p.last).abs() > 1e-12).collect();
    let corners: Vec<crate::topology::VertexId> = all.iter().map(|p| model.add_vertex(p.at(0.0), tolerance)).collect();
    let mut edges = Vec::with_capacity(all.len());
    for (k, p) in all.iter().enumerate() {
        let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
        let edge = model.add_edge(a, b, tolerance)?;
        let (lo, hi) = (p.first.min(p.last), p.first.max(p.last));
        let mut curve = p.curve.clone();
        let (lo, hi) = if p.last < p.first {
            // Run backwards: the reversed curve over the mirrored range.
            curve.reverse();
            (curve.reversed_parameter(hi), curve.reversed_parameter(lo))
        } else {
            (lo, hi)
        };
        model.set_edge_curve(edge.into(), curve, tolerance)?;
        model.set_edge_range(edge.into(), lo, hi)?;
        edges.push(Shape::from(edge));
    }
    let closed = crate::make::wire(model, &edges)?;
    let region = crate::make::planar_face(model, closed)?;
    Ok(RibProfile { region, ends: Some((first_point, last_point)), check })
}

/// A rib's planar region, where its wire meets the base (an open wire's ends) and OCCT's
/// check point on the region's side.
struct RibProfile {
    region: Shape,
    ends: Option<(Point, Point)>,
    check: Point,
}

/// A wire's pieces in the order it runs (as OCCT's wire explorer: a wire's edges need not
/// be stored in order): from the start no other piece ends at, each on from where the
/// last ends; turned round where the wire's orientations do not run on.
fn in_order(mut pieces: Vec<CurvePiece>) -> Vec<CurvePiece> {
    let size = pieces.iter().map(|p| p.at(0.0).distance(p.at(1.0))).fold(0.0f64, f64::max).max(1e-9);
    let close = |a: Point, b: Point| a.distance(b) <= 1e-6 * size.max(1.0);
    let start = (0..pieces.len()).find(|&i| !pieces.iter().enumerate().any(|(j, q)| j != i && close(q.at(1.0), pieces[i].at(0.0)))).unwrap_or(0);
    let mut out = vec![pieces.remove(start)];
    while !pieces.is_empty() {
        let end = out[out.len() - 1].at(1.0);
        if let Some(k) = pieces.iter().position(|p| close(p.at(0.0), end)) {
            out.push(pieces.remove(k));
        } else if let Some(k) = pieces.iter().position(|p| close(p.at(1.0), end)) {
            let p = pieces.remove(k);
            out.push(CurvePiece { first: p.last, last: p.first, ..p });
        } else {
            // Not connected on: the rest as stored.
            out.append(&mut pieces);
        }
    }
    out
}

/// The pieces of a chain from position `from` to `to` (piece index and fraction, `from`
/// first along it).
fn pieces_between(pieces: &[CurvePiece], from: (usize, f64), to: (usize, f64)) -> Vec<CurvePiece> {
    let mut out = Vec::new();
    for k in from.0..=to.0 {
        let a = if k == from.0 { from.1 } else { 0.0 };
        let b = if k == to.0 { to.1 } else { 1.0 };
        if b > a + 1e-12 {
            out.push(pieces[k].part(a, b));
        }
    }
    out
}

/// The pieces of a closed chain from position `from` to `to`, going forwards (along it,
/// wrapping round) or backwards (against it); `from_at`, `to_at` are the positions as
/// sample counts along the chain (`total` of them).
fn pieces_between_cyclic(pieces: &[CurvePiece], from: (usize, f64), to: (usize, f64), from_at: f64, to_at: f64, total: f64, forwards: bool) -> Vec<CurvePiece> {
    if !forwards {
        // Backwards: the forward way from `to` to `from`, reversed.
        return pieces_between_cyclic(pieces, to, from, to_at, from_at, total, true)
            .into_iter()
            .rev()
            .map(|p| CurvePiece { first: p.last, last: p.first, ..p })
            .collect();
    }
    if to_at >= from_at {
        return pieces_between(pieces, from, to);
    }
    // Round the end of the chain and on from its start.
    let last = pieces.len() - 1;
    let mut out = pieces_between(pieces, from, (last, 1.0));
    out.extend(pieces_between(pieces, (0, 0.0), to));
    out
}

/// The greatest distance of `shape`'s edges (sampled) from `axis`.
fn extent_from_axis(model: &Model, shape: Shape, axis: Axis) -> f64 {
    let a = Vector::from(axis.direction);
    let mut far = 0.0f64;
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        if let Some((curve, first, last)) = model.edge_curve(edge) {
            for i in 0..=32 {
                let p = curve.point(first + (last - first) * i as f64 / 32.0) - axis.origin;
                far = far.max((p - a * p.dot(a)).length());
            }
        }
    }
    far.max(1e-300)
}

/// The faces an extent is limited by.
fn limits(extent: FeatureExtent) -> Vec<Shape> {
    match extent {
        FeatureExtent::Until(u) | FeatureExtent::UntilLength(u, _) => vec![u],
        FeatureExtent::FromUntil(f, u) => vec![f, u],
        _ => Vec::new(),
    }
}

/// The unit normal of a planar profile, as the face is used.
fn profile_normal(model: &Model, profile: Shape) -> Result<Vector> {
    let surface = model.surface(profile).ok_or_else(|| Error::InvalidInput("a profile without a surface".into()))?.transformed(&profile.placement);
    let Surface::Plane(plane) = surface else { return Err(Error::InvalidInput("a drafted prism of a profile not on a plane".into())) };
    let n = Vector::from(plane.normal());
    Ok(if profile.orientation == Orientation::Reversed { -n } else { n })
}

/// The tool of a prism (drafted if `draft`): swept from `a` to `b` along `d` (a unit
/// vector; `0` at the profile), as the extent needs, and trimmed by its limiting faces.
#[allow(clippy::too_many_arguments)]
fn prism_tool(
    model: &mut Model,
    base: Shape,
    profile: Shape,
    center: Point,
    d: Vector,
    draft: Option<f64>,
    height: f64,
    mode: FeatureMode,
    extent: FeatureExtent,
) -> Result<Shape> {
    let d = d * (1.0 / d.length());
    let max_side = {
        let mut bbox = crate::properties::bounding_box(model, base);
        bbox.add_box(&crate::properties::bounding_box(model, profile));
        let size = bbox.max - bbox.min;
        size.x.max(size.y).max(size.z)
    };
    // Where a limiting face is along the line through the profile's centre (the nearest
    // crossing), and the side of the profile it is on.
    let base_faces: Vec<crate::topology::NodeId> = model.unique_subshapes(base, NodeKind::Face).iter().map(|f| f.node).collect();
    let hit = |model: &Model, face: Shape| -> Result<f64> {
        line_hit(model, face, center, d, 2.0 * height).ok_or_else(|| Error::Failed("a feature's limiting face is not across its direction".into()))
    };
    // Far enough to pass the limiting faces.
    let height = limits(extent).iter().filter_map(|&f| line_hit(model, f, center, d, 1e3 * height)).fold(height, |h, t| h.max(1.5 * t.abs()));
    let (a, b, trims): (f64, f64, Vec<(Shape, f64)>) = match extent {
        FeatureExtent::Length(l) => (l.min(0.0), l.max(0.0), Vec::new()),
        // (A drafted prism narrows or widens as it goes: only somewhat past its limits, where
        // a long one would narrow to nothing.)
        FeatureExtent::Until(until) => {
            // (A pocket goes to the last crossing of the face on its way, a boss to the first,
            // as OCCT takes them.)
            let t = hit(model, until)?;
            let s = t.signum();
            let t = if mode == FeatureMode::Remove {
                let (lo, hi) = crossing_extent(model, until, center, d);
                if s > 0.0 { hi.max(t) } else { lo.min(t) }
            } else {
                t
            };
            // (Every crossing of the face on that side: a curved face of the base is met again.)
            let reach = |face: Shape| {
                let (lo, hi) = crossing_extent(model, face, center, d);
                if s > 0.0 { hi.max(t) } else { (-lo).max(-t) }
            };
            let far = if draft.is_some() { 1.25 * reach(until).abs() + 1e-3 * height } else { 2.0 * height };
            (if s > 0.0 { 0.0 } else { -far }, if s > 0.0 { far } else { 0.0 }, vec![(until, t)])
        }
        FeatureExtent::FromUntil(from, until) => {
            let (tf, tu) = (hit(model, from)?, hit(model, until)?);
            if draft.is_some() {
                let (ef, eu) = (crossing_extent(model, from, center, d), crossing_extent(model, until, center, d));
                let (lo, hi) = (ef.0.min(eu.0).min(tf.min(tu)), ef.1.max(eu.1).max(tf.max(tu)));
                let margin = 0.25 * (hi - lo) + 1e-3 * height;
                (lo - margin, hi + margin, vec![(from, tf), (until, tu)])
            } else {
                (-2.0 * height, 2.0 * height, vec![(from, tf), (until, tu)])
            }
        }
        // A drafted prism widens or narrows as it goes: as long as OCCT makes it (the largest
        // side of the bounding box, along its slanted sides), not longer.
        FeatureExtent::UntilLength(until, length) => {
            let t = hit(model, until)?;
            let l = length.abs() * t.signum();
            (l.min(0.0), l.max(0.0), vec![(until, t)])
        }
        FeatureExtent::UntilEnd => (0.0, match draft { Some(angle) => max_side * angle.cos(), None => 2.0 * height }, Vec::new()),
        FeatureExtent::ThruAll => match draft {
            Some(angle) => (-max_side * angle.cos(), max_side * angle.cos(), Vec::new()),
            None => (-1.5 * height, 1.5 * height, Vec::new()),
        },
    };
    if !(b > a) {
        return Err(Error::Failed("a feature of no extent".into()));
    }
    let mut tool = swept(model, profile, d, a, b, draft)?;
    // Cut off beyond each limiting face: the side away from the other limit (or from the
    // profile).
    let reference_t = match &trims[..] {
        [(_, t1), (_, t2)] => 0.5 * (t1 + t2),
        _ => 0.0,
    };
    // (A curved limiting face of the base itself, or for a pocket one on the surface of one,
    // is not cut at: the tool's pieces are selected up to its far end instead, as OCCT does;
    // it may be met more than once.)
    let base_surfaces: Vec<Surface> = model
        .unique_subshapes(base, NodeKind::Face)
        .iter()
        .filter_map(|f| model.surface(*f).map(|s| s.transformed(&f.placement)))
        .filter(|s| !matches!(s, Surface::Plane(_)))
        .collect();
    for &(face, t) in &trims {
        let curved_base_face = model.unique_subshapes(face, NodeKind::Face).iter().any(|f| {
            !matches!(model.surface(*f), Some(Surface::Plane(_)))
                && (base_faces.contains(&f.node)
                    || mode == FeatureMode::Remove && model.surface(*f).is_some_and(|s| base_surfaces.contains(&s.transformed(&f.placement))))
        });
        if curved_base_face {
            continue;
        }
        let away = if t >= reference_t { 1.0 } else { -1.0 };
        let half = half_space_beyond(model, face, center + d * t, d * away, height)?;
        tool = crate::boolean::difference(model, tool, half)?;
    }
    // A curved limit can leave pieces of the tool beyond it as well: only those between the
    // limits (or the profile and the limit) are kept.
    if !trims.is_empty() {
        let (lo, hi) = match &trims[..] {
            [(_, t)] => (t.min(0.0), t.max(0.0)),
            [(_, t1), (_, t2)] => (t1.min(*t2), t1.max(*t2)),
            _ => (f64::NEG_INFINITY, f64::INFINITY),
        };
        let solids = model.unique_subshapes(tool, NodeKind::Solid);
        if solids.len() > 1 {
            let delta = 10.0 * crate::tolerance::LINEAR;
            let kept: Vec<Shape> = solids
                .into_iter()
                .filter(|&solid| {
                    let (a, b) = crossing_extent(model, solid, center, d);
                    a <= hi - delta && b >= lo + delta
                })
                .collect();
            if !kept.is_empty() {
                tool = compound_of(model, &kept)?;
            }
        }
    }
    Ok(tool)
}

/// The profile swept from `a` to `b` along unit `d` (drafted by `draft`, about the plane
/// it starts from).
fn swept(model: &mut Model, profile: Shape, d: Vector, a: f64, b: f64, draft: Option<f64>) -> Result<Shape> {
    let start = if a != 0.0 {
        crate::modify::transform(model, profile, &Transform::translation(d * a), Default::default())?
    } else {
        profile
    };
    let swept = crate::sweep::extrude(model, start, d * (b - a))?;
    // (For the history: the profile's edges and face, through their copies where the
    // profile was moved first, sweeping the tool's sides and closing its end.)
    let uid = model.uid();
    if start.node != profile.node {
        for kind in [NodeKind::Face, NodeKind::Edge] {
            for (from, to) in model.unique_subshapes(profile, kind).into_iter().zip(model.unique_subshapes(start, kind)) {
                crate::boolean::record(uid, from.node, to.node);
            }
        }
    }
    for (input, outputs) in swept.history() {
        for out in outputs.iter().filter(|o| matches!(o.kind(), NodeKind::Face | NodeKind::Edge)) {
            crate::boolean::record(uid, *input, out.node);
        }
    }
    for &(input, _, last) in swept.ends() {
        crate::boolean::record(uid, input, last.node);
    }
    let solid = swept.shape;
    let Some(angle) = draft else { return Ok(solid) };
    // The sides tilted about the profile's plane (the faces of the solid other than its two
    // ends, which are parallel to the profile).
    let neutral_surface = model.surface(profile).cloned().ok_or_else(|| Error::InvalidInput("a profile without a surface".into()))?;
    let Surface::Plane(neutral) = neutral_surface.transformed(&profile.placement) else {
        return Err(Error::InvalidInput("a drafted prism of a profile not on a plane".into()));
    };
    let direction = Direction::try_from_vector(d).ok_or_else(|| Error::InvalidInput("a null direction".into()))?;
    let sides: Vec<crate::draft::DraftedFace> = model
        .unique_subshapes(solid, NodeKind::Face)
        .into_iter()
        .filter(|&f| match model.surface(f) {
            Some(Surface::Plane(p)) => Vector::from(p.normal()).cross(d).length() > 1e-9,
            _ => true,
        })
        .map(|face| crate::draft::DraftedFace { face, angle, neutral })
        .collect();
    crate::draft::draft(model, solid, direction, &sides)
}

/// How far `shape` reaches along unit `d` from `origin`: the least and greatest
/// coordinates of its vertices and of points along its edges.
fn extent_along(model: &Model, shape: Shape, origin: Point, d: Vector) -> (f64, f64) {
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    let mut take = |p: Point| {
        let t = (p - origin).dot(d);
        range = (range.0.min(t), range.1.max(t));
    };
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        if let Some((curve, first, last)) = model.edge_curve(edge) {
            let curve = curve.transformed(&edge.placement);
            for i in 0..=32 {
                take(curve.point(first + (last - first) * i as f64 / 32.0));
            }
        }
    }
    range
}

/// The parameter along the line `center + t d` (|t| up to `reach`) where it crosses
/// `face`'s surface, nearest the centre.
fn line_hit(model: &Model, face: Shape, center: Point, d: Vector, reach: f64) -> Option<f64> {
    let face = model.unique_subshapes(face, NodeKind::Face).into_iter().next()?;
    let surface = model.surface(face)?.transformed(&face.placement);
    // Crossings of the surface, kept where they are inside the face (unless none is: a
    // face bounded smaller than the profile still limits it).
    let mut hits: Vec<f64> = Vec::new();
    if let Surface::Plane(plane) = &surface {
        let n = Vector::from(plane.normal());
        let along = d.dot(n);
        if along.abs() < 1e-12 {
            return None;
        }
        hits.push((plane.origin() - center).dot(n) / along);
    } else {
        let segment = Curve::Bezier(BezierCurve::new(vec![center + d * (-reach), center + d * reach]));
        if let CurveSurfaceExtrema::Isolated(found) = curve_surface(&segment, &surface) {
            hits.extend(found.iter().filter(|e| e.distance <= 1e-6 * reach.max(1.0)).map(|e| (e.t * 2.0 - 1.0) * reach));
        }
    }
    let nearest = |ts: &mut dyn Iterator<Item = f64>| ts.min_by(|x: &f64, y: &f64| x.abs().total_cmp(&y.abs()));
    nearest(&mut hits.iter().copied().filter(|&t| inside_face(model, face, &surface, center + d * t))).or_else(|| nearest(&mut hits.iter().copied()))
}

/// The half-space bounded by `face`'s surface (a plane extended well past the tool) on the
/// side `beyond` points to from `at`.
fn half_space_beyond(model: &mut Model, face: Shape, at: Point, beyond: Vector, height: f64) -> Result<Shape> {
    let face = model.unique_subshapes(face, NodeKind::Face).into_iter().next().ok_or_else(|| Error::InvalidInput("a limit without a face".into()))?;
    let surface = model.surface(face).cloned().ok_or_else(|| Error::InvalidInput("a limiting face without a surface".into()))?.transformed(&face.placement);
    let bound = match &surface {
        Surface::Plane(plane) => {
            let (u, v) = plane.parameters(at);
            let r = 4.0 * height;
            crate::make::face_on_domain(model, Surface::Plane(Plane::from_frame(*plane.frame())), u - r, u + r, v - r, v + r)?
        }
        // A face over an unbounded surface (made from a bare surface): one on the surface
        // around the crossing (whole periods round, a reach across), where the cut can be made.
        other if model.face_uv_bounds(face).0.abs().max(model.face_uv_bounds(face).1.abs()).max(model.face_uv_bounds(face).2.abs()).max(model.face_uv_bounds(face).3.abs()) > 1e6 => {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(at, other, None);
            let (su0, su1, sv0, sv1) = other.bounds();
            let r = 4.0 * height;
            let span = |x: f64, lo: f64, hi: f64, periodic: bool| if periodic { (lo, hi) } else { ((x - r).max(lo), (x + r).min(hi)) };
            let (u0, u1) = span(u, su0, su1, other.is_u_periodic());
            let (v0, v1) = span(v, sv0, sv1, other.is_v_periodic());
            crate::make::face_on_domain(model, other.clone(), u0, u1, v0, v1)?
        }
        _ => Shape { orientation: Orientation::Forward, ..face },
    };
    crate::primitives::half_space(model, bound, at + beyond * (1e-3 * height))
}

/// The tool of a revolution: the profile turned about `axis` by the extent's angle, or up
/// to (or between) limiting faces, found by turning the profile's centre.
fn revolution_tool(model: &mut Model, _base: Shape, profile: Shape, center: Point, axis: Axis, height: f64, extent: FeatureExtent) -> Result<Shape> {
    use std::f64::consts::TAU;
    let turned = |p: Point, angle: f64| p.transformed(&Transform::rotation(axis, angle));
    // The first angle (one way, then the other) where the turning centre crosses a face.
    let angle_to = |model: &Model, face: Shape, from: bool| -> Option<f64> {
        let face = model.unique_subshapes(face, NodeKind::Face).into_iter().next()?;
        let surface = model.surface(face)?.transformed(&face.placement);
        let side = |angle: f64| {
            let p = turned(center, angle);
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
            let (q, du, dv) = surface.d1(u, v);
            (p - q).dot(du.cross(dv))
        };
        // Already on the face (a "from" face the profile lies on): no turn.
        if from && side(0.0).abs() <= 1e-7 * height.max(1.0) && inside_face(model, face, &surface, center) {
            return Some(0.0);
        }
        const N: usize = 360;
        // Samples half a step off the whole steps (crossings there are common), from just past
        // the start.
        let at = |sign: f64, k: usize| sign * TAU * (k as f64 + 0.5) / N as f64;
        for sign in [1.0, -1.0] {
            let mut prev = side(at(sign, 0));
            for k in 1..N {
                let cur = side(at(sign, k));
                if (prev < 0.0) != (cur < 0.0) {
                    let (mut lo, mut hi) = (at(sign, k - 1), at(sign, k));
                    for _ in 0..60 {
                        let m = 0.5 * (lo + hi);
                        if (side(m) < 0.0) == (prev < 0.0) { lo = m } else { hi = m }
                    }
                    let found = 0.5 * (lo + hi);
                    // Where the turning centre meets the face itself (not its surface beyond).
                    if inside_face(model, face, &surface, turned(center, found)) {
                        return Some(found);
                    }
                }
                prev = cur;
            }
        }
        None
    };
    let (start, sweep_angle, trims): (f64, f64, Vec<(Shape, f64)>) = match extent {
        FeatureExtent::Length(angle) => (0.0, angle, Vec::new()),
        FeatureExtent::ThruAll | FeatureExtent::UntilEnd => (0.0, TAU, Vec::new()),
        FeatureExtent::Until(until) => {
            // (Never met again, as the face the profile lies on: a whole turn.)
            let a = angle_to(model, until, false).unwrap_or(TAU);
            if a.abs() >= TAU - 1e-9 {
                return Ok(crate::sweep::revolve(model, profile, axis, TAU)?.shape);
            }
            // Past the face, then cut back to it.
            let over = (a.abs() + 0.25 * (TAU - a.abs()).min(std::f64::consts::FRAC_PI_2)).min(TAU) * a.signum();
            (0.0, over, vec![(until, a)])
        }
        FeatureExtent::UntilLength(until, angle) => {
            let a = angle_to(model, until, false).unwrap_or(TAU);
            (0.0, angle.abs() * a.signum(), vec![(until, a)])
        }
        FeatureExtent::FromUntil(from, until) => {
            let (af, au) = (
                angle_to(model, from, true).ok_or_else(|| Error::Failed("a revolution's limiting face is not on its way".into()))?,
                angle_to(model, until, false).unwrap_or(TAU),
            );
            let (lo, hi) = (af.min(au), af.max(au));
            let margin = 0.25 * (TAU - (hi - lo)).min(std::f64::consts::FRAC_PI_2);
            (lo - margin, (hi - lo) + 2.0 * margin, vec![(from, af), (until, au)])
        }
    };
    // Limits on planes through the axis are reached exactly by the angles: no trimming.
    let through_axis = |model: &Model, face: Shape| {
        let face = model.unique_subshapes(face, NodeKind::Face).into_iter().next();
        match face.and_then(|f| model.surface(f).map(|s| s.transformed(&f.placement))) {
            Some(Surface::Plane(plane)) => {
                let n = Vector::from(plane.normal());
                Vector::from(axis.direction).dot(n).abs() <= 1e-9 && (axis.origin - plane.origin()).dot(n).abs() <= 1e-7 * height.max(1.0)
            }
            _ => false,
        }
    };
    let (start, sweep_angle, trims) = if !trims.is_empty() && trims.iter().all(|&(f, _)| through_axis(model, f)) {
        match &trims[..] {
            // Up to the face, or less far if the extent says so.
            [(_, a)] => (0.0, if sweep_angle.abs() < a.abs() { sweep_angle } else { *a }, Vec::new()),
            [(_, af), (_, au)] => (*af, au - af, Vec::new()),
            _ => (start, sweep_angle, trims),
        }
    } else {
        (start, sweep_angle, trims)
    };
    let start_profile = if start != 0.0 {
        crate::modify::transform(model, profile, &Transform::rotation(axis, start), Default::default())?
    } else {
        profile
    };
    let mut tool = crate::sweep::revolve(model, start_profile, axis, sweep_angle)?.shape;
    let middle = match &trims[..] {
        [(_, a1), (_, a2)] => 0.5 * (a1 + a2),
        _ => 0.0,
    };
    for &(face, angle) in &trims {
        let at = turned(center, angle);
        // The way the centre turns there.
        let tangent = Vector::from(axis.direction).cross(at - axis.origin);
        let away = if angle >= middle { 1.0 } else { -1.0 };
        let half = half_space_beyond(model, face, at, tangent * (away / tangent.length().max(1e-300)), height)?;
        tool = crate::boolean::difference(model, tool, half)?;
    }
    Ok(tool)
}

/// Whether `p` (on or near `face`'s surface) is inside the face: its parameters within the
/// face's loops (an even-odd count of crossings of their sampled pcurves).
fn inside_face(model: &Model, face: Shape, surface: &Surface, p: Point) -> bool {
    use crate::geometry::ParametricCurve as _;
    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, surface, None);
    let mut crossings = 0usize;
    for coedge in model.unique_subshapes(face, NodeKind::CoEdge) {
        let Ok(id) = crate::topology::CoEdgeId::try_from(coedge.node) else { continue };
        let Some((pcurve, first, last)) = model.effective_pcurve(id) else { continue };
        const N: usize = 64;
        let points: Vec<crate::math::Point2> = (0..=N).map(|i| pcurve.point(first + (last - first) * i as f64 / N as f64)).collect();
        for w in points.windows(2) {
            let (a, b) = (w[0], w[1]);
            // A ray from (u, v) along +u.
            if (a.y > v) != (b.y > v) {
                let x = a.x + (v - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > u {
                    crossings += 1;
                }
            }
        }
    }
    crossings % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cuboid;
    use crate::properties::volume;

    /// The face of `shape` on the plane `z = z` (facing either way).
    fn face_at_z(model: &Model, shape: Shape, z: f64) -> Shape {
        model
            .unique_subshapes(shape, NodeKind::Face)
            .into_iter()
            .find(|&f| matches!(model.surface(f), Some(Surface::Plane(p)) if (p.origin().z - z).abs() < 1e-9 && p.normal().z().abs() > 0.5))
            .expect("a face at z")
    }

    /// A square profile `[x0, x0 + s]²` (in x and y) at height `z`.
    fn square(model: &mut Model, x0: f64, s: f64, z: f64) -> Shape {
        Cuboid::from_corners(Point::new(x0, x0, z), Point::new(x0 + s, x0 + s, z + 1.0))
            .build(model)
            .map(|b| face_at_z(model, b, z))
            .unwrap()
    }

    #[test]
    fn a_pocket_a_hole_and_a_boss() {
        let mut model = Model::new();
        let base = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        let profile = square(&mut model, 3.0, 4.0, 10.0);
        let down = FeatureSweep::Prism(-Direction::Z);
        // A pocket 2 deep.
        let pocket = feature(&mut model, base, profile, down, FeatureMode::Remove, FeatureExtent::Length(2.0)).unwrap();
        assert!((volume(&model, pocket).mass - (1000.0 - 32.0)).abs() < 1e-6);
        // A hole right through.
        let hole = feature(&mut model, base, profile, down, FeatureMode::Remove, FeatureExtent::ThruAll).unwrap();
        assert!((volume(&model, hole).mass - (1000.0 - 160.0)).abs() < 1e-6);
        // A boss from the top up to a plane at z = 13.
        let lid = Cuboid::from_corners(Point::new(-5.0, -5.0, 13.0), Point::new(15.0, 15.0, 14.0)).build(&mut model).unwrap();
        let until = face_at_z(&model, lid, 13.0);
        let boss = feature(&mut model, base, profile, FeatureSweep::Prism(Direction::Z), FeatureMode::Add, FeatureExtent::Until(until)).unwrap();
        assert!((volume(&model, boss).mass - (1000.0 + 48.0)).abs() < 1e-6, "{}", volume(&model, boss).mass);
    }

    #[test]
    fn a_tent_rib_on_a_box() {
        // OCCT's featlf A3: a roof-shaped wire on the top of a unit box, thickened 0.3.
        let mut model = Model::new();
        let base = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let [a, b, c] = [Point::new(0.7, 0.5, 1.0), Point::new(0.5, 0.5, 1.5), Point::new(0.3, 0.5, 1.0)].map(|p| crate::make::vertex(&mut model, p));
        let e1 = crate::make::edge(&mut model, a, b).unwrap();
        let e2 = crate::make::edge(&mut model, b, c).unwrap();
        let wire = crate::make::wire(&mut model, &[e1, e2]).unwrap();
        let plane = Plane::new(Point::new(0.5, 0.5, 1.0), -Direction::Y);
        let sweep = RibSweep::Linear { side: Vector::new(0.0, 0.3, 0.0), other_side: Vector::ZERO };
        let ribbed = rib(&mut model, base, wire, &plane, sweep, FeatureMode::Add).unwrap();
        // A triangle 0.4 wide and 0.5 high, 0.3 thick: its sides and slopes added, the
        // covered strip of the top removed.
        let slope = (0.2f64 * 0.2 + 0.5 * 0.5).sqrt();
        let expected = 6.0 + 2.0 * 0.1 + 2.0 * slope * 0.3 - 0.4 * 0.3;
        assert!((crate::properties::area(&model, ribbed).mass - expected).abs() < 1e-9, "{}", crate::properties::area(&model, ribbed).mass);
        assert!((volume(&model, ribbed).mass - 1.03).abs() < 1e-9);
    }

    #[test]
    fn a_quarter_turn_boss() {
        let mut model = Model::new();
        let base = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        // A 2 × 4 profile on the top at the edge x = 0, turned about that edge (along -y) up.
        let profile = Cuboid::from_corners(Point::new(0.0, 3.0, 10.0), Point::new(2.0, 7.0, 11.0)).build(&mut model).map(|b| face_at_z(&model, b, 10.0)).unwrap();
        let axis = Axis::new(Point::new(0.0, 0.0, 10.0), -Direction::Y);
        let boss = feature(&mut model, base, profile, FeatureSweep::Revolution(axis), FeatureMode::Add, FeatureExtent::Length(std::f64::consts::FRAC_PI_2)).unwrap();
        // A quarter of a cylinder of radius 2 and length 4 on top.
        let expected = 1000.0 + 0.25 * std::f64::consts::PI * 4.0 * 4.0;
        assert!((volume(&model, boss).mass - expected).abs() < 1e-6, "{}", volume(&model, boss).mass);
    }
}
