//! Fillets made by a ball rolling along any edge between two faces.
//!
//! The ball's centre is found section by section along the edge (continued a little past
//! its ends): in the plane square to the edge, at the ball's radius from both faces. The
//! fillet surface interpolates the arcs between the contact points (exact circles across,
//! a cubic B-spline along the edge, rational). The solid is then rebuilt around it: the
//! two faces lose the strips between the edge and the contact curves, the faces at the
//! edge's ends lose (or gain) the corner cut off by the fillet surface, and the edges
//! there are shortened to where the contact curves meet them.

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{BSplineCurve2, Curve2, Line2};
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{extend, fit, project, ParametricCurve, ParametricSurface};
use crate::math::{Direction2, Point, Point2, Vector};
use crate::topology::{CoEdgeId, EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

#[path = "network.rs"]
mod network;
pub(super) use network::{blend_network, fillet_network};


/// Sections along the edge, and past each of its ends.
const SECTIONS: usize = 48;
const PAST_END: usize = 8;
/// Samples across the fillet where it meets an end face.
const ACROSS: usize = 16;
/// How near (a fraction of the fillet's width) a contact meeting an edge at the fillet's
/// end passes the edge's far vertex for it to end there.
const ONTO_VERTEX: f64 = 2e-3;
/// How near (a fraction of the fillet's width) a contact passes an edge at the fillet's end
/// to meet it.
const ON_EDGE: f64 = 1e-3;

/// A face along the edge: its surface, +1 where the surface's normal points out of the
/// solid (-1 otherwise), and the direction of the edge in its boundary (+1 if the face's
/// boundary runs along the edge's own direction).
#[derive(Clone)]
struct Side {
    face: Shape,
    surface: Surface,
    sign: f64,
    along: f64,
    /// The parameters of the last point projected (where the surface has no closed-form
    /// inverse, the next projection starts there).
    last: std::cell::Cell<Option<(f64, f64)>>,
    /// The face's range of `v`: a surface of revolution's inverse may take a point to its
    /// meridian's half across the axis (the same point, with the normal the other way).
    v_range: (f64, f64),
    /// For a surface of revolution, the face's range of `v` carried on as long as its
    /// meridian stays on its side of the axis (found once).
    own_range: std::cell::Cell<Option<(f64, f64)>>,
}

impl Side {
    /// The signed distance from the surface (positive outside the solid), the foot point
    /// and the outward unit normal there.
    /// The parameters of `p`'s foot on the surface.
    fn uv(&self, p: Point) -> Result<Point2> {
        let (u, v) = self.parameters(p);
        Ok(Point2::new(u, v))
    }

    /// The parameters of `p`'s foot on the surface: in closed form where there is one,
    /// otherwise the nearest point (from the last one found).
    fn parameters(&self, p: Point) -> (f64, f64) {
        self.foot(p).0
    }

    /// [`Self::parameters`], with the foot point and its derivatives where the search found
    /// them (to spare evaluating them again).
    fn foot(&self, p: Point) -> ((f64, f64), Option<(Point, Vector, Vector)>) {
        let Some(uv) = self.surface.parameters(p) else {
            // (Locally from the last foot, the points following one another; the whole
            // surface searched the first time.)
            return match self.last.get() {
                Some(start) => {
                    let foot = crate::geometry::extrema::newton_foot(p, &self.surface, start);
                    self.last.set(Some((foot.u, foot.v)));
                    ((foot.u, foot.v), Some((foot.point, foot.du, foot.dv)))
                }
                None => {
                    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &self.surface, None);
                    self.last.set(Some((u, v)));
                    ((u, v), None)
                }
            };
        };
        (self.closed_form(p, uv), None)
    }

    /// The parameters of `p`'s foot from the surface's closed-form inverse (`uv`), kept
    /// within the face's range of `v` on a surface of revolution.
    fn closed_form(&self, p: Point, uv: (f64, f64)) -> (f64, f64) {
        let (v0, v1) = self.v_range;
        let margin = 0.25 * (v1 - v0);
        let near = 1e-9 * (v1 - v0);
        if !(matches!(self.surface, Surface::Revolution(_)) && (v1 - v0).is_finite()) || (uv.1 >= v0 - near && uv.1 <= v1 + near) {
            return uv;
        }
        // (Found again within the face's range: off it, the point may be on the
        // meridian's half across the axis, with the normal the other way.)
        if let Surface::Revolution(r) = &self.surface {
            // (Within the face's range carried on as long as the meridian stays on its
            // side of the axis.)
            let range = match self.own_range.get() {
                Some(range) => range,
                None => {
                    use crate::geometry::ParametricCurve;
                    let (o, z) = (r.axis().origin, Vector::from(r.axis().direction));
                    let radial = |v: f64| {
                        let d = r.basis().point(v) - o;
                        d - z * d.dot(z)
                    };
                    let e = radial(0.5 * (v0 + v1));
                    let mut range = (v0, v1);
                    for k in 1..=16 {
                        let v = v0 - margin * k as f64 / 16.0;
                        if radial(v).dot(e) <= 0.0 {
                            break;
                        }
                        range.0 = v;
                    }
                    for k in 1..=16 {
                        let v = v1 + margin * k as f64 / 16.0;
                        if radial(v).dot(e) <= 0.0 {
                            break;
                        }
                        range.1 = v;
                    }
                    self.own_range.set(Some(range));
                    range
                }
            };
            return r.parameters_between(p, range.0, range.1);
        }
        uv
    }

    fn distance(&self, p: Point) -> Result<(f64, Point, Vector)> {
        let ((u, v), found) = self.foot(p);
        let (foot, du, dv) = found.unwrap_or_else(|| self.surface.d1(u, v));
        let n = du.cross(dv);
        if n.length() <= 1e-300 {
            return Err(Error::InvalidInput("a fillet section meets a singular point of a face".into()));
        }
        let n = n.normalized() * self.sign;
        Ok(((p - foot).dot(n), foot, n))
    }
}

/// A section of the fillet: the spine point on the edge (or its continuation), the ball's
/// centre and its contact points with the two faces.
struct Section {
    t: f64,
    center: Point,
    contacts: [Point; 2],
    /// The contacts' parameters on the faces.
    contact_uv: [Point2; 2],
    radius: f64,
    /// A chamfer's section: straight from contact to contact (`center` then the spine
    /// point reflected through the chord's middle, on the side a ball's centre would be).
    straight: bool,
}

/// What runs along the edge: a ball of radius `law(s)` (`s` the fraction of the edge's
/// length from its start), or a chamfer cutting back `distance` along the first face
/// (`first`, else the first side found) and along the second the other distance, or as far
/// as a cut at `angle` (radians) to the first face reaches.
pub(super) enum Profile<'a> {
    Ball(&'a dyn Fn(f64) -> f64),
    Chamfer { first: Option<crate::topology::NodeId>, distance: f64, second: ChamferSecond },
}

/// How far a chamfer cuts back along its second face.
#[derive(Clone, Copy)]
pub(super) enum ChamferSecond {
    Distance(f64),
    Angle(f64),
}

/// The ball's sections along `edge` of `shape` (and past its ends), with what the fillet
/// is built from: the edge, its faces, the contact curves and whether it is convex.
fn rolling(model: &Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<Rolling> {
    rolling_past(model, shape, edge, radius, 1)
}

/// The fillet doesn't reach the face at an end of its edge within the usual reach past it
/// (an end face oblique to the edge).
const NOT_REACHED: &str = "the fillet does not reach the face at the end of its edge";

/// [`rolling`] with the sections carried `far` times as far past the edge's ends.
fn rolling_past(model: &Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64, far: usize) -> Result<Rolling> {
    rolling_with(model, shape, edge, &Profile::Ball(radius), far)
}

/// [`rolling_past`] for a ball or a chamfer.
fn rolling_with(model: &Model, shape: Shape, edge: Shape, profile: &Profile, far: usize) -> Result<Rolling> {
    // (On the faces' surfaces carried on past their patches, else, where a section is not
    // found so, on the patches themselves.)
    rolling_on(model, shape, edge, profile, far, true).or_else(|e| rolling_on(model, shape, edge, profile, far, false).map_err(|_| e))
}

/// [`rolling_with`], on the faces' surfaces carried on past their patches (`carry`) or not.
fn rolling_on(model: &Model, shape: Shape, edge: Shape, profile: &Profile, far: usize, carry: bool) -> Result<Rolling> {
    let id = EdgeId::try_from(edge.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not an edge")))?;
    let forward = Shape { orientation: Orientation::Forward, ..edge };
    let (curve, t0, t1) = model.edge_curve(forward).ok_or_else(|| Error::InvalidInput("a degenerate edge has no fillet".into()))?;
    let (va, vb) = model.edge_vertices(forward);
    let closed = va.node == vb.node;
    // (Their surfaces carried on past their patches as far as the sections past the edge's
    // ends reach.)
    let size = match profile {
        Profile::Ball(law) => (0..=8).map(|i| law(i as f64 / 8.0)).fold(0.0, f64::max),
        Profile::Chamfer { distance, second: ChamferSecond::Distance(d), .. } => distance.max(*d),
        Profile::Chamfer { distance, .. } => 2.0 * distance,
    };
    let sides = sides_of(model, shape, id, carry.then_some((4 * far + 1) as f64 * size))?;
    // The fraction of the edge's length from its start, for the radius law.
    let lengths: Vec<f64> = {
        let n = 256;
        let mut acc = vec![0.0];
        for k in 1..=n {
            let (a, b) = (t0 + (t1 - t0) * (k - 1) as f64 / n as f64, t0 + (t1 - t0) * k as f64 / n as f64);
            let last = *acc.last().expect("lengths");
            acc.push(last + curve.point(a).distance(curve.point(0.5 * (a + b))) + curve.point(0.5 * (a + b)).distance(curve.point(b)));
        }
        acc
    };
    let fraction = |t: f64| {
        let x = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0) * (lengths.len() - 1) as f64;
        let k = (x.floor() as usize).min(lengths.len() - 2);
        let l = lengths[k] + (lengths[k + 1] - lengths[k]) * (x - k as f64);
        l / lengths[lengths.len() - 1]
    };

    // Convex where the material lies between the faces: the ball rolls inside it (at
    // distance -r from both faces), else outside (+r).
    let mid = 0.5 * (t0 + t1);
    let (o, d) = curve.d1(mid);
    let tangent = d.normalized();
    let normals = [sides[0].distance(o)?.2, sides[1].distance(o)?.2];
    let into = |i: usize| normals[i].cross(tangent * sides[i].along).normalized();
    if into(0).cross(into(1)).length() <= crate::tolerance::ANGULAR.sqrt() {
        return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
    }
    let convex = normals[0].dot(into(1)) < 0.0;
    let side_of_ball = if convex { -1.0 } else { 1.0 };
    // Ends where the faces meet tangentially (the mitre between fillets of two sides,
    // running up to the face both roll onto): a ball's fillet narrows to nothing there, its
    // last section the vertex itself, none past it.
    let vanishing: [bool; 2] = match profile {
        Profile::Ball(_) if !closed => [t0, t1].map(|t| {
            let p = curve.point(t);
            matches!((sides[0].distance(p), sides[1].distance(p)), (Ok(a), Ok(b)) if a.2.dot(b.2) >= 1.0 - crate::tolerance::ANGULAR.sqrt())
        }),
        _ => [false; 2],
    };
    // A chamfer's distances along the two sides (the second from the angle in the faces'
    // cross-section at the edge's middle, by the sines' law in its triangle).
    let chamfer: Option<[f64; 2]> = match profile {
        Profile::Ball(_) => None,
        Profile::Chamfer { first, distance, second } => {
            let swapped = first.is_some_and(|f| sides[1].face.node == f);
            let (a, b) = if swapped { (1, 0) } else { (0, 1) };
            let other = match *second {
                ChamferSecond::Distance(d) => d,
                ChamferSecond::Angle(angle) => {
                    let between = into(a).dot(into(b)).clamp(-1.0, 1.0).acos();
                    let opposite = (between + angle).sin();
                    if !(opposite > crate::tolerance::ANGULAR) {
                        return Err(Error::InvalidInput(format!("a chamfer at {angle} rad does not meet the other face")));
                    }
                    distance * angle.sin() / opposite
                }
            };
            Some(if swapped { [other, *distance] } else { [*distance, other] })
        }
    };
    // (At an angle: the face it is to, and the angle; the second distance above only a
    // first guess, each section cut at the angle to the first face there.)
    let cut_angle: Option<(usize, f64)> = match profile {
        Profile::Chamfer { first, second: ChamferSecond::Angle(angle), .. } => Some((if first.is_some_and(|f| sides[1].face.node == f) { 1 } else { 0 }, *angle)),
        _ => None,
    };
    let chamfer_size = chamfer.map(|d| d[0].max(d[1]));
    let ball_law = |s: f64| match profile {
        Profile::Ball(law) => law(s),
        Profile::Chamfer { .. } => chamfer_size.expect("a chamfer's size"),
    };
    let radius: &dyn Fn(f64) -> f64 = &ball_law;

    // The spine: points of the edge, continued past its ends (by a few radii) along the
    // faces' intersection.
    let reach = if closed { 0.0 } else { 4.0 * far as f64 * (0..=8).map(|i| radius(i as f64 / 8.0)).fold(0.0, f64::max) };
    let past = |t_end: f64| -> f64 { reach / curve.d1(t_end).1.length() };
    let (t_min, t_max) = (t0 - past(t0), t1 + past(t1));
    // (Along a closed edge whose curve turns a little where it closes on itself, both its
    // ends square to the mean of its directions there: its sections there one, the fillet
    // closing on itself.)
    let seam = closed.then(|| (curve.d1(t0).1.normalized() + curve.d1(t1).1.normalized()).normalized());
    let spine_at = |t: f64| -> Result<(Point, Vector)> {
        if (t0..=t1).contains(&t) {
            let (p, d) = curve.d1(t);
            if let Some(mean) = seam.filter(|_| t == t0 || t == t1) {
                return Ok((p, mean));
            }
            return Ok((p, d.normalized()));
        }
        // Past an end: the edge's curve continued (keeping its parametrization smooth),
        // brought onto both faces' surfaces.
        let (p, d) = curve.d1(t);
        let q = onto_both(&sides, p, d.normalized())?;
        let n = sides[0].distance(q)?.2.cross(sides[1].distance(q)?.2);
        Ok((q, if n.dot(d) < 0.0 { -n.normalized() } else { n.normalized() }))
    };
    // (Whether sections lie in the balls' own planes across, not square to the edge.)
    let characteristic = std::cell::Cell::new(false);
    // (Each section from the centre of the last found, relative to its point of the spine.)
    let last_offset: std::cell::Cell<Option<Vector>> = std::cell::Cell::new(None);
    let last_contacts: std::cell::Cell<Option<[Vector; 2]>> = std::cell::Cell::new(None);
    let section_at = |t: f64| -> Result<Section> {
        let r = radius(fraction(t));
        if !(r > crate::tolerance::LINEAR) {
            return Err(Error::InvalidInput(format!("fillet radius must be positive, got {r}")));
        }
        let (o, tangent) = spine_at(t)?;
        if let Some(distances) = chamfer {
            let near = last_contacts.get();
            let on_edge = (t0.min(t1)..=t0.max(t1)).contains(&t);
            let s = chamfer_section(&sides, o, tangent, distances, cut_angle, near, on_edge).or_else(|e| match near {
                Some(_) => chamfer_section(&sides, o, tangent, distances, cut_angle, None, on_edge),
                None => Err(e),
            })?;
            last_contacts.set(Some([s.contacts[0] - o, s.contacts[1] - o]));
            return Ok(Section { t, ..s });
        }
        let found = match section(&sides, o, tangent, r, side_of_ball, last_offset.get(), characteristic.get()) {
            Ok(found) => Ok((o, Section { t, ..found })),
            // On the edge, where the faces' surfaces meet a seam or an end of theirs: just
            // inside it.
            Err(e) if (t0.min(t1)..=t0.max(t1)).contains(&t) => {
                let inside = t + (0.5 * (t0 + t1) - t).signum() * 1e-4 * (t1 - t0).abs();
                let (o, tangent) = spine_at(inside)?;
                section(&sides, o, tangent, r, side_of_ball, last_offset.get(), characteristic.get()).map(|found| (o, Section { t: inside, ..found })).map_err(|_| e)
            }
            Err(e) => Err(e),
        };
        let (o, s) = found?;
        last_offset.set(Some(s.center - o));
        Ok(s)
    };
    // (At a vanishing end: the ball touching both faces at the vertex.)
    let vanished = |t: f64| -> Result<Section> {
        let p = curve.point(t);
        let n = sides[0].distance(p)?.2;
        let center = p + n * (side_of_ball * radius(fraction(t)));
        last_offset.set(Some(center - p));
        Ok(Section { t, center, contacts: [p, p], contact_uv: [sides[0].uv(p)?, sides[1].uv(p)?], radius: radius(fraction(t)), straight: false })
    };
    // More sections until the contact curves interpolated through them hold to fresh
    // sections midway.
    let size = curve.point(t0).distance(curve.point(0.5 * (t0 + t1))) + curve.point(t1).distance(curve.point(0.5 * (t0 + t1))) + reach;
    let target = 1e-9 * size.max(1.0);
    let sweep = || -> Result<(Vec<Section>, [BSplineCurve; 2])> {
        let mut n = SECTIONS;
        let mut previous_worst = f64::INFINITY;
        let mut previous_poor = false;
        loop {
            let past_n = if closed { 0 } else { far * (PAST_END * n / SECTIONS).max(PAST_END) };
            // Past each end, outwards as far as the ball still fits between the faces (their
            // surfaces may end or degenerate further on).
            // (Not where the spine runs off: a curve carried on past its range, onto the faces'
            // surfaces, may leap away; no farther from the end than twice the reach.)
            // (Nor past where the faces turn tangent and the ball's arc, shrunk to nothing,
            // would open the other way round: a ball there is on the faces' other side.)
            let turn = |s: &Section| (s.contacts[0] - s.center).cross(s.contacts[1] - s.center).dot(curve.d1(s.t).1);
            let past = |end: f64, far: f64| -> Vec<Section> {
                let at = curve.point(end);
                // (The way the first section's arc turns and its width, and the narrowest yet.)
                let mut way: Option<(f64, f64)> = None;
                let mut narrowest = f64::INFINITY;
                (1..=past_n)
                    .map(|j| end + (far - end) * j as f64 / past_n as f64)
                    .map_while(|t| {
                        let near = spine_at(t).ok().is_some_and(|(o, _)| o.distance(at) <= 2.0 * reach.max(crate::tolerance::LINEAR));
                        let s = if near { section_at(t).ok() } else { None }?;
                        if s.straight {
                            return Some(s);
                        }
                        let (w, width) = (turn(&s), s.contacts[0].distance(s.contacts[1]));
                        match way {
                            Some((first, wide)) if first * w < 0.0 && narrowest < 0.25 * wide => None,
                            _ => {
                                way = way.or(Some((w, width)));
                                narrowest = narrowest.min(width);
                                Some(s)
                            }
                        }
                    })
                    .collect()
            };
            // (A chamfer's sections, each run from where it starts: a section is solved from the
            // one before it, which must be its neighbour, else it may find the points on the
            // faces' surfaces carried on past the edge, on its other side. A ball's from the one
            // before, as found.)
            let restart = |from: Option<&Section>| {
                let spine = from.and_then(|s| spine_at(s.t).ok().map(|(o, _)| (s, o)));
                last_contacts.set(spine.map(|(s, o)| [s.contacts[0] - o, s.contacts[1] - o]));
            };
            restart(None);
            let before = if vanishing[0] { Vec::new() } else { past(t0, t_min) };
            restart(before.first());
            // (On from the section past the start nearest it.)
            let mut sections: Vec<Section> = before.into_iter().rev().collect();
            // (Where no ball fits near an end of an open edge (a face's surface singular there,
            // as at a cone's apex), the sections stop short of it: the fillet must end on the
            // face there before.)
            let mut short = [false; 2];
            for k in 0..=n {
                let t = t0 + (t1 - t0) * k as f64 / n as f64;
                let found = match (k, vanishing) {
                    (0, [true, _]) => vanished(t0),
                    (k, [_, true]) if k == n => vanished(t1),
                    _ => section_at(t),
                };
                match found {
                    Ok(s) => sections.push(s),
                    // (Next to a vanishing end, where the faces are all but tangent, a section
                    // may not be found: left out.)
                    Err(_) if (vanishing[0] && 8 * k < n) || (vanishing[1] && 8 * k > 7 * n) => {}
                    Err(_) if !closed && 4 * k < n && (short[0] || k == 0) => {
                        short[0] = true;
                        sections.clear();
                    }
                    Err(_) if !closed && 4 * k > 3 * n => {
                        short[1] = true;
                        break;
                    }
                    Err(e) => return Err(e),
                }
            }
            if !short[1] && !vanishing[1] {
                sections.extend(past(t1, t_max));
            }
            let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();
            let contact = |i: usize| -> Result<BSplineCurve> {
                let points: Vec<Point> = sections.iter().map(|s| s.contacts[i]).collect();
                fit::interpolate_with_parameters(&points, &ts)
            };
            let contacts = [contact(0)?, contact(1)?];
            let mut worst: f64 = 0.0;
            restart(sections.first());
            for w in ts.windows(2) {
                // (Past the edge's ends, where the faces' surfaces may end between two sections
                // found, a section not found midway is left out of the measure.)
                let mid = match section_at(0.5 * (w[0] + w[1])) {
                    Ok(mid) => mid,
                    Err(_) if w[1] <= t0.min(t1) || w[0] >= t0.max(t1) => continue,
                    Err(_) if (vanishing[0] && w[0] == t0) || (vanishing[1] && w[1] == t1) => continue,
                    Err(e) => return Err(e),
                };
                for i in 0..2 {
                    worst = worst.max(contacts[i].point(mid.t).distance(mid.contacts[i]));
                }
            }
            // (Not further where doubling the sections has twice running hardly helped, and the
            // gap is still far from closed: the contact curves are not smooth there (the faces'
            // own seams), and twice as many would not close it. Errors jumping about instead
            // (sections finding other balls) are left to the full count.)
            let ratio = previous_worst / worst;
            let poor = (0.25..4.0).contains(&ratio);
            let stalled = n >= 4 * SECTIONS && worst > 100.0 * target && poor && previous_poor;
            if worst <= target || stalled || n >= 32 * SECTIONS {
                return Ok((sections, contacts));
            }
            previous_worst = worst;
            previous_poor = poor;
            n *= 2;
        }
    };

    // (Where the balls' centres lie beyond the edge's centre of curvature (a ball near as
    // large as the curve's radius, inside its bend), the sections square to it fold back
    // (the centres run back along the edge) or are not found. Then each section in its
    // ball's own plane across (the arc's, square to the path of the centres) through the
    // edge's point instead.)
    let folds = |sections: &[Section]| {
        sections.windows(2).any(|w| {
            let inside = |t: f64| (t0.min(t1)..=t0.max(t1)).contains(&t);
            inside(w[0].t) && inside(w[1].t) && (w[1].center - w[0].center).dot(curve.d1(w[0].t).1) * (w[1].t - w[0].t) < -1e-9 * size * (w[1].t - w[0].t).abs()
        })
    };
    let (sections, contacts) = match sweep() {
        Ok(swept) if chamfer.is_some() || !folds(&swept.0) => swept,
        Err(e) if chamfer.is_some() => return Err(e),
        first => {
            characteristic.set(true);
            match sweep() {
                Ok(swept) if !folds(&swept.0) => swept,
                _ => {
                    characteristic.set(false);
                    first?
                }
            }
        }
    };
    Ok(Rolling { id, va, vb, closed, sides, sections, contacts, convex, target, t0, t1, vanishing })
}

/// What a rolling-ball fillet is made from (see [`rolling`]).
struct Rolling {
    id: EdgeId,
    va: Shape,
    vb: Shape,
    closed: bool,
    sides: [Side; 2],
    sections: Vec<Section>,
    contacts: [BSplineCurve; 2],
    convex: bool,
    target: f64,
    t0: f64,
    t1: f64,
    /// Whether the fillet narrows to nothing at each end, the faces tangent there.
    vanishing: [bool; 2],
}

/// `shape` with `edge` rounded by a ball of radius `radius(s)`, `s` the fraction of the
/// edge's length from its start.
pub(super) fn rolling_ball(model: &mut Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<Shape> {
    // (Further past the ends where the end faces are not reached.)
    match rolling_ball_past(model, shape, edge, radius, 1) {
        Err(Error::InvalidInput(m)) if m == NOT_REACHED => rolling_ball_past(model, shape, edge, radius, 3),
        other => other,
    }
}

fn rolling_ball_past(model: &mut Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64, far: usize) -> Result<Shape> {
    let Rolling { id, va, vb, closed, sides, sections, contacts, convex, target, .. } = rolling_past(model, shape, edge, radius, far)?;
    let surface = fillet_surface(&sections)?;
    let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();

    // Where the fillet meets the faces at the edge's ends.
    let mut ends = Vec::new();
    if !closed {
        for (vertex, at_end) in [(va, false), (vb, true)] {
            ends.push(end_cut(model, shape, id, vertex, &sides, &surface, &ts, at_end, target)?);
        }
    }
    let mut builder = Rebuild { model, id, sides: &sides, surface: &surface, contacts: &contacts, ends: &ends, ts: &ts };
    builder.build(shape, closed, &sections, convex)
}

/// `shape` with `edge` chamfered as `profile` (a [`Profile::Chamfer`]) has it, the faces
/// along it cut back and the chamfer ending on the faces at the edge's ends.
pub(super) fn rolling_chamfer(model: &mut Model, shape: Shape, edge: Shape, profile: &Profile) -> Result<Shape> {
    let attempt = |model: &mut Model, far: usize| -> Result<Shape> {
        let Rolling { id, va, vb, closed, sides, sections, contacts, convex, target, .. } = rolling_with(model, shape, edge, profile, far)?;
        let surface = fillet_surface(&sections)?;
        let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();
        let mut ends = Vec::new();
        if !closed {
            for (vertex, at_end) in [(va, false), (vb, true)] {
                ends.push(end_cut(model, shape, id, vertex, &sides, &surface, &ts, at_end, target)?);
            }
        }
        let mut builder = Rebuild { model, id, sides: &sides, surface: &surface, contacts: &contacts, ends: &ends, ts: &ts };
        builder.build(shape, closed, &sections, convex)
    };
    match attempt(model, 1) {
        Err(Error::InvalidInput(m)) if m == NOT_REACHED => attempt(model, 3),
        other => other,
    }
}

/// The solid a rolling-ball fillet of `edge` takes away (at a convex edge) or adds (at a
/// concave one), for a boolean operation, with whether it is convex: in each section the
/// region between the ball's arc and the edge, bounded by the rays from the ball's centre
/// through the contacts (carried well past the edge) and closed across beyond it. A convex
/// tool runs on past the edge's ends, for the faces there to trim it; a concave one stops
/// at them.
pub(super) fn rolling_tool(model: &mut Model, shape: Shape, edge: Shape, radius: &dyn Fn(f64) -> f64) -> Result<(Shape, bool)> {
    use crate::geometry::curve::Circle;
    use crate::math::{Direction, Frame};
    // (On the faces' own surfaces: the tool goes on past the edge's ends, where a boolean
    // meets it with the faces as they are.)
    let Rolling { sections, convex, closed, t0, t1, .. } = rolling_on(model, shape, edge, &Profile::Ball(radius), 1, false)?;
    if closed {
        return Err(Error::InvalidInput("rolling-ball tools along closed edges are not supported yet".into()));
    }
    let span = (t1 - t0).abs() * 1e-9;
    let kept: Vec<&Section> = sections.iter().filter(|s| convex || (s.t >= t0.min(t1) - span && s.t <= t0.max(t1) + span)).collect();
    let mut wires = Vec::with_capacity(kept.len());
    let kept_radius = kept.iter().map(|s| s.radius).fold(0.0, f64::max);
    for s in kept {
        let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        // (Not where the faces turn tangent: the arc shrinks to nothing there.)
        if half <= 1e-3 || s.contacts[0].distance(s.contacts[1]) <= 1e3 * crate::tolerance::LINEAR {
            continue;
        }
        // Far enough along the rays for the closing segment to pass beyond the edge's
        // corner (which lies within the wedge, at r / sin(half) from the centre).
        let reach = (s.radius / half.sin().max(1e-3) + 0.5 * s.radius) / half.cos().max(1e-3);
        let far = [s.center + a.normalized() * reach, s.center + b.normalized() * reach];
        let normal = a.cross(b);
        let frame = Frame::new(s.center, Direction::from_vector(normal), Direction::from_vector(a));
        // (The lines start where the arc ends, the contacts being only nearly on it.)
        let circle = Circle::new(frame, a.length());
        let ends = [circle.point(0.0), circle.point(2.0 * half)];
        let arc = crate::make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, 2.0 * half)?;
        let edges = [
            arc,
            crate::make::edge_between(model, ends[1], far[1])?,
            crate::make::edge_between(model, far[1], far[0])?,
            crate::make::edge_between(model, far[0], ends[0])?,
        ];
        wires.push(crate::make::wire(model, &edges)?);
    }
    let tool = crate::loft::loft(model, &wires, &crate::loft::LoftOptions { solid: true, ruled: false, ..Default::default() })?.shape;
    let reach = 0.05 * kept_radius;
    cover_sides(model, tool, shape, edge, reach)?;
    Ok((tool, convex))
}

/// Raises the tolerance of each edge of `tool` (and its vertices') lying on one of the faces
/// along `edge` within `reach` to cover how far it strays from that face: the lofted tool
/// only passes through its sections' points on the faces, and an operation meeting it with
/// the shape finds its edges on the faces only within that much.
fn cover_sides(model: &mut Model, tool: Shape, shape: Shape, edge: Shape, reach: f64) -> Result<()> {
    let Ok(id) = EdgeId::try_from(edge.node) else { return Ok(()) };
    let sides = sides_of(model, shape, id, None)?;
    for e in model.unique_subshapes(tool, NodeKind::Edge) {
        let Some((curve, first, last)) = model.edge_curve(e) else { continue };
        let points: Vec<Point> = (0..=32).map(|k| curve.point(first + (last - first) * k as f64 / 32.0)).collect();
        let strays = sides
            .iter()
            .filter_map(|side| points.iter().map(|&p| side.distance(p).map(|d| d.0.abs())).collect::<Result<Vec<f64>>>().ok())
            .map(|d| d.into_iter().fold(0.0, f64::max))
            .filter(|&d| d <= reach)
            .fold(f64::INFINITY, f64::min);
        if !strays.is_finite() {
            continue;
        }
        let needed = 1.5 * strays;
        let Ok(edge_id) = EdgeId::try_from(e.node) else { continue };
        if needed > model.tolerance(e) {
            model.set_edge_tolerance(edge_id, needed)?;
            let (a, b) = model.edge_vertices(e);
            for v in [a, b] {
                if needed > model.tolerance(v) {
                    model.set_vertex_tolerance(VertexId::try_from(v.node).expect("a vertex"), needed)?;
                }
            }
        }
    }
    Ok(())
}

/// The two faces along the edge.
fn sides_of(model: &Model, shape: Shape, id: EdgeId, carry: Option<f64>) -> Result<[Side; 2]> {
    let mut sides = Vec::new();
    for face in model.explore(shape, NodeKind::Face) {
        for wire in model.children(face) {
            for coedge in model.children(wire) {
                let c = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                if model[c].edge != id {
                    continue;
                }
                let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.clone();
                // (A patch carried on past its bounds, the same over them: the
                // sections past the edge's ends then follow the faces on smoothly, as
                // projected onto the patch itself they would stop at its bounds.)
                // (By `carry`, no more than half the face's size: farther, the surface carried
                // on may curl round.)
                let surface = match carry {
                    Some(reach) => {
                        let size = { let b = crate::properties::bounding_box(model, face); (b.max - b.min).length() };
                        let _ = reach;
                        carried_on(&surface, 0.5 * size).unwrap_or(surface)
                    }
                    None => surface,
                };
                let sign = |o: Orientation| if o == Orientation::Reversed { -1.0 } else { 1.0 };
                sides.push(Side {
                    face,
                    surface,
                    sign: sign(face.orientation),
                    along: sign(model[c].orientation.compose(wire.orientation)),
                    last: Default::default(),
                    v_range: { let b = model.face_uv_bounds(face); (b.2, b.3) },
                    own_range: Default::default(),
                });
            }
        }
    }
    match <[Side; 2]>::try_from(sides) {
        Ok(sides) => Ok(sides),
        Err(sides) => Err(Error::InvalidInput(format!("a fillet edge must bound two faces, it bounds {}", sides.len()))),
    }
}

/// The point on both faces' surfaces nearest to `guess` in the plane through it square to
/// `tangent` (Newton).
fn onto_both(sides: &[Side; 2], guess: Point, tangent: Vector) -> Result<Point> {
    newton(guess, |x| Ok([sides[0].distance(x)?.0, sides[1].distance(x)?.0, (x - guess).dot(tangent)]))
}

/// The ball of radius `r` touching both faces, centred in the plane through `o` square to
/// `tangent` (or, `characteristic`, in the plane through `o` of its own arc between the
/// contacts, square to the path of the centres), on the given side of the faces (-1 inside
/// the solid): from the centre at offset `near` from `o` (a neighbouring section's), else
/// from between the faces' tangent planes.
#[allow(clippy::too_many_arguments)]
fn section(sides: &[Side; 2], o: Point, tangent: Vector, r: f64, side: f64, near: Option<Vector>, characteristic: bool) -> Result<Section> {
    // (The plane's normal at a centre: the path of the centres runs square to both faces'
    // normals at the contacts.)
    let across = |q: Point| -> Result<Vector> {
        if !characteristic {
            return Ok(tangent);
        }
        let along = sides[0].distance(q)?.2.cross(sides[1].distance(q)?.2);
        Ok(if along.length() > 1e-9 { along.normalized() * along.normalized().dot(tangent).signum() } else { tangent })
    };
    if let Some(offset) = near {
        let start = o + (offset - tangent * offset.dot(tangent));
        let solved = newton_damped(start, |q| Ok([sides[0].distance(q)?.0 - side * r, sides[1].distance(q)?.0 - side * r, (q - o).dot(across(q)?)]), characteristic);
        if let Ok(center) = solved {
            let contacts = [sides[0].distance(center)?.1, sides[1].distance(center)?.1];
            let contact_uv = [sides[0].uv(center)?, sides[1].uv(center)?];
            return Ok(Section { t: 0.0, center, contacts, contact_uv, radius: r, straight: false });
        }
    }
    // (The faces' normals at `o`, or where a surface is singular there (a cone's apex), just
    // along the edge.)
    let normal = |i: usize| -> Result<Vector> {
        let reach = 1e-6 * (1.0 + (o - Point::ORIGIN).length());
        sides[i].distance(o).or_else(|e| sides[i].distance(o + tangent * reach).or_else(|_| sides[i].distance(o - tangent * reach)).map_err(|_| e)).map(|d| d.2)
    };
    let normals = [normal(0)?, normal(1)?];
    let into: Vec<Vector> = (0..2).map(|i| normals[i].cross(tangent * sides[i].along)).collect();
    let tiny = 1e-9;
    let unit = |v: Vector| if v.length() > tiny { v * (1.0 / v.length()) } else { v };
    let sum = unit(into[0]) + unit(into[1]);
    if into.iter().any(|v| v.length() <= tiny) || sum.length() <= tiny {
        return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
    }
    let into: Vec<Vector> = into.into_iter().map(|v| v.normalized()).collect();
    let half = 0.5 * into[0].dot(into[1]).clamp(-1.0, 1.0).acos();
    let start = o + sum.normalized() * (r / half.sin());
    let center = newton_damped(start, |q| {
        Ok([sides[0].distance(q)?.0 - side * r, sides[1].distance(q)?.0 - side * r, (q - o).dot(across(q)?)])
    }, characteristic)?;
    let contacts = [sides[0].distance(center)?.1, sides[1].distance(center)?.1];
    let contact_uv = [sides[0].uv(center)?, sides[1].uv(center)?];
    Ok(Section { t: 0.0, center, contacts, contact_uv, radius: r, straight: false })
}

/// The chamfer's section in the plane through `o` square to `tangent`: on each face the
/// point `distances[i]` from `o` (in that plane, into the face), from the points at offsets
/// `near` from `o` (a neighbouring section's), else from along the faces' tangent planes.
/// (Off a face's patch, in its parameters, only for a section `on_edge`: past the edge's
/// ends the sections stop where the faces do.)
fn chamfer_section(sides: &[Side; 2], o: Point, tangent: Vector, distances: [f64; 2], angle: Option<(usize, f64)>, near: Option<[Vector; 2]>, on_edge: bool) -> Result<Section> {
    let normal = |i: usize| -> Result<Vector> {
        let reach = 1e-6 * (1.0 + (o - Point::ORIGIN).length());
        sides[i].distance(o).or_else(|e| sides[i].distance(o + tangent * reach).or_else(|_| sides[i].distance(o - tangent * reach)).map_err(|_| e)).map(|d| d.2)
    };
    let mut contacts = [o; 2];
    // (At an angle to the first face, that face's contact first.)
    let order = match angle {
        Some((first, _)) => [first, 1 - first],
        None => [0, 1],
    };
    for i in order {
        let d = distances[i];
        let start = match near {
            Some(offsets) => o + (offsets[i] - tangent * offsets[i].dot(tangent)),
            None => {
                let into = normal(i)?.cross(tangent * sides[i].along);
                if into.length() <= 1e-9 {
                    return Err(Error::InvalidInput("the faces at a fillet edge are tangent".into()));
                }
                o + into.normalized() * d
            }
        };
        // At an angle to the first face: where the cut from its contact meets this face,
        // turned by the angle from the first face's way back to the edge there (from the
        // distance the faces' tangents at the edge give).
        if let Some((first, angle)) = angle.filter(|&(first, _)| first != i) {
            if let Some(x) = cut_at_angle(sides, o, tangent, contacts[first], first, angle, start) {
                contacts[i] = x;
                continue;
            }
        }
        // (Else in the face's parameters: past a patch's bounds a point's distance from the
        // surface jumps; its parametric form carries on.)
        contacts[i] = match newton(start, |x| Ok([sides[i].distance(x)?.0, (x - o).dot(tangent), (x - o).length() - d])) {
            Ok(x) => x,
            Err(e) if on_edge => on_surface_at(&sides[i], start, o, tangent, d).ok_or(e)?,
            Err(e) => return Err(e),
        };
    }
    let middle = Point::barycenter(&[(contacts[0], 1.0), (contacts[1], 1.0)]);
    let center = middle + (middle - o);
    let contact_uv = [sides[0].uv(contacts[0])?, sides[1].uv(contacts[1])?];
    Ok(Section { t: 0.0, center, contacts, contact_uv, radius: distances[0].max(distances[1]), straight: true })
}

/// Where the cut from `from` (on face `first`'s surface, in the plane through `o` square to
/// `tangent`) meets the other face's surface: turned by `angle` from the first face's way
/// back towards `o` there, towards `guess` (a point on the other face near it).
fn cut_at_angle(sides: &[Side; 2], o: Point, tangent: Vector, from: Point, first: usize, angle: f64, guess: Point) -> Option<Point> {
    let n = sides[first].distance(from).ok()?.2;
    let back = n.cross(tangent);
    let back = if back.dot(o - from) >= 0.0 { back } else { -back };
    if back.length() <= 1e-12 {
        return None;
    }
    let back = back.normalized();
    // (Turned about the edge's tangent, either way: the one heading nearer the guess.)
    let turned = |a: f64| back * a.cos() + tangent.cross(back) * a.sin();
    let towards = (guess - from).normalized();
    let way = [turned(angle), turned(-angle)].into_iter().max_by(|x, y| x.dot(towards).total_cmp(&y.dot(towards)))?;
    let other = &sides[1 - first];
    let at = |l: f64| other.distance(from + way * l).map(|d| d.0);
    // Secant steps along the cut from the guess's distance.
    let (mut l0, mut l1) = (0.9 * (guess - from).length(), (guess - from).length());
    let (mut f0, mut f1) = (at(l0).ok()?, at(l1).ok()?);
    for _ in 0..60 {
        if f1.abs() <= 1e-13 * (1.0 + l1) {
            break;
        }
        if f1 == f0 {
            return None;
        }
        let l2 = l1 - f1 * (l1 - l0) / (f1 - f0);
        if !(l2.is_finite() && l2 > 0.0) {
            return None;
        }
        (l0, f0) = (l1, f1);
        l1 = l2;
        f1 = at(l1).ok()?;
    }
    (f1.abs() <= 1e-9 * (1.0 + l1)).then(|| from + way * l1)
}

/// The point of `side`'s surface in the plane through `o` square to `tangent`, `d` from
/// `o`: Newton's method in the surface's parameters, from `start`'s.
fn on_surface_at(side: &Side, start: Point, o: Point, tangent: Vector, d: f64) -> Option<Point> {
    let (mut u, mut v) = side.parameters(start);
    let scale = d.max(1e-9);
    for _ in 0..60 {
        let (p, du, dv) = side.surface.d1(u, v);
        let r = p - o;
        let length = r.length();
        if !(length > 0.0) {
            return None;
        }
        let f = [r.dot(tangent), length - d];
        if f[0].abs() <= 1e-12 * scale && f[1].abs() <= 1e-12 * scale {
            return Some(p);
        }
        let radial = r * (1.0 / length);
        let (a, b, c, e) = (du.dot(tangent), dv.dot(tangent), du.dot(radial), dv.dot(radial));
        let det = a * e - b * c;
        if det.abs() <= 1e-300 {
            return None;
        }
        let (su, sv) = ((f[0] * e - b * f[1]) / det, (a * f[1] - c * f[0]) / det);
        u -= su;
        v -= sv;
        if !(u.is_finite() && v.is_finite()) {
            return None;
        }
    }
    let p = side.surface.point(u, v);
    ((p - o).dot(tangent).abs() <= 1e-8 * scale && ((p - o).length() - d).abs() <= 1e-8 * scale).then_some(p)
}

/// Newton's method on a function of a point with three components, from `start`.
fn newton(start: Point, f: impl Fn(Point) -> Result<[f64; 3]>) -> Result<Point> {
    newton_damped(start, f, false)
}

/// [`newton`], with each step shortened (`damped`) where the whole of it leaves the system
/// further from solved.
fn newton_damped(start: Point, f: impl Fn(Point) -> Result<[f64; 3]>, damped: bool) -> Result<Point> {
    let scale = 1.0 + (start - Point::ORIGIN).length();
    let h = 1e-7 * scale;
    let mut x = start;
    for _ in 0..60 {
        let fx = f(x)?;
        if fx.iter().all(|v| v.abs() <= 1e-12 * scale) {
            return Ok(x);
        }
        let columns: Vec<[f64; 3]> = [Vector::new(h, 0.0, 0.0), Vector::new(0.0, h, 0.0), Vector::new(0.0, 0.0, h)]
            .iter()
            .map(|&d| {
                let g = f(x + d)?;
                Ok([(g[0] - fx[0]) / h, (g[1] - fx[1]) / h, (g[2] - fx[2]) / h])
            })
            .collect::<Result<_>>()?;
        let step = super::solve3([columns[0], columns[1], columns[2]], fx)
            .ok_or_else(|| Error::InvalidInput("no fillet section here (singular system)".into()))?;
        let step = Vector::new(step[0], step[1], step[2]);
        // (Shorter where the whole step leaves the system further from solved, as where the
        // plane of a section turns fast with its centre; the whole step where none helps.)
        if !damped {
            x = x - step;
            continue;
        }
        let norm = |v: [f64; 3]| v.iter().map(|x| x * x).sum::<f64>();
        let size = norm(fx);
        let mut part = 1.0;
        while part > 1e-3 && !f(x - step * part).is_ok_and(|g| norm(g) < size) {
            part *= 0.5;
        }
        x = x - step * if part > 1e-3 { part } else { 1.0 };
    }
    if f(x)?.iter().all(|v| v.abs() <= 1e-8 * scale) {
        Ok(x)
    } else {
        Err(Error::InvalidInput("no fillet section here (no convergence)".into()))
    }
}

/// The fillet surface: across (u in [0, 1]) the arc of each section from the contact on the
/// first face to the one on the second, as a rational quadratic; along the edge (v, the
/// spine's parameter) a cubic interpolation of the arcs' homogeneous control points.
fn fillet_surface(sections: &[Section]) -> Result<Surface> {
    let ts: Vec<f64> = sections.iter().map(|s| s.t).collect();
    if sections.iter().all(|s| s.straight) {
        // A chamfer: ruled between the contact curves (across, u in [0, 1], degree 1).
        let rows = [
            fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[0]).collect::<Vec<_>>(), &ts)?,
            fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[1]).collect::<Vec<_>>(), &ts)?,
        ];
        let u = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let v = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
        return Ok(Surface::BSpline(BSplineSurface::try_new(vec![rows[0].poles().to_vec(), rows[1].poles().to_vec()], u, v)?));
    }
    let mut middle = Vec::with_capacity(sections.len());
    let mut weights = Vec::with_capacity(sections.len());
    for s in sections {
        let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        let w = half.cos();
        // The corner of the arc's control polygon, where its end tangents meet.
        let corner = s.center + (a + b).normalized() * (s.radius / w);
        middle.push(Point::ORIGIN + (corner - Point::ORIGIN) * w);
        weights.push(Point::new(w, 0.0, 0.0));
    }
    let rows = [
        fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[0]).collect::<Vec<_>>(), &ts)?,
        fit::interpolate_with_parameters(&middle, &ts)?,
        fit::interpolate_with_parameters(&sections.iter().map(|s| s.contacts[1]).collect::<Vec<_>>(), &ts)?,
    ];
    let w = fit::interpolate_with_parameters(&weights, &ts)?;
    let w_poles: Vec<f64> = w.poles().iter().map(|p| p.x).collect();
    if w_poles.iter().any(|&w| w <= 0.0) {
        return Err(Error::InvalidInput("the fillet turns too sharply".into()));
    }
    let poles = vec![
        rows[0].poles().to_vec(),
        rows[1].poles().iter().zip(&w_poles).map(|(p, w)| Point::ORIGIN + (*p - Point::ORIGIN) * (1.0 / w)).collect(),
        rows[2].poles().to_vec(),
    ];
    let n = w_poles.len();
    let u = KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2);
    let v = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
    let surface = BSplineSurface::try_new(poles, u, v)?.with_weights(vec![vec![1.0; n], w_poles, vec![1.0; n]]);
    Ok(Surface::BSpline(surface))
}

/// Where the fillet ends at a vertex of the edge: the face there, the curve along which the
/// fillet surface crosses it (from the contact curve on the first face to the one on the
/// second) and the two edges at the vertex it shortens.
struct EndCut {
    vertex: VertexId,
    face: Side,
    /// The crossing as a 3D curve over [0, 1] and as a pcurve on the fillet surface.
    curve: BSplineCurve,
    on_fillet: BSplineCurve2,
    /// The crossing on the end face, over [0, 1].
    on_face: BSplineCurve2,
    /// The fillet's v where it crosses at u = 0 and u = 1.
    v: [f64; 2],
    /// Per face along the edge: the edge at the vertex it shares with the end face, and the
    /// parameter where the contact curve meets it.
    edges: [(EdgeId, f64); 2],
    /// Per face along the edge: whether the contact meets that edge's curve past the vertex
    /// (the edge is kept, and carried on to the contact by a new edge).
    extended: [bool; 2],
    /// An end on the shell's free boundary: no face there (`face` is the first side's, not
    /// rebuilt); the crossing runs straight across the fillet's parameters between where
    /// the contact curves meet the two free edges, and is left free.
    open: bool,
    /// Per face along the edge: the far vertex of the edge at the vertex, where the contact
    /// meets it there (as near as a fraction of the fillet's width): the contact ends at
    /// that vertex and the edge is taken away whole.
    onto: [Option<VertexId>; 2],
    /// Where a contact runs over the end face's edge onto the face beyond: the crossing
    /// (`curve` and the rest) is the end face's piece, from where the two meet.
    beyond: Option<Beyond>,
    /// Whether `face`'s surface is the end face's carried on past its patch's bounds (for
    /// the face rebuilt to take).
    carried: bool,
    /// Where the fillet ends across a fan of faces round the vertex: the crossing (`curve`
    /// and the rest) is the first face's piece, these the others'.
    fan: Option<Fan>,
}

/// The faces of a fan after the first (in the order u runs), the crossing with each, and the
/// edges between each two (from the first) with where the crossing passes them.
struct Fan {
    faces: Vec<Side>,
    pieces: Vec<Piece>,
    splits: Vec<(EdgeId, f64)>,
}

/// The faces round `vertex` from the face across edge `a` from the first side to the face
/// across edge `b` from the second, each sharing an edge at the vertex with the one before
/// (those edges given too); `None` where they do not go round so.
fn fan_round(model: &Model, shape: Shape, vertex: Shape, sides: &[Side; 2], a: EdgeId, b: EdgeId) -> Option<(Vec<Shape>, Vec<EdgeId>)> {
    let faces_of = |e: EdgeId| -> Vec<Shape> { model.explore(shape, NodeKind::Face).filter(|f| model.explore(*f, NodeKind::Edge).any(|x| x.node == crate::topology::NodeId::from(e))).collect() };
    let at_vertex = |e: EdgeId| {
        let (p, q) = model.edge_vertices(Shape::from(e));
        p.node == vertex.node || q.node == vertex.node
    };
    let mut face = faces_of(a).into_iter().find(|f| f.node != sides[0].face.node)?;
    let (mut faces, mut inner) = (vec![face], Vec::new());
    let mut came = a;
    for _ in 0..8 {
        let edges: Vec<EdgeId> = model
            .explore(face, NodeKind::Edge)
            .filter_map(|e| EdgeId::try_from(e.node).ok())
            .filter(|&e| e != came && at_vertex(e) && model.edge_curve(Shape::from(e)).is_some())
            .collect();
        let [next] = edges[..] else { return None };
        if next == b {
            return faces_of(b).iter().any(|f| f.node == sides[1].face.node).then_some((faces, inner));
        }
        face = faces_of(next).into_iter().find(|f| f.node != face.node)?;
        if faces.iter().any(|f| f.node == face.node) {
            return None;
        }
        faces.push(face);
        inner.push(next);
        came = next;
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn end_cut(
    model: &Model,
    shape: Shape,
    id: EdgeId,
    vertex: Shape,
    sides: &[Side; 2],
    fillet: &Surface,
    ts: &[f64],
    at_end: bool,
    target: f64,
) -> Result<EndCut> {
    // (Where the end face grows to reach a contact and its surface is a patch: on its
    // surface carried on past the patch's bounds (the same over them), which the face then
    // takes, as points off the patch are not found on the patch itself. Where no crossing
    // is found, again with each crossing's feet on the end face searched afresh.)
    let on = |carry: bool, fresh: bool| end_cut_on(model, shape, id, vertex, sides, fillet, ts, at_end, target, carry, fresh);
    match on(false, false) {
        Ok(cut) if cut.open || !cut.extended.iter().any(|&x| x) => Ok(cut),
        Ok(cut) => Ok(on(true, false).or_else(|_| on(true, true)).ok().filter(|c| c.carried).unwrap_or(cut)),
        Err(e) => on(false, true).map_err(|_| e),
    }
}

/// Whether `pcurve` over [`lo`, `hi`] leaves the patch of `side`'s face (its surface
/// carried on past it, the face then to take that surface).
fn off_patch(model: &Model, side: &Side, pcurve: &Curve2, lo: f64, hi: f64) -> bool {
    let Some(own) = model.surface(side.face) else { return false };
    if !matches!(own, Surface::Bezier(_) | Surface::BSpline(_)) || own == &side.surface {
        return false;
    }
    let (u0, u1, v0, v1) = own.bounds();
    let (du, dv) = (1e-9 * (u1 - u0), 1e-9 * (v1 - v0));
    (0..=32).any(|k| {
        let q = pcurve.point(lo + (hi - lo) * k as f64 / 32.0);
        q.x < u0 - du || q.x > u1 + du || q.y < v0 - dv || q.y > v1 + dv
    })
}

/// `surface` carried on past each bound by about `length`, if it is a patch (a Bezier or
/// B-spline surface).
fn carried_on(surface: &Surface, length: f64) -> Option<Surface> {
    if !matches!(surface, Surface::Bezier(_) | Surface::BSpline(_)) {
        return None;
    }
    // (Not round a period, which has no bound: carried on across it, it would be closed no
    // longer.)
    let (u_periodic, v_periodic) = (surface.is_u_periodic(), surface.is_v_periodic());
    if u_periodic || v_periodic {
        return None;
    }
    let mut carried = surface.clone();
    for side in [extend::SurfaceSide::UMin, extend::SurfaceSide::UMax, extend::SurfaceSide::VMin, extend::SurfaceSide::VMax] {
        carried = Surface::BSpline(extend::surface_by_length(&carried, length, 2, side)?);
    }
    Some(carried)
}

/// [`end_cut`], on the end face's surface carried on past its patch (`carry`, where it is
/// one), its feet there searched afresh for each crossing (`fresh`).
#[allow(clippy::too_many_arguments)]
fn end_cut_on(
    model: &Model,
    shape: Shape,
    id: EdgeId,
    vertex: Shape,
    sides: &[Side; 2],
    fillet: &Surface,
    ts: &[f64],
    at_end: bool,
    target: f64,
    carry: bool,
    fresh: bool,
) -> Result<EndCut> {
    let v_id = VertexId::try_from(vertex.node).expect("a vertex");
    // The other edges at the vertex: one on each face along the edge.
    let at_vertex = |e: EdgeId| {
        let (a, b) = model.edge_vertices(Shape::from(e));
        a.node == vertex.node || b.node == vertex.node
    };
    // (Not degenerate edges: a cone's apex, say.)
    let face_edges = |face: Shape| -> Vec<EdgeId> {
        model
            .explore(face, NodeKind::Edge)
            .filter_map(|e| EdgeId::try_from(e.node).ok())
            .filter(|&e| e != id && at_vertex(e) && model.edge_curve(Shape::from(e)).is_some())
            .collect()
    };
    let neighbours: Vec<Vec<EdgeId>> = sides.iter().map(|s| face_edges(s.face)).collect();
    let (Some(&a), Some(&b)) = (neighbours[0].first(), neighbours[1].first()) else {
        return Err(Error::InvalidInput("rolling-ball fillets ending where no other edge continues each face are not supported yet".into()));
    };
    if neighbours[0].len() != 1 || neighbours[1].len() != 1 || a == b {
        return Err(Error::InvalidInput("rolling-ball fillets ending at vertices with more than three edges are not supported yet".into()));
    }
    // (Each face going on from the vertex along a free edge: an open end.)
    let uses = |e: EdgeId| model.explore(shape, NodeKind::Face).filter(|f| model.explore(*f, NodeKind::Edge).any(|x| x.node == crate::topology::NodeId::from(e))).count();
    if uses(a) == 1 && uses(b) == 1 {
        return open_end(model, v_id, [a, b], sides, fillet, ts, at_end, target);
    }
    // The end face: the one using both edges.
    let end_faces: Vec<Shape> = model
        .explore(shape, NodeKind::Face)
        .filter(|&f| f.node != sides[0].face.node && f.node != sides[1].face.node)
        .filter(|&f| {
            let edges: Vec<EdgeId> = model.explore(f, NodeKind::Edge).filter_map(|e| EdgeId::try_from(e.node).ok()).collect();
            edges.contains(&a) && edges.contains(&b)
        })
        .collect();
    // (Else the faces round the vertex from the first face's edge to the second's: a fan,
    // the fillet cut by each in turn.)
    let (face, fan) = match end_faces[..] {
        [face] => (face, None),
        _ => {
            let (faces, inner) = fan_round(model, shape, vertex, sides, a, b)
                .ok_or_else(|| Error::InvalidInput("rolling-ball fillets ending at vertices with more than three faces are not supported yet".into()))?;
            (faces[0], Some((faces, inner)))
        }
    };
    let side_of = |face: Shape| Side {
        face,
        surface: model.surface(face).expect("a face's surface").clone(),
        sign: if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 },
        along: 1.0,
        last: Default::default(),
        v_range: { let b = model.face_uv_bounds(face); (b.2, b.3) },
        own_range: Default::default(),
    };
    let fan_sides: Vec<Side> = fan.as_ref().map(|(faces, _)| faces.iter().map(|&f| side_of(f)).collect()).unwrap_or_default();
    let carry = carry && fan.is_none();
    let sign = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    let surface = model.surface(face).expect("a face's surface").clone();
    // (By twice the fillet's width across the middle of the edge.)
    let carried = if carry {
        let v = ts[ts.len() / 2];
        carried_on(&surface, 2.0 * fillet.point(0.0, v).distance(fillet.point(1.0, v)))
    } else {
        None
    };
    let carries = carried.is_some();
    let end = Side { face, surface: carried.unwrap_or(surface), sign, along: 1.0, last: Default::default(), v_range: { let b = model.face_uv_bounds(face); (b.2, b.3) }, own_range: Default::default() };

    // Across the fillet, where it leaves the end face's surface: going from the middle of
    // the edge towards this end, the outward distance turns positive.
    let mid = 0.5 * (ts[0] + ts[ts.len() - 1]);
    let order: Vec<f64> = if at_end {
        ts.iter().copied().filter(|&t| t >= mid).collect()
    } else {
        ts.iter().rev().copied().filter(|&t| t <= mid).collect()
    };
    // (Where the fillet crosses the end face's surface on its way to this end, either way
    // across it: the end face may lie on either side of the fillet. Of several, the one
    // nearest the edge's vertex.)
    let vertex_point = model.point(vertex);
    let crossing_on = |end: &Side, u: f64| -> Result<(f64, f64)> {
        // (Afresh: from a search of the whole surface, then each point from the one before,
        // as a foot found for another crossing may lead to a far one.)
        if fresh {
            end.last.set(None);
        }
        let d = |v: f64| -> Result<f64> { Ok(end.distance(fillet.point(u, v))?.0) };
        // (Into the end face's outside first, as the fillet usually leaves the solid there.)
        let mut best: Option<(f64, f64)> = None;
        for outwards in [true, false] {
        if best.is_some() {
            break;
        }
        // (Each section's distance found once, in order along the edge: the projections onto
        // the end face each start from the foot before, as the crossing's refinement then
        // does from its bracket's end.)
        let mut previous: Option<f64> = None;
        for w in order.windows(2) {
            let da = match previous {
                Some(x) => x,
                None => d(w[0])?,
            };
            let db = d(w[1])?;
            previous = Some(db);
            if (da <= 0.0) != (db <= 0.0) && (da <= 0.0) == outwards {
                let (mut lo, mut hi) = (w[0], w[1]);
                // (Halved until the halves meet in the last bits: no further step moves.)
                for _ in 0..80 {
                    let m = 0.5 * (lo + hi);
                    if m == lo || m == hi {
                        break;
                    }
                    if (d(m)? <= 0.0) == (da <= 0.0) {
                        lo = m;
                    } else {
                        hi = m;
                    }
                }
                let v = 0.5 * (lo + hi);
                // (Outwards the first from the middle; the other way, the one nearest the
                // vertex.)
                if outwards {
                    best = Some((0.0, v));
                    break;
                }
                let gap = fillet.point(u, v).distance(vertex_point);
                if best.is_none_or(|(g, _)| gap < g) {
                    best = Some((gap, v));
                }
            }
        }
        }
        best.map(|(_, v)| (u, v)).ok_or_else(|| Error::InvalidInput(NOT_REACHED.into()))
    };
    let crossing = |u: f64| crossing_on(&end, u);
    // Where the contact curves meet the edges at the vertex.
    // (Where a contact curve would meet the edge past its far end, the fillet is wider than
    // the face there: it would run over onto the next face, which is not built here. Told
    // from the crossing's ends before it is refined.)
    // (On the edge, or on its curve carried on past the vertex: the end face then grows to
    // reach the contact, as where the fillet ends wider than its end face's edge.)
    let last = fan_sides.last().unwrap_or(&end);
    let ends = [crossing(0.0)?, crossing_on(last, 1.0)?];
    let width = fillet.point(ends[0].0, ends[0].1).distance(fillet.point(ends[1].0, ends[1].1));
    // (Where the contact meets the edge at its far vertex, or as good as, a little short
    // or past it: the edge is taken away whole, rather than leave a sliver of it or of the
    // face beyond.)
    let onto = |e: EdgeId, p: Point| -> Option<VertexId> {
        let (a, b) = model.edge_vertices(Shape::from(e));
        let far = if a.node == vertex.node { b } else { a };
        (far.node != vertex.node && model.point(far).distance(p) <= ONTO_VERTEX * width).then(|| VertexId::try_from(far.node).expect("a vertex"))
    };
    let near_to = |p: Point| (1e3 * target).max(1e-6 * (1.0 + (p - Point::ORIGIN).length()));
    let on_edge = |e: EdgeId, p: Point| -> Result<Option<f64>> {
        let (c, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        let t = parameter_on(&c, first, last, p);
        // (Or as good as: the end face's surface and the edge's curve each hold to the faces
        // only so closely.)
        Ok((c.point(t).distance(p) <= near_to(p).max(ON_EDGE * width)).then_some(t))
    };
    // (Where the contact on face `i` passes the far vertex of its edge `e` at the vertex, the
    // fillet runs over the end face's edge from there onto the face beyond: cut by that face
    // from the contact to where it crosses that edge, by the end face from there on.)
    let run_over = |i: usize, e: EdgeId| -> Result<Option<Beyond>> {
        let (ea, eb) = model.edge_vertices(Shape::from(e));
        let far = if ea.node == vertex.node { eb } else { ea };
        if far.node == vertex.node {
            return Ok(None);
        }
        let at_far = |face: Shape| -> Vec<EdgeId> {
            model
                .explore(face, NodeKind::Edge)
                .filter_map(|x| EdgeId::try_from(x.node).ok())
                .filter(|&x| x != e && model.edge_curve(Shape::from(x)).is_some())
                .filter(|&x| {
                    let (a, b) = model.edge_vertices(Shape::from(x));
                    a.node == far.node || b.node == far.node
                })
                .collect()
        };
        let (on_side, on_end) = (at_far(sides[i].face), at_far(end.face));
        let ([cut_side], [cut_end]) = (&on_side[..], &on_end[..]) else { return Ok(None) };
        let beyond_faces: Vec<Shape> = model
            .explore(shape, NodeKind::Face)
            .filter(|&f| f.node != sides[i].face.node && f.node != end.face.node)
            .filter(|&f| {
                let edges: Vec<EdgeId> = model.explore(f, NodeKind::Edge).filter_map(|x| EdgeId::try_from(x.node).ok()).collect();
                edges.contains(cut_side) && edges.contains(cut_end)
            })
            .collect();
        let [g] = beyond_faces[..] else { return Ok(None) };
        let sign = if g.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
        let next = Side { face: g, surface: model.surface(g).expect("a face's surface").clone(), sign, along: 1.0, last: Default::default(), v_range: { let b = model.face_uv_bounds(g); (b.2, b.3) }, own_range: Default::default() };
        // Across the fillet from this contact, the crossing with the end face's surface
        // outside the face beyond, then inside it: where it goes in, the end face's edge.
        let outside = |u: f64| -> Result<f64> {
            let (u, v) = crossing(u)?;
            Ok(next.distance(fillet.point(u, v))?.0)
        };
        let (u_side, u_other) = (i as f64, 1.0 - i as f64);
        if !(outside(u_side)? > 0.0 && outside(u_other)? < 0.0) {
            return Ok(None);
        }
        let (mut lo, mut hi) = (u_side, u_other);
        for _ in 0..80 {
            let m = 0.5 * (lo + hi);
            if outside(m)? > 0.0 { lo = m } else { hi = m }
        }
        let at = 0.5 * (lo + hi);
        let meeting = { let (u, v) = crossing(at)?; fillet.point(u, v) };
        let Some(t_end) = on_edge(*cut_end, meeting)? else { return Ok(None) };
        let contact = { let (u, v) = crossing_on(&next, u_side)?; fillet.point(u, v) };
        let Some(t_side) = on_edge(*cut_side, contact)? else { return Ok(None) };
        Ok(Some(Beyond { side: i, face: next, at, far: VertexId::try_from(far.node).expect("a vertex"), split: (*cut_end, t_end), trimmed: (*cut_side, t_side), removed: e, piece: None }))
    };
    let meet = |e: EdgeId, p: Point| -> Result<(f64, bool)> {
        let (c, first, last) = model.edge_curve(Shape::from(e)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        if onto(e, p).is_some() {
            let (a, _) = model.edge_vertices(Shape::from(e));
            return Ok((if a.node == vertex.node { last } else { first }, false));
        }
        let near = near_to(p);
        if let Some(t) = on_edge(e, p)? {
            return Ok((t, false));
        }
        let (a, _) = model.edge_vertices(Shape::from(e));
        let at = if a.node == vertex.node { first } else { last };
        let beyond = at + (at - if a.node == vertex.node { last } else { first });
        let t = parameter_on(&c, at, beyond, p);
        if c.point(t).distance(p) <= near && (t - at).abs() > 1e-9 * (last - first).abs() {
            return Ok((t, true));
        }
        // (Or along where the side's surface and the end face's meet, from the vertex: the
        // edge's curve does not carry on there, as a cone's generatrix at its apex.)
        if let Some(i) = sides.iter().position(|s| model.explore(s.face, NodeKind::Edge).any(|x| x.node == crate::topology::NodeId::from(e))) {
            if intersection_from(&sides[i], &end, model.point(vertex), p).is_ok() {
                return Ok((f64::NAN, true));
            }
        }
        Err(Error::InvalidInput("a fillet running over onto the next face is not supported yet".into()))
    };
    let (pa, pb) = (fillet.point(ends[0].0, ends[0].1), fillet.point(ends[1].0, ends[1].1));
    let onto = [onto(a, pa), onto(b, pb)];
    // (A run-over where the contact is neither on its edge nor at its far vertex.)
    let mut beyond = None;
    for (i, e, p) in [(0, a, pa), (1, b, pb)] {
        if fan.is_none() && onto[i].is_none() && on_edge(e, p)?.is_none() && beyond.is_none() {
            beyond = run_over(i, e)?;
        }
    }
    let (ma, mb) = match &beyond {
        Some(r) if r.side == 0 => ((f64::NAN, false), meet(b, pb)?),
        Some(_) => (meet(a, pa)?, (f64::NAN, false)),
        None => (meet(a, pa)?, meet(b, pb)?),
    };
    let edges = [(a, ma.0), (b, mb.0)];
    let extended = [ma.1, mb.1];
    // The crossing with `side`'s surface for u from `lo` to `hi`, over [0, 1]: more samples
    // across until the interpolated crossing stays on the face (or doubling them hardly
    // helps).
    let piece = |side: &Side, lo: f64, hi: f64| -> Result<Piece> {
        let mut across = ACROSS;
        let mut previous_worst = f64::INFINITY;
        let (uv, curve, on_fillet) = loop {
            let uv: Vec<(f64, f64)> = (0..=across).map(|j| crossing_on(side, lo + (hi - lo) * j as f64 / across as f64)).collect::<Result<_>>()?;
            let params: Vec<f64> = (0..=across).map(|j| j as f64 / across as f64).collect();
            let points: Vec<Point> = uv.iter().map(|&(u, v)| fillet.point(u, v)).collect();
            let curve = fit::interpolate_with_parameters(&points, &params)?;
            let on_fillet = fit::interpolate_with_parameters(&uv.iter().map(|&(u, v)| Point2::new(u, v)).collect::<Vec<_>>(), &params)?;
            let mut worst: f64 = 0.0;
            for w in params.windows(2) {
                let s = 0.5 * (w[0] + w[1]);
                let p = curve.point(s);
                let q = on_fillet.point(s);
                worst = worst.max(side.distance(p)?.0.abs()).max(fillet.point(q.x, q.y).distance(p));
            }
            let stalled = across >= 4 * ACROSS && worst > 0.25 * previous_worst;
            if worst <= target || stalled || across >= 64 * ACROSS {
                break (uv, curve, on_fillet);
            }
            previous_worst = worst;
            across *= 2;
        };
        let params: Vec<f64> = (0..uv.len()).map(|j| j as f64 / (uv.len() - 1) as f64).collect();
        let face_uv: Vec<Point2> = uv.iter().map(|&(u, v)| side.uv(fillet.point(u, v))).collect::<Result<_>>()?;
        let on_face = fit::interpolate_with_parameters(&unwrapped(face_uv, &side.surface), &params)?;
        Ok(Piece { curve, on_fillet, on_face, v: [uv[0].1, uv[uv.len() - 1].1] })
    };
    if let Some((_, inner)) = &fan {
        // Across the fan: on each face's surface from where it took over from the face before
        // to where the next takes over (the crossing passes the edge between them).
        if extended.iter().any(|&x| x) {
            return Err(Error::InvalidInput("a fillet ending across several faces where a contact misses its edge is not supported yet".into()));
        }
        let mut at = vec![0.0];
        for (k, w) in fan_sides.windows(2).enumerate() {
            let over = |u: f64| -> Result<f64> {
                let (u, v) = crossing_on(&w[0], u)?;
                Ok(w[1].distance(fillet.point(u, v))?.0)
            };
            let from = at[k];
            let n = 64;
            let mut found = None;
            let mut before = over(from + (1.0 - from) * 1e-9)?;
            for j in 1..=n {
                let u = from + (1.0 - from) * j as f64 / n as f64;
                let now = over(u)?;
                if (now > 0.0) != (before > 0.0) {
                    let (mut lo, mut hi) = (from + (1.0 - from) * (j - 1) as f64 / n as f64, u);
                    for _ in 0..80 {
                        let m = 0.5 * (lo + hi);
                        if (over(m)? > 0.0) == (before > 0.0) { lo = m } else { hi = m }
                    }
                    found = Some(0.5 * (lo + hi));
                    break;
                }
                before = now;
            }
            at.push(found.ok_or_else(|| Error::InvalidInput("a fillet's end not found passing from face to face round its vertex".into()))?);
        }
        at.push(1.0);
        let mut pieces = Vec::new();
        for (k, s) in fan_sides.iter().enumerate() {
            pieces.push(piece(s, at[k], at[k + 1])?);
        }
        // Where each passes to the next: on the edge between them.
        let mut splits = Vec::new();
        for (k, &e) in inner.iter().enumerate() {
            let p = pieces[k].curve.point(1.0);
            let t = on_edge(e, p)?.ok_or_else(|| Error::InvalidInput("a fillet's end passing from face to face away from the edge between them".into()))?;
            splits.push((e, t));
        }
        let mut pieces = pieces.into_iter();
        let Piece { curve, on_fillet, on_face, v: v_first } = pieces.next().expect("a piece");
        let rest: Vec<Piece> = pieces.collect();
        let v = [v_first[0], rest.last().map_or(v_first[1], |p| p.v[1])];
        let next = fan_sides[1..].to_vec();
        return Ok(EndCut { vertex: v_id, face: end, curve, on_fillet, on_face, v, edges, extended, open: false, onto, beyond: None, carried: false, fan: Some(Fan { faces: next, pieces: rest, splits }) });
    }
    match beyond {
        None => {
            let Piece { curve, on_fillet, on_face, v } = piece(&end, 0.0, 1.0)?;
            Ok(EndCut { vertex: v_id, face: end, curve, on_fillet, on_face, v, edges, extended, open: false, onto, beyond: None, carried: carries, fan: None })
        }
        Some(mut r) => {
            // (Each piece in the order u runs: the end face's from the meeting on, the face
            // beyond's up to it, or the other way round.)
            let (Piece { curve, on_fillet, on_face, v: v_end }, next) = if r.side == 0 {
                (piece(&end, r.at, 1.0)?, piece(&r.face, 0.0, r.at)?)
            } else {
                (piece(&end, 0.0, r.at)?, piece(&r.face, r.at, 1.0)?)
            };
            let v = if r.side == 0 { [next.v[0], v_end[1]] } else { [v_end[0], next.v[1]] };
            r.piece = Some(next);
            Ok(EndCut { vertex: v_id, face: end, curve, on_fillet, on_face, v, edges, extended, open: false, onto, beyond: Some(r), carried: carries, fan: None })
        }
    }
}

/// A crossing of the fillet with a face's surface, over [0, 1] in the order u runs: as a 3D
/// curve, on the fillet and on the face; and the fillet's v at its ends.
struct Piece {
    curve: BSplineCurve,
    on_fillet: BSplineCurve2,
    on_face: BSplineCurve2,
    v: [f64; 2],
}

/// Where a fillet's contact on face `side` runs over the end face's edge at the far vertex
/// `far` onto the face beyond (`face`): the fillet cut by it up to where it crosses that
/// edge (`split`, at u = `at`), the end face's edge from the vertex (`removed`) taken away,
/// and the face's edge there (`trimmed`) ending at the contact.
struct Beyond {
    side: usize,
    face: Side,
    at: f64,
    far: VertexId,
    split: (EdgeId, f64),
    trimmed: (EdgeId, f64),
    removed: EdgeId,
    /// The crossing with the face beyond.
    piece: Option<Piece>,
}

/// [`end_cut`] where the faces along the edge go on from `vertex` along free edges `free`
/// (one on each face): where each contact curve (the fillet at u = 0 and u = 1) meets its
/// free edge, nearest the vertex, and across the fillet's parameters between those.
#[allow(clippy::too_many_arguments)]
fn open_end(model: &Model, vertex: VertexId, free: [EdgeId; 2], sides: &[Side; 2], fillet: &Surface, ts: &[f64], at_end: bool, target: f64) -> Result<EndCut> {
    let mid = 0.5 * (ts[0] + ts[ts.len() - 1]);
    let far = if at_end { ts[ts.len() - 1] } else { ts[0] };
    let mut v = [0.0; 2];
    let mut edges = [(free[0], 0.0); 2];
    for i in 0..2 {
        let (curve, first, last) = model.edge_curve(Shape::from(free[i])).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
        let u = i as f64;
        let gap = |w: f64| -> (f64, f64) {
            let p = fillet.point(u, w);
            let t = parameter_on(&curve, first, last, p);
            (curve.point(t).distance(p), t)
        };
        // From the middle of the edge towards this end, the first place the contact comes
        // nearest the free edge: sampled, then narrowed.
        let n = 400;
        let ws: Vec<f64> = (0..=n).map(|k| mid + (far - mid) * k as f64 / n as f64).collect();
        let gaps: Vec<f64> = ws.iter().map(|&w| gap(w).0).collect();
        let k = (0..=n).min_by(|&x, &y| gaps[x].total_cmp(&gaps[y])).expect("samples");
        let (mut lo, mut hi) = (ws[k.saturating_sub(1)], ws[(k + 1).min(n)]);
        for _ in 0..100 {
            let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
            if gap(a).0 <= gap(b).0 { hi = b } else { lo = a }
        }
        let w = 0.5 * (lo + hi);
        let (d, t) = gap(w);
        let size = fillet.point(0.0, w).distance(fillet.point(1.0, w)).max(1e-9);
        // (Within the free edge's tolerance, too: its curve strays from its faces that much.)
        if d > (1e3 * target).max(1e-4 * size).max(2.0 * model.tolerance(Shape::from(free[i]))) {
            return Err(Error::InvalidInput("a fillet's contact does not meet the free edge at its end".into()));
        }
        v[i] = w;
        edges[i] = (free[i], t);
    }
    // The crossing: straight in the fillet's parameters, from (0, v0) to (1, v1).
    let params: Vec<f64> = (0..=ACROSS).map(|j| j as f64 / ACROSS as f64).collect();
    let uv: Vec<Point2> = params.iter().map(|&u| Point2::new(u, v[0] + (v[1] - v[0]) * u)).collect();
    let points: Vec<Point> = uv.iter().map(|q| fillet.point(q.x, q.y)).collect();
    let curve = fit::interpolate_with_parameters(&points, &params)?;
    let on_fillet = BSplineCurve2::new(vec![uv[0], uv[ACROSS]], vec![0.0, 1.0], vec![2, 2], 1);
    Ok(EndCut { vertex, face: sides[0].clone(), curve, on_fillet: on_fillet.clone(), on_face: on_fillet, v, edges, extended: [false; 2], open: true, onto: [None; 2], beyond: None, carried: false, fan: None })
}

/// `result`, where edges were carried on past a fillet's end (`extended`), with each joined
/// to the edge it carries on: one edge between the same two faces, without a vertex of
/// theirs alone between (which a fillet along them later could not get past).
fn joined_extensions(model: &mut Model, result: Shape, extended: bool) -> Result<Shape> {
    if !extended {
        return Ok(result);
    }
    let options = crate::unify::UnifyOptions { faces: false, edges: true, same_surface_only: false };
    crate::unify::unify_same_domain_with(model, result, options)
}

/// The new edges of `result` (not `shape`'s) at the vertices `snapped` (contacts ended on
/// them, a little off): no tighter than the vertex, which reaches the contacts' ends.
fn reach_snapped(model: &mut Model, shape: Shape, result: Shape, snapped: &[VertexId]) -> Result<()> {
    if snapped.is_empty() {
        return Ok(());
    }
    let old: std::collections::HashSet<crate::topology::NodeId> = model.unique_subshapes(shape, NodeKind::Edge).iter().map(|e| e.node).collect();
    for e in model.unique_subshapes(result, NodeKind::Edge) {
        if old.contains(&e.node) {
            continue;
        }
        let (a, b) = model.edge_vertices(e);
        let reach = [a, b].iter().filter(|v| snapped.iter().any(|&x| crate::topology::NodeId::from(x) == v.node)).map(|&v| model.tolerance(v)).fold(0.0, f64::max);
        if reach > model.tolerance(e) {
            model.set_edge_tolerance(EdgeId::try_from(e.node).expect("an edge"), reach)?;
        }
    }
    Ok(())
}

/// `vertex`'s tolerance widened to reach `curve` at `t` (an edge shortened to end there).
fn reach_curve(model: &mut Model, vertex: VertexId, curve: &Curve, t: f64) -> Result<()> {
    let gap = curve.point(t).distance(model.point(vertex));
    if gap > model.tolerance(vertex) {
        model.set_vertex_tolerance(vertex, 1.01 * gap)?;
    }
    Ok(())
}

/// The vertex where contact `i` ends at `p`, at an end cut `cut` (if any): the vertex the
/// contact ends on there, its tolerance reaching `p`, else a new one.
fn contact_vertex(model: &mut Model, cut: Option<&EndCut>, i: usize, p: Point, tolerance: f64) -> Result<VertexId> {
    match cut.and_then(|cut| cut.onto[i]) {
        Some(x) => {
            let reach = model.tolerance(x).max(tolerance).max(1.01 * model.point(x).distance(p));
            model.set_vertex_tolerance(x, reach)?;
            Ok(x)
        }
        None => Ok(model.add_vertex(p, tolerance)),
    }
}

/// The edge carrying edge `i` of `cut` (along the fillet's face `side`) on past the vertex to
/// the contact's vertex `to`: over [0, 1] from the vertex, with its pcurves on `side`'s face
/// and on the end face, and its vertices.
fn extension(model: &mut Model, cut: &EndCut, i: usize, side: &Side, to: VertexId, tolerance: f64) -> Result<(EdgeId, Curve2, Curve2, [VertexId; 2])> {
    let (old, t) = cut.edges[i];
    let (c, first, last) = model.edge_curve(Shape::from(old)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet's end".into()))?;
    let (a, _) = model.edge_vertices(Shape::from(old));
    let at = if a.node == crate::topology::NodeId::from(cut.vertex) { first } else { last };
    let n = 16;
    let ws: Vec<f64> = (0..=n).map(|k| k as f64 / n as f64).collect();
    let points: Vec<Point> = if t.is_nan() {
        // (Along where the side's surface and the end face's meet.)
        intersection_from(side, &cut.face, model.point(Shape::from(cut.vertex)), model.point(Shape::from(to)))?
    } else {
        ws.iter().map(|&w| c.point(at + (t - at) * w)).collect()
    };
    let curve = fit::interpolate_with_parameters(&points, &ws)?;
    let e = model.add_edge(cut.vertex, to, tolerance)?;
    model.set_edge_curve(e.into(), Curve::BSpline(curve.clone()), tolerance)?;
    model.set_edge_range(e.into(), 0.0, 1.0)?;
    let on = |s: &Side| -> Result<Curve2> {
        if let Some(p) = project::analytic_pcurve(&Curve::BSpline(curve.clone()), &s.surface) {
            return Ok(p);
        }
        // (From a search of the whole surface, then each from the one before.)
        s.last.set(None);
        let uv: Vec<Point2> = points.iter().map(|&p| s.uv(p)).collect::<Result<_>>()?;
        Ok(Curve2::BSpline(fit::interpolate_with_parameters(&unwrapped(uv, &s.surface), &ws)?))
    };
    let (on_side, on_end) = (on(side)?, on(&cut.face)?);
    // (As loose as its curve and pcurves are apart at its ends, and its vertices with it: the
    // carried curve reaches the contact as nearly as it was found on it.)
    let mut gap: f64 = 0.0;
    for (w, v) in [(0.0, cut.vertex), (1.0, to)] {
        let q = model.point(Shape::from(v));
        let mut far = curve.point(w).distance(q);
        for (pcurve, s) in [(&on_side, side), (&on_end, &cut.face)] {
            let uv = pcurve.point(w);
            far = far.max(s.surface.point(uv.x, uv.y).distance(q));
        }
        gap = gap.max(far);
        if far > model.tolerance(Shape::from(v)) {
            model.set_vertex_tolerance(v, 1.01 * far)?;
        }
    }
    if gap > tolerance {
        model.set_edge_tolerance(e, 1.01 * gap)?;
    }
    Ok((e, on_side, on_end, [cut.vertex, to]))
}

/// Raises `edge`'s tolerance (and its vertices') to cover how far `surface` at `pcurve`
/// strays from `curve` over `[first, last]` (a pcurve interpolated where the surface's
/// parameters crowd, as near a pole).
fn cover_pcurve(model: &mut Model, edge: EdgeId, curve: &Curve, pcurve: &Curve2, surface: &Surface, first: f64, last: f64) -> Result<()> {
    let mut worst: f64 = 0.0;
    for k in 0..=64 {
        let t = first + (last - first) * k as f64 / 64.0;
        let q = pcurve.point(t);
        worst = worst.max(surface.point(q.x, q.y).distance(curve.point(t)));
    }
    let needed = 1.5 * worst;
    if needed > model.tolerance(Shape::from(edge)) {
        model.set_edge_tolerance(edge, needed)?;
        let (a, b) = model.edge_vertices(Shape::from(edge));
        for v in [a, b] {
            if needed > model.tolerance(v) {
                model.set_vertex_tolerance(VertexId::try_from(v.node).expect("a vertex"), needed)?;
            }
        }
    }
    Ok(())
}

/// Points (17, evenly from `from` to `to`) where the surfaces of `a` and `b` meet: each on
/// both, in the plane square to the chord through the chord's point.
fn intersection_from(a: &Side, b: &Side, from: Point, to: Point) -> Result<Vec<Point>> {
    let chord = to - from;
    if chord.length() <= 1e-12 {
        return Err(Error::InvalidInput("no way along the faces to a fillet's end".into()));
    }
    let tangent = chord.normalized();
    let sides = [a.clone(), b.clone()];
    let n = 16;
    let mut points = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let guess = from + chord * (k as f64 / n as f64);
        let p = if k == 0 {
            from
        } else if k == n {
            to
        } else {
            onto_both(&sides, guess, tangent)?
        };
        // (Close to the chord's point: the same branch of the intersection.)
        if p.distance(guess) > 0.5 * chord.length() {
            return Err(Error::InvalidInput("no way along the faces to a fillet's end".into()));
        }
        points.push(p);
    }
    Ok(points)
}

/// The parameter of the point of `curve` over `[first, last]` nearest to `p`.
fn parameter_on(curve: &Curve, first: f64, last: f64, p: Point) -> f64 {
    let n = 200;
    let gap = |t: f64| curve.point(t).distance(p);
    let mut best = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).map(|t| (t, gap(t))).min_by(|a, b| a.1.total_cmp(&b.1)).expect("samples").0;
    let step = (last - first) / n as f64;
    let (mut lo, mut hi) = ((best - step).max(first.min(last)), (best + step).min(first.max(last)));
    let ratio = 0.618_033_988_749_895;
    for _ in 0..100 {
        let (x1, x2) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
        if gap(x1) < gap(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    best = 0.5 * (lo + hi);
    best
}

/// Rebuilds the solid with the fillet.
struct Rebuild<'a> {
    model: &'a mut Model,
    id: EdgeId,
    sides: &'a [Side; 2],
    surface: &'a Surface,
    contacts: &'a [BSplineCurve; 2],
    ends: &'a [EndCut],
    ts: &'a [f64],
}

/// A coedge to build: edge, orientation, pcurve and range.
type Use = (EdgeId, Orientation, Curve2, f64, f64);

impl Rebuild<'_> {
    fn build(&mut self, shape: Shape, closed: bool, sections: &[Section], convex: bool) -> Result<Shape> {
        // Within twice the contact curves' distance from their faces, and no finer than
        // the model's size allows pcurves to be fitted (1e-8 of it).
        let size = {
            let (a, b) = (self.ts[0], *self.ts.last().expect("sections"));
            let (p, q) = (self.contacts[0].point(a), self.contacts[0].point(b));
            p.distance(q) + (p - Point::ORIGIN).length() + 1.0
        };
        // The fillet's v range along each contact curve.
        let (v_start, v_end) = match self.ends {
            [start, end] => ([start.v[0], start.v[1]], [end.v[0], end.v[1]]),
            _ => ([self.ts[0]; 2], [*self.ts.last().expect("sections"); 2]),
        };
        let tolerance = (self.deviation(v_start, v_end)? * 2.0).max(crate::tolerance::LINEAR).max(1e-8 * size);
        // Vertices at the ends of the contact curves.
        let corner = |model: &mut Model, k: usize, i: usize, v: f64| contact_vertex(model, self.ends.get(k), i, self.contacts[i].point(v), tolerance);
        let starts = [corner(self.model, 0, 0, v_start[0])?, corner(self.model, 0, 1, v_start[1])?];
        let finishes = if closed { starts } else { [corner(self.model, 1, 0, v_end[0])?, corner(self.model, 1, 1, v_end[1])?] };
        // Along a closed edge, where each contact curve crosses the other edges of its face
        // at the edge's vertex (a seam): its vertex there.
        let vertex_id = VertexId::try_from(self.model.edge_vertices(Shape::from(self.id)).0.node).expect("a vertex");
        let at_vertex: [Vec<EdgeId>; 2] = [0, 1].map(|i| {
            self.model
                .explore(self.sides[i].face, NodeKind::Edge)
                .filter_map(|e| EdgeId::try_from(e.node).ok())
                .filter(|&e| e != self.id)
                .filter(|&e| {
                    let (a, b) = self.model.edge_vertices(Shape::from(e));
                    let v: crate::topology::NodeId = vertex_id.into();
                    a.node == v || b.node == v
                })
                .collect()
        });
        let mut splits: [Option<(f64, VertexId)>; 2] = [None, None];
        if closed {
            for i in 0..2 {
                let Some(&other) = at_vertex[i].first() else { continue };
                let Some((curve, first, last)) = self.model.edge_curve(Shape::from(other)) else { continue };
                let gap = |v: f64| crate::geometry::extrema::nearest_on_curve(self.contacts[i].point(v), &curve, first, last).distance;
                let n = 400;
                let (mut best, mut at) = (f64::INFINITY, v_start[i]);
                for k in 0..=n {
                    let v = v_start[i] + (v_end[i] - v_start[i]) * k as f64 / n as f64;
                    let d = gap(v);
                    if d < best {
                        (best, at) = (d, v);
                    }
                }
                let step = (v_end[i] - v_start[i]) / n as f64;
                let (mut lo, mut hi) = (at - step, at + step);
                for _ in 0..80 {
                    let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
                    if gap(a) < gap(b) { hi = b } else { lo = a }
                }
                let v = 0.5 * (lo + hi);
                let span = (v_end[i] - v_start[i]).abs();
                if v > v_start[i] + 1e-6 * span && v < v_end[i] - 1e-6 * span {
                    splits[i] = Some((v, self.model.add_vertex(self.contacts[i].point(v), tolerance)));
                }
            }
        }
        // The contact curves as edges (in pieces where split), in order along them.
        let mut contact_edges: [Vec<(EdgeId, f64, f64)>; 2] = [Vec::new(), Vec::new()];
        for i in 0..2 {
            let pieces = match splits[i] {
                Some((v, x)) => vec![(starts[i], x, v_start[i], v), (x, finishes[i], v, v_end[i])],
                None => vec![(starts[i], finishes[i], v_start[i], v_end[i])],
            };
            for (a, b, lo, hi) in pieces {
                let e = self.model.add_edge(a, b, tolerance)?;
                self.model.set_edge_curve(e.into(), Curve::BSpline(self.contacts[i].clone()), tolerance)?;
                self.model.set_edge_range(e.into(), lo, hi)?;
                contact_edges[i].push((e, lo, hi));
            }
        }
        // The fillet's boundary across at each end: the end cut, or for a closed edge the arc
        // of the first section (a seam).
        let mut across = Vec::new();
        let mut shortened: Vec<(EdgeId, EdgeId, f64, f64)> = Vec::new();
        // (face, edge, pcurve, its vertices) carrying edges on past a fillet's end
        let mut side_inserts: Vec<(Shape, (EdgeId, Curve2, [VertexId; 2]))> = Vec::new();
        let mut end_inserts: Vec<(Shape, (EdgeId, Curve2, [VertexId; 2]))> = Vec::new();
        // (edge, its vertex at a fillet's end, the vertex replacing it, the parameter there)
        let mut trims: Vec<(EdgeId, VertexId, VertexId, f64)> = Vec::new();
        // (edges the fillet takes away whole, its contact ending at their far vertex or
        // running past it)
        let mut removed: Vec<EdgeId> = Vec::new();
        // (per end cut, the vertices of its piece on the end face; and the pieces across the
        // faces beyond: the end, the side running over, the edge and its pcurve on the fillet)
        let mut across_ends: Vec<[VertexId; 2]> = Vec::new();
        let mut beyond_across: Vec<(usize, usize, EdgeId, Curve2)> = Vec::new();
        if closed {
            let s = &sections[0];
            let arc = arc_curve(s);
            let e = self.model.add_edge(starts[0], starts[1], tolerance)?;
            self.model.set_edge_curve(e.into(), arc, tolerance)?;
            self.model.set_edge_range(e.into(), 0.0, 1.0)?;
            across.push(e);
            // Other edges of the faces at the edge's vertex (a seam) now end on the contact
            // curves.
            let vertex = self.model.edge_vertices(Shape::from(self.id)).0;
            let v_id = VertexId::try_from(vertex.node).expect("a vertex");
            for i in 0..2 {
                let others: Vec<EdgeId> = self
                    .model
                    .explore(self.sides[i].face, NodeKind::Edge)
                    .filter_map(|e| EdgeId::try_from(e.node).ok())
                    .filter(|&e| e != self.id && !shortened.iter().any(|s| s.0 == e))
                    .collect();
                for old in others {
                    let (a, b) = self.model.edge_vertices(Shape::from(old));
                    let (a_id, b_id) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
                    if a_id != v_id && b_id != v_id {
                        continue;
                    }
                    let (curve, first, last) = self.model.edge_curve(Shape::from(old)).ok_or_else(|| Error::InvalidInput("a degenerate edge at a fillet".into()))?;
                    // Ending where the contact curve crosses it.
                    let (meet, at) = match splits[i] {
                        Some((v, x)) => (x, self.contacts[i].point(v)),
                        None => (starts[i], self.contacts[i].point(v_start[i])),
                    };
                    let t = parameter_on(&curve, first, last, at);
                    let (from, to, lo, hi) = if a_id == v_id { (meet, b_id, t, last) } else { (a_id, meet, first, t) };
                    // (No closer than the edge it replaces holds to its faces.)
                    let kept = tolerance.max(self.model.tolerance(Shape::from(old)));
                    let e = self.model.add_edge(from, to, kept)?;
                    self.model.set_edge_curve(e.into(), curve, kept)?;
                    self.model.set_edge_range(e.into(), lo, hi)?;
                    shortened.push((old, e, lo, hi));
                }
            }
        } else {
            for (k, cut) in self.ends.iter().enumerate() {
                let k_end = k;
                let vs = if k == 0 { starts } else { finishes };
                // (Where the fillet runs over onto the face beyond, the end face's piece runs
                // from where the two meet, the face beyond's piece to it, in the order u runs.)
                let ends = match (&cut.beyond, &cut.fan) {
                    (None, Some(fan)) => {
                        // Across a fan: a vertex where each piece passes to the next (on the edge
                        // between their faces, which now ends there), each piece but the first an
                        // edge of its own, in its face.
                        let splits: Vec<VertexId> = fan.splits.iter().enumerate().map(|(k, &(edge, at))| {
                            let point = if k == 0 { cut.curve.point(1.0) } else { fan.pieces[k - 1].curve.point(1.0) };
                            let s = self.model.add_vertex(point, tolerance);
                            trims.push((edge, cut.vertex, s, at));
                            s
                        }).collect();
                        for (k, (piece, side)) in fan.pieces.iter().zip(&fan.faces).enumerate() {
                            let from = splits[k];
                            let to = if k + 1 < splits.len() { splits[k + 1] } else { vs[1] };
                            let e = self.model.add_edge(from, to, tolerance)?;
                            self.model.set_edge_curve(e.into(), Curve::BSpline(piece.curve.clone()), tolerance)?;
                            self.model.set_edge_range(e.into(), 0.0, 1.0)?;
                            let pcurve = project::analytic_pcurve(&Curve::BSpline(piece.curve.clone()), &side.surface).unwrap_or_else(|| Curve2::BSpline(piece.on_face.clone()));
                            end_inserts.push((side.face, (e, pcurve, [from, to])));
                            beyond_across.push((k_end, 1, e, Curve2::BSpline(piece.on_fillet.clone())));
                        }
                        [vs[0], splits[0]]
                    }
                    (None, None) => vs,
                    (Some(r), _) => {
                        let next = r.piece.as_ref().expect("the piece beyond");
                        let meeting = cut.curve.point(if r.side == 0 { 0.0 } else { 1.0 });
                        let s = self.model.add_vertex(meeting, tolerance);
                        let ends = if r.side == 0 { [s, vs[1]] } else { [vs[0], s] };
                        let next_ends = if r.side == 0 { [vs[0], s] } else { [s, vs[1]] };
                        let e = self.model.add_edge(next_ends[0], next_ends[1], tolerance)?;
                        self.model.set_edge_curve(e.into(), Curve::BSpline(next.curve.clone()), tolerance)?;
                        self.model.set_edge_range(e.into(), 0.0, 1.0)?;
                        let pcurve = project::analytic_pcurve(&Curve::BSpline(next.curve.clone()), &r.face.surface).unwrap_or_else(|| Curve2::BSpline(next.on_face.clone()));
                        end_inserts.push((r.face.face, (e, pcurve, next_ends)));
                        beyond_across.push((k, r.side, e, Curve2::BSpline(next.on_fillet.clone())));
                        // The end face's edge from the vertex goes; the face's edge there and
                        // the end face's edge onto it now end at the contact and the meeting.
                        removed.push(r.removed);
                        trims.push((r.trimmed.0, r.far, vs[r.side], r.trimmed.1));
                        trims.push((r.split.0, r.far, s, r.split.1));
                        for (edge, at, x) in [(r.trimmed.0, r.trimmed.1, vs[r.side]), (r.split.0, r.split.1, s)] {
                            if let Some((c, ..)) = self.model.edge_curve(Shape::from(edge)) {
                                let reach = self.model.tolerance(x).max(1.01 * c.point(at).distance(self.model.point(x)));
                                self.model.set_vertex_tolerance(x, reach)?;
                            }
                        }
                        ends
                    }
                };
                let e = self.model.add_edge(ends[0], ends[1], tolerance)?;
                self.model.set_edge_curve(e.into(), Curve::BSpline(cut.curve.clone()), tolerance)?;
                self.model.set_edge_range(e.into(), 0.0, 1.0)?;
                across.push(e);
                across_ends.push(ends);
                // The edges at the vertex now end where the contact curves meet them (or are
                // carried on to them).
                for (i, &(old, t)) in cut.edges.iter().enumerate() {
                    if cut.beyond.as_ref().is_some_and(|r| r.side == i) {
                        continue;
                    }
                    if cut.onto[i].is_some() {
                        removed.push(old);
                    } else if cut.extended[i] {
                        let (e, on_side, on_end, vs2) = extension(self.model, cut, i, &self.sides[i], vs[i], tolerance)?;
                        side_inserts.push((self.sides[i].face, (e, on_side, vs2)));
                        end_inserts.push((cut.face.face, (e, on_end, vs2)));
                    } else {
                        trims.push((old, cut.vertex, vs[i], t));
                    }
                }
            }
            // Each trimmed edge once, at one end or both.
            let mut done: Vec<EdgeId> = Vec::new();
            for &(old, ..) in &trims {
                if done.contains(&old) || removed.contains(&old) {
                    continue;
                }
                done.push(old);
                let (a, b) = self.model.edge_vertices(Shape::from(old));
                let (curve, first, last) = self.model.edge_curve(Shape::from(old)).expect("a curve");
                let (mut from, mut to) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
                let (mut lo, mut hi) = (first, last);
                let (a_id, b_id) = (from, to);
                for &(_, at, new_vertex, t) in trims.iter().filter(|c| c.0 == old) {
                    if at == a_id {
                        (from, lo) = (new_vertex, t);
                    } else if at == b_id {
                        (to, hi) = (new_vertex, t);
                    }
                    reach_curve(self.model, new_vertex, &curve, t)?;
                }
                // (No closer than the edge it replaces holds to its faces.)
                let kept = tolerance.max(self.model.tolerance(Shape::from(old)));
                let e = self.model.add_edge(from, to, kept)?;
                self.model.set_edge_curve(e.into(), curve, kept)?;
                self.model.set_edge_range(e.into(), lo, hi)?;
                shortened.push((old, e, lo, hi));
            }
        }

        // The faces along the edge, with the contact curves for the edge, and the end faces,
        // with the end cuts.
        let mut rebuilt: Vec<(Shape, Shape)> = Vec::new();
        for (i, side) in self.sides.iter().enumerate() {
            let pcurve = match project::analytic_pcurve(&Curve::BSpline(self.contacts[i].clone()), &side.surface) {
                Some(p) => p,
                None => Curve2::BSpline(fit::interpolate_with_parameters(
                    &unwrapped(sections.iter().map(|s| s.contact_uv[i]).collect(), &side.surface),
                    self.ts,
                )?),
            };
            for &(e, lo, hi) in &contact_edges[i] {
                cover_pcurve(self.model, e, &Curve::BSpline(self.contacts[i].clone()), &pcurve, &side.surface, lo, hi)?;
            }
            let extra: Vec<(EdgeId, Curve2, [VertexId; 2])> = side_inserts.iter().filter(|(f, _)| f.node == side.face.node).map(|(_, x)| x.clone()).collect();
            let (lo, hi) = (v_start[i].min(v_end[i]), v_start[i].max(v_end[i]));
            let carried = off_patch(self.model, side, &pcurve, lo, hi).then_some(&side.surface);
            let face = rebuild_face(self.model, side.face, &[(self.id, &contact_edges[i], pcurve)], &shortened, &removed, &extra, carried)?;
            rebuilt.push((side.face, face));
        }
        // Each end face once, with the cuts of both ends where it is at both.
        let mut end_faces: Vec<(Shape, Vec<(EdgeId, Curve2, [VertexId; 2])>)> = Vec::new();
        for (k, cut) in self.ends.iter().enumerate() {
            if cut.open {
                continue;
            }
            let pcurve = match project::analytic_pcurve(&Curve::BSpline(cut.curve.clone()), &cut.face.surface) {
                Some(p) => p,
                None => Curve2::BSpline(cut.on_face.clone()),
            };
            let insert = (across[k], pcurve, across_ends[k]);
            match end_faces.iter_mut().find(|(f, _)| f.node == cut.face.face.node) {
                Some((_, inserts)) => inserts.push(insert),
                None => end_faces.push((cut.face.face, vec![insert])),
            }
        }
        for (face, insert) in end_inserts {
            match end_faces.iter_mut().find(|(f, _)| f.node == face.node) {
                Some((_, inserts)) => inserts.push(insert),
                None => end_faces.push((face, vec![insert])),
            }
        }
        for (face, inserts) in end_faces {
            if self.sides.iter().any(|s| s.face.node == face.node) {
                todo!("rolling-ball fillets ending on one of their own faces");
            }
            let carried = self.ends.iter().find(|cut| cut.carried && cut.face.face.node == face.node).map(|cut| &cut.face.surface);
            let new = rebuild_face(self.model, face, &[], &shortened, &removed, &inserts, carried)?;
            rebuilt.push((face, new));
        }

        // The fillet face: across at the start (u from 0 to 1), along the second contact,
        // back across at the end, back along the first contact: counter-clockwise in (u, v).
        let line = |o: (f64, f64), d: (f64, f64)| Curve2::Line(Line2::new(Point2::new(o.0, o.1), Direction2::new(d.0, d.1)));
        let across_pcurve = |k: usize| -> Curve2 {
            if closed {
                line((0.0, if k == 0 { v_start[0] } else { v_end[0] }), (1.0, 0.0))
            } else {
                Curve2::BSpline(self.ends[k].on_fillet.clone())
            }
        };
        let last_across = if closed { across[0] } else { across[1] };
        // (Each end's pieces in the order u runs.)
        let pieces = |k: usize, edge: EdgeId| -> Vec<(EdgeId, Curve2)> {
            let mut pieces = vec![(edge, across_pcurve(k))];
            for (_, side, e, p) in beyond_across.iter().filter(|b| b.0 == k) {
                if *side == 0 { pieces.insert(0, (*e, p.clone())) } else { pieces.push((*e, p.clone())) }
            }
            pieces
        };
        let mut uses: Vec<Use> = pieces(0, across[0]).into_iter().map(|(e, p)| (e, Orientation::Forward, p, 0.0, 1.0)).collect();
        for &(e, lo, hi) in &contact_edges[1] {
            uses.push((e, Orientation::Forward, line((1.0, 0.0), (0.0, 1.0)), lo, hi));
        }
        for (e, p) in pieces(1, last_across).into_iter().rev() {
            uses.push((e, Orientation::Reversed, p, 0.0, 1.0));
        }
        for &(e, lo, hi) in contact_edges[0].iter().rev() {
            uses.push((e, Orientation::Reversed, line((0.0, 0.0), (0.0, 1.0)), lo, hi));
        }
        let fillet = new_face(self.model, self.surface.clone(), vec![(uses, Orientation::Forward)], tolerance)?;
        // Outward: away from the ball's centre for a convex edge (the ball is in the solid).
        let mid = &sections[sections.len() / 2];
        let (p, du, dv) = self.surface.d1(0.5, mid.t);
        let away = (p - mid.center).dot(du.cross(dv)) > 0.0;
        let fillet = if away == convex { fillet } else { fillet.reversed() };

        // The solid with the faces replaced and the fillet added.
        let mut solids = Vec::new();
        // (Or its shells, where it has no solids: an open shell, filleted.)
        let kind = if self.model.explore(shape, NodeKind::Solid).next().is_some() { NodeKind::Solid } else { NodeKind::Shell };
        for solid in self.model.explore(shape, kind).collect::<Vec<_>>() {
            let faces: Vec<Shape> = self.model.explore(solid, NodeKind::Face).collect();
            let touched = faces.iter().any(|f| rebuilt.iter().any(|(old, _)| old.node == f.node));
            if !touched {
                solids.push(solid);
                continue;
            }
            let shell = self.model.new_empty(NodeKind::Shell);
            for f in faces {
                match rebuilt.iter().find(|(old, _)| old.node == f.node) {
                    Some(&(_, new)) => {
                        self.model.add(shell, Shape { orientation: f.orientation, ..new })?;
                        if new.node == rebuilt[1].1.node {
                            self.model.add(shell, fillet)?;
                        }
                    }
                    None => self.model.add(shell, f)?,
                }
            }
            if kind == NodeKind::Shell {
                solids.push(shell);
                continue;
            }
            let new_solid = self.model.new_empty(NodeKind::Solid);
            self.model.add(new_solid, shell)?;
            solids.push(new_solid);
        }
        let compound = self.model.new_empty(NodeKind::Compound);
        for s in solids {
            self.model.add(compound, s)?;
        }
        let snapped: Vec<VertexId> = self.ends.iter().flat_map(|cut| cut.onto.into_iter().flatten()).collect();
        reach_snapped(self.model, shape, compound, &snapped)?;
        let extended = self.ends.iter().any(|cut| cut.extended.iter().any(|&x| x));
        joined_extensions(self.model, compound, extended)
    }

    /// The largest distance of the contact curves from their faces' surfaces.
    /// (Between `from` and `to` along each: past the fillet's ends the contact curves run
    /// on through sections past the edge's ends, which may stray from the faces' surfaces.)
    fn deviation(&self, from: [f64; 2], to: [f64; 2]) -> Result<f64> {
        let mut worst: f64 = 0.0;
        for i in 0..2 {
            let (a, b) = (from[i], to[i]);
            for k in 0..=200 {
                let p = self.contacts[i].point(a + (b - a) * k as f64 / 200.0);
                worst = worst.max(self.sides[i].distance(p)?.0.abs());
            }
        }
        Ok(worst)
    }

}

/// A copy of `face` (same surface, forward) whose boundary has `replace`'s edge for the
/// filleted one, the shortened edges for the old ones, without the `removed` ones, and
/// `insert`'s edge between the two shortened edges (which end at `insert`'s vertices); on
/// its surface `carried` on past its patch, if given.
fn rebuild_face(
    model: &mut Model,
    face: Shape,
    replaces: &[(EdgeId, &[(EdgeId, f64, f64)], Curve2)],
    shortened: &[(EdgeId, EdgeId, f64, f64)],
    removed: &[EdgeId],
    inserts: &[(EdgeId, Curve2, [VertexId; 2])],
    carried: Option<&Surface>,
) -> Result<Shape> {
    let forward = Shape { orientation: Orientation::Forward, ..face };
    let surface = match carried {
        Some(s) => s.clone(),
        None => model.surface(forward).expect("a face's surface").clone(),
    };
    let tolerance = model.tolerance(forward);
    let mut wires = Vec::new();
    for wire in model.children(forward).collect::<Vec<_>>() {
        let mut uses: Vec<Use> = Vec::new();
        let mut coedges: Vec<Shape> = model.children_local(wire).collect();
        // (Starting at a coedge kept, for the replacing pieces to follow one.)
        let replaced = |model: &Model, c: &Shape| replaces.iter().any(|r| r.0 == model[CoEdgeId::try_from(c.node).expect("a wire holds coedges")].edge);
        if let Some(k) = coedges.iter().position(|c| !replaced(model, c)) {
            coedges.rotate_left(k);
        }
        for coedge in coedges {
            let c = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
            let (edge, orientation) = (model[c].edge, model[c].orientation);
            let (pcurve, first, last) = model.effective_pcurve(c).ok_or_else(|| Error::InvalidInput("a face without pcurves".into()))?;
            if let Some((_, pieces, p)) = replaces.iter().find(|r| r.0 == edge) {
                {
                    // In the same periods as the edge it replaces; its pieces in the order
                    // the wire runs.
                    let (a, b) = (pieces[0].1, pieces[pieces.len() - 1].2);
                    let target = pcurve.point(0.5 * (first + last));
                    let shifted = align(p, 0.5 * (a + b), target, &surface);
                    let mut ordered: Vec<&(EdgeId, f64, f64)> = if orientation == Orientation::Reversed { pieces.iter().rev().collect() } else { pieces.iter().collect() };
                    // (A closed edge's pieces go round from where the coedge before them
                    // ends: a seam may reach them at any of their vertices.)
                    if let Some(previous) = uses.last() {
                        let (s, e) = model.edge_vertices(Shape::from(previous.0));
                        let arrived = if previous.1 == Orientation::Reversed { s.node } else { e.node };
                        let start_of = |model: &Model, edge: EdgeId| {
                            let (s, e) = model.edge_vertices(Shape::from(edge));
                            if orientation == Orientation::Reversed { e.node } else { s.node }
                        };
                        if let Some(k) = ordered.iter().position(|p| start_of(model, p.0) == arrived) {
                            ordered.rotate_left(k);
                        }
                    }
                    for &&(new, lo, hi) in &ordered {
                        // Each piece in the periods where the coedge before it ends (a
                        // closed edge's pieces lie on either side of a seam).
                        let curve = match uses.last() {
                            Some(previous) if pieces.len() > 1 => {
                                let target = previous.2.point(if previous.1 == Orientation::Reversed { previous.3 } else { previous.4 });
                                let from = if orientation == Orientation::Reversed { hi } else { lo };
                                align(p, from, target, &surface)
                            }
                            _ => shifted.clone(),
                        };
                        uses.push((new, orientation, curve, lo, hi));
                    }
                    continue;
                }
            }
            if removed.contains(&edge) {
                continue;
            }
            if let Some(&(_, new, a, b)) = shortened.iter().find(|s| s.0 == edge) {
                uses.push((new, orientation, pcurve, a, b));
                continue;
            }
            uses.push((edge, orientation, pcurve, first, last));
        }
        // (Without degenerate edges at a vertex the fillet took away: no other edge of the
        // wire ends there now.)
        let ends = |model: &Model, e: EdgeId| {
            let (a, b) = model.edge_vertices(Shape::from(e));
            (a.node, b.node)
        };
        let lone: Vec<usize> = (0..uses.len())
            .filter(|&k| {
                let (a, b) = ends(model, uses[k].0);
                a == b && model.edge_curve(Shape::from(uses[k].0)).is_none() && !uses.iter().enumerate().any(|(j, u)| j != k && { let (c, d) = ends(model, u.0); c == a || d == a })
            })
            .collect();
        for k in lone.into_iter().rev() {
            uses.remove(k);
        }
        // Each between the coedge ending at one of the edge's vertices and the one starting
        // at the other; where none fits so yet, after a coedge ending at one of them (a chain
        // of them inserted one after another).
        let ends = |model: &Model, u: &Use| {
            let (s, e) = model.edge_vertices(Shape::from(u.0));
            let (s, e) = (VertexId::try_from(s.node).expect("vertex"), VertexId::try_from(e.node).expect("vertex"));
            if u.1 == Orientation::Reversed { (e, s) } else { (s, e) }
        };
        let mut pending: Vec<&(EdgeId, Curve2, [VertexId; 2])> = inserts.iter().collect();
        while !pending.is_empty() {
            let n = uses.len();
            let fits = |model: &Model, uses: &[Use], &(_, _, [a, b]): &(EdgeId, Curve2, [VertexId; 2]), loose: bool| {
                (0..n).find(|&k| {
                    let (_, end) = ends(model, &uses[k]);
                    let (start, _) = ends(model, &uses[(k + 1) % n]);
                    (end == a && (loose || start == b)) || (end == b && (loose || start == a))
                })
            };
            let found = (0..pending.len())
                .find_map(|i| fits(model, &uses, pending[i], false).map(|k| (i, k)))
                .or_else(|| (0..pending.len()).find_map(|i| fits(model, &uses, pending[i], true).map(|k| (i, k))));
            let Some((i, k)) = found else { break };
            let &(edge, ref pcurve, [a, _]) = pending.remove(i);
            let (_, end) = ends(model, &uses[k]);
            let orientation = if end == a { Orientation::Forward } else { Orientation::Reversed };
            // Starting where the previous coedge ends, in the surface's periods.
            let previous = &uses[k];
            let target = previous.2.point(if previous.1 == Orientation::Reversed { previous.3 } else { previous.4 });
            let from = if orientation == Orientation::Reversed { 1.0 } else { 0.0 };
            let shifted = align(pcurve, from, target, &surface);
            uses.insert(k + 1, (edge, orientation, shifted, 0.0, 1.0));
        }
        // (A degenerate edge kept (a pole) runs from where the coedge before it ends to where
        // the one after starts, which may have moved round the pole.)
        let n = uses.len();
        for k in 0..n {
            if n < 3 || model.edge_curve(Shape::from(uses[k].0)).is_some() {
                continue;
            }
            let (prev, next) = (&uses[(k + n - 1) % n], &uses[(k + 1) % n]);
            let a = prev.2.point(if prev.1 == Orientation::Reversed { prev.3 } else { prev.4 });
            let b = next.2.point(if next.1 == Orientation::Reversed { next.4 } else { next.3 });
            let (first, last) = (uses[k].3, uses[k].4);
            let (p, q) = if uses[k].1 == Orientation::Reversed { (b, a) } else { (a, b) };
            uses[k].2 = Curve2::BSpline(BSplineCurve2::new(vec![p, q], vec![first, last], vec![2, 2], 1));
        }
        let uses = round_pole(model, &surface, forward, uses, tolerance)?;
        wires.push((uses, wire.orientation));
    }
    new_face(model, surface, wires, tolerance)
}

/// `uses` with the way round a pole of the surface where two consecutive pcurves meet a
/// period apart in u (a cut now running round the pole, as at a fillet ending on a surface of
/// revolution at its axis): up a seam to the pole, along it, and back down the seam.
fn round_pole(model: &mut Model, surface: &Surface, face: Shape, mut uses: Vec<Use>, tolerance: f64) -> Result<Vec<Use>> {
    let Some(period) = surface.u_period() else { return Ok(uses) };
    let (_, _, v0, v1) = model.face_uv_bounds(face);
    let end_of = |u: &Use| u.2.point(if u.1 == Orientation::Reversed { u.3 } else { u.4 });
    let start_of = |u: &Use| u.2.point(if u.1 == Orientation::Reversed { u.4 } else { u.3 });
    let mut k = 0;
    while k < uses.len() {
        let (a, b) = (end_of(&uses[k]), start_of(&uses[(k + 1) % uses.len()]));
        let du = b.x - a.x;
        if !((du.abs() - period).abs() <= 1e-6 * period && (b.y - a.y).abs() <= 1e-9 * (1.0 + a.y.abs())) {
            k += 1;
            continue;
        }
        // The pole: the end of the face's v range where the surface shrinks to a point.
        let collapsed = |v: f64| surface.point(0.0, v).distance(surface.point(0.5 * period, v)) <= tolerance;
        let pole = [v0, v1].into_iter().filter(|&v| collapsed(v)).min_by(|x, y| (x - a.y).abs().total_cmp(&(y - a.y).abs()));
        let Some(pole) = pole else {
            return Err(Error::InvalidInput("a face's boundary goes round the surface where it has no pole".into()));
        };
        let (_, at) = {
            let (s, e) = model.edge_vertices(Shape::from(uses[k].0));
            if uses[k].1 == Orientation::Reversed { (e, s) } else { (s, e) }
        };
        let at = VertexId::try_from(at.node).expect("a vertex");
        let apex = model.add_vertex(surface.point(a.x, pole), tolerance);
        // The seam: from the meeting point up to the pole, along the surface's meridian.
        let n = 16;
        let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let points: Vec<Point> = params.iter().map(|&w| surface.point(a.x, a.y + (pole - a.y) * w)).collect();
        let seam = model.add_edge(at, apex, tolerance)?;
        model.set_edge_curve(seam.into(), Curve::BSpline(fit::interpolate_with_parameters(&points, &params)?), tolerance)?;
        model.set_edge_range(seam.into(), 0.0, 1.0)?;
        let segment = |p: (f64, f64), q: (f64, f64)| Curve2::BSpline(BSplineCurve2::new(vec![Point2::new(p.0, p.1), Point2::new(q.0, q.1)], vec![0.0, 1.0], vec![2, 2], 1));
        let around = model.add_edge(apex, apex, tolerance)?;
        let way = [
            (seam, Orientation::Forward, segment((a.x, a.y), (a.x, pole)), 0.0, 1.0),
            (around, Orientation::Forward, segment((a.x, pole), (b.x, pole)), 0.0, 1.0),
            (seam, Orientation::Reversed, segment((b.x, a.y), (b.x, pole)), 0.0, 1.0),
        ];
        for (j, u) in way.into_iter().enumerate() {
            uses.insert(k + 1 + j, u);
        }
        k += 4;
    }
    Ok(uses)
}

fn new_face(model: &mut Model, surface: Surface, wires: Vec<(Vec<Use>, Orientation)>, tolerance: f64) -> Result<Shape> {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, tolerance)?;
    for (uses, orientation) in wires {
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, o, pcurve, first, last) in uses {
            let c = model.add_coedge(edge, o)?;
            model.set_pcurve_with_range(c, pcurve, first, last)?;
            coedges.push(c);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, Shape { orientation, ..Shape::from(wire) })?;
    }
    Ok(face)
}

/// Parameters along a curve made continuous across the surface's periods.
fn unwrapped(mut uv: Vec<Point2>, surface: &Surface) -> Vec<Point2> {
    for i in 1..uv.len() {
        if let Some(period) = surface.u_period() {
            uv[i].x += ((uv[i - 1].x - uv[i].x) / period).round() * period;
        }
        if let Some(period) = surface.v_period() {
            uv[i].y += ((uv[i - 1].y - uv[i].y) / period).round() * period;
        }
    }
    uv
}

/// `pcurve` moved by whole periods of the surface's parameters so that its point at `t`
/// is nearest to `target`.
fn align(pcurve: &Curve2, t: f64, target: Point2, surface: &Surface) -> Curve2 {
    use crate::math::{Transform2, Transformable2, Vector2};
    let p = pcurve.point(t);
    let shift = |period: Option<f64>, d: f64| period.map_or(0.0, |period| (d / period).round() * period);
    let (du, dv) = (shift(surface.u_period(), target.x - p.x), shift(surface.v_period(), target.y - p.y));
    if du == 0.0 && dv == 0.0 {
        return pcurve.clone();
    }
    pcurve.transformed(&Transform2::translation(Vector2::new(du, dv)))
}

/// The arc of a section, from the first contact to the second, over [0, 1] (a rational
/// quadratic B-spline).
fn arc_curve(s: &Section) -> Curve {
    if s.straight {
        return Curve::BSpline(BSplineCurve::new(vec![s.contacts[0], s.contacts[1]], vec![0.0, 1.0], vec![2, 2], 1));
    }
    let (a, b) = (s.contacts[0] - s.center, s.contacts[1] - s.center);
    let w = (0.5 * a.cross(b).length().atan2(a.dot(b))).cos();
    let corner = s.center + (a + b).normalized() * (s.radius / w);
    let curve = BSplineCurve::new(vec![s.contacts[0], corner, s.contacts[1]], vec![0.0, 1.0], vec![3, 3], 2).with_weights(vec![1.0, w, 1.0]);
    Curve::BSpline(curve)
}
