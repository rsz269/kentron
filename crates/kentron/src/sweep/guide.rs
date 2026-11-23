//! Frames steered by a guide wire: along each spine piece, the normal points from the spine
//! point to a matched point of the guide, and the tangent is the spine's, made square to it
//! (curvilinear matching) or the spine's own (plane matching).
//!
//! - **Curvilinear**: the guide point at the same reduced arc length: each spine edge's share
//!   of the spine's length is mapped onto the same share of the guide's.
//! - **Plane**: where the plane through the spine point, square to the spine's tangent, meets
//!   the guide; the meeting followed is the one a table of meetings at twenty points along
//!   the piece leads a root search to.

use crate::error::{Error, Result};
use crate::geometry::abscissa;
use crate::geometry::adaptor::{CurveAdaptor, GeomAdaptor};
use crate::geometry::curvlin::CurvlinFunc;
use crate::geometry::function_root::function_root;
use crate::geometry::int_curve_plane::{intersect_curve_plane, plane_through};
use crate::geometry::Continuity;
use crate::math::{Direction, Vector};
use crate::topology::comp_curve::{CompCurve, EdgeWire, ShapeKeys};
use crate::topology::{Model, Orientation, Shape};

use super::frenet::{travelled, unit, Frenet};
use super::guide_wire::{compatible_guide, compute_acr, search_origin};
use super::pipe::SpinePiece;
use super::plan_func::PlanFunc;

/// Arc lengths are computed to this.
const LENGTH_TOLERANCE: f64 = 1e-7;

/// Meetings of plane and guide tabulated along each piece.
const TABLE_POINTS: usize = 20;

/// One spine piece of a curvilinear guide law.
pub(super) struct CurvilinearPiece {
    spine: CurvlinFunc<'static>,
    /// The piece's C1 stretches (for the tangent at a break, from the stretch before).
    breaks: Vec<f64>,
    view: GeomAdaptor,
    /// The share of the guide's length the piece is matched with.
    share: (f64, f64),
}

/// One spine piece of a plane guide law.
pub(super) struct PlanePiece {
    frenet: Frenet,
    /// Table parameters along the piece and the guide parameters met there.
    table: Vec<(f64, f64)>,
    missed: bool,
}

/// The frames of every spine piece towards a guide.
pub(super) enum GuideFraming {
    Curvilinear { guide: CurvlinFunc<'static>, pieces: Vec<CurvilinearPiece> },
    Plane { guide: CompCurve, pieces: Vec<PlanePiece> },
}

/// The guide wire prepared to run with the spine (see `guide_wire`), as one curve.
fn prepared_guide(model: &Model, spine: Shape, spine_wire: &EdgeWire, guide: Shape) -> Result<CompCurve> {
    let guide_wire = EdgeWire::from_model(model, Shape { orientation: Orientation::Forward, ..guide }, &mut ShapeKeys::new())?;
    let length = |c: &GeomAdaptor| abscissa::length(c, None);
    let prepared = match (spine_wire.closed, guide_wire.closed) {
        (false, false) => compatible_guide(model, spine, guide, &length)?,
        (false, true) => {
            let (Some(a), Some(b)) = spine_wire.ends() else { return Err(Error::Failed("a spine without ends".into())) };
            let origin = a.point + (b.point - a.point) * 0.5;
            search_origin(model, guide, origin, a.point - b.point, 1e-2)?
        }
        (true, true) => {
            let curve = CompCurve::from_wire(spine_wire.clone())?;
            let (p, d) = curve.d1(0.0);
            search_origin(model, guide, p, d, 1e-2)?
        }
        (true, false) => guide_wire,
    };
    CompCurve::from_wire(prepared)
}

impl GuideFraming {
    /// The laws along the spine `spine` (its edges `pieces`, in connection order) towards
    /// `guide`, by curvilinear abscissa or by planes.
    pub(super) fn new(model: &Model, spine: Shape, pieces: &[SpinePiece], guide: Shape, curvilinear: bool) -> Result<GuideFraming> {
        let spine_wire = EdgeWire::from_model(model, Shape { orientation: Orientation::Forward, ..spine }, &mut ShapeKeys::new())?;
        let curve = prepared_guide(model, spine, &spine_wire, guide)?;
        if curvilinear {
            let length = |c: &GeomAdaptor| abscissa::length(c, None);
            let acr = compute_acr(&spine_wire, &length);
            let guide = CurvlinFunc::owning(Box::new(curve), LENGTH_TOLERANCE);
            let pieces = pieces
                .iter()
                .enumerate()
                .map(|(k, p)| {
                    let (view, _, _) = travelled(&p.curve, p.first, p.last);
                    let breaks = view.intervals(Continuity::C1);
                    let alpha = if k == 0 { 0.0 } else { acr.get(k).copied().unwrap_or(1.0) };
                    let beta = acr.get(k + 1).copied().unwrap_or(1.0);
                    CurvilinearPiece { spine: CurvlinFunc::owning(Box::new(view.clone()), LENGTH_TOLERANCE), breaks, view, share: (alpha, beta) }
                })
                .collect();
            return Ok(GuideFraming::Curvilinear { guide, pieces });
        }
        let (g0, g1) = (curve.first_parameter(), curve.last_parameter());
        let delta = 0.5 * (g1 - g0);
        let extended = curve.trim(g0 - delta / 100.0, g1 + delta / 100.0, 0.0);
        let mut out = Vec::with_capacity(pieces.len());
        for p in pieces {
            let (view, a, b) = travelled(&p.curve, p.first, p.last);
            let frenet = Frenet::without_singular_points(view);
            let mut table = Vec::with_capacity(TABLE_POINTS);
            let mut missed = false;
            for i in 0..TABLE_POINTS {
                let t = if i + 1 == TABLE_POINTS { b } else { a + (b - a) * i as f64 / (TABLE_POINTS - 1) as f64 };
                let point = frenet.curve().value(t);
                let tangent = frenet.tangent(t);
                let plane = plane_through(point, Direction::from_vector(tangent));
                let meetings = intersect_curve_plane(&extended, &plane).map_err(|e| Error::Failed(format!("a guide's plane intersection: {e:?}")))?;
                let nearest = meetings.points.iter().fold(None::<(f64, f64)>, |best, m| {
                    let d = m.point.distance(point);
                    if best.is_none_or(|(bd, _)| d < bd) { Some((d, m.w)) } else { best }
                });
                let w = match nearest {
                    Some((_, w)) => w,
                    None => {
                        missed = true;
                        let previous = table.last().map_or(g0, |&(_, w): &(f64, f64)| w);
                        if (previous - g0).abs() <= (previous - g1).abs() { g0 } else { g1 }
                    }
                };
                table.push((t, w));
            }
            out.push(PlanePiece { frenet, table, missed });
        }
        Ok(GuideFraming::Plane { guide: curve, pieces: out })
    }

    /// `Err` where a plane square to the spine misses the guide (plane matching).
    pub(super) fn status(&self) -> Result<()> {
        match self {
            GuideFraming::Plane { pieces, .. } if pieces.iter().any(|p| p.missed) => Err(Error::Failed("a plane square to the spine misses the guide".into())),
            _ => Ok(()),
        }
    }

    /// (tangent, normal) at fraction `s` of spine piece `k`.
    pub(super) fn at_fraction(&self, k: usize, s: f64) -> Result<(Vector, Vector)> {
        match self {
            GuideFraming::Curvilinear { guide, pieces } => curvilinear_frame(guide, &pieces[k], s),
            GuideFraming::Plane { guide, pieces } => plane_frame(guide, &pieces[k], s),
        }
    }
}

fn degenerate() -> Error {
    Error::Failed("a guided frame with its guide point on the spine or along its tangent".into())
}

/// The frame towards the guide point at the piece's reduced arc length.
fn curvilinear_frame(guide: &CurvlinFunc<'static>, piece: &CurvilinearPiece, s: f64) -> Result<(Vector, Vector)> {
    let (a, b) = (piece.view.first_parameter(), piece.view.last_parameter());
    let u = if s == 1.0 { b } else { a + (b - a) * s };
    let sigma = piece.spine.get_s_parameter(u);
    let (alpha, beta) = piece.share;
    let g = (alpha + sigma * (beta - alpha)).clamp(0.0, 1.0);
    let w = guide.get_u_parameter(g)?;
    let q = guide.curve().value(w);
    // (The tangent from the C1 stretch holding u: the one before where it is a bound.)
    let breaks = &piece.breaks;
    let i = breaks[1..].partition_point(|&x| x < u).min(breaks.len() - 2);
    let view = if breaks.len() == 2 { piece.view.clone() } else { piece.view.trim(breaks[i], breaks[i + 1]) };
    let (p, d1) = view.d1(u);
    let towards = q - p;
    let (nl, bl) = (towards.length(), d1.cross(towards).length());
    if !(nl > 0.0) || !(bl > 0.0) {
        return Err(degenerate());
    }
    let n = towards * (1.0 / nl);
    let bn = unit(d1.cross(n));
    let t = unit(n.cross(bn));
    if !(t.length().is_finite() && bn.length().is_finite()) {
        return Err(degenerate());
    }
    Ok((t, n))
}

/// The frame towards the guide point on the plane square to the spine.
fn plane_frame(guide: &CompCurve, piece: &PlanePiece, s: f64) -> Result<(Vector, Vector)> {
    let frenet = &piece.frenet;
    let u = frenet.at_fraction(s);
    let t = frenet.tangent(u);
    let p = frenet.curve().value(u);
    let table = &piece.table;
    let j = table.partition_point(|e| e.0 < u).clamp(1, table.len() - 1);
    let ((t0, w0), (t1, w1)) = (table[j - 1], table[j]);
    let guess = if (t1 - t0).abs() < 1e-7 { 0.5 * (w0 + w1) } else { w0 + (w1 - w0) * (u - t0) / (t1 - t0) };
    let (g0, g1) = (guide.first_parameter(), guide.last_parameter());
    let delta = 0.5 * (g1 - g0);
    let mut function = PlanFunc::new(p, t, guide);
    let root = function_root(&mut function, guess, guide.resolution(1e-6), g0 - delta, g1 + delta, 50);
    if !root.done {
        return Err(Error::Failed("no point of the guide on the plane square to the spine".into()));
    }
    let towards = guide.value(root.root) - p;
    if !(towards.length() > 0.0) {
        return Err(degenerate());
    }
    Ok((t, unit(towards)))
}

#[cfg(test)]
mod tests {
    use crate::geometry::adaptor::CurveAdaptor;
    use crate::geometry::curve::Curve;
    use crate::geometry::ParametricCurve;
    use crate::make;
    use crate::math::{Point, Vector};
    use crate::topology::Model;

    /// A straight spine up the axis of a helical guide (a fitted B-spline, one and a half
    /// turns): by curvilinear abscissa the normal points at the guide point at the same
    /// fraction of the guide's length, by planes at the guide point at the same height.
    #[test]
    fn a_coaxial_helix_guides_the_normal_round_the_axis() {
        let mut model = Model::new();
        let (r, height, turns) = (2.0, 9.0, 1.5);
        let points: Vec<Point> = (0..=96)
            .map(|i| {
                let a = std::f64::consts::TAU * turns * i as f64 / 96.0;
                Point::new(r * a.cos(), r * a.sin(), height * i as f64 / 96.0)
            })
            .collect();
        let helix = crate::geometry::fit::interpolate(&points, false, 1e-9).unwrap();
        let curve = Curve::BSpline(helix);
        let guide = make::edge_from_curve(&mut model, curve.clone()).unwrap();
        let spine = make::edge_between(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, height)).unwrap();
        let fractions: Vec<f64> = (0..=12).map(|i| i as f64 / 12.0).collect();
        let by_length = super::super::guided_frames(&model, spine, guide, true, &fractions).unwrap();
        let by_planes = super::super::guided_frames(&model, spine, guide, false, &fractions).unwrap();
        let view = crate::geometry::adaptor::GeomAdaptor::new(&curve, ParametricCurve::first_parameter(&curve), ParametricCurve::last_parameter(&curve));
        let total = crate::geometry::abscissa::length(&view, None);
        for (i, &f) in fractions.iter().enumerate() {
            let p = Point::new(0.0, 0.0, height * f);
            // (The guide point at the fraction f of its length.)
            let u = if f == 0.0 {
                view.first_parameter()
            } else {
                crate::geometry::abscissa::travel(&view, view.first_parameter(), f * total, true, 1e-10).unwrap()
            };
            let q = curve.point(u);
            let (t, n) = by_length[0][i];
            let towards = (q - p).normalized();
            assert!((n - towards).length() < 1e-6, "by length at {f}: {n:?} vs {towards:?}");
            assert!(t.dot(n).abs() < 1e-12);
            // (By planes: the guide point at the spine point's height, level with it.)
            let (t, n) = by_planes[0][i];
            assert!((t - Vector::new(0.0, 0.0, 1.0)).length() < 1e-12);
            assert!(n.z.abs() < 1e-6, "by planes at {f}: {n:?}");
            let a = std::f64::consts::TAU * turns * f;
            assert!((n - Vector::new(a.cos(), a.sin(), 0.0)).length() < 1e-4, "by planes at {f}: {n:?}");
        }
    }
}
