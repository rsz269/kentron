//! OCCT's guided pipe trihedra, ported: `GeomFill_GuideTrihedronAC` (the guide's point at the
//! same reduced curvilinear abscissa as the spine's, a pipe shell's `CurvilinearEquivalence`)
//! and `GeomFill_GuideTrihedronPlan` (where the plane square to the spine meets the guide),
//! each under `GeomFill_LocationGuide` without contact (no rotation: its frame is the
//! trihedron's, set into continuity by the location law).

use crate::error::{Error, Result};
use crate::geometry::adaptor::{in_period, CurveAdaptor, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::curvlin::CurvlinFunc;
use crate::geometry::function_root::function_root;
use crate::geometry::int_curve_plane::{intersect_curve_plane, plane_through};
use crate::geometry::Continuity;
use crate::geometry::abscissa;
use crate::math::{Direction, Point, Vector};
use crate::topology::comp_curve::{CompCurve, EdgeWire, ShapeKeys};
use crate::topology::{Model, Shape};

use super::corrected::edge_adaptor;
use super::frenet::{fuse_intervals, Frenet, P_CONFUSION};
use super::guide_wire::{compatible_guide, compute_acr, search_origin};
use super::pipe::SpinePiece;
use super::plan_func::PlanFunc;

/// `gp_Vec::Normalized`.
fn normalized(v: Vector) -> Vector {
    let d = v.length();
    Vector::new(v.x / d, v.y / d, v.z / d)
}

/// A reduced abscissa on the guide brought into `[0, 1]`. Not OCCT's: its
/// `Approx_CurvlinFunc::GetUParameter` throws for one beyond, which a spine's abscissa can be
/// by rounding alone (at its end, an interval's bound found by summing ten steps that lands
/// past the curve's last parameter by its last bit); kentron takes the end.
fn on_guide(s: f64) -> f64 {
    s.clamp(0.0, 1.0)
}

/// The interval of `intervals` holding `param` (the first of two at a bound).
fn interval_holding(intervals: &[f64], param: f64) -> usize {
    intervals.windows(2).position(|w| w[0] <= param && param <= w[1]).unwrap_or(if param < intervals[0] { 0 } else { intervals.len() - 2 })
}

/// `GeomFill_GuideTrihedronAC` under `GeomFill_LocationGuide`, along one spine edge (as
/// `BRepFill_ACRLaw` sets it: the edge's share `[orig1, orig2]` of the spine's reduced
/// abscissa), with `SetInterval` done for each of the law's intervals.
#[doc(alias = "GeomFill_GuideTrihedronAC")]
#[doc(alias = "GeomFill_LocationGuide")]
pub(super) struct CurvilinearGuideLaw {
    /// `myCurve` and `myCurveAC`.
    curve: GeomAdaptor,
    curve_ac: CurvlinFunc<'static>,
    /// `myGuideAC` (holding `myGuide`).
    guide_ac: CurvlinFunc<'static>,
    /// `Orig1`, `Orig2`.
    orig1: f64,
    orig2: f64,
    /// The location law's intervals (`GeomAbs_C2`), and for each the trimmed spine and
    /// guide `SetInterval` makes (`myTrimmed`, `myTrimG`).
    intervals: Vec<f64>,
    trimmed: Vec<(GeomAdaptor, Box<dyn CurveAdaptor>)>,
}

/// `Approx_CurvlinFunc`'s tolerance here.
const AC_TOLERANCE: f64 = 1e-7;

impl CurvilinearGuideLaw {
    /// The law along `curve` from `first` to `last` (`last` the smaller for a reversed edge),
    /// the edge's share of the spine's reduced abscissa `[orig1, orig2]`, the guide made by
    /// `guide` (`BRepAdaptor_CompCurve` of the reoriented guide wire; each copy of the law
    /// has its own).
    pub(super) fn new(curve: &Curve, first: f64, last: f64, orig1: f64, orig2: f64, guide: Box<dyn CurveAdaptor>) -> Result<CurvilinearGuideLaw> {
        let curve = edge_adaptor(curve, first, last);
        let curve_ac = CurvlinFunc::owning(Box::new(curve.clone()), AC_TOLERANCE);
        let guide_ac = CurvlinFunc::owning(guide, AC_TOLERANCE);
        let mut law = CurvilinearGuideLaw { curve, curve_ac, guide_ac, orig1, orig2, intervals: Vec::new(), trimmed: Vec::new() };
        law.intervals = law.location_intervals(Continuity::C2)?;
        for k in 0..law.intervals.len() - 1 {
            let pair = law.set_interval(law.intervals[k], law.intervals[k + 1])?;
            law.trimmed.push(pair);
        }
        Ok(law)
    }

    /// `GeomFill_GuideTrihedronAC::Intervals`: the spine's and the guide's intervals in
    /// reduced abscissa fused (the guide's taken as they are, whatever the edge's share), as
    /// parameters of the spine.
    fn trihedron_intervals(&self, continuity: Continuity) -> Result<Vec<f64>> {
        let disc_c = self.curve_ac.intervals(continuity);
        let disc_g = self.guide_ac.intervals(continuity);
        let seq = fuse_intervals(&disc_c, &disc_g, P_CONFUSION, false);
        seq.iter().map(|&s| self.curve_ac.get_u_parameter_on(&self.curve, s)).collect()
    }

    /// `GeomFill_LocationGuide::Intervals`: the curve's and the trihedron's fused (either
    /// alone where the other has one interval).
    fn location_intervals(&self, continuity: Continuity) -> Result<Vec<f64>> {
        let of_curve = self.curve.intervals(continuity);
        let of_law = self.trihedron_intervals(continuity)?;
        if of_curve.len() == 2 {
            return Ok(of_law);
        }
        if of_law.len() == 2 {
            return Ok(of_curve);
        }
        Ok(fuse_intervals(&of_curve, &of_law, P_CONFUSION * 0.99, false))
    }

    /// `SetInterval(First, Last)`: the spine trimmed, and the guide over the stretch matched
    /// to it.
    fn set_interval(&self, first: f64, last: f64) -> Result<(GeomAdaptor, Box<dyn CurveAdaptor>)> {
        let trimmed = self.curve.trim(first, last);
        let (sf, sl) = (self.curve_ac.get_s_parameter(first), self.curve_ac.get_s_parameter(last));
        let guide = self.guide_ac.curve();
        let u = on_guide(self.orig1 + sf * (self.orig2 - self.orig1));
        let sf = self.guide_ac.get_u_parameter_on(guide, u)?;
        let u = on_guide(self.orig1 + sl * (self.orig2 - self.orig1));
        let sl = self.guide_ac.get_u_parameter_on(guide, u)?;
        Ok((trimmed, guide.trimmed(sf, sl, P_CONFUSION)))
    }

    /// `D0` in the `k`th interval: tangent and normal (`LocationGuide`'s third and first
    /// columns).
    fn d0(&self, k: usize, param: f64) -> Result<(Vector, Vector)> {
        let s = self.curve_ac.get_s_parameter(param);
        let orig_g = on_guide(self.orig1 + s * (self.orig2 - self.orig1));
        let t_g = self.guide_ac.get_u_parameter_on(self.guide_ac.curve(), orig_g)?;
        let (trimmed, trim_g) = &self.trimmed[k];
        let (p, to) = trimmed.d1(param);
        let pg = trim_g.value(t_g);
        // (OCCT guards none of these against vanishing.)
        let normal = normalized(pg - p);
        let b = to.cross(normal);
        let binormal = normalized(b);
        let tangent = normalized(normal.cross(binormal));
        Ok((tangent, normal))
    }

    /// The frame's tangent and normal at fraction `fraction` of the way along.
    pub(super) fn at_fraction(&self, fraction: f64) -> Result<(Vector, Vector)> {
        let (first, last) = (self.curve.first_parameter(), self.curve.last_parameter());
        let param = first + (last - first) * fraction;
        self.d0(interval_holding(&self.intervals, param), param)
    }
}

/// `GeomFill_GuideTrihedronPlan` under `GeomFill_LocationGuide`, along one spine edge (as
/// `BRepFill_Edge3DLaw` sets it): the Frenet tangent, the normal towards where the plane
/// square to it meets the guide (a root of `GeomFill_PlanFunc`, from a guess interpolated
/// between the meetings found at twenty points along the edge).
#[doc(alias = "GeomFill_GuideTrihedronPlan")]
pub(super) struct PlaneGuideLaw {
    /// `myCurve` and its Frenet law (never trimmed).
    curve: GeomAdaptor,
    frenet: Frenet,
    /// `myGuide` and `myTrimG` (the guide a little longer, unless periodic).
    guide: Box<dyn CurveAdaptor>,
    trim_g: Box<dyn CurveAdaptor>,
    /// `Pole`: the spine's parameter and the guide's at each of the twenty points.
    poles: Vec<(f64, f64)>,
    /// `XTol(1)`, `Inf(1)`, `Sup(1)`.
    x_tol: f64,
    inf: f64,
    sup: f64,
    /// `myStatus` after `Init`: whether a plane missed the guide.
    missed: bool,
}

/// `myNbPts`.
const NB_PTS: usize = 20;

/// `InGoodPeriod`: `current` shifted by whole periods to within half a period of `prec`.
fn in_good_period(prec: f64, period: f64, current: f64) -> f64 {
    let mut current = current;
    let diff = current - prec;
    let nb = (diff / period).trunc();
    current -= nb * period;
    let diff = current - prec;
    if diff > period / 2.0 {
        current -= period;
    } else if diff < -period / 2.0 {
        current += period;
    }
    current
}

impl PlaneGuideLaw {
    /// The law along `curve` from `first` to `last` (`last` the smaller for a reversed edge),
    /// the guide `guide` (`BRepAdaptor_CompCurve` of the guide wire): `SetCurve` and `Init`
    /// (a plane missing the guide noted, as OCCT's status notes it).
    pub(super) fn new(curve: &Curve, first: f64, last: f64, guide: Box<dyn CurveAdaptor>) -> Result<PlaneGuideLaw> {
        let curve = edge_adaptor(curve, first, last);
        let frenet = Frenet::new(curve.clone());
        let (g0, g1) = (guide.first_parameter(), guide.last_parameter());
        let x_tol = guide.resolution(1e-6);
        let mut delta_g = (g1 - g0) / 2.0;
        let (inf, sup) = (g0 - delta_g, g1 + delta_g);
        let trim_g = if !guide.is_periodic() { guide.trimmed(g0 - delta_g / 100.0, g1 + delta_g / 100.0, delta_g * 1e-7) } else { guide.trimmed(g0, g1, 0.0) };
        delta_g /= 3.0;
        let (f, l) = (curve.first_parameter(), curve.last_parameter());
        let mut poles: Vec<(f64, f64)> = Vec::with_capacity(NB_PTS);
        let mut missed = false;
        let mut w = 0.0;
        for ii in 1..=NB_PTS {
            let t = ((NB_PTS - ii) as f64 * f + (ii - 1) as f64 * l) / (NB_PTS - 1) as f64;
            let p = curve.value(t);
            let (_, (tangent, _, _)) = frenet.d0(frenet.curve(), t);
            let plane = plane_through(p, Direction::from_vector(tangent));
            let found = intersect_curve_plane(trim_g.as_ref(), &plane).map_err(|_| Error::Failed("a guide's intersection with a plane failed (Standard_NumericError)".into()))?;
            if found.points.is_empty() {
                w = if (g1 - w).abs() > (g0 - w).abs() { g0 } else { g1 };
                missed = true;
            } else {
                // (The meeting nearest the spine's point, the first of equals.)
                let mut best = &found.points[0];
                let mut d_min = p.distance(best.point);
                for point in &found.points[1..] {
                    if p.distance(point.point) < d_min {
                        best = point;
                        d_min = p.distance(point.point);
                    }
                }
                w = best.w;
            }
            if let Some(&(_, previous)) = poles.last() {
                if (w - previous).abs() > delta_g && guide.is_periodic() {
                    w = in_good_period(previous, guide.period(), w);
                }
            }
            poles.push((t, w));
        }
        Ok(PlaneGuideLaw { curve, frenet, guide, trim_g, poles, x_tol, inf, sup, missed })
    }

    /// `InitX(Param)`: the guess, the guide's parameter interpolated between the poles about
    /// `param`.
    fn init_x(&self, param: f64) -> f64 {
        let poles = &self.poles;
        let (mut ideb, mut ifin) = (0, poles.len() - 1);
        if param == poles[ideb].0 {
            ifin = ideb + 1;
        }
        if param == poles[ifin].0 {
            ideb = ifin - 1;
        }
        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            let valeur = poles[idemi].0;
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                ideb = idemi;
                ifin = ideb + 1;
            }
        }
        let (t1, t2) = (poles[ideb].0, poles[ifin].0);
        let diff = t2 - t1;
        let mut x = if diff > 1e-7 {
            let (b, a) = ((param - t1) / diff, (t2 - param) / diff);
            poles[ideb].1 * a + poles[ifin].1 * b
        } else {
            (poles[ideb].1 + poles[ifin].1) / 2.0
        };
        if self.guide.is_periodic() {
            x = in_period(x, self.guide.first_parameter(), self.guide.last_parameter());
        }
        x
    }

    /// `D0`: tangent and normal; an error where the root is not found
    /// (`GeomFill_PlaneNotIntersectGuide`).
    fn d0(&self, param: f64) -> Result<(Vector, Vector)> {
        let p = self.curve.value(param);
        let (_, (tangent, _, _)) = self.frenet.d0(self.frenet.curve(), param);
        let x = self.init_x(param);
        let mut function = PlanFunc::new(p, tangent, self.guide.as_ref());
        let root = function_root(&mut function, x, self.x_tol, self.inf, self.sup, 50);
        if !root.done {
            return Err(Error::Failed("a plane square to the spine misses the guide (GeomFill_PlaneNotIntersectGuide)".into()));
        }
        let pprime = self.trim_g.value(root.root);
        Ok((tangent, normalized(pprime - p)))
    }

    /// The frame's tangent and normal at fraction `fraction` of the way along.
    pub(super) fn at_fraction(&self, fraction: f64) -> Result<(Vector, Vector)> {
        let (first, last) = (self.curve.first_parameter(), self.curve.last_parameter());
        self.d0(first + (last - first) * fraction)
    }
}

/// A guided pipe's trihedron laws, one per spine piece, as `BRepFill_PipeShell::Set` makes
/// them for a guide wire without contact: the guide reoriented to run with the spine
/// (`BRepFill_CompatibleWires` for an open spine and guide, `BRepFill::SearchOrigin` for a
/// closed guide), taken as one curve (`BRepAdaptor_CompCurve`), and the laws by curvilinear
/// abscissa (`BRepFill_ACRLaw`, each edge its share of the spine's length) or by planes
/// (`BRepFill_Edge3DLaw`).
pub(super) enum GuideFraming {
    Curvilinear(Vec<CurvilinearGuideLaw>),
    Plane(Vec<PlaneGuideLaw>),
}

/// `BRepFill_PipeShell`'s `myTol3d` (its default).
const TOL_3D: f64 = 1e-4;

impl GuideFraming {
    pub(super) fn new(model: &Model, spine: Shape, pieces: &[SpinePiece], guide: Shape, curvilinear: bool) -> Result<GuideFraming> {
        let length = |c: &GeomAdaptor| abscissa::length(c, None);
        let mut keys = ShapeKeys::new();
        let spine_wire = EdgeWire::from_model(model, spine, &mut keys)?;
        let guide_wire = EdgeWire::from_model(model, guide, &mut keys)?;
        let the_guide = if !spine_wire.closed && !guide_wire.closed {
            compatible_guide(model, spine, guide, &length)?
        } else if guide_wire.closed {
            // (The guide's origin searched for about the middle of the spine's ends, along the
            // way back from its end to its start; or at the start of a closed spine, along it.)
            let (origin, dir) = if !spine_wire.closed {
                let (vf, vl) = spine_wire.ends();
                let missing = || Error::InvalidInput("a pipe's spine without ends".into());
                let (sp_or, p) = (vf.ok_or_else(missing)?.point, vl.ok_or_else(missing)?.point);
                (Point::new(0.5 * sp_or.x + 0.5 * p.x, 0.5 * sp_or.y + 0.5 * p.y, 0.5 * sp_or.z + 0.5 * p.z), sp_or - p)
            } else {
                CompCurve::from_wire(spine_wire.clone())?.d1(0.0)
            };
            search_origin(model, guide, origin, dir, 100.0 * TOL_3D)?
        } else {
            guide_wire
        };
        let guide_curve = CompCurve::from_wire(the_guide)?;
        if curvilinear {
            let acr = compute_acr(&spine_wire, &length);
            let laws = pieces
                .iter()
                .enumerate()
                .map(|(k, p)| {
                    let t1 = if k == 0 { 0.0 } else { acr[k] };
                    CurvilinearGuideLaw::new(&p.curve, p.first, p.last, t1, acr[k + 1], Box::new(guide_curve.clone()))
                })
                .collect::<Result<_>>()?;
            Ok(GuideFraming::Curvilinear(laws))
        } else {
            let laws = pieces.iter().map(|p| PlaneGuideLaw::new(&p.curve, p.first, p.last, Box::new(guide_curve.clone()))).collect::<Result<_>>()?;
            Ok(GuideFraming::Plane(laws))
        }
    }

    /// `ErrorStatus()` after `SetCurve`: an error where a plane missed the guide
    /// (`GeomFill_PlaneNotIntersectGuide`, which fails the pipe).
    pub(super) fn status(&self) -> Result<()> {
        if let GuideFraming::Plane(laws) = self {
            if laws.iter().any(|l| l.missed) {
                return Err(Error::Failed("a plane square to the spine misses the guide (GeomFill_PlaneNotIntersectGuide)".into()));
            }
        }
        Ok(())
    }

    /// The frame's tangent and normal at fraction `s` of piece `k`.
    pub(super) fn at_fraction(&self, k: usize, s: f64) -> Result<(Vector, Vector)> {
        match self {
            GuideFraming::Curvilinear(laws) => laws[k].at_fraction(s),
            GuideFraming::Plane(laws) => laws[k].at_fraction(s),
        }
    }
}
