//! Points onto surfaces and curves, the places where a surface collapses, and whether it
//! closes on itself: what building a face's parameter-space curve needs to know of the face's
//! surface (see [`pcurve_projection`](super::pcurve_projection)).
//!
//! - Point inversion: the parameters of the surface point nearest a 3D point. Elementary
//!   surfaces (plane, cylinder, cone, sphere, torus) answer in closed form; others by a
//!   search of the whole patch (a grid of samples and Newton steps on the squared distance,
//!   [`extrema::nearest_on_surface`]), or, continuing from the parameters of a nearby point,
//!   by Newton steps from there, checked to have reached a true foot.
//! - Singular places: boundary isos that the surface shrinks to (nearly) one point, such as
//!   a sphere's poles or a cone's apex, with an upper estimate of their 3D size.
//! - Closure: whether the surface meets itself across its bounds in `u` or in `v`.
//! - Feet on curves: the nearest point of a curve over a range, or the distance extremum
//!   nearest a previous parameter (to follow successive points along an iso).

use crate::math::{Point, Point2, Vector};
use crate::tolerance::LINEAR;

use super::adaptor::GeomAdaptor;
use super::curve::Curve;
use super::extrema;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// The value to add to `value` (a whole multiple of `period`) so that it lies within half a
/// period of `target`: 0 when `|value - target| <= period/2`. For a vanishing period, the
/// whole difference `target - value`.
pub(crate) fn period_shift(value: f64, target: f64, period: f64) -> f64 {
    let period = period.abs();
    if period <= 1e-100 {
        return target - value;
    }
    let d = target - value;
    if !(d.abs() > 0.5 * period) {
        return 0.0;
    }
    let k = (d / period).round();
    k * period
}

/// A foot of a point on a curve.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CurveFoot {
    pub parameter: f64,
    /// The curve's point at `parameter`.
    pub point: Point,
    /// The distance from the point projected.
    pub distance: f64,
}

impl CurveFoot {
    fn at(curve: &Curve, p: Point, parameter: f64) -> CurveFoot {
        let point = curve.point(parameter);
        CurveFoot { parameter, point, distance: point.distance(p) }
    }
}

/// Whether a curve has ends of its own (rather than being a line or conic).
fn bounded(curve: &Curve) -> bool {
    matches!(curve, Curve::BSpline(_) | Curve::Bezier(_) | Curve::Trimmed(_))
}

/// The range `[first, last]` in order, widened at both ends for a curve that does not close
/// (feet slightly beyond its ends are allowed), within the curve's own domain.
fn search_range(curve: &Curve, preci: f64, first: f64, last: f64) -> (f64, f64) {
    let (a, b) = (first.min(last), first.max(last));
    let closed = curve.is_periodic() || (bounded(curve) && curve.is_closed());
    if closed || !(b - a).is_finite() {
        return (a, b);
    }
    let widen = super::resolution::curve_resolution(curve, preci.max(LINEAR)).min(0.1 * (b - a));
    let widen = if widen.is_finite() && widen > 0.0 { widen } else { 0.0 };
    let (lo, hi) = (a - widen, b + widen);
    if bounded(curve) {
        (lo.max(curve.first_parameter().min(a)), hi.min(curve.last_parameter().max(b)))
    } else {
        (lo, hi)
    }
}

/// An end of a bounded curve within `tol` of `p`, the first end first.
fn snapped_end(curve: &Curve, p: Point, a: f64, b: f64, tol: f64) -> Option<CurveFoot> {
    if !bounded(curve) {
        return None;
    }
    [a, b].into_iter().map(|t| CurveFoot::at(curve, p, t)).find(|f| f.distance <= tol)
}

/// The nearest point of `curve` over `[first, last]` (order of the bounds free). A bounded
/// curve's end within `preci` of `p` (`snap_to_ends`, else within `tolerance::LINEAR`) is
/// taken as it is; otherwise the range is widened a little beyond the ends of a curve that
/// does not close, and the global minimum of the distance over it found (in closed form on
/// lines and circles).
pub(crate) fn foot_on_curve(curve: &Curve, p: Point, preci: f64, first: f64, last: f64, snap_to_ends: bool) -> CurveFoot {
    let (a, b) = (first.min(last), first.max(last));
    let tol = if snap_to_ends { preci.max(0.0) } else { LINEAR };
    if let Some(end) = snapped_end(curve, p, a, b, tol) {
        return end;
    }
    let (lo, hi) = search_range(curve, preci, a, b);
    nearest_over(curve, p, lo, hi)
}

/// The distance extremum of `curve` nearest the parameter `previous` (a continuation of a
/// previous foot) over `[first, last]` widened as in [`foot_on_curve`] (a bounded curve's end
/// within `tolerance::LINEAR` of `p` taken as it is); where the range holds no extremum, the
/// nearest point.
pub(crate) fn foot_on_curve_from(previous: f64, curve: &Curve, p: Point, preci: f64, first: f64, last: f64) -> CurveFoot {
    let (a, b) = (first.min(last), first.max(last));
    if let Some(end) = snapped_end(curve, p, a, b, LINEAR) {
        return end;
    }
    let (lo, hi) = search_range(curve, preci, a, b);
    if lo.is_finite() && hi.is_finite() && previous.is_finite() {
        let view = GeomAdaptor::new(curve, lo, hi);
        if let Some(found) = super::locate_ext_pc::locate_ext_pc(p, &view, previous, lo, hi, 1e-10) {
            return CurveFoot::at(curve, p, found.parameter);
        }
    }
    nearest_over(curve, p, lo, hi)
}

/// The global nearest point of `curve` over `[lo, hi]` (ends included).
fn nearest_over(curve: &Curve, p: Point, lo: f64, hi: f64) -> CurveFoot {
    let inf = |t: f64| t.abs() >= 1e100;
    match curve {
        Curve::Line(line) => {
            let d = Vector::from(line.direction);
            let t = (p - line.origin).dot(d);
            let t = if inf(lo) { t } else { t.max(lo) };
            let t = if inf(hi) { t } else { t.min(hi) };
            return CurveFoot::at(curve, p, t);
        }
        Curve::Circle(circle) => {
            let f = circle.frame();
            let q = p - f.origin();
            let (x, y) = (q.dot(Vector::from(f.x())), q.dot(Vector::from(f.y())));
            let angle = if x == 0.0 && y == 0.0 { lo } else { y.atan2(x) };
            let middle = 0.5 * (lo + hi);
            let t = angle + period_shift(angle, middle, std::f64::consts::TAU);
            if t >= lo && t <= hi {
                return CurveFoot::at(curve, p, t);
            }
            let (s, e) = (CurveFoot::at(curve, p, lo), CurveFoot::at(curve, p, hi));
            return if e.distance < s.distance { e } else { s };
        }
        _ => {}
    }
    if !(lo.is_finite() && hi.is_finite()) || inf(lo) || inf(hi) {
        let (a, b) = (if inf(lo) { -1e3 } else { lo }, if inf(hi) { 1e3 } else { hi });
        let found = extrema::nearest_on_curve(p, curve, a, b);
        return CurveFoot { parameter: found.parameter, point: found.point, distance: found.distance };
    }
    let mut best = {
        let (s, e) = (CurveFoot::at(curve, p, lo), CurveFoot::at(curve, p, hi));
        if e.distance < s.distance { e } else { s }
    };
    let view = GeomAdaptor::new(curve, lo, hi);
    match super::ext_pc::ext_pc_other(&view, p, lo, hi, 1e-10) {
        Some(found) => {
            for e in found.iter().filter(|e| e.is_min) {
                let foot = CurveFoot::at(curve, p, e.parameter.clamp(lo, hi));
                if foot.distance < best.distance {
                    best = foot;
                }
            }
        }
        None => {
            let found = extrema::nearest_on_curve(p, curve, lo, hi);
            if found.distance < best.distance {
                best = CurveFoot { parameter: found.parameter, point: found.point, distance: found.distance };
            }
        }
    }
    best
}

/// A place where the surface collapses one of its boundary isos to (nearly) a point.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SingularPlace {
    /// A representative 3D point of the place.
    pub point: Point,
    /// Its 3D size: 0 for a true pole, an upper estimate otherwise.
    pub spread: f64,
    /// True: the place is the iso `u = fixed` (`v` runs along it); false: the iso `v = fixed`.
    pub fixed_is_u: bool,
    /// The constant parameter of that iso.
    pub fixed: f64,
    /// The range of the other parameter along it.
    pub free_range: (f64, f64),
}

/// What a surface is, for inversion.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Kind {
    Plane,
    /// A cylinder, cone or sphere: `u` an angle.
    Quadric,
    Torus,
    /// Anything inverted numerically.
    FreeForm,
}

/// The measurement behind a closure verdict: the largest gap between the two boundary isos,
/// and the distance from the boundary iso to the middle one where it is attained (`None`:
/// decided without measuring).
#[derive(Clone, Copy, Debug)]
enum Closure {
    Decided(bool),
    Measured { gap: f64, middle: f64 },
}

/// A surface prepared for repeated point inversion (caches what it learns about the surface;
/// nothing about the points).
pub(crate) struct SurfaceInverse<'a> {
    surface: &'a Surface,
    /// The elementary surface a trimmed one is part of, else the surface itself.
    basis: &'a Surface,
    kind: Kind,
    bounds: (f64, f64, f64, f64),
    places: Option<Vec<SingularPlace>>,
    closure_u: Option<Closure>,
    closure_v: Option<Closure>,
    /// The knots of reduced (C0) continuity of a B-spline surface, in `u` and in `v`.
    c0_knots: Option<(Vec<f64>, Vec<f64>)>,
}

/// The surface under any rectangular trims.
fn untrimmed(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => untrimmed(t.basis()),
        s => s,
    }
}

fn infinite(x: f64) -> bool {
    !(x.abs() < 1e100)
}

impl<'a> SurfaceInverse<'a> {
    pub(crate) fn new(surface: &'a Surface) -> Self {
        // (An offset of an elementary surface shares its normal lines, so its feet: inverted
        // by the basis's closed form, but a torus turned through its tube's middle.)
        let basis = match untrimmed(surface) {
            Surface::Offset(o) if !matches!(o.basis(), Surface::Torus(_)) => untrimmed(o.basis()),
            s => s,
        };
        let kind = match basis {
            Surface::Plane(_) => Kind::Plane,
            Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) => Kind::Quadric,
            Surface::Torus(_) => Kind::Torus,
            _ => Kind::FreeForm,
        };
        let basis = if kind == Kind::FreeForm { surface } else { basis };
        SurfaceInverse { surface, basis, kind, bounds: surface.bounds(), places: None, closure_u: None, closure_v: None, c0_knots: None }
    }

    /// Whether the surface is inverted in closed form (a plane, cylinder, cone, sphere or
    /// torus, or a part of one).
    pub(crate) fn is_elementary(&self) -> bool {
        self.kind != Kind::FreeForm
    }

    pub(crate) fn surface(&self) -> &'a Surface {
        self.surface
    }

    /// `surface.bounds()`.
    pub(crate) fn bounds(&self) -> (f64, f64, f64, f64) {
        self.bounds
    }

    /// The parameters `(u, v)` and the gap of `p` in closed form, angles brought within half
    /// a turn of `(ut, vt)` (`None` for free-form surfaces).
    fn closed_form(&self, p: Point, ut: f64, vt: f64) -> Option<(Point2, f64)> {
        if self.kind == Kind::FreeForm {
            return None;
        }
        let (u, v) = self.basis.parameters(p)?;
        let turn = std::f64::consts::TAU;
        let (u, v) = match self.kind {
            Kind::Plane => (u, v),
            Kind::Quadric => (u + period_shift(u, ut, turn), v),
            _ => (u + period_shift(u, ut, turn), v + period_shift(v, vt, turn)),
        };
        Some((Point2::new(u, v), self.surface.point(u, v).distance(p)))
    }

    /// The middles of the bounds (0 where infinite).
    fn middles(&self) -> (f64, f64) {
        let (u0, u1, v0, v1) = self.bounds;
        let mid = |a: f64, b: f64| if infinite(a) || infinite(b) { 0.0 } else { 0.5 * (a + b) };
        (mid(u0, u1), mid(v0, v1))
    }

    /// Global inversion: the parameters of the point of the surface nearest `p` (within the
    /// bounds, angles of elementary surfaces within half a turn of the bounds' middles), and
    /// the gap.
    pub(crate) fn nearest(&mut self, p: Point, preci: f64) -> (Point2, f64) {
        let _ = preci;
        let (um, vm) = self.middles();
        if let Some(found) = self.closed_form(p, um, vm) {
            return found;
        }
        self.global_free_form(p)
    }

    /// The global foot on a free-form surface.
    fn global_free_form(&self, p: Point) -> (Point2, f64) {
        let s = self.surface;
        let (u0, u1, v0, v1) = self.bounds;
        // (Surfaces of revolution, extrusions and turned tori have closed forms of their own;
        // taken where they land within the bounds.)
        if matches!(s, Surface::Revolution(_) | Surface::Extrusion(_) | Surface::Offset(_)) {
            if let Some((u, v)) = s.parameters(p) {
                let u = match s.u_period() {
                    Some(period) if !infinite(u0) && !infinite(u1) => u + period_shift(u, 0.5 * (u0 + u1), period),
                    _ => u,
                };
                let slack = 1e-9 * (1.0 + u.abs().max(v.abs()));
                let inside = |x: f64, a: f64, b: f64| x >= a - slack && x <= b + slack;
                if u.is_finite() && v.is_finite() && inside(u, u0, u1) && inside(v, v0, v1) {
                    let found = (Point2::new(u, v), s.point(u, v).distance(p));
                    // (Checked against a search of the patch: a closed form of a curve's
                    // nearest point may pick a far branch.)
                    let searched = extrema::nearest_on_surface(p, &self.searchable(), None);
                    let searched = self.clamped(searched);
                    return if searched.1 + LINEAR < found.1 { searched } else { found };
                }
            }
            let searched = extrema::nearest_on_surface(p, &self.searchable(), None);
            return self.clamped(searched);
        }
        let (u, v, d) = extrema::nearest_on_surface(p, s, None);
        self.across_the_bounds(p, (Point2::new(u, v), d))
    }

    /// A foot on a bound, compared with Newton steps from the opposite bound (where a closed
    /// surface meets itself, a search may stop on the wrong side of its seam): the better.
    fn across_the_bounds(&self, p: Point, found: (Point2, f64)) -> (Point2, f64) {
        let s = self.surface;
        let (u0, u1, v0, v1) = self.bounds;
        let (uv, mut best) = (found.0, found);
        let near = |x: f64, a: f64, b: f64| (x - a).abs() <= 1e-9 * (1.0 + (b - a).abs());
        let mut starts = Vec::new();
        if !s.is_u_periodic() && !infinite(u0) && !infinite(u1) {
            if near(uv.x, u0, u1) {
                starts.push((u1, uv.y));
            } else if near(uv.x, u1, u0) {
                starts.push((u0, uv.y));
            }
        }
        if !s.is_v_periodic() && !infinite(v0) && !infinite(v1) {
            if near(uv.y, v0, v1) {
                starts.push((uv.x, v1));
            } else if near(uv.y, v1, v0) {
                starts.push((uv.x, v0));
            }
        }
        for start in starts {
            let foot = extrema::newton_foot(p, s, start);
            let d = foot.point.distance(p);
            if d + LINEAR < best.1 && foot.u.is_finite() && foot.v.is_finite() {
                best = (Point2::new(foot.u, foot.v), d);
            }
        }
        best
    }

    /// The surface as searched without closed forms (a trim over the whole bounds for those
    /// with a closed form, so that the search runs over the patch).
    fn searchable(&self) -> Surface {
        let (u0, u1, v0, v1) = self.bounds;
        let finite = |a: f64, b: f64| if infinite(a) || infinite(b) { (-1e3, 1e3) } else { (a, b) };
        let ((u0, u1), (v0, v1)) = (finite(u0, u1), finite(v0, v1));
        Surface::Trimmed(super::surface::TrimmedSurface::new(self.surface.clone(), u0, u1, v0, v1))
    }

    fn clamped(&self, found: (f64, f64, f64)) -> (Point2, f64) {
        (Point2::new(found.0, found.1), found.2)
    }

    /// Inversion continued from `hint` (the parameters of a nearby point): on elementary
    /// surfaces the closed form with every angle within half a turn of the hint's; on others
    /// Newton steps from the hint, kept where they reach a true foot no farther than the
    /// hint's own point, else (or where the gap exceeds `max_gap > 0`, the surface is nearly
    /// degenerate there, or the hint lies on a knot line where the surface is only C0) the
    /// better of that and the global answer.
    pub(crate) fn nearest_from(&mut self, hint: Point2, p: Point, preci: f64, max_gap: f64) -> (Point2, f64) {
        if let Some(found) = self.closed_form(p, hint.x, hint.y) {
            return found;
        }
        if !(hint.x.is_finite() && hint.y.is_finite()) {
            return self.nearest(p, preci);
        }
        let local = self.local_foot(hint, p);
        let on_kink = self.on_c0_knot(hint);
        match local {
            Some((uv, gap, degenerate)) => {
                let suspicious = degenerate || on_kink || (max_gap > 0.0 && gap > max_gap + LINEAR);
                if !suspicious {
                    return (uv, gap);
                }
                let global = self.nearest(p, preci);
                if global.1 + LINEAR < gap { global } else { (uv, gap) }
            }
            None => self.nearest(p, preci),
        }
    }

    /// Newton steps from `hint` to a foot of `p`: its parameters, gap and whether the surface
    /// is nearly degenerate there; `None` where the steps do not reach a true foot (the
    /// residual not orthogonal to the surface, allowing for a bound) or end farther away than
    /// the hint's point.
    fn local_foot(&self, hint: Point2, p: Point) -> Option<(Point2, f64, bool)> {
        let s = self.surface;
        let foot = extrema::newton_foot(p, s, (hint.x, hint.y));
        let (u, v) = (foot.u, foot.v);
        if !(u.is_finite() && v.is_finite()) {
            return None;
        }
        let r = p - foot.point;
        let gap = r.length();
        let at_hint = s.point(hint.x, hint.y).distance(p);
        if !(gap <= at_hint + LINEAR) {
            return None;
        }
        let (u0, u1, v0, v1) = self.bounds;
        let (lu, lv) = (foot.du.length(), foot.dv.length());
        let tiny = 1e-300;
        // The residual's component along each direction: zero at an interior foot, pointing
        // out of the bounds at a foot on a bound.
        let settled = |g: f64, len: f64, x: f64, lo: f64, hi: f64, periodic: bool| {
            if len <= tiny {
                return true;
            }
            let along = g / len;
            if along.abs() <= LINEAR {
                return true;
            }
            !periodic && ((x <= lo && along < 0.0) || (x >= hi && along > 0.0))
        };
        let ok = settled(r.dot(foot.du), lu, u, u0, u1, s.is_u_periodic()) && settled(r.dot(foot.dv), lv, v, v0, v1, s.is_v_periodic());
        if !ok {
            return None;
        }
        let degenerate = !(foot.du.cross(foot.dv).length() > 1e-6 * lu * lv);
        Some((Point2::new(u, v), gap, degenerate))
    }

    /// Whether `uv` lies (within `tolerance::LINEAR` in parameter) on a knot line where the
    /// B-spline surface is only C0.
    fn on_c0_knot(&mut self, uv: Point2) -> bool {
        if self.c0_knots.is_none() {
            let mut knots = (Vec::new(), Vec::new());
            if let Surface::BSpline(b) = untrimmed(self.surface) {
                let inner = |ks: &[f64], ms: &[usize], degree: usize| -> Vec<f64> {
                    let n = ks.len();
                    ks.iter().zip(ms).enumerate().filter(|&(i, (_, &m))| i > 0 && i + 1 < n && m >= degree).map(|(_, (&k, _))| k).collect()
                };
                knots = (inner(b.u_knots(), b.u_multiplicities(), b.u_degree()), inner(b.v_knots(), b.v_multiplicities(), b.v_degree()));
            }
            self.c0_knots = Some(knots);
        }
        let (us, vs) = self.c0_knots.as_ref().expect("knots");
        us.iter().any(|k| (uv.x - k).abs() <= LINEAR) || vs.iter().any(|k| (uv.y - k).abs() <= LINEAR)
    }

    /// All singular places, ascending `spread` (computed once).
    pub(crate) fn singular_places(&mut self) -> &[SingularPlace] {
        if self.places.is_none() {
            let mut places = self.find_places();
            places.sort_by(|a, b| a.spread.total_cmp(&b.spread));
            self.places = Some(places);
        }
        self.places.as_deref().expect("places")
    }

    fn find_places(&self) -> Vec<SingularPlace> {
        let (u0, u1, v0, v1) = self.bounds;
        let s = self.surface;
        let along_v = |fixed: f64, spread: f64, point: Point| SingularPlace { point, spread, fixed_is_u: false, fixed, free_range: (u0, u1) };
        match s {
            Surface::Cone(cone) => {
                let v = -cone.radius() / cone.semi_angle_sin_cos().0;
                vec![along_v(v, 0.0, cone.apex())]
            }
            Surface::Sphere(_) => vec![along_v(v1, 0.0, s.point(u0, v1)), along_v(v0, 0.0, s.point(u0, v0))],
            Surface::Torus(torus) => {
                let (big, small) = (torus.major_radius(), torus.minor_radius());
                let pi = std::f64::consts::PI;
                let ratio = if small > 0.0 { (big / small).min(1.0) } else { 1.0 };
                let inner = pi - ratio.acos();
                let mut out = vec![along_v(inner, (big - small).max(0.0), s.point(u0, inner))];
                if big <= small {
                    let outer = pi + ratio.acos();
                    out.push(along_v(outer, 0.0, s.point(u0, outer)));
                }
                out
            }
            Surface::Plane(_) | Surface::Cylinder(_) | Surface::Extrusion(_) => Vec::new(),
            _ => {
                let mut out = Vec::new();
                for (fixed_is_u, fixed, lo, hi, other_lo, other_hi) in [(true, u0, v0, v1, u0, u1), (true, u1, v0, v1, u0, u1), (false, v0, u0, u1, v0, v1), (false, v1, u0, u1, v0, v1)] {
                    if infinite(fixed) || infinite(lo) || infinite(hi) || infinite(other_lo) || infinite(other_hi) {
                        continue;
                    }
                    let at = |t: f64| if fixed_is_u { s.point(fixed, t) } else { s.point(t, fixed) };
                    let spread = self.pole_row_spread(fixed_is_u, fixed == other_lo).unwrap_or_else(|| {
                        let points: Vec<Point> = (0..=SAMPLES).map(|i| at(lo + (hi - lo) * i as f64 / SAMPLES as f64)).collect();
                        diameter(&points)
                    });
                    let point = at(0.5 * (lo + hi));
                    out.push(SingularPlace { point, spread, fixed_is_u, fixed, free_range: (lo, hi) });
                }
                out
            }
        }
    }

    /// For a B-spline (clamped in that direction, not trimmed) or Bézier surface: the
    /// diameter of the boundary row of poles, which bounds the boundary iso's.
    fn pole_row_spread(&self, fixed_is_u: bool, at_first: bool) -> Option<f64> {
        let row = |poles: &[Vec<Point>]| -> Vec<Point> {
            if fixed_is_u {
                let i = if at_first { 0 } else { poles.len() - 1 };
                poles[i].clone()
            } else {
                poles.iter().map(|r| if at_first { r[0] } else { r[r.len() - 1] }).collect()
            }
        };
        match self.surface {
            Surface::Bezier(b) => Some(diameter(&row(b.poles()))),
            Surface::BSpline(b) => {
                let (periodic, mults, degree) = if fixed_is_u { (b.is_u_periodic(), b.u_multiplicities(), b.u_degree()) } else { (b.is_v_periodic(), b.v_multiplicities(), b.v_degree()) };
                let clamped = !periodic && mults.first().is_some_and(|&m| m > degree) && mults.last().is_some_and(|&m| m > degree);
                clamped.then(|| diameter(&row(b.poles())))
            }
            _ => None,
        }
    }

    /// Geometric closure in `u` (computed once; `preci` used only for the verdict).
    pub(crate) fn closed_in_u(&mut self, preci: f64) -> bool {
        if self.closure_u.is_none() {
            self.closure_u = Some(self.closure(true));
        }
        verdict(self.closure_u.expect("closure"), preci)
    }

    /// Geometric closure in `v` (computed once; `preci` used only for the verdict).
    pub(crate) fn closed_in_v(&mut self, preci: f64) -> bool {
        if self.closure_v.is_none() {
            self.closure_v = Some(self.closure(false));
        }
        verdict(self.closure_v.expect("closure"), preci)
    }

    fn closure(&self, in_u: bool) -> Closure {
        let s = self.surface;
        let trimmed = matches!(s, Surface::Trimmed(_));
        let (u0, u1, v0, v1) = self.bounds;
        let (a, b, c, d) = if in_u { (u0, u1, v0, v1) } else { (v0, v1, u0, u1) };
        if infinite(a) || infinite(b) {
            return Closure::Decided(false);
        }
        if !trimmed {
            let (closed, periodic) = if in_u { (s.is_u_closed(), s.is_u_periodic()) } else { (s.is_v_closed(), s.is_v_periodic()) };
            if closed || periodic {
                return Closure::Decided(true);
            }
        }
        match untrimmed(s) {
            Surface::Plane(_) => return Closure::Decided(false),
            Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Extrusion(_) if !in_u => return Closure::Decided(false),
            _ => {}
        }
        if !trimmed {
            match s {
                Surface::BSpline(bs) => {
                    let poles = if in_u { bs.u_pole_count() } else { bs.v_pole_count() };
                    if poles < 3 {
                        return Closure::Decided(false);
                    }
                }
                Surface::Bezier(bz) => {
                    let poles = if in_u { bz.poles().len() } else { bz.poles()[0].len() };
                    if poles < 3 {
                        return Closure::Decided(false);
                    }
                }
                // (A generatrix or basis curve whose ends meet.)
                Surface::Revolution(r) if !in_u => {
                    let gap = r.basis().point(a).distance(r.basis().point(b));
                    return Closure::Measured { gap, middle: f64::INFINITY };
                }
                Surface::Extrusion(e) if in_u => {
                    let gap = e.basis().point(a).distance(e.basis().point(b));
                    return Closure::Measured { gap, middle: f64::INFINITY };
                }
                _ => {}
            }
        }
        // Rows across the other direction (its bounds kept finite), and a B-spline's knots.
        let (c, d) = if infinite(c) || infinite(d) { (-1.0, 1.0) } else { (c, d) };
        let mut rows: Vec<f64> = (0..=ROWS).map(|i| c + (d - c) * i as f64 / ROWS as f64).collect();
        if let Surface::BSpline(bs) = untrimmed(s) {
            let knots = if in_u { bs.v_knots() } else { bs.u_knots() };
            rows.extend(knots.iter().copied().filter(|&k| k > c && k < d));
        }
        let at = |t: f64, row: f64| if in_u { s.point(t, row) } else { s.point(row, t) };
        let mut worst = (f64::NEG_INFINITY, c);
        for &row in &rows {
            let gap = at(a, row).distance(at(b, row));
            let gap = if gap.is_nan() { f64::INFINITY } else { gap };
            if gap > worst.0 {
                worst = (gap, row);
            }
        }
        let middle = at(a, worst.1).distance(at(0.5 * (a + b), worst.1));
        Closure::Measured { gap: worst.0, middle }
    }
}

/// Samples along an iso for its diameter.
const SAMPLES: usize = 16;

/// Rows measured for closure.
const ROWS: usize = 32;

fn verdict(closure: Closure, preci: f64) -> bool {
    match closure {
        Closure::Decided(closed) => closed,
        Closure::Measured { gap, middle } => gap <= preci.max(LINEAR) && gap < middle,
    }
}

/// The largest distance between two of the points.
fn diameter(points: &[Point]) -> f64 {
    let mut d: f64 = 0.0;
    for (i, a) in points.iter().enumerate() {
        for b in &points[i + 1..] {
            let x = a.distance(*b);
            d = if x.is_nan() { f64::INFINITY } else { d.max(x) };
        }
    }
    d
}

#[cfg(test)]
mod tests;
