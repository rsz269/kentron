//! OCCT's `ShapeAnalysis_Surface` (projecting points onto a surface: `ValueOfUV`,
//! `NextValueOfUV`, its singularities, closedness and bound isos) and the point projections
//! of `ShapeAnalysis_Curve` (`Project`, `ProjectAct`, `NextProject`), ported for
//! `ShapeConstruct_ProjectCurveOnSurface` ([`super::pcurve_projection`]).
//!
//! OCCT's `Extrema_ExtPS` (the global search `ValueOfUV` starts from) is
//! [`super::extrema::nearest_on_surface`] here: both find the nearest point of the surface.

use std::f64::consts::PI;

use crate::math::{BoundingBox, Point, Point2, Vector};

use super::curve::Curve;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;
use super::Continuity;

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = CONFUSION * 0.01;
/// `Precision::Infinite()`.
const INFINITE: f64 = 2.0e100;
/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;
/// `gp::Resolution()`.
const RESOLUTION: f64 = f64::MIN_POSITIVE;

/// `Precision::IsInfinite`.
pub(crate) fn is_infinite(value: f64) -> bool {
    value.abs() >= 0.5 * INFINITE
}

/// `Precision::IsPositiveInfinite`.
fn is_positive_infinite(value: f64) -> bool {
    value >= 0.5 * INFINITE
}

/// `Precision::IsNegativeInfinite`.
fn is_negative_infinite(value: f64) -> bool {
    value <= -0.5 * INFINITE
}

/// `RestrictBounds`: infinite ends replaced (by ±1000, or 2000 from a finite end).
fn restrict_bounds(first: &mut f64, last: &mut f64) {
    let f_inf = is_negative_infinite(*first);
    let l_inf = is_positive_infinite(*last);
    if f_inf || l_inf {
        if f_inf && l_inf {
            *first = -1000.0;
            *last = 1000.0;
        } else if f_inf {
            *first = *last - 2000.0;
        } else {
            *last = *first + 2000.0;
        }
    }
}

/// `ShapeAnalysis::AdjustByPeriod`: the multiple of the period to add to `value` to bring
/// it within half a period of `to_value`.
#[doc(alias = "ShapeAnalysis::AdjustByPeriod")]
pub(crate) fn adjust_by_period(value: f64, to_value: f64, period: f64) -> f64 {
    let diff = value - to_value;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

/// `ShapeAnalysis::AdjustToPeriod`: as [`adjust_by_period`], into `[min, max]`.
#[doc(alias = "ShapeAnalysis::AdjustToPeriod")]
pub(crate) fn adjust_to_period(value: f64, min: f64, max: f64) -> f64 {
    adjust_by_period(value, 0.5 * (min + max), max - min)
}

pub(crate) use super::adaptor::in_period;

/// `normalizeAngle` of ElCLib: into `[0, 2π)`.
fn normalize_angle(mut angle: f64) -> f64 {
    while angle < 0.0 {
        angle += 2.0 * PI;
    }
    while angle >= 2.0 * PI {
        angle -= 2.0 * PI;
    }
    angle
}

/// `gp_Vec::Angle`.
fn vector_angle(a: Vector, b: Vector) -> f64 {
    let cosinus = a.dot(b) / (a.length() * b.length());
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = a.cross(b).length() / (a.length() * b.length());
        if cosinus < 0.0 { PI - sinus.asin() } else { sinus.asin() }
    }
}

/// `gp_Vec::AngleWithRef`.
fn angle_with_ref(a: Vector, b: Vector, reference: Vector) -> f64 {
    let angle = vector_angle(a, b);
    if a.cross(b).dot(reference) >= 0.0 { angle } else { -angle }
}

/// `ElCLib::Parameter` for the analytic curves: the parameter of the point of the curve
/// nearest `p` as ElCLib computes it; `None` for other curves.
#[doc(alias = "ElCLib::Parameter")]
fn elclib_parameter(curve: &Curve, p: Point) -> Option<f64> {
    Some(match curve {
        Curve::Line(l) => (p - l.origin).dot(Vector::from(l.direction)),
        Curve::Circle(c) => {
            let f = c.frame();
            let v = p - f.origin();
            if v.dot(v) < RESOLUTION {
                return Some(0.0);
            }
            let dir = Vector::from(f.z());
            let projected = dir.cross(v).cross(dir);
            if projected.dot(projected) < RESOLUTION {
                return Some(0.0);
            }
            normalize_angle(angle_with_ref(Vector::from(f.x()), projected, dir))
        }
        Curve::Ellipse(e) => {
            let f = e.frame();
            let op = p - f.origin();
            let (x, y) = (Vector::from(f.x()), Vector::from(f.y()));
            let (nx, ny) = (op.dot(x), op.dot(y));
            if nx.abs() <= RESOLUTION && ny.abs() <= RESOLUTION {
                return Some(0.0);
            }
            let om = x * nx + y * (ny * (e.major_radius() / e.minor_radius()));
            normalize_angle(angle_with_ref(x, om, Vector::from(f.z())))
        }
        Curve::Hyperbola(h) => ((p - h.frame().origin()).dot(Vector::from(h.frame().y())) / h.minor_radius()).asinh(),
        Curve::Parabola(c) => (p - c.frame().origin()).dot(Vector::from(c.frame().y())),
        _ => return None,
    })
}

/// The curve `GeomAdaptor_Curve` sees (the basis of a trimmed curve).
fn adapted(curve: &Curve) -> &Curve {
    match curve {
        Curve::Trimmed(t) => adapted(t.basis()),
        c => c,
    }
}

/// `IsKind(STANDARD_TYPE(Geom_BoundedCurve))`.
fn is_bounded(curve: &Curve) -> bool {
    matches!(curve, Curve::BSpline(_) | Curve::Bezier(_) | Curve::Trimmed(_))
}

/// `GeomAdaptor_Curve::IsClosed` over `[first, last]`: its ends within confusion.
fn adaptor_is_closed(curve: &Curve, first: f64, last: f64) -> bool {
    if !is_positive_infinite(last) && !is_negative_infinite(first) {
        return curve.point(first).distance(curve.point(last)) <= CONFUSION;
    }
    false
}

/// `Bnd_Box::Distance` from a point.
fn box_distance(b: &BoundingBox, p: Point) -> f64 {
    let gap = |lo: f64, hi: f64, x: f64| if x < lo { lo - x } else if x > hi { x - hi } else { 0.0 };
    let (dx, dy, dz) = (gap(b.min.x, b.max.x, p.x), gap(b.min.y, b.max.y, p.y), gap(b.min.z, b.max.z, p.z));
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// `ProjectOnSegments`: the curve sampled at `segments + 1` points over `[start, end]`, the
/// projection moved to the nearest closer than `distance`, the range narrowed round it.
fn project_on_segments(curve: &Curve, p: Point, segments: usize, start: &mut f64, end: &mut f64, distance: &mut f64, proj: &mut Point, param: &mut f64) {
    if segments == 0 {
        return;
    }
    let step = (*end - *start) / segments as f64;
    let mut min_square = *distance * *distance;
    let mut changed = false;
    for i in 0..=segments {
        let t = *start + step * i as f64;
        let q = curve.point(t);
        let d = q.distance_squared(p);
        if d < min_square {
            min_square = d;
            *proj = q;
            *param = t;
            changed = true;
        }
    }
    if changed {
        *distance = min_square.sqrt();
    }
    *end = end.min(*param + step);
    *start = start.max(*param - step);
}

/// A point projected onto a curve: its distance, the point and its parameter.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CurveProjection {
    pub distance: f64,
    pub point: Point,
    pub parameter: f64,
}

/// `ShapeAnalysis_Curve::Project` (the `Geom_Curve` form) over `[cf, cl]`.
#[doc(alias = "ShapeAnalysis_Curve::Project")]
pub(crate) fn project_on_curve(curve: &Curve, p: Point, preci: f64, cf: f64, cl: f64, adjust_to_ends: bool) -> CurveProjection {
    let mut u_min = if cf < cl { cf } else { cl };
    let mut u_max = if cf < cl { cl } else { cf };
    if is_bounded(curve) {
        let prec = if adjust_to_ends { preci } else { CONFUSION };
        let low = curve.point(u_min);
        let high = curve.point(u_max);
        let d = low.distance(p);
        if d <= prec {
            return CurveProjection { distance: d, point: low, parameter: u_min };
        }
        let d = high.distance(p);
        if d <= prec {
            return CurveProjection { distance: d, point: high, parameter: u_max };
        }
    }
    if !curve.is_closed() {
        let delta = super::resolution::curve_resolution(curve, preci).min((u_max - u_min) * 0.1);
        u_min -= delta;
        u_max += delta;
    }
    project_act(curve, u_min, u_max, p, preci)
}

/// `ShapeAnalysis_Curve::Project` (the `Adaptor3d_Curve` form) on the curve over
/// `[first, last]`.
fn project_on_adaptor(curve: &Curve, first: f64, last: f64, p: Point, preci: f64, adjust_to_ends: bool) -> CurveProjection {
    if is_infinite(first) && is_infinite(last) {
        return project_act(curve, first, last, p, preci);
    }
    let prec = if adjust_to_ends { preci } else { CONFUSION };
    let low = curve.point(first);
    let high = curve.point(last);
    let (dl, dh) = (low.distance(p), high.distance(p));
    if dl <= prec {
        return CurveProjection { distance: dl, point: low, parameter: first };
    }
    if dh <= prec {
        return CurveProjection { distance: dh, point: high, parameter: last };
    }
    let projected = project_act(curve, first, last, p, preci);
    if projected.distance < dl + CONFUSION && projected.distance < dh + CONFUSION {
        return projected;
    }
    if dl < dh {
        return CurveProjection { distance: dl, point: low, parameter: first };
    }
    CurveProjection { distance: dh, point: high, parameter: last }
}

/// `ShapeAnalysis_Curve::ProjectAct` on the curve over `[u_min, u_max]`.
#[doc(alias = "ShapeAnalysis_Curve::ProjectAct")]
fn project_act(curve: &Curve, mut u_min: f64, mut u_max: f64, p: Point, tolerance: f64) -> CurveProjection {
    // (The adaptor, over the range given.)
    let adaptor = super::adaptor::GeomAdaptor::new(curve, u_min, u_max);
    let mut ok = false;
    let mut param = 0.0;
    let mut proj = Point::ORIGIN;
    // (Extrema_ExtPC on the adaptor: its nearest minimum.)
    let search = super::extrema::PointCurveSearch::with_domain(curve, u_min, u_max);
    if let Ok(found) = search.find(p, 1e-10) {
        let mut min = REAL_LAST;
        let mut best = None;
        for e in found.iter() {
            if !e.is_minimum {
                continue;
            }
            if e.square_distance < min {
                min = e.square_distance;
                best = Some((e.parameter, e.point));
            }
        }
        if let Some((t, q)) = best {
            param = t;
            proj = q;
            ok = true;
        }
    }
    let computed = (param, proj);
    let mut closed = false;
    let mut period = 0.0;
    let mut distance;
    let mut mod_min = INFINITE;
    let mut old = None;
    if ok {
        old = Some((param, proj));
        distance = proj.distance(p);
        mod_min = distance;
        if distance > tolerance {
            ok = false;
        }
        if adaptor_is_closed(curve, u_min, u_max) {
            closed = true;
            period = u_max - u_min;
        }
    }
    if !ok {
        param = 0.0;
        let basis = adapted(curve);
        match basis {
            Curve::Circle(c) => {
                proj = c.frame().origin();
                if c.radius() <= RESOLUTION || p.distance_squared(proj) <= RESOLUTION {
                    param = u_min;
                    proj = proj + Vector::from(c.frame().x()) * c.radius();
                } else {
                    param = elclib_parameter(basis, p).expect("a circle");
                    proj = basis.point(param);
                }
                closed = true;
                period = 2.0 * PI;
            }
            Curve::Hyperbola(_) | Curve::Parabola(_) | Curve::Line(_) => {
                param = elclib_parameter(basis, p).expect("a conic");
                proj = basis.point(param);
            }
            Curve::Ellipse(_) => {
                param = elclib_parameter(basis, p).expect("an ellipse");
                proj = basis.point(param);
                closed = true;
                period = 2.0 * PI;
            }
            _ => {
                distance = INFINITE;
                project_on_segments(curve, p, 25, &mut u_min, &mut u_max, &mut distance, &mut proj, &mut param);
                if distance <= tolerance {
                    return CurveProjection { distance, point: proj, parameter: param };
                }
                if let Some(e) = super::locate_ext_pc::locate_ext_pc(p, &adaptor, param, u_min, u_max, tolerance) {
                    param = e.parameter;
                    proj = e.point;
                    let newton = p.distance(proj);
                    if newton < mod_min {
                        return CurveProjection { distance: newton, point: proj, parameter: param };
                    }
                }
                for segments in [40, 20, 25, 40] {
                    project_on_segments(curve, p, segments, &mut u_min, &mut u_max, &mut distance, &mut proj, &mut param);
                    if distance <= tolerance {
                        return CurveProjection { distance, point: proj, parameter: param };
                    }
                }
                if distance > mod_min {
                    distance = mod_min;
                    param = computed.0;
                    proj = computed.1;
                }
                return CurveProjection { distance, point: proj, parameter: param };
            }
        }
    }
    if closed && (param < u_min || param > u_max) {
        param += adjust_by_period(param, 0.5 * (u_min + u_max), period);
    }
    if let Some((old_param, old_proj)) = old {
        if old_proj.distance_squared(p) < proj.distance_squared(p) {
            proj = old_proj;
            param = old_param;
        }
    }
    CurveProjection { distance: proj.distance(p), point: proj, parameter: param }
}

/// `ShapeAnalysis_Curve::NextProject` (the `Geom_Curve` form): the projection searched
/// from `previous` first.
#[doc(alias = "ShapeAnalysis_Curve::NextProject")]
pub(crate) fn next_project_on_curve(previous: f64, curve: &Curve, p: Point, preci: f64, cf: f64, cl: f64, adjust_to_ends: bool) -> CurveProjection {
    let mut u_min = if cf < cl { cf } else { cl };
    let mut u_max = if cf < cl { cl } else { cf };
    if is_bounded(curve) {
        let prec = if adjust_to_ends { preci } else { CONFUSION };
        let low = curve.point(u_min);
        let high = curve.point(u_max);
        let d = low.distance(p);
        if d <= prec {
            return CurveProjection { distance: d, point: low, parameter: u_min };
        }
        let d = high.distance(p);
        if d <= prec {
            return CurveProjection { distance: d, point: high, parameter: u_max };
        }
    }
    if !curve.is_closed() {
        let delta = super::resolution::curve_resolution(curve, preci).min((u_max - u_min) * 0.1);
        u_min -= delta;
        u_max += delta;
    }
    let adaptor = super::adaptor::GeomAdaptor::new(curve, u_min, u_max);
    if let Some(e) = super::locate_ext_pc::locate_ext_pc(p, &adaptor, previous, u_min, u_max, preci) {
        return CurveProjection { distance: p.distance(e.point), point: e.point, parameter: e.parameter };
    }
    project_on_adaptor(curve, u_min, u_max, p, preci, false)
}

/// The kind of surface `GeomAdaptor_Surface` sees (`GetType`): a trimmed surface's basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Bezier,
    BSpline,
    Extrusion,
    Revolution,
    Offset,
}

/// The surface `GeomAdaptor_Surface` evaluates: the basis of a trimmed surface.
fn adapted_surface(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => adapted_surface(t.basis()),
        s => s,
    }
}

fn kind(surface: &Surface) -> Kind {
    match adapted_surface(surface) {
        Surface::Plane(_) => Kind::Plane,
        Surface::Cylinder(_) => Kind::Cylinder,
        Surface::Cone(_) => Kind::Cone,
        Surface::Sphere(_) => Kind::Sphere,
        Surface::Torus(_) => Kind::Torus,
        Surface::Bezier(_) => Kind::Bezier,
        Surface::BSpline(_) => Kind::BSpline,
        Surface::Extrusion(_) => Kind::Extrusion,
        Surface::Revolution(_) => Kind::Revolution,
        Surface::Offset(_) => Kind::Offset,
        Surface::Trimmed(_) => unreachable!("a trimmed surface's basis"),
    }
}

/// A singular place of a surface (a pole, an apex, a degenerate side): `ShapeAnalysis_Surface`'s
/// `Singularity`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Singularity {
    /// How far the place spreads in 3D (0 for a true pole).
    pub preci: f64,
    pub point: Point,
    pub first_uv: Point2,
    pub last_uv: Point2,
    pub first_parameter: f64,
    pub last_parameter: f64,
    /// True if the place is a u-iso (the parameter running along it is v).
    pub u_iso: bool,
}

/// OCCT's `ShapeAnalysis_Surface`: a surface with what is learned about it while points are
/// projected onto it (its singularities, closedness, bound isos), and the gap of the last
/// projection.
#[doc(alias = "ShapeAnalysis_Surface")]
pub(crate) struct SurfaceAnalysis<'a> {
    surface: &'a Surface,
    uf: f64,
    ul: f64,
    vf: f64,
    vl: f64,
    gap: f64,
    singularities: Option<Vec<Singularity>>,
    u_close_value: f64,
    v_close_value: f64,
    u_delta: f64,
    v_delta: f64,
    /// The isos at the four bounds (uf, ul, vf, vl), once computed.
    isos: Option<[Option<Curve>; 4]>,
    /// Their boxes.
    iso_boxes: Option<[Option<BoundingBox>; 4]>,
}

impl<'a> SurfaceAnalysis<'a> {
    pub(crate) fn new(surface: &'a Surface) -> Self {
        let (uf, ul, vf, vl) = surface.bounds();
        SurfaceAnalysis {
            surface,
            uf,
            ul,
            vf,
            vl,
            gap: 0.0,
            singularities: None,
            u_close_value: -1.0,
            v_close_value: -1.0,
            u_delta: 0.01,
            v_delta: 0.01,
            isos: None,
            iso_boxes: None,
        }
    }

    pub(crate) fn surface(&self) -> &'a Surface {
        self.surface
    }

    /// `Bounds`.
    pub(crate) fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.uf, self.ul, self.vf, self.vl)
    }

    /// The gap of the last projection (`Gap`).
    pub(crate) fn gap(&self) -> f64 {
        self.gap
    }

    fn value(&self, uv: Point2) -> Point {
        self.surface.point(uv.x, uv.y)
    }

    fn u_resolution(&self, r3d: f64) -> f64 {
        super::resolution::u_resolution(self.surface, self.surface.bounds(), r3d)
    }

    fn v_resolution(&self, r3d: f64) -> f64 {
        super::resolution::v_resolution(self.surface, self.surface.bounds(), r3d)
    }

    /// `ComputeSingularities` (once).
    fn singularities(&mut self) -> &[Singularity] {
        if self.singularities.is_none() {
            let found = self.compute_singularities();
            self.singularities = Some(found);
        }
        self.singularities.as_deref().expect("computed")
    }

    #[doc(alias = "ShapeAnalysis_Surface::ComputeSingularities")]
    fn compute_singularities(&self) -> Vec<Singularity> {
        let (su1, su2, sv1, sv2) = self.bounds();
        let s = self.surface;
        let mut out: Vec<Singularity> = Vec::new();
        let p2 = Point2::new;
        match s {
            Surface::Cone(cone) => {
                let v_apex = -cone.radius() / cone.semi_angle().sin();
                out.push(Singularity { preci: 0.0, point: cone.apex(), first_uv: p2(su1, v_apex), last_uv: p2(su2, v_apex), first_parameter: su1, last_parameter: su2, u_iso: false });
            }
            Surface::Torus(torus) => {
                let minor = torus.minor_radius();
                let major = torus.major_radius();
                let angle = (1.0f64.min(major / minor)).acos();
                let preci = 0.0f64.max(major - minor);
                let first = Singularity { preci, point: s.point(0.0, PI - angle), first_uv: p2(su1, PI - angle), last_uv: p2(su2, PI - angle), first_parameter: su1, last_parameter: su2, u_iso: false };
                let second = Singularity { preci, point: s.point(0.0, PI + angle), first_uv: p2(su2, PI + angle), last_uv: p2(su1, PI + angle), first_parameter: su1, last_parameter: su2, u_iso: false };
                out.push(first);
                if major <= minor {
                    out.push(second);
                }
            }
            Surface::Sphere(_) => {
                // (The northern pole first.)
                out.push(Singularity { preci: 0.0, point: s.point(su1, sv2), first_uv: p2(su2, sv2), last_uv: p2(su1, sv2), first_parameter: su1, last_parameter: su2, u_iso: false });
                out.push(Singularity { preci: 0.0, point: s.point(su1, sv1), first_uv: p2(su1, sv1), last_uv: p2(su2, sv1), first_parameter: su1, last_parameter: su2, u_iso: false });
            }
            Surface::Bezier(_) | Surface::BSpline(_) | Surface::Trimmed(_) | Surface::Revolution(_) | Surface::Offset(_) => {
                let p0 = s.point(su1, 0.5 * (sv1 + sv2));
                let p1 = s.point(su2, 0.5 * (sv1 + sv2));
                let p2_ = s.point(0.5 * (su1 + su2), sv1);
                let p3 = s.point(0.5 * (su1 + su2), sv2);
                let c1 = s.point(su1, sv1);
                let c2 = s.point(su1, sv2);
                let c3 = s.point(su2, sv1);
                let c4 = s.point(su2, sv2);
                let preci0 = c1.distance(c2).max(p0.distance(c1).max(p0.distance(c2)));
                let preci1 = c3.distance(c4).max(p1.distance(c3).max(p1.distance(c4)));
                let preci2 = c1.distance(c3).max(p2_.distance(c1).max(p2_.distance(c3)));
                let preci3 = c2.distance(c4).max(p3.distance(c2).max(p3.distance(c4)));
                out.push(Singularity { preci: preci0, point: p0, first_uv: p2(su1, sv2), last_uv: p2(su1, sv1), first_parameter: sv1, last_parameter: sv2, u_iso: true });
                out.push(Singularity { preci: preci1, point: p1, first_uv: p2(su2, sv1), last_uv: p2(su2, sv2), first_parameter: sv1, last_parameter: sv2, u_iso: true });
                out.push(Singularity { preci: preci2, point: p2_, first_uv: p2(su1, sv1), last_uv: p2(su2, sv1), first_parameter: su1, last_parameter: su2, u_iso: false });
                out.push(Singularity { preci: preci3, point: p3, first_uv: p2(su2, sv2), last_uv: p2(su1, sv2), first_parameter: su1, last_parameter: su2, u_iso: false });
            }
            _ => {}
        }
        // `SortSingularities`: a selection sort by `preci`.
        let n = out.len();
        for i in 0..n.saturating_sub(1) {
            let mut min_preci = out[i].preci;
            let mut min_index = i;
            for j in i + 1..n {
                if min_preci > out[j].preci {
                    min_preci = out[j].preci;
                    min_index = j;
                }
            }
            if min_index != i {
                out.swap(min_index, i);
            }
        }
        out
    }

    /// `NbSingularities`.
    pub(crate) fn nb_singularities(&mut self, preci: f64) -> usize {
        self.singularities().iter().filter(|s| s.preci <= preci).count()
    }

    /// `Singularity(num)` (from 0).
    pub(crate) fn singularity(&mut self, index: usize) -> Option<Singularity> {
        self.singularities().get(index).copied()
    }

    /// `IsDegenerated(P3d, preci)`.
    pub(crate) fn is_degenerated(&mut self, p: Point, preci: f64) -> bool {
        let singularities = self.singularities().to_vec();
        for s in singularities.iter().take_while(|s| s.preci <= preci) {
            self.gap = s.point.distance(p);
            if self.gap <= preci {
                return true;
            }
        }
        false
    }

    /// `ProjectDegenerated(nbrPnt, points, pnt2d, preci, direct)`: the points at an end of a
    /// run lying on a singular place given the other parameter of the first point off it.
    pub(crate) fn project_degenerated(&mut self, points: &[Point], uvs: &mut [Point2], preci: f64, direct: bool) -> bool {
        let n = points.len();
        let singularities = self.singularities().to_vec();
        let step: isize = if direct { 1 } else { -1 };
        let mut index = None;
        let mut gap_min = REAL_LAST;
        let prec2 = preci * preci;
        let j = if direct { 0 } else { n - 1 };
        for (i, s) in singularities.iter().enumerate().take_while(|(_, s)| s.preci <= preci) {
            let mut gap2 = s.point.distance_squared(points[j]);
            if gap2 > prec2 {
                gap2 = gap2.min(s.point.distance_squared(self.value(uvs[j])));
            }
            if gap2 <= prec2 && gap_min > gap2 {
                gap_min = gap2;
                index = Some(i);
            }
        }
        let Some(index) = index else { return false };
        self.gap = gap_min.sqrt();
        let s = singularities[index];
        let mut k = j as isize + step;
        let mut pk = Point2::new(0.0, 0.0);
        while k >= 0 && (k as usize) < n {
            pk = uvs[k as usize];
            if s.point.distance_squared(points[k as usize]) > prec2 && s.point.distance_squared(self.value(pk)) > prec2 {
                break;
            }
            k += step;
        }
        if k < 0 || k as usize >= n {
            // (The whole pcurve degenerate: spread evenly.)
            let x1 = if s.u_iso { uvs[0].y } else { uvs[0].x };
            let x2 = if s.u_iso { uvs[n - 1].y } else { uvs[n - 1].x };
            for (jj, uv) in uvs.iter_mut().enumerate() {
                let x = (x1 * (n - 1 - jj) as f64 + x2 * jj as f64) / (n - 1) as f64;
                if !s.u_iso {
                    uv.x = x;
                } else {
                    uv.y = x;
                }
            }
            return true;
        }
        let mut jj = k - step;
        while jj >= 0 && (jj as usize) < n {
            if !s.u_iso {
                uvs[jj as usize].x = pk.x;
            } else {
                uvs[jj as usize].y = pk.y;
            }
            jj -= step;
        }
        true
    }

    /// `ComputeBoundIsos` (once).
    fn bound_isos(&mut self) -> &[Option<Curve>; 4] {
        if self.isos.is_none() {
            let s = self.surface;
            let iso = |u_type: bool, at: f64| -> Option<Curve> {
                // (Infinite bounds give no iso.)
                if is_infinite(at) {
                    return None;
                }
                Some(geom_iso(s, u_type, at))
            };
            self.isos = Some([iso(true, self.uf), iso(true, self.ul), iso(false, self.vf), iso(false, self.vl)]);
        }
        self.isos.as_ref().expect("computed")
    }

    /// `UIso(U)`.
    pub(crate) fn u_iso(&mut self, u: f64) -> Option<Curve> {
        if u == self.uf {
            return self.bound_isos()[0].clone();
        }
        if u == self.ul {
            return self.bound_isos()[1].clone();
        }
        Some(geom_iso(self.surface, true, u))
    }

    /// `VIso(V)`.
    pub(crate) fn v_iso(&mut self, v: f64) -> Option<Curve> {
        if v == self.vf {
            return self.bound_isos()[2].clone();
        }
        if v == self.vl {
            return self.bound_isos()[3].clone();
        }
        Some(geom_iso(self.surface, false, v))
    }

    /// `ComputeBoxes`: the boxes of the bound isos (uf, ul, vf, vl).
    pub(crate) fn iso_boxes(&mut self) -> [Option<BoundingBox>; 4] {
        if self.iso_boxes.is_none() {
            let isos = self.bound_isos().clone();
            let boxes = isos.map(|iso| iso.map(|c| super::bounds::curve(&c, CONFUSION)));
            self.iso_boxes = Some(boxes);
        }
        self.iso_boxes.clone().expect("computed")
    }

    /// `IsUClosed(preci)`.
    #[doc(alias = "ShapeAnalysis_Surface::IsUClosed")]
    pub(crate) fn is_u_closed(&mut self, preci: f64) -> bool {
        let prec = preci.max(CONFUSION);
        let mut mid = -1.0;
        if self.u_close_value < 0.0 {
            let (mut uf, mut ul, mut vf, mut vl) = self.bounds();
            restrict_bounds(&mut uf, &mut ul);
            restrict_bounds(&mut vf, &mut vl);
            self.u_delta = (ul - uf).abs() / 20.0;
            if self.surface.is_u_closed() {
                self.u_close_value = 0.0;
                self.u_delta = 0.0;
                self.gap = 0.0;
                return true;
            }
            let s = self.surface;
            let trimmed = matches!(s, Surface::Trimmed(_));
            match (trimmed, adapted_surface(s)) {
                (false, Surface::Plane(_)) => self.u_close_value = REAL_LAST,
                (false, Surface::Extrusion(e)) => {
                    let c = e.basis();
                    let (f, l) = (c.first_parameter(), c.last_parameter());
                    if !is_infinite(f) && !is_infinite(l) {
                        let p1 = c.point(f);
                        let p2 = c.point(l);
                        self.u_close_value = p1.distance_squared(p2);
                        let pm = c.point((f + l) / 2.0);
                        mid = p1.distance_squared(pm);
                    } else {
                        self.u_close_value = REAL_LAST;
                    }
                }
                (false, Surface::BSpline(bs)) => {
                    let nbup = bs.poles().len();
                    let mut distmin;
                    if bs.is_u_periodic() {
                        self.u_close_value = 0.0;
                        self.u_delta = 0.0;
                    } else if nbup < 3 {
                        self.u_close_value = REAL_LAST;
                    } else if bs.is_u_rational()
                        || bs.u_multiplicities()[0] != bs.u_degree() + 1
                        || *bs.u_multiplicities().last().expect("knots") != bs.u_degree() + 1
                    {
                        let vk = bs.v_knots();
                        let mut v = vk[0];
                        let p1 = s.point(uf, v);
                        let p2 = s.point(ul, v);
                        self.u_close_value = p1.distance_squared(p2);
                        let pm = s.point((uf + ul) / 2.0, v);
                        mid = p1.distance_squared(pm);
                        distmin = self.u_close_value;
                        for i in 1..vk.len() {
                            v = 0.5 * (vk[i - 1] + vk[i]);
                            let p1 = s.point(uf, v);
                            let p2 = s.point(ul, v);
                            let d = p1.distance_squared(p2);
                            if d > self.u_close_value {
                                self.u_close_value = d;
                                let pm = s.point((uf + ul) / 2.0, v);
                                mid = p1.distance_squared(pm);
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.u_delta = self.u_delta.min(0.5 * self.u_resolution(distmin));
                    } else {
                        let nbvp = bs.poles()[0].len();
                        self.u_close_value = bs.pole(0, 0).distance_squared(bs.pole(nbup - 1, 0));
                        mid = bs.pole(0, 0).distance_squared(bs.pole(nbup / 2, 0));
                        distmin = self.u_close_value;
                        for i in 1..nbvp {
                            let d = bs.pole(0, i).distance_squared(bs.pole(nbup - 1, i));
                            if d > self.u_close_value {
                                self.u_close_value = d;
                                mid = bs.pole(0, i).distance_squared(bs.pole(nbup / 2, i));
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.u_delta = self.u_delta.min(0.5 * self.u_resolution(distmin));
                    }
                }
                (false, Surface::Bezier(bz)) => {
                    let nbup = bz.poles().len();
                    if nbup < 3 {
                        self.u_close_value = REAL_LAST;
                    } else {
                        let nbvp = bz.poles()[0].len();
                        self.u_close_value = bz.pole(0, 0).distance_squared(bz.pole(nbup - 1, 0));
                        mid = bz.pole(0, 0).distance_squared(bz.pole(nbup / 2, 0));
                        let mut distmin = self.u_close_value;
                        for i in 0..nbvp {
                            let d = bz.pole(0, i).distance_squared(bz.pole(nbup - 1, i));
                            if d > self.u_close_value {
                                self.u_close_value = d;
                                mid = bz.pole(0, i).distance_squared(bz.pole(nbup / 2, i));
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.u_delta = self.u_delta.min(0.5 * self.u_resolution(distmin));
                    }
                }
                _ => {
                    // (Trimmed and offset surfaces, and the rest: sampled.)
                    let points = 101;
                    let p1 = s.point(uf, vf);
                    let p2 = s.point(ul, vf);
                    self.u_close_value = p1.distance_squared(p2);
                    let pm = s.point((uf + ul) / 2.0, vf);
                    mid = p1.distance_squared(pm);
                    let mut distmin = self.u_close_value;
                    for i in 1..points {
                        let v = vf + (vl - vf) * i as f64 / (points - 1) as f64;
                        let p1 = s.point(uf, v);
                        let p2 = s.point(ul, v);
                        let d = p1.distance_squared(p2);
                        if d > self.u_close_value {
                            self.u_close_value = d;
                            let pm = s.point((uf + ul) / 2.0, v);
                            mid = p1.distance_squared(pm);
                        } else {
                            distmin = distmin.min(d);
                        }
                    }
                    distmin = distmin.sqrt();
                    self.u_delta = self.u_delta.min(0.5 * self.u_resolution(distmin));
                }
            }
            self.gap = self.u_close_value.sqrt();
            self.u_close_value = self.gap;
        }
        if mid > 0.0 && self.u_close_value > mid.sqrt() {
            self.u_close_value = REAL_LAST;
            return false;
        }
        self.u_close_value <= prec
    }

    /// `IsVClosed(preci)`.
    #[doc(alias = "ShapeAnalysis_Surface::IsVClosed")]
    pub(crate) fn is_v_closed(&mut self, preci: f64) -> bool {
        let prec = preci.max(CONFUSION);
        let mut mid = -1.0;
        if self.v_close_value < 0.0 {
            let (mut uf, mut ul, mut vf, mut vl) = self.bounds();
            restrict_bounds(&mut uf, &mut ul);
            restrict_bounds(&mut vf, &mut vl);
            self.v_delta = (vl - vf).abs() / 20.0;
            if self.surface.is_v_closed() {
                self.v_close_value = 0.0;
                self.v_delta = 0.0;
                self.gap = 0.0;
                return true;
            }
            let s = self.surface;
            let trimmed = matches!(s, Surface::Trimmed(_));
            match (trimmed, adapted_surface(s)) {
                (false, Surface::Plane(_) | Surface::Cone(_) | Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Extrusion(_)) => {
                    self.v_close_value = REAL_LAST;
                }
                (false, Surface::Revolution(r)) => {
                    let c = r.basis();
                    let p1 = c.point(c.first_parameter());
                    let p2 = c.point(c.last_parameter());
                    self.v_close_value = p1.distance_squared(p2);
                }
                (false, Surface::BSpline(bs)) => {
                    let nbvp = bs.poles()[0].len();
                    let mut distmin;
                    if bs.is_v_periodic() {
                        self.v_close_value = 0.0;
                        self.v_delta = 0.0;
                    } else if nbvp < 3 {
                        self.v_close_value = REAL_LAST;
                    } else if bs.is_v_rational()
                        || bs.v_multiplicities()[0] != bs.v_degree() + 1
                        || *bs.v_multiplicities().last().expect("knots") != bs.v_degree() + 1
                    {
                        let uk = bs.u_knots();
                        let mut u = uk[0];
                        let p1 = s.point(u, vf);
                        let p2 = s.point(u, vl);
                        self.v_close_value = p1.distance_squared(p2);
                        let pm = s.point(u, (vf + vl) / 2.0);
                        mid = p1.distance_squared(pm);
                        distmin = self.v_close_value;
                        for i in 1..uk.len() {
                            u = 0.5 * (uk[i - 1] + uk[i]);
                            let p1 = s.point(u, vf);
                            let p2 = s.point(u, vl);
                            let d = p1.distance_squared(p2);
                            if d > self.v_close_value {
                                self.v_close_value = d;
                                let pm = s.point(u, (vf + vl) / 2.0);
                                mid = p1.distance_squared(pm);
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.v_delta = self.v_delta.min(0.5 * self.v_resolution(distmin));
                    } else {
                        let nbup = bs.poles().len();
                        self.v_close_value = bs.pole(0, 0).distance_squared(bs.pole(0, nbvp - 1));
                        mid = bs.pole(0, 0).distance_squared(bs.pole(0, nbvp / 2));
                        distmin = self.v_close_value;
                        for i in 1..nbup {
                            let d = bs.pole(i, 0).distance_squared(bs.pole(i, nbvp - 1));
                            if d > self.v_close_value {
                                self.v_close_value = d;
                                mid = bs.pole(i, 0).distance_squared(bs.pole(i, nbvp / 2));
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.v_delta = self.v_delta.min(0.5 * self.v_resolution(distmin));
                    }
                }
                (false, Surface::Bezier(bz)) => {
                    let nbvp = bz.poles()[0].len();
                    if nbvp < 3 {
                        self.v_close_value = REAL_LAST;
                    } else {
                        let nbup = bz.poles().len();
                        self.v_close_value = bz.pole(0, 0).distance_squared(bz.pole(0, nbvp - 1));
                        mid = bz.pole(0, 0).distance_squared(bz.pole(0, nbvp / 2));
                        let mut distmin = self.v_close_value;
                        for i in 1..nbup {
                            let d = bz.pole(i, 0).distance_squared(bz.pole(i, nbvp - 1));
                            if d > self.v_close_value {
                                self.v_close_value = d;
                                mid = bz.pole(i, 0).distance_squared(bz.pole(i, nbvp / 2));
                            } else {
                                distmin = distmin.min(d);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.v_delta = self.v_delta.min(0.5 * self.v_resolution(distmin));
                    }
                }
                _ => {
                    let points = 101;
                    let p1 = s.point(uf, vf);
                    let p2 = s.point(uf, vl);
                    let pm = s.point(uf, (vf + vl) / 2.0);
                    self.v_close_value = p1.distance_squared(p2);
                    mid = p1.distance_squared(pm);
                    let mut distmin = self.v_close_value;
                    for i in 1..points {
                        let u = uf + (ul - uf) * i as f64 / (points - 1) as f64;
                        let p1 = s.point(u, vf);
                        let p2 = s.point(u, vl);
                        let d = p1.distance_squared(p2);
                        if d > self.v_close_value {
                            self.v_close_value = d;
                            let pm = s.point(u, (vf + vl) / 2.0);
                            mid = p1.distance_squared(pm);
                        } else {
                            distmin = distmin.min(d);
                        }
                    }
                    distmin = distmin.sqrt();
                    self.v_delta = self.v_delta.min(0.5 * self.v_resolution(distmin));
                }
            }
            self.gap = self.v_close_value.sqrt();
            self.v_close_value = self.gap;
        }
        if mid > 0.0 && self.v_close_value > mid.sqrt() {
            self.v_close_value = REAL_LAST;
            return false;
        }
        self.v_close_value <= prec
    }

    /// `SurfaceNewton`: Newton steps from `previous` toward the foot of `p`; 0 if they fail,
    /// 2 where the surface's normal is nearly degenerate there, else 1, with the solution.
    #[doc(alias = "ShapeAnalysis_Surface::SurfaceNewton")]
    fn surface_newton(&self, previous: Point2, p: Point, preci: f64) -> (i32, Point2) {
        let (uf, ul, vf, vl) = self.bounds();
        let du = self.u_resolution(preci);
        let dv = self.v_resolution(preci);
        let (u_lo, u_hi, v_lo, v_hi) = (uf - du, ul + du, vf - dv, vl + dv);
        const TOL: f64 = CONFUSION;
        const TOL2: f64 = TOL * TOL;
        let (mut u, mut v) = (previous.x, previous.y);
        let rs_first = p - self.surface.point(u, v);
        for _ in 0..25 {
            let d2 = self.surface.d2(u, v);
            let (ru, rv, ruu, rvv, ruv) = (d2.du, d2.dv, d2.duu, d2.dvv, d2.duv);
            let ru2 = ru.dot(ru);
            let rv2 = rv.dot(rv);
            let n = ru.cross(rv);
            let nrm2 = n.dot(n);
            if nrm2 < 1e-10 || is_positive_infinite(nrm2) {
                break;
            }
            let rs = p - self.surface.point(u, v);
            let r_suu = rs.dot(ruu);
            let r_svv = rs.dot(rvv);
            let r_suv = rs.dot(ruv);
            let d = -nrm2 + rv2 * r_suu + ru2 * r_svv - 2.0 * r_suv * ru.dot(rv) + r_suv * r_suv - r_suu * r_svv;
            if d.abs() < 1e-10 {
                break;
            }
            let fract = 1.0 / d;
            let du = rs.dot(n.cross(rv) + ru * r_svv - rv * r_suv) * fract;
            let dv = rs.dot(ru.cross(n) + rv * r_suu - ru * r_suv) * fract;
            u += du;
            v += dv;
            if u < u_lo || u > u_hi || v < v_lo || v > v_hi {
                break;
            }
            let resolution = 1e-12f64.max((u + v) * 10e-16);
            if du.abs() + dv.abs() > resolution {
                continue;
            }
            let rs2 = rs.dot(rs);
            if rs2 > rs_first.dot(rs_first) {
                break;
            }
            let rsn = rs.dot(n);
            if rs2 - rsn * rsn / nrm2 > TOL2 {
                break;
            }
            return (if nrm2 < 0.01 * ru2 * rv2 { 2 } else { 1 }, Point2::new(u, v));
        }
        (0, Point2::new(u, v))
    }

    /// `NextValueOfUV`: the parameters of `p`, searched from `previous` first (on free-form
    /// surfaces).
    #[doc(alias = "ShapeAnalysis_Surface::NextValueOfUV")]
    pub(crate) fn next_value_of_uv(&mut self, previous: Point2, p: Point, preci: f64, max_preci: f64) -> Point2 {
        let k = kind(self.surface);
        if matches!(k, Kind::Bezier | Kind::BSpline | Kind::Extrusion | Kind::Revolution | Kind::Offset) {
            if let Surface::BSpline(bs) = adapted_surface(self.surface) {
                // (Near a C0 knot: the full search.)
                if u_continuity_is_c0(bs, self.uf, self.ul, true) && bs.u_knots().iter().any(|&k| (k - previous.x).abs() < CONFUSION) {
                    return self.value_of_uv(p, preci);
                }
                if u_continuity_is_c0(bs, self.vf, self.vl, false) && bs.v_knots().iter().any(|&k| (k - previous.y).abs() < CONFUSION) {
                    return self.value_of_uv(p, preci);
                }
            }
            let (res, sol) = self.surface_newton(previous, p, preci);
            if res != 0 {
                let gap = p.distance(self.value(sol));
                if res == 2 || (max_preci > 0.0 && gap - max_preci > CONFUSION) {
                    let (mut u, mut v) = (sol.x, sol.y);
                    self.gap = self.uv_from_iso(p, preci, &mut u, &mut v);
                    if gap >= self.gap {
                        return Point2::new(u, v);
                    }
                }
                self.gap = gap;
                return sol;
            }
        }
        self.value_of_uv(p, preci)
    }

    /// `ValueOfUV`: the parameters of `p` on the surface, its gap kept.
    #[doc(alias = "ShapeAnalysis_Surface::ValueOfUV")]
    pub(crate) fn value_of_uv(&mut self, p: Point, preci: f64) -> Point2 {
        let (mut s, mut t);
        self.gap = -1.0;
        let (mut uf, mut ul, mut vf, mut vl) = self.bounds();
        let k = kind(self.surface);
        let basis = adapted_surface(self.surface);
        match k {
            Kind::Plane => {
                (s, t) = basis.parameters(p).expect("a plane");
            }
            Kind::Cylinder | Kind::Cone | Kind::Sphere => {
                (s, t) = basis.parameters(p).expect("an elementary surface");
                s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * PI);
            }
            Kind::Torus => {
                (s, t) = basis.parameters(p).expect("a torus");
                s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * PI);
                t += adjust_by_period(t, 0.5 * (vf + vl), 2.0 * PI);
            }
            Kind::Bezier | Kind::BSpline | Kind::Extrusion | Kind::Revolution | Kind::Offset => {
                s = (uf + ul) / 2.0;
                t = (vf + vl) / 2.0;
                if k == Kind::Extrusion && is_infinite(uf) && is_infinite(ul) {
                    let (res, solution) = self.surface_newton(Point2::new(s, t), p, preci);
                    if res != 0 {
                        return solution;
                    }
                    uf = -500.0;
                    ul = 500.0;
                }
                restrict_bounds(&mut uf, &mut ul);
                restrict_bounds(&mut vf, &mut vl);
                if self.u_close_value < 0.0 {
                    self.is_u_closed(0.0);
                }
                if self.v_close_value < 0.0 {
                    self.is_v_closed(0.0);
                }
                // (Extrema_ExtPS over the bounds, widened by the resolution: the nearest point.)
                let (es, et, ed) = super::extrema::nearest_on_surface(p, self.surface, None);
                if es.is_finite() && et.is_finite() && ed.is_finite() {
                    s = es;
                    t = et;
                    let check = self.surface.point(s, t);
                    let mut dis_surf = p.distance_squared(check).sqrt();
                    let (mut uu, mut vv) = (s, t);
                    let mut dist_min_on_iso = REAL_LAST;
                    let mut possibly_local = false;
                    if dis_surf > preci {
                        let (res, pp) = self.surface_newton(Point2::new(uu, vv), p, preci);
                        if res != 0 {
                            let dist = p.distance(self.value(pp));
                            if dist < dis_surf {
                                dis_surf = dist;
                                s = pp.x;
                                uu = pp.x;
                                t = pp.y;
                                vv = pp.y;
                            }
                        }
                        if dis_surf < 10.0 * preci && self.surface.continuity() != Continuity::C0 {
                            const TOL: f64 = CONFUSION;
                            let (pnt, d1u, d1v) = self.surface.d1(uu, vv);
                            let b = d1u.cross(d1v);
                            let a = p - pnt;
                            let ab = a.dot(b);
                            let nrm2 = b.dot(b);
                            if nrm2 > 1e-10 {
                                let dist = a.dot(a) - (ab * ab) / nrm2;
                                possibly_local = dist < TOL * TOL;
                            }
                        }
                        if !possibly_local {
                            dist_min_on_iso = self.uv_from_iso(p, preci, &mut uu, &mut vv);
                        }
                    }
                    if dis_surf > dist_min_on_iso {
                        s = uu;
                        t = vv;
                        self.gap = dist_min_on_iso;
                    } else {
                        self.gap = dis_surf;
                    }
                } else {
                    let (mut uu, mut vv) = (s, t);
                    self.gap = self.uv_from_iso(p, preci, &mut uu, &mut vv);
                    s = uu;
                    t = vv;
                }
            }
        }
        if self.gap <= 0.0 {
            self.gap = p.distance(self.surface.point(s, t));
        }
        Point2::new(s, t)
    }

    /// `UVFromIso`: the parameters improved by projecting onto the bound isos and the isos
    /// through `(u, v)`; the distance reached.
    #[doc(alias = "ShapeAnalysis_Surface::UVFromIso")]
    fn uv_from_iso(&mut self, p: Point, preci: f64, u: &mut f64, v: &mut f64) -> f64 {
        let (mut uu, mut vv) = (*u, *v);
        let start = self.surface.point(*u, *v);
        let mut the_min = start.distance(p);
        if the_min < preci / 10.0 {
            return the_min;
        }
        let boxes = self.iso_boxes();
        // (OCCT's isos at infinite bounds lie infinitely far: none here, skipped below.)
        let isos = self.bound_isos().clone();
        let offset = kind(self.surface) == Kind::Offset;
        let mut uv_flag = true;
        let (mut par, mut other);
        for num in 0..6 {
            uv_flag = num < 3;
            let (iso, iso_box) = match num {
                0 => {
                    par = self.uf;
                    (isos[0].clone(), boxes[0].clone())
                }
                1 => {
                    par = self.ul;
                    (isos[1].clone(), boxes[1].clone())
                }
                2 => {
                    par = *u;
                    (self.u_iso(*u), None)
                }
                3 => {
                    par = self.vf;
                    (isos[2].clone(), boxes[2].clone())
                }
                4 => {
                    par = self.vl;
                    (isos[3].clone(), boxes[3].clone())
                }
                _ => {
                    par = *v;
                    (self.v_iso(*v), None)
                }
            };
            if !offset {
                let Some(iso) = iso else { continue };
                if is_infinite(par) {
                    continue;
                }
                if iso_box.as_ref().is_some_and(|b| box_distance(b, p) > the_min) {
                    continue;
                }
                let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
                restrict_bounds(&mut cf, &mut cl);
                let projected = project_on_curve(&iso, p, preci, cf, cl, false);
                other = projected.parameter;
                if projected.distance < the_min {
                    the_min = projected.distance;
                    uu = if uv_flag { par } else { other };
                    vv = if uv_flag { other } else { par };
                }
            } else {
                // (On an offset surface the isos are the adaptor's.)
                if iso_box.as_ref().is_some_and(|b| box_distance(b, p) > the_min) {
                    continue;
                }
                let Some(iso) = iso else { continue };
                let projected = project_on_adaptor(&iso, iso.first_parameter(), iso.last_parameter(), p, preci, true);
                other = projected.parameter;
                if projected.distance < the_min {
                    the_min = projected.distance;
                    uu = if uv_flag { par } else { other };
                    vv = if uv_flag { other } else { par };
                }
            }
        }
        // (Then alternately along the u- and v-isos through the point found.)
        let (mut prev_u, mut prev_v) = (*u, *v);
        let max_iters = 5;
        let mut iters = 0;
        while (prev_u != uu || prev_v != vv) && iters < max_iters && the_min > preci {
            prev_u = uu;
            prev_v = vv;
            for second in [false, true] {
                if second {
                    uv_flag = !uv_flag;
                }
                let iso = if uv_flag { self.u_iso(uu) } else { self.v_iso(vv) };
                if let Some(iso) = iso {
                    let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
                    restrict_bounds(&mut cf, &mut cl);
                    let projected = if !offset {
                        project_on_curve(&iso, p, preci, cf, cl, false)
                    } else if !second {
                        project_on_adaptor(&iso, iso.first_parameter(), iso.last_parameter(), p, preci, true)
                    } else {
                        project_act(&iso, iso.first_parameter(), iso.last_parameter(), p, preci)
                    };
                    if projected.distance < the_min {
                        the_min = projected.distance;
                        if uv_flag {
                            vv = projected.parameter;
                        } else {
                            uu = projected.parameter;
                        }
                    }
                }
            }
            uv_flag = !uv_flag;
            iters += 1;
        }
        *u = uu;
        *v = vv;
        the_min
    }
}

/// `Geom_Surface::UIso` (`VIso` for `!u_type`) as OCCT makes it: a sphere's meridian is the
/// half circle over its latitudes, a trimmed surface's iso is trimmed to its range (kentron's
/// own isos run over the whole curve).
fn geom_iso(surface: &Surface, u_type: bool, at: f64) -> Curve {
    // (`Geom_TrimmedCurve`: on a periodic basis the range moved into its first period.)
    let trimmed = |c: Curve, first: f64, last: f64| -> Curve {
        let (first, last) = match c.period() {
            Some(period) if c.is_periodic() => {
                let (lo, hi) = (c.first_parameter(), c.first_parameter() + period);
                super::adaptor::adjust_periodic(lo, hi, ((last - first).abs() / 2.0).min(P_CONFUSION), first, last)
            }
            _ => (first, last),
        };
        Curve::Trimmed(super::curve::TrimmedCurve::new(c, first, last))
    };
    match surface {
        Surface::Sphere(_) if u_type => trimmed(surface.u_iso(at), -PI / 2.0, PI / 2.0),
        Surface::Trimmed(t) => {
            let (u1, u2, v1, v2) = t.trim();
            let iso = geom_iso(t.basis(), u_type, at);
            if u_type { trimmed(iso, v1, v2) } else { trimmed(iso, u1, u2) }
        }
        _ => if u_type { surface.u_iso(at) } else { surface.v_iso(at) },
    }
}

/// `GeomAdaptor_Surface::UContinuity() == GeomAbs_C0` (`VContinuity` for `!along_u`) for a
/// B-spline over `[first, last]`: a knot inside of multiplicity reaching the degree.
fn u_continuity_is_c0(bs: &super::bspline::BSplineSurface, first: f64, last: f64, along_u: bool) -> bool {
    let (knots, multiplicities, degree, periodic) = if along_u {
        (bs.u_knots(), bs.u_multiplicities(), bs.u_degree(), bs.is_u_periodic())
    } else {
        (bs.v_knots(), bs.v_multiplicities(), bs.v_degree(), bs.is_v_periodic())
    };
    let n = knots.len();
    let interior = (0..n).filter(|&i| periodic || (i > 0 && i + 1 < n)).filter(|&i| knots[i] > first + P_CONFUSION && knots[i] < last - P_CONFUSION);
    interior.map(|i| multiplicities[i]).max().is_some_and(|m| m >= degree)
}
