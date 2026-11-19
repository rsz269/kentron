//! OCCT's `ShapeConstruct_ProjectCurveOnSurface`, ported: the pcurve of a 3D curve on a
//! surface as OCCT's shape healing (and so its STEP and IGES readers) makes it, the curve's
//! points projected one by one ([`super::surface_analysis`]) and interpolated, isolines and
//! straight runs recognised, seams, periods and singular points dealt with.
//!
//! `ProjLib_ProjectedCurve`, which OCCT uses on planes and on B-splines of very uneven
//! parametrization, is kentron's [`super::project`] here.

use crate::math::{Point, Point2, Vector2};

use super::curve2::BSplineCurve2;
use super::curve::Curve;
use super::curve2::{Curve2, Line2};
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;
use super::surface_analysis::{adjust_to_period, in_period, is_infinite, next_project_on_curve, project_on_curve, SurfaceAnalysis};

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = CONFUSION * 0.01;
/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = CONFUSION * CONFUSION;
/// `gp::Resolution()`.
const RESOLUTION: f64 = f64::MIN_POSITIVE;
/// `THE_NCONTROL`: the number of points a curve is sampled at, at least.
const N_CONTROL: usize = 23;

/// A point projected before, kept from one edge to the next (`CachePoint`).
type CachePoint = (Point, Point2);

/// OCCT's `ShapeConstruct_ProjectCurveOnSurface`: projects curves onto one surface, keeping
/// the ends of the last one projected (so that the next, sharing an end, starts in the same
/// period).
#[doc(alias = "ShapeConstruct_ProjectCurveOnSurface")]
pub(crate) struct PCurveProjector<'a> {
    surface: SurfaceAnalysis<'a>,
    preci: f64,
    /// `AdjustOverDegenMode`: -1 leaves runs over a singular point alone; 1 (the default)
    /// and 0 put them on the first (last) bound.
    pub adjust_over_degen: i32,
    cache: Vec<CachePoint>,
    /// `ShapeExtend_DONE4`: a run over a singular point adjusted.
    pub adjusted_over_degen: bool,
}

impl<'a> PCurveProjector<'a> {
    /// `Init(surface, preci)`.
    pub(crate) fn new(surface: &'a Surface, preci: f64) -> Self {
        PCurveProjector { surface: SurfaceAnalysis::new(surface), preci, adjust_over_degen: 1, cache: Vec::new(), adjusted_over_degen: false }
    }

    /// `SetPrecision`.
    pub(crate) fn set_precision(&mut self, preci: f64) {
        self.preci = preci;
    }

    /// The surface's analysis (shared with what uses the projector).
    pub(crate) fn analysis(&mut self) -> &mut SurfaceAnalysis<'a> {
        &mut self.surface
    }

    /// `Perform`: the pcurve of `curve` over `[first, last]` (parametrized as the curve),
    /// `tol_first` and `tol_last` the tolerances of its end vertices (negative: none).
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::Perform")]
    pub(crate) fn perform(&mut self, curve: &Curve, first: f64, last: f64, tol_first: f64, tol_last: f64) -> Option<Curve2> {
        self.adjusted_over_degen = false;
        if let Some(c2d) = self.project_analytic(curve, first, last) {
            return Some(c2d);
        }
        if is_bspline_curve_invalid(curve, first, last) {
            if let Some(c2d) = self.perform_by_projlib(curve, first, last) {
                return Some(c2d);
            }
        }
        let (mut points, mut params) = generate_curve_points(curve, first, last, N_CONTROL);
        let n = points.len();
        let mut points2d = vec![Point2::new(0.0, 0.0); n];
        if let Some(c2d) = self.approx_pcurve(n, curve, tol_first, tol_last, &mut points, &mut params, &mut points2d) {
            return Some(c2d);
        }
        self.interpolate_pcurve(&points2d, &params)
    }

    /// `PerformByProjLib`.
    fn perform_by_projlib(&self, curve: &Curve, first: f64, last: f64) -> Option<Curve2> {
        super::project::pcurve_on_range(curve, self.surface.surface(), first, last, CONFUSION).map(Curve2::BSpline)
    }

    /// `projectAnalytic`: on a plane (or a trimmed or offset plane), the curve projected
    /// along the normal.
    fn project_analytic(&self, curve: &Curve, first: f64, last: f64) -> Option<Curve2> {
        let plane = match self.surface.surface() {
            Surface::Plane(p) => *p,
            Surface::Trimmed(t) => match t.basis() {
                Surface::Plane(p) => *p,
                _ => return None,
            },
            Surface::Offset(o) => match o.basis() {
                Surface::Plane(p) => *p,
                _ => return None,
            },
            _ => return None,
        };
        let _ = (first, last);
        super::project::analytic_pcurve(curve, &Surface::Plane(plane))
    }

    /// `getLine`: the pcurve as a straight segment where the curve's points all lie on it
    /// (the ends projected, the ends' parameters filled in either way).
    fn get_line(&mut self, points: &[Point], params: &[f64], points2d: &mut [Point2], tol: f64, recompute: &mut bool, from_cache: &mut bool) -> Option<Curve2> {
        let n = points.len();
        let p = [points[0], points[1], points[n - 2], points[n - 1]];
        let mut p2d = [Point2::new(0.0, 0.0); 4];
        let mut tol2 = tol * tol;
        let mut tol_working = tol;
        let surface = self.surface.surface();
        let periodic_u = surface.is_u_periodic();
        let periodic_v = surface.is_v_periodic();
        if tol2 > 1.0 {
            tol_working = CONFUSION;
            tol2 = tol_working * tol_working;
        }
        if tol2 < SQUARE_CONFUSION {
            tol2 = SQUARE_CONFUSION;
        }
        let old_tol2 = tol2;
        let mut saved_index: i32 = -1;
        let mut saved_point = Point2::new(0.0, 0.0);
        let corners = corner_cache(surface);
        let cache = self.cache.clone();
        // (OCCT's lambda sees the tolerance as it grows: passed in.)
        let mut project = |sas: &mut SurfaceAnalysis<'a>, point: Point, index: usize, tol2: f64| -> Point2 {
            for c in &cache {
                if c.0.distance_squared(point) < tol2 {
                    saved_index = index as i32;
                    saved_point = c.1;
                    if index == 0 {
                        *from_cache = true;
                    }
                    return projector_next_value(sas, &corners, c.1, point, tol_working, tol2, tol_working);
                }
            }
            projector_value(sas, &corners, point, tol_working, tol2)
        };
        for i in [0, 3] {
            p2d[i] = project(&mut self.surface, p[i], i, tol2);
            let d = self.surface.gap();
            if tol2 < d * d {
                tol2 = d * d;
            }
        }
        if periodic_u || periodic_v {
            for i in [1, 2] {
                p2d[i] = project(&mut self.surface, p[i], i, tol2);
                let d = self.surface.gap();
                if tol2 < d * d {
                    tol2 = d * d;
                }
            }
            if periodic_u {
                *recompute = fix_periodicity_troubles(&mut p2d, 0, surface.u_period().expect("periodic"), saved_index, saved_point.x);
            }
            if periodic_v {
                *recompute = fix_periodicity_troubles(&mut p2d, 1, surface.v_period().expect("periodic"), saved_index, saved_point.y);
            }
        }
        if *recompute && matches!(surface, Surface::Sphere(_)) {
            return None;
        }
        points2d[0] = p2d[0];
        points2d[n - 1] = p2d[3];
        tol2 = old_tol2;
        let d_par = params[n - 1] - params[0];
        if d_par.abs() < P_CONFUSION {
            return None;
        }
        let vec0 = p2d[3] - p2d[0];
        let vec = vec0 / d_par;
        let mut normal_check = surface.is_cn_u(1) && surface.is_cn_v(1);
        if normal_check {
            for i in 0..n {
                let cur = p2d[0] + vec * (params[i] - params[0]);
                let (q, du, dv) = surface.d1(cur.x, cur.y);
                let normal = du.cross(dv);
                if normal.dot(normal) < SQUARE_CONFUSION {
                    normal_check = false;
                    break;
                }
                // (The distance of the point from the normal line there.)
                let nd = normal / normal.length();
                let w = points[i] - q;
                let dist = (w - nd * w.dot(nd)).length();
                if dist > tol_working {
                    return None;
                }
            }
        }
        if !normal_check {
            let first_dist = surface.point(p2d[0].x, p2d[0].y).distance_squared(points[0]);
            tol2 = tol2.max(tol2 * 2.0 * first_dist);
            for i in 1..n - 1 {
                let cur = p2d[0] + vec * (params[i] - params[0]);
                let d1 = surface.point(cur.x, cur.y).distance_squared(points[i]);
                if (first_dist - d1).abs() > tol2 {
                    return None;
                }
            }
        }
        let length = vec0.length();
        if (length - d_par).abs() <= P_CONFUSION {
            let dir = vec0 / length;
            let origin = p2d[0] + dir * (-params[0]);
            return Some(Curve2::Line(Line2::new(origin, crate::math::Direction2::new(dir.x, dir.y))));
        }
        Some(Curve2::BSpline(BSplineCurve2::new(vec![p2d[0], p2d[3]], vec![params[0], params[n - 1]], vec![2, 2], 1)))
    }

    /// The ends kept for the next curve (`myCache`), in the order OCCT keeps them.
    fn fill_cache(&mut self, points: &[Point], points2d: &[Point2], change_cycle: bool) {
        let n = points.len();
        self.cache = if change_cycle {
            vec![(points[0], points2d[0]), (points[n - 1], points2d[n - 1])]
        } else {
            vec![(points[n - 1], points2d[n - 1]), (points[0], points2d[0])]
        };
    }

    /// `approxPCurve`: the curve's points projected (`points2d`), with periods, seams and
    /// singular points dealt with; a straight pcurve where `getLine` finds one.
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::approxPCurve")]
    fn approx_pcurve(&mut self, nb_pnt: usize, c3d: &Curve, tol_first: f64, tol_last: f64, points: &mut Vec<Point>, params: &mut Vec<f64>, points2d: &mut Vec<Point2>) -> Option<Curve2> {
        let mut recompute = false;
        let mut from_cache_line = false;
        let preci = self.preci;
        if let Some(line) = self.get_line(points, params, points2d, preci, &mut recompute, &mut from_cache_line) {
            let change_cycle = !self.cache.is_empty()
                && self.cache[0].0.distance(points[0]) > self.cache[0].0.distance(points[nb_pnt - 1])
                && self.cache[0].0.distance(points[nb_pnt - 1]) < CONFUSION;
            self.fill_cache(points, points2d, change_cycle);
            return Some(line);
        }
        let surface = self.surface.surface();
        let analytic = !matches!(surface, Surface::Bezier(_) | Surface::BSpline(_));
        let (uf, ul, vf, vl) = surface.bounds();
        let iso_closed;
        let mut pout = vec![0.0; nb_pnt];
        let iso = self.is_an_isoparametric(nb_pnt, points, params, &mut pout);
        let (mut iso_value, mut iso_par1, mut iso_par2) = (0.0, 0.0, 0.0);
        let mut t_par = 0.0;
        let mut value_p1 = iso.value_p1;
        let mut value_p2 = iso.value_p2;
        if iso.iso_param {
            let (parf, parl);
            if iso.type_u {
                iso_value = value_p1.x;
                iso_par1 = value_p1.y;
                iso_par2 = value_p2.y;
                iso_closed = self.surface.is_v_closed(preci);
                parf = vf;
                parl = vl;
            } else {
                iso_value = value_p1.y;
                iso_par1 = value_p1.x;
                iso_par2 = value_p2.x;
                iso_closed = self.surface.is_u_closed(preci);
                parf = uf;
                parl = ul;
            }
            if !iso.par_2d_3d && !analytic {
                let tdeb = pout[1];
                if iso_closed && (iso_par1 == parf || iso_par1 == parl) {
                    iso_par1 = if (tdeb - parf).abs() < (tdeb - parl).abs() { parf } else { parl };
                    if iso.type_u {
                        value_p1.y = iso_par1;
                    } else {
                        value_p1.x = iso_par1;
                    }
                }
                if iso_closed && (iso_par2 == parf || iso_par2 == parl) {
                    let tfin = pout[nb_pnt - 2];
                    iso_par2 = if (tfin - parf).abs() < (tfin - parl).abs() { parf } else { parl };
                    if iso.type_u {
                        value_p2.y = iso_par2;
                    } else {
                        value_p2.x = iso_par2;
                    }
                }
                if !iso_closed && (tdeb - iso_par1).abs() > (tdeb - iso_par2).abs() && (pout[nb_pnt - 2] - iso_par2).abs() > (pout[nb_pnt - 2] - iso_par1).abs() {
                    std::mem::swap(&mut value_p1, &mut value_p2);
                    if iso.type_u {
                        iso_value = value_p1.x;
                        iso_par1 = value_p1.y;
                        iso_par2 = value_p2.y;
                    } else {
                        iso_value = value_p1.y;
                        iso_par1 = value_p1.x;
                        iso_par2 = value_p2.x;
                    }
                }
            }
        }
        let up = ul - uf;
        let vp = vl - vf;
        let mut gap = preci;
        let mut change_cycle = false;
        let mut from_cache = false;
        let mut saved_point = Point2::new(0.0, 0.0);
        if !self.cache.is_empty() && self.cache[0].0.distance(points[0]) > self.cache[0].0.distance(points[nb_pnt - 1]) && self.cache[0].0.distance(points[nb_pnt - 1]) < CONFUSION {
            change_cycle = true;
        }
        let mut resolve_u_jump = false;
        let mut resolve_v_jump = false;
        let mut prev_p3d = Point::ORIGIN;
        let mut prev_p2d = Point2::new(0.0, 0.0);
        let mut p2d = Point2::new(0.0, 0.0);
        for ii in 1..=nb_pnt {
            let index = if change_cycle { nb_pnt - ii } else { ii - 1 };
            let p3d = points[index];
            if iso.iso_param {
                if iso.par_2d_3d {
                    t_par = if iso_par2 > iso_par1 { params[index] } else { iso.t1 + iso.t2 - params[index] };
                } else if !analytic {
                    t_par = if index == 0 {
                        iso_par1
                    } else if index == nb_pnt - 1 {
                        iso_par2
                    } else {
                        pout[index]
                    };
                }
                if !iso.par_2d_3d && analytic {
                    if index == 0 {
                        p2d = value_p1;
                    } else if index == nb_pnt - 1 {
                        p2d = value_p2;
                    } else {
                        p2d = self.surface.next_value_of_uv(p2d, p3d, preci, CONFUSION + 1000.0 * gap);
                        gap = self.surface.gap();
                    }
                } else if iso.type_u {
                    p2d = Point2::new(iso_value, t_par);
                } else {
                    p2d = Point2::new(t_par, iso_value);
                }
            } else if index == 0 && iso.p1_on_iso {
                p2d = value_p1;
            } else if index == nb_pnt - 1 && iso.p2_on_iso {
                p2d = value_p2;
            } else {
                if index == 0 || index == nb_pnt - 1 {
                    if !recompute {
                        p2d = points2d[index];
                        gap = self.surface.gap();
                        if index == 0 {
                            from_cache = from_cache_line;
                            saved_point = p2d;
                        }
                        // (OCCT's `continue`: no jump check, no extrapolation for this one.)
                        continue;
                    }
                    let mut found = false;
                    for c in self.cache.clone() {
                        if c.0.distance_squared(p3d) < preci * preci {
                            p2d = self.surface.next_value_of_uv(c.1, p3d, preci, CONFUSION + gap);
                            if index == 0 {
                                from_cache = true;
                                saved_point = c.1;
                            }
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        p2d = self.surface.value_of_uv(p3d, preci);
                    }
                } else {
                    p2d = self.surface.next_value_of_uv(p2d, p3d, preci, CONFUSION + 1000.0 * gap);
                }
                gap = self.surface.gap();
            }
            points2d[index] = p2d;
            if nb_pnt > 23 && ii > 2 && ii < nb_pnt {
                let bspline = matches!(surface, Surface::BSpline(_));
                if (p2d.x - prev_p2d.x).abs() > 0.95 * up && prev_p3d.distance(p3d) < preci && !self.surface.is_u_closed(preci) && self.surface.nb_singularities(preci) > 0 && bspline {
                    resolve_u_jump = true;
                }
                if (p2d.y - prev_p2d.y).abs() > 0.95 * vp && prev_p3d.distance(p3d) < preci && !self.surface.is_v_closed(preci) && self.surface.nb_singularities(preci) > 0 && bspline {
                    resolve_v_jump = true;
                }
            }
            prev_p3d = p3d;
            prev_p2d = p2d;
            if ii > 1 {
                // (The next search starts from the point extrapolated.)
                let other = if change_cycle { points2d[index + 1] } else { points2d[index - 1] };
                p2d = Point2::new(2.0 * p2d.x - other.x, 2.0 * p2d.y - other.y);
            }
        }
        if !iso.par_2d_3d {
            self.surface.project_degenerated(points, points2d, preci, true);
            self.surface.project_degenerated(points, points2d, preci, false);
        }
        // Ends at singular points.
        let point_first = points[0];
        let point_last = points[points.len() - 1];
        let tol_first = if tol_first < 0.0 { CONFUSION } else { tol_first };
        let tol_last = if tol_last < 0.0 { CONFUSION } else { tol_last };
        let mut i = 0;
        while let Some(s) = self.surface.singularity(i) {
            if s.preci <= CONFUSION && point_first.distance(s.point) <= tol_first {
                self.correct_extremity(c3d, params, points2d, true, s.first_uv, s.u_iso);
            }
            if s.preci <= CONFUSION && point_last.distance(s.point) <= tol_last {
                self.correct_extremity(c3d, params, points2d, false, s.first_uv, s.u_iso);
            }
            i += 1;
        }
        let tol_on_u_period = CONFUSION * up;
        let tol_on_v_period = CONFUSION * vp;
        // U-closed surfaces.
        if self.surface.is_u_closed(preci) || resolve_u_jump {
            self.follow_period(0, uf, ul, up, tol_on_u_period, from_cache, saved_point.x, c3d, points, params, points2d);
        }
        // V-closed surfaces.
        let sphere = matches!(surface, Surface::Sphere(_));
        if self.surface.is_v_closed(preci) || resolve_v_jump || sphere {
            self.follow_period(1, vf, vl, vp, tol_on_v_period, from_cache, saved_point.y, c3d, points, params, points2d);
        }
        // V-closed: a jump across the seam put on the bound nearer.
        if self.surface.is_v_closed(preci) || sphere {
            for k in 1..points2d.len() {
                let dist2d = (points2d[k].y - points2d[k - 1].y).abs();
                if dist2d > vp / 2.0 {
                    // (The nearest of the four distances, the first on a tie, as OCCT's chain
                    // of comparisons picks it.)
                    let distances = [
                        (points2d[k - 1].y - vf).abs(),
                        (points2d[k - 1].y - vl).abs(),
                        (points2d[k].y - vf).abs(),
                        (points2d[k].y - vl).abs(),
                    ];
                    let mut nearest = 0;
                    for (i, &d) in distances.iter().enumerate().skip(1) {
                        if d < distances[nearest] {
                            nearest = i;
                        }
                    }
                    match nearest {
                        0 => points2d[k - 1].y = vf,
                        1 => points2d[k - 1].y = vl,
                        2 => points2d[k].y = vf,
                        _ => points2d[k].y = vl,
                    }
                }
            }
        }
        // Runs over a singular point.
        if self.adjust_over_degen != -1 {
            if self.surface.is_u_closed(preci) {
                self.surface.is_degenerated(Point::ORIGIN, preci);
                if self.surface.nb_singularities(preci) > 0 {
                    self.adjust_over_degenerated(0, uf, ul, up, points, points2d);
                }
            } else if self.surface.is_v_closed(preci) {
                self.surface.is_degenerated(Point::ORIGIN, preci);
                if self.surface.nb_singularities(preci) > 0 {
                    self.adjust_over_degenerated(1, vf, vl, vp, points, points2d);
                }
            }
        }
        self.fill_cache(points, points2d, change_cycle);
        None
    }

    /// The `Handle U-closed (V-closed) surfaces` part of `approxPCurve`: the first point
    /// brought into the period, each next one within half a period of the one before (or a
    /// point added between them where the curve crosses the seam), all moved by a period
    /// where their middle lies outside.
    #[allow(clippy::too_many_arguments)]
    fn follow_period(&mut self, coord: usize, first: f64, last: f64, period: f64, tol_on_period: f64, from_cache: bool, saved: f64, c3d: &Curve, points: &mut Vec<Point>, params: &mut Vec<f64>, points2d: &mut Vec<Point2>) {
        let get = |p: &Point2| if coord == 0 { p.x } else { p.y };
        let set = |p: &mut Point2, x: f64| if coord == 0 { p.x = x } else { p.y = x };
        let mut first_x = get(&points2d[0]);
        if !from_cache {
            while first_x < first {
                first_x += period;
                set(&mut points2d[0], first_x);
            }
            while first_x > last {
                first_x -= period;
                set(&mut points2d[0], first_x);
            }
        }
        let periodic = if coord == 0 { self.surface.surface().is_u_periodic() } else { self.surface.surface().is_v_periodic() };
        if periodic && from_cache {
            let (mut min, mut max) = (first, last);
            while min > saved {
                min -= period;
                max -= period;
            }
            while max < saved {
                min += period;
                max += period;
            }
            first_x += adjust_to_period(first_x, min, max);
            set(&mut points2d[0], first_x);
        }
        let mut prev_x = first_x;
        let (mut min_x, mut max_x) = (first_x, first_x);
        let mut to_adjust = false;
        let mut k = 1;
        while k < points2d.len() {
            let mut cur_x = get(&points2d[k]);
            let dist2d = (cur_x - prev_x).abs();
            if dist2d > period / 2.0 {
                self.insert_additional_point_or_adjust(&mut to_adjust, coord, period, tol_on_period, &mut cur_x, prev_x, c3d, &mut k, points, params, points2d);
            }
            prev_x = cur_x;
            if min_x > cur_x {
                min_x = cur_x;
            } else if max_x < cur_x {
                max_x = cur_x;
            }
            k += 1;
        }
        if !from_cache {
            let mid = 0.5 * (min_x + max_x);
            let shift = if mid > last { -period } else if mid < first { period } else { 0.0 };
            if shift != 0.0 {
                for p in points2d.iter_mut() {
                    let x = get(p) + shift;
                    set(p, x);
                }
            }
        }
    }

    /// The `Handle AdjustOverDegen` part of `approxPCurve` for coordinate `coord` (0: u).
    fn adjust_over_degenerated(&mut self, coord: usize, first: f64, last: f64, period: f64, points: &[Point], points2d: &mut [Point2]) {
        let get = |p: &Point2| if coord == 0 { p.x } else { p.y };
        let set = |p: &mut Point2, x: f64| if coord == 0 { p.x = x } else { p.y = x };
        let mut prev_x = 0.0;
        let (mut on_bound, mut prev_on_bound) = (false, false);
        let mut start = true;
        let mut ind = 0;
        while ind < points2d.len() {
            let cur_x = get(&points2d[ind]);
            if self.surface.is_degenerated(points[ind], CONFUSION) {
                ind += 1;
                continue;
            }
            on_bound = ((cur_x - 0.5 * (last + first)).abs() - period / 2.0).abs() <= P_CONFUSION;
            if !start && ((cur_x - prev_x).abs() - period / 2.0).abs() <= 0.01 * period {
                break;
            }
            start = false;
            prev_x = cur_x;
            prev_on_bound = on_bound;
            ind += 1;
        }
        if ind < points2d.len() {
            let prev = if self.adjust_over_degen != 0 { first } else { last };
            let d = period / 2.0 + P_CONFUSION;
            let bring = |p: &mut Point2| {
                let mut x = get(p);
                while x < prev - d {
                    x += period;
                    set(p, x);
                }
                while x > prev + d {
                    x -= period;
                    set(p, x);
                }
            };
            if prev_on_bound {
                set(&mut points2d[ind - 1], prev);
                for j in (0..ind.saturating_sub(1)).rev() {
                    bring(&mut points2d[j]);
                }
            } else if on_bound {
                set(&mut points2d[ind], prev);
                for p in points2d.iter_mut().skip(ind + 1) {
                    bring(p);
                }
            }
            self.adjusted_over_degen = true;
        }
    }

    /// `correctExtremity`: an end at a singular point moved along its iso to where the
    /// pcurve's last segments, extended, meet it.
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::correctExtremity")]
    fn correct_extremity(&mut self, c3d: &Curve, params: &[f64], points2d: &mut [Point2], is_first: bool, point_on_iso: Point2, u_iso: bool) {
        let nb = points2d.len();
        let ind = if u_iso { 1 } else { 0 };
        let other = 1 - ind;
        let coord = |p: Point2, i: usize| if i == 0 { p.x } else { p.y };
        let singularity_coord = coord(point_on_iso, other);
        let end_point = if is_first { points2d[0] } else { points2d[nb - 1] };
        let finish_coord = coord(end_point, other);
        let iso_dir = if u_iso { Vector2::new(0.0, 1.0) } else { Vector2::new(1.0, 0.0) };
        let periodic = if u_iso { self.surface.surface().is_v_periodic() } else { self.surface.surface().is_u_periodic() };
        let (mut first_of_line, mut second_of_line, finish_param, mut second_param) = if is_first {
            (points2d[2], points2d[1], params[0], params[1])
        } else {
            (points2d[nb - 3], points2d[nb - 2], params[nb - 1], params[nb - 2])
        };
        if singularity_coord > finish_coord && coord(second_of_line, other) > finish_coord {
            return;
        }
        if singularity_coord < finish_coord && coord(second_of_line, other) < finish_coord {
            return;
        }
        {
            let prev_dist = (coord(second_of_line, ind) - coord(first_of_line, ind)).abs();
            let cur_dist = (coord(end_point, ind) - coord(second_of_line, ind)).abs();
            if cur_dist <= 2.0 * prev_dist {
                return;
            }
        }
        let on_finish = |p: Point2| if u_iso { Point2::new(finish_coord, p.y) } else { Point2::new(p.x, finish_coord) };
        let mut finish_point = on_finish(second_of_line);
        loop {
            if (coord(second_of_line, other) - finish_coord).abs() <= 2.0 * P_CONFUSION {
                break;
            }
            let v = second_of_line - first_of_line;
            if v.dot(v) <= 1.0e-32 {
                break;
            }
            // (IntCurve_IntConicConic of the iso through the end and the line through the last
            // two points.)
            finish_point = match intersect_lines(end_point, iso_dir, first_of_line, v) {
                Some(p) => p,
                None => on_finish(second_of_line),
            };
            let prev_point = first_of_line;
            first_of_line = second_of_line;
            let first_param = second_param;
            second_param = (first_param + finish_param) / 2.0;
            if (second_param - first_param).abs() <= 2.0 * P_CONFUSION {
                break;
            }
            let p3d = c3d.point(second_param);
            second_of_line = self.surface.next_value_of_uv(first_of_line, p3d, self.preci, CONFUSION);
            if periodic {
                adjust_second_to_first_point(first_of_line, &mut second_of_line, self.surface.surface());
            }
            let prev_dist = (coord(first_of_line, ind) - coord(prev_point, ind)).abs();
            let cur_dist = (coord(second_of_line, ind) - coord(first_of_line, ind)).abs();
            if cur_dist > 2.0 * prev_dist {
                break;
            }
        }
        if is_first {
            points2d[0] = finish_point;
        } else {
            points2d[nb - 1] = finish_point;
        }
    }

    /// `insertAdditionalPointOrAdjust`: where two points lie more than half a period apart,
    /// either the second moved by a period to follow the first, or a point inserted between
    /// them where the curve crosses the seam.
    #[allow(clippy::too_many_arguments)]
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::insertAdditionalPointOrAdjust")]
    fn insert_additional_point_or_adjust(&mut self, to_adjust: &mut bool, coord: usize, period: f64, tol_on_period: f64, cur: &mut f64, prev: f64, c3d: &Curve, index: &mut usize, points: &mut Vec<Point>, params: &mut Vec<f64>, points2d: &mut Vec<Point2>) {
        let get = |p: Point2| if coord == 0 { p.x } else { p.y };
        let corrected = in_period(*cur, prev - period / 2.0, prev + period / 2.0);
        if !*to_adjust {
            let cur_par = params[*index];
            let prev_par = params[*index - 1];
            let mut mid_par = (prev_par + cur_par) / 2.0;
            let mut mid_p3d = c3d.point(mid_par);
            let mut mid_p2d = self.surface.value_of_uv(mid_p3d, self.preci);
            let mut mid = in_period(get(mid_p2d), prev - period / 2.0, prev + period / 2.0);
            let (mut first_coord, mut last_coord) = (prev, corrected);
            if last_coord < first_coord {
                std::mem::swap(&mut first_coord, &mut last_coord);
            }
            if last_coord - first_coord <= tol_on_period {
                *to_adjust = true;
            } else if first_coord <= mid && mid <= last_coord {
                *to_adjust = true;
            } else {
                let mut success = true;
                let (mut first_t, mut last_t) = (prev_par, cur_par);
                mid = get(mid_p2d);
                while (mid - prev).abs() >= period / 2.0 - tol_on_period || (*cur - mid).abs() >= period / 2.0 - tol_on_period {
                    if mid_par - first_t <= P_CONFUSION || last_t - mid_par <= P_CONFUSION {
                        success = false;
                        break;
                    }
                    if (mid - prev).abs() >= period / 2.0 - tol_on_period {
                        last_t = (first_t + last_t) / 2.0;
                    } else {
                        first_t = (first_t + last_t) / 2.0;
                    }
                    mid_par = (first_t + last_t) / 2.0;
                    mid_p3d = c3d.point(mid_par);
                    mid_p2d = self.surface.value_of_uv(mid_p3d, self.preci);
                    mid = get(mid_p2d);
                }
                if success {
                    points.insert(*index, mid_p3d);
                    params.insert(*index, mid_par);
                    points2d.insert(*index, mid_p2d);
                    *index += 1;
                } else {
                    *to_adjust = true;
                }
            }
        }
        if *to_adjust {
            *cur = corrected;
            if coord == 0 {
                points2d[*index].x = corrected;
            } else {
                points2d[*index].y = corrected;
            }
        }
    }

    /// `interpolatePCurve`: the cubic through the points at the curve's parameters
    /// (`Geom2dAPI_Interpolate`), points too close together dropped.
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::interpolatePCurve")]
    fn interpolate_pcurve(&self, points2d: &[Point2], params: &[f64]) -> Option<Curve2> {
        let n = points2d.len();
        let tolerance2d = self.preci / (100.0 * n as f64);
        let mut points = points2d.to_vec();
        let mut params = params.to_vec();
        let mut preci = tolerance2d;
        check_points_2d(&mut points, &mut params, &mut preci);
        if points.len() < 2 {
            return None;
        }
        super::fit::interpolate_with_tangents(&points, Some(&params), &[], false, false, preci).ok().map(Curve2::BSpline)
    }

    /// `isAnIsoparametric`: whether the curve runs along a bound iso of the surface.
    #[doc(alias = "ShapeConstruct_ProjectCurveOnSurface::isAnIsoparametric")]
    fn is_an_isoparametric(&mut self, nb_pnt: usize, points: &[Point], params: &[f64], params_out: &mut [f64]) -> IsoParametric {
        let prec = CONFUSION;
        let mut out = IsoParametric::default();
        let surface = self.surface.surface();
        let (u1, u2, v1, v2) = match surface {
            Surface::Trimmed(t) => {
                let (a, b, c, d) = t.trim();
                (a, b, c, d)
            }
            _ => self.surface.bounds(),
        };
        let mut mpt = [0u8; 2];
        let mut tpar = [0.0; 2];
        let mut iso_value = 0.0;
        let mut mindist2 = 4.0 * prec * prec;
        let mut mind2 = [mindist2; 2];
        let boxes = self.surface.iso_boxes();
        let mut c_iso: Option<Curve> = None;
        for j in 1..=4 {
            let (iso_val, iso_u, iso_box) = match j {
                1 => (u1, true, &boxes[0]),
                2 => (u2, true, &boxes[1]),
                3 => (v1, false, &boxes[2]),
                _ => (v2, false, &boxes[3]),
            };
            if is_infinite(iso_val) {
                continue;
            }
            let ci = if iso_u { self.surface.u_iso(iso_val) } else { self.surface.v_iso(iso_val) };
            let Some(ci) = ci else { continue };
            let (tt1, tt2) = if iso_u { (v1, v2) } else { (u1, u2) };
            let ext1 = ci.point(tt1);
            let ext2 = ci.point(tt2);
            let extmi = ci.point((tt1 + tt2) / 2.0);
            if ext1.distance(ext2) <= prec && ext1.distance(extmi) <= prec {
                continue;
            }
            let (mut eq1, mut eq2) = (false, false);
            let mut currd2 = [0.0; 2];
            let mut tp = [0.0; 2];
            let mut mp = [0u8; 2];
            for i in 0..2 {
                let k = if i == 0 { 0 } else { nb_pnt - 1 };
                currd2[i] = points[k].distance_squared(ext1);
                if currd2[i] <= prec * prec && !eq1 {
                    mp[i] = 1;
                    tp[i] = tt1;
                    eq1 = true;
                    continue;
                }
                currd2[i] = points[k].distance_squared(ext2);
                if currd2[i] <= prec * prec && !eq2 {
                    mp[i] = 2;
                    tp[i] = tt2;
                    eq2 = true;
                    continue;
                }
                if matches!(surface, Surface::Sphere(_)) && !iso_u {
                    continue;
                }
                // (`Bnd_Box::IsOut` of a box of OCCT's iso, which an infinite bound's never has.)
                if iso_box.as_ref().is_some_and(|b| !b.contains(points[k])) {
                    continue;
                }
                let (mut cf, mut cl) = (ci.first_parameter(), ci.last_parameter());
                if is_infinite(cf) {
                    cf = -1000.0;
                }
                if is_infinite(cl) {
                    cl = 1000.0;
                }
                let projected = project_on_curve(&ci, points[k], prec, cf, cl, true);
                currd2[i] = projected.distance * projected.distance;
                if projected.distance <= prec && projected.parameter >= cf && projected.parameter <= cl {
                    mp[i] = 3;
                    tp[i] = projected.parameter;
                }
            }
            if mp[0] > 0 && mp[1] > 0 && (tp[0] - tp[1]).abs() < P_CONFUSION {
                continue;
            }
            if mp[0] > 0 && (!out.p1_on_iso || currd2[0] < mind2[0]) {
                out.p1_on_iso = true;
                mind2[0] = currd2[0];
                out.value_p1 = if iso_u { Point2::new(iso_val, tp[0]) } else { Point2::new(tp[0], iso_val) };
            }
            if mp[1] > 0 && (!out.p2_on_iso || currd2[1] < mind2[1]) {
                out.p2_on_iso = true;
                mind2[1] = currd2[1];
                out.value_p2 = if iso_u { Point2::new(iso_val, tp[1]) } else { Point2::new(tp[1], iso_val) };
            }
            if mp[0] == 0 || mp[1] == 0 {
                continue;
            }
            let md2 = currd2[0] + currd2[1];
            if mindist2 <= md2 {
                continue;
            }
            mindist2 = md2;
            mpt = mp;
            tpar = tp;
            out.type_u = iso_u;
            iso_value = iso_val;
            c_iso = Some(ci);
            out.t1 = tt1;
            out.t2 = tt2;
        }
        if mpt[0] > 0 && mpt[1] > 0 {
            let c_iso = c_iso.expect("an iso found");
            out.p1_on_iso = true;
            out.p2_on_iso = true;
            if out.type_u {
                out.value_p1 = Point2::new(iso_value, tpar[0]);
                out.value_p2 = Point2::new(iso_value, tpar[1]);
            } else {
                out.value_p1 = Point2::new(tpar[0], iso_value);
                out.value_p2 = Point2::new(tpar[1], iso_value);
            }
            if mpt[0] != 3 && mpt[1] != 3 {
                out.par_2d_3d = true;
                for i in 1..nb_pnt - 1 {
                    if !out.par_2d_3d {
                        break;
                    }
                    let t = if tpar[1] > tpar[0] { params[i] } else { out.t1 + out.t2 - params[i] };
                    if points[i].distance(c_iso.point(t)) > prec {
                        out.par_2d_3d = false;
                    }
                }
            }
            if out.par_2d_3d {
                out.iso_param = true;
            } else {
                let mut prev_param = tpar[0];
                let mut by_distance = true;
                let (mut cf, mut cl) = (c_iso.first_parameter(), c_iso.last_parameter());
                if is_infinite(cf) {
                    cf = -1000.0;
                }
                if is_infinite(cl) {
                    cl = 1000.0;
                }
                for i in 1..nb_pnt - 1 {
                    if !by_distance {
                        break;
                    }
                    let projected = next_project_on_curve(prev_param, &c_iso, points[i], prec, cf, cl, false);
                    prev_param = projected.parameter;
                    params_out[i] = projected.parameter;
                    if projected.distance > prec || projected.parameter < cf || projected.parameter > cl {
                        by_distance = false;
                    }
                }
                if by_distance {
                    out.iso_param = true;
                }
            }
        }
        out
    }
}

/// What `isAnIsoparametric` finds.
#[derive(Debug, Clone, Copy, Default)]
struct IsoParametric {
    iso_param: bool,
    type_u: bool,
    p1_on_iso: bool,
    value_p1: Point2,
    p2_on_iso: bool,
    value_p2: Point2,
    par_2d_3d: bool,
    t1: f64,
    t2: f64,
}

/// The corners a clamped B-spline surface passes through, with their parameters
/// (`SurfaceProjectorWithCache`'s corner cache).
fn corner_cache(surface: &Surface) -> Vec<CachePoint> {
    let Surface::BSpline(bs) = surface else { return Vec::new() };
    let (nu, nv) = (bs.poles().len(), bs.poles()[0].len());
    let (du, dv) = (bs.u_degree(), bs.v_degree());
    let (um, vm) = (bs.u_multiplicities(), bs.v_multiplicities());
    let u_first_clamped = um[0] >= du + 1;
    let u_last_clamped = um[um.len() - 1] >= du + 1;
    let v_first_clamped = vm[0] >= dv + 1;
    let v_last_clamped = vm[vm.len() - 1] >= dv + 1;
    let (uk, vk) = (bs.u_knots(), bs.v_knots());
    let (u_first, u_last, v_first, v_last) = (uk[0], uk[uk.len() - 1], vk[0], vk[vk.len() - 1]);
    let mut out = Vec::new();
    if u_first_clamped && v_first_clamped {
        out.push((bs.pole(0, 0), Point2::new(u_first, v_first)));
    }
    if u_last_clamped && v_first_clamped {
        out.push((bs.pole(nu - 1, 0), Point2::new(u_last, v_first)));
    }
    if u_first_clamped && v_last_clamped {
        out.push((bs.pole(0, nv - 1), Point2::new(u_first, v_last)));
    }
    if u_last_clamped && v_last_clamped {
        out.push((bs.pole(nu - 1, nv - 1), Point2::new(u_last, v_last)));
    }
    out
}

/// `SurfaceProjectorWithCache::ValueOfUV`.
fn projector_value(sas: &mut SurfaceAnalysis<'_>, corners: &[CachePoint], p: Point, tol: f64, tol_sq: f64) -> Point2 {
    if let Some(c) = corners.iter().find(|c| c.0.distance_squared(p) < tol_sq) {
        return sas.next_value_of_uv(c.1, p, tol, tol);
    }
    sas.value_of_uv(p, tol)
}

/// `SurfaceProjectorWithCache::NextValueOfUV`.
fn projector_next_value(sas: &mut SurfaceAnalysis<'_>, corners: &[CachePoint], hint: Point2, p: Point, tol: f64, tol_sq: f64, step: f64) -> Point2 {
    if let Some(c) = corners.iter().find(|c| c.0.distance_squared(p) < tol_sq) {
        return sas.next_value_of_uv(c.1, p, tol, tol);
    }
    sas.next_value_of_uv(hint, p, tol, step)
}

/// The B-spline under any trims (`extractBSplineCurve`).
fn extract_bspline(curve: &Curve) -> Option<&super::bspline::BSplineCurve> {
    match curve {
        Curve::Trimmed(t) => extract_bspline(t.basis()),
        Curve::BSpline(b) => Some(b),
        _ => None,
    }
}

/// `adjustSecondToFirstPoint`: the second point within half a period of the first.
fn adjust_second_to_first_point(first: Point2, second: &mut Point2, surface: &Surface) {
    if let Some(period) = surface.u_period().filter(|_| surface.is_u_periodic()) {
        second.x = in_period(second.x, first.x - period / 2.0, first.x + period / 2.0);
    }
    if let Some(period) = surface.v_period().filter(|_| surface.is_v_periodic()) {
        second.y = in_period(second.y, first.y - period / 2.0, first.y + period / 2.0);
    }
}

/// `fixPeriodicityTroubles`: the four points (ends and their neighbours) brought into one
/// period of coordinate `index` (0: u), and a jump between them undone; true if there was
/// one.
fn fix_periodicity_troubles(p: &mut [Point2; 4], index: usize, period: f64, saved_point: i32, saved_param: f64) -> bool {
    let get = |p: &Point2| if index == 0 { p.x } else { p.y };
    let set = |p: &mut Point2, x: f64| if index == 0 { p.x = x } else { p.y = x };
    let (saved_param, saved_point) = if saved_point < 0 { (0.5 * period, 0usize) } else { (saved_param, saved_point as usize) };
    let (mut min_param, mut max_param) = (0.0, period);
    // (With no point saved, the middle of [0, period] stays in it.)
    while min_param > saved_param {
        min_param -= period;
        max_param -= period;
    }
    while max_param < saved_param {
        min_param += period;
        max_param += period;
    }
    let mut fix_iso_param = min_param;
    let mut iso_line = false;
    if max_param - saved_param < P_CONFUSION || saved_param - min_param < P_CONFUSION {
        fix_iso_param = saved_param;
        iso_line = true;
    }
    for q in p.iter_mut() {
        let mut param = get(q);
        param += adjust_to_period(param, min_param, max_param);
        if iso_line {
            if max_param - param < P_CONFUSION || param - min_param < P_CONFUSION {
                param = fix_iso_param;
            }
        } else {
            if max_param - param < P_CONFUSION {
                param = max_param;
            }
            if param - min_param < P_CONFUSION {
                param = min_param;
            }
        }
        set(q, param);
    }
    let mut jump = false;
    let mut prev_diff = 0.0;
    let mut sum_diff = 1.0;
    for i in 0..3 {
        let diff = get(&p[i + 1]) - get(&p[i]);
        if diff < -P_CONFUSION {
            sum_diff *= -1.0;
        }
        if diff * prev_diff < -P_CONFUSION {
            jump = true;
        }
        prev_diff = diff;
    }
    if !jump {
        return false;
    }
    if sum_diff > 0.0 {
        for i in (1..=saved_point).rev() {
            if get(&p[i]) > get(&p[i - 1]) {
                let x = get(&p[i - 1]) + period;
                set(&mut p[i - 1], x);
            }
        }
        for i in saved_point..3 {
            if get(&p[i]) < get(&p[i + 1]) {
                let x = get(&p[i + 1]) - period;
                set(&mut p[i + 1], x);
            }
        }
    } else {
        for i in (1..=saved_point).rev() {
            if get(&p[i]) < get(&p[i - 1]) {
                let x = get(&p[i - 1]) - period;
                set(&mut p[i - 1], x);
            }
        }
        for i in saved_point..3 {
            if get(&p[i]) > get(&p[i + 1]) {
                let x = get(&p[i + 1]) + period;
                set(&mut p[i + 1], x);
            }
        }
    }
    true
}

/// `isBSplineCurveInvalid`: a B-spline whose speed varies over its knot spans by more than
/// 1500 times (projected by ProjLib instead).
fn is_bspline_curve_invalid(curve: &Curve, first: f64, last: f64) -> bool {
    let Some(bs) = extract_bspline(curve) else { return false };
    let knots = bs.knots();
    let nb = knots.len();
    let mut first_param = first;
    let mut idx = 0;
    while idx < nb && first_param < last {
        if knots[idx] > first {
            break;
        }
        idx += 1;
    }
    let mut min_speed = 2.0e100;
    let mut coeffs: Vec<f64> = Vec::new();
    while idx < nb && first_param < last {
        let last_param = last.min(knots[idx]);
        let mut nb_int = N_CONTROL;
        if idx > 0 {
            let len_thres = 1.0e-2;
            let len_ratio = (last_param - first_param) / (knots[idx] - knots[idx - 1]);
            if len_ratio < len_thres {
                nb_int = ((len_ratio / len_thres * nb_int as f64) as usize).max(2);
            }
        }
        let step = (last_param - first_param) / (nb_int - 1) as f64;
        let mut p1 = curve.point(first_param);
        let mut length = 0.0;
        for k in 1..nb_int {
            let p2 = curve.point(first_param + step * k as f64);
            let d = p2.distance(p1);
            length += d;
            p1 = p2;
            min_speed = f64::min(min_speed, d / step);
        }
        let coeff = length / (last_param - first_param);
        if coeff.abs() > RESOLUTION {
            coeffs.push(coeff);
        }
        first_param = last_param;
        idx += 1;
    }
    let mut evenly = 0.0;
    if !coeffs.is_empty() {
        let max = coeffs.iter().copied().fold(f64::MIN, f64::max);
        let min = coeffs.iter().copied().fold(f64::MAX, f64::min);
        evenly = max / min;
    }
    evenly > 1500.0 && min_speed > CONFUSION
}

/// `generateCurvePoints`: the curve sampled evenly in its parameter, at least 23 points and
/// on a B-spline at least (degree + 1) per span (more by 22 at a time).
fn generate_curve_points(curve: &Curve, first: f64, last: f64, n_control: usize) -> (Vec<Point>, Vec<f64>) {
    let mut n = n_control;
    if let Some(bs) = extract_bspline(curve) {
        let knots = bs.knots();
        let used = (0..knots.len() - 1).filter(|&i| knots[i + 1] > first && knots[i] < last).count();
        let min = used * (bs.degree() + 1);
        while n < min {
            n += N_CONTROL - 1;
        }
    }
    let delta = (last - first) / (n - 1) as f64;
    let mut points = Vec::with_capacity(n);
    let mut params = Vec::with_capacity(n);
    for i in 1..=n {
        let t = if i == 1 {
            first
        } else if i == n {
            last
        } else {
            first + (i - 1) as f64 * delta
        };
        points.push(curve.point(t));
        params.push(t);
    }
    (points, params)
}

/// `checkPoints2d`: points closer to the one before than `gp::Resolution` dropped (the last
/// kept, the one before it dropped instead); `preci` set to 0.9 of the shortest step left.
fn check_points_2d(points: &mut Vec<Point2>, params: &mut Vec<f64>, preci: &mut f64) {
    let n = points.len();
    let mut keep = vec![true; n];
    let mut dropped = 0;
    let mut last_valid = 0;
    let mut dist_min2 = f64::MAX;
    let mut prev = points[0];
    for i in 1..n {
        let cur = points[i];
        let d2 = prev.distance_squared(cur);
        if d2 < RESOLUTION {
            dropped += 1;
            if i == n - 1 {
                keep[last_valid] = false;
            } else {
                keep[i] = false;
            }
        } else {
            if d2 < dist_min2 {
                dist_min2 = d2;
            }
            last_valid = i;
            prev = cur;
        }
    }
    if dist_min2 < f64::MAX {
        *preci = 0.9 * dist_min2.sqrt();
    }
    if dropped == 0 {
        return;
    }
    if n - dropped < 2 {
        // (A pcurve of minimal length.)
        keep[0] = true;
        keep[n - 1] = true;
        points[n - 1] = Point2::new(points[n - 1].x + *preci, points[n - 1].y + *preci);
        let new: Vec<(Point2, f64)> = vec![(points[0], params[0]), (points[n - 1], params[n - 1])];
        *points = new.iter().map(|x| x.0).collect();
        *params = new.iter().map(|x| x.1).collect();
        return;
    }
    let mut new_points = Vec::with_capacity(n - dropped);
    let mut new_params = Vec::with_capacity(n - dropped);
    for i in 0..n {
        if keep[i] {
            new_points.push(points[i]);
            new_params.push(params[i]);
        }
    }
    *points = new_points;
    *params = new_params;
}

/// Where the line through `p1` along `d1` meets the one through `p2` along `d2`
/// (`IntCurve_IntConicConic` of two lines, tolerance 1e-10); `None` if parallel.
fn intersect_lines(p1: Point2, d1: Vector2, p2: Point2, d2: Vector2) -> Option<Point2> {
    let d1 = d1 / d1.length();
    let d2 = d2 / d2.length();
    let cross = d1.x * d2.y - d1.y * d2.x;
    if cross.abs() <= 1.0e-10 {
        return None;
    }
    let w = p2 - p1;
    let t = (w.x * d2.y - w.y * d2.x) / cross;
    Some(p1 + d1 * t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::{BSplineCurve, BSplineSurface, KnotVector};
    use crate::geometry::curve::{Circle, Ellipse};
    use crate::geometry::surface::{Cone, Cylinder, Sphere, Torus};
    use crate::math::{Direction, Frame};
    use std::f64::consts::PI;

    /// OCCT 8.1's pcurves of the cases (its curve's kind, then nine points), from a driver
    /// running `ShapeConstruct_ProjectCurveOnSurface` on the same curves and surfaces.
    const OCCT: &[(&str, &str, [(f64, f64); 9])] = &[
        ("circle_on_cylinder", "line", [
            (0.0, 1.0), (0.78539816339744828, 1.0), (1.5707963267948966, 1.0),
            (2.3561944901923448, 1.0), (3.1415926535897931, 1.0), (3.9269908169872414, 1.0),
            (4.7123889803846897, 1.0), (5.497787143782138, 1.0), (6.2831853071795862, 1.0),
        ]),
        ("ellipse_on_cylinder", "bspline poles 25", [
            (1.5707963267948966, 0.84558643747632356), (2.3561944901923448, 0.59792100095845058), (3.1415926535897931, -2.8732900868994049e-07),
            (3.9269908169872414, -0.59791302757845866), (4.7123889803846897, -0.84558643747632356), (5.4977871437821371, -0.59791302757845877),
            (6.2831853071795862, -2.8732900886341284e-07), (7.0685834705770336, 0.59792100095845047), (7.8539816339744819, 0.84558643747632356),
        ]),
        ("ellipse_arc_on_cylinder", "bspline poles 25", [
            (2.5707963267948966, 0.45687230197928325), (2.9457963267948961, 0.16450696196313253), (3.320796326794897, -0.15072231712251716),
            (3.6957963267948966, -0.44500342099456525), (4.0707963267948957, -0.67743617577723314), (4.4457963267948966, -0.81571496858139325),
            (4.8207963267948957, -0.84062179841421081), (5.1957963267948966, -0.74869602003395364), (5.5707963267948966, -0.55271218074518635),
        ]),
        ("arc_across_seam", "line", [
            (5.0, 0.5), (5.375, 0.5), (5.75, 0.5),
            (6.125, 0.5), (6.5, 0.5), (6.875, 0.5),
            (7.25, 0.5), (7.625, 0.5), (8.0, 0.5),
        ]),
        ("latitude_on_sphere", "line", [
            (0.0, 0.52359877559829882), (0.78539816339744828, 0.52359877559829882), (1.5707963267948966, 0.52359877559829882),
            (2.3561944901923448, 0.52359877559829882), (3.1415926535897931, 0.52359877559829882), (3.9269908169872414, 0.52359877559829882),
            (4.7123889803846897, 0.52359877559829882), (5.497787143782138, 0.52359877559829882), (6.2831853071795862, 0.52359877559829882),
        ]),
        ("meridian_over_pole", "bspline poles 25", [
            (0.0, 0.49999999999999911), (-8.1007772312179015e-06, 0.75000000464868899), (0.00043600917995507542, 0.99999974979301176),
            (-0.01209115396652214, 1.250006938595235), (-5.5511151231257827e-17, 1.5), (3.0964678526646359, 1.3901107544373092),
            (3.14321986200195, 1.1416460911389334), (3.1415624210775852, 0.89159166075363006), (3.1415926535897931, 0.64159265358979334),
        ]),
        ("meridian_to_pole", "line", [
            (0.0, 0.5), (0.0, 0.63384954084936207), (0.0, 0.76769908169872414),
            (0.0, 0.90154862254808621), (0.0, 1.0353981633974483), (0.0, 1.1692477042468103),
            (0.0, 1.3030972450961724), (0.0, 1.4369467859455345), (0.0, 1.5707963267948966),
        ]),
        ("meridian_on_torus", "line", [
            (0.0, 0.0), (0.0, 0.78539816339744828), (0.0, 1.5707963267948966),
            (0.0, 2.3561944901923448), (0.0, 3.1415926535897931), (0.0, 3.9269908169872414),
            (0.0, 4.7123889803846897), (0.0, 5.497787143782138), (0.0, 6.2831853071795862),
        ]),
        ("circle_on_cone", "line", [
            (0.0, 2.0), (0.78539816339744828, 2.0), (1.5707963267948966, 2.0),
            (2.3561944901923448, 2.0), (3.1415926535897931, 2.0), (3.9269908169872414, 2.0),
            (4.7123889803846897, 2.0), (5.497787143782138, 2.0), (6.2831853071795862, 2.0),
        ]),
        ("curve_near_bspline", "bspline poles 25", [
            (0.099999999999999978, 0.19999999999999996), (0.24030457746014502, 0.23735033528006341), (0.33712169728006847, 0.28001227347008273),
            (0.41029446477625192, 0.32881657945357123), (0.47825158652311706, 0.38507448353312201), (0.54939103616863938, 0.45376009089842073),
            (0.62758683111287139, 0.54147873604631536), (0.71656669927611272, 0.65473725919794712), (0.82000000000000006, 0.80000000000000004),
        ]),
        ("curve_near_bspline_part", "bspline poles 25", [
            (0.30234128653042519, 0.26225450709740811), (0.36065503047995168, 0.2939723738656228), (0.41029119201739683, 0.32881775681296388),
            (0.4577105827164456, 0.36719253888139863), (0.50610596139256692, 0.41068672968963221), (0.55683020287861618, 0.46156700089623459),
            (0.6112016616715954, 0.52209801722692406), (0.6705022528423239, 0.5945138132545833), (0.73597551917389559, 0.6810193574631388),
        ]),
        ("bound_iso_on_bspline", "line", [
            (-9.1445772232641205e-19, -1.5295455196264383e-18), (-8.0015050703561052e-19, 0.125), (-6.8584329174480899e-19, 0.25),
            (-5.7153607645400756e-19, 0.375), (-4.5722886116320603e-19, 0.5), (-3.429216458724045e-19, 0.625),
            (-2.2861443058160306e-19, 0.75), (-1.1430721529080153e-19, 0.875), (0.0, 1.0),
        ]),
        ("inner_iso_on_bspline", "line", [
            (0.099999999999999964, 0.37000000000000005), (0.18749999999999994, 0.37000000000000005), (0.27499999999999997, 0.37000000000000005),
            (0.36249999999999999, 0.37), (0.45000000000000001, 0.37), (0.53750000000000009, 0.37),
            (0.625, 0.37), (0.71250000000000002, 0.37), (0.80000000000000004, 0.37),
        ]),
    ];

    fn frame(o: (f64, f64, f64), z: (f64, f64, f64), x: (f64, f64, f64)) -> Frame {
        Frame::new(Point::new(o.0, o.1, o.2), Direction::new(z.0, z.1, z.2), Direction::new(x.0, x.1, x.2))
    }

    /// The pcurve of `curve` checked against OCCT's: the same kind, the same points at
    /// nine parameters.
    fn run(name: &str, surface: &Surface, curve: &Curve, f: f64, l: f64, preci: f64) {
        let mut projector = PCurveProjector::new(surface, preci);
        let c2d = projector.perform(curve, f, l, -1.0, -1.0).expect(name);
        let kind = match &c2d {
            Curve2::Line(_) => "line".to_string(),
            Curve2::BSpline(b) => format!("bspline poles {}", b.poles().len()),
            _ => "other".to_string(),
        };
        let (_, occt_kind, occt) = OCCT.iter().find(|c| c.0 == name).expect("an OCCT case");
        assert_eq!(kind, *occt_kind, "{name}");
        for (i, &(x, y)) in occt.iter().enumerate() {
            let q = c2d.point(f + (l - f) * i as f64 / 8.0);
            assert!((q.x - x).abs() <= 1e-9 && (q.y - y).abs() <= 1e-9, "{name} point {i}: ({}, {}) against OCCT's ({x}, {y})", q.x, q.y);
        }
    }

    #[test]
    fn cases_against_occt() {
        let z = frame((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let cyl = Surface::Cylinder(Cylinder::new(z, 2.0));
        run("circle_on_cylinder", &cyl, &Curve::Circle(Circle::new(frame((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 2.0)), 0.0, 2.0 * PI, 1e-7);
        let a: f64 = 0.4;
        let e = Curve::Ellipse(Ellipse::new(frame((0.0, 0.0, 0.0), (0.0, -a.sin(), a.cos()), (0.0, a.cos(), a.sin())), 2.0 / a.cos(), 2.0));
        run("ellipse_on_cylinder", &cyl, &e, 0.0, 2.0 * PI, 1e-7);
        run("ellipse_arc_on_cylinder", &cyl, &e, 1.0, 4.0, 1e-7);
        run("arc_across_seam", &cyl, &Curve::Circle(Circle::new(frame((0.0, 0.0, 0.5), (0.0, 0.0, 1.0), (5.0f64.cos(), 5.0f64.sin(), 0.0)), 2.0)), 0.0, 3.0, 1e-7);
        let sph = Surface::Sphere(Sphere::new(z, 3.0));
        run("latitude_on_sphere", &sph, &Curve::Circle(Circle::new(frame((0.0, 0.0, 1.5), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), (9.0f64 - 2.25).sqrt())), 0.0, 2.0 * PI, 1e-7);
        let meridian = Curve::Circle(Circle::new(frame((0.0, 0.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0)), 3.0));
        run("meridian_over_pole", &sph, &meridian, 0.5, 2.5, 1e-7);
        run("meridian_to_pole", &sph, &meridian, 0.5, PI / 2.0, 1e-7);
        let tor = Surface::Torus(Torus::new(z, 5.0, 1.0));
        run("meridian_on_torus", &tor, &Curve::Circle(Circle::new(frame((5.0, 0.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0)), 1.0)), 0.0, 2.0 * PI, 1e-7);
        let cone = Surface::Cone(Cone::new(z, 0.3, 1.0));
        run("circle_on_cone", &cone, &Curve::Circle(Circle::new(frame((0.0, 0.0, 2.0 * 0.3f64.cos()), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)), 1.0 + 2.0 * 0.3f64.sin())), 0.0, 2.0 * PI, 1e-7);
        let poles: Vec<Vec<Point>> = (1..=4)
            .map(|i| (1..=4).map(|j| Point::new((i - 1) as f64, (j - 1) as f64, 0.3 * (1.3 * i as f64).sin() * (0.7 * j as f64).cos())).collect())
            .collect();
        let bs = Surface::BSpline(BSplineSurface::new(poles, KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3), KnotVector::new(vec![0.0, 1.0], vec![4, 4], 3)));
        let cp: Vec<Point> = (1..=5)
            .map(|k| {
                let u = 0.1 + 0.18 * (k - 1) as f64;
                let v = 0.2 + 0.6 * ((k - 1) * (k - 1)) as f64 / 16.0;
                bs.point(u, v)
            })
            .collect();
        let c = Curve::BSpline(BSplineCurve::new(cp, vec![0.0, 0.4, 1.0], vec![4, 1, 4], 3));
        run("curve_near_bspline", &bs, &c, 0.0, 1.0, 1e-4);
        run("curve_near_bspline_part", &bs, &c, 0.2, 0.9, 1e-4);
        run("bound_iso_on_bspline", &bs, &bs.u_iso(0.0), 0.0, 1.0, 1e-7);
        run("inner_iso_on_bspline", &bs, &bs.v_iso(0.37), 0.1, 0.8, 1e-7);
    }
}
