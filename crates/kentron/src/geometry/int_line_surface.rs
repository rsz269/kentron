//! OCCT's intersection of a line with the surface of a face, as `IntCurvesFace_Intersector`
//! runs it before it classifies the points against the face (`InternalCall`): the line as an
//! infinite `gp_Lin`, the surface as the face's `BRepAdaptor_Surface` (its parameters
//! restricted to the face's UV box).
//!
//! The constructor's dispatch is OCCT's: planes, cylinders, cones, spheres and tori are
//! solved in closed form when the line is cast (`IntCurveSurface_HInter::Perform(curve,
//! surface)`, which reaches `PerformConicSurfLine`: `IntAna_IntConicQuad` of the line and the
//! plane or the quadric, `IntAna_IntLinTorus` for a torus); any other surface is sampled once
//! into a polyhedron (`IntCurveSurface_ThePolyhedronOfHInter`, its grid sized as the
//! constructor sizes it from `BRepTopAdaptor_TopolTool::NbSamplesU/V`, split at the surface's
//! C2 breaks by `ComputeSamplePars`), each cast clips the line to the polyhedron's box
//! (`Intf_Tool::LinBox`), samples it into a polygon of five points, finds where the
//! polygon's segments cross the polyhedron's triangles or pass within its deflection of their
//! sides (`Intf_InterferencePolygonPolyhedron` through `Bnd_BoundSortBox`), turns each such
//! section point into starting parameters (`SectionPointToParameters`), sorts them and solves
//! `S(u, v) = L(w)` from each with `math_FunctionSetRoot` (`IntImp_IntCS`,
//! `IntImp_ZerCSParFunc`). Every solution is filtered against the line's and the face's
//! bounds and given its transition (`ComputeAppendPoint`), and appended unless one with the
//! same parameters and transition is there (`IntCurveSurface_Intersection::Append`).
//!
//! One departure: a plane or quadric is solved with the line taken to the surface's own
//! origin (see `conic_surf_line`), for precision far from the origin.
//!
//! What is not ported: the polyhedral fallback `PerformConicSurfLine` takes for a torus
//! whose every quartic root `IntAna_IntLinTorus` rejects (then no point is reported), and for
//! a cone whose half-angle is within `1e-7` of a right angle (refused: kentron has no such
//! cones); `Bnd_BoundSortBox`'s grid is a bounding-volume tree here, which finds the same
//! triangles in another order (the start points are sorted after).

use crate::math::{BoundingBox, Bvh, Point, Vector};

use super::adaptor::{epsilon, first_u_knot_index, in_period, CurveAdaptor, GeomAdaptor};
use super::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use super::int_curve_plane::{line_plane, Ax3, ConicQuad, Quadric, Transition, TOLERANCE_ANGULAIRE};
use super::polynomial_roots::PolynomialRoots;
use super::resolution::{u_resolution, v_resolution};
use super::surface::Surface;
use super::{Continuity, ParametricSurface};

/// `THE_TOLTANGENCY`: the tolerance on the line's and the face's parameter bounds.
const TOLTANGENCY: f64 = 0.00000001;
/// `Precision::Infinite()`: the bounds of the infinite line.
pub(crate) const INFINITE: f64 = 2.0e100;
/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1.0e-9;
/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = 1.0e-14;
/// `Precision::Computational()`.
const COMPUTATIONAL: f64 = f64::EPSILON;
/// `gp::Resolution()`.
const RESOLUTION: f64 = f64::MIN_POSITIVE;
/// `Precision::Angular()`.
const ANGULAR: f64 = 1.0e-12;
/// `PARAMEQUAL`'s tolerance in `IntCurveSurface_Intersection::Append`.
const PARAM_EQUAL: f64 = 1e-8;
/// `THE_MIN_EDGE_LENGTH_SQUARED` of the polyhedron.
const MIN_EDGE_LENGTH_SQUARED: f64 = 1e-15;

/// One point where the line meets the surface (`IntCurveSurface_IntersectionPoint`).
#[doc(alias = "IntCurveSurface_IntersectionPoint")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct LinePoint {
    /// `Pnt()`: the line's point at `w`.
    pub point: Point,
    /// `U()`, `V()`: the surface's parameters.
    pub u: f64,
    pub v: f64,
    /// `W()`: the distance along the line from its origin.
    pub w: f64,
    pub transition: Transition,
}

/// A face's surface made ready for lines (`IntCurvesFace_Intersector`'s constructor).
#[doc(alias = "IntCurvesFace_Intersector")]
pub(crate) struct LineFace {
    surface: Surface,
    /// The face's UV box, `(u0, u1, v0, v1)`.
    bounds: (f64, f64, f64, f64),
    /// `myPolyhedron`, with its `Bnd_BoundSortBox`, for a surface not solved in closed form.
    polyhedron: Option<Polyhedron>,
    /// `myReady`: false for a sampled surface without extent along one of its directions.
    ready: bool,
}

/// The surface `GeomAdaptor_Surface` sees: a trimmed surface's basis.
fn adapted(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => adapted(t.basis()),
        other => other,
    }
}

/// Whether OCCT solves a line against the surface in closed form.
fn analytic(surface: &Surface) -> bool {
    matches!(adapted(surface), Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_))
}

impl LineFace {
    /// `IntCurvesFace_Intersector(Face, Tol)` for a face whose surface (placed) is `surface`
    /// and whose UV box is `bounds`.
    pub(crate) fn new(surface: Surface, bounds: (f64, f64, f64, f64)) -> LineFace {
        let mut face = LineFace { surface, bounds, polyhedron: None, ready: false };
        if !analytic(&face.surface) {
            let (u0, u1, v0, v1) = bounds;
            let (mut nbsu, mut nbsv) = topol_samples(&face.surface, bounds);
            let u_res = u_resolution(&face.surface, bounds, 1.0);
            let v_res = v_resolution(&face.surface, bounds, 1.0);
            // Checking the correlation between the number of samples and the length of the
            // face along each axis.
            const TRESH: f64 = 100.0;
            let mut min_samples = 20;
            const MAX_SAMPLES: i64 = 40;
            const MAX_SAMPLES2: i64 = MAX_SAMPLES * MAX_SAMPLES;
            let du = (u1 - u0) / u_res;
            let dv = (v1 - v0) / v_res;
            nbsu = nbsu.max(min_samples).min(MAX_SAMPLES);
            nbsv = nbsv.max(min_samples).min(MAX_SAMPLES);
            if du > CONFUSION && dv > CONFUSION {
                if du.max(dv) > du.min(dv) * TRESH {
                    min_samples = 10;
                    nbsu = (((du / dv).sqrt() * MAX_SAMPLES as f64) as i64).max(min_samples);
                    nbsv = MAX_SAMPLES2 / nbsu;
                    if nbsv < min_samples {
                        nbsv = min_samples;
                        nbsu = MAX_SAMPLES2 / min_samples;
                    }
                }
            } else {
                // The surface has no extension along one of its directions.
                return face;
            }
            let u_breaks = intervals(&face.surface, bounds, true);
            let v_breaks = intervals(&face.surface, bounds, false);
            let (us, vs) = if u_breaks.len() > 2 || v_breaks.len() > 2 {
                (sample_pars(&u_breaks, nbsu), sample_pars(&v_breaks, nbsv))
            } else {
                // (The uniform constructor: at least three subdivisions each way.)
                let (nu, nv) = (nbsu.max(3) as usize, nbsv.max(3) as usize);
                let du = (u1 - u0) / nu as f64;
                let dv = (v1 - v0) / nv as f64;
                ((0..=nu).map(|i| u0 + i as f64 * du).collect(), (0..=nv).map(|j| v0 + j as f64 * dv).collect())
            };
            face.polyhedron = Some(Polyhedron::new(&face.surface, &us, &vs));
        }
        face.ready = true;
        face
    }

    /// `Perform(gp_Lin(origin, direction), ParMin, ParMax)` up to `InternalCall`'s
    /// classification: the points the surface meets the line at, and the range of `w`
    /// `InternalCall` keeps them in. `None` where OCCT is not done (`myReady` false).
    /// `direction` is a unit vector.
    pub(crate) fn perform(&self, origin: Point, direction: Vector, par_min: f64, par_max: f64) -> Option<(Vec<LinePoint>, f64, f64)> {
        if !self.ready {
            return None;
        }
        let line = LineCurve { origin, direction };
        let Some(polyhedron) = &self.polyhedron else {
            return Some((self.conic_surf_line(&line), par_min, par_max));
        };
        let empty = Some((Vec::new(), par_min, par_max));
        let Some((begin, end)) = lin_box(origin, direction, &polyhedron.bnd) else { return empty };
        let mut pinf = begin;
        let mut psup = end;
        let pppp = 0.05 * (psup - pinf);
        pinf -= pppp;
        psup += pppp;
        if psup - pinf < 1e-10 {
            pinf -= 1e-10;
            psup += 1e-10;
        }
        let (mut parinf, mut parsup) = (pinf, psup);
        if parinf > par_max || parsup < par_min {
            return empty;
        }
        parinf = parinf.max(par_min);
        parsup = parsup.min(par_max);
        if parinf > parsup - 1e-9 {
            return empty;
        }
        let polygon = Polygon::new(&line, parinf, parsup, 2);
        Some((self.polygon_polyhedron(&line, &polygon, polyhedron), parinf, parsup))
    }

    /// `PerformConicSurfLine` for the surfaces it solves in closed form.
    fn conic_surf_line(&self, line: &LineCurve) -> Vec<LinePoint> {
        let d = [line.direction.x, line.direction.y, line.direction.z];
        // (A departure: the line and the plane or quadric are taken to the surface's own
        // origin first, as `IntAna_IntLinTorus` takes the line to the torus's, so that a
        // surface far from the origin is solved as precisely as one near it; OCCT's
        // coefficients, of the size of the squared distance to the origin, lose the roots'
        // digits there. The parameters on the line are unchanged by the move.)
        let (frame, mut out) = match adapted(&self.surface) {
            Surface::Plane(p) => (*p.frame(), Vec::new()),
            Surface::Cylinder(c) => (*c.frame(), Vec::new()),
            Surface::Sphere(s) => (*s.frame(), Vec::new()),
            Surface::Cone(c) => (*c.frame(), Vec::new()),
            _ => (crate::math::Frame::WORLD, Vec::new()),
        };
        let loc = frame.origin();
        let rel = line.origin - loc;
        let o = [rel.x, rel.y, rel.z];
        let mut pos = Ax3::from_frame(&frame);
        pos.location = [0.0, 0.0, 0.0];
        match adapted(&self.surface) {
            Surface::Plane(_) => {
                let ana = line_plane(o, d, &pos, TOLERANCE_ANGULAIRE, 0.0, 0.0);
                self.append_int_ana(line, loc, &ana, &mut out);
            }
            Surface::Cylinder(c) => {
                let ana = line_quadric(o, d, &cylinder_quadric(&pos, c.radius()));
                self.append_int_ana(line, loc, &ana, &mut out);
            }
            Surface::Sphere(s) => {
                let ana = line_quadric(o, d, &sphere_quadric(&pos, s.radius()));
                self.append_int_ana(line, loc, &ana, &mut out);
            }
            Surface::Torus(t) => {
                // (A torus whose roots `IntAna_IntLinTorus` all rejects falls to a polyhedron
                // in OCCT: not ported, no point.)
                for (w, u, v) in line_torus(line, t) {
                    if let Some(p) = self.compute_append_point(line, w, u, v) {
                        append(&mut out, p);
                    }
                }
            }
            Surface::Cone(c) => {
                let correction = 1.0e5 * ANGULAR;
                if c.semi_angle().abs() < std::f64::consts::FRAC_PI_2 - correction {
                    let ana = line_quadric(o, d, &cone_quadric(&pos, c.semi_angle(), c.radius()));
                    self.append_int_ana(line, loc, &ana, &mut out);
                }
            }
            _ => unreachable!("an analytic surface"),
        }
        out
    }

    /// `AppendIntAna` (`ProcessIntAna`): a line parallel to the surface or in it gives
    /// nothing; otherwise each solution's parameters are found (`ComputeParamsOnQuadric`) and
    /// it is filtered and appended.
    fn append_int_ana(&self, line: &LineCurve, loc: Point, ana: &ConicQuad, out: &mut Vec<LinePoint>) {
        if !ana.done || ana.inquadric || ana.parallel {
            return;
        }
        let mut points = Vec::new();
        for &(p, w) in &ana.points {
            let (u, v) = quadric_parameters(adapted(&self.surface), loc + Vector::new(p[0], p[1], p[2]));
            if let Some(point) = self.compute_append_point(line, w, u, v) {
                points.push(point);
            }
        }
        for p in points {
            append(out, p);
        }
    }

    /// `ComputeAppendPoint`: the solution brought into the surface's periods, rejected
    /// outside the line's or the face's bounds, given its transition and the line's point.
    fn compute_append_point(&self, line: &LineCurve, lw: f64, su: f64, sv: f64) -> Option<LinePoint> {
        let (w0, w1) = (-INFINITE, INFINITE);
        let (u0, u1, v0, v1) = self.bounds;
        let (w, mut u, mut v) = (lw, su, sv);
        if w0 - w >= TOLTANGENCY || w - w1 >= TOLTANGENCY {
            return None;
        }
        let s = adapted(&self.surface);
        let angular_u = matches!(s, Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_));
        if self.surface.is_u_periodic() || angular_u {
            let period = self.surface.u_period().unwrap_or(2.0 * std::f64::consts::PI);
            u = in_period(u, u0, u0 + period);
        }
        if self.surface.is_v_periodic() {
            let period = self.surface.v_period().unwrap_or(2.0 * std::f64::consts::PI);
            v = in_period(v, v0, v0 + period);
        }
        if u0 - u >= TOLTANGENCY || u - u1 >= TOLTANGENCY {
            return None;
        }
        if v0 - v >= TOLTANGENCY || v - v1 >= TOLTANGENCY {
            return None;
        }
        let transition = compute_transition(line, &self.surface, u, v);
        Some(LinePoint { point: line.value(w), u, v, w, transition })
    }

    /// `HICS.Perform(HLL, polygon, Hsurface, *myPolyhedron, *myBndBounding)`
    /// (`PerformPolygonPolyhedronBSB`, `InternalPerformBSB`).
    fn polygon_polyhedron(&self, line: &LineCurve, polygon: &Polygon, polyhedron: &Polyhedron) -> Vec<LinePoint> {
        let (u0, u1, v0, v1) = self.bounds;
        let sections = interference(polygon, polyhedron);
        let mut starts: Vec<(f64, f64, f64)> = sections.iter().map(|sp| section_point_parameters(sp, polyhedron, polygon)).collect();
        sort_start_points(&mut starts);
        let mut function = CsFunction { surface: &self.surface, line, f: 0.0 };
        // `IntImp_IntCS(F, TolTangency)`, `math_FunctionSetRoot(F)` (100 iterations).
        let tol = (TOLTANGENCY * TOLTANGENCY).max(SQUARE_CONFUSION);
        let mut rsnld = FunctionSetRoot::new(3, 100);
        let (winf, wsup) = (polygon.binf, polygon.bsup);
        let mut out = Vec::new();
        // `ProcessSortedPoints`.
        let ptol = 10.0 * P_CONFUSION;
        let (mut su, mut sv, mut sw) = (0.0, 0.0, 0.0);
        for (i, &(u, v, w)) in starts.iter().enumerate() {
            if i == 0 {
                su = u - 1.0;
            }
            if (u - su).abs() > ptol || (v - sv).abs() > ptol || (w - sw).abs() > ptol {
                if let Some((w, u, v)) = exact_perform(&mut function, &mut rsnld, tol, [u, v, w], [u0, v0, winf], [u1, v1, wsup], self.bounds) {
                    if let Some(p) = self.compute_append_point(line, w, u, v) {
                        out.push(p);
                    }
                }
            }
            (su, sv, sw) = (u, v, w);
        }
        let mut appended = Vec::new();
        for p in out {
            append(&mut appended, p);
        }
        appended
    }
}

/// The line `origin + w · direction` (a `GeomAdaptor_Curve` of a `Geom_Line`).
struct LineCurve {
    origin: Point,
    direction: Vector,
}

impl LineCurve {
    fn value(&self, w: f64) -> Point {
        self.origin + self.direction * w
    }
}

/// `IntCurveSurface_Intersection::Append(Point)`: skipped when a point with the same
/// parameters (within 1e-8) and transition is there already.
fn append(points: &mut Vec<LinePoint>, point: LinePoint) {
    let equal = |a: f64, b: f64| (a - b).abs() < PARAM_EQUAL;
    if points.iter().any(|o| equal(point.u, o.u) && equal(point.v, o.v) && equal(point.w, o.w) && o.transition == point.transition) {
        return;
    }
    points.push(point);
}

/// `ComputeTransitions`: the line's direction against the surface's normal `D1U ^ D1V`.
fn compute_transition(line: &LineCurve, surface: &Surface, u: f64, v: f64) -> Transition {
    let (_, d1u, d1v) = surface.d1(u, v);
    let nsurf = d1u.cross(d1v);
    let tangent = line.direction;
    let norm = nsurf.length();
    if norm > TOLERANCE_ANGULAIRE && tangent.length_squared() > TOLERANCE_ANGULAIRE {
        let t = tangent * (1.0 / tangent.length());
        let cos_dir = nsurf.dot(t) / norm;
        if -cos_dir > TOLERANCE_ANGULAIRE {
            Transition::In
        } else if cos_dir > TOLERANCE_ANGULAIRE {
            Transition::Out
        } else {
            Transition::Tangent
        }
    } else {
        Transition::Tangent
    }
}

/// `ComputeParamsOnQuadric` (`ElSLib::Parameters`).
fn quadric_parameters(surface: &Surface, p: Point) -> (f64, f64) {
    match surface {
        Surface::Plane(s) => s.parameters(p),
        Surface::Cylinder(s) => s.parameters(p),
        Surface::Cone(s) => s.parameters(p),
        Surface::Sphere(s) => s.parameters(p),
        _ => (0.0, 0.0),
    }
}

/// The rows of `gp_Trsf::SetTransformation(pos)`: `T(i, 1..=3)` and `T(i, 4)`.
fn trsf_rows(pos: &Ax3) -> [[f64; 4]; 3] {
    let (m, loc) = pos.transformation();
    [0, 1, 2].map(|i| [m[i][0], m[i][1], m[i][2], loc[i]])
}

/// `IntAna_Quadric(gp_Cylinder)` (`gp_Cylinder::Coefficients`).
fn cylinder_quadric(pos: &Ax3, radius: f64) -> Quadric {
    let t = trsf_rows(pos);
    let ([t11, t12, t13, t14], [t21, t22, t23, t24]) = (t[0], t[1]);
    Quadric {
        cxx: t11 * t11 + t21 * t21,
        cyy: t12 * t12 + t22 * t22,
        czz: t13 * t13 + t23 * t23,
        cxy: t11 * t12 + t21 * t22,
        cxz: t11 * t13 + t21 * t23,
        cyz: t12 * t13 + t22 * t23,
        cx: t11 * t14 + t21 * t24,
        cy: t12 * t14 + t22 * t24,
        cz: t13 * t14 + t23 * t24,
        ccte: t14 * t14 + t24 * t24 - radius * radius,
    }
}

/// `IntAna_Quadric(gp_Cone)` (`gp_Cone::Coefficients`).
fn cone_quadric(pos: &Ax3, semi_angle: f64, radius: f64) -> Quadric {
    let t = trsf_rows(pos);
    let k = semi_angle.tan();
    let ([t11, t12, t13, t14], [t21, t22, t23, t24]) = (t[0], t[1]);
    let [t31, t32, t33, t34] = t[2].map(|x| x * k);
    Quadric {
        cxx: t11 * t11 + t21 * t21 - t31 * t31,
        cyy: t12 * t12 + t22 * t22 - t32 * t32,
        czz: t13 * t13 + t23 * t23 - t33 * t33,
        cxy: t11 * t12 + t21 * t22 - t31 * t32,
        cxz: t11 * t13 + t21 * t23 - t31 * t33,
        cyz: t12 * t13 + t22 * t23 - t32 * t33,
        cx: t11 * t14 + t21 * t24 - t31 * (radius + t34),
        cy: t12 * t14 + t22 * t24 - t32 * (radius + t34),
        cz: t13 * t14 + t23 * t24 - t33 * (radius + t34),
        ccte: t14 * t14 + t24 * t24 - radius * radius - t34 * t34 - 2.0 * radius * t34,
    }
}

/// `IntAna_Quadric(gp_Sphere)` (`gp_Sphere::Coefficients`).
fn sphere_quadric(pos: &Ax3, radius: f64) -> Quadric {
    let t = trsf_rows(pos);
    let ([t11, t12, t13, t14], [t21, t22, t23, t24], [t31, t32, t33, t34]) = (t[0], t[1], t[2]);
    Quadric {
        cxx: t11 * t11 + t21 * t21 + t31 * t31,
        cyy: t12 * t12 + t22 * t22 + t32 * t32,
        czz: t13 * t13 + t23 * t23 + t33 * t33,
        cxy: t11 * t12 + t21 * t22 + t31 * t32,
        cxz: t11 * t13 + t21 * t23 + t31 * t33,
        cyz: t12 * t13 + t22 * t23 + t32 * t33,
        cx: t11 * t14 + t21 * t24 + t31 * t34,
        cy: t12 * t14 + t22 * t24 + t32 * t34,
        cz: t13 * t14 + t23 * t24 + t33 * t34,
        ccte: t14 * t14 + t24 * t24 + t34 * t34 - radius * radius,
    }
}

/// `IntAna_IntConicQuad::Perform(gp_Lin, IntAna_Quadric)`: the line substituted into the
/// quadric, a polynomial of degree two in its parameter.
#[doc(alias = "IntAna_IntConicQuad")]
fn line_quadric(origin: [f64; 3], direction: [f64; 3], quad: &Quadric) -> ConicQuad {
    let Quadric { cxx: qxx, cyy: qyy, czz: qzz, cxy: qxy, cxz: qxz, cyz: qyz, cx: qx, cy: qy, cz: qz, ccte: qcte } = *quad;
    let [lx, ly, lz] = direction;
    let [lx0, ly0, lz0] = origin;
    let a0 = qcte + qxx * lx0 * lx0 + qyy * ly0 * ly0 + qzz * lz0 * lz0 + 2.0 * (lx0 * (qx + qxy * ly0 + qxz * lz0) + ly0 * (qy + qyz * lz0) + qz * lz0);
    let a1 = 2.0 * (lx * (qx + qxx * lx0 + qxy * ly0 + qxz * lz0) + ly * (qy + qxy * lx0 + qyy * ly0 + qyz * lz0) + lz * (qz + qxz * lx0 + qyz * ly0 + qzz * lz0));
    let a2 = qxx * lx * lx + qyy * ly * ly + qzz * lz * lz + 2.0 * (lx * (qxy * ly + qxz * lz) + qyz * ly * lz);
    let roots = PolynomialRoots::quadratic(a2, a1, a0);
    let mut out = ConicQuad::default();
    if roots.done {
        out.done = true;
        if roots.infinite {
            out.inquadric = true;
        } else {
            for &t in &roots.roots {
                out.points.push(([lx0 + lx * t, ly0 + ly * t, lz0 + lz * t], t));
            }
        }
    }
    out
}

/// `IntAna_IntLinTorus::Perform(L, T)` and `ProcessLinTorus`'s reading of it: the line's
/// parameter and the torus's at each root of the quartic that lies on the torus (within
/// `1e-5`). Empty where it is not done.
#[doc(alias = "IntAna_IntLinTorus")]
fn line_torus(line: &LineCurve, torus: &super::surface::Torus) -> Vec<(f64, f64, f64)> {
    let pl = line.origin;
    let dl = line.direction;
    // Reparametrize the line: its location nearest the torus's.
    let tor_loc = torus.frame().origin();
    let param_of_new_pl = (tor_loc - pl).dot(dl);
    let new_pl = pl + dl * param_of_new_pl;
    // The line in the torus's frame.
    let pos = Ax3::from_frame(torus.frame());
    let (m, loc) = pos.transformation();
    let q = [new_pl.x, new_pl.y, new_pl.z];
    let d = [dl.x, dl.y, dl.z];
    let row = |i: usize, v: [f64; 3]| m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    let (x0, y0, z0) = (row(0, q) + loc[0], row(1, q) + loc[1], row(2, q) + loc[2]);
    let (x1, y1, z1) = (row(0, d), row(1, d), row(2, d));
    let r_major = torus.major_radius();
    let r2 = r_major * r_major;
    let r = torus.minor_radius();
    let rr2 = r * r;
    let a = x1 * x1 + y1 * y1 + z1 * z1;
    let b = 2.0 * (x1 * x0 + y1 * y0 + z1 * z0);
    let c = x0 * x0 + y0 * y0 + z0 * z0 - (r2 + rr2);
    let a4 = a * a;
    let a3 = 2.0 * a * b;
    let a2 = 2.0 * a * c + 4.0 * r2 * z1 * z1 + b * b;
    let a1 = 2.0 * b * c + 8.0 * r2 * z1 * z0;
    let a0 = c * c + 4.0 * r2 * (z0 * z0 - rr2);
    let roots = PolynomialRoots::quartic(a4, a3, a2, a1, a0);
    if !roots.done {
        return Vec::new();
    }
    let mut out = Vec::new();
    for &t in &roots.roots {
        let t = t + param_of_new_pl;
        let psol_l = line.value(t);
        let (u, v) = torus.parameters(psol_l);
        let psol_t = Surface::Torus(torus.clone()).point(u, v);
        if psol_t.distance_squared(psol_l) <= 0.0000000001 {
            out.push((t, u, v));
        }
    }
    out
}

/// `BRepTopAdaptor_TopolTool::NbSamplesU()`, `NbSamplesV()` (`ComputeSamplePoints`) for a
/// surface OCCT samples into a polyhedron.
#[doc(alias = "BRepTopAdaptor_TopolTool::ComputeSamplePoints")]
fn topol_samples(surface: &Surface, _bounds: (f64, f64, f64, f64)) -> (i64, i64) {
    let s = adapted(surface);
    let (mut nbsu, mut nbsv): (i64, i64) = match s {
        Surface::Bezier(b) => (3 + b.u_pole_count() as i64, 3 + b.v_pole_count() as i64),
        Surface::BSpline(b) => {
            let nbsv = ((b.v_knots().len() * b.v_degree()) as i64).max(4);
            let nbsu = ((b.u_knots().len() * b.u_degree()) as i64).max(4);
            (nbsu, nbsv)
        }
        Surface::Revolution(_) | Surface::Extrusion(_) => (25, 15),
        _ => (10, 10),
    };
    nbsu = nbsu.max(10);
    nbsv = nbsv.max(10);
    if nbsu > 10 || nbsv > 10 {
        let poles = match s {
            Surface::BSpline(b) => Some(b.poles()),
            Surface::Bezier(b) => Some(b.poles()),
            _ => None,
        };
        if let Some(poles) = poles {
            (nbsu, nbsv) = analyse(poles);
        }
    }
    (nbsu.max(10), nbsv.max(10))
}

/// `Analyse(array2, nbup, nbvp, ...)`: from the turns of the second differences of the
/// poles along each direction, `(NbSamplesU, NbSamplesV)`.
fn analyse(poles: &[Vec<Point>]) -> (i64, i64) {
    let nbup = poles.len();
    let nbvp = poles.first().map_or(0, |r| r.len());
    let at = |i: usize, j: usize| poles[i - 1][j - 1];
    let second = |a: Point, b: Point, c: Point| Vector::new(c.x - b.x - b.x + a.x, c.y - b.y - b.y + a.y, c.z - b.z - b.z + a.z);
    let mut sh = 1;
    let mut nbch = 0;
    if nbvp > 2 {
        for i in 2..nbup {
            let mut vi = second(at(i, 1), at(i, 2), at(i, 3));
            let mut locnbch = 0;
            for j in 3..nbvp {
                let vip1 = second(at(i, j - 1), at(i, j), at(i, j + 1));
                let pd = vi.dot(vip1);
                vi = vip1;
                if !(-1.0e-7..=1.0e-7).contains(&pd) {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            nbch = nbch.max(locnbch);
        }
    }
    let samples_v = nbch + 5;
    nbch = 0;
    if nbup > 2 {
        for j in 2..nbvp {
            let mut vi = second(at(1, j), at(2, j), at(3, j));
            let mut locnbch = 0;
            for i in 3..nbup {
                let vip1 = second(at(i - 1, j), at(i, j), at(i + 1, j));
                let pd = vi.dot(vip1);
                vi = vip1;
                if !(-1.0e-7..=1.0e-7).contains(&pd) {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            nbch = nbch.max(locnbch);
        }
    }
    (nbch + 5, samples_v)
}

/// The bounds of the face adaptor's C2 intervals along u (`along_u`) or v
/// (`GeomAdaptor_Surface::UIntervals` / `VIntervals`), the first and last bounds included.
fn intervals(surface: &Surface, bounds: (f64, f64, f64, f64), along_u: bool) -> Vec<f64> {
    intervals_of(surface, bounds, along_u, Continuity::C2)
}

fn intervals_of(surface: &Surface, bounds: (f64, f64, f64, f64), along_u: bool, continuity: Continuity) -> Vec<f64> {
    let (u0, u1, v0, v1) = bounds;
    let (first, last) = if along_u { (u0, u1) } else { (v0, v1) };
    let whole = vec![first, last];
    let bspline_curve = |c: &super::curve::Curve| matches!(c, super::curve::Curve::BSpline(_));
    match surface {
        Surface::Trimmed(t) => intervals_of(t.basis(), bounds, along_u, continuity),
        Surface::BSpline(b) => {
            // `GeomAdaptor_Curve` on the iso at the first knot of the other direction.
            let iso = if along_u {
                let k = b.v_knot_vector();
                b.v_iso(k.knots[first_u_knot_index(k.degree, &k.multiplicities) - 1])
            } else {
                let k = b.u_knot_vector();
                b.u_iso(k.knots[first_u_knot_index(k.degree, &k.multiplicities) - 1])
            };
            GeomAdaptor::new(&iso, first, last).intervals(continuity)
        }
        Surface::Extrusion(e) if along_u && bspline_curve(e.basis()) => GeomAdaptor::new(e.basis(), first, last).intervals(continuity),
        Surface::Revolution(r) if !along_u && bspline_curve(r.basis()) => GeomAdaptor::new(r.basis(), first, last).intervals(continuity),
        Surface::Offset(o) => {
            let base = match continuity {
                Continuity::C0 => Continuity::C1,
                Continuity::C1 => Continuity::C2,
                Continuity::C2 => Continuity::C3,
                other => other,
            };
            intervals_of(o.basis(), bounds, along_u, base)
        }
        _ => whole,
    }
}

/// `ComputeSamplePars`: each C2 interval sampled in proportion to its length, `nbs` over
/// the whole range.
fn sample_pars(breaks: &[f64], nbs: i64) -> Vec<f64> {
    let n = breaks.len() - 1;
    let t = 1.0 / (breaks[n] - breaks[0]);
    let mut out = Vec::new();
    for i in 0..n {
        let dt = breaks[i + 1] - breaks[i];
        // (`RealToInt`: rounds half away from zero.)
        let sub = (nbs as f64 * dt * t).round() as i64 + 1;
        out.push(breaks[i]);
        let step = dt / sub as f64;
        let mut x = breaks[i];
        for _ in 1..sub {
            x += step;
            out.push(x);
        }
    }
    out.push(breaks[n]);
    out
}

/// A `Bnd_Box`: its extent, its gap already added.
#[derive(Clone, Copy, Debug)]
struct Bnd {
    lo: Point,
    hi: Point,
}

impl Bnd {
    fn of(points: impl IntoIterator<Item = Point>, gap: f64) -> Bnd {
        let b = BoundingBox::from_points(points);
        let g = Vector::new(gap, gap, gap);
        Bnd { lo: b.min - g, hi: b.max + g }
    }

    /// `IsOut(Other)` false.
    fn meets(&self, other: &Bnd) -> bool {
        !(self.lo.x > other.hi.x || other.lo.x > self.hi.x || self.lo.y > other.hi.y || other.lo.y > self.hi.y || self.lo.z > other.hi.z || other.lo.z > self.hi.z)
    }
}

/// `Intf_Tool::LinBox(L, domain)`: the line's parameters where it enters and leaves the
/// box, slab by slab (`None` where a slab the line runs along does not hold it). As OCCT's,
/// the segment is reported even when it is empty (enters after it leaves).
#[doc(alias = "Intf_Tool::LinBox")]
fn lin_box(origin: Point, direction: Vector, domain: &Bnd) -> Option<(f64, f64)> {
    let (mut parmin, mut parmax) = (-INFINITE, INFINITE);
    let o = [origin.x, origin.y, origin.z];
    let d = [direction.x, direction.y, direction.z];
    let lo = [domain.lo.x, domain.lo.y, domain.lo.z];
    let hi = [domain.hi.x, domain.hi.y, domain.hi.z];
    for k in 0..3 {
        let (a, b) = if d[k] > 0.0 {
            ((lo[k] - o[k]) / d[k], (hi[k] - o[k]) / d[k])
        } else if d[k] < 0.0 {
            ((hi[k] - o[k]) / d[k], (lo[k] - o[k]) / d[k])
        } else {
            if o[k] < lo[k] || hi[k] < o[k] {
                return None;
            }
            continue;
        };
        if k == 0 {
            (parmin, parmax) = (a, b);
        } else {
            parmin = parmin.max(a);
            parmax = parmax.min(b);
        }
    }
    Some((parmin, parmax))
}

/// `IntCurveSurface_ThePolygonOfHInter(Curve, U1, U2, NbPnt)` of the line.
struct Polygon {
    pnts: Vec<Point>,
    binf: f64,
    bsup: f64,
    deflection: f64,
    bnd: Bnd,
}

impl Polygon {
    fn new(line: &LineCurve, u1: f64, u2: f64, nb_pnt: usize) -> Polygon {
        let n = nb_pnt.max(5);
        // `InitUniform`.
        let du = (u2 - u1) / (n - 1) as f64;
        let mut u = u1;
        let mut pnts = Vec::with_capacity(n);
        for _ in 0..n {
            pnts.push(line.value(u));
            u += du;
        }
        let mut deflection = 0.0f64;
        let mut u = u1 + du * 0.5;
        for i in 0..n - 1 {
            let pm = line.value(u);
            let (p1, p2) = (pnts[i], pnts[i + 1]);
            let dir = (p2 - p1) * (1.0 / (p2 - p1).length());
            let w = pm - p1;
            deflection = deflection.max((w - dir * w.dot(dir)).length());
            u += du;
        }
        let bnd = Bnd::of(pnts.iter().copied(), 1.5 * deflection);
        Polygon { pnts, binf: u1, bsup: u2, deflection, bnd }
    }

    fn nb_segments(&self) -> usize {
        self.pnts.len() - 1
    }

    /// `ApproxParamOnCurve(Index, ParamOnLine)`.
    fn approx_param_on_curve(&self, index: usize, param_on_line: f64) -> f64 {
        let n = self.pnts.len();
        if !(0.0..=1.0).contains(&param_on_line) {
            return self.binf + (param_on_line * (self.bsup - self.binf)) / (n - 1) as f64;
        }
        let (mut index, mut param) = (index, param_on_line);
        if index == n && param == 0.0 {
            index -= 1;
            param = 1.0;
        }
        let du = (self.bsup - self.binf) / (n - 1) as f64;
        let u = self.binf + du * (index - 1) as f64;
        u + du * param
    }
}

/// `IntCurveSurface_ThePolyhedronOfHInter` over a grid of parameters, with its triangles'
/// boxes in a tree (`Bnd_BoundSortBox`). Points and triangles are numbered from 1.
struct Polyhedron {
    nbdu: usize,
    nbdv: usize,
    pnts: Vec<Point>,
    u: Vec<f64>,
    v: Vec<f64>,
    on_bounds: Vec<bool>,
    deflection: f64,
    border_deflection: f64,
    bnd: Bnd,
    /// The triangles with a box (degenerate ones have none), by the tree's index.
    boxed: Vec<usize>,
    tree: Bvh,
}

impl Polyhedron {
    /// `Init(Surface, Upars, Vpars)`.
    fn new(surface: &Surface, us: &[f64], vs: &[f64]) -> Polyhedron {
        let (nbdu, nbdv) = (us.len() - 1, vs.len() - 1);
        let mut pnts = vec![Point::new(0.0, 0.0, 0.0)];
        let (mut uu, mut vv, mut on_bounds) = (vec![0.0], vec![0.0], vec![false]);
        let grid = surface.points_on_grid(us, vs);
        for (i1, &u) in us.iter().enumerate() {
            for (i2, &v) in vs.iter().enumerate() {
                pnts.push(grid[i2][i1]);
                uu.push(u);
                vv.push(v);
                on_bounds.push(i1 == 0 || i1 == nbdu || i2 == 0 || i2 == nbdv);
            }
        }
        let mut poly = Polyhedron { nbdu, nbdv, pnts, u: uu, v: vv, on_bounds, deflection: epsilon(100.0), border_deflection: 0.0, bnd: Bnd { lo: Point::new(0.0, 0.0, 0.0), hi: Point::new(0.0, 0.0, 0.0) }, boxed: Vec::new(), tree: Bvh::new(&[]) };
        // `ComputeMaxDeflection`, `DeflectionOverEstimation(tol * 1.2)`.
        let mut tol = 0.0f64;
        for t in 1..=poly.nb_triangles() {
            let [i1, i2, i3] = poly.triangle(t);
            let uc = (poly.u[i1] + poly.u[i2] + poly.u[i3]) / 3.0;
            let vc = (poly.v[i1] + poly.v[i2] + poly.v[i3]) / 3.0;
            let center = surface.point(uc, vc);
            tol = tol.max(deflection_with_center(poly.pnts[i1], poly.pnts[i2], poly.pnts[i3], center));
        }
        const MIN_DEFLECTION: f64 = 0.0001;
        let flec = tol * 1.2;
        poly.deflection = if flec < MIN_DEFLECTION { MIN_DEFLECTION } else { flec };
        poly.bnd = Bnd::of(poly.pnts[1..].iter().copied(), poly.deflection);
        // `FillBounding`.
        let mut boxes = Vec::new();
        for t in 1..=poly.nb_triangles() {
            let [a, b, c] = poly.triangle(t).map(|i| poly.pnts[i]);
            if a.distance_squared(b) > MIN_EDGE_LENGTH_SQUARED && a.distance_squared(c) > MIN_EDGE_LENGTH_SQUARED && b.distance_squared(c) > MIN_EDGE_LENGTH_SQUARED {
                let bx = Bnd::of([a, b, c], poly.deflection);
                boxes.push(BoundingBox::from_points([bx.lo, bx.hi]));
                poly.boxed.push(t);
            }
        }
        poly.tree = Bvh::new(&boxes);
        // `ComputeMaxBorderDeflection`.
        let (u0, u1, v0, v1) = (us[0], us[nbdu], vs[0], vs[nbdv]);
        let mut border = f64::MIN;
        border = border.max(border_deflection(surface, u0, v0, v1, true, nbdv));
        border = border.max(border_deflection(surface, u1, v0, v1, true, nbdv));
        border = border.max(border_deflection(surface, v0, u0, u1, false, nbdu));
        border = border.max(border_deflection(surface, v1, u0, u1, false, nbdu));
        poly.border_deflection = border;
        poly
    }

    fn nb_triangles(&self) -> usize {
        self.nbdu * self.nbdv * 2
    }

    /// `Triangle(Index, P1, P2, P3)`.
    fn triangle(&self, index: usize) -> [usize; 3] {
        let nbdv = self.nbdv;
        let line = 1 + (index - 1) / (nbdv * 2);
        let colon = 1 + (index - 1) % (nbdv * 2);
        let colpnt = colon.div_ceil(2);
        let p1 = (line - 1) * (nbdv + 1) + colpnt;
        let p2 = line * (nbdv + 1) + colpnt + (colon - 1) % 2;
        let p3 = (line - 1 + colon % 2) * (nbdv + 1) + colpnt + 1;
        [p1, p2, p3]
    }

    /// `IsOnBound(Index1, Index2)`.
    fn is_on_bound(&self, i1: usize, i2: usize) -> bool {
        let diff = i1.abs_diff(i2);
        let nbdv = self.nbdv;
        if diff != 1 && diff != nbdv + 1 {
            return false;
        }
        for i in 0..=self.nbdu {
            if i1 == 1 + i * (nbdv + 1) && i2 + 1 == i1 {
                return false;
            }
            if i1 == (1 + i) * (nbdv + 1) && i2 == i1 + 1 {
                return false;
            }
        }
        self.on_bounds[i1] && self.on_bounds[i2]
    }
}

/// `ComputeDeflectionWithCenter`.
fn deflection_with_center(p1: Point, p2: Point, p3: Point, center: Point) -> f64 {
    if p1.distance_squared(p2) <= MIN_EDGE_LENGTH_SQUARED || p1.distance_squared(p3) <= MIN_EDGE_LENGTH_SQUARED || p2.distance_squared(p3) <= MIN_EDGE_LENGTH_SQUARED {
        return 0.0;
    }
    let (x1, x2, x3) = (p2 - p1, p3 - p2, p1 - p3);
    let n = x1.cross(x2) + x2.cross(x3) + x3.cross(x1);
    let len = n.length();
    if len < RESOLUTION {
        return 0.0;
    }
    ((center - p1).dot(n) / len).abs()
}

/// `ComputeBorderDeflection`: the largest distance of an isoline's midpoints from the chords
/// of its samples.
fn border_deflection(surface: &Surface, parameter: f64, pmin: f64, pmax: f64, u_iso: bool, nb_samples: usize) -> f64 {
    if nb_samples == 0 {
        return 0.0;
    }
    let delta = (pmax - pmin) / nb_samples as f64;
    let at = |t: f64| if u_iso { surface.point(parameter, t) } else { surface.point(t, parameter) };
    let mut deflection = f64::MIN;
    for i in 0..nb_samples {
        let p1 = at(pmin + i as f64 * delta);
        let p2 = at(pmin + (i + 1) as f64 * delta);
        let pm = at(pmin + (i as f64 + 0.5) * delta);
        let mid = Point::new((p2.x + p1.x) / 2.0, (p2.y + p1.y) / 2.0, (p2.z + p1.z) / 2.0);
        deflection = deflection.max(mid.distance(pm));
    }
    deflection
}

/// `Intf::PlaneEquation`: the triangle's unit normal and polar distance.
fn plane_equation(p1: Point, p2: Point, p3: Point) -> (Vector, f64) {
    let (v1, v2, v3) = (p2 - p1, p3 - p2, p1 - p3);
    if v1.length_squared() <= MIN_EDGE_LENGTH_SQUARED || v2.length_squared() <= MIN_EDGE_LENGTH_SQUARED || v3.length_squared() <= MIN_EDGE_LENGTH_SQUARED {
        return (Vector::new(1.0, 0.0, 0.0), 0.0);
    }
    let n = v1.cross(v2) + v2.cross(v3) + v3.cross(v1);
    let len = n.length();
    if len < RESOLUTION {
        (n, 0.0)
    } else {
        let n = n * (1.0 / len);
        (n, n.dot(p1 - Point::new(0.0, 0.0, 0.0)))
    }
}

/// What an `Intf_SectionPoint` says about the polyhedron (`InfoSecond`).
#[derive(Clone, Copy, Debug)]
enum OnPolyhedron {
    Vertex(usize),
    /// An edge between two points (the parameter along it is always 0 here).
    Edge(usize, usize),
    Face(usize),
}

/// An `Intf_SectionPoint`: the segment of the polygon and the parameter along it, and where
/// on the polyhedron.
#[derive(Clone, Copy, Debug)]
struct SectionPoint {
    /// `Pnt()`.
    point: Point,
    segment: usize,
    param: f64,
    on: OnPolyhedron,
}

const POURCENT3: [usize; 4] = [0, 1, 2, 0];

/// `Intf_InterferencePolygonPolyhedron(Polygon, Polyhedron, BoundSortBox)`.
#[doc(alias = "Intf_InterferencePolygonPolyhedron")]
fn interference(polygon: &Polygon, polyhedron: &Polyhedron) -> Vec<SectionPoint> {
    let mut tolerance = polygon.deflection + polyhedron.deflection;
    if tolerance == 0.0 {
        tolerance = epsilon(1000.0);
    }
    let mut out = Vec::new();
    if !polygon.bnd.meets(&polyhedron.bnd) {
        return out;
    }
    // `Interference(thePolyg, thePolyh, PolyhGrid)`.
    let nb_segments = polygon.nb_segments();
    let def_ph = polyhedron.deflection;
    for i_lin in 1..=nb_segments {
        let (p1, p2) = (polygon.pnts[i_lin - 1], polygon.pnts[i_lin]);
        let seg = Bnd::of([p1, p2], polygon.deflection);
        let list: Vec<usize> = polyhedron.tree.overlapping(&BoundingBox::from_points([seg.lo, seg.hi])).into_iter().map(|k| polyhedron.boxed[k]).collect();
        let (mut beg0, mut end0) = (p1, p2);
        for tri in list {
            let [a, b, c] = polyhedron.triangle(tri).map(|i| polyhedron.pnts[i]);
            let (tri_nor, tri_dp) = plane_equation(a, b, c);
            // The boundary segments are lengthened.
            if i_lin == 1 {
                let dif = p1 - p2;
                let dist = dif.length();
                if dist > RESOLUTION {
                    let dif = dif * (1.0 / dist);
                    let cos = dif.dot(tri_nor).abs();
                    if cos > RESOLUTION {
                        beg0 = p1 + dif * (def_ph / cos);
                    }
                }
            } else if i_lin == nb_segments {
                let dif = p2 - p1;
                let dist = dif.length();
                if dist > RESOLUTION {
                    let dif = dif * (1.0 / dist);
                    let cos = dif.dot(tri_nor).abs();
                    if cos > RESOLUTION {
                        end0 = p2 + dif * (def_ph / cos);
                    }
                }
            }
            let origin = Point::new(0.0, 0.0, 0.0);
            let d_beg = tri_nor.dot(beg0 - origin) - tri_dp;
            let d_end = tri_nor.dot(end0 - origin) - tri_dp;
            intersect(beg0, end0, i_lin, tri, polyhedron, tri_nor, d_beg, d_end, tolerance, &mut out);
        }
    }
    out
}

/// `Intersect(BegO, EndO, Infinite = false, TTri, thePolyh, TriNormal, TriDp, dBegTri,
/// dEndTri)`: where the segment crosses the triangle, and where it passes within the
/// tolerance of its sides.
#[allow(clippy::too_many_arguments)]
fn intersect(beg: Point, end: Point, i_lin: usize, tri: usize, polyhedron: &Polyhedron, tri_nor: Vector, d_beg: f64, d_end: f64, tolerance: f64, out: &mut Vec<SectionPoint>) {
    let p_tri = polyhedron.triangle(tri);
    let pt = |k: usize| polyhedron.pnts[p_tri[k]];
    let t = d_beg - d_end;
    let mut param = if !(-1.0e-16..1.0e-16).contains(&t) { d_beg / t } else { d_beg };
    let floatgap = epsilon(1000.0);
    let mut no_intersection = false;
    if (-floatgap..=floatgap).contains(&d_beg) {
        param = 0.0;
    } else if (-floatgap..=floatgap).contains(&d_end) {
        param = 1.0;
    }
    if !(0.0..=1.0).contains(&param) {
        no_intersection = true;
    }
    if !no_intersection {
        let sp_lieu = beg + (end - beg) * param;
        let mut d_pi_e = [0.0f64; 3];
        let mut is = 0;
        let mut s_edge: Option<usize> = None;
        let mut s_vertex: Option<usize> = None;
        let sign = |x: f64| if x > floatgap { 1.0 } else if x < -floatgap { -1.0 } else { 0.0 };
        for k in 0..3 {
            let seg_t = pt((k + 1) % 3) - pt(k);
            let vec_p = sp_lieu - pt(k);
            if vec_p.length() <= floatgap {
                s_vertex = Some(k);
                is = k;
                // (At the first two vertices the search stops; at the last the side is
                // still measured, as in OCCT.)
                if k < 2 {
                    break;
                }
            }
            let cross = seg_t.cross(vec_p);
            d_pi_e[k] = sign(cross.dot(tri_nor)) * (cross.length() / seg_t.length());
            if (-floatgap..=floatgap).contains(&d_pi_e[k]) {
                s_edge = Some(k);
                is = k;
                if k < 2 {
                    break;
                }
            }
        }
        if let Some(_v) = s_vertex {
            out.push(SectionPoint { point: sp_lieu, segment: i_lin, param, on: OnPolyhedron::Vertex(p_tri[is]) });
        } else if let Some(e) = s_edge {
            let (a, b) = (p_tri[e], p_tri[POURCENT3[e + 1]]);
            out.push(SectionPoint { point: sp_lieu, segment: i_lin, param, on: OnPolyhedron::Edge(a.min(b), a.max(b)) });
        } else if d_pi_e[0] > 0.0 && d_pi_e[1] > 0.0 && d_pi_e[2] > 0.0 {
            out.push(SectionPoint { point: sp_lieu, segment: i_lin, param, on: OnPolyhedron::Face(tri) });
        } else {
            // The triangulation may not cover the whole face: the deflection between the
            // boundary isolines and the boundary triangles is taken into account.
            for i in 1..=3 {
                let ind_p1 = if i == 3 { p_tri[0] } else { p_tri[i] };
                let ind_p2 = p_tri[i - 1];
                if polyhedron.is_on_bound(ind_p1, ind_p2) {
                    let beg_p = polyhedron.pnts[ind_p1];
                    let end_p = polyhedron.pnts[ind_p2];
                    let vec_tri = end_p - beg_p;
                    let dir_tri = vec_tri * (1.0 / vec_tri.length());
                    let w = sp_lieu - beg_p;
                    let dist = (w - dir_tri * w.dot(dir_tri)).length();
                    if dist <= polyhedron.border_deflection {
                        let par = w.dot(dir_tri);
                        if par >= 0.0 && par <= vec_tri.length() {
                            out.push(SectionPoint { point: sp_lieu, segment: i_lin, param, on: OnPolyhedron::Face(tri) });
                        }
                    }
                }
            }
        }
    }
    // The distance between the triangle's sides and the polygon: within the tolerance, a
    // section point.
    let vec_pol = end - beg;
    let n_vec_pol = vec_pol.length();
    if n_vec_pol < COMPUTATIONAL {
        return;
    }
    let dir_pol = vec_pol * (1.0 / n_vec_pol);
    for i in 0..3 {
        let p_ip1 = p_tri[POURCENT3[i + 1]];
        let p_i = p_tri[i];
        let beg_t = polyhedron.pnts[p_ip1];
        let end_t = polyhedron.pnts[p_i];
        let vec_tri = end_t - beg_t;
        let n_vec_tri = vec_tri.length();
        if n_vec_tri < COMPUTATIONAL {
            continue;
        }
        let dir_tri = vec_tri * (1.0 / n_vec_tri);
        let Some((po, ptt, dist2)) = ext_lin_lin(beg, dir_pol, beg_t, dir_tri) else { continue };
        if dist2 <= tolerance * tolerance {
            if let Some(param_on_o) = in_segment(vec_pol, po - beg, n_vec_pol, tolerance) {
                if let Some(param_on_t) = in_segment(vec_tri, ptt - beg_t, n_vec_tri, tolerance) {
                    let sp_lieu = beg_t + (end_t - beg_t) * param_on_t;
                    out.push(SectionPoint { point: sp_lieu, segment: i_lin, param: param_on_o, on: OnPolyhedron::Edge(p_i.min(p_ip1), p_i.max(p_ip1)) });
                }
            }
        }
    }
}

/// `IsInSegment`: the parameter (in `[0, 1]`) along the segment of the projection of a point
/// within the tolerance of it.
fn in_segment(p1p2: Vector, p1p: Vector, np1p2: f64, tolerance: f64) -> Option<f64> {
    if np1p2 < COMPUTATIONAL {
        return None;
    }
    let mut param = p1p2.dot(p1p) / np1p2;
    if param > np1p2 + tolerance || param < -tolerance {
        return None;
    }
    param /= np1p2;
    Some(param.clamp(0.0, 1.0))
}

/// `Extrema_ExtElC(gp_Lin, gp_Lin, AngTol)` for lines that are not parallel: the nearest
/// points and their squared distance; `None` for parallel lines.
#[doc(alias = "Extrema_ExtElC")]
fn ext_lin_lin(l1: Point, d1: Vector, l2: Point, d2: Vector) -> Option<(Point, Point, f64)> {
    let cos_a = d1.dot(d2);
    let sq_sin_a = 1.0 - cos_a * cos_a;
    // (`IsParallel(D2, Precision::Angular())`: the angle within it, either way.)
    let angle = d1.cross(d2).length().atan2(cos_a);
    if sq_sin_a < RESOLUTION || angle <= ANGULAR || std::f64::consts::PI - angle <= ANGULAR {
        return None;
    }
    let l1l2 = l2 - l1;
    let (d1l, d2l) = (d1.dot(l1l2), d2.dot(l1l2));
    let u1 = (d1l - cos_a * d2l) / sq_sin_a;
    let u2 = (cos_a * d1l - d2l) / sq_sin_a;
    if u1.abs() >= 0.5 * INFINITE || u2.abs() >= 0.5 * INFINITE {
        return None;
    }
    let (p1, p2) = (l1 + d1 * u1, l2 + d2 * u2);
    Some((p1, p2, p1.distance_squared(p2)))
}

/// `SectionPointToParameters`: the section point's approximate parameters `(u, v, w)`.
fn section_point_parameters(sp: &SectionPoint, polyhedron: &Polyhedron, polygon: &Polygon) -> (f64, f64, f64) {
    let (u1, v1) = match sp.on {
        OnPolyhedron::Vertex(a) => (polyhedron.u[a], polyhedron.v[a]),
        // (`u1 += Param * (u - u1)` with the parameter 0.)
        OnPolyhedron::Edge(a, _) => (polyhedron.u[a], polyhedron.v[a]),
        OnPolyhedron::Face(tri) => {
            let [pt1, pt2, pt3] = polyhedron.triangle(tri);
            let (pa, pb, pc) = (polyhedron.pnts[pt1], polyhedron.pnts[pt2], polyhedron.pnts[pt3]);
            let (ua, va, ub, vb, uc, vc) = (polyhedron.u[pt1], polyhedron.v[pt1], polyhedron.u[pt2], polyhedron.v[pt2], polyhedron.u[pt3], polyhedron.v[pt3]);
            let p = sp.point;
            let normale = (pb - pa).cross(pc - pa);
            let mut cc = (pb - pa).cross(p - pa).dot(normale);
            let mut ca = (pc - pb).cross(p - pb).dot(normale);
            let mut cb = (pa - pc).cross(p - pc).dot(normale);
            let cabc = ca + cb + cc;
            if cabc.abs() > COMPUTATIONAL {
                ca /= cabc;
                cb /= cabc;
                cc /= cabc;
                (ca * ua + cb * ub + cc * uc, ca * va + cb * vb + cc * vc)
            } else {
                let (ab, bc, ca_) = (pb - pa, pc - pb, pa - pc);
                let (sq_ab, sq_bc, sq_ca) = (ab.length_squared(), bc.length_squared(), ca_.length_squared());
                if sq_ab >= sq_bc && sq_ab >= sq_ca && sq_ab > SQUARE_CONFUSION {
                    let t = ((p - pa).dot(ab) / sq_ab).clamp(0.0, 1.0);
                    (ua + t * (ub - ua), va + t * (vb - va))
                } else if sq_bc >= sq_ca && sq_bc > SQUARE_CONFUSION {
                    let t = ((p - pb).dot(bc) / sq_bc).clamp(0.0, 1.0);
                    (ub + t * (uc - ub), vb + t * (vc - vb))
                } else if sq_ca > SQUARE_CONFUSION {
                    let t = ((p - pc).dot(ca_) / sq_ca).clamp(0.0, 1.0);
                    (uc + t * (ua - uc), vc + t * (va - vc))
                } else {
                    (ua, va)
                }
            }
        }
    };
    let w = polygon.approx_param_on_curve(sp.segment, sp.param);
    (u1, v1, w)
}

/// `SortStartPoints`: by `w`, then `u` for the same `w`, then `v` for the same `w` and `u`
/// (bubble sorts, near-equal values made equal).
fn sort_start_points(points: &mut [(f64, f64, f64)]) {
    let n = points.len();
    if n == 0 {
        return;
    }
    let ptol = 10.0 * P_CONFUSION;
    loop {
        let mut ok = true;
        for i in 1..n {
            if points[i].2 < points[i - 1].2 {
                points.swap(i, i - 1);
                ok = false;
            }
        }
        if ok {
            break;
        }
    }
    loop {
        let mut ok = true;
        for i in 1..n {
            if points[i].2 - points[i - 1].2 < ptol {
                points[i].2 = points[i - 1].2;
                if points[i].0 < points[i - 1].0 {
                    let (a, b) = (points[i], points[i - 1]);
                    points[i].0 = b.0;
                    points[i].1 = b.1;
                    points[i - 1].0 = a.0;
                    points[i - 1].1 = a.1;
                    ok = false;
                }
            }
        }
        if ok {
            break;
        }
    }
    loop {
        let mut ok = true;
        for i in 1..n {
            if points[i].2 - points[i - 1].2 < ptol && points[i].0 - points[i - 1].0 < ptol {
                points[i].0 = points[i - 1].0;
                if points[i].1 < points[i - 1].1 {
                    let (a, b) = (points[i].1, points[i - 1].1);
                    points[i].1 = b;
                    points[i - 1].1 = a;
                    ok = false;
                }
            }
        }
        if ok {
            break;
        }
    }
}

/// `IntImp_ZerCSParFunc`: `S(u, v) - L(w)`, remembering the square of its last value.
struct CsFunction<'a> {
    surface: &'a Surface,
    line: &'a LineCurve,
    f: f64,
}

impl FunctionSetWithDerivatives for CsFunction<'_> {
    fn nb_variables(&self) -> usize {
        3
    }

    fn nb_equations(&self) -> usize {
        3
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool {
        let (ps, du, dv) = self.surface.d1(x[0], x[1]);
        let pc = self.line.value(x[2]);
        let dw = self.line.direction;
        d.copy_from_slice(&[du.x, dv.x, -dw.x, du.y, dv.y, -dw.y, du.z, dv.z, -dw.z]);
        f[0] = ps.x - pc.x;
        f[1] = ps.y - pc.y;
        f[2] = ps.z - pc.z;
        self.f = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        true
    }
}

/// `IntImp_IntCS::Perform(U, V, W, Rsnld, u0, u1, v0, v1, w0, w1)`: the solution `(w, u, v)`
/// from the start, else from the start moved to either end of the line's range.
fn exact_perform(function: &mut CsFunction, rsnld: &mut FunctionSetRoot, tol: f64, start: [f64; 3], inf: [f64; 3], sup: [f64; 3], bounds: (f64, f64, f64, f64)) -> Option<(f64, f64, f64)> {
    let mut uvap = start;
    let tolerance = [u_resolution(function.surface, bounds, CONFUSION), v_resolution(function.surface, bounds, CONFUSION), CONFUSION];
    rsnld.set_tolerance(&tolerance);
    for attempt in 0..3 {
        if attempt == 1 {
            uvap[2] = inf[2];
        } else if attempt == 2 {
            uvap[2] = sup[2];
        }
        rsnld.perform(function, &uvap, &inf, &sup);
        if rsnld.is_done() && function.f.abs() <= tol {
            let root = rsnld.root();
            return Some((root[2], root[0], root[1]));
        }
    }
    None
}

#[cfg(test)]
mod tests;
