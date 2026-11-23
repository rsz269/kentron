//! Where a curve meets an infinite plane: every transversal crossing and tangential contact
//! over the curve's range, with the curve parameter, the plane parameters, the point and the
//! direction of the crossing relative to the plane's normal.
//!
//! Lines and conics are solved in closed form by substituting their parametrisation into the
//! plane's equation `N · (q − O) = 0`: a line gives a linear equation, a circle or an ellipse
//! `A cos t + B sin t + C = 0`, a parabola a quadratic in `t`, a hyperbola a quadratic in
//! `eᵗ`. Any other curve is searched numerically: the signed distance to the plane is
//! sampled along each C1 span, sign changes are refined by Brent's method, and dips of the
//! distance between samples (seen from the samples' derivatives) are examined for touching
//! points and for pairs of crossings closer together than the sampling.
//!
//! The module also holds a few small analytic records shared with the line–surface
//! intersection: a frame as coordinate arrays ([`Ax3`]), an implicit quadric's coefficients
//! ([`Quadric`]), a conic–quadric result ([`ConicQuad`]) and the line–plane solver
//! ([`line_plane`]).

use std::f64::consts::PI;

use crate::geometry::adaptor::{in_period, CurveAdaptor, CurveKind};
use crate::geometry::curve::Curve;
use crate::geometry::function_roots::{bracketed_root, Function, NumericError};
use crate::geometry::surface::Plane;
use crate::geometry::Continuity;
use crate::math::{Direction, Frame, Point, Vector};

/// The angular tolerance of the analytic intersections (and of the transition test).
pub(super) const TOLERANCE_ANGULAIRE: f64 = 1.0e-12;

/// A solution this far outside the curve's range is dropped; nearer ones are kept.
const RANGE_TOLERANCE: f64 = 1e-8;

/// Two solutions this close in every parameter (and of the same transition) are one.
const DUPLICATE: f64 = 1e-8;

/// A distance from the plane up to this counts as a contact (and stretches of it as lying in
/// the plane).
const CONTACT: f64 = 1e-8;

/// How a curve passes the plane at a meeting, relative to the plane's normal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Transition {
    /// The curve touches the plane, or runs along it.
    Tangent,
    /// The curve goes against the normal (from the side the normal points to).
    In,
    /// The curve goes along the normal.
    Out,
}

/// One meeting of a curve with a plane.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CurvePlanePoint {
    /// The curve point at `w`.
    pub point: Point,
    /// The meeting's coordinates in the plane (along its X and Y axes from its origin).
    pub u: f64,
    pub v: f64,
    /// The curve parameter.
    pub w: f64,
    pub transition: Transition,
}

/// All meetings of a curve with a plane, by ascending curve parameter.
#[derive(Clone, Debug, Default)]
pub(crate) struct CurvePlaneIntersection {
    pub points: Vec<CurvePlanePoint>,
    /// The curve (a line or a conic) is parallel to the plane, or lies in it.
    pub parallel: bool,
}

/// A frame as coordinate arrays: origin and three unit axes.
#[derive(Clone, Copy, Debug)]
pub(super) struct Ax3 {
    pub(super) location: [f64; 3],
    pub(super) x: [f64; 3],
    pub(super) y: [f64; 3],
    pub(super) z: [f64; 3],
}

fn coordinates(v: Vector) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl Ax3 {
    pub(super) fn from_frame(frame: &Frame) -> Ax3 {
        let o = frame.origin();
        Ax3 {
            location: [o.x, o.y, o.z],
            x: coordinates(Vector::from(frame.x())),
            y: coordinates(Vector::from(frame.y())),
            z: coordinates(Vector::from(frame.z())),
        }
    }

    /// `(m, t)` taking a global point `p` to its coordinates `m·p + t` in the frame: the rows
    /// of `m` are the axes, and `t = −m·location`.
    pub(super) fn transformation(&self) -> ([[f64; 3]; 3], [f64; 3]) {
        let m = [self.x, self.y, self.z];
        let t = [-dot3(m[0], self.location), -dot3(m[1], self.location), -dot3(m[2], self.location)];
        (m, t)
    }
}

/// The quadric `cxx x² + cyy y² + czz z² + 2(cxy xy + cxz xz + cyz yz) + 2(cx x + cy y + cz z)
/// + ccte = 0`.
#[derive(Clone, Copy, Debug)]
pub(super) struct Quadric {
    pub(super) cxx: f64,
    pub(super) cyy: f64,
    pub(super) czz: f64,
    pub(super) cxy: f64,
    pub(super) cxz: f64,
    pub(super) cyz: f64,
    pub(super) cx: f64,
    pub(super) cy: f64,
    pub(super) cz: f64,
    pub(super) ccte: f64,
}

/// The meeting of a conic (or line) with a quadric (or plane).
#[derive(Clone, Debug, Default)]
pub(super) struct ConicQuad {
    /// The computation completed.
    pub(super) done: bool,
    /// The conic is parallel to the surface: no isolated points.
    pub(super) parallel: bool,
    /// The conic lies in the surface.
    pub(super) inquadric: bool,
    /// Isolated meetings: coordinates and the parameter on the conic.
    pub(super) points: Vec<([f64; 3], f64)>,
}

/// The line `origin + t·direction` against the plane through `plane.location` with normal
/// `plane.z`. Parallel when the direction's normal component is below `tolang`, unless
/// (`len ≠ 0`) going `len` along the line from the foot of `origin` leaves the plane by more
/// than `tol`; then `inquadric` when `origin` is within `tolang` of the plane.
pub(super) fn line_plane(origin: [f64; 3], direction: [f64; 3], plane: &Ax3, tolang: f64, tol: f64, len: f64) -> ConicQuad {
    let nz = plane.z;
    let norm = dot3(nz, nz).sqrt();
    let n = [nz[0] / norm, nz[1] / norm, nz[2] / norm];
    let delta = dot3(n, direction);
    let rel = [origin[0] - plane.location[0], origin[1] - plane.location[1], origin[2] - plane.location[2]];
    let s = dot3(n, rel);
    let mut parallel = delta.abs() < tolang;
    if parallel && len != 0.0 && delta != 0.0 {
        // The foot of `origin`, moved `len` along the line: its distance from the plane.
        let dlen = dot3(direction, direction).sqrt();
        let moved = (len / dlen) * delta;
        if moved.abs() > tol {
            parallel = false;
        }
    }
    if parallel {
        return ConicQuad { done: true, parallel: true, inquadric: s.abs() < tolang, points: Vec::new() };
    }
    let t = -s / delta;
    let p = [origin[0] + t * direction[0], origin[1] + t * direction[1], origin[2] + t * direction[2]];
    ConicQuad { done: true, parallel: false, inquadric: false, points: vec![(p, t)] }
}

/// The plane through `p` with unit normal `v` and the conventional default X axis (that of
/// [`Frame::from_z`]: the smallest component of `v` dropped, the other two swapped with one
/// sign change).
pub(crate) fn plane_through(p: Point, v: Direction) -> Plane {
    let (a, b, c) = (v.x(), v.y(), v.z());
    let (aa, ab, ac) = (a.abs(), b.abs(), c.abs());
    let hint = if ab <= aa && ab <= ac {
        if aa > ac { (-c, 0.0, a) } else { (c, 0.0, -a) }
    } else if aa <= ab && aa <= ac {
        if ab > ac { (0.0, -c, b) } else { (0.0, c, -b) }
    } else if aa > ab {
        (-b, a, 0.0)
    } else {
        (b, -a, 0.0)
    };
    let x = Direction::new(hint.0, hint.1, hint.2);
    let y = Direction::from_vector(Vector::from(v).cross(Vector::from(x)));
    Plane::from_frame(Frame::from_axes_unchecked(p, x, y, v))
}

/// A closure as a [`Function`].
struct Closure<F: FnMut(f64) -> f64>(F);

impl<F: FnMut(f64) -> f64> Function for Closure<F> {
    fn value(&mut self, x: f64, value: &mut f64) -> bool {
        *value = (self.0)(x);
        value.is_finite()
    }
}

/// The plane's normal and origin as plain vectors.
struct PlaneData {
    origin: Point,
    x: Vector,
    y: Vector,
    n: Vector,
}

impl PlaneData {
    fn new(plane: &Plane) -> PlaneData {
        let f = plane.frame();
        let (x, y) = (Vector::from(f.x()), Vector::from(f.y()));
        PlaneData { origin: f.origin(), x, y, n: x.cross(y) }
    }

    fn distance(&self, q: Point) -> f64 {
        self.n.dot(q - self.origin)
    }
}

/// Every meeting of `curve` (over its range) with `plane`, by ascending curve parameter.
pub(crate) fn intersect_curve_plane(curve: &dyn CurveAdaptor, plane: &Plane) -> Result<CurvePlaneIntersection, NumericError> {
    let pd = PlaneData::new(plane);
    let (w0, w1) = (curve.first_parameter(), curve.last_parameter());
    let kind = curve.curve_type();
    // Solutions as (curve parameter, meeting point).
    let mut parallel = false;
    let solutions: Vec<(f64, Point)> = match (kind, curve.basis_curve()) {
        (Some(CurveKind::Line), Some(Curve::Line(l))) => {
            let d = Vector::from(l.direction);
            let delta = pd.n.dot(d);
            if delta.abs() < TOLERANCE_ANGULAIRE {
                parallel = true;
                Vec::new()
            } else {
                let t = -pd.distance(l.origin) / delta;
                vec![(t, l.origin + d * t)]
            }
        }
        (Some(CurveKind::Circle), Some(Curve::Circle(c))) => {
            let f = c.frame();
            let r = c.radius();
            if Vector::from(f.z()).cross(pd.n).length() <= TOLERANCE_ANGULAIRE {
                parallel = true;
                Vec::new()
            } else {
                trig_solutions(&pd, f, r, r, true).into_iter().map(|t| (t, ellipse_point(f, r, r, t))).collect()
            }
        }
        (Some(CurveKind::Ellipse), Some(Curve::Ellipse(e))) => {
            let f = e.frame();
            let (a, b) = (e.major_radius(), e.minor_radius());
            let (ca, cb, cc) = (a * pd.n.dot(Vector::from(f.x())), b * pd.n.dot(Vector::from(f.y())), pd.distance(f.origin()));
            if ca.hypot(cb) <= TOLERANCE_ANGULAIRE * a.max(b).max(1.0) {
                // (In a plane parallel to the cutting one: in it, or no meeting at all.)
                parallel = cc.abs() <= TOLERANCE_ANGULAIRE * a.max(1.0);
                Vec::new()
            } else {
                trig_solutions(&pd, f, a, b, false).into_iter().map(|t| (t, ellipse_point(f, a, b, t))).collect()
            }
        }
        (Some(CurveKind::Parabola), Some(Curve::Parabola(p))) => {
            let f = p.frame();
            let focal = p.focal();
            let (nx, ny) = (pd.n.dot(Vector::from(f.x())), pd.n.dot(Vector::from(f.y())));
            let cc = pd.distance(f.origin());
            let point = |t: f64| f.origin() + Vector::from(f.x()) * (t * t / (4.0 * focal)) + Vector::from(f.y()) * t;
            if nx.abs() <= TOLERANCE_ANGULAIRE && ny.abs() <= TOLERANCE_ANGULAIRE {
                parallel = cc.abs() <= TOLERANCE_ANGULAIRE * focal.max(1.0);
                Vec::new()
            } else {
                // nx t²/(4f) + ny t + c = 0.
                let g = |t: f64| nx * t * t / (4.0 * focal) + ny * t + cc;
                let dg = |t: f64| nx * t / (2.0 * focal) + ny;
                let roots = crate::geometry::polynomial_roots::PolynomialRoots::quadratic(nx / (4.0 * focal), ny, cc);
                let mut ts: Vec<f64> = roots.roots.iter().map(|&t| polish(&g, &dg, t)).collect();
                dedup_close(&mut ts);
                ts.into_iter().map(|t| (t, point(t))).collect()
            }
        }
        (Some(CurveKind::Hyperbola), Some(Curve::Hyperbola(h))) => {
            let f = h.frame();
            let (a, b) = (h.major_radius(), h.minor_radius());
            let (ca, cb) = (a * pd.n.dot(Vector::from(f.x())), b * pd.n.dot(Vector::from(f.y())));
            let cc = pd.distance(f.origin());
            let point = |t: f64| f.origin() + Vector::from(f.x()) * (a * t.cosh()) + Vector::from(f.y()) * (b * t.sinh());
            if ca.abs() <= TOLERANCE_ANGULAIRE * a.max(1.0) && cb.abs() <= TOLERANCE_ANGULAIRE * b.max(1.0) {
                parallel = cc.abs() <= TOLERANCE_ANGULAIRE * a.max(1.0);
                Vec::new()
            } else {
                // A cosh t + B sinh t + C = 0; with v = eᵗ: (A + B) v² + 2C v + (A − B) = 0.
                let g = |t: f64| ca * t.cosh() + cb * t.sinh() + cc;
                let dg = |t: f64| ca * t.sinh() + cb * t.cosh();
                let roots = crate::geometry::polynomial_roots::PolynomialRoots::quadratic(ca + cb, 2.0 * cc, ca - cb);
                let mut ts: Vec<f64> = roots.roots.iter().filter(|&&v| v > 0.0 && v.is_finite()).map(|&v| polish(&g, &dg, v.ln())).collect();
                dedup_close(&mut ts);
                ts.into_iter().map(|t| (t, point(t))).collect()
            }
        }
        _ => numerical_solutions(curve, &pd)?.into_iter().map(|w| (w, curve.value(w))).collect(),
    };

    let periodic = curve.is_periodic() || matches!(kind, Some(CurveKind::Circle | CurveKind::Ellipse));
    let mut accepted: Vec<CurvePlanePoint> = Vec::new();
    let mut decorated: Vec<CurvePlanePoint> = Vec::new();
    for (w, q) in solutions {
        let w = if periodic { in_period(w, w0, w0 + curve.period()) } else { w };
        if w0 - w >= RANGE_TOLERANCE || w - w1 >= RANGE_TOLERANCE {
            continue;
        }
        let rel = q - pd.origin;
        let (u, v) = (rel.dot(pd.x), rel.dot(pd.y));
        let point = curve.value(w);
        let d = curve.d1(w).1;
        let transition = transition_of(pd.n, d);
        decorated.push(CurvePlanePoint { point, u, v, w, transition });
    }
    decorated.sort_by(|a, b| a.w.total_cmp(&b.w));
    for p in decorated {
        let duplicate = accepted.iter().any(|a| {
            (a.u - p.u).abs() < DUPLICATE && (a.v - p.v).abs() < DUPLICATE && (a.w - p.w).abs() < DUPLICATE && a.transition == p.transition
        });
        if !duplicate {
            accepted.push(p);
        }
    }
    Ok(CurvePlaneIntersection { points: accepted, parallel })
}

/// The transition of a curve going along `d` through a plane of normal `n`.
fn transition_of(n: Vector, d: Vector) -> Transition {
    let nn = n.length();
    let dd = d.length_squared();
    if nn <= TOLERANCE_ANGULAIRE || dd <= TOLERANCE_ANGULAIRE {
        return Transition::Tangent;
    }
    let c = n.dot(d) / dd.sqrt() / nn;
    if c < -TOLERANCE_ANGULAIRE {
        Transition::In
    } else if c > TOLERANCE_ANGULAIRE {
        Transition::Out
    } else {
        Transition::Tangent
    }
}

/// `c + a cos t X + b sin t Y` in the frame `f`.
fn ellipse_point(f: &Frame, a: f64, b: f64, t: f64) -> Point {
    let (s, c) = t.sin_cos();
    f.origin() + Vector::from(f.x()) * (a * c) + Vector::from(f.y()) * (b * s)
}

/// Newton steps on `g` from `t`, kept while they lower `|g|`.
fn polish(g: &dyn Fn(f64) -> f64, dg: &dyn Fn(f64) -> f64, mut t: f64) -> f64 {
    let mut gt = g(t);
    for _ in 0..8 {
        let d = dg(t);
        if gt == 0.0 || d == 0.0 || !d.is_finite() {
            break;
        }
        let next = t - gt / d;
        let gn = g(next);
        if !(gn.abs() < gt.abs()) {
            break;
        }
        t = next;
        gt = gn;
    }
    t
}

/// Sorts `ts` and merges values closer than `1e-12` (relative).
fn dedup_close(ts: &mut Vec<f64>) {
    ts.sort_by(f64::total_cmp);
    ts.dedup_by(|b, a| (*b - *a).abs() <= 1e-12 * a.abs().max(1.0));
}

/// The angles where the ellipse (or circle, `a = b`) of frame `f` meets the plane:
/// `A cos t + B sin t + C = 0` with `A = a N·X`, `B = b N·Y`, `C` the centre's distance. A
/// circle tangent to the plane to within a few units in the last place of its radius gives
/// one angle.
fn trig_solutions(pd: &PlaneData, f: &Frame, a: f64, b: f64, circle: bool) -> Vec<f64> {
    let ca = a * pd.n.dot(Vector::from(f.x()));
    let cb = b * pd.n.dot(Vector::from(f.y()));
    let cc = pd.distance(f.origin());
    let r = ca.hypot(cb);
    let phi = cb.atan2(ca);
    let k = -cc / r;
    // (`|k|` is the centre's distance from the common line over the radius, for a circle.)
    let slack = if circle { 4.0 * f64::EPSILON } else { 8.0 * f64::EPSILON };
    if !k.is_finite() || k.abs() > 1.0 + slack {
        return Vec::new();
    }
    let g = |t: f64| {
        let (s, c) = t.sin_cos();
        ca * c + cb * s + cc
    };
    let dg = |t: f64| {
        let (s, c) = t.sin_cos();
        -ca * s + cb * c
    };
    if k.abs() >= 1.0 - slack {
        let t = if k > 0.0 { phi } else { phi + PI };
        return vec![t];
    }
    let alpha = k.acos();
    let mut ts = vec![polish(&g, &dg, phi - alpha), polish(&g, &dg, phi + alpha)];
    dedup_close(&mut ts);
    ts
}

/// The roots of the distance to the plane along a general curve: crossings, touching
/// points, ascending (duplicates at span ends possible; the caller removes them).
fn numerical_solutions(curve: &dyn CurveAdaptor, pd: &PlaneData) -> Result<Vec<f64>, NumericError> {
    let (w0, w1) = (curve.first_parameter(), curve.last_parameter());
    if !(w1 >= w0) {
        return Ok(Vec::new());
    }
    if w1 == w0 {
        let d = pd.distance(curve.value(w0));
        return Ok(if d.abs() <= CONTACT { vec![w0] } else { Vec::new() });
    }
    let breaks = curve.intervals(Continuity::C1);
    let per_span = samples_per_span(curve, breaks.len().saturating_sub(1).max(1));
    let mut out = Vec::new();
    // Spans lying in the plane (their ends are not reported for themselves).
    let mut flat: Vec<(f64, f64)> = Vec::new();
    for pair in breaks.windows(2) {
        let (a, b) = (pair[0].max(w0), pair[1].min(w1));
        if !(b > a) {
            continue;
        }
        if !span_solutions(curve, pd, a, b, per_span, &mut out)? {
            flat.push((a, b));
        }
    }
    // A range end within the contact distance whose Newton step stays within it: a crossing
    // or contact just at (or a rounding step beyond) the end.
    for end in [w0, w1] {
        let (c, d) = curve.d1(end);
        let dist = pd.distance(c);
        let slope = pd.n.dot(d);
        if flat.iter().any(|&(a, b)| a == end || b == end) {
            continue;
        }
        if dist.abs() <= CONTACT && (dist.abs() <= 1e-12 || (dist / slope).abs() <= RANGE_TOLERANCE) && !out.iter().any(|&w: &f64| (w - end).abs() < DUPLICATE) {
            out.push(end);
        }
    }
    out.sort_by(f64::total_cmp);
    Ok(out)
}

/// Samples per span: enough for the oscillation of a polynomial piece of the curve's degree,
/// with a modest cap overall.
fn samples_per_span(curve: &dyn CurveAdaptor, spans: usize) -> usize {
    let degree = match curve.basis_curve() {
        Some(Curve::BSpline(b)) => b.degree(),
        Some(Curve::Bezier(b)) => b.degree(),
        _ => 6,
    };
    let wanted = (4 * (degree + 1)).max(16);
    let cap = (4000 / spans).max(8);
    wanted.min(cap)
}

/// The roots of the distance on the span `[a, b]` sampled at `n + 1` points; `false` (and
/// none) when the span lies in the plane.
fn span_solutions(curve: &dyn CurveAdaptor, pd: &PlaneData, a: f64, b: f64, n: usize, out: &mut Vec<f64>) -> Result<bool, NumericError> {
    let xs: Vec<f64> = (0..=n).map(|i| if i == n { b } else { a + (b - a) * i as f64 / n as f64 }).collect();
    let mut fs = Vec::with_capacity(n + 1);
    let mut ds = Vec::with_capacity(n + 1);
    for &x in &xs {
        let (c, d) = curve.d1(x);
        let (f, s) = (pd.distance(c), pd.n.dot(d));
        if !f.is_finite() || !s.is_finite() {
            return Err(NumericError);
        }
        fs.push(f);
        ds.push(s);
    }
    // (A span lying in the plane gives no points.)
    if fs.iter().all(|f| f.abs() <= CONTACT) {
        return Ok(false);
    }
    let dist = |w: f64| pd.distance(curve.value(w));
    let slope = |w: f64| pd.n.dot(curve.d1(w).1);
    for i in 0..=n {
        if fs[i] == 0.0 {
            out.push(xs[i]);
        }
    }
    let mut dips: Vec<(f64, f64, f64)> = Vec::new();
    for i in 0..n {
        let (x0, x1, f0, f1) = (xs[i], xs[i + 1], fs[i], fs[i + 1]);
        if f0 == 0.0 || f1 == 0.0 {
            continue;
        }
        if (f0 < 0.0) != (f1 < 0.0) {
            out.push(root(&dist, x0, x1)?);
            continue;
        }
        // Same sign at both ends: does the distance dip towards the plane in between (both
        // ends' slopes heading towards it)?
        let sign = f0.signum();
        if sign * ds[i] < 0.0 && sign * ds[i + 1] > 0.0 {
            dips.push((x0, x1, sign));
        }
    }
    // Samples nearest the plane among their neighbours (of one sign): the dip is around them.
    for i in 0..=n {
        let f = fs[i];
        if f == 0.0 {
            continue;
        }
        let sign = f.signum();
        let lower = |j: usize| fs[j] != 0.0 && fs[j].signum() == sign && sign * fs[j] >= sign * f;
        let left = i == 0 || lower(i - 1);
        let right = i == n || lower(i + 1);
        if !(left && right) {
            continue;
        }
        let (lo, hi) = (xs[i.saturating_sub(1)], xs[(i + 1).min(n)]);
        if (i == 0 || i == n) && f.abs() > CONTACT {
            continue;
        }
        // (Dips found twice give the same points, removed later.)
        dips.retain(|&(a, b, _)| !(a >= lo && b <= hi));
        dips.push((lo, hi, sign));
        if (i == 0 || i == n) && f.abs() <= CONTACT {
            // A contact at the span's end unless the slope has a root next to it.
            let (s0, s1) = (slope(lo), slope(hi));
            if (s0 < 0.0) == (s1 < 0.0) && s0 != 0.0 && s1 != 0.0 {
                out.push(xs[i]);
            }
        }
    }
    for (lo, hi, sign) in dips {
        // The stationary point of the distance in the bracket where the slope changes sign.
        let (s0, s1) = (slope(lo), slope(hi));
        if s0 != 0.0 && s1 != 0.0 && (s0 < 0.0) == (s1 < 0.0) {
            continue;
        }
        let m = root(&slope, lo, hi)?;
        let fm = dist(m);
        if fm == 0.0 || (fm.abs() <= CONTACT && fm.signum() == sign) {
            // A touching point.
            out.push(m);
        } else if fm.signum() != sign {
            // Two crossings closer than the sampling.
            let (flo, fhi) = (dist(lo), dist(hi));
            if flo.signum() == sign {
                out.push(root(&dist, lo, m)?);
            }
            if fhi.signum() == sign {
                out.push(root(&dist, m, hi)?);
            }
        }
    }
    Ok(true)
}

/// A root of `f` bracketed by `[a, b]` to full precision.
fn root(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> Result<f64, NumericError> {
    let mut function = Closure(|x: f64| f(x));
    let r = bracketed_root(&mut function, a, b, 0.0, 200, f64::EPSILON);
    if r.done { Ok(r.root) } else if r.value.is_finite() { Ok(r.root) } else { Err(NumericError) }
}

#[cfg(test)]
pub(crate) mod tests;
