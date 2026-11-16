//! OCCT's Frenet trihedron (`GeomFill_Frenet`), ported: the tangent, principal normal and
//! binormal of a curve, with the points where its curvature vanishes (`D1 × D2` nought) found
//! as OCCT finds them (`Init`, by the extrema of `GeomFill_SnglrFunc` from the origin) and
//! the trihedron there taken from higher derivatives (`DoSingular`, `SingularD0`) or from a
//! point nearby (`D0` where the first derivative vanishes).

use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::ext_pc::ext_pc_other;
use crate::geometry::Continuity;
use crate::math::{Point, Vector};

/// `Precision::Confusion()`, `Precision::PConfusion()`, `Precision::Angular()`.
pub(super) const CONFUSION: f64 = 1e-7;
pub(super) const P_CONFUSION: f64 = 1e-9;
pub(super) const ANGULAR: f64 = 1e-12;
/// `gp::Resolution()`.
pub(super) const RESOLUTION: f64 = f64::MIN_POSITIVE;

/// The file's constants (`NullTol`, `MaxSingular`, `maxDerivOrder`).
const NULL_TOL: f64 = 1e-10;
const MAX_SINGULAR: f64 = 1e-5;
const MAX_DERIV_ORDER: u32 = 3;

/// `gp_Dir(v)`: `v` divided by its length, unless its square length is 1 within
/// `gp::Resolution` (kept as it is then).
pub(super) fn direction(v: Vector) -> Vector {
    let square = v.x * v.x + v.y * v.y + v.z * v.z;
    if (1.0 - RESOLUTION..=1.0 + RESOLUTION).contains(&square) {
        return v;
    }
    let d = square.sqrt();
    Vector::new(v.x / d, v.y / d, v.z / d)
}

/// `gp_Vec::Angle` (by way of `gp_Dir::Angle`, with its constants): the angle between two
/// vectors, in `[0, π]`.
#[allow(clippy::approx_constant)]
pub(super) fn angle(a: Vector, b: Vector) -> f64 {
    // (OCCT's constant, not quite 1/√2.)
    let (a, b) = (direction(a), direction(b));
    let cos = a.dot(b);
    if cos > -0.70710678118655 && cos < 0.70710678118655 {
        cos.acos()
    } else {
        let sin = a.cross(b).length();
        if cos < 0.0 { std::f64::consts::PI - sin.asin() } else { sin.asin() }
    }
}

/// `gp_Ax2(P, V).YDirection()`: the Y direction OCCT gives an axis system about `v`.
pub(super) fn default_y(v: Vector) -> Vector {
    let (a, b, c) = (v.x.abs(), v.y.abs(), v.z.abs());
    // (gp_Ax2's X direction: the axis of the smallest component dropped.)
    let x = if b <= a && b <= c {
        if a > c { Vector::new(-v.z, 0.0, v.x) } else { Vector::new(v.z, 0.0, -v.x) }
    } else if a <= b && a <= c {
        if b > c { Vector::new(0.0, -v.z, v.y) } else { Vector::new(0.0, v.z, -v.y) }
    } else if a > b {
        Vector::new(-v.y, v.x, 0.0)
    } else {
        Vector::new(v.y, -v.x, 0.0)
    };
    let x = x * (1.0 / x.length());
    let z = v * (1.0 / v.length());
    z.cross(x)
}

/// `CosAngle`: the cosine between two vectors, 1 where either is nought.
fn cos_angle(a: Vector, b: Vector) -> f64 {
    let (m1, m2) = (a.length(), b.length());
    if m1 <= RESOLUTION || m2 <= RESOLUTION {
        return 1.0;
    }
    (a.dot(b) / (m1 * m2)).clamp(-1.0, 1.0)
}

/// `GeomFill_SnglrFunc`: the curve `ratio · D1 × D2` of a curve, whose zeros are where the
/// curve's curvature vanishes.
pub(super) struct SingularFunction<'a> {
    pub curve: &'a GeomAdaptor,
    pub ratio: f64,
}

impl CurveAdaptor for SingularFunction<'_> {
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let of_curve = match continuity {
            Continuity::C0 => Continuity::C2,
            Continuity::C1 => Continuity::C3,
            Continuity::C2 | Continuity::C3 | Continuity::CN => Continuity::CN,
            // (G1 and G2: OCCT leaves C0.)
            Continuity::G1 | Continuity::G2 => Continuity::C0,
        };
        self.curve.intervals(of_curve)
    }

    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.curve.period()
    }

    fn value(&self, u: f64) -> Point {
        let (_, d1, d2) = self.curve.d2(u);
        let v = (d1 * self.ratio).cross(d2);
        Point::new(v.x, v.y, v.z)
    }

    fn d1(&self, u: f64) -> (Point, Vector) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let dc = d1 * self.ratio;
        let v = dc.cross(d2);
        (Point::new(v.x, v.y, v.z), dc.cross(d3))
    }

    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let d4 = self.curve.dn(u, 4);
        let v = d1.cross(d2) * self.ratio;
        (Point::new(v.x, v.y, v.z), d1.cross(d3) * self.ratio, (d2.cross(d3) + d1.cross(d4)) * self.ratio)
    }

    fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
        let (_, d1, d2, d3) = self.curve.d3(u);
        let (d4, d5) = (self.curve.dn(u, 4), self.curve.dn(u, 5));
        let v = d1.cross(d2) * self.ratio;
        (Point::new(v.x, v.y, v.z), d1.cross(d3) * self.ratio, (d2.cross(d3) + d1.cross(d4)) * self.ratio, (d1.cross(d5) + d2.cross(d4) * 2.0) * self.ratio)
    }

    fn dn(&self, u: f64, n: u32) -> Vector {
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            3 => self.d3(u).3,
            // (OCCT throws Standard_NotImplemented; nothing here asks for more.)
            _ => Vector::new(f64::NAN, f64::NAN, f64::NAN),
        }
    }

    /// `Precision::Parametric(R3D)`.
    fn resolution(&self, r3d: f64) -> f64 {
        r3d * 0.01
    }
}

/// `GeomFill_Frenet` on one curve.
#[derive(Clone, Debug)]
pub(super) struct Frenet {
    /// `myCurve`.
    curve: GeomAdaptor,
    /// `mySngl` and `mySnglLen`: the singular points and the half-lengths of the intervals
    /// about them where the trihedron is taken from higher derivatives (`isSngl` where any).
    singular: Vec<(f64, f64)>,
}

/// A trihedron: tangent, normal, binormal.
pub(super) type Trihedron = (Vector, Vector, Vector);

impl Frenet {
    /// `SetCurve`: the singular points searched for unless the curve is a line or a conic.
    pub(super) fn new(curve: GeomAdaptor) -> Frenet {
        let mut frenet = Frenet { curve, singular: Vec::new() };
        if !matches!(frenet.curve.kind(), CurveKind::Line | CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola) {
            frenet.init();
        }
        frenet
    }

    pub(super) fn curve(&self) -> &GeomAdaptor {
        &self.curve
    }

    /// `Init`: the singular points.
    fn init(&mut self) {
        const TOL_F: f64 = 1.0e-10;
        const TOL: f64 = 10.0 * TOL_F;
        const TOL2: f64 = TOL * TOL;
        let mut func = SingularFunction { curve: &self.curve, ratio: 1.0 };
        // (Whether the curve has straight stretches.)
        let c2_disc = self.curve.intervals(Continuity::C2);
        let nb = c2_disc.len() - 1;
        let mut is_lin = vec![true; nb];
        let mut is_const = vec![true; nb];
        let mut ave_func = vec![0.0; nb];
        const NB_CONTROL: usize = 10;
        // (Summed over all the intervals so far, as OCCT does.)
        let mut average = 0.0;
        for i in 0..nb {
            let step = (c2_disc[i + 1] - c2_disc[i]) / NB_CONTROL as f64;
            let mut c1 = Point::new(0.0, 0.0, 0.0);
            for j in 1..=NB_CONTROL {
                let c = func.value(c2_disc[i] + (j - 1) as f64 * step);
                if j == 1 {
                    c1 = c;
                }
                let modulus = Vector::new(c.x, c.y, c.z).length();
                if modulus > TOL {
                    is_lin[i] = false;
                }
                average += modulus;
                if is_const[i] && ((c.x - c1.x).abs() > TOL || (c.y - c1.y).abs() > TOL || (c.z - c1.z).abs() > TOL) {
                    is_const[i] = false;
                }
            }
            ave_func[i] = average / NB_CONTROL as f64;
        }
        // The singularities.
        let origin = Point::new(0.0, 0.0, 0.0);
        let mut seq: Vec<Vec<f64>> = vec![Vec::new(); nb];
        let mut sngl: Vec<f64> = Vec::new();
        for i in 0..nb {
            if !is_lin[i] && !is_const[i] {
                func.ratio = 1.0 / ave_func[i];
                if let Some(ext) = ext_pc_other(&func, origin, c2_disc[i], c2_disc[i + 1], TOL_F) {
                    for e in ext {
                        if e.square_distance < TOL2 {
                            seq[i].push(e.parameter);
                        }
                    }
                }
                seq[i].sort_by(f64::total_cmp);
            }
        }
        for s in &seq {
            sngl.extend_from_slice(s);
        }
        // ("Extrema works bad, need to pass second time".)
        for i in 0..nb {
            if seq[i].is_empty() {
                continue;
            }
            let mut bounds = vec![c2_disc[i]];
            bounds.extend_from_slice(&seq[i]);
            bounds.push(c2_disc[i + 1]);
            func.ratio = 1.0 / ave_func[i];
            for j in 0..bounds.len() - 1 {
                if bounds[j + 1] - bounds[j] > P_CONFUSION {
                    if let Some(ext) = ext_pc_other(&func, origin, bounds[j], bounds[j + 1], TOL_F) {
                        for e in ext {
                            let t = e.parameter;
                            if e.square_distance < TOL2 && t - bounds[j] > P_CONFUSION && bounds[j + 1] - t > P_CONFUSION {
                                sngl.push(t);
                            }
                        }
                    }
                }
            }
        }
        if sngl.is_empty() {
            return;
        }
        sngl.sort_by(f64::total_cmp);
        // (Repeats dropped.)
        let mut j = 0;
        loop {
            let mut found = false;
            let mut i = j;
            while i + 1 < sngl.len() {
                if sngl[i + 1] - sngl[i] <= P_CONFUSION {
                    sngl.remove(i + 1);
                    j = i;
                    found = true;
                    break;
                }
                i += 1;
            }
            if !found {
                break;
            }
        }
        // The lengths of the singular intervals (with the ratio last set).
        let mut singular: Vec<(f64, f64)> = sngl
            .iter()
            .map(|&s| {
                let (_, der, der2) = func.d2(s);
                let length = if der.length() > RESOLUTION {
                    (NULL_TOL / der.length()).min(MAX_SINGULAR)
                } else if der2.length() > RESOLUTION {
                    (2.0 * NULL_TOL / der2.length()).sqrt().min(MAX_SINGULAR)
                } else {
                    MAX_SINGULAR
                };
                (s, length)
            })
            .collect();
        if singular.len() > 1 {
            // (Singular points whose intervals overlap merged.)
            let mut merged: Vec<(f64, f64)> = vec![singular[0]];
            for &(s, l) in &singular[1..] {
                let last = *merged.last().expect("one at least");
                let u12 = last.0 + last.1;
                let u21 = s - l;
                if u12 >= u21 {
                    let u11 = last.0 - last.1;
                    let u22 = s + l;
                    *merged.last_mut().expect("one at least") = ((u11 + u22) / 2.0, (u22 - u11) / 2.0);
                } else {
                    merged.push((s, l));
                }
            }
            singular = merged;
        }
        self.singular = singular;
    }

    /// `NbIntervals(S)`, `Intervals(T, S)`: the curve's intervals of two orders more, the
    /// singular points added (`GeomLib::FuseIntervals`, kept to the first's values).
    pub(super) fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let of_curve = match continuity {
            Continuity::C0 => Continuity::C2,
            Continuity::C1 => Continuity::C3,
            // (OCCT throws Standard_OutOfRange for G1 and G2.)
            _ => Continuity::CN,
        };
        let trimmed = self.curve.intervals(of_curve);
        if self.singular.is_empty() {
            return trimmed;
        }
        let points: Vec<f64> = self.singular.iter().map(|s| s.0).collect();
        fuse_intervals(&trimmed, &points, P_CONFUSION, true)
    }

    /// `IsSingular`: the singular point whose interval holds `u`.
    fn is_singular(&self, u: f64) -> Option<usize> {
        self.singular.iter().position(|&(s, l)| (u - s).abs() < l)
    }

    /// `D0` with `trimmed` as `myTrimmed` (the interval `SetInterval` last set): the
    /// trihedron at `param`, and whether OCCT reports it found.
    pub(super) fn d0(&self, trimmed: &GeomAdaptor, param: f64) -> (bool, Trihedron) {
        let mut delta = 0.0;
        if let Some(index) = self.is_singular(param) {
            let (ok, trihedron, d) = self.singular_d0(trimmed, param, index);
            delta = d;
            if ok {
                return (true, trihedron);
            }
        }
        let a_param = param + delta;
        let (point, tangent, d2) = trimmed.d2(a_param);
        let _ = point;
        let division_factor = 1.0e-3;
        let (u_inf, u_sup) = (trimmed.first_parameter(), trimmed.last_parameter());
        let a_delta = (u_sup - u_inf) * division_factor;
        let ndu = tangent.length();
        if ndu <= RESOLUTION {
            // The derivative approximated by the Taylor series.
            let mut index = 1;
            let (mut tn, found) = loop {
                index += 1;
                let tn = trimmed.dn(param, index);
                let found = tn.length() > RESOLUTION;
                if found || index >= MAX_DERIV_ORDER {
                    break (tn, found);
                }
            };
            if found {
                let u = if param - u_inf < a_delta { param + a_delta } else { param - a_delta };
                let p1 = trimmed.value(param.min(u));
                let p2 = trimmed.value(param.max(u));
                if tn.dot(p2 - p1) < 0.0 {
                    tn = -tn;
                }
            } else {
                // By three points.
                let (p1, p2, p3, grown) = if param - u_inf < 2.0 * a_delta {
                    (trimmed.value(param), trimmed.value(param + a_delta), trimmed.value(param + 2.0 * a_delta), true)
                } else {
                    (trimmed.value(param - 2.0 * a_delta), trimmed.value(param - a_delta), trimmed.value(param), false)
                };
                let v = |p: Point| Vector::new(p.x, p.y, p.z);
                tn = if grown { v(p1) * -3.0 + v(p2) * 4.0 - v(p3) } else { v(p1) - v(p2) * 4.0 + v(p3) * 3.0 };
            }
            // (The trihedron of a point nearby, turned onto that tangent.)
            let d_par = 10.0 * a_delta;
            let (ok, near) = if param - u_inf < d_par { self.d0(trimmed, a_param + d_par) } else { self.d0(trimmed, a_param - d_par) };
            if !ok {
                return (false, near);
            }
            let (ok, turned) = rotate_trihedron(near, tn);
            return (ok, turned);
        }
        let tangent = tangent / ndu;
        let binormal = tangent.cross(d2);
        let norm = binormal.length();
        let binormal = if norm <= RESOLUTION { default_y(tangent) } else { binormal / norm };
        let normal = binormal.cross(tangent);
        (true, (tangent, normal, binormal))
    }

    /// `DoSingular`: at the singular point `index`, the tangent and binormal from the first
    /// derivatives that do not vanish, signed as at a point beside it; `delta` the step off
    /// the point where none serve.
    fn do_singular(&self, trimmed: &GeomAdaptor, u: f64, index: usize) -> (bool, Vector, Vector, i32, i32, f64) {
        const MAX_N: u32 = 20;
        let mut h = 2.0 * self.singular[index].1;
        let (a, b) = (trimmed.first_parameter(), trimmed.last_parameter());
        if u >= (a + b) / 2.0 {
            h = -h;
        }
        let mut tangent = Vector::new(0.0, 0.0, 0.0);
        let mut binormal = Vector::new(0.0, 0.0, 0.0);
        let mut i = 1;
        while i <= MAX_N {
            tangent = trimmed.dn(u, i);
            if tangent.length() > CONFUSION {
                break;
            }
            i += 1;
        }
        if i > MAX_N {
            return (false, tangent, binormal, 1, 1, 0.0);
        }
        tangent = tangent / tangent.length();
        i += 1;
        while i <= MAX_N {
            binormal = tangent.cross(trimmed.dn(u, i));
            let magnitude = binormal.length();
            if magnitude > CONFUSION {
                // (OCC605.)
                let next = tangent.cross(trimmed.dn(u, i + 1));
                if next.length() > magnitude {
                    i += 1;
                    binormal = next;
                }
                break;
            }
            i += 1;
        }
        if i > MAX_N {
            return (false, tangent, binormal, 1, 1, h);
        }
        binormal = binormal / binormal.length();
        let (_, (t, _, bn)) = self.d0(trimmed, u + h);
        let t_flag = if angle(tangent, t) > std::f64::consts::FRAC_PI_2 { -1 } else { 1 };
        let bn_flag = if angle(binormal, bn) > std::f64::consts::FRAC_PI_2 { -1 } else { 1 };
        (true, tangent, binormal, t_flag, bn_flag, 0.0)
    }

    /// `SingularD0`.
    fn singular_d0(&self, trimmed: &GeomAdaptor, param: f64, index: usize) -> (bool, Trihedron, f64) {
        let (ok, tangent, binormal, t_flag, bn_flag, delta) = self.do_singular(trimmed, param, index);
        if !ok {
            return (false, (tangent, Vector::new(0.0, 0.0, 0.0), binormal), delta);
        }
        let tangent = tangent * t_flag as f64;
        let binormal = binormal * bn_flag as f64;
        (true, (tangent, binormal.cross(tangent), binormal), delta)
    }
}

/// `RotateTrihedron`: the trihedron turned to bring its tangent onto `new_tangent`, and
/// whether it got there.
fn rotate_trihedron((tangent, normal, binormal): Trihedron, new_tangent: Vector) -> (bool, Trihedron) {
    let inf_cos = ANGULAR.cos();
    let axis = tangent.cross(new_tangent);
    let nt = axis.length();
    if nt <= RESOLUTION {
        return (true, (tangent, normal, binormal));
    }
    let axis = axis / nt;
    let (px, py, pz) = (axis.x, axis.y, axis.z);
    let c = cos_angle(tangent, new_tangent);
    let add = 1.0 - c;
    let s = (1.0 - c * c).sqrt();
    let v11 = Vector::new(add * px * px + c, add * px * py - pz * s, add * px * pz + py * s);
    let v12 = Vector::new(add * px * px + c, add * px * py + pz * s, add * px * pz - py * s);
    let v21 = Vector::new(add * px * py + pz * s, add * py * py + c, add * py * pz - px * s);
    let v22 = Vector::new(add * px * py - pz * s, add * py * py + c, add * py * pz + px * s);
    let v31 = Vector::new(add * px * pz - py * s, add * py * pz + px * s, add * pz * pz + c);
    let v32 = Vector::new(add * px * pz + py * s, add * py * pz - px * s, add * pz * pz + c);
    let by = |v: Vector, a: Vector, b: Vector, c: Vector| Vector::new(v.dot(a), v.dot(b), v.dot(c));
    let t1 = by(tangent, v11, v21, v31);
    let t2 = by(tangent, v12, v22, v32);
    let turned = if cos_angle(t1, new_tangent) >= cos_angle(t2, new_tangent) {
        (t1, by(normal, v11, v21, v31), by(binormal, v11, v21, v31))
    } else {
        (t2, by(normal, v12, v22, v32), by(binormal, v12, v22, v32))
    };
    (cos_angle(turned.0, new_tangent) >= inf_cos, turned)
}

/// `GeomLib::FuseIntervals`: two increasing lists merged, values within `confusion` taken
/// once (the first's, or their mean); what is left of either appended (the second's first).
pub(super) fn fuse_intervals(first: &[f64], second: &[f64], confusion: f64, adjust_to_first: bool) -> Vec<f64> {
    let mut out = Vec::with_capacity(first.len() + second.len());
    let (mut i, mut j) = (0, 0);
    while i < first.len() && j < second.len() {
        let (a, b) = (first[i], second[j]);
        if (a - b).abs() <= confusion {
            out.push(if adjust_to_first { a } else { 0.5 * (a + b) });
            i += 1;
            j += 1;
        } else if a < b {
            out.push(a);
            i += 1;
        } else {
            out.push(b);
            j += 1;
        }
    }
    out.extend_from_slice(&second[j..]);
    out.extend_from_slice(&first[i..]);
    out
}
