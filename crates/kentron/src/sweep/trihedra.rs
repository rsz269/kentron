//! OCCT's other pipe trihedra, ported: the constant binormal one (`GeomFill_ConstantBiNormal`,
//! the Frenet frame turned about a fixed binormal) and the discrete one
//! (`GeomFill_DiscreteTrihedron`, a frame carried from sample to sample along the curve by
//! the least rotations taking each tangent onto the next).

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;
use crate::math::Vector;

use super::corrected::edge_adaptor;
use super::frenet::{direction, Frenet, CONFUSION, P_CONFUSION, RESOLUTION};

/// `gp_Vec::Normalized`: `v` divided by its length.
fn normalized(v: Vector) -> Vector {
    let d = v.length();
    Vector::new(v.x / d, v.y / d, v.z / d)
}

/// `gp_Dir::Crossed`: the cross product made a unit vector (as `gp_Dir` makes one).
fn dir_cross(a: Vector, b: Vector) -> Vector {
    direction(a.cross(b))
}

/// `gp_Dir::AngleWithRef`: the angle from `a` to `b` (unit vectors), in `(-π, π]`, positive
/// where `a × b` points the way `reference` does.
#[allow(clippy::approx_constant)]
fn angle_with_ref(a: Vector, b: Vector, reference: Vector) -> f64 {
    let cross = a.cross(b);
    let cos = a.dot(b);
    let sin = cross.length();
    // (OCCT's constant, not quite 1/√2.)
    let angle = if cos > -0.70710678118655 && cos < 0.70710678118655 {
        cos.acos()
    } else if cos < 0.0 {
        std::f64::consts::PI - sin.asin()
    } else {
        sin.asin()
    };
    if cross.dot(reference) >= 0.0 { angle } else { -angle }
}

/// `gp_Mat::SetRotation(Axis, Ang)` applied to `v` (`gp_Dir::Rotate` by way of `gp_Trsf`): `v`
/// turned by `angle` about `axis` (a unit vector), not renormalized.
fn rotated(v: Vector, axis: Vector, angle: f64) -> Vector {
    let (a, b, c) = (axis.x, axis.y, axis.z);
    let (cos, sin) = (angle.cos(), angle.sin());
    let om_cos = 1.0 - cos;
    let (a2, b2, c2, ab, ac, bc) = (a * a, b * b, c * c, a * b, a * c, b * c);
    let m = [
        [1.0 + om_cos * (-(b2 + c2)), om_cos * ab - sin * c, om_cos * ac + sin * b],
        [om_cos * ab + sin * c, 1.0 + om_cos * (-(a2 + c2)), om_cos * bc - sin * a],
        [om_cos * ac - sin * b, om_cos * bc + sin * a, 1.0 + om_cos * (-(a2 + b2))],
    ];
    Vector::new(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z, m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z, m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z)
}

/// The directions of a `gp_Ax2` (its origin is always the world's here): main, X and Y.
#[doc(alias = "gp_Ax2")]
#[derive(Clone, Copy, Debug)]
struct Axes {
    main: Vector,
    x: Vector,
    y: Vector,
}

impl Axes {
    /// `gp_Ax2(P, V)`: the X direction square to `v` (a unit vector) with the smallest of its
    /// components dropped (`SetXDirection`).
    fn new(v: Vector) -> Axes {
        let (a, b, c) = (v.x, v.y, v.z);
        let (a_abs, b_abs, c_abs) = (a.abs(), b.abs(), c.abs());
        let d = if b_abs <= a_abs && b_abs <= c_abs {
            if a_abs > c_abs { Vector::new(-c, 0.0, a) } else { Vector::new(c, 0.0, -a) }
        } else if a_abs <= b_abs && a_abs <= c_abs {
            if b_abs > c_abs { Vector::new(0.0, -c, b) } else { Vector::new(0.0, c, -b) }
        } else if a_abs > b_abs {
            Vector::new(-b, a, 0.0)
        } else {
            Vector::new(b, -a, 0.0)
        };
        let d = direction(d);
        // (`SetXDirection`: X = V ^ (D ^ V), Y = V ^ X.)
        let x = direction(v.cross(d.cross(v)));
        Axes { main: v, x, y: dir_cross(v, x) }
    }

    /// `Rotated(gp_Ax1(Origin, axis), angle)`: X and Y turned, the main direction their cross
    /// product.
    fn rotated(&self, axis: Vector, angle: f64) -> Axes {
        let axis = direction(axis);
        let (x, y) = (rotated(self.x, axis, angle), rotated(self.y, axis, angle));
        Axes { main: dir_cross(x, y), x, y }
    }

    /// `SetDirection(V)`: the main direction `v` (a unit vector), X brought square to it.
    fn set_direction(&mut self, v: Vector) {
        const ANGULAR: f64 = 1e-12;
        let a = v.dot(self.x);
        if (a.abs() - 1.0).abs() <= ANGULAR {
            if a > 0.0 {
                self.x = self.y;
                self.y = self.main;
            } else {
                self.x = self.main;
            }
            self.main = v;
        } else {
            self.main = v;
            self.x = direction(v.cross(self.x.cross(v)));
            self.y = dir_cross(v, self.x);
        }
    }
}

/// The interval of `intervals` holding `param` (the first of two at a bound).
fn interval_holding(intervals: &[f64], param: f64) -> usize {
    intervals.windows(2).position(|w| w[0] <= param && param <= w[1]).unwrap_or(if param < intervals[0] { 0 } else { intervals.len() - 2 })
}

/// `GeomFill_ConstantBiNormal` (a pipe shell's binormal mode), along `curve` run from `first`
/// to `last`: the Frenet frame's tangent turned square to a fixed binormal.
#[doc(alias = "GeomFill_ConstantBiNormal")]
#[derive(Clone, Debug)]
pub(super) struct ConstantBinormalLaw {
    frenet: Frenet,
    binormal: Vector,
}

impl ConstantBinormalLaw {
    pub(super) fn new(curve: &Curve, first: f64, last: f64, binormal: Vector) -> ConstantBinormalLaw {
        ConstantBinormalLaw { frenet: Frenet::new(edge_adaptor(curve, first, last)), binormal: direction(binormal) }
    }

    /// `D0`: tangent and normal. (Its Frenet law is never trimmed: `SetInterval` trims the
    /// law's own curve, not the Frenet law's.)
    pub(super) fn d0(&self, param: f64) -> (Vector, Vector) {
        let (_, (tangent, normal, _)) = self.frenet.d0(self.frenet.curve(), param);
        let b = self.binormal;
        if b.cross(tangent).length() > CONFUSION {
            let normal = normalized(b.cross(tangent));
            (normal.cross(b), normal)
        } else {
            let tangent = normalized(normal.cross(b));
            (tangent, b.cross(tangent))
        }
    }

    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let curve = self.frenet.curve();
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        self.d0(first + (last - first) * fraction)
    }
}

/// `GeomFill_DiscreteTrihedron` (a pipe shell's discrete mode), along `curve` run from
/// `first` to `last`: lines and conics keep their Frenet frame; along other curves, a frame
/// set at ten samples per interval of the curve's (infinite) continuity, each sample's the
/// last one's turned by the least rotation taking its tangent onto this one's, and between
/// samples the frame at the sample before turned likewise.
#[doc(alias = "GeomFill_DiscreteTrihedron")]
#[derive(Clone, Debug)]
pub(super) struct DiscreteLaw {
    curve: GeomAdaptor,
    /// `myUseFrenet`: the Frenet law, for lines and conics.
    frenet: Option<Frenet>,
    /// `myKnots` and `myTrihedrons`.
    knots: Vec<f64>,
    axes: Vec<Axes>,
}

/// `NbSamples`.
const NB_SAMPLES: usize = 10;

impl DiscreteLaw {
    /// `SetCurve`, with `Init` for curves other than lines and conics; an error where OCCT
    /// throws (tangents turning right round between samples, a sample where the tangent
    /// vanishes at the end).
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> Result<DiscreteLaw> {
        let curve = edge_adaptor(curve, first, last);
        if matches!(curve.kind(), CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola | CurveKind::Line) {
            return Ok(DiscreteLaw { frenet: Some(Frenet::new(curve.clone())), curve, knots: Vec::new(), axes: Vec::new() });
        }
        let mut law = DiscreteLaw { curve, frenet: None, knots: Vec::new(), axes: Vec::new() };
        law.init()?;
        Ok(law)
    }

    /// `Init`.
    fn init(&mut self) -> Result<()> {
        let intervals = self.curve.intervals(Continuity::CN);
        for w in intervals.windows(2) {
            let delta = (w[1] - w[0]) / NB_SAMPLES as f64;
            for j in 0..NB_SAMPLES {
                self.knots.push(w[0] + j as f64 * delta);
            }
        }
        self.knots.push(*intervals.last().expect("intervals"));
        let mut i = 0;
        while i < self.knots.len() {
            let param = self.knots[i];
            let (point, mut tangent) = self.curve.d1(param);
            if tangent.length() < CONFUSION {
                // (OCCT reads past its knots here at the last one: Standard_OutOfRange.)
                let next = *self.knots.get(i + 1).ok_or_else(|| Error::Failed("a discrete trihedron's path stops at its end".into()))?;
                let mut subdelta = (next - param) / NB_SAMPLES as f64;
                if subdelta < P_CONFUSION {
                    subdelta = next - param;
                }
                tangent = self.curve.value(param + subdelta) - point;
            }
            // (`gp_Dir(Tangent)`: Standard_ConstructionError where it is nought.)
            if tangent.length() <= RESOLUTION {
                return Err(Error::Failed("a discrete trihedron's path with no tangent".into()));
            }
            let tangent = direction(tangent);
            if i == 0 {
                self.axes.push(Axes::new(tangent));
            } else {
                let last = *self.axes.last().expect("a trihedron");
                let axis = last.main.cross(tangent);
                if axis.length() <= RESOLUTION {
                    if last.main.dot(tangent) > 0.0 {
                        self.axes.push(last);
                    } else {
                        // (Opposite tangents: a sample put halfway back, and this one again.)
                        let new_param = (self.knots[i - 1] + self.knots[i]) / 2.0;
                        if new_param - self.knots[i - 1] < RESOLUTION {
                            return Err(Error::Failed("GeomFill_DiscreteTrihedron : impassable singularities on path curve".into()));
                        }
                        self.knots.insert(i, new_param);
                        continue;
                    }
                } else {
                    let angle = angle_with_ref(last.main, tangent, direction(axis));
                    let mut axes = last.rotated(axis, angle);
                    axes.set_direction(tangent);
                    self.axes.push(axes);
                }
            }
            i += 1;
        }
        Ok(())
    }

    /// `Intervals`: the curve's of infinite continuity, whatever is asked.
    pub(super) fn intervals(&self) -> Vec<f64> {
        self.curve.intervals(Continuity::CN)
    }

    /// `D0` with `trimmed` as `myTrimmed`: tangent and normal.
    fn d0(&self, trimmed: &GeomAdaptor, param: f64) -> (Vector, Vector) {
        if let Some(frenet) = &self.frenet {
            let (_, (t, n, _)) = frenet.d0(frenet.curve(), param);
            return (t, n);
        }
        let knots = &self.knots;
        let (mut i1, mut i2) = (0, knots.len() - 1);
        loop {
            let i = (i1 + i2) / 2;
            if param <= knots[i] {
                i2 = i;
            } else {
                i1 = i;
            }
            if i2 - i1 <= 1 {
                break;
            }
        }
        let index = if (param - knots[i2]).abs() < P_CONFUSION { i2 } else { i1 };
        let previous = self.axes[index];
        let axes = if (param - knots[index]).abs() < P_CONFUSION {
            previous
        } else {
            let (point, mut tangent) = trimmed.d1(param);
            if tangent.length() < CONFUSION {
                let mut subdelta = (knots[index + 1] - param) / NB_SAMPLES as f64;
                if subdelta < P_CONFUSION {
                    subdelta = knots[index + 1] - param;
                }
                tangent = trimmed.value(param + subdelta) - point;
            }
            let tangent = direction(tangent);
            let axis = previous.main.cross(tangent);
            let mut axes = if axis.length() <= RESOLUTION {
                previous
            } else {
                let angle = angle_with_ref(previous.main, tangent, direction(axis));
                previous.rotated(axis, angle)
            };
            axes.set_direction(tangent);
            axes
        };
        (axes.main, axes.x)
    }

    /// The frame's tangent and normal at fraction `fraction` of the way along (evaluated in
    /// the law's interval holding it, as a sweep sets it).
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let (first, last) = (self.curve.first_parameter(), self.curve.last_parameter());
        let param = first + (last - first) * fraction;
        let intervals = self.intervals();
        let k = interval_holding(&intervals, param);
        self.d0(&self.curve.trim(intervals[k], intervals[k + 1]), param)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Line, OffsetCurve};
    use crate::math::{Direction, Frame, Point};
    use crate::topology::{Model, NodeKind};

    /// The binormal and discrete trihedra match OCCT 8.1's (`GeomFill_ConstantBiNormal`,
    /// `GeomFill_DiscreteTrihedron`, set on the curve as `BRepFill_Edge3DLaw` sets them, each
    /// evaluated in its interval) at the start, a fifth, half, 85% and the end: B-splines
    /// planar and twisted (run backwards), a Bézier, one with a C1 knot both ways, a quintic
    /// and a periodic closed curve run backwards, an offset, a circle and a line (the binormal
    /// along it and across it).
    #[test]
    fn trihedra_match_occts() {
        let p = Point::new;
        let v = Vector::new;
        let cubic = |poles: Vec<Point>, knots: Vec<f64>, mults: Vec<usize>| Curve::BSpline(BSplineCurve::new(poles, knots, mults, 3));
        let planar_s = cubic(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 2.0, 0.0), p(3.0, -2.0, 0.0), p(4.0, -2.0, 0.0), p(5.0, 0.0, 0.0)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
        let twisted_s = cubic(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.5), p(2.0, 2.0, 1.0), p(3.0, -2.0, 1.5), p(4.0, -2.0, 2.0), p(5.0, 0.0, 2.5)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
        let bezier = Curve::Bezier(BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 1.0, 2.0), p(2.0, -1.0, 3.0), p(0.0, -1.0, 4.0)]));
        let c1_knot = cubic(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 2.0, 1.0), p(3.0, 0.0, 2.0), p(4.0, -1.0, 1.0), p(5.0, 1.0, 0.0), p(6.0, 2.0, 2.0)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 2, 1, 4]);
        let quintic = Curve::BSpline(BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, 0.0, 1.0), p(3.0, 1.0, 1.0), p(4.0, 0.0, 0.0), p(5.0, 1.0, 2.0), p(6.0, 0.0, 1.0), p(7.0, 1.0, 0.0)], vec![0.0, 1.0, 2.0, 3.0], vec![6, 1, 1, 6], 5));
        let offset = Curve::Offset(OffsetCurve::new(twisted_s.clone(), 0.3, Direction::Z));
        let normal = Direction::from_vector(v(1.0, 1.0, 1.0));
        let circle = Curve::Circle(Circle::new(Frame::from_z(p(1.0, 2.0, 3.0), normal), 2.0));
        let line = Curve::Line(Line::new(p(1.0, 0.0, 0.0), Direction::Z));
        let mut model = Model::new();
        let trefoil = crate::io::brep::read(&mut model, &mut include_str!("tests/occt_periodic_trefoil.brep").as_bytes()).unwrap();
        let (trefoil, t0, t1) = model.edge_curve(model.unique_subshapes(trefoil, NodeKind::Edge)[0]).unwrap();
        let (z, w) = (v(0.0, 0.0, 1.0), v(0.3, 0.2, 1.0));
        let binormal = |c: &Curve, a: f64, b: f64, d: Vector| -> Box<dyn Fn(f64) -> (Vector, Vector)> {
            let law = ConstantBinormalLaw::new(c, a, b, d);
            Box::new(move |f| law.at_fraction(f))
        };
        let discrete = |c: &Curve, a: f64, b: f64| -> Box<dyn Fn(f64) -> (Vector, Vector)> {
            let law = DiscreteLaw::new(c, a, b).unwrap();
            Box::new(move |f| law.at_fraction(f))
        };
        let laws = [
            binormal(&planar_s, 0.0, 3.0, z),
            binormal(&twisted_s, 3.0, 0.0, w),
            binormal(&bezier, 0.0, 1.0, w),
            binormal(&line, 0.0, 5.0, z),
            binormal(&line, 0.0, 5.0, v(1.0, 0.0, 0.0)),
            binormal(&circle, 0.0, 4.0, w),
            discrete(&planar_s, 0.0, 3.0),
            discrete(&twisted_s, 3.0, 0.0),
            discrete(&bezier, 0.0, 1.0),
            discrete(&c1_knot, 0.0, 3.0),
            discrete(&c1_knot, 3.0, 0.0),
            discrete(&quintic, 3.0, 0.0),
            discrete(&trefoil, t1, t0),
            discrete(&offset, 0.0, 3.0),
            discrete(&circle, 0.0, 4.0),
        ];
        // (OCCT's tangent and normal at each fraction.)
        let occt: [(&str, [(Vector, Vector); 5]); 15] = [
            ("bn_planar", [(v(0.447213595499958, 0.894427190999916, -0.0), v(-0.894427190999916, 0.447213595499958, 0.0)), (v(0.9895864392845115, 0.14393984571411012, -0.0), v(-0.14393984571411012, 0.9895864392845115, 0.0)), (v(0.3511234415883917, -0.9363291775690445, 0.0), v(0.9363291775690445, 0.3511234415883917, -0.0)), (v(0.8035010957679819, 0.5953032748941102, -0.0), v(-0.5953032748941102, 0.8035010957679819, 0.0)), (v(0.447213595499958, 0.894427190999916, -0.0), v(-0.894427190999916, 0.447213595499958, 0.0))]),
            ("bn_twisted_rev", [(v(-0.3417491637537273, -0.8965367672500377, 0.28183210257612573), v(0.8964159142800361, -0.4010281721779109, -0.18871913984842864)), (v(-0.9579120738331143, 0.0015810391150544376, 0.2870574143269234), v(0.052520863527052236, -0.9821401479558918, 0.1806717705330627)), (v(-0.32256660689238764, 0.9420992963206242, -0.09164987719640855), v(-0.9034958585510046, -0.2775800529283207, 0.3265647681509655)), (v(-0.7470390662998554, -0.5721238333518617, 0.3385364865603289), v(0.6019024968570287, -0.798295740451041, -0.020911600966900434)), (v(-0.3417491637537273, -0.8965367672500377, 0.28183210257612573), v(0.8964159142800361, -0.4010281721779109, -0.18871913984842864))]),
            ("bn_bezier", [(v(0.38099892893717247, 0.8779540536378322, -0.2898904894087182), v(-0.8804509063256238, 0.4402254531628119, 0.17609018126512474)), (v(0.8702120912343603, 0.36256985777239614, -0.33357759892478733), v(-0.40383771316013406, 0.9127676967040459, -0.06140222539276898)), (v(-0.05006951972142579, -0.9763556345678034, 0.2102919828299884), v(0.9580433317422291, -0.10644925908246987, -0.2661231477061747)), (v(-0.8039422189321461, -0.4886826095728545, 0.33891918759421474), v(0.5234795993344367, -0.8519337280637325, 0.013342865812415495)), (v(-0.955545921346429, -0.03211919063349341, 0.29308761453062737), v(0.08535791995526398, -0.9816160794855358, 0.17071583991052797))]),
            ("bn_line_along", [(v(0.0, -1.0, 0.0), v(1.0, 0.0, -0.0)), (v(0.0, -1.0, 0.0), v(1.0, 0.0, -0.0)), (v(0.0, -1.0, 0.0), v(1.0, 0.0, -0.0)), (v(0.0, -1.0, 0.0), v(1.0, 0.0, -0.0)), (v(0.0, -1.0, 0.0), v(1.0, 0.0, -0.0))]),
            ("bn_line", [(v(-0.0, 0.0, 1.0), v(0.0, -1.0, 0.0)), (v(-0.0, 0.0, 1.0), v(0.0, -1.0, 0.0)), (v(-0.0, 0.0, 1.0), v(0.0, -1.0, 0.0)), (v(-0.0, 0.0, 1.0), v(0.0, -1.0, 0.0)), (v(-0.0, 0.0, 1.0), v(0.0, -1.0, 0.0))]),
            ("bn_circle", [(v(-0.33110982210466783, 0.939427867366732, -0.08855262684194605), v(-0.9004000599770637, -0.28649092817452027, 0.3274182036280232)), (v(-0.821559812100556, 0.5537320094114848, 0.13572154174786982), v(-0.4953720394230829, -0.811161285851725, 0.31084386899726985)), (v(-0.7688169986471152, -0.5421679326611458, 0.3390786861263637), v(0.5738243676353653, -0.8189357134143337, -0.008360167607742828)), (v(0.527006951462552, -0.8497783039484753, 0.011853575350929457), v(0.8016343651634655, 0.49242116531477126, -0.33897454261199383)), (v(0.8444596476875046, -0.5140312515636614, -0.1505316439935191), v(0.45523827357591606, 0.836883290744534, -0.3039481402216816))]),
            ("dt_planar", [(v(0.4472135954999579, 0.8944271909999159, 0.0), v(0.894427190999916, -0.447213595499958, 0.0)), (v(0.9895864392845115, 0.14393984571411012, 0.0), v(0.14393984571411012, -0.9895864392845115, 0.0)), (v(0.35112344158839165, -0.9363291775690445, 0.0), v(-0.9363291775690445, -0.35112344158839165, 0.0)), (v(0.8035010957679819, 0.5953032748941102, 0.0), v(0.5953032748941102, -0.8035010957679819, 0.0)), (v(0.4472135954999579, 0.8944271909999159, 0.0), v(0.894427190999916, -0.447213595499958, 0.0))]),
            ("dt_twisted_rev", [(v(-0.4364357804719847, -0.8728715609439694, -0.21821789023599236), v(-0.8944271909999159, 0.4472135954999579, 0.0)), (v(-0.8869525753639098, -0.12901128368929538, -0.4434762876819549), v(-0.2846432199773797, 0.9088559371495054, 0.30489198551126806)), (v(-0.3458342506225409, 0.9222246683267759, -0.17291712531127046), v(0.5771138235198484, 0.354374653931171, 0.7357705072598821)), (v(-0.7455811295485929, -0.5523911422862476, -0.37279056477429645), v(-0.6317114252847321, 0.7639933126512977, 0.13135788285759162)), (v(-0.4364357804719847, -0.8728715609439694, -0.21821789023599236), v(-0.894427190999916, 0.44721359549995804, -2.4980018054066027e-16))]),
            ("dt_bezier", [(v(0.4472135954999579, 0.8944271909999159, 0.0), v(0.894427190999916, -0.447213595499958, 0.0)), (v(0.6679929279245405, 0.3400691269434024, 0.6619202649434084), v(0.5842683086035676, -0.7905488302053639, -0.18347504087056926)), (v(0.16169041669088866, -0.5659164584181103, 0.8084520834544433), v(0.32677391179761855, -0.7423111670135168, -0.5849725992689856)), (v(-0.6780005385445668, -0.4070180963189924, 0.6120878523562944), v(0.19268764647091977, -0.9019962182525079, -0.3863603151925314)), (v(-0.8944271909999159, 0.0, 0.4472135954999579), v(-0.09685170270883323, -0.9762677595885255, -0.19370340541766645))]),
            ("dt_c1", [(v(0.4472135954999579, 0.8944271909999159, 0.0), v(0.894427190999916, -0.447213595499958, 0.0)), (v(0.7784498303861735, -0.03797316245786258, 0.6265571805547251), v(0.24636202185244163, -0.8995934482480552, -0.36060696340746184)), (v(0.7043607250604991, -0.6163156344279367, -0.35218036253024954), v(-0.5330208087264644, -0.1315601036746288, -0.8358114360223283)), (v(0.592264133228619, 0.7713075132353514, 0.2330405898500321), v(-0.25146419493755295, 0.4517178225305931, -0.8559887659733788)), (v(0.408248290463863, 0.408248290463863, 0.816496580927726), v(0.09432217510523663, 0.8707754888763015, -0.48254883199076903))]),
            ("dt_c1_rev", [(v(-0.408248290463863, -0.408248290463863, -0.816496580927726), v(-0.894427190999916, 0.0, 0.447213595499958)), (v(-0.6185281318792845, -0.7784923039170305, 0.10664278135849758), v(-0.49814761756522263, 0.4934577701179351, 0.7129855399826338)), (v(-0.7043607250604991, 0.6163156344279367, 0.35218036253024954), v(0.6071684129270545, 0.7801183667335125, -0.15087031592953817)), (v(-0.7807040146620593, -0.3496337868166658, -0.5179357649427014), v(0.21188372049908905, 0.6316235814655236, -0.7457593045508244)), (v(-0.4472135954999579, -0.8944271909999159, 0.0), v(-0.2684772728804972, 0.1342386364402486, -0.9538867555602781))]),
            ("dt_quintic_rev", [(v(-0.5773502691896257, -0.5773502691896257, 0.5773502691896257), v(0.7071067811865475, 0.0, 0.7071067811865475)), (v(-0.9845459045431045, 0.16778585099940516, -0.0501724032891127), v(-0.05217975506145581, -0.007573649457912455, 0.9986089890420651)), (v(-0.9986506538179852, -0.023224433809720425, -0.04644886761944138), v(-0.0472902199436763, 0.03709514019681476, 0.9981921586906287)), (v(-0.8739328047646481, 0.15350852324523281, -0.4611687175505982), v(-0.4572253154645048, 0.062221868893858764, 0.8871716011740567)), (v(-0.7071067811865476, -0.7071067811865474, 3.925231146709438e-17), v(-0.1701439258891902, 0.1701439258891903, 0.9706194357038331))]),
            ("dt_trefoil_rev", [(v(-0.6183515251469651, 0.2522130710898149, 0.7443318870771735), v(0.7691987562340656, 0.0, 0.6390096035334418)), (v(0.8032403754789572, 0.2581630090654324, -0.5368023471918739), v(0.05614760653621304, 0.8643684141378611, 0.49971460947330093)), (v(-0.7424142544184067, 0.00025715594339414125, -0.6699410486805403), v(-0.226674889153592, -0.9411165194928954, 0.2508349883185866)), (v(0.06936972458795142, -0.7855450707338075, -0.6149038812338143), v(0.7445377480416856, 0.45101801706272976, -0.49218521922729314)), (v(-0.6183515251469555, 0.25221307108982677, 0.7443318870771776), v(-0.6682858631770889, -0.6671387096108883, -0.329119958702556))]),
            ("dt_offset", [(v(0.4360542907453469, 0.8721085814906938, 0.22199837299598038), v(0.894427190999916, -0.447213595499958, 0.0)), (v(0.4525488820381729, 0.06582529193282469, 0.8893067751389045), v(0.5656388395524011, -0.7921596977651723, -0.22920671112863172)), (v(0.3458342506225409, -0.9222246683267759, 0.17291712531127046), v(-0.0803999428987475, -0.21273675484113547, -0.9737961400218943)), (v(0.7677422073708129, 0.5688099900914785, 0.29500355623877467), v(0.5434676927118929, -0.3341633677848155, -0.7700504597811247)), (v(0.43679826006250766, 0.8735965201250153, 0.2145609471265207), v(0.6855534168799582, -0.16884504279034315, -0.708172199489179))]),
            ("dt_circle", [(v(-0.4082482904638631, 0.8164965809277261, -0.4082482904638631), v(-0.7071067811865477, 0.0, 0.7071067811865477)), (v(-0.7916766794462098, 0.5688586460913676, 0.22281803335484218), v(-0.19978664091397608, -0.5857187955271395, 0.7855054364411154)), (v(-0.47307914202160245, -0.339782469204631, 0.8128616112262335), v(0.6654793701165883, -0.7424382400495483, 0.07695886993295994)), (v(0.5753885554591733, -0.7893874186882351, 0.21399886322906186), v(0.5793053399661292, 0.20864843609141087, -0.7879537760575401)), (v(0.8019890672413155, -0.5336977815803581, -0.2682912856609574), v(0.15323251186374132, 0.617926649856495, -0.7711591617202362))]),
        ];
        let fractions = [0.0, 0.2, 0.5, 0.85, 1.0];
        let mut worst = 0.0f64;
        for (law, (name, frames)) in laws.iter().zip(occt) {
            for (f, (t, n)) in fractions.iter().zip(frames) {
                let (lt, ln) = law(*f);
                let off = (lt - t).length().max((ln - n).length());
                worst = worst.max(off);
                assert!(off < 1e-9, "{name} at {f}: {lt:?}, {ln:?} vs OCCT's {t:?}, {n:?}");
            }
        }
        eprintln!("largest difference from OCCT: {worst:e}");
    }
}
