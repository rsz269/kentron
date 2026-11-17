//! OCCT's `GeomFill_PlanFunc`, ported: the signed distance, scaled by the length of `V`, from
//! the plane through `P` normal to `V` to the point of a curve, as a function of the curve's
//! parameter. `GeomFill_GuideTrihedronPlan` solves it for where a guide curve crosses the
//! plane normal to the path.

use crate::geometry::adaptor::CurveAdaptor;
use crate::geometry::function_roots::FunctionWithDerivative;
use crate::math::{Point, Vector};

/// `GeomFill_PlanFunc(P, V, C)`: `F(X) = V · (C(X) - P)`.
#[doc(alias = "GeomFill_PlanFunc")]
pub(crate) struct PlanFunc<'a> {
    /// `myPnt`, `myVec`.
    pnt: [f64; 3],
    vec: [f64; 3],
    curve: &'a dyn CurveAdaptor,
}

/// `gp_XYZ::Dot`.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `V.SetLinearForm(-1, P, G)`: `-1 · P + G`, coordinate by coordinate.
fn from_pnt(pnt: [f64; 3], g: Point) -> [f64; 3] {
    [-pnt[0] + g.x, -pnt[1] + g.y, -pnt[2] + g.z]
}

impl<'a> PlanFunc<'a> {
    pub(crate) fn new(p: Point, v: Vector, curve: &'a dyn CurveAdaptor) -> Self {
        PlanFunc { pnt: [p.x, p.y, p.z], vec: [v.x, v.y, v.z], curve }
    }

    /// `DEDT(X, DPnt, DVec, DFDT)`: the derivative of `F(X)` with respect to a parameter
    /// moving `P` and `V` at the rates `dpnt` and `dvec`, at fixed `X`.
    #[doc(alias = "GeomFill_PlanFunc::DEDT")]
    pub(crate) fn dedt(&self, x: f64, dpnt: Vector, dvec: Vector) -> f64 {
        let g = self.curve.value(x);
        let v = from_pnt(self.pnt, g);
        dot([dvec.x, dvec.y, dvec.z], v) - dot(self.vec, [dpnt.x, dpnt.y, dpnt.z])
    }

    /// `D2E(X, DP, D2P, DV, D2V, DFDT, D2FDT2, D2FDTDX)`: OCCT computes only `DFDT`,
    /// `DV · (C(X) - P) + V · (C(X) - DP)` (the two second derivatives are left unset).
    #[doc(alias = "GeomFill_PlanFunc::D2E")]
    pub(crate) fn d2e(&self, x: f64, dp: Vector, dv: Vector) -> f64 {
        let (g, _) = self.curve.d1(x);
        let v = from_pnt(self.pnt, g);
        let dvdt = [-dp.x + g.x, -dp.y + g.y, -dp.z + g.z];
        dot([dv.x, dv.y, dv.z], v) + dot(self.vec, dvdt)
    }
}

impl FunctionWithDerivative for PlanFunc<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let g = self.curve.value(x);
        let v = from_pnt(self.pnt, g);
        *f = dot(self.vec, v);
        true
    }

    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (_, dg) = self.curve.d1(x);
        *d = dot(self.vec, [dg.x, dg.y, dg.z]);
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (g, dg) = self.curve.d1(x);
        let v = from_pnt(self.pnt, g);
        *f = dot(self.vec, v);
        *d = dot(self.vec, [dg.x, dg.y, dg.z]);
        true
    }

    /// (`math_Function::GetStateNumber`'s default.)
    fn get_state_number(&mut self) -> i32 {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Curve;
    use crate::geometry::function_root::function_root;
    use crate::geometry::int_curve_plane::tests::{helix_bspline, tangent_bspline, wavy_bspline, whole};

    /// `math_FunctionRoot(GeomFill_PlanFunc(P, V, C), Guess, Tol, A, B, 50)` against OCCT
    /// 8.1: root, value there, derivative and iteration count. (Searches started at a knot
    /// of the helix, or stalling at a tangency where the derivative is ~1e-9, take other
    /// paths than OCCT's on last-bit differences of the B-spline's derivative: left out.)
    #[test]
    #[allow(clippy::excessive_precision)]
    fn roots_match_occt() {
        let wavy = whole(Curve::BSpline(wavy_bspline()));
        let helix = whole(Curve::BSpline(helix_bspline()));
        let tangent = whole(Curve::BSpline(tangent_bspline()));
        type Case<'a> = (&'a str, &'a dyn CurveAdaptor, [f64; 3], [f64; 3], f64, f64, f64, f64, (f64, f64, f64, usize));
        let cases: [Case; 9] = [
            ("bs_x_from_03", &wavy, [2.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3, 1e-9, -0.5, 1.5, (0.40141929031152657, 0.0, 4.0237614053708723, 4)),
            ("bs_x_from_09", &wavy, [2.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.9, 1e-9, -0.5, 1.5, (0.40141929031152657, 0.0, 4.0237614053708723, 5)),
            ("bs_wavy_05", &wavy, [2.0, 0.0, 0.05], [0.1, 0.05, 1.0], 0.5, 1e-9, -0.5, 1.5, (0.19514912648254937, -3.4694469519536142e-18, 5.0769308392495569, 5)),
            ("bs_wavy_02", &wavy, [2.0, 0.0, 0.05], [0.1, 0.05, 1.0], 0.2, 1e-6, -0.5, 1.5, (0.1951491264825494, 6.2450045135165055e-17, 5.0769308392495542, 3)),
            ("bs_far", &wavy, [100.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1e-9, 0.0, 1.0, (1.0, -94.0, 15.0, 2)),
            ("bs_tight", &wavy, [2.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1e-9, 0.0, 0.2, (0.20000000000000001, -0.95481481481481478, 5.8444444444444441, 3)),
            ("helix_z", &helix, [0.0, 0.0, 1.1], [0.0, 0.0, 1.0], 2.5, 1e-9, -7.0, 21.0, (4.5, 0.0, 0.20000000000000009, 1)),
            ("helix_tilted", &helix, [0.5, 0.5, 1.2], [0.3, -0.2, 2.0], 5.0, 1e-9, -7.0, 21.0, (4.6383166031326208, -1.6653345369377348e-16, 0.88826731396570946, 4)),
            ("tangent_min", &tangent, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.2, 1e-9, -0.5, 1.5, (0.49999999999999989, 0.0, -8.8817841970012523e-16, 1)),
        ];
        let mut deviation: f64 = 0.0;
        for (name, curve, p, v, guess, tol, a, b, (root, value, derivative, iterations)) in cases {
            let mut f = PlanFunc::new(Point::new(p[0], p[1], p[2]), Vector::new(v[0], v[1], v[2]), curve);
            let r = function_root(&mut f, guess, tol, a, b, 50);
            assert!(r.done, "{name}");
            assert_eq!(r.iterations, iterations, "{name}");
            let mut f_root = 0.0;
            f.value(r.root, &mut f_root);
            for (x, y) in [(r.root, root), (f_root, value), (r.derivative, derivative)] {
                deviation = deviation.max((x - y).abs());
            }
        }
        assert!(deviation < 1e-12, "deviation {deviation:e} from OCCT");
        eprintln!("PlanFunc roots: max deviation from OCCT {deviation:e}");
    }

    /// `DEDT` and `D2E` by their formulas (OCCT's `D2E` takes `C(X) - DP` as written).
    #[test]
    fn time_derivatives() {
        let wavy = whole(Curve::BSpline(wavy_bspline()));
        let f = PlanFunc::new(Point::new(2.0, 0.0, 0.05), Vector::new(0.1, 0.05, 1.0), &wavy);
        let (dp, dv) = (Vector::new(0.3, -0.2, 0.1), Vector::new(-0.5, 0.25, 0.75));
        let g = wavy.value(0.4);
        let d = [g.x - 2.0, g.y - 0.0, g.z - 0.05];
        let n = [0.1, 0.05, 1.0];
        let expected = dot([dv.x, dv.y, dv.z], d) - dot(n, [dp.x, dp.y, dp.z]);
        assert!((f.dedt(0.4, dp, dv) - expected).abs() < 1e-15);
        let expected = dot([dv.x, dv.y, dv.z], d) + dot(n, [g.x - dp.x, g.y - dp.y, g.z - dp.z]);
        assert!((f.d2e(0.4, dp, dv) - expected).abs() < 1e-15);
    }
}
