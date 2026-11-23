//! Where a curve meets a plane, as a scalar function of the curve's parameter: the plane
//! through `P` with normal `V` (not necessarily unit), and `F(x) = V · (C(x) − P)`.
//!
//! The plane-guided trihedron solves `F(x) = 0` to find where the plane square to the spine
//! meets the guide. The guide may be evaluated beyond its own range (the search there is
//! widened by half the guide's range on each side), so `x` is passed to the curve as it is.

use crate::geometry::adaptor::CurveAdaptor;
use crate::geometry::function_roots::FunctionWithDerivative;
use crate::math::{Point, Vector};

/// `F(x) = V · (C(x) − P)` and its derivative `V · C′(x)`.
pub(crate) struct PlanFunc<'a> {
    p: Point,
    v: Vector,
    curve: &'a dyn CurveAdaptor,
}

impl<'a> PlanFunc<'a> {
    pub(crate) fn new(p: Point, v: Vector, curve: &'a dyn CurveAdaptor) -> Self {
        PlanFunc { p, v, curve }
    }

    /// `V · (c − P)`, coordinate-wise in a fixed order (so that `value` and `values` agree
    /// to the last bit).
    fn offset(&self, c: Point) -> f64 {
        let (dx, dy, dz) = (c.x - self.p.x, c.y - self.p.y, c.z - self.p.z);
        self.v.x * dx + self.v.y * dy + self.v.z * dz
    }

    /// The rate of change of `F(x)` at fixed `x` when `P` moves at `dpnt` and `V` at `dvec`:
    /// `dvec · (C(x) − P) − V · dpnt`.
    pub(crate) fn dedt(&self, x: f64, dpnt: Vector, dvec: Vector) -> f64 {
        let w = self.curve.value(x) - self.p;
        dvec.dot(w) - self.v.dot(dpnt)
    }

    /// `dv · (C(x) − P) + V · (C(x) − dp)` (the second term as written: `C(x) − dp`, a
    /// formula kept from the start and pinned by a test).
    pub(crate) fn d2e(&self, x: f64, dp: Vector, dv: Vector) -> f64 {
        let c = self.curve.value(x);
        let w = c - self.p;
        let q = Vector::new(c.x - dp.x, c.y - dp.y, c.z - dp.z);
        dv.dot(w) + self.v.dot(q)
    }
}

impl FunctionWithDerivative for PlanFunc<'_> {
    fn value(&mut self, x: f64, value: &mut f64) -> bool {
        *value = self.offset(self.curve.value(x));
        true
    }

    fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool {
        let (_, d) = self.curve.d1(x);
        *derivative = self.v.dot(d);
        true
    }

    fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool {
        let (c, d) = self.curve.d1(x);
        *value = self.offset(c);
        *derivative = self.v.dot(d);
        true
    }

    fn get_state_number(&mut self) -> i32 {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::adaptor::GeomAdaptor;
    use crate::geometry::curve::{Curve, Line};
    use crate::geometry::function_root::function_root;
    use crate::geometry::int_curve_plane::tests::{helix_bspline, tangent_bspline, wavy_bspline, whole};
    use crate::geometry::int_curve_plane::{intersect_curve_plane, plane_through};
    use crate::math::Direction;

    /// Roots found from guesses near the crossings that the curve–plane intersection reports
    /// are those crossings.
    #[test]
    fn roots_agree_with_the_intersection() {
        let cases: Vec<(Curve, Point, Vector)> = vec![
            (wavy_bspline(), Point::new(0.0, 0.0, 0.1), Vector::new(0.0, 0.0, 1.0)),
            (wavy_bspline(), Point::new(1.5, 0.0, 0.0), Vector::new(1.0, 0.2, 0.3)),
            (helix_bspline(), Point::new(0.0, 0.0, 1.0), Vector::new(0.0, 0.3, 2.0)),
            (helix_bspline(), Point::new(0.5, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
            (tangent_bspline(), Point::new(0.0, 0.0, 0.05), Vector::new(0.0, 0.0, 1.0)),
        ];
        let mut checked = 0;
        for (curve, p, v) in cases {
            let a = whole(&curve);
            let (f, l) = (a.first_parameter(), a.last_parameter());
            let found = intersect_curve_plane(&a, &plane_through(p, Direction::from_vector(v))).unwrap();
            for q in &found.points {
                let guess = (q.w + 1e-3 * (l - f)).min(l);
                let mut function = PlanFunc::new(p, v, &a);
                let root = function_root(&mut function, guess, 1e-12, f, l, 50);
                assert!(root.done && root.iterations <= 50);
                let mut value = 0.0;
                function.value(root.root, &mut value);
                assert!(value.abs() < 1e-12 * v.length(), "F({}) = {value}", root.root);
                // (A root of a smooth function near a simple crossing is that crossing.)
                assert!((root.root - q.w).abs() < 1e-9, "{} vs {}", root.root, q.w);
                checked += 1;
            }
        }
        assert!(checked >= 8, "{checked}");
    }

    /// A line: `F` is affine, with slope `V·D` and root `V·(P − O)/(V·D)`; evaluation beyond
    /// the range passes through.
    #[test]
    fn a_line_and_its_root() {
        let line = Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::new(1.0, 1.0, 0.0)));
        let a = GeomAdaptor::new(&line, 0.0, 1.0);
        let (p, v) = (Point::new(4.0, -1.0, 2.0), Vector::new(0.5, 2.0, -1.0));
        let mut function = PlanFunc::new(p, v, &a);
        let d = Vector::from(Direction::new(1.0, 1.0, 0.0));
        let slope = v.dot(d);
        let expected = v.dot(p - Point::new(1.0, 2.0, 3.0)) / slope;
        for x in [-7.0, 0.0, 0.25, 30.0] {
            let (mut f, mut df, mut g) = (0.0, 0.0, 0.0);
            assert!(function.values(x, &mut f, &mut df));
            assert!(function.value(x, &mut g));
            assert_eq!(f.to_bits(), g.to_bits());
            assert!((f - slope * (x - expected)).abs() < 1e-12);
            assert!((df - slope).abs() < 1e-14);
            let mut d1 = 0.0;
            function.derivative(x, &mut d1);
            assert_eq!(d1, df);
        }
        assert_eq!(function.get_state_number(), 0);
        let root = function_root(&mut function, 0.5, 1e-12, -100.0, 100.0, 50);
        assert!(root.done && (root.root - expected).abs() < 1e-12);
    }

    /// The two derivatives with respect to the plane's motion, against their formulas.
    #[test]
    fn time_derivatives() {
        let curve = helix_bspline();
        let a = whole(&curve);
        let (p, v) = (Point::new(0.3, -0.2, 0.7), Vector::new(0.4, 1.1, -0.6));
        let function = PlanFunc::new(p, v, &a);
        let (dp, dv) = (Vector::new(0.2, -0.5, 1.5), Vector::new(-1.0, 0.25, 0.75));
        for x in [a.first_parameter(), 0.37 * a.last_parameter(), a.last_parameter()] {
            let c = a.value(x);
            let w = c - p;
            let dedt = dv.x * w.x + dv.y * w.y + dv.z * w.z - (v.x * dp.x + v.y * dp.y + v.z * dp.z);
            assert!((function.dedt(x, dp, dv) - dedt).abs() < 1e-15 * 10.0);
            let q = Vector::new(c.x - dp.x, c.y - dp.y, c.z - dp.z);
            let d2e = dv.x * w.x + dv.y * w.y + dv.z * w.z + (v.x * q.x + v.y * q.y + v.z * q.z);
            assert!((function.d2e(x, dp, dv) - d2e).abs() < 1e-15 * 10.0);
        }
    }
}
