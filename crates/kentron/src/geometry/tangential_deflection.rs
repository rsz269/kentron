//! OCCT's `GCPnts_TangentialDeflection` for any curve (its `GeomAbs_OtherCurve` path), ported.
//!
//! Only the path a curve of type `GeomAbs_OtherCurve` takes is translated: `initialize`
//! sends such a curve to `PerformCurve` (and from there, for a straight curve, to
//! `PerformLinear`); the line, circle, Bezier and B-spline branches of OCCT's type switches
//! (`PerformCircular`, `ArcAngularStep`, the pole-count and degree checks) are not reached
//! and not ported, nor is `AddPoint`, which nothing in this path calls.

use crate::math::{Point, Vector};

use super::Continuity;
use super::adaptor::CurveAdaptor;
use super::brent_minimum::BrentMinimum;
use super::particle_swarm::{ParticleSwarm, std_max, std_min};

/// `Us3`, as OCCT writes it (the double nearest 1/3).
#[allow(clippy::excessive_precision)]
const ONE_THIRD: f64 = 0.3333333333333333333333333333;

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;

/// `Precision::Angular()`.
const ANGULAR: f64 = 1.0e-12;

/// The parameters of `GCPnts_TangentialDeflection(C, AngularDeflection, CurvatureDeflection,
/// MinimumOfPoints, UTol, theMinLen)` over the curve's whole range, in order
/// (`Parameter(1..=NbPoints())`).
#[doc(alias = "GCPnts_TangentialDeflection")]
pub(crate) fn tangential_deflection(curve: &dyn CurveAdaptor, angular_deflection: f64, curvature_deflection: f64, minimum_of_points: usize, u_tol: f64, min_len: f64) -> Vec<f64> {
    let mut sampler = TangentialDeflection::new(curve, curve.first_parameter(), curve.last_parameter(), angular_deflection, curvature_deflection, minimum_of_points, u_tol, min_len);
    sampler.perform_curve();
    sampler.parameters
}

/// The state of a `GCPnts_TangentialDeflection` computation.
struct TangentialDeflection<'a> {
    curve: &'a dyn CurveAdaptor,
    /// `myAngularDeflection`.
    angular_deflection: f64,
    /// `myCurvatureDeflection`.
    curvature_deflection: f64,
    /// `myUTol`.
    u_tol: f64,
    /// `myMinNbPnts`.
    min_nb_points: usize,
    /// `myMinLen`.
    min_len: f64,
    /// `myFirstu`, `myLastU`.
    first_u: f64,
    last_u: f64,
    /// `myParameters`, `myPoints` (OCCT's 1-based sequences, kept 0-based here).
    parameters: Vec<f64>,
    points: Vec<Point>,
}

/// `EstimAngl`: the sine of the angle at `middle` of the polyline `p1`, `middle`, `p2`.
#[doc(alias = "EstimAngl")]
fn estimate_angle(p1: Point, middle: Point, p2: Point) -> f64 {
    let v1 = middle - p1;
    let v2 = p2 - middle;
    let l = v1.length() * v2.length();
    if l > f64::MIN_POSITIVE { v1.cross(v2).length() / l } else { 0.0 }
}

/// `getIntervalIdx`: the 1-based index of the interval of `intervals` (OCCT's 1-based
/// `Intervs`) holding `param`, searched from `previous_index`; `intervals.len()` (OCCT's
/// `Upper()`) where none does.
#[doc(alias = "getIntervalIdx")]
fn interval_index(param: f64, intervals: &[f64], previous_index: usize) -> usize {
    let upper = intervals.len();
    let mut index = previous_index;
    while index < upper {
        if param >= intervals[index - 1] && param <= intervals[index] {
            break;
        }
        index += 1;
    }
    index
}

/// `GCPnts_DistFunction`: minus the square distance from the curve to the chord between
/// two parameters.
#[doc(alias = "GCPnts_DistFunction")]
struct ChordDistance<'a> {
    curve: &'a dyn CurveAdaptor,
    /// `myLin`: its location and unit direction.
    line_location: Point,
    line_direction: Vector,
    u1: f64,
    u2: f64,
}

impl<'a> ChordDistance<'a> {
    fn new(curve: &'a dyn CurveAdaptor, u1: f64, u2: f64) -> Self {
        let p1 = curve.value(u1);
        let mut p2 = curve.value(u2);
        #[allow(clippy::neg_cmp_op_on_partial_ord)] // OCCT's test, NaN included
        if !(p1.distance_squared(p2) > f64::MIN_POSITIVE) {
            // For #28812
            p2 = curve.value(u1 + 0.01 * (u2 - u1));
        }
        ChordDistance { curve, line_location: p1, line_direction: occt_direction(p2 - p1), u1, u2 }
    }

    /// `GCPnts_DistFunction::Value`: `-gp_Lin::SquareDistance(C.Value(X))`, `None` outside
    /// the chord's parameter range.
    fn value(&self, x: f64) -> Option<f64> {
        if x < self.u1 || x > self.u2 {
            return None;
        }
        let v = self.curve.value(x) - self.line_location;
        Some(-v.cross(self.line_direction).length_squared())
    }
}

/// `gp_Dir(gp_XYZ)`: `v` divided by its length (OCCT's fast path for a vector of square
/// length within `gp::Resolution()` of 1 keeps the coordinates, which is the same thing for
/// doubles). Panics on a zero vector, where OCCT raises `Standard_ConstructionError`.
fn occt_direction(v: Vector) -> Vector {
    let square_modulus = v.x * v.x + v.y * v.y + v.z * v.z;
    #[allow(clippy::neg_cmp_op_on_partial_ord)] // OCCT's test, NaN included
    let zero_norm = square_modulus <= f64::MIN_POSITIVE * f64::MIN_POSITIVE;
    assert!(!zero_norm, "gp_Dir() - input vector has zero norm");
    let modulus = square_modulus.sqrt();
    Vector::new(v.x / modulus, v.y / modulus, v.z / modulus)
}

impl<'a> TangentialDeflection<'a> {
    /// `GCPnts_TangentialDeflection::initialize` up to its switch on the curve type.
    #[allow(clippy::too_many_arguments)]
    fn new(curve: &'a dyn CurveAdaptor, first_parameter: f64, last_parameter: f64, angular_deflection: f64, curvature_deflection: f64, minimum_of_points: usize, u_tol: f64, min_len: f64) -> Self {
        assert!(!(curvature_deflection < CONFUSION || angular_deflection < ANGULAR), "GCPnts_TangentialDeflection::Initialize - Zero Deflection");
        let (first_u, last_u) = if first_parameter < last_parameter { (first_parameter, last_parameter) } else { (last_parameter, first_parameter) };
        TangentialDeflection {
            curve,
            angular_deflection,
            curvature_deflection,
            u_tol,
            min_nb_points: minimum_of_points.max(2),
            min_len: std_max(min_len, CONFUSION),
            first_u,
            last_u,
            parameters: Vec::new(),
            points: Vec::new(),
        }
    }

    /// `EvaluateDu`: the point at `u` into `point`, and, where the curvature there is not
    /// negligible, the parameter step for the curvature deflection into `du` with
    /// `not_done` cleared (left as they are otherwise).
    #[doc(alias = "GCPnts_TangentialDeflection::EvaluateDu")]
    fn evaluate_du(&self, u: f64, point: &mut Point, du: &mut f64, not_done: &mut bool) {
        let (p, t, n) = self.curve.d2(u);
        *point = p;
        let lt = t.length();
        const LTOL: f64 = CONFUSION;
        if lt > LTOL && n.length() > LTOL {
            let lc = n.cross(t).length();
            let ln = lc / lt;
            if ln > LTOL {
                *du = (8.0 * std_max(self.curvature_deflection, self.min_len) / ln).sqrt();
                *not_done = false;
            }
        }
    }

    /// `PerformLinear`.
    #[doc(alias = "GCPnts_TangentialDeflection::PerformLinear")]
    fn perform_linear(&mut self) {
        let p = self.curve.value(self.first_u);
        self.parameters.push(self.first_u);
        self.points.push(p);
        if self.min_nb_points > 2 {
            let du = (self.last_u - self.first_u) / self.min_nb_points as f64;
            let mut u = self.first_u + du;
            for _ in 2..self.min_nb_points {
                let p = self.curve.value(u);
                self.parameters.push(u);
                self.points.push(p);
                u += du;
            }
        }
        let p = self.curve.value(self.last_u);
        self.parameters.push(self.last_u);
        self.points.push(p);
    }

    /// `PerformCurve`.
    #[doc(alias = "GCPnts_TangentialDeflection::PerformCurve")]
    fn perform_curve(&mut self) {
        let curve = self.curve;
        let mut middle_point = Point::ORIGIN;
        let mut current_point = Point::ORIGIN;
        let mut middle_u;
        let mut v1: Vector;
        let mut v2: Vector;
        let mut l1: f64;
        let mut l2: f64;

        let mut u1 = self.first_u;
        const LTOL: f64 = CONFUSION; // protection against zero length
        let mut atol = 1.0e-2 * self.angular_deflection;
        #[allow(clippy::manual_clamp)] // f64::clamp differs for NaN
        if atol > 1.0e-2 {
            atol = 1.0e-2;
        } else if atol < 1.0e-7 {
            atol = 1.0e-7;
        }

        let last_point = curve.value(self.last_u);

        // Initialization of the computation
        let mut not_done = true;
        let mut dusave = (self.last_u - self.first_u) * ONE_THIRD;
        let mut du = dusave;
        self.evaluate_du(u1, &mut current_point, &mut du, &mut not_done);
        self.parameters.push(u1);
        self.points.push(current_point);

        // Used to detect "isLine" current bspline and in Du computation in general handling.
        // OCCT's 1-based `Intervs(i)` is `intervals[i - 1]`.
        let mut intervals = curve.intervals(Continuity::CN);
        let nb_interv = intervals.len() - 1;

        if not_done || du > 5.0 * dusave {
            // It is either a line or a singularity:
            v1 = last_point - current_point;
            l1 = v1.length();
            if l1 > LTOL {
                // If it is a line, we verify by computing minNbPoints (the curve type is
                // GeomAbs_OtherCurve, so no degree raises the count):
                let mut is_line = true;
                let nb_points = if self.min_nb_points > 3 { self.min_nb_points } else { 3 };
                let mut param = 0.0;
                let mut i = 1;
                while i <= nb_interv && is_line {
                    // Avoid usage intervals out of [myFirstu, myLastU].
                    if intervals[i] < self.first_u || intervals[i - 1] > self.last_u {
                        i += 1;
                        continue;
                    }
                    // Fix border points in applicable intervals, to avoid be out of target interval.
                    if intervals[i - 1] < self.first_u && intervals[i] > self.first_u {
                        intervals[i - 1] = self.first_u;
                    }
                    if intervals[i - 1] < self.last_u && intervals[i] > self.last_u {
                        intervals[i] = self.last_u;
                    }
                    let delta = (intervals[i] - intervals[i - 1]) / nb_points as f64;
                    let mut j = 1;
                    while j <= nb_points && is_line {
                        param = intervals[i - 1] + j as f64 * delta;
                        middle_point = curve.value(param);
                        v2 = middle_point - current_point;
                        l2 = v2.length();
                        if l2 > LTOL {
                            let angle = v2.cross(v1).length() / (l1 * l2);
                            is_line = angle < atol;
                        }
                        j += 1;
                    }
                    i += 1;
                }

                if is_line {
                    self.parameters.clear();
                    self.points.clear();
                    self.perform_linear();
                    return;
                } else {
                    // It was a singularity, continue:
                    self.evaluate_du(param, &mut middle_point, &mut du, &mut not_done);
                }
            } else {
                du = (self.last_u - self.first_u) / 2.1;
                middle_u = self.first_u + du;
                middle_point = curve.value(middle_u);
                v1 = middle_point - current_point;
                l1 = v1.length();
                if l1 < LTOL {
                    // L1 < LTol: this is a zero-length curve, computation is finished:
                    // return a segment of 2 points (protection)
                    self.parameters.push(self.last_u);
                    self.points.push(last_point);
                    return;
                }
            }
        }

        if du > dusave {
            du = dusave;
        } else {
            dusave = du;
        }

        if du < self.u_tol {
            du = self.last_u - self.first_u;
            if du < self.u_tol {
                self.parameters.push(self.last_u);
                self.points.push(last_point);
                return;
            }
        }

        // Normal processing for a curve
        let mut more_points = true;
        let mut u2 = self.first_u;
        let angle_max = self.angular_deflection * 0.5; // because we take the midpoint
        // Indexes of intervals of U1 and U2, used to handle non-uniform case.
        let mut interval_indices = [1usize, 1usize];
        let mut is_need_to_check = false;
        let mut previous_point = *self.points.last().unwrap();

        while more_points {
            interval_indices[0] = interval_index(u1, &intervals, interval_indices[0]);
            u2 += du;

            if u2 >= self.last_u {
                // End of curve
                u2 = self.last_u;
                current_point = last_point;
                du = u2 - u1;
                dusave = du;
            } else {
                current_point = curve.value(u2); // Next point
            }

            let mut coef;
            let mut acoef = 0.0;
            let mut fcoef = 0.0;
            let mut too_large = false;
            let mut correction = true;
            let mut too_small = false;

            while correction {
                // Adjustment of Du
                if is_need_to_check {
                    interval_indices[1] = interval_index(u2, &intervals, interval_indices[0]);
                    if interval_indices[1] > interval_indices[0] {
                        // Jump to another polynom.
                        // Set Du to the smallest value and check deflection on it.
                        let k = interval_indices[0];
                        if du > (intervals[k] - intervals[k - 1]) * ONE_THIRD {
                            du = (intervals[k] - intervals[k - 1]) * ONE_THIRD;
                            u2 = u1 + du;
                            if u2 > self.last_u {
                                u2 = self.last_u;
                            }
                            current_point = curve.value(u2);
                        }
                    }
                }
                middle_u = (u1 + u2) * 0.5; // Verification at the midpoint
                middle_point = curve.value(middle_u);

                v1 = current_point - previous_point; // Deflection criterion
                v2 = middle_point - previous_point;
                l1 = v1.length();

                fcoef = if l1 > self.min_len { v1.cross(v2).length() / (l1 * self.curvature_deflection) } else { 0.0 };

                v1 = current_point - middle_point; // Angular criterion
                l1 = v1.length();
                l2 = v2.length();
                if l1 > self.min_len && l2 > self.min_len {
                    let angg = v1.cross(v2).length() / (l1 * l2);
                    acoef = angg / angle_max;
                } else {
                    acoef = 0.0;
                }

                // Keep the most penalizing coefficient
                coef = std_max(acoef, fcoef);

                if is_need_to_check && coef < 0.55 {
                    is_need_to_check = false;
                    du = dusave;
                    u2 = u1 + du;
                    if u2 > self.last_u {
                        u2 = self.last_u;
                    }
                    current_point = curve.value(u2);
                    continue;
                }

                if coef <= 1.0 {
                    if (self.last_u - u2).abs() < self.u_tol {
                        self.parameters.push(self.last_u);
                        self.points.push(last_point);
                        more_points = false;
                        correction = false;
                    } else if coef >= 0.55 || too_large {
                        self.parameters.push(u2);
                        self.points.push(current_point);
                        previous_point = current_point;
                        correction = false;
                        is_need_to_check = true;
                    } else if too_small {
                        correction = false;
                        previous_point = current_point;
                    } else {
                        too_small = true;
                        du += std_min((u2 - u1) * (1.0 - coef), du * ONE_THIRD);

                        u2 = u1 + du;
                        if u2 > self.last_u {
                            u2 = self.last_u;
                        }
                        current_point = curve.value(u2);
                    }
                } else if coef >= 1.5 {
                    if !previous_point.approx_eq(*self.points.last().unwrap(), CONFUSION) {
                        self.parameters.push(u1);
                        self.points.push(previous_point);
                    }
                    u2 = middle_u;
                    du = u2 - u1;
                    current_point = middle_point;
                } else {
                    du *= 0.9;
                    u2 = u1 + du;
                    current_point = curve.value(u2);
                    too_large = true;
                }
            }

            du = u2 - u1;

            if more_points {
                if u1 > self.first_u {
                    if fcoef > acoef {
                        // Deflection is the splitting criterion
                        self.evaluate_du(u2, &mut current_point, &mut du, &mut not_done);
                        if not_done {
                            du += (du - dusave) * (du / dusave);
                            if du > 1.5 * dusave {
                                du = 1.5 * dusave;
                            }
                            if du < 0.75 * dusave {
                                du = 0.75 * dusave;
                            }
                        }
                    } else {
                        // Angle is the splitting criterion
                        du += (du - dusave) * (du / dusave);
                        if du > 1.5 * dusave {
                            du = 1.5 * dusave;
                        }
                        if du < 0.75 * dusave {
                            du = 0.75 * dusave;
                        }
                    }
                }

                if du < self.u_tol {
                    du = self.last_u - u2;
                    if du < self.u_tol {
                        self.parameters.push(self.last_u);
                        self.points.push(last_point);
                        more_points = false;
                    } else if du * ONE_THIRD > self.u_tol {
                        du *= ONE_THIRD;
                    }
                }
                u1 = u2;
                dusave = du;
            }
        }

        // Readjustment of the second to last point (1-based index `i`):
        let i = self.points.len() - 1;
        if i >= 2 {
            middle_u = self.parameters[i - 2];
            middle_u = (self.last_u + middle_u) * 0.5;
            middle_point = curve.value(middle_u);
            self.parameters[i - 1] = middle_u;
            self.points[i - 1] = middle_point;
        }

        // Add points at midpoints of segments if the minimum
        // number of points is not yet reached.
        let mut nbp = self.points.len();
        while nbp < self.min_nb_points {
            let mut i = 2;
            while i <= nbp {
                middle_u = (self.parameters[i - 2] + self.parameters[i - 1]) * 0.5;
                middle_point = curve.value(middle_u);
                self.parameters.insert(i - 1, middle_u);
                self.points.insert(i - 1, middle_point);
                nbp += 1;
                i += 2;
            }
        }

        // Additional check for intervals
        let min_len2 = self.min_len * self.min_len;
        let max_nbp = 10 * nbp;
        let mut i = 1;
        while i < nbp {
            let u1 = self.parameters[i - 1];
            let u2 = self.parameters[i];

            if u2 - u1 <= self.u_tol {
                i += 1;
                continue;
            }

            // Check maximal deflection on interval;
            let (dmax, umax) = self.estimate_deflection(u1, u2);
            let p1 = self.points[i - 1];
            let p2 = self.points[i];
            middle_point = curve.value(umax);
            let amax = estimate_angle(p1, middle_point, p2);
            if (dmax > self.curvature_deflection || amax > angle_max)
                && umax - u1 > self.u_tol
                && u2 - umax > self.u_tol
                && p1.distance_squared(middle_point) > min_len2
                && p2.distance_squared(middle_point) > min_len2
            {
                self.parameters.insert(i, umax);
                self.points.insert(i, middle_point);
                nbp += 1;
                // i must point to first part of split interval: not advanced.
                if nbp > max_nbp {
                    break;
                }
                continue;
            }
            i += 1;
        }
    }

    /// `EstimDefl`: the largest distance from the curve between `u1` and `u2` to their chord
    /// and its parameter, by Brent's method, or, where that fails, a particle swarm polished
    /// by Brent's method again.
    #[doc(alias = "GCPnts_TangentialDeflection::EstimDefl")]
    fn estimate_deflection(&self, u1: f64, u2: f64) -> (f64, f64) {
        let distance = ChordDistance::new(self.curve, u1, u2);
        let mut function = |x: f64| distance.value(x);

        let nb_iter = 100;
        let rel_tol = std_max(1.0e-3, 2.0 * self.u_tol / (u1.abs() + u2.abs()));

        let mut optimizer = BrentMinimum::new(rel_tol, nb_iter, self.u_tol);
        optimizer.perform(&mut function, u1, (u1 + u2) / 2.0, u2);
        if optimizer.is_done() {
            return ((-optimizer.minimum()).sqrt(), optimizer.location());
        }

        self.estimate_deflection_by_swarm(&distance, u1, u2, &mut optimizer)
    }

    /// The second half of `EstimDefl`, where Brent's method over the whole interval failed:
    /// the best point of a particle swarm over the interval, polished by `optimizer` again
    /// over a step around it.
    fn estimate_deflection_by_swarm(&self, distance: &ChordDistance, u1: f64, u2: f64, optimizer: &mut BrentMinimum) -> (f64, f64) {
        let du = self.last_u - self.first_u;
        let step = std_max(0.1 * du, 100.0 * self.u_tol);
        // `RealToInt`: truncation, clamped to the int range (as `as` does).
        let nb_particles = 8.max((32.0 * (u2 - u1) / du) as i32) as usize;
        let mut swarm = ParticleSwarm::new(vec![u1], vec![u2], vec![step], nb_particles);
        let mut function_mv = |x: &[f64]| distance.value(x[0]);
        let (value, t) = swarm.perform(&mut function_mv, &[step], 100);
        let t = t[0];

        let mut function = |x: f64| distance.value(x);
        optimizer.perform(&mut function, std_max(t - step, u1), t, std_min(t + step, u2));
        if optimizer.is_done() {
            return ((-optimizer.minimum()).sqrt(), optimizer.location());
        }

        ((-value).sqrt(), t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::ParametricCurve;
    use crate::geometry::curve::{BSplineCurve, BezierCurve};

    /// `GeomFill_SnglrFunc` of a spine (D1 x D2 and its derivatives) over `first..last`,
    /// with the CN intervals OCCT gives it.
    struct SingularFunction<C: ParametricCurve<Point = Point, Vector = Vector>> {
        spine: C,
        first: f64,
        last: f64,
        intervals: Vec<f64>,
    }

    impl<C: ParametricCurve<Point = Point, Vector = Vector>> SingularFunction<C> {
        /// Over a single-span spine's whole range.
        fn new(spine: C) -> Self {
            let (first, last) = (spine.first_parameter(), spine.last_parameter());
            SingularFunction { spine, first, last, intervals: vec![first, last] }
        }
    }

    impl<C: ParametricCurve<Point = Point, Vector = Vector>> CurveAdaptor for SingularFunction<C> {
        fn first_parameter(&self) -> f64 {
            self.first
        }
        fn last_parameter(&self) -> f64 {
            self.last
        }
        fn intervals(&self, _continuity: Continuity) -> Vec<f64> {
            self.intervals.clone()
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            unreachable!()
        }
        fn value(&self, u: f64) -> Point {
            let (_, d1, d2) = self.spine.d2(u);
            Point::ORIGIN + d1.cross(d2)
        }
        fn d1(&self, u: f64) -> (Point, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            (Point::ORIGIN + d1.cross(d2), d1.cross(d3))
        }
        fn d2(&self, u: f64) -> (Point, Vector, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            let d4 = self.spine.dn(u, 4);
            (Point::ORIGIN + d1.cross(d2), d1.cross(d3), d2.cross(d3) + d1.cross(d4))
        }
        fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            let d4 = self.spine.dn(u, 4);
            let d5 = self.spine.dn(u, 5);
            (Point::ORIGIN + d1.cross(d2), d1.cross(d3), d2.cross(d3) + d1.cross(d4), d1.cross(d5) + d2.cross(d4) * 2.0)
        }
        fn dn(&self, u: f64, n: u32) -> Vector {
            match n {
                1 => self.d1(u).1,
                2 => self.d2(u).2,
                3 => self.d3(u).3,
                _ => unimplemented!(),
            }
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d / 100.0
        }
    }

    fn points(poles: &[[f64; 3]]) -> Vec<Point> {
        poles.iter().map(|p| Point::new(p[0], p[1], p[2])).collect()
    }

    fn bspline(poles: &[[f64; 3]], first: f64, last: f64) -> BSplineCurve {
        let degree = poles.len() - 1;
        BSplineCurve::new(points(poles), vec![first, last], vec![degree + 1, degree + 1], degree)
    }

    /// The arguments `Extrema_CurveTool::DeflCurvIntervals` passes (angular deflection
    /// π/6, two points at least), with OCCT's (cadquery-ocp 7.9.3)
    /// `GCPnts_TangentialDeflection` parameters for the same singular function.
    fn check(curve: &dyn CurveAdaptor, curvature_deflection: f64, min_len: f64, u_tol: f64, expected: &[f64]) {
        check_with_minimum(curve, curvature_deflection, min_len, u_tol, 2, expected);
    }

    fn check_with_minimum(curve: &dyn CurveAdaptor, curvature_deflection: f64, min_len: f64, u_tol: f64, minimum_of_points: usize, expected: &[f64]) {
        let parameters = tangential_deflection(curve, std::f64::consts::PI / 6.0, curvature_deflection, minimum_of_points, u_tol, min_len);
        assert_eq!(parameters.len(), expected.len(), "{parameters:?} vs {expected:?}");
        for (got, want) in parameters.iter().zip(expected) {
            assert!((got - want).abs() < 1.0e-9, "{parameters:?} vs {expected:?}");
        }
    }

    #[test]
    fn planar_s_cubic_is_a_line() {
        let spine = BezierCurve::new(points(&[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, -2.0, 0.0], [3.0, 0.0, 0.0]]));
        check(&SingularFunction::new(spine), 3.95340878640268, 0.024839999999999997, 1e-05, &[0.0, 1.0]);
    }

    #[test]
    fn twisted_cubic() {
        let spine = BezierCurve::new(points(&[[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 1.0], [3.0, 1.0, 1.0]]));
        check(&SingularFunction::new(spine), 1.6139723785537097, 0.010140887535122359, 1e-05, &[0.0, 1.0]);
    }

    #[test]
    fn helix_quintic() {
        let spine = BezierCurve::new(points(&[[1.0, 0.0, 0.0], [1.0, 0.8, 0.2], [0.2, 1.4, 0.5], [-1.0, 1.2, 0.8], [-1.3, 0.1, 1.1], [-0.8, -0.9, 1.4]]));
        check(
            &SingularFunction::new(spine),
            0.9694735455120522,
            0.006091381936860627,
            1e-05,
            &[0.0, 0.17837073295125083, 0.26755609942687625, 0.3650802976679727, 0.4654138349530513, 0.5657473722381299, 0.699525421951568, 0.8497627109757839, 1.0],
        );
    }

    #[test]
    fn near_planar_s_quintic() {
        let spine = BezierCurve::new(points(&[[0.0, 0.0, 0.0], [1.0, 1.5, 0.01], [2.0, 1.0, -0.02], [3.0, -1.0, 0.02], [4.0, -1.5, -0.01], [5.0, 0.0, 0.0]]));
        check(&SingularFunction::new(spine), 7.325928935898821, 0.04603016905148126, 1e-05, &[0.0, 1.0]);
    }

    #[test]
    fn wavy_quintic_bspline() {
        let spine = bspline(&[[0.0, 0.0, 0.0], [1.0, 3.0, 1.0], [2.0, -3.0, -1.0], [3.0, 3.0, 2.0], [4.0, -3.0, 0.0], [5.0, 0.0, 1.0]], 1.5, 4.0);
        check(
            &SingularFunction::new(spine),
            2.2235481044109067,
            0.01397096477944163,
            2.5e-05,
            &[
                1.5,
                1.8770428491165752,
                2.322851650612239,
                2.3734246606147984,
                2.398711165616078,
                2.4113544181167175,
                2.4239976706173576,
                2.4502774312372493,
                2.476557191857141,
                2.5341248356031447,
                2.9453309216135377,
                3.0598576870326974,
                3.079184078697181,
                3.098204147053782,
                3.1213135301070523,
                3.560656765053526,
                4.0,
            ],
        );
    }

    #[test]
    fn twisted_cubic_bspline() {
        let spine = bspline(&[[0.0, 0.0, 0.0], [2.0, 1.0, -1.0], [1.0, 3.0, 2.0], [4.0, 2.0, 0.0]], -2.0, 3.0);
        check(&SingularFunction::new(spine), 0.048663755332934294, 0.00030576339250007505, 5e-05, &[-2.0, 0.22222222222222188, 1.611111111111111, 3.0]);
    }

    /// A multi-span B-spline spine over `range`, with the singular function's CN `intervals`
    /// as OCCT 8.1 gives them and OCCT 8.1's parameters for the deflection arguments
    /// `(curvature deflection, min length, u tolerance, minimum of points)`.
    #[allow(clippy::too_many_arguments)]
    fn check_multi_span(poles: &[[f64; 3]], knots: &[f64], multiplicities: &[usize], degree: usize, range: (f64, f64), intervals: &[f64], arguments: (f64, f64, f64, usize), expected: &[f64]) {
        let spine = BSplineCurve::new(points(poles), knots.to_vec(), multiplicities.to_vec(), degree);
        let curve = SingularFunction { spine, first: range.0, last: range.1, intervals: intervals.to_vec() };
        let (curvature_deflection, min_len, u_tol, minimum_of_points) = arguments;
        check_with_minimum(&curve, curvature_deflection, min_len, u_tol, minimum_of_points, expected);
    }

    /// Three-span quartic, nine points at least: jumps between polynomials, midpoints added, splits.
    #[test]
    fn three_span_quartic_nine_points() {
        check_multi_span(
            &[[-0.412711653700475, 0.7778398390937324, -0.0034132620107640045], [1.0067802060415274, 0.7234892808564224, 0.008277609672392399], [2.0139404653944633, 1.2501206253293784, -0.00642182402102649], [2.7416832443245096, -1.661715165489479, 0.008738243451305831], [4.1094431952569455, -0.11069194405554095, -7.754726016749418e-05], [4.652921311122854, -0.7369276743480895, 0.0046390265546195925], [6.389525826891221, 1.741038953526548, 0.009263851126104278]],
            &[0.0, 0.6444440200447117, 1.1634677522949326, 1.5373851690910998],
            &[5, 1, 1, 5],
            4,
            (0.0, 1.5373851690910998),
            &[0.0, 0.6444440200447117, 1.1634677522949326, 1.5373851690910998],
            (4.759746707679074, 0.029906370579585572, 1.5373851690911e-05, 9),
            &[0.0, 0.19217314613638747, 0.2810257642799804, 0.38434629227277495, 0.5765194384091624, 0.7686925845455499, 0.9608657306819374, 1.1530388768183248, 1.3452120229547124, 1.5373851690910998],
        );
    }

    /// Three-span quartic, seventeen points at least: all the step corrections.
    #[test]
    fn three_span_quartic_seventeen_points() {
        check_multi_span(
            &[[-0.23505173517535216, 0.8318475068774345, 1.4495949894895084], [0.704506650656431, 0.32752595455972733, -0.43340939019426594], [1.7827317660207163, 1.127888744922695, 1.9946930288463527], [2.8619764919854016, 2.901307597816132, -1.6123336091668694], [4.185695952891123, 1.6714837226545822, 1.2437436370943287], [5.442934671559797, 1.0405930680660114, 1.01904403974439], [5.7731516738661375, -0.1732203298296473, -0.5945146427751946]],
            &[-0.4345402770194111, 0.6635725616641295, 2.6125439970329216, 3.9830929514040534],
            &[5, 1, 1, 5],
            4,
            (-0.4345402770194111, 3.9830929514040534),
            &[-0.4345402770194111, 0.6635725616641295, 2.6125439970329216, 3.9830929514040534],
            (2.0693203932826165, 0.013001923490920421, 4.417633228423465e-05, 17),
            &[-0.4345402770194111, -0.25320038622849006, -0.07186049543756895, 0.10947939535335213, 0.2908192861442732, 0.41171254667155394, 0.47215917693519427, 0.5326058071988347, 0.5449992553600767, 0.5573927035213186, 0.5688681184854316, 0.5803435334495446, 0.6185949166632545, 0.6568462998769644, 0.7333490663043842, 1.2178665870113763, 1.7023841077183681, 1.7669333168915036, 1.8314825260646388, 1.8637571306512064, 1.896031735237774, 2.297529561318788, 2.699027387399802, 2.896594430168617, 3.094161472937432, 3.316394342554087, 3.5386272121707427, 3.6497436469790703, 3.7608600817873983, 3.871976516595726, 3.9830929514040534],
        );
    }

    /// Three-span cubic starting at a singularity.
    #[test]
    fn three_span_cubic_from_a_singularity() {
        check_multi_span(
            &[[-0.26098639611236185, -0.3502500308176191, -0.003112620777590034], [1.0648493173385951, -0.08403984153263444, 0.005181125179746605], [1.5417386780657494, 0.0673353008517128, 0.003128611662670493], [3.0147343815735113, -0.3958827400923379, -0.005883868869879707], [4.178710685885465, 0.7134342242759177, 0.003928667522195221], [4.580582703575115, -1.4289456595696914, -0.0002498255933535351]],
            &[0.1005449931205149, 1.6941232562168684, 2.5033485495342527, 4.044901643877793],
            &[4, 1, 1, 4],
            3,
            (0.1005449931205149, 4.044901643877793),
            &[0.1005449931205149, 1.6941232562168684, 2.5033485495342527, 4.044901643877793],
            (0.07589549353285013, 0.0004768654498467473, 3.944356650757279e-05, 17),
            &[0.1005449931205149, 0.3470672837928448, 0.5935895744651747, 0.8401118651375046, 1.0866341558098345, 1.3331564464821644, 1.3639717328162055, 1.3947870191502467, 1.4564175918183293, 1.5796787371544943, 1.8262010278268241, 2.072723318499154, 2.319245609171484, 2.565767899843814, 2.8122901905161437, 3.0588124811884736, 3.3053347718608035, 3.5518570625331334, 3.7983793532054633, 4.044901643877793],
        );
    }

    /// Three-span quartic starting at a singularity.
    #[test]
    fn three_span_quartic_from_a_singularity() {
        check_multi_span(
            &[[0.03974907845094788, 1.1811270516608463, -0.00764652428404812], [0.7810947040023741, 1.7282422356297666, 0.0013933308553280774], [2.438138572977868, 0.9247509612272262, -0.00660345594171565], [3.3202953048542816, 1.3017520512334095, 0.0015552988238679583], [3.591524331516508, -1.9722610913033005, 0.0073169251400770855], [5.3666318059866205, -0.32240390528612295, -0.005566957007099986], [5.522643510250447, -1.7377320000518806, -0.0025455086254815495]],
            &[-1.0761213468540838, -0.2358310329340776, 0.7744882572297322, 1.737307790733266],
            &[5, 1, 1, 5],
            4,
            (-1.0761213468540838, 1.737307790733266),
            &[-1.0761213468540838, -0.2358310329340776, 0.7744882572297322, 1.737307790733266],
            (1.660331510223194, 0.010432170550081665, 2.81342913758735e-05, 5),
            &[-1.0761213468540838, -0.37276406245724636, 0.33059322193959106, 1.0339505063364285, 1.737307790733266],
        );
    }

    /// Two-span cubic whose singular function is straight.
    #[test]
    fn two_span_cubic_with_straight_singular_function() {
        check_multi_span(
            &[[-0.09995099244303396, -0.39697108595480834, -0.0076785802654261714], [1.407846650899622, -0.6087691105966941, -0.0017786957228815379], [1.9815181269217415, 1.4636799472711637, 0.004803689902083557], [2.5969486358561085, 0.9763916705707336, 0.005505556719493647], [4.17079300207153, -1.4073057783290666, -0.003203793781149644]],
            &[1.9357283006415593, 2.4389260044686654, 4.232327416589332],
            &[4, 1, 4],
            3,
            (1.9357283006415593, 4.232327416589332),
            &[1.9357283006415593, 2.4389260044686654, 4.232327416589332],
            (3.328766070332656, 0.020915254064152074, 2.2965991159477726e-05, 9),
            &[1.9357283006415593, 2.190905980191312, 2.4460836597410642, 2.7012613392908165, 2.9564390188405687, 3.211616698390321, 3.466794377940073, 3.7219720574898254, 4.232327416589332],
        );
    }

    /// Three-span sextic (OCCT gives one CN interval).
    #[test]
    fn three_span_sextic() {
        check_multi_span(
            &[[-0.32200633874940343, 0.25777448365901146, -1.6461781747147994], [1.036197337460754, 1.8058624410679829, 0.9822612423851029], [1.8063993169945118, 1.5453137028464243, 1.7235068476206563], [2.6238869767970634, 1.9454292087127842, 1.9661927274653253], [4.476282235020886, 13.359161788852902, -1.8646906910468632], [4.767083971397948, 0.3462045455686753, -0.9412180886500257], [6.303361295938601, -1.9572220394460538, 0.7941315209853053], [7.15371598797266, 0.5690846889461403, -0.4574378691669989], [7.898989576435468, -0.00897563050190442, -0.37605996661476837]],
            &[1.4754303328668357, 1.935671389340006, 3.001168207209597, 3.319721831858709],
            &[7, 1, 1, 7],
            6,
            (1.4754303328668357, 3.319721831858709),
            &[1.4754303328668357, 3.319721831858709],
            (54.57128529375701, 0.34288149795163947, 1.8442914989918733e-05, 2),
            &[1.4754303328668357, 1.6357132902203413, 1.6757840295587179, 1.7058370840625001, 1.7423515452845957, 1.7867166156694416, 1.8199904184580762, 1.8477185874486048, 1.876245510278367, 2.1115240627421823, 2.3778439033683547, 2.404922193398168, 2.4252309109205283, 2.445217735270159, 2.480623795000799, 2.5491657166387642, 2.9832645536792106, 3.060374215521921, 3.0820613079151835, 3.088648762229637, 3.095131805393206, 3.1030087028369433, 3.2113652673478263, 3.319721831858709],
        );
    }

    /// Three-span quintic trimmed inside its knots.
    #[test]
    fn trimmed_three_span_quintic() {
        check_multi_span(
            &[[-0.30794928544757183, 0.717442230898405, -0.4766013359012442], [1.1044228521085468, 0.09200430844337593, -1.13206350616404], [1.7833106333515414, 0.17234339696640788, -1.0999863966074135], [2.683826854841059, 1.2183779246895168, -1.356157105017838], [4.064115589602215, -0.21826779200188717, -0.6283021721029178], [5.362008987882685, -1.2457604477315019, -1.7647038905622225], [5.658694758900803, -0.10794010606544324, 1.762133201520819], [6.79871240508376, 0.9074670997746619, 0.43947917497081823], [8.00929102623023, 1.4873097470669454, 1.411866385365684]],
            &[0.0, 0.3230358282079954, 0.7633671017415435, 1.7301231524615805],
            &[6, 1, 2, 6],
            5,
            (0.7004309757825861, 1.2533722206983544),
            &[0.7004309757825861, 0.7633671017415435, 1.2533722206983544],
            (1.8433077582022335, 0.011581844222946416, 5.529412449157684e-06, 2),
            &[0.7004309757825861, 0.7454077456849562, 0.8053767722214497, 0.8681581057425256, 0.8917011058129292, 0.9045732411014222, 0.9178161786410242, 0.9393486310486772, 0.96088108345633, 1.0058807440191946, 1.0855941996696892, 1.1706582228066549, 1.2533722206983544],
        );
    }

    /// The particle swarm half of `EstimDefl` (which Brent's method failing over a whole
    /// interval leads to, unseen in the curves above), against OCCT 8.1's
    /// `GCPnts_DistFunction`, `math_PSO` and `math_BrentMinimum` run the same way:
    /// `(u tolerance, u1, u2, deflection, its parameter)`.
    fn check_swarm_deflection(curve: &dyn CurveAdaptor, cases: &[(f64, f64, f64, f64, f64)]) {
        for &(u_tol, u1, u2, deflection, parameter) in cases {
            let sampler = TangentialDeflection::new(curve, curve.first_parameter(), curve.last_parameter(), std::f64::consts::PI / 6.0, 1.0, 2, u_tol, 1.0e-3);
            let distance = ChordDistance::new(curve, u1, u2);
            let mut optimizer = BrentMinimum::new(std_max(1.0e-3, 2.0 * u_tol / (u1.abs() + u2.abs())), 100, u_tol);
            let (got_deflection, got_parameter) = sampler.estimate_deflection_by_swarm(&distance, u1, u2, &mut optimizer);
            assert!((got_deflection - deflection).abs() < 1.0e-12 * deflection.max(1.0), "{u1} {u2}: {got_deflection} vs {deflection}");
            assert!((got_parameter - parameter).abs() < 1.0e-9, "{u1} {u2}: {got_parameter} vs {parameter}");
        }
    }

    #[test]
    fn swarm_deflection_helix_quintic() {
        let spine = BezierCurve::new(points(&[[1.0, 0.0, 0.0], [1.0, 0.8, 0.2], [0.2, 1.4, 0.5], [-1.0, 1.2, 0.8], [-1.3, 0.1, 1.1], [-0.8, -0.9, 1.4]]));
        check_swarm_deflection(
            &SingularFunction::new(spine),
            &[
                (1e-05, 0.0, 1.0, 25.891143345683847, 0.5186395637014166),
                (1e-05, 0.0, 0.5, 8.525612982708306, 0.24735681944220464),
                (1e-05, 0.2, 0.9, 12.355407426199166, 0.5664303310962348),
                (1e-05, 0.6, 0.61, 0.0025240246876251544, 0.6050124120855797),
                (1e-05, 0.17837073295125083, 0.26755609942687625, 0.27608162294969585, 0.2229064685191759),
            ],
        );
    }

    #[test]
    fn swarm_deflection_wavy_quintic_bspline() {
        let spine = bspline(&[[0.0, 0.0, 0.0], [1.0, 3.0, 1.0], [2.0, -3.0, -1.0], [3.0, 3.0, 2.0], [4.0, -3.0, 0.0], [5.0, 0.0, 1.0]], 1.5, 4.0);
        check_swarm_deflection(
            &SingularFunction::new(spine),
            &[
                (2.5e-05, 1.5, 4.0, 2.319666072994962, 3.7641169982125517),
                // Not (2.0, 3.5): the curve is as far from the chord at both ends of the swarm's
                // box, OCCT's evaluation putting the start ahead by 2e-17 and kentron's the end,
                // so the swarms part.
                (2.5e-05, 1.7, 3.3, 4.066600609027816, 2.2918833380620036),
                (2.5e-05, 2.1, 3.9, 3.694734071132065, 3.2317778789786384),
                (2.5e-05, 2.322851650612239, 2.9453309216135377, 0.5308432933241741, 2.5350496096526083),
                (2.5e-05, 3.0, 3.01, 0.00034043879987279396, 3.005037673667905),
            ],
        );
    }
}
