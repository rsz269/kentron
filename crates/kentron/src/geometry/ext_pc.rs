//! OCCT's point-to-curve extrema for any curve (`Extrema_ExtPC`, that is
//! `Extrema_GGExtPC` on its `GeomAbs_OtherCurve` path, with `Extrema_GGenExtPC`,
//! `Extrema_GFuncExtPC` and `math_FunctionRoots`), ported.

use std::f64::consts::PI;

use crate::math::{Point, Vector};

use super::Continuity;
use super::adaptor::CurveAdaptor;
use super::function_roots::{FunctionWithDerivative, function_roots};
use super::tangential_deflection::tangential_deflection;

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = CONFUSION * 0.01;
/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = CONFUSION * CONFUSION;
/// `Precision::Infinite()`.
const INFINITE: f64 = 2.0e100;
/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()`.
const REAL_FIRST: f64 = -f64::MAX;

/// `Precision::IsInfinite`.
fn is_infinite(value: f64) -> bool {
    value.abs() >= 0.5 * INFINITE
}

/// `std::max(a, b)`: `b` only when `a < b`.
fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `std::min(a, b)`: `b` only when `b < a`.
fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// `Epsilon(x)` (Standard_Real): the gap from `x` to the next double away from zero.
fn epsilon(value: f64) -> f64 {
    // `std::nextafter(x, RealLast())` stays at `RealLast()` (`next_up` would overflow).
    if value >= 0.0 {
        let next = if value == REAL_LAST { value } else { value.next_up() };
        next - value
    } else {
        let next = if value == REAL_FIRST { value } else { value.next_down() };
        value - next
    }
}

/// `ElCLib::InPeriod(U, UFirst, ULast)`: `u` shifted by whole periods into
/// `[u_first, u_last)`.
#[doc(alias = "ElCLib::InPeriod")]
fn in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    if is_infinite(u) || is_infinite(u_first) || is_infinite(u_last) {
        // In order to avoid FLT_Overflow exception.
        return u;
    }
    let period = u_last - u_first;
    if period < epsilon(u_last) {
        return u;
    }
    std_max(u_first, u + period * ((u_first - u) / period).ceil())
}

/// One extremum of the distance (`Extrema_POnCurv` with its square distance and kind).
#[derive(Clone, Copy, Debug)]
pub(crate) struct ExtPCPoint {
    pub parameter: f64,
    pub point: Point,
    pub square_distance: f64,
    pub is_min: bool,
}

/// `Extrema_GFuncExtPC::TolFactor`.
const TOL_FACTOR: f64 = 1.0e-12;
/// `Extrema_GFuncExtPC::MinTol`.
const MIN_TOL: f64 = 1.0e-20;
/// `Extrema_GFuncExtPC::MinStep`.
const MIN_STEP: f64 = 1.0e-7;
/// `Extrema_GFuncExtPC::MaxOrder`.
const MAX_ORDER: u32 = 3;

/// The function `F(u) = (C(u) - P).C'(u) / |C'(u)|` whose roots are the extrema of the
/// distance from `P` to the curve `C`, keeping the extrema found (`GetStateNumber`).
///
/// OCCT's `myPinit`/`myCinit` checks (exceptions) are left out: the port is always built
/// with its curve and given its point before evaluation. `myD1Init` is left out: OCCT never
/// reads it.
#[doc(alias = "Extrema_GFuncExtPC")]
#[doc(alias = "Extrema_PCFOfEPCOfExtPC")]
struct DistanceFunction<'a> {
    /// `myC`.
    curve: &'a dyn CurveAdaptor,
    /// `myP`.
    point: Point,
    /// `myU`.
    u: f64,
    /// `myPc`.
    curve_point: Point,
    /// `myD1f`.
    d1f: f64,
    /// `mySqDist`, `myIsMin`, `myPoint`.
    extrema: Vec<ExtPCPoint>,
    /// `myTol`.
    tol: f64,
    /// `myMaxDerivOrder`.
    max_derivative_order: u32,
    /// `myUinfium`.
    u_infimum: f64,
    /// `myUsupremum`.
    u_supremum: f64,
}

impl<'a> DistanceFunction<'a> {
    /// The default constructor followed by `Initialize(C)` for a curve of type
    /// `GeomAbs_OtherCurve`.
    fn new(curve: &'a dyn CurveAdaptor) -> Self {
        let mut function = DistanceFunction {
            curve,
            point: Point::ORIGIN,
            u: 0.0,
            curve_point: Point::ORIGIN,
            d1f: 0.0,
            extrema: Vec::new(),
            tol: MIN_TOL,
            max_derivative_order: 0,
            u_infimum: 0.0,
            u_supremum: 0.0,
        };
        function.sub_interval_initialize(curve.first_parameter(), curve.last_parameter());
        // GeomAbs_OtherCurve.
        function.max_derivative_order = MAX_ORDER;
        function.tol = function.search_of_tolerance();
        function
    }

    /// `SetPoint`.
    fn set_point(&mut self, point: Point) {
        self.point = point;
        self.extrema.clear();
    }

    /// `SubIntervalInitialize`.
    fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.u_infimum = u_first;
        self.u_supremum = u_last;
    }

    /// `SearchOfTolerance`: the tolerance under which the first derivative counts as zero,
    /// `1e-12` of its largest length at eleven samples (at least `1e-20`).
    fn search_of_tolerance(&self) -> f64 {
        const SAMPLE_COUNT: i32 = 10;
        let step = (self.u_supremum - self.u_infimum) / f64::from(SAMPLE_COUNT);
        let mut index = 0;
        let mut max = -INFINITE;
        loop {
            let mut u = self.u_infimum + f64::from(index) * step;
            if u > self.u_supremum {
                u = self.u_supremum;
            }
            let (_, derivative) = self.curve.d1(u);
            // OCCT's `continue` goes to the loop condition.
            if !(is_infinite(derivative.x) || is_infinite(derivative.y)) {
                let length = derivative.length();
                if length > max {
                    max = length;
                }
            }
            index += 1;
            if index > SAMPLE_COUNT {
                break;
            }
        }
        std_max(max * TOL_FACTOR, MIN_TOL)
    }

    /// `aDelta` of the singular cases: `division_factor` of the interval, at least `MinStep`.
    fn singular_delta(&self, division_factor: f64) -> f64 {
        let du = if self.u_supremum >= REAL_LAST || self.u_infimum <= REAL_FIRST { 0.0 } else { self.u_supremum - self.u_infimum };
        std_max(du * division_factor, MIN_STEP)
    }
}

impl FunctionWithDerivative for DistanceFunction<'_> {
    /// `Value`: where the first derivative vanishes (below `myTol`), the direction of the
    /// first non-zero higher derivative (oriented along the curve), else a three-point
    /// difference.
    fn value(&mut self, u: f64, value: &mut f64) -> bool {
        self.u = u;
        let (curve_point, mut d1c) = self.curve.d1(self.u);
        self.curve_point = curve_point;
        if is_infinite(d1c.x) || is_infinite(d1c.y) {
            *value = INFINITE;
            return false;
        }
        let mut ndu = d1c.length();
        if self.max_derivative_order != 0 && ndu <= self.tol {
            // Singular case.
            let delta = self.singular_delta(1.0e-3);
            // The derivative approximated by the Taylor series.
            let mut order = 1;
            let mut derivative;
            let mut derivative_found;
            loop {
                order += 1;
                derivative = self.curve.dn(self.u, order);
                ndu = derivative.length();
                derivative_found = ndu > self.tol;
                if derivative_found || order >= self.max_derivative_order {
                    break;
                }
            }
            if derivative_found {
                let other = if self.u - self.u_infimum < delta { self.u + delta } else { self.u - delta };
                let p1 = self.curve.value(std_min(self.u, other));
                let p2 = self.curve.value(std_max(self.u, other));
                let chord = p2 - p1;
                let direction_factor = derivative.dot(chord);
                d1c = if direction_factor < 0.0 { -derivative } else { derivative };
            } else {
                // The derivative approximated by three points (as vectors from the origin).
                let (p1, p2, p3);
                let parameter_grown;
                if self.u - self.u_infimum < 2.0 * delta {
                    p1 = self.curve.value(self.u);
                    p2 = self.curve.value(self.u + delta);
                    p3 = self.curve.value(self.u + 2.0 * delta);
                    parameter_grown = true;
                } else {
                    p1 = self.curve.value(self.u - 2.0 * delta);
                    p2 = self.curve.value(self.u - delta);
                    p3 = self.curve.value(self.u);
                    parameter_grown = false;
                }
                let (v1, v2, v3) = (p1 - Point::ORIGIN, p2 - Point::ORIGIN, p3 - Point::ORIGIN);
                d1c = if parameter_grown { v1 * -3.0 + v2 * 4.0 - v3 } else { v1 - v2 * 4.0 + v3 * 3.0 };
            }
            ndu = d1c.length();
        }
        if ndu <= MIN_TOL {
            // The first derivative is zero.
            return false;
        }
        let to_curve: Vector = self.curve_point - self.point;
        *value = to_curve.dot(d1c) / ndu;
        true
    }

    /// `Derivative`.
    fn derivative(&mut self, u: f64, derivative: &mut f64) -> bool {
        let mut value = 0.0;
        self.values(u, &mut value, derivative)
    }

    /// `Values`: where the first derivative vanishes, `F'` by a three-point difference of
    /// `F`. Keeps OCCT's restoring of `myPc` to its value before the call (not at `u`)
    /// after that difference.
    fn values(&mut self, u: f64, value: &mut f64, derivative: &mut f64) -> bool {
        let curve_point_old = self.curve_point;
        let point_old = self.point;
        if !self.value(u, value) {
            return false;
        }
        self.u = u;
        self.curve_point = curve_point_old;
        self.point = point_old;
        let (curve_point, d1c, d2c) = self.curve.d2(self.u);
        self.curve_point = curve_point;
        let ndu = d1c.length();
        if ndu <= self.tol {
            // Singular case: the derivative approximated by three points.
            let delta = self.singular_delta(0.01);
            if self.u - self.u_infimum < 2.0 * delta {
                let f1 = *value;
                let u2 = self.u + delta;
                let u3 = self.u + delta * 2.0;
                let mut f2 = 0.0;
                let mut f3 = 0.0;
                if !(self.value(u2, &mut f2) && self.value(u3, &mut f3)) {
                    return false;
                }
                *derivative = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * delta);
            } else {
                let f3 = *value;
                let u1 = self.u - delta * 2.0;
                let u2 = self.u - delta;
                let mut f1 = 0.0;
                let mut f2 = 0.0;
                if !(self.value(u2, &mut f2) && self.value(u1, &mut f1)) {
                    return false;
                }
                *derivative = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * delta);
            }
            self.u = u;
            self.curve_point = curve_point_old;
            self.point = point_old;
        } else {
            let to_curve = self.curve_point - self.point;
            *derivative = ndu + (to_curve.dot(d2c) / ndu) - *value * (d1c.dot(d2c)) / (ndu * ndu);
        }
        self.d1f = *derivative;
        true
    }

    /// `GetStateNumber`: saves the extremum at the last evaluated parameter (square
    /// distance before, point after recomputing `F'`, as OCCT).
    fn get_state_number(&mut self) -> i32 {
        let square_distance = self.curve_point.distance_squared(self.point);
        let mut value = 0.0;
        let mut derivative = 0.0;
        self.values(self.u, &mut value, &mut derivative);
        let is_min = self.d1f > 0.0;
        self.extrema.push(ExtPCPoint { parameter: self.u, point: self.curve_point, square_distance, is_min });
        0
    }
}

/// The extrema of the distance from a point to a curve over a parameter interval, as the
/// roots of [`DistanceFunction`].
#[doc(alias = "Extrema_GGenExtPC")]
#[doc(alias = "Extrema_EPCOfExtPC")]
struct IntervalExtrema<'a> {
    /// `myDone`.
    done: bool,
    /// `mynbsample`.
    sample_count: i32,
    /// `myumin`.
    umin: f64,
    /// `myusup`.
    usup: f64,
    /// `mytolu`.
    tolu: f64,
    /// `mytolF`.
    tolf: f64,
    /// `myF`.
    function: DistanceFunction<'a>,
}

impl<'a> IntervalExtrema<'a> {
    /// The default constructor followed by `Initialize(C)`.
    fn new(curve: &'a dyn CurveAdaptor) -> Self {
        IntervalExtrema { done: false, sample_count: 0, umin: 0.0, usup: 0.0, tolu: 0.0, tolf: 0.0, function: DistanceFunction::new(curve) }
    }

    /// `Initialize(NbU, Umin, Usup, TolU, TolF)`.
    fn initialize(&mut self, sample_count: i32, umin: f64, usup: f64, tolu: f64, tolf: f64) {
        self.sample_count = sample_count;
        self.tolu = tolu;
        self.tolf = tolf;
        self.umin = umin;
        self.usup = usup;
    }

    /// `Perform(P)`.
    fn perform(&mut self, point: Point) {
        self.function.set_point(point);
        self.function.sub_interval_initialize(self.umin, self.usup);
        self.done = false;
        let roots = function_roots(&mut self.function, self.umin, self.usup, self.sample_count, self.tolu, self.tolf, self.tolf, 0.0);
        if !roots.done || roots.all_null {
            return;
        }
        self.done = true;
    }
}

/// `Extrema_CurveTool::DeflCurvIntervals(C)`: parameters splitting a curve of one
/// `C2` interval, from `GCPnts_TangentialDeflection` (just the ends for short, slow or long
/// parameter ranges). The length estimate measures every sample from the first point, as
/// OCCT's (it never advances its start point).
#[doc(alias = "Extrema_CurveTool::DeflCurvIntervals")]
fn defl_curv_intervals(curve: &dyn CurveAdaptor) -> Vec<f64> {
    const EPSD: f64 = 1.0e-3;
    const MAX_DEFLECTION: f64 = 1.0e3;
    const MIN_DEFLECTION: f64 = 1.0e-3;
    const POINT_COUNT: i32 = 23;
    let mut length = 0.0;
    let tf = curve.first_parameter();
    let tl = curve.last_parameter();
    let start = curve.value(tf);
    for i in 2..=POINT_COUNT {
        let t = (tf * f64::from(POINT_COUNT - i) + f64::from(i - 1) * tl) / f64::from(POINT_COUNT - 1);
        let p1 = curve.value(t);
        length += start.distance(p1);
    }
    let dl_dt = length / (tl - tf);
    // `(tl - tf) > 10000.`: to avoid problems with GCPnts.
    if length <= CONFUSION || dl_dt < EPSD || (tl - tf) > 10000.0 {
        return vec![tf, tl];
    }
    let deflection = std_max(0.01 * length / (2.0 * PI), MIN_DEFLECTION);
    if deflection > MAX_DEFLECTION {
        return vec![tf, tl];
    }
    let min_length = std_max(0.00001 * length, CONFUSION);
    let tolerance = std_max(0.00001 * (tl - tf), P_CONFUSION);
    tangential_deflection(curve, PI / 6.0, deflection, 2, tolerance, min_length)
}

/// `RealToInt`: truncation toward zero, saturated.
fn real_to_int(value: f64) -> i32 {
    value as i32
}

/// The state of `Extrema_GGExtPC` for a curve of type `GeomAbs_OtherCurve`.
#[doc(alias = "Extrema_GGExtPC")]
#[doc(alias = "Extrema_ExtPC")]
struct PointCurveExtrema<'a> {
    /// `myC`.
    curve: &'a dyn CurveAdaptor,
    /// `mypoint`, `mySqDist`, `myismin`.
    extrema: Vec<ExtPCPoint>,
    /// `mydone`.
    done: bool,
    /// `mytolu`.
    tolu: f64,
    /// `mytolf`.
    tolf: f64,
    /// `mysample`.
    sample: i32,
    /// `myintuinf`.
    int_uinf: f64,
    /// `myintusup`.
    int_usup: f64,
    /// `myuinf`.
    uinf: f64,
    /// `myusup`.
    usup: f64,
}

impl<'a> PointCurveExtrema<'a> {
    /// `Initialize(C, Uinf, Usup, TolF)`.
    fn initialize(curve: &'a dyn CurveAdaptor, uinf: f64, usup: f64, tolf: f64) -> Self {
        PointCurveExtrema {
            curve,
            extrema: Vec::new(),
            done: false,
            tolu: curve.resolution(CONFUSION),
            tolf,
            sample: 17,
            int_uinf: uinf,
            int_usup: usup,
            uinf,
            usup,
        }
    }

    /// `Perform(P)`: the `default:` (any curve) branch, then the endpoint post-processing of
    /// B-spline, offset and other curves.
    fn perform(&mut self, point: Point) {
        self.extrema.clear();
        self.sample = 17;
        // `Pf`, `Pl` (default-constructed members until set).
        let mut first_point = Point::ORIGIN;
        let mut last_point = Point::ORIGIN;
        let first_distance = if is_infinite(self.uinf) {
            REAL_LAST
        } else {
            first_point = self.curve.value(self.uinf);
            point.distance_squared(first_point)
        };
        let last_distance = if is_infinite(self.usup) {
            REAL_LAST
        } else {
            last_point = self.curve.value(self.usup);
            point.distance_squared(last_point)
        };

        // `default:`
        const MAX_SAMPLE: i32 = 17;
        let mut interval_extrema_done = false;
        let mut bounds = self.curve.intervals(Continuity::C2);
        let mut n = bounds.len() as i32 - 1;
        if n <= 1 {
            bounds = defl_curv_intervals(self.curve);
            n = bounds.len() as i32 - 1;
        }
        self.sample = (self.sample / n).max(MAX_SAMPLE);
        let mut max_interval = 0.0;
        for i in 1..=n as usize {
            let dt = bounds[i] - bounds[i - 1];
            if max_interval < dt {
                max_interval = dt;
            }
        }
        let is_periodic = self.curve.is_periodic();
        let mut s2 = 0.0;
        let mut interval_extrema = IntervalExtrema::new(self.curve);
        for i in 1..=n {
            self.int_uinf = bounds[i as usize - 1];
            self.int_usup = bounds[i as usize];
            self.sample = real_to_int(f64::from(MAX_SAMPLE) * (self.int_usup - self.int_uinf) / max_interval).max(3);
            let mut inf_to_check = self.int_uinf;
            let mut sup_to_check = self.int_usup;
            if is_periodic {
                let period = self.curve.period();
                inf_to_check = in_period(self.int_uinf, self.uinf, self.uinf + period);
                sup_to_check = self.int_usup + (inf_to_check - self.int_uinf);
            }
            let interval_is_not_valid = (self.uinf > sup_to_check) || (self.usup < inf_to_check);
            if interval_is_not_valid {
                continue;
            }
            if self.uinf >= inf_to_check {
                inf_to_check = self.uinf;
            }
            if self.usup <= sup_to_check {
                sup_to_check = self.usup;
            }
            if (sup_to_check - inf_to_check) <= self.tolu {
                continue;
            }
            if i != 1 {
                let (curve_point, derivative) = self.curve.d1(self.int_uinf);
                let s1 = (curve_point - point).dot(derivative);
                if s1 * s2 < 0.0 {
                    self.extrema.push(ExtPCPoint { parameter: self.int_uinf, point: curve_point, square_distance: curve_point.distance_squared(point), is_min: s1 < 0.0 });
                }
            }
            if i != n {
                let (curve_point, derivative) = self.curve.d1(self.int_usup);
                s2 = (curve_point - point).dot(derivative);
            }
            self.interval_perform(point, &mut interval_extrema);
            interval_extrema_done = interval_extrema_done || self.done;
        }
        self.done = interval_extrema_done;

        // GeomAbs_OtherCurve: the end points at zero distance.
        if first_distance < SQUARE_CONFUSION || last_distance < SQUARE_CONFUSION {
            let mut is_first_added = false;
            let mut is_last_added = false;
            for extremum in &self.extrema {
                let u = extremum.parameter;
                if (u - self.uinf).abs() < self.tolu {
                    is_first_added = true;
                } else if (self.usup - u).abs() < self.tolu {
                    is_last_added = true;
                }
            }
            if !is_first_added && first_distance < SQUARE_CONFUSION {
                self.extrema.insert(0, ExtPCPoint { parameter: self.uinf, point: first_point, square_distance: first_distance, is_min: true });
            }
            if !is_last_added && last_distance < SQUARE_CONFUSION {
                self.extrema.push(ExtPCPoint { parameter: self.usup, point: last_point, square_distance: last_distance, is_min: true });
            }
            self.done = true;
        }
    }

    /// `IntervalPerform(P)`: the extrema over `[myintuinf, myintusup]`, added by
    /// [`Self::add_sol`] when within `[myuinf, myusup]`.
    fn interval_perform(&mut self, point: Point, interval_extrema: &mut IntervalExtrema<'a>) {
        interval_extrema.initialize(self.sample, self.int_uinf, self.int_usup, self.tolu, self.tolf);
        interval_extrema.perform(point);
        self.done = interval_extrema.done;
        if self.done {
            for extremum in &interval_extrema.function.extrema {
                let mut u = extremum.parameter;
                if self.curve.is_periodic() {
                    u = in_period(u, self.uinf, self.uinf + self.curve.period());
                }
                if u >= self.uinf - self.tolu && u <= self.usup + self.tolu {
                    self.add_sol(u, extremum.point, extremum.square_distance, extremum.is_min);
                }
            }
        }
    }

    /// `AddSol`: appends the extremum unless one is already within `mytolu` of it.
    fn add_sol(&mut self, u: f64, point: Point, square_distance: f64, is_min: bool) {
        if self.extrema.iter().any(|extremum| (extremum.parameter - u).abs() <= self.tolu) {
            return;
        }
        self.extrema.push(ExtPCPoint { parameter: u, point, square_distance, is_min });
    }
}

/// `Extrema_ExtPC ext; ext.Initialize(curve, uinf, usup, tolf); ext.Perform(p)` for a curve
/// of type `GeomAbs_OtherCurve`: the extrema in OCCT's order, or `None` where `IsDone()` is
/// false.
pub(crate) fn ext_pc_other(curve: &dyn CurveAdaptor, p: Point, uinf: f64, usup: f64, tolf: f64) -> Option<Vec<ExtPCPoint>> {
    let mut extrema = PointCurveExtrema::initialize(curve, uinf, usup, tolf);
    extrema.perform(p);
    extrema.done.then_some(extrema.extrema)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::BSplineCurve;
    use crate::geometry::parametric::ParametricCurve;

    /// `GeomFill_SnglrFunc` over a B-spline: `ratio * D1 x D2`, with OCCT 8.1's order of
    /// operations.
    struct SingularFunction {
        spine: BSplineCurve,
        ratio: f64,
    }

    impl SingularFunction {
        /// The spine's intervals of continuity `continuity` (`GeomAdaptor_Curve::Intervals`
        /// for a non-periodic B-spline): breaks at the interior knots of multiplicity above
        /// `degree - k` (every interior knot for `CN`), none when the spine is that smooth.
        fn spine_intervals(&self, continuity: Continuity) -> Vec<f64> {
            let degree = self.spine.degree();
            let knots = self.spine.knots();
            let multiplicities = self.spine.multiplicities();
            let order = match continuity {
                Continuity::C0 => return vec![knots[0], knots[knots.len() - 1]],
                Continuity::C1 => 1,
                Continuity::C2 => 2,
                Continuity::C3 => 3,
                Continuity::CN => degree,
                Continuity::G1 | Continuity::G2 => panic!("GeomAdaptor_Curve::NbIntervals(): DomainError"),
            };
            let mut bounds = vec![knots[0]];
            for i in 1..knots.len() - 1 {
                if multiplicities[i] as i64 > degree as i64 - order as i64 {
                    bounds.push(knots[i]);
                }
            }
            bounds.push(knots[knots.len() - 1]);
            bounds
        }
    }

    impl CurveAdaptor for SingularFunction {
        fn first_parameter(&self) -> f64 {
            self.spine.first_parameter()
        }
        fn last_parameter(&self) -> f64 {
            self.spine.last_parameter()
        }
        fn intervals(&self, continuity: Continuity) -> Vec<f64> {
            let spine_continuity = match continuity {
                Continuity::C0 => Continuity::C2,
                Continuity::C1 => Continuity::C3,
                Continuity::C2 | Continuity::C3 | Continuity::CN => Continuity::CN,
                Continuity::G1 | Continuity::G2 => Continuity::C0,
            };
            self.spine_intervals(spine_continuity)
        }
        fn is_periodic(&self) -> bool {
            self.spine.is_periodic()
        }
        fn period(&self) -> f64 {
            self.spine.period().expect("periodic")
        }
        fn value(&self, u: f64) -> Point {
            let (_, d1, d2) = self.spine.d2(u);
            Point::ORIGIN + (d1 * self.ratio).cross(d2)
        }
        fn d1(&self, u: f64) -> (Point, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            let scaled = d1 * self.ratio;
            (Point::ORIGIN + scaled.cross(d2), scaled.cross(d3))
        }
        fn d2(&self, u: f64) -> (Point, Vector, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            let d4 = self.spine.dn(u, 4);
            (Point::ORIGIN + d1.cross(d2) * self.ratio, d1.cross(d3) * self.ratio, (d2.cross(d3) + d1.cross(d4)) * self.ratio)
        }
        fn d3(&self, u: f64) -> (Point, Vector, Vector, Vector) {
            let (_, d1, d2, d3) = self.spine.d3(u);
            let d4 = self.spine.dn(u, 4);
            let d5 = self.spine.dn(u, 5);
            (
                Point::ORIGIN + d1.cross(d2) * self.ratio,
                d1.cross(d3) * self.ratio,
                (d2.cross(d3) + d1.cross(d4)) * self.ratio,
                (d1.cross(d5) + d2.cross(d4) * 2.0) * self.ratio,
            )
        }
        fn dn(&self, u: f64, n: u32) -> Vector {
            match n {
                1 => self.d1(u).1,
                2 => self.d2(u).2,
                3 => self.d3(u).3,
                _ => panic!("GeomFill_SnglrFunc::DN: derivative order greater than 3"),
            }
        }
        fn resolution(&self, r3d: f64) -> f64 {
            // Precision::Parametric.
            r3d * 0.01
        }
    }

    fn singular_function(poles: &[(f64, f64, f64)], knots: &[f64], multiplicities: &[usize], degree: usize) -> SingularFunction {
        let poles = poles.iter().map(|&(x, y, z)| Point::new(x, y, z)).collect();
        // The ratio of the OCP references.
        SingularFunction { spine: BSplineCurve::new(poles, knots.to_vec(), multiplicities.to_vec(), degree), ratio: 0.37 }
    }

    /// Compares with `(parameter, square distance, is_min)` from OCP (OCCT 7.9.3):
    /// `Extrema_ExtPC.Initialize(GeomFill_SnglrFunc(GeomAdaptor_Curve(spine)) with
    /// SetRatio(0.37), uinf, usup, 1e-10)`, `Perform(origin)`.
    fn check(function: &SingularFunction, uinf: f64, usup: f64, expected: &[(f64, f64, bool)]) {
        let extrema = ext_pc_other(function, Point::ORIGIN, uinf, usup, 1e-10).expect("done");
        let found: Vec<_> = extrema.iter().map(|e| (e.parameter, e.square_distance, e.is_min)).collect();
        assert_eq!(extrema.len(), expected.len(), "on [{uinf}, {usup}]: {found:?}");
        for (extremum, &(parameter, square_distance, is_min)) in extrema.iter().zip(expected) {
            assert!((extremum.parameter - parameter).abs() < 1e-9, "on [{uinf}, {usup}]: {found:?}");
            assert!((extremum.square_distance - square_distance).abs() <= 1e-9 * square_distance.max(1.0), "on [{uinf}, {usup}]: {found:?}");
            assert_eq!(extremum.is_min, is_min, "on [{uinf}, {usup}]: {found:?}");
            assert!((extremum.point.distance_squared(Point::ORIGIN) - extremum.square_distance).abs() <= 1e-12 * square_distance.max(1.0));
        }
    }

    #[test]
    fn lifted_s_curve_near_inflection() {
        let function = singular_function(&[(0.0, 0.0, 0.0), (1.0, 1.0, 0.05), (2.0, 1.0, 0.1), (3.0, -1.0, 0.12), (4.0, -1.0, 0.1), (5.0, 0.0, 0.0)], &[0.0, 1.0, 2.0, 3.0], &[4, 1, 1, 4], 3);
        assert_eq!(function.intervals(Continuity::C2).len() - 1, 3);
        check(&function, 0.0, 3.0, &[(1.500172161578981, 0.0013264068143910548, true)]);
        check(&function, 0.8999999999999999, 2.25, &[(1.500172161578981, 0.0013264068143910548, true)]);
    }

    #[test]
    fn planar_s_curve_vanishing() {
        let function = singular_function(&[(0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 0.0), (3.0, -1.0, 0.0), (4.0, -1.0, 0.0), (5.0, 0.0, 0.0)], &[0.0, 1.0, 2.0, 3.0], &[4, 1, 1, 4], 3);
        check(&function, 0.0, 3.0, &[(1.5, 1.5186805020668887e-32, true)]);
        check(&function, 0.8999999999999999, 2.25, &[(1.5, 1.5186805020668887e-32, true)]);
    }

    #[test]
    fn helix_like_spline_without_extrema() {
        let function = singular_function(
            &[(1.0, 0.0, 0.0), (1.0, 1.0, 0.2), (0.0, 1.3, 0.4), (-1.0, 1.0, 0.6), (-1.2, 0.0, 0.8), (-1.0, -1.0, 1.0), (0.0, -1.3, 1.2), (1.0, -1.0, 1.4), (1.2, 0.0, 1.6)],
            &[0.0, 0.5, 1.0, 2.0, 2.5, 2.8, 3.0],
            &[4, 1, 1, 1, 1, 1, 4],
            3,
        );
        assert_eq!(function.intervals(Continuity::C2).len() - 1, 6);
        check(&function, 0.0, 3.0, &[]);
        check(&function, 0.8999999999999999, 2.25, &[]);
    }

    #[test]
    fn quartic_wiggle_minima_and_maxima() {
        let function = singular_function(
            &[(0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (2.0, -2.0, 0.01), (3.0, 2.0, -0.01), (4.0, -2.0, 0.0), (5.0, 1.0, 0.02), (6.0, 0.0, 0.0)],
            &[0.0, 0.4, 1.1, 2.0],
            &[5, 1, 1, 5],
            4,
        );
        check(
            &function,
            0.0,
            2.0,
            &[
                (0.27281424998693254, 0.04494692891695159, true),
                (0.3708547069790087, 883.1722464272625, false),
                (0.6781598772328843, 0.0001285629818813302, true),
                (0.8913568883177169, 34.35745405390091, false),
                (1.1089974317101317, 0.0038790530635522775, true),
                (1.3687704798462108, 46.18397283836634, false),
                (1.6029368141836051, 0.020553439966697562, true),
            ],
        );
        check(
            &function,
            0.6,
            1.5,
            &[
                (0.6781598772328843, 0.0001285629818813302, true),
                (0.8913568883177169, 34.35745405390091, false),
                (1.1089974317101317, 0.0038790530635522775, true),
                (1.3687704798462108, 46.18397283836634, false),
            ],
        );
    }

    #[test]
    fn asymmetric_s_curve_with_extremum_at_a_knot() {
        let function = singular_function(
            &[(0.0, 0.0, 0.0), (0.7, 1.3, 0.02), (2.2, 0.8, 0.05), (2.9, -1.4, 0.03), (4.1, -0.6, -0.04), (5.0, 0.4, 0.0), (5.5, 1.5, 0.1)],
            &[0.0, 0.6, 1.0, 1.7, 3.0],
            &[4, 1, 1, 1, 4],
            3,
        );
        check(&function, 0.0, 3.0, &[(0.8605897808675906, 0.06108254926874029, true), (0.999999999362116, 22.614818927442656, true)]);
        check(&function, 0.8999999999999999, 2.25, &[(0.999999999362116, 22.614818927442656, true)]);
    }

    /// Spines of one span: the singular function has one `C2` interval, so the intervals
    /// come from `DeflCurvIntervals` (`GCPnts_TangentialDeflection`).
    #[test]
    fn single_span_spines_through_deflection_intervals() {
        let function = singular_function(&[(0.0, 0.0, 0.0), (1.0, 1.0, 0.1), (2.0, -1.0, 0.1), (3.0, 0.0, 0.0)], &[0.0, 1.0], &[4, 4], 3);
        assert_eq!(function.intervals(Continuity::C2).len() - 1, 1);
        check(&function, 0.0, 1.0, &[(0.5, 0.5544450000000001, true)]);
        check(&function, 0.3, 0.75, &[(0.5, 0.5544450000000001, true)]);

        let function = singular_function(&[(0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (2.0, -2.0, 0.05), (3.0, 2.0, -0.05), (4.0, -2.0, 0.0), (5.0, 0.0, 0.1)], &[0.0, 2.0], &[6, 6], 5);
        check(
            &function,
            0.0,
            2.0,
            &[
                (0.5528179235468673, 0.01017396083906832, true),
                (0.7414723526606957, 1.4298464802584225, false),
                (0.9993273758742189, 0.003529259100741573, true),
                (1.2598803724218255, 1.46723656045831, false),
                (1.4463787814870708, 0.08993783078053545, true),
            ],
        );
        check(
            &function,
            0.6,
            1.5,
            &[
                (0.7414723526606957, 1.4298464802584225, false),
                (0.9993273758742189, 0.003529259100741573, true),
                (1.2598803724218255, 1.46723656045831, false),
                (1.4463787814870708, 0.08993783078053545, true),
            ],
        );
    }

    #[test]
    fn in_period_shifts_into_the_period() {
        assert_eq!(in_period(7.0, 0.0, 2.0 * PI), 7.0 - 2.0 * PI);
        assert_eq!(in_period(-1.0, 0.0, 2.0), 1.0);
        assert_eq!(in_period(0.5, 0.0, 0.0), 0.5);
    }
}
