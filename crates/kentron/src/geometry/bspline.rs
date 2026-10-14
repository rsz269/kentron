//! B-spline (NURBS) curves (2D and 3D) and surfaces.
//!
//! A B-spline is given by its degree, poles, optional weights, and a knot vector stored as
//! distinct increasing knots with multiplicities. Indices (poles, knots) start at 0. A
//! periodic B-spline wraps around: its pole and knot lists describe one period.

use crate::error::{Error, Result};
use crate::math::{
    AffineTransform, AffineTransformable, Point, Point2, Transform, Transform2, Transformable,
    Transformable2,
};

use super::control::ControlPoint;
use super::{Continuity, ParametricCurve};

/// How knots are spaced.
#[doc(alias = "GeomAbs_BSplKnotDistribution")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KnotDistribution {
    NonUniform,
    /// Equally spaced knots, all of multiplicity 1.
    Uniform,
    /// Equally spaced knots of multiplicity 1, except at the ends.
    QuasiUniform,
    /// Every knot has full multiplicity: the curve is a chain of Bezier pieces.
    PiecewiseBezier,
}

/// A B-spline curve. The type parameter is the point type: [`Point`] (3D, the default) or
/// [`Point2`] ([`BSplineCurve2`](super::curve2::BSplineCurve2)).
#[doc(alias = "Geom_BSplineCurve")]
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineCurve<P = Point> {
    degree: usize,
    poles: Vec<P>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    multiplicities: Vec<usize>,
    periodic: bool,
    /// The knots repeated by multiplicity, kept for evaluation.
    flat: Vec<f64>,
    /// The curve as Bézier segments, for evaluating it (reset when it changes).
    segments: super::patches::SegmentCache,
}

/// The knots repeated by multiplicity.
/// Flat knot `index` (counting each knot by its multiplicity), without spelling them out.
fn flat_knot(kv: &KnotVector, index: usize) -> f64 {
    let mut seen = 0;
    for (&k, &m) in kv.knots.iter().zip(&kv.multiplicities) {
        seen += m;
        if index < seen {
            return k;
        }
    }
    kv.knots[kv.knots.len() - 1]
}

fn flatten(knots: &[f64], multiplicities: &[usize]) -> Vec<f64> {
    knots.iter().zip(multiplicities).flat_map(|(&k, &m)| std::iter::repeat(k).take(m)).collect()
}

/// The flat knots of a periodic curve over one period, extended by whole periods on both
/// sides so that the `n` poles, repeated cyclically (pole `i mod n`), make an ordinary
/// B-spline of `n + degree` poles valid over the period (as OCCT's `BSplCLib::KnotSequence`).
/// For reversing a periodic B-spline of `n` poles, degree `p` and end multiplicity `m0`:
/// the old index of new pole `j` (the one whose basis function has the reflected support).
fn reversed_pole(n: usize, p: usize, m0: usize) -> impl Fn(usize) -> usize {
    move |j| (2 * n + p - j - m0) % n
}

pub(crate) fn periodic_flatten(knots: &[f64], multiplicities: &[usize], degree: usize) -> Vec<f64> {
    let period = knots[knots.len() - 1] - knots[0];
    // One period's knots (the last knot is the first of the next period).
    let one: Vec<f64> = flatten(&knots[..knots.len() - 1], &multiplicities[..multiplicities.len() - 1]);
    let n = one.len();
    let at = |j: isize| -> f64 {
        let (q, r) = (j.div_euclid(n as isize), j.rem_euclid(n as isize));
        one[r as usize] + q as f64 * period
    };
    // Index 0 of `one` (the first copy of the first knot) sits at position
    // `degree + 1 - m0`, so that position `degree` holds the last copy of the first knot.
    let shift = (degree + 1 - multiplicities[0]) as isize;
    (0..(n + 2 * degree + 1) as isize).map(|i| at(i - shift)).collect()
}

impl<P: ControlPoint> BSplineCurve<P> {
    /// The largest supported degree.
    pub const MAX_DEGREE: usize = 25;

    /// A non-periodic polynomial B-spline. Panics on inconsistent data (see
    /// [`BSplineCurve::try_new`]).
    pub fn new(poles: Vec<P>, knots: Vec<f64>, multiplicities: Vec<usize>, degree: usize) -> Self {
        Self::try_new(poles, knots, multiplicities, degree).expect("invalid B-spline curve")
    }

    /// Fails unless `1 <= degree <= MAX_DEGREE`, knots are strictly increasing,
    /// multiplicities are at most `degree` inside and `degree + 1` at the ends, and the
    /// pole count is `sum(multiplicities) - degree - 1`.
    pub fn try_new(
        poles: Vec<P>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
        degree: usize,
    ) -> Result<Self> {
        let invalid = |msg: String| Err(Error::InvalidInput(msg));
        if !(1..=Self::MAX_DEGREE).contains(&degree) {
            return invalid(format!("B-spline degree {degree} out of range"));
        }
        if knots.len() < 2 || knots.len() != multiplicities.len() {
            return invalid("a B-spline needs at least two knots, each with a multiplicity".into());
        }
        if knots.windows(2).any(|w| !(w[0] < w[1])) {
            return invalid("B-spline knots must be strictly increasing".into());
        }
        let last = multiplicities.len() - 1;
        let bad_multiplicity = multiplicities.iter().enumerate().any(|(i, &m)| {
            let max = if i == 0 || i == last { degree + 1 } else { degree };
            m == 0 || m > max
        });
        if bad_multiplicity {
            return invalid(format!("invalid B-spline multiplicities {multiplicities:?} for degree {degree}"));
        }
        let total: usize = multiplicities.iter().sum();
        if total < degree + 2 || poles.len() != total - degree - 1 {
            return invalid(format!("a B-spline of degree {degree} with these knots needs {} poles, got {}", total.saturating_sub(degree + 1), poles.len()));
        }
        let flat = flatten(&knots, &multiplicities);
        Ok(Self { degree, poles, weights: None, knots, multiplicities, periodic: false, flat, segments: Default::default() })
    }

    /// A periodic polynomial B-spline: the pole count is `sum(multiplicities)` minus the
    /// last knot's multiplicity, and the curve wraps from the last pole to the first.
    pub fn new_periodic(poles: Vec<P>, knots: Vec<f64>, multiplicities: Vec<usize>, degree: usize) -> Self {
        Self::try_new_periodic(poles, knots, multiplicities, degree).expect("invalid periodic B-spline curve")
    }

    pub fn try_new_periodic(
        poles: Vec<P>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
        degree: usize,
    ) -> Result<Self> {
        let invalid = |msg: String| Err(Error::InvalidInput(msg));
        if !(1..=Self::MAX_DEGREE).contains(&degree) {
            return invalid(format!("B-spline degree {degree} out of range"));
        }
        if knots.len() < 2 || knots.len() != multiplicities.len() {
            return invalid("a B-spline needs at least two knots, each with a multiplicity".into());
        }
        if knots.windows(2).any(|w| !(w[0] < w[1])) {
            return invalid("B-spline knots must be strictly increasing".into());
        }
        let last = multiplicities.len() - 1;
        if multiplicities.iter().any(|&m| m == 0 || m > degree) || multiplicities[0] != multiplicities[last] {
            return invalid(format!("invalid periodic B-spline multiplicities {multiplicities:?} for degree {degree}"));
        }
        let count: usize = multiplicities[..last].iter().sum();
        if poles.len() != count || count < 2 {
            return invalid(format!("a periodic B-spline with these knots needs {count} poles, got {}", poles.len()));
        }
        let flat = periodic_flatten(&knots, &multiplicities, degree);
        Ok(Self { degree, poles, weights: None, knots, multiplicities, periodic: true, flat, segments: Default::default() })
    }

    /// Makes the curve rational with the given weights. Panics unless there is one
    /// positive weight per pole.
    pub fn with_weights(mut self, weights: Vec<f64>) -> Self {
        assert!(
            weights.len() == self.poles.len() && weights.iter().all(|&w| w > 0.0),
            "invalid B-spline weights"
        );
        self.weights = Some(weights);
        self.segments = Default::default();
        self
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn poles(&self) -> &[P] {
        &self.poles
    }

    pub fn pole(&self, i: usize) -> P {
        self.poles[i]
    }

    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// The weight of pole `i` (1 for a polynomial curve).
    pub fn weight(&self, i: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |w| w[i])
    }

    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// The distinct knots, increasing.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    pub fn knot(&self, i: usize) -> f64 {
        self.knots[i]
    }

    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }

    pub fn multiplicity(&self, i: usize) -> usize {
        self.multiplicities[i]
    }

    /// The knots repeated by multiplicity (the "flat" knot vector).
    #[doc(alias = "KnotSequence")]
    pub fn flat_knots(&self) -> Vec<f64> {
        self.flat.clone()
    }

    /// The flat knots evaluation uses (extended by a period on both sides for a periodic
    /// curve).
    fn evaluation_knots(&self) -> &[f64] {
        &self.flat
    }

    /// Derivatives of orders `0..=order` at `t`, as coordinates: outside the parameter
    /// range those of the end spans' polynomials (like OCCT), wrapped round for a periodic
    /// curve.
    fn derivatives(&self, t: f64, order: usize) -> Vec<[f64; 3]> {
        let flat = self.evaluation_knots();
        let p = self.degree;
        // A periodic curve's poles repeat: its first `p` poles again after the last.
        let poles = self.poles.len();
        let n = if self.periodic { poles + p } else { poles };
        let (lo, hi) = (flat[p], flat[n]);
        let t = if self.periodic && !(lo..=hi).contains(&t) { lo + (t - lo).rem_euclid(hi - lo) } else { t };
        // The span: flat[span] <= t < flat[span + 1], with p <= span < n.
        let span = (flat[..n].partition_point(|&k| k <= t)).saturating_sub(1).clamp(p, n - 1);
        let basis = basis_derivatives(flat, span, t, p, order.min(p));
        // Derivatives of the homogeneous curve (w·P, w).
        let mut homogeneous = vec![[0.0; 4]; order + 1];
        for (k, row) in basis.iter().enumerate() {
            for (j, &b) in row.iter().enumerate() {
                let i = (span - p + j) % poles;
                let w = self.weight(i);
                let c = self.poles[i].coords();
                for d in 0..3 {
                    homogeneous[k][d] += b * w * c[d];
                }
                homogeneous[k][3] += b * w;
            }
        }
        // Back to Cartesian: C(k) = (A(k) - Σ binom(k, i)·w(i)·C(k - i)) / w.
        let mut out: Vec<[f64; 3]> = Vec::with_capacity(order + 1);
        for k in 0..=order {
            let mut c = [homogeneous[k][0], homogeneous[k][1], homogeneous[k][2]];
            let mut binom = 1.0;
            for i in 1..=k {
                binom = binom * (k - i + 1) as f64 / i as f64;
                for d in 0..3 {
                    c[d] -= binom * homogeneous[i][3] * out[k - i][d];
                }
            }
            for v in &mut c {
                *v /= homogeneous[0][3];
            }
            out.push(c);
        }
        out
    }

    /// Derivatives of orders `0..=order` (at most 2) as [`Self::derivatives`] gives them, on
    /// the stack.
    fn fast_derivatives<const N: usize>(&self, t: f64, order: usize) -> [[f64; 3]; 3] {
        let flat = self.evaluation_knots();
        let p = self.degree;
        let poles = self.poles.len();
        let n = if self.periodic { poles + p } else { poles };
        let (lo, hi) = (flat[p], flat[n]);
        let t = if self.periodic && !(lo..=hi).contains(&t) { lo + (t - lo).rem_euclid(hi - lo) } else { t };
        let span = (flat[..n].partition_point(|&k| k <= t)).saturating_sub(1).clamp(p, n - 1);
        // Values and first derivatives always the same way, so that point, d1 and d2 agree.
        let (b, d) = basis_and_first::<N>(flat, 0, span, t, p, order >= 1);
        let second = if order >= 2 { basis_to_second::<N>(flat, 0, span, t, p)[2] } else { [0.0; N] };
        let basis = [b, d, second];
        let mut homogeneous = [[0.0f64; 4]; 3];
        for j in 0..=p {
            let i = (span - p + j) % poles;
            let w = self.weight(i);
            let c = self.poles[i].coords();
            let h = [w * c[0], w * c[1], w * c[2], w];
            for k in 0..=order {
                for d in 0..4 {
                    homogeneous[k][d] += basis[k][j] * h[d];
                }
            }
        }
        let w = homogeneous[0][3];
        let mut out = [[0.0f64; 3]; 3];
        for d in 0..3 {
            out[0][d] = homogeneous[0][d] / w;
        }
        if order >= 1 {
            for d in 0..3 {
                out[1][d] = (homogeneous[1][d] - homogeneous[1][3] * out[0][d]) / w;
            }
        }
        if order >= 2 {
            for d in 0..3 {
                out[2][d] = (homogeneous[2][d] - 2.0 * homogeneous[1][3] * out[1][d] - homogeneous[2][3] * out[0][d]) / w;
            }
        }
        out
    }

    /// [`Self::fast_derivatives`] with arrays sized for the degree (they are zeroed on every
    /// evaluation); `None` beyond the largest degree. From the Bézier segments where there
    /// are.
    fn fast_derivatives_by_degree(&self, t: f64, order: usize) -> Option<[[f64; 3]; 3]> {
        if let Some(d) = self.segments().and_then(|s| s.evaluate(t, order)) {
            return Some(d);
        }
        Some(match self.degree {
            0..=3 => self.fast_derivatives::<4>(t, order),
            4..=8 => self.fast_derivatives::<9>(t, order),
            9..=25 => self.fast_derivatives::<26>(t, order),
            _ => return None,
        })
    }

    /// The curve as Bézier segments (made on first use); `None` for one with too many.
    fn segments(&self) -> Option<&super::patches::Segments> {
        self.segments.get(|| {
            let p = self.degree;
            let count = self.flat.len() - p - 1;
            let n = self.poles.len();
            // (A periodic curve's poles repeat cyclically.)
            let points: Vec<[f64; 4]> = (0..count)
                .map(|i| {
                    let w = self.weight(i % n);
                    let c = self.poles[i % n].coords();
                    [w * c[0], w * c[1], w * c[2], w]
                })
                .collect();
            super::patches::Segments::new(&self.flat, p, &points, self.periodic)
        })
    }

    /// The poles as homogeneous rows `w·x, w·y, w·z, w`.
    fn rows(&self) -> Vec<Vec<f64>> {
        (0..self.poles.len())
            .map(|i| {
                let (w, c) = (self.weight(i), self.poles[i].coords());
                vec![w * c[0], w * c[1], w * c[2], w]
            })
            .collect()
    }

    /// Poles and weights (`None` unless `rational`) from homogeneous rows.
    fn from_rows(rows: &[Vec<f64>], rational: bool) -> (Vec<P>, Option<Vec<f64>>) {
        let poles = rows.iter().map(|r| P::from_coords([r[0] / r[3], r[1] / r[3], r[2] / r[3]])).collect();
        (poles, rational.then(|| rows.iter().map(|r| r[3]).collect()))
    }

    /// The (non-periodic) curve as a [`Spline`](super::knots::Spline); a periodic one
    /// unrolled over its period.
    pub(crate) fn spline(&self) -> super::knots::Spline {
        match self.periodic {
            true => {
                let (lo, hi) = (self.knots[0], self.knots[self.knots.len() - 1]);
                self.periodic_form().unroll(lo, hi)
            }
            false => super::knots::Spline { degree: self.degree, flat: self.flat.clone(), poles: self.rows() },
        }
    }

    fn periodic_form(&self) -> super::knots::Periodic {
        super::knots::Periodic {
            degree: self.degree,
            knots: self.knots.clone(),
            multiplicities: self.multiplicities.clone(),
            poles: self.rows(),
        }
    }

    /// Becomes the non-periodic curve of `spline`.
    pub(crate) fn set_spline(&mut self, spline: super::knots::Spline) {
        let (knots, multiplicities) = super::knots::group(&spline.flat);
        let (poles, weights) = Self::from_rows(&spline.poles, self.is_rational());
        *self = Self { degree: spline.degree, poles, weights, knots, multiplicities, periodic: false, flat: spline.flat, segments: Default::default() };
    }

    /// Becomes the periodic curve of `periodic`.
    fn set_periodic_form(&mut self, periodic: super::knots::Periodic) {
        let (poles, weights) = Self::from_rows(&periodic.poles, self.is_rational());
        let flat = periodic_flatten(&periodic.knots, &periodic.multiplicities, periodic.degree);
        *self = Self {
            degree: periodic.degree,
            poles,
            weights,
            knots: periodic.knots,
            multiplicities: periodic.multiplicities,
            periodic: true,
            flat,
            segments: Default::default(),
        };
    }

    /// `u` moved onto a knot within the parametric tolerance of it.
    fn snapped(&self, u: f64) -> f64 {
        self.knots.iter().copied().find(|k| (k - u).abs() <= crate::tolerance::PARAMETRIC).unwrap_or(u)
    }

    /// The same curve, non-periodic and clamped (end knots of full multiplicity).
    pub(crate) fn clamped(&self) -> Self {
        let mut spline = self.spline();
        let (lo, hi) = spline.domain();
        spline.segment(lo, hi);
        let mut out = self.clone();
        out.set_spline(spline);
        out
    }

    /// This clamped curve followed by `next` (clamped, of the same degree, starting where
    /// this one ends): `next`'s parameters continue this one's, the junction is one pole
    /// (the mean of the two) with a knot of multiplicity `degree`.
    pub(crate) fn joined(&self, next: &Self) -> Self {
        let p = self.degree;
        let rational = self.is_rational() || next.is_rational();
        let (a, mut b) = (self.rows(), next.rows());
        // Weights scaled so that the junction poles weigh the same (the curve is unchanged).
        let scale = a[a.len() - 1][3] / b[0][3];
        for row in &mut b {
            row.iter_mut().for_each(|x| *x *= scale);
        }
        let junction: Vec<f64> = a[a.len() - 1].iter().zip(&b[0]).map(|(x, y)| 0.5 * (x + y)).collect();
        let mut rows: Vec<Vec<f64>> = a[..a.len() - 1].to_vec();
        rows.push(junction);
        rows.extend_from_slice(&b[1..]);
        let shift = self.knots[self.knots.len() - 1] - next.knots[0];
        let mut flat: Vec<f64> = self.flat[..self.flat.len() - 1].to_vec();
        flat.extend(next.flat[p + 1..].iter().map(|k| k + shift));
        let mut out = self.clone();
        out.weights = rational.then(|| vec![1.0; rows.len()]);
        out.set_spline(super::knots::Spline { degree: p, flat, poles: rows });
        if !rational {
            out.weights = None;
        }
        out
    }

    #[doc(alias = "KnotDistribution")]
    pub fn knot_distribution(&self) -> KnotDistribution {
        let (p, last) = (self.degree, self.knots.len() - 1);
        let inner = &self.multiplicities[1..last];
        let clamped = !self.periodic && self.multiplicities[0] == p + 1 && self.multiplicities[last] == p + 1;
        if clamped && inner.iter().all(|&m| m == p) {
            return KnotDistribution::PiecewiseBezier;
        }
        let step = (self.knots[last] - self.knots[0]) / last as f64;
        let even = self.knots.windows(2).all(|w| ((w[1] - w[0]) - step).abs() <= 1e-9 * step.abs().max(1.0));
        match (even, inner.iter().all(|&m| m == 1)) {
            (true, true) if self.multiplicities[0] == 1 && self.multiplicities[last] == 1 => KnotDistribution::Uniform,
            (true, true) if clamped => KnotDistribution::QuasiUniform,
            _ => KnotDistribution::NonUniform,
        }
    }

    pub fn set_pole(&mut self, i: usize, p: P) {
        self.poles[i] = p;
        self.segments = Default::default();
    }

    pub fn set_pole_with_weight(&mut self, i: usize, p: P, weight: f64) {
        self.set_pole(i, p);
        self.set_weight(i, weight);
    }

    /// Sets the weight of pole `i`, making the curve rational. Panics if `weight <= 0`.
    pub fn set_weight(&mut self, i: usize, weight: f64) {
        assert!(weight > 0.0, "non-positive weight {weight}");
        let n = self.poles.len();
        self.weights.get_or_insert_with(|| vec![1.0; n])[i] = weight;
        self.segments = Default::default();
    }

    /// Moves knot `i` to `value`; panics if that breaks the knots' increasing order.
    pub fn set_knot(&mut self, i: usize, value: f64) {
        let last = self.knots.len() - 1;
        assert!(
            (i == 0 || self.knots[i - 1] < value) && (i == last || value < self.knots[i + 1]),
            "knot {i} moved to {value} would break the knots' order"
        );
        self.knots[i] = value;
        self.flat = match self.periodic {
            true => periodic_flatten(&self.knots, &self.multiplicities, self.degree),
            false => flatten(&self.knots, &self.multiplicities),
        };
        self.segments = Default::default();
    }

    /// Raises the degree without changing the shape.
    pub fn increase_degree(&mut self, degree: usize) {
        if degree <= self.degree {
            return;
        }
        assert!(degree <= Self::MAX_DEGREE, "B-spline degree {degree} out of range");
        let periodic = self.periodic;
        let mut spline = self.spline();
        let (lo, hi) = spline.domain();
        spline.segment(lo, hi);
        spline.elevate(degree - self.degree);
        self.set_spline(spline);
        if periodic {
            self.set_periodic();
        }
    }

    /// Inserts knot `u` with multiplicity `multiplicity` (or raises an existing knot's
    /// multiplicity to it), without changing the shape.
    #[doc(alias = "InsertKnot")]
    pub fn insert_knot(&mut self, u: f64, multiplicity: usize) {
        let u = self.snapped(u);
        let current = self.knots.iter().position(|&k| k == u).map_or(0, |i| self.multiplicities[i]);
        let multiplicity = multiplicity.min(self.degree);
        if multiplicity <= current {
            return;
        }
        if self.periodic {
            let mut periodic = self.periodic_form();
            for _ in current..multiplicity {
                periodic.insert_once(u);
            }
            self.set_periodic_form(periodic);
        } else {
            let mut spline = self.spline();
            spline.insert(u, multiplicity - current);
            self.set_spline(spline);
        }
    }

    /// Inserts several `(knot, multiplicity)` pairs.
    #[doc(alias = "InsertKnots")]
    pub fn insert_knots(&mut self, knots: &[(f64, usize)]) {
        for &(u, multiplicity) in knots {
            self.insert_knot(u, multiplicity);
        }
    }

    /// Raises the multiplicity of knot `i` to `multiplicity`.
    #[doc(alias = "IncreaseMultiplicity")]
    pub fn increase_multiplicity(&mut self, i: usize, multiplicity: usize) {
        self.insert_knot(self.knots[i], multiplicity);
    }

    /// Lowers the multiplicity of knot `i` to `multiplicity` (removing the knot at 0) if
    /// the shape changes by at most `tolerance`. Returns whether it did.
    #[doc(alias = "RemoveKnot")]
    pub fn remove_knot(&mut self, i: usize, multiplicity: usize, tolerance: f64) -> bool {
        let last = self.knots.len() - 1;
        if self.periodic || i == 0 || i >= last {
            return false;
        }
        if multiplicity >= self.multiplicities[i] {
            return true;
        }
        let u = self.knots[i];
        let mut spline = self.spline();
        for _ in multiplicity..self.multiplicities[i] {
            if !spline.remove_once(u, tolerance) {
                return false;
            }
        }
        self.set_spline(spline);
        true
    }

    /// Restricts the curve to `[u1, u2]`, keeping the parametrization.
    #[doc(alias = "Segment")]
    pub fn segment(&mut self, u1: f64, u2: f64) {
        assert!(u1 < u2, "segment [{u1}, {u2}] is empty");
        let (u1, u2) = (self.snapped(u1), self.snapped(u2));
        let mut spline = match self.periodic {
            true => self.periodic_form().unroll(u1, u2),
            false => self.spline(),
        };
        let (lo, hi) = spline.domain();
        spline.segment(u1.max(lo), u2.min(hi));
        self.set_spline(spline);
    }

    /// Makes a periodic curve non-periodic, without changing its shape.
    #[doc(alias = "SetNotPeriodic")]
    pub fn set_not_periodic(&mut self) {
        if self.periodic {
            let mut spline = self.spline();
            let (lo, hi) = spline.domain();
            spline.segment(lo, hi);
            self.set_spline(spline);
        }
    }

    /// Makes a closed curve periodic.
    ///
    /// The ends' multiplicities become the degree (the curve stays the same); panics unless
    /// the curve is clamped and its end poles (and weights) coincide.
    #[doc(alias = "SetPeriodic")]
    pub fn set_periodic(&mut self) {
        if self.periodic {
            return;
        }
        let (p, last, n) = (self.degree, self.knots.len() - 1, self.poles.len());
        assert!(
            self.multiplicities[0] == p + 1 && self.multiplicities[last] == p + 1,
            "only a clamped B-spline can be made periodic"
        );
        let rows = self.rows();
        let gap = rows[0].iter().zip(&rows[n - 1]).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(gap <= crate::tolerance::LINEAR, "only a closed B-spline can be made periodic");
        let mut multiplicities = self.multiplicities.clone();
        multiplicities[0] = p;
        multiplicities[last] = p;
        let periodic = super::knots::Periodic { degree: p, knots: self.knots.clone(), multiplicities, poles: rows[..n - 1].to_vec() };
        self.set_periodic_form(periodic);
    }

    /// For a periodic curve: makes knot `i` the start of the parameter range.
    #[doc(alias = "SetOrigin")]
    pub fn set_origin(&mut self, i: usize) {
        assert!(self.periodic, "only a periodic B-spline has an origin to move");
        let last = self.knots.len() - 1;
        if i == 0 || i == last {
            return;
        }
        let period = self.knots[last] - self.knots[0];
        let (m0, mi) = (self.multiplicities[0], self.multiplicities[i]);
        // Pole j of the new origin is old pole j + (flat knots before knot i) + mi - m0.
        let before: usize = self.multiplicities[..i].iter().sum();
        let n = self.poles.len();
        let rotation = (before + mi + n - m0) % n;
        let mut knots: Vec<f64> = self.knots[i..].to_vec();
        knots.extend(self.knots[1..=i].iter().map(|k| k + period));
        let mut multiplicities: Vec<usize> = self.multiplicities[i..].to_vec();
        multiplicities.extend_from_slice(&self.multiplicities[1..=i]);
        let rows = self.rows();
        let poles = (0..n).map(|j| rows[(j + rotation) % n].clone()).collect();
        self.set_periodic_form(super::knots::Periodic { degree: self.degree, knots, multiplicities, poles });
    }

    /// Moves the curve so that the point at `u` becomes `target`, changing only poles
    /// `first_pole..=last_pole`. Returns the range of poles actually modified, or `None`
    /// if that was not possible.
    pub fn move_point(&mut self, u: f64, target: P, first_pole: usize, last_pole: usize) -> Option<(usize, usize)> {
        // As OCCT's BSplCLib::MovePoint: the poles round the most influential one move,
        // less the further they are from it.
        if self.periodic || first_pole > last_pole || last_pole >= self.poles.len() {
            return None;
        }
        let (start, basis) = nonzero_basis(&self.flat, self.degree, u)?;
        let along = Influence::new(start, &basis, first_pole, last_pole)?;
        let (mut d1, mut d2) = (0.0, 0.0);
        for (k, &b) in basis.iter().enumerate() {
            let i = start + k;
            let h = if self.is_rational() { self.weight(i) * b } else { b };
            d2 += h;
            if along.covers(i) {
                d1 += h / (along.distance(i) + 1.0);
            }
        }
        if !(d1 > 0.0) {
            return None;
        }
        let coef = if self.is_rational() { d2 / d1 } else { 1.0 / d1 };
        let (now, goal) = (self.point(u).coords(), target.coords());
        for i in along.first..=along.last {
            let f = coef / (along.distance(i) + 1.0);
            let c = self.poles[i].coords();
            self.set_pole(i, P::from_coords([0, 1, 2].map(|a| c[a] + f * (goal[a] - now[a]))));
        }
        Some((along.first, along.last))
    }

    /// Changes the poles so that the curve passes through `target` at `u` with derivative
    /// `tangent` there, as OCCT's `MovePointAndTangent`: two smooth bumps over the poles
    /// (cubic in the poles' Greville abscissae, peaking at the abscissae nearest `u`)
    /// carry the corrections. The first `starting + 1` and last `ending + 1` poles stay
    /// put (each from -1, none, up to the degree); `u` must lie within `tolerance` inside
    /// the ends unless they are free. `None` (the curve unchanged) where there is not
    /// enough freedom, or for a periodic curve.
    #[doc(alias = "Geom_BSplineCurve::MovePointAndTangent")]
    pub fn move_point_and_tangent(&mut self, u: f64, target: P, tangent: P::Vector, tolerance: f64, starting: i32, ending: i32) -> Option<()> {
        if self.periodic {
            return None;
        }
        let degree = self.degree as i32;
        let flat = self.flat.clone();
        let n = self.poles.len();
        if !((-1..=degree).contains(&starting) && (-1..=degree).contains(&ending) && (starting + ending + 4) as usize <= n) {
            return None;
        }
        // (1-based indices, as OCCT's.)
        let knot = |i: usize| flat[i - 1];
        let (start, end) = (self.degree + 1, flat.len() - self.degree);
        let inside_start = if starting == -1 { knot(start) <= u } else { knot(start) + tolerance < u };
        let inside_end = if ending == -1 { knot(end) >= u } else { knot(end) - tolerance > u };
        if !(inside_start && inside_end) {
            return None;
        }
        // The Greville (Schoenberg) abscissae of the poles.
        let s: Vec<f64> = (1..=n).map(|i| (1..=self.degree).map(|j| knot(i + j)).sum::<f64>() / self.degree as f64).collect();
        let at = |i: usize| s[i - 1];
        let (first, last) = ((starting + 2) as usize, n - ending as usize - 1);
        let index = locate_in(&s, u, first, last);
        let other = if index == first {
            index + 1
        } else if index == last {
            index - 1
        } else if u - knot(index) < knot(index + 1) - u {
            index - 1
        } else {
            index + 1
        };
        let start_value = if first == 1 { at(1) - (at(n) - at(1)) } else { at(first - 1) };
        let end_value = if last == n { at(n) + (at(n) - at(1)) } else { at(last + 1) };
        // A bump over the poles first..=last, 1 at the pole `peak`.
        let bump = |peak: usize| {
            let mut f = vec![0.0; n + 1];
            let rise = 1.0 / (at(peak) - start_value);
            for (i, fi) in f.iter_mut().enumerate().take(peak + 1).skip(first) {
                *fi = ((at(i) - start_value) * rise).powi(3);
            }
            let fall = 1.0 / (end_value - at(peak));
            for (i, fi) in f.iter_mut().enumerate().take(last + 1).skip(peak) {
                *fi = ((end_value - at(i)) * fall).powi(3);
            }
            f
        };
        let (f1, f2) = (bump(index), bump(other));
        // Each bump's value and derivative at u, as a function on the curve's basis (and
        // weights).
        let as_function = |f: &[f64]| -> Option<(f64, f64)> {
            let poles: Vec<crate::math::Point2> = f[1..].iter().map(|&x| crate::math::Point2::new(x, 0.0)).collect();
            let mut c = BSplineCurve::<crate::math::Point2>::try_new(poles, self.knots.clone(), self.multiplicities.clone(), self.degree).ok()?;
            if let Some(w) = &self.weights {
                c = c.with_weights(w.clone());
            }
            let (p, d) = c.d1(u);
            Some((p.x, d.x))
        };
        let ((v1, d1), (v2, d2)) = (as_function(&f1)?, as_function(&f2)?);
        // [v1 d1; v2 d2]⁻¹ applied to (Δ point, Δ derivative).
        let det = v1 * d2 - d1 * v2;
        if det.abs() <= 1e-300 {
            return None;
        }
        let (now, slope) = self.d1(u);
        let (now, slope, goal, tangent) = (now.coords(), P::vector_coords(slope), target.coords(), P::vector_coords(tangent));
        let delta: [f64; 3] = std::array::from_fn(|k| goal[k] - now[k]);
        let delta_d: [f64; 3] = std::array::from_fn(|k| tangent[k] - slope[k]);
        let a: [f64; 3] = std::array::from_fn(|k| (d2 * delta[k] - v2 * delta_d[k]) / det);
        let b: [f64; 3] = std::array::from_fn(|k| (-d1 * delta[k] + v1 * delta_d[k]) / det);
        for i in 0..n {
            let c = self.poles[i].coords();
            self.set_pole(i, P::from_coords(std::array::from_fn(|k| c[k] + f1[i + 1] * a[k] + f2[i + 1] * b[k])));
        }
        Some(())
    }

    /// The knots bounding `u`: `(i, i)` if `u` is knot `i` within `tolerance`, otherwise
    /// `(i, i + 1)` with `knot(i) < u < knot(i + 1)`.
    #[doc(alias = "LocateU")]
    pub fn locate(&self, u: f64, tolerance: f64) -> (usize, usize) {
        if let Some(i) = self.knots.iter().position(|k| (k - u).abs() <= tolerance) {
            return (i, i);
        }
        let i = self.knots.partition_point(|&k| k < u);
        match i {
            0 => (0, 0),
            i if i == self.knots.len() => (i - 1, i - 1),
            i => (i - 1, i),
        }
    }

    /// The point at `u` computed on the span between knots `from` and `to` only.
    #[doc(alias = "LocalValue")]
    pub fn local_point(&self, u: f64, from: usize, to: usize) -> P {
        self.local(from, to).point(u)
    }

    /// The curve's piece between knots `from` and `to`, whose end polynomials extend it.
    fn local(&self, from: usize, to: usize) -> Self {
        let mut piece = self.clone();
        piece.segment(self.knots[from.min(to)], self.knots[from.max(to)]);
        piece
    }

    #[doc(alias = "LocalD1")]
    pub fn local_d1(&self, u: f64, from: usize, to: usize) -> (P, P::Vector) {
        self.local(from, to).d1(u)
    }

    /// True if both curves have the same degree, knots, multiplicities and periodicity,
    /// and poles (and weights) within `tolerance`.
    pub fn is_equal(&self, other: &Self, tolerance: f64) -> bool {
        let near = |a: [f64; 3], b: [f64; 3]| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt() <= tolerance;
        self.degree == other.degree
            && self.periodic == other.periodic
            && self.multiplicities == other.multiplicities
            && self.knots.len() == other.knots.len()
            && self.knots.iter().zip(&other.knots).all(|(a, b)| (a - b).abs() <= crate::tolerance::PARAMETRIC)
            && self.poles.len() == other.poles.len()
            && self.poles.iter().zip(&other.poles).all(|(a, b)| near(a.coords(), b.coords()))
            && (0..self.poles.len()).all(|i| (self.weight(i) - other.weight(i)).abs() <= crate::tolerance::PARAMETRIC)
    }
}

/// The non-zero basis functions of degree `p` on span `span` and their derivatives up to
/// `order`: `result[k][j]` is the `k`-th derivative of `N(span - p + j, p)` at `t` (The NURBS
/// Book, algorithm A2.3).

/// The span of `t` among `count` poles (as in a flat knot vector: the last flat knot <= t,
/// kept within [p, count - 1]) and the flat knots from `span - p + 1` to `span + p`, found
/// from the distinct knots and multiplicities without spelling the flat vector out, or by
/// bisection in the `flat` knots when given (for long knot vectors).
fn knot_window<const N: usize, const W: usize>(kv: &KnotVector, count: usize, t: f64, flat: Option<&[f64]>) -> ([f64; W], usize) {
    let p = kv.degree;
    if let Some(flat) = flat {
        let at_most = flat[..count.min(flat.len())].partition_point(|&k| k <= t);
        let span = at_most.saturating_sub(1).clamp(p, count - 1);
        let mut window = [0.0f64; W];
        let (lo, hi) = (span + 1 - p, (span + p).min(flat.len() - 1));
        window[..=hi - lo].copy_from_slice(&flat[lo..=hi]);
        return (window, span);
    }
    // How many of the first `count` flat knots are <= t.
    let mut at_most = 0;
    let mut seen = 0;
    for (&k, &m) in kv.knots.iter().zip(&kv.multiplicities) {
        if seen >= count || k > t {
            break;
        }
        at_most += m.min(count - seen);
        seen += m;
    }
    let span = at_most.saturating_sub(1).clamp(p, count - 1);
    let (lo, hi) = (span + 1 - p, span + p);
    let mut window = [0.0f64; W];
    let mut index = 0;
    for (&k, &m) in kv.knots.iter().zip(&kv.multiplicities) {
        let (from, to) = (index.max(lo), (index + m).min(hi + 1));
        for i in from..to {
            window[i - lo] = k;
        }
        index += m;
        if index > hi {
            break;
        }
    }
    (window, span)
}

/// The span of `t` among `count` poles, the non-zero basis functions of degree `p` there and
/// (with `derivatives`) their first derivatives (The NURBS Book, A2.2 and eq. 2.7), all on
/// the stack (arrays of `N`, `W = 2 N`, more than the degree); `None` for a higher degree.
fn fast_basis<const N: usize, const W: usize>(kv: &KnotVector, count: usize, t: f64, derivatives: bool, flat: Option<&[f64]>) -> Option<([f64; N], [f64; N], usize)> {
    if kv.degree >= N {
        return None;
    }
    let (window, span) = knot_window::<N, W>(kv, count, t, flat);
    let (n, d) = basis_and_first::<N>(&window, span + 1 - kv.degree, span, t, kv.degree, derivatives);
    Some((n, d, span))
}

/// The basis functions of degree `p` non-zero at `t` in `span` and (with `derivatives`) their
/// first derivatives.
pub(crate) fn basis_and_first<const N: usize>(flat: &[f64], offset: usize, span: usize, t: f64, p: usize, derivatives: bool) -> ([f64; N], [f64; N]) {
    let flat = |i: usize| flat[i - offset];
    let (mut n, mut previous) = ([0.0f64; N], [0.0f64; N]);
    let (mut left, mut right) = ([0.0f64; N], [0.0f64; N]);
    n[0] = 1.0;
    for j in 1..=p {
        if j == p {
            previous = n;
        }
        left[j] = t - flat(span + 1 - j);
        right[j] = flat(span + j) - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denominator = right[r + 1] + left[j - r];
            let temp = if denominator != 0.0 { n[r] / denominator } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    let mut d = [0.0f64; N];
    if derivatives {
        // N'_{i,p} = p (N_{i,p-1} / (u_{i+p} - u_i) - N_{i+1,p-1} / (u_{i+p+1} - u_{i+1})).
        for i in 0..=p {
            let mut value = 0.0;
            if i >= 1 {
                let span_width = flat(span + i) - flat(span + i - p);
                if span_width != 0.0 {
                    value += previous[i - 1] / span_width;
                }
            }
            if i < p {
                let span_width = flat(span + i + 1) - flat(span + i + 1 - p);
                if span_width != 0.0 {
                    value -= previous[i] / span_width;
                }
            }
            d[i] = p as f64 * value;
        }
    }
    (n, d)
}

/// The span of `t` and the non-zero basis functions there with their first and second
/// derivatives (The NURBS Book, A2.3), on the stack (as [`fast_basis`]); `None` for a
/// higher degree.
fn fast_basis_derivatives<const N: usize, const W: usize>(kv: &KnotVector, count: usize, t: f64, flat: Option<&[f64]>) -> Option<([[f64; N]; 3], usize)> {
    if kv.degree >= N {
        return None;
    }
    let (window, span) = knot_window::<N, W>(kv, count, t, flat);
    Some((basis_to_second::<N>(&window, span + 1 - kv.degree, span, t, kv.degree), span))
}

/// The basis functions of degree `p` non-zero at `t` in `span`, with their first and second
/// derivatives.
fn basis_to_second<const N: usize>(flat: &[f64], offset: usize, span: usize, t: f64, p: usize) -> [[f64; N]; 3] {
    let flat = |i: usize| flat[i - offset];
    let order = 2.min(p);
    let mut ndu = [[0.0f64; N]; N];
    let (mut left, mut right) = ([0.0f64; N], [0.0f64; N]);
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - flat(span + 1 - j);
        right[j] = flat(span + j) - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r] != 0.0 { ndu[r][j - 1] / ndu[j][r] } else { 0.0 };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut ders = [[0.0f64; N]; 3];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = [[0.0f64; N]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=order {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                a[s2][0] = if ndu[pk + 1][rk as usize] != 0.0 { a[s1][0] / ndu[pk + 1][rk as usize] } else { 0.0 };
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize { k - 1 } else { p - r };
            for j in j1..=j2 {
                let idx = (rk + j as isize) as usize;
                a[s2][j] = if ndu[pk + 1][idx] != 0.0 { (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx] } else { 0.0 };
                d += a[s2][j] * ndu[idx][pk];
            }
            if r <= pk {
                a[s2][k] = if ndu[pk + 1][r] != 0.0 { -a[s1][k - 1] / ndu[pk + 1][r] } else { 0.0 };
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for k in 1..=order {
        for v in &mut ders[k][..=p] {
            *v *= factor;
        }
        factor *= (p - k) as f64;
    }
    ders
}

pub(crate) fn basis_derivatives(flat: &[f64], span: usize, t: f64, p: usize, order: usize) -> Vec<Vec<f64>> {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - flat[span + 1 - j];
        right[j] = flat[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r] != 0.0 { ndu[r][j - 1] / ndu[j][r] } else { 0.0 };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut ders = vec![vec![0.0; p + 1]; order + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = vec![vec![0.0; p + 1]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=order {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                a[s2][0] = if ndu[pk + 1][rk as usize] != 0.0 { a[s1][0] / ndu[pk + 1][rk as usize] } else { 0.0 };
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize { k - 1 } else { p - r };
            for j in j1..=j2 {
                let idx = (rk + j as isize) as usize;
                a[s2][j] = if ndu[pk + 1][idx] != 0.0 { (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx] } else { 0.0 };
                d += a[s2][j] * ndu[idx][pk];
            }
            if r <= pk {
                a[s2][k] = if ndu[pk + 1][r] != 0.0 { -a[s1][k - 1] / ndu[pk + 1][r] } else { 0.0 };
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for k in 1..=order {
        for v in &mut ders[k] {
            *v *= factor;
        }
        factor *= (p - k) as f64;
    }
    ders
}

impl<P: ControlPoint> ParametricCurve for BSplineCurve<P> {
    type Point = P;
    type Vector = P::Vector;
    fn first_parameter(&self) -> f64 {
        if self.periodic { self.knots[0] } else { self.flat[self.degree] }
    }
    fn last_parameter(&self) -> f64 {
        if self.periodic {
            *self.knots.last().expect("knots")
        } else {
            self.flat[self.poles.len()]
        }
    }
    fn is_closed(&self) -> bool {
        let (a, b) = (self.point(self.first_parameter()).coords(), self.point(self.last_parameter()).coords());
        self.periodic || (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt() <= crate::tolerance::LINEAR
    }
    fn is_periodic(&self) -> bool {
        self.periodic
    }
    fn period(&self) -> Option<f64> {
        self.periodic.then(|| self.last_parameter() - self.first_parameter())
    }
    fn continuity(&self) -> Continuity {
        let n = self.multiplicities.len();
        let inner = self.multiplicities[1..n - 1].iter().copied().max().unwrap_or(0);
        match self.degree.saturating_sub(inner) {
            _ if inner == 0 => Continuity::CN,
            0 => Continuity::C0,
            1 => Continuity::C1,
            2 => Continuity::C2,
            _ => Continuity::C3,
        }
    }
    fn is_cn(&self, n: u32) -> bool {
        let m = self.multiplicities.len();
        let inner = self.multiplicities[1..m - 1].iter().copied().max().unwrap_or(0);
        inner == 0 || (self.degree - inner) as u32 >= n
    }
    fn point(&self, t: f64) -> P {
        if let Some(d) = self.fast_derivatives_by_degree(t, 0) {
            return P::from_coords(d[0]);
        }
        P::from_coords(self.derivatives(t, 0)[0])
    }
    fn d1(&self, t: f64) -> (P, P::Vector) {
        if let Some(d) = self.fast_derivatives_by_degree(t, 1) {
            return (P::from_coords(d[0]), P::vector_from_coords(d[1]));
        }
        let d = self.derivatives(t, 1);
        (P::from_coords(d[0]), P::vector_from_coords(d[1]))
    }
    fn d2(&self, t: f64) -> (P, P::Vector, P::Vector) {
        if let Some(d) = self.fast_derivatives_by_degree(t, 2) {
            return (P::from_coords(d[0]), P::vector_from_coords(d[1]), P::vector_from_coords(d[2]));
        }
        let d = self.derivatives(t, 2);
        (P::from_coords(d[0]), P::vector_from_coords(d[1]), P::vector_from_coords(d[2]))
    }
    fn d3(&self, t: f64) -> (P, P::Vector, P::Vector, P::Vector) {
        let d = self.derivatives(t, 3);
        (P::from_coords(d[0]), P::vector_from_coords(d[1]), P::vector_from_coords(d[2]), P::vector_from_coords(d[3]))
    }
    fn dn(&self, t: f64, n: u32) -> P::Vector {
        assert!(n > 0, "derivative of order 0");
        P::vector_from_coords(self.derivatives(t, n as usize)[n as usize])
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        self.first_parameter() + self.last_parameter() - t
    }
    fn reverse(&mut self) {
        let (a, b) = (self.knots[0], *self.knots.last().expect("knots"));
        if self.periodic {
            let (n, p, m0) = (self.poles.len(), self.degree, self.multiplicities[0]);
            let old = reversed_pole(n, p, m0);
            self.knots = self.knots.iter().rev().map(|&k| a + b - k).collect();
            self.multiplicities.reverse();
            self.flat = periodic_flatten(&self.knots, &self.multiplicities, p);
            self.poles = (0..n).map(|j| self.poles[old(j)]).collect();
            if let Some(w) = &self.weights {
                self.weights = Some((0..n).map(|j| w[old(j)]).collect());
            }
            self.segments = Default::default();
            return;
        }
        self.knots = self.knots.iter().rev().map(|&k| a + b - k).collect();
        self.multiplicities.reverse();
        self.flat = flatten(&self.knots, &self.multiplicities);
        self.poles.reverse();
        self.segments = Default::default();
        if let Some(w) = &mut self.weights {
            w.reverse();
        }
    }
    fn breaks(&self) -> Vec<f64> {
        let last = self.knots.len() - 1;
        self.knots[1..last].to_vec()
    }
    fn resolution(&self, tolerance: f64) -> f64 {
        // The parameter step moving the curve by at most `tolerance`: tolerance over a bound
        // of the first derivative (degree times the largest pole step over its knot span).
        let flat = &self.flat;
        let p = self.degree;
        let mut bound: f64 = 0.0;
        for i in 0..self.poles.len().saturating_sub(1) {
            let (a, b) = (self.poles[i].coords(), self.poles[i + 1].coords());
            let step = (0..3).map(|k| (b[k] - a[k]).powi(2)).sum::<f64>().sqrt();
            let span = flat.get(i + p + 1).copied().unwrap_or(0.0) - flat.get(i + 1).copied().unwrap_or(0.0);
            if span > 0.0 {
                bound = bound.max(p as f64 * step / span);
            }
        }
        if bound > 0.0 { tolerance / bound } else { f64::INFINITY }
    }
}

// (The transformed curves and surfaces keep the Bézier pieces made so far, transformed:
// cheaper than making them again.)

impl Transformable for BSplineCurve<Point> {
    fn transformed(&self, t: &Transform) -> Self {
        if *t == Transform::IDENTITY {
            return self.clone();
        }
        let map = |c: [f64; 3]| Point::new(c[0], c[1], c[2]).transformed(t).coords();
        let segments = self.segments.mapped(|s| s.map_points(|q| super::patches::map_homogeneous(q, map)));
        Self { poles: self.poles.iter().map(|p| p.transformed(t)).collect(), segments, ..self.clone() }
    }
}

impl Transformable2 for BSplineCurve<Point2> {
    fn transformed(&self, t: &Transform2) -> Self {
        if *t == Transform2::IDENTITY {
            return self.clone();
        }
        let map = |c: [f64; 3]| Point2::new(c[0], c[1]).transformed(t).coords();
        let segments = self.segments.mapped(|s| s.map_points(|q| super::patches::map_homogeneous(q, map)));
        Self { poles: self.poles.iter().map(|p| p.transformed(t)).collect(), segments, ..self.clone() }
    }
}

impl AffineTransformable for BSplineCurve<Point> {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        // Affine maps commute with (rational) B-spline evaluation: map the poles.
        let map = |c: [f64; 3]| Point::new(c[0], c[1], c[2]).affine_transformed(t).coords();
        let segments = self.segments.mapped(|s| s.map_points(|q| super::patches::map_homogeneous(q, map)));
        Self { poles: self.poles.iter().map(|p| p.affine_transformed(t)).collect(), segments, ..self.clone() }
    }
}

/// The knot vector of one parameter direction of a [`BSplineSurface`].
#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub multiplicities: Vec<usize>,
    pub periodic: bool,
}

impl KnotVector {
    /// A non-periodic knot vector.
    pub fn new(knots: Vec<f64>, multiplicities: Vec<usize>, degree: usize) -> Self {
        Self { degree, knots, multiplicities, periodic: false }
    }

    /// A periodic knot vector.
    pub fn periodic(knots: Vec<f64>, multiplicities: Vec<usize>, degree: usize) -> Self {
        Self { degree, knots, multiplicities, periodic: true }
    }
}

/// A B-spline surface: poles `poles[i][j]` (`i` along u, `j` along v), optional weights, and
/// a knot vector per direction.
#[doc(alias = "Geom_BSplineSurface")]
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineSurface {
    poles: Vec<Vec<Point>>,
    weights: Option<Vec<Vec<f64>>>,
    u: KnotVector,
    v: KnotVector,
    /// The surface as Bézier patches, for evaluating it (reset when it changes).
    patches: super::patches::PatchCache,
    /// A periodic surface as a non-periodic one over one period, for evaluating it (reset
    /// when it changes).
    unrolled: super::patches::Shared<BSplineSurface>,
    /// The flat knot vectors along u and v, for finding spans by bisection on a surface
    /// with many knots (reset when it changes).
    flat_knots: super::patches::Shared<[Vec<f64>; 2]>,
}

impl BSplineSurface {
    pub const MAX_DEGREE: usize = 25;

    /// True if the boundary curves `u = u_min` and `u = u_max` coincide within `tolerance`
    /// for `v` in `[v1, v2]`, comparing points (so rational surfaces with different but
    /// projectively equal boundary poles count as closed).
    #[doc(alias = "GeomLib::IsBSplUClosed")]
    pub fn is_u_closed_between(&self, v1: f64, v2: f64, tolerance: f64) -> bool {
        use crate::geometry::ParametricSurface;
        let (u0, u1, _, _) = self.bounds();
        (0..=20).all(|i| {
            let v = v1 + (v2 - v1) * i as f64 / 20.0;
            self.point(u0, v).distance(self.point(u1, v)) <= tolerance
        })
    }

    /// [`is_u_closed_between`](Self::is_u_closed_between) in the other direction.
    #[doc(alias = "GeomLib::IsBSplVClosed")]
    pub fn is_v_closed_between(&self, u1: f64, u2: f64, tolerance: f64) -> bool {
        use crate::geometry::ParametricSurface;
        let (_, _, v0, v1) = self.bounds();
        (0..=20).all(|i| {
            let u = u1 + (u2 - u1) * i as f64 / 20.0;
            self.point(u, v0).distance(self.point(u, v1)) <= tolerance
        })
    }

    /// Panics on inconsistent data (see [`BSplineSurface::try_new`]).
    pub fn new(poles: Vec<Vec<Point>>, u: KnotVector, v: KnotVector) -> Self {
        Self::try_new(poles, u, v).expect("invalid B-spline surface")
    }

    /// Fails unless each knot vector is valid (as for [`BSplineCurve::try_new`]) and the
    /// pole grid has the matching size.
    pub fn try_new(poles: Vec<Vec<Point>>, u: KnotVector, v: KnotVector) -> Result<Self> {
        let invalid = |msg: String| Err(Error::InvalidInput(msg));
        for kv in [&u, &v] {
            let (p, knots, mults) = (kv.degree, &kv.knots, &kv.multiplicities);
            if kv.periodic {
                let last = mults.len().saturating_sub(1);
                if !(1..=Self::MAX_DEGREE).contains(&p) || knots.len() < 2 || knots.len() != mults.len() || knots.windows(2).any(|w| !(w[0] < w[1])) {
                    return invalid(format!("invalid periodic B-spline surface knot vector {kv:?}"));
                }
                if mults.iter().any(|&m| m == 0 || m > p) || mults[0] != mults[last] {
                    return invalid(format!("invalid periodic B-spline surface multiplicities {mults:?}"));
                }
                continue;
            }
            if !(1..=Self::MAX_DEGREE).contains(&p) || knots.len() < 2 || knots.len() != mults.len() {
                return invalid(format!("invalid B-spline surface knot vector {kv:?}"));
            }
            if knots.windows(2).any(|w| !(w[0] < w[1])) {
                return invalid("B-spline surface knots must be strictly increasing".into());
            }
            let last = mults.len() - 1;
            if mults.iter().enumerate().any(|(i, &m)| m == 0 || m > if i == 0 || i == last { p + 1 } else { p }) {
                return invalid(format!("invalid B-spline surface multiplicities {mults:?}"));
            }
        }
        let count = |kv: &KnotVector| match kv.periodic {
            true => kv.multiplicities[..kv.multiplicities.len() - 1].iter().sum::<usize>(),
            false => kv.multiplicities.iter().sum::<usize>().saturating_sub(kv.degree + 1),
        };
        if poles.len() != count(&u) || poles.iter().any(|row| row.len() != count(&v)) {
            return invalid(format!("a B-spline surface with these knots needs {} x {} poles", count(&u), count(&v)));
        }
        Ok(Self { poles, weights: None, u, v, patches: Default::default(), unrolled: Default::default(), flat_knots: Default::default() })
    }

    /// What identifies this surface and its unchanged clones (for caches keyed by surface).
    pub(crate) fn id(&self) -> super::SurfaceId {
        self.patches.id()
    }

    /// The surface as Bézier patches (made on first use); `None` for one with too many.
    fn patches(&self) -> Option<&super::patches::Patches> {
        self.patches.get(|| {
            let net: Vec<Vec<[f64; 4]>> = self
                .poles
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(j, p)| {
                            let w = self.weights.as_ref().map_or(1.0, |w| w[i][j]);
                            [w * p.x, w * p.y, w * p.z, w]
                        })
                        .collect()
                })
                .collect();
            super::patches::Patches::new(&self.u, &self.v, &net)
        })
    }

    /// The points and first derivatives at `(u, v)` for each `u` of `us` (best ascending),
    /// at once: bit for bit [`ParametricSurface::d1`] at each, for less work.
    pub fn d1_along_u(&self, us: &[f64], v: f64) -> Vec<(Point, crate::math::Vector, crate::math::Vector)> {
        use super::ParametricSurface;
        if let Some(s) = self.unrolled() {
            let wrapped: Vec<f64> = us.iter().map(|&u| self.wrapped(u, v).0).collect();
            return s.d1_along_u(&wrapped, self.wrapped(us.first().copied().unwrap_or(0.0), v).1);
        }
        let mut found = Vec::new();
        if !self.patches().is_some_and(|p| p.evaluate_along_u(us, v, 1, &mut found)) {
            return us.iter().map(|&u| self.d1(u, v)).collect();
        }
        us.iter().zip(found).map(|(&u, d)| d.map_or_else(|| self.d1(u, v), |d| (d.point, d.du, d.dv))).collect()
    }

    /// The points at `(u, v)` for each `u` of `us` (best ascending), bit for bit
    /// [`ParametricSurface::point`] at each.
    pub fn points_along_u(&self, us: &[f64], v: f64) -> Vec<Point> {
        use super::ParametricSurface;
        if let Some(s) = self.unrolled() {
            let wrapped: Vec<f64> = us.iter().map(|&u| self.wrapped(u, v).0).collect();
            return s.points_along_u(&wrapped, self.wrapped(us.first().copied().unwrap_or(0.0), v).1);
        }
        let mut found = Vec::new();
        if !self.patches().is_some_and(|p| p.evaluate_along_u(us, v, 0, &mut found)) {
            return us.iter().map(|&u| self.point(u, v)).collect();
        }
        us.iter().zip(found).map(|(&u, d)| d.map_or_else(|| self.point(u, v), |d| d.point)).collect()
    }

    /// The flat knot vectors, kept for surfaces with many knots (where scanning the knots
    /// for a span would cost more than evaluating).
    fn flat_knots(&self) -> Option<&[Vec<f64>; 2]> {
        let flat = |kv: &KnotVector| -> Vec<f64> { kv.knots.iter().zip(&kv.multiplicities).flat_map(|(&k, &m)| std::iter::repeat_n(k, m)).collect() };
        self.flat_knots.get(|| (self.u.knots.len() + self.v.knots.len() > 64).then(|| [flat(&self.u), flat(&self.v)]))
    }

    /// The point and derivatives up to `order` from the Bézier patches.
    fn patch_derivatives(&self, u: f64, v: f64, order: usize) -> Option<super::patches::Derivatives> {
        self.patches()?.evaluate(u, v, order)
    }

    /// [`Self::fast_d1`] with arrays sized for the degrees (they are zeroed on every
    /// evaluation).
    fn fast_d1_by_degree(&self, u: f64, v: f64, derivatives: bool) -> Option<(Point, crate::math::Vector, crate::math::Vector)> {
        match self.u.degree.max(self.v.degree) {
            0..=3 => self.fast_d1::<4, 8>(u, v, derivatives),
            4..=8 => self.fast_d1::<9, 18>(u, v, derivatives),
            9..=25 => self.fast_d1::<26, 52>(u, v, derivatives),
            _ => None,
        }
    }

    /// [`Self::fast_d2`] with arrays sized for the degrees.
    fn fast_d2_by_degree(&self, u: f64, v: f64) -> Option<crate::geometry::SurfaceD2> {
        match self.u.degree.max(self.v.degree) {
            0..=3 => self.fast_d2::<4, 8>(u, v),
            4..=8 => self.fast_d2::<9, 18>(u, v),
            9..=25 => self.fast_d2::<26, 52>(u, v),
            _ => None,
        }
    }

    /// The point and first derivatives at `(u, v)` without allocating (the frequent case):
    /// basis functions in fixed arrays from the flat knots gathered on the stack. `None` for
    /// knot vectors too long for the buffer.
    fn fast_d1<const N: usize, const W: usize>(&self, u: f64, v: f64, derivatives: bool) -> Option<(Point, crate::math::Vector, crate::math::Vector)> {
        let (nu, nv) = (self.poles.len(), self.poles[0].len());
        let flat = self.flat_knots();
        let (bu, du, su) = fast_basis::<N, W>(&self.u, nu, u, derivatives, flat.map(|f| &f[0][..]))?;
        let (bv, dv, sv) = fast_basis::<N, W>(&self.v, nv, v, derivatives, flat.map(|f| &f[1][..]))?;
        let (pu, pv) = (self.u.degree, self.v.degree);
        // Homogeneous sums: (w P, w) with the basis, and its u and v derivatives.
        let (mut a, mut au, mut av) = ([0.0f64; 4], [0.0f64; 4], [0.0f64; 4]);
        for i in 0..=pu {
            let row = &self.poles[su - pu + i];
            let weights = self.weights.as_ref().map(|w| &w[su - pu + i]);
            for j in 0..=pv {
                let jj = sv - pv + j;
                let w = weights.map_or(1.0, |w| w[jj]);
                let pt = row[jj];
                let h = [w * pt.x, w * pt.y, w * pt.z, w];
                let (f, fu, fv) = (bu[i] * bv[j], du[i] * bv[j], bu[i] * dv[j]);
                for k in 0..4 {
                    a[k] += f * h[k];
                    au[k] += fu * h[k];
                    av[k] += fv * h[k];
                }
            }
        }
        let w = a[3];
        let p = [a[0] / w, a[1] / w, a[2] / w];
        let point = Point::new(p[0], p[1], p[2]);
        if !derivatives {
            return Some((point, crate::math::Vector::new(0.0, 0.0, 0.0), crate::math::Vector::new(0.0, 0.0, 0.0)));
        }
        let along = |d: [f64; 4]| crate::math::Vector::new((d[0] - d[3] * p[0]) / w, (d[1] - d[3] * p[1]) / w, (d[2] - d[3] * p[2]) / w);
        Some((point, along(au), along(av)))
    }

    /// [`Self::fast_d1`] with the second derivatives.
    fn fast_d2<const N: usize, const W: usize>(&self, u: f64, v: f64) -> Option<crate::geometry::SurfaceD2> {
        let (nu, nv) = (self.poles.len(), self.poles[0].len());
        let flat = self.flat_knots();
        let (bu, su) = fast_basis_derivatives::<N, W>(&self.u, nu, u, flat.map(|f| &f[0][..]))?;
        let (bv, sv) = fast_basis_derivatives::<N, W>(&self.v, nv, v, flat.map(|f| &f[1][..]))?;
        let (pu, pv) = (self.u.degree, self.v.degree);
        // Homogeneous sums a[k][l] = Σ N_i^(k) N_j^(l) (w P, w) for k + l <= 2.
        let mut a = [[[0.0f64; 4]; 3]; 3];
        for i in 0..=pu {
            let row = &self.poles[su - pu + i];
            let weights = self.weights.as_ref().map(|w| &w[su - pu + i]);
            for j in 0..=pv {
                let jj = sv - pv + j;
                let w = weights.map_or(1.0, |w| w[jj]);
                let pt = row[jj];
                let h = [w * pt.x, w * pt.y, w * pt.z, w];
                for (k, l) in [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)] {
                    let f = bu[k][i] * bv[l][j];
                    for c in 0..4 {
                        a[k][l][c] += f * h[c];
                    }
                }
            }
        }
        let w = a[0][0][3];
        let vector = |c: [f64; 3]| crate::math::Vector::new(c[0], c[1], c[2]);
        let s = [a[0][0][0] / w, a[0][0][1] / w, a[0][0][2] / w];
        let first = |k: usize, l: usize| std::array::from_fn::<f64, 3, _>(|c| (a[k][l][c] - a[k][l][3] * s[c]) / w);
        let (s_u, s_v) = (first(1, 0), first(0, 1));
        let s_uu: [f64; 3] = std::array::from_fn(|c| (a[2][0][c] - 2.0 * a[1][0][3] * s_u[c] - a[2][0][3] * s[c]) / w);
        let s_vv: [f64; 3] = std::array::from_fn(|c| (a[0][2][c] - 2.0 * a[0][1][3] * s_v[c] - a[0][2][3] * s[c]) / w);
        let s_uv: [f64; 3] = std::array::from_fn(|c| (a[1][1][c] - a[1][0][3] * s_v[c] - a[0][1][3] * s_u[c] - a[1][1][3] * s[c]) / w);
        Some(crate::geometry::SurfaceD2 {
            point: Point::new(s[0], s[1], s[2]),
            du: vector(s_u),
            dv: vector(s_v),
            duu: vector(s_uu),
            dvv: vector(s_vv),
            duv: vector(s_uv),
        })
    }

    /// Derivatives `S(k, l)` for `k <= ku`, `l <= kv` at `(u, v)` (outside the bounds,
    /// those of the end spans' polynomials), as coordinates (The NURBS Book, A3.6 and A4.4).
    fn derivatives(&self, u: f64, v: f64, ku: usize, kv: usize) -> Vec<Vec<[f64; 3]>> {
        let along = |kvec: &KnotVector, count: usize, t: f64, order: usize| {
            let flat = flatten(&kvec.knots, &kvec.multiplicities);
            let p = kvec.degree;
            let span = flat[..count].partition_point(|&k| k <= t).saturating_sub(1).clamp(p, count - 1);
            (span, basis_derivatives(&flat, span, t, p, order.min(p)))
        };
        let (nu, nv) = (self.poles.len(), self.poles[0].len());
        let (su, bu) = along(&self.u, nu, u, ku);
        let (sv, bv) = along(&self.v, nv, v, kv);
        let (pu, pv) = (self.u.degree, self.v.degree);
        // Homogeneous derivatives A(k, l) = Σ N_i^(k)(u) N_j^(l)(v) (w P, w).
        let mut a = vec![vec![[0.0; 4]; kv + 1]; ku + 1];
        for (k, ru) in bu.iter().enumerate() {
            for (l, rv) in bv.iter().enumerate() {
                for (i, &nui) in ru.iter().enumerate() {
                    for (j, &nvj) in rv.iter().enumerate() {
                        let (ii, jj) = (su - pu + i, sv - pv + j);
                        let w = self.weight(ii, jj);
                        let pt = self.poles[ii][jj];
                        let f = nui * nvj * w;
                        a[k][l][0] += f * pt.x;
                        a[k][l][1] += f * pt.y;
                        a[k][l][2] += f * pt.z;
                        a[k][l][3] += f;
                    }
                }
            }
        }
        let binom = |n: usize, r: usize| (1..=r).fold(1.0, |acc, i| acc * (n - r + i) as f64 / i as f64);
        let mut s = vec![vec![[0.0; 3]; kv + 1]; ku + 1];
        for k in 0..=ku {
            for l in 0..=kv {
                let mut c = [a[k][l][0], a[k][l][1], a[k][l][2]];
                for j in 1..=l {
                    let f = binom(l, j) * a[0][j][3];
                    for d in 0..3 {
                        c[d] -= f * s[k][l - j][d];
                    }
                }
                for i in 1..=k {
                    let f = binom(k, i) * a[i][0][3];
                    for d in 0..3 {
                        c[d] -= f * s[k - i][l][d];
                    }
                    for j in 1..=l {
                        let g = binom(k, i) * binom(l, j) * a[i][j][3];
                        for d in 0..3 {
                            c[d] -= g * s[k - i][l - j][d];
                        }
                    }
                }
                for d in 0..3 {
                    c[d] /= a[0][0][3];
                }
                s[k][l] = c;
            }
        }
        s
    }

    /// Makes the surface rational. Panics unless the weights match the poles and are
    /// positive.
    pub fn with_weights(mut self, weights: Vec<Vec<f64>>) -> Self {
        assert!(
            weights.len() == self.poles.len()
                && weights.iter().all(|row| row.len() == self.poles[0].len())
                && weights.iter().flatten().all(|&w| w > 0.0),
            "invalid B-spline surface weights"
        );
        self.weights = Some(weights);
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
        self
    }

    pub fn u_degree(&self) -> usize {
        self.u.degree
    }

    pub fn v_degree(&self) -> usize {
        self.v.degree
    }

    pub fn u_knot_vector(&self) -> &KnotVector {
        &self.u
    }

    pub fn v_knot_vector(&self) -> &KnotVector {
        &self.v
    }

    pub fn u_knots(&self) -> &[f64] {
        &self.u.knots
    }

    pub fn v_knots(&self) -> &[f64] {
        &self.v.knots
    }

    pub fn u_knot(&self, i: usize) -> f64 {
        self.u.knots[i]
    }

    pub fn v_knot(&self, i: usize) -> f64 {
        self.v.knots[i]
    }

    pub fn u_multiplicity(&self, i: usize) -> usize {
        self.u.multiplicities[i]
    }

    pub fn v_multiplicity(&self, i: usize) -> usize {
        self.v.multiplicities[i]
    }

    pub fn u_multiplicities(&self) -> &[usize] {
        &self.u.multiplicities
    }

    pub fn v_multiplicities(&self) -> &[usize] {
        &self.v.multiplicities
    }

    pub fn u_pole_count(&self) -> usize {
        self.poles.len()
    }

    pub fn v_pole_count(&self) -> usize {
        self.poles[0].len()
    }

    pub fn poles(&self) -> &[Vec<Point>] {
        &self.poles
    }

    pub fn pole(&self, i: usize, j: usize) -> Point {
        self.poles[i][j]
    }

    pub fn weights(&self) -> Option<&[Vec<f64>]> {
        self.weights.as_deref()
    }

    pub fn weight(&self, i: usize, j: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |w| w[i][j])
    }

    /// True if the weights vary along u.
    pub fn is_u_rational(&self) -> bool {
        // Some row of weights differs from another.
        self.weights.as_ref().is_some_and(|w| w.iter().any(|row| row.iter().zip(&w[0]).any(|(a, b)| (a - b).abs() > f64::EPSILON * a.abs().max(1.0))))
    }

    /// True if the weights vary along v.
    pub fn is_v_rational(&self) -> bool {
        self.weights.as_ref().is_some_and(|w| w.iter().any(|row| row.iter().any(|&x| (x - row[0]).abs() > f64::EPSILON * x.abs().max(1.0))))
    }

    pub fn set_pole(&mut self, i: usize, j: usize, p: Point) {
        self.poles[i][j] = p;
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }

    /// Panics if `weight <= 0`.
    pub fn set_weight(&mut self, i: usize, j: usize, weight: f64) {
        assert!(weight > 0.0, "non-positive weight {weight}");
        let (nu, nv) = (self.u_pole_count(), self.v_pole_count());
        self.weights.get_or_insert_with(|| vec![vec![1.0; nv]; nu])[i][j] = weight;
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }

    pub fn set_pole_row(&mut self, i: usize, poles: &[Point]) {
        for (j, &p) in poles.iter().enumerate() {
            self.set_pole(i, j, p);
        }
    }

    pub fn set_pole_col(&mut self, j: usize, poles: &[Point]) {
        for (i, &p) in poles.iter().enumerate() {
            self.set_pole(i, j, p);
        }
    }

    pub fn set_weight_row(&mut self, i: usize, weights: &[f64]) {
        for (j, &w) in weights.iter().enumerate() {
            self.set_weight(i, j, w);
        }
    }

    pub fn set_weight_col(&mut self, j: usize, weights: &[f64]) {
        for (i, &w) in weights.iter().enumerate() {
            self.set_weight(i, j, w);
        }
    }

    pub fn set_u_knot(&mut self, i: usize, value: f64) {
        let last = self.u.knots.len() - 1;
        assert!(
            (i == 0 || self.u.knots[i - 1] < value) && (i == last || value < self.u.knots[i + 1]),
            "u knot {i} moved to {value} would break the knots' order"
        );
        self.u.knots[i] = value;
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }

    pub fn set_v_knot(&mut self, i: usize, value: f64) {
        self.along_v(|s| s.set_u_knot(i, value));
    }

    /// Raises the degrees without changing the shape.
    pub fn increase_degree(&mut self, u_degree: usize, v_degree: usize) {
        if u_degree > self.u.degree {
            self.set_u_not_periodic();
            let lo_hi = (flatten(&self.u.knots, &self.u.multiplicities), self.poles.len());
            let (lo, hi) = (lo_hi.0[self.u.degree], lo_hi.0[lo_hi.1]);
            let t = u_degree - self.u.degree;
            self.edit_u(|s| {
                s.segment(lo, hi);
                s.elevate(t);
            });
        }
        if v_degree > self.v.degree {
            self.along_v(|s| s.increase_degree(v_degree, s.v.degree));
        }
    }

    pub fn insert_u_knot(&mut self, u: f64, multiplicity: usize) {
        let u = self.u.knots.iter().copied().find(|k| (k - u).abs() <= crate::tolerance::PARAMETRIC).unwrap_or(u);
        let current = self.u.knots.iter().position(|&k| k == u).map_or(0, |i| self.u.multiplicities[i]);
        let multiplicity = multiplicity.min(self.u.degree);
        if multiplicity <= current {
            return;
        }
        if self.u.periodic {
            let mut periodic = super::knots::Periodic {
                degree: self.u.degree,
                knots: self.u.knots.clone(),
                multiplicities: self.u.multiplicities.clone(),
                poles: self.u_rows(),
            };
            for _ in current..multiplicity {
                periodic.insert_once(u);
            }
            let kv = KnotVector::periodic(periodic.knots, periodic.multiplicities, periodic.degree);
            self.set_u_rows(periodic.poles, kv);
        } else {
            self.edit_u(|s| s.insert(u, multiplicity - current));
        }
    }

    pub fn insert_v_knot(&mut self, v: f64, multiplicity: usize) {
        self.along_v(|s| s.insert_u_knot(v, multiplicity));
    }

    pub fn increase_u_multiplicity(&mut self, i: usize, multiplicity: usize) {
        self.insert_u_knot(self.u.knots[i], multiplicity);
    }

    /// See [`BSplineCurve::remove_knot`].
    pub fn remove_u_knot(&mut self, i: usize, multiplicity: usize, tolerance: f64) -> bool {
        let last = self.u.knots.len() - 1;
        if self.u.periodic || i == 0 || i >= last {
            return false;
        }
        if multiplicity >= self.u.multiplicities[i] {
            return true;
        }
        let u = self.u.knots[i];
        let times = self.u.multiplicities[i] - multiplicity;
        let mut trial = self.clone();
        let mut ok = true;
        trial.edit_u(|s| {
            for _ in 0..times {
                if !s.remove_once(u, tolerance) {
                    ok = false;
                    break;
                }
            }
        });
        if ok {
            *self = trial;
        }
        ok
    }

    pub fn remove_v_knot(&mut self, i: usize, multiplicity: usize, tolerance: f64) -> bool {
        self.along_v(|s| s.remove_u_knot(i, multiplicity, tolerance))
    }

    /// Restricts the surface to `[u1, u2] x [v1, v2]`, keeping the parametrization.
    pub fn segment(&mut self, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.segment_u(u1, u2);
        self.along_v(|s| s.segment_u(v1, v2));
    }

    /// Restricts the surface to `[u1, u2]` along u (non-periodic and clamped there).
    fn segment_u(&mut self, u1: f64, u2: f64) {
        let snap = |kv: &KnotVector, x: f64| kv.knots.iter().copied().find(|k| (k - x).abs() <= crate::tolerance::PARAMETRIC).unwrap_or(x);
        let (u1, u2) = (snap(&self.u, u1), snap(&self.u, u2));
        if self.u.periodic {
            let periodic = super::knots::Periodic {
                degree: self.u.degree,
                knots: self.u.knots.clone(),
                multiplicities: self.u.multiplicities.clone(),
                poles: self.u_rows(),
            };
            let spline = periodic.unroll(u1, u2);
            let (knots, multiplicities) = super::knots::group(&spline.flat);
            let kv = KnotVector::new(knots, multiplicities, spline.degree);
            self.set_u_rows(spline.poles, kv);
        }
        let flat = flatten(&self.u.knots, &self.u.multiplicities);
        let (lo, hi) = (flat[self.u.degree], flat[self.poles.len()]);
        self.edit_u(|s| s.segment(u1.max(lo), u2.min(hi)));
    }

    pub fn exchange_uv(&mut self) {
        *self = self.swapped();
    }

    pub fn set_u_not_periodic(&mut self) {
        if !self.u.periodic {
            return;
        }
        let (lo, hi) = (self.u.knots[0], self.u.knots[self.u.knots.len() - 1]);
        let periodic = super::knots::Periodic {
            degree: self.u.degree,
            knots: self.u.knots.clone(),
            multiplicities: self.u.multiplicities.clone(),
            poles: self.u_rows(),
        };
        let mut spline = periodic.unroll(lo, hi);
        spline.segment(lo, hi);
        let (knots, multiplicities) = super::knots::group(&spline.flat);
        let kv = KnotVector::new(knots, multiplicities, spline.degree);
        self.set_u_rows(spline.poles, kv);
    }

    pub fn set_v_not_periodic(&mut self) {
        if self.v.periodic {
            self.along_v(|s| s.set_u_not_periodic());
        }
    }

    /// See [`BSplineCurve::move_point`]; pole ranges are given per direction and the
    /// modified ranges `(u_first, u_last, v_first, v_last)` are returned.
    pub fn move_point(
        &mut self,
        u: f64,
        v: f64,
        target: Point,
        u_poles: (usize, usize),
        v_poles: (usize, usize),
    ) -> Option<(usize, usize, usize, usize)> {
        // As OCCT's BSplSLib::MovePoint.
        let (rows, cols) = (self.poles.len(), self.poles[0].len());
        if self.u.periodic || self.v.periodic || u_poles.0 > u_poles.1 || v_poles.0 > v_poles.1 || u_poles.1 >= rows || v_poles.1 >= cols {
            return None;
        }
        let (u_start, u_basis) = nonzero_basis(&flatten(&self.u.knots, &self.u.multiplicities), self.u.degree, u)?;
        let (v_start, v_basis) = nonzero_basis(&flatten(&self.v.knots, &self.v.multiplicities), self.v.degree, v)?;
        let along_u = Influence::new(u_start, &u_basis, u_poles.0, u_poles.1)?;
        let along_v = Influence::new(v_start, &v_basis, v_poles.0, v_poles.1)?;
        let rational = self.weights.is_some();
        let (mut d1, mut d2) = (0.0, 0.0);
        for (a, &bu) in u_basis.iter().enumerate() {
            for (b, &bv) in v_basis.iter().enumerate() {
                let (i, j) = (u_start + a, v_start + b);
                let h = if rational { self.weight(i, j) * bu * bv } else { bu * bv };
                d2 += h;
                if along_u.covers(i) && along_v.covers(j) {
                    d1 += h / (along_u.distance(i) + along_v.distance(j) + 1.0);
                }
            }
        }
        if !(d1 > 0.0) {
            return None;
        }
        let coef = if rational { d2 / d1 } else { 1.0 / d1 };
        let displacement = target - super::ParametricSurface::point(self, u, v);
        for i in along_u.first..=along_u.last {
            for j in along_v.first..=along_v.last {
                let f = coef / (along_u.distance(i) + along_v.distance(j) + 1.0);
                let moved = self.poles[i][j] + displacement * f;
                self.set_pole(i, j, moved);
            }
        }
        Some((along_u.first, along_u.last, along_v.first, along_v.last))
    }

    /// See [`BSplineCurve::locate`].
    pub fn locate_u(&self, u: f64, tolerance: f64) -> (usize, usize) {
        locate_knot(&self.u.knots, u, tolerance)
    }

    pub fn locate_v(&self, v: f64, tolerance: f64) -> (usize, usize) {
        locate_knot(&self.v.knots, v, tolerance)
    }

    /// The point at `(u, v)` computed on the given knot spans only.
    #[doc(alias = "LocalD0")]
    pub fn local_point(&self, u: f64, v: f64, u_span: (usize, usize), v_span: (usize, usize)) -> Point {
        use crate::geometry::ParametricSurface;
        self.local(u_span, v_span).point(u, v)
    }

    /// The derivatives at `(u, v)` computed on the given knot spans only.
    #[doc(alias = "LocalD1")]
    pub fn local_d1(&self, u: f64, v: f64, u_span: (usize, usize), v_span: (usize, usize)) -> (Point, crate::math::Vector, crate::math::Vector) {
        use crate::geometry::ParametricSurface;
        self.local(u_span, v_span).d1(u, v)
    }
}

impl BSplineSurface {
    /// For a periodic surface, the same surface non-periodic over its period (made on first
    /// use); `None` for a non-periodic one.
    fn unrolled(&self) -> Option<&BSplineSurface> {
        if !self.u.periodic && !self.v.periodic {
            return None;
        }
        self.unrolled.get(|| {
            let mut s = self.clone();
            s.unrolled = Default::default();
            s.set_u_not_periodic();
            s.set_v_not_periodic();
            Some(s)
        })
    }

    /// `(u, v)` moved into the period in each periodic direction.
    fn wrapped(&self, u: f64, v: f64) -> (f64, f64) {
        let wrap = |x: f64, kv: &KnotVector| match kv.periodic {
            true => {
                let (lo, hi) = (kv.knots[0], kv.knots[kv.knots.len() - 1]);
                lo + (x - lo).rem_euclid(hi - lo)
            }
            false => x,
        };
        (wrap(u, &self.u), wrap(v, &self.v))
    }

    /// The rows of poles along u as homogeneous numbers: row `i` holds `w·x, w·y, w·z, w`
    /// of each pole of it.
    fn u_rows(&self) -> Vec<Vec<f64>> {
        (0..self.poles.len())
            .map(|i| {
                (0..self.poles[i].len())
                    .flat_map(|j| {
                        let (p, w) = (self.poles[i][j], self.weight(i, j));
                        [w * p.x, w * p.y, w * p.z, w]
                    })
                    .collect()
            })
            .collect()
    }

    /// Takes poles and weights from homogeneous rows along u, and the u knot vector.
    fn set_u_rows(&mut self, rows: Vec<Vec<f64>>, u: KnotVector) {
        let rational = self.weights.is_some();
        let poles: Vec<Vec<Point>> = rows.iter().map(|r| r.chunks(4).map(|c| Point::new(c[0] / c[3], c[1] / c[3], c[2] / c[3])).collect()).collect();
        let weights: Vec<Vec<f64>> = rows.iter().map(|r| r.chunks(4).map(|c| c[3]).collect()).collect();
        self.poles = poles;
        self.weights = rational.then_some(weights);
        self.u = u;
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }

    /// Runs a knot operation along u on the (non-periodic) surface.
    pub(crate) fn edit_u(&mut self, f: impl FnOnce(&mut super::knots::Spline)) {
        let mut spline = super::knots::Spline { degree: self.u.degree, flat: flatten(&self.u.knots, &self.u.multiplicities), poles: self.u_rows() };
        f(&mut spline);
        let (knots, multiplicities) = super::knots::group(&spline.flat);
        let u = KnotVector::new(knots, multiplicities, spline.degree);
        self.set_u_rows(spline.poles, u);
    }

    /// Runs `f` with u and v exchanged.
    pub(crate) fn along_v<R>(&mut self, f: impl FnOnce(&mut BSplineSurface) -> R) -> R {
        let mut swapped = self.swapped();
        let r = f(&mut swapped);
        *self = swapped.swapped();
        r
    }
}

/// The first pole index whose basis function of degree `p` is non-zero at `u`, and the
/// `p + 1` values from there.
fn nonzero_basis(flat: &[f64], p: usize, u: f64) -> Option<(usize, Vec<f64>)> {
    let m = flat.len();
    if m < 2 * p + 2 || !u.is_finite() {
        return None;
    }
    // The span: flat[s] <= u < flat[s + 1], within the domain's spans.
    let s = flat.partition_point(|&k| k <= u).saturating_sub(1).clamp(p, m - p - 2);
    let mut n = vec![0.0; p + 1];
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = u - flat[s + 1 - j];
        right[j] = flat[s + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let t = if denom != 0.0 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * t;
            saved = left[j - r] * t;
        }
        n[j] = saved;
    }
    Some((s - p, n))
}

/// Which poles a point move shifts (OCCT's BSplCLib::MovePoint): those allowed and
/// non-zero at the point, by their index distance from the one or two with the largest
/// basis value.
struct Influence {
    first: usize,
    last: usize,
    peak: (usize, usize),
}

impl Influence {
    fn new(start: usize, basis: &[f64], allowed_first: usize, allowed_last: usize) -> Option<Self> {
        let first = start.max(allowed_first);
        let last = (start + basis.len() - 1).min(allowed_last);
        if first > last {
            return None;
        }
        let (mut peak, mut best) = (first, 0.0);
        for i in first..=last {
            if basis[i - start] > best {
                (peak, best) = (i, basis[i - start]);
            }
        }
        let second = if peak + 1 <= last && (basis[peak + 1 - start] - best).abs() < 1e-10 { peak + 1 } else { peak };
        Some(Self { first, last, peak: (peak, second) })
    }

    fn covers(&self, i: usize) -> bool {
        (self.first..=self.last).contains(&i)
    }

    fn distance(&self, i: usize) -> f64 {
        if i < self.peak.0 {
            (self.peak.0 - i) as f64
        } else if i > self.peak.1 {
            (i - self.peak.1) as f64
        } else {
            0.0
        }
    }
}

/// The knots bounding `x` (see [`BSplineCurve::locate`]).
fn locate_knot(knots: &[f64], x: f64, tolerance: f64) -> (usize, usize) {
    if let Some(i) = knots.iter().position(|k| (k - x).abs() <= tolerance) {
        return (i, i);
    }
    match knots.partition_point(|&k| k < x) {
        0 => (0, 0),
        i if i == knots.len() => (i - 1, i - 1),
        i => (i - 1, i),
    }
}

impl BSplineSurface {
    /// The piece between the given knots in each direction (its end polynomials extend it).
    fn local(&self, (a, b): (usize, usize), (c, d): (usize, usize)) -> BSplineSurface {
        let mut piece = self.clone();
        let (u1, u2) = (self.u.knots[a.min(b)], self.u.knots[a.max(b)]);
        let (v1, v2) = (self.v.knots[c.min(d)], self.v.knots[c.max(d)]);
        if u1 < u2 && v1 < v2 {
            piece.segment(u1, u2, v1, v2);
        }
        piece
    }
}

impl crate::geometry::ParametricSurface for BSplineSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        let (u, v) = (&self.u, &self.v);
        let range = |kv: &KnotVector, count: usize| match kv.periodic {
            true => (kv.knots[0], kv.knots[kv.knots.len() - 1]),
            false => (flat_knot(kv, kv.degree), flat_knot(kv, count)),
        };
        let ((a, b), (c, d)) = (range(u, self.poles.len()), range(v, self.poles[0].len()));
        (a, b, c, d)
    }
    fn is_u_closed(&self) -> bool {
        if self.u.periodic {
            return true;
        }
        let (a, b) = (&self.poles[0], &self.poles[self.poles.len() - 1]);
        a.iter().zip(b).all(|(p, q)| p.distance(*q) <= crate::tolerance::LINEAR)
    }
    fn is_v_closed(&self) -> bool {
        if self.v.periodic {
            return true;
        }
        let n = self.poles[0].len();
        self.poles.iter().all(|row| row[0].distance(row[n - 1]) <= crate::tolerance::LINEAR)
    }
    fn is_u_periodic(&self) -> bool {
        self.u.periodic
    }
    fn is_v_periodic(&self) -> bool {
        self.v.periodic
    }
    fn u_period(&self) -> Option<f64> {
        let (a, b, _, _) = self.bounds();
        self.u.periodic.then_some(b - a)
    }
    fn v_period(&self) -> Option<f64> {
        let (_, _, a, b) = self.bounds();
        self.v.periodic.then_some(b - a)
    }
    fn continuity(&self) -> Continuity {
        let inner = |kv: &KnotVector| {
            let n = kv.multiplicities.len();
            let m = kv.multiplicities[1..n - 1].iter().copied().max().unwrap_or(0);
            if m == 0 { usize::MAX } else { kv.degree - m }
        };
        match inner(&self.u).min(inner(&self.v)) {
            0 => Continuity::C0,
            1 => Continuity::C1,
            2 => Continuity::C2,
            usize::MAX => Continuity::CN,
            _ => Continuity::C3,
        }
    }
    fn is_cn_u(&self, n: u32) -> bool {
        let m = self.u.multiplicities.len();
        let inner = self.u.multiplicities[1..m - 1].iter().copied().max().unwrap_or(0);
        inner == 0 || (self.u.degree - inner) as u32 >= n
    }
    fn is_cn_v(&self, n: u32) -> bool {
        let m = self.v.multiplicities.len();
        let inner = self.v.multiplicities[1..m - 1].iter().copied().max().unwrap_or(0);
        inner == 0 || (self.v.degree - inner) as u32 >= n
    }
    fn point(&self, u: f64, v: f64) -> Point {
        if let Some(s) = self.unrolled() {
            let (u, v) = self.wrapped(u, v);
            return s.point(u, v);
        }
        if let Some(d) = self.patch_derivatives(u, v, 0) {
            return d.point;
        }
        if let Some((p, _, _)) = self.fast_d1_by_degree(u, v, false) {
            return p;
        }
        let c = self.derivatives(u, v, 0, 0)[0][0];
        Point::new(c[0], c[1], c[2])
    }
    fn d1(&self, u: f64, v: f64) -> (Point, crate::math::Vector, crate::math::Vector) {
        if let Some(s) = self.unrolled() {
            let (u, v) = self.wrapped(u, v);
            return s.d1(u, v);
        }
        if let Some(d) = self.patch_derivatives(u, v, 1) {
            return (d.point, d.du, d.dv);
        }
        if let Some(d) = self.fast_d1_by_degree(u, v, true) {
            return d;
        }
        let d = self.derivatives(u, v, 1, 1);
        let vec = |c: [f64; 3]| crate::math::Vector::new(c[0], c[1], c[2]);
        (Point::new(d[0][0][0], d[0][0][1], d[0][0][2]), vec(d[1][0]), vec(d[0][1]))
    }
    fn d2(&self, u: f64, v: f64) -> crate::geometry::SurfaceD2 {
        if let Some(s) = self.unrolled() {
            let (u, v) = self.wrapped(u, v);
            return s.d2(u, v);
        }
        if let Some(d) = self.patch_derivatives(u, v, 2) {
            return crate::geometry::SurfaceD2 { point: d.point, du: d.du, dv: d.dv, duu: d.duu, dvv: d.dvv, duv: d.duv };
        }
        if let Some(d) = self.fast_d2_by_degree(u, v) {
            return d;
        }
        let d = self.derivatives(u, v, 2, 2);
        let vec = |c: [f64; 3]| crate::math::Vector::new(c[0], c[1], c[2]);
        crate::geometry::SurfaceD2 {
            point: Point::new(d[0][0][0], d[0][0][1], d[0][0][2]),
            du: vec(d[1][0]),
            dv: vec(d[0][1]),
            duu: vec(d[2][0]),
            dvv: vec(d[0][2]),
            duv: vec(d[1][1]),
        }
    }
    fn d3(&self, u: f64, v: f64) -> crate::geometry::SurfaceD3 {
        if let Some(s) = self.unrolled() {
            let (u, v) = self.wrapped(u, v);
            return s.d3(u, v);
        }
        let d = self.derivatives(u, v, 3, 3);
        let vec = |c: [f64; 3]| crate::math::Vector::new(c[0], c[1], c[2]);
        crate::geometry::SurfaceD3 {
            point: Point::new(d[0][0][0], d[0][0][1], d[0][0][2]),
            du: vec(d[1][0]),
            dv: vec(d[0][1]),
            duu: vec(d[2][0]),
            dvv: vec(d[0][2]),
            duv: vec(d[1][1]),
            duuu: vec(d[3][0]),
            dvvv: vec(d[0][3]),
            duuv: vec(d[2][1]),
            duvv: vec(d[1][2]),
        }
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> crate::math::Vector {
        assert!(nu + nv > 0, "derivative of order 0");
        if let Some(s) = self.unrolled() {
            let (u, v) = self.wrapped(u, v);
            return s.dn(u, v, nu, nv);
        }
        let c = self.derivatives(u, v, nu as usize, nv as usize)[nu as usize][nv as usize];
        crate::math::Vector::new(c[0], c[1], c[2])
    }
    fn u_iso(&self, u: f64) -> crate::geometry::curve::Curve {
        if let Some(s) = self.unrolled() {
            let (u, _) = self.wrapped(u, 0.0);
            return s.u_iso(u);
        }
        // Poles along v: the u-combinations of each column, in homogeneous form.
        let flat = flatten(&self.u.knots, &self.u.multiplicities);
        let (n, p) = (self.poles.len(), self.u.degree);
        let t = u.clamp(flat[p], flat[n]);
        let span = flat[..n].partition_point(|&k| k <= t).saturating_sub(1).clamp(p, n - 1);
        let basis = basis_derivatives(&flat, span, t, p, 0);
        let m = self.poles[0].len();
        let (mut poles, mut weights) = (Vec::with_capacity(m), Vec::with_capacity(m));
        for j in 0..m {
            let (mut c, mut w) = ([0.0; 3], 0.0);
            for (i, &b) in basis[0].iter().enumerate() {
                let ii = span - p + i;
                let wi = self.weight(ii, j) * b;
                let q = self.poles[ii][j];
                c[0] += wi * q.x;
                c[1] += wi * q.y;
                c[2] += wi * q.z;
                w += wi;
            }
            poles.push(Point::new(c[0] / w, c[1] / w, c[2] / w));
            weights.push(w);
        }
        let curve = BSplineCurve::new(poles, self.v.knots.clone(), self.v.multiplicities.clone(), self.v.degree);
        crate::geometry::curve::Curve::BSpline(if self.weights.is_some() { curve.with_weights(weights) } else { curve })
    }
    fn v_iso(&self, v: f64) -> crate::geometry::curve::Curve {
        if let Some(s) = self.unrolled() {
            let (_, v) = self.wrapped(0.0, v);
            return s.v_iso(v);
        }
        let swapped = self.swapped();
        swapped.u_iso(v)
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        let (a, b, _, _) = self.bounds();
        a + b - u
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        let (_, _, a, b) = self.bounds();
        a + b - v
    }
    fn u_reverse(&mut self) {
        let (a, b) = (self.u.knots[0], *self.u.knots.last().expect("knots"));
        let (n, p, m0) = (self.poles.len(), self.u.degree, self.u.multiplicities[0]);
        self.u.knots = self.u.knots.iter().rev().map(|&k| a + b - k).collect();
        self.u.multiplicities.reverse();
        if self.u.periodic {
            let old = reversed_pole(n, p, m0);
            self.poles = (0..n).map(|j| self.poles[old(j)].clone()).collect();
            if let Some(w) = &mut self.weights {
                *w = (0..n).map(|j| w[old(j)].clone()).collect();
            }
        } else {
            self.poles.reverse();
            if let Some(w) = &mut self.weights {
                w.reverse();
            }
        }
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }
    fn v_reverse(&mut self) {
        let (a, b) = (self.v.knots[0], *self.v.knots.last().expect("knots"));
        let (n, p, m0) = (self.poles[0].len(), self.v.degree, self.v.multiplicities[0]);
        self.v.knots = self.v.knots.iter().rev().map(|&k| a + b - k).collect();
        self.v.multiplicities.reverse();
        let periodic = self.v.periodic;
        let old = reversed_pole(n, p, m0);
        let flip = |row: &mut Vec<_>| {
            if periodic {
                *row = (0..n).map(|j| row[old(j)]).collect();
            } else {
                row.reverse();
            }
        };
        for row in &mut self.poles {
            flip(row);
        }
        if let Some(w) = &mut self.weights {
            for row in w {
                if periodic {
                    *row = (0..n).map(|j| row[old(j)]).collect();
                } else {
                    row.reverse();
                }
            }
        }
        self.patches = Default::default();
        self.flat_knots = Default::default();
        self.unrolled = Default::default();
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        let (a, b, c, d) = self.bounds();
        let size = self.poles.iter().flatten().fold(0.0f64, |m, p| m.max((*p - Point::ORIGIN).length())) * 2.0;
        let k = (self.u.degree.max(self.v.degree) as f64) * size.max(1e-300);
        (tolerance * (b - a) / k, tolerance * (d - c) / k)
    }
}

impl BSplineSurface {
    /// The same surface with u and v exchanged.
    fn swapped(&self) -> BSplineSurface {
        let (n, m) = (self.poles.len(), self.poles[0].len());
        let poles = (0..m).map(|j| (0..n).map(|i| self.poles[i][j]).collect()).collect();
        let weights = self.weights.as_ref().map(|w| (0..m).map(|j| (0..n).map(|i| w[i][j]).collect()).collect());
        BSplineSurface { poles, weights, u: self.v.clone(), v: self.u.clone(), patches: Default::default(), unrolled: Default::default(), flat_knots: Default::default() }
    }
}


impl Transformable for BSplineSurface {
    fn transformed(&self, t: &Transform) -> Self {
        if *t == Transform::IDENTITY {
            return self.clone();
        }
        let poles = self.poles.iter().map(|row| row.iter().map(|p| p.transformed(t)).collect()).collect();
        let map = |c: [f64; 3]| Point::new(c[0], c[1], c[2]).transformed(t).coords();
        let patches = self.patches.mapped(|s| s.map_points(|q| super::patches::map_homogeneous(q, map)));
        BSplineSurface { poles, weights: self.weights.clone(), u: self.u.clone(), v: self.v.clone(), patches, unrolled: Default::default(), flat_knots: Default::default() }
    }
}

impl AffineTransformable for BSplineSurface {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        // Affine maps commute with (rational) B-spline evaluation: map the poles.
        let poles = self.poles.iter().map(|row| row.iter().map(|p| p.affine_transformed(t)).collect()).collect();
        BSplineSurface { poles, weights: self.weights.clone(), u: self.u.clone(), v: self.v.clone(), patches: Default::default(), unrolled: Default::default(), flat_knots: Default::default() }
    }
}

#[cfg(test)]
mod fast_evaluation_tests {
    use super::*;
    use crate::geometry::{ParametricCurve, ParametricSurface};

    #[test]
    fn moving_a_point_and_its_tangent() {
        let poles: Vec<Point> = (0..7).map(|i| Point::new(i as f64, (i as f64 * 0.7).sin(), 0.0)).collect();
        let plain = BSplineCurve::new(poles, vec![0.0, 0.3, 0.5, 1.0], vec![4, 2, 1, 4], 3);
        let rational = plain.clone().with_weights((0..7).map(|i| 1.0 + 0.2 * (i % 3) as f64).collect());
        for mut curve in [plain, rational] {
            let (target, tangent) = (Point::new(3.0, 2.0, 1.0), crate::math::Vector::new(5.0, -1.0, 2.0));
            let (first, last) = (curve.poles()[0], curve.poles()[6]);
            curve.move_point_and_tangent(0.4, target, tangent, 1e-5, 0, 0).unwrap();
            let (p, d) = curve.d1(0.4);
            assert!(p.distance(target) < 1e-9 && (d - tangent).length() < 1e-9, "{p:?} {d:?}");
            // The end poles stay put.
            assert_eq!((curve.poles()[0], curve.poles()[6]), (first, last));
        }
        // Too near an end that must stay put.
        let mut curve = BSplineCurve::new((0..4).map(|i| Point::new(i as f64, 0.0, 0.0)).collect(), vec![0.0, 1.0], vec![4, 4], 3);
        assert!(curve.move_point_and_tangent(1e-7, Point::ORIGIN, crate::math::Vector::new(1.0, 1.0, 0.0), 1e-5, 0, 0).is_none());
    }

    #[test]
    fn fast_curve_derivatives_match_the_general_ones() {
        let poles: Vec<Point> = (0..7).map(|i| Point::new(i as f64, (i as f64 * 0.7).sin(), (i * i) as f64 * 0.1)).collect();
        let plain = BSplineCurve::new(poles, vec![0.0, 0.3, 0.5, 1.0], vec![4, 2, 1, 4], 3);
        let rational = plain.clone().with_weights((0..7).map(|i| 1.0 + 0.2 * (i % 3) as f64).collect());
        for curve in [plain, rational] {
            for i in -3..=43 {
                let t = i as f64 / 40.0;
                let general = curve.derivatives(t, 2);
                let (p, d1, d2) = curve.d2(t);
                let close = |a: [f64; 3], b: [f64; 3]| (0..3).all(|k| (a[k] - b[k]).abs() <= 1e-12 * (1.0 + b[k].abs()));
                assert!(close([p.x, p.y, p.z], general[0]) && close([d1.x, d1.y, d1.z], general[1]) && close([d2.x, d2.y, d2.z], general[2]), "at {t}");
                assert_eq!(curve.point(t), p);
                assert_eq!(curve.d1(t).1, d1);
            }
        }
    }

    #[test]
    fn segments_match_the_general_derivatives() {
        // Periodic, unclamped ends, degree 5, a knot of full multiplicity; rational and not;
        // 2D too; inside and outside the range.
        let poles: Vec<Point> = (0..9).map(|i| Point::new((i as f64 * 0.8).cos() * 3.0, (i as f64 * 0.8).sin() * 2.0, (i % 4) as f64 * 0.3)).collect();
        // (With whether to check outside the range too: not for the fitted-like curve, whose
        // short end span, extrapolated 100 times its length, magnifies rounding.)
        let curves = [
            BSplineCurve::new_periodic(poles.clone(), (0..=9).map(|i| i as f64 * 0.5).collect(), vec![1; 10], 3),
            BSplineCurve::new(poles.clone(), vec![0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0], vec![1, 1, 2, 3, 2, 1, 3], 3),
            BSplineCurve::new(poles[..6].to_vec(), vec![0.0, 1.0], vec![6, 6], 5),
            BSplineCurve::new(poles[..4].to_vec(), vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], 1),
            // Like a fitted curve: many poles, uneven knots, two nearly together.
            {
                let n = 200;
                let poles: Vec<Point> = (0..n).map(|i| Point::new(i as f64 * 0.1, (i as f64 * 0.37).sin(), (i as f64 * 0.11).cos())).collect();
                let mut knots: Vec<f64> = (0..=n - 3).map(|i| (i as f64 / (n - 3) as f64).powf(1.3)).collect();
                knots[50] = knots[49] + 1e-9;
                let mut mults = vec![1; knots.len()];
                mults[0] = 4;
                *mults.last_mut().unwrap() = 4;
                BSplineCurve::new(poles, knots, mults, 3)
            },
        ];
        let count = curves.len();
        for (k, plain) in curves.into_iter().enumerate() {
            let range = if k + 1 == count { 0..=40 } else { -5..=45 };
            let rational = plain.clone().with_weights((0..plain.poles().len()).map(|i| 0.7 + 0.2 * (i % 3) as f64).collect());
            for curve in [plain, rational] {
                let (a, b) = (curve.first_parameter(), curve.last_parameter());
                for i in range.clone() {
                    let t = a + (b - a) * i as f64 / 40.0;
                    let general = curve.derivatives(t, 2);
                    let (p, d1, d2) = curve.d2(t);
                    let scale = 1.0 + general[0].iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    // (Each order to its own size: a short span's derivatives are large.)
                    let close = |a: [f64; 3], b: [f64; 3]| {
                        let size = 1.0 + b.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                        (0..3).all(|k| (a[k] - b[k]).abs() <= 1e-10 * size)
                    };
                    assert!(close([p.x, p.y, p.z], general[0]) && close([d1.x, d1.y, d1.z], general[1]) && close([d2.x, d2.y, d2.z], general[2]), "at {t}: {p:?} {d1:?} {d2:?} vs {general:?}");
                    assert!(curve.point(t).distance(p) <= 1e-14 * scale);
                }
                let flat: BSplineCurve<Point2> = BSplineCurve::new(
                    curve.poles().iter().take(4).map(|p| Point2::new(p.x, p.y)).collect(),
                    vec![0.0, 1.0],
                    vec![4, 4],
                    3,
                );
                let q = flat.point(0.3);
                let g = flat.derivatives(0.3, 0)[0];
                assert!((q.x - g[0]).abs() <= 1e-12 && (q.y - g[1]).abs() <= 1e-12);
            }
        }
    }

    #[test]
    fn fast_point_and_derivatives_match_the_general_ones() {
        // Degrees 3 x 2, an inner knot of multiplicity 2 along u, rational and not.
        let u = KnotVector::new(vec![0.0, 0.4, 1.0], vec![4, 2, 4], 3);
        let v = KnotVector::new(vec![-1.0, 0.0, 2.0], vec![3, 1, 3], 2);
        let poles: Vec<Vec<Point>> = (0..6)
            .map(|i| (0..4).map(|j| Point::new(i as f64, j as f64 + 0.3 * (i as f64).sin(), ((i * j) as f64).cos())).collect())
            .collect();
        let plain = BSplineSurface::new(poles, u, v);
        let weights: Vec<Vec<f64>> = (0..6).map(|i| (0..4).map(|j| 1.0 + 0.1 * ((i + 2 * j) % 5) as f64).collect()).collect();
        let rational = plain.clone().with_weights(weights);
        for surface in [plain, rational] {
            for i in -2..=22 {
                for j in -2..=22 {
                    let (s, t) = (i as f64 / 20.0, -1.0 + 3.0 * j as f64 / 20.0);
                    let general = surface.derivatives(s, t, 1, 1);
                    let (p, du, dv) = surface.d1(s, t);
                    let close = |a: [f64; 3], b: [f64; 3]| (0..3).all(|k| (a[k] - b[k]).abs() <= 1e-12 * (1.0 + b[k].abs()));
                    assert!(close([p.x, p.y, p.z], general[0][0]), "point at {s} {t}");
                    assert!(close([du.x, du.y, du.z], general[1][0]), "du at {s} {t}");
                    assert!(close([dv.x, dv.y, dv.z], general[0][1]), "dv at {s} {t}");
                    assert!(surface.point(s, t).distance(p) <= 1e-14 * (1.0 + p.distance(Point::ORIGIN)));
                    let general = surface.derivatives(s, t, 2, 2);
                    let d2 = surface.d2(s, t);
                    let c = |v: crate::math::Vector| [v.x, v.y, v.z];
                    assert!(close(c(d2.duu), general[2][0]) && close(c(d2.dvv), general[0][2]) && close(c(d2.duv), general[1][1]), "second derivatives at {s} {t}");
                    assert!(close(c(d2.du), general[1][0]) && close(c(d2.dv), general[0][1]) && close([d2.point.x, d2.point.y, d2.point.z], general[0][0]));
                }
            }
        }
    }

    #[test]
    fn patches_match_the_general_derivatives() {
        // Unclamped ends, a knot of full multiplicity, degree 5 (the larger arrays), a
        // single span; rational and not; inside and outside the bounds.
        let cases = [
            (KnotVector::new(vec![0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0], vec![1, 1, 2, 3, 1, 1, 2], 3), KnotVector::new(vec![0.0, 1.0], vec![6, 6], 5)),
            (KnotVector::new(vec![2.0, 3.0], vec![2, 2], 1), KnotVector::new(vec![0.0, 0.25, 0.5, 1.0], vec![3, 2, 1, 3], 2)),
        ];
        for (u, v) in cases {
            let count = |kv: &KnotVector| kv.multiplicities.iter().sum::<usize>() - kv.degree - 1;
            let (nu, nv) = (count(&u), count(&v));
            let poles: Vec<Vec<Point>> = (0..nu)
                .map(|i| (0..nv).map(|j| Point::new(i as f64 + 0.2 * (j as f64).sin(), j as f64, ((i * j) as f64 * 0.7).cos())).collect())
                .collect();
            let plain = BSplineSurface::new(poles, u, v);
            let weights: Vec<Vec<f64>> = (0..nu).map(|i| (0..nv).map(|j| 0.6 + 0.15 * ((i + 3 * j) % 5) as f64).collect()).collect();
            let rational = plain.clone().with_weights(weights);
            for surface in [plain, rational] {
                let (u0, u1, v0, v1) = surface.bounds();
                for i in -3..=33 {
                    for j in -3..=33 {
                        let (s, t) = (u0 + (u1 - u0) * i as f64 / 30.0, v0 + (v1 - v0) * j as f64 / 30.0);
                        let general = surface.derivatives(s, t, 2, 2);
                        let d2 = surface.d2(s, t);
                        let scale = 1.0 + general[0][0].iter().fold(0.0f64, |m, x| m.max(x.abs()));
                        let close = |a: crate::math::Vector, b: [f64; 3]| (a.x - b[0]).abs().max((a.y - b[1]).abs()).max((a.z - b[2]).abs()) <= 1e-10 * scale;
                        let p = d2.point;
                        assert!(close(crate::math::Vector::new(p.x, p.y, p.z), general[0][0]), "point at {s} {t}");
                        assert!(close(d2.du, general[1][0]) && close(d2.dv, general[0][1]), "first derivatives at {s} {t}");
                        assert!(close(d2.duu, general[2][0]) && close(d2.dvv, general[0][2]) && close(d2.duv, general[1][1]), "second derivatives at {s} {t}");
                        let (q, du, _) = surface.d1(s, t);
                        assert!(q.distance(p) <= 1e-14 * scale && (du - d2.du).length() <= 1e-12 * scale);
                    }
                }
            }
        }
    }

    /// Rows evaluated at once agree with single evaluations to the last bit (in and out of
    /// bounds, rational and not, low and high degree, periodic).
    #[test]
    fn rows_match_single_evaluations() {
        let cases = [
            (KnotVector::new(vec![0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0], vec![1, 1, 2, 3, 1, 1, 2], 3), KnotVector::new(vec![0.0, 1.0], vec![6, 6], 5)),
            (KnotVector::new(vec![2.0, 3.0], vec![2, 2], 1), KnotVector::new(vec![0.0, 0.25, 0.5, 1.0], vec![3, 2, 1, 3], 2)),
            (KnotVector::new(vec![0.0, 0.5, 1.0], vec![4, 2, 4], 3), KnotVector::new(vec![0.0, 0.3, 0.6, 1.0], vec![4, 1, 1, 4], 3)),
        ];
        let bits = |p: Point| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        let vbits = |p: crate::math::Vector| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        for (u, v) in cases {
            let count = |kv: &KnotVector| kv.multiplicities.iter().sum::<usize>() - kv.degree - 1;
            let (nu, nv) = (count(&u), count(&v));
            let poles: Vec<Vec<Point>> = (0..nu)
                .map(|i| (0..nv).map(|j| Point::new(i as f64 + 0.2 * (j as f64).sin(), j as f64, ((i * j) as f64 * 0.7).cos())).collect())
                .collect();
            let plain = BSplineSurface::new(poles, u, v);
            let weights: Vec<Vec<f64>> = (0..nu).map(|i| (0..nv).map(|j| 0.6 + 0.15 * ((i + 3 * j) % 5) as f64).collect()).collect();
            let rational = plain.clone().with_weights(weights);
            // Periodic in u: 6 poles along u, the net's rows reused.
            let around = KnotVector::periodic(vec![0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0], vec![1; 7], 3);
            let periodic_poles: Vec<Vec<Point>> = (0..6).map(|i| (0..nv).map(|j| Point::new((i as f64).cos() * 3.0, (i as f64).sin() * 3.0, j as f64 + 0.1 * i as f64)).collect()).collect();
            let periodic = BSplineSurface::new(periodic_poles, around, plain.v.clone());
            for surface in [plain, rational, periodic] {
                let (u0, u1, v0, v1) = surface.bounds();
                let us: Vec<f64> = (-3..=33).map(|i| u0 + (u1 - u0) * i as f64 / 30.0).collect();
                for j in -3..=33 {
                    let t = v0 + (v1 - v0) * j as f64 / 30.0;
                    let rows = surface.d1_along_u(&us, t);
                    let points = surface.points_along_u(&us, t);
                    for ((&s, (p, du, dv)), q) in us.iter().zip(rows).zip(points) {
                        let (p1, du1, dv1) = surface.d1(s, t);
                        assert_eq!((bits(p), vbits(du), vbits(dv)), (bits(p1), vbits(du1), vbits(dv1)), "d1 at {s} {t}");
                        assert_eq!(bits(q), bits(surface.point(s, t)), "point at {s} {t}");
                    }
                }
            }
        }
    }
}

/// The interval of `u` among the increasing `values` (1-based, as OCCT's
/// `BSplCLib::LocateParameter`): the last `i` with `values[i] <= u` (a value equal to `u`
/// within rounding counting as below), kept within `first..last` and past empty intervals.
fn locate_in(values: &[f64], u: f64, first: usize, last: usize) -> usize {
    let at = |i: usize| values[i - 1];
    let eps = (values[values.len() - 1].abs().min(u.abs()) * f64::EPSILON).max(f64::MIN_POSITIVE);
    let mut index = values.partition_point(|&x| x <= u);
    if index < values.len() && (u - at(index + 1)).abs() <= eps {
        index += 1;
    }
    index = index.clamp(first, last - 1);
    while index != last - 1 && index < values.len() && (at(index + 1) - at(index)).abs() <= eps {
        index += 1;
    }
    index
}
