//! A curve seen through its normalised arc length `s = length(u_first, u)/L ∈ [0, 1]`.
//!
//! On construction the curve's smooth spans are each cut into a few pieces and the arc
//! length at every cut is tabulated (to the tolerance); `s(u)` then integrates only from the
//! cut below `u`, and `u(s)` finds the piece holding `s·L` in the table and solves inside it
//! by Newton's method on the arc length. Results depend on the arguments only.

use crate::error::{Error, Result};
use crate::geometry::abscissa::signed_length;
use crate::geometry::adaptor::CurveAdaptor;
use crate::geometry::Continuity;

/// Pieces each smooth span is cut into for the table.
const PIECES_PER_SPAN: usize = 4;

/// The curve held: borrowed or owned.
enum Held<'a> {
    Borrowed(&'a dyn CurveAdaptor),
    Owned(Box<dyn CurveAdaptor>),
}

impl Held<'_> {
    fn get(&self) -> &dyn CurveAdaptor {
        match self {
            Held::Borrowed(c) => *c,
            Held::Owned(c) => c.as_ref(),
        }
    }
}

/// A curve reparametrised by normalised arc length, possibly restricted to a sub-range of
/// `s` (the abscissa stays that of the whole curve).
pub(crate) struct CurvlinFunc<'a> {
    original: Held<'a>,
    trimmed: Option<Box<dyn CurveAdaptor>>,
    tolerance: f64,
    length: f64,
    first: f64,
    last: f64,
    /// Parameters of the cuts of the original curve, increasing.
    cuts: Vec<f64>,
    /// Arc length from the first parameter to each cut.
    lengths: Vec<f64>,
}

impl<'a> CurvlinFunc<'a> {
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, tolerance: f64) -> Self {
        Self::build(Held::Borrowed(curve), tolerance)
    }

    pub(crate) fn owning(curve: Box<dyn CurveAdaptor>, tolerance: f64) -> CurvlinFunc<'static> {
        CurvlinFunc::build(Held::Owned(curve), tolerance)
    }

    fn build(original: Held<'a>, tolerance: f64) -> Self {
        let mut f = CurvlinFunc { original, trimmed: None, tolerance, length: 0.0, first: 0.0, last: 1.0, cuts: Vec::new(), lengths: Vec::new() };
        f.tabulate();
        f
    }

    /// The table of arc lengths at the cuts, and the total length.
    fn tabulate(&mut self) {
        let curve = self.original.get();
        let spans = curve.intervals(Continuity::CN);
        let mut cuts = Vec::with_capacity(spans.len() * PIECES_PER_SPAN);
        for w in spans.windows(2) {
            for k in 0..PIECES_PER_SPAN {
                cuts.push(w[0] + (w[1] - w[0]) * k as f64 / PIECES_PER_SPAN as f64);
            }
        }
        cuts.push(*spans.last().unwrap());
        // (Each piece to a share of the tolerance, so that the sum keeps to it.)
        let piece_tolerance = self.tolerance / (cuts.len() as f64);
        let mut lengths = Vec::with_capacity(cuts.len());
        let mut total = 0.0;
        lengths.push(0.0);
        for w in cuts.windows(2) {
            total += signed_length(curve, w[0], w[1], piece_tolerance).abs();
            lengths.push(total);
        }
        self.cuts = cuts;
        self.lengths = lengths;
        self.length = total;
    }

    /// Later lengths and parameters are computed to `tolerance` (the total length is kept).
    pub(crate) fn set_tolerance(&mut self, tolerance: f64) {
        self.tolerance = tolerance;
    }

    /// The first `s` of the current range.
    pub(crate) fn first_parameter(&self) -> f64 {
        self.first
    }

    /// The last `s` of the current range.
    pub(crate) fn last_parameter(&self) -> f64 {
        self.last
    }

    /// The curve currently viewed (the original, or its restriction by the last `trim`).
    pub(crate) fn curve(&self) -> &dyn CurveAdaptor {
        match &self.trimmed {
            Some(t) => t.as_ref(),
            None => self.original.get(),
        }
    }

    pub(crate) fn nb_intervals(&self, continuity: Continuity) -> usize {
        self.curve().intervals(continuity).len() - 1
    }

    /// The current curve's breakpoints of `continuity` as abscissae (the ends exactly the
    /// current range's).
    pub(crate) fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        let breaks = self.curve().intervals(continuity);
        let n = breaks.len();
        breaks
            .iter()
            .enumerate()
            .map(|(i, &u)| if i == 0 { self.first } else if i + 1 == n { self.last } else { self.get_s_parameter(u) })
            .collect()
    }

    /// Restricts the view to the abscissae `[first, last]` (relative to the original curve);
    /// nothing happens for a range shorter than `tolerance`.
    pub(crate) fn trim(&mut self, first: f64, last: f64, tolerance: f64) -> Result<()> {
        if first < 0.0 || last > 1.0 {
            return Err(Error::InvalidInput(format!("abscissae [{first}, {last}] outside [0, 1]")));
        }
        if last - first < tolerance {
            return Ok(());
        }
        let original = self.original.get();
        let u1 = self.get_u_parameter_on(original, first)?;
        let u2 = self.get_u_parameter_on(original, last)?;
        self.trimmed = Some(original.trimmed(u1, u2, tolerance));
        self.first = first;
        self.last = last;
        Ok(())
    }

    /// The total length of the original curve.
    pub(crate) fn get_length(&self) -> f64 {
        self.length
    }

    /// The abscissa of the parameter `u`.
    pub(crate) fn get_s_parameter(&self, u: f64) -> f64 {
        let curve = self.original.get();
        // The cut at or below `u` (the first one for a `u` before them).
        let i = self.cuts.partition_point(|&c| c <= u).saturating_sub(1).min(self.cuts.len() - 1);
        let base = self.cuts[i];
        let along = if u == base { 0.0 } else { signed_length(curve, base, u, 0.5 * self.tolerance) };
        (self.lengths[i] + along) / self.length
    }

    /// The parameter of the abscissa `s ∈ [0, 1]` of the original curve.
    pub(crate) fn get_u_parameter(&self, s: f64) -> Result<f64> {
        self.get_u_parameter_on(self.curve(), s)
    }

    /// The parameter of the abscissa `s ∈ [0, 1]`, evaluated along `curve` (which shares the
    /// original curve's parametrisation).
    pub(crate) fn get_u_parameter_on(&self, curve: &dyn CurveAdaptor, s: f64) -> Result<f64> {
        if !(0.0..=1.0).contains(&s) {
            return Err(Error::InvalidInput(format!("abscissa {s} outside [0, 1]")));
        }
        let n = self.cuts.len();
        if s == 0.0 {
            return Ok(self.cuts[0]);
        }
        if s == 1.0 {
            return Ok(self.cuts[n - 1]);
        }
        let target = s * self.length;
        // The piece holding the target.
        let i = self.lengths.partition_point(|&l| l <= target).clamp(1, n - 1) - 1;
        let (a, b) = (self.cuts[i], self.cuts[i + 1]);
        let remaining = target - self.lengths[i];
        let piece = self.lengths[i + 1] - self.lengths[i];
        if remaining <= 0.0 {
            return Ok(a);
        }
        if remaining >= piece {
            return Ok(b);
        }
        solve(curve, a, b, remaining, piece, 0.5 * self.tolerance).ok_or_else(|| Error::Failed(format!("no parameter found at the abscissa {s}")))
    }

    /// The point (`order` 0) or a derivative with respect to the abscissa (1, 2) at `s`;
    /// `None` for higher orders.
    pub(crate) fn eval_case1(&self, s: f64, order: u32) -> Result<Option<[f64; 3]>> {
        if order > 2 {
            return Ok(None);
        }
        let u = self.get_u_parameter(s)?;
        let curve = self.curve();
        let (p, d1, d2) = curve.d2(u);
        let l = self.length;
        let out = match order {
            0 => [p.x, p.y, p.z],
            1 => {
                let v = d1 * (l / d1.length());
                [v.x, v.y, v.z]
            }
            _ => {
                let n = d1.length();
                let du = l / n;
                let ddu = -l * d1.dot(d2) * du / (n * n * n);
                let v = d2 * (du * du) + d1 * ddu;
                [v.x, v.y, v.z]
            }
        };
        Ok(Some(out))
    }
}

/// The parameter between `a < b` at the arc length `remaining` from `a` (the piece being
/// `piece` long), to `tolerance`.
fn solve(curve: &dyn CurveAdaptor, a: f64, b: f64, remaining: f64, piece: f64, tolerance: f64) -> Option<f64> {
    let quad = 0.25 * tolerance;
    let (mut lo, mut hi) = (a, b);
    let mut u = a + (b - a) * (remaining / piece);
    let mut along = signed_length(curve, a, u, quad);
    for _ in 0..100 {
        let h = along - remaining;
        if h.abs() <= 0.5 * tolerance {
            return Some(u);
        }
        if h < 0.0 {
            lo = u;
        } else {
            hi = u;
        }
        let speed = curve.d1(u).1.length();
        let mut next = if speed > 0.0 && speed.is_finite() { u - h / speed } else { f64::NAN };
        if !next.is_finite() || next <= lo || next >= hi {
            next = 0.5 * (lo + hi);
        }
        if next == u || next == lo || next == hi {
            return Some(u);
        }
        along += signed_length(curve, u, next, quad);
        u = next;
    }
    Some(u)
}
