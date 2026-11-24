//! The plane section of what a moving solid leaves of a blank: a region of the plane bounded
//! by pieces of the blank's own outline and of curves on which the swept set's boundary must
//! lie (the candidates), each piece kept where it bounds the region.
//!
//! The candidates are given as functions of a parameter; a point of one lies on the boundary
//! of the swept set unless the set covers it, which `Cover::depth` says (how deep the point
//! lies in the moving solid at its deepest, positive inside). So each candidate is split where
//! that changes or where it leaves the blank, at known junctions (`breaks`: parameters where
//! another candidate ends on it, touching it, where no crossing shows the change), and the
//! pieces kept are chained end to end into loops, the region on their left.

use crate::error::{Error, Result};

pub(super) type P2 = [f64; 2];

pub(super) fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

pub(super) fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

pub(super) fn norm(a: P2) -> f64 {
    a[0].hypot(a[1])
}

/// A curve on which the swept set's boundary may lie, `point(s)` for s in `range`; `closed`
/// if it comes back to where it starts (a periodic parameter). `breaks` are parameters where
/// another candidate ends on it.
pub(super) struct Candidate<'a> {
    pub point: Box<dyn Fn(f64) -> P2 + Sync + 'a>,
    pub range: (f64, f64),
    pub closed: bool,
    pub breaks: Vec<f64>,
}

impl Candidate<'_> {
    /// The point at `s`, a closed candidate's parameter taken round its period.
    pub fn at(&self, s: f64) -> P2 {
        if self.closed {
            let (a, b) = self.range;
            let period = b - a;
            let mut s = s;
            while s > b {
                s -= period;
            }
            while s < a {
                s += period;
            }
            (self.point)(s)
        } else {
            (self.point)(s)
        }
    }
}

/// How deep a point lies in the moving solid at its deepest: positive where the swept set
/// covers it (within its interior), zero on its boundary, negative outside.
pub(super) trait Cover: Sync {
    fn depth(&self, p: P2) -> f64;
}

/// The blank's outline in the section's plane: a disc, or a rectangle of the plane's two
/// coordinates.
#[derive(Debug, Clone, Copy)]
pub(super) enum Outline {
    Disc { centre: P2, radius: f64 },
    Rectangle { min: P2, max: P2 },
}

impl Outline {
    /// How far inside the outline `p` lies (negative outside).
    pub fn inside(&self, p: P2) -> f64 {
        match *self {
            Outline::Disc { centre, radius } => radius - norm(sub(p, centre)),
            Outline::Rectangle { min, max } => (p[0] - min[0]).min(max[0] - p[0]).min(p[1] - min[1]).min(max[1] - p[1]),
        }
    }

    fn size(&self) -> f64 {
        match *self {
            Outline::Disc { radius, .. } => 2.0 * radius,
            Outline::Rectangle { min, max } => norm(sub(max, min)),
        }
    }
}

/// One piece of a loop of the region's boundary, run with the region on its left.
#[derive(Debug, Clone, Copy)]
pub(super) enum Piece {
    /// Candidate `index` from parameter `from` to `to` (either the larger; past its range for
    /// a closed one, taken round).
    Candidate { index: usize, from: f64, to: f64 },
    /// The blank's circle from angle `from` to `to` (radians about its centre from the
    /// plane's first axis; decreasing, clockwise).
    Arc { centre: P2, radius: f64, from: f64, to: f64 },
    /// A straight piece of the blank's rectangle.
    Segment { from: P2, to: P2 },
}

impl Piece {
    pub fn point(&self, candidates: &[Candidate], f: f64) -> P2 {
        match *self {
            Piece::Candidate { index, from, to } => candidates[index].at(from + (to - from) * f),
            Piece::Arc { centre, radius, from, to } => {
                let a = from + (to - from) * f;
                [centre[0] + radius * a.cos(), centre[1] + radius * a.sin()]
            }
            Piece::Segment { from, to } => [from[0] + (to[0] - from[0]) * f, from[1] + (to[1] - from[1]) * f],
        }
    }
}

/// A loop of the region's boundary: its pieces in order, each one's end the next one's
/// start (`junctions[i]` where piece i starts).
pub(super) struct Loop {
    pub pieces: Vec<Piece>,
    pub junctions: Vec<P2>,
}

impl Loop {
    /// Points along the loop, in order, for areas and containment.
    pub fn polygon(&self, candidates: &[Candidate]) -> Vec<P2> {
        let mut points = Vec::new();
        for piece in &self.pieces {
            for i in 0..64 {
                points.push(piece.point(candidates, i as f64 / 64.0));
            }
        }
        points
    }
}

pub(super) fn signed_area(polygon: &[P2]) -> f64 {
    let n = polygon.len();
    (0..n).map(|i| cross(polygon[i], polygon[(i + 1) % n])).sum::<f64>() / 2.0
}

pub(super) fn contains(polygon: &[P2], p: P2) -> bool {
    let n = polygon.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
            if x > p[0] {
                inside = !inside;
            }
        }
    }
    inside
}

/// The region's boundary loops: of the blank within `outline`, what the swept set does not
/// cover, its boundary lying on the outline and on `candidates`.
pub(super) fn region(candidates: &[Candidate], cover: &dyn Cover, outline: Outline) -> Result<Vec<Loop>> {
    let size = outline.size();
    // (How deep a point of a candidate must lie to count as covered: its own pose leaves it
    // at depth zero.)
    let eps = 1e-11 * size;
    let in_region = |p: P2| outline.inside(p) > 0.0 && cover.depth(p) < 0.0;
    let mut pieces: Vec<Piece> = Vec::new();
    let kept = |p: P2| outline.inside(p) > 0.0 && cover.depth(p) <= eps;
    // (Each candidate's on its own: taken together.)
    let kept_on = crate::parallel::map(candidates, |candidate| kept_intervals(candidate, &kept, size));
    for (index, intervals) in kept_on.into_iter().enumerate() {
        for (from, to) in intervals {
            let piece = Piece::Candidate { index, from, to };
            // Run with the region on its left (dropped if it has the region on neither side
            // or both: where the swept set only touches itself).
            let (p, q) = (piece.point(candidates, 0.5 - 1e-3), piece.point(candidates, 0.5 + 1e-3));
            let m = piece.point(candidates, 0.5);
            let tangent = sub(q, p);
            let length = norm(tangent);
            if !(length > 0.0) {
                continue;
            }
            let delta = 1e-6 * size;
            let left = [m[0] - tangent[1] / length * delta, m[1] + tangent[0] / length * delta];
            let right = [m[0] + tangent[1] / length * delta, m[1] - tangent[0] / length * delta];
            match (in_region(left), in_region(right)) {
                (true, false) => pieces.push(piece),
                (false, true) => pieces.push(Piece::Candidate { index, from: to, to: from }),
                _ => {}
            }
        }
    }
    // (Candidates may coincide where the tool's symmetry makes two of its corners or sides
    // run along the same curve: such a piece is kept once.)
    let mut unique: Vec<Piece> = Vec::with_capacity(pieces.len());
    let mut samples: Vec<Samples> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        if !unique.iter().zip(&samples).any(|(kept, along)| lies_on(candidates, &piece, kept, along, 1e-9 * size)) {
            samples.push(Samples::of(candidates, &piece));
            unique.push(piece);
        }
    }
    let mut pieces = unique;
    // The outline where the swept set does not cover it, run counter-clockwise.
    let outside = |p: P2| cover.depth(p) < 0.0;
    match outline {
        Outline::Disc { centre, radius } => {
            let circle = Candidate {
                point: Box::new(move |a: f64| [centre[0] + radius * a.cos(), centre[1] + radius * a.sin()]),
                range: (0.0, std::f64::consts::TAU),
                closed: true,
                breaks: Vec::new(),
            };
            for (from, to) in kept_intervals(&circle, &outside, size) {
                pieces.push(Piece::Arc { centre, radius, from, to });
            }
        }
        Outline::Rectangle { min, max } => {
            let corners = [min, [max[0], min[1]], max, [min[0], max[1]]];
            for k in 0..4 {
                let (a, b) = (corners[k], corners[(k + 1) % 4]);
                let side = Candidate {
                    point: Box::new(move |f: f64| [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]),
                    range: (0.0, 1.0),
                    closed: false,
                    breaks: Vec::new(),
                };
                for (from, to) in kept_intervals(&side, &outside, size) {
                    pieces.push(Piece::Segment { from: side.at(from), to: side.at(to) });
                }
            }
        }
    }
    chain(candidates, pieces, size)
}

/// A piece's points, even in its parameter, and the box holding them: what `lies_on` searches
/// for a nearest point.
struct Samples {
    points: Vec<P2>,
    lo: P2,
    hi: P2,
}

impl Samples {
    const COUNT: usize = 2048;

    fn of(candidates: &[Candidate], piece: &Piece) -> Samples {
        let n = Self::COUNT;
        let points: Vec<P2> = (0..=n).map(|i| piece.point(candidates, i as f64 / n as f64)).collect();
        let (lo, hi) = points.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| {
            ([lo[0].min(p[0]), lo[1].min(p[1])], [hi[0].max(p[0]), hi[1].max(p[1])])
        });
        Samples { points, lo, hi }
    }
}

/// True if all of piece `a` lies within `tolerance` of piece `b` (whose samples are `along`):
/// sampled, each sample's nearest point on `b` found among `b`'s samples and then by golden
/// section.
fn lies_on(candidates: &[Candidate], a: &Piece, b: &Piece, along: &Samples, tolerance: f64) -> bool {
    let n = Samples::COUNT;
    let (lo, hi) = (along.lo, along.hi);
    (0..=16).all(|k| {
        let p = a.point(candidates, k as f64 / 16.0);
        if p[0] < lo[0] - tolerance || p[1] < lo[1] - tolerance || p[0] > hi[0] + tolerance || p[1] > hi[1] + tolerance {
            return false;
        }
        let distances: Vec<f64> = along.points.iter().map(|&q| norm(sub(q, p))).collect();
        let nearest = (0..=n).min_by(|&i, &j| distances[i].total_cmp(&distances[j])).expect("samples");
        let distance = |f: f64| norm(sub(b.point(candidates, f), p));
        let (mut x, mut y) = ((nearest.max(1) - 1) as f64 / n as f64, ((nearest + 1).min(n)) as f64 / n as f64);
        let g = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..60 {
            let (x1, x2) = (y - g * (y - x), x + g * (y - x));
            if distance(x1) < distance(x2) {
                y = x2;
            } else {
                x = x1;
            }
        }
        distance(0.5 * (x + y)) <= tolerance
    })
}

/// The parameter intervals of `candidate` where `kept` holds, split at its breaks.
fn kept_intervals(candidate: &Candidate, kept: &(dyn Fn(P2) -> bool + Sync), size: f64) -> Vec<(f64, f64)> {
    let (a, b) = candidate.range;
    let span = b - a;
    // Samples: even in the parameter, then halved where they lie far apart, the breaks
    // among them (and, on an open candidate, its ends, which are junctions too).
    let mut params: Vec<(f64, bool)> = (0..=256).map(|i| (a + span * i as f64 / 256.0, false)).collect();
    if !candidate.closed {
        params[0].1 = true;
        params[256].1 = true;
    } else {
        params.pop();
    }
    for &s in &candidate.breaks {
        if s > a && s < b {
            params.push((s, true));
        } else if candidate.closed && (s == a || s == b) {
            params[0] = (a, true);
        }
    }
    params.sort_by(|x, y| x.0.total_cmp(&y.0));
    params.dedup_by(|x, y| {
        if (x.0 - y.0).abs() <= 1e-14 * span.abs().max(1.0) {
            y.1 |= x.1;
            true
        } else {
            false
        }
    });
    let h = size / 400.0;
    for _ in 0..12 {
        let points: Vec<P2> = params.iter().map(|&(s, _)| candidate.at(s)).collect();
        let n = params.len();
        let count = if candidate.closed { n } else { n - 1 };
        let mut refined = Vec::with_capacity(2 * n);
        let mut any = false;
        for i in 0..count {
            refined.push(params[i]);
            let j = (i + 1) % n;
            let next = if j == 0 { params[0].0 + span } else { params[j].0 };
            let both_breaks = params[i].1 && params[j].1;
            if norm(sub(points[j], points[i])) > h || both_breaks {
                refined.push((0.5 * (params[i].0 + next), false));
                any = true;
            }
        }
        if !candidate.closed {
            refined.push(params[n - 1]);
        }
        params = refined;
        if !any {
            break;
        }
    }
    // Each sample's state, the breaks taking their neighbours' (a junction lies on both).
    let states: Vec<Option<bool>> =
        crate::parallel::map(&params, |&(s, is_break)| if is_break { None } else { Some(kept(candidate.at(s))) });
    let n = params.len();
    let defined: Vec<usize> = (0..n).filter(|&i| states[i].is_some()).collect();
    if defined.is_empty() {
        return Vec::new();
    }
    let state = |i: usize| states[i].expect("defined");
    // The transition between defined samples i and j (in order, j after i, maybe round the
    // end of a closed candidate): at a break between them, else where `kept` changes.
    let transition = |i: usize, j: usize| -> f64 {
        let (si, mut sj) = (params[i].0, params[j].0);
        if sj <= si {
            sj += span;
        }
        let mut k = (i + 1) % n;
        while k != j {
            if params[k].1 {
                let mut s = params[k].0;
                if s < si {
                    s += span;
                }
                return s;
            }
            k = (k + 1) % n;
        }
        let (mut lo, mut hi) = (si, sj);
        let at_lo = state(i);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if kept(candidate.at(mid)) == at_lo {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    // The transitions between neighbouring defined samples whose states differ (each found
    // on its own: taken together, keyed by the first sample).
    let pairs: Vec<(usize, usize)> = (0..defined.len())
        .filter(|&w| candidate.closed || w + 1 < defined.len())
        .map(|w| (defined[w], defined[(w + 1) % defined.len()]))
        .filter(|&(i, j)| state(i) != state(j))
        .collect();
    let mut found = vec![f64::NAN; n];
    for (&(i, _), s) in pairs.iter().zip(crate::parallel::map(&pairs, |&(i, j)| transition(i, j))) {
        found[i] = s;
    }
    let transition = |i: usize, _: usize| found[i];
    let mut intervals = Vec::new();
    if candidate.closed {
        // Start from a sample not kept, if any.
        let Some(start) = defined.iter().position(|&i| !state(i)) else {
            return vec![(a, b)];
        };
        let order: Vec<usize> = (0..defined.len()).map(|k| defined[(start + k) % defined.len()]).collect();
        let mut open: Option<f64> = None;
        let first = params[order[0]].0;
        for w in 0..order.len() {
            let (i, j) = (order[w], order[(w + 1) % order.len()]);
            if state(i) != state(j) {
                // (Measured from the start of the walk, increasing.)
                let mut s = transition(i, j);
                while s < first {
                    s += span;
                }
                if state(j) {
                    open = Some(s);
                } else if let Some(from) = open.take() {
                    let mut to = s;
                    while to < from {
                        to += span;
                    }
                    intervals.push((from, to));
                }
            }
        }
    } else {
        let mut open: Option<f64> = if state(defined[0]) { Some(a) } else { None };
        for w in 0..defined.len() - 1 {
            let (i, j) = (defined[w], defined[w + 1]);
            if state(i) != state(j) {
                let s = transition(i, j);
                if state(j) {
                    open = Some(s);
                } else if let Some(from) = open.take() {
                    intervals.push((from, s));
                }
            }
        }
        if let Some(from) = open {
            intervals.push((from, b));
        }
    }
    intervals.retain(|&(from, to)| to - from > 1e-12 * span.abs());
    intervals
}

/// The pieces chained end to end into loops.
fn chain(candidates: &[Candidate], pieces: Vec<Piece>, size: f64) -> Result<Vec<Loop>> {
    let tolerance = 1e-7 * size;
    let ends: Vec<(P2, P2)> = pieces.iter().map(|p| (p.point(candidates, 0.0), p.point(candidates, 1.0))).collect();
    // The direction a piece leaves its start and arrives at its end.
    let leaving = |k: usize| {
        let d = sub(pieces[k].point(candidates, 1e-4), ends[k].0);
        let l = norm(d);
        [d[0] / l, d[1] / l]
    };
    let arriving = |k: usize| {
        let d = sub(ends[k].1, pieces[k].point(candidates, 1.0 - 1e-4));
        let l = norm(d);
        [d[0] / l, d[1] / l]
    };
    let mut used = vec![false; pieces.len()];
    let mut loops = Vec::new();
    for first in 0..pieces.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut order = vec![first];
        let mut current = first;
        loop {
            let end = ends[current].1;
            // What may follow: the unused pieces starting here, and the first (closing the
            // loop); where several meet, the one turning most sharply right, which keeps the
            // region on the left.
            let mut next: Vec<usize> = (0..pieces.len()).filter(|&k| !used[k] && norm(sub(ends[k].0, end)) <= tolerance).collect();
            if norm(sub(ends[first].0, end)) <= tolerance {
                next.push(first);
            }
            let into = arriving(current);
            let turn = |k: usize| {
                let out = leaving(k);
                cross(into, out).atan2(into[0] * out[0] + into[1] * out[1])
            };
            let Some(chosen) = next.iter().copied().min_by(|&x, &y| turn(x).total_cmp(&turn(y))) else {
                return Err(Error::Failed(format!(
                    "the section's boundary does not close: a piece ends at ({:.9}, {:.9}) and none starts there",
                    end[0], end[1]
                )));
            };
            if chosen == first {
                break;
            }
            used[chosen] = true;
            order.push(chosen);
            current = chosen;
        }
        let n = order.len();
        // Each junction where two pieces meet: halfway between their ends.
        let junctions: Vec<P2> = (0..n)
            .map(|i| {
                let (a, b) = (ends[order[(i + n - 1) % n]].1, ends[order[i]].0);
                [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]
            })
            .collect();
        loops.push(Loop { pieces: order.iter().map(|&k| pieces[k]).collect(), junctions });
    }
    Ok(loops)
}

/// A sure bound on how fast a function of the motion's parameter changes: `over(t0, t1)` no
/// less than |f′| anywhere on [t0, t1] (not a sampled estimate: derived from the motion), and
/// `margin` no less than the rounding in a value of f. With it `largest` passes over spans
/// where f provably stays below zero.
pub(super) struct Rate<'a> {
    pub over: &'a dyn Fn(f64, f64) -> f64,
    pub margin: f64,
}

/// The largest of `f` over `[a, b]` (or round a period from `a` if `periodic`): `f` sampled
/// at spacing `step` (within which it changes by no more than `slope` a unit), then each
/// sample's neighbourhood that could hold more than the best found searched by golden
/// section. Stops early once more than `enough` is found.
///
/// With a sure `rate`, the samples whose neighbourhoods `f` provably keeps below zero are
/// passed over (found by halving the samples' range, each half bounded from the value at its
/// middle sample): every value the search would take there is negative, so passing them over
/// changes neither whether the result is negative nor whether it is more than `enough` (which
/// is all a depth is asked); a negative result is then only a bound.
#[allow(clippy::too_many_arguments)]
pub(super) fn largest(f: &dyn Fn(f64) -> f64, a: f64, b: f64, step: f64, slope: f64, periodic: bool, enough: f64, rate: Option<&Rate>) -> f64 {
    let n = (((b - a) / step).ceil() as usize).max(4);
    let h = (b - a) / n as f64;
    let count = if periodic { n } else { n + 1 };
    let time = |i: usize| a + h * i as f64;
    // (NaN where not taken.)
    let mut values = vec![f64::NAN; count];
    let mut live = vec![true; count];
    let mut best = f64::NEG_INFINITY;
    if let Some(rate) = rate {
        live.fill(false);
        let mut ceiling = f64::NEG_INFINITY;
        // The span of what the search takes for samples lo..=hi: their neighbourhoods.
        let span = |lo: usize, hi: usize| {
            (if periodic || lo > 0 { time(lo) - h } else { time(lo) }, if periodic || hi + 1 < count { time(hi) + h } else { time(hi) })
        };
        let mut stack = vec![(0usize, count - 1)];
        while let Some((lo, hi)) = stack.pop() {
            let mid = (lo + hi) / 2;
            if values[mid].is_nan() {
                values[mid] = f(time(mid));
            }
            let v = values[mid];
            best = best.max(v);
            if v > enough {
                // (A sample's: so would the samples' largest be.)
                return v;
            }
            let (t0, t1) = span(lo, hi);
            let reach = (time(mid) - t0).max(t1 - time(mid));
            let cap = v + (rate.over)(t0, t1) * reach + rate.margin;
            if cap < 0.0 {
                ceiling = ceiling.max(cap);
                continue;
            }
            if lo == hi {
                live[mid] = true;
                continue;
            }
            stack.push((mid + 1, hi));
            stack.push((lo, mid));
        }
        if !live.iter().any(|&l| l) {
            return ceiling;
        }
        // The live samples' neighbours, for telling their local largest.
        for i in 0..count {
            if !live[i] {
                continue;
            }
            let neighbours = if periodic { [Some((i + count - 1) % count), Some((i + 1) % count)] } else { [i.checked_sub(1), Some(i + 1).filter(|&j| j < count)] };
            for j in neighbours.into_iter().flatten().chain([i]) {
                if values[j].is_nan() {
                    values[j] = f(time(j));
                    best = best.max(values[j]);
                }
            }
        }
    } else {
        for i in 0..count {
            values[i] = f(time(i));
            best = best.max(values[i]);
        }
    }
    if best > enough {
        return best;
    }
    let reach = slope * h;
    for i in 0..count {
        if !live[i] {
            continue;
        }
        let v = values[i];
        let (before, after) = if periodic {
            (values[(i + count - 1) % count], values[(i + 1) % count])
        } else {
            (if i > 0 { values[i - 1] } else { f64::NEG_INFINITY }, if i + 1 < count { values[i + 1] } else { f64::NEG_INFINITY })
        };
        if v < before || v < after || v + reach < best {
            continue;
        }
        let s = time(i);
        let (mut lo, mut hi) = (if periodic || i > 0 { s - h } else { s }, if periodic || i + 1 < count { s + h } else { s });
        let g = 0.5 * (5f64.sqrt() - 1.0);
        let (mut x1, mut x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        let (mut f1, mut f2) = (f(x1), f(x2));
        for _ in 0..80 {
            if f1 < f2 {
                lo = x1;
                x1 = x2;
                f1 = f2;
                x2 = lo + g * (hi - lo);
                f2 = f(x2);
            } else {
                hi = x2;
                x2 = x1;
                f2 = f1;
                x1 = hi - g * (hi - lo);
                f1 = f(x1);
            }
            if hi - lo <= 1e-15 * (1.0 + s.abs()) {
                break;
            }
        }
        best = best.max(f1).max(f2);
        if best > enough {
            return best;
        }
    }
    best
}
