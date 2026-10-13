//! Knot operations on B-splines (inserting and removing knots, cutting out a segment,
//! raising the degree), shared by curves and surfaces.
//!
//! They work on a [`Spline`]: a non-periodic B-spline whose "poles" are rows of numbers of
//! any length, homogeneous (`w·x, w·y, w·z, w`) so that rational splines need nothing more,
//! and with a surface's whole row or column of poles as one row. A periodic spline is
//! [unrolled](unroll) into one first.

/// A non-periodic B-spline: degree, flat knots (each repeated by its multiplicity) and
/// `flat.len() - degree - 1` poles.
#[derive(Debug, Clone)]
pub(crate) struct Spline {
    pub degree: usize,
    pub flat: Vec<f64>,
    pub poles: Vec<Vec<f64>>,
}

/// `a·x + b·y`, element by element.
fn mix(a: f64, x: &[f64], b: f64, y: &[f64]) -> Vec<f64> {
    x.iter().zip(y).map(|(x, y)| a * x + b * y).collect()
}

/// The largest distance between two rows, by blocks of 4 (one homogeneous point each).
fn distance(x: &[f64], y: &[f64]) -> f64 {
    x.chunks(4)
        .zip(y.chunks(4))
        .map(|(a, b)| a.iter().zip(b).map(|(a, b)| (a - b) * (a - b)).sum::<f64>().sqrt())
        .fold(0.0, f64::max)
}

/// The distinct knots of flat knots and their multiplicities.
pub(crate) fn group(flat: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut multiplicities: Vec<usize> = Vec::new();
    for &k in flat {
        match knots.last() {
            Some(&last) if last == k => *multiplicities.last_mut().expect("a multiplicity") += 1,
            _ => {
                knots.push(k);
                multiplicities.push(1);
            }
        }
    }
    (knots, multiplicities)
}

/// The flat knots of distinct knots with multiplicities.
pub(crate) fn flatten(knots: &[f64], multiplicities: &[usize]) -> Vec<f64> {
    knots.iter().zip(multiplicities).flat_map(|(&k, &m)| std::iter::repeat(k).take(m)).collect()
}

impl Spline {
    /// The parameter range: from flat knot `degree` to flat knot `poles`.
    pub fn domain(&self) -> (f64, f64) {
        (self.flat[self.degree], self.flat[self.poles.len()])
    }

    /// How many flat knots equal `u`.
    pub fn multiplicity(&self, u: f64) -> usize {
        self.flat.iter().filter(|&&k| k == u).count()
    }

    /// Inserts `u` (within the domain) `times` times more, without changing the shape, at
    /// most up to multiplicity `degree + 1` (where the spline splits in two).
    pub fn insert(&mut self, u: f64, times: usize) {
        let p = self.degree;
        let times = times.min((p + 1).saturating_sub(self.multiplicity(u)));
        for _ in 0..times {
            let n = self.poles.len();
            // The span: flat[k] <= u < flat[k + 1], p <= k < n (u at the domain's end in the
            // last span).
            let k = self.flat[..n].partition_point(|&t| t <= u).saturating_sub(1).clamp(p, n - 1);
            let mut poles = Vec::with_capacity(n + 1);
            poles.extend_from_slice(&self.poles[..=k - p]);
            for i in k - p + 1..=k {
                let span = self.flat[i + p] - self.flat[i];
                let alpha = if span > 0.0 { (u - self.flat[i]) / span } else { 0.0 };
                poles.push(mix(alpha, &self.poles[i], 1.0 - alpha, &self.poles[i - 1]));
            }
            poles.extend_from_slice(&self.poles[k..]);
            self.poles = poles;
            self.flat.insert(k + 1, u);
        }
    }

    /// Removes one copy of the interior knot `u` if the shape moves by at most `tolerance`
    /// (in homogeneous coordinates); returns whether it did.
    pub fn remove_once(&mut self, u: f64, tolerance: f64) -> bool {
        let p = self.degree;
        let n = self.poles.len();
        let (lo, hi) = self.domain();
        let s = self.multiplicity(u);
        if s == 0 || !(lo < u && u < hi) {
            return false;
        }
        // The last flat index of u; poles first..=last are replaced by last - first new ones.
        let r = self.flat.iter().rposition(|&k| k == u).expect("the knot");
        let (first, last) = (r - p, r - s);
        // P[i] = a(i)·Q[i] + (1 - a(i))·Q[i - 1] for i in first..=last, with Q[first - 1] =
        // P[first - 1] and Q[last] = P[last + 1]: solved from both ends, the middle
        // equation checked.
        let alpha = |i: usize| (u - self.flat[i]) / (self.flat[i + p + 1] - self.flat[i]);
        let count = last - first;
        let from_left = count.div_ceil(2);
        let mut q: Vec<Vec<f64>> = vec![Vec::new(); count];
        let below = |q: &Vec<Vec<f64>>, j: usize| if j == 0 { self.poles[first - 1].clone() } else { q[j - 1].clone() };
        for j in 0..from_left {
            let a = alpha(first + j);
            let previous = below(&q, j);
            q[j] = mix(1.0 / a, &self.poles[first + j], -(1.0 - a) / a, &previous);
        }
        for j in (from_left..count).rev() {
            let a = alpha(first + j + 1);
            let next = if j + 1 == count { self.poles[last + 1].clone() } else { q[j + 1].clone() };
            q[j] = mix(1.0 / (1.0 - a), &self.poles[first + j + 1], -a / (1.0 - a), &next);
        }
        let i = first + from_left;
        let a = alpha(i);
        let upper = if from_left == count { self.poles[last + 1].clone() } else { q[from_left].clone() };
        let lower = below(&q, from_left);
        if distance(&self.poles[i], &mix(a, &upper, 1.0 - a, &lower)) > tolerance {
            return false;
        }
        let mut poles = Vec::with_capacity(n - 1);
        poles.extend_from_slice(&self.poles[..first]);
        poles.extend(q);
        poles.extend_from_slice(&self.poles[last + 1..]);
        self.poles = poles;
        self.flat.remove(r);
        true
    }

    /// Restricts the spline to `[u1, u2]` (within the domain), clamped there: both ends get
    /// multiplicity `degree + 1`.
    pub fn segment(&mut self, u1: f64, u2: f64) {
        let p = self.degree;
        self.insert(u1, p + 1);
        self.insert(u2, p + 1);
        let a = self.flat.iter().position(|&k| k == u1).expect("the first knot");
        let b = self.flat.iter().rposition(|&k| k == u2).expect("the last knot");
        self.poles = self.poles[a..b - p].to_vec();
        self.flat = self.flat[a..=b].to_vec();
    }

    /// The same spline traversed backwards over the same domain.
    pub fn reversed(&self) -> Spline {
        let (lo, hi) = self.domain();
        let flat = self.flat.iter().rev().map(|&k| lo + hi - k).collect();
        let poles = self.poles.iter().rev().cloned().collect();
        Spline { degree: self.degree, flat, poles }
    }

    /// Derivatives `0..=k` (as rows) at the end of the domain. The spline must be clamped.
    pub fn end_derivatives(&self, k: usize) -> Vec<Vec<f64>> {
        let mut out = vec![self.poles[self.poles.len() - 1].clone()];
        let (mut flat, mut poles, mut p) = (self.flat.clone(), self.poles.clone(), self.degree);
        for _ in 0..k {
            if p == 0 {
                out.push(vec![0.0; out[0].len()]);
                continue;
            }
            // The derivative's poles: p (P[i+1] - P[i]) / (u[i+p+1] - u[i+1]).
            poles = (0..poles.len() - 1)
                .map(|i| {
                    let span = flat[i + p + 1] - flat[i + 1];
                    let c = if span > 0.0 { p as f64 / span } else { 0.0 };
                    mix(c, &poles[i + 1], -c, &poles[i])
                })
                .collect();
            flat = flat[1..flat.len() - 1].to_vec();
            p -= 1;
            out.push(poles[poles.len() - 1].clone());
        }
        out
    }

    /// Appends a Bézier piece over `[end, end + delta]` with the end's derivatives `0..=k`
    /// (so the junction is C`k`) and last pole `last`; the poles between go evenly to
    /// `last`. The spline must be clamped, with degree above `k`.
    pub fn append_piece(&mut self, delta: f64, k: usize, last: Vec<f64>) {
        let n = self.degree;
        debug_assert!(n > k);
        let derivatives = self.end_derivatives(k);
        let binomial = |n: usize, k: usize| (0..k).fold(1.0, |b, i| b * (n - i) as f64 / (i + 1) as f64);
        // Bézier derivatives at the start: D_r = n!/(n-r)! / delta^r · Σ (-1)^(r-i) C(r,i) P_i.
        let mut piece: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
        for (r, d) in derivatives.iter().enumerate() {
            let falling = (0..r).fold(1.0, |f, i| f * (n - i) as f64);
            let mut row: Vec<f64> = d.iter().map(|x| x * delta.powi(r as i32) / falling).collect();
            for (i, p) in piece.iter().enumerate() {
                let c = -binomial(r, i) * if (r - i) % 2 == 0 { 1.0 } else { -1.0 };
                row = mix(1.0, &row, c, p);
            }
            piece.push(row);
        }
        let from = piece[k].clone();
        for i in k + 1..=n {
            let f = (i - k) as f64 / (n - k) as f64;
            piece.push(mix(1.0 - f, &from, f, &last));
        }
        let (_, end) = self.domain();
        self.flat.pop();
        self.flat.extend(std::iter::repeat(end + delta).take(n + 1));
        self.poles.extend(piece.into_iter().skip(1));
    }

    /// Prepends a piece as [`Spline::append_piece`] over `[start - delta, start]`, keeping
    /// the parameters of the rest.
    pub fn prepend_piece(&mut self, delta: f64, k: usize, first: Vec<f64>) {
        let mut reversed = self.reversed();
        reversed.append_piece(delta, k, first);
        let mut out = reversed.reversed();
        // Reversing back shifted the old part by `delta`.
        out.flat.iter_mut().for_each(|x| *x -= delta);
        *self = out;
    }

    /// Raises the degree by `t` without changing the shape. The spline must be clamped.
    pub fn elevate(&mut self, t: usize) {
        if t == 0 {
            return;
        }
        let p = self.degree;
        let (knots, multiplicities) = group(&self.flat);
        // Bézier pieces: every interior knot to multiplicity p.
        for (&k, &m) in knots[1..knots.len() - 1].iter().zip(&multiplicities[1..]) {
            self.insert(k, p - m.min(p));
        }
        let pieces = knots.len() - 1;
        let binomial = |n: usize, k: usize| (0..k).fold(1.0, |b, i| b * (n - i) as f64 / (i + 1) as f64);
        let mut poles: Vec<Vec<f64>> = Vec::with_capacity(pieces * (p + t) + 1);
        for piece in 0..pieces {
            let old = &self.poles[piece * p..=piece * p + p];
            let start = if piece == 0 { 0 } else { 1 };
            for i in start..=p + t {
                let mut row = vec![0.0; old[0].len()];
                for j in i.saturating_sub(t)..=i.min(p) {
                    let c = binomial(p, j) * binomial(t, i - j) / binomial(p + t, i);
                    for (r, x) in row.iter_mut().zip(&old[j]) {
                        *r += c * x;
                    }
                }
                poles.push(row);
            }
        }
        let ends = p + t + 1;
        let inner: Vec<usize> = knots[1..knots.len() - 1].iter().map(|_| p + t).collect();
        let multiplicities_now: Vec<usize> = std::iter::once(ends).chain(inner).chain(std::iter::once(ends)).collect();
        self.flat = flatten(&knots, &multiplicities_now);
        self.poles = poles;
        self.degree = p + t;
        // Back to the original continuity: each interior knot to its multiplicity + t (exact
        // removals; one numerically refused only keeps a knot more).
        let scale = self.poles.iter().flatten().fold(1.0f64, |s, x| s.max(x.abs()));
        for (&k, &m) in knots[1..knots.len() - 1].iter().zip(&multiplicities[1..]) {
            for _ in 0..p - m.min(p) {
                if !self.remove_once(k, 1e-9 * scale) {
                    break;
                }
            }
        }
    }
}

/// A periodic B-spline as the knots of one period (the last knot the first one period on,
/// with the same multiplicity) and one period's poles.
#[derive(Debug, Clone)]
pub(crate) struct Periodic {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub multiplicities: Vec<usize>,
    pub poles: Vec<Vec<f64>>,
}

impl Periodic {
    fn period(&self) -> f64 {
        self.knots[self.knots.len() - 1] - self.knots[0]
    }

    /// One period's flat knots (without the last knot).
    fn one(&self) -> Vec<f64> {
        let last = self.knots.len() - 1;
        flatten(&self.knots[..last], &self.multiplicities[..last])
    }

    /// Flat knot at extended position `e`: positions continue the period's flat knots both
    /// ways, laid out so that pole `e mod n` goes with position `e` (as evaluation does).
    fn knot_at(&self, one: &[f64], e: isize) -> f64 {
        let n = one.len() as isize;
        let shift = (self.degree + 1 - self.multiplicities[0]) as isize;
        let j = e - shift;
        one[j.rem_euclid(n) as usize] + j.div_euclid(n) as f64 * self.period()
    }

    /// The same curve as a non-periodic (unclamped) spline whose domain covers `[u1, u2]`.
    pub fn unroll(&self, u1: f64, u2: f64) -> Spline {
        let p = self.degree as isize;
        let one = self.one();
        let n = one.len() as isize;
        // Extended positions e0..=e1 hold the poles; the domain is [knot(e0 + p), knot(e1 + 1)].
        let mut e0 = 0;
        while self.knot_at(&one, e0 + p) > u1 {
            e0 -= n;
        }
        let mut e1 = n + p - 1;
        while self.knot_at(&one, e1 + 1) < u2 {
            e1 += n;
        }
        let flat = (e0..=e1 + p + 1).map(|e| self.knot_at(&one, e)).collect();
        let poles = (e0..=e1).map(|e| self.poles[e.rem_euclid(n) as usize].clone()).collect();
        Spline { degree: self.degree, flat, poles }
    }

    /// Inserts knot `u` once (moved into the first period), without changing the shape.
    pub fn insert_once(&mut self, u: f64) {
        let p = self.degree as isize;
        let (lo, period) = (self.knots[0], self.period());
        let u = lo + (u - lo).rem_euclid(period);
        let one = self.one();
        let n = one.len() as isize;
        // The span of u among extended positions p..n + p.
        let mut k = p;
        while self.knot_at(&one, k + 1) <= u {
            k += 1;
        }
        let pole = |e: isize| &self.poles[e.rem_euclid(n) as usize];
        // New poles at positions k - p + 1..=k - p + 1 + n (all n + 1 of a period).
        let start = k - p + 1;
        let mut placed: Vec<(isize, Vec<f64>)> = Vec::with_capacity(n as usize + 1);
        for e in start..=start + n {
            let row = if e <= k {
                let (a, b) = (self.knot_at(&one, e), self.knot_at(&one, e + p));
                let alpha = (u - a) / (b - a);
                mix(alpha, pole(e), 1.0 - alpha, pole(e - 1))
            } else {
                pole(e - 1).clone()
            };
            placed.push((e, row));
        }
        // Inserting at the start of the period shifts the period's first flat knot one
        // position on.
        let at_start = u == lo;
        let offset = if at_start { 1 } else { 0 };
        let mut poles = vec![Vec::new(); n as usize + 1];
        for (e, row) in placed {
            poles[(e - offset).rem_euclid(n + 1) as usize] = row;
        }
        self.poles = poles;
        match self.knots.iter().position(|&x| x == u) {
            Some(0) => {
                let last = self.multiplicities.len() - 1;
                self.multiplicities[0] += 1;
                self.multiplicities[last] += 1;
            }
            Some(i) => self.multiplicities[i] += 1,
            None => {
                let i = self.knots.partition_point(|&x| x < u);
                self.knots.insert(i, u);
                self.multiplicities.insert(i, 1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cubic over [0, 3] with one interior knot of each multiplicity 1, its value at `u`
    /// (de Boor, by brute force on the basis functions).
    fn value(s: &Spline, u: f64) -> Vec<f64> {
        fn basis(flat: &[f64], i: usize, p: usize, u: f64, last: bool) -> f64 {
            if p == 0 {
                let inside = flat[i] <= u && u < flat[i + 1];
                let at_end = last && u == flat[i + 1] && flat[i] < flat[i + 1];
                return if inside || at_end { 1.0 } else { 0.0 };
            }
            let mut b = 0.0;
            let d1 = flat[i + p] - flat[i];
            if d1 > 0.0 {
                b += (u - flat[i]) / d1 * basis(flat, i, p - 1, u, last);
            }
            let d2 = flat[i + p + 1] - flat[i + 1];
            if d2 > 0.0 {
                b += (flat[i + p + 1] - u) / d2 * basis(flat, i + 1, p - 1, u, last);
            }
            b
        }
        let (_, hi) = s.domain();
        let last = u == hi;
        // At the domain's end only the last span counts.
        let mut out = vec![0.0; s.poles[0].len()];
        for (i, pole) in s.poles.iter().enumerate() {
            let b = basis(&s.flat, i, s.degree, u, last);
            for (o, x) in out.iter_mut().zip(pole) {
                *o += b * x;
            }
        }
        out
    }

    fn cubic() -> Spline {
        Spline {
            degree: 3,
            flat: vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0],
            poles: [[0.0, 0.0], [1.0, 2.0], [2.0, -1.0], [3.0, 3.0], [4.0, 0.0], [5.0, 1.0]]
                .iter()
                .map(|p| vec![p[0], p[1], 0.0, 1.0])
                .collect(),
        }
    }

    fn same_shape(a: &Spline, b: &Spline, from: f64, to: f64) {
        for i in 0..=20 {
            let u = from + (to - from) * i as f64 / 20.0;
            assert!(distance(&value(a, u), &value(b, u)) < 1e-9, "at {u}: {:?} vs {:?}", value(a, u), value(b, u));
        }
    }

    #[test]
    fn insertion_and_removal_keep_the_shape() {
        let original = cubic();
        let mut s = original.clone();
        s.insert(1.5, 2);
        assert_eq!(s.poles.len(), 8);
        same_shape(&original, &s, 0.0, 3.0);
        assert!(s.remove_once(1.5, 1e-9) && s.remove_once(1.5, 1e-9));
        assert_eq!(s.poles.len(), 6);
        same_shape(&original, &s, 0.0, 3.0);
        // A knot the shape needs stays.
        assert!(!s.remove_once(1.0, 1e-9));
    }

    #[test]
    fn a_segment_keeps_its_part_of_the_shape() {
        let original = cubic();
        let mut s = original.clone();
        s.segment(0.5, 2.25);
        assert_eq!(s.domain(), (0.5, 2.25));
        same_shape(&original, &s, 0.5, 2.25);
    }

    #[test]
    fn elevation_keeps_the_shape_and_continuity() {
        let original = cubic();
        let mut s = original.clone();
        s.elevate(2);
        assert_eq!(s.degree, 5);
        let (_, multiplicities) = group(&s.flat);
        assert_eq!(multiplicities, vec![6, 3, 3, 6]);
        same_shape(&original, &s, 0.0, 3.0);
    }

    #[test]
    fn a_periodic_insertion_keeps_the_shape() {
        // A periodic quadratic: four poles round a square.
        let periodic = Periodic {
            degree: 2,
            knots: vec![0.0, 1.0, 2.0, 3.0, 4.0],
            multiplicities: vec![1, 1, 1, 1, 1],
            poles: [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]].iter().map(|p| vec![p[0], p[1], 0.0, 1.0]).collect(),
        };
        for u in [0.0, 0.5, 2.5, 3.75] {
            let mut inserted = periodic.clone();
            inserted.insert_once(u);
            assert_eq!(inserted.poles.len(), 5);
            same_shape(&periodic.unroll(0.0, 4.0), &inserted.unroll(0.0, 4.0), 0.0, 4.0);
        }
    }
}
