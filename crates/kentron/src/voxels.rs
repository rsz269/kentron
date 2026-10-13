//! Boolean grids of voxels (occupancy, as [`crate::mesh::voxelize`] gives it): erosion,
//! dilation, opening and closing by a structuring element, as image morphology has them,
//! with cells outside the grid empty.
//!
//! A grid of `dims = [nx, ny, nz]` is a slice of `nx * ny * nz` cells, x-major
//! (`(i * ny + j) * nz + k`). A structuring element is a grid of odd size centred on its
//! middle cell. A ball (the offsets within a radius) is applied through the squared
//! distances to the nearest cells of a kind, in time linear in the grid whatever its
//! radius; another element cell by cell.

/// A structuring element: the offsets it holds, from its centre.
pub struct Element {
    offsets: Vec<[isize; 3]>,
    /// For a ball, the squared radius in cells: it holds exactly the offsets `(i, j, k)` with
    /// `i² + j² + k² <= ball`.
    ball: Option<u32>,
}

impl Element {
    /// The element of the cells set in `cells`, a grid of odd sizes `dims` centred on its
    /// middle cell; `None` for even sizes.
    pub fn new(cells: &[bool], dims: [usize; 3]) -> Option<Element> {
        if dims.iter().any(|&n| n % 2 == 0) || cells.len() != dims.iter().product::<usize>() {
            return None;
        }
        let half = dims.map(|n| (n / 2) as isize);
        let mut offsets = Vec::new();
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    if cells[(i * dims[1] + j) * dims[2] + k] {
                        offsets.push([i as isize - half[0], j as isize - half[1], k as isize - half[2]]);
                    }
                }
            }
        }
        // A ball if it holds every offset (of the grid) up to its largest squared distance and
        // none beyond.
        let squared = |o: &[isize; 3]| (o[0] * o[0] + o[1] * o[1] + o[2] * o[2]) as u32;
        let ball = offsets.iter().map(squared).max().filter(|&s| {
            let fits = |n: isize| n * n <= s as isize;
            let reach_in_grid = half.iter().all(|&h| !fits(h + 1));
            let all_in = (0..dims[0]).all(|i| {
                (0..dims[1]).all(|j| {
                    (0..dims[2]).all(|k| {
                        let o = [i as isize - half[0], j as isize - half[1], k as isize - half[2]];
                        (squared(&o) <= s) == cells[(i * dims[1] + j) * dims[2] + k]
                    })
                })
            });
            reach_in_grid && all_in
        });
        Some(Element { offsets, ball })
    }
}

/// The cells whose element, placed on them, lies wholly in `grid` (cells outside it empty).
pub fn erode(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    by_density(grid, dims, element, erode_all, |g, set, d, e| {
        let mut out = vec![false; g.len()];
        for c in erode_listed(g, set, d, e) {
            out[c] = true;
        }
        out
    })
}

/// The cells whose element, placed on them, meets `grid`.
pub fn dilate(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    by_density(grid, dims, element, dilate_all, |g, set, d, e| dilate_listed(set, g.len(), d, e).0)
}

/// Erosion, then dilation: `grid` without the parts the element does not fit in.
#[doc(alias = "binary_opening")]
pub fn open(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    by_density(grid, dims, element, |g, d, e| dilate_all(&erode_all(g, d, e), d, e), |g, set, d, e| dilate_listed(&erode_listed(g, set, d, e), g.len(), d, e).0)
}

/// Dilation, then erosion: `grid` with the gaps the element does not fit in filled.
#[doc(alias = "binary_closing")]
pub fn close(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    by_density(grid, dims, element, |g, d, e| erode_all(&dilate_all(g, d, e), d, e), |g, set, d, e| {
        let (dilated, listed) = dilate_listed(set, g.len(), d, e);
        let mut out = vec![false; g.len()];
        for c in erode_listed(&dilated, &listed, d, e) {
            out[c] = true;
        }
        out
    })
}

/// `operation` over the whole grid, unless the grid has few cells set for the element's
/// size: then `sparse`, from the list of those cells (cell by cell, the same result).
fn by_density(
    grid: &[bool],
    dims: [usize; 3],
    element: &Element,
    operation: impl Fn(&[bool], [usize; 3], &Element) -> Vec<bool>,
    sparse: impl Fn(&[bool], &[usize], [usize; 3], &Element) -> Vec<bool>,
) -> Vec<bool> {
    let count = grid.iter().filter(|&&v| v).count();
    if count.saturating_mul(element.offsets.len()) > grid.len() / 4 {
        return operation(grid, dims, element);
    }
    let set: Vec<usize> = grid.iter().enumerate().filter_map(|(c, &v)| v.then_some(c)).collect();
    sparse(grid, &set, dims, element)
}

/// The cell `c` moved by `offset`, if still in the grid.
fn moved(c: usize, offset: &[isize; 3], dims: [usize; 3]) -> Option<usize> {
    let [_, ny, nz] = dims;
    let at = [c / (ny * nz), (c / nz) % ny, c % nz];
    let mut to = [0usize; 3];
    for a in 0..3 {
        let p = at[a] as isize + offset[a];
        if p < 0 || p >= dims[a] as isize {
            return None;
        }
        to[a] = p as usize;
    }
    Some((to[0] * ny + to[1]) * nz + to[2])
}

/// The cells of `set` (the set cells of `grid`) whose element lies wholly in the grid.
fn erode_listed(grid: &[bool], set: &[usize], dims: [usize; 3], element: &Element) -> Vec<usize> {
    set.iter().copied().filter(|&c| element.offsets.iter().all(|o| moved(c, o, dims).is_some_and(|m| grid[m]))).collect()
}

/// The grid of the cells the element, placed on cells of `set`, covers (reflected, as
/// dilation has it), and the list of them.
fn dilate_listed(set: &[usize], cells: usize, dims: [usize; 3], element: &Element) -> (Vec<bool>, Vec<usize>) {
    let mut out = vec![false; cells];
    let mut listed = Vec::new();
    for &c in set {
        for o in &element.offsets {
            if let Some(m) = moved(c, &[-o[0], -o[1], -o[2]], dims) {
                if !out[m] {
                    out[m] = true;
                    listed.push(m);
                }
            }
        }
    }
    (out, listed)
}

/// [`erode`] over the whole grid.
fn erode_all(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    match element.ball {
        // Farther than the ball from every empty cell, those just outside included.
        Some(s) if s < NEAR => near(grid, false, dims, true, s).into_iter().map(|d| u32::from(d) > s).collect(),
        Some(s) => squared_distances_to(dims, |c| !grid[c], true).into_iter().map(|d| d > s).collect(),
        None => cell_by_cell(grid, dims, element, true),
    }
}

/// [`dilate`] over the whole grid.
fn dilate_all(grid: &[bool], dims: [usize; 3], element: &Element) -> Vec<bool> {
    match element.ball {
        Some(s) if s < NEAR => near(grid, true, dims, false, s).into_iter().map(|d| u32::from(d) <= s).collect(),
        Some(s) => squared_distances_to(dims, |c| grid[c], false).into_iter().map(|d| d <= s).collect(),
        None => cell_by_cell(grid, dims, element, false),
    }
}

/// The grid's connected parts (cells sharing a face, as `scipy.ndimage.label` joins them by
/// default), in the order of their first cells: each one's number of cells and its centre
/// (the mean of its cells' indices).
#[doc(alias = "label")]
pub fn components(grid: &[bool], dims: [usize; 3]) -> Vec<(usize, [f64; 3])> {
    let [nx, ny, nz] = dims;
    let mut seen = vec![false; grid.len()];
    let mut stack: Vec<usize> = Vec::new();
    let mut out = Vec::new();
    for start in 0..grid.len() {
        if !grid[start] || seen[start] {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        let (mut count, mut sums) = (0usize, [0u64; 3]);
        while let Some(c) = stack.pop() {
            let (i, j, k) = (c / (ny * nz), (c / nz) % ny, c % nz);
            count += 1;
            sums[0] += i as u64;
            sums[1] += j as u64;
            sums[2] += k as u64;
            let neighbours = [
                (i > 0).then(|| c - ny * nz),
                (i + 1 < nx).then(|| c + ny * nz),
                (j > 0).then(|| c - nz),
                (j + 1 < ny).then(|| c + nz),
                (k > 0).then(|| c - 1),
                (k + 1 < nz).then(|| c + 1),
            ];
            for n in neighbours.into_iter().flatten() {
                if grid[n] && !seen[n] {
                    seen[n] = true;
                    stack.push(n);
                }
            }
        }
        // (The sums are exact in f64 well beyond any grid's size.)
        out.push((count, sums.map(|s| s as f64 / count as f64)));
    }
    out
}

/// Erosion (`all`) or dilation (any) offset by offset: a cell's result from the cells at its
/// element's offsets (reflected for dilation; outside the grid, empty).
fn cell_by_cell(grid: &[bool], dims: [usize; 3], element: &Element, all: bool) -> Vec<bool> {
    let [nx, ny, nz] = dims.map(|n| n as isize);
    let at = |i: isize, j: isize, k: isize| i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz && grid[((i * ny + j) * nz + k) as usize];
    let mut out = vec![false; grid.len()];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let cell = |o: &[isize; 3]| if all { at(i + o[0], j + o[1], k + o[2]) } else { at(i - o[0], j - o[1], k - o[2]) };
                out[((i * ny + j) * nz + k) as usize] = if all { element.offsets.iter().all(cell) } else { element.offsets.iter().any(cell) };
            }
        }
    }
    out
}

/// No cell of the kind: farther than any distance in a grid.
const FAR: u32 = u32::MAX;

/// Squared radii below this are applied by [`near`] (values in a byte).
const NEAR: u32 = 250;

/// Each cell's squared distance (in cells) to the nearest cell of `grid` equal to `kind`,
/// as [`squared_distances_to`] gives it, but only up to `s` (`s + 1` for any beyond): axis by
/// axis, the least of the values within `⌊√s⌋` cells plus their squared offsets. (Exact as
/// far as it goes: a cell within `s` of its nearest is within that reach of it along each
/// axis, through values within `s`.) Each pass runs along the grid's memory, whole rows or
/// the whole grid at a time.
fn near(grid: &[bool], kind: bool, dims: [usize; 3], outside: bool, s: u32) -> Vec<u8> {
    let beyond = (s + 1) as u8;
    let reach = (s as f64).sqrt().floor() as usize;
    // (Outside the grid: of the kind if `outside`, else none near.)
    let edge = if outside { 0 } else { beyond };
    let [nx, ny, nz] = dims;
    let n = nx * ny * nz;
    let mut a: Vec<u8> = grid.iter().map(|&g| if g == kind { 0 } else { beyond }).collect();
    let squares: Vec<u8> = (0..=reach).map(|r| (r * r) as u8).collect();
    let add = |v: u8, r: usize| v.saturating_add(squares[r]).min(beyond);
    // Along z: the whole grid shifted by each offset (a cell's z neighbours are the next
    // and previous in memory), then the cells near the ends of rows, whose neighbours
    // there are in other rows, done again within their own rows.
    let mut b = a.clone();
    for r in 1..=reach.min(n) {
        let square = squares[r];
        for (t, &f) in b[r..].iter_mut().zip(&a[..n - r]) {
            *t = (*t).min(f.saturating_add(square).min(beyond));
        }
        for (t, &f) in b[..n - r].iter_mut().zip(&a[r..]) {
            *t = (*t).min(f.saturating_add(square).min(beyond));
        }
    }
    let ends = reach.min(nz);
    for row in 0..nx * ny {
        let base = row * nz;
        for k in (0..ends).chain(nz.saturating_sub(ends).max(ends)..nz) {
            let mut m = a[base + k];
            for r in 1..=reach {
                let lo = if k >= r { a[base + k - r] } else { edge };
                let hi = if k + r < nz { a[base + k + r] } else { edge };
                m = m.min(add(lo, r)).min(add(hi, r));
            }
            b[base + k] = m;
        }
    }
    // Along y, then x: rows combined whole (the inner loop runs along z).
    let pass = |from: &[u8], to: &mut [u8], count: usize, step: usize| {
        let outer = from.len() / (count * step);
        for o in 0..outer {
            for q in 0..count {
                let at = |p: usize| (o * count + p) * step;
                let target = &mut to[at(q)..at(q) + step];
                target.copy_from_slice(&from[at(q)..at(q) + step]);
                for r in 1..=reach {
                    for side in [q.checked_sub(r), (q + r < count).then_some(q + r)] {
                        match side {
                            Some(p) => {
                                let source = &from[at(p)..at(p) + step];
                                for (t, &v) in target.iter_mut().zip(source) {
                                    *t = (*t).min(add(v, r));
                                }
                            }
                            None => {
                                let v = add(edge, r);
                                for t in target.iter_mut() {
                                    *t = (*t).min(v);
                                }
                            }
                        }
                    }
                }
            }
        }
    };
    pass(&b, &mut a, ny, nz);
    pass(&a, &mut b, nx, ny * nz);
    b
}

/// Each cell's squared distance (in cells) to the nearest cell of a kind (`of`, by index),
/// with the cells just outside the grid of that kind too if `outside`: the exact Euclidean
/// distance transform, axis by axis (Felzenszwalb and Huttenlocher's lower envelopes of
/// parabolas).
fn squared_distances_to(dims: [usize; 3], of: impl Fn(usize) -> bool, outside: bool) -> Vec<u32> {
    // With the outside: the grid inside a layer of such cells.
    let pad = usize::from(outside);
    let [px, py, pz] = dims.map(|n| n + 2 * pad);
    let inside = |p: usize, n: usize| p >= pad && p < n + pad;
    let mut d = vec![FAR; px * py * pz];
    for i in 0..px {
        for j in 0..py {
            for k in 0..pz {
                let kind = if inside(i, dims[0]) && inside(j, dims[1]) && inside(k, dims[2]) {
                    of(((i - pad) * dims[1] + (j - pad)) * dims[2] + (k - pad))
                } else {
                    true
                };
                if kind {
                    d[(i * py + j) * pz + k] = 0;
                }
            }
        }
    }
    // Along z, y, then x: each line's values replaced by the lower envelope of the parabolas
    // rooted at them.
    let mut line = Vec::new();
    let mut out = Vec::new();
    let strides = [(py * pz, px), (pz, py), (1, pz)];
    for axis in [2, 1, 0] {
        let (stride, n) = strides[axis];
        let starts: Vec<usize> = (0..px * py * pz).filter(|&c| (c / stride) % n == 0).collect();
        for start in starts {
            line.clear();
            line.extend((0..n).map(|q| d[start + q * stride]));
            envelope(&line, &mut out);
            for (q, &v) in out.iter().enumerate() {
                d[start + q * stride] = v;
            }
        }
    }
    if !outside {
        return d;
    }
    let mut cropped = Vec::with_capacity(dims.iter().product());
    for i in 0..dims[0] {
        for j in 0..dims[1] {
            let from = ((i + 1) * py + (j + 1)) * pz + 1;
            cropped.extend_from_slice(&d[from..from + dims[2]]);
        }
    }
    cropped
}

/// `out[q] = min over p of (q - p)² + f[p]` (with `FAR` for none).
fn envelope(f: &[u32], out: &mut Vec<u32>) {
    out.clear();
    let sites: Vec<usize> = (0..f.len()).filter(|&p| f[p] != FAR).collect();
    if sites.is_empty() {
        out.resize(f.len(), FAR);
        return;
    }
    // The envelope's parabolas (their roots) and where each starts to be the lowest.
    let height = |p: usize| f[p] as f64 + (p * p) as f64;
    let mut roots: Vec<usize> = vec![sites[0]];
    let mut from: Vec<f64> = vec![f64::NEG_INFINITY];
    for &q in &sites[1..] {
        // Parabolas the new one is below from where they start on are off the envelope (the
        // first, starting at -inf, never is).
        loop {
            let v = *roots.last().expect("a parabola");
            let s = (height(q) - height(v)) / (2.0 * (q as f64 - v as f64));
            if s <= *from.last().expect("a start") {
                roots.pop();
                from.pop();
                continue;
            }
            roots.push(q);
            from.push(s);
            break;
        }
    }
    let mut k = 0;
    for q in 0..f.len() {
        while k + 1 < roots.len() && from[k + 1] < q as f64 {
            k += 1;
        }
        let p = roots[k];
        let dq = q.abs_diff(p) as u64;
        out.push((dq * dq + f[p] as u64).min(FAR as u64 - 1) as u32);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ball of squared radius `s` in a grid of `2 n + 1` cells a side.
    fn ball(n: usize, s: isize) -> (Vec<bool>, [usize; 3]) {
        let m = 2 * n + 1;
        let h = n as isize;
        let cells = (0..m * m * m)
            .map(|c| {
                let (i, j, k) = ((c / (m * m)) as isize - h, ((c / m) % m) as isize - h, (c % m) as isize - h);
                i * i + j * j + k * k <= s
            })
            .collect();
        (cells, [m, m, m])
    }

    #[test]
    fn balls_by_distances_match_cell_by_cell() {
        // A lumpy grid: a few overlapping blocks and holes.
        let dims = [23, 19, 17];
        let grid: Vec<bool> = (0..dims.iter().product::<usize>())
            .map(|c| {
                let (i, j, k) = (c / (dims[1] * dims[2]), (c / dims[2]) % dims[1], c % dims[2]);
                let block = (2..15).contains(&i) && (1..12).contains(&j) && (3..14).contains(&k);
                let other = (9..22).contains(&i) && (6..18).contains(&j) && (0..9).contains(&k);
                let hole = (5..7).contains(&i) && (4..9).contains(&j) && (5..11).contains(&k);
                (block || other) && !hole || (i * 7 + j * 3 + k) % 31 == 0
            })
            .collect();
        for (n, s) in [(1, 1), (2, 4), (3, 6), (3, 9)] {
            let (cells, edims) = ball(n, s);
            let element = Element::new(&cells, edims).expect("odd sizes");
            assert!(element.ball.is_some(), "a ball of {s}");
            let plain = Element { offsets: element.offsets.clone(), ball: None };
            assert_eq!(erode(&grid, dims, &element), erode(&grid, dims, &plain), "erosion by {s}");
            assert_eq!(dilate(&grid, dims, &element), dilate(&grid, dims, &plain), "dilation by {s}");
            assert_eq!(open(&grid, dims, &element), open(&grid, dims, &plain));
            assert_eq!(close(&grid, dims, &element), close(&grid, dims, &plain));
            // The general transform (for large balls) agrees too.
            let s = s as u32;
            let by_transform: Vec<bool> = squared_distances_to(dims, |c| !grid[c], true).into_iter().map(|d| d > s).collect();
            assert_eq!(by_transform, erode(&grid, dims, &plain), "erosion by the transform, {s}");
            let by_transform: Vec<bool> = squared_distances_to(dims, |c| grid[c], false).into_iter().map(|d| d <= s).collect();
            assert_eq!(by_transform, dilate(&grid, dims, &plain), "dilation by the transform, {s}");
        }
    }

    #[test]
    fn sparse_grids_by_their_cells_match_whole() {
        // Blobs well inside and at an edge of a larger grid, a ball and a cross.
        let dims = [30, 26, 22];
        let blob = |i: usize, j: usize, k: usize| {
            let inner = (10..15).contains(&i) && (8..13).contains(&j) && (6..12).contains(&k) && (i + j + k) % 7 != 0;
            let edge = i < 3 && (2..6).contains(&j) && k + 2 >= 22;
            inner || edge
        };
        let grid: Vec<bool> = (0..dims.iter().product::<usize>()).map(|c| blob(c / (26 * 22), (c / 22) % 26, c % 22)).collect();
        let (ball_cells, ball_dims) = ball(2, 5);
        let mut cross = vec![false; 27];
        for c in [4, 10, 12, 13, 14, 16, 22] {
            cross[c] = true;
        }
        for (cells, edims) in [(ball_cells, ball_dims), (cross, [3, 3, 3])] {
            let element = Element::new(&cells, edims).expect("odd sizes");
            let plain = Element { offsets: element.offsets.clone(), ball: None };
            assert_eq!(erode(&grid, dims, &element), erode_all(&grid, dims, &plain));
            assert_eq!(dilate(&grid, dims, &element), dilate_all(&grid, dims, &plain));
            assert_eq!(open(&grid, dims, &element), dilate_all(&erode_all(&grid, dims, &plain), dims, &plain));
            assert_eq!(close(&grid, dims, &element), erode_all(&dilate_all(&grid, dims, &plain), dims, &plain));
        }
        let empty = vec![false; grid.len()];
        let (cells, edims) = ball(1, 1);
        assert_eq!(close(&empty, dims, &Element::new(&cells, edims).expect("odd")), empty);
    }

    #[test]
    fn components_of_two_blocks_and_a_cell() {
        let dims = [6, 5, 4];
        let mut grid = vec![false; 120];
        let at = |i: usize, j: usize, k: usize| (i * 5 + j) * 4 + k;
        for (i, j, k) in [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)] {
            grid[at(i, j, k)] = true;
        }
        grid[at(5, 4, 3)] = true;
        // Touching only along an edge: apart.
        grid[at(2, 2, 1)] = true;
        let found = components(&grid, dims);
        assert_eq!(found.len(), 3);
        assert_eq!(found[0], (4, [0.25, 0.5, 0.75]));
        assert_eq!(found[1], (1, [2.0, 2.0, 1.0]));
        assert_eq!(found[2], (1, [5.0, 4.0, 3.0]));
    }

    #[test]
    fn other_elements_are_not_balls() {
        let mut cells = vec![false; 27];
        cells[13] = true;
        cells[14] = true;
        assert!(Element::new(&cells, [3, 3, 3]).expect("odd").ball.is_none());
    }
}
