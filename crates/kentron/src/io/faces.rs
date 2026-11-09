//! Faces built from loops of edges with curves on the face's surface, as the file readers
//! (STEP, IGES) make them: pcurves checked against the edges or projected, put in the right
//! periods and on the right sides of seams, edges added at poles, outer boundaries turned
//! counter-clockwise; tolerances raised to cover the gaps the files leave.

use crate::hashing::FastMap as HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::{BSplineCurve2, Curve2};
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2};
use crate::topology::{EdgeId, Model, Orientation, Shape, VertexId};

fn invalid(what: String) -> Error {
    Error::InvalidInput(what)
}

/// An edge's curve over its range.
pub(crate) struct CurveSpan<'a> {
    pub curve: &'a Curve,
    pub first: f64,
    pub last: f64,
}

/// The tolerances edges and vertices need to cover the gaps between their curves, pcurves
/// and points, as noted face by face; applied at the end.
#[derive(Default)]
pub(crate) struct Gaps {
    edge: HashMap<EdgeId, f64>,
    vertex: HashMap<VertexId, f64>,
    /// Edges split while reading (for a missing seam), each into two halves at a parameter:
    /// faces read after use the halves.
    splits: HashMap<EdgeId, (EdgeId, EdgeId, f64)>,
    /// The faces read, each with its sense against its surface.
    sense: HashMap<crate::topology::FaceId, bool>,
}

impl Gaps {
    /// The loop with each edge split while reading replaced by its halves.
    fn expand(&self, coedges: Vec<Coedge>) -> Vec<Coedge> {
        if self.splits.is_empty() {
            return coedges;
        }
        let mut out = Vec::with_capacity(coedges.len());
        let mut stack: Vec<Coedge> = coedges.into_iter().rev().collect();
        while let Some(c) = stack.pop() {
            match self.splits.get(&c.0) {
                Some(&(h1, h2, t)) => {
                    let (first, second) = ((h1, c.1, c.2.clone(), c.3, t), (h2, c.1, c.2, t, c.4));
                    // (Pushed in reverse: the half the coedge runs along first on top.)
                    if c.1 {
                        stack.push(second);
                        stack.push(first);
                    } else {
                        stack.push(first);
                        stack.push(second);
                    }
                }
                None => out.push(c),
            }
        }
        out
    }

    /// Splits `edge` at `parameter` (by [`Model::split_edge`], in the faces read so far),
    /// noting it for the faces read after.
    fn split(&mut self, model: &mut Model, edge: EdgeId, parameter: f64) -> Result<(EdgeId, EdgeId)> {
        let (curve, _, _) = model.edge_curve(Shape::from(edge)).ok_or_else(|| invalid(format!("{edge:?} has no curve")))?;
        let point = curve.point(parameter);
        let vertex = model.add_vertex(point, crate::tolerance::LINEAR);
        let (h1, h2) = model.split_edge(edge, vertex, parameter)?;
        if let Some(tol) = self.edge.remove(&edge) {
            self.edge.insert(h1, tol);
            self.edge.insert(h2, tol);
        }
        self.splits.insert(edge, (h1, h2, parameter));
        Ok((h1, h2))
    }

    /// Widens the tolerances of an edge and its vertices to cover how far its pcurve on a
    /// surface and its vertices are from its curve.
    fn note(&mut self, model: &Model, edge: EdgeId, pcurve: &Curve2, first: f64, last: f64, surface: &Surface) {
        let Some((curve, _, _)) = model.edge_curve(Shape::from(edge)) else { return };
        let mut deviation: f64 = 0.0;
        for i in 0..=8 {
            let t = first + (last - first) * i as f64 / 8.0;
            let uv = pcurve.point(t);
            deviation = deviation.max(surface.point(uv.x, uv.y).distance(curve.point(t)));
        }
        let (a, b) = model.edge_vertices(Shape::from(edge));
        let ends = [(a, first), (b, last)];
        for (v, t) in ends {
            let id = VertexId::try_from(v.node).expect("a vertex");
            let p = model.point(v);
            let uv = pcurve.point(t);
            let gap = p.distance(curve.point(t)).max(p.distance(surface.point(uv.x, uv.y)));
            let tol = self.vertex.entry(id).or_insert(crate::tolerance::LINEAR);
            *tol = tol.max(gap * 1.1).max(deviation * 1.1);
            deviation = deviation.max(gap);
        }
        let tol = self.edge.entry(edge).or_insert(crate::tolerance::LINEAR);
        *tol = tol.max(deviation * 1.1);
    }


    /// Sets the tolerances noted while reading.
    pub(crate) fn apply(&self, model: &mut Model) -> Result<()> {
        for (&edge, &tol) in &self.edge {
            if tol > crate::tolerance::LINEAR {
                model.set_edge_tolerance(edge, tol)?;
            }
        }
        for (&v, &tol) in &self.vertex {
            if tol > crate::tolerance::LINEAR {
                let p = model.point(Shape::from(v));
                model.set_vertex_point(v.into(), p, tol)?;
            }
        }
        Ok(())
    }
}


/// The face from its loops (outer?, coedges with pcurves): seams and poles placed,
/// outer boundaries counter-clockwise, holes in the outer one's period.
pub(crate) fn finish_face(model: &mut Model, gaps: &mut Gaps, face: Shape, surface: &Surface, loops: Vec<(bool, Vec<(EdgeId, bool, Curve2, f64, f64)>)>, same_sense: bool) -> Result<Shape> {
    let surface = surface.clone();
    // (Seams placed by what all the face's loops span: a hole may tell which side.)
    let spans = face_spans(&loops.iter().map(|(_, c)| c.as_slice()).collect::<Vec<_>>(), &surface);
    let mut sided = Vec::with_capacity(loops.len());
    let mut ties = Vec::new();
    for (k, (outer, coedges)) in loops.into_iter().enumerate() {
        let mut coedges = gaps.expand(coedges);
        ties.extend(seams_to_sides(&mut coedges, &surface, spans).into_iter().map(|(i, shift)| (k, i, shift)));
        ties.extend(period_ties(&coedges, &surface).into_iter().map(|(i, shift)| (k, i, shift)));
        sided.push((outer, coedges));
    }
    // A pcurve along the seam with nothing in the face to tell its side (a half sphere
    // between two meridians): the side the faces already read across its edges tell.
    if !ties.is_empty() {
        let mut other = sided.clone();
        for &(k, i, shift) in &ties {
            use crate::math::Transformable2;
            other[k].1[i].2 = other[k].1[i].2.transformed(&crate::math::Transform2::translation(shift));
        }
        if agreement(model, gaps, &other, &surface, same_sense) > agreement(model, gaps, &sided, &surface, same_sense) {
            sided = other;
        }
    }
    let mut loops = Vec::with_capacity(sided.len());
    for (outer, coedges) in sided {
        loops.push((outer, close_at_poles(model, coedges, &surface)?));
    }
    // Outer boundaries run counter-clockwise in the surface's parameters.
    let (outer, open, clockwise) = orientation(&loops, &surface, same_sense);
    if clockwise {
        for (_, coedges) in &mut loops {
            coedges.reverse();
            for c in coedges.iter_mut() {
                c.1 = !c.1;
            }
        }
    }
    let joined = match open {
        Some(axis) => add_missing_seam(model, gaps, &mut loops, &surface, axis)?,
        None => None,
    };
    // Holes in the outer boundary's period of a periodic surface.
    if let (Some(i), None) = (outer, open) {
        let outer_box = uv_box(&loops[i].1);
        for (j, (_, coedges)) in loops.iter_mut().enumerate() {
            if j != i {
                into_period(coedges, outer_box, &surface);
            }
        }
    }
    // The outer boundary first (the model's outer wire; files may list a hole first).
    if let Some(i) = if open.is_some() { joined } else { outer } {
        let first = loops.remove(i);
        loops.insert(0, first);
    }
    for (_, coedges) in loops {
        let mut ids = Vec::with_capacity(coedges.len());
        for (edge, forward, pcurve, first, last) in coedges {
            gaps.note(model, edge, &pcurve, first, last, &surface);
            let c = model.add_coedge(edge, if forward { Orientation::Forward } else { Orientation::Reversed })?;
            model.set_pcurve_with_range(c, pcurve, first, last)?;
            ids.push(c);
        }
        let wire = model.add_wire(&ids)?;
        model.add(face, wire.into())?;
    }
    if let Ok(id) = crate::topology::FaceId::try_from(face.node) {
        gaps.sense.insert(id, same_sense);
    }
    Ok(if same_sense { face } else { face.reversed() })
}


/// On a periodic surface, the pcurves along an iso-line of the period's parameter (a
/// sphere's meridian) as near their neighbours a period away (across a pole): with the
/// shift to the other place.
fn period_ties(coedges: &[Coedge], surface: &Surface) -> Vec<(usize, crate::math::Vector2)> {
    let n = coedges.len();
    let mut ties = Vec::new();
    if n < 2 {
        return ties;
    }
    for (axis, period) in [(0, surface.u_period()), (1, surface.v_period())] {
        let Some(period) = period else { continue };
        let at = 1e-9 * period;
        for i in 0..n {
            let c = &coedges[i];
            if coedges.iter().filter(|d| d.0 == c.0).count() > 1 {
                continue; // a seam
            }
            let xs: Vec<f64> = (0..=4).map(|k| component(c.2.point(c.3 + (c.4 - c.3) * k as f64 / 4.0), axis)).collect();
            if xs.iter().any(|x| (x - xs[0]).abs() > at) {
                continue;
            }
            let (start, end) = coedge_ends(c);
            let (before, after) = (coedge_ends(&coedges[(i + n - 1) % n]).1, coedge_ends(&coedges[(i + 1) % n]).0);
            let gap = |shift: f64| {
                let moved = |p: Point2| if axis == 0 { Point2::new(p.x + shift, p.y) } else { Point2::new(p.x, p.y + shift) };
                moved(start).distance(before) + moved(end).distance(after)
            };
            for shift in [-period, period] {
                if (gap(shift) - gap(0.0)).abs() <= 1e-9 * (1.0 + gap(0.0)) {
                    ties.push((i, if axis == 0 { crate::math::Vector2::new(shift, 0.0) } else { crate::math::Vector2::new(0.0, shift) }));
                }
            }
        }
    }
    ties
}


/// The face's outer loop, the parameter whose period loops go round open (if one), and
/// whether the loops run clockwise in the surface's parameters (to be turned). Loops round
/// a period (a band of a cylinder between two circles, its seam left out) are open in the
/// parameters: then the face's area from all its loops at once, by an integral seams add
/// nothing to.
fn orientation(loops: &[(bool, Vec<Coedge>)], surface: &Surface, same_sense: bool) -> (Option<usize>, Option<usize>, bool) {
    let areas: Vec<f64> = loops.iter().map(|(_, c)| uv_area(c)).collect();
    let open = round_period(loops, surface);
    let outer = loops.iter().position(|(o, _)| *o).or_else(|| (0..areas.len()).max_by(|&i, &j| areas[i].abs().total_cmp(&areas[j].abs())));
    let clockwise = outer.is_some_and(|i| {
        let scale = areas.iter().map(|a| a.abs()).fold(0.0, f64::max);
        match open {
            // (Loops going round the period both ways, as a band's; a cap's one loop, closed
            // by a pole, says nothing of its side: the file's word for it then.)
            Some(axis) if paired(loops, surface, axis) => loops.iter().map(|(_, c)| band_area(c, axis)).sum::<f64>() < 0.0,
            Some(_) => !same_sense,
            None => areas[i] < -1e-9 * scale.max(1e-300),
        }
    });
    (outer, open, clockwise)
}


/// How far the face's loops, once turned counter-clockwise, run their edges against the
/// faces already read on them (as a closed shell's faces do): each edge of another face run
/// the other way counts one, the same way minus one.
fn agreement(model: &Model, gaps: &Gaps, loops: &[(bool, Vec<Coedge>)], surface: &Surface, same_sense: bool) -> i64 {
    let (_, _, clockwise) = orientation(loops, surface, same_sense);
    let mut score = 0;
    for (_, coedges) in loops {
        for c in coedges {
            // Its direction in space: along the edge's curve or against it.
            let along = (c.1 != clockwise) == same_sense;
            for f in model.ancestors_of::<crate::topology::entities::Face>(c.0) {
                let Some(&sense) = gaps.sense.get(&f) else { continue };
                let uses: Vec<bool> = model
                    .face_coedges(f)
                    .into_iter()
                    .filter_map(|u| model.get(u).filter(|u| u.edge == c.0).map(|u| (u.orientation == Orientation::Forward) == sense))
                    .collect();
                if let [theirs] = uses[..] {
                    score += if theirs != along { 1 } else { -1 };
                }
            }
        }
    }
    score
}


/// The loop with degenerate edges added where consecutive pcurves end at the same point
/// of a singular part of the surface (a cone's apex, a sphere's pole) but apart in its
/// parameters: files leave those edges out.
fn close_at_poles(model: &mut Model, coedges: Vec<(EdgeId, bool, Curve2, f64, f64)>, surface: &Surface) -> Result<Vec<(EdgeId, bool, Curve2, f64, f64)>> {
    let n = coedges.len();
    let ends = |c: &(EdgeId, bool, Curve2, f64, f64)| {
        let (a, b) = (c.2.point(c.3), c.2.point(c.4));
        if c.1 { (a, b) } else { (b, a) }
    };
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..n {
        let (_, end) = ends(&coedges[i]);
        let (start, _) = ends(&coedges[(i + 1) % n]);
        out.push(coedges[i].clone());
        if end.distance(start) <= 1e-9 {
            continue;
        }
        let (p, q) = (surface.point(end.x, end.y), surface.point(start.x, start.y));
        let (_, du, dv) = surface.d1(end.x, end.y);
        let singular = du.length().min(dv.length()) <= 1e-9 * du.length().max(dv.length()).max(1e-300);
        if !singular || p.distance(q) > 1e-6 * (1.0 + p.distance(Point::ORIGIN)) {
            continue;
        }
        // At the vertex the previous coedge ends at.
        let (a, b) = model.edge_vertices(Shape::from(coedges[i].0));
        let vertex = if coedges[i].1 { b } else { a };
        let v = VertexId::try_from(vertex.node).expect("a vertex");
        let edge = model.add_edge(v, v, crate::tolerance::LINEAR)?;
        let pcurve = Curve2::BSpline(BSplineCurve2::new(vec![end, start], vec![0.0, 1.0], vec![2, 2], 1));
        out.push((edge, true, pcurve, 0.0, 1.0));
    }
    Ok(out)
}


/// True if the pcurve follows the edge's curve on the surface.
pub(crate) fn fits(pcurve: &Curve2, edge: &CurveSpan, surface: &Surface) -> bool {
    let length: f64 = (0..8)
        .map(|i| {
            let t = |k: f64| edge.first + (edge.last - edge.first) * k / 8.0;
            edge.curve.point(t(i as f64)).distance(edge.curve.point(t(i as f64 + 1.0)))
        })
        .sum();
    let allowed = (2e-3 * length).clamp(1e-4, 1e-2);
    (0..=4).all(|i| {
        let t = edge.first + (edge.last - edge.first) * i as f64 / 4.0;
        let uv = pcurve.point(t);
        let p = surface.point(uv.x, uv.y);
        let q = edge.curve.point(t);
        // Files are often accurate to only 1e-3 or so; a pcurve parametrized otherwise than
        // the curve (a rational circle) or on the other side of a seam is off by more.
        p.distance(q) <= allowed
    })
}


/// True if the pcurve, mapped onto the surface, runs along the edge's curve (whatever its
/// parametrization): its points near the curve, its ends at the edge's.
pub(crate) fn follows(pcurve: &Curve2, edge: &CurveSpan, surface: &Surface) -> bool {
    let length: f64 = (0..8)
        .map(|i| {
            let t = |k: f64| edge.first + (edge.last - edge.first) * k / 8.0;
            edge.curve.point(t(i as f64)).distance(edge.curve.point(t(i as f64 + 1.0)))
        })
        .sum();
    let allowed = (2e-3 * length).clamp(1e-4, 1e-2);
    let on_surface = |t: f64| {
        let uv = pcurve.point(t);
        surface.point(uv.x, uv.y)
    };
    let ends_match = on_surface(edge.first).distance(edge.curve.point(edge.first)) <= allowed
        && on_surface(edge.last).distance(edge.curve.point(edge.last)) <= allowed;
    ends_match
        && (1..8).all(|i| {
            let p = on_surface(edge.first + (edge.last - edge.first) * i as f64 / 8.0);
            crate::geometry::extrema::nearest_on_curve(p, edge.curve, edge.first, edge.last).distance <= allowed
        })
}


/// The edge's curve projected onto the surface (`has_pcurve`: the file gives a pcurve to fall
/// back on).
pub(crate) fn project(edge: &CurveSpan, surface: &Surface, has_pcurve: bool) -> Result<Curve2> {
    use crate::geometry::project as pr;
    // (An edge running its curve backwards: the pcurve over the range the other way round, at
    // the same parameters.)
    if edge.first > edge.last {
        return project(&CurveSpan { curve: edge.curve, first: edge.last, last: edge.first }, surface, has_pcurve);
    }
    if let Some(pc) = pr::analytic_pcurve(edge.curve, surface) {
        return Ok(pc);
    }
    // An edge well off the surface (bad data) fails at once, before costly fitting. (Each
    // point's foot from the one before, the whole surface searched only to confirm one off.)
    // (The first point's foot, from the whole surface, starts the fitting too.)
    let mut first_foot = None;
    if surface.parameters(edge.curve.point(edge.first)).is_none() {
        use crate::geometry::extrema::{nearest_on_surface, newton_on_surface};
        use crate::geometry::ParametricSurface;
        // (On a cylinder or cone written as a B-spline, in closed form.)
        let circular = pr::Circular::of(surface);
        let mut start = None;
        for i in 0..=4 {
            let p = edge.curve.point(edge.first + (edge.last - edge.first) * i as f64 / 4.0);
            let exact = circular.as_ref().and_then(|c| c.parameters(p)).map(|(u, v)| (u, v, surface.point(u, v).distance(p)));
            let near = exact.or_else(|| start.map(|s| newton_on_surface(p, surface, s)).filter(|f| f.2 <= 1e-2));
            let (u, v, d) = near.unwrap_or_else(|| nearest_on_surface(p, surface, start));
            start = Some((u, v));
            if i == 0 {
                first_foot = start;
            }
            if d > 1e-2 {
                return Err(invalid("an edge that does not project onto its face".into()));
            }
        }
    }
    // On a free-form surface (inverted numerically) a coarser search: the file's pcurve
    // serves where this fails.
    let freeform = surface.parameters(edge.curve.point(edge.first)).is_none();
    // (Files put edges up to a hundredth off their free-form faces, as OCCT's tolerances
    // after reading show: the edge's tolerance then covers the gap.)
    let tiers: &[f64] = if freeform { &[1e-5, 1e-3, 1e-2] } else { &[1e-7, 1e-5, 1e-3] };
    // The polyline (fine, adaptive, costly) is no use for a curve off the surface by more
    // than the loosest fit allows: the file's pcurve serves there if it gives one.
    let off_surface = || {
        (0..=8).any(|i| {
            let p = edge.curve.point(edge.first + (edge.last - edge.first) * i as f64 / 8.0);
            surface.parameters(p).is_some_and(|(u, v)| surface.point(u, v).distance(p) > 1e-3)
        })
    };
    // Last, for an edge off its face by more than the fits allow (files accurate to a
    // hundredth or so, which OCCT reads with tolerances to match): fitted at the gap measured,
    // the edge's tolerance then covering it. (Not for one off by a tenth of its length: no
    // edge of the face.)
    let at_gap = || {
        if has_pcurve {
            return None;
        }
        let gap = gap_to(edge, surface, first_foot)?;
        let length = (0..16)
            .map(|i| {
                let t = |k: f64| edge.first + (edge.last - edge.first) * k / 16.0;
                edge.curve.point(t(i as f64)).distance(edge.curve.point(t(i as f64 + 1.0)))
            })
            .sum::<f64>();
        if gap > 0.1 * length {
            return None;
        }
        let tolerance = 1.5 * gap + 1e-7;
        pr::pcurve_on_range_from(edge.curve, surface, edge.first, edge.last, &[tolerance], first_foot)
            .or_else(|| pr::polyline_pcurve_on_range(edge.curve, surface, edge.first, edge.last, tolerance))
    };
    pr::pcurve_on_range_from(edge.curve, surface, edge.first, edge.last, tiers, first_foot)
        .or_else(|| if freeform || (has_pcurve && off_surface()) { None } else { pr::polyline_pcurve_on_range(edge.curve, surface, edge.first, edge.last, 1e-6) })
        .or_else(at_gap)
        .map(Curve2::BSpline)
        .ok_or_else(|| invalid("an edge that does not project onto its face".into()))
}


/// How far the edge strays from the surface: the largest distance of its points (64, each
/// foot found from the one before) from their feet; `None` where a foot is not found.
fn gap_to(edge: &CurveSpan, surface: &Surface, start: Option<(f64, f64)>) -> Option<f64> {
    use crate::geometry::extrema::nearest_on_surface;
    use crate::geometry::ParametricSurface;
    let mut previous = start;
    let mut gap: f64 = 0.0;
    for i in 0..=64 {
        let p = edge.curve.point(edge.first + (edge.last - edge.first) * i as f64 / 64.0);
        let (u, v, d) = nearest_on_surface(p, surface, previous);
        if !d.is_finite() {
            return None;
        }
        previous = Some((u, v));
        gap = gap.max(surface.point(u, v).distance(p));
    }
    Some(gap)
}

/// Pcurves moved by whole periods of a periodic surface so that each starts where the one
/// before it ends (the file's and projected pcurves may lie in different periods; a seam's
/// two uses end up on its two sides).
pub(crate) fn align_periods(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface) {
    latitudes_into_range(coedges, surface);
    let (pu, pv) = (surface.u_period(), surface.v_period());
    if (pu.is_none() && pv.is_none()) || coedges.is_empty() {
        return;
    }
    let ends = |c: &(EdgeId, bool, Curve2, f64, f64)| {
        let (a, b) = (c.2.point(c.3), c.2.point(c.4));
        if c.1 { (a, b) } else { (b, a) }
    };
    let steps = |p: Option<f64>| -> Vec<f64> { p.map_or(vec![0.0], |p| (-3..=3).map(|k| k as f64 * p).collect()) };
    // (Coedge i moved by whole periods to start where the one before it ends.)
    let follow = |coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], i: usize| {
        let (_, previous_end) = ends(&coedges[i - 1]);
        let (start, _) = ends(&coedges[i]);
        let mut best = (0.0, 0.0, start.distance(previous_end));
        for du in steps(pu) {
            for dv in steps(pv) {
                let d = Point2::new(start.x + du, start.y + dv).distance(previous_end);
                if d < best.2 - 1e-12 {
                    best = (du, dv, d);
                }
            }
        }
        if best.0 != 0.0 || best.1 != 0.0 {
            let t = crate::math::Transform2::translation(crate::math::Vector2::new(best.0, best.1));
            use crate::math::Transformable2;
            coedges[i].2 = coedges[i].2.transformed(&t);
        }
    };
    for i in 1..coedges.len() {
        follow(coedges, i);
    }
    // A seam's two uses on its two sides: where the loop leaves out the edge collapsed at a
    // pole between them (a meridian's end on the axis), the second joined the first across
    // the pole, on the same side; moved a period over, whichever way meets its neighbours
    // better (the pole's edge then added across the gap).
    // (Or two edges along one meridian, where a file models the seam so.)
    let n = coedges.len();
    let middle = |c: &(EdgeId, bool, Curve2, f64, f64)| c.2.point(0.5 * (c.3 + c.4));
    for i in 0..n {
        let Some(_) = (0..i).find(|&j| middle(&coedges[i]).distance(middle(&coedges[j])) <= 1e-9) else { continue };
        let (start, end) = ends(&coedges[i]);
        let (previous, next) = (ends(&coedges[(i + n - 1) % n]).1, ends(&coedges[(i + 1) % n]).0);
        let shifts: Vec<crate::math::Vector2> = [pu.map(|p| crate::math::Vector2::new(p, 0.0)), pv.map(|p| crate::math::Vector2::new(0.0, p))]
            .into_iter()
            .flatten()
            .flat_map(|d| [d, -d])
            .collect();
        let gap = |d: crate::math::Vector2| (start + d).distance(previous).min((end + d).distance(next));
        if let Some(d) = shifts.into_iter().min_by(|a, b| gap(*a).total_cmp(&gap(*b))) {
            use crate::math::Transformable2;
            coedges[i].2 = coedges[i].2.transformed(&crate::math::Transform2::translation(d));
            // (Those after it following it there.)
            for k in i + 1..n {
                follow(coedges, k);
            }
        }
    }
    // A loop still ending a whole period from where it starts, round an edge along a pole
    // (where every u is one point, so that a pcurve there may run either way): that edge's
    // pcurve run on by the period, those after it moved over by it.
    if n > 1 {
        let (start, _) = ends(&coedges[0]);
        let (_, end) = ends(&coedges[n - 1]);
        let gap = start - end;
        let whole = |d: f64, period: Option<f64>| period.is_some_and(|p| ((d / p).round() * p - d).abs() <= 1e-9 * p.max(1.0) && (d / p).round() != 0.0);
        let along_pole = |c: &(EdgeId, bool, Curve2, f64, f64)| {
            use crate::geometry::ParametricSurface;
            let (a, b) = (c.2.point(c.3), c.2.point(c.4));
            let m = c.2.point(0.5 * (c.3 + c.4));
            let (pa, pb, pm) = (surface.point(a.x, a.y), surface.point(b.x, b.y), surface.point(m.x, m.y));
            let size = 1.0 + pa.x.abs().max(pa.y.abs()).max(pa.z.abs());
            a.distance(b) > 1e-9 && pa.distance(pb) <= 1e-9 * size && pa.distance(pm) <= 1e-9 * size
        };
        let shifted = (whole(gap.x, pu) && gap.y.abs() <= 1e-9) || (whole(gap.y, pv) && gap.x.abs() <= 1e-9);
        // (Or two consecutive edges meeting at a pole, the edge along it left out: those after
        // the pole moved over, the pole's edge added across the gap then running round.)
        let at_pole = |i: usize| {
            use crate::geometry::ParametricSurface;
            let (_, end) = ends(&coedges[i]);
            let p = surface.point(end.x, end.y);
            let size = 1.0 + p.x.abs().max(p.y.abs()).max(p.z.abs());
            let step = |d: f64, period: Option<f64>| period.map_or(d, |p| 0.25 * p);
            surface.point(end.x + step(1.0, pu), end.y).distance(p) <= 1e-9 * size
        };
        if shifted && !(0..n).any(|i| along_pole(&coedges[i])) {
            if let Some(i) = (0..n - 1).find(|&i| at_pole(i)) {
                use crate::math::Transformable2;
                let t = crate::math::Transform2::translation(gap);
                for k in i + 1..n {
                    coedges[k].2 = coedges[k].2.transformed(&t);
                }
            }
        } else if shifted {
            if let Some(i) = (0..n).find(|&i| along_pole(&coedges[i])) {
                use crate::math::Transformable2;
                let (c, first, last) = (&coedges[i].2, coedges[i].3, coedges[i].4);
                let (mut a, mut b) = (c.point(first), c.point(last));
                // (Its end where the loop leaves it: the curve's last point run forward, its
                // first run backward.)
                if coedges[i].1 {
                    b = b + gap;
                } else {
                    a = a + gap;
                }
                coedges[i].2 = Curve2::BSpline(crate::geometry::curve2::BSplineCurve2::new(vec![a, b], vec![first, last], vec![2, 2], 1));
                let t = crate::math::Transform2::translation(gap);
                for k in i + 1..n {
                    coedges[k].2 = coedges[k].2.transformed(&t);
                }
            }
        }
    }
}


/// The ranges of `u` and `v` the pcurves of a face's loops span, leaving out those along the
/// seam of a surface closed but not periodic and those of edges collapsed to a point (whose
/// place [`seams_to_sides`] decides).
pub(crate) fn face_spans(loops: &[&[(EdgeId, bool, Curve2, f64, f64)]], surface: &Surface) -> [(f64, f64); 2] {
    use crate::geometry::ParametricSurface;
    let (u0, u1, v0, v1) = surface.bounds();
    let size = surface.point(u0.max(-1e3), v0.max(-1e3)).distance(surface.point(u1.min(1e3), v1.min(1e3))) + 1.0;
    let mut spans = [(f64::INFINITY, f64::NEG_INFINITY); 2];
    for coedges in loops {
        for c in coedges.iter() {
            let pts: Vec<Point2> = [0.0, 0.5, 1.0].iter().map(|s| c.2.point(c.3 + (c.4 - c.3) * s)).collect();
            let images: Vec<Point> = pts.iter().map(|q| surface.point(q.x, q.y)).collect();
            if images.iter().all(|q| q.distance(images[0]) <= 1e-9 * size) {
                continue;
            }
            for (k, (lo, hi)) in [(0, (u0, u1)), (1, (v0, v1))] {
                let at = 1e-7 * (hi - lo).abs().max(1.0);
                let xs: Vec<f64> = pts.iter().map(|q| if k == 0 { q.x } else { q.y }).collect();
                let on_seam = lo.is_finite() && hi.is_finite() && (xs.iter().all(|x| (x - lo).abs() <= at) || xs.iter().all(|x| (x - hi).abs() <= at));
                if on_seam {
                    continue;
                }
                for x in xs {
                    spans[k] = (spans[k].0.min(x), spans[k].1.max(x));
                }
            }
        }
    }
    spans
}


/// On a surface closed but not periodic (its parameter range wrapping round once), a pcurve
/// along the seam (on one end of the range) put on the end the loop lies on: by the range
/// the loop's other pcurves span, or (a seam used twice) where the pcurves before and after
/// it meet it. Pcurves of edges collapsed to a point of the surface (a pole), whose
/// parameters the projection leaves arbitrary, then join their neighbours.
pub(crate) fn seams_to_sides(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface, spans: [(f64, f64); 2]) -> Vec<(usize, crate::math::Vector2)> {
    use crate::geometry::ParametricSurface;
    use crate::math::Transformable2;
    let mut ties = Vec::new();
    let n = coedges.len();
    if n < 2 {
        return ties;
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let closed = [
        surface.u_period().is_none() && surface.is_u_closed() && u0.is_finite() && u1.is_finite(),
        surface.v_period().is_none() && surface.is_v_closed() && v0.is_finite() && v1.is_finite(),
    ];
    if !closed[0] && !closed[1] {
        return ties;
    }
    let ends = |c: &(EdgeId, bool, Curve2, f64, f64)| {
        let (a, b) = (c.2.point(c.3), c.2.point(c.4));
        if c.1 { (a, b) } else { (b, a) }
    };
    let samples = |c: &(EdgeId, bool, Curve2, f64, f64)| -> Vec<Point2> { [0.0, 0.5, 1.0].iter().map(|s| c.2.point(c.3 + (c.4 - c.3) * s)).collect() };
    // Edges collapsed to a point on the surface.
    let size = surface.point(u0.max(-1e3), v0.max(-1e3)).distance(surface.point(u1.min(1e3), v1.min(1e3))) + 1.0;
    let collapsed: Vec<bool> = coedges
        .iter()
        .map(|c| {
            let pts: Vec<Point> = samples(c).into_iter().map(|q| surface.point(q.x, q.y)).collect();
            pts.iter().all(|q| q.distance(pts[0]) <= 1e-9 * size)
        })
        .collect();
    for (k, (lo, hi)) in [(0, (u0, u1)), (1, (v0, v1))] {
        if !closed[k] {
            continue;
        }
        let at = 1e-7 * (hi - lo).abs().max(1.0);
        let coord = |p: Point2| if k == 0 { p.x } else { p.y };
        let side_of = |c: &(EdgeId, bool, Curve2, f64, f64)| {
            let xs: Vec<f64> = samples(c).into_iter().map(coord).collect();
            if xs.iter().all(|x| (x - lo).abs() <= at) {
                Some(lo)
            } else if xs.iter().all(|x| (x - hi).abs() <= at) {
                Some(hi)
            } else {
                None
            }
        };
        let seams: Vec<usize> = (0..n).filter(|&i| !collapsed[i] && side_of(&coedges[i]).is_some()).collect();
        let span = spans[k];
        let same_edge_twice = seams.iter().any(|&i| seams.iter().any(|&j| j != i && coedges[j].0 == coedges[i].0));
        for &i in &seams {
            let side = side_of(&coedges[i]).expect("a seam");
            let other = if side == lo { hi } else { lo };
            let shift = if k == 0 { crate::math::Vector2::new(other - side, 0.0) } else { crate::math::Vector2::new(0.0, other - side) };
            let better = if span.0 <= span.1 && !same_edge_twice {
                let width = |x: f64| span.1.max(x) - span.0.min(x);
                if (width(other) - width(side)).abs() <= at {
                    ties.push((i, shift));
                }
                width(other) < width(side) - at
            } else {
                // By the neighbours that are neither seams (a seam's other use, which may
                // move too) nor collapsed.
                let fixed = |j: usize| !collapsed[j] && !seams.contains(&j);
                let (start, end) = ends(&coedges[i]);
                let (previous, next) = ((i + n - 1) % n, (i + 1) % n);
                let gap = |shift: f64| {
                    let moved = |p: Point2| if k == 0 { Point2::new(p.x + shift, p.y) } else { Point2::new(p.x, p.y + shift) };
                    let before = if fixed(previous) { moved(start).distance(ends(&coedges[previous]).1) } else { 0.0 };
                    let after = if fixed(next) { moved(end).distance(ends(&coedges[next]).0) } else { 0.0 };
                    before + after
                };
                gap(other - side) < gap(0.0) - 1e-9
            };
            if better {
                coedges[i].2 = coedges[i].2.transformed(&crate::math::Transform2::translation(shift));
            }
        }
    }
    // Collapsed edges from the end of the pcurve before to the start of the one after.
    for i in 0..n {
        if !collapsed[i] || collapsed[(i + n - 1) % n] || collapsed[(i + 1) % n] {
            continue;
        }
        let (_, from) = ends(&coedges[(i + n - 1) % n]);
        let (to, _) = ends(&coedges[(i + 1) % n]);
        let (first, last) = (coedges[i].3, coedges[i].4);
        // (Along the edge's parameter: from where the coedge starts if it runs along it.)
        let (a, b) = if coedges[i].1 { (from, to) } else { (to, from) };
        if let Ok(line) = BSplineCurve2::try_new(vec![a, b], vec![first, last], vec![2, 2], 1) {
            coedges[i].2 = Curve2::BSpline(line);
        }
    }
    ties
}


/// The (u, v) box of a loop's pcurves (by their ends and middles).
pub(crate) fn uv_box(coedges: &[(EdgeId, bool, Curve2, f64, f64)]) -> (Point2, Point2) {
    let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
    for (_, _, pcurve, first, last) in coedges {
        for s in [0.0, 0.5, 1.0] {
            let p = pcurve.point(first + (last - first) * s);
            lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
        }
    }
    (lo, hi)
}


/// Pcurves on a sphere moved by whole turns of latitude into its range (-π/2..π/2): files
/// write some a turn or more away (-2π..-3π/2 for 0..π/2), the same points.
fn latitudes_into_range(coedges: &mut [Coedge], surface: &Surface) {
    if !matches!(surface, Surface::Sphere(_)) {
        return;
    }
    use crate::math::Transformable2;
    let turn = std::f64::consts::TAU;
    let half = std::f64::consts::FRAC_PI_2 + 1e-6;
    for c in coedges.iter_mut() {
        let (a, b) = (c.2.point(c.3), c.2.point(c.4));
        let middle = c.2.point(0.5 * (c.3 + c.4)).y;
        let (lo, hi) = (a.y.min(b.y).min(middle), a.y.max(b.y).max(middle));
        if lo >= -half && hi <= half {
            continue;
        }
        let k = (-middle / turn).round();
        if k != 0.0 && lo + k * turn >= -half && hi + k * turn <= half {
            c.2 = c.2.transformed(&crate::math::Transform2::translation(crate::math::Vector2::new(0.0, k * turn)));
        }
    }
}


/// Moves a loop by whole periods so that its middle lies within `target`'s range (a hole
/// read in another period than its face's outer boundary).
pub(crate) fn into_period(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], target: (Point2, Point2), surface: &Surface) {
    let (lo, hi) = uv_box(coedges);
    let middle = Point2::new(0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y));
    let shift = |period: Option<f64>, at: f64, from: f64, to: f64| {
        period.map_or(0.0, |p| if at >= from && at <= to { 0.0 } else { (((from + to) * 0.5 - at) / p).round() * p })
    };
    let du = shift(surface.u_period(), middle.x, target.0.x, target.1.x);
    let dv = shift(surface.v_period(), middle.y, target.0.y, target.1.y);
    if du != 0.0 || dv != 0.0 {
        use crate::math::Transformable2;
        let t = crate::math::Transform2::translation(crate::math::Vector2::new(du, dv));
        for c in coedges.iter_mut() {
            c.2 = c.2.transformed(&t);
        }
    }
}


type Coedge = (EdgeId, bool, Curve2, f64, f64);

/// A coedge's ends in the surface's parameters, in its direction.
fn coedge_ends(c: &Coedge) -> (Point2, Point2) {
    let (a, b) = (c.2.point(c.3), c.2.point(c.4));
    if c.1 { (a, b) } else { (b, a) }
}

/// How far a loop's end is from its start in the surface's parameters.
fn loop_shift(coedges: &[Coedge]) -> Point2 {
    match (coedges.first(), coedges.last()) {
        (Some(first), Some(last)) => {
            let (a, b) = (coedge_ends(first).0, coedge_ends(last).1);
            Point2::new(b.x - a.x, b.y - a.y)
        }
        _ => Point2::new(0.0, 0.0),
    }
}

fn component(p: Point2, axis: usize) -> f64 {
    if axis == 0 { p.x } else { p.y }
}

/// The parameter (0 for u, 1 for v) whose period some of the loops go round, open in the
/// parameters, if only one is.
fn round_period(loops: &[(bool, Vec<Coedge>)], surface: &Surface) -> Option<usize> {
    let periods = [surface.u_period(), surface.v_period()];
    let open: Vec<usize> = (0..2)
        .filter(|&axis| periods[axis].is_some_and(|p| loops.iter().any(|(_, c)| component(loop_shift(c), axis).abs() > 0.5 * p)))
        .collect();
    (open.len() == 1).then(|| open[0])
}

/// True if the loops going round the period along `axis` go round it as often each way.
fn paired(loops: &[(bool, Vec<Coedge>)], surface: &Surface, axis: usize) -> bool {
    let period = if axis == 0 { surface.u_period() } else { surface.v_period() };
    let Some(period) = period else { return false };
    let turns: f64 = loops.iter().map(|(_, c)| (component(loop_shift(c), axis) / period).round()).sum();
    turns == 0.0
}

/// The area a loop adds to its face's in the surface's parameters, positive
/// counter-clockwise: -∮v du along a u period (`axis` 0), ∮u dv along a v period. Iso-lines
/// across that period (seams) add nothing, so the loops of a band left without its seam
/// sum to the band's area.
fn band_area(coedges: &[Coedge], axis: usize) -> f64 {
    let mut points = Vec::new();
    for (_, forward, pcurve, first, last) in coedges {
        for i in 0..=16 {
            let s = i as f64 / 16.0;
            let t = if *forward { first + (last - first) * s } else { last - (last - first) * s };
            points.push(pcurve.point(t));
        }
    }
    points
        .windows(2)
        .map(|w| if axis == 0 { -0.5 * (w[0].y + w[1].y) * (w[1].x - w[0].x) } else { 0.5 * (w[0].x + w[1].x) * (w[1].y - w[0].y) })
        .sum()
}

/// Two loops going round a period opposite ways (a band of a cylinder between two circles)
/// joined into one by the seam the file left out, as OCCT's `ShapeFix_Face::FixMissingSeam`:
/// the first loop, the seam across to the second's vertex at the same place on the period,
/// the second loop, the seam back. Left as they are unless the second loop has a vertex
/// there: one of its edges is split there (in the faces read before as well). Returns
/// the joined loop's index.
fn add_missing_seam(model: &mut Model, gaps: &mut Gaps, loops: &mut Vec<(bool, Vec<Coedge>)>, surface: &Surface, axis: usize) -> Result<Option<usize>> {
    let Some(period) = (if axis == 0 { surface.u_period() } else { surface.v_period() }) else { return Ok(None) };
    let open: Vec<usize> = (0..loops.len()).filter(|&i| component(loop_shift(&loops[i].1), axis).abs() > 0.5 * period).collect();
    let [i, j] = open[..] else { return Ok(None) };
    let (sa, sb) = (loop_shift(&loops[i].1), loop_shift(&loops[j].1));
    if component(sa, axis) * component(sb, axis) >= 0.0 {
        return Ok(None);
    }
    let a0 = coedge_ends(&loops[i].1[0]).0;
    let a1 = Point2::new(a0.x + sa.x, a0.y + sa.y);
    // The second loop started at its vertex across from the first's start.
    let tolerance = 1e-7 * period.max(1.0);
    let target = component(a1, axis);
    let off = |x: f64| {
        let d = x - target;
        d - (d / period).round() * period
    };
    let at_vertex = |b: &[Coedge]| (0..b.len()).find(|&k| off(component(coedge_ends(&b[k]).0, axis)).abs() <= tolerance);
    if at_vertex(&loops[j].1).is_none() {
        // A vertex made there, splitting the edge that passes it.
        let crossing = loops[j].1.iter().enumerate().find_map(|(k, c)| {
            let f = |t: f64| off(component(c.2.point(t), axis));
            let n = 32;
            (0..n).find_map(|i| {
                let (mut lo, mut hi) = (c.3 + (c.4 - c.3) * i as f64 / n as f64, c.3 + (c.4 - c.3) * (i + 1) as f64 / n as f64);
                let (flo, fhi) = (f(lo), f(hi));
                // (A change of sign, not the jump half a period away.)
                if flo * fhi > 0.0 || (flo - fhi).abs() > 0.25 * period {
                    return None;
                }
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if f(mid) * flo > 0.0 { lo = mid } else { hi = mid }
                }
                let t = 0.5 * (lo + hi);
                let span = (c.4 - c.3).abs();
                ((t - c.3).abs() > 1e-9 * span && (c.4 - t).abs() > 1e-9 * span).then_some((k, t))
            })
        });
        let Some((k, t)) = crossing else { return Ok(None) };
        let c = loops[j].1[k].clone();
        let (h1, h2) = gaps.split(model, c.0, t)?;
        let halves = if c.1 { [(h1, c.1, c.2.clone(), c.3, t), (h2, c.1, c.2, t, c.4)] } else { [(h2, c.1, c.2.clone(), t, c.4), (h1, c.1, c.2, c.3, t)] };
        loops[j].1.splice(k..=k, halves);
    }
    let b = &loops[j].1;
    let Some(k) = at_vertex(b) else { return Ok(None) };
    use crate::math::Transformable2;
    let shift = |p: Point2| crate::math::Transform2::translation(crate::math::Vector2::new(p.x, p.y));
    let mut second: Vec<Coedge> = b[k..].to_vec();
    second.extend(b[..k].iter().map(|c| (c.0, c.1, c.2.transformed(&shift(sb)), c.3, c.4)));
    let b0 = coedge_ends(&second[0]).0;
    let d = component(b0, axis) - component(a1, axis);
    let along = -(d / period).round() * period;
    let offset = if axis == 0 { Point2::new(along, 0.0) } else { Point2::new(0.0, along) };
    for c in &mut second {
        c.2 = c.2.transformed(&shift(offset));
    }
    let b0 = coedge_ends(&second[0]).0;
    let (ta, tb) = (component(a0, 1 - axis), component(b0, 1 - axis));
    if (ta - tb).abs() <= tolerance {
        return Ok(None);
    }
    // The vertices the seam joins.
    let start_vertex = |c: &Coedge| {
        let (s, e) = model.edge_vertices(Shape::from(c.0));
        VertexId::try_from(if c.1 { s } else { e }.node).expect("a vertex")
    };
    let (va, vb) = (start_vertex(&loops[i].1[0]), start_vertex(&second[0]));
    let low_first = ta < tb;
    let edge = if low_first { model.add_edge(va, vb, crate::tolerance::LINEAR)? } else { model.add_edge(vb, va, crate::tolerance::LINEAR)? };
    let at = component(a0, axis);
    let cut = crate::divide::Cut::iso(surface, axis == 0, at, ta.min(tb), ta.max(tb));
    model.set_edge_curve(edge.into(), cut.curve, crate::tolerance::LINEAR)?;
    model.set_edge_range(edge.into(), cut.first, cut.last)?;
    let iso = |at: f64| crate::divide::Cut::iso(surface, axis == 0, at, 0.0, 1.0).pcurve;
    // Across from the first loop's end to the second's start, and back at the first's start.
    let across = (edge, low_first, iso(component(a1, axis)), cut.first, cut.last);
    let back = (edge, !low_first, iso(at), cut.first, cut.last);
    let mut joined = std::mem::take(&mut loops[i].1);
    joined.push(across);
    joined.extend(second);
    joined.push(back);
    let outer = loops[i].0 || loops[j].0;
    loops[i] = (outer, joined);
    loops.remove(j);
    Ok(Some(i))
}

/// The signed area a loop encloses in the surface's parameters (positive counter-clockwise).
pub(crate) fn uv_area(coedges: &[(EdgeId, bool, Curve2, f64, f64)]) -> f64 {
    let mut points = Vec::new();
    for (_, forward, pcurve, first, last) in coedges {
        for i in 0..16 {
            let s = i as f64 / 16.0;
            let t = if *forward { first + (last - first) * s } else { last - (last - first) * s };
            points.push(pcurve.point(t));
        }
    }
    let n = points.len();
    (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum::<f64>() * 0.5
}

