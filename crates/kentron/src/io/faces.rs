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
}

impl Gaps {

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
    let mut placed = Vec::with_capacity(loops.len());
    for (outer, mut coedges) in loops {
        seams_to_sides(&mut coedges, &surface, spans);
        placed.push((outer, close_at_poles(model, coedges, &surface)?));
    }
    let mut loops = placed;
    // Outer boundaries run counter-clockwise in the surface's parameters.
    let areas: Vec<f64> = loops.iter().map(|(_, c)| uv_area(c)).collect();
    let outer = loops.iter().position(|(o, _)| *o).or_else(|| (0..areas.len()).max_by(|&i, &j| areas[i].abs().total_cmp(&areas[j].abs())));
    if let Some(i) = outer {
        let scale = areas.iter().map(|a| a.abs()).fold(0.0, f64::max);
        if areas[i] < -1e-9 * scale.max(1e-300) {
            for (_, coedges) in &mut loops {
                coedges.reverse();
                for c in coedges.iter_mut() {
                    c.1 = !c.1;
                }
            }
        }
    }
    // Holes in the outer boundary's period of a periodic surface.
    if let Some(i) = outer {
        let outer_box = uv_box(&loops[i].1);
        for (j, (_, coedges)) in loops.iter_mut().enumerate() {
            if j != i {
                into_period(coedges, outer_box, &surface);
            }
        }
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
    Ok(if same_sense { face } else { face.reversed() })
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
    if let Some(pc) = pr::analytic_pcurve(edge.curve, surface) {
        return Ok(pc);
    }
    // An edge well off the surface (bad data) fails at once, before costly fitting. (Each
    // point's foot from the one before, the whole surface searched only to confirm one off.)
    // (The first point's foot, from the whole surface, starts the fitting too.)
    let mut first_foot = None;
    if surface.parameters(edge.curve.point(edge.first)).is_none() {
        use crate::geometry::extrema::{nearest_on_surface, newton_on_surface};
        let mut start = None;
        for i in 0..=4 {
            let p = edge.curve.point(edge.first + (edge.last - edge.first) * i as f64 / 4.0);
            let near = start.map(|s| newton_on_surface(p, surface, s)).filter(|f| f.2 <= 1e-2);
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
    pr::pcurve_on_range_from(edge.curve, surface, edge.first, edge.last, tiers, first_foot)
        .or_else(|| if freeform || (has_pcurve && off_surface()) { None } else { pr::polyline_pcurve_on_range(edge.curve, surface, edge.first, edge.last, 1e-6) })
        .map(Curve2::BSpline)
        .ok_or_else(|| invalid("an edge that does not project onto its face".into()))
}


/// Pcurves moved by whole periods of a periodic surface so that each starts where the one
/// before it ends (the file's and projected pcurves may lie in different periods; a seam's
/// two uses end up on its two sides).
pub(crate) fn align_periods(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface) {
    let (pu, pv) = (surface.u_period(), surface.v_period());
    if (pu.is_none() && pv.is_none()) || coedges.is_empty() {
        return;
    }
    let ends = |c: &(EdgeId, bool, Curve2, f64, f64)| {
        let (a, b) = (c.2.point(c.3), c.2.point(c.4));
        if c.1 { (a, b) } else { (b, a) }
    };
    let steps = |p: Option<f64>| -> Vec<f64> { p.map_or(vec![0.0], |p| (-3..=3).map(|k| k as f64 * p).collect()) };
    for i in 1..coedges.len() {
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
pub(crate) fn seams_to_sides(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface, spans: [(f64, f64); 2]) {
    use crate::geometry::ParametricSurface;
    use crate::math::Transformable2;
    let n = coedges.len();
    if n < 2 {
        return;
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let closed = [
        surface.u_period().is_none() && surface.is_u_closed() && u0.is_finite() && u1.is_finite(),
        surface.v_period().is_none() && surface.is_v_closed() && v0.is_finite() && v1.is_finite(),
    ];
    if !closed[0] && !closed[1] {
        return;
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
            let better = if span.0 <= span.1 && !same_edge_twice {
                let width = |x: f64| span.1.max(x) - span.0.min(x);
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
                let shift = if k == 0 { crate::math::Vector2::new(other - side, 0.0) } else { crate::math::Vector2::new(0.0, other - side) };
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

