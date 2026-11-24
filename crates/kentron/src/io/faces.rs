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
            // (A band round a surface closed but not periodic (a closed B-spline's): by its
            // loops' areas round that range, as a periodic band's.)
            None => match closed_band(loops, surface) {
                Some(axis) => loops.iter().map(|(_, c)| band_area(c, axis)).sum::<f64>() < 0.0,
                None => areas[i] < -1e-9 * scale.max(1e-300),
            },
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


/// The loop on a sphere with each edge along a great circle through a pole split there (in
/// the faces read before as well), as OCCT's `ShapeFix` splits it: no one pcurve runs over a
/// pole, its longitude jumping there. Each half is a meridian's arc, its pcurve a line of
/// fixed longitude; [`close_at_poles`] then joins the halves round the pole. Files leave
/// such edges whole (a corner of rounds meeting at a sphere's pole).
pub(crate) fn split_at_poles(model: &mut Model, gaps: &mut Gaps, coedges: Vec<Coedge>, surface: &Surface, precision: f64) -> Result<Vec<Coedge>> {
    let Surface::Sphere(sphere) = surface else { return Ok(coedges) };
    let (center, r) = (sphere.frame().origin(), sphere.radius());
    let z = crate::math::Vector::from(sphere.frame().z());
    let poles = [center + z * r, center - z * r];
    let tolerance = 1e-9 * (1.0 + r + (center - Point::ORIGIN).length());
    let mut out = Vec::with_capacity(coedges.len());
    for c in coedges {
        let Some((curve, _, _)) = model.edge_curve(Shape::from(c.0)) else {
            out.push(c);
            continue;
        };
        let (lo, hi) = (c.3.min(c.4), c.3.max(c.4));
        // (A great circle: centred on the sphere's, of its radius; the latitude then moves as
        // the circle's angle does. Another curve (a circle written as a rational B-spline)
        // split too where it passes a pole, its pieces projected.)
        let great = matches!(&curve, Curve::Circle(circle) if circle.center().distance(center) <= tolerance && (circle.radius() - r).abs() <= tolerance);
        // (Within the tolerance a file's B-spline keeps to the sphere.)
        let near = if great { tolerance } else { 1e-6 * (1.0 + r) };
        let mut at: Vec<f64> = poles
            .iter()
            .filter_map(|&pole| {
                let foot = crate::geometry::extrema::nearest_on_curve(pole, &curve, lo, hi);
                let inside = foot.parameter - lo > 1e-9 * (hi - lo) && hi - foot.parameter > 1e-9 * (hi - lo);
                (foot.distance <= near && inside).then_some(foot.parameter)
            })
            .collect();
        at.sort_by(f64::total_cmp);
        // The pieces between, each a meridian's arc: its pcurve the line of its longitude.
        let meridian = |a: f64, b: f64| -> Option<Curve2> {
            let (_, va) = sphere.parameters(curve.point(a));
            let (_, vb) = sphere.parameters(curve.point(b));
            let (u, _) = sphere.parameters(curve.point(0.5 * (a + b)));
            let sense = (vb - va).signum();
            let along = ((vb - va).abs() - (b - a).abs()).abs() <= 1e-6 * (b - a).abs().max(1.0);
            along.then(|| Curve2::Line(crate::geometry::curve2::Line2::new(Point2::new(u, va - sense * a), crate::math::Direction2::new(0.0, sense))))
        };
        let bounds: Vec<f64> = std::iter::once(lo).chain(at.iter().copied()).chain(std::iter::once(hi)).collect();
        let pcurves: Option<Vec<Curve2>> = if great {
            bounds.windows(2).map(|w| meridian(w[0], w[1])).collect()
        } else {
            bounds.windows(2).map(|w| project(&CurveSpan { curve: &curve, first: w[0], last: w[1] }, surface, false, precision).ok()).collect()
        };
        let (Some(pcurves), false) = (pcurves, at.is_empty()) else {
            out.push(c);
            continue;
        };
        // (Split one point after another, each time the piece past the last split.)
        let mut pieces = Vec::with_capacity(pcurves.len());
        let mut rest = c.0;
        for (k, &t) in at.iter().enumerate() {
            let (h1, h2) = gaps.split(model, rest, t)?;
            pieces.push((h1, c.1, pcurves[k].clone(), bounds[k], t));
            rest = h2;
        }
        pieces.push((rest, c.1, pcurves[at.len()].clone(), at[at.len() - 1], hi));
        if !c.1 {
            pieces.reverse();
        }
        out.extend(pieces);
    }
    Ok(out)
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
        // (Its tolerance, and its vertex's, covering how far the vertex is from the point on
        // the surface, as OCCT's FixVertexTolerance raises them.)
        let point = model.point(vertex);
        let gap = 1.0000001 * point.distance(p).max(point.distance(q));
        let tolerance = gap.max(crate::tolerance::LINEAR);
        if tolerance > model[v].tolerance {
            model.set_vertex_point(v.into(), point, tolerance)?;
        }
        let edge = model.add_edge(v, v, tolerance)?;
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
/// back on), `precision` the file's (OCCT's reading precision). An edge off the surface is
/// projected all the same, as OCCT's readers do, its tolerance then covering the gap.
///
/// Off planes, as OCCT's readers add pcurves (`ShapeFix_Edge::FixAddPCurve`): by its
/// `ShapeConstruct_ProjectCurveOnSurface`, ported, the edge's vertices' tolerances the
/// precision; then remade where it strays from the curve at the same parameters
/// ([`same_parameter`](crate::geometry::same_parameter::same_parameter)).
pub(crate) fn project(edge: &CurveSpan, surface: &Surface, has_pcurve: bool, precision: f64) -> Result<Curve2> {
    // (An edge running its curve backwards: the pcurve over the range the other way round, at
    // the same parameters.)
    if edge.first > edge.last {
        return project(&CurveSpan { curve: edge.curve, first: edge.last, last: edge.first }, surface, has_pcurve, precision);
    }
    if !matches!(surface, Surface::Plane(_)) {
        let mut projector = crate::geometry::pcurve_projection::PCurveProjector::new(surface, precision);
        let pcurve = projector.perform(edge.curve, edge.first, edge.last, precision, precision).ok_or_else(|| invalid("an edge that does not project onto its face".into()))?;
        // (Then made to follow the curve parameter for parameter, where OCCT's readers run
        // `ShapeFix_Edge::FixSameParameter`: kentron's own way.)
        return Ok(crate::geometry::same_parameter::same_parameter(edge.curve, edge.first, edge.last, surface, &pcurve, crate::tolerance::LINEAR).0);
    }
    project_fitted(edge, surface, has_pcurve)
}

/// The edge's curve projected onto the surface by kentron's own fitting: in closed form
/// where it can be, else fitted to the projected points within tolerances tried in turn
/// (`has_pcurve`: the file gives a pcurve to fall back on).
pub(crate) fn project_fitted(edge: &CurveSpan, surface: &Surface, has_pcurve: bool) -> Result<Curve2> {
    use crate::geometry::project as pr;
    if edge.first > edge.last {
        return project_fitted(&CurveSpan { curve: edge.curve, first: edge.last, last: edge.first }, surface, has_pcurve);
    }
    // (On a sphere, a meridian's arc on the half of its great circle its range lies on.)
    if let Some(pc) = pr::analytic_pcurve_on_range(edge.curve, surface, edge.first, edge.last).or_else(|| pr::analytic_pcurve(edge.curve, surface)) {
        return Ok(pc);
    }
    // An edge well off the surface (bad data, which OCCT reads with tolerances to cover the
    // gap) goes straight to the fit at the gap, before costly fitting. (Each point's foot from
    // the one before, the whole surface searched only to confirm one off.) (The first point's
    // foot, from the whole surface, starts the fitting too.)
    let mut first_foot = None;
    let mut far = false;
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
                far = true;
                break;
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
    // hundredth or so, or worse, which OCCT reads with tolerances to match): fitted at the
    // gap measured, the edge's tolerance then covering it.
    let at_gap = || {
        if has_pcurve {
            return None;
        }
        let gap = gap_to(edge, surface, first_foot)?;
        let tolerance = 1.5 * gap + 1e-7;
        // (The gap is measured at samples; between them the curve may stray further, a few
        // times as far: the fit then allowed that, the edge's tolerance covering it.)
        pr::pcurve_on_range_from(edge.curve, surface, edge.first, edge.last, &[tolerance, 4.0 * tolerance], first_foot)
            .or_else(|| pr::polyline_pcurve_on_range(edge.curve, surface, edge.first, edge.last, tolerance))
    };
    let fitted = if far { None } else {
        pr::pcurve_on_range_from(edge.curve, surface, edge.first, edge.last, tiers, first_foot)
            .or_else(|| if freeform || (has_pcurve && off_surface()) { None } else { pr::polyline_pcurve_on_range(edge.curve, surface, edge.first, edge.last, 1e-6) })
    };
    fitted
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
        // (First the way that closes the loop (its end where it starts, those after following
        // it): round a pole the seam's two uses run the latitudes' period one way and the
        // pole's edge the other; then the way that meets the neighbours better.)
        let closure = |d: crate::math::Vector2| {
            use crate::math::Transformable2;
            let mut trial = coedges.to_vec();
            trial[i].2 = trial[i].2.transformed(&crate::math::Transform2::translation(d));
            for k in i + 1..n {
                follow(&mut trial, k);
            }
            ends(&trial[0]).0.distance(ends(&trial[n - 1]).1)
        };
        let scored: Vec<(crate::math::Vector2, f64, f64)> = shifts.into_iter().map(|d| (d, closure(d), gap(d))).collect();
        if let Some(&(d, _, _)) = scored.iter().min_by(|a, b| {
            let close = |c: f64| (c * 1e6).round();
            close(a.1).total_cmp(&close(b.1)).then(a.2.total_cmp(&b.2))
        }) {
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
        // (The other parameter's gap at rounding level of its values: a generator's range
        // in a file may end a few digits short of the circle it meets.)
        let level = |p: Point2| 1e-9 * (1.0 + p.x.abs().max(p.y.abs()));
        let near = level(start).max(level(end));
        let shifted = (whole(gap.x, pu) && gap.y.abs() <= near) || (whole(gap.y, pv) && gap.x.abs() <= near);
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
        // (Not where the loop itself ends at a pole, a period from its start, as a whole
        // sphere's seam used both ways: the pole's edge added there closes it.)
        if shifted && !(0..n).any(|i| along_pole(&coedges[i])) && !at_pole(n - 1) {
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
            // By the neighbours that are neither seams (a seam's other use, which may move
            // too) nor collapsed.
            let fixed = |j: usize| !collapsed[j] && !seams.contains(&j);
            let (start, end) = ends(&coedges[i]);
            let (previous, next) = ((i + n - 1) % n, (i + 1) % n);
            let gap = |shift: f64| {
                let moved = |p: Point2| if k == 0 { Point2::new(p.x + shift, p.y) } else { Point2::new(p.x, p.y + shift) };
                let before = if fixed(previous) { moved(start).distance(ends(&coedges[previous]).1) } else { 0.0 };
                let after = if fixed(next) { moved(end).distance(ends(&coedges[next]).0) } else { 0.0 };
                before + after
            };
            let better = if span.0 <= span.1 && !same_edge_twice {
                let width = |x: f64| span.1.max(x) - span.0.min(x);
                if (width(other) - width(side)).abs() <= at {
                    // (A tie by the span (the face across the whole range): the side its
                    // neighbours in the loop meet it on, else left to the faces around.)
                    let (moved, stay) = (gap(other - side), gap(0.0));
                    if (moved - stay).abs() > at {
                        moved < stay
                    } else {
                        ties.push((i, shift));
                        false
                    }
                } else {
                    width(other) < width(side) - at
                }
            } else {
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
        // (By more than rounding: a loop between two meridians a half turn apart, run round a
        // pole and back, ends half a period from its start either way.)
        .filter(|&axis| periods[axis].is_some_and(|p| loops.iter().any(|(_, c)| component(loop_shift(c), axis).abs() > 0.5 * p * (1.0 + 1e-9))))
        .collect();
    (open.len() == 1).then(|| open[0])
}

/// The period of a parameter, or for one closed but not periodic (a closed B-spline
/// surface's) its range: what a loop going round the surface that way spans.
fn closure_period(surface: &Surface, axis: usize) -> Option<f64> {
    let (u0, u1, v0, v1) = surface.bounds();
    let (period, closed, lo, hi) = if axis == 0 { (surface.u_period(), surface.is_u_closed(), u0, u1) } else { (surface.v_period(), surface.is_v_closed(), v0, v1) };
    period.or_else(|| (closed && lo.is_finite() && hi.is_finite() && hi > lo).then_some(hi - lo))
}

/// The parameter, closed but not periodic, that two of the loops go round opposite ways
/// (a band), if one is.
fn closed_band(loops: &[(bool, Vec<Coedge>)], surface: &Surface) -> Option<usize> {
    (0..2).find(|&axis| {
        let periodic = if axis == 0 { surface.u_period() } else { surface.v_period() };
        let Some(range) = closure_period(surface, axis).filter(|_| periodic.is_none()) else { return false };
        let turns: Vec<f64> = loops.iter().map(|(_, c)| (component(loop_shift(c), axis) / range).round()).filter(|t| *t != 0.0).collect();
        turns.len() >= 2 && turns.iter().sum::<f64>() == 0.0
    })
}

/// True if the loops going round the period along `axis` go round it as often each way.
fn paired(loops: &[(bool, Vec<Coedge>)], surface: &Surface, axis: usize) -> bool {
    let Some(period) = closure_period(surface, axis) else { return false };
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
    let mut open: Vec<usize> = (0..loops.len()).filter(|&i| component(loop_shift(&loops[i].1), axis).abs() > 0.5 * period).collect();
    // One loop round a sphere's period (a cap with its pole inside): as FixMissingSeam, the
    // second loop made, a degenerate edge at the pole the loop's way round faces (north for
    // a loop running u up), running the period back; the two then joined as a band's are.
    if let ([i], Surface::Sphere(sphere), 0) = (&open[..], surface, axis) {
        let sign = component(loop_shift(&loops[*i].1), 0).signum();
        let start = coedge_ends(&loops[*i].1[0]).0;
        let after = start.x + sign * period;
        let v = sign * std::f64::consts::FRAC_PI_2;
        let vertex = model.add_vertex(sphere.point(after, v), crate::tolerance::LINEAR);
        let edge = model.add_edge(vertex, vertex, crate::tolerance::LINEAR)?;
        let pcurve = Curve2::BSpline(BSplineCurve2::new(vec![Point2::new(after, v), Point2::new(after - sign * period, v)], vec![0.0, 1.0], vec![2, 2], 1));
        loops.push((false, vec![(edge, true, pcurve, 0.0, 1.0)]));
        open.push(loops.len() - 1);
    }
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Line;
    use crate::geometry::surface::Torus;
    use crate::math::{Direction, Frame};

    /// A short edge (0.1 long) off a small torus by 0.011 to 0.013, more than a tenth of its
    /// length: projected all the same (as OCCT does), its tolerance to cover the gap.
    #[test]
    fn a_short_edge_well_off_its_face() {
        let torus = Surface::Torus(Torus::new(Frame::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 0.0, 0.0)), 0.755, 0.663));
        let line = Curve::Line(Line::new(Point::new(1.429, 0.0, -0.05), Direction::new(0.0, 0.0, 1.0)));
        let edge = CurveSpan { curve: &line, first: 0.0, last: 0.1 };
        let pcurve = project(&edge, &torus, false, 1e-4).unwrap();
        // (Along the tube's meridian at the outer equator: u 0, v about 0, in whichever
        // period, rising across it.)
        let (a, b) = (pcurve.point(0.0), pcurve.point(0.1));
        assert!(a.x.abs() < 1e-3 && b.x.abs() < 1e-3 && a.y.sin() < 0.0 && b.y.sin() > 0.0 && (b.y - a.y - 0.148).abs() < 0.01, "{a:?} {b:?}");
    }
}
