//! Fillets whose section across the edge is not a circle: a conic of a given shoulder (rho), the
//! elliptic section of a round with two distances, and the curvature-continuous blend section of
//! order 2.
//!
//! OCCT has no counterpart: its fillet sections (`ChFi3d_FilletShape`: rational, quasi-angular,
//! polynomial) are all ways of making the circle. These are built from the sections' own
//! definitions on kentron's machinery for edges whose faces cross the edge's cross-section in
//! straight lines (straight edges between planes, circular edges between planes and the
//! cylinders and cones about the circle's axis, as the chamfer's): each edge's cross-section
//! between the edge and the section curve, swept along the edge (extruded or revolved), run past
//! its ends, mitred where two edges meet, and cut away (or added at a concave edge).

use super::*;
use crate::geometry::curve::BSplineCurve;

/// The shape of a conic fillet's section across its edge: a conic arc from one contact to the
/// other, tangent to the faces there, with its corner (the point where the faces' tangents at
/// the two contacts meet) on the edge. As a rational quadratic Bézier its poles are the two
/// contacts and the corner; the two kinds differ in the corner's weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConicSection {
    /// The conic whose shoulder, the point where it crosses the line from the corner to the
    /// middle of the chord between the contacts, lies `rho` of the way from the chord to the
    /// corner (`0 < rho < 1`): the corner's weight is `rho/(1 − rho)`. `rho = ½` is a parabola,
    /// below an ellipse, above a hyperbola; with equal distances at a square edge, `√2 − 1` is
    /// the quarter circle.
    Rho(f64),
    /// The ellipse arc of a round with two distances: the circular round's section (the weight
    /// `sin(α/2)` of the circle's arc between faces at `α`, as the faces' traces meet in the
    /// cross-section) stretched along the faces to meet each its own distance from the edge.
    /// Between faces square to each other it is the quarter ellipse with the two distances as
    /// its semi-axes, centred where the lines through the contacts square to the faces cross;
    /// with equal distances `d` it is the rolling ball's round of radius `d·tan(α/2)`.
    Elliptic,
}

/// An edge to round with a [`ConicSection`], with how far back along each face the round meets
/// it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConicFilletEdge {
    pub edge: Shape,
    /// The distance from the edge to the contact on `face`, and to the contact on the other
    /// face, measured along the faces square to the edge.
    pub distances: (f64, f64),
    /// One of the two faces at the edge: the face of the first distance (either if `None`).
    pub face: Option<Shape>,
    pub section: ConicSection,
}

/// The section of a blend across its edge, in the cross-section at the edge's start.
#[derive(Clone, Copy, Debug)]
enum Kind {
    /// A [`ConicSection`] meeting the faces `distances` from the edge: the first along the
    /// chosen face, if one is.
    Conic { distances: [f64; 2], section: ConicSection },
    /// The blend section of order 2 between the contacts of the ball of `offset`, its apex on
    /// the edge, at the thumbweight whose section has radius `center` at its middle.
    G2 { offset: f64, center: f64 },
}

/// An edge to blend: the edge, the face its first distance is measured on, and the section.
#[derive(Clone, Copy, Debug)]
struct Wanted {
    edge: Shape,
    face: Option<Shape>,
    kind: Kind,
}

/// `shape` with `edges` rounded by a conic instead of a circle: each edge replaced by a blend
/// meeting both faces `setback` from the edge, tangent to them there, whose section across the
/// edge is the conic arc from one contact to the other with its corner on the edge and its
/// shoulder `rho` of the way from the chord's middle to the corner: the rational quadratic with
/// the contacts and the corner as poles and the middle weight `rho/(1 − rho)`. `rho = ½` is a
/// parabola; at a square edge `√2 − 1` is the quarter circle. [`fillet_conic_edges`] with
/// [`ConicSection::Rho`] and equal distances.
///
/// Supported where both faces cross the edge's cross-section in straight lines: straight edges
/// between planes, and circular edges between planes, cylinders and cones about the circle's
/// axis (the blend then a surface of revolution, the section turned about the axis). Edges
/// carrying a named edge on tangentially are blended too, and two blended edges meeting at a
/// vertex meet in a mitre, as chamfers do. Three or more meeting at a vertex, or a convex and a
/// concave one meeting, are refused. `setback` must be positive, both contacts on their faces
/// (on a circular edge, on the edge's side of its axis) and `0 < rho < 1`.
pub fn fillet_conic(model: &mut Model, shape: Shape, edges: &[Shape], setback: f64, rho: f64) -> Result<Shape> {
    let edges: Vec<ConicFilletEdge> = edges.iter().map(|&edge| ConicFilletEdge { edge, distances: (setback, setback), face: None, section: ConicSection::Rho(rho) }).collect();
    fillet_conic_edges(model, shape, &edges)
}

/// `shape` with each of `edges` rounded by its [`ConicSection`]: replaced by a blend meeting the
/// edge's two faces at its two distances from it (the first on its chosen face), tangent to them
/// there, its section across the edge the conic arc between those contacts with its corner on
/// the edge. With [`ConicSection::Elliptic`] and two distances it is an unequal round of
/// elliptic section; with [`ConicSection::Rho`], a conic blend of that shoulder.
///
/// Supported as [`fillet_conic`] is. The distances must be positive and a rho strictly between 0
/// and 1.
pub fn fillet_conic_edges(model: &mut Model, shape: Shape, edges: &[ConicFilletEdge]) -> Result<Shape> {
    let mut wanted = Vec::with_capacity(edges.len());
    for e in edges {
        let (a, b) = e.distances;
        if !(a > tolerance::LINEAR) || !(b > tolerance::LINEAR) {
            return Err(Error::InvalidInput(format!("a conic fillet's distances must be positive, got {a} and {b}")));
        }
        if let ConicSection::Rho(rho) = e.section {
            if !(rho > 0.0 && rho < 1.0) {
                return Err(Error::InvalidInput(format!("a conic fillet's rho must lie strictly between 0 and 1, got {rho}")));
            }
        }
        wanted.push(Wanted { edge: e.edge, face: e.face, kind: Kind::Conic { distances: [a, b], section: e.section } });
    }
    let reach = edges.iter().map(|e| e.distances.0.max(e.distances.1)).fold(0.0, f64::max);
    shaped(model, shape, &wanted, reach)
}

/// `shape` with `edges` rounded by a curvature-continuous blend: it meets each face where a
/// ball of `offset` rolling along the edge touches it, and its section across the edge is the
/// blend section of order 2 between the two contacts with its apex on the edge (E. Hartmann,
/// "Parametric Gⁿ blending of curves and surfaces", The Visual Computer 17, 2001),
/// `apex + α(t)(a − apex) + β(t)(b − apex)` with `α = f·(1 − p₁)`, `β = (1 − f)·p₂`,
/// `f = t̄³/(t̄³ + t³)`, `p₁ = t/(1 − w + wt)`, `p₂ = t(1 − w)/(1 − wt)`, at the thumbweight `w`
/// below 1 that gives the section's middle the radius of curvature `center`. That curve is
/// exactly a rational Bézier of degree 5 (its numerators and its denominator
/// `(t̄³ + t³)(1 − w + wt)(1 − wt)` are quintics with positive Bernstein coefficients), so the
/// blend is exact. Its first three poles lie on the leg: the blend starts as straight as the
/// face it leaves.
///
/// Supported as [`fillet_conic`] is. `offset` and `center` must be positive, the ball of
/// `offset` must touch both faces on them, and `center` must be a radius the section takes at
/// some thumbweight below 1.
pub fn fillet_g2(model: &mut Model, shape: Shape, edges: &[Shape], offset: f64, center: f64) -> Result<Shape> {
    if !(offset > tolerance::LINEAR) || !(center > tolerance::LINEAR) {
        return Err(Error::InvalidInput(format!("a G2 fillet's offset and center radius must be positive, got {offset} and {center}")));
    }
    let wanted: Vec<Wanted> = edges.iter().map(|&edge| Wanted { edge, face: None, kind: Kind::G2 { offset, center } }).collect();
    shaped(model, shape, &wanted, offset)
}

fn shaped(model: &mut Model, shape: Shape, edges: &[Wanted], reach: f64) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    let edges = sharp(model, shape, edges, |e| e.edge)?;
    // (An edge whose distances are measured from a chosen face is not carried onto a copy,
    // where that face is not.)
    let with_edge = |e: &Wanted, c: Shape| e.face.is_none().then_some(Wanted { edge: c, ..*e });
    on_solids(model, shape, reach, |model, shape, found| {
        let edges = remapped(model, shape, &edges, |e| e.edge, with_edge, found)?;
        each_solid(model, shape, &edges, |e| e.edge, with_edge, shaped_solid_edges)
    })
}

/// A section curve: a rational Bézier from one contact to the other.
struct Section {
    poles: Vec<Point>,
    weights: Vec<f64>,
}

impl Section {
    fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    /// The conic arc from `a` to `b` with its corner at `apex` and middle weight `w`.
    fn conic(a: Point, apex: Point, b: Point, w: f64) -> Section {
        Section { poles: vec![a, apex, b], weights: vec![1.0, w, 1.0] }
    }

    /// The order-2 blend section of thumbweight `w` from `a` to `b` on the legs from `apex`, as
    /// the rational quintic it is: with `A = a − apex`, `B = b − apex`, the point times the
    /// denominator `D = (t̄³ + t³)(1 − w + wt)(1 − wt)` is
    /// `apex·D + (1 − w)t̄⁴(1 − wt)·A + (1 − w)t⁴(1 − w + wt)·B`.
    fn g2(a: Point, apex: Point, b: Point, w: f64) -> Section {
        let k = 1.0 - w;
        let (d1, d2) = ([k, 1.0], [1.0, k]);
        let x: Vec<f64> = bernstein_product(&[1.0, 0.0, 0.0, 0.0, 0.0], &d2).into_iter().map(|c| k * c).collect();
        let y: Vec<f64> = bernstein_product(&[0.0, 0.0, 0.0, 0.0, 1.0], &d1).into_iter().map(|c| k * c).collect();
        let d = bernstein_product(&bernstein_product(&[1.0, 0.0, 0.0, 1.0], &d1), &d2);
        let (legs_a, legs_b) = (a - apex, b - apex);
        let poles = (0..6).map(|i| apex + legs_a * (x[i] / d[i]) + legs_b * (y[i] / d[i])).collect();
        Section { poles, weights: d.into_iter().map(|c| c / k).collect() }
    }

    /// The radius of curvature at `t`, from the homogeneous curve's derivatives.
    fn radius_at(&self, t: f64) -> f64 {
        let n = self.degree();
        let homogeneous: Vec<[f64; 4]> = self.poles.iter().zip(&self.weights).map(|(p, &w)| [p.x * w, p.y * w, p.z * w, w]).collect();
        let at = |c: &[[f64; 4]]| -> [f64; 4] { std::array::from_fn(|j| bernstein(&c.iter().map(|q| q[j]).collect::<Vec<_>>(), t)) };
        let hodograph = |c: &[[f64; 4]]| -> Vec<[f64; 4]> { c.windows(2).map(|q| std::array::from_fn(|j| (c.len() - 1) as f64 * (q[1][j] - q[0][j]))).collect() };
        let first = hodograph(&homogeneous);
        let second = if n >= 2 { hodograph(&first) } else { vec![[0.0; 4]] };
        let (h0, h1, h2) = (at(&homogeneous), at(&first), at(&second));
        let v = |h: [f64; 4]| Vector::new(h[0], h[1], h[2]);
        let c = v(h0) * (1.0 / h0[3]);
        let d1 = (v(h1) - c * h1[3]) * (1.0 / h0[3]);
        let d2 = (v(h2) - d1 * (2.0 * h1[3]) - c * h2[3]) * (1.0 / h0[3]);
        d1.length().powi(3) / d1.cross(d2).length()
    }

    fn edge(&self, model: &mut Model) -> Result<Shape> {
        let degree = self.degree();
        let mut curve = BSplineCurve::try_new(self.poles.clone(), vec![0.0, 1.0], vec![degree + 1, degree + 1], degree)?.with_weights(self.weights.clone());
        // (The quintic in eight spans, the same curve: its speed changes much along it, most
        // at the ends, and the faces' integrals, taken span by span, would miss by 1e-4 of
        // the blend's area over one span.)
        if degree > 2 {
            curve.insert_knots(&(1..8).map(|k| (k as f64 / 8.0, 1)).collect::<Vec<_>>());
        }
        make::edge_from_curve(model, Curve::BSpline(curve))
    }
}

/// The corner weight of the circle's arc tangent to two lines from a corner at equal distances
/// from it, the lines leaving it along the unit directions `into`: `cos` of half the arc's
/// turn, `π − α` between lines meeting at `α`, so `sin(α/2)`, half the chord between `into`.
fn elliptic_weight(into: [Vector; 2]) -> f64 {
    0.5 * (into[0] - into[1]).length()
}

/// The value at `t` of the polynomial with Bernstein coefficients `c` (de Casteljau).
fn bernstein(c: &[f64], t: f64) -> f64 {
    let mut b = c.to_vec();
    for r in 1..b.len() {
        for i in 0..b.len() - r {
            b[i] = (1.0 - t) * b[i] + t * b[i + 1];
        }
    }
    b[0]
}

/// The Bernstein coefficients of the product of two polynomials given by theirs.
fn bernstein_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let binomial = |n: usize, k: usize| (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64);
    let (m, n) = (a.len() - 1, b.len() - 1);
    (0..=m + n)
        .map(|k| {
            let mut sum = 0.0;
            for i in k.saturating_sub(n)..=k.min(m) {
                sum += binomial(m, i) * binomial(n, k - i) * a[i] * b[k - i];
            }
            sum / binomial(m + n, k)
        })
        .collect()
}

/// The thumbweight (below 1) at which the order-2 section on the legs from `apex` to `a` and
/// `b` has radius of curvature `center` at its middle (`t = ½`, by symmetry). The middle radius
/// falls as `w` rises, from unbounded to zero at 1: bisected on a bracket stepped down from 0.
/// (The section is sound for every `w` below 1: its weights stay positive. The equal-area
/// thumbweight, which the spec counts among the sections, lies below 0, about −0.525 at a
/// square edge, so the spec's "between 0 and 1" cannot be meant as a bound on `w`.)
fn g2_thumbweight(a: Point, apex: Point, b: Point, center: f64) -> Result<f64> {
    let radius = |w: f64| Section::g2(a, apex, b, w).radius_at(0.5);
    let (mut lo, mut hi) = (0.0f64, 1.0 - 1e-12);
    while !(radius(lo) > center) {
        if !(lo > -1e9) || !radius(lo).is_finite() {
            return Err(Error::InvalidInput(format!("a G2 fillet's center radius {center} is no radius its section takes")));
        }
        hi = lo;
        lo = 2.0 * lo - 1.0;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if radius(mid) > center {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 1e-15 * (1.0 - lo) {
            break;
        }
    }
    Ok(0.5 * (lo + hi))
}

impl Corner {
    /// The contacts and section curve of a `kind` blend, in the cross-section at `origin`: a
    /// conic's first distance along the second face if `first_on_second`.
    fn shaped_curve(&self, kind: Kind, first_on_second: bool) -> Result<Section> {
        let o = self.origin;
        Ok(match kind {
            Kind::Conic { distances: [a, b], section } => {
                let [d0, d1] = if first_on_second { [b, a] } else { [a, b] };
                let w = match section {
                    ConicSection::Rho(rho) => rho / (1.0 - rho),
                    ConicSection::Elliptic => elliptic_weight(self.into),
                };
                Section::conic(o + self.into[0] * d0, o, o + self.into[1] * d1, w)
            }
            Kind::G2 { offset, center } => {
                let (_, t1, t2) = self.ball(offset)?;
                Section::g2(t1, o, t2, g2_thumbweight(t1, o, t2, center)?)
            }
        })
    }


    /// The cross-section of the tool: the faces' traces from the edge to the contacts, and the
    /// section curve between them.
    fn shaped_section(&self, model: &mut Model, section: &Section) -> Result<Shape> {
        let o = self.origin;
        let (t1, t2) = (section.poles[0], section.poles[section.poles.len() - 1]);
        let leg1 = make::edge_between(model, o, t1)?;
        let curve = section.edge(model)?;
        let leg2 = make::edge_between(model, t2, o)?;
        let wire = make::wire(model, &[leg1, curve, leg2])?;
        make::planar_face(model, wire)
    }

    /// Whether a blend touching the faces `faces` at `contacts` in the cross-section at `origin`
    /// touches them off their bounds somewhere along the edge (on a circular edge, past its
    /// axis anywhere).
    fn shaped_rolls_off(&self, model: &Model, faces: [Shape; 2], contacts: [Point; 2], known: &runover::KnownFaces) -> Result<bool> {
        let ends = [self.origin, self.end];
        Ok(match self.path {
            Path::Line { length } => runover::rolls_off(model, faces, self.tangent, length, contacts, ends, known, false),
            Path::Arc { axis, radius, angle } => {
                let z = Vector::from(axis.direction);
                let foot = axis.origin + z * (self.origin - axis.origin).dot(z);
                let outward = (self.origin - foot).normalized();
                if contacts.iter().any(|t| !((*t - foot).dot(outward) > tolerance::LINEAR)) {
                    return Err(Error::InvalidInput("a conic or G2 fillet's contact lies past its circular edge's axis".into()));
                }
                const SAMPLES: usize = 32;
                // (Not at the very ends of an open arc, where the contact meets the faces there.)
                let (from, span) = if angle >= TAU - tolerance::ANGULAR { (0.0, 1.0) } else { (0.02, 0.96) };
                faces.iter().zip(contacts).any(|(&face, contact)| {
                    (0..=SAMPLES).any(|i| {
                        let turn = Transform::rotation(axis, angle * (from + span * i as f64 / SAMPLES as f64));
                        runover::off_face(model, face, contact.transformed(&turn), radius * angle, ends, known, false)
                    })
                })
            }
        })
    }
}

/// [`shaped`] on a shape with solids.
fn shaped_solid_edges(model: &mut Model, shape: Shape, edges: &[Wanted]) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    // The edges continuing them tangentially (indices into `edges` as their "radius").
    let indexed: Vec<(Shape, f64)> = edges.iter().enumerate().map(|(i, e)| (e.edge, i as f64)).collect();
    let chains = tangent_chains(model, shape, &indexed, 0.0);
    let chosen: Vec<Option<Shape>> = edges.iter().map(|e| e.face).collect();
    let firsts = chain_faces(model, shape, &chains, &chosen);
    let faces = model.unique_subshapes(shape, NodeKind::Face);
    let known = runover::KnownFaces::default();
    let mut corners: Vec<(Corner, f64)> = Vec::new();
    let mut sections = Vec::new();
    for (k, (edge, index)) in chains.into_iter().enumerate() {
        let wanted = edges[index as usize];
        let corner = Corner::of(model, shape, edge, true)?.filter(|c| c.traces.iter().all(|t| matches!(t, Trace::Line)));
        let Some(corner) = corner else {
            return Err(Error::InvalidInput(
                "a conic or G2 fillet is supported only where both faces cross the edge's cross-section in straight lines: straight edges between planes, circular edges between planes, cylinders and cones about the circle's axis".into(),
            ));
        };
        // (The first distance along the chosen face: for an edge carrying the chosen one on,
        // the face on the same side.)
        let first_on_second = firsts[k].is_some_and(|f| corner.faces[1].0 == f);
        let section = corner.shaped_curve(wanted.kind, first_on_second)?;
        let contacts = [section.poles[0], section.poles[section.poles.len() - 1]];
        // Both contacts on their faces all along the edge.
        let find = |node: NodeId| faces.iter().copied().find(|f| f.node == node);
        if let (Some(a), Some(b)) = (find(corner.faces[0].0), find(corner.faces[1].0)) {
            if corner.shaped_rolls_off(model, [a, b], contacts, &known)? {
                return Err(Error::InvalidInput("a fillet's contact lies off its face".into()));
            }
        }
        let size = contacts.iter().map(|&c| c.distance(corner.origin)).fold(0.0, f64::max);
        corners.push((corner, size));
        sections.push(section);
    }
    // Where blended edges meet: two in a mitre; more, or a convex and a concave one, refused.
    let mut meetings: Vec<(Point, Vec<usize>)> = Vec::new();
    for (i, (corner, _)) in corners.iter().enumerate() {
        for p in corner.ends() {
            match meetings.iter_mut().find(|(q, _)| q.distance(p) <= tolerance::LINEAR) {
                Some((_, at)) => at.push(i),
                None => meetings.push((p, vec![i])),
            }
        }
    }
    for (_, at) in &meetings {
        let tangent = |i: usize, j: usize| corners[i].0.tangent.cross(corners[j].0.tangent).length() <= tolerance::ANGULAR;
        match at[..] {
            [_] => {}
            [i, j] if tangent(i, j) => {}
            [i, j] if corners[i].0.convex == corners[j].0.convex => {}
            [_, _] => return Err(Error::InvalidInput("a conic or G2 fillet's convex and concave edges meeting at a vertex".into())),
            _ => return Err(Error::InvalidInput("a conic or G2 fillet's corner of three edges or more".into())),
        }
    }
    let mitred = mitres(&corners);
    let extended = extensions(model, shape, &corners);
    let overhangs = overhanging_ends(model, shape, &corners, &extended, &mitred);
    let mut tools = Vec::with_capacity(corners.len());
    for (i, (corner, size)) in corners.iter().enumerate() {
        let section = corner.shaped_section(model, &sections[i])?;
        let tool = corner.sweep_section(model, section, mitre_runs(corner, *size, extended[i], &mitred[i]))?;
        let tool = mitred_tool(model, tool, &mitred[i])?;
        let tool = mitred_tool(model, tool, &overhangs[i])?;
        record_tool(model, tool, &corner.edges);
        tools.push((tool, corner.convex));
    }
    let blended = with_tools(model, shape, &tools)?;
    // (The faces the tools' square ends flush with a face at their edge's end split, one again.)
    let merged = crate::unify::unify_same_domain_with(model, blended, crate::unify::UnifyOptions { faces: true, edges: false, same_surface_only: false })?;
    crate::unify::unify_same_domain_with(model, merged, crate::unify::UnifyOptions { faces: false, edges: true, same_surface_only: true })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::ParametricCurve;

    /// The conic of `rho = √2 − 1` at a square corner is the quarter circle: its radius of
    /// curvature is the setback all along.
    #[test]
    fn a_conic_at_a_square_corner_with_rho_root_two_less_one_is_a_quarter_circle() {
        let rho = 2f64.sqrt() - 1.0;
        let s = Section::conic(Point::new(2.0, 0.0, 0.0), Point::ORIGIN, Point::new(0.0, 2.0, 0.0), rho / (1.0 - rho));
        for t in [0.1, 0.3, 0.5, 0.8] {
            assert!((s.radius_at(t) - 2.0).abs() < 1e-9, "{t}: {}", s.radius_at(t));
        }
    }

    /// The order-2 section's poles: the point times its denominator is as the definition has
    /// it, and the curve starts as straight as each leg.
    #[test]
    fn the_g2_section_is_its_definition_and_flat_at_its_ends() {
        let (a, apex, b) = (Point::new(3.0, 0.0, 0.0), Point::new(0.0, 0.0, 1.0), Point::new(0.0, 2.0, 0.0));
        for w in [0.05, 0.4, 0.9] {
            let s = Section::g2(a, apex, b, w);
            let curve = BSplineCurve::new(s.poles.clone(), vec![0.0, 1.0], vec![6, 6], 5).with_weights(s.weights.clone());
            for t in [0.0f64, 0.2, 0.5, 0.7, 1.0] {
                let u = 1.0 - t;
                let f = u.powi(3) / (u.powi(3) + t.powi(3));
                let (p1, p2) = (t / (1.0 - w + w * t), t * (1.0 - w) / (1.0 - w * t));
                let want = apex + (a - apex) * (f * (1.0 - p1)) + (b - apex) * ((1.0 - f) * p2);
                assert!(Curve::BSpline(curve.clone()).point(t).distance(want) < 1e-12, "{w} {t}");
            }
            assert!(s.radius_at(1e-4) > 1e3 && s.radius_at(1.0 - 1e-4) > 1e3);
        }
    }

    /// The parabola (`rho = ½`) between a square box edge's faces: a parabolic segment is two
    /// thirds of the triangle on its chord and tangents, so the cut takes the remaining third
    /// of the triangle, `d²/6`, per unit length.
    #[test]
    fn a_parabolic_round_of_a_box_edge_takes_a_sixth_of_the_setback_squared() {
        let mut model = Model::new();
        let solid = crate::primitives::Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let edge = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                p.x.abs() < 1e-9 && p.y.abs() < 1e-9 && q.x.abs() < 1e-9 && q.y.abs() < 1e-9
            })
            .unwrap();
        let d = 2.0;
        let rounded = fillet_conic(&mut model, solid, &[edge], d, 0.5).unwrap();
        let v = crate::properties::volume(&model, rounded).mass;
        assert!((v - (6000.0 - d * d / 6.0 * 30.0)).abs() < 1e-6, "{v}");
        assert!(crate::validity::check(&model, rounded).is_valid());
        assert!(fillet_conic(&mut model, solid, &[edge], d, 1.0).is_err());
        assert!(fillet_conic(&mut model, solid, &[edge], 0.0, 0.5).is_err());
        assert!(fillet_conic(&mut model, solid, &[edge], 11.0, 0.5).is_err());
    }

    /// A G2 round of a box edge takes off the area its section cuts from the triangle on its
    /// legs, `½∫(αβ' − α'β)` with legs of the offset (at a square edge the ball's contacts lie
    /// `offset` from the edge), here by the midpoint rule on the definition; its thumbweight
    /// gives the middle the asked radius of curvature.
    #[test]
    fn a_g2_round_of_a_box_edge_takes_off_its_sections_area() {
        let mut model = Model::new();
        let solid = crate::primitives::Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let edge = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                p.x.abs() < 1e-9 && p.y.abs() < 1e-9 && q.x.abs() < 1e-9 && q.y.abs() < 1e-9
            })
            .unwrap();
        let (offset, center) = (2.0, 1.2);
        let (a, o, b) = (Point::new(offset, 0.0, 0.0), Point::ORIGIN, Point::new(0.0, offset, 0.0));
        let w = g2_thumbweight(a, o, b, center).unwrap();
        assert!((Section::g2(a, o, b, w).radius_at(0.5) - center).abs() < 1e-9);
        let weights = |t: f64| {
            let u = 1.0 - t;
            let f = u.powi(3) / (u.powi(3) + t.powi(3));
            (f * (1.0 - t / (1.0 - w + w * t)), (1.0 - f) * t * (1.0 - w) / (1.0 - w * t))
        };
        let n = 20000;
        let (h, mut share) = (1e-7, 0.0);
        for k in 0..n {
            let t = (k as f64 + 0.5) / n as f64;
            let ((x, y), (x0, y0), (x1, y1)) = (weights(t), weights(t - h), weights(t + h));
            share += (x * (y1 - y0) - (x1 - x0) * y) / (2.0 * h) / n as f64;
        }
        let rounded = fillet_g2(&mut model, solid, &[edge], offset, center).unwrap();
        let removed = 6000.0 - crate::properties::volume(&model, rounded).mass;
        assert!((removed - 30.0 * 0.5 * share * offset * offset).abs() < 1e-6, "{removed} {share}");
        assert!(crate::validity::check(&model, rounded).is_valid());
        assert!(fillet_g2(&mut model, solid, &[edge], 11.0, center).is_err());
        assert!(fillet_g2(&mut model, solid, &[edge], offset, 0.0).is_err());
    }

    /// The edge of a 10 × 20 × 30 box along z through the origin, and the box.
    fn box_and_edge(model: &mut Model) -> (Shape, Shape) {
        let solid = crate::primitives::Cuboid::new(10.0, 20.0, 30.0).build(model).unwrap();
        let edge = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (a, b) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(b));
                p.x.abs() < 1e-9 && p.y.abs() < 1e-9 && q.x.abs() < 1e-9 && q.y.abs() < 1e-9
            })
            .unwrap();
        (solid, edge)
    }

    /// The face of `solid` on the plane through `p` square to `normal`.
    fn face_on(model: &Model, solid: Shape, p: Point, normal: Vector) -> Shape {
        model
            .unique_subshapes(solid, NodeKind::Face)
            .into_iter()
            .find(|&f| match model.surface(f) {
                Some(Surface::Plane(plane)) => {
                    let n = Vector::from(plane.normal());
                    n.cross(normal).length() < 1e-9 && (p - plane.origin()).dot(n).abs() < 1e-9
                }
                _ => false,
            })
            .unwrap()
    }

    /// A conic's shoulder, its point at `t = ½`, lies `rho` of the way from the chord's middle
    /// to the corner; at `rho = ½` it is the parabola, the polynomial quadratic Bézier.
    #[test]
    fn a_conic_s_shoulder_lies_rho_of_the_way_to_its_corner_and_rho_one_half_is_the_parabola() {
        let (a, o, b) = (Point::new(3.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new(0.5, 2.0, 0.0));
        let chord_middle = a + (b - a) * 0.5;
        for rho in [0.2, 2f64.sqrt() - 1.0, 0.5, 0.8] {
            let s = Section::conic(a, o, b, rho / (1.0 - rho));
            let curve = Curve::BSpline(BSplineCurve::new(s.poles.clone(), vec![0.0, 1.0], vec![3, 3], 2).with_weights(s.weights.clone()));
            assert!(curve.point(0.5).distance(chord_middle + (o - chord_middle) * rho) < 1e-12, "{rho}");
            if rho == 0.5 {
                for t in [0.0f64, 0.15, 0.4, 0.75, 1.0] {
                    let u = 1.0 - t;
                    let want = Point::ORIGIN + (a - Point::ORIGIN) * (u * u) + (o - Point::ORIGIN) * (2.0 * t * u) + (b - Point::ORIGIN) * (t * t);
                    assert!(curve.point(t).distance(want) < 1e-12, "{t}");
                }
            }
        }
    }

    /// With equal legs, the elliptic section's weight makes the circle's arc tangent to both
    /// legs, of radius `d·tan(α/2)` between legs at `α`.
    #[test]
    fn the_elliptic_section_with_equal_distances_is_the_circle() {
        let d = 2.0;
        for degrees in [35.0f64, 60.0, 90.0, 125.0] {
            let alpha = degrees.to_radians();
            let (u, v) = (Vector::new(1.0, 0.0, 0.0), Vector::new(alpha.cos(), alpha.sin(), 0.0));
            let s = Section::conic(Point::ORIGIN + u * d, Point::ORIGIN, Point::ORIGIN + v * d, elliptic_weight([u, v]));
            for t in [0.0, 0.2, 0.5, 0.9] {
                assert!((s.radius_at(t) - d * (0.5 * alpha).tan()).abs() < 1e-9, "{degrees} {t}: {}", s.radius_at(t));
            }
        }
    }

    /// An elliptic round of a square edge with equal distances is the rolling ball's round of
    /// that radius, measure for measure.
    #[test]
    fn an_elliptic_round_with_equal_distances_of_a_box_edge_is_its_circular_fillet() {
        let mut model = Model::new();
        let (solid, edge) = box_and_edge(&mut model);
        let round = ConicFilletEdge { edge, distances: (2.0, 2.0), face: None, section: ConicSection::Elliptic };
        let elliptic = fillet_conic_edges(&mut model, solid, &[round]).unwrap();
        let circular = super::super::fillet(&mut model, solid, &[edge], 2.0).unwrap();
        let measures = |s: Shape| (crate::properties::volume(&model, s).mass, crate::properties::area(&model, s).mass);
        let ((v, a), (v0, a0)) = (measures(elliptic), measures(circular));
        assert!((v - v0).abs() < 1e-9 * v0 && (a - a0).abs() < 1e-9 * a0, "{v} {v0} {a} {a0}");
        assert!(crate::validity::check(&model, elliptic).is_valid());
    }

    /// An elliptic round of a box edge, 1 along one face and 3 along the other, takes off the
    /// corner outside the quarter ellipse, `ab(1 − π/4)` per unit length, each distance on its
    /// own face.
    #[test]
    fn an_elliptic_round_of_a_box_edge_takes_the_corner_outside_its_quarter_ellipse() {
        let mut model = Model::new();
        let (solid, edge) = box_and_edge(&mut model);
        let side_x = face_on(&model, solid, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0));
        let side_y = face_on(&model, solid, Point::ORIGIN, Vector::new(0.0, 1.0, 0.0));
        // (A point near the face x = 0, 2 from the edge: within a round of 3 along that face,
        // not one of 1.)
        let probe = Point::new(0.02, 2.0, 15.0);
        for (face, removed_probe) in [(side_x, true), (side_y, false)] {
            let round = ConicFilletEdge { edge, distances: (3.0, 1.0), face: Some(face), section: ConicSection::Elliptic };
            let rounded = fillet_conic_edges(&mut model, solid, &[round]).unwrap();
            let v = crate::properties::volume(&model, rounded).mass;
            assert!((v - (6000.0 - 30.0 * 3.0 * (1.0 - std::f64::consts::FRAC_PI_4))).abs() < 1e-6, "{v}");
            assert!(crate::validity::check(&model, rounded).is_valid());
            let at = crate::boolean::classify_point(&model, rounded, probe).unwrap();
            assert_eq!(at == crate::boolean::PointLocation::Outside, removed_probe, "{at:?}");
        }
    }

    /// Round a circular edge: the rim of a cylinder of radius 5 and height 10, by a parabola 2
    /// from the edge on both faces and by a quarter ellipse 1 across the top and 3 down the
    /// side, each the section turned about the axis; volumes by Pappus, areas by the sections'
    /// arc lengths times the circles they run round (midpoint rule).
    #[test]
    fn conic_rounds_of_a_cylinder_s_rim_are_their_sections_turned_about_its_axis() {
        use std::f64::consts::PI;
        let mut model = Model::new();
        let solid = crate::primitives::Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
        let rim = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .find(|&e| matches!(model.edge_curve(e), Some((Curve::Circle(_), _, _))) && model.point(model.edge_vertices(e).0).z > 5.0)
            .unwrap();
        let top = face_on(&model, solid, Point::new(0.0, 0.0, 10.0), Vector::new(0.0, 0.0, 1.0));
        // ∫ ρ |c'| over the section `c(t) = (ρ, z)`, t in 0..1.
        let swept = |c: &dyn Fn(f64) -> (f64, f64)| {
            let (n, h) = (200000, 1e-7);
            (0..n)
                .map(|k| {
                    let t = (k as f64 + 0.5) / n as f64;
                    let ((r, _), (r0, z0), (r1, z1)) = (c(t), c(t - h), c(t + h));
                    r * ((r1 - r0).powi(2) + (z1 - z0).powi(2)).sqrt() / (2.0 * h) / n as f64
                })
                .sum::<f64>()
        };
        let parabola = fillet_conic(&mut model, solid, &[rim], 2.0, 0.5).unwrap();
        let (v, a) = (crate::properties::volume(&model, parabola).mass, crate::properties::area(&model, parabola).mass);
        let want = 114.0 * PI + 2.0 * PI * swept(&|t| (5.0 - 2.0 * (1.0 - t).powi(2), 10.0 - 2.0 * t * t));
        assert!((v - (250.0 * PI - 92.0 * PI / 15.0)).abs() < 1e-8 * v, "{v}");
        assert!((a - want).abs() < 1e-8 * a, "{a} {want}");
        assert!(crate::validity::check(&model, parabola).is_valid());
        let round = ConicFilletEdge { edge: rim, distances: (1.0, 3.0), face: Some(top), section: ConicSection::Elliptic };
        let elliptic = fillet_conic_edges(&mut model, solid, &[round]).unwrap();
        let (v, a) = (crate::properties::volume(&model, elliptic).mass, crate::properties::area(&model, elliptic).mass);
        let (ea, eb) = (1.0, 3.0);
        let want_v = 250.0 * PI - 2.0 * PI * (5.0 * ea * eb * (1.0 - PI / 4.0) - ea * ea * eb * (5.0 / 6.0 - PI / 4.0));
        let want_a = 111.0 * PI + 2.0 * PI * swept(&|t| {
            let theta = 0.5 * PI * t;
            (4.0 + theta.cos(), 7.0 + 3.0 * theta.sin())
        });
        assert!((v - want_v).abs() < 1e-8 * v, "{v} {want_v}");
        assert!((a - want_a).abs() < 1e-8 * a, "{a} {want_a}");
        assert!(crate::validity::check(&model, elliptic).is_valid());
        // Refused: a contact past the axis, a rho out of range, a distance not positive.
        assert!(fillet_conic(&mut model, solid, &[rim], 6.0, 0.5).is_err());
        assert!(fillet_conic(&mut model, solid, &[rim], 2.0, 0.0).is_err());
        let round = ConicFilletEdge { edge: rim, distances: (1.0, -3.0), face: Some(top), section: ConicSection::Elliptic };
        assert!(fillet_conic_edges(&mut model, solid, &[round]).is_err());
    }

    /// A parabolic blend at a boss's root, a concave circular edge, adds the region between the
    /// corner and the parabola, `d²/6` with moment `d³/30` about the boss's side, turned about
    /// the axis (Pappus).
    #[test]
    fn a_parabolic_blend_at_a_boss_s_root_adds_its_section_turned_about_the_axis() {
        use crate::math::Direction;
        use std::f64::consts::PI;
        let mut model = Model::new();
        let plate = crate::primitives::Cylinder::new(10.0, 2.0).build(&mut model).unwrap();
        let frame = crate::math::Frame::new(Point::new(0.0, 0.0, 2.0), Direction::Z, Direction::X);
        let boss = crate::primitives::Cylinder::new(3.0, 5.0).with_frame(frame).build(&mut model).unwrap();
        let solid = crate::boolean::union(&mut model, plate, boss).unwrap();
        let roots: Vec<Shape> = model
            .unique_subshapes(solid, NodeKind::Edge)
            .into_iter()
            .filter(|&e| matches!(model.edge_curve(e), Some((Curve::Circle(c), _, _)) if (c.radius() - 3.0).abs() < 1e-9 && (c.center().z - 2.0).abs() < 1e-9))
            .collect();
        let d = 1.0;
        let blended = fillet_conic(&mut model, solid, &roots, d, 0.5).unwrap();
        let v = crate::properties::volume(&model, blended).mass;
        let want = 200.0 * PI + 45.0 * PI + 2.0 * PI * (3.0 * d * d / 6.0 + d * d * d / 30.0);
        assert!((v - want).abs() < 1e-8 * want, "{v} {want}");
        assert!(crate::validity::check(&model, blended).is_valid());
    }
}
