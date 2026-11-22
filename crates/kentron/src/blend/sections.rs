//! Fillets whose section across the edge is not a circle: a conic of a given shoulder, and the
//! curvature-continuous blend section of order 2.
//!
//! OCCT has no counterpart: its fillet sections (`ChFi3d_FilletShape`: rational, quasi-angular,
//! polynomial) are all ways of making the circle. These are built from the spec's own definition
//! of the two sections on kentron's machinery for straight edges between planes (as the
//! chamfer's): each edge's cross-section between the edge and the section curve, swept along
//! the edge, run past its ends, mitred where two edges meet, and cut away (or added at a
//! concave edge).

use super::*;
use crate::geometry::curve::BSplineCurve;

/// The section of a blend across its edge, in the cross-section at the edge's start.
#[derive(Clone, Copy, Debug)]
enum Kind {
    /// The conic arc meeting each face `setback` from the edge, tangent there, its corner on
    /// the edge and the middle weight `w` (shoulder `w/(1 + w)`).
    Conic { setback: f64, w: f64 },
    /// The blend section of order 2 between the contacts of the ball of `offset`, its apex on
    /// the edge, at the thumbweight whose section has radius `center` at its middle.
    G2 { offset: f64, center: f64 },
}

/// `shape` with `edges` rounded by a conic instead of a circle: each edge replaced by a blend
/// meeting both faces `setback` from the edge, tangent to them there, whose section across the
/// edge is the conic arc from one contact to the other with its corner on the edge and its
/// shoulder `rho` of the way from the chord's middle to the corner: the rational quadratic with
/// the contacts and the corner as poles and the middle weight `rho/(1 − rho)`. `rho = ½` is a
/// parabola; at a square edge `√2 − 1` is the quarter circle.
///
/// Supported for straight edges between planes (the section is defined there); edges carrying
/// a named edge on tangentially are blended too, and two blended edges meeting at a vertex
/// meet in a mitre, as chamfers do. Three or more meeting at a vertex, or a convex and a
/// concave one meeting, are refused. `setback` must be positive, both contacts on their faces
/// and `0 < rho < 1`.
pub fn fillet_conic(model: &mut Model, shape: Shape, edges: &[Shape], setback: f64, rho: f64) -> Result<Shape> {
    if !(setback > tolerance::LINEAR) {
        return Err(Error::InvalidInput(format!("a conic fillet's setback must be positive, got {setback}")));
    }
    if !(rho > 0.0 && rho < 1.0) {
        return Err(Error::InvalidInput(format!("a conic fillet's rho must lie strictly between 0 and 1, got {rho}")));
    }
    shaped(model, shape, edges, Kind::Conic { setback, w: rho / (1.0 - rho) }, setback)
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
    shaped(model, shape, edges, Kind::G2 { offset, center }, offset)
}

fn shaped(model: &mut Model, shape: Shape, edges: &[Shape], kind: Kind, reach: f64) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    let edges = sharp(model, shape, edges, |e| *e)?;
    on_solids(model, shape, reach, |model, shape, found| {
        let edges = remapped(model, shape, &edges, |e| *e, |_, c| Some(c), found)?;
        each_solid(model, shape, &edges, |e| *e, |_, c| Some(c), |model, shape, edges| shaped_solid_edges(model, shape, edges, kind))
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
    /// The contacts and section curve of a `kind` blend, in the cross-section at `origin`.
    fn shaped_curve(&self, kind: Kind) -> Result<Section> {
        let o = self.origin;
        Ok(match kind {
            Kind::Conic { setback, w } => Section::conic(o + self.into[0] * setback, o, o + self.into[1] * setback, w),
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
}

/// [`shaped`] on a shape with solids.
fn shaped_solid_edges(model: &mut Model, shape: Shape, edges: &[Shape], kind: Kind) -> Result<Shape> {
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    let indexed: Vec<(Shape, f64)> = edges.iter().map(|&e| (e, 0.0)).collect();
    let chains = tangent_chains(model, shape, &indexed, 0.0);
    let faces = model.unique_subshapes(shape, NodeKind::Face);
    let known = runover::KnownFaces::default();
    let mut corners: Vec<(Corner, f64)> = Vec::new();
    let mut sections = Vec::new();
    for (edge, _) in chains {
        let corner = Corner::of(model, shape, edge, true)?.filter(|c| matches!(c.path, Path::Line { .. }) && c.faces.iter().all(|f| matches!(untrimmed(&f.1), Surface::Plane(_))));
        let Some(corner) = corner else {
            return Err(Error::InvalidInput("a conic or G2 fillet is supported only along straight edges between planes".into()));
        };
        let section = corner.shaped_curve(kind)?;
        let contacts = [section.poles[0], section.poles[section.poles.len() - 1]];
        // Both contacts on their faces all along the edge.
        let Path::Line { length } = corner.path else { unreachable!() };
        let find = |node: NodeId| faces.iter().copied().find(|f| f.node == node);
        if let (Some(a), Some(b)) = (find(corner.faces[0].0), find(corner.faces[1].0)) {
            if runover::rolls_off(model, [a, b], corner.tangent, length, contacts, [corner.origin, corner.end], &known, false) {
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
}
