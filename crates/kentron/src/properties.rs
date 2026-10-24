//! Global properties of shapes: length, area, volume and their centers.

use crate::geometry::curve::Line;
use crate::geometry::curve2::Curve2;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::geometry::surface::{Plane, Surface};
use crate::math::{BoundingBox, Direction, Matrix3, Point, Transformable, Vector};
use crate::topology::{FaceId, Model, NodeKind, Orientation, Shape};

/// The measure of a shape and its centroid, assuming unit density.
#[doc(alias = "GProp_GProps")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// Length, area or volume, depending on the query.
    pub mass: f64,
    pub center_of_mass: Point,
    /// Matrix of inertia about the center of mass.
    pub inertia: Matrix3,
}

impl MassProperties {
    /// The principal moments of inertia (eigenvalues of `inertia`), in ascending order.
    #[doc(alias = "GProp_PrincipalProps::Moments")]
    pub fn principal_moments(&self) -> [f64; 3] {
        symmetric_eigen(&self.inertia).0
    }
}

/// How properties are computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Options {
    /// Count sub-shapes shared by several parents only once.
    pub skip_shared: bool,
    /// Integrate over the stored meshes instead of the exact geometry.
    pub use_meshes: bool,
}

/// Total length of the edges of `shape`. An edge counts once per use (twice for an
/// edge between two faces) unless `Options::skip_shared` is set.
#[doc(alias = "BRepGProp::LinearProperties")]
pub fn length(model: &Model, shape: Shape) -> MassProperties {
    length_with(model, shape, Options::default())
}

/// Total area of the faces of `shape`, counting each use of a face (see [`length`]).
#[doc(alias = "BRepGProp::SurfaceProperties")]
pub fn area(model: &Model, shape: Shape) -> MassProperties {
    area_with(model, shape, Options::default())
}

/// Volume enclosed by the solids of `shape`.
#[doc(alias = "BRepGProp::VolumeProperties")]
pub fn volume(model: &Model, shape: Shape) -> MassProperties {
    volume_with(model, shape, Options::default())
}

pub fn length_with(model: &Model, shape: Shape, options: Options) -> MassProperties {
    use crate::geometry::curve::Curve;
    use crate::geometry::ParametricCurve;
    if options.use_meshes {
        todo!("length from meshes");
    }
    let edges: Vec<Shape> = match options.skip_shared {
        true => model.unique_subshapes(shape, NodeKind::Edge),
        false => model.explore(shape, NodeKind::Edge).collect(),
    };
    let mut moments = Moments::default();
    for edge in edges {
        let Some((curve, first, last)) = model.edge_curve(edge) else { continue };
        // Pieces each Gauss rule integrates: between a B-spline's knots, short in angle on
        // a conic.
        let mut breaks = vec![first];
        match &curve {
            Curve::Line(_) => {}
            // Each span in pieces: the speed |C'| is not a polynomial.
            Curve::BSpline(b) => {
                let mut knots: Vec<f64> = vec![first];
                knots.extend(b.knots().iter().copied().filter(|&k| k > first && k < last));
                knots.push(last);
                // (Eight pieces a span on a curve of few, fewer on one of many: at least 64
                // in all, as many as its spans where they are more.)
                let per_span = (64 / (knots.len() - 1)).clamp(1, 8);
                for w in knots.windows(2) {
                    breaks.extend((1..=per_span).map(|i| w[0] + (w[1] - w[0]) * i as f64 / per_span as f64));
                }
                breaks.pop();
            }
            Curve::Circle(_) | Curve::Ellipse(_) => {
                let n = ((last - first).abs() / MAX_ANGLE_STEP).ceil().max(1.0) as usize;
                breaks.extend((1..n).map(|i| first + (last - first) * i as f64 / n as f64));
            }
            _ => breaks.extend((1..16).map(|i| first + (last - first) * i as f64 / 16.0)),
        }
        breaks.push(last);
        for w in breaks.windows(2) {
            let (a, b) = (w[0], w[1]);
            for (node, weight) in GAUSS_LEGENDRE {
                let t = a + (b - a) * (node + 1.0) / 2.0;
                let (p, d) = curve.d1(t);
                moments.add(p, d.length() * weight * (b - a).abs() / 2.0, 1.0, 1.0);
            }
        }
    }
    moments.properties()
}

pub fn area_with(model: &Model, shape: Shape, options: Options) -> MassProperties {
    if options.use_meshes {
        todo!("area from meshes");
    }
    let mut total = Moments::default();
    let mut plane: Option<Option<crate::geometry::surface::Plane>> = None;
    for face in faces(model, shape, options) {
        let moments = face_area_moments(model, face);
        // The sign follows the boundary's direction of travel; an area is positive.
        total.merge_area(&moments, moments.measure.signum());
        // (Whether all the faces lie in one plane.)
        let own = match model.surface(face) {
            Some(crate::geometry::surface::Surface::Plane(p)) => Some(p.transformed(&face.placement)),
            _ => None,
        };
        plane = Some(match (plane, own) {
            (None, own) => own,
            (Some(Some(p)), Some(q)) if same_plane(&p, &q) => Some(p),
            _ => None,
        });
    }
    let mut properties = total.properties();
    // (The centre of faces in one plane lies in it: exactly, not as near as the sum comes,
    // so that faces in one plane sort by their centres' height as they stand.)
    if let Some(Some(p)) = plane {
        let n = crate::math::Vector::from(p.normal());
        let c = properties.center_of_mass;
        properties.center_of_mass = c - n * (c - p.origin()).dot(n);
    }
    properties
}

/// Whether two planes are the same, either way round.
fn same_plane(a: &crate::geometry::surface::Plane, b: &crate::geometry::surface::Plane) -> bool {
    let (n, m) = (crate::math::Vector::from(a.normal()), crate::math::Vector::from(b.normal()));
    n.cross(m).length() <= 1e-12 && (b.origin() - a.origin()).dot(n).abs() <= 1e-12 * (1.0 + (a.origin() - Point::ORIGIN).length())
}

pub fn volume_with(model: &Model, shape: Shape, options: Options) -> MassProperties {
    if options.use_meshes {
        todo!("volume from meshes");
    }
    // Sum of the cones from a fixed apex to each surface element: a cone over the element
    // `N du dv` at `q` (from the apex) has volume `q·N/3`, its centroid at `3/4 q` and its
    // second moment `3/5 q qᵀ` times its volume.
    let mut moments = Moments::default();
    for face in faces(model, shape, options) {
        integrate_face(model, face, |p, n, weight| {
            let q = p - *moments.origin.get_or_insert(p);
            moments.add(p, q.dot(n) * weight / 3.0, 0.75, 0.6);
        });
    }
    moments.properties()
}

/// A face's terms of the volume a closed shell bounds and of its first moment: `∫∫ p·N / 3`
/// about the origin and the same weighted by the centroid of each cone from the origin
/// (`3/4 p`), signed by the face's orientation (a shell's volume and moment are the sums
/// over its faces). Kept in the model for an unplaced face.
pub(crate) fn face_volume_terms(model: &Model, face: Shape) -> (f64, Vector) {
    let (volume, moment, _) = face_terms(model, face);
    (volume, moment)
}

/// A face's volume terms (see [`face_volume_terms`]) and its area's moments (about its first
/// sample: signed by the face's orientation, as is its area), found in one pass. Kept in the
/// model for an unplaced face.
fn face_terms(model: &Model, face: Shape) -> (f64, Vector, Moments) {
    struct FaceTerms(f64, Vector, Moments);
    let compute = |f: Shape| {
        let (mut volume, mut moment, mut area) = (0.0, Vector::new(0.0, 0.0, 0.0), Moments::default());
        integrate_face(model, f, |p, n, weight| {
            let q = p - Point::ORIGIN;
            let m = q.dot(n) * weight / 3.0;
            volume += m;
            moment = moment + q * (0.75 * m);
            area.add(p, n.length() * weight, 1.0, 1.0);
        });
        (volume, moment, area)
    };
    match FaceId::try_from(face.node) {
        Ok(id) if face.placement == crate::math::Transform::IDENTITY => {
            let forward = Shape { orientation: Orientation::Forward, ..face };
            let sign = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
            let kept = model.derived(face.node, model.face_fingerprint(id), || {
                let (v, m, a) = compute(forward);
                FaceTerms(v, m, a)
            });
            (sign * kept.0, kept.1 * sign, kept.2.signed(sign))
        }
        _ => compute(face),
    }
}

/// A face occurrence's area moments (see [`face_terms`]): of a face placed without
/// reflecting or scaling it, its unplaced face's (kept in the model), moved.
fn face_area_moments(model: &Model, face: Shape) -> Moments {
    let t = face.placement;
    if t == crate::math::Transform::IDENTITY || (t.scale_factor() - 1.0).abs() > 1e-12 || t.is_negative() {
        return face_terms(model, face).2;
    }
    let local = face_terms(model, Shape { placement: crate::math::Transform::IDENTITY, ..face }).2;
    // About the moved origin, a point's q turns with it: first R f, second R S Rᵀ.
    let r: [Vector; 3] = std::array::from_fn(|c| Vector::new(f64::from(c == 0), f64::from(c == 1), f64::from(c == 2)).transformed(&t));
    let rot = |row: usize, col: usize| r[col][row];
    let mut second = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            second[i][j] = (0..3).map(|a| (0..3).map(|b| rot(i, a) * local.second[a][b] * rot(j, b)).sum::<f64>()).sum();
        }
    }
    Moments { origin: local.origin.map(|o| o.transformed(&t)), measure: local.measure, first: local.first.transformed(&t), second }
}

/// The volume a closed shell bounds, from its faces' terms (kept in the model).
pub(crate) fn closed_shell_volume(model: &Model, shell: Shape) -> f64 {
    model.explore(shell, NodeKind::Face).map(|f| face_volume_term(model, f)).sum()
}

/// A face's term of the volume a closed shell bounds (see [`face_volume_terms`]).
pub(crate) fn face_volume_term(model: &Model, face: Shape) -> f64 {
    face_volume_terms(model, face).0
}

/// The volume and centre of mass of a shape made of solids (every face of it on one), from
/// its faces' terms, which the model keeps: much quicker than [`volume`] for faces met
/// before. `None` for other shapes.
pub fn solid_volume(model: &Model, shape: Shape) -> Option<(f64, Point)> {
    let faces: Vec<Shape> = model.explore(shape, NodeKind::Face).collect();
    let on_solids: usize = model.explore(shape, NodeKind::Solid).map(|s| model.explore(s, NodeKind::Face).count()).sum();
    if faces.is_empty() || on_solids != faces.len() {
        return None;
    }
    // A solid placed without reflecting or scaling it bounds the volume it bounds unplaced,
    // its centre moved: found from its unplaced faces' kept terms (a lattice's cells, each a
    // placed copy of one, at the cost of one). Others from their placed faces' terms.
    let terms = |faces: &mut dyn Iterator<Item = Shape>| {
        faces.fold((0.0, Vector::new(0.0, 0.0, 0.0)), |(v, m), f| {
            let (fv, fm) = face_volume_terms(model, f);
            (v + fv, m + fm)
        })
    };
    let rigid = |t: &crate::math::Transform| (t.scale_factor() - 1.0).abs() <= 1e-12 && !t.is_negative();
    let (volume, moment) = model.explore(shape, NodeKind::Solid).fold((0.0, Vector::new(0.0, 0.0, 0.0)), |(v, m), s| {
        let (sv, sm) = if s.placement != crate::math::Transform::IDENTITY && rigid(&s.placement) {
            let (sv, sm) = terms(&mut model.explore(Shape { placement: crate::math::Transform::IDENTITY, ..s }, NodeKind::Face));
            (sv, sm.transformed(&s.placement) + s.placement.translation_part() * sv)
        } else {
            terms(&mut model.explore(s, NodeKind::Face))
        };
        (v + sv, m + sm)
    });
    let center = if volume != 0.0 { Point::ORIGIN + moment * (1.0 / volume) } else { Point::ORIGIN };
    Some((volume, center))
}

/// The face occurrences of `shape` (each distinct face once with `skip_shared`).
fn faces(model: &Model, shape: Shape, options: Options) -> impl Iterator<Item = Shape> + '_ {
    let all: Vec<Shape> = match options.skip_shared {
        true => model.unique_subshapes(shape, NodeKind::Face),
        false => model.explore(shape, NodeKind::Face).collect(),
    };
    all.into_iter()
}

/// Calls `sample(point, normal, weight)` for quadrature samples of a face occurrence such
/// that `Σ weight·f(point, normal)` approximates `∫∫ f dudv` over the face's parameter
/// domain, with `normal = Su × Sv` and points and normals in assembly coordinates. The
/// weights are signed by the direction of travel of the boundary (composed with the
/// face's orientation), so a reversed face integrates with the opposite sign.
///
/// The domain integral becomes a boundary integral by Green's theorem,
/// `∫∫ g dudv = ∮ G dv` with `G(u, v) = ∫_{u0}^{u} g(s, v) ds`, both evaluated by
/// Gauss–Legendre quadrature over short pieces of the pcurves and of the inner span.
fn integrate_face(model: &Model, face: Shape, mut sample: impl FnMut(Point, Vector, f64)) {
    let id = FaceId::from_index(face.node.index);
    let record = model.get(id).expect("face in model");
    let surface = &model.get(record.surface).expect("face with a surface").surface;
    // Gauss pieces per unit of u and of v: enough for an angle (a piece spans at most
    // MAX_ANGLE_STEP), for a curve's knot spans, none for a straight direction.
    let angular = 1.0 / MAX_ANGLE_STEP;
    let along = |c: &crate::geometry::curve::Curve| -> f64 {
        use crate::geometry::curve::Curve;
        match c {
            Curve::Line(_) => 0.0,
            Curve::Circle(_) | Curve::Ellipse(_) => angular,
            Curve::Bezier(b) => 4.0 * b.degree() as f64,
            Curve::BSpline(b) => {
                let k = b.knots();
                let shortest = k.windows(2).map(|w| w[1] - w[0]).filter(|&h| h > 0.0).fold(f64::INFINITY, f64::min);
                if shortest.is_finite() { 2.0 * b.degree() as f64 / shortest } else { 0.0 }
            }
            // Unbounded conics and other curves: by parameter steps of the order of one.
            _ => 4.0,
        }
    };
    fn densities(surface: &Surface, angular: f64, along: &dyn Fn(&crate::geometry::curve::Curve) -> f64) -> (f64, f64) {
        match surface {
            Surface::Plane(_) => (0.0, 0.0),
            Surface::Cylinder(_) | Surface::Cone(_) => (angular, 0.0),
            Surface::Sphere(_) | Surface::Torus(_) => (angular, angular),
            // About a piece per knot span (at most 8 across the surface), so that each
            // Gauss rule sees a few polynomial spans, not dozens; more where the surface
            // turns further (a helical sweep winding many times), each piece short in angle.
            Surface::BSpline(b) => {
                use crate::geometry::ParametricSurface;
                let angle = |x: Vector, y: Vector| {
                    let (lx, ly) = (x.length(), y.length());
                    if lx > 0.0 && ly > 0.0 { (x.dot(y) / (lx * ly)).clamp(-1.0, 1.0).acos() } else { 0.0 }
                };
                // (Along the middle iso-lines, over chords long enough that a fitted surface's
                // ripples between its knots do not add up.)
                let (u0, u1, v0, v1) = b.bounds();
                const CHORDS: usize = 128;
                let turned = |at: &dyn Fn(f64) -> Point| {
                    let points: Vec<Point> = (0..=CHORDS).map(|i| at(i as f64 / CHORDS as f64)).collect();
                    points.windows(3).map(|w| angle(w[1] - w[0], w[2] - w[1])).sum::<f64>()
                };
                let u_turn = turned(&|s| b.point(u0 + (u1 - u0) * s, 0.5 * (v0 + v1)));
                let v_turn = turned(&|s| b.point(0.5 * (u0 + u1), v0 + (v1 - v0) * s));
                let density = |knots: &[f64], turn: f64| {
                    let spans = knots.len().saturating_sub(1);
                    let range = knots[knots.len() - 1] - knots[0];
                    let pieces = (spans.min(8) as f64).max((turn / MAX_ANGLE_STEP).ceil().min(4.0 * spans as f64));
                    if spans > 1 && range > 0.0 { pieces / range } else { 0.0 }
                };
                (density(b.u_knots(), u_turn), density(b.v_knots(), v_turn))
            }
            Surface::Bezier(_) => (0.0, 0.0),
            Surface::Revolution(r) => (angular, along(r.basis())),
            Surface::Extrusion(e) => (along(e.basis()), 0.0),
            // Trimmed and offset surfaces are parametrized as their bases.
            Surface::Trimmed(t) => densities(t.basis(), angular, along),
            Surface::Offset(o) => densities(o.basis(), angular, along),
        }
    }
    let (u_density, v_density) = densities(surface, angular, &along);
    // (Along a pcurve, v counts too where the surface is straight across it (a ruled face's
    // single span), a few pieces across its range: the integrand along u changes with v all
    // the same, and a pcurve turning sharply across it (where a face narrows to nothing) is
    // integrated piece by piece.)
    let v_weight = {
        let (_, _, v0, v1) = surface.bounds();
        if matches!(surface, Surface::Plane(_)) || !(v1 - v0).is_finite() || v1 - v0 >= 1e100 || !(v1 > v0) { v_density } else { v_density.max(4.0 / (v1 - v0)) }
    };
    // An extrusion of a curve without a density of its own (a parabola, a hyperbola, ...):
    // pieces along u where its tangent turns, graded (by the parameter, a parabola's turn
    // is all near its vertex, however long it is).
    fn turning_basis(surface: &Surface) -> Option<&crate::geometry::curve::Curve> {
        use crate::geometry::curve::Curve;
        match surface {
            Surface::Extrusion(e) => (!matches!(e.basis(), Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) | Curve::Bezier(_) | Curve::BSpline(_))).then(|| e.basis()),
            Surface::Trimmed(t) => turning_basis(t.basis()),
            Surface::Offset(o) => turning_basis(o.basis()),
            _ => None,
        }
    }
    let u_turning = turning_basis(surface);
    // The angle travelled along u from `a` to `b`.
    let mut scratch = Vec::new();
    let mut u_span = |a: f64, b: f64| match u_turning {
        Some(c) => {
            scratch.clear();
            turning_breaks(c, a, b, &mut scratch)
        }
        None => (b - a).abs() * u_density * MAX_ANGLE_STEP,
    };
    // Knots of a B-spline surface where it is at most C1 (as at the joins of a piecewise
    // Bézier surface): no Gauss rule may span one.
    fn kinks(surface: &Surface) -> (Vec<f64>, Vec<f64>) {
        match surface {
            Surface::BSpline(b) => {
                let low = |knots: &[f64], mults: &[usize], degree: usize| -> Vec<f64> {
                    let n = knots.len();
                    (1..n.saturating_sub(1)).filter(|&i| mults[i] + 1 >= degree).map(|i| knots[i]).collect()
                };
                (low(b.u_knots(), b.u_multiplicities(), b.u_degree()), low(b.v_knots(), b.v_multiplicities(), b.v_degree()))
            }
            Surface::Trimmed(t) => kinks(t.basis()),
            Surface::Offset(o) => kinks(o.basis()),
            _ => (Vec::new(), Vec::new()),
        }
    }
    let (u_kinks, v_kinks) = kinks(surface);
    // The parameters in (a, b) where `f` crosses one of `values`, by bisection.
    let crossings = |f: &dyn Fn(f64) -> f64, a: f64, b: f64, values: &[f64]| -> Vec<f64> {
        let mut out = Vec::new();
        if values.is_empty() {
            return out;
        }
        let n = 16;
        for i in 0..n {
            let (s0, s1) = (a + (b - a) * i as f64 / n as f64, a + (b - a) * (i + 1) as f64 / n as f64);
            let (f0, f1) = (f(s0), f(s1));
            for &k in values {
                if (f0 - k) * (f1 - k) < 0.0 {
                    let (mut lo, mut hi) = (s0, s1);
                    // (f(lo) kept from where lo was found: the same value.)
                    let mut at_lo = f(lo);
                    for _ in 0..60 {
                        let m = 0.5 * (lo + hi);
                        let at_m = f(m);
                        if (at_lo - k) * (at_m - k) <= 0.0 {
                            hi = m
                        } else {
                            (lo, at_lo) = (m, at_m)
                        }
                    }
                    out.push(0.5 * (lo + hi));
                }
            }
        }
        out
    };
    // The length of a parameter step that keeps each Gauss piece short in angle.
    let pieces = |angular_span: f64| ((angular_span / MAX_ANGLE_STEP).ceil() as usize).clamp(1, 256);
    // The normal from the placed derivatives: a scaling stretches the area element by its
    // square.
    let place = |p: Point, su: Vector, sv: Vector| {
        (p.transformed(&face.placement), su.transformed(&face.placement).cross(sv.transformed(&face.placement)))
    };

    let planar = matches!(surface, Surface::Plane(_));
    // (On a B-spline face straight along u, its normal keeping its direction along each line
    // of constant v, as files write cylinders and cones, the integrands along u have degree
    // 3 at most, as on a plane: there too the 2-point rule is exact, a piece between kinks.)
    let straight_u = match surface {
        Surface::BSpline(b) => b.straight_along_u(),
        Surface::Trimmed(t) => matches!(t.basis(), Surface::BSpline(b) if b.straight_along_u()),
        _ => false,
    };
    let inner_rule: &[(f64, f64)] = if planar || straight_u { &GAUSS_LEGENDRE_2 } else { &GAUSS_LEGENDRE };
    // On a surface of revolution (u its angle about the axis) a point, its derivatives and
    // so the integrands of the moments (of degree 3 at most in them) are trigonometric
    // polynomials of degree 3 at most in u: their integral along u is exact from 7 samples
    // a seventh of a turn apart, each the point at the pcurve turned (one evaluation).
    fn revolution_axis(surface: &Surface) -> Option<crate::math::Axis> {
        let about = |f: &crate::math::Frame| {
            let z = Vector::from(f.x()).cross(Vector::from(f.y()));
            Some(crate::math::Axis::new(f.origin(), Direction::new(z.x, z.y, z.z)))
        };
        match surface {
            Surface::Cylinder(c) => about(c.frame()),
            Surface::Cone(c) => about(c.frame()),
            Surface::Sphere(s) => about(s.frame()),
            Surface::Torus(t) => about(t.frame()),
            Surface::Revolution(r) => Some(r.axis()),
            Surface::Trimmed(t) => revolution_axis(t.basis()),
            Surface::Offset(o) => revolution_axis(o.basis()),
            _ => None,
        }
    }
    let turns: Option<[crate::math::Transform; 7]> =
        revolution_axis(surface).map(|axis| std::array::from_fn(|j| crate::math::Transform::rotation(axis, std::f64::consts::TAU * j as f64 / 7.0)));
    let mut u0 = None;
    // (Buffers reused from node to node.)
    let (mut ends, mut nodes, mut us, mut graded): (Vec<f64>, Vec<(f64, f64, f64)>, Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for &wire_use_id in &record.wires {
        let wire_use = model.get(wire_use_id).expect("wire use");
        let wire = model.get(wire_use.wire).expect("wire in model");
        let loop_orientation = wire_use.orientation.compose(face.orientation);
        for &coedge_id in &wire.coedges {
            let coedge = model.get(coedge_id).expect("coedge in model");
            let (pcurve, first, last) = model
                .effective_pcurve(coedge_id)
                .unwrap_or_else(|| todo!("properties of a face without pcurves ({coedge_id:?})"));
            let sign = match coedge.orientation.compose(loop_orientation) {
                Orientation::Reversed => -1.0,
                _ => 1.0,
            };
            // Stretches that each Gauss rule integrates: between knots of a B-spline (where
            // it is only finitely smooth; with many knots, groups of spans, the joins being
            // smooth enough), short in angle.
            let mut breaks = vec![first];
            if let Curve2::BSpline(b) = &pcurve {
                let (lo, hi) = (first.min(last), first.max(last));
                let mut inner: Vec<f64> = b.knots().iter().copied().filter(|&k| k > lo && k < hi).collect();
                if last < first {
                    inner.reverse();
                }
                // Group spans only where they join smoothly (not at a polyline's kinks).
                let group = if b.degree() >= 2 { inner.len().div_ceil(32).max(1) } else { 1 };
                breaks.extend(inner.into_iter().skip(group - 1).step_by(group));
            }
            breaks.push(last);
            if !u_kinks.is_empty() || !v_kinks.is_empty() {
                let mut more = Vec::new();
                for w in breaks.windows(2) {
                    more.extend(crossings(&|t| pcurve.point(t).x, w[0], w[1], &u_kinks));
                    more.extend(crossings(&|t| pcurve.point(t).y, w[0], w[1], &v_kinks));
                }
                breaks.extend(more);
                if last < first {
                    breaks.sort_by(|a, b| b.total_cmp(a));
                } else {
                    breaks.sort_by(f64::total_cmp);
                }
                breaks.dedup_by(|a, b| (*a - *b).abs() <= 1e-15 * (1.0 + a.abs()));
            }
            let mut stretches: Vec<(f64, f64)> = Vec::new();
            for w in breaks.windows(2) {
                let (s0, s1) = (w[0], w[1]);
                let (a, b) = (pcurve.point(s0), pcurve.point(s1));
                let span = match pcurve {
                    Curve2::Circle(_) => (s1 - s0).abs(),
                    Curve2::Line(_) => u_span(a.x, b.x) + (b.y - a.y).abs() * v_density * MAX_ANGLE_STEP,
                    _ => {
                        // The angle travelled, from a few samples; where it is travelled
                        // unevenly (most of it in one part, as round a sharp turn), halved
                        // until it is even or short, so that the pieces gather there.
                        let mut todo = vec![(s0, s1, 0u32)];
                        while let Some((a, b, depth)) = todo.pop() {
                            let parts: Vec<f64> = (0..8)
                                .map(|i| {
                                    let (p, q) = (pcurve.point(a + (b - a) * i as f64 / 8.0), pcurve.point(a + (b - a) * (i + 1) as f64 / 8.0));
                                    u_span(p.x, q.x) + (q.y - p.y).abs() * v_weight * MAX_ANGLE_STEP
                                })
                                .collect();
                            let span: f64 = parts.iter().sum();
                            let most = parts.iter().copied().fold(0.0, f64::max);
                            if span > MAX_ANGLE_STEP && most > 0.5 * span && depth < 16 {
                                let m = 0.5 * (a + b);
                                todo.push((m, b, depth + 1));
                                todo.push((a, m, depth + 1));
                                continue;
                            }
                            let n = pieces(span);
                            let step = (b - a) / n as f64;
                            stretches.extend((0..n).map(|k| (a + step * k as f64, step)));
                        }
                        continue;
                    }
                };
                let n = pieces(span);
                let step = (s1 - s0) / n as f64;
                stretches.extend((0..n).map(|k| (s0 + step * k as f64, step)));
            }
            // On a plane the integrands are polynomials: of degree 3 at most in u (a cone's
            // second moments), and with a straight pcurve of degree 4 at most along it, which
            // Gauss rules of 2 and 3 points integrate exactly.
            let outer_rule: &[(f64, f64)] = if planar && matches!(pcurve, Curve2::Line(_)) { &GAUSS_LEGENDRE_3 } else { &GAUSS_LEGENDRE };
            for (start, step) in stretches {
                for &(node, w) in outer_rule {
                    let t = start + step * (node + 1.0) / 2.0;
                    let (uv, duv) = pcurve.d1(t);
                    let outer_weight = sign * w * step / 2.0 * duv.y;
                    if outer_weight == 0.0 {
                        continue;
                    }
                    let u_start = *u0.get_or_insert(uv.x);
                    if let Some(turns) = &turns {
                        // ∫ from u_start to u of the samples' interpolant (the Dirichlet
                        // kernel D(s) = 1 + 2 Σ cos ks about each), in a form exact for
                        // short spans.
                        let span = uv.x - u_start;
                        let half: [(f64, f64); 3] = std::array::from_fn(|k| ((k + 1) as f64 * span / 2.0).sin_cos());
                        let (p, su, sv) = surface.d1(uv.x, uv.y);
                        for (j, turn) in turns.iter().enumerate() {
                            // (cos k(s/2 + θ) by the sum rule, k θ's sine and cosine tabled.)
                            let kernel: f64 = (0..3).map(|k| {
                                let (sin_t, cos_t) = SEVENTHS[(k + 1) * j % 7];
                                half[k].0 * (half[k].1 * cos_t - half[k].0 * sin_t) / (k + 1) as f64
                            }).sum();
                            let (p, n) = place(p.transformed(turn), su.transformed(turn), sv.transformed(turn));
                            sample(p, n, outer_weight * (span + 4.0 * kernel) / 7.0);
                        }
                        continue;
                    }
                    // From u_start to the point, in stretches between the surface's kinks.
                    ends.clear();
                    ends.push(u_start);
                    let (lo, hi) = (u_start.min(uv.x), u_start.max(uv.x));
                    let inside = u_kinks.iter().copied().filter(|&k| k > lo && k < hi);
                    if uv.x >= u_start { ends.extend(inside) } else { ends.extend(inside.rev()) }
                    ends.push(uv.x);
                    // The nodes along u at this v (evaluated together on a B-spline).
                    nodes.clear();
                    for e in ends.windows(2) {
                        if let Some(c) = u_turning {
                            graded.clear();
                            turning_breaks(c, e[0], e[1], &mut graded);
                            let mut s0 = e[0];
                            for &s1 in &graded {
                                let inner_step = s1 - s0;
                                for &(inner_node, inner_w) in inner_rule {
                                    nodes.push((s0 + inner_step * (inner_node + 1.0) / 2.0, inner_w, inner_step));
                                }
                                s0 = s1;
                            }
                            continue;
                        }
                        let inner = if straight_u { 1 } else { pieces((e[1] - e[0]).abs() * u_density * MAX_ANGLE_STEP) };
                        let inner_step = (e[1] - e[0]) / inner as f64;
                        for inner_piece in 0..inner {
                            let s0 = e[0] + inner_step * inner_piece as f64;
                            for &(inner_node, inner_w) in inner_rule {
                                nodes.push((s0 + inner_step * (inner_node + 1.0) / 2.0, inner_w, inner_step));
                            }
                        }
                    }
                    if let Surface::BSpline(b) = surface {
                        us.clear();
                        us.extend(nodes.iter().map(|n| n.0));
                        for ((p, su, sv), &(_, inner_w, inner_step)) in b.d1_along_u(&us, uv.y).into_iter().zip(&nodes) {
                            let (p, n) = place(p, su, sv);
                            sample(p, n, outer_weight * inner_w * inner_step / 2.0);
                        }
                    } else {
                        for &(u, inner_w, inner_step) in &nodes {
                            let (p, su, sv) = surface.d1(u, uv.y);
                            let (p, n) = place(p, su, sv);
                            sample(p, n, outer_weight * inner_w * inner_step / 2.0);
                        }
                    }
                }
            }
        }
    }
}

/// The longest angular parameter span integrated by one Gauss–Legendre rule.
/// Splits `a`..`b` (either way) where `curve`'s tangent turns by at most [`MAX_ANGLE_STEP`]
/// on each piece, the pieces' ends (after `a`, `b` last) pushed to `out`; the angle turned.
fn turning_breaks(curve: &crate::geometry::curve::Curve, a: f64, b: f64, out: &mut Vec<f64>) -> f64 {
    let direction = |t: f64| curve.d1(t).1;
    let angle = |x: Vector, y: Vector| {
        let (lx, ly) = (x.length(), y.length());
        if lx > 0.0 && ly > 0.0 { (x.dot(y) / (lx * ly)).clamp(-1.0, 1.0).acos() } else { 0.0 }
    };
    // (Halved until the tangents at the ends and the middle turn by at most the step, up to
    // 4096 pieces.)
    let mut stack = vec![(b, direction(b), 0u32)];
    let (mut at, mut d_at) = (a, direction(a));
    let mut turned = 0.0;
    while let Some(&(end, d_end, depth)) = stack.last() {
        let middle = 0.5 * (at + end);
        let d_middle = direction(middle);
        let turn = angle(d_at, d_middle) + angle(d_middle, d_end);
        if turn <= MAX_ANGLE_STEP || depth >= 12 || middle == at || middle == end {
            stack.pop();
            out.push(end);
            turned += turn;
            (at, d_at) = (end, d_end);
        } else {
            stack.last_mut().expect("the piece's end").2 = depth + 1;
            stack.push((middle, d_middle, depth + 1));
        }
    }
    turned
}

const MAX_ANGLE_STEP: f64 = std::f64::consts::PI / 4.0;

/// The sine and cosine of `2π j / 7` for `j` in `0..7`.
const SEVENTHS: [(f64, f64); 7] = [
    (0.0, 1.0),
    (0.781_831_482_468_029_8, 0.623_489_801_858_733_5),
    (0.974_927_912_181_823_6, -0.222_520_933_956_314_4),
    (0.433_883_739_117_558_1, -0.900_968_867_902_419_1),
    (-0.433_883_739_117_558_1, -0.900_968_867_902_419_1),
    (-0.974_927_912_181_823_6, -0.222_520_933_956_314_4),
    (-0.781_831_482_468_029_8, 0.623_489_801_858_733_5),
];

/// The 2-point Gauss–Legendre rule on [-1, 1] (exact for cubics).
const GAUSS_LEGENDRE_2: [(f64, f64); 2] = [(-0.577_350_269_189_625_8, 1.0), (0.577_350_269_189_625_8, 1.0)];

/// The 3-point Gauss–Legendre rule on [-1, 1] (exact for polynomials of degree 5).
const GAUSS_LEGENDRE_3: [(f64, f64); 3] = [(-0.774_596_669_241_483_4, 0.555_555_555_555_555_6), (0.0, 0.888_888_888_888_888_9), (0.774_596_669_241_483_4, 0.555_555_555_555_555_6)];

/// The 10-point Gauss–Legendre rule on [-1, 1]: (node, weight).
const GAUSS_LEGENDRE: [(f64, f64); 10] = [
    (-0.973_906_528_517_171_7, 0.066_671_344_308_688_14),
    (-0.865_063_366_688_984_5, 0.149_451_349_150_580_6),
    (-0.679_409_568_299_024_4, 0.219_086_362_515_982_04),
    (-0.433_395_394_129_247_2, 0.269_266_719_309_996_35),
    (-0.148_874_338_981_631_2, 0.295_524_224_714_752_87),
    (0.148_874_338_981_631_2, 0.295_524_224_714_752_87),
    (0.433_395_394_129_247_2, 0.269_266_719_309_996_35),
    (0.679_409_568_299_024_4, 0.219_086_362_515_982_04),
    (0.865_063_366_688_984_5, 0.149_451_349_150_580_6),
    (0.973_906_528_517_171_7, 0.066_671_344_308_688_14),
];

/// Accumulates a measure and its first and second moments from weighted samples,
/// relative to the first sample's point.
#[derive(Default, Clone)]
struct Moments {
    origin: Option<Point>,
    measure: f64,
    first: Vector,
    second: [[f64; 3]; 3],
}

impl Moments {
    /// Adds the measure `dm` at `p`, its centroid at `centroid · q` and its second moment
    /// `second · q qᵀ · dm`, where `q` is `p` relative to the origin.
    fn add(&mut self, p: Point, dm: f64, centroid: f64, second: f64) {
        let q = p - *self.origin.get_or_insert(p);
        self.measure += dm;
        self.first = self.first + q * (centroid * dm);
        for r in 0..3 {
            for c in 0..3 {
                self.second[r][c] += second * dm * q[r] * q[c];
            }
        }
    }

    /// The accumulation times `sign`.
    fn signed(&self, sign: f64) -> Moments {
        let mut out = Moments { origin: self.origin, ..Moments::default() };
        out.merge(self, sign);
        out
    }

    /// Adds `sign` times another accumulation of an area (its centroids and second moments
    /// those of its points, as [`Self::add`] with factors 1), moved to this one's origin.
    fn merge_area(&mut self, other: &Moments, sign: f64) {
        let (Some(origin), Some(theirs)) = (self.origin, other.origin) else {
            return self.merge(other, sign);
        };
        // About this origin, a point's q is its own q + d.
        let d = theirs - origin;
        let (m, first) = (sign * other.measure, other.first * sign);
        self.measure += m;
        self.first = self.first + first + d * m;
        for r in 0..3 {
            for c in 0..3 {
                self.second[r][c] += sign * other.second[r][c] + first[r] * d[c] + d[r] * first[c] + m * d[r] * d[c];
            }
        }
    }

    /// Adds `sign` times another accumulation (which must share the origin).
    fn merge(&mut self, other: &Moments, sign: f64) {
        if self.origin.is_none() {
            self.origin = other.origin;
        }
        assert!(other.measure == 0.0 || other.origin == self.origin, "moments with different origins");
        self.measure += sign * other.measure;
        self.first = self.first + other.first * sign;
        for r in 0..3 {
            for c in 0..3 {
                self.second[r][c] += sign * other.second[r][c];
            }
        }
    }

    /// The accumulated measure, centroid and inertia about the centroid.
    fn properties(&self) -> MassProperties {
        let origin = self.origin.unwrap_or(Point::ORIGIN);
        if self.measure == 0.0 {
            return MassProperties { mass: 0.0, center_of_mass: origin, inertia: Matrix3::ZERO };
        }
        let g = self.first / self.measure;
        let mut central = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                central[r][c] = self.second[r][c] - self.measure * g[r] * g[c];
            }
        }
        let trace = central[0][0] + central[1][1] + central[2][2];
        let inertia = std::array::from_fn(|r| {
            std::array::from_fn(|c| if r == c { trace } else { 0.0 } - central[r][c])
        });
        MassProperties { mass: self.measure, center_of_mass: origin + g, inertia: Matrix3::from_rows(inertia) }
    }
}

/// A box enclosing `shape` (not necessarily the tightest one).
#[doc(alias = "BRepBndLib::Add")]
pub fn bounding_box(model: &Model, shape: Shape) -> BoundingBox {
    let mut b = BoundingBox::EMPTY;
    // A face moved without turning: its whole box (surface, edges and vertices), found once
    // unplaced and kept with it, then moved. Other faces, edges and vertices each.
    let mut faces: Vec<Shape> = Vec::new();
    for s in std::iter::once(shape).chain(model.explore_all(shape).stop_at(NodeKind::Face)) {
        match s.node.kind {
            NodeKind::Face => faces.push(s),
            NodeKind::Vertex => b.add_box(&vertex_box(model, s)),
            NodeKind::Edge => b.add_box(&edge_box(model, s)),
            _ => {}
        }
    }
    let boxes = crate::parallel::map(&faces, |&f| {
        let moved_only = crate::math::Transform::translation(f.placement.translation_part()) == f.placement;
        if moved_only {
            return whole_face_box(model, f);
        }
        let mut b = BoundingBox::EMPTY;
        for v in model.explore(f, NodeKind::Vertex) {
            b.add_box(&vertex_box(model, v));
        }
        for e in model.explore(f, NodeKind::Edge) {
            b.add_box(&edge_box(model, e));
        }
        if let Some(surface) = face_box(model, f) {
            b.add_box(&surface);
        }
        b
    });
    for face_box in boxes {
        b.add_box(&face_box);
    }
    b
}

/// A vertex's part of [`bounding_box`]: its point, enlarged by its tolerance.
fn vertex_box(model: &Model, v: Shape) -> BoundingBox {
    let p = model.point(v);
    BoundingBox::new(p, p).enlarged(model.tolerance(v))
}

/// An edge's part of [`bounding_box`]: points along its curve, enlarged to cover the curve
/// between them (found once unplaced, kept with the edge, then placed).
fn edge_box(model: &Model, edge: Shape) -> BoundingBox {
    struct EdgeBox(Option<BoundingBox>);
    let Ok(id) = crate::topology::EdgeId::try_from(edge.node) else { return BoundingBox::EMPTY };
    let local = Shape { orientation: Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..edge };
    let kept = model.derived(edge.node, model.edge_fingerprint(id), || {
        EdgeBox(model.edge_curve(local).map(|(curve, first, last)| {
            let n = if matches!(curve, crate::geometry::curve::Curve::Line(_)) { 1 } else { 64 };
            let points: Vec<Point> = (0..=n).map(|i| curve.point(first + (last - first) * i as f64 / n as f64)).collect();
            // A curve bulges from a chord by at most a quarter of it (arcs up to half a turn).
            let chord = points.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max);
            let margin = if n == 1 { 0.0 } else { 0.25 * chord };
            BoundingBox::from_points(points).enlarged(margin + model.tolerance(local))
        }))
    });
    kept.0.map_or(BoundingBox::EMPTY, |e| placed_box(&e, &edge.placement))
}

/// A face's whole part of [`bounding_box`] (its vertices, edges and, for curved ones, the
/// surface within them), found once unplaced and kept with the face, then placed.
fn whole_face_box(model: &Model, f: Shape) -> BoundingBox {
    struct WholeFaceBox(BoundingBox);
    let Ok(id) = FaceId::try_from(f.node) else { return BoundingBox::EMPTY };
    let local = Shape { orientation: Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..f };
    let kept = model.derived(f.node, model.face_fingerprint(id), || {
        let mut b = BoundingBox::EMPTY;
        for v in model.explore(local, NodeKind::Vertex) {
            b.add_box(&vertex_box(model, v));
        }
        for e in model.explore(local, NodeKind::Edge) {
            b.add_box(&edge_box(model, e));
        }
        if let Some(surface) = face_box(model, local) {
            b.add_box(&surface);
        }
        WholeFaceBox(b)
    });
    if kept.0.is_empty() { kept.0 } else { placed_box(&kept.0, &f.placement) }
}

/// A box along axes of its own (an oriented bounding box): `center ± half_sizes[i] *
/// axes[i]`.
#[doc(alias = "Bnd_OBB")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedBox {
    pub center: Point,
    pub axes: [Vector; 3],
    pub half_sizes: [f64; 3],
}

impl OrientedBox {
    /// Whether the boxes overlap (no plane separates them: none across their 15
    /// candidate axes, the boxes' own and the cross products of a pair).
    #[doc(alias = "Bnd_OBB::IsOut")]
    pub fn intersects(&self, other: &OrientedBox) -> bool {
        let gap = other.center - self.center;
        let reach = |b: &OrientedBox, axis: Vector| (0..3).map(|i| b.half_sizes[i] * b.axes[i].dot(axis).abs()).sum::<f64>();
        let mut axes: Vec<Vector> = self.axes.iter().chain(other.axes.iter()).copied().collect();
        for a in self.axes {
            for b in other.axes {
                let c = a.cross(b);
                if c.length() > 1e-9 {
                    axes.push(c.normalized());
                }
            }
        }
        axes.iter().all(|&axis| gap.dot(axis).abs() <= reach(self, axis) + reach(other, axis))
    }
}

/// A box along axes fitting `shape` (not necessarily the smallest one): the smallest of
/// boxes along the frames of its planar faces (their normal and an edge's direction) and
/// along the principal axes of points spread over it, around those points enlarged by the
/// tolerances.
#[doc(alias = "BRepBndLib::AddOBB")]
pub fn oriented_bounding_box(model: &Model, shape: Shape) -> OrientedBox {
    // Points: the vertices, along the edges and over the curved faces.
    let mut points: Vec<Point> = model.explore(shape, NodeKind::Vertex).map(|v| model.point(v)).collect();
    let mut margin = model.explore(shape, NodeKind::Vertex).map(|v| model.tolerance(v)).fold(0.0, f64::max);
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        margin = margin.max(model.tolerance(edge));
        if let Some((curve, first, last)) = model.edge_curve(edge) {
            let n = if matches!(curve, crate::geometry::curve::Curve::Line(_)) { 1 } else { 64 };
            points.extend((0..=n).map(|i| curve.point(first + (last - first) * i as f64 / n as f64)));
        }
    }
    let faces = model.unique_subshapes(shape, NodeKind::Face);
    let mut frames: Vec<[Vector; 3]> = Vec::new();
    for &f in &faces {
        if let Some(b) = face_box(model, f) {
            // (A curved face: its box's corners stand for it.)
            for i in 0..8 {
                points.push(Point::new(
                    if i & 1 == 0 { b.min.x } else { b.max.x },
                    if i & 2 == 0 { b.min.y } else { b.max.y },
                    if i & 4 == 0 { b.min.z } else { b.max.z },
                ));
            }
        }
        if let Some(Surface::Plane(plane)) = model.surface(f).map(|s| s.transformed(&f.placement)) {
            let n = Vector::from(plane.normal());
            for e in model.explore(f, NodeKind::Edge).take(4) {
                let Some((Some(curve), first, last)) = model.edge_curve(e).map(|(c, a, b)| (Some(c), a, b)) else { continue };
                let d = curve.point(last) - curve.point(first);
                let x = d - n * d.dot(n);
                if x.length() > 1e-9 {
                    let x = x.normalized();
                    frames.push([n, x, n.cross(x)]);
                }
            }
        }
    }
    if points.is_empty() {
        return OrientedBox { center: Point::ORIGIN, axes: [Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, 1.0)], half_sizes: [0.0; 3] };
    }
    // The principal axes of the points, and the world's.
    let mean = points.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + (*p - Point::ORIGIN)) * (1.0 / points.len() as f64);
    let mut covariance = [[0.0; 3]; 3];
    for p in &points {
        let d = (*p - Point::ORIGIN) - mean;
        let d = [d.x, d.y, d.z];
        for r in 0..3 {
            for c in 0..3 {
                covariance[r][c] += d[r] * d[c];
            }
        }
    }
    let (_, principal) = symmetric_eigen(&Matrix3::from_rows(covariance));
    frames.push(principal.map(|v| v.normalized()));
    frames.push([Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0), Vector::new(0.0, 0.0, 1.0)]);
    let along = |axes: &[Vector; 3]| {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &points {
            let d = *p - Point::ORIGIN;
            for i in 0..3 {
                let x = d.dot(axes[i]);
                (lo[i], hi[i]) = (lo[i].min(x), hi[i].max(x));
            }
        }
        (lo, hi)
    };
    let mut best: Option<(f64, OrientedBox)> = None;
    for axes in frames {
        let (lo, hi) = along(&axes);
        let half_sizes: [f64; 3] = std::array::from_fn(|i| 0.5 * (hi[i] - lo[i]) + margin);
        let center = Point::ORIGIN + (0..3).fold(Vector::new(0.0, 0.0, 0.0), |s, i| s + axes[i] * (0.5 * (lo[i] + hi[i])));
        let volume = half_sizes.iter().product::<f64>();
        if best.as_ref().is_none_or(|(v, _)| volume < *v) {
            best = Some((volume, OrientedBox { center, axes, half_sizes }));
        }
    }
    best.expect("a frame").1
}

/// A curved face's part of [`bounding_box`]: points of its surface over the part of the
/// parameter space its boundary spans, enlarged to cover the surface between them.
fn face_box(model: &Model, f: Shape) -> Option<BoundingBox> {
    // (Found once in the face's own coordinates, kept with it, then placed: a shape moved
    // about keeps its faces.)
    struct FaceBox(Option<BoundingBox>);
    let id = FaceId::try_from(f.node).ok()?;
    let local = Shape { orientation: Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..f };
    let kept = model.derived(f.node, model.face_fingerprint(id), || FaceBox(local_face_box(model, local)));
    kept.0.map(|b| placed_box(&b, &f.placement))
}

/// The box of `b`'s corners placed by `t` (exact under a translation, around it otherwise).
fn placed_box(b: &BoundingBox, t: &crate::math::Transform) -> BoundingBox {
    if t.is_identity() {
        return *b;
    }
    BoundingBox::from_points((0..8).map(|i| {
        Point::new(
            if i & 1 == 0 { b.min.x } else { b.max.x },
            if i & 2 == 0 { b.min.y } else { b.max.y },
            if i & 4 == 0 { b.min.z } else { b.max.z },
        )
        .transformed(t)
    }))
}

/// [`face_box`] of an unplaced face.
fn local_face_box(model: &Model, f: Shape) -> Option<BoundingBox> {
    let surface = model.surface(f)?;
    let surface = surface.transformed(&f.placement);
    if matches!(surface, Surface::Plane(_)) {
        return None;
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for wire in model.children(Shape { orientation: Orientation::Forward, ..f }) {
        for coedge in model.children(wire) {
            let Ok(id) = crate::topology::CoEdgeId::try_from(coedge.node) else { continue };
            let Some((pcurve, a, c)) = model.effective_pcurve(id) else { continue };
            for i in 0..=32 {
                let uv = pcurve.point(a + (c - a) * i as f64 / 32.0);
                lo = [lo[0].min(uv.x), lo[1].min(uv.y)];
                hi = [hi[0].max(uv.x), hi[1].max(uv.y)];
            }
        }
    }
    if lo[0] > hi[0] {
        return None;
    }
    let n = 32;
    let grid: Vec<Point> = (0..=n)
        .flat_map(|i| (0..=n).map(move |j| (i, j)))
        .map(|(i, j)| surface.point(lo[0] + (hi[0] - lo[0]) * i as f64 / n as f64, lo[1] + (hi[1] - lo[1]) * j as f64 / n as f64))
        .collect();
    // The largest step between neighbours of the grid (along u or v, not from the end of a
    // row to the start of the next).
    let at = |i: usize, j: usize| grid[i * (n + 1) + j];
    let step = (0..=n).flat_map(|i| (0..=n).map(move |j| (i, j))).fold(0.0f64, |m, (i, j)| {
        let along = if j < n { at(i, j).distance(at(i, j + 1)) } else { 0.0 };
        let across = if i < n { at(i, j).distance(at(i + 1, j)) } else { 0.0 };
        m.max(along).max(across)
    });
    Some(BoundingBox::from_points(grid).enlarged(0.25 * step + model.tolerance(f)))
}

/// A box around `shape` like [`tight_bounding_box`], except that a face carrying a mesh
/// counts by its mesh's nodes, widened by the mesh's deviation from the surface and the
/// face's tolerance (as OCCT bounds triangulated faces).
#[doc(alias = "BRepBndLib::AddOptimal")]
pub fn bounding_box_using_meshes(model: &Model, shape: Shape) -> BoundingBox {
    let faces: Vec<Shape> = model.explore(shape, NodeKind::Face).collect();
    if faces.is_empty() {
        return tight_bounding_box(model, shape);
    }
    let mut b = BoundingBox::EMPTY;
    for face in faces {
        let local = Shape { orientation: Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..face };
        match model.face_mesh(local) {
            Some(mesh) if !mesh.nodes().is_empty() => {
                let nodes = BoundingBox::from_points(mesh.nodes().iter().map(|p| p.transformed(&face.placement)));
                b.add_box(&nodes.enlarged(mesh.deflection() + model.tolerance(face)));
            }
            _ => b.add_box(&tight_bounding_box(model, face)),
        }
    }
    b
}

/// The table of lengths along `edge`'s curve over its range (kept with the edge while it
/// is unchanged); `None` for an edge without a curve.
pub fn edge_length_table(model: &Model, edge: Shape) -> Option<std::sync::Arc<crate::geometry::measure::LengthTable>> {
    let id = crate::topology::EdgeId::try_from(edge.node).ok()?;
    let (curve, first, last) = model.edge_curve(Shape { orientation: crate::topology::Orientation::Forward, placement: crate::math::Transform::IDENTITY, ..edge })?;
    Some(model.derived(edge.node, model.edge_fingerprint(id), || crate::geometry::measure::LengthTable::new(&curve, first, last)))
}

/// The tightest box around `shape`'s geometry (without tolerances): its vertices, the extremes
/// of its edges' curves and of its curved faces' surfaces within their boundaries.
#[doc(alias = "BRepBndLib::AddOptimal")]
pub fn tight_bounding_box(model: &Model, shape: Shape) -> BoundingBox {
    let mut b = BoundingBox::EMPTY;
    // A compound's children bounded each alone: copies of one shape placed apart without
    // turning (a lattice's cells) bounded once, the box moved to each.
    if shape.node.kind == NodeKind::Compound {
        let mut bounded: std::collections::HashMap<crate::topology::NodeId, Vec<(crate::math::Transform, BoundingBox)>> = Default::default();
        for child in model.children(shape) {
            let by = child.placement.translation_part();
            let turned = crate::math::Transform::translation(-by) * child.placement;
            let kept = bounded.entry(child.node).or_default();
            let local = match kept.iter().find(|(t, _)| *t == turned) {
                Some(&(_, local)) => local,
                None => {
                    let local = tight_bounding_box(model, Shape { placement: turned, ..child });
                    kept.push((turned, local));
                    local
                }
            };
            if !local.is_empty() {
                b.add_box(&BoundingBox::new(local.min + by, local.max + by));
            }
        }
        return b;
    }
    for v in model.explore(shape, NodeKind::Vertex) {
        b.add_point(model.point(v));
    }
    for edge in model.unique_subshapes(shape, NodeKind::Edge) {
        let (points, by) = edge_extreme_points(model, edge);
        for &p in points.iter() {
            b.add_point(p + by);
        }
    }
    for f in model.unique_subshapes(shape, NodeKind::Face) {
        let Some(surface) = model.surface(f) else { continue };
        if matches!(surface, Surface::Plane(_)) {
            continue; // bounded by its edges
        }
        let surface = surface.transformed(&f.placement);
        // A face whose surface lies wholly within the box so far cannot widen it.
        if control_box(&surface).is_some_and(|c| b.contains_box(&c)) {
            continue;
        }
        let (points, by) = face_extreme_points(model, f, &surface);
        for &p in points.iter() {
            b.add_point(p + by);
        }
    }
    b
}

/// [`tight_bounding_box`] given the shape's distinct edges and, for each face as
/// `model.explore` lists them, a box known to hold it: a face whose box lies within the box found so far cannot widen it, and its
/// extremes are not sought.
pub(crate) fn tight_bounding_box_within(model: &Model, shape: Shape, edges: &[Shape], face_bounds: &[BoundingBox]) -> BoundingBox {
    let mut b = BoundingBox::EMPTY;
    // The edges' extreme points include their ends (a degenerate edge's, its vertex): the
    // vertices count alone only without edges.
    if edges.is_empty() {
        for v in model.explore(shape, NodeKind::Vertex) {
            b.add_point(model.point(v));
        }
    }
    for &edge in edges {
        let (points, by) = edge_extreme_points(model, edge);
        for &p in points.iter() {
            b.add_point(p + by);
        }
    }
    let faces: Vec<Shape> = model.explore(shape, NodeKind::Face).collect();
    // The faces most likely to widen the box first.
    let mut order: Vec<usize> = (0..faces.len()).collect();
    let size = |k: usize| face_bounds.get(k).map_or(f64::INFINITY, |c| c.min.distance(c.max));
    order.sort_by(|&x, &y| size(y).total_cmp(&size(x)));
    for k in order {
        let f = faces[k];
        let Some(surface) = model.surface(f) else { continue };
        if matches!(surface, Surface::Plane(_)) {
            continue; // bounded by its edges
        }
        if face_bounds.get(k).is_some_and(|c| b.contains_box(c)) {
            continue;
        }
        let surface = surface.transformed(&f.placement);
        if control_box(&surface).is_some_and(|c| b.contains_box(&c)) {
            continue;
        }
        let (points, by) = face_extreme_points(model, f, &surface);
        for &p in points.iter() {
            b.add_point(p + by);
        }
    }
    b
}

/// An edge's points where a coordinate is extreme (its ends included), kept in the model
/// for an unplaced edge: the points, and the translation to add to them.
fn edge_extreme_points(model: &Model, edge: Shape) -> (std::sync::Arc<Vec<Point>>, Vector) {
    struct EdgeExtremes(std::sync::Arc<Vec<Point>>);
    // (Found unplaced, kept with the edge, then moved where it is only moved.)
    let by = edge.placement.translation_part();
    match crate::topology::EdgeId::try_from(edge.node) {
        Ok(id) if crate::math::Transform::translation(by) == edge.placement => {
            let local = Shape { placement: crate::math::Transform::IDENTITY, ..edge };
            (model.derived(edge.node, model.edge_fingerprint(id), || EdgeExtremes(std::sync::Arc::new(edge_extreme_points_of(model, local)))).0.clone(), by)
        }
        _ => (std::sync::Arc::new(edge_extreme_points_of(model, edge)), Vector::ZERO),
    }
}

/// [`edge_extreme_points`], found anew.
fn edge_extreme_points_of(model: &Model, edge: Shape) -> Vec<Point> {
    let Some((curve, first, last)) = model.edge_curve(edge) else {
        let (a, b) = model.edge_vertices(edge);
        return vec![model.point(a), model.point(b)];
    };
    if matches!(curve, crate::geometry::curve::Curve::Line(_)) {
        return vec![curve.point(first), curve.point(last)];
    }
    curve_extremes(&curve, first, last).into_iter().map(|t| curve.point(t)).collect()
}

/// A curved face's points where a coordinate is extreme, kept in the model for an unplaced
/// face; `surface` is the face's, placed.
fn face_extreme_points(model: &Model, face: Shape, surface: &Surface) -> (std::sync::Arc<Vec<Point>>, Vector) {
    struct FaceExtremes(std::sync::Arc<Vec<Point>>);
    // (Found unplaced, kept with the face, then moved where it is only moved.)
    let by = face.placement.translation_part();
    match FaceId::try_from(face.node) {
        Ok(id) if crate::math::Transform::translation(by) == face.placement => {
            let local = Shape { placement: crate::math::Transform::IDENTITY, ..face };
            let kept = model
                .derived(face.node, model.face_fingerprint(id), || {
                    let unplaced = if by == Vector::ZERO { surface.clone() } else { model.surface(local).cloned().unwrap_or_else(|| surface.transformed(&crate::math::Transform::translation(-by))) };
                    FaceExtremes(std::sync::Arc::new(face_extremes(model, local, &unplaced)))
                })
                .0
                .clone();
            (kept, by)
        }
        _ => (std::sync::Arc::new(face_extremes(model, face, surface)), Vector::ZERO),
    }
}

/// The box of a free-form surface's poles, which holds the whole surface (its weights being
/// positive); none for other surfaces.
pub(crate) fn control_box(surface: &Surface) -> Option<BoundingBox> {
    let (poles, weights) = match surface {
        Surface::BSpline(s) => (s.poles(), s.weights()),
        Surface::Bezier(s) => (s.poles(), s.weights()),
        _ => return None,
    };
    if weights.is_some_and(|w| w.iter().flatten().any(|&w| w <= 0.0)) {
        return None;
    }
    Some(BoundingBox::from_points(poles.iter().flatten().copied()))
}

/// Parameters in `[first, last]` where each coordinate of `curve` is extreme: the ends, and
/// the stationary points found from samples by Newton steps.
fn curve_extremes(curve: &crate::geometry::curve::Curve, first: f64, last: f64) -> Vec<f64> {
    let (lo, hi) = (first.min(last), first.max(last));
    let n = 64;
    let ts: Vec<f64> = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).collect();
    let points: Vec<Point> = ts.iter().map(|&t| curve.point(t)).collect();
    let mut found = vec![lo, hi];
    for k in 0..3 {
        let c = |p: Point| [p.x, p.y, p.z][k];
        for i in 1..n {
            let (a, m, z) = (c(points[i - 1]), c(points[i]), c(points[i + 1]));
            if (m - a) * (z - m) > 0.0 {
                continue; // monotonic here
            }
            // Newton on the coordinate's derivative.
            let mut t = ts[i];
            for _ in 0..30 {
                let (_, d1, d2) = curve.d2(t);
                let (g, h) = ([d1.x, d1.y, d1.z][k], [d2.x, d2.y, d2.z][k]);
                if h.abs() <= 1e-300 {
                    break;
                }
                let next = (t - g / h).clamp(ts[i - 1], ts[i + 1]);
                let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
                t = next;
                if done {
                    break;
                }
            }
            found.push(t);
        }
    }
    found
}

/// Points of a curved face where a coordinate is extreme inside its boundary (its extremes on
/// the boundary are its edges'): from a grid over the boundary's parameter range, refined by
/// Newton steps that stay inside.
fn face_extremes(model: &Model, face: Shape, surface: &Surface) -> Vec<Point> {
    // On a surface straight along v (a cylinder, a cone, an extrusion) each coordinate is
    // linear along a v-line, and the v-line through a point of the face meets the face's
    // boundary at its ends: the face's extremes are its boundary's, which its edges give.
    if matches!(surface, Surface::Cylinder(_) | Surface::Cone(_) | Surface::Extrusion(_)) {
        return Vec::new();
    }
    // The boundary in parameter space.
    let mut loops: Vec<Vec<[f64; 2]>> = Vec::new();
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    let (u0, u1, v0, v1) = surface.bounds();
    let clamp_v = |v: f64| v.clamp(-1e6, 1e6);
    for wire in model.children(Shape { orientation: Orientation::Forward, ..face }) {
        let mut points: Vec<[f64; 2]> = Vec::new();
        for coedge in model.children(wire) {
            let Ok(id) = crate::topology::CoEdgeId::try_from(coedge.node) else { continue };
            let Some((pcurve, a, c)) = model.effective_pcurve(id) else { continue };
            // Along the loop's direction (the coedge's own orientation in the wire, the
            // wire's on the face).
            let reversed = model[id].orientation.compose(wire.orientation) == Orientation::Reversed;
            let (a, c) = if reversed { (c, a) } else { (a, c) };
            for i in 0..32 {
                let uv = pcurve.point(a + (c - a) * i as f64 / 32.0);
                points.push([uv.x, uv.y]);
            }
        }
        // A loop running round a whole period without a seam (a sphere's cap bounded by
        // one circle) encloses the face with the surface's bound on its left.
        if let (Some(period), Some(&first), Some(&last)) = (surface.u_period(), points.first(), points.last()) {
            let travel = last[0] - first[0];
            if (travel.abs() - period).abs() <= 1e-6 * period && (last[1] - first[1]).abs() <= 1e-6 {
                let bound = if travel > 0.0 { clamp_v(v1) } else { clamp_v(v0) };
                points.push([last[0], bound]);
                points.push([first[0], bound]);
            }
        }
        for p in &points {
            lo = [lo[0].min(p[0]), lo[1].min(p[1])];
            hi = [hi[0].max(p[0]), hi[1].max(p[1])];
        }
        if !points.is_empty() {
            loops.push(points);
        }
    }
    let _ = (u0, u1);
    if lo[0] > hi[0] {
        return Vec::new();
    }
    let inside = |uv: [f64; 2]| {
        let mut crossings = 0;
        for points in &loops {
            for i in 0..points.len() {
                let (a, b) = (points[i], points[(i + 1) % points.len()]);
                if (a[1] > uv[1]) != (b[1] > uv[1]) && uv[0] < a[0] + (uv[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    };
    // (On a sphere a coordinate is stationary only at the points of its axes: those the face
    // has, found in closed form below, and its boundary's, its edges' extremes. No search.)
    let sphere = matches!(surface, Surface::Sphere(_));
    let n = if sphere { 0 } else { 24 };
    let grid: Vec<[f64; 2]> = (0..=n)
        .flat_map(|i| (0..=n).map(move |j| (i, j)))
        .map(|(i, j)| [lo[0] + (hi[0] - lo[0]) * i as f64 / n as f64, lo[1] + (hi[1] - lo[1]) * j as f64 / n as f64])
        .filter(|&uv| !sphere && inside(uv))
        .collect();
    if grid.is_empty() && !sphere {
        return Vec::new();
    }
    // The grid's points, found once for all six directions.
    let points: Vec<Point> = grid.iter().map(|uv| surface.point(uv[0], uv[1])).collect();
    let mut found = Vec::new();
    for k in 0..3 {
        let c = |p: Point| [p.x, p.y, p.z][k];
        for sign in [1.0, -1.0] {
            if grid.is_empty() {
                continue;
            }
            // The best grid point, then Newton on the coordinate's gradient.
            let best = (0..grid.len()).max_by(|&a, &b| (sign * c(points[a])).total_cmp(&(sign * c(points[b])))).expect("grid");
            let mut uv = grid[best];
            for _ in 0..20 {
                let d = surface.d2(uv[0], uv[1]);
                let g = [[d.du.x, d.du.y, d.du.z][k], [d.dv.x, d.dv.y, d.dv.z][k]];
                let h = [[d.duu.x, d.duu.y, d.duu.z][k], [d.duv.x, d.duv.y, d.duv.z][k], [d.dvv.x, d.dvv.y, d.dvv.z][k]];
                let det = h[0] * h[2] - h[1] * h[1];
                if det.abs() <= 1e-300 {
                    break;
                }
                let step = [-(g[0] * h[2] - g[1] * h[1]) / det, -(g[1] * h[0] - g[0] * h[1]) / det];
                let next = [uv[0] + step[0], uv[1] + step[1]];
                if !inside(next) || sign * c(surface.point(next[0], next[1])) < sign * c(surface.point(uv[0], uv[1])) {
                    break;
                }
                uv = next;
                if step[0].abs() + step[1].abs() <= 1e-14 {
                    break;
                }
            }
            found.push(surface.point(uv[0], uv[1]));
        }
    }
    // A sphere's extremes in closed form (often on the border of its parameters, at a
    // pole or the seam, where the search above stops short), where the face has them.
    if let Surface::Sphere(sphere) = surface {
        let near_border = |uv: [f64; 2]| {
            loops.iter().any(|points| {
                (0..points.len()).any(|i| {
                    let (a, b) = (points[i], points[(i + 1) % points.len()]);
                    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                    let len2 = dx * dx + dy * dy;
                    let t = if len2 > 0.0 { (((uv[0] - a[0]) * dx + (uv[1] - a[1]) * dy) / len2).clamp(0.0, 1.0) } else { 0.0 };
                    (uv[0] - a[0] - t * dx).hypot(uv[1] - a[1] - t * dy) <= 1e-9
                })
            })
        };
        let c = sphere.frame().origin();
        for k in 0..3 {
            for sign in [1.0, -1.0] {
                let mut e = [0.0; 3];
                e[k] = sign * sphere.radius();
                let p = c + Vector::new(e[0], e[1], e[2]);
                let Some((u, v)) = surface.parameters(p) else { continue };
                let at_pole = (v.abs() - std::f64::consts::FRAC_PI_2).abs() <= 1e-9;
                let tau = std::f64::consts::TAU;
                let us: Vec<f64> = if at_pole { vec![0.5 * (lo[0] + hi[0])] } else { vec![u, u - tau, u + tau] };
                if us.iter().any(|&u| inside([u, v]) || near_border([u, v])) {
                    found.push(p);
                }
            }
        }
    }
    found
}

/// Mass properties of point masses `(point, mass)`; the mass of an empty set is 0 (with
/// the origin as center).
#[doc(alias = "GProp_PGProps")]
pub fn point_masses(points: &[(Point, f64)]) -> MassProperties {
    let mass: f64 = points.iter().map(|p| p.1).sum();
    if points.is_empty() || mass == 0.0 {
        return MassProperties { mass, center_of_mass: Point::ORIGIN, inertia: Matrix3::ZERO };
    }
    let center = points.iter().fold(Vector::ZERO, |acc, &(p, m)| acc + (p - Point::ORIGIN) * m) * (1.0 / mass);
    let center_of_mass = Point::ORIGIN + center;
    // About the center: Σ m (|r|² I - r rᵀ).
    let mut inertia = Matrix3::ZERO;
    for &(p, m) in points {
        let r = p - center_of_mass;
        let c = [r.x, r.y, r.z];
        let r2 = r.dot(r);
        for i in 0..3 {
            for j in 0..3 {
                inertia[(i, j)] += m * (if i == j { r2 } else { 0.0 } - c[i] * c[j]);
            }
        }
    }
    MassProperties { mass, center_of_mass, inertia }
}

/// Mass properties of unit point masses.
#[doc(alias = "GProp_PGProps")]
pub fn points(points: &[Point]) -> MassProperties {
    point_masses(&points.iter().map(|&p| (p, 1.0)).collect::<Vec<_>>())
}

/// The average of `points`; panics if there are none.
#[doc(alias = "GProp_PGProps::Barycentre")]
pub fn barycenter(points: &[Point]) -> Point {
    assert!(!points.is_empty(), "the barycenter of no points");
    self::points(points).center_of_mass
}

/// The smallest-dimensional shape that holds a point set within a tolerance.
#[doc(alias = "GProp_PEquation")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PointSetShape {
    /// All points coincide.
    Point(Point),
    /// All points are on a line.
    Line(Line),
    /// All points are on a plane.
    Plane(Plane),
    /// The points span space; they fit in the box at `corner` with edges `edges`
    /// (mutually perpendicular, along the principal axes of the set).
    Space { corner: Point, edges: [Vector; 3] },
}

/// Classifies `points` by the dimension of the space they span within `tolerance`.
/// Panics if there are no points.
#[doc(alias = "GProp_PEquation")]
pub fn point_set_shape(points: &[Point], tolerance: f64) -> PointSetShape {
    assert!(!points.is_empty(), "the shape of no points");
    let props = self::points(points);
    let g = props.center_of_mass;
    // The extents along the principal axes of the set, smallest moment first (the axis
    // of largest spread).
    let (_, axes) = symmetric_eigen(&props.inertia);
    let extents: Vec<(f64, f64, Vector)> = axes
        .iter()
        .map(|&d| {
            let (lo, hi) = points.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &p| {
                let t = (p - g).dot(d);
                (lo.min(t), hi.max(t))
            });
            (lo, hi, d)
        })
        .collect();
    let flat: Vec<bool> = extents.iter().map(|e| e.1 - e.0 <= tolerance).collect();
    match flat.iter().filter(|&&f| f).count() {
        3 => PointSetShape::Point(g),
        2 => {
            let k = flat.iter().position(|&f| !f).expect("one spread axis");
            PointSetShape::Line(Line::new(g, Direction::from_vector(extents[k].2)))
        }
        1 => {
            let k = flat.iter().position(|&f| f).expect("one flat axis");
            PointSetShape::Plane(Plane::new(g, Direction::from_vector(extents[k].2)))
        }
        _ => {
            let corner = extents.iter().fold(g, |c, e| c + e.2 * e.0);
            let edges = [0, 1, 2].map(|k| extents[k].2 * (extents[k].1 - extents[k].0));
            PointSetShape::Space { corner, edges }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Direction;
    use crate::primitives::Cuboid;
    use std::f64::consts::TAU;

    /// The area of a helicoid strip winding twenty times in 63 spans (a B-spline ruled
    /// between two helices, as a thread's flank): its Gauss rules each spanning only a little
    /// of a turn, as a fine sum over the surface finds.
    #[test]
    fn areas_of_surfaces_winding_many_times() {
        use crate::geometry::surface::{BSplineSurface, KnotVector};
        use crate::geometry::ParametricSurface;
        let (turns, n) = (20.0, 64);
        let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let helix = |r: f64| {
            let points: Vec<Point> = params.iter().map(|&t| Point::new(r * (TAU * turns * t).cos(), r * (TAU * turns * t).sin(), turns * t)).collect();
            crate::geometry::fit::interpolate_with_parameters(&points, &params).unwrap()
        };
        let (inner, outer) = (helix(2.0), helix(3.0));
        let poles: Vec<Vec<Point>> = inner.poles().iter().zip(outer.poles()).map(|(&a, &b)| vec![a, b]).collect();
        let u = KnotVector::new(inner.knots().to_vec(), inner.multiplicities().to_vec(), inner.degree());
        let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let strip = BSplineSurface::new(poles, u, v);
        let (nu, nv) = (20000, 8);
        let fine: f64 = (0..nu)
            .flat_map(|i| (0..nv).map(move |j| ((i as f64 + 0.5) / nu as f64, (j as f64 + 0.5) / nv as f64)))
            .map(|(u, v)| {
                let (_, du, dv) = strip.d1(u, v);
                du.cross(dv).length()
            })
            .sum::<f64>()
            / (nu * nv) as f64;
        let mut model = Model::new();
        let face = crate::make::face_on_surface(&mut model, crate::geometry::surface::Surface::BSpline(strip)).unwrap();
        let a = area(&model, face).mass;
        assert!((a - fine).abs() <= 1e-5 * fine, "area {a} != {fine}");
    }

    #[test]
    fn oriented_boxes_fit_turned_blocks() {
        let mut model = Model::new();
        let turn = |s: Shape, axis: Direction, angle: f64| s.moved(&crate::math::Transform::rotation(crate::math::Axis::new(Point::ORIGIN, axis), angle));
        let block = Cuboid::new(100.0, 150.0, 200.0).build(&mut model).unwrap();
        let a = turn(block, Direction::new(1.0, 2.0, 3.0), 0.7);
        let b = oriented_bounding_box(&model, a);
        let volume = 8.0 * b.half_sizes.iter().product::<f64>();
        assert!((volume - 100.0 * 150.0 * 200.0).abs() < 1e-6 * volume, "{volume}");
        // Apart though their axis-aligned boxes overlap; and overlapping.
        let thin = Cuboid::new(300.0, 10.0, 10.0).build(&mut model).unwrap();
        let (c, d) = (turn(thin, Direction::Z, 0.785), turn(thin, Direction::Z, 0.785).moved(&crate::math::Transform::translation(crate::math::Vector::new(20.0, -20.0, 0.0))));
        assert!(!oriented_bounding_box(&model, c).intersects(&oriented_bounding_box(&model, d)));
        assert!(oriented_bounding_box(&model, c).intersects(&oriented_bounding_box(&model, c)));
    }

    #[test]
    fn tight_boxes_of_tilted_cylinders_and_cones_hold_their_sides() {
        let mut model = Model::new();
        let tilt = crate::math::Frame::new(Point::new(1.0, -2.0, 0.5), Direction::new(0.3, -0.5, 0.8), Direction::new(0.9, 0.4, -0.1));
        let cylinder = crate::primitives::Cylinder::new(1.5, 4.0).with_frame(tilt).build(&mut model).unwrap();
        let cone = crate::primitives::Cone::new(2.0, 0.5, 3.0).with_frame(tilt).build(&mut model).unwrap();
        for solid in [cylinder, cone] {
            let b = tight_bounding_box(&model, solid);
            // The sides' points, densely: within the box, and reaching its faces.
            let mut sampled = BoundingBox::EMPTY;
            for face in model.explore(solid, NodeKind::Face).collect::<Vec<_>>() {
                let surface = model.surface(face).unwrap().transformed(&face.placement);
                if matches!(surface, Surface::Plane(_)) {
                    continue;
                }
                // Over the face's parameter range, from its boundary.
                let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
                for wire in model.children(face).collect::<Vec<_>>() {
                    for coedge in model.children(wire).collect::<Vec<_>>() {
                        let id = crate::topology::CoEdgeId::try_from(coedge.node).unwrap();
                        let (pcurve, a, c) = model.effective_pcurve(id).unwrap();
                        for k in 0..=50 {
                            let q = pcurve.point(a + (c - a) * k as f64 / 50.0);
                            lo = [lo[0].min(q.x), lo[1].min(q.y)];
                            hi = [hi[0].max(q.x), hi[1].max(q.y)];
                        }
                    }
                }
                for i in 0..=400 {
                    for j in 0..=40 {
                        let (u, v) = (lo[0] + (hi[0] - lo[0]) * i as f64 / 400.0, lo[1] + (hi[1] - lo[1]) * j as f64 / 40.0);
                        sampled.add_point(surface.point(u, v));
                    }
                }
            }
            for (got, near) in [(b.min.x, sampled.min.x), (b.min.y, sampled.min.y), (b.min.z, sampled.min.z)] {
                assert!(got <= near + 1e-9 && near - got <= 1e-3, "{got} against {near}");
            }
            for (got, near) in [(b.max.x, sampled.max.x), (b.max.y, sampled.max.y), (b.max.z, sampled.max.z)] {
                assert!(got >= near - 1e-9 && got - near <= 1e-3, "{got} against {near}");
            }
        }
    }

    #[test]
    fn face_volume_terms_add_up_to_the_volume() {
        let mut model = Model::new();
        let frame = crate::math::Frame::new(Point::new(3.0, -2.0, 5.0), Direction::Z, Direction::X);
        let cylinder = crate::primitives::Cylinder::new(1.5, 4.0).with_frame(frame).build(&mut model).unwrap();
        let cuboid = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let holed = crate::boolean::difference(&mut model, cuboid, cylinder).unwrap();
        for shape in [cylinder, cuboid, holed] {
            let terms: f64 = model.explore(shape, NodeKind::Face).collect::<Vec<_>>().into_iter().map(|f| face_volume_term(&model, f)).sum();
            let exact = volume(&model, shape);
            let volume = exact.mass;
            assert!((terms - volume).abs() <= 1e-9 * volume.abs().max(1.0), "{terms} against {volume}");
            let (v, center) = solid_volume(&model, shape).expect("solids");
            assert!((v - volume).abs() <= 1e-9 * volume.abs().max(1.0));
            assert!(center.distance(exact.center_of_mass) <= 1e-9, "{center:?} against {:?}", exact.center_of_mass);
            // Kept: the same again, and a reversed face counts negatively.
            let face = model.explore(shape, NodeKind::Face).next().unwrap();
            assert_eq!(face_volume_term(&model, face), face_volume_term(&model, face));
            assert_eq!(face_volume_term(&model, face.reversed()), -face_volume_term(&model, face));
        }
    }

    fn near(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * b.abs().max(1.0)
    }

    #[test]
    fn box_mass_properties() {
        let mut model = Model::new();
        let frame = crate::math::Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::new(1.0, 1.0, 0.0),
            Direction::new(0.0, 0.0, 1.0),
        );
        let (a, b, c) = (3.0, 4.0, 5.0);
        let shape = Cuboid::new(a, b, c).with_frame(frame).build(&mut model).unwrap();
        let props = volume(&model, shape);
        assert!(near(props.mass, 60.0), "{props:?}");

        let expected_center = frame.origin()
            + Vector::from(frame.x()) * (a / 2.0)
            + Vector::from(frame.y()) * (b / 2.0)
            + Vector::from(frame.z()) * (c / 2.0);
        assert!(props.center_of_mass.distance(expected_center) <= 1e-9, "{props:?}");

        // Inertia about the frame's X axis through the center: m (b² + c²) / 12.
        let x = Vector::from(frame.x());
        let ix = (0..3).map(|r| (0..3).map(|k| x[r] * props.inertia[(r, k)] * x[k]).sum::<f64>()).sum::<f64>();
        assert!(near(ix, 60.0 * (b * b + c * c) / 12.0), "{ix}");
    }

    #[test]
    fn box_from_corners_and_reversed_solid() {
        let mut model = Model::new();
        let shape = Cuboid::from_corners(Point::new(4.0, 6.0, 8.0), Point::new(1.0, 2.0, 3.0))
            .build(&mut model)
            .unwrap();
        assert!(near(volume(&model, shape).mass, 60.0));
        let mut reversed = shape;
        reversed.orientation = Orientation::Reversed;
        assert!(near(volume(&model, reversed).mass, -60.0));
    }

    #[test]
    fn box_area_and_shared_faces() {
        let mut model = Model::new();
        let shape = Cuboid::new(3.0, 4.0, 5.0).build(&mut model).unwrap();
        let props = area(&model, shape);
        assert!(near(props.mass, 2.0 * (12.0 + 15.0 + 20.0)), "{props:?}");
        assert!(props.center_of_mass.distance(Point::new(1.5, 2.0, 2.5)) <= 1e-9, "{props:?}");

        // Each box face is reached once, so skipping shared faces changes nothing.
        let options = Options { skip_shared: true, use_meshes: false };
        assert!(near(area_with(&model, shape, options).mass, props.mass));
    }

    fn inertia_about(props: &MassProperties, axis: Vector) -> f64 {
        (0..3).map(|r| (0..3).map(|c| axis[r] * props.inertia[(r, c)] * axis[c]).sum::<f64>()).sum()
    }

    #[test]
    fn cylinder_properties() {
        use crate::primitives::Cylinder;
        use std::f64::consts::PI;
        let (r, h) = (2.0, 7.0);
        let mut model = Model::new();
        let full = Cylinder::new(r, h).build(&mut model).unwrap();
        let props = volume(&model, full);
        let mass = PI * r * r * h;
        assert!(near(props.mass, mass), "{props:?}");
        assert!(props.center_of_mass.distance(Point::new(0.0, 0.0, h / 2.0)) <= 1e-9, "{props:?}");
        assert!(near(inertia_about(&props, Vector::new(0.0, 0.0, 1.0)), mass * r * r / 2.0));
        assert!(near(inertia_about(&props, Vector::new(1.0, 0.0, 0.0)), mass * (3.0 * r * r + h * h) / 12.0));
        assert!(near(area(&model, full).mass, 2.0 * PI * r * h + 2.0 * PI * r * r));

        // Half cylinder: half the volume, centroid 4r/(3π) off the axis towards +y.
        let half = Cylinder::new(r, h).with_angle(PI).build(&mut model).unwrap();
        let props = volume(&model, half);
        assert!(near(props.mass, mass / 2.0), "{props:?}");
        let expected = Point::new(0.0, 4.0 * r / (3.0 * PI), h / 2.0);
        assert!(props.center_of_mass.distance(expected) <= 1e-9, "{props:?}");
        assert!(near(area(&model, half).mass, PI * r * h + PI * r * r + 2.0 * r * h));
    }

    /// A compound of one cylinder placed three ways (moved, and turned two ways and moved): its
    /// volume, centre, area and box those of the three cylinders, the placed copies measured
    /// from the unplaced one.
    #[test]
    fn copies_of_a_cylinder_placed_in_a_compound() {
        use crate::math::{Axis, Transform};
        use crate::primitives::Cylinder;
        use std::f64::consts::PI;
        let (r, h) = (2.0, 7.0);
        let mut model = Model::new();
        let cylinder = Cylinder::new(r, h).build(&mut model).unwrap();
        let x = Direction::from_vector(Vector::new(1.0, 0.0, 0.0));
        let placements = [
            Transform::translation(Vector::new(0.0, 20.0, 0.0)),
            Transform::translation(Vector::new(10.0, 0.0, 0.0)) * Transform::rotation(Axis::new(Point::ORIGIN, x), PI / 2.0),
            Transform::translation(Vector::new(0.0, 0.0, -17.0)) * Transform::rotation(Axis::new(Point::ORIGIN, Direction::from_vector(Vector::new(0.0, 0.0, 1.0))), PI / 3.0),
        ];
        let compound = model.new_empty(NodeKind::Compound);
        for t in &placements {
            model.add(compound, cylinder.moved(t)).unwrap();
        }
        let props = volume(&model, compound);
        assert!(near(props.mass, 3.0 * PI * r * r * h), "{props:?}");
        let centre = Point::new(10.0 / 3.0, (20.0 - 3.5) / 3.0, (3.5 - 13.5) / 3.0);
        assert!(props.center_of_mass.distance(centre) <= 1e-9, "{props:?}");
        let props = area(&model, compound);
        assert!(near(props.mass, 3.0 * (2.0 * PI * r * h + 2.0 * PI * r * r)), "{props:?}");
        assert!(props.center_of_mass.distance(centre) <= 1e-9, "{props:?}");
        let b = tight_bounding_box(&model, compound);
        assert!(b.min.distance(Point::new(-2.0, -7.0, -17.0)) <= 1e-9 && b.max.distance(Point::new(12.0, 22.0, 7.0)) <= 1e-9, "{b:?}");
    }

    #[test]
    fn sphere_properties() {
        use crate::primitives::Sphere;
        use std::f64::consts::PI;
        let r = 3.0;
        let mut model = Model::new();
        let full = Sphere::new(r).build(&mut model).unwrap();
        let props = volume(&model, full);
        let mass = 4.0 / 3.0 * PI * r * r * r;
        assert!(near(props.mass, mass), "{props:?}");
        assert!(props.center_of_mass.distance(Point::ORIGIN) <= 1e-9, "{props:?}");
        for axis in [Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)] {
            assert!(near(inertia_about(&props, axis), 0.4 * mass * r * r));
        }
        assert!(near(area(&model, full).mass, 4.0 * PI * r * r));

        // Quarter sphere (a lune of π/2): its centroid is 3r/8 from the center along the
        // bisector of the lune, scaled by the lune's centroid factor sin(α/2)/(α/2).
        let alpha = PI / 2.0;
        let lune = Sphere::new(r).with_angle(alpha).build(&mut model).unwrap();
        let props = volume(&model, lune);
        assert!(near(props.mass, mass * alpha / (2.0 * PI)), "{props:?}");
        let distance = 3.0 * PI * r / 16.0 * (alpha / 2.0).sin() / (alpha / 2.0);
        let bisector = Vector::new((alpha / 2.0).cos(), (alpha / 2.0).sin(), 0.0);
        assert!(props.center_of_mass.distance(Point::ORIGIN + bisector * distance) <= 1e-9, "{props:?}");
        assert!(near(area(&model, lune).mass, 4.0 * PI * r * r * alpha / (2.0 * PI) + PI * r * r));
    }

    #[test]
    fn cone_properties() {
        use crate::primitives::Cone;
        use std::f64::consts::PI;
        let h = 6.0;
        let mut model = Model::new();
        // Widening, narrowing and pointed frustums (an apex at either end).
        for (r1, r2) in [(1.0, 3.0), (3.0, 1.0), (2.5, 0.0), (0.0, 2.5)] {
            let shape = Cone::new(r1, r2, h).build(&mut model).unwrap();
            let props = volume(&model, shape);
            let sum = r1 * r1 + r1 * r2 + r2 * r2;
            let mass = PI * h / 3.0 * sum;
            assert!(near(props.mass, mass), "{r1} {r2}: {props:?}");
            let height = h * (r1 * r1 + 2.0 * r1 * r2 + 3.0 * r2 * r2) / (4.0 * sum);
            assert!(props.center_of_mass.distance(Point::new(0.0, 0.0, height)) <= 1e-9, "{r1} {r2}: {props:?}");
            let axial = 0.3 * mass * (r1.powi(5) - r2.powi(5)) / (r1.powi(3) - r2.powi(3));
            assert!(near(inertia_about(&props, Vector::new(0.0, 0.0, 1.0)), axial), "{r1} {r2}");
            let slant = (h * h + (r1 - r2) * (r1 - r2)).sqrt();
            let expected_area = PI * (r1 + r2) * slant + PI * (r1 * r1 + r2 * r2);
            assert!(near(area(&model, shape).mass, expected_area), "{r1} {r2}");
        }

        // A quarter of a pointed cone: its centroid lies on the bisector of the sector at
        // (2/3) sin(α/2) ∫ρ³dz / ((α/2) ∫ρ²dz) from the axis.
        let (r, alpha) = (2.0, PI / 2.0);
        let shape = Cone::new(r, 0.0, h).with_angle(alpha).build(&mut model).unwrap();
        let props = volume(&model, shape);
        assert!(near(props.mass, PI * r * r * h / 3.0 * alpha / TAU), "{props:?}");
        let off_axis = (2.0 / 3.0) * (alpha / 2.0).sin() * (r.powi(3) * h / 4.0) / ((alpha / 2.0) * (r * r * h / 3.0));
        let bisector = Vector::new((alpha / 2.0).cos(), (alpha / 2.0).sin(), 0.0);
        let expected = Point::new(0.0, 0.0, h / 4.0) + bisector * off_axis;
        assert!(props.center_of_mass.distance(expected) <= 1e-9, "{props:?}");
    }

    #[test]
    fn sphere_zone_properties() {
        use crate::primitives::Sphere;
        use std::f64::consts::{FRAC_PI_2, PI};
        let r = 2.0;
        let mut model = Model::new();
        for (v1, v2) in [(0.0, FRAC_PI_2), (-0.4, 0.9), (-FRAC_PI_2, -0.2)] {
            let shape = Sphere::new(r).with_latitudes(v1, v2).build(&mut model).unwrap();
            let (z1, z2) = (r * f64::sin(v1), r * f64::sin(v2));
            let mass = PI * ((r * r * z2 - z2.powi(3) / 3.0) - (r * r * z1 - z1.powi(3) / 3.0));
            let moment = PI * ((r * r * z2 * z2 / 2.0 - z2.powi(4) / 4.0) - (r * r * z1 * z1 / 2.0 - z1.powi(4) / 4.0));
            let props = volume(&model, shape);
            assert!(near(props.mass, mass), "{v1} {v2}: {props:?}");
            assert!(props.center_of_mass.distance(Point::new(0.0, 0.0, moment / mass)) <= 1e-9, "{v1} {v2}: {props:?}");
            let cap = |z: f64| if (z.abs() - r).abs() <= 1e-12 { 0.0 } else { PI * (r * r - z * z) };
            let expected_area = 2.0 * PI * r * (z2 - z1) + cap(z1) + cap(z2);
            assert!(near(area(&model, shape).mass, expected_area), "{v1} {v2}");

            // A third of the zone has a third of its volume.
            let wedge = Sphere::new(r).with_latitudes(v1, v2).with_angle(TAU / 3.0).build(&mut model).unwrap();
            assert!(near(volume(&model, wedge).mass, mass / 3.0), "{v1} {v2}");
        }
    }

    #[test]
    fn non_positive_box_is_rejected() {
        let mut model = Model::new();
        assert!(Cuboid::new(1.0, 0.0, 1.0).build(&mut model).is_err());
    }
}

/// The eigenvalues of a symmetric matrix, ascending, with unit eigenvectors (Jacobi
/// rotations).
fn symmetric_eigen(m: &Matrix3) -> ([f64; 3], [Vector; 3]) {
    let mut a = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = 0.5 * (m[(i, j)] + m[(j, i)]);
        }
    }
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..100 {
        // The largest off-diagonal entry.
        let (mut p, mut q) = (0, 1);
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            if a[i][j].abs() > a[p][q].abs() {
                (p, q) = (i, j);
            }
        }
        let scale = a.iter().flatten().fold(0.0f64, |s, x| s.max(x.abs()));
        if a[p][q].abs() <= 1e-15 * scale.max(1e-300) {
            break;
        }
        let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
        let t = if theta == 0.0 { 1.0 } else { t };
        let c = 1.0 / (t * t + 1.0).sqrt();
        let s = t * c;
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
        for row in &mut v {
            let (vp, vq) = (row[p], row[q]);
            row[p] = c * vp - s * vq;
            row[q] = s * vp + c * vq;
        }
    }
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.map(|i| a[i][i]);
    let vectors = order.map(|i| Vector::new(v[0][i], v[1][i], v[2][i]));
    (values, vectors)
}
