//! What CadQuery's layer needs beyond the basic bindings: curves and faces read as its
//! OCCT adaptors read them, the edges, solids and sweeps it makes, and shapes split by one
//! another.

use pyo3::prelude::*;

use kentron::geometry::curve::{BezierCurve, Curve, Ellipse};
use kentron::geometry::extrema::nearest_on_curve;
use kentron::geometry::fit;
use kentron::geometry::measure;
use kentron::geometry::surface::{Plane, Surface};
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Point, Vector};
use kentron::topology::{Model, NodeKind, Orientation, Shape};
use kentron::{boolean, make, primitives, properties, sweep};

use super::{direction, frame, handle, nearest_parameters, point, vector, with_model, xyz, Handle, Placed, Xyz};

// --- curves as adaptors see them ---------------------------------------------------------

/// One edge of a [`Path`]: its curve over `[first, last]`, run backwards if `reversed`,
/// starting at path parameter `start`.
struct Piece {
    curve: Curve,
    first: f64,
    last: f64,
    reversed: bool,
    start: f64,
    /// Lengths along the curve over its range (kept with the edge), and the factor the
    /// edge's placement scales them by.
    lengths: Option<(std::sync::Arc<kentron::geometry::measure::LengthTable>, f64)>,
}

impl Piece {
    fn span(&self) -> f64 {
        self.last - self.first
    }

    /// The curve parameter at path parameter `u`.
    fn local(&self, u: f64) -> f64 {
        if self.reversed { self.last - (u - self.start) } else { self.first + (u - self.start) }
    }

    /// The path parameter at curve parameter `t`.
    fn global(&self, t: f64) -> f64 {
        if self.reversed { self.start + (self.last - t) } else { self.start + (t - self.first) }
    }
}

/// A curve as CadQuery's adaptors see it: an edge's own curve over its range, whatever the
/// edge's orientation (BRepAdaptor_Curve), or a wire's edges in order, each run along the
/// wire, their ranges laid end to end from 0 (BRepAdaptor_CompCurve).
struct Path {
    pieces: Vec<Piece>,
    first: f64,
    last: f64,
}

fn path(m: &Model, s: Shape) -> kentron::Result<Path> {
    let invalid = |msg: &str| kentron::Error::InvalidInput(msg.into());
    match s.kind() {
        NodeKind::Edge => {
            let (curve, first, last) = m.edge_curve(s).ok_or_else(|| invalid("a degenerate edge"))?;
            Ok(Path { pieces: vec![Piece { curve, first, last, reversed: false, start: first, lengths: lengths(m, s) }], first, last })
        }
        NodeKind::Wire => {
            let mut pieces = Vec::new();
            let mut at = 0.0;
            for e in m.ordered_edges(s) {
                let Some((curve, first, last)) = m.edge_curve(e) else { continue };
                let piece = Piece { curve, first, last, reversed: e.orientation == Orientation::Reversed, start: at, lengths: lengths(m, e) };
                at += piece.span();
                pieces.push(piece);
            }
            if pieces.is_empty() {
                return Err(invalid("a wire without curves"));
            }
            Ok(Path { pieces, first: 0.0, last: at })
        }
        _ => Err(invalid("an edge or wire is needed")),
    }
}

/// An edge's table of lengths (see [`Piece::lengths`]).
fn lengths(m: &Model, edge: Shape) -> Option<(std::sync::Arc<kentron::geometry::measure::LengthTable>, f64)> {
    Some((kentron::properties::edge_length_table(m, edge)?, edge.placement.scale_factor().abs()))
}

impl Piece {
    /// The piece's whole length.
    fn length(&self) -> f64 {
        match &self.lengths {
            Some((table, scale)) => table.total() * scale,
            None => measure::length_between(&self.curve, self.first, self.last),
        }
    }

    /// The curve parameter `distance` along the piece from its start (as the path runs it).
    fn parameter_at_length(&self, distance: f64) -> f64 {
        match &self.lengths {
            Some((table, scale)) => {
                let along = if self.reversed { self.length() - distance } else { distance };
                table.parameter_at(&self.curve, along / scale)
            }
            None if self.reversed => measure::parameter_at_length(&self.curve, self.last, -distance),
            None => measure::parameter_at_length(&self.curve, self.first, distance),
        }
    }
}

impl Path {
    /// The piece holding path parameter `u` (the end ones beyond the ends).
    fn piece(&self, u: f64) -> &Piece {
        self.pieces.iter().find(|p| u <= p.start + p.span()).unwrap_or_else(|| self.pieces.last().expect("a piece"))
    }

    /// The point and first two derivatives at `u`.
    fn d2(&self, u: f64) -> (Point, Vector, Vector) {
        let piece = self.piece(u);
        let (p, d1, d2) = piece.curve.d2(piece.local(u));
        if piece.reversed { (p, -d1, d2) } else { (p, d1, d2) }
    }

    fn length(&self) -> f64 {
        self.pieces.iter().map(Piece::length).sum()
    }

    /// The parameter `distance` along from the start.
    fn parameter_at_length(&self, distance: f64) -> f64 {
        let mut left = distance;
        for (i, p) in self.pieces.iter().enumerate() {
            let length = p.length();
            if left <= length || i + 1 == self.pieces.len() {
                // (At the piece's end, or past it: its end, not searched for.)
                let t = if left >= length * (1.0 - 1e-14) {
                    if p.reversed { p.first } else { p.last }
                } else {
                    p.parameter_at_length(left)
                };
                return p.global(t);
            }
            left -= length;
        }
        self.last
    }

    /// The parameter of the point nearest `q`.
    fn project(&self, q: Point) -> f64 {
        self.pieces
            .iter()
            .map(|p| {
                let found = nearest_on_curve(q, &p.curve, p.first, p.last);
                (found.distance, p.global(found.parameter))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .map_or(self.first, |(_, u)| u)
    }
}

/// (first, last, closed, periodic, period) of an edge's or wire's curve, as CadQuery's
/// adaptors see it.
#[pyfunction]
pub fn curve_info(shape: Handle) -> PyResult<(f64, f64, bool, bool, f64)> {
    with_model(|m| {
        let path = path(m, shape.shape)?;
        let (a, b) = (path.d2(path.first).0, path.d2(path.last).0);
        let closed = a.distance(b) <= kentron::tolerance::LINEAR;
        let (periodic, period) = match &path.pieces[..] {
            [p] if p.curve.is_periodic() => (true, p.curve.period().unwrap_or(0.0)),
            _ => (false, 0.0),
        };
        Ok((path.first, path.last, closed, periodic, period))
    })
}

/// (point, first derivative, second derivative) at each of the curve parameters `us`.
#[pyfunction]
pub fn curve_derivatives(shape: Handle, us: Vec<f64>) -> PyResult<Vec<(Xyz, Xyz, Xyz)>> {
    with_model(|m| {
        let path = path(m, shape.shape)?;
        let v = |v: Vector| (v.x, v.y, v.z);
        Ok(us.iter().map(|&u| path.d2(u)).map(|(p, d1, d2)| (xyz(p), v(d1), v(d2))).collect())
    })
}

/// The curve's length (over its parameter range).
#[pyfunction]
pub fn curve_length(shape: Handle) -> PyResult<f64> {
    with_model(|m| Ok(path(m, shape.shape)?.length()))
}

/// The curve parameters at `fractions` of its length from its start.
#[pyfunction]
pub fn curve_parameters_at_lengths(shape: Handle, fractions: Vec<f64>) -> PyResult<Vec<f64>> {
    with_model(|m| {
        let path = path(m, shape.shape)?;
        let length = path.length();
        Ok(fractions.iter().map(|&f| path.parameter_at_length(f * length)).collect())
    })
}

/// The curve parameters of the points nearest `points`.
#[pyfunction]
pub fn curve_parameters_of_points(shape: Handle, points: Vec<Xyz>) -> PyResult<Vec<f64>> {
    with_model(|m| {
        let path = path(m, shape.shape)?;
        Ok(points.into_iter().map(|q| path.project(point(q))).collect())
    })
}

/// The normal of the plane a curve lies in: the axis of a circle or ellipse, else that of
/// the plane through its points; `None` if it is not planar.
#[pyfunction]
pub fn curve_plane_normal(shape: Handle) -> PyResult<Option<Xyz>> {
    with_model(|m| {
        let path = path(m, shape.shape)?;
        if let [p] = &path.pieces[..] {
            match &p.curve {
                Curve::Circle(c) => return Ok(Some(xyz(Point::ORIGIN + Vector::from(c.frame().z())))),
                Curve::Ellipse(e) => return Ok(Some(xyz(Point::ORIGIN + Vector::from(e.frame().z())))),
                Curve::Line(_) => return Ok(None),
                _ => {}
            }
        }
        // The plane through samples: the normal of largest cross product from the first.
        let n = 64;
        let pts: Vec<Point> = (0..=n).map(|i| path.d2(path.first + (path.last - path.first) * i as f64 / n as f64).0).collect();
        let size = pts.iter().map(|p| p.distance(pts[0])).fold(0.0, f64::max);
        let far = pts.iter().copied().max_by(|a, b| a.distance(pts[0]).total_cmp(&b.distance(pts[0]))).expect("samples");
        let u = far - pts[0];
        let Some(normal) = pts.iter().map(|&p| u.cross(p - pts[0])).max_by(|a, b| a.length().total_cmp(&b.length())) else { return Ok(None) };
        if normal.length() <= 1e-12 * size * size {
            return Ok(None);
        }
        let normal = normal * (1.0 / normal.length());
        let planar = pts.iter().all(|&p| (p - pts[0]).dot(normal).abs() <= 1e-7 * size.max(1.0));
        Ok(planar.then(|| xyz(Point::ORIGIN + normal)))
    })
}

/// (kind, centre, X, Y and Z directions, major and minor radii) of a circular ("CIRCLE",
/// radii equal) or elliptic ("ELLIPSE") edge's curve; `None` for other curves.
#[pyfunction]
pub fn edge_conic(edge: Handle) -> PyResult<Option<(&'static str, Xyz, Xyz, Xyz, Xyz, f64, f64)>> {
    with_model(|m| {
        let d = |d: kentron::math::Direction| xyz(Point::ORIGIN + Vector::from(d));
        Ok(match m.edge_curve(edge.shape) {
            Some((Curve::Circle(c), _, _)) => {
                let f = c.frame();
                Some(("CIRCLE", xyz(c.center()), d(f.x()), d(f.y()), d(f.z()), c.radius(), c.radius()))
            }
            Some((Curve::Ellipse(e), _, _)) => {
                let f = e.frame();
                Some(("ELLIPSE", xyz(e.center()), d(f.x()), d(f.y()), d(f.z()), e.major_radius(), e.minor_radius()))
            }
            _ => None,
        })
    })
}

/// (degree, poles, flat knots, weights if rational, periodic) of the (clamped) B-spline of
/// an edge's curve over its range.
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn edge_bspline(edge: Handle) -> PyResult<(usize, Vec<Xyz>, Vec<f64>, Option<Vec<f64>>, bool)> {
    with_model(|m| {
        let invalid = |msg: &str| kentron::Error::InvalidInput(msg.into());
        let (curve, first, last) = m.edge_curve(edge.shape).ok_or_else(|| invalid("a degenerate edge"))?;
        // (Over the edge's range: clamped, periodic curves too.)
        let trimmed = Curve::Trimmed(kentron::geometry::curve::TrimmedCurve::new(curve, first, last));
        let spline = kentron::geometry::convert::curve_to_bspline(&trimmed).ok_or_else(|| invalid("a curve without a B-spline form"))?;
        Ok((
            spline.degree(),
            spline.poles().iter().map(|&p| xyz(p)).collect(),
            spline.flat_knots(),
            spline.weights().map(|w| w.to_vec()),
            spline.is_periodic(),
        ))
    })
}

/// The B-spline edge with these poles (and `weights`, rational), `knots` of
/// `multiplicities`, `degree`, periodic or not.
#[pyfunction]
#[pyo3(signature = (poles, weights, knots, multiplicities, degree, periodic=false))]
pub fn make_bspline(poles: Vec<Xyz>, weights: Option<Vec<f64>>, knots: Vec<f64>, multiplicities: Vec<usize>, degree: usize, periodic: bool) -> PyResult<Handle> {
    let poles: Vec<Point> = poles.into_iter().map(point).collect();
    with_model(|m| {
        let mut spline = if periodic {
            kentron::geometry::curve::BSplineCurve::try_new_periodic(poles, knots, multiplicities, degree)?
        } else {
            kentron::geometry::curve::BSplineCurve::try_new(poles, knots, multiplicities, degree)?
        };
        if let Some(w) = weights {
            if w.len() != spline.poles().len() || w.iter().any(|&x| !(x > 0.0)) {
                return Err(kentron::Error::InvalidInput("a B-spline needs one positive weight per pole".into()));
            }
            spline = spline.with_weights(w);
        }
        make::edge_from_curve(m, Curve::BSpline(spline))
    })
    .map(handle)
}

// --- edges -------------------------------------------------------------------------------

/// The spline edge through `points`, as OCCT interpolates it (GeomAPI_Interpolate): with
/// `params` (one more, periodic), and tangents (`None` where free; two for the ends only).
#[pyfunction]
#[pyo3(signature = (points, tangents=None, params=None, periodic=false, scale=true, tolerance=1e-6))]
pub fn make_spline(points: Vec<Xyz>, tangents: Option<Vec<Option<Xyz>>>, params: Option<Vec<f64>>, periodic: bool, scale: bool, tolerance: f64) -> PyResult<Handle> {
    let pts: Vec<Point> = points.into_iter().map(point).collect();
    let mut tangents: Vec<Option<Vector>> = tangents.unwrap_or_default().into_iter().map(|t| t.map(vector)).collect();
    if tangents.len() == 2 && pts.len() != 2 {
        let (a, b) = (tangents[0], tangents[1]);
        tangents = vec![None; pts.len()];
        tangents[0] = a;
        *tangents.last_mut().expect("points") = b;
    }
    with_model(|m| {
        let curve = fit::interpolate_with_tangents(&pts, params.as_deref(), &tangents, periodic, scale, tolerance)?;
        make::edge_from_curve(m, Curve::BSpline(curve))
    })
    .map(handle)
}

/// A spline edge near `points` (within `tolerance`, degree `min_degree` to `max_degree`;
/// with `smoothing`, the weights of its length, bending and twisting kept small besides).
#[pyfunction]
#[pyo3(signature = (points, tolerance=1e-3, min_degree=1, max_degree=6, smoothing=None))]
pub fn make_spline_approx(points: Vec<Xyz>, tolerance: f64, min_degree: usize, max_degree: usize, smoothing: Option<[f64; 3]>) -> PyResult<Handle> {
    let pts: Vec<Point> = points.into_iter().map(point).collect();
    with_model(|m| {
        let options = fit::ApproximationOptions { min_degree, max_degree, tolerance, smoothing, ..Default::default() };
        let curve = fit::approximate_with(&pts, options)?;
        make::edge_from_curve(m, Curve::BSpline(curve))
    })
    .map(handle)
}

/// A spline face near the grid of `points` (rows along u), within `tolerance`, of degree
/// `min_degree` to `max_degree`.
#[pyfunction]
#[pyo3(signature = (points, tolerance=1e-2, min_degree=1, max_degree=3))]
pub fn make_spline_surface_approx(points: Vec<Vec<Xyz>>, tolerance: f64, min_degree: usize, max_degree: usize) -> PyResult<Handle> {
    let grid: Vec<Vec<Point>> = points.into_iter().map(|row| row.into_iter().map(point).collect()).collect();
    with_model(|m| {
        let options = fit::ApproximationOptions { min_degree, max_degree, tolerance, ..Default::default() };
        let surface = fit::approximate_surface(&grid, options)?;
        make::face_on_surface(m, Surface::BSpline(surface))
    })
    .map(handle)
}

/// The face spanning the closed loop of `edges` (edges or wires, in any order), through
/// `points`, of least bending energy of `order` (OCCT's filling `Degree`).
#[pyfunction]
#[pyo3(signature = (edges, points=Vec::new(), order=3, points_per_edge=15))]
pub fn make_filling(edges: Vec<Handle>, points: Vec<Xyz>, order: usize, points_per_edge: usize) -> PyResult<Handle> {
    let edges: Vec<Shape> = edges.iter().map(|h| h.shape).collect();
    let options = kentron::fill::FillingOptions { order, points: points.into_iter().map(point).collect(), points_per_edge };
    with_model(|m| kentron::fill::n_sided(m, &edges, &options)).map(handle)
}

/// `shape` with `faces` drafted by `angle` (radians) against the pull `direction`, each
/// turned about where it meets the neutral plane through `origin` with `normal` (OCCT's
/// draft angle; the faces tangent to them drafted alike).
#[pyfunction]
pub fn draft_faces(shape: Handle, faces: Vec<Handle>, direction_: Xyz, angle: f64, origin: Xyz, normal: Xyz) -> PyResult<Handle> {
    let pull = direction(direction_)?;
    let neutral = Plane::new(point(origin), direction(normal)?);
    let faces: Vec<kentron::draft::DraftedFace> = faces.iter().map(|h| kentron::draft::DraftedFace { face: h.shape, angle, neutral }).collect();
    with_model(|m| kentron::draft::draft(m, shape.shape, pull, &faces)).map(handle)
}

/// A prism feature (OCCT's BRepFeat_MakePrism, or MakeDPrism with an `angle`): `profile`
/// swept along `direction` (or, with a taper `angle` in radians, its normal), fused with
/// `base` (`additive`) or cut from it; as far as `length`, up to the face `until` (from the
/// face `start`), or, with neither, through all of `base`.
#[pyfunction]
#[pyo3(signature = (base, profile, direction=None, angle=0.0, additive=true, length=None, until=None, start=None))]
#[allow(clippy::too_many_arguments)]
pub fn feature_prism(base: Handle, profile: Handle, direction: Option<Xyz>, angle: f64, additive: bool, length: Option<f64>, until: Option<Handle>, start: Option<Handle>) -> PyResult<Handle> {
    use kentron::feature::{FeatureExtent, FeatureMode, FeatureSweep};
    let sweep = if angle != 0.0 {
        FeatureSweep::DraftedPrism(angle)
    } else {
        let d = direction.ok_or_else(|| pyo3::exceptions::PyValueError::new_err("a prism feature without a direction"))?;
        FeatureSweep::Prism(super::direction(d)?)
    };
    let extent = match (length, until, start) {
        (Some(l), None, None) => FeatureExtent::Length(l),
        (None, Some(u), None) => FeatureExtent::Until(u.shape),
        (None, Some(u), Some(s)) => FeatureExtent::FromUntil(s.shape, u.shape),
        (None, None, None) => FeatureExtent::ThruAll,
        _ => return Err(pyo3::exceptions::PyValueError::new_err("a prism feature's extent: a length, a face up to, or faces from and up to")),
    };
    let mode = if additive { FeatureMode::Add } else { FeatureMode::Remove };
    with_model(|m| kentron::feature::feature(m, base.shape, profile.shape, sweep, mode, extent)).map(handle)
}

/// `face` with the closed wires `holes` (on its surface, within it) cut out.
#[pyfunction]
pub fn face_with_holes(face: Handle, holes: Vec<Handle>) -> PyResult<Handle> {
    let holes: Vec<Shape> = holes.iter().map(|h| h.shape).collect();
    with_model(|m| make::face_with_holes(m, face.shape, &holes)).map(handle)
}

/// Where the infinite line through `origin` along `direction` crosses the faces of `shape`
/// (within `tolerance`): each point, the face, and the line's parameter there (its signed
/// distance from `origin`), in the line's order.
#[pyfunction]
#[pyo3(signature = (shape, origin, direction, tolerance=1e-4))]
pub fn line_intersections(shape: Handle, origin: Xyz, direction: Xyz, tolerance: f64) -> PyResult<Vec<(Xyz, Handle, f64)>> {
    use kentron::geometry::curve::Line;
    let d = kentron::math::Direction::try_from_vector(Vector::new(direction.0, direction.1, direction.2))
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("a line along a zero direction"))?;
    let line = Curve::Line(Line::new(point(origin), d));
    let found = with_model(|m| {
        let mut out = Vec::new();
        for face in m.unique_subshapes(shape.shape, NodeKind::Face) {
            let Some(surface) = m.surface(face).map(|s| kentron::math::Transformable::transformed(s, &face.placement)) else { continue };
            for hit in kentron::geometry::intersect::curve_surface(&line, &surface, tolerance).points {
                // (Within the face: no further from it than the tolerance.)
                let v = make::vertex(m, hit.point);
                if kentron::distance::distance(m, v, face)?.value <= tolerance {
                    out.push((hit.point, face, hit.t));
                }
            }
        }
        out.sort_by(|a, b| a.2.total_cmp(&b.2));
        Ok::<_, kentron::Error>(out)
    })?;
    Ok(found.into_iter().map(|(p, f, t)| ((p.x, p.y, p.z), handle(f), t)).collect())
}

/// The Bézier edge with poles `points`.
#[pyfunction]
pub fn make_bezier(points: Vec<Xyz>) -> PyResult<Handle> {
    let pts: Vec<Point> = points.into_iter().map(point).collect();
    with_model(|m| {
        let curve = BezierCurve::try_new(pts).ok_or_else(|| kentron::Error::InvalidInput("a Bézier curve needs two to 26 poles".into()))?;
        make::edge_from_curve(m, Curve::Bezier(curve))
    })
    .map(handle)
}

/// An elliptic edge (an arc from `angle1` to `angle2` degrees from the frame's X, which
/// is along the `major` radius).
#[pyfunction]
#[pyo3(signature = (major, minor, centre, normal, x=None, angle1=0.0, angle2=360.0))]
pub fn make_ellipse(major: f64, minor: f64, centre: Xyz, normal: Xyz, x: Option<Xyz>, angle1: f64, angle2: f64) -> PyResult<Handle> {
    let f = frame(centre, normal, x)?;
    with_model(|m| {
        // (Of equal radii: the circle, parametrized alike, so that what is made from it stays
        // analytic (a cylinder extruded, a circle offset).)
        let curve = if (major - minor).abs() <= 1e-12 * major.abs().max(minor.abs()) {
            Curve::Circle(kentron::geometry::curve::Circle::new(f, major))
        } else {
            Curve::Ellipse(Ellipse::new(f, major, minor))
        };
        if (angle2 - angle1).abs() >= 360.0 - 1e-9 {
            make::edge_from_curve(m, curve)
        } else {
            make::edge_from_curve_range(m, curve, angle1.to_radians(), angle2.to_radians())
        }
    })
    .map(handle)
}

// --- solids and sweeps -----------------------------------------------------------------

/// A wedge: the box `dx` × `dy` × `dz` whose top (at y = `dy`) spans `[xmin, xmax]` ×
/// `[zmin, zmax]`, in the frame at `origin` with Z along `z` (and X along `x`).
#[pyfunction]
#[pyo3(signature = (dx, dy, dz, xmin, zmin, xmax, zmax, origin, z, x=None))]
#[allow(clippy::too_many_arguments)]
pub fn make_wedge(dx: f64, dy: f64, dz: f64, xmin: f64, zmin: f64, xmax: f64, zmax: f64, origin: Xyz, z: Xyz, x: Option<Xyz>) -> PyResult<Handle> {
    let f = frame(origin, z, x)?;
    with_model(|m| primitives::Wedge::with_top(dx, dy, dz, xmin, zmin, xmax, zmax).with_frame(f).build(m)).map(handle)
}

/// `profile` swept along `spine` (an edge or wire), moved onto the spine first if
/// `contact`, and squared to it if `correction`.
/// At the spine's corners, the profile carried on unturned ("transformed"), or turned with
/// the spine and the sweeps mitred ("right") or rounded ("round"), as CadQuery's transitions;
/// with `frenet`, carried in the spine's Frenet frame.
#[pyfunction]
#[pyo3(signature = (profile, spine, contact=false, correction=false, corner="transformed", frenet=false))]
pub fn pipe(profile: Handle, spine: Handle, contact: bool, correction: bool, corner: &str, frenet: bool) -> PyResult<Handle> {
    let corner = match corner {
        "transformed" => sweep::Corner::Unturned,
        "right" => sweep::Corner::Mitred,
        "round" => sweep::Corner::Rounded,
        other => return Err(pyo3::exceptions::PyValueError::new_err(format!("unknown transition {other}"))),
    };
    with_model(|m| {
        let options = sweep::PipeOptions { contact, correction, frenet };
        sweep::pipe_cornered(m, profile.shape, spine.shape, &options, corner)
    })
    .map(handle)
}

/// `sections` (wires or vertices) swept along `spine`, blended between: a shell, or with
/// `solid` a solid.
#[pyfunction]
#[pyo3(signature = (sections, spine, solid=true, contact=false, correction=false))]
pub fn pipe_through(sections: Vec<Handle>, spine: Handle, solid: bool, contact: bool, correction: bool) -> PyResult<Handle> {
    let options = sweep::PipeOptions { contact, correction, ..Default::default() };
    let sections: Vec<(Shape, sweep::PipeOptions)> = sections.iter().map(|h| (h.shape, options)).collect();
    with_model(|m| sweep::pipe_through(m, &sections, spine.shape, solid)).map(handle)
}

/// The matrix of inertia (rows) about the centre of mass, of the volume, area or length
/// of `shape` (`dimension` 3, 2 or 1).
#[pyfunction]
pub fn inertia(shape: Handle, dimension: u32) -> PyResult<[[f64; 3]; 3]> {
    with_model(|m| {
        let p = match dimension {
            3 => properties::volume(m, shape.shape),
            2 => properties::area(m, shape.shape),
            _ => properties::length(m, shape.shape),
        };
        Ok(std::array::from_fn(|i| std::array::from_fn(|j| p.inertia[(i, j)])))
    })
}

/// `wire` (in one plane) with its corners at `vertices` rounded by `radius`, its edges lines
/// or curves.
#[pyfunction]
pub fn fillet_wire_2d(wire: Handle, radius: f64, vertices: Vec<Handle>) -> PyResult<Handle> {
    let vertices: Vec<Shape> = vertices.iter().map(|h| h.shape).collect();
    with_model(|m| kentron::blend::fillet_wire_2d(m, wire.shape, radius, Some(&vertices))).map(handle)
}

/// `shape` with `edges` chamfered: each `(edge, d1, d2, face)` cut back by `d1` along `face`
/// (either face if `None`) and `d2` along the other.
#[pyfunction]
pub fn chamfer(shape: Handle, edges: Vec<(Handle, f64, f64, Option<Handle>)>) -> PyResult<Handle> {
    use kentron::blend::{ChamferSize, ChamferedEdge};
    let edges: Vec<ChamferedEdge> = edges
        .into_iter()
        .map(|(e, d1, d2, f)| ChamferedEdge { edge: e.shape, size: ChamferSize::Distances(d1, d2), face: f.map(|f| f.shape) })
        .collect();
    with_model(|m| kentron::blend::chamfer_edges(m, shape.shape, &edges)).map(handle)
}

// --- faces -----------------------------------------------------------------------------

/// The face's parameter bounds (u0, u1, v0, v1).
#[pyfunction]
pub fn face_uv_bounds(face: Handle) -> PyResult<(f64, f64, f64, f64)> {
    with_model(|m| Ok(m.face_uv_bounds(face.shape)))
}

/// (point, unit normal as used, u derivative, v derivative) at each of the parameters
/// `uvs`.
#[pyfunction]
pub fn face_points(face: Handle, uvs: Vec<(f64, f64)>) -> PyResult<Vec<(Xyz, Xyz, Xyz, Xyz)>> {
    with_model(|m| {
        let s = m.surface(face.shape).ok_or_else(|| kentron::Error::InvalidInput("a face without a surface".into()))?.clone();
        let s = s.transformed_by(&face.shape.placement);
        let flip = face.shape.orientation == Orientation::Reversed;
        let v = |v: Vector| (v.x, v.y, v.z);
        Ok(uvs
            .into_iter()
            .map(|(a, b)| {
                let (p, du, dv) = s.d1(a, b);
                let n = du.cross(dv);
                let n = if n.length() > 0.0 { n * (1.0 / n.length()) } else { n };
                let n = if flip { -n } else { n };
                (xyz(p), v(n), v(du), v(dv))
            })
            .collect())
    })
}

/// The parameters of the face's surface points nearest `points`.
#[pyfunction]
pub fn face_parameters(face: Handle, points: Vec<Xyz>) -> PyResult<Vec<(f64, f64)>> {
    with_model(|m| {
        let s = m.surface(face.shape).ok_or_else(|| kentron::Error::InvalidInput("a face without a surface".into()))?.clone();
        let s = s.transformed_by(&face.shape.placement);
        Ok(points.into_iter().map(|p| nearest_parameters(&s, point(p))).collect())
    })
}

// --- topology ------------------------------------------------------------------------

/// The edges of a wire in the order they connect, each oriented along it.
#[pyfunction]
pub fn ordered_edges(wire: Handle) -> PyResult<Vec<Handle>> {
    with_model(|m| Ok(m.ordered_edges(wire.shape).into_iter().map(handle).collect()))
}

/// A new edge on `edge`'s curve over `[u0, u1]`.
#[pyfunction]
pub fn edge_trimmed(edge: Handle, u0: f64, u1: f64) -> PyResult<Handle> {
    with_model(|m| {
        let (curve, _, _) = m.edge_curve(edge.shape).ok_or_else(|| kentron::Error::InvalidInput("a degenerate edge".into()))?;
        make::edge_from_curve_range(m, curve, u0, u1)
    })
    .map(handle)
}

/// A new face on `face`'s surface over `[u0, u1] × [v0, v1]`.
#[pyfunction]
pub fn face_trimmed(face: Handle, u0: f64, u1: f64, v0: f64, v1: f64) -> PyResult<Handle> {
    with_model(|m| {
        let s = m.surface(face.shape).ok_or_else(|| kentron::Error::InvalidInput("a face without a surface".into()))?.clone();
        let s = s.transformed_by(&face.shape.placement);
        make::face_on_domain(m, s, u0, u1, v0, v1)
    })
    .map(handle)
}

/// The plane face over `[u0, u1] × [v0, v1]` of the plane at `origin` with normal `normal`
/// (and X along `x`).
#[pyfunction]
#[pyo3(signature = (origin, normal, u0, u1, v0, v1, x=None))]
pub fn make_plane_face(origin: Xyz, normal: Xyz, u0: f64, u1: f64, v0: f64, v1: f64, x: Option<Xyz>) -> PyResult<Handle> {
    let f = frame(origin, normal, x)?;
    with_model(|m| make::face_on_domain(m, Surface::Plane(Plane::from_frame(f)), u0, u1, v0, v1)).map(handle)
}

/// `shapes`' faces sewn together where their boundaries meet within `tolerance`.
#[pyfunction]
#[pyo3(signature = (shapes, tolerance=1e-6))]
pub fn sew(shapes: Vec<Handle>, tolerance: f64) -> PyResult<Handle> {
    let shapes: Vec<Shape> = shapes.iter().map(|h| h.shape).collect();
    with_model(|m| kentron::sew::sew(m, &shapes, &kentron::sew::SewOptions { tolerance, ..Default::default() })).map(handle)
}

/// A shell of `faces` as they are (not sewn).
#[pyfunction]
pub fn make_shell(faces: Vec<Handle>) -> PyResult<Handle> {
    with_model(|m| {
        let shell = m.new_empty(NodeKind::Shell);
        for f in faces {
            let face = unplaced(m, f.shape)?;
            m.add(shell, face)?;
        }
        Ok(shell)
    })
    .map(handle)
}

/// `shape` copied where its placement puts it, where it has one (a shell or solid holds its
/// children unplaced).
fn unplaced(m: &mut Model, shape: Shape) -> kentron::Result<Shape> {
    if shape.placement == kentron::math::Transform::IDENTITY {
        return Ok(shape);
    }
    kentron::modify::copy(m, shape, false)
}

/// The solid bounded by `shells` (the first the outer one), turned so that it has a
/// positive volume.
#[pyfunction]
pub fn make_solid(shells: Vec<Handle>) -> PyResult<Handle> {
    with_model(|m| {
        // (Placed shells, as those of a solid on a workplane, moved where they are.)
        let shells: Vec<Shape> = shells.iter().map(|s| unplaced(m, s.shape)).collect::<kentron::Result<_>>()?;
        let build = |m: &mut Model, flip: bool| -> kentron::Result<Shape> {
            let solid = m.new_empty(NodeKind::Solid);
            for (i, &s) in shells.iter().enumerate() {
                m.add(solid, if flip && i == 0 { s.reversed() } else { s })?;
            }
            Ok(solid)
        };
        let solid = build(m, false)?;
        if properties::volume(m, solid).mass < 0.0 { build(m, true) } else { Ok(solid) }
    })
    .map(handle)
}

// --- cells -----------------------------------------------------------------------------

/// `shapes` split by one another (OCCT's general fuse): each cell, with whether it lies
/// inside each of the shapes.
#[pyfunction]
pub fn cells(shapes: Vec<Handle>) -> PyResult<Vec<(Vec<bool>, Handle)>> {
    let shapes: Vec<Shape> = shapes.iter().map(|h| h.shape).collect();
    with_model(|m| Ok(boolean::cells(m, &shapes)?.into_iter().map(|c| (c.inside, handle(c.shape))).collect()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cells, m)?)?;
    m.add_function(wrap_pyfunction!(fillet_wire_2d, m)?)?;
    m.add_function(wrap_pyfunction!(ordered_edges, m)?)?;
    m.add_function(wrap_pyfunction!(edge_trimmed, m)?)?;
    m.add_function(wrap_pyfunction!(face_trimmed, m)?)?;
    m.add_function(wrap_pyfunction!(make_plane_face, m)?)?;
    m.add_function(wrap_pyfunction!(sew, m)?)?;
    m.add_function(wrap_pyfunction!(make_solid, m)?)?;
    m.add_function(wrap_pyfunction!(make_shell, m)?)?;
    m.add_function(wrap_pyfunction!(edge_conic, m)?)?;
    m.add_function(wrap_pyfunction!(edge_bspline, m)?)?;
    m.add_function(wrap_pyfunction!(make_bspline, m)?)?;
    let _ = direction;
    m.add_function(wrap_pyfunction!(curve_info, m)?)?;
    m.add_function(wrap_pyfunction!(curve_derivatives, m)?)?;
    m.add_function(wrap_pyfunction!(curve_length, m)?)?;
    m.add_function(wrap_pyfunction!(curve_parameters_at_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(curve_parameters_of_points, m)?)?;
    m.add_function(wrap_pyfunction!(curve_plane_normal, m)?)?;
    m.add_function(wrap_pyfunction!(make_spline, m)?)?;
    m.add_function(wrap_pyfunction!(make_spline_approx, m)?)?;
    m.add_function(wrap_pyfunction!(make_spline_surface_approx, m)?)?;
    m.add_function(wrap_pyfunction!(make_filling, m)?)?;
    m.add_function(wrap_pyfunction!(draft_faces, m)?)?;
    m.add_function(wrap_pyfunction!(feature_prism, m)?)?;
    m.add_function(wrap_pyfunction!(face_with_holes, m)?)?;
    m.add_function(wrap_pyfunction!(line_intersections, m)?)?;
    m.add_function(wrap_pyfunction!(make_bezier, m)?)?;
    m.add_function(wrap_pyfunction!(make_ellipse, m)?)?;
    m.add_function(wrap_pyfunction!(make_wedge, m)?)?;
    m.add_function(wrap_pyfunction!(pipe, m)?)?;
    m.add_function(wrap_pyfunction!(pipe_through, m)?)?;
    m.add_function(wrap_pyfunction!(inertia, m)?)?;
    m.add_function(wrap_pyfunction!(chamfer, m)?)?;
    m.add_function(wrap_pyfunction!(face_uv_bounds, m)?)?;
    m.add_function(wrap_pyfunction!(face_points, m)?)?;
    m.add_function(wrap_pyfunction!(face_parameters, m)?)?;
    Ok(())
}
