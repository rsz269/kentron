//! OCCT's native ASCII BRep format (`.brep` files, as used by DRAW's `restore`/`save`).
//!
//! A file holds tables, each element referred to by its 1-based index: locations
//! (transformations), 2D curves, 3D curves, 3D polygons, polygons on triangulations,
//! surfaces, triangulations, then the shapes ("TShapes"), sub-shapes before the shapes using
//! them, and last a reference to the root shape. A reference to a shape is its orientation
//! sign, its number counted back from the end of the table, and a location.
//!
//! Reading builds the shapes in a [`Model`]. Locations of sub-shapes are applied to their
//! geometry, except those of a compound's children, which become the children's placements.
//! Meshes (polygons and triangulations) are skipped.

use std::collections::HashMap;
use std::io::{Read, Write};

use crate::error::{Error, Result};
use crate::geometry::curve::{
    BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, OffsetCurve, Parabola, TrimmedCurve,
};
use crate::geometry::curve2::{Circle2, Curve2, Ellipse2, Hyperbola2, Line2, OffsetCurve2, Parabola2};
use crate::geometry::surface::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ExtrusionSurface, KnotVector, OffsetSurface, Plane,
    RevolutionSurface, Sphere, Surface, Torus, TrimmedSurface,
};
use crate::geometry::{ControlPoint, ParametricCurve};
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Transformable, Vector};
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape, ShellId, VertexId, WireId};

/// Revisions of the BRep format.
#[doc(alias = "TopTools_FormatVersion")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FormatVersion {
    /// Without the (u, v) end points of curves on surfaces.
    V1,
    /// Adds the (u, v) end points of curves on surfaces.
    V2,
    /// Adds per-node normals of triangulation-only faces.
    #[default]
    V3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WriteOptions {
    pub version: FormatVersion,
    /// Also write face meshes.
    pub with_meshes: bool,
    /// Also write mesh normals.
    pub with_normals: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self { version: FormatVersion::default(), with_meshes: true, with_normals: false }
    }
}

/// Writes `shape` and everything below it.
#[doc(alias = "BRepTools::Write")]
pub fn write(model: &Model, shape: Shape, out: &mut dyn Write, options: &WriteOptions) -> Result<()> {
    let text = Writer::new(model, options.version).file(shape)?;
    out.write_all(text.as_bytes()).map_err(|e| Error::InvalidInput(format!("writing BRep: {e}")))
}

/// Reads a shape into `model`; any format version is accepted.
#[doc(alias = "BRepTools::Read")]
pub fn read(model: &mut Model, input: &mut dyn Read) -> Result<Shape> {
    let mut text = String::new();
    input.read_to_string(&mut text).map_err(|e| Error::InvalidInput(format!("reading BRep: {e}")))?;
    let file = parse(&text)?;
    Builder::new(model, &file).root()
}

/// A point, curve or surface saved by OCCT's DRAW.
#[derive(Debug, Clone, PartialEq)]
pub enum SavedGeometry {
    Point(Point),
    Point2(Point2),
    Curve(Curve),
    Curve2(Curve2),
    Surface(Surface),
}

/// Reads a point, curve or surface saved by OCCT's DRAW (`DrawTrSurf_Point`,
/// `DrawTrSurf_Curve`, `DrawTrSurf_Curve2d`, `DrawTrSurf_Surface` and their B-spline and
/// Bézier forms: the record names the kind, then holds the geometry as in a BRep file).
/// `None` if the input holds something else (a shape, a triangulation).
#[doc(alias = "DrawTrSurf")]
pub fn read_saved_geometry(input: &mut dyn Read) -> Result<Option<SavedGeometry>> {
    let mut text = String::new();
    input.read_to_string(&mut text).map_err(|e| Error::InvalidInput(format!("reading DRAW geometry: {e}")))?;
    let mut w = Words::new(&text);
    let Ok(kind) = w.word() else { return Ok(None) };
    Ok(Some(match kind {
        "DrawTrSurf_Point" => match w.int()? {
            1 => SavedGeometry::Point(w.point()?),
            _ => SavedGeometry::Point2(w.point2()?),
        },
        "DrawTrSurf_Curve" | "DrawTrSurf_BSplineCurve" | "DrawTrSurf_BezierCurve" => {
            SavedGeometry::Curve(RawCurve::read(&mut w)?.build()?)
        }
        "DrawTrSurf_Curve2d" | "DrawTrSurf_BSplineCurve2d" | "DrawTrSurf_BezierCurve2d" => {
            SavedGeometry::Curve2(RawCurve2::read(&mut w)?.build()?)
        }
        "DrawTrSurf_Surface" | "DrawTrSurf_BSplineSurface" | "DrawTrSurf_BezierSurface" => {
            SavedGeometry::Surface(RawSurface::read(&mut w)?.build()?)
        }
        _ => return Ok(None),
    }))
}

// --- parsing ------------------------------------------------------------------------------

const VERSIONS: [(&str, FormatVersion); 3] = [
    ("CASCADE Topology V1, (c) Matra-Datavision", FormatVersion::V1),
    ("CASCADE Topology V2, (c) Matra-Datavision", FormatVersion::V2),
    ("CASCADE Topology V3, (c) Open Cascade", FormatVersion::V3),
];

fn bad(msg: impl Into<String>) -> Error {
    Error::InvalidInput(format!("BRep: {}", msg.into()))
}

/// Whitespace-separated words.
struct Words<'a> {
    words: std::str::SplitWhitespace<'a>,
    /// The rest of a word an integer was read from (C++ streams stop at the first character
    /// that is not a digit: old files write `2CN` for `2` and `CN`).
    rest: Option<&'a str>,
}

impl<'a> Words<'a> {
    fn new(text: &'a str) -> Self {
        Self { words: text.split_whitespace(), rest: None }
    }

    fn word(&mut self) -> Result<&'a str> {
        if let Some(rest) = self.rest.take() {
            return Ok(rest);
        }
        self.words.next().ok_or_else(|| bad("unexpected end of file"))
    }

    fn expect(&mut self, keyword: &str) -> Result<()> {
        let w = self.word()?;
        if w == keyword {
            Ok(())
        } else {
            Err(bad(format!("expected {keyword}, found {w}")))
        }
    }

    /// A real number; like OCCT's `strtod`, the longest numeric prefix (0 if none).
    fn real(&mut self) -> Result<f64> {
        let w = self.word()?;
        if let Ok(v) = w.parse() {
            return Ok(v);
        }
        Ok((1..w.len()).rev().find_map(|n| w.get(..n)?.parse().ok()).unwrap_or(0.0))
    }

    fn int(&mut self) -> Result<i64> {
        let w = self.word()?;
        let digits = w.char_indices().find(|&(i, c)| !(c.is_ascii_digit() || (i == 0 && (c == '-' || c == '+'))));
        let (number, rest) = match digits {
            Some((i, _)) => (&w[..i], Some(&w[i..])),
            None => (w, None),
        };
        let value = number.parse().map_err(|_| bad(format!("expected an integer, found {w}")))?;
        self.rest = rest;
        Ok(value)
    }

    fn count(&mut self) -> Result<usize> {
        usize::try_from(self.int()?).map_err(|_| bad("negative count"))
    }

    fn flag(&mut self) -> Result<bool> {
        Ok(self.int()? != 0)
    }

    fn point(&mut self) -> Result<Point> {
        Ok(Point::new(self.real()?, self.real()?, self.real()?))
    }

    fn point2(&mut self) -> Result<Point2> {
        Ok(Point2::new(self.real()?, self.real()?))
    }

    fn direction(&mut self) -> Result<Direction> {
        let (x, y, z) = (self.real()?, self.real()?, self.real()?);
        Direction::try_new(x, y, z).ok_or_else(|| bad("a null direction"))
    }

    fn direction2(&mut self) -> Result<Direction2> {
        let (x, y) = (self.real()?, self.real()?);
        Direction2::try_new(x, y).ok_or_else(|| bad("a null direction"))
    }

    /// `gp_Ax2`: origin, main direction, X direction, Y direction (right-handed; the Y
    /// direction is not used).
    fn frame(&mut self) -> Result<Frame> {
        let (origin, z, x) = (self.point()?, self.direction()?, self.direction()?);
        self.direction()?;
        Ok(Frame::new(origin, z, x))
    }

    /// `gp_Ax3`: as `gp_Ax2`, left-handed if the Y direction says so.
    fn frame3(&mut self) -> Result<Frame> {
        let (origin, z, x, y) = (self.point()?, self.direction()?, self.direction()?, self.direction()?);
        let mut frame = Frame::new(origin, z, x);
        if Vector::from(y).dot(Vector::from(z).cross(Vector::from(x))) < 0.0 {
            frame.flip_y();
        }
        Ok(frame)
    }

    /// `gp_Ax22d`: origin, X and Y directions.
    fn frame2(&mut self) -> Result<Frame2> {
        let (origin, x, y) = (self.point2()?, self.direction2()?, self.direction2()?);
        Ok(Frame2::new(origin, x, x.x() * y.y() - x.y() * y.x() > 0.0))
    }

    fn knots(&mut self, count: usize) -> Result<(Vec<f64>, Vec<usize>)> {
        let mut knots = Vec::with_capacity(count);
        let mut multiplicities = Vec::with_capacity(count);
        for _ in 0..count {
            knots.push(self.real()?);
            multiplicities.push(self.count()?);
        }
        Ok((knots, multiplicities))
    }
}

fn pole3(w: &mut Words) -> Result<Point> {
    w.point()
}

fn pole2(w: &mut Words) -> Result<Point2> {
    w.point2()
}

/// A B-spline curve's data, built when used: periodic B-splines are not supported yet.
#[derive(Debug, Clone)]
struct Spline<P> {
    poles: Vec<P>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    multiplicities: Vec<usize>,
    degree: usize,
    periodic: bool,
}

impl<P: ControlPoint> Spline<P> {
    fn read(w: &mut Words, pole: impl Fn(&mut Words) -> Result<P>) -> Result<Self> {
        let (rational, periodic) = (w.flag()?, w.flag()?);
        let (degree, pole_count, knot_count) = (w.count()?, w.count()?, w.count()?);
        let mut poles = Vec::with_capacity(pole_count);
        let mut weights = Vec::with_capacity(pole_count);
        for _ in 0..pole_count {
            poles.push(pole(w)?);
            if rational {
                weights.push(w.real()?);
            }
        }
        let (knots, multiplicities) = w.knots(knot_count)?;
        Ok(Self { poles, weights: rational.then_some(weights), knots, multiplicities, degree, periodic })
    }

    fn build(&self) -> Result<BSplineCurve<P>> {
        let (poles, knots, mults) = (self.poles.clone(), self.knots.clone(), self.multiplicities.clone());
        let curve = if self.periodic {
            BSplineCurve::try_new_periodic(poles, knots, mults, self.degree)?
        } else {
            BSplineCurve::try_new(poles, knots, mults, self.degree)?
        };
        Ok(match &self.weights {
            Some(w) => curve.with_weights(w.clone()),
            None => curve,
        })
    }
}

fn bezier<P: ControlPoint>(w: &mut Words, pole: impl Fn(&mut Words) -> Result<P>) -> Result<BezierCurve<P>> {
    let rational = w.flag()?;
    let degree = w.count()?;
    let mut poles = Vec::with_capacity(degree + 1);
    let mut weights = Vec::with_capacity(degree + 1);
    for _ in 0..=degree {
        poles.push(pole(w)?);
        if rational {
            weights.push(w.real()?);
        }
    }
    Ok(if rational { BezierCurve::rational(poles, weights) } else { BezierCurve::new(poles) })
}

/// A 3D curve as read.
#[derive(Debug, Clone)]
enum RawCurve {
    Ready(Curve),
    BSpline(Spline<Point>),
    Trimmed(Box<RawCurve>, f64, f64),
    Offset(Box<RawCurve>, f64, Direction),
}

impl RawCurve {
    fn read(w: &mut Words) -> Result<Self> {
        Ok(match w.int()? {
            1 => RawCurve::Ready(Curve::Line(Line::new(w.point()?, w.direction()?))),
            2 => RawCurve::Ready(Curve::Circle(Circle::new(w.frame()?, w.real()?))),
            3 => RawCurve::Ready(Curve::Ellipse(Ellipse::new(w.frame()?, w.real()?, w.real()?))),
            4 => RawCurve::Ready(Curve::Parabola(Parabola::new(w.frame()?, w.real()?))),
            5 => RawCurve::Ready(Curve::Hyperbola(Hyperbola::new(w.frame()?, w.real()?, w.real()?))),
            6 => RawCurve::Ready(Curve::Bezier(bezier(w, pole3)?)),
            7 => RawCurve::BSpline(Spline::read(w, pole3)?),
            8 => {
                let (u1, u2) = (w.real()?, w.real()?);
                RawCurve::Trimmed(Box::new(RawCurve::read(w)?), u1, u2)
            }
            9 => {
                let (offset, direction) = (w.real()?, w.direction()?);
                RawCurve::Offset(Box::new(RawCurve::read(w)?), offset, direction)
            }
            t => return Err(bad(format!("unknown curve type {t}"))),
        })
    }

    fn build(&self) -> Result<Curve> {
        Ok(match self {
            RawCurve::Ready(c) => c.clone(),
            RawCurve::BSpline(s) => Curve::BSpline(s.build()?),
            RawCurve::Trimmed(basis, u1, u2) => Curve::Trimmed(trimmed(basis.build()?, *u1, *u2)),
            RawCurve::Offset(basis, offset, direction) => {
                Curve::Offset(OffsetCurve::new(basis.build()?, *offset, *direction))
            }
        })
    }
}

/// `Geom_TrimmedCurve`: a reversed range on a non-periodic curve trims the reversed curve.
fn trimmed<C: ParametricCurve>(mut basis: C, u1: f64, u2: f64) -> TrimmedCurve<C> {
    if u1 > u2 && !basis.is_periodic() {
        let (r1, r2) = (basis.reversed_parameter(u1), basis.reversed_parameter(u2));
        basis.reverse();
        return TrimmedCurve::new(basis, r1, r2);
    }
    TrimmedCurve::new(basis, u1, u2)
}

/// A 2D curve as read.
#[derive(Debug, Clone)]
enum RawCurve2 {
    Ready(Curve2),
    BSpline(Spline<Point2>),
    Trimmed(Box<RawCurve2>, f64, f64),
    Offset(Box<RawCurve2>, f64),
}

impl RawCurve2 {
    fn read(w: &mut Words) -> Result<Self> {
        Ok(match w.int()? {
            1 => RawCurve2::Ready(Curve2::Line(Line2::new(w.point2()?, w.direction2()?))),
            2 => RawCurve2::Ready(Curve2::Circle(Circle2::new(w.frame2()?, w.real()?))),
            3 => RawCurve2::Ready(Curve2::Ellipse(Ellipse2::new(w.frame2()?, w.real()?, w.real()?))),
            4 => RawCurve2::Ready(Curve2::Parabola(Parabola2::new(w.frame2()?, w.real()?))),
            5 => RawCurve2::Ready(Curve2::Hyperbola(Hyperbola2::new(w.frame2()?, w.real()?, w.real()?))),
            6 => RawCurve2::Ready(Curve2::Bezier(bezier(w, pole2)?)),
            7 => RawCurve2::BSpline(Spline::read(w, pole2)?),
            8 => {
                let (u1, u2) = (w.real()?, w.real()?);
                RawCurve2::Trimmed(Box::new(RawCurve2::read(w)?), u1, u2)
            }
            9 => {
                let offset = w.real()?;
                RawCurve2::Offset(Box::new(RawCurve2::read(w)?), offset)
            }
            t => return Err(bad(format!("unknown 2D curve type {t}"))),
        })
    }

    fn build(&self) -> Result<Curve2> {
        Ok(match self {
            RawCurve2::Ready(c) => c.clone(),
            RawCurve2::BSpline(s) => Curve2::BSpline(s.build()?),
            RawCurve2::Trimmed(basis, u1, u2) => Curve2::Trimmed(trimmed(basis.build()?, *u1, *u2)),
            RawCurve2::Offset(basis, offset) => Curve2::Offset(OffsetCurve2::new(basis.build()?, *offset)),
        })
    }
}

/// A surface as read.
#[derive(Debug, Clone)]
enum RawSurface {
    Ready(Surface),
    Extrusion(RawCurve, Direction),
    Revolution(RawCurve, Axis),
    BSpline { poles: Vec<Vec<Point>>, weights: Option<Vec<Vec<f64>>>, u: KnotVector, v: KnotVector },
    Trimmed(Box<RawSurface>, [f64; 4]),
    Offset(Box<RawSurface>, f64),
}

impl RawSurface {
    fn read(w: &mut Words) -> Result<Self> {
        Ok(match w.int()? {
            1 => RawSurface::Ready(Surface::Plane(Plane::from_frame(w.frame3()?))),
            2 => RawSurface::Ready(Surface::Cylinder(Cylinder::new(w.frame3()?, w.real()?))),
            3 => {
                let (frame, radius, angle) = (w.frame3()?, w.real()?, w.real()?);
                RawSurface::Ready(Surface::Cone(Cone::new(frame, angle, radius)))
            }
            4 => RawSurface::Ready(Surface::Sphere(Sphere::new(w.frame3()?, w.real()?))),
            5 => RawSurface::Ready(Surface::Torus(Torus::new(w.frame3()?, w.real()?, w.real()?))),
            6 => {
                let direction = w.direction()?;
                RawSurface::Extrusion(RawCurve::read(w)?, direction)
            }
            7 => {
                let axis = Axis::new(w.point()?, w.direction()?);
                RawSurface::Revolution(RawCurve::read(w)?, axis)
            }
            8 => {
                let (u_rational, v_rational) = (w.flag()?, w.flag()?);
                let (u_degree, v_degree) = (w.count()?, w.count()?);
                let (poles, weights) = grid(w, u_degree + 1, v_degree + 1, u_rational || v_rational)?;
                RawSurface::Ready(Surface::Bezier(match weights {
                    Some(weights) => BezierSurface::rational(poles, weights),
                    None => BezierSurface::new(poles),
                }))
            }
            9 => {
                let (u_rational, v_rational) = (w.flag()?, w.flag()?);
                let (u_periodic, v_periodic) = (w.flag()?, w.flag()?);
                let (u_degree, v_degree) = (w.count()?, w.count()?);
                let (u_poles, v_poles) = (w.count()?, w.count()?);
                let (u_knots, v_knots) = (w.count()?, w.count()?);
                let (poles, weights) = grid(w, u_poles, v_poles, u_rational || v_rational)?;
                let knot_vector = |w: &mut Words, count, degree, periodic| -> Result<KnotVector> {
                    let (knots, mults) = w.knots(count)?;
                    Ok(if periodic {
                        KnotVector::periodic(knots, mults, degree)
                    } else {
                        KnotVector::new(knots, mults, degree)
                    })
                };
                let u = knot_vector(w, u_knots, u_degree, u_periodic)?;
                let v = knot_vector(w, v_knots, v_degree, v_periodic)?;
                RawSurface::BSpline { poles, weights, u, v }
            }
            10 => {
                let bounds = [w.real()?, w.real()?, w.real()?, w.real()?];
                RawSurface::Trimmed(Box::new(RawSurface::read(w)?), bounds)
            }
            11 => {
                let offset = w.real()?;
                RawSurface::Offset(Box::new(RawSurface::read(w)?), offset)
            }
            t => return Err(bad(format!("unknown surface type {t}"))),
        })
    }

    fn build(&self) -> Result<Surface> {
        Ok(match self {
            RawSurface::Ready(s) => s.clone(),
            RawSurface::Extrusion(c, d) => Surface::Extrusion(ExtrusionSurface::new(c.build()?, *d)),
            RawSurface::Revolution(c, a) => Surface::Revolution(RevolutionSurface::new(c.build()?, *a)),
            RawSurface::BSpline { poles, weights, u, v } => {
                let surface = BSplineSurface::try_new(poles.clone(), u.clone(), v.clone())?;
                Surface::BSpline(match weights {
                    Some(w) => surface.with_weights(w.clone()),
                    None => surface,
                })
            }
            RawSurface::Trimmed(basis, [u1, u2, v1, v2]) => {
                Surface::Trimmed(TrimmedSurface::new(basis.build()?, *u1, *u2, *v1, *v2))
            }
            RawSurface::Offset(basis, offset) => Surface::Offset(OffsetSurface::new(basis.build()?, *offset)),
        })
    }
}

/// A grid of poles `[i][j]` (with weights if `rational`).
#[allow(clippy::type_complexity)]
fn grid(w: &mut Words, nu: usize, nv: usize, rational: bool) -> Result<(Vec<Vec<Point>>, Option<Vec<Vec<f64>>>)> {
    let mut poles = vec![Vec::with_capacity(nv); nu];
    let mut weights = vec![Vec::with_capacity(nv); nu];
    for i in 0..nu {
        for _ in 0..nv {
            poles[i].push(w.point()?);
            if rational {
                weights[i].push(w.real()?);
            }
        }
    }
    Ok((poles, rational.then_some(weights)))
}

/// A reference to a shape: its index in the table (0-based), orientation and location.
#[derive(Debug, Clone, Copy)]
struct Ref {
    index: usize,
    orientation: Orientation,
    location: usize,
}

/// A curve on a surface: one 2D curve, or two for a seam (forward and reversed uses).
#[derive(Debug, Clone)]
struct OnSurface {
    curve: usize,
    reversed_curve: Option<usize>,
    surface: usize,
    location: usize,
    first: f64,
    last: f64,
}

#[derive(Debug, Clone)]
enum Geometry {
    None,
    Vertex {
        point: Point,
        tolerance: f64,
    },
    Edge {
        tolerance: f64,
        degenerate: bool,
        /// Curve, location, range.
        curve: Option<(usize, usize, f64, f64)>,
        on_surfaces: Vec<OnSurface>,
    },
    Face {
        tolerance: f64,
        /// Surface and location; `None` for a face with only a triangulation.
        surface: Option<(usize, usize)>,
    },
}

#[derive(Debug, Clone)]
struct RawShape {
    kind: NodeKind,
    geometry: Geometry,
    children: Vec<Ref>,
}

struct File {
    /// Locations as transformations; index 0 is the identity.
    locations: Vec<Transform>,
    curves2: Vec<RawCurve2>,
    curves: Vec<RawCurve>,
    surfaces: Vec<RawSurface>,
    shapes: Vec<RawShape>,
    root: Ref,
}

fn parse(text: &str) -> Result<File> {
    // The words after the header line naming the version.
    let (start, version) = VERSIONS
        .iter()
        .filter_map(|&(header, v)| text.find(header).map(|i| (i + header.len(), v)))
        .min_by_key(|&(i, _)| i)
        .ok_or_else(|| bad("no CASCADE Topology header"))?;
    let mut w = Words::new(&text[start..]);
    let locations = read_locations(&mut w)?;
    w.expect("Curve2ds")?;
    let curves2 = (0..w.count()?).map(|_| RawCurve2::read(&mut w)).collect::<Result<Vec<_>>>()?;
    w.expect("Curves")?;
    let curves = (0..w.count()?).map(|_| RawCurve::read(&mut w)).collect::<Result<Vec<_>>>()?;
    skip_meshes(&mut w, version, "Polygon3D")?;
    w.expect("Surfaces")?;
    let surfaces = (0..w.count()?).map(|_| RawSurface::read(&mut w)).collect::<Result<Vec<_>>>()?;
    skip_meshes(&mut w, version, "Triangulations")?;
    w.expect("TShapes")?;
    let count = w.count()?;
    let mut shapes = Vec::with_capacity(count);
    for _ in 0..count {
        shapes.push(read_shape(&mut w, version, count)?);
    }
    let root = read_ref(&mut w, count)?.ok_or_else(|| bad("no root shape"))?;
    Ok(File { locations, curves2, curves, surfaces, shapes, root })
}

/// The 3D curve of `pcurve` on `surface` over `[first, last]`, at the same parameters:
/// exact on a plane, otherwise interpolated through points of it (doubled until within
/// 1e-9, relative, between them).
fn curve_on_surface(pcurve: &Curve2, surface: &Surface, first: f64, last: f64) -> Result<Curve> {
    use crate::geometry::ParametricSurface;
    if let Surface::Plane(plane) = surface {
        return Ok(pcurve.to_3d(plane.frame()));
    }
    let at = |t: f64| {
        let uv = pcurve.point(t);
        surface.point(uv.x, uv.y)
    };
    let mut n = 16;
    loop {
        let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
        let points: Vec<Point> = params.iter().map(|&t| at(t)).collect();
        let curve = crate::geometry::fit::interpolate_with_parameters(&points, &params)
            .map_err(|_| bad("an edge whose curve on a surface cannot be followed"))?;
        let size = points.iter().fold(1.0f64, |m, p| m.max((*p - Point::ORIGIN).length()));
        let error = params.windows(2).map(|w| {
            let t = 0.5 * (w[0] + w[1]);
            curve.point(t).distance(at(t))
        }).fold(0.0, f64::max);
        if error <= 1e-9 * size || n >= 4096 {
            return Ok(Curve::BSpline(curve));
        }
        n *= 2;
    }
}

/// The location table: elementary transformations (`1` and a 3×4 matrix) and products of
/// powers of earlier ones (`2`, then pairs of index and power ending with `0`). Identities
/// are not numbered.
fn read_locations(w: &mut Words) -> Result<Vec<Transform>> {
    w.expect("Locations")?;
    let count = w.count()?;
    let mut locations = vec![Transform::IDENTITY];
    for _ in 0..count {
        let t = match w.int()? {
            1 => {
                let mut rows = [[0.0; 4]; 3];
                for row in &mut rows {
                    for v in row.iter_mut() {
                        *v = w.real()?;
                    }
                }
                Some(Transform::from_rows(rows).ok_or_else(|| bad("a singular location"))?)
            }
            2 => {
                let mut product = Transform::IDENTITY;
                let mut empty = true;
                loop {
                    let index = w.count()?;
                    if index == 0 {
                        break;
                    }
                    let power = w.int()? as i32;
                    let l = locations.get(index).ok_or_else(|| bad(format!("no location {index}")))?;
                    product = l.powi(power) * product;
                    empty &= power == 0;
                }
                (!empty).then_some(product)
            }
            t => return Err(bad(format!("unknown location type {t}"))),
        };
        if let Some(t) = t {
            locations.push(t);
        }
    }
    Ok(locations)
}

/// Skips the mesh tables starting at `first`: 3D polygons and polygons on triangulations,
/// or triangulations.
fn skip_meshes(w: &mut Words, version: FormatVersion, first: &str) -> Result<()> {
    if first == "Polygon3D" {
        w.expect("Polygon3D")?;
        for _ in 0..w.count()? {
            let (nodes, has_parameters) = (w.count()?, w.flag()?);
            w.real()?;
            for _ in 0..nodes * 3 + if has_parameters { nodes } else { 0 } {
                w.word()?;
            }
        }
        w.expect("PolygonOnTriangulations")?;
        for _ in 0..w.count()? {
            let nodes = w.count()?;
            for _ in 0..nodes {
                w.word()?;
            }
            w.word()?; // "p"
            w.real()?;
            if w.flag()? {
                for _ in 0..nodes {
                    w.word()?;
                }
            }
        }
    } else {
        w.expect("Triangulations")?;
        for _ in 0..w.count()? {
            let (nodes, triangles, has_uv) = (w.count()?, w.count()?, w.flag()?);
            let has_normals = version == FormatVersion::V3 && w.flag()?;
            w.real()?;
            let words = nodes * 3
                + if has_uv { nodes * 2 } else { 0 }
                + triangles * 3
                + if has_normals { nodes * 3 } else { 0 };
            for _ in 0..words {
                w.word()?;
            }
        }
    }
    Ok(())
}

/// A shape reference (`+3 0`), or `None` for the end of a list (`*`).
fn read_ref(w: &mut Words, count: usize) -> Result<Option<Ref>> {
    let word = w.word()?;
    if word.starts_with('*') {
        return Ok(None);
    }
    let orientation = match word.as_bytes()[0] {
        b'+' => Orientation::Forward,
        b'-' => Orientation::Reversed,
        b'i' => Orientation::Internal,
        b'e' => Orientation::External,
        _ => return Err(bad(format!("bad shape reference {word}"))),
    };
    let number: usize = word[1..].parse().map_err(|_| bad(format!("bad shape reference {word}")))?;
    if number == 0 || number > count {
        return Err(bad(format!("shape reference {word} out of range")));
    }
    let location = w.count()?;
    Ok(Some(Ref { index: count - number, orientation, location }))
}

fn read_shape(w: &mut Words, version: FormatVersion, count: usize) -> Result<RawShape> {
    let kind = match w.word()? {
        "Ve" => NodeKind::Vertex,
        "Ed" => NodeKind::Edge,
        "Wi" => NodeKind::Wire,
        "Fa" => NodeKind::Face,
        "Sh" => NodeKind::Shell,
        "So" => NodeKind::Solid,
        "CS" => NodeKind::CompSolid,
        "Co" => NodeKind::Compound,
        k => return Err(bad(format!("unknown shape type {k}"))),
    };
    let geometry = match kind {
        NodeKind::Vertex => {
            let (tolerance, point) = (w.real()?, w.point()?);
            // Parameters of the point on curves and surfaces: not needed.
            loop {
                w.real()?;
                match w.int()? {
                    0 => break,
                    1 => {
                        w.count()?;
                    }
                    2 => {
                        w.count()?;
                        w.count()?;
                    }
                    3 => {
                        w.real()?;
                        w.count()?;
                    }
                    t => return Err(bad(format!("unknown point representation {t}"))),
                }
                w.count()?;
            }
            Geometry::Vertex { point, tolerance }
        }
        NodeKind::Edge => {
            let tolerance = w.real()?;
            let (_same_parameter, _same_range, degenerate) = (w.flag()?, w.flag()?, w.flag()?);
            let (mut curve, mut on_surfaces) = (None, Vec::new());
            loop {
                match w.int()? {
                    0 => break,
                    1 => {
                        let (c, l) = (w.count()?, w.count()?);
                        curve = Some((c, l, w.real()?, w.real()?));
                    }
                    rep @ (2 | 3) => {
                        let curve = w.count()?;
                        let reversed_curve = if rep == 3 {
                            let c = w.count()?;
                            w.word()?; // continuity
                            Some(c)
                        } else {
                            None
                        };
                        let (surface, location) = (w.count()?, w.count()?);
                        let (first, last) = (w.real()?, w.real()?);
                        if version == FormatVersion::V2 {
                            for _ in 0..4 {
                                w.real()?;
                            }
                        }
                        on_surfaces.push(OnSurface { curve, reversed_curve, surface, location, first, last });
                    }
                    4 => {
                        // Continuity across two faces.
                        w.word()?;
                        for _ in 0..4 {
                            w.count()?;
                        }
                    }
                    5 => {
                        w.count()?;
                        w.count()?;
                    }
                    6 => {
                        for _ in 0..3 {
                            w.count()?;
                        }
                    }
                    7 => {
                        for _ in 0..4 {
                            w.count()?;
                        }
                    }
                    t => return Err(bad(format!("unknown edge representation {t}"))),
                }
            }
            Geometry::Edge { tolerance, degenerate, curve, on_surfaces }
        }
        NodeKind::Face => {
            let restriction = w.int()?;
            if restriction == 2 {
                w.count()?;
                Geometry::Face { tolerance: 0.0, surface: None }
            } else {
                let tolerance = w.real()?;
                let surface = (w.count()?, w.count()?);
                Geometry::Face { tolerance, surface: Some(surface) }
            }
        }
        _ => Geometry::None,
    };
    // An optional triangulation of the face (`2 index`), then the flags.
    let mut flags = w.word()?;
    if kind == NodeKind::Face && flags == "2" {
        w.count()?;
        flags = w.word()?;
    }
    if flags.len() < 7 || !flags.bytes().all(|b| b == b'0' || b == b'1') {
        return Err(bad(format!("bad shape flags {flags}")));
    }
    let mut children = Vec::new();
    while let Some(r) = read_ref(w, count)? {
        children.push(r);
    }
    Ok(RawShape { kind, geometry, children })
}

// --- building -----------------------------------------------------------------------------

/// A face's surface (indices of the surface and its location) and the location of a wire
/// on it.
#[derive(Debug, Clone, Copy)]
struct OnFace {
    surface: usize,
    location: usize,
    wire_location: usize,
}

/// True if two transformations move points alike.
fn same_transform(a: &Transform, b: &Transform) -> bool {
    [Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), Point::new(0.0, 0.0, 1.0)].iter().all(|p| {
        let (q, r) = (p.transformed(a), p.transformed(b));
        q.distance(r) <= 1e-9 * (1.0 + (q - Point::ORIGIN).length())
    })
}

/// The locations met from the root down to a shape (the indices of those other than the
/// identity).
type Path = Vec<usize>;

/// Shapes built, by index and placement: a shape reached along two paths of locations that
/// compose to the same placement is built once.
struct Memo<T>(HashMap<usize, Vec<(Transform, T)>>);

impl<T> Default for Memo<T> {
    fn default() -> Self {
        Memo(HashMap::new())
    }
}

impl<T: Copy> Memo<T> {
    fn get(&self, index: usize, t: &Transform) -> Option<T> {
        self.0.get(&index)?.iter().find(|(u, _)| same_transform(u, t)).map(|&(_, x)| x)
    }

    fn insert(&mut self, index: usize, t: Transform, x: T) {
        self.0.entry(index).or_default().push((t, x));
    }
}

struct Builder<'a> {
    model: &'a mut Model,
    file: &'a File,
    vertices: Memo<VertexId>,
    edges: Memo<EdgeId>,
    shapes: Memo<Shape>,
}

impl<'a> Builder<'a> {
    fn new(model: &'a mut Model, file: &'a File) -> Self {
        Self { model, file, vertices: Memo::default(), edges: Memo::default(), shapes: Memo::default() }
    }

    fn location(&self, index: usize) -> Result<Transform> {
        self.file.locations.get(index).copied().ok_or_else(|| bad(format!("no location {index}")))
    }

    fn transform(&self, path: &[usize]) -> Result<Transform> {
        path.iter().try_fold(Transform::IDENTITY, |t, &l| Ok(t * self.location(l)?))
    }

    fn raw(&self, index: usize) -> &'a RawShape {
        &self.file.shapes[index]
    }

    fn root(&mut self) -> Result<Shape> {
        let root = self.file.root;
        let shape = self.shape(Ref { location: 0, ..root }, &[])?;
        Ok(shape.with_placement(self.location(root.location)?))
    }

    /// The shape `r` refers to, below a parent reached by `path`.
    fn shape(&mut self, r: Ref, path: &[usize]) -> Result<Shape> {
        let mut path = path.to_vec();
        if r.location != 0 {
            path.push(r.location);
        }
        let shape = match self.raw(r.index).kind {
            NodeKind::Vertex => Shape::new(self.vertex(r.index, &path)?.into()),
            NodeKind::Edge => Shape::new(self.edge(r.index, &path)?.into()),
            NodeKind::Wire => self.wire_on_face(r.index, &path, Orientation::Forward, None)?,
            _ => {
                let t = self.transform(&path)?;
                match self.shapes.get(r.index, &t) {
                    Some(s) => s,
                    None => {
                        let s = self.container(r.index, &path)?;
                        self.shapes.insert(r.index, t, s);
                        s
                    }
                }
            }
        };
        Ok(Shape { orientation: r.orientation, ..shape })
    }

    fn vertex(&mut self, index: usize, path: &Path) -> Result<VertexId> {
        let t = self.transform(path)?;
        if let Some(v) = self.vertices.get(index, &t) {
            return Ok(v);
        }
        let Geometry::Vertex { point, tolerance } = self.raw(index).geometry else {
            return Err(bad("a vertex without a point"));
        };
        let v = self.model.add_vertex(point.transformed(&self.transform(path)?), tolerance);
        self.vertices.insert(index, t, v);
        Ok(v)
    }

    fn edge(&mut self, index: usize, path: &Path) -> Result<EdgeId> {
        let placed = self.transform(path)?;
        if let Some(e) = self.edges.get(index, &placed) {
            return Ok(e);
        }
        let raw = self.raw(index);
        let Geometry::Edge { tolerance, degenerate, curve, on_surfaces } = &raw.geometry else {
            return Err(bad("an edge without geometry"));
        };
        let (mut start, mut end) = (None, None);
        for child in &raw.children {
            let mut child_path = path.clone();
            if child.location != 0 {
                child_path.push(child.location);
            }
            match child.orientation {
                Orientation::Forward => start = Some(self.vertex(child.index, &child_path)?),
                Orientation::Reversed => end = Some(self.vertex(child.index, &child_path)?),
                // Vertices inside the edge: not kept.
                _ => {}
            }
        }
        let (start, end) = match (start, end) {
            (Some(s), Some(e)) => (s, e),
            (Some(v), None) | (None, Some(v)) => (v, v),
            (None, None) => todo!("BRep edges without vertices"),
        };
        let edge = self.model.add_edge(start, end, *tolerance)?;
        match curve {
            Some((c, l, first, last)) if !degenerate => {
                let t = self.transform(path)? * self.location(*l)?;
                let raw_curve = self.file.curves.get(c.wrapping_sub(1)).ok_or_else(|| bad(format!("no curve {c}")))?;
                self.model.set_edge_curve(edge.into(), raw_curve.build()?.transformed(&t), *tolerance)?;
                self.model.set_edge_range(edge.into(), *first, *last)?;
            }
            _ if *degenerate => {}
            _ => {
                // Only curves on surfaces: the 3D curve made from the first (as OCCT's
                // BRepLib::BuildCurves3d does).
                let o = on_surfaces.first().ok_or_else(|| bad("an edge without curves"))?;
                let t = self.transform(path)? * self.location(o.location)?;
                let raw_surface = self.file.surfaces.get(o.surface.wrapping_sub(1)).ok_or_else(|| bad(format!("no surface {}", o.surface)))?;
                let raw_curve = self.file.curves2.get(o.curve.wrapping_sub(1)).ok_or_else(|| bad(format!("no 2D curve {}", o.curve)))?;
                let curve = curve_on_surface(&raw_curve.build()?, &raw_surface.build()?.transformed(&t), o.first, o.last)?;
                self.model.set_edge_curve(edge.into(), curve, *tolerance)?;
                self.model.set_edge_range(edge.into(), o.first, o.last)?;
            }
        }
        self.edges.insert(index, placed, edge);
        Ok(edge)
    }

    /// The curve of an edge on a face's surface; for a seam, the one for `orientation`.
    fn set_pcurve(
        &mut self,
        coedge: CoEdgeId,
        edge: usize,
        orientation: Orientation,
        face: OnFace,
        expected: &Transform,
    ) -> Result<()> {
        let Geometry::Edge { on_surfaces, .. } = &self.raw(edge).geometry else {
            return Ok(());
        };
        let candidates: Vec<&OnSurface> = on_surfaces.iter().filter(|o| o.surface == face.surface).collect();
        let at = |o: &&&OnSurface| self.location(o.location).is_ok_and(|t| same_transform(&t, expected));
        // (On a plane, only a curve for the surface placed as expected is taken: OCCT finds
        // the curve on a plane from the edge's own where none is, as kentron does.)
        let plane = self.file.surfaces.get(face.surface.wrapping_sub(1)).and_then(|r| r.build().ok()).is_some_and(|s| matches!(s, Surface::Plane(_)));
        let Some(on) = candidates
            .iter()
            .find(at)
            .or_else(|| if plane { None } else { candidates.iter().find(|o| o.location == face.location).or(candidates.first()) })
        else {
            return Ok(());
        };
        let index = match (orientation, on.reversed_curve) {
            (Orientation::Reversed, Some(c)) => c,
            _ => on.curve,
        };
        let curve = self.file.curves2.get(index.wrapping_sub(1)).ok_or_else(|| bad(format!("no 2D curve {index}")))?;
        self.model.set_pcurve_with_range(coedge, curve.build()?, on.first, on.last)?;
        Ok(())
    }

    /// A face, shell, solid, compsolid or compound.
    fn container(&mut self, index: usize, path: &Path) -> Result<Shape> {
        let raw = self.raw(index);
        let node = self.model.new_empty(raw.kind);
        match (&raw.kind, &raw.geometry) {
            (NodeKind::Face, Geometry::Face { tolerance, surface }) => {
                let Some((s, l)) = *surface else { todo!("BRep faces with only a triangulation") };
                let t = self.transform(path)? * self.location(l)?;
                let raw_surface =
                    self.file.surfaces.get(s.wrapping_sub(1)).ok_or_else(|| bad(format!("no surface {s}")))?;
                self.model.set_face_surface(node, raw_surface.build()?.transformed(&t), *tolerance)?;
                for child in &raw.children {
                    let mut wire_path = path.clone();
                    if child.location != 0 {
                        wire_path.push(child.location);
                    }
                    // A wire's own orientation reverses its edges' uses on the face.
                    let on = OnFace { surface: s, location: l, wire_location: child.location };
                    let wire = self.wire_on_face(child.index, &wire_path, child.orientation, Some(on))?;
                    self.model.add(node, Shape { orientation: child.orientation, ..wire })?;
                }
            }
            (NodeKind::Compound, _) => {
                // Children keep their locations as placements.
                for child in &raw.children {
                    let shape = self.shape(Ref { location: 0, ..*child }, &[])?;
                    self.model.add(node, shape.with_placement(self.location(child.location)?))?;
                }
            }
            _ => {
                for child in &raw.children {
                    let shape = self.shape(*child, path)?;
                    self.model.add(node, shape)?;
                }
            }
        }
        Ok(node)
    }

    /// A wire; on a face (given the face's surface and the wire's orientation on it), with
    /// the edges' curves on the surface.
    fn wire_on_face(
        &mut self,
        index: usize,
        path: &Path,
        orientation: Orientation,
        surface: Option<OnFace>,
    ) -> Result<Shape> {
        let raw = self.raw(index);
        let mut coedges: Vec<CoEdgeId> = Vec::new();
        for child in &raw.children {
            let mut edge_path = path.clone();
            if child.location != 0 {
                edge_path.push(child.location);
            }
            let edge = self.edge(child.index, &edge_path)?;
            let coedge = self.model.add_coedge(edge, child.orientation)?;
            if let Some(on) = surface {
                // As OCCT: the curve on the surface placed where the face's surface is, seen
                // from the edge (both below the wire).
                let relative = self.location(on.wire_location)? * self.location(child.location)?;
                let expected = relative.inverse() * self.location(on.location)?;
                self.set_pcurve(coedge, child.index, child.orientation.compose(orientation), on, &expected)?;
            }
            coedges.push(coedge);
        }
        // OCCT lists a wire's edges in any order.
        let wire = self.model.add_wire(&coedges)?;
        self.model.set_coedge_order(wire, &coedges)?;
        Ok(Shape::new(wire.into()))
    }
}

// --- writing ------------------------------------------------------------------------------

/// Coordinates of a pole, as written.
trait Coordinates {
    fn coordinates(&self) -> Vec<f64>;
}

impl Coordinates for Point {
    fn coordinates(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }
}

impl Coordinates for Point2 {
    fn coordinates(&self) -> Vec<f64> {
        vec![self.x, self.y]
    }
}

/// Text in the BRep format's word-separated style.
#[derive(Default)]
struct Text(String);

impl Text {
    fn word(&mut self, w: impl std::fmt::Display) -> &mut Self {
        use std::fmt::Write as _;
        let _ = write!(self.0, "{w} ");
        self
    }

    fn line(&mut self) -> &mut Self {
        self.0.push('\n');
        self
    }

    /// Reals, in exponent form when plain digits would be long: OCCT reads each number
    /// into a short buffer.
    fn reals(&mut self, values: &[f64]) -> &mut Self {
        for &v in values {
            self.real(v);
        }
        self
    }

    fn real(&mut self, v: f64) -> &mut Self {
        if v == 0.0 || (1e-4..1e15).contains(&v.abs()) {
            self.word(v)
        } else {
            self.word(format_args!("{v:e}"))
        }
    }

    fn point(&mut self, p: impl Coordinates) -> &mut Self {
        self.reals(&p.coordinates())
    }

    fn direction(&mut self, d: Direction) -> &mut Self {
        self.reals(&[d.x(), d.y(), d.z()])
    }

    fn direction2(&mut self, d: Direction2) -> &mut Self {
        self.reals(&[d.x(), d.y()])
    }

    /// `gp_Ax3` (of surfaces): origin, main, X and Y directions, left-handed if the frame is.
    fn frame(&mut self, f: &Frame) -> &mut Self {
        self.point(f.origin()).direction(f.z()).direction(f.x()).direction(f.y())
    }

    /// `gp_Ax2` (of curves), always right-handed: a left-handed frame (a reflected conic's)
    /// is written with the opposite main direction, which keeps the parametrization.
    fn right_handed_frame(&mut self, f: &Frame) -> &mut Self {
        let z = Direction::from_vector(Vector::from(f.x()).cross(Vector::from(f.y())));
        self.point(f.origin()).direction(z).direction(f.x()).direction(f.y())
    }

    fn frame2(&mut self, f: &Frame2) -> &mut Self {
        self.point(f.origin()).direction2(f.x()).direction2(f.y())
    }

    fn flag(&mut self, b: bool) -> &mut Self {
        self.word(u8::from(b))
    }

    fn knots(&mut self, knots: &[f64], multiplicities: &[usize]) -> &mut Self {
        for (k, m) in knots.iter().zip(multiplicities) {
            self.real(*k).word(m);
        }
        self
    }

    fn poles<P: Coordinates + Copy>(&mut self, poles: &[P], weights: Option<&[f64]>) -> &mut Self {
        for (i, p) in poles.iter().enumerate() {
            self.point(*p);
            if let Some(w) = weights {
                self.real(w[i]);
            }
        }
        self
    }

    fn bezier<P: ControlPoint + Coordinates>(&mut self, c: &BezierCurve<P>) -> &mut Self {
        self.flag(c.is_rational()).word(c.degree()).poles(c.poles(), c.weights())
    }

    fn bspline<P: ControlPoint + Coordinates>(&mut self, c: &BSplineCurve<P>) -> &mut Self {
        self.flag(c.is_rational()).flag(c.is_periodic()).word(c.degree());
        self.word(c.poles().len()).word(c.knots().len());
        self.poles(c.poles(), c.weights()).knots(c.knots(), c.multiplicities())
    }

    fn curve(&mut self, c: &Curve) -> &mut Self {
        match c {
            Curve::Line(l) => self.word(1).point(l.origin).direction(l.direction),
            Curve::Circle(c) => self.word(2).right_handed_frame(c.frame()).real(c.radius()),
            Curve::Ellipse(e) => self.word(3).right_handed_frame(e.frame()).real(e.major_radius()).real(e.minor_radius()),
            Curve::Parabola(p) => self.word(4).right_handed_frame(p.frame()).real(p.focal()),
            Curve::Hyperbola(h) => self.word(5).right_handed_frame(h.frame()).real(h.major_radius()).real(h.minor_radius()),
            Curve::Bezier(b) => self.word(6).bezier(b),
            Curve::BSpline(b) => self.word(7).bspline(b),
            Curve::Trimmed(t) => {
                let (u1, u2) = t.range();
                self.word(8).real(u1).real(u2).line().curve(t.basis())
            }
            Curve::Offset(o) => self.word(9).real(o.offset()).direction(o.direction()).line().curve(o.basis()),
        }
    }

    fn curve2(&mut self, c: &Curve2) -> &mut Self {
        match c {
            Curve2::Line(l) => self.word(1).point(l.origin).direction2(l.direction),
            Curve2::Circle(c) => self.word(2).frame2(c.frame()).real(c.radius()),
            Curve2::Ellipse(e) => self.word(3).frame2(e.frame()).real(e.major_radius()).real(e.minor_radius()),
            Curve2::Parabola(p) => self.word(4).frame2(p.frame()).real(p.focal()),
            Curve2::Hyperbola(h) => self.word(5).frame2(h.frame()).real(h.major_radius()).real(h.minor_radius()),
            Curve2::Bezier(b) => self.word(6).bezier(b),
            Curve2::BSpline(b) => self.word(7).bspline(b),
            Curve2::Trimmed(t) => {
                let (u1, u2) = t.range();
                self.word(8).real(u1).real(u2).line().curve2(t.basis())
            }
            Curve2::Offset(o) => self.word(9).real(o.offset()).line().curve2(o.basis()),
        }
    }

    fn grid(&mut self, poles: &[Vec<Point>], weights: Option<&[Vec<f64>]>) -> &mut Self {
        for (i, row) in poles.iter().enumerate() {
            self.poles(row, weights.map(|w| w[i].as_slice()));
        }
        self
    }

    fn surface(&mut self, s: &Surface) -> &mut Self {
        match s {
            Surface::Plane(p) => self.word(1).frame(p.frame()),
            Surface::Cylinder(c) => self.word(2).frame(c.frame()).real(c.radius()),
            Surface::Cone(c) => self.word(3).frame(c.frame()).real(c.radius()).real(c.semi_angle()),
            Surface::Sphere(s) => self.word(4).frame(s.frame()).real(s.radius()),
            Surface::Torus(t) => self.word(5).frame(t.frame()).real(t.major_radius()).real(t.minor_radius()),
            Surface::Extrusion(e) => self.word(6).direction(e.direction()).line().curve(e.basis()),
            Surface::Revolution(r) => {
                let axis = r.axis();
                self.word(7).point(axis.origin).direction(axis.direction).line().curve(r.basis())
            }
            Surface::Bezier(b) => {
                let rational = b.weights().is_some();
                self.word(8).flag(rational).flag(rational).word(b.u_degree()).word(b.v_degree());
                self.grid(b.poles(), b.weights())
            }
            Surface::BSpline(b) => {
                let (u, v) = (b.u_knot_vector(), b.v_knot_vector());
                let rational = b.weights().is_some();
                self.word(9).flag(rational).flag(rational).flag(u.periodic).flag(v.periodic);
                self.word(u.degree).word(v.degree).word(b.u_pole_count()).word(b.v_pole_count());
                self.word(u.knots.len()).word(v.knots.len());
                self.grid(b.poles(), b.weights());
                self.knots(&u.knots, &u.multiplicities).knots(&v.knots, &v.multiplicities)
            }
            Surface::Trimmed(t) => {
                let (u1, u2, v1, v2) = t.trim();
                self.word(10).reals(&[u1, u2, v1, v2]).line().surface(t.basis())
            }
            Surface::Offset(o) => self.word(11).real(o.offset()).line().surface(o.basis()),
        }
    }

    fn reference(&mut self, orientation: Orientation, number: usize, location: usize) -> &mut Self {
        let sign = match orientation {
            Orientation::Forward => '+',
            Orientation::Reversed => '-',
            Orientation::Internal => 'i',
            Orientation::External => 'e',
        };
        let _ = std::fmt::Write::write_fmt(&mut self.0, format_args!("{sign}{number} {location} "));
        self
    }
}

/// A child as written: the node, its orientation and the index of its location.
type Child = (NodeId, Orientation, usize);

/// Gathers the tables of a file, then writes them.
struct Writer<'a> {
    model: &'a Model,
    version: FormatVersion,
    /// Locations other than the identity; location `i` is `locations[i - 1]`.
    locations: Vec<Transform>,
    curves2: Text,
    curve2_count: usize,
    curves: Text,
    curve_count: usize,
    surfaces: Text,
    surface_count: usize,
    /// Shapes, each after its sub-shapes, with their children.
    shapes: Vec<(NodeId, Vec<Child>)>,
    numbers: HashMap<NodeId, usize>,
    face_surfaces: HashMap<FaceId, usize>,
}

impl<'a> Writer<'a> {
    fn new(model: &'a Model, version: FormatVersion) -> Self {
        Self {
            model,
            version,
            locations: Vec::new(),
            curves2: Text::default(),
            curve2_count: 0,
            curves: Text::default(),
            curve_count: 0,
            surfaces: Text::default(),
            surface_count: 0,
            shapes: Vec::new(),
            numbers: HashMap::new(),
            face_surfaces: HashMap::new(),
        }
    }

    fn location(&mut self, t: Transform) -> usize {
        if t == Transform::IDENTITY {
            return 0;
        }
        match self.locations.iter().position(|l| *l == t) {
            Some(i) => i + 1,
            None => {
                self.locations.push(t);
                self.locations.len()
            }
        }
    }

    fn add_curve(&mut self, c: &Curve) -> usize {
        self.curves.curve(c).line();
        self.curve_count += 1;
        self.curve_count
    }

    fn add_curve2(&mut self, c: &Curve2) -> usize {
        self.curves2.curve2(c).line();
        self.curve2_count += 1;
        self.curve2_count
    }

    /// Adds `node` and everything below it to the shape table.
    fn collect(&mut self, node: NodeId) -> Result<()> {
        if self.numbers.contains_key(&node) {
            return Ok(());
        }
        let model = self.model;
        let mut children = Vec::new();
        match node.kind {
            NodeKind::Wire => {
                let wire = model.get(WireId::try_from(node).map_err(|_| bad("not a wire"))?).ok_or_else(|| bad("no wire"))?;
                for &c in &wire.coedges {
                    let coedge = model.get(c).ok_or_else(|| bad("no coedge"))?;
                    children.push((coedge.edge.into(), coedge.orientation, 0));
                }
            }
            NodeKind::CoEdge | NodeKind::Product | NodeKind::Instance => {
                return Err(Error::InvalidInput(format!("BRep files cannot hold a {:?}", node.kind)));
            }
            _ => {
                for child in model.children_local(Shape::new(node)) {
                    let location = self.location(child.placement);
                    children.push((child.node, child.orientation, location));
                }
            }
        }
        for &(child, _, _) in &children {
            self.collect(child)?;
        }
        self.numbers.insert(node, self.shapes.len());
        self.shapes.push((node, children));
        Ok(())
    }

    fn file(mut self, root: Shape) -> Result<String> {
        self.collect(root.node)?;
        let root_location = self.location(root.placement);
        let model = self.model;
        // Surfaces, then the curves of edges and of their uses on faces.
        let faces: Vec<FaceId> =
            self.shapes.iter().filter_map(|(n, _)| FaceId::try_from(*n).ok()).collect();
        // Per edge: (face, the coedge's orientation on the face, coedge).
        let mut uses: HashMap<EdgeId, Vec<(FaceId, Orientation, CoEdgeId)>> = HashMap::new();
        for &face in &faces {
            let surface = model.surface(Shape::new(face.into())).ok_or_else(|| bad("a face without a surface"))?;
            self.surfaces.surface(surface).line();
            self.surface_count += 1;
            self.face_surfaces.insert(face, self.surface_count);
            let record = model.get(face).ok_or_else(|| bad("no face"))?;
            for &use_id in &record.wires {
                let wire_use = model.get(use_id).ok_or_else(|| bad("no wire use"))?;
                let wire = model.get(wire_use.wire).ok_or_else(|| bad("no wire"))?;
                for &c in &wire.coedges {
                    let coedge = model.get(c).ok_or_else(|| bad("no coedge"))?;
                    uses.entry(coedge.edge).or_default().push((face, coedge.orientation.compose(wire_use.orientation), c));
                }
            }
        }
        let mut shapes = Text::default();
        shapes.word("TShapes").word(self.shapes.len()).line();
        let count = self.shapes.len();
        let shape_list = std::mem::take(&mut self.shapes);
        for (node, children) in &shape_list {
            let shape = Shape::new(*node);
            match node.kind {
                NodeKind::Vertex => {
                    let point = model.point(shape);
                    shapes.word("Ve").line().real(model.tolerance(shape)).line().point(point).line();
                    shapes.word(0).word(0).line().line().word("0101101").line();
                }
                NodeKind::Edge => {
                    let degenerate = model.is_degenerate(shape);
                    let edge = EdgeId::try_from(*node).map_err(|_| bad("not an edge"))?;
                    let edge_uses = uses.remove(&edge).unwrap_or_default();
                    let same_parameter = edge_uses.iter().all(|u| model.is_same_parameter(u.2));
                    let same_range = edge_uses.iter().all(|u| model.is_same_range(u.2));
                    shapes.word("Ed").line().real(model.tolerance(shape));
                    shapes.flag(same_parameter).flag(same_range).flag(degenerate).line();
                    if let Some((curve, first, last)) = model.edge_curve(shape) {
                        let c = self.add_curve(&curve);
                        shapes.word(1).word(c).word(0).real(first).real(last).line();
                    }
                    let mut done = vec![false; edge_uses.len()];
                    for i in 0..edge_uses.len() {
                        if done[i] {
                            continue;
                        }
                        let (face, orientation, coedge) = edge_uses[i];
                        let surface = self.face_surfaces[&face];
                        // A seam: the other use on the same face, the other way.
                        let seam = (i + 1..edge_uses.len())
                            .find(|&j| !done[j] && edge_uses[j].0 == face && edge_uses[j].1 != orientation);
                        let Some((pcurve, first, last)) = model.effective_pcurve(coedge) else { continue };
                        match seam {
                            Some(j) => {
                                done[j] = true;
                                let Some((other, _, _)) = model.effective_pcurve(edge_uses[j].2) else { continue };
                                let (forward, reversed) =
                                    if orientation == Orientation::Reversed { (other, pcurve.clone()) } else { (pcurve.clone(), other) };
                                let (c1, c2) = (self.add_curve2(&forward), self.add_curve2(&reversed));
                                shapes.word(3).word(c1).word(c2).word("CN").word(surface).word(0).real(first).real(last);
                            }
                            None => {
                                let c = self.add_curve2(&pcurve);
                                shapes.word(2).word(c).word(surface).word(0).real(first).real(last);
                            }
                        }
                        if self.version == FormatVersion::V2 {
                            shapes.point(pcurve.point(first)).point(pcurve.point(last));
                        }
                        shapes.line();
                    }
                    let closed = model.is_closed_edge(shape);
                    shapes.word(0).line().line().word(if closed { "0101100" } else { "0101000" }).line();
                }
                NodeKind::Wire => {
                    let closed = model.is_closed_wire(WireId::try_from(*node).map_err(|_| bad("not a wire"))?);
                    shapes.word("Wi").line().line().word(if closed { "0101100" } else { "0101000" }).line();
                }
                NodeKind::Face => {
                    let face = FaceId::try_from(*node).map_err(|_| bad("not a face"))?;
                    let surface = self.face_surfaces[&face];
                    shapes.word("Fa").line().word(0).real(model.tolerance(shape)).word(surface).word(0).line();
                    shapes.line().word("0101000").line();
                }
                NodeKind::Shell => {
                    let closed = model.is_closed_shell(ShellId::try_from(*node).map_err(|_| bad("not a shell"))?);
                    shapes.word("Sh").line().line().word(if closed { "0101100" } else { "0101000" }).line();
                }
                NodeKind::Solid => {
                    shapes.word("So").line().line().word("0100000").line();
                }
                NodeKind::CompSolid => {
                    shapes.word("CS").line().line().word("0101000").line();
                }
                NodeKind::Compound => {
                    shapes.word("Co").line().line().word("0101000").line();
                }
                _ => unreachable!("not collected"),
            }
            for &(child, orientation, location) in children {
                shapes.reference(orientation, count - self.numbers[&child], location);
            }
            shapes.word("*").line().line();
        }
        shapes.reference(root.orientation, count - self.numbers[&root.node], root_location).line();

        let header = VERSIONS.iter().find(|(_, v)| *v == self.version).map(|(h, _)| *h).unwrap_or(VERSIONS[2].0);
        let mut text = Text::default();
        text.0.push_str(header);
        text.line().word("Locations").word(self.locations.len()).line();
        for l in &self.locations {
            text.word(1).line();
            for row in l.rows() {
                text.reals(&row).line();
            }
        }
        text.word("Curve2ds").word(self.curve2_count).line();
        text.0.push_str(&self.curves2.0);
        text.word("Curves").word(self.curve_count).line();
        text.0.push_str(&self.curves.0);
        text.word("Polygon3D").word(0).line().word("PolygonOnTriangulations").word(0).line();
        text.word("Surfaces").word(self.surface_count).line();
        text.0.push_str(&self.surfaces.0);
        text.word("Triangulations").word(0).line().line();
        text.0.push_str(&shapes.0);
        Ok(text.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A line segment between two vertices, as OCCT writes it.
    #[test]
    fn reads_an_edge() {
        let text = "DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\nLocations 1\n1\n\
            1 0 0 5\n0 1 0 0\n0 0 1 0\nCurve2ds 0\nCurves 1\n1 0 0 0 1 0 0 \nPolygon3D 0\n\
            PolygonOnTriangulations 0\nSurfaces 0\nTriangulations 0\n\nTShapes 3\n\
            Ve\n1e-07\n0 0 0\n0 0\n\n0101101\n*\n\
            Ve\n1e-07\n10 0 0\n0 0\n\n0101101\n*\n\
            Ed\n 1e-07 1 1 0\n1  1 0 0 10\n0\n\n0101000\n+3 0 -2 0 *\n\n+1 1 ";
        let mut model = Model::new();
        let edge = read(&mut model, &mut text.as_bytes()).unwrap();
        assert_eq!(edge.kind(), NodeKind::Edge);
        let (a, b) = model.edge_vertices(edge);
        assert_eq!(model[crate::topology::VertexId::try_from(a.node).unwrap()].point, Point::ORIGIN);
        assert_eq!(model[crate::topology::VertexId::try_from(b.node).unwrap()].point, Point::new(10.0, 0.0, 0.0));
        assert!(Point::ORIGIN.transformed(&edge.placement).distance(Point::new(5.0, 0.0, 0.0)) < 1e-12);
    }
    /// Solids written and read back keep their volumes and stay valid (a cylinder has a
    /// seam; a compound child keeps its placement).
    #[test]
    fn writes_solids_that_read_back() {
        use crate::primitives::{Cuboid, Cylinder};
        let mut model = Model::new();
        let cuboid = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let cylinder = Cylinder::new(1.0, 2.0).build(&mut model).unwrap();
        let moved = cylinder.with_placement(Transform::translation(Vector::new(5.0, 0.0, 0.0)));
        let compound = model.new_empty(NodeKind::Compound);
        model.add(compound, cuboid).unwrap();
        model.add(compound, moved).unwrap();
        let mut bytes = Vec::new();
        write(&model, compound, &mut bytes, &WriteOptions::default()).unwrap();
        let mut other = Model::new();
        let back = read(&mut other, &mut bytes.as_slice()).unwrap();
        assert!(crate::validity::check(&other, back).is_valid());
        let volume = crate::properties::volume(&other, back);
        assert!((volume.mass - (6.0 + 2.0 * std::f64::consts::PI)).abs() < 1e-9, "{}", volume.mass);
        assert!(volume.center_of_mass.x > 2.0);
    }
}
