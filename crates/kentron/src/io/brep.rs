//! OCCT's native ASCII BRep format (`.brep` files, as used by DRAW's `restore`/`save`).
//!
//! A file holds tables, each element referred to by its 1-based index: locations
//! (transformations), 2D curves, 3D curves, 3D polygons, polygons on triangulations,
//! surfaces, triangulations, then the shapes ("TShapes"), sub-shapes before the shapes using
//! them, and last a reference to the root shape. A reference to a shape is its orientation
//! sign, its number counted back from the end of the table, and a location.
//!
//! OCCT's binary BRep format (`BinTools`, DRAW's `binsave`/`binrestore`) holds the same
//! tables, its numbers in binary ([`read_binary`]).
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
    let text = Writer::<Text>::new(model, options.version).file(shape)?;
    out.write_all(text.as_bytes()).map_err(|e| Error::InvalidInput(format!("writing BRep: {e}")))
}

/// Writes `shape` and everything below it in OCCT's binary BRep format (`BinTools`, version
/// 3, without meshes).
#[doc(alias = "BinTools::Write")]
pub fn write_binary(model: &Model, shape: Shape, out: &mut dyn Write) -> Result<()> {
    let bytes = Writer::<Bin>::new(model, FormatVersion::V3).file(shape)?;
    out.write_all(&bytes).map_err(|e| Error::InvalidInput(format!("writing BRep: {e}")))
}

/// Reads a shape into `model`; any format version is accepted.
#[doc(alias = "BRepTools::Read")]
pub fn read(model: &mut Model, input: &mut dyn Read) -> Result<Shape> {
    let mut text = String::new();
    input.read_to_string(&mut text).map_err(|e| Error::InvalidInput(format!("reading BRep: {e}")))?;
    let file = parse(&text)?;
    Builder::new(model, &file).root()
}

/// Reads a shape saved in OCCT's binary BRep format (`BinTools::Write`, DRAW's `binsave`)
/// into `model` (any version).
#[doc(alias = "BinTools::Read")]
pub fn read_binary(model: &mut Model, input: &mut dyn Read) -> Result<Shape> {
    let mut data = Vec::new();
    input.read_to_end(&mut data).map_err(|e| Error::InvalidInput(format!("reading BRep: {e}")))?;
    let file = parse_binary(&data)?;
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
}

/// What the tables are read from: the text format's words, or the binary format's numbers
/// ([`Bytes`]: type codes and flags there single bytes, degrees two).
trait Tokens {
    fn real(&mut self) -> Result<f64>;

    fn int(&mut self) -> Result<i64>;

    /// A geometry's type code.
    fn code(&mut self) -> Result<i64>;

    fn flag(&mut self) -> Result<bool>;

    /// A curve's or surface's degree.
    fn degree(&mut self) -> Result<usize>;

    fn count(&mut self) -> Result<usize> {
        usize::try_from(self.int()?).map_err(|_| bad("negative count"))
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

impl Tokens for Words<'_> {
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

    fn code(&mut self) -> Result<i64> {
        self.int()
    }

    fn flag(&mut self) -> Result<bool> {
        Ok(self.int()? != 0)
    }

    fn degree(&mut self) -> Result<usize> {
        self.count()
    }
}

fn pole3<T: Tokens>(w: &mut T) -> Result<Point> {
    w.point()
}

fn pole2<T: Tokens>(w: &mut T) -> Result<Point2> {
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
    fn read<T: Tokens>(w: &mut T, pole: impl Fn(&mut T) -> Result<P>) -> Result<Self> {
        let (rational, periodic) = (w.flag()?, w.flag()?);
        let (degree, pole_count, knot_count) = (w.degree()?, w.count()?, w.count()?);
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

fn bezier<P: ControlPoint, T: Tokens>(w: &mut T, pole: impl Fn(&mut T) -> Result<P>) -> Result<BezierCurve<P>> {
    let rational = w.flag()?;
    let degree = w.degree()?;
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
    fn read<T: Tokens>(w: &mut T) -> Result<Self> {
        Ok(match w.code()? {
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
    fn read<T: Tokens>(w: &mut T) -> Result<Self> {
        Ok(match w.code()? {
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
    fn read<T: Tokens>(w: &mut T) -> Result<Self> {
        Ok(match w.code()? {
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
                let (u_degree, v_degree) = (w.degree()?, w.degree()?);
                let (poles, weights) = grid(w, u_degree + 1, v_degree + 1, u_rational || v_rational)?;
                RawSurface::Ready(Surface::Bezier(match weights {
                    Some(weights) => BezierSurface::rational(poles, weights),
                    None => BezierSurface::new(poles),
                }))
            }
            9 => {
                let (u_rational, v_rational) = (w.flag()?, w.flag()?);
                let (u_periodic, v_periodic) = (w.flag()?, w.flag()?);
                let (u_degree, v_degree) = (w.degree()?, w.degree()?);
                let (u_poles, v_poles) = (w.count()?, w.count()?);
                let (u_knots, v_knots) = (w.count()?, w.count()?);
                let (poles, weights) = grid(w, u_poles, v_poles, u_rational || v_rational)?;
                fn knot_vector<T: Tokens>(w: &mut T, count: usize, degree: usize, periodic: bool) -> Result<KnotVector> {
                    let (knots, mults) = w.knots(count)?;
                    Ok(if periodic {
                        KnotVector::periodic(knots, mults, degree)
                    } else {
                        KnotVector::new(knots, mults, degree)
                    })
                }
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
fn grid<T: Tokens>(w: &mut T, nu: usize, nv: usize, rational: bool) -> Result<(Vec<Vec<Point>>, Option<Vec<Vec<f64>>>)> {
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
pub(crate) fn curve_on_surface(pcurve: &Curve2, surface: &Surface, first: f64, last: f64) -> Result<Curve> {
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

// --- the binary format --------------------------------------------------------------------

/// The binary format's header lines, by version (`BinTools_FormatVersion`: 2 adds the
/// (u, v) end points of curves on surfaces, 3 writes points on curves and surfaces
/// unambiguously, 4 adds meshes' normals and drops the end points again).
const BINARY_VERSIONS: [(&str, u8); 4] = [
    ("Open CASCADE Topology V1 (c)", 1),
    ("Open CASCADE Topology V2 (c)", 2),
    ("Open CASCADE Topology V3 (c)", 3),
    ("Open CASCADE Topology V4, (c) Open Cascade", 4),
];

/// The binary format's contents: little-endian numbers (reals of 8 bytes, integers of 4,
/// degrees of 2, type codes and flags of 1), each table headed by its name and size in
/// words, as the text format's.
struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
}

impl<'a> Bytes<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let bytes = self.data.get(self.at..self.at + n).ok_or_else(|| bad("unexpected end of file"))?;
        self.at += n;
        Ok(bytes)
    }

    fn byte(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn skip_space(&mut self) {
        while self.data.get(self.at).is_some_and(|b| b.is_ascii_whitespace()) {
            self.at += 1;
        }
    }

    /// A word, as C++ streams read one (`>>`): after any white space, up to the next.
    fn word(&mut self) -> Result<&'a str> {
        self.skip_space();
        let start = self.at;
        while self.data.get(self.at).is_some_and(|b| !b.is_ascii_whitespace()) {
            self.at += 1;
        }
        if start == self.at {
            return Err(bad("unexpected end of file"));
        }
        std::str::from_utf8(&self.data[start..self.at]).map_err(|_| bad("a table header that is not text"))
    }

    /// A table's header: its name and size, then the end of the line.
    fn table(&mut self, name: &str) -> Result<usize> {
        let word = self.word()?;
        if word != name {
            return Err(bad(format!("expected {name}, found {word}")));
        }
        let count = self.word()?.parse().map_err(|_| bad(format!("a bad size of the {name} table")))?;
        self.byte()?;
        Ok(count)
    }

    fn skip(&mut self, n: usize) -> Result<()> {
        self.take(n).map(|_| ())
    }
}

impl Tokens for Bytes<'_> {
    fn real(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().expect("8 bytes")))
    }

    fn int(&mut self) -> Result<i64> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes")) as i64)
    }

    fn code(&mut self) -> Result<i64> {
        Ok(self.byte()? as i64)
    }

    fn flag(&mut self) -> Result<bool> {
        Ok(self.byte()? != 0)
    }

    fn degree(&mut self) -> Result<usize> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().expect("2 bytes")) as usize)
    }
}

fn parse_binary(data: &[u8]) -> Result<File> {
    // The header line naming the version (the format's lines before it skipped).
    let mut at = 0;
    let version = loop {
        let end = data[at..].iter().position(|&b| b == b'\n').map_or(data.len(), |i| at + i);
        let line = String::from_utf8_lossy(&data[at..end]);
        let line = line.trim_end_matches('\r');
        if let Some(&(_, v)) = BINARY_VERSIONS.iter().find(|(h, _)| *h == line) {
            at = end + 1;
            break v;
        }
        if end >= data.len() {
            return Err(bad("no Open CASCADE Topology header"));
        }
        at = end + 1;
    };
    let mut b = Bytes { data, at };
    // The locations: elementary (`1` and a 3×4 matrix) or products of powers of earlier ones
    // (`2`, then index and power pairs ending with index 0); identities not numbered.
    let mut locations = vec![Transform::IDENTITY];
    for _ in 0..b.table("Locations")? {
        let t = match b.byte()? {
            1 => {
                let mut rows = [[0.0; 4]; 3];
                for row in &mut rows {
                    for v in row.iter_mut() {
                        *v = b.real()?;
                    }
                }
                Some(Transform::from_rows(rows).ok_or_else(|| bad("a singular location"))?)
            }
            2 => {
                let mut product = Transform::IDENTITY;
                let mut empty = true;
                loop {
                    let index = b.count()?;
                    if index == 0 {
                        break;
                    }
                    let power = b.int()? as i32;
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
    let curves2 = (0..b.table("Curve2ds")?).map(|_| RawCurve2::read(&mut b)).collect::<Result<Vec<_>>>()?;
    let curves = (0..b.table("Curves")?).map(|_| RawCurve::read(&mut b)).collect::<Result<Vec<_>>>()?;
    // Meshes, skipped: 3D polygons, polygons on triangulations.
    for _ in 0..b.table("Polygon3D")? {
        let (nodes, has_parameters) = (b.count()?, b.flag()?);
        b.skip(8 * (1 + 3 * nodes + if has_parameters { nodes } else { 0 }))?;
    }
    for _ in 0..b.table("PolygonOnTriangulations")? {
        let nodes = b.count()?;
        b.skip(4 * nodes + 8)?;
        if b.flag()? {
            b.skip(8 * nodes)?;
        }
    }
    let surfaces = (0..b.table("Surfaces")?).map(|_| RawSurface::read(&mut b)).collect::<Result<Vec<_>>>()?;
    for _ in 0..b.table("Triangulations")? {
        let (nodes, triangles, has_uv) = (b.count()?, b.count()?, b.flag()?);
        let has_normals = version >= 4 && b.flag()?;
        b.skip(8 + 8 * 3 * nodes + if has_uv { 8 * 2 * nodes } else { 0 } + 4 * 3 * triangles + if has_normals { 4 * 3 * nodes } else { 0 })?;
    }
    let count = b.table("TShapes")?;
    let mut shapes = Vec::with_capacity(count);
    for _ in 0..count {
        shapes.push(read_binary_shape(&mut b, version, count)?);
    }
    let root = read_binary_ref(&mut b, count)?.ok_or_else(|| bad("no root shape"))?;
    Ok(File { locations, curves2, curves, surfaces, shapes, root })
}

/// A shape reference (an orientation byte, the number counted back from the end of the
/// table and a location), or `None` for the end of a list (`*`).
fn read_binary_ref(b: &mut Bytes, count: usize) -> Result<Option<Ref>> {
    // (Read as a C++ stream reads a character: after any white space.)
    b.skip_space();
    let orientation = match b.byte()? {
        b'*' => return Ok(None),
        0 => Orientation::Forward,
        1 => Orientation::Reversed,
        2 => Orientation::Internal,
        3 => Orientation::External,
        o => return Err(bad(format!("bad shape orientation {o}"))),
    };
    let number = b.count()?;
    if number == 0 || number > count {
        return Err(bad(format!("shape reference {number} out of range")));
    }
    let location = b.count()?;
    Ok(Some(Ref { index: count - number, orientation, location }))
}

fn read_binary_shape(b: &mut Bytes, version: u8, count: usize) -> Result<RawShape> {
    // (`TopAbs_ShapeEnum`.)
    let kind = match b.byte()? {
        0 => NodeKind::Compound,
        1 => NodeKind::CompSolid,
        2 => NodeKind::Solid,
        3 => NodeKind::Shell,
        4 => NodeKind::Face,
        5 => NodeKind::Wire,
        6 => NodeKind::Edge,
        7 => NodeKind::Vertex,
        k => return Err(bad(format!("unknown shape type {k}"))),
    };
    let geometry = match kind {
        NodeKind::Vertex => {
            let (tolerance, point) = (b.real()?, b.point()?);
            // Parameters of the point on curves and surfaces: not needed. (Versions 3 and 4
            // give the kind first; older ones the parameter, or nothing before a kind of 0,
            // told apart as OCCT does.)
            loop {
                let kind = if version >= 3 {
                    let kind = b.byte()?;
                    if (1..=3).contains(&kind) {
                        b.real()?;
                    }
                    kind
                } else {
                    let at = b.at;
                    b.real()?;
                    let kind = b.byte()?;
                    if (1..=3).contains(&kind) {
                        kind
                    } else {
                        b.at = at;
                        let kind = b.byte()?;
                        if (1..=3).contains(&kind) {
                            b.real()?;
                        }
                        kind
                    }
                };
                match kind {
                    0 => break,
                    1 => b.skip(4)?,
                    2 => b.skip(8)?,
                    3 => b.skip(8 + 4)?,
                    t => return Err(bad(format!("unknown point representation {t}"))),
                }
                b.skip(4)?; // location
            }
            Geometry::Vertex { point, tolerance }
        }
        NodeKind::Edge => {
            let tolerance = b.real()?;
            let (_same_parameter, _same_range, degenerate) = (b.flag()?, b.flag()?, b.flag()?);
            let (mut curve, mut on_surfaces) = (None, Vec::new());
            loop {
                match b.byte()? {
                    0 => break,
                    1 => {
                        let (c, l) = (b.count()?, b.count()?);
                        curve = Some((c, l, b.real()?, b.real()?));
                    }
                    rep @ (2 | 3) => {
                        let curve = b.count()?;
                        let reversed_curve = if rep == 3 {
                            let c = b.count()?;
                            b.byte()?; // continuity
                            Some(c)
                        } else {
                            None
                        };
                        let (surface, location) = (b.count()?, b.count()?);
                        let (first, last) = (b.real()?, b.real()?);
                        if version == 2 || version == 3 {
                            b.skip(4 * 8)?; // the (u, v) end points
                        }
                        on_surfaces.push(OnSurface { curve, reversed_curve, surface, location, first, last });
                    }
                    // Continuity across two faces.
                    4 => b.skip(1 + 4 * 4)?,
                    5 => b.skip(4 * 2)?,
                    6 => b.skip(4 * 3)?,
                    7 => b.skip(4 * 4)?,
                    t => return Err(bad(format!("unknown edge representation {t}"))),
                }
            }
            Geometry::Edge { tolerance, degenerate, curve, on_surfaces }
        }
        NodeKind::Face => {
            let _natural_restriction = b.flag()?;
            let tolerance = b.real()?;
            let (surface, location) = (b.count()?, b.count()?);
            // An optional triangulation (`2` and its index).
            if b.byte()? == 2 {
                b.skip(4)?;
            }
            Geometry::Face { tolerance, surface: (surface > 0).then_some((surface, location)) }
        }
        _ => Geometry::None,
    };
    // The flags (free, modified, checked, orientable, closed, infinite, convex), then the
    // sub-shapes.
    b.skip(7)?;
    let mut children = Vec::new();
    while let Some(r) = read_binary_ref(b, count)? {
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
        // The 3D curve: the edge's own, else made from its first curve on a surface (as
        // OCCT's BRepLib::BuildCurves3d does); none for a degenerate edge.
        let geometry: Option<(Curve, f64, f64)> = match curve {
            Some((c, l, first, last)) if !degenerate => {
                let t = self.transform(path)? * self.location(*l)?;
                let raw_curve = self.file.curves.get(c.wrapping_sub(1)).ok_or_else(|| bad(format!("no curve {c}")))?;
                Some((raw_curve.build()?.transformed(&t), *first, *last))
            }
            _ if *degenerate => None,
            _ => {
                let o = on_surfaces.first().ok_or_else(|| bad("an edge without curves"))?;
                let t = self.transform(path)? * self.location(o.location)?;
                let raw_surface = self.file.surfaces.get(o.surface.wrapping_sub(1)).ok_or_else(|| bad(format!("no surface {}", o.surface)))?;
                let raw_curve = self.file.curves2.get(o.curve.wrapping_sub(1)).ok_or_else(|| bad(format!("no 2D curve {}", o.curve)))?;
                let curve = curve_on_surface(&raw_curve.build()?, &raw_surface.build()?.transformed(&t), o.first, o.last)?;
                Some((curve, o.first, o.last))
            }
        };
        let (start, end) = match (start, end) {
            (Some(s), Some(e)) => (s, e),
            (Some(v), None) | (None, Some(v)) => (v, v),
            // (An edge without vertices (OCCT allows them, on an infinite line, say): vertices
            // made at the ends of its range, one where they meet.)
            (None, None) => {
                let Some((c, first, last)) = &geometry else { return Err(bad("a degenerate edge without vertices")) };
                let (a, b) = (c.point(*first), c.point(*last));
                let start = self.model.add_vertex(a, *tolerance);
                let end = if a.distance(b) <= *tolerance { start } else { self.model.add_vertex(b, *tolerance) };
                (start, end)
            }
        };
        let edge = self.model.add_edge(start, end, *tolerance)?;
        if let Some((c, first, last)) = geometry {
            self.model.set_edge_curve(edge.into(), c, *tolerance)?;
            self.model.set_edge_range(edge.into(), first, last)?;
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

/// What the BRep formats write, typed: OCCT's text format writes each as a word, its binary
/// format each in its own width (reals 8 bytes, integers 4, degrees 2, type codes and flags 1,
/// little-endian). The geometry is written alike in both.
trait Out: Sized {
    fn real(&mut self, v: f64) -> &mut Self;
    fn int(&mut self, v: i64) -> &mut Self;
    /// A geometry's type code.
    fn code(&mut self, v: u8) -> &mut Self;
    /// A curve's or surface's degree.
    fn degree(&mut self, v: usize) -> &mut Self;
    fn flag(&mut self, b: bool) -> &mut Self;
    /// The end of a line (of the text format only).
    fn line(&mut self) -> &mut Self;

    fn count(&mut self, v: usize) -> &mut Self {
        self.int(v as i64)
    }

    fn reals(&mut self, values: &[f64]) -> &mut Self {
        for &v in values {
            self.real(v);
        }
        self
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

    fn knots(&mut self, knots: &[f64], multiplicities: &[usize]) -> &mut Self {
        for (k, m) in knots.iter().zip(multiplicities) {
            self.real(*k).count(*m);
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
        self.flag(c.is_rational()).degree(c.degree()).poles(c.poles(), c.weights())
    }

    fn bspline<P: ControlPoint + Coordinates>(&mut self, c: &BSplineCurve<P>) -> &mut Self {
        self.flag(c.is_rational()).flag(c.is_periodic()).degree(c.degree());
        self.count(c.poles().len()).count(c.knots().len());
        self.poles(c.poles(), c.weights()).knots(c.knots(), c.multiplicities())
    }

    fn curve(&mut self, c: &Curve) -> &mut Self {
        match c {
            Curve::Line(l) => self.code(1).point(l.origin).direction(l.direction),
            Curve::Circle(c) => self.code(2).right_handed_frame(c.frame()).real(c.radius()),
            Curve::Ellipse(e) => self.code(3).right_handed_frame(e.frame()).real(e.major_radius()).real(e.minor_radius()),
            Curve::Parabola(p) => self.code(4).right_handed_frame(p.frame()).real(p.focal()),
            Curve::Hyperbola(h) => self.code(5).right_handed_frame(h.frame()).real(h.major_radius()).real(h.minor_radius()),
            Curve::Bezier(b) => self.code(6).bezier(b),
            Curve::BSpline(b) => self.code(7).bspline(b),
            Curve::Trimmed(t) => {
                let (u1, u2) = t.range();
                self.code(8).real(u1).real(u2).line().curve(t.basis())
            }
            Curve::Offset(o) => self.code(9).real(o.offset()).direction(o.direction()).line().curve(o.basis()),
        }
    }

    fn curve2(&mut self, c: &Curve2) -> &mut Self {
        match c {
            Curve2::Line(l) => self.code(1).point(l.origin).direction2(l.direction),
            Curve2::Circle(c) => self.code(2).frame2(c.frame()).real(c.radius()),
            Curve2::Ellipse(e) => self.code(3).frame2(e.frame()).real(e.major_radius()).real(e.minor_radius()),
            Curve2::Parabola(p) => self.code(4).frame2(p.frame()).real(p.focal()),
            Curve2::Hyperbola(h) => self.code(5).frame2(h.frame()).real(h.major_radius()).real(h.minor_radius()),
            Curve2::Bezier(b) => self.code(6).bezier(b),
            Curve2::BSpline(b) => self.code(7).bspline(b),
            Curve2::Trimmed(t) => {
                let (u1, u2) = t.range();
                self.code(8).real(u1).real(u2).line().curve2(t.basis())
            }
            Curve2::Offset(o) => self.code(9).real(o.offset()).line().curve2(o.basis()),
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
            Surface::Plane(p) => self.code(1).frame(p.frame()),
            Surface::Cylinder(c) => self.code(2).frame(c.frame()).real(c.radius()),
            Surface::Cone(c) => self.code(3).frame(c.frame()).real(c.radius()).real(c.semi_angle()),
            Surface::Sphere(s) => self.code(4).frame(s.frame()).real(s.radius()),
            Surface::Torus(t) => self.code(5).frame(t.frame()).real(t.major_radius()).real(t.minor_radius()),
            Surface::Extrusion(e) => self.code(6).direction(e.direction()).line().curve(e.basis()),
            Surface::Revolution(r) => {
                let axis = r.axis();
                self.code(7).point(axis.origin).direction(axis.direction).line().curve(r.basis())
            }
            Surface::Bezier(b) => {
                let rational = b.weights().is_some();
                self.code(8).flag(rational).flag(rational).degree(b.u_degree()).degree(b.v_degree());
                self.grid(b.poles(), b.weights())
            }
            Surface::BSpline(b) => {
                let (u, v) = (b.u_knot_vector(), b.v_knot_vector());
                let rational = b.weights().is_some();
                self.code(9).flag(rational).flag(rational).flag(u.periodic).flag(v.periodic);
                self.degree(u.degree).degree(v.degree).count(b.u_pole_count()).count(b.v_pole_count());
                self.count(u.knots.len()).count(v.knots.len());
                self.grid(b.poles(), b.weights());
                self.knots(&u.knots, &u.multiplicities).knots(&v.knots, &v.multiplicities)
            }
            Surface::Trimmed(t) => {
                let (u1, u2, v1, v2) = t.trim();
                self.code(10).reals(&[u1, u2, v1, v2]).line().surface(t.basis())
            }
            Surface::Offset(o) => self.code(11).real(o.offset()).line().surface(o.basis()),
        }
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

impl Out for Text {
    /// Reals, in exponent form when plain digits would be long: OCCT reads each number
    /// into a short buffer.
    fn real(&mut self, v: f64) -> &mut Self {
        if v == 0.0 || (1e-4..1e15).contains(&v.abs()) {
            self.word(v)
        } else {
            self.word(format_args!("{v:e}"))
        }
    }

    fn int(&mut self, v: i64) -> &mut Self {
        self.word(v)
    }

    fn code(&mut self, v: u8) -> &mut Self {
        self.word(v)
    }

    fn degree(&mut self, v: usize) -> &mut Self {
        self.word(v)
    }

    fn flag(&mut self, b: bool) -> &mut Self {
        self.word(u8::from(b))
    }

    fn line(&mut self) -> &mut Self {
        self.0.push('\n');
        self
    }
}

/// Bytes in OCCT's binary BRep format (`BinTools`).
#[derive(Default)]
struct Bin(Vec<u8>);

impl Bin {
    fn byte(&mut self, b: u8) -> &mut Self {
        self.0.push(b);
        self
    }

    /// A table's header: its name and size as words, then a line's end.
    fn table(&mut self, name: &str, count: usize) -> &mut Self {
        self.0.extend_from_slice(format!("{name} {count}\n").as_bytes());
        self
    }

    fn reference(&mut self, orientation: Orientation, number: usize, location: usize) -> &mut Self {
        let o = match orientation {
            Orientation::Forward => 0,
            Orientation::Reversed => 1,
            Orientation::Internal => 2,
            Orientation::External => 3,
        };
        self.byte(o).count(number).count(location)
    }

    /// A shape's flags (free, modified, checked, orientable, closed, infinite, convex), as
    /// the text format's digits.
    fn flags(&mut self, digits: &str) -> &mut Self {
        for c in digits.bytes() {
            self.byte(c - b'0');
        }
        self
    }
}

impl Out for Bin {
    fn real(&mut self, v: f64) -> &mut Self {
        self.0.extend_from_slice(&v.to_le_bytes());
        self
    }

    fn int(&mut self, v: i64) -> &mut Self {
        self.0.extend_from_slice(&(v as i32).to_le_bytes());
        self
    }

    fn code(&mut self, v: u8) -> &mut Self {
        self.byte(v)
    }

    fn degree(&mut self, v: usize) -> &mut Self {
        self.0.extend_from_slice(&(v as u16).to_le_bytes());
        self
    }

    fn flag(&mut self, b: bool) -> &mut Self {
        self.byte(u8::from(b))
    }

    fn line(&mut self) -> &mut Self {
        self
    }
}

/// A child as written: the node, its orientation and the index of its location.
type Child = (NodeId, Orientation, usize);

/// Gathers the tables of a file, then writes them (as text or binary, by `O`).
struct Writer<'a, O> {
    model: &'a Model,
    version: FormatVersion,
    /// Locations other than the identity; location `i` is `locations[i - 1]`.
    locations: Vec<Transform>,
    curves2: O,
    curve2_count: usize,
    curves: O,
    curve_count: usize,
    surfaces: O,
    surface_count: usize,
    /// Shapes, each after its sub-shapes, with their children.
    shapes: Vec<(NodeId, Vec<Child>)>,
    numbers: HashMap<NodeId, usize>,
    face_surfaces: HashMap<FaceId, usize>,
}

impl<'a, O: Out + Default> Writer<'a, O> {
    fn new(model: &'a Model, version: FormatVersion) -> Self {
        Self {
            model,
            version,
            locations: Vec::new(),
            curves2: O::default(),
            curve2_count: 0,
            curves: O::default(),
            curve_count: 0,
            surfaces: O::default(),
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

    /// The faces' surfaces written, and per edge its uses on faces: (face, the coedge's
    /// orientation on the face, coedge).
    fn surfaces_and_uses(&mut self) -> Result<HashMap<EdgeId, Vec<(FaceId, Orientation, CoEdgeId)>>> {
        let model = self.model;
        let faces: Vec<FaceId> = self.shapes.iter().filter_map(|(n, _)| FaceId::try_from(*n).ok()).collect();
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
        Ok(uses)
    }
}

impl Writer<'_, Text> {
    fn file(mut self, root: Shape) -> Result<String> {
        self.collect(root.node)?;
        let root_location = self.location(root.placement);
        let model = self.model;
        // Surfaces, then the curves of edges and of their uses on faces.
        let mut uses = self.surfaces_and_uses()?;
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

impl Writer<'_, Bin> {
    /// The shape in OCCT's binary BRep format, version 3 (as `BinTools::Write` of OCCT 7.5
    /// writes it; later versions read it): its tables as the text format's, each number in
    /// its own width; no meshes.
    fn file(mut self, root: Shape) -> Result<Vec<u8>> {
        self.collect(root.node)?;
        let root_location = self.location(root.placement);
        let model = self.model;
        let mut uses = self.surfaces_and_uses()?;
        let mut shapes = Bin::default();
        let count = self.shapes.len();
        shapes.table("TShapes", count);
        let shape_list = std::mem::take(&mut self.shapes);
        for (node, children) in &shape_list {
            let shape = Shape::new(*node);
            // (`TopAbs_ShapeEnum`.)
            let kind = match node.kind {
                NodeKind::Compound => 0,
                NodeKind::CompSolid => 1,
                NodeKind::Solid => 2,
                NodeKind::Shell => 3,
                NodeKind::Face => 4,
                NodeKind::Wire => 5,
                NodeKind::Edge => 6,
                NodeKind::Vertex => 7,
                _ => unreachable!("not collected"),
            };
            shapes.byte(kind);
            let flags = match node.kind {
                NodeKind::Vertex => {
                    // (No representations of the point on curves or surfaces.)
                    shapes.real(model.tolerance(shape)).point(model.point(shape)).byte(0);
                    "0101101"
                }
                NodeKind::Edge => {
                    let edge = EdgeId::try_from(*node).map_err(|_| bad("not an edge"))?;
                    let edge_uses = uses.remove(&edge).unwrap_or_default();
                    let same_parameter = edge_uses.iter().all(|u| model.is_same_parameter(u.2));
                    let same_range = edge_uses.iter().all(|u| model.is_same_range(u.2));
                    shapes.real(model.tolerance(shape)).flag(same_parameter).flag(same_range).flag(model.is_degenerate(shape));
                    if let Some((curve, first, last)) = model.edge_curve(shape) {
                        let c = self.add_curve(&curve);
                        shapes.byte(1).count(c).count(0).real(first).real(last);
                    }
                    let mut done = vec![false; edge_uses.len()];
                    for i in 0..edge_uses.len() {
                        if done[i] {
                            continue;
                        }
                        let (face, orientation, coedge) = edge_uses[i];
                        let surface = self.face_surfaces[&face];
                        let seam = (i + 1..edge_uses.len()).find(|&j| !done[j] && edge_uses[j].0 == face && edge_uses[j].1 != orientation);
                        let Some((pcurve, first, last)) = model.effective_pcurve(coedge) else { continue };
                        match seam {
                            Some(j) => {
                                done[j] = true;
                                let Some((other, _, _)) = model.effective_pcurve(edge_uses[j].2) else { continue };
                                let (forward, reversed) = if orientation == Orientation::Reversed { (other, pcurve.clone()) } else { (pcurve.clone(), other) };
                                let (c1, c2) = (self.add_curve2(&forward), self.add_curve2(&reversed));
                                // (Continuity across the seam: `GeomAbs_CN`.)
                                shapes.byte(3).count(c1).count(c2).byte(6).count(surface).count(0).real(first).real(last);
                                shapes.point(forward.point(first)).point(forward.point(last));
                            }
                            None => {
                                let c = self.add_curve2(&pcurve);
                                shapes.byte(2).count(c).count(surface).count(0).real(first).real(last);
                                shapes.point(pcurve.point(first)).point(pcurve.point(last));
                            }
                        }
                    }
                    shapes.byte(0);
                    if model.is_closed_edge(shape) { "0101100" } else { "0101000" }
                }
                NodeKind::Wire => {
                    if model.is_closed_wire(WireId::try_from(*node).map_err(|_| bad("not a wire"))?) { "0101100" } else { "0101000" }
                }
                NodeKind::Face => {
                    let face = FaceId::try_from(*node).map_err(|_| bad("not a face"))?;
                    let surface = self.face_surfaces[&face];
                    // (Natural restriction no; no triangulation, 1 (2 would be one).)
                    shapes.flag(false).real(model.tolerance(shape)).count(surface).count(0).byte(1);
                    "0101000"
                }
                NodeKind::Shell => {
                    if model.is_closed_shell(ShellId::try_from(*node).map_err(|_| bad("not a shell"))?) { "0101100" } else { "0101000" }
                }
                NodeKind::Solid => "0100000",
                _ => "0101000",
            };
            shapes.flags(flags);
            for &(child, orientation, location) in children {
                shapes.reference(orientation, count - self.numbers[&child], location);
            }
            shapes.byte(b'*');
        }
        shapes.reference(root.orientation, count - self.numbers[&root.node], root_location);

        let mut out = Bin::default();
        out.0.extend_from_slice(b"\nOpen CASCADE Topology V3 (c)\n");
        out.table("Locations", self.locations.len());
        for l in &self.locations {
            out.byte(1);
            for row in l.rows() {
                out.reals(&row);
            }
        }
        out.table("Curve2ds", self.curve2_count);
        out.0.extend_from_slice(&self.curves2.0);
        out.table("Curves", self.curve_count);
        out.0.extend_from_slice(&self.curves.0);
        out.table("Polygon3D", 0).table("PolygonOnTriangulations", 0);
        out.table("Surfaces", self.surface_count);
        out.0.extend_from_slice(&self.surfaces.0);
        out.table("Triangulations", 0);
        out.0.extend_from_slice(&shapes.0);
        Ok(out.0)
    }
}

// --- OCCT's shape streams (BinTools_ShapeWriter / ShapeReader) -----------------------------

/// OCCT's shape streams, as its binary documents (OCAF's XBF, format version 12) keep a
/// shape in each attribute that holds one: each shape, location and surface written where
/// first met, a later use of it a reference back to that place in the file (its distance,
/// in 1, 2, 4 or 8 bytes). One writer serves a whole document, its streams written in file
/// order, each at the absolute position given; geometry bodies as the binary BRep format's.
#[derive(Default)]
pub struct ShapeStreamWriter {
    shapes: HashMap<NodeId, u64>,
    locations: Vec<(Transform, u64)>,
    surfaces: HashMap<FaceId, u64>,
    /// Per edge its uses on faces: (face, the coedge's orientation on it, coedge).
    uses: HashMap<EdgeId, Vec<(FaceId, Orientation, CoEdgeId)>>,
    faces: std::collections::HashSet<FaceId>,
}

/// (`BinTools_ObjectType`.)
const STREAM_LOCATION: u8 = 5;
const STREAM_SIMPLE_LOCATION: u8 = 6;
const STREAM_IDENTITY: u8 = 7;
const STREAM_LOCATION_END: u8 = 8;
const STREAM_CURVE: u8 = 9;
const STREAM_NO_CURVE: u8 = 10;
const STREAM_CURVE2: u8 = 11;
const STREAM_NO_CURVE2: u8 = 12;
const STREAM_SURFACE: u8 = 13;
const STREAM_NO_SURFACE: u8 = 14;
const STREAM_POLYGON3: u8 = 15;
const STREAM_POLYGON_ON_MESH: u8 = 17;
const STREAM_MESH: u8 = 19;
const STREAM_NO_SHAPE: u8 = 198;
const STREAM_SHAPE_END: u8 = 199;

/// A shape's flags as a byte (free, modified, checked, orientable, closed, infinite, convex:
/// bits 0 to 6), from the text format's digits.
fn flag_bits(digits: &str) -> u8 {
    digits.bytes().enumerate().fold(0, |b, (i, c)| b | ((c - b'0') << i))
}

fn stream_code(kind: NodeKind) -> u8 {
    match kind {
        NodeKind::Compound => 0,
        NodeKind::CompSolid => 1,
        NodeKind::Solid => 2,
        NodeKind::Shell => 3,
        NodeKind::Face => 4,
        NodeKind::Wire => 5,
        NodeKind::Edge => 6,
        _ => 7,
    }
}

fn orientation_code(o: Orientation) -> u8 {
    match o {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    }
}

impl ShapeStreamWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// The stream of `shape` (its placement and orientation the occurrence's), its first byte
    /// at file position `base`.
    pub fn shape(&mut self, model: &Model, shape: Shape, base: u64) -> Result<Vec<u8>> {
        for face in model.unique_subshapes(shape, NodeKind::Face) {
            let Ok(face) = FaceId::try_from(face.node) else { continue };
            if !self.faces.insert(face) {
                continue;
            }
            let record = model.get(face).ok_or_else(|| bad("no face"))?;
            for &use_id in &record.wires {
                let wire_use = model.get(use_id).ok_or_else(|| bad("no wire use"))?;
                let wire = model.get(wire_use.wire).ok_or_else(|| bad("no wire"))?;
                for &c in &wire.coedges {
                    let coedge = model.get(c).ok_or_else(|| bad("no coedge"))?;
                    self.uses.entry(coedge.edge).or_default().push((face, coedge.orientation.compose(wire_use.orientation), c));
                }
            }
        }
        let mut out = Bin::default();
        self.write_shape(model, &mut out, base, shape.node, shape.orientation, shape.placement)?;
        Ok(out.0)
    }

    /// The stream of a location (`7` for the identity), its first byte at `base`.
    pub fn location(&mut self, t: Transform, base: u64) -> Vec<u8> {
        let mut out = Bin::default();
        self.write_location(&mut out, base, t);
        out.0
    }

    fn reference(out: &mut Bin, base: u64, target: u64) {
        let delta = base + out.0.len() as u64 - target;
        if delta <= 0xff {
            out.byte(1).byte(delta as u8);
        } else if delta <= 0xffff {
            out.byte(2);
            out.0.extend_from_slice(&(delta as u16).to_le_bytes());
        } else if delta <= 0xffff_ffff {
            out.byte(3);
            out.0.extend_from_slice(&(delta as u32).to_le_bytes());
        } else {
            out.byte(4);
            out.0.extend_from_slice(&delta.to_le_bytes());
        }
    }

    fn write_location(&mut self, out: &mut Bin, base: u64, t: Transform) {
        if t == Transform::IDENTITY {
            out.byte(STREAM_IDENTITY);
            return;
        }
        if let Some(&(_, at)) = self.locations.iter().find(|(l, _)| *l == t) {
            Self::reference(out, base, at);
            return;
        }
        self.locations.push((t, base + out.0.len() as u64));
        out.byte(STREAM_SIMPLE_LOCATION);
        for row in t.rows() {
            out.reals(&row);
        }
    }

    fn write_surface(&mut self, model: &Model, out: &mut Bin, base: u64, face: FaceId) -> Result<()> {
        if let Some(&at) = self.surfaces.get(&face) {
            Self::reference(out, base, at);
            return Ok(());
        }
        let surface = model.surface(Shape::new(face.into())).ok_or_else(|| bad("a face without a surface"))?;
        self.surfaces.insert(face, base + out.0.len() as u64);
        out.byte(STREAM_SURFACE).surface(surface);
        Ok(())
    }

    fn write_shape(&mut self, model: &Model, out: &mut Bin, base: u64, node: NodeId, orientation: Orientation, placement: Transform) -> Result<()> {
        if let Some(&at) = self.shapes.get(&node) {
            Self::reference(out, base, at);
            self.write_location(out, base, placement);
            out.byte(orientation_code(orientation));
            return Ok(());
        }
        if matches!(node.kind, NodeKind::CoEdge | NodeKind::Product | NodeKind::Instance) {
            return Err(Error::InvalidInput(format!("shape streams cannot hold a {:?}", node.kind)));
        }
        self.shapes.insert(node, base + out.0.len() as u64);
        out.byte(200 + 4 * stream_code(node.kind) + orientation_code(orientation));
        self.write_location(out, base, placement);
        let shape = Shape::new(node);
        let flags = match node.kind {
            NodeKind::Vertex => {
                out.real(model.tolerance(shape)).point(model.point(shape)).byte(0);
                "0101101"
            }
            NodeKind::Edge => {
                let edge = EdgeId::try_from(node).map_err(|_| bad("not an edge"))?;
                let edge_uses = self.uses.get(&edge).cloned().unwrap_or_default();
                let same_parameter = edge_uses.iter().all(|u| model.is_same_parameter(u.2));
                let same_range = edge_uses.iter().all(|u| model.is_same_range(u.2));
                out.real(model.tolerance(shape)).byte(u8::from(same_parameter) | (u8::from(same_range) << 1) | (u8::from(model.is_degenerate(shape)) << 2));
                if let Some((curve, first, last)) = model.edge_curve(shape) {
                    out.byte(1).byte(STREAM_CURVE).curve(&curve);
                    out.byte(STREAM_IDENTITY).real(first).real(last);
                }
                let mut done = vec![false; edge_uses.len()];
                for i in 0..edge_uses.len() {
                    if done[i] {
                        continue;
                    }
                    let (face, use_orientation, coedge) = edge_uses[i];
                    let seam = (i + 1..edge_uses.len()).find(|&j| !done[j] && edge_uses[j].0 == face && edge_uses[j].1 != use_orientation);
                    let Some((pcurve, first, last)) = model.effective_pcurve(coedge) else { continue };
                    match seam {
                        Some(j) => {
                            done[j] = true;
                            let Some((other, _, _)) = model.effective_pcurve(edge_uses[j].2) else { continue };
                            let (forward, reversed) = if use_orientation == Orientation::Reversed { (other, pcurve) } else { (pcurve, other) };
                            out.byte(3).byte(STREAM_CURVE2).curve2(&forward).byte(STREAM_CURVE2).curve2(&reversed).byte(6);
                        }
                        None => {
                            out.byte(2).byte(STREAM_CURVE2).curve2(&pcurve);
                        }
                    }
                    self.write_surface(model, out, base, face)?;
                    out.byte(STREAM_IDENTITY).real(first).real(last);
                }
                out.byte(0);
                if model.is_closed_edge(shape) { "0101100" } else { "0101000" }
            }
            NodeKind::Face => {
                let face = FaceId::try_from(node).map_err(|_| bad("not a face"))?;
                out.flag(false).real(model.tolerance(shape));
                self.write_surface(model, out, base, face)?;
                // (Identity location; no mesh written.)
                out.byte(STREAM_IDENTITY).byte(0);
                "0101000"
            }
            NodeKind::Wire => {
                if model.is_closed_wire(WireId::try_from(node).map_err(|_| bad("not a wire"))?) { "0101100" } else { "0101000" }
            }
            NodeKind::Shell => {
                if model.is_closed_shell(ShellId::try_from(node).map_err(|_| bad("not a shell"))?) { "0101100" } else { "0101000" }
            }
            NodeKind::Solid => "0100000",
            _ => "0101000",
        };
        out.byte(flag_bits(flags));
        let children: Vec<(NodeId, Orientation, Transform)> = if node.kind == NodeKind::Wire {
            let wire = model.get(WireId::try_from(node).map_err(|_| bad("not a wire"))?).ok_or_else(|| bad("no wire"))?;
            let mut out = Vec::with_capacity(wire.coedges.len());
            for &c in &wire.coedges {
                let coedge = model.get(c).ok_or_else(|| bad("no coedge"))?;
                out.push((coedge.edge.into(), coedge.orientation, Transform::IDENTITY));
            }
            out
        } else {
            model.children_local(shape).map(|c| (c.node, c.orientation, c.placement)).collect()
        };
        for (child, o, t) in children {
            self.write_shape(model, out, base, child, o, t)?;
        }
        out.byte(STREAM_SHAPE_END);
        Ok(())
    }
}

/// Reads OCCT's shape streams (see [`ShapeStreamWriter`]) out of a document's bytes: the
/// shapes whose streams start at `shapes` and the locations at `locations` (positions in
/// file order, every stream of the document among them, as later ones may refer back into
/// earlier ones). A shape met several times is one shape of `model`.
pub fn read_shape_streams(model: &mut Model, data: &[u8], shapes: &[u64], locations: &[u64]) -> Result<(Vec<Option<Shape>>, Vec<Transform>)> {
    let mut reader = StreamReader { file: File { locations: vec![Transform::IDENTITY], curves2: Vec::new(), curves: Vec::new(), surfaces: Vec::new(), shapes: Vec::new(), root: Ref { index: 0, orientation: Orientation::Forward, location: 0 } }, at: HashMap::new() };
    // All the streams parsed in file order, then built.
    let mut starts: Vec<(u64, bool)> = shapes.iter().map(|&p| (p, true)).chain(locations.iter().map(|&p| (p, false))).collect();
    starts.sort_by_key(|s| s.0);
    let mut shape_refs: HashMap<u64, Option<Ref>> = HashMap::new();
    let mut location_refs: HashMap<u64, usize> = HashMap::new();
    for (p, is_shape) in starts {
        let mut b = Bytes { data, at: p as usize };
        if is_shape {
            let r = reader.shape(&mut b)?;
            shape_refs.insert(p, r);
        } else {
            let l = reader.location(&mut b)?;
            location_refs.insert(p, l);
        }
    }
    let file = reader.file;
    let mut builder = Builder::new(model, &file);
    let mut out = Vec::with_capacity(shapes.len());
    for p in shapes {
        out.push(match shape_refs[p] {
            Some(r) => {
                let s = builder.shape(Ref { location: 0, ..r }, &[])?;
                Some(s.with_placement(builder.location(r.location)?))
            }
            None => None,
        });
    }
    let transforms = locations.iter().map(|p| file.locations[location_refs[p]]).collect();
    Ok((out, transforms))
}

/// The tables a stream's definitions fill (a binary BRep file's), by where each was defined.
struct StreamReader {
    file: File,
    /// (kind, position) -> index in its table.
    at: HashMap<(u8, u64), usize>,
}

impl StreamReader {
    /// A reference's target: the position of its type byte less the distance it gives.
    fn target(b: &mut Bytes, code: u8, start: usize) -> Result<u64> {
        let delta = match code {
            1 => b.byte()? as u64,
            2 => u16::from_le_bytes(b.take(2)?.try_into().expect("2 bytes")) as u64,
            3 => u32::from_le_bytes(b.take(4)?.try_into().expect("4 bytes")) as u64,
            _ => u64::from_le_bytes(b.take(8)?.try_into().expect("8 bytes")),
        };
        (start as u64).checked_sub(delta).ok_or_else(|| bad("a reference before the file"))
    }

    fn known(&self, kind: u8, at: u64) -> Result<usize> {
        self.at.get(&(kind, at)).copied().ok_or_else(|| bad(format!("a reference to nothing read at {at}")))
    }

    fn location(&mut self, b: &mut Bytes) -> Result<usize> {
        let start = b.at;
        match b.byte()? {
            STREAM_IDENTITY => Ok(0),
            c @ 1..=4 => {
                let at = Self::target(b, c, start)?;
                self.known(STREAM_LOCATION, at)
            }
            STREAM_SIMPLE_LOCATION => {
                let mut rows = [[0.0; 4]; 3];
                for row in &mut rows {
                    for v in row.iter_mut() {
                        *v = b.real()?;
                    }
                }
                let t = Transform::from_rows(rows).ok_or_else(|| bad("a singular location"))?;
                self.file.locations.push(t);
                let index = self.file.locations.len() - 1;
                self.at.insert((STREAM_LOCATION, start as u64), index);
                Ok(index)
            }
            STREAM_LOCATION => {
                let mut product = Transform::IDENTITY;
                while b.data.get(b.at) != Some(&STREAM_LOCATION_END) {
                    let l = self.location(b)?;
                    let power = b.int()? as i32;
                    product = self.file.locations[l].powi(power) * product;
                }
                b.byte()?;
                self.file.locations.push(product);
                let index = self.file.locations.len() - 1;
                self.at.insert((STREAM_LOCATION, start as u64), index);
                Ok(index)
            }
            c => Err(bad(format!("a location of type {c}"))),
        }
    }

    /// A curve, 2D curve or surface: `None` for an empty one; its index (from 1) otherwise.
    fn geometry(&mut self, b: &mut Bytes, define: u8, empty: u8) -> Result<Option<usize>> {
        let start = b.at;
        let code = b.byte()?;
        if code == empty {
            return Ok(None);
        }
        if (1..=4).contains(&code) {
            let at = Self::target(b, code, start)?;
            return self.known(define, at).map(Some);
        }
        if code != define {
            return Err(bad(format!("geometry of type {code}")));
        }
        let index = match define {
            STREAM_CURVE => {
                self.file.curves.push(RawCurve::read(b)?);
                self.file.curves.len()
            }
            STREAM_CURVE2 => {
                self.file.curves2.push(RawCurve2::read(b)?);
                self.file.curves2.len()
            }
            _ => {
                self.file.surfaces.push(RawSurface::read(b)?);
                self.file.surfaces.len()
            }
        };
        self.at.insert((define, start as u64), index);
        Ok(Some(index))
    }

    /// Meshes and polygons, skipped (their definitions kept for references to them).
    fn mesh(&mut self, b: &mut Bytes, define: u8) -> Result<()> {
        let start = b.at;
        let code = b.byte()?;
        if (1..=4).contains(&code) {
            Self::target(b, code, start)?;
            return Ok(());
        }
        if code != define {
            return Ok(()); // empty
        }
        match define {
            STREAM_POLYGON3 => {
                let (nodes, has_parameters) = (b.count()?, b.flag()?);
                b.skip(8 + 24 * nodes + if has_parameters { 8 * nodes } else { 0 })?;
            }
            STREAM_POLYGON_ON_MESH => {
                let nodes = b.count()?;
                b.skip(4 * nodes + 8)?;
                if b.flag()? {
                    b.skip(8 * nodes)?;
                }
            }
            _ => {
                let (nodes, triangles, has_uv, has_normals) = (b.count()?, b.count()?, b.flag()?, b.flag()?);
                b.skip(8 + 24 * nodes + if has_uv { 16 * nodes } else { 0 } + 12 * triangles + if has_normals { 12 * nodes } else { 0 })?;
            }
        }
        self.at.insert((define, start as u64), 0);
        Ok(())
    }

    fn shape(&mut self, b: &mut Bytes) -> Result<Option<Ref>> {
        let start = b.at;
        let code = b.byte()?;
        if code == STREAM_NO_SHAPE {
            return Ok(None);
        }
        let orientation = |o: u8| match o & 3 {
            0 => Orientation::Forward,
            1 => Orientation::Reversed,
            2 => Orientation::Internal,
            _ => Orientation::External,
        };
        if (1..=4).contains(&code) {
            let at = Self::target(b, code, start)?;
            let index = self.known(STREAM_SHAPE_END, at)?;
            let location = self.location(b)?;
            let o = b.byte()?;
            return Ok(Some(Ref { index, orientation: orientation(o), location }));
        }
        if code < 200 {
            return Err(bad(format!("a shape of type {code}")));
        }
        let kind = match (code - 200) >> 2 {
            0 => NodeKind::Compound,
            1 => NodeKind::CompSolid,
            2 => NodeKind::Solid,
            3 => NodeKind::Shell,
            4 => NodeKind::Face,
            5 => NodeKind::Wire,
            6 => NodeKind::Edge,
            7 => NodeKind::Vertex,
            k => return Err(bad(format!("unknown shape type {k}"))),
        };
        let location = self.location(b)?;
        let geometry = match kind {
            NodeKind::Vertex => {
                let (tolerance, point) = (b.real()?, b.point()?);
                loop {
                    match b.byte()? {
                        0 => break,
                        1 => {
                            b.real()?;
                            self.geometry(b, STREAM_CURVE, STREAM_NO_CURVE)?;
                        }
                        2 => {
                            b.real()?;
                            self.geometry(b, STREAM_CURVE2, STREAM_NO_CURVE2)?;
                            self.geometry(b, STREAM_SURFACE, STREAM_NO_SURFACE)?;
                        }
                        3 => {
                            b.real()?;
                            b.real()?;
                            self.geometry(b, STREAM_SURFACE, STREAM_NO_SURFACE)?;
                        }
                        t => return Err(bad(format!("unknown point representation {t}"))),
                    }
                    self.location(b)?;
                }
                Geometry::Vertex { point, tolerance }
            }
            NodeKind::Edge => {
                let tolerance = b.real()?;
                let flags = b.byte()?;
                let (mut curve, mut on_surfaces) = (None, Vec::new());
                loop {
                    match b.byte()? {
                        0 => break,
                        1 => {
                            let c = self.geometry(b, STREAM_CURVE, STREAM_NO_CURVE)?;
                            let l = self.location(b)?;
                            let (first, last) = (b.real()?, b.real()?);
                            if let Some(c) = c {
                                curve = Some((c, l, first, last));
                            }
                        }
                        rep @ (2 | 3) => {
                            let c = self.geometry(b, STREAM_CURVE2, STREAM_NO_CURVE2)?;
                            let reversed = if rep == 3 {
                                let r = self.geometry(b, STREAM_CURVE2, STREAM_NO_CURVE2)?;
                                b.byte()?; // continuity
                                r
                            } else {
                                None
                            };
                            let s = self.geometry(b, STREAM_SURFACE, STREAM_NO_SURFACE)?;
                            let l = self.location(b)?;
                            let (first, last) = (b.real()?, b.real()?);
                            if let (Some(curve), Some(surface)) = (c, s) {
                                on_surfaces.push(OnSurface { curve, reversed_curve: reversed, surface, location: l, first, last });
                            }
                        }
                        4 => {
                            b.byte()?;
                            for _ in 0..2 {
                                self.geometry(b, STREAM_SURFACE, STREAM_NO_SURFACE)?;
                                self.location(b)?;
                            }
                        }
                        5 => {
                            self.mesh(b, STREAM_POLYGON3)?;
                            self.location(b)?;
                        }
                        rep @ (6 | 7) => {
                            self.mesh(b, STREAM_POLYGON_ON_MESH)?;
                            if rep == 7 {
                                self.mesh(b, STREAM_POLYGON_ON_MESH)?;
                            }
                            self.mesh(b, STREAM_MESH)?;
                            self.location(b)?;
                        }
                        t => return Err(bad(format!("unknown edge representation {t}"))),
                    }
                }
                Geometry::Edge { tolerance, degenerate: flags & 4 != 0, curve, on_surfaces }
            }
            NodeKind::Face => {
                b.byte()?; // natural restriction
                let tolerance = b.real()?;
                let surface = self.geometry(b, STREAM_SURFACE, STREAM_NO_SURFACE)?;
                let l = self.location(b)?;
                if b.byte()? == 2 {
                    self.mesh(b, STREAM_MESH)?;
                }
                Geometry::Face { tolerance, surface: surface.map(|s| (s, l)) }
            }
            _ => Geometry::None,
        };
        b.byte()?; // flags
        // (Its place in the table taken before its sub-shapes', as they may refer to it.)
        let index = self.file.shapes.len();
        self.file.shapes.push(RawShape { kind, geometry, children: Vec::new() });
        self.at.insert((STREAM_SHAPE_END, start as u64), index);
        let mut children = Vec::new();
        while b.data.get(b.at) != Some(&STREAM_SHAPE_END) {
            if let Some(r) = self.shape(b)? {
                children.push(r);
            }
        }
        b.byte()?;
        self.file.shapes[index].children = children;
        Ok(Some(Ref { index, orientation: orientation(code - 200), location }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT's binary BRep (BinTools), versions 3 and 4 (CadQuery 2.8's OCCT 7.9 writes 4):
    /// a box of 10 with a hole of radius 2 and its vertical edges rounded by 1, turned by 30°
    /// about z and moved to (1, 2, 3) by locations, and meshed (tables read past). Its
    /// volume, area and centre as OCCT's; its vertices farthest along x where the rounding
    /// leaves the side at (5, -4), turned: at 1 + 5 cos 30° + 4 sin 30°.
    #[test]
    fn reads_binary_files() {
        for data in [&include_bytes!("tests/filleted_box_v3.bin")[..], &include_bytes!("tests/filleted_box_v4.bin")[..]] {
            let mut model = Model::new();
            let solid = read_binary(&mut model, &mut &data[..]).unwrap();
            assert!(crate::validity::check(&model, solid).is_valid());
            let volume = crate::properties::volume(&model, solid);
            let area = crate::properties::area(&model, solid).mass;
            assert!((volume.mass - 865.752220392).abs() < 1e-6, "{}", volume.mass);
            assert!((area - 681.646003294).abs() < 1e-6, "{area}");
            assert!(volume.center_of_mass.distance(Point::new(1.0, 2.0, 3.0)) < 1e-6, "{:?}", volume.center_of_mass);
            let x = model.unique_subshapes(solid, NodeKind::Vertex).into_iter().map(|v| model.point(v).x).fold(f64::MIN, f64::max);
            let (sin, cos) = 30f64.to_radians().sin_cos();
            assert!((x - (1.0 + 5.0 * cos + 4.0 * sin)).abs() < 1e-9, "{x}");
        }
    }

    /// A box 1 x 2 x 3 OCCT wrote unmeshed (version 4): each face's surface followed by 1, no
    /// triangulation.
    #[test]
    fn reads_binary_files_without_meshes() {
        let mut model = Model::new();
        let solid = read_binary(&mut model, &mut &include_bytes!("tests/box_unmeshed_v4.bin")[..]).unwrap();
        assert!(crate::validity::check(&model, solid).is_valid());
        assert!((crate::properties::volume(&model, solid).mass - 6.0).abs() < 1e-12);
    }

    /// The filleted box of `reads_binary_files`, written in the binary format and read back:
    /// the same volume, area and centre (its fillets' seams on two pcurves each).
    #[test]
    fn writes_binary_files_that_read_back() {
        let mut model = Model::new();
        let solid = read_binary(&mut model, &mut &include_bytes!("tests/filleted_box_v4.bin")[..]).unwrap();
        let mut bytes = Vec::new();
        write_binary(&model, solid, &mut bytes).unwrap();
        let mut other = Model::new();
        let back = read_binary(&mut other, &mut bytes.as_slice()).unwrap();
        assert!(crate::validity::check(&other, back).is_valid());
        let volume = crate::properties::volume(&other, back);
        assert!((volume.mass - 865.752220392).abs() < 1e-6, "{}", volume.mass);
        assert!((crate::properties::area(&other, back).mass - 681.646003294).abs() < 1e-6);
        assert!(volume.center_of_mass.distance(Point::new(1.0, 2.0, 3.0)) < 1e-6, "{:?}", volume.center_of_mass);
    }

    /// The filleted box of `reads_binary_files` in a shape stream, then the same box moved in
    /// a second stream after it (written as a reference back into the first, with its
    /// placement), and a placement stream: read back, the same box twice, sharing its faces,
    /// and the placement.
    #[test]
    fn shape_streams_read_back() {
        let mut model = Model::new();
        let solid = read_binary(&mut model, &mut &include_bytes!("tests/filleted_box_v4.bin")[..]).unwrap();
        let moved = solid.moved(&Transform::translation(Vector::new(10.0, 0.0, 0.0)));
        let mut writer = ShapeStreamWriter::new();
        let mut data = vec![0u8; 5];
        let first = data.len() as u64;
        data.extend(writer.shape(&model, solid, first).unwrap());
        let second = data.len() as u64;
        let again = writer.shape(&model, moved, second).unwrap();
        // (A reference, its placement (a type byte and 12 reals) and orientation.)
        assert!(again.len() <= 1 + 8 + 1 + 96 + 1, "the second stream refers back: {} bytes", again.len());
        data.extend(again);
        let third = data.len() as u64;
        data.extend(writer.location(Transform::translation(Vector::new(1.0, 2.0, 3.0)), third));
        let mut other = Model::new();
        let (shapes, locations) = read_shape_streams(&mut other, &data, &[first, second], &[third]).unwrap();
        let (a, b) = (shapes[0].unwrap(), shapes[1].unwrap());
        assert!(crate::validity::check(&other, a).is_valid());
        for (s, x) in [(a, 1.0), (b, 11.0)] {
            let volume = crate::properties::volume(&other, s);
            assert!((volume.mass - 865.752220392).abs() < 1e-6, "{}", volume.mass);
            assert!((volume.center_of_mass.x - x).abs() < 1e-6, "{:?}", volume.center_of_mass);
        }
        assert_eq!(a.node, b.node);
        assert!(locations[0] == Transform::translation(Vector::new(1.0, 2.0, 3.0)));
    }

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
