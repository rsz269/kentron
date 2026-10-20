//! Reading STEP files (ISO 10303, application protocols 203, 214 and 242): boundary
//! representations of solids and shells, their assembly structure (products and placed
//! occurrences of them) and the colours of their solids.
//!
//! Geometry keeps the file's parametrization, so that the pcurves it gives fit the edges'
//! curves; where a pcurve is missing or does not fit, the edge's curve is projected onto the
//! face's surface instead. Faces' boundaries are turned, where needed, to run counter-clockwise
//! in the surface's parameters (outer boundaries; holes clockwise), the face itself reversed
//! where the file says it faces against its surface.

mod parse;
mod write;

pub use write::{write, write_assembly, write_assembly_in, write_in};

use crate::hashing::FastMap as HashMap;
use std::io::Read;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Circle, Curve, Ellipse, Line};
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Line2};
use crate::geometry::extrema::nearest_on_curve;
use crate::geometry::surface::{
    BSplineSurface, Cone, Cylinder, ExtrusionSurface, KnotVector, OffsetSurface, Plane, RevolutionSurface, Sphere, Surface, Torus,
};
use crate::geometry::ParametricCurve;
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Vector};
use crate::topology::{EdgeId, Model, NodeKind, Shape, VertexId};

use super::faces::{align_periods, finish_face, fits, follows, project, CurveSpan, Gaps};

pub mod binstep;
pub use parse::{parse, Exchange, Instance, Value};

/// A node of a file's assembly structure: a product (part or assembly) as placed in its
/// parent.
#[derive(Debug, Clone, Default)]
pub struct Node {
    /// The occurrence's name (or the product's, for a root or an unnamed occurrence).
    pub name: String,
    /// The product's own shape (a solid, or a compound of its solids and shells); `None`
    /// for an assembly of other products only.
    pub shape: Option<Shape>,
    /// The colour of all its solids, if they share one.
    pub color: Option<[f64; 3]>,
    /// That colour's opacity, where it is not opaque (1 opaque, 0 clear).
    pub alpha: Option<f64>,
    /// Where the product sits in its parent.
    pub placement: Transform,
    pub children: Vec<Node>,
    /// Faces and edges of the shape with names, colours or layers of their own (a wire's
    /// layer is its edges').
    pub subshapes: Vec<Subshape>,
    /// What the product is made of.
    pub material: Option<Material>,
}

/// A face or edge of a product's shape with a name, colour or layer of its own.
#[derive(Debug, Clone, PartialEq)]
pub struct Subshape {
    /// A face, edge or wire of the product's shape (unplaced within it).
    pub shape: Shape,
    pub name: Option<String>,
    /// Red, green, blue and opacity.
    pub color: Option<[f64; 4]>,
    pub layer: Option<String>,
}

/// A product's material: its name and description, and its density (in its unit, a name
/// such as "g/cm^3").
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub name: String,
    pub description: String,
    pub density: f64,
    pub density_unit: String,
}

/// A length unit of a STEP file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit {
    Millimetre,
    Centimetre,
    Metre,
    Kilometre,
    Micrometre,
    Nanometre,
    Inch,
    Foot,
    Mile,
}

impl LengthUnit {
    /// Millimetres per unit.
    pub fn millimetres(self) -> f64 {
        match self {
            LengthUnit::Millimetre => 1.0,
            LengthUnit::Centimetre => 10.0,
            LengthUnit::Metre => 1e3,
            LengthUnit::Kilometre => 1e6,
            LengthUnit::Micrometre => 1e-3,
            LengthUnit::Nanometre => 1e-6,
            LengthUnit::Inch => 25.4,
            LengthUnit::Foot => 304.8,
            LengthUnit::Mile => 1_609_344.0,
        }
    }

    /// The unit's name in a STEP conversion-based unit (for those not SI).
    fn name(self) -> &'static str {
        match self {
            LengthUnit::Inch => "INCH",
            LengthUnit::Foot => "FOOT",
            LengthUnit::Mile => "MILE",
            _ => "",
        }
    }

    /// The SI prefix of a metric unit as written in `SI_UNIT` (`$` for the metre itself);
    /// `None` for units defined by conversion.
    fn si_prefix(self) -> Option<&'static str> {
        match self {
            LengthUnit::Millimetre => Some(".MILLI."),
            LengthUnit::Centimetre => Some(".CENTI."),
            LengthUnit::Metre => Some("$"),
            LengthUnit::Kilometre => Some(".KILO."),
            LengthUnit::Micrometre => Some(".MICRO."),
            LengthUnit::Nanometre => Some(".NANO."),
            _ => None,
        }
    }
}

/// Reads a STEP file's products as a tree: one node per root product. The file may be STEP
/// text or `.binstep` (told apart by its first bytes).
#[doc(alias = "STEPCAFControl_Reader")]
pub fn read_assembly(model: &mut Model, input: &mut dyn Read) -> Result<Vec<Node>> {
    let mut bytes = Vec::new();
    input.read_to_end(&mut bytes).map_err(|e| Error::InvalidInput(format!("reading STEP: {e}")))?;
    let exchange = exchange_of(&bytes)?;
    let mut reader = Reader::new(&exchange, model)?;
    let roots = reader.assembly()?;
    reader.gaps.apply(reader.model)?;
    Ok(roots)
}

/// Reads a STEP file's shapes as one: a compound of its products' shapes, placed.
#[doc(alias = "STEPControl_Reader")]
pub fn read(model: &mut Model, input: &mut dyn Read) -> Result<Shape> {
    let roots = read_assembly(model, input)?;
    let mut shapes = Vec::new();
    fn collect(node: &Node, placement: Transform, out: &mut Vec<Shape>) {
        let here = placement * node.placement;
        if let Some(s) = node.shape {
            out.push(s.moved(&here));
        }
        for child in &node.children {
            collect(child, here, out);
        }
    }
    for root in &roots {
        collect(root, Transform::IDENTITY, &mut shapes);
    }
    if shapes.len() == 1 {
        return Ok(shapes[0]);
    }
    let compound = model.new_empty(NodeKind::Compound);
    for s in shapes {
        model.add(compound, s)?;
    }
    Ok(compound)
}

/// The exchange structure of a STEP file's bytes: `.binstep` decoded, text parsed.
pub fn exchange_of(bytes: &[u8]) -> Result<Exchange> {
    if binstep::is_binstep(bytes) { binstep::decode(bytes) } else { parse(bytes) }
}

/// A STEP file's bytes (text or `.binstep`) in the other form or the same: `.binstep`
/// if `binary`, else text.
pub fn convert(bytes: &[u8], binary: bool) -> Result<Vec<u8>> {
    let exchange = exchange_of(bytes)?;
    Ok(if binary { binstep::encode(&exchange) } else { exchange.to_text().into_bytes() })
}

/// [`write`], as `.binstep`.
pub fn write_binstep(model: &Model, shape: Shape, name: &str, color: Option<[f64; 3]>, out: &mut dyn std::io::Write) -> Result<()> {
    let mut text = Vec::new();
    write(model, shape, name, color, &mut text)?;
    out.write_all(&convert(&text, true)?).map_err(|e| Error::InvalidInput(format!("writing .binstep: {e}")))
}

/// [`write_assembly`], as `.binstep`.
pub fn write_assembly_binstep(model: &Model, roots: &[Node], out: &mut dyn std::io::Write) -> Result<()> {
    let mut text = Vec::new();
    write_assembly(model, roots, &mut text)?;
    out.write_all(&convert(&text, true)?).map_err(|e| Error::InvalidInput(format!("writing .binstep: {e}")))
}

fn invalid(what: String) -> Error {
    Error::InvalidInput(format!("STEP: {what}"))
}

/// An edge of the file as made in the model: the model's edge runs along its curve's
/// parameter; `flipped` if that is against the file's edge (from its second vertex).
impl EdgeMade {
    fn span(&self) -> CurveSpan<'_> {
        CurveSpan { curve: &self.curve, first: self.first, last: self.last }
    }
}

#[derive(Clone)]
struct EdgeMade {
    id: EdgeId,
    flipped: bool,
    curve: Curve,
    first: f64,
    last: f64,
    /// For a line, the file's parameter per unit of the model's (its vector's length).
    scale: f64,
    /// Pcurves given by the file: (surface instance, curve instance).
    pcurves: Vec<(usize, usize)>,
}

struct Reader<'a> {
    x: &'a Exchange,
    model: &'a mut Model,
    /// Millimetres per length unit of the file.
    length: f64,
    /// Radians per angle unit of the file.
    angle: f64,
    vertices: HashMap<usize, VertexId>,
    /// Vertices made at the ends of boundary curves (which the file gives no vertices),
    /// shared by position.
    point_vertices: Vec<(Point, VertexId)>,
    edges: HashMap<usize, EdgeMade>,
    surfaces: HashMap<usize, Surface>,
    /// Solids and shells by their instance, made once.
    bodies: HashMap<usize, Shape>,
    colors: HashMap<usize, [f64; 4]>,
    layers: HashMap<usize, String>,
    materials: HashMap<usize, Material>,
    /// The file's instance of each face and edge read (for their names, colours and layers).
    face_instances: HashMap<crate::topology::NodeId, usize>,
    edge_instances: HashMap<crate::topology::NodeId, usize>,
    /// Tolerances needed to cover the gaps between the file's curves, pcurves and points.
    gaps: Gaps,
}

impl<'a> Reader<'a> {
    fn new(x: &'a Exchange, model: &'a mut Model) -> Result<Self> {
        let (length, angle) = units(x);
        let mut reader = Reader {
            x,
            model,
            length,
            angle,
            vertices: HashMap::default(),
            point_vertices: Vec::new(),
            edges: HashMap::default(),
            surfaces: HashMap::default(),
            bodies: HashMap::default(),
            colors: HashMap::default(),
            layers: HashMap::default(),
            materials: HashMap::default(),
            face_instances: HashMap::default(),
            edge_instances: HashMap::default(),
            gaps: Gaps::default(),
        };
        reader.colors = colors(x);
        reader.layers = layers(x);
        reader.materials = materials(x);
        Ok(reader)
    }

    fn get(&self, id: usize) -> Result<Instance<'a>> {
        self.x.get(id)
    }

    fn reference(v: &Value) -> Result<usize> {
        v.as_reference().ok_or_else(|| invalid(format!("expected a reference, got {v:?}")))
    }

    /// A surface curve's (or seam or intersection curve's) 3D curve and its list of pcurves:
    /// a complex instance's record has no name before them.
    fn surface_curve(&self, inst: &Instance<'a>) -> Option<(Value, Value)> {
        let p = inst.record("SURFACE_CURVE").or(inst.record("SEAM_CURVE")).or(inst.record("INTERSECTION_CURVE"))?;
        let k = if inst.is_complex() { 0 } else { 1 };
        Some((p.get(k)?.clone(), p.get(k + 1)?.clone()))
    }

    fn real(&self, v: &Value) -> Result<f64> {
        self.x.real(v).ok_or_else(|| invalid(format!("expected a number, got {v:?}")))
    }

    fn list(&self, v: &Value) -> Result<&'a [Value]> {
        let x: &'a Exchange = self.x;
        x.list(v).ok_or_else(|| invalid(format!("expected a list, got {v:?}")))
    }

    // --- placements and points -----------------------------------------------------------

    fn coordinates(&self, id: usize) -> Result<Vec<f64>> {
        let inst = self.get(id)?;
        let p = inst.record("CARTESIAN_POINT").ok_or_else(|| invalid(format!("#{id} is not a point")))?;
        self.list(&p[1])?.iter().map(|v| self.real(v)).collect()
    }

    fn point(&self, id: usize) -> Result<Point> {
        let c = self.coordinates(id)?;
        let s = self.length;
        Ok(Point::new(c[0] * s, c.get(1).copied().unwrap_or(0.0) * s, c.get(2).copied().unwrap_or(0.0) * s))
    }

    fn point2(&self, id: usize) -> Result<Point2> {
        // Parameter space: no unit.
        let c = self.coordinates(id)?;
        Ok(Point2::new(c[0], c[1]))
    }

    fn ratios(&self, id: usize) -> Result<Vec<f64>> {
        let inst = self.get(id)?;
        let d = inst.record("DIRECTION").ok_or_else(|| invalid(format!("#{id} is not a direction")))?;
        self.list(&d[1])?.iter().map(|v| self.real(v)).collect()
    }

    fn direction(&self, id: usize) -> Result<Direction> {
        let r = self.ratios(id)?;
        Direction::try_from_vector(Vector::new(r[0], r[1], r.get(2).copied().unwrap_or(0.0))).ok_or_else(|| invalid(format!("#{id}: a zero direction")))
    }

    fn direction2(&self, id: usize) -> Result<Direction2> {
        let r = self.ratios(id)?;
        let v = crate::math::Vector2::new(r[0], r[1]);
        if v.length() == 0.0 {
            return Err(invalid(format!("#{id}: a zero direction")));
        }
        Ok(Direction2::from_vector(v))
    }

    /// `VECTOR(name, direction, magnitude)`: the unit direction and the length.
    fn vector(&self, id: usize) -> Result<(usize, f64)> {
        let p = self.get(id)?.record("VECTOR").ok_or_else(|| invalid(format!("#{id} is not a vector")))?;
        Ok((Self::reference(&p[1])?, self.real(&p[2])?))
    }

    fn frame(&self, id: usize) -> Result<Frame> {
        let inst = self.get(id)?;
        let p = inst.record("AXIS2_PLACEMENT_3D").ok_or_else(|| invalid(format!("#{id} is not a placement")))?;
        let origin = self.point(Self::reference(&p[1])?)?;
        let z = match p.get(2).and_then(Value::as_reference) {
            Some(d) => self.direction(d)?,
            None => Direction::Z,
        };
        let x = match p.get(3).and_then(Value::as_reference) {
            Some(d) => self.direction(d)?,
            None => {
                if z.z().abs() < 0.9 {
                    Direction::X
                } else {
                    Direction::from_vector(Vector::from(z).cross(Vector::new(1.0, 0.0, 0.0)).cross(Vector::from(z)))
                }
            }
        };
        // (A reference direction along the axis says nothing of X: as OCCT, any one then.)
        if Vector::from(z).cross(Vector::from(x)).length() <= 1e-9 {
            return Ok(Frame::from_z(origin, z));
        }
        Ok(Frame::new(origin, z, x))
    }

    fn frame2(&self, id: usize) -> Result<Frame2> {
        let p = self.get(id)?.record("AXIS2_PLACEMENT_2D").ok_or_else(|| invalid(format!("#{id} is not a 2D placement")))?;
        let origin = self.point2(Self::reference(&p[1])?)?;
        let x = match p.get(2).and_then(Value::as_reference) {
            Some(d) => self.direction2(d)?,
            None => Direction2::X,
        };
        Ok(Frame2::new(origin, x, true))
    }

    // --- curves --------------------------------------------------------------------------

    /// A 3D curve, and for a line the file's parameter per unit of length.
    fn curve(&self, id: usize) -> Result<(Curve, f64)> {
        let inst = self.get(id)?;
        let s = self.length;
        if let Some(p) = inst.record("LINE") {
            let origin = self.point(Self::reference(&p[1])?)?;
            let (d, magnitude) = self.vector(Self::reference(&p[2])?)?;
            return Ok((Curve::Line(Line { origin, direction: self.direction(d)? }), magnitude * s));
        }
        if let Some(p) = inst.record("CIRCLE") {
            return Ok((Curve::Circle(Circle::new(self.frame(Self::reference(&p[1])?)?, self.real(&p[2])? * s)), 1.0));
        }
        if let Some(p) = inst.record("ELLIPSE") {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok((Curve::Ellipse(Ellipse::new(frame, self.real(&p[2])? * s, self.real(&p[3])? * s)), 1.0));
        }
        if inst.has("B_SPLINE_CURVE_WITH_KNOTS") || inst.has("B_SPLINE_CURVE") {
            return Ok((Curve::BSpline(self.bspline(&inst, false)?.0), 1.0));
        }
        if let Some((curve, _)) = self.surface_curve(&inst) {
            return self.curve(Self::reference(&curve)?);
        }
        if let Some(p) = inst.record("TRIMMED_CURVE") {
            // The basis (edges set their own range from their vertices).
            return self.curve(Self::reference(&p[1])?);
        }
        Err(invalid(format!("#{id}: curves of type {} are not read yet", inst.name())))
    }

    /// The poles, weights, knots, multiplicities and degree of a B-spline curve (3D, or 2D
    /// if `planar`).
    fn bspline_data(&self, inst: &Instance, planar: bool) -> Result<(Vec<Vec<f64>>, Option<Vec<f64>>, Vec<f64>, Vec<usize>, usize)> {
        // Simple instance: (name, degree, poles, form, closed, self_intersect, mults, knots,
        // spec); complex: B_SPLINE_CURVE(degree, poles, ...) and B_SPLINE_CURVE_WITH_KNOTS(
        // mults, knots, spec).
        let (degree, poles, knots_record) = if !inst.is_complex() && inst.has("B_SPLINE_CURVE_WITH_KNOTS") {
            let p = inst.params();
            (&p[1], &p[2], Some(&p[6..]))
        } else {
            let b = inst.record("B_SPLINE_CURVE").ok_or_else(|| invalid("a B-spline curve without its record".into()))?;
            (&b[0], &b[1], inst.record("B_SPLINE_CURVE_WITH_KNOTS"))
        };
        let degree = degree.as_integer().ok_or_else(|| invalid("a B-spline degree".into()))? as usize;
        let poles: Vec<Vec<f64>> = self.list(poles)?
            .iter()
            .map(|v| {
                let c = self.coordinates(Self::reference(v)?)?;
                Ok(if planar { c } else { c.iter().map(|x| x * self.length).collect() })
            })
            .collect::<Result<_>>()?;
        let (mults, knots) = match knots_record {
            Some(k) => (
                self.list(&k[0])?.iter().map(|v| v.as_integer().map(|i| i as usize).ok_or_else(|| invalid("a multiplicity".into()))).collect::<Result<Vec<_>>>()?,
                self.list(&k[1])?.iter().map(|v| self.real(v)).collect::<Result<Vec<_>>>()?,
            ),
            None => uniform_knots(inst, poles.len(), degree),
        };
        let weights = match inst.record("RATIONAL_B_SPLINE_CURVE") {
            Some(w) => Some(self.list(&w[0])?.iter().map(|v| self.real(v)).collect::<Result<Vec<_>>>()?),
            None => None,
        };
        Ok((poles, weights, knots, mults, degree))
    }

    fn bspline(&self, inst: &Instance, _planar: bool) -> Result<(BSplineCurve, ())> {
        let (poles, weights, knots, mults, degree) = self.bspline_data(inst, false)?;
        let poles: Vec<Point> = poles.into_iter().map(|c| Point::new(c[0], c[1], c.get(2).copied().unwrap_or(0.0))).collect();
        let curve = BSplineCurve::try_new(poles, knots, mults, degree)?;
        Ok((match weights {
            Some(w) => curve.with_weights(w),
            None => curve,
        }, ()))
    }

    /// A curve in a surface's parameters, of the edge whose model parameter is the file's
    /// times `scale`, over `[first, last]`.
    fn curve2(&self, id: usize, scale: f64, first: f64, last: f64) -> Result<Curve2> {
        let inst = self.get(id)?;
        if let Some(p) = inst.record("PCURVE") {
            let rep = self.get(Self::reference(&p[2])?)?;
            let items = rep.record("DEFINITIONAL_REPRESENTATION").ok_or_else(|| invalid("a pcurve without its representation".into()))?;
            let curve = Self::reference(&self.list(&items[1])?[0])?;
            return self.curve2(curve, scale, first, last);
        }
        if let Some(p) = inst.record("LINE") {
            let origin = self.point2(Self::reference(&p[1])?)?;
            let (d, magnitude) = self.vector(Self::reference(&p[2])?)?;
            let direction = self.direction2(d)?;
            // The file's parameter t = s / scale moves magnitude · t along the line.
            let speed = magnitude / scale;
            if (speed - 1.0).abs() <= 1e-12 {
                return Ok(Curve2::Line(Line2::new(origin, direction)));
            }
            let at = |s: f64| origin + crate::math::Vector2::from(direction) * (s * speed);
            return Ok(Curve2::BSpline(BSplineCurve2::new(vec![at(first), at(last)], vec![first, last], vec![2, 2], 1)));
        }
        if (scale - 1.0).abs() > 1e-12 {
            if inst.has("B_SPLINE_CURVE_WITH_KNOTS") || inst.has("B_SPLINE_CURVE") {
                let (poles, weights, knots, mults, degree) = self.bspline_data(&inst, true)?;
                let poles = poles.into_iter().map(|c| Point2::new(c[0], c[1])).collect();
                let knots = knots.iter().map(|k| k * scale).collect();
                let c = BSplineCurve2::try_new(poles, knots, mults, degree)?;
                return Ok(Curve2::BSpline(match weights {
                    Some(w) => c.with_weights(w),
                    None => c,
                }));
            }
            return Err(invalid(format!("#{id}: a pcurve of a rescaled line edge")));
        }
        if let Some(p) = inst.record("CIRCLE") {
            return Ok(Curve2::Circle(Circle2::new(self.frame2(Self::reference(&p[1])?)?, self.real(&p[2])?)));
        }
        if let Some(p) = inst.record("ELLIPSE") {
            return Ok(Curve2::Ellipse(Ellipse2::new(self.frame2(Self::reference(&p[1])?)?, self.real(&p[2])?, self.real(&p[3])?)));
        }
        if inst.has("B_SPLINE_CURVE_WITH_KNOTS") || inst.has("B_SPLINE_CURVE") {
            let (poles, weights, knots, mults, degree) = self.bspline_data(&inst, true)?;
            let poles = poles.into_iter().map(|c| Point2::new(c[0], c[1])).collect();
            let c = BSplineCurve2::try_new(poles, knots, mults, degree)?;
            return Ok(Curve2::BSpline(match weights {
                Some(w) => c.with_weights(w),
                None => c,
            }));
        }
        Err(invalid(format!("#{id}: 2D curves of type {} are not read yet", inst.name())))
    }

    // --- surfaces ------------------------------------------------------------------------

    fn surface(&mut self, id: usize) -> Result<Surface> {
        if let Some(s) = self.surfaces.get(&id) {
            return Ok(s.clone());
        }
        let s = self.make_surface(id)?;
        self.surfaces.insert(id, s.clone());
        Ok(s)
    }

    fn make_surface(&mut self, id: usize) -> Result<Surface> {
        let inst = self.get(id)?;
        let l = self.length;
        if let Some(p) = inst.record("PLANE") {
            return Ok(Surface::Plane(Plane::from_frame(self.frame(Self::reference(&p[1])?)?)));
        }
        if let Some(p) = inst.record("CYLINDRICAL_SURFACE") {
            return Ok(Surface::Cylinder(Cylinder::new(self.frame(Self::reference(&p[1])?)?, self.real(&p[2])? * l)));
        }
        if let Some(p) = inst.record("CONICAL_SURFACE") {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok(Surface::Cone(Cone::new(frame, self.real(&p[3])? * self.angle, self.real(&p[2])? * l)));
        }
        if let Some(p) = inst.record("SPHERICAL_SURFACE") {
            return Ok(Surface::Sphere(Sphere::new(self.frame(Self::reference(&p[1])?)?, self.real(&p[2])? * l)));
        }
        if let Some(p) = inst.record("TOROIDAL_SURFACE") {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok(Surface::Torus(Torus::new(frame, self.real(&p[2])? * l, self.real(&p[3])? * l)));
        }
        if let Some(p) = inst.record("SURFACE_OF_LINEAR_EXTRUSION") {
            let (curve, _) = self.curve(Self::reference(&p[1])?)?;
            let (d, _) = self.vector(Self::reference(&p[2])?)?;
            return Ok(Surface::Extrusion(ExtrusionSurface::new(curve, self.direction(d)?)));
        }
        if let Some(p) = inst.record("SURFACE_OF_REVOLUTION") {
            let (curve, _) = self.curve(Self::reference(&p[1])?)?;
            let a = self.get(Self::reference(&p[2])?)?.record("AXIS1_PLACEMENT").ok_or_else(|| invalid("a revolution axis".into()))?.to_vec();
            let origin = self.point(Self::reference(&a[1])?)?;
            let dir = match a[2].as_reference() {
                Some(d) => self.direction(d)?,
                None => Direction::Z,
            };
            return Ok(Surface::Revolution(RevolutionSurface::new(curve, Axis::new(origin, dir))));
        }
        if let Some(p) = inst.record("OFFSET_SURFACE") {
            let basis = self.surface(Self::reference(&p[1])?)?;
            return Ok(Surface::Offset(OffsetSurface::new(basis, self.real(&p[2])? * l)));
        }
        if inst.has("B_SPLINE_SURFACE_WITH_KNOTS") || inst.has("B_SPLINE_SURFACE") {
            return self.bspline_surface(&inst);
        }
        Err(invalid(format!("#{id}: surfaces of type {} are not read yet", inst.name())))
    }

    fn bspline_surface(&self, inst: &Instance) -> Result<Surface> {
        let (ud, vd, grid, knots) = if !inst.is_complex() {
            let p = inst.params();
            (&p[1], &p[2], &p[3], &p[8..])
        } else {
            let b = inst.record("B_SPLINE_SURFACE").ok_or_else(|| invalid("a B-spline surface without its record".into()))?;
            let k = inst.record("B_SPLINE_SURFACE_WITH_KNOTS").ok_or_else(|| invalid("B-spline surfaces without knots are not read yet".into()))?;
            (&b[0], &b[1], &b[2], k)
        };
        let int = |v: &Value| v.as_integer().map(|i| i as usize).ok_or_else(|| invalid("an integer".into()));
        let (ud, vd) = (int(ud)?, int(vd)?);
        let poles: Vec<Vec<Point>> = self.list(grid)?
            .iter()
            .map(|row| self.list(row)?.iter().map(|v| self.point(Self::reference(v)?)).collect())
            .collect::<Result<_>>()?;
        let ints = |v: &Value| -> Result<Vec<usize>> { self.list(v)?.iter().map(int).collect() };
        let reals = |v: &Value| -> Result<Vec<f64>> { self.list(v)?.iter().map(|v| self.real(v)).collect() };
        let u = KnotVector::new(reals(&knots[2])?, ints(&knots[0])?, ud);
        let v = KnotVector::new(reals(&knots[3])?, ints(&knots[1])?, vd);
        let surface = BSplineSurface::try_new(poles, u, v)?;
        Ok(Surface::BSpline(match inst.record("RATIONAL_B_SPLINE_SURFACE") {
            Some(w) => surface.with_weights(self.list(&w[0])?.iter().map(reals).collect::<Result<_>>()?),
            None => surface,
        }))
    }

    // --- topology ------------------------------------------------------------------------

    fn vertex(&mut self, id: usize) -> Result<VertexId> {
        if let Some(&v) = self.vertices.get(&id) {
            return Ok(v);
        }
        let p = self.get(id)?.record("VERTEX_POINT").ok_or_else(|| invalid(format!("#{id} is not a vertex")))?;
        let point = self.point(Self::reference(&p[1])?)?;
        let v = self.model.add_vertex(point, crate::tolerance::LINEAR);
        self.vertices.insert(id, v);
        Ok(v)
    }

    fn edge(&mut self, id: usize) -> Result<EdgeMade> {
        if let Some(e) = self.edges.get(&id) {
            return Ok(e.clone());
        }
        let p = self.get(id)?.record("EDGE_CURVE").ok_or_else(|| invalid(format!("#{id} is not an edge")))?.to_vec();
        let (v1, v2) = (Self::reference(&p[1])?, Self::reference(&p[2])?);
        let geometry = Self::reference(&p[3])?;
        let same_sense = p[4].as_bool().unwrap_or(true);
        let (curve, scale) = self.curve(geometry)?;
        // The pcurves the file gives with the curve.
        let mut pcurves = Vec::new();
        let g = self.get(geometry)?;
        if let Some((_, list)) = self.surface_curve(&g) {
            for item in self.list(&list)? {
                let item = Self::reference(item)?;
                if let Some(pc) = self.get(item)?.record("PCURVE") {
                    pcurves.push((Self::reference(&pc[1])?, item));
                }
            }
        }
        let (a, b) = (self.vertex(v1)?, self.vertex(v2)?);
        let (pa, pb) = (self.model.point(Shape::from(a)), self.model.point(Shape::from(b)));
        // The model's edge runs along the curve: from the first vertex if the file's edge
        // does, else from the second.
        let (start, end, ps, pe, flipped) = if same_sense { (a, b, pa, pb, false) } else { (b, a, pb, pa, true) };
        let (first, last) = range(&curve, ps, pe, start == end);
        let edge = self.model.add_edge(start, end, crate::tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), curve.clone(), crate::tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), first, last)?;
        let made = EdgeMade { id: edge, flipped, curve, first, last, scale, pcurves };
        self.edges.insert(id, made.clone());
        self.edge_instances.insert(edge.into(), id);
        Ok(made)
    }

    fn face(&mut self, id: usize) -> Result<Shape> {
        let inst = self.get(id)?;
        let p = inst
            .record("ADVANCED_FACE")
            .or(inst.record("FACE_SURFACE"))
            .ok_or_else(|| invalid(format!("#{id}: faces of type {} are not read yet", inst.name())))?
            .to_vec();
        let surface_id = Self::reference(&p[2])?;
        let surface = self.surface(surface_id)?;
        let same_sense = p[3].as_bool().unwrap_or(true);
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface.clone(), crate::tolerance::LINEAR)?;
        self.face_instances.insert(face.node, id);
        // Boundaries: (outer?, coedges with pcurves).
        let mut loops: Vec<(bool, Vec<(EdgeId, bool, Curve2, f64, f64)>)> = Vec::new();
        for bound in self.list(&p[1])? {
            let bound = Self::reference(bound)?;
            let b = self.get(bound)?;
            let outer = b.has("FACE_OUTER_BOUND");
            let bp = b.record("FACE_OUTER_BOUND").or(b.record("FACE_BOUND")).ok_or_else(|| invalid(format!("#{bound} is not a bound")))?.to_vec();
            let along = bp[2].as_bool().unwrap_or(true);
            let lp = self.get(Self::reference(&bp[1])?)?;
            let Some(edges) = lp.record("EDGE_LOOP") else {
                continue; // a vertex loop (a pole): nothing to bound
            };
            let mut coedges = Vec::new();
            let mut seen: HashMap<EdgeId, usize> = HashMap::default();
            for oe in self.list(&edges[1])? {
                let oe = self.get(Self::reference(oe)?)?.record("ORIENTED_EDGE").ok_or_else(|| invalid("an oriented edge".into()))?.to_vec();
                let forward = oe[4].as_bool().unwrap_or(true);
                let made = self.edge(Self::reference(&oe[3])?)?;
                // Which of the file's pcurves on this surface (a seam has two).
                let k = seen.entry(made.id).or_insert(0);
                let pcurve = self.coedge_pcurve(&made, surface_id, &surface, *k).map_err(|e| {
                    invalid(format!("face #{id}: edge #{} ({}): {e}", Self::reference(&oe[3]).unwrap_or(0), curve_kind(&made.curve)))
                })?;
                *k += 1;
                coedges.push((made.id, forward != made.flipped, pcurve, made.first, made.last));
            }
            if !along {
                coedges.reverse();
                for c in &mut coedges {
                    c.1 = !c.1;
                }
            }
            align_periods(&mut coedges, &surface);
            loops.push((outer, coedges));
        }
        finish_face(self.model, &mut self.gaps, face, &surface, loops, same_sense)
    }

    /// A face of a surface bounded by curves on it (no edges given: each boundary a
    /// composite curve of segments, each a curve with pcurves).
    fn curve_bounded_face(&mut self, id: usize) -> Result<Shape> {
        let p = self.get(id)?.record("CURVE_BOUNDED_SURFACE").ok_or_else(|| invalid(format!("#{id} is not a curve-bounded surface")))?.to_vec();
        let surface_id = Self::reference(&p[1])?;
        let surface = self.surface(surface_id)?;
        let boundaries = self.list(&p[2])?.to_vec();
        if boundaries.is_empty() {
            return crate::make::face_on_surface(self.model, surface);
        }
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface.clone(), crate::tolerance::LINEAR)?;
        let mut loops = Vec::new();
        for b in boundaries {
            let b = Self::reference(&b)?;
            let inst = self.get(b)?;
            let outer = inst.has("OUTER_BOUNDARY_CURVE");
            let bp = inst
                .record("OUTER_BOUNDARY_CURVE")
                .or(inst.record("BOUNDARY_CURVE"))
                .or(inst.record("COMPOSITE_CURVE_ON_SURFACE"))
                .or(inst.record("COMPOSITE_CURVE"))
                .ok_or_else(|| invalid(format!("#{b}: boundaries of type {} are not read yet", inst.name())))?
                .to_vec();
            let mut coedges = Vec::new();
            for segment in self.list(&bp[1])?.to_vec() {
                let sp = self.get(Self::reference(&segment)?)?.record("COMPOSITE_CURVE_SEGMENT").ok_or_else(|| invalid("a composite curve's segment".into()))?.to_vec();
                let along = sp[1].as_bool().unwrap_or(true);
                let made = self.segment_edge(Self::reference(&sp[2])?)?;
                let pcurve = self.coedge_pcurve(&made, surface_id, &surface, 0).map_err(|e| invalid(format!("curve-bounded surface #{id}: {e}")))?;
                coedges.push((made.id, along, pcurve, made.first, made.last));
            }
            align_periods(&mut coedges, &surface);
            loops.push((outer, coedges));
        }
        finish_face(self.model, &mut self.gaps, face, &surface, loops, true)
    }

    /// An edge along a boundary curve's segment (a trimmed curve, possibly with pcurves),
    /// between vertices at its ends shared with the other segments meeting there.
    fn segment_edge(&mut self, id: usize) -> Result<EdgeMade> {
        if let Some(e) = self.edges.get(&id) {
            return Ok(e.clone());
        }
        let (curve, scale) = self.curve(id)?;
        let mut pcurves = Vec::new();
        let g = self.get(id)?;
        if let Some((_, list)) = self.surface_curve(&g) {
            for item in self.list(&list)? {
                let item = Self::reference(item)?;
                if let Some(pc) = self.get(item)?.record("PCURVE") {
                    pcurves.push((Self::reference(&pc[1])?, item));
                }
            }
        }
        let (first, last) = self.trim_range(id, &curve, scale)?;
        let (ps, pe) = (curve.point(first), curve.point(last));
        let (start, end) = (self.vertex_at(ps), self.vertex_at(pe));
        let edge = self.model.add_edge(start, end, crate::tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), curve.clone(), crate::tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), first, last)?;
        let made = EdgeMade { id: edge, flipped: false, curve, first, last, scale, pcurves };
        self.edges.insert(id, made.clone());
        Ok(made)
    }

    /// The model's range of a trimmed curve (inside a surface curve or not), increasing; the
    /// whole curve where it is bounded and not trimmed. (Which way a boundary runs along it
    /// is the segment's to say.)
    fn trim_range(&self, id: usize, curve: &Curve, scale: f64) -> Result<(f64, f64)> {
        let inst = self.get(id)?;
        if let Some((basis, _)) = self.surface_curve(&inst) {
            return self.trim_range(Self::reference(&basis)?, curve, scale);
        }
        if let Some(p) = inst.record("TRIMMED_CURVE") {
            // Each trim a parameter value (else a point, found on the curve).
            let trim = |v: &Value| -> Result<f64> {
                for item in self.list(v)? {
                    if let Some(("PARAMETER_VALUE", _)) = self.x.typed(item) {
                        // (A conic's parameter is an angle, in the file's angle unit.)
                        let unit = if matches!(curve, Curve::Circle(_) | Curve::Ellipse(_)) { self.angle } else { scale };
                        return Ok(self.real(item)? * unit);
                    }
                }
                for item in self.list(v)? {
                    if let Some(r) = item.as_reference() {
                        let q = self.point(r)?;
                        let (lo, hi) = (curve.first_parameter().max(-1e6), curve.last_parameter().min(1e6));
                        return Ok(nearest_on_curve(q, curve, lo, hi).parameter);
                    }
                }
                Err(invalid(format!("#{id}: a trim of no parameter or point")))
            };
            // From the first trim to the second: on a periodic curve round the way its sense
            // says.
            let (t1, t2) = (trim(&p[2])?, trim(&p[3])?);
            let sense = p[4].as_bool().unwrap_or(true);
            return Ok(match curve.period() {
                Some(period) => {
                    let (a, b) = if sense { (t1, t2) } else { (t2, t1) };
                    (a, if b <= a { b + period * ((a - b) / period).ceil().max(1.0) } else { b })
                }
                None => (t1.min(t2), t1.max(t2)),
            });
        }
        let (lo, hi) = (curve.first_parameter(), curve.last_parameter());
        if lo.is_finite() && hi.is_finite() { Ok((lo, hi)) } else { Err(invalid(format!("#{id}: an unbounded curve bounds a surface"))) }
    }

    /// The vertex at `p` (made once for points within the linear tolerance).
    fn vertex_at(&mut self, p: Point) -> VertexId {
        let near = 1e-6 * (1.0 + (p - Point::ORIGIN).length());
        if let Some(&(_, v)) = self.point_vertices.iter().find(|(q, _)| q.distance(p) <= near) {
            return v;
        }
        let v = self.model.add_vertex(p, crate::tolerance::LINEAR);
        self.point_vertices.push((p, v));
        v
    }

    /// The pcurve of an edge on a face's surface: the file's (the `k`-th on it: a seam has
    /// two) where it follows the curve, else the curve projected, else the file's after all
    /// (the tolerances then cover the gap).
    fn coedge_pcurve(&mut self, made: &EdgeMade, surface_id: usize, surface: &Surface, k: usize) -> Result<Curve2> {
        let on: Vec<usize> = made.pcurves.iter().filter(|(s, _)| *s == surface_id).map(|&(_, c)| c).collect();
        let given = match on.get(k % on.len().max(1)) {
            Some(&c) => self.curve2(c, made.scale, made.first, made.last).ok(),
            None => None,
        };
        Ok(match given {
            Some(pc) if fits(&pc, &made.span(), surface) => pc,
            // Parametrized otherwise but along the edge on a free-form surface (where
            // projecting is costly): kept, the tolerances covering the difference.
            Some(pc) if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) && follows(&pc, &made.span(), surface) => pc,
            given => match (project(&made.span(), surface, given.is_some()), given) {
                (Ok(pc), _) => pc,
                (Err(_), Some(pc)) => pc,
                (Err(e), None) => return Err(e),
            },
        })
    }

    fn shell(&mut self, id: usize) -> Result<Shape> {
        let inst = self.get(id)?;
        let p = inst
            .record("CLOSED_SHELL")
            .or(inst.record("OPEN_SHELL"))
            .or(inst.record("ORIENTED_CLOSED_SHELL"))
            .ok_or_else(|| invalid(format!("#{id} is not a shell")))?
            .to_vec();
        if inst.has("ORIENTED_CLOSED_SHELL") {
            let s = self.shell(Self::reference(&p[2])?)?;
            return Ok(if p[3].as_bool().unwrap_or(true) { s } else { s.reversed() });
        }
        let shell = self.model.new_empty(NodeKind::Shell);
        for f in self.list(&p[1])? {
            let face = self.face(Self::reference(f)?)?;
            self.model.add(shell, face)?;
        }
        Ok(shell)
    }

    /// A solid or shell of a representation (`None` for other items, such as placements).
    fn body(&mut self, id: usize) -> Result<Option<Shape>> {
        if let Some(&s) = self.bodies.get(&id) {
            return Ok(Some(s));
        }
        let inst = self.get(id)?;
        let made = if let Some(p) = inst.record("MANIFOLD_SOLID_BREP").or(inst.record("FACETED_BREP")) {
            let p = p.to_vec();
            let solid = self.model.new_empty(NodeKind::Solid);
            let outer = self.shell(Self::reference(&p[1])?)?;
            self.model.add(solid, outer)?;
            if let Some(voids) = inst.record("BREP_WITH_VOIDS") {
                for v in self.list(&voids[0])?.to_vec() {
                    let void = self.shell(Self::reference(&v)?)?;
                    self.model.add(solid, void)?;
                }
            }
            solid
        } else if let Some(p) = inst.record("SHELL_BASED_SURFACE_MODEL") {
            let shells = self.list(&p[1])?.to_vec();
            if shells.len() == 1 {
                self.shell(Self::reference(&shells[0])?)?
            } else {
                let c = self.model.new_empty(NodeKind::Compound);
                for s in shells {
                    let shell = self.shell(Self::reference(&s)?)?;
                    self.model.add(c, shell)?;
                }
                c
            }
        } else if let Some(p) = inst.record("GEOMETRIC_SET").or(inst.record("GEOMETRIC_CURVE_SET")) {
            // Its surfaces as faces (curves and points left out): bounded ones whole, those
            // bounded by curves on them by those.
            let c = self.model.new_empty(NodeKind::Compound);
            let mut any = false;
            for item in self.list(&p[1])?.to_vec() {
                let item = Self::reference(&item)?;
                let it = self.get(item)?;
                let face = if it.has("CURVE_BOUNDED_SURFACE") {
                    self.curve_bounded_face(item)?
                } else if it.has("B_SPLINE_SURFACE") || it.has("B_SPLINE_SURFACE_WITH_KNOTS") || it.has("RECTANGULAR_TRIMMED_SURFACE") || it.has("BEZIER_SURFACE") {
                    let surface = self.surface(item)?;
                    crate::make::face_on_surface(self.model, surface)?
                } else {
                    continue;
                };
                self.model.add(c, face)?;
                any = true;
            }
            if !any {
                return Ok(None);
            }
            c
        } else {
            return Ok(None);
        };
        self.bodies.insert(id, made);
        Ok(Some(made))
    }

    // --- products ------------------------------------------------------------------------

    fn assembly(&mut self) -> Result<Vec<Node>> {
        let x = self.x;
        // Product definitions: their names and representations.
        let mut reps_of: HashMap<usize, Vec<usize>> = HashMap::default();
        for (_, inst) in x.of_type("SHAPE_DEFINITION_REPRESENTATION") {
            let p = inst.params();
            let pds = self.get(Self::reference(&p[0])?)?;
            if let Some(d) = pds.record("PRODUCT_DEFINITION_SHAPE") {
                if let Some(pd) = d[2].as_reference() {
                    reps_of.entry(pd).or_default().push(Self::reference(&p[1])?);
                }
            }
        }
        // Representations linked without a placement (a part's shape representation and its
        // B-rep one).
        let mut linked: HashMap<usize, Vec<usize>> = HashMap::default();
        for (_, inst) in x.of_type("SHAPE_REPRESENTATION_RELATIONSHIP") {
            if inst.has("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION") {
                continue;
            }
            let p = inst.record("REPRESENTATION_RELATIONSHIP").or(inst.record("SHAPE_REPRESENTATION_RELATIONSHIP")).unwrap_or(&[]);
            if p.len() >= 4 {
                let (a, b) = (Self::reference(&p[2])?, Self::reference(&p[3])?);
                linked.entry(a).or_default().push(b);
                linked.entry(b).or_default().push(a);
            }
        }
        // Occurrences: parent definition -> (child definition, occurrence, name).
        let mut children: HashMap<usize, Vec<(usize, usize, String)>> = HashMap::default();
        let mut is_child: crate::hashing::FastSet<usize> = Default::default();
        for (id, inst) in x.of_type("NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            let p = inst.params();
            let (parent, child) = (Self::reference(&p[3])?, Self::reference(&p[4])?);
            let name = p[1].as_str().unwrap_or("").to_string();
            children.entry(parent).or_default().push((child, id, name));
            is_child.insert(child);
        }
        // Placements of occurrences.
        let mut placement_of: HashMap<usize, Transform> = HashMap::default();
        for (_, inst) in x.of_type("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") {
            let p = inst.params();
            let rel = self.get(Self::reference(&p[0])?)?;
            let pds = self.get(Self::reference(&p[1])?)?;
            let Some(occurrence) = pds.record("PRODUCT_DEFINITION_SHAPE").and_then(|d| d[2].as_reference()) else { continue };
            let Some(t) = rel.record("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION") else { continue };
            let idt = self.get(Self::reference(&t[0])?)?;
            let Some(items) = idt.record("ITEM_DEFINED_TRANSFORMATION") else { continue };
            let (from, to) = (self.frame(Self::reference(&items[2])?)?, self.frame(Self::reference(&items[3])?)?);
            placement_of.insert(occurrence, Transform::frame_to_frame(&from, &to));
        }
        let definitions: Vec<usize> = x.of_type("PRODUCT_DEFINITION").map(|(id, _)| id).collect();
        let roots: Vec<usize> = definitions.iter().copied().filter(|d| !is_child.contains(d)).collect();
        let mut nodes = Vec::new();
        for root in roots {
            let name = self.product_name(root)?;
            nodes.push(self.node(root, name, Transform::IDENTITY, &reps_of, &linked, &children, &placement_of, 0)?);
        }
        // Files without products: every body.
        if nodes.is_empty() {
            let mut shapes = Vec::new();
            for id in x.ids() {
                if let Some(s) = self.body(id)? {
                    shapes.push(s);
                }
            }
            let shape = self.compound_of(shapes)?;
            nodes.push(Node { shape, ..Default::default() });
        }
        Ok(nodes)
    }

    fn product_name(&self, definition: usize) -> Result<String> {
        let pd = self.get(definition)?.params();
        let formation = self.get(Self::reference(&pd[2])?)?.params();
        let product = self.get(Self::reference(&formation[2])?)?.params();
        Ok(product[0].as_str().or(product[1].as_str()).unwrap_or("").to_string())
    }

    #[allow(clippy::too_many_arguments)]
    fn node(
        &mut self,
        definition: usize,
        name: String,
        placement: Transform,
        reps_of: &HashMap<usize, Vec<usize>>,
        linked: &HashMap<usize, Vec<usize>>,
        children: &HashMap<usize, Vec<(usize, usize, String)>>,
        placement_of: &HashMap<usize, Transform>,
        depth: usize,
    ) -> Result<Node> {
        if depth > 64 {
            return Err(invalid("an assembly nested too deep (a cycle?)".into()));
        }
        // The product's bodies: items of its representations and of those linked to them.
        let mut reps: Vec<usize> = reps_of.get(&definition).cloned().unwrap_or_default();
        let mut i = 0;
        while i < reps.len() {
            for &r in linked.get(&reps[i]).map(Vec::as_slice).unwrap_or(&[]) {
                if !reps.contains(&r) {
                    reps.push(r);
                }
            }
            i += 1;
        }
        let mut bodies = Vec::new();
        let mut colors = Vec::new();
        for rep in reps {
            let inst = self.get(rep)?;
            let Some(items) = inst.records().find(|(_, v)| v.len() >= 2 && self.x.list(&v[1]).is_some()).map(|(_, v)| v[1].clone()) else { continue };
            for item in self.list(&items)? {
                let item = Self::reference(item)?;
                if let Some(body) = self.body(item)? {
                    if !bodies.iter().any(|b: &Shape| b.node == body.node) {
                        bodies.push(body);
                        colors.push(self.colors.get(&item).copied());
                    }
                }
            }
        }
        let shared = match colors.first() {
            Some(&Some(c)) if colors.iter().all(|&k| k == Some(c)) => Some(c),
            _ => None,
        };
        let (color, alpha) = (shared.map(|[r, g, b, _]| [r, g, b]), shared.map(|c| c[3]).filter(|&a| a < 1.0));
        // The faces and edges with names, colours or layers of their own.
        let mut subshapes = Vec::new();
        for &body in &bodies {
            for (kind, instances) in [(NodeKind::Face, &self.face_instances), (NodeKind::Edge, &self.edge_instances)] {
                for s in self.model.unique_subshapes(body, kind) {
                    let Some(&inst) = instances.get(&s.node) else { continue };
                    let name = self.x.get(inst).ok().and_then(|i| i.params().first().and_then(Value::as_str).map(str::to_string)).filter(|n| !n.is_empty() && n != "NONE");
                    let color = if kind == NodeKind::Face { self.colors.get(&inst).copied() } else { None };
                    let layer = self.layers.get(&inst).cloned();
                    if name.is_some() || color.is_some() || layer.is_some() {
                        subshapes.push(Subshape { shape: Shape { orientation: crate::topology::Orientation::Forward, ..s }, name, color, layer });
                    }
                }
            }
        }
        let material = self.materials.get(&definition).cloned();
        let shape = self.compound_of(bodies)?;
        let mut kids = Vec::new();
        for (child, occurrence, occurrence_name) in children.get(&definition).cloned().unwrap_or_default() {
            let child_name = if occurrence_name.is_empty() { self.product_name(child)? } else { occurrence_name };
            let at = placement_of.get(&occurrence).copied().unwrap_or(Transform::IDENTITY);
            kids.push(self.node(child, child_name, at, reps_of, linked, children, placement_of, depth + 1)?);
        }
        Ok(Node { name, shape, color, alpha, placement, children: kids, subshapes, material })
    }

    fn compound_of(&mut self, shapes: Vec<Shape>) -> Result<Option<Shape>> {
        Ok(match shapes.len() {
            0 => None,
            1 => Some(shapes[0]),
            _ => {
                let c = self.model.new_empty(NodeKind::Compound);
                for s in shapes {
                    self.model.add(c, s)?;
                }
                Some(c)
            }
        })
    }
}

/// Knots of the uniform and quasi-uniform B-spline types (which leave them implicit).
fn uniform_knots(inst: &Instance, poles: usize, degree: usize) -> (Vec<usize>, Vec<f64>) {
    if inst.has("QUASI_UNIFORM_CURVE") {
        let spans = poles - degree;
        let mut mults = vec![1; spans + 1];
        mults[0] = degree + 1;
        mults[spans] = degree + 1;
        (mults, (0..=spans).map(|i| i as f64).collect())
    } else {
        let n = poles + degree + 1;
        (vec![1; n], (0..n).map(|i| i as f64 - degree as f64).collect())
    }
}

/// The length and angle units of the file's geometry: millimetres and radians per unit.
fn units(x: &Exchange) -> (f64, f64) {
    let mut length = 1.0;
    let mut angle = 1.0;
    for (_, inst) in x.of_type("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
        let Some(list) = inst.record("GLOBAL_UNIT_ASSIGNED_CONTEXT").and_then(|p| x.list(&p[0])) else { continue };
        for unit in list.iter().filter_map(Value::as_reference) {
            let Ok(u) = x.get(unit) else { continue };
            if u.has("LENGTH_UNIT") {
                length = unit_factor(x, &u, true);
            } else if u.has("PLANE_ANGLE_UNIT") {
                angle = unit_factor(x, &u, false);
            }
        }
        break;
    }
    (length, angle)
}

/// Millimetres (lengths) or radians (angles) per unit.
fn unit_factor(x: &Exchange, u: &Instance, is_length: bool) -> f64 {
    if let Some(si) = u.record("SI_UNIT") {
        let prefix = match si[0].as_enumeration() {
            Some("MILLI") => 1e-3,
            Some("CENTI") => 1e-2,
            Some("DECI") => 1e-1,
            Some("KILO") => 1e3,
            Some("MICRO") => 1e-6,
            Some("NANO") => 1e-9,
            _ => 1.0,
        };
        return if is_length { prefix * 1e3 } else { prefix };
    }
    if let Some(c) = u.record("CONVERSION_BASED_UNIT") {
        // CONVERSION_BASED_UNIT(name, MEASURE_WITH_UNIT(value, unit)).
        let Some(m) = c[1].as_reference().and_then(|r| x.get(r).ok()) else { return 1.0 };
        let Some(mw) = m.records().find(|(n, _)| n.ends_with("MEASURE_WITH_UNIT")).map(|(_, v)| v) else { return 1.0 };
        let value = x.real(&mw[0]).unwrap_or(1.0);
        let base = mw[1].as_reference().and_then(|r| x.get(r).ok()).map_or(1.0, |b| unit_factor(x, &b, is_length));
        return value * base;
    }
    1.0
}

/// The colours of styled items: item instance -> RGB.
fn colors(x: &Exchange) -> HashMap<usize, [f64; 4]> {
    let mut out = HashMap::default();
    for (_, inst) in x.of_type("STYLED_ITEM") {
        let p = inst.record("STYLED_ITEM").expect("record");
        let Some(item) = p.get(2).and_then(Value::as_reference) else { continue };
        let mut stack: Vec<&Value> = vec![&p[1]];
        let mut depth = 0;
        let mut found = None;
        // (The transparency of its rendering, where it has one.)
        let mut transparency: Option<f64> = None;
        while let Some(v) = stack.pop() {
            depth += 1;
            if depth > 200 {
                break;
            }
            match v {
                Value::List(_) => stack.extend(x.list(v).unwrap_or(&[])),
                Value::Typed(_) => stack.extend(x.typed(v).map_or(&[][..], |(_, params)| params)),
                Value::Reference(r) => {
                    let Ok(i) = x.get(*r) else { continue };
                    if let Some(t) = i.record("SURFACE_STYLE_TRANSPARENT") {
                        transparency = transparency.or(t.first().and_then(|v| x.real(v)));
                    }
                    if let Some(c) = i.record("COLOUR_RGB") {
                        let rgb = [x.real(&c[1]), x.real(&c[2]), x.real(&c[3])];
                        if let [Some(r), Some(g), Some(b)] = rgb {
                            found = found.or(Some([r, g, b]));
                            continue;
                        }
                    }
                    if let Some(c) = i.record("DRAUGHTING_PRE_DEFINED_COLOUR") {
                        if let Some(c) = c[0].as_str().and_then(named_colour) {
                            found = found.or(Some(c));
                            continue;
                        }
                    }
                    for (_, params) in i.records() {
                        stack.extend(params.iter());
                    }
                }
                _ => {}
            }
        }
        if let Some([r, g, b]) = found {
            out.entry(item).or_insert([r, g, b, 1.0 - transparency.unwrap_or(0.0).clamp(0.0, 1.0)]);
        }
    }
    out
}

/// The layers items are on: item -> the layer's name.
fn layers(x: &Exchange) -> HashMap<usize, String> {
    let mut out = HashMap::default();
    for (_, inst) in x.of_type("PRESENTATION_LAYER_ASSIGNMENT") {
        let p = inst.params();
        let Some(name) = p.first().and_then(Value::as_str) else { continue };
        for item in p.get(2).and_then(|v| x.list(v)).unwrap_or(&[]) {
            if let Some(r) = item.as_reference() {
                out.entry(r).or_insert_with(|| name.to_string());
            }
        }
    }
    out
}

/// The products' materials: product definition -> material (its designation's name; its
/// description, density and density's unit where a property gives them).
fn materials(x: &Exchange) -> HashMap<usize, Material> {
    let definition_of = |pds: usize| x.get(pds).ok().and_then(|i| i.record("PRODUCT_DEFINITION_SHAPE").and_then(|d| d.get(2)?.as_reference()));
    let mut out: HashMap<usize, Material> = HashMap::default();
    for (_, inst) in x.of_type("MATERIAL_DESIGNATION") {
        let p = inst.params();
        let Some(name) = p.first().and_then(Value::as_str) else { continue };
        // (Of one definition, or a set of them.)
        let targets: Vec<usize> = match p.get(1) {
            Some(v) if x.list(v).is_some() => x.list(v).unwrap_or(&[]).iter().filter_map(Value::as_reference).collect(),
            Some(v) => v.as_reference().into_iter().collect(),
            None => Vec::new(),
        };
        for pds in targets {
            if let Some(pd) = definition_of(pds) {
                out.insert(pd, Material { name: name.to_string(), description: String::new(), density: 0.0, density_unit: String::new() });
            }
        }
    }
    // Properties of a material: description, density and its unit.
    for (_, inst) in x.of_type("PROPERTY_DEFINITION_REPRESENTATION") {
        let p = inst.params();
        let (Some(prop), Some(rep)) = (p.first().and_then(Value::as_reference), p.get(1).and_then(Value::as_reference)) else { continue };
        let Ok(prop) = x.get(prop) else { continue };
        let Some(pp) = prop.record("PROPERTY_DEFINITION") else { continue };
        if pp.first().and_then(Value::as_str) != Some("material property") {
            continue;
        }
        let Some(pd) = pp.get(2).and_then(Value::as_reference).and_then(definition_of) else { continue };
        let Some(material) = out.get_mut(&pd) else { continue };
        let Ok(rep) = x.get(rep) else { continue };
        let Some(items) = rep.records().find(|(_, v)| v.len() >= 2 && x.list(&v[1]).is_some()).and_then(|(_, v)| x.list(&v[1])) else { continue };
        for item in items.iter().filter_map(Value::as_reference) {
            let Ok(i) = x.get(item) else { continue };
            if let Some(d) = i.record("DESCRIPTIVE_REPRESENTATION_ITEM") {
                let (key, value) = (d.first().and_then(Value::as_str).unwrap_or(""), d.get(1).and_then(Value::as_str).unwrap_or(""));
                match key {
                    "description" => material.description = value.to_string(),
                    "density" => material.density = value.trim().parse().unwrap_or(material.density),
                    "density unit" => material.density_unit = value.to_string(),
                    _ => {}
                }
            } else if let Some(m) = i.record("MEASURE_REPRESENTATION_ITEM") {
                // (As OCCT writes a density: a measure with a unit.)
                if m.first().and_then(Value::as_str).is_some_and(|k| k.contains("density")) {
                    if let Some(v) = m.get(1).and_then(|v| x.real(v)) {
                        material.density = v;
                    }
                }
            }
        }
    }
    out
}

fn named_colour(name: &str) -> Option<[f64; 3]> {
    Some(match name {
        "red" => [1.0, 0.0, 0.0],
        "green" => [0.0, 1.0, 0.0],
        "blue" => [0.0, 0.0, 1.0],
        "yellow" => [1.0, 1.0, 0.0],
        "magenta" => [1.0, 0.0, 1.0],
        "cyan" => [0.0, 1.0, 1.0],
        "black" => [0.0, 0.0, 0.0],
        "white" => [1.0, 1.0, 1.0],
        _ => return None,
    })
}

/// The parameter range along `curve` from `start` to `end` (a whole period when the edge
/// closes on itself).
fn range(curve: &Curve, start: Point, end: Point, closed: bool) -> (f64, f64) {
    let period = curve.period();
    let param = |p: Point| match curve {
        Curve::Line(l) => (p - l.origin).dot(Vector::from(l.direction)),
        _ => {
            let (lo, hi) = match period {
                Some(t) => (curve.first_parameter(), curve.first_parameter() + t),
                None => (curve.first_parameter(), curve.last_parameter()),
            };
            // Mostly the vertices are the curve's ends.
            let exact = 1e-12 * (1.0 + p.x.abs().max(p.y.abs()).max(p.z.abs()));
            if lo.is_finite() && curve.point(lo).distance(p) <= exact {
                return lo;
            }
            if period.is_none() && hi.is_finite() && curve.point(hi).distance(p) <= exact {
                return hi;
            }
            nearest_on_curve(p, curve, lo, hi).parameter
        }
    };
    let mut first = param(start);
    // Into [0, period): a point at the curve's origin is found at either end of the period,
    // and the file's pcurves expect 0.
    if let Some(t) = period {
        let origin = curve.first_parameter();
        first = origin + (first - origin).rem_euclid(t);
        if t - (first - origin) <= 1e-9 * t {
            first = origin;
        }
    }
    let mut last = if closed {
        match period {
            Some(t) => first + t,
            None => curve.last_parameter(),
        }
    } else {
        param(end)
    };
    if closed && period.is_none() {
        return (curve.first_parameter(), curve.last_parameter());
    }
    if let Some(t) = period {
        while last <= first + 1e-12 {
            last += t;
        }
        while last - first > t + 1e-12 {
            last -= t;
        }
    }
    (first, last)
}

fn curve_kind(c: &Curve) -> &'static str {
    match c {
        Curve::Line(_) => "line",
        Curve::Circle(_) => "circle",
        Curve::Ellipse(_) => "ellipse",
        Curve::BSpline(_) => "B-spline",
        _ => "curve",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CUBE: &str = include_str!("tests/cube.step");

    #[test]
    fn a_cube_with_colour_and_name() {
        let mut model = Model::new();
        // Written by CadQuery: a root product holding one occurrence of `cube`.
        let roots = read_assembly(&mut model, &mut CUBE.as_bytes()).unwrap();
        assert_eq!(roots.len(), 1);
        let cube = &roots[0].children[0];
        assert_eq!(cube.name, "cube");
        let shape = cube.shape.unwrap();
        assert_eq!(shape.kind(), NodeKind::Solid);
        let volume = crate::properties::volume(&model, shape).mass;
        assert!((volume - 1000.0).abs() < 1e-6, "volume {volume}");
        assert!(crate::validity::check(&model, shape).errors().next().is_none());
        assert_eq!(cube.color, Some([1.0, 0.0, 0.0]));
    }
}
