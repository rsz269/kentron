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
use crate::geometry::curve::{BSplineCurve, Circle, Curve, Ellipse, Hyperbola, Line, Parabola};
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Hyperbola2, Line2, Parabola2};
use crate::geometry::extrema::nearest_on_curve;
use crate::geometry::surface::{
    BSplineSurface, Cone, Cylinder, ExtrusionSurface, KnotVector, OffsetSurface, Plane, RevolutionSurface, Sphere, Surface, Torus,
};
use crate::geometry::ParametricCurve;
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Vector, Vector2};
use crate::topology::{EdgeId, Model, NodeKind, Shape, VertexId};

use super::faces::{align_periods, finish_face, fits, follows, project, split_at_poles, CurveSpan, Gaps};

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
    /// Parts of the shape with names, colours or layers of their own: faces and edges, and
    /// the shells and solids of a shape of several (a wire's layer is its edges').
    pub subshapes: Vec<Subshape>,
    /// What the product is made of.
    pub material: Option<Material>,
    /// The layer the product's shape is on as a whole, if any.
    pub layer: Option<String>,
    /// The validation properties the file gives the product.
    pub validation: Validation,
    /// Those it gives this occurrence of it (in its parent).
    pub occurrence_validation: Validation,
    /// The file of the external reference the file gives the product, by name (not read).
    pub external: Option<String>,
    /// This occurrence's own colour (red, green, blue, opacity), over the product's: a style
    /// the file gives it in the context of its parent.
    pub occurrence_color: Option<[f64; 4]>,
    /// This occurrence's own layer, if the file puts it on one.
    pub occurrence_layer: Option<String>,
    /// The file's product definition this node is an occurrence of (its entity number): the
    /// same for every occurrence of one product. `None` for shapes read without products, and
    /// for nodes made otherwise than by reading.
    pub definition: Option<usize>,
}

/// Properties a file gives to check its reading by (geometric validation properties): the
/// centroid, volume and area of a product's shape, or of an occurrence of one.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Validation {
    pub centroid: Option<Point>,
    pub volume: Option<f64>,
    pub area: Option<f64>,
}

/// A face, wire, edge or vertex of a product's shape with a name, colour or layer of its own.
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
    // (Edges meeting at a shallow angle and crossing just beside their vertex, within the
    // reading's largest tolerance (read.maxprecision.val, 1): the vertex's tolerance raised
    // to cover it, as OCCT's ShapeFix_Wire::FixIntersectingEdges does.)
    let mut faces: Vec<Shape> = reader.face_instances.keys().map(|&node| Shape { node, orientation: crate::topology::Orientation::Forward, placement: Transform::IDENTITY }).collect();
    faces.sort_by_key(|f| f.node.index);
    crate::validity::cover_crossings_at_vertices(reader.model, &faces, 1.0);
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

/// A triangulated face or surface set of the file.
struct Triangulation {
    /// Its coordinates list's instance.
    list: usize,
    /// Each node's point in that list (0-based), and the point.
    at: Vec<usize>,
    nodes: Vec<Point>,
    /// Over the nodes, 0-based.
    triangles: Vec<[u32; 3]>,
    normals: Option<Vec<Direction>>,
    /// The face it is linked to.
    link: Option<usize>,
}

struct Reader<'a> {
    x: &'a Exchange,
    model: &'a mut Model,
    /// Millimetres per length unit of the file.
    length: f64,
    /// Radians per angle unit of the file.
    angle: f64,
    /// The file's distance uncertainty, in millimetres (OCCT's reading precision).
    precision: f64,
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
    /// Validation properties by the product definition or occurrence they are given to.
    validations: HashMap<usize, Validation>,
    /// The file's instance of each face and edge read (for their names, colours and layers).
    face_instances: HashMap<crate::topology::NodeId, usize>,
    edge_instances: HashMap<crate::topology::NodeId, usize>,
    shell_instances: HashMap<crate::topology::NodeId, usize>,
    /// The file's instance of each wire (a composite curve) and vertex (a point) of a
    /// geometric set.
    item_instances: HashMap<crate::topology::NodeId, usize>,
    /// External references: product definition -> the referenced file's name.
    externals: HashMap<usize, String>,
    /// Occurrences' own colours (styles by context, or on an occurrence's mapped item).
    occurrence_colors: HashMap<usize, [f64; 4]>,
    /// Occurrences' own layers (those of an occurrence's mapped item).
    occurrence_layers: HashMap<usize, String>,
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
            precision: uncertainty(x),
            vertices: HashMap::default(),
            point_vertices: Vec::new(),
            edges: HashMap::default(),
            surfaces: HashMap::default(),
            bodies: HashMap::default(),
            colors: HashMap::default(),
            layers: HashMap::default(),
            materials: HashMap::default(),
            validations: HashMap::default(),
            face_instances: HashMap::default(),
            edge_instances: HashMap::default(),
            shell_instances: HashMap::default(),
            item_instances: HashMap::default(),
            externals: HashMap::default(),
            occurrence_colors: HashMap::default(),
            occurrence_layers: HashMap::default(),
            gaps: Gaps::default(),
        };
        reader.colors = colors(x);
        reader.externals = external_references(x);
        reader.occurrence_colors = occurrence_colors(x);
        reader.layers = layers(x);
        // (An occurrence shaped by a representation of its own holding the mapped item that
        // places its product: what the file gives that item is the occurrence's, as OCCT
        // puts it on the component's label.)
        for (item, occurrence) in occurrence_items(x) {
            if let Some(&c) = reader.colors.get(&item) {
                reader.occurrence_colors.entry(occurrence).or_insert(c);
            }
            if let Some(l) = reader.layers.get(&item) {
                reader.occurrence_layers.entry(occurrence).or_insert_with(|| l.clone());
            }
        }
        reader.materials = materials(x);
        reader.validations = validations(x, length);
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

    /// Whether the geometry of `id` is in 3D: no point it is built on has other than three
    /// coordinates (OCCT's StepToGeom makes nothing of a 2D point, nor of what is built on
    /// one).
    fn three_d(&self, id: usize) -> bool {
        let mut stack = vec![id];
        let mut seen: crate::hashing::FastSet<usize> = Default::default();
        while let Some(i) = stack.pop() {
            if !seen.insert(i) || seen.len() > 10_000 {
                continue;
            }
            let Ok(inst) = self.get(i) else { continue };
            if inst.has("CARTESIAN_POINT") {
                if self.coordinates(i).map_or(true, |c| c.len() != 3) {
                    return false;
                }
                continue;
            }
            let mut values: Vec<&Value> = inst.records().flat_map(|(_, v)| v.iter()).collect();
            while let Some(v) = values.pop() {
                match v {
                    Value::Reference(r) => stack.push(*r),
                    Value::List(_) => values.extend(self.x.list(v).unwrap_or(&[])),
                    _ => {}
                }
            }
        }
        true
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
        // (Conics' parameters as kentron's: edges set their ranges from their vertices.)
        if let Some(p) = inst.record("HYPERBOLA") {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok((Curve::Hyperbola(Hyperbola::new(frame, self.real(&p[2])? * s, self.real(&p[3])? * s)), 1.0));
        }
        if let Some(p) = inst.record("PARABOLA") {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok((Curve::Parabola(Parabola::new(frame, self.real(&p[2])? * s)), 1.0));
        }
        if is_bspline_curve(&inst) {
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
        } else if !inst.is_complex() {
            // (A simple uniform, quasi-uniform or Bézier curve, or a B-spline one without knots:
            // name, degree, poles, form, closed, self_intersect; its knots implicit.)
            let p = inst.params();
            (&p[1], &p[2], None)
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
            // (Simple, or complex with its items in the REPRESENTATION part.)
            let items = rep
                .record("DEFINITIONAL_REPRESENTATION")
                .filter(|p| p.len() >= 2)
                .or(rep.record("REPRESENTATION"))
                .ok_or_else(|| invalid("a pcurve without its representation".into()))?;
            let curve = items.get(1).map(|i| self.list(i)).transpose()?.and_then(|l| l.first()).ok_or_else(|| invalid("a pcurve's representation without a curve".into()))?;
            return self.curve2(Self::reference(curve)?, scale, first, last);
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
            if is_bspline_curve(&inst) {
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
        if let Some(p) = inst.record("HYPERBOLA") {
            return Ok(Curve2::Hyperbola(Hyperbola2::new(self.frame2(Self::reference(&p[1])?)?, self.real(&p[2])?, self.real(&p[3])?)));
        }
        if let Some(p) = inst.record("PARABOLA") {
            return Ok(Curve2::Parabola(Parabola2::new(self.frame2(Self::reference(&p[1])?)?, self.real(&p[2])?)));
        }
        if is_bspline_curve(&inst) {
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
        // (A degenerate torus, its flag choosing a spindle torus' inner or outer part, read as
        // the torus, as OCCT does: StepGeom_DegenerateToroidalSurface is a ToroidalSurface.)
        if let Some(p) = inst.record("TOROIDAL_SURFACE").or(inst.record("DEGENERATE_TOROIDAL_SURFACE")) {
            let frame = self.frame(Self::reference(&p[1])?)?;
            return Ok(Surface::Torus(Torus::new(frame, self.real(&p[2])? * l, self.real(&p[3])? * l)));
        }
        if let Some(p) = inst.record("RECTANGULAR_TRIMMED_SURFACE") {
            // As OCCT's StepToGeom::MakeRectangularTrimmedSurface: the trims in the basis'
            // parameters, scaled by its kind of parameter (angles, lengths, along a cone's
            // generatrix).
            let p = p.to_vec();
            let basis = self.surface(Self::reference(&p[1])?)?;
            let a = self.angle;
            let (u_factor, v_factor) = match &basis {
                Surface::Sphere(_) | Surface::Torus(_) => (a, a),
                Surface::Cylinder(_) => (a, l),
                Surface::Revolution(_) => (a, 1.0),
                Surface::Cone(c) => (a, l / c.semi_angle().cos()),
                Surface::Plane(_) => (l, l),
                _ => (1.0, 1.0),
            };
            let (u1, u2, v1, v2) = (self.real(&p[2])? * u_factor, self.real(&p[3])? * u_factor, self.real(&p[4])? * v_factor, self.real(&p[5])? * v_factor);
            let (u_sense, v_sense) = (p[6].as_bool().unwrap_or(true), p[7].as_bool().unwrap_or(true));
            return crate::geometry::surface::TrimmedSurface::with_senses(basis, u1, u2, v1, v2, u_sense, v_sense)
                .map(Surface::Trimmed)
                .ok_or_else(|| invalid(format!("#{id}: a trim outside its surface")));
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
        if ["B_SPLINE_SURFACE_WITH_KNOTS", "B_SPLINE_SURFACE", "UNIFORM_SURFACE", "QUASI_UNIFORM_SURFACE", "BEZIER_SURFACE"].iter().any(|t| inst.has(t)) {
            return self.bspline_surface(&inst);
        }
        Err(invalid(format!("#{id}: surfaces of type {} are not read yet", inst.name())))
    }

    fn bspline_surface(&self, inst: &Instance) -> Result<Surface> {
        // (Simple: name, u degree, v degree, poles, form, u closed, v closed, self_intersect
        // and, with knots, their multiplicities and values in u and v, and their spec; complex:
        // B_SPLINE_SURFACE(u degree, v degree, poles, ...) and B_SPLINE_SURFACE_WITH_KNOTS(...).
        // Uniform, quasi-uniform and Bézier surfaces leave their knots implicit.)
        let (ud, vd, grid, knots) = if !inst.is_complex() {
            let p = inst.params();
            (&p[1], &p[2], &p[3], inst.has("B_SPLINE_SURFACE_WITH_KNOTS").then(|| &p[8..]))
        } else {
            let b = inst.record("B_SPLINE_SURFACE").ok_or_else(|| invalid("a B-spline surface without its record".into()))?;
            (&b[0], &b[1], &b[2], inst.record("B_SPLINE_SURFACE_WITH_KNOTS"))
        };
        let int = |v: &Value| v.as_integer().map(|i| i as usize).ok_or_else(|| invalid("an integer".into()));
        let (ud, vd) = (int(ud)?, int(vd)?);
        let poles: Vec<Vec<Point>> = self.list(grid)?
            .iter()
            .map(|row| self.list(row)?.iter().map(|v| self.point(Self::reference(v)?)).collect())
            .collect::<Result<_>>()?;
        let ints = |v: &Value| -> Result<Vec<usize>> { self.list(v)?.iter().map(int).collect() };
        let reals = |v: &Value| -> Result<Vec<f64>> { self.list(v)?.iter().map(|v| self.real(v)).collect() };
        let (u, v) = match knots {
            Some(knots) => (KnotVector::new(reals(&knots[2])?, ints(&knots[0])?, ud), KnotVector::new(reals(&knots[3])?, ints(&knots[1])?, vd)),
            None => {
                let quasi = !inst.has("UNIFORM_SURFACE") || inst.has("QUASI_UNIFORM_SURFACE");
                let rows = poles.len();
                let columns = poles.first().map_or(0, Vec::len);
                let (um, uk) = implicit_knots(quasi, rows, ud);
                let (vm, vk) = implicit_knots(quasi, columns, vd);
                (KnotVector::new(uk, um, ud), KnotVector::new(vk, vm, vd))
            }
        };
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
        // (A face bounded by one vertex loop only, on a sphere, B-spline or revolution
        // surface: the whole surface within its natural bounds, as OCCT's
        // StepToTopoDS_TranslateFace makes it, ISO 10303-42 calling it the surface's domain.
        // On a torus OCCT leaves the face without bounds, which is the whole torus too.)
        let bounds = self.list(&p[1])?.to_vec();
        if let [only] = &bounds[..] {
            let b = self.get(Self::reference(only)?)?;
            let lp = b.record("FACE_OUTER_BOUND").or(b.record("FACE_BOUND")).and_then(|bp| bp.get(1)).and_then(Value::as_reference);
            let vertex_loop = lp.and_then(|l| self.get(l).ok()).is_some_and(|l| l.has("VERTEX_LOOP"));
            if vertex_loop && matches!(surface, Surface::Sphere(_) | Surface::BSpline(_) | Surface::Revolution(_) | Surface::Torus(_)) {
                let face = crate::make::face_on_surface(self.model, surface)?;
                self.face_instances.insert(face.node, id);
                return Ok(if same_sense { face } else { face.reversed() });
            }
        }
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
            let mut coedges = split_at_poles(self.model, &mut self.gaps, coedges, &surface, self.precision)?;
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
            let mut coedges = split_at_poles(self.model, &mut self.gaps, coedges, &surface, self.precision)?;
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
                // (A conic's parameter is an angle, in the file's angle unit.)
                let unit = if matches!(curve, Curve::Circle(_) | Curve::Ellipse(_)) { self.angle } else { scale };
                for item in self.list(v)? {
                    if let Some(("PARAMETER_VALUE", _)) = self.x.typed(item) {
                        return Ok(self.real(item)? * unit);
                    }
                }
                // (A bare number, as some files write it: a parameter value too, as OCCT
                // reads it.)
                for item in self.list(v)? {
                    if item.as_reference().is_none() && self.x.typed(item).is_none() {
                        if let Ok(t) = self.real(item) {
                            return Ok(t * unit);
                        }
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
            Some(&c) => self.curve2(c, made.scale, made.first, made.last).ok().map(|pc| within_bounds(pc, made, surface)),
            None => None,
        };
        Ok(match given {
            Some(pc) if fits(&pc, &made.span(), surface) => pc,
            // Parametrized otherwise but along the edge on a free-form surface (where
            // projecting is costly): kept, the tolerances covering the difference.
            Some(pc) if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) && follows(&pc, &made.span(), surface) => pc,
            given => match (project(&made.span(), surface, given.is_some(), self.precision), given) {
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
        self.shell_instances.insert(shell.node, id);
        Ok(shell)
    }

    // --- tessellated geometry (AP242) ----------------------------------------------------

    /// The triangulation of a triangulated face or surface set, plain or complex (`None` for
    /// other instances), as OCCT's `StepToTopoDS_TranslateFace` makes it
    /// (`CreatePolyTriangulation`): a node per index of the file's `pnindex` into its list of
    /// coordinates, or per point of that list where `pnindex` is empty; the triangles given
    /// as such, or those of the strips and fans (strips' and fans' triangles taken in OCCT's
    /// order and turn, those naming a point twice left out).
    fn triangulation(&self, id: usize) -> Result<Option<Triangulation>> {
        let inst = self.get(id)?;
        // (The record; whether it has a geometric link, which surface sets have not; whether
        // its triangles come as strips and fans.)
        let (p, linked, complex) = if let Some(p) = inst.record("TRIANGULATED_FACE") {
            (p, true, false)
        } else if let Some(p) = inst.record("COMPLEX_TRIANGULATED_FACE") {
            (p, true, true)
        } else if let Some(p) = inst.record("TRIANGULATED_SURFACE_SET") {
            (p, false, false)
        } else if let Some(p) = inst.record("COMPLEX_TRIANGULATED_SURFACE_SET") {
            (p, false, true)
        } else {
            return Ok(None);
        };
        let k = if linked { 5 } else { 4 };
        if p.len() < k + if complex { 3 } else { 2 } {
            return Err(invalid(format!("#{id}: too few parameters for a {}", inst.name())));
        }
        let list = Self::reference(&p[1])?;
        let coordinates = self.get(list)?.record("COORDINATES_LIST").ok_or_else(|| invalid(format!("#{list} is not a coordinates list")))?;
        let s = self.length;
        let points = self
            .list(coordinates.get(2).ok_or_else(|| invalid(format!("#{list}: no coordinates")))?)?
            .iter()
            .map(|c| match self.list(c)? {
                [x, y, z] => Ok(Point::new(self.real(x)? * s, self.real(y)? * s, self.real(z)? * s)),
                _ => Err(invalid(format!("#{list}: a point without three coordinates"))),
            })
            .collect::<Result<Vec<Point>>>()?;
        let int = |v: &Value| v.as_integer().ok_or_else(|| invalid(format!("#{id}: expected an integer, got {v:?}")));
        let ints = |v: &Value| -> Result<Vec<i64>> { self.x.list(v).unwrap_or(&[]).iter().map(int).collect() };
        // Each node's point in the list (0-based).
        let pnindex = ints(&p[k])?;
        let at: Vec<usize> = if pnindex.is_empty() {
            (0..points.len()).collect()
        } else {
            pnindex
                .iter()
                .map(|&i| if (1..=points.len() as i64).contains(&i) { Ok(i as usize - 1) } else { Err(invalid(format!("#{id} names point {i} of {}", points.len()))) })
                .collect::<Result<_>>()?
        };
        let mut triangles: Vec<[i64; 3]> = Vec::new();
        if !complex {
            for t in self.x.list(&p[k + 1]).unwrap_or(&[]) {
                match self.list(t)? {
                    [a, b, c] => triangles.push([int(a)?, int(b)?, int(c)?]),
                    _ => return Err(invalid(format!("#{id}: a triangle without three nodes"))),
                }
            }
        } else {
            // (OCCT's SetTriangles, 1-based as there: a strip's triangles at its odd
            // positions, then at its even ones; a fan's about its first node.)
            for strip in self.x.list(&p[k + 1]).unwrap_or(&[]) {
                let v = ints(strip)?;
                let at = |j: usize| v[j - 1];
                for j in (3..=v.len()).step_by(2) {
                    if at(j) != at(j - 2) && at(j) != at(j - 1) {
                        triangles.push([at(j - 2), at(j), at(j - 1)]);
                    }
                }
                for j in (4..=v.len()).step_by(2) {
                    if at(j) != at(j - 2) && at(j) != at(j - 1) {
                        triangles.push([at(j - 2), at(j - 1), at(j)]);
                    }
                }
            }
            for fan in self.x.list(&p[k + 2]).unwrap_or(&[]) {
                let v = ints(fan)?;
                let at = |j: usize| v[j - 1];
                for j in 3..=v.len() {
                    if at(j) != at(j - 2) && at(j - 1) != at(j - 2) {
                        triangles.push([at(1), at(j), at(j - 1)]);
                    }
                }
            }
        }
        let n = at.len() as i64;
        let triangles = triangles
            .into_iter()
            .map(|t| if t.iter().all(|i| (1..=n).contains(i)) { Ok(t.map(|i| (i - 1) as u32)) } else { Err(invalid(format!("#{id}: a triangle names a node beyond its {n}"))) })
            .collect::<Result<Vec<[u32; 3]>>>()?;
        // (OCCT's SetNormals: one for every node, or one per node.)
        let normals = self
            .x
            .list(&p[3])
            .unwrap_or(&[])
            .iter()
            .map(|v| match self.list(v)? {
                [x, y, z] => Ok(Direction::try_from_vector(Vector::new(self.real(x)?, self.real(y)?, self.real(z)?))),
                _ => Ok(None),
            })
            .collect::<Result<Option<Vec<Direction>>>>()?
            .and_then(|d| match d.len() {
                1 => Some(vec![d[0]; at.len()]),
                m if m == at.len() && m > 0 => Some(d),
                _ => None,
            });
        let link = if linked { p[4].as_reference() } else { None };
        Ok(Some(Triangulation { list, nodes: at.iter().map(|&i| points[i]).collect(), at, triangles, normals, link }))
    }

    /// The planar faces of triangulations, in a compound (OCCT's
    /// `BRepBuilderAPI_MakeShapeOnMesh` on one mesh of all of them): nodes at the same
    /// point of the same coordinates list one vertex, so that the faces of one solid or
    /// shell share their edges across the file's faces. OCCT reads a triangulated face as a
    /// face carrying its triangulation and no surface; kentron's measures, classification
    /// and checks work on faces' surfaces, so the triangles are read as the planes they span.
    fn mesh_faces(&mut self, triangulations: &[Triangulation]) -> Result<Shape> {
        let mut index: HashMap<(usize, usize), u32> = HashMap::default();
        let (mut nodes, mut triangles) = (Vec::new(), Vec::new());
        for t in triangulations {
            let mut node = Vec::with_capacity(t.at.len());
            for (k, &i) in t.at.iter().enumerate() {
                node.push(*index.entry((t.list, i)).or_insert_with(|| {
                    nodes.push(t.nodes[k]);
                    nodes.len() as u32 - 1
                }));
            }
            triangles.extend(t.triangles.iter().map(|tri| tri.map(|i| node[i as usize])));
        }
        crate::make::shape_on_mesh(self.model, &crate::mesh::TriangleMesh::new(nodes, triangles))
    }

    /// A tessellated solid or shell, or a triangulated face or surface set standing alone
    /// (`None` for other instances), as OCCT's `StepToTopoDS_TranslateSolid`,
    /// `TranslateShell` and `TranslateFace` read them (with `read.step.tessellated` On, its
    /// default): a solid or shell already read from the B-rep the item is linked to is that
    /// solid or shell, its faces linked to faces read given their triangulations; otherwise
    /// a new solid of one shell, or a new shell, of the items' faces ([`Self::mesh_faces`]).
    /// (A face linked to a B-rep face read, in a solid or shell not linked, is read from its
    /// triangles too, where OCCT would put the B-rep face in the new shell.)
    fn tessellated(&mut self, id: usize) -> Result<Option<Shape>> {
        let inst = self.get(id)?;
        if let Some(t) = self.triangulation(id)? {
            return Ok(Some(self.mesh_faces(&[t])?));
        }
        let (p, solid) = match (inst.record("TESSELLATED_SOLID"), inst.record("TESSELLATED_SHELL")) {
            (Some(p), _) => (p.to_vec(), true),
            (None, Some(p)) => (p.to_vec(), false),
            _ => return Ok(None),
        };
        let mut triangulations = Vec::new();
        for item in self.x.list(&p[1]).unwrap_or(&[]) {
            let item = Self::reference(item)?;
            // (Only tessellated faces: a surface set is no item of a solid or shell.)
            let face = self.get(item)?;
            if face.has("TRIANGULATED_FACE") || face.has("COMPLEX_TRIANGULATED_FACE") {
                triangulations.extend(self.triangulation(item)?);
            }
        }
        let link = p.get(2).and_then(Value::as_reference);
        let linked = link.and_then(|l| {
            if solid {
                self.bodies.get(&l).copied()
            } else {
                self.shell_instances.iter().find(|&(_, &i)| i == l).map(|(&node, _)| Shape { node, orientation: crate::topology::Orientation::Forward, placement: Transform::IDENTITY })
            }
        });
        if let Some(linked) = linked {
            for t in triangulations {
                let Some(face_link) = t.link else { continue };
                let Some(node) = self.face_instances.iter().find(|&(_, &i)| i == face_link).map(|(&node, _)| node) else { continue };
                let mesh = crate::mesh::TriangleMesh::new(t.nodes, t.triangles);
                let mesh = match t.normals {
                    Some(normals) => mesh.with_normals(normals),
                    None => mesh,
                };
                self.model.set_face_mesh(Shape { node, orientation: crate::topology::Orientation::Forward, placement: Transform::IDENTITY }, Some(mesh))?;
            }
            return Ok(Some(linked));
        }
        let faces = self.mesh_faces(&triangulations)?;
        let shell = self.model.new_empty(NodeKind::Shell);
        for face in self.model.children(faces).collect::<Vec<_>>() {
            self.model.add(shell, face)?;
        }
        if !solid {
            self.shell_instances.insert(shell.node, id);
            return Ok(Some(shell));
        }
        let made = self.model.new_empty(NodeKind::Solid);
        self.model.add(made, shell)?;
        Ok(Some(made))
    }

    /// A solid or shell of a representation (`None` for other items, such as placements).
    fn body(&mut self, id: usize) -> Result<Option<Shape>> {
        if let Some(&s) = self.bodies.get(&id) {
            return Ok(Some(s));
        }
        if let Some(made) = self.tessellated(id)? {
            self.bodies.insert(id, made);
            return Ok(Some(made));
        }
        let inst = self.get(id)?;
        let made = if let Some(p) = inst.record("MANIFOLD_SOLID_BREP").or(inst.record("FACETED_BREP")).or(inst.record("BREP_WITH_VOIDS").filter(|_| !inst.is_complex())) {
            let p = p.to_vec();
            let solid = self.model.new_empty(NodeKind::Solid);
            let outer = self.shell(Self::reference(&p[1])?)?;
            self.model.add(solid, outer)?;
            // (Its voids: the third parameter of a simple BREP_WITH_VOIDS, the only one of the
            // complex entity's part.)
            let voids = if inst.is_complex() { inst.record("BREP_WITH_VOIDS").and_then(|v| v.first().cloned()) } else { p.get(2).cloned() };
            if let Some(voids) = voids {
                for v in self.list(&voids)?.to_vec() {
                    let void = self.shell(Self::reference(&v)?)?;
                    self.model.add(solid, void)?;
                }
            }
            solid
        } else if let Some(p) = inst.record("EDGE_BASED_WIREFRAME_MODEL") {
            // Its connected edge sets' edges, as a compound.
            let c = self.model.new_empty(NodeKind::Compound);
            for set in self.list(&p[1])?.to_vec() {
                let set = self.get(Self::reference(&set)?)?;
                let Some(edges) = set.record("CONNECTED_EDGE_SET").map(|e| e.to_vec()) else { continue };
                for e in self.list(&edges[1])?.to_vec() {
                    let made = self.edge(Self::reference(&e)?)?;
                    self.model.add(c, Shape::from(made.id))?;
                }
            }
            c
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
            // As OCCT's StepToTopoDS_Builder::Init for a geometric set, element by element:
            // a composite curve as its segments' edges, another curve as an edge over its
            // range (an unbounded one left out), a point as a vertex, a curve-bounded surface
            // as its face, another bounded one as a face over its bounds, any other item as
            // the shape it reads as.
            let c = self.model.new_empty(NodeKind::Compound);
            let mut any = false;
            for item in self.list(&p[1])?.to_vec() {
                let item = Self::reference(&item)?;
                let it = self.get(item)?;
                let element = if let Some(cp) = it.record("COMPOSITE_CURVE").map(<[Value]>::to_vec) {
                    // (As StepToTopoDS_TranslateCompositeCurve: a wire of its segments, each
                    // the way it runs; their edges loose where they do not join.)
                    let mut edges = Vec::new();
                    for segment in self.list(&cp[1])?.to_vec() {
                        let sp = self.get(Self::reference(&segment)?)?.record("COMPOSITE_CURVE_SEGMENT").ok_or_else(|| invalid("a composite curve's segment".into()))?.to_vec();
                        let made = self.segment_edge(Self::reference(&sp[2])?)?;
                        let edge = Shape::from(made.id);
                        edges.push(if sp[1].as_bool().unwrap_or(true) { edge } else { edge.reversed() });
                    }
                    match crate::make::wire(self.model, &edges) {
                        Ok(wire) => {
                            self.item_instances.insert(wire.node, item);
                            Some(wire)
                        }
                        Err(_) => {
                            let loose = self.model.new_empty(NodeKind::Compound);
                            for e in edges {
                                self.model.add(loose, e)?;
                            }
                            Some(loose)
                        }
                    }
                } else if (is_curve(&it) || it.has("CARTESIAN_POINT")) && !self.three_d(item) {
                    None
                } else if is_curve(&it) {
                    match self.curve(item).and_then(|(curve, scale)| Ok((self.trim_range(item, &curve, scale)?, curve))) {
                        // (An unbounded curve (Precision::IsInfinite) left out: OCCT makes an
                        // infinite edge of it, kentron none.)
                        Ok(((a, b), _)) if a.abs() >= 1e100 || b.abs() >= 1e100 => None,
                        Ok(((a, b), curve)) => match crate::make::edge_from_curve_range(self.model, curve, a, b) {
                            Ok(edge) => {
                                // (The curve's styles and layers are the edge's.)
                                self.edge_instances.insert(edge.node, item);
                                Some(edge)
                            }
                            Err(_) => None,
                        },
                        Err(_) => None,
                    }
                } else if it.has("CARTESIAN_POINT") {
                    let vertex = crate::make::vertex(self.model, self.point(item)?);
                    self.item_instances.insert(vertex.node, item);
                    Some(vertex)
                } else if it.has("CURVE_BOUNDED_SURFACE") {
                    let face = self.curve_bounded_face(item)?;
                    // (The surface's styles and layers are the face's.)
                    self.face_instances.insert(face.node, item);
                    Some(face)
                } else if it.has("B_SPLINE_SURFACE") || it.has("B_SPLINE_SURFACE_WITH_KNOTS") || it.has("RECTANGULAR_TRIMMED_SURFACE") || it.has("BEZIER_SURFACE") {
                    let surface = self.surface(item)?;
                    let face = crate::make::face_on_surface(self.model, surface)?;
                    self.face_instances.insert(face.node, item);
                    Some(face)
                } else {
                    self.body(item)?
                };
                if let Some(element) = element {
                    self.model.add(c, element)?;
                    any = true;
                }
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
            // (As OCCT's ReadNames: the occurrence's description, else its name, else its id.)
            let text = |k: usize| p.get(k).and_then(Value::as_str).filter(|s| !s.trim().is_empty());
            let name = text(2).or(text(1)).or(p.first().and_then(Value::as_str)).unwrap_or("").to_string();
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
        // (With its subtype AP203 writes, as OCCT's StepBasic_ProductDefinitionWithAssociated-
        // Documents is a product definition.)
        let mut definitions: Vec<usize> = x.of_type("PRODUCT_DEFINITION").chain(x.of_type("PRODUCT_DEFINITION_WITH_ASSOCIATED_DOCUMENTS")).map(|(id, _)| id).collect();
        definitions.sort_unstable();
        definitions.dedup();
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
        // (As OCCT's STEPCAFControl_Reader::ReadNames: the product's name, else its id.)
        let text = |k: usize| product.get(k).and_then(Value::as_str).filter(|s| !s.trim().is_empty());
        Ok(text(1).or(product.first().and_then(Value::as_str)).unwrap_or("").to_string())
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
        // Each representation's bodies (each body once), with the items they are made of.
        let mut groups: Vec<(usize, Vec<(usize, Shape)>)> = Vec::new();
        let mut seen: crate::hashing::FastSet<crate::topology::NodeId> = Default::default();
        for rep in reps {
            let inst = self.get(rep)?;
            let Some(items) = inst.records().find(|(_, v)| v.len() >= 2 && self.x.list(&v[1]).is_some()).map(|(_, v)| v[1].clone()) else { continue };
            let mut group = Vec::new();
            for item in self.list(&items)? {
                let item = Self::reference(item)?;
                if let Some(body) = self.body(item)? {
                    if seen.insert(body.node) {
                        group.push((item, body));
                    }
                }
            }
            if !group.is_empty() {
                groups.push((rep, group));
            }
        }
        let material = self.materials.get(&definition).cloned();
        let mut kids = Vec::new();
        // A shape spread over several representations linked to the product's (as a surface
        // model and a wireframe): an assembly of them, as OCCT reads it, each a part of its
        // own (named as the representation is, else as the product).
        let (shape, color, alpha, layer, subshapes) = if groups.len() > 1 && !children.contains_key(&definition) {
            for (rep, group) in &groups {
                let (shape, color, alpha, layer, subshapes) = self.attributes(group)?;
                let rep_name = self.get(*rep)?.records().find_map(|(_, v)| v.first().and_then(Value::as_str).map(str::to_string)).filter(|n| !n.is_empty());
                kids.push(Node { name: rep_name.unwrap_or_else(|| name.clone()), shape, color, alpha, placement: Transform::IDENTITY, children: Vec::new(), subshapes, material: None, layer, definition: None, ..Default::default() });
            }
            (None, None, None, None, Vec::new())
        } else {
            let all: Vec<(usize, Shape)> = groups.into_iter().flat_map(|(_, g)| g).collect();
            self.attributes(&all)?
        };
        for (child, occurrence, occurrence_name) in children.get(&definition).cloned().unwrap_or_default() {
            let child_name = if occurrence_name.is_empty() { self.product_name(child)? } else { occurrence_name };
            let at = placement_of.get(&occurrence).copied().unwrap_or(Transform::IDENTITY);
            let mut kid = self.node(child, child_name, at, reps_of, linked, children, placement_of, depth + 1)?;
            kid.occurrence_validation = self.validations.get(&occurrence).copied().unwrap_or_default();
            kid.occurrence_color = self.occurrence_colors.get(&occurrence).copied();
            kid.occurrence_layer = self.occurrence_layers.get(&occurrence).cloned();
            kids.push(kid);
        }
        let validation = self.validations.get(&definition).copied().unwrap_or_default();
        // (An external reference only for a product with a shape representation, as OCCT
        // binds one only to a definition it transfers, and whose shape is a compound (an
        // assembly's, or one of several bodies), as STEPCAFControl_Reader::AddShape records
        // one only for a compound.)
        let compound = shape.is_none_or(|s| s.kind() == NodeKind::Compound);
        let external = self.externals.get(&definition).filter(|_| compound && reps_of.contains_key(&definition)).cloned();
        Ok(Node { name, shape, color, alpha, placement, children: kids, subshapes, material, layer, validation, occurrence_validation: Validation::default(), definition: Some(definition), external, occurrence_color: None, occurrence_layer: None })
    }

    /// The shape of `bodies` (each with the item it is made of: a compound of several) and
    /// what the file gives them: the colour they share (and its opacity), the layer of the
    /// shape as a whole, and the parts with names, colours or layers of their own (each body
    /// of several, the shells of bodies (a surface model's), faces and edges).
    #[allow(clippy::type_complexity)]
    fn attributes(&mut self, bodies: &[(usize, Shape)]) -> Result<(Option<Shape>, Option<[f64; 3]>, Option<f64>, Option<String>, Vec<Subshape>)> {
        // (A body that is a shell, as a surface model of one shell reads: what the file gives
        // the shell is the body's, as OCCT puts a styled shape's colour on the shape's label.)
        let shell_of = |body: Shape| self.shell_instances.get(&body.node).copied();
        let body_color = |item: usize, body: Shape| self.colors.get(&item).or_else(|| shell_of(body).and_then(|i| self.colors.get(&i))).copied();
        let body_layer = |item: usize, body: Shape| self.layers.get(&item).or_else(|| shell_of(body).and_then(|i| self.layers.get(&i))).cloned();
        let colors: Vec<Option<[f64; 4]>> = bodies.iter().map(|&(item, body)| body_color(item, body)).collect();
        let shared = match colors.first() {
            Some(&Some(c)) if colors.iter().all(|&k| k == Some(c)) => Some(c),
            _ => None,
        };
        let (color, alpha) = (shared.map(|[r, g, b, _]| [r, g, b]), shared.map(|c| c[3]).filter(|&a| a < 1.0));
        let layer = match bodies {
            [(item, body)] => body_layer(*item, *body),
            _ => None,
        };
        let forward = |s: Shape| Shape { orientation: crate::topology::Orientation::Forward, ..s };
        let mut subshapes = Vec::new();
        if bodies.len() > 1 {
            for &(item, body) in bodies {
                let color = body_color(item, body);
                let layer = body_layer(item, body);
                if color.is_some() || layer.is_some() {
                    subshapes.push(Subshape { shape: forward(body), name: None, color, layer });
                }
            }
        }
        for &(_, body) in bodies {
            for s in self.model.unique_subshapes(body, NodeKind::Shell) {
                if bodies.iter().any(|(_, b)| b.node == s.node) {
                    continue;
                }
                let Some(&inst) = self.shell_instances.get(&s.node) else { continue };
                let color = self.colors.get(&inst).copied();
                let layer = self.layers.get(&inst).cloned();
                if color.is_some() || layer.is_some() {
                    subshapes.push(Subshape { shape: forward(s), name: None, color, layer });
                }
            }
        }
        for &(_, body) in bodies {
            for (kind, instances) in [(NodeKind::Face, &self.face_instances), (NodeKind::Wire, &self.item_instances), (NodeKind::Edge, &self.edge_instances), (NodeKind::Vertex, &self.item_instances)] {
                for s in self.model.unique_subshapes(body, kind) {
                    let Some(&inst) = instances.get(&s.node) else { continue };
                    let name = self.x.get(inst).ok().and_then(|i| i.params().first().and_then(Value::as_str).map(str::to_string)).filter(|n| !n.is_empty() && n != "NONE");
                    let color = self.colors.get(&inst).copied();
                    let layer = self.layers.get(&inst).cloned();
                    if name.is_some() || color.is_some() || layer.is_some() {
                        subshapes.push(Subshape { shape: forward(s), name, color, layer });
                    }
                }
            }
        }
        let shape = self.compound_of(bodies.iter().map(|&(_, b)| b).collect())?;
        Ok((shape, color, alpha, layer, subshapes))
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

/// Whether `inst` is a B-spline curve: with knots, or of one of the types leaving them
/// implicit (uniform, quasi-uniform, Bézier).
/// Whether an instance is a curve `Reader::curve` reads.
fn is_curve(inst: &Instance) -> bool {
    ["LINE", "CIRCLE", "ELLIPSE", "HYPERBOLA", "PARABOLA", "TRIMMED_CURVE", "SURFACE_CURVE", "SEAM_CURVE", "INTERSECTION_CURVE"].iter().any(|t| inst.has(t)) || is_bspline_curve(inst)
}

fn is_bspline_curve(inst: &Instance) -> bool {
    ["B_SPLINE_CURVE_WITH_KNOTS", "B_SPLINE_CURVE", "UNIFORM_CURVE", "QUASI_UNIFORM_CURVE", "BEZIER_CURVE"].iter().any(|t| inst.has(t))
}

/// Knots of the uniform and quasi-uniform B-spline types (which leave them implicit; a
/// Bézier curve's are a quasi-uniform one's).
fn uniform_knots(inst: &Instance, poles: usize, degree: usize) -> (Vec<usize>, Vec<f64>) {
    implicit_knots(!inst.has("UNIFORM_CURVE") || inst.has("QUASI_UNIFORM_CURVE"), poles, degree)
}

/// The multiplicities and knots of a quasi-uniform B-spline (`quasi`: clamped at both ends,
/// uniform between) or a uniform one, of `poles` poles and degree `degree`.
fn implicit_knots(quasi: bool, poles: usize, degree: usize) -> (Vec<usize>, Vec<f64>) {
    if quasi {
        let spans = poles.saturating_sub(degree).max(1);
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


/// The file's distance uncertainty in millimetres, OCCT's reading precision
/// (`STEPControl_ActorRead::PrepareUnits`, `read.precision.mode` "File"): the global
/// uncertainty's length measure, else `read.precision.val` (1e-4).
fn uncertainty(x: &Exchange) -> f64 {
    for (_, inst) in x.of_type("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT") {
        let Some(list) = inst.record("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT").and_then(|p| x.list(&p[0])) else { continue };
        for u in list.iter().filter_map(Value::as_reference) {
            let Ok(m) = x.get(u) else { continue };
            let Some(mw) = m.record("UNCERTAINTY_MEASURE_WITH_UNIT") else { continue };
            let Some(value) = x.real(&mw[0]) else { continue };
            let factor = mw[1].as_reference().and_then(|r| x.get(r).ok()).map_or(1.0, |unit| unit_factor(x, &unit, true));
            if value > 0.0 {
                return value * factor;
            }
        }
    }
    1e-4
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

/// The colours of styled items: item instance -> RGB and opacity. As OCCT's
/// STEPConstruct_Styles::GetColors reads a style and STEPCAFControl_Reader applies it: of
/// each presentation style, a surface style usage's fill colour, its rendering's colour (and
/// transparency) over it, and its boundary's curve colour; a curve style's colour. kentron
/// keeps one colour an item: the surface's, else the boundary's or curve's (parameter-line
/// styles and the like give none).
fn colors(x: &Exchange) -> HashMap<usize, [f64; 4]> {
    let mut out = HashMap::default();
    // (Over-riding styled items too: the same, with the style they override after. Styles
    // by context are an occurrence's: `occurrence_colors`.)
    for (_, inst) in x.of_type("STYLED_ITEM").chain(x.of_type("OVER_RIDING_STYLED_ITEM")) {
        let Some(p) = inst.record("STYLED_ITEM").or(inst.record("OVER_RIDING_STYLED_ITEM")) else { continue };
        let Some(item) = p.get(2).and_then(Value::as_reference) else { continue };
        let assignments: Vec<Value> = x.list(&p[1]).map(<[Value]>::to_vec).unwrap_or_default();
        let plain: Vec<Value> = assignments.into_iter().filter(|a| a.as_reference().and_then(|r| x.get(r).ok()).is_some_and(|i| i.has("PRESENTATION_STYLE_ASSIGNMENT") && !i.has("PRESENTATION_STYLE_BY_CONTEXT"))).collect();
        if let Some(c) = style_colour(x, &plain) {
            out.entry(item).or_insert(c);
        }
    }
    out
}

/// The mapped items of occurrences' own shape representations: item -> the occurrence (a
/// next assembly usage) whose product definition shape the representation is the shape of.
fn occurrence_items(x: &Exchange) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for (_, inst) in x.of_type("SHAPE_DEFINITION_REPRESENTATION") {
        let p = inst.params();
        let (Some(pds), Some(rep)) = (p.first().and_then(Value::as_reference), p.get(1).and_then(Value::as_reference)) else { continue };
        let Some(definition) = x.get(pds).ok().and_then(|i| i.record("PRODUCT_DEFINITION_SHAPE").and_then(|d| d.get(2).and_then(Value::as_reference))) else { continue };
        if !x.get(definition).is_ok_and(|d| d.has("NEXT_ASSEMBLY_USAGE_OCCURRENCE")) {
            continue;
        }
        let Some(items) = x.get(rep).ok().and_then(|r| r.records().find_map(|(_, v)| v.get(1).and_then(|l| x.list(l)).map(<[Value]>::to_vec))) else { continue };
        for item in items.iter().filter_map(Value::as_reference) {
            if x.get(item).is_ok_and(|i| i.has("MAPPED_ITEM")) {
                out.push((item, definition));
            }
        }
    }
    out
}

/// The colours styles by context give occurrences: the occurrence (a next assembly usage)
/// whose product definition shape the context's shape representation is the shape of ->
/// RGB and opacity, as OCCT reads a component's colour.
fn occurrence_colors(x: &Exchange) -> HashMap<usize, [f64; 4]> {
    let mut out = HashMap::default();
    // Shape representation -> the occurrence it is a shape of.
    let mut occurrence_of: HashMap<usize, usize> = HashMap::default();
    for (_, inst) in x.of_type("SHAPE_DEFINITION_REPRESENTATION") {
        let p = inst.params();
        let (Some(pds), Some(rep)) = (p.first().and_then(Value::as_reference), p.get(1).and_then(Value::as_reference)) else { continue };
        let Some(definition) = x.get(pds).ok().and_then(|i| i.record("PRODUCT_DEFINITION_SHAPE").and_then(|d| d.get(2).and_then(Value::as_reference))) else { continue };
        if x.get(definition).is_ok_and(|d| d.has("NEXT_ASSEMBLY_USAGE_OCCURRENCE")) {
            occurrence_of.insert(rep, definition);
        }
    }
    if occurrence_of.is_empty() {
        return out;
    }
    for (_, inst) in x.of_type("STYLED_ITEM").chain(x.of_type("OVER_RIDING_STYLED_ITEM")) {
        let Some(p) = inst.record("STYLED_ITEM").or(inst.record("OVER_RIDING_STYLED_ITEM")) else { continue };
        for a in x.list(&p[1]).unwrap_or(&[]) {
            let Some(psbc) = a.as_reference().and_then(|r| x.get(r).ok()).and_then(|i| i.record("PRESENTATION_STYLE_BY_CONTEXT").map(<[Value]>::to_vec)) else { continue };
            let Some(&occurrence) = psbc.get(1).and_then(Value::as_reference).and_then(|c| occurrence_of.get(&c)) else { continue };
            if let Some(c) = style_colour(x, std::slice::from_ref(a)) {
                out.entry(occurrence).or_insert(c);
            }
        }
    }
    out
}

/// The colour a styled item's presentation style assignments give, as OCCT's
/// STEPConstruct_Styles::GetColors reads them and STEPCAFControl_Reader applies them: of
/// each, a surface style usage's fill colour, its rendering's colour (and transparency) over
/// it, and its boundary's curve colour; a curve style's colour. kentron keeps one colour an
/// item: the surface's, else the boundary's or curve's (parameter-line styles and the like
/// give none).
fn style_colour(x: &Exchange, assignments: &[Value]) -> Option<[f64; 4]> {
    let get = |v: &Value| v.as_reference().and_then(|r| x.get(r).ok());
    let list = |v: &Value| -> Vec<Value> { x.list(v).map(<[Value]>::to_vec).unwrap_or_default() };
    // A colour entity's components.
    let colour = |v: &Value| -> Option<[f64; 3]> {
        let c = get(v)?;
        if let Some(c) = c.record("COLOUR_RGB") {
            return Some([x.real(&c[1])?, x.real(&c[2])?, x.real(&c[3])?]);
        }
        if let Some(c) = c.record("DRAUGHTING_PRE_DEFINED_COLOUR") {
            return c[0].as_str().and_then(named_colour);
        }
        None
    };
    let (mut surface, mut boundary, mut curve, mut rendering) = (None, None, None, None);
    let mut transparency: Option<f64> = None;
    for psa in assignments {
        let Some(psa) = get(psa) else { continue };
        let Some(styles) = psa.record("PRESENTATION_STYLE_ASSIGNMENT").or(psa.record("PRESENTATION_STYLE_BY_CONTEXT")).and_then(|r| r.first().cloned()) else { continue };
        for select in list(&styles) {
            let Some(style) = get(&select) else { continue };
            if let Some(usage) = style.record("SURFACE_STYLE_USAGE") {
                let Some(side) = usage.get(1).and_then(get) else { continue };
                let Some(elements) = side.record("SURFACE_SIDE_STYLE").and_then(|r| r.get(1).cloned()) else { continue };
                for element in list(&elements) {
                    let Some(e) = get(&element) else { continue };
                    if let Some(fill) = e.record("SURFACE_STYLE_FILL_AREA") {
                        // FILL_AREA_STYLE(name, (FILL_AREA_STYLE_COLOUR(name, colour), ...))
                        let fas = fill.first().and_then(get);
                        let found = fas.and_then(|f| f.record("FILL_AREA_STYLE").and_then(|r| r.get(1).cloned())).map(|fs| {
                            list(&fs).iter().filter_map(get).find_map(|c| c.record("FILL_AREA_STYLE_COLOUR").and_then(|r| r.get(1).and_then(colour)))
                        });
                        surface = surface.or(found.flatten());
                    } else if let Some(b) = e.record("SURFACE_STYLE_BOUNDARY") {
                        boundary = boundary.or(b.first().and_then(get).and_then(|cs| cs.record("CURVE_STYLE").and_then(|r| r.get(3).and_then(colour))));
                    } else if let Some(r) = e.record("SURFACE_STYLE_RENDERING_WITH_PROPERTIES").or(e.record("SURFACE_STYLE_RENDERING")) {
                        rendering = rendering.or(r.get(1).and_then(colour));
                        if let Some(props) = r.get(2) {
                            for prop in list(props) {
                                if let Some(t) = get(&prop).and_then(|t| t.record("SURFACE_STYLE_TRANSPARENT").and_then(|v| v.first().and_then(|v| x.real(v)))) {
                                    transparency = transparency.or(Some(t));
                                }
                            }
                        }
                    }
                }
            } else if let Some(cs) = style.record("CURVE_STYLE") {
                curve = curve.or(cs.get(3).and_then(colour));
            }
        }
    }
    let [r, g, b] = rendering.or(surface).or(boundary).or(curve)?;
    Some([r, g, b, 1.0 - transparency.unwrap_or(0.0).clamp(0.0, 1.0)])
}

/// The external references of product definitions: definition -> the referenced file's name,
/// as OCCT's STEPConstruct_ExternRefs::LoadExternRefs and FileName find them (an applied
/// document reference's, AP214; else an AP203 product definition with associated documents'),
/// those of a format other than STEP AP203 or AP214 left out, as STEPCAFControl_Reader does.
/// (The file is not looked for: kentron reads from bytes. OCCT, not finding it, keeps the
/// name the external identification gives, else the document file's.)
fn external_references(x: &Exchange) -> HashMap<usize, String> {
    let mut out: HashMap<usize, String> = HashMap::default();
    let adrs: Vec<usize> = x.of_type("APPLIED_DOCUMENT_REFERENCE").map(|(id, _)| id).collect();
    let mut pdwads: Vec<usize> = x.of_type("PRODUCT_DEFINITION_WITH_ASSOCIATED_DOCUMENTS").map(|(id, _)| id).collect();
    if adrs.is_empty() && pdwads.is_empty() {
        return out;
    }
    // (Graph().Shareds: the instances an instance refers to; Sharings: those referring to it.)
    let shareds = |id: usize| -> Vec<usize> {
        let mut found = Vec::new();
        let Ok(inst) = x.get(id) else { return found };
        let mut values: Vec<&Value> = inst.records().flat_map(|(_, v)| v.iter()).collect();
        while let Some(v) = values.pop() {
            match v {
                Value::Reference(r) => found.push(*r),
                Value::List(_) => values.extend(x.list(v).unwrap_or(&[])),
                Value::Typed(_) => values.extend(x.typed(v).map_or(&[][..], |(_, p)| p)),
                _ => {}
            }
        }
        found.reverse();
        found
    };
    let mut sharings: HashMap<usize, Vec<usize>> = HashMap::default();
    for id in x.ids() {
        for r in shareds(id) {
            sharings.entry(r).or_default().push(id);
        }
    }
    let sharing = |id: usize| sharings.get(&id).cloned().unwrap_or_default();
    let is = |id: usize, t: &str| x.get(id).is_ok_and(|i| i.has(t));
    let is_definition = |id: usize| is(id, "PRODUCT_DEFINITION") || is(id, "PRODUCT_DEFINITION_WITH_ASSOCIATED_DOCUMENTS");
    let is_formation = |id: usize| is(id, "PRODUCT_DEFINITION_FORMATION") || is(id, "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE");
    let text = |id: usize, k: usize| x.get(id).ok().and_then(|i| i.params().get(k).and_then(Value::as_str).map(str::to_string)).filter(|s| !s.is_empty());
    for adr in adrs {
        // findPDWADandExcludeExcess: the definition with documents of the product the
        // referenced document is equivalent to, taken off the AP203 ones.
        let mut pdwad = None;
        for doc in shareds(adr).into_iter().filter(|&d| is(d, "DOCUMENT") || is(d, "DOCUMENT_FILE")) {
            for dpe in sharing(doc).into_iter().filter(|&e| is(e, "DOCUMENT_PRODUCT_EQUIVALENCE")) {
                for pdf in shareds(dpe).into_iter().filter(|&f| is_formation(f)) {
                    for s in sharing(pdf).into_iter().filter(|&s| is(s, "PRODUCT_DEFINITION_WITH_ASSOCIATED_DOCUMENTS")) {
                        pdwad = Some(s);
                    }
                    if let Some(p) = pdwad {
                        pdwads.retain(|&q| q != p);
                    }
                }
            }
        }
        let shape = shareds(adr).into_iter().filter(|&d| is_definition(d)).last();
        // The document file: the definition's (with documents), else the reference's; with the
        // format a shape representation of it describes.
        let doc_file = shareds(pdwad.unwrap_or(adr)).into_iter().filter(|&d| is(d, "DOCUMENT_FILE")).last();
        let Some(doc_file) = doc_file else { continue };
        let mut format = None;
        for pd in sharing(doc_file).into_iter().filter(|&p| is(p, "PROPERTY_DEFINITION")) {
            for pdr in sharing(pd).into_iter().filter(|&r| is(r, "PROPERTY_DEFINITION_REPRESENTATION")) {
                let Some(rep) = x.get(pdr).ok().and_then(|i| i.params().get(1).and_then(Value::as_reference)) else { continue };
                if !(is(rep, "SHAPE_REPRESENTATION") || x.get(rep).is_ok_and(|i| i.name().ends_with("SHAPE_REPRESENTATION"))) {
                    continue;
                }
                let items = x.get(rep).ok().and_then(|i| i.params().get(1).and_then(|v| x.list(v)).map(<[Value]>::to_vec)).unwrap_or_default();
                format = items.iter().filter_map(Value::as_reference).find(|&i| is(i, "DESCRIPTIVE_REPRESENTATION_ITEM")).map(|i| text(i, 1).unwrap_or_default());
            }
        }
        if let Some(f) = &format {
            if !f.eq_ignore_ascii_case("STEP AP203") && !f.eq_ignore_ascii_case("STEP AP214") {
                continue;
            }
        }
        // FileName: the external identification's assigned id (not an entity number), or
        // its source's; else the document file's id, else its name.
        let aeia = sharing(doc_file).into_iter().find(|&a| is(a, "APPLIED_EXTERNAL_IDENTIFICATION_ASSIGNMENT"));
        let from_aeia = aeia.and_then(|a| {
            text(a, 0).filter(|n| !n.starts_with('#')).or_else(|| {
                let source = x.get(a).ok()?.params().get(2)?.as_reference()?;
                let id = x.get(source).ok()?.params().first()?.clone();
                x.typed(&id).and_then(|(_, p)| p.first().and_then(Value::as_str)).or(id.as_str()).map(str::to_string).filter(|s| !s.is_empty())
            })
        });
        let name = from_aeia.or_else(|| text(doc_file, 0)).or_else(|| text(doc_file, 1));
        if let (Some(shape), Some(name)) = (shape, name) {
            out.insert(shape, name);
        }
    }
    // The AP203 ones left: the first of the definition's documents with a name.
    for pdwad in pdwads {
        let docs = x.get(pdwad).ok().and_then(|i| i.params().get(4).and_then(|v| x.list(v)).map(<[Value]>::to_vec)).unwrap_or_default();
        if let Some(name) = docs.iter().filter_map(Value::as_reference).find_map(|d| text(d, 1)) {
            out.insert(pdwad, name);
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
/// The geometric validation properties: by the product definition (or occurrence) they are
/// given to (directly, or through a shape aspect of its shape), with lengths in millimetres
/// (`length` per unit).
fn validations(x: &Exchange, length: f64) -> HashMap<usize, Validation> {
    let mut out: HashMap<usize, Validation> = HashMap::default();
    // (The definition or occurrence a product definition shape, or a shape aspect of one, is
    // of.)
    let target = |id: usize| -> Option<usize> {
        let i = x.get(id).ok()?;
        let pds = match i.record("SHAPE_ASPECT") {
            Some(a) => x.get(a.get(2)?.as_reference()?).ok()?,
            None => i,
        };
        pds.record("PRODUCT_DEFINITION_SHAPE")?.get(2)?.as_reference()
    };
    for (_, inst) in x.of_type("PROPERTY_DEFINITION_REPRESENTATION") {
        let p = inst.params();
        let (Some(prop), Some(rep)) = (p.first().and_then(Value::as_reference), p.get(1).and_then(Value::as_reference)) else { continue };
        let Ok(prop) = x.get(prop) else { continue };
        let Some(pp) = prop.record("PROPERTY_DEFINITION") else { continue };
        let name = pp.first().and_then(Value::as_str).unwrap_or("").to_ascii_lowercase().replace('_', " ");
        if !name.contains("validation") {
            continue;
        }
        let Some(of) = pp.get(2).and_then(Value::as_reference).and_then(target) else { continue };
        let Ok(rep) = x.get(rep) else { continue };
        let Some(items) = rep.records().find(|(_, v)| v.len() >= 2 && x.list(&v[1]).is_some()).and_then(|(_, v)| x.list(&v[1])) else { continue };
        let v = out.entry(of).or_default();
        for item in items.iter().filter_map(Value::as_reference) {
            let Ok(i) = x.get(item) else { continue };
            if let Some(c) = i.record("CARTESIAN_POINT") {
                let coords: Vec<f64> = c.get(1).and_then(|l| x.list(l)).unwrap_or(&[]).iter().filter_map(|v| x.real(v)).collect();
                if coords.len() == 3 {
                    v.centroid = Some(Point::new(coords[0] * length, coords[1] * length, coords[2] * length));
                }
            } else if let Some(m) = i.record("MEASURE_REPRESENTATION_ITEM") {
                let Some((kind, values)) = m.get(1).and_then(|t| x.typed(t)) else { continue };
                let Some(value) = values.first().and_then(|v| x.real(v)) else { continue };
                match kind {
                    "VOLUME_MEASURE" => v.volume = Some(value * length.powi(3)),
                    "AREA_MEASURE" => v.area = Some(value * length.powi(2)),
                    _ => {}
                }
            }
        }
    }
    out
}

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

/// A straight pcurve of an edge on a periodic curve (a meridian of a sphere), taken a
/// period of the curve along where it leaves the surface's range in a direction that is not
/// periodic: the edge's range is put in the curve's first period, the file's may be another,
/// and the surface's formula, periodic in that direction too, hides the difference.
fn within_bounds(pcurve: Curve2, made: &EdgeMade, surface: &Surface) -> Curve2 {
    use crate::geometry::ParametricSurface;
    let (Curve2::Line(line), Some(period)) = (&pcurve, made.curve.period()) else { return pcurve };
    let (u0, u1, v0, v1) = surface.bounds();
    let slack = 1e-9 * (1.0 + period);
    let outside = |l: &Line2| {
        [made.first, made.last].iter().any(|&t| {
            let p = l.origin + Vector2::from(l.direction) * t;
            (surface.u_period().is_none() && (p.x < u0 - slack || p.x > u1 + slack)) || (surface.v_period().is_none() && (p.y < v0 - slack || p.y > v1 + slack))
        })
    };
    if !outside(line) {
        return pcurve;
    }
    for k in [1.0, -1.0, 2.0, -2.0] {
        let moved = Line2 { origin: line.origin + Vector2::from(line.direction) * (k * period), direction: line.direction };
        if !outside(&moved) {
            return Curve2::Line(moved);
        }
    }
    pcurve
}

/// The parameter range along `curve` from `start` to `end` (a whole period when the edge
/// closes on itself).
fn range(curve: &Curve, start: Point, end: Point, closed: bool) -> (f64, f64) {
    let period = curve.period();
    let param = |p: Point| match curve {
        Curve::Line(l) => (p - l.origin).dot(Vector::from(l.direction)),
        // (Unbounded conics in closed form: from the point's place along the frame's y.)
        Curve::Hyperbola(h) if h.minor_radius() > 0.0 => ((p - h.frame().origin()).dot(Vector::from(h.frame().y())) / h.minor_radius()).asinh(),
        Curve::Parabola(c) => (p - c.frame().origin()).dot(Vector::from(c.frame().y())),
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
    // (An edge between two vertices found at one place on a curve not periodic, the curve
    // coming back about there: the edge runs it round, to the end where its end vertex lies
    // near, as a curve closed only to the file's accuracy.)
    if period.is_none() && !closed {
        let (lo, hi) = (curve.first_parameter(), curve.last_parameter());
        let span = (hi - lo).abs();
        let apart = start.distance(end);
        if lo.is_finite() && hi.is_finite() && (last - first).abs() <= 1e-6 * span && apart > 0.0 {
            let near = |q: Point| q.distance(end) <= 4.0 * apart + 1e-7;
            if (first - lo).abs() <= 1e-6 * span && near(curve.point(hi)) {
                return (lo, hi);
            }
            if (first - hi).abs() <= 1e-6 * span && near(curve.point(lo)) {
                return (lo, hi);
            }
        }
    }
    // A closed curve (not periodic) run through its closing point: a vertex there is at the
    // end of the curve's range where the edge ends, at its start where the edge starts.
    if period.is_none() && last < first {
        let (lo, hi) = (curve.first_parameter(), curve.last_parameter());
        if lo.is_finite() && hi.is_finite() {
            let (a, b) = (curve.point(lo), curve.point(hi));
            let size = 1.0 + a.x.abs().max(a.y.abs()).max(a.z.abs());
            let span = (hi - lo).abs().max(1e-300);
            if a.distance(b) <= 1e-9 * size {
                if (last - lo).abs() <= 1e-9 * span {
                    last = hi;
                } else if (first - hi).abs() <= 1e-9 * span {
                    first = lo;
                }
            }
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

    /// A cone whose side, a surface of revolution of a line, closes at the apex without an
    /// edge there (the seam's two uses meet across it): its seam on its two sides, the
    /// apex's edge added, the cone whole.
    #[test]
    fn a_side_closing_at_a_pole_without_its_edge() {
        let mut model = Model::new();
        let shape = read(&mut model, &mut include_str!("tests/cone_without_apex_edge.step").as_bytes()).unwrap();
        let report = crate::validity::check(&model, shape);
        assert!(report.is_valid(), "{report:?}");
        let volume = crate::properties::volume(&model, shape).mass;
        let expected = std::f64::consts::PI * 100.0 * 5.0 / 3.0;
        assert!((volume - expected).abs() < 1e-6 * expected, "volume {volume}");
    }

    /// A quarter ball whose sphere face is bounded by a meridian running over the pole, as
    /// files write a corner of rounds meeting at a sphere's pole (no seam, no edge at the
    /// pole): the meridian split at the pole, the halves joined round it, the ball whole.
    #[test]
    fn a_sphere_bounded_by_a_meridian_over_its_pole() {
        let mut model = Model::new();
        let shape = read(&mut model, &mut include_str!("tests/quarter_ball_over_its_pole.step").as_bytes()).unwrap();
        let report = crate::validity::check(&model, shape);
        assert!(report.is_valid(), "{report:?}");
        let volume = crate::properties::volume(&model, shape).mass;
        let expected = std::f64::consts::PI * 1000.0 / 3.0;
        assert!((volume - expected).abs() < 1e-6 * expected, "volume {volume}");
        let area = crate::properties::area(&model, shape).mass;
        let expected = std::f64::consts::PI * 100.0 * 2.0;
        assert!((area - expected).abs() < 1e-6 * expected, "area {area}");
    }

    /// An AP242 file of one part whose shape is the tessellated items `items` (entities
    /// numbered from #20, the last the item of the representation).
    fn tessellated_file(items: &[String]) -> String {
        let mut data = String::from(
            "#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n\
             #2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2011,#1);\n\
             #3=PRODUCT_CONTEXT('',#1,'mechanical');\n#4=PRODUCT('p','p','',(#3));\n\
             #5=PRODUCT_DEFINITION_FORMATION('','',#4);\n#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n\
             #7=PRODUCT_DEFINITION('design','',#5,#6);\n#8=PRODUCT_DEFINITION_SHAPE('','',#7);\n\
             #9=SHAPE_DEFINITION_REPRESENTATION(#8,#19);\n#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
             #11=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n#12=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n\
             #13=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#10,'distance_accuracy_value','');\n\
             #14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''));\n",
        );
        let last = 19 + items.len();
        data += &format!("#19=TESSELLATED_SHAPE_REPRESENTATION('',(#{last}),#14);\n");
        for (i, item) in items.iter().enumerate() {
            data += &format!("#{}={item};\n", 20 + i);
        }
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('t','',(''),(''),'','','');\n\
             FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF {{ 1 0 10303 442 1 1 4 }}'));\nENDSEC;\nDATA;\n{data}ENDSEC;\nEND-ISO-10303-21;\n"
        )
    }

    /// The corners of a 2 × 3 × 4 box, point k (1-based) at bits k − 1 = (x, y, z), as a
    /// coordinates list.
    fn box_corners() -> String {
        let points: Vec<String> = (0..8).map(|k| format!("({:.1},{:.1},{:.1})", 2.0 * (k & 1) as f64, 3.0 * ((k >> 1) & 1) as f64, 4.0 * ((k >> 2) & 1) as f64)).collect();
        format!("COORDINATES_LIST('',8,({}))", points.join(","))
    }

    /// The box's sides, each counter-clockwise seen from outside.
    const BOX_SIDES: [[usize; 4]; 6] = [[1, 3, 4, 2], [5, 6, 8, 7], [1, 2, 6, 5], [3, 7, 8, 4], [1, 5, 7, 3], [2, 4, 8, 6]];

    fn triangles(sides: &[[usize; 4]], node: impl Fn(usize) -> usize) -> String {
        let t: Vec<String> = sides.iter().flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]]).map(|t| format!("({},{},{})", node(t[0]), node(t[1]), node(t[2]))).collect();
        t.join(",")
    }

    fn assert_box(model: &Model, shape: Shape, faces: usize) {
        assert_eq!(shape.kind(), NodeKind::Solid);
        assert_eq!(model.explore(shape, NodeKind::Face).count(), faces);
        let report = crate::validity::check(model, shape);
        assert!(report.is_valid(), "{report:?}");
        let volume = crate::properties::volume(model, shape).mass;
        assert!((volume - 24.0).abs() < 1e-9, "volume {volume}");
        let area = crate::properties::area(model, shape).mass;
        assert!((area - 52.0).abs() < 1e-9, "area {area}");
    }

    /// A tessellated solid of one triangulated face: a solid of the triangles' planar faces,
    /// sharing their edges, the box's volume and area.
    #[test]
    fn a_tessellated_solid_of_one_triangulated_face() {
        let file = tessellated_file(&[
            box_corners(),
            format!("TRIANGULATED_FACE('',#20,8,(),$,(),({}))", triangles(&BOX_SIDES, |k| k)),
            "TESSELLATED_SOLID('',(#21),$)".into(),
        ]);
        let mut model = Model::new();
        let shape = read(&mut model, &mut file.as_bytes()).unwrap();
        assert_box(&model, shape, 12);
    }

    /// Two triangulated faces over one coordinates list, one naming its points through
    /// `pnindex` in reverse: their nodes at the same points one vertex, the shell closed.
    #[test]
    fn triangulated_faces_sharing_their_points() {
        let file = tessellated_file(&[
            box_corners(),
            format!("TRIANGULATED_FACE('',#20,8,(),$,(),({}))", triangles(&BOX_SIDES[..3], |k| k)),
            format!("TRIANGULATED_FACE('',#20,8,(),$,(8,7,6,5,4,3,2,1),({}))", triangles(&BOX_SIDES[3..], |k| 9 - k)),
            "TESSELLATED_SOLID('',(#21,#22),$)".into(),
        ]);
        let mut model = Model::new();
        let shape = read(&mut model, &mut file.as_bytes()).unwrap();
        assert_box(&model, shape, 12);
        assert_eq!(model.unique_subshapes(shape, NodeKind::Vertex).len(), 8);
    }

    /// A complex triangulated face of fans, one per side, in a tessellated shell: each fan
    /// taken about its first node, a fan (a, b, c, …) giving (a, c, b), (a, d, c), … as
    /// OCCT reads it.
    #[test]
    fn a_complex_triangulated_face_of_fans() {
        let fans: Vec<String> = BOX_SIDES.iter().map(|&[a, b, c, d]| format!("({a},{d},{c},{b})")).collect();
        let file = tessellated_file(&[
            box_corners(),
            format!("COMPLEX_TRIANGULATED_FACE('',#20,8,(),$,(),(),({}))", fans.join(",")),
            "TESSELLATED_SHELL('',(#21),$)".into(),
        ]);
        let mut model = Model::new();
        let shell = read(&mut model, &mut file.as_bytes()).unwrap();
        assert_eq!(shell.kind(), NodeKind::Shell);
        let solid = model.new_empty(NodeKind::Solid);
        model.add(solid, shell).unwrap();
        assert_box(&model, solid, 12);
    }

    /// A triangle naming a node the face has not is refused.
    #[test]
    fn a_triangle_beyond_the_nodes_refused() {
        let file = tessellated_file(&[box_corners(), "TRIANGULATED_FACE('',#20,8,(),$,(),((1,2,9)))".into(), "TESSELLATED_SOLID('',(#21),$)".into()]);
        assert!(read(&mut Model::new(), &mut file.as_bytes()).is_err());
    }

    /// Unbounded conics' edges get their parameters from their vertices in closed form.
    #[test]
    fn hyperbola_and_parabola_edges_ranged_by_their_vertices() {
        use crate::geometry::ParametricCurve;
        let frame = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 0.0, 0.0));
        for curve in [Curve::Hyperbola(Hyperbola::new(frame, 2.0, 0.5)), Curve::Parabola(Parabola::new(frame, 0.75))] {
            let (a, b) = (-0.7, 1.3);
            let (first, last) = range(&curve, curve.point(a), curve.point(b), false);
            assert!((first - a).abs() < 1e-12 && (last - b).abs() < 1e-12, "{curve:?}: {first} {last}");
        }
    }
}
