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

pub use write::{write, write_assembly};

use crate::hashing::FastMap as HashMap;
use std::io::Read;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Circle, Curve, Ellipse, Line};
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Line2};
use crate::geometry::extrema::nearest_on_curve;
use crate::geometry::surface::{
    BSplineSurface, Cone, Cylinder, ExtrusionSurface, KnotVector, OffsetSurface, Plane, RevolutionSurface, Sphere, Surface, Torus,
};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Vector};
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

pub mod binstep;
pub use parse::{parse, Exchange, Instance, Value};

/// A node of a file's assembly structure: a product (part or assembly) as placed in its
/// parent.
#[derive(Debug, Clone)]
pub struct Node {
    /// The occurrence's name (or the product's, for a root or an unnamed occurrence).
    pub name: String,
    /// The product's own shape (a solid, or a compound of its solids and shells); `None`
    /// for an assembly of other products only.
    pub shape: Option<Shape>,
    /// The colour of all its solids, if they share one.
    pub color: Option<[f64; 3]>,
    /// Where the product sits in its parent.
    pub placement: Transform,
    pub children: Vec<Node>,
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
    reader.apply_tolerances()?;
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
    colors: HashMap<usize, [f64; 3]>,
    /// Tolerances needed to cover the gaps between the file's curves, pcurves and points.
    edge_tolerance: HashMap<EdgeId, f64>,
    vertex_tolerance: HashMap<VertexId, f64>,
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
            edge_tolerance: HashMap::default(),
            vertex_tolerance: HashMap::default(),
        };
        reader.colors = colors(x);
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
        self.finish_face(face, &surface, loops, same_sense)
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
        self.finish_face(face, &surface, loops, true)
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
            Some(pc) if fits(&pc, made, surface) => pc,
            // Parametrized otherwise but along the edge on a free-form surface (where
            // projecting is costly): kept, the tolerances covering the difference.
            Some(pc) if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) && follows(&pc, made, surface) => pc,
            given => match (project(made, surface, given.is_some()), given) {
                (Ok(pc), _) => pc,
                (Err(_), Some(pc)) => pc,
                (Err(e), None) => return Err(e),
            },
        })
    }

    /// The face from its loops (outer?, coedges with pcurves): seams and poles placed,
    /// outer boundaries counter-clockwise, holes in the outer one's period.
    fn finish_face(&mut self, face: Shape, surface: &Surface, loops: Vec<(bool, Vec<(EdgeId, bool, Curve2, f64, f64)>)>, same_sense: bool) -> Result<Shape> {
        let surface = surface.clone();
        // (Seams placed by what all the face's loops span: a hole may tell which side.)
        let spans = face_spans(&loops.iter().map(|(_, c)| c.as_slice()).collect::<Vec<_>>(), &surface);
        let mut placed = Vec::with_capacity(loops.len());
        for (outer, mut coedges) in loops {
            seams_to_sides(&mut coedges, &surface, spans);
            placed.push((outer, self.close_at_poles(coedges, &surface)?));
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
                self.note_gaps(edge, &pcurve, first, last, &surface);
                let c = self.model.add_coedge(edge, if forward { Orientation::Forward } else { Orientation::Reversed })?;
                self.model.set_pcurve_with_range(c, pcurve, first, last)?;
                ids.push(c);
            }
            let wire = self.model.add_wire(&ids)?;
            self.model.add(face, wire.into())?;
        }
        Ok(if same_sense { face } else { face.reversed() })
    }

    /// The loop with degenerate edges added where consecutive pcurves end at the same point
    /// of a singular part of the surface (a cone's apex, a sphere's pole) but apart in its
    /// parameters: files leave those edges out.
    fn close_at_poles(&mut self, coedges: Vec<(EdgeId, bool, Curve2, f64, f64)>, surface: &Surface) -> Result<Vec<(EdgeId, bool, Curve2, f64, f64)>> {
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
            let (a, b) = self.model.edge_vertices(Shape::from(coedges[i].0));
            let vertex = if coedges[i].1 { b } else { a };
            let v = VertexId::try_from(vertex.node).expect("a vertex");
            let edge = self.model.add_edge(v, v, crate::tolerance::LINEAR)?;
            let pcurve = Curve2::BSpline(BSplineCurve2::new(vec![end, start], vec![0.0, 1.0], vec![2, 2], 1));
            out.push((edge, true, pcurve, 0.0, 1.0));
        }
        Ok(out)
    }

    /// Widens the tolerances of an edge and its vertices to cover how far its pcurve on a
    /// surface and its vertices are from its curve.
    fn note_gaps(&mut self, edge: EdgeId, pcurve: &Curve2, first: f64, last: f64, surface: &Surface) {
        let Some((curve, _, _)) = self.model.edge_curve(Shape::from(edge)) else { return };
        let mut deviation: f64 = 0.0;
        for i in 0..=8 {
            let t = first + (last - first) * i as f64 / 8.0;
            let uv = pcurve.point(t);
            deviation = deviation.max(surface.point(uv.x, uv.y).distance(curve.point(t)));
        }
        let (a, b) = self.model.edge_vertices(Shape::from(edge));
        let ends = [(a, first), (b, last)];
        for (v, t) in ends {
            let id = VertexId::try_from(v.node).expect("a vertex");
            let p = self.model.point(v);
            let uv = pcurve.point(t);
            let gap = p.distance(curve.point(t)).max(p.distance(surface.point(uv.x, uv.y)));
            let tol = self.vertex_tolerance.entry(id).or_insert(crate::tolerance::LINEAR);
            *tol = tol.max(gap * 1.1).max(deviation * 1.1);
            deviation = deviation.max(gap);
        }
        let tol = self.edge_tolerance.entry(edge).or_insert(crate::tolerance::LINEAR);
        *tol = tol.max(deviation * 1.1);
    }

    /// Sets the tolerances noted while reading.
    fn apply_tolerances(&mut self) -> Result<()> {
        for (&edge, &tol) in &self.edge_tolerance {
            if tol > crate::tolerance::LINEAR {
                self.model.set_edge_tolerance(edge, tol)?;
            }
        }
        for (&v, &tol) in &self.vertex_tolerance {
            if tol > crate::tolerance::LINEAR {
                let p = self.model.point(Shape::from(v));
                self.model.set_vertex_point(v.into(), p, tol)?;
            }
        }
        Ok(())
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
            nodes.push(Node { name: String::new(), shape, color: None, placement: Transform::IDENTITY, children: Vec::new() });
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
        let color = match colors.first() {
            Some(&Some(c)) if colors.iter().all(|&k| k == Some(c)) => Some(c),
            _ => None,
        };
        let shape = self.compound_of(bodies)?;
        let mut kids = Vec::new();
        for (child, occurrence, occurrence_name) in children.get(&definition).cloned().unwrap_or_default() {
            let child_name = if occurrence_name.is_empty() { self.product_name(child)? } else { occurrence_name };
            let at = placement_of.get(&occurrence).copied().unwrap_or(Transform::IDENTITY);
            kids.push(self.node(child, child_name, at, reps_of, linked, children, placement_of, depth + 1)?);
        }
        Ok(Node { name, shape, color, placement, children: kids })
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
fn colors(x: &Exchange) -> HashMap<usize, [f64; 3]> {
    let mut out = HashMap::default();
    for (_, inst) in x.of_type("STYLED_ITEM") {
        let p = inst.record("STYLED_ITEM").expect("record");
        let Some(item) = p.get(2).and_then(Value::as_reference) else { continue };
        let mut stack: Vec<&Value> = vec![&p[1]];
        let mut depth = 0;
        let mut found = None;
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
                    if let Some(c) = i.record("COLOUR_RGB") {
                        let rgb = [x.real(&c[1]), x.real(&c[2]), x.real(&c[3])];
                        if let [Some(r), Some(g), Some(b)] = rgb {
                            found = Some([r, g, b]);
                            break;
                        }
                    }
                    if let Some(c) = i.record("DRAUGHTING_PRE_DEFINED_COLOUR") {
                        found = c[0].as_str().and_then(named_colour);
                        if found.is_some() {
                            break;
                        }
                    }
                    for (_, params) in i.records() {
                        stack.extend(params.iter());
                    }
                }
                _ => {}
            }
        }
        if let Some(c) = found {
            out.entry(item).or_insert(c);
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

/// True if the pcurve follows the edge's curve on the surface.
fn fits(pcurve: &Curve2, edge: &EdgeMade, surface: &Surface) -> bool {
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
fn follows(pcurve: &Curve2, edge: &EdgeMade, surface: &Surface) -> bool {
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
            crate::geometry::extrema::nearest_on_curve(p, &edge.curve, edge.first, edge.last).distance <= allowed
        })
}

/// The edge's curve projected onto the surface (`has_pcurve`: the file gives a pcurve to fall
/// back on).
fn project(edge: &EdgeMade, surface: &Surface, has_pcurve: bool) -> Result<Curve2> {
    use crate::geometry::project as pr;
    if let Some(pc) = pr::analytic_pcurve(&edge.curve, surface) {
        return Ok(pc);
    }
    // An edge well off the surface (bad data) fails at once, before costly fitting.
    if surface.parameters(edge.curve.point(edge.first)).is_none() {
        let mut start = None;
        for i in 0..=4 {
            let p = edge.curve.point(edge.first + (edge.last - edge.first) * i as f64 / 4.0);
            let (u, v, d) = crate::geometry::extrema::nearest_on_surface(p, surface, start);
            start = Some((u, v));
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
    pr::pcurve_on_range_within(&edge.curve, surface, edge.first, edge.last, tiers)
        .or_else(|| if freeform || (has_pcurve && off_surface()) { None } else { pr::polyline_pcurve_on_range(&edge.curve, surface, edge.first, edge.last, 1e-6) })
        .map(Curve2::BSpline)
        .ok_or_else(|| invalid("an edge that does not project onto its face".into()))
}

/// Pcurves moved by whole periods of a periodic surface so that each starts where the one
/// before it ends (the file's and projected pcurves may lie in different periods; a seam's
/// two uses end up on its two sides).
fn align_periods(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface) {
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
fn face_spans(loops: &[&[(EdgeId, bool, Curve2, f64, f64)]], surface: &Surface) -> [(f64, f64); 2] {
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
fn seams_to_sides(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], surface: &Surface, spans: [(f64, f64); 2]) {
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
fn uv_box(coedges: &[(EdgeId, bool, Curve2, f64, f64)]) -> (Point2, Point2) {
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
fn into_period(coedges: &mut [(EdgeId, bool, Curve2, f64, f64)], target: (Point2, Point2), surface: &Surface) {
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
fn uv_area(coedges: &[(EdgeId, bool, Curve2, f64, f64)]) -> f64 {
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
