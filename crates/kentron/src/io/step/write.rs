//! Writing STEP files (AP214): solids and shells as boundary representations, products
//! placed in assemblies, and their colours.
//!
//! Each product's shape is written in its own frame, its placement in its parent as an
//! item-defined transformation. Edges carry their 3D curves only (readers project them onto
//! the faces' surfaces). Frames are written right-handed, as STEP has them: where a
//! reflected (left-handed) surface turns its normal round, the face says so by its sense.

use crate::hashing::FastMap as HashMap;
use std::fmt::Write as _;
use std::io::Write;

use super::parse::real;
use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::surface::Surface;
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Frame, Point, Transform, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::{LengthUnit, Node};

/// Writes `shape` as one product named `name`, with an optional colour (in millimetres).
#[doc(alias = "STEPControl_Writer")]
pub fn write(model: &Model, shape: Shape, name: &str, color: Option<[f64; 3]>, out: &mut dyn Write) -> Result<()> {
    write_in(model, shape, name, color, LengthUnit::Millimetre, out)
}

/// [`write`], the shape's lengths taken to be in `unit`: written as they are, the file
/// declaring that unit.
pub fn write_in(model: &Model, shape: Shape, name: &str, color: Option<[f64; 3]>, unit: LengthUnit, out: &mut dyn Write) -> Result<()> {
    let root = Node { name: name.to_string(), shape: Some(shape), color, ..Default::default() };
    write_assembly_in(model, &[root], unit, out)
}

/// Writes products: each node a product (an assembly of its children, with its own shape if
/// any), placed in its parent (in millimetres).
#[doc(alias = "STEPCAFControl_Writer")]
pub fn write_assembly(model: &Model, roots: &[Node], out: &mut dyn Write) -> Result<()> {
    write_assembly_in(model, roots, LengthUnit::Millimetre, out)
}

/// [`write_assembly`], lengths taken to be in `unit` (see [`write_in`]).
pub fn write_assembly_in(model: &Model, roots: &[Node], unit: LengthUnit, out: &mut dyn Write) -> Result<()> {
    let mut w = Writer::new(model);
    w.unit = unit;
    let (application, product_context, definition_context) = w.contexts();
    for root in roots {
        w.product(root, &application, (product_context, definition_context))?;
    }
    // The layers, each with the faces and edges on it.
    for (name, items) in std::mem::take(&mut w.layers) {
        let list = items.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",");
        let name = w.escape(&name);
        w.add(format!("PRESENTATION_LAYER_ASSIGNMENT('{name}','visible',({list}))"));
    }
    let mut text = String::new();
    text.push_str("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('kentron model'),'2;1');\n");
    let _ = writeln!(text, "FILE_NAME('{}','',(''),(''),'kentron','kentron','');", w.escape(roots.first().map_or("", |r| r.name.as_str())));
    text.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n");
    for (i, line) in w.lines.iter().enumerate() {
        let _ = writeln!(text, "#{} = {};", i + 1, line);
    }
    text.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    out.write_all(text.as_bytes()).map_err(|e| Error::InvalidInput(format!("writing STEP: {e}")))
}


/// A face's or edge's own name, colour (with its opacity) and layer.
#[derive(Debug, Clone, Default)]
struct Mark {
    name: Option<String>,
    color: Option<[f64; 4]>,
    layer: Option<String>,
}

struct Writer<'a> {
    model: &'a Model,
    lines: Vec<String>,
    /// The marks of the product being written's faces and edges, by node.
    marks: HashMap<NodeId, Mark>,
    /// Styles of the product being written's faces.
    face_styles: Vec<usize>,
    /// Layers and the items on each, in the order met.
    layers: Vec<(String, Vec<usize>)>,
    /// Shapes written, by node and placement within the product being written.
    vertices: HashMap<(NodeId, [u64; 12]), usize>,
    edges: HashMap<(NodeId, [u64; 12]), usize>,
    context: Option<usize>,
    unit: LengthUnit,
}

fn key(t: &Transform) -> [u64; 12] {
    let r = t.rows();
    let mut k = [0u64; 12];
    for i in 0..3 {
        for j in 0..4 {
            k[i * 4 + j] = r[i][j].to_bits();
        }
    }
    k
}

impl<'a> Writer<'a> {
    fn new(model: &'a Model) -> Self {
        Writer {
            model,
            lines: Vec::new(),
            marks: HashMap::default(),
            face_styles: Vec::new(),
            layers: Vec::new(),
            vertices: HashMap::default(),
            edges: HashMap::default(),
            context: None,
            unit: LengthUnit::Millimetre,
        }
    }

    /// Adds an entity; its number.
    fn add(&mut self, line: String) -> usize {
        self.lines.push(line);
        self.lines.len()
    }

    /// A string in STEP's form: quotes doubled, non-ASCII characters encoded (`\X2\..\X0\`).
    fn escape(&self, s: &str) -> String {
        let mut out = String::new();
        for c in s.chars() {
            match c {
                '\'' => out.push_str("''"),
                '\\' => out.push_str("\\\\"),
                c if c.is_ascii() && !c.is_ascii_control() => out.push(c),
                c => {
                    let mut units = [0u16; 2];
                    out.push_str("\\X2\\");
                    for u in c.encode_utf16(&mut units) {
                        let _ = write!(out, "{:04X}", u);
                    }
                    out.push_str("\\X0\\");
                }
            }
        }
        out
    }

    fn contexts(&mut self) -> (usize, usize, usize) {
        let application = self.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')".into());
        self.add(format!("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{application})"));
        let product = self.add(format!("PRODUCT_CONTEXT('',#{application},'mechanical')"));
        let definition = self.add(format!("PRODUCT_DEFINITION_CONTEXT('part definition',#{application},'design')"));
        (application, product, definition)
    }

    /// The geometric context with millimetres, radians and the uncertainty (shared).
    fn geometric_context(&mut self) -> usize {
        if let Some(c) = self.context {
            return c;
        }
        let length = match self.unit.si_prefix() {
            Some(prefix) => self.add(format!("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT({prefix},.METRE.) )")),
            // A unit defined as so many millimetres (as OCCT writes inches).
            None => {
                let mm = self.add("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )".into());
                let measure = self.add(format!("LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE({}),#{mm})", real(self.unit.millimetres())));
                let exponents = self.add("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)".into());
                self.add(format!("( CONVERSION_BASED_UNIT('{}',#{measure}) LENGTH_UNIT() NAMED_UNIT(#{exponents}) )", self.unit.name()))
            }
        };
        let angle = self.add("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".into());
        let solid = self.add("( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )".into());
        let uncertainty = self.add(format!("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{length},'distance_accuracy_value','confusion accuracy')"));
        let c = self.add(format!(
            "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{uncertainty})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{length},#{angle},#{solid})) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') )"
        ));
        self.context = Some(c);
        c
    }

    /// A product for `node` and its children; its product definition and representation.
    fn product(&mut self, node: &Node, application: &usize, (product_context, definition_context): (usize, usize)) -> Result<(usize, usize)> {
        let _ = application;
        let name = self.escape(&node.name);
        let product = self.add(format!("PRODUCT('{name}','{name}','',(#{product_context}))"));
        self.add(format!("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{product}))"));
        let formation = self.add(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let definition = self.add(format!("PRODUCT_DEFINITION('design','',#{formation},#{definition_context})"));
        let definition_shape = self.add(format!("PRODUCT_DEFINITION_SHAPE('','',#{definition})"));
        if let Some(material) = &node.material {
            self.material(material, definition_shape);
        }
        let context = self.geometric_context();
        let origin = self.axis(&Frame::WORLD);
        // The faces' and edges' own names, colours and layers (a wire's are its edges', a
        // solid's or shell's its faces').
        self.marks.clear();
        for sub in &node.subshapes {
            let targets: Vec<NodeId> = match sub.shape.kind() {
                NodeKind::Face | NodeKind::Edge => vec![sub.shape.node],
                NodeKind::Wire => self.model.unique_subshapes(sub.shape, NodeKind::Edge).iter().map(|e| e.node).collect(),
                _ => self.model.unique_subshapes(sub.shape, NodeKind::Face).iter().map(|f| f.node).collect(),
            };
            for t in targets {
                let mark = self.marks.entry(t).or_default();
                mark.name = sub.name.clone().or(mark.name.take());
                mark.color = sub.color.or(mark.color);
                mark.layer = sub.layer.clone().or(mark.layer.take());
            }
        }
        // The product's own geometry.
        let (mut items, mut styled) = (vec![origin], Vec::new());
        let mut has_shells = false;
        if let Some(shape) = node.shape {
            let local = Shape { orientation: shape.orientation, ..shape };
            for (item, shell) in self.bodies(local)? {
                items.push(item);
                has_shells |= shell;
                if let Some(c) = node.color {
                    styled.push(self.style(item, c, node.alpha));
                }
            }
        }
        styled.append(&mut self.face_styles);
        // Placements of the children.
        let mut child_axes = Vec::new();
        for child in &node.children {
            let frame = Frame::WORLD.transformed(&child.placement);
            let axis = self.axis(&frame);
            items.push(axis);
            child_axes.push(axis);
        }
        let list = items.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",");
        let kind = if node.shape.is_some() && !has_shells && node.children.is_empty() { "ADVANCED_BREP_SHAPE_REPRESENTATION" } else { "SHAPE_REPRESENTATION" };
        let representation = self.add(format!("{kind}('',({list}),#{context})"));
        self.add(format!("SHAPE_DEFINITION_REPRESENTATION(#{definition_shape},#{representation})"));
        if !styled.is_empty() {
            let list = styled.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",");
            self.add(format!("MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('',({list}),#{context})"));
        }
        for (k, (child, axis)) in node.children.iter().zip(child_axes).enumerate() {
            let (child_definition, child_representation) = self.product(child, application, (product_context, definition_context))?;
            let transformation = self.add(format!("ITEM_DEFINED_TRANSFORMATION('','',#{origin},#{axis})"));
            let relationship = self.add(format!(
                "( REPRESENTATION_RELATIONSHIP('','',#{child_representation},#{representation}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{transformation}) SHAPE_REPRESENTATION_RELATIONSHIP() )"
            ));
            let child_name = self.escape(&child.name);
            let usage = self.add(format!("NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{child_name}','',#{definition},#{child_definition},$)", k + 1));
            let placement = self.add(format!("PRODUCT_DEFINITION_SHAPE('Placement','Placement of an item',#{usage})"));
            self.add(format!("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{relationship},#{placement})"));
        }
        Ok((definition, representation))
    }

    /// A colour on an item, with its opacity where it is not opaque (written as the
    /// rendering's transparency).
    fn style(&mut self, item: usize, rgb: [f64; 3], alpha: Option<f64>) -> usize {
        let colour = self.colour(rgb);
        let fill_colour = self.add(format!("FILL_AREA_STYLE_COLOUR('',#{colour})"));
        let fill = self.add(format!("FILL_AREA_STYLE('',(#{fill_colour}))"));
        let area = self.add(format!("SURFACE_STYLE_FILL_AREA(#{fill})"));
        let side = match alpha.filter(|&a| a < 1.0) {
            Some(a) => {
                let transparent = self.add(format!("SURFACE_STYLE_TRANSPARENT({})", real((1.0 - a).clamp(0.0, 1.0))));
                let rendering = self.add(format!("SURFACE_STYLE_RENDERING_WITH_PROPERTIES(.NORMAL_SHADING.,#{colour},(#{transparent}))"));
                self.add(format!("SURFACE_SIDE_STYLE('',(#{area},#{rendering}))"))
            }
            None => self.add(format!("SURFACE_SIDE_STYLE('',(#{area}))")),
        };
        let usage = self.add(format!("SURFACE_STYLE_USAGE(.BOTH.,#{side})"));
        let assignment = self.add(format!("PRESENTATION_STYLE_ASSIGNMENT((#{usage}))"));
        self.add(format!("STYLED_ITEM('color',(#{assignment}),#{item})"))
    }

    /// A colour: one of the predefined ones by its name where it is exactly that (as OCCT
    /// writes them), else by its components.
    fn colour(&mut self, [r, g, b]: [f64; 3]) -> usize {
        const NAMED: [(&str, [f64; 3]); 8] = [
            ("black", [0.0, 0.0, 0.0]),
            ("white", [1.0, 1.0, 1.0]),
            ("red", [1.0, 0.0, 0.0]),
            ("green", [0.0, 1.0, 0.0]),
            ("blue", [0.0, 0.0, 1.0]),
            ("yellow", [1.0, 1.0, 0.0]),
            ("magenta", [1.0, 0.0, 1.0]),
            ("cyan", [0.0, 1.0, 1.0]),
        ];
        match NAMED.iter().find(|(_, c)| (c[0] - r).abs() <= 1e-9 && (c[1] - g).abs() <= 1e-9 && (c[2] - b).abs() <= 1e-9) {
            Some((name, _)) => self.add(format!("DRAUGHTING_PRE_DEFINED_COLOUR('{name}')")),
            None => self.add(format!("COLOUR_RGB('',{},{},{})", real(r), real(g), real(b))),
        }
    }

    /// A product's material: its designation, and a property of its description, density
    /// and the density's unit.
    fn material(&mut self, material: &super::Material, definition_shape: usize) {
        let name = self.escape(&material.name);
        self.add(format!("MATERIAL_DESIGNATION('{name}',#{definition_shape})"));
        let property = self.add(format!("PROPERTY_DEFINITION('material property','{name}',#{definition_shape})"));
        let description = self.escape(&material.description);
        let unit = self.escape(&material.density_unit);
        let items = [
            self.add(format!("DESCRIPTIVE_REPRESENTATION_ITEM('description','{description}')")),
            self.add(format!("DESCRIPTIVE_REPRESENTATION_ITEM('density','{}')", real(material.density))),
            self.add(format!("DESCRIPTIVE_REPRESENTATION_ITEM('density unit','{unit}')")),
        ];
        let list = items.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",");
        let context = self.geometric_context();
        let representation = self.add(format!("REPRESENTATION('material property',({list}),#{context})"));
        self.add(format!("PROPERTY_DEFINITION_REPRESENTATION(#{property},#{representation})"));
    }

    /// Puts `item` on the layer `name`.
    fn layer(&mut self, name: &str, item: usize) {
        match self.layers.iter_mut().find(|(n, _)| n == name) {
            Some((_, items)) => items.push(item),
            None => self.layers.push((name.to_string(), vec![item])),
        }
    }

    /// The solids (as manifold solid B-reps) and the shells outside solids (as shell-based
    /// surface models) of a shape; each item with whether it is a shell model.
    fn bodies(&mut self, shape: Shape) -> Result<Vec<(usize, bool)>> {
        let model = self.model;
        let mut out = Vec::new();
        let mut in_solids: Vec<(NodeId, [u64; 12])> = Vec::new();
        for solid in model.explore(shape, NodeKind::Solid).collect::<Vec<_>>() {
            let mut shells = Vec::new();
            for shell in model.children(solid).collect::<Vec<_>>() {
                in_solids.push((shell.node, key(&shell.placement)));
                shells.push(self.shell(shell, true)?);
            }
            let Some((&outer, voids)) = shells.split_first() else { continue };
            let item = if voids.is_empty() {
                self.add(format!("MANIFOLD_SOLID_BREP('',#{outer})"))
            } else {
                let list = voids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
                self.add(format!("BREP_WITH_VOIDS('',#{outer},({list}))"))
            };
            out.push((item, false));
        }
        let free: Vec<Shape> = model.explore(shape, NodeKind::Shell).filter(|s| !in_solids.contains(&(s.node, key(&s.placement)))).collect();
        if !free.is_empty() {
            let mut shells = Vec::new();
            for shell in free {
                shells.push(self.shell(shell, false)?);
            }
            let list = shells.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
            out.push((self.add(format!("SHELL_BASED_SURFACE_MODEL('',({list}))")), true));
        }
        if out.is_empty() && model.explore(shape, NodeKind::Face).next().is_some() {
            // Loose faces: an open shell of them.
            let faces: Vec<Shape> = model.explore(shape, NodeKind::Face).collect();
            let mut ids = Vec::new();
            for f in faces {
                ids.push(self.face(f)?);
            }
            let list = ids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
            let shell = self.add(format!("OPEN_SHELL('',({list}))"));
            out.push((self.add(format!("SHELL_BASED_SURFACE_MODEL('',(#{shell}))")), true));
        }
        Ok(out)
    }

    fn shell(&mut self, shell: Shape, closed: bool) -> Result<usize> {
        let faces: Vec<Shape> = self.model.children(shell).collect();
        let mut ids = Vec::with_capacity(faces.len());
        for f in faces {
            ids.push(self.face(f)?);
        }
        let list = ids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
        let kind = if closed { "CLOSED_SHELL" } else { "OPEN_SHELL" };
        Ok(self.add(format!("{kind}('',({list}))")))
    }

    /// A face (placed, used as `face.orientation` says).
    fn face(&mut self, face: Shape) -> Result<usize> {
        let model = self.model;
        let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
        let surface = model.surface(local).ok_or_else(|| Error::InvalidInput("writing STEP: a face without a surface".into()))?.transformed(&face.placement);
        let (surface_id, flipped) = self.surface(&surface)?;
        let forward = face.orientation != Orientation::Reversed;
        let mut bounds = Vec::new();
        for (i, wire) in model.children_local(local).collect::<Vec<_>>().into_iter().enumerate() {
            // The edges in the order they connect (as STEP's loops have them).
            let mut oriented = Vec::new();
            for e in model.ordered_edges(Shape { placement: Transform::IDENTITY, ..wire }) {
                // Degenerate edges (at poles) are left out, as OCCT does: readers rebuild
                // them where the loop's curves on the surface end apart.
                if model.edge_curve(Shape { orientation: Orientation::Forward, ..e }).is_none() {
                    continue;
                }
                let edge = self.edge(Shape { orientation: Orientation::Forward, placement: face.placement, ..e })?;
                let sense = e.orientation != Orientation::Reversed;
                oriented.push(self.add(format!("ORIENTED_EDGE('',*,*,#{edge},{})", if sense { ".T." } else { ".F." })));
            }
            let list = oriented.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
            let lp = self.add(format!("EDGE_LOOP('',({list}))"));
            let kind = if i == 0 { "FACE_OUTER_BOUND" } else { "FACE_BOUND" };
            bounds.push(self.add(format!("{kind}('',#{lp},.T.)")));
        }
        let list = bounds.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
        // The face's sense against its written surface: the model's (forward or reversed),
        // turned again where the written surface's normal runs the other way.
        let same_sense = forward != flipped;
        let mark = self.marks.get(&face.node).cloned().unwrap_or_default();
        let name = mark.name.as_deref().map(|n| self.escape(n)).unwrap_or_default();
        let id = self.add(format!("ADVANCED_FACE('{name}',({list}),#{surface_id},{})", if same_sense { ".T." } else { ".F." }));
        if let Some([r, g, b, a]) = mark.color {
            let styled = self.style(id, [r, g, b], Some(a));
            self.face_styles.push(styled);
        }
        if let Some(layer) = &mark.layer {
            self.layer(layer, id);
        }
        Ok(id)
    }

    fn vertex(&mut self, vertex: Shape) -> usize {
        let k = (vertex.node, key(&vertex.placement));
        if let Some(&v) = self.vertices.get(&k) {
            return v;
        }
        let p = self.model.point(vertex);
        let point = self.point(p);
        let v = self.add(format!("VERTEX_POINT('',#{point})"));
        self.vertices.insert(k, v);
        v
    }

    fn edge(&mut self, edge: Shape) -> Result<usize> {
        let k = (edge.node, key(&edge.placement));
        if let Some(&e) = self.edges.get(&k) {
            return Ok(e);
        }
        let model = self.model;
        let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..edge });
        let (va, vb) = (self.vertex(a), self.vertex(b));
        let (curve, first, last) = model
            .edge_curve(Shape { orientation: Orientation::Forward, ..edge })
            .ok_or_else(|| Error::InvalidInput("writing STEP: an edge without a curve".into()))?;
        let curve = self.curve(&curve, first, last)?;
        let mark = self.marks.get(&edge.node).cloned().unwrap_or_default();
        let name = mark.name.as_deref().map(|n| self.escape(n)).unwrap_or_default();
        let e = self.add(format!("EDGE_CURVE('{name}',#{va},#{vb},#{curve},.T.)"));
        if let Some(layer) = &mark.layer {
            self.layer(layer, e);
        }
        self.edges.insert(k, e);
        Ok(e)
    }

    fn point(&mut self, p: Point) -> usize {
        self.add(format!("CARTESIAN_POINT('',({},{},{}))", real(p.x), real(p.y), real(p.z)))
    }

    fn direction(&mut self, d: Direction) -> usize {
        self.add(format!("DIRECTION('',({},{},{}))", real(d.x()), real(d.y()), real(d.z())))
    }

    /// `AXIS2_PLACEMENT_3D` of a frame, right-handed: the main direction is `x × y`.
    fn axis(&mut self, f: &Frame) -> usize {
        let z = Direction::from_vector(Vector::from(f.x()).cross(Vector::from(f.y())));
        let (o, zd, xd) = (self.point(f.origin()), self.direction(z), self.direction(f.x()));
        self.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{zd},#{xd})"))
    }

    /// A curve, its parameter increasing from the edge's start to its end over `[first,
    /// last]`.
    fn curve(&mut self, curve: &Curve, first: f64, last: f64) -> Result<usize> {
        Ok(match curve {
            Curve::Line(l) => {
                let p = self.point(l.origin);
                let d = self.direction(l.direction);
                let v = self.add(format!("VECTOR('',#{d},1.)"));
                self.add(format!("LINE('',#{p},#{v})"))
            }
            Curve::Circle(c) => {
                let a = self.axis(c.frame());
                self.add(format!("CIRCLE('',#{a},{})", real(c.radius())))
            }
            Curve::Ellipse(e) => {
                let a = self.axis(e.frame());
                self.add(format!("ELLIPSE('',#{a},{},{})", real(e.major_radius()), real(e.minor_radius())))
            }
            Curve::BSpline(b) if !b.is_periodic() => self.bspline_curve(b.degree(), b.poles(), b.weights(), b.knots(), b.multiplicities()),
            Curve::Bezier(b) => {
                let s = b.to_bspline();
                self.bspline_curve(s.degree(), s.poles(), s.weights(), s.knots(), s.multiplicities())
            }
            Curve::Trimmed(t) => self.curve(t.basis(), first, last)?,
            // Others (parabolas, hyperbolas, offsets, periodic B-splines): interpolated through
            // points of the edge.
            other => {
                let n = 64;
                let points: Vec<Point> = (0..=n).map(|i| other.point(first + (last - first) * i as f64 / n as f64)).collect();
                let b = crate::geometry::fit::interpolate(&points, false, 0.0)
                    .map_err(|e| Error::InvalidInput(format!("writing STEP: approximating a curve: {e}")))?;
                self.bspline_curve(b.degree(), b.poles(), b.weights(), b.knots(), b.multiplicities())
            }
        })
    }

    fn bspline_curve(&mut self, degree: usize, poles: &[Point], weights: Option<&[f64]>, knots: &[f64], multiplicities: &[usize]) -> usize {
        let ids: Vec<usize> = poles.iter().map(|&p| self.point(p)).collect();
        let list = ids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
        let mults = multiplicities.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",");
        let ks = knots.iter().map(|&k| real(k)).collect::<Vec<_>>().join(",");
        match weights {
            None => self.add(format!("B_SPLINE_CURVE_WITH_KNOTS('',{degree},({list}),.UNSPECIFIED.,.F.,.F.,({mults}),({ks}),.UNSPECIFIED.)")),
            Some(w) => {
                let ws = w.iter().map(|&x| real(x)).collect::<Vec<_>>().join(",");
                self.add(format!(
                    "( BOUNDED_CURVE() B_SPLINE_CURVE({degree},({list}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({mults}),({ks}),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({ws})) REPRESENTATION_ITEM('') )"
                ))
            }
        }
    }

    /// A surface; with whether its written normal runs against the model's `∂u × ∂v` (a
    /// reflected elementary surface, written right-handed).
    fn surface(&mut self, surface: &Surface) -> Result<(usize, bool)> {
        let left_handed = |f: &Frame| Vector::from(f.x()).cross(Vector::from(f.y())).dot(Vector::from(f.z())) < 0.0;
        Ok(match surface {
            Surface::Plane(p) => {
                let a = self.axis(p.frame());
                (self.add(format!("PLANE('',#{a})")), false)
            }
            Surface::Cylinder(c) => {
                let a = self.axis(c.frame());
                (self.add(format!("CYLINDRICAL_SURFACE('',#{a},{})", real(c.radius()))), left_handed(c.frame()))
            }
            Surface::Cone(c) => {
                // STEP's cones widen along their axis: a narrowing one is written about the
                // reversed axis (the same surface; its normal still away from the axis).
                let f = c.frame();
                let z = Vector::from(f.x()).cross(Vector::from(f.y()));
                let angle = if left_handed(f) { -c.semi_angle() } else { c.semi_angle() };
                let axis = if angle < 0.0 { -z } else { z };
                let (o, zd, xd) = (self.point(f.origin()), self.direction(Direction::from_vector(axis)), self.direction(f.x()));
                let a = self.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{zd},#{xd})"));
                (self.add(format!("CONICAL_SURFACE('',#{a},{},{})", real(c.radius()), real(angle.abs()))), left_handed(f))
            }
            Surface::Sphere(s) => {
                let a = self.axis(s.frame());
                (self.add(format!("SPHERICAL_SURFACE('',#{a},{})", real(s.radius()))), left_handed(s.frame()))
            }
            Surface::Torus(t) => {
                let a = self.axis(t.frame());
                (self.add(format!("TOROIDAL_SURFACE('',#{a},{},{})", real(t.major_radius()), real(t.minor_radius()))), left_handed(t.frame()))
            }
            Surface::BSpline(b) if !b.u_knot_vector().periodic && !b.v_knot_vector().periodic => {
                let (u, v) = (b.u_knot_vector(), b.v_knot_vector());
                (self.bspline_surface(u.degree, v.degree, b.poles(), b.weights(), (&u.knots, &u.multiplicities), (&v.knots, &v.multiplicities)), false)
            }
            Surface::Bezier(_) => {
                let b = crate::geometry::convert::surface_to_bspline(surface).ok_or_else(|| Error::InvalidInput("writing STEP: a Bezier surface".into()))?;
                let (u, v) = (b.u_knot_vector().clone(), b.v_knot_vector().clone());
                (self.bspline_surface(u.degree, v.degree, b.poles(), b.weights(), (&u.knots, &u.multiplicities), (&v.knots, &v.multiplicities)), false)
            }
            Surface::Revolution(r) => {
                let (lo, hi) = (r.basis().first_parameter(), r.basis().last_parameter());
                let c = self.curve(r.basis(), lo, hi)?;
                let axis = r.axis();
                let (o, d) = (self.point(axis.origin), self.direction(axis.direction));
                let a = self.add(format!("AXIS1_PLACEMENT('',#{o},#{d})"));
                (self.add(format!("SURFACE_OF_REVOLUTION('',#{c},#{a})")), false)
            }
            Surface::Extrusion(e) => {
                let (lo, hi) = (e.basis().first_parameter(), e.basis().last_parameter());
                let c = self.curve(e.basis(), lo, hi)?;
                let d = self.direction(e.direction());
                let v = self.add(format!("VECTOR('',#{d},1.)"));
                (self.add(format!("SURFACE_OF_LINEAR_EXTRUSION('',#{c},#{v})")), false)
            }
            Surface::Offset(o) => {
                let (basis, flipped) = self.surface(o.basis())?;
                let distance = if flipped { -o.offset() } else { o.offset() };
                (self.add(format!("OFFSET_SURFACE('',#{basis},{},.F.)", real(distance))), flipped)
            }
            Surface::Trimmed(t) => self.surface(t.basis())?,
            other => return Err(Error::InvalidInput(format!("writing STEP: a {:?} surface", std::mem::discriminant(other)))),
        })
    }

    fn bspline_surface(
        &mut self,
        u_degree: usize,
        v_degree: usize,
        poles: &[Vec<Point>],
        weights: Option<&[Vec<f64>]>,
        (uk, um): (&[f64], &[usize]),
        (vk, vm): (&[f64], &[usize]),
    ) -> usize {
        let rows: Vec<String> = poles
            .iter()
            .map(|row| {
                let ids: Vec<usize> = row.iter().map(|&p| self.point(p)).collect();
                format!("({})", ids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(","))
            })
            .collect();
        let grid = rows.join(",");
        let join_u = |m: &[usize]| m.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let join_r = |k: &[f64]| k.iter().map(|&x| real(x)).collect::<Vec<_>>().join(",");
        match weights {
            None => self.add(format!(
                "B_SPLINE_SURFACE_WITH_KNOTS('',{u_degree},{v_degree},({grid}),.UNSPECIFIED.,.F.,.F.,.F.,({}),({}),({}),({}),.UNSPECIFIED.)",
                join_u(um),
                join_u(vm),
                join_r(uk),
                join_r(vk)
            )),
            Some(w) => {
                let ws = w.iter().map(|row| format!("({})", join_r(row))).collect::<Vec<_>>().join(",");
                self.add(format!(
                    "( BOUNDED_SURFACE() B_SPLINE_SURFACE({u_degree},{v_degree},({grid}),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS(({}),({}),({}),({}),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({ws})) REPRESENTATION_ITEM('') SURFACE() )",
                    join_u(um),
                    join_u(vm),
                    join_r(uk),
                    join_r(vk)
                ))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Cylinder};

    /// A unit box written as inches or metres reads back (in millimetres) scaled to match;
    /// the file declares the unit.
    #[test]
    fn units_other_than_millimetres() {
        for (unit, mm, declared) in [(LengthUnit::Inch, 25.4, "CONVERSION_BASED_UNIT('INCH'"), (LengthUnit::Metre, 1000.0, "SI_UNIT($,.METRE.)")] {
            let mut model = Model::new();
            let cube = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
            let mut bytes = Vec::new();
            write_in(&model, cube, "cube", None, unit, &mut bytes).unwrap();
            let text = String::from_utf8(bytes.clone()).unwrap();
            assert!(text.contains(declared), "{unit:?} declared");
            let mut read_model = Model::new();
            let read = super::super::read(&mut read_model, &mut bytes.as_slice()).unwrap();
            let volume = crate::properties::volume(&read_model, read).mass;
            assert!((volume - mm * mm * mm).abs() <= 1e-9 * mm * mm * mm, "{unit:?}: volume {volume}");
        }
    }

    /// A box and a cylinder, written and read back: same volumes, names and colours.
    #[test]
    fn written_assemblies_read_back() {
        let mut model = Model::new();
        let cube = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let rod = Cylinder::new(0.5, 2.0).build(&mut model).unwrap();
        let roots = [Node {
            name: "top".into(),
            shape: None,
            color: None,
            placement: Transform::IDENTITY,
            children: vec![
                Node { name: "cube".into(), shape: Some(cube), color: Some([1.0, 0.0, 0.0]), placement: Transform::translation(Vector::new(5.0, 0.0, 0.0)), ..Default::default() },
                Node { name: "rod".into(), shape: Some(rod), color: Some([0.0, 0.0, 1.0]), ..Default::default() },
            ],
            ..Default::default()
        }];
        let mut bytes = Vec::new();
        write_assembly(&model, &roots, &mut bytes).unwrap();
        let mut other = Model::new();
        let back = super::super::read_assembly(&mut other, &mut bytes.as_slice()).unwrap();
        assert_eq!(back.len(), 1);
        let children = &back[0].children;
        assert_eq!(children.iter().map(|c| c.name.as_str()).collect::<Vec<_>>(), ["cube", "rod"]);
        assert_eq!(children[0].color, Some([1.0, 0.0, 0.0]));
        let volume = |n: &Node, m: &Model| crate::properties::volume(m, n.shape.unwrap()).mass;
        assert!((volume(&children[0], &other) - 6.0).abs() < 1e-9);
        assert!((volume(&children[1], &other) - 0.5 * std::f64::consts::PI).abs() < 1e-6);
        assert!(crate::validity::check(&other, children[1].shape.unwrap()).is_valid());
        assert!((Point::ORIGIN.transformed(&children[0].placement).x - 5.0).abs() < 1e-12);
    }

    /// A face's name, colour and layer, an edge's layer, a part's translucent colour and
    /// its material: written (named colours by name, as OCCT writes them) and read back.
    #[test]
    fn subshapes_and_materials_read_back() {
        use super::super::{Material, Subshape};
        let mut model = Model::new();
        let cube = Cuboid::new(1.0, 2.0, 3.0).build(&mut model).unwrap();
        let top = model.unique_subshapes(cube, NodeKind::Face).into_iter().max_by(|a, b| {
            let z = |f: &Shape| crate::properties::area(&model, *f).center_of_mass.z;
            z(a).total_cmp(&z(b))
        }).unwrap();
        let edge = model.unique_subshapes(cube, NodeKind::Edge)[0];
        let material = Material { name: "copper".into(), description: "a metal".into(), density: 8.96, density_unit: "g/cm^3".into() };
        let roots = [Node {
            name: "top".into(),
            children: vec![Node {
                name: "cube".into(),
                shape: Some(cube),
                color: Some([0.2, 0.4, 0.6]),
                alpha: Some(0.5),
                subshapes: vec![
                    Subshape { shape: top, name: Some("lid".into()), color: Some([1.0, 0.0, 0.0, 1.0]), layer: Some("lids".into()) },
                    Subshape { shape: edge, name: None, color: None, layer: Some("edges".into()) },
                ],
                material: Some(material.clone()),
                ..Default::default()
            }],
            ..Default::default()
        }];
        let mut bytes = Vec::new();
        write_assembly(&model, &roots, &mut bytes).unwrap();
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert!(text.contains("ADVANCED_FACE('lid'") && text.contains("DRAUGHTING_PRE_DEFINED_COLOUR('red')"));
        assert!(text.contains("PRESENTATION_LAYER_ASSIGNMENT('lids','visible'") && text.contains("MATERIAL_DESIGNATION('copper'"));
        let mut other = Model::new();
        let back = super::super::read_assembly(&mut other, &mut bytes.as_slice()).unwrap();
        let cube = &back[0].children[0];
        assert_eq!(cube.color, Some([0.2, 0.4, 0.6]));
        assert!((cube.alpha.unwrap() - 0.5).abs() < 1e-12);
        assert_eq!(cube.material.as_ref(), Some(&material));
        let lid = cube.subshapes.iter().find(|s| s.name.as_deref() == Some("lid")).unwrap();
        assert_eq!((lid.shape.kind(), lid.color, lid.layer.as_deref()), (NodeKind::Face, Some([1.0, 0.0, 0.0, 1.0]), Some("lids")));
        let height = crate::properties::area(&model, top).center_of_mass.z;
        assert!((crate::properties::area(&other, lid.shape).center_of_mass.z - height).abs() < 1e-9);
        let edges: Vec<&Subshape> = cube.subshapes.iter().filter(|s| s.layer.as_deref() == Some("edges")).collect();
        assert_eq!(edges.len(), 1);
        assert_eq!(edges[0].shape.kind(), NodeKind::Edge);
    }
}
