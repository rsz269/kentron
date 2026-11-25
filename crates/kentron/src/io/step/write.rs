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

use super::parse::{push_real, real};
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
    // (Written as it goes rather than gathered first: a file of large B-splines runs to
    // hundreds of megabytes.)
    let mut out = std::io::BufWriter::with_capacity(1 << 16, out);
    let mut written = || -> std::io::Result<()> {
        out.write_all(b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('kentron model'),'2;1');\n")?;
        writeln!(out, "FILE_NAME('{}','',(''),(''),'kentron','kentron','');", w.escape(roots.first().map_or("", |r| r.name.as_str())))?;
        out.write_all(b"FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n")?;
        for (i, line) in w.lines.iter().enumerate() {
            writeln!(out, "#{} = {};", i + 1, line)?;
        }
        out.write_all(b"ENDSEC;\nEND-ISO-10303-21;\n")?;
        out.flush()
    };
    written().map_err(|e| Error::InvalidInput(format!("writing STEP: {e}")))
}

/// `head` and three reals, closed: `head1.,2.,3.))`.
fn triple(head: &str, v: [f64; 3]) -> String {
    let mut line = String::with_capacity(head.len() + 72);
    line.push_str(head);
    for (i, &x) in v.iter().enumerate() {
        if i > 0 {
            line.push(',');
        }
        push_real(&mut line, x);
    }
    line.push_str("))");
    line
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
    /// The marks of its solids, shells and compounds of their own (styled as wholes).
    body_marks: HashMap<NodeId, Mark>,
    /// Products written (their definitions and representations), by the definition they
    /// were read from: each written once, however often it occurs.
    products: HashMap<usize, (usize, usize)>,
    /// Each shape representation written's body items (to style an occurrence of it).
    bodies_of: HashMap<usize, Vec<usize>>,
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
            body_marks: HashMap::default(),
            products: HashMap::default(),
            bodies_of: HashMap::default(),
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

    /// Entity `id` (a curve's or surface's, written unnamed) named `name` (already escaped),
    /// as OCCT names the geometry of named faces and edges.
    fn name_geometry(&mut self, id: usize, name: &str) {
        if name.is_empty() {
            return;
        }
        let line = &mut self.lines[id - 1];
        if line.contains("REPRESENTATION_ITEM('')") {
            *line = line.replacen("REPRESENTATION_ITEM('')", &format!("REPRESENTATION_ITEM('{name}')"), 1);
        } else if let Some(at) = line.find("(''") {
            line.replace_range(at..at + 3, &format!("('{name}'"));
        }
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
        if let Some(written) = node.definition.and_then(|d| self.products.get(&d)) {
            return Ok(*written);
        }
        let name = self.escape(&node.name);
        let product = self.add(format!("PRODUCT('{name}','{name}','',(#{product_context}))"));
        self.add(format!("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{product}))"));
        let formation = self.add(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let definition = self.add(format!("PRODUCT_DEFINITION('design','',#{formation},#{definition_context})"));
        let definition_shape = self.add(format!("PRODUCT_DEFINITION_SHAPE('','',#{definition})"));
        if let Some(material) = &node.material {
            self.material(material, definition_shape);
        }
        self.validation(&node.validation, definition_shape);
        let context = self.geometric_context();
        let origin = self.axis(&Frame::WORLD);
        // The faces' and edges' own names, colours and layers (a wire's are its edges', a
        // solid's or shell's its faces').
        self.marks.clear();
        self.body_marks.clear();
        for sub in &node.subshapes {
            // (A solid's, shell's or compound's as a whole: on its item; a wire's on its
            // edges.)
            if matches!(sub.shape.kind(), NodeKind::Solid | NodeKind::Shell | NodeKind::Compound) {
                let mark = self.body_marks.entry(sub.shape.node).or_default();
                mark.name = sub.name.clone().or(mark.name.take());
                mark.color = sub.color.or(mark.color);
                mark.layer = sub.layer.clone().or(mark.layer.take());
                continue;
            }
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
                if let Some(layer) = &node.layer {
                    self.layer(layer, item);
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
        self.bodies_of.insert(representation, items[1..items.len() - child_axes.len()].to_vec());
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
            self.validation(&child.occurrence_validation, placement);
            // The occurrence's own colour: its bodies styled in the context of a shape
            // representation of the occurrence, as OCCT's STEPCAFControl_Writer writes a
            // component's colour.
            if let Some([r, g, b, a]) = child.occurrence_color {
                let own = self.add(format!("SHAPE_REPRESENTATION('',(#{origin}),#{context})"));
                self.add(format!("SHAPE_DEFINITION_REPRESENTATION(#{placement},#{own})"));
                let mut styles = Vec::new();
                for item in self.bodies_of.get(&child_representation).cloned().unwrap_or_default() {
                    let usage = self.surface_usage([r, g, b], Some(a));
                    let by_context = self.add(format!("PRESENTATION_STYLE_BY_CONTEXT((#{usage}),#{own})"));
                    styles.push(self.add(format!("STYLED_ITEM('overriding color',(#{by_context}),#{item})")));
                }
                if !styles.is_empty() {
                    let list = styles.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",");
                    self.add(format!("MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('',({list}),#{context})"));
                }
            }
        }
        if let Some(d) = node.definition {
            self.products.insert(d, (definition, representation));
        }
        Ok((definition, representation))
    }

    /// A colour on an item, with its opacity where it is not opaque (written as the
    /// rendering's transparency).
    fn style(&mut self, item: usize, rgb: [f64; 3], alpha: Option<f64>) -> usize {
        let usage = self.surface_usage(rgb, alpha);
        let assignment = self.add(format!("PRESENTATION_STYLE_ASSIGNMENT((#{usage}))"));
        self.add(format!("STYLED_ITEM('color',(#{assignment}),#{item})"))
    }

    /// A surface style usage of a colour, with its opacity where it is not opaque.
    fn surface_usage(&mut self, rgb: [f64; 3], alpha: Option<f64>) -> usize {
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
        self.add(format!("SURFACE_STYLE_USAGE(.BOTH.,#{side})"))
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

    /// Geometric validation properties of the product definition shape `of` (each in a
    /// property of its own, as OCCT writes them; lengths in the file's unit, as written).
    fn validation(&mut self, v: &super::Validation, of: usize) {
        if v == &super::Validation::default() {
            return;
        }
        let context = self.geometric_context();
        // (Units of area and volume: the length unit squared and cubed.)
        let length = self.lines.iter().position(|l| l.contains("LENGTH_UNIT()")).map_or(0, |i| i + 1);
        let mut items: Vec<(&str, String)> = Vec::new();
        if let Some(c) = v.centroid {
            items.push(("centroid", format!("CARTESIAN_POINT('centre point',({},{},{}))", real(c.x), real(c.y), real(c.z))));
        }
        if let Some(volume) = v.volume {
            let element = self.add(format!("DERIVED_UNIT_ELEMENT(#{length},3.)"));
            let unit = self.add(format!("( DERIVED_UNIT((#{element})) VOLUME_UNIT() )"));
            items.push(("volume", format!("MEASURE_REPRESENTATION_ITEM('volume measure',VOLUME_MEASURE({}),#{unit})", real(volume))));
        }
        if let Some(area) = v.area {
            let element = self.add(format!("DERIVED_UNIT_ELEMENT(#{length},2.)"));
            let unit = self.add(format!("( AREA_UNIT() DERIVED_UNIT((#{element})) )"));
            items.push(("area", format!("MEASURE_REPRESENTATION_ITEM('surface area measure',AREA_MEASURE({}),#{unit})", real(area))));
        }
        for (what, item) in items {
            let item = self.add(item);
            let property = self.add(format!("PROPERTY_DEFINITION('geometric validation property','{what}',#{of})"));
            let representation = self.add(format!("REPRESENTATION('{what}',(#{item}),#{context})"));
            self.add(format!("PROPERTY_DEFINITION_REPRESENTATION(#{property},#{representation})"));
        }
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
            self.mark_whole(solid.node, item);
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
        // Loose faces (outside shells, beside them or not): each an open shell of its own in a
        // model of its own, as OCCT's TopoDSToStep_MakeShellBasedSurfaceModel writes a face.
        let in_shells: crate::hashing::FastSet<NodeId> = model.explore(shape, NodeKind::Shell).flat_map(|s| model.explore(s, NodeKind::Face).map(|f| f.node).collect::<Vec<_>>()).collect();
        let loose_faces: Vec<Shape> = model.explore(shape, NodeKind::Face).filter(|f| !in_shells.contains(&f.node)).collect();
        if !loose_faces.is_empty() {
            let faces = loose_faces;
            for f in faces {
                let face = self.face(f)?;
                let shell = self.add(format!("OPEN_SHELL('',(#{face}))"));
                out.push((self.add(format!("SHELL_BASED_SURFACE_MODEL('',(#{shell}))")), true));
            }
        }
        // Loose edges (outside faces, beside them or not): a set of their curves, each trimmed
        // to the edge (at its ends' points, as well as parameters), named as the edges are
        // (the set as the first).
        let mut loose_edges: Vec<Shape> = Vec::new();
        let in_faces: crate::hashing::FastSet<NodeId> = model.explore(shape, NodeKind::Face).flat_map(|f| model.explore(f, NodeKind::Edge).map(|e| e.node).collect::<Vec<_>>()).collect();
        for e in model.explore(shape, NodeKind::Edge) {
            if !in_faces.contains(&e.node) && !loose_edges.iter().any(|l| l.node == e.node) {
                loose_edges.push(e);
            }
        }
        // (A set for each compound of its own (with a colour or layer as a whole) holding
        // loose edges, the innermost one an edge is in; the shape's own set for the rest.)
        let mut owned: Vec<(NodeId, crate::hashing::FastSet<NodeId>)> = Vec::new();
        for c in model.explore_all(shape).filter(|c| c.kind() == NodeKind::Compound) {
            if c.node != shape.node && self.body_marks.contains_key(&c.node) && !owned.iter().any(|(o, _)| *o == c.node) {
                owned.push((c.node, model.explore(c, NodeKind::Edge).map(|e| e.node).collect()));
            }
        }
        owned.sort_by_key(|(_, edges)| edges.len());
        let mut groups: Vec<(NodeId, Vec<Shape>)> = Vec::new();
        for e in loose_edges {
            let owner = owned.iter().find(|(_, edges)| edges.contains(&e.node)).map_or(shape.node, |(c, _)| *c);
            match groups.iter_mut().find(|(c, _)| *c == owner) {
                Some((_, members)) => members.push(e),
                None => groups.push((owner, vec![e])),
            }
        }
        for (compound, loose_edges) in groups {
            let mut ids = Vec::new();
            let mut set_name = String::new();
            for e in loose_edges {
                let Some((curve, first, last)) = model.edge_curve(Shape { orientation: Orientation::Forward, ..e }) else { continue };
                let basis = self.curve(&curve, first, last)?;
                let mark = self.marks.get(&e.node).cloned().unwrap_or_default();
                let name = mark.name.as_deref().map(|n| self.escape(n)).unwrap_or_default();
                self.name_geometry(basis, &name);
                if set_name.is_empty() {
                    set_name = name.clone();
                }
                let (a, b) = (self.point(curve.point(first)), self.point(curve.point(last)));
                // (Parameters too where the curve is written on its own (not interpolated
                // anew, over a range of its own).)
                let own = match &curve {
                    Curve::Trimmed(t) => !matches!(t.basis(), Curve::BSpline(b) if !written_exactly(b, first, last)) && matches!(t.basis(), Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) | Curve::BSpline(_)),
                    Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) => true,
                    Curve::BSpline(b) => written_exactly(b, first, last),
                    _ => false,
                };
                let trims = if own {
                    format!("(#{a},PARAMETER_VALUE({})),(#{b},PARAMETER_VALUE({}))", real(first), real(last))
                } else {
                    format!("(#{a}),(#{b})")
                };
                let trimmed = self.add(format!("TRIMMED_CURVE('{name}',#{basis},{trims},.T.,.CARTESIAN.)"));
                // (The edge's colour and layer on its curve, as on an edge of a face.)
                if let Some([r, g, b, a]) = mark.color {
                    let styled = self.style(trimmed, [r, g, b], Some(a));
                    self.face_styles.push(styled);
                }
                if let Some(layer) = &mark.layer {
                    self.layer(layer, trimmed);
                }
                ids.push(trimmed);
            }
            if !ids.is_empty() {
                let list = ids.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
                let set = self.add(format!("GEOMETRIC_CURVE_SET('{set_name}',({list}))"));
                self.mark_whole(compound, set);
                out.push((set, true));
            }
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
        let id = self.add(format!("{kind}('',({list}))"));
        self.mark_whole(shell.node, id);
        Ok(id)
    }

    /// The colour and layer of a solid, shell or compound of its own (`node`) on the item
    /// written for it.
    fn mark_whole(&mut self, node: NodeId, item: usize) {
        let Some(mark) = self.body_marks.get(&node).cloned() else { return };
        if let Some([r, g, b, a]) = mark.color {
            let styled = self.style(item, [r, g, b], Some(a));
            self.face_styles.push(styled);
        }
        if let Some(layer) = &mark.layer {
            self.layer(layer, item);
        }
    }

    /// A face (placed, used as `face.orientation` says).
    fn face(&mut self, face: Shape) -> Result<usize> {
        let model = self.model;
        let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
        let surface = model.surface(local).ok_or_else(|| Error::InvalidInput("writing STEP: a face without a surface".into()))?.transformed(&face.placement);
        let (surface_id, flipped) = self.surface(&surface, Some(model.face_uv_bounds(local)))?;
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
            // (The loop runs anticlockwise about the surface's normal as the model has it; ISO
            // 10303-42 has loops anticlockwise about the face's normal, so on a face used
            // reversed the bound is taken against its loop, as OCCT writes it.)
            bounds.push(self.add(format!("{kind}('',#{lp},{})", if forward { ".T." } else { ".F." })));
        }
        let list = bounds.iter().map(|v| format!("#{v}")).collect::<Vec<_>>().join(",");
        // The face's sense against its written surface: the model's (forward or reversed),
        // turned again where the written surface's normal runs the other way.
        let same_sense = forward != flipped;
        let mark = self.marks.get(&face.node).cloned().unwrap_or_default();
        let name = mark.name.as_deref().map(|n| self.escape(n)).unwrap_or_default();
        self.name_geometry(surface_id, &name);
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
        self.name_geometry(curve, &name);
        let e = self.add(format!("EDGE_CURVE('{name}',#{va},#{vb},#{curve},.T.)"));
        if let Some([r, g, b, a]) = mark.color {
            let styled = self.style(e, [r, g, b], Some(a));
            self.face_styles.push(styled);
        }
        if let Some(layer) = &mark.layer {
            self.layer(layer, e);
        }
        self.edges.insert(k, e);
        Ok(e)
    }

    fn point(&mut self, p: Point) -> usize {
        self.add(triple("CARTESIAN_POINT('',(", [p.x, p.y, p.z]))
    }

    /// The points' entities, one after another (as [`Self::point`] of each in turn): their
    /// ids. (Formatted on the pool: a large B-spline's poles run to millions of reals.)
    fn points(&mut self, points: &[Point]) -> Vec<usize> {
        const CHUNK: usize = 4096;
        let chunks: Vec<&[Point]> = points.chunks(CHUNK).collect();
        let lines = crate::parallel::map(&chunks, |chunk| chunk.iter().map(|p| triple("CARTESIAN_POINT('',(", [p.x, p.y, p.z])).collect::<Vec<String>>());
        lines.into_iter().flatten().map(|line| self.add(line)).collect()
    }

    fn direction(&mut self, d: Direction) -> usize {
        self.add(triple("DIRECTION('',(", [d.x(), d.y(), d.z()]))
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
        // A periodic B-spline along part of a period: the same curve made non-periodic from
        // the knot before the edge that leaves it the most room, its knots moved by whole
        // periods to hold the edge's parameters. (Its ends then inside the curve's range,
        // not at its ends: there a reader finding the edge's ends on the curve, closed,
        // may take either end for the other, and the edge for the rest of the curve; and
        // over the curve's origin, the edge runs on.)
        if let Curve::BSpline(b) = curve {
            if b.is_periodic() && !whole_period(b, first, last) {
                if let Some((knots, b)) = from_knot_before(b, first, last) {
                    return Ok(self.bspline_curve(b.degree(), b.poles(), b.weights(), &knots, b.multiplicities()));
                }
            }
        }
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
            // (A parabola and a hyperbola as themselves: STEP's are kentron's, a parabola's
            // parameter only scaled (λ(u) = C + f(u²·x + 2u·y)), which nothing written here
            // gives.)
            Curve::Parabola(p) if p.focal() > 0.0 => {
                let a = self.axis(p.frame());
                self.add(format!("PARABOLA('',#{a},{})", real(p.focal())))
            }
            Curve::Hyperbola(h) => {
                let a = self.axis(h.frame());
                self.add(format!("HYPERBOLA('',#{a},{},{})", real(h.major_radius()), real(h.minor_radius())))
            }
            Curve::BSpline(b) if !b.is_periodic() => self.bspline_curve(b.degree(), b.poles(), b.weights(), b.knots(), b.multiplicities()),
            Curve::Bezier(b) => {
                let s = b.to_bspline();
                self.bspline_curve(s.degree(), s.poles(), s.weights(), s.knots(), s.multiplicities())
            }
            // (A periodic one as the same curve made non-periodic, where the edge keeps within
            // one period: STEP's knots clamp.)
            Curve::BSpline(b) if within_period(b, first, last) => {
                let mut b = b.clone();
                b.set_not_periodic();
                self.bspline_curve(b.degree(), b.poles(), b.weights(), b.knots(), b.multiplicities())
            }
            Curve::Trimmed(t) => self.curve(t.basis(), first, last)?,
            // Others (degenerate parabolas, offsets, periodic B-splines over more than a
            // period): interpolated through points of the edge.
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
        let ids = self.points(poles);
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
    /// A surface, for a face over `bounds` (its `(u0, u1, v0, v1)`, where known): a surface
    /// of extrusion of a periodic curve is written over the curve's period the face's range
    /// of `u` keeps within.
    fn surface(&mut self, surface: &Surface, bounds: Option<(f64, f64, f64, f64)>) -> Result<(usize, bool)> {
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
            // (A periodic one as the same surface made non-periodic: STEP's knots clamp.)
            Surface::BSpline(b) => {
                let mut b = b.clone();
                b.set_u_not_periodic();
                b.set_v_not_periodic();
                let (u, v) = (b.u_knot_vector().clone(), b.v_knot_vector().clone());
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
                let whole = (e.basis().first_parameter(), e.basis().last_parameter());
                let (lo, hi) = match (e.basis(), bounds) {
                    (Curve::BSpline(b), Some((u0, u1, _, _))) if b.is_periodic() && u1 - u0 < whole.1 - whole.0 => (u0, u1),
                    _ => whole,
                };
                let c = self.curve(e.basis(), lo, hi)?;
                let d = self.direction(e.direction());
                let v = self.add(format!("VECTOR('',#{d},1.)"));
                (self.add(format!("SURFACE_OF_LINEAR_EXTRUSION('',#{c},#{v})")), false)
            }
            Surface::Offset(o) => {
                let (basis, flipped) = self.surface(o.basis(), bounds)?;
                let distance = if flipped { -o.offset() } else { o.offset() };
                (self.add(format!("OFFSET_SURFACE('',#{basis},{},.F.)", real(distance))), flipped)
            }
            Surface::Trimmed(t) => self.surface(t.basis(), bounds)?,
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
        let all: Vec<Point> = poles.iter().flatten().copied().collect();
        let ids = self.points(&all);
        let mut grid = String::with_capacity(ids.len() * 9);
        let mut next = ids.iter();
        for (i, row) in poles.iter().enumerate() {
            grid.push_str(if i == 0 { "(" } else { ",(" });
            for (j, id) in next.by_ref().take(row.len()).enumerate() {
                if j > 0 {
                    grid.push(',');
                }
                let _ = write!(grid, "#{id}");
            }
            grid.push(')');
        }
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

/// The periodic B-spline `b` made non-periodic from the knot before `[first, last]` that
/// leaves it the most room either side, with its knots (moved by the whole periods that bring
/// the range within them): the same curve at the same parameters over the range. `None`
/// where the range is longer than a period.
fn from_knot_before(b: &crate::geometry::curve::BSplineCurve, first: f64, last: f64) -> Option<(Vec<f64>, crate::geometry::curve::BSplineCurve)> {
    if !b.is_periodic() {
        return None;
    }
    let knots = b.knots();
    let period = knots[knots.len() - 1] - knots[0];
    let (lo, hi) = (first.min(last), first.max(last));
    let slack = 1e-12 * period.abs().max(1.0);
    if !(period > 0.0) || hi - lo > period + slack {
        return None;
    }
    // Each knot of a period, moved by whole periods to the last place not past `lo`; its
    // room: the least of the way from it to `lo` and from `hi` to a period on.
    let (i, start) = (0..knots.len() - 1)
        .map(|i| (i, knots[i] + ((lo + slack - knots[i]) / period).floor() * period))
        .filter(|&(_, start)| hi <= start + period + slack)
        .max_by(|x, y| (lo - x.1).min(x.1 + period - hi).total_cmp(&(lo - y.1).min(y.1 + period - hi)))?;
    let mut curve = b.clone();
    curve.set_origin(i);
    curve.set_not_periodic();
    let shift = start - curve.knots()[0];
    Some((curve.knots().iter().map(|k| k + shift).collect(), curve))
}

/// Whether `[first, last]` is a whole period of the periodic B-spline `b` (to rounding).
fn whole_period(b: &crate::geometry::curve::BSplineCurve, first: f64, last: f64) -> bool {
    let knots = b.knots();
    let period = knots[knots.len() - 1] - knots[0];
    (last - first).abs() >= period - 1e-12 * period.abs().max(1.0)
}

/// Whether a B-spline is written as itself over `[first, last]`, at the same parameters (a
/// periodic one over no more than a period).
fn written_exactly(b: &crate::geometry::curve::BSplineCurve, first: f64, last: f64) -> bool {
    !b.is_periodic() || within_period(b, first, last) || from_knot_before(b, first, last).is_some()
}

/// Whether `[first, last]` keeps within the knots of the periodic B-spline `b` (so that the
/// curve made non-periodic has the edge at the same parameters).
fn within_period(b: &crate::geometry::curve::BSplineCurve, first: f64, last: f64) -> bool {
    let knots = b.knots();
    let (lo, hi) = (knots[0], knots[knots.len() - 1]);
    let slack = 1e-12 * (hi - lo).abs().max(1.0);
    b.is_periodic() && first.min(last) >= lo - slack && first.max(last) <= hi + slack
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Cylinder};

    /// A B-spline surface of more poles than one batch of points formatted together: its
    /// points written one after another row by row, each as it would be alone, and the surface
    /// naming them in that order, row by row.
    #[test]
    fn large_bspline_surfaces_name_their_points_in_order() {
        let model = Model::new();
        let (rows, columns) = (90, 70);
        let poles: Vec<Vec<Point>> = (0..rows).map(|i| (0..columns).map(|j| Point::new(i as f64 * 0.1, j as f64 / 3.0, ((i * j) as f64).sin())).collect()).collect();
        let knots = |n: usize| -> (Vec<f64>, Vec<usize>) {
            let inner = n - 3;
            ((0..=inner).map(|k| k as f64).collect(), (0..=inner).map(|k| if k == 0 || k == inner { 4 } else { 1 }).collect())
        };
        let ((uk, um), (vk, vm)) = (knots(rows), knots(columns));
        let mut w = Writer::new(&model);
        let id = w.bspline_surface(3, 3, &poles, None, (&uk, &um), (&vk, &vm));
        assert_eq!(id, rows * columns + 1);
        for (k, p) in poles.iter().flatten().enumerate() {
            assert_eq!(w.lines[k], triple("CARTESIAN_POINT('',(", [p.x, p.y, p.z]));
        }
        let grid: Vec<String> = (0..rows).map(|i| format!("({})", (0..columns).map(|j| format!("#{}", i * columns + j + 1)).collect::<Vec<_>>().join(","))).collect();
        assert!(w.lines[rows * columns].starts_with(&format!("B_SPLINE_SURFACE_WITH_KNOTS('',3,3,({}),", grid.join(","))));
    }

    /// Every face's loops run anticlockwise about the face's normal (its surface's turned by
    /// `same_sense`), as ISO 10303-42 directs, and the normals point out of the solid: a cube's
    /// planes, read from the file by hand.
    #[test]
    fn loops_anticlockwise_about_face_normals() {
        let mut model = Model::new();
        let cube = Cuboid::new(2.0, 3.0, 4.0).build(&mut model).unwrap();
        loops_checked(&model, cube, true);
        // (And a face used reversed: the square z = 0 seen from below.)
        let corners = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let sides: Vec<Shape> = (0..4).map(|i| crate::make::edge_between(&mut model, corners[i], corners[(i + 1) % 4]).unwrap()).collect();
        let square = crate::make::wire(&mut model, &sides).unwrap();
        let face = crate::make::planar_face(&mut model, square).unwrap();
        let below = Shape { orientation: Orientation::Reversed, ..face };
        loops_checked(&model, below, false);
        loops_checked(&model, face, false);
    }

    /// The loops of `shape`'s faces, written and read from the file by hand, anticlockwise
    /// about the faces' normals (and, `solid`, those pointing out of the cube about (1, 1.5, 2)).
    fn loops_checked(model: &Model, shape: Shape, solid: bool) {
        let mut bytes = Vec::new();
        write(model, shape, "cube", None, &mut bytes).unwrap();
        let text = String::from_utf8(bytes).unwrap().replace('\n', "");
        let entities: HashMap<usize, String> = text
            .split(';')
            .filter_map(|e| {
                let (id, body) = e.trim().strip_prefix('#')?.split_once('=')?;
                Some((id.trim().parse().ok()?, body.trim().to_string()))
            })
            .collect();
        let refs = |id: usize| -> Vec<usize> { entities[&id].split('#').skip(1).map(|r| r.chars().take_while(char::is_ascii_digit).collect::<String>().parse().unwrap()).collect() };
        let flag = |id: usize| entities[&id].trim_end_matches(')').ends_with(".T.");
        let numbers = |id: usize| -> Vec<f64> {
            let body = &entities[&id];
            let inner = &body[body.rfind('(').unwrap() + 1..body.find(')').unwrap()];
            inner.split(',').map(|x| x.trim().parse().unwrap()).collect()
        };
        let point = |vertex: usize| { let c = numbers(refs(vertex)[0]); Point::new(c[0], c[1], c[2]) };
        let centre = Point::new(1.0, 1.5, 2.0);
        let mut faces = 0;
        for (&id, body) in &entities {
            if !body.starts_with("ADVANCED_FACE") {
                continue;
            }
            faces += 1;
            let r = refs(id);
            let (bounds, plane) = (&r[..r.len() - 1], r[r.len() - 1]);
            let axis = refs(refs(plane)[0]);
            let z = numbers(axis[1]);
            let sense = if flag(id) { 1.0 } else { -1.0 };
            let normal = Vector::new(z[0], z[1], z[2]) * sense;
            for &b in bounds {
                let mut points = Vec::new();
                for oe in refs(refs(b)[0]) {
                    let ec = refs(oe)[0];
                    let (start, end) = (refs(ec)[0], refs(ec)[1]);
                    points.push(point(if flag(oe) { start } else { end }));
                }
                if !flag(b) {
                    points.reverse();
                }
                // (The loop's turning about the normal: its vector area.)
                let mut area = Vector::new(0.0, 0.0, 0.0);
                for i in 0..points.len() {
                    area = area + (points[i] - centre).cross(points[(i + 1) % points.len()] - centre);
                }
                assert!(area.dot(normal) > 0.0, "face #{id}: its loop runs clockwise about its normal");
                assert!(!solid || (points[0] - centre).dot(normal) > 0.0, "face #{id}: its normal points into the cube");
            }
        }
        assert_eq!(faces, if solid { 6 } else { 1 });
    }

    /// A periodic B-spline curve is written as itself, made non-periodic (not fitted anew).
    #[test]
    fn periodic_curves_written_exactly() {
        use crate::geometry::curve::BSplineCurve;
        let mut model = Model::new();
        let poles = [(14.0, 0.0), (12.0, 8.0), (4.0, 6.0), (0.0, -2.625), (-4.0, 6.0), (-12.0, 8.0), (-14.0, 0.0), (-12.0, -8.0), (-4.0, -6.0), (0.0, 2.625), (4.0, -6.0), (12.0, -8.0)];
        let curve = BSplineCurve::try_new_periodic(poles.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect(), (0..=12).map(f64::from).collect(), vec![1; 13], 3).unwrap();
        let rim = crate::make::edge_from_curve(&mut model, Curve::BSpline(curve.clone())).unwrap();
        let rim = crate::make::wire(&mut model, &[rim]).unwrap();
        let face = crate::make::planar_face(&mut model, rim).unwrap();
        let mut bytes = Vec::new();
        write(&model, face, "peanut", None, &mut bytes).unwrap();
        let mut read_model = Model::new();
        let read = super::super::read(&mut read_model, &mut bytes.as_slice()).unwrap();
        let edge = read_model.unique_subshapes(read, NodeKind::Edge)[0];
        let (written, first, last) = read_model.edge_curve(edge).unwrap();
        let Curve::BSpline(b) = &written else { panic!("a B-spline: {written:?}") };
        assert_eq!(b.poles().len(), 15, "the twelve poles, three repeated where it closes");
        for i in 0..=48 {
            let t = first + (last - first) * f64::from(i) / 48.0;
            assert!(written.point(t).distance(curve.point(t)) <= 1e-12, "at {t}");
        }
    }

    /// Two prisms on a closed periodic B-spline outline (eleven poles about a circle, its knots
    /// from 3), the second's moved across by 2e-4 and up by a quarter of their height: the
    /// union's sides each run all round the outline below the other's end, their edges
    /// along the outlines' pieces between where the outlines cross, some ending at the
    /// curves' origins (where their sides have seams) and some running across them. Written
    /// with each piece's curve made non-periodic from a knot outside it, and the side's
    /// surface from one outside its range, the union reads back with the same volume and
    /// valid. (Written whole, made non-periodic at their origins, the pieces were read back
    /// with their ends there taken at the other end of the curve, and the union came back
    /// an eighth larger.)
    #[test]
    fn prisms_on_a_periodic_outline_and_the_same_moved_a_hair_written() {
        use crate::geometry::curve::BSplineCurve;
        let poles = [
            (2.0287421860930874, 0.0),
            (1.9082164839169744, 1.2263362702027996),
            (0.870576301772249, 1.9062961946583574),
            (-0.31242701423038965, 2.172977614001407),
            (-1.4013433565211482, 1.6172364448156276),
            (-1.7425155631923277, 0.5116487337106861),
            (-1.627670835096333, -0.4779272789673983),
            (-1.2325251501555514, -1.4224098488836912),
            (-0.34547186371883687, -2.4028095905138143),
            (1.0531481476137061, -2.3060727728523887),
            (1.881460500340129, -1.2091412436522042),
        ];
        let height = 7.87355609995029;
        let prism = |model: &mut Model, (dx, dy, z): (f64, f64, f64)| {
            let curve = BSplineCurve::try_new_periodic(poles.iter().map(|&(x, y)| Point::new(x + dx, y + dy, z)).collect(), (3..=14).map(f64::from).collect(), vec![1; 12], 3).unwrap();
            let rim = crate::make::edge_from_curve(model, Curve::BSpline(curve)).unwrap();
            let rim = crate::make::wire(model, &[rim]).unwrap();
            let face = crate::make::planar_face(model, rim).unwrap();
            crate::sweep::extrude(model, face, crate::math::Vector::new(0.0, 0.0, height)).unwrap().shape
        };
        let mut model = Model::new();
        let a = prism(&mut model, (0.0, 0.0, 0.0));
        let b = prism(&mut model, (-1.1018005558799e-4, 1.8549186073198794e-4, 1.9838604188755724));
        let union = crate::boolean::union(&mut model, a, b).unwrap();
        let before = crate::properties::volume(&model, union).mass;
        let mut bytes = Vec::new();
        write(&model, union, "union", None, &mut bytes).unwrap();
        let mut read_model = Model::new();
        let read = super::super::read(&mut read_model, &mut bytes.as_slice()).unwrap();
        let report = crate::validity::check(&read_model, read);
        assert!(report.problems.is_empty(), "{report:?}");
        let after = crate::properties::volume(&read_model, read).mass;
        assert!((after - before).abs() < 1e-10 * before, "{after} against {before}");
    }

    /// A paraboloid of revolution (a parabola turned about its axis, the basis trimmed to the
    /// edge) is written on the parabola itself, its frame and focal length, not on points of
    /// it fitted anew (which strayed from the tip, a ten-thousandth of the way up a parabola
    /// running on 8 high, by more than its height): the surface and the edge, both on it.
    #[test]
    fn parabolas_written_exactly() {
        use crate::geometry::curve::Parabola;
        let mut model = Model::new();
        let frame = Frame::new(Point::new(0.0, 0.0, 0.0), Direction::Y, Direction::Z);
        let rim = 8f64.sqrt() * 1e-2;
        let side = crate::make::edge_from_curve_range(&mut model, Curve::Parabola(Parabola::new(frame, 2.0)), 0.0, rim).unwrap();
        let top = crate::make::edge_between(&mut model, Point::new(rim, 0.0, 1e-4), Point::new(0.0, 0.0, 1e-4)).unwrap();
        let axis = crate::make::edge_between(&mut model, Point::new(0.0, 0.0, 1e-4), Point::new(0.0, 0.0, 0.0)).unwrap();
        let wire = crate::make::wire(&mut model, &[side, top, axis]).unwrap();
        let section = crate::make::planar_face(&mut model, wire).unwrap();
        let tip = crate::sweep::revolve(&mut model, section, crate::math::Axis::Z, std::f64::consts::TAU).unwrap().shape;
        let mut bytes = Vec::new();
        write(&model, tip, "tip", None, &mut bytes).unwrap();
        let text = String::from_utf8(bytes).unwrap();
        assert!(!text.contains("B_SPLINE_CURVE"), "no curve fitted anew");
        let parabola = text.lines().find(|l| l.contains("PARABOLA(")).expect("a parabola");
        assert!(parabola.trim_end().ends_with(",2.);"), "{parabola}");
        let id = parabola.trim().split(' ').next().unwrap();
        assert!(text.lines().any(|l| l.contains("SURFACE_OF_REVOLUTION('',") && l.contains(&format!("{id},"))), "the surface on it");
    }

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

    /// A ring lofted round through squares, periodic across them (B-spline faces periodic in
    /// one direction): written as the same faces made non-periodic, read back whole.
    #[test]
    fn periodic_bspline_faces_written() {
        use crate::geometry::surface::Surface;
        let mut model = Model::new();
        let square = |model: &mut Model, angle: f64| {
            let turn = Transform::rotation(crate::math::Axis::new(Point::ORIGIN, crate::math::Direction::Z), angle);
            let corners: Vec<Point> = [(9.0, -1.0), (11.0, -1.0), (11.0, 1.0), (9.0, 1.0)].iter().map(|&(x, z)| Point::new(x, 0.0, z).transformed(&turn)).collect();
            let edges: Vec<Shape> = (0..4).map(|i| crate::make::edge_between(model, corners[i], corners[(i + 1) % 4]).unwrap()).collect();
            crate::make::wire(model, &edges).unwrap()
        };
        let mut ring: Vec<Shape> = (0..8).map(|k| square(&mut model, k as f64 * std::f64::consts::TAU / 8.0)).collect();
        ring.push(ring[0]);
        let options = crate::loft::LoftOptions { solid: false, ..Default::default() };
        let shell = crate::loft::loft_periodic(&mut model, &ring, &options).unwrap().shape;
        let shell = model.unique_subshapes(shell, NodeKind::Shell)[0];
        let mut solid = model.new_empty(NodeKind::Solid);
        model.add(solid, shell).unwrap();
        if crate::properties::volume(&model, solid).mass < 0.0 {
            solid = model.new_empty(NodeKind::Solid);
            model.add(solid, shell.reversed()).unwrap();
        }
        let periodic = model.unique_subshapes(solid, NodeKind::Face).into_iter().any(|f| {
            matches!(model.surface(f), Some(Surface::BSpline(b)) if b.u_knot_vector().periodic || b.v_knot_vector().periodic)
        });
        assert!(periodic, "the loft's faces are periodic");
        let volume = crate::properties::volume(&model, solid).mass;
        let mut bytes = Vec::new();
        write_in(&model, solid, "ring", None, LengthUnit::Millimetre, &mut bytes).unwrap();
        let mut other = Model::new();
        let back = super::super::read(&mut other, &mut bytes.as_slice()).unwrap();
        assert!(crate::validity::check(&other, back).is_valid(), "{:?}", crate::validity::check(&other, back));
        let read = crate::properties::volume(&other, back).mass;
        assert!((read - volume).abs() < 1e-6 * volume, "{read} vs {volume}");
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
    /// A named face's surface and a named loose edge's curve (a segment, as CadQuery's
    /// `name_geometries` names them) carry the names too, as OCCT writes them; the segment,
    /// alone in its product, written as a curve set and read back as an edge.
    #[test]
    fn named_geometry_and_loose_edges() {
        use super::super::Subshape;
        let mut model = Model::new();
        let cube = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        let top = model.unique_subshapes(cube, NodeKind::Face)[5];
        let segment = crate::make::edge_between(&mut model, Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)).unwrap();
        let named = |shape: Shape, name: &str| Subshape { shape, name: Some(name.into()), color: None, layer: None };
        let roots = [Node {
            name: "top".into(),
            children: vec![
                Node { name: "cube".into(), shape: Some(cube), subshapes: vec![named(top, "lid")], ..Default::default() },
                Node { name: "rod".into(), shape: Some(segment), subshapes: vec![named(segment, "seg")], ..Default::default() },
            ],
            ..Default::default()
        }];
        let mut bytes = Vec::new();
        write_assembly(&model, &roots, &mut bytes).unwrap();
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert_eq!(text.matches("'lid'").count(), 2, "{text}");
        assert!(text.contains("PLANE('lid'"));
        assert_eq!(text.matches("'seg'").count(), 3);
        assert!(text.contains("LINE('seg'") && text.contains("TRIMMED_CURVE('seg'") && text.contains("GEOMETRIC_CURVE_SET('seg'"));
        let mut other = Model::new();
        let back = super::super::read_assembly(&mut other, &mut bytes.as_slice()).unwrap();
        let rod = back[0].children.iter().find(|c| c.name == "rod").unwrap();
        let edges = other.unique_subshapes(rod.shape.unwrap(), NodeKind::Edge);
        assert_eq!(edges.len(), 1);
        let (a, b) = other.edge_vertices(edges[0]);
        assert!(other.point(a).distance(Point::new(0.0, 0.0, 0.0)) < 1e-9 && other.point(b).distance(Point::new(0.0, 1.0, 0.0)) < 1e-9);
    }

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
