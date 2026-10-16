//! A minimal DRAW-like test session.
//!
//! OCCT's DRAW dispatches each Tcl command to a C++ function taking `argc/argv`. We mirror
//! that: every DRAW command is a method on [`Draw`] taking its (already substituted)
//! arguments as strings and returning the Tcl result string. Command methods live in
//! `commands.rs`; `tools/draw2rs.py` appends shims for commands it has not seen before.

use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use kentron::geometry::surface::Surface;
use kentron::topology::{Model, Shape};

mod commands;
pub mod dispatch;
pub mod tcl;

/// Evaluates DRAW numeric expressions as `Draw::Atof` does: numbers, `+ - * /`, parentheses,
/// numeric DRAW variables (`pi` predefined) and the functions `sin cos tan asin acos atan2
/// sqrt`. Like DRAW, characters after a complete expression are ignored (`136.06,` is 136.06).
struct Expression<'a> {
    text: &'a [u8],
    at: usize,
    vars: &'a HashMap<String, DrawObject>,
}

impl Expression<'_> {
    fn skip_blanks(&mut self) {
        while matches!(self.text.get(self.at), Some(b' ' | b'\t')) {
            self.at += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    /// Sums and differences of products.
    fn sum(&mut self) -> Option<f64> {
        let mut x = self.product()?;
        loop {
            match self.peek() {
                Some(b'+') => {
                    self.at += 1;
                    x += self.product()?;
                }
                Some(b'-') => {
                    self.at += 1;
                    x -= self.product()?;
                }
                _ => return Some(x),
            }
        }
    }

    fn product(&mut self) -> Option<f64> {
        let mut x = self.value()?;
        loop {
            match self.peek() {
                Some(b'*') => {
                    self.at += 1;
                    x *= self.value()?;
                }
                Some(b'/') => {
                    self.at += 1;
                    x /= self.value()?;
                }
                _ => return Some(x),
            }
        }
    }

    fn value(&mut self) -> Option<f64> {
        self.skip_blanks();
        let x = match self.peek()? {
            b'(' => {
                self.at += 1;
                let x = self.sum()?;
                (self.peek() == Some(b')')).then_some(())?;
                self.at += 1;
                x
            }
            b'+' => {
                self.at += 1;
                self.value()?
            }
            b'-' => {
                self.at += 1;
                -self.value()?
            }
            c if c.is_ascii_digit() || c == b'.' => {
                let start = self.at;
                while matches!(self.peek(), Some(c) if c.is_ascii_digit() || c == b'.') {
                    self.at += 1;
                }
                if matches!(self.peek(), Some(b'e' | b'E'))
                    && matches!(self.text.get(self.at + 1), Some(c) if c.is_ascii_digit() || *c == b'+' || *c == b'-')
                {
                    self.at += 2;
                    while matches!(self.peek(), Some(c) if c.is_ascii_digit()) {
                        self.at += 1;
                    }
                }
                std::str::from_utf8(&self.text[start..self.at]).ok()?.parse().ok()?
            }
            c if c.is_ascii_alphabetic() || c == b'_' => {
                let start = self.at;
                while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_') {
                    self.at += 1;
                }
                let name = std::str::from_utf8(&self.text[start..self.at]).ok()?;
                match self.vars.get(name) {
                    Some(DrawObject::Number(x)) => *x,
                    Some(_) => return None,
                    None if name == "pi" => std::f64::consts::PI,
                    None => self.call(name)?,
                }
            }
            _ => return None,
        };
        self.skip_blanks();
        Some(x)
    }

    /// A function call `name(arg, ...)`.
    fn call(&mut self, name: &str) -> Option<f64> {
        self.skip_blanks();
        (self.peek() == Some(b'(')).then_some(())?;
        self.at += 1;
        let mut args = vec![self.sum()?];
        while self.peek() == Some(b',') {
            self.at += 1;
            args.push(self.sum()?);
        }
        (self.peek() == Some(b')')).then_some(())?;
        self.at += 1;
        let x = args[0];
        Some(match (name.to_ascii_lowercase().as_str(), args.len()) {
            ("sin", 1) => x.sin(),
            ("cos", 1) => x.cos(),
            ("tan", 1) => x.tan(),
            ("asin", 1) => x.asin(),
            ("acos", 1) => x.acos(),
            ("sqrt", 1) => x.sqrt(),
            ("atan2", 2) => x.atan2(args[1]),
            _ => return None,
        })
    }
}

/// Fails with `message` unless `count` arguments or more were given.
pub(crate) fn need(args: &[&str], count: usize, usage: &str) -> DrawResult<()> {
    match args.len() >= count {
        true => Ok(()),
        false => Err(DrawError::Failed(format!("usage: {usage}"))),
    }
}

impl From<kentron::Error> for DrawError {
    fn from(e: kentron::Error) -> Self {
        DrawError::Failed(e.to_string())
    }
}

/// An object stored in a DRAW variable (shapes, curves, surfaces, numbers, ...).
///
/// Kernel types get added here as they are implemented. Shapes live in the session's
/// [`Model`].
#[derive(Debug, Clone)]
pub enum DrawObject {
    Number(f64),
    Shape(Shape),
    Curve(kentron::geometry::curve::Curve),
    Curve2(kentron::geometry::curve2::Curve2),
    Surface(Surface),
    Point(kentron::math::Point),
    Point2(kentron::math::Point2),
    /// An operation's history (`savehistory`): the shapes generated from each input.
    History(Vec<(kentron::topology::NodeId, Vec<Shape>)>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum DrawError {
    /// The command has no implementation yet.
    NotImplemented { command: String, args: Vec<String> },
    /// The command ran and reported an error (bad arguments, failed check, ...).
    Failed(String),
}

impl fmt::Display for DrawError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DrawError::NotImplemented { command, args } => {
                write!(f, "DRAW command not implemented: {command} {}", args.join(" "))
            }
            DrawError::Failed(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for DrawError {}

pub type DrawResult<T = String> = Result<T, DrawError>;

/// Runs the body of a test expected to fail (its script is in `expected_failures.txt`; the
/// test is `#[should_panic]`): a failure becomes a panic, so the test passes while it fails
/// and fails ("did not panic") once it passes. With `KENTRON_SKIP_EXPECTED_FAILURES` set,
/// it panics without running.
pub fn expected_failure(test: impl FnOnce() -> DrawResult<()>) {
    if crate::skip_expected_failures() {
        panic!("expected failure, not run (KENTRON_SKIP_EXPECTED_FAILURES)");
    }
    if let Err(e) = test() {
        panic!("expected failure: {e}");
    }
}

#[derive(Debug, Default)]
pub struct Draw {
    vars: HashMap<String, DrawObject>,
    /// Owns every shape the session creates.
    pub model: Model,
    /// The object and tool of the last `bop`, for `bopfuse`, `bopcut`, ...
    pub(crate) bop_operands: Option<(Shape, Shape)>,
    /// The evolving fillet being described by `mkevol` and `updatevol`: the result's name,
    /// the shape and the edges with their `(u, r)` pairs.
    pub(crate) evolving: Option<(String, Shape, Vec<(Shape, Vec<(f64, f64)>)>)>,
    /// The offset being described by `offsetparameter` and `offsetload`: the join, and the
    /// shape, distance and faces left out.
    pub(crate) offset_join: kentron::offset::Join,
    pub(crate) offset_load: Option<(Shape, f64, Vec<Shape>)>,
    /// Faces moved by distances of their own (`offsetonface`).
    pub(crate) offset_on_faces: Vec<(Shape, f64)>,
    /// Whether the offset's faces are merged where they share a surface (`offsetparameter
    /// ... r`).
    pub(crate) offset_merge: bool,
    /// The form features defined by `featprism`, `featrevol` and `featdprism` (by kind:
    /// `prism`, `revol`, `dprism`): base, profile, sweep and mode.
    pub(crate) features: HashMap<String, (Shape, Shape, kentron::feature::FeatureSweep, kentron::feature::FeatureMode)>,
    /// The ribs and slots defined by `featlf` and `featrf` (by kind: `lf`, `rf`): base, wire,
    /// sweep and mode.
    pub(crate) ribs: HashMap<String, (Shape, Shape, kentron::geometry::surface::Plane, kentron::feature::RibSweep, kentron::feature::FeatureMode)>,
    /// The arguments of the general fuse (`baddobjects`, `baddtools`), and their cells once
    /// split (`bfillds` and the commands building on it).
    pub(crate) fuse_objects: Vec<Shape>,
    pub(crate) fuse_tools: Vec<Shape>,
    pub(crate) fuse_cells: Option<Vec<kentron::boolean::Cell>>,
    /// The cells builder's result (`bcadd`, `bcremove`, ...): cells taken, by index, with
    /// their material (0 for none), and the materials whose cells are merged.
    pub(crate) cells_taken: Vec<(usize, i64)>,
    pub(crate) cells_merged: std::collections::HashSet<i64>,
    /// The sweep being described by `mksweep`, `setsweep` and `addsweep`: the spine, the
    /// trihedron mode (`setsweep`'s option, `-CF` by default) and the profiles with their
    /// contact and correction flags.
    pub(crate) sweep_spine: Option<Shape>,
    pub(crate) sweep_mode: Vec<String>,
    pub(crate) sweep_profiles: Vec<(Shape, bool, bool)>,
    /// Whether `bapibop` unifies its result's edges and faces (`bsimplify -e`, `-f`).
    pub(crate) simplify: (bool, bool),
    /// The history of the last operation that keeps one (`savehistory` saves it).
    pub(crate) last_history: Vec<(kentron::topology::NodeId, Vec<Shape>)>,
    /// The variables created since last taken (the Tcl runner makes each a Tcl variable
    /// too, as DRAW does).
    pub(crate) created: Vec<String>,
    /// The solids flagged closed, as OCCT's builders flag them (`thrusections`'s); `whatis`
    /// reports other shapes closed by their topology.
    pub closed_solids: std::collections::HashSet<kentron::topology::NodeId>,
    /// The axis of the helices made (`setaxis`; the world's by default).
    pub(crate) helix_frame: Option<kentron::math::Frame>,
}

impl Draw {
    pub fn new() -> Self {
        Self::default()
    }

    /// A compound of `parts`.
    pub(crate) fn compound_of(&mut self, parts: &[Shape]) -> DrawResult<Shape> {
        let compound = self.model.new_empty(kentron::topology::NodeKind::Compound);
        for &part in parts {
            self.model.add(compound, part)?;
        }
        Ok(compound)
    }

    /// The cells of the general fuse's arguments (objects, then tools), split on first use.
    pub(crate) fn fuse_cells(&mut self) -> DrawResult<&Vec<kentron::boolean::Cell>> {
        if self.fuse_cells.is_none() {
            let shapes: Vec<Shape> = self.fuse_objects.iter().chain(&self.fuse_tools).copied().collect();
            if shapes.is_empty() {
                return Err(DrawError::Failed("no arguments for the general fuse".into()));
            }
            self.fuse_cells = Some(kentron::boolean::cells(&mut self.model, &shapes)?);
        }
        Ok(self.fuse_cells.as_ref().expect("split"))
    }

    /// The cells picked by `result s1 1 s2 0 ... [-m material] [-u]` (inside the shapes
    /// marked 1, outside those marked 0), the material and whether to merge.
    pub(crate) fn cells_selection(&mut self, command: &str, args: &[&str]) -> DrawResult<(Vec<usize>, i64, bool)> {
        if args.is_empty() {
            return Err(DrawError::Failed(format!("{command}: no result name")));
        }
        let arguments: Vec<Shape> = self.fuse_objects.iter().chain(&self.fuse_tools).copied().collect();
        let (mut take, mut avoid) = (Vec::new(), Vec::new());
        let (mut material, mut update) = (0, false);
        let mut k = 1;
        while k < args.len() {
            match args[k] {
                "-m" => {
                    material = args.get(k + 1).and_then(|m| m.parse().ok()).ok_or_else(|| DrawError::Failed(format!("{command}: bad material")))?;
                    k += 2;
                }
                "-u" => {
                    update = true;
                    k += 1;
                }
                name => {
                    let shape = self.shape_var(name)?;
                    let index = arguments.iter().position(|&a| a == shape).ok_or_else(|| DrawError::Failed(format!("{command}: {name} is not an argument")))?;
                    let flag = args.get(k + 1).copied().unwrap_or("1");
                    if flag == "0" { avoid.push(index) } else { take.push(index) }
                    k += 2;
                }
            }
        }
        let selected = self
            .fuse_cells()?
            .iter()
            .enumerate()
            .filter(|(_, c)| take.iter().all(|&i| c.inside[i]) && avoid.iter().all(|&i| !c.inside[i]))
            .map(|(i, _)| i)
            .collect();
        Ok((selected, material, update))
    }

    /// The cells taken, as the variable `name`: each on its own, those of a merged
    /// material fused.
    pub(crate) fn set_cells_result(&mut self, name: &str) -> DrawResult {
        let cells: Vec<Shape> = self.fuse_cells()?.iter().map(|c| c.shape).collect();
        let mut parts = Vec::new();
        let mut merged: HashMap<i64, Shape> = HashMap::new();
        // (Merged, a material's cells lose their boundaries between them: pieces of lower
        // dimension inside its solids go with them.)
        let solid = |model: &kentron::topology::Model, s: Shape| !model.unique_subshapes(s, kentron::topology::NodeKind::Solid).is_empty();
        let has_solids: std::collections::HashSet<i64> =
            self.cells_taken.iter().filter(|&&(i, m)| m != 0 && solid(&self.model, cells[i])).map(|&(_, m)| m).collect();
        for &(i, material) in &self.cells_taken.clone() {
            if material != 0 && self.cells_merged.contains(&material) && has_solids.contains(&material) && !solid(&self.model, cells[i]) {
                continue;
            }
            if material != 0 && self.cells_merged.contains(&material) {
                let fused = match merged.get(&material) {
                    Some(&so_far) => kentron::boolean::union(&mut self.model, so_far, cells[i])?,
                    None => cells[i],
                };
                merged.insert(material, fused);
            } else {
                parts.push(cells[i]);
            }
        }
        let mut materials: Vec<i64> = merged.keys().copied().collect();
        materials.sort_unstable();
        for m in materials {
            let fused = merged[&m];
            parts.extend(self.model.unique_subshapes(fused, kentron::topology::NodeKind::Solid));
        }
        let result = self.compound_of(&parts)?;
        self.set_shape_var(name, result)
    }

    /// The names of the variables set, in no particular order.
    pub fn variable_names(&self) -> Vec<String> {
        self.vars.keys().cloned().collect()
    }

    pub fn get(&self, name: &str) -> Option<&DrawObject> {
        self.vars.get(name)
    }

    pub fn set(&mut self, name: &str, value: DrawObject) {
        if self.vars.insert(name.to_string(), value).is_none() {
            self.created.push(name.to_string());
        }
    }

    /// Whether `name` is a DRAW variable.
    pub fn has(&self, name: &str) -> bool {
        self.vars.contains_key(name)
    }

    /// The variable `name`, removed.
    pub(crate) fn take_var(&mut self, name: &str) -> DrawResult<DrawObject> {
        self.vars.remove(name).ok_or_else(|| DrawError::Failed(format!("{name} is not a DRAW variable")))
    }

    /// Fails unless the DRAW variable `name` exists (models the `isdraw result` checks in
    /// OCCT's group `end` scripts).
    pub fn require_var(&self, name: &str) -> DrawResult<()> {
        match self.vars.contains_key(name) {
            true => Ok(()),
            false => Err(DrawError::Failed(format!("Error : variable {name} was not created"))),
        }
    }

    /// Fails unless `bopcheck` finds the named shape's sub-shapes meeting only where they
    /// share vertices (the end script of OCCT's volumemaker grid).
    pub fn require_no_interference(&mut self, name: &str) -> DrawResult<()> {
        let report = self.bopcheck(&[name])?;
        match report.contains("This shape seems to be OK.") {
            true => Ok(()),
            false => Err(DrawError::Failed(format!("Error: bopcheck failed: {}", report.trim()))),
        }
    }

    /// The shape in variable `name`.
    pub fn shape_var(&self, name: &str) -> DrawResult<Shape> {
        match self.vars.get(name) {
            Some(DrawObject::Shape(s)) => Ok(*s),
            Some(other) => Err(DrawError::Failed(format!("{name} is not a shape: {other:?}"))),
            None => Err(DrawError::Failed(format!("{name} is not a DRAW variable"))),
        }
    }

    /// The 3D curve in variable `name`.
    pub fn curve_var(&self, name: &str) -> DrawResult<kentron::geometry::curve::Curve> {
        match self.vars.get(name) {
            Some(DrawObject::Curve(c)) => Ok(c.clone()),
            Some(other) => Err(DrawError::Failed(format!("{name} is not a 3D curve: {other:?}"))),
            None => Err(DrawError::Failed(format!("{name} is not a DRAW variable"))),
        }
    }

    /// The frame of an analytic curve's arguments: `x y z` (axis Z), `x y z nx ny nz` (X axis
    /// chosen like gp_Ax2) or `x y z nx ny nz ux uy uz`.
    pub(crate) fn curve_frame(v: &[f64]) -> DrawResult<kentron::math::Frame> {
        use kentron::math::{Direction, Frame, Point, Vector};
        let dir = |a: &[f64]| {
            Direction::try_from_vector(Vector::new(a[0], a[1], a[2])).ok_or_else(|| DrawError::Failed("null direction".into()))
        };
        let origin = Point::new(v[0], v[1], v[2]);
        Ok(match v.len() {
            3 => Frame::from_z(origin, Direction::Z),
            6 => Frame::from_z(origin, dir(&v[3..6])?),
            9 => Frame::new(origin, dir(&v[3..6])?, dir(&v[6..9])?),
            _ => return Err(DrawError::Failed("bad frame arguments".into())),
        })
    }

    /// The 2D frame of an analytic curve's arguments: `x y` (axis X) or `x y dx dy`.
    pub(crate) fn curve_frame2(v: &[f64]) -> DrawResult<kentron::math::Frame2> {
        use kentron::math::{Direction2, Frame2, Point2};
        let x = match v.len() {
            2 => Direction2::X,
            4 => Direction2::try_new(v[2], v[3]).ok_or_else(|| DrawError::Failed("null direction".into()))?,
            _ => return Err(DrawError::Failed("bad frame arguments".into())),
        };
        Ok(Frame2::new(Point2::new(v[0], v[1]), x, true))
    }

    /// The surface in variable `name`.
    pub fn surface_var(&self, name: &str) -> DrawResult<Surface> {
        match self.vars.get(name) {
            Some(DrawObject::Surface(s)) => Ok(s.clone()),
            Some(other) => Err(DrawError::Failed(format!("{name} is not a surface: {other:?}"))),
            None => Err(DrawError::Failed(format!("{name} is not a DRAW variable"))),
        }
    }

    /// The value of a numeric DRAW argument: a number or an expression (see
    /// [`Expression`]).
    #[doc(alias = "Draw::Atof")]
    pub fn number(&self, arg: &str) -> DrawResult<f64> {
        let mut e = Expression { text: arg.as_bytes(), at: 0, vars: &self.vars };
        e.sum().ok_or_else(|| DrawError::Failed(format!("{arg} is not a number")))
    }

    /// Kentron's own command (not OCCT's), used where ported scripts name sub-shapes by
    /// OCCT's enumeration order: `nearest var shape E|F|V x y z` binds `var` to the edge, face
    /// or vertex of `shape` whose centroid (of its length, area, or its point) is nearest to
    /// `(x, y, z)`; with `measure`, among those sharing that centroid, the one whose length or
    /// area it is. Fails unless exactly one fits.
    pub fn nearest(&mut self, args: &[&str]) -> DrawResult {
        use kentron::topology::NodeKind;
        if args.len() != 6 && args.len() != 7 {
            return Err(DrawError::Failed("usage: nearest var shape E|F|V x y z [measure]".into()));
        }
        let shape = self.shape_var(args[1])?;
        let kind = match args[2] {
            "E" => NodeKind::Edge,
            "F" => NodeKind::Face,
            "V" => NodeKind::Vertex,
            "W" => NodeKind::Wire,
            "SH" => NodeKind::Shell,
            "SO" => NodeKind::Solid,
            k => return Err(DrawError::Failed(format!("nearest: unknown type {k}"))),
        };
        let v = self.numbers(&args[3..])?;
        let props = |model: &Model, s: Shape| match kind {
            NodeKind::Edge | NodeKind::Wire => kentron::properties::length(model, s),
            NodeKind::Face | NodeKind::Shell => kentron::properties::area(model, s),
            _ => kentron::properties::volume(model, s),
        };
        let measure = |model: &Model, s: Shape| if kind == NodeKind::Vertex { 0.0 } else { props(model, s).mass };
        let target = kentron::math::Point::new(v[0], v[1], v[2]);
        let centroid = |model: &Model, s: Shape| if kind == NodeKind::Vertex { model.point(s) } else { props(model, s).center_of_mass };
        let tolerance = 1e-5 * (1.0 + (target - kentron::math::Point::ORIGIN).length());
        let mut found: Vec<(f64, Shape)> = self.model.unique_subshapes(shape, kind).into_iter().map(|s| (centroid(&self.model, s).distance(target), s)).collect();
        found.sort_by(|a, b| a.0.total_cmp(&b.0));
        if let Some(&m) = v.get(3) {
            // Among those at the centroid, by length or area.
            let close: Vec<Shape> = found.iter().filter(|f| f.0 <= tolerance).map(|f| f.1).collect();
            if close.len() > 1 {
                let fits: Vec<Shape> = close.into_iter().filter(|&s| (measure(&self.model, s) - m).abs() <= 1e-5 * (1.0 + m.abs())).collect();
                found = fits.into_iter().map(|s| (0.0, s)).collect();
            }
        }
        match found[..] {
            [(d, s), ..] if d <= tolerance && found.get(1).is_none_or(|n| n.0 > tolerance) => {
                self.set(args[0], DrawObject::Shape(s));
                Ok(args[0].to_string())
            }
            [(d, _), ..] => Err(DrawError::Failed(format!("nearest: no single {k} of {} at ({}, {}, {}); the closest is {d} away", args[1], v[0], v[1], v[2], k = args[2]))),
            [] => Err(DrawError::Failed(format!("nearest: {} has no {}", args[1], args[2]))),
        }
    }

    /// The values of numeric arguments.
    pub fn numbers(&self, args: &[&str]) -> DrawResult<Vec<f64>> {
        args.iter().map(|a| self.number(a)).collect()
    }

    /// Moves the shapes in the named variables by `t` (DRAW's `ttranslate`, `trotate`, ...).
    pub(crate) fn move_shapes(&mut self, names: &[&str], t: &kentron::math::Transform) -> DrawResult {
        for name in names {
            let shape = self.shape_var(name)?;
            self.set_shape_var(name, shape.moved(t))?;
        }
        Ok(String::new())
    }

    /// The optional leading plane variable of a primitive command (its frame, else the world
    /// frame) and the numeric arguments after it.
    pub(crate) fn primitive_frame(&self, args: &[&str]) -> DrawResult<(kentron::math::Frame, Vec<f64>)> {
        let (frame, rest) = match args.first().and_then(|a| self.get(a)) {
            Some(DrawObject::Surface(Surface::Plane(plane))) => (*plane.frame(), &args[1..]),
            _ => (kentron::math::Frame::WORLD, args),
        };
        let values = rest.iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        Ok((frame, values))
    }

    /// Transforms the curves, surfaces and shapes in the named variables by `t` (DRAW's
    /// `translate`, `rotate`, `pmirror`, ...).
    pub(crate) fn transform_vars(&mut self, names: &[&str], t: &kentron::math::Transform) -> DrawResult {
        use kentron::math::Transformable;
        for name in names {
            let moved = match self.vars.get(*name) {
                Some(DrawObject::Shape(s)) => DrawObject::Shape(s.moved(t)),
                Some(DrawObject::Curve(c)) => DrawObject::Curve(c.transformed(t)),
                Some(DrawObject::Surface(s)) => DrawObject::Surface(s.transformed(t)),
                Some(DrawObject::Point(p)) => DrawObject::Point(p.transformed(t)),
                Some(other) => return Err(DrawError::Failed(format!("{name} cannot be transformed in 3D: {other:?}"))),
                // (A name that is not a variable is passed over, as DRAW does.)
                None => continue,
            };
            self.set(name, moved);
        }
        Ok(String::new())
    }

    /// Transforms the 2D curves in the named variables by `t` (DRAW's `2dtranslate`, ...).
    pub(crate) fn transform_vars2(&mut self, names: &[&str], t: &kentron::math::Transform2) -> DrawResult {
        use kentron::math::Transformable2;
        for name in names {
            let moved = match self.vars.get(*name) {
                Some(DrawObject::Curve2(c)) => DrawObject::Curve2(c.transformed(t)),
                Some(DrawObject::Point2(p)) => DrawObject::Point2(p.transformed(t)),
                Some(other) => return Err(DrawError::Failed(format!("{name} is not a 2D curve: {other:?}"))),
                // (A name that is not a variable is passed over, as DRAW does.)
                None => continue,
            };
            self.set(name, moved);
        }
        Ok(String::new())
    }

    /// The 3D curve of variable `name`: a curve, or the curve of an edge (trimmed to the
    /// edge's range).
    pub(crate) fn curve_or_edge(&self, name: &str) -> DrawResult<kentron::geometry::curve::Curve> {
        use kentron::geometry::curve::{Curve, TrimmedCurve};
        match self.vars.get(name) {
            Some(DrawObject::Shape(s)) if s.kind() == kentron::topology::NodeKind::Edge => {
                let (curve, first, last) =
                    self.model.edge_curve(*s).ok_or_else(|| DrawError::Failed(format!("{name} has no 3D curve")))?;
                Ok(Curve::Trimmed(TrimmedCurve::new(curve, first, last)))
            }
            _ => self.curve_var(name),
        }
    }

    /// The fields OCCT's `dump` lists for a line or conic, by name (Origin, Center, Axis,
    /// XAxis, YAxis, Radii, Focal), for `checkdump`.
    pub(crate) fn dump_fields(&self, name: &str) -> DrawResult<Vec<(&'static str, Vec<f64>)>> {
        use kentron::geometry::curve::Curve;
        use kentron::geometry::curve2::Curve2;
        let d3 = |d: kentron::math::Direction| vec![d.x(), d.y(), d.z()];
        let p3 = |p: kentron::math::Point| vec![p.x, p.y, p.z];
        let d2 = |d: kentron::math::Direction2| vec![d.x(), d.y()];
        let p2 = |p: kentron::math::Point2| vec![p.x, p.y];
        let frame3 = |f: &kentron::math::Frame| vec![("Center", p3(f.origin())), ("Axis", d3(f.z())), ("XAxis", d3(f.x())), ("YAxis", d3(f.y()))];
        let frame2 = |f: &kentron::math::Frame2| vec![("Center", p2(f.origin())), ("Axis", d2(f.x())), ("XAxis", d2(f.x())), ("YAxis", d2(f.y()))];
        let mut fields = match self.vars.get(name) {
            Some(DrawObject::Curve(c)) => match c {
                Curve::Line(l) => vec![("Origin", p3(l.origin)), ("Axis", d3(l.direction))],
                Curve::Circle(c) => [frame3(c.frame()), vec![("Radius", vec![c.radius()])]].concat(),
                Curve::Ellipse(c) => [frame3(c.frame()), vec![("Radii", vec![c.major_radius(), c.minor_radius()])]].concat(),
                Curve::Hyperbola(c) => [frame3(c.frame()), vec![("Radii", vec![c.major_radius(), c.minor_radius()])]].concat(),
                Curve::Parabola(c) => [frame3(c.frame()), vec![("Focal", vec![c.focal()])]].concat(),
                _ => Vec::new(),
            },
            Some(DrawObject::Curve2(c)) => match c {
                Curve2::Line(l) => vec![("Origin", p2(l.origin)), ("Axis", d2(l.direction))],
                Curve2::Circle(c) => [frame2(c.frame()), vec![("Radius", vec![c.radius()])]].concat(),
                Curve2::Ellipse(c) => [frame2(c.frame()), vec![("Radii", vec![c.major_radius(), c.minor_radius()])]].concat(),
                Curve2::Hyperbola(c) => [frame2(c.frame()), vec![("Radii", vec![c.major_radius(), c.minor_radius()])]].concat(),
                Curve2::Parabola(c) => [frame2(c.frame()), vec![("Focal", vec![c.focal()])]].concat(),
                _ => Vec::new(),
            },
            Some(other) => return Err(DrawError::Failed(format!("{name}: no dump for {other:?}"))),
            None => return Err(DrawError::Failed(format!("{name} is not a DRAW variable"))),
        };
        fields.retain(|f| !f.1.is_empty());
        Ok(fields)
    }

    pub(crate) fn set_shape_var(&mut self, name: &str, shape: Shape) -> DrawResult {
        self.set(name, DrawObject::Shape(shape));
        Ok(String::new())
    }

    pub(crate) fn not_implemented(&mut self, command: &str, args: &[&str]) -> DrawResult {
        Err(DrawError::NotImplemented {
            command: command.to_string(),
            args: args.iter().map(|a| a.to_string()).collect(),
        })
    }
}

/// Resolves a test data file the way DRAW's `locate_data_file` does: by file name, searched
/// recursively in the directories listed in `OCCT_TEST_DATA` (path-list syntax), then in
/// `~/occt-data` (OCCT's test dataset, an external download: see `crates/occt-tests/README.md`),
/// `$OCCT_ROOT/data` and `~/OCCT/data`, then in a BRL-CAD source checkout's `db` and
/// `src/libbrep/tests` (`$BRLCAD_SRC`, else `~/brlcad`).
///
/// Panics if the file can't be found.
pub fn data_file(name: &str) -> String {
    find_data_file(name).unwrap_or_else(|| panic!("OCCT test data file not found: {name} (set OCCT_TEST_DATA)"))
}

/// [`data_file`], or `None` if it can't be found.
pub fn find_data_file(name: &str) -> Option<String> {
    static INDEX: OnceLock<HashMap<String, PathBuf>> = OnceLock::new();
    let index = INDEX.get_or_init(|| {
        let mut dirs: Vec<PathBuf> = std::env::var_os("OCCT_TEST_DATA")
            .map(|v| std::env::split_paths(&v).collect())
            .unwrap_or_default();
        let home = std::env::var_os("HOME").map(PathBuf::from);
        if let Some(home) = &home {
            dirs.push(home.join("occt-data"));
        }
        if let Some(root) = std::env::var_os("OCCT_ROOT") {
            dirs.push(Path::new(&root).join("data"));
        }
        if let Some(home) = &home {
            dirs.push(home.join("OCCT").join("data"));
        }
        // BRL-CAD's sample models and test files (a source checkout).
        let brlcad = std::env::var_os("BRLCAD_SRC").map(PathBuf::from).or_else(|| home.as_ref().map(|h| h.join("brlcad")));
        if let Some(root) = brlcad {
            dirs.push(root.join("db"));
            dirs.push(root.join("src").join("libbrep").join("tests"));
        }
        let mut index = HashMap::new();
        for dir in dirs {
            index_dir(&dir, &mut index);
        }
        index
    });
    index.get(name).map(|p| p.to_string_lossy().into_owned())
}

fn index_dir(dir: &Path, index: &mut HashMap<String, PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            index_dir(&path, index);
        } else if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            index.entry(name.to_string()).or_insert(path);
        }
    }
}
