//! Runs test cases written in YAML (see `docs/test-cases.md`) against kentron.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use kentron::geometry::curve::{BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, Parabola};
use kentron::geometry::curve2::{Circle2, Curve2, Ellipse2, Hyperbola2, Line2, Parabola2};
use kentron::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Transformable, Vector};
use kentron::primitives::{Cone, Cuboid, Cylinder, Sphere, Torus, Wedge};
use kentron::topology::{Model, NodeKind, Orientation, Shape};
use kentron::{blend, boolean, convert, make, modify, properties, sweep, validity};
use kentron::geometry::surface::Surface;
use kentron::sketch::{self, Outcome, Step};
use yaml_rust2::{Yaml, YamlLoader};

/// A case read from a file.
#[derive(Debug, Clone)]
pub struct Case {
    pub name: String,
    pub file: PathBuf,
    /// Why the case is not run, if it is not.
    pub skip: Option<String>,
    /// Known to fail for want of a feature.
    pub not_done: bool,
    /// Its own limit on processor time, in seconds (a case known to run long and fail).
    pub timeout: Option<f64>,
    doc: Yaml,
}

/// The cases of a YAML file (one per document).
pub fn load(path: &Path) -> Result<Vec<Case>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let docs = YamlLoader::load_from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut cases = Vec::new();
    for doc in docs {
        if doc.is_null() {
            continue;
        }
        let name = doc["case"].as_str().ok_or_else(|| format!("{}: a case without a name", path.display()))?.to_string();
        let mut skip = doc["skip"].as_str().map(str::to_string);
        if skip.is_none() {
            // OCCT's test data is an external download: without it, the cases using it don't run.
            let missing: Vec<&str> = seq(&doc["build"])
                .into_iter()
                .filter_map(|step| {
                    let load = &step.as_hash()?.values().next()?["load"];
                    load["brep"].as_str().or(load["step"].as_str())
                })
                .filter(|f| crate::draw::find_data_file(f).is_none())
                .collect();
            if !missing.is_empty() {
                skip = Some(format!("test data not found: {} (set OCCT_TEST_DATA or BRLCAD_SRC)", missing.join(", ")));
            }
        }
        let not_done = doc["expect"].as_str() == Some("not_done");
        let timeout = doc["timeout"].as_f64().or_else(|| doc["timeout"].as_i64().map(|t| t as f64));
        cases.push(Case { name, file: path.to_path_buf(), skip, not_done, timeout, doc });
    }
    Ok(cases)
}

/// All case files below `root`, sorted.
pub fn files(root: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else { continue };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "yaml") {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

/// Runs a case: `Err` with the reason it failed.
/// Runs a case's build steps (not its checks) as far as they go: the model and its named
/// shapes, for looking into a case.
pub fn build(case: &Case) -> Result<(Model, Vec<(String, Shape)>), String> {
    let (env, _) = built(case)?;
    let shapes = env.vars.iter().filter_map(|(k, v)| match v { Value::Shape(s) => Some((k.clone(), *s)), _ => None }).collect();
    Ok((env.model, shapes))
}

pub fn run(case: &Case) -> Result<(), String> {
    let (env, failed) = built(case)?;
    if let Some(e) = failed {
        return Err(e);
    }
    for check in seq(&case.doc["checks"]) {
        env.check(check)?;
    }
    Ok(())
}

/// The case's build steps run up to the first one that fails (with its error).
fn built(case: &Case) -> Result<(Env, Option<String>), String> {
    let mut env = Env::default();
    if let Yaml::Hash(params) = &case.doc["params"] {
        for (k, v) in params {
            let name = k.as_str().ok_or("a parameter needs a name")?;
            let value = env.number(v)?;
            env.params.insert(name.to_string(), value);
        }
    }
    for step in seq(&case.doc["build"]) {
        let Yaml::Hash(map) = step else { return Err(format!("a build step is a map: {step:?}")) };
        let (key, value) = map.iter().next().ok_or("an empty build step")?;
        let key = key.as_str().ok_or("a build step's key is a name")?;
        if key == "check" {
            env.check(value)?;
        } else {
            match env.operation(value) {
                Ok(v) => {
                    env.vars.insert(key.to_string(), v);
                }
                Err(e) => return Ok((env, Some(format!("{key}: {e}")))),
            }
        }
    }
    Ok((env, None))
}

fn seq(y: &Yaml) -> Vec<&Yaml> {
    match y {
        Yaml::Array(items) => items.iter().collect(),
        Yaml::BadValue | Yaml::Null => Vec::new(),
        other => vec![other],
    }
}

#[derive(Debug, Clone)]
enum Value {
    #[allow(dead_code)] // for numeric results (none yet)
    Number(f64),
    Shape(Shape),
    Curve(Curve),
    #[allow(dead_code)] // 2D curves are built but not yet used by any operation
    Curve2(Curve2),
    Plane(Frame),
}

type R<T> = Result<T, String>;

#[derive(Default)]
struct Env {
    model: Model,
    params: HashMap<String, f64>,
    vars: HashMap<String, Value>,
}

thread_local! {
    /// Where the last panic on this thread happened.
    static PANIC_AT: std::cell::RefCell<Option<String>> = const { std::cell::RefCell::new(None) };
}

/// A panic hook recording the panic's location for [`guarded`], printing nothing.
pub fn quiet_panic_hook(info: &std::panic::PanicHookInfo) {
    let at = info.location().map(|l| format!("{}:{}", l.file(), l.line()));
    PANIC_AT.with(|p| *p.borrow_mut() = at);
}

/// Panics in kentron (`todo!`) become failures with their message and location.
fn guarded<T>(f: impl FnOnce() -> R<T>) -> R<T> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(r) => r,
        Err(p) => {
            let message = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "panic".into());
            Err(match PANIC_AT.with(|p| p.borrow_mut().take()) {
                Some(at) => format!("{message} ({at})"),
                None => message,
            })
        }
    }
}

fn err(e: kentron::Error) -> String {
    e.to_string()
}

impl Env {
    // --- values ---------------------------------------------------------------------------

    fn number(&self, y: &Yaml) -> R<f64> {
        match y {
            Yaml::Real(s) => s.parse::<f64>().map_err(|e| e.to_string()),
            Yaml::Integer(i) => Ok(*i as f64),
            Yaml::String(s) => match self.vars.get(s) {
                Some(Value::Number(v)) => Ok(*v),
                _ => Expr { text: s.as_bytes(), at: 0, params: &self.params }.parse(),
            },
            other => Err(format!("not a number: {other:?}")),
        }
    }

    fn numbers(&self, y: &Yaml) -> R<Vec<f64>> {
        seq(y).into_iter().map(|v| self.number(v)).collect()
    }

    fn point(&self, y: &Yaml) -> R<Point> {
        match self.numbers(y)?[..] {
            [x, yy, z] => Ok(Point::new(x, yy, z)),
            [x, yy] => Ok(Point::new(x, yy, 0.0)),
            _ => Err(format!("not a point: {y:?}")),
        }
    }

    fn vector(&self, y: &Yaml) -> R<Vector> {
        let p = self.point(y)?;
        Ok(Vector::new(p.x, p.y, p.z))
    }

    fn direction(&self, y: &Yaml) -> R<Direction> {
        if let Some(s) = y.as_str() {
            let named = match s {
                "x" => Some(Direction::X),
                "y" => Some(Direction::Y),
                "z" => Some(Direction::Z),
                "-x" => Some(-Direction::X),
                "-y" => Some(-Direction::Y),
                "-z" => Some(-Direction::Z),
                _ => None,
            };
            if let Some(d) = named {
                return Ok(d);
            }
        }
        Direction::try_from_vector(self.vector(y)?).ok_or_else(|| "a null direction".into())
    }

    fn frame(&self, y: &Yaml) -> R<Frame> {
        match y {
            Yaml::BadValue | Yaml::Null => Ok(Frame::WORLD),
            Yaml::String(name) => match self.vars.get(name) {
                Some(Value::Plane(f)) => Ok(*f),
                _ => Err(format!("{name} is not a plane")),
            },
            _ => {
                let origin = if y["origin"].is_badvalue() { Point::ORIGIN } else { self.point(&y["origin"])? };
                let z = if y["z"].is_badvalue() { Direction::Z } else { self.direction(&y["z"])? };
                Ok(match &y["x"] {
                    Yaml::BadValue => Frame::from_z(origin, z),
                    x => Frame::new(origin, z, self.direction(x)?),
                })
            }
        }
    }

    /// The frame of a curve: `frame`, or `center` (and `normal`).
    fn curve_frame(&self, y: &Yaml) -> R<Frame> {
        if !y["frame"].is_badvalue() {
            return self.frame(&y["frame"]);
        }
        let origin = if y["center"].is_badvalue() { Point::ORIGIN } else { self.point(&y["center"])? };
        let z = if y["normal"].is_badvalue() { Direction::Z } else { self.direction(&y["normal"])? };
        Ok(match &y["x"] {
            Yaml::BadValue => Frame::from_z(origin, z),
            x => Frame::new(origin, z, self.direction(x)?),
        })
    }

    /// A 2D frame from `center` (2 numbers) and an optional `x`.
    fn frame2(&self, y: &Yaml) -> R<Option<Frame2>> {
        let center = self.numbers(&y["center"]).unwrap_or_default();
        if center.len() != 2 {
            return Ok(None);
        }
        let x = match &y["x"] {
            Yaml::BadValue => Direction2::X,
            x => match self.numbers(x)?[..] {
                [a, b] => Direction2::try_new(a, b).ok_or("a null direction")?,
                _ => return Err("a 2D direction has two numbers".into()),
            },
        };
        Ok(Some(Frame2::new(Point2::new(center[0], center[1]), x, true)))
    }

    /// A step of a `profile`: a one-key map.
    fn profile_step(&self, y: &Yaml) -> R<Step> {
        let (key, v) = match y.as_hash().map(|h| h.iter().collect::<Vec<_>>()).as_deref() {
            Some([(k, v)]) => (k.as_str().unwrap_or_default(), *v),
            _ => return Err(format!("a profile step is a one-key map: {y:?}")),
        };
        let two = || -> R<(f64, f64)> {
            match self.numbers(v)?[..] {
                [a, b] => Ok((a, b)),
                _ => Err(format!("{key} takes two numbers")),
            }
        };
        Ok(match key {
            "start" => two().map(|(a, b)| Step::Start(a, b))?,
            "origin" => Step::Origin(self.point(v)?),
            "plane" => Step::Plane { normal: self.vector(&v["normal"])?, x: self.vector(&v["x"])? },
            "on" => match self.get(v)? {
                Value::Plane(f) => Step::Frame(f),
                Value::Shape(f) => match self.model.surface(f) {
                    Some(Surface::Plane(p)) if f.kind() == NodeKind::Face => Step::Frame(*p.transformed(&f.placement).frame()),
                    _ => return Err(format!("{v:?} is not a planar face")),
                },
                _ => return Err(format!("{v:?} is not a plane")),
            },
            "along_x" => Step::AlongX(self.number(v)?),
            "to_x" => Step::ToX(self.number(v)?),
            "along_y" => Step::AlongY(self.number(v)?),
            "to_y" => Step::ToY(self.number(v)?),
            "forward" => Step::Forward(self.number(v)?),
            "forward_to_x" => Step::ForwardToX(self.number(v)?),
            "forward_to_y" => Step::ForwardToY(self.number(v)?),
            "by" => two().map(|(a, b)| Step::By(a, b))?,
            "to" => two().map(|(a, b)| Step::To(a, b))?,
            "turn" => Step::Turn(self.number(v)?),
            "heading" => Step::Heading(self.number(v)?),
            "direction" => two().map(|(a, b)| Step::Direction(a, b))?,
            "arc" => Step::Arc { radius: self.number(&v["radius"])?, angle: self.number(&v["angle"])? },
            _ => return Err(format!("unknown profile step {key}")),
        })
    }

    fn get(&self, y: &Yaml) -> R<Value> {
        let name = y.as_str().ok_or_else(|| format!("not a name: {y:?}"))?;
        self.vars.get(name).cloned().ok_or_else(|| format!("{name} is not bound"))
    }

    fn shape(&self, y: &Yaml) -> R<Shape> {
        match self.get(y)? {
            Value::Shape(s) => Ok(s),
            other => Err(format!("{y:?} is not a shape: {other:?}")),
        }
    }

    /// An operand of a boolean: a shape, or a list of them (as a compound).
    fn operand(&mut self, y: &Yaml) -> R<Shape> {
        let shapes = self.shapes(y)?;
        match shapes[..] {
            [s] => Ok(s),
            _ => {
                let c = self.model.new_empty(NodeKind::Compound);
                for s in shapes {
                    self.model.add(c, s).map_err(err)?;
                }
                Ok(c)
            }
        }
    }

    fn shapes(&self, y: &Yaml) -> R<Vec<Shape>> {
        seq(y).into_iter().map(|v| self.shape(v)).collect()
    }

    fn curve(&self, y: &Yaml) -> R<Curve> {
        match self.get(y)? {
            Value::Curve(c) => Ok(c),
            other => Err(format!("{y:?} is not a curve: {other:?}")),
        }
    }

    fn angle(&self, y: &Yaml) -> R<Option<f64>> {
        if y.is_badvalue() {
            return Ok(None);
        }
        Ok(Some(self.number(y)?.to_radians()))
    }

    // --- operations -----------------------------------------------------------------------

    fn operation(&mut self, y: &Yaml) -> R<Value> {
        let Yaml::Hash(map) = y else {
            // A bare name: an alias.
            return self.get(y);
        };
        let (op, a) = map.iter().next().ok_or("an empty operation")?;
        let op = op.as_str().ok_or("an operation's key is its name")?;
        guarded(|| self.apply(op, a))
    }

    fn apply(&mut self, op: &str, a: &Yaml) -> R<Value> {
        let shape = |s: Shape| Ok(Value::Shape(s));
        match op {
            // primitives
            "box" => {
                let corner = if a["corner"].is_badvalue() { Point::ORIGIN } else { self.point(&a["corner"])? };
                let size = self.vector(&a["size"])?;
                if [size.x, size.y, size.z].iter().any(|d| d.abs() <= kentron::tolerance::LINEAR) {
                    return Err("box: null dimension".into());
                }
                let frame = self.frame(&a["frame"])?;
                let b = if a["frame"].is_badvalue() {
                    Cuboid::from_corners(corner, corner + size)
                } else {
                    Cuboid::new(size.x, size.y, size.z).with_frame(frame)
                };
                shape(b.build(&mut self.model).map_err(err)?)
            }
            "wedge" => {
                let s = self.numbers(&a["size"])?;
                let [dx, dy, dz] = s[..] else { return Err("wedge: size has three numbers".into()) };
                let w = if !a["top_dx"].is_badvalue() {
                    Wedge::new(dx, dy, dz, self.number(&a["top_dx"])?)
                } else {
                    let t = self.numbers(&a["top"])?;
                    let [x1, z1, x2, z2] = t[..] else { return Err("wedge: top has four numbers".into()) };
                    Wedge::with_top(dx, dy, dz, x1, z1, x2, z2)
                };
                shape(w.with_frame(self.frame(&a["frame"])?).build(&mut self.model).map_err(err)?)
            }
            "cylinder" => {
                let mut c = Cylinder::new(self.number(&a["radius"])?, self.number(&a["height"])?).with_frame(self.frame(&a["frame"])?);
                if let Some(angle) = self.angle(&a["angle"])? {
                    c = c.with_angle(angle);
                }
                shape(c.build(&mut self.model).map_err(err)?)
            }
            "cone" => {
                let mut c = Cone::new(self.number(&a["r1"])?, self.number(&a["r2"])?, self.number(&a["height"])?).with_frame(self.frame(&a["frame"])?);
                if let Some(angle) = self.angle(&a["angle"])? {
                    c = c.with_angle(angle);
                }
                shape(c.build(&mut self.model).map_err(err)?)
            }
            "sphere" => {
                let mut s = Sphere::new(self.number(&a["radius"])?).with_frame(self.frame(&a["frame"])?);
                if !a["latitudes"].is_badvalue() {
                    let l = self.numbers(&a["latitudes"])?;
                    s = s.with_latitudes(l[0].to_radians(), l[1].to_radians());
                }
                if let Some(angle) = self.angle(&a["angle"])? {
                    s = s.with_angle(angle);
                }
                shape(s.build(&mut self.model).map_err(err)?)
            }
            "torus" => {
                let mut t = Torus::new(self.number(&a["major"])?, self.number(&a["minor"])?).with_frame(self.frame(&a["frame"])?);
                if !a["tube"].is_badvalue() {
                    let l = self.numbers(&a["tube"])?;
                    t = t.with_tube_angles(l[0].to_radians(), l[1].to_radians());
                }
                if let Some(angle) = self.angle(&a["angle"])? {
                    t = t.with_angle(angle);
                }
                shape(t.build(&mut self.model).map_err(err)?)
            }
            "plane" => Ok(Value::Plane(self.frame(&a["frame"])?)),
            "half_space" => {
                let face = self.shape(&a["face"])?;
                let toward = self.point(&a["toward"])?;
                shape(kentron::primitives::half_space(&mut self.model, face, toward).map_err(err)?)
            }
            // curves
            "circle" | "ellipse" | "parabola" | "hyperbola" => {
                if let Some(f) = self.frame2(a)? {
                    let c = match op {
                        "circle" => Curve2::Circle(Circle2::new(f, self.number(&a["radius"])?)),
                        "ellipse" => Curve2::Ellipse(Ellipse2::new(f, self.number(&a["major"])?, self.number(&a["minor"])?)),
                        "hyperbola" => Curve2::Hyperbola(Hyperbola2::new(f, self.number(&a["major"])?, self.number(&a["minor"])?)),
                        _ => Curve2::Parabola(Parabola2::new(f, self.number(&a["focal"])?)),
                    };
                    return Ok(Value::Curve2(c));
                }
                let f = self.curve_frame(a)?;
                Ok(Value::Curve(match op {
                    "circle" => Curve::Circle(Circle::new(f, self.number(&a["radius"])?)),
                    "ellipse" => Curve::Ellipse(Ellipse::new(f, self.number(&a["major"])?, self.number(&a["minor"])?)),
                    "hyperbola" => Curve::Hyperbola(Hyperbola::new(f, self.number(&a["major"])?, self.number(&a["minor"])?)),
                    _ => Curve::Parabola(Parabola::new(f, self.number(&a["focal"])?)),
                }))
            }
            "line" => {
                let origin = self.numbers(&a["origin"])?;
                let d = self.numbers(&a["direction"])?;
                if origin.len() == 2 {
                    let dir = Direction2::try_new(d[0], d[1]).ok_or("a null direction")?;
                    return Ok(Value::Curve2(Curve2::Line(Line2::new(Point2::new(origin[0], origin[1]), dir))));
                }
                Ok(Value::Curve(Curve::Line(Line::new(self.point(&a["origin"])?, self.direction(&a["direction"])?))))
            }
            "bezier" => {
                let poles: Vec<Point> = seq(&a["poles"]).into_iter().map(|p| self.point(p)).collect::<R<_>>()?;
                let curve = match &a["weights"] {
                    Yaml::BadValue => BezierCurve::new(poles),
                    w => BezierCurve::rational(poles, self.numbers(w)?),
                };
                Ok(Value::Curve(Curve::Bezier(curve)))
            }
            "bspline" => {
                let poles: Vec<Point> = seq(&a["poles"]).into_iter().map(|p| self.point(p)).collect::<R<_>>()?;
                let mut knots = Vec::new();
                let mut mults = Vec::new();
                for k in seq(&a["knots"]) {
                    let v = self.numbers(k)?;
                    knots.push(v[0]);
                    mults.push(v[1] as usize);
                }
                let degree = self.number(&a["degree"])? as usize;
                let mut curve = if a["periodic"].as_bool() == Some(true) {
                    BSplineCurve::new_periodic(poles, knots, mults, degree)
                } else {
                    BSplineCurve::new(poles, knots, mults, degree)
                };
                if !a["weights"].is_badvalue() {
                    curve = curve.with_weights(self.numbers(&a["weights"])?);
                }
                Ok(Value::Curve(Curve::BSpline(curve)))
            }
            // topology
            "vertex" => {
                let p = self.point(&a["at"])?;
                shape(make::vertex(&mut self.model, p))
            }
            "edge" => {
                if !a["curve"].is_badvalue() {
                    let c = self.curve(&a["curve"])?;
                    let e = match &a["range"] {
                        Yaml::BadValue => make::edge_from_curve(&mut self.model, c),
                        r => {
                            let r = self.numbers(r)?;
                            make::edge_from_curve_range(&mut self.model, c, r[0], r[1])
                        }
                    };
                    return shape(e.map_err(err)?);
                }
                let end = |env: &mut Env, y: &Yaml| -> R<Shape> {
                    match y.as_str() {
                        Some(_) => env.shape(y),
                        None => {
                            let p = env.point(y)?;
                            Ok(make::vertex(&mut env.model, p))
                        }
                    }
                };
                let (from, to) = (end(self, &a["from"])?, end(self, &a["to"])?);
                shape(make::edge(&mut self.model, from, to).map_err(err)?)
            }
            "polyline" => {
                let points: Vec<Point> = seq(&a["points"]).into_iter().map(|p| self.point(p)).collect::<R<_>>()?;
                let closed = points.len() > 2 && points[0].distance(*points.last().unwrap()) <= kentron::tolerance::LINEAR;
                let distinct = if closed { &points[..points.len() - 1] } else { &points[..] };
                let vertices: Vec<Shape> = distinct.iter().map(|&p| make::vertex(&mut self.model, p)).collect();
                let count = if closed { vertices.len() } else { vertices.len() - 1 };
                let edges = (0..count)
                    .map(|i| make::edge(&mut self.model, vertices[i], vertices[(i + 1) % vertices.len()]))
                    .collect::<kentron::Result<Vec<_>>>()
                    .map_err(err)?;
                shape(make::wire(&mut self.model, &edges).map_err(err)?)
            }
            "wire" => {
                let mut edges = Vec::new();
                for s in self.shapes(&a["edges"])? {
                    match s.kind() {
                        NodeKind::Edge => edges.push(s),
                        NodeKind::Wire => edges.extend(self.model.explore(s, NodeKind::Edge)),
                        k => return Err(format!("wire: a {k:?} is not an edge or a wire")),
                    }
                }
                shape(make::wire(&mut self.model, &edges).map_err(err)?)
            }
            "face" => {
                let wire = self.shape(&a["wire"])?;
                shape(make::planar_face(&mut self.model, wire).map_err(err)?)
            }
            "compound" => {
                let c = self.model.new_empty(NodeKind::Compound);
                for s in self.shapes(&a["of"])? {
                    self.model.add(c, s).map_err(err)?;
                }
                shape(c)
            }
            "empty" => {
                let kind = match a.as_str() {
                    Some("solid") => NodeKind::Solid,
                    Some("shell") => NodeKind::Shell,
                    Some("compound") => NodeKind::Compound,
                    Some("compsolid") => NodeKind::CompSolid,
                    Some("wire") => NodeKind::Wire,
                    Some("face") => NodeKind::Face,
                    other => return Err(format!("empty: unknown kind {other:?}")),
                };
                shape(self.model.new_empty(kind))
            }
            "add" => {
                let (to, what) = (self.shape(&a["to"])?, self.shape(&a["shape"])?);
                self.model.add(to, what).map_err(err)?;
                shape(to)
            }
            "load" => {
                // A BRep or STEP file of the test data.
                let (file, step) = match (a["brep"].as_str(), a["step"].as_str()) {
                    (Some(f), _) => (f, false),
                    (None, Some(f)) => (f, true),
                    _ => return Err("load: a brep or step file name".into()),
                };
                let path = crate::draw::data_file(file);
                let mut f = std::fs::File::open(&path).map_err(|e| format!("{path}: {e}"))?;
                shape(if step { kentron::io::step::read(&mut self.model, &mut f) } else { kentron::io::brep::read(&mut self.model, &mut f) }.map_err(err)?)
            }
            "select" => self.select(a),
            // transforms
            "translate" | "rotate" | "scale" | "reverse" | "orient" | "copy" => {
                let s = self.shape(&a["shape"])?;
                let moved = match op {
                    "translate" => s.moved(&Transform::translation(self.vector(&a["by"])?)),
                    "rotate" => {
                        let axis = &a["axis"];
                        let axis = Axis::new(self.point(&axis["origin"])?, self.direction(&axis["direction"])?);
                        s.moved(&Transform::rotation(axis, self.number(&a["angle"])?.to_radians()))
                    }
                    "scale" => {
                        let t = Transform::scaling(self.point(&a["center"])?, self.number(&a["factor"])?);
                        modify::transform(&mut self.model, s, &t, Default::default()).map_err(err)?
                    }
                    "reverse" => s.reversed(),
                    "orient" => {
                        let orientation = match a["as"].as_str() {
                            Some("forward") => Orientation::Forward,
                            Some("reversed") => Orientation::Reversed,
                            Some("internal") => Orientation::Internal,
                            Some("external") => Orientation::External,
                            other => return Err(format!("orient: unknown orientation {other:?}")),
                        };
                        Shape { orientation, ..s }
                    }
                    _ => modify::copy(&mut self.model, s, false).map_err(err)?,
                };
                shape(moved)
            }
            // sweeps and conversions
            "extrude" => {
                let p = self.shape(&a["profile"])?;
                let by = self.vector(&a["by"])?;
                shape(sweep::extrude(&mut self.model, p, by).map_err(err)?.shape)
            }
            "revolve" => {
                let p = self.shape(&a["profile"])?;
                let axis = Axis::new(self.point(&a["axis"]["origin"])?, self.direction(&a["axis"]["direction"])?);
                let angle = self.number(&a["angle"])?.to_radians();
                shape(sweep::revolve(&mut self.model, p, axis, angle).map_err(err)?.shape)
            }
            "pipe" => {
                let (p, spine) = (self.shape(&a["profile"])?, self.shape(&a["spine"])?);
                shape(sweep::pipe(&mut self.model, p, spine).map_err(err)?.shape)
            }
            "to_nurbs" => {
                let s = self.shape(&a["shape"])?;
                shape(convert::to_nurbs(&mut self.model, s).map_err(err)?)
            }
            // booleans and fillets
            "union" | "difference" | "intersection" | "section" => {
                let (o, t) = (self.operand(&a["object"])?, self.operand(&a["tool"])?);
                let r = match op {
                    "union" => boolean::union(&mut self.model, o, t),
                    "difference" => boolean::difference(&mut self.model, o, t),
                    "intersection" => boolean::intersection(&mut self.model, o, t),
                    _ => boolean::section(&mut self.model, o, t),
                };
                shape(r.map_err(err)?)
            }
            "sew" => {
                let shapes = seq(&a["shapes"]).into_iter().map(|s| self.shape(s)).collect::<R<Vec<Shape>>>()?;
                let mut options = kentron::sew::SewOptions::default();
                // (A missing key reads as a bad value.)
                if !a["tolerance"].is_null() && !a["tolerance"].is_badvalue() {
                    options.tolerance = self.number(&a["tolerance"])?;
                }
                options.cutting = a["cutting"].as_bool().unwrap_or(true);
                options.non_manifold = a["non_manifold"].as_bool().unwrap_or(false);
                // The other switches (face analysis, floating edges, same-parameter
                // processing, the minimal tolerance) have no counterpart: kentron always
                // analyses faces, keeps floating edges out and reparametrizes pcurves.
                shape(kentron::sew::sew(&mut self.model, &shapes, &options).map_err(err)?)
            }
            "loft" => {
                let sections = seq(&a["sections"]).into_iter().map(|s| self.shape(s)).collect::<R<Vec<Shape>>>()?;
                let options = kentron::loft::LoftOptions {
                    solid: a["solid"].as_bool().unwrap_or(false),
                    ruled: a["ruled"].as_bool().unwrap_or(false),
                    // As OCCT's thrusections matches the sections' edges, unless told.
                    matching: match a["matching"].as_str() {
                        Some("by_length") => kentron::loft::Matching::ByLength,
                        Some("as_given") => kentron::loft::Matching::AsGiven,
                        _ => kentron::loft::Matching::ByRays,
                    },
                };
                shape(kentron::loft::loft(&mut self.model, &sections, &options).map_err(err)?.shape)
            }
            "offset_wire" => {
                let s = self.shape(&a["shape"])?;
                // A distance, or a fraction of the wire's length.
                let distance = if a["length_fraction"].is_badvalue() {
                    self.number(&a["distance"])?
                } else {
                    properties::length(&self.model, s).mass * self.number(&a["length_fraction"])?
                };
                shape(kentron::offset::planar_wire(&mut self.model, s, distance, Default::default()).map_err(err)?)
            }
            "offset_shape" | "thick_solid" => {
                let s = self.shape(&a["shape"])?;
                let join = match a["join"].as_str() {
                    Some("arc") | None => kentron::offset::Join::Arc,
                    Some("intersection") => kentron::offset::Join::Intersection,
                    Some(other) => return Err(format!("{op}: unknown join {other}")),
                };
                let r = if op == "offset_shape" {
                    let distance = self.number(&a["distance"])?;
                    kentron::offset::shape(&mut self.model, s, distance, join)
                } else {
                    let faces = seq(&a["faces"]).into_iter().map(|f| self.shape(f)).collect::<R<Vec<Shape>>>()?;
                    let thickness = self.number(&a["thickness"])?;
                    kentron::offset::thick_solid(&mut self.model, s, &faces, thickness, join)
                };
                shape(r.map_err(err)?)
            }
            // The general fuse: every argument split by the others, the pieces together; a
            // split: the objects' pieces only.
            "partition" | "split" => {
                let (shapes, objects) = if op == "partition" {
                    let all = seq(&a["of"]).into_iter().map(|x| self.shape(x)).collect::<R<Vec<Shape>>>()?;
                    let n = all.len();
                    (all, n)
                } else {
                    let objects = seq(&a["objects"]).into_iter().map(|x| self.shape(x)).collect::<R<Vec<Shape>>>()?;
                    let tools = seq(&a["tools"]).into_iter().map(|x| self.shape(x)).collect::<R<Vec<Shape>>>()?;
                    let n = objects.len();
                    (objects.into_iter().chain(tools).collect(), n)
                };
                let cells = kentron::boolean::cells(&mut self.model, &shapes).map_err(err)?;
                let result = self.model.new_empty(NodeKind::Compound);
                for cell in cells.iter().filter(|c| c.inside[..objects].iter().any(|&i| i)) {
                    self.model.add(result, cell.shape).map_err(err)?;
                }
                shape(result)
            }
            "profile" => {
                let frame = self.frame(&a["frame"])?;
                let steps = seq(&a["steps"]).into_iter().map(|s| self.profile_step(s)).collect::<R<Vec<Step>>>()?;
                let outcome = match a["outcome"].as_str().unwrap_or("face") {
                    "face" => Outcome::Face,
                    "closed_wire" => Outcome::ClosedWire,
                    "open_wire" => Outcome::OpenWire,
                    o => return Err(format!("profile outcome {o}?")),
                };
                shape(sketch::profile(&mut self.model, frame, &steps, outcome).map_err(err)?)
            }
            "filleted_union" | "filleted_difference" => {
                let (o, t) = (self.shape(&a["object"])?, self.shape(&a["tool"])?);
                let radius = self.number(&a["radius"])?;
                let (r, sections) = if op == "filleted_union" {
                    boolean::union_with_sections(&mut self.model, o, t)
                } else {
                    boolean::difference_with_sections(&mut self.model, o, t)
                }
                .map_err(err)?;
                if sections.is_empty() {
                    return shape(r);
                }
                let edges: Vec<(Shape, f64)> = sections.into_iter().map(|e| (e, radius)).collect();
                // Like OCCT's bfuseblend, the boolean's result where the fillet fails.
                shape(blend::fillet_edges(&mut self.model, r, &edges).unwrap_or(r))
            }
            "fillet" => {
                let s = self.shape(&a["shape"])?;
                let mut radii = Vec::new();
                let mut laws = Vec::new();
                for e in seq(&a["edges"]) {
                    if e.as_str().is_some() {
                        radii.push((self.shape(e)?, self.number(&a["radius"])?));
                    } else if !e["law"].is_badvalue() {
                        let law: Vec<(f64, f64)> = seq(&e["law"]).into_iter().map(|p| self.numbers(p).map(|v| (v[0], v[1]))).collect::<R<_>>()?;
                        laws.push((self.shape(&e["edge"])?, law));
                    } else {
                        let r = if e["radius"].is_badvalue() { self.number(&a["radius"])? } else { self.number(&e["radius"])? };
                        radii.push((self.shape(&e["edge"])?, r));
                    }
                }
                if !laws.is_empty() && !radii.is_empty() {
                    return Err("fillet: constant and varying radii together".into());
                }
                let r = if laws.is_empty() {
                    blend::fillet_edges(&mut self.model, s, &radii)
                } else {
                    blend::fillet_evolving(&mut self.model, s, &laws)
                };
                shape(r.map_err(err)?)
            }
            other => Err(format!("unknown operation {other}")),
        }
    }

    /// `select: {of, edge|face|vertex: P, length|area: m}`.
    fn select(&self, a: &Yaml) -> R<Value> {
        let of = self.shape(&a["of"])?;
        let kinds = [
            ("edge", NodeKind::Edge),
            ("face", NodeKind::Face),
            ("vertex", NodeKind::Vertex),
            ("wire", NodeKind::Wire),
            ("shell", NodeKind::Shell),
            ("solid", NodeKind::Solid),
        ];
        let (kind, at) = kinds
            .into_iter()
            .find(|(k, _)| !a[*k].is_badvalue())
            .map(|(k, kind)| (kind, &a[k]))
            .ok_or("select: edge, face, vertex, wire, shell or solid")?;
        let target = self.point(at)?;
        // Of the length of edges and wires, the area of faces and shells, the volume of solids.
        let props = |s: Shape| match kind {
            NodeKind::Edge | NodeKind::Wire => properties::length(&self.model, s),
            NodeKind::Face | NodeKind::Shell => properties::area(&self.model, s),
            _ => properties::volume(&self.model, s),
        };
        let measure_of = |s: Shape| if kind == NodeKind::Vertex { 0.0 } else { props(s).mass };
        let centroid = |s: Shape| if kind == NodeKind::Vertex { self.model.point(s) } else { props(s).center_of_mass };
        let tolerance = 1e-5 * (1.0 + (target - Point::ORIGIN).length());
        let mut close: Vec<Shape> =
            self.model.unique_subshapes(of, kind).into_iter().filter(|&s| centroid(s).distance(target) <= tolerance).collect();
        let measure = [&a["length"], &a["area"], &a["volume"]].into_iter().find(|m| !m.is_badvalue());
        if close.len() > 1 {
            if let Some(m) = measure {
                let m = self.number(m)?;
                close.retain(|&s| (measure_of(s) - m).abs() <= 1e-5 * (1.0 + m.abs()));
            }
        }
        // (A shell or solid, where the shape holds only one, as DRAW's explode takes it: its
        // centroid need not match OCCT's to the digit; the checks after it measure it.)
        if close.is_empty() && matches!(kind, NodeKind::Shell | NodeKind::Solid) {
            if let [only] = self.model.unique_subshapes(of, kind)[..] {
                close.push(only);
            }
        }
        match close[..] {
            [s] => Ok(Value::Shape(s)),
            [] => Err(format!("select: no {kind:?} of {:?} at {target:?}", a["of"])),
            _ => Err(format!("select: {} {kind:?}s of {:?} at {target:?}", close.len(), a["of"])),
        }
    }

    // --- checks ---------------------------------------------------------------------------

    fn check(&self, y: &Yaml) -> R<()> {
        guarded(|| self.check_inner(y))
    }

    fn check_inner(&self, y: &Yaml) -> R<()> {
        if let Some(name) = y["valid"].as_str() {
            let s = self.shape(&y["valid"])?;
            let report = validity::check(&self.model, s);
            return match report.is_valid() {
                true => Ok(()),
                false => Err(format!(
                    "{name} is not valid: {}",
                    report.errors().map(|p| p.description.clone()).collect::<Vec<_>>().join("; ")
                )),
            };
        }
        if let Some(name) = y["exists"].as_str() {
            return match self.vars.contains_key(name) {
                true => Ok(()),
                false => Err(format!("{name} is not bound")),
            };
        }
        if !y["count"].is_badvalue() {
            let s = self.shape(&y["count"])?;
            let (cmp, expected) = ["eq", "lt", "le", "gt", "ge"]
                .into_iter()
                .find_map(|c| match &y[c] {
                    Yaml::Hash(h) => Some((c, h)),
                    _ => None,
                })
                .ok_or("count: eq (or lt, le, gt, ge) is a map of counts")?;
            let kinds = [
                ("vertices", NodeKind::Vertex),
                ("edges", NodeKind::Edge),
                ("wires", NodeKind::Wire),
                ("faces", NodeKind::Face),
                ("shells", NodeKind::Shell),
                ("solids", NodeKind::Solid),
                ("compsolids", NodeKind::CompSolid),
                ("compounds", NodeKind::Compound),
            ];
            // Distinct sub-shapes (one placed twice counts once), or with `placements: true`
            // distinct placed sub-shapes, or with `occurrences: true` each use (an explorer).
            let occurrences = y["occurrences"].as_bool() == Some(true);
            let placements = y["placements"].as_bool() == Some(true);
            let count = |k: NodeKind| {
                if occurrences {
                    self.model.explore(s, k).count()
                } else if placements {
                    self.model.unique_subshapes(s, k).len()
                } else {
                    let mut nodes: Vec<kentron::topology::NodeId> = self.model.explore(s, k).map(|x| x.node).collect();
                    nodes.sort_unstable();
                    nodes.dedup();
                    nodes.len()
                }
            };
            let mut wrong = Vec::new();
            for (k, v) in expected {
                let k = k.as_str().unwrap_or_default();
                let want = self.number(v)? as usize;
                let have = match k {
                    "shapes" => kinds.iter().map(|&(_, kind)| count(kind)).sum(),
                    _ => count(kinds.iter().find(|(n, _)| *n == k).ok_or_else(|| format!("count: unknown kind {k}"))?.1),
                };
                let ok = match cmp {
                    "eq" => have == want,
                    "lt" => have < want,
                    "le" => have <= want,
                    "gt" => have > want,
                    _ => have >= want,
                };
                if !ok {
                    wrong.push(format!("{k}: {have}, expected {cmp} {want}"));
                }
            }
            return match wrong.is_empty() {
                true => Ok(()),
                false => Err(format!("wrong numbers of sub-shapes of {:?}: {}", y["count"], wrong.join(", "))),
            };
        }
        if !y["centroid"].is_badvalue() {
            // Of the volume of solids, else of the area of faces, else of the length of edges.
            let s = self.shape(&y["centroid"])?;
            let props = if self.model.explore(s, NodeKind::Solid).next().is_some() {
                properties::volume(&self.model, s)
            } else if self.model.explore(s, NodeKind::Face).next().is_some() {
                properties::area(&self.model, s)
            } else {
                properties::length(&self.model, s)
            };
            let want = self.point(&y["eq"])?;
            let tol = if y["abs"].is_badvalue() { 1e-7 } else { self.number(&y["abs"])? };
            let c = props.center_of_mass;
            return if (c.x - want.x).abs() <= tol && (c.y - want.y).abs() <= tol && (c.z - want.z).abs() <= tol {
                Ok(())
            } else {
                Err(format!("the centroid is {c:?}, expected {want:?}"))
            };
        }
        if !y["kind"].is_badvalue() {
            let s = self.shape(&y["kind"])?;
            let want = y["is"].as_str().ok_or("kind: is")?;
            let have = format!("{:?}", s.kind()).to_lowercase();
            return if have == want { Ok(()) } else { Err(format!("{:?} is a {have}, not a {want}", y["kind"])) };
        }
        if !y["closed"].is_badvalue() {
            let s = self.shape(&y["closed"])?;
            let shells: Vec<Shape> = if s.kind() == NodeKind::Shell { vec![s] } else { self.model.explore(s, NodeKind::Shell).collect() };
            for shell in &shells {
                let id = kentron::topology::ShellId::try_from(shell.node).map_err(|n| format!("{n:?} is not a shell"))?;
                if !self.model.is_closed_shell(id) {
                    return Err(format!("{:?} has an open shell", y["closed"]));
                }
            }
            return if shells.is_empty() { Err(format!("{:?} has no shell", y["closed"])) } else { Ok(()) };
        }
        if !y["free_vertices"].is_badvalue() {
            // Vertices met an odd number of times along the edges (OCCT's checksection).
            let s = self.shape(&y["free_vertices"])?;
            let mut odd: Vec<kentron::topology::NodeId> = Vec::new();
            for v in self.model.explore(s, NodeKind::Vertex) {
                match odd.iter().position(|&n| n == v.node) {
                    Some(i) => {
                        odd.swap_remove(i);
                    }
                    None => odd.push(v.node),
                }
            }
            return self.compare("number of free vertices", odd.len() as f64, y);
        }
        if !y["faults"].is_badvalue() {
            // Faulty sub-shapes (the shape included), each counted once.
            let faults = |name: &Yaml| -> R<usize> {
                let report = validity::check(&self.model, self.shape(name)?);
                let mut nodes: Vec<kentron::topology::NodeId> = report.errors().map(|p| p.node).collect();
                nodes.sort();
                nodes.dedup();
                Ok(nodes.len())
            };
            let actual = faults(&y["faults"])?;
            if y["source"].is_badvalue() {
                return self.compare("number of faults", actual as f64, y);
            }
            let source = faults(&y["source"])?;
            return match actual <= source {
                true => Ok(()),
                false => Err(format!("{actual} faults, more than the source's {source}")),
            };
        }
        if !y["free_edges"].is_badvalue() {
            // Edges of exactly one face, on closed (or open) chains of such edges: OCCT's
            // checkfreebounds.
            let s = self.shape(&y["free_edges"])?;
            let closed = match y["bounds"].as_str().unwrap_or("closed") {
                "closed" => true,
                "open" => false,
                b => return Err(format!("bounds: {b}?")),
            };
            let mut uses: Vec<(kentron::topology::NodeId, Shape, usize)> = Vec::new();
            for f in self.model.explore(s, NodeKind::Face) {
                for e in self.model.explore(f, NodeKind::Edge) {
                    match uses.iter_mut().find(|u| u.0 == e.node) {
                        Some(u) => u.2 += 1,
                        None => uses.push((e.node, e, 1)),
                    }
                }
            }
            let free: Vec<(kentron::topology::NodeId, kentron::topology::NodeId)> = uses
                .iter()
                .filter(|u| u.2 == 1 && !self.model.is_degenerate(u.1))
                .map(|u| {
                    let (a, b) = self.model.edge_vertices(u.1);
                    (a.node, b.node)
                })
                .collect();
            // A chain is closed when each of its vertices ends an even number of its edges.
            let mut component: Vec<usize> = (0..free.len()).collect();
            fn root(c: &mut [usize], mut i: usize) -> usize {
                while c[i] != i {
                    c[i] = c[c[i]];
                    i = c[i];
                }
                i
            }
            for i in 0..free.len() {
                for j in 0..i {
                    let (a, b) = (free[i], free[j]);
                    if a.0 == b.0 || a.0 == b.1 || a.1 == b.0 || a.1 == b.1 {
                        let (ri, rj) = (root(&mut component, i), root(&mut component, j));
                        component[ri] = rj;
                    }
                }
            }
            let mut count = 0;
            for i in 0..free.len() {
                let r = root(&mut component, i);
                let mut ends: Vec<(kentron::topology::NodeId, usize)> = Vec::new();
                for j in (0..free.len()).filter(|&j| root(&mut component, j) == r) {
                    for v in [free[j].0, free[j].1] {
                        match ends.iter_mut().find(|e| e.0 == v) {
                            Some(e) => e.1 += 1,
                            None => ends.push((v, 1)),
                        }
                    }
                }
                if ends.iter().all(|e| e.1 % 2 == 0) == closed {
                    count += 1;
                }
            }
            return self.compare("number of free edges", count as f64, y);
        }
        for quantity in ["volume", "area", "length", "tolerance"] {
            if y[quantity].is_badvalue() {
                continue;
            }
            let s = self.shape(&y[quantity])?;
            let actual = match quantity {
                "volume" => properties::volume(&self.model, s).mass,
                "area" => properties::area(&self.model, s).mass,
                "length" => properties::length(&self.model, s).mass,
                // The largest of the vertices', edges' and faces' tolerances.
                _ => [NodeKind::Vertex, NodeKind::Edge, NodeKind::Face]
                    .into_iter()
                    .flat_map(|k| self.model.unique_subshapes(s, k))
                    .map(|x| self.model.tolerance(x))
                    .fold(0.0, f64::max),
            };
            return self.compare(quantity, actual, y);
        }
        Err(format!("unknown check {y:?}"))
    }

    fn compare(&self, what: &str, actual: f64, y: &Yaml) -> R<()> {
        // `abs`, `rel` (relative to the expected value), or their sum when both are given.
        let tol = |expected: f64| -> R<f64> {
            let (abs, rel) = (&y["abs"], &y["rel"]);
            Ok(match (abs.is_badvalue(), rel.is_badvalue()) {
                (true, true) => 1e-6 * expected.abs(),
                (false, true) => self.number(abs)?,
                (true, false) => self.number(rel)? * expected.abs(),
                (false, false) => self.number(abs)? + self.number(rel)? * expected.abs(),
            })
        };
        let fail = |cond: &str, v: f64| Err(format!("the {what} is {actual}, expected {cond} {v}"));
        if !y["eq"].is_badvalue() {
            let v = self.number(&y["eq"])?;
            if v == 0.0 && y["abs"].is_badvalue() {
                return if actual == 0.0 { Ok(()) } else { fail("exactly", v) };
            }
            return if (actual - v).abs() <= tol(v)? { Ok(()) } else { fail("", v) };
        }
        for (key, ok) in [
            ("lt", (|a: f64, v: f64| a < v) as fn(f64, f64) -> bool),
            ("le", |a, v| a <= v),
            ("gt", |a, v| a > v),
            ("ge", |a, v| a >= v),
        ] {
            if !y[key].is_badvalue() {
                let v = self.number(&y[key])?;
                return if ok(actual, v) { Ok(()) } else { fail(key, v) };
            }
        }
        Err(format!("{what}: no comparison"))
    }
}

/// Expressions over named parameters: `+ - * /`, `^` or `**`, parentheses, `pi` and a few
/// functions.
struct Expr<'a> {
    text: &'a [u8],
    at: usize,
    params: &'a HashMap<String, f64>,
}

impl Expr<'_> {
    fn parse(&mut self) -> R<f64> {
        let v = self.sum()?;
        self.space();
        if self.at != self.text.len() {
            return Err(format!("bad expression {:?}", String::from_utf8_lossy(self.text)));
        }
        Ok(v)
    }

    fn space(&mut self) {
        while self.at < self.text.len() && self.text[self.at].is_ascii_whitespace() {
            self.at += 1;
        }
    }

    fn eat(&mut self, s: &str) -> bool {
        self.space();
        if self.text[self.at..].starts_with(s.as_bytes()) {
            self.at += s.len();
            true
        } else {
            false
        }
    }

    fn sum(&mut self) -> R<f64> {
        let mut v = self.product()?;
        loop {
            if self.eat("+") {
                v += self.product()?;
            } else if self.eat("-") {
                v -= self.product()?;
            } else {
                return Ok(v);
            }
        }
    }

    fn product(&mut self) -> R<f64> {
        let mut v = self.power()?;
        loop {
            self.space();
            if self.text[self.at..].starts_with(b"**") {
                return Ok(v);
            }
            if self.eat("*") {
                v *= self.power()?;
            } else if self.eat("/") {
                v /= self.power()?;
            } else {
                return Ok(v);
            }
        }
    }

    fn power(&mut self) -> R<f64> {
        let base = self.unary()?;
        if self.eat("^") || self.eat("**") {
            return Ok(base.powf(self.power()?));
        }
        Ok(base)
    }

    fn unary(&mut self) -> R<f64> {
        if self.eat("-") {
            return Ok(-self.unary()?);
        }
        if self.eat("+") {
            return self.unary();
        }
        self.atom()
    }

    fn atom(&mut self) -> R<f64> {
        self.space();
        if self.eat("(") {
            let v = self.sum()?;
            if !self.eat(")") {
                return Err("missing )".into());
            }
            return Ok(v);
        }
        let start = self.at;
        let c = *self.text.get(self.at).ok_or("an expression ends early")?;
        if c.is_ascii_digit() || c == b'.' {
            while self.at < self.text.len() {
                let c = self.text[self.at];
                let exponent = (c == b'e' || c == b'E') && self.at + 1 < self.text.len();
                if c.is_ascii_digit() || c == b'.' {
                    self.at += 1;
                } else if exponent {
                    self.at += 1;
                    if matches!(self.text[self.at], b'+' | b'-') {
                        self.at += 1;
                    }
                } else {
                    break;
                }
            }
            let s = std::str::from_utf8(&self.text[start..self.at]).map_err(|e| e.to_string())?;
            return s.parse::<f64>().map_err(|e| format!("{s}: {e}"));
        }
        while self.at < self.text.len() && (self.text[self.at].is_ascii_alphanumeric() || self.text[self.at] == b'_') {
            self.at += 1;
        }
        let name = std::str::from_utf8(&self.text[start..self.at]).map_err(|e| e.to_string())?;
        if name.is_empty() {
            return Err(format!("bad expression {:?}", String::from_utf8_lossy(self.text)));
        }
        if self.eat("(") {
            let x = self.sum()?;
            if self.eat(",") {
                let y = self.sum()?;
                if !self.eat(")") {
                    return Err("missing )".into());
                }
                return match name {
                    "atan2" => Ok(x.atan2(y)),
                    "pow" => Ok(x.powf(y)),
                    "min" => Ok(x.min(y)),
                    "max" => Ok(x.max(y)),
                    other => Err(format!("unknown function {other}/2")),
                };
            }
            if !self.eat(")") {
                return Err("missing )".into());
            }
            let f: fn(f64) -> f64 = match name {
                "sqrt" => f64::sqrt,
                "sin" => f64::sin,
                "cos" => f64::cos,
                "tan" => f64::tan,
                "asin" => f64::asin,
                "acos" => f64::acos,
                "atan" => f64::atan,
                "abs" => f64::abs,
                "exp" => f64::exp,
                "log" => f64::ln,
                other => return Err(format!("unknown function {other}")),
            };
            return Ok(f(x));
        }
        match name {
            "pi" | "PI" => Ok(std::f64::consts::PI),
            _ => self.params.get(name).copied().ok_or_else(|| format!("unknown parameter {name}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expressions() {
        let params = HashMap::from([("SCALE1".to_string(), 10.0)]);
        let e = |s: &str| Expr { text: s.as_bytes(), at: 0, params: &params }.parse().unwrap();
        assert_eq!(e("1*SCALE1"), 10.0);
        assert_eq!(e("2^3^2"), 512.0);
        assert_eq!(e("1.e-7*SCALE1"), 1e-6);
        assert!((e("cos(pi)") + 1.0).abs() < 1e-15);
        assert_eq!(e("10*10*(10+5)/2"), 750.0);
    }
}
