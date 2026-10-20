//! Reading IGES files (version 5.3): curves, surfaces, trimmed and bounded surfaces, and
//! boundary representations (manifold solid B-rep objects), as OCCT's `IGESControl_Reader`
//! transfers them.
//!
//! The file's independent entities are transferred, as a compound: curves as edges (a
//! composite curve's pieces in a wire), points as vertices, surfaces as faces on their whole
//! domain, trimmed and bounded surfaces as faces (each with edges of its own, as OCCT's
//! faces of such entities have), B-rep objects as solids sharing their edges; groups as
//! compounds of their members, subfigure instances as their definitions placed. Lengths are
//! converted to millimetres.
//!
//! Faces' boundaries: on B-spline surfaces the file's curves in the surface's parameters are
//! used where they fit the edges; elsewhere, and where they do not, the edges' curves are
//! projected onto the surface (see [`super::faces`]).

use crate::hashing::FastMap as HashMap;
use std::io::Read;

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Circle, Curve, Ellipse, Line};
use crate::geometry::curve2::{BSplineCurve2, Circle2, Curve2, Ellipse2, Line2};
use crate::geometry::surface::{BSplineSurface, ExtrusionSurface, KnotVector, Plane, RevolutionSurface, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Axis, Direction, Direction2, Frame, Frame2, Point, Point2, Transform, Transformable, Vector};
use crate::topology::{EdgeId, Model, NodeKind, Shape, VertexId};

use super::faces::{align_periods, finish_face, fits, follows, project, CurveSpan, Gaps};

fn invalid(what: impl Into<String>) -> Error {
    Error::InvalidInput(format!("reading IGES: {}", what.into()))
}

/// Reads an IGES file's independent entities as one compound (in millimetres).
#[doc(alias = "IGESControl_Reader")]
pub fn read(model: &mut Model, input: &mut dyn Read) -> Result<Shape> {
    let mut bytes = Vec::new();
    input.read_to_end(&mut bytes).map_err(|e| invalid(e.to_string()))?;
    let file = File::parse(&String::from_utf8_lossy(&bytes))?;
    let mut reader = Reader { file: &file, model, edges_made: HashMap::default(), vertices_made: HashMap::default(), gaps: Gaps::default() };
    let out = reader.model.new_empty(NodeKind::Compound);
    // (Where no entity is marked independent, the B-rep objects at the top of the file's
    // structure: its solids, else its shells, else its faces.)
    let marked = file.entries.iter().any(Entry::is_root);
    let top = [186, 514, 510].into_iter().find(|&k| file.entries.iter().any(|e| e.kind == k));
    for (k, entry) in file.entries.iter().enumerate() {
        if !(entry.is_root() || (!marked && Some(entry.kind) == top)) {
            continue;
        }
        // (An entity that cannot be read is left out, as OCCT's reader leaves it.)
        if let Ok(Some(shape)) = reader.transfer(k) {
            reader.model.add(out, shape)?;
        }
    }
    reader.gaps.apply(reader.model)?;
    if (file.millimetres - 1.0).abs() > 1e-12 {
        let scaled = crate::modify::transform(reader.model, out, &Transform::scaling(Point::ORIGIN, file.millimetres), Default::default())?;
        return Ok(scaled);
    }
    Ok(out)
}

// --- the file -----------------------------------------------------------------------------

/// A parameter of an entity.
#[derive(Debug, Clone, PartialEq)]
enum Value {
    Number(f64),
    Text(String),
    /// Left empty: its default.
    Empty,
}

/// A directory entry.
#[derive(Debug, Clone)]
struct Entry {
    kind: i32,
    /// Index into the parameter data lines, and how many.
    params: Vec<Value>,
    /// The directory entry of its transformation matrix (0 for none).
    transform: usize,
    /// Subordinate switch and use flag of its status number.
    subordinate: u32,
    use_flag: u32,
    form: i32,
}

impl Entry {
    /// Whether it is transferred on its own: independent, neither a definition nor a curve
    /// in some surface's parameters, of a kind that makes a shape.
    fn is_root(&self) -> bool {
        self.subordinate == 0
            && self.use_flag != 2
            && self.use_flag != 5
            && matches!(self.kind, 100 | 102 | 104 | 106 | 110 | 112 | 116 | 120 | 122 | 126 | 128 | 140 | 108 | 143 | 144 | 186 | 190 | 402 | 408 | 514 | 510)
    }
}

struct File {
    entries: Vec<Entry>,
    /// Millimetres per unit of the file.
    millimetres: f64,
    /// The least distance the file tells apart (in its units).
    resolution: f64,
}

impl File {
    fn parse(text: &str) -> Result<File> {
        let mut global = String::new();
        let mut directory: Vec<String> = Vec::new();
        let mut parameter: HashMap<usize, String> = HashMap::default();
        // (Lines ended by carriage returns alone, or none: the 80-column records cut out.)
        let text = text.replace("\r\n", "\n").replace('\r', "\n");
        let records: Vec<String> = text
            .lines()
            .flat_map(|l| {
                let chars: Vec<char> = l.chars().collect();
                if chars.len() > 80 && chars.len() % 80 == 0 { chars.chunks(80).map(|c| c.iter().collect()).collect::<Vec<String>>() } else { vec![l.to_string()] }
            })
            .collect();
        for raw in &records {
            let line: Vec<char> = raw.chars().collect();
            if line.len() < 73 {
                continue;
            }
            let data: String = line[..72].iter().collect();
            match line[72] {
                'G' => global.push_str(&data),
                'D' => directory.push(data),
                'P' => {
                    let seq: String = line[73..line.len().min(80)].iter().collect();
                    if let Ok(n) = seq.trim().parse::<usize>() {
                        parameter.insert(n, line[..64].iter().collect());
                    }
                }
                _ => {}
            }
        }
        // The global section: its delimiters (given as Hollerith strings, else ',' and ';'),
        // and the unit.
        let (delimiter, record) = {
            let g = global.trim_start();
            let first = if let Some(rest) = g.strip_prefix("1H") { rest.chars().next().unwrap_or(',') } else { ',' };
            let after = g.find(first).map_or("", |i| &g[i + first.len_utf8()..]);
            let second = if let Some(rest) = after.trim_start().strip_prefix("1H") { rest.chars().next().unwrap_or(';') } else { ';' };
            (first, second)
        };
        let gp = params(&global, delimiter, record);
        let unit_flag = gp.get(13).and_then(as_number).unwrap_or(2.0) as i32;
        let unit_name = gp.get(14).and_then(as_text).unwrap_or_default().to_ascii_uppercase();
        let millimetres = match unit_flag {
            1 => 25.4,
            2 => 1.0,
            4 => 304.8,
            5 => 1_609_344.0,
            6 => 1000.0,
            7 => 1_000_000.0,
            8 => 0.0254,
            9 => 0.001,
            10 => 10.0,
            11 => 0.0000254,
            3 => match unit_name.trim_start_matches("4H").trim_start_matches("2H") {
                "IN" | "INCH" => 25.4,
                "FT" => 304.8,
                "M" => 1000.0,
                "CM" => 10.0,
                _ => 1.0,
            },
            _ => 1.0,
        };
        let field = |line: &str, k: usize| -> String { line.chars().skip(8 * k).take(8).collect::<String>() };
        let int = |s: &str| s.trim().parse::<i64>().unwrap_or(0);
        let mut entries = Vec::with_capacity(directory.len() / 2);
        for pair in directory.chunks(2) {
            let [a, b] = pair else { break };
            let kind = int(&field(a, 0)) as i32;
            let start = int(&field(a, 1)).max(0) as usize;
            let transform = int(&field(a, 6)).max(0) as usize;
            let status = format!("{:0>8}", field(a, 8).trim());
            let digits: Vec<u32> = status.chars().map(|c| c.to_digit(10).unwrap_or(0)).collect();
            let subordinate = digits.get(2).copied().unwrap_or(0) * 10 + digits.get(3).copied().unwrap_or(0);
            let use_flag = digits.get(4).copied().unwrap_or(0) * 10 + digits.get(5).copied().unwrap_or(0);
            let count = int(&field(b, 3)).max(1) as usize;
            let form = int(&field(b, 4)) as i32;
            let mut text = String::new();
            for k in start..start + count {
                if let Some(l) = parameter.get(&k) {
                    text.push_str(l);
                }
            }
            let mut ps = params(&text, delimiter, record);
            // (The first parameter is the entity's type.)
            if !ps.is_empty() {
                ps.remove(0);
            }
            entries.push(Entry { kind, params: ps, transform, subordinate, use_flag, form });
        }
        let resolution = gp.get(18).and_then(as_number).filter(|&r| r > 0.0).unwrap_or(1e-6);
        Ok(File { entries, millimetres, resolution })
    }

    /// The entry of a directory pointer (the sequence number of its first line).
    fn entry(&self, pointer: usize) -> Result<&Entry> {
        if pointer == 0 || pointer % 2 == 0 {
            return Err(invalid(format!("no entity at directory line {pointer}")));
        }
        self.entries.get((pointer - 1) / 2).ok_or_else(|| invalid(format!("no entity at directory line {pointer}")))
    }
}

fn as_number(v: &Value) -> Option<f64> {
    match v {
        Value::Number(x) => Some(*x),
        _ => None,
    }
}

fn as_text(v: &Value) -> Option<String> {
    match v {
        Value::Text(s) => Some(s.clone()),
        _ => None,
    }
}

/// The parameters of a record: delimited fields up to the record delimiter, Hollerith
/// strings (`nH...`) taken whole, numbers with `D` exponents read as `E`.
fn params(text: &str, delimiter: char, record: char) -> Vec<Value> {
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    let mut field = String::new();
    let mut text_value: Option<String> = None;
    let flush = |field: &mut String, text_value: &mut Option<String>, out: &mut Vec<Value>| {
        let f = field.trim();
        out.push(match text_value.take() {
            Some(t) => Value::Text(t),
            None if f.is_empty() => Value::Empty,
            None => f.replace(['D', 'd'], "E").parse::<f64>().map(Value::Number).unwrap_or_else(|_| Value::Text(f.to_string())),
        });
        field.clear();
    };
    while i < chars.len() {
        let c = chars[i];
        if c == delimiter {
            flush(&mut field, &mut text_value, &mut out);
            i += 1;
            continue;
        }
        if c == record {
            flush(&mut field, &mut text_value, &mut out);
            break;
        }
        if c == 'H' && text_value.is_none() && !field.trim().is_empty() && field.trim().chars().all(|d| d.is_ascii_digit()) {
            let n: usize = field.trim().parse().unwrap_or(0);
            let s: String = chars.iter().skip(i + 1).take(n).collect();
            i += 1 + n;
            text_value = Some(s);
            field.clear();
            continue;
        }
        field.push(c);
        i += 1;
    }
    out
}

// --- transferring --------------------------------------------------------------------------

/// A piece of a curve: the curve over its range.
type Piece = (Curve, f64, f64);

struct Reader<'a> {
    file: &'a File,
    model: &'a mut Model,
    /// B-rep edges made, by (edge list, index), and vertices by (vertex list, index).
    edges_made: HashMap<(usize, usize), (EdgeId, Piece)>,
    vertices_made: HashMap<(usize, usize), VertexId>,
    gaps: Gaps,
}

impl Reader<'_> {
    fn entry(&self, pointer: usize) -> Result<&Entry> {
        self.file.entry(pointer)
    }

    fn real(e: &Entry, k: usize) -> f64 {
        e.params.get(k).and_then(as_number).unwrap_or(0.0)
    }

    fn int(e: &Entry, k: usize) -> i64 {
        Self::real(e, k).round() as i64
    }

    fn pointer(e: &Entry, k: usize) -> usize {
        Self::int(e, k).unsigned_abs() as usize
    }

    /// The placement an entity's transformation matrix (and those it is placed by) gives.
    fn transform_of(&self, e: &Entry) -> Result<Transform> {
        let mut t = Transform::IDENTITY;
        let mut pointer = e.transform;
        let mut depth = 0;
        while pointer != 0 && depth < 32 {
            let m = self.entry(pointer)?;
            if m.kind != 124 {
                break;
            }
            let r = |k: usize| Self::real(m, k);
            let rows = [[r(0), r(1), r(2), r(3)], [r(4), r(5), r(6), r(7)], [r(8), r(9), r(10), r(11)]];
            let this = Transform::from_rows(rows).ok_or_else(|| invalid("a singular transformation matrix"))?;
            t = this * t;
            pointer = m.transform;
            depth += 1;
        }
        Ok(t)
    }

    // --- curves ---

    /// A curve entity's pieces, in model space.
    fn curve(&self, pointer: usize) -> Result<Vec<Piece>> {
        let e = self.entry(pointer)?;
        let t = self.transform_of(e)?;
        let pieces = self.curve_local(e)?;
        Ok(if t == Transform::IDENTITY { pieces } else { pieces.into_iter().map(|(c, a, b)| (c.transformed(&t), a, b)).collect() })
    }

    /// A curve entity's pieces in its own definition space.
    fn curve_local(&self, e: &Entry) -> Result<Vec<Piece>> {
        let r = |k: usize| Self::real(e, k);
        let p3 = |k: usize| Point::new(r(k), r(k + 1), r(k + 2));
        match e.kind {
            110 => {
                let (a, b) = (p3(0), p3(3));
                let d = Direction::try_from_vector(b - a).ok_or_else(|| invalid("a line of no length"))?;
                Ok(vec![(Curve::Line(Line::new(a, d)), 0.0, a.distance(b))])
            }
            100 => {
                let z = r(0);
                let (c, s, f) = (Point::new(r(1), r(2), z), Point::new(r(3), r(4), z), Point::new(r(5), r(6), z));
                let radius = c.distance(s);
                let x = Direction::try_from_vector(s - c).ok_or_else(|| invalid("an arc of no radius"))?;
                let frame = Frame::new(c, Direction::Z, x);
                let mut end = (f.y - c.y).atan2(f.x - c.x) - (s.y - c.y).atan2(s.x - c.x);
                end = end.rem_euclid(std::f64::consts::TAU);
                if end <= 1e-12 {
                    end = std::f64::consts::TAU;
                }
                Ok(vec![(Curve::Circle(Circle::new(frame, radius)), 0.0, end)])
            }
            104 => {
                // An ellipse in its standard position (form 1): A x² + C y² + F = 0.
                let (a, c, f, z) = (r(0), r(2), r(5), r(6));
                if e.form != 1 || a * f >= 0.0 || c * f >= 0.0 {
                    return Err(invalid("conic arcs other than ellipses"));
                }
                let (ra, rb) = ((-f / a).sqrt(), (-f / c).sqrt());
                let (s, g) = (Point::new(r(7), r(8), z), Point::new(r(9), r(10), z));
                let along_x = ra >= rb;
                let frame = if along_x { Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X) } else { Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::Y) };
                let (major, minor) = if along_x { (ra, rb) } else { (rb, ra) };
                let angle = |p: Point| if along_x { (p.y / minor).atan2(p.x / major) } else { (-p.x / minor).atan2(p.y / major) };
                let start = angle(s);
                let mut end = angle(g);
                while end <= start + 1e-12 {
                    end += std::f64::consts::TAU;
                }
                Ok(vec![(Curve::Ellipse(Ellipse::new(frame, major, minor)), start, end)])
            }
            106 => {
                // Piecewise linear (forms 11, 12, 13): its points' segments.
                let (ip, n) = (Self::int(e, 0), Self::int(e, 1).max(0) as usize);
                let points: Vec<Point> = match ip {
                    1 => (0..n).map(|k| Point::new(r(3 + 2 * k), r(4 + 2 * k), r(2))).collect(),
                    2 => (0..n).map(|k| p3(2 + 3 * k)).collect(),
                    _ => (0..n).map(|k| p3(2 + 6 * k)).collect(),
                };
                let mut out = Vec::new();
                for w in points.windows(2) {
                    if let Some(d) = Direction::try_from_vector(w[1] - w[0]) {
                        out.push((Curve::Line(Line::new(w[0], d)), 0.0, w[0].distance(w[1])));
                    }
                }
                Ok(out)
            }
            112 => {
                // A parametric spline: each segment's cubic as a Bézier, joined.
                let n = Self::int(e, 3).max(0) as usize;
                let breaks: Vec<f64> = (0..=n).map(|k| r(4 + k)).collect();
                let base = 5 + n;
                let mut out = Vec::new();
                for seg in 0..n {
                    let at = |axis: usize, k: usize| r(base + 12 * seg + 4 * axis + k);
                    let h = breaks[seg + 1] - breaks[seg];
                    // x(s) = a + b s + c s² + d s³, s from 0 to h: Bézier poles.
                    let pole = |i: usize| {
                        let coord = |axis: usize| {
                            let (a, b, c, d) = (at(axis, 0), at(axis, 1), at(axis, 2), at(axis, 3));
                            match i {
                                0 => a,
                                1 => a + b * h / 3.0,
                                2 => a + 2.0 * b * h / 3.0 + c * h * h / 3.0,
                                _ => a + b * h + c * h * h + d * h * h * h,
                            }
                        };
                        Point::new(coord(0), coord(1), coord(2))
                    };
                    let bezier = BSplineCurve::new((0..4).map(pole).collect(), vec![breaks[seg], breaks[seg + 1]], vec![4, 4], 3);
                    out.push((Curve::BSpline(bezier), breaks[seg], breaks[seg + 1]));
                }
                Ok(out)
            }
            126 => {
                let (k, m) = (Self::int(e, 0).max(0) as usize, Self::int(e, 1).max(0) as usize);
                let flat: Vec<f64> = (0..k + m + 2).map(|i| r(6 + i)).collect();
                let w0 = 6 + k + m + 2;
                let weights: Vec<f64> = (0..=k).map(|i| r(w0 + i)).collect();
                let poles: Vec<Point> = (0..=k).map(|i| p3(w0 + k + 1 + 3 * i)).collect();
                let v = w0 + k + 1 + 3 * (k + 1);
                // (An empty range given: the knots'.)
                let (first, last) = if r(v + 1) > r(v) { (r(v), r(v + 1)) } else { (f64::NEG_INFINITY, f64::INFINITY) };
                let rational = weights.iter().any(|&w| (w - weights[0]).abs() > 1e-15);
                let mut out = Vec::new();
                for (poles, weights, flat) in curve_parts(poles, weights, flat, m)? {
                    let (lo, hi) = (flat[m].max(first), flat[flat.len() - 1 - m].min(last));
                    if !(hi > lo) {
                        continue;
                    }
                    let (knots, mults) = distinct(&flat);
                    let mut curve = BSplineCurve::new(poles, knots, mults, m);
                    if rational {
                        curve = curve.with_weights(weights);
                    }
                    out.push((Curve::BSpline(curve), lo, hi));
                }
                Ok(out)
            }
            102 => {
                let n = Self::int(e, 0).max(0) as usize;
                let mut out = Vec::new();
                for k in 0..n {
                    let member = Self::pointer(e, 1 + k);
                    // (Points in a composite curve bound nothing.)
                    if self.entry(member).map(|m| m.kind == 116).unwrap_or(false) {
                        continue;
                    }
                    out.extend(self.curve(member)?);
                }
                Ok(out)
            }
            142 => match Self::pointer(e, 3) {
                0 => Err(invalid("a curve on a surface without its curve in space")),
                c => self.curve(c),
            },
            other => Err(invalid(format!("curves of type {other}"))),
        }
    }

    /// A curve entity's pieces in a surface's parameters (its x and y as u and v).
    fn curve2(&self, pointer: usize) -> Result<Vec<(Curve2, f64, f64)>> {
        let pieces = self.curve(pointer)?;
        pieces.into_iter().map(|(c, a, b)| Ok((to_2d(&c).ok_or_else(|| invalid("a curve in a surface's parameters off its plane"))?, a, b))).collect()
    }

    // --- surfaces ---

    /// A surface entity in model space, its own parameter bounds where it has them, and
    /// whether its normal runs against the file's (a surface read as an elementary one whose
    /// parametrization turns it round).
    fn surface(&self, pointer: usize) -> Result<(Surface, Option<(f64, f64, f64, f64)>, bool)> {
        let e = self.entry(pointer)?;
        let t = self.transform_of(e)?;
        let r = |k: usize| Self::real(e, k);
        let p3 = |k: usize| Point::new(r(k), r(k + 1), r(k + 2));
        let mut flipped = false;
        let (surface, bounds) = match e.kind {
            128 => {
                let (k1, k2, m1, m2) = (Self::int(e, 0).max(0) as usize, Self::int(e, 1).max(0) as usize, Self::int(e, 2).max(0) as usize, Self::int(e, 3).max(0) as usize);
                let s0 = 9;
                let su: Vec<f64> = (0..k1 + m1 + 2).map(|i| r(s0 + i)).collect();
                let t0 = s0 + k1 + m1 + 2;
                let sv: Vec<f64> = (0..k2 + m2 + 2).map(|i| r(t0 + i)).collect();
                let w0 = t0 + k2 + m2 + 2;
                let count = (k1 + 1) * (k2 + 1);
                let p0 = w0 + count;
                let mut poles = vec![vec![Point::ORIGIN; k2 + 1]; k1 + 1];
                let mut weights = vec![vec![1.0; k2 + 1]; k1 + 1];
                for j in 0..=k2 {
                    for i in 0..=k1 {
                        let n = j * (k1 + 1) + i;
                        weights[i][j] = r(w0 + n);
                        poles[i][j] = p3(p0 + 3 * n);
                    }
                }
                let b0 = p0 + 3 * count;
                let bounds = (r(b0), r(b0 + 1), r(b0 + 2), r(b0 + 3));
                // (Knots repeated past the degree inside: the poles there merged, where they
                // coincide, along either direction.)
                let (su, sv) = surface_knots(&mut poles, &mut weights, su, sv, m1, m2)?;
                let (uk, um) = distinct(&su);
                let (vk, vm) = distinct(&sv);
                let mut surface = BSplineSurface::new(poles, KnotVector::new(uk, um, m1), KnotVector::new(vk, vm, m2));
                let w = weights[0][0];
                if weights.iter().flatten().any(|&x| (x - w).abs() > 1e-15) {
                    surface = surface.with_weights(weights);
                }
                (Surface::BSpline(surface), Some(bounds))
            }
            108 => {
                let plane = Plane::from_equation(r(0), r(1), r(2), -r(3)).ok_or_else(|| invalid("a plane without a normal"))?;
                (Surface::Plane(plane), None)
            }
            190 => {
                let at = self.point(Self::pointer(e, 0))?;
                let normal = self.direction(Self::pointer(e, 1))?;
                let frame = match Self::pointer(e, 2) {
                    0 => Frame::from_z(at, normal),
                    reference => Frame::new(at, normal, self.direction(reference)?),
                };
                (Surface::Plane(Plane::from_frame(frame)), None)
            }
            120 => {
                let axis = self.curve(Self::pointer(e, 0))?;
                let [(Curve::Line(line), _, _)] = &axis[..] else { return Err(invalid("a surface of revolution whose axis is not a line")) };
                let generatrix = self.single_curve(Self::pointer(e, 1))?;
                let (sa, ta) = (r(2), r(3));
                let axis = Axis::new(line.origin, line.direction);
                // (A line revolved, as OCCT reads it: a cylinder or a cone.)
                match elementary_revolution(&generatrix, &axis) {
                    Some((surface, (v0, v1), turned)) => {
                        flipped = turned;
                        (surface, Some((sa, ta, v0, v1)))
                    }
                    None => {
                        let (_, c0, c1) = generatrix;
                        let surface = RevolutionSurface::new(generatrix.0, axis);
                        // (Its parameters: the angle about the axis, then the generatrix's.)
                        (Surface::Revolution(surface), Some((sa, ta, c0, c1)))
                    }
                }
            }
            122 => {
                let directrix = self.single_curve(Self::pointer(e, 0))?;
                let start = directrix.0.point(directrix.1);
                let end = p3(1);
                let d = end - start;
                let direction = Direction::try_from_vector(d).ok_or_else(|| invalid("a tabulated cylinder without a generatrix"))?;
                let (c0, c1) = (directrix.1, directrix.2);
                (Surface::Extrusion(ExtrusionSurface::new(directrix.0, direction)), Some((c0, c1, 0.0, d.length())))
            }
            140 => {
                // Offset by D along the side of the surface's normal its indicator points to.
                let indicator = Vector::new(r(0), r(1), r(2));
                let (basis, bounds, turned) = self.surface(Self::pointer(e, 4))?;
                flipped = turned;
                let (u0, u1, v0, v1) = bounds.unwrap_or_else(|| basis.bounds());
                let (um, vm) = (if (u1 - u0).is_finite() { 0.5 * (u0 + u1) } else { 0.0 }, if (v1 - v0).is_finite() { 0.5 * (v0 + v1) } else { 0.0 });
                let (_, du, dv) = basis.d1(um, vm);
                let sign = if du.cross(dv).dot(indicator) < 0.0 { -1.0 } else { 1.0 };
                (Surface::Offset(crate::geometry::surface::OffsetSurface::new(basis, sign * r(3))), bounds)
            }
            other => return Err(invalid(format!("surfaces of type {other}"))),
        };
        Ok((if t == Transform::IDENTITY { surface } else { surface.transformed(&t) }, bounds, flipped))
    }

    /// A curve entity of one piece (a composite of one, or a curve joined into one B-spline).
    fn single_curve(&self, pointer: usize) -> Result<Piece> {
        let mut pieces = self.curve(pointer)?;
        match pieces.len() {
            1 => Ok(pieces.remove(0)),
            0 => Err(invalid("a curve of no pieces")),
            _ => {
                // Joined: interpolated through its pieces' points.
                let mut points = Vec::new();
                for (c, a, b) in &pieces {
                    for k in 0..16 {
                        points.push(c.point(a + (b - a) * k as f64 / 16.0));
                    }
                }
                let (c, a, b) = pieces.last().expect("a piece");
                points.push(c.point(*b));
                let _ = a;
                let fitted = crate::geometry::fit::interpolate(&points, false, 0.0).map_err(|_| invalid("a composite generatrix"))?;
                let (f, l) = (fitted.first_parameter(), fitted.last_parameter());
                Ok((Curve::BSpline(fitted), f, l))
            }
        }
    }

    fn point(&self, pointer: usize) -> Result<Point> {
        let e = self.entry(pointer)?;
        let t = self.transform_of(e)?;
        Ok(Point::new(Self::real(e, 0), Self::real(e, 1), Self::real(e, 2)).transformed(&t))
    }

    fn direction(&self, pointer: usize) -> Result<Direction> {
        let e = self.entry(pointer)?;
        let t = self.transform_of(e)?;
        Direction::try_from_vector(Vector::new(Self::real(e, 0), Self::real(e, 1), Self::real(e, 2)).transformed(&t)).ok_or_else(|| invalid("a direction of no length"))
    }

    // --- shapes ---

    /// An independent entity as a shape (`None` for entities that make none).
    fn transfer(&mut self, index: usize) -> Result<Option<Shape>> {
        let e = self.file.entries[index].clone();
        let pointer = 2 * index + 1;
        Ok(Some(match e.kind {
            116 => {
                let p = self.point(pointer)?;
                let v = self.model.add_vertex(p, crate::tolerance::LINEAR);
                Shape::from(v)
            }
            100 | 102 | 104 | 106 | 110 | 112 | 126 => {
                let pieces = self.curve(pointer)?;
                let edges = self.chain(&pieces)?;
                if edges.len() == 1 {
                    Shape::from(edges[0].0)
                } else {
                    let uses: Vec<Shape> = edges.iter().map(|&(e, _)| Shape::from(e)).collect();
                    crate::make::wire(self.model, &uses)?
                }
            }
            108 | 120 | 122 | 128 | 140 | 190 => {
                let (surface, bounds, flipped) = self.surface(pointer)?;
                let (u0, u1, v0, v1) = bounds.unwrap_or_else(|| surface.bounds());
                let face = crate::make::face_on_domain(self.model, surface, u0, u1, v0, v1)?;
                if flipped { face.reversed() } else { face }
            }
            144 => self.trimmed(&e)?,
            143 => self.bounded(&e)?,
            510 => self.brep_face(&e)?,
            514 => self.brep_shell(&e)?,
            186 => self.brep_solid(&e)?,
            402 => {
                let n = Self::int(&e, 0).max(0) as usize;
                let c = self.model.new_empty(NodeKind::Compound);
                for k in 0..n {
                    let member = Self::pointer(&e, 1 + k);
                    if let Ok(Some(s)) = self.transfer((member.max(1) - 1) / 2) {
                        self.model.add(c, s)?;
                    }
                }
                c
            }
            408 => {
                let definition = self.entry(Self::pointer(&e, 0))?.clone();
                let n = Self::int(&definition, 2).max(0) as usize;
                let c = self.model.new_empty(NodeKind::Compound);
                for k in 0..n {
                    let member = Self::pointer(&definition, 3 + k);
                    if let Ok(Some(s)) = self.transfer((member.max(1) - 1) / 2) {
                        self.model.add(c, s)?;
                    }
                }
                let scale = match e.params.get(4).and_then(as_number) {
                    Some(s) if s != 0.0 => s,
                    _ => 1.0,
                };
                let t = self.transform_of(&e)? * Transform::translation(Vector::new(Self::real(&e, 1), Self::real(&e, 2), Self::real(&e, 3))) * Transform::scaling(Point::ORIGIN, scale);
                crate::modify::transform(self.model, c, &t, Default::default())?
            }
            _ => return Ok(None),
        }))
    }

    /// Edges along pieces of curve, each starting at the vertex the one before ends at
    /// (where they meet), the last ending at the first's start where the chain closes.
    fn chain(&mut self, pieces: &[Piece]) -> Result<Vec<(EdgeId, Piece)>> {
        let mut out: Vec<(EdgeId, Piece)> = Vec::new();
        // (Ends meet within the file's resolution, or a ten-thousandth of the pieces' size:
        // files are often only so accurate. The vertex's tolerance then covers the gap.)
        let size = pieces.iter().map(|(c, a, b)| c.point(*a).distance(c.point(0.5 * (a + b))) + c.point(0.5 * (a + b)).distance(c.point(*b))).fold(0.0, f64::max);
        let reach = (10.0 * self.file.resolution).max(1e-4 * size);
        let close = |a: Point, b: Point| a.distance(b) <= reach.max(1e-6 * (1.0 + a.distance(Point::ORIGIN)));
        let covering = |model: &mut Model, v: VertexId, p: Point| -> Result<()> {
            let gap = model.point(Shape::from(v)).distance(p);
            if gap > model.tolerance(Shape::from(v)) {
                let at = model.point(Shape::from(v));
                model.set_vertex_point(v.into(), at, 1.1 * gap)?;
            }
            Ok(())
        };
        let mut first_vertex: Option<(VertexId, Point)> = None;
        let mut previous: Option<(VertexId, Point)> = None;
        for (k, (curve, a, b)) in pieces.iter().enumerate() {
            let (p, q) = (curve.point(*a), curve.point(*b));
            let start = match previous {
                Some((v, at)) if close(at, p) => {
                    covering(self.model, v, p)?;
                    v
                }
                _ => self.model.add_vertex(p, crate::tolerance::LINEAR),
            };
            if first_vertex.is_none() {
                first_vertex = Some((start, p));
            }
            let last = k + 1 == pieces.len();
            let end = match first_vertex {
                Some((v, at)) if last && close(at, q) => {
                    covering(self.model, v, q)?;
                    v
                }
                _ if close(p, q) && pieces.len() == 1 => start,
                _ => self.model.add_vertex(q, crate::tolerance::LINEAR),
            };
            let tol = p.distance(curve.point(*a)).max(crate::tolerance::LINEAR);
            let edge = self.model.add_edge(start, end, tol)?;
            self.model.set_edge_curve(edge.into(), curve.clone(), crate::tolerance::LINEAR)?;
            self.model.set_edge_range(edge.into(), *a, *b)?;
            out.push((edge, (curve.clone(), *a, *b)));
            previous = Some((end, q));
        }
        Ok(out)
    }

    /// The curve on a surface of an edge: the file's (in the surface's parameters, over the
    /// edge's range) where it fits, else the edge's curve projected, else its points' nearest
    /// on the surface interpolated (files are often only so accurate: the tolerances then
    /// cover the gaps).
    fn pcurve(&self, piece: &Piece, given: Option<&(Curve2, f64, f64)>, surface: &Surface) -> Result<Curve2> {
        let span = CurveSpan { curve: &piece.0, first: piece.1, last: piece.2 };
        // (The file's parameters are the model's only for B-spline surfaces.)
        let given = given.filter(|_| matches!(surface, Surface::BSpline(_))).map(|g| over_range(g, piece.1, piece.2));
        if let Some(pc) = &given {
            if fits(pc, &span, surface) || follows(pc, &span, surface) {
                return Ok(pc.clone());
            }
        }
        if let Ok(pc) = project(&span, surface, given.is_some()) {
            return Ok(pc);
        }
        if let Some(pc) = given {
            return Ok(pc);
        }
        nearest_pcurve(piece, surface).ok_or_else(|| invalid("an edge off its face"))
    }

    /// A loop of curves (in space, and in the surface's parameters where given) as coedges
    /// with their curves on the surface.
    fn surface_loop(&mut self, pieces: &[Piece], params: Option<&[(Curve2, f64, f64)]>, surface: &Surface) -> Result<Vec<(EdgeId, bool, Curve2, f64, f64)>> {
        let edges = self.chain(pieces)?;
        let params = params.filter(|p| p.len() == edges.len());
        let mut coedges = Vec::with_capacity(edges.len());
        for (k, (edge, piece)) in edges.iter().enumerate() {
            let pc = self.pcurve(piece, params.map(|p| &p[k]), surface)?;
            coedges.push((*edge, true, pc, piece.1, piece.2));
        }
        align_periods(&mut coedges, surface);
        Ok(coedges)
    }

    /// A curve on a parametric surface (142): its pieces in space (computed from those in
    /// the surface's parameters where it gives no curve in space) and in the parameters.
    fn curve_on_surface(&self, pointer: usize, surface: &Surface) -> Result<(Vec<Piece>, Option<Vec<(Curve2, f64, f64)>>)> {
        let e = self.entry(pointer)?;
        if e.kind != 142 {
            return Ok((self.curve(pointer)?, None));
        }
        let (b, c) = (Self::pointer(e, 2), Self::pointer(e, 3));
        let params = if b != 0 { self.curve2(b).ok() } else { None };
        let space = if c != 0 {
            self.curve(c)?
        } else {
            // In space through the surface: each piece of the curve in its parameters mapped.
            let params = params.as_ref().ok_or_else(|| invalid("a curve on a surface with neither of its curves"))?;
            params
                .iter()
                .map(|(pc, a, b)| {
                    let points: Vec<Point> = (0..=32).map(|k| { let uv = pc.point(a + (b - a) * k as f64 / 32.0); surface.point(uv.x, uv.y) }).collect();
                    let params: Vec<f64> = (0..=32).map(|k| a + (b - a) * k as f64 / 32.0).collect();
                    let fitted = crate::geometry::fit::interpolate_with_parameters(&points, &params).map_err(|_| invalid("a curve on a surface"))?;
                    Ok((Curve::BSpline(fitted), *a, *b))
                })
                .collect::<Result<Vec<_>>>()?
        };
        Ok((space, params))
    }

    /// A trimmed surface (144) as a face.
    fn trimmed(&mut self, e: &Entry) -> Result<Shape> {
        let (surface, bounds, flipped) = self.surface(Self::pointer(e, 0))?;
        let (n1, n2) = (Self::int(e, 1), Self::int(e, 2).max(0) as usize);
        let inner: Vec<usize> = (0..n2).map(|k| Self::pointer(e, 4 + k)).collect();
        let face = if n1 == 0 {
            // Bounded by the surface's own boundary; the holes cut out of it.
            let (u0, u1, v0, v1) = bounds.unwrap_or_else(|| surface.bounds());
            let face = crate::make::face_on_domain(self.model, surface.clone(), u0, u1, v0, v1)?;
            if inner.is_empty() {
                face
            } else {
                let mut holes = Vec::new();
                for p in inner {
                    let (space, _) = self.curve_on_surface(p, &surface)?;
                    let edges = self.chain(&space)?;
                    let uses: Vec<Shape> = edges.iter().map(|&(e, _)| Shape::from(e)).collect();
                    holes.push(crate::make::wire(self.model, &uses)?);
                }
                crate::make::face_with_holes(self.model, face, &holes)?
            }
        } else {
            let mut loops = Vec::new();
            for (k, p) in std::iter::once(Self::pointer(e, 3)).chain(inner).enumerate() {
                let (space, params) = self.curve_on_surface(p, &surface)?;
                let coedges = self.surface_loop(&space, params.as_deref(), &surface)?;
                loops.push((k == 0, coedges));
            }
            let face = self.model.new_empty(NodeKind::Face);
            self.model.set_face_surface(face, surface.clone(), crate::tolerance::LINEAR)?;
            finish_face(self.model, &mut self.gaps, face, &surface, loops, true)?
        };
        let face = if flipped { face.reversed() } else { face };
        self.placed(e, face)
    }

    /// A bounded surface (143) as a face: its boundaries (141), the first the outer one.
    fn bounded(&mut self, e: &Entry) -> Result<Shape> {
        let (surface, _, flipped) = self.surface(Self::pointer(e, 1))?;
        let n = Self::int(e, 2).max(0) as usize;
        let mut loops = Vec::new();
        for k in 0..n {
            let b = self.entry(Self::pointer(e, 3 + k))?.clone();
            let curves = Self::int(&b, 3).max(0) as usize;
            let (mut space, mut params): (Vec<Piece>, Vec<(Curve2, f64, f64)>) = (Vec::new(), Vec::new());
            let mut at = 4;
            let mut with_params = true;
            for _ in 0..curves {
                let (curve, sense, count) = (Self::pointer(&b, at), Self::int(&b, at + 1), Self::int(&b, at + 2).max(0) as usize);
                let mut pieces = self.curve(curve)?;
                let mut pcs = if count >= 1 { self.curve2(Self::pointer(&b, at + 3)).unwrap_or_default() } else { Vec::new() };
                if pcs.len() != pieces.len() {
                    with_params = false;
                }
                if sense == 2 {
                    pieces.reverse();
                    pieces = pieces.into_iter().map(reversed_piece).collect();
                    pcs.reverse();
                    pcs = pcs.into_iter().map(reversed_piece2).collect();
                }
                space.extend(pieces);
                params.extend(pcs);
                at += 3 + count;
            }
            let coedges = self.surface_loop(&space, with_params.then_some(params.as_slice()), &surface)?;
            loops.push((k == 0, coedges));
        }
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface.clone(), crate::tolerance::LINEAR)?;
        let face = finish_face(self.model, &mut self.gaps, face, &surface, loops, !flipped)?;
        self.placed(e, face)
    }

    /// `shape` placed by the entity's transformation matrix, where it has one.
    fn placed(&mut self, e: &Entry, shape: Shape) -> Result<Shape> {
        let t = self.transform_of(e)?;
        if t == Transform::IDENTITY {
            return Ok(shape);
        }
        crate::modify::transform(self.model, shape, &t, Default::default())
    }

    // --- boundary representations ---

    fn brep_vertex(&mut self, list: usize, index: usize) -> Result<VertexId> {
        if let Some(&v) = self.vertices_made.get(&(list, index)) {
            return Ok(v);
        }
        let e = self.entry(list)?.clone();
        let t = self.transform_of(&e)?;
        let k = 1 + 3 * (index.max(1) - 1);
        let p = Point::new(Self::real(&e, k), Self::real(&e, k + 1), Self::real(&e, k + 2)).transformed(&t);
        let v = self.model.add_vertex(p, crate::tolerance::LINEAR);
        self.vertices_made.insert((list, index), v);
        Ok(v)
    }

    fn brep_edge(&mut self, list: usize, index: usize) -> Result<(EdgeId, Piece)> {
        if let Some(found) = self.edges_made.get(&(list, index)) {
            return Ok(found.clone());
        }
        let e = self.entry(list)?.clone();
        let k = 1 + 5 * (index.max(1) - 1);
        let curve = Self::pointer(&e, k);
        let (sv, si, tv, ti) = (Self::pointer(&e, k + 1), Self::pointer(&e, k + 2), Self::pointer(&e, k + 3), Self::pointer(&e, k + 4));
        let piece = self.single_curve(curve)?;
        let (a, b) = (self.brep_vertex(sv, si)?, self.brep_vertex(tv, ti)?);
        let edge = self.model.add_edge(a, b, crate::tolerance::LINEAR)?;
        self.model.set_edge_curve(edge.into(), piece.0.clone(), crate::tolerance::LINEAR)?;
        self.model.set_edge_range(edge.into(), piece.1, piece.2)?;
        self.edges_made.insert((list, index), (edge, piece.clone()));
        Ok((edge, piece))
    }

    fn brep_face(&mut self, e: &Entry) -> Result<Shape> {
        let (surface, _, flipped) = self.surface(Self::pointer(e, 0))?;
        let n = Self::int(e, 1).max(0) as usize;
        let outer_first = Self::int(e, 2) == 1;
        let mut loops = Vec::new();
        for k in 0..n {
            let l = self.entry(Self::pointer(e, 3 + k))?.clone();
            let count = Self::int(&l, 0).max(0) as usize;
            let mut coedges = Vec::new();
            let mut at = 1;
            for _ in 0..count {
                let (kind, list, index, forward, curves) = (Self::int(&l, at), Self::pointer(&l, at + 1), Self::pointer(&l, at + 2), Self::int(&l, at + 3) == 1, Self::int(&l, at + 4).max(0) as usize);
                let given = if curves >= 1 { self.curve2(Self::pointer(&l, at + 6)).ok().and_then(|mut p| (p.len() == 1).then(|| p.remove(0))) } else { None };
                at += 5 + 2 * curves;
                if kind != 0 {
                    continue; // a vertex loop
                }
                let (edge, piece) = self.brep_edge(list, index)?;
                // (A curve given along the use, not the edge, runs the other way.)
                let given = given.map(|g| if forward { g } else { reversed_piece2(g) });
                let pc = self.pcurve(&piece, given.as_ref(), &surface)?;
                coedges.push((edge, forward, pc, piece.1, piece.2));
            }
            align_periods(&mut coedges, &surface);
            loops.push((outer_first && k == 0, coedges));
        }
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface.clone(), crate::tolerance::LINEAR)?;
        finish_face(self.model, &mut self.gaps, face, &surface, loops, !flipped)
    }

    fn brep_shell(&mut self, e: &Entry) -> Result<Shape> {
        let n = Self::int(e, 0).max(0) as usize;
        let shell = self.model.new_empty(NodeKind::Shell);
        for k in 0..n {
            let f = self.entry(Self::pointer(e, 1 + 2 * k))?.clone();
            let agrees = Self::int(e, 2 + 2 * k) != 0;
            let face = self.brep_face(&f)?;
            self.model.add(shell, if agrees { face } else { face.reversed() })?;
        }
        Ok(shell)
    }

    fn brep_solid(&mut self, e: &Entry) -> Result<Shape> {
        let solid = self.model.new_empty(NodeKind::Solid);
        let outer = self.entry(Self::pointer(e, 0))?.clone();
        let agrees = Self::int(e, 1) != 0;
        let shell = self.brep_shell(&outer)?;
        self.model.add(solid, if agrees { shell } else { shell.reversed() })?;
        let voids = Self::int(e, 2).max(0) as usize;
        for k in 0..voids {
            let v = self.entry(Self::pointer(e, 3 + 2 * k))?.clone();
            let agrees = Self::int(e, 4 + 2 * k) != 0;
            let shell = self.brep_shell(&v)?;
            self.model.add(solid, if agrees { shell } else { shell.reversed() })?;
        }
        let solid = self.placed(e, solid)?;
        // (Outward: a positive volume.)
        if crate::properties::volume(self.model, solid).mass < 0.0 {
            return Ok(solid.reversed());
        }
        Ok(solid)
    }
}

/// The cylinder or cone a line revolved about `axis` sweeps (its angle starting at the line,
/// its frame at the line's end farther from the axis), with the range of its second parameter
/// the line covers; `None` for other generatrices, and for lines square to the axis or askew
/// to it.
fn elementary_revolution(piece: &Piece, axis: &Axis) -> Option<(Surface, (f64, f64), bool)> {
    let (surface, range) = elementary_revolution_of(piece, axis)?;
    // (Whether its normal runs against the surface of revolution's, at the line's start: as
    // OCCT makes that, about the axis run from its end to its start, by the angle about it
    // and then the generatrix's parameter, its normal the angle's derivative across the
    // generatrix's.)
    let (curve, a, _) = piece;
    let (p, tangent) = curve.d1(*a);
    let d = -Vector::from(axis.direction);
    let theirs = d.cross(p - axis.origin).cross(tangent);
    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
    let (_, du, dv) = surface.d1(u, v);
    let turned = du.cross(dv).dot(theirs) < 0.0;
    Some((surface, range, turned))
}

fn elementary_revolution_of((curve, a, b): &Piece, axis: &Axis) -> Option<(Surface, (f64, f64))> {
    use crate::geometry::surface::{Cone, Cylinder};
    let Curve::Line(_) = curve else { return None };
    // (From the end farther from the axis: a cone's radius there leaves room for offsets.)
    let (p, q) = (curve.point(*a), curve.point(*b));
    let distance = |x: Point| (x - axis.origin).cross(Vector::from(axis.direction)).length();
    let (p, q) = if distance(q) > distance(p) { (q, p) } else { (p, q) };
    let d = Vector::from(axis.direction);
    let foot = |x: Point| axis.origin + d * (x - axis.origin).dot(d);
    let radial = |x: Point| x - foot(x);
    let (rp, rq) = (radial(p), radial(q));
    let (r0, r1) = (rp.length(), rq.length());
    let length = p.distance(q);
    let size = 1.0 + length + r0 + r1;
    // (In a plane through the axis.)
    if rp.cross(rq).dot(d).abs() > 1e-9 * size * size {
        return None;
    }
    let along = (q - p).dot(d);
    let x = Direction::try_from_vector(if r0 > 1e-9 * size { rp } else { rq })?;
    let frame = Frame::new(foot(p), axis.direction, x);
    if (r1 - r0).abs() <= 1e-9 * size {
        // Parallel to the axis: a cylinder, its height from the line's start.
        if r0 <= 1e-9 * size {
            return None;
        }
        return Some((Surface::Cylinder(Cylinder::new(frame, r0)), (along.min(0.0), along.max(0.0))));
    }
    if along.abs() <= 1e-9 * size {
        return None; // square to the axis: a plane
    }
    // A cone: its generator's length from the line's start (the radius there r0), the sign of
    // its axial component.
    let v1 = if along > 0.0 { length } else { -length };
    let semi = ((r1 - r0) / v1).atan2(along / v1);
    if semi == 0.0 || semi.abs() >= std::f64::consts::FRAC_PI_2 {
        return None;
    }
    Some((Surface::Cone(Cone::new(frame, semi, r0)), (v1.min(0.0), v1.max(0.0))))
}

/// A curve in a surface's parameters run over `[first, last]` instead of its own range
/// (linearly): sampled there and interpolated, unless its range is that already.
fn over_range((curve, a, b): &(Curve2, f64, f64), first: f64, last: f64) -> Curve2 {
    if (a - first).abs() <= 1e-12 * (1.0 + first.abs()) && (b - last).abs() <= 1e-12 * (1.0 + last.abs()) {
        return curve.clone();
    }
    let n = 32;
    let ts: Vec<f64> = (0..=n).map(|k| first + (last - first) * k as f64 / n as f64).collect();
    let uv: Vec<Point2> = (0..=n).map(|k| curve.point(a + (b - a) * k as f64 / n as f64)).collect();
    match crate::geometry::fit::interpolate_with_parameters(&uv, &ts) {
        Ok(c) => Curve2::BSpline(c),
        Err(_) => curve.clone(),
    }
}

/// The curve on `surface` through the nearest points of it to the piece's (each found from
/// the last), interpolated over the piece's range.
fn nearest_pcurve((curve, a, b): &Piece, surface: &Surface) -> Option<Curve2> {
    let n = 64;
    let ts: Vec<f64> = (0..=n).map(|k| a + (b - a) * k as f64 / n as f64).collect();
    let mut from: Option<(f64, f64)> = None;
    let uv: Vec<Point2> = ts
        .iter()
        .map(|&t| {
            let (u, v, _) = crate::geometry::extrema::nearest_on_surface(curve.point(t), surface, from);
            from = Some((u, v));
            Point2::new(u, v)
        })
        .collect();
    crate::geometry::fit::interpolate_with_parameters(&uv, &ts).ok().map(Curve2::BSpline)
}

/// Whether flat knots of `degree` suit `poles` poles: as many as the poles and the degree
/// need, non-decreasing, none repeated more than the degree inside or once more at the ends.
fn sound_knots(flat: &[f64], poles: usize, degree: usize) -> bool {
    if flat.len() != poles + degree + 1 || flat.windows(2).any(|w| w[1] < w[0]) || poles <= degree {
        return false;
    }
    let (_, mults) = distinct(flat);
    let n = mults.len();
    n >= 2 && mults[0] <= degree + 1 && mults[n - 1] <= degree + 1 && mults[1..n - 1].iter().all(|&m| m <= degree)
}

/// The interior run of a knot repeated more than `degree` times: its first flat index and
/// length.
fn full_run(flat: &[f64], degree: usize) -> Option<(usize, usize)> {
    let mut i = degree + 1;
    while i < flat.len().saturating_sub(degree + 1) {
        let mut j = i;
        while j + 1 < flat.len() && flat[j + 1] == flat[i] {
            j += 1;
        }
        let m = j - i + 1;
        if m > degree && j + 1 < flat.len() - degree {
            return Some((i, m));
        }
        i = j + 1;
    }
    None
}

/// A curve's poles, weights and flat knots in parts with knots of at most `degree`
/// multiplicity inside: at a knot repeated once more, the curve's two sides' poles merged
/// where they coincide, else the curve split there.
fn curve_parts(poles: Vec<Point>, weights: Vec<f64>, flat: Vec<f64>, degree: usize) -> Result<Vec<(Vec<Point>, Vec<f64>, Vec<f64>)>> {
    let (mut poles, mut weights, mut flat) = (poles, weights, flat);
    let mut out = Vec::new();
    loop {
        let Some((s, m)) = full_run(&flat, degree) else { break };
        if m > degree + 1 {
            return Err(invalid("a B-spline curve with a knot repeated more than once past its degree"));
        }
        let (left, right) = (s - 1, s);
        let close = poles[left].distance(poles[right]) <= 1e-9 * (1.0 + poles[left].distance(Point::ORIGIN)) && (weights[left] - weights[right]).abs() <= 1e-12;
        if close {
            poles.remove(right);
            weights.remove(right);
            flat.remove(s);
        } else {
            out.push((poles[..s].to_vec(), weights[..s].to_vec(), flat[..s + degree + 1].to_vec()));
            poles = poles[s..].to_vec();
            weights = weights[s..].to_vec();
            flat = flat[s..].to_vec();
        }
    }
    out.push((poles, weights, flat));
    for (poles, _, flat) in &out {
        if !sound_knots(flat, poles.len(), degree) {
            return Err(invalid("a B-spline curve whose knots do not suit its poles"));
        }
    }
    Ok(out)
}

/// A surface's flat knots along u and v with the poles (rows along u, columns along v)
/// merged where a knot repeats past the degree and the rows (columns) on both sides of it
/// coincide; an error where they do not.
fn surface_knots(poles: &mut Vec<Vec<Point>>, weights: &mut Vec<Vec<f64>>, su: Vec<f64>, sv: Vec<f64>, du: usize, dv: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let (mut su, mut sv) = (su, sv);
    while let Some((s, m)) = full_run(&su, du) {
        let coincide = m == du + 1 && (0..poles[0].len()).all(|j| poles[s - 1][j].distance(poles[s][j]) <= 1e-9 * (1.0 + poles[s][j].distance(Point::ORIGIN)));
        if !coincide {
            return Err(invalid("a B-spline surface broken across a knot"));
        }
        poles.remove(s);
        weights.remove(s);
        su.remove(s);
    }
    while let Some((s, m)) = full_run(&sv, dv) {
        let coincide = m == dv + 1 && (0..poles.len()).all(|i| poles[i][s - 1].distance(poles[i][s]) <= 1e-9 * (1.0 + poles[i][s].distance(Point::ORIGIN)));
        if !coincide {
            return Err(invalid("a B-spline surface broken across a knot"));
        }
        for row in poles.iter_mut() {
            row.remove(s);
        }
        for row in weights.iter_mut() {
            row.remove(s);
        }
        sv.remove(s);
    }
    if !sound_knots(&su, poles.len(), du) || !sound_knots(&sv, poles[0].len(), dv) {
        return Err(invalid("a B-spline surface whose knots do not suit its poles"));
    }
    Ok((su, sv))
}

/// Flat knots as distinct knots and their multiplicities.
fn distinct(flat: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let (mut knots, mut mults): (Vec<f64>, Vec<usize>) = (Vec::new(), Vec::new());
    for &k in flat {
        match knots.last() {
            Some(&last) if (k - last).abs() <= 1e-12 * (1.0 + k.abs()) => *mults.last_mut().expect("a multiplicity") += 1,
            _ => {
                knots.push(k);
                mults.push(1);
            }
        }
    }
    (knots, mults)
}

/// A curve in a plane z = constant as one in (x, y): lines, circles, ellipses and B-splines.
fn to_2d(curve: &Curve) -> Option<Curve2> {
    let flat = |p: Point| Point2::new(p.x, p.y);
    Some(match curve {
        Curve::Line(l) => Curve2::Line(Line2::new(flat(l.origin), Direction2::try_new(l.direction.x(), l.direction.y())?)),
        Curve::Circle(c) => {
            let f = c.frame();
            let ccw = f.z().z() > 0.0;
            Curve2::Circle(Circle2::new(Frame2::new(flat(f.origin()), Direction2::try_new(f.x().x(), f.x().y())?, ccw), c.radius()))
        }
        Curve::Ellipse(c) => {
            let f = c.frame();
            let ccw = f.z().z() > 0.0;
            Curve2::Ellipse(Ellipse2::new(Frame2::new(flat(f.origin()), Direction2::try_new(f.x().x(), f.x().y())?, ccw), c.major_radius(), c.minor_radius()))
        }
        Curve::BSpline(b) => {
            let poles: Vec<Point2> = b.poles().iter().map(|&p| flat(p)).collect();
            let mut c = BSplineCurve2::new(poles, b.knots().to_vec(), b.multiplicities().to_vec(), b.degree());
            if let Some(w) = b.weights() {
                c = c.with_weights(w.to_vec());
            }
            Curve2::BSpline(c)
        }
        _ => return None,
    })
}

/// A piece run the other way.
fn reversed_piece((curve, a, b): Piece) -> Piece {
    let mut c = curve;
    let (ra, rb) = (c.reversed_parameter(b), c.reversed_parameter(a));
    c.reverse();
    (c, ra, rb)
}

fn reversed_piece2((curve, a, b): (Curve2, f64, f64)) -> (Curve2, f64, f64) {
    let mut c = curve;
    let (ra, rb) = (c.reversed_parameter(b), c.reversed_parameter(a));
    c.reverse();
    (c, ra, rb)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Lines of an IGES file: data in the first 72 columns, the section and a sequence number.
    fn record(data: &str, section: char, seq: usize) -> String {
        format!("{data:<72}{section}{seq:>7}\n")
    }

    fn directory(kind: i32, param: usize, count: usize, status: &str, seq: usize) -> String {
        let a = format!("{kind:>8}{param:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{status:>8}", 0, 0, 0, 0, 0, 0);
        let b = format!("{kind:>8}{:>8}{:>8}{count:>8}{:>8}", 0, 0, 0);
        record(&a, 'D', seq) + &record(&b, 'D', seq + 1)
    }

    /// A unit square trimmed out of a B-spline plane (degree 1, 2 × 2 poles), in inches: a
    /// face of the square's area, in millimetres.
    #[test]
    fn a_trimmed_bilinear_patch_in_inches() {
        let mut text = record("", 'S', 1);
        text += &record("1H,,1H;,,,,,,,,,,,1.,1,2HIN;", 'G', 1);
        // 1: the surface (128); 3: the outer boundary (142) on it, of 5: a composite of 7..13
        // lines (110); 15: the trimmed surface (144).
        text += &directory(128, 1, 2, "00010000", 1);
        text += &directory(142, 3, 1, "00010000", 3);
        text += &directory(102, 4, 1, "00010000", 5);
        let mut seq = 7;
        for k in 0..4 {
            text += &directory(110, 5 + k, 1, "00010000", seq);
            seq += 2;
        }
        text += &directory(144, 9, 1, "00000000", 15);
        let mut p = Vec::new();
        p.push(record("128,1,1,1,1,0,0,1,0,0,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,", 'P', 1));
        p.push(record("0.,0.,0.,2.,0.,0.,0.,2.,0.,2.,2.,0.,0.,1.,0.,1.;", 'P', 2));
        p.push(record("142,1,1,0,5,2;", 'P', 3));
        p.push(record("102,4,7,9,11,13;", 'P', 4));
        let corners = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)];
        for k in 0..4 {
            let (a, b) = (corners[k], corners[(k + 1) % 4]);
            p.push(record(&format!("110,{},{},0.,{},{},0.;", a.0, a.1, b.0, b.1), 'P', 5 + k));
        }
        p.push(record("144,1,1,0,3;", 'P', 9));
        // (The parameter lines name their directory entries in columns 66 to 72: not read.)
        for line in p {
            text += &line;
        }
        let mut model = Model::new();
        let shape = read(&mut model, &mut text.as_bytes()).unwrap();
        let faces = model.unique_subshapes(shape, NodeKind::Face);
        assert_eq!(faces.len(), 1);
        let area = crate::properties::area(&model, shape).mass;
        assert!((area - 25.4 * 25.4).abs() < 1e-6, "{area}");
        assert!(crate::validity::check(&model, shape).is_valid());
    }
}
