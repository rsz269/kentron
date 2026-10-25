//! Text as planar faces, from the glyph outlines of a TrueType or OpenType font.
//!
//! As OCCT's `Font_BRepTextBuilder` lays text out and builds it: each character's glyph at
//! the pen, which moves on by the glyph's advance (and the font's kerning of the pair, from
//! its `kern` table); a line feed starts a new line a line spacing lower, a tab moves on to
//! the next of every 8 spaces. Each contour of a glyph is a closed wire of an edge per
//! segment of the outline (lines, and quadratic or cubic Bézier curves); contours inside an
//! even number of others bound faces, those inside an odd number are their holes. The faces
//! lie in the XY plane, facing +Z, the em `size` high. Characters the font has no glyph for
//! are set in a fallback font for their script where one is given, else as the font's
//! missing-glyph glyph.

use crate::error::{Error, Result};
use crate::geometry::curve::BezierCurve;
use crate::geometry::curve::{Curve, Line};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Point, Vector};
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// Where a text's lines stand across: their pen's start at x = 0 (`Left`), their middle
/// (`Center`) or their end (`Right`), each by the line's advance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HorizontalAlignment {
    #[default]
    Left,
    Center,
    Right,
}

/// Where a text stands up and down: the font's ascender of the first line at y = 0
/// (`Top`), the last line's baseline there (`Bottom`), or the middle of its lines (`Center`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VerticalAlignment {
    #[default]
    Top,
    Center,
    Bottom,
}

/// The scripts text set in a font can take glyphs the font lacks from a fallback font for,
/// as OCCT's `Font_UnicodeSubset` divides characters.
#[doc(alias = "Font_UnicodeSubset")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Script {
    Western,
    /// Chinese, Japanese and Korean ideographs, with Japanese kana.
    Cjk,
    /// Korean Hangul.
    Korean,
    Arabic,
}

impl Script {
    /// The script `c` is set in.
    pub fn of(c: char) -> Script {
        let u = u32::from(c);
        let within = |ranges: &[(u32, u32)]| ranges.iter().any(|&(a, b)| (a..=b).contains(&u));
        if within(&[(0x3400, 0x4DFF), (0x4E00, 0x9FFF), (0xF900, 0xFAFF), (0x20000, 0x2A6DF), (0x2F800, 0x2FA1F), (0x3040, 0x309F), (0x30A0, 0x30FF)]) {
            Script::Cjk
        } else if within(&[(0x1100, 0x11FF), (0x3130, 0x318F), (0xAC00, 0xD7A3)]) {
            Script::Korean
        } else if within(&[(0x600, 0x6FF)]) {
            Script::Arabic
        } else {
            Script::Western
        }
    }
}

/// A font for the characters of a script that the text's own font has no glyph for.
#[derive(Debug, Clone, Copy)]
pub struct Fallback<'a> {
    pub script: Script,
    /// The bytes of a TrueType or OpenType file (or collection).
    pub font: &'a [u8],
    /// Which face of the file.
    pub index: u32,
}

/// The faces of `text` set in the font `font` (the bytes of a TrueType or OpenType file,
/// its first face) at em height `size`, aligned as given: a compound of faces in the XY
/// plane facing +Z (empty for text with no glyph outlines, as spaces).
#[doc(alias = "Font_BRepTextBuilder")]
#[doc(alias = "StdPrs_BRepFont")]
pub fn text(model: &mut Model, font: &[u8], text: &str, size: f64, horizontal: HorizontalAlignment, vertical: VerticalAlignment) -> Result<Shape> {
    text_with_fallbacks(model, font, &[], text, size, horizontal, vertical)
}

/// [`text`], with the glyphs `font` lacks taken from the fallback for their characters'
/// script where it has them (at the same em height), and the font's missing-glyph glyph
/// (`.notdef`) for the others, as OCCT sets them. The lines stand by `font`'s metrics.
pub fn text_with_fallbacks(model: &mut Model, font: &[u8], fallbacks: &[Fallback], text: &str, size: f64, horizontal: HorizontalAlignment, vertical: VerticalAlignment) -> Result<Shape> {
    if !(size > 0.0 && size.is_finite()) {
        return Err(Error::InvalidInput(format!("text size must be positive, got {size}")));
    }
    let face = ttf_parser::Face::parse(font, 0).map_err(|e| Error::InvalidInput(format!("not a font: {e}")))?;
    let mut faces = vec![(face, Script::Western)];
    for f in fallbacks {
        let face = ttf_parser::Face::parse(f.font, f.index).map_err(|e| Error::InvalidInput(format!("not a font: {e}")))?;
        faces.push((face, f.script));
    }
    let scales: Vec<f64> = faces.iter().map(|(f, _)| size / f64::from(f.units_per_em())).collect();
    // Each character's face (an index into `faces`) and glyph.
    let glyph = |c: char| -> (usize, ttf_parser::GlyphId) {
        if let Some(id) = faces[0].0.glyph_index(c) {
            return (0, id);
        }
        let script = Script::of(c);
        for (k, (f, s)) in faces.iter().enumerate().skip(1) {
            if *s == script {
                if let Some(id) = f.glyph_index(c) {
                    return (k, id);
                }
            }
        }
        (0, ttf_parser::GlyphId(0))
    };
    let face = &faces[0].0;
    let scale = scales[0];
    // (A font with CFF outlines stands as high and low as its glyphs' bounding box, with no
    // gap between lines, as OCCT sets it; one with TrueType outlines by its own metrics.)
    let tables = face.tables();
    let (ascender, descender, gap) = if tables.cff.is_some() || tables.cff2.is_some() {
        let b = face.global_bounding_box();
        (f64::from(b.y_max), f64::from(b.y_min), 0.0)
    } else {
        (f64::from(face.ascender()), f64::from(face.descender()), f64::from(face.line_gap()))
    };
    let line_spacing = (ascender - descender + gap) * scale;
    let ascender = ascender * scale;
    let advance_of = |c: char, next: Option<char>| {
        let (k, id) = glyph(c);
        // (Kerning only between glyphs of one face.)
        let next = next.map(glyph).filter(|&(n, _)| n == k).map(|(_, id)| id);
        advance(&faces[k].0, id, next) * scales[k]
    };
    // The glyphs' pen positions, line by line.
    let chars: Vec<char> = text.chars().collect();
    let mut lines: Vec<(Vec<(char, f64)>, f64)> = vec![(Vec::new(), 0.0)];
    let mut column = 0;
    for (k, &c) in chars.iter().enumerate() {
        let next = chars.get(k + 1).copied();
        let (glyphs, pen) = lines.last_mut().expect("a line");
        match c {
            '\n' => {
                lines.push((Vec::new(), 0.0));
                column = 0;
            }
            '\t' => {
                let spaces = 8 - column % 8;
                *pen += advance_of(' ', next) * spaces as f64;
                column += spaces;
            }
            c if c.is_control() => {}
            c => {
                glyphs.push((c, *pen));
                *pen += advance_of(c, next);
                column += 1;
            }
        }
    }
    let count = lines.len() as f64;
    let shift_y = match vertical {
        VerticalAlignment::Top => 0.0,
        VerticalAlignment::Bottom => ascender + (count - 1.0) * line_spacing,
        VerticalAlignment::Center => 0.5 * line_spacing * count,
    };
    let compound = model.new_empty(NodeKind::Compound);
    let tolerance = 1e-9 * size.max(1.0);
    for (row, (glyphs, width)) in lines.iter().enumerate() {
        let shift_x = match horizontal {
            HorizontalAlignment::Left => 0.0,
            HorizontalAlignment::Center => -0.5 * width,
            HorizontalAlignment::Right => -width,
        };
        let baseline = -ascender - row as f64 * line_spacing + shift_y;
        for &(c, pen) in glyphs {
            let (k, id) = glyph(c);
            let mut outline = Outline::default();
            if faces[k].0.outline_glyph(id, &mut outline).is_none() {
                continue;
            }
            let origin = Vector::new(pen + shift_x, baseline, 0.0);
            let scale = scales[k];
            let at = |(x, y): (f32, f32)| Point::new(f64::from(x) * scale, f64::from(y) * scale, 0.0) + origin;
            let contours: Vec<Vec<Segment<Point>>> = outline.contours.iter().map(|c| c.iter().map(|s| s.map(at)).collect()).collect();
            for f in glyph_faces(model, &contours, tolerance)? {
                model.add(compound, f)?;
            }
        }
    }
    Ok(compound)
}

/// A glyph's advance in font units, with the font's kerning before the glyph `next`.
fn advance(face: &ttf_parser::Face, id: ttf_parser::GlyphId, next: Option<ttf_parser::GlyphId>) -> f64 {
    let mut a = f64::from(face.glyph_hor_advance(id).unwrap_or(0));
    if let (Some(next), Some(kern)) = (next, face.tables().kern) {
        for table in kern.subtables {
            if table.horizontal && !table.variable {
                if let Some(k) = table.glyphs_kerning(id, next) {
                    a += f64::from(k);
                    break;
                }
            }
        }
    }
    a
}

/// A segment of a glyph's outline: its points (the ends, and the control points between).
#[derive(Debug, Clone)]
enum Segment<P = (f32, f32)> {
    Line(P, P),
    Quadratic(P, P, P),
    Cubic(P, P, P, P),
}

impl<P: Copy> Segment<P> {
    fn map<Q>(&self, f: impl Fn(P) -> Q) -> Segment<Q> {
        match *self {
            Segment::Line(a, b) => Segment::Line(f(a), f(b)),
            Segment::Quadratic(a, b, c) => Segment::Quadratic(f(a), f(b), f(c)),
            Segment::Cubic(a, b, c, d) => Segment::Cubic(f(a), f(b), f(c), f(d)),
        }
    }

    fn start(&self) -> P {
        match *self {
            Segment::Line(a, _) | Segment::Quadratic(a, _, _) | Segment::Cubic(a, _, _, _) => a,
        }
    }

    fn end(&self) -> P {
        match *self {
            Segment::Line(_, b) | Segment::Quadratic(_, _, b) | Segment::Cubic(_, _, _, b) => b,
        }
    }
}

/// A glyph's outline as `ttf_parser` reads it: contours of segments.
#[derive(Default)]
struct Outline {
    contours: Vec<Vec<Segment>>,
    at: (f32, f32),
    start: (f32, f32),
}

impl ttf_parser::OutlineBuilder for Outline {
    fn move_to(&mut self, x: f32, y: f32) {
        self.contours.push(Vec::new());
        self.at = (x, y);
        self.start = (x, y);
    }

    fn line_to(&mut self, x: f32, y: f32) {
        self.push(Segment::Line(self.at, (x, y)));
    }

    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        self.push(Segment::Quadratic(self.at, (x1, y1), (x, y)));
    }

    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        self.push(Segment::Cubic(self.at, (x1, y1), (x2, y2), (x, y)));
    }

    fn close(&mut self) {
        if self.at != self.start {
            self.push(Segment::Line(self.at, self.start));
        }
    }
}

impl Outline {
    fn push(&mut self, s: Segment) {
        self.at = s.end();
        if let Some(c) = self.contours.last_mut() {
            c.push(s);
        }
    }
}

/// The faces of a glyph's contours: each contour inside an even number of the others bounds
/// a face, with those directly inside it as its holes.
fn glyph_faces(model: &mut Model, contours: &[Vec<Segment<Point>>], tolerance: f64) -> Result<Vec<Shape>> {
    // Each contour without its segments too short to matter, with a polygon along it.
    let contours: Vec<Vec<Segment<Point>>> = contours
        .iter()
        .map(|c| c.iter().filter(|s| s.start().distance(s.end()) > tolerance || !matches!(s, Segment::Line(..))).cloned().collect::<Vec<_>>())
        .filter(|c: &Vec<Segment<Point>>| !c.is_empty())
        .collect();
    let polygons: Vec<Vec<Point>> = contours.iter().map(|c| polygon(c)).collect();
    let areas: Vec<f64> = polygons.iter().map(|p| signed_area(p)).collect();
    let keep: Vec<usize> = (0..contours.len()).filter(|&i| areas[i].abs() > tolerance * tolerance).collect();
    // How many others hold each.
    let inside = |i: usize, j: usize| i != j && contains(&polygons[j], polygons[i][0]);
    let depth: Vec<usize> = (0..contours.len()).map(|i| keep.iter().filter(|&&j| inside(i, j)).count()).collect();
    let mut faces = Vec::new();
    for &i in keep.iter().filter(|&&i| depth[i] % 2 == 0) {
        // (Outer contours counter-clockwise, holes clockwise: the face facing +Z.)
        let outer = wire(model, &oriented(&contours[i], areas[i] > 0.0), tolerance)?;
        let face = crate::make::planar_face(model, outer)?;
        for &j in keep.iter().filter(|&&j| depth[j] == depth[i] + 1 && inside(j, i)) {
            let hole = wire(model, &oriented(&contours[j], areas[j] < 0.0), tolerance)?;
            model.add(face, hole)?;
        }
        faces.push(face);
    }
    Ok(faces)
}

/// A contour as it is (`keep`), or run the other way round.
fn oriented(contour: &[Segment<Point>], keep: bool) -> Vec<Segment<Point>> {
    if keep {
        return contour.to_vec();
    }
    contour
        .iter()
        .rev()
        .map(|s| match *s {
            Segment::Line(a, b) => Segment::Line(b, a),
            Segment::Quadratic(a, b, c) => Segment::Quadratic(c, b, a),
            Segment::Cubic(a, b, c, d) => Segment::Cubic(d, c, b, a),
        })
        .collect()
}

/// The closed wire of a contour: an edge per segment, between shared vertices.
fn wire(model: &mut Model, contour: &[Segment<Point>], tolerance: f64) -> Result<Shape> {
    let n = contour.len();
    let vertices: Vec<_> = contour.iter().map(|s| model.add_vertex(s.start(), tolerance)).collect();
    let mut coedges = Vec::with_capacity(n);
    for (k, s) in contour.iter().enumerate() {
        let (a, b) = (vertices[k], vertices[(k + 1) % n]);
        let end = model.point(Shape::from(b));
        let (curve, last) = match *s {
            Segment::Line(p, _) => {
                let d = end - p;
                let direction = Direction::try_from_vector(d).ok_or_else(|| Error::InvalidInput("a glyph's segment of no length".into()))?;
                (Curve::Line(Line::new(p, direction)), d.length())
            }
            Segment::Quadratic(p, c, _) => (Curve::Bezier(BezierCurve::new(vec![p, c, end])), 1.0),
            Segment::Cubic(p, c1, c2, _) => (Curve::Bezier(BezierCurve::new(vec![p, c1, c2, end])), 1.0),
        };
        let edge = model.add_edge(a, b, tolerance)?;
        model.set_edge_curve(edge.into(), curve, tolerance)?;
        model.set_edge_range(edge.into(), 0.0, last)?;
        coedges.push(model.add_coedge(edge, Orientation::Forward)?);
    }
    Ok(Shape::from(model.add_wire(&coedges)?))
}

/// Points along a contour (each curve sampled).
fn polygon(contour: &[Segment<Point>]) -> Vec<Point> {
    let mut points = Vec::new();
    for s in contour {
        match *s {
            Segment::Line(a, _) => points.push(a),
            Segment::Quadratic(a, b, c) => {
                let curve = BezierCurve::new(vec![a, b, c]);
                points.extend((0..8).map(|k| curve.point(k as f64 / 8.0)));
            }
            Segment::Cubic(a, b, c, d) => {
                let curve = BezierCurve::new(vec![a, b, c, d]);
                points.extend((0..8).map(|k| curve.point(k as f64 / 8.0)));
            }
        }
    }
    points
}

/// The area a polygon encloses in the XY plane, positive counter-clockwise.
fn signed_area(points: &[Point]) -> f64 {
    let n = points.len();
    0.5 * (0..n).map(|i| {
        let (p, q) = (points[i], points[(i + 1) % n]);
        p.x * q.y - q.x * p.y
    }).sum::<f64>()
}

/// Whether `p` lies inside the polygon (even-odd crossings).
fn contains(points: &[Point], p: Point) -> bool {
    let n = points.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (points[i], points[(i + 1) % n]);
        if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
            inside = !inside;
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A font the tests can rely on: one of the system's, if any (they are skipped otherwise).
    fn some_font() -> Option<Vec<u8>> {
        system_font("ttf", |_| true)
    }

    /// One of the system's fonts with the extension `extension` (and an 'O'), for which
    /// `wanted` holds.
    fn system_font(extension: &str, wanted: impl Fn(&ttf_parser::Face) -> bool) -> Option<Vec<u8>> {
        let dirs = ["/usr/share/fonts", "/usr/local/share/fonts", "/System/Library/Fonts", "C:\\Windows\\Fonts"];
        let mut stack: Vec<std::path::PathBuf> = dirs.iter().map(std::path::PathBuf::from).collect();
        while let Some(dir) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&dir) else { continue };
            for e in entries.flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case(extension)) {
                    if let Ok(bytes) = std::fs::read(&p) {
                        if ttf_parser::Face::parse(&bytes, 0).is_ok_and(|f| f.glyph_index('O').is_some() && f.glyph_index('H').is_some() && wanted(&f)) {
                            return Some(bytes);
                        }
                    }
                }
            }
        }
        None
    }

    #[test]
    fn letters_are_faces_with_holes() {
        let Some(font) = some_font() else { return };
        let mut model = Model::new();
        let t = text(&mut model, &font, "O", 10.0, HorizontalAlignment::Left, VerticalAlignment::Bottom).unwrap();
        let faces = model.unique_subshapes(t, NodeKind::Face);
        assert_eq!(faces.len(), 1);
        assert_eq!(model.unique_subshapes(t, NodeKind::Wire).len(), 2);
        assert!(crate::validity::check(&model, t).is_valid());
        let area = crate::properties::area(&model, t).mass;
        assert!(area > 0.0);
        // Facing +Z: extruded along +Z, the solid has a positive volume.
        let solid = crate::sweep::extrude(&mut model, faces[0], Vector::new(0.0, 0.0, 1.0)).unwrap().shape;
        assert!((crate::properties::volume(&model, solid).mass - area).abs() < 1e-6 * area);
        let normal = match model.surface(faces[0]) { Some(crate::geometry::surface::Surface::Plane(p)) => Vector::from(p.normal()), _ => panic!() };
        let normal = if faces[0].orientation == Orientation::Reversed { -normal } else { normal };
        assert!(normal.z > 0.99, "{normal:?}");
        // Facing +Z, on the baseline (bottom), starting past the pen (left).
        let b = crate::properties::bounding_box(&model, t);
        assert!(b.min.y > -0.5 && b.min.x >= 0.0 && b.max.y < 10.0);
    }

    /// A character the font lacks is set in the fallback for its script, as the fallback
    /// sets it; one no font has as the font's missing-glyph glyph, which a fallback for
    /// another script does not stand in for.
    #[test]
    fn characters_a_font_lacks_come_from_the_fallback_for_their_script() {
        let ideograph = '\u{590F}';
        let Some(font) = system_font("ttf", |f| f.glyph_index(ideograph).is_none()) else { return };
        let Some(cjk) = system_font("ttf", |f| f.glyph_index(ideograph).is_some()) else { return };
        let mut model = Model::new();
        let area = |model: &mut Model, font: &[u8], fallbacks: &[Fallback], s: &str| {
            let t = text_with_fallbacks(model, font, fallbacks, s, 10.0, HorizontalAlignment::Left, VerticalAlignment::Bottom).unwrap();
            crate::properties::area(model, t).mass
        };
        let fallback = [Fallback { script: Script::Cjk, font: &cjk, index: 0 }];
        let own = area(&mut model, &cjk, &[], &ideograph.to_string());
        assert!(own > 0.0);
        let borrowed = area(&mut model, &font, &fallback, &ideograph.to_string());
        assert!((borrowed - own).abs() < 1e-9 * own, "{borrowed} {own}");
        // Without a fallback (or with one for another script only), the font's .notdef.
        let notdef = area(&mut model, &font, &[], &ideograph.to_string());
        let other = [Fallback { script: Script::Korean, font: &cjk, index: 0 }];
        assert_eq!(area(&mut model, &font, &other, &ideograph.to_string()), notdef);
        assert_eq!(Script::of(ideograph), Script::Cjk);
        assert_eq!(Script::of('\u{3010}'), Script::Western);
    }

    /// A font with CFF outlines stands as high and low as its glyphs' bounding box: aligned
    /// at the top, its baseline that far down; its lines that box's height apart.
    #[test]
    fn cff_fonts_stand_as_their_bounding_box() {
        let Some(font) = system_font("otf", |f| f.tables().cff.is_some() && f.global_bounding_box().y_max > f.ascender()) else { return };
        let face = ttf_parser::Face::parse(&font, 0).unwrap();
        let scale = 10.0 / f64::from(face.units_per_em());
        let b = face.global_bounding_box();
        let mut model = Model::new();
        let top = |model: &mut Model, s: &str, v: VerticalAlignment| {
            let t = text(model, &font, s, 10.0, HorizontalAlignment::Left, v).unwrap();
            crate::properties::bounding_box(model, t).max.y
        };
        let cap = top(&mut model, "H", VerticalAlignment::Bottom);
        let at_top = top(&mut model, "H", VerticalAlignment::Top);
        assert!((at_top - (cap - f64::from(b.y_max) * scale)).abs() < 1e-6, "{at_top} {cap}");
        let two = top(&mut model, "H\nH", VerticalAlignment::Bottom);
        assert!((two - cap - f64::from(b.y_max - b.y_min) * scale).abs() < 1e-6, "{two} {cap}");
    }
}
