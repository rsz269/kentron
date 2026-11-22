//! Hidden-line drawings: a shape's sharp edges and outlines seen along a direction, each cut
//! into the stretches the viewer sees and the stretches the shape hides.
//!
//! As OCCT's exact hidden-line removal does (`HLRBRep_Algo`, read out by
//! `HLRBRep_HLRToShape`), for a parallel projection:
//!
//! * The data (`HLRBRep_ShapeToHLR`, `HLRBRep_Data::Update`): every face with its surface,
//!   every edge with the faces it bounds, an edge between two faces that meet tangent-plane
//!   continuously along it being a smooth one (`Rg1`, as `BRepLib::EncodeRegularity` finds
//!   it) and an edge a face meets itself along a seam; the scene's size (`myBigSize`, the
//!   widest of its sixteen min-max spans); each face's flags: whether it has outlines, is
//!   simple (no outline, an elementary surface), seen edge-on (`IsSide`), seen from behind
//!   (`Back`) and so hiding.
//! * The outlines (`HLRTopoBRep_DSFiller` through `Contap_Contour`): on a cylinder, cone or
//!   sphere, the lines or circle where the line of sight grazes the surface
//!   (`Contap_ContAna`), cut where the face's boundary arcs cross them, the arcs searched for
//!   where the surface's normal turns square to the line of sight (`Contap_ArcFunction`),
//!   and kept where they lie in the face (`Contap_Contour::LineConstructor`). A boundary
//!   arc lying on the contour all along (a seam the line of sight grazes) is an outline
//!   itself (`Contap_Restriction`).
//! * The hiding (`HLRBRep_InternalAlgo::HideAll`, `HLRBRep_Hider::Hide`): each edge and
//!   outline is cut where its projection meets the projection of a hiding face's boundary
//!   or outlines (the interferences `HLRBRep_Intersector` finds, at crossings and where an
//!   end of one touches the other), and each stretch between cuts is hidden where
//!   `HLRBRep_Data::Classify` finds a hiding face in front of its middle: the line of
//!   sight through the middle meets the face's surface nearer the viewer by more than
//!   `TolZ`, at a point the face's classifier does not find outside it. A simple face's
//!   own edges are not tested against it; another face's own lines are, with its larger
//!   `TolZ`.
//!
//! What is not ported: OCCT sorts the interferences into transitions (`HLRBRep_EdgeIList`,
//! `HLRBRep_EdgeBuilder`) to find which stretches lie inside a face's projection and
//! classifies only those (always, `IsSuspicion` being set); here every stretch is
//! classified, which `Classify` answers alike for a stretch outside the projection. Outlines
//! of other surfaces (a torus, a B-spline, which OCCT follows by walking) are refused, as
//! are perspective views, iso-lines and the smooth (`Rg1`, `RgN`) categories, which are
//! not drawn.

use crate::boolean::{PointClassifier, PointLocation};
use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Line, TrimmedCurve};
use crate::geometry::intersect::curve_surface;
use crate::geometry::surface::Surface;
use crate::geometry::{sample, ParametricCurve, ParametricSurface};
use crate::math::{Direction, Frame, Point, Point2, Transform, Transformable, Vector, Vector2};
use crate::topology::{EdgeId, FaceId, Model, NodeKind, Orientation, Shape};

/// What a drawn line is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineKind {
    /// A stretch of an edge between faces that meet at an angle, or of a free edge.
    Sharp,
    /// A stretch of where the line of sight grazes a curved face.
    Outline,
}

/// One stretch of a [`HiddenLines`] drawing: points on it in the view plane, every chord
/// within the tolerance asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct DrawnLine {
    pub kind: LineKind,
    pub points: Vec<Point2>,
}

/// A hidden-line drawing (see [`hidden_lines`]).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct HiddenLines {
    pub visible: Vec<DrawnLine>,
    pub hidden: Vec<DrawnLine>,
}

/// `BRepLib::EncodeRegularity`'s default angle (rad) within which two faces' normals along
/// an edge make it a smooth one.
const SMOOTH_ANGLE: f64 = 1e-10;
/// `HLRBRep_Data::myToler`.
const TOLER: f64 = 1e-5;

/// Draws `shape` seen along `-view.z()`, the viewer on `view`'s positive z side: each sharp
/// edge and each outline (where the line of sight grazes a curved face), cut into the
/// stretches seen and those hidden behind the shape's faces, as polylines in `view`'s
/// xy plane (a point `p` drawing at `((p - o)·x, (p - o)·y)`) whose vertices lie on the
/// stretches and whose chords lie within `tolerance` of them. Edges between faces meeting
/// smoothly, seams, degenerate edges and edges seen end on are not drawn.
///
/// Refuses a `tolerance` not positive, and a shape with a curved face other than a
/// cylinder, cone or sphere (whose outlines are not found yet).
#[doc(alias = "HLRBRep_Algo")]
#[doc(alias = "HLRBRep_HLRToShape")]
pub fn hidden_lines(model: &Model, shape: Shape, view: &Frame, tolerance: f64) -> Result<HiddenLines> {
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(Error::InvalidInput(format!("a drawing tolerance of {tolerance}")));
    }
    let view = View::new(view);
    let mut data = Data::load(model, shape, &view)?;
    data.outlines(model, &view)?;
    data.flag_faces(model, &view);
    data.hide_all(model, &view, tolerance)
}

/// The projection: `HLRAlgo_Projector` for a parallel view.
struct View {
    origin: Point,
    x: Vector,
    y: Vector,
    /// Toward the viewer.
    z: Vector,
}

impl View {
    fn new(frame: &Frame) -> View {
        View { origin: frame.origin(), x: Vector::from(frame.x()), y: Vector::from(frame.y()), z: Vector::from(frame.z()) }
    }

    fn flat(&self, p: Point) -> Point2 {
        let w = p - self.origin;
        Point2::new(w.dot(self.x), w.dot(self.y))
    }

    fn depth(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.z)
    }

    fn flat_vector(&self, v: Vector) -> Vector2 {
        Vector2::new(v.dot(self.x), v.dot(self.y))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    Sharp,
    /// Between faces meeting tangent-plane continuously (`Rg1Line`).
    Smooth,
    /// Where a face meets itself.
    Seam,
    Outline,
}

/// An edge or an outline (`HLRBRep_EdgeData` with its `HLRBRep_Curve`).
struct HlrLine {
    curve: Curve,
    first: f64,
    last: f64,
    role: Role,
    /// The faces it bounds or, an outline, lies on.
    faces: Vec<usize>,
    track: Track,
    /// Seen end on: its projection a point (`HLRBRep_EdgeData::Vertical`).
    vertical: bool,
}

/// A line's projection sampled: parameters, points and how far each chord may stray from
/// the projected curve.
struct Track {
    params: Vec<f64>,
    points: Vec<Point2>,
    sag: Vec<f64>,
    lo: Point2,
    hi: Point2,
    depth: (f64, f64),
}

impl HlrLine {
    fn new(curve: Curve, first: f64, last: f64, role: Role, faces: Vec<usize>, view: &View, tolerance: f64) -> HlrLine {
        let n = match &curve {
            Curve::Line(_) => 1,
            Curve::Circle(_) | Curve::Ellipse(_) => 128,
            Curve::BSpline(b) => (16 * b.knots().len()).clamp(128, 4096),
            _ => 256,
        };
        let params: Vec<f64> = (0..=n).map(|i| if i == n { last } else { first + (last - first) * i as f64 / n as f64 }).collect();
        let points3: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
        let points: Vec<Point2> = points3.iter().map(|&p| view.flat(p)).collect();
        let sag = (0..n)
            .map(|i| {
                let mid = view.flat(curve.point(0.5 * (params[i] + params[i + 1])));
                let chord = Point2::new(0.5 * (points[i].x + points[i + 1].x), 0.5 * (points[i].y + points[i + 1].y));
                2.0 * mid.distance(chord)
            })
            .collect();
        let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut depth = (f64::INFINITY, f64::NEG_INFINITY);
        for (p, q) in points.iter().zip(&points3) {
            (lo.x, lo.y, hi.x, hi.y) = (lo.x.min(p.x), lo.y.min(p.y), hi.x.max(p.x), hi.y.max(p.y));
            let z = view.depth(*q);
            depth = (depth.0.min(z), depth.1.max(z));
        }
        let vertical = hi.x - lo.x < tolerance && hi.y - lo.y < tolerance;
        HlrLine { curve, first, last, role, faces, track: Track { params, points, sag, lo, hi, depth }, vertical }
    }

    fn flat(&self, view: &View, t: f64) -> Point2 {
        view.flat(self.curve.point(t))
    }

    fn flat_d1(&self, view: &View, t: f64) -> (Point2, Vector2) {
        let (p, d) = self.curve.d1(t);
        (view.flat(p), view.flat_vector(d))
    }

    fn drawn(&self) -> bool {
        matches!(self.role, Role::Sharp | Role::Outline) && !self.vertical
    }
}

/// A face (`HLRBRep_FaceData` with its `HLRBRep_Surface`).
struct HlrFace {
    shape: Shape,
    /// The surface in the face's own coordinates, and placed.
    local: Surface,
    surface: Surface,
    classifier: PointClassifier,
    /// Its boundary edges and outlines, as indices into the lines.
    lines: Vec<usize>,
    with_outlines: bool,
    simple: bool,
    hiding: bool,
    /// Its lines' projected box and nearest depth, and a margin for their chords.
    lo: Point2,
    hi: Point2,
    depth: f64,
    margin: f64,
}

/// `HLRBRep_Data`: the faces and lines, and the scene's size.
struct Data {
    faces: Vec<HlrFace>,
    lines: Vec<HlrLine>,
    big: f64,
}

/// The surface's unit normal at `p`, a point of it; `None` where it has none.
fn normal_at(surface: &Surface, p: Point) -> Option<Vector> {
    let (u, v) = surface.parameters(p).unwrap_or_else(|| {
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, surface, None);
        (u, v)
    });
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    (n.length() > 1e-300).then(|| n.normalized())
}

impl Data {
    /// `HLRBRep_ShapeToHLR::Load`: the faces, and the edges with their regularity.
    fn load(model: &Model, shape: Shape, view: &View) -> Result<Data> {
        let mut faces = Vec::new();
        for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
            let local_shape = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
            let Some(local) = model.surface(local_shape).cloned() else { continue };
            let surface = local.transformed(&face.placement);
            let classifier = PointClassifier::new(model, face)?;
            faces.push(HlrFace { shape: face, local, surface, classifier, lines: Vec::new(), with_outlines: false, simple: false, hiding: false, lo: Point2::new(0.0, 0.0), hi: Point2::new(0.0, 0.0), depth: 0.0, margin: 0.0 });
        }
        // Each edge's uses by the faces (an edge placed twice is two edges).
        let mut uses: Vec<(EdgeId, Transform, Vec<usize>)> = Vec::new();
        for (f, face) in faces.iter().enumerate() {
            let id = FaceId::try_from(face.shape.node).expect("a face");
            for coedge in model.face_coedges(id) {
                let edge = model[coedge].edge;
                match uses.iter_mut().find(|u| u.0 == edge && u.1 == face.shape.placement) {
                    Some(u) => u.2.push(f),
                    None => uses.push((edge, face.shape.placement, vec![f])),
                }
            }
        }
        for edge in model.explore(shape, NodeKind::Edge).collect::<Vec<_>>() {
            let id = EdgeId::try_from(edge.node).expect("an edge");
            if !uses.iter().any(|u| u.0 == id && u.1 == edge.placement) {
                uses.push((id, edge.placement, Vec::new()));
            }
        }
        let mut lines = Vec::new();
        let tolerance = |e: EdgeId| model.tolerance(Shape::from(e)).max(crate::tolerance::LINEAR);
        for (edge, placement, used) in uses {
            let placed = Shape { placement, ..Shape::from(edge) };
            if model.is_degenerate(placed) {
                continue;
            }
            let Some((curve, first, last)) = model.edge_curve(placed) else { continue };
            let role = if used.len() >= 2 && used.iter().all(|&f| f == used[0]) {
                Role::Seam
            } else if used.len() == 2 && smooth(&faces[used[0]], &faces[used[1]], &curve, first, last) {
                Role::Smooth
            } else {
                Role::Sharp
            };
            let mut on: Vec<usize> = used.clone();
            on.dedup();
            let index = lines.len();
            for &f in &on {
                faces[f].lines.push(index);
            }
            lines.push(HlrLine::new(curve, first, last, role, on, view, tolerance(edge)));
        }
        let big = big_size(&lines, view);
        Ok(Data { faces, lines, big })
    }

    /// `HLRTopoBRep_DSFiller::Insert` with `Contap_Contour::PerformAna`: each face's outlines.
    fn outlines(&mut self, model: &Model, view: &View) -> Result<()> {
        let d = Direction::try_from_vector(view.z).expect("a view direction");
        for f in 0..self.faces.len() {
            let face = &self.faces[f];
            // `Contap_ContAna`: lines through a point along a direction, or a circle.
            let mut straight: Vec<(Point, Direction)> = Vec::new();
            let mut round: Option<Circle> = None;
            match &face.surface {
                Surface::Plane(_) => continue,
                Surface::Cylinder(c) => {
                    let axis = c.frame().z();
                    let normale = axis.cross(d);
                    if normale.length() > 1e-15 {
                        let n = normale.normalized();
                        let o = c.frame().origin();
                        straight.push((o + n * c.radius(), axis));
                        straight.push((o - n * c.radius(), axis));
                    }
                }
                Surface::Cone(c) => {
                    let frame = c.frame();
                    let tgtalpha = c.semi_angle().tan();
                    let coefcos = d.dot(frame.x());
                    let coefsin = d.dot(frame.y());
                    let coefcst = d.dot(frame.z()) * tgtalpha;
                    let norm1 = coefcos * coefcos + coefsin * coefsin;
                    if coefcst.abs() < norm1.sqrt() {
                        let prm = (norm1 - coefcst * coefcst).sqrt();
                        let cost0 = (coefcos * coefcst - coefsin * prm) / norm1;
                        let cost1 = (coefcos * coefcst + coefsin * prm) / norm1;
                        let sint0 = (coefcos * prm + coefsin * coefcst) / norm1;
                        let sint1 = (-coefcos * prm + coefsin * coefcst) / norm1;
                        let (x, y, z) = (Vector::from(frame.x()), Vector::from(frame.y()), Vector::from(frame.z()));
                        let apex = c.apex();
                        for (cos, sin) in [(cost0, sint0), (cost1, sint1)] {
                            if let Some(dir) = Direction::try_from_vector(x * cos + y * sin + z * (1.0 / tgtalpha)) {
                                straight.push((apex, dir));
                            }
                        }
                    }
                }
                Surface::Sphere(s) => {
                    let frame = s.frame();
                    let dir2 = if d.dot(frame.x()).abs() < 0.9999999999999 { d.cross(frame.x()) } else { d.cross(frame.y()) };
                    let x = Direction::try_from_vector(dir2).expect("square to the view");
                    round = Some(Circle::new(Frame::new(frame.origin(), d, x), s.radius()));
                }
                other => {
                    return Err(Error::Failed(format!(
                        "hidden lines: the outlines of a {} face are not found yet",
                        surface_name(other)
                    )))
                }
            }
            if straight.is_empty() && round.is_none() {
                continue;
            }
            // `Contap_TheSearch` with `Contap_ArcFunction`: where the boundary arcs cross the
            // outline.
            let on = 1e-7 * (1.0 + self.big);
            let mut crossings: Vec<Point> = Vec::new();
            // An arc lying on the contour all along is an outline itself (`Contap_Restriction`,
            // from `Contap_Contour::ProcessSegments`), its ends the outline's vertices.
            let mut restricted: Vec<usize> = Vec::new();
            for &l in &face.lines {
                let line = &self.lines[l];
                if on_contour(&face.surface, &line.curve, line.first, line.last, d) {
                    restricted.push(l);
                    for p in [line.curve.point(line.first), line.curve.point(line.last)] {
                        if !crossings.iter().any(|q| q.distance(p) <= on) {
                            crossings.push(p);
                        }
                    }
                    continue;
                }
                for p in arc_zeros(&face.surface, &line.curve, line.first, line.last, d) {
                    if !crossings.iter().any(|q| q.distance(p) <= on) {
                        crossings.push(p);
                    }
                }
            }
            // `Contap_Contour::LineConstructor`: the outlines between the crossings, kept where
            // they lie in the face.
            let inside = |p: Point| face.classifier.locate_on_faces(p) == PointLocation::Inside;
            let mut found: Vec<(Curve, f64, f64)> = Vec::new();
            for (origin, dir) in straight {
                let line = Line::new(origin, dir);
                let mut params: Vec<f64> = crossings
                    .iter()
                    .filter(|&&p| line_distance(&line, p) <= on)
                    .map(|&p| (p - origin).dot(Vector::from(dir)))
                    .collect();
                params.sort_by(f64::total_cmp);
                params.dedup_by(|a, b| (*a - *b).abs() <= on);
                for w in params.windows(2) {
                    let mid = 0.5 * (w[0] + w[1]);
                    if w[1] - w[0] > on && inside(origin + Vector::from(dir) * mid) {
                        found.push((Curve::Line(line), w[0], w[1]));
                    }
                }
            }
            if let Some(circle) = round {
                let curve = Curve::Circle(circle);
                let mut params: Vec<f64> = crossings
                    .iter()
                    .filter(|&&p| circle_distance(&circle, p) <= on)
                    .map(|&p| circle_parameter(&circle, p))
                    .collect();
                params.sort_by(f64::total_cmp);
                params.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
                if params.is_empty() {
                    if inside(curve.point(std::f64::consts::PI)) {
                        found.push((curve.clone(), 0.0, std::f64::consts::TAU));
                    }
                } else {
                    let n = params.len();
                    for i in 0..n {
                        let (a, b) = (params[i], if i + 1 < n { params[i + 1] } else { params[0] + std::f64::consts::TAU });
                        if b - a > 1e-12 && inside(curve.point(0.5 * (a + b))) {
                            found.push((curve.clone(), a, b));
                        }
                    }
                }
            }
            for (curve, a, b) in found {
                let index = self.lines.len();
                self.lines.push(HlrLine::new(curve, a, b, Role::Outline, vec![f], view, crate::tolerance::LINEAR));
                self.faces[f].lines.push(index);
                self.faces[f].with_outlines = true;
            }
            for l in restricted {
                // `HLRTopoBRep_Data::AddOutL`: drawn as an outline, a seam or smooth edge too.
                if matches!(self.lines[l].role, Role::Seam | Role::Smooth) {
                    self.lines[l].role = Role::Outline;
                }
                self.faces[f].with_outlines = true;
            }
        }
        Ok(())
    }

    /// `HLRBRep_Data::Update`, the faces: whether simple, seen edge-on, from behind, hiding.
    fn flag_faces(&mut self, model: &Model, view: &View) {
        for face in &mut self.faces {
            let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
            let (mut depth, mut margin) = (f64::NEG_INFINITY, 1e-7 * (1.0 + self.big));
            for &l in &face.lines {
                let t = &self.lines[l].track;
                (lo.x, lo.y, hi.x, hi.y) = (lo.x.min(t.lo.x), lo.y.min(t.lo.y), hi.x.max(t.hi.x), hi.y.max(t.hi.y));
                depth = depth.max(t.depth.1);
                margin = margin.max(t.sag.iter().copied().fold(0.0, f64::max));
            }
            (face.lo, face.hi, face.depth, face.margin) = (lo, hi, depth, margin);
            let elementary = matches!(face.surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_));
            face.simple = !face.with_outlines && elementary;
            let side = match &face.surface {
                Surface::Plane(p) => Vector::from(p.normal()).dot(view.z).abs() < TOLER * 10.0,
                Surface::Cylinder(c) => {
                    let a = Vector::from(c.frame().z());
                    a.dot(view.x).hypot(a.dot(view.y)) < TOLER * 10.0
                }
                _ => false,
            };
            face.hiding = if side {
                false
            } else if face.with_outlines {
                true
            } else {
                // `Back` by the surface's normal at the middle of the first edge where it is
                // not edge-on; a face seen edge-on all round is a side one.
                let mut back = None;
                for &l in &face.lines {
                    let line = &self.lines[l];
                    if line.role == Role::Seam {
                        continue;
                    }
                    let p = line.curve.point(0.5 * (line.first + line.last));
                    let local = p.transformed(&face.shape.placement.inverse());
                    let Some(n) = normal_at(&face.local, local) else { continue };
                    let r = n.transformed(&face.shape.placement).dot(view.z);
                    if r.abs() > TOLER * 10.0 {
                        back = Some(r < 0.0);
                        break;
                    }
                }
                match back {
                    None => false,
                    Some(back) => {
                        let in_solid = model.explore_up(face.shape.node).any(|n| n.kind() == NodeKind::Solid);
                        if !in_solid {
                            true
                        } else if face.shape.orientation == Orientation::Reversed {
                            back
                        } else {
                            !back
                        }
                    }
                }
            };
        }
    }

    /// `HLRBRep_InternalAlgo::HideAll` and `HLRBRep_HLRToShape`: each drawn line cut where
    /// the hiding faces' contours meet it, each stretch classified, and drawn.
    fn hide_all(&self, model: &Model, view: &View, tolerance: f64) -> Result<HiddenLines> {
        let mut out = HiddenLines::default();
        let touch = 1e-7 * (1.0 + self.big);
        for (e, line) in self.lines.iter().enumerate() {
            if !line.drawn() {
                continue;
            }
            // The interferences: crossings with every other line's projection, and the ends of
            // other lines lying on it.
            let mut cuts = vec![line.first, line.last];
            for (k, other) in self.lines.iter().enumerate() {
                if k == e || other.vertical || !boxes_meet(&line.track, &other.track, touch) {
                    continue;
                }
                cuts.extend(crossings(line, other, view, self.big));
                for t in [other.first, other.last] {
                    if let Some(s) = touching(line, other.flat(view, t), view, touch) {
                        cuts.push(s);
                    }
                }
            }
            let span = line.last - line.first;
            cuts.retain(|&t| t >= line.first && t <= line.last);
            cuts.sort_by(f64::total_cmp);
            cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * span.abs().max(1.0));
            // Each stretch hidden or seen, joined with its neighbours of the same state.
            let mut stretches: Vec<(f64, f64, bool)> = Vec::new();
            for w in cuts.windows(2) {
                let hidden = self.faces.iter().enumerate().any(|(f, face)| self.hides(view, f, face, e, 0.5 * (w[0] + w[1])));
                match stretches.last_mut() {
                    Some(last) if last.2 == hidden => last.1 = w[1],
                    _ => stretches.push((w[0], w[1], hidden)),
                }
            }
            let kind = if line.role == Role::Outline { LineKind::Outline } else { LineKind::Sharp };
            for (a, b, hidden) in stretches {
                if !(b > a) {
                    continue;
                }
                let piece = Curve::Trimmed(TrimmedCurve::new(line.curve.clone(), a, b));
                let sampled = sample::quasi_uniform_deflection(&piece, tolerance)
                    .ok_or_else(|| Error::Failed("hidden lines: a stretch could not be sampled".into()))?;
                let points = sampled.points.into_iter().map(|p| view.flat(p)).collect();
                let drawn = DrawnLine { kind, points };
                if hidden {
                    out.hidden.push(drawn);
                } else {
                    out.visible.push(drawn);
                }
            }
        }
        Ok(out)
    }

    /// `HLRBRep_Data::Classify` of line `e` at `t` against face `f`: whether the face hides
    /// the point.
    fn hides(&self, view: &View, f: usize, face: &HlrFace, e: usize, t: f64) -> bool {
        if !face.hiding {
            return false;
        }
        let line = &self.lines[e];
        let own = line.faces.contains(&f);
        // A simple face's own edges are jumped (`HLRBRep_Data::InitEdge`).
        if own && face.simple {
            return false;
        }
        let p = line.curve.point(t);
        let tol_z = if !own {
            self.big * 1e-6
        } else if line.role == Role::Outline {
            self.big * 1e-2
        } else {
            self.big * 1e-3
        };
        // The rejections by the face's min-max: off its projection, or in front of it all.
        let q = view.flat(p);
        let m = face.margin;
        if q.x < face.lo.x - m || q.x > face.hi.x + m || q.y < face.lo.y - m || q.y > face.hi.y + m || view.depth(p) > face.depth + m {
            return false;
        }
        // The line of sight through the point, met by the face's surface nearer the viewer
        // (`HLRBRep_InterCSurf`), at a point the face's classifier does not find outside.
        let toward = Direction::try_from_vector(view.z).expect("a view direction");
        let sight = Line::new(p, toward);
        let curve = match face.surface {
            Surface::Plane(_) | Surface::Cylinder(_) | Surface::Sphere(_) => Curve::Line(sight),
            _ => Curve::Trimmed(TrimmedCurve::new(Curve::Line(sight), -2.0 * self.big - 1.0, 2.0 * self.big + 1.0)),
        };
        let met = curve_surface(&curve, &face.surface, crate::tolerance::LINEAR);
        met.points.iter().any(|h| h.t > tol_z && face.classifier.locate_on_faces(h.point) != PointLocation::Outside)
    }
}

fn surface_name(s: &Surface) -> &'static str {
    match s {
        Surface::Plane(_) => "plane",
        Surface::Cylinder(_) => "cylinder",
        Surface::Cone(_) => "cone",
        Surface::Sphere(_) => "sphere",
        Surface::Torus(_) => "torus",
        Surface::Bezier(_) => "Bézier",
        Surface::BSpline(_) => "B-spline",
        Surface::Extrusion(_) => "extrusion",
        Surface::Revolution(_) => "revolution",
        Surface::Offset(_) => "offset",
        Surface::Trimmed(_) => "trimmed",
    }
}

/// The widest of the sixteen min-max spans of the lines' projections (`myBigSize`, as
/// `HLRAlgo::UpdateMinMax` measures them).
fn big_size(lines: &[HlrLine], view: &View) -> f64 {
    let mut lo = [f64::INFINITY; 15];
    let mut hi = [f64::NEG_INFINITY; 15];
    for line in lines {
        for (p, &t) in line.track.points.iter().zip(&line.track.params) {
            let mut d = [0.0; 15];
            for k in 0..7 {
                let (s, c) = (k as f64 * std::f64::consts::PI / 14.0).sin_cos();
                d[2 * k] = c * p.x + s * p.y;
                d[2 * k + 1] = s * p.x - c * p.y;
            }
            d[14] = view.depth(line.curve.point(t));
            for i in 0..15 {
                lo[i] = lo[i].min(d[i]);
                hi[i] = hi[i].max(d[i]);
            }
        }
    }
    (0..15).map(|i| hi[i] - lo[i]).filter(|d| d.is_finite()).fold(0.0, f64::max)
}

fn line_distance(line: &Line, p: Point) -> f64 {
    let w = p - line.origin;
    let d = Vector::from(line.direction);
    (w - d * w.dot(d)).length()
}

fn circle_distance(circle: &Circle, p: Point) -> f64 {
    let f = circle.frame();
    let w = p - f.origin();
    let n = Vector::from(f.z());
    let h = w.dot(n);
    let r = (w - n * h).length();
    (r - circle.radius()).hypot(h)
}

fn circle_parameter(circle: &Circle, p: Point) -> f64 {
    let f = circle.frame();
    let w = p - f.origin();
    let a = w.dot(Vector::from(f.y())).atan2(w.dot(Vector::from(f.x())));
    if a < 0.0 {
        a + std::f64::consts::TAU
    } else {
        a
    }
}

/// `Contap_ArcFunction` searched along a boundary arc: the points where the surface's
/// normal turns square to `d`.
fn arc_zeros(surface: &Surface, curve: &Curve, first: f64, last: f64, d: Direction) -> Vec<Point> {
    let d = Vector::from(d);
    let g = |t: f64| normal_at(surface, curve.point(t)).map(|n| n.dot(d));
    let n = 64;
    let at = |i: usize| if i == n { last } else { first + (last - first) * i as f64 / n as f64 };
    let mut out = Vec::new();
    let mut prev: Option<(f64, f64)> = None;
    for i in 0..=n {
        let t = at(i);
        let Some(v) = g(t) else {
            prev = None;
            continue;
        };
        if v == 0.0 {
            out.push(curve.point(t));
        } else if let Some((pt, pv)) = prev {
            if pv != 0.0 && (pv < 0.0) != (v < 0.0) {
                let (mut a, mut b, mut ga) = (pt, t, pv);
                for _ in 0..200 {
                    let m = 0.5 * (a + b);
                    if m <= a.min(b) || m >= a.max(b) {
                        break;
                    }
                    let Some(gm) = g(m) else { break };
                    if gm == 0.0 {
                        (a, b) = (m, m);
                        break;
                    }
                    if (gm < 0.0) == (ga < 0.0) {
                        (a, ga) = (m, gm);
                    } else {
                        b = m;
                    }
                }
                out.push(curve.point(0.5 * (a + b)));
            }
        }
        prev = Some((t, v));
    }
    out
}

/// Whether the surface's normal is square to `d` all along a boundary arc (to within
/// [`SMOOTH_ANGLE`]), the arc lying on the contour.
fn on_contour(surface: &Surface, curve: &Curve, first: f64, last: f64, d: Direction) -> bool {
    let d = Vector::from(d);
    (0..=8).all(|i| {
        let t = first + (last - first) * i as f64 / 8.0;
        normal_at(surface, curve.point(t)).is_some_and(|n| n.dot(d).abs() <= SMOOTH_ANGLE)
    })
}

fn boxes_meet(a: &Track, b: &Track, margin: f64) -> bool {
    let sa = a.sag.iter().copied().fold(0.0, f64::max);
    let sb = b.sag.iter().copied().fold(0.0, f64::max);
    let m = margin + sa + sb;
    !(a.lo.x - m > b.hi.x || b.lo.x - m > a.hi.x || a.lo.y - m > b.hi.y || b.lo.y - m > a.hi.y)
}

/// The parameters on `a` where the projections of `a` and `b` cross (`HLRBRep_Intersector`):
/// close pairs of their sampled chords, refined by Newton's method on the curves.
fn crossings(a: &HlrLine, b: &HlrLine, view: &View, big: f64) -> Vec<f64> {
    let (ta, tb) = (&a.track, &b.track);
    let mut out = Vec::new();
    let slack = 1e-9 * (1.0 + big);
    for i in 0..ta.points.len() - 1 {
        let (p0, p1) = (ta.points[i], ta.points[i + 1]);
        for j in 0..tb.points.len() - 1 {
            let (q0, q1) = (tb.points[j], tb.points[j + 1]);
            let margin = ta.sag[i] + tb.sag[j] + slack;
            if p0.x.min(p1.x) - margin > q0.x.max(q1.x)
                || q0.x.min(q1.x) - margin > p0.x.max(p1.x)
                || p0.y.min(p1.y) - margin > q0.y.max(q1.y)
                || q0.y.min(q1.y) - margin > p0.y.max(p1.y)
            {
                continue;
            }
            let (u, v, gap) = segment_closest(p0, p1, q0, q1);
            if gap > margin {
                continue;
            }
            let s = ta.params[i] + u * (ta.params[i + 1] - ta.params[i]);
            let t = tb.params[j] + v * (tb.params[j + 1] - tb.params[j]);
            if let Some(s) = refine_crossing(a, b, s, t, view, big) {
                out.push(s);
            }
        }
    }
    out
}

/// Newton's method on `a(s) = b(t)` in the view plane, within both ranges; `None` where it
/// does not converge (the curves tangent there, or not meeting).
fn refine_crossing(a: &HlrLine, b: &HlrLine, mut s: f64, mut t: f64, view: &View, big: f64) -> Option<f64> {
    let within = |x: f64, lo: f64, hi: f64| x >= lo.min(hi) - 1e-12 * (hi - lo).abs() && x <= lo.max(hi) + 1e-12 * (hi - lo).abs();
    for _ in 0..50 {
        let (pa, da) = a.flat_d1(view, s);
        let (pb, db) = b.flat_d1(view, t);
        let f = pa - pb;
        let det = db.cross(da);
        if det.abs() <= 1e-14 * da.length() * db.length() {
            return None;
        }
        // da·Δs − db·Δt = −f
        let ds = (f.x * db.y - f.y * db.x) / det;
        let dt = (f.x * da.y - f.y * da.x) / det;
        s += ds;
        t += dt;
        if !s.is_finite() || !t.is_finite() {
            return None;
        }
        if ds.abs() <= 1e-15 * (1.0 + s.abs()) && dt.abs() <= 1e-15 * (1.0 + t.abs()) {
            break;
        }
    }
    let gap = a.flat(view, s).distance(b.flat(view, t));
    (gap <= 1e-9 * (1.0 + big) && within(s, a.first, a.last) && within(t, b.first, b.last)).then(|| s.clamp(a.first.min(a.last), a.first.max(a.last)))
}

/// The parameter of the point of `a`'s projection nearest `q`, where that is within
/// `reach` of it.
fn touching(a: &HlrLine, q: Point2, view: &View, reach: f64) -> Option<f64> {
    let tr = &a.track;
    let mut best = (f64::INFINITY, a.first);
    for i in 0..tr.points.len() - 1 {
        let (p0, p1) = (tr.points[i], tr.points[i + 1]);
        let d = p1 - p0;
        let len = d.length_squared();
        let u = if len > 0.0 { ((q - p0).dot(d) / len).clamp(0.0, 1.0) } else { 0.0 };
        let gap = (p0 + d * u).distance(q) - tr.sag[i];
        if gap < best.0 {
            best = (gap, tr.params[i] + u * (tr.params[i + 1] - tr.params[i]));
        }
    }
    if best.0 > reach {
        return None;
    }
    let (lo, hi) = (a.first.min(a.last), a.first.max(a.last));
    let mut t = best.1;
    for _ in 0..50 {
        let (p, d) = a.flat_d1(view, t);
        let len = d.length_squared();
        if len <= 0.0 {
            break;
        }
        let step = (q - p).dot(d) / len;
        let next = (t + step).clamp(lo, hi);
        let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
        t = next;
        if done {
            break;
        }
    }
    (a.flat(view, t).distance(q) <= reach).then_some(t)
}

/// The closest points of two segments, as fractions along each, and their distance.
fn segment_closest(p0: Point2, p1: Point2, q0: Point2, q1: Point2) -> (f64, f64, f64) {
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.length_squared(), d2.length_squared(), d2.dot(r));
    let (s, t) = if a <= 1e-300 && e <= 1e-300 {
        (0.0, 0.0)
    } else if a <= 1e-300 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(r);
        if e <= 1e-300 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            let mut s = if denom > 0.0 { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
            let mut t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
            (s, t)
        }
    };
    let gap = (p0 + d1 * s).distance(q0 + d2 * t);
    (s, t, gap)
}

/// Whether two faces meet tangent-plane continuously along an edge (`BRepLib::
/// ContinuityOfFaces` to G1): their outward normals within [`SMOOTH_ANGLE`] at 21 points
/// along it, short of its ends by a hundredth of its range.
fn smooth(a: &HlrFace, b: &HlrFace, curve: &Curve, first: f64, last: f64) -> bool {
    let outward = |face: &HlrFace, p: Point| {
        let inverse = face.shape.placement.inverse();
        let n = normal_at(&face.local, p.transformed(&inverse))?.transformed(&face.shape.placement);
        Some(if face.shape.orientation == Orientation::Reversed { -n } else { n })
    };
    let eps = (last - first) / 100.0;
    let (f, l) = (first + eps, last - eps);
    for i in 0..=20 {
        let p = curve.point(f + (l - f) * i as f64 / 20.0);
        let (Some(n1), Some(n2)) = (outward(a, p), outward(b, p)) else { return false };
        if n1.cross(n2).length() > SMOOTH_ANGLE || n1.dot(n2) < 0.0 {
            return false;
        }
    }
    true
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Cuboid, Cylinder, Sphere};
    use std::f64::consts::PI;

    const TOLERANCE: f64 = 0.01;

    fn length(lines: &[DrawnLine], kind: LineKind) -> f64 {
        lines.iter().filter(|l| l.kind == kind).flat_map(|l| l.points.windows(2).map(|w| w[0].distance(w[1]))).sum()
    }

    /// The perimeter of an ellipse of semi-axes `a` and `b`, by the midpoint rule.
    fn perimeter(a: f64, b: f64) -> f64 {
        let n = 200_000;
        (0..n).map(|i| (i as f64 + 0.5) * 2.0 * PI / n as f64).map(|t| (a * t.sin()).hypot(b * t.cos())).sum::<f64>() * 2.0 * PI / n as f64
    }

    /// Looking down at `elevation` above the horizontal, toward +y: x across, the viewer
    /// at -y and above.
    fn tilted(elevation: f64) -> Frame {
        let toward = Direction::new(0.0, -elevation.cos(), elevation.sin());
        Frame::new(Point::ORIGIN, toward, Direction::X)
    }

    #[test]
    fn a_box_seen_from_a_corner_shows_nine_edges_and_hides_the_three_at_the_far_corner() {
        let mut model = Model::new();
        let (a, b, c) = (2.0, 3.0, 4.0);
        let solid = Cuboid::new(a, b, c).build(&mut model).unwrap();
        let toward = Direction::new(1.0, 1.0, 1.0);
        let view = Frame::new(Point::ORIGIN, toward, Direction::new(1.0, -1.0, 0.0));
        let drawn = hidden_lines(&model, solid, &view, TOLERANCE).unwrap();
        assert_eq!((drawn.visible.len(), drawn.hidden.len()), (9, 3));
        // Each edge projects to √(2/3) of its length; the far corner, the origin, draws at 0.
        let k = (2.0f64 / 3.0).sqrt();
        assert!((length(&drawn.hidden, LineKind::Sharp) - (a + b + c) * k).abs() < 1e-9);
        assert!((length(&drawn.visible, LineKind::Sharp) - 3.0 * (a + b + c) * k).abs() < 1e-9);
        for line in &drawn.hidden {
            let ends = [line.points[0], *line.points.last().unwrap()];
            assert!(ends.iter().any(|p| p.distance(Point2::new(0.0, 0.0)) < 1e-12), "{ends:?}");
        }
    }

    #[test]
    fn a_cylinder_seen_from_above_hides_the_far_half_of_its_bottom_rim() {
        let mut model = Model::new();
        let (r, h, elevation) = (2.0, 5.0, 0.6f64);
        let solid = Cylinder::new(r, h).build(&mut model).unwrap();
        let drawn = hidden_lines(&model, solid, &tilted(elevation), TOLERANCE).unwrap();
        // A rim draws as the ellipse of semi-axes r and r·sin(elevation); its far half, y > 0,
        // up the drawing, is hidden between the outlines' feet at (±r, 0).
        let minor = r * elevation.sin();
        let rim = perimeter(r, minor);
        assert!(drawn.hidden.iter().all(|l| l.kind == LineKind::Sharp));
        for p in drawn.hidden.iter().flat_map(|l| &l.points) {
            assert!(((p.x / r).powi(2) + (p.y / minor).powi(2) - 1.0).abs() < 1e-12 && p.y > -1e-12, "{p:?}");
        }
        let hidden = length(&drawn.hidden, LineKind::Sharp);
        assert!(hidden < rim / 2.0 && hidden > rim / 2.0 * (1.0 - TOLERANCE / minor), "{hidden} of {}", rim / 2.0);
        // Seen: the top rim, the near half of the bottom one, and two outlines x = ±r from
        // the bottom rim's to the top rim's, h·cos(elevation) up.
        let seen = length(&drawn.visible, LineKind::Sharp);
        assert!(seen < 1.5 * rim && seen > 1.5 * rim * (1.0 - TOLERANCE / minor), "{seen} of {}", 1.5 * rim);
        let outlines: Vec<&DrawnLine> = drawn.visible.iter().filter(|l| l.kind == LineKind::Outline).collect();
        assert_eq!(outlines.len(), 2);
        for line in outlines {
            assert!(line.points.iter().all(|p| (p.x.abs() - r).abs() < 1e-12));
            let (lo, hi) = line.points.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p.y), hi.max(p.y)));
            assert!(lo.abs() < 1e-12 && (hi - h * elevation.cos()).abs() < 1e-12, "{lo} {hi}");
        }
    }

    #[test]
    fn a_ball_draws_as_one_seen_circle() {
        let mut model = Model::new();
        let centre = Point::new(1.0, -2.0, 0.5);
        let solid = Sphere::new(3.0).with_frame(Frame::from_z(centre, Direction::Z)).build(&mut model).unwrap();
        let view = Frame::new(Point::ORIGIN, Direction::new(0.3, -0.5, 0.8), Direction::new(1.0, 0.0, -0.375));
        let drawn = hidden_lines(&model, solid, &view, TOLERANCE).unwrap();
        assert!(drawn.hidden.is_empty());
        assert!(drawn.visible.iter().all(|l| l.kind == LineKind::Outline));
        let w = centre - Point::ORIGIN;
        let c = Point2::new(w.dot(Vector::from(view.x())), w.dot(Vector::from(view.y())));
        assert!(drawn.visible.iter().flat_map(|l| &l.points).all(|p| (p.distance(c) - 3.0).abs() < 1e-12));
        let seen = length(&drawn.visible, LineKind::Outline);
        // A chord of sagitta s shortens its arc by about s/(3ρ) of it.
        assert!(seen < 6.0 * PI && seen > 6.0 * PI * (1.0 - TOLERANCE / 3.0), "{seen}");
    }

    #[test]
    fn a_holes_outlines_are_hidden_by_the_block_round_it() {
        let mut model = Model::new();
        let (side, height, radius, elevation) = (10.0, 2.0, 1.5, 0.9f64);
        let block = Cuboid::new(side, side, height).build(&mut model).unwrap();
        let drill = Cylinder::new(radius, height + 2.0).with_frame(Frame::from_z(Point::new(5.0, 5.0, -1.0), Direction::Z)).build(&mut model).unwrap();
        let drilled = crate::boolean::difference(&mut model, block, drill).unwrap();
        let drawn = hidden_lines(&model, drilled, &tilted(elevation), TOLERANCE).unwrap();
        // The line of sight from the hole's wall where it grazes it runs into the block, out
        // through its top: the two outlines, x = 5 ± r, hidden their whole height.
        assert!(length(&drawn.visible, LineKind::Outline) < 1e-4);
        let outlines: Vec<&DrawnLine> = drawn.hidden.iter().filter(|l| l.kind == LineKind::Outline).collect();
        assert_eq!(outlines.len(), 2);
        assert!(outlines.iter().flat_map(|l| &l.points).all(|p| ((p.x - 5.0).abs() - radius).abs() < 1e-12));
        assert!((length(&drawn.hidden, LineKind::Outline) - 2.0 * height * elevation.cos()).abs() < 1e-4);
    }

    #[test]
    fn a_tolerance_not_positive_is_refused() {
        let mut model = Model::new();
        let solid = Cuboid::new(1.0, 1.0, 1.0).build(&mut model).unwrap();
        assert!(matches!(hidden_lines(&model, solid, &tilted(0.5), 0.0), Err(Error::InvalidInput(_))));
    }
}
