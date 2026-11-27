//! Sheets thickened into solids: the space a face or an open shell sweeps as it moves along its
//! normals.
//!
//! A sheet thickened by a distance is every point p + s·n(p), for p on the sheet, n its unit
//! normal and s from 0 to the distance (to the other side for a negative one), with the gaps
//! its moved faces leave where they part at creases closed as the join says. Its boundary lies
//! on the sheet itself, on the sheet's offset (its faces moved one by one as [`super::shape`]
//! moves them, their creases rounded joined by arcs), on the walls its free edges sweep (the normals along each, from the
//! edge to its image) and, at a free edge's end where the offset is rounded, on the sector
//! between the end's vertex and the round's end. Those make one closed surface, oriented as the
//! boundary of the sweep: the sheet against the way it moves, the offset with it, the walls
//! outwards. Where the sweep folds over itself the surface passes through itself: a face moved
//! past its centre of curvature turns inside out (its offset still the image of its points, now
//! facing the other way), and a wall along it folds where its normals cross. A face moved exactly
//! to its centre or axis has no offset: its free edges sweep walls to the centre or axis.
//!
//! The surface's faces split space into cells (`boolean::enclosed_solids`). How many times the
//! surface winds round a cell is how many times the sweep covers it, each covering counted by
//! the way it turns: the solid is the cells it winds round, merged (the space inside a closed
//! sheet, which the sweep moves away from, it winds round no times).

use crate::error::{Error, Result};
use crate::geometry::surface::{OffsetSurface, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::hashing::FastMap;
use crate::math::{Point, Point2, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::Join;

/// The solid `sheet` (a face, or a shell of faces that is not closed, or a closed one) sweeps
/// moving `distance` along its faces' normals (against them for a negative distance), the
/// gaps the moved faces leave where they part joined as `join` says: by rounds about the
/// creases joined by arcs, by the moved faces extended until they meet joined by
/// intersection. Offsets past a face's centre of curvature (a bend thickened inside beyond
/// its radius) and offsets moving faces to their centre or axis (a hemisphere thickened to its
/// centre) give the space swept as well.
#[doc(alias = "BRepOffset_MakeOffset::Initialize(thickening)")]
#[doc(alias = "BRepOffsetAPI_MakeThickSolid")]
pub fn thicken(model: &mut Model, sheet: Shape, distance: f64, join: Join) -> Result<Shape> {
    if !(distance != 0.0 && distance.is_finite()) {
        return Err(Error::InvalidInput("a sheet thickened by no distance".into()));
    }
    let faces: Vec<Shape> = model.explore(sheet, NodeKind::Face).collect();
    if faces.is_empty() {
        return Err(Error::InvalidInput("a thickening of no faces".into()));
    }
    let shell = if sheet.kind() == NodeKind::Shell {
        sheet
    } else {
        let shell = model.new_empty(NodeKind::Shell);
        for &f in &faces {
            model.add(shell, f)?;
        }
        shell
    };
    let scale = crate::properties::bounding_box(model, shell).size().length().max(distance.abs()).max(1.0);
    let tolerance = 1e-7 * scale;
    let free = free_edges(model, shell);
    // Faces whose offset collapses onto their centre or axis.
    let mut foci = Vec::with_capacity(faces.len());
    for &f in &faces {
        foci.push(focus(model, f, distance)?);
    }
    let collapsing = foci.iter().filter(|f| f.is_some()).count();
    let mut surface: Vec<Shape> = faces.iter().map(|&f| if distance > 0.0 { f.reversed() } else { f }).collect();
    let mut builder = Walls { rulings: FastMap::default(), tolerance };
    if collapsing == 0 {
        // (Its faces moved one by one, those moved past their centres of curvature turned
        // inside out and kept so: the image of the sheet's points, whatever space that
        // bounds. The offset's other ways, finding the points within the distance, give the
        // boundary of a different space.)
        let (moved, _) = super::shape::moved_and_rounded(model, shell, distance, &[], join)?;
        let moved_faces: Vec<Shape> = model.explore(moved, NodeKind::Face).collect();
        surface.extend(moved_faces.iter().map(|&f| if distance < 0.0 { f.reversed() } else { f }));
        surface.extend(builder.to_offset(model, &free, moved, distance)?);
    } else if collapsing == faces.len() {
        let foci: Vec<Focus> = foci.into_iter().flatten().collect();
        surface.extend(builder.to_foci(model, &free, &faces, &foci, distance)?);
    } else {
        return Err(Error::InvalidInput("a thickening moving some of a sheet's faces onto their centres or axes, and others not".into()));
    }
    swept_cells(model, &surface, tolerance)
}

/// A free edge of a sheet: the edge run as its face uses it, and the face.
#[derive(Clone, Copy)]
struct FreeEdge {
    edge: Shape,
    face: Shape,
}

/// The edges of `shape` that bound it: each used by one face once (not shared by two faces, nor
/// a seam), not degenerate; run as their faces use them.
fn free_edges(model: &Model, shape: Shape) -> Vec<FreeEdge> {
    let mut uses: FastMap<NodeId, (usize, FreeEdge)> = FastMap::default();
    let mut order = Vec::new();
    for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
        for edge in model.explore(face, NodeKind::Edge).collect::<Vec<_>>() {
            let entry = uses.entry(edge.node).or_insert_with(|| {
                order.push(edge.node);
                (0, FreeEdge { edge, face })
            });
            entry.0 += 1;
        }
    }
    order.into_iter().filter_map(|n| uses.get(&n).copied()).filter(|&(n, e)| n == 1 && !model.is_degenerate(e.edge)).map(|(_, e)| e).collect()
}

/// The unit normal of `face` (as it is used) at its point nearest `p`.
fn normal_at(model: &Model, face: Shape, p: Point) -> Result<Vector> {
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face of the sheet without a surface".into()))?.transformed(&face.placement);
    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &surface, None);
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    if n.length() <= 1e-300 {
        return Err(Error::InvalidInput("a sheet's normal where its surface has none".into()));
    }
    Ok(n.normalized() * if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 })
}

/// Where a face moved by a distance collapses: a sphere's centre, a cylinder's axis.
#[derive(Clone, Copy)]
enum Focus {
    Centre(Point),
    Axis(Point, Vector),
}

impl Focus {
    /// The point of the focus that `p` (on the face) moves to.
    fn image(&self, p: Point) -> Point {
        match *self {
            Focus::Centre(c) => c,
            Focus::Axis(o, z) => o + z * (p - o).dot(z),
        }
    }
}

/// Where `face` moved `distance` along its normal collapses, if it does: a sphere moved by its
/// radius to its centre, a cylinder to its axis. (Other faces whose offsets collapse are not
/// thickened.)
fn focus(model: &Model, face: Shape, distance: f64) -> Result<Option<Focus>> {
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face of the sheet without a surface".into()))?.transformed(&face.placement);
    let along = if face.orientation == Orientation::Reversed { -distance } else { distance };
    let basis = match &surface {
        Surface::Trimmed(t) => t.basis().clone(),
        s => s.clone(),
    };
    let elementary = matches!(basis, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_));
    if !elementary || OffsetSurface::new(basis.clone(), along).equivalent_surface().is_some() {
        return Ok(None);
    }
    match basis {
        // (A torus moved through its tube's middle circle, or past its axis: no torus, but not
        // collapsing unless its tube vanishes.)
        Surface::Torus(_) if super::shape::moved_tube(&basis, along).is_some_and(|(minor, _)| minor.abs() > crate::tolerance::LINEAR) => Ok(None),
        Surface::Sphere(s) => Ok(Some(Focus::Centre(s.frame().origin()))),
        Surface::Cylinder(c) => Ok(Some(Focus::Axis(c.frame().origin(), Vector::from(c.frame().z())))),
        _ => Err(Error::InvalidInput("a thickening collapsing a face onto a circle or point other than a sphere's centre or a cylinder's axis".into())),
    }
}

/// Points along `edge` as it runs, `n` + 1 of them, ends included.
fn samples(model: &Model, edge: Shape, n: usize) -> Result<Vec<Point>> {
    let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("an edge without a curve".into()))?;
    let reversed = edge.orientation == Orientation::Reversed;
    Ok((0..=n)
        .map(|i| {
            let s = i as f64 / n as f64;
            curve.point(if reversed { last - (last - first) * s } else { first + (last - first) * s })
        })
        .collect())
}

/// `edge` run the other way when `flip`.
fn run(edge: Shape, flip: bool) -> Shape {
    if flip { edge.reversed() } else { edge }
}

/// The vertices `edge` runs from and to (forward).
fn ends(model: &Model, edge: Shape) -> (Shape, Shape) {
    let (a, b) = model.edge_vertices(edge);
    (Shape { orientation: Orientation::Forward, ..a }, Shape { orientation: Orientation::Forward, ..b })
}

/// The normal of the plane through `points`, if they lie on one (within `tolerance`), its
/// sense arbitrary.
fn plane_normal(points: &[Point], tolerance: f64) -> Option<Vector> {
    let origin = *points.first()?;
    let mut best = Vector::ZERO;
    for i in 1..points.len() {
        for j in i + 1..points.len() {
            let n = (points[i] - origin).cross(points[j] - origin);
            if n.length() > best.length() {
                best = n;
            }
        }
    }
    let extent = points.iter().map(|p| p.distance(origin)).fold(0.0, f64::max);
    if best.length() <= 1e-9 * extent * extent {
        return None;
    }
    let n = best.normalized();
    points.iter().all(|p| (*p - origin).dot(n).abs() <= tolerance).then_some(n)
}

/// The signed area normal of the closed polygon `points` (Newell's), pointing so that the
/// polygon runs counter-clockwise about it.
fn newell(points: &[Point]) -> Vector {
    let n = points.len();
    (0..n).fold(Vector::ZERO, |acc, i| acc + (points[i] - Point::ORIGIN).cross(points[(i + 1) % n] - Point::ORIGIN) * 0.5)
}

/// Where the segments `a`→`b` and `c`→`d` cross strictly inside both, if they do.
fn crossing(a: Point, b: Point, c: Point, d: Point, tolerance: f64) -> Option<Point> {
    let (u, v, w) = (b - a, d - c, a - c);
    let (uu, uv, vv, uw, vw) = (u.dot(u), u.dot(v), v.dot(v), u.dot(w), v.dot(w));
    let det = uu * vv - uv * uv;
    if det <= 1e-12 * uu * vv {
        return None;
    }
    let (s, t) = ((uv * vw - vv * uw) / det, (uu * vw - uv * uw) / det);
    let (p, q) = (a + u * s, c + v * t);
    let margin = |x: f64, length: f64| x * length > tolerance && (1.0 - x) * length > tolerance;
    (p.distance(q) <= tolerance && margin(s, uu.sqrt()) && margin(t, vv.sqrt())).then_some(p)
}

/// The walls a sheet's free edges sweep, and the rulings between them: the straight edges each
/// free vertex of the sheet sweeps to its image, by the two vertices, from the sheet's (each
/// run of edges along it, from the sheet's vertex, split where walls fold across it).
struct Walls {
    rulings: FastMap<(NodeId, NodeId), Vec<Shape>>,
    tolerance: f64,
}

/// A wall to be made: the edge, its image (or the focus point it collapses to), and the
/// vertices at both ends of each, the edge run from `from` to `to`, the image back from
/// `to`'s image to `from`'s.
struct Plan {
    edge: Shape,
    image: Option<Shape>,
    from: Shape,
    to: Shape,
    from_image: Shape,
    to_image: Shape,
    /// A planar wall's normal (its sense arbitrary), `None` for a wall on no one plane.
    plane: Option<Vector>,
    /// Where the wall's end rulings cross (a planar wall folding across a centre of
    /// curvature, its normals meeting there).
    fold: Option<Point>,
}

impl Walls {
    /// The ruling from the sheet's vertex `from` to the image vertex `to`, made if it is not
    /// yet (through `through`, a vertex where it is split, if any).
    fn ruling(&mut self, model: &mut Model, from: Shape, to: Shape, through: Option<Shape>) -> Result<Vec<Shape>> {
        if let Some(found) = self.rulings.get(&(from.node, to.node)) {
            return Ok(found.clone());
        }
        let pieces = match through {
            Some(c) => vec![crate::make::edge(model, from, c)?, crate::make::edge(model, c, to)?],
            None => vec![crate::make::edge(model, from, to)?],
        };
        self.rulings.insert((from.node, to.node), pieces.clone());
        Ok(pieces)
    }

    /// The ruling's edges run from the image vertex back to the sheet's.
    fn ruling_back(&mut self, model: &mut Model, from: Shape, to: Shape, through: Option<Shape>) -> Result<Vec<Shape>> {
        Ok(self.ruling(model, from, to, through)?.into_iter().rev().map(|e| e.reversed()).collect())
    }

    /// The walls between the free edges of the sheet (`free`) and their images on `moved`, the
    /// sheet moved `distance`, and the sectors closing them at rounds' ends: faces oriented as
    /// the sweep's boundary.
    fn to_offset(&mut self, model: &mut Model, free: &[FreeEdge], moved: Shape, distance: f64) -> Result<Vec<Shape>> {
        let images = free_edges(model, moved);
        let mut used = vec![false; images.len()];
        let mut plans = Vec::with_capacity(free.len());
        for e in free {
            // The image: the moved sheet's free edge through the point the edge's middle moves
            // to.
            let points = samples(model, e.edge, 2)?;
            let target = points[1] + normal_at(model, e.face, points[1])? * distance;
            let at = crate::make::vertex(model, target);
            let mut nearest: Option<(f64, usize)> = None;
            for (k, image) in images.iter().enumerate() {
                let d = crate::distance::distance(model, at, image.edge)?.value;
                if nearest.is_none_or(|(best, _)| d < best) {
                    nearest = Some((d, k));
                }
            }
            let Some((gap, k)) = nearest else { return Err(Error::Failed("the sheet's offset has no free edge where the sheet has".into())) };
            if gap > 100.0 * self.tolerance {
                return Err(Error::Failed("a free edge of the sheet whose image its offset lacks".into()));
            }
            used[k] = true;
            // Run as the sweep's boundary has them: the edge and its image each the other way
            // from the faces holding them there (the sheet faces against the sweep, the offset
            // with it).
            let edge = run(e.edge, distance < 0.0);
            let image = run(images[k].edge, distance > 0.0);
            let (from, to) = ends(model, edge);
            // (The image runs back: from the image of the edge's end to that of its start.)
            let (image_start, image_end) = ends(model, image);
            let (p_from, p_to) = (model.point(from), model.point(to));
            let (target_from, target_to) = (p_from + normal_at(model, e.face, p_from)? * distance, p_to + normal_at(model, e.face, p_to)? * distance);
            let back = target_to.distance(model.point(image_start)) + target_from.distance(model.point(image_end));
            let along = target_to.distance(model.point(image_end)) + target_from.distance(model.point(image_start));
            if from.node != to.node && along < back {
                return Err(Error::Failed("the sheet's offset runs against the sheet along a free edge".into()));
            }
            plans.push(self.plan(model, edge, Some(image), from, to, image_start, image_end)?);
        }
        let mut walls = self.build(model, &plans)?;
        // The moved sheet's free edges that are no image: the ends of rounds (in one edge or
        // more), each closed by a sector to the vertex of the sheet whose rulings reach both
        // its ends.
        let rest: Vec<(Shape, Shape, Shape)> = images
            .iter()
            .zip(&used)
            .filter(|(_, &u)| !u)
            .map(|(image, _)| {
                let end = run(image.edge, distance > 0.0);
                let (a, b) = ends(model, end);
                (end, a, b)
            })
            .collect();
        let reached = |rulings: &FastMap<(NodeId, NodeId), Vec<Shape>>, m: NodeId| rulings.keys().filter(|k| k.1 == m).map(|k| k.0).collect::<Vec<_>>();
        let mut taken = vec![false; rest.len()];
        for k in 0..rest.len() {
            if taken[k] || reached(&self.rulings, rest[k].1.node).is_empty() {
                continue;
            }
            // From a ruling's end along the round's end to the next.
            let mut chain = vec![k];
            taken[k] = true;
            while reached(&self.rulings, rest[*chain.last().expect("a chain")].2.node).is_empty() {
                let at = rest[*chain.last().expect("a chain")].2.node;
                let Some(next) = (0..rest.len()).find(|&j| !taken[j] && rest[j].1.node == at) else {
                    return Err(Error::Failed("a round's end at a free edge of the sheet that does not reach the sheet's vertex".into()));
                };
                taken[next] = true;
                chain.push(next);
            }
            let (a, b) = (rest[chain[0]].1, rest[*chain.last().expect("a chain")].2);
            let from_a = reached(&self.rulings, a.node);
            let Some(v) = reached(&self.rulings, b.node).into_iter().find(|v| from_a.contains(v)) else {
                return Err(Error::Failed("a round's end at a free edge of the sheet whose ends the rulings of no one vertex reach".into()));
            };
            let v = Shape::new(v);
            let mut loop_: Vec<Shape> = chain.iter().map(|&j| rest[j].0).collect();
            let mut points = Vec::new();
            for &edge in &loop_ {
                points.extend(samples(model, edge, 16)?);
            }
            points.push(model.point(v));
            if plane_normal(&points, self.tolerance).is_none() {
                return Err(Error::InvalidInput("a round's end at a free edge of the sheet that lies on no one plane with the edge's vertex".into()));
            }
            loop_.extend(self.ruling_back(model, v, b, None)?);
            loop_.extend(self.ruling(model, v, a, None)?);
            let wire = crate::make::wire(model, &loop_)?;
            walls.push(crate::make::planar_face(model, wire)?);
        }
        if let Some(k) = (0..rest.len()).find(|&k| !taken[k]) {
            let _ = k;
            return Err(Error::Failed("the sheet's offset has a free edge that is neither an image of the sheet's nor part of a round's end at one of its vertices".into()));
        }
        Ok(walls)
    }

    /// The walls between the free edges of the sheet (`free`) and the foci its faces (`faces`)
    /// collapse onto (`foci`, the faces' own), moved `distance`.
    fn to_foci(&mut self, model: &mut Model, free: &[FreeEdge], faces: &[Shape], foci: &[Focus], distance: f64) -> Result<Vec<Shape>> {
        let mut images: Vec<(Point, Shape)> = Vec::new();
        let tolerance = self.tolerance;
        let mut image_vertex = |model: &mut Model, p: Point| -> Shape {
            match images.iter().find(|(q, _)| q.distance(p) <= tolerance) {
                Some(&(_, v)) => v,
                None => {
                    let v = crate::make::vertex(model, p);
                    images.push((p, v));
                    v
                }
            }
        };
        let mut axes: FastMap<(NodeId, NodeId), Shape> = FastMap::default();
        let mut plans = Vec::with_capacity(free.len());
        for e in free {
            let k = faces.iter().position(|f| f.node == e.face.node).expect("a free edge's face is the sheet's");
            let edge = run(e.edge, distance < 0.0);
            let (from, to) = ends(model, edge);
            let (from_image, to_image) = (image_vertex(model, foci[k].image(model.point(from))), image_vertex(model, foci[k].image(model.point(to))));
            // (The image on an axis: the segment of it, one edge for both walls along it.)
            let image = if from_image.node == to_image.node {
                None
            } else {
                let key = if from_image.node < to_image.node { (from_image.node, to_image.node) } else { (to_image.node, from_image.node) };
                let segment = match axes.get(&key) {
                    Some(&s) => s,
                    None => {
                        let s = crate::make::edge(model, Shape::new(key.0), Shape::new(key.1))?;
                        axes.insert(key, s);
                        s
                    }
                };
                Some(if ends(model, segment).0.node == to_image.node { segment } else { segment.reversed() })
            };
            plans.push(self.plan(model, edge, image, from, to, to_image, from_image)?);
        }
        self.build(model, &plans)
    }

    /// What the wall of `edge` (run `from` → `to`) and its image (run `to_image` →
    /// `from_image`, or the point they both are) is: planar or not, folded or not.
    #[allow(clippy::too_many_arguments)]
    fn plan(&mut self, model: &mut Model, edge: Shape, image: Option<Shape>, from: Shape, to: Shape, to_image: Shape, from_image: Shape) -> Result<Plan> {
        let n = 16;
        let mut points = samples(model, edge, n)?;
        let image_points = match image {
            Some(i) => samples(model, i, n)?,
            None => vec![model.point(to_image)],
        };
        points.extend(image_points.iter().copied());
        let plane = plane_normal(&points, self.tolerance);
        // The end rulings crossing inside both: a fold, which a planar wall's rulings make
        // all through one point (an arc moved past its centre); elsewhere the wall is not
        // made.
        let (a, b, a2, b2) = (model.point(from), model.point(to), model.point(from_image), model.point(to_image));
        let fold = if from.node == to.node { None } else { crossing(a, a2, b, b2, self.tolerance) };
        if let Some(c) = fold {
            let edge_points = samples(model, edge, n)?;
            let mut image_back = image_points.clone();
            image_back.reverse();
            let through = plane.is_some()
                && edge_points.iter().zip(&image_back).all(|(&p, &q)| {
                    let d = q - p;
                    let length = d.length();
                    length > 0.0 && (c - p).cross(d).length() / length <= 1e-6 * length.max(1.0)
                });
            if !through {
                return Err(Error::InvalidInput("a wall at a free edge of the sheet folding other than through one point".into()));
            }
        }
        Ok(Plan { edge, image, from, to, from_image, to_image, plane, fold })
    }

    /// The walls `plans` describe, their rulings made first (split where walls fold).
    fn build(&mut self, model: &mut Model, plans: &[Plan]) -> Result<Vec<Shape>> {
        // (Where a ruling is split, every wall along it uses the pieces.)
        let mut splits: FastMap<(NodeId, NodeId), Shape> = FastMap::default();
        for plan in plans {
            if let Some(c) = plan.fold {
                let c = crate::make::vertex(model, c);
                splits.insert((plan.from.node, plan.from_image.node), c);
                splits.insert((plan.to.node, plan.to_image.node), c);
            }
        }
        let mut walls = Vec::new();
        for plan in plans {
            let through = |from: Shape, to: Shape| splits.get(&(from.node, to.node)).copied();
            let closed = plan.from.node == plan.to.node;
            if plan.plane.is_none() {
                walls.push(self.ruled(model, plan)?);
                continue;
            }
            if closed {
                walls.extend(self.planar_closed(model, plan)?);
                continue;
            }
            let to_ruling = self.ruling(model, plan.to, plan.to_image, through(plan.to, plan.to_image))?;
            let from_ruling = self.ruling_back(model, plan.from, plan.from_image, through(plan.from, plan.from_image))?;
            let loops: Vec<Vec<Shape>> = match plan.fold {
                // A figure eight through the fold: the edge's side of it, and the image's.
                Some(_) => vec![
                    vec![plan.edge, to_ruling[0], from_ruling[1]],
                    std::iter::once(to_ruling[1]).chain(plan.image).chain([from_ruling[0]]).collect(),
                ],
                None => vec![std::iter::once(plan.edge).chain(to_ruling).chain(plan.image).chain(from_ruling).collect()],
            };
            for edges in loops {
                if edges.len() < 2 {
                    continue;
                }
                let wire = crate::make::wire(model, &edges)?;
                walls.push(crate::make::planar_face(model, wire)?);
            }
        }
        Ok(walls)
    }

    /// The planar wall of a closed edge: the ring between it and its image, or where they wind
    /// round the same way (the wall folding through their common centre) or the image is a
    /// point, the disks within each.
    fn planar_closed(&mut self, model: &mut Model, plan: &Plan) -> Result<Vec<Shape>> {
        let normal = plan.plane.expect("a planar wall");
        let outer = crate::make::wire(model, &[plan.edge])?;
        let Some(image) = plan.image else {
            return Ok(vec![crate::make::planar_face(model, outer)?]);
        };
        let inner = crate::make::wire(model, &[image])?;
        let (edge_area, image_area) = (newell(&samples(model, plan.edge, 64)?).dot(normal), newell(&samples(model, image, 64)?).dot(normal));
        if (edge_area > 0.0) == (image_area > 0.0) {
            return Ok(vec![crate::make::planar_face(model, outer)?, crate::make::planar_face(model, inner)?]);
        }
        let (outside, hole) = if edge_area.abs() >= image_area.abs() { (outer, inner) } else { (inner, outer) };
        let face = crate::make::planar_face(model, outside)?;
        let surface = model.surface(face).expect("a planar face's plane").clone();
        Ok(vec![crate::make::face_bounded_on(model, surface, &[outside, hole])?])
    }

    /// The wall of a plan on no one plane: the ruled surface between the edge and its image,
    /// facing as the edge's run and the rulings turn.
    fn ruled(&mut self, model: &mut Model, plan: &Plan) -> Result<Shape> {
        let Some(image) = plan.image else {
            return Err(Error::InvalidInput("a free edge of the sheet sweeping a cone to the point it collapses onto".into()));
        };
        if plan.fold.is_some() {
            return Err(Error::InvalidInput("a wall at a free edge of the sheet folding on no one plane".into()));
        }
        let from = crate::make::wire(model, &[plan.edge])?;
        let to = crate::make::wire(model, &[image.reversed()])?;
        let loft = crate::loft::loft(model, &[from, to], &crate::loft::LoftOptions { ruled: true, ..Default::default() })?;
        let faces = model.unique_subshapes(loft.shape, NodeKind::Face);
        let [face] = faces[..] else { return Err(Error::Failed("a wall lofted into more than one face".into())) };
        // Facing out of the sweep: as the edge runs and then turns along the ruling from its
        // middle.
        let points = samples(model, plan.edge, 2)?;
        let image_points = samples(model, image, 2)?;
        let tangent = points[2] - points[0];
        let expected = if plan.from.node == plan.to.node {
            let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..plan.edge }).expect("a free edge's curve");
            let t = 0.5 * (first + last);
            let d = curve.d1(t).1 * if plan.edge.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
            d.cross(image_points[1] - points[1])
        } else {
            tangent.cross(image_points[1] - points[1])
        };
        let normal = normal_at(model, face, points[1])?;
        Ok(if normal.dot(expected) < 0.0 { face.reversed() } else { face })
    }
}

/// The solid of the cells `surface` (faces oriented as the boundary of a sweep, which may pass
/// through itself) winds round: the space it splits into cells, each cell's winding number
/// counted across the faces between it and its neighbours from the space outside them all.
fn swept_cells(model: &mut Model, surface: &[Shape], tolerance: f64) -> Result<Shape> {
    let enclosed = crate::boolean::enclosed_solids(model, surface)?;
    let cells = model.unique_subshapes(enclosed, NodeKind::Solid);
    if cells.is_empty() {
        return Err(Error::Failed("the sheet, its offset and the walls between them enclose nothing".into()));
    }
    // Each face of a cell, and the cells holding it.
    let mut holders: FastMap<NodeId, Vec<(usize, Shape)>> = FastMap::default();
    let mut faces_of: Vec<Vec<Shape>> = Vec::with_capacity(cells.len());
    for (i, &cell) in cells.iter().enumerate() {
        let faces: Vec<Shape> = model.explore(cell, NodeKind::Face).collect();
        for &f in &faces {
            holders.entry(f.node).or_default().push((i, f));
        }
        faces_of.push(faces);
    }
    // How much the winding number rises crossing each face into the cell holding it first
    // (against the face's normal there): the turns of the surface's faces it lies on.
    let boxes: Vec<crate::math::BoundingBox> = surface.iter().map(|&f| crate::properties::bounding_box(model, f).enlarged(10.0 * tolerance)).collect();
    let mut rise: FastMap<NodeId, i32> = FastMap::default();
    for (&node, held) in &holders {
        let face = held[0].1;
        let (p, normal) = inside(model, face)?;
        let mut turns = 0;
        for (k, &sigma) in surface.iter().enumerate() {
            if !boxes[k].contains(p) {
                continue;
            }
            let at = crate::make::vertex(model, p);
            if crate::distance::distance(model, at, sigma)?.value > 10.0 * tolerance {
                continue;
            }
            turns += if normal_at(model, sigma, p)?.dot(normal) > 0.0 { 1 } else { -1 };
        }
        rise.insert(node, turns);
    }
    // From outside (winding 0) across the faces held by one cell, then from cell to cell: the
    // winding number falls by the turns of a face crossed along its normal, out of a cell.
    let turns = |face: Shape| -> i32 {
        let first = holders[&face.node][0].1;
        if first.orientation == face.orientation { rise[&face.node] } else { -rise[&face.node] }
    };
    let mut winding: Vec<Option<i32>> = vec![None; cells.len()];
    let mut queue = std::collections::VecDeque::new();
    let reach = |winding: &mut Vec<Option<i32>>, queue: &mut std::collections::VecDeque<usize>, cell: usize, w: i32| -> Result<()> {
        match winding[cell] {
            None => {
                winding[cell] = Some(w);
                queue.push_back(cell);
                Ok(())
            }
            Some(known) if known != w => Err(Error::Failed("the sweep's boundary winds round a cell inconsistently".into())),
            _ => Ok(()),
        }
    };
    for (i, faces) in faces_of.iter().enumerate() {
        for &f in faces {
            if holders[&f.node].len() == 1 {
                reach(&mut winding, &mut queue, i, turns(f))?;
            }
        }
    }
    while let Some(i) = queue.pop_front() {
        let w = winding[i].expect("a cell reached");
        for &f in &faces_of[i] {
            for &(j, _) in &holders[&f.node] {
                if j != i {
                    reach(&mut winding, &mut queue, j, w - turns(f))?;
                }
            }
        }
    }
    let kept: Vec<usize> = (0..cells.len()).filter(|&i| winding[i].is_some_and(|w| w != 0)).collect();
    match kept[..] {
        [] => Err(Error::Failed("the sweep of the sheet covers none of the space its boundary encloses".into())),
        [one] => Ok(cells[one]),
        _ => {
            // The kept cells' faces but those between two of them.
            let mut faces = Vec::new();
            for &i in &kept {
                for &f in &faces_of[i] {
                    if holders[&f.node].iter().filter(|(j, _)| kept.contains(j)).count() == 1 {
                        faces.push(f);
                    }
                }
            }
            let merged = crate::boolean::enclosed_solids_unsplit(model, &faces)?;
            let solids = model.unique_subshapes(merged, NodeKind::Solid);
            match solids[..] {
                [one] => Ok(one),
                [] => Err(Error::Failed("the cells the sweep covers do not merge into a solid".into())),
                _ => Ok(merged),
            }
        }
    }
}

/// A point inside `face` (as it is used: away from its edges) and its unit normal there.
fn inside(model: &Model, face: Shape) -> Result<(Point, Vector)> {
    let loops = crate::heal::face_loops(model, Shape { orientation: Orientation::Forward, ..face }).ok_or_else(|| Error::Failed("a cell's face without pcurves".into()))?;
    let polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect();
    let outer = polygons.first().ok_or_else(|| Error::Failed("a cell's face without loops".into()))?;
    let size = outer.iter().fold(0.0f64, |m, p| m.max((*p - outer[0]).length()));
    let near_boundary = |q: Point2| {
        polygons.iter().flat_map(|p| (0..p.len()).map(move |i| (p[i], p[(i + 1) % p.len()]))).map(|(a, b)| {
            let ab = b - a;
            let t = ((q - a).dot(ab) / ab.dot(ab).max(1e-300)).clamp(0.0, 1.0);
            (a + ab * t - q).length()
        }).fold(f64::INFINITY, f64::min)
    };
    // Just inside the middle of each side of the outer loop, deeper first: the deepest found.
    let mut best: Option<(f64, Point2)> = None;
    for depth in [0.25, 0.1, 0.03, 0.01, 1e-3] {
        for i in 0..outer.len() {
            let (a, b) = (outer[i], outer[(i + 1) % outer.len()]);
            let along = b - a;
            if along.length() <= 1e-12 * size {
                continue;
            }
            let k = depth * size / along.length();
            let q = Point2::new(0.5 * (a.x + b.x) - along.y * k, 0.5 * (a.y + b.y) + along.x * k);
            if super::regularize::winding(&polygons, q) == 0 {
                continue;
            }
            let clear = near_boundary(q);
            if best.is_none_or(|(c, _)| clear > c) {
                best = Some((clear, q));
            }
        }
        if best.is_some_and(|(c, _)| c >= 0.5 * depth * size) {
            break;
        }
    }
    let (_, q) = best.ok_or_else(|| Error::Failed("no point found inside a cell's face".into()))?;
    let surface = model.surface(face).ok_or_else(|| Error::Failed("a cell's face without a surface".into()))?.transformed(&face.placement);
    let (p, du, dv) = surface.d1(q.x, q.y);
    let n = du.cross(dv);
    if n.length() <= 1e-300 {
        return Err(Error::Failed("a cell's face without a normal inside it".into()));
    }
    Ok((p, n.normalized() * if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::{Circle, Curve};
    use crate::geometry::surface::{Cone, Sphere};
    use crate::math::{Direction, Frame};
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    fn point(x: f64, y: f64, z: f64) -> Point {
        Point::new(x, y, z)
    }

    /// `solid` closed and sound, its volume and area those given (to 1e-9, relative).
    fn assert_measures(model: &Model, solid: Shape, volume: f64, area: f64) {
        let report = crate::validity::check(model, solid);
        assert!(report.is_valid(), "{report:?}");
        assert!(crate::boolean::closed_shells(model, solid));
        let v = crate::properties::volume(model, solid).mass;
        let a = crate::properties::area(model, solid).mass;
        assert!((v - volume).abs() <= 1e-9 * volume, "volume {v}, not {volume}");
        assert!((a - area).abs() <= 1e-9 * area, "area {a}, not {area}");
    }

    /// The polyline through `points` swept by `along`: a sheet of flat faces.
    fn swept_polyline(model: &mut Model, points: &[Point], along: Vector) -> Shape {
        let edges: Vec<Shape> = points.windows(2).map(|w| crate::make::edge_between(model, w[0], w[1]).unwrap()).collect();
        let wire = crate::make::wire(model, &edges).unwrap();
        crate::sweep::extrude(model, wire, along).unwrap().shape
    }

    /// A creased sheet (the floor from (0, 0, 0) to (10, 8, 0), the wall from (0, 0, 0) to
    /// (0, 8, 10)) thickened 1 to the outside of the crease, joined by arcs: the floor's and
    /// the wall's slabs and the quarter cylinder about the crease between them, closed at each
    /// end by a quarter disk from the crease's end to the round's end.
    #[test]
    fn a_creased_sheet_thickened_outside_its_crease_is_rounded_there() {
        let mut model = Model::new();
        let sheet = swept_polyline(&mut model, &[point(10.0, 0.0, 0.0), point(0.0, 0.0, 0.0), point(0.0, 0.0, 10.0)], Vector::new(0.0, 8.0, 0.0));
        let solid = thicken(&mut model, sheet, -1.0, Join::Arc).unwrap();
        let (l, t, w) = (10.0, 1.0, 8.0);
        let section = 2.0 * l * t + PI * t * t / 4.0;
        assert_measures(&model, solid, w * section, 2.0 * section + w * (4.0 * l + 2.0 * t + PI * t / 2.0));
    }

    /// An open box's five faces thickened 1 outwards, joined by arcs: the slabs, a quarter
    /// cylinder along each of the eight edges, an eighth ball at each of the floor's corners,
    /// its rim flat (quarter disks at the upright edges' tops).
    #[test]
    fn an_open_box_thickened_outwards_is_rounded_at_its_edges_and_corners() {
        let mut model = Model::new();
        let cube = crate::primitives::Cuboid::from_corners(point(0.0, 0.0, 0.0), point(10.0, 10.0, 10.0)).build(&mut model).unwrap();
        let shell = model.new_empty(NodeKind::Shell);
        for face in model.unique_subshapes(cube, NodeKind::Face) {
            if crate::properties::area(&model, face).center_of_mass.z < 9.9 {
                model.add(shell, face).unwrap();
            }
        }
        let solid = thicken(&mut model, shell, 1.0, Join::Arc).unwrap();
        let (a, t) = (10.0, 1.0);
        let volume = 5.0 * a * a * t + 8.0 * PI * t * t * a / 4.0 + 4.0 * PI * t * t * t / 6.0;
        let area = 10.0 * a * a + 8.0 * PI * t * a / 2.0 + 4.0 * PI * t * t / 2.0 + 4.0 * a * t + PI * t * t;
        assert_measures(&model, solid, volume, area);
    }

    /// A sheet bent round a quarter cylinder of radius 2 thickened 3 inside the bend: the
    /// bend's offset turns inside out within the solid, the floor's and the wall's offsets
    /// meet in a concave edge, and the bend fills the quarter disk behind it.
    #[test]
    fn a_bent_sheet_thickened_past_its_bend() {
        let mut model = Model::new();
        let floor = crate::make::edge_between(&mut model, point(-10.0, 0.0, 0.0), point(0.0, 0.0, 0.0)).unwrap();
        let frame = Frame::new(point(0.0, 0.0, 2.0), Direction::try_from_vector(Vector::new(0.0, -1.0, 0.0)).unwrap(), Direction::try_from_vector(Vector::new(0.0, 0.0, -1.0)).unwrap());
        let bend = crate::make::edge_from_curve_range(&mut model, Curve::Circle(Circle::new(frame, 2.0)), 0.0, FRAC_PI_2).unwrap();
        let wall = crate::make::edge_between(&mut model, point(2.0, 0.0, 2.0), point(2.0, 0.0, 12.0)).unwrap();
        let profile = crate::make::wire(&mut model, &[floor, bend, wall]).unwrap();
        let sheet = crate::sweep::extrude(&mut model, profile, Vector::new(0.0, 8.0, 0.0)).unwrap().shape;
        let solid = thicken(&mut model, sheet, -3.0, Join::Arc).unwrap();
        let (l, r, t, w) = (10.0, 2.0, 3.0, 8.0);
        let section = 2.0 * l * t - (t - r) * (t - r) + PI * r * r / 4.0;
        assert_measures(&model, solid, w * section, 2.0 * section + w * (4.0 * l + 2.0 * r + PI * r / 2.0));
    }

    /// A hemisphere of radius 5 thickened 5 towards its centre: its offset collapses to the
    /// centre, and the solid is the half ball, its equator sweeping the disk.
    #[test]
    fn a_hemisphere_thickened_to_its_centre_is_the_half_ball() {
        let mut model = Model::new();
        let sheet = crate::make::face_on_domain(&mut model, Surface::Sphere(Sphere::new(Frame::WORLD, 5.0)), 0.0, TAU, 0.0, FRAC_PI_2).unwrap();
        let solid = thicken(&mut model, sheet, -5.0, Join::Arc).unwrap();
        assert_measures(&model, solid, 2.0 * PI * 125.0 / 3.0, 3.0 * PI * 25.0);
    }

    /// A hemisphere moved 1 inwards (an offset of the face) and that thickened 2 outwards: the
    /// dome between the spheres of radius 4 and 6.
    #[test]
    fn an_offset_hemisphere_thickened_outwards() {
        let mut model = Model::new();
        let sheet = crate::make::face_on_domain(&mut model, Surface::Sphere(Sphere::new(Frame::WORLD, 5.0)), 0.0, TAU, 0.0, FRAC_PI_2).unwrap();
        let inner = super::super::shape(&mut model, sheet, -1.0, Join::Arc).unwrap();
        let solid = thicken(&mut model, inner, 2.0, Join::Arc).unwrap();
        assert_measures(&model, solid, 2.0 * PI * (216.0 - 64.0) / 3.0, 2.0 * PI * (36.0 + 16.0) + PI * (36.0 - 16.0));
    }

    /// A cone frustum's side (radius 6 at z = 0 to 3 at z = 4) thickened 1 along its normal:
    /// the section a 5 by 1 rectangle, its rims cones; by Pappus, 2π x̄ · 5, x̄ = 4.9.
    #[test]
    fn a_frustum_thickened_along_its_slanting_normal() {
        let mut model = Model::new();
        let cone = Surface::Cone(Cone::new(Frame::WORLD, -(3.0f64).atan2(4.0), 6.0));
        let sheet = crate::make::face_on_domain(&mut model, cone, 0.0, TAU, 0.0, 5.0).unwrap();
        let solid = thicken(&mut model, sheet, 1.0, Join::Arc).unwrap();
        assert_measures(&model, solid, TAU * 4.9 * 5.0, TAU * (4.5 * 5.0 + 5.3 * 5.0 + 6.4 + 3.4));
    }

    /// A sphere's face (closed) thickened 1 inwards: the shell between radii 4 and 5; the ball
    /// the offset bounds inside, which the sweep moves away from, left out.
    #[test]
    fn a_closed_sheet_thickened_leaves_out_the_space_inside_its_offset() {
        let mut model = Model::new();
        let ball = crate::primitives::Sphere::new(5.0).build(&mut model).unwrap();
        let face = model.unique_subshapes(ball, NodeKind::Face)[0];
        let solid = thicken(&mut model, face, -1.0, Join::Arc).unwrap();
        assert_eq!(model.unique_subshapes(solid, NodeKind::Shell).len(), 2);
        assert_measures(&model, solid, 4.0 * PI * (125.0 - 64.0) / 3.0, 4.0 * PI * (25.0 + 16.0));
    }

    /// No distance, or no faces: refused.
    #[test]
    fn a_thickening_by_nothing_is_refused() {
        let mut model = Model::new();
        let sheet = swept_polyline(&mut model, &[point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0)], Vector::new(0.0, 1.0, 0.0));
        assert!(thicken(&mut model, sheet, 0.0, Join::Arc).is_err());
        let empty = model.new_empty(NodeKind::Shell);
        assert!(thicken(&mut model, empty, 1.0, Join::Arc).is_err());
    }

    /// A cup of revolution about the z axis (the floor out to radius `c`, a quarter torus of
    /// radius `r` bending up to the wall of radius `c + r`, the wall up to `h`), its faces
    /// tangent.
    fn cup(model: &mut Model, c: f64, r: f64, h: f64) -> Shape {
        let floor = crate::make::edge_between(model, point(0.0, 0.0, 0.0), point(c, 0.0, 0.0)).unwrap();
        let frame = Frame::new(point(c, 0.0, r), Direction::try_from_vector(Vector::new(0.0, -1.0, 0.0)).unwrap(), Direction::try_from_vector(Vector::new(0.0, 0.0, -1.0)).unwrap());
        let bend = crate::make::edge_from_curve_range(model, Curve::Circle(Circle::new(frame, r)), 0.0, FRAC_PI_2).unwrap();
        let wall = crate::make::edge_between(model, point(c + r, 0.0, r), point(c + r, 0.0, h)).unwrap();
        let profile = crate::make::wire(model, &[floor, bend, wall]).unwrap();
        crate::sweep::revolve(model, profile, crate::math::Axis::new(Point::ORIGIN, Direction::Z), TAU).unwrap().shape
    }

    /// Cups thickened inside by more than their bends' radii (by 3 past a bend of 2, 4 past 2, 3
    /// past 1): the bend's offset turns inside out through its tube's middle circle (a torus
    /// moved so is no torus with its parameters, but the offset surface), the moved floor and
    /// wall meet in a concave edge, and the bend fills the quarter disk behind it.
    #[test]
    fn cups_thickened_inside_past_their_bends() {
        for (c, r, h, t) in [(4.0, 2.0, 10.0, 3.0), (3.0, 2.0, 10.0, 4.0), (5.0, 1.0, 8.0, 3.0)] {
            let mut model = Model::new();
            let sheet = cup(&mut model, c, r, h);
            let solid = thicken(&mut model, sheet, -t, Join::Arc).unwrap();
            let (big, w) = (c + r, c + r - t);
            let volume = PI * c * c * t + PI * (big * big - w * w) * (h - r) - PI * (c * c - w * w) * (t - r) + PI * PI * r * r / 2.0 * (c + 4.0 * r / (3.0 * PI));
            let area = PI * c * c + PI * PI * r * (c + 2.0 * r / PI) + 2.0 * PI * big * (h - r) + PI * w * w + 2.0 * PI * w * (h - t) + PI * (big * big - w * w);
            assert_measures(&model, solid, volume, area);
        }
    }

    /// A whole torus thickened, closed with no free edges: the solid between it and its
    /// offset, a cavity inside for an inward thickness.
    #[test]
    fn a_whole_torus_thickened() {
        let (big, small, t) = (3.0, 1.8, 1.3);
        let mut model = Model::new();
        let torus = crate::primitives::Torus::new(big, small).build(&mut model).unwrap();
        let face = model.unique_subshapes(torus, NodeKind::Face)[0];
        let solid = thicken(&mut model, face, -t, Join::Arc).unwrap();
        let inner = small - t;
        assert_measures(&model, solid, 2.0 * PI * PI * big * (small * small - inner * inner), 4.0 * PI * PI * big * (small + inner));
    }
}
