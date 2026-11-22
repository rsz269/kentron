//! Drawing a shape as a viewer shows it: each face meshed on its own with the surface's
//! outward normal at each node, and each edge as the polyline the faces' meshes share, with
//! the faces it bounds and how sharply they meet along it.
//!
//! As OCCT's shaded presentation does (`StdPrs_ShadedShape` on the faces' `Poly_Triangulation`,
//! their normals from `BRepLib_ToolTriangulatedShape::ComputeNormals`, which takes each node's
//! from `GeomLib::NormEstim` at its (u, v)), and its wireframe (`StdPrs_WFShape`, drawing each
//! edge by its polygon on the triangulation, seams apart). Where `NormEstim` finds the surface
//! has no normal at a node (a cone's apex), OCCT gives the node the mean of its triangles'
//! normals; here each triangle gets a node of its own there, with the normal the surface has
//! halfway across that triangle, so that the shading follows the surface on each side of the
//! point.

use crate::hashing::FastMap as HashMap;

use crate::error::Result;
use crate::geometry::props::estimated_normal;
use crate::geometry::surface::Surface;
use crate::math::{Direction, Point, Point2, Transform, Transformable, Vector};
use crate::tolerance::LINEAR;
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeKind, Orientation, Shape};

use super::{mesh_with, MeshOptions};

/// A shape as a viewer draws it (see [`display`]).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Display {
    /// Each face meshed on its own, in the order the shape lists its faces.
    pub faces: Vec<DisplayFace>,
    /// Each edge between faces (not of no length, not a seam).
    pub edges: Vec<DisplayEdge>,
    /// Each seam, where a face meets itself, as a polyline.
    pub seams: Vec<Vec<Point>>,
}

/// One face of a [`Display`]: its own nodes, with the face's outward unit normal at each,
/// and its triangles counter-clockwise seen from outside.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DisplayFace {
    pub positions: Vec<Point>,
    pub normals: Vec<Direction>,
    pub triangles: Vec<[u32; 3]>,
}

/// One edge of a [`Display`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DisplayEdge {
    /// The edge's polyline: the nodes the meshes of its faces share along it.
    pub points: Vec<Point>,
    /// Each face the edge bounds, as an index into [`Display::faces`], once for each use.
    pub faces: Vec<u32>,
    /// For an edge between two distinct faces, the least and greatest angle between their
    /// outward normals over `points` (where both faces have a normal), in radians: positive
    /// where the edge is convex (the solid filling less than half the space about it),
    /// negative where it is concave, 0 where the faces meet smoothly.
    pub turn: Option<[f64; 2]>,
    /// The two faces lie on one surface, so that only a split made the edge.
    pub split: bool,
}

/// A face as drawn, with what the edges need of it.
struct Drawn {
    face: Shape,
    id: FaceId,
    /// The surface in the shape's coordinates.
    surface: Surface,
    /// The face's outward normal at each node of its mesh, in the shape's coordinates;
    /// `None` where the surface has none.
    normals: Vec<Option<Direction>>,
    /// The mesh's nodes, in the face's own coordinates.
    local_nodes: Vec<Point>,
}

/// Draws `shape` as a viewer shows it, meshed within `tolerance` (as [`mesh`](super::mesh)
/// meshes it, and stored there): each face meshed on its own, with the surface's outward
/// normal at each node; each edge between faces as the polyline their meshes share, with
/// its faces and their turn along it, and whether the faces lie on one surface; the seams
/// apart. Degenerate edges (a pole, a cone's apex) are not drawn.
#[doc(alias = "StdPrs_ShadedShape")]
#[doc(alias = "StdPrs_WFShape")]
#[doc(alias = "BRepLib_ToolTriangulatedShape::ComputeNormals")]
pub fn display(model: &mut Model, shape: Shape, tolerance: f64) -> Result<Display> {
    mesh_with(model, shape, MeshOptions { deflection: tolerance, ..MeshOptions::default() })?;
    let model: &Model = model;
    let mut out = Display::default();
    let mut drawn: Vec<Drawn> = Vec::new();
    for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
        let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
        let (Some(mesh), Some(surface)) = (model.face_mesh(local), model.surface(local)) else { continue };
        let id = FaceId::try_from(face.node).expect("a face");
        let reversed = face.orientation == Orientation::Reversed;
        let placement = face.placement;
        let nodes = mesh.nodes();
        // Each node's (u, v): the mesh's, or the surface's own for a mesh made otherwise.
        let uv: Vec<Point2> = match mesh.uv_nodes() {
            Some(uv) => uv.to_vec(),
            None => nodes
                .iter()
                .map(|&p| {
                    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, surface, None);
                    Point2::new(u, v)
                })
                .collect(),
        };
        // The outward normal, in the face's own coordinates, as `NormEstim` finds it.
        let outward = |u: f64, v: f64| estimated_normal(surface, u, v, LINEAR).map(|n| if reversed { -n } else { n });
        let local_normals: Vec<Option<Direction>> = uv.iter().map(|p| outward(p.x, p.y)).collect();
        // Counter-clockwise seen from outside: the mesh's triangles run counter-clockwise
        // about the surface's normal, and a mirror turns them round again.
        let turned = reversed != placement.is_negative();
        let mut positions: Vec<Point> = Vec::new();
        let mut normals: Vec<Direction> = Vec::new();
        let mut vertex_of: Vec<Option<u32>> = vec![None; nodes.len()];
        let mut triangles = Vec::with_capacity(mesh.triangles().len());
        for t in mesh.triangles() {
            let mut corners = [0u32; 3];
            for k in 0..3 {
                let i = t[k] as usize;
                if let Some(n) = local_normals[i] {
                    corners[k] = *vertex_of[i].get_or_insert_with(|| {
                        positions.push(nodes[i].transformed(&placement));
                        normals.push(n.transformed(&placement));
                        (positions.len() - 1) as u32
                    });
                    continue;
                }
                // No one normal here: a node of this triangle's own, with the normal the
                // surface has halfway out from it across the triangle.
                let (b, c) = (t[(k + 1) % 3] as usize, t[(k + 2) % 3] as usize);
                let opposite = nodes[b] + (nodes[c] - nodes[b]) * 0.5;
                let halfway = nodes[i] + (opposite - nodes[i]) * 0.5;
                let estimate = (0.5 * (uv[b].x + uv[c].x), 0.5 * (uv[b].y + uv[c].y));
                let (u, v, _) = crate::geometry::extrema::nearest_on_surface(halfway, surface, Some(estimate));
                let n = outward(u, v).or_else(|| {
                    // (Else the triangle's own, as OCCT's flat normals.)
                    let [p, q, r] = t.map(|j| nodes[j as usize]);
                    let n = (q - p).cross(r - p);
                    Direction::try_from_vector(if reversed { -n } else { n })
                });
                let Some(n) = n else { continue };
                positions.push(nodes[i].transformed(&placement));
                normals.push(n.transformed(&placement));
                corners[k] = (positions.len() - 1) as u32;
            }
            triangles.push(if turned { [corners[0], corners[2], corners[1]] } else { corners });
        }
        out.faces.push(DisplayFace { positions, normals, triangles });
        drawn.push(Drawn {
            face,
            id,
            surface: surface.transformed(&placement),
            normals: local_normals.iter().map(|n| n.map(|n| n.transformed(&placement))).collect(),
            local_nodes: nodes.to_vec(),
        });
    }

    // Each edge's uses by the drawn faces, edges in order of first use (an edge placed
    // twice is two edges).
    let mut uses: Vec<(EdgeId, Transform, Vec<(usize, CoEdgeId)>)> = Vec::new();
    let mut by_edge: HashMap<EdgeId, Vec<usize>> = HashMap::default();
    for (f, d) in drawn.iter().enumerate() {
        for coedge in model.face_coedges(d.id) {
            let edge = model[coedge].edge;
            let entries = by_edge.entry(edge).or_default();
            let found = entries.iter().copied().find(|&j| uses[j].1 == d.face.placement);
            let j = found.unwrap_or_else(|| {
                uses.push((edge, d.face.placement, Vec::new()));
                entries.push(uses.len() - 1);
                uses.len() - 1
            });
            uses[j].2.push((f, coedge));
        }
    }
    for (edge, placement, used) in &uses {
        let edge_shape = Shape::from(*edge);
        if model.is_degenerate(edge_shape) {
            continue;
        }
        let Some(polyline) = model.edge_polyline(edge_shape) else { continue };
        let local = polyline.nodes();
        let points: Vec<Point> = local.iter().map(|p| p.transformed(placement)).collect();
        // Of no length: a pole, an apex.
        if points.len() < 2 || points.iter().all(|p| p.distance(points[0]) <= LINEAR) {
            continue;
        }
        let faces: Vec<u32> = used.iter().map(|&(f, _)| f as u32).collect();
        if faces.len() >= 2 && faces.iter().all(|&f| f == faces[0]) {
            out.seams.push(points);
            continue;
        }
        let distinct = used.len() == 2 && used[0].0 != used[1].0;
        let (mut turn, mut split) = (None, false);
        if distinct {
            let sides = used.iter().map(|&(f, c)| along(model, &drawn[f], c, local)).collect::<Option<Vec<_>>>();
            if let Some(sides) = sides {
                let (one, other) = (&sides[0], &sides[1]);
                let mut range: Option<[f64; 2]> = None;
                let mut agree = true;
                let n = points.len();
                for i in 0..n {
                    let (Some(n1), Some(n2)) = (one.normals[i], other.normals[i]) else { continue };
                    // The edge's direction where it runs with the first face on its left,
                    // seen from outside.
                    let chord = points[(i + 1).min(n - 1)] - points[i.saturating_sub(1)];
                    let along = chord * one.sense;
                    let across = n1.cross(n2);
                    let angle = across.length().atan2(n1.dot(n2));
                    let angle = if across.dot(along) < 0.0 { -angle } else { angle };
                    agree &= n1.dot(n2) > 0.0;
                    range = Some(range.map_or([angle, angle], |[lo, hi]| [lo.min(angle), hi.max(angle)]));
                }
                turn = range;
                let (a, b) = (&drawn[used[0].0], &drawn[used[1].0]);
                let same_slot = model.get(a.id).map(|f| f.surface) == model.get(b.id).map(|f| f.surface) && a.face.placement == b.face.placement;
                split = range.is_some() && agree && (same_slot || crate::unify::same_domain(&a.surface, &b.surface));
            }
        }
        out.edges.push(DisplayEdge { points, faces, turn, split });
    }
    Ok(out)
}

/// A face along one of its edges: its outward normal at each node of the edge's polyline,
/// and the sense (±1) in which the polyline runs with the face on its left seen from outside.
struct Along {
    normals: Vec<Option<Direction>>,
    sense: f64,
}

/// `d` along the edge of `coedge`, whose polyline has `local` nodes: from the coedge's chain
/// of nodes of the face's mesh. `None` where the mesh has no such chain.
fn along(model: &Model, d: &Drawn, coedge: CoEdgeId, local: &[Point]) -> Option<Along> {
    let chain = model.coedge_mesh_polyline(coedge)?.nodes();
    let n = local.len();
    if chain.len() != n || chain.iter().any(|&i| i as usize >= d.local_nodes.len()) {
        return None;
    }
    // The chain runs the way the face's boundary does about the surface's normal, the face
    // to its left: along the polyline, or against it.
    let node = |j: usize| d.local_nodes[chain[j] as usize];
    let forward = node(1) == local[1] && node(0) == local[0];
    let backward = node(n - 2) == local[1] && node(n - 1) == local[0];
    if forward == backward {
        return None;
    }
    let normals = (0..n).map(|i| d.normals[chain[if forward { i } else { n - 1 - i }] as usize]).collect();
    // About the outward normal: the other way round for a reversed face, and a mirror turns
    // it round again.
    let flips = [!forward, d.face.orientation == Orientation::Reversed, d.face.placement.is_negative()];
    let sense = if flips.iter().filter(|&&f| f).count() % 2 == 1 { -1.0 } else { 1.0 };
    Some(Along { normals, sense })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Frame;
    use crate::primitives::{Cone, Cuboid, Cylinder, Sphere};
    use std::f64::consts::{FRAC_PI_2, PI};

    const TOLERANCE: f64 = 0.01;

    /// Whether the faces, their nodes joined where they coincide, close: every side of a
    /// triangle used once each way round. And the volume they bound.
    fn closed(d: &Display) -> (bool, f64) {
        let mut ids: HashMap<[u64; 3], u32> = HashMap::default();
        let mut points: Vec<Point> = Vec::new();
        let mut triangles: Vec<[u32; 3]> = Vec::new();
        for f in &d.faces {
            let map: Vec<u32> = f
                .positions
                .iter()
                .map(|p| {
                    *ids.entry([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]).or_insert_with(|| {
                        points.push(*p);
                        (points.len() - 1) as u32
                    })
                })
                .collect();
            triangles.extend(f.triangles.iter().map(|t| t.map(|i| map[i as usize])).filter(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]));
        }
        let mut count: HashMap<(u32, u32), i32> = HashMap::default();
        let mut volume = 0.0;
        for t in &triangles {
            for k in 0..3 {
                *count.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
            }
            let [a, b, c] = t.map(|i| points[i as usize] - Point::ORIGIN);
            volume += a.dot(b.cross(c)) / 6.0;
        }
        let shut = count.iter().all(|(&(a, b), &n)| n == 1 && count.get(&(b, a)) == Some(&1));
        (shut, volume)
    }

    /// Every vertex lies on its face's surface with the surface's outward normal there (`at`
    /// gives both, `None` where the vertex is to be judged by its triangles), and every
    /// triangle faces the way its vertices' normals do.
    fn check_faces(d: &Display, at: impl Fn(Point) -> Option<(f64, Direction)>) {
        for f in &d.faces {
            assert_eq!(f.positions.len(), f.normals.len());
            for (p, n) in f.positions.iter().zip(&f.normals) {
                assert!((Vector::from(*n).length() - 1.0).abs() < 1e-12);
                if let Some((off, want)) = at(*p) {
                    assert!(off < 1e-9, "{p:?} off the surface by {off}");
                    assert!(n.dot(want) > 1.0 - 1e-9, "normal {n:?} at {p:?}, not {want:?}");
                }
            }
            for t in &f.triangles {
                let [a, b, c] = t.map(|i| f.positions[i as usize]);
                let n = (b - a).cross(c - a);
                if n.length() <= 1e-12 {
                    continue;
                }
                for &i in t {
                    assert!(n.dot(Vector::from(f.normals[i as usize])) > 0.0, "a triangle facing inwards: {:?} {:?} {:?} {n:?} {:?}", a, b, c, f.normals[i as usize]);
                }
            }
        }
    }

    fn turns(d: &Display) -> Vec<[f64; 2]> {
        d.edges.iter().map(|e| e.turn.expect("a turn")).collect()
    }

    #[test]
    fn a_box_has_six_faces_and_twelve_square_edges() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
        let d = display(&mut model, b, TOLERANCE).unwrap();
        assert_eq!((d.faces.len(), d.edges.len(), d.seams.len()), (6, 12, 0));
        let (shut, volume) = closed(&d);
        assert!(shut && (volume - 6000.0).abs() < 1e-9, "{shut} {volume}");
        check_faces(&d, |_| None);
        // Each face's normal is its side's, out of the box.
        let middle = Point::new(5.0, 10.0, 15.0);
        for f in &d.faces {
            let n = f.normals[0];
            assert!(f.normals.iter().all(|&m| m == n));
            let level = |p: Point| (p - middle).dot(Vector::from(n));
            assert!(level(f.positions[0]) > 0.0);
            assert!(f.positions.iter().all(|&p| (level(p) - level(f.positions[0])).abs() < 1e-12));
        }
        for e in &d.edges {
            assert_eq!(e.points.len(), 2);
            assert_eq!(e.faces.len(), 2);
            assert_ne!(e.faces[0], e.faces[1]);
            assert!(!e.split);
            let [lo, hi] = e.turn.unwrap();
            assert!((lo - FRAC_PI_2).abs() < 1e-12 && (hi - FRAC_PI_2).abs() < 1e-12, "{lo} {hi}");
        }
    }

    #[test]
    fn a_cylinder_lists_its_seam_apart_and_its_rims_turn_square() {
        let mut model = Model::new();
        let c = Cylinder::new(5.0, 8.0).build(&mut model).unwrap();
        let d = display(&mut model, c, TOLERANCE).unwrap();
        assert_eq!((d.faces.len(), d.edges.len(), d.seams.len()), (3, 2, 1));
        let (shut, volume) = closed(&d);
        let exact = PI * 25.0 * 8.0;
        assert!(shut && (volume - exact).abs() <= 2.0 * TOLERANCE * (2.0 * PI * 25.0 + 2.0 * PI * 5.0 * 8.0), "{shut} {volume}");
        check_faces(&d, |p| {
            let r = p.x.hypot(p.y);
            if p.z.abs() < 1e-12 && r < 5.0 - 1e-9 {
                return Some((0.0, Direction::NEG_Z));
            }
            if (p.z - 8.0).abs() < 1e-12 && r < 5.0 - 1e-9 {
                return Some((0.0, Direction::Z));
            }
            // On the side (and its rims, where either normal is the face's).
            let side = Direction::try_from_vector(Vector::new(p.x, p.y, 0.0)).unwrap();
            ((r - 5.0).abs() < 1e-9 && p.z > 1e-12 && p.z < 8.0 - 1e-12).then_some((0.0, side))
        });
        // The seam is a ruling from rim to rim.
        let seam = &d.seams[0];
        assert!(seam.iter().all(|p| (p.x.hypot(p.y) - 5.0).abs() < 1e-9));
        for e in &d.edges {
            // Within the tolerance of the circle, every vertex on it.
            assert!(e.points.iter().all(|p| (p.x.hypot(p.y) - 5.0).abs() < 1e-9));
            for w in e.points.windows(2) {
                let m = w[0] + (w[1] - w[0]) * 0.5;
                assert!(5.0 - m.x.hypot(m.y) <= TOLERANCE);
            }
            let [lo, hi] = e.turn.unwrap();
            assert!((lo - FRAC_PI_2).abs() < 1e-9 && (hi - FRAC_PI_2).abs() < 1e-9, "{lo} {hi}");
            assert!(!e.split);
        }
    }

    #[test]
    fn a_hole_through_a_box_has_square_rims_and_a_blind_ones_floor_is_concave() {
        let mut model = Model::new();
        let block = Cuboid::new(20.0, 20.0, 10.0).build(&mut model).unwrap();
        let drill = Cylinder::new(4.0, 20.0).with_frame(Frame::from_z(Point::new(10.0, 10.0, -5.0), Direction::Z)).build(&mut model).unwrap();
        let drilled = crate::boolean::difference(&mut model, block, drill).unwrap();
        let d = display(&mut model, drilled, TOLERANCE).unwrap();
        let (shut, volume) = closed(&d);
        let exact = 4000.0 - PI * 16.0 * 10.0;
        let area = 1600.0 - 2.0 * PI * 16.0 + 2.0 * PI * 4.0 * 10.0;
        assert!(shut && (volume - exact).abs() <= 2.0 * TOLERANCE * area, "{shut} {volume}");
        assert_eq!((d.edges.len(), d.seams.len()), (14, 1));
        // Through the box, the hole's rims are convex, square.
        for t in turns(&d) {
            assert!((t[0] - FRAC_PI_2).abs() < 1e-9 && (t[1] - FRAC_PI_2).abs() < 1e-9, "{t:?}");
        }

        // A blind hole: its floor meets its wall in a concave circle.
        let block = Cuboid::new(20.0, 20.0, 10.0).build(&mut model).unwrap();
        let drill = Cylinder::new(4.0, 20.0).with_frame(Frame::from_z(Point::new(10.0, 10.0, 4.0), Direction::Z)).build(&mut model).unwrap();
        let pocketed = crate::boolean::difference(&mut model, block, drill).unwrap();
        let d = display(&mut model, pocketed, TOLERANCE).unwrap();
        let (shut, _) = closed(&d);
        assert!(shut);
        let floor: Vec<&DisplayEdge> = d.edges.iter().filter(|e| e.points.iter().all(|p| (p.z - 4.0).abs() < 1e-9)).collect();
        assert_eq!(floor.len(), 1);
        let [lo, hi] = floor[0].turn.unwrap();
        assert!((lo + FRAC_PI_2).abs() < 1e-9 && (hi + FRAC_PI_2).abs() < 1e-9, "{lo} {hi}");
        let rest = d.edges.iter().filter(|e| !std::ptr::eq(*e, floor[0]));
        for e in rest {
            let [lo, hi] = e.turn.unwrap();
            assert!((lo - FRAC_PI_2).abs() < 1e-9 && (hi - FRAC_PI_2).abs() < 1e-9, "{lo} {hi}");
        }
    }

    #[test]
    fn a_rounded_edge_meets_its_sides_smoothly() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        // The upright edge at x = y = 0.
        let edge = model
            .unique_subshapes(b, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (a, z) = model.edge_vertices(e);
                let (p, q) = (model.point(a), model.point(z));
                p.x.abs() + p.y.abs() + q.x.abs() + q.y.abs() < 1e-12
            })
            .unwrap();
        let rounded = crate::blend::fillet(&mut model, b, &[edge], 2.0).unwrap();
        let d = display(&mut model, rounded, TOLERANCE).unwrap();
        let (shut, volume) = closed(&d);
        let exact = 1000.0 - (4.0 - PI) * 10.0;
        assert!(shut && (volume - exact).abs() <= 2.0 * TOLERANCE * 600.0, "{shut} {volume}");
        assert_eq!(d.faces.len(), 7);
        let smooth = turns(&d).into_iter().filter(|t| t[0].abs() < 1e-6 && t[1].abs() < 1e-6).count();
        let square = turns(&d).into_iter().filter(|t| (t[0] - FRAC_PI_2).abs() < 1e-6 && (t[1] - FRAC_PI_2).abs() < 1e-6).count();
        // The two where the round meets the sides; the rest square, the round's arcs too.
        assert_eq!((d.edges.len(), smooth, square), (15, 2, 13), "{:?}", turns(&d));
        assert!(d.edges.iter().all(|e| !e.split));
    }

    #[test]
    fn a_ball_draws_no_edge_and_has_the_normal_at_its_poles() {
        let mut model = Model::new();
        let s = Sphere::new(3.0).build(&mut model).unwrap();
        let d = display(&mut model, s, TOLERANCE).unwrap();
        assert!(d.edges.is_empty());
        let (shut, volume) = closed(&d);
        let exact = 4.0 / 3.0 * PI * 27.0;
        assert!(shut && (volume - exact).abs() <= 2.0 * TOLERANCE * 4.0 * PI * 9.0, "{shut} {volume}");
        check_faces(&d, |p| Some(((p - Point::ORIGIN).length() - 3.0).abs()).map(|off| (off, Direction::try_from_vector(p - Point::ORIGIN).unwrap())));
        // The poles are drawn as vertices of the face, with the axis as their normal.
        let poles = d.faces[0].positions.iter().filter(|p| p.x.hypot(p.y) < 1e-12).count();
        assert!(poles >= 2);
        // No seam has no length.
        assert!(d.seams.iter().all(|s| s.len() >= 2));
    }

    #[test]
    fn a_cones_apex_has_a_vertex_for_each_triangle_with_its_rulings_normal() {
        let mut model = Model::new();
        let c = Cone::new(4.0, 0.0, 6.0).build(&mut model).unwrap();
        let d = display(&mut model, c, TOLERANCE).unwrap();
        let (shut, volume) = closed(&d);
        let exact = PI * 16.0 * 6.0 / 3.0;
        assert!(shut && (volume - exact).abs() < 0.1, "{shut} {volume}");
        // One edge (the rim; the apex is not drawn), one seam.
        assert_eq!((d.edges.len(), d.seams.len()), (1, 1));
        let half = (4.0f64 / 6.0).atan();
        let [lo, hi] = d.edges[0].turn.unwrap();
        assert!((lo - (FRAC_PI_2 + half)).abs() < 1e-9 && (hi - (FRAC_PI_2 + half)).abs() < 1e-9, "{lo} {hi}");
        // The cone's outward normal at a point (off the apex): out and up.
        let ruling = |p: Point| {
            let r = Vector::new(p.x, p.y, 0.0).normalized();
            Direction::try_from_vector(r * half.cos() + Vector::from(Direction::Z) * half.sin()).unwrap()
        };
        let side = d.faces.iter().find(|f| f.normals.iter().any(|n| n.z() > 0.0)).unwrap();
        let apex = Point::new(0.0, 0.0, 6.0);
        let mut at_apex = 0;
        for t in &side.triangles {
            for (k, &i) in t.iter().enumerate() {
                let p = side.positions[i as usize];
                if p.distance(apex) > 1e-12 {
                    assert!(side.normals[i as usize].dot(ruling(p)) > 1.0 - 1e-9);
                    continue;
                }
                at_apex += 1;
                // Its own vertex, with the normal of the ruling halfway across.
                assert_eq!(side.triangles.iter().flatten().filter(|&&j| j == i).count(), 1);
                let (b, c) = (side.positions[t[(k + 1) % 3] as usize], side.positions[t[(k + 2) % 3] as usize]);
                let halfway = b + (c - b) * 0.5;
                assert!(side.normals[i as usize].dot(ruling(halfway)) > 1.0 - 1e-9);
            }
        }
        assert!(at_apex >= 3);
        check_faces(&d, |_| None);
    }

    #[test]
    fn a_face_cut_in_two_meets_its_other_half_on_a_split() {
        let mut model = Model::new();
        let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
        let other = Cuboid::new(10.0, 10.0, 10.0).with_frame(Frame::from_z(Point::new(10.0, 0.0, 0.0), Direction::Z)).build(&mut model).unwrap();
        // Joined, the long box's sides come in two faces each, unless the union merges them.
        let joined = crate::boolean::union(&mut model, b, other).unwrap();
        let d = display(&mut model, joined, TOLERANCE).unwrap();
        let (shut, volume) = closed(&d);
        assert!(shut && (volume - 2000.0).abs() < 1e-9);
        for e in &d.edges {
            let [lo, hi] = e.turn.unwrap();
            if e.split {
                assert!(lo.abs() < 1e-12 && hi.abs() < 1e-12);
                assert!(e.points.iter().all(|p| (p.x - 10.0).abs() < 1e-9));
            } else {
                assert!((lo - FRAC_PI_2).abs() < 1e-9 && (hi - FRAC_PI_2).abs() < 1e-9);
            }
        }
        // Where the sides stay in two, the long box's edges are drawn in two too.
        let splits = d.edges.iter().filter(|e| e.split).count();
        assert!(splits == 0 || splits == 4, "{splits} splits");
        let length: f64 = d.edges.iter().filter(|e| !e.split).flat_map(|e| e.points.windows(2).map(|w| w[0].distance(w[1]))).sum();
        assert!((length - 160.0).abs() < 1e-9, "{length}");
    }
}

