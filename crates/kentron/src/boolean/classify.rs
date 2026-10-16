//! Where pieces of one operand's boundary lie with respect to the other operand.

use crate::math::{Point, Vector};

use super::geom::{self, on_tol};
use super::graph::{Graph, Location};

/// Where a piece of one operand's boundary lies with respect to the other operand.
#[doc(alias = "TopAbs_State")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum State {
    Inside,
    Outside,
    /// On the other operand's boundary, with the same outward normal.
    OnSame,
    /// On the other operand's boundary, with the opposite outward normal.
    OnOpposite,
}

/// Directions for rays, chosen to avoid axis-aligned geometry; a ray grazing an edge or
/// a surface is ambiguous and replaced by the next.
const DIRECTIONS: [(f64, f64, f64); 7] = [
    (0.5773, 0.3141, 0.7536),
    (-0.2718, 0.8414, 0.4142),
    (0.6931, -0.5403, 0.1732),
    (-0.1414, -0.3606, -0.9220),
    (0.8660, 0.1050, -0.4890),
    (0.1234, -0.9123, 0.3901),
    (-0.7071, 0.2236, -0.6708),
];

impl Graph {
    /// Where `p`, inside a piece with outward normal `outward`, lies with respect to
    /// operand `operand`.
    pub fn classify(&self, p: Point, outward: Vector, operand: usize) -> State {
        let faces = self.faces_of(operand);
        for g in self.index().faces_containing(operand, p) {
            let face = &self.faces[g];
            if geom::signed_distance(&face.surface, p).abs() <= on_tol()
                && face.locate(p, &self.edges, &self.vertices) != Location::Outside
            {
                return if outward.dot(face.outward_normal(p)) > 0.0 { State::OnSame } else { State::OnOpposite };
            }
        }
        let inside = if self.open[operand] { self.behind(faces, p) } else { self.contains(operand, faces, p) };
        if inside { State::Inside } else { State::Outside }
    }

    /// For a half-space (an open operand): true if `p` is behind its boundary, on the side
    /// opposite the outward normal of the nearest face's surface.
    pub(super) fn behind(&self, faces: &[usize], p: Point) -> bool {
        let nearest = faces
            .iter()
            .map(|&g| {
                let face = &self.faces[g];
                let foot = geom::surface_point(&face.surface, geom::uv_of(&face.surface, p));
                (foot.distance(p), g, foot)
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match nearest {
            Some((_, g, foot)) => (p - foot).dot(self.faces[g].outward_normal(foot)) < 0.0,
            None => false,
        }
    }

    /// How many times the faces of `operand` wind round `p` (off them): along a ray from it,
    /// the faces it leaves through less those it enters by (1 inside a solid, 0 outside it;
    /// more inside a shell passing through itself). `None` where every ray grazes an edge,
    /// or the operand has far faces.
    pub(super) fn winding(&self, operand: usize, p: Point) -> Option<i32> {
        if self.far.iter().any(|f| f.operand == operand) {
            return None;
        }
        let faces = self.faces_of(operand);
        'rays: for (x, y, z) in DIRECTIONS {
            let d = Vector::new(x, y, z).normalized();
            let mut winding = 0;
            for g in self.index().faces_along_ray(operand, p, d).into_iter().filter(|g| faces.binary_search(g).is_ok()) {
                let face = &self.faces[g];
                let Some(hits) = geom::ray_hits(&face.surface, p, d, &face.bbox) else { continue 'rays };
                for t in hits {
                    let at = p + d * t;
                    match face.locate(at, &self.edges, &self.vertices) {
                        Location::Inside => {
                            let along = face.outward_normal(at).dot(d);
                            // (Grazing the face: ambiguous.)
                            if along.abs() <= 1e-9 {
                                continue 'rays;
                            }
                            winding += if along > 0.0 { 1 } else { -1 };
                        }
                        Location::Boundary => continue 'rays,
                        Location::Outside => {}
                    }
                }
            }
            return Some(winding);
        }
        None
    }

    /// True if `p` (off their boundary) is inside the solid bounded by `faces` (faces of
    /// `operand`) and the operand's far faces.
    pub(super) fn contains(&self, operand: usize, faces: &[usize], p: Point) -> bool {
        // (Its faces may all be far ones, left out of the graph.)
        if faces.is_empty() && !self.far.iter().any(|f| f.operand == operand) {
            return false;
        }
        // Outside the operand's (tight) box, outside the solid: rays from just beside a
        // nearly coincident surface may miscount.
        if let (Some((lo, hi)), 3, false) = (self.tight_boxes[operand], self.dimension[operand], self.open[operand]) {
            let margin = on_tol();
            if p.x < lo.x - margin || p.x > hi.x + margin || p.y < lo.y - margin || p.y > hi.y + margin || p.z < lo.z - margin || p.z > hi.z + margin {
                return false;
            }
        }
        'rays: for (x, y, z) in DIRECTIONS {
            let d = Vector::new(x, y, z).normalized();
            let mut crossings = 0;
            // Only faces whose box the ray crosses.
            let along = self.index().faces_along_ray(operand, p, d);
            for g in along.into_iter().filter(|g| faces.binary_search(g).is_ok()) {
                let face = &self.faces[g];
                let Some(hits) = geom::ray_hits(&face.surface, p, d, &face.bbox) else { continue 'rays };
                for t in hits {
                    match face.locate(p + d * t, &self.edges, &self.vertices) {
                        Location::Inside => crossings += 1,
                        Location::Boundary => continue 'rays,
                        Location::Outside => {}
                    }
                }
            }
            // The operand's faces left out of the graph count too.
            for g in self.index().far_along_ray(operand, p, d) {
                let far = &self.far[g];
                let data = &far.data;
                // (The ray where the face is unplaced: a placed far face's data is unplaced.)
                let (p, d) = if far.placement == crate::math::Transform::IDENTITY {
                    (p, d)
                } else {
                    use crate::math::Transformable;
                    let inverse = far.placement.inverse();
                    (p.transformed(&inverse), d.transformed(&inverse))
                };
                let Some(hits) = geom::ray_hits(&data.surface, p, d, &data.sample.bbox) else { continue 'rays };
                for t in hits {
                    match data.locate(p + d * t) {
                        Location::Inside => crossings += 1,
                        Location::Boundary => continue 'rays,
                        Location::Outside => {}
                    }
                }
            }
            return crossings % 2 == 1;
        }
        false
    }
}
