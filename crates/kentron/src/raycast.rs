//! Rays cast at a shape: where a ray meets the shape's faces, as distances along it.
//!
//! As OCCT's `IntCurvesFace_ShapeIntersector` does (`Load` the shape once, `Perform` a
//! `gp_Lin` over a range of its parameter): every face of the shape is made ready for lines
//! once (`IntCurvesFace_Intersector`: a closed-form solution for a plane, cylinder, cone,
//! sphere or torus, a polyhedron sampled from any other surface, see
//! `geometry::int_line_surface`), and each ray is intersected with each face's surface over
//! the face's UV box, the points the face's classifier finds inside it or on its boundary
//! kept (`InternalCall`), all faces' points sorted by their distance (`SortResult`).
//!
//! Two departures. The classifier is kentron's (the boolean's face classifier,
//! `PointClassifier::locate_on_faces`, which tells a point of the face's surface inside,
//! on its boundary within `1e-6`, or outside, as `BRepClass_FaceClassifier` does) where OCCT
//! classifies the point's parameters with `BRepTopAdaptor_TopolTool::Classify` at the tolerance
//! of the face's edges. And the faces a ray cannot reach are passed over by their boxes (a
//! tree of the faces' bounding boxes, enlarged), which OCCT tries one by one; a face whose
//! box the ray misses has no point on it.

use crate::boolean::{PointClassifier, PointLocation};
use crate::error::{Error, Result};
use crate::geometry::int_line_surface::{LineFace, INFINITE};
use crate::math::{Bvh, Point, Transform, Transformable, Vector};
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// How far the faces' boxes are enlarged for the rays: beyond the classifier's boundary
/// tolerance.
const BOX_GAP: f64 = 1e-5;

/// A shape made ready for casting rays at (`IntCurvesFace_ShapeIntersector::Load`).
#[doc(alias = "IntCurvesFace_ShapeIntersector")]
pub struct Raycaster {
    faces: Vec<RayFace>,
    tree: Bvh,
}

struct RayFace {
    line: LineFace,
    classifier: PointClassifier,
}

impl Raycaster {
    /// The faces of `shape`, each made ready for lines.
    pub fn new(model: &Model, shape: Shape) -> Result<Raycaster> {
        let mut faces = Vec::new();
        let mut boxes = Vec::new();
        for face in model.explore(shape, NodeKind::Face).collect::<Vec<_>>() {
            let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..face };
            let Some(surface) = model.surface(local).cloned() else { continue };
            let surface = surface.transformed(&face.placement);
            let bounds = model.face_uv_bounds(face);
            let classifier = PointClassifier::new(model, face)?;
            boxes.push(crate::properties::bounding_box(model, face).enlarged(BOX_GAP));
            faces.push(RayFace { line: LineFace::new(surface, bounds), classifier });
        }
        Ok(Raycaster { tree: Bvh::new(&boxes), faces })
    }

    /// `Perform(gp_Lin(origin, direction), 0, Precision::Infinite())`: the distances from
    /// `origin` along `direction` (which need not be a unit vector) to every point where the
    /// ray meets a face, in increasing order. A ray through an edge meets every face there,
    /// and lists the point once for each. Errs for a zero direction, and where a face cannot
    /// be intersected (a sampled surface with no extent along one of its directions, which
    /// OCCT's intersector is not done with).
    pub fn cast(&self, origin: Point, direction: Vector) -> Result<Vec<f64>> {
        let len = direction.length();
        if !(len > 0.0) || !len.is_finite() {
            return Err(Error::InvalidInput("a ray without a direction".into()));
        }
        let unit = direction * (1.0 / len);
        let mut found = Vec::new();
        for f in self.tree.along_ray(origin, unit) {
            let face = &self.faces[f];
            let Some((points, parinf, parsup)) = face.line.perform(origin, unit, 0.0, INFINITE) else {
                return Err(Error::Failed("a face without extent along one of its directions".into()));
            };
            // `InternalCall`: the points the face holds, within the range.
            for p in points {
                if face.classifier.locate_on_faces(p.point) != PointLocation::Outside && p.w >= parinf && p.w <= parsup {
                    found.push(p.w);
                }
            }
        }
        found.sort_by(f64::total_cmp);
        Ok(found)
    }

    /// [`Raycaster::cast`] for each ray, `(origin, direction)`, in the rays' order.
    pub fn cast_all(&self, rays: &[(Point, Vector)]) -> Result<Vec<Vec<f64>>> {
        crate::parallel::map(rays, |&(o, d)| self.cast(o, d)).into_iter().collect()
    }
}

/// The distances from `origin` along `direction` (not necessarily a unit vector) to every
/// point where the ray meets a face of `shape`, in increasing order (see [`Raycaster`], which
/// prepares the shape once for many rays).
#[doc(alias = "IntCurvesFace_ShapeIntersector")]
pub fn raycast(model: &Model, shape: Shape, origin: Point, direction: Vector) -> Result<Vec<f64>> {
    Raycaster::new(model, shape)?.cast(origin, direction)
}

/// [`raycast`] for each ray, `(origin, direction)`, in the rays' order: the shape prepared
/// once.
#[doc(alias = "IntCurvesFace_ShapeIntersector")]
pub fn raycast_all(model: &Model, shape: Shape, rays: &[(Point, Vector)]) -> Result<Vec<Vec<f64>>> {
    Raycaster::new(model, shape)?.cast_all(rays)
}

#[cfg(test)]
mod tests;
