//! Geometry and mesh slots of topology entities.
//!
//! Each slot belongs to exactly one entity and holds its geometry by value; geometry is
//! never shared between entities, so nothing is reference counted. Identical geometry used
//! by several entities is duplicated (see [`Model::deduplicate`](super::Model::deduplicate)).

use crate::geometry::{curve::Curve, curve2::Curve2, surface::Surface};
use crate::mesh::{MeshPolyline, Polyline, Polyline2, TriangleMesh};

use super::entities::{CoEdgeId, EdgeId, FaceId};

/// The 3D curve of an edge, restricted to `[first, last]`.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeCurve {
    pub edge: EdgeId,
    pub curve: Curve,
    pub first: f64,
    pub last: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EdgePolyline {
    pub edge: EdgeId,
    pub polyline: Polyline,
}

/// The curve of a coedge in the (u, v) parameter space of its face's surface, restricted
/// to `[first, last]`.
#[derive(Debug, Clone, PartialEq)]
pub struct PCurve {
    pub coedge: CoEdgeId,
    pub curve: Curve2,
    pub first: f64,
    pub last: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PCurvePolyline {
    pub coedge: CoEdgeId,
    pub polyline: Polyline2,
}

/// A coedge's discretization as nodes of its face's mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeMeshPolyline {
    pub coedge: CoEdgeId,
    pub polyline: MeshPolyline,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FaceSurface {
    pub face: FaceId,
    pub surface: Surface,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FaceMesh {
    pub face: FaceId,
    pub mesh: TriangleMesh,
}
