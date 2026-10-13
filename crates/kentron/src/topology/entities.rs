//! Topology entities and the "use" entries that place them in their parents.
//!
//! An entity holds what is intrinsic to it (its geometry slot, its tolerance). How a
//! parent uses a child (orientation, placement) lives on a separate use entry, so a child
//! shared by several parents (an edge between two faces) is stored once.

use crate::math::{Point, Transform};

use super::id::{Id, NodeId};
use super::slots::{EdgeCurve, EdgeMeshPolyline, EdgePolyline, FaceMesh, FaceSurface, PCurve,
    PCurvePolyline};

pub type VertexId = Id<Vertex>;
pub type EdgeId = Id<Edge>;
pub type CoEdgeId = Id<CoEdge>;
pub type WireId = Id<Wire>;
pub type FaceId = Id<Face>;
pub type ShellId = Id<Shell>;
pub type SolidId = Id<Solid>;
pub type CompoundId = Id<Compound>;
pub type CompSolidId = Id<CompSolid>;
pub type ProductId = Id<Product>;
pub type InstanceId = Id<Instance>;

pub type VertexUseId = Id<VertexUse>;
pub type WireUseId = Id<WireUse>;
pub type FaceUseId = Id<FaceUse>;
pub type ShellUseId = Id<ShellUse>;
pub type SolidUseId = Id<SolidUse>;
pub type CompoundChildId = Id<CompoundChild>;
pub type ProductChildId = Id<ProductChild>;

/// How a child is oriented within its parent.
#[doc(alias = "TopAbs_Orientation")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Orientation {
    #[default]
    Forward,
    Reversed,
    /// Lies inside the parent's material (e.g. an edge embedded in a face).
    Internal,
    /// Lies outside the parent's material.
    External,
}

impl Orientation {
    /// Forward and reversed swap; internal and external are unchanged.
    pub fn reversed(self) -> Orientation {
        match self {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            other => other,
        }
    }

    /// The orientation of a child oriented `self` within a parent that is itself
    /// oriented `parent` (OCCT's `TopAbs::Compose`).
    #[doc(alias = "TopAbs::Compose")]
    pub fn compose(self, parent: Orientation) -> Orientation {
        match (parent, self) {
            (Orientation::Forward, o) => o,
            (Orientation::Reversed, o) => o.reversed(),
            (Orientation::Internal, _) => Orientation::Internal,
            (Orientation::External, _) => Orientation::External,
        }
    }
}

// --- entities -----------------------------------------------------------------------------

#[doc(alias = "TopoDS_Vertex")]
#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    pub point: Point,
    pub tolerance: f64,
}

#[doc(alias = "TopoDS_Edge")]
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    /// Invalid for a degenerate edge (collapsed to a point, e.g. at a sphere's pole).
    pub curve: Id<EdgeCurve>,
    pub polyline: Id<EdgePolyline>,
    pub tolerance: f64,
    pub start: VertexUseId,
    pub end: VertexUseId,
}

/// One side of an edge as used by a wire: the half-edge of a face boundary. Owns the edge's
/// curve in the face's parameter space; a seam edge (where a closed surface meets itself)
/// has two coedges on the same face, with opposite orientations.
#[derive(Debug, Clone, PartialEq)]
pub struct CoEdge {
    pub wire: WireId,
    pub edge: EdgeId,
    /// Invalid for a wire that does not bound a face.
    pub face: FaceId,
    /// Orientation relative to the edge's curve.
    pub orientation: Orientation,
    pub pcurve: Id<PCurve>,
    pub pcurve_polyline: Id<PCurvePolyline>,
    pub mesh_polyline: Id<EdgeMeshPolyline>,
}

/// A connected chain of coedges; on a face, a boundary loop.
#[doc(alias = "TopoDS_Wire")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Wire {
    /// Coedges in order along the wire.
    pub coedges: Vec<CoEdgeId>,
}

/// A bounded region of a surface.
#[doc(alias = "TopoDS_Face")]
#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub surface: Id<FaceSurface>,
    pub mesh: Id<FaceMesh>,
    pub tolerance: f64,
    pub wires: Vec<WireUseId>,
}

/// Faces connected along shared edges.
#[doc(alias = "TopoDS_Shell")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Shell {
    pub faces: Vec<FaceUseId>,
}

/// A region of space bounded by shells.
#[doc(alias = "TopoDS_Solid")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Solid {
    pub shells: Vec<ShellUseId>,
}

/// An arbitrary collection of nodes.
#[doc(alias = "TopoDS_Compound")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Compound {
    pub children: Vec<CompoundChildId>,
}

/// Solids sharing faces.
#[doc(alias = "TopoDS_CompSolid")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CompSolid {
    pub solids: Vec<SolidUseId>,
}

/// A reusable definition: a part (with a topology root) or an assembly of placed
/// instances of other products.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Product {
    /// The topology of a part (its root, oriented and placed within the part); `None` for
    /// an assembly.
    pub shape: Option<super::Shape>,
    pub children: Vec<ProductChildId>,
}

impl Product {
    /// True for a product with a shape.
    pub fn is_part(&self) -> bool {
        self.shape.is_some()
    }

    /// True for a product holding instances. A new empty product is neither a part nor an
    /// assembly.
    pub fn is_assembly(&self) -> bool {
        !self.children.is_empty()
    }
}

/// A placed occurrence of a product (or of a topology root) inside a parent product.
#[derive(Debug, Clone, PartialEq)]
pub struct Instance {
    pub target: NodeId,
}

// --- uses ---------------------------------------------------------------------------------

/// A vertex at one end of an edge.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexUse {
    pub edge: EdgeId,
    pub vertex: VertexId,
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WireUse {
    pub face: FaceId,
    pub wire: WireId,
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FaceUse {
    pub shell: ShellId,
    pub face: FaceId,
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ShellUse {
    pub solid: SolidId,
    pub shell: ShellId,
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SolidUse {
    pub comp_solid: CompSolidId,
    pub solid: SolidId,
    pub orientation: Orientation,
}

/// A node placed in a compound.
#[derive(Debug, Clone, PartialEq)]
pub struct CompoundChild {
    pub compound: CompoundId,
    pub child: NodeId,
    pub orientation: Orientation,
    pub placement: Transform,
}

/// An instance placed in a product.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductChild {
    pub product: ProductId,
    pub instance: InstanceId,
    pub placement: Transform,
}
