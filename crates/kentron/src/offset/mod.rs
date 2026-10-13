//! Offsets: curves and shapes moved a distance away from themselves.
//!
//! * [`planar_wire`]: a wire offset within its plane.
//! * [`shape`]: a solid or shell whose faces are moved along their normals.
//! * [`thick_solid`]: a solid hollowed out, its walls of a given thickness, with some faces
//!   removed to open it.

pub(crate) mod shape;
mod wire;

pub use shape::{shape, shape_with_faces, thick_solid};
pub use wire::planar_wire;

/// How offsets are joined where they part (at corners of a wire turning away from the offset
/// side, at edges of a shape between faces whose offsets move apart).
#[doc(alias = "GeomAbs_JoinType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Join {
    /// By an arc around the corner (a rolling-ball face along the edge). Wires only, and
    /// shapes where no offset faces part.
    #[default]
    Arc,
    /// By extending the offsets until they meet (along their tangents, for wires).
    Intersection,
}

