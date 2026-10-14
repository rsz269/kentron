//! kentron: a CAD kernel.
//!
//! Modules, from the bottom up:
//!
//! * [`math`]: points, vectors, directions, axes, frames, matrices, transformations.
//! * [`tolerance`]: default tolerances for geometric comparisons.
//! * [`geometry`]: curves and surfaces as plain values.
//! * [`mesh`]: polygonal approximations.
//! * [`topology`]: the boundary representation, stored in a [`topology::Model`] arena.
//! * Modeling, on a model: [`make`] (vertices, edges, wires, faces), [`primitives`]
//!   (boxes, cylinders, ...), [`sweep`] (extrusion, revolution), [`loft`] (shells and solids
//!   through sections), [`offset`] (offset wires and shapes, thick solids), [`boolean`] (union, difference,
//!   intersection), [`sew`] (joining faces into shells), [`blend`] (fillets), [`draft`] (draft
//!   angles), [`feature`] (bosses and pockets), [`convert`], [`modify`] (transforming shapes),
//!   [`sketch`] (profiles drawn by a moving pen).
//! * Analysis: [`properties`] (length, area, volume, bounding box), [`validity`]; repair:
//!   [`heal`].
//! * [`meshing`]: triangulation; [`io`]: reading and writing shapes.
//!
//! The functionality is modelled on Open CASCADE Technology (OCCT), whose tests are being
//! ported in `crates/occt-tests`, but names and structure are kentron's own. Types carry
//! `#[doc(alias)]`es for their OCCT counterparts, and `docs/occt-mapping.md` maps OCCT
//! concepts to kentron ones.
//!
//! # Conventions
//!
//! * Small value types (`Point`, `Direction`, `Frame`, `Line`, ...) are `Copy`.
//! * Constructors that can fail on input return `Option` (`Direction::try_new`), with a
//!   panicking variant for inputs known to be valid (`Direction::new`).
//! * Everything movable implements [`math::Transformable`] (or `Transformable2`).
//! * Operations that can fail return [`Result`] with the crate's [`Error`].
//! * No reference counting: shared data lives in a [`topology::Model`] and is referenced
//!   by typed ids. See [`topology`] for the rationale.

// Shim parameters and fields are unused until implemented.
#![allow(unused_variables, dead_code)]

pub mod blend;
pub mod boolean;
pub mod convert;
pub mod distance;
pub mod divide;
pub mod draft;
pub mod feature;
mod error;
pub mod geometry;
pub mod hashing;
pub mod heal;
pub mod helix;
pub mod io;
pub mod loft;
pub mod make;
pub mod offset;
pub mod math;
mod parallel;
pub mod mesh;
pub mod meshing;
pub mod modify;
pub mod primitives;
pub mod project;
pub mod properties;
pub mod sew;
pub mod sketch;
pub mod sweep;
#[cfg(feature = "text")]
pub mod text;
pub mod tolerance;
pub mod topology;
pub mod unify;
pub mod validity;
pub mod voxels;

pub use error::{Error, Result};
