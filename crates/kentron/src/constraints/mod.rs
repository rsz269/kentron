//! Geometric constraint solving: placements of assembly parts, and the lines and arcs of a
//! sketch, that best satisfy constraints between them (CadQuery's `Assembly.solve` and
//! `Sketch.solve`), as CadQuery's solvers solve them: IPOPT's Newton steps for assemblies,
//! nlopt's SLSQP for sketches.

pub mod assembly;
pub mod sketch;
mod dual;
mod least_squares;
mod newton;
pub(crate) mod slsqp;
mod sketch_costs;
