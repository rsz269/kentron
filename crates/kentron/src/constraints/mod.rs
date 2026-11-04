//! Geometric constraint solving: placements of assembly parts, and the lines and arcs of a
//! sketch, that best satisfy constraints between them (CadQuery's `Assembly.solve` and
//! `Sketch.solve`), by least squares with exact derivatives.

pub mod assembly;
pub mod sketch;
mod dual;
mod least_squares;
mod newton;
