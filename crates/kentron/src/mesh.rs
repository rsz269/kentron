//! Polygonal approximations of geometry. Node indices are 0-based.

use crate::math::{AffineTransform, AffineTransformable, Direction, Point, Point2};

/// A triangle mesh approximating a surface.
#[doc(alias = "Poly_Triangulation")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TriangleMesh {
    nodes: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    normals: Option<Vec<Direction>>,
    deflection: f64,
    /// The deflection it was made for, where known (it may come out coarser).
    requested: Option<f64>,
}

impl TriangleMesh {
    /// Panics if a triangle refers to a missing node.
    pub fn new(nodes: Vec<Point>, triangles: Vec<[u32; 3]>) -> Self {
        assert!(
            triangles.iter().flatten().all(|&i| (i as usize) < nodes.len()),
            "triangle node index out of range"
        );
        Self { nodes, triangles, normals: None, deflection: 0.0, requested: None }
    }

    /// Sets one normal per node; panics on a count mismatch.
    pub fn with_normals(self, normals: Vec<Direction>) -> Self {
        assert_eq!(normals.len(), self.nodes.len(), "one normal per node");
        Self { normals: Some(normals), ..self }
    }

    /// Sets the maximal distance between the mesh and the exact geometry.
    pub fn with_deflection(self, deflection: f64) -> Self {
        Self { deflection, ..self }
    }

    pub fn nodes(&self) -> &[Point] {
        &self.nodes
    }

    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The total area of the triangles.
    #[doc(alias = "bg_trimesh_area")]
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| crate::math::triangle_area(self.nodes[t[0] as usize], self.nodes[t[1] as usize], self.nodes[t[2] as usize]))
            .sum()
    }

    pub fn normals(&self) -> Option<&[Direction]> {
        self.normals.as_deref()
    }

    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Records the deflection the mesh was made for (which it may not reach).
    pub fn with_requested(self, deflection: f64) -> Self {
        Self { requested: Some(deflection), ..self }
    }

    /// The deflection the mesh was made for, where known.
    pub fn requested(&self) -> Option<f64> {
        self.requested
    }
}

/// A 3D polyline approximating a curve, optionally with the curve parameter of each node.
#[doc(alias = "Poly_Polygon3D")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polyline {
    nodes: Vec<Point>,
    parameters: Option<Vec<f64>>,
    deflection: f64,
}

impl Polyline {
    pub fn new(nodes: Vec<Point>) -> Self {
        Self { nodes, parameters: None, deflection: 0.0 }
    }

    /// Panics on a count mismatch.
    pub fn with_parameters(self, parameters: Vec<f64>) -> Self {
        assert_eq!(parameters.len(), self.nodes.len(), "one parameter per node");
        Self { parameters: Some(parameters), ..self }
    }

    pub fn with_deflection(self, deflection: f64) -> Self {
        Self { deflection, ..self }
    }

    pub fn nodes(&self) -> &[Point] {
        &self.nodes
    }

    pub fn parameters(&self) -> Option<&[f64]> {
        self.parameters.as_deref()
    }

    pub fn deflection(&self) -> f64 {
        self.deflection
    }
}

/// A 2D polyline approximating a curve in a surface's parameter space.
#[doc(alias = "Poly_Polygon2D")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polyline2 {
    nodes: Vec<Point2>,
    deflection: f64,
}

impl Polyline2 {
    pub fn new(nodes: Vec<Point2>) -> Self {
        Self { nodes, deflection: 0.0 }
    }

    pub fn nodes(&self) -> &[Point2] {
        &self.nodes
    }
}

/// An edge's discretization given as node indices into a face's [`TriangleMesh`],
/// optionally with the curve parameter of each node.
#[doc(alias = "Poly_PolygonOnTriangulation")]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MeshPolyline {
    nodes: Vec<u32>,
    parameters: Option<Vec<f64>>,
    deflection: f64,
}

impl MeshPolyline {
    pub fn new(nodes: Vec<u32>) -> Self {
        Self { nodes, parameters: None, deflection: 0.0 }
    }

    pub fn with_deflection(self, deflection: f64) -> Self {
        Self { deflection, ..self }
    }

    pub fn nodes(&self) -> &[u32] {
        &self.nodes
    }

    pub fn parameters(&self) -> Option<&[f64]> {
        self.parameters.as_deref()
    }

    pub fn deflection(&self) -> f64 {
        self.deflection
    }
}

/// Nodes are mapped; normals are mapped by the inverse transpose (and kept when the map
/// collapses them); a reflection reverses each triangle's node order so that the winding
/// still agrees with the normals; the deflection grows by the map's largest stretch.
impl AffineTransformable for TriangleMesh {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        let reflects = t.matrix().determinant() < 0.0;
        let triangles = match reflects {
            true => self.triangles.iter().map(|&[a, b, c]| [a, c, b]).collect(),
            false => self.triangles.clone(),
        };
        let normals = self.normals.as_ref().map(|n| n.iter().map(|&d| t.transform_normal(d).unwrap_or(d)).collect());
        Self {
            nodes: self.nodes.iter().map(|&p| t.transform_point(p)).collect(),
            triangles,
            normals,
            deflection: self.deflection * t.max_stretch(),
            // (What it was made for is of the shape it was made on.)
            requested: None,
        }
    }
}

/// Nodes are mapped, parameters kept, the deflection grows by the largest stretch.
impl AffineTransformable for Polyline {
    fn affine_transformed(&self, t: &AffineTransform) -> Self {
        Self {
            nodes: self.nodes.iter().map(|&p| t.transform_point(p)).collect(),
            parameters: self.parameters.clone(),
            deflection: self.deflection * t.max_stretch(),
        }
    }
}

impl Polyline2 {
    /// The polyline with each node mapped by `f`.
    pub(crate) fn mapped(&self, f: impl Fn(Point2) -> Point2) -> Self {
        Self { nodes: self.nodes.iter().map(|&p| f(p)).collect(), deflection: self.deflection }
    }
}

/// Which points of a grid lie inside a closed mesh, by the parity of its crossings along
/// rays in +Z: `out[(i * ys.len() + j) * zs.len() + k]` for the point `(xs[i], ys[j],
/// zs[k])`. Each column's crossings are where the triangles over it (their projections
/// closed, edges included) meet it; between successive pairs of them, sorted, the points
/// from the first (included) to the second (excluded) are inside. `xs`, `ys` and `zs` are
/// increasing; columns through edges or vertices of the mesh count a crossing for each
/// triangle there, so grids are best offset from the mesh's features. Degenerate triangles
/// (seen edge-on from +Z) are passed over.
pub fn voxelize(nodes: &[Point], triangles: &[[u32; 3]], xs: &[f64], ys: &[f64], zs: &[f64]) -> Vec<bool> {
    let (nx, ny, nz) = (xs.len(), ys.len(), zs.len());
    // The first index of `values` at or above `v`.
    let at_or_above = |values: &[f64], v: f64| values.partition_point(|&x| x < v);
    let mut hits: Vec<Vec<f64>> = vec![Vec::new(); nx * ny];
    for t in triangles {
        let [a, b, c] = t.map(|i| nodes[i as usize]);
        let lo = [a.x.min(b.x).min(c.x), a.y.min(b.y).min(c.y)];
        let hi = [a.x.max(b.x).max(c.x), a.y.max(b.y).max(c.y)];
        let (i0, i1) = (at_or_above(xs, lo[0]), at_or_above(xs, hi[0]));
        let (j0, j1) = (at_or_above(ys, lo[1]), at_or_above(ys, hi[1]));
        if i0 >= i1 || j0 >= j1 {
            continue;
        }
        let d = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
        if d.abs() < 1e-12 {
            continue;
        }
        for i in i0..i1 {
            for j in j0..j1 {
                let (x, y) = (xs[i], ys[j]);
                let l1 = ((b.y - c.y) * (x - c.x) + (c.x - b.x) * (y - c.y)) / d;
                let l2 = ((c.y - a.y) * (x - c.x) + (a.x - c.x) * (y - c.y)) / d;
                let l3 = 1.0 - l1 - l2;
                if l1 >= 0.0 && l2 >= 0.0 && l3 >= 0.0 {
                    hits[i * ny + j].push(l1 * a.z + l2 * b.z + l3 * c.z);
                }
            }
        }
    }
    let mut out = vec![false; nx * ny * nz];
    for (column, h) in hits.iter_mut().enumerate() {
        h.sort_by(f64::total_cmp);
        for pair in h.chunks_exact(2) {
            let (k0, k1) = (at_or_above(zs, pair[0]), at_or_above(zs, pair[1]));
            if k0 < k1 {
                out[column * nz + k0..column * nz + k1].fill(true);
            }
        }
    }
    out
}

#[cfg(test)]
mod voxelize_tests {
    use super::*;

    #[test]
    fn a_box_fills_the_points_inside_it() {
        // A unit cube's mesh (12 triangles), on a grid offset from its features.
        let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let nodes = vec![p(0., 0., 0.), p(1., 0., 0.), p(1., 1., 0.), p(0., 1., 0.), p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)];
        let triangles = vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
        // (x and y offset differently: a column with x = y would run along the faces'
        // diagonals.)
        let grid = |offset: f64| -> Vec<f64> { (0..12).map(|i| -0.4 + offset + 0.163 * i as f64).collect() };
        let (xs, ys, zs) = (grid(0.022), grid(0.011), grid(0.0));
        let out = voxelize(&nodes, &triangles, &xs, &ys, &zs);
        let n = xs.len();
        for (i, &x) in xs.iter().enumerate() {
            for (j, &y) in ys.iter().enumerate() {
                for (k, &z) in zs.iter().enumerate() {
                    let inside = [x, y].iter().all(|&c| c > 0.0 && c < 1.0) && z >= 0.0 && z < 1.0;
                    assert_eq!(out[(i * n + j) * n + k], inside, "at {x} {y} {z}");
                }
            }
        }
    }
}
