//! `src/libbg/tests/trimesh_area.c` (`bg_trimesh_area`): a square of two triangles has the
//! square's area, and a bent one does not.

use crate::p;
use kentron::mesh::TriangleMesh;

#[test]
fn bg_trimesh_area() {
    let faces = vec![[0, 1, 2], [1, 2, 3]];
    for i in 0..10 {
        let k = i as f64;
        let square = TriangleMesh::new(vec![p(0.0, 0.0, 7.0), p(0.0, k, 7.0), p(k, k, 7.0), p(k, 0.0, 7.0)], faces.clone());
        assert!((square.area() - k * k).abs() <= 1e-12, "{} != {}", square.area(), k * k);
        if i == 0 {
            continue;
        }
        let bent = TriangleMesh::new(vec![p(0.0, 0.0, 7.0), p(0.0, k, 7.0), p(k, k, 7.0), p(k, 0.0, 10.0)], faces.clone());
        assert!((bent.area() - k * k).abs() > 1e-12);
    }
}
