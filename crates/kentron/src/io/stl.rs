//! STL files: triangles, each with its normal, in binary or ASCII form.

use std::io::Write;

use crate::error::{Error, Result};
use crate::math::Point;

/// Writes triangles (indices into `nodes`, counter-clockwise seen from outside).
#[doc(alias = "StlAPI_Writer")]
#[doc(alias = "writestl")]
pub fn write(nodes: &[Point], triangles: &[[u32; 3]], out: &mut dyn Write, ascii: bool) -> Result<()> {
    let io = |e: std::io::Error| Error::InvalidInput(format!("writing STL: {e}"));
    let normal = |t: &[u32; 3]| {
        let [a, b, c] = t.map(|i| nodes[i as usize]);
        let n = (b - a).cross(c - a);
        let l = n.length();
        if l > 0.0 { n * (1.0 / l) } else { n }
    };
    if ascii {
        let mut text = String::from("solid kentron\n");
        for t in triangles {
            let n = normal(t);
            text.push_str(&format!(" facet normal {:e} {:e} {:e}\n  outer loop\n", n.x, n.y, n.z));
            for &i in t {
                let p = nodes[i as usize];
                text.push_str(&format!("   vertex {:e} {:e} {:e}\n", p.x, p.y, p.z));
            }
            text.push_str("  endloop\n endfacet\n");
        }
        text.push_str("endsolid kentron\n");
        return out.write_all(text.as_bytes()).map_err(io);
    }
    let mut bytes = Vec::with_capacity(84 + 50 * triangles.len());
    let mut header = [0u8; 80];
    let name = b"kentron";
    header[..name.len()].copy_from_slice(name);
    bytes.extend_from_slice(&header);
    bytes.extend_from_slice(&(triangles.len() as u32).to_le_bytes());
    for t in triangles {
        let n = normal(t);
        for v in [n.x, n.y, n.z] {
            bytes.extend_from_slice(&(v as f32).to_le_bytes());
        }
        for &i in t {
            let p = nodes[i as usize];
            for v in [p.x, p.y, p.z] {
                bytes.extend_from_slice(&(v as f32).to_le_bytes());
            }
        }
        bytes.extend_from_slice(&[0, 0]);
    }
    out.write_all(&bytes).map_err(io)
}

#[cfg(test)]
mod tests {
    use super::write;
    use crate::math::Point;

    #[test]
    fn binary_and_ascii_files_hold_the_triangles() {
        let nodes = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let mut binary = Vec::new();
        write(&nodes, &[[0, 1, 2]], &mut binary, false).unwrap();
        assert_eq!(binary.len(), 80 + 4 + 50);
        assert!(binary.starts_with(b"kentron"));
        assert_eq!(u32::from_le_bytes(binary[80..84].try_into().unwrap()), 1);

        let mut ascii = Vec::new();
        write(&nodes, &[[0, 1, 2]], &mut ascii, true).unwrap();
        let text = String::from_utf8(ascii).unwrap();
        assert!(text.starts_with("solid kentron\n") && text.ends_with("endsolid kentron\n"));
        assert_eq!(text.matches("facet normal").count(), 1);
    }
}
