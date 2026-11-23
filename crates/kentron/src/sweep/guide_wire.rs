//! Preparing a guided pipe's guide wire to run with the spine: an open guide along an open
//! spine turned to run the spine's way, a closed guide opened at the point nearest a given
//! one and given a direction, and the spine's reduced arc lengths (each edge's share of the
//! wire's length) that match spine edges to stretches of the guide.

use crate::error::{Error, Result};
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::ext_pc::ext_pc_other;
use crate::math::{Point, Vector};
use crate::topology::comp_curve::{EdgeWire, ShapeKeys, WireEdge, WireVertex};
use crate::topology::{Model, Orientation, Shape};

/// The reduced arc lengths of `wire`'s edges: entry `i ≥ 1` is the length up to the end of
/// edge `i` over the whole length; entry 0 is the whole length itself. Degenerate edges
/// count nothing. A wire shorter than `1e-7` keeps its lengths unreduced, but for the last
/// entry, 1.
pub(crate) fn compute_acr(wire: &EdgeWire, length: &dyn Fn(&GeomAdaptor) -> f64) -> Vec<f64> {
    let n = wire.edges.len();
    let mut acr = vec![0.0; n + 1];
    let mut total = 0.0;
    for (i, edge) in wire.edges.iter().enumerate() {
        if let (false, Some(c)) = (edge.degenerate, &edge.curve) {
            total += length(c);
        }
        acr[i + 1] = total;
    }
    acr[0] = total;
    if total > 1e-7 {
        for a in &mut acr[1..] {
            *a /= total;
        }
    }
    acr[n] = 1.0;
    acr
}

/// Where on a wire the point nearest a given one lies.
enum Nearest {
    /// At the start vertex of edge `i`.
    Vertex(usize),
    /// Inside edge `i`, at its curve parameter `t`.
    Edge(usize, f64, Point),
}

/// The point of `wire` nearest `p`: vertices first, in wire order, then edge interiors in
/// edge order, among points equally near within `1e-7`.
fn nearest(wire: &EdgeWire, p: Point) -> Result<Nearest> {
    let mut vertices: Vec<(f64, usize)> = Vec::new();
    let mut interiors: Vec<(f64, usize, f64, Point)> = Vec::new();
    for (i, edge) in wire.edges.iter().enumerate() {
        let Some(v) = edge.start else { return Err(Error::Failed("a guide edge without vertices".into())) };
        vertices.push((v.point.distance(p), i));
        let Some(curve) = &edge.curve else {
            if edge.degenerate {
                continue;
            }
            return Err(Error::Failed("a guide edge without a curve".into()));
        };
        if let Some(extrema) = ext_pc_other(curve, p, edge.first, edge.last, 1e-10) {
            for e in extrema.into_iter().filter(|e| e.is_min) {
                if e.parameter > edge.first && e.parameter < edge.last {
                    interiors.push((e.square_distance.sqrt(), i, e.parameter, e.point));
                }
            }
        }
    }
    let best = vertices.iter().map(|v| v.0).chain(interiors.iter().map(|e| e.0)).fold(f64::INFINITY, f64::min);
    if let Some(&(_, i)) = vertices.iter().find(|v| v.0 <= best + 1e-7) {
        return Ok(Nearest::Vertex(i));
    }
    let &(_, i, t, x) = interiors.iter().find(|e| e.0 <= best + 1e-7).expect("a nearest point");
    Ok(Nearest::Edge(i, t, x))
}

/// The closed wire `wire` of `model`, opened at its point nearest `p` and run the way of
/// `dir` there; `tolerance` is how near `p` must be to an edge's vertex to start there.
pub(crate) fn search_origin(model: &Model, wire: Shape, p: Point, dir: Vector, tolerance: f64) -> Result<EdgeWire> {
    let mut keys = ShapeKeys::new();
    let wire = EdgeWire::from_model(model, Shape { orientation: Orientation::Forward, ..wire }, &mut keys)?;
    search_origin_of(&wire, p, dir, tolerance, &mut keys)
}

/// [`search_origin`] on a wire already read.
///
/// Where the nearest point lies inside an edge, the edge's own (unrun) start and end decide
/// first: `p` within `tolerance` of the start makes the end the origin, and then `p` within
/// `tolerance` of the end, and nearer it than the start, makes the start the origin; else
/// the edge is cut at the nearest point. Run backwards, the wire begins with the start edge
/// reversed, so it begins at that edge's far end, not at the origin.
pub(crate) fn search_origin_of(wire: &EdgeWire, p: Point, dir: Vector, tolerance: f64, keys: &mut ShapeKeys) -> Result<EdgeWire> {
    if !wire.closed {
        return Err(Error::InvalidInput("searching the origin of an open wire".into()));
    }
    if !(dir.length() > 0.0) {
        return Err(Error::Failed("a null direction for a guide's origin".into()));
    }
    let mut edges: Vec<WireEdge> = wire.edges.clone();
    let start_at_vertex = |edges: &[WireEdge], key: usize| -> usize {
        edges.iter().position(|e| e.start.is_some_and(|v| v.key == key)).unwrap_or(edges.len() - 1)
    };
    let start = match nearest(wire, p)? {
        Nearest::Vertex(i) => i,
        Nearest::Edge(j, t, x) => {
            let edge = &edges[j];
            let (Some(run_start), Some(run_end)) = (edge.start, edge.end) else { return Err(Error::Failed("a guide edge without vertices".into())) };
            // (The edge's own start and end, whichever way the wire runs it.)
            let (vf, vl) = if edge.reversed { (run_end, run_start) } else { (run_start, run_end) };
            let mut origin: Option<WireVertex> = None;
            let (df, dl) = (p.distance(vf.point), p.distance(vl.point));
            if df < tolerance {
                origin = Some(vl);
            }
            if dl < tolerance && dl < df {
                origin = Some(vf);
            }
            match origin {
                Some(v) => start_at_vertex(&edges, v.key),
                None => {
                    let at = WireVertex { key: keys.fresh(), point: x, tolerance: edge.start.map_or(1e-7, |v| v.tolerance) };
                    let (before, after) = edge.split(t, at);
                    edges.splice(j..=j, [before, after]);
                    j + 1
                }
            }
        }
    };
    let first = &edges[start];
    let (t0, _) = first.run();
    let d = first.run_derivative(t0).ok_or_else(|| Error::Failed("a guide edge without a curve".into()))?;
    if !(d.length() > 0.0) {
        return Err(Error::Failed("a guide with a null derivative at its origin".into()));
    }
    let forward = d.dot(dir) >= 0.0;
    let n = edges.len();
    let order: Vec<WireEdge> = if forward {
        (0..n).map(|k| edges[(start + k) % n].clone()).collect()
    } else {
        (0..=start).rev().chain((start + 1..n).rev()).map(|k| edges[k].reversed()).collect()
    };
    Ok(EdgeWire { edges: order, closed: true })
}

/// Whether the wire counts as straight: on no one plane, and two of its principal radii of
/// gyration equal to the largest within `1e-7`.
fn straight(model: &Model, wire: Shape) -> bool {
    if crate::make::plane_of_wire(model, wire).is_ok() {
        return false;
    }
    let props = crate::properties::length(model, wire);
    if !(props.mass > 0.0) {
        return false;
    }
    let radii = props.principal_moments().map(|m| (m.max(0.0) / props.mass).sqrt());
    let largest = radii[2];
    radii.iter().filter(|&&r| (largest - r).abs() <= 1e-7).count() >= 2
}

/// The ends of an open wire, as points.
fn end_points(wire: &EdgeWire) -> Result<(Point, Point)> {
    match wire.ends() {
        (Some(a), Some(b)) if a.key != b.key => Ok((a.point, b.point)),
        (Some(_), Some(_)) => Err(Error::Failed("a closed wire where an open one is needed".into())),
        _ => Err(Error::Failed("a wire without ends".into())),
    }
}

/// The open `guide` of `model` run the way of the open `spine`: kept when its ends lie
/// nearer the spine's ends taken in order than crosswise (if either wire is straight), or
/// when its chord makes an acute angle with the spine's; else turned.
pub(crate) fn compatible_guide(model: &Model, spine: Shape, guide: Shape, length: &dyn Fn(&GeomAdaptor) -> f64) -> Result<EdgeWire> {
    let _ = length;
    let mut keys = ShapeKeys::new();
    let forward = |s: Shape| Shape { orientation: Orientation::Forward, ..s };
    let spine_wire = EdgeWire::from_model(model, forward(spine), &mut keys)?;
    let guide_wire = EdgeWire::from_model(model, forward(guide), &mut keys)?;
    if spine_wire.closed || guide_wire.closed {
        return Err(Error::Failed("making an open guide run with a closed wire".into()));
    }
    let (sf, sl) = end_points(&spine_wire)?;
    let (gf, gl) = end_points(&guide_wire)?;
    let keep = if straight(model, forward(spine)) || straight(model, forward(guide)) {
        sf.distance(gf) + sl.distance(gl) <= sf.distance(gl) + sl.distance(gf)
    } else {
        let chord = |w: &EdgeWire, f: Point, l: Point| -> Vector {
            if f.distance(l) >= 1e-7 {
                return l - f;
            }
            // (A chord too short to point anywhere: from the first edge's start a hair on.)
            let e = &w.edges[0];
            match e.point(e.first + 1e-7) {
                Some(q) => q - f,
                None => l - f,
            }
        };
        let (s, g) = (chord(&spine_wire, sf, sl), chord(&guide_wire, gf, gl));
        let angle = s.cross(g).length().atan2(s.dot(g));
        angle < std::f64::consts::FRAC_PI_2
    };
    Ok(if keep { guide_wire } else { guide_wire.reversed() })
}

#[cfg(test)]
#[path = "guide_wire_tests.rs"]
mod tests;
